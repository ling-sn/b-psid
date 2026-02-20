#!/usr/bin/env python
'''
Created December 2025
@author: chasew and chat GPT 5.2
Runs HISAT2/STAR aligner, optionally filters against contaminant index, and optionally dedups with UMIs
requires HISAT2 or STAR aligner, samtools, bowtie2, and umi tools

'''

import argparse
import shlex
import subprocess
import sys
from pathlib import Path

def run_cmd(cmd, stdout=None, stderr=None):
	print(">>> " + " ".join(shlex.quote(c) for c in cmd))
	subprocess.run(cmd, check=True, stdout=stdout, stderr=stderr)

def pigz_compress(path, threads=4, keep=False, force=True):
	"""
	Compress a file using pigz.

	Args:
		path (Path or str): file to compress (e.g. .fastq)
		threads (int): number of pigz threads
		keep (bool): keep original file (-k)
		force (bool): overwrite existing .gz (-f)
	"""
	path = Path(path)

	if path.suffix == ".gz":
		print(f">>> pigz: {path.name} already compressed, skipping")
		return

	cmd = ["pigz", "-p", str(threads)]

	if keep:
		cmd.append("-k")
	if force:
		cmd.append("-f")

	cmd.append(str(path))

	run_cmd(cmd)

def list_sample_dirs(input_dir: Path):
	return sorted(d for d in input_dir.iterdir() if d.is_dir())


def select_fastqs(sample_dir: Path):
	"""
	Select paired FASTQs for alignment.

	Priority:
	  1) *_filtered_R1/R2.fastq.gz
	  2) *_fastp_R1/R2.fastq.gz
	  3) *_cutadapt_R1/R2.fastq.gz
	"""
	for tag in ("_filtered_", "_fastp_", "_cutadapt_"):
		r1 = next(sample_dir.glob(f"*{tag}R1.fastq.gz"), None)
		r2 = next(sample_dir.glob(f"*{tag}R2.fastq.gz"), None)
		if r1 and r2:
			return r1, r2, tag.strip("_")

	raise RuntimeError(f"No paired FASTQs found in {sample_dir}")

def infer_history_tag(sample_dir: Path) -> str:
	"""
	Infer preprocessing history from FASTQs present in the sample directory.
	Order is always: cutadapt -> fastp -> filtered
	"""
	history = []

	if list(sample_dir.glob("*_cutadapt_R1.fastq.gz")):
		history.append("cutadapt")

	if list(sample_dir.glob("*_fastp_R1.fastq.gz")):
		history.append("fastp")

	if list(sample_dir.glob("*_filtered_R1.fastq.gz")):
		history.append("filtered")

	return "_".join(history) if history else "raw"

def bowtie2_filter(r1, r2, index, threads, out_bam, unaligned_prefix=None):
	"""
	Run Bowtie2 in filtering mode only.
	Writes:
	  - aligned reads -> out_bam
	  - unaligned paired reads -> <prefix>.1.gz / <prefix>.2.gz
	"""
	index = Path(index).expanduser()

	cmd = [
		"bowtie2",
		"-x", index,
		"-1", str(r1),
		"-2", str(r2),
		"-p", str(threads),
	]

	if unaligned_prefix is not None:
		cmd += ["--un-conc", str(unaligned_prefix)]

	sort_cmd = [
		"samtools", "sort",
		"-@", str(threads),
		"-o", str(out_bam)
	]

	index_cmd = ["samtools", "index", str(out_bam)]

	with subprocess.Popen(cmd, stdout=subprocess.PIPE) as p1:
		with subprocess.Popen(sort_cmd, stdin=p1.stdout) as p2:
			p1.stdout.close()
			p2.communicate()
			if p1.wait() != 0 or p2.returncode != 0:
				raise RuntimeError(f"Bowtie2 filtering failed for {out_bam}")

	run_cmd(index_cmd)

def hisat2_run(r1, r2, index, threads, library_type, out_bam, extra_args=None):
	index = Path(index).expanduser()

	cmd = [
		"hisat2",
		"-x", index,
		"-1", str(r1),
		"-2", str(r2),
		"-p", str(threads),
	]

	if library_type != "unstranded":
		cmd += ["--rna-strandness", library_type]

	if extra_args:
		cmd += extra_args

	sort_cmd = [
		"samtools", "sort",
		"-@", str(threads),
		"-o", str(out_bam)
	]

	index_cmd = ["samtools", "index", str(out_bam)]

	with subprocess.Popen(cmd, stdout=subprocess.PIPE) as p1:
		with subprocess.Popen(sort_cmd, stdin=p1.stdout) as p2:
			p1.stdout.close()
			p2.communicate()
			if p1.wait() != 0 or p2.returncode != 0:
				raise RuntimeError(f"HISAT2 failed for {out_bam}")

	run_cmd(index_cmd)


def star_run(r1, r2, index, threads, out_bam, twopass=False, discovery=False, extra_args=None):
	index = Path(index).expanduser()

	prefix = out_bam.parent / (out_bam.stem + ".")

	cmd = [
		"STAR",
		"--genomeDir", index,
		"--runThreadN", str(threads),
		"--readFilesIn", str(r1), str(r2),
		"--readFilesCommand", "zcat",
		"--outFileNamePrefix", str(prefix),
		"--outSAMtype", "BAM", "SortedByCoordinate",
		"--outSAMstrandField", "intronMotif",
		"--outSAMattributes", "NH", "HI", "AS", "nM", "XS",
		]

	if twopass:
		cmd += ["--twopassMode", "Basic"]

	if discovery:
		cmd += [
			"--alignSJoverhangMin", "6",
			"--outFilterMultimapNmax", "20"
		]

	if extra_args:
		cmd += extra_args

	run_cmd(cmd)

	produced = prefix.with_name(prefix.name + "Aligned.sortedByCoord.out.bam")

	if not produced.exists():
		raise RuntimeError("STAR did not produce expected BAM")

	produced.rename(out_bam)
	run_cmd(["samtools", "index", str(out_bam)])

def write_dedup_slurm(out_file, bam_files, threads, email, account):
	out_file = Path(out_file)
	mem = threads * 3000

	# Convert BAMs into task lines
	new_tasks = [f"\"python parallel_dedup.py --bam {bam} -C {threads}\"" for bam in bam_files]

	# Case 1: file does NOT exist → create fresh
	if not out_file.exists():
		n = len(new_tasks)

		with open(out_file, "w") as fh:
			fh.write(f"""#!/usr/bin/env bash
#SBATCH --job-name=DEDUP
#SBATCH --mail-user={email}
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=DEDUP_%u_%A_%a.out
#SBATCH --array=0-{n-1}
#SBATCH --account={account}
#SBATCH --time=1:30:00
#SBATCH --mem={mem}m
#SBATCH --cpus-per-task={threads}

module purge
eval "$(conda shell.bash hook)"
conda activate ~/miniconda3/envs/RNA-SEQ

declare -a tasks=(
""")

			for t in new_tasks:
				fh.write(t + "\n")

			fh.write(""")
eval ${tasks[$SLURM_ARRAY_TASK_ID]}
""")
		return

	# Case 2: file EXISTS → edit in place
	with open(out_file) as fh:
		lines = fh.readlines()

	header = []
	tasks = []
	footer = []

	in_tasks = False

	for line in lines:
		if line.startswith("declare -a tasks=("):
			in_tasks = True
			continue
		elif in_tasks and line.strip() == ")":
			in_tasks = False
			continue

		if in_tasks:
			tasks.append(line.rstrip())
		elif line.strip().startswith("eval ${tasks"):
			footer.append(line)
		else:
			header.append(line)

	# Append new tasks if missing
	existing = set(tasks)
	for t in new_tasks:
		if t not in existing:
			tasks.append(t)

	# Update array line
	ntasks = len(tasks)
	new_header = []
	for line in header:
		if line.startswith("#SBATCH --array="):
			new_header.append(f"#SBATCH --array=0-{ntasks-1}\n")
		else:
			new_header.append(line)

	# Rewrite file
	with open(out_file, "w") as fh:
		fh.writelines(new_header)
		fh.write("declare -a tasks=(\n")
		for t in tasks:
			fh.write(t + "\n")
		fh.write(")\n")
		fh.writelines(footer)

if __name__ == '__main__': # allows another python script to import the functions
	parser = argparse.ArgumentParser(description="Align per-sample paired FASTQs (fastp preferred, cutadapt fallback) using HISAT2 or STAR.")

	required = parser.add_argument_group('Required Input', 'Specify input and output directories.')

	required.add_argument("--input", required=True, help="Directory containing per-sample folders.")
	required.add_argument("--output", required=True, help="Output directory for BAM files.")
	required.add_argument("--aligner", required=True, choices=["hisat2", "star"])
	required.add_argument("--index", required=True, help="HISAT2 index prefix or STAR genomeDir.")

	data_opt = parser.add_argument_group('Basic Arguments', 'These options can be used to change how the script is generally run.')

	data_opt.add_argument("-C","--threads", type=int, default=4, help="Specify the number of cpus. Default is 4")
	data_opt.add_argument("-L","--library_type", choices=["RF", "FR", "unstranded"], default="unstranded", help="Specify which strand is on read 1 (RF is AS, FR is S). Default is unstranded (but NEBNext needs RF)")
	data_opt.add_argument("-O","--overwrite", action="store_true", help="Specify to overwrite existing outputs")
	data_opt.add_argument("-S","--sample", help="Process only a single sample directory within the input directory")
	data_opt.add_argument("--extra-aligner-args", default="", help=argparse.SUPPRESS)

	filter_opt = parser.add_argument_group("Alignment-based filtering", "Align against a limited reference and keep only unaligned reads for downstream processing.")

	filter_opt.add_argument("-F","--filter_index", help="If set, align reads to this bowtie2 index and retain unaligned reads for downstream processing.")
	filter_opt.add_argument("-f","--filter_only", action="store_true", help="Only perform alignment-based filtering and do not run final alignment.")

	star_opt = parser.add_argument_group("STAR options", "Options specific to STAR alignment.")

	star_opt.add_argument("-T","--star_twopass", action="store_true", help="Enable STAR 2-pass alignment (recommended for transcript discovery)")
	star_opt.add_argument("-D","--star_discovery", action="store_true", help="Enable more discovery-friendly STAR settings (at expense of more ambiguous mapped reads)")

	dedup_opt = parser.add_argument_group('Dedup scripting', 'These options can be used to help setup dedup with umi tools.')
	parser.add_argument("--emit_dedup_slurm", metavar="FILE", help="Write a SLURM array script to deduplicate all produced BAMs")
	parser.add_argument("--slurm_email", default="<uniqname>@umich.edu", help=argparse.SUPPRESS)
	parser.add_argument("--slurm_account", default="<account>", help=argparse.SUPPRESS)
	parser.add_argument("--dedup_threads", type=int, default=8, help=argparse.SUPPRESS)

	args = parser.parse_args()

	if args.filter_only and not args.filter_index:
		parser.error("--filter_only requires --filter_index")

	input_dir = Path(args.input).resolve()
	output_dir = Path(args.output).resolve()
	output_dir.mkdir(parents=True, exist_ok=True)

	extra_args = shlex.split(args.extra_aligner_args) if args.extra_aligner_args else []

	final_bams = []

	requested_sample = args.sample

	for sample_dir in list_sample_dirs(input_dir):
		sample = sample_dir.name

		if requested_sample and sample != requested_sample:
			continue

		try:
			r1, r2, stage = select_fastqs(sample_dir)
		except RuntimeError as e:
			print(f">>> Skipping {sample}: {e}")
			continue

		if args.filter_index:
			# Look for existing filtered FASTQs
			filtered_r1 = next(sample_dir.glob(f"{sample}*_filtered_R1.fastq.gz"), None)
			filtered_r2 = next(sample_dir.glob(f"{sample}*_filtered_R2.fastq.gz"), None)

			if filtered_r1 and filtered_r2 and not args.overwrite:
				print(f">>> Found existing filtered FASTQs for {sample}, skipping contaminant filter")
				r1, r2, stage = filtered_r1, filtered_r2, "filtered"

			else:
				pre_filter_tag = infer_history_tag(sample_dir)
				print(f">>> Filtering contaminants for {sample}")

				filtered_prefix = sample_dir / f"{sample}_{pre_filter_tag}_filtered"
				contam_bam = sample_dir / f"{sample}_{pre_filter_tag}_contaminants.bam"

				bowtie2_filter(
					r1=r1,
					r2=r2,
					index=args.filter_index,
					threads=args.threads,
					out_bam=contam_bam,
					unaligned_prefix=filtered_prefix
				)

				# Rename Bowtie2 outputs
				(filtered_prefix.with_suffix(".1")).rename(sample_dir / f"{sample}_{pre_filter_tag}_filtered_R1.fastq")
				(filtered_prefix.with_suffix(".2")).rename(sample_dir / f"{sample}_{pre_filter_tag}_filtered_R2.fastq")

				pigz_compress(sample_dir / f"{sample}_{pre_filter_tag}_filtered_R1.fastq", threads=args.threads)
				pigz_compress(sample_dir / f"{sample}_{pre_filter_tag}_filtered_R2.fastq", threads=args.threads)

				# Re-detect filtered FASTQs after creation
				r1, r2, stage = select_fastqs(sample_dir)

				if args.filter_only:
					print(f">>> Filter-only mode: skipping final alignment for {sample}")
					continue

		history_tag = infer_history_tag(sample_dir)
		sample_outdir = output_dir / sample
		sample_outdir.mkdir(parents=True, exist_ok=True)

		out_bam = sample_outdir / f"{sample}_{history_tag}_{args.aligner}.bam"

		if out_bam.exists() and not args.overwrite:
			final_bams.append(out_bam)
			print(f">>> Skipping {sample} alignment: BAM exists")
			continue

		print(f"\n=== {sample} ===")
		print(f">>> Stage: {stage}")
		print(f">>> R1: {r1.name}")
		print(f">>> R2: {r2.name}")
		print(f">>> Aligner: {args.aligner}")
		print(f">>> Output: {out_bam}")

		if args.aligner == "hisat2":
			hisat2_run(
				r1=r1,
				r2=r2,
				index=args.index,
				threads=args.threads,
				library_type=args.library_type,
				out_bam=out_bam,
				extra_args=extra_args
			)

		else:
			if args.star_twopass:
				print(">>> STAR 2-pass mode enabled")
			if args.star_discovery:
				print(">>> STAR discovery mode enabled")
			star_run(
				r1=r1,
				r2=r2,
				index=args.index,
				threads=args.threads,
				out_bam=out_bam,
				twopass=args.star_twopass,
				discovery=args.star_discovery,
				extra_args=extra_args
			)

		final_bams.append(out_bam)

	if args.sample and not final_bams:
		sys.exit(f"ERROR: sample '{args.sample}' not found in {input_dir}")

	if args.emit_dedup_slurm:
		write_dedup_slurm(
			out_file=args.emit_dedup_slurm,
			bam_files=final_bams,
			threads=args.dedup_threads,
			email=args.slurm_email,
			account=args.slurm_account
		)
		print(f">>> Wrote dedup SLURM script to {args.emit_dedup_slurm}")

	print("\n>>> Done.")
