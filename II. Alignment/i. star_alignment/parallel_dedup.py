#!/usr/bin/env python

"""
Created December 2025
@author: chasew and ChatGPT 5.2

Runs umi_tools deduplication on an aligned, coordinate-sorted BAM
with UMIs encoded in read names. Splits process into chromosomes.
Requires samtools.
"""

import argparse
import shlex
import subprocess
import sys
import shutil
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import re
import json


UMI_STATS_RE = {
	"input_reads": re.compile(r"Input Reads:\s*(\d+)"),
	"read_pairs": re.compile(r"Read pairs:\s*(\d+)"),
	"output_reads": re.compile(r"Number of reads out:\s*(\d+)"),
	"positions_deduped": re.compile(r"Total number of positions deduplicated:\s*(\d+)"),
	"mean_umis": re.compile(r"Mean number of unique UMIs per position:\s*([\d.]+)"),
	"max_umis": re.compile(r"Max\. number of unique UMIs per position:\s*(\d+)")
}

def parse_umi_tools_stats(text: str):
	stats = {
		"input_reads": None,
		"read_pairs": None,
		"output_reads": None,
		"positions_deduped": None,
		"mean_umis": None,
		"max_umis": None,
	}

	for key, rx in UMI_STATS_RE.items():
		m = rx.search(text)
		if m:
			val = m.group(1)
			stats[key] = float(val) if "." in val else int(val)

	# Derived fields
	if stats["input_reads"] is not None and stats["output_reads"] is not None:
		stats["duplicates_removed"] = stats["input_reads"] - stats["output_reads"]
		stats["duplication_rate"] = stats["duplicates_removed"] / stats["input_reads"]
	else:
		stats["duplicates_removed"] = None
		stats["duplication_rate"] = None

	return stats

def run_cmd(cmd, stdout=None, stderr=None):
	print(">>> " + " ".join(shlex.quote(c) for c in cmd))
	subprocess.run(cmd, check=True, stdout=stdout, stderr=stderr)

def get_chromosomes(bam: Path):
	cmd = ["samtools", "idxstats", str(bam)]
	out = subprocess.check_output(cmd, text=True)
	chroms = []
	for line in out.strip().splitlines():
		chrom, length, *_ = line.split("\t")
		if chrom != "*" and int(length) > 0:
			chroms.append(chrom)
	return chroms

def split_bam_by_chrom(bam: Path, chroms, tmpdir: Path):
	tmpdir.mkdir(parents=True, exist_ok=True)
	bam_files = []

	for chrom in chroms:
		out = tmpdir / f"{chrom}.bam"
		run_cmd([
			"samtools", "view",
			"-b", str(bam),
			chrom,
			"-o", str(out)
		])

		count = subprocess.check_output(["samtools", "view", "-c", str(out)], text=True).strip()
		if count == "0":
			continue

		run_cmd(["samtools", "index", str(out)])
		bam_files.append(out)

	return bam_files

def umi_dedup(in_bam: Path, out_bam: Path, stats_dir: Path):
	"""
	Run umi_tools dedup and capture UMI summary stats.
	"""
	stats_dir.mkdir(parents=True, exist_ok=True)

	chrom = in_bam.stem
	stats_file = stats_dir / f"{chrom}.json"

	cmd = [
		"umi_tools", "dedup",
		"--paired",
		"--extract-umi-method=read_id",
		"--umi-separator=:UMI_",
		"-I", str(in_bam),
		"-S", str(out_bam),
	]

	print(">>> " + " ".join(shlex.quote(c) for c in cmd))
	proc = subprocess.run(
		cmd,
		check=True,
		stdout=subprocess.PIPE,
		stderr=subprocess.PIPE,
		text=True
	)

	combined_log = proc.stderr + "\n" + proc.stdout
	stats = parse_umi_tools_stats(combined_log)
	stats["chrom"] = chrom

	with open(stats_file, "w") as fh:
		json.dump(stats, fh)

def dedup_one_chrom(args):
	bam, stats_dir = args
	out = bam.with_name(bam.stem + ".dedup.bam")
	umi_dedup(bam, out, stats_dir)
	return out

def parallel_dedup(bam_files, threads, stats_dir):
	deduped = []
	max_workers = min(threads, len(bam_files))
	with ProcessPoolExecutor(max_workers=max_workers) as exe:
		args = [(bam, stats_dir) for bam in bam_files]
		for out in exe.map(dedup_one_chrom, args):
			deduped.append(out)
	return deduped    

def merge_bams(bams, out_bam):
	run_cmd([
		"samtools", "merge",
		"-c",
		"-p",
		str(out_bam),
		*map(str, bams)
	])

def summarize_umi_stats(stats_dir: Path, out_tsv: Path):
	rows = []
	total_input = 0
	total_output = 0

	for f in sorted(stats_dir.glob("*.json")):
		with open(f) as fh:
			s = json.load(fh)
		rows.append(s)

		if s["input_reads"] is not None:
			total_input += s["input_reads"]
		if s["output_reads"] is not None:
			total_output += s["output_reads"]

	total_duplicates = total_input - total_output
	total_dup_rate = total_duplicates / total_input if total_input else 0.0

	with open(out_tsv, "w") as out:
		out.write(
			"chrom\tinput_reads\toutput_reads\tduplicates_removed\tduplication_rate\t"
			"positions_deduped\tmean_umis\tmax_umis\n"
		)

		for r in rows:
			out.write(
				f"{r['chrom']}\t"
				f"{r['input_reads']}\t"
				f"{r['output_reads']}\t"
				f"{r['duplicates_removed']}\t"
				f"{r['duplication_rate']:.4f}\t"
				f"{r['positions_deduped']}\t"
				f"{r['mean_umis']}\t"
				f"{r['max_umis']}\n"
			)

		out.write(
			f"TOTAL\t{total_input}\t{total_output}\t"
			f"{total_duplicates}\t{total_dup_rate:.4f}\t.\t.\t.\n"
		)

def umi_dedup_parallel(in_bam, out_bam, threads=4, tmpdir=None):
	tmpdir = tmpdir or in_bam.parent / f"dedup_tmp_{in_bam.stem}"
	stats_dir = tmpdir / "umi_stats"

	chroms = get_chromosomes(in_bam)
	chrom_bams = split_bam_by_chrom(in_bam, chroms, tmpdir)
	deduped = parallel_dedup(chrom_bams, threads, stats_dir)
	merge_bams(deduped, out_bam)

	summarize_umi_stats(stats_dir, out_bam.with_suffix(".umi_stats.tsv"))

	shutil.rmtree(tmpdir)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Deduplicate an aligned BAM using umi_tools")

	required = parser.add_argument_group('Required Input', 'Specify input and output directories.')

	required.add_argument("--bam", required=True, help="Path to coordinate-sorted BAM file with UMIs in read names")

	data_opt = parser.add_argument_group('Basic Arguments', 'These options can be used to change how the script is generally run.')

	data_opt.add_argument("-C","--threads", type=int, default=4, help="Specify the number of cpus. Default is 4")

	args = parser.parse_args()

	bam = Path(args.bam).resolve()

	if not bam.exists():
		sys.exit(f"ERROR: BAM file not found: {bam}")

	dedup_bam = bam.with_name(bam.stem + "_dedup.bam")

	print(f">>> Deduplicating {bam.name}")

	umi_dedup_parallel(in_bam=bam, out_bam=dedup_bam, threads=args.threads)

	tmp_sorted = dedup_bam.with_suffix(".sorted.bam")
	run_cmd(["samtools", "sort", "-o", str(tmp_sorted), str(dedup_bam)])
	tmp_sorted.replace(dedup_bam)
	run_cmd(["samtools", "index", str(dedup_bam)])