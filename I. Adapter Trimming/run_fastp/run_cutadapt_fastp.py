#!/usr/bin/env python
'''
Created December 2025
@author: chasew and chat GPT 5.2
Runs cutadapt and fastP

'''
import os, sys, argparse
import subprocess
import shutil
import tempfile
import re

FASTQ_EXT_RE = r"(?:\.fastq|\.fq)(?:\.gz)?$"

R1_RE = re.compile(rf"(.*?)(_R?1(?:_\d+)?){FASTQ_EXT_RE}", re.IGNORECASE)
R2_RE = re.compile(rf"(.*?)(_R?2(?:_\d+)?){FASTQ_EXT_RE}", re.IGNORECASE)

def parse_fastq_read(filename):
	if R1_RE.match(filename):
		return 1
	if R2_RE.match(filename):
		return 2
	return None
	
def get_sample_id(filename):
	return filename.split("_", 1)[0]

def collect_sample_pairs(input_dir):
	pairs = {}

	for fname in os.listdir(input_dir):
		read = parse_fastq_read(fname)
		if read is None:
			continue

		sample = get_sample_id(fname)
		pairs.setdefault(sample, {})[read] = os.path.join(input_dir, fname)

	return pairs

def run_cutadapt_single_sample(sample, r1, r2, output_dir, adapter_fasta, threads=2, error_rate=0.1, minimum_length=30, overlap=8, extra_args=None):
	os.makedirs(output_dir, exist_ok=True)

	r1_out = os.path.join(output_dir, f"{sample}_cutadapt_R1.fastq.gz")
	r2_out = os.path.join(output_dir, f"{sample}_cutadapt_R2.fastq.gz")
	json_out = os.path.join(output_dir, f"{sample}_cutadapt.json")

	cmd = [
		"cutadapt",
		"-j", str(threads),
		"-a", f"file:{adapter_fasta}",
		"-A", f"file:{adapter_fasta}",
		"--overlap", str(overlap),
		"-e", str(error_rate),
		"--minimum-length", str(minimum_length),
		"--json", json_out,
		"-o", r1_out,
		"-p", r2_out,
		r1, r2
	]

	if extra_args:
		cmd += extra_args

	subprocess.run(cmd, check=True)

def run_fastp_single_sample(sample, r1, r2, output_dir, umi=True, umi_loc="per_read", umi_len=12, threads=2, window_trim=4, meanQ_trim=15, min_read_length=50, meanQ_total=15, pcent_passQ=30, adapters_already_trimmed=True):
	sample_dir = os.path.join(output_dir, sample)
	os.makedirs(sample_dir, exist_ok=True)

	r1_out = os.path.join(sample_dir, f"{sample}_fastp_R1.fastq.gz")
	r2_out = os.path.join(sample_dir, f"{sample}_fastp_R2.fastq.gz")
	json_out = os.path.join(sample_dir, f"{sample}_fastp.json")
	html_out = os.path.join(sample_dir, f"{sample}_fastp.html")

	cmd = [
		"fastp",
		"-i", r1,
		"-I", r2,
		"-o", r1_out,
		"-O", r2_out,
		"-j", json_out,
		"-h", html_out,
		"-w", str(threads),
		"-l", str(min_read_length),
		"-5", "-3",
		"-W", str(window_trim),
		"-M", str(meanQ_trim),
		"-q", str(meanQ_total),
		"-u", str(pcent_passQ),
		"-c",
		"-p"
	]

	if umi:
		cmd += [
			"-U",
			"--umi_loc", umi_loc,
			"--umi_len", str(umi_len),
			"--umi_prefix", "UMI"
		]

	if adapters_already_trimmed:
		cmd += ["-A", "-G"]

	subprocess.run(cmd, check=True)

def emit_fastqs(sample, source_r1, source_r2, output_dir, prefix):
	os.makedirs(output_dir, exist_ok=True)

	destination_r1 = os.path.join(output_dir, f"{sample}_{prefix}_R1.fastq.gz")
	destination_r2 = os.path.join(output_dir, f"{sample}_{prefix}_R2.fastq.gz")

	shutil.copy2(source_r1, destination_r1)
	shutil.copy2(source_r2, destination_r2)

def emit_json(src, dst):
	os.makedirs(os.path.dirname(dst), exist_ok=True)
	shutil.copy2(src, dst)

if __name__ == '__main__': # allows another python script to import the functions

	parser = argparse.ArgumentParser(description="Runs fastP", usage='%(prog)s [options] --input in_dir --output out_dir',add_help=False) # a description of the function

	required = parser.add_argument_group('Required Input', 'Specify input and output directories.')

	required.add_argument("--input",help="path to directory with paired fastqs.", required=True)
	required.add_argument("--output",help="path to new directory to place UMI_labeled fastqs.", required=True)

	data_opt = parser.add_argument_group('Basic Arguments', 'These options can be used to change how the script is generally run.')

	data_opt.add_argument("-X","--cutadapt_mode", action="count",help="Control Cutadapt behavior (default = Cutadapt first): use -X for Cutadapt only, -XX to skip Cutadapt.")
	data_opt.add_argument("-C","--threads",type=int,default=2,help="Specify number of cpus. Default is 2.")
	data_opt.add_argument("-S","--sample", help="Process only a single sample prefix within the input directory")

	cut_opt = parser.add_argument_group('Cutadapt Options', 'These options can be used to change how cutadapt runs.')
	cut_opt.add_argument("-b","--barcode_fasta",default="ARL_adapters.fa",help="Specify barcode primers fasta. Default is ARL_adapters.fa")
	cut_opt.add_argument("-f","--Tn_fasta",default="Tn_adapters.fa",help="Specify Tn mosaic ends fasta. Default is Tn_adapters.fa")
	cut_opt.add_argument("-e","--cut_error",type=float,default=0.1,help="Specify max error rate when finding adapters to trim. Default = 0.1")
	cut_opt.add_argument("-m","--min_cutadapt_read",type=int,default=30,help="Specify minimum read length to keep after adapter trimming. Default = 30")
	cut_opt.add_argument("-o","--overlap",type=int,default=8,help="Specify minimum read-adapter overlap for trimming. Default = 8")
	# Hidden expert escape hatch for arbitrary cutadapt flags
	parser.add_argument("--cutadapt-extra",nargs="+",help=argparse.SUPPRESS)

	fastp_opt = parser.add_argument_group('Fastp Options', 'These options can be used to change how fastp runs.')
	fastp_opt.add_argument("-W","--window_trim",type=int,default=4,help="Specify the window size when end trimming. Default is 4")
	fastp_opt.add_argument("-T","--meanQ_trim",type=int,default=15,help="Specify the mean Q of a window when trimming necessary to keep. Default is 15")
	fastp_opt.add_argument("-L","--min_length",type=int,default=30,help="Specify minimum read length to keep. Default = 30")
	fastp_opt.add_argument("-M","--meanQ_total",type=int,default=15,help="Specify Q score to serve as keep threshold. Default is 15")
	fastp_opt.add_argument("-P","--pcent_passQ",type=int,default=30,help="Specify percent of read nucleotides that must pass meanQ to keep. Default is 30")
	fastp_opt.add_argument("-U","--umi_len",type=int,default=None,help="Length of UMI at 5' end of reads. Enables UMI extraction if set.")

	help_opt = parser.add_argument_group('Help')
	help_opt.add_argument('-h', '--help', action="help", help="show this help message and exit")

	args = parser.parse_args()

	mode = args.cutadapt_mode or 0
	umi_enabled = args.umi_len is not None

	pairs = collect_sample_pairs(args.input)

	if args.sample:
		if args.sample not in pairs:
			sys.exit(f">>> Requested sample '{args.sample}' not found in {args.input}")
		pairs = {args.sample: pairs[args.sample]}

	for sample, reads in pairs.items():
		if 1 not in reads or 2 not in reads:
			print(f">>> Skipping {sample}: missing R1 or R2")
			continue

		print(f">>> Processing sample {sample}")

		with tempfile.TemporaryDirectory(prefix=f"cutadapt_{sample}_") as tmp:
			pass1_dir = os.path.join(tmp, "pass1")
			pass2_dir = os.path.join(tmp, "pass2")

			# -------------------------
			# MODE 0 or 1: run cutadapt
			# -------------------------

			if mode in (0, 1):
				# --- Cutadapt pass 1 ---
				run_cutadapt_single_sample(
					sample,
					reads[1],
					reads[2],
					pass1_dir,
					args.barcode_fasta,
					threads=args.threads,
					error_rate=args.cut_error,
					minimum_length=args.min_cutadapt_read,
					overlap=args.overlap,
					extra_args=args.cutadapt_extra
				)

				# --- Cutadapt pass 2 ---
				r1_p1 = os.path.join(pass1_dir, f"{sample}_cutadapt_R1.fastq.gz")
				r2_p1 = os.path.join(pass1_dir, f"{sample}_cutadapt_R2.fastq.gz")

				run_cutadapt_single_sample(
					sample,
					r1_p1,
					r2_p1,
					pass2_dir,
					args.Tn_fasta,
					threads=args.threads,
					error_rate=args.cut_error,
					minimum_length=args.min_cutadapt_read,
					overlap=args.overlap,
					extra_args=args.cutadapt_extra
				)
				
				r1_p2 = os.path.join(pass2_dir, f"{sample}_cutadapt_R1.fastq.gz")
				r2_p2 = os.path.join(pass2_dir, f"{sample}_cutadapt_R2.fastq.gz")

				sample_dir = os.path.join(args.output, sample)
				emit_fastqs(sample, r1_p2, r2_p2, sample_dir, prefix="cutadapt")
				emit_json(os.path.join(pass1_dir, f"{sample}_cutadapt.json"), os.path.join(sample_dir, f"{sample}_cutadapt_pass1.json"))
				emit_json(os.path.join(pass2_dir, f"{sample}_cutadapt.json"), os.path.join(sample_dir, f"{sample}_cutadapt_pass2.json"))		

			# -------------------------
			# MODE 1: cutadapt only
			# -------------------------
			if mode == 1:
				print(f">>> Cutadapt-only mode for {sample}, skipping fastp")
				continue

			# -------------------------
			# MODE 2: fastp only
			# -------------------------
			if mode == 2:
				r1_in = reads[1]
				r2_in = reads[2]
				adapters_trimmed = False

			else:
				r1_in = os.path.join(pass2_dir, f"{sample}_cutadapt_R1.fastq.gz")
				r2_in = os.path.join(pass2_dir, f"{sample}_cutadapt_R2.fastq.gz")
				adapters_trimmed = True

			# -------------------------
			# Fastp (mode 0 or 2)
			# -------------------------
			
			# --- Fastp ---
			run_fastp_single_sample(
				sample,
				r1_in,
				r2_in,
				args.output,
				umi=umi_enabled,
				umi_loc="per_read",
				umi_len=args.umi_len,
				threads=args.threads,
				window_trim=args.window_trim,
				meanQ_trim=args.meanQ_trim,
				min_read_length=args.min_length,
				meanQ_total=args.meanQ_total,
				pcent_passQ=args.pcent_passQ,
				adapters_already_trimmed=adapters_trimmed
			)

		# Rename fastp outputs in default mode to reflect full history
		if mode == 0:
			sample_dir = os.path.join(args.output, sample)

			os.rename(
				os.path.join(sample_dir, f"{sample}_fastp_R1.fastq.gz"),
				os.path.join(sample_dir, f"{sample}_cutadapt_fastp_R1.fastq.gz")
			)
			os.rename(
				os.path.join(sample_dir, f"{sample}_fastp_R2.fastq.gz"),
				os.path.join(sample_dir, f"{sample}_cutadapt_fastp_R2.fastq.gz")
			)
