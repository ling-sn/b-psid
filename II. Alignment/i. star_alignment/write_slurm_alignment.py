## Use RNA-SEQ conda environment
from pathlib import Path
import traceback
import argparse
import textwrap
import os

def main(input_folder: str, output_folder: str, aligner_type: str,
         genome_idx: str, filter_idx: str, library: str,
         two_pass: bool, emit_dedup: str, email: str, 
         slurm_acct: str, walltime: str, mem: int):
    '''
    PURPOSE:
    * Goes into folder containing all subfolders for trimmed reads,
      then obtains all sample names.
    * Automatically appends new task containing sample name to SBATCH script,
      which then aligns reads using HISAT2 or STAR.
    '''
    current_path = Path.cwd()
    output = current_path/"SBATCHSubArr-Align-STAR.sbatch"
    
    start_dir = Path(current_path/input_folder)
    if not start_dir.exists():
        raise FileNotFoundError(
            "Please input a folder name that exists "
            "in your current working directory."
        )
    
    '''
    1. Count all subfolders in start_dir
    2. Subtract 1 so count is 0-based
    '''
    num_jobs = len(next(os.walk(start_dir))[1]) - 1
    
    ## Create SBATCH file if it doesn't exist
    if not output.exists():
        with open(output, "w") as f:
            template_start = textwrap.dedent(f"""\
                                            #!/usr/bin/env bash
                                            #SBATCH --job-name=ALIGN
                                            #SBATCH --mail-user={email}
                                            #SBATCH --mail-type=BEGIN,END,FAIL
                                            #SBATCH --output=ALIGN_%u_%A_%a.out
                                            #SBATCH --array=0-{num_jobs}
                                            #SBATCH --account={slurm_acct}
                                            #SBATCH --time={walltime}
                                            #SBATCH --mem={mem}m
                                            #SBATCH --partition=standard
                                            #SBATCH --cpus-per-task=8
                                            ################################################################################
                                            # Edit the strings under 'declare -a tasks=(' to match your experiments.
                                            #
                                            # Recommend 1.5 hours per 30M read sample, 2.5 for 2-pass STAR.
                                            #
                                            # This requires a conda environment for genome alignment, edit to your named
                                            # version in the activate command.
                                            # 
                                            # To call this script:
                                            # sbatch SBATCHSubArr-Align-STAR.sbatch
                                            ################################################################################

                                            module purge
                                            eval "$(conda shell.bash hook)"
                                            conda activate ~/miniconda3/envs/RNA-SEQ

                                            declare -a tasks=(
                                            """)

            f.write(template_start)

    try:
        for subfolder in start_dir.iterdir():
            if subfolder.is_dir():
                ## Append new tasks to SBATCH
                sample_name = str(subfolder.stem)
                with open(output, "a") as f:
                    task = (f'\n"python3 -u run_align.py --input {input_folder} --output {output_folder} '
                            f'--aligner {aligner_type} --index {genome_idx} '
                            f'-C 8 -L {library} -S {sample_name}')
                    
                    if (filter_idx):
                        task += f' --filter_index {filter_idx}'

                    if (two_pass):
                        task += f' -T --emit_dedup_slurm {emit_dedup}'
                    
                    task += '"'

                    f.write(task)

        ## Once all tasks have been added, finish up SBATCH template
        with open(output, "a") as f:
            f.write("\n)\neval ${tasks[$SLURM_ARRAY_TASK_ID]}")
    except Exception as e:
        print("Failed to create SBATCH file: {e}")
        traceback.print_exc()
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Writes SBATCH script for STAR alignment.")
    parser.add_argument("--input_folder", help = "Name of folder containing all trimmed reads", 
                        required = True)
    parser.add_argument("--output_folder", help = "Name of output folder for aligned reads (after running run_align)", 
                        default = "star_aligned")
    parser.add_argument("--aligner_type", help = "Specify aligner type", choices = ["hisat2", "star"], required = True)
    parser.add_argument("--genome_idx", help = "Directory to genome index", required = True)
    parser.add_argument("--filter_idx", help = "Directory to contaminants index, including prefix of index (e.g., everything up to .1.bt2)")
    parser.add_argument("--library", help = "Specify which strand is on Read 1. Default is unstranded, but NEBNext needs RF",
                        choices = ["RF", "FR", "unstranded"], required = True)
    parser.add_argument("--two_pass", action = "store_true")
    parser.add_argument("--emit_dedup", help = "Name of output SBATCH file for deduplication (next step of pipeline)")
    parser.add_argument("--email", default = "<uniqname>@umich.edu")
    parser.add_argument("--slurm_acct", default = "<account>")
    parser.add_argument("--walltime", default = "<time>")
    parser.add_argument("--mem", help = "Memory for job (in MB)", default = "<memory>")

    args = parser.parse_args()

    print("Writing SBATCH script...")
    main(args.input_folder, args.output_folder, args.aligner_type, args.genome_idx,
         args.filter_idx, args.library, args.two_pass, args.emit_dedup,
         args.email, args.slurm_acct, 
         args.walltime, args.mem)
    print("Process finished.")