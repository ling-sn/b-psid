## Use RNA-STAR conda environment
from pathlib import Path
import traceback
import argparse
import textwrap
import os

def main(star_folder: Path, fasta: Path, discard: bool, 
         email: str, slurm_acct: str, walltime: str, mem: int):
    '''
    PURPOSE:
    * Goes into folder containing all subfolders for deduplicated reads,
      then obtains all sample names.
    * Automatically appends new task containing sample name to SBATCH script,
      which then realigns reads to be more deletion-aware.
    '''
    current_path = Path.cwd()
    output = current_path/"realignGap.sbatch"
    
    start_dir = Path(current_path/star_folder)
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
                                            #SBATCH --job-name=REALIGN
                                            #SBATCH --mail-user={email}
                                            #SBATCH --mail-type=BEGIN,END,FAIL
                                            #SBATCH --output=REALIGN_%u_%A_%a.out
                                            #SBATCH --array=0-{num_jobs}
                                            #SBATCH --account={slurm_acct}
                                            #SBATCH --time={walltime}
                                            #SBATCH --mem={mem}m
                                            #SBATCH --partition=standard
                                            #SBATCH --ntasks-per-node=1
                                            #SBATCH --nodes=1
                                            ################################################################################
                                            # Edit the strings under 'declare -a tasks=(' to match your experiments.
                                            #
                                            # This requires a conda environment with samtools, pysam, and parasail (RNA-STAR)
                                            # 
                                            # To call this script:
                                            # sbatch realignGap.sbatch
                                            ################################################################################

                                            module purge
                                            eval "$(conda shell.bash hook)"
                                            conda activate ~/miniconda3/envs/RNA-STAR

                                            declare -a tasks=(
                                            """)

            f.write(template_start)

    try:
        for subfolder in start_dir.iterdir():
            if subfolder.is_dir():
                ## Append new tasks to SBATCH
                folder = str(subfolder.stem)
                with open(output, "a") as f:
                    task = (f'\n"python3 realignGap.py --star_folder {star_folder} --subf {folder} --fasta_dir {fasta}')

                    if (discard):
                        task += f' --discard'
                    
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
    parser = argparse.ArgumentParser(description = "Writes SBATCH script for STAR realignment.")
    parser.add_argument("--star_folder", help = "Path to folder created from STAR alignment", 
                        required = True)
    parser.add_argument("--fasta", help = "Path of reference fasta file", required = True)
    parser.add_argument("--discard", 
                        help = "Write discarded reads into file for debugging. Use '--no-discard' to disable and '--discard' to enable.",
                        action = "store_true", default = False)
    parser.add_argument("--email", default = "<uniqname>@umich.edu")
    parser.add_argument("--slurm_acct", default = "<account>")
    parser.add_argument("--walltime", default = "<time>")
    parser.add_argument("--mem", help = "Memory for job (in MB)", default = "<memory>")

    args = parser.parse_args()

    print("Writing SBATCH script...")
    main(args.star_folder, args.fasta, args.discard,
         args.email, args.slurm_acct, args.walltime, args.mem)
    print("Process finished.")