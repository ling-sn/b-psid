## Use RNA-STAR conda environment
from pathlib import Path
import traceback
import argparse
import textwrap

def main(email, slurm_acct, walltime, mem, fa):
    """
    PURPOSE:
    1. Goes into each sample group folder and obtains all subfolder names.
    2. Automatically adds new task containing folder and subfolder names
       to SBATCH script, allowing for BAM parallelization downstream 
       (each subfolder = 1 BAM).
    """
    current_path = Path.cwd()
    output = Path("run_calculate_dr.sbatch")
    
    start_dir = current_path/"realignments"
    if not start_dir.exists():
        raise FileNotFoundError(
            "Realignments folder doesn't exist. Did you run realignGap.py?"
        )

    """
    - Recursive count of directories = Number of jobs in array
    - Subtract 1 so it's 0-based
    """
    num_jobs = len(list(start_dir).rglob("*/")) - 1
    
    ## Create SBATCH file if it doesn't exist
    if not output.exists():
        with open(output, "w") as f:
            template_start = textwrap.dedent(f"""\
                                            #!/usr/bin/env bash
                                            #SBATCH --job-name=calculate_dr
                                            #SBATCH --mail-user={email}
                                            #SBATCH --mail-type=BEGIN,END,FAIL
                                            #SBATCH --output=calculate_dr_%u_%A_%a.out
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
                                            # The #SBATCH --array variable above creates an array [0,1,2,3]. Change it so that the length
                                            # is how many jobs you need (same as number of strings under $tasks).
                                            #
                                            # This script is submitted that many times, but only one line from $tasks is
                                            # evaluated each time.
                                            #
                                            # For more info on #SBATCH variables, see https://arc.umich.edu/greatlakes/slurm-user-guide/
                                            # and https://slurm.schedmd.com/sbatch.html
                                            #
                                            # This requires a conda environment with samtools and pysam (RNA-STAR)
                                            # 
                                            # To call this script:
                                            # sbatch run_calculate_dr.sbatch
                                            ################################################################################

                                            module purge
                                            eval "$(conda shell.bash hook)"
                                            conda activate ~/miniconda3/envs/RNA-STAR

                                            declare -a tasks=(
                                            """)

            f.write(template_start)

    try:
        for folder in start_dir.iterdir():
            ## Obtain all subfolder names in each folder within 'realignments'
            subfolders = [subf.stem for subf in folder.iterdir() if subf.is_dir()]

            ## Append new tasks to SBATCH
            with open(output, "a") as f:
                for subf in subfolders:
                    task = f'\n"python3 run_calculate_dr.py --folder_name {folder} --subf_name {subf} --fasta {fa}"'
                    f.write(task)
        
        ## Once all tasks have been added, finish up SBATCH template
        with open(output, "a") as f:
            f.write("\n)\neval ${tasks[$SLURM_ARRAY_TASK_ID]}")
    except Exception as e:
        print("Failed to create SBATCH file: {e}")
        traceback.print_exc()
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Writes SBATCH script for deletion rate calculation.")
    parser.add_argument("--email", default = "<uniqname>@umich.edu")
    parser.add_argument("--slurm_acct", default = "<account>")
    parser.add_argument("--walltime", required = True)
    parser.add_argument("--mem", help = "Memory for job", required = True)
    parser.add_argument("--fa", help = "Directory to FASTA file of reference genome", required = True)

    args = parser.parse_args()

    print("Writing SBATCH script...")
    main(args.email, args.slurm_acct, args.walltime, args.mem, args.fa)
    print("Process finished.")