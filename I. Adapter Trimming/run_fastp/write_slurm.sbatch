#!/usr/bin/env bash
#SBATCH --job-name=write_slurm
#SBATCH --mail-user=lingsn@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=write_slurm_%u_%A_%a.out
#SBATCH --account=cweidman99
#SBATCH --time=1:00
#SBATCH --mem=200m
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
# This requires the RNA-SEQ conda environment.
# 
# To call this script:
# sbatch write_slurm.sbatch
################################################################################

module purge
eval "$(conda shell.bash hook)"
conda activate ~/miniconda3/envs/RNA-SEQ

python3 write_slurm.py --input_folder raw_fastqs --output_folder trimmed_reads --email lingsn@umich.edu \
    --slurm_acct cweidman99 --walltime 1:00:00 --mem 10000
