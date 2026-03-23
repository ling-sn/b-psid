#!/usr/bin/env bash
#SBATCH --job-name=clean_tsv
#SBATCH --mail-user=YOUR_UNIQNAME@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=clean_tsv_%u_%A_%a.out
#SBATCH --array=0-0
#SBATCH --account=cweidman99
#SBATCH --time=30:00
#SBATCH --mem=16000m
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
# sbatch clean_tsv.sbatch
################################################################################

module purge
eval "$(conda shell.bash hook)"
conda activate ~/miniconda3/envs/RNA-STAR

declare -a tasks=(
"python3 clean_tsv.py"
)
eval ${tasks[$SLURM_ARRAY_TASK_ID]}
