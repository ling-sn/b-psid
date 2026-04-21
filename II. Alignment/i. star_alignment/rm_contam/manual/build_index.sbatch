#!/usr/bin/env bash
#SBATCH --job-name=build_index
#SBATCH --mail-user=YOUR_UNIQNAME@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/dev/null
#SBATCH --account=cweidman99
#SBATCH --time=1:00
#SBATCH --mem=1000m
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
# This requires a conda environment with samtools and bowtie2 (B-PSID), run create_env.sbatch 
# before running this file
# 
# To call this script:
# sbatch build_index.sbatch
################################################################################

module purge
eval "$(conda shell.bash hook)"
conda activate ~/miniconda3/envs/B-PSID
python3 build_index.py
