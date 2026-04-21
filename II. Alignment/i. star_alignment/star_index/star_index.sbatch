#!/usr/bin/env bash
#SBATCH --job-name=build_index
#SBATCH --mail-user=YOUR_UNIQNAME@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=build_index_%u_%A_%a.out
#SBATCH --account=cweidman99
#SBATCH --time=9:00:00
#SBATCH --mem=45000m
#SBATCH --partition=standard
#SBATCH --ntasks-per-node=2
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
# To call this script:
# sbatch star_index.sbatch
################################################################################

module purge
eval "$(conda shell.bash hook)"
conda activate ~/miniconda3/envs/B-PSID
python3 star_index.py --overhang=100
