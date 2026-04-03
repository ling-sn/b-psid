#!/usr/bin/env bash
#SBATCH --job-name=create_env
#SBATCH --mail-user=YOUR_UNIQNAME@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=create_env_%u_%A_%a.out
#SBATCH --account=cweidman99
#SBATCH --time=5:00
#SBATCH --mem=3500m
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
# To call this script:
# sbatch create_env.sbatch
################################################################################

module purge
eval "$(conda shell.bash hook)"
conda create -n B-PSID python=3.10 samtools=1.17 star=2.7.11b pysam=0.21.0 parasail-python=1.3.4 numpy=1.26 pandas=2.0 openpyxl=3.1 scipy=1.11 matplotlib=3.7.2 seaborn=0.12.2 python_abi=3.10 pysamstats=1.1.2 bs4=4.12.3 requests=2.32.5

