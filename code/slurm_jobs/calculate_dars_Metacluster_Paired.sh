#!/bin/bash
#SBATCH --job-name=calculate_dars_Metacluster_Paired
#SBATCH --cpus-per-task=48
#SBATCH --mem=256G
#SBATCH --time=21-00:00:00
#SBATCH --output=calculate_dars_Metacluster_Paired_%j.out
#SBATCH --error=calculate_dars_Metacluster_Paired_%j.err
#SBATCH --chdir="/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/slurm_jobs"

# Load R module
module load R

Rscript /n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/prepare_inputs/calculate_dars_P24_P48_Adult_Metacluster_Paired.R &&

Rscript /n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/prepare_inputs/calculate_dars_P24_P48_Metacluster_Paired.R
