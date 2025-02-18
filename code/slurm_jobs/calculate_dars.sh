#!/bin/bash
#SBATCH --job-name=calculate_dars
#SBATCH --cpus-per-task=48
#SBATCH --mem=256G
#SBATCH --time=21-00:00:00
#SBATCH --output=calculate_dars_%j.out
#SBATCH --error=calculate_dars_%j.err

# Load R module
module load R

Rscript /n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/prepare_inputs/calculate_dars_P24_P48_Adult.R &&

Rscript /n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/prepare_inputs/calculate_dars_P24_P48.R
