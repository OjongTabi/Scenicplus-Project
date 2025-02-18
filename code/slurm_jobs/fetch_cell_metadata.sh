#!/bin/bash
#SBATCH --job-name=fetch_cell_metadata
#SBATCH --cpus-per-task=48
#SBATCH --mem=256G
#SBATCH --time=7-00:00:00
#SBATCH --output=fetch_cell_metadata_%j.out
#SBATCH --error=fetch_cell_metadata_%j.err

# Load R module
module load R

Rscript /n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/prepare_inputs/fetch_cell_metadata.R
