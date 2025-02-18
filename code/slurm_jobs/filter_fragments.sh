#!/bin/bash
#SBATCH --job-name=filter_fragments
#SBATCH --cpus-per-task=48
#SBATCH --mem=256G
#SBATCH --time=7-00:00:00
#SBATCH --output=filter_fragments_%j.out
#SBATCH --error=filter_fragments_%j.err

# Load R module
module load R

Rscript /n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/prepare_inputs/filter_fragments.R
