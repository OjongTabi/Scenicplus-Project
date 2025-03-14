#!/bin/bash
#SBATCH --job-name=filter_fragments_Metacluster_Paired
#SBATCH --cpus-per-task=48
#SBATCH --mem=256G
#SBATCH --time=7-00:00:00
#SBATCH --output=filter_fragments_Metacluster_Paired_%j.out
#SBATCH --error=filter_fragments_Metacluster_Paired_%j.err
#SBATCH --chdir="/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/slurm_jobs/"

# Load R module
module load R

Rscript /n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/prepare_inputs/filter_fragments_Metacluster_Paired.R
