#!/bin/bash
#SBATCH --job-name=P0_Dm_M24_Without_NB_NE_prepare_inputs
#SBATCH --cpus-per-task=48
#SBATCH --mem=256G
#SBATCH --time=21-00:00:00
#SBATCH --output=P0_Dm_M24_Without_NB_NE_prepare_inputs_%j.out
#SBATCH --error=P0_Dm_M24_Without_NB_NE_prepare_inputs_%j.err
#SBATCH --chdir="/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/slurm_jobs"

# Load R module
module load R

Rscript /n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/prepare_inputs/calculate_dars_P0_Dm_M24_Without_NB_NE.R &&

Rscript /n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/prepare_inputs/convert_rds_to_h5ad_P0_Dm_M24_Without_NB_NE.R &&

Rscript /n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/prepare_inputs/filter_fragments_P0_Dm_M24_Without_NB_NE.R
