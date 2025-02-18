#!/bin/bash
#SBATCH --job-name=convert_rds_to_h5ad_comprehensive
#SBATCH --cpus-per-task=48
#SBATCH --mem=256G
#SBATCH --time=7-00:00:00
#SBATCH --output=convert_rds_to_h5ad_comprehensive_%j.out
#SBATCH --error=convert_rds_to_h5ad_comprehensive_%j.err

# Load R module
module load R

Rscript /n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/prepare_inputs/convert_rds_to_h5ad_Comprehensive_RNA.R &&

Rscript /n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/prepare_inputs/convert_rds_to_h5ad_Comprehensive_RNA_P24_P48.R