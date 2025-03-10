#!/bin/bash
#SBATCH --cpus-per-task=48
#SBATCH --mem=256G
#SBATCH --time=7-00:00:00
#SBATCH --job-name="process_four_approaches_M1toM42"
#SBATCH --output="process_four_approaches_M1toM42_%j.out"
#SBATCH --error="process_four_approaches_M1toM42_%j.err"
#SBATCH --chdir="/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/slurm_jobs/"

source ~/miniforge3/etc/profile.d/conda.sh
conda activate scenicplus_v1.0a2
export PYTHONHASHSEED=42

python3 /n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/tools/process_four_approaches_M1toM42.py