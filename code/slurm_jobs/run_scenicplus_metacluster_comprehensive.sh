#!/bin/bash

# Read metaclusters from command-line arguments
METACLUSTERS=("$@")

# Ensure at least one metacluster is provided
if [[ ${#METACLUSTERS[@]} -eq 0 ]]; then
    echo "Error: No metaclusters provided!"
    exit 1
fi

# Loop through each metacluster and submit a job
for metacluster in "${METACLUSTERS[@]}"; do
    echo "Processing $metacluster"

    # Define paths with the metacluster variable
    RUN_SCRIPT="/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/run_scenicplus/run_scenicplus_comprehensive_RNA.sh"
    MULTIOME_RUN_PATH="/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/${metacluster}/P24_P48_Adult/"
    COMPREHENSIVE_DATA_PATH="/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/data/GEX_h5ad/Comprehensive/P24_P48_Adult/${metacluster}_Comprehensive_RNA.h5ad"
    OUTPUT_PATH="/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/${metacluster}/P24_P48_Adult_Comprehensive/"
    metacluster=$metacluster
    sbatch <<EOF
#!/bin/bash
#SBATCH --cpus-per-task=48
#SBATCH --mem=256G
#SBATCH --time=7-00:00:00
#SBATCH --job-name="run_scenicplus_comprehensive_${metacluster}"
#SBATCH --output="run_scenicplus_comprehensive_${metacluster}_%j.out"
#SBATCH --error="run_scenicplus_comprehensive_${metacluster}_%j.err"

# Load Conda
source ~/miniforge3/etc/profile.d/conda.sh
conda activate scenicplus_v1.0a2

# Run the SCENIC+ script
$RUN_SCRIPT $MULTIOME_RUN_PATH $COMPREHENSIVE_DATA_PATH $OUTPUT_PATH $metacluster

EOF
done
