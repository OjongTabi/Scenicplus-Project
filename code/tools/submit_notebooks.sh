#!/bin/bash
#SBATCH --cpus-per-task=48
#SBATCH --mem=256G
#SBATCH --time=21-00:00:00
#SBATCH --job-name=scenicplus_nb
#SBATCH --output=scenicplus_nb_%j.out
#SBATCH --error=scenicplus_nb_%j.err

module load papermill  # Ensure Papermill is available
module load jupyter  # Ensure Jupyter is available

# Loop over metaclusters from M1 to M42
for metacluster in {3..42}
do
    metacluster_name="M${metacluster}"
    echo "Running notebook for ${metacluster_name}"

    # Run the Jupyter notebook for the current metacluster using Papermill
    papermill /n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/tools/SCENICPLUS_Visualization_Template_Comprehensive.ipynb               /n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/${metacluster_name}/P24_P48_Adult_Comprehensive/SCENICPLUS_Visualization_${metacluster_name}.ipynb               -p metacluster ${metacluster_name}

    # Convert to HTML
    jupyter nbconvert --to html /n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/${metacluster_name}/P24_P48_Adult_Comprehensive/SCENICPLUS_Visualization_${metacluster_name}.ipynb

    echo "Notebook for ${metacluster_name} completed"
done
