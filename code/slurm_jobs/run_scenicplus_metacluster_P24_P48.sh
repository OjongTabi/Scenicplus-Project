#!/bin/bash

# Read metaclusters from command-line arguments
metaclusters=("$@")

# Ensure at least one metacluster is provided
if [[ ${#metaclusters[@]} -eq 0 ]]; then
    echo "Error: No metaclusters provided!"
    exit 1
fi


# Loop through each metacluster and submit as a separate job
for metacluster in "${metaclusters[@]}"; do
    echo "Submitting job for $metacluster"
    # Output directory
    path_to_output="/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/${metacluster}/P24_P48"

    sbatch <<EOF
#!/bin/bash
#SBATCH --cpus-per-task=48
#SBATCH --mem=256G
#SBATCH --time=7-00:00:00
#SBATCH --job-name="run_scenicplus_${metacluster}_P24_P48"
#SBATCH --output="run_scenicplus_${metacluster}_P24_P48_%j.out"
#SBATCH --error="run_scenicplus_${metacluster}_P24_P48_%j.err"

source ~/miniforge3/etc/profile.d/conda.sh
conda activate scenicplus_v1.0a2

echo "Processing $metacluster"

echo $path_to_output

mkdir -p "$path_to_output"

# Run pycisTopic
python3 /n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/run_scenicplus/run_pycistopic_P24_P48.py \
    --path_to_fragments "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/data/Fragments/atac_fragments_${metacluster}.tsv.gz" \
    --path_to_regions "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/data/dars/P24_P48/${metacluster}/P24_P48_peaks.bed" \
    --path_to_cell "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/data/cell_metadata/${metacluster}_cell_data.tsv" \
    --path_to_output "$path_to_output"

# Copy Signac DARs
cp -r "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/data/dars/P24_P48/${metacluster}/DARs_cell_type" "$path_to_output/region_sets/"

# Run pycisTarget
/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/run_scenicplus/creating_custom_cistarget_database.sh \
    "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/data/dars/P24_P48/${metacluster}/P24_P48_peaks.bed" \
    "${metacluster}" \
    "$path_to_output"

# Run SCENIC+
/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/run_scenicplus/run_scenicplus.sh \
    "$path_to_output" \
    "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/data/GEX_h5ad/${metacluster}_RNA.h5ad" \
    "${metacluster}"

echo "Finished processing $metacluster"
EOF
done