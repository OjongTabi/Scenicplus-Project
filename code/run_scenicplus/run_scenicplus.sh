#!/bin/bash

# ---------------------------------------------------------------------
# Provide arguments in the following order when running the script:
# 1. Output file path from step 2 (pycisTarget)
# 2. GEX object path
# 3. Database prefix
# ---------------------------------------------------------------------

# Set error handling and output display
set -e
echo "Starting SCENIC+ pipeline setup and data preprocessing..."

# Define paths
WORK_DIR=$1
SCPLUS_PIPELINE_DIR="${WORK_DIR}/scplus_pipeline"
cisTopic_path="${WORK_DIR}/cistopic_obj.pkl"
RAW_DATA_PATH=$2
PROCESSED_DATA_PATH="${WORK_DIR}/processed_RNA.h5ad"
NORMALIZED_DATA_PATH="${WORK_DIR}/normalized_adata.h5ad"
CONFIG_PATH="${SCPLUS_PIPELINE_DIR}/Snakemake/config/config.yaml"
GENOME_ANNOTATION_PATH="/n/sci/SCI-004375-NYUDATA/Ojong/Archive/code/scenicplus_train_M21_9_4_2024/scplus_pipeline/Snakemake"
MOTIF_ANNOTATION_PATH="/n/sci/SCI-004375-NYUDATA/Ojong/Archive/code/aertslab_motif_colleciton/v10nr_clust_public/snapshots/motifs-v10-nr.flybase-m0.00001-o0.0_updated.tbl"
CTX_DB_PATH="${WORK_DIR}/$3.regions_vs_motifs.rankings.feather"
DEM_DB_PATH="${WORK_DIR}/$3.regions_vs_motifs.scores.feather"

# Create pipeline directory
mkdir -p "${SCPLUS_PIPELINE_DIR}"
scenicplus init_snakemake --out_dir "${SCPLUS_PIPELINE_DIR}"

# Display directory structure
tree "${SCPLUS_PIPELINE_DIR}"

# Load and preprocess the raw data using Python
python3 <<EOF
import scanpy as sc

# Load raw data
adata = sc.read('${RAW_DATA_PATH}')

# Rename index columns
if '_index' in adata.var.columns:
    adata.var.rename(columns={'_index': 'index'}, inplace=True)

if adata.raw is not None and '_index' in adata.raw.var.columns:
    adata.raw.var.rename(columns={'_index': 'index'}, inplace=True)

# Set raw data if available, or use the current data
new_adata = adata.raw.to_adata() if adata.raw is not None else adata
if adata.raw is not None and 'index' in adata.raw.var.columns:
    new_adata.var_names = adata.raw.var['index']
else:
    new_adata.var_names = adata.var_names  # Use existing variable names

new_adata.var_names.name = None

# Save processed data
new_adata.write('${PROCESSED_DATA_PATH}')

# Reload processed data
adata = sc.read('${PROCESSED_DATA_PATH}')
if adata.raw is None:
    adata.raw = adata

# Normalize and log-transform data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.write('${NORMALIZED_DATA_PATH}')
EOF

# Update configuration file
python3 <<EOF
import os
import yaml

config_path = '${CONFIG_PATH}'
with open(config_path, 'r') as file:
    config_data = yaml.safe_load(file)

config_str = yaml.dump(config_data)
updated_config_str = config_str.replace("search_space_downstream: 1000 150000", "search_space_downstream: 1000 50000")
updated_config_str = updated_config_str.replace("search_space_upstream: 1000 150000", "search_space_upstream: 1000 50000")
updated_config_str = updated_config_str.replace("min_target_genes: 10", "min_target_genes: 5")
updated_config_str = updated_config_str.replace("species: homo_sapiens", "species: drosophila_melanogaster")
updated_config_str = updated_config_str.replace("species: hsapiens", "species: dmelanogaster")
updated_config_str = updated_config_str.replace("cisTopic_obj_fname: ''", f"cisTopic_obj_fname: '${cisTopic_path}'")
updated_config_str = updated_config_str.replace("GEX_anndata_fname: ''", f"GEX_anndata_fname: '${NORMALIZED_DATA_PATH}'")
updated_config_str = updated_config_str.replace("region_set_folder: ''", f"region_set_folder: '${WORK_DIR}/region_sets/'")
updated_config_str = updated_config_str.replace("ctx_db_fname: ''", f"ctx_db_fname: '${CTX_DB_PATH}'")
updated_config_str = updated_config_str.replace("dem_db_fname: ''", f"dem_db_fname: '${DEM_DB_PATH}'")
updated_config_str = updated_config_str.replace("path_to_motif_annotations: ''", f"path_to_motif_annotations: '${MOTIF_ANNOTATION_PATH}'")
updated_config_str = updated_config_str.replace("temp_dir: ''", "temp_dir: '/tmp/'")

with open(config_path, 'w') as file:
    file.write(updated_config_str)
EOF

# Display updated configuration
cat "${CONFIG_PATH}"

# Copy required files
cp "${GENOME_ANNOTATION_PATH}/genome_annotation.tsv" "${SCPLUS_PIPELINE_DIR}/Snakemake/"
cp "${GENOME_ANNOTATION_PATH}/chromsizes.tsv" "${SCPLUS_PIPELINE_DIR}/Snakemake/"

cd "${WORK_DIR}/scplus_pipeline/Snakemake"

snakemake all --cores 20 --latency-wait 60

