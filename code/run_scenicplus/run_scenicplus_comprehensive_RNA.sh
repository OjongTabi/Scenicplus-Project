#!/bin/bash

# ---------------------------------------------------------------------
# Provide arguments in the following order when running the script:
# 1. Output file path from multiome run
# 2. Comprehensive GEX object path
# 3. Output path for this run
# 4. Metacluster
# ---------------------------------------------------------------------

# Load Conda
source ~/miniforge3/etc/profile.d/conda.sh
conda activate scenicplus_v1.0a2

# Set error handling and output display
set -e

# Define paths
MULTIOME_RUN_PATH="$1/scplus_pipeline/Snakemake"
GEX_PATH=$2
OUT_DIR=$3
cisTopic_path="$1cistopic_obj.pkl"
PROCESSED_DATA_PATH="$OUT_DIR/processed_RNA.h5ad"
NORMALIZED_DATA_PATH="$OUT_DIR/normalized_adata.h5ad"
metacluster=$4

mkdir -p "$OUT_DIR"

# Load and preprocess the raw data using Python
python3 <<EOF
import scanpy as sc

# Load raw data
adata = sc.read('${GEX_PATH}')

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

scenicplus prepare_data prepare_GEX_ACC \
    --cisTopic_obj_fname "$cisTopic_path" \
    --GEX_anndata_fname "$NORMALIZED_DATA_PATH" \
    --out_file "$OUT_DIR/ACC_GEX.h5mu" \
    --is_not_multiome \
    --key_to_group_by AnnotatedJul24

scenicplus grn_inference TF_to_gene \
    --multiome_mudata_fname "$OUT_DIR/ACC_GEX.h5mu" \
    --tf_names "$MULTIOME_RUN_PATH/tf_names.txt" \
    --temp_dir "dir/" \
    --out_tf_to_gene_adjacencies "$OUT_DIR/tf_to_gene_adj.tsv" \
    --method GBM \
    --n_cpu 40 \
    --seed 666

scenicplus grn_inference eGRN \
    --is_extended \
    --TF_to_gene_adj_fname "$OUT_DIR/tf_to_gene_adj.tsv" \
    --region_to_gene_adj_fname "$MULTIOME_RUN_PATH/region_to_gene_adj.tsv" \
    --cistromes_fname "$MULTIOME_RUN_PATH/cistromes_extended.h5ad" \
    --ranking_db_fname "$1/$metacluster.regions_vs_motifs.rankings.feather" \
    --eRegulon_out_fname "$OUT_DIR/eRegulons_extended.tsv" \
    --temp_dir "dir/" \
    --order_regions_to_genes_by importance \
    --order_TFs_to_genes_by importance \
    --gsea_n_perm 1000 \
    --quantiles 0.85 0.90 0.95 \
    --top_n_regionTogenes_per_gene 5 10 15 \
    --top_n_regionTogenes_per_region  \
    --min_regions_per_gene 0 \
    --rho_threshold 0.05 \
    --min_target_genes 5 \
    --n_cpu 40

scenicplus grn_inference AUCell \
    --eRegulon_fname "$OUT_DIR/eRegulons_extended.tsv" \
    --multiome_mudata_fname "$MULTIOME_RUN_PATH/ACC_GEX.h5mu" \
    --aucell_out_fname "$OUT_DIR/AUCell_extended.h5mu" \
    --n_cpu 40

scenicplus grn_inference eGRN \
    --TF_to_gene_adj_fname "$OUT_DIR/tf_to_gene_adj.tsv" \
    --region_to_gene_adj_fname "$MULTIOME_RUN_PATH/region_to_gene_adj.tsv" \
    --cistromes_fname "$MULTIOME_RUN_PATH/cistromes_direct.h5ad" \
    --ranking_db_fname "$1/$metacluster.regions_vs_motifs.rankings.feather" \
    --eRegulon_out_fname "$OUT_DIR/eRegulon_direct.tsv" \
    --temp_dir "dir/" \
    --order_regions_to_genes_by importance \
    --order_TFs_to_genes_by importance \
    --gsea_n_perm 1000 \
    --quantiles 0.85 0.90 0.95 \
    --top_n_regionTogenes_per_gene 5 10 15 \
    --top_n_regionTogenes_per_region  \
    --min_regions_per_gene 0 \
    --rho_threshold 0.05 \
    --min_target_genes 5 \
    --n_cpu 40

scenicplus grn_inference AUCell \
    --eRegulon_fname "$OUT_DIR/eRegulon_direct.tsv" \
    --multiome_mudata_fname "$MULTIOME_RUN_PATH/ACC_GEX.h5mu" \
    --aucell_out_fname "$OUT_DIR/AUCell_direct.h5mu" \
    --n_cpu 40



scenicplus grn_inference create_scplus_mudata \
    --multiome_mudata_fname "$MULTIOME_RUN_PATH/ACC_GEX.h5mu" \
    --e_regulon_auc_direct_mudata_fname "$OUT_DIR/AUCell_direct.h5mu" \
    --e_regulon_auc_extended_mudata_fname "$OUT_DIR/AUCell_extended.h5mu" \
    --e_regulon_metadata_direct_fname "$OUT_DIR/eRegulon_direct.tsv" \
    --e_regulon_metadata_extended_fname "$OUT_DIR/eRegulons_extended.tsv" \
    --out_file "$OUT_DIR/scplusmdata.h5mu"
