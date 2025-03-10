#!/usr/bin/env python3

import pathlib
from typing import List
from scenicplus.cli.commands import (
    prepare_GEX_ACC,
    run_motif_enrichment_cistarget,
    run_motif_enrichment_dem,
    prepare_motif_enrichment_results,
    get_search_space_command,
    infer_region_to_gene,
    infer_TF_to_gene,
    infer_grn,
    calculate_auc,
    create_scplus_mudata
)

# ------------------------------------------------------------------------------
# Step 1: Prepare GEX_ACC
# ------------------------------------------------------------------------------
cisTopic_obj_fname = pathlib.Path(
    "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M1/P24_P48_Adult/cistopic_obj.pkl"
)
GEX_anndata_fname = pathlib.Path(
    "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M1/P24_P48_Adult/normalized_adata.h5ad"
)
out_file = pathlib.Path("ACC_GEX.h5mu")
bc_transform_func = lambda x: f"{x}"

prepare_GEX_ACC(
    cisTopic_obj_fname=cisTopic_obj_fname,
    GEX_anndata_fname=GEX_anndata_fname,
    out_file=out_file,
    use_raw_for_GEX_anndata=True,
    is_multiome=True,
    bc_transform_func=bc_transform_func,
    key_to_group_by=None,
    nr_metacells=None,
    nr_cells_per_metacells=10
)

# ------------------------------------------------------------------------------
# Step 2: Run Motif Enrichment (CisTarget & DEM)
# ------------------------------------------------------------------------------
region_set_folder = "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M1/P24_P48_Adult/region_sets"

cistarget_db_fname = "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M1/P24_P48_Adult/M1.regions_vs_motifs.rankings.feather"

output_fname_cistarget_result = "ctx_results.hdf5"
path_to_motif_annotations ="/n/sci/SCI-004375-NYUDATA/Ojong/Archive/code/aertslab_motif_colleciton/v10nr_clust_public/snapshots/motifs-v10-nr.flybase-m0.00001-o0.0_updated.tbl"


run_motif_enrichment_cistarget(
    region_set_folder=region_set_folder,
    cistarget_db_fname=cistarget_db_fname,
    output_fname_cistarget_result=output_fname_cistarget_result,
    n_cpu=1,
    fraction_overlap_w_cistarget_database=0.4,
    auc_threshold=0.005,
    nes_threshold=3.0,
    rank_threshold=0.05,
    path_to_motif_annotations=path_to_motif_annotations,
    annotation_version="v10nr_clust",
    motif_similarity_fdr=0.001,
    orthologous_identity_threshold=0.0,
    temp_dir="/tmp/",
    species="drosophila_melanogaster",
    annotations_to_use=["Direct_annot", "Orthology_annot"]
)

run_motif_enrichment_dem(
    region_set_folder=region_set_folder,
    dem_db_fname="/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M1/P24_P48_Adult/M1.regions_vs_motifs.scores.feather",
    output_fname_dem_result="dem_results.hdf5",
    n_cpu=1,
    temp_dir="/tmp/",
    species="drosophila_melanogaster",
    fraction_overlap_w_dem_database=0.4,
    max_bg_regions=500,
    path_to_genome_annotation="/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/data/genome_annotation.tsv",
    balance_number_of_promoters=True,
    promoter_space=1000,
    adjpval_thr=0.05,
    log2fc_thr=1.0,
    mean_fg_thr=0.0,
    motif_hit_thr=3.0,
    path_to_motif_annotations=path_to_motif_annotations,
    annotation_version="v10nr_clust",
    annotations_to_use=["Direct_annot","Orthology_annot"],
    motif_similarity_fdr=0.001,
    orthologous_identity_threshold=0.0,
    seed=666,
    write_html=True,
    output_fname_dem_html=pathlib.Path("dem_results.html")
)

# ------------------------------------------------------------------------------
# Step 3: Prepare Motif Enrichment Results
# ------------------------------------------------------------------------------
prepare_motif_enrichment_results(
    paths_to_motif_enrichment_results=["dem_results.hdf5","ctx_results.hdf5"],
    multiome_mudata_fname="ACC_GEX.h5mu",
    out_file_direct_annotation="cistromes_direct.h5ad",
    out_file_extended_annotation="cistromes_extended.h5ad",
    out_file_tf_names="tf_names.txt",
    direct_annotation=["Direct_annot"],
    extended_annotation=["Orthology_annot"]
)

# ------------------------------------------------------------------------------
# Step 4: Get Search Space
# ------------------------------------------------------------------------------
get_search_space_command(
    multiome_mudata_fname=pathlib.Path("ACC_GEX.h5mu"),
    gene_annotation_fname=pathlib.Path("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/data/genome_annotation.tsv"),
    chromsizes_fname=pathlib.Path("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/data/chromsizes.tsv"),
    out_fname=pathlib.Path("search_space.tsv"),
    use_gene_boundaries=True,
    upstream=(1000, 50000),
    downstream=(1000, 50000),
    extend_tss=(10, 10),
    remove_promoters=False
)

# ------------------------------------------------------------------------------
# Step 5: Infer Region to Gene & TF to Gene
# ------------------------------------------------------------------------------
infer_region_to_gene(
    multiome_mudata_fname=pathlib.Path("ACC_GEX.h5mu"),
    search_space_fname=pathlib.Path("search_space.tsv"),
    temp_dir=pathlib.Path("/tmp/"),
    adj_out_fname=pathlib.Path("region_to_gene_adj.tsv"),
    importance_scoring_method="GBM",
    correlation_scoring_method="SR",
    mask_expr_dropout=False,
    n_cpu=1
)

infer_TF_to_gene(
    multiome_mudata_fname=pathlib.Path("ACC_GEX.h5mu"),
    tf_names_fname=pathlib.Path("tf_names.txt"),
    temp_dir=pathlib.Path("/tmp/"),
    adj_out_fname=pathlib.Path("tf_to_gene_adj.tsv"),
    method="GBM",
    n_cpu=1,
    seed=666
)

# ------------------------------------------------------------------------------
# Step 6: Infer GRN
# ------------------------------------------------------------------------------
infer_grn(
    TF_to_gene_adj_fname="tf_to_gene_adj.tsv",
    region_to_gene_adj_fname="region_to_gene_adj.tsv",
    cistromes_fname="cistromes_direct.h5ad",
    eRegulon_out_fname=pathlib.Path("eRegulon_direct.tsv"),
    ranking_db_fname=pathlib.Path(
        "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M1/P24_P48_Adult/M1.regions_vs_motifs.scores.feather"
    ),
    is_extended=False,
    temp_dir=pathlib.Path("/tmp/"),
    order_regions_to_genes_by="importance",
    order_TFs_to_genes_by="importance",
    gsea_n_perm=1000,
    quantiles=[0.85, 0.90, 0.95],
    top_n_regionTogenes_per_gene=[5, 10, 15],
    top_n_regionTogenes_per_region=[],
    binarize_using_basc=False,
    min_regions_per_gene=0,
    rho_dichotomize_tf2g=False,
    rho_dichotomize_r2g=False,
    rho_dichotomize_eregulon=False,
    keep_only_activating=False,
    rho_threshold=0.05,
    min_target_genes=5,
    n_cpu=1,
    seed=666
)

infer_grn(
    TF_to_gene_adj_fname="tf_to_gene_adj.tsv",
    region_to_gene_adj_fname="region_to_gene_adj.tsv",
    cistromes_fname="cistromes_extended.h5ad",
    eRegulon_out_fname=pathlib.Path("eRegulon_extended.tsv"),
    ranking_db_fname=pathlib.Path(
        "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M1/P24_P48_Adult/M1.regions_vs_motifs.scores.feather"
    ),
    is_extended=True,
    temp_dir=pathlib.Path("/tmp/"),
    order_regions_to_genes_by="importance",
    order_TFs_to_genes_by="importance",
    gsea_n_perm=1000,
    quantiles=[0.85, 0.90, 0.95],
    top_n_regionTogenes_per_gene=[5, 10, 15],
    top_n_regionTogenes_per_region=[],
    binarize_using_basc=False,
    min_regions_per_gene=0,
    rho_dichotomize_tf2g=False,
    rho_dichotomize_r2g=False,
    rho_dichotomize_eregulon=False,
    keep_only_activating=False,
    rho_threshold=0.05,
    min_target_genes=5,
    n_cpu=1,
    seed=666
)

# ------------------------------------------------------------------------------
# Step 7: Calculate AUC
# ------------------------------------------------------------------------------
calculate_auc(
    eRegulons_fname="eRegulon_direct.tsv",
    multiome_mudata_fname="ACC_GEX.h5mu",
    out_file="AUCell_direct.h5mu",
    n_cpu=1
)

calculate_auc(
    eRegulons_fname="eRegulon_extended.tsv",
    multiome_mudata_fname="ACC_GEX.h5mu",
    out_file="AUCell_extended.h5mu",
    n_cpu=1
)

# ------------------------------------------------------------------------------
# Step 8: Create SCENIC+ MuData Object
# ------------------------------------------------------------------------------
create_scplus_mudata(
    multiome_mudata_fname=pathlib.Path("ACC_GEX.h5mu"),
    e_regulon_auc_direct_mudata_fname=pathlib.Path("AUCell_direct.h5mu"),
    e_regulon_auc_extended_mudata_fname=pathlib.Path("AUCell_extended.h5mu"),
    e_regulon_metadata_direct_fname=pathlib.Path("eRegulon_direct.tsv"),
    e_regulon_metadata_extended_fname=pathlib.Path("eRegulon_extended.tsv"),
    out_file=pathlib.Path("scplusmdata.h5mu")
)
