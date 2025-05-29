#!/usr/bin/env python3

import os
import mudata as md
import anndata as ad
import pandas as pd
import mudata
from scenicplus.plotting.dotplot import heatmap_dotplot

for i in range(1, 43):  
    M_ID = f"M{i}"  

    os.chdir(f"/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/{M_ID}/")

    scplus_run1 = mudata.read(f"/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/{M_ID}/P24_P48_Adult/scplus_pipeline/Snakemake/scplusmdata.h5mu")
    scplus_run2 = mudata.read(f"/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/{M_ID}/P24_P48/scplus_pipeline/Snakemake/scplusmdata.h5mu")
    scplus_run3 = mudata.read(f"/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/{M_ID}/P24_P48_Adult_Comprehensive/scplusmdata.h5mu")
    scplus_run4 = mudata.read(f"/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/{M_ID}/P24_P48_Comprehensive/scplusmdata.h5mu")

    scplus_run1.mod["scATAC_counts"].obs["Annotation_celltype_stage_analysis_type"] = scplus_run1.mod["scATAC_counts"].obs["Annotation_celltype_stage"].astype(str) + "_P24_P48_Adult"
    scplus_run2.mod["scATAC_counts"].obs["Annotation_celltype_stage_analysis_type"] = scplus_run2.mod["scATAC_counts"].obs["Annotation_celltype_stage"].astype(str) + "_P24_P48"
    scplus_run3.mod["scATAC_counts"].obs["Annotation_celltype_stage_analysis_type"] = scplus_run3.mod["scATAC_counts"].obs["Annotation_celltype_stage"].astype(str) + "_P24_P48_Adult_Comprehensive"
    scplus_run4.mod["scATAC_counts"].obs["Annotation_celltype_stage_analysis_type"] = scplus_run4.mod["scATAC_counts"].obs["Annotation_celltype_stage"].astype(str) + "_P24_P48_Comprehensive"

    modalities_1 = scplus_run1.mod
    modalities_2 = scplus_run2.mod
    modalities_3 = scplus_run3.mod
    modalities_4 = scplus_run4.mod

    merged_modalities = {}

    for modality in modalities_1.keys():
        merged_modalities[modality] = ad.concat(
            [modalities_1[modality], modalities_2[modality], modalities_3[modality], modalities_4[modality]],
            axis=0, 
            join='outer',
            label="dataset",  
            keys=["P24_P48_Adult", "P24_P48", "P24_P48_Adult_Comprehensive", "P24_P48_Comprehensive"],  
            index_unique="-", 
    )

    new_mdata = md.MuData(merged_modalities)

    new_mdata.uns = {key: pd.concat([scplus_run1.uns[key], scplus_run2.uns[key], scplus_run3.uns[key], scplus_run4.uns[key]], axis=0, ignore_index=True) for key in scplus_run1.uns}

    new_mdata.write(f"{M_ID}_FourApproaches_mudata.h5mu")

    y = list(new_mdata.mod["scATAC_counts"].obs['Annotation_celltype_stage_analysis_type'].unique())

    # Custom sort: Place items starting with 'P' first, then 'A'
    sorted_y = sorted(y, key=lambda x: (not x.startswith("P"), x))

    heatmap_dotplot(
        scplus_mudata=new_mdata,
        size_modality="direct_region_based_AUC",
        color_modality="direct_gene_based_AUC",
        group_variable="scATAC_counts:Annotation_celltype_stage_analysis_type",
        eRegulon_metadata_key="direct_e_regulon_metadata",
        color_feature_key="Gene_signature_name",
        size_feature_key="Region_signature_name",
        feature_name_key="eRegulon_name",
        sort_data_by="direct_gene_based_AUC",
        orientation="vertical",
        group_variable_order=sorted_y,
        figsize=(15.0, 25.0),
        save=f"{M_ID}_HeatmapDotplot_FourApproaches.png"
    )
