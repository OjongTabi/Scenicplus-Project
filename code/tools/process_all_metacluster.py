#!/usr/bin/env python3

import os
import mudata as md
import anndata as ad
import pandas as pd
import mudata
from scenicplus.plotting.dotplot import heatmap_dotplot



os.chdir("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/")


#Read all 42 MuData files dynamically
##for i in range(1, 43):
    ##globals()[f"scplus_run{i}"] = md.read(f"/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M{i}/M{i}_FourApproaches_mudata.h5mu")  # Dynamically creates global variables
scplus_run1 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M1/M1_FourApproaches_mudata.h5mu")
scplus_run2 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M2/M2_FourApproaches_mudata.h5mu")
scplus_run3 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M3/M3_FourApproaches_mudata.h5mu")
scplus_run4 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M4/M4_FourApproaches_mudata.h5mu")
scplus_run5 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M5/M5_FourApproaches_mudata.h5mu")
scplus_run6 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M6/M6_FourApproaches_mudata.h5mu")
scplus_run7 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M7/M7_FourApproaches_mudata.h5mu")
scplus_run8 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M8/M8_FourApproaches_mudata.h5mu")
scplus_run9 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M9/M9_FourApproaches_mudata.h5mu")
scplus_run10 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M10/M10_FourApproaches_mudata.h5mu")
scplus_run11 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M11/M11_FourApproaches_mudata.h5mu")
scplus_run12 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M12/M12_FourApproaches_mudata.h5mu")
scplus_run13 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M13/M13_FourApproaches_mudata.h5mu")
scplus_run14 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M14/M14_FourApproaches_mudata.h5mu")
scplus_run15 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M15/M15_FourApproaches_mudata.h5mu")
scplus_run16 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M16/M16_FourApproaches_mudata.h5mu")
scplus_run17 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M17/M17_FourApproaches_mudata.h5mu")
scplus_run18 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M18/M18_FourApproaches_mudata.h5mu")
scplus_run19 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M19/M19_FourApproaches_mudata.h5mu")
scplus_run20 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M20/M20_FourApproaches_mudata.h5mu")
scplus_run21 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M21/M21_FourApproaches_mudata.h5mu")
scplus_run22 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M22/M22_FourApproaches_mudata.h5mu")
scplus_run23 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M23/M23_FourApproaches_mudata.h5mu")
scplus_run24 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M24/M24_FourApproaches_mudata.h5mu")
scplus_run25 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M25/M25_FourApproaches_mudata.h5mu")
scplus_run26 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M26/M26_FourApproaches_mudata.h5mu")
scplus_run27 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M27/M27_FourApproaches_mudata.h5mu")
scplus_run28 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M28/M28_FourApproaches_mudata.h5mu")
scplus_run29 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M29/M29_FourApproaches_mudata.h5mu")
scplus_run30 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M30/M30_FourApproaches_mudata.h5mu")
scplus_run31 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M31/M31_FourApproaches_mudata.h5mu")
scplus_run32 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M32/M32_FourApproaches_mudata.h5mu")
scplus_run33 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M33/M33_FourApproaches_mudata.h5mu")
scplus_run34 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M34/M34_FourApproaches_mudata.h5mu")
scplus_run35 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M35/M35_FourApproaches_mudata.h5mu")
scplus_run36 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M36/M36_FourApproaches_mudata.h5mu")
scplus_run37 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M37/M37_FourApproaches_mudata.h5mu")
scplus_run38 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M38/M38_FourApproaches_mudata.h5mu")
scplus_run39 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M39/M39_FourApproaches_mudata.h5mu")
scplus_run40 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M40/M40_FourApproaches_mudata.h5mu")
scplus_run41 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M41/M41_FourApproaches_mudata.h5mu")
scplus_run42 = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M42/M42_FourApproaches_mudata.h5mu")

#Extract modalities from each MuData object
for i in range(1, 43):
    globals()[f"modalities_{i}"] = globals()[f"scplus_run{i}"].mod

#Merge modalities across all runs
merged_modalities = {}


modalities_list = list(globals()["modalities_1"].keys())  # Get modality names from the first object

for modality in modalities_list:
    merged_modalities[modality] = ad.concat(
        [globals()[f"modalities_{i}"][modality] for i in range(1, 43)],  # Concatenate all 42 modalities
        axis=0, 
        join='outer',
        label="dataset",  
        keys=[f"M{i}" for i in range(1, 43)],  # Dynamic keys for each dataset
        index_unique="-"
    )

#Create new MuData object
new_mdata = md.MuData(merged_modalities)

#Merge `uns` attributes from all runs
new_mdata.uns = {
    key: pd.concat(
        [globals()[f"scplus_run{i}"].uns[key] for i in range(1, 43)],  
        axis=0, 
        ignore_index=True
    ) for key in globals()["scplus_run1"].uns  
}

#Save the final merged MuData object
new_mdata.write("All_Metacluster_Combined_mudata.h5mu")


heatmap_dotplot(
        scplus_mudata=new_mdata,
        size_modality="direct_region_based_AUC",
        color_modality="direct_gene_based_AUC",
        group_variable="scATAC_counts:Metacluster.idents",
        eRegulon_metadata_key="direct_e_regulon_metadata",
        color_feature_key="Gene_signature_name",
        size_feature_key="Region_signature_name",
        feature_name_key="eRegulon_name",
        sort_data_by="direct_gene_based_AUC",
        orientation="vertical",
        #group_variable_order=,
        figsize=(25.0, 25.0),
        save="HeatmapDotplot_All_Metacluster_Combined.png"
    )   
    
