import os
import pandas as pd
import mudata

scplus_mdata = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M28/P24_P48_Adult/scplus_pipeline/Snakemake/scplusmdata.h5mu")

RNA_Counts = pd.DataFrame(
    data=scplus_mdata.mod['scRNA_counts'].X.toarray(),
    index=scplus_mdata.mod['scRNA_counts'].obs.index,  # Rows: cells
    columns=scplus_mdata.mod['scRNA_counts'].var.index  # Columns: genes
)
unique_tfs = scplus_mdata.uns['direct_e_regulon_metadata'].TF.unique()
RNA_Counts_tfs = RNA_Counts.loc[:, RNA_Counts.columns.intersection(unique_tfs)]

annotations = scplus_mdata.mod['scATAC_counts'].obs['Annotation_celltype_stage']

unique_categories = annotations.unique()

output_dir = "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M28/P24_P48_Adult/scplus_pipeline/Snakemake"

for i in unique_categories:
    
    cell_indices = annotations[annotations == i].index
    
    RNA_Counts_tfs_cell_type = RNA_Counts_tfs.loc[cell_indices]

    file_path = os.path.join(output_dir, f"RNA_Counts_tfs_{i}.tsv")

    RNA_Counts_tfs_cell_type.to_csv(file_path, sep='\t')