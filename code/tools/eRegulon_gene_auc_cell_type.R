import os
import pandas as pd
import mudata

scplus_mdata = mudata.read("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M28/P24_P48_Adult/scplus_pipeline/Snakemake/scplusmdata.h5mu")

auc_matrix = scplus_mdata.mod['direct_gene_based_AUC'].X

auc_df = pd.DataFrame(
    data=auc_matrix.toarray() if hasattr(auc_matrix, "toarray") else auc_matrix,
    index=scplus_mdata.mod['direct_gene_based_AUC'].obs_names,
    columns=scplus_mdata.mod['direct_gene_based_AUC'].var_names
)

annotations = scplus_mdata.mod['scATAC_counts'].obs['Annotation_celltype_stage']

unique_categories = annotations.unique()

output_dir = "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M28/P24_P48_Adult/scplus_pipeline/Snakemake"

for i in unique_categories:
    
    cell_indices = annotations[annotations == i].index
    
    eRegulon_gene_auc_cell_type = RNA_Counts_tfs.loc[cell_indices]

    file_path = os.path.join(output_dir, f"eRegulon_gene_auc_{i}.tsv")

    eRegulon_gene_auc_cell_type.to_csv(file_path, sep='\t')
