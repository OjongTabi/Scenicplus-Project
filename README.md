SCENIC+ Project: Rewiring Gene Regulatory Networks in Neural Circuits

This repository showcases my work utilizing SCENIC+ to infer gene regulatory networks (GRNs) within neural circuits. By integrating single-cell RNA sequencing (scRNA-seq) and ATAC sequencing (scATAC-seq) data, the project aims to understand and potentially rewire these networks to influence neural circuit function.

Repository Structure

code/: Contains scripts and workflows for data processing and analysis.

data/: Includes raw and processed datasets used in the study.

results/: Stores output files, figures, and findings from the analyses.


Key Highlights

Data Integration: Merging scRNA-seq and scATAC-seq data to construct comprehensive GRNs.

Regulatory Analysis: Identifying critical transcription factors and their target genes within neural circuits.

Network Modulation: Exploring strategies to rewire GRNs to alter neural circuit behavior.

For more details or inquiries, please contact me at ojongtabi@ymail.com.




SCENIC+ Pipeline Quick Guide
============================

To run SCENIC+ with each of the four different approaches, a separate .sh script is provided.
Each script handles all SCENIC+ steps: pycisTopic, pycisTarget, and the core SCENIC+ analysis and takes as input the selected metacluster(s) to be processed.

Run Scripts:
------------

1. Multiome P24-P48:
   /n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/slurm_jobs/run_scenicplus_metacluster_P24_P48.sh

2. Multiome P24-P48-Adult:
   /n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/slurm_jobs/run_scenicplus_metacluster.sh

3. Comprehensive P24-P48:
   /n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/slurm_jobs/run_scenicplus_metacluster_comprehensive_P24_P48.sh

4. Comprehensive P24-P48-Adult:
   /n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/slurm_jobs/run_scenicplus_metacluster_comprehensive.sh


SCENIC+ Input Files:
--------------------
SCENIC+ run requires the following four inputs:
- Cell metadata 
- ATAC fragments file 
- Differentially Accessible Regions 
- Gene expression `.h5ad` file

These inputs were generated using the scripts below :

Preparation Scripts:
--------------------

Located in:  
/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/prepare_inputs/


- Cell metadata:
  └── fetch_cell_metadata.R  
      (/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/slurm_jobs/fetch_cell_metadata.sh)

- ATAC fragments file:
  └── filter_fragments.R  
      (/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/slurm_jobs/filter_fragments.sh)

- Differentially Accessible Regions (DARs):
  ├── calculate_dars_P24_P48.R  
  │   (Multiome – P24-P48)
  └── calculate_dars_P24_P48_Adult.R  
      (Multiome – P24-P48-Adult)

- Gene expression `.h5ad` file:
  ├── convert_rds_to_h5ad.R  
  │   (Multiome – P24-P48)
  ├── convert_rds_to_h5ad_Comprehensive_RNA_P24_P48.R  
  │   (Comprehensive – P24-P48)
  ├── convert_rds_to_h5ad_Comprehensive_RNA.R  
     (Comprehensive – P24-P48-Adult)
  


Each of these prepares a subset of the required files, depending on whether the run is multiome or comprehensive, and which stages (P24, P48, Adult) are included.

