#!/usr/bin/env python

# ---------------------------------------------------------------------
# 1. Path to fragments (--fragments)  
# 2. Path to common peaks (--regions)  
# 3. Path to cell metadata (--cell)  
# 4. Output file path (--output)  
# ---------------------------------------------------------------------


# Import libraries for cisTopic analysis and visualization.
import os
import pycisTopic
import pandas as pd
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
import polars as pl
import pickle
from pycisTopic.clust_vis import (
    find_clusters,
    run_umap,
    run_tsne,
    plot_metadata,
    plot_topic,
    cell_topic_heatmap
)
import matplotlib.colors as mcolors
import random
from pycisTopic.lda_models import run_cgs_models_mallet
from pycisTopic.lda_models import evaluate_models
from pycisTopic.topic_binarization import binarize_topics
from pycisTopic.topic_qc import compute_topic_metrics, plot_topic_qc, topic_annotation
import matplotlib.pyplot as plt
from pycisTopic.utils import fig2img
from pycisTopic.topic_qc import compute_topic_metrics, plot_topic_qc, topic_annotation
import matplotlib.pyplot as plt
from pycisTopic.utils import fig2img
from pycisTopic.utils import region_names_to_coordinates
from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features
)
import numpy as np
from pycisTopic.lda_models import evaluate_models
import time
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--path_to_fragments")
parser.add_argument("--path_to_regions")
parser.add_argument("--path_to_cell")
parser.add_argument("--path_to_output")
args = parser.parse_args()

# Assign arguments to variables
path_to_fragments = args.path_to_fragments
path_to_regions = args.path_to_regions
path_to_cell = args.path_to_cell
out_dir = args.path_to_output


# Create a cistopic object for the current sample
cistopic_obj = create_cistopic_object_from_fragments(
    path_to_fragments=path_to_fragments,
    path_to_regions=path_to_regions,
    path_to_blacklist="/n/sci/SCI-004375-NYUDATA/Hassan/working_dir/data/input_files/blacklist.bed",
    n_cpu=1,
    project="",
    split_pattern=''
)

# Load cell metadata
cell_data = pd.read_table(path_to_cell, index_col = 0, dtype = str)
print(cell_data['AnnotatedJul24'].unique())

# Add cell metadata to the cisTopic object
cistopic_obj.add_cell_data(cell_data)
cistopic_obj.cell_data['Annotation_celltype_stage'] = cistopic_obj.cell_data['stage'] + "_" + cistopic_obj.cell_data['AnnotatedJul24']


# Ensure the output directory exists
os.makedirs(out_dir, exist_ok=True)  # Create the directory if it does not exist
os.chdir(out_dir)

# Save the cistopic object
output_file = os.path.join(out_dir, "cistopic_obj.pkl")
with open(output_file, "wb") as f:
    pickle.dump(cistopic_obj, f)

# Load the cistopic object
with open(os.path.join(out_dir, "cistopic_obj.pkl"), 'rb') as file:
    cistopic_obj = pickle.load(file)



# Download and extract Mallet
url = "https://github.com/mimno/Mallet/releases/download/v202108/Mallet-202108-bin.tar.gz"
filename = "Mallet-202108-bin.tar.gz"

os.system(f"wget {url}")

# Extract the tar.gz file
os.system(f"tar -xf {filename}")

# Create the tutorial directory
os.makedirs("tutorial", exist_ok=True)

# Set the MALLET memory allocation to 200GB.
os.environ['MALLET_MEMORY'] = '200G'

# Configure path Mallet
mallet_path='Mallet-202108/bin/mallet'


# Run models
models = run_cgs_models_mallet(
    cistopic_obj,
    n_topics=[2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 52, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 102, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150], 
    n_cpu=15,
    n_iter=500,
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    tmp_path='tutorial',
    save_path='tutorial',
    mallet_path=mallet_path,
)

pickle.dump(
    models,
    open(os.path.join(out_dir, "models.pkl"), "wb")
)

file_path = os.path.join(out_dir, 'models.pkl')
while not os.path.exists(file_path):
    time.sleep(10)  # Check every 10 seconds

# Open the file in binary read mode and load the object
with open(file_path, 'rb') as file:
    models = pickle.load(file)

# Evaluate the models and return the best-performing model.
model = evaluate_models(
    models,
    return_model = True
)

# Add the selected LDA model to the cisTopic object.
cistopic_obj.add_LDA_model(model)

# Save the cisTopic object to a pickle file
pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj.pkl"), "wb")
)

# Print the selected model from the cisTopic object.
print(cistopic_obj.selected_model)

# Perform clustering on the cisTopic object with specified parameters for resolution and scaling.
find_clusters(
    cistopic_obj,  
    target='cell',  
    k=10,  
    res=[0.6, 1.2, 3],  
    prefix='pycisTopic_',  
    scale=True,  
    split_pattern=''  
)

# Run UMAP and t-SNE dimensionality reduction on the cisTopic object for cell-level analysis with scaling.
run_umap(
    cistopic_obj,  
    target='cell', 
    scale=True  
)


run_tsne(
    cistopic_obj,
    target  = 'cell', scale=True)

# Plot metadata variables on the UMAP reduction for cell-level analysis with customized text and dot sizes.
plot_metadata(
    cistopic_obj,  
    reduction_name='UMAP',  
    variables=['Annotation_celltype_stage', 'pycisTopic_leiden_10_0.6', 'pycisTopic_leiden_10_1.2', 'pycisTopic_leiden_10_3'],  # List of metadata variables to be plotted.
    target='cell',  
    num_columns=2, 
    text_size=10,  
    dot_size=5  
)


plot_metadata(
    cistopic_obj,  
    reduction_name='UMAP', 
    variables=['Annotation_celltype_stage', 'pycisTopic_leiden_10_0.6', 'pycisTopic_leiden_10_1.2', 'pycisTopic_leiden_10_3'],  # List of metadata variables to be plotted.
    target='cell',  
    num_columns=2,  
    text_size=10,  
    dot_size=10  
)

# Create a dictionary to annotate clusters at different resolutions with the most frequent cell type and stage.
annot_dict = {}


for resolution in [0.6, 1.2, 3]:
    
    annot_dict[f"pycisTopic_leiden_10_{resolution}"] = {}
    
    
    for cluster in set(cistopic_obj.cell_data[f"pycisTopic_leiden_10_{resolution}"]):
        
        counts = cistopic_obj.cell_data.loc[
            cistopic_obj.cell_data.loc[cistopic_obj.cell_data[f"pycisTopic_leiden_10_{resolution}"] == cluster].index,
            "Annotation_celltype_stage"].value_counts()
        
        
        annot_dict[f"pycisTopic_leiden_10_{resolution}"][cluster] = f"{counts.index[counts.argmax()]}({cluster})"


annot_dict

# Update cluster annotations in the cisTopic object using the predefined annotation dictionary for each resolution.
for resolution in [0.6, 1.2, 3]:
    
    cistopic_obj.cell_data[f'pycisTopic_leiden_10_{resolution}'] = [
        annot_dict[f'pycisTopic_leiden_10_{resolution}'][x] for x in cistopic_obj.cell_data[f'pycisTopic_leiden_10_{resolution}'].tolist()
    ]

# Visualize metadata variables on the UMAP reduction for cell-level analysis with specified layout and style parameters.
plot_metadata(
    cistopic_obj,
    reduction_name='UMAP',
    variables=['Annotation_celltype_stage', 'pycisTopic_leiden_10_0.6', 'pycisTopic_leiden_10_1.2', 'pycisTopic_leiden_10_3'],
    target='cell', num_columns=2,
    text_size=10,
    dot_size=5)


# Plot topics on the UMAP reduction for cell-level visualization with a specified number of columns.
plot_topic(
    cistopic_obj,
    reduction_name = 'UMAP',
    target = 'cell',
    num_columns=3
)


# Generate a random color dictionary for unique cell type-stage annotations and plot a cell-topic heatmap with customized legend and figure settings.
unique_levels = cistopic_obj.cell_data['Annotation_celltype_stage'].unique()


color_dictionary = {'Annotation_celltype_stage': {}}


random.seed(42)  
for level in unique_levels:
    
    color = mcolors.to_rgb("#" + "%06x" % random.randint(0, 0xFFFFFF))
    color_hex = mcolors.to_hex(color)
    
    color_dictionary['Annotation_celltype_stage'][level] = color_hex


cell_topic_heatmap(
    cistopic_obj,
    variables=['Annotation_celltype_stage'],
    scale=False,
    legend_loc_x=1.0,
    legend_loc_y=-1.2,
    legend_dist_y=-1,
    figsize=(10, 10),
    color_dictionary=color_dictionary
)



# Binarize topics using the top 3,000 regions per topic and visualize the results with plots arranged in 3 columns.
##region_bin_topics_top_3k = binarize_topics(
    ##cistopic_obj,    
    ##method='ntop',   
    ##ntop=3_000,      
    ##plot=True,       
    ##num_columns=3    
##)

# Set ntop to the available regions if fewer than 3,000.
ntop = 3_000  
while ntop > 0:
    try:
        region_bin_topics_top_3k = binarize_topics(cistopic_obj, method='ntop', ntop=ntop, plot=True, num_columns=3)
        print(f"Using ntop={ntop}")
        break  
    except IndexError:
        ntop -= 1  

if ntop <= 0:
    print("Binarization failed: No valid ntop found")

# Binarize topics using the Otsu method and visualize the results with plots arranged in 5 columns.
region_bin_topics_otsu = binarize_topics(
    cistopic_obj, method='otsu',
    plot=True, num_columns=5
)

# Binarize cell-topic associations using the Li method with 100 bins and visualize the results with plots arranged in 5 columns.
binarized_cell_topic = binarize_topics(
    cistopic_obj,
    target='cell',
    method='li',
    plot=True,
    num_columns=5, nbins=100)



# Compute topic quality control (QC) metrics and create visualizations for various topic relationships.
# Store QC plots in a dictionary and arrange them in a grid layout for visualization.

# Compute topic QC metrics
topic_qc_metrics = compute_topic_metrics(cistopic_obj)


# Generate QC plots for different variable pairs
fig_dict={}
fig_dict['CoherenceVSAssignments']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Log10_Assignments', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['AssignmentsVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Log10_Assignments', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSRegions_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Regions_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSMarginal_dist']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Marginal_topic_dist', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSGini_index']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Gini_index', var_color='Gini_index', plot=False, return_fig=True)


fig=plt.figure(figsize=(40, 43))
i = 1
for fig_ in fig_dict.keys():
    plt.subplot(2, 3, i)
    img = fig2img(fig_dict[fig_])
    plt.imshow(img)
    plt.axis('off')
    i += 1
plt.subplots_adjust(wspace=0, hspace=-0.70)
plt.show()


# Annotate topics based on the 'Annotation_celltype_stage' variable and binarized cell-topic associations,
# using a general topic threshold of 0.2 for annotation.
topic_annot = topic_annotation(
    cistopic_obj,
    annot_var='Annotation_celltype_stage',
    binarized_cell_topic=binarized_cell_topic,
    general_topic_thr = 0.2
)


# Impute accessibility for the cisTopic object across all cells and regions, scaling the results by a factor of 10^6.
imputed_acc_obj = impute_accessibility(
    cistopic_obj,
    selected_cells=None,
    selected_regions=None,
    scale_factor=10**6
)

# Normalize the imputed accessibility scores with a scaling factor of 10^4.
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)

# Identify highly variable regions from the normalized imputed accessibility object using specified thresholds,
# and visualize the results with a plot.
variable_regions = find_highly_variable_features(
    normalized_imputed_acc_obj,
    min_disp = 0.05,
    min_mean = 0.0125,
    max_mean = 3,
    max_disp = np.inf,
    n_bins=20,
    n_top_features=None,
    plot=True
)


len(variable_regions)

# Identify differential features for each level of 'Annotation_celltype_stage' using imputed accessibility data
# and variable regions, applying thresholds for adjusted p-value and log2 fold change.
markers_dict = find_diff_features(
    cistopic_obj,           
    imputed_acc_obj,        
    variable='Annotation_celltype_stage', 
    var_features=variable_regions,  
    contrasts=None,         
    adjpval_thr=0.05,       
    log2fc_thr=0.5, 
    n_cpu=1,                
    _temp_dir='/scratch/leuven/330/vsc33053/ray_spill',  
    split_pattern=''       
)

# Print the number of DARs for each category, create necessary directories for region sets and DARs,
# and export region sets and DARs to BED files, organizing them by topics (Otsu and top 3k) and cell types.
print("Number of DARs found:")
print("---------------------")
for x in markers_dict:
    print(f"  {x}: {len(markers_dict[x])}")


os.makedirs(os.path.join(out_dir, "region_sets"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "Topics_otsu"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "Topics_top_3k"), exist_ok = True)
os.makedirs(os.path.join(out_dir,"DARs_cell_type"), exist_ok = True)


for topic in region_bin_topics_otsu:
    region_names_to_coordinates(
        region_bin_topics_otsu[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(out_dir, "region_sets", "Topics_otsu", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )


for topic in region_bin_topics_top_3k:
    region_names_to_coordinates(
        region_bin_topics_top_3k[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(out_dir, "region_sets", "Topics_top_3k", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )



for cell_type in markers_dict:
    region_names_to_coordinates(
        markers_dict[cell_type].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(out_dir, "DARs_cell_type", f"{cell_type}.bed"),
        sep = "\t",
        header = False, index = False
    )

