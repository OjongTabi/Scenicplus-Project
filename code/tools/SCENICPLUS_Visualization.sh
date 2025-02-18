#!/bin/bash

# Loop over metaclusters from M1 to M42
for metacluster in {1..42}  
do
    metacluster_name="M${metacluster}"  
    echo "Running notebook for $metacluster_name"

    # Run the Jupyter notebook for the current metacluster using Papermill
    papermill /n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/tools/SCENICPLUS_Visualization_Template.ipynb /n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/${metacluster_name}/P24_P48_Adult/SCENICPLUS_Visualization_${metacluster_name}.ipynb -p metacluster $metacluster_name
    jupyter nbconvert --to html /n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/${metacluster_name}/P24_P48_Adult/SCENICPLUS_Visualization_${metacluster_name}.ipynb
    echo "Notebook for $metacluster_name completed"
done


# Loop over metaclusters from M1 to M42
for metacluster in {1..42}  
do
    metacluster_name="M${metacluster}"  
    echo "Running notebook for $metacluster_name"

    # Run the Jupyter notebook for the current metacluster using Papermill
    papermill /n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/tools/SCENICPLUS_Visualization_Template_P24_P48.ipynb /n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/${metacluster_name}/P24_P48/SCENICPLUS_Visualization_${metacluster_name}.ipynb -p metacluster $metacluster_name
    jupyter nbconvert --to html /n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/${metacluster_name}/P24_P48/SCENICPLUS_Visualization_${metacluster_name}.ipynb
    echo "Notebook for $metacluster_name completed"
done