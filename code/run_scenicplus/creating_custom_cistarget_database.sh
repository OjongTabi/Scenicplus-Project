#!/bin/bash

# ---------------------------------------------------------------------
# Provide arguments in the following order when running the script:
# 1. Path to regions (first argument)
# 2. Database prefix (second argument)
# 3. Output file path (third argument)
# ---------------------------------------------------------------------

# Add current directory to PATH for easy access to local executables
export PATH=$PWD:$PATH

source ~/miniforge3/etc/profile.d/conda.sh 
conda activate scenicplus_v1.0a2 
export PYTHONHASHSEED=42 

# Set paths and variables
PATH_TO_REGIONS=$1
GENOME_FASTA="/n/sci/SCI-004375-NYUDATA/working/10XDrosophila-ARC/fasta/genome.fa"
CHROMSIZES="/n/sci/SCI-004375-NYUDATA/Ojong/Archive/code/scenicplus_train_M21_8_20_2024/create_cisTarget_databases/genome.chrom.sizes"
DATABASE_PREFIX=$2
SCRIPT_DIR="/n/sci/SCI-004375-NYUDATA/Ojong/Archive/code/create_cisTarget_databases"
CBDIR="/n/sci/SCI-004375-NYUDATA/Ojong/Archive/code/aertslab_motif_colleciton/v10nr_clust_public/singletons"
OUT_DIR=$3
FASTA_FILE="${OUT_DIR}/dm6.${DATABASE_PREFIX}.with_1kb_bg_padding.fa"
MOTIF_LIST="${OUT_DIR}/motifs.txt"
CBUST_DIR="/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/tools"

module --ignore_cache load scenicplus/v1.0a2
module --ignore_cache load bedtools/2.30.0

# Step 1: Create FASTA file with padded background
"${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh" \
    "${GENOME_FASTA}" \
    "${CHROMSIZES}" \
    "${PATH_TO_REGIONS}" \
    "${FASTA_FILE}" \
    1000 \
    yes 

# Step 2: Verify FASTA file

head -n 10 "${FASTA_FILE}"

# Step 3: Load motif collection
ls "${CBDIR}" > "${MOTIF_LIST}" 

# Step 4: Prepare Cluster-Buster binary (cbust)

# Add cbust directory to PATH
export PATH=$CBUST_DIR:$PATH

# Step 5: Run cistarget motif database creation
python "${SCRIPT_DIR}/create_cistarget_motif_databases.py" \
    -f "${FASTA_FILE}" \
    -M "${CBDIR}" \
    -m "${MOTIF_LIST}" \
    -o "${OUT_DIR}/${DATABASE_PREFIX}" \
    --bgpadding 1000 \
    -c "$CBUST_DIR/cbust" \
    -t 10
    -s 666 

