#!/bin/bash
#SBATCH --job-name=pgMI_array
#SBATCH --output=pgMI_%A_%a.log
#SBATCH --error=pgMI_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=1-20

# ----------------------------
# Base directory and scripts
# ----------------------------
BASE_DIR="/fh/fast/berger_a/Siobhan/pgMI_analyses/010525_batching/Full_pDNA"

SCRIPTS=("01_qc_batch_pDNA_010725.R" "02_annotations_batch2_pDNA_010526.R" "03_filter_LFC_batch_pDNA_010526.R" "04_calc_gi_rep_early_pDNA_010625.R")

# List of 20 cell lines (order matters for SLURM array)
CELL_LINES=("A549" "PC9" "AALE" "CALU6" "HCC15" "HCC827" "HCC4006" "HT1197" "JHOS4" "NCIH650" \
            "NCIH1373" "NCIH1975" "NCIH2110" "NCIH2228" "NCIH3122" "OVCAR8" "PANC0203" "SALE" "SKOV3" "T24")

# Select cell line based on SLURM_ARRAY_TASK_ID
CELL_LINE=${CELL_LINES[$SLURM_ARRAY_TASK_ID-1]}
echo "Array task $SLURM_ARRAY_TASK_ID running cell line: $CELL_LINE"

# Load R
module load R  # adjust to your environment

# ----------------------------
# Run scripts sequentially
# ----------------------------
for SCRIPT in "${SCRIPTS[@]}"; do
    echo "Running $SCRIPT for $CELL_LINE"
    Rscript "$BASE_DIR/$SCRIPT" "$BASE_DIR" "$CELL_LINE"
    if [ $? -ne 0 ]; then
        echo "Error: $SCRIPT failed for $CELL_LINE. Stopping this cell line."
        exit 1
    fi
done

echo "Completed cell line: $CELL_LINE"
