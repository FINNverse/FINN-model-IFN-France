#!/bin/bash
# ===========================================================================
# FINN-model-IFN-France | run ONE calibration script on a SLURM CPU node
# sbatch calibration/run_fit.sbatch calibration/fit.R
#
# Resources only. Every scientific choice (data set, model, epochs, losses)
# is written in the R script itself. Expects the repository at
# $HOME/FINN-zoo/frnfi and the data at $HOME/FINN-zoo-data/frnfi (the
# ../../FINN-zoo-data/frnfi layout the scripts use), FINN 0.2.0 (0848b60) in
# $HOME/R/lib-finn-0848b60 and the remaining packages in $R_LIBS_USER.
# ===========================================================================
#SBATCH --job-name=finn-ifn-fr
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/fit_%j.out
#SBATCH --error=logs/fit_%j.err
set -euo pipefail
module purge 2>/dev/null || true
module load math/R/4.4.2 2>/dev/null || true
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/site-library}"
export LD_LIBRARY_PATH="$R_LIBS_USER/torch/lib:${LD_LIBRARY_PATH:-}"
export TORCH_INSTALL=0
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}" MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
cd "${SLURM_SUBMIT_DIR:-$HOME/FINN-zoo/frnfi}"
mkdir -p logs
SCRIPT="${1:?usage: sbatch calibration/run_fit.sbatch calibration/<script>.R}"
echo "[job] $SCRIPT, $OMP_NUM_THREADS threads, $(date)"
Rscript "$SCRIPT"
echo "[job] done $(date)"
