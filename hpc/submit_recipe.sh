#!/bin/bash
# ===========================================================================
# FINN-model-IFN-France | submit the zoo v0.2 recipe chain on the cluster
#   bash hpc/submit_recipe.sh            (from the repository root)
#
# Three jobs, each one self-contained R script run by hpc/run_pilot.sbatch
# (resources only; every scientific choice is in the script; the wrapper puts
# the libtorch of lib-finn-0848b60 on LD_LIBRARY_PATH):
#   stage A  calibration/v02/recipe/fit_stageA.R   demography, observed recruitment
#   stage B  calibration/v02/recipe/fit_stageB.R   regeneration (afterok: stage A)
#   report   calibration/v02/recipe/validation.R   holdout report (afterok: stage B)
# Time limits from the measured cost (laptop, one thread, 2026-11-10): stage A
# 0.3 s/epoch on the 200 smoke plots beyond the fixed cost -> ~55 s/epoch on the
# 35,783 training plots, 300 epochs ~4.6 h (A) and 150 epochs ~2.3 h (B); the
# cluster's cores ran the v0.2 pilots ~3x slower than the laptop (~14 h / ~7 h),
# hence 22 h and 12 h below. Validation simulates all 12,093 holdout plots and
# 1,500 training plots over 5 years plus six 500-year bare-ground runs.
# ===========================================================================
set -euo pipefail
cd "$(dirname "$0")/.."
A=$(sbatch --parsable --job-name=zoo-fr-A --time=22:00:00 --mem=16G hpc/run_pilot.sbatch calibration/v02/recipe/fit_stageA.R)
B=$(sbatch --parsable --job-name=zoo-fr-B --time=12:00:00 --mem=16G --dependency=afterok:"$A" hpc/run_pilot.sbatch calibration/v02/recipe/fit_stageB.R)
V=$(sbatch --parsable --job-name=zoo-fr-V --time=04:00:00 --mem=16G --cpus-per-task=4 --dependency=afterok:"$B" hpc/run_pilot.sbatch calibration/v02/recipe/validation.R)
echo "submitted: stage A $A -> stage B $B (afterok) -> validation $V (afterok); logs in logs/pilot_<jobid>.out|err"
