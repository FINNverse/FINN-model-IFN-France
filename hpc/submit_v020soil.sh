#!/bin/bash
# ===========================================================================
# FINN-model-IFN-France | submit the v020-soil production chain on the cluster
#   bash hpc/submit_v020soil.sh        (from the repository root)
#
# The chain of hpc/submit_v020.sh (same shape and limits) for calibration/v020-soil/
# (the v020 joint fit + soil_wet + soil_nutrient_poor, data fiveyear-v3-soil, FINN
# release-0.2.0 at 20b184f). Six jobs, each one self-contained R script run by
# hpc/run_v020soil.sbatch (resources only; every scientific choice is in the script):
#   fit        calibration/v020-soil/fit.R           24 h segment 1 (refuses to start over existing checkpoints)
#   continue   calibration/v020-soil/fit_continue.R  24 h segments 2-5 (afterany; each exits at once when model.pt exists)
#   report     calibration/v020-soil/validation.R    holdout report, 4 CPUs (afterok: segment 5)
# Sizing: the soil columns add two inputs to the first layer of each 20-20 env network and
# two weights per species to the linear regeneration response, so the cost per epoch is
# that of the v020 fit at batch 100 (358 steps per epoch). The recipe's stage A at batch 100
# ran 151 s/epoch at the start and ~700 s/epoch on average on a loaded cluster core (112
# epochs in 22 h); 500 epochs therefore take ~21-97 h, i.e. 1-5 segments of 24 h, which is
# why the chain has the fit + four continuations (120 h; a finished fit makes the remaining
# continuations exit at once). Laptop smoke of calibration/v020-soil (2026-11-14, one
# thread on efficiency cores beside another fit, batch 100): 515 s/epoch, max RSS 1.8 GB
# (peak footprint 3.9 GB); 16 G requested as for v020. Validation: 4 CPUs, 4 h as v020
# (the bare-ground run now has 24 sites = 6 climates x 4 soils instead of 6).
# Before submitting: the data folder prepared/fiveyear-v3-soil/{train,holdout_blocked}
# must be on the cluster (rsync from the laptop; manifest prepared/fiveyear-v3-soil.md5).
# ===========================================================================
set -euo pipefail
cd "$(dirname "$0")/.."
F=$(sbatch --parsable --job-name=zoo-fr-v020soil-fit --time=24:00:00 --mem=16G hpc/run_v020soil.sbatch calibration/v020-soil/fit.R)
C1=$(sbatch --parsable --job-name=zoo-fr-v020soil-c1 --time=24:00:00 --mem=16G --dependency=afterany:"$F" hpc/run_v020soil.sbatch calibration/v020-soil/fit_continue.R)
C2=$(sbatch --parsable --job-name=zoo-fr-v020soil-c2 --time=24:00:00 --mem=16G --dependency=afterany:"$C1" hpc/run_v020soil.sbatch calibration/v020-soil/fit_continue.R)
C3=$(sbatch --parsable --job-name=zoo-fr-v020soil-c3 --time=24:00:00 --mem=16G --dependency=afterany:"$C2" hpc/run_v020soil.sbatch calibration/v020-soil/fit_continue.R)
C4=$(sbatch --parsable --job-name=zoo-fr-v020soil-c4 --time=24:00:00 --mem=16G --dependency=afterany:"$C3" hpc/run_v020soil.sbatch calibration/v020-soil/fit_continue.R)
V=$(sbatch --parsable --job-name=zoo-fr-v020soil-val --time=04:00:00 --mem=16G --cpus-per-task=4 --dependency=afterok:"$C4" hpc/run_v020soil.sbatch calibration/v020-soil/validation.R)
echo "submitted: fit $F -> continue $C1 -> $C2 -> $C3 -> $C4 (afterany) -> validation $V (afterok); logs in logs/v020soil_<jobid>.out|err"
