#!/bin/bash
# ===========================================================================
# FINN-model-IFN-France | submit the v020-defaults production chain on the cluster
#   bash hpc/submit_v020.sh            (from the repository root, after the smoke job)
#
# Four jobs, each one self-contained R script run by hpc/run_v020.sbatch
# (resources only; every scientific choice is in the script):
#   fit        calibration/v020/fit.R           24 h segment 1 (refuses to start over existing checkpoints)
#   continue   calibration/v020/fit_continue.R  24 h segment 2 (afterany: fit; exits at once when model.pt exists)
#   continue   calibration/v020/fit_continue.R  24 h segment 3 (afterany: segment 2; same)
#   report     calibration/v020/validation.R    holdout report, 4 CPUs (afterok: segment 3)
# Sizing from the smoke (fit_smoke.R: 26 s/epoch on a laptop core at batch 500;
# the cluster's cores ran the v0.2 pilots ~3x slower): 500 epochs ~ 11-14 h on
# the cluster, so segment 1 should finish the fit and the two continuations
# exit at once; they are the safety net for a slower node (the recipe's stage A
# timed out at 22 h). Memory: 2.4 GB peak on the laptop smoke, 16 G requested.
# Smoke first: sbatch --time=02:00:00 --mem=16G hpc/run_v020.sbatch calibration/v020/fit_smoke.R
# ===========================================================================
set -euo pipefail
cd "$(dirname "$0")/.."
F=$(sbatch --parsable --job-name=zoo-fr-v020-fit --time=24:00:00 --mem=16G hpc/run_v020.sbatch calibration/v020/fit.R)
C1=$(sbatch --parsable --job-name=zoo-fr-v020-c1 --time=24:00:00 --mem=16G --dependency=afterany:"$F" hpc/run_v020.sbatch calibration/v020/fit_continue.R)
C2=$(sbatch --parsable --job-name=zoo-fr-v020-c2 --time=24:00:00 --mem=16G --dependency=afterany:"$C1" hpc/run_v020.sbatch calibration/v020/fit_continue.R)
# batch 100 (as the recipe) can take up to ~700 s/epoch on a loaded core: two more
# segments, so the chain covers 120 h and the fit runs as long as it needs
C3=$(sbatch --parsable --job-name=zoo-fr-v020-c3 --time=24:00:00 --mem=16G --dependency=afterany:"$C2" hpc/run_v020.sbatch calibration/v020/fit_continue.R)
C4=$(sbatch --parsable --job-name=zoo-fr-v020-c4 --time=24:00:00 --mem=16G --dependency=afterany:"$C3" hpc/run_v020.sbatch calibration/v020/fit_continue.R)
V=$(sbatch --parsable --job-name=zoo-fr-v020-val --time=04:00:00 --mem=16G --cpus-per-task=4 --dependency=afterok:"$C4" hpc/run_v020.sbatch calibration/v020/validation.R)
echo "submitted: fit $F -> continue $C1 -> $C2 -> $C3 -> $C4 (afterany) -> validation $V (afterok); logs in logs/v020_<jobid>.out|err"
