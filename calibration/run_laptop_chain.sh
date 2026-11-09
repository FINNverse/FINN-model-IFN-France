#!/bin/bash
# ===========================================================================
# FINN-model-IFN-France | calibration/run_laptop_chain.sh
#
# Sequencing only (no scientific parameters): waits until the machine's
# priority fits are gone (processes matching WAIT_PATTERN), then runs the
# production fit and its validation, single-threaded (set inside the scripts).
#   nohup bash calibration/run_laptop_chain.sh > logs/laptop_chain.log 2>&1 &
# from the repository root.
# ===========================================================================
set -euo pipefail
WAIT_PATTERN="${1:-fits/2026-11-09}"
mkdir -p logs
echo "[chain] $(date -u +%FT%TZ) waiting for '$WAIT_PATTERN' processes to finish"
while pgrep -f "$WAIT_PATTERN" > /dev/null; do sleep 120; done
echo "[chain] $(date -u +%FT%TZ) starting calibration/fit.R"
Rscript calibration/fit.R > logs/fit.log 2>&1
echo "[chain] $(date -u +%FT%TZ) fit done, starting calibration/validation.R"
Rscript calibration/validation.R > logs/validation.log 2>&1
echo "[chain] $(date -u +%FT%TZ) done"
