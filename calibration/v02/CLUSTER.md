# Running the v0.2 pilots on a SLURM cluster

Everything a pilot needs, so that the coordinator can sync and submit. The
scripts resolve the data as `../../FINN-zoo-data/frnfi` relative to the
repository root, so keep this layout on the cluster:

```
~/FINN-zoo/frnfi/                      <- this repository (calibration/, hpc/, R/, inst/)
~/FINN-zoo-data/frnfi/prepared/fiveyear-v2/   <- data (see below)
~/FINN-zoo-data/frnfi/fits/v02-pilots/        <- written by the pilots
```

## Files to sync (repository)
- `calibration/v02/pilot_A_v01recipe.R`, `pilot_B_thetafix.R`, `pilot_C_linear.R`,
  `pilot_D_treeloss.R` (one self-contained fit each; one torch thread; the
  optimizer falls back to `torch::optim_adam` when the installed torch has no
  `optim_ignite_adam`)
- `calibration/v02/pilot_eval.R` (scores every pilot found in `fits/v02-pilots/`
  on holdout1000 and pilot2000; run after the pilots, one core, ~2 min per pilot)
- `hpc/run_pilot.sbatch` (resource-only wrapper: 1 CPU, 8 GB, 2 h)
- `R/`, `inst/model/species.csv` are not needed by the pilots.

## Data to sync (from `FINN-zoo-data/frnfi/prepared/fiveyear-v2/` only)
- `pilot2000/`  : obs_dt.csv, env_dt.csv, init_cohort_dt.csv, tree_obs_dt.csv, species_dt.csv, siteID_dt.csv  (~30 MB)
- `holdout1000/`: the same six files (~15 MB)
- nothing else (no raw data, no national/train sets, no `.pt`)

## Submit
```
cd ~/FINN-zoo/frnfi && mkdir -p logs
for p in A_v01recipe B_thetafix C_linear D_treeloss; do
  sbatch hpc/run_pilot.sbatch calibration/v02/pilot_${p}.R
done
# after all four finished:
sbatch hpc/run_pilot.sbatch calibration/v02/pilot_eval.R
```
Expected cost per pilot: 2000 sites, 200 epochs, ~4-5 s/epoch on one core
(laptop pilot: 1.0 s/epoch per 500 sites) -> 15-20 min; pilot D (per-tree
losses) somewhat more. Memory well below 8 GB (the 500-site smoke test used
< 2 GB).

## Results to sync back
- `~/FINN-zoo-data/frnfi/fits/v02-pilots/<name>/{model.pt, loss_history.csv, settings.md, regeneration_fixed.csv, eval_predictions.csv}`
- `~/FINN-zoo/frnfi/calibration/v02/pilot_table.csv`, `pilot_table.md`, `logs/pilot_*.out`

## Library
`~/R/lib-finn-0848b60` (FINN 0.2.0, main 0848b60) first on `.libPaths()`, as in
every script; data.table and torch from `$R_LIBS_USER`. Nothing is installed by
the scripts.
