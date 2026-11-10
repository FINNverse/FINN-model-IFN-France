# Running the v0.2 pilots on a SLURM cluster

Everything a pilot needs, so that the coordinator can sync and submit. The
scripts resolve the data as `../../FINN-zoo-data/frnfi` relative to the
repository root; keep this layout on the cluster:

```
~/FINN-zoo/frnfi/                               <- this repository (calibration/, hpc/)
~/FINN-zoo-data/frnfi/prepared/fiveyear-v2/     <- data (below)
~/FINN-zoo-data/frnfi/fits/v02-pilots/          <- written by the pilots
```

## Files to sync (repository)
- `calibration/v02/pilot_baseline.R`, `pilot_cosine.R`, `pilot_thetafix.R`,
  `pilot_treemort.R`, `pilot_floor.R`, `pilot_shade.R` (one self-contained fit
  each, one torch thread; the optimizer falls back to `torch::optim_adam` when
  the installed torch has no `optim_ignite_adam`)
- `calibration/v02/pilot_eval.R` (scores every pilot found in `fits/v02-pilots/`
  on holdout1000 and pilot2000, stratified; one core, ~3 min per pilot)
- `hpc/run_pilot.sbatch` (resource-only wrapper: 1 CPU, 8 GB, 2 h)

## Data to sync (from `FINN-zoo-data/frnfi/prepared/fiveyear-v2/` only)
- `pilot2000/`  : obs_dt.csv, env_dt.csv, init_cohort_dt.csv, tree_obs_dt.csv,
  species_dt.csv, siteID_dt.csv, obs_dt_shade.csv, init_cohort_shade_dt.csv (~30 MB)
- `holdout1000/`: the same eight files (~15 MB)
- nothing else (no raw data, no national/train sets, no `.pt`)

## Libraries
- `~/R/lib-finn-0848b60` (FINN 0.2.0, main 0848b60) first on `.libPaths()` in
  every pilot except `pilot_floor.R`; data.table and torch from `$R_LIBS_USER`.
- `pilot_floor.R` and `pilot_eval.R` put `~/R/lib-finn-zoo-v02` FIRST (the FINN
  branch `zoo-v02-mortfloor`, `createProcess(mortality_floor=, mortality_draw=)`;
  installed by the coordinator, never by the scripts). `pilot_eval.R` only reads
  models with it (a superset of 0.2.0), so one library serves all pilots.

## Submit
```
cd ~/FINN-zoo/frnfi && mkdir -p logs
for p in baseline cosine thetafix treemort floor shade; do
  sbatch hpc/run_pilot.sbatch calibration/v02/pilot_${p}.R
done
# after all six finished:
sbatch hpc/run_pilot.sbatch calibration/v02/pilot_eval.R
```
Expected cost per pilot: 2000 sites, 200 epochs, ~4-5 s/epoch on one core
(laptop smoke test: 1.0 s/epoch per 500 sites) -> 15-20 min; `treemort` and
`floor` (per-tree loss) somewhat more. Memory well below 8 GB.

## Results to sync back
- `~/FINN-zoo-data/frnfi/fits/v02-pilots/<name>/{model.pt, loss_history.csv, settings.md, regeneration_fixed.csv, eval_predictions.csv}`
- `~/FINN-zoo/frnfi/calibration/v02/pilot_table.csv`, `pilot_table_strata.csv`, `pilot_table.md`, `logs/pilot_*.out`
