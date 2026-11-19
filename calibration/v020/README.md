# calibration/v020: the refit with the FINN 0.2.0 release defaults ("v020-defaults")

One joint `fit()` of the French NFI model (IGN, 5-year remeasurement, data
`fiveyear-v2`, 30 species) with every FINN 0.2.0 default left at its default, as
decided for the whole zoo on 2026-11-13.

## Set-up

- **FINN**: release candidate `release-0.2.0` at b6244e5 (PR FINNverse/FINN#68) in
  `~/R/lib-finn-release-b6244e5` (laptop and cluster; the only library these scripts load).
- **Calibration**: ONE joint fit of competition, growth, mortality and regeneration
  against all six responses (dbh, ba, trees, growth, mortality, regeneration),
  `loss_weights = "auto"`, `update_step = 1`, `loss_recruit = "simulated"`.
- **Defaults in use** (none passed, read back from the built model and recorded in
  `fit_summary.md`): `regeneration_capped` (cap "all", prior "max_obs_recruits" =
  max(reg / period_length) per ha and year = 353.7, prior_lambda 0.1), theta trained
  with `theta_floor = 0.5`, `mortality_draw = "binomial"`, `regeneration_sampler =
  "refined"`, mortality limits 0.001 / 0.999, `reg_floor = 1e-3`, `trees_offset = 0.01`,
  `rb_gradient = TRUE`, `light_mask = TRUE`, competition height bound log(2.3).
- **Model structure** (as the v0.2 recipe fit, `calibration/v02/recipe/`): U-shaped
  mortality; growth and mortality env = SELU network 20-20 (dropout 0.1) on the
  WorldClim normals `tavg_c + prec_mm` plus one linear per-species weight on the IGN
  management trace `gest` (0/1/2, on its own scale) outside the network (module
  `env_dnn_linmgmt`, pure-shift check in the script); static climate; recruit
  observation operator `recruit_obs_weight` = the 6 m circle (0.0113 ha); recruits
  enter at 8.3 cm; one 0.0707 ha patch per plot; `growth_rate = "annual"`; growth and
  mortality scored per tallied tree (`tree_obs_dt`, `period_mode = "per_site"`).
- **Schedule**: 300 epochs at lr 0.01, then a cosine decay to 1e-4 over 200 epochs
  (500 in total), plain Adam (`torch::optim_adam`), batch 100 (as the recipe), `clip_norm`
  2, checkpoints (parameter snapshot + model file) every 25 epochs, seed 1; 24 h cluster
  segments chained through `fit_continue.R` four times (120 h), so the fit runs as long
  as it needs (decided on 2026-11-13).
- **Data**: `FINN-zoo-data/frnfi/prepared/fiveyear-v2/train` (35,783 plots), unchanged;
  `holdout_blocked` (12,093 plots in 100 km blocks) for the validation. The inventory
  problems of the IGN data and how the prepared data handle them are listed in the
  header of `fit.R` (rare-species rule without an "other" class, five-year interval,
  total removal incl. cut trees with GEST as the management covariate, annualisation
  trap -> per-tree mortality, ingrowth vs recruitment, concentric circles and the
  6 m recruit circle, dispersion collapse -> theta floor, fuzzed coordinates).

## Differences to the v0.2 recipe (`calibration/v02/recipe/`)

| | v0.2 recipe | v020-defaults |
|---|---|---|
| stages | A: demography with the observed recruits fed in; B: regeneration alone with frozen demography | one joint fit |
| regeneration process | `regeneration_saturation`, cap K per species started at 50, theta frozen at its per-species ML value in count space | `regeneration_capped` (0.2.0 default): one cap K shared by all species, started at and pulled towards max(reg / period_length) (prior_lambda 0.1), theta trained above the floor 0.5 |
| death draw / mortality limits | gamma draw, floor 0.01 / ceiling 0.99 | binomial draw, 0.001 / 0.999 |
| recruit draw | gaussian sampler | refined sampler |
| trees loss | offset 1e-5 | offset 0.01, `rb_gradient` |
| schedule | cosine 0.005 -> 5e-5 over 300 (A) + 150 (B) epochs | constant 0.01 for 300, cosine to 1e-4 over 200 |
| batch | 100 | 100 (the first submission used 500 for the cluster budget; on 2026-11-13 I decided to keep the recipe's 100 and let the fit run as long as it needs, five 24 h segments) |
| FINN | main 0848b60 in lib-finn-0848b60 | release-0.2.0 b6244e5 in lib-finn-release-b6244e5 |
| API | `weights`, `lr_scheduler`, `optimizer` of `fit()`; `finn(reg_floor =)` | `loss_weights`, `control = finnControl()`, `loss_recruit`; `finn(growth_rate =)` |

Everything else (data, U-shape, env module, management term, recruit operator,
per-tree losses, light mask) is the recipe's.

## Files

- `fit.R`: the production fit (segment 1); refuses to start over existing checkpoints.
- `fit_continue.R`: resumes from the newest checkpoint with the schedule continued;
  exits at once when `model.pt` exists.
- `fit_smoke.R`: `fit.R` with the schedule cut to 3 epochs and its own output folder
  (`fits/fiveyear-v2-v020-smoke`); measures s/epoch and exercises every check.
- `validation.R`: the recipe validation adapted to one stage (holdout report with the
  persistence null, change metrics, strata, per-tree AUC incl. the natural-removal AUC
  without cut trees, bare-ground succession, parameter trajectories).
- `../../hpc/run_v020.sbatch`: the resource-only wrapper (one CPU, one thread; the
  release library's torch 0.17 on the library path); `../../hpc/submit_v020.sh`: the
  chain fit -> continue -> continue -> validation.

## Running

Laptop (one thread): `Rscript calibration/v020/fit_smoke.R`, then
`Rscript calibration/v020/fit.R` and `Rscript calibration/v020/validation.R` from the
repository root. Cluster: the smoke job, then `bash hpc/submit_v020.sh`.

Outputs in `FINN-zoo-data/frnfi/fits/fiveyear-v2-v020/` (outside git): `model_init.pt`,
`model.pt`, `checkpoints/`, `loss_history.csv`, `lr_schedule.csv`, `species_parameters.csv`
(incl. theta and the U-shape), `regeneration_parameters.csv` (predicted / observed
recruitment per species, K), `management_weights.csv`, `module_checks.csv`,
`grad_norms*.csv`, `cap_trajectory.csv`, `convergence.csv`, `fit_summary.md`,
`report.md`, `validation/`.

Convergence is judged on the parameter drift between checkpoints (`convergence.csv`,
`validation/fig6_parameters.png`), never on the loss alone.

## Smoke (laptop, 2026-11-13, one thread, torch 0.15.1)

3 epochs on the full training set: 26.0 s/epoch at batch 500 (71 steps per epoch),
peak RSS 2.4 GB, every check passed, every output written, checkpoint files ~114 MB
each. Baseline-scaled losses after 3 epochs: dbh 0.10, ba 0.07, trees 0.08, growth 0.68,
mortality 0.86, regeneration 0.83. Predicted / observed recruitment on train 3.20 at the
start (K = prior) -> 0.99; theta 1 -> 0.50-0.65 (the floor binds for some species
within 3 epochs); K 353.7 -> 353.3. `validation.R` and `fit_continue.R` were exercised
on the smoke folder.
