# FINNmodelIFNFrance 0.1.0

First fitted release (tag `v0.1.0`), variant `process`.

* Data: IGN raw data 2005-2025 (export 2026-10-07, Licence Ouverte / Etalab
  2.0), prepared data `fiveyear-v1` (58,731 plots with one 5-year interval,
  first visits 2010-2020, 35 species + other, WorldClim 2.1 `tavg_c` and
  `prec_mm`), 43,964 training plots in 66 of 82 100-km blocks, 14,767 holdout
  plots in the other blocks.
* Model: FINN process model (growth, mortality, saturating regeneration,
  competition), annual timestep, each process with a two-layer SELU network
  (20, 20) on the two climate normals and per-species mechanistic parameters;
  one patch of 0.0707 ha (the IGN plot). Natural dynamics only (no harvest).
* Fit: `calibration/fit.R`, 150 epochs, 1.63 h on one core, FINN 0.2.0 (main
  0848b60), torch 0.15.1, `weights = "auto"` (baseline-scaled losses; no
  numeric weight vector, so `weight_scale` is moot).
* Release assets (GitHub release `v0.1.0`): `parameters.pt` (named list of 28
  tensors incl. env networks and env standardisation, what `load_model()`
  rebuilds from) and `model.pt` (the whole fitted module for
  `torch::torch_load()` under the same FINN version); sha256 of both in
  `inst/model/weights.json`. `build_model(parameters.pt)` reproduces the fitted
  model bit for bit. Pinned for simulation: one torch thread and
  `batchsize = 100`; a different batchsize changes the random draws but not the
  model.
* Holdout (14,767 plots in blocks never seen, `calibration/report.md`): stand
  basal area R2 0.94 (persistence of the first-visit stand 0.92), stand stems
  R2 0.78 (persistence 0.91; the model adds about 170 stems/ha too many),
  species-level basal area R2 0.98, species-mean dbh R2 0.85, stems R2 0.89,
  growth R2 0.38 (r 0.63), mortality AUC 0.66, recruitment R2 negative (about
  four times too many recruits, mostly in "other", *Abies alba* and the
  maples).

Known limitations:

* The states are usable (basal area at the level of persistence or slightly
  above), the rates are weak: recruitment is over-predicted about fourfold and
  not ranked (r 0.03); mortality is over-predicted. The report's site-level
  comparison (5-year death share 7.7 % predicted vs 2.1 % observed) overstates
  this: counted per tallied tree, natural mortality is over-predicted about
  1.5x (5-year 5.6 % vs 3.6 %), because the site-level target annualised
  complete removals of one-tree cohorts.
* Bare-ground succession does not settle at a plausible basal area (89-158
  m2/ha after 1000 years, still drifting); use the model for short projections
  from an inventory, not for succession.
* FINN's simulated mortality has a floor of 0.01 per model step
  (FINNverse/FINN#51), here 1 % per year, which the fit did not see; where the
  fitted annual death probability lies below 1 % (several oaks and beech on
  many plots, judging by the holdout rates) the simulation runs at the floor.
* Natural dynamics only: plots with cut trees entered the calibration only
  through growth, mortality and recruitment, not through their end state; the
  "other" class (176 minor species lumped) carries most of the recruitment
  error.
* The calibration range is metropolitan France; outside it the climate
  networks extrapolate without constraint.

Licences: code and fitted parameters GPL (>= 3). IGN raw data: Licence Ouverte
/ Open Licence Etalab 2.0, citation "IGN - Inventaire forestier national
francais, Donnees brutes, Campagnes annuelles 2005 et suivantes,
https://inventaire-forestier.ign.fr/dataIFN/"; this repository does not
redistribute them, `calibration/` downloads them from the source. WorldClim
2.1: CC BY-SA 4.0 (Fick & Hijmans 2017, doi:10.1002/joc.5086). Details in
`calibration/PROVENANCE.md`.

# FINNmodelIFNFrance (development)

* Repository created from the FINN model template.
