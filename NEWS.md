# FINNmodelIFNFrance (development, towards v0.1.0)

* First production fit of the variant "process" (`calibration/fit.R`,
  2026-11-10): IGN raw data 2005-2025 (export 2026-10-07), prepared data
  `fiveyear-v1` (58,731 plots with one 5-year interval, 35 species + other,
  WorldClim 2.1 `tavg_c` and `prec_mm`), 43,964 training plots in 66 of 82
  100-km blocks, 150 epochs, FINN 0.2.0 (main 0848b60), torch 0.15.1.
* Release assets (tag `v0.1.0`, to be published): `parameters.pt` (named list of
  28 tensors incl. env networks and env scaling, sha256
  `f278626c…1a0218`) and `model.pt` (whole module, sha256 `e3c52b4d…91db05`);
  `build_model(parameters.pt)` reproduces the fitted model bit for bit.
* Holdout skill (14,767 plots in blocks never seen, `calibration/report.md`):
  stand basal area R2 0.94 (persistence of the first-visit stand 0.92), stand
  stems R2 0.78 (persistence 0.91; the model adds ~170 stems/ha too many),
  species-level basal area R2 0.98, growth R2 0.38 (r 0.63), mortality AUC 0.66
  with the 5-year death share overpredicted (7.7 % vs 2.1 % observed),
  recruitment R2 negative (about four times too many recruits, mostly in
  "other", *Abies alba* and the maples). The states are usable, the rates are
  weak; see the model card's verdict and limitations.
* `build_model()` rebuilds the model from the released named list of parameter
  tensors, failing loudly on missing or unexpected names; `inst/model/species.csv`
  and `covariates.csv` describe the inputs; quickstart vignette.
