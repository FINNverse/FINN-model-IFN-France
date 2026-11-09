# FINNmodelIFNFrance (development)

* Calibration pipeline for the first release (variant "process"): IGN raw data
  2005-2025 (export 2026-10-07), prepared data `fiveyear-v1` (58,731 plots with
  one 5-year interval, 35 species + other, WorldClim 2.1 `tavg_c` and `prec_mm`),
  `fit.R` on 43,964 training plots in 66 of 82 100-km blocks, `validation.R` on
  the 14,767 holdout plots. FINN 0.2.0 (main 0848b60).
* `build_model()` rebuilds the model from the released named list of parameter
  tensors (incl. the env standardisation constants); `inst/model/species.csv`
  and `covariates.csv` describe the inputs.
* No fitted release yet: the production fit is running; `inst/model/weights.json`
  is filled by `calibration/fit.R` (tag `v0.1.0` planned) and the holdout metrics
  will be listed here with the release.
