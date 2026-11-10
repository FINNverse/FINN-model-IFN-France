First release of FINNmodelIFNFrance: a fitted FINN forest model for France, calibrated to the 5-year tree-by-tree remeasurement of the IGN national forest inventory (first visits 2010-2020) and packaged as an R simulation tool.

Variant `process`: 35 species + "other", annual step on one patch of 0.0707 ha (the IGN plot), FINN process model (growth, mortality, saturating regeneration, competition) with two-layer SELU climate networks on WorldClim 2.1 `tavg_c` and `prec_mm`; natural dynamics only. Fitted with FINN 0.2.0 (commit 0848b60), torch 0.15.1, one thread, 150 epochs, baseline-scaled losses (`weights = "auto"`).

Install and use:

```r
remotes::install_github("FINNverse/FINN@0848b60")
remotes::install_github("FINNverse/FINN-model-IFN-France@v0.1.0")
library(FINNmodelIFNFrance)
torch::torch_set_num_threads(1L)
m <- load_model()   # downloads parameters.pt once, checks its sha256, rebuilds the model with the installed FINN
```

Simulations are pinned to one torch thread and `batchsize = 100`; a different batchsize changes the random draws but not the model.

Assets:

| file | what | sha256 |
|---|---|---|
| parameters.pt | named list of 28 tensors (env networks and env standardisation included); `load_model()` rebuilds from it | f278626cc0a3d5625a8e6437d98d974d4150311bc105e7e042717391351a0218 |
| model.pt | the whole fitted module for `torch::torch_load()` under FINN 0848b60 | e3c52b4d830118ab3fce3511c2f68cf931adc98420f7cdca665af6c93991db05 |

Holdout skill (14,767 plots in 100-km blocks the fit never saw; the persistence of the first visit in brackets): stand basal area R2 0.94 (0.92), stand stems R2 0.78 (0.91; about 170 stems/ha too many), species-level basal area R2 0.98, growth R2 0.38, mortality AUC 0.66 (over-predicted), recruitment R2 below zero (about four times too many recruits).

Limitations: the states are usable, the rates are weak; bare-ground succession does not settle at a plausible basal area, so use the model for short projections from an inventory; FINN's simulated mortality has a floor of 0.01 per annual step (FINNverse/FINN#51), which the fit did not see. Details in NEWS.md and `calibration/report.md`.

Data: `fiveyear-v1`, rebuilt from the IGN raw-data export by `calibration/` (provenance and checksums in `calibration/PROVENANCE.md`). Licences: code and parameters GPL (>= 3); IGN raw data Licence Ouverte / Etalab 2.0 (citation "IGN - Inventaire forestier national francais, Donnees brutes, Campagnes annuelles 2005 et suivantes, https://inventaire-forestier.ign.fr/dataIFN/"), not redistributed here; WorldClim 2.1 CC BY-SA 4.0.
