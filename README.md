# FINN-model-IFN-France

A fitted [FINN](https://github.com/FINNverse/FINN) forest model calibrated to the French National Forest Inventory
(IGN, Inventaire forestier national, raw data under the Licence Ouverte / Etalab 2.0), packaged as a versioned
simulation tool (R package `FINNmodelIFNFrance`). The model card is in the
[FINN-zoo](https://finnverse.github.io/FINN-zoo/models.html).

**Current release: v0.1.0** (see NEWS.md). The fitted parameters are release assets of this repository;
`load_model()` downloads them once and checks their sha256.

```r
remotes::install_github("FINNverse/FINN@0848b60")               # FINN 0.2.0, the version the model was fitted with
remotes::install_github("FINNverse/FINN-model-IFN-France@v0.1.0")
library(FINNmodelIFNFrance)
torch::torch_set_num_threads(1L)   # pinned: the fit and the validation ran single-threaded
m = load_model()          # downloads the released parameters once, checks sha256, rebuilds with the installed FINN
species_table()           # 35 species + "other": index, Latin name, IGN ESPAR code, French name
model_info()              # release tag, files, sha256, FINN commit, patch settings, pinned threads and batchsize
vignette("quickstart", package = "FINNmodelIFNFrance")
```

Simulations are run with one torch thread and an explicit `batchsize` (100 sites per `simulateForest()` batch,
the value pinned in `model_info()`; `simulate_model()` uses it by default). A different batchsize changes the
random draws but not the model.

What the model is: a FINN process model (growth, mortality, saturating regeneration, competition) for 35 French
tree species + "other", annual timestep, two WorldClim 2.1 covariates (`tavg_c`, `prec_mm`; see
`inst/model/covariates.csv` for how to derive them), calibrated on 43,964 IGN plots with a tree-by-tree 5-year
remeasurement (first visits 2010-2020) and validated on 14,767 plots in spatially blocked holdout squares.
One patch = 0.0707 ha (the IGN plot). Natural dynamics only (no harvest).

This repository is a tool, not the code of a manuscript: it is versioned for
users (see NEWS.md) and may differ from the models in our papers. The
calibration pipeline that produced each release is in `calibration/`
(`calibration/README.md`; `run_all.R` runs it end to end; `PROVENANCE.md`
pins every raw layer).

License: GPL (>= 3) for code and fitted parameters; the calibration data keep
the license of their source (IGN raw data: Licence Ouverte / Open Licence Etalab 2.0,
citation "IGN – Inventaire forestier national français, Données brutes, Campagnes
annuelles 2005 et suivantes, https://inventaire-forestier.ign.fr/dataIFN/").
Climate normals: WorldClim 2.1 (Fick & Hijmans 2017, doi:10.1002/joc.5086), free for
academic and other non-commercial use; the rasters are not redistributed, only values
extracted at the plot locations are included.
