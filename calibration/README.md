# Calibration pipeline

Self-contained R scripts that produce the released model, run from the
repository root with `Rscript <script>` and no arguments (house rules:
r-pipeline-conventions, finn-standalone-fit). `run_all.R` runs them in order.

| order | script | what | writes |
|---|---|---|---|
| 1 | `01_download_ign_raw.R` | IGN raw data 2005-2025 and documentation, SHA-256 pinned | `../../FINN-zoo-data/frnfi/raw/`, `PROVENANCE.md` |
| 2 | `02_prepare_fiveyear_v1.R` | FINN tables for the 5-year interval, WorldClim covariates, blocked split | `.../prepared/fiveyear-v1/{national,train,holdout_blocked,pilot}`, `PREPARATION.md`, `expected_preparation.csv`, climate rows of `PROVENANCE.md` |
| 3 | `fit.R` | the production fit (variant "process") | `.../fits/fiveyear-v1-process/{model.pt, parameters.pt}`, `loss_history.csv`, `fit_summary.md`, `../inst/model/weights.json` |
| 4 | `validation.R` | holdout report | `report.md`, `validation/*.csv`, `figures/*.png` (and copies into the zoo card folder) |
| – | `fit_pilot.R` | smoke test (500 sites, 3 epochs, one thread) and epoch timing | `.../fits/pilot-smoke/` |

Runners (resources and sequencing only): `run_laptop_chain.sh` waits for the
machine to be free, then runs fit and validation single-threaded;
`run_fit.sbatch` runs one script on a SLURM CPU node.

Data layout: raw and prepared data, `*.pt` and logs live outside git in
`../../FINN-zoo-data/frnfi/` relative to the repository root (the scripts'
`DATA_DIR`). Committed records: `PROVENANCE.md` (every raw layer with URL,
licence, date and checksum), `PREPARATION.md` (what the tables mean, counts,
PARAMS), `expected_preparation.csv` (summary of the prepared tables for
diffing a rerun), `fit_summary.md`, `loss_history.csv`, `report.md`.
