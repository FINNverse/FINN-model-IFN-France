#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-IFN-France | calibration/v02/pilot_baseline.R
##
## WHAT   v0.2 pilot BASELINE: the v0.1 fitting recipe on the fiveyear-v2 data
##        (species set S of 30, no "other", total-removal target, management
##        trace `gest` as covariate): hidden env networks on tavg_c + prec_mm +
##        gest, free recruitment dispersion, trainable Beverton-Holt cap K
##        initialised at 50 stems/ha/yr, NB per hectare, constant learning rate,
##        site x species losses. Every other pilot changes ONE thing against this
##        script.
## DATA   <DATA_DIR>/prepared/fiveyear-v2/pilot2000 (fit); holdout1000 (pilot_eval.R)
## RUN    Rscript calibration/v02/pilot_baseline.R     (one thread, ~20 min)
## OUTPUT <DATA_DIR>/fits/v02-pilots/baseline/{model.pt, loss_history.csv, settings.md}
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-0848b60"), .libPaths()))   # FINN 0.2.0 (0848b60)
suppressPackageStartupMessages({ library(FINN); library(data.table); library(torch) })
torch_set_num_threads(1L)
DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "frnfi")
FIT_DIR  <- file.path(DATA_DIR, "fits", "v02-pilots", "baseline"); dir.create(FIT_DIR, showWarnings = FALSE, recursive = TRUE)

PARAMS <- list(
  prep_version = "fiveyear-v2", train_set = "pilot2000",
  patch_size_ha = pi * 15^2 / 1e4, n_patches = 1L,
  env_vars      = c("tavg_c", "prec_mm", "gest"),
  recruits_dbh  = 8.3, hidden = c(20L, 20L), dropout = 0.1, reg_floor = 1e-3,
  reg_K_init    = 50,                                   # trainable cap, init 50 stems/ha/yr (v0.1)
  epochs = 200L, batchsize = 200L, lr = 0.01, lr_scheduler = "none",
  weights = "auto",
  loss_family = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality = "binomial", regeneration = "nbinom"),
  checkpoints = 10L, device = "cpu", seed = 1L
)
d <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$train_set)
obs_all <- fread(file.path(d, "obs_dt.csv")); env_all <- fread(file.path(d, "env_dt.csv")); init <- fread(file.path(d, "init_cohort_dt.csv")); species <- fread(file.path(d, "species_dt.csv"))[order(species)]
n_species <- nrow(species)
obs <- obs_all[, .(siteID, year, species, ba, dbh, trees, growth, growth_n, mort, n_at_risk, n_died, reg, period_length)]
env <- env_all[, c("siteID", "year", PARAMS$env_vars), with = FALSE]
init_cohort <- makeInitCohorts(init[, .(siteID, patchID, species, dbh, trees, treeName)], Nspecies = n_species)

model <- finn(
  N_species = n_species, recruits_dbh = PARAMS$recruits_dbh, reg_floor = PARAMS$reg_floor,
  competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
  growth_process       = createProcess(~ tavg_c + prec_mm + gest, FINN::growth, hidden = PARAMS$hidden, dropout = PARAMS$dropout, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(~ tavg_c + prec_mm + gest, FINN::mortality, hidden = PARAMS$hidden, dropout = PARAMS$dropout, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(~ tavg_c + prec_mm + gest, FINN::regeneration_saturation, hidden = PARAMS$hidden, dropout = PARAMS$dropout,
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE, custom_parameters = list(reg_logK = rep(log(PARAMS$reg_K_init), n_species)))
)
## torch >= 0.14 has the faster ignite Adam (FINN's default); older torch (the cluster library) falls back to optim_adam
optimizer <- if (exists("optim_ignite_adam", envir = asNamespace("torch"))) torch::optim_ignite_adam else torch::optim_adam
FINN.seed(PARAMS$seed); t0 <- Sys.time()
fit(model, env = env, data = obs, init_cohort = init_cohort,
    epochs = PARAMS$epochs, batchsize = PARAMS$batchsize, lr = PARAMS$lr, lr_scheduler = PARAMS$lr_scheduler,
    weights = PARAMS$weights, loss_family = PARAMS$loss_family,
    patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, optimizer = optimizer, device = PARAMS$device,
    folder = NULL, checkpoints = PARAMS$checkpoints, plot_progress = FALSE)
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
torch_save(model, file.path(FIT_DIR, "model.pt"))
hist <- as.data.table(do.call(rbind, model$history)); setnames(hist, paste0("V", 1:7), c("dbh", "ba", "trees", "growth", "mortality", "regeneration", "unused"), skip_absent = TRUE); hist[, epoch := .I]
fwrite(hist, file.path(FIT_DIR, "loss_history.csv"))
th <- as.numeric(model$par_theta_recruits$detach()); K <- exp(as.numeric(model$reg_logK$detach()))
writeLines(c("# pilot baseline", sprintf("%d sites, %d epochs, %.0f s (%.1f s/epoch); v0.1 recipe on fiveyear-v2 + gest", uniqueN(obs$siteID), PARAMS$epochs, secs, secs / PARAMS$epochs),
             sprintf("final theta: median %.3f (min %.3f, max %.3f); K/ha/yr: median %.1f (min %.2f, max %.0f)", median(th), min(th), max(th), median(K), min(K), max(K)),
             "", "```r", capture.output(str(PARAMS, give.attr = FALSE)), "```"), file.path(FIT_DIR, "settings.md"))
message(sprintf("done baseline: %.1f s/epoch; theta median %.3f", secs / PARAMS$epochs, median(th))); print(hist[c(1, 50, 100, 150, 200)])
