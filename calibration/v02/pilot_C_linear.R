#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-IFN-France | calibration/v02/pilot_C_linear.R
##
## WHAT   v0.2 pilot C: pilot B (frozen theta, count-space NB, fixed K, gest,
##        cosine lr) with LINEAR environmental responses instead of the hidden
##        networks: hidden = NULL, so each process gets one linear layer over
##        the standardised covariates, with quadratic terms for the two climate
##        normals so that a species can still have an optimum
##        (~ tavg_c + I(tavg_c^2) + prec_mm + I(prec_mm^2) + gest).
## HYPOTHESIS  With two climate axes a quadratic response is enough; the linear
##        model settles (parameter change falls), is less prone to the loss
##        competition that moved the v0.1 networks without improving the loss
##        (BWI lesson 9), and transfers to the holdout blocks at least as well.
## DATA   <DATA_DIR>/prepared/fiveyear-v2/pilot2000 (fit), holdout1000 (pilot_eval.R)
## RUN    Rscript calibration/v02/pilot_C_linear.R     (one thread, ~20 min)
## OUTPUT <DATA_DIR>/fits/v02-pilots/C_linear/{model.pt, loss_history.csv, settings.md}
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-0848b60"), .libPaths()))
suppressPackageStartupMessages({ library(FINN); library(data.table); library(torch) })
torch_set_num_threads(1L)
DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "frnfi")
FIT_DIR  <- file.path(DATA_DIR, "fits", "v02-pilots", "C_linear"); dir.create(FIT_DIR, showWarnings = FALSE, recursive = TRUE)

PARAMS <- list(
  prep_version = "fiveyear-v2", train_set = "pilot2000",
  patch_size_ha = pi * 15^2 / 1e4, n_patches = 1L,
  env_vars      = c("tavg_c", "prec_mm", "gest"),
  recruits_dbh  = 8.3, hidden = NULL, reg_floor = 1e-3,   # linear env responses
  recruit_obs_weight = pi * 6^2 / 1e4, theta_floor = 0.05, theta_cap = 10, K_quantile = 0.99, K_floor_ha_yr = 2,
  epochs = 200L, batchsize = 200L, lr = 0.01, lr_scheduler = "cosine",
  weights = "auto",
  loss_family = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality = "binomial", regeneration = "nbinom"),
  checkpoints = 10L, device = "cpu", seed = 1L
)
P <- 5L
d <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$train_set)
obs_all <- fread(file.path(d, "obs_dt.csv")); env_all <- fread(file.path(d, "env_dt.csv")); init <- fread(file.path(d, "init_cohort_dt.csv")); species <- fread(file.path(d, "species_dt.csv"))[order(species)]
n_species <- nrow(species)
obs <- obs_all[, .(siteID, year, species, ba, dbh, trees, growth, growth_n, mort, n_at_risk, n_died, reg, period_length)]
env <- env_all[, c("siteID", "year", PARAMS$env_vars), with = FALSE]
init_cohort <- makeInitCohorts(init[, .(siteID, patchID, species, dbh, trees, treeName)], Nspecies = n_species)
cnt <- obs[, .(species, k = reg * PARAMS$recruit_obs_weight)]
mm <- cnt[, .(m = mean(k), v = var(k)), by = species][order(species)]
mm[, theta := pmin(pmax(fifelse(v > m & m > 0, m^2 / (v - m), PARAMS$theta_cap), PARAMS$theta_floor), PARAMS$theta_cap)]
Kq <- obs[, .(K = max(quantile(reg, PARAMS$K_quantile) / P, PARAMS$K_floor_ha_yr)), by = species][order(species)]
reg_tab <- merge(merge(mm, Kq, by = "species"), species, by = "species")
fwrite(reg_tab, file.path(FIT_DIR, "regeneration_fixed.csv"))

model <- finn(
  N_species = n_species, recruits_dbh = PARAMS$recruits_dbh, reg_floor = PARAMS$reg_floor, recruit_obs_weight = PARAMS$recruit_obs_weight,
  competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
  growth_process       = createProcess(~ tavg_c + I(tavg_c^2) + prec_mm + I(prec_mm^2) + gest, FINN::growth, hidden = NULL, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(~ tavg_c + I(tavg_c^2) + prec_mm + I(prec_mm^2) + gest, FINN::mortality, hidden = NULL, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(~ tavg_c + I(tavg_c^2) + prec_mm + I(prec_mm^2) + gest, FINN::regeneration_saturation, hidden = NULL,
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE, dispersion_parameter = reg_tab$theta,
                                       custom_parameters = list(reg_logK = list(init = log(reg_tab$K), trainable = FALSE)))
)
invisible(model$par_theta_recruits_raw$requires_grad_(FALSE))   # theta frozen (as pilot B)

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
writeLines(c("# pilot C_linear", sprintf("%d sites, %d epochs, %.0f s (%.1f s/epoch); pilot B with linear + quadratic climate responses (hidden = NULL)", uniqueN(obs$siteID), PARAMS$epochs, secs, secs / PARAMS$epochs),
             "", "```r", capture.output(str(PARAMS, give.attr = FALSE)), "```"), file.path(FIT_DIR, "settings.md"))
message(sprintf("done C_linear: %.1f s/epoch", secs / PARAMS$epochs)); print(hist[c(1, 50, 100, 150, 200)])
