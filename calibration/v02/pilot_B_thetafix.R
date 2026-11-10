#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-IFN-France | calibration/v02/pilot_B_thetafix.R
##
## WHAT   v0.2 pilot B: the v0.1 model structure (hidden env networks) on the
##        v2 data with the regeneration fixed where the audit found it broken:
##        - recruitment dispersion theta FROZEN at a method-of-moments estimate
##          per species from the observed tallied recruit counts (floored at
##          0.05, capped at 10): in v0.1 theta fell to 0.016 for 33 of 36
##          species while the loss stayed flat (BWI lesson 7);
##        - the negative binomial scored in TALLIED-RECRUIT COUNT space,
##          finn(recruit_obs_weight = 0.0113) = the 6 m circle (ha) on which a
##          recruit is tallied, so observed 88.4 recruits/ha is 1 count (v0.1
##          scored per hectare, where theta has no meaning);
##        - the Beverton-Holt cap K FIXED per species (not trainable) at the
##          observed 99th percentile of recruitment per ha and year on the
##          training plots, floored at 2/ha/yr: in v0.1 K ran to 1e4 for Salix
##          and to 0.1 for twelve species;
##        - total-removal mortality with the management trace `gest` as a
##          covariate in every process (v2 data);
##        - cosine learning-rate decay (BWI lesson 10).
## HYPOTHESIS  recruitment ratio -> ~1, stems no longer drift, bare ground
##        no longer "other"-dominated; growth / mortality skill unchanged.
## DATA   <DATA_DIR>/prepared/fiveyear-v2/pilot2000 (fit), holdout1000 (pilot_eval.R)
## RUN    Rscript calibration/v02/pilot_B_thetafix.R     (one thread, ~20 min)
## OUTPUT <DATA_DIR>/fits/v02-pilots/B_thetafix/{model.pt, loss_history.csv, settings.md}
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-0848b60"), .libPaths()))
suppressPackageStartupMessages({ library(FINN); library(data.table); library(torch) })
torch_set_num_threads(1L)
DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "frnfi")
FIT_DIR  <- file.path(DATA_DIR, "fits", "v02-pilots", "B_thetafix"); dir.create(FIT_DIR, showWarnings = FALSE, recursive = TRUE)

PARAMS <- list(
  prep_version = "fiveyear-v2", train_set = "pilot2000",
  patch_size_ha = pi * 15^2 / 1e4, n_patches = 1L,
  env_vars      = c("tavg_c", "prec_mm", "gest"),      # v2: + management trace (0/1/2)
  recruits_dbh  = 8.3, hidden = c(20L, 20L), dropout = 0.1, reg_floor = 1e-3,
  recruit_obs_weight = pi * 6^2 / 1e4,                  # 0.0113 ha: NB on tallied recruits (one count = 88.4/ha)
  theta_floor = 0.05, theta_cap = 10,                   # method-of-moments theta per species, bounded
  K_quantile = 0.99, K_floor_ha_yr = 2,                 # fixed per-species cap: observed 99 % quantile of recruits/ha/yr
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

## ---- recruitment dispersion and cap from the observed counts ------------------------------
cnt <- obs[, .(species, k = reg * PARAMS$recruit_obs_weight)]          # tallied recruits per plot and 5 yr
mm <- cnt[, .(m = mean(k), v = var(k)), by = species][order(species)]
mm[, theta := fifelse(v > m & m > 0, m^2 / (v - m), PARAMS$theta_cap)]
mm[, theta := pmin(pmax(theta, PARAMS$theta_floor), PARAMS$theta_cap)]
Kq <- obs[, .(K = max(quantile(reg, PARAMS$K_quantile) / P, PARAMS$K_floor_ha_yr)), by = species][order(species)]
reg_tab <- merge(merge(mm, Kq, by = "species"), species, by = "species")
fwrite(reg_tab, file.path(FIT_DIR, "regeneration_fixed.csv")); print(reg_tab[, .(species_name, m, theta, K)], digits = 3)

## ---- model --------------------------------------------------------------------------------
model <- finn(
  N_species = n_species, recruits_dbh = PARAMS$recruits_dbh, reg_floor = PARAMS$reg_floor,
  recruit_obs_weight = PARAMS$recruit_obs_weight,
  competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
  growth_process       = createProcess(~ tavg_c + prec_mm + gest, FINN::growth, hidden = PARAMS$hidden, dropout = PARAMS$dropout, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(~ tavg_c + prec_mm + gest, FINN::mortality, hidden = PARAMS$hidden, dropout = PARAMS$dropout, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(~ tavg_c + prec_mm + gest, FINN::regeneration_saturation, hidden = PARAMS$hidden, dropout = PARAMS$dropout,
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE, dispersion_parameter = reg_tab$theta,
                                       custom_parameters = list(reg_logK = list(init = log(reg_tab$K), trainable = FALSE)))
)
## FINN always registers theta as trainable; take it out of the gradient (deviation
## from "configure only through finn()/fit()": the same call the BWI stage-3 fits use)
invisible(model$par_theta_recruits_raw$requires_grad_(FALSE))
stopifnot(max(abs(as.numeric(model$par_theta_recruits$detach()) - reg_tab$theta)) < 1e-3)

## ---- fit ----------------------------------------------------------------------------------
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
writeLines(c("# pilot B_thetafix", sprintf("%d sites, %d epochs, %.0f s (%.1f s/epoch), theta frozen, K fixed, NB in count space (w %.4f), cosine lr, gest covariate, hidden env",
                                            uniqueN(obs$siteID), PARAMS$epochs, secs, secs / PARAMS$epochs, PARAMS$recruit_obs_weight),
             sprintf("theta unchanged after the fit: %s", max(abs(as.numeric(model$par_theta_recruits$detach()) - reg_tab$theta)) < 1e-3),
             "", "```r", capture.output(str(PARAMS, give.attr = FALSE)), "```"), file.path(FIT_DIR, "settings.md"))
message(sprintf("done B_thetafix: %.1f s/epoch", secs / PARAMS$epochs)); print(hist[c(1, 50, 100, 150, 200)])
