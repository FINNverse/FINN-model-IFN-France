#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-IFN-France | calibration/v02/pilot_thetafix.R
##
## WHAT   v0.2 pilot THETAFIX = pilot_baseline.R with ONE change, the
##        regeneration fixed where the audit found it broken:
##        - recruitment dispersion theta FROZEN at a method-of-moments estimate
##          per species from the observed tallied recruit counts (bounded to
##          [0.05, 10]): in v0.1 theta fell to 0.016 for 33 of 36 species while
##          the loss stayed flat (BWI lesson 7);
##        - the negative binomial scored in TALLIED-RECRUIT COUNT space,
##          finn(recruit_obs_weight = 0.0113) = the 6 m circle (ha) a recruit is
##          tallied on, so 88.4 recruits/ha is one count;
##        - the Beverton-Holt cap K FIXED per species (custom parameter, not
##          trainable) at the observed 99 % quantile of recruitment per ha and
##          year on the training plots (floor 2/ha/yr): in v0.1 K ran to 1e4 for
##          Salix caprea and to 0.1 for twelve species.
## HYPOTHESIS  recruitment ratio -> ~1, stems stop drifting, bare ground is no
##        longer filled by a few species; growth / removal skill unchanged.
## DATA   <DATA_DIR>/prepared/fiveyear-v2/pilot2000; holdout1000 (pilot_eval.R)
## RUN    Rscript calibration/v02/pilot_thetafix.R     (one thread, ~20 min)
## OUTPUT <DATA_DIR>/fits/v02-pilots/thetafix/{model.pt, loss_history.csv, regeneration_fixed.csv, settings.md}
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-0848b60"), .libPaths()))
suppressPackageStartupMessages({ library(FINN); library(data.table); library(torch) })
torch_set_num_threads(1L)
DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "frnfi")
FIT_DIR  <- file.path(DATA_DIR, "fits", "v02-pilots", "thetafix"); dir.create(FIT_DIR, showWarnings = FALSE, recursive = TRUE)

PARAMS <- list(
  prep_version = "fiveyear-v2", train_set = "pilot2000",
  patch_size_ha = pi * 15^2 / 1e4, n_patches = 1L,
  env_vars      = c("tavg_c", "prec_mm", "gest"),
  recruits_dbh  = 8.3, hidden = c(20L, 20L), dropout = 0.1, reg_floor = 1e-3,
  recruit_obs_weight = pi * 6^2 / 1e4,                  # THE change (1/3): NB on tallied recruits (0.0113 ha)
  theta_floor = 0.05, theta_cap = 10,                   # THE change (2/3): theta frozen at the method-of-moments estimate, bounded
  K_quantile = 0.99, K_floor_ha_yr = 2,                 # THE change (3/3): K fixed at the observed 99 % quantile per species
  epochs = 200L, batchsize = 200L, lr = 0.01, lr_scheduler = "none",
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
cnt <- obs[, .(species, k = reg * PARAMS$recruit_obs_weight)]                          # tallied recruits per plot and 5 yr
mm <- cnt[, .(m = mean(k), v = var(k)), by = species][order(species)]
mm[, theta := pmin(pmax(fifelse(v > m & m > 0, m^2 / (v - m), PARAMS$theta_cap), PARAMS$theta_floor), PARAMS$theta_cap)]
Kq <- obs[, .(K = max(quantile(reg, PARAMS$K_quantile) / P, PARAMS$K_floor_ha_yr)), by = species][order(species)]
reg_tab <- merge(merge(mm, Kq, by = "species"), species, by = "species")
fwrite(reg_tab, file.path(FIT_DIR, "regeneration_fixed.csv")); print(reg_tab[, .(species_name, m, theta, K)], digits = 3)

model <- finn(
  N_species = n_species, recruits_dbh = PARAMS$recruits_dbh, reg_floor = PARAMS$reg_floor, recruit_obs_weight = PARAMS$recruit_obs_weight,
  competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
  growth_process       = createProcess(~ tavg_c + prec_mm + gest, FINN::growth, hidden = PARAMS$hidden, dropout = PARAMS$dropout, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(~ tavg_c + prec_mm + gest, FINN::mortality, hidden = PARAMS$hidden, dropout = PARAMS$dropout, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(~ tavg_c + prec_mm + gest, FINN::regeneration_saturation, hidden = PARAMS$hidden, dropout = PARAMS$dropout,
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE, dispersion_parameter = reg_tab$theta,
                                       custom_parameters = list(reg_logK = list(init = log(reg_tab$K), trainable = FALSE)))
)
## FINN always registers theta as trainable; take it out of the gradient (deviation from
## "configure only through finn()/fit()": the same torch call the FINN-bwi stage-3 fits use)
invisible(model$par_theta_recruits_raw$requires_grad_(FALSE))
stopifnot(max(abs(as.numeric(model$par_theta_recruits$detach()) - reg_tab$theta)) < 1e-3)
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
writeLines(c("# pilot thetafix", sprintf("%d sites, %d epochs, %.0f s (%.1f s/epoch); baseline + frozen theta, count-space NB (w %.4f), fixed K", uniqueN(obs$siteID), PARAMS$epochs, secs, secs / PARAMS$epochs, PARAMS$recruit_obs_weight),
             sprintf("theta unchanged after the fit: %s", max(abs(as.numeric(model$par_theta_recruits$detach()) - reg_tab$theta)) < 1e-3),
             "", "```r", capture.output(str(PARAMS, give.attr = FALSE)), "```"), file.path(FIT_DIR, "settings.md"))
message(sprintf("done thetafix: %.1f s/epoch", secs / PARAMS$epochs)); print(hist[c(1, 50, 100, 150, 200)])
