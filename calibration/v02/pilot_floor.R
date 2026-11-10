#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-IFN-France | calibration/v02/pilot_floor.R
##
## WHAT   v0.2 pilot FLOOR: the expected best candidate (pilot_thetafix.R's
##        regeneration + pilot_treemort.R's per-tree mortality loss) with the
##        simulation's mortality floor removed:
##          createProcess(..., FINN::mortality, mortality_floor = 1e-6, mortality_draw = "poisson")
##        LIBRARY: ~/R/lib-finn-zoo-v02 FIRST on .libPaths() (FINN branch
##        zoo-v02-mortfloor, FINNverse/FINN#51 / PR #52), then lib-finn-0848b60.
##        Released FINN 0.2.0 clamps the annual death draw at >= 1 %/yr
##        (finn-class.R, torch_clamp(m, 0.01, 0.99)); the fitted rate and the
##        per-tree loss are not clamped, every simulation is.
## WHY    On the v0.1 holdout 62 % of the tallied trees had a fitted annual
##        removal rate below 1 %/yr (beech 94 %, Quercus pubescens / ilex 100 %),
##        so the simulated stand loses more trees than the fitted rate says and
##        the states push the rate parameters (audit section 1b).
## HYPOTHESIS  Stem numbers and basal area at year 5 improve where removal is
##        low; the per-tree AUC is unchanged; bare ground no longer thins the
##        low-mortality species artificially.
## DATA   <DATA_DIR>/prepared/fiveyear-v2/pilot2000; holdout1000
## RUN    Rscript calibration/v02/pilot_floor.R     (one thread, ~25 min)
## OUTPUT <DATA_DIR>/fits/v02-pilots/floor/{model.pt, loss_history.csv, regeneration_fixed.csv, settings.md}
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-zoo-v02"), path.expand("~/R/lib-finn-0848b60"), .libPaths()))   # branch build with mortality_floor FIRST
suppressPackageStartupMessages({ library(FINN); library(data.table); library(torch) })
if (!"mortality_floor" %in% names(formals(FINN::createProcess))) stop("this FINN has no createProcess(mortality_floor=): install the zoo-v02-mortfloor build in ~/R/lib-finn-zoo-v02")
torch_set_num_threads(1L)
DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "frnfi")
FIT_DIR  <- file.path(DATA_DIR, "fits", "v02-pilots", "floor"); dir.create(FIT_DIR, showWarnings = FALSE, recursive = TRUE)

PARAMS <- list(
  prep_version = "fiveyear-v2", train_set = "pilot2000",
  patch_size_ha = pi * 15^2 / 1e4, n_patches = 1L,
  env_vars      = c("tavg_c", "prec_mm", "gest"),
  recruits_dbh  = 8.3, hidden = c(20L, 20L), dropout = 0.1, reg_floor = 1e-3,
  recruit_obs_weight = pi * 6^2 / 1e4, theta_floor = 0.05, theta_cap = 10, K_quantile = 0.99, K_floor_ha_yr = 2,   # as pilot_thetafix
  loss_aggregation = c(mortality = "tree"),                                                                        # as pilot_treemort
  mortality_floor = 1e-6, mortality_draw = "poisson",                                                              # THE change
  epochs = 200L, batchsize = 200L, lr = 0.01, lr_scheduler = "none",
  weights = "auto",
  loss_family = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality = "binomial", regeneration = "nbinom"),
  checkpoints = 10L, device = "cpu", seed = 1L
)
P <- 5L
d <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$train_set)
obs_all <- fread(file.path(d, "obs_dt.csv")); env_all <- fread(file.path(d, "env_dt.csv")); init <- fread(file.path(d, "init_cohort_dt.csv")); species <- fread(file.path(d, "species_dt.csv"))[order(species)]
tree_obs <- fread(file.path(d, "tree_obs_dt.csv"))[, .(siteID, year, patchID, species, dbh, trees, growth, died, period_length)]
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
  growth_process       = createProcess(~ tavg_c + prec_mm + gest, FINN::growth, hidden = PARAMS$hidden, dropout = PARAMS$dropout, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(~ tavg_c + prec_mm + gest, FINN::mortality, hidden = PARAMS$hidden, dropout = PARAMS$dropout, optimizeSpecies = TRUE, optimizeEnv = TRUE,
                                       mortality_floor = PARAMS$mortality_floor, mortality_draw = PARAMS$mortality_draw),
  regeneration_process = createProcess(~ tavg_c + prec_mm + gest, FINN::regeneration_saturation, hidden = PARAMS$hidden, dropout = PARAMS$dropout,
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE, dispersion_parameter = reg_tab$theta,
                                       custom_parameters = list(reg_logK = list(init = log(reg_tab$K), trainable = FALSE)))
)
invisible(model$par_theta_recruits_raw$requires_grad_(FALSE))
optimizer <- if (exists("optim_ignite_adam", envir = asNamespace("torch"))) torch::optim_ignite_adam else torch::optim_adam
FINN.seed(PARAMS$seed); t0 <- Sys.time()
fit(model, env = env, data = obs, init_cohort = init_cohort, tree_data = tree_obs, loss_aggregation = PARAMS$loss_aggregation,
    epochs = PARAMS$epochs, batchsize = PARAMS$batchsize, lr = PARAMS$lr, lr_scheduler = PARAMS$lr_scheduler,
    weights = PARAMS$weights, loss_family = PARAMS$loss_family,
    patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, optimizer = optimizer, device = PARAMS$device,
    folder = NULL, checkpoints = PARAMS$checkpoints, plot_progress = FALSE)
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
torch_save(model, file.path(FIT_DIR, "model.pt"))
hist <- as.data.table(do.call(rbind, model$history)); setnames(hist, paste0("V", 1:7), c("dbh", "ba", "trees", "growth", "mortality", "regeneration", "unused"), skip_absent = TRUE); hist[, epoch := .I]
fwrite(hist, file.path(FIT_DIR, "loss_history.csv"))
writeLines(c("# pilot floor", sprintf("%d sites, %d epochs, %.0f s (%.1f s/epoch); thetafix + treemort + mortality_floor %g, draw %s; FINN %s from %s",
                                       uniqueN(obs$siteID), PARAMS$epochs, secs, secs / PARAMS$epochs, PARAMS$mortality_floor, PARAMS$mortality_draw, packageVersion("FINN"), system.file(package = "FINN")),
             "", "```r", capture.output(str(PARAMS, give.attr = FALSE)), "```"), file.path(FIT_DIR, "settings.md"))
message(sprintf("done floor: %.1f s/epoch", secs / PARAMS$epochs)); print(hist[c(1, 50, 100, 150, 200)])
