#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-IFN-France | calibration/v02/pilot_treemort.R
##
## WHAT   v0.2 pilot TREEMORT = pilot_baseline.R with ONE change: mortality is
##        scored PER TALLIED TREE, fit(loss_aggregation = c(mortality = "tree"),
##        tree_data = tree_obs_dt). Each tallied tree enters with its diameter
##        at the first visit, its virtual copies as `trees` (light), and its
##        removal 0/1 over the five years; FINN compares the compounded interval
##        survival with the fate.
## HYPOTHESIS  The site x species target annualises complete removals of
##        small cohorts to 100 %/yr (pooled 6.0 %/yr against 2.5 %/yr counted
##        per tree; audit section 3); per tree the level and the size gradient
##        (v0.1 under-predicted removal of trees above 40 cm by half) come right.
## DATA   <DATA_DIR>/prepared/fiveyear-v2/pilot2000 (incl. tree_obs_dt.csv); holdout1000
## RUN    Rscript calibration/v02/pilot_treemort.R     (one thread, ~25 min)
## OUTPUT <DATA_DIR>/fits/v02-pilots/treemort/{model.pt, loss_history.csv, settings.md}
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-0848b60"), .libPaths()))
suppressPackageStartupMessages({ library(FINN); library(data.table); library(torch) })
torch_set_num_threads(1L)
DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "frnfi")
FIT_DIR  <- file.path(DATA_DIR, "fits", "v02-pilots", "treemort"); dir.create(FIT_DIR, showWarnings = FALSE, recursive = TRUE)

PARAMS <- list(
  prep_version = "fiveyear-v2", train_set = "pilot2000",
  patch_size_ha = pi * 15^2 / 1e4, n_patches = 1L,
  env_vars      = c("tavg_c", "prec_mm", "gest"),
  recruits_dbh  = 8.3, hidden = c(20L, 20L), dropout = 0.1, reg_floor = 1e-3,
  reg_K_init    = 50,
  loss_aggregation = c(mortality = "tree"),             # THE change
  epochs = 200L, batchsize = 200L, lr = 0.01, lr_scheduler = "none",
  weights = "auto",
  loss_family = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality = "binomial", regeneration = "nbinom"),
  checkpoints = 10L, device = "cpu", seed = 1L
)
d <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$train_set)
obs_all <- fread(file.path(d, "obs_dt.csv")); env_all <- fread(file.path(d, "env_dt.csv")); init <- fread(file.path(d, "init_cohort_dt.csv")); species <- fread(file.path(d, "species_dt.csv"))[order(species)]
tree_obs <- fread(file.path(d, "tree_obs_dt.csv"))[, .(siteID, year, patchID, species, dbh, trees, growth, died, period_length)]
n_species <- nrow(species)
obs <- obs_all[, .(siteID, year, species, ba, dbh, trees, growth, growth_n, mort, n_at_risk, n_died, reg, period_length)]
env <- env_all[, c("siteID", "year", PARAMS$env_vars), with = FALSE]
init_cohort <- makeInitCohorts(init[, .(siteID, patchID, species, dbh, trees, treeName)], Nspecies = n_species)
message(sprintf("tree table: %d tallied trees, %d removed", nrow(tree_obs), sum(tree_obs$died)))

model <- finn(
  N_species = n_species, recruits_dbh = PARAMS$recruits_dbh, reg_floor = PARAMS$reg_floor,
  competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
  growth_process       = createProcess(~ tavg_c + prec_mm + gest, FINN::growth, hidden = PARAMS$hidden, dropout = PARAMS$dropout, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(~ tavg_c + prec_mm + gest, FINN::mortality, hidden = PARAMS$hidden, dropout = PARAMS$dropout, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(~ tavg_c + prec_mm + gest, FINN::regeneration_saturation, hidden = PARAMS$hidden, dropout = PARAMS$dropout,
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE, custom_parameters = list(reg_logK = rep(log(PARAMS$reg_K_init), n_species)))
)
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
th <- as.numeric(model$par_theta_recruits$detach()); K <- exp(as.numeric(model$reg_logK$detach()))
writeLines(c("# pilot treemort", sprintf("%d sites, %d epochs, %.0f s (%.1f s/epoch); baseline + per-tree mortality loss", uniqueN(obs$siteID), PARAMS$epochs, secs, secs / PARAMS$epochs),
             sprintf("final theta: median %.3f (min %.3f, max %.3f); K/ha/yr: median %.1f (min %.2f, max %.0f)", median(th), min(th), max(th), median(K), min(K), max(K)),
             "", "```r", capture.output(str(PARAMS, give.attr = FALSE)), "```"), file.path(FIT_DIR, "settings.md"))
message(sprintf("done treemort: %.1f s/epoch", secs / PARAMS$epochs)); print(hist[c(1, 50, 100, 150, 200)])
