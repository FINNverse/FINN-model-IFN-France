#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-IFN-France | calibration/v02/pilot_D_treeloss.R
##
## WHAT   v0.2 pilot D: pilot B (hidden env, frozen theta, count-space NB,
##        fixed K, gest, cosine lr) with growth and mortality scored PER TALLIED
##        TREE instead of as site x species means:
##        fit(loss_aggregation = c(growth = "tree", mortality = "tree"), tree_data = ...).
##        Each tallied tree enters with its diameter at the first visit, its
##        virtual copies as `trees` (light), its annual relative growth (survivors)
##        and its removal 0/1 over the five years; FINN compares the compounded
##        interval survival with the fate, so no annualisation of a share is
##        needed and a plot with one tree counts as one tree.
## HYPOTHESIS  The size gradient of removal (v0.1 under-predicted mortality of
##        trees above 40 cm by half) and the species levels improve because the
##        loss sees every tree with its own size; recruitment as in B.
## DATA   <DATA_DIR>/prepared/fiveyear-v2/pilot2000 (fit incl. tree_obs_dt.csv), holdout1000 (pilot_eval.R)
## RUN    Rscript calibration/v02/pilot_D_treeloss.R     (one thread, ~25 min)
## OUTPUT <DATA_DIR>/fits/v02-pilots/D_treeloss/{model.pt, loss_history.csv, settings.md}
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-0848b60"), .libPaths()))
suppressPackageStartupMessages({ library(FINN); library(data.table); library(torch) })
torch_set_num_threads(1L)
DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "frnfi")
FIT_DIR  <- file.path(DATA_DIR, "fits", "v02-pilots", "D_treeloss"); dir.create(FIT_DIR, showWarnings = FALSE, recursive = TRUE)

PARAMS <- list(
  prep_version = "fiveyear-v2", train_set = "pilot2000",
  patch_size_ha = pi * 15^2 / 1e4, n_patches = 1L,
  env_vars      = c("tavg_c", "prec_mm", "gest"),
  recruits_dbh  = 8.3, hidden = c(20L, 20L), dropout = 0.1, reg_floor = 1e-3,
  recruit_obs_weight = pi * 6^2 / 1e4, theta_floor = 0.05, theta_cap = 10, K_quantile = 0.99, K_floor_ha_yr = 2,
  loss_aggregation = c(growth = "tree", mortality = "tree"),   # the one difference to pilot B
  epochs = 200L, batchsize = 200L, lr = 0.01, lr_scheduler = "cosine",
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
message(sprintf("tree table: %d tallied trees, %d removed, %d with growth", nrow(tree_obs), sum(tree_obs$died), sum(is.finite(tree_obs$growth))))
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
  mortality_process    = createProcess(~ tavg_c + prec_mm + gest, FINN::mortality, hidden = PARAMS$hidden, dropout = PARAMS$dropout, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(~ tavg_c + prec_mm + gest, FINN::regeneration_saturation, hidden = PARAMS$hidden, dropout = PARAMS$dropout,
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE, dispersion_parameter = reg_tab$theta,
                                       custom_parameters = list(reg_logK = list(init = log(reg_tab$K), trainable = FALSE)))
)
invisible(model$par_theta_recruits_raw$requires_grad_(FALSE))

## torch >= 0.14 has the faster ignite Adam (FINN's default); older torch (the cluster library) falls back to optim_adam
optimizer <- if (exists("optim_ignite_adam", envir = asNamespace("torch"))) torch::optim_ignite_adam else torch::optim_adam
FINN.seed(PARAMS$seed); t0 <- Sys.time()
fit(model, env = env, data = obs, init_cohort = init_cohort,
    tree_data = tree_obs, loss_aggregation = PARAMS$loss_aggregation,
    epochs = PARAMS$epochs, batchsize = PARAMS$batchsize, lr = PARAMS$lr, lr_scheduler = PARAMS$lr_scheduler,
    weights = PARAMS$weights, loss_family = PARAMS$loss_family,
    patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, optimizer = optimizer, device = PARAMS$device,
    folder = NULL, checkpoints = PARAMS$checkpoints, plot_progress = FALSE)
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
torch_save(model, file.path(FIT_DIR, "model.pt"))
hist <- as.data.table(do.call(rbind, model$history)); setnames(hist, paste0("V", 1:7), c("dbh", "ba", "trees", "growth", "mortality", "regeneration", "unused"), skip_absent = TRUE); hist[, epoch := .I]
fwrite(hist, file.path(FIT_DIR, "loss_history.csv"))
writeLines(c("# pilot D_treeloss", sprintf("%d sites, %d epochs, %.0f s (%.1f s/epoch); pilot B with growth and mortality scored per tallied tree", uniqueN(obs$siteID), PARAMS$epochs, secs, secs / PARAMS$epochs),
             "", "```r", capture.output(str(PARAMS, give.attr = FALSE)), "```"), file.path(FIT_DIR, "settings.md"))
message(sprintf("done D_treeloss: %.1f s/epoch", secs / PARAMS$epochs)); print(hist[c(1, 50, 100, 150, 200)])
