#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-IFN-France | calibration/v02/pilot_shade.R
##
## WHAT   v0.2 pilot SHADE = pilot_baseline.R with ONE change, option (b) of the
##        agreed rare-species rule: the out-of-S trees of the kept plots enter
##        as ONE extra class (species index 31, "shade class") that shades but
##        claims no demography: its observations are NA (never in the loss) and
##        its recruitment cap K is fixed at 1e-6 stems/ha/yr, so it cannot
##        recruit in fitting or simulation. Data: obs_dt_shade.csv and
##        init_cohort_shade_dt.csv of fiveyear-v2.
## WHAT THE INTERFACE DOES NOT ALLOW (deviation, disclosed): createProcess()
##        freezes a process for ALL species (optimizeSpecies / optimizeEnv), not
##        for one class, so the shade class's growth and mortality parameters
##        stay trainable; they receive gradient only indirectly (through the
##        light it casts on S), start at FINN's defaults and are reported with
##        the pilot. A per-species freeze would need a FINN change (proposal in
##        the notes if this arm wins).
## HYPOTHESIS  Keeping the out-of-S basal area as competitors improves the
##        year-5 state and growth on plots with 1-2 % out-of-S share without
##        touching recruitment; if not, option (a) (drop) stays.
## DATA   <DATA_DIR>/prepared/fiveyear-v2/pilot2000 (obs_dt_shade.csv, init_cohort_shade_dt.csv); holdout1000
## RUN    Rscript calibration/v02/pilot_shade.R     (one thread, ~20 min)
## OUTPUT <DATA_DIR>/fits/v02-pilots/shade/{model.pt, loss_history.csv, settings.md}
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-0848b60"), .libPaths()))
suppressPackageStartupMessages({ library(FINN); library(data.table); library(torch) })
torch_set_num_threads(1L)
DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "frnfi")
FIT_DIR  <- file.path(DATA_DIR, "fits", "v02-pilots", "shade"); dir.create(FIT_DIR, showWarnings = FALSE, recursive = TRUE)

PARAMS <- list(
  prep_version = "fiveyear-v2", train_set = "pilot2000",
  patch_size_ha = pi * 15^2 / 1e4, n_patches = 1L,
  env_vars      = c("tavg_c", "prec_mm", "gest"),
  recruits_dbh  = 8.3, hidden = c(20L, 20L), dropout = 0.1, reg_floor = 1e-3,
  reg_K_init    = 50, shade_K = 1e-6,                   # THE change: class 31 = out-of-S trees, NA observations, no recruitment
  epochs = 200L, batchsize = 200L, lr = 0.01, lr_scheduler = "none",
  weights = "auto",
  loss_family = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality = "binomial", regeneration = "nbinom"),
  checkpoints = 10L, device = "cpu", seed = 1L
)
d <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$train_set)
obs_all <- fread(file.path(d, "obs_dt_shade.csv")); env_all <- fread(file.path(d, "env_dt.csv")); init <- fread(file.path(d, "init_cohort_shade_dt.csv")); species <- fread(file.path(d, "species_dt.csv"))[order(species)]
n_species <- nrow(species) + 1L                        # S + the shade class
shade_id <- n_species
stopifnot(max(obs_all$species) == shade_id, max(init$species) == shade_id)
obs <- obs_all[, .(siteID, year, species, ba, dbh, trees, growth, growth_n, mort, n_at_risk, n_died, reg, period_length)]
env <- env_all[, c("siteID", "year", PARAMS$env_vars), with = FALSE]
init_cohort <- makeInitCohorts(init[, .(siteID, patchID, species, dbh, trees, treeName)], Nspecies = n_species)
message(sprintf("shade class: %d cohorts, %.1f %% of the initial basal area", init[species == shade_id, .N], 100 * init[, sum(dbh2ba(dbh) * trees * (species == shade_id)) / sum(dbh2ba(dbh) * trees)]))

K_init <- c(rep(PARAMS$reg_K_init, n_species - 1L), PARAMS$shade_K)
model <- finn(
  N_species = n_species, recruits_dbh = PARAMS$recruits_dbh, reg_floor = PARAMS$reg_floor,
  competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
  growth_process       = createProcess(~ tavg_c + prec_mm + gest, FINN::growth, hidden = PARAMS$hidden, dropout = PARAMS$dropout, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(~ tavg_c + prec_mm + gest, FINN::mortality, hidden = PARAMS$hidden, dropout = PARAMS$dropout, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(~ tavg_c + prec_mm + gest, FINN::regeneration_saturation, hidden = PARAMS$hidden, dropout = PARAMS$dropout,
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE,
                                       custom_parameters = list(reg_logK = list(init = log(K_init), trainable = FALSE)))   # K fixed: S at 50 (as the baseline init), shade at ~0
)
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
th <- as.numeric(model$par_theta_recruits$detach())
pg <- as.matrix(model$par_growth$detach()); pm <- as.matrix(model$par_mortality$detach())
writeLines(c("# pilot shade", sprintf("%d sites, %d epochs, %.0f s (%.1f s/epoch); baseline + shade class (index %d, NA obs, K %g), K of S fixed at %g (so the one change is the class)", uniqueN(obs$siteID), PARAMS$epochs, secs, secs / PARAMS$epochs, shade_id, PARAMS$shade_K, PARAMS$reg_K_init),
             sprintf("shade class parameters after the fit (indirectly identified): growth %s; mortality %s; S medians growth %s, mortality %s",
                     paste(signif(pg[shade_id, ], 3), collapse = "/"), paste(signif(pm[shade_id, ], 3), collapse = "/"),
                     paste(signif(apply(pg[-shade_id, , drop = FALSE], 2, median), 3), collapse = "/"), paste(signif(apply(pm[-shade_id, , drop = FALSE], 2, median), 3), collapse = "/")),
             "", "```r", capture.output(str(PARAMS, give.attr = FALSE)), "```"), file.path(FIT_DIR, "settings.md"))
message(sprintf("done shade: %.1f s/epoch", secs / PARAMS$epochs)); print(hist[c(1, 50, 100, 150, 200)])
