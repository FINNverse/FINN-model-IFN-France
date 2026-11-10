#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-IFN-France | calibration/v02/fit.R
##
## WHAT   The v0.2 production fit (variant "process", fits/fiveyear-v2-process):
##        FINN process model for the 30 species of S on the fiveyear-v2 plots
##        (species outside S < 2 % of the start-census basal area; no "other"),
##        total-removal mortality with the management trace `gest` as a
##        covariate, and the regeneration and mortality settings the v0.2 audit
##        asked for. PREPARED BEFORE THE PILOT TABLE: the settings below are the
##        expected best arm (thetafix + treemort + cosine); the coordinator
##        checks them against calibration/v02/pilot_table.md before starting,
##        and the mortality floor arm (pilot_floor.R) replaces the mortality
##        createProcess() if it wins.
## WHY THIS SHAPE (each item is tied to an audit finding, AUDIT_v02.md)
##        - species set S + plot rule instead of "other": the v0.1 "other" took
##          154 of 193 excess recruits per ha and every bare-ground run (A1);
##        - total removal + gest: removal is 1.8 / 2.4 / 2.7 %/yr for GEST 0/1/2;
##          censoring cuts biased natural mortality and masked 17 % of states (A3);
##        - recruits by back-modelled start diameter (A2);
##        - theta frozen at the method-of-moments estimate, NB in tallied-recruit
##          count space (recruit_obs_weight = 0.0113 ha), K fixed at the observed
##          99 % quantile: v0.1 theta collapsed to 0.016, K ran to 1e4 / 0.1 (A4);
##        - mortality scored per tallied tree: the site x species annualisation
##          turns a cut single tree into 100 %/yr (pooled 6.0 vs 2.5 %/yr) (A3);
##        - cosine lr decay: the v0.1 loss was flat from epoch 50 while the
##          parameters kept moving (A6);
##        - hidden env networks on tavg_c + prec_mm + gest as v0.1 (the pilots
##          test no alternative structure); recruits enter at 8.3 cm.
##        Losses: dbh/ba mse, stems poisson, growth mse (site x species),
##        removal binomial per tree, recruitment nbinom; weights = "auto"
##        (baseline-scaled; no numeric vector, weight_scale unused).
## DATA   <DATA_DIR>/prepared/fiveyear-v2/train (fit), holdout_blocked (validation.R)
## RUN    Rscript calibration/v02/fit.R      (repository root; one thread; ~37k sites,
##        ~35 s/epoch (v0.1: 39 s for 44k) -> ~1.5 h for 150 epochs; or
##        sbatch calibration/run_fit.sbatch calibration/v02/fit.R)
## OUTPUT <DATA_DIR>/fits/fiveyear-v2-process/{model.pt, parameters.pt, loss_history.csv}
##        calibration/v02/loss_history.csv, fit_summary.md, inst/model/weights.json
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-0848b60"), .libPaths()))   # FINN 0.2.0, main commit 0848b60
suppressPackageStartupMessages({ library(FINN); library(data.table); library(torch) })

DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "frnfi")
FIT_DIR  <- file.path(DATA_DIR, "fits", "fiveyear-v2-process")
CAL_DIR  <- file.path("calibration", "v02")
dir.create(FIT_DIR, showWarnings = FALSE, recursive = TRUE)

PARAMS <- list(
  prep_version  = "fiveyear-v2", train_set = "train", holdout_set = "holdout_blocked",
  patch_size_ha = pi * 15^2 / 1e4, n_patches = 1L,
  env_vars      = c("tavg_c", "prec_mm", "gest"),     # WorldClim 2.1 bio1, bio12; IGN management trace 0/1/2
  recruits_dbh  = 8.3, hidden = c(20L, 20L), dropout = 0.1, reg_floor = 1e-3,
  recruit_obs_weight = pi * 6^2 / 1e4,                 # NB on tallied recruits: the 6 m circle (0.0113 ha)
  theta_floor = 0.05, theta_cap = 10,                  # theta frozen at the method-of-moments estimate per species, bounded
  K_quantile = 0.99, K_floor_ha_yr = 2,                # K fixed at the observed 99 % quantile of recruits/ha/yr per species
  loss_aggregation = c(mortality = "tree"),            # removal scored per tallied tree
  epochs = 150L, batchsize = 500L, lr = 0.01, lr_scheduler = "cosine",
  weights = "auto",
  loss_family = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality = "binomial", regeneration = "nbinom"),
  checkpoints = 10L, threads = 1L, device = "cpu", seed = 1L
)
P <- 5L
torch_set_num_threads(PARAMS$threads)

## ---- 1. data ------------------------------------------------------------------------
d <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$train_set)
train <- list(obs = fread(file.path(d, "obs_dt.csv")), env = fread(file.path(d, "env_dt.csv")), init = fread(file.path(d, "init_cohort_dt.csv")),
              trees = fread(file.path(d, "tree_obs_dt.csv")), species = fread(file.path(d, "species_dt.csv"))[order(species)])
n_species <- nrow(train$species)
obs <- train$obs[, .(siteID, year, species, ba, dbh, trees, growth, growth_n, mort, n_at_risk, n_died, reg, period_length)]
env <- train$env[, c("siteID", "year", PARAMS$env_vars), with = FALSE]
tree_obs <- train$trees[, .(siteID, year, patchID, species, dbh, trees, growth, died, period_length)]
init_cohort <- makeInitCohorts(train$init[, .(siteID, patchID, species, dbh, trees, treeName)], Nspecies = n_species)
message(sprintf("train: %d sites, %d species, %d obs rows, %d cohorts, %d tallied trees (%d removed)", uniqueN(obs$siteID), n_species, nrow(obs), nrow(train$init), nrow(tree_obs), sum(tree_obs$died)))
## recruitment dispersion and cap from the observed tallied counts (frozen / fixed in the model)
cnt <- obs[, .(species, k = reg * PARAMS$recruit_obs_weight)]
mm <- cnt[, .(m = mean(k), v = var(k)), by = species][order(species)]
mm[, theta := pmin(pmax(fifelse(v > m & m > 0, m^2 / (v - m), PARAMS$theta_cap), PARAMS$theta_floor), PARAMS$theta_cap)]
Kq <- obs[, .(K = max(quantile(reg, PARAMS$K_quantile) / P, PARAMS$K_floor_ha_yr)), by = species][order(species)]
reg_tab <- merge(merge(mm, Kq, by = "species"), train$species, by = "species")
fwrite(reg_tab, file.path(CAL_DIR, "regeneration_fixed.csv"))

## ---- 2. the model, written out ----------------------------------------------------------
model <- finn(
  N_species = n_species, recruits_dbh = PARAMS$recruits_dbh, reg_floor = PARAMS$reg_floor, recruit_obs_weight = PARAMS$recruit_obs_weight,
  competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
  growth_process       = createProcess(~ tavg_c + prec_mm + gest, FINN::growth, hidden = PARAMS$hidden, dropout = PARAMS$dropout, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(~ tavg_c + prec_mm + gest, FINN::mortality, hidden = PARAMS$hidden, dropout = PARAMS$dropout, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(~ tavg_c + prec_mm + gest, FINN::regeneration_saturation, hidden = PARAMS$hidden, dropout = PARAMS$dropout,
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE, dispersion_parameter = reg_tab$theta,
                                       custom_parameters = list(reg_logK = list(init = log(reg_tab$K), trainable = FALSE)))
)
invisible(model$par_theta_recruits_raw$requires_grad_(FALSE))   # theta frozen (deviation: a torch call on a public parameter, as FINN-bwi stage 3)
stopifnot(max(abs(as.numeric(model$par_theta_recruits$detach()) - reg_tab$theta)) < 1e-3)

## ---- 3. fit -------------------------------------------------------------------------------
optimizer <- if (exists("optim_ignite_adam", envir = asNamespace("torch"))) torch::optim_ignite_adam else torch::optim_adam
FINN.seed(PARAMS$seed); t_start <- Sys.time()
fit(model, env = env, data = obs, init_cohort = init_cohort, tree_data = tree_obs, loss_aggregation = PARAMS$loss_aggregation,
    epochs = PARAMS$epochs, batchsize = PARAMS$batchsize, lr = PARAMS$lr, lr_scheduler = PARAMS$lr_scheduler,
    weights = PARAMS$weights, loss_family = PARAMS$loss_family,
    patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, optimizer = optimizer, device = PARAMS$device,
    folder = NULL, checkpoints = PARAMS$checkpoints, plot_progress = FALSE)
runtime_h <- as.numeric(difftime(Sys.time(), t_start, units = "hours"))
torch_save(model, file.path(FIT_DIR, "model.pt"))
## release object: every parameter AND buffer by name (reg_logK is a fixed buffer here) + env scaling
parameters <- lapply(model$parameters, function(p) p$detach()$cpu()$clone())
for (nm in names(model$buffers)) parameters[[paste0("buffer.", nm)]] <- model$buffers[[nm]]$detach()$cpu()$clone()
for (i in seq_len(nrow(model$env_scaling))) {
  v <- model$env_scaling$variable[i]
  parameters[[paste0("env_scaling.", v, ".center")]] <- torch_tensor(model$env_scaling$center[i], dtype = torch_float64())
  parameters[[paste0("env_scaling.", v, ".scale")]]  <- torch_tensor(model$env_scaling$scale[i],  dtype = torch_float64())
}
torch_save(parameters, file.path(FIT_DIR, "parameters.pt"))
hist <- as.data.table(do.call(rbind, model$history)); setnames(hist, paste0("V", 1:7), c("dbh", "ba", "trees", "growth", "mortality", "regeneration", "unused"), skip_absent = TRUE); hist[, epoch := .I]
fwrite(hist, file.path(CAL_DIR, "loss_history.csv")); fwrite(hist, file.path(FIT_DIR, "loss_history.csv"))

## ---- 3b. round trip through the package's build_model() (bit-identical) ------------------
species_table <- function() utils::read.csv(file.path("inst", "model", "species.csv"))
source(file.path("R", "definition.R"))
rebuilt <- build_model("process", torch_load(file.path(FIT_DIR, "parameters.pt")))
check_ids <- sort(unique(obs$siteID))[seq_len(200L)]
roundtrip <- function(m) { m$eval(); FINN.seed(PARAMS$seed)
  as.data.table(simulateForest(m, env = env[siteID %in% check_ids], init_cohort = makeInitCohorts(train$init[siteID %in% check_ids, .(siteID, patchID, species, dbh, trees, treeName)], Nspecies = n_species),
                               patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, device = PARAMS$device, batchsize = 200L)$wide$site) }
stopifnot(identical(roundtrip(model), roundtrip(rebuilt)))
message("round trip: build_model(parameters.pt) reproduces the fitted model bit for bit on ", length(check_ids), " sites")

## ---- 4. release metadata and fit summary -------------------------------------------------
sha <- unname(tools::sha256sum(file.path(FIT_DIR, "parameters.pt"))); sha_model <- unname(tools::sha256sum(file.path(FIT_DIR, "model.pt")))
release_tag <- "v0.2.0"
weights <- list(variants = list(process = list(
  tag = release_tag, file = "parameters.pt", url = sprintf("https://github.com/FINNverse/FINN-model-IFN-France/releases/download/%s/parameters.pt", release_tag), sha256 = sha,
  model_file = "model.pt", model_url = sprintf("https://github.com/FINNverse/FINN-model-IFN-France/releases/download/%s/model.pt", release_tag), model_sha256 = sha_model,
  finn_commit = "0848b60", finn_version = as.character(packageVersion("FINN")), torch_version = as.character(packageVersion("torch")),
  patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, data_version = PARAMS$prep_version, train_set = PARAMS$train_set, epochs = PARAMS$epochs,
  covariates = PARAMS$env_vars, fitted = format(Sys.time(), "%Y-%m-%d"))))
writeLines(jsonlite::toJSON(weights, auto_unbox = TRUE, pretty = TRUE, digits = NA), file.path("inst", "model", "weights.json"))
last <- hist[.N]
writeLines(c("# Fit summary: fiveyear-v2-process (v0.2)", "",
  sprintf("Fitted %s by calibration/v02/fit.R; FINN %s (main 0848b60), torch %s, R %s; %d thread(s); runtime %.2f h for %d epochs (%.1f s/epoch).",
          format(Sys.time(), "%Y-%m-%d %H:%M %Z"), packageVersion("FINN"), packageVersion("torch"), getRversion(), torch_get_num_threads(), runtime_h, PARAMS$epochs, runtime_h * 3600 / PARAMS$epochs), "",
  sprintf("Data: %s/%s, %d sites, %d species. Removal per tallied tree, theta frozen, K fixed, NB in count space (w %.4f), cosine lr. Weights `auto` (no numeric vector; weight_scale unused).",
          PARAMS$prep_version, PARAMS$train_set, uniqueN(obs$siteID), n_species, PARAMS$recruit_obs_weight), "",
  sprintf("Release assets: parameters.pt (%d tensors incl. buffers and env scaling) sha256 %s; model.pt sha256 %s; round trip bit-identical on %d sites.", length(parameters), sha, sha_model, length(check_ids)), "",
  "Final losses (baseline-scaled, last epoch):", "", paste0("- ", names(last), ": ", vapply(last, function(x) format(signif(as.numeric(x), 4)), character(1))), "",
  "## PARAMS", "", "```r", capture.output(str(PARAMS, give.attr = FALSE)), "```"), file.path(CAL_DIR, "fit_summary.md"))
message("done: ", FIT_DIR)
