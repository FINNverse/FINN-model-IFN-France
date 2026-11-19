#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-zoo | FINN-model-IFN-France/calibration/v020-soil/fit.R
##
## VARIANT v020-soil (2026-11-14): the v020-defaults joint fit
##        (calibration/v020/fit.R) on the data version fiveyear-v3-soil
##        (calibration/v020-soil/prepare_fiveyear_v3_soil.R = fiveyear-v2 plus
##        two static soil covariates per plot from the IGN soil pit, ECOLOGIE.csv,
##        on the soil scales shared by the zoo inventories
##        and defined in prepare_fiveyear_v3_soil.R) with TWO soil covariates added to env_vars: soil_wet
##        (hydromorphy 0 none within 1 m .. 6 peat) and soil_nutrient_poor
##        (1 eutrophic .. 6 dystrophic). They enter every process that gets the
##        climate in v020, through env_formula: the growth and mortality
##        networks (env_dnn_linmgmt: the climate network now sees tavg_c,
##        prec_mm, soil_wet, soil_nutrient_poor; gest stays the linear term
##        outside it) and the linear regeneration response; competition has no
##        env (~0), as in v020. The soil columns are z-scaled like the climate
##        (env_autoscale "auto" over the training rows). Why: no fiveyear-v2
##        covariate describes the soil, while French forest types are sorted
##        by water and nutrients (alder and ash on wet soils, maritime pine on
##        poor sands), and the general
##        model will pool the inventories on the same two scales. On the
##        national plots the basal-area-weighted soil_wet orders alder 3.0 >
##        aspen 2.0 > pedunculate oak 1.5 > ... > beech 0.37, and
##        soil_nutrient_poor puts maritime pine poorest (4.8), ash 1.9 and
##        Aleppo pine 1.7 richest (PREPARATION.md of fiveyear-v3-soil).
##        Library: FINN release-0.2.0 at 20b184f (~/R/lib-finn-release-20b184f,
##        laptop and cluster) instead of b6244e5; between the two: the help
##        texts of ?recruitment, regeneration_adult_effect() (not used here), a
##        never-fitted model draws deaths as a fitted one does, and
##        simulateForest(init_cohort = NULL) starts from bare ground (used only
##        by build_networks below, two years on two sites). The fit API and
##        every default are the same (formals compared on 2026-11-14).
##        Everything else is identical to calibration/v020/fit.R (seed,
##        schedule, batch 100, defaults, checks, outputs); the rest of this
##        header is v020's and still holds, except DATA / RUN / OUTPUT.
##
## WHAT   The zoo refit of the French NFI model (IGN, 5-year remeasurement,
##        data fiveyear-v2, 30 species) with the FINN 0.2.0 release defaults
##        ("v020-defaults", the set-up shared by every zoo model): ONE
##        joint fit() of competition, growth, mortality and regeneration against
##        all six responses (dbh, ba, trees, growth, mortality, regeneration),
##        loss_weights "auto", update_step 1. It replaces the two-stage v0.2 recipe
##        (calibration/v02/recipe: stage A with the observed recruits fed in,
##        stage B regeneration_saturation with frozen demography and a frozen
##        ML theta); the France stage A of 2026-11-10 timed out on the cluster
##        at epoch 112 of 300 (22 h) and no recipe model exists.
## WHY    The zoo uses the 0.2.0 default set-up (decided on 2026-11-13):
##        every model of the zoo is refitted with it, then the zoo report is
##        updated. The dataset-specific inventory handling stays as prepared.
## FINN   Release candidate release-0.2.0 at b6244e5 (PR FINNverse/FINN#68),
##        installed in ~/R/lib-finn-release-b6244e5 (laptop and cluster).
##        Every 0.2.0 default is left at its default and recorded in
##        fit_summary.md from the built model and formals(): regeneration_capped
##        (cap "all", prior "max_obs_recruits" = max(reg / period_length) per ha
##        and year, prior_lambda 0.1), theta trained with theta_floor 0.5,
##        mortality_draw "binomial", regeneration_sampler "refined", mortality
##        limits 0.001 / 0.999, reg_floor 1e-3, trees_offset 0.01, rb_gradient
##        TRUE, light_mask TRUE, competition height bound log(2.3). No custom
##        floors, no gamma / poisson draw, no weight_scale "raw", no frozen
##        stage, no fed recruits.
## MODEL STRUCTURE (as the v0.2 recipe fit, the only change being the
##        regeneration process):
##   - U-shaped mortality `mortality_ushape` (sen in [0, sen_max] on the logit at
##     dbh 1 m; growth enters as growth / 0.01);
##   - growth and mortality env = module env_dnn_linmgmt: SELU network hidden
##     20-20, dropout 0.1, on the WorldClim normals (tavg_c, prec_mm) PLUS one
##     linear per-species weight on the IGN management trace GEST (0 none, 1
##     old, 2 recent, recorded at visit 1) OUTSIDE the network, so forcing gest
##     shifts the growth / mortality logit by exactly beta_s per GEST step
##     (pure-shift check below); GEST on its 0/1/2 scale (env_autoscale identity);
##   - static climate (WorldClim 2.1 normals, env rows 1..5 identical per plot);
##   - regeneration = FINN::regeneration_capped with a LINEAR env response on
##     tavg_c + prec_mm + gest (the 0.2.0 default process instead of the
##     recipe's frozen regeneration_saturation);
##   - recruit_obs_weight = the 6 m IGN circle (0.0113 ha), recruits_dbh 8.3 cm,
##     one 0.0707 ha patch per plot, growth_rate "annual", growth and mortality
##     scored per tallied tree (tree_obs_dt, period_mode per_site).
## INVENTORY PROBLEMS of the IGN data and how this fit handles them (all in the
##        prepared data fiveyear-v2, unchanged here; calibration/AUDIT_v02.md,
##        calibration/v02/PREPARATION_v2.md):
##   1. Rare species / the old "other" class (A5): the v1 "other" (126 codes, 5 %
##      of stems) took 80 % of the excess recruitment. fiveyear-v2 has NO lumped
##      class: S = the 30 most frequent species, a plot is kept when out-of-S
##      species hold < 2 % of its START-census basal area (48,190 of 59,960
##      plots), out-of-S trees are dropped; the END-census out-of-S share is a
##      validation stratum (end_stratum), never a filter.
##   2. Five-year remeasurement (A3, trap 1): obs at year 5, period_length 5,
##      env rows 1..5 (static climate); rates annualised in the preparation;
##      fit() compares the observation after exactly five steps; recruitment is
##      summed over the window, the cap prior is max(reg / 5) per ha and year.
##   3. Felling / harvest (A1, A3): the mortality target is TOTAL REMOVAL (dead,
##      fallen, broken or cut; cut trees are 48,254 of 747,100 visit-1 trees,
##      19 % of the plots), only not-found trees censored, NO year-5 state masked
##      (v1 censored the cut trees and masked their plots' end state). GEST is
##      the covariate that separates managed from unmanaged plots (removal 1.8 /
##      2.4 / 2.7 %/yr by GEST 0/1/2); the tool simulates with gest = 0 (no
##      trace) or 1 / 2. validation.R reports the natural-removal AUC without
##      the cut trees (`cut` in tree_obs_dt).
##   4. Annualisation trap (A3): 30 % of the site x species rows have one tree
##      at risk, so the site-level `mort` is inflated 2-2.5x (pooled 6.0 %/yr
##      vs 2.47 %/yr per tree). Mortality (and growth) are therefore scored PER
##      TALLIED TREE (loss_aggregation "tree", tree_obs_dt), which has no such
##      bias; the site-level `mort` stays in obs_dt for comparability only.
##   5. Ingrowth vs recruitment (A2): a "new" tree is a recruit only if its
##      back-modelled start diameter (dbh at visit 2 minus the species q90 of the
##      small-tree 5-yr increment) is below the 7.5 cm tally threshold; 3,295 of
##      36,616 new trees were dropped as missed at visit 1.
##   6. Concentric-circle design and expansion (A3): every tallied tree is
##      expanded to its design weight (PB 6.25, MB 2.78, GB 1.0 virtual copies on
##      the 0.0707 ha patch, largest-remainder rounding per plot x species;
##      basal area preserved to -0.2 %); sample sizes (n_at_risk, growth_n) are
##      tallied trees. Recruits are tallied on the 6 m circle, so the recruitment
##      observation enters the negative binomial as reg x 0.0113 ha (whole
##      counts to +-0.1, the virtual-copy rounding; count_tol below).
##   7. Recruit dispersion collapse (A4): v0.1's theta fell to 0.016 for 33 of
##      36 species. Here theta is trained above FINN's default floor 0.5 (the
##      0.2.0 guard rail); the fitted values and the predicted / observed
##      recruitment are reported, and `recruitment` (FINN help) warns that a
##      floored theta does not calibrate the recruitment level.
##   8. Coordinates are the 1 km grid-cell centre (fuzzed <= 700 m): no
##      coordinates as predictors; holdout = 100 km blocks.
##   9. Plots with planting / sowing / clearing in the interval, poplar
##      plantations and plots without GEST are dropped in the preparation.
## SCHEDULE 300 epochs at lr 0.01, then a cosine decay to 1e-4 over 200 epochs
##        (500 in total), plain Adam, batch 100 (as the recipe), checkpoints
##        (parameter snapshot + model file) every 25 epochs; 24 h cluster
##        segments continued by fit_continue.R for as long as the fit needs
##        (decided on 2026-11-13: no compromise on the batch for the
##        cluster budget; the recipe's stage A at batch 100 ran up to ~700
##        s/epoch on a loaded cluster core). Convergence is judged on the
##        parameter drift between checkpoints (convergence.csv, validation.R
##        fig6), never on the loss.
## DATA   ../../FINN-zoo-data/frnfi/prepared/fiveyear-v3-soil/train
##        (holdout_blocked for validation.R)
## RUN    Rscript calibration/v020-soil/fit.R     (from the repo root, one thread)
##        cluster: bash hpc/submit_v020soil.sh (fit, fit_continue x4, validation, chained)
## OUTPUT FIT_DIR = ../../FINN-zoo-data/frnfi/fits/fiveyear-v3-soil-v020 (outside git):
##        model_init.pt, model.pt, checkpoints/epoch_<n>model.pt, loss_history.csv,
##        lr_schedule.csv, species_parameters.csv, regeneration_parameters.csv,
##        management_weights.csv, module_checks.csv, grad_norms.csv,
##        grad_norms_by_epoch.csv, cap_trajectory.csv, convergence.csv, fit_summary.md
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-release-20b184f"), .libPaths()))   # VARIANT: FINN release-0.2.0 at 20b184f (v020: b6244e5; same API and defaults)
suppressPackageStartupMessages({
  library(FINN)
  library(data.table)
  library(torch)
})

DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "frnfi")
FIT_DIR  <- file.path(DATA_DIR, "fits", "fiveyear-v3-soil-v020")   # VARIANT: its own folder (v020: fiveyear-v2-v020)

PARAMS <- list(
  ## ---- per model: what changes from data set to data set ----
  prep_version  = "fiveyear-v3-soil",      # VARIANT: fiveyear-v2 (30 species, no "other", start-census site rule x = 2 %, 5-yr interval, plots first visited 2010-2020) + static soil columns
  train_set     = "train",                 # 35,783 training plots in the 66 retained 100 km blocks
  holdout_set   = "holdout_blocked",       # not read here; recorded for validation.R
  env_vars      = c("tavg_c", "prec_mm", "soil_wet", "soil_nutrient_poor", "gest"),   # WorldClim 2.1 normals (bio1, bio12) + VARIANT: the two harmonised soil scores (static per plot) + the IGN management trace GEST at visit 1 (0/1/2)
  linear_var    = "gest",                  # outside the network: one linear per-species weight per GEST step; the climate network sees the other columns
  mgmt_levels   = c(0, 1, 2),              # the values the covariate takes (0 none, 1 old trace, 2 recent); the pure-shift check forces each of them
  patch_size_ha = pi * 15^2 / 1e4,         # one FINN patch = the 15 m IGN circle (0.0707 ha); every tallied tree expanded to its design weight (PB 6.25, MB 2.78, GB 1 copies)
  n_patches     = 1L,
  period        = 5L,                      # the IGN remeasurement interval: obs at year 5, env rows 1..5
  recruits_dbh  = 8.3,                     # recruits enter just above the 7.5 cm tally threshold (median dbh of the fiveyear-v2 recruits at visit 2)
  growth_rate   = "annual",                # `growth` in obs_dt and tree_obs_dt is the annualised relative dbh increment (PREPARATION_v2.md)
  recruit_obs_weight_ha = pi * 6^2 / 1e4,  # the area one observed recruit stands for: the 6 m circle (0.0113 ha) on which IGN tallies the 7.5-22.5 cm trees, so reg x w is the tallied count
  count_tol     = 0.1,                     # |reg x w - round(reg x w)| allowed: the virtual-copy rounding of the preparation (6.25 copies per recruit, largest remainder)
  check_rows    = 200L,                    # env rows of the module checks (forward, gradients, pure shift)
  check_tol     = 1e-5,                    # |output change - beta * delta| / (1 + |output|) allowed in the pure-shift check (float32 rounding)
  ## ---- the v020 set-up (same in every model) ----
  hidden        = c(20L, 20L),             # climate network of growth and mortality (SELU), as the recipe
  dropout       = 0.1,                     # active while training, off in simulateForest (eval mode)
  beta_init     = 0,                       # management weight start
  sen_max       = 5,                       # U-shape cap: at most +5 on the mortality logit at dbh 1 m
  sen_init      = 0.05,                    # U-shape start: +0.05 on the logit at 1 m
  growth_scale  = 0.01,                    # growth enters mortality as growth / 0.01 (FINN 0.2.0's own default, asserted)
  n_quantiles   = 10L,                     # FINN's default height-quantile competition, written out
  light_mask    = TRUE,                    # 0.2.0's masked competition (FINN's default, written out)
  loss_family   = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse",
                    mortality = "binomial", regeneration = "nbinom"),   # FINN's default families, all six responses scored
  loss_weights  = "auto",                  # each term / its intercept-only baseline (FINN's default)
  loss_recruit  = "simulated",             # the regeneration process recruits and is estimated with the rest (FINN's default; no fed recruits)
  loss_aggregation = list(growth = "tree", mortality = "tree"),   # per tallied tree on the observed stand (inventory problem 4); needs period_mode per_site
  clip_norm     = 2,                       # FINN's default per-group clip budget, written out
  update_step   = 1L,                      # FINN's default
  epochs_const  = 300L,                    # epochs at constant lr, then the cosine decay
  anneal_epochs = 200L,                    # epochs of the cosine decay after them
  lr            = 0.01,                    # the constant lr = the start of the cosine
  lr_min        = 1e-4,                    # the end of the cosine (1 % of the start)
  checkpoints   = 25L,                     # a parameter snapshot and a model file (resume point for fit_continue.R) every 25 epochs
  batchsize     = 100L,                    # as the recipe (~358 steps per epoch); the fit runs as long as it needs (header)
  threads       = 1L,                      # one torch thread (one CPU per fit job)
  device        = "cpu",
  seed          = 1L                       # the recipe's seed
)
EPOCHS <- PARAMS$epochs_const + PARAMS$anneal_epochs
torch_set_num_threads(PARAMS$threads)
if (file.exists(file.path(FIT_DIR, "model.pt"))) {
  message("model.pt exists: the fit has finished"); quit(save = "no")
}
dir.create(file.path(FIT_DIR, "checkpoints"), recursive = TRUE, showWarnings = FALSE)
if (length(list.files(file.path(FIT_DIR, "checkpoints"), pattern = "model\\.pt$", recursive = TRUE)))
  stop("checkpoints exist already: resume with fit_continue.R instead of starting over")
## the installed FINN must be the release-0.2.0 API (defaults are read, not set)
stopifnot(packageVersion("FINN") == "0.2.0",
          basename(dirname(find.package("FINN"))) == "lib-finn-release-20b184f",   # VARIANT: the library this fit names, not another 0.2.0 build
          identical(formals(FINN::mortality)$growth_scale, PARAMS$growth_scale),
          all(c("loss_weights", "loss_recruit", "control", "rb_gradient", "trees_offset", "weight_scale") %in% names(formals(fit))),
          !any(c("weights", "recruits", "lr_scheduler", "optimizer", "clip_norm", "deaths") %in% names(formals(fit))),
          all(c("cap", "prior", "prior_lambda", "theta_floor", "optimizeTheta", "mortality_draw", "mortality_floor",
                "mortality_ceiling", "regeneration_sampler", "light_mask", "reg_floor") %in% names(formals(createProcess))),
          all(c("growth_rate", "recruit_obs_weight") %in% names(formals(finn))),
          is.function(FINN::regeneration_capped), is.function(FINN::finnControl))
DEFAULTS <- list(trees_offset = eval(formals(fit)$trees_offset), rb_gradient = eval(formals(fit)$rb_gradient))
PS    <- PARAMS$patch_size_ha
w_obs <- PARAMS$recruit_obs_weight_ha

## ---- the learning-rate schedule over the whole run --------------------------------
## lr_of(t), t = epochs completed; FINN steps the scheduler once per epoch, so epoch
## t + 1 runs at lr_of(t). torch::lr_lambda multiplies the base lr by lr_lambda(e),
## e = epochs stepped in THIS fit() call (fit_continue.R sets k_done).
lr_of <- function(t) ifelse(t < PARAMS$epochs_const, PARAMS$lr,
  PARAMS$lr_min + (PARAMS$lr - PARAMS$lr_min) * (1 + cos(pi * (t - PARAMS$epochs_const) / PARAMS$anneal_epochs)) / 2)
k_done <- 0L
schedule <- function(optimizer) torch::lr_lambda(optimizer, lr_lambda = function(e) lr_of(k_done + e) / PARAMS$lr)
SCHED <- data.table(epoch = seq_len(EPOCHS))[, lr := lr_of(epoch - 1L)]

## ---- 1. data -----------------------------------------------------------------
train_dir <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$train_set)
obs   <- fread(file.path(train_dir, "obs_dt.csv"))[, .(siteID, year, species, ba, dbh, trees, growth, growth_n, mort, n_at_risk, n_died, reg, period_length)]
env   <- fread(file.path(train_dir, "env_dt.csv"))[, c("siteID", "year", PARAMS$env_vars), with = FALSE]
init  <- fread(file.path(train_dir, "init_cohort_dt.csv"))[, .(siteID, patchID, species, dbh, trees, treeName)]
trees_all <- fread(file.path(train_dir, "tree_obs_dt.csv"))
trees <- trees_all[, .(siteID, year, patchID, species, dbh, trees, growth, died, period_length)]   # the columns FINN's per-tree losses read; `cut` / `died_natural` stay diagnostics
species_dt <- fread(file.path(train_dir, "species_dt.csv"))
n_species <- nrow(species_dt)
stopifnot(max(obs$species) == n_species, all(obs$year == PARAMS$period), all(obs$period_length == PARAMS$period),
          max(env$year) == PARAMS$period, uniqueN(init$patchID) == PARAMS$n_patches,
          all(env[[PARAMS$linear_var]] %in% PARAMS$mgmt_levels), !anyNA(obs$reg), all(obs$reg >= 0))
soil_vars <- grep("^soil_", PARAMS$env_vars, value = TRUE)   # VARIANT: the soil columns are complete and static per plot
stopifnot(length(soil_vars) == 2L, !anyNA(env[, ..soil_vars]),
          env[, all(vapply(.SD, uniqueN, integer(1)) == 1L), by = siteID, .SDcols = soil_vars][, all(V1)])
init_cohort <- makeInitCohorts(init, Nspecies = n_species)
## the recruit observation space: reg (per ha per 5 yr) x the 6 m circle = tallied recruits
count_dev <- max(abs(obs$reg * w_obs - round(obs$reg * w_obs)))
stopifnot(count_dev < PARAMS$count_tol)
K_prior_expected <- max(obs$reg / obs$period_length)   # what prior = "max_obs_recruits" resolves to (per ha per year)
mgmt_share <- env[year == 1, 100 * prop.table(table(factor(get(PARAMS$linear_var), levels = PARAMS$mgmt_levels)))]
message(sprintf("train: %d plots, %d species, interval %d yr, %d tallied trees (%d removed, %d of them cut); gest 0/1/2 = %s %% of the plots; reg > 0 in %d of %d rows; reg x %.5f ha whole to %.3f; max reg / period %.1f per ha per year",
                uniqueN(obs$siteID), n_species, PARAMS$period, nrow(trees), sum(trees$died), sum(trees_all$cut),
                paste(sprintf("%.1f", mgmt_share), collapse = " / "), sum(obs$reg > 0), nrow(obs), w_obs, count_dev, K_prior_expected))
env_formula <- as.formula(paste("~", paste(PARAMS$env_vars, collapse = " + ")))
## scaling per column: climate "auto" (z-scaled), the management trace "identity" (beta = the shift of one GEST step)
autoscale_spec <- c(setNames(rep("auto", length(setdiff(PARAMS$env_vars, PARAMS$linear_var))), setdiff(PARAMS$env_vars, PARAMS$linear_var)),
                    setNames("identity", PARAMS$linear_var))

## ---- 2. the env module of growth and mortality: MLP(climate) + linear management ----
## FINN's env design matrix is stats::model.matrix(formula, env): the intercept
## column first, then the formula's terms in order. The management column's
## position is read from it, not assumed. The module: input [..., n_inputs] (FINN
## passes [sites, years, n] or [sites, n]), output [..., n_species]. mlp = FINN's
## own network layout (linear, SELU, dropout after every hidden activation, linear
## output with bias) on every column except management (the constant intercept
## column stays in, as in FINN's networks); beta = one linear weight per species
## on the management column, added outside the network.
X_names <- colnames(model.matrix(env_formula, env[1:5, PARAMS$env_vars, with = FALSE]))
mgmt_col <- match(PARAMS$linear_var, X_names)
stopifnot(identical(X_names, c("(Intercept)", PARAMS$env_vars)), !is.na(mgmt_col))
env_dnn_linmgmt <- torch::nn_module("env_dnn_linmgmt",
  initialize = function(n_inputs, n_species, mgmt_col, hidden, dropout, beta_init) {
    self$mgmt_col  <- as.integer(mgmt_col)
    self$clim_cols <- setdiff(seq_len(n_inputs), self$mgmt_col)
    layers <- list(); n_in <- length(self$clim_cols)
    for (h in hidden) {
      layers <- c(layers, torch::nn_linear(n_in, h), torch::nn_selu())
      if (dropout > 0) layers <- c(layers, torch::nn_dropout(p = dropout))
      n_in <- h
    }
    self$mlp  <- do.call(torch::nn_sequential, c(layers, torch::nn_linear(n_in, n_species)))
    self$beta <- torch::nn_parameter(torch::torch_tensor(beta_init, dtype = torch::torch_float32()))
  },
  forward = function(x) {
    self$mlp(x[.., self$clim_cols]) + x[.., self$mgmt_col, drop = FALSE] * self$beta
  })
new_env_module <- function() env_dnn_linmgmt(n_inputs = length(X_names), n_species = n_species, mgmt_col = mgmt_col,
                                             hidden = PARAMS$hidden, dropout = PARAMS$dropout,
                                             beta_init = rep(PARAMS$beta_init, n_species))
## the checks on check_rows real design-matrix rows (climate z-scaled as FINN will
## see it): (1) forward pass in both call forms FINN uses; (2) gradients reach every
## MLP tensor and beta; (3) forcing management to each of its levels changes the
## output by exactly beta_s * delta on every row (eval mode = without dropout)
scaling_chk <- compute_env_scaling(env, autoscale_spec)
X_chk <- model.matrix(env_formula, apply_env_scaling(env, scaling_chk)[round(seq(1, .N, length.out = PARAMS$check_rows)), PARAMS$env_vars, with = FALSE])
stopifnot(nrow(X_chk) %% 20L == 0L)
check_module <- function(nn, label, grads = TRUE) {
  x <- torch::torch_tensor(unname(X_chk), dtype = torch::torch_float32())
  was_training <- nn$training
  out <- list(module = label, rows = nrow(X_chk))
  if (grads) {
    nn$train()
    y3 <- nn(x$view(c(20L, -1L, ncol(X_chk))))
    stopifnot(identical(as.integer(y3$shape), c(20L, nrow(X_chk) %/% 20L, n_species)))
    g <- torch::autograd_grad(nn(x)$pow(2)$sum(), nn$parameters)
    g <- setNames(vapply(g, function(t) as.numeric(t$abs()$sum()), numeric(1)), names(nn$parameters))
    stopifnot(all(is.finite(g)), all(g[names(g) != "beta"] > 0))
    out$grad_mlp_min_abs_sum <- min(g[names(g) != "beta"]); out$grad_beta_abs_sum <- g[["beta"]]
  }
  nn$eval()
  torch::with_no_grad({
    y2 <- as_array(nn(x))
    stopifnot(identical(dim(y2), c(nrow(X_chk), n_species)))
    beta <- as_array(nn$beta)
    dev <- 0
    for (t in PARAMS$mgmt_levels) {
      xt <- X_chk; xt[, mgmt_col] <- t
      yt <- as_array(nn(torch::torch_tensor(unname(xt), dtype = torch::torch_float32())))
      shift <- outer(xt[, mgmt_col] - X_chk[, mgmt_col], beta)
      dev <- max(dev, max(abs(yt - y2 - shift) / (1 + abs(y2))))
    }
  })
  if (was_training) nn$train()
  stopifnot(dev < PARAMS$check_tol)
  out$shift_max_rel_dev <- dev
  out$max_abs_output <- max(abs(y2))
  as.data.table(out)
}

## ---- 3. the model ------------------------------------------------------------------
## U-shaped mortality (the recipe's function; growth_scale and sen_max are buffers
## so that nothing depends on a formal's default)
mortality_ushape <- function(dbh, species, trees, parMort, pred, light, base_steepness = 5, debug = FALSE, growth = NULL) {
  if (is.null(growth)) growth <- self$g
  d   <- dbh / 100
  sen <- (self$mort_sen_max * torch::torch_sigmoid(self$mort_sen_raw))[species]
  torch::torch_sigmoid(pred + growth / self$mort_growth_scale * parMort[, 3][species] + light * parMort[, 1][species] +
                       parMort[, 2][species] * d + sen * d * d)
}
sen_raw_init <- rep(qlogis(PARAMS$sen_init / PARAMS$sen_max), n_species)
build_model <- function() finn(
  N_species    = n_species,
  recruits_dbh = PARAMS$recruits_dbh,
  growth_rate  = PARAMS$growth_rate,
  recruit_obs_weight = w_obs,              # the inventory observation operator (header)
  competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE,
                                       n_quantiles = PARAMS$n_quantiles, light_mask = PARAMS$light_mask),
  growth_process       = createProcess(env_formula, FINN::growth, NN = new_env_module(),
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(env_formula, mortality_ushape, NN = new_env_module(),   # death draw and limits: FINN's defaults
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE,
                                       custom_parameters = list(mort_sen_raw = list(init = sen_raw_init, trainable = TRUE),
                                                                mort_sen_max = list(init = PARAMS$sen_max, trainable = FALSE),
                                                                mort_growth_scale = list(init = PARAMS$growth_scale, trainable = FALSE))),
  regeneration_process = createProcess(env_formula, FINN::regeneration_capped,   # FINN's default process; cap, prior, sampler, theta, floor: defaults
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE)
)
## FINN builds the env networks lazily; a two-site, two-year bare-ground
## simulation builds them so that model_init.pt survives torch_load. The climate
## is z-scaled here as fit() will scale it (the model carries no scaling yet):
## on the raw normals the untrained linear regeneration response overflows.
build_networks <- function(m) {
  e <- apply_env_scaling(env, scaling_chk)[siteID %in% sort(unique(siteID))[1:2] & year <= 2L]
  e[, siteID := match(siteID, sort(unique(siteID)))]
  invisible(simulateForest(m, env = e, init_cohort = NULL, patches = 1L, patch_size = PS, device = "cpu"))
}
FINN.seed(PARAMS$seed)                     # species parameters and the networks are initialised here
model <- build_model()
build_networks(model)
stopifnot(inherits(model$nn_growth, "env_dnn_linmgmt"), inherits(model$nn_mortality, "env_dnn_linmgmt"),
          !any(c("mort_sen_max", "mort_growth_scale") %in% names(model$parameters)),
          all(c("par_theta_recruits_raw", "reg_logK") %in% names(model$parameters)),
          identical(model$reg_cap, "all"), model$reg_logK$numel() == 1L,
          identical(model$theta_floor, 0.5),   # FINN's default floor on theta (header)
          abs(model$recruit_obs_weight - w_obs) < 1e-12)
CHECKS <- rbind(cbind(stage = "init", check_module(model$nn_growth, "growth")),
                cbind(stage = "init", check_module(model$nn_mortality, "mortality")))
print(CHECKS)
## what the installed FINN made of the defaults (recorded in fit_summary.md)
USED <- list(mortality_draw = model$mortality_draw, mortality_floor = model$mortality_floor, mortality_ceiling = model$mortality_ceiling,
             theta_floor = model$theta_floor, theta_optimized = model$theta_optimized, regeneration_sampler = model$regeneration_sampler,
             sample_regeneration = model$sample_regeneration, reg_floor = model$reg_floor, cap = model$reg_cap, prior = model$reg_prior,
             prior_lambda = model$reg_prior_lambda, light_mask = model$light_mask, trees_offset = DEFAULTS$trees_offset, rb_gradient = DEFAULTS$rb_gradient)
message("FINN defaults in use: ", paste(names(USED), vapply(USED, function(x) paste(format(x), collapse = ","), ""), sep = " = ", collapse = "; "))
torch_save(model, file.path(FIT_DIR, "model_init.pt"))   # reg_logK holds FINN's placeholder until fit() resolves the prior
START <- list(theta = as.numeric(model$par_theta_recruits$detach()$cpu()))
model$period_mode <- "per_site"            # FINN scores loss_aggregation = "tree" on its per-site path only

## ---- 4. the optimizer: plain Adam that logs the clipped gradient norms ------------
## FINN clips each parameter group (mechanistic, nn, loss) to clip_norm and then
## calls optimizer$step(); the subclass records, per step and parameter, the L2
## norm of the gradient it receives, and the cap K after the step (cap_trajectory.csv),
## then runs Adam's step unchanged.
GRADLOG <- new.env(); GRADLOG$rows <- list(); GRADLOG$K <- numeric(0)
optim_adam_gradlog <- torch::optimizer("optim_adam_gradlog", inherit = torch::optim_adam,
  initialize = function(params, lr, ...) {
    GRADLOG$names <- names(params)
    super$initialize(params, lr = lr, ...)
  },
  step = function(closure = NULL) {
    ps <- self$param_groups[[1]]$params
    GRADLOG$rows[[length(GRADLOG$rows) + 1L]] <- vapply(ps, function(p) tryCatch({
      g <- p$grad
      if (is.null(g) || !g$numel()) 0 else as.numeric(g$norm()$item())
    }, error = function(e) 0), numeric(1))
    out <- super$step(closure)
    GRADLOG$K <- c(GRADLOG$K, exp(as.numeric(ps[[match("reg_logK", GRADLOG$names)]]$detach()$cpu())))
    out
  })

## ---- 5. the joint fit ----------------------------------------------------------------
## nothing is meant to be fixed: every parameter must require gradient (header)
stopifnot(all(vapply(model$parameters, function(p) isTRUE(p$requires_grad), logical(1))))
FINN.seed(PARAMS$seed)
t0 <- Sys.time()
fit(model,
    env              = env,
    data             = obs,                # all six responses; max(reg / period_length) resolves the cap's prior
    tree_data        = trees,
    loss_aggregation = PARAMS$loss_aggregation,
    init_cohort      = init_cohort,
    epochs           = EPOCHS,
    batchsize        = PARAMS$batchsize,
    update_step      = PARAMS$update_step,
    lr               = PARAMS$lr,          # the base lr the schedule multiplies
    loss_weights     = PARAMS$loss_weights,
    loss_family      = PARAMS$loss_family,
    loss_recruit     = PARAMS$loss_recruit,
    patches          = PARAMS$n_patches,
    patch_size       = PS,
    control          = finnControl(optimizer = optim_adam_gradlog,   # plain Adam + the gradient-norm log
                                   lr_scheduler = schedule,          # constant, then cosine; stepped once per epoch by FINN
                                   clip_norm = PARAMS$clip_norm),
    device           = PARAMS$device,
    env_autoscale    = autoscale_spec,
    folder           = file.path(FIT_DIR, "checkpoints"),
    checkpoints      = PARAMS$checkpoints,
    plot_progress    = FALSE)              # trees_offset, rb_gradient, weight_scale: FINN's defaults (header)
min_fit <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
lr_end <- model$optimizer$param_groups[[1]]$lr
torch_save(model, file.path(FIT_DIR, "model.pt"))
hist <- as.data.table(do.call(rbind, model$history))[, 1:6]
setnames(hist, c("dbh", "ba", "trees", "growth", "mort", "reg"))
hist[, epoch := .I]
fwrite(hist, file.path(FIT_DIR, "loss_history.csv"))
fwrite(SCHED, file.path(FIT_DIR, "lr_schedule.csv"))
K_prior <- model$reg_prior_value
K_fit   <- exp(as.numeric(model$reg_logK$detach()$cpu()))
stopifnot(nrow(hist) == EPOCHS, abs(lr_end - lr_of(EPOCHS)) < 1e-9,
          abs(K_prior - K_prior_expected) < 1e-6 * K_prior_expected,   # max_obs_recruits = max(reg / period_length), per ha per year
          abs(as.numeric(model$mort_sen_max$cpu()) - PARAMS$sen_max) < 1e-6,
          abs(as.numeric(model$mort_growth_scale$cpu()) - PARAMS$growth_scale) < 1e-9)

## ---- 6. diagnostics: gradient norms and the cap ------------------------------------
G <- as.data.table(do.call(rbind, GRADLOG$rows)); setnames(G, GRADLOG$names)
G[, step := .I]
steps_per_epoch <- nrow(G) / EPOCHS
G[, epoch := as.integer(ceiling(step / steps_per_epoch))]
GL <- melt(G, id.vars = c("step", "epoch"), variable.name = "parameter", value.name = "norm")
GL[, group := fifelse(grepl("^par_loss_", parameter), "loss", fifelse(grepl("^nn_", parameter), "nn", "mechanistic"))]
GL[, regeneration := grepl("regeneration|theta_recruits|reg_logK", parameter)]
GS <- GL[, .(group_norm = sqrt(sum(norm^2)), reg_share = sum(norm[regeneration]^2) / max(sum(norm^2), 1e-30)),
         by = .(step, epoch, group)]
GS[, at_budget := group_norm > 0.999 * PARAMS$clip_norm]
fwrite(GL, file.path(FIT_DIR, "grad_norms.csv"))
GE <- GS[, .(group_norm_mean = mean(group_norm), share_clipped = mean(at_budget), reg_share_mean = mean(reg_share)),
         by = .(epoch, group)]
fwrite(GE, file.path(FIT_DIR, "grad_norms_by_epoch.csv"))
KT <- data.table(step = seq_along(GRADLOG$K), K = GRADLOG$K)[, epoch := as.integer(ceiling(step / steps_per_epoch))][
  , .(K_end = K[.N], K_min = min(K), K_max = max(K)), by = epoch][, K_prior := K_prior]
fwrite(KT, file.path(FIT_DIR, "cap_trajectory.csv"))

## ---- 7. checks and tables -------------------------------------------------------------
## the trained modules reloaded from model.pt: still the module class, still a pure shift
chk <- torch_load(file.path(FIT_DIR, "model.pt"))
stopifnot(inherits(chk$nn_growth, "env_dnn_linmgmt"), inherits(chk$nn_mortality, "env_dnn_linmgmt"))
CHECKS <- rbind(CHECKS, cbind(stage = "trained, reloaded", check_module(chk$nn_growth, "growth", grads = FALSE)),
                cbind(stage = "trained, reloaded", check_module(chk$nn_mortality, "mortality", grads = FALSE)), fill = TRUE)
fwrite(CHECKS, file.path(FIT_DIR, "module_checks.csv"))
rm(chk)
spar <- function(block, nm) `colnames<-`(as.data.table(as.matrix(model[[block]]$detach()$cpu())), nm)
SPP <- cbind(species_dt[order(species)], spar("par_competition", c("comp_height", "comp_strength")),
             spar("par_growth", c("growth_light", "growth_size")), spar("par_mortality", c("mort_light", "mort_size", "mort_growth")),
             spar("par_regeneration", "reg_light"))
SPP[, `:=`(sen = PARAMS$sen_max * plogis(as.numeric(model$mort_sen_raw$detach()$cpu())),
           theta_start = START$theta, theta_fit = as.numeric(model$par_theta_recruits$detach()$cpu()))]
SPP[, dbh_min_mortality_cm := fifelse(sen > 1e-3 & mort_size < 0, -mort_size / (2 * sen) * 100, Inf)]
fwrite(SPP, file.path(FIT_DIR, "species_parameters.csv"))
MGMT <- merge(data.table(species = seq_len(n_species), growth_beta = as_array(model$nn_growth$beta$detach()),
                         mortality_beta = as_array(model$nn_mortality$beta$detach())), species_dt, by = "species")
fwrite(MGMT, file.path(FIT_DIR, "management_weights.csv"))   # the shift of log growth / the mortality logit per GEST step
## one simulation of train over fit()'s windows (observation at year p <- sim years
## 1..p): predicted / observed recruitment, and the vanished stem predictions
## (simulated 0, or below FINN's trees_offset, against a positive count)
sim_check <- function(m) {
  m$eval()
  FINN.seed(PARAMS$seed)
  S <- as.data.table(simulateForest(m, env = env, init_cohort = init_cohort, patches = PARAMS$n_patches, patch_size = PS,
                                    device = "cpu", batchsize = PARAMS$batchsize)$wide$site)
  for (cc in c("siteID", "year", "species")) S[, (cc) := as.integer(as.character(get(cc)))]
  W <- unique(obs[, .(siteID, obs_year = year, period_length)])[
    , .(year = (obs_year - period_length + 1L):obs_year), by = .(siteID, obs_year)]
  R <- merge(W, S[, .(siteID, year, species, r_mean_ha)], by = c("siteID", "year"))[, .(pred = sum(r_mean_ha)), by = species]
  R <- merge(R, obs[, .(obs = sum(reg)), by = species], by = "species")
  Z <- merge(obs[trees > 0, .(siteID, year, species)], S[, .(siteID, year, species, trees_m = trees)],
             by = c("siteID", "year", "species"), all.x = TRUE)[is.na(trees_m), trees_m := 0]
  list(total = sum(R$pred) / sum(R$obs), species = R[, .(species, ratio = pred / obs)],
       zero = Z[, .(n_pos = .N, n_zero = sum(trees_m == 0), n_below_offset = sum(trees_m < DEFAULTS$trees_offset))])
}
m_start <- torch_load(file.path(FIT_DIR, "model_init.pt"))
invisible(torch::with_no_grad(m_start$reg_logK$fill_(log(K_prior))))   # the start fit() used: K = the resolved prior
m_start$env_scaling <- model$env_scaling   # and the env scaling fit() learned from the same data before its first step
ck_start <- sim_check(m_start); rm(m_start)
ck_end   <- sim_check(model)
REGT <- merge(merge(SPP[, .(species, species_name, reg_light, theta_start, theta_fit)],
                    ck_start$species[, .(species, rec_ratio_start = ratio)], by = "species"),
              ck_end$species[, .(species, rec_ratio_fit = ratio)], by = "species")
REGT <- merge(REGT, obs[, .(obs_max_reg_per_period = max(reg)), by = species], by = "species")
REGT[, `:=`(K_fit = K_fit, K_prior = K_prior, K_over_obs_max_annual = K_fit / (obs_max_reg_per_period / PARAMS$period))]
fwrite(REGT, file.path(FIT_DIR, "regeneration_parameters.csv"))
ph <- model$param_history
conv <- if (length(ph) >= 2) rbindlist(lapply(names(ph[[length(ph)]]), function(nm)
  data.table(block = nm, last_change = mean(abs(ph[[length(ph)]][[nm]] - ph[[length(ph) - 1]][[nm]]))))) else data.table(block = NA_character_, last_change = NA_real_)
fwrite(conv, file.path(FIT_DIR, "convergence.csv"))

## ---- 8. what was fitted --------------------------------------------------------
commit <- suppressWarnings(tryCatch(system("git rev-parse --short HEAD", intern = TRUE, ignore.stderr = TRUE),
                                    error = function(e) character(0)))
if (!length(commit)) commit <- "unknown (no git checkout)"
sum_l <- function(rows) { rows <- unique(pmin(pmax(rows, 1L), nrow(hist)))
  round(colMeans(as.matrix(hist[rows, .(dbh, ba, trees, growth, mort, reg)])), 4) }
fmt <- function(v) paste(sprintf("%s %.4g", names(v), v), collapse = ", ")
writeLines(c(
  sprintf("# FINN-model-IFN-France, v020-defaults joint fit + soil (%s)", PARAMS$prep_version), "",
  sprintf("VARIANT v020-soil: calibration/v020/fit.R with `soil_wet` and `soil_nutrient_poor` (static, harmonised scales; z-scaled like the climate) added to env_vars, so to the growth and mortality networks and the linear regeneration response; data %s; FINN release-0.2.0 at 20b184f.", PARAMS$prep_version), "",
  sprintf("Fitted %s by `calibration/v020-soil/fit.R` (FINN %s from library `%s`, release-0.2.0; torch %s; repository %s): %.1f min for %d epochs (%.1f s/epoch), update_step %d, %d thread(s) on %s.",
          format(Sys.time(), "%Y-%m-%d %H:%M"), packageVersion("FINN"), basename(dirname(find.package("FINN"))),
          packageVersion("torch"), commit, min_fit, EPOCHS, 60 * min_fit / EPOCHS, PARAMS$update_step, PARAMS$threads, PARAMS$device), "",
  sprintf("FINN defaults in use (none passed): %s.", paste(names(USED), vapply(USED, function(x) paste(format(x), collapse = ","), ""), sep = " = ", collapse = "; ")), "",
  sprintf("Data: `prepared/%s/%s` (%d plots, %d species, interval %d yr, %d tallied trees of which %d removed (%d cut); gest 0 / 1 / 2 on %s %% of the plots). Data convention kept: total removal (dead, fallen, broken or cut) is the loss target, not-found trees censored, no state masked. All six responses scored; recruitment through recruit_obs_weight = %.5f ha (the 6 m circle; reg x w whole to %.3f). Holdout for validation: `%s`.",
          PARAMS$prep_version, PARAMS$train_set, uniqueN(obs$siteID), n_species, PARAMS$period, nrow(trees), sum(trees$died), sum(trees_all$cut),
          paste(sprintf("%.1f", mgmt_share), collapse = " / "), w_obs, count_dev, PARAMS$holdout_set), "",
  sprintf("Plain joint fit; FINN default initialisation (seed %d), theta trained above FINN's default floor %g, no warm start, no L2 other than the cap prior, no fed recruits. Growth and mortality env `%s` through env_dnn_linmgmt (hidden %s, dropout %g, + linear `%s` on its 0/1/2 scale, start %g); U-shaped mortality (sen in [0, %g], start %g at 1 m; growth / %g); light_mask %s, %d quantiles; regeneration FINN::regeneration_capped with a linear env response on the same formula; growth and mortality loss_aggregation = tree (scored on the observed tallied trees, period_mode per_site); growth_rate %s; recruits_dbh %g; %d patch of %.4f ha.",
          PARAMS$seed, model$theta_floor, paste(deparse(env_formula), collapse = ""), paste(PARAMS$hidden, collapse = "-"), PARAMS$dropout,
          PARAMS$linear_var, PARAMS$beta_init, PARAMS$sen_max, PARAMS$sen_init, PARAMS$growth_scale, PARAMS$light_mask, PARAMS$n_quantiles,
          PARAMS$growth_rate, PARAMS$recruits_dbh, PARAMS$n_patches, PS), "",
  sprintf("Module checks (`module_checks.csv`, %d rows): forward pass in both call forms ok; gradients reach every MLP tensor and beta at the start; forcing %s to each of %s moves the output by beta_s x delta with max relative deviation %.2g (start) and %.2g (trained, reloaded) = a pure shift.",
          nrow(X_chk), PARAMS$linear_var, paste(PARAMS$mgmt_levels, collapse = "/"), max(CHECKS[stage == "init", shift_max_rel_dev]), max(CHECKS[stage != "init", shift_max_rel_dev])), "",
  sprintf("Management weights (shift per GEST step; `management_weights.csv`): log growth %s; mortality logit %s.",
          paste(sprintf("%s %.2f", MGMT$species_name, MGMT$growth_beta), collapse = ", "),
          paste(sprintf("%s %.2f", MGMT$species_name, MGMT$mortality_beta), collapse = ", ")), "",
  sprintf("Recruitment cap: prior (max_obs_recruits, per ha per year) %.1f = start K; fitted K %.1f per ha per year (K / prior %.3f).", K_prior, K_fit, K_fit / K_prior), "",
  sprintf("Weights \"auto\": baselines %s; weights %s.", fmt(model$loss_baseline), fmt(model$loss_weights)), "",
  sprintf("Schedule (`lr_schedule.csv`): epochs 1-%d at lr %g, epochs %d-%d cosine %g -> %g; the optimizer ended at lr %.3g. Plain Adam (torch::optim_adam), clip_norm %g, batch %d (%d steps per epoch), checkpoints every %d epochs.",
          PARAMS$epochs_const, PARAMS$lr, PARAMS$epochs_const + 1L, EPOCHS, PARAMS$lr, PARAMS$lr_min, lr_end,
          PARAMS$clip_norm, PARAMS$batchsize, as.integer(round(steps_per_epoch)), PARAMS$checkpoints), "",
  sprintf("Baseline-scaled losses (mean of 10 epochs): start %s; end %s.", fmt(sum_l(1:10)), fmt(sum_l((nrow(hist) - 9):nrow(hist)))), "",
  sprintf("Predicted / observed recruitment on train over the fit windows: %.2f at the start (K = prior), %.2f after the fit (`regeneration_parameters.csv`). Theta: start %s, fit %s (floor %g).",
          ck_start$total, ck_end$total, paste(signif(range(START$theta), 3), collapse = "-"), paste(signif(range(SPP$theta_fit), 3), collapse = "-"), model$theta_floor), "",
  sprintf("Vanished stem predictions on train (positive observed count, simulated 0 / below trees_offset %g): start %d / %d, end %d / %d of %d.",
          DEFAULTS$trees_offset, ck_start$zero$n_zero, ck_start$zero$n_below_offset, ck_end$zero$n_zero, ck_end$zero$n_below_offset, ck_end$zero$n_pos), "",
  sprintf("U-shape: sen %s; dbh of minimum mortality (cm): %s (`species_parameters.csv`).",
          paste(sprintf("%.3f", range(SPP$sen)), collapse = "-"),
          paste(sprintf("%s %.0f", SPP$species_name, SPP$dbh_min_mortality_cm), collapse = ", ")), "",
  sprintf("Clipped gradients (`grad_norms_by_epoch.csv`, last 10 epochs): %s.",
          paste(GE[epoch > EPOCHS - 10, .(s = sprintf("%s clipped %.0f %% of steps, regeneration share %.2f", group[1], 100 * mean(share_clipped), mean(reg_share_mean))), by = group]$s, collapse = "; ")), "",
  sprintf("Largest parameter changes between the last two snapshots (every %d epochs; `convergence.csv`): %s.", PARAMS$checkpoints,
          paste(sprintf("%s %.2g", conv[order(-last_change)][1:5, block], conv[order(-last_change)][1:5, last_change]), collapse = ", ")), "",
  "Validation report: `Rscript calibration/v020-soil/validation.R`. Convergence: judge by the parameter drift between checkpoints, not by the loss."
), file.path(FIT_DIR, "fit_summary.md"))
message(sprintf("done: %.1f min (%.1f s/epoch); final losses: %s", min_fit, 60 * min_fit / EPOCHS,
                paste(sprintf("%s=%.3f", names(hist)[1:6], as.numeric(hist[.N, 1:6])), collapse = " ")))
