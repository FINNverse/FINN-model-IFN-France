#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-IFN-France | calibration/validation.R
##
## WHAT   The validation report of the model fitted by calibration/fit.R
##        (fiveyear-v1-process). Runs on the saved model, so it can be redone
##        without refitting. Two sets, both simulated from their own first-visit
##        stands over the 5-year interval and aggregated as fit() scores them
##        (state at year 5; growth and mortality = mean of the annual rates over
##        years 1..5; recruitment = sum):
##          TEST    = fiveyear-v1/train, the sites the model was fitted on
##          HOLDOUT = fiveyear-v1/holdout_blocked, 100 km blocks never seen
##        The model has no coordinates, so the holdout is predicted as is.
##        Sections: metrics per response and set (species-level rows pooled, and
##        stand totals), per-species R2, mortality AUC (weighted by the tallied
##        counts), persistence of the first-visit stand as the model-free
##        reference for the states; figures: loss curves, observed vs predicted,
##        ROC, map, bare-ground succession across the climate space, parameter
##        trajectories; provenance table.
## RUN    Rscript calibration/validation.R      (from the repository root, ~20 min on one core)
## OUTPUT calibration/report.md, calibration/validation/{metrics.csv, metrics_species.csv,
##        mortality_auc.csv, bare_ground.csv, parameter_change.csv}, calibration/figures/*.png;
##        <DATA_DIR>/fits/fiveyear-v1-process/validation/predictions.csv (large);
##        copies of metrics.csv, species.csv and the figures into the zoo folder if present
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-0848b60"), .libPaths()))   # FINN 0.2.0, main commit 0848b60
suppressPackageStartupMessages({ library(FINN); library(data.table); library(torch); library(ggplot2) })

DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "frnfi")
FIT_DIR  <- file.path(DATA_DIR, "fits", "fiveyear-v1-process")
CAL_DIR  <- "calibration"
VAL      <- file.path(CAL_DIR, "validation"); FIG <- file.path(CAL_DIR, "figures")
ZOO_DIR  <- file.path("..", "..", "FINN-zoo-worktrees", "frnfi", "IGN-FR_process")   # the model card's folder, if checked out
for (d in c(VAL, FIG, file.path(FIT_DIR, "validation"))) dir.create(d, showWarnings = FALSE, recursive = TRUE)

PARAMS <- list(
  prep_version   = "fiveyear-v1",
  sets           = c(test = "train", holdout = "holdout_blocked"),
  patch_size_ha  = pi * 15^2 / 1e4, n_patches = 1L,       # as fit.R
  env_vars       = c("tavg_c", "prec_mm"),
  period         = 5L,
  sim_batchsize  = 2000L,                                  # sites per simulateForest batch
  bare_locations = 5L,                                     # climates for the bare-ground runs: quantiles of the training tavg_c
  bare_quantiles = c(0.05, 0.275, 0.5, 0.725, 0.95),
  bare_years     = 500L, bare_patches = 15L,
  n_top_species  = 6L,                                     # coloured in the succession figure; the rest pale grey
  threads        = 1L,                                     # forward simulations only; one core keeps a shared machine usable
  seed           = 1L
)
P <- PARAMS$period
torch_set_num_threads(PARAMS$threads)

## ---- helpers ---------------------------------------------------------------------
metric <- function(o, p) {
  ok <- is.finite(o) & is.finite(p); o <- o[ok]; p <- p[ok]
  if (!length(o)) return(data.table(n = 0L, rmse = NA_real_, r2 = NA_real_, r = NA_real_, bias = NA_real_, obs_mean = NA_real_))
  sst <- sum((o - mean(o))^2)
  data.table(n = length(o), rmse = sqrt(mean((p - o)^2)), r2 = if (sst > 0) 1 - sum((p - o)^2) / sst else NA_real_,
             r = if (sd(o) > 0 && sd(p) > 0) cor(o, p) else NA_real_, bias = mean(p) - mean(o), obs_mean = mean(o))
}
## AUC for grouped counts: each row contributes n_pos deaths and n_neg survivors at one predicted probability; ties count half
weighted_auc <- function(p, n_pos, n_neg) {
  dt <- data.table(p = p, n_pos = n_pos, n_neg = n_neg)[is.finite(p) & (n_pos + n_neg) > 0][, .(n_pos = sum(n_pos), n_neg = sum(n_neg)), by = p][order(p)]
  if (!nrow(dt) || sum(dt$n_pos) == 0 || sum(dt$n_neg) == 0) return(NA_real_)
  before <- cumsum(dt$n_neg) - dt$n_neg
  sum(dt$n_pos * (before + 0.5 * dt$n_neg)) / (sum(dt$n_pos) * sum(dt$n_neg))
}
roc_points <- function(p, n_pos, n_neg) {
  dt <- data.table(p = p, n_pos = n_pos, n_neg = n_neg)[is.finite(p) & (n_pos + n_neg) > 0][, .(n_pos = sum(n_pos), n_neg = sum(n_neg)), by = p][order(-p)]
  data.table(fpr = c(0, cumsum(dt$n_neg) / sum(dt$n_neg)), tpr = c(0, cumsum(dt$n_pos) / sum(dt$n_pos)))
}
md_table <- function(dt) {
  dt <- copy(dt); num <- names(dt)[vapply(dt, is.numeric, logical(1))]
  dt[, (num) := lapply(.SD, signif, 3), .SDcols = num]
  c(paste0("| ", paste(names(dt), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(dt)), collapse = "|"), "|"),
    apply(dt, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
read_set <- function(set) {
  d <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, set)
  list(obs = fread(file.path(d, "obs_dt.csv")), env = fread(file.path(d, "env_dt.csv")), init = fread(file.path(d, "init_cohort_dt.csv")),
       species = fread(file.path(d, "species_dt.csv")), sites = fread(file.path(d, "siteID_dt.csv")))
}

## ---- load ------------------------------------------------------------------------
model <- torch_load(file.path(FIT_DIR, "model.pt")); model$eval()
loss_hist <- fread(file.path(FIT_DIR, "loss_history.csv"))
test <- read_set(PARAMS$sets[["test"]]); hold <- read_set(PARAMS$sets[["holdout"]])
species <- test$species[order(species)]
n_species <- nrow(species)
message(sprintf("test %d sites | holdout %d sites | %d species", uniqueN(test$obs$siteID), uniqueN(hold$obs$siteID), n_species))

## ---- 1. site-level predictions, aggregated as fit() scores them ----------------
predict_set <- function(d) {
  env <- d$env[, c("siteID", "year", PARAMS$env_vars), with = FALSE]
  ids <- sort(unique(d$obs$siteID)); chunks <- split(ids, ceiling(seq_along(ids) / PARAMS$sim_batchsize))
  S <- rbindlist(lapply(chunks, function(ch) {
    FINN.seed(PARAMS$seed)
    sim <- simulateForest(model, env = env[siteID %in% ch],
                          init_cohort = makeInitCohorts(d$init[siteID %in% ch, .(siteID, patchID, species, dbh, trees, treeName)], Nspecies = n_species),
                          patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, device = "cpu")
    W <- as.data.table(sim$wide$site)
    for (cc in c("siteID", "year", "species")) W[, (cc) := as.integer(as.character(get(cc)))]
    W
  }))
  ## state at year 5; mean annual growth / mortality and summed recruitment over years 1..5
  A <- S[year >= 1L & year <= P, .(growth_m = mean(growth), mort_m = mean(mort), reg_m = sum(r_mean_ha)), by = .(siteID, species)]
  St <- S[year == P, .(siteID, species, dbh_m = dbh, ba_m = ba, trees_m = trees)]
  Pr <- merge(St, A, by = c("siteID", "species"))
  M <- merge(d$obs, Pr, by = c("siteID", "species"))
  M[, mort_period_m := 1 - (1 - pmin(pmax(mort_m, 0), 1))^P]
  ## persistence: the first-visit stand as the model-free reference for the states
  I <- d$init[, .(ba_0 = sum(dbh2ba(dbh) * trees), trees_0 = sum(trees), dbh_0 = weighted.mean(dbh, trees)), by = .(siteID, species)]
  M <- merge(M, I, by = c("siteID", "species"), all.x = TRUE)
  M[is.na(ba_0), `:=`(ba_0 = 0, trees_0 = 0)]
  M <- merge(M, d$sites[, .(siteID, x, y)], by = "siteID")
  M
}
M <- rbind(cbind(set = "test", predict_set(test)), cbind(set = "holdout", predict_set(hold)))
M[, set := factor(set, levels = c("test", "holdout"))]
M <- merge(M, species, by = "species", suffixes = c("", "_sp"))
fwrite(M, file.path(FIT_DIR, "validation", "predictions.csv"))

## ---- 2. metrics ------------------------------------------------------------------
vars <- list(
  list("dbh (cm, species mean)",              quote(trees > 0 & state_masked == FALSE), "dbh",    "dbh_m"),
  list("basal area (m2 per 0.0707 ha patch)", quote(state_masked == FALSE),             "ba",     "ba_m"),
  list("trees (per 0.0707 ha patch)",         quote(state_masked == FALSE),             "trees",  "trees_m"),
  list("growth (rel. dbh increment / yr)",    quote(growth_n > 0),                      "growth", "growth_m"),
  list("mortality (annual rate)",             quote(n_at_risk > 0),                     "mort",   "mort_m"),
  list("recruitment (stems / ha / 5 yr)",     quote(rep(TRUE, .N)),                     "reg",    "reg_m"))
metrics <- rbindlist(lapply(levels(M$set), function(s) rbindlist(lapply(vars, function(v) {
  d <- M[set == s][eval(v[[2]])]
  cbind(set = s, variable = v[[1]], level = "species x site", metric(d[[v[[3]]]], d[[v[[4]]]]))
}))))
## stand totals over species (basal area, stems, recruitment), unmasked sites
TOT <- M[state_masked == FALSE, .(ba = sum(ba), ba_m = sum(ba_m), ba_0 = sum(ba_0), trees = sum(trees), trees_m = sum(trees_m), trees_0 = sum(trees_0),
                                  reg = sum(reg), reg_m = sum(reg_m), x = x[1], y = y[1]), by = .(set, siteID)]
metrics <- rbind(metrics, rbindlist(lapply(levels(M$set), function(s) rbind(
  cbind(set = s, variable = "basal area (m2 per patch)", level = "stand total", metric(TOT[set == s]$ba, TOT[set == s]$ba_m)),
  cbind(set = s, variable = "trees (per patch)",         level = "stand total", metric(TOT[set == s]$trees, TOT[set == s]$trees_m)),
  cbind(set = s, variable = "recruitment (stems / ha / 5 yr)", level = "stand total", metric(TOT[set == s]$reg, TOT[set == s]$reg_m)),
  cbind(set = s, variable = "basal area (m2 per patch)", level = "stand total, PERSISTENCE reference", metric(TOT[set == s]$ba, TOT[set == s]$ba_0)),
  cbind(set = s, variable = "trees (per patch)",         level = "stand total, PERSISTENCE reference", metric(TOT[set == s]$trees, TOT[set == s]$trees_0))))))
auc_tab <- rbindlist(lapply(levels(M$set), function(s) {
  d <- M[set == s & n_at_risk > 0]
  data.table(set = s, n_trees_at_risk = sum(d$n_at_risk), n_died = sum(d$n_died),
             auc = weighted_auc(d$mort_period_m, d$n_died, d$n_at_risk - d$n_died),
             obs_death_share = sum(d$n_died) / sum(d$n_at_risk), pred_death_share = sum(d$n_at_risk * d$mort_period_m) / sum(d$n_at_risk))
}))
metrics[variable == "mortality (annual rate)", auc := auc_tab$auc[match(set, auc_tab$set)]]
fwrite(metrics, file.path(VAL, "metrics.csv")); fwrite(auc_tab, file.path(VAL, "mortality_auc.csv"))
print(metrics, digits = 3); print(auc_tab, digits = 3)
## per species (holdout and test): basal area, stems, growth, mortality AUC, recruitment
SP <- M[, {
  u <- state_masked == FALSE
  .(sites_present = sum(trees > 0, na.rm = TRUE),
    ba_r2 = metric(ba[u], ba_m[u])$r2, ba_rmse = metric(ba[u], ba_m[u])$rmse,
    trees_r2 = metric(trees[u], trees_m[u])$r2,
    growth_r = metric(growth[growth_n > 0], growth_m[growth_n > 0])$r, growth_bias = metric(growth[growth_n > 0], growth_m[growth_n > 0])$bias,
    mort_auc = weighted_auc(mort_period_m[n_at_risk > 0], n_died[n_at_risk > 0], (n_at_risk - n_died)[n_at_risk > 0]),
    mort_obs = sum(n_died) / sum(n_at_risk), mort_pred = sum(n_at_risk * mort_period_m) / sum(n_at_risk),
    reg_r = metric(reg, reg_m)$r, reg_obs = mean(reg), reg_pred = mean(reg_m))
}, by = .(set, species, species_name)][order(set, species)]
fwrite(SP, file.path(VAL, "metrics_species.csv"))

## ---- 3. figures -------------------------------------------------------------------
## fig1 loss per response over epochs
loss_cols <- setdiff(names(loss_hist), c("epoch", "unused"))
LH <- melt(loss_hist[, c("epoch", loss_cols), with = FALSE][, (loss_cols) := lapply(.SD, as.numeric), .SDcols = loss_cols],
           id.vars = "epoch", variable.name = "response", value.name = "loss")[is.finite(loss)]
g1 <- ggplot(LH, aes(epoch, loss)) + geom_line() + facet_wrap(~ response, scales = "free_y") + theme_bw() +
  labs(title = "Training loss per response (baseline-scaled)", y = "loss")
ggsave(file.path(FIG, "fig1_loss_epochs.png"), g1, width = 9, height = 5, dpi = 150, bg = "white")
## fig2 observed vs predicted, every response, both sets
OP <- rbind(
  M[state_masked == FALSE & trees > 0, .(set, variable = "dbh (cm)", obs = dbh, pred = dbh_m)],
  TOT[, .(set, variable = "basal area, stand (m2/ha)", obs = ba / PARAMS$patch_size_ha, pred = ba_m / PARAMS$patch_size_ha)],
  TOT[, .(set, variable = "trees, stand (per ha)", obs = trees / PARAMS$patch_size_ha, pred = trees_m / PARAMS$patch_size_ha)],
  M[growth_n > 0, .(set, variable = "growth (rel. / yr)", obs = growth, pred = growth_m)],
  M[n_at_risk > 0, .(set, variable = "mortality (annual rate)", obs = mort, pred = mort_m)],
  TOT[, .(set, variable = "recruitment, stand (per ha / 5 yr)", obs = reg, pred = reg_m)])
OP[, panel := factor(paste(set, variable, sep = ": "), levels = unique(paste(rep(levels(M$set), each = 6), unique(variable), sep = ": ")))]
g2 <- ggplot(OP[is.finite(obs) & is.finite(pred)], aes(obs, pred)) + geom_hex(bins = 60) + geom_abline(colour = "red") +
  facet_wrap(~ panel, scales = "free", nrow = 2) + scale_fill_viridis_c(trans = "log10") + theme_bw() +
  labs(title = "Observed vs predicted at year 5", x = "observed", y = "predicted")
ggsave(file.path(FIG, "fig2_obs_pred.png"), g2, width = 15, height = 6, dpi = 150, bg = "white")
## fig3 ROC
ROC <- M[n_at_risk > 0, roc_points(mort_period_m, n_died, n_at_risk - n_died), by = set]
g3 <- ggplot(ROC, aes(fpr, tpr, colour = set)) + geom_path() + geom_abline(linetype = 2) + coord_equal() + theme_bw() +
  labs(title = "Mortality ROC (5-year death probability vs tallied fates)", subtitle = paste(auc_tab[, sprintf("%s AUC %.3f", set, auc)], collapse = " | "),
       x = "false positive rate", y = "true positive rate")
ggsave(file.path(FIG, "fig3_mortality_roc.png"), g3, width = 6, height = 6, dpi = 150, bg = "white")
## fig4 map: observed and predicted stand basal area, both sets
MAP <- melt(TOT[, .(set, x, y, observed = ba / PARAMS$patch_size_ha, predicted = ba_m / PARAMS$patch_size_ha)], id.vars = c("set", "x", "y"), variable.name = "what", value.name = "ba_ha")
g4 <- ggplot(MAP, aes(x, y, colour = pmin(ba_ha, 60))) + geom_point(size = 0.3) + facet_grid(what ~ set) + coord_quickmap() +
  scale_colour_viridis_c(name = "BA m2/ha\n(capped 60)") + theme_bw() + labs(title = "Stand basal area at year 5: test (fitted blocks) and holdout blocks", x = "longitude", y = "latitude")
ggsave(file.path(FIG, "fig4_map.png"), g4, width = 9, height = 9, dpi = 150, bg = "white")
## fig5 bare-ground succession at five climates (training tavg_c quantiles, with the matching median prec_mm)
TE <- unique(test$env[, c("siteID", PARAMS$env_vars), with = FALSE])
qs <- quantile(TE$tavg_c, PARAMS$bare_quantiles)
loc <- rbindlist(lapply(seq_along(qs), function(i) {
  nb <- TE[abs(tavg_c - qs[i]) <= 0.25]
  data.table(siteID = i, tavg_c = qs[[i]], prec_mm = median(nb$prec_mm))
}))
env_bare <- loc[, .(year = seq_len(PARAMS$bare_years)), by = .(siteID, tavg_c, prec_mm)][, .(siteID, year, tavg_c, prec_mm)]
FINN.seed(PARAMS$seed)
sim_bare <- simulateForest(model, env = env_bare, init_cohort = NULL, patches = PARAMS$bare_patches, patch_size = PARAMS$patch_size_ha, device = "cpu")
B <- as.data.table(sim_bare$wide$site)
for (cc in c("siteID", "year", "species")) B[, (cc) := as.integer(as.character(get(cc)))]
B <- merge(B, species, by = "species"); B <- merge(B, loc[, .(siteID, label = sprintf("%.1f degC, %.0f mm", tavg_c, prec_mm))], by = "siteID")
B[, label := factor(label, levels = loc[order(tavg_c), sprintf("%.1f degC, %.0f mm", tavg_c, prec_mm)])]   # panels from cold to warm
B[, ba_ha := ba / PARAMS$patch_size_ha]
fwrite(B[, .(siteID, label, year, species, species_name, ba_ha, trees_ha = trees / PARAMS$patch_size_ha)], file.path(VAL, "bare_ground.csv"))
top <- B[, .(m = max(ba_ha)), by = species_name][order(-m)][seq_len(PARAMS$n_top_species), species_name]
B[, shown := fifelse(species_name %in% top, species_name, "other species")]
g5 <- ggplot() +
  geom_line(data = B[shown == "other species"], aes(year, ba_ha, group = species_name), colour = "grey80", linewidth = 0.3) +
  geom_line(data = B[shown != "other species"], aes(year, ba_ha, colour = shown), linewidth = 0.6) +
  facet_wrap(~ label, nrow = 1) + theme_bw() + labs(title = sprintf("Succession from bare ground, %d patches of 0.0707 ha, constant climate", PARAMS$bare_patches),
                                                     y = "basal area (m2/ha)", colour = "species")
ggsave(file.path(FIG, "fig5_bare_ground.png"), g5, width = 14, height = 4.5, dpi = 150, bg = "white")
EQ <- B[, .(ba_ha = sum(ba_ha), trees_ha = sum(trees) / PARAMS$patch_size_ha), by = .(label, year)]
eq_tab <- EQ[, .(ba_last100 = mean(ba_ha[year > max(year) - 100]), ba_prev100 = mean(ba_ha[year > max(year) - 200 & year <= max(year) - 100]),
                 trees_last100 = mean(trees_ha[year > max(year) - 100])), by = label][, change_pct := 100 * (ba_last100 / ba_prev100 - 1)]
## fig6 parameter trajectories (unconstrained species parameters, one line per species) and settling
PH <- model$param_history
if (length(PH) > 1) {
  sp_blocks <- grep("^par_|^reg_logK", names(PH[[1]]), value = TRUE)
  TR <- rbindlist(lapply(seq_along(PH), function(k) rbindlist(lapply(sp_blocks, function(b) {
    m <- as.matrix(PH[[k]][[b]]); if (is.null(dim(m)) || ncol(m) == 1) m <- matrix(m, ncol = 1)
    data.table(snapshot = k, block = b, species = rep(seq_len(nrow(m)), ncol(m)), column = rep(seq_len(ncol(m)), each = nrow(m)), value = as.vector(m))
  }))))
  TR[, panel := paste0(block, "[", column, "]")]
  g6a <- ggplot(TR, aes(snapshot, value, group = species, colour = factor(species))) + geom_line(linewidth = 0.3, show.legend = FALSE) +
    facet_wrap(~ panel, scales = "free_y") + theme_bw() + labs(title = "Species parameters over the fit (unconstrained scale, one line per species)", x = "snapshot (every 10 epochs)")
  ggsave(file.path(FIG, "fig6_parameter_trajectories.png"), g6a, width = 12, height = 7, dpi = 150, bg = "white")
  CH <- rbindlist(lapply(2:length(PH), function(k) rbindlist(lapply(names(PH[[k]]), function(b)
    data.table(snapshot = k, block = b, mean_abs_change = mean(abs(as.matrix(PH[[k]][[b]]) - as.matrix(PH[[k - 1]][[b]]))))))))
  fwrite(CH, file.path(VAL, "parameter_change.csv"))
  g6b <- ggplot(CH[mean_abs_change > 0], aes(snapshot, mean_abs_change, colour = block)) + geom_line() + scale_y_log10() + theme_bw() +
    labs(title = "Mean absolute parameter change between snapshots", x = "snapshot (every 10 epochs)")
  ggsave(file.path(FIG, "fig6b_parameter_change.png"), g6b, width = 9, height = 5, dpi = 150, bg = "white")
}

## ---- 4. report.md -------------------------------------------------------------------
fit_summary <- if (file.exists(file.path(CAL_DIR, "fit_summary.md"))) readLines(file.path(CAL_DIR, "fit_summary.md")) else character(0)
sha <- if (file.exists(file.path(FIT_DIR, "parameters.pt"))) unname(tools::sha256sum(file.path(FIT_DIR, "parameters.pt"))) else NA
hold_sp <- SP[set == "holdout"][order(-sites_present)]
report <- c(
  "# Validation report: fiveyear-v1-process", "",
  sprintf("Written %s by calibration/validation.R; FINN %s (main 0848b60), torch %s. Model: %s/model.pt; release object parameters.pt sha256 %s.",
          format(Sys.time(), "%Y-%m-%d %H:%M %Z"), packageVersion("FINN"), packageVersion("torch"), FIT_DIR, sha), "",
  "## Sets", "",
  sprintf("- TEST = fiveyear-v1/train: %d sites the model was fitted on (%d with the year-5 state masked because of a cut tree).", uniqueN(test$obs$siteID), test$obs[state_masked == TRUE, uniqueN(siteID)]),
  sprintf("- HOLDOUT = fiveyear-v1/holdout_blocked: %d sites in 100 km blocks the model never saw (%d state-masked). No coordinates in the model, predicted as is.", uniqueN(hold$obs$siteID), hold$obs[state_masked == TRUE, uniqueN(siteID)]),
  "- Aggregation as fit() scores it: state at year 5; growth and mortality = mean of the annual rates over years 1..5; recruitment = sum over years 1..5.", "",
  "## Metrics", "", md_table(metrics[, .(set, variable, level, n, rmse, r2, r, bias, obs_mean, auc)]), "",
  "A negative R2 means worse than predicting the mean. The persistence rows predict the year-5 state by the first-visit state (model-free reference).", "",
  "## Mortality AUC (5-year death probability against the tallied fates; each site x species row weighted by its trees at risk)", "", md_table(auc_tab), "",
  "## Per species, holdout", "", md_table(hold_sp[, .(species_name, sites_present, ba_r2, trees_r2, growth_r, growth_bias, mort_auc, mort_obs, mort_pred, reg_r, reg_obs, reg_pred)]), "",
  "## Bare-ground equilibrium (stand totals, last 100 vs previous 100 years)", "", md_table(eq_tab), "",
  "## Figures", "",
  "- figures/fig1_loss_epochs.png: loss per response over epochs", "- figures/fig2_obs_pred.png: observed vs predicted, both sets",
  "- figures/fig3_mortality_roc.png: mortality ROC", "- figures/fig4_map.png: observed and predicted stand basal area, both sets",
  "- figures/fig5_bare_ground.png: succession from bare ground at five climates", "- figures/fig6_parameter_trajectories.png, fig6b_parameter_change.png: parameter trajectories and settling", "",
  "## Provenance", "", "| step | script | inputs | outputs |", "|---|---|---|---|",
  "| raw data | calibration/01_download_ign_raw.R | IGN DataIFN export 2005-2025 (sha256 in PROVENANCE.md) | raw/ign_data/unz, PROVENANCE.md |",
  "| preparation | calibration/02_prepare_fiveyear_v1.R | raw tables, WorldClim 2.1 30 s | prepared/fiveyear-v1/{national,train,holdout_blocked,pilot}, PREPARATION.md, expected_preparation.csv |",
  "| fit | calibration/fit.R | fiveyear-v1/train | fits/fiveyear-v1-process/{model.pt, parameters.pt}, loss_history.csv, fit_summary.md |",
  "| validation | calibration/validation.R | model.pt, train, holdout_blocked | report.md, validation/*.csv, figures/*.png |", "",
  "## Fit summary (from fit.R)", "", fit_summary)
writeLines(report, file.path(CAL_DIR, "report.md"))

## ---- 5. the zoo folder: metrics.csv, species.csv, figures --------------------------
if (dir.exists(ZOO_DIR)) {
  dir.create(file.path(ZOO_DIR, "figures"), showWarnings = FALSE)
  fwrite(metrics, file.path(ZOO_DIR, "metrics.csv"))
  file.copy(file.path(DATA_DIR, "prepared", PARAMS$prep_version, "species.csv"), file.path(ZOO_DIR, "species.csv"), overwrite = TRUE)
  for (f in list.files(FIG, pattern = "\\.png$", full.names = TRUE)) file.copy(f, file.path(ZOO_DIR, "figures", basename(f)), overwrite = TRUE)
  message("copied metrics, species and figures to ", ZOO_DIR)
}
message("done: ", file.path(CAL_DIR, "report.md"))
