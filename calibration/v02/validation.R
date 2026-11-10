#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-IFN-France | calibration/v02/validation.R
##
## WHAT   The validation report of the v0.2 production fit (calibration/v02/fit.R,
##        fits/fiveyear-v2-process). Copy of calibration/validation.R (v0.1) with:
##        - data fiveyear-v2 (train / holdout_blocked), covariates tavg_c, prec_mm, gest;
##        - the mortality target is TOTAL REMOVAL; the report shows the removal
##          level on ONE scale (annual observed vs annual predicted; the v0.1
##          report compared annual observed deaths with the 5-year predicted
##          probability, which read as a 3.6-fold over-prediction where the
##          fitted rate was in fact 0.78x), plus the natural-only share for
##          reference;
##        - the fitted per-tree 5-year removal probability from FINN::predictTrees
##          against the tallied fates (not subject to the simulation's 1 %/yr
##          mortality floor), next to the simulated stand;
##        - change metrics (dBA, dN) against persistence;
##        - bare-ground runs at gest 0 (no management trace) and gest 2.
##        TEST = fiveyear-v2/train (fitted), HOLDOUT = fiveyear-v2/holdout_blocked.
## RUN    Rscript calibration/v02/validation.R   (repository root, one thread, ~25 min)
## OUTPUT calibration/v02/report.md, calibration/v02/validation/*.csv, calibration/v02/figures/*.png;
##        <DATA_DIR>/fits/fiveyear-v2-process/validation/predictions.csv; copies of
##        metrics.csv and the figures into the zoo card folder if present
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-0848b60"), .libPaths()))
suppressPackageStartupMessages({ library(FINN); library(data.table); library(torch); library(ggplot2) })

DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "frnfi")
FIT_DIR  <- file.path(DATA_DIR, "fits", "fiveyear-v2-process")
CAL_DIR  <- file.path("calibration", "v02")
VAL      <- file.path(CAL_DIR, "validation"); FIG <- file.path(CAL_DIR, "figures")
ZOO_DIR  <- file.path("..", "..", "FINN-zoo-worktrees", "frnfi", "IGN-FR_process")
for (d in c(VAL, FIG, file.path(FIT_DIR, "validation"))) dir.create(d, showWarnings = FALSE, recursive = TRUE)

PARAMS <- list(
  prep_version   = "fiveyear-v2",
  sets           = c(test = "train", holdout = "holdout_blocked"),
  patch_size_ha  = pi * 15^2 / 1e4, n_patches = 1L,
  env_vars       = c("tavg_c", "prec_mm", "gest"),
  period         = 5L,
  sim_batchsize  = 2000L,
  bare_quantiles = c(0.05, 0.275, 0.5, 0.725, 0.95),  # tavg_c quantiles of the training sites
  bare_gest      = c(0, 2),                             # bare ground without and with clear management traces
  bare_years     = 500L, bare_patches = 15L,
  n_top_species  = 6L,
  threads        = 1L,
  seed           = 1L
)
P <- PARAMS$period
torch_set_num_threads(PARAMS$threads)

## ---- helpers ------------------------------------------------------------------------
metric <- function(o, p) {
  ok <- is.finite(o) & is.finite(p); o <- o[ok]; p <- p[ok]
  if (!length(o)) return(data.table(n = 0L, rmse = NA_real_, r2 = NA_real_, r = NA_real_, bias = NA_real_, obs_mean = NA_real_))
  sst <- sum((o - mean(o))^2)
  data.table(n = length(o), rmse = sqrt(mean((p - o)^2)), r2 = if (sst > 0) 1 - sum((p - o)^2) / sst else NA_real_,
             r = if (sd(o) > 0 && sd(p) > 0) cor(o, p) else NA_real_, bias = mean(p) - mean(o), obs_mean = mean(o))
}
weighted_auc <- function(p, n_pos, n_neg) {
  dt <- data.table(p = p, n_pos = n_pos, n_neg = n_neg)[is.finite(p) & (n_pos + n_neg) > 0][, .(n_pos = sum(n_pos), n_neg = sum(n_neg)), by = p][order(p)]
  if (!nrow(dt) || sum(dt$n_pos) == 0 || sum(dt$n_neg) == 0) return(NA_real_)
  before <- cumsum(dt$n_neg) - dt$n_neg; sum(dt$n_pos * (before + 0.5 * dt$n_neg)) / (sum(dt$n_pos) * sum(dt$n_neg))
}
roc_points <- function(p, n_pos, n_neg) {
  dt <- data.table(p = p, n_pos = n_pos, n_neg = n_neg)[is.finite(p) & (n_pos + n_neg) > 0][, .(n_pos = sum(n_pos), n_neg = sum(n_neg)), by = p][order(-p)]
  data.table(fpr = c(0, cumsum(dt$n_neg) / sum(dt$n_neg)), tpr = c(0, cumsum(dt$n_pos) / sum(dt$n_pos)))
}
md_table <- function(dt) {
  dt <- copy(dt); num <- names(dt)[vapply(dt, is.numeric, logical(1))]; dt[, (num) := lapply(.SD, signif, 3), .SDcols = num]
  c(paste0("| ", paste(names(dt), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(dt)), collapse = "|"), "|"), apply(dt, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
read_set <- function(set) {
  d <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, set)
  list(obs = fread(file.path(d, "obs_dt.csv")), env = fread(file.path(d, "env_dt.csv")), init = fread(file.path(d, "init_cohort_dt.csv")),
       trees = fread(file.path(d, "tree_obs_dt.csv")), species = fread(file.path(d, "species_dt.csv")), sites = fread(file.path(d, "siteID_dt.csv")))
}

## ---- load ----------------------------------------------------------------------------
model <- torch_load(file.path(FIT_DIR, "model.pt")); model$eval()
loss_hist <- fread(file.path(FIT_DIR, "loss_history.csv"))
test <- read_set(PARAMS$sets[["test"]]); hold <- read_set(PARAMS$sets[["holdout"]])
species <- test$species[order(species)]; n_species <- nrow(species)
message(sprintf("test %d sites | holdout %d sites | %d species classes", uniqueN(test$obs$siteID), uniqueN(hold$obs$siteID), n_species))

## ---- 1. site-level predictions, aggregated as fit() scores them --------------------------
predict_set <- function(d) {
  env <- d$env[, c("siteID", "year", PARAMS$env_vars), with = FALSE]
  ids <- sort(unique(d$obs$siteID)); chunks <- split(ids, ceiling(seq_along(ids) / PARAMS$sim_batchsize))
  S <- rbindlist(lapply(chunks, function(ch) {
    FINN.seed(PARAMS$seed)
    sim <- simulateForest(model, env = env[siteID %in% ch],
                          init_cohort = makeInitCohorts(d$init[siteID %in% ch, .(siteID, patchID, species, dbh, trees, treeName)], Nspecies = n_species),
                          patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, device = "cpu")
    W <- as.data.table(sim$wide$site); for (cc in c("siteID", "year", "species")) W[, (cc) := as.integer(as.character(get(cc)))]; W
  }))
  A <- S[year >= 1L & year <= P, .(growth_m = mean(growth), mort_m = mean(mort), reg_m = sum(r_mean_ha)), by = .(siteID, species)]
  M <- merge(d$obs, merge(S[year == P, .(siteID, species, dbh_m = dbh, ba_m = ba, trees_m = trees)], A, by = c("siteID", "species")), by = c("siteID", "species"))
  M[, mort_period_m := 1 - (1 - pmin(pmax(mort_m, 0), 1))^P]
  I <- d$init[, .(ba_0 = sum(dbh2ba(dbh) * trees), trees_0 = sum(trees), dbh_0 = weighted.mean(dbh, trees)), by = .(siteID, species)]
  M <- merge(M, I, by = c("siteID", "species"), all.x = TRUE); M[is.na(ba_0), `:=`(ba_0 = 0, trees_0 = 0)]
  merge(M, d$sites[, .(siteID, x, y, gest, end_stratum, end_share_out)], by = "siteID")
}
M <- rbind(cbind(set = "test", predict_set(test)), cbind(set = "holdout", predict_set(hold)))
M[, set := factor(set, levels = c("test", "holdout"))]
M <- merge(M, species, by = "species", suffixes = c("", "_sp"))
fwrite(M, file.path(FIT_DIR, "validation", "predictions.csv"))

## ---- 1b. per-tree removal probability (not floored) ---------------------------------------
tree_pred <- function(d) {
  td <- d$trees[, .(siteID, year, patchID, species, dbh, trees, growth, died, period_length)]
  pt <- predictTrees(model, tree_data = td, env = d$env[, c("siteID", "year", PARAMS$env_vars), with = FALSE], patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, device = "cpu")
  stopifnot(nrow(pt) == nrow(td))
  cbind(d$trees[, .(siteID, species, species_name, dbh, trees, growth, died, died_natural, cut)], pt[, .(growth_pred, mort_pred)])
}
TR <- rbind(cbind(set = "test", tree_pred(test)), cbind(set = "holdout", tree_pred(hold)))
TR[, set := factor(set, levels = c("test", "holdout"))]
TR[, dcls := cut(dbh, c(0, 15, 25, 40, 60, Inf))]
fwrite(TR[, .(set, siteID, species, dbh, trees, growth, died, died_natural, cut, growth_pred, mort_pred)], file.path(FIT_DIR, "validation", "tree_predictions.csv"))

## ---- 2. metrics ------------------------------------------------------------------------------
vars <- list(
  list("dbh (cm, species mean)",              quote(trees > 0),     "dbh",    "dbh_m"),
  list("basal area (m2 per 0.0707 ha patch)", quote(rep(TRUE, .N)), "ba",     "ba_m"),
  list("trees (per 0.0707 ha patch)",         quote(rep(TRUE, .N)), "trees",  "trees_m"),
  list("growth (rel. dbh increment / yr)",    quote(growth_n > 0),  "growth", "growth_m"),
  list("removal (annual rate, dead + cut)",   quote(n_at_risk > 0), "mort",   "mort_m"),
  list("recruitment (stems / ha / 5 yr)",     quote(rep(TRUE, .N)), "reg",    "reg_m"))
metrics <- rbindlist(lapply(levels(M$set), function(s) rbindlist(lapply(vars, function(v) {
  d <- M[set == s][eval(v[[2]])]; cbind(set = s, variable = v[[1]], level = "species x site", metric(d[[v[[3]]]], d[[v[[4]]]])) }))))
TOT <- M[, .(ba = sum(ba), ba_m = sum(ba_m), ba_0 = sum(ba_0), trees = sum(trees), trees_m = sum(trees_m), trees_0 = sum(trees_0), reg = sum(reg), reg_m = sum(reg_m), x = x[1], y = y[1], gest = gest[1]), by = .(set, siteID)]
metrics <- rbind(metrics, rbindlist(lapply(levels(M$set), function(s) rbind(
  cbind(set = s, variable = "basal area (m2 per patch)", level = "stand total", metric(TOT[set == s]$ba, TOT[set == s]$ba_m)),
  cbind(set = s, variable = "trees (per patch)",         level = "stand total", metric(TOT[set == s]$trees, TOT[set == s]$trees_m)),
  cbind(set = s, variable = "recruitment (stems / ha / 5 yr)", level = "stand total", metric(TOT[set == s]$reg, TOT[set == s]$reg_m)),
  cbind(set = s, variable = "basal area (m2 per patch)", level = "stand total, PERSISTENCE reference", metric(TOT[set == s]$ba, TOT[set == s]$ba_0)),
  cbind(set = s, variable = "trees (per patch)",         level = "stand total, PERSISTENCE reference", metric(TOT[set == s]$trees, TOT[set == s]$trees_0)),
  cbind(set = s, variable = "5-yr change of basal area (m2 per patch)", level = "stand total", metric(TOT[set == s, ba - ba_0], TOT[set == s, ba_m - ba_0])),
  cbind(set = s, variable = "5-yr change of trees (per patch)",         level = "stand total", metric(TOT[set == s, trees - trees_0], TOT[set == s, trees_m - trees_0]))))))
## removal: one scale (annual), level ratio, AUC on the simulated rate and on the fitted per-tree probability
rem <- rbindlist(lapply(levels(M$set), function(s) {
  d <- M[set == s & n_at_risk > 0]; t <- TR[set == s & is.finite(mort_pred)]
  data.table(set = s, n_trees_at_risk = sum(d$n_at_risk),
             obs_annual_removal = sum(d$n_died) / sum(d$n_at_risk), pred_annual_removal = sum(d$n_at_risk * d$mort_m) / sum(d$n_at_risk),
             obs_annual_natural = sum(d$n_died_natural) / sum(d$n_at_risk),
             auc_site_rate = weighted_auc(d$mort_period_m, d$n_died, d$n_at_risk - d$n_died),
             tree_obs_5yr_removal = mean(t$died), tree_pred_5yr_removal = mean(t$mort_pred), auc_per_tree = weighted_auc(t$mort_pred, t$died, 1 - t$died),
             share_trees_below_floor = mean(1 - (1 - t$mort_pred)^(1 / P) < 0.01))
}))[, `:=`(ratio_annual = pred_annual_removal / obs_annual_removal, ratio_tree = tree_pred_5yr_removal / tree_obs_5yr_removal)]
metrics[variable == "removal (annual rate, dead + cut)", auc := rem$auc_site_rate[match(set, rem$set)]]
fwrite(metrics, file.path(VAL, "metrics.csv")); fwrite(rem, file.path(VAL, "removal_level_auc.csv"))
print(metrics, digits = 3); print(rem, digits = 3)
SP <- M[, {
  .(sites_present = sum(trees > 0, na.rm = TRUE), ba_r2 = metric(ba, ba_m)$r2, trees_r2 = metric(trees, trees_m)$r2,
    growth_r = metric(growth[growth_n > 0], growth_m[growth_n > 0])$r, growth_bias = metric(growth[growth_n > 0], growth_m[growth_n > 0])$bias,
    removal_auc = weighted_auc(mort_period_m[n_at_risk > 0], n_died[n_at_risk > 0], (n_at_risk - n_died)[n_at_risk > 0]),
    removal_obs_annual = sum(n_died) / sum(n_at_risk), removal_pred_annual = sum(n_at_risk * mort_m) / sum(n_at_risk), natural_obs_annual = sum(n_died_natural) / sum(n_at_risk),
    reg_r = metric(reg, reg_m)$r, reg_obs = mean(reg), reg_pred = mean(reg_m))
}, by = .(set, species, species_name)][order(set, species)]
fwrite(SP, file.path(VAL, "metrics_species.csv"))
TS <- TR[, .(trees = .N, obs_5yr = mean(died), pred_5yr = mean(mort_pred), growth_obs = mean(growth, na.rm = TRUE), growth_pred = mean(growth_pred[is.finite(growth)], na.rm = TRUE)), by = .(set, dcls)][order(set, dcls)]
fwrite(TS, file.path(VAL, "tree_by_size.csv"))
GS <- M[, .(gest = gest[1], obs_annual = sum(n_died) / sum(n_at_risk), pred_annual = sum(n_at_risk * mort_m) / sum(n_at_risk)), by = .(set, siteID)][, .(sites = .N, obs_annual = mean(obs_annual, na.rm = TRUE), pred_annual = mean(pred_annual, na.rm = TRUE)), by = .(set, gest)][order(set, gest)]
fwrite(GS, file.path(VAL, "removal_by_gest.csv"))
## holdout errors by the END-census out-of-S basal-area stratum (agreed rule 3): where does the model stop applying?
TOTs <- merge(TOT, unique(M[, .(set, siteID, end_stratum)]), by = c("set", "siteID"))
STR <- rbindlist(lapply(c("<1%", "1-5%", ">5%"), function(s) {
  t <- TOTs[set == "holdout" & end_stratum == s]; m <- M[set == "holdout" & end_stratum == s]
  data.table(end_stratum = s, sites = nrow(t), ba_r2 = metric(t$ba, t$ba_m)$r2, ba_r2_persist = metric(t$ba, t$ba_0)$r2, ba_rmse = metric(t$ba, t$ba_m)$rmse,
             trees_r2 = metric(t$trees, t$trees_m)$r2, trees_r2_persist = metric(t$trees, t$trees_0)$r2,
             growth_r2 = metric(m[growth_n > 0]$growth, m[growth_n > 0]$growth_m)$r2,
             removal_auc = weighted_auc(m[n_at_risk > 0]$mort_period_m, m[n_at_risk > 0]$n_died, m[n_at_risk > 0, n_at_risk - n_died]),
             reg_ratio = sum(t$reg_m) / sum(t$reg))
}))
fwrite(STR, file.path(VAL, "metrics_by_end_stratum.csv")); print(STR, digits = 3)

## ---- 3. figures ---------------------------------------------------------------------------------
loss_cols <- setdiff(names(loss_hist), c("epoch", "unused"))
LH <- melt(loss_hist[, c("epoch", loss_cols), with = FALSE][, (loss_cols) := lapply(.SD, as.numeric), .SDcols = loss_cols], id.vars = "epoch", variable.name = "response", value.name = "loss")[is.finite(loss)]
g1 <- ggplot(LH, aes(epoch, loss)) + geom_line() + facet_wrap(~ response, scales = "free_y") + theme_bw() + labs(title = "Training loss per response (baseline-scaled)", y = "loss")
ggsave(file.path(FIG, "fig1_loss_epochs.png"), g1, width = 9, height = 5, dpi = 150, bg = "white")
OP <- rbind(M[trees > 0, .(set, variable = "dbh (cm)", obs = dbh, pred = dbh_m)],
            TOT[, .(set, variable = "basal area, stand (m2/ha)", obs = ba / PARAMS$patch_size_ha, pred = ba_m / PARAMS$patch_size_ha)],
            TOT[, .(set, variable = "trees, stand (per ha)", obs = trees / PARAMS$patch_size_ha, pred = trees_m / PARAMS$patch_size_ha)],
            M[growth_n > 0, .(set, variable = "growth (rel. / yr)", obs = growth, pred = growth_m)],
            M[n_at_risk > 0, .(set, variable = "removal (annual rate)", obs = mort, pred = mort_m)],
            TOT[, .(set, variable = "recruitment, stand (per ha / 5 yr)", obs = reg, pred = reg_m)])
OP[, panel := factor(paste(set, variable, sep = ": "), levels = unique(paste(rep(levels(M$set), each = 6), unique(variable), sep = ": ")))]
g2 <- ggplot(OP[is.finite(obs) & is.finite(pred)], aes(obs, pred)) + geom_hex(bins = 60) + geom_abline(colour = "red") + facet_wrap(~ panel, scales = "free", nrow = 2) +
  scale_fill_viridis_c(trans = "log10") + theme_bw() + labs(title = "Observed vs predicted at year 5", x = "observed", y = "predicted")
ggsave(file.path(FIG, "fig2_obs_pred.png"), g2, width = 15, height = 6, dpi = 150, bg = "white")
ROC <- rbind(M[n_at_risk > 0, roc_points(mort_period_m, n_died, n_at_risk - n_died), by = set][, what := "simulated site rate"],
             TR[is.finite(mort_pred), roc_points(mort_pred, died, 1 - died), by = set][, what := "fitted per-tree probability"])
g3 <- ggplot(ROC, aes(fpr, tpr, colour = set, linetype = what)) + geom_path() + geom_abline(linetype = 3) + coord_equal() + theme_bw() +
  labs(title = "Removal ROC (5-yr probability vs tallied fates)", subtitle = paste(rem[, sprintf("%s: site AUC %.3f, per-tree AUC %.3f", set, auc_site_rate, auc_per_tree)], collapse = " | "), x = "false positive rate", y = "true positive rate")
ggsave(file.path(FIG, "fig3_removal_roc.png"), g3, width = 7, height = 6, dpi = 150, bg = "white")
MAP <- melt(TOT[, .(set, x, y, observed = ba / PARAMS$patch_size_ha, predicted = ba_m / PARAMS$patch_size_ha)], id.vars = c("set", "x", "y"), variable.name = "what", value.name = "ba_ha")
g4 <- ggplot(MAP, aes(x, y, colour = pmin(ba_ha, 60))) + geom_point(size = 0.3) + facet_grid(what ~ set) + coord_quickmap() + scale_colour_viridis_c(name = "BA m2/ha\n(capped 60)") + theme_bw() +
  labs(title = "Stand basal area at year 5: test and holdout blocks", x = "longitude", y = "latitude")
ggsave(file.path(FIG, "fig4_map.png"), g4, width = 9, height = 9, dpi = 150, bg = "white")
## bare ground at five climates x gest 0 / 2
TE <- unique(test$env[, .(siteID, tavg_c, prec_mm)]); qs <- quantile(TE$tavg_c, PARAMS$bare_quantiles)
loc <- rbindlist(lapply(seq_along(qs), function(i) data.table(clim = i, tavg_c = qs[[i]], prec_mm = median(TE[abs(tavg_c - qs[i]) <= 0.25, prec_mm]))))
loc <- loc[, .(gest = PARAMS$bare_gest), by = .(clim, tavg_c, prec_mm)][, siteID := .I]
env_bare <- loc[, .(year = seq_len(PARAMS$bare_years)), by = .(siteID, tavg_c, prec_mm, gest)][, .(siteID, year, tavg_c, prec_mm, gest)]
FINN.seed(PARAMS$seed)
B <- as.data.table(simulateForest(model, env = env_bare, init_cohort = NULL, patches = PARAMS$bare_patches, patch_size = PARAMS$patch_size_ha, device = "cpu")$wide$site)
for (cc in c("siteID", "year", "species")) B[, (cc) := as.integer(as.character(get(cc)))]
B <- merge(merge(B, species, by = "species"), loc[, .(siteID, label = sprintf("%.1f degC, %.0f mm", tavg_c, prec_mm), gest, tavg_c)], by = "siteID")
B[, label := factor(label, levels = unique(loc[order(tavg_c), sprintf("%.1f degC, %.0f mm", tavg_c, prec_mm)]))][, ba_ha := ba / PARAMS$patch_size_ha]
fwrite(B[, .(siteID, label, gest, year, species, species_name, ba_ha, trees_ha = trees / PARAMS$patch_size_ha)], file.path(VAL, "bare_ground.csv"))
top <- B[, .(m = max(ba_ha)), by = species_name][order(-m)][seq_len(PARAMS$n_top_species), species_name]
B[, shown := fifelse(species_name %in% top, species_name, "other species")]
g5 <- ggplot() + geom_line(data = B[shown == "other species"], aes(year, ba_ha, group = species_name), colour = "grey80", linewidth = 0.3) +
  geom_line(data = B[shown != "other species"], aes(year, ba_ha, colour = shown), linewidth = 0.6) + facet_grid(paste("gest", gest) ~ label) + theme_bw() +
  labs(title = sprintf("Succession from bare ground, %d patches of 0.0707 ha, constant climate, management trace 0 and 2", PARAMS$bare_patches), y = "basal area (m2/ha)", colour = "species")
ggsave(file.path(FIG, "fig5_bare_ground.png"), g5, width = 14, height = 7, dpi = 150, bg = "white")
EQ <- B[, .(ba_ha = sum(ba_ha), trees_ha = sum(trees) / PARAMS$patch_size_ha), by = .(label, gest, year)]
eq_tab <- EQ[, .(ba_last100 = mean(ba_ha[year > max(year) - 100]), ba_prev100 = mean(ba_ha[year > max(year) - 200 & year <= max(year) - 100]), trees_last100 = mean(trees_ha[year > max(year) - 100])), by = .(label, gest)][, change_pct := 100 * (ba_last100 / ba_prev100 - 1)]
PH <- model$param_history
if (length(PH) > 1) {
  CH <- rbindlist(lapply(2:length(PH), function(k) rbindlist(lapply(names(PH[[k]]), function(b) data.table(snapshot = k, block = b, mean_abs_change = mean(abs(as.matrix(PH[[k]][[b]]) - as.matrix(PH[[k - 1]][[b]]))))))))
  fwrite(CH, file.path(VAL, "parameter_change.csv"))
  g6 <- ggplot(CH[mean_abs_change > 0], aes(snapshot, mean_abs_change, colour = block)) + geom_line() + scale_y_log10() + theme_bw() + labs(title = "Mean absolute parameter change between snapshots", x = "snapshot (every 10 epochs)")
  ggsave(file.path(FIG, "fig6b_parameter_change.png"), g6, width = 9, height = 5, dpi = 150, bg = "white")
}

## ---- 4. report.md ---------------------------------------------------------------------------------
fit_summary <- if (file.exists(file.path(CAL_DIR, "fit_summary.md"))) readLines(file.path(CAL_DIR, "fit_summary.md")) else character(0)
sha <- if (file.exists(file.path(FIT_DIR, "parameters.pt"))) unname(tools::sha256sum(file.path(FIT_DIR, "parameters.pt"))) else NA
report <- c("# Validation report: fiveyear-v2-process (v0.2)", "",
  sprintf("Written %s by calibration/v02/validation.R; FINN %s (main 0848b60), torch %s. Model: %s/model.pt; parameters.pt sha256 %s.", format(Sys.time(), "%Y-%m-%d %H:%M %Z"), packageVersion("FINN"), packageVersion("torch"), FIT_DIR, sha), "",
  "## Sets", "",
  sprintf("- TEST = fiveyear-v2/train: %d sites the model was fitted on.", uniqueN(test$obs$siteID)),
  sprintf("- HOLDOUT = fiveyear-v2/holdout_blocked: %d sites in 100 km blocks never seen. No coordinates in the model.", uniqueN(hold$obs$siteID)),
  "- Aggregation as fit() scores it: state at year 5; growth and removal = mean annual rate over years 1..5; recruitment = sum. Removal = dead + cut (total removal); natural-only for reference.", "",
  "## Metrics", "", md_table(metrics[, .(set, variable, level, n, rmse, r2, r, bias, obs_mean, auc)]), "",
  "## Removal level and AUC on one scale", "", md_table(rem), "",
  "`share_trees_below_floor`: share of tallied trees whose fitted annual removal probability is below the 1 %/yr floor FINN applies when it draws the dying trees in simulation; for them the simulated stand loses more trees than the fitted rate says.", "",
  "## Per-tree removal and growth by start diameter", "", md_table(TS), "",
  "## Removal by management trace (site means)", "", md_table(GS), "",
  "## Holdout by end-census out-of-S basal-area stratum (the model's domain: species outside S < 2 % at the START census; the end-census share says where it stops applying)", "", md_table(STR), "",
  "## Per species, holdout", "", md_table(SP[set == "holdout"][order(-sites_present), .(species_name, sites_present, ba_r2, trees_r2, growth_r, growth_bias, removal_auc, removal_obs_annual, removal_pred_annual, natural_obs_annual, reg_r, reg_obs, reg_pred)]), "",
  "## Bare-ground equilibrium (stand totals, last 100 vs previous 100 years)", "", md_table(eq_tab), "",
  "## Figures", "", "- figures/fig1_loss_epochs.png", "- figures/fig2_obs_pred.png", "- figures/fig3_removal_roc.png", "- figures/fig4_map.png", "- figures/fig5_bare_ground.png", "- figures/fig6b_parameter_change.png", "",
  "## Provenance", "", "| step | script | inputs | outputs |", "|---|---|---|---|",
  "| raw data | calibration/01_download_ign_raw.R | IGN DataIFN export 2005-2025 | raw/ign_data/unz, PROVENANCE.md |",
  "| preparation v2 | calibration/v02/02_prepare_fiveyear_v2.R | raw tables, WorldClim 2.1 | prepared/fiveyear-v2/*, PREPARATION_v2.md, expected_preparation_v2.csv |",
  "| audit | calibration/v02/audit_v01.R, pilots pilot_*.R, pilot_eval.R | v0.1 fit, fiveyear-v2 pilot sets | AUDIT_v02.md, pilot_table.csv |",
  "| fit | calibration/v02/fit.R | fiveyear-v2/train | fits/fiveyear-v2-process/{model.pt, parameters.pt}, loss_history.csv, fit_summary.md |",
  "| validation | calibration/v02/validation.R | model.pt, train, holdout_blocked | report.md, validation/*.csv, figures/*.png |", "",
  "## Fit summary (from fit.R)", "", fit_summary)
writeLines(report, file.path(CAL_DIR, "report.md"))
if (dir.exists(ZOO_DIR)) {
  dir.create(file.path(ZOO_DIR, "figures"), showWarnings = FALSE)
  fwrite(metrics, file.path(ZOO_DIR, "metrics.csv"))
  file.copy(file.path(DATA_DIR, "prepared", PARAMS$prep_version, "species.csv"), file.path(ZOO_DIR, "species.csv"), overwrite = TRUE)
  for (f in list.files(FIG, pattern = "\\.png$", full.names = TRUE)) file.copy(f, file.path(ZOO_DIR, "figures", basename(f)), overwrite = TRUE)
}
message("done: ", file.path(CAL_DIR, "report.md"))
