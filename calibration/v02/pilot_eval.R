#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-IFN-France | calibration/v02/pilot_eval.R
##
## WHAT   Scores every v0.2 pilot in <DATA_DIR>/fits/v02-pilots/<name>/model.pt on
##        the SAME holdout subset (fiveyear-v2/holdout1000) and on its training
##        set (pilot2000), the way the v0.1 validation does (simulate the 5-year
##        interval from the first-visit stand; state at year 5; mean annual
##        rates; summed recruitment), and against persistence. One table:
##          BA / stems R2 (stand totals) and the persistence R2, change metrics
##          (r of dBA and dN), growth R2, removal AUC (weighted by trees at
##          risk) and the annual removal ratio pred/obs, recruitment ratio and
##          r, bare-ground stand basal area after 300 years at two climates
##          (median training tavg_c, and the 90 % quantile) with gest 0,
##          and the loss trend over the last third of the epochs.
## RUN    Rscript calibration/v02/pilot_eval.R   (one thread, ~2 min per pilot)
## OUTPUT calibration/v02/pilot_table.csv, pilot_table.md; per pilot
##        <DATA_DIR>/fits/v02-pilots/<name>/eval_predictions.csv
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-0848b60"), .libPaths()))
suppressPackageStartupMessages({ library(FINN); library(data.table); library(torch) })
torch_set_num_threads(1L)
DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "frnfi")
PIL_DIR  <- file.path(DATA_DIR, "fits", "v02-pilots")
PARAMS <- list(prep_version = "fiveyear-v2", sets = c(train = "pilot2000", holdout = "holdout1000"),
               patch_size_ha = pi * 15^2 / 1e4, period = 5L, env_vars = c("tavg_c", "prec_mm", "gest"),
               bare_years = 300L, bare_patches = 10L, bare_quantiles = c(0.5, 0.9), seed = 1L)
P <- PARAMS$period
r2 <- function(o, p) { ok <- is.finite(o) & is.finite(p); o <- o[ok]; p <- p[ok]; 1 - sum((o - p)^2) / sum((o - mean(o))^2) }
weighted_auc <- function(p, n_pos, n_neg) {
  dt <- data.table(p = p, n_pos = n_pos, n_neg = n_neg)[is.finite(p) & (n_pos + n_neg) > 0][, .(n_pos = sum(n_pos), n_neg = sum(n_neg)), by = p][order(p)]
  if (!nrow(dt) || sum(dt$n_pos) == 0 || sum(dt$n_neg) == 0) return(NA_real_)
  before <- cumsum(dt$n_neg) - dt$n_neg; sum(dt$n_pos * (before + 0.5 * dt$n_neg)) / (sum(dt$n_pos) * sum(dt$n_neg))
}
read_set <- function(set) { d <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, set)
  list(obs = fread(file.path(d, "obs_dt.csv")), env = fread(file.path(d, "env_dt.csv")), init = fread(file.path(d, "init_cohort_dt.csv")), species = fread(file.path(d, "species_dt.csv"))) }
sets <- lapply(PARAMS$sets, read_set); n_species <- nrow(sets$train$species)
## the env formula may contain I(x^2): simulateForest only needs the raw columns
predict_set <- function(model, d) {
  env <- d$env[, c("siteID", "year", PARAMS$env_vars), with = FALSE]
  FINN.seed(PARAMS$seed)
  S <- as.data.table(simulateForest(model, env = env, init_cohort = makeInitCohorts(d$init[, .(siteID, patchID, species, dbh, trees, treeName)], Nspecies = n_species),
                                    patches = 1L, patch_size = PARAMS$patch_size_ha, device = "cpu", batchsize = 500L)$wide$site)
  for (cc in c("siteID", "year", "species")) S[, (cc) := as.integer(as.character(get(cc)))]
  A <- S[year <= P, .(growth_m = mean(growth), mort_m = mean(mort), reg_m = sum(r_mean_ha)), by = .(siteID, species)]
  M <- merge(d$obs, merge(S[year == P, .(siteID, species, ba_m = ba, trees_m = trees)], A, by = c("siteID", "species")), by = c("siteID", "species"))
  I <- d$init[, .(ba_0 = sum(dbh2ba(dbh) * trees), trees_0 = sum(trees)), by = .(siteID, species)]
  M <- merge(M, I, by = c("siteID", "species"), all.x = TRUE); M[is.na(ba_0), `:=`(ba_0 = 0, trees_0 = 0)]
  M
}
score <- function(M) {
  T <- M[, .(ba = sum(ba), ba_m = sum(ba_m), ba_0 = sum(ba_0), n = sum(trees), n_m = sum(trees_m), n_0 = sum(trees_0), reg = sum(reg), reg_m = sum(reg_m)), by = siteID]
  d <- M[n_at_risk > 0]; g <- M[growth_n > 0]
  data.table(ba_r2 = r2(T$ba, T$ba_m), ba_r2_persist = r2(T$ba, T$ba_0), dba_r = cor(T$ba - T$ba_0, T$ba_m - T$ba_0), dba_bias = mean(T$ba_m - T$ba),
             n_r2 = r2(T$n, T$n_m), n_r2_persist = r2(T$n, T$n_0), dn_r = cor(T$n - T$n_0, T$n_m - T$n_0), dn_bias = mean(T$n_m - T$n),
             sp_ba_r2 = r2(M$ba, M$ba_m), growth_r2 = r2(g$growth, g$growth_m), growth_ratio = weighted.mean(g$growth_m, g$growth_n) / weighted.mean(g$growth, g$growth_n),
             removal_auc = weighted_auc(d$mort_m, d$n_died, d$n_at_risk - d$n_died), removal_ratio = sum(d$n_at_risk * d$mort_m) / sum(d$n_died),
             reg_ratio = sum(T$reg_m) / sum(T$reg), reg_r = cor(T$reg, T$reg_m), reg_sp_r = cor(M$reg, M$reg_m))
}
bare <- function(model) {
  TE <- unique(sets$train$env[, .(siteID, tavg_c, prec_mm)]); qs <- quantile(TE$tavg_c, PARAMS$bare_quantiles)
  loc <- rbindlist(lapply(seq_along(qs), function(i) data.table(siteID = i, tavg_c = qs[[i]], prec_mm = median(TE[abs(tavg_c - qs[i]) < 0.25, prec_mm]), gest = 0)))
  env <- loc[, .(year = seq_len(PARAMS$bare_years)), by = .(siteID, tavg_c, prec_mm, gest)][, .(siteID, year, tavg_c, prec_mm, gest)]
  FINN.seed(PARAMS$seed)
  B <- as.data.table(simulateForest(model, env = env, init_cohort = NULL, patches = PARAMS$bare_patches, patch_size = PARAMS$patch_size_ha, device = "cpu")$wide$site)
  for (cc in c("siteID", "year", "species")) B[, (cc) := as.integer(as.character(get(cc)))]
  B <- merge(B, sets$train$species, by = "species")
  tot <- B[, .(ba_ha = sum(ba) / PARAMS$patch_size_ha), by = .(siteID, year)]
  dom <- B[year > PARAMS$bare_years - 50, .(ba = sum(ba)), by = .(siteID, species_name)][order(siteID, -ba)][, .SD[1], by = siteID]
  data.table(bare_ba_median_clim = tot[siteID == 1 & year > PARAMS$bare_years - 50, mean(ba_ha)], bare_ba_warm_clim = tot[siteID == 2 & year > PARAMS$bare_years - 50, mean(ba_ha)],
             bare_dominant_median = dom[siteID == 1, species_name], bare_dominant_warm = dom[siteID == 2, species_name],
             bare_drift_pct = 100 * (tot[siteID == 1 & year > PARAMS$bare_years - 50, mean(ba_ha)] / tot[siteID == 1 & year > PARAMS$bare_years - 100 & year <= PARAMS$bare_years - 50, mean(ba_ha)] - 1))
}
pilots <- list.dirs(PIL_DIR, recursive = FALSE, full.names = FALSE)
pilots <- pilots[file.exists(file.path(PIL_DIR, pilots, "model.pt"))]
rows <- rbindlist(lapply(pilots, function(nm) {
  message("== ", nm)
  model <- torch_load(file.path(PIL_DIR, nm, "model.pt")); model$eval()
  Mh <- predict_set(model, sets$holdout); Mt <- predict_set(model, sets$train)
  fwrite(rbind(cbind(set = "holdout", Mh), cbind(set = "train", Mt)), file.path(PIL_DIR, nm, "eval_predictions.csv"))
  lh <- fread(file.path(PIL_DIR, nm, "loss_history.csv")); n <- nrow(lh); last <- lh[epoch > 2 * n / 3]
  tot_loss <- rowSums(last[, .(dbh, ba, trees, growth, mortality, regeneration)])
  trend <- unname(coef(lm(tot_loss ~ seq_along(tot_loss)))[2]) * length(tot_loss) / mean(tot_loss)
  th <- as.numeric(model$par_theta_recruits$detach())
  cbind(pilot = nm, rbind(cbind(set = "holdout", score(Mh)), cbind(set = "train", score(Mt))),
        theta_median = median(th), loss_trend_last_third = trend, bare(model)[rep(1, 2)])
}), fill = TRUE)
fwrite(rows, file.path("calibration", "v02", "pilot_table.csv"))
num <- names(rows)[vapply(rows, is.numeric, logical(1))]; show <- copy(rows)[, (num) := lapply(.SD, signif, 3), .SDcols = num]
md <- c("# v0.2 pilots: holdout1000 (and the pilot2000 training set) against persistence", "",
        sprintf("Written %s by calibration/v02/pilot_eval.R. Columns: *_r2 = R2 of the stand totals at year 5 (persistence = first-visit stand as the prediction); dba_r / dn_r = correlation of the predicted with the observed 5-year change; growth_r2 on species x site rows; removal_auc weighted by trees at risk; removal_ratio = predicted / observed annual removal (total removal, incl. cuts); reg_ratio = predicted / observed recruits; bare_* = stand basal area (m2/ha) in years 250-300 from bare ground at the median and the 90 %% quantile of the training tavg_c (gest 0), with the dominant species and the drift between the last two 50-year windows.", format(Sys.time(), "%Y-%m-%d %H:%M")), "",
        paste0("| ", paste(names(show), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(show)), collapse = "|"), "|"),
        apply(show, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
writeLines(md, file.path("calibration", "v02", "pilot_table.md"))
print(rows[, .(pilot, set, ba_r2, ba_r2_persist, n_r2, n_r2_persist, dba_r, dn_r, growth_r2, removal_auc, removal_ratio, reg_ratio, reg_r, bare_ba_median_clim, bare_dominant_median)], digits = 3)
