#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-IFN-France | calibration/v02/pilot_eval.R
##
## WHAT   Scores every v0.2 pilot in <DATA_DIR>/fits/v02-pilots/<name>/model.pt on
##        the SAME holdout subset (fiveyear-v2/holdout1000) and on its training
##        set (pilot2000): simulate the 5-year interval from the first-visit
##        stand (state at year 5, mean annual rates, summed recruitment), and
##        the fitted PER-TREE removal probability from FINN::predictTrees (not
##        subject to the simulation's 1 %/yr floor). Against persistence.
##        One row per pilot x set, plus the holdout rows STRATIFIED by the
##        end-census out-of-S basal-area share (< 1 %, 1-5 %, > 5 %; agreed
##        rule 3): if the error grows with the stratum, the card says where the
##        model stops applying.
##        Columns: BA / stems R2 of the stand totals and the persistence R2,
##        correlation of the 5-yr changes (dba_r, dn_r) and their bias,
##        species-level BA R2, growth R2 and level ratio, removal: per-tree AUC
##        and 5-yr level ratio (fitted / observed) and the SIMULATED annual rate
##        ratio (floored stand), recruitment ratio and r, bare-ground stand
##        basal area after 300 years at two climates (gest 0), dominant species
##        and drift, median theta, loss trend over the last third.
##        A pilot whose model has |S| + 1 species (the shade arm) is simulated
##        from the shade initial cohorts; its scores refer to S only.
## RUN    Rscript calibration/v02/pilot_eval.R   (one thread, ~3 min per pilot)
## OUTPUT calibration/v02/pilot_table.csv, pilot_table.md, pilot_table_strata.csv;
##        <DATA_DIR>/fits/v02-pilots/<name>/eval_predictions.csv
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-zoo-v02"), path.expand("~/R/lib-finn-0848b60"), .libPaths()))   # the branch build reads every pilot (superset of 0.2.0)
suppressPackageStartupMessages({ library(FINN); library(data.table); library(torch) })
torch_set_num_threads(1L)
DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "frnfi")
PIL_DIR  <- file.path(DATA_DIR, "fits", "v02-pilots")
PARAMS <- list(prep_version = "fiveyear-v2", sets = c(train = "pilot2000", holdout = "holdout1000"),
               patch_size_ha = pi * 15^2 / 1e4, period = 5L, env_vars = c("tavg_c", "prec_mm", "gest"),
               bare_years = 300L, bare_patches = 10L, bare_quantiles = c(0.5, 0.9), seed = 1L)
P <- PARAMS$period
r2 <- function(o, p) { ok <- is.finite(o) & is.finite(p); o <- o[ok]; p <- p[ok]; if (length(o) < 3) return(NA_real_); 1 - sum((o - p)^2) / sum((o - mean(o))^2) }
weighted_auc <- function(p, n_pos, n_neg) {
  dt <- data.table(p = p, n_pos = n_pos, n_neg = n_neg)[is.finite(p) & (n_pos + n_neg) > 0][, .(n_pos = sum(n_pos), n_neg = sum(n_neg)), by = p][order(p)]
  if (!nrow(dt) || sum(dt$n_pos) == 0 || sum(dt$n_neg) == 0) return(NA_real_)
  before <- cumsum(dt$n_neg) - dt$n_neg; sum(dt$n_pos * (before + 0.5 * dt$n_neg)) / (sum(dt$n_pos) * sum(dt$n_neg))
}
read_set <- function(set) { d <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, set)
  list(obs = fread(file.path(d, "obs_dt.csv")), env = fread(file.path(d, "env_dt.csv")), init = fread(file.path(d, "init_cohort_dt.csv")),
       init_shade = fread(file.path(d, "init_cohort_shade_dt.csv")), trees = fread(file.path(d, "tree_obs_dt.csv")),
       species = fread(file.path(d, "species_dt.csv")), sites = fread(file.path(d, "siteID_dt.csv"))) }
sets <- lapply(PARAMS$sets, read_set); n_S <- nrow(sets$train$species)
predict_set <- function(model, d, n_model) {
  env <- d$env[, c("siteID", "year", PARAMS$env_vars), with = FALSE]
  init <- if (n_model > n_S) d$init_shade else d$init
  FINN.seed(PARAMS$seed)
  S <- as.data.table(simulateForest(model, env = env, init_cohort = makeInitCohorts(init[, .(siteID, patchID, species, dbh, trees, treeName)], Nspecies = n_model),
                                    patches = 1L, patch_size = PARAMS$patch_size_ha, device = "cpu", batchsize = 500L)$wide$site)
  for (cc in c("siteID", "year", "species")) S[, (cc) := as.integer(as.character(get(cc)))]
  S <- S[species <= n_S]
  A <- S[year <= P, .(growth_m = mean(growth), mort_m = mean(mort), reg_m = sum(r_mean_ha)), by = .(siteID, species)]
  M <- merge(d$obs, merge(S[year == P, .(siteID, species, ba_m = ba, trees_m = trees)], A, by = c("siteID", "species")), by = c("siteID", "species"))
  I <- d$init[, .(ba_0 = sum(dbh2ba(dbh) * trees), trees_0 = sum(trees)), by = .(siteID, species)]
  M <- merge(M, I, by = c("siteID", "species"), all.x = TRUE); M[is.na(ba_0), `:=`(ba_0 = 0, trees_0 = 0)]
  merge(M, d$sites[, .(siteID, end_stratum, end_share_out)], by = "siteID")
}
tree_pred <- function(model, d) {
  td <- d$trees[, .(siteID, year, patchID, species, dbh, trees, growth, died, period_length)]
  pt <- predictTrees(model, tree_data = td, env = d$env[, c("siteID", "year", PARAMS$env_vars), with = FALSE], patches = 1L, patch_size = PARAMS$patch_size_ha, device = "cpu")
  stopifnot(nrow(pt) == nrow(td))
  merge(cbind(td, pt[, .(growth_pred, mort_pred)]), d$sites[, .(siteID, end_stratum)], by = "siteID")
}
score <- function(M, T2) {
  T <- M[, .(ba = sum(ba), ba_m = sum(ba_m), ba_0 = sum(ba_0), n = sum(trees), n_m = sum(trees_m), n_0 = sum(trees_0), reg = sum(reg), reg_m = sum(reg_m)), by = siteID]
  d <- M[n_at_risk > 0]; g <- M[growth_n > 0]; t <- T2[is.finite(mort_pred)]
  data.table(sites = nrow(T), ba_r2 = r2(T$ba, T$ba_m), ba_r2_persist = r2(T$ba, T$ba_0), dba_r = cor(T$ba - T$ba_0, T$ba_m - T$ba_0), dba_bias = mean(T$ba_m - T$ba),
             n_r2 = r2(T$n, T$n_m), n_r2_persist = r2(T$n, T$n_0), dn_r = cor(T$n - T$n_0, T$n_m - T$n_0), dn_bias = mean(T$n_m - T$n),
             sp_ba_r2 = r2(M$ba, M$ba_m), growth_r2 = r2(g$growth, g$growth_m), growth_ratio = weighted.mean(g$growth_m, g$growth_n) / weighted.mean(g$growth, g$growth_n),
             tree_growth_r2 = r2(T2$growth, T2$growth_pred),
             removal_auc_tree = weighted_auc(t$mort_pred, t$died, 1 - t$died), removal_ratio_tree = mean(t$mort_pred) / mean(t$died),
             removal_obs_5yr = mean(t$died), removal_ratio_sim_annual = sum(d$n_at_risk * d$mort_m) / sum(d$n_died),
             reg_ratio = sum(T$reg_m) / sum(T$reg), reg_r = cor(T$reg, T$reg_m), reg_sp_r = cor(M$reg, M$reg_m))
}
bare <- function(model, n_model) {
  TE <- unique(sets$train$env[, .(siteID, tavg_c, prec_mm)]); qs <- quantile(TE$tavg_c, PARAMS$bare_quantiles)
  loc <- rbindlist(lapply(seq_along(qs), function(i) data.table(siteID = i, tavg_c = qs[[i]], prec_mm = median(TE[abs(tavg_c - qs[i]) < 0.25, prec_mm]), gest = 0)))
  env <- loc[, .(year = seq_len(PARAMS$bare_years)), by = .(siteID, tavg_c, prec_mm, gest)][, .(siteID, year, tavg_c, prec_mm, gest)]
  FINN.seed(PARAMS$seed)
  B <- as.data.table(simulateForest(model, env = env, init_cohort = NULL, patches = PARAMS$bare_patches, patch_size = PARAMS$patch_size_ha, device = "cpu")$wide$site)
  for (cc in c("siteID", "year", "species")) B[, (cc) := as.integer(as.character(get(cc)))]
  B <- merge(B, rbind(sets$train$species, data.table(species = n_S + 1L, species_name = "shade class")), by = "species")
  tot <- B[, .(ba_ha = sum(ba) / PARAMS$patch_size_ha), by = .(siteID, year)]
  dom <- B[year > PARAMS$bare_years - 50, .(ba = sum(ba)), by = .(siteID, species_name)][order(siteID, -ba)][, .SD[1], by = siteID]
  data.table(bare_ba_median_clim = tot[siteID == 1 & year > PARAMS$bare_years - 50, mean(ba_ha)], bare_ba_warm_clim = tot[siteID == 2 & year > PARAMS$bare_years - 50, mean(ba_ha)],
             bare_dominant_median = dom[siteID == 1, species_name], bare_dominant_warm = dom[siteID == 2, species_name],
             bare_drift_pct = 100 * (tot[siteID == 1 & year > PARAMS$bare_years - 50, mean(ba_ha)] / tot[siteID == 1 & year > PARAMS$bare_years - 100 & year <= PARAMS$bare_years - 50, mean(ba_ha)] - 1),
             bare_shade_share = if (n_model > n_S) B[year > PARAMS$bare_years - 50, sum(ba[species == n_S + 1L]) / sum(ba)] else NA_real_)
}
pilots <- list.dirs(PIL_DIR, recursive = FALSE, full.names = FALSE); pilots <- pilots[file.exists(file.path(PIL_DIR, pilots, "model.pt"))]
rows <- list(); strata <- list()
for (nm in pilots) {
  message("== ", nm)
  model <- torch_load(file.path(PIL_DIR, nm, "model.pt")); model$eval(); n_model <- model$N_species
  Mh <- predict_set(model, sets$holdout, n_model); Mt <- predict_set(model, sets$train, n_model)
  Th <- tree_pred(model, sets$holdout); Tt <- tree_pred(model, sets$train)
  fwrite(rbind(cbind(set = "holdout", Mh), cbind(set = "train", Mt)), file.path(PIL_DIR, nm, "eval_predictions.csv"))
  lh <- fread(file.path(PIL_DIR, nm, "loss_history.csv")); n <- nrow(lh); last <- lh[epoch > 2 * n / 3]
  tot_loss <- rowSums(last[, .(dbh, ba, trees, growth, mortality, regeneration)]); trend <- unname(coef(lm(tot_loss ~ seq_along(tot_loss)))[2]) * length(tot_loss) / mean(tot_loss)
  th <- as.numeric(model$par_theta_recruits$detach())
  rows[[nm]] <- cbind(pilot = nm, rbind(cbind(set = "holdout", score(Mh, Th)), cbind(set = "train", score(Mt, Tt))), theta_median = median(th), loss_trend_last_third = trend, bare(model, n_model)[rep(1, 2)])
  strata[[nm]] <- rbindlist(lapply(levels(factor(Mh$end_stratum, levels = c("<1%", "1-5%", ">5%"))), function(s) cbind(pilot = nm, set = "holdout", end_stratum = s, score(Mh[end_stratum == s], Th[end_stratum == s]))))
}
rows <- rbindlist(rows, fill = TRUE); strata <- rbindlist(strata, fill = TRUE)
fwrite(rows, file.path("calibration", "v02", "pilot_table.csv")); fwrite(strata, file.path("calibration", "v02", "pilot_table_strata.csv"))
fmt <- function(dt) { num <- names(dt)[vapply(dt, is.numeric, logical(1))]; s <- copy(dt)[, (num) := lapply(.SD, signif, 3), .SDcols = num]
  c(paste0("| ", paste(names(s), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(s)), collapse = "|"), "|"), apply(s, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))) }
md <- c("# v0.2 pilots: holdout1000 (and the pilot2000 training set) against persistence", "",
        sprintf("Written %s by calibration/v02/pilot_eval.R. *_r2 = R2 of the stand totals at year 5 (persistence = first-visit stand as the prediction); dba_r / dn_r = correlation of predicted with observed 5-year change; growth_r2 on species x site rows, tree_growth_r2 per tallied tree; removal_auc_tree / removal_ratio_tree = fitted per-tree 5-yr removal probability (predictTrees, not floored) against the tallied fates; removal_ratio_sim_annual = simulated annual rate / site-level target (both inflated by the site-level annualisation, see AUDIT); reg_ratio = predicted / observed recruits; bare_* = stand basal area (m2/ha) in years 250-300 from bare ground at the median and the 90 %% quantile of the training tavg_c (gest 0).", format(Sys.time(), "%Y-%m-%d %H:%M")), "",
        fmt(rows), "", "## Holdout by end-census out-of-S basal-area stratum", "", fmt(strata))
writeLines(md, file.path("calibration", "v02", "pilot_table.md"))
print(rows[, .(pilot, set, ba_r2, ba_r2_persist, n_r2, n_r2_persist, dba_r, dn_r, growth_r2, removal_auc_tree, removal_ratio_tree, reg_ratio, reg_r, bare_ba_median_clim, bare_dominant_median)], digits = 3)
print(strata[, .(pilot, end_stratum, sites, ba_r2, ba_r2_persist, n_r2, growth_r2, removal_auc_tree, reg_ratio)], digits = 3)
