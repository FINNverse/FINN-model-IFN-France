#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-IFN-France | calibration/v02/audit_v01.R
##
## WHAT   The evidence behind calibration/v02/AUDIT_v02.md: reads the v0.1 fit
##        (fits/fiveyear-v1-process: model.pt, loss_history.csv, validation
##        predictions) and the v1 data and writes tables and small figures to
##        calibration/v02/audit/. Nothing is fitted or simulated anew except
##        reading parameter values out of the saved model.
##        Sections: (1) mortality level on ONE scale (the v0.1 validation mixed
##        annual observed with 5-year predicted shares); by size, species and
##        management trace; (2) recruitment: predicted vs observed per species,
##        the recruitment dispersion theta and the cap K per species over the
##        fit; (3) change metrics (dBA, dN) vs persistence; (4) growth by size;
##        (5) loss and parameter trajectories; (6) data: fates by management
##        trace, new-tree diameters, back-modelled start diameters.
## INPUT  <DATA_DIR>/fits/fiveyear-v1-process/{model.pt, loss_history.csv, validation/predictions.csv}
##        <DATA_DIR>/prepared/fiveyear-v1/{site_table.csv, national/obs_dt.csv}
##        <DATA_DIR>/raw/ign_data/unz/{ARBRE.csv, PLACETTE.csv}
## OUTPUT calibration/v02/audit/*.csv, *.png
## RUN    Rscript calibration/v02/audit_v01.R   (from the repository root, ~3 min, one thread)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-0848b60"), .libPaths()))
suppressPackageStartupMessages({ library(FINN); library(data.table); library(torch); library(ggplot2) })
torch_set_num_threads(1L)
DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "frnfi")
FIT_DIR  <- file.path(DATA_DIR, "fits", "fiveyear-v1-process")
RAW      <- file.path(DATA_DIR, "raw", "ign_data", "unz")
OUT      <- file.path("calibration", "v02", "audit"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
PARAMS <- list(period = 5L, patch_size_ha = pi * 15^2 / 1e4, dbh_breaks = c(15, 25, 40, 60),
               fate_dead = c("M", "5", "C", "2", "A", "T", "1"), fate_cut = c("6", "7"))
P <- PARAMS$period
P1 <- fread(file.path(FIT_DIR, "validation", "predictions.csv"))
site <- fread(file.path(DATA_DIR, "prepared", "fiveyear-v1", "site_table.csv"))
P1 <- merge(P1, site[, .(siteName = as.character(IDP), GEST)], by = "siteName", all.x = TRUE)
H <- P1[set == "holdout"]
H[, dcls := cut(dbh_0, c(0, PARAMS$dbh_breaks, Inf))]

## ---- 1. mortality level on one scale ---------------------------------------------------
## the v0.1 report compared observed ANNUAL deaths (n_died = annual rate x n_at_risk) with
## the predicted 5-YEAR death probability; here both are annual (pred = mean annual rate)
lvl <- P1[n_at_risk > 0, .(n_at_risk = sum(n_at_risk), obs_annual = sum(n_died) / sum(n_at_risk),
                            pred_annual = sum(n_at_risk * mort_m) / sum(n_at_risk)), by = set][, ratio := pred_annual / obs_annual]
fwrite(lvl, file.path(OUT, "mortality_level.csv")); print(lvl, digits = 3)
by_size <- H[n_at_risk > 0, .(n_at_risk = sum(n_at_risk), obs = sum(n_died) / sum(n_at_risk), pred = sum(n_at_risk * mort_m) / sum(n_at_risk)), by = dcls][order(dcls)][, ratio := pred / obs]
by_sp <- H[n_at_risk > 0, .(n_at_risk = sum(n_at_risk), obs = sum(n_died) / sum(n_at_risk), pred = sum(n_at_risk * mort_m) / sum(n_at_risk)), by = species_name][order(-n_at_risk)][, ratio := pred / obs]
by_gest <- H[n_at_risk > 0, .(n_at_risk = sum(n_at_risk), obs = sum(n_died) / sum(n_at_risk), pred = sum(n_at_risk * mort_m) / sum(n_at_risk)), by = GEST][order(GEST)][, ratio := pred / obs]
fwrite(by_size, file.path(OUT, "mortality_by_size.csv")); fwrite(by_sp, file.path(OUT, "mortality_by_species.csv")); fwrite(by_gest, file.path(OUT, "mortality_by_gest.csv"))
print(by_size, digits = 3); print(by_gest, digits = 3)
g1 <- ggplot(melt(by_size, id.vars = c("dcls", "n_at_risk", "ratio")), aes(dcls, value, fill = variable)) + geom_col(position = "dodge") + theme_bw() +
  labs(title = "v0.1 holdout: annual natural mortality by start diameter", x = "dbh at the first visit (cm)", y = "annual rate", fill = "")
ggsave(file.path(OUT, "fig_mortality_by_size.png"), g1, width = 7, height = 4, dpi = 130, bg = "white")

## ---- 1b. the simulation's 1 %/yr mortality floor vs the fitted rate -------------------------
## FINN draws the dying trees with the annual rate clamped to [0.01, 0.99]
## (finn-class.R, binomial_from_gamma(..., torch_clamp(m, 0.01, 0.99))), so the
## SIMULATED stand never loses fewer than 1 %/yr even where the fitted rate is
## lower; the rate itself (mort in the output, the site-level loss, predictTrees)
## is not clamped. Separate the two: the fitted per-tree 5-year death probability
## (FINN::predictTrees on the v2 tallied-tree table of the holdout1000 plots, v1
## species index) against the observed NATURAL deaths (the v0.1 target), by
## species; and how much of the stand sits below the floor.
v2d <- file.path(DATA_DIR, "prepared", "fiveyear-v2", "holdout1000")
if (file.exists(file.path(v2d, "tree_obs_dt.csv"))) {
  model <- torch_load(file.path(FIT_DIR, "model.pt")); model$eval()
  tr <- fread(file.path(v2d, "tree_obs_dt.csv")); ev <- fread(file.path(v2d, "env_dt.csv"))
  sp1 <- fread(file.path(DATA_DIR, "prepared", "fiveyear-v1", "species.csv"))
  tr[, species_v1 := sp1$species[match(fifelse(species_name %in% sp1$species_name, species_name, "other"), sp1$species_name)]]
  td <- tr[, .(siteID, year, patchID, species = species_v1, dbh, trees, growth, died = died_natural, period_length)]
  pt <- predictTrees(model, tree_data = td, env = ev[, .(siteID, year, tavg_c, prec_mm)], patches = 1L, patch_size = PARAMS$patch_size_ha, device = "cpu")
  td[, `:=`(p5 = pt$mort_pred, g_pred = pt$growth_pred)]
  td[, annual_pred := 1 - (1 - p5)^(1 / P)]
  floor_tab <- td[, .(trees = .N, obs_5yr_natural = mean(died), pred_5yr = mean(p5), pred_annual = mean(annual_pred),
                      share_trees_below_1pct = mean(annual_pred < 0.01)), by = .(species_name = sp1$species_name[species])][order(-trees)]
  fwrite(floor_tab, file.path(OUT, "mortality_floor_per_tree.csv"))
  cat(sprintf("per tree (holdout1000, v0.1 model): observed 5-yr natural death share %.4f, fitted %.4f (ratio %.2f); %.0f %% of the trees have a fitted annual rate below the 1 %% simulation floor\n",
              mean(td$died), mean(td$p5), mean(td$p5) / mean(td$died), 100 * mean(td$annual_pred < 0.01)))
  print(floor_tab[1:12], digits = 3)
  g1b <- ggplot(floor_tab[trees > 200], aes(obs_5yr_natural, pred_5yr, label = species_name)) + geom_abline(colour = "grey60") +
    geom_point(aes(size = trees)) + geom_text(size = 2.5, vjust = -0.8) + geom_vline(xintercept = 1 - 0.99^5, linetype = 2, colour = "red") + theme_bw() +
    labs(title = "v0.1: fitted per-tree 5-yr death probability vs observed natural death share (holdout1000)",
         subtitle = "red: 1 %/yr simulation floor over 5 yr; species left of it lose more trees in simulation than the fitted rate says", x = "observed", y = "fitted (not floored)")
  ggsave(file.path(OUT, "fig_mortality_floor.png"), g1b, width = 7, height = 5, dpi = 130, bg = "white")
}

## ---- 2. recruitment: species, theta, K ------------------------------------------------------
rec_sp <- H[, .(obs = mean(reg), pred = mean(reg_m), r = if (sd(reg) > 0 && sd(reg_m) > 0) cor(reg, reg_m) else NA_real_), by = species_name][order(-pred)][, ratio := pred / obs]
fwrite(rec_sp, file.path(OUT, "recruitment_by_species.csv")); print(rec_sp[1:10], digits = 3)
model <- torch_load(file.path(FIT_DIR, "model.pt"))
sp <- fread(file.path(DATA_DIR, "prepared", "fiveyear-v1", "species.csv"))[order(species)]
par_tab <- data.table(species = sp$species, species_name = sp$species_name,
                      theta = as.numeric(model$par_theta_recruits$detach()), K_ha_yr = exp(as.numeric(model$reg_logK$detach())))
par_tab <- merge(par_tab, rec_sp[, .(species_name, reg_obs = obs, reg_pred = pred)], by = "species_name")[order(species)]
fwrite(par_tab, file.path(OUT, "regeneration_parameters.csv")); print(par_tab[order(-reg_pred)][1:8], digits = 3)
ph <- model$param_history
traj <- rbindlist(lapply(seq_along(ph), function(k) data.table(snapshot = k, epoch = 10L * k,
  theta_median = median(1 / (1 + exp(-as.numeric(ph[[k]]$par_theta_recruits_raw))) * 0 + as.numeric(ph[[k]]$par_theta_recruits_raw)),
  logK_median = median(as.numeric(ph[[k]]$reg_logK)), logK_max = max(as.numeric(ph[[k]]$reg_logK)), logK_min = min(as.numeric(ph[[k]]$reg_logK)))))
fwrite(traj, file.path(OUT, "regeneration_trajectory.csv"))
g2 <- ggplot(melt(traj[, .(epoch, theta_raw_median = theta_median, logK_median, logK_max, logK_min)], id.vars = "epoch"), aes(epoch, value, colour = variable)) + geom_line() + theme_bw() +
  labs(title = "v0.1: recruitment dispersion (raw) and log K over the fit", y = "value (unconstrained / log)")
ggsave(file.path(OUT, "fig_regeneration_trajectory.png"), g2, width = 7, height = 4, dpi = 130, bg = "white")

## ---- 3. change metrics vs persistence --------------------------------------------------------
TOT <- H[state_masked == FALSE, .(ba = sum(ba), ba_m = sum(ba_m), ba_0 = sum(ba_0), n = sum(trees), n_m = sum(trees_m), n_0 = sum(trees_0), reg = sum(reg), reg_m = sum(reg_m)), by = siteID]
r2 <- function(o, p) 1 - sum((o - p)^2) / sum((o - mean(o))^2)
chg <- data.table(quantity = c("dBA (m2/patch)", "dN (stems/patch)", "BA year 5 (model)", "BA year 5 (persistence)", "N year 5 (model)", "N year 5 (persistence)", "recruits / ha / 5 yr"),
                  obs_mean = c(mean(TOT$ba - TOT$ba_0), mean(TOT$n - TOT$n_0), mean(TOT$ba), mean(TOT$ba), mean(TOT$n), mean(TOT$n), mean(TOT$reg)),
                  pred_mean = c(mean(TOT$ba_m - TOT$ba_0), mean(TOT$n_m - TOT$n_0), mean(TOT$ba_m), mean(TOT$ba_0), mean(TOT$n_m), mean(TOT$n_0), mean(TOT$reg_m)),
                  r = c(cor(TOT$ba - TOT$ba_0, TOT$ba_m - TOT$ba_0), cor(TOT$n - TOT$n_0, TOT$n_m - TOT$n_0), cor(TOT$ba, TOT$ba_m), cor(TOT$ba, TOT$ba_0), cor(TOT$n, TOT$n_m), cor(TOT$n, TOT$n_0), cor(TOT$reg, TOT$reg_m)),
                  r2 = c(r2(TOT$ba - TOT$ba_0, TOT$ba_m - TOT$ba_0), r2(TOT$n - TOT$n_0, TOT$n_m - TOT$n_0), r2(TOT$ba, TOT$ba_m), r2(TOT$ba, TOT$ba_0), r2(TOT$n, TOT$n_m), r2(TOT$n, TOT$n_0), r2(TOT$reg, TOT$reg_m)))
fwrite(chg, file.path(OUT, "change_metrics.csv")); print(chg, digits = 3)

## ---- 4. growth by size ----------------------------------------------------------------------
gr <- H[growth_n > 0, .(n = sum(growth_n), obs = weighted.mean(growth, growth_n), pred = weighted.mean(growth_m, growth_n)), by = dcls][order(dcls)][, ratio := pred / obs]
fwrite(gr, file.path(OUT, "growth_by_size.csv")); print(gr, digits = 3)

## ---- 5. loss and parameter trajectories -------------------------------------------------------
lh <- fread(file.path(FIT_DIR, "loss_history.csv"))
fwrite(lh[epoch %in% c(1, 10, 25, 50, 100, 150)], file.path(OUT, "loss_milestones.csv"))
CH <- rbindlist(lapply(2:length(ph), function(k) rbindlist(lapply(names(ph[[k]]), function(b)
  data.table(epoch = 10L * k, block = b, mean_abs_change = mean(abs(as.matrix(ph[[k]][[b]]) - as.matrix(ph[[k - 1]][[b]]))))))))
fwrite(dcast(CH, block ~ epoch, value.var = "mean_abs_change"), file.path(OUT, "parameter_change_by_block.csv"))

## ---- 6. data: fates by management trace, new-tree diameters, back-modelled start ---------------
pl <- fread(file.path(RAW, "PLACETTE.csv"), sep = ";", encoding = "UTF-8", colClasses = list(character = c("GEST", "CSA")), select = c("CAMPAGNE", "IDP", "VISITE", "GEST"))
ar <- fread(file.path(RAW, "ARBRE.csv"), sep = ";", encoding = "UTF-8", colClasses = list(character = c("ESPAR", "VEGET", "VEGET5")), select = c("CAMPAGNE", "IDP", "A", "ESPAR", "VEGET", "VEGET5", "C13", "W"))
ids <- as.integer(site$IDP); ar <- ar[IDP %in% ids]; ar[, dbh := C13 * 100 / pi]
ar <- merge(ar, site[, .(IDP, campaign_t0, campaign_t1)], by = "IDP")
t0 <- ar[CAMPAGNE == campaign_t0 & VEGET == "0" & dbh >= 7.5]; t1 <- ar[CAMPAGNE == campaign_t1]
f <- merge(t0[, .(IDP, A, ESPAR, dbh0 = dbh, W)], t1[, .(IDP, A, VEGET5, dbh1 = dbh)], by = c("IDP", "A"), all.x = TRUE)
f <- merge(f, pl[VISITE == 1L, .(IDP, GEST)], by = "IDP", all.x = TRUE)
f[, fate := fcase(VEGET5 == "0", "survivor", VEGET5 %in% PARAMS$fate_dead, "dead", VEGET5 %in% PARAMS$fate_cut, "cut", default = "other")]
f[is.finite(W) == FALSE, W := c(88.42, 39.29, 14.15)[findInterval(dbh0, c(22.5, 37.5)) + 1L]]
ann <- function(p) 1 - (1 - p)^(1 / P)
fg <- f[fate != "other", .(trees = .N, natural_annual = ann(sum(W * (fate == "dead")) / sum(W)), cut_annual = ann(sum(W * (fate == "cut")) / sum(W)),
                           removal_annual = ann(sum(W * (fate %in% c("dead", "cut"))) / sum(W))), by = GEST][order(GEST)]
fwrite(fg, file.path(OUT, "removal_by_gest.csv")); print(fg, digits = 3)
fs <- f[fate != "other", .(trees = .N, natural_annual = ann(sum(W * (fate == "dead")) / sum(W)), cut_annual = ann(sum(W * (fate == "cut")) / sum(W))), by = .(dcls = cut(dbh0, c(0, PARAMS$dbh_breaks, Inf)))][order(dcls)]
fwrite(fs, file.path(OUT, "removal_by_size.csv")); print(fs, digits = 3)
new1 <- t1[VEGET5 == "" & VEGET == "0" & is.finite(dbh) & dbh >= 7.5]
inc <- f[fate == "survivor" & dbh0 < 12.5, .(ESPAR, inc = dbh1 - dbh0)]
q90 <- inc[, .(n = .N, q90 = quantile(inc, 0.9, na.rm = TRUE)), by = ESPAR]; q_pool <- quantile(inc$inc, 0.9, na.rm = TRUE); q90[n < 30, q90 := q_pool]
new1 <- merge(new1, q90[, .(ESPAR, q90)], by = "ESPAR", all.x = TRUE); new1[is.na(q90), q90 := q_pool]
nt <- new1[, .(n = .N, recruit_by_backmodel = sum(dbh - q90 < 7.5)), by = .(dcls = cut(dbh, c(7.5, 10, 12.5, 15, 22.5, 37.5, Inf), include.lowest = TRUE))][order(dcls)]
fwrite(nt, file.path(OUT, "new_trees_by_dbh.csv")); print(nt)
message("done: ", OUT)
