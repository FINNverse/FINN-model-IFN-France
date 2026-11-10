#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-IFN-France | calibration/v02/02_prepare_fiveyear_v2.R
##
## WHAT   Prepared-data version "fiveyear-v2": the v1 preparation
##        (calibration/02_prepare_fiveyear_v1.R, same raw export, same design
##        expansion, same blocks) with the changes the v0.2 audit
##        (calibration/v02/AUDIT_v02.md) asked for. v1 stays untouched.
##
## WHAT CHANGES AGAINST v1
##  1. TOTAL REMOVAL as the mortality target (the user's rule, BWI lesson 4):
##     a tree alive at the first visit that is dead, fallen, broken OR CUT at the
##     second visit counts as removed; only "not found" (N) is censored. In v1 cut
##     trees were censored and the year-5 state of plots with cuts was masked:
##     that asked FINN for a natural mortality the data cannot separate
##     (pre-emptive and sanitary felling is booked as a cut) and threw away 17 %
##     of the states. Now no state is masked. The management trace recorded at
##     the first visit, GEST (0 none, 1 weak/old, 2 clear signs of exploitation),
##     becomes the covariate `gest`; the audit shows removal 2.0 / 2.7 / 3.2 %/yr
##     for GEST 0 / 1 / 2 against natural death 1.5 / 1.6 / 1.1. Plots without a
##     GEST record are dropped (their cut share is 20 %, they are not missing at
##     random). Natural-only columns (`mort_natural`, `n_died_natural`) and the
##     observed cut share ride along for diagnostics, not for fitting.
##  2. RECRUITS split on the back-modelled start-of-period diameter (BWI lesson
##     2): a tree first tallied at the second visit is a recruit only if it can
##     have crossed the 7.5 cm threshold during the five years, i.e. if
##     dbh_t1 - q90_species < 7.5, with q90_species the 90th percentile of the
##     5-year increment of surviving trees that were below 12.5 cm at the first
##     visit (2.6 cm overall). "New" trees above that were standing trees the
##     first visit missed (or crossed a circle's diameter threshold); they are
##     neither recruits nor part of the start stand and are dropped (counted in
##     PREPARATION.md). v1 used a flat 15 cm cap.
##  3. "OTHER" split in two: "other conifer" (ESPAR codes 5x-7x) and "other
##     broadleaf" (the rest). The v1 "other" (126 codes, 5 % of stems) produced
##     154 of the 193 excess recruits per ha and dominated every bare-ground run.
##  4. tree_obs_dt.csv: one row per TALLIED tree at risk (dbh at the first visit,
##     `trees` = its virtual copies for the light computation, annual relative
##     growth of survivors, removed 0/1), for FINN's per-tree mortality and
##     growth losses (fit(loss_aggregation = ..., tree_data = ...)).
##  5. Pilot sets for the v0.2 pilots: pilot2000 (training blocks, seed 1) and
##     holdout1000 (holdout blocks, seed 11), so every pilot is scored on the
##     same sites.
##  Everything else is v1: campaigns, land-cover and planting/clearing filters,
##  0.0707 ha patch with the W x 0.0707 expansion (largest remainder), tallied
##  sample sizes, annualised rates, WorldClim covariates, 100 km blocks (seed 21).
##
## INPUT  <DATA_DIR>/raw/ign_data/unz/{PLACETTE.csv, ARBRE.csv, espar-cdref13.csv}
##        <WORLDCLIM_DIR>/wc2.1_30s_bio/wc2.1_30s_bio_{1,5,6,12,15}.tif, wc2.1_30s_elev.tif
## OUTPUT <DATA_DIR>/prepared/fiveyear-v2/PREPARATION.md, site_table.csv, site_blocks.csv,
##        holdout_separation.csv, species.csv and, in national/, train/,
##        holdout_blocked/, pilot2000/, holdout1000/ (same layout each): obs_dt.csv,
##        env_dt.csv, init_cohort_dt.csv, tree_obs_dt.csv, species_dt.csv, siteID_dt.csv;
##        calibration/v02/PREPARATION_v2.md, calibration/v02/expected_preparation_v2.csv
## RUN    Rscript calibration/v02/02_prepare_fiveyear_v2.R   (from the repository root, ~6 min)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-0848b60"), .libPaths()))   # FINN 0.2.0, main commit 0848b60
suppressPackageStartupMessages({ library(FINN); library(data.table); library(terra) })

DATA_DIR      <- file.path("..", "..", "FINN-zoo-data", "frnfi")
WORLDCLIM_DIR <- file.path("..", "..", "worldclim")
V02_DIR       <- file.path("calibration", "v02")

PARAMS <- list(
  prep_version    = "fiveyear-v2",
  period          = 5L,
  first_campaigns = 2010:2020,
  csa_keep        = c("1", "2", "3"),
  instp5_drop     = c("P", "S"),
  def5_drop       = c("1", "2", "3", "4", "5"),
  gest_keep       = c("0", "1", "2"),       # management trace at visit 1 = the covariate; blank GEST dropped (cut share 20 %, not at random)
  min_dbh_cm      = 7.5,
  class_breaks_cm = c(22.5, 37.5),
  nominal_w       = c(1e4 / (pi * 6^2), 1e4 / (pi * 9^2), 1e4 / (pi * 15^2)),
  patch_size_ha   = pi * 15^2 / 1e4,
  recruit_small_dbh = 12.5,                 # survivors below this at visit 1 define the small-tree increment distribution
  recruit_quantile  = 0.9,                  # a new tree is a recruit if dbh_t1 - q90(increment) < 7.5 (back-modelled start below the threshold)
  recruit_min_trees = 30L,                  # species with fewer small survivors use the pooled q90
  dbh_growth_thresh = c(-5, 30),
  fate_dead       = c("M", "5", "C", "2", "A", "T", "1"),   # natural loss from the growing stock
  fate_cut        = c("6", "7"),            # cut: REMOVAL, counted in the target (total removal), not censored
  fate_lost       = c("N"),                 # not found: censored (unknown fate)
  species_quantile = 0.95,
  conifer_first_digit = c("5", "6", "7"),   # IGN ESPAR: conifer codes start with 5, 6 or 7
  env_rasters = c(tavg_c = "wc2.1_30s_bio/wc2.1_30s_bio_1.tif", prec_mm = "wc2.1_30s_bio/wc2.1_30s_bio_12.tif",
                  tmax_warm_c = "wc2.1_30s_bio/wc2.1_30s_bio_5.tif", tmin_cold_c = "wc2.1_30s_bio/wc2.1_30s_bio_6.tif",
                  prec_season = "wc2.1_30s_bio/wc2.1_30s_bio_15.tif", elev_m = "wc2.1_30s_elev.tif"),
  block_km        = 100, holdout_block_share = 0.2, block_seed = 21L,
  pilot_n         = 2000L, pilot_seed = 1L,       # pilot2000: the v0.2 pilot training set (~4 s/epoch on one thread)
  holdout_pilot_n = 1000L, holdout_pilot_seed = 11L   # holdout1000: the common pilot evaluation set
)
P   <- PARAMS$period
OUT <- file.path(DATA_DIR, "prepared", PARAMS$prep_version)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE); dir.create(V02_DIR, showWarnings = FALSE, recursive = TRUE)
RAW <- file.path(DATA_DIR, "raw", "ign_data", "unz")
log_lines <- character(0)
note <- function(...) { s <- sprintf(...); message(s); log_lines <<- c(log_lines, s); invisible(s) }

## ---- 1. plots ---------------------------------------------------------------------
pl <- fread(file.path(RAW, "PLACETTE.csv"), sep = ";", encoding = "UTF-8",
            colClasses = list(character = c("CSA", "GEST", "DEF5", "INSTP5", "GRECO", "SER", "DEP")),
            select = c("CAMPAGNE", "IDP", "VISITE", "XL", "YL", "CSA", "DEP", "GRECO", "SER", "GEST", "DEF5", "INSTP5"))
v1 <- pl[VISITE == 1L & CAMPAGNE %in% PARAMS$first_campaigns, .(IDP, campaign_t0 = CAMPAGNE, XL, YL, CSA, DEP, GRECO, SER, GEST)]
v2 <- pl[VISITE == 2L, .(IDP, campaign_t1 = CAMPAGNE, DEF5, INSTP5)]
plots <- merge(v1, v2, by = "IDP")
note("plots first-visited %d-%d with a second visit: %d", min(PARAMS$first_campaigns), max(PARAMS$first_campaigns), nrow(plots))
stopifnot(all(plots$campaign_t1 - plots$campaign_t0 == P))
plots <- plots[CSA %in% PARAMS$csa_keep];          note("  land cover 1/2/3 at visit 1: %d", nrow(plots))
plots <- plots[!INSTP5 %in% PARAMS$instp5_drop];   note("  no planting / sowing in the interval: %d", nrow(plots))
plots <- plots[!DEF5 %in% PARAMS$def5_drop];       note("  no clearing in the interval: %d", nrow(plots))
plots <- plots[GEST %in% PARAMS$gest_keep];        note("  management trace GEST recorded (0/1/2): %d  [v2]", nrow(plots))
plots[, gest := as.integer(GEST)]

## ---- 2. trees ---------------------------------------------------------------------
ar <- fread(file.path(RAW, "ARBRE.csv"), sep = ";", encoding = "UTF-8",
            colClasses = list(character = c("ESPAR", "VEGET", "VEGET5")),
            select = c("CAMPAGNE", "IDP", "A", "ESPAR", "VEGET", "VEGET5", "C13", "W"))
ar <- ar[IDP %in% plots$IDP]
ar <- merge(ar, plots[, .(IDP, campaign_t0, campaign_t1)], by = "IDP")
ar[, dbh := C13 * 100 / pi]
size_class <- function(d) findInterval(d, PARAMS$class_breaks_cm) + 1L
t0 <- ar[CAMPAGNE == campaign_t0]; t1 <- ar[CAMPAGNE == campaign_t1]
stand <- t0[VEGET == "0" & is.finite(dbh)]
n_small <- stand[dbh < PARAMS$min_dbh_cm, .N]
stand <- stand[dbh >= PARAMS$min_dbh_cm]
note("  living visit-1 trees >= 7.5 cm: %d (%d below the threshold dropped)", nrow(stand), n_small)
stand[, W := fifelse(is.finite(W) & W > 0, W, PARAMS$nominal_w[size_class(dbh)])]
stand <- merge(stand[, .(IDP, A, ESPAR, dbh_t0 = dbh, W)], t1[, .(IDP, A, VEGET5, dbh_t1 = dbh)], by = c("IDP", "A"), all.x = TRUE)
stand[, fate := fcase(VEGET5 == "0" & is.finite(dbh_t1), "survivor",
                      VEGET5 %in% PARAMS$fate_dead, "dead",
                      VEGET5 %in% PARAMS$fate_cut,  "cut",
                      VEGET5 %in% PARAMS$fate_lost, "lost",
                      default = "unknown")]
fates <- stand[, .N, by = fate][order(-N)]
note("  fates: %s", paste(sprintf("%s %d", fates$fate, fates$N), collapse = ", "))
## v2: the back-modelled start diameter decides what a "new" tree is
inc <- stand[fate == "survivor" & dbh_t0 < PARAMS$recruit_small_dbh, .(ESPAR, inc = dbh_t1 - dbh_t0)]
q_pool <- quantile(inc$inc, PARAMS$recruit_quantile, na.rm = TRUE)
q_sp <- inc[, .(n = .N, q90 = quantile(inc, PARAMS$recruit_quantile, na.rm = TRUE)), by = ESPAR]
q_sp[n < PARAMS$recruit_min_trees, q90 := q_pool]
new1 <- t1[VEGET5 == "" & VEGET == "0" & is.finite(dbh) & dbh >= PARAMS$min_dbh_cm]
new1 <- merge(new1, q_sp[, .(ESPAR, q90)], by = "ESPAR", all.x = TRUE)
new1[is.na(q90), q90 := q_pool]
new1[, dbh_t0_back := dbh - q90]
rec <- new1[dbh_t0_back < PARAMS$min_dbh_cm]
n_missed <- nrow(new1) - nrow(rec)
note("  new trees at visit 2 (alive, >= 7.5 cm): %d; recruits by the back-modelled start diameter (dbh_t1 - species q%d of the small-tree 5-yr increment, pooled %.2f cm, < 7.5 cm): %d; %d dropped as trees missed at visit 1  [v2]",
     nrow(new1), round(100 * PARAMS$recruit_quantile), q_pool, nrow(rec), n_missed)
note("  dropped new trees by dbh: %s", paste(new1[dbh_t0_back >= PARAMS$min_dbh_cm, .N, by = .(cls = cut(dbh, c(7.5, 10, 12.5, 15, 22.5, Inf)))][order(cls), sprintf("%s %d", cls, N)], collapse = ", "))
rec[, W := fifelse(is.finite(W) & W > 0, W, PARAMS$nominal_w[size_class(dbh)])]

## ---- 3. species: 95 % set, "other" split into conifer / broadleaf -----------------
sp <- fread(file.path(RAW, "espar-cdref13.csv"), sep = ";", encoding = "UTF-8", header = FALSE, skip = 1L,
            col.names = c("espar", "lib_espar", "cd_ref", "lib_cdref"), colClasses = "character")
sp[, espar := sub("^// ", "", espar)]
sp[, ESPAR := fifelse(nchar(espar) == 1L, paste0("0", espar), espar)]
sp[, latin := trimws(sub("\\s+(subsp\\.|var\\.|x)\\s.*$", "", lib_cdref))]
sp[, latin := vapply(strsplit(latin, " "), function(w) paste(w[seq_len(min(2L, length(w)))], collapse = " "), character(1))]
tally <- stand[, .(stems_t0 = .N, plots_t0 = uniqueN(IDP)), by = ESPAR][order(-stems_t0)]
tally <- merge(tally, sp[, .(ESPAR, lib_espar, latin)], by = "ESPAR", all.x = TRUE)[order(-stems_t0)]
tally[is.na(latin), `:=`(latin = paste0("ESPAR_", ESPAR), lib_espar = "(not in espar-cdref13)")]
tally[, cum_share := cumsum(stems_t0) / sum(stems_t0)]
n_sel <- which(tally$cum_share >= PARAMS$species_quantile)[1]
tally[, selected := seq_len(.N) <= n_sel & !startsWith(latin, "ESPAR_")]
tally[, other_class := fifelse(substr(ESPAR, 1, 1) %in% PARAMS$conifer_first_digit, "other conifer", "other broadleaf")]
species_map <- tally[, .(ESPAR, species_name = fifelse(selected, latin, other_class))]
note("species: %d ESPAR codes; %d species cover %.1f %% of stems; the rest -> 'other conifer' (%d codes, %d stems) and 'other broadleaf' (%d codes, %d stems)  [v2]",
     nrow(tally), n_sel, 100 * tally$cum_share[n_sel],
     tally[selected == FALSE & other_class == "other conifer", .N], tally[selected == FALSE & other_class == "other conifer", sum(stems_t0)],
     tally[selected == FALSE & other_class == "other broadleaf", .N], tally[selected == FALSE & other_class == "other broadleaf", sum(stems_t0)])
stand <- merge(stand, species_map, by = "ESPAR"); rec <- merge(rec, species_map, by = "ESPAR", all.x = TRUE)
rec[is.na(species_name), species_name := fifelse(substr(ESPAR, 1, 1) %in% PARAMS$conifer_first_digit, "other conifer", "other broadleaf")]

## ---- 4. virtual trees on the 0.0707 ha patch ---------------------------------------
largest_remainder <- function(x) {
  x <- as.numeric(x); x[!is.finite(x) | x < 0] <- 0
  fl <- floor(x); k <- as.integer(round(sum(x)) - sum(fl))
  if (k > 0L) { ord <- order(x - fl, decreasing = TRUE); idx <- ord[seq_len(min(k, length(x)))]; fl[idx] <- fl[idx] + 1 }
  as.integer(fl)
}
stand[, n_trees := largest_remainder(W * PARAMS$patch_size_ha), by = .(IDP, species_name)]
rec[,   n_trees := largest_remainder(W * PARAMS$patch_size_ha), by = .(IDP, species_name)]
ba_design <- stand[, sum(W * dbh2ba(dbh_t0))]; ba_expanded <- stand[, sum(n_trees * dbh2ba(dbh_t0))] / PARAMS$patch_size_ha
note("expansion: basal area %.1f (design) vs %.1f m2/ha (virtual), %.2f %% off; %d tallied trees expand to zero", ba_design / uniqueN(stand$IDP), ba_expanded / uniqueN(stand$IDP), 100 * (ba_expanded / ba_design - 1), stand[n_trees == 0L, .N])
stand[, treeKey := paste(IDP, A, sep = "_")]; rec[, treeKey := paste(IDP, A, sep = "_")]
expand <- function(dt, year, status, living, dbh_col) {
  out <- dt[n_trees > 0L, .SD[rep(seq_len(.N), n_trees)], .SDcols = c("IDP", "treeKey", "species_name", dbh_col)]
  out[, rep_i := seq_len(.N), by = treeKey]
  out[, .(siteName = as.character(IDP), patchName = paste0(IDP, "_1"), treeName = paste0(treeKey, "#", rep_i),
          year = year, species_name, dbh = get(dbh_col), status = status, living = living)]
}
## v2: cut trees are "dead" rows at year 5 (removed from the stand and counted in the closed cohort); lost trees have no row
tree_dt <- rbind(expand(stand, 0L, "alive", TRUE, "dbh_t0"),
                 expand(stand[fate == "survivor"], P, "alive", TRUE, "dbh_t1"),
                 expand(stand[fate %in% c("dead", "cut")], P, "dead", FALSE, "dbh_t0"),
                 expand(rec, P, "new", TRUE, "dbh"))

## ---- 5. environment (WorldClim + gest) -----------------------------------------------
pts <- unique(plots[, .(IDP, XL, YL, gest)])
v   <- project(vect(pts, geom = c("XL", "YL"), crs = "EPSG:2154"), "EPSG:4326")
pts[, c("x", "y") := as.data.table(crds(v))]
for (nm in names(PARAMS$env_rasters)) { r <- rast(file.path(WORLDCLIM_DIR, PARAMS$env_rasters[[nm]])); pts[, (nm) := extract(r, v, ID = FALSE)[[1]]] }
env_ok <- pts[complete.cases(pts)]
note("env: WorldClim for %d of %d plots", nrow(env_ok), nrow(pts))
env_vars <- c(names(PARAMS$env_rasters), "gest", "x", "y")
tree_dt <- tree_dt[siteName %in% as.character(env_ok$IDP)]
env_dt0 <- merge(unique(tree_dt[, .(siteName, year)]), env_ok[, c(list(siteName = as.character(IDP)), .SD), .SDcols = env_vars], by = "siteName")

## ---- 6. makeObsData -> resolveSiteIDs -------------------------------------------------
obs_list <- makeObsData(tree_dt = tree_dt, plotsize = PARAMS$patch_size_ha, aggregate_by_site = FALSE, minNyears = 2L,
                        fix_period_length = P, dbh_growth_thresh = PARAMS$dbh_growth_thresh, Npatches = NULL, Nspecies = NULL, NspeciesQuantile = NULL)
obs_raw <- obs_list$obs_dt; tree_obs <- obs_list$tree_dt[siteName %in% unique(obs_raw$siteName)]
inputs <- resolveSiteIDs(tree_dt = tree_obs, obs_dt = obs_raw, env_dt = env_dt0, createInitCohorts = FALSE)
obs_dt <- inputs$obs_dt; siteID_dt <- unique(inputs$siteID_dt[, .(siteID, siteName)]); species_dt <- inputs$species_dt[order(species)]
note("makeObsData / resolveSiteIDs: %d sites, %d species classes", uniqueN(obs_dt$siteID), nrow(species_dt))

## ---- 7. year-5 observation: tallied counts, total removal, annualised -------------------
obs_dt <- obs_dt[year == 1L][, year := P][, period_length := P]
obs_dt <- merge(obs_dt, siteID_dt, by = "siteID")
eff <- stand[, {
  at_risk <- fate %in% c("survivor", "dead", "cut"); removed <- fate %in% c("dead", "cut"); died_nat <- fate == "dead"; surv <- fate == "survivor"
  g <- fifelse(surv & is.finite(dbh_t1) & dbh_t0 > 0, dbh_t1 / dbh_t0 - 1, NA_real_)
  .(n_at_risk_t = sum(at_risk),
    removal_w   = if (sum(W[at_risk]) > 0) sum(W[removed]) / sum(W[at_risk]) else NA_real_,
    natural_w   = if (sum(W[at_risk]) > 0) sum(W[died_nat]) / sum(W[at_risk]) else NA_real_,
    cut_w       = if (sum(W[at_risk]) > 0) sum(W[fate == "cut"]) / sum(W[at_risk]) else NA_real_,
    growth_n_t  = sum(is.finite(g)),
    growth_w    = if (any(is.finite(g))) sum(W[is.finite(g)] * g[is.finite(g)]) / sum(W[is.finite(g)]) else NA_real_)
}, by = .(siteName = as.character(IDP), species_name)]
obs_dt <- merge(obs_dt, eff, by = c("siteName", "species_name"), all.x = TRUE)
obs_dt[is.na(n_at_risk_t), n_at_risk_t := 0L]; obs_dt[is.na(growth_n_t), growth_n_t := 0L]
ann <- function(p5) 1 - (1 - pmin(pmax(p5, 0), 1))^(1 / P)
obs_dt[, `:=`(n_at_risk = n_at_risk_t,
              mort         = fifelse(n_at_risk_t > 0 & is.finite(removal_w), ann(removal_w), NA_real_),     # TOTAL REMOVAL, annual
              mort_natural = fifelse(n_at_risk_t > 0 & is.finite(natural_w), ann(natural_w), NA_real_),     # diagnostics only
              cut_share_5yr = fifelse(n_at_risk_t > 0, cut_w, NA_real_),
              growth   = fifelse(growth_n_t > 0 & is.finite(growth_w) & growth_w > -1, (1 + growth_w)^(1 / P) - 1, NA_real_),
              growth_n = growth_n_t)]
obs_dt[, n_died := fifelse(is.finite(mort), mort * n_at_risk, 0)]
obs_dt[, n_died_natural := fifelse(is.finite(mort_natural), mort_natural * n_at_risk, 0)]
obs_cols <- c("siteID", "year", "species", "ba", "dbh", "trees", "growth", "growth_n", "mort", "n_at_risk", "n_died", "reg",
              "species_name", "period_length", "mort_natural", "n_died_natural", "cut_share_5yr")
obs_dt <- obs_dt[order(siteID, species), ..obs_cols]
note("removal target: pooled annual removal %.4f (natural %.4f, cut share of the at-risk weight over 5 yr %.3f)",
     obs_dt[, sum(n_died) / sum(n_at_risk)], obs_dt[, sum(n_died_natural) / sum(n_at_risk)], stand[fate != "lost", sum(W * (fate == "cut")) / sum(W * (fate != "unknown"))])

## ---- 8. env rows years 1..5, initial cohorts, per-tree table --------------------------
env_site <- merge(siteID_dt, env_ok[, c(list(siteName = as.character(IDP)), .SD), .SDcols = env_vars], by = "siteName")
env_dt <- rbindlist(lapply(seq_len(P), function(k) copy(env_site)[, year := k]))[order(siteID, year), c("siteID", "year", env_vars), with = FALSE]
init0 <- inputs$tree_dt[year == 0L & living == TRUE]; init0[, treeKey := sub("#[0-9]+$", "", treeName)]
init_cohort_dt <- init0[, .(trees = .N, dbh = dbh[1], species_name = species_name[1]), by = .(siteID, patchID, species, treeName = treeKey)][, year := 0L]
setcolorder(init_cohort_dt, c("siteID", "patchID", "year", "treeName", "species", "species_name", "dbh", "trees"))
## per-tree table: every tallied tree at risk (survivor, dead, cut), by name, aligned to the set's species index
tree_obs_dt <- merge(stand[fate %in% c("survivor", "dead", "cut") & n_trees > 0L,
                           .(siteName = as.character(IDP), treeName = treeKey, species_name, dbh = dbh_t0, trees = n_trees, fate, dbh_t1)],
                     siteID_dt, by = "siteName")
tree_obs_dt <- merge(tree_obs_dt, species_dt, by = "species_name")
tree_obs_dt[, `:=`(year = P, patchID = 1L, period_length = P,
                   growth = fifelse(fate == "survivor" & is.finite(dbh_t1) & dbh > 0, (dbh_t1 / dbh)^(1 / P) - 1, NA_real_),
                   died = as.numeric(fate %in% c("dead", "cut")), died_natural = as.numeric(fate == "dead"), cut = as.numeric(fate == "cut"))]
tree_obs_dt <- tree_obs_dt[order(siteID, species, treeName), .(siteID, year, patchID, species, species_name, treeName, dbh, trees, growth, died, died_natural, cut, period_length)]
stopifnot(setequal(unique(init_cohort_dt$siteID), unique(obs_dt$siteID)), setequal(unique(env_dt$siteID), unique(obs_dt$siteID)))

## ---- 9. site table, blocks, sets --------------------------------------------------------
site_tab <- merge(siteID_dt, plots[, .(siteName = as.character(IDP), IDP, campaign_t0, campaign_t1, XL, YL, CSA, DEP, GRECO, SER, GEST, gest)], by = "siteName")
site_tab <- merge(site_tab, env_ok[, c(list(siteName = as.character(IDP)), .SD), .SDcols = setdiff(env_vars, "gest")], by = "siteName")
site_tab <- merge(site_tab, stand[, .(n_tallied_t0 = .N, n_cut = sum(fate == "cut"), n_dead = sum(fate == "dead"), n_lost = sum(fate == "lost"),
                                      cut_share_w = sum(W * (fate == "cut")) / sum(W * (fate != "lost"))), by = .(siteName = as.character(IDP))], by = "siteName")
bm <- PARAMS$block_km * 1000
site_tab[, block := sprintf("E%02d_N%02d", floor(XL / bm), floor(YL / bm))]
blocks <- site_tab[, .(n_sites = .N), by = block][order(block)]
set.seed(PARAMS$block_seed); blocks[, role := "train"]; blocks[sample(.N, round(PARAMS$holdout_block_share * .N)), role := "holdout"]
site_tab <- merge(site_tab, blocks[, .(block, block_role = role)], by = "block"); site_tab[, set := block_role]
set.seed(PARAMS$pilot_seed);         pilot_ids <- site_tab[set == "train", sample(siteID, PARAMS$pilot_n)]
set.seed(PARAMS$holdout_pilot_seed); hpilot_ids <- site_tab[set == "holdout", sample(siteID, PARAMS$holdout_pilot_n)]
site_tab[, `:=`(pilot2000 = siteID %in% pilot_ids, holdout1000 = siteID %in% hpilot_ids)]
setorder(site_tab, siteID)
note("blocks: %d of %d held out; %d training, %d holdout sites; pilot2000 / holdout1000 drawn", blocks[role == "holdout", .N], nrow(blocks), site_tab[set == "train", .N], site_tab[set == "holdout", .N])
tr_xy <- as.matrix(site_tab[set == "train", .(XL, YL)]); ho_xy <- as.matrix(site_tab[set == "holdout", .(XL, YL)])
nearest <- unlist(lapply(split(seq_len(nrow(ho_xy)), ceiling(seq_len(nrow(ho_xy)) / 500)), function(idx) {
  d2 <- outer(ho_xy[idx, 1], tr_xy[, 1], "-")^2 + outer(ho_xy[idx, 2], tr_xy[, 2], "-")^2; sqrt(apply(d2, 1, min)) }))
sep <- data.table(stat = c("median_km", "min_km", "q10_km", "q90_km"), value = c(median(nearest), min(nearest), quantile(nearest, 0.1), quantile(nearest, 0.9)) / 1000)
note("holdout separation: median %.1f km, min %.1f km, 10 %% below %.1f km", sep$value[1], sep$value[2], sep$value[3])
fwrite(sep, file.path(OUT, "holdout_separation.csv")); fwrite(site_tab, file.path(OUT, "site_table.csv"))
fwrite(site_tab[, .(siteID, siteName, x, y, block, block_role, set, pilot2000, holdout1000)], file.path(OUT, "site_blocks.csv"))
species_out <- merge(species_dt, tally[selected == TRUE, .(species_name = latin, espar = paste(ESPAR, collapse = "+"), french_name = paste(unique(lib_espar), collapse = " / "),
                                                            stems_t0 = sum(stems_t0), plots_t0 = sum(plots_t0)), by = .(latin)][, latin := NULL], by = "species_name", all.x = TRUE)[order(species)]
for (oc in c("other conifer", "other broadleaf"))
  species_out[species_name == oc, `:=`(espar = paste(tally[selected == FALSE & other_class == oc, ESPAR], collapse = "+"),
                                       french_name = paste("all other", sub("other ", "", oc), "species"),
                                       stems_t0 = tally[selected == FALSE & other_class == oc, sum(stems_t0)], plots_t0 = NA_integer_)]
fwrite(species_out, file.path(OUT, "species.csv"))

write_set <- function(ids, name) {
  d <- file.path(OUT, name); dir.create(d, showWarnings = FALSE, recursive = TRUE)
  ids <- sort(ids); map <- data.table(siteID_national = ids, siteID_new = seq_along(ids))
  renum <- function(dt) { x <- merge(dt[siteID %in% ids], map, by.x = "siteID", by.y = "siteID_national"); x[, siteID := siteID_new][, siteID_new := NULL]; x }
  o <- renum(obs_dt); e <- renum(env_dt); i <- renum(init_cohort_dt); t <- renum(tree_obs_dt)
  s <- merge(site_tab[siteID %in% ids, .(siteID, siteName, IDP, x, y, block, set, gest, cut_share_w)], map, by.x = "siteID", by.y = "siteID_national")
  s[, siteID_national := siteID][, siteID := siteID_new][, siteID_new := NULL]
  fwrite(o[order(siteID, species), ..obs_cols], file.path(d, "obs_dt.csv"))
  fwrite(e[order(siteID, year), c("siteID", "year", env_vars), with = FALSE], file.path(d, "env_dt.csv"))
  fwrite(i[order(siteID, patchID, species), c("siteID", "patchID", "year", "treeName", "species", "species_name", "dbh", "trees"), with = FALSE], file.path(d, "init_cohort_dt.csv"))
  fwrite(t[order(siteID, species, treeName)], file.path(d, "tree_obs_dt.csv"))
  fwrite(species_dt, file.path(d, "species_dt.csv"))
  fwrite(s[order(siteID), .(siteID, siteID_national, siteName, IDP, x, y, block, set, gest, cut_share_w)], file.path(d, "siteID_dt.csv"))
  note("set %-16s %6d sites, %7d obs rows, %7d cohorts, %7d tallied trees", name, length(ids), nrow(o), nrow(i), nrow(t))
}
write_set(site_tab$siteID, "national"); write_set(site_tab[set == "train", siteID], "train"); write_set(site_tab[set == "holdout", siteID], "holdout_blocked")
write_set(site_tab[pilot2000 == TRUE, siteID], "pilot2000"); write_set(site_tab[holdout1000 == TRUE, siteID], "holdout1000")

## ---- 10. records: PREPARATION, expected_preparation -------------------------------------
site_totals <- obs_dt[, .(ba = sum(ba), trees = sum(trees), reg = sum(reg)), by = siteID]
md <- c(sprintf("# Prepared data: %s", PARAMS$prep_version), "",
  sprintf("Written %s by calibration/v02/02_prepare_fiveyear_v2.R (FINN %s). Raw layers as v1 (calibration/PROVENANCE.md).", format(Sys.time(), "%Y-%m-%d %H:%M %Z"), packageVersion("FINN")), "",
  "## Changes against fiveyear-v1", "",
  "1. Mortality target = TOTAL REMOVAL (dead, fallen, broken or cut); not-found censored; no state masked. Covariate `gest` (IGN GEST at visit 1: 0 none, 1 weak/old, 2 clear signs of management). Plots without GEST dropped.",
  "   Natural-only diagnostics: `mort_natural`, `n_died_natural`, `cut_share_5yr` (W-weighted share of the at-risk trees cut).",
  "2. Recruits = new trees whose back-modelled start diameter (dbh_t1 - species q90 of the small-tree 5-yr increment) is below 7.5 cm; other new trees dropped as missed at visit 1.",
  "3. 'other' split into 'other conifer' and 'other broadleaf'.",
  "4. tree_obs_dt.csv per set: tallied trees at risk (siteID, year 5, patchID 1, species, dbh at visit 1, trees = virtual copies, annual relative growth of survivors, died = removal 0/1, died_natural, cut, period_length 5).",
  "5. Sets pilot2000 (seed 1, training blocks) and holdout1000 (seed 11, holdout blocks) for comparable pilots.", "",
  "Everything else as v1 (see calibration/PREPARATION.md): one 0.0707 ha patch per plot, W x 0.0707 virtual trees, tallied n_at_risk / growth_n, annualised rates, reg per ha per 5 yr, env years 1..5, WorldClim 2.1 normals in natural units, 100 km blocks (seed 21).", "",
  "## Processing log", "", paste0("- ", log_lines), "",
  "## Summary (national)", "",
  sprintf("- %d sites; mean basal area at year 5 %.1f m2/ha, %.0f stems/ha; pooled annual removal %.4f (natural %.4f); mean recruitment %.1f stems/ha per 5 yr",
          uniqueN(obs_dt$siteID), mean(site_totals$ba) / PARAMS$patch_size_ha, mean(site_totals$trees) / PARAMS$patch_size_ha,
          obs_dt[, sum(n_died) / sum(n_at_risk)], obs_dt[, sum(n_died_natural) / sum(n_at_risk)], mean(site_totals$reg)), "",
  "## Species", "", "| species | name | ESPAR | French | stems at visit 1 | plots |", "|---|---|---|---|---|---|",
  species_out[, sprintf("| %d | %s | %s | %s | %d | %s |", species, species_name, espar, french_name, stems_t0, ifelse(is.na(plots_t0), "", plots_t0))], "",
  "## PARAMS", "", "```r", capture.output(str(PARAMS, give.attr = FALSE)), "```")
writeLines(md, file.path(OUT, "PREPARATION.md")); writeLines(md, file.path(V02_DIR, "PREPARATION_v2.md"))
per_set <- rbindlist(lapply(c("national", "train", "holdout_blocked", "pilot2000", "holdout1000"), function(s) {
  o <- fread(file.path(OUT, s, "obs_dt.csv")); i <- fread(file.path(OUT, s, "init_cohort_dt.csv")); t <- fread(file.path(OUT, s, "tree_obs_dt.csv"))
  st <- o[, .(ba = sum(ba), trees = sum(trees), reg = sum(reg)), by = siteID]
  data.table(set = s, sites = uniqueN(o$siteID), obs_rows = nrow(o), cohorts = nrow(i), tallied_trees = nrow(t), virtual_trees_t0 = sum(i$trees),
             ba_ha_t5_mean = round(mean(st$ba) / PARAMS$patch_size_ha, 3), trees_ha_t5_mean = round(mean(st$trees) / PARAMS$patch_size_ha, 2),
             growth_annual_mean = round(weighted.mean(o$growth, o$growth_n, na.rm = TRUE), 6), n_at_risk = sum(o$n_at_risk),
             n_died = round(sum(o$n_died), 2), removal_annual_pooled = round(sum(o$n_died) / sum(o$n_at_risk), 6),
             natural_annual_pooled = round(sum(o$n_died_natural) / sum(o$n_at_risk), 6), reg_ha_5yr_mean = round(mean(st$reg), 4))
}))
per_species <- obs_dt[, .(set = "national", sites_present_t5 = sum(trees > 0), ba_sum = round(sum(ba), 3), n_at_risk = sum(n_at_risk),
                          n_died = round(sum(n_died), 2), n_died_natural = round(sum(n_died_natural), 2), reg_sum = round(sum(reg), 2)), by = .(species, species_name)]
expected <- rbind(melt(per_set, id.vars = "set", variable.name = "quantity", value.name = "value", variable.factor = FALSE)[, `:=`(species = NA_integer_, species_name = NA_character_)],
                  melt(per_species, id.vars = c("set", "species", "species_name"), variable.name = "quantity", value.name = "value", variable.factor = FALSE), use.names = TRUE)
fwrite(expected[, .(set, species, species_name, quantity, value)], file.path(V02_DIR, "expected_preparation_v2.csv"))
message("done: ", OUT)
