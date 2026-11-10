#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-zoo | IGN-FR_process/scripts/02_prepare_fiveyear_v1.R
##
## WHAT   Turns the IGN raw tables (scripts/01) into FINN's calibration tables
##        for the 5-year remeasurement interval of the French NFI:
##        obs_dt (site x year x species: ba, dbh, trees, growth, mortality
##        counts, recruitment), env_dt (WorldClim 2.1 normals per site and
##        year), init_cohort_dt (the first-visit stand), species_dt, siteID_dt,
##        a spatially blocked train / holdout split and PREPARATION.md.
##        Prepared-data version "fiveyear-v1".
##
## THE DATA
##   Since 2010 every IGN plot is revisited exactly five campaigns after its
##   first visit (same IDP, same tree id A). 2010-2014 revisits only recorded
##   whether a tree was cut; from the 2015 campaign on every first-visit tree is
##   remeasured (C13) or booked dead / cut / not found, and new trees >= 7.5 cm
##   dbh are tallied. So the usable cohort is first visits 2010-2020, second
##   visits 2015-2025 (one interval per plot).
##
## THE DESIGN AND HOW IT BECOMES A FINN PATCH
##   Trees are tallied on three concentric circles: 6 m (all trees, dbh >= 7.5
##   cm), 9 m (dbh >= 22.5), 15 m (dbh >= 37.5). There are no tree positions, but
##   every tree carries W = stems per hectare it stands for (10000 / circle
##   area, with IGN's edge corrections). FINN needs fixed-area patches with
##   integer tree counts, so -- as FINN-bwi does for angle-count plots -- the
##   plot is declared ONE patch of 0.0707 ha (the 15 m circle) and every tallied
##   tree is expanded into round(W * 0.0707) identical "virtual" trees (6.25 /
##   2.78 / 1.0 for the three size classes; largest-remainder rounding within
##   plot x species, so species totals are preserved). Density and basal area
##   per hectare are thereby exact; the sample sizes behind the rates are NOT
##   (one measured tree is not six observations), so n_at_risk and growth_n are
##   reset to the number of TALLIED trees while the rates stay design-weighted.
##
## FATES (VEGET5 at the second visit) AND WHAT FINN SEES
##   0        survivor: growth from C13 at both visits
##   M 5 C    dead standing; 2 A dead/fallen (chablis); T broken below 1.3 m;
##   1        fallen but alive (chablis vivant) -> all counted as NATURAL DEATH
##            (the tree has left the growing stock)
##   6 7      cut -> CENSORED from the mortality cohort (a competing risk FINN
##            has no process for); the plot's year-5 STATE targets (ba, dbh,
##            trees) are set to NA because the model cannot reproduce a harvest
##   N        not found -> censored, state masked like a cut
##   (none)   first-visit tree without a second-visit row -> censored
##   new trees (VEGET5 empty, alive): recruitment; a "new" tree larger than
##   recruit_max_dbh cannot have crossed 7.5 cm in five years, it was missed at
##   the first visit, and is dropped from everything (counted in PREPARATION.md)
##
## TIME   One interval per plot: year 0 = first visit, year 5 = second visit,
##        period_length = 5. Growth and mortality are annualised (FINN scores
##        the mean of its annual rates over the window), recruitment is the
##        5-year total per hectare (FINN sums). Env rows for years 1..5 (static
##        WorldClim normals; no year-0 row, as in FINN-fia).
##
## INPUT  <DATA_DIR>/raw/ign_data/unz/{PLACETTE.csv, ARBRE.csv, espar-cdref13.csv}
##        <WORLDCLIM_DIR>/wc2.1_30s_bio/wc2.1_30s_bio_{1,5,6,12,15}.tif, wc2.1_30s_elev.tif
## OUTPUT <DATA_DIR>/prepared/fiveyear-v1/PREPARATION.md, site_table.csv, site_blocks.csv,
##        holdout_separation.csv, species.csv and, in national/, train/,
##        holdout_blocked/, pilot/ (same layout each): obs_dt.csv, env_dt.csv,
##        init_cohort_dt.csv, species_dt.csv, siteID_dt.csv
## RUN    Rscript IGN-FR_process/scripts/02_prepare_fiveyear_v1.R   (from the repo root, ~10 min)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c(path.expand("~/R/lib-finn-0848b60"), .libPaths()))   # FINN 0.2.0, main commit 0848b60
suppressPackageStartupMessages({ library(FINN); library(data.table); library(terra) })

DATA_DIR      <- file.path("..", "..", "FINN-zoo-data", "frnfi")   # zoo layout: data outside git
WORLDCLIM_DIR <- file.path("..", "..", "worldclim")                # WorldClim 2.1 30 s rasters

PARAMS <- list(
  prep_version    = "fiveyear-v1",
  period          = 5L,                     # IGN revisit interval, campaigns
  first_campaigns = 2010:2020,              # first visits whose revisit (2015-2025) remeasured C13
  csa_keep        = c("1", "2", "3"),       # land cover at visit 1: closed forest, bosquet, open forest; drops poplar plantations (5)
  instp5_drop     = c("P", "S"),            # planting / sowing between the visits = a regenerated stand, not natural dynamics
  def5_drop       = c("1", "2", "3", "4", "5"),   # any clearing of the plot between the visits
  min_dbh_cm      = 7.5,                    # IGN tally threshold (C13 >= 0.235 m)
  class_breaks_cm = c(22.5, 37.5),          # PB / MB / GB size classes of the concentric design
  nominal_w       = c(1e4 / (pi * 6^2), 1e4 / (pi * 9^2), 1e4 / (pi * 15^2)),   # stems/ha per tree when W is missing: 88.42, 39.29, 14.15
  patch_size_ha   = pi * 15^2 / 1e4,        # 0.0707 ha, the 15 m circle = the plot as one FINN patch
  recruit_max_dbh = 15,                     # a new tree above this (> 1.5 cm/yr from 7.5 cm) was missed at visit 1, not ingrowth
  dbh_growth_thresh = c(-5, 30),            # cm over 5 yr; a site with any tree outside is dropped (makeObsData)
  fate_dead       = c("M", "5", "C", "2", "A", "T", "1"),   # natural loss from the growing stock (see header)
  fate_cut        = c("6", "7"),            # harvested: censored, plot state masked
  fate_lost       = c("N"),                 # not found: censored, plot state masked
  species_quantile = 0.95,                  # species covering 95 % of living tallied stems at visit 1; the rest -> "other"
  ## env: WorldClim 2.1 normals 1970-2000, 30 arc-seconds. Model covariates are
  ## tavg_c and prec_mm (the FINN-fia names, so regions can be merged); the
  ## others ride along for other zoo models.
  env_rasters = c(tavg_c = "wc2.1_30s_bio/wc2.1_30s_bio_1.tif", prec_mm = "wc2.1_30s_bio/wc2.1_30s_bio_12.tif",
                  tmax_warm_c = "wc2.1_30s_bio/wc2.1_30s_bio_5.tif", tmin_cold_c = "wc2.1_30s_bio/wc2.1_30s_bio_6.tif",
                  prec_season = "wc2.1_30s_bio/wc2.1_30s_bio_15.tif", elev_m = "wc2.1_30s_elev.tif"),
  ## spatially blocked holdout: 100 km squares on the Lambert-93 grid (metres),
  ## a fifth of the blocks held out; training sites come only from the rest
  block_km        = 100, holdout_block_share = 0.2, block_seed = 21L,
  pilot_n         = 500L, pilot_seed = 1L   # a small training subset for smoke tests
)
P   <- PARAMS$period
OUT <- file.path(DATA_DIR, "prepared", PARAMS$prep_version)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
RAW <- file.path(DATA_DIR, "raw", "ign_data", "unz")
log_lines <- character(0)
note <- function(...) { s <- sprintf(...); message(s); log_lines <<- c(log_lines, s); invisible(s) }

## ---- 1. plots: one interval per IDP ----------------------------------------------
pl <- fread(file.path(RAW, "PLACETTE.csv"), sep = ";", encoding = "UTF-8",
            colClasses = list(character = c("CSA", "GEST", "DEF5", "INSTP5", "GRECO", "SER", "DEP")),
            select = c("CAMPAGNE", "IDP", "VISITE", "XL", "YL", "CSA", "DEP", "GRECO", "SER", "GEST", "DEF5", "INSTP5"))
v1 <- pl[VISITE == 1L & CAMPAGNE %in% PARAMS$first_campaigns,
         .(IDP, campaign_t0 = CAMPAGNE, XL, YL, CSA, DEP, GRECO, SER, GEST)]
v2 <- pl[VISITE == 2L, .(IDP, campaign_t1 = CAMPAGNE, DEF5, INSTP5)]
plots <- merge(v1, v2, by = "IDP")
note("plots first-visited %d-%d with a second visit: %d", min(PARAMS$first_campaigns), max(PARAMS$first_campaigns), nrow(plots))
stopifnot(all(plots$campaign_t1 - plots$campaign_t0 == P))
plots <- plots[CSA %in% PARAMS$csa_keep];           note("  land cover 1/2/3 at visit 1 (no poplar plantations): %d", nrow(plots))
plots <- plots[!INSTP5 %in% PARAMS$instp5_drop];    note("  no planting / sowing in the interval: %d", nrow(plots))
plots <- plots[!DEF5 %in% PARAMS$def5_drop];        note("  no clearing in the interval: %d", nrow(plots))

## ---- 2. trees: visit-1 stand, visit-2 fates and recruits ------------------------
ar <- fread(file.path(RAW, "ARBRE.csv"), sep = ";", encoding = "UTF-8",
            colClasses = list(character = c("ESPAR", "VEGET", "VEGET5")),
            select = c("CAMPAGNE", "IDP", "A", "ESPAR", "VEGET", "VEGET5", "C13", "W"))
ar <- ar[IDP %in% plots$IDP]
ar <- merge(ar, plots[, .(IDP, campaign_t0, campaign_t1)], by = "IDP")
ar[, dbh := C13 * 100 / pi]                              # circumference (m) -> diameter (cm)
size_class <- function(d) findInterval(d, PARAMS$class_breaks_cm) + 1L   # 1 PB, 2 MB, 3 GB

t0 <- ar[CAMPAGNE == campaign_t0]
t1 <- ar[CAMPAGNE == campaign_t1]
note("tree rows at visit 1: %d; at visit 2: %d", nrow(t0), nrow(t1))
## the first-visit stand: living trees at or above the tally threshold
stand <- t0[VEGET == "0" & is.finite(dbh)]
note("  living at visit 1: %d (dropped: %d dead / fallen, %d without C13)",
     nrow(stand), nrow(t0[VEGET != "0"]), nrow(t0[VEGET == "0" & !is.finite(dbh)]))
n_small <- stand[dbh < PARAMS$min_dbh_cm, .N]
stand <- stand[dbh >= PARAMS$min_dbh_cm];  note("  below the 7.5 cm threshold, dropped: %d", n_small)
stand[, W := fifelse(is.finite(W) & W > 0, W, PARAMS$nominal_w[size_class(dbh)])]
## fates from the second-visit row of the same tree (IDP, A)
stand <- merge(stand[, .(IDP, A, ESPAR, dbh_t0 = dbh, W)],
               t1[, .(IDP, A, VEGET5, dbh_t1 = dbh)], by = c("IDP", "A"), all.x = TRUE)
stand[, fate := fcase(
  VEGET5 == "0" & is.finite(dbh_t1), "survivor",
  VEGET5 %in% PARAMS$fate_dead,      "dead",
  VEGET5 %in% PARAMS$fate_cut,       "cut",
  VEGET5 %in% PARAMS$fate_lost,      "lost",
  default = "unknown")]                                  # no visit-2 row, or alive without a remeasured C13
fates <- stand[, .N, by = fate][order(-N)]
note("  fates: %s", paste(sprintf("%s %d", fates$fate, fates$N), collapse = ", "))
## recruits: trees first tallied at the second visit, alive, within the ingrowth range
rec <- t1[VEGET5 == "" & VEGET == "0" & is.finite(dbh) & dbh >= PARAMS$min_dbh_cm]
n_rec_big <- rec[dbh > PARAMS$recruit_max_dbh, .N]
rec <- rec[dbh <= PARAMS$recruit_max_dbh]
rec[, W := fifelse(is.finite(W) & W > 0, W, PARAMS$nominal_w[size_class(dbh)])]
note("  new trees at visit 2 (alive, >= 7.5 cm): %d kept as recruits, %d above %g cm dropped as missed at visit 1",
     nrow(rec), n_rec_big, PARAMS$recruit_max_dbh)

## ---- 3. species: ESPAR -> Latin name, 95 % set ----------------------------------
sp <- fread(file.path(RAW, "espar-cdref13.csv"), sep = ";", encoding = "UTF-8", header = FALSE, skip = 1L,
            col.names = c("espar", "lib_espar", "cd_ref", "lib_cdref"), colClasses = "character")
sp[, espar := sub("^// ", "", espar)]
## ARBRE pads one-digit codes with a zero ("09"), the lookup does not ("9")
sp[, ESPAR := fifelse(nchar(espar) == 1L, paste0("0", espar), espar)]
## a Latin binomial: drop subspecies / variety suffixes so e.g. "Picea abies subsp. abies" -> "Picea abies"
sp[, latin := trimws(sub("\\s+(subsp\\.|var\\.|x)\\s.*$", "", lib_cdref))]
sp[, latin := vapply(strsplit(latin, " "), function(w) paste(w[seq_len(min(2L, length(w)))], collapse = " "), character(1))]
tally <- stand[, .(stems_t0 = .N, plots_t0 = uniqueN(IDP)), by = ESPAR][order(-stems_t0)]
tally <- merge(tally, sp[, .(ESPAR, lib_espar, latin)], by = "ESPAR", all.x = TRUE)[order(-stems_t0)]
## a few rare codes are not in the lookup (e.g. 29AF, 68CE): they can only be "other"
n_unnamed <- tally[is.na(latin), .N]
tally[is.na(latin), `:=`(latin = paste0("ESPAR_", ESPAR), lib_espar = "(not in espar-cdref13)")]
tally[, cum_share := cumsum(stems_t0) / sum(stems_t0)]
n_sel <- which(tally$cum_share >= PARAMS$species_quantile)[1]
tally[, selected := seq_len(.N) <= n_sel & !startsWith(latin, "ESPAR_")]
note("  %d ESPAR codes without a lookup entry (%d stems) lumped into 'other'", n_unnamed, tally[startsWith(latin, "ESPAR_"), sum(stems_t0)])
species_map <- tally[, .(ESPAR, species_name = fifelse(selected, latin, "other"))]
note("species: %d ESPAR codes among living visit-1 trees; %d species cover %.1f %% of stems, %d codes lumped into 'other'",
     nrow(tally), n_sel, 100 * tally$cum_share[n_sel], nrow(tally) - n_sel)
## two ESPAR codes may share a binomial (e.g. Pinus nigra var. corsicana / subsp. nigra): both map to it
stand <- merge(stand, species_map, by = "ESPAR"); rec <- merge(rec, species_map, by = "ESPAR", all.x = TRUE)
rec[is.na(species_name), species_name := "other"]

## ---- 4. virtual trees on a 0.0707 ha patch --------------------------------------
largest_remainder <- function(x) {      # integers that preserve round(sum(x)) within the group
  x <- as.numeric(x); x[!is.finite(x) | x < 0] <- 0
  fl <- floor(x); k <- as.integer(round(sum(x)) - sum(fl))
  if (k > 0L) { ord <- order(x - fl, decreasing = TRUE); fl[ord[seq_len(min(k, length(x)))]] <- fl[ord[seq_len(min(k, length(x)))]] + 1 }
  as.integer(fl)
}
stand[, n_trees := largest_remainder(W * PARAMS$patch_size_ha), by = .(IDP, species_name)]
rec[,   n_trees := largest_remainder(W * PARAMS$patch_size_ha), by = .(IDP, species_name)]
ba_design   <- stand[, sum(W * dbh2ba(dbh_t0))]
ba_expanded <- stand[, sum(n_trees * dbh2ba(dbh_t0))] / PARAMS$patch_size_ha
note("expansion: visit-1 basal area %.1f m2/ha design-weighted vs %.1f m2/ha from the virtual trees (%.2f %% off); %d of %d tallied trees expand to zero",
     ba_design / uniqueN(stand$IDP), ba_expanded / uniqueN(stand$IDP), 100 * (ba_expanded / ba_design - 1), stand[n_trees == 0L, .N], nrow(stand))
stand[, treeKey := paste(IDP, A, sep = "_")]
rec[,   treeKey := paste(IDP, A, sep = "_")]
expand <- function(dt, year, status, living, dbh_col) {
  out <- dt[n_trees > 0L, .SD[rep(seq_len(.N), n_trees)], .SDcols = c("IDP", "treeKey", "species_name", dbh_col)]
  out[, rep_i := seq_len(.N), by = treeKey]
  out[, .(siteName = as.character(IDP), patchName = paste0(IDP, "_1"), treeName = paste0(treeKey, "#", rep_i),
          year = year, species_name, dbh = get(dbh_col), status = status, living = living)]
}
tree_dt <- rbind(
  expand(stand,                       0L, "alive", TRUE,  "dbh_t0"),
  expand(stand[fate == "survivor"],   P,  "alive", TRUE,  "dbh_t1"),
  expand(stand[fate == "dead"],       P,  "dead",  FALSE, "dbh_t0"),
  expand(rec,                         P,  "new",   TRUE,  "dbh"))
## cut, lost and unknown trees have no year-5 row: they leave the closed cohort
## (censored) and the year-5 stand. A plot whose every tree is censored has no
## year-5 row at all and drops out in makeObsData (minNyears).

## ---- 5. environment: WorldClim normals at the (fuzzed) plot centre --------------
pts <- unique(plots[, .(IDP, XL, YL)])
v   <- project(vect(pts, geom = c("XL", "YL"), crs = "EPSG:2154"), "EPSG:4326")
pts[, c("x", "y") := as.data.table(crds(v))]
for (nm in names(PARAMS$env_rasters)) {
  r <- rast(file.path(WORLDCLIM_DIR, PARAMS$env_rasters[[nm]]))
  pts[, (nm) := extract(r, v, ID = FALSE)[[1]]]
}
env_ok <- pts[complete.cases(pts)]
note("env: WorldClim extracted for %d of %d plots (%d without a land pixel dropped)", nrow(env_ok), nrow(pts), nrow(pts) - nrow(env_ok))
env_vars <- c(names(PARAMS$env_rasters), "x", "y")
tree_dt <- tree_dt[siteName %in% as.character(env_ok$IDP)]
env_dt0 <- merge(unique(tree_dt[, .(siteName, year)]),
                 env_ok[, c(list(siteName = as.character(IDP)), .SD), .SDcols = env_vars], by = "siteName")

## ---- 6. FINN's makeObsData -> resolveSiteIDs --------------------------------------
obs_list <- makeObsData(tree_dt = tree_dt, plotsize = PARAMS$patch_size_ha, aggregate_by_site = FALSE,
                        minNyears = 2L, fix_period_length = P, dbh_growth_thresh = PARAMS$dbh_growth_thresh,
                        Npatches = NULL, Nspecies = NULL, NspeciesQuantile = NULL)
obs_raw  <- obs_list$obs_dt
tree_obs <- obs_list$tree_dt[siteName %in% unique(obs_raw$siteName)]
inputs   <- resolveSiteIDs(tree_dt = tree_obs, obs_dt = obs_raw, env_dt = env_dt0, createInitCohorts = FALSE)
obs_dt    <- inputs$obs_dt
siteID_dt <- unique(inputs$siteID_dt[, .(siteID, siteName)])
species_dt <- inputs$species_dt[order(species)]
note("makeObsData / resolveSiteIDs: %d sites, %d species incl. 'other'", uniqueN(obs_dt$siteID), nrow(species_dt))

## ---- 7. the year-5 observation: tallied sample sizes, annualised rates, masks ---
obs_dt <- obs_dt[year == 1L][, year := P]          # resolveSiteIDs numbers inventories 0, 1; the second visit is year 5
obs_dt[, period_length := P]
obs_dt <- merge(obs_dt, siteID_dt, by = "siteID")
## design-weighted rates and TALLIED counts per plot x species
eff <- stand[, {
  at_risk <- fate %in% c("survivor", "dead"); died <- fate == "dead"; surv <- fate == "survivor"
  g <- fifelse(surv & is.finite(dbh_t1) & dbh_t0 > 0, dbh_t1 / dbh_t0 - 1, NA_real_)
  .(n_at_risk_t = sum(at_risk),
    mort_w      = if (sum(W[at_risk]) > 0) sum(W[died]) / sum(W[at_risk]) else NA_real_,
    growth_n_t  = sum(is.finite(g)),
    growth_w    = if (any(is.finite(g))) sum(W[is.finite(g)] * g[is.finite(g)]) / sum(W[is.finite(g)]) else NA_real_)
}, by = .(siteName = as.character(IDP), species_name)]
obs_dt <- merge(obs_dt, eff, by = c("siteName", "species_name"), all.x = TRUE)
obs_dt[is.na(n_at_risk_t), n_at_risk_t := 0L]; obs_dt[is.na(growth_n_t), growth_n_t := 0L]
## annualise: FINN compares the mean of its annual rates over the 5-year window
obs_dt[, `:=`(
  n_at_risk = n_at_risk_t,
  mort      = fifelse(n_at_risk_t > 0 & is.finite(mort_w), 1 - (1 - pmin(pmax(mort_w, 0), 1))^(1 / P), NA_real_),
  growth    = fifelse(growth_n_t > 0 & is.finite(growth_w) & growth_w > -1, (1 + growth_w)^(1 / P) - 1, NA_real_),
  growth_n  = growth_n_t)]
obs_dt[, n_died := fifelse(is.finite(mort), mort * n_at_risk, 0)]
## plots with a cut or a lost tree: the year-5 state is not what FINN can simulate
masked_plots <- stand[fate %in% c("cut", "lost"), unique(as.character(IDP))]
obs_dt[, state_masked := siteName %in% masked_plots]
obs_dt[state_masked == TRUE, `:=`(ba = NA_real_, dbh = NA_real_, trees = NA_real_)]
note("year-5 state (ba, dbh, trees) masked on %d of %d sites with a cut or not-found tree; growth, mortality and recruitment kept",
     obs_dt[state_masked == TRUE, uniqueN(siteID)], uniqueN(obs_dt$siteID))
obs_cols <- c("siteID", "year", "species", "ba", "dbh", "trees", "growth", "growth_n", "mort", "n_at_risk", "n_died", "reg",
              "species_name", "period_length", "state_masked")
obs_dt <- obs_dt[order(siteID, species), ..obs_cols]

## ---- 8. env rows years 1..5, initial cohorts -------------------------------------
env_site <- merge(siteID_dt, env_ok[, c(list(siteName = as.character(IDP)), .SD), .SDcols = env_vars], by = "siteName")
env_dt <- rbindlist(lapply(seq_len(P), function(k) copy(env_site)[, year := k]))
env_dt <- env_dt[order(siteID, year), c("siteID", "year", env_vars), with = FALSE]
## one cohort per TALLIED tree with trees = its virtual copies (cheaper than one cohort per copy, same stand)
init0 <- inputs$tree_dt[year == 0L & living == TRUE]
init0[, treeKey := sub("#[0-9]+$", "", treeName)]
init_cohort_dt <- init0[, .(trees = .N, dbh = dbh[1], species_name = species_name[1]),
                        by = .(siteID, patchID, species, treeName = treeKey)][, year := 0L]
setcolorder(init_cohort_dt, c("siteID", "patchID", "year", "treeName", "species", "species_name", "dbh", "trees"))
stopifnot(setequal(unique(init_cohort_dt$siteID), unique(obs_dt$siteID)), setequal(unique(env_dt$siteID), unique(obs_dt$siteID)))

## ---- 9. site table and the spatially blocked split ------------------------------
site_tab <- merge(siteID_dt, plots[, .(siteName = as.character(IDP), IDP, campaign_t0, campaign_t1, XL, YL, CSA, DEP, GRECO, SER, GEST)], by = "siteName")
site_tab <- merge(site_tab, env_ok[, c(list(siteName = as.character(IDP)), .SD), .SDcols = env_vars], by = "siteName")
site_tab <- merge(site_tab, stand[, .(n_tallied_t0 = .N, n_cut = sum(fate == "cut"), n_dead = sum(fate == "dead")), by = .(siteName = as.character(IDP))], by = "siteName")
site_tab[, state_masked := siteName %in% masked_plots]
bm <- PARAMS$block_km * 1000
site_tab[, block := sprintf("E%02d_N%02d", floor(XL / bm), floor(YL / bm))]
blocks <- site_tab[, .(n_sites = .N), by = block][order(block)]
set.seed(PARAMS$block_seed)
blocks[, role := "train"]
blocks[sample(.N, round(PARAMS$holdout_block_share * .N)), role := "holdout"]
site_tab <- merge(site_tab, blocks[, .(block, block_role = role)], by = "block")
site_tab[, set := block_role]
set.seed(PARAMS$pilot_seed)
pilot_ids <- site_tab[set == "train", sample(siteID, PARAMS$pilot_n)]
site_tab[, pilot := siteID %in% pilot_ids]
setorder(site_tab, siteID)
note("blocks: %d of %d %d-km blocks held out; %d training sites, %d holdout sites",
     blocks[role == "holdout", .N], nrow(blocks), PARAMS$block_km, site_tab[set == "train", .N], site_tab[set == "holdout", .N])
## separation: distance from every holdout site to the nearest training site (chunked, metres on Lambert-93)
tr_xy <- as.matrix(site_tab[set == "train", .(XL, YL)]); ho_xy <- as.matrix(site_tab[set == "holdout", .(XL, YL)])
nearest <- unlist(lapply(split(seq_len(nrow(ho_xy)), ceiling(seq_len(nrow(ho_xy)) / 500)), function(idx) {
  d2 <- outer(ho_xy[idx, 1], tr_xy[, 1], "-")^2 + outer(ho_xy[idx, 2], tr_xy[, 2], "-")^2
  sqrt(apply(d2, 1, min))
}))
sep <- data.table(stat = c("median_km", "min_km", "q10_km", "q90_km"),
                  value = c(median(nearest), min(nearest), quantile(nearest, 0.1), quantile(nearest, 0.9)) / 1000)
note("holdout separation: nearest training site median %.1f km, minimum %.1f km, 10 %% below %.1f km", sep$value[1], sep$value[2], sep$value[3])
fwrite(sep, file.path(OUT, "holdout_separation.csv"))
fwrite(site_tab, file.path(OUT, "site_table.csv"))
fwrite(site_tab[, .(siteID, siteName, x, y, block, block_role, set, pilot)], file.path(OUT, "site_blocks.csv"))
species_out <- merge(species_dt, tally[selected == TRUE, .(species_name = latin, espar = paste(ESPAR, collapse = "+"),
                                                            french_name = paste(unique(lib_espar), collapse = " / "),
                                                            stems_t0 = sum(stems_t0), plots_t0 = sum(plots_t0)), by = .(latin)][, latin := NULL],
                     by = "species_name", all.x = TRUE)[order(species)]
species_out[species_name == "other", `:=`(espar = paste(tally[selected == FALSE, ESPAR], collapse = "+"),
                                          french_name = "all other species", stems_t0 = tally[selected == FALSE, sum(stems_t0)],
                                          plots_t0 = NA_integer_)]
fwrite(species_out, file.path(OUT, "species.csv"))

## ---- 10. write the sets (siteIDs renumbered 1..n within each set) ---------------
write_set <- function(ids, name) {
  d <- file.path(OUT, name); dir.create(d, showWarnings = FALSE, recursive = TRUE)
  ids <- sort(ids); map <- data.table(siteID_national = ids, siteID = seq_along(ids))
  renum <- function(dt) { dt <- merge(dt, map, by.x = "siteID", by.y = "siteID_national"); dt[, siteID := NULL]; setnames(dt, "siteID.y", "siteID") ; dt }
  o <- merge(obs_dt[siteID %in% ids], map, by.x = "siteID", by.y = "siteID_national")[, siteID := siteID.y][, siteID.y := NULL]
  e <- merge(env_dt[siteID %in% ids], map, by.x = "siteID", by.y = "siteID_national")[, siteID := siteID.y][, siteID.y := NULL]
  i <- merge(init_cohort_dt[siteID %in% ids], map, by.x = "siteID", by.y = "siteID_national")[, siteID := siteID.y][, siteID.y := NULL]
  s <- merge(site_tab[siteID %in% ids, .(siteID, siteName, IDP, x, y, block, set, state_masked)], map, by.x = "siteID", by.y = "siteID_national")
  s[, siteID_national := siteID][, siteID := siteID.y][, siteID.y := NULL]
  fwrite(o[order(siteID, species), ..obs_cols], file.path(d, "obs_dt.csv"))
  fwrite(e[order(siteID, year), c("siteID", "year", env_vars), with = FALSE], file.path(d, "env_dt.csv"))
  fwrite(i[order(siteID, patchID, species), c("siteID", "patchID", "year", "treeName", "species", "species_name", "dbh", "trees"), with = FALSE],
         file.path(d, "init_cohort_dt.csv"))
  fwrite(species_dt, file.path(d, "species_dt.csv"))
  fwrite(s[order(siteID), .(siteID, siteID_national, siteName, IDP, x, y, block, set, state_masked)], file.path(d, "siteID_dt.csv"))
  note("set %-16s %6d sites, %7d obs rows, %7d cohorts", name, length(ids), nrow(o), nrow(i))
}
write_set(site_tab$siteID, "national")
write_set(site_tab[set == "train", siteID], "train")
write_set(site_tab[set == "holdout", siteID], "holdout_blocked")
write_set(site_tab[pilot == TRUE, siteID], "pilot")

## ---- 11. PREPARATION.md ------------------------------------------------------------
summ <- obs_dt[, .(
  sites = uniqueN(siteID),
  ba_ha_mean = mean(ba, na.rm = TRUE) / PARAMS$patch_size_ha,
  growth_mean = weighted.mean(growth, growth_n, na.rm = TRUE),
  mort_annual = sum(n_died) / sum(n_at_risk),
  reg_ha_5yr = mean(reg))]
site_totals <- obs_dt[, .(ba = sum(ba), trees = sum(trees), reg = sum(reg)), by = siteID]
md <- c(
  sprintf("# Prepared data: %s", PARAMS$prep_version), "",
  sprintf("Written %s by IGN-FR_process/scripts/02_prepare_fiveyear_v1.R (FINN %s, data.table %s, terra %s).",
          format(Sys.time(), "%Y-%m-%d %H:%M %Z"), packageVersion("FINN"), packageVersion("data.table"), packageVersion("terra")), "",
  "Source: IGN, Inventaire forestier national, données brutes, campaigns 2005-2025 (export dated 2026-10-07,",
  "downloaded 2026-11-09, checksums in raw/SHA256SUMS.txt and raw/DOWNLOAD.md), Licence Ouverte / Open Licence Etalab 2.0.",
  "Climate: WorldClim 2.1 (1970-2000 normals, 30 arc-seconds) at the fuzzed plot centre (1 km grid cell).", "",
  "## What the tables mean", "",
  "- One site = one IGN plot (IDP), one patch of 0.0707 ha (the 15 m circle). Every tallied tree stands for W stems/ha",
  "  and is expanded to round(W x 0.0707) virtual trees (largest remainder within plot x species). `ba` (m2), `trees` and",
  "  `dbh` (cm) in obs_dt and `trees` in init_cohort_dt refer to that patch; divide ba and trees by 0.0707 for per-hectare values.",
  "- `year` 0 = first visit (campaign_t0 in site_table.csv), 5 = second visit; `period_length` = 5. Env rows years 1..5.",
  "- `growth` = annualised relative dbh increment of survivors, design-weighted; `growth_n` = tallied survivors behind it.",
  "- `n_at_risk` = tallied first-visit trees that survived or died naturally (cut / not-found censored); `mort` = annual",
  "  death probability from the design-weighted 5-year share; `n_died` = mort x n_at_risk (real-valued).",
  "- `reg` = new trees >= 7.5 cm (<= 15 cm) per hectare over the 5 years.",
  "- `state_masked` = plot with a cut or not-found tree: ba, dbh, trees at year 5 are NA there (FINN has no harvest).",
  "- env (natural units; FINN z-scales): tavg_c (bio1, degC), prec_mm (bio12, mm/yr), tmax_warm_c (bio5), tmin_cold_c (bio6),",
  "  prec_season (bio15, CV %), elev_m, x = longitude, y = latitude (WGS84, fuzzed plot centre).", "",
  "## Processing log", "", paste0("- ", log_lines), "",
  "## Summary (national)", "",
  sprintf("- %d sites; mean basal area at the first visit %.1f m2/ha (unmasked plots), mean annual relative growth %.4f,",
          summ$sites, mean(site_totals$ba, na.rm = TRUE) / PARAMS$patch_size_ha, summ$growth_mean),
  sprintf("  pooled annual mortality %.4f, mean recruitment %.0f stems/ha per 5 yr", summ$mort_annual, mean(site_totals$reg)), "",
  "## Sets", "",
  "- national: every site passing the filters", "- train: every site in the training blocks",
  "- holdout_blocked: every site in the held-out blocks (never for fitting)",
  sprintf("- pilot: %d training sites drawn with seed %d (smoke tests only)", PARAMS$pilot_n, PARAMS$pilot_seed), "",
  sprintf("Spatial blocking: %d km squares on EPSG:2154, %d of %d blocks held out (seed %d). Holdout site to nearest training site:",
          PARAMS$block_km, blocks[role == "holdout", .N], nrow(blocks), PARAMS$block_seed),
  sprintf("median %.1f km, minimum %.1f km, 10 %% below %.1f km, 90 %% below %.1f km.", sep$value[1], sep$value[2], sep$value[3], sep$value[4]), "",
  "Files per set: obs_dt.csv, env_dt.csv, init_cohort_dt.csv, species_dt.csv, siteID_dt.csv (siteID within the set,",
  "siteID_national, IDP, lon/lat, block, set, state_masked). Site attributes: site_table.csv; species list with ESPAR",
  "codes, French and Latin names: species.csv.", "",
  "## Species", "", "| species | name | ESPAR | French | stems at visit 1 | plots |", "|---|---|---|---|---|---|",
  species_out[, sprintf("| %d | %s | %s | %s | %d | %s |", species, species_name, espar, french_name, stems_t0, ifelse(is.na(plots_t0), "", plots_t0))], "",
  "## PARAMS", "", "```r", capture.output(str(PARAMS, give.attr = FALSE)), "```")
writeLines(md, file.path(OUT, "PREPARATION.md"))
writeLines(md, file.path("calibration", "PREPARATION.md"))          # the committed copy

## ---- 12. committed records: expected_preparation.csv, climate provenance ---------
## A summary of the final tables so a rerun (same vintage, same PARAMS) can be diffed.
per_set <- rbindlist(lapply(c("national", "train", "holdout_blocked", "pilot"), function(s) {
  o <- fread(file.path(OUT, s, "obs_dt.csv")); i <- fread(file.path(OUT, s, "init_cohort_dt.csv"))
  st <- o[, .(ba = sum(ba), trees = sum(trees), reg = sum(reg)), by = siteID]
  data.table(set = s, sites = uniqueN(o$siteID), obs_rows = nrow(o), cohorts = nrow(i), virtual_trees_t0 = sum(i$trees),
             sites_state_masked = o[state_masked == TRUE, uniqueN(siteID)],
             ba_ha_t5_mean = round(mean(st$ba, na.rm = TRUE) / PARAMS$patch_size_ha, 3),
             trees_ha_t5_mean = round(mean(st$trees, na.rm = TRUE) / PARAMS$patch_size_ha, 2),
             growth_annual_mean = round(weighted.mean(o$growth, o$growth_n, na.rm = TRUE), 6),
             n_at_risk = sum(o$n_at_risk), n_died = round(sum(o$n_died), 2),
             mort_annual_pooled = round(sum(o$n_died) / sum(o$n_at_risk), 6),
             reg_ha_5yr_mean = round(mean(st$reg), 4))
}))
per_species <- obs_dt[, .(set = "national", sites_present_t5 = sum(trees > 0, na.rm = TRUE), ba_sum = round(sum(ba, na.rm = TRUE), 3),
                          n_at_risk = sum(n_at_risk), n_died = round(sum(n_died), 2), reg_sum = round(sum(reg), 2)),
                      by = .(species, species_name)]
expected <- rbind(melt(per_set, id.vars = "set", variable.name = "quantity", value.name = "value")[, `:=`(species = NA_integer_, species_name = NA_character_)],
                  melt(per_species, id.vars = c("set", "species", "species_name"), variable.name = "quantity", value.name = "value"),
                  use.names = TRUE)
fwrite(expected[, .(set, species, species_name, quantity, value)], file.path("calibration", "expected_preparation.csv"))
prov_file <- file.path("calibration", "PROVENANCE.md")
prov <- if (file.exists(prov_file)) readLines(prov_file) else "# Provenance of every raw layer"
prov <- prov[!grepl("^## Climate", prov) & seq_along(prov) < (if (any(grepl("^## Climate", prov))) which(grepl("^## Climate", prov))[1] else length(prov) + 1L)]
wc <- rbindlist(lapply(names(PARAMS$env_rasters), function(nm) {
  f <- file.path(WORLDCLIM_DIR, PARAMS$env_rasters[[nm]])
  data.table(variable = nm, file = basename(f), sha256 = unname(tools::sha256sum(f)), size = file.size(f))
}))
prov <- c(prov, "## Climate: WorldClim 2.1 (02_prepare_fiveyear_v1.R)", "",
          "Source: https://www.worldclim.org/data/worldclim21.html (Fick & Hijmans 2017), 1970-2000 normals, 30 arc-seconds,",
          "files wc2.1_30s_bio.zip and wc2.1_30s_elev.zip; licence: WorldClim's own terms, free for academic and other non-commercial use, no redistribution without permission (the rasters are not redistributed, only values extracted at the plot locations are included); local copies pulled 2019-09-09 (file dates),",
          "values extracted at the fuzzed IGN plot centres (1 km cell) on the date of PREPARATION.md.", "",
          "| variable | file | SHA-256 | size (bytes) |", "|---|---|---|---|",
          wc[, sprintf("| %s | %s | %s | %d |", variable, file, sha256, size)], "")
writeLines(prov, prov_file)
message("done: ", OUT, "; calibration/PREPARATION.md, expected_preparation.csv, PROVENANCE.md")
