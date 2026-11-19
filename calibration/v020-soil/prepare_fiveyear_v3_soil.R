#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-IFN-France | calibration/v020-soil/prepare_fiveyear_v3_soil.R
##
## WHAT   Prepared-data version "fiveyear-v3-soil" = fiveyear-v2 plus two static
##        soil covariates per plot, soil_wet and soil_nutrient_poor, built from
##        the IGN soil pit (ECOLOGIE.csv) on the soil scales shared by the
##        zoo inventories (defined below). Nothing
##        else changes: every v2 file is a byte copy except env_dt.csv, which is
##        v2's table (same rows, same order, same columns) with the two soil
##        columns appended (the same value in the five env rows of a plot).
## WHY    The v020 model sees climate (WorldClim normals) and the management
##        trace only; nothing describes the soil, while French forest types are
##        sorted by water and nutrients. The general model will pool the
##        inventories later, so France uses the harmonised definitions, not its own.
##
## SOURCE IGN, Inventaire forestier national, raw data export 2005-2025 (DataIFN,
##        exported 2026-10-07), ECOLOGIE.csv: one soil pit per plot at the first
##        visit; codes in metadonnees.csv, documentation IGN_DB_doc_ecologie.pdf
##        (version 2.3, 14/10/2025). Licence Ouverte / Open Licence Etalab 2.0.
##        Matched to the plots by IDP (= siteName of fiveyear-v2); every v2 plot
##        has an ECOLOGIE row from its first-visit campaign (2010-2020), so the
##        humus and soil-type codes are the 2005-2021 lists (HUMUS, TSOL), not
##        the 2022 lists.
##
## READING OF THE DEPTH COLUMNS PGLEY / PPSEUDO / POX / PCALC (decided from the
##        IGN documentation): a number d is the depth class of the feature in dm
##        (d = 0..9: d*10 +-5 cm, 10 = >95 cm). The documentation says "En cas
##        d'absence, cette donnee n'est pas renseignee" and that the columns are
##        empty when no pit could be dug (OBSPEDO = 5). From 2017 the code list
##        1/10CM_X adds "X" = "phenomene absent". Hence: X = absent; empty =
##        absent when the pit was described, unknown when it was not (OBSPEDO =
##        5 and no depth recorded). Checked on the v2 plots: before 2017 no X
##        exists and empty cells are the absent features; from 2017 nearly every
##        described pit has X or a number and the empty cells are the OBSPEDO = 5
##        plots (the few exceptions are counted in PREPARATION.md, read as absent).
##
## soil_wet (0 no hydromorphy within 1 m .. 6 permanently wet, peat), wetness
##        crosswalk, France column:
##          depth score = the wettest of: PGLEY < 5 dm -> 5; PGLEY 5-9 dm -> 3;
##            PPSEUDO 0-2 dm -> 4; PPSEUDO 3-5 dm -> 3; PPSEUDO 6-9 dm -> 2;
##            POX 0-5 dm -> 2; POX 6-9 dm -> 1 (bins: "0-3 dm" = classes 0-2,
##            "3-6 dm" = 3-5, "within 6 dm" = 0-5, "within 5 dm" = 0-4)
##          type score  = TSOL 83 or HUMUS 80 (peat) -> 6; TSOL 82/85/89 (gley,
##            stagnogley, humic gley) or HUMUS 18/28/48/81 (hydromorphic humus)
##            -> 5; TSOL 81/87 (pseudogley) -> 4; TSOL 28/38/48/58 (hydromorphic
##            variants) -> 3; anything else -> 0
##          combination (harmonised rule: score from the shallowest hydromorphic
##            depth, TSOL / HUMUS alone where the depth columns are empty):
##            peat -> 6; else a type score of 5 -> max(depth score, 5); else the
##            depth score where any depth column records a feature; else the
##            type score (0 on a described pit without hydromorphy).
##        CHOICES the crosswalk leaves open: PGLEY deeper than 5 dm -> 3 ("groundwater
##        in the subsoil"); TSOL 28 (calcaire hydromorphe = calcarisol redoxique)
##        joins the hydromorphic variants 38/48/58 at 3; the soil-type criteria of
##        scores 3 and 4 do not override a recorded depth (the pseudogley soil
##        type TSOL 81 means mottles within 65 cm, which PPSEUDO places exactly),
##        those of scores 5-6 do (water at the surface, which mottle depths do
##        not measure). OBSHYDR is not used (kept in soil_plots.csv).
##
## soil_nutrient_poor (1 eutrophic .. 6 dystrophic), nutrient crosswalk:
##          humus score: 50 eumull, 55 / 45 carbonated mull -> 1; 40 mesomull -> 2;
##            31 oligomull, 30 dysmull -> 3; 22 hemimoder, 20 moder, 42 amphimull
##            -> 4; 21 dysmoder -> 5; 10 mor -> 6. Not in the table (CHOICES, by
##            the family the IGN key assigns them to): 15 carbonated mor -> 6 and
##            25 carbonated moder -> 4 (then the carbonate rule), 18 hydromor -> 6,
##            28 hydromoder -> 4, 48 hydromull -> 2, 81 anmoor -> 2 (a mull-type,
##            biologically active A horizon), 80 peat -> 4 (as the Dutch peat
##            row); 99 none -> missing.
##          carbonate rule (nutrient crosswalk: where a pH exists the pH decides and the
##            humus form moves the score by at most one class; France via PCALC):
##            active carbonate within 5 dm (PCALC 0-4), a carbonated humus (15,
##            25, 45, 55) or a calcareous soil type (TSOL 21, 23, 24, 25, 27, 28,
##            29, 31: rendzines, humocalcaric, calcaric brown soils) = pH > 6.5,
##            class 1 -> score = min(humus score, 2), and 1 where the humus is
##            missing.
##          podzol rule (harmonised rule): podzolic soils (TSOL 45, 51, 52, 53, 58,
##            59, 87) -> +1, at most 6.
##
## GAPS   A value is missing where the pit was not described and no soil type
##        was coded (soil_wet), or where the humus is 99 / empty without
##        carbonate evidence (soil_nutrient_poor). Every gap is filled with the
##        MEDIAN of the observed values of the same variable among the national
##        plots of the same sylvo-ecoregion (SER) and soil type (TSOL) when at
##        least min_group plots hold it, else of the same SER, else of the same
##        GRECO, else nationally. soil_plots.csv and <set>/soil_dt.csv record the
##        source of every value ("pit" or the fill level). No plot is dropped.
##
## SCALING env_dt.csv keeps the unscaled ordinal scores (as v2 keeps the climate
##        in natural units); the fit scripts z-scale them through FINN's
##        env_autoscale "auto" over their own training rows (as for tavg_c,
##        prec_mm), so the scores stay comparable across inventories.
##
## INPUT  <DATA_DIR>/prepared/fiveyear-v2/ (every set and top-level file; never
##        modified); <DATA_DIR>/raw/ign_data/unz/ECOLOGIE.csv
## OUTPUT <DATA_DIR>/prepared/fiveyear-v3-soil/
##          PREPARATION.md, soil_plots.csv (every national plot: raw codes,
##          scores, sources), soil_species_means.csv, soil_climate_correlation.csv,
##          v2_v3_check.csv; every v2 top-level file as a byte copy;
##          <set>/env_dt.csv  v2 + soil_wet, soil_nutrient_poor
##          <set>/soil_dt.csv siteID, siteID_national, IDP, the two scores, their sources
##          <set>/<other>     byte copy of v2
##        <DATA_DIR>/prepared/fiveyear-v3-soil.md5  manifest of every file
## RUN    Rscript calibration/v020-soil/prepare_fiveyear_v3_soil.R   (from the repository root, ~2 min)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
suppressPackageStartupMessages(library(data.table))

DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "frnfi")

PARAMS <- list(
  version      = "fiveyear-v3-soil",       # the new version (v2 + static soil)
  base_version = "fiveyear-v2",            # the base; never modified
  ecologie     = file.path("raw", "ign_data", "unz", "ECOLOGIE.csv"),   # the IGN soil pit, one row per plot at its first visit
  soil_vars    = c("soil_wet", "soil_nutrient_poor"),   # the two harmonised scores appended to env_dt.csv (the fit uses both)
  ## soil_wet (wetness crosswalk); depth classes in dm
  gley_shallow_dm   = 5L,                  # PGLEY below this class ("within 5 dm") -> 5, deeper -> 3
  pseudo_breaks_dm  = c(3L, 6L),           # PPSEUDO classes 0-2 -> 4, 3-5 -> 3, 6-9 -> 2
  ox_shallow_dm     = 6L,                  # POX below this class ("within 6 dm") -> 2, deeper -> 1
  tsol_peat    = "83",                     # Tourbe (Histosol) -> 6
  humus_peat   = "80",                     # Tourbe -> 6
  tsol_gley    = c("82", "85", "89"),      # gley, stagnogley, humic gley -> 5
  humus_hydro  = c("18", "28", "48", "81"),# hydromor, hydromoder, hydromull, anmoor -> 5
  tsol_pseudo  = c("81", "87"),            # pseudogley, podzolic pseudogley -> 4
  tsol_redox   = c("28", "38", "48", "58"),# calcareous / brown / leached / podzolic soils with redox features -> 3 (28 added to the crosswalk's 38/48/58)
  ## soil_nutrient_poor (nutrient crosswalk)
  humus_score  = c("50" = 1, "55" = 1, "45" = 1, "40" = 2, "31" = 3, "30" = 3, "22" = 4, "20" = 4, "42" = 4,
                   "21" = 5, "10" = 6, "15" = 6, "25" = 4, "18" = 6, "28" = 4, "48" = 2, "81" = 2, "80" = 4),   # Jabiol humus form -> class; 99 / empty missing
  carb_depth_dm = 5L,                      # active carbonate (PCALC) above this class ("within 5 dm") = carbonate present
  humus_carb   = c("15", "25", "45", "55"),# carbonated humus forms (A1 effervesces with HCl)
  tsol_calc    = c("21", "23", "24", "25", "27", "28", "29", "31"),   # rendzines, humocalcaric, calcaric brown soils: carbonate in the profile
  carb_cap     = 2,                        # carbonate = pH class 1; the humus moves it by at most one class -> at most 2
  tsol_podzol  = c("45", "51", "52", "53", "58", "59", "87"),         # podzolic soils: +1
  ## gaps
  min_group    = 10L,                      # plots with an observed value an SER x TSOL (or SER, GRECO) group needs to fill a gap
  climate_vars = c("tavg_c", "prec_mm", "tmax_warm_c", "tmin_cold_c", "prec_season", "elev_m", "gest")   # v2 env columns the soil scores are correlated with
)
BASE <- file.path(DATA_DIR, "prepared", PARAMS$base_version)
OUT  <- file.path(DATA_DIR, "prepared", PARAMS$version)
stopifnot(dir.exists(BASE), !dir.exists(OUT), file.exists(file.path(DATA_DIR, PARAMS$ecologie)))   # a version is never modified
sets <- list.dirs(BASE, full.names = FALSE, recursive = FALSE)
message("sets: ", paste(sets, collapse = ", "))

## ---- 1. the soil pit of every national plot -------------------------------------------
eco <- fread(file.path(DATA_DIR, PARAMS$ecologie), sep = ";", colClasses = "character", encoding = "UTF-8")   # empty cells stay ""
setnames(eco, sub("^﻿", "", names(eco)))
eco <- eco[, c("CAMPAGNE", "IDP", "OBSPEDO", "HUMUS", "TSOL", "PGLEY", "PPSEUDO", "POX", "OBSHYDR", "PCALC",
               "TEXT1", "TEXT2", "PROF1", "PROF2", "CAI40", "ROCHE"), with = FALSE]
stopifnot(!anyDuplicated(eco$IDP))
nat  <- fread(file.path(BASE, "national", "siteID_dt.csv"))[, .(siteID_national = siteID, IDP = as.character(IDP))]
st   <- fread(file.path(BASE, "site_table.csv"))[, .(IDP = as.character(IDP), campaign_t0, SER, GRECO)]
SP <- merge(merge(nat, st, by = "IDP"), eco, by = "IDP", all.x = TRUE)
stopifnot(nrow(SP) == nrow(nat), !anyNA(SP$CAMPAGNE), all(SP$CAMPAGNE == as.character(SP$campaign_t0)),
          all(as.integer(SP$CAMPAGNE) <= 2021L))   # the 2005-2021 HUMUS / TSOL code lists apply
setorder(SP, siteID_national)
dm <- function(x) suppressWarnings(as.integer(x))   # a depth class 0..9 (10 = >95 cm), NA for "X" (absent) and empty
SP[, `:=`(pgley = dm(PGLEY), ppseudo = dm(PPSEUDO), pox = dm(POX), pcalc = dm(PCALC))]
## a pit counts as described unless OBSPEDO = 5 and no depth column holds a number (before
## 2017 OBSPEDO 5 = "other case, soil not coded", and a few such plots still carry depths)
SP[, pit_described := OBSPEDO != "5" | !is.na(pgley) | !is.na(ppseudo) | !is.na(pox)]
stopifnot(SP[as.integer(CAMPAGNE) < 2017L, !any(c(PGLEY, PPSEUDO, POX) == "X")])   # X exists from 2017 only
DEPTH_READING <- SP[, .(plots = .N, pit_not_described = sum(!pit_described),
                        empty_depth_cell_on_described_pit_2017plus = sum(as.integer(CAMPAGNE) >= 2017L & pit_described & (PGLEY == "" | PPSEUDO == "" | POX == ""))),
                    by = .(era = fifelse(as.integer(CAMPAGNE) >= 2017L, "2017-2020 (X = absent)", "2010-2016 (empty = absent)"))]
print(DEPTH_READING)   # 2017+: an empty cell on a described pit is read as absent too (a handful of plots)
tsol_known <- function(x) !(x %in% c("", "99"))

## soil_wet
SP[, depth_score := pmax(fifelse(!is.na(pgley), fifelse(pgley < PARAMS$gley_shallow_dm, 5, 3), NA_real_),
                         fifelse(!is.na(ppseudo), fifelse(ppseudo < PARAMS$pseudo_breaks_dm[1], 4,
                                                          fifelse(ppseudo < PARAMS$pseudo_breaks_dm[2], 3, 2)), NA_real_),
                         fifelse(!is.na(pox), fifelse(pox < PARAMS$ox_shallow_dm, 2, 1), NA_real_), na.rm = TRUE)]
SP[, type_score := fcase(TSOL == PARAMS$tsol_peat | HUMUS == PARAMS$humus_peat, 6,
                         TSOL %in% PARAMS$tsol_gley | HUMUS %in% PARAMS$humus_hydro, 5,
                         TSOL %in% PARAMS$tsol_pseudo, 4,
                         TSOL %in% PARAMS$tsol_redox, 3,
                         default = 0)]
SP[, c("soil_wet", "soil_wet_source") := list(NA_real_, NA_character_)]
SP[type_score == 6, c("soil_wet", "soil_wet_source") := list(6, "pit: peat")]
SP[is.na(soil_wet) & type_score == 5, c("soil_wet", "soil_wet_source") := list(pmax(depth_score, 5, na.rm = TRUE), "pit: gley type / hydromorphic humus")]
SP[is.na(soil_wet) & !is.na(depth_score), c("soil_wet", "soil_wet_source") := list(depth_score, "pit: depth")]
SP[is.na(soil_wet) & (pit_described | tsol_known(TSOL)), c("soil_wet", "soil_wet_source") := list(type_score, "pit: soil type")]

## soil_nutrient_poor
SP[, humus_score := unname(PARAMS$humus_score[HUMUS])]
SP[, carbonate := (!is.na(pcalc) & pcalc < PARAMS$carb_depth_dm) | HUMUS %in% PARAMS$humus_carb | TSOL %in% PARAMS$tsol_calc]
SP[, podzol := TSOL %in% PARAMS$tsol_podzol]
SP[, soil_nutrient_poor := humus_score]
SP[carbonate == TRUE, soil_nutrient_poor := fifelse(is.na(humus_score), 1, pmin(humus_score, PARAMS$carb_cap))]
SP[podzol == TRUE & !is.na(soil_nutrient_poor), soil_nutrient_poor := pmin(soil_nutrient_poor + 1, 6)]
SP[, soil_nutrient_poor_source := fifelse(is.na(soil_nutrient_poor), NA_character_,
                                          fifelse(carbonate, "pit: humus + carbonate", fifelse(podzol, "pit: humus + podzol", "pit: humus")))]
stopifnot(SP[!is.na(soil_wet), all(soil_wet %in% 0:6)], SP[!is.na(soil_nutrient_poor), all(soil_nutrient_poor %in% 1:6)],
          all(SP[!HUMUS %in% c("", "99"), HUMUS] %in% names(PARAMS$humus_score)))   # every coded humus form has a class
OBSERVED <- SP[, .(variable = PARAMS$soil_vars, observed = c(sum(!is.na(soil_wet)), sum(!is.na(soil_nutrient_poor))), plots = .N)]
print(OBSERVED)

## ---- 2. gaps: median of the same SER x TSOL, SER, GRECO, nation ---------------------------
for (v in PARAMS$soil_vars) {
  src <- paste0(v, "_source")
  obs_v <- SP[!is.na(get(v))]
  g1 <- obs_v[tsol_known(TSOL), .(m = median(get(v)), n = .N), by = .(SER, TSOL)][n >= PARAMS$min_group]
  g2 <- obs_v[, .(m = median(get(v)), n = .N), by = .(SER)][n >= PARAMS$min_group]
  g3 <- obs_v[, .(m = median(get(v)), n = .N), by = .(GRECO)][n >= PARAMS$min_group]
  g0 <- median(obs_v[[v]])
  SP[g1, on = .(SER, TSOL), f1 := i.m][g2, on = .(SER), f2 := i.m][g3, on = .(GRECO), f3 := i.m]
  SP[is.na(get(v)) & tsol_known(TSOL) & !is.na(f1), c(v, src) := list(f1, "fill: SER x TSOL median")]
  SP[is.na(get(v)) & !is.na(f2), c(v, src) := list(f2, "fill: SER median")]
  SP[is.na(get(v)) & !is.na(f3), c(v, src) := list(f3, "fill: GRECO median")]
  SP[is.na(get(v)), c(v, src) := list(g0, "fill: national median")]
  SP[, c("f1", "f2", "f3") := NULL]
}
stopifnot(!anyNA(SP[, PARAMS$soil_vars, with = FALSE]))
COVER <- rbindlist(lapply(PARAMS$soil_vars, function(v) SP[, .N, by = c(source = paste0(v, "_source"))][, `:=`(variable = v, share = N / sum(N))]))
print(COVER[order(variable, -N)])
DIST <- rbindlist(lapply(PARAMS$soil_vars, function(v) SP[, .N, by = .(value = get(v))][, `:=`(variable = v, share = N / sum(N))]))[order(variable, value)]
print(dcast(DIST, value ~ variable, value.var = "N", fill = 0))

## ---- 3. plausibility: correlation with the climate, basal-area-weighted means per species --
env_nat <- fread(file.path(BASE, "national", "env_dt.csv"))[year == 1L, c("siteID", PARAMS$climate_vars), with = FALSE]
CC <- merge(env_nat, SP[, c("siteID_national", PARAMS$soil_vars), with = FALSE], by.x = "siteID", by.y = "siteID_national")
stopifnot(nrow(CC) == nrow(SP))
CORR <- CJ(soil = PARAMS$soil_vars, covariate = c(PARAMS$climate_vars, PARAMS$soil_vars), sorted = FALSE)[soil != covariate]
CORR[, `:=`(pearson = mapply(function(a, b) cor(CC[[a]], CC[[b]]), soil, covariate),
            spearman = mapply(function(a, b) cor(CC[[a]], CC[[b]], method = "spearman"), soil, covariate))]
print(CORR)
init_nat <- fread(file.path(BASE, "national", "init_cohort_dt.csv"))[, .(siteID, species_name, ba = trees * pi * (dbh / 200)^2)]
SPM <- merge(init_nat[, .(ba = sum(ba)), by = .(siteID, species_name)], CC[, c("siteID", PARAMS$soil_vars), with = FALSE], by = "siteID")[
  , .(plots = .N, ba_weighted_soil_wet = weighted.mean(soil_wet, ba), ba_weighted_soil_nutrient_poor = weighted.mean(soil_nutrient_poor, ba)), by = species_name][order(-ba_weighted_soil_wet)]
print(SPM)

## ---- 4. the new version: byte copies + env_dt with the soil columns --------------------
dir.create(OUT)
top_files <- setdiff(list.files(BASE), c(sets, "PREPARATION.md"))
stopifnot(all(file.copy(file.path(BASE, top_files), OUT, copy.date = TRUE)))
fwrite(SP[, c("siteID_national", "IDP", "CAMPAGNE", "SER", "GRECO", "OBSPEDO", "HUMUS", "TSOL", "PGLEY", "PPSEUDO", "POX", "OBSHYDR", "PCALC",
              "TEXT1", "TEXT2", "PROF1", "PROF2", "CAI40", "ROCHE", "depth_score", "type_score", "humus_score", "carbonate", "podzol",
              "soil_wet", "soil_wet_source", "soil_nutrient_poor", "soil_nutrient_poor_source"), with = FALSE],
       file.path(OUT, "soil_plots.csv"))
fwrite(SPM, file.path(OUT, "soil_species_means.csv"))
fwrite(CORR, file.path(OUT, "soil_climate_correlation.csv"))

CHECKS <- list()
for (s in sets) {
  dir.create(file.path(OUT, s))
  files <- list.files(file.path(BASE, s))
  copy <- setdiff(files, "env_dt.csv")
  stopifnot("env_dt.csv" %in% files, all(file.copy(file.path(BASE, s, copy), file.path(OUT, s), copy.date = TRUE)))
  sites <- fread(file.path(BASE, s, "siteID_dt.csv"))[, .(siteID, siteID_national, IDP = as.character(IDP))]
  site_soil <- merge(sites, SP[, c("IDP", "siteID_national", PARAMS$soil_vars, paste0(PARAMS$soil_vars, "_source")), with = FALSE],
                     by = c("IDP", "siteID_national"), all.x = TRUE)
  stopifnot(nrow(site_soil) == nrow(sites), !anyNA(site_soil[, PARAMS$soil_vars, with = FALSE]), !anyDuplicated(site_soil$siteID))
  fwrite(site_soil[order(siteID), c("siteID", "siteID_national", "IDP", PARAMS$soil_vars, paste0(PARAMS$soil_vars, "_source")), with = FALSE],
         file.path(OUT, s, "soil_dt.csv"))
  env <- fread(file.path(BASE, s, "env_dt.csv"))
  stopifnot(!any(PARAMS$soil_vars %in% names(env)))
  new <- copy(env)   # an update join keeps v2's row order
  new[site_soil, on = "siteID", (PARAMS$soil_vars) := mget(paste0("i.", PARAMS$soil_vars))]
  stopifnot(nrow(new) == nrow(env), !anyNA(new[, PARAMS$soil_vars, with = FALSE]),
            new[, all(vapply(.SD, uniqueN, integer(1)) == 1L), by = siteID, .SDcols = PARAMS$soil_vars][, all(V1)])
  fwrite(new, file.path(OUT, s, "env_dt.csv"))
  back <- fread(file.path(OUT, s, "env_dt.csv"))
  md5_old <- tools::md5sum(file.path(BASE, s, copy)); md5_new <- tools::md5sum(file.path(OUT, s, copy))
  CHECKS[[s]] <- data.table(
    set = s, env_rows = nrow(back), sites = uniqueN(back$siteID),
    v2_env_columns_and_rows_identical = isTRUE(all.equal(env, back[, names(env), with = FALSE], check.attributes = FALSE)),
    env_columns = paste(names(back), collapse = " "),
    other_files_byte_identical = all(unname(md5_old) == unname(md5_new)), n_other_files = length(copy),
    share_soil_wet_from_pit = mean(grepl("^pit", site_soil$soil_wet_source)),
    share_soil_nutrient_poor_from_pit = mean(grepl("^pit", site_soil$soil_nutrient_poor_source)))
}
CHECKS <- rbindlist(CHECKS)
top_ok <- all(unname(tools::md5sum(file.path(BASE, top_files))) == unname(tools::md5sum(file.path(OUT, top_files))))
print(CHECKS[, !"env_columns"])
stopifnot(all(CHECKS$v2_env_columns_and_rows_identical), all(CHECKS$other_files_byte_identical), top_ok)
fwrite(CHECKS, file.path(OUT, "v2_v3_check.csv"))

## ---- 5. PREPARATION.md and the manifest ----------------------------------------------------
md_table <- function(x) {
  x <- copy(x)
  for (cc in names(x)) if (is.numeric(x[[cc]])) x[, (cc) := if (all(get(cc) == round(get(cc)))) format(get(cc), big.mark = ",", trim = TRUE) else formatC(get(cc), digits = 3, format = "g")]
  c(paste("|", paste(names(x), collapse = " | "), "|"), paste("|", paste(rep("---", ncol(x)), collapse = " | "), "|"),
    apply(x, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}
commit <- suppressWarnings(tryCatch(system("git rev-parse --short HEAD", intern = TRUE, ignore.stderr = TRUE), error = function(e) "unknown"))
SUMM <- rbindlist(lapply(PARAMS$soil_vars, function(v) SP[, .(variable = v, mean = mean(get(v)), sd = sd(get(v)),
                                                              pit_mean = mean(get(v)[grepl("^pit", get(paste0(v, "_source")))]),
                                                              filled_mean = mean(get(v)[grepl("^fill", get(paste0(v, "_source")))]),
                                                              q05 = quantile(get(v), 0.05), q50 = median(get(v)), q95 = quantile(get(v), 0.95))]))
CORRw <- dcast(CORR[, .(soil, covariate, spearman = round(spearman, 2))], soil ~ covariate, value.var = "spearman")
prep <- c(
  "# Prepared data: fiveyear-v3-soil", "",
  sprintf("Written %s by `calibration/v020-soil/prepare_fiveyear_v3_soil.R` (FINN-model-IFN-France at commit `%s`).", format(Sys.time(), "%Y-%m-%d %H:%M %Z"), commit), "",
  "fiveyear-v3-soil = fiveyear-v2 plus two static soil covariates per plot, `soil_wet` and `soil_nutrient_poor`, on the soil scales shared by the zoo inventories. Nothing else changes.", "",
  "**Source.** IGN, Inventaire forestier national, raw data export 2005-2025 (DataIFN, exported 2026-10-07), `ECOLOGIE.csv` (soil pit at the first visit; codes `metadonnees.csv`; documentation `IGN_DB_doc_ecologie.pdf` v2.3). Licence Ouverte / Open Licence Etalab 2.0. Matched by IDP; every v2 plot has its first-visit row (campaigns 2010-2020, the 2005-2021 HUMUS / TSOL code lists).", "",
  "**Depth columns.** PGLEY / PPSEUDO / POX / PCALC: a number = depth class in dm (0..9); `X` (code list 1/10CM_X, from 2017) = feature absent; empty = absent when the pit was described (the IGN documentation: \"En cas d'absence, cette donnee n'est pas renseignee\"), unknown when no pit could be dug (OBSPEDO = 5). Before 2017 no `X` occurs; from 2017 the empty cells are the OBSPEDO = 5 plots, with the exceptions counted here (read as absent):", "",
  md_table(DEPTH_READING), "",
  "**soil_wet** (0 no hydromorphy within 1 m .. 6 peat). Depth score = the wettest of PGLEY < 5 dm -> 5, PGLEY 5-9 dm -> 3, PPSEUDO 0-2 dm -> 4, 3-5 dm -> 3, 6-9 dm -> 2, POX 0-5 dm -> 2, 6-9 dm -> 1. Type score = TSOL 83 / HUMUS 80 (peat) -> 6; TSOL 82/85/89 / HUMUS 18/28/48/81 -> 5; TSOL 81/87 -> 4; TSOL 28/38/48/58 -> 3; else 0. Combination: peat 6; a type score of 5 -> max(depth score, 5); else the depth score where a depth column records a feature; else the type score. Choices beyond the crosswalk: PGLEY deeper than 5 dm -> 3; TSOL 28 joins 38/48/58; the soil-type criteria of scores 3-4 never override a recorded depth, those of 5-6 do. OBSHYDR not used.", "",
  "**soil_nutrient_poor** (1 eutrophic .. 6 dystrophic). Humus class: 50/55/45 -> 1, 40 -> 2, 31/30 -> 3, 22/20/42 -> 4, 21 -> 5, 10 -> 6; beyond the crosswalk's table: 15 -> 6, 25 -> 4, 18 -> 6, 28 -> 4, 48 -> 2, 81 -> 2, 80 -> 4; 99 missing. Carbonate (PCALC within 5 dm, carbonated humus 15/25/45/55, or TSOL 21/23/24/25/27/28/29/31) = pH class 1, the humus moves it at most one class: min(humus class, 2), 1 without a humus. Podzolic soils (TSOL 45/51/52/53/58/59/87): +1, at most 6.", "",
  sprintf("**Gaps.** Filled with the median of the observed values of the same variable among the national plots of the same SER x TSOL (>= %d plots), else SER, else GRECO, else national; the source of every value is in `soil_plots.csv` and `<set>/soil_dt.csv`. No plot dropped.", PARAMS$min_group), "",
  md_table(dcast(COVER, variable ~ source, value.var = "share", fill = 0)), "",
  "**Distribution** over the national plots (count per value; fills by a median can be half classes):", "",
  md_table(dcast(DIST, value ~ variable, value.var = "N", fill = 0)), "",
  md_table(SUMM), "",
  "**Spearman correlation** with the v2 env covariates over the national plots (`soil_climate_correlation.csv` has Pearson too):", "",
  md_table(CORRw), "",
  "**Basal-area-weighted mean per species** at the first visit (national; `soil_species_means.csv`):", "",
  md_table(SPM[, .(species_name, plots, soil_wet = round(ba_weighted_soil_wet, 2), soil_nutrient_poor = round(ba_weighted_soil_nutrient_poor, 2))]), "",
  "**Scaling.** `env_dt.csv` keeps the unscaled ordinal scores (v2 keeps the climate in natural units too); the fit scripts z-scale them with FINN's `env_autoscale = \"auto\"` over their training rows.", "",
  "**Files.** `<set>/env_dt.csv` = v2's table with `soil_wet`, `soil_nutrient_poor` appended (v2 columns and row order identical: `v2_v3_check.csv`); new `<set>/soil_dt.csv`; every other file, in the sets and at the top level, a byte copy of v2 (md5 checked). Manifest: `../fiveyear-v3-soil.md5`.", "",
  md_table(CHECKS[, !"env_columns"]), "",
  "## PARAMS", "",
  sprintf("- `%s` = `%s`", names(PARAMS), vapply(PARAMS, function(x) paste(if (is.null(names(x))) x else paste(names(x), x, sep = ":"), collapse = ", "), character(1)))
)
writeLines(prep, file.path(OUT, "PREPARATION.md"))
all_files <- list.files(OUT, recursive = TRUE)
writeLines(sprintf("%s  ./%s", unname(tools::md5sum(file.path(OUT, all_files))), all_files), paste0(OUT, ".md5"))
message("done: ", OUT)
