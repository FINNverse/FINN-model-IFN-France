#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-IFN-France | calibration/v02/species_subset_tradeoff.R
##
## WHAT   The trade-off behind the v0.2 species rule (user decision 2026-11-10:
##        no "other" class; keep only plots where the species outside the
##        modelled set S hold less than x % of the basal area at BOTH visits).
##        For S = the N most frequent species (by living tallied stems at the
##        first visit, N = 15 ... 50) and x = 1, 2, 3, 5 %, the table gives:
##        plots retained (of the v1-filtered cohort with a GEST record), holdout
##        blocks still populated (>= 20 plots), the climate envelope retained
##        (share of the 2.5-97.5 % range of tavg_c and prec_mm of all plots that
##        the retained plots still span), the stems of S lost with the dropped
##        plots, the END-census out-of-S share on the retained plots (the
##        evaluation stratum; the rule itself uses the START census only) and
##        the identifiability columns of the rarest species in S: plots, trees
##        at risk, holdout blocks with presence, climate span.
## INPUT  <DATA_DIR>/raw/ign_data/unz/{PLACETTE.csv, ARBRE.csv, espar-cdref13.csv},
##        <DATA_DIR>/prepared/fiveyear-v1/site_table.csv (plots, blocks, climate; the
##        v1 plot filters minus GEST are reproduced from it)
## OUTPUT calibration/v02/species_subset_tradeoff.csv, species_subset_tradeoff.md,
##        species_subset_candidates.csv (which species enter S for each N)
## RUN    Rscript calibration/v02/species_subset_tradeoff.R   (repository root, ~2 min)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
suppressPackageStartupMessages(library(data.table))
DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "frnfi")
RAW <- file.path(DATA_DIR, "raw", "ign_data", "unz")
OUT <- file.path("calibration", "v02")
PARAMS <- list(N_candidates = c(15L, 20L, 25L, 30L, 35L, 40L, 50L), x_candidates = c(1, 2, 3, 5),
               min_block_plots = 20L, env_quantiles = c(0.025, 0.975), gest_keep = c("0", "1", "2"))

site <- fread(file.path(DATA_DIR, "prepared", "fiveyear-v1", "site_table.csv"))   # the v1 plot cohort (filters applied), blocks, climate
site <- site[GEST %in% PARAMS$gest_keep]
pl_ids <- site$IDP
ar <- fread(file.path(RAW, "ARBRE.csv"), sep = ";", encoding = "UTF-8", colClasses = list(character = c("ESPAR", "VEGET", "VEGET5")), select = c("CAMPAGNE", "IDP", "A", "ESPAR", "VEGET", "VEGET5", "C13", "W"))
ar <- ar[IDP %in% pl_ids]; ar <- merge(ar, site[, .(IDP, campaign_t0, campaign_t1)], by = "IDP"); ar[, dbh := C13 * 100 / pi]
nomW <- function(d) c(1e4 / (pi * 6^2), 1e4 / (pi * 9^2), 1e4 / (pi * 15^2))[findInterval(d, c(22.5, 37.5)) + 1L]
t0 <- ar[CAMPAGNE == campaign_t0 & VEGET == "0" & is.finite(dbh) & dbh >= 7.5]; t0[!is.finite(W) | W <= 0, W := nomW(dbh)]
t1 <- ar[CAMPAGNE == campaign_t1]
## living trees at the second visit: survivors (remeasured) and new trees; W from the first visit where the second has none
## a second-visit row carries neither the species nor (2021+) the weight of a remeasured tree: take both from its first-visit row
w0 <- t0[, .(IDP, A, W0 = W, ESPAR0 = ESPAR)]
alive1 <- t1[(VEGET5 == "0" | (VEGET5 == "" & VEGET == "0")) & is.finite(dbh) & dbh >= 7.5]
alive1 <- merge(alive1, w0, by = c("IDP", "A"), all.x = TRUE)
alive1[ESPAR == "" & !is.na(ESPAR0), ESPAR := ESPAR0]
alive1[!is.finite(W) | W <= 0, W := W0]; alive1[!is.finite(W) | W <= 0, W := nomW(dbh)]
ba0 <- t0[, .(ba = sum(W * (dbh / 200)^2 * pi)), by = .(IDP, ESPAR)]; ba1 <- alive1[, .(ba = sum(W * (dbh / 200)^2 * pi)), by = .(IDP, ESPAR)]
sp <- fread(file.path(RAW, "espar-cdref13.csv"), sep = ";", encoding = "UTF-8", header = FALSE, skip = 1L, col.names = c("espar", "lib_espar", "cd_ref", "lib_cdref"), colClasses = "character")
sp[, espar := sub("^// ", "", espar)]; sp[, ESPAR := fifelse(nchar(espar) == 1L, paste0("0", espar), espar)]
sp[, latin := trimws(sub("\\s+(subsp\\.|var\\.|x)\\s.*$", "", lib_cdref))]
sp[, latin := vapply(strsplit(latin, " "), function(w) paste(w[seq_len(min(2L, length(w)))], collapse = " "), character(1))]
rank <- t0[, .(stems = .N, plots = uniqueN(IDP)), by = ESPAR][order(-stems)]
rank <- merge(rank, sp[, .(ESPAR, latin)], by = "ESPAR", all.x = TRUE)[order(-stems)]
rank[is.na(latin), latin := paste0("ESPAR_", ESPAR)]
## the same binomial under two codes (Pinus nigra) counts as one species; rank by binomial
rank_sp <- rank[!startsWith(latin, "ESPAR_"), .(stems = sum(stems), plots = max(plots), codes = paste(ESPAR, collapse = "+")), by = latin][order(-stems)]
rank_sp[, N_rank := .I]
fwrite(rank_sp, file.path(OUT, "species_subset_candidates.csv"))
code2sp <- rank[, .(ESPAR, latin)]
ba0 <- merge(ba0, code2sp, by = "ESPAR", all.x = TRUE); ba1 <- merge(ba1, code2sp, by = "ESPAR", all.x = TRUE)
env_rng <- function(x) quantile(x, PARAMS$env_quantiles)
rng_all <- list(t = env_rng(site$tavg_c), p = env_rng(site$prec_mm))
## the plot rule is applied to the START census only (agreed rule: never condition on the
## end census); the end-census out-of-S share is reported as the evaluation stratum
t0 <- merge(t0, code2sp, by = "ESPAR", all.x = TRUE)
## per-species identifiability counts on the retained plots: plots present, trees at
## risk (living at visit 1), holdout blocks with presence, climate span (2.5-97.5 % of
## the species' plots relative to the cohort's range)
fate <- merge(t0[, .(IDP, A)], t1[, .(IDP, A, VEGET5)], by = c("IDP", "A"), all.x = TRUE)
t0[, at_risk := fate$VEGET5 %in% c("0", "M", "5", "C", "2", "A", "T", "1", "6", "7")]
res <- rbindlist(lapply(PARAMS$N_candidates, function(N) rbindlist(lapply(PARAMS$x_candidates, function(x) {
  S <- rank_sp$latin[seq_len(N)]
  s0 <- ba0[, .(share_out = sum(ba[!latin %in% S]) / sum(ba)), by = IDP]
  keep <- s0[share_out < x / 100, IDP]
  k <- site[IDP %in% keep]
  s1 <- ba1[IDP %in% keep, .(share_out = sum(ba[!latin %in% S]) / sum(ba)), by = IDP]
  rng_k <- list(t = env_rng(k$tavg_c), p = env_rng(k$prec_mm))
  tk <- merge(t0[IDP %in% keep & latin %in% S], site[, .(IDP, block, block_role, tavg_c)], by = "IDP")
  ident <- tk[, .(plots = uniqueN(IDP), trees_at_risk = sum(at_risk), holdout_blocks = uniqueN(block[block_role == "holdout"]),
                  tavg_span = diff(env_rng(tavg_c)) / diff(rng_all$t)), by = latin]
  data.table(N = N, x_pct = x, plots = length(keep), plots_share = length(keep) / nrow(site),
             holdout_plots = k[block_role == "holdout", .N], holdout_blocks_populated = k[block_role == "holdout", .N, by = block][N >= PARAMS$min_block_plots, .N],
             holdout_blocks_total = site[block_role == "holdout", uniqueN(block)],
             tavg_envelope_share = diff(rng_k$t) / diff(rng_all$t), prec_envelope_share = diff(rng_k$p) / diff(rng_all$p),
             end_share_out_median = median(s1$share_out), end_share_gt5pct = mean(s1$share_out > 0.05),
             stems_S_lost_share = 1 - t0[IDP %in% keep & latin %in% S, .N] / t0[latin %in% S, .N],
             species_min_plots = min(ident$plots), species_min_trees_at_risk = min(ident$trees_at_risk),
             species_min_holdout_blocks = min(ident$holdout_blocks), species_min_tavg_span = min(ident$tavg_span),
             rarest_species = ident[which.min(plots), latin])
}))))
## the per-species identifiability table for the chosen candidates is written for every N (x = 2 %)
ident_all <- rbindlist(lapply(PARAMS$N_candidates, function(N) {
  S <- rank_sp$latin[seq_len(N)]; keep <- ba0[, .(share_out = sum(ba[!latin %in% S]) / sum(ba)), by = IDP][share_out < 0.02, IDP]
  tk <- merge(t0[IDP %in% keep & latin %in% S], site[, .(IDP, block, block_role, tavg_c)], by = "IDP")
  cbind(N = N, tk[, .(plots = uniqueN(IDP), trees_at_risk = sum(at_risk), holdout_blocks = uniqueN(block[block_role == "holdout"]), tavg_span = diff(env_rng(tavg_c)) / diff(rng_all$t)), by = latin])
}))
fwrite(ident_all, file.path(OUT, "species_subset_identifiability_x2.csv"))
fwrite(res, file.path(OUT, "species_subset_tradeoff.csv"))
show <- copy(res); num <- names(show)[vapply(show, is.numeric, logical(1))]; show[, (num) := lapply(.SD, signif, 3), .SDcols = num]
md <- c("# Species set S and plot purity x: the trade-off", "",
        sprintf("Written %s by calibration/v02/species_subset_tradeoff.R. Cohort: the v1 plot filters plus a GEST record, %d plots, %d holdout blocks. S = the N most frequent species (living tallied stems at the first visit; Pinus nigra's two codes as one species). A plot is kept when species outside S hold < x %% of its basal area at BOTH visits (W-weighted). `holdout_blocks_populated` = held-out 100 km blocks with >= %d retained plots; `*_envelope_share` = width of the retained 2.5-97.5 %% climate range relative to all plots; `stems_S_lost_share` = stems of S on dropped plots; `species_min_plots` = plots of the rarest species of S among the retained plots.",
                format(Sys.time(), "%Y-%m-%d %H:%M"), nrow(site), site[block_role == "holdout", uniqueN(block)], PARAMS$min_block_plots), "",
        paste0("| ", paste(names(show), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(show)), collapse = "|"), "|"),
        apply(show, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")), "",
        "## Species by rank", "", "| rank | species | stems | plots | ESPAR |", "|---|---|---|---|---|",
        rank_sp[1:50, sprintf("| %d | %s | %d | %d | %s |", N_rank, latin, stems, plots, codes)])
writeLines(md, file.path(OUT, "species_subset_tradeoff.md"))
print(res[, .(N, x_pct, plots, plots_share = round(plots_share, 3), holdout_plots, holdout_blocks_populated, tavg_env = round(tavg_envelope_share, 3), prec_env = round(prec_envelope_share, 3),
              end_gt5 = round(end_share_gt5pct, 3), S_lost = round(stems_S_lost_share, 3), min_plots = species_min_plots, min_risk = species_min_trees_at_risk, min_hblocks = species_min_holdout_blocks, min_span = round(species_min_tavg_span, 2), rarest = rarest_species)])
