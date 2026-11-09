#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-zoo | IGN-FR_process/scripts/01_download_ign_raw.R
##
## WHAT   Downloads the raw data of the French national forest inventory
##        (IGN, "Inventaire forestier national", données brutes, all campaigns
##        2005-2025) and their documentation, verifies the SHA-256 checksums
##        recorded on the day the model was built, unpacks both archives and
##        writes raw/DOWNLOAD.md (URL, licence, access date, checksums).
## WHY    The zoo model must be reproducible from an official open source. The
##        checksums pin the export the model was built from: IGN replaces the
##        archive every year (and sometimes in between), so a later download
##        can differ; the script then stops and says so instead of silently
##        preparing different data.
## SOURCE https://inventaire-forestier.ign.fr/dataifn/ ("Télécharger les
##        données brutes"), Licence Ouverte / Open Licence Etalab 2.0. No
##        account, no terms click-through. Required citation:
##        "IGN - Inventaire forestier national français, Données brutes,
##        Campagnes annuelles 2005 et suivantes,
##        https://inventaire-forestier.ign.fr/dataIFN/, site consulté le 09/11/2026".
## INPUT  none (network)
## OUTPUT <DATA_DIR>/raw/ign_data/export_dataifn_2005_2025.zip + unz/ (CSV tables,
##        table documentation PDFs, metadonnees.csv, espar-cdref13.csv)
##        <DATA_DIR>/raw/ign_doc/export_dataifn_doc_2025.zip + unz/
##        <DATA_DIR>/raw/SHA256SUMS.txt, <DATA_DIR>/raw/DOWNLOAD.md
## RUN    Rscript IGN-FR_process/scripts/01_download_ign_raw.R   (from the repo root)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

## Where the data live, outside git: ../../FINN-zoo-data/frnfi relative to the
## repo root (the zoo convention). Change this one line for another layout.
DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "frnfi")

PARAMS <- list(
  ## the two archives IGN offers on the DataIFN start page (2026-11-09)
  files = list(
    data = list(
      url    = "https://inventaire-forestier.ign.fr/dataifn/data/export_dataifn_2005_2025.zip",
      dest   = file.path(DATA_DIR, "raw", "ign_data", "export_dataifn_2005_2025.zip"),
      ## SHA-256 of the archive downloaded 2026-11-09 (export dated 2026-10-07
      ## inside the archive); a different value means IGN changed the export
      sha256 = "be583e7eb01d43d4e74cc46097b8900092f1dcbfbe355c0ee3080001157c5b70"),
    doc = list(
      url    = "https://inventaire-forestier.ign.fr/dataifn/data/export_dataifn_doc_2025.zip",
      dest   = file.path(DATA_DIR, "raw", "ign_doc", "export_dataifn_doc_2025.zip"),
      sha256 = "8af9111099a2d59804bab8b696bb6deacc7aa1fced7725600a4e012b930f2480")
  ),
  licence  = "Licence Ouverte / Open Licence Etalab 2.0 (https://www.etalab.gouv.fr/wp-content/uploads/2017/04/ETALAB-Licence-Ouverte-v2.0.pdf)",
  citation = "IGN - Inventaire forestier national français, Données brutes, Campagnes annuelles 2005 et suivantes, https://inventaire-forestier.ign.fr/dataIFN/, site consulté le 09/11/2026",
  timeout_s = 1800        # a 71 MB archive over a slow line
)

options(timeout = PARAMS$timeout_s)
sha256 <- function(path) unname(tools::sha256sum(path))

lines <- c("# IGN raw data: download record", "",
           sprintf("Written %s by IGN-FR_process/scripts/01_download_ign_raw.R.", format(Sys.time(), "%Y-%m-%d %H:%M %Z")), "",
           paste("Licence:", PARAMS$licence), "",
           paste("Citation required by IGN:", PARAMS$citation), "",
           "| file | URL | SHA-256 | size (bytes) | downloaded |", "|---|---|---|---|---|")
sums <- character(0)
for (nm in names(PARAMS$files)) {
  f <- PARAMS$files[[nm]]
  dir.create(dirname(f$dest), showWarnings = FALSE, recursive = TRUE)
  fresh <- FALSE
  if (!file.exists(f$dest)) {
    message("downloading ", f$url)
    download.file(f$url, f$dest, mode = "wb", quiet = TRUE)
    fresh <- TRUE
  } else message("present: ", f$dest)
  got <- sha256(f$dest)
  if (!identical(got, f$sha256))
    stop(sprintf("checksum mismatch for %s:\n  expected %s\n  got      %s\nIGN has replaced the export; record the new checksum and re-check the preparation.",
                 basename(f$dest), f$sha256, got))
  unz_dir <- file.path(dirname(f$dest), "unz")
  if (!dir.exists(unz_dir) || !length(list.files(unz_dir, recursive = TRUE))) {
    dir.create(unz_dir, showWarnings = FALSE)
    unzip(f$dest, exdir = unz_dir)
  }
  sums <- c(sums, sprintf("%s  %s/%s", got, basename(dirname(f$dest)), basename(f$dest)))
  lines <- c(lines, sprintf("| %s | %s | %s | %d | %s |", basename(f$dest), f$url, got,
                            file.size(f$dest), if (fresh) format(Sys.Date()) else "2026-11-09 (kept)"))
}
writeLines(sums, file.path(DATA_DIR, "raw", "SHA256SUMS.txt"))
lines <- c(lines, "",
           "Unpacked tables (ign_data/unz): ARBRE.csv (trees, one row per tree x campaign), PLACETTE.csv (plots),",
           "BOIS_MORT.csv, COUVERT.csv, ECOLOGIE.csv, FLORE.csv, HABITAT.csv; metadonnees.csv (codes and labels),",
           "espar-cdref13.csv (ESPAR species code -> TAXREF v13 name), IGN_DB_doc_*.pdf (table documentation),",
           "IGN_DB_documentation_generale.pdf, precision-donnees_2020-2024.csv. CSV separator ';', decimal '.'.", "",
           "Plot coordinates XL, YL are Lambert-93 (EPSG:2154) centres of the 1 km sampling-grid cell nearest",
           "to the plot (true position up to ~700 m away), as released by IGN.")
writeLines(lines, file.path(DATA_DIR, "raw", "DOWNLOAD.md"))
## the committed provenance record of the calibration (calibration/PROVENANCE.md):
## the IGN layers from this script; 02_prepare_fiveyear_v1.R appends the climate layers
prov <- c("# Provenance of every raw layer", "",
          "Pinned vintage: IGN export dated 2026-10-07 (campaigns 2005-2025), downloaded 2026-11-09. The IGN archive is a",
          "live download replaced every year; a later archive with another checksum is another vintage and needs a new",
          "prepared-data version.", "",
          "## IGN raw data (01_download_ign_raw.R)", "",
          paste("Licence:", PARAMS$licence), "", paste("Citation required by IGN:", PARAMS$citation), "",
          lines[grep("^\\| ", lines)], "")
writeLines(prov, file.path("calibration", "PROVENANCE.md"))
message("done: ", file.path(DATA_DIR, "raw", "DOWNLOAD.md"), " and calibration/PROVENANCE.md")
