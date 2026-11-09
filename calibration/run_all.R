#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-IFN-France | calibration/run_all.R
##
## WHAT   Runs the whole calibration chain in its fixed order, each step as its
##        own self-contained script (every scientific choice is inside the step
##        scripts, this file only sequences them):
##          1. 01_download_ign_raw.R      IGN raw data, checksums -> PROVENANCE.md
##          2. 02_prepare_fiveyear_v1.R   FINN tables fiveyear-v1, PREPARATION.md,
##                                        expected_preparation.csv
##          3. fit.R                      the production fit (hours; a cluster node)
##          4. validation.R               holdout report, metrics, figures
##        The smoke test fit_pilot.R is not part of the chain.
## RUN    Rscript calibration/run_all.R   (from the repository root)
##        Set `steps` below to rerun a part, e.g. only the validation.
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
steps <- c("01_download_ign_raw.R", "02_prepare_fiveyear_v1.R", "fit.R", "validation.R")
for (s in steps) {
  message(sprintf("==== %s  %s", format(Sys.time(), "%Y-%m-%d %H:%M"), s))
  status <- system2("Rscript", file.path("calibration", s))
  if (status != 0L) stop(sprintf("%s failed with status %d", s, status), call. = FALSE)
}
message("==== calibration chain complete ", format(Sys.time(), "%Y-%m-%d %H:%M"))
