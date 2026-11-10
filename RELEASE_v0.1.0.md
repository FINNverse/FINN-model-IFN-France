# Release checklist: FINNmodelIFNFrance v0.1.0

Prepared 2026-11-10. Nothing here has been executed; tag, push and release only
after the go. All commands from the repository root.

## Commit to tag: `399541e` on branch `release-v0.1.0` (NOT on `main`)

`399541e` ("WorldClim licence stated correctly ...") is the tip of the release
branch: the release commit `cbbd214`, the checklist commit `a23a2ac` and the
licence correction (README, NEWS, PROVENANCE and the script writing it,
release notes: WorldClim is not CC BY-SA 4.0, it is free for academic and other
non-commercial use, rasters not redistributed). The same correction is on
`main` as `3df3c3d`.

`main` cannot be tagged for v0.1.0: commit `7f6564c` ("v0.2 audit, species rule
trade-off and fiveyear-v2 preparation, six pilots, production candidate")
rewrote `R/definition.R` for the v0.2 structure (30 species, covariates
`tavg_c + prec_mm + gest`, fixed cap as a buffer) while `inst/model/` still
describes v0.1 (36 species, two covariates). Under that `build_model()` the
v0.1 parameter file fails: `parameters.pt lacks env_scaling.gest.center`
(checked). The v0.2 scripts in `calibration/v02/fit.R` source
`R/definition.R`, so reverting it on `main` would break the v0.2 work.

Therefore the release lives on `release-v0.1.0`, branched from `bd30ddd` (the
last v0.1 commit: production fit, release assets, NEWS) plus the release commit
`cbbd214` ("FINNmodelIFNFrance 0.1.0: version, FINN pinned at 0848b60, ...")
and the two commits after it (checklist, licence correction).
On that branch `R/`, `inst/model/`, `DESCRIPTION` 0.1.0, vignette and
`calibration/` are the v0.1 package and pipeline; `calibration/`, `hpc/`,
`logs/` and this file are excluded from the build by `.Rbuildignore`. `main`
keeps the v0.2 work untouched. When v0.2 is released from `main`, the v0.1
package stays reachable through the tag.

## Assets (outside git), sha256 verified 2026-11-10 with `shasum -a 256`

| asset | local path (relative to this repo) | sha256 | size |
|---|---|---|---|
| parameters.pt | ../../FINN-zoo-data/frnfi/fits/fiveyear-v1-process/parameters.pt | f278626cc0a3d5625a8e6437d98d974d4150311bc105e7e042717391351a0218 | 21,732 B |
| model.pt | ../../FINN-zoo-data/frnfi/fits/fiveyear-v1-process/model.pt | e3c52b4d830118ab3fce3511c2f68cf931adc98420f7cdca665af6c93991db05 | 100 MB |

The same values are in `inst/model/weights.json`; `load_model()` refuses a
file whose sha256 differs. (`fits/pilot-smoke/` holds a smoke test, not a
release asset.)

## Checks done (2026-11-10, one thread, FINN 0.2.0 = 0848b60, torch 0.15.1, R 4.5.0, on `release-v0.1.0`)

- `R CMD build` with the vignette (chunks `eval = FALSE`: they need the release
  download) and `R CMD INSTALL` into a temporary library: OK.
- `load_model(path = <parameters.pt>)`: sha256 match, model rebuilt (24
  parameter tensors), 3-year bare-ground simulation with the pinned batchsize
  runs. `model.pt` sha256 match.
- Vignette renders.

## Steps

```sh
# 0. sanity: the sha256 of the files to upload equal weights.json
shasum -a 256 ../../FINN-zoo-data/frnfi/fits/fiveyear-v1-process/parameters.pt ../../FINN-zoo-data/frnfi/fits/fiveyear-v1-process/model.pt

# 1. tag the release commit on the release branch
git tag -a v0.1.0 399541e -m "FINNmodelIFNFrance 0.1.0: first fitted release (variant process), FINN 0.2.0 commit 0848b60"

# 2. push the release branch, main (unchanged by the release) and the tag
git push origin release-v0.1.0
git push origin main
git push origin v0.1.0

# 3. create the GitHub release with both assets (the tag already exists, so
#    gh uses it; --target is not needed)
gh release create v0.1.0 --repo FINNverse/FINN-model-IFN-France \
  --title "FINNmodelIFNFrance 0.1.0" \
  --notes-file RELEASE_NOTES_v0.1.0.md \
  ../../FINN-zoo-data/frnfi/fits/fiveyear-v1-process/parameters.pt \
  ../../FINN-zoo-data/frnfi/fits/fiveyear-v1-process/model.pt

# 4. verify the uploaded assets byte for byte
tmp=$(mktemp -d) && gh release download v0.1.0 --repo FINNverse/FINN-model-IFN-France -D "$tmp" && shasum -a 256 "$tmp"/*.pt

# 5. (user's go) make the repository public
gh repo edit FINNverse/FINN-model-IFN-France --visibility public --accept-visibility-change-consequences
```

Note: `RELEASE_NOTES_v0.1.0.md` and this file live on `release-v0.1.0` (not
on `main`); run the `gh release create` from that branch's checkout, or point
`--notes-file` at the file there.

## After the release

```r
remotes::install_github("FINNverse/FINN-model-IFN-France@v0.1.0")
library(FINNmodelIFNFrance); torch::torch_set_num_threads(1L)
m <- load_model()             # downloads from the release URL in weights.json and checks the sha256
```

Then the FINN-zoo card (`IGN-FR_process/model-card.qmd`) cites tag `v0.1.0`
and these sha256 values. The README and install line on `main` still say
"in preparation"; they get updated with the v0.2 release (or cherry-pick the
README hunk of `cbbd214` onto `main` if wanted earlier).

## Release notes

The text of `RELEASE_NOTES_v0.1.0.md` (used by `--notes-file` above).
