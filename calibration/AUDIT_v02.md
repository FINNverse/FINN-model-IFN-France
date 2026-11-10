# Audit of the v0.1 fit and the IGN data for the v0.2 round

2026-11-10. Evidence: `calibration/v02/audit/*.csv|png` (written by
`calibration/v02/audit_v01.R` from the v0.1 fit `fits/fiveyear-v1-process`, its
holdout predictions and the raw IGN tables), `calibration/v02/species_subset_tradeoff.md`
(`species_subset_tradeoff.R`), and the v2 preparation log (`calibration/v02/PREPARATION_v2.md`).
Numbers are holdout (14,767 plots in unseen 100 km blocks) unless stated.

## A0. The v0.1 verdict had one wrong number

The v0.1 report said "mortality overpredicted 3.6-fold". It compared the observed
*annual* deaths (`n_died` = annual rate x trees at risk) with the predicted
*5-year* death probability. On one scale (`audit/mortality_level.csv`) the
v0.1 model predicts an annual natural mortality of 1.66 % against the site-level
target of 2.12 % (ratio 0.78; test 0.83). The card, NEWS and notes carry the
wrong statement and are corrected with v0.2. The AUC (0.66) is unaffected
(rank-based).

But the site-level target itself is inflated (A3): counted per tallied tree the
5-year natural death share is 5.2 % = 1.07 %/yr, and the fitted per-tree
probability is 1.55x that (A1b). So v0.1 mortality is *over*-predicted after
all, by 1.5x, not 3.6x, and for a reason that sits in the target, not in the
model.

## A1. Mortality: the fitted rate, the floored simulation, the target

**A1a. Level and gradients (fitted rate from the simulation, annual, not floored).**
By start diameter (`mortality_by_size.csv`): < 15 cm 0.0258 obs / 0.0237 pred,
15-25 0.0212 / 0.0162, 25-40 0.0181 / 0.0102, 40-60 0.0117 / 0.0053, > 60 cm
0.0074 / 0.0031: the fitted rate falls with size twice as fast as the data.
By species (`mortality_by_species.csv`): oaks, beech, hornbeam, fir within
10 %; Quercus pubescens 0.0024 vs 0.0081, Picea 0.019 vs 0.030, Castanea 0.027
vs 0.044, Fraxinus 0.039 vs 0.063 (ash dieback). By management trace: the
natural-death target does not vary much with GEST (1.5 / 1.7 / 1.1 %/yr).

**A1b. Fitted per-tree probability vs the 1 %/yr simulation floor.** FINN draws
the dying trees with the annual rate clamped to [0.01, 0.99]
(`finn-class.R`, `torch_clamp(m, 0.01, 0.99)`); the rate itself, the per-tree
loss and `predictTrees()` are not clamped. `FINN::predictTrees` of the v0.1
model on the v2 tallied trees of holdout1000 (`mortality_floor_per_tree.csv`,
`fig_mortality_floor.png`): observed 5-yr natural death share 3.6 % (weighted),
fitted 5.6 % (ratio 1.55); **62 % of the tallied trees have a fitted annual
rate below 1 %/yr** (Fagus 94 %, Pinus pinaster 89 %, Quercus pubescens and
Q. ilex 100 %, Abies 69 %, Q. petraea 63 %), so for most of the French stand the
simulation kills more trees than the fitted rate says. The floor binds on
exactly the species whose fitted rate is *over* the per-tree truth: the state
losses see the floored stand losing too many trees and have no way to lower it,
while the per-tree-equivalent site target (A3) pulls the rate up. Both effects
act on low-mortality species and sites and on every equilibrium run (the
bare-ground stands of v0.1 lost >= 1 %/yr of every cohort regardless of the
fitted rate). Treated as a model-side bias in v0.2: one pilot arm removes the
floor (`pilot_floor.R`, branch build `~/R/lib-finn-zoo-v02`); nothing is
compensated with weights.

**A1c. Removal, not natural mortality (v2 target).** Per tallied tree, 5-yr
removal 11.7 % (natural 5.2 %, cut 6.5 %) = 2.47 %/yr. By GEST
(`removal_by_gest.csv`): removal 1.82 / 2.39 / 2.71 %/yr for GEST 0 / 1 / 2,
cut 0.34 / 0.68 / 1.50, natural 1.46 / 1.66 / 1.14. By size
(`removal_by_size.csv`): natural death falls from 1.8 %/yr (< 15 cm) to
0.43 %/yr (> 60 cm), cutting rises from 0.9 to 1.3 %/yr, so total removal is
much flatter in size than natural death. v1 censored cuts and masked the
year-5 state on 17 % of the plots; the user's rule (keep total removal, management
as covariate; pre-emptive and sanitary felling is booked as a cut) is applied in
v2 with GEST as `gest`.

## A2. Recruits and ingrowth

New trees at the second visit, by dbh (`new_trees_by_dbh.csv`, v1 cohort):
7.5-10 cm 28,042 (27,345 recruits by the back-model), 10-12.5 cm 2,959 (1,490),
12.5-15 cm 634 (311), 15-22.5 cm 186 (34), >= 22.5 cm 6 (0). The circle thresholds
(22.5 / 37.5 cm) produce almost no "new" trees because the 9 and 15 m rings
only add trees that were already large; the real ingrowth confound is the
10-15 cm band: a tree of 12 cm cannot have been below 7.5 cm five years earlier
(species q90 of the small-tree increment 2.6 cm). v2 keeps a new tree as a
recruit only if dbh_t1 - q90_species < 7.5 cm (91 % of the new trees); the rest
are standing trees the first visit missed and are dropped. The recruitment
target changes little in level (v1 cap at 15 cm kept 99 % of them), but the
ingrowth-sized trees no longer enter the end stand.

## A3. Expansion, sample sizes and the annualisation trap

- Virtual-tree expansion: W x 0.0707 preserves basal area to -0.2 %; no tallied
  tree expands to zero (PB 6.25, MB 2.78, GB 1.0 copies). Sample sizes are
  tallied (v1 already; `n_at_risk` median 1-2 per species x site).
- **Annualising a site x species share is not the same as annualising the
  pooled share.** `mort = 1 - (1 - share_5yr)^(1/5)` per row: a row with one
  tree at risk that was cut has share 1 and annual rate 1. 30 % of the
  species x site rows have one tree at risk. Pooled over rows the v2 removal
  target is 6.0 %/yr, the natural-only v1 target 2.15 %/yr, against 2.47 and
  1.07 %/yr counted per tree (`PREPARATION_v2.md`, log). The binomial site loss
  is therefore biased upward by a factor 2-2.5, mostly through small cohorts,
  which is also where the v0.1 model's rate exceeds the per-tree truth. The
  per-tree loss (`fit(loss_aggregation = c(mortality = "tree"))`) scores each
  tallied tree's compounded 5-year survival and has no such bias: v2 ships
  `tree_obs_dt.csv` and the production fit uses it. The site-level `mort` stays
  in obs_dt for the state terms and for comparability, with this caveat.
- Intervals (trap 1): env rows 1..5, observation at year 5, `period_length` 5;
  `fit()` compares the observation at loop index 5 after exactly five steps
  (`year_sequence` = 5, window steps 1..5); simulateForest with env years 1..5
  returns years 1..5. Checked on the pilot model; no off-by-one.

## A4. Recruitment: the dispersion collapsed, the cap ran away

`regeneration_parameters.csv`, `regeneration_trajectory.csv`: the per-species
NB dispersion theta fell to 0.016-0.018 for 33 of 36 species (the raw
parameter drifted monotonically through all 15 snapshots, median 5 -> 57); the
cap K went to 1e4 stems/ha/yr for Salix caprea, 8e3 for Acer campestre, 5e3
for Populus tremula, and to 0.1-0.8 for twelve species. Predicted recruitment
per site sat in a flat band of 158-354 /ha/5 yr (observed mean 49), 154 of it
"other", 35 Abies alba, 18 + 16 the two maples, while beech, hornbeam, chestnut
and birch were *under*-predicted (ratios 0.3-0.7). This is BWI lesson 7: with
theta -> 0 the NB loss barely penalises over-prediction, and the loss kept
falling. v0.1 scored the NB per hectare (`recruit_obs_weight = 1`), where a
count of 1 tallied recruit is 88 /ha and no theta is meaningful. v0.2: theta
frozen at a method-of-moments estimate per species in tallied-count space
(`recruit_obs_weight` = 0.0113 ha, the 6 m circle), K fixed at the observed
99 % quantile per species, and no "other" class (A5).

## A5. "Other" and the species rule (agreed 2026-11-10)

The v1 "other" (126 codes, 5 % of stems) took 80 % of the excess recruitment
and dominated every bare-ground run. Rule: no lumped class; a species set S and
a plot rule on the START census only (species outside S < x % of basal area);
the END-census out-of-S share is an evaluation stratum, never a filter.
Trade-off (`species_subset_tradeoff.md`, cohort 58,528 plots with a GEST record,
16 holdout blocks):

| N | x % | plots kept | share | holdout plots | holdout blocks populated (>= 20) | tavg envelope | prec envelope | end share > 5 % | S stems lost | rarest species: plots / trees at risk / holdout blocks / tavg span |
|---|---|---|---|---|---|---|---|---|---|---|
| 20 | 2 | 35,265 | 0.60 | 8,866 | 14 | 0.91 | 1.01 | 1.2 % | 0.33 | Alnus glutinosa 747 / 4,300 / 3 / 0.37 |
| 25 | 2 | 42,551 | 0.73 | 10,868 | 14 | 0.94 | 0.99 | 0.8 % | 0.22 | Pinus halepensis 793 / 2,838 / 3 / 0.37 |
| **30** | **1** | 46,335 | 0.79 | 11,698 | 15 | 0.95 | 1.00 | 0.6 % | 0.17 | Larix decidua 499 / 2,116 / 3 / 0.37 |
| **30** | **2** | **47,013** | **0.80** | **11,881** | **15** | **0.95** | **1.00** | **0.6 %** | **0.15** | **Larix decidua 509 / 2,160 / 3 / 0.37** |
| 30 | 3 | 47,683 | 0.81 | 12,061 | 15 | 0.95 | 1.00 | 0.6 % | 0.14 | 519 / 2,217 / 3 / 0.37 |
| 30 | 5 | 48,754 | 0.83 | 12,304 | 15 | 0.95 | 1.00 | 1.2 % | 0.11 | 535 / 2,296 / 3 / 0.38 |
| 35 | 2 | 49,309 | 0.84 | 12,446 | 15 | 0.95 | 0.99 | 0.5 % | 0.12 | Arbutus unedo 286 / 1,613 / 3 / 0.37 |
| 40 | 2 | 51,410 | 0.88 | 13,052 | 15 | 0.97 | 1.00 | 0.4 % | 0.09 | Pinus mugo 277 / 1,458 / 3 / 0.37 |
| 50 | 2 | 54,016 | 0.92 | 13,675 | 15 | 0.98 | 1.00 | 0.3 % | 0.06 | Fraxinus angustifolia 140 / 643 / 2 / 0.37 |

(the "rarest" columns give the minimum over S of each count; at N = 30 the
minima are Larix decidua for plots and trees, Pinus halepensis for holdout
blocks (3) and climate span (0.40); full per-species table
`species_subset_identifiability_x2.csv`.)

**Choice: N = 30, x = 2 %.** Every species of S is on >= 500 plots, has
>= 2,000 trees at risk, is present in >= 3 holdout blocks and spans >= 40 % of
the cohort's temperature range (k = 500, m = 2,000, n = 3, span 0.4; these are
coverage choices, not measured identifiability limits, the planned simulation
study having been dropped). N = 35 adds Arbutus unedo (286 plots, Mediterranean
corner), Quercus rubra and cultivated Populus (plantations) for +2,300 plots;
N = 25 drops Pinus nigra, Robinia, Salix caprea, Sorbus torminalis and Crataegus
and loses 4,500 plots. x: from 1 % to 5 % the plots grow by 5 % while the plots
whose end-census out-of-S share exceeds 5 % double (0.6 -> 1.2 %); 2 % keeps the
competition of dropped trees negligible (776 out-of-S trees among 595,291 on
the kept plots) and leaves option (b) little to add. Population claim for the
card: valid for stands where species outside S hold < 2 % of the basal area at
the start of the interval = 80 % of the inventory's forest plots in the 2010-2020
cohort. Built as `fiveyear-v2`: 47,876 plots after the WorldClim join and the
diameter screen (train 35,783 in 66 blocks, holdout_blocked 12,093 in 16 blocks,
pilot2000 / holdout1000); per tallied tree the 5-year removal share is 11.8 %
(2.49 %/yr; natural 4.9 %); end-census strata < 1 % 46,991 / 1-5 % 991 / > 5 % 208
plots.
Out-of-S trees on kept plots are dropped (option a); the option-(b) tables
(one shade-only class, index 31, NA observations) are written for the pilot.

## A6. Convergence, loss competition, growth

- Loss flat from epoch ~50 (`loss_milestones.csv`: total baseline-scaled loss
  2.82 -> 2.21 at epoch 10 -> 2.27 at 50 -> 2.31 at 150) while every parameter
  block kept moving at a constant rate (`parameter_change_by_block.csv`: the
  mean absolute change between snapshots does not decrease; theta raw grows).
  v0.2 pilots test a cosine decay; the production fit uses it.
- Change metrics (`change_metrics.csv`): dBA observed +0.167 m2/patch, predicted
  +0.211, r 0.26; dN observed -0.2 stems/patch, predicted +12.1, r 0.16. State
  R2 (BA 0.94, N 0.78) is persistence plus a little (BA) or minus a lot (N).
- Growth (`growth_by_size.csv`): ratio pred/obs 0.98 (< 15 cm), 1.01, 0.93, 0.75
  (40-60 cm), 0.55 (> 60 cm): the size decline of growth is too steep, same
  sign as the mortality size gradient. Both point at the size terms being pulled
  by the state losses in the joint fit (BWI lesson 9); the per-tree losses and
  the removed recruitment explosion are the v0.2 answer; a staged recipe is the
  fallback if the pilots do not fix it.

## A7. What v0.2 changes (data) and what the pilots test (settings)

Data `fiveyear-v2` (`calibration/v02/02_prepare_fiveyear_v2.R`): S of 30 and
the 2 % start-census rule, no "other"; total removal + `gest`; recruits by the
back-modelled start diameter; per-tree table; pilot2000 / holdout1000; end-census
strata. Pilots (`calibration/v02/pilot_*.R`, one change each against the
baseline = v0.1 recipe on v2): cosine; thetafix (frozen theta, count-space NB,
fixed K); treemort (per-tree removal loss); floor (thetafix + treemort with
`mortality_floor = 1e-6`, poisson draw, branch build); shade (option b).
Evaluation `pilot_eval.R`: holdout1000 against persistence, change metrics,
per-tree removal AUC and level, recruitment ratio, bare ground, stratified by
the end-census out-of-S share. Production candidate `calibration/v02/fit.R`
(thetafix + treemort + cosine), to be confirmed against the pilot table.
