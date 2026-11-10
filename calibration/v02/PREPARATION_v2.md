# Prepared data: fiveyear-v2

Written 2026-11-10 07:14 CET by calibration/v02/02_prepare_fiveyear_v2.R (FINN 0.2.0). Raw layers as v1 (calibration/PROVENANCE.md).

## Changes against fiveyear-v1

1. No 'other' class. S = the 30 most frequent species; a plot is kept when species outside S hold < 2 % of its basal area at the START census (48190 of 59960 plots). Out-of-S trees on kept plots are dropped from the model tables (option a) and written as one shade-only class (index 31, NA observations) in obs_dt_shade.csv / init_cohort_shade_dt.csv (option b). The END-census out-of-S share is recorded per plot (`end_share_out`, `end_stratum` < 1 % / 1-5 % / > 5 %) for stratified evaluation; it is not a filter.
2. Mortality target = TOTAL REMOVAL (dead, fallen, broken or cut); not-found censored; no state masked. Covariate `gest` (IGN GEST at visit 1). Plots without GEST dropped. Natural-only diagnostics: `mort_natural`, `n_died_natural`, `cut_share_5yr`. The site x species annualisation inflates complete removals (see the log); the per-tree table carries the unbiased fates.
3. Recruits = new trees whose back-modelled start diameter (dbh_t1 - species q90 of the small-tree 5-yr increment) is below 7.5 cm; other new trees dropped.
4. tree_obs_dt.csv per set: tallied trees at risk (siteID, year 5, patchID 1, species, dbh at visit 1, trees = virtual copies, annual relative growth of survivors, died = removal 0/1, died_natural, cut, period_length 5).
5. Sets pilot2000 (seed 1, training blocks) and holdout1000 (seed 11, holdout blocks).

Everything else as v1: one 0.0707 ha patch per plot, W x 0.0707 virtual trees, tallied n_at_risk / growth_n, annualised rates, reg per ha per 5 yr, env years 1..5, WorldClim 2.1 normals in natural units, 100 km blocks (seed 21).

## Processing log

- plots first-visited 2010-2020 with a second visit: 66021
-   land cover 1/2/3 at visit 1: 63907
-   no planting / sowing in the interval: 63285
-   no clearing in the interval: 63170
-   management trace GEST recorded (0/1/2): 62919  [v2]
-   living visit-1 trees >= 7.5 cm: 747100 (436 below the threshold dropped)
-   fates: survivor 659139, cut 48254, dead 39165, lost 542
-   new trees at visit 2 (alive, >= 7.5 cm): 36616; recruits by the back-modelled start diameter (pooled q90 2.61 cm): 33321; 3295 dropped as missed at visit 1  [v2]
- species: S = 30 species (Fagus sylvatica ... Ulmus minor); plot rule out-of-S basal area < 2 % at the START census: 48190 of 59960 plots kept  [v2]
-   end-census out-of-S share on the kept plots: <1% 46991, 1-5% 991, >5% 208
-   out-of-S trees on the kept plots: 776 of 595291 living visit-1 trees dropped from the model tables (kept as the shade class in the option-(b) tables); out-of-S recruits dropped
- expansion: basal area 24.8 (design) vs 24.8 m2/ha (virtual), -0.18 % off
- env: WorldClim for 48180 of 48190 plots
- makeObsData / resolveSiteIDs: 47876 sites, 30 species
- removal target: pooled site-level annual removal 0.0599 (natural 0.0181); per tallied tree: 5-yr removal share 0.1183 (annual 0.0249), natural 0.0492
- blocks: 16 of 82 held out; 35783 training, 12093 holdout sites
- holdout separation: median 21.2 km, min 1.0 km, 10 % below 5.1 km
- set national          47876 sites, 1436280 obs rows,  590024 cohorts (+766 shade),  589714 tallied trees
- set train             35783 sites, 1073490 obs rows,  442939 cohorts (+565 shade),  442709 tallied trees
- set holdout_blocked   12093 sites,  362790 obs rows,  147085 cohorts (+201 shade),  147005 tallied trees
- set pilot2000          2000 sites,   60000 obs rows,   24367 cohorts (+26 shade),   24357 tallied trees
- set holdout1000        1000 sites,   30000 obs rows,   12113 cohorts (+22 shade),   12099 tallied trees

## Summary (national)

- 47876 sites; mean basal area at year 5 25.3 m2/ha, 679 stems/ha; per-tree 5-yr removal share 0.1185; mean recruitment 39.6 stems/ha per 5 yr

## Species (S) with identifiability counts on the kept plots

| species | name | ESPAR | French | plots | stems at visit 1 | trees at risk | holdout blocks |
|---|---|---|---|---|---|---|---|
| 1 | Abies alba | 61 | Sapin pectin<U+00E9> | 4495 | 33590 | 33577 | 13 |
| 2 | Acer campestre | 21C | <U+00C9>rable champ<U+00EA>tre | 3276 | 8757 | 8753 | 15 |
| 3 | Acer pseudoplatanus | 15S | <U+00C9>rable sycomore | 2367 | 6754 | 6750 | 14 |
| 4 | Alnus glutinosa | 13G | Aulne glutineux | 1112 | 5740 | 5734 | 13 |
| 5 | Betula pendula | 12V | Bouleau verruqueux | 4556 | 15051 | 15032 | 13 |
| 6 | Carpinus betulus | 11 | Charme | 10597 | 57583 | 57554 | 12 |
| 7 | Castanea sativa | 10 | Ch<U+00E2>taignier | 6143 | 41536 | 41519 | 13 |
| 8 | Corylus avellana | 31 | Noisetier coudrier | 2631 | 9789 | 9770 | 13 |
| 9 | Crataegus monogyna | 49AM | Aub<U+00E9>pine monogyne | 1710 | 3490 | 3484 | 14 |
| 10 | Fagus sylvatica | 09 | H<U+00EA>tre | 11579 | 59824 | 59813 | 14 |
| 11 | Fraxinus excelsior | 17C | Fr<U+00EA>ne commun | 5772 | 21463 | 21406 | 14 |
| 12 | Larix decidua | 63 | M<U+00E9>l<U+00E8>ze d'Europe | 514 | 3048 | 3047 | 8 |
| 13 | Picea abies | 62 | <U+00C9>pic<U+00E9>a commun | 4423 | 33485 | 33474 | 14 |
| 14 | Pinus halepensis | 57A | Pin d'Alep | 808 | 4212 | 4210 | 3 |
| 15 | Pinus nigra | 53CA+53CO+53S+54 | Pin laricio de Calabre / Pin laricio de Corse / Pin de Salzmann / Pin noir d'Autriche | 1695 | 12521 | 12509 | 13 |
| 16 | Pinus pinaster | 51 | Pin maritime | 3177 | 23393 | 23385 | 10 |
| 17 | Pinus sylvestris | 52 | Pin sylvestre | 5634 | 33236 | 33215 | 15 |
| 18 | Populus tremula | 24 | Tremble | 2328 | 7326 | 7314 | 14 |
| 19 | Prunus avium | 22M | Merisier | 3201 | 6456 | 6450 | 13 |
| 20 | Pseudotsuga menziesii | 64 | Douglas | 2250 | 18579 | 18577 | 11 |
| 21 | Quercus ilex | 06 | Ch<U+00EA>ne vert | 1494 | 14359 | 14350 | 4 |
| 22 | Quercus petraea | 03 | Ch<U+00EA>ne sessile | 11646 | 59094 | 59079 | 13 |
| 23 | Quercus pubescens | 05 | Ch<U+00EA>ne pubescent | 4977 | 33809 | 33787 | 10 |
| 24 | Quercus robur | 02 | Ch<U+00EA>ne p<U+00E9>doncul<U+00E9> | 13876 | 55801 | 55788 | 12 |
| 25 | Robinia pseudoacacia | 14 | Robinier faux acacia | 1581 | 9030 | 9027 | 12 |
| 26 | Salix caprea | 25M | Saule marsault | 1691 | 5622 | 5597 | 14 |
| 27 | Sorbus aria | 23AB | Alisier blanc | 1174 | 2587 | 2586 | 8 |
| 28 | Sorbus torminalis | 41 | Alisier torminal | 1623 | 3013 | 3013 | 12 |
| 29 | Tilia cordata | 20P | Tilleul <U+00E0> petites feuilles | 869 | 3149 | 3149 | 12 |
| 30 | Ulmus minor | 18C | Orme champ<U+00EA>tre | 941 | 2218 | 2216 | 13 |

## PARAMS

```r
List of 30
 $ prep_version       : chr "fiveyear-v2"
 $ period             : int 5
 $ first_campaigns    : int [1:11] 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 ...
 $ csa_keep           : chr [1:3] "1" "2" "3"
 $ instp5_drop        : chr [1:2] "P" "S"
 $ def5_drop          : chr [1:5] "1" "2" "3" "4" ...
 $ gest_keep          : chr [1:3] "0" "1" "2"
 $ species_N          : int 30
 $ purity_x           : num 0.02
 $ shade_class        : logi TRUE
 $ min_dbh_cm         : num 7.5
 $ class_breaks_cm    : num [1:2] 22.5 37.5
 $ nominal_w          : num [1:3] 88.4 39.3 14.1
 $ patch_size_ha      : num 0.0707
 $ recruit_small_dbh  : num 12.5
 $ recruit_quantile   : num 0.9
 $ recruit_min_trees  : int 30
 $ dbh_growth_thresh  : num [1:2] -5 30
 $ fate_dead          : chr [1:7] "M" "5" "C" "2" ...
 $ fate_cut           : chr [1:2] "6" "7"
 $ fate_lost          : chr "N"
 $ env_rasters        : Named chr [1:6] "wc2.1_30s_bio/wc2.1_30s_bio_1.tif" "wc2.1_30s_bio/wc2.1_30s_bio_12.tif" "wc2.1_30s_bio/wc2.1_30s_bio_5.tif" "wc2.1_30s_bio/wc2.1_30s_bio_6.tif" ...
 $ block_km           : num 100
 $ holdout_block_share: num 0.2
 $ block_seed         : int 21
 $ pilot_n            : int 2000
 $ pilot_seed         : int 1
 $ holdout_pilot_n    : int 1000
 $ holdout_pilot_seed : int 11
 $ strata_breaks      : num [1:2] 0.01 0.05
```
