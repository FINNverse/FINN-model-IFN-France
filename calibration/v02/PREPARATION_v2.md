# Prepared data: fiveyear-v2

Written 2026-11-10 06:27 CET by calibration/v02/02_prepare_fiveyear_v2.R (FINN 0.2.0). Raw layers as v1 (calibration/PROVENANCE.md).

## Changes against fiveyear-v1

1. Mortality target = TOTAL REMOVAL (dead, fallen, broken or cut); not-found censored; no state masked. Covariate `gest` (IGN GEST at visit 1: 0 none, 1 weak/old, 2 clear signs of management). Plots without GEST dropped.
   Natural-only diagnostics: `mort_natural`, `n_died_natural`, `cut_share_5yr` (W-weighted share of the at-risk trees cut).
2. Recruits = new trees whose back-modelled start diameter (dbh_t1 - species q90 of the small-tree 5-yr increment) is below 7.5 cm; other new trees dropped as missed at visit 1.
3. 'other' split into 'other conifer' and 'other broadleaf'.
4. tree_obs_dt.csv per set: tallied trees at risk (siteID, year 5, patchID 1, species, dbh at visit 1, trees = virtual copies, annual relative growth of survivors, died = removal 0/1, died_natural, cut, period_length 5).
5. Sets pilot2000 (seed 1, training blocks) and holdout1000 (seed 11, holdout blocks) for comparable pilots.

Everything else as v1 (see calibration/PREPARATION.md): one 0.0707 ha patch per plot, W x 0.0707 virtual trees, tallied n_at_risk / growth_n, annualised rates, reg per ha per 5 yr, env years 1..5, WorldClim 2.1 normals in natural units, 100 km blocks (seed 21).

## Processing log

- plots first-visited 2010-2020 with a second visit: 66021
-   land cover 1/2/3 at visit 1: 63907
-   no planting / sowing in the interval: 63285
-   no clearing in the interval: 63170
-   management trace GEST recorded (0/1/2): 62919  [v2]
-   living visit-1 trees >= 7.5 cm: 747100 (436 below the threshold dropped)
-   fates: survivor 659139, cut 48254, dead 39165, lost 542
-   new trees at visit 2 (alive, >= 7.5 cm): 36616; recruits by the back-modelled start diameter (dbh_t1 - species q90 of the small-tree 5-yr increment, pooled 2.61 cm, < 7.5 cm): 33321; 3295 dropped as trees missed at visit 1  [v2]
-   dropped new trees by dbh: (7.5,10] 724, (10,12.5] 1774, (12.5,15] 459, (15,22.5] 330, (22.5,Inf] 8
- species: 162 ESPAR codes; 36 species cover 95.2 % of stems; the rest -> 'other conifer' (35 codes, 6999 stems) and 'other broadleaf' (91 codes, 28505 stems)  [v2]
- expansion: basal area 24.6 (design) vs 24.5 m2/ha (virtual), -0.21 % off; 0 tallied trees expand to zero
- env: WorldClim for 62904 of 62919 plots
- makeObsData / resolveSiteIDs: 59543 sites, 37 species classes
- removal target: pooled annual removal 0.0600 (natural 0.0198, cut share of the at-risk weight over 5 yr 0.062)
- blocks: 16 of 82 held out; 44580 training, 14963 holdout sites; pilot2000 / holdout1000 drawn
- holdout separation: median 21.0 km, min 1.0 km, 10 % below 5.0 km
- set national          59543 sites, 2203091 obs rows,  741175 cohorts,  740701 tallied trees
- set train             44580 sites, 1649460 obs rows,  557848 cohorts,  557484 tallied trees
- set holdout_blocked   14963 sites,  553631 obs rows,  183327 cohorts,  183217 tallied trees
- set pilot2000          2000 sites,   74000 obs rows,   25362 cohorts,   25339 tallied trees
- set holdout1000        1000 sites,   37000 obs rows,   12308 cohorts,   12300 tallied trees

## Summary (national)

- 59543 sites; mean basal area at year 5 25.1 m2/ha, 704 stems/ha; pooled annual removal 0.0600 (natural 0.0198); mean recruitment 42.2 stems/ha per 5 yr

## Species

| species | name | ESPAR | French | stems at visit 1 | plots |
|---|---|---|---|---|---|
| 1 | Abies alba | 61 | Sapin pectin<U+00E9> | 37216 | 5234 |
| 2 | Acer campestre | 21C | <U+00C9>rable champ<U+00EA>tre | 10845 | 4091 |
| 3 | Acer pseudoplatanus | 15S | <U+00C9>rable sycomore | 8793 | 3086 |
| 4 | Alnus glutinosa | 13G | Aulne glutineux | 8580 | 1679 |
| 5 | Arbutus unedo | 40 | Arbousier | 2891 | 497 |
| 6 | Betula pendula | 12V | Bouleau verruqueux | 17652 | 5347 |
| 7 | Carpinus betulus | 11 | Charme | 62692 | 11851 |
| 8 | Castanea sativa | 10 | Ch<U+00E2>taignier | 48026 | 7435 |
| 9 | Corylus avellana | 31 | Noisetier coudrier | 13240 | 3461 |
| 10 | Crataegus monogyna | 49AM | Aub<U+00E9>pine monogyne | 5535 | 2574 |
| 11 | Fagus sylvatica | 09 | H<U+00EA>tre | 67741 | 13298 |
| 12 | Fraxinus excelsior | 17C | Fr<U+00EA>ne commun | 29275 | 7883 |
| 13 | Ilex aquifolium | 49IA | Houx | 2248 | 1023 |
| 14 | Larix decidua | 63 | M<U+00E9>l<U+00E8>ze d'Europe | 3991 | 678 |
| 15 | Picea abies | 62 | <U+00C9>pic<U+00E9>a commun | 36962 | 5122 |
| 16 | Pinus halepensis | 57A | Pin d'Alep | 4700 | 939 |
| 17 | Pinus mugo | 58 | Pin <U+00E0> crochets | 2128 | 296 |
| 18 | Pinus nigra | 53CO+54 | Pin laricio de Corse / Pin noir d'Autriche | 13782 | 1985 |
| 19 | Pinus pinaster | 51 | Pin maritime | 26153 | 3759 |
| 20 | Pinus sylvestris | 52 | Pin sylvestre | 38121 | 6646 |
| 21 | Populus tremula | 24 | Tremble | 8361 | 2722 |
| 22 | Prunus avium | 22M | Merisier | 8054 | 3943 |
| 23 | Pseudotsuga menziesii | 64 | Douglas | 20182 | 2554 |
| 24 | Quercus ilex | 06 | Ch<U+00EA>ne vert | 22165 | 2333 |
| 25 | Quercus petraea | 03 | Ch<U+00EA>ne sessile | 64879 | 13015 |
| 26 | Quercus pubescens | 05 | Ch<U+00EA>ne pubescent | 43803 | 6647 |
| 27 | Quercus pyrenaica | 07 | Ch<U+00EA>ne tauzin | 2129 | 693 |
| 28 | Quercus robur | 02 | Ch<U+00EA>ne p<U+00E9>doncul<U+00E9> | 65611 | 16668 |
| 29 | Robinia pseudoacacia | 14 | Robinier faux acacia | 11297 | 2039 |
| 30 | Salix caprea | 25M | Saule marsault | 7242 | 2168 |
| 31 | Sorbus aria | 23AB | Alisier blanc | 3899 | 1712 |
| 32 | Sorbus torminalis | 41 | Alisier torminal | 3711 | 1980 |
| 33 | Tilia cordata | 20P | Tilleul <U+00E0> petites feuilles | 3656 | 1021 |
| 34 | Tilia platyphyllos | 20G | Tilleul <U+00E0> grandes feuilles | 2738 | 854 |
| 35 | Ulmus minor | 18C | Orme champ<U+00EA>tre | 3298 | 1352 |
| 36 | other broadleaf | 04+12P+25C+21O+19+21M+23SO+08S+33N+25B+49SN+15P+49BS+17O+23PM+49PS+33G+25R+49FL+23PC+27C+23C+13B+18M+49PM+17F+38AU+49EA+32+34+22G+23PD+25FR+33B+28+39+49RC+29AF+49BO+29MA+49EV+29EN+22S+13C+49LN+49AE+36+26OC+49RA+23PA+27N+49PT+29MI+26E+23AM+29AI+18D+29CM+23F+16+42+23AF+49PL+22C+37+20X+38AL+49RP+49PC+49SR+23A+25D+25E3+23PF+25V+49FA+49AA+49PB+23SS+25XR+49CS+49MB+49RT+25E5+25FD+29FI+29PT+49C+49MN+49TA+49VV | all other broadleaf species | 28505 |  |
| 37 | other conifer | 73+69JO+69JC+65+74J+72V+55+56+67+59+77+71+68TG+53S+68PM+68CL+57B+68CE+74H+68TH+68CJ+68SV+53CA+66+68SC+68PC+68CM+68CA+69+70SC+70SE+76+68CC+68CB+72N | all other conifer species | 6999 |  |

## PARAMS

```r
List of 28
 $ prep_version       : chr "fiveyear-v2"
 $ period             : int 5
 $ first_campaigns    : int [1:11] 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 ...
 $ csa_keep           : chr [1:3] "1" "2" "3"
 $ instp5_drop        : chr [1:2] "P" "S"
 $ def5_drop          : chr [1:5] "1" "2" "3" "4" ...
 $ gest_keep          : chr [1:3] "0" "1" "2"
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
 $ species_quantile   : num 0.95
 $ conifer_first_digit: chr [1:3] "5" "6" "7"
 $ env_rasters        : Named chr [1:6] "wc2.1_30s_bio/wc2.1_30s_bio_1.tif" "wc2.1_30s_bio/wc2.1_30s_bio_12.tif" "wc2.1_30s_bio/wc2.1_30s_bio_5.tif" "wc2.1_30s_bio/wc2.1_30s_bio_6.tif" ...
 $ block_km           : num 100
 $ holdout_block_share: num 0.2
 $ block_seed         : int 21
 $ pilot_n            : int 2000
 $ pilot_seed         : int 1
 $ holdout_pilot_n    : int 1000
 $ holdout_pilot_seed : int 11
```
