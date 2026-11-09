# Prepared data: fiveyear-v1

Written 2026-11-09 22:20 CET by IGN-FR_process/scripts/02_prepare_fiveyear_v1.R (FINN 0.2.0, data.table 1.18.6.1, terra 1.8.60).

Source: IGN, Inventaire forestier national, données brutes, campaigns 2005-2025 (export dated 2026-10-07,
downloaded 2026-11-09, checksums in raw/SHA256SUMS.txt and raw/DOWNLOAD.md), Licence Ouverte / Open Licence Etalab 2.0.
Climate: WorldClim 2.1 (1970-2000 normals, 30 arc-seconds) at the fuzzed plot centre (1 km grid cell).

## What the tables mean

- One site = one IGN plot (IDP), one patch of 0.0707 ha (the 15 m circle). Every tallied tree stands for W stems/ha
  and is expanded to round(W x 0.0707) virtual trees (largest remainder within plot x species). `ba` (m2), `trees` and
  `dbh` (cm) in obs_dt and `trees` in init_cohort_dt refer to that patch; divide ba and trees by 0.0707 for per-hectare values.
- `year` 0 = first visit (campaign_t0 in site_table.csv), 5 = second visit; `period_length` = 5. Env rows years 1..5.
- `growth` = annualised relative dbh increment of survivors, design-weighted; `growth_n` = tallied survivors behind it.
- `n_at_risk` = tallied first-visit trees that survived or died naturally (cut / not-found censored); `mort` = annual
  death probability from the design-weighted 5-year share; `n_died` = mort x n_at_risk (real-valued).
- `reg` = new trees >= 7.5 cm (<= 15 cm) per hectare over the 5 years.
- `state_masked` = plot with a cut or not-found tree: ba, dbh, trees at year 5 are NA there (FINN has no harvest).
- env (natural units; FINN z-scales): tavg_c (bio1, degC), prec_mm (bio12, mm/yr), tmax_warm_c (bio5), tmin_cold_c (bio6),
  prec_season (bio15, CV %), elev_m, x = longitude, y = latitude (WGS84, fuzzed plot centre).

## Processing log

- plots first-visited 2010-2020 with a second visit: 66021
-   land cover 1/2/3 at visit 1 (no poplar plantations): 63907
-   no planting / sowing in the interval: 63285
-   no clearing in the interval: 63170
- tree rows at visit 1: 803893; at visit 2: 841512
-   living at visit 1: 749160 (dropped: 54733 dead / fallen, 0 without C13)
-   below the 7.5 cm threshold, dropped: 436
-   fates: survivor 660677, cut 48256, dead 39240, lost 551
-   new trees at visit 2 (alive, >= 7.5 cm): 36281 kept as recruits, 432 above 15 cm dropped as missed at visit 1
-   4 ESPAR codes without a lookup entry (194 stems) lumped into 'other'
- species: 162 ESPAR codes among living visit-1 trees; 36 species cover 95.2 % of stems, 126 codes lumped into 'other'
- expansion: visit-1 basal area 24.5 m2/ha design-weighted vs 24.5 m2/ha from the virtual trees (-0.21 % off); 0 of 748724 tallied trees expand to zero
- env: WorldClim extracted for 63155 of 63170 plots (15 without a land pixel dropped)
- makeObsData / resolveSiteIDs: 58731 sites, 36 species incl. 'other'
- year-5 state (ba, dbh, trees) masked on 9880 of 58731 sites with a cut or not-found tree; growth, mortality and recruitment kept
- blocks: 16 of 82 100-km blocks held out; 43964 training sites, 14767 holdout sites
- holdout separation: nearest training site median 21.0 km, minimum 1.0 km, 10 % below 5.0 km
- set national          58731 sites, 2114316 obs rows,  730773 cohorts
- set train             43964 sites, 1582704 obs rows,  549952 cohorts
- set holdout_blocked   14767 sites,  531612 obs rows,  180821 cohorts
- set pilot               500 sites,   18000 obs rows,    6355 cohorts

## Summary (national)

- 58731 sites; mean basal area at the first visit 25.9 m2/ha (unmasked plots), mean annual relative growth 0.0155,
  pooled annual mortality 0.0215, mean recruitment 47 stems/ha per 5 yr

## Sets

- national: every site passing the filters
- train: every site in the training blocks
- holdout_blocked: every site in the held-out blocks (never for fitting)
- pilot: 500 training sites drawn with seed 1 (smoke tests only)

Spatial blocking: 100 km squares on EPSG:2154, 16 of 82 blocks held out (seed 21). Holdout site to nearest training site:
median 21.0 km, minimum 1.0 km, 10 % below 5.0 km, 90 % below 42.9 km.

Files per set: obs_dt.csv, env_dt.csv, init_cohort_dt.csv, species_dt.csv, siteID_dt.csv (siteID within the set,
siteID_national, IDP, lon/lat, block, set, state_masked). Site attributes: site_table.csv; species list with ESPAR
codes, French and Latin names: species.csv.

## Species

| species | name | ESPAR | French | stems at visit 1 | plots |
|---|---|---|---|---|---|
| 1 | Abies alba | 61 | Sapin pectin<U+00E9> | 37285 | 5243 |
| 2 | Acer campestre | 21C | <U+00C9>rable champ<U+00EA>tre | 10850 | 4095 |
| 3 | Acer pseudoplatanus | 15S | <U+00C9>rable sycomore | 8851 | 3092 |
| 4 | Alnus glutinosa | 13G | Aulne glutineux | 8594 | 1682 |
| 5 | Arbutus unedo | 40 | Arbousier | 2902 | 502 |
| 6 | Betula pendula | 12V | Bouleau verruqueux | 17689 | 5354 |
| 7 | Carpinus betulus | 11 | Charme | 62718 | 11856 |
| 8 | Castanea sativa | 10 | Ch<U+00E2>taignier | 48031 | 7437 |
| 9 | Corylus avellana | 31 | Noisetier coudrier | 13250 | 3464 |
| 10 | Crataegus monogyna | 49AM | Aub<U+00E9>pine monogyne | 5563 | 2584 |
| 11 | Fagus sylvatica | 09 | H<U+00EA>tre | 67843 | 13314 |
| 12 | Fraxinus excelsior | 17C | Fr<U+00EA>ne commun | 29326 | 7904 |
| 13 | Ilex aquifolium | 49IA | Houx | 2249 | 1024 |
| 14 | Larix decidua | 63 | M<U+00E9>l<U+00E8>ze d'Europe | 4091 | 698 |
| 15 | Picea abies | 62 | <U+00C9>pic<U+00E9>a commun | 37021 | 5137 |
| 16 | Pinus halepensis | 57A | Pin d'Alep | 4712 | 944 |
| 17 | Pinus mugo | 58 | Pin <U+00E0> crochets | 2286 | 319 |
| 18 | Pinus nigra | 53CO+54 | Pin laricio de Corse / Pin noir d'Autriche | 13795 | 1991 |
| 19 | Pinus pinaster | 51 | Pin maritime | 26179 | 3761 |
| 20 | Pinus sylvestris | 52 | Pin sylvestre | 38237 | 6676 |
| 21 | Populus tremula | 24 | Tremble | 8363 | 2723 |
| 22 | Prunus avium | 22M | Merisier | 8061 | 3947 |
| 23 | Pseudotsuga menziesii | 64 | Douglas | 20182 | 2554 |
| 24 | Quercus ilex | 06 | Ch<U+00EA>ne vert | 22276 | 2351 |
| 25 | Quercus petraea | 03 | Ch<U+00EA>ne sessile | 64923 | 13028 |
| 26 | Quercus pubescens | 05 | Ch<U+00EA>ne pubescent | 43922 | 6677 |
| 27 | Quercus pyrenaica | 07 | Ch<U+00EA>ne tauzin | 2129 | 693 |
| 28 | Quercus robur | 02 | Ch<U+00EA>ne p<U+00E9>doncul<U+00E9> | 65614 | 16670 |
| 29 | Robinia pseudoacacia | 14 | Robinier faux acacia | 11297 | 2039 |
| 30 | Salix caprea | 25M | Saule marsault | 7251 | 2171 |
| 31 | Sorbus aria | 23AB | Alisier blanc | 3943 | 1728 |
| 32 | Sorbus torminalis | 41 | Alisier torminal | 3711 | 1980 |
| 33 | Tilia cordata | 20P | Tilleul <U+00E0> petites feuilles | 3677 | 1027 |
| 34 | Tilia platyphyllos | 20G | Tilleul <U+00E0> grandes feuilles | 2754 | 858 |
| 35 | Ulmus minor | 18C | Orme champ<U+00EA>tre | 3302 | 1354 |
| 36 | other | 04+12P+25C+73+21O+19+21M+23SO+08S+33N+25B+49SN+15P+49BS+17O+23PM+69JO+69JC+65+49PS+33G+74J+72V+25R+49FL+23PC+27C+55+56+23C+13B+18M+49PM+17F+38AU+49EA+32+34+22G+23PD+25FR+28+33B+67+59+39+49RC+77+29AF+49BO+71+29MA+49EV+29EN+22S+13C+68TG+53S+49LN+68PM+68CL+57B+68CE+49AE+36+74H+68TH+26OC+49RA+23PA+27N+49PT+29MI+26E+68CJ+23AM+29AI+18D+29CM+23F+16+37+68SV+42+23AF+66+49PL+53CA+25V+22C+68SC+20X+38AL+49RP+68PC+68CM+49PC+49SR+68CA+23A+25D+25E3+69+70SC+49PB+70SE+23PF+49FA+76+49AA+68CC+23SS+25FD+25XR+49CS+49MB+49RT+25E5+29FI+29PT+49C+49MN+49TA+49VV+68CB+72N | all other species | 35847 |  |

## PARAMS

```r
List of 22
 $ prep_version       : chr "fiveyear-v1"
 $ period             : int 5
 $ first_campaigns    : int [1:11] 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 ...
 $ csa_keep           : chr [1:3] "1" "2" "3"
 $ instp5_drop        : chr [1:2] "P" "S"
 $ def5_drop          : chr [1:5] "1" "2" "3" "4" ...
 $ min_dbh_cm         : num 7.5
 $ class_breaks_cm    : num [1:2] 22.5 37.5
 $ nominal_w          : num [1:3] 88.4 39.3 14.1
 $ patch_size_ha      : num 0.0707
 $ recruit_max_dbh    : num 15
 $ dbh_growth_thresh  : num [1:2] -5 30
 $ fate_dead          : chr [1:7] "M" "5" "C" "2" ...
 $ fate_cut           : chr [1:2] "6" "7"
 $ fate_lost          : chr "N"
 $ species_quantile   : num 0.95
 $ env_rasters        : Named chr [1:6] "wc2.1_30s_bio/wc2.1_30s_bio_1.tif" "wc2.1_30s_bio/wc2.1_30s_bio_12.tif" "wc2.1_30s_bio/wc2.1_30s_bio_5.tif" "wc2.1_30s_bio/wc2.1_30s_bio_6.tif" ...
 $ block_km           : num 100
 $ holdout_block_share: num 0.2
 $ block_seed         : int 21
 $ pilot_n            : int 500
 $ pilot_seed         : int 1
```
