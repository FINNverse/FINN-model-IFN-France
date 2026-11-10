# Validation report: fiveyear-v1-process

Written 2026-11-10 03:05 CET by calibration/validation.R; FINN 0.2.0 (main 0848b60), torch 0.15.1. Model: ../../FINN-zoo-data/frnfi/fits/fiveyear-v1-process/model.pt; release object parameters.pt sha256 f278626cc0a3d5625a8e6437d98d974d4150311bc105e7e042717391351a0218.

## Sets

- TEST = fiveyear-v1/train: 43964 sites the model was fitted on (7357 with the year-5 state masked because of a cut tree).
- HOLDOUT = fiveyear-v1/holdout_blocked: 14767 sites in 100 km blocks the model never saw (2523 state-masked). No coordinates in the model, predicted as is.
- Aggregation as fit() scores it: state at year 5; growth and mortality = mean of the annual rates over years 1..5; recruitment = sum over years 1..5.

## Metrics

| set | variable | level | n | rmse | r2 | r | bias | obs_mean | auc |
|---|---|---|---|---|---|---|---|---|---|
| test | dbh (cm, species mean) | species x site |   96600 |   6.1600 |  0.84700 | 0.9250 |  -1.43000 | 24.9000 | NA |
| test | basal area (m2 per 0.0707 ha patch) | species x site | 1320000 |   0.0419 |  0.97700 | 0.9890 |   0.00109 |  0.0513 | NA |
| test | trees (per 0.0707 ha patch) | species x site | 1320000 |   2.9100 |  0.89500 | 0.9470 |   0.33400 |  1.5000 | NA |
| test | growth (rel. dbh increment / yr) | species x site |  110000 |   0.0137 |  0.40700 | 0.6500 |  -0.00216 |  0.0174 | NA |
| test | mortality (annual rate) | species x site |  113000 |   0.1610 |  0.00634 | 0.1310 |  -0.01690 |  0.0355 | 0.673 |
| test | recruitment (stems / ha / 5 yr) | species x site | 1580000 |  32.9000 | -2.08000 | 0.0323 |   5.49000 |  1.3000 | NA |
| holdout | dbh (cm, species mean) | species x site |   31800 |   6.0900 |  0.85100 | 0.9280 |  -1.46000 | 25.1000 | NA |
| holdout | basal area (m2 per 0.0707 ha patch) | species x site |  441000 |   0.0410 |  0.97600 | 0.9880 |   0.00123 |  0.0496 | NA |
| holdout | trees (per 0.0707 ha patch) | species x site |  441000 |   2.8600 |  0.89400 | 0.9470 |   0.34200 |  1.4600 | NA |
| holdout | growth (rel. dbh increment / yr) | species x site |   36300 |   0.0138 |  0.37500 | 0.6250 |  -0.00218 |  0.0174 | NA |
| holdout | mortality (annual rate) | species x site |   37300 |   0.1580 |  0.00105 | 0.1140 |  -0.01730 |  0.0344 | 0.664 |
| holdout | recruitment (stems / ha / 5 yr) | species x site |  532000 |  32.5000 | -2.34000 | 0.0335 |   5.50000 |  1.2700 | NA |
| test | basal area (m2 per patch) | stand total |   36600 |   0.2510 |  0.94300 | 0.9720 |   0.03920 |  1.8500 | NA |
| test | trees (per patch) | stand total |   36600 |  19.0000 |  0.77900 | 0.9320 |  12.00000 | 54.0000 | NA |
| test | recruitment (stems / ha / 5 yr) | stand total |   36600 | 234.0000 | -2.77000 | 0.0114 | 194.00000 | 49.6000 | NA |
| test | basal area (m2 per patch) | stand total, PERSISTENCE reference |   36600 |   0.2930 |  0.92200 | 0.9740 |  -0.16900 |  1.8500 | NA |
| test | trees (per patch) | stand total, PERSISTENCE reference |   36600 |  12.0000 |  0.91100 | 0.9580 |   0.39200 | 54.0000 | NA |
| holdout | basal area (m2 per patch) | stand total |   12200 |   0.2480 |  0.93800 | 0.9700 |   0.04430 |  1.7900 | NA |
| holdout | trees (per patch) | stand total |   12200 |  18.8000 |  0.77600 | 0.9340 |  12.30000 | 52.6000 | NA |
| holdout | recruitment (stems / ha / 5 yr) | stand total |   12200 | 229.0000 | -3.06000 | 0.0203 | 193.00000 | 48.9000 | NA |
| holdout | basal area (m2 per patch) | stand total, PERSISTENCE reference |   12200 |   0.2870 |  0.91800 | 0.9720 |  -0.16700 |  1.7900 | NA |
| holdout | trees (per patch) | stand total, PERSISTENCE reference |   12200 |  11.7000 |  0.91400 | 0.9590 |   0.20900 | 52.6000 | NA |

A negative R2 means worse than predicting the mean. The persistence rows predict the year-5 state by the first-visit state (model-free reference).

## Mortality AUC (5-year death probability against the tallied fates; each site x species row weighted by its trees at risk)

| set | n_trees_at_risk | n_died | auc | obs_death_share | pred_death_share |
|---|---|---|---|---|---|
| test | 523000 | 11300 | 0.673 | 0.0215 | 0.0823 |
| holdout | 172000 |  3640 | 0.664 | 0.0212 | 0.0773 |

## Per species, holdout

| species_name | sites_present | ba_r2 | trees_r2 | growth_r | growth_bias | mort_auc | mort_obs | mort_pred | reg_r | reg_obs | reg_pred |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Quercus robur | 3540 | 0.983 |  0.957 | 0.6560 | -0.001920 | 0.611 | 0.01340 | 0.063600 |  0.07160 | 1.840 |   0.5000 |
| Quercus petraea | 2620 | 0.986 |  0.965 | 0.6770 | -0.001670 | 0.691 | 0.00959 | 0.056700 |  0.09320 | 1.790 |   0.4590 |
| Fagus sylvatica | 2500 | 0.990 |  0.978 | 0.6310 | -0.002540 | 0.520 | 0.00707 | 0.030500 |  0.11600 | 3.380 |   0.9750 |
| Carpinus betulus | 2280 | 0.982 |  0.961 | 0.5730 | -0.003760 | 0.568 | 0.01330 | 0.067100 |  0.08830 | 5.510 |   2.4600 |
| other | 2080 | 0.902 | -0.246 | 0.3860 |  0.000411 | 0.609 | 0.04290 | 0.216000 | -0.00614 | 3.800 | 154.0000 |
| Pinus sylvestris | 1520 | 0.983 |  0.973 | 0.6630 | -0.003190 | 0.482 | 0.01630 | 0.085600 |  0.11100 | 1.380 |   0.6200 |
| Quercus pubescens | 1490 | 0.992 |  0.984 | 0.4580 | -0.001830 | 0.557 | 0.00813 | 0.012200 |  0.10300 | 2.050 |   0.6440 |
| Fraxinus excelsior | 1480 | 0.920 |  0.877 | 0.4870 | -0.001580 | 0.462 | 0.06250 | 0.171000 |  0.08450 | 1.440 |   0.9240 |
| Castanea sativa | 1320 | 0.962 |  0.928 | 0.5020 | -0.003990 | 0.553 | 0.04440 | 0.127000 |  0.09160 | 2.540 |   1.6800 |
| Abies alba | 1230 | 0.982 |  0.741 | 0.6290 | -0.003610 | 0.542 | 0.01070 | 0.047600 |  0.14700 | 1.830 |  35.3000 |
| Betula pendula | 1040 | 0.947 |  0.911 | 0.6320 | -0.001730 | 0.457 | 0.04870 | 0.062800 |  0.12900 | 1.990 |   1.0000 |
| Picea abies |  807 | 0.965 |  0.952 | 0.6860 | -0.002160 | 0.602 | 0.03020 | 0.089400 |  0.12900 | 0.974 |   0.4270 |
| Acer campestre |  741 | 0.965 |  0.797 | 0.4280 | -0.003930 | 0.450 | 0.00968 | 0.018900 |  0.03670 | 1.170 |  17.8000 |
| Prunus avium |  739 | 0.936 |  0.934 | 0.5320 | -0.002050 | 0.452 | 0.03330 | 0.119000 |  0.03330 | 0.619 |   0.4130 |
| Pinus pinaster |  706 | 0.981 |  0.909 | 0.8640 |  0.000752 | 0.503 | 0.01500 | 0.029500 |  0.09390 | 1.480 |   0.6990 |
| Corylus avellana |  663 | 0.878 |  0.854 | 0.0738 | -0.000987 | 0.405 | 0.06600 | 0.265000 |  0.10100 | 2.940 |   2.7100 |
| Quercus ilex |  641 | 0.994 |  0.991 | 0.2420 | -0.001530 | 0.515 | 0.00382 | 0.024500 |  0.11500 | 1.430 |   0.9160 |
| Populus tremula |  562 | 0.938 |  0.865 | 0.6470 | -0.002190 | 0.454 | 0.06460 | 0.078300 |  0.05890 | 0.728 |   0.6350 |
| Pinus nigra |  534 | 0.988 |  0.986 | 0.6980 | -0.000794 | 0.585 | 0.00751 | 0.071100 |  0.05950 | 0.467 |   0.2440 |
| Crataegus monogyna |  529 | 0.895 |  0.901 | 0.2440 | -0.002820 | 0.467 | 0.03840 | 0.219000 |  0.04050 | 0.857 |   0.8600 |
| Acer pseudoplatanus |  491 | 0.967 |  0.621 | 0.6830 | -0.003410 | 0.553 | 0.01440 | 0.050500 |  0.04610 | 0.701 |  15.6000 |
| Pinus halepensis |  490 | 0.986 |  0.958 | 0.8050 | -0.004020 | 0.461 | 0.01080 | 0.000907 |  0.12600 | 0.572 |   0.0806 |
| Robinia pseudoacacia |  458 | 0.958 |  0.926 | 0.6500 | -0.005130 | 0.471 | 0.02640 | 0.144000 |  0.02770 | 1.440 |   0.6520 |
| Pseudotsuga menziesii |  394 | 0.991 |  0.967 | 0.7340 |  0.000699 | 0.691 | 0.00754 | 0.027100 |  0.08420 | 0.361 |   0.4150 |
| Alnus glutinosa |  366 | 0.972 |  0.883 | 0.6440 | -0.002590 | 0.538 | 0.02860 | 0.154000 |  0.05500 | 0.751 |   0.2660 |
| Sorbus aria |  360 | 0.940 |  0.938 | 0.3710 | -0.001270 | 0.446 | 0.03930 | 0.160000 |  0.06820 | 0.358 |   0.2400 |
| Sorbus torminalis |  358 | 0.963 |  0.949 | 0.4370 | -0.001210 | 0.593 | 0.02130 | 0.085200 |  0.03710 | 0.328 |   0.3910 |
| Salix caprea |  357 | 0.854 |  0.784 | 0.4510 | -0.003690 | 0.344 | 0.11500 | 0.273000 |  0.04880 | 0.976 |   0.5870 |
| Ilex aquifolium |  274 | 0.982 |  0.950 | 0.3530 | -0.002980 | 0.472 | 0.01370 | 0.053000 |  0.02860 | 0.481 |   0.3680 |
| Ulmus minor |  217 | 0.779 |  0.764 | 0.2480 | -0.003760 | 0.657 | 0.10800 | 0.037000 |  0.04590 | 0.594 |   0.6130 |
| Larix decidua |  215 | 0.993 |  0.986 | 0.7580 | -0.001090 | 0.683 | 0.00837 | 0.027400 |  0.09990 | 0.105 |   0.0385 |
| Quercus pyrenaica |  215 | 0.988 |  0.978 | 0.6170 |  0.002300 | 0.543 | 0.01120 | 0.028900 |  0.04990 | 0.192 |   0.0468 |
| Tilia platyphyllos |  208 | 0.987 |  0.975 | 0.5070 | -0.003830 | 0.472 | 0.00986 | 0.010500 |  0.01580 | 0.242 |   0.1370 |
| Tilia cordata |  197 | 0.985 |  0.971 | 0.6420 | -0.001310 | 0.476 | 0.00521 | 0.043100 |  0.02990 | 0.258 |   0.1920 |
| Arbutus unedo |   92 | 0.974 |  0.969 | 0.3210 | -0.002650 | 0.454 | 0.00941 | 0.104000 |  0.06340 | 0.198 |   0.2720 |
| Pinus mugo |   43 | 0.988 |  0.989 | 0.3190 |  0.000683 | 0.325 | 0.01750 | 0.069200 |  0.05080 | 0.023 |   0.0633 |

## Bare-ground equilibrium (stand totals, last 100 vs previous 100 years)

| label | ba_last100 | ba_prev100 | trees_last100 | change_pct |
|---|---|---|---|---|
| 7.0 degC, 1162 mm | 158.0 | 171 | 3500 |  -7.60 |
| 9.5 degC, 827 mm | 104.0 | 133 | 2090 | -21.80 |
| 10.3 degC, 752 mm | 103.0 | 111 | 1690 |  -6.73 |
| 11.3 degC, 744 mm |  88.7 | 108 | 1370 | -17.80 |
| 13.0 degC, 900 mm | 130.0 | 149 | 1840 | -13.10 |

## Figures

- figures/fig1_loss_epochs.png: loss per response over epochs
- figures/fig2_obs_pred.png: observed vs predicted, both sets
- figures/fig3_mortality_roc.png: mortality ROC
- figures/fig4_map.png: observed and predicted stand basal area, both sets
- figures/fig5_bare_ground.png: succession from bare ground at five climates
- figures/fig6_parameter_trajectories.png, fig6b_parameter_change.png: parameter trajectories and settling

## Provenance

| step | script | inputs | outputs |
|---|---|---|---|
| raw data | calibration/01_download_ign_raw.R | IGN DataIFN export 2005-2025 (sha256 in PROVENANCE.md) | raw/ign_data/unz, PROVENANCE.md |
| preparation | calibration/02_prepare_fiveyear_v1.R | raw tables, WorldClim 2.1 30 s | prepared/fiveyear-v1/{national,train,holdout_blocked,pilot}, PREPARATION.md, expected_preparation.csv |
| fit | calibration/fit.R | fiveyear-v1/train | fits/fiveyear-v1-process/{model.pt, parameters.pt}, loss_history.csv, fit_summary.md |
| validation | calibration/validation.R | model.pt, train, holdout_blocked | report.md, validation/*.csv, figures/*.png |

## Fit summary (from fit.R)

# Fit summary: fiveyear-v1-process

Fitted 2026-11-10 03:03 CET by calibration/fit.R; FINN 0.2.0 (main 0848b60), torch 0.15.1, R 4.5.0; 1 threads; runtime 1.63 h for 150 epochs (39.0 s/epoch).

Data: fiveyear-v1/train, 43964 sites, 36 species (incl. other). Weights: `auto` (baseline-scaled; no numeric vector, weight_scale not used).

Release assets: parameters.pt (named list of 28 tensors incl. env scaling), sha256 f278626cc0a3d5625a8e6437d98d974d4150311bc105e7e042717391351a0218; model.pt (whole fitted module), sha256 e3c52b4d830118ab3fce3511c2f68cf931adc98420f7cdca665af6c93991db05. build_model(parameters.pt) reproduced the fitted model bit for bit on 200 sites.

Final losses (baseline-scaled, last epoch):

- dbh: 0.167
- ba: 0.02212
- trees: 0.1811
- growth: 0.5887
- mortality: 1.097
- regeneration: 0.2567
- unused: 0
- epoch: 150

## PARAMS

```r
List of 20
 $ prep_version : chr "fiveyear-v1"
 $ train_set    : chr "train"
 $ holdout_set  : chr "holdout_blocked"
 $ patch_size_ha: num 0.0707
 $ n_patches    : int 1
 $ env_vars     : chr [1:2] "tavg_c" "prec_mm"
 $ recruits_dbh : num 8.3
 $ hidden       : int [1:2] 20 20
 $ dropout      : num 0.1
 $ reg_K_init   : num 50
 $ reg_floor    : num 0.001
 $ epochs       : int 150
 $ batchsize    : int 500
 $ lr           : num 0.01
 $ threads      : int 1
 $ weights      : chr "auto"
 $ loss_family  : Named chr [1:6] "mse" "mse" "poisson" "mse" ...
 $ checkpoints  : int 10
 $ device       : chr "cpu"
 $ seed         : int 1
```
