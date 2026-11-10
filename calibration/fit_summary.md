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
