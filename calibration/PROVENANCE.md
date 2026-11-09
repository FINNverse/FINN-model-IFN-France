# Provenance of every raw layer

Pinned vintage: IGN export dated 2026-10-07 (campaigns 2005-2025), downloaded 2026-11-09. The IGN archive is a
live download replaced every year; a later archive with another checksum is another vintage and needs a new
prepared-data version.

## IGN raw data (01_download_ign_raw.R)

Licence: Licence Ouverte / Open Licence Etalab 2.0 (https://www.etalab.gouv.fr/wp-content/uploads/2017/04/ETALAB-Licence-Ouverte-v2.0.pdf)

Citation required by IGN: IGN - Inventaire forestier national français, Données brutes, Campagnes annuelles 2005 et suivantes, https://inventaire-forestier.ign.fr/dataIFN/, site consulté le 09/11/2026

| file | URL | SHA-256 | size (bytes) | downloaded |
| export_dataifn_2005_2025.zip | https://inventaire-forestier.ign.fr/dataifn/data/export_dataifn_2005_2025.zip | be583e7eb01d43d4e74cc46097b8900092f1dcbfbe355c0ee3080001157c5b70 | 71136669 | 2026-11-09 (kept) |
| export_dataifn_doc_2025.zip | https://inventaire-forestier.ign.fr/dataifn/data/export_dataifn_doc_2025.zip | 8af9111099a2d59804bab8b696bb6deacc7aa1fced7725600a4e012b930f2480 | 6207768 | 2026-11-09 (kept) |

## Climate: WorldClim 2.1 (02_prepare_fiveyear_v1.R)

Source: https://www.worldclim.org/data/worldclim21.html (Fick & Hijmans 2017), 1970-2000 normals, 30 arc-seconds,
files wc2.1_30s_bio.zip and wc2.1_30s_elev.zip, licence CC BY-SA 4.0; local copies pulled 2019-09-09 (file dates),
values extracted at the fuzzed IGN plot centres (1 km cell) on the date of PREPARATION.md.

| variable | file | SHA-256 | size (bytes) |
|---|---|---|---|
| tavg_c | wc2.1_30s_bio_1.tif | 382c991cea4a20b42e345fcd5a897788feee05bc0d613dec032ff0cad1f4c375 | 753795566 |
| prec_mm | wc2.1_30s_bio_12.tif | 6df0ef2f57f900733a569fcc4c8debd27b490978c5d82c2970485f0e66abbf28 | 487101349 |
| tmax_warm_c | wc2.1_30s_bio_5.tif | dd770c326664dcd16d8799fdfe28adb65f5e054c7c3692452e1382eba56a37c3 | 378786176 |
| tmin_cold_c | wc2.1_30s_bio_6.tif | fb5532167aa42fed22f03e552a65d1055acd35a946b9c114c903fb5db1218c82 | 398541228 |
| prec_season | wc2.1_30s_bio_15.tif | 9bfde1c593b72b707f24b52d32da292f29a77c5cc5d7bbbb81832df7440296c9 | 1006627399 |
| elev_m | wc2.1_30s_elev.tif | 1d57ce9fd3584cea6a227396bc75ef88522fd5346e8487c17370885c0499716e | 342195708 |

