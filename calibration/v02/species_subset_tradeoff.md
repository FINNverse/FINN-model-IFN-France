# Species set S and plot purity x: the trade-off

Written 2026-11-10 06:57 by calibration/v02/species_subset_tradeoff.R. Cohort: the v1 plot filters plus a GEST record, 58528 plots, 16 holdout blocks. S = the N most frequent species (living tallied stems at the first visit; Pinus nigra's two codes as one species). A plot is kept when species outside S hold < x % of its basal area at BOTH visits (W-weighted). `holdout_blocks_populated` = held-out 100 km blocks with >= 20 retained plots; `*_envelope_share` = width of the retained 2.5-97.5 % climate range relative to all plots; `stems_S_lost_share` = stems of S on dropped plots; `species_min_plots` = plots of the rarest species of S among the retained plots.

| N | x_pct | plots | plots_share | holdout_plots | holdout_blocks_populated | holdout_blocks_total | tavg_envelope_share | prec_envelope_share | end_share_out_median | end_share_gt5pct | stems_S_lost_share | species_min_plots | species_min_trees_at_risk | species_min_holdout_blocks | species_min_tavg_span | rarest_species |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 15 | 1 | 28800 | 0.492 |  7370 | 14 | 16 | 0.912 | 0.996 | 0 | 0.01580 | 0.4140 | 1040 | 5660 | 3 | 0.369 | Quercus ilex |
| 15 | 2 | 29700 | 0.508 |  7620 | 14 | 16 | 0.913 | 0.999 | 0 | 0.01590 | 0.3880 | 1080 | 6220 | 3 | 0.369 | Quercus ilex |
| 15 | 3 | 30800 | 0.526 |  7860 | 14 | 16 | 0.913 | 1.000 | 0 | 0.01850 | 0.3630 | 1140 | 6810 | 3 | 0.370 | Quercus ilex |
| 15 | 5 | 32600 | 0.557 |  8280 | 14 | 16 | 0.912 | 1.000 | 0 | 0.03310 | 0.3200 | 1200 | 7940 | 3 | 0.371 | Quercus ilex |
| 20 | 1 | 34300 | 0.586 |  8620 | 14 | 16 | 0.906 | 1.000 | 0 | 0.01150 | 0.3520 |  716 | 4080 | 3 | 0.372 | Alnus glutinosa |
| 20 | 2 | 35300 | 0.603 |  8870 | 14 | 16 | 0.906 | 1.010 | 0 | 0.01160 | 0.3260 |  747 | 4300 | 3 | 0.372 | Alnus glutinosa |
| 20 | 3 | 36300 | 0.621 |  9110 | 14 | 16 | 0.905 | 1.010 | 0 | 0.01340 | 0.3010 |  770 | 4470 | 3 | 0.370 | Alnus glutinosa |
| 20 | 5 | 38200 | 0.653 |  9520 | 14 | 16 | 0.904 | 1.000 | 0 | 0.02630 | 0.2580 |  839 | 4900 | 3 | 0.370 | Alnus glutinosa |
| 25 | 1 | 41700 | 0.712 | 10600 | 14 | 16 | 0.935 | 0.993 | 0 | 0.00733 | 0.2420 |  786 | 2740 | 3 | 0.372 | Pinus halepensis |
| 25 | 2 | 42600 | 0.727 | 10900 | 14 | 16 | 0.935 | 0.994 | 0 | 0.00762 | 0.2200 |  793 | 2840 | 3 | 0.373 | Pinus halepensis |
| 25 | 3 | 43500 | 0.744 | 11100 | 14 | 16 | 0.934 | 0.996 | 0 | 0.00895 | 0.1980 |  805 | 2960 | 3 | 0.372 | Pinus halepensis |
| 25 | 5 | 45100 | 0.770 | 11400 | 14 | 16 | 0.934 | 0.996 | 0 | 0.01770 | 0.1640 |  823 | 3200 | 3 | 0.372 | Pinus halepensis |
| 30 | 1 | 46300 | 0.792 | 11700 | 15 | 16 | 0.947 | 0.996 | 0 | 0.00553 | 0.1680 |  499 | 2120 | 3 | 0.372 | Larix decidua |
| 30 | 2 | 47000 | 0.803 | 11900 | 15 | 16 | 0.947 | 0.996 | 0 | 0.00567 | 0.1510 |  509 | 2160 | 3 | 0.370 | Larix decidua |
| 30 | 3 | 47700 | 0.815 | 12100 | 15 | 16 | 0.949 | 0.996 | 0 | 0.00647 | 0.1360 |  519 | 2220 | 3 | 0.370 | Larix decidua |
| 30 | 5 | 48800 | 0.833 | 12300 | 15 | 16 | 0.949 | 0.996 | 0 | 0.01180 | 0.1140 |  535 | 2300 | 3 | 0.372 | Larix decidua |
| 35 | 1 | 48800 | 0.833 | 12300 | 15 | 16 | 0.947 | 0.991 | 0 | 0.00499 | 0.1320 |  262 | 1520 | 3 | 0.369 | Arbutus unedo |
| 35 | 2 | 49300 | 0.842 | 12400 | 15 | 16 | 0.950 | 0.993 | 0 | 0.00508 | 0.1190 |  286 | 1610 | 3 | 0.373 | Arbutus unedo |
| 35 | 3 | 49800 | 0.852 | 12600 | 15 | 16 | 0.951 | 0.994 | 0 | 0.00553 | 0.1070 |  293 | 1640 | 3 | 0.377 | Arbutus unedo |
| 35 | 5 | 50700 | 0.867 | 12800 | 15 | 16 | 0.952 | 0.994 | 0 | 0.00949 | 0.0888 |  320 | 1680 | 3 | 0.377 | Arbutus unedo |
| 40 | 1 | 50900 | 0.870 | 12900 | 15 | 16 | 0.966 | 0.996 | 0 | 0.00427 | 0.1060 |  263 | 1380 | 3 | 0.369 | Arbutus unedo |
| 40 | 2 | 51400 | 0.878 | 13100 | 15 | 16 | 0.968 | 0.997 | 0 | 0.00436 | 0.0925 |  277 | 1460 | 3 | 0.372 | Pinus mugo |
| 40 | 3 | 51900 | 0.887 | 13200 | 15 | 16 | 0.969 | 0.997 | 0 | 0.00478 | 0.0815 |  278 | 1500 | 3 | 0.372 | Pinus mugo |
| 40 | 5 | 52700 | 0.901 | 13400 | 15 | 16 | 0.970 | 0.997 | 0 | 0.00800 | 0.0657 |  280 | 1600 | 3 | 0.377 | Pinus mugo |
| 50 | 1 | 53700 | 0.917 | 13600 | 15 | 16 | 0.984 | 1.000 | 0 | 0.00328 | 0.0672 |  137 |  636 | 2 | 0.375 | Fraxinus angustifolia |
| 50 | 2 | 54000 | 0.923 | 13700 | 15 | 16 | 0.984 | 1.000 | 0 | 0.00334 | 0.0585 |  140 |  643 | 2 | 0.375 | Fraxinus angustifolia |
| 50 | 3 | 54400 | 0.929 | 13800 | 15 | 16 | 0.984 | 1.000 | 0 | 0.00357 | 0.0507 |  143 |  651 | 2 | 0.377 | Fraxinus angustifolia |
| 50 | 5 | 55000 | 0.939 | 13900 | 15 | 16 | 0.985 | 1.000 | 0 | 0.00607 | 0.0395 |  146 |  660 | 2 | 0.376 | Fraxinus angustifolia |

## Species by rank

| rank | species | stems | plots | ESPAR |
|---|---|---|---|---|
| 1 | Fagus sylvatica | 66532 | 13075 | 09 |
| 2 | Quercus robur | 64399 | 16344 | 02 |
| 3 | Quercus petraea | 64118 | 12819 | 03 |
| 4 | Carpinus betulus | 61745 | 11666 | 11 |
| 5 | Castanea sativa | 46083 | 7219 | 10 |
| 6 | Quercus pubescens | 43192 | 6556 | 05 |
| 7 | Pinus sylvestris | 37457 | 6529 | 52 |
| 8 | Abies alba | 36471 | 5138 | 61 |
| 9 | Picea abies | 34523 | 4901 | 62 |
| 10 | Fraxinus excelsior | 28724 | 7747 | 17C |
| 11 | Pinus pinaster | 24532 | 3538 | 51 |
| 12 | Quercus ilex | 21794 | 2293 | 06 |
| 13 | Pseudotsuga menziesii | 19519 | 2481 | 64 |
| 14 | Betula pendula | 17091 | 5195 | 12V |
| 15 | Pinus nigra | 13586 | 1038 | 53CO+54+53S+53CA |
| 16 | Corylus avellana | 12937 | 3394 | 31 |
| 17 | Robinia pseudoacacia | 10896 | 1978 | 14 |
| 18 | Acer campestre | 10681 | 4033 | 21C |
| 19 | Acer pseudoplatanus | 8691 | 3046 | 15S |
| 20 | Alnus glutinosa | 8384 | 1648 | 13G |
| 21 | Populus tremula | 8187 | 2659 | 24 |
| 22 | Prunus avium | 7933 | 3874 | 22M |
| 23 | Salix caprea | 7129 | 2125 | 25M |
| 24 | Crataegus monogyna | 5432 | 2528 | 49AM |
| 25 | Pinus halepensis | 4623 | 926 | 57A |
| 26 | Larix decidua | 3930 | 671 | 63 |
| 27 | Sorbus aria | 3869 | 1694 | 23AB |
| 28 | Sorbus torminalis | 3653 | 1950 | 41 |
| 29 | Tilia cordata | 3608 | 1006 | 20P |
| 30 | Ulmus minor | 3231 | 1321 | 18C |
| 31 | Arbutus unedo | 2810 | 488 | 40 |
| 32 | Tilia platyphyllos | 2724 | 846 | 20G |
| 33 | Ilex aquifolium | 2206 | 996 | 49IA |
| 34 | Populus | 2192 | 457 | 19+33G |
| 35 | Quercus rubra | 2096 | 445 | 04 |
| 36 | Pinus mugo | 2093 | 292 | 58 |
| 37 | Quercus pyrenaica | 2074 | 675 | 07 |
| 38 | Betula pubescens | 2032 | 418 | 12P |
| 39 | Salix cinerea | 1965 | 497 | 25C |
| 40 | Acer opalus | 1732 | 581 | 21O |
| 41 | Picea sitchensis | 1594 | 252 | 73 |
| 42 | Acer monspessulanum | 1518 | 560 | 21M |
| 43 | Sorbus aucuparia | 1364 | 564 | 23SO |
| 44 | Quercus suber | 1237 | 296 | 08S |
| 45 | Salix alba | 905 | 233 | 25B |
| 46 | Populus nigra | 904 | 281 | 33N |
| 47 | Sambucus nigra | 850 | 419 | 49SN |
| 48 | Acer platanoides | 824 | 401 | 15P |
| 49 | Buxus sempervirens | 805 | 246 | 49BS |
| 50 | Fraxinus angustifolia | 772 | 182 | 17O |
