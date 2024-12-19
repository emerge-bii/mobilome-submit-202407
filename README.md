
# Description
Scripts from "Mobilome in Mire" project. This repo consolidates code and scripts throughout the project. Note that many of them are not directly portable, ie., they need minor modification to run on a different computer. Here we provide them for the purpose of disclosure and documenting the steps performed in this analysis. Below are steps to reproduce the figures.

# Setup

## Create conda environment

```bash
git clone https://github.com/emerge-bii/mobilome-submit-202407.git
cd mobilome-submit-202407
mamba env create -f requirements.yaml
conda activate mobilome_submit_v0.1
```

## Download data package

You can download the data package from this [url](https://zenodo.org/records/14287562) and
name it as `som-data`.
```bash
curl -o som-data.tar.gz  https://zenodo.org/records/14287562/files/som-data.tar.gz?download=1
tar -xzf som-data.tar.gz
```

# Reproduce figs

```bash
./generate-fig.sh
```

You can see figs main figures: Fig.1, Fig.2, Fig.4 and supplementary figures: fig.S2 - fig.S12 in `fig.outdir`. Note that Fig.3, Fig.5, fig.S1 were manually created so are not included here.


# Package versions
See R package used and versions details below (from `sessionInfo()`):

```
R version 4.1.3 (2022-03-10)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: Ubuntu 22.04.4 LTS

Matrix products: default
BLAS/LAPACK: /mnt/ufs18/home-111/guojiaro/miniforge3/envs/mobilome_submit_v0.1/lib/libopenblasp-r0.3.28.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
 [1] ggpubr_0.6.0     hrbrthemes_0.8.0 patchwork_1.1.2  forcats_1.0.0
 [5] stringr_1.5.0    dplyr_1.1.2      purrr_1.0.1      readr_2.1.4
 [9] tidyr_1.3.0      tibble_3.2.1     ggplot2_3.4.2    tidyverse_1.3.2
[13] here_1.0.1

loaded via a namespace (and not attached):
 [1] httr_1.4.6              jsonlite_1.8.5          carData_3.0-5
 [4] modelr_0.1.11           shiny_1.7.4             fontLiberation_0.1.0
 [7] googlesheets4_1.1.0     cellranger_1.1.0        gdtools_0.3.3
[10] Rttf2pt1_1.3.12         pillar_1.9.0            backports_1.4.1
[13] glue_1.6.2              extrafontdb_1.0         digest_0.6.31
[16] promises_1.2.0.1        ggsignif_0.6.4          rvest_1.0.3
[19] colorspace_2.1-0        htmltools_0.5.5         httpuv_1.6.11
[22] gfonts_0.2.0            fontBitstreamVera_0.1.1 pkgconfig_2.0.3
[25] httpcode_0.3.0          broom_1.0.5             haven_2.5.2
[28] xtable_1.8-4            scales_1.2.1            later_1.3.1
[31] fontquiver_0.2.1        tzdb_0.4.0              timechange_0.2.0
[34] googledrive_2.1.0       car_3.1-2               generics_0.1.3
[37] ellipsis_0.3.2          withr_2.5.0             cli_3.6.1
[40] magrittr_2.0.3          crayon_1.5.2            readxl_1.4.2
[43] mime_0.12               evaluate_0.21           fs_1.6.2
[46] fansi_1.0.4             rstatix_0.7.2           xml2_1.3.4
[49] tools_4.1.3             hms_1.1.3               gargle_1.4.0
[52] lifecycle_1.0.3         munsell_0.5.0           reprex_2.0.2
[55] compiler_4.1.3          systemfonts_1.0.4       rlang_1.1.1
[58] grid_4.1.3              rmarkdown_2.22          gtable_0.3.3
[61] abind_1.4-5             DBI_1.1.3               curl_4.3.3
[64] R6_2.5.1                lubridate_1.9.2         knitr_1.43
[67] fastmap_1.1.1           extrafont_0.18          utf8_1.2.3
[70] rprojroot_2.0.3         stringi_1.7.12          crul_1.4.0
[73] Rcpp_1.0.10             vctrs_0.6.2             dbplyr_2.3.2
[76] tidyselect_1.2.0        xfun_0.39
```
