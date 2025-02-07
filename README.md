# Analysis for 
# _Multivariate signals of population collapse in a high-throughput ecological experiment_
#### Cerini, F., Jackson, J., O'Brien, D., Childs, D.Z., & Clements, C.F.

#### 2025-02-07
#### Repository created by John Jackson and Francesco Cerini

---

[![DOI](https://zenodo.org/badge/721133613.svg)](https://zenodo.org/doi/10.5281/zenodo.10160252)

This directory contains scripts and analysis data for our work on an empirical demonstration of a sequence of signals preceding a population collapse in the ciliate _Paramecium caudatum_. For the manuscript please see [the Authorea entry](XXXX.XXXX.XXXX). Package version info for this analysis is given below.

Analysis scripts can be found in the `scripts/` sub-repository, figures in the `output/` sub-repository, and analysis data in the `data/` sub-repository. 

Scripts are labeled A and B in order of the analysis, and are as follows:

1. `A_data_cleaning.R` - Data cleaning, seasonal decomposition and supplementary figures.
2. `B_autocorrelation_k_selection.R` - Preparatory analysis selecting for k (basis dimension in GAM) and exploring autocorrelation
3. `C_gam_ews.R` - Additive model analysis for P.caudatum, include timeline component models and EWS analysis
4. `D_piecewise_regression.R`- Secondary analysis approach using piecewise Bayesian linear regression/

## System Information and Package Versions

<details>
  <summary>Click here to expand</summary>
  
```
R version 4.3.1 (2023-06-16)
Platform: x86_64-apple-darwin20 (64-bit)
Running under: macOS Monterey 12.6.9

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggpubr_0.6.0     patchwork_1.1.2  MASS_7.3-60      mgcv_1.9-0       nlme_3.1-162     EWSmethods_1.1.2 lubridate_1.9.2  forcats_1.0.0   
 [9] stringr_1.5.0    dplyr_1.1.2      purrr_1.0.1      readr_2.1.4      tidyr_1.3.0      tibble_3.2.1     ggplot2_3.4.2    tidyverse_2.0.0 

loaded via a namespace (and not attached):
 [1] gtable_0.3.3       rstatix_0.7.2      lattice_0.21-8     tzdb_0.4.0         quadprog_1.5-8     vctrs_0.6.3        tools_4.3.1       
 [8] generics_0.1.3     curl_5.0.2         parallel_4.3.1     fansi_1.0.4        xts_0.13.1         pkgconfig_2.0.3    Matrix_1.6-1.1    
[15] lifecycle_1.0.3    compiler_4.3.1     munsell_0.5.0      codetools_0.2-19   carData_3.0-5      pillar_1.9.0       car_3.1-2         
[22] seasonal_1.9.0     iterators_1.0.14   abind_1.4-5        foreach_1.5.2      fracdiff_1.5-2     tidyselect_1.2.0   stringi_1.7.12    
[29] tseries_0.10-54    splines_4.3.1      grid_4.3.1         colorspace_2.1-0   cli_3.6.1          magrittr_2.0.3     utf8_1.2.3        
[36] broom_1.0.5        withr_2.5.0        scales_1.2.1       backports_1.4.1    forecast_8.21      x13binary_1.1.57-3 timechange_0.2.0  
[43] TTR_0.24.3         infotheo_1.2.0.1   quantmod_0.4.24    nnet_7.3-19        timeDate_4022.108  ggsignif_0.6.4     zoo_1.8-12        
[50] hms_1.1.3          urca_1.3-3         lmtest_0.9-40      rlang_1.1.1        Rcpp_1.0.11        glue_1.6.2         rstudioapi_0.15.0 
[57] R6_2.5.1                
```

</details>
