# 03_GOA_salmon_case_study

This folder contains data and code that supports Appendix 3: Gulf of Alaska salmon analysis and all associated data processing, analyses, and figure production used in the main text that applies to the Gulf of Alaska case study.

**Data Folder**

The data folder contains the processed and raw data files. **NOTE** the raw SST data that is processed from the process.SST.R script is too large for github and the .nc is not containted on this repository. SST data is from the [National Centers for Environmental Information Extended Reconstructed Sea Surface Temperature]((https://www.ncei.noaa.gov/products/extended-reconstructed-sst)). More information of data files are contained within that folder.

**Scripts**

File: 00_process.SST.R

Description: post process sea surface temperature netcdf to generate 3-year running mean in SST that is exported as the "winterSST_3yr_running_mean.csv" file that is used as the temperature covariate

File: 01_Appendix3_GoASalmon.qmd

Description: This script produces all analyses for Gulf of Alaska salmon analysis from the post-processed files "winterSST_3yr_running_mean.rds" and "dfa_trend.rds" to produce figures, tables, and scripts in Appendix 3
