# 03_GOA_salmon_case_study

## Data Folder

File: dfa_loadings.csv

Description: Dynamic factor analysis (DFA) loadings for the latent trends used in this analysis

*Variables*
- names: common names of each of the four species of salmon included in this analysis
- mean: mean value of the DFA loading for each species
- upCI: upper bound of the confidence interval for the mean loading estimate
- lowCI: lower bound of the confidence interval for the mean loading estimate

File: dfa_model_selection_table.csv

Description: Model comparison for four DFA models with different parameterization of the observation error variance-covariance matices (R)

*Variables*
- R: observation error variance-covariance matrix description
- m: number of latent states
- loglik: log-likelihood of the fitted DFA model
- K: number of estimated parameters
- AICc: Akaike Information Criterion (AIC) corrected for small sample size for the fitted DFA model
- delAIC: difference in AIC between the fitted model and the best model

File: dfa_trend.csv

Description: latent state of DFA model

*Variables*
- t: year
- estimate: estimate of latent state representing salmon abundance
- conf.low: lower bound of the confidence interval of the latent state
- conf.high: upper bound of the confidence interval of the latent state

File: GOA_salmon_catch.csv

Description: Salmon catch estimates using a PCA

*Variables*
- catcjh_pc1: catch estimate from the first principle component
- sst_3yr_running: 3-year running mean of SST
- year: year

File: winterSST_3yr_running_mean.csv

Description: 3-year winter running mean of SST data generated from 00_process.SST.R

*Variables*
year: year
sst_3yr_running: running mean of 3 year winter SST
