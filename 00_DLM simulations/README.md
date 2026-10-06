# 00_DLM simulations

This folder contains data and code that supports Appendix 4: DLM Simulations. The file "simulations.R" contains code to simulate data which is exported as "simulation_pars.rds". The quarto file "Appendix4_DLMsimulations.qmd" uses "simulation_pars.rds" to produce "Appendix4_DLMsimulations.pdf"

## Scripts

Script: simulation.R
This script simulates two random walk datasets simulation_pars.rds and simulation_pars_ar1.rds. simulation_pars.rds is used to generate Appendix 4 

Script: Appendix4_DLMsimulations.qmd
This script produces all figures, tables, and scripts in Appendix 4 using simulation_pars.rds.

## Data files

### File: simulation_pars.rds & simulation_pars_ar1.rds
Description: Simulated time series data representing white noise (simulation_pars.rds) and first-order autoregressive (simulation_pars_ar1.rds) time series processes.

*Variables*
- seed
- n_years
- covar_trend
- coef_sd
- rho
- R
- q_alpha
- q_beta
- convergence
-  rho_cov
