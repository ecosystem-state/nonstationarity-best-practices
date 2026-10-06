# 00_DLM simulations

This folder contains data and code that supports Appendix 4: DLM Simulations. The file "simulations.R" contains code to simulate data which is exported as "simulation_pars.rds". The quarto file "Appendix4_DLMsimulations.qmd" uses "simulation_pars.rds" to produce "Appendix4_DLMsimulations.pdf"

**Scripts**

File: simulation.R
Description: This script simulates two random walk datasets simulation_pars.rds and simulation_pars_ar1.rds. simulation_pars.rds is used to generate Appendix 4 

File: Appendix4_DLMsimulations.qmd
Description: This script produces all figures, tables, and scripts in Appendix 4 using simulation_pars.rds.

**Data Files**

File: simulation_pars.rds & simulation_pars_ar1.rds
Description: Simulated time series data representing white noise (simulation_pars.rds) and first-order autoregressive (simulation_pars_ar1.rds) time series processes.

*Variables*
- seed: seed for random process
- n_years: number of time steps
- covar_trend: covariate trend
- coef_sd: variability of random walk of the time varying covariate
- rho: autocorrelation of cov time series
- R: observation error covariance
- q_alpha: process error covariance for intercept
- q_beta: process error covariance for slope
- convergence: convergence of model. 0 indicates converged successfully; 2 indicates convergence issues
- rho_cov: latent state of the time-varying slope coefficient
