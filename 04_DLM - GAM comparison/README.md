# 04_DLM - GAM comparison

We included an analysis that compared DLM and GAM output across case studies (GoA salmon and Lake Washington) for Figure 5 of the main text. This folder contains the data that was combined and the resulting output for the combined analysis. For more context on these case studies, see the associated case study folders. 

**Scripts**

File: 01_dlm_gam_comparison.qmd

Description: Script to produce side by side comparisons of DLMs and GAMs for both Lake Washington and Gulf of Alaska case studies with time-varying intercepts and time-varying intercepts and slopes. This script produces Figure 5.

**Data Files**

File: data_for_AK_salmon_example.rds

Description: Post processed data for Gulf of Alaska salmon GAMs and DLMs

*Variables*
- year: calendar year
- salmon_DFA: latent trend of salmon catch from salmon DFA
- conf.low: lower bound of the confidence interval of the salmon DFA latent trend
- conf.high: upper bound of the confidence interval of the salmon DFA latent trend
- sst_3yr_running_mean: 3-year running mean of winter SST
- z_sst: z-scored sst_3yr_running_mean

File: data_for_lakeWA_example.rds

Description: Post processed data for Lake Washington GAMs and DLMs

*Variables*
- Year: calendar year
- Month: month of year
- Temp: water temperature
- TP: total phosphorus
- pH: pH of water
- Cryptomonas: counts of Cryptomonas (cells per mL)
- Diatoms: counts of Diatoms (cells per mL)
- Greens: counts of green algae (Chlorophyta) (cells per mL)
- Bluegreens: counts of blue-green algae (cyanobacteria) (cells per mL)
- Unicells: counts of unicellular algae (cells per mL)
- Other.algae: counts of other algae (cells per mL)
- Cyclops: counts of copepod in family Cyclopidae (organisms per L)
- Daphnia: counts of daphnia (organisms per L)
- Epischura: counts of Epischura genus (organisms per L)
- Leptodora: counts of Leptodora genus (organisms per L)
- Neomysis: counts of Neomysis (organisms per L)
- Non.daphnid.cladocerans: counts of caldocerans that are not Daphnia (organisms per L)
- Non.colonial.rotifers: counts of rotifers that are not colonial (organisms per L)
- date: date denoted as yyyy-mm-dd
- y: log-transformed counts of blue-green algae (cyanobacteria) (cells per mL)
- zp: z-scored total phosphorus
- y_interp: log-transformed counts of blue-green algae (cyanobacteria) (cells per mL) wth interpolated missing data
- month_mean: mean value  log-transformed counts of blue-green algae (cyanobacteria) (cells per mL) for a given month
- y_adj: log-transformed counts of blue-green algae (cyanobacteria) with the monthly mean removed (deseasoned)
- lagged_zp: z-scored total phosphorus with a 6 month temporal lag
