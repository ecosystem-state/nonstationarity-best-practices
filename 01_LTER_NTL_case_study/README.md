# 01_LTER_NTL_case_study

This folder contains data and code that supports Appendix 1: Trout Lake and all associated data processing, analyses, and figure production used in the main text that applies to the Trout Lake case study.

**Scripts**

File: 00_LTER-NTL_DataProcessing.R

Description: Process the data that was accessed at [North Temperate Lakes Long-Term Ecological Research](https://lter.limnology.wisc.edu/core-datasets/) and contained in the "NTL-LTER-TroutLake" folder which were accessed on July, 30 2024 (see folder for all additional information). This script produces the six "NTL_LTER_TR_.rds" objects which are used by "02_Appendix1_TroutLake.qmd" for the Trout Lake portion of the analysis and Appendix 1.

File: 01_Appendix1_TroutLake.qmd

Description: This script produces all analyses for Trout Lake from the post-processed file "NTL_LTER_TL_data.rds" to produce figures, tables, and scripts in Appendix 1 in addition to Figure 3 included in the main text.

**Data Files**

File: NTL_LTER_TL_data.rds

Description: Post processed data from [North Temperate Lakes Long-Term Ecological Research](https://lter.limnology.wisc.edu/core-datasets/) to generate annual means and associated confidence intervals. 

*Variables*
- year: year of observation remaned from original name of foles 'year4'
- mean_sec: mean annual secchi depth
- sd_sec: standard deviation of annual secchi depth
- ci_sec_lwr: lower bound of confidence interval of annual secchi depth
- ci_sec_upr: upper bound of confidence interval of annual secchi depth
- mean_drsif: mean annual dissolved reactive sillica 
- sd_drsif: standard deviation of annual filtered dissolved reactive sillica concentration
- ci_drsif_lwr: lower bound of confidence interval of annual filtered dissolved reactive sillica
- ci_drsif_upr: upper bound of confidence interval of annual filtered dissolved reactive sillica
- mean_chl: mean annual chlorophyll a concetration
- sd_chl: standard deviation of annual chlorophyll a concetration
- ci_chl_lwr: lower bound of confidence interval of annual chlorophyll a concetration
- ci_chl_upr: upper bound of confidence interval of annual chlorophyll a concetration
- mean_totnf: mean annual filtered total nitrogen concetration
- mean_no3no2: mean annual nitrate plus nitrite concetration
- mean_nh4: mean annual ammonium concetration
- mean_totpuf: mean annual total phosphorus unfiltered concentration
- Large: annual mean density of large zooplankton
- Small: annual mean density of small zooplankton
- Calanoid: annual mean density of *Calanoid* species
- Daphnia: annual mean density of *Daphnia* species
- CISCO: annual mean density of Cisco based on acoustic data
- LAKETROUT: annual mean density of lake trout based on acoustic data
- period: period assignment for the 3 major eras, zooplankton regime (1) piscivore regime (2), and the novel regime
- Group: species used for Mean_pred
- Mean_pred: mean annual abundance of predatory zooplanton, *Bythotrophes*
- sd_pred: standard deviation of annual abundance of predatory zooplanton, *Bythotrophes*
- year4: original year variable name from the NTL_LTER_TroutLake_Secchi.csv data file 
