
library(MARSS)
library(mgcv)
library(dplyr)
library(forecast)
library(ggplot2)
library(lubridate)
library(marssTMB)
library(tidyverse)
library(strucchange)
library(data.table)


## Salmon Data ##

# Load salmon timeseries collapsed into DFA trend
salm.dat <- readRDS(here::here("03_GOA_salmon_case_study/Data/dfa_trend.rds")) 

# Load 3-year running mean winter SST timeseries
sst.dat <- readRDS(here::here("03_GOA_salmon_case_study/Data/winterSST_3yr_running_mean.rds")) 

# Bind
right_join(salm.dat, sst.dat) -> dat


# Scale predictor sst
dat %>%
  mutate(z_sst = scale(sst_3yr_running_mean)) %>%
  na.omit() -> dat2

Fit_slm<-ggplot(data=dat2, aes(group=period,col=period, y=salmon_DFA,x=z_sst)) +
  stat_smooth(method = "lm")+
  scale_colour_discrete(labels=c("before 1988/89", "after 1988/89"))+
  geom_point()+
  ylab("Salmon trend")+
  xlab("SST")

## LTER Data ##

full_data <- readRDS("01_LTER_NTL_case_study/NTL_LTER_TL_data.rds")
analysis_data<-full_data%>%
  select(year4, mean_sec,mean_chl, Large,mean_totpuf,period)%>%
  mutate(mean_sec_scale=scale(mean_sec),mean_chl_scale=scale(mean_chl), 
         Large_scale=scale(Large),mean_totpuf_scale=scale(mean_totpuf))
analysis_long<-analysis_data%>%
  pivot_longer(!c(year4,period), names_to = 'Dataset', values_to='values')
dataset_names <- list(
  'mean_sec_scale'="Water Clarity",
  'mean_totpuf_scale'="Total Phosphorus",
  'mean_chl_scale'="Chlorophyll",
  'Large_scale'="Large Zooplankton"
)
facet_labeller <- function(variable,value){
  return(dataset_names[value])
}

## Lake Washington ##
dat_lakewa <- readRDS("04_DLM - GAM comparison/data_for_lakeWA_example.rds")

# renaming
dat_lakewa <- dplyr::rename(dat_lakewa, 
                            response = y_adj,
                            driver = lagged_zp)%>%
  select(Year, Month, Bluegreens, TP, date)%>%
  pivot_longer(cols=-c(date, Year, Month),names_to = "TimeSeries", values_to = "Values")%>%
  group_by(TimeSeries)%>%
  mutate(Anomalies=scale(Values))



### Plots ###
annotateWA<-data.frame(label=c("Sewage Inputs", "Water Treatment"),
                       date=c(dat_lakewa$date[50],dat_lakewa$date[200]),
                       y=c(5,5), TimeSeries=c("TP","TP"))
lakeWA<- ggplot(data=dat_lakewa, aes(y=Anomalies, x=date, group=TimeSeries,linetype=TimeSeries, col=TimeSeries)) +
  geom_hline(yintercept=0)+
  geom_line(size = 1) +
  #ylim(c(-2,5))+
 
  labs(color = "Time Series")+
  scale_colour_manual(values=c(col[2],col[3]),name="Time Series",labels=c("Cyanobacteria \n (Biological Response)", "Total Phosphorus \n (Ecological Driver)"))+
  scale_linetype_manual(values=c("solid","dashed"),name="Time Series",labels=c("Cyanobacteria \n (Biological Response)", "Total Phosphorus \n (Ecological Driver)"))+
  theme_classic()+
  ylab("Scaled Anomaly")+
  geom_vline(xintercept=as.numeric(dat_lakewa$date[111]), lty=2,col='grey')+
  geom_text(data=annotateWA, aes(x = date, y = y, 
                                label = label),  col='black',fontface = "italic")+
  xlab("")
lakeWA


annotateGoA<-data.frame(label=c("High Aleutian \n Low Variability",
                                "Low Aleutian \n Low Variability",
                                "Marine \n Heatwaves"),
                       date=c(1977, 2000,2020),
                       y=c(3.25,3.25,3.25), Anomalies=c("z_sst","z_sst","z_sst"))
GoA_TS<- ggplot(data=dat2%>%
                  select(year, salmon_DFA, z_sst)%>%
                  mutate(salmon=scale(salmon_DFA))%>%
                  pivot_longer(-c(year,salmon_DFA),values_to = 'values', names_to = 'Anomalies'),
                aes(y=values,x=year,col=Anomalies, linetype=Anomalies))+
  geom_line(size = 1) +
  geom_point() +
 # labs(color = "Time Series")+
  scale_linetype_manual(values=c("solid","dashed"), name="Time Series",labels=c("Salmon\n (Biological Response)", "SST\n (Ecological Driver)"))+
  scale_colour_manual(values=c(col[2],col[3]), name="Time Series",labels=c("Salmon\n (Biological Response)", "SST\n (Ecological Driver)"))+
  geom_vline(xintercept=c(1988,2014), lty=c(2,2),col='grey')+
  theme_classic()+
  geom_text(data=annotateGoA, aes(x = date, y = y, 
                                 label = label),  
            col='black',fontface = "italic")+
  geom_hline(yintercept=0)+
  ylab("")+
  ylim(c(-2.5,3.5))+
  xlim(c(1970,2024))+
  xlab("Year")

GoA_TS

annotateLTER<-data.frame(label=c("Zooplankton \n regime",
                                "Piscivore \n regime",
                                "Novel \n regime"),
                        date=c(1992, 2011,2020),
                        y=c(3.35,3.35,3.35), Dataset=c("Large_scale","Large_scale","Large_scale"))
LTER<-ggplot(data=analysis_long%>%
         filter(Dataset=="mean_sec_scale"| Dataset=="Large_scale"), aes(y=values, x=year4, group=Dataset,linetype=Dataset, col=Dataset)) +
  geom_hline(yintercept=0)+
  geom_line(size = 1) +
  ylim(c(-2.5,3.7))+
  geom_point() +
  labs(color = "Time Series")+
  scale_colour_manual(values=c(col[2],col[3]),name="Time Series",labels=c("Large Zooplankton \n (Biological Response)", "Water Clarity \n (Ecological Driver)"))+
  scale_linetype_manual(values=c("solid","dashed"),name="Time Series",labels=c("Large Zooplankton \n (Biological Response)", "Water Clarity \n (Ecological Driver)"))+
   theme_classic()+
  geom_text(data=annotateLTER, aes(x = date, y = y, 
                                  label = label),  
            col='black',fontface = "italic")+
  ylab("")+
  geom_vline(xintercept=2007, lty=2,col='grey')+
  geom_vline(xintercept=2014.5, lty=2,col='grey')+
  xlab("")
LTER

ggarrange(LTER, lakeWA, GoA_TS, labels=c("A.", "B.", "C."), nrow=3)
ggsave("Figures/Figure_2_timeseries.png", height = 8, width = 7)
