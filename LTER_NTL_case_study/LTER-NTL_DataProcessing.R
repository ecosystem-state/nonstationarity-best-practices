library(ggpubr)
library(nord)
library(ggplot2)
library(dplyr)
library(stringr)
library(dplyr)
library(ggplot2)
library(nord)
library(tidyr)
library(data.table)
library(mgcv)

##### Functions #####
cis <- NA
ci_func <- function(data,column) {
  for(i in 1:length(years)){
    temp<-data%>%filter(year4==years[i])
    ci_annual_lwr=t.test(temp$column)$conf.int[1]
    ci_annual_upr=t.test(temp$column)$conf.int[2]
    cis_temp<-cbind(year=years[i],ci_annual_lwr,ci_annual_upr)
    cis<-rbind(cis,cis_temp)
  }
}

###### Secchi data #####
secchi <- read.csv("LTER_NTL_case_study/NTL-LTER-TroutLake/NTL_LTER_TroutLake_Secchi.csv")
secchi_TR<-secchi%>%
  filter(lakeid=="TR")%>%
  select(year4, daynum, sampledate, sta, secview, secnview)%>% 
  filter_at(vars(secview, secnview),any_vars(!is.na(.))) #removing where there are no observations (both NAs))
secchi_TR<-secchi_TR%>%
  mutate(secviewMean=(secview+secnview)/2)%>%
  mutate(secviewMean = coalesce(secviewMean, secview))%>% 
  mutate(secviewMean = coalesce(secviewMean, secnview))

###### removing seasonality ######
gam_seas<-gam(secviewMean~s(daynum, bs = "cc") + s(year4),
       data = secchi_TR)
years<-unique(secchi_TR$year4)
newdata<-data.frame(year4=years, daynum = rep(0, length(years)))
pred_secchi<-predict(gam_seas, newdata = newdata, interval = "confidence", se.fit = TRUE)
pred_secchi_dat<-data.frame(cbind(year4=years,fit=as.numeric(pred_secchi$fit),se.fit=as.numeric(pred_secchi$se.fit)))%>%
  mutate(ci_upr=fit + (2 * se.fit), ci_lwr=fit - (2 * se.fit))
ggplot(data = pred_secchi_dat, aes(x = year4, y = fit)) +
  #  geom_ribbon(alpha=0.2,linetype = 0,aes(ymin=ln.sr-sr.CV*ln.sr, ymax=ln.sr+sr.CV*ln.sr,lwd=0))+
  geom_point(size=0.75)+
  geom_errorbar(aes(ymin=ci_lwr,ymax=ci_upr))+
  geom_line()

###### plotting raw data ######
cis <- NA
for(i in 1:length(years)){
  secchi_temp<-secchi_TR%>%filter(year4==years[i])
  ci_annual_lwr=t.test(secchi_temp$secviewMean)$conf.int[1]
  ci_annual_upr=t.test(secchi_temp$secviewMean)$conf.int[2]
  cis_temp<-cbind(year=years[i],ci_annual_lwr,ci_annual_upr)
  cis<-rbind(cis,cis_temp)
}

secchi_TR_summary<-secchi_TR%>%  
  group_by(year4)%>%
  summarise(mean_annual=mean(secviewMean), sd_annual=sd(secviewMean))%>%
  cbind(na.omit(cis))
secchi.plot <- ggplot(data =secchi_TR_summary, aes(x = year4, y = mean_annual)) +
  #  geom_ribbon(alpha=0.2,linetype = 0,aes(ymin=ln.sr-sr.CV*ln.sr, ymax=ln.sr+sr.CV*ln.sr,lwd=0))+
  geom_point(size=0.75)+
  ylim(c(3,7.5))+
  ggtitle("Water Clarity")+
  ylab("Secchi depth (m)")+
  xlab("Year")+
  geom_smooth(span=0.2, col='black', lwd=.25)+
  geom_errorbar(aes(ymin=ci_annual_lwr,ymax=ci_annual_upr))
secchi.plot

##### Dissolved Reactive Silica #####
chemistry <- read.csv("LTER_NTL_case_study/NTL-LTER-TroutLake/NTL_LTER_TroutLake_Chemistry.csv")
chemistry_TR<-chemistry%>%
  filter(lakeid=="TR", depth==0.0)%>%
  select(year4, daynum, sampledate, depth, rep, sta,drsif)%>% 
  filter_at(vars(drsif),all_vars(!is.na(.))) #removing where there are no observations (both NAs))
###### plotting raw data ######
cis <- NA
years<-unique(chemistry_TR$year4)
for(i in 2:length(years)){
  chemistry_temp<-chemistry_TR%>%filter(year4==years[i])
  ci_annual_lwr=t.test(chemistry_temp$drsif)$conf.int[1]
  ci_annual_upr=t.test(chemistry_temp$drsif)$conf.int[2]
  cis_temp<-cbind(year=years[i],ci_annual_lwr,ci_annual_upr)
  cis<-rbind(cis,cis_temp)
}

chemistry_TR_summary<-chemistry_TR%>%  
  group_by(year4)%>%
  summarise(mean_annual=mean(drsif), sd_annual=sd(drsif))%>%
  cbind(cis)
drsf.plot <- ggplot(data =chemistry_TR_summary, aes(x = year4, y = mean_annual)) +
  #  geom_ribbon(alpha=0.2,linetype = 0,aes(ymin=ln.sr-sr.CV*ln.sr, ymax=ln.sr+sr.CV*ln.sr,lwd=0))+
  geom_point(size=0.75)+
  #ylim(c(3,7.5))+
  ggtitle("Silica (Datiom Indicator)")+
  ylab(expression(paste('Dissolved Reactive Silica (filtered) ('~mu, "g L"^"-1",")")))+
  xlab("Year")+
  geom_smooth(span=0.2, col='black', lwd=.25)+
  geom_errorbar(aes(ymin=ci_annual_lwr,ymax=ci_annual_upr))
drsf.plot


##### Chlorophyll #####