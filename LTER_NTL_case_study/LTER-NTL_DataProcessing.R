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
    ci_sec_lwr=t.test(temp$column)$conf.int[1]
    ci_sec_upr=t.test(temp$column)$conf.int[2]
    cis_temp<-cbind(year=years[i],ci_sec_lwr,ci_sec_upr)
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
  #geom_errorbar(aes(ymin=ci_lwr,ymax=ci_upr))+
  geom_line()

###### plotting raw data ######
cis <- NA
for(i in 1:length(years)){
  secchi_temp<-secchi_TR%>%filter(year4==years[i])
  ci_sec_lwr=t.test(secchi_temp$secviewMean)$conf.int[1]
  ci_sec_upr=t.test(secchi_temp$secviewMean)$conf.int[2]
  cis_temp<-cbind(year=years[i],ci_sec_lwr,ci_sec_upr)
  cis<-rbind(cis,cis_temp)
}

secchi_TR_summary<-secchi_TR%>%  
  group_by(year4)%>%
  summarise(mean_sec=mean(secviewMean), sd_sec=sd(secviewMean))%>%
  cbind(na.omit(cis))
secchi.plot <- ggplot(data =secchi_TR_summary, aes(x = year4, y = mean_sec)) +
  #  geom_ribbon(alpha=0.2,linetype = 0,aes(ymin=ln.sr-sr.CV*ln.sr, ymax=ln.sr+sr.CV*ln.sr,lwd=0))+
  geom_point(size=0.75)+
  ylim(c(3,7.5))+
  ggtitle("Water Clarity")+
  ylab("Secchi depth (m)")+
  xlab("Year")+
  geom_smooth(span=0.2, col='black', lwd=.25)+
  #geom_errorbar(aes(ymin=ci_sec_lwr,ymax=ci_sec_upr))+
  geom_vline(xintercept=2007, lty=2)+
  geom_vline(xintercept=2014, lty=2)
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
  ci_drsif_lwr=t.test(chemistry_temp$drsif)$conf.int[1]
  ci_drsif_upr=t.test(chemistry_temp$drsif)$conf.int[2]
  cis_temp<-cbind(year=years[i],ci_drsif_lwr,ci_drsif_upr)
  cis<-rbind(cis,cis_temp)
}

chemistry_TR_summary<-chemistry_TR%>%  
  group_by(year4)%>%
  summarise(mean_drsif=mean(drsif), sd_drsif=sd(drsif))%>%
  cbind(cis)
drsf.plot <- ggplot(data =chemistry_TR_summary, aes(x = year4, y = mean_drsif)) +
  #  geom_ribbon(alpha=0.2,linetype = 0,aes(ymin=ln.sr-sr.CV*ln.sr, ymax=ln.sr+sr.CV*ln.sr,lwd=0))+
  geom_point(size=0.75)+
  #ylim(c(3,7.5))+
  ggtitle("Silica (Datiom Indicator)")+
  ylab(expression(paste('Dissolved Reactive Silica (filtered) ('~mu, "g L"^"-1",")")))+
  xlab("Year")+
  geom_smooth(span=0.2, col='black', lwd=.25)+
 # geom_errorbar(aes(ymin=ci_drsif_lwr,ymax=ci_drsif_upr))+
  geom_vline(xintercept=2007, lty=2)+
  geom_vline(xintercept=2014, lty=2)
drsf.plot


##### Chlorophyll #####
chl <- read.csv("LTER_NTL_case_study/NTL-LTER-TroutLake/NTL_LTER_TroutLake_ChlA.csv")
chl_TR<-chl%>%
  filter(lakeid=="TR", depth==0.0)%>%
  select(year4, daynum, sampledate, depth, rep, sta,chlor)%>% 
  filter_at(vars(chlor),all_vars(!is.na(.))) #removing where there are no observations (both NAs))
###### plotting raw data ######
cis <- NA
years<-unique(chl_TR$year4)
for(i in 1:length(years)){
  chl_temp<-chl_TR%>%filter(year4==years[i])
  ci_chl_lwr=t.test(chl_temp$chlor)$conf.int[1]
  ci_chl_upr=t.test(chl_temp$chlor)$conf.int[2]
  cis_temp<-cbind(year=years[i],ci_chl_lwr,ci_chl_upr)
  cis<-rbind(cis,cis_temp)
}

chl_TR_summary<-chl_TR%>%  
  group_by(year4)%>%
  summarise(mean_chl=mean(chlor), sd_chl=sd(chlor))%>%
  cbind(na.omit(cis))
chl.plot <- ggplot(data =chl_TR_summary, aes(x = year4, y = mean_chl)) +
  #  geom_ribbon(alpha=0.2,linetype = 0,aes(ymin=ln.sr-sr.CV*ln.sr, ymax=ln.sr+sr.CV*ln.sr,lwd=0))+
  geom_point(size=0.75)+
  #ylim(c(3,7.5))+
  ggtitle("Chlorophyll")+
  ylab(expression(paste('Chlorophyll'~alpha," concentrations ("~mu, "g L"^"-1",")")))+
  xlab("Year")+
  geom_smooth(span=0.2, col='black', lwd=.25)+
 # geom_errorbar(aes(ymin=ci_chl_lwr,ymax=ci_chl_upr))+
  geom_vline(xintercept=2007, lty=2)+
  geom_vline(xintercept=2014, lty=2)
chl.plot


##### Predatory Zooplankton #####
pred <- read.csv("LTER_NTL_case_study/NTL-LTER-TroutLake/NTL_LTER_PredZoop.csv")
colnames(pred)
unique(pred$taxon)
pred_TR<-pred%>%
  filter(lakeid=="TR")%>%
  select(year4, sta, depth, avg_ind_per_m3,stdev_ind_per_m3,taxon)%>%
  mutate(Group=ifelse(taxon=="CHAOBORUS LARVAE"|taxon=="CHAOBORUS PUPAE","CHAOBORUS",taxon))

pred_TR_summary<-na.omit(pred_TR)%>%
  group_by(year4,Group)%>%
  summarise(mean_pred=mean(avg_ind_per_m3), sd_pred=mean(stdev_ind_per_m3))

pred.plot<- ggplot(data =pred_TR_summary, aes(x = year4, y = mean_pred, group=Group)) +
  #  geom_ribbon(alpha=0.2,linetype = 0,aes(ymin=ln.sr-sr.CV*ln.sr, ymax=ln.sr+sr.CV*ln.sr,lwd=0))+
  #facet_grid(rows = vars(Type),scales='free')+
  geom_point(size=2.5, aes(pch=Group))+
  #ylim(c(3,7.5))+
   geom_line()+
  ggtitle("Predatory Zooplankton")+
  ylab(expression(paste('Density (Individuals m' ^"-3",")")))+
  xlab("Year")+
   geom_vline(xintercept=2007, lty=2)+
   geom_vline(xintercept=2014, lty=2)
 # geom_smooth(span=0.2, col='black', lwd=.25)+
 # geom_errorbar(aes(ymin=ci_chl_lwr,ymax=ci_chl_upr))


##### Predatory Zooplankton #####
pred <- read.csv("LTER_NTL_case_study/NTL-LTER-TroutLake/NTL_LTER_PredZoop.csv")
colnames(pred)
unique(pred$taxon)
pred_TR<-pred%>%
  filter(lakeid=="TR")%>%
  select(year4, sta, depth, avg_ind_per_m3,stdev_ind_per_m3,taxon)%>%
  mutate(Group=ifelse(taxon=="CHAOBORUS LARVAE"|taxon=="CHAOBORUS PUPAE","CHAOBORUS",taxon))

pred_TR_summary<-na.omit(pred_TR)%>%
  group_by(year4,Group)%>%
  summarise(mean_pred=mean(avg_ind_per_m3), sd_pred=mean(stdev_ind_per_m3))

pred.plot<- ggplot(data =pred_TR_summary, aes(x = year4, y = mean_pred, group=Group)) +
  #  geom_ribbon(alpha=0.2,linetype = 0,aes(ymin=ln.sr-sr.CV*ln.sr, ymax=ln.sr+sr.CV*ln.sr,lwd=0))+
  #facet_grid(rows = vars(Type),scales='free')+
  geom_point(size=2.5, aes(pch=Group))+
  #ylim(c(3,7.5))+
   geom_line()+
  ggtitle("Predatory Zooplankton")+
  ylab(expression(paste('Density (Individuals m' ^"-3",")")))+
  xlab("Year")+
   geom_vline(xintercept=2007, lty=2)+
   geom_vline(xintercept=2014, lty=2)
 # geom_smooth(span=0.2, col='black', lwd=.25)+
 # geom_errorbar(aes(ymin=ci_chl_lwr,ymax=ci_chl_upr))

##### Zooplankton #####
zoop <- read.csv("LTER_NTL_case_study/NTL-LTER-TroutLake/NTL_LTER_TroutLake_Zooplankton.csv")
zoopsCODE <- read.csv("LTER_NTL_case_study/NTL-LTER-TroutLake/zoopsCODE.csv")

zoop_TR<-zoop%>%
  filter(lakeid=="TR")%>%
  select(year4, sample_date, species_code, species_name, density)%>% #removing where there are no observations (both NAs))
  left_join(zoopsCODE)

zoop_TR_summary<-na.omit(zoop_TR)%>%
  group_by(year4,Group, Type)%>%
  summarise(mean_zoop=mean(density), sd_zoop=sd(density))

zoop.plot<- ggplot(data =zoop_TR_summary, aes(x = year4, y = mean_zoop, group=Group)) +
  #  geom_ribbon(alpha=0.2,linetype = 0,aes(ymin=ln.sr-sr.CV*ln.sr, ymax=ln.sr+sr.CV*ln.sr,lwd=0))+
  facet_grid(rows = vars(Type),scales='free')+
  geom_point(size=2.5, aes(pch=Group))+
  #ylim(c(3,7.5))+
  geom_line()+
  ggtitle("Zooplankton taxa")+
  ylab(expression(paste('Density (Individuals L' ^"-1",")")))+
  xlab("Year")+
  geom_vline(xintercept=2007, lty=2)+
  geom_vline(xintercept=2014, lty=2)
# geom_smooth(span=0.2, col='black', lwd=.25)+
# geom_errorbar(aes(ymin=ci_chl_lwr,ymax=ci_chl_upr))



##### Cisco and Lake Trout #####
acoustic <- read.csv("LTER_NTL_case_study/NTL-LTER-TroutLake/NTL_LTER_TroutLake_FishDensity.csv")
CPUE <- read.csv("LTER_NTL_case_study/NTL-LTER-TroutLake/NTL_LTER_TroutLake_CPUE.csv")
Lengths <- read.csv("LTER_NTL_case_study/NTL-LTER-TroutLake/NTL_LTER_Lengths.csv")

acoustic_TR<-acoustic%>%
  filter(lakeid=="TR"&species=="CISCO"|species=="LAKETROUT")%>%
  select(year4, sampledate,species, density)
acoustic_TR_summary<-na.omit(acoustic_TR)%>%
  group_by(year4,species)%>%
  summarise(mean_acoustic=mean(density))

acoustic_TR_summary2<-left_join(data.frame(year4 = rep(c(2001:2022), 2),
          species=rep(c("CISCO","LAKETROUT"), each=22)), 
          acoustic_TR_summary)%>%
  mutate(hectare=mean_acoustic/343)
acoustic_TR_summary2[is.na(acoustic_TR_summary2)] <- 0
acoustic.plot<- ggplot(data =acoustic_TR_summary2, aes(x = year4, y = hectare, group=species)) +
  facet_grid(rows = vars(species),scales='free')+
  geom_bar(stat = "identity")+
  ggtitle("Hydroacoustic Density")+
  ylab(expression(paste('Population Density (# per hectare)')))+
  xlab("Year")+
  geom_vline(xintercept=2007, lty=2)+
  geom_vline(xintercept=2014, lty=2)

TR_Cisco_Lengths <-Lengths%>%
  filter(lakeid=="TR"&spname=="CISCO")%>%
  filter(gearid=="VGN019"|gearid== "VGN032"|gearid==  "VGN038"|
           gearid== "VGN051" |gearid== "VGN089" |gearid== "VGN025"|
           gearid==  "VGN064")%>%
  filter(length>250)%>%
  select(year4, length)%>%
  group_by(year4)%>%
  summarise(n=n())

CPUE_TR<-CPUE%>%
  filter(lakeid=="TR"&spname=="CISCO")%>%
  filter(gearid=="VGN019"|gearid== "VGN032"|gearid==  "VGN038"|
           gearid== "VGN051" |gearid== "VGN089" |gearid== "VGN025"|
           gearid==  "VGN064")%>%
  group_by(year4)%>%
  summarise(effort2=max(effort))%>% 
  left_join(TR_Cisco_Lengths)%>%
  mutate(CPUE=n/effort2)
CPUE_TR[is.na(CPUE_TR)] <- 0

CPUE.plot<-ggplot(data =CPUE_TR, aes(x = year4, y = CPUE)) +
  geom_bar(stat = "identity")+
  ggtitle("Cisco > 250 mm")+
  ylab(expression(paste('Catch per Unit Effort')))+
  xlab("Year")+
  geom_vline(xintercept=2007, lty=2)+
  geom_vline(xintercept=2014, lty=2)



##### writing full plots and datasets ####
pdf(file = "LTER_NTL_case_study/NTL_LTER_TSplots.pdf",  
    width = 6,
    height = 8) 

ggarrange(secchi.plot, drsf.plot,
          chl.plot, nrow = 3)
ggarrange(zoop.plot, pred.plot, nrow = 2)
ggarrange( acoustic.plot,
           CPUE.plot, nrow = 2)
dev.off()


