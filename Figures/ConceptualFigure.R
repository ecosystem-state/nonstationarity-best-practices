library(MARSS)
library(mgcv)
library(dplyr)
library(forecast)
library(ggplot2)
library(lubridate)
library(viridis)
library(paletteer)
library(ggpubr)
set.seed(1035)
set.seed(9450)
col<-paletteer_d("nationalparkcolors::Everglades")
col2<-col[2:4]
# defining the mean and varince of the different parameter distributions
sds=c(0.1,0.1,0.35)
means= c(0.2,-0.15,-0.15)
lengths=c(17,20,25)
periods=c(1,2,3)
parameter_stationary=c(1,1,1)
# simulating the distributions
ndist<-3 #number of distributions to simulate
simulated<-data.frame() # setting up empty dataframe to fill

#looping over means/sds/timeseries lengths for 3 distinct distributions
for(i in 1:ndist){
  parameter=arima.sim(n = lengths[i],list(ar=.5),rand.gen = rnorm,
            sd = sds[i], mean=means[i])
  temp=data.frame(parameter=parameter,
             period=as.factor(periods[i]),
             mean=means[i],
             sd=sds[i],
             length=lengths[i],
             parameter_stat=parameter_stationary[i])
  simulated=simulated%>%bind_rows(temp)
}

#Adding additional columns for plotting
ystart<-1950

simulated<-simulated%>%
  mutate(Year=seq(ystart,ystart+sum(lengths)-1,1))%>%# adding a year column for plots
  mutate(xdata=rnorm(sum(lengths),0,0.2))%>%
  mutate(data1 = parameter*xdata+parameter_stat)%>%
  mutate(data2 =parameter_stat*xdata+parameter)%>%
  mutate(data3 =parameter*xdata+parameter)%>%
  mutate(sim=arima.sim(n = sum(lengths),list(ar=.5),rand.gen = rnorm,
                       sd =mean(sds), mean=mean(means)))

simulated_results<-simulated%>%
  group_by(period)%>%
  summarise(mean=mean(parameter),sd=sd(parameter))%>%
  mutate(ymax=mean+sd,ymin=mean-sd,parameter_stationary=parameter_stationary)%>%
  mutate(xmin=c(1950, 1967,1987), xmax=c(1966,1988,2012))

simulated_results2<-simulated%>%
  group_by(period)%>%
  summarise(mean=mean(sim),sd=sd(sim))%>%
  mutate(ymax=mean+sd,ymin=mean-sd,parameter_stationary=parameter_stationary)%>%
  mutate(xmin=c(1950, 1967,1987), xmax=c(1966,1988,2012))

aplot<-ggplot(data=simulated,aes(x=Year,y=parameter,color=period,fill=period)) +
  geom_point(aes(group=period))+
  geom_line(aes(group=period))+
  geom_segment(data=simulated_results,aes(x = xmin, y = mean,
        xend = xmax, yend = mean, group=period,color=period),lwd=1)+ 
  geom_segment(data=simulated_results,aes(x = xmin, y = ymax,
        xend = xmax, yend = ymax, group=period,color=period),lwd=1,lty=2)+
  geom_segment(data=simulated_results,aes(x = xmin, y = ymin,
                                          xend = xmax, yend = ymin, group=period,color=period),lwd=1,lty=2)+
  scale_colour_manual(values=col2,labels = c("1950 - 1966", "1966 - 1986", "1986 - 2011"))+  
  scale_fill_manual(values=col2,labels = c("1950 - 1966", "1966 - 1986", "1986 - 2011"))+
  theme_bw()+ 
  ylab("Parameter")+
  theme(legend.position="none")
aplot

bplot<-ggplot(data=simulated, aes(x=parameter,group=period,fill=period)) +
  geom_density(alpha=0.8)+
  scale_fill_manual(values=col2,labels = c("1950 - 1966", "1966 - 1986", "1986 - 2011"))+
  xlim(c(-1,1))+
  xlab("Parameter")+
  theme_bw()+ 
  theme(legend.position="none")
bplot


cplot<-ggplot(data=simulated_results,) +
  geom_point(data=simulated, aes(x=xdata,y=data1,group = period,color=period))+
  geom_abline(aes(intercept = parameter_stationary, slope = mean, color=factor(period)),lwd=1)+
 # geom_smooth(data=simulated, aes(x=xdata,y=data1,group = period,color=period,fill=period),method='lm',fullrange=T,alpha=0.2)+
  scale_colour_manual(values=col2,labels = c("1950 - 1966", "1966 - 1986", "1986 - 2011"))+
  ggtitle("Time-Varying Slope")+
  scale_fill_manual(values=col2,labels = c("1950 - 1966", "1966 - 1986", "1986 - 2011"))+
 xlab("")+
  ylab("Response")+
  #geom_ribbon(aes(y=mean,ymax=ymax,ymin=ymin))
  theme_bw()+
 # ylim(c(-1,1))+
  theme(plot.title = element_text(hjust = 0.5))+ 
  theme(legend.position="none")
cplot

dplot<-ggplot(data=simulated_results,) +
  geom_point(data=simulated, aes(x=xdata,y=data2,group = period,color=period))+
 geom_abline(aes(intercept = mean, slope = parameter_stationary, color=factor(period)),lwd=1)+
  #geom_ribbon(aes(group = period,y=mean,ymax=ymax,ymin=ymin, xmin=0.5, xmax=1.5), fill='grey')+
 #geom_smooth(data=simulated, aes(x=xdata,y=data2,group = period,color=period,fill=period),method='lm',fullrange=T,alpha=0.2)+
  scale_colour_manual(values=col2,labels = c("1950 - 1966", "1966 - 1986", "1986 - 2011"))+
  ggtitle("Time-Varying Intercept")+
  scale_fill_manual(values=col2,labels = c("1950 - 1966", "1966 - 1986", "1986 - 2011"))+
  xlab("Driver")+
  ylab("")+
 # ylim(c(-1,1))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+ 
  theme(legend.position="none")
dplot

eplot<-ggplot(data=simulated_results) +
  geom_point(data=simulated, aes(x=xdata,y=data3,group = period,color=period))+
 # geom_smooth(data=simulated, aes(x=xdata,y=data3,group = period,color=period),method='lm')+
  geom_abline(aes(intercept = mean, slope = mean, color=factor(period)),lwd=1)+
 # geom_smooth(data=simulated, aes(x=xdata,y=data3,group = period,color=period,fill=period),method='lm',fullrange=T,alpha=0.2)+
  scale_colour_manual(values=col2,labels = c("1950 - 1966", "1966 - 1986", "1986 - 2011"))+
  scale_fill_manual(values=col2,labels = c("1950 - 1966", "1966 - 1986", "1986 - 2011"))+
  xlab("")+
  ylab("")+
  ggtitle("Time-Varying Slope & Intercept")+
  theme_bw()+ 
 # ylim(c(-1,1))+
 # theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5))
eplot

fplot<-ggplot(data=simulated,aes(x=Year,y=sim,color=period,fill=period)) +
  geom_point(aes(group=period))+
  geom_line(aes(group=period))+
  geom_segment(data=simulated_results2,aes(x = xmin, y = mean,
                                          xend = xmax, yend = mean, group=period,color=period),lwd=1)+ 
  geom_segment(data=simulated_results2,aes(x = xmin, y = ymax,
                                          xend = xmax, yend = ymax, group=period,color=period),lwd=1,lty=2)+
  geom_segment(data=simulated_results2,aes(x = xmin, y = ymin,
                                          xend = xmax, yend = ymin, group=period,color=period),lwd=1,lty=2)+
  
  scale_colour_manual(values=col2,labels = c("1950 - 1966", "1966 - 1986", "1986 - 2011"))+
  scale_fill_manual(values=col2,labels = c("1950 - 1966", "1966 - 1986", "1986 - 2011"))+
  #geom_smooth(aes(group=period),method='lm',alpha=0.2)+
  ylab("Parameter")+
  theme_bw()+ 
  theme(legend.position="none")
fplot

 gplot<-ggplot(data=simulated, aes(x=sim,group=period,fill=period)) +
  geom_density(alpha=0.8)+
  xlim(c(-1,1))+
  scale_fill_manual(values=col2,labels = c("1950 - 1966", "1966 - 1986", "1986 - 2011"))+
  theme_bw()+ 
  xlab("Parameter")+
  theme(legend.position="none")
gplot

row1<-ggarrange(aplot,bplot, labels=c("A.", "B."), nrow=1,ncol=2)
row1a<-annotate_figure(row1, top = text_grob("Non-stationary Parameter",face = "bold", size = 14)
                       )
row2<-ggarrange(cplot,dplot,eplot, labels=c("C.", "D.","E."), nrow=1,ncol=3,widths = c(1,1,1.4))
row3<-ggarrange(fplot,gplot, labels=c("F.", "G."), nrow=1,ncol=2)
row3a<-annotate_figure(row3, top = text_grob("Stationary Parameter",face = "bold", size = 14))

ggarrange(row1,row2,row3,nrow=3,ncol=1)
pdf(file = "Figure1_Theorhetical.pdf",   # The directory you want to save the file in
        width = 9, # The width of the plot in inches
       height = 9)
ggarrange(row1a,row2,row3a,nrow=3,ncol=1)
dev.off()
