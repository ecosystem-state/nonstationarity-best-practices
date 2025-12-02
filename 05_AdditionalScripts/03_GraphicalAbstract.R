library(MARSS)
library(mgcv)
library(dplyr)
library(forecast)
library(ggplot2)
library(lubridate)
library(viridis)
library(paletteer)
library(ggpubr)
#set.seed(1035)
col<-paletteer_d("nationalparkcolors::Everglades")
paletteer_d("nationalparkcolors::Everglades")
col2<-col[2:4]

set.seed(31)

col3<-col2
# defining the mean and varince of the different parameter distributions
sds=c(0.05,0.05)
means= c(-0.6,-0.1)
means2= c(1250,750)
means3=c(1000,1800)
lengths=c(20,30)
ran= c(0,300)
ransd=c(50,50)
periods=c(1,2)
parameter_stationary=c(1,1)
# simulating the distributions
ndist<-2 #number of distributions to simulate
simulated<-data.frame() # setting up empty dataframe to fill

#looping over means/sds/timeseries lengths for 3 distinct distributions
for(i in 1:ndist){
  parameter= arima.sim(n = lengths[i],list(ar=.2),rand.gen = rnorm,
                       sd = 0.1, mean=means[i])
  driver =  arima.sim(n = lengths[i],list(ar=.2),rand.gen = rnorm,
                                sd = 200, mean=means2[i])
  error =arima.sim(n = lengths[i],list(ar=.1),rand.gen = rnorm,
                      sd = ransd[i], mean=ran[i])
  response = driver*means[i]+error+means3[i]
  temp=data.frame(parameter=parameter,
                  period=as.factor(periods[i]),
                  mean=means[i],
                  sd=sds[i],
                  length=lengths[i],
                  parameter_stat=parameter_stationary[i],
                  driver = driver,
                  response)
  simulated=simulated%>%bind_rows(temp)
}

ystart<-1950

simulated<-simulated%>%
  mutate(Year=seq(ystart,ystart+sum(lengths)-1,1))%>%# adding a year column for plots
  mutate(period2=as.factor(ifelse(Year<=1969,1,2)))%>%
  mutate(driver_stand = scale(driver), response_stand=scale(response))

ggplot(data=simulated, aes(x=Year))+
  geom_line(aes(y=driver),color=col[3],lwd=1.25)+
  geom_line(aes(y=response),color=col[4],lwd=1.25)+
  ylab("Biological Condition \n (e.g. abundance)")+
  #ylim(c(0,3000))
  ggtitle("Change in Ecosystem")+
  geom_vline(xintercept=1968, lty=2)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))
ggsave("GraphicalAbstract_tsplot2.png", height =2.5, width = 3.5)
ggplot(data=simulated, aes(x=Year))+
  #geom_line(aes(y=mean),color=col[2],lwd=1.25)+
  geom_point(aes(y=parameter),color=col[1])+
  ylab("Parameter \n (e.g. slope)")+
  #ylim(c(0,1))+
  geom_vline(xintercept=1969, lty=2)+
  ggtitle("Observed Driver/Response \n Relationship")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))
ggsave("GraphicalAbstract_parameter2.png", height =2.5, width = 3.5)

sim1<-simulated%>%filter(period==1)%>%
  mutate(response=scale(response), driver=scale(driver))
sim2<-simulated%>%filter(period==2)%>%
  mutate(response=scale(response), driver=scale(driver))

sim<-bind_rows(sim1,sim2)

gam <- gam(response ~ Year+s(Year, by =driver),  data =sim)
gam_extract <- plot(gam, seWithMean = TRUE, n = nrow(sim))

coef_gam <- data.frame(date = simulated$Year, 
                               int_est = gam_extract[[1]]$fit,
                               int_se = gam_extract[[1]]$se,
                               slope_est = gam_extract[[1]]$fit,
                               slope_se = gam_extract[[1]]$se,
                       Model="Generalized Additive \n Model")

lm <- lm(response ~ period*driver, data =sim)
summ_lm<-summary(lm)
coef_lm<-bind_rows(coef_gam%>%filter(date<1970)%>%
  mutate(int_est=summ_lm$coefficients[1,1],
         int_se=summ_lm$coefficients[1,2],
         slope_est=summ_lm$coefficients[3,1],
         slope_se=summ_lm$coefficients[3,2]),
coef_gam%>%filter(date>1969)%>%
  mutate(int_est=summ_lm$coefficients[2,1]+summ_lm$coefficients[1,1],
         int_se=summ_lm$coefficients[2,2],
         slope_est=summ_lm$coefficients[4,1]+summ_lm$coefficients[3,1],
         slope_se=summ_lm$coefficients[4,2]))%>%
  mutate(Model="Linear Model")

dat <- sim
m <- 2
TT <- nrow(dat)
B <- diag(m)  ## 2x2; Identity
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now
diag(Q)[1] <- c("q.alpha")
diag(Q)[2] <- c("q.beta")
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- scale(dat$driver)  ## Nx1; predictor variable
A <- matrix(0)  ## 1x1; scalar = 0
R <- matrix("r")  ## 1x1; scalar = r
## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))
## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)
# convert response to matrix
dat_mat <- matrix(scale(dat$response), nrow = 1)
# fit the model -- crank up the maxit to ensure convergence
dlm_aksalmon_1 <- MARSS(dat_mat, inits = inits_list, model = mod_list,
                        control = list(maxit=4000), method="TMB")

# Fit a second model with just time-varying slope
diag(Q)[1] <- 1e-10
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)
dlm <- MARSS(dat_mat, inits = inits_list, model = mod_list,
                        control = list(maxit=4000), method="TMB")

coef_dlm <- data.frame(date = dat$Year,
                           int_est = 0,
                           int_se = 0,
                           slope_est = dlm$states[2,],
                           slope_se = dlm$states.se[2,],
                           #time_varying = "Slope",
                           Model = "Dynamic Linear \n Model")


coefs <- rbind(coef_gam, coef_dlm,coef_lm)
  

ggplot(data=coefs, aes(x=date, y=slope_est))+
geom_line(lwd=1, aes(linetype=Model),col=col[2])+
  geom_point(data=simulated, aes(x=Year, y=parameter-0.1), col=col[1])+
  ylim(c(-2, 1))+
  ylab("Parameter (e.g. slope)")+
  xlab("Year")+
  ggtitle("Model Estimated \n Driver/Response Relationship")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  geom_vline(xintercept=1969, lty=2)+
  theme(plot.title = element_text(hjust = 0.5))
ggsave("GraphicalAbstract_modeModel2.png", height =3, width = 6)

ggplot()+
  geom_point(data=simulated%>%filter(period==2),aes(y=response, x=driver),color=col3[2])+
  geom_point(data=simulated%>%filter(period==1),aes(y=response, x=driver),color=col3[1])+
  geom_smooth(data=simulated%>%filter(period==2),aes(y=response, x=driver),method='lm',color=col3[2])+
  geom_smooth(data=simulated%>%filter(period==1),aes(y=response, x=driver),method='lm',color=col3[1])+
  
   # ylab("Parameter (e.g. slope)")+
  #ylim(c(0,1))+
  #geom_vline(xintercept=1969, lty=2)+
  theme_classic()
  
cor(simulated%>%filter(period==2)%>%select(response),simulated%>%filter(period==2)%>%select(driver))
cor(simulated%>%filter(period==1)%>%select(response),simulated%>%filter(period==1)%>%select(driver))
