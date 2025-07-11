library(tidyverse)
library(reshape2)
library(MARSS)
library(ggplot2)
library(viridis)


dat=read.csv("data/wa_oyster_condition.csv")
names(dat)=tolower(names(dat))
dat.test=arrange(dat,year)
dfa.dat=dat[which(dat$year>=1955&dat$year<=2023),] #same length as covariate data

melt.dat=melt(dfa.dat, id.vars = c("year","site"))
names(melt.dat)[2:3]=c("code","month")


####### Choose seasonal period to average oyster condition index (CI)

####Annual spring-summer means based May to Sept
subyear.dat <- melt.dat[melt.dat$month %in% c("may","jun","jul","aug","sep"),]   
ann.dat=subyear.dat %>%
  group_by(year,code) %>%
  summarize(mean.val= mean(value, na.rm=T)
  )      
names(ann.dat)=c("year","code","mean.val")


##Annual winter means based Oct to April data
subyear.dat <- melt.dat[melt.dat$month %in% c("oct", "nov","dec","jan", "feb", "mar", "apr"),]
subyear.dat$win.yr <- ifelse(subyear.dat$month %in% c("oct", "nov","dec"), subyear.dat$year+1, subyear.dat$year)
ann.dat=subyear.dat %>%
  group_by(win.yr,code) %>%
  summarize(mean.val= mean(value, na.rm=T)
  )
names(ann.dat)=c("year","code","mean.val")
ann.dat=ann.dat[which(ann.dat$year>=1955&ann.dat$year<=2023),] 


#Plot raw data
ggplot(ann.dat, aes(year,mean.val)) +
  geom_line() +
  facet_wrap(~code,scale="free_y") +
  theme_bw() +
  xlab("Year") + ylab("Oyster condition (April-Oct") +
  theme(strip.background =element_rect(fill="white")) +
  xlim(c(1955,2023)) +
  theme(strip.text.x = element_text(size = 7))
#ggsave("sprsum_rawdata.pdf",width = 11, height = 7)
ggsave("fallwin_rawdata.pdf",width = 11, height = 7)

Y <- reshape2::dcast(ann.dat, code ~ year)
names = Y$code
Y = as.matrix(Y[,-which(names(Y) == "code")])


## Run DFA models (MARSS)

# # set up forms of R matrices
levels.R = c("diagonal and equal",
             "diagonal and unequal")
             #"equalvarcov") 
model.data = data.frame()


## changing convergence criterion to ensure convergence
cntl.list = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)
# fit models & store results
for(R in levels.R) {
  for(m in 1:2) {  # allowing up to 1 trend
    dfa.model = list(A="zero", R=R, m=m)
    kemz = MARSS(Y, model=dfa.model,
                 form="dfa", z.score=TRUE, control=cntl.list)
    model.data = rbind(model.data,
                       data.frame(R=R,
                                  m=m,
                                  logLik=kemz$logLik,
                                  K=kemz$num.params,
                                  AICc=kemz$AICc,
                                  stringsAsFactors=FALSE))
    assign(paste("kemz", m, R, sep="."), kemz)
  } # end m loop
} # end R loop

###################################################
### makemodeltable
###################################################
# calculate delta-AICc
model.data$delta.AICc = model.data$AICc - min(model.data$AICc)
# calculate Akaike weights
wt = exp(-0.5*model.data$delta.AICc)
model.data$Ak.wt = wt/sum(wt)
# sort results
model.tbl = model.data[order(model.data$AICc),-4]
# drop AICc from table
# calculate cumulative wts
model.tbl$Ak.wt.cum = cumsum(model.tbl$Ak.wt)
model.tbl = model.tbl[,-6]
print(model.tbl)

#write.csv(model.tbl,'sprsumDFA_AIC.csv')
write.csv(model.tbl,'fallwinDFA_AIC.csv')

model.list = list(A="zero", m=1, R= "diagonal and unequal")
dfa.mod = MARSS(Y, model=model.list, z.score=TRUE, form="dfa")

#saveRDS(dfa.mod,"sprsum_trend.rds")
saveRDS(dfa.mod,"fallwin_trend.rds")

#### To plot 1 trend model
# get CI and plot loadings...
modCI <- MARSSparamCIs(dfa.mod)
modCI 

loadings <- data.frame(names = names,
                       loading = modCI$par$Z,
                       upCI = modCI$par.upCI$Z,
                       lowCI = modCI$par.lowCI$Z)

#loadings$names <- reorder(loadings$names, loadings$order)

#write.csv(loadings,'sprsum_loadings.csv')
write.csv(loadings,'fallwin_loadings.csv')

#quartz()
loadings.plot=ggplot(loadings, aes(names, loading)) +
  geom_bar(stat="identity", fill="light grey") +
  geom_errorbar(aes(ymin=lowCI, ymax=upCI), width=0.2) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=57, vjust=1, hjust=1))+
  ylab("Loading")

#ggsave("sprsum_loadings.pdf",width = 8, height=6)
ggsave("fallwin_loadings.pdf",width = 8, height=6)


trend <- data.frame(year = 1955:2023,
                    trend = as.vector(dfa.mod$states),
                    ymin = as.vector(dfa.mod$states-1.96*dfa.mod$states.se),
                    ymax = as.vector(dfa.mod$states+1.96*dfa.mod$states.se))

#write.csv(trend,'sprsum_model.trend.csv')
write.csv(trend,'fallwin_model.trend.csv')

#quartz()
trend.plot <- ggplot(trend, aes(year, trend)) +
  geom_ribbon(aes(ymin=ymin, ymax=ymax), fill="grey90") +
  geom_line(color="black") +
  geom_point(color="black") +
  geom_hline(yintercept = 0) +
  theme(axis.title.x = element_blank()) +
  ylab("Trend") +
  scale_x_continuous(breaks = seq(1955, 2023, 5))
#ggsave("sprsum_trend.pdf",width = 8, height=6)
ggsave("fallwin_trend.pdf",width = 8, height=6)



# load PDO data and average May-Sept
dat<-read.csv("pdo.csv") 
dat<-dat[dat$year>=1955&dat$year<=2023,]
pdo.dat=dat[dat$month %in% c("5","6","7","8","9"),]

ann.pdo=pdo.dat %>%
  group_by(year) %>%
  summarize(mean.val= mean(pdo, na.rm=T)
  )
write.csv(ann.pdo,'sprsum.pdo.csv')

# load PDO data and average Oct to April
dat<-read.csv("pdo.csv") 
dat<-dat[dat$year>=1955&dat$year<=2023,]
pdo.dat=dat[dat$month %in% c("10","11","12","1","2","3","4"),]

pdo.dat$win.yr <- ifelse(pdo.dat$month %in% c("10", "11","12"), pdo.dat$year+1, pdo.dat$year)
ann.pdo=pdo.dat %>%
  group_by(win.yr) %>%
  summarize(mean.val= mean(pdo, na.rm=T)
  )
write.csv(ann.pdo,'win.pdo.csv')


# load upwelling data data and average May-Sept
dat<-read.csv("upw.csv") 
dat=dat[dat$lat %in% c("45N"),]   
dat=dat[-c(1:2)]
melt.dat=melt(dat, id.vars = c("year"))
names(melt.dat)[2:3]=c("month","value")

upw.dat <- melt.dat[melt.dat$month %in% c("may","jun","jul","aug","sep"),]   
upw.dat<-upw.dat[upw.dat$year>=1955&upw.dat$year<=2023,]

ann.upw=upw.dat %>%
  group_by(year) %>%
  summarize(mean.val= mean(value, na.rm=T)
  )
write.csv(ann.upw,'sprsum.upw.csv')

# load upwelling data and average Oct to April
dat<-read.csv("upw.csv") 
dat=dat[dat$lat %in% c("45N"),]   
dat=dat[-c(1:2)]
melt.dat=melt(dat, id.vars = c("year"))
names(melt.dat)[2:3]=c("month","value")

upw.dat=melt.dat[melt.dat$month %in% c("oct", "nov","dec","jan","feb","mar","apr"),]
upw.dat<-upw.dat[upw.dat$year>=1955&upw.dat$year<=2023,]
upw.dat$win.yr <- ifelse(upw.dat$month %in% c("oct", "nov","dec"), upw.dat$year+1, upw.dat$year)
ann.upw=upw.dat %>%
  group_by(win.yr) %>%
  summarize(mean.val= mean(value, na.rm=T)
  )
write.csv(ann.upw,'win.upw.csv')

