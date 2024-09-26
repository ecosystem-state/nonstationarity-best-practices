# PURPOSE: to process ERSST SST data for use in DLM/GAMs

# Load libs

library(tidyverse)
library(ncdf4)
library(zoo)
library(maps)
library(mapdata)
library(chron)
library(fields)
library(oce)

# Load full North Pacific sst, 1854-2022
nc <- nc_open("./GOA_salmon_case_study/Data/nceiErsstv5_ee08_74ee_6f8f.nc")

# process

ncvar_get(nc, "time")   # seconds since 1-1-1970
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))
m <- months(d)
yr <- years(d)

x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")

SST <- ncvar_get(nc, "sst", verbose = F)

SST <- aperm(SST, 3:1)

SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))
lon <- rep(x, each = length(y))
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# identify polygons for each area
# extract study area
# GOA polygon
goa.x <- c(201, 201, 205, 208, 225, 231, 201)
goa.y <- c(55, 56.5, 59, 61, 61, 55, 55)


# define cells within each polygon and plot to check
goa.sst <- as.data.frame(SST)

xp <- cbind(goa.x, goa.y)
loc=cbind(lon, lat)
check <- in.poly(loc, xp=xp)

goa.sst[,!check] <- NA

# plot to check
temp.mean <- colMeans(goa.sst, na.rm=T)
z <- t(matrix(temp.mean,length(y)))
image.plot(x,y,z, col=oceColorsPalette(64), xlab = "", ylab = "", xlim = c(170, 240), ylim = c(40, 66))
contour(x, y, z, add=T)
map('world2Hires',c('Canada', 'usa', 'USSR'), fill=T,add=T, lwd=1, col="lightyellow3")

# save names of each df for processing
sst.data.names <- c("goa_sst")

# and a vector of clean names
sst.clean.names <- c("Gulf_of_Alaska")

# loop through each df, process, summarize, combine, save

# create vector of latitudes to weight mean sst by cell area

cell.weight <- sqrt(cos(lat*pi/180))

# plot to check
hist(cell.weight, breaks = 50)
unique(cell.weight) # looks right

# create a weighted mean function
weighted.cell.mean <- function(x) weighted.mean(x, w = cell.weight, na.rm = T)

# create a function to compute monthly anomalies
monthly.anomalies <- function(x) tapply(x, m, mean) 

# # short year and month vectors for 1950-2021
# yr.1950.2021 <- yr[yr %in% 1950:2021]
# m.1950.2021 <- m[yr %in% 1950:2021]

# and define winter year
winter.year <- if_else(m  %in% c("Nov", "Dec"), as.numeric(as.character(yr)) +1, as.numeric(as.character(yr)))
winter.year <- winter.year[m  %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]

# create blank data frame for catching results
temp.time.series <- temp.anomaly.time.series <- data.frame()


# Data frame to catch all outputs
sst.data.frames <- list()
sst.data.frames[[1]] <- goa.sst

# processing loop

for(i in 1: length(sst.data.names)){
  
  # pull out sst dataframe of interest
  temp.dat <- sst.data.frames[[i]]
  
  # first, calculate monthly anomalies
  mu <- apply(temp.dat, 2, monthly.anomalies)	# compute monthly means for each time series (cell)
  mu <- mu[rep(1:12, length(m)/12),]  # replicate means matrix for each year at each location
  
  temp.anom <- temp.dat - mu   # compute matrix of anomalies
  
  # calculate weighted monthly means
  temp.monthly.anom <- apply(temp.anom, 1, weighted.cell.mean)
  
  ## now create time series of annual means
  ## unsmoothed, and two- and three-year running means
  
  # calculate monthly mean temp weighted by area  
  temp.monthly.sst <- apply(temp.dat, 1, weighted.cell.mean) 
  
  # calculate annual means
  temp.annual <- tapply(temp.monthly.sst, yr, mean) 
  
  temp.2yr <- rollmean(temp.annual, 2, fill = NA, align = "left") # for salmon - year of and year after ocean entry
  
  temp.3yr <- rollmean(temp.annual, 3, fill = NA, align = "center") # for salmon - year before, year of, and year after ocean entry
  
  # calculate winter means
  temp.winter.monthly.sst <- temp.monthly.sst[m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]
  
  temp.winter <- tapply(temp.winter.monthly.sst, winter.year, mean)
  
  # change leading and trailing years to NA b/c these are incomplete
  leading <- min(names(temp.winter))
  trailing <- max(names(temp.winter))
  
  temp.winter[names(temp.winter) %in% c(leading, trailing)] <- NA
  
  temp.winter.2yr <- rollmean(temp.winter, 2, fill = NA, align = "left")
  
  temp.winter.3yr <- rollmean(temp.winter, 3, fill = NA, align = "center")
  
  
  # combine into data frame of time series by region
  temp.time.series <- rbind(temp.time.series,
                            data.frame(region = sst.clean.names[i],
                                       year = 1854:2024,
                                       annual.unsmoothed = temp.annual,
                                       annual.two.yr.running.mean = temp.2yr,
                                       annual.three.yr.running.mean = temp.3yr,
                                       winter.unsmoothed = temp.winter[names(temp.winter) %in% 1854:2024],
                                       winter.two.yr.running.mean = temp.winter.2yr[names(temp.winter) %in% 1854:2024],
                                       winter.three.yr.running.mean = temp.winter.3yr[names(temp.winter) %in% 1854:2024]))
  
}


# Compare with previous data from Litzow et al. 2018
old.sst <- read.csv("./GOA_salmon_case_study/Data/GOA_salmon_catch.csv") %>%
  dplyr::select(sst_3yr_running_mean, year)

temp.time.series %>%
  mutate(year = rownames(.)) %>%
  dplyr::select(winter.three.yr.running.mean, year) %>%
  filter(year > 1964) %>%
  rename(sst_3yr_running_mean = winter.three.yr.running.mean) -> new.sst

rbind(old.sst %>% mutate(type = "old"), new.sst %>% mutate(type = "new")) -> plot.dat

ggplot(plot.dat, aes(as.numeric(as.character(year)), sst_3yr_running_mean, color = type))+
  geom_line(group = 1, linewidth = 1)+
  theme_bw()


# Save data
write.csv(new.sst, "./GOA_salmon_case_study/Data/winterSST_3yr_running_mean.csv")

read.csv("./GOA_salmon_case_study/Data/winterSST_3yr_running_mean.csv")
