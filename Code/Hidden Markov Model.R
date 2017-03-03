# clean everything first
rm(list=ls())

#libraries
library(adehabitatLT)
library(geosphere)
library(moveHMM)
library(rworldmap)
library(rworldxtra) # creates high resolution maps 
# Plot Lat/Long Points On A Google Map
# load the ggplot2 package
library(ggplot2)
# load the ggmap package
library(ggmap)

setwd("C:\\Users\\akane\\Desktop\\Science\\Manuscripts\\White-backed Vulture Tracking Data\\AWBV-Tracking-Data\\Data")

data <- read.table("476.csv", header=T,sep=",")
head(data)
data$id<-"ID1"
length(data$DateTime)
#clean data (Adam's code)
data$DateTime<-as.POSIXct(data$DateTime, format= "%d-%m-%y %H:%M", tz = "UTC")# tz="Africa/Johannesburg")
date1 <- strptime("2015-04-08","%Y-%m-%d")
data<-subset(data, data$DateTime  > date1)
names(data)[names(data) == 'Longitude_E'] <- 'lon'
names(data)[names(data) == 'Latitude_N'] <- 'lat'
names(data)[names(data) == 'Altitude_m'] <- 'alt'
# remove duplicated rows
data<-data[!duplicated(data$DateTime), ] 
length(data$DateTime)

# remove points that have extreme changes in latitude
data <- transform(data, latDiff = c(lat[1], abs(diff(lat)) ))
head(data)
hist(data$latDiff)
length(data$lat)
data<-data[data$latDiff < 0.04, ]
length(data[data$latDiff < 0.04, ]$lat)
head(data)
length(data$lat)

# remove points that have extreme changes in longitude
data <- transform(data, lonDiff = c(lon[1], abs(diff(lon)) ))
head(data)
hist(data$lonDiff)
length(data$lon)
data<-data[data$lonDiff < 0.04, ]
length(data[data$lonDiff < 0.04, ]$lon)
head(data)
length(data$lon)

# plot the data
newmap <- getMap(resolution = "high")
plot(newmap,
     xlim = c(30, 32),
     ylim = c(-30, -21.5),
     asp = 1)
points(data$lon,data$lat, cex = .6)
#---------------------------------------------
# alternative plotting option with google maps 
#---------------------------------------------
# load the data points to be plotted
gps <- data

# get the map from google maps, centered on the median long/lat. 
mapImageData <- get_googlemap(
  center = c(lon = median(gps$lon), lat = median(gps$lat)),
  zoom = 6,
  maptype = c("terrain")
)

# plot the points on the map in red (the aes has some problems inheriting from previous models, so that is why we've FALSE'd the inherit)
ggmap(mapImageData, extent = "device") +
  geom_point(inherit.aes = FALSE, aes(x = gps$lon, y = gps$lat),
             data = gps,
             colour = gps,
             size = 0.01,
             pch = 16
  )

#

#calculate turning angle (relative, in radians)
tr<-as.ltraj(data.frame(X=data$lon,Y=data$lat),id=data$id,typeII=F) #create trajectory
head(tr[[1]])
data$theta<-tr[[1]]$rel.angle    #extract turning angle
hist(data$theta)
summary(data$theta)

#calculate step length in meters
data$prevlat<-c(NA,data$lat[-length(data$lat)])
data$prevlon<-c(NA,data$lon[-length(data$lon)])
data$step<-distMeeus(cbind(data$lon,data$lat),cbind(data$prevlon,data$prevlat), a=6378137, f=1/298.257223563)
hist(data$step)
summary(data$step)
hist(data$step[data$step<2000])
hist(data$step[data$step>2000 &data$step<10000],breaks=1000)
hist(data$step[data$step>2000],breaks=1000)

#calculate time difference in minutes
data$timediff<-NA
data$timediff[2:nrow(data)]<-diff(data$DateTime)/60
hist(data$timediff)
data$datetime2<-strptime(data$DateTime, format= "%Y-%m-%d %H:%M:%S")
data$daynight<-0
data$daynight[data$datetime$hour>=6 & data$datetime$hour<18]<-1
hist(data$timediff[data$daynight==1])

#prepare data with moveHMM
trackData <- data[,c(11,3,2)]
colnames(trackData)[1] <- c("ID")
data2 <- prepData(trackData,type="LL",coordNames=c("lon","lat"))
plot(data2,compact=T)

#apply two state HMM
## initial parameters for gamma and von Mises distributions
mu0 <- c(0.1,1) # step mean (two parameters: one for each state)
sigma0 <- c(0.1,1) # step SD
zeromass0 <- c(0.1,0.05) # step zero-mass
stepPar0 <- c(mu0,sigma0,zeromass0)
angleMean0 <- c(pi,0) # angle mean
kappa0 <- c(1,1) # angle concentration
anglePar0 <- c(angleMean0,kappa0)

m1 <- fitHMM(data=data2,nbStates=2,stepPar0=stepPar0,anglePar0=anglePar0,
             formula=~1) # no covariate
m1
plot(m1)

states <- viterbi(m1)
states[1:25]

sp <- stateProbs(m1)
head(sp)
plotStates(m1)
# ----------------------------------------------
# Try with interpolated data
# ----------------------------------------------
# interpolate data 
# n = n normally, but can change to any value
n <- nrow(data)
data <- data.frame(
  lat = approx(x = data$DateTime, y = data$lat, n = n)$y,
  lon = approx(x = data$DateTime, y = data$lon, n = n)$y,
  time = as.POSIXct(approx(data$DateTime, 1:nrow(data), n = n)$x, origin = "1970-01-01"),
  type = "corrected"
)

# plot the interpolated data
ggplot(data, aes(x = time, y = lat, color = type)) +
  geom_line() +
  geom_point() +
  scale_x_datetime(breaks = data$time, minor_breaks = NULL, labels = format(data$time, format = "%H:%M") ) +
  theme_minimal()

#prepare data with moveHMM
trackData <- data[,c(1,2,3)]
#colnames(trackData)[4] <- c("ID")
data2 <- prepData(trackData,type="LL",coordNames=c("lon","lat"))
plot(data2,compact=T)

# plot the data
newmap <- getMap(resolution = "low")
plot(newmap,
     xlim = c(28, 32),
     ylim = c(-30, -24.5),
     asp = 1)
points(data2$x,data2$y, cex = .6)

#apply two state HMM
## initial parameters for gamma and von Mises distributions
mu0 <- c(0.1,1) # step mean (two parameters: one for each state)
sigma0 <- c(0.1,1) # step SD
zeromass0 <- c(0.1,0.05) # step zero-mass
stepPar0 <- c(mu0,sigma0,zeromass0)
angleMean0 <- c(pi,0) # angle mean
kappa0 <- c(1,1) # angle concentration
anglePar0 <- c(angleMean0,kappa0)

m1 <- fitHMM(data=data2,nbStates=2,stepPar0=stepPar0,anglePar0=anglePar0,
             formula=~1) # no covariate

m1
plot(m1)

states <- viterbi(m1)
states[1:25]

sp <- stateProbs(m1)
head(sp)
plotStates(m1)
# ---------------------------------------------------------------------------
# Partition by daytime only - 2 options 
# ---------------------------------------------------------------------------
# option 1 using lubridate
library(lubridate)
time <- dmy_hm(tracks$V1)
tracks[!(hour(time) > 18 | hour(time)< 6),]
# option 2 using base R
hm <- strftime(as.POSIXct(tracks$V1, format="%m/%d/%Y %H:%M"), "%H:%M")
tracks <- tracks["06:00" < hm & hm < "18:00",]

# ---------------------------------------------------------------------------
# linear interpolation of the data, using adehabitatLT
# ---------------------------------------------------------------------------
tr2<-as.ltraj(data.frame(X=data$lon,Y=data$lat),date=data$DateTime,id=data$id,typeII=T) #create trajectory
tstep<-60 #time step we want for the interpolation, in seconds
newtr<-redisltraj(tr2, u=tstep, type = "time")
head(newtr)
head(newtr[[1]])

# convert object of class ltraj to a dataframe 
df<-ld(newtr)
names(df)[names(df) == 'x'] <- 'lon'
names(df)[names(df) == 'y'] <- 'lat'

#calculate turning angle (relative, in radians)
#extract turning angle
hist(df$rel.angle)
summary(df$rel.angle)

#calculate step length in meters
df$prevlat<-c(NA,df$lat[-length(df$lat)])
df$prevlon<-c(NA,df$lon[-length(df$lon)])
df$step<-distMeeus(cbind(df$lon,df$lat),cbind(df$prevlon,df$prevlat), a=6378137, f=1/298.257223563)
hist(df$step)
summary(df$step)
hist(df$step[df$step<2000])
hist(df$step[df$step>2000 &df$step<10000],breaks=1000)
hist(df$step[df$step>2000],breaks=1000)

#prepare data with moveHMM
trackData <- df[,c(1,2)]
#colnames(trackData)[4] <- c("ID")
data2 <- prepData(trackData,type="LL",coordNames=c("lon","lat"))
plot(data2,compact=T)

#apply two state HMM
## initial parameters for gamma and von Mises distributions
mu0 <- c(0.1,1) # step mean (two parameters: one for each state)
sigma0 <- c(0.1,1) # step SD
zeromass0 <- c(0.1,0.05) # step zero-mass
stepPar0 <- c(mu0,sigma0,zeromass0)
angleMean0 <- c(pi,0) # angle mean
kappa0 <- c(1,1) # angle concentration
anglePar0 <- c(angleMean0,kappa0)

m1 <- fitHMM(data=data2,nbStates=2,stepPar0=stepPar0,anglePar0=anglePar0,
             formula=~1) # no covariate

m1
plot(m1)
