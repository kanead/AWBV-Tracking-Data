##################################################################################
## African White-backed Vulture Tracking Data
##################################################################################
rm(list=ls()) 
library("plotKML")
library("chron")
library("sp")
library(OpenStreetMap)
library(rgdal)
library(move)
library(rworldmap)
## https://pakillo.github.io/R-GIS-tutorial/
##################################################################################
setwd("C:/Users/akane/Desktop/Science/Manuscripts/White-backed Vulture Tracking Data/AWBV-Tracking-Data/Data")
data<-read.csv("475June.csv", header=T, sep=",")
head(data)
class(data$DateTime)
##################################################################################
## make the date column into a date-time class object
##################################################################################
data$DateTime<-as.POSIXct(data$DateTime, format= "%d-%m-%y %H:%M", tz="UTC")
data$DateTime
length(data$DateTime)
##################################################################################
## rename columns
##################################################################################
names(data)[names(data) == 'Longitude_E'] <- 'lon'
names(data)[names(data) == 'Latitude_N'] <- 'lat'
names(data)[names(data) == 'Altitude_m'] <- 'alt'
names(data)
##################################################################################
## can identify outliers by plotting the data
##################################################################################
newmap <- getMap(resolution = "low")
plot(newmap,
     xlim = c(28, 32),
     ylim = c(-30, -25),
     asp = 1)
points(data$lon,data$lat, cex = .6)
##################################################################################
## clean the data 
##################################################################################
date1 <- strptime("2015-04-08","%Y-%m-%d")
data<-subset(data, data$DateTime  > date1 
             & data$alt  !="NegAlt" & data$alt != "No Fix"
                & data$lon > 29.5 & data$lon < 32.5   
                  &  data$lat > -29 & data$lat < -25)
## these coords should be checked before removal! 
##################################################################################
## check cleaned data and plot again 
##################################################################################
length(data$DateTime)
min(data$lon)
max(data$lon)
min(data$lat)
max(data$lat)

newmap <- getMap(resolution = "low")
plot(newmap,
     xlim = c(28, 32),
     ylim = c(-30, -25),
     asp = 1)
points(data$lon,data$lat, cex = .6)
##################################################################################
## add the coordinates so we can remove remaining niggley outliers 
##################################################################################
text(data$lon, data$lat, round(data$lon, 6), cex=0.8)
# -26.39311; 32.10486
# -26.01529; 31.93312
data<-subset(data, data$lat != -26.39311 & data$lat != -26.01529)
##################################################################################
## Creating a Spatial Points Dataframe 
##################################################################################
coordinates(data) <- ~lon+lat
proj4string(data) <- CRS("+proj=longlat +datum=WGS84")
data
plot(data)
##################################################################################
## Alternatively create a Move object from move package 
##################################################################################
setwd("C:/Users/akane/Desktop/Science/Manuscripts/White-backed Vulture Tracking Data/AWBV-Tracking-Data/Data")
data <- read.table("475June.csv", header=T,sep=",")
length(data$lon)
length(data$lon)
## subset the data so that we remove the duplicate DateTime values
data<-subset(data, !duplicated(DateTime))
length(data$lon)

data <- move(x=data$lon, y=data$lat,
             time=as.POSIXct(data$DateTime,format="%d-%m-%y %H:%M", tz="UTC"),
             data=data, proj=CRS("+proj=longlat +ellps=WGS84"), animal="unknown")

head(data)
class(data)
data

## distance is given in metres
summary(data)
##################################################################################
## Eagle Model - Enrico 
##################################################################################
model
{
#States:
# 1 ascending
# 2 descending
# 3 convoluted, change in Altitude, small steps
# 4 convoluted on the ground

#priors by state:
a.mu[1] ~ djl.dnorm.trunc(100,0.0005,0,10000)
a.mu[2] <- -a.mu[1]
a.mu[3] ~ dnorm(100,0.0005)
a.mu[4] <- 0 #this is not used (when on the ground, altitude is elevation)

a.sd[1] ~ dunif(0,200)
a.sd[2] <- a.sd[1]
a.sd[3] ~ dunif(a.sd[4],a.sd[1])
a.sd[4] <- 50/1.96 #manufacturer estimate of altitude error + error in elevation + error in interpolation (CI: +- 15+2.44+30m) (bird is on ground, so no process/model error)
for (i in 1:nstates){
  a.tau[i]<-1/a.sd[i]/a.sd[i]
}

rho[1] ~ dbeta(2, 1.5)	
rho[2] <- rho[1]
dev2 ~ dbeta(1, 1)	
rho[3] <- rho[1]*dev2
dev ~ dbeta(1, 1)			
rho[4] <- rho[1]*dev	

mu <- 0

logalpha[1] ~ dunif(0,log.maxalpha) 
logalpha[2] <- logalpha[1]
logalpha[3] ~ dunif(0,logalpha[1]) 
logalpha[4] ~ dunif(0,logalpha[1]) 
logbeta[1] ~ dunif(0,log.maxbeta)
logbeta[2] <- logbeta[1]
logbeta[3] ~ dunif(0,log.maxbeta)
logbeta[4] ~ dunif(0,log.maxbeta)

for (i in 1:nstates){
  alpha[i] <- exp(logalpha[i]) #scale parameter
  beta[i] <- exp(logbeta[i]) #shape parameter
  #JAGS uses different Weibull parameterization than R
  lambda[i] <- 1/pow(alpha[i],beta[i]) 
}

phi1[1:nstates] ~ ddirch(phiprior[1:nstates])   #initial state probabilities
b[1] ~ dcat(phi1[1:nstates]) #initial state

for (i in 1:nstates){
  phi[i,1:nstates] ~ ddirch(phiprior[1:nstates]) 
}

#for (k in 1:N){ #individual loop (when multiple inds)
for (t in 2:steps){ #steps[k]
  
  b[t] ~ dcat(phi[b[t-1],1:nstates])
  
  g[t] <- step(b[t]-4)   #equal to 1 if bird in state 4 (on the ground), 0 otherwise
  
  alt.mu[t] <- (alt[t-1] + a.mu[b[t]])*(1-g[t]) + elev[t]*g[t] #when bird flying (g=0), avg altitude is previous altitude + drift. When on ground (g=1), avg altitude is the elevation
  alt[t] ~ dnorm(alt.mu[t], a.tau[b[t]]) # altitude (includes process/model error and instrument error)
  
  step[t]~dweib(beta[b[t]],lambda[b[t]])
  
  #"ones" trick to sample from the Wrapped Cauchy distribution
  ones[t] <- 1
  ones[t] ~ dbern(wC[t])
  wC[t] <- ( 1/(2*Pi)*(1-rho[b[t]]*rho[b[t]])/(1+rho[b[t]]*rho[b[t]]-2*rho[b[t]]*cos(turn[t]-mu)) )/500
  
  state.cnt[1,t]<-equals(b[t],1)
  state.cnt[2,t]<-equals(b[t],2)
  state.cnt[3,t]<-equals(b[t],3)
  state.cnt[4,t]<-equals(b[t],4)
  
} #close temporal loop

#monitor convergence of states (looking at proportions)
sumb[1]<-sum(state.cnt[1,2:steps])/steps
sumb[2]<-sum(state.cnt[2,2:steps])/steps
sumb[3]<-sum(state.cnt[3,2:steps])/steps
sumb[4]<-sum(state.cnt[4,2:steps])/steps


# } #close individual loop
} #end
