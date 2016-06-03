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
##################################################################################
setwd("C:/Users/kanead/Desktop/PhD/White back Tracking Data/September/00000475")
data<-read.csv("475.csv", header=T, sep=",")
head(data)
class(data$DateTime)
##################################################################################
## make the date column into a date-time class object
##################################################################################
data$DateTime<-as.POSIXct(data$DateTime, format= "%d-%m-%y %H:%M")
data$DateTime
length(data$DateTime)
##################################################################################
## rename columns
##################################################################################
names(data)[names(data) == 'Longitude_E'] <- 'lon'
names(data)[names(data) == 'Latitude_N'] <- 'lat'
names(data)
##################################################################################
## can identify outliers by plotting the data
##################################################################################
plot(data$lon,data$lat)
##################################################################################
## clean the data 
##################################################################################
date1 <- strptime("2015-04-08","%Y-%m-%d")
data<-subset(data, data$DateTime  > date1)
data<-subset(data, data$lon < 32.5 & data$lon > 31.4)
data<-subset(data, data$lat > -28 & data$lat < -25 )
##################################################################################
## check cleaned data
##################################################################################
length(data$DateTime)
min(data$lon)
max(data$lon)
min(data$lat)
max(data$lat)

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
## Creating a Move object from move package
##################################################################################
setwd("C:/Users/kanead/Desktop/PhD/White back Tracking Data/September/00000475")
data <- read.table("475cleaned.csv", header=T,sep=",")
length(data$lon)
length(data$lon)
## subset the data so that we remove the duplicate DateTime values
data<-subset(data, !duplicated(DateTime))
length(data$lon)

data <- move(x=data$lon, y=data$lat,
             time=as.POSIXct(data$DateTime,
                             format="%d-%m-%y %H:%M", tz="UTC"),
             data=data, proj=CRS("+proj=longlat +ellps=WGS84"), animal="unknown")

head(data)
class(data)
data

## distance is given in metres
summary(data)