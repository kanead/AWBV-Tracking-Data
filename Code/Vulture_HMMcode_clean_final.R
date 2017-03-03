##################################################################################
##                African White-backed Vulture Tracking Data                    ##
##################################################################################

#libraries
library(adehabitatLT)
library(geosphere)
library(moveHMM)

#set working directory
setwd("C:\\Users\\akane\\Desktop\\Science\\Manuscripts\\White-backed Vulture Tracking Data\\AWBV-Tracking-Data\\Data")

#load data
data <- read.table("475_3Dec2016.csv", header=T,sep=",")
data$id<-"ID1"
data2 <- read.table("476_3Dec2016.csv", header=T,sep=",")
data2$id<-"ID2"
data<-rbind(data,data2)
table(data$id)

#clean data (Adam's code)
data$DateTime<-as.POSIXct(data$DateTime, format= "%d-%m-%Y %H:%M:%S", tz="UTC")
date1 <- strptime("2015-04-08","%Y-%m-%d")
names(data)[names(data) == 'Longitude_E'] <- 'lon'
names(data)[names(data) == 'Latitude_N'] <- 'lat'
names(data)[names(data) == 'Altitude_m'] <- 'alt'
tapply(data$lat,data$id,summary)
tapply(data$lon,data$id,summary)
data<-subset(data, data$DateTime  > date1
             & data$alt  !="NegAlt" & data$alt != "No Fix" & data$alt  !="2D Fix"
                & data$lon > 29 & data$lon < 34
                  &  data$lat > -29 & data$lat < -21)

#calculate step length in meters
data$prevlat<-NA
data$prevlon<-NA
data$step<-NA
for (i in unique(data$id)){
 data$prevlat[data$id==i]<-c(NA,data$lat[data$id==i][-length(data$lat[data$id==i])])
 data$prevlon[data$id==i]<-c(NA,data$lon[data$id==i][-length(data$lon[data$id==i])])
 data$step[data$id==i]<-distMeeus(cbind(data$lon[data$id==i],data$lat[data$id==i]),cbind(data$prevlon[data$id==i],data$prevlat[data$id==i]), a=6378137, f=1/298.257223563)
 }

#calculate time difference in minutes
data$timediff<-NA
for (i in unique(data$id)){
 data$timediff[data$id==i][2:nrow(data[data$id==i,])]<-diff(data$DateTime[data$id==i])/60
 }

#filter out unrealistic locations
data<-data[data$timediff!=0 | is.na(data$timediff),]
data$speedkmh<-data$step/(data$timediff*60)*60*60/1000
plot(10,length(data$speedkmh[data$speedkmh>10]),xlim=c(0,500),ylim=c(1,18000))  #at which speed in km/h the number of locations reaches an asymptote (minimum)?
for (i in seq(20,500,by=10)){
 points(i,length(data$speedkmh[data$speedkmh>i]))
 }
tmp<-data[(data$speedkmh<100 | is.na(data$speedkmh)) & !is.na(data$DateTime),] #exclude unrealistic data points

#explore time difference
tmp$gap<-NA  #in minutes
for (i in unique(tmp$id)){
 tmp$gap[tmp$id==i]<-as.numeric(c(NA,diff(tmp$DateTime[tmp$id==i])))
 }
head(tmp)#check that it is in minutes
tapply(tmp$gap,tmp$id,summary)
summary(tmp$gap)

#choose time step to use
for (i in c(0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.99)){
 print(c(i,quantile(tmp$gap,probs=i,na.rm=T)))
 }
#use 10 minute time step, which is ~90 quantile
sum(tmp$gap<10,na.rm=T)/length(tmp$gap) #88%

tmp$daynight<-0
tmp$DateTime<-as.POSIXlt(tmp$DateTime)
tmp$daynight[tmp$DateTime$hour>=6 & tmp$DateTime$hour<18]<-1
tmp<-tmp[tmp$daynight==1,]

#recalculate gap and separate where greater than 100 minutes
tmp$gap<-NA  #in minutes
for (i in unique(tmp$id)){
 tmp$gap[tmp$id==i]<-as.numeric(c(NA,diff(tmp$DateTime[tmp$id==i])))
 }
head(tmp)#check that it is in minutes
tapply(tmp$gap,tmp$id,summary)
k<-1
tmp$id_new<-NA
tmp$id_new[1]<-1
for (i in 2:nrow(tmp)){
 if(is.na(tmp$gap[i]) | tmp$gap[i]>=100){
  k<-k+1
  tmp$id_new[i]<-k
  }else{
  tmp$id_new[i]<-k
  }
 }
table(tmp$id_new)
length(which(table(tmp$id_new)>100)) #how many portions of tracks with more than 100 locations
tmp2<-tmp[tmp$id_new %in% which(table(tmp$id_new)>100),]

#regularise track (every 10 minutes)
tmp2$DateTime<-as.POSIXct(tmp2$DateTime)
tr2<-as.ltraj(data.frame(X=tmp2$lon,Y=tmp2$lat),date=tmp2$DateTime,id=tmp2$id_new,typeII=T) #create trajectory
tstep<-10*60
newtr<-redisltraj(tr2, u=tstep, type = "time")
ndat<-c()
for (i in 1:length(unique(tmp2$id_new))){
 idi<-cbind(id_new=unique(tmp2[tmp2$id_new==sort(unique(tmp2$id_new))[i],19]),newtr[[i]])
 ndat<-rbind(ndat,idi)
 }
table(ndat$id_new)

#prepare data for HMM analysis
trackData <- ndat[,c(2:3,1)]
names(trackData)[3]<-"ID"
datareg <- prepData(trackData,type="LL",coordNames=c("x","y"))
head(datareg)
table(datareg$ID)
#plot(datareg)

#initial parameters for distributions (model with 3 states)
mu0 <- c(0.1,1,1)                     #step mean
sigma0 <- c(0.1,1,1)                  #step SD
zeromass0 <- c(0.2,0.1,0.1)           #step zero-mass
stepPar0 <- c(mu0,sigma0,zeromass0)
angleMean0 <- c(pi,0,pi)              #angle mean
kappa0 <- c(1,2,2)                    #angle concentration
anglePar0 <- c(angleMean0,kappa0)

#HMM fitting 3 states
m23 <- fitHMM(data=datareg,nbStates=3,stepPar0=stepPar0,anglePar0=anglePar0, stepDist="gamma",angleDist="vm") #stepDist="weibull",angleDist="wrpcauchy",  #stepDist="gamma",angleDist="vm",
m23
AIC(m23)

#plot results 3 states
plot(m23)
states23 <- viterbi(m23)
seqs<-round(seq(1,nrow(ndat),length.out=100)) #plot subsets of data points
for (i in 2:length(seqs)){
 plot(ndat$x[seqs[i-1]:seqs[i]], ndat$y[seqs[i-1]:seqs[i]],type="l")
 points(ndat$x[seqs[i-1]:seqs[i]], ndat$y[seqs[i-1]:seqs[i]], col=states23[seqs[i-1]:seqs[i]],pch=16,cex=0.5)
 }

 plot(ndat$x, ndat$y,type="l",xlab="Longitude",ylab="Latitude")
 points(ndat$x[states23==2], ndat$y[states23==2], col=states23[states23==2],pch=16,cex=0.5)
 points(ndat$x[states23==1], ndat$y[states23==1], col=states23[states23==1],pch=16,cex=0.5)
 points(ndat$x[states23==3], ndat$y[states23==3], col=states23[states23==3],pch=16,cex=0.2)

# Add covariate to HMM
# read data back in 
 df<-read.csv("African White backed Vultures 10 min regularised-1251491397440348311.csv",header=T,sep=",")
 
 # rename columns
 names(df)[names(df) == 'MODIS.Land.Terra.Vegetation.Indices.500m.16d.NDVI'] <- 'NDVI'
 names(df)[names(df) == 'location.long'] <- 'lon'
 names(df)[names(df) == 'location.lat'] <- 'lat'
 names(df)[names(df) == 'tag.local.identifier'] <- 'ID'
 
 # prepare data with moveHMM
 datareg <- df[,c(4,5,8,11)]
 head(datareg)
 datareg <- prepData(datareg,type="LL",coordNames=c("lon","lat"))
 
 #initial parameters for distributions (model with 3 states)
 mu0 <- c(0.1,1,1)                     #step mean
 sigma0 <- c(0.1,1,1)                  #step SD
 zeromass0 <- c(0.2,0.1,0.1)           #step zero-mass
 stepPar0 <- c(mu0,sigma0,zeromass0)
 angleMean0 <- c(pi,0,pi)              #angle mean
 kappa0 <- c(1,2,2)                    #angle concentration
 anglePar0 <- c(angleMean0,kappa0)
 
 #HMM fitting 3 states
 m24 <- fitHMM(data=datareg,nbStates=3,stepPar0=stepPar0,anglePar0=anglePar0, stepDist="gamma",angleDist="vm",
               formula=~NDVI) 
m24
AIC(m24)
plot(m24)

save(m24, file = "m24_withNDVI.RData")
