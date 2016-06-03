#######################################################################################################################################
#################################### HENDERSON SKINK DENSITY ESTIMATION ####################################################
#######################################################################################################################################


# written by steffen.oppel@rspb.org.uk
# adapted from Aquatic Warbler monitoring script in 2013
# separate models for each species
# last update 17 Dec 2015



############### LOAD THE REQUIRED PACKAGES #########################################################

library(unmarked)
library(Hmisc)
library(reshape)
library(slam)
library(readxl)
library(RODBC)



#####################################################################################################################################################
#############    LOAD RAW DATA AND MANIPULATE FOR UNMARKED FORMAT      ##############################################################################
#####################################################################################################################################################


############### SET THE WORKING DIRECTORY #########################################################

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Henderson\\Data")
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Henderson\\Data")

rawdatXL<-read_excel("REPTILE_DATA_ENTRY.xlsx", sheet = "Pointcounts")
rawdatXL$Date<-as.Date(rawdatXL$Date, format="%Y-%m-%d")
head(rawdatXL)
dim(rawdatXL)
dat<-rawdatXL[,1:19]



### INPUT DATA ####
library(RODBC)
SP<-odbcConnectAccess2007('HendersonData_2015.accdb')
siteCov <- sqlQuery(SP, "SELECT * FROM siteCov") 
panCov <- sqlQuery(SP, "SELECT * FROM pandanus_cover")  
odbcClose(SP)

siteCov$pandanus<-panCov$cover[match(siteCov$LocationID, panCov$LocationID)]
siteCov$pandanus[is.na(siteCov$pandanus)]<-0

head(dat)
head(siteCov)


### FORMAT TIME AND DATES ####

dat$dtime<-as.numeric(difftime(dat$Time,min(dat$Time))/3600)			### daytime in hrs since first observation
dat$DAY<-as.numeric(dat$Date)-min(as.numeric(dat$Date))				### day as number from first day of observation




### DEFINE NUMBER OF SPECIES AND POINTS ###
species<-c("EC","CP","DS","LN")
locs<-unique(dat$Point)



#####################################################################################################################################################
#############    CREATE UNMARKED DATA FRAME       ###################################################################################################
#####################################################################################################################################################

distdat1<-dat[,c(2,1,4,5,8,9,21,20,10,12,14,16)]
distdat1$dist<-1
names(distdat1)[9:12]<-species
distdat2<-dat[,c(2,1,4,5,8,9,21,20,11,13,15,17)]
distdat2$dist<-2
names(distdat2)[9:12]<-species
dist_bands<-rbind(distdat1,distdat2)





######################################################################################################
# 
# 3. RUN MODELS FOR brown-tailed copper-striped skink Emoia cyanura 
# 
#######################################################################################################

s=1

OBSDAT<-distdat1[,c(1:8)]
names(OBSDAT)[c(3,6)]<-c("morning","cloud")

SITEDAT<-aggregate(DAY~Point, OBSDAT, FUN=mean)			#### needs some meaningful covariates!!

SITEDAT$canopy<-siteCov$canopy_cover[match(siteCov$LocationID,SITEDAT$Point)]
SITEDAT$density<-siteCov$understory_dens[match(siteCov$LocationID,SITEDAT$Point)]
SITEDAT$substrate<-siteCov$substrate[match(siteCov$LocationID,SITEDAT$Point)]
SITEDAT$litter<-siteCov$leaf_litter[match(siteCov$LocationID,SITEDAT$Point)]
SITEDAT$treesize<-siteCov$AvgOfDBH[match(siteCov$LocationID,SITEDAT$Point)]
SITEDAT$pandanus<-siteCov$pandanus[match(siteCov$LocationID,SITEDAT$Point)]


obs.y<-cast(dist_bands, Point~Count+dist, value=species[s])


#### PREPARE NUMERIC COVARIATES AND MAKE SURE ALL IS IN SAME ORDER ######

OBSDAT<-as.data.frame(OBSDAT[order(OBSDAT$Point, OBSDAT$Count),])
OBSDAT$Wind<-as.numeric(ifelse(OBSDAT$Wind=="Moderate",1,0))
OBSDAT$morning<-as.numeric(ifelse(OBSDAT$morning=="am",1,0))
OBSDAT$Observer<-as.numeric(ifelse(OBSDAT$Observer=="SH",1,0))

SITEDAT<-SITEDAT[order(SITEDAT$Point),]
obs.y<-obs.y[order(obs.y$Point),]

BIRD.d<-as.matrix(obs.y[,2:13], dimnames=NULL)							      # convert to matrix which is required for 'unmarked'
head(BIRD.d)
dim(BIRD.d)
str(SITEDAT)
str(OBSDAT)



##### COMBINE RESPONSE AND OBSERVATION COVARIATES TO UNMARKED FRAME AND STANDARDIZE NUMERIC COVARIATES #######

REPdis<-unmarkedFrameGDS(y=BIRD.d, siteCovs=SITEDAT, yearlySiteCovs=OBSDAT, numPrimary=6, dist.breaks=c(0, 2, 6), survey='point', unitsIn='m')
siteCovs(REPdis)[,c(3,6,7,8)] <- scale(siteCovs(REPdis)[,c(3,6,7,8)])
yearlySiteCovs(REPdis)[,c(6:8)] <- scale(yearlySiteCovs(REPdis)[,c(6:8)])
summary(REPdis)



###################################################################################################
######## ANALYSIS OF DATA #########################################################################
###################################################################################################

### in gdistsamp, first formula is for abundance, second formula is for availability, third formula is for detection  ####

#####  BIOLOGICALLY PLAUSIBLE MODELS FOR DETECTION  #####
null <- gdistsamp(~ 1, ~ 1, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
am_obs <- gdistsamp(~ 1, ~ morning, ~Observer, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
day_obs <- gdistsamp(~ 1, ~ DAY, ~Observer, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
time_obs <- gdistsamp(~ 1, ~ dtime, ~Observer, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
cloud_obs <- gdistsamp(~ 1, ~ cloud, ~Observer, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)

am_windobs <- gdistsamp(~ 1, ~ morning, ~Observer+Wind, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
day_windobs <- gdistsamp(~ 1, ~ DAY, ~Observer+Wind, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
time_windobs <- gdistsamp(~ 1, ~ dtime, ~Observer+Wind, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
cloud_windobs <- gdistsamp(~ 1, ~ cloud, ~Observer+Wind, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)

am <- gdistsamp(~ 1, ~ morning, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
day <- gdistsamp(~ 1, ~ DAY, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
time <- gdistsamp(~ 1, ~ dtime, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
cloud <- gdistsamp(~ 1, ~ cloud, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)

obs <- gdistsamp(~ 1, ~ 1, ~Observer, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
wind <- gdistsamp(~ 1, ~ 1, ~Wind, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)

fl <- fitList(null,am_obs,day_obs,time_obs,cloud_obs,am_windobs,cloud_windobs,time_windobs,day_windobs, am, day, time, cloud, obs, wind) 
ms <- modSel(fl, nullmod="null")
ms



### USE BEST OBSERVATION MODEL TO CONSTRUCT PLAUSIBLE DENSITY MODELS

canop <- gdistsamp(~ canopy, ~ DAY, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
litt <- gdistsamp(~ litter, ~ DAY, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
tree <- gdistsamp(~ treesize, ~ DAY, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
subs <- gdistsamp(~ substrate, ~ DAY, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
dens <- gdistsamp(~ density, ~ DAY, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
pand <- gdistsamp(~ pandanus, ~ DAY, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)


fl <- fitList(day,canop,litt,tree,subs,dens,pand) 
ms <- modSel(fl, nullmod="day")
ms


###################################################################################################
######## SUMMARY OF TOP MODEL #########################################################################
###################################################################################################

out<-summary(day)

out_tab<-rbind(out$lambda, out$phi, out$det)
write.table(out_tab, "clipboard", sep="\t")


###################################################################################################
######## GOF TEST OF TOP MODEL #########################################################################
###################################################################################################

# Function adopted from Sillet et al. 2012
freeTuke <- function(fm) {
    observed <- getY(fm@data)
    expected <- fitted(fm)
    sum((sqrt(observed) - sqrt(expected))^2)
    }

# should set nsim=100 or more, but nsim=200 takes a few hours to run
pb <- parboot(top, freeTuke, nsim=200, report=1)



########################################################
#### DENSITY ESTIMATES FOR EACH SPECIES ####
########################################################
# type='state' MUST BE type='lambda' for gdistsamp!

new<-data.frame(DAY=0, obs=1)
RESULT<-predict(day,type='lambda',newdata=new,appendData=T)
names(RESULT)[1:4]<-c('density','se','lcl','ucl')
RESULT$Species<-species[s]
write.table(RESULT, "clipboard", sep="\t", row.names=F)






######################################################################################################
# 
# 4. RUN MODELS FOR oceanic snake-eyed skink Cryptoblepharus poecilopleurus
# 
#######################################################################################################

s=2

OBSDAT<-distdat1[,c(1:8)]
names(OBSDAT)[c(3,6)]<-c("morning","cloud")

SITEDAT<-aggregate(DAY~Point, OBSDAT, FUN=mean)			#### needs some meaningful covariates!!

SITEDAT$canopy<-siteCov$canopy_cover[match(siteCov$LocationID,SITEDAT$Point)]
SITEDAT$density<-siteCov$understory_dens[match(siteCov$LocationID,SITEDAT$Point)]
SITEDAT$substrate<-siteCov$substrate[match(siteCov$LocationID,SITEDAT$Point)]
SITEDAT$litter<-siteCov$leaf_litter[match(siteCov$LocationID,SITEDAT$Point)]
SITEDAT$pandanus<-siteCov$pandanus[match(siteCov$LocationID,SITEDAT$Point)]
SITEDAT$treesize<-siteCov$AvgOfDBH[match(siteCov$LocationID,SITEDAT$Point)]

obs.y<-cast(dist_bands, Point~Count+dist, value=species[s])


#### PREPARE NUMERIC COVARIATES AND MAKE SURE ALL IS IN SAME ORDER ######

OBSDAT<-as.data.frame(OBSDAT[order(OBSDAT$Point, OBSDAT$Count),])
OBSDAT$Wind<-as.numeric(ifelse(OBSDAT$Wind=="Moderate",1,0))
OBSDAT$morning<-as.numeric(ifelse(OBSDAT$morning=="am",1,0))
OBSDAT$Observer<-as.numeric(ifelse(OBSDAT$Observer=="SH",1,0))

SITEDAT<-SITEDAT[order(SITEDAT$Point),]
obs.y<-obs.y[order(obs.y$Point),]

BIRD.d<-as.matrix(obs.y[,2:13], dimnames=NULL)							      # convert to matrix which is required for 'unmarked'
head(BIRD.d)
dim(BIRD.d)
str(SITEDAT)
str(OBSDAT)



##### COMBINE RESPONSE AND OBSERVATION COVARIATES TO UNMARKED FRAME AND STANDARDIZE NUMERIC COVARIATES #######

REPdis<-unmarkedFrameGDS(y=BIRD.d, siteCovs=SITEDAT, yearlySiteCovs=OBSDAT, numPrimary=6, dist.breaks=c(0, 2, 6), survey='point', unitsIn='m')
siteCovs(REPdis)[,c(3,6,7,8)] <- scale(siteCovs(REPdis)[,c(3,6,7,8)])
yearlySiteCovs(REPdis)[,c(6:8)] <- scale(yearlySiteCovs(REPdis)[,c(6:8)])
summary(REPdis)



###################################################################################################
######## ANALYSIS OF DATA #########################################################################
###################################################################################################

### in gdistsamp, first formula is for abundance, second formula is for availability, third formula is for detection  ####

#####  BIOLOGICALLY PLAUSIBLE MODELS FOR DETECTION  #####
null <- gdistsamp(~ 1, ~ 1, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
am_obs <- gdistsamp(~ 1, ~ morning, ~Observer, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
day_obs <- gdistsamp(~ 1, ~ DAY, ~Observer, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
time_obs <- gdistsamp(~ 1, ~ dtime, ~Observer, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
cloud_obs <- gdistsamp(~ 1, ~ cloud, ~Observer, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)

am_windobs <- gdistsamp(~ 1, ~ morning, ~Observer+Wind, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
day_windobs <- gdistsamp(~ 1, ~ DAY, ~Observer+Wind, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
time_windobs <- gdistsamp(~ 1, ~ dtime, ~Observer+Wind, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
cloud_windobs <- gdistsamp(~ 1, ~ cloud, ~Observer+Wind, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)

am <- gdistsamp(~ 1, ~ morning, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
day <- gdistsamp(~ 1, ~ DAY, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
time <- gdistsamp(~ 1, ~ dtime, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
cloud <- gdistsamp(~ 1, ~ cloud, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)

obs <- gdistsamp(~ 1, ~ 1, ~Observer, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
wind <- gdistsamp(~ 1, ~ 1, ~Wind, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)

fl <- fitList(null,am_obs,day_obs,time_obs,cloud_obs,am_windobs,cloud_windobs,time_windobs,day_windobs, am, day, time, cloud, obs, wind) 
ms <- modSel(fl, nullmod="null")
ms



### USE BEST OBSERVATION MODEL TO CONSTRUCT PLAUSIBLE DENSITY MODELS

canop <- gdistsamp(~ canopy, ~ cloud, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
litt <- gdistsamp(~ litter, ~ cloud, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
tree <- gdistsamp(~ treesize, ~ cloud, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
subs <- gdistsamp(~ substrate, ~ cloud, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
dens <- gdistsamp(~ density, ~ cloud, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
pand <- gdistsamp(~ pandanus, ~ cloud, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)


fl <- fitList(cloud,canop,litt,tree,subs,dens,pand) 
ms <- modSel(fl, nullmod="cloud")
ms


###################################################################################################
######## SUMMARY OF TOP MODEL #########################################################################
###################################################################################################

out<-summary(subs)

out_tab<-rbind(out$lambda, out$phi, out$det)
write.table(out_tab, "clipboard", sep="\t")


###################################################################################################
######## GOF TEST OF TOP MODEL #########################################################################
###################################################################################################

# Function adopted from Sillet et al. 2012
freeTuke <- function(fm) {
    observed <- getY(fm@data)
    expected <- fitted(fm)
    sum((sqrt(observed) - sqrt(expected))^2)
    }

# should set nsim=100 or more, but nsim=200 takes a few hours to run
pb <- parboot(top, freeTuke, nsim=200, report=1)



########################################################
#### DENSITY ESTIMATES FOR EACH SPECIES ####
########################################################
# type='state' MUST BE type='lambda' for gdistsamp!

new<-data.frame(cloud=0, obs=1, substrate=unique(SITEDAT$substrate))
RESULT<-predict(subs,type='lambda',newdata=new,appendData=T)
names(RESULT)[1:4]<-c('density','se','lcl','ucl')
RESULT$Species<-species[s]
write.table(RESULT, "clipboard", sep="\t", row.names=F)







######################################################################################################
# 
# 5. RUN MODELS FOR DS
# 
#######################################################################################################

s=3

OBSDAT<-distdat1[,c(1:8)]
names(OBSDAT)[c(3,6)]<-c("morning","cloud")

SITEDAT<-aggregate(DAY~Point, OBSDAT, FUN=mean)			#### needs some meaningful covariates!!

SITEDAT$canopy<-siteCov$canopy_cover[match(siteCov$LocationID,SITEDAT$Point)]
SITEDAT$density<-siteCov$understory_dens[match(siteCov$LocationID,SITEDAT$Point)]
SITEDAT$substrate<-siteCov$substrate[match(siteCov$LocationID,SITEDAT$Point)]
SITEDAT$litter<-siteCov$leaf_litter[match(siteCov$LocationID,SITEDAT$Point)]
SITEDAT$pandanus<-siteCov$pandanus[match(siteCov$LocationID,SITEDAT$Point)]
SITEDAT$treesize<-siteCov$AvgOfDBH[match(siteCov$LocationID,SITEDAT$Point)]

obs.y<-cast(dist_bands, Point~Count+dist, value=species[s])


#### PREPARE NUMERIC COVARIATES AND MAKE SURE ALL IS IN SAME ORDER ######

OBSDAT<-as.data.frame(OBSDAT[order(OBSDAT$Point, OBSDAT$Count),])
OBSDAT$Wind<-as.numeric(ifelse(OBSDAT$Wind=="Moderate",1,0))
OBSDAT$morning<-as.numeric(ifelse(OBSDAT$morning=="am",1,0))
OBSDAT$Observer<-as.numeric(ifelse(OBSDAT$Observer=="SH",1,0))

SITEDAT<-SITEDAT[order(SITEDAT$Point),]
obs.y<-obs.y[order(obs.y$Point),]

BIRD.d<-as.matrix(obs.y[,2:13], dimnames=NULL)							      # convert to matrix which is required for 'unmarked'
head(BIRD.d)
dim(BIRD.d)
str(SITEDAT)
str(OBSDAT)



##### COMBINE RESPONSE AND OBSERVATION COVARIATES TO UNMARKED FRAME AND STANDARDIZE NUMERIC COVARIATES #######

REPdis<-unmarkedFrameGDS(y=BIRD.d, siteCovs=SITEDAT, yearlySiteCovs=OBSDAT, numPrimary=6, dist.breaks=c(0, 2, 6), survey='point', unitsIn='m')
siteCovs(REPdis)[,c(3,6,7,8)] <- scale(siteCovs(REPdis)[,c(3,6,7,8)])
yearlySiteCovs(REPdis)[,c(6:8)] <- scale(yearlySiteCovs(REPdis)[,c(6:8)])
summary(REPdis)

###################################################################################################
######## ANALYSIS OF DATA #########################################################################
###################################################################################################

### in gdistsamp, first formula is for abundance, second formula is for availability, third formula is for detection  ####

#####  BIOLOGICALLY PLAUSIBLE MODELS FOR DETECTION  #####
null <- gdistsamp(~ 1, ~ 1, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
am_obs <- gdistsamp(~ 1, ~ morning, ~Observer, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
day_obs <- gdistsamp(~ 1, ~ DAY, ~Observer, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
time_obs <- gdistsamp(~ 1, ~ dtime, ~Observer, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
cloud_obs <- gdistsamp(~ 1, ~ cloud, ~Observer, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)

am_windobs <- gdistsamp(~ 1, ~ morning, ~Observer+Wind, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
#day_windobs <- gdistsamp(~ 1, ~ DAY, ~Observer+Wind, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
time_windobs <- gdistsamp(~ 1, ~ dtime, ~Observer+Wind, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
cloud_windobs <- gdistsamp(~ 1, ~ cloud, ~Observer+Wind, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)


am <- gdistsamp(~ 1, ~ morning, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
day <- gdistsamp(~ 1, ~ DAY, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
time <- gdistsamp(~ 1, ~ dtime, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
cloud <- gdistsamp(~ 1, ~ cloud, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)

obs <- gdistsamp(~ 1, ~ 1, ~Observer, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
wind <- gdistsamp(~ 1, ~ 1, ~Wind, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)

fl <- fitList(null,am_obs,day_obs,time_obs,cloud_obs,am_windobs,cloud_windobs,time_windobs,am, day, time, cloud, obs, wind) 
ms <- modSel(fl, nullmod="null")
ms



### USE BEST OBSERVATION MODEL TO CONSTRUCT PLAUSIBLE DENSITY MODELS

canop <- gdistsamp(~ canopy, ~ dtime, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
litt <- gdistsamp(~ litter, ~ dtime, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
tree <- gdistsamp(~ treesize, ~ dtime, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
subs <- gdistsamp(~ substrate, ~ dtime, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
dens <- gdistsamp(~ density, ~ dtime, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
pand <- gdistsamp(~ pandanus, ~ dtime, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)


fl <- fitList(time,canop,litt,tree,subs,dens,pand) 
ms <- modSel(fl, nullmod="time")
ms


###################################################################################################
######## SUMMARY OF TOP MODEL #########################################################################
###################################################################################################

out<-summary(pand)

out_tab<-rbind(out$lambda, out$phi, out$det)
write.table(out_tab, "clipboard", sep="\t")

out<-summary(tree)

out_tab<-rbind(out$lambda, out$phi, out$det)
write.table(out_tab, "clipboard", sep="\t")

out<-summary(litt)

out_tab<-rbind(out$lambda, out$phi, out$det)
write.table(out_tab, "clipboard", sep="\t")

###################################################################################################
######## GOF TEST OF TOP MODEL #########################################################################
###################################################################################################

# Function adopted from Sillet et al. 2012
freeTuke <- function(fm) {
    observed <- getY(fm@data)
    expected <- fitted(fm)
    sum((sqrt(observed) - sqrt(expected))^2)
    }

# should set nsim=100 or more, but nsim=200 takes a few hours to run
pb <- parboot(top, freeTuke, nsim=200, report=1)



########################################################
#### DENSITY ESTIMATES FOR EACH SPECIES ####
########################################################
# type='state' MUST BE type='lambda' for gdistsamp!


new<-data.frame(dtime=0, obs=1, pandanus=seq(min(REPdis@siteCovs$pandanus),max(REPdis@siteCovs$pandanus),1))
RESULT<-predict(pand,type='lambda',newdata=new,appendData=T)
names(RESULT)[1:4]<-c('density','se','lcl','ucl')
RESULT$Species<-species[s]
RESULT$pandanus<-(RESULT$pandanus*sd(SITEDAT$treesize))+mean(SITEDAT$pandanus)
write.table(RESULT, "clipboard", sep="\t", row.names=F)

new<-data.frame(dtime=0, obs=1, treesize=seq(min(REPdis@siteCovs$treesize),max(REPdis@siteCovs$treesize),1))
RESULT<-predict(tree,type='lambda',newdata=new,appendData=T)
names(RESULT)[1:4]<-c('density','se','lcl','ucl')
RESULT$Species<-species[s]
RESULT$treesize<-(RESULT$treesize*sd(SITEDAT$treesize))+mean(SITEDAT$treesize)
write.table(RESULT, "clipboard", sep="\t", row.names=F)

new<-data.frame(dtime=0, obs=1,litter=seq(min(REPdis@siteCovs$litter),max(REPdis@siteCovs$litter),1))
RESULT<-predict(litt,type='lambda',newdata=new,appendData=T)
names(RESULT)[1:4]<-c('density','se','lcl','ucl')
RESULT$Species<-species[s]
RESULT$litter<-(RESULT$litter*sd(SITEDAT$litter))+mean(SITEDAT$litter)
write.table(RESULT, "clipboard", sep="\t", row.names=F)







######################################################################################################
# 
# 6. RUN MODELS FOR Lipinia noctua
# 
#######################################################################################################

s=4

OBSDAT<-distdat1[,c(1:8)]
names(OBSDAT)[c(3,6)]<-c("morning","cloud")

SITEDAT<-aggregate(DAY~Point, OBSDAT, FUN=mean)			#### needs some meaningful covariates!!

SITEDAT$canopy<-siteCov$canopy_cover[match(siteCov$LocationID,SITEDAT$Point)]
SITEDAT$density<-siteCov$understory_dens[match(siteCov$LocationID,SITEDAT$Point)]
SITEDAT$substrate<-siteCov$substrate[match(siteCov$LocationID,SITEDAT$Point)]
SITEDAT$litter<-siteCov$leaf_litter[match(siteCov$LocationID,SITEDAT$Point)]
SITEDAT$pandanus<-siteCov$pandanus[match(siteCov$LocationID,SITEDAT$Point)]
SITEDAT$treesize<-siteCov$AvgOfDBH[match(siteCov$LocationID,SITEDAT$Point)]

obs.y<-cast(dist_bands, Point~Count+dist, value=species[s])


#### PREPARE NUMERIC COVARIATES AND MAKE SURE ALL IS IN SAME ORDER ######

OBSDAT<-as.data.frame(OBSDAT[order(OBSDAT$Point, OBSDAT$Count),])
OBSDAT$Wind<-as.numeric(ifelse(OBSDAT$Wind=="Moderate",1,0))
OBSDAT$morning<-as.numeric(ifelse(OBSDAT$morning=="am",1,0))
OBSDAT$Observer<-as.numeric(ifelse(OBSDAT$Observer=="SH",1,0))

SITEDAT<-SITEDAT[order(SITEDAT$Point),]
obs.y<-obs.y[order(obs.y$Point),]

BIRD.d<-as.matrix(obs.y[,2:13], dimnames=NULL)							      # convert to matrix which is required for 'unmarked'
head(BIRD.d)
dim(BIRD.d)
str(SITEDAT)
str(OBSDAT)



##### COMBINE RESPONSE AND OBSERVATION COVARIATES TO UNMARKED FRAME AND STANDARDIZE NUMERIC COVARIATES #######

REPdis<-unmarkedFrameGDS(y=BIRD.d, siteCovs=SITEDAT, yearlySiteCovs=OBSDAT, numPrimary=6, dist.breaks=c(0, 2, 6), survey='point', unitsIn='m')
siteCovs(REPdis)[,c(3,6,7,8)] <- scale(siteCovs(REPdis)[,c(3,6,7,8)])
yearlySiteCovs(REPdis)[,c(6:8)] <- scale(yearlySiteCovs(REPdis)[,c(6:8)])
summary(REPdis)




###################################################################################################
######## ANALYSIS OF DATA #########################################################################
###################################################################################################

### in gdistsamp, first formula is for abundance, second formula is for availability, third formula is for detection  ####

#####  BIOLOGICALLY PLAUSIBLE MODELS FOR DETECTION  #####
null <- gdistsamp(~ 1, ~ 1, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
am_obs <- gdistsamp(~ 1, ~ morning, ~Observer, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
day_obs <- gdistsamp(~ 1, ~ DAY, ~Observer, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
time_obs <- gdistsamp(~ 1, ~ dtime, ~Observer, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
cloud_obs <- gdistsamp(~ 1, ~ cloud, ~Observer, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)

am <- gdistsamp(~ 1, ~ morning, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
day <- gdistsamp(~ 1, ~ DAY, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
time <- gdistsamp(~ 1, ~ dtime, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
cloud <- gdistsamp(~ 1, ~ cloud, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)

obs <- gdistsamp(~ 1, ~ 1, ~Observer, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)

fl <- fitList(null,am_obs,day_obs,time_obs,cloud_obs,am, day, time, cloud, obs) 
ms <- modSel(fl, nullmod="null")
ms



### USE BEST OBSERVATION MODEL TO CONSTRUCT PLAUSIBLE DENSITY MODELS

canop <- gdistsamp(~ canopy, ~ 1, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
litt <- gdistsamp(~ litter, ~ 1, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
tree <- gdistsamp(~ treesize, ~ 1, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
subs <- gdistsamp(~ substrate, ~ 1, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
dens <- gdistsamp(~ density, ~ 1, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)
pand <- gdistsamp(~ pandanus, ~ 1, ~1, data= REPdis, keyfun="halfnorm", output="density", unitsOut="ha", mixture = "P", K=25)


fl <- fitList(null,canop,litt,tree,subs,dens,pand) 
ms <- modSel(fl, nullmod="null")
ms


###################################################################################################
######## SUMMARY OF TOP MODEL #########################################################################
###################################################################################################

out<-summary(null)

out_tab<-rbind(out$lambda, out$phi, out$det)
write.table(out_tab, "clipboard", sep="\t")


###################################################################################################
######## GOF TEST OF TOP MODEL #########################################################################
###################################################################################################

# Function adopted from Sillet et al. 2012
freeTuke <- function(fm) {
    observed <- getY(fm@data)
    expected <- fitted(fm)
    sum((sqrt(observed) - sqrt(expected))^2)
    }

# should set nsim=100 or more, but nsim=200 takes a few hours to run
pb <- parboot(top, freeTuke, nsim=200, report=1)



########################################################
#### DENSITY ESTIMATES FOR EACH SPECIES ####
########################################################
# type='state' MUST BE type='lambda' for gdistsamp!

new<-data.frame(day=0, obs=1)
RESULT<-predict(null,type='lambda',newdata=new,appendData=T)
names(RESULT)[1:4]<-c('density','se','lcl','ucl')
RESULT$Species<-species[s]
write.table(RESULT, "clipboard", sep="\t", row.names=F)


