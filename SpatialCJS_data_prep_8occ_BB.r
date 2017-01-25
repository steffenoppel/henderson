#######################################################################################
##  RAT SURVIVAL ESTIMATION ON HENDERSON ISLAND  - DATA PREPARATION FOR SPATIAL CJS MODEL ##
#######################################################################################
## written by steffen.oppel@rspb.org.uk on 22 July 2016
## based on book chapter 16.3 in Royle/Gardner/Sollmann book
## uses the primary periods as encounter occasions, and the number of secondary (daily) captures for each rat at each trap
## trapping effort is the number of active trap nights at each location
## UPDATE 23 July 2016: included sex and age of killed rats because number of individuals was too large for computer - REMOVED 28 NOVEMBER
## UPDATE 26 July 2016: added trap-specific recapture probability

## SHORTENED 21 Oct 2016 - temporarily abandoned on 24 OCT 2016

## RESUMED ON 25 NOV 2016 because all approaches using single nights as occasions (single session or robuist design models) produced ludicrous estimates.
## 25 NOV 2016: REMOVED AGE AND ADDED INDIVIDUAL CAPTURE HETEROGENEITY EFFECTS

## MODIFIED on 7 DEC 2016 for BeachBack

## UPDATE 25 Jan 2017: INCLUDED RAINFALL AND SNAP TRAPPING


## DATA ANALYSIS MOVED TO SEPARATE SCRIPT

library(RODBC)
library(sp)
library(rgdal)
library(maptools)
library(reshape)




################################################################################################################
######################### INPUT DATA #############################################
################################################################################################################



setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Henderson\\Data")
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Henderson\\Data")
SP<-odbcConnectAccess2007('HendersonData_2015.accdb')
rats <- sqlQuery(SP, "SELECT * FROM SECR_rat_input_Gardner2010")
ratinds <- sqlQuery(SP, "SELECT * FROM Rat_individuals") 
trap <- sqlQuery(SP, "SELECT * FROM SECR_input_traps_Gardner2010")
kills<-sqlQuery(SP, "SELECT * FROM Kill_rat_last_locs")
sessions <- sqlQuery(SP, "SELECT * FROM Sessions")
traptypes <- sqlQuery(SP, "SELECT * FROM trap_deployments")
odbcClose(SP)

head(rats)
head(trap)







################################################################################################################
### REMOVE PLATEAU FROM ALL DATA FRAMES ###
################################################################################################################

rats$area<-trap$area[match(rats$Traploc_ID, trap$Traploc_ID)]
trap<-trap[trap$area=="Beachback",]
rats<-rats[rats$area=="Beachback",]
sessions<-sessions[sessions$habitat=="Beachback",]
rats$area<-NULL
traptypes<-traptypes[traptypes$Traploc_ID %in% trap$Traploc_ID,]

T <- max(rats$TrapSession)-min(rats$TrapSession)+1   				### number of primary study periods (of 10 trap nights each)
ntraps<-dim(trap)[1]					### number of trap locations	






#############################################################################
##   READ IN SNAP TRAPPING DATA FOR SESSION 8  ####
#############################################################################

setwd("A:\\RSPB\\UKOT\\Henderson\\Data\\Rats")
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Henderson\\Data\\Rats")
files<-list.files(pattern = "\\.xlsx$")	
files[5]

tl<-read_excel(files[5], sheet="TrapData")
tl<-tl[tl$Habitat=="Beach-back",]


## NEED TO MATCH LOCATIONS IN tl TO TRAP LOCATIONS 
tl$Traploc_ID<-NA
for (l in 1: dim(tl)[1]){
matr<-spDistsN1(pts=as.matrix(trap[,4:3]),pt=as.numeric(tl[l,8:7]), longlat = T)
tl$Traploc_ID[l]<-as.character(trap$Traploc_ID[which(matr == min(matr))])
}
head(tl)


## READ IN RAT DATA AND MATCH TO TRAPLOC

krats<-read_excel(files[5], sheet="TrapCaptures")
head(krats)
names(krats)[6]<-"Rat_ID"
krats<-krats[!is.na(krats$Rat_ID),c(1,2,4:6)]
krats<-merge(krats,tl,by=c("TrapLine","TrapNumber"), all.x=T)
dim(krats)


## ADD TO general rats data set
head(rats)
krats$TrapSession<-8
krats$Sex<-ifelse(krats$Condition=="Male",'m','f')
krats$mass<-75			### randomly made up - not used in analysis
krats$Age<-'ad'			### randomly made up - not used in analysis
krats$tail<-55			### randomly made up - not used in analysis
krats$body<-105			### randomly made up - not used in analysis
head(krats)
names(krats)

rats<-rbind(rats,krats[,c(5,3,13,12,14:18)])


### KILLED RATS THAT ARE NOT AVAILABLE IN SPREADSHEED ##

kills<-kills[-(kills$Rat_ID %in% krats$Rat_ID),]
kills<-kills[kills$area=="Beachback",]
kills$TrapSession<-8
names(kills)[3]<-"Traploc_ID"
head(kills)
rats<-rbind(rats,kills[,c(1,2,10,3:8)])






################################################################################################################
##################### ARRANGE TRAP INPUT FILE #########################################################
################################################################################################################

head(trap)
trap<-trap[order(trap$Traploc_ID, decreasing=F),]
trap[is.na(trap)]<-0


### PLOT LIVE AND SNAP TRAPS
plot(trap[,4:3])
points(tl[,8:7], pch=16, cex=0.5)

### convert lat and long into UTM grid
trapSP<-SpatialPoints(trap[,4:3], proj4string=CRS("+proj=longlat +datum=WGS84"))
trapTransformed<-spTransform(trapSP, CRSobj =CRS("+proj=aeqd"))
trap[,4:3]<-trapTransformed@coords
trap$Longitude<- (trap$Longitude)*(-1)					### because long is on other side of the world, looks mirror-imaged on plot
trap$Longitude<- trap$Longitude-min(trap$Longitude)			### reduce the dimensions of the coordinates to avoid numerical problems
trap$Latitude<- trap$Latitude-min(trap$Latitude)			### reduce the dimensions of the coordinates to avoid numerical problems

### create unique number for each trap
trap$ID<-seq(1:dim(trap)[1])
head(trap)
grid<-as.matrix(trap[,4:3])
plot(grid)


### create matrix with number of trap nights for each trap in each primary period
year<-as.matrix(trap[,8:10])
year<-cbind(year,rep(2,dim(year)[1]))
names(year)<-NULL





################################################################################################################
##################### ARRANGE TRAPTYPE INPUT MATRIX #########################################################
################################################################################################################

ttyp<-year
ttyp[,]<-1			### same traps for all sessions except session 8 (=4)


### ADD SNAP TRAPS FOR SESSION 8
ttyp[,4]<-2


################################################################################################################
########################## ARRANGE CAPTURE DATA ###################################
################################################################################################################

# remove rat captures without ID
rats<-rats[!is.na(rats$Rat_ID),]
dim(rats)

# remove all rats that were only captured once (transients)
rats$capt<-1
ncapt<-aggregate(capt~Rat_ID, rats, FUN=sum)
hist(ncapt$capt,30)
exclude<-ncapt$Rat_ID[ncapt$capt==1]
dim(rats)
rats<-rats[!(rats$Rat_ID %in% exclude),]
dim(rats)

# update sex info
rats$Sex[is.na(rats$Sex)]<-'f'
rats$sexnum<-ifelse(rats$Sex=="m",1,2)




#### RE-NUMBER SESSIONS
head(rats)
rats$sess<-rats$TrapSession-3
rats$sess[rats$TrapSession==8]<-4
T <- 4   				### number of primary study periods (of 10 trap nights each)




################################################################################################################
### CALCULATE INTERVAL BETWEEN PRIMARY SESSIONS ####
################################################################################################################


sessmidtime<-aggregate(Date~TrapSession, sessions, FUN=mean)
interval<-vector()
str(sessmidtime)
for (t in 1:(T-1)){
interval[t]<-as.numeric(sessmidtime$Date[t+1]-sessmidtime$Date[t])
}




#############################################################################
##   READ IN WEATHER DATA FROM PITCAIRN IN 2015  ####
#############################################################################
## these data are provided by Pitcairn in an Excel file with two sheets per month
library(readxl)
setwd("A:\\RSPB\\UKOT\\Henderson\\Data\\Weather\\Pitcairn")
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Henderson\\Data\\Weather\\Pitcairn")
files<-list.files(pattern = "\\.xlsx$")	
files[10]

rain<-data.frame()
for (sh in seq(5,10,1)){
shit<-read_excel(files[10], sh,col_names=F,skip=7)
shit<-as.data.frame(shit)
shit<-shit[,1:2]
del<-as.numeric(min(row.names(shit[is.na(shit[,1]),])))
delete<-seq(del,dim(shit)[1],1)
shit<-shit[-(delete),]
shit[,2]<-as.numeric(shit[,2])
shit$Month=sh
rain<-rbind(rain,shit)
}
names(rain)<-c("Day","mm","Month")
rain$Date<-as.POSIXct(paste("2015",rain$Month, rain$Day, sep="-"), format="%Y-%m-%d")


### CALCULATE SUM OF RAINFALL FOR EACH OCCASION AND HABITAT - AS SURVIVAL COVARIATE
survrain<-matrix(0,ncol=2,nrow=3)
for (sess in 4:5){
	x<-sessions[sessions$TrapSession==sess,]
	x1<-sessions[sessions$TrapSession==sess+1,]
	xa<-x[x$habitat==habs[1],]
	xa1<-x1[x1$habitat==habs[1],]
	sessrange<-seq.Date(as.Date(min(xa$Date, na.rm=T)),as.Date(min(xa1$Date, na.rm=T)),1)
	survrain[sess-3,1]<-sum(rain$mm[as.Date(rain$Date) %in% sessrange], na.rm=T)/10			## rainfall is provided in micrometers
	}

sessrange<-seq.Date(as.Date(min(xa1$Date, na.rm=T)),as.Date(max(rats$Date, na.rm=T)),1)
interval[3]<-as.numeric(difftime(max(sessrange),min(sessrange),'days'))
survrain[3,1]<-sum(rain$mm[as.Date(rain$Date) %in% sessrange], na.rm=T)/10	
survrain

### Standardise for JAGS ###
raininput<-(survrain[,1]-mean(survrain[,1]))/(sd(survrain[,1]))






################################################################################################################
### CREATE CAPTURE ARRAY FOR TRAPS x INDIVIDUALS x OCCASIONS  ##
################################################################################################################

# sort data and use numeric RatID
rats<-rats[order(rats$Rat_ID, decreasing=F),]
ratids<-unique(rats$Rat_ID)
rats$ID<-match(rats$Rat_ID,ratids)

# match the trap IDs with numeric trap ID
rats$TrapID<-trap$ID[match(rats$Traploc_ID, trap$Traploc_ID)]

# create the data input matrix, which is an array that lists the number of captures of each rat at each trap in each session
rats$capt<-1

y<-array(0,dim=c(max(rats$ID),ntraps,T))						# array for the rat detections
K<-array(0,dim=c(ntraps,T))								# array for the trap deployments

for (p in 1:T){			## start loop over all primary sessions
K[,p]<-trap[,p+7]			## write the activity for all traps into K vector
yp<-cast(rats[rats$sess==p,], formula= ID~TrapID, value='capt', fun.aggregate=sum)		## extract and cast rat captures
trapsused<-seq(1:dim(yp)[2])	## sequential number of the traps in which rats were caught (all others are not in this matrix)

for (m in unique(yp$ID)){			## start loop over all rats captured in this session
ypp<-yp[yp$ID==m,]				## extract all capture locations
ypp2<-as.numeric(ypp)

for (t in names(ypp)[2:length(ypp)]){	## start loop over all capture locations
tn<-as.numeric(t)					## this is the trap number and the index into which the observation has to go into the y array
vecpos<-trapsused[match(t,names(ypp))]	## this is the position in the capture number vector for that individual, matched to the TrapID by vector position
y[m,tn,p]<-ypp2[vecpos]
}							## end loop over capture locations
}							## end loop over individuals
}							## end loop over primary sessions




### DEAL WITH SNAP TRAP SESSION ###
K[,4]<-0
head(tl)
head(trap)
tl$ID<-trap$ID[match(tl$Traploc_ID, trap$Traploc_ID)]
K[tl$ID,4]<-2


### CHECK THAT NO RAT CAUGHT IN SESSION 8 WAS ATA TRAP WITH 0 nights
dim(rats[rats$TrapSession==8,])
length(which(y[,,4]!=0))
test<-melt(y[,,4], fun.aggregate=max)
names(test)<-c("ID","Trap_ID","capt")
testcast<-aggregate(capt~Trap_ID, test,FUN=max)
testcast$K<-K[,4]
testcast$test<-testcast$capt-testcast$K
K[,4]<-apply(testcast[,2:3],1,max)





################################################################################################################
########################## ARRANGE SIZE AND SEX INFORMATION ###################################
################################################################################################################

ALT<-aggregate(sexnum~ID, rats, FUN=min)
body<-aggregate(body~ID, rats, FUN=mean)

## standardize body size ##
body$stand<-(body$body-mean(body$body))/(sd(body$body))
hist(body$stand)

## merge the two

ALT<-merge(ALT, body, by="ID", all.x=T)
ALT$stand[is.na(ALT$stand)]<-0		### set rats with no body length info to 0

## create Animal Lookup Table (ALT) for JAGS TO FIND THE RIGHT INDEX OF SIGMA ###
ALT<-as.matrix(ALT[,c(4,2)])			
dimnames(ALT)<-NULL
dim(ALT)
dim(y)





################################################################################################################
########### SPECIFY INDIVIDUAL HOME RANGE CENTERS IN FIRST CAPTURE OCCASION ###################################
################################################################################################################

### set the upper and lower x and y coordinates for the state space
xl<- min(trap$Longitude)-150
yl<- min(trap$Latitude)-150
xu<- max(trap$Longitude)+150
yu<- max(trap$Latitude)+150
area<-((xu-xl)*(yu-yl))/10000


### create vector of first capture session and first Home-Range centre
first<-as.numeric()
xl<-as.numeric()
yl<-as.numeric()
xu<-as.numeric()
yu<-as.numeric()
for (i in 1:max(rats$ID)){
first[i]<-min(rats$TrapSession[rats$ID==i])-3
firstcaploc<-rats$Traploc_ID[rats$ID==i & rats$Date==min(rats$Date[rats$ID==i])]
xl[i]<- trap$Longitude[match(firstcaploc,trap$Traploc_ID)]-50
yl[i]<- trap$Latitude[match(firstcaploc,trap$Traploc_ID)]-50
xu[i]<- trap$Longitude[match(firstcaploc,trap$Traploc_ID)]+50
yu[i]<- trap$Latitude[match(firstcaploc,trap$Traploc_ID)]+50
}








################################################################################################################
########################## ARRANGE JAGS INPUT DATA ###################################
################################################################################################################


### Set up a data input for JAGS
CJS_data<-list(y=y, M=max(rats$ID), K=K, grid=grid, T=T,ntraps=ntraps, xl=xl, yl=yl, xu=xu, yu=yu, first=first, interval=interval, ttyp=ttyp, ALT=ALT, rain=raininput)

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Bayes_SECR")
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Bayes_SECR")

save.image("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Bayes_SECR\\spatialCJS_input_8occ_BB.RData")
save.image("C:\\STEFFEN\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Bayes_SECR\\spatialCJS_input_8occ_BB.RData")





