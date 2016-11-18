#######################################################################################
##  RAT DENSITY ESTIMATION ON HENDERSON ISLAND  - DATA PREPARATION FOR RD SECR MODEL ##
#######################################################################################
## written by steffen.oppel@rspb.org.uk on 30 May 2015
## updated on 11 Dez 2015, after incorporating Phase 2 data
## updated on 8 March 2016 to prepare data for RD SECR model of Ergon and Gardner (2014)
## updated on 4 April 2016: solved the problem of NA in J by filling with 1s
## updated on 4 April 2016: reference to 'tod' results in index of 0 for non-captured rats because it refers to H-1
## updated on 4 April 2016: 'tod' reference SOLVED BY CREATING NEW 'Htt2' MATRIX WITH ACTUAL TRAP NUMBERS
## updated 3 June 2016: fixed some database problem to delete NA in 'Htt2'
## noticed critical error in gr group assignment - now fixed

## CHANGED ACCESS QUERY TO EXCLUDE ALL TRAPS AND CAPTURES NOT ON PLATEAU

## UPDATE 25 OCT 2016: updated data preparation to account for changes in 'Sessions' table in Access database
## re-introduced Beachback data and split N into two groups (censored and normal).

## UPDATE 31 OCTOBER 2016: continued problem with convergence of sigma and survival estimates - removed BB data
## UPDATE 7 NOV 2016: included body size data, and new variable 'seen' based on Adam Butler's advice.
## UPDATED 14 NOV 2016 - removed rats with LARGE MOVEMENTS to stabilise sigma
## UPDATED 15 NOV 2016 - removed rats never recaptured


## DATA ANALYSIS MOVED TO SEPARATE SCRIPT

library(RODBC)
library(sp)
library(rgdal)
library(maptools)
library(reshape)


######################### INPUT DATA #############################################

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Henderson\\Data")
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Henderson\\Data")
SP<-odbcConnectAccess2007('HendersonData_2015.accdb')
rats <- sqlQuery(SP, "SELECT * FROM SECR_rat_input")
censoredrats <- sqlQuery(SP, "SELECT * FROM censoredrats")  
trap <- sqlQuery(SP, "SELECT * FROM SECR_trap_input")
sessions <- sqlQuery(SP, "SELECT * FROM Sessions")
traptypes <- sqlQuery(SP, "SELECT * FROM trap_deployments")
odbcClose(SP)

head(rats)
head(trap)
head(traptypes)


### REMOVE BEACHBACK FROM ALL DATA FRAMES ###
rats$area<-trap$area[match(rats$Traploc_ID, trap$Traploc_ID)]
trap<-trap[trap$area=="Plateau",]
rats<-rats[rats$area=="Plateau",]
sessions<-sessions[sessions$habitat=="Plateau",]
rats$area<-NULL
traptypes<-traptypes[traptypes$Traploc_ID %in% trap$Traploc_ID,]




##################### ARRANGE TRAP INPUT FILE #########################################################

trap<-trap[order(trap$Traploc_ID, decreasing=F),]
trap[is.na(trap)]<-0

### convert lat and long into UTM grid
trapSP<-SpatialPoints(trap[,4:3], proj4string=CRS("+proj=longlat +datum=WGS84"))
trapTransformed<-spTransform(trapSP, CRSobj =CRS("+proj=aeqd"))
trap[,4:3]<-trapTransformed@coords
trap$Longitude<- (trap$Longitude)*-1
trap$Longitude<- trap$Longitude-min(trap$Longitude)			### reduce the dimensions of the coordinates to avoid numerical problems
trap$Latitude<- trap$Latitude-min(trap$Latitude)			### reduce the dimensions of the coordinates to avoid numerical problems

### create unique number for each trap
trap$ID<-seq(1:dim(trap)[1])
head(trap)
grid<-as.matrix(trap[,4:3])
plot(grid)






########################## ARRANGE CAPTURE DATA ###################################

# remove rat captures without ID
rats<-rats[!is.na(rats$Rat_ID),]
dim(rats)

# update sex info
rats$Sex[is.na(rats$Sex)]<-'f'
rats$sexnum<-ifelse(rats$Sex=="m",1,2)

# remove all rats that were only captured in the last session, as they are useless for survival estimation
#firstcapt<-aggregate(TrapSession~Rat_ID, rats, FUN=min)
#exclude<-firstcapt$Rat_ID[firstcapt$TrapSession==7]
#rats<-rats[!(rats$Rat_ID %in% exclude),]
#dim(rats)

# sort data and use numeric RatID
rats<-rats[order(rats$Rat_ID, decreasing=F),]
ratids<-unique(rats$Rat_ID)
rats$ID<-match(rats$Rat_ID,ratids)

# match the trap IDs with numeric trap ID
rats$TrapID<-trap$ID[match(rats$Traploc_ID, trap$Traploc_ID)]

# create the data input matrix, which is an array that lists the number of captures of each rat at each trap in each session
rats$capt<-1




################################################################################################################
###################### REMOVE CRAZY MOVERS WITH FEW RECAPTURES  ################################################
################################################################################################################

nrecap<-aggregate(capt~Rat_ID, rats, sum)
totdists<-read.table("A:\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Movement_analysis\\Rat_Ind_move_dist_total.csv", sep=",", header=T)
nrecap<-merge(nrecap,totdists, by="Rat_ID", all.x=T)
nrecap$ratio<-nrecap$distance/nrecap$capt
head(nrecap)
plot(distance~capt, nrecap)
plot(ratio~capt, nrecap)

## remove rats with a recap distance >200 m
move200rats<-subset(nrecap,subset=ratio>200)
## remove rats with a total movement distance >4000 m
#move4krats<-subset(nrecap,subset=distance>4000)
sum(move200rats$capt)

dim(rats)
rats<-rats[!(rats$Rat_ID %in% move200rats$Rat_ID),]
dim(rats)



################################################################################################################
###################### REMOVE RATS NEVER RECAPTURED  ################################################
################################################################################################################
dim(nrecap)
length(nrecap$Rat_ID[nrecap$capt>1])

recaps<-nrecap$Rat_ID[nrecap$capt>1]

dim(rats)
rats<-rats[rats$Rat_ID %in% recaps,]
dim(rats)



###################################################################
### SPECIFIC DATA PREPARATION FOR RD SECR MODEL			    ###
### adopted from Ergon and Gardner				          ###
###################################################################
# DATA:
# - N[1]: Number of individuals that are only seen in the last primary session
#         or the session they are censored (these individuals must be sorted
#         first in the data)
# - N[2]:     Total number of individuals



### Determine rats that are immediately 'censored' (=killed) in the same primary session as they were marked

immcens<-censoredrats$Rat_ID[censoredrats$TrapSession==censoredrats$TrapSession.1]

### Determine rats that are marked in Trap Session 7

firstmarked<-aggregate(TrapSession~Rat_ID, rats, FUN=min)
latemark<-firstmarked$Rat_ID[firstmarked$TrapSession==7]
censrats<-c(immcens,latemark)


### ORDER THE NUMBER OF RATS AND GIVE THEM SEQUENTIAL NUMBERS
rats$sort<-ifelse(rats$Rat_ID %in% censrats,"first","second")
rats<-rats[order(rats$sort, rats$Rat_ID),]
ratids<-unique(rats$Rat_ID)
censrats<-ratids[ratids %in% censrats]
rats$ID<-match(rats$Rat_ID,ratids)
unique(rats$ID)

## for troubleshooting rat ID
rats[rats$ID==174,]

### COUNT THE NUMBER OF RAT INDIVIDUALS
N<-vector()
N[1]<-length(censrats)
N[2]<-length(ratids)



# - n.prim:   Number of primary sessions
n.prim<-max(sessions$TrapSession)
n.sec<-aggregate(secsess~TrapSession,sessions[sessions$habitat=="Plateau",],FUN=max)
n.sec<-n.sec$secsess

# - dt[k]:    Length of interval k (primary session)
sessmidtime<-aggregate(Date~TrapSession, sessions[sessions$habitat=="Plateau",], FUN=mean)
dt<-vector()
str(sessmidtime)
for (t in 1:(n.prim-1)){
dt[t]<-as.numeric(sessmidtime$Date[t+1]-sessmidtime$Date[t])
}

# - R:        Number of traps
R<-length(trap$Traploc_ID)

# - tod[k,j]: Category of trapping session k,j; two categories (e.g., time of day in original model, for us this is TRAP TYPE)
# FOR HENDERSON THIS VARIES BY TRAP LOCATION AND PRIMARY PERIOD - but not by secondary period
# SIMPLE MATRIX DOES NOT WORK BECAUSE INDEX from H WILL BE 0
# NEED COMPLEX ARRAY THAT MATCHES ARRAY FOR H

#tod<-array(NA, dim=c(length(trap$Traploc_ID),n.prim))
trap$Nr<-seq(2,R+1,1)
traptypes$Nr<-trap$Nr[match(traptypes$Traploc_ID,trap$Traploc_ID)]
traptypes$type<-ifelse(traptypes$TrapType=="Tomahawk",1,2)
tod.raw<-cast(traptypes, Nr~TrapSession, value='type')
#cast(traptypes, Traploc_ID~TrapSession, value='TrapType')			## for troubleshooting to see which trap type was missing
tod<-as.matrix(tod.raw[,2:8])


# - first[i]: First primary session of individual i
# may need to re-order individuals?
firstmarked<-aggregate(TrapSession~ID, rats, FUN=min)
first<-firstmarked$TrapSession

# - K[i]:     Last primary session for individual i (allows censoring)
lastseen<-aggregate(TrapSession~ID, rats, FUN=max)
K<-lastseen$TrapSession


# - J[i,k]:   Last secondary session for individual i in primary session k (allows censoring)
lastsecsess<-cast(rats, ID~TrapSession, value='secsess', fun.aggregate=max)
J<-as.matrix(lastsecsess[,c(2:8)])
J[J<0]<-NA


### counter J strikes an error if nothing was caught in a primary session
### need to fill in J with numbers between first and last primary session

for(i in 1:dim(J)[1]){
start<-first[i]
end<-K[i]
capts<-J[i,start:end]
capts[is.na(capts)]<-1			### replace all NA instances with 10
J[i,start:end]<-capts
}



# - gr[i]:    Group (sex) of individual i
rats$group<-ifelse(rats$Sex=="m",1,2)
#unique(rats$ID)
gr<-cast(rats, ID~., value='group', fun.aggregate=max)[,2]			## 
gr[is.na(gr)]<-2					## model will not run with uncertain sex assignment, so we arbitrarily assign all unknowns as males or females to check whether it affects output


# - size[i]:    Size (body length) of individual i
malesize<-mean(rats$body[rats$group==1], na.rm=T)
femalesize<-mean(rats$body[rats$group==2], na.rm=T)
rats$body<-ifelse(is.na(rats$body),ifelse(rats$group==1,malesize,femalesize),rats$body)
meanbody<-mean(rats$body, na.rm=T)
sdbody<-sd(rats$body, na.rm=T)
rats$bodystand<-(rats$body-meanbody)/sdbody
size<-cast(rats, ID~., value='bodystand', fun.aggregate=max)[,2]			## 




# - X[r,]:    Location of trap r = 1..R (coordinate)
X<-as.matrix(trap[,3:4])



# - H[i,j,k]: Trap index for capture of individual i in secondary session j within primary session k. 1 = not captured, other values are r+1
# requires a consecutive number for traps from 1:R
trap$Nr<-seq(2,R+1,1)
rats$TrapNr<-trap$Nr[match(rats$Traploc_ID, trap$Traploc_ID)]

H<-array(1, dim=c(N[2],max(J,na.rm=T),n.prim))
seen<-array(1, dim=c(N[2],max(J,na.rm=T),n.prim))

for(i in 1:N[2]){ 						# loop over each ind
ir<-subset(rats, ID==i)

     for(k in 1:K[i]){ 						# primary session
	irk<-subset(ir, TrapSession==k)

        for(j in 1:max(J,na.rm=T)){					# secondary session
		irkj<-subset(irk, secsess==j)
         H[i,j,k] <-ifelse(length(irkj$TrapNr)>0,irkj$TrapNr,1)
         seen[i,j,k] <-ifelse(length(irkj$TrapNr)>0,1,0)
        }		## close second sess
      }		## close prim sess
  } 			## close ind


### CHECK FOR NA
H[is.na(H)]



### HENDERSON ADD ON - instead of 'tod'-matrix we need a 3-dimensional array to specify the trap type for each location and session

Htt2<-array(1, dim=c(N[2],max(J,na.rm=T),n.prim))	# array for trap type for individuals in N[2]
for(i in 1:N[2]){ 						# loop over each ind
ir<-subset(rats, ID==i)
     for(k in 1:K[i]){ 						# primary session
	irk<-subset(ir, TrapSession==k)
        for(j in 1:max(J,na.rm=T)){					# secondary session
	   irkj<-subset(irk, secsess==j)
	   if (length(irkj$TrapNr)==1){
         Htt2[i,j,k] <-tod[irkj$TrapNr-1,k]}else{Htt2[i,j,k] <-1}
        }		## close second sess
      }		## close prim sess
  } 			## close ind


### CHECK NA
Htt2[is.na(Htt2)]<-1


# - Ones[i,j,k]: An array of all ones used in the "Bernoulli-trick" (see 
#                OpenBUGS user manual) in the observation likelihood. This saves
#                computation time as it is not necessary to compute the complete 
#                capture probabilities for every individual*trap*session (it is
#                sufficient to compute the g-variable for every
#                ind*trap*primary)

Ones<-array(1, dim=c(N[2],max(J,na.rm=T),n.prim))


# - xlow[i]: Lower bound in uniform prior for first x-location of individual i
# - xupp[i]: Upper bound in uniform prior for first x-location of individual i
# - ylow[i]: Lower bound in uniform prior for first y-location of individual i
# - yupp[i]: Upper bound in uniform prior for first y-location of individual i

xlow<-vector()
xupp<-vector()
ylow<-vector()
yupp<-vector()


for (i in 1:N[2]){
ir<-subset(rats, ID==i)
irf<-ir[ir$EncOcc==min(ir$EncOcc),18]			## should refer to 'TrapNr', was 13, changed on 7 Nov 2016

xlow[i]<-trap[trap$Nr==irf,4]-100				## x=longitude
ylow[i]<-trap[trap$Nr==irf,3]-100				## y=latitude
xupp[i]<-trap[trap$Nr==irf,4]+100				## x=longitude
yupp[i]<-trap[trap$Nr==irf,3]+100				## y=latitude
}



#### Vector of males and females
males<-which(x=(gr==1))
females<-which(x=(gr==2))


###################################################################
### MODEL FORMULATION FOR RD SECR MODEL			    ###
### adopted from Ergon and Gardner				          ###
###################################################################

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Bayes_SECR")
#setwd("C:\\STEFFEN\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Bayes_SECR")
save.image("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Bayes_SECR\\RD_CJS_SECR_input_size.RData")




###################################################################
### TROUBLE-SHOOTING MODEL FORMULATION	    ###
###################################################################

### counter J strikes an error if nothing was caught in a primary session
### need to fill in J with numbers between first and last primary session

for(i in 1:dim(J)[1]){

start<-first[i]
end<-K[i]
capts<-J[i,start:end]
capts[is.na(capts)]<-1			### replace all NA instances with 10
J[i,start:end]<-capts

}