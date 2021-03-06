#######################################################################################
##  RAT SURVIVAL ESTIMATION ON HENDERSON ISLAND  - RUNNING THE spatialCJS MODEL 	###
#######################################################################################

## model based on Ch 16.3 in Royle/Gardner/Sollmann book
## modified by steffen.oppel@rspb.org.uk on 22 July 2016
## agesex model run with all data was not very useful - unrealistic survival and home range estimates

## revised on 24 Oct 2016 - updated model to include time-dependent survival and Markovian transience
## survival estimates nearly constant and movement ranges are daft
## abandoned 24 Oct 2016 to use all 69 occasions rather than just 7 primary occasions

## revisited on 25 Nov 2016 because all other approaches using 69 occasions or robust-design failed

## revised 2 Dec because single sigma.ar does not converge "Rat_spatialCJS_trans_captHet.jags"
## 2 Dec: removed individual capture heterogeneity and introduced individual movement heterogeneity 

## FINALISED GRAPHS ON 5 DEC - MODEL CONVERGED!


library(reshape)
library(jagsUI)


##############################################
### LOADING WORKSPACE WITH PREPARED DATA   ###
##############################################
## data preparation in script "SpatialCJS_data_prep_7occ.r"

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Bayes_SECR")
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Bayes_SECR")
load("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Bayes_SECR\\spatialCJS_input_7occ.RData")
load("C:\\STEFFEN\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Bayes_SECR\\spatialCJS_input_7occ.RData")





##############################################
### VISUALISE PREPARED DATA   ###
##############################################
## this plots a 50 m buffer around each trap from where the home range center for each animal is drawn
plot(c(CJS_data$xl, CJS_data$xu), c(CJS_data$yl, CJS_data$yu))
points(CJS_data$grid, pch=16, cex=0.7)

plot(CJS_data$grid, pch=16, cex=0.7)









###################################################################
### MODEL FORMULATION FOR SPATIAL CJS MODEL			    ###
### adopted from Ch 16.3 and Royle 2016					          ###
###################################################################

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Bayes_SECR")
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Bayes_SECR")




# Writing the model specification to a text file:
sink("Rat_spatialCJS_trans_moveHet.jags")
cat("

model {


#### PRIORS AND CONSTRAINTS ####

## Movement around home-range centre affects capture probability ###
for(t in 1:T){						# allow different sigma for each session
	for (s in 1:2){					# allow different sigma for two sexes
	sigma[t,s] ~ dunif(5, 90)     		
	sigma2[t,s] <- sigma[t,s]*sigma[t,s]

	}
}



## Survival probability varies by sex and season - specified as annual survival###
for (s in 1:2){					# allow different phi for two sexes
	for(t in 1:(T-1)){
		daily.phi[t,s] ~ dunif(0.9,1)				### 0.997,0.9994
 		annual.surv[t,s]<-pow(daily.phi[t,s], 365)
 		monthly.surv[t,s]<-pow(daily.phi[t,s], 30)
    		phi[t,s] <- pow(daily.phi[t,s], interval[t])
	}
}




#### TRAP-DEPENDENT CAPTURE PROBABILITY ###

## Basic capture probability depends on trap type and season###
for(t in 1:T){
	for (u in 1:2){				## capture probability differs by traptype
	lam0[t,u]~dgamma(0.1,0.1)
	}
}




#### INDIVIDUAL MOVEMENT HETEROGENEITY ###

## Movement between primary sessions ###
for (i in 1:M){
	sigma.ar[i]~dunif(0,200)
	tau[i]<-1/(sigma.ar[i]*sigma.ar[i])
}



#### MODEL LIKELIHOOD ####  


### LOOP OVER ALL INDIVIDUALS ##     
   for (i in 1:M){    			# loop through the population of M trapped rats
       z[i,first[i]] <- 1       	# CJS is conditional on capture, so alive state in first capture period must be 1   
       SX[i,first[i]]~dunif(xl[i], xu[i])    	# priors for the activity centers for each individual at first capture
       SY[i,first[i]]~dunif(yl[i], yu[i])    	# xl is the lower x coordinate, xu is the upper x value

## OBSERVATION MODEL FOR FIRST SESSION ##
     for(j in 1:ntraps) {     	#loop through the J trap locations
           D2[i,j,first[i]] <- pow(SX[i,first[i]]-grid[j,1], 2) + pow(SY[i,first[i]]-grid[j,2],2)		## calculate distance from home range center to each trap location		
           lam[i,j,first[i]] <- lam0[first[i],ttyp[j,first[i]]]*exp(-D2[i,j,first[i]]/(2*sigma2[first[i],ALT[i,2]]))							## exponential detection function based on distance
           tmp[i,j,first[i]] <- lam[i,j,first[i]]*K[j,first[i]]									## correct by number of nights this trap was actually available
           y[i,j,first[i]] ~ dpois(tmp[i,j,first[i]])										## count of observations at trap j is a poisson draw of capture probability
				}	# close the trap loop for first capture period

## OBSERVATION MODEL FOR SUBSEQUENT SESSIONS ##       
	for(t in (first[i]+1):T){   			#loop through the primary periods after the first capture session
       	#SX[i,t]~dunif(xl[i], xu[i])    	# priors for the activity centers completely random
       	#SY[i,t]~dunif(yl[i], yu[i])    	# xl is the lower x coordinate, xu is the upper x value
       	SX[i,t]~dnorm(SX[i,t-1], tau[i])#T(xl[i], xu[i])    	# tau[t,ALT[i,2]]; priors for the activity centers based on movement from first capture location, truncated by state space
       	SY[i,t]~dnorm(SY[i,t-1], tau[i])#T(yl[i], yu[i])   	# tau[t,ALT[i,2]]; xl is the lower x coordinate, xu is the upper x value

	#HR_shift[i,t] <- pow(pow(SX[i,t]-SX[i,t-1], 2) + pow(SY[i,t]-SY[i,t-1],2),0.5)	## shift of home range centre from one prim occ to next

     for(j in 1:ntraps) {     	#loop through the J trap locations
           D2[i,j,t] <- pow(pow(SX[i,t]-grid[j,1], 2) + pow(SY[i,t]-grid[j,2],2),0.5)		## calculate distance from home range center to each trap location		
           lam[i,j,t] <- lam0[t,ttyp[j,t]]*exp(-D2[i,j,t]/(2*sigma2[t,ALT[i,2]]))		## exponential detection function based on distance
           tmp[i,j,t] <- lam[i,j,t]*K[j,t]*z[i,t]								## correct by number of nights this trap was actually available and by alive state of individual i at time t
           y[i,j,t] ~ dpois(tmp[i,j,t])										## count of observations at trap j is a poisson draw of 
				}	# close the trap loop for subsequent capture periods

### SURVIVAL MODEL FOR TIME PERIOD T ##
		phiUP[i,t]<-z[i,t-1]*phi[t-1,ALT[i,2]]				# change phi when not constant
		z[i,t]~dbern(phiUP[i,t])

		} 		## close the primary session loop			

          }		## close the individual loop



### DERIVED PARAMETERS FOR SURVIVAL AND CAPTURE PROBABILITY ###

capt.prob.Sherman<- mean(lam0[1:6,1])
capt.prob.Tomahawk<- mean(lam0[6:7,2])



} ## close the model loop

",fill = TRUE)
sink()










##################################################################
### FITTING THE MODEL IN JAGS ###
##################################################################

#Set the initial values
z<-matrix(NA,nrow=CJS_data$M,ncol=T)
for (i in 1:CJS_data$M){
z[i,first[i]]<- NA
if(first[i]<7){
	for(t in (first[i]+1):T){
 z[i,t]<- 1
}}
}

geninits =  function() {list(z=z,daily.phi=matrix(runif((T-1)*2,0.9,1),ncol=2), lam0=matrix(runif(T*2,0,1),ncol=2)) }
geninits =  function() {list(z=z,daily.phi=matrix(runif((T-1)*2,0.9,1),ncol=2), lam0=matrix(rgamma(T*2,0.1,0.1),ncol=2)) }

# Parameters to monitor
params = c('monthly.surv','capt.prob.Sherman','capt.prob.Tomahawk','sigma','sigma.ar')

# MCMC settings for the model run 
ni <- 50000		## number of iterations
nt <- 5		## thinning rate
nb <- 25000		## 'burn-in' - the number of iterations that are discarded at the start of each chain
nc <- 3		## number of chains

# Generating initial values for all chains
inits = list(geninits(),geninits(),geninits())

# Fitting the model in JAGS
ratmodel <- jags(CJS_data, inits, params, "A:\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Bayes_SECR\\Rat_spatialCJS_trans_moveHet.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=T, n.cores=nc)


# Summary:
summary(ratmodel)

# Export Results
out<-as.data.frame(ratmodel$summary)
out$parameter<-row.names(ratmodel$summary)
write.table(out,"Henderson_spatialCJS_output_7occ.csv", sep=",", row.names=F)

# TRouble-shoot model:
traceplot(ratmodel, c("sigma","capt.prob.Sherman","capt.prob.Tomahawk","sigma.ar"))






########################################################################################################################################
############### ARRANGING SUMMARY TABLE FOR PLOTTING ESTIMATES #################################################################
########################################################################################################################################
head(out)
dim(out)

surv.dat<-out[1:12,]			## annual survival
hr.dat<-out[15:28,]
move.dat<-out[29:552,]


### add estimates from BB:
outBB<-read.table("Henderson_spatialCJS_output_7occ_BB.csv", sep=",", header=T)
surv.datBB<-outBB[1:4,]			## annual survival
hr.datBB<-outBB[8:13,]






##########################################################################################################
################ PLOT SURVIVAL PROBABILITY BETWEEN SESSIONS #####################################
##########################################################################################################

surv.dat$Sex<-rep(c('male','female'), each=6)
names(surv.dat)[c(3,7)]<-c("lcl","ucl")

surv.dat$session<-as.numeric(substr(surv.dat$parameter,14,14))			### change to 13 for constant survival
surv.dat$sessname<-"June"
surv.dat$sessname<-ifelse(surv.dat$session==1,"June",surv.dat$sessname)
surv.dat$sessname<-ifelse(surv.dat$session==2,"early July",surv.dat$sessname)
surv.dat$sessname<-ifelse(surv.dat$session==3,"late July",surv.dat$sessname)
surv.dat$sessname<-ifelse(surv.dat$session==4,"early Aug",surv.dat$sessname)
surv.dat$sessname<-ifelse(surv.dat$session==5,"late Aug",surv.dat$sessname)
surv.dat$sessname<-ifelse(surv.dat$session==6,"Sept",surv.dat$sessname)
plotF<-subset(surv.dat, Sex=="female")
plotF$x<-as.numeric(plotF$session)-0.2
plotM<-subset(surv.dat, Sex=="male")
plotM$x<-as.numeric(plotM$session)+0.2


#### ADD BB DATA ###
surv.datBB$Sex<-rep(c('male','female'), each=2)
names(surv.datBB)[c(3,7)]<-c("lcl","ucl")

surv.datBB$session<-as.numeric(substr(surv.datBB$parameter,14,14))			### change to 13 for constant survival
surv.datBB$sessname<-"June"
surv.datBB$sessname<-ifelse(surv.datBB$session==1,"early Aug",surv.datBB$sessname)
surv.datBB$sessname<-ifelse(surv.datBB$session==2,"late Aug",surv.datBB$sessname)
plotFBB<-subset(surv.datBB, Sex=="female")
plotFBB$x<-as.numeric(plotFBB$session)+2.9
plotMBB<-subset(surv.datBB, Sex=="male")
plotMBB$x<-as.numeric(plotMBB$session)+3.3




pdf("A:\\MANUSCRIPTS\\in_prep\\Henderson\\Rat_homerange\\Rat_surv_prob.pdf", width=11, height=8)
par(mar=c(4,5,1,1), oma=c(1,1,0,0))
errbar(plotF$x,plotF$mean,plotF$lcl,plotF$ucl,type='p', las=1,ylab="Monthly survival probability", xlab="", main="", xlim=c(0.5,6.5), ylim=c(0,1), cex=1.5, cex.lab=1.5, axes=F, pch=16,col='darkred', frame=F)
par(new=T)
errbar(plotM$x,plotM$mean,plotM$lcl,plotM$ucl,type='p', axes=F, las=1,ylab="", xlab="", main="", xlim=c(0.5,6.5),ylim=c(0,1), cex=1.5, pch=16,col='navyblue', frame=F)
par(new=T)
errbar(plotFBB$x,plotFBB$mean,plotFBB$lcl,plotFBB$ucl,type='p', las=1,ylab="", xlab="", main="", xlim=c(0.5,6.5), ylim=c(0,1), cex=1.5, cex.lab=1.5, axes=F, pch=1,col='darkred', frame=F)
par(new=T)
errbar(plotMBB$x,plotMBB$mean,plotMBB$lcl,plotMBB$ucl,type='p', axes=F, las=1,ylab="", xlab="", main="", xlim=c(0.5,6.5),ylim=c(0,1), cex=1.5, pch=1,col='navyblue', frame=F)

axis(1, at=c(1:6), labels=plotF$sessname, cex.axis=1.5, cex.lab=1.5)
axis(2, at=seq(0,1,0.2), las=1, cex.axis=1.5, cex.lab=1.5)
legend('bottomleft', pch=16, col=c('darkred','navyblue'),legend=c("female","male"), cex=1.4, bty='n')
dev.off()





##########################################################################################################
################ PLOT HOME RANGE RADIUS BETWEEN SESSIONS #####################################
##########################################################################################################


hr.dat$Sex<-rep(c('male','female'),each=7)
names(hr.dat)[c(3,7)]<-c("lcl","ucl")

hr.dat$session<-as.numeric(substr(hr.dat$parameter,7,7))
hr.dat$sessname<-"June"
hr.dat$sessname<-ifelse(hr.dat$session==1,"June",hr.dat$sessname)
hr.dat$sessname<-ifelse(hr.dat$session==2,"early July",hr.dat$sessname)
hr.dat$sessname<-ifelse(hr.dat$session==3,"late July",hr.dat$sessname)
hr.dat$sessname<-ifelse(hr.dat$session==4,"early Aug",hr.dat$sessname)
hr.dat$sessname<-ifelse(hr.dat$session==5,"late Aug",hr.dat$sessname)
hr.dat$sessname<-ifelse(hr.dat$session==6,"Sept",hr.dat$sessname)
hr.dat$sessname<-ifelse(hr.dat$session==7,"Oct",hr.dat$sessname)

### CALCULATE HOME RANGE RADIUS
library(secr)
for (l in 1:dim(hr.dat)[1]){
hr.dat$HR_rad[l]<-circular.r (p = 0.95, detectfn = 'EX', sigma = hr.dat$mean[l], hazard = TRUE)		### this is the mean home-range radius of rats in m
hr.dat$HR_lcl[l]<-circular.r (p = 0.95, detectfn = 'EX', sigma = hr.dat$lcl[l], hazard = TRUE)		### this is the mean home-range radius of rats in m
hr.dat$HR_ucl[l]<-circular.r (p = 0.95, detectfn = 'EX', sigma = hr.dat$ucl[l], hazard = TRUE)		### this is the mean home-range radius of rats in m
}


plotF<-subset(hr.dat, Sex=="female")
plotF$x<-as.numeric(plotF$session)-0.2
plotM<-subset(hr.dat, Sex=="male")
plotM$x<-as.numeric(plotM$session)+0.1


### ADD BEACHBACK DATA

hr.datBB$Sex<-rep(c('male','female'),each=3)
names(hr.datBB)[c(3,7)]<-c("lcl","ucl")

hr.datBB$session<-as.numeric(substr(hr.datBB$parameter,7,7))
hr.datBB$sessname<-"June"
hr.datBB$sessname<-ifelse(hr.datBB$session==1,"early Aug",hr.datBB$sessname)
hr.datBB$sessname<-ifelse(hr.datBB$session==2,"late Aug",hr.datBB$sessname)
hr.datBB$sessname<-ifelse(hr.datBB$session==3,"Sept",hr.datBB$sessname)


### CALCULATE HOME RANGE RADIUS
library(secr)
for (l in 1:dim(hr.datBB)[1]){
hr.datBB$HR_rad[l]<-circular.r (p = 0.95, detectfn = 'EX', sigma = hr.datBB$mean[l], hazard = TRUE)		### this is the mean home-range radius of rats in m
hr.datBB$HR_lcl[l]<-circular.r (p = 0.95, detectfn = 'EX', sigma = hr.datBB$lcl[l], hazard = TRUE)		### this is the mean home-range radius of rats in m
hr.datBB$HR_ucl[l]<-circular.r (p = 0.95, detectfn = 'EX', sigma = hr.datBB$ucl[l], hazard = TRUE)		### this is the mean home-range radius of rats in m
}


plotFBB<-subset(hr.datBB, Sex=="female")
plotFBB$x<-as.numeric(plotFBB$session)+2.9
plotMBB<-subset(hr.datBB, Sex=="male")
plotMBB$x<-as.numeric(plotMBB$session)+3.2





#pdf("A:\\MANUSCRIPTS\\in_prep\\Henderson\\Rat_homerange\\Rat_HR_size.pdf", width=11, height=8)
par(mar=c(2,5,1,1), oma=c(1,1,0,0), mgp=c(4,1,0))
errbar(plotF$x,plotF$HR_rad,plotF$HR_lcl,plotF$HR_ucl,type='p', las=1,ylab="Home range radius (m)", xlab="", main="", xlim=c(0.5,7.5), ylim=c(0,400), cex=1.5, cex.lab=1.5, axes=F, pch=16,col='darkred', frame=F)
par(new=T)
errbar(plotM$x,plotM$HR_rad,plotM$HR_lcl,plotM$HR_ucl,type='p', axes=F, las=1,ylab="", xlab="", main="", xlim=c(0.5,7.5),ylim=c(0,400), cex=1.5, pch=16,col='navyblue', frame=F)
par(new=T)
errbar(plotMBB$x,plotMBB$HR_rad,plotMBB$HR_lcl,plotMBB$HR_ucl,type='p', axes=F, las=1,ylab="", xlab="", main="", xlim=c(0.5,7.5),ylim=c(0,400), cex=1.5, pch=1,col='navyblue', frame=F)
par(new=T)
errbar(plotFBB$x,plotFBB$HR_rad,plotFBB$HR_lcl,plotFBB$HR_ucl,type='p', las=1,ylab="", xlab="", main="", xlim=c(0.5,7.5), ylim=c(0,400), cex=1.5, cex.lab=1.5, axes=F, pch=1,col='darkred', frame=F)
axis(1, at=c(1:7), labels=plotF$sessname, cex.axis=1.5, cex.lab=1.5)
axis(2, at=seq(0,400,100), las=1, cex.axis=1.5, cex.lab=1.5)
legend('topleft', pch=16, col=c('darkred','navyblue'),legend=c("female","male"), cex=1.4, bty='n')
dev.off()




### CALCULATE HOME RANGE AREA
hr.dat$AREA<-((hr.dat$HR_rad^2)*pi)/10000
hr.datBB$AREA<-((hr.datBB$HR_rad^2)*pi)/10000

hr.dat$AREA_lcl<-((hr.dat$HR_lcl^2)*pi)/10000
hr.datBB$AREA_lcl<-((hr.datBB$HR_lcl^2)*pi)/10000


### CALCULATE BAIT IN MIN HOME RANGE
platBait<-10000			## g/ha
beachBait<-40000			## g/ha
pellet<-2.5				## g for one bait pellet

platBaitDensity<-platBait/pellet
beachBaitDensity<-beachBait/pellet

min(hr.dat$AREA_lcl)*platBaitDensity
min(hr.datBB$AREA_lcl)*beachBaitDensity






##########################################################################################################
################ PLOT HOME RANGE SHIFT BETWEEN SESSIONS #####################################
##########################################################################################################
move.dat$Sex<-ifelse(ALT[,2]==1,'male','female')
names(move.dat)[c(3,7)]<-c("lcl","ucl")
head(move.dat)





## HOME RANGE CENTRE SHIFT ##

#pdf("A:\\MANUSCRIPTS\\in_prep\\Henderson\\Rat_homerange\\Rat_HR_shift.pdf", width=10, height=15)
par(mfrow=c(2,1),mar=c(4,5,1,1), oma=c(1,1,0,0))
hist(move.dat$mean[move.dat$Sex=="female"],30, col='darkred', xlim=c(0,200),main="", xlab="", ylab="N females", cex.axis=1.5, cex.lab=1.5, las=1)
hist(move.dat$mean[move.dat$Sex=="male"],30, col='navyblue', xlim=c(0,200),main="", xlab="Home range centre shift (m)", ylab="N males", cex.axis=1.5, cex.lab=1.5, las=1)
dev.off()






## for session-specific shifts
#move.dat$ncar<-nchar(move.dat$parameter)
#move.dat$ID<-as.numeric(substr(move.dat$parameter,10,(move.dat$ncar-3)))
#move.dat$session<-substr(move.dat$parameter,(move.dat$ncar-1),(move.dat$ncar-1))
plotdat<-aggregate(mean~Sex+session, move.dat, FUN=mean)
plotdat$ucl<-aggregate(ucl~Sex+session, move.dat, FUN=mean)[,3]
plotdat$lcl<-aggregate(lcl~Sex+session, move.dat, FUN=mean)[,3]
plotdat$session<-as.numeric(plotdat$session)
plotdat$sessname<-"June"
plotdat$sessname<-ifelse(plotdat$session==3,"early July",plotdat$sessname)
plotdat$sessname<-ifelse(plotdat$session==4,"late July",plotdat$sessname)
plotdat$sessname<-ifelse(plotdat$session==5,"early Aug",plotdat$sessname)
plotdat$sessname<-ifelse(plotdat$session==6,"late Aug",plotdat$sessname)
plotdat$sessname<-ifelse(plotdat$session==7,"Sept",plotdat$sessname)
plotF<-subset(plotdat, Sex=="female")
plotF$x<-as.numeric(plotF$session)-0.2
plotM<-subset(plotdat, Sex=="male")
plotM$x<-as.numeric(plotM$session)+0.2


## HOME RANGE CENTRE SHIFT ##

par(mar=c(4,5,1,1), oma=c(1,1,0,0))
errbar(plotF$x,plotF$mean,plotF$lcl,plotF$ucl,type='p', las=1,ylab="Home range centre shift (m)", xlab="", main="", xlim=c(1.5,7.5), ylim=c(0,500), cex=1.5, cex.lab=1.5, axes=F, pch=16,col='darkred', frame=F)
par(new=T)
errbar(plotM$x,plotM$mean,plotM$lcl,plotM$ucl,type='p', axes=F, las=1,ylab="", xlab="", main="", xlim=c(1.5,7.5),ylim=c(0,500), cex=1.5, pch=16,col='navyblue', frame=F)
axis(1, at=c(2:7), labels=plotF$sessname, cex.axis=1.5, cex.lab=1.5)
axis(2, at=seq(0,500,100), las=1, cex.axis=1.5, cex.lab=1.5)
legend('topleft', pch=16, col=c('darkred','navyblue'),legend=c("female","male"), cex=1.4, bty='n')










##########################################################################################################
################ CORRELATION BETWEEN SURVIVAL AND MOVEMENT PARAMETERS #####################################
##########################################################################################################

str(ratmodel$sims.list$sigma)
str(ratmodel$sims.list$monthly.surv)

correlations<-data.frame(session=rep(c(1:6), each=2), sex=rep(c(1:2),6), sigma=0,surv=0)
#par(mfrow=c(6,2))
for (sess in 1:6){
for (sex in 1:2){
title<-cor.test(ratmodel$sims.list$sigma[,sess,sex],ratmodel$sims.list$monthly.surv[,sess,sex])
#plot(ratmodel$sims.list$sigma[,sess,sex]~ratmodel$sims.list$monthly.surv[,sess,sex], pch=16, cex=0.5, xlab="", ylab="")
correlations$sigma[correlations$sex==sex & correlations$session==sess]<-mean(ratmodel$sims.list$sigma[,sess,sex])
correlations$surv[correlations$sex==sex & correlations$session==sess]<-mean(ratmodel$sims.list$monthly.surv[,sess,sex])
}
}

plot(sigma~surv, correlations)
cor.test(correlations$sigma,correlations$surv, method='spearman')
cor.test(correlations$sigma[correlations$sex==1],correlations$surv[correlations$sex==1], method='spearman')
cor.test(correlations$sigma[correlations$sex==2],correlations$surv[correlations$sex==2], method='spearman')