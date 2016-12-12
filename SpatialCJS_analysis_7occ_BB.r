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

## ADJUSTED ON 7 DEC 2016 FOR BEACHBACK



###################################################################
### MODEL FORMULATION FOR SPATIAL CJS MODEL			    ###
### adopted from Ch 16.3 and Royle 2016					          ###
###################################################################

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Bayes_SECR")
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Bayes_SECR")




# Writing the model specification to a text file:
sink("Rat_spatialCJS_trans_moveHet_BB.jags")
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
	lam0[t]~dgamma(0.1,0.1)
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
           lam[i,j,first[i]] <- lam0[first[i]]*exp(-D2[i,j,first[i]]/(2*sigma2[first[i],ALT[i,2]]))		## exponential detection function based on distance
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
           lam[i,j,t] <- lam0[t]*exp(-D2[i,j,t]/(2*sigma2[t,ALT[i,2]]))		## exponential detection function based on distance
           tmp[i,j,t] <- lam[i,j,t]*K[j,t]*z[i,t]								## correct by number of nights this trap was actually available and by alive state of individual i at time t
           y[i,j,t] ~ dpois(tmp[i,j,t])										## count of observations at trap j is a poisson draw of 
				}	# close the trap loop for subsequent capture periods

### SURVIVAL MODEL FOR TIME PERIOD T ##
		phiUP[i,t]<-z[i,t-1]*phi[t-1,ALT[i,2]]				# change phi when not constant
		z[i,t]~dbern(phiUP[i,t])

		} 		## close the primary session loop			

          }		## close the individual loop


} ## close the model loop

",fill = TRUE)
sink()





library(reshape)
library(jagsUI)


##############################################
### LOADING WORKSPACE WITH PREPARED DATA   ###
##############################################
## data preparation in script "SpatialCJS_data_prep_7occ.r"

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Bayes_SECR")
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Bayes_SECR")
load("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Bayes_SECR\\spatialCJS_input_7occ_BB.RData")
load("C:\\STEFFEN\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Bayes_SECR\\spatialCJS_input_7occ_BB.RData")





##############################################
### VISUALISE PREPARED DATA   ###
##############################################
## this plots a 50 m buffer around each trap from where the home range center for each animal is drawn
plot(c(CJS_data$xl, CJS_data$xu), c(CJS_data$yl, CJS_data$yu))
points(CJS_data$grid, pch=16, cex=0.7)

plot(CJS_data$grid, pch=16, cex=0.7)







##################################################################
### FITTING THE MODEL IN JAGS ###
##################################################################

#Set the initial values
z<-matrix(NA,nrow=CJS_data$M,ncol=T)
for (i in 1:CJS_data$M){
z[i,first[i]]<- NA
if(first[i]<3){
	for(t in (first[i]+1):T){
 z[i,t]<- 1
}}
}

geninits =  function() {list(z=z,daily.phi=matrix(runif((T-1)*2,0.9,1),ncol=2), lam0=rgamma(T,0.1,0.1)) }

# Parameters to monitor
params = c('monthly.surv','lam0','sigma','sigma.ar')

# MCMC settings for the model run 
ni <- 50000		## number of iterations
nt <- 5		## thinning rate
nb <- 25000		## 'burn-in' - the number of iterations that are discarded at the start of each chain
nc <- 3		## number of chains

# Generating initial values for all chains
inits = list(geninits(),geninits(),geninits())

# Fitting the model in JAGS
ratmodel <- jags(CJS_data, inits, params, "S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Bayes_SECR\\Rat_spatialCJS_trans_moveHet_BB.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=T, n.cores=nc)


# Summary:
summary(ratmodel)

# Export Results
out<-as.data.frame(ratmodel$summary)
out$parameter<-row.names(ratmodel$summary)
write.table(out,"Henderson_spatialCJS_output_7occ_BB.csv", sep=",", row.names=F)






########################################################################################################################################
############### ARRANGING SUMMARY TABLE FOR PLOTTING ESTIMATES #################################################################
########################################################################################################################################
head(out)
dim(out)

surv.dat<-out[1:4,]			## annual survival
hr.dat<-out[8:13,]






##########################################################################################################
################ PLOT SURVIVAL PROBABILITY BETWEEN SESSIONS #####################################
##########################################################################################################

surv.dat$Sex<-rep(c('male','female'), each=2)
names(surv.dat)[c(3,7)]<-c("lcl","ucl")

surv.dat$session<-as.numeric(substr(surv.dat$parameter,14,14))			### change to 13 for constant survival
surv.dat$sessname<-"June"
surv.dat$sessname<-ifelse(surv.dat$session==1,"early Aug",surv.dat$sessname)
surv.dat$sessname<-ifelse(surv.dat$session==2,"late Aug",surv.dat$sessname)
plotF<-subset(surv.dat, Sex=="female")
plotF$x<-as.numeric(plotF$session)-0.2
plotM<-subset(surv.dat, Sex=="male")
plotM$x<-as.numeric(plotM$session)+0.2


pdf("S:\\ConSci\\DptShare\\SteffenOppel\\MANUSCRIPTS\\in_prep\\Henderson\\Rat_homerange\\Rat_surv_prob_BB.pdf", width=9, height=8)
par(mar=c(4,5,1,1), oma=c(1,1,0,0))
errbar(plotF$x,plotF$mean,plotF$lcl,plotF$ucl,type='p', las=1,ylab="Monthly survival probability", xlab="", main="", xlim=c(0.5,2.5), ylim=c(0,1), cex=1.5, cex.lab=1.5, axes=F, pch=16,col='darkred', frame=F)
par(new=T)
errbar(plotM$x,plotM$mean,plotM$lcl,plotM$ucl,type='p', axes=F, las=1,ylab="", xlab="", main="", xlim=c(0.5,2.5),ylim=c(0,1), cex=1.5, pch=16,col='navyblue', frame=F)
axis(1, at=c(1:2), labels=plotF$sessname, cex.axis=1.5, cex.lab=1.5)
axis(2, at=seq(0,1,0.2), las=1, cex.axis=1.5, cex.lab=1.5)
legend('bottomleft', pch=16, col=c('darkred','navyblue'),legend=c("female","male"), cex=1.4, bty='n')
dev.off()





##########################################################################################################
################ PLOT HOME RANGE SIZE BETWEEN SESSIONS #####################################
##########################################################################################################


hr.dat$Sex<-rep(c('male','female'),each=3)
names(hr.dat)[c(3,7)]<-c("lcl","ucl")

hr.dat$session<-as.numeric(substr(hr.dat$parameter,7,7))
hr.dat$sessname<-"June"
hr.dat$sessname<-ifelse(hr.dat$session==1,"early Aug",hr.dat$sessname)
hr.dat$sessname<-ifelse(hr.dat$session==2,"late Aug",hr.dat$sessname)
hr.dat$sessname<-ifelse(hr.dat$session==3,"Sept",hr.dat$sessname)


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
plotM$x<-as.numeric(plotM$session)+0.2



#pdf("A:\\MANUSCRIPTS\\in_prep\\Henderson\\Rat_homerange\\Rat_HR_size.pdf", width=11, height=8)
par(mar=c(4,5,1,1), oma=c(1,1,0,0))
errbar(plotF$x,plotF$HR_rad,plotF$HR_lcl,plotF$HR_ucl,type='p', las=1,ylab="Home range radius (m)", xlab="", main="", xlim=c(0.5,3.5), ylim=c(0,400), cex=1.5, cex.lab=1.5, axes=F, pch=16,col='darkred', frame=F)
par(new=T)
errbar(plotM$x,plotM$HR_rad,plotM$HR_lcl,plotM$HR_ucl,type='p', axes=F, las=1,ylab="", xlab="", main="", xlim=c(0.5,3.5),ylim=c(0,400), cex=1.5, pch=16,col='navyblue', frame=F)
axis(1, at=c(1:3), labels=plotF$sessname, cex.axis=1.5, cex.lab=1.5)
axis(2, at=seq(0,400,100), las=1, cex.axis=1.5, cex.lab=1.5)
legend('topleft', pch=16, col=c('darkred','navyblue'),legend=c("female","male"), cex=1.4, bty='n')
dev.off()






