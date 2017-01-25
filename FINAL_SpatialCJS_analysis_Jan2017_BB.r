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
## ADJUSTED ON 7 DEC 2016 FOR SIMPLER MODEL WITHOUT SEX DIFFERENCES
## UPDATED on 4 Jan 2017 to include rain covariate effect on survival and capture prob
## DOES NOT WORK WITH COVARIATES ON CAPTURE PROBABILITY
## UPDATED on 6 Jan 2017 to use mixture model for sigma.ar

## UPDATED 9 Jan to use mixture for sigma, and revert to simple sigma.ar

## NEW DATA ON 11 JAN 2015 - added snap trapping as session 8 and trap type 3

## UPDATED 16 JAN 2017 to re-insert 3 trap types and fix K input vector

## UPDATED 25 JAN 2017 TO ADAPT FOR BEACHBACK



###################################################################
### MODEL FORMULATION FOR SPATIAL CJS MODEL			    ###
###################################################################

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Bayes_SECR")
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Bayes_SECR")




# Writing the model specification to a text file:
sink("Rat_spatialCJS_FINAL_BB.jags")
cat("

model {


#### PRIORS AND CONSTRAINTS ####

## Movement around home-range centre affects capture probability ###
## INCLUDED two part mixture for the normal rats and the long-moving ones

prop.longmoves ~ dbeta(0.5,2)			## this sets up a prior for drawing the sigma 

for (i in 1:M){
	move.selector[i] ~ dbern(prop.longmoves)		## draws a random number of 0 or 1, reflecting the two movement distributions

		for(t in first[i]:T){						# allow different sigma for each session
			#for (s in 1:2){					# allow different sigma for two sexes
			#sigma[t] ~ dunif(0, 200)     		

		# prior for very long movements
		sigma[i,t,2]~dunif(30,400)				## index 2 is for the crazy moving rats
		# prior for normal movements
		sigma[i,t,1]~dnorm(20,0.001)				## index 1 is for the normal moving rats - what we expect with sigma 20-40

		sigma2[i,t] <- (sigma[i,t,move.selector[i]+1]*sigma[i,t,move.selector[i]+1])
		}
}



## Survival probability varies by sex and season - specified as annual survival###
	for(t in 1:(T-1)){
		daily.phi[t] ~ dunif(0.9,1)				### 0.997,0.9994
 		annual.surv[t]<-pow(daily.phi[t], 365)
 		monthly.surv[t]<-pow(daily.phi[t], 30)
    		phi[t] <- pow(daily.phi[t], interval[t])
  		lphi[t] <- log(phi[t]/(1-phi[t]))			### transform to logit scale to allow inclusion of continuous covariates
	}





#### TRAP-DEPENDENT CAPTURE PROBABILITY ###

## Basic capture probability depends on trap type and season###
for(t in 1:T){
	for(f in 1:2){					# trap types Sherman, Tomahawk and snap trap
		lam0[t,f]~dgamma(0.1,0.1)		# for poisson draw
	}
}



#### INDIVIDUAL MOVEMENT HETEROGENEITY ###
## rather than define per individual, make this a mixture

## Movement between primary sessions ###
for (i in 1:M){
sigma.ar[i]~dunif(0,200)
tau[i]<-1/(sigma.ar[i]*sigma.ar[i])
}


#### COVARIATE EFFECTS #####

b.sex.surv ~ dnorm(0,1)
b.rain.surv ~ dnorm(0,1)



#### MODEL LIKELIHOOD ####  


### LOOP OVER ALL INDIVIDUALS ##     
   for (i in 1:M){    			# loop through the population of M trapped rats
       z[i,first[i]] <- 1       	# CJS is conditional on capture, so alive state in first capture period must be 1   
       SX[i,first[i]]~dunif(xl[i], xu[i])    	# priors for the activity centers for each individual at first capture
       SY[i,first[i]]~dunif(yl[i], yu[i])    	# xl is the lower x coordinate, xu is the upper x value

## OBSERVATION MODEL FOR FIRST SESSION ##
     for(j in 1:ntraps) {     	#loop through the J trap locations
           D2[i,j,first[i]] <- pow(SX[i,first[i]]-grid[j,1], 2) + pow(SY[i,first[i]]-grid[j,2],2)		## calculate distance from home range center to each trap location		
           lam[i,j,first[i]] <- lam0[first[i],ttyp[j,first[i]]]*exp(-D2[i,j,first[i]]/(2*sigma2[i,first[i]]))		## exponential detection function based on distance
           #log(lam[i,j,first[i]]) <- (lp[first[i]]*(-D2[i,j,first[i]]/(2*sigma2[i,first[i]])))+ 1/sqrt(tau.capt)*eta[i] + b.rain.capt*rain[first[i]] + b.sex.capt*ALT[i,2]		## exponential detection function based on distance
           tmp[i,j,first[i]] <- lam[i,j,first[i]]*K[j,first[i]]									## correct by number of nights this trap was actually available
           y[i,j,first[i]] ~ dpois(tmp[i,j,first[i]])										## count of observations at trap j is a poisson draw of capture probability
				}	# close the trap loop for first capture period
### DERIVED PARAMETERS FOR HOME RANGE IN FIRST PERIOD ##
		hrradius[i,first[i]]<-pow(sigma2[i,first[i]],0.5)						## home range radius as 'sigma' for conversion by secr

## OBSERVATION MODEL FOR SUBSEQUENT SESSIONS ##       
	for(t in (first[i]+1):T){   			#loop through the primary periods after the first capture session
       	SX[i,t]~dnorm(SX[i,t-1], tau[i])#T(xl[i], xu[i])    	# tau[t,ALT[i,2]]; priors for the activity centers based on movement from first capture location, truncated by state space
       	SY[i,t]~dnorm(SY[i,t-1], tau[i])#T(yl[i], yu[i])   	# tau[t,ALT[i,2]]; xl is the lower x coordinate, xu is the upper x value

     for(j in 1:ntraps) {     	#loop through the J trap locations
           D2[i,j,t] <- pow(pow(SX[i,t]-grid[j,1], 2) + pow(SY[i,t]-grid[j,2],2),0.5)		## calculate distance from home range center to each trap location		
           lam[i,j,t] <- lam0[t,ttyp[j,t]]*exp(-D2[i,j,t]/(2*sigma2[i,t]))		## exponential detection function based on distance
           tmp[i,j,t] <- lam[i,j,t]*K[j,t]*z[i,t]								## correct by number of nights this trap was actually available and by alive state of individual i at time t
           y[i,j,t] ~ dpois(tmp[i,j,t])										## count of observations at trap j is a poisson draw of 
				}	# close the trap loop for subsequent capture periods

### SURVIVAL MODEL FOR TIME PERIOD T ##
    		logit(PHI[i,t-1]) <- lphi[t-1] + b.sex.surv*ALT[i,2] + b.rain.surv*rain[t-1]
		phiUP[i,t]<-z[i,t-1]*PHI[i,t-1]				# changed from phi phi when not constant
		z[i,t]~dbern(phiUP[i,t])

### DERIVED PARAMETERS FOR HOME RANGE AND MOVEMENTS ##
		hrradius[i,t]<-z[i,t]*pow(sigma2[i,t],0.5)						## home range radius as 'sigma' for conversion by secr
		hrshift[i,t] <- pow(pow(SX[i,t]-SX[i,t-1], 2) + pow(SY[i,t]-SY[i,t-1],2),0.5)	## shift of home range centre from one prim occ to next

		} 		## close the primary session loop			

          }		## close the individual loop


} ## close the model loop

",fill = TRUE)
sink()





library(reshape)
library(jagsUI)
library(secr)
library(plyr)

##############################################
### LOADING WORKSPACE WITH PREPARED DATA   ###
##############################################
## data preparation in script "SpatialCJS_data_prep_8occ.r"

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Bayes_SECR")
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Bayes_SECR")
load("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Bayes_SECR\\spatialCJS_input_8occ_BB.RData")
load("C:\\STEFFEN\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Bayes_SECR\\spatialCJS_input_8occ_BB.RData")




##############################################
### VISUALISE PREPARED DATA   ###
##############################################
## this plots a 50 m buffer around each trap from where the home range center for each animal is drawn
plot(c(CJS_data$xl, CJS_data$xu), c(CJS_data$yl, CJS_data$yu))
points(CJS_data$grid, pch=16, cex=0.7)







##################################################################
### FITTING THE MODEL IN JAGS ###
##################################################################

#Set the initial values
z<-matrix(NA,nrow=CJS_data$M,ncol=T)
for (i in 1:CJS_data$M){
z[i,first[i]]<- NA
if(first[i]<4){
	for(t in (first[i]+1):T){
 z[i,t]<- 1
}}
}

geninits =  function() {list(z=z,daily.phi=runif((T-1),0.9,1), lam0=matrix(rgamma(T*2,0.1,0.1),ncol=max(ttyp))) }

# Parameters to monitor
params = c('monthly.surv','b.sex.surv','b.rain.surv','prop.longmoves','move.selector','hrradius','hrshift')


# MCMC settings for the model run 
ni <- 30000		## number of iterations
nt <- 1		## thinning rate
nb <- 2000		## 'burn-in' - the number of iterations that are discarded at the start of each chain
nc <- 3		## number of chains
na <- 5000		## adaptation samples, additional to burnin, default is 100

# Generating initial values for all chains
inits = list(geninits(),geninits(),geninits())


# Fitting the model in JAGS
ratmodel <- jags(CJS_data, inits, params, "S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Bayes_SECR\\Rat_spatialCJS_FINAL_BB.jags", n.chains = nc, n.thin = nt, n.adapt=na, n.iter = ni, n.burnin = nb, parallel=T, n.cores=nc)


# Summary:
summary(ratmodel)

# Export Results
out<-as.data.frame(ratmodel$summary)
out$parameter<-row.names(ratmodel$summary)
write.table(out,"Henderson_spatialCJS_output_8occ_FINAL_BB.csv", sep=",", row.names=F)
dim(out)




##########################################################################################################
################ DIAGNOSTICS FOR CONVERGENCE #####################################
##########################################################################################################
hist(out$Rhat,30)
unique(substr(out$parameter[(out$Rhat>1.1)],1,7))
traceplot(ratmodel)




##########################################################################################################
################ SAVING THE WORKSPACE #####################################
##########################################################################################################

OUT<-ratmodel
OUT$sims.list<-NULL
OUT$mcmc.info<-NULL
OUT$model<-NULL
OUT$samples<-NULL
OUT$f<-NULL
OUT$Rhat<-NULL
OUT$n.eff<-NULL
OUT$overlap0<-NULL
rm(ratmodel)




########################################################################################################################################
############### EXTRACTING INDIVIDUAL INFO ON MOVEMENT  #################################################################
########################################################################################################################################
head(out)
dim(out)
move.sel<-out[grep("move.selector",out$parameter),]
head(move.sel)
ind.move<-as.data.frame(CJS_data$ALT)
names(ind.move)<-c('body','sex')
ind.move$crazy<-as.numeric(round(move.sel$mean))
table(ind.move$sex,ind.move$crazy)
mean(ind.move$crazy)
aggregate(crazy~sex, ind.move,mean)
ind.move$ID<-seq(1,60,1)


########################################################################################################################################
############### EXTRACTING PARAMETER ESTIMATES FOR SURVIVAL COVARIATES  #################################################################
########################################################################################################################################

OUT$mean$b.sex.surv
OUT$q2.5$b.sex.surv
OUT$q97.5$b.sex.surv

OUT$mean$b.rain.surv
OUT$q2.5$b.rain.surv
OUT$q97.5$b.rain.surv










##########################################################################################################
################ PLOT HOME RANGE SIZE BETWEEN SESSIONS #####################################
##########################################################################################################

######### CREATE MATRIX OF ACTUAL DETECTIONS ###########

det<-OUT$mean$hrradius
det[,]<-0
dim(det)
for (p in 1:T){
det[,p]<-apply(y[,,p],1,max)
}
#det[det>1]<-1
pick<-melt(det)
names(pick)[1:3]<-c("ID","session","SELECT")
pick$parameter<-paste("hrradius[",pick[,1],",",pick[,2],"]",sep="")

hist(pick$SELECT)
hist(OUT$mean$hrradius,50)
hr.dat<-pick[pick$SELECT>0,]
hr.dat<-hr.dat[hr.dat$session<8,]			### no data on the 8th trapping session because only kill traps used


######### CREATE MATRIX OF ACTUAL DETECTIONS ###########
hr.dat$hr_radius<-0
hr.dat$lcl<-0
hr.dat$ucl<-0

for (i in 1:dim(hr.dat)[1]){
		
			hr.dat$hr_radius[i]<-circular.r (p = 0.95, detectfn = 'EX', sigma = OUT$mean$hrradius[hr.dat$ID[i],hr.dat$session[i]], hazard = TRUE)		### this is the mean home-range radius of rats in m
			hr.dat$ucl[i]<-try(circular.r (p = 0.95, detectfn = 'EX', sigma = OUT$q97.5$hrradius[hr.dat$ID[i],hr.dat$session[i]], hazard = TRUE),silent = TRUE)		### this is the mean home-range radius of rats in m
			hr.dat$lcl[i]<-try(circular.r (p = 0.95, detectfn = 'EX', sigma = OUT$q2.5$hrradius[hr.dat$ID[i],hr.dat$session[i]], hazard = TRUE),silent = TRUE)		### this is the mean home-range radius of rats in m

	}

head(hr.dat)		
hr.dat$lcl<-as.numeric(as.character(hr.dat$lcl))
hr.dat$ucl<-as.numeric(as.character(hr.dat$ucl))
dim(hr.dat)
head(hr.dat)

hr.dat$sessname<-"June"
hr.dat$sessname<-ifelse(hr.dat$session==1,"early Aug",hr.dat$sessname)
hr.dat$sessname<-ifelse(hr.dat$session==2,"late Aug",hr.dat$sessname)
hr.dat$sessname<-ifelse(hr.dat$session==3,"Sept",hr.dat$sessname)
hr.dat$sessname<-factor(hr.dat$sessname, levels=c("early Aug","late Aug","Sept"))



######## INSPECT THE NA radii ######################
hr.dat[is.na(hr.dat$lcl),]





###### CALCULATE PROPORTION OF CRAZY MOVERS PER SESSION ########
hr.dat$crazy<-ind.move$crazy[match(hr.dat$ID, ind.move$ID)]

ddply(hr.dat, c('session','crazy'), summarize,mean=mean(hr_radius), lcl=mean(lcl, na.rm=T), ucl=mean(ucl, na.rm=T))
aggregate(crazy~session, hr.dat, mean)




######### ALTERNATIVE FIGURE 3: PLOT OF MEANS ###########

hr.means.all <- ddply(hr.dat, c('crazy','session','sessname'), summarize,mean=mean(hr_radius), lcl=mean(lcl, na.rm=T), ucl=mean(ucl, na.rm=T))
hr.means<-hr.means.all[hr.means.all$crazy==0,]
roam.means<- hr.means.all[hr.means.all$crazy==1,]
roam.means$prop<-aggregate(crazy~session, hr.dat, mean)[3:4,2]
hr.means$prop<-1-(aggregate(crazy~session, hr.dat, mean)[,2])

par(mfrow=c(2,1),mar=c(4,5,1,1), oma=c(1,1,0,0))
errbar(hr.means$session,hr.means$mean,hr.means$lcl,hr.means$ucl,type='p', las=1,ylab="Home range radius (m)", xlab="", main="", xlim=c(0.5,3.5), ylim=c(0,600), cex=1.5, cex.lab=1.5, axes=F, pch=16,frame=F)
axis(1, at=c(1:3), labels=hr.means$sessname, cex.axis=1.5, cex.lab=1.5)
axis(2, at=seq(0,600,200), las=1, cex.axis=1.5, cex.lab=1.5)

errbar(roam.means$session,roam.means$mean,roam.means$lcl,roam.means$ucl,type='p', las=1,ylab="Home range radius (m)", xlab="", main="", xlim=c(0.5,3.5), ylim=c(0,1500), cex=1.5, cex.lab=1.5, axes=F, pch=16,frame=F)
axis(1, at=c(1:3), labels=hr.means$sessname, cex.axis=1.5, cex.lab=1.5)
axis(2, at=seq(0,1600,200), las=1, cex.axis=1.5, cex.lab=1.5)


write.table(hr.means,"clipboard", sep="\t", row.names=F)
write.table(roam.means,"clipboard", sep="\t", row.names=F)



########## CALCULATE BAIT COVERAGE ######


### CALCULATE BAIT IN MIN HOME RANGE
platBait<-10000			## g/ha
beachBait<-40000			## g/ha
pellet<-1.8				## g for one bait pellet


beachBaitDensity<-beachBait/pellet

((min(hr.dat$lcl[hr.dat$crazy==0])^2*pi)/10000)*beachBaitDensity
((min(hr.dat$ucl[hr.dat$crazy==0])^2*pi)/10000)*beachBaitDensity






######### FIGURE 3 FOR PAPER ###########

pdf("S:\\ConSci\\DptShare\\SteffenOppel\\MANUSCRIPTS\\in_prep\\Henderson\\Rat_homerange\\FIGURE3_BB.pdf", width=12, height=8)
ggplot(hr.dat, aes(x=hr_radius)) +   facet_wrap(~sessname, ncol=2)+
    geom_histogram(binwidth = 30) +      # Thinner lines
    xlab("Home range radius (m)") +
    ylab("N individual rats") +
    #theme_bw()

  	geom_vline(data=hr.means, aes(xintercept = mean), colour="black") +
  	geom_vline(data=hr.means, aes(xintercept = lcl), colour="black", linetype = "dashed") +
  	geom_vline(data=hr.means, aes(xintercept = ucl), colour="black", linetype = "dashed") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank())
dev.off()









##########################################################################################################
################ PLOT SURVIVAL PROBABILITY BETWEEN SESSIONS #####################################
##########################################################################################################
surv.dat<-out[1:3,]			## annual survival
names(surv.dat)[c(3,7)]<-c("lcl","ucl")

surv.dat$session<-as.numeric(substr(surv.dat$parameter,14,14))			### change to 13 for constant survival
surv.dat$sessname<-"June"
surv.dat$sessname<-ifelse(surv.dat$session==1,"early Aug",surv.dat$sessname)
surv.dat$sessname<-ifelse(surv.dat$session==2,"late Aug",surv.dat$sessname)
surv.dat$sessname<-ifelse(surv.dat$session==3,"Sept",surv.dat$sessname)



## FOR INFERENCE, REMOVE THE LAST SURVIVAL ESTIMATE ###
#surv.dat<-surv.dat[1:5,]

pdf("S:\\ConSci\\DptShare\\SteffenOppel\\MANUSCRIPTS\\in_prep\\Henderson\\Rat_homerange\\Rat_surv_prob_MarkovSECR_moveMix.pdf", width=11, height=8)
par(mar=c(4,5,1,1), oma=c(1,1,0,0))
errbar(surv.dat$session,surv.dat$mean,surv.dat$lcl,surv.dat$ucl,type='p', las=1,ylab="Monthly survival probability", xlab="", main="", xlim=c(0.5,7.5), ylim=c(0,1), cex=1.5, cex.lab=1.5, axes=F, pch=16,col='black', frame=F)
axis(1, at=c(1:6), labels=surv.dat$sessname, cex.axis=1.5, cex.lab=1.5)
axis(2, at=seq(0,1,0.2), las=1, cex.axis=1.5, cex.lab=1.5)
dev.off()





##########################################################################################################
################ PLOT HOME RANGE SIZE BETWEEN SESSIONS #####################################
##########################################################################################################
### NEED TO DO: WEED OUT THE 'NON-RAT' ESTIMATES OF INDIVIDUALSxSESSION where the rat wasn't caught
hr.dat<-out[grep("hrradius",out$parameter),]
dim(hr.dat)
names(hr.dat)[c(3,7)]<-c("lcl","ucl")
hr.dat<-hr.dat[!hr.dat$ucl==0,]
hr.dat<-hr.dat[!hr.dat$lcl==0,]
dim(hr.dat)

hr.dat$session<-as.numeric(substr(hr.dat$parameter,nchar(hr.dat$parameter)-1,nchar(hr.dat$parameter)-1))
hr.dat$sessname<-"June"
hr.dat$sessname<-ifelse(hr.dat$session==1,"June",hr.dat$sessname)
hr.dat$sessname<-ifelse(hr.dat$session==2,"early July",hr.dat$sessname)
hr.dat$sessname<-ifelse(hr.dat$session==3,"late July",hr.dat$sessname)
hr.dat$sessname<-ifelse(hr.dat$session==4,"early Aug",hr.dat$sessname)
hr.dat$sessname<-ifelse(hr.dat$session==5,"late Aug",hr.dat$sessname)
hr.dat$sessname<-ifelse(hr.dat$session==6,"Sept",hr.dat$sessname)
hr.dat$sessname<-ifelse(hr.dat$session==7,"Oct",hr.dat$sessname)
hr.dat$sessname<-factor(hr.dat$sessname, levels=c("June","early July","late July","early Aug","late Aug","Sept","Oct"))

### CALCULATE HOME RANGE RADIUS
library(secr)
hr.dat$HR_lcl<-0
for (l in 1:dim(hr.dat)[1]){
hr.dat$HR_rad[l]<-circular.r (p = 0.95, detectfn = 'EX', sigma = hr.dat$mean[l], hazard = TRUE)		### this is the mean home-range radius of rats in m
hr.dat$HR_lcl[l]<-try(circular.r (p = 0.95, detectfn = 'EX', sigma = hr.dat$lcl[l], hazard = TRUE),silent = TRUE)		### this is the mean home-range radius of rats in m
hr.dat$HR_ucl[l]<-circular.r (p = 0.95, detectfn = 'EX', sigma = hr.dat$ucl[l], hazard = TRUE)		### this is the mean home-range radius of rats in m
}


head(hr.dat)


#pdf("S:\\ConSci\\DptShare\\SteffenOppel\\MANUSCRIPTS\\in_prep\\Henderson\\Rat_homerange\\Rat_HR_size_covariates.pdf", width=11, height=8)
par(mar=c(4,5,1,1), oma=c(1,1,0,0))
errbar(hr.dat$session,hr.dat$HR_rad,hr.dat$HR_lcl,hr.dat$HR_ucl,type='p', las=1,ylab="Home range radius (m)", xlab="", main="", xlim=c(0.5,7.5), ylim=c(0,1000), cex=1.5, cex.lab=1.5, axes=F, pch=16,col='darkred', frame=F)
axis(1, at=c(1:7), labels=hr.dat$sessname, cex.axis=1.5, cex.lab=1.5)
axis(2, at=seq(0,1000,200), las=1, cex.axis=1.5, cex.lab=1.5)
dev.off()


library(plyr)
hr.means <- ddply(hr.dat, c('session','sessname'), summarize,mean=mean(HR_rad), lcl=mean(HR_lcl, na.rm=T), ucl=mean(HR_ucl))




#pdf("S:\\ConSci\\DptShare\\SteffenOppel\\MANUSCRIPTS\\in_prep\\Henderson\\Rat_homerange\\Rat_HR_size_histograms.pdf", width=12, height=8)
ggplot(hr.dat, aes(x=HR_rad)) +   facet_wrap(~sessname, ncol=2)+
    geom_histogram(binwidth = 30) +      # Thinner lines
    xlab("Home range radius (m)") +
    ylab("N individual rats") +
    #theme_bw()

  	#geom_vline(data=hr.means, aes(xintercept = mean), colour="black") +
  	#geom_vline(data=hr.means, aes(xintercept = lcl), colour="black", linetype = "dashed") +
  	#geom_vline(data=hr.means, aes(xintercept = ucl), colour="black", linetype = "dashed") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank())
dev.off()



##########################################################################################################
################ PLOT HOME RANGE SHIFT BETWEEN SESSIONS #####################################
##########################################################################################################
### NEED TO DO: WEED OUT THE 'NON-RAT' ESTIMATES OF INDIVIDUALSxSESSION where the rat wasn't caught

move.dat<-out[grep("hrshift",out$parameter),]
dim(move.dat)
names(move.dat)[c(3,7)]<-c("lcl","ucl")
#move.dat<-move.dat[!move.dat$ucl==0,]
head(move.dat)
dim(move.dat)


move.dat$session<-as.numeric(substr(move.dat$parameter,nchar(move.dat$parameter)-1,nchar(move.dat$parameter)-1))
move.dat$sessname<-"June"
move.dat$sessname<-ifelse(move.dat$session==2,"June",move.dat$sessname)
move.dat$sessname<-ifelse(move.dat$session==3,"early July",move.dat$sessname)
move.dat$sessname<-ifelse(move.dat$session==4,"late July",move.dat$sessname)
move.dat$sessname<-ifelse(move.dat$session==5,"early Aug",move.dat$sessname)
move.dat$sessname<-ifelse(move.dat$session==6,"late Aug",move.dat$sessname)
move.dat$sessname<-ifelse(move.dat$session==7,"Sept",move.dat$sessname)
move.dat$sessname<-factor(move.dat$sessname, levels=c("June","early July","late July","early Aug","late Aug","Sept"))


#pdf("S:\\ConSci\\DptShare\\SteffenOppel\\MANUSCRIPTS\\in_prep\\Henderson\\Rat_homerange\\Rat_HR_shift_histograms.pdf", width=12, height=8)
ggplot(move.dat, aes(x=mean)) +   facet_wrap(~sessname, ncol=2)+
    geom_histogram(binwidth = 30) +      # Thinner lines
    xlab("Shift of virtual home range centre (m)") +
    ylab("N individual rats") +
    #theme_bw()
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank())
dev.off()













##########################################################################################################
################ CORRELATION BETWEEN SURVIVAL AND MOVEMENT PARAMETERS #####################################
##########################################################################################################

str(ratmodel$sims.list$sigma)
str(ratmodel$sims.list$monthly.surv)

correlations<-data.frame(session=c(1:6), sigma=0,surv=0)
for (sess in 1:6){
title<-cor.test(ratmodel$sims.list$sigma[,sess],ratmodel$sims.list$monthly.surv[,sess])
#plot(ratmodel$sims.list$sigma[,sess]~ratmodel$sims.list$monthly.surv[,sess], pch=16, cex=0.5, xlab="", ylab="")
correlations$sigma[correlations$session==sess]<-mean(ratmodel$sims.list$sigma[,sess])
correlations$surv[correlations$session==sess]<-mean(ratmodel$sims.list$monthly.surv[,sess])
}
}

plot(sigma~surv, correlations)
cor.test(correlations$sigma,correlations$surv, method='spearman')

plot(ratmodel$sims.list$sigma[,1:6]~ratmodel$sims.list$monthly.surv)
cor.test(ratmodel$sims.list$sigma[,1:6],ratmodel$sims.list$monthly.surv, method='spearman')