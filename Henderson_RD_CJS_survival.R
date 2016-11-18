#######################################################################################
##  RAT SURVIVAL ESTIMATION ON HENDERSON ISLAND  - RUNNING THE SIMPLE CJS MODEL FOR EACH SESSION		###
#######################################################################################
## ADAPTED FROM RD SECR MODEL THAT WOULD NOT CONVERGE

## UPDATE 16 Nov 2016: WROTE SIMPLE CJS MODEL FOR EACH SESSION





###################################################################
### MODEL FORMULATION FOR SINGLE SESSION CJS SECR MODEL	###
###################################################################


# Writing the model specification to a text file:
sink("CJS_singlesess.jags")
cat("


model {


#### PRIORS ####

#  Space use and recapture probability parameters:
#beta.size ~ dunif(-10,10)## adjust sigma according to size
#beta.size.sq ~ dunif(-10,10)## adjust sigma according to size, quadratic effect

for (s in 1:2){			# movement varies by sex
sigma[s] ~ dunif(0, 200)     	
}

for(t in 1:T){			# capture probability varies by day
#for (u in 1:2){## capture probability differs by traptype
lam0[t]~dunif(0,1)
#}
}

for(s in 1:2){			# survival probability varies by sex
phi[s]~dunif(0, 1)
}


#### MODEL ####       
   for (i in 1:M){    # loop through the population of M trapped rats
       z[i,first[i]] <- 1       # CJS is conditional on capture, so alive state in first capture period must be 1   
       SX[i,first[i]]~dunif(xl[i], xu[i])    # priors for the activity centers for each individual at first capture
       SY[i,first[i]]~dunif(yl[i], yu[i])    # xl is the lower x coordinate, xu is the upper x value

## OBSERVATION MODEL FOR FIRST SESSION ##
     for(j in 1:ntraps) {     #loop through the J trap locations
           D2[i,j,first[i]] <- sqrt(pow(SX[i,first[i]]-grid[j,1], 2) + pow(SY[i,first[i]]-grid[j,2],2))## calculate distance from home range center to each trap location
           lam[i,j,first[i]] <- lam0[first[i]]*exp(-D2[i,j,first[i]]/(sigma[gr[i]]))## exponential detection function based on distance
           tmp[i,j,first[i]] <- lam[i,j,first[i]]*ON[j,first[i]]				## correct by activity index
           y[i,j,first[i]] ~ dbern(tmp[i,j,first[i]])						## bernoulli draw of capture 
}# close the trap loop for first capture period

## OBSERVATION MODEL FOR SUBSEQUENT SESSIONS ##       
for(t in (first[i]+1):T){   #loop through the primary periods after the first capture session
       SX[i,t]~dunif(xl[i], xu[i])    # priors for the activity centers for each individual at first capture
       SY[i,t]~dunif(yl[i], yu[i])    # xl is the lower x coordinate, xu is the upper x value

     for(j in 1:ntraps) {     #loop through the J trap locations
           D2[i,j,t] <- sqrt(pow(SX[i,t]-grid[j,1], 2) + pow(SY[i,t]-grid[j,2],2))## calculate distance from home range center to each trap location
           lam[i,j,t] <- lam0[t]*exp(-D2[i,j,t]/(sigma[gr[i]]))				## exponential detection function based on distance
           tmp[i,j,t] <- lam[i,j,t]*ON[j,t]*z[i,t]						## correct by activity index
           y[i,j,t] ~ dbern(tmp[i,j,t])								## bernoulli draw of capture probability 
				}# close the trap loop for subsequent capture periods

### SURVIVAL MODEL FOR TIME PERIOD T ##
	phiUP[i,t]<-z[i,t-1]*phi[gr[i]]			## sex-specific survival probability
	z[i,t]~dbern(phiUP[i,t])

				} ## close the primary session loop

}## close the individual loop


### DERIVED PARAMETERS FOR SURVIVAL AND CAPTURE PROBABILITY ###

for(s in 1:2){
month.surv[s]<- pow(phi[s],30)
}


} ## close the model loop



",fill = TRUE)
sink()










##############################################
### LOADING WORKSPACE WITH PREPARED DATA   ###
##############################################
## data preparation in script "Henderson_CJS_survival_data_prep_singlesession.R"


library(jagsUI)
library(secr)			## needed for home-range radius conversion from sigma
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Bayes_SECR")
#setwd("C:\\STEFFEN\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Bayes_SECR")
#load("RD_CJS_SECR_input_data.RData")

ALL_OUT<-data.frame()
for (sess in 1:7){
#sess=3

load(sprintf("CJS_input_sess%s.RData",sess))


#################################################################
### FUNCTION FOR GENERATING INITIAL VALUES FOR SIMPLE MODEL  ###
#################################################################


### Set up a data input for JAGS
CJS_data<-list(y=y, M=M, grid=grid, T=T,ntraps=ntrap, xl=xlow, yl=ylow, xu=xupp, yu=yupp, first=first, gr=gr, ON=ON)
#CJS_data<-list(y=y, M=max(rats$ID), K=K, grid=grid, T=T,ntraps=ntraps, xl=xl, yl=yl, xu=xu, yu=yu, first=first, gr=gr, ON=ON)


### Set the initial values
z<-matrix(NA,nrow=CJS_data$M,ncol=T)
for (i in 1:CJS_data$M){
z[i,first[i]]<- NA
if(first[i]<CJS_data$T){
	for(t in (first[i]+1):T){
 z[i,t]<- 1
}}}

geninits =  function() {list(z=z,phi=runif(2,0,1), lam0=runif(T,0,1))}

# Parameters to monitor
params = c('month.surv','lam0','sigma')


# MCMC settings for the model run  (1.5 hrs for 2000 iterations)
ni <- 4000		## number of iterations
nt <- 2		## thinning rate
nb <- 1000		## 'burn-in' - the number of iterations that are discarded at the start of each chain
nc <- 4		## number of chains

# Generating initial values for all chains
inits = list(geninits(),geninits(),geninits(),geninits())

# Fitting the model in JAGS
ratmodel <- jags(CJS_data, inits, params, "A:\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Bayes_SECR\\CJS_singlesess.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=T, n.cores=nc)



# Export Results
out<-as.data.frame(ratmodel$summary)
out$parameter<-row.names(ratmodel$summary)
out$session<-sess
ALL_OUT<-rbind(ALL_OUT,out)
write.table(out,sprintf("Henderson_CJS_surv_est_sess%s.csv",sess), sep=",", row.names=F)
#save.image("Henderson_RD_CJS_survival_estimates_sigma_ranef.RData")



} ### end loop over all 7 sessions



########################################################################################################################################
############### ARRANGING SUMMARY TABLE FOR HR SHIFT #################################################################
########################################################################################################################################
head(out)
dim(out)


surv.dat<-ALL_OUT[ALL_OUT$parameter %in% c("month.surv[1]","month.surv[2]"),]			## annual survival
hr.dat<-ALL_OUT[ALL_OUT$parameter %in% c("sigma[1]","sigma[2]"),]	





##########################################################################################################
################ PLOT SURVIVAL PROBABILITY BETWEEN SESSIONS #####################################
##########################################################################################################

surv.dat$sessname<-"June"
surv.dat$sessname<-ifelse(surv.dat$session==1,"June",surv.dat$sessname)
surv.dat$sessname<-ifelse(surv.dat$session==2,"early July",surv.dat$sessname)
surv.dat$sessname<-ifelse(surv.dat$session==3,"late July",surv.dat$sessname)
surv.dat$sessname<-ifelse(surv.dat$session==4,"early Aug",surv.dat$sessname)
surv.dat$sessname<-ifelse(surv.dat$session==5,"late Aug",surv.dat$sessname)
surv.dat$sessname<-ifelse(surv.dat$session==6,"Sept",surv.dat$sessname)
surv.dat$sex<-c('m','f')
plotF<-surv.dat[surv.dat$sex=='f',c(13,3,5,7,14)]
names(plotF)[c(2:4)]<-c("lcl","mean","ucl")
plotF$session<-c(1:7)
plotF$x<-as.numeric(plotF$session)-0.2
plotM<-surv.dat[surv.dat$sex=='m',c(13,3,5,7,14)]
names(plotM)[c(2:4)]<-c("lcl","mean","ucl")
plotM$session<-c(1:7)
plotM$x<-as.numeric(plotM$session)+0.2


par(mar=c(4,5,1,1), oma=c(1,1,0,0))
errbar(plotF$x,plotF$mean,plotF$lcl,plotF$ucl,type='p', las=1,ylab="Monthly survival probability", xlab="", main="", xlim=c(0.5,7.5), ylim=c(0,1), cex=1.5, cex.lab=1.5, axes=F, pch=16,col='darkred', frame=F)
par(new=T)
errbar(plotM$x,plotM$mean,plotM$lcl,plotM$ucl,type='p', axes=F, las=1,ylab="", xlab="", main="", xlim=c(0.5,7.5),ylim=c(0,1), cex=1.5, pch=16,col='navyblue', frame=F)
axis(1, at=c(1:7), labels=plotF$sessname, cex.axis=1.5, cex.lab=1.5)
axis(2, at=seq(0,1,0.2), las=1, cex.axis=1.5, cex.lab=1.5)
legend('topright', pch=16, col=c('darkred','navyblue'),legend=c("female","male"), cex=1.4, bty='n')






##########################################################################################################
################ PLOT HOME RANGE SIZE BETWEEN SESSIONS #####################################
##########################################################################################################


hr.dat$Sex<-rep(c('male','female'),7)
names(hr.dat)[c(3,7)]<-c("lcl","ucl")

hr.dat$session<-as.numeric(substr(hr.dat$parameter,9,9))
hr.dat$sessname<-"June"
hr.dat$sessname<-ifelse(hr.dat$session==1,"June",hr.dat$sessname)
hr.dat$sessname<-ifelse(hr.dat$session==2,"early July",hr.dat$sessname)
hr.dat$sessname<-ifelse(hr.dat$session==3,"late July",hr.dat$sessname)
hr.dat$sessname<-ifelse(hr.dat$session==4,"early Aug",hr.dat$sessname)
hr.dat$sessname<-ifelse(hr.dat$session==5,"late Aug",hr.dat$sessname)
hr.dat$sessname<-ifelse(hr.dat$session==6,"Sept",hr.dat$sessname)
hr.dat$sessname<-ifelse(hr.dat$session==7,"Oct",hr.dat$sessname)

### CALCULATE HOME RANGE RADIUS
for (l in 1:dim(hr.dat)[1]){
hr.dat$HR_rad[l]<-circular.r (p = 0.95, detectfn = 'EX', sigma = hr.dat$mean[l], hazard = TRUE)		### this is the mean home-range radius of rats in m
hr.dat$HR_lcl[l]<-circular.r (p = 0.95, detectfn = 'EX', sigma = hr.dat$lcl[l], hazard = TRUE)		### this is the mean home-range radius of rats in m
hr.dat$HR_ucl[l]<-circular.r (p = 0.95, detectfn = 'EX', sigma = hr.dat$ucl[l], hazard = TRUE)		### this is the mean home-range radius of rats in m
}


plotF<-subset(hr.dat, Sex=="female")
plotF$x<-as.numeric(plotF$session)-0.2
plotM<-subset(hr.dat, Sex=="male")
plotM$x<-as.numeric(plotM$session)+0.2


par(mar=c(4,5,1,1), oma=c(1,1,0,0))
errbar(plotF$x,plotF$HR_rad,plotF$HR_lcl,plotF$HR_ucl,type='p', las=1,ylab="Home range radius (m)", xlab="", main="", xlim=c(0.5,6.5), ylim=c(0,200), cex=1.5, cex.lab=1.5, axes=F, pch=16,col='darkred', frame=F)
par(new=T)
errbar(plotM$x,plotM$HR_rad,plotM$HR_lcl,plotM$HR_ucl,type='p', axes=F, las=1,ylab="", xlab="", main="", xlim=c(0.5,6.5),ylim=c(0,200), cex=1.5, pch=16,col='navyblue', frame=F)
axis(1, at=c(1:7), labels=plotF$sessname, cex.axis=1.5, cex.lab=1.5)
axis(2, at=seq(0,200,40), labels=T, las=1, cex.axis=1.5, cex.lab=1.5)
legend('topleft', pch=16, col=c('darkred','navyblue'),legend=c("female","male"), cex=1.4, bty='n')







