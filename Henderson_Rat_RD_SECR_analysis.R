#######################################################################################
##  RAT DENSITY ESTIMATION ON HENDERSON ISLAND  - RUNNING THE RD SECR MODEL 		###
#######################################################################################

## model based on Ergon, T., Gardner, B., 2014. Separating mortality and emigration: modelling space use, dispersal and survival with robust-design spatial capture–recapture data. Meth. Ecol. Evol. 5, 1327–1336.
## modified by steffen.oppel@rspb.org.uk on 30 May 2015
## updated on 11 Dez 2015, after incorporating Phase 2 data
## updated on 4 April 2016: introduced a secondary-occasion and location specific recapture probability depending on trap type
## altered model on 6 June 2016 to avoid error when PI ==0 by introducing temporary fudge (advise from Adam Butler)





###################################################################
### MODEL FORMULATION FOR RD SECR MODEL				    ###
### adopted from Ergon and Gardner				          ###
###################################################################

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Henderson\\Analyses\\Rats")
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Henderson\\Analyses\\Rats")

save.image("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\RD_SECR_input_data.RData")


# PARAMETERS:
# - kappa, sigma and lambda: Space use and recapture probability parameters
#                            (eq. 5 in paper)
# - beta: additive effects on log(lambda):
#           - beta[1]: effect of tod==2
#           - beta[2]: effect of gr==2
# - Phi[gr[i]]:   Survival for group gr[i] per time-unit 
# - phi[gr[i],k]: Survival probability for individual i belonging to group (sex) 
#                 gr[i] over the k'th interval given that the individual is
#                 alive at the beginning of the interval.
# - dmean[gr[i]]: Mean dispersal distance (exponential distribution) given
#                 dispersal for group gr[i]
# - Psi[gr[i]]:   Probability of not emigrating for group gr[i] per time-unit
#                 (only for zero-inflated model, commented out)
# - psi[gr[i],k]: Probability of dispersal (1-zero-inflation probability) for
#                 group gr[i] in interval k (only for zero-inflated model, 
#                 commented out)
# - d.offset[gr[i]]: Minimum dispersal distance given dispersal for group gr[i]
#                    (only for zero-inflated model, commented out)                 
# STATE VARIABLES:
# - z[i,k]: Alive state of individual i in primary session k (1=alive, 0=dead)
# - S[i,,k]: Centre of activity of individual i in primary session k (x- and
#            y-coordinate)
# - u[i,k]: Dispersal state for individual i in interval k (1=dispersal,
#           0=no dispersal)(only for zero-inflated model, commented out)



# Writing the model specification to a text file:
sink("Rat_RD_SECR.jags")
cat("

model{

## PRIORS AND CONSTRAINTS: ##

#  Space use and recapture probability parameters:
for(sex in 1:2){
  kappa[sex] ~ dunif(0,50)
  sigma[sex] ~ dunif(0.01,20)
}
for(sex in 1:2){
  for(TOD in 1:2){
    lambda[TOD, sex] <- exp(log.lambda0) * pow(beta[1],(TOD-1)) * pow(beta[2],(sex-1))
  }
}
PL ~ dunif(0.01,0.99)
log.lambda0 <- log(-log(1-PL)) 
beta[1] ~ dunif(0.1,10)
beta[2] ~ dunif(0.1,10)

# Survival parameters:
for(sex in 1:2){
  Phi[sex] ~ dunif(0,1)
  for(k in 1:(n.prim-1)){
    phi[sex,k] <- pow(Phi[sex], dt[k])
  }
}

# Dispersal parameters:
for(sex in 1:2){
  dmean[sex] ~ dunif(0,1000)
  dlambda[sex] <- 1/dmean[sex]
# For zero-inflated model:  
#  d.offset[sex] ~ dunif(5,100)
#  Psi0[sex] ~ dunif(0,1)
#  for(k in 1:(n.prim-1)){
#    psi[sex,k] <- 1 - pow(Psi0[sex], dt[k])
#  }
}

## MODEL: ##

# Loop over individuals that are only seen in the last primary session or the
# session they are censored
for(i in 1:N[1]){
  z[i,first[i]] ~ dbern(1)
  S[i,1,first[i]] ~ dunif(xlow[i], xupp[i]) # Prior for the first x coordinate
  S[i,2,first[i]] ~ dunif(ylow[i], yupp[i]) # Prior for the first y coordinate
  g[i,first[i],1] <- 0
  for(r in 1:R){ # trap
      D[i,r,first[i]] <- sqrt(pow(S[i,1,first[i]]-X[r,1],2) + pow(S[i,2,first[i]]-X[r,2],2))
      g[i,first[i],r+1] <- exp(-pow(D[i,r,first[i]]/sigma[gr[i]], kappa[gr[i]])) # Trap exposure
  }
  G[i,first[i]] <- sum(g[i,first[i],]) # Total trap exposure
  for(j in 1:J[i,first[i]]){
      P[i,j,first[i]] <- 1 - exp(-lambda[Htt2[i,j,first[i]],gr[i]]*G[i,first[i]]) ### Probability of being captured, CHANGED tod[first[i],j] to trap location and primary period tod[H[i,j,first[i]]-1,first[i]] but then created new array Htt2
    PI[i,first[i],j] <- step(H[i,j,first[i]]-2)*(g[i,first[i],H[i,j,first[i]]]/(G[i,first[i]]+ 0.000000001))*P[i,j,first[i]] + (1-step(H[i,j,first[i]]-2))*(1-P[i,j,first[i]])
    PIdum[i,first[i],j]<- PI[i,first[i],j]+equals(PI[i,first[i],j],0)*0.0001
    Ones[i,j,first[i]] ~ dbern(PIdum[i,first[i],j])
  }
}

# Loop over all other individuals
for(i in (N[1]+1):N[2]){
  z[i,first[i]] ~ dbern(1)
  S[i,1,first[i]] ~ dunif(xlow[i], xupp[i]) # Prior for the first x coordinate
  S[i,2,first[i]] ~ dunif(ylow[i], yupp[i]) # Prior for the first y coordinate
  # First primary session:
  g[i,first[i],1] <- 0
  for(r in 1:R){ # trap
      D[i,r,first[i]] <- sqrt(pow(S[i,1,first[i]]-X[r,1],2) + pow(S[i,2,first[i]]-X[r,2],2))
      g[i,first[i],r+1] <- exp(-pow(D[i,r,first[i]]/sigma[gr[i]], kappa[gr[i]])) # Trap exposure
  }
  G[i,first[i]] <- sum(g[i,first[i],]) # Total trap exposure
  for(j in 1:J[i,first[i]]){
      P[i,j,first[i]] <- 1 - exp(-lambda[Htt2[i,j,first[i]],gr[i]]*G[i,first[i]]) # Probability of being captured
    PI[i,first[i],j] <- step(H[i,j,first[i]]-2)*(g[i,first[i],H[i,j,first[i]]]/(G[i,first[i]]+ 0.000000001))*P[i,j,first[i]] + (1-step(H[i,j,first[i]]-2))*(1-P[i,j,first[i]])
	PIdum[i,first[i],j]<- PI[i,first[i],j]+equals(PI[i,first[i],j],0)*0.0001    
	Ones[i,j,first[i]] ~ dbern(PIdum[i,first[i],j])
  }
  ## Later primary sessions
  for(k in (first[i]+1):K[i]){ # primary session
    theta[i,k-1] ~ dunif(-3.141593,3.141593) # Prior for dispersal direction 
    z[i,k] ~ dbern(Palive[i,k-1])
    Palive[i,k-1] <- z[i,k-1]*phi[gr[i],k-1] # Pr(alive in primary session k) gr[i] = sex
    d[i,k-1] ~ dexp(dlambda[gr[i]])
    # For zero-inflated model, replace line above with:
    #u[i,k-1] ~ dbern(psi[gr[i],k-1])
    #dd[i,k-1] ~ dexp(dlambda[gr[i]])
    #d[i,k-1] <- u[i,k-1]*(d.offset[gr[i]] + dd[i,k-1])
    S[i,1,k] <- S[i,1,k-1] + d[i,k-1]*cos(theta[i,k-1])
    S[i,2,k] <- S[i,2,k-1] + d[i,k-1]*sin(theta[i,k-1])
    g[i,k,1] <- 0
    for(r in 1:R){ # trap
      D[i,r,k] <- sqrt(pow(S[i,1,k]-X[r,1],2) + pow(S[i,2,k]-X[r,2],2))  # Squared distance to trap
      g[i,k,r+1] <- exp(-pow(D[i,r,k]/sigma[gr[i]], kappa[gr[i]])) # Trap exposure
    }
    G[i,k] <- sum(g[i,k,]) # Total trap exposure
    for(j in 1:J[i,k]){
      P[i,j,k] <- (1 - exp(-lambda[Htt2[i,j,first[i]],gr[i]]*G[i,k]))*z[i,k] # Probability of being captured
    PI[i,k,j] <- step(H[i,j,k]-2)*(g[i,k,H[i,j,k]]/(G[i,k] + 0.000000001))*P[i,j,k] + (1-step(H[i,j,k]-2))*(1-P[i,j,k])
		PIdum[i,k,j]<- PI[i,k,j]+equals(PI[i,k,j],0)*0.0001  
    Ones[i,j,k] ~ dbern(PIdum[i,k,j])
    }
  }
}
}

",fill = TRUE)
sink()









##############################################
### LOADING WORKSPACE WITH PREPARED DATA   ###
##############################################
## data preparation in script "Henderson_RD_SECR_Ergon_Gardner2014.r"

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Henderson\\Analyses\\Rats")
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Henderson\\Analyses\\Rats")
load("RD_SECR_input_data.RData")

library(reshape)
library(slam)
library(doParallel)
library(snow)
library(dclone)
library(R2jags)
library(jagsUI)
library(Hmisc)





##############################################
### FUNCTION FOR GENERATING INITIAL VALUES ###
##############################################

Inits = function(H,X,K){
  n = dim(H)[1]
  n.prim = dim(H)[3]
  mean.x = apply(H, c(1,3), function(i) mean(X[i-1,1], na.rm=T))
  mean.y = apply(H, c(1,3), function(i) mean(X[i-1,2], na.rm=T))

  # For initial values of dispersal distance:
  for(i in (n.prim-1):1){
    mean.x[,i][is.na(mean.x[,i])] = mean.x[,i+1][is.na(mean.x[,i])]
    mean.y[,i][is.na(mean.y[,i])] = mean.y[,i+1][is.na(mean.y[,i])]
  }
  dx = mean.x[,2:n.prim] - mean.x[,1:(n.prim-1)]
  dy = mean.y[,2:n.prim] - mean.y[,1:(n.prim-1)]
  d.emp = sqrt(dx^2 + dy^2)

  ch = apply(H,c(1,3), function(i) any(i!=1))
  first.last = apply(ch, 1, function(i) range(which(i)))
  z = ch

  theta = atan2(dy,dx)
  theta[is.na(theta)] = runif(sum(is.na(theta)), -pi, pi)

  d = d.emp
  d[is.na(d)] = 0.001

  S = array(NA, c(n,2,n.prim)) # For initial values og FIRST location

  for(i in 1:n){
    S[i,1,first.last[1,i]] = mean.x[i,first.last[1,i]]
    S[i,2,first.last[1,i]] = mean.y[i,first.last[1,i]]
    z[i, first.last[1,i]:first.last[2,i]] = 1   # 1 when known to be alive, 0 otherwise
    if(first.last[1,i] != 1){
      theta[i,1:(first.last[1,i]-1)] = NA
      z[i,1:(first.last[1,i]-1)] = NA
      d[i,1:(first.last[1,i]-1)] = NA
    }
    if(K[i]!=n.prim){ # Adding NA's for censored individuals
      theta[i,K[i]:(n.prim-1)] = NA
      z[i, (K[i]+1):n.prim] = NA # after being removed
      d[i, K[i]:(n.prim-1)] = NA
    }
  }
  list(
    kappa = runif(2, 1.9, 2.1),
    sigma = runif(2, 4, 11), # Do not set too high if kappa is low
    PL = runif(1,0.7,0.9),
    beta = runif(2,0.5,2),
  	dmean = c(10,15) + runif(2,-3,3),
    Phi = runif(2,0.5,0.8),
    theta = theta,
    d = d,
    z = z,
    S = S
  )
}




##################################################################
### PREPARING THE DATA FILE TO PASS TO JAGS			 ###
##################################################################

# Priors for first centre of activity:
# Assuming that first centre of activity is always within +/- maxd
# from the mean capture location in both x and y direction.
trap.dist<-20
maxd = 2*trap.dist
mean.x = apply(H, c(1,3), function(i) mean(X[i-1,1], na.rm=T))
mean.y = apply(H, c(1,3), function(i) mean(X[i-1,2], na.rm=T))
first.mean.x = apply(mean.x,1, function(i) i[min(which(!is.na(i)))])
first.mean.y = apply(mean.y,1, function(i) i[min(which(!is.na(i)))])
xlow = first.mean.x - maxd
xupp = first.mean.x + maxd
ylow = first.mean.y - maxd
yupp = first.mean.y + maxd



JAGS.data = list(
  N = N,					### removed this c(sum((K-first)==0), length(first))
  K = K,
  R = nrow(X),
  J = J,
  #tod = tod,				### this used to be time of day in original model, replaced by Htt2 to allow occasion and location specific capt prob
  Htt2 = Htt2,
  first = first,
  X = as.matrix(X),
  H = H,
  n.prim = dim(H)[3],
  dt = dt,
  gr = gr,
  Ones = array(1, dim(H)),
  # Prior parameters for initial centre of activity:
  xlow = xlow,
  xupp = xupp,
  ylow = ylow,
  yupp = yupp
)



##################################################################
### FITTING THE MODEL IN JAGS ###
##################################################################


# Parameters to monitor
params = c("kappa", "sigma", "log.lambda0", "beta", "dmean", "phi", "Phi")

# MCMC settings for the model run 
ni <- 2500		## number of iterations
nt <- 1		## thinning rate
nb <- 1000		## 'burn-in' - the number of iterations that are discarded at the start of each chain
nc <- 8		## number of chains

# Generating initial values for all chains
inits = list(Inits(H,X=X,K=K),Inits(H,X=X,K=K),Inits(H,X=X,K=K),Inits(H,X=X,K=K),Inits(H,X=X,K=K),Inits(H,X=X,K=K))

# Fitting the model in JAGS (iter+burnin = 3500 on 6 chains took 8 hrs)
ratmodel <- jags(JAGS.data, inits, params, "S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Henderson\\Analyses\\Rats\\Rat_RD_SECR.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=T, n.cores=nc)

# Summary:
summary(ratmodel)

# Export Results
out<-as.data.frame(ratmodel$summary)
out$parameter<-row.names(ratmodel$summary)
write.table(out,"Henderson_rat_RD_SECR_output.csv", sep=",", row.names=F)

