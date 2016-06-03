###################################################################################################################
######### HENDERSON ISLAND RAIL BAIT CONSUMPTION TRIAL ############################################################
###################################################################################################################

## written by steffen.oppel@rspb.org on 21 June 2015
## modified 14 March 2016 to revise analysis based on Wildlife Research comments

## NOTE THAT PELLET CONSUMPTION IS SET TO 0 FOR SOME DAYS AND BIRDS WHEN MASSES INDICATE PELLET LOSS -> The 0 is correct because pellets were lost in cages due to rail rage etc. and NOT due to consumption


#### load packages
library(RODBC)
library(lme4)
library(MuMIn)
library(AICcmodavg)
library(piecewiseSEM)			## for R2 of mixed models

#### set working directory
setwd("A:\\RSPB\\UKOT\\Henderson\\Data")


#### load data
SP<-odbcConnectAccess2007('HendersonData_2015.accdb')
rails <- sqlQuery(SP, "SELECT * FROM Rail_pellet_consumption")  
mass <- sqlQuery(SP, "SELECT * FROM Rail_mass_change") 
odbcClose(SP)
rails$mass<-mass$Mass[match(rails$Bird_ID, mass$Bird_ID)]			### insert capture mass of birds into feeding table


############################################# ANALYSE BODY MASS CHANGE #####################################
head(mass)
mass<-mass[!is.na(mass$Mass.1),]
mass$diff<-mass$Mass.1-mass$Mass
mass$diffprop<-(mass$diff/mass$Mass)*100
aggregate(diff~sex, mass, FUN=mean)
aggregate(diff~sex, mass, FUN=sd)
aggregate(diffprop~sex, mass, FUN=mean)
aggregate(diffprop~sex, mass, FUN=sd)

t.test(mass$Mass.1[mass$sex=="male"],mass$Mass[mass$sex=="male"], paired=T)

t.test(mass$Mass.1[mass$sex=="female"],mass$Mass[mass$sex=="female"], paired=T)





############################################# ANALYSE PELLET CONSUMPTION #####################################


#### manipulate data [remove NA and data from the control trays]
head(rails)
dim(rails)


#### calculate consumption only for rows where it is not set to 0
rails$consumed[is.na(rails$consumed)]<-rails$Pellet_in[is.na(rails$consumed)]-rails$Pellet_out[is.na(rails$consumed)]

rails<-rails[!is.na(rails$consumed),]
rails$foodconsumption<-rails$food_in-rails$food_out
rails$foodconsumption[is.na(rails$foodconsumption)]<-mean(rails$foodconsumption, na.rm=T)		## assume average consumption if no data
dim(rails)
unique(rails$Bird_ID)
rails<-subset(rails,Bird_ID>0)
dim(rails)
rails$Bird_ID<-as.factor(as.character(rails$Bird_ID))
rails$Date<-as.factor(as.character(rails$Date))
str(rails)

#### proportional consumption
rails$prop_consumed<-(rails$Pellet_in-rails$Pellet_out)/rails$Pellet_in
rails$prop_consumed<-ifelse(rails$consumed==0,0,rails$prop_consumed)


#### fill in blanks for sex and pool immature and adult  [makes no difference whether male or female - outcome is identical for bait colour]
rails$sex<-as.character(rails$sex)
#rails$sex<-ifelse(is.na(rails$sex)==TRUE,"male",rails$sex)
rails$sex<-ifelse(is.na(rails$sex)==TRUE,"female",rails$sex)
#rails<-rails[!is.na(rails$sex),]		## exclude the birds without sex id
rails$sex<-as.factor(rails$sex)

rails$age<-as.character(rails$age)
rails$age<-ifelse(rails$age=="immature","adult",rails$age)
rails$age<-as.factor(rails$age)




##### CLEAN DATA AND REMOVE DAYS WHERE RAIN WASHED OUT BOWLS #####
rails$consumed[rails$consumed<0]<-0
rails$prop_consumed[rails$prop_consumed<0]<-0



#### inspect and visualise data
summary(rails[rails$bait_colour=="blue",])
summary(rails[rails$bait_colour=="green",])

simplesummary<-aggregate(consumed~bait_colour+sex,rails[rails$bait_state=="wet",],FUN=mean)
simplesummary$sd<-aggregate(consumed~bait_colour+sex,rails[rails$bait_state=="wet",],FUN=sd)[,3]
simplesummary
write.table(simplesummary, "clipboard", sep="\t")

#par(mfrow=c(2,1))
hist(rails$consumed[rails$bait_colour=="blue"], breaks=20, col="blue", xlim=c(-10,40), ylim=c(0,10), main="", xlab="amount of pellets consumed (g)")
par(new=T)
hist(rails$consumed[rails$bait_colour=="green"], breaks=20, col="green", xlim=c(-10,40), ylim=c(0,10), main="", xlab="")


#### plot summary
par(mar=c(3,5,0,0))
errbar(c(0.9,1.1,1.9,2.1), simplesummary$consumed, (simplesummary$consumed-0.5*simplesummary$sd), (simplesummary$consumed+0.5*simplesummary$sd), xlim=c(0,3), ylim=c(0,20),pch=16, cex=2, col="darkgreen", axes=F, xlab="",ylab="Pellets consumed", cex.lab=1.7)
axis(1, at=c(0:3), labels=c("","female","male",""), cex.axis=1.5)
axis(2, at=seq(0,20,5), labels=T,las=1, cex.axis=1.5)
points(c(0.9,1.9),simplesummary$consumed[c(1,3)],pch=16, cex=2, col="darkblue")
points(c(1.1,2.1),simplesummary$consumed[c(2,4)],pch=16, cex=2, col="darkgreen")

#### simple analysis
t.test(consumed~bait_colour, data=rails)

t.test(consumed~bait_colour, data=rails[rails$bait_state=="wet",])


####
cor.test(rails$consumed,rails$foodconsumption)
plot(rails$consumed~rails$foodconsumption)



#### INDIVIDUAL SUMMARY
indsummary<-aggregate(prop_consumed~bait_colour+Bird_ID+sex,rails[rails$bait_state=="wet",],FUN=mean)
indsummary$sd<-aggregate(prop_consumed~bait_colour+Bird_ID+sex,rails[rails$bait_state=="wet",],FUN=sd)[,4]
indsummary$min<-aggregate(prop_consumed~bait_colour+Bird_ID+sex,rails[rails$bait_state=="wet",],FUN=min)[,4]
indsummary$max<-aggregate(prop_consumed~bait_colour+Bird_ID+sex,rails[rails$bait_state=="wet",],FUN=max)[,4]
indsummary$mass<-aggregate(mass~bait_colour+Bird_ID+sex,rails[rails$bait_state=="wet",],FUN=mean)[,4]
indsummary[order(indsummary$bait_colour, indsummary$sex),]
write.table(indsummary[order(indsummary$bait_colour, indsummary$sex),], "clipboard", sep="\t", row.names=F)


######################## ANALYSE DATA ############################
## finalised on 11 Aug after trialling various options!
## multimodel inference
## updated 14 March 2016: tried including body mass


### EITHER WITH BODY MASS AS FIXED EFFECT
#m0<-lmer(consumed~mass+(1|Bird_ID), data=rails[rails$bait_state=="wet",], na.action=na.fail)		### consider age, sex, bait state
#m1<-lmer(consumed~mass+bait_colour+(1|Bird_ID), data=rails[rails$bait_state=="wet",], na.action=na.fail)
#m2<-lmer(consumed~mass+sex+(1|Bird_ID), data=rails[rails$bait_state=="wet",], na.action=na.fail)
#m3<-lmer(consumed~mass+bait_colour+sex+(1|Bird_ID), data=rails[rails$bait_state=="wet",], na.action=na.fail)
#m4<-lmer(consumed~mass+bait_colour:sex+(1|Bird_ID), data=rails[rails$bait_state=="wet",], na.action=na.fail)


### BECAUSE BODY MASS IS ACCOUNTED FOR IN IND RANDOM EFFECT WE INCORPORATED TRIAL_DAY
m0<-lmer(consumed~Trial_Day+(1|Bird_ID), data=rails[rails$bait_state=="wet",], na.action=na.fail)		### consider age, sex, bait state
m1<-lmer(consumed~Trial_Day+bait_colour+(1|Bird_ID), data=rails[rails$bait_state=="wet",], na.action=na.fail)
m2<-lmer(consumed~Trial_Day+sex+(1|Bird_ID), data=rails[rails$bait_state=="wet",], na.action=na.fail)
m3<-lmer(consumed~Trial_Day+bait_colour+sex+(1|Bird_ID), data=rails[rails$bait_state=="wet",], na.action=na.fail)
m4<-lmer(consumed~bait_colour:sex+Trial_Day+(1|Bird_ID), data=rails[rails$bait_state=="wet",], na.action=na.fail)
#anova(m0,m3)			### abandoned likelihood ratio test

modtab<-aictab(list(m0,m1, m2,m3,m4), modnames = c("null", "colour", "sex", "sex+colour", "sex x colour"))
out<-summary(m4)
sem.model.fits(m4)
write.table(modtab, "clipboard", sep="\t")

write.table(rails, "clipboard", sep="\t", row.names=F)


##### SUMMARISE DIFFERENCES IN CONSUMPTION BY PREDICTING TO NEW DATA #########
newdat<-data.frame(sex=rep(c("female","male"),2), bait_colour=rep(c("blue","green"), each=2), Trial_Day=rep(13,4), Bird_ID=rep(rails$Bird_ID[1],4))
newdat$consumed<-predict(m4, newdat, type = "response")

mm <- model.matrix(terms(m4),newdat)
pvar1 <- diag(mm %*% tcrossprod(vcov(m4),mm))
tvar1 <- pvar1+VarCorr(m4)[[1]][1]  			## must be adapted for more complex models
newdat$se<-sqrt(tvar1)
newdat$consumed_lcl<-newdat$consumed-2*sqrt(tvar1)
newdat$consumed_ucl<- newdat$consumed+2*sqrt(tvar1)
newdat
write.table(newdat, "clipboard", sep="\t")




##### SUMMARISE DIFFERENCES IN CONSUMPTION #########
consump<-out$coefficients

1-(consump[1,1]/consump[2,1])

1-(consump[3,1]/consump[4,1])



### model dredging - some support for bait colour
m1<-lmer(consumed~foodconsumption+mass+sex+bait_state+bait_colour+AvgOftemperature+AvgOfhumidity+(1|Bird_ID)+(1|Date), data=rails, na.action=na.fail)
dredge(m1)
summary(m1)




