###################################################################
#### MURPHYS PETREL INCUBATION TRIP ANALYSIS			 ####
###################################################################
## written 7 September 2015 by steffen.oppel@rspb.org.uk
## modified to insert simple map on 23 Feb 2016
## updated 12 March 2016 to include analysis for seabird conference
## updated 12 May to make poster graph and table
## updated 12 May to include Tommy's proportion wet data
## 27 May 2016: added plot of departure directions

# Load necessary library
library(maptools)
library(rgdal)
require(maps)
require(mapdata)
require(maptools)
library(ggmap)
library(move)
library(RODBC)
require(geosphere)
library(rworldmap)
data(countriesLow)

source("C:\\STEFFEN\\RSPB\\Marine\\IBA\\Analysis\\tripSplit.r")
source("C:\\STEFFEN\\RSPB\\Marine\\IBA\\Analysis\\TripSummary.r")			## does not work for reasons I don't understand


source("A:\\RSPB\\Marine\\IBA\\Analysis\\tripSplit.r")
source("A:\\RSPB\\Marine\\IBA\\Analysis\\TripSummary.r")			## does not work for reasons I don't understand




#### read in deployment data from database
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Henderson\\Data")
setwd("A:\\RSPB\\UKOT\\Henderson\\Data")
SP<-odbcConnectAccess2007('HendersonData_2015.accdb')
birds <- sqlQuery(SP, "SELECT * FROM Logger_deployments")
odbcClose(SP)
birds<-birds[birds$species=="MUPE",]

birds$deployed<-as.POSIXct(paste(birds$MinOfDate, format(birds$MinOfRelease_time,format="%H:%M:%S")), format = "%Y-%m-%d %H:%M:%S")
birds$retrieved<-as.POSIXct(paste(birds$MaxOfDate, format(birds$MaxOfCapture_time,format="%H:%M:%S")), format = "%Y-%m-%d %H:%M:%S")
names(birds)[c(3,11:14)]<-c("Nest","Nest_Lat","Nest_Long","deployment_mass","retrieval_mass")
birds<-birds[,c(1:4,11:16)]


### INPUT LOCATION OF OUTPUT FILES HERE ####
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Henderson\\Data\\Seabirds\\PathTrack_GPS_data")
setwd("A:\\RSPB\\UKOT\\Henderson\\Data\\Seabirds\\PathTrack_GPS_data")

### THIS SHOULD RUN WITHOUT ANY CHANGES ####
targetfiles<-list.files(pattern = "\\.pos$")

output<-data.frame()

for (f in 3:length(targetfiles)){			### start at 2 to exclude HEPE

bird_track<-read.table(targetfiles[f], sep=",", header=F, skip=5)

names(bird_track)<-c("Day", "Month", "Year","Hour","Minute","Second","elapsed","SVs","Latitude","Longitude","Altitude","Accuracy_ind","Clock_offset","Voltage")
bird_track$Year<-bird_track$Year+2000


flength<-nchar(targetfiles[f])
bird_track<-bird_track[!is.na(bird_track$Latitude),]
bird_track$Species<-substr(targetfiles[f],1,4)
bird_track$Nest<-ifelse(is.na(as.numeric(substr(targetfiles[f],6,6)))==T,substr(targetfiles[f],1,5),substr(targetfiles[f],1,6))
bird_track$GPS_ID<-as.numeric(substr(targetfiles[f],flength-8,flength-4))
bird_track$date<-as.POSIXct(paste(bird_track$Year, bird_track$Month,bird_track$Day, sep="-"), format = "%Y-%m-%d")
bird_track$time<-as.POSIXct(paste(bird_track$Hour, bird_track$Minute,bird_track$Second, sep=":"), format = "%H:%M:%S")
bird_track$Date<-as.Date(bird_track$date, format = "%Y-%m-%d")
bird_track$Time<-format(bird_track$time,format="%H:%M:%S")
bird_track$DateTime<-as.POSIXct(paste(bird_track$Date, bird_track$Time), format = "%Y-%m-%d %H:%M:%S")
bird_track$TrackTime <- as.double(bird_track$DateTime)
bird_track<-bird_track[order(bird_track$DateTime),]
bird_track$Sequence<-seq(1,dim(bird_track)[1],1)
bird_track$ID<-as.integer(ifelse(is.na(as.numeric(substr(targetfiles[f],6,6)))==T,substr(targetfiles[f],5,5),substr(targetfiles[f],5,6)))

head(bird_track)

output<-rbind(output, bird_track)
}
head(output)
dim(output)


#### REMOVE NON-LOCATIONS ####
output<-output[!output$Latitude==0,]
dim(output)




#### CHECK FOR DUPLICATE TIMESTAMPS ###
#output$count<-1
#check<-aggregate(count~Nest+DateTime, output, FUN=sum)
#check[check$count>1,]
head(output)


#### REDUCE FOR FURTHER PROCESSING 
tracks<-output[,c(25,17,20,21,22,9:10,23,24)]
#names(tracks)[1]<-"ID"
head(tracks)
#write.table(tracks, "HEPE_GPS_positions_Henderson2015raw.csv", sep=",", row.names=F)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE A MOVE OBJECT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
animals<-unique(output$ID)
output<-output[order(output$ID, output$DateTime),]
#output$type<-outsum$type[match(output$ID, as.integer(substr(outsum$Nest,5,nchar(as.character(outsum$Nest)))))]
MOBJ_all<-move(x=output$Longitude, y=output$Latitude, time=output$DateTime, proj=CRS("+proj=longlat +ellps=WGS84"), animal=output$ID)
MOBJ_all_df <- as(MOBJ_all, "data.frame")
names(MOBJ_all_df)[10] <- "ID"
#MOBJ_all<-spTransform(x=MOBJ_all, CRSobj="+proj=aeqd", center=TRUE)
proj4string(MOBJ_all)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOT A SIMPLE MAP WITHOUT GOOGLE EARTH
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

jpeg("Fig28_MUPE_track_map.jpg", width=600, height=600, quality = 100)
par(mar=c(4.2,4,0.2,0.2), oma=c(0,0,0,0))
plot(MOBJ_all, type='l', xlab="Longitude", ylab="Latitude", xlim=c(-145,-75), ylim=c(-50,-10))
plot(countriesLow, col='darkgrey', add=T)
points(birds[,6],birds[,5], pch=16, col='darkred', cex=2)
map.scale(x=-143, y=-55, ratio=FALSE, relwidth=0.2)
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPLIT INTO FORAGING TRIPS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Trips <- NULL
loc<-aggregate(Nest_Lat~species, data=birds, FUN=mean)		## Colony location is mean of all nest locations
loc$Longitude<-aggregate(Nest_Long~species, data=birds, FUN=mean)[,2]
names(loc)[2]<-"Latitude"
for(i in 1:length(unique(tracks$ID)))
  {
  Temp <- subset(tracks, ID == unique(tracks$ID)[i])
  Trip <- tripSplit(Track=Temp, Colony=loc[loc$species=='MUPE',2:3], InnerBuff=3, ReturnBuff=25, Duration=2.5, plotit=T)
  if(i == 1) {Trips <- Trip} else
  Trips <- spRbind(Trips,Trip)

  }

#str(Trip)
dim(Trips@data[Trips@data$trip_id>0,])
head(Trips)






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SUMMARISE FORAGING DISTANCES FOR EACH TRIP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nests<-birds[,c(3,5,6)]
nests$ID<-as.integer(ifelse(nchar(as.character(nests$Nest))==5,substr(nests$Nest,5,5),substr(nests$Nest,5,6)))
names(nests)[2:3]<-c("Latitude","Longitude")


### SUMMARISE MAX DIST FROM COLONY AND TRIP TRAVELLING TIME FOR EACH TRIP

trip_distances<-data.frame(trip=as.numeric(unique(Trips@data$trip_id[Trips@data$trip_id>0])), max_dist=0, duration=0, total_dist=0, type=NA, departure=min(birds$deployed), return=max(birds$retrieved))
trip_distances<-trip_distances[trip_distances$trip>0,]	
trip_distances$ID<-Trips@data$ID[match(trip_distances$trip, Trips@data$trip_id)]

for (i in as.numeric(unique(trip_distances$trip))){
x<-Trips@data[Trips@data$trip_id==i,]
maxdist<-c(x$Longitude[x$ColDist==max(x$ColDist)], x$Latitude[x$ColDist==max(x$ColDist)])	### this used to be column 5:4, but unstable by numeric indexing
trip_distances[trip_distances$trip==i,2]<-max(Trips@data$ColDist[Trips@data$trip_id==i,])/1000
trip_distances[trip_distances$trip==i,3]<-(max(Trips@data$TrackTime[Trips@data$trip_id==i])-min(Trips@data$TrackTime[Trips@data$trip_id==i]))/3600


## Calculate distances from one point to the next and total trip distance
x$Dist[1]<-x$ColDist[1]/1000				### distance to first point is assumed a straight line from the nest/colony
for (p in 2:dim(x)[1]){
p1<-c(x$Longitude[p-1],x$Latitude[p-1])
p2<-c(x$Longitude[p],x$Latitude[p])
#x$Dist[p]<-pointDistance(p1,p2, lonlat=T, allpairs=FALSE)/1000			### no longer works in geosphere
x$Dist[p]<-distMeeus(p1,p2)/1000						### great circle distance according to Meeus, converted to km

}
trip_distances[trip_distances$trip==i,4]<-sum(x$Dist)+(x$ColDist[p]/1000)				## total trip distance is the sum of all steps plus the dist from the nest of the last location
trip_distances[trip_distances$trip==i,5]<-ifelse(max(x$Longitude, na.rm=T)>-120,'long','short')	## divide into long and short trips depending on whether they go far enough east

trip_distances[trip_distances$trip==i,6]<-min(x$DateTime)	## departure time of trip
trip_distances[trip_distances$trip==i,7]<-max(x$DateTime)	## return time of trip
trip_distances$n_locs[trip_distances$trip==i]<-dim(x)[1]		## number of locations per trip

origin<- nests[match(unique(x$ID), nests$ID),]
trip_distances$bearing[trip_distances$trip==i]<-bearing(origin[,3:2],maxdist)			## great circle route bearing of trip
}
trip_distances$Nest<-paste("MUPE",trip_distances$ID, sep="")

outsum<-merge(birds,trip_distances[,c(11,1,6,7,9,5,10,3,2,4)], by="Nest", all.x=T)


### AVERAGE TRAVEL SPEED FOR MUPES ###
outsum$travel_speed<-outsum$total_dist/outsum$duration



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE A FANCY MAP ON GOOGLE EARTH
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
output$type<-outsum$type[match(output$ID, as.integer(substr(outsum$Nest,5,nchar(as.character(outsum$Nest)))))]
output$colour<-ifelse(output$type=="long", "tomato3", "royalblue3")
map.scale.black<-map.scale
fix(map.scale.black)			
## change lines(linexy) to lines(linexy, lwd=6), 
## change text(x + ats/perdeg, y + dy - 0.5 * cxy[2], to text(x + ats/perdeg, y + 2+ dy - 0.5 * cxy[2],
## change linexy[3 * i + 2, ] <- c(x + i * dx, y + dy) to linexy[3 * i + 2, ] <- c(x + i * dx, y - 3*dy) 
## change linexy[2, ] <- c(x, y + dy) to linexy[2, ] <- c(x, y - 3*dy) 


jpeg("MUPE_poster_map.jpg", width=6000, height=4000, quality = 100)
par(mar=c(10,10,0.2,0.2), oma=c(0,0,0,0))
plot(Latitude~Longitude,data=output, type='l', xlab="", ylab="", xlim=c(-145,-60), ylim=c(-50,0), col=colour, lwd=10, axes=F)
plot(countriesLow, col='darkseagreen', add=T)
lines(Latitude~Longitude,data=output[output$type=="long",], type='l', col=colour, lwd=10)
points(birds[,6],birds[,5], pch=16, col='black', cex=20)
map.scale.black(x=-143, y=-49.5, ratio=FALSE, relwidth=0.4, cex=15)
dev.off()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE MASS GAIN DURING TRIP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## assume daily mass loss during incubation of 5.25g
str(outsum)
outsum$predep<-as.numeric(difftime(outsum$departure,outsum$deployed, units="days"))
outsum$postarr<-as.numeric(difftime(outsum$retrieved,outsum$return, units="days"))

outsum$dept_mass<-outsum$deployment_mass-(outsum$predep*5.25)
outsum$return_mass<-outsum$retrieval_mass+(outsum$postarr*5.25)
outsum$mass_gain<-outsum$return_mass-outsum$dept_mass
outsum$prop_mass_gain<-(outsum$mass_gain/outsum$dept_mass)*100

outsum$col<-ifelse(outsum$type=="long", "tomato3", "royalblue3")






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOT DEPARTURE DIRECTIONS FOR POSTER 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(circular)

#### FOR EACH TRIP, GET DIRECTION TO LOC AFTER 24 hrs ###

trip_directions<-data.frame(trip=as.numeric(unique(Trips@data$trip_id[Trips@data$trip_id>0])), dept_direction12=0, dept_direction24=0, dept_direction36=0)
trip_directions$Nest<-paste("MUPE",Trips@data$ID[match(trip_directions$trip, Trips@data$trip_id)], sep="")

for (i in as.numeric(unique(trip_distances$trip))){
x<-Trips@data[Trips@data$trip_id==i,]
origin<-birds[birds$Nest==trip_directions$Nest[match(i,trip_directions$trip)],6:5]			## extract long and lat of nest
dept_tim<-min(x$DateTime)			# departure time
dept12<-dept_tim+12*3600			# 12 hrs after departure
dept24<-dept_tim+24*3600			# 24 hrs after departure
dept36<-dept_tim+36*3600			# 36 hrs after departure
x<-x[x$DateTime>dept12,]
trip_directions[trip_distances$trip==i,2]<-ifelse(bearing(origin,x[1,7:6])<0,360+bearing(origin,x[1,7:6]),bearing(origin,x[1,7:6]))		## calculate direction and scale to 0-360
x<-x[x$DateTime>dept24,]
trip_directions[trip_distances$trip==i,3]<-ifelse(bearing(origin,x[1,7:6])<0,360+bearing(origin,x[1,7:6]),bearing(origin,x[1,7:6]))
x<-x[x$DateTime>dept36,]
trip_directions[trip_distances$trip==i,4]<-ifelse(bearing(origin,x[1,7:6])<0,360+bearing(origin,x[1,7:6]),bearing(origin,x[1,7:6]))
}	# end loop for individual trip


### MERGE with output summary
outsum<-merge(outsum, trip_directions, by=c("Nest","trip"), al.x=T)



### ANALYSE TRIP DIRECTIONS BETWEEN SHORT AND LONG TRIPS
directions<-circular(outsum$dept_direction12,units = "degrees", template =  "geographics", zero = 0, rotation = "clock")
aov.circular(directions, outsum$type)
directions<-circular(outsum$dept_direction24,units = "degrees", template =  "geographics", zero = 0, rotation = "clock")
aov.circular(directions, outsum$type)
directions<-circular(outsum$dept_direction36,units = "degrees", template =  "geographics", zero = 0, rotation = "clock")
aov.circular(directions, outsum$type)


### PLOT DEPARTURE DIRECTIONS
library(plotrix)
pdf("Petrel_departure_directions.pdf", width=8, height=8)
par(mar=c(10,10,10,10))
polar.plot(c(outsum$duration,250),c(outsum$dept_direction24,39),main="",start=90,clockwise=TRUE,lwd=3,line.col=c(outsum$col,"darkgreen"),
	labels=c("N","NE","E","SE","S","SW","W","NW"),label.pos=c(0,45,90,135,180,225,270,315),show.grid.labels=F)
legend(280,550,lwd=3,col=c("tomato3", "royalblue3", "darkgreen"), legend=c("MUPE long", "MUPE short", "HEPE"), bty = "n")
dev.off()


### abandoned approaches
plot(circular(outsum$dept_direction12[outsum$type=='short'],type = "directions", units = "degrees", rotation = "clock"), zero=pi/2, rotation = "clock",pch=16, col=outsum$col[outsum$type=='short'],cex=2,shrink=1.5, tol=0, tcl.text=-0.4)
par(new=T)
plot(circular(outsum$dept_direction12[outsum$type=='long'],type = "directions", units = "degrees", rotation = "clock"), zero=pi/2, rotation = "clock",pch=16, col=outsum$col[outsum$type=='long'], cex=2,shrink=1.4, tol=0, tcl.text=-0.4, axes=F, control.circle=circle.control(col = 0))
legend(-1.8,1.6,pch=16, col=c("tomato3", "royalblue3", "darkgreen"), legend=c("MUPE long", "MUPE short", "HEPE"), cex=1.3, bty = "n")

#### rose diagram does not allow overlay of different colours
rose.diag(circular(outsum$dept_direction12[outsum$type=='short'],type = "directions", units = "degrees", rotation = "clock"), bins = 16,shrink=1.5,prop=2,col="lightblue",axes=T,tcl.text=-0.2, cex.lab=1.5,zero=c(rad(90)),ticks=TRUE)		#homing direction








#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RELATE MASS GAIN TO TRIP DISTANCE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Henderson\\Analyses\\Seabirds")
setwd("A:\\RSPB\\UKOT\\Henderson\\Analyses\\Seabirds")

pdf("MUPE_mass_trip_dist.pdf", width=10, height=8)
par(mfrow=c(2,2), cex=1.3, mar=c(4,4,1,0), oma=c(2,0,0,0))
plot(prop_mass_gain~total_dist, outsum, pch=16, col=outsum$col,frame=F, xlab="", ylab="prop mass gain (%)")
abline(lm(prop_mass_gain~total_dist, outsum))
plot(mass_gain~total_dist, outsum, pch=16, col=outsum$col,frame=F, xlab="", ylab="mass gain (g)")
abline(lm(mass_gain~total_dist, outsum))
plot(dept_mass~total_dist, outsum, pch=16, col=outsum$col,frame=F, xlab="", ylab="departure mass (g)")
abline(lm(dept_mass~total_dist, outsum))
plot(return_mass~total_dist, outsum, pch=16, col=outsum$col,frame=F, xlab="", ylab="return mass (g)")
abline(lm(return_mass~total_dist, outsum))
mtext("total trip distance (km)", side=1, outer=T, cex=1.5)
dev.off()




pdf("MUPE_mass_trip_duration.pdf", width=10, height=8)
par(mfrow=c(2,2), cex=1.3, mar=c(4,4,1,0), oma=c(2,0,0,0))
plot(prop_mass_gain~duration, outsum, pch=16, col=outsum$col,frame=F, xlab="", ylab="prop mass gain (%)")
abline(lm(prop_mass_gain~duration, outsum))
plot(mass_gain~duration, outsum, pch=16, col=outsum$col,frame=F, xlab="", ylab="mass gain (g)")
abline(lm(mass_gain~duration, outsum))
plot(dept_mass~duration, outsum, pch=16, col=outsum$col,frame=F, xlab="", ylab="departure mass (g)")
abline(lm(dept_mass~duration, outsum))
plot(return_mass~duration, outsum, pch=16, col=outsum$col,frame=F, xlab="", ylab="return mass (g)")
abline(lm(return_mass~duration, outsum))
mtext("total trip duration (hrs)", side=1, outer=T, cex=1.5)
dev.off()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# POSTER GRAPH
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
jpeg("MUPE_poster_plot.jpg", width=2000, height=1200, quality = 100)
par(mar=c(20,20,6,5), oma=c(1,1,1,1))
plot(return_mass~total_dist, outsum, pch=16, col=outsum$col,frame=F, ylim=c(350,550), xlim=c(4000,16000),xlab="total trip distance (km)", ylab="return mass (g)", cex=7, axes=F,cex.lab=8, mgp=c(15,7,0))
axis(1,at=c(4000,7000,10000,13000,16000) , cex.axis=7, cex.lab=8, lwd=8, mgp=c(15,7,0))
axis(2,at=c(350,400,450,500,550), cex.axis=8, cex.lab=7, lwd=8, mgp=c(15,3,0), las=1)
abline(lm(return_mass~total_dist, outsum), lwd=8)
dev.off()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RELATE SPEED TO TRIP DISTANCE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pdf("MUPE_travel_speed.pdf", width=10, height=6)
par(mfrow=c(1,2), cex=1.3, mar=c(4,4,1,0), oma=c(0,2,0,0))
plot(travel_speed~total_dist, outsum, pch=16, col=outsum$col,frame=F, xlab="trip distance (km)", ylab="")
abline(lm(travel_speed~total_dist, outsum))
plot(travel_speed~duration, outsum, pch=16, col=outsum$col,frame=F, xlab="trip duration (hrs)", ylab="")
abline(lm(travel_speed~duration, outsum))
mtext("travel speed (km/h)", side=2, outer=T, cex=1.5)
dev.off()













#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE TRAVEL DISTANCES, SPEED AND TIME INTERVALS FOR EACH ANIMAL TRACK
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
head(output)

alltrips<-as.numeric(unique(Trips@data$trip_id[Trips@data$trip_id>0]))
migration<-data.frame()

for (a in alltrips){

input<-Trips@data[Trips@data$trip_id == a,]
input<-input[order(input$DateTime),]
head(input)
Nest<-paste("MUPE",input$ID[1],sep="")
input$Nest<-Nest
input$step_dist<-0
input$home_dist<-0
input$cumul_dist<-0
input$time_diff<-0
input$speed<-0
#first<-SpatialPoints(data.frame(birds$MaxOfLatitude[birds$Nest_Nr==a], birds$MaxOfLatitude[birds$Nest_Nr==a]), proj4string=CRS("+proj=longlat + datum=wgs84"))
first<-SpatialPoints(data.frame(birds$Nest_Lat[birds$Nest==Nest], birds$Nest_Long[birds$Nest==Nest]), proj4string=CRS("+proj=longlat + datum=wgs84"))

for (l in 2: dim(input)[1]){
input$time_diff[l]<-as.numeric(difftime(input$DateTime[l],input$DateTime[l-1], units="hours"))
fromloc<-SpatialPoints(data.frame(input$Longitude[l-1], input$Latitude[l-1]), proj4string=CRS("+proj=longlat + datum=wgs84"))
toloc<-SpatialPoints(data.frame(input$Longitude[l], input$Latitude[l]), proj4string=CRS("+proj=longlat + datum=wgs84"))
input$step_dist[l]<-spDistsN1(fromloc, toloc, longlat=T)
input$home_dist[l]<-spDistsN1(first, toloc, longlat=T)
input$cumul_dist[l]<-sum(input$step_dist)
input$speed[l]<-input$step_dist[l]/input$time_diff[l]
}

migration<-rbind(migration, input)

}




#### PLOT THE INSTANTANEOUS TRAVEL SPEED OVER TIME ###
head(migration)
pdf("MUPE_instantaneous_speed.pdf", width=12, height=18)
par(mfrow=c(11,4), mar=c(3,3,0,0), oma=c(2,2,0,0))
for (a in alltrips){
input<-migration[migration$trip_id == a,]
plot(speed~home_dist, input, type="l", ylim=c(0,80), frame=F, xlab="", ylab="")
text(a, x=0, y=75)
plot(speed~DateTime, input, type="l", ylim=c(0,80), frame=F, xlab="", ylab="")
plot(speed~cumul_dist, input, type="l", ylim=c(0,80), frame=F, xlab="", ylab="")
plot(speed~Longitude, input, type="l", ylim=c(0,80), frame=F, xlab="", ylab="")
}
mtext("instantaneous travel speed (km/h)", side=2, outer=T, cex=1.5)
mtext("Dist to nest                      Date                         Cumul distance                     Longitude", side=1, outer=T, cex=1.5)
dev.off()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IDENTIFY AREAS OF FORAGING BEHAVIOUR ALONG EACH TRIP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(EMbC)

## does not work on MoveStack but only on List of 'move' objects
migration<-migration[order(migration$ID, migration$DateTime),]
MOBJ_all<-list()
for (a in 1:length(alltrips)){
input<-migration[migration$trip_id == alltrips[a],]
MOBJ_all[[a]]<-move(x=input$Longitude, y=input$Latitude, time=input$DateTime, proj=CRS("+proj=longlat +ellps=WGS84"))
#MOBJ_all<-spTransform(x=MOBJ_all, CRSobj="+proj=aeqd", center=TRUE)			### not clear whether EMbC works on lat/long info or needs UTM
}


### BEHAVIOUR ASSIGNMENT FOR A SINGLE TRACK ###
x<-move(x=input$Longitude, y=input$Latitude, time=input$DateTime, proj=CRS("+proj=longlat +ellps=WGS84"))
mybcp <-stbc(x,info=-1)
stts(mybcp)
sctr(mybcp)
view(mybcp, lims=c(0,5000))
pkml(mybcp, display=T)


### BEHAVIOUR ASSIGNMENT FOR ALL TRACKS TOGETHER ###
mybcp <-stbc(MOBJ_all,info=-1)
STATES<-stts(mybcp)
sctr(mybcp)



### CAPTURE STATE ASSIGNMENT ###
migration$STATE<-NA
pdf("MUPE_foraging.pdf", width=16, height=12)
par (mfrow=c(3,4))
for(a in 1:length(alltrips)){
pkml(mybcp@bCS[[a]], display=F)
view(mybcp@bCS[[a]], lims=c(0,5000))
migration$STATE[migration$trip_id == alltrips[a]]<-mybcp@bCS[[a]]@A
}
dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MATCH BEHAVIOURAL ASSIGNMENT WITH GLS-IMMERSION DATA 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
getwd()
immersion<-read.table("A:\\RSPB\\UKOT\\Henderson\\Data\\Seabirds\\PathTrack_GPS_data\\MUPE_GPS_activity.csv", sep=",", header=T)
migration<-merge(migration,immersion[,c(2,3,7:9,16:18)], by=c("ID", "GPS_ID","Latitude","Longitude","TrackTime"), all.x=T)
migration$behav<-ifelse(migration$STATE==1,"rest",ifelse(migration$STATE==2,"forage",ifelse(migration$STATE==3,"travel","other")))
migration$behav<-as.factor(migration$behav)
boxplot(PropWet~behav, migration, range=0, outline=F)
plotdat<-aggregate(PropWet~behav, migration, FUN=mean)
plotdat$sd<-aggregate(PropWet~behav, migration, FUN=sd)[,2]
ucl<-plotdat$PropWet+0.5*plotdat$sd
lcl<-plotdat$PropWet-0.5*plotdat$sd
par(mar=c(5,5,0,2.5))
errbar(1:4,plotdat$PropWet, lcl, ucl, xlab="inferred behaviour",ylab="Prop Wet", axes=F, cex=1.5, cex.lab=1.8)
axis(1, at=c(0,1,2,3,4,5,6), labels=c("",plotdat$behav,""), cex.axis=1.5, cex=1.5, cex.lab=1.5)
axis(2, at=seq(0,1,0.1), labels=T, las=1, cex=1.5, cex.lab=1.5, cex.axis=1.5)
kruskal.test(PropWet~behav, migration)
kruskal.test(NoWetBout~behav, migration)
hist(migration$PropWet, 50)
hist(migration$NoWetBoutHr, 50)
hist(migration$NoWetBout, 50)
NoWetBoutHr

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DOES THE BEHAVIOUR DIFFER AMONG LONG AND SHORT TRIPS?
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sink("MUPE_behav_states.txt")
stts(mybcp)
sink()

STATES<-read.table("MUPE_behav_states.txt", header=T, sep="\t")
str(STATES)


behav<-data.frame()

for(a in 1:length(alltrips)){
input<-migration[migration$trip_id == alltrips[a],]
input$count<-1
bx<-aggregate(count~STATE, input, FUN=sum)
bx$time<-aggregate(time_diff~STATE, input, FUN=sum)[,2]
bx$Nest<-input$Nest[1]
bx$trip<-alltrips[a]
behav<-rbind(behav,bx)
}


behav<-merge(behav, outsum[,c(1,11,3,7,8,14,15,17:20, 25,26,27)], by=c("Nest","trip"))
behav$prop_time<-behav$time/behav$duration
head(behav)


### TEST FOR DIFFERENCES IN TIME SPENT IN CERTAIN BEHAVIOUR ###

forage<-behav[behav$STATE==2,]
wilcox.test(prop_time~type, forage)


pdf("MUPE_prop_foraging.pdf", width=10, height=8)

par(mfrow=c(2,2), cex=1.3, mar=c(4,4,1,0), oma=c(0,0,0,0))
plot(prop_time~total_dist, forage, pch=16, col=forage$col,frame=F, xlab="trip distance", ylab="proportion of trip foraging")
abline(lm(prop_time~total_dist, forage))
boxplot(prop_time~type, forage)

plot(time~total_dist, forage, pch=16, col=forage$col,frame=F, xlab="trip distance", ylab="total foraging time (hrs)")
abline(lm(time~total_dist, forage))
boxplot(time~type, forage)

dev.off()




### TEST FOR DIFFERENCES IN TIME SPENT IN CERTAIN BEHAVIOUR ###


par(mfrow=c(2,2), cex=1.3, mar=c(4,4,1,0), oma=c(0,0,0,0))
plot(prop_time~prop_mass_gain, forage, pch=16, col=forage$col,frame=F, xlab="trip distance", ylab="proportion of trip foraging")
abline(lm(prop_time~prop_mass_gain, forage))

plot(time~prop_mass_gain, forage, pch=16, col=forage$col,frame=F, xlab="trip distance", ylab="total foraging time (hrs)")
abline(lm(time~prop_mass_gain, forage))

plot(prop_time~prop_mass_gain, forage, pch=16, col=forage$col,frame=F, xlab="mass gain", ylab="")
abline(lm(prop_time~prop_mass_gain, forage))

plot(prop_time~duration, forage, pch=16, col=forage$col,frame=F, xlab="duration", ylab="")
abline(lm(prop_time~duration, forage))

mtext("proportion of time spent foraging", side=2, outer=T, cex=1.5)
dev.off()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SIMPLE SUMMARY FOR POSTER
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

postertab<-aggregate(duration~type, outsum, FUN=mean)
postertab$duration<-postertab$duration/24
postertab$duration.sd<-aggregate(duration~type, outsum, FUN=sd)[,2]/24

postertab$max_dist<-aggregate(max_dist~type, outsum, FUN=mean)[,2]
postertab$max_dist.sd<-aggregate(max_dist~type, outsum, FUN=sd)[,2]

postertab$tot_dist<-aggregate(total_dist~type, outsum, FUN=mean)[,2]
postertab$tot_dist.sd<-aggregate(total_dist~type, outsum, FUN=sd)[,2]

postertab$prop_mass_gain<-aggregate(prop_mass_gain~type, outsum, FUN=mean)[,2]
postertab$prop_mass_gain.sd<-aggregate(prop_mass_gain~type, outsum, FUN=sd)[,2]

write.table(postertab,"clipboard", sep="\t", row.names=F)

aggregate(dept_mass~type, outsum, FUN=mean)
aggregate(dept_mass~type, outsum, FUN=sd)

aggregate(return_mass~type, outsum, FUN=mean)
aggregate(return_mass~type, outsum, FUN=sd)

aggregate(time~type, forage, FUN=mean)
aggregate(time~type, forage, FUN=sd)

aggregate(prop_time~type, forage, FUN=mean)
aggregate(prop_time~type, forage, FUN=sd)




