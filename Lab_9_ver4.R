####################################################################################################################################
# WILD 562 - R code for LAB 8 - Step Selection Functions
####################################################################################################################################

### Mark Hebblewhite
### Daniel Eacker
### Adapted from code written
#'author: Björn Reineking 
#'email: bjoern.reineking@irstea.fr

####' Step selection functions: Literature
#' ========================================================
#' - Forester et al. (2009) Ecology 90: 3554-3565
#' - Potts et al. (2014) J. R. Soc. Interface 11: 20140333
#' - Avgar et al. (2016) Methods in Ecology and Evolution 7: 619-630

#function to install and load required packages
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

#load or install these packages:
packages <- c("sp","raster","rgdal","ggmap","survival","TwoStepCLogit","pbs","dplyr","lubridate","move","MASS","fdrtool","circular","CircStats", "coxme", "mclogit","data.table","scales")


#run function to install packages
ipak(packages)

##### 0.1 Preliminaries: setting working directory #######################################################################

## define working directory on each computer
wd_laptop <- "C:/Users/danea/OneDrive/Documents/Archive (1)/Lab8_new"
wd_desktop <- "C:/Users/danea/OneDrive/Documents/Archive (1)/Lab8_new"

## automatically set working directory depending which computer you're on
ifelse(file.exists(wd_laptop), setwd(wd_laptop), setwd(wd_desktop)) 
getwd()
list.files() ## handy command to see what is inside the working directory

####################################################################################################################################
# Let's just work with elk 29 to start out with, and then we will do 
#the entire dataset

# read in csv file
yl29 <- read.csv("yl29.csv", header=T)

# inspect the dataset
head(yl29)

## Plot the data
ggplot(yl29, aes(x=UTMX, y = UTMY)) + geom_point()
ggplot(yl29, aes(x=Date, y= UTMY)) + geom_point()
ggplot(yl29, aes(x=Date, y= UTMX)) + geom_point()
### What date did this elk migrate?

####################################################################################################################################
# load custom functions from BjÃ¶rn Reineking 
source("ssf_fun_2016_final.R", verbose = FALSE)

####################################################################################################################################
# create "move" object using move package

# first format time stamp to bring date and time together

# format Date as date class
yl29$Date <- as.Date(yl29$Date, format="%m/%d/%y")

# make Time variable to standardize number of digits
yl29$Time <- as.character(yl29$Time)

new.time <- numeric(length(yl29$Time))

for(i in 1:length(yl29$Time)){
  if(nchar(substr(yl29$Time[i], 1, 5))==4) {new.time[i] <- paste("0",yl29$Time[i],sep="")}else
    if(nchar(substr(yl29$Time[i], 1, 5))==5) {new.time[i] <- yl29$Time[i]}
}

# now paste together Date and Time to make timestamp and convert to POSIXct format
yl29$rtime <- as.POSIXct(paste(yl29$Date,new.time, sep=" "),  format="%Y-%m-%d %H:%M", tz="America/Denver")

### Make the yl29 data a Move object or MoveStack (multiple Move objects, i.e., animals)
### ?move
yl29.mv <- move(x=yl29$UTMX, y=yl29$UTMY, 
                time=yl29$rtime, 
                data=yl29, proj=CRS("+proj=utm +zone=11 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"), animal=yl29$ELKUID)

# visualize the movement data
par(mfrow=c(1,1))
plot(yl29.mv) ## note you can do this with the Move Object

# look at sumary of move object
summary(yl29.mv)

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# THIS PART OF THE ANALYSIS IS FOR THINNING THE DATA BASED ON A MAXIMUM TIME STEP AND CREATING A moveBurst object that
# can be used in a step selection function = I have modified some of the code Björn Reineking to make this work
# as his code only deals with "Move" or "MoveStack" objects

# set maximum step duration
step <- 3

# burst track and create burstId variable with long (i.e., >=3 hour timelag) and normal (i.e., <3 hour timelag) segments
bursted<-move::burst(yl29.mv, c('normal','long')[1+(timeLag(yl29.mv, units='hours')>step)])

head(bursted@burstId)
#[1] normal normal normal normal normal normal
#Levels: long normal

# The projection seems to matter just slightly whether you work in UTMs or latlong, but the angles are the same when calculated
#bursted2 <- spTransform(bursted, CRS("+proj=longlat +ellps=GRS80 +datum=NAD83 +towgs84=0,0,0"))

# create empty list										
burstlyl29 <- vector("list",(nrow(bursted)-1))

# create list of bursted segments (note you need two locations for this)
for(i in 2:nrow(bursted)){

		burstlyl29[[(i-1)]] <- bursted[(i-1):i]
									
				
}

# create empty index to hold indicator variable of whether a track was of normal or long timestep
# We use this to filter the data latter
index1 <- numeric(nrow(bursted)-1)

# loop over index1 using ifelse statements to create indicator variable
for(i in 1:(nrow(bursted)-1)){

	if(burstlyl29[[i]]@burstId=="normal"){index1[i] <- 1}else
	if(burstlyl29[[i]]@burstId=="long"){index1[i] <- 0}
	}


#This function seems correct (a negative of an angle should be the same
#radians as the positive of the angle just with an opposite sign.)
## The function is currently defined as
degrees2radians <- function (x) 
{
  radians <- (x * (pi/180))
  return(radians)
}

# calculate distances and turn angles (relative as well in radians) for ALL data
	dist <- distance(yl29.mv)
	yl29.mv <- spTransform(yl29.mv, CRS = "+proj=longlat +ellps=GRS80 +datum=NAD83 +towgs84=0,0,0")
	abs.angle <- angle(yl29.mv)
    abs.angle3 <- degrees2radians(abs.angle)

# calculate turn angles by using diff function to find difference between consecutive angles
ta.data <- data.frame(index1, c(NA, diff(abs.angle3)))
names(ta.data) <- c("index1", "turn.angle3")

# set long bursts to NA using index1 created above where index1==0
for(i in 1:(length(ta.data$turn.angle3))){

 if(ta.data$index1[i]==0){ta.data$turn.angle3[i] <- NA}else
 if(ta.data$index1[i]==1){ta.data$turn.angle3[i] <- ta.data$turn.angle3[i]}

}

# set NA for turn angle for  after missing burst (the next turn angle cannot be calcuted with a missing before between it)
indext1 <- which(is.na(ta.data$turn.angle3)==T)+1
indext1 <- indext1[2:207]

ta.data$turn.angle3[indext0] <- NA
ta.data$turn.angle3[indext1] <- NA

# put a place holder of 0 to hold first row of data
ta.data$turn.angle3[1] <- 0

# set distances to NA where turn angle is not available
ta.data$dist <- ifelse(is.na(ta.data$turn.angle3)==F,dist,NA)

# set absolute angles to NA where turn angle is not available
ta.data$abs.angle3 <- ifelse(is.na(ta.data$turn.angle3)==F,abs.angle3,NA)



# remove na's before fitting distributions to distances and relative turn angles
rel.angle3 <- ifelse(abs(na.omit(ta.data$turn.angle3)) > pi, sign(na.omit(ta.data$turn.angle3)) * (abs(na.omit(ta.data$turn.angle3)) - 2 * pi), na.omit(ta.data$turn.angle3))
dist <- na.omit(ta.data$dist)
abs.angle3 <- na.omit(ta.data$abs.angle3)


# note that 1 radian is just under 57.3 degrees
# note that this function does not seem to work right
#degrees2radians <-function(x) {
#  -2*pi*(x-90)/360
#}

# I also think this is unnecesary; the resulting degrees2 radians is within
#within pi and -pi
# wrap_rel_angle <- function(x) {
# wraps relative angles to the range from -pi to pi
#  ifelse(abs(x) > pi, sign(x) * (abs(x) - 2 * pi), x)
#}


# basic summary of angles (in degrees) and distances (in meters) between successive locations
angle_dist_summary <- data.frame(cbind(summary(dist), summary(rel.angle3)))
names(angle_dist_summary) <- c("distance", "relative angle (radians)")
angle_dist_summary

# examine empirical distances and turning angles
par(mfrow = c(1, 2))
hist(dist, breaks = 20, main = "", xlab = "Distance (m)")
hist(rel.angle3, breaks = seq(-pi, pi,len=11), main="", xlab="Relative angle (radians)")


# Fit distributions to distances
fexp <- fitdistr(dist, "exponential")
fgam <- fitdistr(dist, "gamma")
flogn <- fitdistr(dist, "lognormal")

AICout <- AIC(flogn, fgam, fexp)
AICout <- cbind(AICout, c(0, diff(AICout$AIC)))
names(AICout) <- c("df", "AIC", "deltaAIC")
AICout

## So logn seems to be the best fitting function

par(mfrow = c(1,1))
hist(dist, breaks = 50, prob = TRUE, xlim = c(0, 8000), ylim = c(0, 2e-3),xlab = "Step length (m)", main = "")
plot(function(x) dexp(x, rate = fexp$estimate), add = TRUE, from = 0, to = 5000)
plot(function(x) dexp(x, rate = fexp$estimate/2), add = TRUE, from = 0, to = 5000, col="purple")
plot(function(x) dgamma(x, shape = fgam$estimate["shape"], rate = fgam$estimate["rate"]), add = TRUE, 
     from = 0, to = 5000, col = "red")
plot(function(x) dlnorm(x, meanlog = flogn$estimate["meanlog"], 
                        sdlog = flogn$estimate["sdlog"]), add = TRUE, 
     from = 0, to = 5000, col = "blue")
plot(function(x) dhalfnorm(x, theta = 1/mean(dist)), add = TRUE, 
     from = 0, to = 5000, col = "green")
legend("topright", lty = 1, col = c("black", "red", "blue", "green"),
       legend = c("exp","exp/2" ,"gamma", "lnorm", "halfnorm"))

# Fit distributions to angles
fkappa <- est.kappa(rel.angle3)
fkappa

hist(rel.angle3, prob = TRUE, main = "", xlab = "Turning angle")
# use distribution von Mi
plot(function(x) dvonmises(x, circular(0), kappa=fkappa), add = TRUE, from = -pi, to = pi, col = "red")

####################################################################################################################################
# subset bursted list for bursts that have turn angle data

# index for turn angle does not equal NA
sub.index <- which(is.na(ta.data$turn.angle3)==F)

# return first burst turn angle to NA (it cannot be estimated)
ta.data$turn.angle3[1] <- NA
rel.angle3[1] <- NA

# subset bursted list
burstlyl29.2 <- burstlyl29[sub.index]

#check timeLags
timeLags <- numeric(length(burstlyl29.2))

# calculate timeLags by looping over 
for(i in 1:length(burstlyl29.2)){
	timeLags[i] <- timeLag(burstlyl29.2[[i]], units="hours")

}	

# inspect timeLags to make sure < than maximum step
timeLags

# look at time lag between successive locations
hist(timeLags)
summary(timeLags)


####################################################################################################################################
# Objective 2.0 Create match case control data points 

utm_crs <-CRS("+proj=utm +zone=11 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

#Control steps: sample with one of the following distance distributions: exponential, 
#gamma, lognormal, halfnormal, potentially drawing turning angles from a von Mises 
#distribution. We will use an lognormal distribution here because it seems to fit the 
#data well (the gamma and exponential seem to fit pretty close as well). We will use a radially symmetric proposal 
#distribution because it is easier to interpret parameter estimates for 
#cos(rel.angle), and because in an earlier analysis on these data it showed less 
#residual autocorrelation (note that this was not checked for the elk dataset).	 

# create empty data frame to hold moveBursts in 
newdata <- as.data.frame(matrix(nrow=(length(burstlyl29.2)),ncol=9))
names(newdata) <- names(burstlyl29.2[[1]]@data)
newdata$Date <- as.Date(newdata$Date)
newdata$rtime <- as.POSIXct(newdata$rtime,  format="%Y-%m-%d %H:%M:%S", tz="America/Denver")
newdata$individual.local.identifier <- as.character(newdata$individual.local.identifier)
newdata$sensor <- as.character(newdata$sensor)
newdata$Time <- as.character(newdata$Time)

#make new dataset with just 1st coordinate of each moveBurst
for(i in 1:length(burstlyl29.2)){

newdata[i,] <-    as.data.frame(burstlyl29.2[[i]]@data[1,])

}

# Mannual add data from colums 8 and 9
newdata[,8] <- rep("X29", nrow(newdata))
newdata[,9] <- rep("unknown", nrow(newdata))

# add in turn angles and distances
newdata$dist <- dist
newdata$rel.angle3 <- rel.angle3
newdata$abs.angle3 <- abs.angle3 

####################################################################################################################################
# split data and prepare sampling loop to create match case

# split data into list based on ID2 variable
data_split <- split(newdata, newdata$ID2)

# function to create ssf with MoveBursts; data input is a list of locations; note that the variables "ID2", "SORTID", "UTMX", "UTMY",
#"abs.angle3 are used in this function that is specific to the dataset. You will have to customize the function a bit
# if these variables do not exist or change names to these varaibles in the function below

prepare_ssf_burst <- function(data=NULL, K=NULL, theta=NULL, kappa=NULL, distr = "lognormal"){
		
		# create empty lists to hold results of for loop
		ssf_data1 <- vector("list",(length(data)))
		temp_data <- vector("list",(length(data)))
		temp_data2 <- vector("list",(length(data)))
		temp_data3 <- vector("list",(length(data)))
		data2 <- data

for(i in 2:length(data2)){  
			
			#sample turning angle
			r <- switch(distr,
            exponential = rexp(K, theta),
            gamma = rgamma(K, theta[1], theta[2]),
            lognormal = rlnorm(K, theta[1], theta[2]),
            halfnormal = abs(rnorm(K, 0, sqrt(pi/2)/theta)))
			start_angles <- rep(data2[[i]]$abs.angle3,K)
			theta2 <- rvonmises(K, circular(0), kappa) + start_angles
			dx <- r * cos(theta2)
			dy <- r * sin(theta2)
			temp_data[[i]] <- data.frame(cbind(dx, dy) )
			temp_data2[[i]] <- rbind(cbind(data2[[i]]$UTMX, data2[[i]]$UTMY), 
			cbind(temp_data[[i]][,1] + data2[[i]]$UTMX, temp_data[[i]][,2] + data2[[i]]$UTMY))
			temp_data3[[i]] <- cbind(temp_data[[i]][,1] + data2[[i]]$UTMX, temp_data[[i]][,2] + data2[[i]]$UTMY)
			abs.angle <- atan2(temp_data3[[i]][,2] - data2[[i]]$UTMY, temp_data3[[i]][,1] - data2[[i]]$UTMX)
			abs.angle2 <- (abs.angle - data2[[i-1]]$abs.angle3) # make i-1
			rel.angle3 <- ifelse(abs(abs.angle2) > pi, sign(abs.angle2) * (abs(abs.angle2) - 2 * pi), abs.angle2)
			ssf_data1[[i]]<- data.frame(rep(data2[[i]]$individual.local.identifier,K+1),temp_data2[[i]][,1],temp_data2[[i]][,2],
			rep(data2[[i]]$rtime,K+1),
			rep(paste(data2[[i]]$individual.local.identifier,"_",data2[[i]]$ID2,"_",
			data2[[i]]$SORTID,sep=""),K+1),
			c(1, rep(0, K)),
			c(data2[[i]]$dist, r),
			c(data2[[i]]$rel.angle3,rel.angle3))
			
			}
		
			ssf_data2 <- rbindlist(ssf_data1)	
			class(ssf_data2) <- "data.frame"
			names(ssf_data2) <- c("id","Easting","Northing","rtime","stratum","used","dist","rel.angle")
			ssf_data<-SpatialPointsDataFrame(ssf_data2[,c("Easting","Northing")], ssf_data2[,c(1:ncol(ssf_data2))])
			proj4string(ssf_data) <- utm_crs
			return(ssf_data)
}

set.seed(5)
# note that data much be in a list of MoveBursts
ssf_out <- prepare_ssf_burst(data=data_split, K=5, theta=flogn$estimate, kappa=fkappa, dist = "lognormal")


####################################################################################################################################
# check that available points make sense

dat0 <- subset(ssf_out, stratum == "X29_903_2" & used == 0)
dat1 <- subset(ssf_out, stratum == "X29_903_2" & used == 1)

plot(dat0, axes=TRUE)
plot(dat1, add=T, col="green")
points(subset(ssf_out, stratum == "X29_904_3" &
                used == 1), col ="blue")
points(subset(ssf_out, stratum == "X29_905_4" &
                used == 1), col ="red")


dat2 <- as.data.frame(subset(ssf_out, stratum %in% c("X29_903_2", "X29_904_3") &
                               used == 1))
lines(dat2$x, dat2$y, col ="green")
legend("topright", col = c("green", "blue", "red"), pch = 1, 
       legend = c("t-1", "t", "t+1"))


####################################################################################################################################
####################################################################################################################################
# Create ssf for Move object using all of the data for comparison

angle_dist <- prepare_angle_dist(yl29.mv)
# Inspect the data
summary(angle_dist)

dist <- angle_dist[,1]
rel.angle3 <- angle_dist[,2]

# examine empirical distances and turning angles
par(mfrow = c(1, 2))
hist(dist, breaks = 20, main = "", xlab = "Distance (m)")
hist(rel.angle3, breaks = seq(-pi, pi,len=11), main="", xlab="Relative angle (radians)")

# Fit distributions to distances
fexp <- fitdistr(dist, "exponential")
fgam <- fitdistr(dist, "gamma")
flogn <- fitdistr(dist, "lognormal")

AICout <- AIC(flogn, fgam, fexp)
AICout <- cbind(AICout, c(0, diff(AICout$AIC)))
names(AICout) <- c("df", "AIC", "deltaAIC")
AICout

## So logn seems to be the best fitting function

par(mfrow = c(1,1))
hist(dist, breaks = 50, prob = TRUE, xlim = c(0, 8000), ylim = c(0, 2e-3),xlab = "Step length (m)", main = "")
plot(function(x) dexp(x, rate = fexp$estimate), add = TRUE, from = 0, to = 5000)
plot(function(x) dexp(x, rate = fexp$estimate/2), add = TRUE, from = 0, to = 5000, col="purple")
plot(function(x) dgamma(x, shape = fgam$estimate["shape"], rate = fgam$estimate["rate"]), add = TRUE, 
     from = 0, to = 5000, col = "red")
plot(function(x) dlnorm(x, meanlog = flogn$estimate["meanlog"], 
                        sdlog = flogn$estimate["sdlog"]), add = TRUE, 
     from = 0, to = 5000, col = "blue")
plot(function(x) dhalfnorm(x, theta = 1/mean(dist)), add = TRUE, 
     from = 0, to = 5000, col = "green")
legend("topright", lty = 1, col = c("black", "red", "blue", "green"),
       legend = c("exp","exp/2" ,"gamma", "lnorm", "halfnorm"))

# Fit distributions to angles
fkappa <- est.kappa(rel.angle3)
fkappa

hist(rel.angle3, prob = TRUE, main = "", xlab = "Turning angle")
# use distribution von Mi
plot(function(x) dvonmises(x, circular(0), kappa=fkappa), add = TRUE, from = -pi, to = pi, col = "red")

# set projection for estimation
utm_crs <- "+proj=utm +zone=11 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"

####################################################################################################################################
# run prepare_ssf_steps on elk yl29.mv

# set the number of control locations
my_K <- 5

# set random number generating seed
set.seed(5)

# use prepare_ssf_steps to create the step selection function data
ssf_data3 <- prepare_ssf_steps(yl29.mv, 
                                    method = "lognormal", K = my_K, 
                                    theta = flogn$estimate, kappa = fkappa,
                                    crs = utm_crs)
	

#set projection to UTMs
ssf_data4 <- spTransform(ssf_data3 , CRS("+proj=utm +zone=11 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
	

####################################################################################################################################
# make plot of thinned and unthinned data for elk yl29			
							
tiff("ssf_yl29_bursted_all4.tiff", res=600, compression = "lzw", height=5, width=10, units="in")
par(mfrow=c(1,2))

used1 <- ssf_out[ssf_out$used==1,]
used0 <- ssf_out[ssf_out$used==0,]
plot(used0, col="green",pch = 1,main="thinned",axes=T)
plot(used1, col="red", add=T, ,pch = 1)
legend("topleft", legend=c("available", "used"), col = c("green","red"),pch = 1) 

used1 <- ssf_data4[ssf_data4$used==1,]
used0 <- ssf_data4[ssf_data4$used==0,]
plot(used0, col="green",pch = 1, main="unthinned",axes=T)
plot(used1, col="red", add=T,pch = 1)
legend("topleft", legend=c("available", "used"), col = c("green","red"),pch = 1) 

dev.off()
	


	
####################################################################################################################################	####################################################################################################################################
####################################################################################################################################
## Create ssf for MoveStack - multiple individuals (note that this should follow the Björn Reineking code more closely since
# that code was created for a MoveStack

# Read in data for all elk
yl.all <- read.csv("yl_all.csv",header=T)

# look at data
head(yl.all)

#look at elk ids
levels(as.factor(yl.all$ElkID))

#create time stamp
# first format time stamp to bring date and time together

# format Date as date class
yl.all$Date <- as.Date(yl.all$Date, format="%m/%d/%Y")

# make Time variable to standardize number of digits
new.hour <- numeric(length(yl.all$Hour))

#first fix hours
for(i in 1:length(yl.all$Hour)){
  if(nchar(substr(yl.all$Hour[i], 1, 2))==1) {new.hour[i] <- paste("0",yl.all$Hour[i],sep="")}else
    if(nchar(substr(yl.all$Hour[i], 1, 2))==2) {new.hour[i] <- yl.all$Hour[i]}
}

# make Time variable to standardize number of digits
new.minute <- numeric(length(yl.all$Minute))

#first fix hours
for(i in 1:length(yl.all$Minute)){
  if(nchar(substr(yl.all$Minute[i], 1, 2))==1) {new.minute[i] <- paste("0",yl.all$Minute[i],sep="")}else
    if(nchar(substr(yl.all$Minute[i], 1, 2))==2) {new.minute[i] <- yl.all$Minute[i]}
}

hr_min <- paste(new.hour,":",new.minute,sep="")

# now paste together Date and Time to make timestamp and convert to POSIXct format
yl.all$rtime <- as.POSIXct(paste(yl.all$Date,hr_min,sep=" "),  format="%Y-%m-%d %H:%M", tz="America/Denver")

# inspect first 6 elements of rtime
head(yl.all$rtime)

# order by timestamps
yl.all <- yl.all[order(yl.all$rtime),]

# order by individuals
yl.all <- yl.all[order(yl.all$ElkID),]

# create move stack 
all.mv <- move(x=yl.all$UTMX, y=yl.all$UTMY, 
                time=yl.all$rtime, 
                data=yl.all, proj=CRS("+proj=utm +zone=11 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"), animal=yl.all$ElkID, removeDuplicatedTimestamps=F)

####################################################################################################################################
# get distance step length and turn angle summaries
angle_dist <- prepare_angle_dist(all.mv)
# Inspect the data
summary(angle_dist)

dist <- angle_dist[,1]
rel.angle3 <- angle_dist[,2]

# examine empirical distances and turning angles
par(mfrow = c(1, 2))
hist(dist, breaks = 20, main = "", xlab = "Distance (m)")
hist(rel.angle3, breaks = seq(-pi, pi,len=11), main="", xlab="Relative angle (radians)")

# Fit distributions to distances
fexp <- fitdistr(dist, "exponential")
fgam <- fitdistr(dist, "gamma")
flogn <- fitdistr(dist, "lognormal")

AICout <- AIC(flogn, fgam, fexp)
AICout <- cbind(AICout, c(0, diff(AICout$AIC)))
names(AICout) <- c("df", "AIC", "deltaAIC")
AICout

## So logn seems to be the best fitting function

par(mfrow = c(1,1))
hist(dist, breaks = 50, prob = TRUE, xlim = c(0, 8000), ylim = c(0, 2e-3),xlab = "Step length (m)", main = "")
plot(function(x) dexp(x, rate = fexp$estimate), add = TRUE, from = 0, to = 5000)
plot(function(x) dexp(x, rate = fexp$estimate/2), add = TRUE, from = 0, to = 5000, col="purple")
plot(function(x) dgamma(x, shape = fgam$estimate["shape"], rate = fgam$estimate["rate"]), add = TRUE, 
     from = 0, to = 5000, col = "red")
plot(function(x) dlnorm(x, meanlog = flogn$estimate["meanlog"], 
                        sdlog = flogn$estimate["sdlog"]), add = TRUE, 
     from = 0, to = 5000, col = "blue")
plot(function(x) dhalfnorm(x, theta = 1/mean(dist)), add = TRUE, 
     from = 0, to = 5000, col = "green")
legend("topright", lty = 1, col = c("black", "red", "blue", "green"),
       legend = c("exp","exp/2" ,"gamma", "lnorm", "halfnorm"))

# Fit distributions to angles
fkappa <- est.kappa(rel.angle3)
fkappa

hist(rel.angle3, prob = TRUE, main = "", xlab = "Turning angle")
# use distribution von Mi
plot(function(x) dvonmises(x, circular(0), kappa=fkappa), add = TRUE, from = -pi, to = pi, col = "red")

# set projection for estimation
utm_crs <- "+proj=utm +zone=11 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"

####################################################################################################################################
# run prepare_ssf_steps on elk yl29.mv

# set the number of control locations
my_K <- 5

# set random number generating seed
set.seed(5)

# use prepare_ssf_steps to create the step selection function data
ssf_data5 <- prepare_ssf_steps(all.mv, 
                                    method = "lognormal", K = my_K, 
                                    theta = flogn$estimate, kappa = fkappa,
                                    crs = utm_crs)
	

#set projection to UTMs
ssf_data6 <- spTransform(ssf_data5 , CRS("+proj=utm +zone=11 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))				

# create plot for all individuals

tiff("ssf_all.tiff", res=600, compression = "lzw", height=5, width=5, units="in")
par(mfrow=c(1,1))

used1 <- ssf_data6[ssf_data6$used==1,]
used0 <- ssf_data6[ssf_data6$used==0,]
plot(used0, col="green",pch = 1, main="all individuals",axes=T)
plot(used1, col="red", add=T,pch = 1)
legend("topleft", legend=c("available", "used"), col = c("green","red"),pch = 1) 

dev.off()

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
######## Objective 3.0 Build covariate dataset


#create an empty raster
mask.raster <- raster()

#set extent (note that I customized this extent so it covered both elc_habitat and humanacess)
extent(mask.raster) <- c(xmin=443680.6, xmax=650430.4, ymin=5618405, ymax=5789236) 	

#set the resolution to 30 m 
res(mask.raster)<-30
			
#match projection to elc_habitat shapefile
projection(mask.raster)<- utm_crs
			
#set all values of mask.raster to zero
mask.raster[]<-0

# check working directory
getwd()

# set working directory to read in other raster layers (I needed this to get the projection from a layer)
setwd("C:/Users/danea/OneDrive/Documents/Archive (1)/Lab8_new")
elevation<-raster("Elevation2.tif") #resampled
disthumanaccess2<-raster("DistFromHumanAccess2.tif") #resampled
disthhu<-raster("DistFromHighHumanAccess2.tif")
landcover<-raster("landcover16.tif") 

#  set the resolution to utm_crs above (same as ssf_data6)
elevation <- projectRaster(elevation, mask.raster)
disthumanaccess2 <- projectRaster(disthumanaccess2, mask.raster)
disthhu <- projectRaster(disthhu, mask.raster)
landcover <- projectRaster(landcover, mask.raster)

# check projection of elevation
elevation@crs@projargs

#  MAKE NDVI_ID variable
# check min and max dates of timestamp for ssf_data6 (a SpatialPointsDataFrame for all individuals ssf)
min(ssf_data6@data$date)
#"2003-04-15 02:00:00 MDT"

max(ssf_data6@data$date)
# "2003-10-14 20:00:00 MDT"

# extract Month from ssf_data6
ssf_data6@data$Month <- as.numeric(substr(ssf_data6@data$date, 6, 7))

# extract Day from ssf_data6
ssf_data6@data$Day <- as.numeric(substr(ssf_data6@data$date, 9, 10))

# create NDVI_ID 
ssf_data6@data$NDVI_ID<-with(ssf_data6@data, 

	ifelse(Month=="4" & Day<23, 7,
	ifelse(Month=="4" & Day>=23 | Month=="5" & Day<9, 8,
	ifelse(Month=="5" & Day>=9  & Day<25, 9,
	ifelse(Month=="5" & Day>=25 | Month=="6" & Day<10, 10,
	ifelse(Month=="6" & Day>=10  & Day<26, 11,
	ifelse(Month=="6" & Day>=26 | Month=="7" & Day<12, 12,
	ifelse(Month=="7" & Day>=12  & Day<28, 13,
	ifelse(Month=="7" & Day>=28 | Month=="8" & Day<13, 14,
	ifelse(Month=="8" & Day>=13  & Day<29, 15,
	ifelse(Month=="8" & Day>=29 | Month=="9" & Day<14, 16,
	ifelse(Month=="9" & Day>=14 & Day<30, 17,
	ifelse(Month=="9" & Day>=30 | Month=="10" & Day<16, 18,
	ifelse(Month=="10" & Day>=16 & Day<=31, 19,NA))))))))))))))

# examine NDVI_ID
table(ssf_data6@data$NDVI_ID)

ssf_data6@data$NDVI_ID2 <- as.character(ssf_data6@data$NDVI_ID)

#first NDVI_ID to make 2 digit
for(i in 1:length(ssf_data6@data$NDVI_ID)){
  if(nchar(substr(ssf_data6@data$NDVI_ID[i], 1, 2))==1) {ssf_data6@data$NDVI_ID2[i] <- paste("0",ssf_data6@data$NDVI_ID[i],sep="")}else
    if(nchar(substr(ssf_data6@data$NDVI_ID[i], 1, 2))==2) {ssf_data6@data$NDVI_ID2[i] <- ssf_data6@data$NDVI_ID[i]}
}


############################################################
#read in 2003 NDVI Raster layers (note these have been resampled to the resolution and extent and crs of the other rasters)

# MAKE SURE WORKING DIRECTORY IS POINTED AT NDVI FILES
getwd()

# read in NDVI rasters
NDVI_2003_07<-raster("NDVI_2003_07.tif") 
NDVI_2003_08<-raster("NDVI_2003_08.tif") 
NDVI_2003_09<-raster("NDVI_2003_09.tif") 
NDVI_2003_10<-raster("NDVI_2003_10.tif") 
NDVI_2003_11<-raster("NDVI_2003_11.tif") 
NDVI_2003_12<-raster("NDVI_2003_12.tif") 
NDVI_2003_13<-raster("NDVI_2003_13.tif") 
NDVI_2003_14<-raster("NDVI_2003_14.tif") 
NDVI_2003_15<-raster("NDVI_2003_15.tif") 
NDVI_2003_16<-raster("NDVI_2003_16.tif") 
NDVI_2003_17<-raster("NDVI_2003_17.tif") 
NDVI_2003_18<-raster("NDVI_2003_18.tif") 


# mask raster stack of NDVI tiles
NDVIstk <- stack(NDVI_2003_07, NDVI_2003_08, NDVI_2003_09, NDVI_2003_10, NDVI_2003_11, NDVI_2003_12, NDVI_2003_13, NDVI_2003_14,
	      NDVI_2003_15, NDVI_2003_16, NDVI_2003_17, NDVI_2003_18)

NDVIstk <- NDVIstk *0.0001		  
		  
# examine CRS of raster files 
NDVIstk@crs@projargs

# resample NDVI stack to extent and resolution and projection of other layers
NDVIstk <- projectRaster(NDVIstk, mask.raster)

# average over stacked NDVI tiles
meanNDVI <- mean(NDVIstk, na.rm=T)

# get values for meanNDVI
meanNDVI@data@values <- getValues(meanNDVI)

# look at mean NDVI values plot
plot(meanNDVI)

#now try stacking them
all_base_rasters <- stack(elevation, disthumanaccess2, disthhu, landcover2, meanNDVI)

# extract covariates
cov.ext<-extract(all_base_rasters, ssf_data6)
head(cov.ext)
str(cov.ext)

# add extracted covariates to ssf_data
ssf_data6@data$elevation <- cov.ext[,1]
ssf_data6@data$distha <- cov.ext[,2]
ssf_data6@data$disthha <- cov.ext[,3]
ssf_data6@data$landcover <- cov.ext[,4]
ssf_data6@data$meanNDVI <- cov.ext[,5]

# look at first 6 rows of data
head(ssf_data6@data)

# look at histogram of mean NDVI values
hist(ssf_data6@data$meanNDVI)
		  
###########################################################
# Extract NDVI values for time-dependent	

# split data by NDVI_ID
temp_data <- split(ssf_data6, ssf_data6@data$NDVI_ID)

# create empty list to hold output
timeNDVI <- vector("list",length(temp_data))

# extract time dependent variables by looping over data listed by NDVI period
for(i in 1:length(temp_data)){
	
	timeNDVI[[i]] <- extract(NDVIstk[[i]], temp_data[[i]])
	timeNDVI[[i]] <- as.data.frame(timeNDVI[[i]])

}

# okay now row bind list and create variable
temp <- as.vector(as.data.frame(rbindlist(timeNDVI)))
ssf_data6@data <- cbind(ssf_data6@data,temp[,1])
names(ssf_data6@data)[19] <- "timeNDVI"

# change ssf_data6 to ssf_data
ssf_data <- ssf_data6

### Add landcover categories to ssf_data # note that you shouldn't need this landcover.legend file
#setwd("/Users/mark.hebblewhite/Dropbox/WILD 562/Spring2017/Lab8/")
#landcover.legend <- read.csv("landcover_legend.csv", header = TRUE, sep = ",")
#landcover.legend 

ssf_data@data$landcover.cat = ifelse(ssf_data@data$landcover == 0, "NA", 
       ifelse(ssf_data@data$landcover == 1, "Open Conifer", 
          ifelse(ssf_data@data$landcover == 2, "Moderate Conifer", 
             ifelse(ssf_data@data$landcover == 3, "Closed Conifer", 
                ifelse(ssf_data@data$landcover == 4, "Deciduous", 
                   ifelse(ssf_data@data$landcover == 5, "Mixed", 
                      ifelse(ssf_data@data$landcover == 6, "Regen", 
                        ifelse(ssf_data@data$landcover == 7, "Herbaceous",                                   
                           ifelse(ssf_data@data$landcover == 8, "Shrub",                               
                              ifelse(ssf_data@data$landcover == 9, "Water", 
                                 ifelse(ssf_data@data$landcover == 10, "Rock-Ice", 
                                    ifelse(ssf_data@data$landcover == 11, "Cloud", 
                                         ifelse(ssf_data@data$landcover == 12, "Burn-Forest",                               
                                               ifelse(ssf_data@data$landcover == 13, "Burn-Grassland", 
                                                    ifelse(ssf_data@data$landcover == 14, "Burn-Shrub", 
                                                           ifelse(ssf_data@data$landcover == 15, "Alpine Herb", "Alpine Shrub"))))))))))))))))
table(ssf_data$landcover.cat, ssf_data$used)
head(ssf_data)

##################################################################################################################################

# create regular data frame with data
ssf_data <- data.frame(ssf_data@data)

# subset NAs out of data for model fitting and selection

# check for NA

length(which(is.na(ssf_data$meanNDVI)==T))
length(which(is.na(ssf_data$elevation)==T))
length(which(is.na(ssf_data$distha)==T))
length(which(is.na(ssf_data$disthha)==T))
length(which(is.na(ssf_data$landcover.cat)==T))
length(which(is.na(ssf_data$timeNDVI)==T))

# looks like just landcover is an issue
# subset out NA's for landcover (i.e., remove rows so data frames are consistent across models)

indexes.to.keep <- which(is.na(ssf_data$landcover.cat)==F)
ssf_data <- ssf_data[indexes.to.keep, ]

#alternatively we could have done
# ssf_data <- ssf_data[is.na(ssf_data$landcover.cat)==F,]

# check that NA's are removed from landcover
length(which(is.na(ssf_data$landcover.cat)==T))

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
################################################################################################################
#### Objective 4.0 FIT step selection function Clogit model

## first make sure we understand the stratum field
str(ssf_data)
head(ssf_data)
head(sort(ssf_data$stratum))
head(table(ssf_data$stratum, ssf_data$used))

### So, there are 5 matched case availabile locations for each 1 used location. It is this stratum that distinguishes clogit from logistic regression. 

## start with univariate models 
m.elev <- clogit(used ~ elevation + strata(stratum), data = ssf_data)
summary(m.elev)	   

#### Compare to naive logistic regression model ignoring structure
elev <- glm(used ~ elevation, data = ssf_data, family = binomial(logit))
summary(elev)	

### Compare clogit and naive logit for next variables. 
### Distance to human access
m.distha <- clogit(used ~ distha + strata(stratum), data = ssf_data)
summary(m.distha)

distha <- glm(used ~ distha, data = ssf_data, family = binomial(logit))
summary(distha)

#### Distance to high human access
m.disthha <- clogit(used ~ disthha + strata(stratum), data = ssf_data)
summary(m.disthha)

disthha <- glm(used ~ disthha, data = ssf_data, family = binomial(logit))
summary(disthha)

#### Landcover
m.landcov <- clogit(used ~ I(landcover.cat) + strata(stratum), data = ssf_data)
summary(m.landcov)

landcov  <- glm(used ~ I(landcover.cat), data = ssf_data, family = binomial(logit))
summary(landcov)

#### Discuss interpretation in class, but fair enough comparison is a bit lame because nothing is turning out significant. 

#### Compare Predictions from full models

m.full  <- clogit(used ~ elevation + distha + disthha + I(landcover.cat) +strata(stratum), data = ssf_data)
summary(m.full)

full  <- glm(used ~ elevation + distha + disthha + I(landcover.cat), data = ssf_data, family = binomial(logit))
summary(full)

#### You can now get predictions from a Clogit model, but the predictions mean something very different here
#### The Cox model is a relative risk model; predictions of type "linear predictor", "risk", and "terms" are all relative to the sample from which they came. By default, the reference value for each of these is the mean covariate within strata. The primary underlying reason is statistical: a Cox model only predicts relative risks between pairs of subjects within the same strata, and hence the addition of a constant to any covariate, either overall or only within a particular stratum, has no effect on the fitted results. Using the reference="strata" option causes this to be true for predictions as well.
ssf_data$m.full.pred <- predict(m.full, type="expected")
hist(ssf_data$m.full.pred)

ssf_data$full.pred <- predict(full, type="response")

plot(ssf_data$full.pred, ssf_data$m.full.pred)
abline(lm(ssf_data$full.pred~ ssf_data$m.full.pred))
cor(ssf_data$full.pred, ssf_data$m.full.pred)

## Discussion: what to make of the different interpretations of Clogit versus Logit for the same dataset?



#### MOdel Selection for cLogit Models

## fit null model
m.null  <- clogit(used ~ strata(stratum), data = ssf_data)
summary(m.null)

AIC(m.full, m.landcov, m.disthha, m.distha, m.elev, m.null)


stepAIC(m.full, direction = "backward")

## So - you cannot fit models with 0 parameters in Coxph model
AIC(m.full, m.landcov, m.disthha, m.distha, m.elev, full)
## note you get an error message, cannot compare to the Naive Logit model for 2 reasons
## 1) differnet # of data points (why?)
## 2) different calculation of the log-likelihood - confiditonal likelihood vs. maximum likelihood. apples and oranges. 
## Cannot solve whether or not to add conditioning based on AIC. 


####### Objective 5.0 Matched-case control over multiple individual - Mixed-effects Clogit

##### 1) Here are the first 2 papers that figured out how to add a random intercept for each individual animal (e.g.), however, it did so in MATLAB. So, its mostly inaccessible to biologists. 
### Craiu, R. V., T. Duchesne, D. Fortin, and S. Baillargeon. 2011. Conditional Logistic Regression With Longitudinal Follow-up and Individual-Level Random Coefficients: A Stable and Efficient Two-Step Estimation Method. Journal of Computational and Graphical Statistics 20:767-784.
### Duchesne, T., D. Fortin, and N. Courbin. 2010. Mixed conditional logistic regression for habitat selection studies. Journal of Animal Ecology 79:548-555.

### 2) However, there have been a few big breakthrough’s lately with the mclogit package http://cran.r-project.org/web/packages/mclogit/mclogit.pdf  or the coxme package here http://cran.r-project.org/web/packages/coxme/coxme.pdf I just played around with both of these packages and they are actually. 
### I figured out the mclogit package just now, and, for example, provide example code here. Imagine instead of using just elk29’s data in our analysis, we had 2 elk, identified by an elkid field. I created 1 at random in the wolf id field to keep track of the 2 wolf packs in the original mccwolf data. 

install.packages(c("coxme", "mclogit"))
library(coxme)
library(mclogit)

##### Bring in Bison dataset 

##### This data set was collected in order to study habitat selection by groups of free-ranging bison. For each observed group, 
#### two individuals (dyad) equipped with GPS radio-collars were followed simultaneously. A cluster is defined here as a pair of bison. 
#### This data set contains 20 clusters. The number of strata per cluster varies between 13 and 345 for a total of 1410 strata. 
### A stratum is composed of two visited GPS locations (one for each individual) gathered at the same time, together with 10 random locations 
### (five drawn within 700 m of each of the two focal bison). Therefore, there are 12 observations per stratum, with 2 cases (Y=1) and 10 controls 
### (Y=0). However, due to problems in the data collection, 17 of the 1410 strata have only 6 observations (1 case and 5 controls).

#install.packages("TwoStepCLogit")
#library(TwoStepCLogit)
head(bison)
str(bison)
head(table(bison$Strata,bison$Cluster))
hist(bison$biomass)
hist(bison$meadow)
boxplot(bison$biomass~ bison$meadow)

bison.mcclogit <- mclogit(cbind(Y, Strata) ~pmeadow + biomass, random ~1|Cluster, data=bison)
summary(bison.mcclogit)


bison.mcfixed <- mclogit(cbind(Y, Strata) ~pmeadow + biomass, data=bison)
summary(bison.mcfixed)

bison.naive <- glm(Y ~ pmeadow + biomass, data = bison, family = binomial(logit))
summary(bison.naive)

AIC(bison.mcclogit, bison.mcfixed)

AIC(bison.naive)
#### Note that the AIC's are not comparable!!

bison$mcclogit <- predict(bison.mcclogit, response="expected")
bison$mcfixed <- predict(bison.mcfixed, response="expected")
bison$naive <- predict(bison.naive, type="response")
str(bison)

plot(bison.mcclogit, bison.naive)
