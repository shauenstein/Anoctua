# run calibrationScript.r up to Line 174
set.seed(01022016)

# load packages
library(sp)
library(raster)
library(ranger)
library(LOSim)
library(abctools)
library(BayesianTools) 

# source functions
source("../code/CalibDataPrep.r")
source("../code/sumStats.r")
source("../code/getSumStats.r")
source("../code/createPriors.r")
source("../code/runABC.r")
source("../code/ABCsummaryChoice.r")
source("../code/plotSummariesCorrelation.r")
source("../code/ABCgetEstimate.r")
source("../code/mergeLOSim.r")

##############################
# load/prepare model input   #
#--------##--------##--------##--------##--------##--------##--------##--------#
# load telemetry data
cleaned <- read.csv("../data/natal_dispersal_cleaned.csv",
                    stringsAsFactors = FALSE)
# date as POSIXct
cleaned$date <- as.POSIXct(cleaned$date)
cleaned$date <- lubridate::with_tz(cleaned$date, "GMT")
# load IDs of training and validation set
load("../data/trainTestIDs.RData")

indiv = 1
#----------------------------------------------------#

# classify dispersal locations # date of first outing (400m buffer)
firstOut <- 1:3
for(i in seq_along(unique(cleaned$id))){
  idIndex <- which(cleaned$id == unique(cleaned$id)[i])
  birthIndex <- idIndex[which.min(cleaned$date[cleaned$id == 
                                                 unique(cleaned$id)[i]])]
  distToBL <- sapply(idIndex, function(x) 
    sqrt((cleaned$posX[x] - cleaned$posX[birthIndex])^2
         + (cleaned$posY[x] - cleaned$posY[birthIndex])^2))
  
  search = TRUE
  j = 1
  while(search){
    if(distToBL[j] > 400){
      search = FALSE
      firstOut[i] = idIndex[j]
    }
    j = j + 1
  }
} 

firstOut <- cleaned[firstOut,]


# select locations within dispersal period, i.e. 
dispersal <- cleaned[-c(1:NROW(cleaned)),]
for(x in unique(cleaned$id)){
  subs <-  cleaned[cleaned$id == x 
                   & cleaned$date >= firstOut$date[firstOut$id == x] 
                   & cleaned$date <= (firstOut$date[firstOut$id == x] + 3600*24*60),]
  dispersal <- rbind(dispersal, subs)
}
rm(search, x, j, i, idIndex, distToBL, birthIndex, subs)

# # spatialpointsdataframe
dispSpdf <- SpatialPointsDataFrame(coords = dispersal[,c("posX","posY")], 
                                   data = dispersal) 
# add correct projection: DHDN / Gauss-Kruger zone 3
proj4string(dispSpdf) <- CRS("+init=epsg:31467")

#----------------------------------------------------#
# get start points from dispersal data, i.e. first recorded locations in dispersal period
startIndex <- 1:3
for(i in seq_along(unique(dispSpdf$id))){
  idIndex <- which(dispSpdf$id == unique(dispSpdf$id)[i])
  startIndex[i] <- idIndex[which(dispSpdf$date[idIndex] == min(dispSpdf$date[idIndex]))]
} 
startPoints <- dispSpdf[startIndex, ]

# telemetry data extent + 50 km buffer
calibrationExtent <- extent(startPoints) + rep(c(-30000,30000), 2)

#----------------------------------------------------#
# load habitat suitability
habitatSuitability <- raster("../envVars/envOnStudyExtent/habitatSuitabilityScaled.tif")
# crop habitat suitability map
habSuitCalib <- crop(habitatSuitability, calibrationExtent)
habSuitCalib[is.na(values(habSuitCalib))] <- 0 # set NAs to zero <- reflection
HSvals <- as.matrix(habSuitCalib)

#----------------------------------------------------#
# priors for observation error
Accuracy <- sapply(unique(dispSpdf$id), function(x) dispersal$accuracy[dispersal$id == x])
accLookUp <- c(0,5,20,50,100,200)
AccMetric <- sapply(unique(dispSpdf$id), function(x) accLookUp[Accuracy[[x]]+1])
AccMetric <- unlist(AccMetric)

#----------------------------------------------------#
# prepare telemetry data
calibrationData <- calibTransform(x = dispersal$posX, y = dispersal$posY,
                                  date = dispersal$date, id = dispersal$id,
                                  habitatSuitability = habSuitCalib)

# filter individuals with more than 50 recorded locations
nRelocs <- sapply(seq_along(calibrationData), function(x) NROW(calibrationData[[x]]$data))
calibID <- names(calibrationData)[nRelocs  > 50 & names(calibrationData) %in% trainingID]
validID <- names(calibrationData)[nRelocs  > 50 & names(calibrationData) %in% testID]

#save(calibID,validID, file = "../data/trainTestIDs_g50rl.RData")

####################
# parameter sample #
#--------##--------##--------##--------##--------##--------##--------##--------#
# number of simulations per individual
n = 1

# parameter names
pnames <- c("habitatPreference", 
            "stepShape", 
            "stepScale", 
            "directionalBias", 
            "dispRestMean", 
            "dispRestSd", 
            "maximumEffort", 
            "observationError")

# prior distributions
distr <- list(runif, runif, NULL, runif, NULL, NULL, runif, rnorm)

# prior ranges
ranges <- list(c(0.1, 5), 
               c(0.1, 4), 
               2, 
               c(0.1, 3.5), 
               1.5,  
               0.7, 
               c(500, 8000),
               c(mean(AccMetric),sd(AccMetric))
)

# which parameter are fixed
fixed = c(3,5,6)

# create parameters
Ps <- createPriors(pnames, distr, ranges, n, redraw = 8, fixed)

# careful: throw away first x in y. run, throw away first 30 000 (r.seed issue)
# Ps <- Ps[100001:180000,]
####################
# model setup      #
#--------##--------##--------##--------##--------##--------##--------##--------#

# environment input
HS <- list(values = HSvals, upperleft = calibrationExtent[c(1,4)], resolution = 20)


# runABC(Ps, HS, calibrationData[calibID[1]], parallel = 13, 
#        maxIterations = 40000, randomSeed = NULL, 
#        sumStatsSampleSize = 5000, sumStats, getSumStats)

# initiate runABC arguments
parameters = Ps 
habitatSuitability = HS 
observations = calibrationData[calibID[indiv]]
parallel = FALSE 
maxIterations = 40000 
randomSeed = NULL
sumStatsSampleSize = 5000

# number of parameter combinations
nParameters = NROW(parameters)

# prepare parameters
parMat <- cbind(rep(observations[[1]]$data[1 , "xObserved"], nParameters),
                rep(observations[[1]]$data[1 , "yObserved"], nParameters),
                parameters
)

# run simulation
simData <- LOSim::runSimulation(habitatSuitability$values, 
                                maxIterations, 
                                parMat, 
                                randomSeed,
                                "observed", 
                                observations[[1]]$data[, "timestamp"], 
                                NULL,
                                habitatSuitability$upperleft,
                                habitatSuitability$resolution,
                                observations[[1]]$startTime,
                                observations[[1]]$startDoy,
                                observations[[1]]$centroid[2],
                                observations[[1]]$centroid[1])

sims <- simData[[1]][[1]]


# plot timestamps, observed and simulated
x11(15,5)
plot(jitter(sims[,1], 100), rep(0.5, nrow(sims)), type = "b", pch = "|")
points(observations[[1]]$data[, "timestamp"], rep(.6,nrow(observations[[1]]$data)), col = 2, pch = "|")

# plot paths
x11(10,10)
# observed
plot(observations[[1]]$data[, c("xObserved", "yObserved")], type = "b", 
     col= paste0("grey",ceiling(seq(1,80, len = nrow(observations[[1]]$data)))),
     xlim=range(observations[[1]]$data[,"xObserved"])+c(-5000,5000),
     ylim=range(observations[[1]]$data[,"yObserved"])+c(-5000,5000))

for(i in 1:n){
  points(simData$simulation[[i]][, c("xObserved", "yObserved")], type = "b", col = i +1)
}

# for all individuals
plot(calibrationData[[1]]$data[, c("xObserved", "yObserved")], type = "b", 
     xlim=calibrationExtent[1:2],
     ylim=calibrationExtent[3:4])

for(i in seq_along(calibrationData)[-1]){
  points(calibrationData[[i]]$data[, c("xObserved", "yObserved")], type = "b")
}

## check processes
library(plot3D)

# test binomial model for roost selection
exp(- 1 / 0.98 * 0.7)
curve(plogis(exp(- 1 / x * 1))-0.5, xlim = c(0,1), ylim = c(0,1), xlab = "energy", ylab = "prob of roosting")
curve(plogis(exp(- 1 / x * 0.1))-0.5, add = T)


# energy has negative effect, HS positive effect on prob of roosting
# HS = 0
intercept = 0; a = -5;b = 0; c = 5; d = 0
# curve(plogis(0 - a * x + b * HS + c * HS^2 + d*x^2), xlim = c(0,1), ylim = c(0,1), xlab = "energy", ylab = "prob of roosting")
# for(hs in seq(0.1,1,0.1)) curve(plogis(0 - a * x + b * hs  + c * hs^2 + d*x^2), add = T)
# abline(h = 0.5)
# 3D plot
hs = en = seq(0,1,length.out = 300)
M <- mesh(hs, en)
x <- M$x; y <- M$y
surf3D(x,y, z = plogis(intercept + a * x + b * y  + c * y^2 + d*x^2), bty = "b2", xlab = "energy", ylab = "habitat suit.", zlab = "probability of roosting", theta = 50, phi = 10, zlim=c(0,1))

# plogis = 1 / (1 + exp(-(x-m)/s)); with m = 0, s = 1
curve(1/(1+exp(-(intercept + a * x + b * HS))), xlim = c(0,1), ylim = c(0,1), xlab = "energy", ylab = "prob of roosting")
for(hs in seq(0.1,1,0.1)) curve(1/(1+exp(-(intercept + a * x + b * hs))), add = T)


# energy uptake
# energy += ((1/(1 + exp(- pow(cumEnv/20, recoveryRate) * ((timestamp.at(iter + 1) - timestamp.at(iter)))))) 
#   - 0.5) / 0.5
rR = 2.5
HS = 0
curve(((1/(1 + exp(- (HS/20)^rR * ((x)))))  - 0.5) / 0.5, xlim = c(0,3*24*60), ylim = c(0,1), xlab = "energy", ylab = "prob of roosting")
for(hs in seq(0.1,1,0.1)) curve(((1/(1 + exp(- (hs/20)^rR * ((x)))))  - 0.5) / 0.5, add = T)


# roosting time
#timestamp.at(iter) + exp(- energy / pow(cumEnv*2, 2)) * maxRoostTime
energy = 0.1
HS = 0.5
curve(exp(- energy / (2*HS)^2)*x, ylim = c(0,3*24*60), xlim = c(0,3*24*60), xlab = "maximum roosting time", ylab = "roosting time")


