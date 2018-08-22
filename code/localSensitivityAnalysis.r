# SH 23-02-2017
#
# Local Sensitivity Analysis
#
#-#-#-#-#-#-#-#-#-#-#

##########
# Setup  #
#--------##--------##--------##--------##--------##--------##--------##--------#
#set seed
set.seed(01022016) # makes problem at n = 1303

# load packages
library(sp)
library(raster)
library(ranger)
library(LOSim)
library(abctools)
library(BayesianTools) 
library(lubridate)


# source functions
source("../code/CalibDataPrep.r")
source("../code/sumStats.r")
source("../code/getSumStats.r")
source("../code/createPriors.r")
source("../code/runABC.r")


# parameter sample function
spar <- function(tar, tseq, fvals, nr, names){
  pm <- matrix(NA, nrow = nr, ncol = length(fvals))
  colnames(pm) <- names
  for(i in seq(length(fvals))){
    if(i == tar){
      pm[,i] <- tseq
    }else{
      pm[,i] <- rep(fvals[i], nr)
    } 
  }
  return(pm)
}

# calculate adj.r.squared
Simp <- function(tseq, Ssim, plot.it=FALSE, imp.level = 0.2){
  imp <- rep(NA, NCOL(Ssim))
  for(i in seq(NCOL(Ssim))){
    fm <- lm(Ssim[,i]~tseq)
    imp[i] <- ifelse(summary(fm)$adj.r.squared > imp.level, 1, 0)
    if(plot.it){
      plot(Ssim[,i]~tseq, xlab = "", ylab = "", main = paste0("S",i))
      abline(fm, col=2)
      mtext(paste("adj.R2 =",summary(fm)$adj.r.squared), 1, line = 3)
      readline(prompt="Press [enter] to continue")
    }
  }
  return(which(imp == 1))
}


# performs sensitivity analysis
SAS <- function(n, by.rep = 20, tar, t.range, fvals, pnames, HS, obs, 
                plotSAS = FALSE, imp.level = 0.5){
  
  tseq <- rep(seq(t.range[1],t.range[2], len = n/by.rep), each = by.rep)
  Ps <- spar(tar, tseq, fvals, n, pnames) 
  
  Ssim <- runABC(Ps, HS, obs, parallel = FALSE, 
                 maxIterations = 40000, randomSeed = NULL, 
                 sumStatsSampleSize = 5000, sumStats, getSumStats)[["summaries"]][[1]]
  
  return(Simp(tseq, Ssim, plot.it = plotSAS, imp.level = imp.level))
  
}

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
calibrationExtent <- extent(startPoints) + rep(c(-40000,40000), 2)

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

# choose observed data
obs <- calibrationData[1]



# environment input
HS <- list(values = HSvals, upperleft = calibrationExtent[c(1,4)], resolution = 20)

####################
# parameter sample #
#--------##--------##--------##--------##--------##--------##--------##--------#
# number of simulations
# parameter names
pnames <- c("habitatPreference", 
            "stepShape", 
            "stepScale", 
            "directionalBias", 
            "dispRestMean", 
            "dispRestSd", 
            "maximumEffort", 
            "roostLambda",
            "observationError",
            "rscRange")

#--------##--------##--------##--------##--------##--------##--------##--------#


t1 <- SAS(n=500, by.rep = 20, tar = 1, t.range = c(-10,4), 
          fvals = c(NA,2,2,1.5,1.5,0.7,4000,2,mean(AccMetric), 100),
          pnames = pnames, HS = HS, obs = obs)

t2 <- SAS(n=500, by.rep = 20, tar = 2, t.range = c(0.1,4), 
          fvals = c(0,NA,2,1.5,1.5,0.7,4000,2,mean(AccMetric), 100),
          pnames = pnames, HS = HS, obs = obs)

t4 <- SAS(n=500, by.rep = 20, tar = 4, t.range = c(0.1,4), 
          fvals = c(0,2,2,NA,1.5,0.7,4000, 2,mean(AccMetric), 100),
          pnames = pnames, HS = HS, obs = obs)

t8 <- SAS(n=500, by.rep = 20, tar = 8, t.range = c(-10,4), 
          fvals = c(0,2,2,1.5,1.5,0.7,4000,NA,mean(AccMetric), 100),
          pnames = pnames, HS = HS, obs = obs)

#--------##--------##--------##--------##--------##--------##--------##--------#


# 
n=200; 
by.rep = 25; 
tar = 1; 
t.range = c(-10,4); 
fvals = c(NA,2,2,1.5,1.5,0.7,1000,8,mean(AccMetric), 100)

tseq <- rep(seq(t.range[1],t.range[2], len = n/by.rep), each = by.rep)
Ps <- spar(tar, tseq, fvals, n, pnames) 

# number of parameter combinations
nParameters = NROW(Ps)

# prepare parameters
parMat <- cbind(rep(obs[[1]]$data[1 , "xObserved"], nParameters),
                rep(obs[[1]]$data[1 , "yObserved"], nParameters),
                Ps
)

# run simulation
simData <- LOSim::runSimulation(HS$values, 
                                40000, 
                                parMat, 
                                NULL,
                                "observed",
                                obs[[1]]$data[, "timestamp"], 
                                NULL,
                                HS$upperleft,
                                HS$resolution,
                                obs[[1]]$startTime,
                                obs[[1]]$startDoy,
                                obs[[1]]$centroid[1],
                                obs[[1]]$centroid[2]
)



x11()
plot(0, type = "n", xlim=parMat[1,1] + c(-20000,20000),
     ylim=parMat[1,2] + c(-20000,20000))
cols = rep(1:5, each = 4)
for(i in 1:20){
  lines(simData$simulation[[i]][,c("xObserved","yObserved")], col= cols[i])
}

  