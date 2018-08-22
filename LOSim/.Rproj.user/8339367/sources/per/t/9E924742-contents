################################################################################
#SH25-09-2016                                                                  #
# re-fit parameters                                                            #
#                                                                              #
################################################################################

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##
#set seed
set.seed(25012017)

# load packages
library(LOSim)
library(sp)
library(raster)
library(ranger)
library(foreach)
library(doSNOW)
library(ggplot2)

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

# load rcpp function to calculate distance between observed and simulated data
Rcpp::sourceCpp("../code/c++/distanceObservedSimulated.cpp")

# individual to use
gp <- 69

# if predictions should consecutively be included as predictors
predictionPredictor <- TRUE
parallel <- FALSE
regr.adj <- TRUE
comp.MAP <- TRUE
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##
# in sample refitting                                                         ##
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##

# load simulations + remove all but for the first individual
load("../data/summaries.RData")
sumOut <- summaries; rm(summaries)
# upper and lower bounds for truncated multivariate normal distribution
lbounds <- c(habitatPreference = 0.1, stepShape = 0.1, directionalBias = 0.1, roostLambda = 0.01)
ubounds <- c(habitatPreference = 5, stepShape = 4, directionalBias = 3.5, roostLambda = 30)
sumOut$summarySelection$simulation <- sumOut$summarySelection$simulation[[gp]]
sumOut$summarySelection$observed <- sumOut$validity <- sumOut$observations <- NULL
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##

# sample parameters from prior sample
n = 1000
cutoff <- 0
Pwhich <- sample(which(rowSums(sapply(seq_along(sumOut$summarySelection$targetParameters),
                                      function(x) sumOut$parameters[,sumOut$summarySelection$targetParameters[x]] > 
                                        (lbounds[x] + (ubounds[x]-lbounds[x])*cutoff) &
                                        sumOut$parameters[,sumOut$summarySelection$targetParameters[x]] < 
                                        (ubounds[x] - (ubounds[x]-lbounds[x])*cutoff))) == 4),
                 size = n)
# Pwhich <- sample(NROW(sumOut$parameters), size = n)
PReFit <- sumOut$parameters[Pwhich,]

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##

# rejection sampling
proportionFiltered <- 1/500
numberParametersFiltered <- ceiling(nrow(sumOut$parameters) * proportionFiltered)

distances <- filtered <- pMedian <- pQuantiles <- pMEst <- vector("list", n)   
if(regr.adj) pMedian.regrAdj <- pQuantiles.regrAdj <- filtered.regrAdj <- vector("list", n)

for(i in seq(n)){
  
  # standardise by range
  reference <- unlist(diff(apply(sumOut$summarySelection$simulation[-Pwhich[i],], 2, range)))
  
  distances[[i]] <- dObsSim(simulated = as.matrix(sumOut$summarySelection$simulation[-Pwhich[i],]),
                            observed = unlist(sumOut$summarySelection$simulation[Pwhich[i],]),
                            reference = reference)
  
  # accepted indices
  filtered_indices <- sort(distances[[i]], method = "quick",
                           index.return = TRUE)$ix[1:numberParametersFiltered]
  
  filtered[[i]] <- sumOut$parameters[-Pwhich[i],][filtered_indices,
                                                  sumOut$summarySelection$targetParameters]
  
  if(regr.adj){
    # copy filtered
    filtered.regrAdj[[i]] <- filtered[[i]]
    # S_sim of accepted sample
    filtered_ss <- sumOut$summarySelection$simulation[-Pwhich[i],][filtered_indices,]
    for(c in seq(NCOL(filtered_ss))) # center around S_obs
      filtered_ss[,c] <- filtered_ss[,c] - sumOut$summarySelection$simulation[Pwhich[i],c]
    # post-sampling regression adjustment, multivariate linear model
    psr.fm <- lm(as.formula(paste("filtered.regrAdj[[i]]~",paste(colnames(filtered_ss), collapse = "+"))), 
                 data = filtered_ss)$coefficients[-1,] # remove intercept
    # correct for epsilon != 0
    filtered.regrAdj[[i]] <- sapply(seq(NCOL(filtered_ss)), 
                            function(x) filtered.regrAdj[[i]][,x] - rowSums(t(psr.fm[,x] * t(filtered_ss))))
    # restrict corrections to prior bounds
    filtered.regrAdj[[i]] <- sapply(seq(NCOL(filtered_ss)),
                            function(x) ifelse(filtered.regrAdj[[i]][,x] < lbounds[x], lbounds[x], 
                                               ifelse(filtered.regrAdj[[i]][,x] > ubounds[x], ubounds[x], 
                                                      filtered.regrAdj[[i]][,x])))
    colnames(filtered.regrAdj[[i]]) <- colnames(sumOut$parameters)[sumOut$summarySelection$targetParameters]
    
    pMedian.regrAdj[[i]] <- apply(filtered.regrAdj[[i]], 2, median)
    pQuantiles.regrAdj[[i]] <- apply(filtered.regrAdj[[i]], 2, quantile, probs = c(0.025,0.975))
  }
  
  pMedian[[i]] <- apply(filtered[[i]], 2, median)
  pQuantiles[[i]] <- apply(filtered[[i]], 2, quantile, probs = c(0.025,0.975))
  pMEst[[i]] <- unlist(sumOut$summarySelection$simulation[Pwhich[i],])
}

if(comp.MAP){
  if(!parallel){
    library(tmvtnorm)
    pMAP <- vector("list", n)
    if(regr.adj) pMAP.regrAdj <- vector("list", n)
    # progress bar
    pb <- txtProgressBar(min = 1, max = n, style = 3)
    for(i in seq.int(n)){
      pMAP[[i]] <- mle.tmvnorm(filtered[[i]], method = "L-BFGS-B", 
                               lower.bounds = lbounds, 
                               upper.bounds = ubounds)@coef[1:NCOL(filtered[[i]])]
      pMAP[[i]] <- ifelse(pMAP[[i]] < lbounds, lbounds,
                          ifelse(pMAP[[i]] > ubounds, ubounds,
                                 pMAP[[i]]))
      names(pMAP[[i]]) <- colnames(sumOut$parameters)[sumOut$summarySelection$targetParameters]
      if(regr.adj){
        pMAP.regrAdj[[i]] <- mle.tmvnorm(filtered.regrAdj[[i]], method = "L-BFGS-B", 
                                 lower.bounds = lbounds, 
                                 upper.bounds = ubounds)@coef[1:NCOL(filtered.regrAdj[[i]])]
        pMAP.regrAdj[[i]] <- ifelse(pMAP.regrAdj[[i]] < lbounds, lbounds,
                            ifelse(pMAP.regrAdj[[i]] > ubounds, ubounds,
                                   pMAP.regrAdj[[i]]))
        names(pMAP.regrAdj[[i]]) <- colnames(sumOut$parameters)[sumOut$summarySelection$targetParameters]
      }
      # progress bar
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }else{
    # parallel execution
    library(foreach)
    # library(doParallel)
    library(doSNOW)
    
    if(parallel == T | parallel == "auto"){
      cores <- parallel::detectCores() - 1
      message("parallel, set cores automatically to ", cores)
    }else if (is.numeric(parallel)){
      cores <- parallel
      message("parallel, set number of cores manually to ", cores)
    }else stop("wrong argument to parallel")
    
    
    cl <- parallel::makeCluster(cores)
    doSNOW::registerDoSNOW(cl)
    
    # set up progress bar
    pb <- txtProgressBar(min = 1, max = n, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    pMAP <- foreach::foreach(i = seq.int(n), 
                             .options.snow=opts, .packages = "tmvtnorm") %dopar% {
                               maps <- mle.tmvnorm(filtered[[i]], method = "L-BFGS-B",
                                                   lower.bounds = lbounds, 
                                                   upper.bounds = ubounds)@coef[1:NCOL(filtered[[i]])]
                               maps <- ifelse(maps < lbounds, lbounds,
                                              ifelse(maps > ubounds, ubounds,
                                                     maps))
                               names(maps) <- colnames(sumOut$parameters)[sumOut$summarySelection$targetParameters]
                               maps
                             }
    if(regr.adj){
      pMAP.regrAdj <- foreach::foreach(i = seq.int(n), 
                               .options.snow=opts, .packages = "tmvtnorm") %dopar% {
                                 maps <- mle.tmvnorm(filtered.regrAdj[[i]], method = "L-BFGS-B",
                                                     lower.bounds = lbounds, 
                                                     upper.bounds = ubounds)@coef[1:NCOL(filtered.regrAdj[[i]])]
                                 maps <- ifelse(maps < lbounds, lbounds,
                                                ifelse(maps > ubounds, ubounds,
                                                       maps))
                                 names(maps) <- colnames(sumOut$parameters)[sumOut$summarySelection$targetParameters]
                                 maps
                               }
    }
    
    # free memory from workers
    close(pb)
    stopCluster(cl)
  }
}


posteriorInS = NULL
posteriorInS$parameters = filtered
posteriorInS$median = pMedian
posteriorInS$quantiles = pQuantiles
if(comp.MAP) posteriorInS$MAP = pMAP
if(regr.adj){
  posteriorInS$parameters.regrAdj = filtered.regrAdj
  posteriorInS$median.regrAdj = pMedian.regrAdj
  posteriorInS$quantiles.regrAdj = pQuantiles.regrAdj
  if(comp.MAP) posteriorInS$MAP.regrAdj = pMAP.regrAdj
} 
posteriorInS$model.estimate = pMEst
posteriorInS$targetParameters = sumOut$summarySelection$targetParameters
posteriorInS$proportion = proportionFiltered


# save(posteriorInS, file = "../data/posteriorInS_ParameterReFit.RData")
# load("../data/posteriorInS_ParameterReFit.RData")

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##
# out of sample refitting                                                     ##
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##

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
calibrationExtent <- extent(startPoints) + rep(c(-50000,50000), 2)

#----------------------------------------------------#
# load habitat suitability
habitatSuitability <- raster("../envVars/envOnStudyExtent/habitatSuitabilityScaled.tif")
# crop habitat suitability map
habSuitCalib <- crop(habitatSuitability, calibrationExtent)
habSuitCalib[is.na(values(habSuitCalib))] <- 0 # set NAs to zero <- reflection
HSvals <- as.matrix(habSuitCalib)
# update calibrationExtent
calibrationExtent <- extent(habSuitCalib)

#----------------------------------------------------#
# priors for observation error
Accuracy <- sapply(unique(dispSpdf$id), function(x) dispersal$accuracy[dispersal$id == x])
accLookUp <- c(0,5,20,50,100,200)
AccMetric <- sapply(unique(dispSpdf$id), function(x) accLookUp[Accuracy[[x]]+1])
AccMetric <- unlist(AccMetric)

#----------------------------------------------------#
# prepare telemetry data
rscRange = 200
calibrationData <- calibTransform(x = dispersal$posX, y = dispersal$posY,
                                  date = dispersal$date, id = dispersal$id,
                                  habitatSuitability = habSuitCalib, rscExtent = rscRange)

# filter individuals with more than 50 recorded locations
nRelocs <- sapply(seq_along(calibrationData), function(x) NROW(calibrationData[[x]]$data))
calibID <- names(calibrationData)[nRelocs  > 50 & names(calibrationData) %in% trainingID]

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##
# simulate observed data for gp
# use parameters from above: PReFit


# create parameters
PsObs <- cbind(calibrationData[[calibID[gp]]]$data[1,"xObserved"],
               calibrationData[[calibID[gp]]]$data[1,"yObserved"],
               PReFit)

# simulate observations
simObs <- runSimulation(environment = HSvals, 
                        iterations = 40000, 
                        PsObs, 
                        randomSeed = NULL,
                        "observed", 
                        TimeAtObservation = calibrationData[[calibID[gp]]]$data[, "timestamp"],
                        maxDuration = NULL,
                        upperleft = calibrationExtent[c(1,4)],
                        environmentResolution = 20,
                        calibrationData[[calibID[gp]]]$startTime,
                        calibrationData[[calibID[gp]]]$startDoy,
                        calibrationData[[calibID[gp]]]$centroid[2],
                        calibrationData[[calibID[gp]]]$centroid[1])

names(simObs[[1]]) <- paste0("reFit_",seq.int(length(simObs[[1]])))

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##
# load simulations
load("../data/simulation.RData")

# reorder simulation output
simOut <- simulation; rm(simulation)
simOut$summaries <- list(simOut$summaries[[gp]]) # here 1 = gp
simOut$observations <- simObs[[1]]

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##

#  manual summary selection, take chuncks from LOSim::abcCreateSummaries

# number of available cores for parallelisation of ranger random forest
nCPU = 4

# compute summary statistics for simulated observations
ssSimObs <- getSumStats(simulations = simOut$observations, 
                        sampleSize = 5000, 
                        sumStats = sumStats)

# don't select any summary statistics + choose target parameters
summarySelection = 1:NCOL(simOut$summaries[[1]])
targetParameters = c(1,2,4,8)

tempSum = simOut$summaries[[1]]
tempSum[is.na( tempSum)] = 0

simObsPreds <- vector("list", n)
rfsummaries <- vector("list", length(targetParameters))
simRFPreds <- as.data.frame(matrix(0, ncol = length(targetParameters), nrow = NROW(tempSum)))
colnames(simRFPreds) <- paste0("S", seq_along(targetParameters))
simObsRFPreds <- as.data.frame(matrix(0, ncol = length(targetParameters), nrow = n))
colnames(simObsRFPreds) <- paste0("S", seq_along(targetParameters))

index = 1
for(j in targetParameters){
  if(predictionPredictor){
    TrainData <- data.frame(y = simOut$parameters[,j], 
                            tempSum[,summarySelection], 
                            simRFPreds)
  }else{
    TrainData <- data.frame(y = simOut$parameters[,j], 
                            tempSum[,summarySelection])
  }
  
  rfsummaries[[index]] <- ranger::ranger(y ~ ., data = TrainData, 
                                         num.trees = 500, write.forest = TRUE, 
                                         num.threads = nCPU, verbose = FALSE,
                                         importance =  "impurity")
  
  simRFPreds[,index] <- predict(rfsummaries[[index]], data = TrainData, num.threads = nCPU)$predictions
  
  # predict for simulated observations
  for(t in seq(NROW(ssSimObs[[1]]))){
    if(predictionPredictor){
      PredSimObs <- as.data.frame(matrix(c(ssSimObs[[1]][t, summarySelection], 
                                           simObsRFPreds[t,]), nrow = 1))
    }else{
      PredSimObs <- as.data.frame(matrix(ssSimObs[[1]][t, summarySelection], nrow = 1))
    }
    
    colnames(PredSimObs) <- colnames(TrainData)[-1]
    simObsRFPreds[t,index] <- predict(rfsummaries[[index]], data = PredSimObs, num.threads = nCPU)$predictions
  }
  cat("targetParameter",targetParameters[index], "\n")
  index <- index + 1
}

# re-organise output
for(t in seq(NROW(ssSimObs[[1]]))){
  simObsPreds[[t]] <- simObsRFPreds[t,]
}


# collect results
# simOut$summarySelection$predict = predict
simOut$summarySelection$predictionPredictor = predictionPredictor
simOut$summarySelection$simulation = simRFPreds
simOut$summarySelection$observed = simObsPreds
simOut$summarySelection$method = "rF"
simOut$summarySelection$targetParameters = targetParameters
simOut$summarySelection$summarySelection = summarySelection

# save(simOut, file = "../data/simObs_ParameterReFit.RData")
# load("../data/simObs_ParameterReFit.RData")


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##
# copy simulations 
simOut$summarySelection$simulation <- lapply(seq(n), function(x) simOut$summarySelection$simulation)
# rejection
posteriorOutS <- getEstimate(simOut, proportionFiltered = proportionFiltered, 
                             parallel = FALSE, regr.adj = regr.adj, comp.MAP = comp.MAP)

# save(posteriorOutS, file = "../data/posteriorOutS_ParameterReFit.RData")
# load("../data/posteriorOutS_ParameterReFit.RData")

# save all in one file
save(posteriorInS,
     posteriorOutS,
     PReFit,
     Pwhich, file = "../data/parameterReFit.RData")
# load("../data/parameterReFit.RData")
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##

pNames <- c("habitatPreference", "stepShape", "directionalBias", "roostLambda")
for(pInt in pNames){
  ggdf <- data.frame(truep = PReFit[,pInt],
                     estMedian = do.call(rbind, posteriorOutS$median.regrAdj)[,pInt],
                     upper = sapply(seq_along(posteriorOutS$quantiles.regrAdj), function(x) posteriorOutS$quantiles.regrAdj[[x]][2,pInt]),
                     lower = sapply(seq_along(posteriorOutS$quantiles.regrAdj), function(x) posteriorOutS$quantiles.regrAdj[[x]][1,pInt]))
  
  ggplot(ggdf, aes(truep, estMedian)) +
    geom_linerange(aes(ymin = lower, ymax = upper), col = "grey70") + 
    geom_point(size = 1.4, shape = 15) +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    coord_fixed() +
    xlim(lbounds[pInt], ubounds[pInt]) +
    ylim(lbounds[pInt], ubounds[pInt]) +
    xlab("") + 
    ylab("") +
    theme(
      panel.border = element_rect(fill = NA),
      axis.title = element_text(family = "Linux Libertine", size = 60),
      axis.text = element_text(family = "Linux Libertine", size = 60),
      legend.background = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  ggsave(paste0("refit_rf_",pInt,".pdf"),
         device = cairo_pdf,
         path = "../figures/results/refit_indiv_plots",
         width = 25, height = 25, units = "cm", dpi = 400)
}


# plot for in sample
# real vs. fitted

# plotParameterReFit <- function(posteriorObj, parameterNames, truePar, regr.adj = TRUE,
#                                plotMedian = FALSE, plotMAP = TRUE, plotMEst = FALSE, idLine = TRUE, ...){
#   # edit par
#   par(...)
#   
#   for(p in parameterNames){
#     plot(1:3, type = "n", xlim = range(truePar[,p]), 
#          ylim = range(truePar[,p]), ylab = "", xlab = "")
#     # 95% CI as line
#     if(regr.adj){
#       for(i in 1:NROW(truePar)) lines(rep(truePar[i,p],2), 
#                                       posteriorObj$quantiles.regrAdj[[i]][,p], col = rgb(.5,.5,.5,.5))
#       if(plotMedian){
#         # median as point
#         for(i in 1:NROW(truePar)) points(truePar[i,p], 
#                                          posteriorObj$median.regrAdj[[i]][p], cex = 0.5,pch = 16) 
#       }
#       if(plotMAP){
#         # median as point
#         for(i in 1:NROW(truePar)) points(truePar[i,p], 
#                                          posteriorObj$MAP.regrAdj[[i]][p], cex = 0.5,pch = 16) 
#       }
#     }else{
#       for(i in 1:NROW(truePar)) lines(rep(truePar[i,p],2), 
#                                       posteriorObj$quantiles[[i]][,p], col = rgb(.5,.5,.5,.5))
#       if(plotMedian){
#         # median as point
#         for(i in 1:NROW(truePar)) points(truePar[i,p], 
#                                          posteriorObj$median[[i]][p], cex = 0.5,pch = 16) 
#       }
#       if(plotMAP){
#         # median as point
#         for(i in 1:NROW(truePar)) points(truePar[i,p], 
#                                          posteriorObj$MAP[[i]][p], cex = 0.5,pch = 16) 
#       }
#     }
#     if(plotMEst){
#       # median as point
#       
#       for(i in 1:NROW(truePar)) points(truePar[i,p], 
#                                        posteriorObj$model.estimate[[i]][parameterNames == p], cex = 0.5,pch = 16) 
#     }
#     if(idLine) abline(0,1, col = 2, xpd = FALSE)
#   }
# }
# 
# pNames <- c("habitatPreference", "stepShape", "directionalBias", "roostLambda")
# pdf("../figures/results/refit_rf_inSample.pdf", width = 10, height = 10)
# plotParameterReFit(posteriorInS, pNames, PReFit, regr.adj = F, plotMedian = T, plotMAP = F, plotMEst = F, idLine = TRUE, 
#                    mfrow = c(2,2), pty = "s", mar = c(2,1,1,1), oma = c(1,1.5,0,0), las = 1, 
#                    mgp = c(1.7,0.5,0), tcl = .2, cex.axis = 2, cex.lab = 2, xpd = NA)
# dev.off()
# 
# pdf("../figures/results/refit_rf_outSample.pdf", width = 10, height = 10)
# plotParameterReFit(posteriorOutS, pNames, PReFit, regr.adj = F, plotMedian = T, plotMAP = FALSE, plotMEst = F, idLine = TRUE, 
#                    mfrow = c(2,2), pty = "s", mar = c(2,2,1,1), oma = c(1,1.5,0,0), las = 1, 
#                    mgp = c(1.7,0.8,0), tcl = .2, cex.axis = 2.3, cex.lab = 2, xpd = NA)
# dev.off()
# 


mestIs <- do.call(rbind, posteriorInS$model.estimate)
mestOs <- do.call(rbind, posteriorOutS$model.estimate)
medIs <- do.call(rbind, posteriorInS$median)
medOs <- do.call(rbind, posteriorOutS$median)

cor(medIs, mestIs)
cor(medOs, mestOs)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# plot for out of sample
# real vs. fitted
# 
# x11(10,10)
# par(mfrow=c(2,2), pty = "s")
# # habitat preference
# plot(PReFit[,"habitatPreference"], sapply(1:500, function(x) posteriorOutS$median[[x]]["habitatPreference"]))
# # habitat preference
# plot(PReFit[,"stepShape"], sapply(1:500, function(x) posteriorOutS$median[[x]]["stepShape"]))
# # habitat preference
# plot(PReFit[,"directionalBias"], sapply(1:500, function(x) posteriorOutS$median[[x]]["directionalBias"]))
# # habitat preference
# plot(PReFit[,"roostLambda"], sapply(1:500, function(x) posteriorOutS$median[[x]]["roostLambda"]))



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##