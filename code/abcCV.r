################################################################################
#SH15-02-2017                                                                  #
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


abcCV <- function(summaryChoice, n = 500, method = "out", individual = NULL, predictionPredictor = TRUE, proportionFiltered = 1/200, nCPU = 16, sumStats){
  if(is.null(individual)) individual <- sample(length(summaryChoice$summarySelection$simulation), size = 1)
  if(individual > length(summaryChoice$summarySelection$simulation)) stop("The selected individual does not exist.")
  
  # reduce summaryChoice to chosen individual
  summaryChoice$summarySelection$simulation <- summaryChoice$summarySelection$simulation[[individual]]
  summaryChoice$summarySelection$observed <- summaryChoice$summarySelection$observed[[individual]]
  summaryChoice$observations <- list(summaryChoice$observations[[individual]])
  summaryChoice$summaries <- list(summaryChoice$summaries[[individual]])
  summaryChoice$validity <- list(summaryChoice$validity[[individual]])
  summaryChoice$nSimSteps  <- list(summaryChoice$nSimSteps[[individual]])
  
  
  # sample parameters from parameter priors
  Pwhich <- sample(NROW(summaryChoice$parameters), size = n)
  PReFit <- summaryChoice$parameters[Pwhich,]
  
  numberParametersFiltered <- ceiling(nrow(summaryChoice$parameters) * proportionFiltered)
  
  
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##
  # in sample refitting                                                         ##
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##
  if(method == "in"){
    distances <- filtered <- pMedian <- pQuantiles <- vector("list", n)   
    
    for(i in seq(n)){
      
      # standardise by range
      reference <- unlist(diff(apply(summaryChoice$summarySelection$simulation[-Pwhich[i],], 2, range)))
      
      distances[[i]] <- dObsSim(simulated = as.matrix(summaryChoice$summarySelection$simulation[-Pwhich[i],]),
                                observed = unlist(summaryChoice$summarySelection$simulation[Pwhich[i],]),
                                reference = reference)
      
      filtered[[i]] <- summaryChoice$parameters[-Pwhich[i],][sort(distances[[i]], method = "quick",
                                                           index.return = TRUE)$ix[1:numberParametersFiltered],
                                                           summaryChoice$summarySelection$targetParameters]
      
      pMedian[[i]] <- apply(filtered[[i]], 2, median)
      pQuantiles[[i]] <- apply(filtered[[i]], 2, quantile, probs = c(0.025,0.975))
      
      cat("replication ",i, "\n")
    }
    
    posteriorInS = NULL
    posteriorInS$parameters = filtered
    posteriorInS$median = pMedian
    posteriorInS$quantiles = pQuantiles
    posteriorInS$targetParameters = sumOut$summarySelection$targetParameters
    posteriorInS$proportion = proportionFiltered
  
  }else{
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##
    # out of sample refitting                                                     ##
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##
    
    # simulate new observations
    # create parameters
    PsObs <- cbind(summaryChoice$observations$data[1,"xObserved"],
                   summaryChoice$observations$data[1,"yObserved"],
                   PReFit)
    
    # simulate cv observations
    simObs <- runSimulation(environment = summaryChoice$environment$values, 
                            iterations = summaryChoice$settings$maxIterations, 
                            PsObs, 
                            randomSeed = summaryChoice$settings$randomSeed,
                            "observed", 
                            TimeAtObservation = summaryChoice$observations$data[, "timestamp"],
                            maxDuration = NULL,
                            upperleft = summaryChoice$environment$upperleft,
                            environmentResolution = summaryChoice$environment$resolution,
                            summaryChoice$observations$startTime,
                            summaryChoice$observations$startDoy,
                            summaryChoice$observations$centroid[2],
                            summaryChoice$observations$centroid[1])
    
    names(simObs[[1]]) <- paste0("reFit_",seq.int(length(simObs[[1]])))
    
    # overwrite observations
    summaryChoice$observations <- simObs[[1]]
    
    # compute summary statistics for simulated observations
    ssSimObs <- getSumStats(simulations = summaryChoice$observations, 
                            sampleSize = summaryChoice$settings$sumStatSampleSize, 
                            sumStats = sumStats)
    
    # don't select any summary statistics + choose target parameters
    summarySelection = summaryChoice$summarySelection$summarySelection
    targetParameters = summaryChoice$summarySelection$targetParameters
    
    tempSum = summaryChoice$summaries[[1]]
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
        TrainData <- data.frame(y = summaryChoice$parameters[,j], 
                                tempSum[,summarySelection], 
                                simRFPreds)
      }else{
        TrainData <- data.frame(y = summaryChoice$parameters[,j], 
                                tempSum[,summarySelection])
      }
      
      rfsummaries[[index]] <- ranger::ranger(y ~ ., data = TrainData, 
                                             num.trees = 500, write.forest = TRUE, 
                                             num.threads = nCPU, verbose = FALSE,
                                             importance =  "impurity")
      
      simRFPreds[,index] <- predict(rfsummaries[[index]], data = TrainData)$predictions
      
      # predict for simulated observations
      for(t in seq(NROW(ssSimObs[[1]]))){
        if(predictionPredictor){
          PredSimObs <- as.data.frame(matrix(c(ssSimObs[[1]][t, summarySelection], 
                                               simObsRFPreds[t,]), nrow = 1))
        }else{
          PredSimObs <- as.data.frame(matrix(ssSimObs[[1]][t, summarySelection], nrow = 1))
        }
        
        colnames(PredSimObs) <- colnames(TrainData)[-1]
        simObsRFPreds[t,index] <- predict(rfsummaries[[index]], data = PredSimObs)$predictions
      }
      cat("targetParameter",targetParameters[index], "\n")
      index <- index + 1
    }
    
    # re-organise output
    for(t in seq(NROW(ssSimObs[[1]]))){
      simObsPreds[[t]] <- simObsRFPreds[t,]
    }
    
    # collect results
    summaryChoice$summarySelection$simulation = simRFPreds
    summaryChoice$summarySelection$observed = simObsPreds
    
    # copy simulations 
    summaryChoice$summarySelection$simulation <- lapply(seq(n), function(x) summaryChoice$summarySelection$simulation)
  }
  
  return(summaryChoice)
}








# save(simOut, file = "../data/simObs_ParameterReFit.RData")
# load("../data/simObs_ParameterReFit.RData")


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##

# rejection
posteriorOutS <- getEstimate(simOut, proportionFiltered = proportionFiltered, parallel = FALSE)

# save(posteriorOutS, file = "../data/posteriorOutS_ParameterReFit.RData")
# load("../data/posteriorOutS_ParameterReFit.RData")

# save all in one file
# save(posteriorInS, 
#      posteriorOutS, 
#      PReFit,
#      Pwhich, file = "../data/parameterReFit.RData")

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##

# plot for in sample
# real vs. fitted

plotParameterReFit <- function(posteriorObj, parameterNames, truePar, plotMedian = TRUE, ...){
  # edit par
  par(...)
  
  for(p in parameterNames){
    plot(1:3, type = "n", xlim = range(truePar[,p]), 
         ylim = range(truePar[,p]))
    # 95% CI as line
    for(i in 1:NROW(truePar)) lines(rep(truePar[i,p],2), 
                                    posteriorObj$quantiles[[i]][,p], col = rgb(.5,.5,.5,.5))
    if(plotMedian){
      # median as point
      for(i in 1:NROW(truePar)) points(truePar[i,p], 
                                       posteriorObj$median[[i]][p], cex = 0.5,pch = 16) 
    }
  }
}

pNames <- c("habitatPreference", "stepShape", "directionalBias", "maximumEffort")
x11(10,10)
plotParameterReFit(posteriorInS, pNames, PReFit, TRUE, mfrow = c(2,2), pty = "s")
x11(10,10)
plotParameterReFit(posteriorOutS, pNames, PReFit, TRUE, mfrow = c(2,2), pty = "s")

# # habitat preference
# plot(PReFit[,"stepShape"], sapply(1:500, function(x) posteriorInS$median[[x]]["stepShape"]))
# # habitat preference
# plot(PReFit[,"directionalBias"], sapply(1:500, function(x) posteriorInS$median[[x]]["directionalBias"]))
# # habitat preference
# plot(PReFit[,"maximumEffort"], sapply(1:500, function(x) posteriorInS$median[[x]]["maximumEffort"]))

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# plot for out of sample
# real vs. fitted

x11(10,10)
par(mfrow=c(2,2), pty = "s")
# habitat preference
plot(PReFit[,"habitatPreference"], sapply(1:500, function(x) posteriorOutS$median[[x]]["habitatPreference"]))
# habitat preference
plot(PReFit[,"stepShape"], sapply(1:500, function(x) posteriorOutS$median[[x]]["stepShape"]))
# habitat preference
plot(PReFit[,"directionalBias"], sapply(1:500, function(x) posteriorOutS$median[[x]]["directionalBias"]))
# habitat preference
plot(PReFit[,"maximumEffort"], sapply(1:500, function(x) posteriorOutS$median[[x]]["maximumEffort"]))



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##