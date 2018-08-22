# run ABC

runABC <- function(parameters, 
                   habitatSuitability, 
                   observations, 
                   parallel = FALSE, 
                   maxIterations, 
                   randomSeed = NULL,
                   sumStatsSampleSize,
                   sumStats,
                   getSumStats,
                   clusterType = "SOCK" # or "MPI"
                   ){
  # start timing 
  ptm <- proc.time()

  # number of parameter combinations
  nParameters = NROW(parameters)
  
  if(!parallel){
    library(LOSim); library(moments)
    ssSim <- vector("list", length(observations))
    # progress bar
    pb <- txtProgressBar(min = 0, max = length(observations), style = 3)
    for(i in seq_along(observations)){

      # prepare parameters
      parMat <- cbind(rep(observations[[i]]$data[1 , "xObserved"], nParameters),
                      rep(observations[[i]]$data[1 , "yObserved"], nParameters),
                      parameters
      )
      
      # run simulation
      simData <- LOSim::runSimulation(habitatSuitability$values, 
                                      maxIterations, 
                                      parMat, 
                                      randomSeed,
                                      "observed", 
                                      observations[[i]]$data[, "timestamp"], 
                                      NULL,
                                      habitatSuitability$upperleft,
                                      habitatSuitability$resolution,
                                      observations[[i]]$startTime,
                                      observations[[i]]$startDoy,
                                      observations[[i]]$centroid[1],
                                      observations[[i]]$centroid[2]
                                      )
      
      
      # compute summary statistics
      ssSim[[i]] <- getSumStats(simData[[1]], sumStatsSampleSize, simData[[2]], simData[[3]], sumStats)
      
      # progress bar
      setTxtProgressBar(pb, i)
    }
    # close progress bar
    close(pb)
  }else{
    library(doSNOW)
    # parallel execution
    if(clusterType == "MPI"){
      library(doMPI)
      cl <- doMPI::startMPIcluster()
      doMPI::registerDoMPI(cl)
    }else{
      if(parallel == T | parallel == "auto"){
        cores <- parallel::detectCores() - 1
        message("parallel, set cores automatically to ", cores)
      }else if (is.numeric(parallel)){
        cores <- parallel
        message("parallel, set number of cores manually to ", cores)
      }else stop("wrong argument to parallel")
      
      cl <- parallel::makeCluster(cores, type = clusterType)
      doSNOW::registerDoSNOW(cl)
      
    }  
    
    # set up progress bar
    pb <- txtProgressBar(min = 1, max = length(observations), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    ssSim <- foreach::foreach(i = seq_along(observations), 
                              .packages=c("LOSim", "moments"), .options.snow=opts) %dopar% {

                                # prepare parameters
                                parMat <- cbind(rep(observations[[i]]$data[1 , "xObserved"], nParameters),
                                                rep(observations[[i]]$data[1 , "yObserved"], nParameters),
                                                parameters
                                )
                                
                                # run simulation
                                simData <- LOSim::runSimulation(habitatSuitability$values, 
                                                                maxIterations, 
                                                                parMat, 
                                                                randomSeed,
                                                                "observed", 
                                                                observations[[i]]$data[, "timestamp"], 
                                                                NULL,
                                                                habitatSuitability$upperleft,
                                                                habitatSuitability$resolution,
                                                                observations[[i]]$startTime,
                                                                observations[[i]]$startDoy,
                                                                observations[[i]]$centroid[2],
                                                                observations[[i]]$centroid[1])
                                
                                
                                
                                # compute summary statistics
                                getSumStats(simData[[1]], sumStatsSampleSize, simData[[2]], simData[[3]], sumStats)
                              }
    # free memory from workers
    close(pb)
    if(clusterType == "MPI"){
      doMPI::closeCluster(cl)
      Rmpi::mpi.finalize()
    } else{
      parallel::stopCluster(cl)
    }
  } # end parallel
  
  # note settings
  settings <- vector("list", 0L)
  settings$parallel <- parallel
  settings$maxIterations <- maxIterations
  settings$randomSeed <- randomSeed
  settings$sumStatSampleSize <- sumStatsSampleSize
  
  # results
  output <- list(parameters = parameters, 
                 summaries = lapply(seq_along(ssSim), function(x) ssSim[[x]]$summaryStatistics),
                 validity = lapply(seq_along(ssSim), function(x) ssSim[[x]]$validity),
                 #nSimSteps = lapply(seq_along(ssSim), function(x) ssSim[[x]]$nSimSteps),
                 observations = observations,
                 #environment = habitatSuitability,
                 settings = settings
  )
  
  
  # stop timing
  print(proc.time() - ptm)
  
  return(output) 
  
}