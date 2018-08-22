# forward simulation function 
# samples parameter combinations from LOSim posterior object
LOSim <- function(posteriorObj, 
                  habitatSuitability,
                  individuals = NULL, 
                  startLocations, 
                  fixParameters,
                  maxPeriod = 60, 
                  maxIterations = 40000, 
                  generations = 10,
                  randomSeed = NULL,
                  reflection = 20,
                  reflectionValues = 1e-20,
                  parallel = FALSE,
                  clusterType = "SOCK" # MPI
                  ){
  
  # start timing 
  ptm <- proc.time()
  
  library(LOSim)
  
  # add ring of very small values around habitatSuitability map to enable reflection
  if(is.numeric(reflection)){
    # add ring of reflection cells with zeros as buffer zone to enable reflection
    rowZeros <- matrix(reflectionValues, nrow = reflection, ncol = NCOL(habitatSuitability$values))
    colP2Zeros <- matrix(reflectionValues, nrow = NROW(habitatSuitability$values)+2*reflection, ncol = reflection)
    #
    habitatSuitability$values <- rbind(rowZeros, habitatSuitability$values, rowZeros)
    habitatSuitability$values <- cbind(colP2Zeros, habitatSuitability$values, colP2Zeros)
    
    # edit upperleft
    habitatSuitability$upperleft <- habitatSuitability$upperleft + c(-1,1) * reflection*habitatSuitability$resolution
    
  }
  
  
  # number of parameter combinations per individual
  posteriorSize <- NROW(posteriorObj$parameters[[1]])
  # if only for a subset of individuals
  if(!is.null(individuals) & is.numeric(individuals)) posteriorObj$parameters <- posteriorObj$parameters[individuals]
  
  # number of virtual individuals to create
  n <- NROW(startLocations)
  
  #--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#
  
  
  if(!parallel){
    # initialisation
    SimList <- vector("list", generations)
    names(SimList) <- paste0("generation",1:generations)
    # start generation loop
    # progress bar
    pb <- txtProgressBar(min = 0, max = generations, style = 3)
    for(g in seq.int(generations)){
      
      # sample individuals
      if(!is.null(individuals)){
        ISample <- sample(individuals, 
                          size = n, replace = TRUE)
      }else{
        ISample <- sample(length(posteriorObj$parameters), 
                          size = n, replace = TRUE)
      }
      
      # sample from individual specific posterior distributions
      PSample <- t(sapply(ISample, 
                          function(x) posteriorObj$parameters[[x]][sample(posteriorSize, size = 1),]))
      
      #interspace with fix parameters
      Pmat <- matrix(NA, nrow = n, 
                           ncol = 2+length(fixParameters)+length(posteriorObj$targetParameters))
      
      Pmat[,2+posteriorObj$targetParameters] <- PSample
      
      fIndex <- 1
      for(f in seq.int(2+length(fixParameters)+
                       length(posteriorObj$targetParameters))[-c(1,2,posteriorObj$targetParameters+2)]){
        Pmat[,f] <- rep(fixParameters[fIndex], n)
        fIndex <- fIndex + 1
      }
      
      # first generation
      if(g == 1){
        # add initial start coordinates to PMat
        Pmat[1:n,1:2] <- as.matrix(startLocations[,c("x","y")])
      }else{
        # use end coordinates of last run as start coordinates
        Pmat[1:n,1:2] <- t(sapply(seq_along(SimList[[g-1]]), 
                                  function(x) SimList[[g-1]][[x]][NROW(SimList[[g-1]][[x]]), 
                                                                  c("x", "y")]))
      }
      
      
  
      SimList[[g]] <- LOSim::runSimulation(habitatSuitability$values, 
                                           maxIterations, 
                                           Pmat, 
                                           randomSeed,
                                           "raw", 
                                           NULL, 
                                           maxPeriod,
                                           habitatSuitability$upperleft,
                                           habitatSuitability$resolution,
                                           startTime = runif(1,0,23.9), # somewhere between 0 and 24 
                                           startDayOfYear = sample(205:352, size = 1), # empiric take off doys
                                           solarLatitude = round(mean(startLocations[,2])), 
                                           solarLongitude = round(mean(startLocations[,1])))[["simulation"]]
      
      # progress bar
      setTxtProgressBar(pb, g)
    }
    close(pb)
  }else{
    library(doSNOW)
    # parallel execution
    if(clusterType == "MPI"){
      library(doMPI)
      cl <- doMPI::startMPIcluster()
      cores <- cl$workerCount
      print(cores)
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
    
    # split parameter matrix for parallel processing
    chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))
    rows <- chunk(1:n, cores)
    
    mergeList <- function(...) {
      mapply(c, ..., SIMPLIFY=FALSE)
    }
    
    # set up progress bar
    pb <- txtProgressBar(min = 1, max = length(rows), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
  # simulate in parallel
    SimList <- foreach(i = seq_along(rows), .combine = mergeList,
                       .packages=c("LOSim"),.options.snow=opts) %dopar% {
                         prelist <- vector("list", generations)
                         names(prelist) <- paste0("generation",1:generations)
                         for(g in seq.int(generations)){
                           # sample individuals
                           if(!is.null(individuals)){
                             ISample <- sample(individuals, 
                                               size = length(rows[[i]]), replace = TRUE)
                           }else{
                             ISample <- sample(length(posteriorObj$parameters), 
                                               size = length(rows[[i]]), replace = TRUE)
                           }
                           
                           # sample from individual specific posterior distributions
                           PSample <- t(sapply(ISample, 
                                               function(x) posteriorObj$parameters[[x]][sample(posteriorSize, 
                                                                                               size = 1),]))
                           
                           #interspace with fix parameters
                           Pmat <- matrix(NA, nrow = length(rows[[i]]), 
                                          ncol = 2+length(fixParameters)+length(posteriorObj$targetParameters))
                           
                           Pmat[,2+posteriorObj$targetParameters] <- PSample
                           
                           fIndex <- 1
                           for(f in seq.int(2+length(fixParameters)+
                                            length(posteriorObj$targetParameters))[-c(1,2,
                                                                                      posteriorObj$targetParameters+2)]){
                             Pmat[,f] <- rep(fixParameters[fIndex], length(rows[[i]]))
                             fIndex <- fIndex + 1
                           }
                           
                           # first generation
                           if(g == 1){
                             # add initial start coordinates to PMat
                             Pmat[,1:2] <- as.matrix(startLocations[rows[[i]],c("x", "y")])
                           }else{
                             # use end coordinates of last run as start coordinates
                             Pmat[,1:2] <- t(sapply(seq_along(prelist[[g-1]]), 
                                                       function(x) prelist[[g-1]][[x]][NROW(prelist[[g-1]][[x]]), 
                                                                                       c("x", "y")]))
                           }
                           prelist[[g]] <- LOSim::runSimulation(habitatSuitability$values, 
                                                maxIterations, 
                                                Pmat, 
                                                randomSeed,
                                                "raw", 
                                                NULL, 
                                                maxPeriod,
                                                habitatSuitability$upperleft,
                                                habitatSuitability$resolution,
                                                startTime = runif(1,0,23.9), # somewhere between 0 and 24 
                                                startDayOfYear = sample(205:352, size = 1), # empiric take off doys
                                                solarLatitude = round(mean(startLocations[,2])), 
                                                solarLongitude = round(mean(startLocations[,1])))[["simulation"]]
                         }
                         prelist
                        }
    
    # close cluster
    close(pb) # and progress bar
    if(clusterType == "MPI"){
      doMPI::closeCluster(cl)
      Rmpi::mpi.finalize()
    } else{
      parallel::stopCluster(cl)
    }
    
    # reorganise simList
    
  }

  
  res <- NULL
  # res$parameters <- PList
  res$simulations <- SimList
  res$startLocations <- startLocations
  #res$environment <- habitatSuitability
  res$dispersalduration <- maxPeriod
  res$individuals <- individuals
  res$reflection$ncell <- reflection
  res$reflection$vals <- reflectionValues
  
  # stop timing
  print(proc.time() - ptm)
  
  return(res)
}


