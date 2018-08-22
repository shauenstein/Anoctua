abcCreateSummaries <- function(simulationOutput, 
                               method = "saABC", 
                               summarySelection = NULL, 
                               targetParameters = NULL, 
                               getSumStats,
                               sumStats,
                               sumStatsSampleSize,
                               parallel = FALSE,
                               predictionPredictor = TRUE,
                               clusterType = "SOCK" # or "MPI"
                               ){
  
  # start timing 
  ptm <- proc.time()
  
  # simulationOutput$info = returnInfo(simulationOutput$info, note = "abcCreateSummaries")
  ssData <- getSumStats(simulations = lapply(seq_along(simulationOutput$observations),
                                             function(x) simulationOutput$observations[[x]]$data),
                        sampleSize = sumStatsSampleSize, 
                        sumStats = sumStats)
  
  if(is.null(summarySelection)) summarySelection = 1:NCOL(simulationOutput$summaries[[1]])
  if(is.null(targetParameters)) targetParameters = 1:NCOL(simulationOutput$parameters)
  
  # output lists
  simPreds <- dataPreds <- vector("list", length(simulationOutput$summaries))
  
  # saABC: linear regression
  if(method == "saABC"){
    # predict function
    predict <- function(summaries){
      if(is.vector(summaries)) summaries = matrix(summaries, nrow = 1)
      summaries[is.na(summaries)] = 0
      summaries = cbind(summaries, summaries^2)
      res = summaries %*% t(out$B)
      colnames(res) = paste0("S", seq_along(targetParameters))
      return(res)
    }
    
    # loop through observed individuals
    # in parallel
    if(parallel){
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
      pb <- txtProgressBar(min = 1, max = length(simulationOutput$summaries), style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      
      preds <- foreach::foreach(i = seq_along(simulationOutput$summaries), 
                                .packages=c("abctools"), .options.snow=opts) %dopar% {
                                  tempSum = simulationOutput$summaries[[i]]
                                  tempSum[is.na( tempSum)] = 0
                                  
                                  out = abctools::saABC(theta = simulationOutput$parameters[,targetParameters], 
                                                        X = cbind(tempSum[,summarySelection], tempSum[,summarySelection]^2),
                                                        plot = F)
                                  
                                  out$B[is.na(out$B) ] = 0
                                  
                                  
                                  list(simPreds = predict(simulationOutput$summaries[[i]][, summarySelection]),
                                       dataPreds = predict(ssData[[1]][i, summarySelection]))
                        
                                }
      # free memory from workers
      close(pb)
      if(clusterType == "MPI"){
        doMPI::closeCluster(cl)
        Rmpi::mpi.finalize()
      } else{
        parallel::stopCluster(cl)
      }
      for(i in seq_along(preds)){
        simPreds[[i]] <- preds$simPreds[[i]]
        dataPreds[[i]] <- preds$dataPreds[[i]]
      }
    }else{
      for(i in seq_along(simulationOutput$summaries)){
        
        tempSum = simulationOutput$summaries[[i]]
        tempSum[is.na( tempSum)] = 0
        
        out = abctools::saABC(theta = simulationOutput$parameters[,targetParameters], 
                              X = cbind(tempSum[,summarySelection], tempSum[,summarySelection]^2), 
                              plot = F)
        
        out$B[is.na(out$B) ] = 0
        
        simPreds[[i]] <- predict(simulationOutput$summaries[[i]][, summarySelection])
        dataPreds[[i]] <- predict(ssData[[1]][i, summarySelection])
        cat("individual No. ",i, "\n")
      }
    }
  }
    
  
  # ranger randomForest
  if(method == "rF"){
    library(ranger)
    if(!parallel) parallel <- 1
    if(parallel == TRUE) parallel <- parallel::detectCores() - 1
    for(i in seq_along(simulationOutput$summaries)){
      tempSum = simulationOutput$summaries[[i]]
      tempSum[is.na( tempSum)] = 0
      
      rfsummaries <- vector("list", length(targetParameters))
      simRFPreds <- as.data.frame(matrix(0, ncol = length(targetParameters), 
                                         nrow = NROW(tempSum)))
      colnames(simRFPreds) <- paste0("S", seq_along(targetParameters))
      dataRFPreds <- as.data.frame(matrix(0, ncol = length(targetParameters), 
                                          nrow = 1))
      colnames(dataRFPreds) <- paste0("S", seq_along(targetParameters))
      
      index = 1
      for(j in targetParameters){
        if(predictionPredictor){
          TrainData <- data.frame(y = simulationOutput$parameters[,j], 
                                  tempSum[,summarySelection], 
                                  simRFPreds)
          PredDataTrue <- as.data.frame(matrix(c(ssData[[1]][i, summarySelection],
                                                 dataRFPreds), nrow = 1))
        }else{
          TrainData <- data.frame(y = simulationOutput$parameters[,j], 
                                  tempSum[,summarySelection])
          PredDataTrue <- as.data.frame(matrix(ssData[[1]][i, summarySelection], nrow = 1))
        }
        
        colnames(PredDataTrue) <- colnames(TrainData)[-1]
        rfsummaries[[index]] <- ranger::ranger(y ~ ., data = TrainData, 
                                               num.trees = 500, write.forest = TRUE, 
                                               num.threads = parallel, verbose = FALSE)
        
        simRFPreds[,index] <- predict(rfsummaries[[index]], 
                                      data = TrainData, 
                                      num.threads = parallel)$predictions
        # simRFPreds[,index] <- rfsummaries[[index]]$predictions
        dataRFPreds[,index] <- predict(rfsummaries[[index]], 
                                                      data = PredDataTrue, 
                                                      num.threads = parallel)$predictions  
        
        index <- index + 1
      }
      simPreds[[i]] <- simRFPreds
      dataPreds[[i]] <- dataRFPreds
      
      cat("individual No. ",i, ": ", proc.time()[3] - ptm[3], "sec.","\n")
    }
  }  
  
  # collect results
  output <- list()
  output$parameters = simulationOutput$parameters
  output$observations = simulationOutput$observations
  #output$summarySelection$predict = predict # exports the function enviornment, clamps up disc memory
  output$summarySelection$predictionPredictor = predictionPredictor
  output$summarySelection$simulation = simPreds
  output$summarySelection$observed = dataPreds
  output$summarySelection$method = method
  output$summarySelection$targetParameters = targetParameters
  output$summarySelection$summarySelection = summarySelection
  
  
  
  # stop timing
  print(proc.time() - ptm)
  
  return(output)
}

