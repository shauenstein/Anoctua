calculateDelta <- function(parMatrix, modArgs, trueParameters, summaryStatistics, parNames){
  # run model
  
  dataTrue <- LOSim::runSimulation(modArgs$environment, modArgs$iterations, trueParameters,
                                   modArgs$randomSeed)
  
  summariesTrue = t(sapply(dataTrue, summaryStatistics))
  
  sdTrue <- apply(summariesTrue, 2, sd)
  
  
  data <- LOSim::runSimulation(modArgs$environment, modArgs$iterations, parMatrix, 
                               modArgs$randomSeed)
  
  summaries = t(sapply(data, summaryStatistics))
  
  #sdSummaries <- apply(summaries, 2, sd)
  
  # compute delta / standardise by standard deviation for fixed parameters
  meanSummaryTrue <- colMeans(summariesTrue)
  distance <- apply(summaries, 1, function(x) sqrt(sum((x - meanSummaryTrue)^2 / sdTrue^2)))
  
  
  out = list(parameters = parMatrix, 
             parameterNames = parNames,
             summaries =  summaries, 
             distance = distance, 
             modArgs = modArgs, 
             summaryStatistics = summaryStatistics, 
             trueParameters = trueParameters, 
             summariesTrue = summariesTrue)
  
  # calculate summaries
  return(out)
}
