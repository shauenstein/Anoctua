# get summary statistics for simulated observations
getSumStats <- function(simulations, sampleSize, validityCode = NULL, nSteps = NULL, sumStats){
  if(mode(simulations) != "list") stop("Please provide the simulations as a list.")
  ssOut <- t(sapply(seq_along(simulations), function(x) sumStats(simulations[[x]], n = sampleSize)))
  return(list(summaryStatistics = ssOut, validity = validityCode, nSimSteps = nSteps))
}