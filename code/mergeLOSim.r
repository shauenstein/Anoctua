# merge simulation outputs
mergeLOSim <- function(sim1, sim2){
  if(sum(!names(sim1$observations) %in% names(sim2$observations)) != 0 | 
     length(sim1$observations) != length(sim2$observations)){
    stop("Only merge simulation outputs from the same individuals!") 
  } 
  
  sim1$parameters <- rbind(sim1$parameters, sim2$parameters)
  for(i in seq_along(sim1$summaries)){
    sim1$summaries[[i]] <- rbind(sim1$summaries[[i]], sim2$summaries[[i]])
    sim1$validity[[i]] <- c(sim1$validity[[i]], sim2$validity[[i]])
    sim1$nSimSteps[[i]] <- c(sim1$nSimSteps[[i]], sim2$nSimSteps[[i]])
  }
  
  return(sim1)
}