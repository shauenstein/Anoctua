pt <- proc.time()
nSteps <- 200
data <- randomWalk(startValues = c(0,0,0), param = c(2,1), steps = nSteps)
data <- cbind(data, rep(NA, nSteps+1),rep(NA, nSteps+1))

## Observation model

observationModel <- function(realData, sd=1){
  realData[,3] = rnorm(NROW(realData), mean = realData[,1], sd = sd)
  realData[,4] = rnorm(NROW(realData), mean = realData[,2], sd = sd)
  realData[realData[,3] - floor(realData[,3]) > 0.7 , 3] = NA
  return(realData)
}

obsdata <- observationModel(data)
plot(data, type = "l")
points(obsdata[,3], obsdata[,4], col = "red", pch = 4)



## Summary statistics
summaryStatistics <- function(dat){
  meandisplacement = mean(sqrt(diff(dat[,3])^2 + diff(dat[,4])^2), na.rm = T)
  
  meandisplacement10 = mean(sqrt(diff(dat[,3], lag = 2)^2 + diff(dat[,4], lag = 2)^2), na.rm = T)/3
  
  #meanturning = mean(abs(diff(atan2(diff(dat[,4]),diff(dat[,3])))), na.rm = T) 
  
  return(c(meandisplacement, meandisplacement10))
} 

dataSummary <- summaryStatistics(obsdata)
dataSummary


## ABC rejection algorithm

n = 10000
fit = data.frame(movementLength = runif(n, 0, 5), error = runif(n, 0,5), distance = rep(NA, n))

for (i in 1:n){
  simulatedPrediction <- randomWalk(startValues = c(0,0,0), param = c(fit[i,1], 1), steps = nSteps)
  simulatedPrediction <- cbind(simulatedPrediction, rep(NA, nSteps+1),rep(NA, nSteps+1))
  simulatedObservation<- observationModel(simulatedPrediction, fit[i,2])
  simulatedSummary <- summaryStatistics(simulatedObservation)
  simulatedSummary
  #deviation = max( simulatedSummary - dataSummary)
  deviation = sqrt(sum((simulatedSummary - dataSummary)^2))
  fit[i,3] = deviation
}

## I had already calculated the euclidian distance between observed and simulated summaries. We now plot parameters for different acceptance intervals

plot(fit[fit[,3] < 1, 1:2], xlim = c(0,5), ylim = c(0,2), col = "lightgrey", main = "Accepted parameters for \n different values of epsilon")
points(fit[fit[,3] < 0.2, 1:2],  pch = 18, col = "gray")
points(fit[fit[,3] < 0.1, 1:2],  pch = 8, col = "red")

legend("topright", c("< 1", "< 0.2", "< 0.1"), pch = c(1,18,8), col = c("lightgrey", "gray", "red"), bty = "n")

abline(v = 2)
abline(h = 1) 


## ABC-MCMC Algorithm

n = 10000
fit = data.frame(movementLength = rep(NA, n), error = rep(NA, n), distance = rep(NA, n))


currentPar = c(2,0.9)
for (i in 1:n){
  newPar = rnorm(2,currentPar, sd = c(0.2,0.2))
  if (min(newPar) < 0 ) fit[i,] = c(currentPar, deviation)
  
  else{
    simulatedPrediction <- randomWalk(startValues = c(0,0,0), param = c(newPar[1], 1), steps = 200)
    simulatedPrediction <- cbind(simulatedPrediction, rep(NA, nSteps+1),rep(NA, nSteps+1))
    simulatedObservation<- observationModel(simulatedPrediction, newPar[2])
    simulatedSummary <- summaryStatistics(simulatedObservation)
    deviation = sqrt(sum( simulatedSummary - dataSummary)^2)
    
    if (deviation < 0.2){
      fit[i,] = c(newPar, deviation)
      currentPar = newPar
    } 
    else fit[i,] = c(currentPar, deviation)
  }
}

plot(fit[, 1:2], xlim = c(0,9), ylim = c(0,2), col = "#00000022", main = "Accepted parameters")

abline(v = 2)
abline(h = 1) 
pt <- proc.time() - pt