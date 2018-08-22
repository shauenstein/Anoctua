
####################################
# add buffer of zeros function
addZ <- function(mat, k=1){
  emat <- matrix(0, nrow = NROW(mat)+k*2, ncol = NCOL(mat)+k*2)
  emat[(k+1):(NROW(mat)+k),(k+1):(NCOL(mat)+k)] <- mat
  return(emat)
}
####################################
# set randomisation seed
set.seed(14)

# load packages
library(LOSim)
library(RandomFields)
library(plotrix)
library(randomForest)
library(raster)

lat <- 200; lon = 200 # number of gridcells in latitudinal and longitudinal direction

#----------------------------------------------------------------------------------------------#

# gaussian random fields
expCov <- RMexp(var = 1, scale = .05) # exponential covariance model
env <- RFsimulate(expCov, x = seq(0,1,length.out = lon), # simulate random fields
                  y = seq(0,1,length.out = lat), spConform = FALSE) # using the exp-cov model

# scale to [0;1]
env <- (env - min(env)) / (max(env) - min(env))

#----------------------------------------------------------------------------------------------#

# circle
env <- matrix(NA, nrow = lat, ncol = lon)
for(i in seq.int(lat)){
  for(j in seq.int(lon)){
    env[i,j] <- dnorm(i, mean = lat/2, sd = lat/3) *
      dnorm(j, mean = lon/2, sd = lon/3)
  }
}

# scale to [0;1]
env <- (env - min(env)) / (max(env) - min(env))

#----------------------------------------------------------------------------------------------#

# random uniform
env <- matrix(runif(lat*lon, 0, 0.5), nrow = lat, ncol = lon)
env[c(10,12,14,14,12,10),c(10,12,14,14,12,10)] = runif(6,0.8,1)
#----------------------------------------------------------------------------------------------#

#add buffer of zeros
k = 10
env <- addZ(mat = env, k)

# uniform env
# env <- matrix(runif(lat*lon, 0, 1), lat, lon)


#plot env # green is high, white is low
# dats <- data.frame(expand.grid(seq(lat + k*2), seq(lon + k*2)), c(env)) # transform to 2 col matrix
# colnames(dats) <- c("lat", "lon", "env")
# plot(dats$lon, dats$lat, pch=15, cex=1.1, col=terrain.colors(101)[ceiling(dats$env*100)+1], xlab = "Longitude", ylab = "Latitude")
# color.legend(xl=lat, yb=lon-100, xr=lat + 15, yt=lon, legend=c(0,1), rect.col=terrain.colors(101), gradient="y", align="rb", xpd = NA)

# as raster
rast <- raster(env, xmn = 1, xmx = lat + 2*k, ymn = 1, ymx = lon + 2*k)
plot(rast)

# parameters
param <- as.matrix(data.frame(xstart = (lon + k*2)/2 - 90,
                              ystart = (lat + k*2)/2 + 90,
                              habitatPreference=.3, 
                              stepShape=1, 
                              stepScale=2, 
                              directionalBias=.01, 
                              dispRestMean=1.5, 
                              dispRestSd=.7, 
                              maximumEffort=6000, 
                              roostLambda =2
                              ))


envir <- list(values = env, upperleft = c(1,NROW(env)), resolution = 1)

# run simulation
# plot(1:10, type = "n",xlim=c(0,1),ylim=c(0,20))
for(i in 1:10){
  testsim <- LOSim::runSimulation(environment = envir$values, 
                                  iterations = 1000, 
                                  parameters =param, 
                                  randomSeed = NULL,
                                  output = "raw", 
                                  TimeAtObservation = NULL,
                                  maxDuration = 60,
                                  upperleft = envir$upperleft,
                                  environmentResolution = envir$resolution,
                                  startTime = 12,
                                  startDayOfYear = 250,
                                  solarLatitude = 9,
                                  solarLongitude = 48
  )
  
  # lines(density(testsim$simulation[[1]][,"habitatSuitability"]), col =3)
  lines(testsim$simulation[[1]][,c("x","y")], col =6)
}


#plot(1,type = "n", ylim = c(0,1000), xlim = c(0,1000), xlab = "x", ylab = "y")
# plot(test[[1]][,2], test[[1]][,3], col = 1, type = "l")
# points(test[[1]][test[[1]][,8] == 1 ,2], test[[1]][test[[1]][,8] == 1,3], col = 2)
for(i in seq(n)) lines(test[[i]][,2], test[[i]][,3], col = rgb(0,0,0,0.2))
points(500,500, col = 2)

####################################
# check summaries correlation
windows(10,10)
plot(Hmisc::varclus(trialRun$summaries[[2]][, c(1,2,5,7,9,10,12,13,14,16,17,18,19,26,27,28,29,30,36,37,38,39,50,59,60)])) # yielding the figure below
abline(h=0.5, col="darkred", lwd=2, lty=2) 
####################################

# habitat preference
hP = seq(0.1,3,0.02)
curve(x^(0.1*exp(0.1)), xlim = c(0,1), ylim = c(0,1))
for(i in hP) curve(x^(i*exp(i)), xlim = c(0,1),add = T)


# time at roost site
envValAtIterp1 = 0.1
maxRoostTime = 3*24*60
curve(exp(- x / envValAtIterp1) * maxRoostTime, from = 0, to = 1)
for(i in seq(0.2,1,0.1)) curve(exp(- x / i) * maxRoostTime, from = 0, to = 1, add = T)


# energy gain: ((1/(1 + exp(- pow(cumEnv, recoveryRate) * ((timestamp.at(iter + 1) - timestamp.at(iter)))))) // energy uptake
#- 0.5) / 0.5
curve(((1/(1 + exp(- (0.5/20)^2 * (x*24*60)))) - 0.5) / 0.5, from = 0, to = 5)
# energy loss
curve(exp( - 5*20 / x), from = 0, to = 100)
# not exponential
curve(x/6000, from = 0, to = 6000, ylim = c(0,1))
abline(h = 0.5)


# step length
curve(pgamma(x, shape = 0.5, scale = 2), from = 0, to = 20, ylim=c(0,1))
curve(dgamma(x, shape = 4, scale = 2), from = 0, to = 20, ylim=c(0,1))
abline(h = 0.95)

# directional bias
curve(dnorm(x, sd = 1/(exp(2)-1)), xlim=c(-pi,pi))
for(i in rev(seq(0.1,2,0.1))) curve(dnorm(x, sd = 1/(exp(i)-1)), add=T)

# prob to roost
curve(exp(-1/x*(1*0.1)), from = 0, to = 1, ylim = c(0,1))
for(i in seq(0.2,1,0.1)){
  curve(exp(-1/x*(1*i)), from = 0, to = 1, add = T)
  Sys.sleep(2)
}


######
iterations <- 200
startX <- 500
startY <- 500
pRanges <- 8
nOpts <- 0.8
nRanges <- 0.1
steplengthScale <- exp(seq(0.1,5, 0.1))-1
steplengthLocationFactor <- seq(0.01,1, 0.01) # factor of preception range # 0.5 = half pR 
obsErrorReal <- 2
parReal <- as.matrix(expand.grid(startX, startY, pRanges, nOpts, nRanges, steplengthScale, steplengthLocationFactor, obsErrorReal))




test <- LOSim::runSimulation(env, iterations, parReal, 14)

# plot environment
windows(10,10)
dats <- data.frame(expand.grid(seq(lat), seq(lon)), c(env))
colnames(dats) <- c("lat", "lon", "env")
plot(dats$lon, dats$lat, pch=15, cex=1.1, col=terrain.colors(101)[ceiling(dats$env*100)+1], xlab = "Longitude", ylab = "Latitude")
color.legend(xl=lat +60, yb=lon -200, xr=lat + 75, yt=lon, legend=c(0,1), rect.col=terrain.colors(101), gradient="y", align="rb", xpd = NA)


for(i in seq_along(test)){
  lines(test[[i]]$x, test[[i]]$y, col = 1)
  Sys.sleep(2)
} 


######


dispMat <- matrix(NA, nrow = iterations, ncol = length(test))
for(i in seq_along(test)){
  dispMat[,i] <- sqrt(diff(test[[i]]$xObs)^2 + diff(test[[i]]$yObs)^2)
}

sd <- means <- skew <- 1:3
for(i in seq_along(test)){
  sd[i] <- sd(dispMat[,i])
  means[i] <- mean(dispMat[,i])
#   geary[i] <- moments::geary(dispMat[,i])
  # kurt[i] <- e1071::kurtosis(dispMat[,i])
  skew[i] <- e1071::skewness(dispMat[,i])
}

par(mfrow = c(1,2))
plot(log(parReal[,6]), (sd) )
plot(parReal[,7], means)
plot(parReal[,7], skew )

plot(density(dispMat[,1]), ylim = c(0,0.15))
for(i in seq(NCOL(dispMat))[-1]){
  lines(density(dispMat[,i]))
  Sys.sleep(0.5)
}

####

## sample parameters

iterations <- 300
startX <- 500
startY <- 500
pRanges <- 8
nOpts <- 0.8
nRanges <- 0.2
steplengthFactor <- 0.5
obsErrorReal <- 2
parReal <- matrix(rep(c(startX, startY, pRanges, nOpts, nRanges, 
                        steplengthFactor, obsErrorReal), 1000), 
                  nr = 1000, byrow = TRUE)


n = 1000
parSample <- cbind(rep(500, n), # starX
                   rep(500, n), # startY
                   sample(1:15, size = n, replace = TRUE), # perception range
                   rep(nOpts, n), # niche optimum
                   rep(nRanges, n),#runif(n, 0.05, 0.5), # niche range // cannot be 0
                   runif(n, -1, 1), # steplength factor
                   runif(n, 0, 5)) # observation error

dataTrue <- LOSim::runSimulation(env, iterations, parReal, 1)

summariesTrue = t(sapply(dataTrue, summaryStatistics))
sdTrue <- apply(summariesTrue, 2, sd)

data <- LOSim::runSimulation(env, iterations, parSample, 14)

summaries = t(sapply(data, summaryStatistics))
meanSummaryTrue <- colMeans(summariesTrue)
summariesScaled <- summaries / sdTrue
summariesTrueScaled <- meanSummaryTrue / sdTrue


rfPar1 <- randomForest(y = parSample[,3], x = summariesScaled)
rfPar2 <- randomForest(y = parSample[,6], x = summariesScaled)
rfPar3 <- randomForest(y = parSample[,7], x = summariesScaled)

predsP1 <- predict(rfPar1)
predsP2 <- predict(rfPar2)
predsP3 <- predict(rfPar3)

preds <- cbind(predsP1,predsP2,predsP3)
distance <- apply(preds, 1, function(x) sqrt(sum((x - summariesTrueScaled)^2)))



# test cpp distance function
Rcpp::sourceCpp('../code/c++/distanceObservedSimulated.cpp')

set.seed(1)
sim <- data.frame(S1 = runif(10000), S2 = runif(10000), S3 = runif(10000), S4 = runif(10000))
obs <- data.frame(S1 = runif(1), S2 = runif(1), S3 = runif(1), S4 = runif(1))
reference <- diff(apply(sim, 2, range))

system.time(dcpp <- dObsSim(as.matrix(sim), unlist(obs), unlist(reference)))
system.time(dr <- apply(sim, 1, function(x) sqrt(sum(((x - obs[1,]) 
                                              / reference)^2))))


## multinomial selection for day-rest, 1 night, 2 nights
x11(8,4); par(mfrow=c(1,2))
#energy
h = 1
curve(0.02 + 0.85*x ,xlim=c(0,1),ylim=c(0,1))
curve(0.5 + 0.48*x, add = T)
# curve(plogis(-6 + 9*x),xlim=c(0,1),ylim=c(0,1))
# curve(plogis(-2 + 7*x),xlim=c(0,1),ylim=c(0,1), add = T)
# abline(v=0.5)
# #habitat
# curve(plogis(6 - 9*x),xlim=c(0,1),ylim=c(0,1))
# curve(plogis(2 - 7*x),xlim=c(0,1),ylim=c(0,1), add = T)
# abline(v=0.5)


library(plot3D)
# energy has negative effect, HS positive effect on prob of roosting
# HS = 0
# day rest
b01 = 0.22; b11 = 0.7;b21 = -0.2
# night rest
b02 = 0.5; b12 = 0.48;b22 = -0.15
hs = en = seq(0,1,length.out = 300)
M <- mesh(hs, en)
x <- M$x; y <- M$y
x11()
surf3D(x,y, z = b01 + b11 * x + b21 * y, bty = "b2", xlab = "energy", ylab = "habitat suit.", zlab = "", theta = -30, phi = 0, zlim=c(0,1), col = "blue")
surf3D(x,y, z = b02 + b12 * x + b22 * y, col = "red", add = T)

