# simulate environment with gaussian random fields (as in simSAC)
set.seed(1)
lat <- 1000; lon = 1000
expCov <- RandomFields::RMexp(var = 0.1, scale = 0.1)
env <- RandomFields::RFsimulate(expCov, x = seq(0,1,length.out = lon), y = seq(0,1,length.out = lat), spConform = FALSE)
# standardise to 0 ; 1
env <- (env - min(env)) / (max(env) - min(env))

# plot environment
windows(10,10)
dats <- data.frame(expand.grid(seq(lat), seq(lon)), c(env))
colnames(dats) <- c("lat", "lon", "env")
plot(dats$lon, dats$lat, pch=15, cex=1.1, col=terrain.colors(101)[ceiling(dats$env*100)+1], xlab = "Longitude", ylab = "Latitude")
plotrix::color.legend(xl=lat +60, yb=lon -200, xr=lat + 75, yt=lon, legend=c(0,1), rect.col=terrain.colors(101), gradient="y", align="rb", xpd = NA)


# run model on environment
iterations <- 400
startX <- rep(500, 500)
startY <- c(500)
pRanges <- c(15)
nOpts <- c(1)
nRanges <- c(0.4)
param <- as.matrix(expand.grid(startX, startY, pRanges, nOpts, nRanges))
ranSeed <- 1

test <- runSimulation(env,iterations, param, ranSeed)


for(i in seq_along(test)){
  lines(test[[i]][,1], test[[i]][,2], col = rgb(0,0,0,0.2))  # 
}
points(startX[1], startY, pch = "*", cex = 2, col = 2)







###
lat <- 1000; lon = 1000
EnvMat <- matrix(NA, lat, lon)
for(i in seq.int(lat)){
  for(j in seq.int(lon)){
    EnvMat[i,j] <- dnorm(i, lat/2, 1, TRUE) + dnorm(j, lon/2, 1, TRUE)
  }
}
EnvMat <- (EnvMat - min(EnvMat)) / (max(EnvMat) - min(EnvMat))

windows(10,10)
dats <- data.frame(expand.grid(seq(lon), seq(lat)), c(EnvMat))
colnames(dats) <- c("lon", "lat", "env")
plot(dats$lat, dats$lon, pch=15, cex=1.1, col=terrain.colors(101)[ceiling(dats$env*100)+1])
plotrix::color.legend(xl=lat +50, yb=lon -100, xr=lat + 60, yt=lon, legend=c(0,1), rect.col=terrain.colors(101), gradient="y", align="rb", xpd = NA)

#points(dats$lat[dats$env < 0.85 & dats$env > 0.75], dats$lon[dats$env < 0.85 & dats$env > 0.75])

# nicheOpt = 1
# nicheRange = 0.5
# curve(exp(-(nicheOpt - x)^2/nicheRange^2), from = 0, to = 1, n = 200) # check implemented function
# curve(dnorm(x, mean = nicheOpt, sd = nicheRange), from = 0, to = 1, n = 200) # check normal distribution and parameter effect

iterations <- 20
startX <- rep(500, 1)
startY <- c(500)
pRanges <- c(10)
nOpts <- c(1)
nRanges <- c(0.5)
param <- as.matrix(expand.grid(startX, startY, pRanges, nOpts, nRanges))
ranSeed <- 1

test <- runSimulation(EnvMat,iterations, param, ranSeed)

for(i in seq_along(test)){
  lines(test[[i]][,1], test[[i]][,2], col = rgb(0,0,0,0.2))  
}
points(startX[1], startY, pch = "*", cex = 2, col = 2)


# for(l in seq(iterations)){
#   lines(test[[1]][l:(l+1),1], test[[1]][l:(l+1),2], type = "l", col = 2)
#   Sys.sleep(0.8)
#   lines(test[[1]][l:(l+1),1], test[[1]][l:(l+1),2], type = "l", col = 1)
#   
# }


### binary environment
lat <- 100; lon = 100
EnvMat <- matrix(0, lat, lon)
# change values in the center
EnvMat[(0.25*lat):(0.75*lat),(0.25*lon):(0.75*lon)] <- 1

dats <- data.frame(expand.grid(seq(lon), seq(lat)), c(EnvMat))
colnames(dats) <- c("lat", "lon", "env")
plot(dats$lon, dats$lat, pch=15, cex=1.1, col=ifelse(EnvMat == 1, "orangered", "yellow"))


iterations <- 200
startX <- c(50)
startY <- c(50)
pRanges <- c(5)
nOpts <- c(1)
nRanges <- c(0.01)
param <- as.matrix(expand.grid(startX, startY, pRanges, nOpts, nRanges))
ranSeed <- 1

test <- runSimulation(EnvMat,iterations, param, ranSeed)

points(startX, startY, col = "White", pch = "+", cex = 1.5)
for(i in seq_along(test)){
  lines(test[[i]][,1], test[[i]][,2], col =i)  
}
