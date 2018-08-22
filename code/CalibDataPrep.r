# transform telemetry data to suitable format for calibration
calibTransform <- function(x, y, date, id, habitatSuitability, 
                           epsg.input = "+init=epsg:31467", 
                           epsg.output = "+init=epsg:4326",
                           rscExtent = 100){
  #library(raster); library(lubridate)
  
  # load rcpp function to calculate distance between observed and simulated data
  Rcpp::sourceCpp("../code/c++/rscCalculateForObs.cpp")
  
  # lists matrices for each id
  outList <- vector(mode = "list", length = length(unique(id)))
  
  # progress bar
  pb <- txtProgressBar(min = 0, max = length(unique(id)), style = 3)
  
  for(i in seq_along(unique(id))){
    idStartT <- lubridate::hour(date[id == unique(id)[i]][1]) + 
      lubridate::minute(date[id == unique(id)[i]][1])/60 + 
      lubridate::second(date[id == unique(id)[i]][1])/3600
    idDoy <- lubridate::yday(date[id == unique(id)[i]][1])
    idLatLon <- matrix(c(mean(x), mean(y)),nrow=1)
    idLatLon <- as.vector(coordinates(sp::spTransform(sp::SpatialPoints(idLatLon, CRS(epsg.input)), CRS(epsg.output))))
    
    idMat <- matrix(NA, nrow = sum(id == unique(id)[i]), ncol = 7)
    
    idMat[,1] <- (as.numeric(date[id == unique(id)[i]]) - 
                    as.numeric(min(date[id == unique(id)[i]])))/60 # in min
    idMat[,2] <- x[id == unique(id)[i]]
    idMat[,3] <- y[id == unique(id)[i]]
    for(j in seq(NROW(idMat))){
      if(j == 1){
        idMat[j,4] <- NA
      }else{
        idMat[j,4] <- sqrt((idMat[j,2] - idMat[(j-1),2])^2 + 
                             (idMat[j,3] - idMat[(j-1),3])^2)
      }
      if(j < 3){
        idMat[j,5] <- NA
      }else{
        prev <- c((idMat[(j-1),2] - idMat[(j-2),2]),
                  (idMat[(j-1),3] - idMat[(j-2),3]))
        
        curr <- c((idMat[j,2] - idMat[(j-1),2]),
                  (idMat[j,3] - idMat[(j-1),3]))
        
        # if no movement in one of the steps
        if((prev[1] == 0 & prev[2] == 0) |
           (curr[1] == 0 & curr[2] == 0)){
          idMat[j,5]  <- NA
        }else{
          # if turning back exactly 180 degree, i.e. pi
          if(sum(prev + curr) == 0 & diff(prev + curr) == 0){
            idMat[j,5]  <- pi
          }else{
            idMat[j,5] <-  acos((curr[1]*prev[1] + curr[2]*prev[2]) /
                                  (sqrt(curr[1]^2 + curr[2]^2) * sqrt(prev[1]^2 + prev[2]^2))) 
          }
        }
      }
    }
    idMat[,6] <- raster::extract(habitatSuitability, idMat[,c(2,3)])
    upperleft <- raster::extent(habitatSuitability)[c(1,4)]
    envRes <- raster::res(habitatSuitability)[1]
    cellcoordinates <- idMat[,c(2,3)]
    cellcoordinates[,1] <- (cellcoordinates[,1] - upperleft[1]) / envRes
    cellcoordinates[,2] <- -(cellcoordinates[,2] - upperleft[2]) / envRes
    idMat[,7] <- rsc_calc_obs(as.matrix(habitatSuitability), cellcoordinates, 
                              floor(rscExtent / envRes))
    
    colnames(idMat) <- c("timestamp", "xObserved", 
                         "yObserved", "stepDistanceObserved", 
                         "turningAngleObserved", "habitatSuitabilityObserved", "rsc")
    
    outList[[i]] <- list(data = idMat, startTime = idStartT, 
                         startDoy = idDoy, centroid = idLatLon)
    
    setTxtProgressBar(pb, i)
  }
  names(outList) <- unique(id)
  # close progress bar
  close(pb)
  return(outList)
}