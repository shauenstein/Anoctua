
# summary statistics
sumStats <- function(observed, n){
  # habitat preference
  ss1 <- quantile(observed[, "habitatSuitabilityObserved"], probs = c(0.1,0.9), na.rm = TRUE) # ,0.2,0.8
  ss2 <- sd(observed[, "habitatSuitabilityObserved"], na.rm = TRUE)
  ss3 <- mean(observed[, "habitatSuitabilityObserved"], na.rm = TRUE)
  
  # gamma distribution: step length, perception range
  from <- sample(NROW(observed), size = n, replace = TRUE)
  to <- sample(NROW(observed), size = n, replace = TRUE)
  deltaT <- abs(observed[to, "timestamp"] - observed[from, "timestamp"])
  distCum <- sapply(1:n, function(x) sum(observed[from[x]:to[x], "stepDistanceObserved"], 
                                         na.rm = TRUE)) 
  distDir <- sqrt((observed[to, "xObserved"] - observed[from, "xObserved"])^2 + 
                    (observed[to, "yObserved"] - observed[from, "yObserved"])^2)
  
  stepFm <- lm(deltaT ~ distCum + I(distCum^2) + distDir)
  ss4 <- stepFm$coefficients
  
  # points less than 30 min apart
  timeDiff <- diff(observed[, "timestamp"])
  shortSteps <- which(timeDiff < 30) + 1
  
  # distances traveled in short time intervals
  shortDist <- observed[shortSteps, "stepDistanceObserved"]
  ss5 <- sd(shortDist, na.rm = TRUE)
  ss6 <- mean(shortDist, na.rm = TRUE)
  ss7 <- quantile(shortDist,  probs = c(0.1,0.9), na.rm = TRUE) # ,0.2,0.8
  
  # angles during short time intervals
  shortAngle <- observed[shortSteps, "turningAngleObserved"]
  ss8 <- sd(shortAngle, na.rm = TRUE)
  ss9 <- mean(shortAngle, na.rm = TRUE)
  ss10 <- quantile(shortAngle,  probs = c(0.1,0.2,0.8,0.9), na.rm = TRUE)
  
  # start-end distance
  ss11 <- sqrt((observed[NROW(observed), "xObserved"] - observed[1, "xObserved"])^2 + 
                 (observed[NROW(observed), "yObserved"] - observed[1, "yObserved"])^2)
  
  # max. effort
  qfm <- quantreg::rq(distCum ~ deltaT, tau = c(0.1, 0.9), model = FALSE) # , 0.2, 0.8
  ss12 <- qfm$coefficients
  
  # distCum moments
  ss13 <- quantile(distCum, probs = c(0.1,0.9), na.rm = TRUE) # ,0.2,0.8
  ss14 <- quantile(distDir, probs = c(0.1,0.9), na.rm = TRUE) # ,0.2,0.8
  ss15 <- moments::skewness(distCum, na.rm = TRUE)
  ss16 <- moments::skewness(distDir, na.rm = TRUE)
  
  # extent of movement for max effort
  ss17 <- (max(observed[, "xObserved"]) - min(observed[, "xObserved"])) *
    (max(observed[, "yObserved"]) - min(observed[, "yObserved"]))
  # cumulated step distance
  ss18 <- sum(observed[, "stepDistanceObserved"], na.rm = TRUE)
  ss19 <- sum(distCum, na.rm = TRUE)
  ss20 <- sum(distDir, na.rm = TRUE)
  # convex hull area
  chpts <- grDevices::chull(observed[,c("xObserved","yObserved")])
  chpts <- c(chpts, chpts[1])
  ss21 <- sp::Polygon(observed[chpts ,c("xObserved","yObserved")], hole = F)@area
  
  # cluster in distance matrix
  dd = sp::spDists(observed[,c("xObserved","yObserved")])
  hdd <- hist(as.vector(dd), breaks = nrow(observed)/2, plot = FALSE)
  ss22 <- sd(hdd$counts)
  ss23 <- sum(diff(sign(diff(hdd$counts)))==-2)
  
  # 
  qfm <- quantreg::rq(distDir ~ deltaT, tau = c(0.1,0.8,0.9), model = FALSE) # , 0.2
  ss24 <- qfm$coefficients
  
  # skewness of hs
  ss25 <- moments::skewness(observed[, "habitatSuitabilityObserved"], na.rm = TRUE)
  
  # resource selection function
  ss26 <- mean(observed[, "rsc"], na.rm = TRUE)
  ss27 <- sd(observed[, "rsc"], na.rm = TRUE)
  ss28 <- quantile(observed[, "rsc"], probs=c(.1,.2,.8,.9),na.rm = TRUE)
  
  out <- unlist(sapply(seq.int(28), function(x) get(paste0("ss",x))))
  names(out) <- NULL
  return(out)
  
}
