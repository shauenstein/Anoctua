# SH 07-09-2016
#
# IBM calibration
#
#-#-#-#-#-#-#-#-#-#-#

##########
# Setup  #
#--------##--------##--------##--------##--------##--------##--------##--------#
#set seed
set.seed(01022016)

# load packages
library(sp)
library(raster)
library(ranger)
library(LOSim)
library(abctools)
library(BayesianTools) 
library(lubridate)


# source functions
source("../code/CalibDataPrep.r")
source("../code/sumStats.r")
source("../code/getSumStats.r")
source("../code/createPriors.r")
source("../code/runABC.r")
source("../code/ABCsummaryChoice.r")
source("../code/plotSummariesCorrelation.r")
source("../code/ABCgetEstimate.r")
source("../code/mergeLOSim.r")

##############################
# load/prepare model input   #
#--------##--------##--------##--------##--------##--------##--------##--------#
# load telemetry data
cleaned <- read.csv("../data/natal_dispersal_cleaned.csv",
                    stringsAsFactors = FALSE)
# date as POSIXct
cleaned$date <- as.POSIXct(cleaned$date)
cleaned$date <- lubridate::with_tz(cleaned$date, "GMT")
# load IDs of training and validation set
load("../data/trainTestIDs.RData")


#----------------------------------------------------#

# classify dispersal locations # date of first outing (400m buffer)
firstOut <- 1:3
for(i in seq_along(unique(cleaned$id))){
  idIndex <- which(cleaned$id == unique(cleaned$id)[i])
  birthIndex <- idIndex[which.min(cleaned$date[cleaned$id == 
                                                 unique(cleaned$id)[i]])]
  distToBL <- sapply(idIndex, function(x) 
    sqrt((cleaned$posX[x] - cleaned$posX[birthIndex])^2
         + (cleaned$posY[x] - cleaned$posY[birthIndex])^2))
  
  search = TRUE
  j = 1
  while(search){
    if(distToBL[j] > 400){
      search = FALSE
      firstOut[i] = idIndex[j]
    }
    j = j + 1
  }
} 

firstOut <- cleaned[firstOut,]


# select locations within dispersal period, i.e. 
dispersal <- cleaned[-c(1:NROW(cleaned)),]
for(x in unique(cleaned$id)){
  subs <-  cleaned[cleaned$id == x 
                   & cleaned$date >= firstOut$date[firstOut$id == x] 
                   & cleaned$date <= (firstOut$date[firstOut$id == x] + 3600*24*60),]
  dispersal <- rbind(dispersal, subs)
}
rm(search, x, j, i, idIndex, distToBL, birthIndex, subs)

# # spatialpointsdataframe
dispSpdf <- SpatialPointsDataFrame(coords = dispersal[,c("posX","posY")], 
                                   data = dispersal) 
# add correct projection: DHDN / Gauss-Kruger zone 3
proj4string(dispSpdf) <- CRS("+init=epsg:31467")

#----------------------------------------------------#
# get start points from dispersal data, i.e. first recorded locations in dispersal period
startIndex <- 1:3
for(i in seq_along(unique(dispSpdf$id))){
  idIndex <- which(dispSpdf$id == unique(dispSpdf$id)[i])
  startIndex[i] <- idIndex[which(dispSpdf$date[idIndex] == min(dispSpdf$date[idIndex]))]
} 
startPoints <- dispSpdf[startIndex, ]

# telemetry data extent + 50 km buffer
calibrationExtent <- extent(startPoints) + rep(c(-40000,40000), 2)

#----------------------------------------------------#
# load habitat suitability
habitatSuitability <- raster("../envVars/envOnStudyExtent/habitatSuitability_scaled_20.tif")
# habitatSuitability <- raster("../envVars/JF/S123_bins.tif")
# proj4string(habitatSuitability) <- CRS("+init=epsg:31467")
# crop habitat suitability map
habSuitCalib <- crop(habitatSuitability, calibrationExtent)
habSuitCalib[is.na(values(habSuitCalib))] <- 0 # set NAs to zero <- reflection
HSvals <- as.matrix(habSuitCalib)

#----------------------------------------------------#
# priors for observation error
Accuracy <- sapply(unique(dispSpdf$id), function(x) dispersal$accuracy[dispersal$id == x])
accLookUp <- c(0,5,20,50,100,200)
AccMetric <- sapply(unique(dispSpdf$id), function(x) accLookUp[Accuracy[[x]]+1])
AccMetric <- unlist(AccMetric)

#----------------------------------------------------#
# prepare telemetry data
rscRange = 200
calibrationData <- calibTransform(x = dispersal$posX, y = dispersal$posY,
                                  date = dispersal$date, id = dispersal$id,
                                  habitatSuitability = habSuitCalib, rscExtent = rscRange)

# filter individuals with more than 50 recorded locations
nRelocs <- sapply(seq_along(calibrationData), function(x) NROW(calibrationData[[x]]$data))
# calibID <- names(calibrationData)[nRelocs  > 50] # & names(calibrationData) %in% trainingID]
calibID <- names(calibrationData)[nRelocs  > 50]

#save(calibID,validID, file = "../data/trainTestIDs_g50rl.RData")

####################
# parameter sample #
#--------##--------##--------##--------##--------##--------##--------##--------#
# number of simulations per individual
n = 500000

# parameter names
pnames <- c("habitatPreference", 
            "stepShape", 
            "stepScale", 
            "directionalBias", 
            "dispRestMean", 
            "dispRestSd", 
            "maximumEffort", 
            "roostLambda",
            "observationError",
            "rscRange")

# prior distributions
distr <- list(runif, runif, NULL, runif, NULL, NULL, NULL, runif, rnorm, NULL)

# prior ranges
ranges <- list(c(0.1, 5), 
               c(0.1, 4), 
               2, 
               c(0.1, 3.5), 
               1.5,  
               0.7, 
               6000,
               c(.01,30),
               c(mean(AccMetric),sd(AccMetric)),
               rscRange
)

# which parameters are fixed
fixed = c(3,5,6,7,10)

# create parameters
Ps <- createPriors(pnames, distr, ranges, n, redraw = 9, fixed)

# careful: throw away first x in y. run, throw away first 30 000 (r.seed issue)
# Ps <- Ps[30001:n,]
####################
# model setup      #
#--------##--------##--------##--------##--------##--------##--------##--------#

# environment input
HS <- list(values = HSvals, upperleft = extent(habSuitCalib)[c(1,4)], resolution = raster::res(habSuitCalib)[1])

save(Ps, HS, calibrationData, calibID, sumStats, getSumStats, runABC, file = "../data/input.RData")
#--------##--------##--------##--------##--------##--------##--------##--------#

simulation <- runABC(Ps, HS, calibrationData[calibID], parallel = TRUE, 
                      maxIterations = 40000, randomSeed = NULL, 
                      sumStatsSampleSize = 5000, sumStats, getSumStats)


# save(simulation, file = "../data/simulation.RData")
# load("../data/simulation.RData")
# load("../data/simulation01.RData")
# load("../data/simulation02.RData")


# merge simulations
# simulation <- mergeLOSim(simulation01, simulation02)
# save(simulation, file = "../data/simulation.RData")
# load("../data/simulation.RData")

# selected summary estimation
summaries <- abcCreateSummaries(simulation, method = "rF", summarySelection = NULL, 
                                  targetParameters = c(1,2,4,8), getSumStats = getSumStats, 
                                  sumStats = sumStats, sumStatsSampleSize = 5000, parallel = TRUE,
                                  predictionPredictor = TRUE)
# save(summaries, file = "../data/summaries.RData")
# load("../data/summaries.RData")
# load("../data/summaries_valid.RData")

# rejection
posterior <- getEstimate(summaries, proportionFiltered = 1/500, parallel = FALSE, regr.adj = TRUE, comp.MAP = FALSE)
posterior_valid <- getEstimate(summaries_valid, proportionFiltered = 1/500, parallel = FALSE, regr.adj = TRUE, comp.MAP = FALSE)
# merge to one posterior object
for(i in names(posterior)[1:7]){
  posterior[[i]] <- append(posterior[[i]],posterior_valid[[i]])
}
# save(posterior, file = "../data/posterior.RData")
# load("../data/posterior.RData")


# plot posterior distributions
addInfo_all <- read.csv("../data/Additional_data_juveniles.csv")

# load calibration IDs
load("../data/trainTestIDs_g50rl.RData") # reload original calib, valid indices, to organise add. info in the same order
# use only calibID information
addInfo <- addInfo_all[addInfo_all$Name1 %in% calibID,]
addInfo <- rbind(addInfo, addInfo_all[addInfo_all$Name1 %in% validID,])
# pool by Sex

# post-sampling regression adjustment or not
regrAdj <- "parameters.regrAdj"

# habitat preference
fHP <- unlist(lapply(which(addInfo$Sex == "f"), function(x) posterior[[regrAdj]][[x]][,1])); quantile(fHP, probs = c(0.1,0.5,0.9))
mHP <- unlist(lapply(which(addInfo$Sex == "m"), function(x) posterior[[regrAdj]][[x]][,1])); quantile(mHP, probs = c(0.1,0.5,0.9))
# t.test(fHP,mHP)

# step shape
fSS <- unlist(lapply(which(addInfo$Sex == "f"), function(x) posterior[[regrAdj]][[x]][,2])); quantile(fSS, probs = c(0.1,0.5,0.9))
mSS <- unlist(lapply(which(addInfo$Sex == "m"), function(x) posterior[[regrAdj]][[x]][,2])); quantile(mSS, probs = c(0.1,0.5,0.9))
# seq(0,10, len = 1000000)[which.max(dgamma(seq(0,10, len = 1000000), shape = quantile(fSS, probs = c(0.1,0.5,0.9))[3], scale = 2))]*20
# seq(0,10, len = 1000000)[which.max(dgamma(seq(0,10, len = 1000000), shape = quantile(mSS, probs = c(0.1,0.5,0.9))[3], scale = 2))]*20
# t.test(fSS,mSS)

# directional bias
fDB <- unlist(lapply(which(addInfo$Sex == "f"), function(x) posterior[[regrAdj]][[x]][,3])); quantile(fDB, probs = c(0.1,0.5,0.9))
mDB <- unlist(lapply(which(addInfo$Sex == "m"), function(x) posterior[[regrAdj]][[x]][,3])); quantile(mDB, probs = c(0.1,0.5,0.9))
# t.test(fDB,mDB)

# maximum effort
fME <- unlist(lapply(which(addInfo$Sex == "f"), function(x) posterior[[regrAdj]][[x]][,4])); quantile(fME, probs = c(0.1,0.5,0.9))
mME <- unlist(lapply(which(addInfo$Sex == "m"), function(x) posterior[[regrAdj]][[x]][,4])); quantile(mME, probs = c(0.1,0.5,0.9))
# t.test(fME,mME)



# plot posterior distributions colour coded by sex
pNames <- c("habitatPreference", "stepShape", "directionalBias", "roostLambda")
for(pInt in pNames){
  ggdf <- data.frame(xval = unlist(lapply(seq_along(posterior[[regrAdj]]), function(x) density(posterior[[regrAdj]][[x]][,pInt], n = 500)$x)),
                     dens = unlist(lapply(seq_along(posterior[[regrAdj]]), function(x) density(posterior[[regrAdj]][[x]][,pInt], n = 500)$y)),
                     id = rep(seq_along(posterior[[regrAdj]]), each = 500),
                     sex = rep(ifelse(is.na(addInfo$Sex),
                                      rgb(0,0,0,0.4),
                                      ifelse(addInfo$Sex == "f",
                                             rgb(0.55,0,0,0.4),
                                             rgb(0,0,0.55,0.4))), each = 500)
  )
  ggdf_f <- data.frame(xval = unlist(lapply(which(addInfo$Sex == "f"), function(x) posterior[[regrAdj]][[x]][,pInt])))
  ggdf_m <- data.frame(xval = unlist(lapply(which(addInfo$Sex == "m"), function(x) posterior[[regrAdj]][[x]][,pInt])))
  
  ggplot() +
    geom_line(data = ggdf, aes(x=xval,y=dens,group=id,colour=sex)) +
    scale_colour_identity() +
    stat_density(data = ggdf_f, aes(xval), colour = rgb(0.55,0,0), size = 1.5, geom = "line") +
    stat_density(data = ggdf_m, aes(xval), colour = rgb(0,0,0.55), size = 1.5, geom = "line") +
    xlab("") + 
    ylab("") +
    theme(
      panel.border = element_rect(fill = NA),
      axis.text.x = element_text(family = "Linux Libertine", size = 60, margin = margin(t = 10, r = 0, b = 0, l = 0)),
      axis.text.y = element_blank(),
      axis.ticks.length = unit(.25, "cm"),
      axis.ticks.y = element_blank(),
      legend.background = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  ggsave(paste0("posterior_",pInt,".pdf"),
         device = cairo_pdf,
         path = "../figures/results/posterior_indiv_plots",
         width = 25, height = 25, units = "cm", dpi = 400)
}


# plot posterior distributions colour coded by fed/not-fed
for(pInt in pNames){
  ggdf <- data.frame(xval = unlist(lapply(seq_along(posterior[[regrAdj]]), function(x) density(posterior[[regrAdj]][[x]][,pInt], n = 500)$x)),
                     dens = unlist(lapply(seq_along(posterior[[regrAdj]]), function(x) density(posterior[[regrAdj]][[x]][,pInt], n = 500)$y)),
                     id = rep(seq_along(posterior[[regrAdj]]), each = 500),
                     fed = rep(ifelse(is.na(addInfo$fed),
                                      rgb(0,0,0,0.4),
                                      ifelse(addInfo$fed == 0,
                                             rgb(0.8,0,0,0.4),
                                             rgb(0,.6,0,0.4))), each = 500)
  )
  ggdf_nfed <- data.frame(xval = unlist(lapply(which(addInfo$fed == 0), function(x) posterior[[regrAdj]][[x]][,pInt])))
  ggdf_fed <- data.frame(xval = unlist(lapply(which(addInfo$fed == 1), function(x) posterior[[regrAdj]][[x]][,pInt])))
  
  ggplot() +
    geom_line(data = ggdf, aes(x=xval,y=dens,group=id,colour=fed)) +
    scale_colour_identity() +
    stat_density(data = ggdf_nfed, aes(xval), colour = rgb(0.8,0,0), size = 1.5, geom = "line") +
    stat_density(data = ggdf_fed, aes(xval), colour = rgb(0,0.6,0), size = 1.5, geom = "line") +
    xlab("") + 
    ylab("") +
    theme(
      panel.border = element_rect(fill = NA),
      axis.text.x = element_text(family = "Linux Libertine", size = 60, margin = margin(t = 10, r = 0, b = 0, l = 0)),
      axis.text.y = element_blank(),
      axis.ticks.length = unit(.25, "cm"),
      axis.ticks.y = element_blank(),
      legend.background = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  ggsave(paste0("posterior_",pInt,"_feeding.pdf"),
         device = cairo_pdf,
         path = "../figures/results/posterior_indiv_plots",
         width = 25, height = 25, units = "cm", dpi = 400)
}
# pdf("../figures/posterior.pdf", width = 10, height = 10)
# par(mfrow=c(2,2), mar=c(1.5,.5,.5,0), oma=c(0,0,0,1), cex.lab = 1.5, cex.axis = 1.3, las = 1, mgp = c(1.5,0.4,0), tcl = -0.2, pty = "s")
# plot(density(posterior[[regrAdj]][[1]][,1]), xlim=c(0.1,5), 
#      main = "", xlab = "", ylab = "", 
#      ylim =c(0,2.8), type = "n", yaxt = "n",lwd=0.5)
# for(i in 1:length(posterior[[regrAdj]])) 
#   lines(density(posterior[[regrAdj]][[i]][,1]), col = ifelse(is.na(addInfo$Sex[i]),
#                                                                 rgb(0.7,0.7,0.7),
#                                                                 ifelse(addInfo$Sex[i] == "f",
#                                                                        rgb(0.55,0,0,0.4),
#                                                                        rgb(0,0,0.55,0.4))),
#         lwd=0.5)
# # pooled
# lines(density(fHP), col = rgb(0.55,0,0), lwd=2) 
# lines(density(mHP), col = rgb(0,0,0.55), lwd=2)
# 
# plot(density(posterior[[regrAdj]][[1]][,2]), xlim=c(0.1,4), 
#      main = "", xlab = "", ylab = "", 
#      ylim =c(0,3.6), type = "n", yaxt = "n",lwd=0.5)
# for(i in 1:length(posterior[[regrAdj]])) 
#   lines(density(posterior[[regrAdj]][[i]][,2]), col = ifelse(is.na(addInfo$Sex[i]),
#                                                                 rgb(0,0,0,0.4),
#                                                                 ifelse(addInfo$Sex[i] == "f",
#                                                                        rgb(0.55,0,0,0.4),
#                                                                        rgb(0,0,0.55,0.4))),
#         lwd=0.5)
# # pooled
# lines(density(fSS), col = rgb(0.55,0,0), lwd=2) 
# lines(density(mSS), col = rgb(0,0,0.55), lwd=2)
# 
# plot(density(posterior[[regrAdj]][[1]][,3]), xlim=c(0.1,3.5), 
#      main = "", xlab = "", ylab = "",
#      ylim =c(0,4.1), type = "n", yaxt = "n",lwd=0.5)
# for(i in 1:length(posterior[[regrAdj]])) 
#   lines(density(posterior[[regrAdj]][[i]][,3]), col = ifelse(is.na(addInfo$Sex[i]),
#                                                                 rgb(0,0,0,0.4),
#                                                                 ifelse(addInfo$Sex[i] == "f",
#                                                                        rgb(0.55,0,0,0.4),
#                                                                        rgb(0,0,0.55,0.4))),
#         lwd=0.5)
# # pooled
# lines(density(fDB), col = rgb(0.55,0,0), lwd=2) 
# lines(density(mDB), col = rgb(0,0,0.55), lwd=2)
# 
# plot(density(posterior[[regrAdj]][[1]][,4]), xlim=c(0.01,30), 
#      main = "", xlab = "", ylab = "",
#      ylim =c(0,0.9), type = "n", yaxt = "n",lwd=0.5)
# for(i in 1:length(posterior[[regrAdj]])) 
#   lines(density(posterior[[regrAdj]][[i]][,4]), col = ifelse(is.na(addInfo$Sex[i]),
#                                                                 rgb(0,0,0,0.4),
#                                                                 ifelse(addInfo$Sex[i] == "f",
#                                                                        rgb(0.55,0,0,0.4),
#                                                                        rgb(0,0,0.55,0.4))),
#         lwd=0.5)
# # pooled
# lines(density(fME), col = rgb(0.55,0,0), lwd=2) 
# lines(density(mME), col = rgb(0,0,0.55), lwd=2)
# 
# dev.off()


#--------##--------##--------##--------##--------##--------##--------##--------#
x11(10,10)
correlationPlot(posterior[[regrAdj]][[69]], density = "smooth")
png("../figures/results/summaries_rf_corPlot.png", width = 1000, height = 1000, units = "px")
plotSC(summaries, plot.which = 1, cex.axis = 2)
dev.off()

load("../data/summariesLr.RData")
png("../figures/results/summaries_lr_corPlot.png", width = 1000, height = 1000, units = "px")
plotSC(summariesLr, plot.which = 1, cex.axis = 2)
dev.off()

#--------##--------##--------##--------##--------##--------##--------##--------#

# plot posterior median against weight
pdf("../figures/results/appendix/posterior_weight.pdf", width = 10.5, height = 10)
par(mfrow=c(2,2), mar=c(1.5,.5,.5,0), oma=c(1,1.5,0,.1), cex.lab = 2, cex.axis = 2, las = 1, mgp = c(1.5,0.6,0), tcl = 0.2, pty = "s")
plot(do.call(rbind,posterior$median)[,"habitatPreference"]~addInfo$weight,
     main = "", xlab = "", ylab = "", pch = "*",
     ylim =c(0,5),lwd=0.5)
plot(do.call(rbind,posterior$median)[,"stepShape"]~addInfo$weight,
     main = "", xlab = "", ylab = "", pch = "*",
     ylim =c(0,4),lwd=0.5)
plot(do.call(rbind,posterior$median)[,"directionalBias"]~addInfo$weight,
     main = "", xlab = "", ylab = "", pch = "*",
     ylim =c(0,3.5),lwd=0.5)
plot(do.call(rbind,posterior$median)[,"roostLambda"]~addInfo$weight,
     main = "", xlab = "", ylab = "", pch = "*",
     ylim =c(0,30),lwd=0.5)
dev.off()

# against tarsus
pdf("../figures/results/appendix/posterior_tarsus.pdf", width = 10.5, height = 10)
par(mfrow=c(2,2), mar=c(1.5,.5,.5,0), oma=c(1,1.5,0,.1), cex.lab = 2, cex.axis = 2, las = 1, mgp = c(1.5,0.6,0), tcl = 0.2, pty = "s")
plot(do.call(rbind,posterior$median)[,"habitatPreference"]~addInfo$tarsus,
     main = "", xlab = "", ylab = "", pch = "*",
     ylim =c(0,5),lwd=0.5)
plot(do.call(rbind,posterior$median)[,"stepShape"]~addInfo$tarsus,
     main = "", xlab = "", ylab = "", pch = "*",
     ylim =c(0,4),lwd=0.5)
plot(do.call(rbind,posterior$median)[,"directionalBias"]~addInfo$tarsus,
     main = "", xlab = "", ylab = "", pch = "*",
     ylim =c(0,3.5),lwd=0.5)
plot(do.call(rbind,posterior$median)[,"roostLambda"]~addInfo$tarsus,
     main = "", xlab = "", ylab = "", pch = "*",
     ylim =c(0,30),lwd=0.5)
dev.off()

# against wing
pdf("../figures/results/appendix/posterior_wing.pdf", width = 10.5, height = 10)
par(mfrow=c(2,2), mar=c(1.5,.5,.5,0), oma=c(1,1.5,0,.1), cex.lab = 2, cex.axis = 2, las = 1, mgp = c(1.5,0.6,0), tcl = 0.2, pty = "s")
plot(do.call(rbind,posterior$median)[,"habitatPreference"]~addInfo$wing,
     main = "", xlab = "", ylab = "", pch = "*",
     ylim =c(0,5),lwd=0.5)
plot(do.call(rbind,posterior$median)[,"stepShape"]~addInfo$wing,
     main = "", xlab = "", ylab = "", pch = "*",
     ylim =c(0,4),lwd=0.5)
plot(do.call(rbind,posterior$median)[,"directionalBias"]~addInfo$wing,
     main = "", xlab = "", ylab = "", pch = "*",
     ylim =c(0,3.5),lwd=0.5)
plot(do.call(rbind,posterior$median)[,"roostLambda"]~addInfo$wing,
     main = "", xlab = "", ylab = "", pch = "*",
     ylim =c(0,30),lwd=0.5)
dev.off()

# against fed/non-fed
pdf("../figures/results/appendix/posterior_fed.pdf", width = 10.5, height = 10)
par(mfrow=c(2,2), mar=c(1.5,.5,.5,0), oma=c(1,1.5,0,.1), cex.lab = 2, cex.axis = 2, las = 1, mgp = c(1.5,0.6,0), tcl = 0.2, pty = "s")
plot(do.call(rbind,posterior$median)[,"habitatPreference"]~factor(addInfo$fed, levels = c(0,1), labels = c("not fed", "fed")),
     main = "", xlab = "", ylab = "", ylim =c(0,5),lwd=0.5)
plot(do.call(rbind,posterior$median)[,"stepShape"]~factor(addInfo$fed, levels = c(0,1), labels = c("not fed", "fed")),
     main = "", xlab = "", ylab = "", ylim =c(0,4),lwd=0.5)
plot(do.call(rbind,posterior$median)[,"directionalBias"]~factor(addInfo$fed, levels = c(0,1), labels = c("not fed", "fed")),
     main = "", xlab = "", ylab = "", ylim =c(0,3.5),lwd=0.5)
plot(do.call(rbind,posterior$median)[,"roostLambda"]~factor(addInfo$fed, levels = c(0,1), labels = c("not fed", "fed")),
     main = "", xlab = "", ylab = "", ylim =c(0,30),lwd=0.5)
dev.off()

# against sex
pdf("../figures/results/appendix/posterior_sex.pdf", width = 10.5, height = 10)
par(mfrow=c(2,2), mar=c(1.5,.5,.5,0), oma=c(1,1.5,0,.1), cex.lab = 2, cex.axis = 2, las = 1, mgp = c(1.5,0.6,0), tcl = 0.2, pty = "s")
plot(do.call(rbind,posterior$median)[,"habitatPreference"]~factor(addInfo$Sex, levels = c("f","m"), labels = c("female", "male")),
     main = "", xlab = "", ylab = "", ylim =c(0,5),lwd=0.5)
plot(do.call(rbind,posterior$median)[,"stepShape"]~factor(addInfo$Sex, levels = c("f","m"), labels = c("female", "male")),
     main = "", xlab = "", ylab = "", ylim =c(0,4),lwd=0.5)
plot(do.call(rbind,posterior$median)[,"directionalBias"]~factor(addInfo$Sex, levels = c("f","m"), labels = c("female", "male")),
     main = "", xlab = "", ylab = "", ylim =c(0,3.5),lwd=0.5)
plot(do.call(rbind,posterior$median)[,"roostLambda"]~factor(addInfo$Sex, levels = c("f","m"), labels = c("female", "male")),
     main = "", xlab = "", ylab = "", ylim =c(0,30),lwd=0.5)
dev.off()

#--------##--------##--------##--------##--------##--------##--------##--------#
# estimate summary statistics using partial least squre regression

# load("../data/simulation.RData")
# 
# # selected summary estimation
# plsr_ssEstimates <- abcCreateSummaries(simulation, method = "saABC", summarySelection = NULL, 
#                                        targetParameters = c(1,2,4,8), getSumStats = getSumStats, 
#                                        sumStats = sumStats, sumStatsSampleSize = 5000, nCPU = NULL,
#                                        predictionPredictor = FALSE)

# save(plsr_ssEstimates, file = "../data/plsrEstimates.RData")

#--------##--------##--------##--------##--------##--------##--------##--------#
#--------##--------##--------##--------##--------##--------##--------##--------#
#--------##--------##--------##--------##--------##--------##--------##--------#
#--------##--------##--------##--------##--------##--------##--------##--------#
#--------##--------##--------##--------##--------##--------##--------##--------#
#--------##--------##--------##--------##--------##--------##--------##--------#
#--------##--------##--------##--------##--------##--------##--------##--------#
#--------##--------##--------##--------##--------##--------##--------##--------#
#--------##--------##--------##--------##--------##--------##--------##--------#
#--------##--------##--------##--------##--------##--------##--------##--------#