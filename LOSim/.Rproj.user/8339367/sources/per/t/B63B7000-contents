# SH 05-10-2016
#
# forward simulation script
#
#-#-#-#-#-#-#-#-#-#-#

##########
# Setup  #
#--------##--------##--------##--------##--------##--------##--------##--------#
#set seed
set.seed(06102016); 


# load packages
library(sp)
library(rgdal)
library(rgeos)
library(raster)
library(LOSim)
library(spatstat)
library(ggplot2)


# source functions
source("../code/forwardSimulation.r")
source("../code/plotLOSim.r")


#--------##--------##--------##--------##--------##--------##--------##--------#

# load posterior
load("../data/posterior.RData")

# load habitat suitability
habitatSuitability <- raster("../envVars/envOnStudyExtent/habitatSuitability_scaled_20.tif")
habitatSuitability[is.na(values(habitatSuitability))] <- 0 # set NAs to zero <- reflection

HSvals <- as.matrix(habitatSuitability)

HS <- list(values = HSvals, upperleft = extent(habitatSuitability)[c(1,4)], resolution = raster::res(habitatSuitability)[1])


#--------##--------##--------##--------##--------##--------##--------##--------#

# load ring data to provide start locations
ringRaw <- read.csv("../data/Steinkauz_Beringungsdaten_ab_2005.csv")
ringData <- data.frame(id = as.factor(ringRaw$STRRINGNR),
                       date = as.POSIXct(strptime(ringRaw$DTMDATE, "%d/%m/%Y %H:%M", 
                                                  tz = "Europe/Berlin")),
                       longitude = ringRaw$LNGLONG,
                       latitude = ringRaw$LNGLAT,
                       geoReg = ringRaw$STRTEXT.2) 

# as spatial data frame
ringSDF <- SpatialPointsDataFrame(coords = ringData[,c("longitude","latitude")], 
                                  data = ringData) 

# add correct projection: wgs84
proj4string(ringSDF) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

## Load administrative boundary BW (http://www.gadm.org)
bw <- readOGR("../envVars/administrativeBoundary", layer = "DEU_adm1",
              encoding = "ESRI Shapefile", verbose = FALSE)
bw <- bw[bw$NAME_1 == "Baden-Württemberg",]

## crop ringing data to bw
ringSDF <- crop(ringSDF, extent(bw))
ringBW <- gIntersects(ringSDF, bw, byid = TRUE)
ringBW <- ringSDF[c(ringBW),]


# reproject to GK3
ringBW <- spTransform(ringBW, proj4string(habitatSuitability))

# select locations from 2012 to 2015
startLocations <- as.data.frame(coordinates(ringBW)[format(ringBW$date, format="%Y") %in% 2011:2015,])
colnames(startLocations) <- c("x","y")
# assign population id
startLocations$pop.id <- ifelse(startLocations$y < 5300000, 1, 
                                ifelse(startLocations$y < 5400000 & startLocations$x < 3440000, 2, 
                                       3))


# target 5000 individuals for each population
for(pop in 1:3){
  slSample <- sample(which(startLocations$pop.id == pop), size = 5000, replace = TRUE)
  startLocs <- data.frame(x = as.integer(startLocations$x[slSample]), y = as.integer(startLocations$y[slSample]))
  
  save(HS, posterior, startLocs, pop, slSample, file = paste0("../data/input_forward5000pop",pop,".RData"))
}


#--------##--------##--------##--------##--------##--------##--------##--------#
# add option to define simulation period in LOSim::runSimulation

# forward simulations
# sim55 <- LOSim(posteriorObj = posterior,
#                habitatSuitability = HS,
#                individuals = NULL,
#                startLocations = startLocs,
#                fixParameters = c(2,1.5,0.7,6000),
#                maxPeriod = 60,
#                maxIterations = 40000,
#                generations = 55,
#                randomSeed = NULL,
#                reflection = 30,
#                reflectionValues = 0,
#                parallel = T,
#                clusterType = "SOCK")


# save(sim55, file = "../data/prediction.RData")
# load("../data/prediction.RData")
# plot results
# x11(8,10)
# plot(habitatSuitability, xlab = "", ylab = "")
# plotLOSim(sim55, plotHS = "onTop", col = rgb(0,0,0,0.1))



#--------##--------##--------##--------##--------##--------##--------##--------#
# order by population
# pop1 <- lapply(seq_along(sim55$simulations), 
#                function(x) sim55$simulations[[x]][which(sim55$startLocactions$pop.id == 1)])
# pop2 <- lapply(seq_along(sim55$simulations), 
#                function(x) sim55$simulations[[x]][which(sim55$startLocactions$pop.id == 2)])
# pop3 <- lapply(seq_along(sim55$simulations), 
#                function(x) sim55$simulations[[x]][which(sim55$startLocactions$pop.id == 3)])
# 
# for(g in seq_along(pop1)){
#   pop1[[g]] <- do.call(rbind, pop1[[g]])
#   # pop2[[g]] <- do.call(rbind, pop2[[g]])
#   # pop3[[g]] <- do.call(rbind, pop3[[g]])
# }
# # load("../data/preds_pop1.RData")
# # load("../data/preds_pop2.RData")
# # load("../data/preds_pop3.RData")
# 
# rast <- raster(ext = extent(habitatSuitability), resolution = c(1000,1000), 
#                crs = proj4string(habitatSuitability))
# 
# p1_1 <- rasterize(pop1[[1]][,c("x","y")], rast, fun = "count")
# p1_7 <- rasterize(do.call(rbind, pop1[1:7])[,c("x","y")], rast, fun = "count")
# p1_20 <- rasterize(do.call(rbind, pop1[1:20])[,c("x","y")], rast, fun = "count")
# p1_55 <- rasterize(do.call(rbind, pop1)[,c("x","y")], rast, fun = "count")
# par(mfrow = c(2,2))
# plot(p1_1); plot(p1_7); plot(p1_20); plot(p1_55);
# par(mfrow = c(1,1))

#--------##--------##--------##--------##--------##--------##--------##--------#
# plot dispersal map, coloured by population for generation 1, 3, 7, 20, 55

rast <- raster(xmn=3308355, xmx=3610195, ymn=5182816, ymx=5517416, 
               crs = CRS("+proj=tmerc +lat_0=0 +lon_0=9 +k=1 +x_0=3500000 +y_0=0 +ellps=bessel +towgs84=598.1,73.7,418.2,0.202,0.045,-2.455,6.7 +units=m +no_defs"), resolution = 500, vals=NA)
load("../data/preds5000_pop1.RData");load("../data/preds5000_pop2.RData");load("../data/preds5000_pop3.RData");

for(g in seq_along(sim5000pop1$simulations)){
  sim5000pop1$simulations[[g]] <- do.call(rbind, sim5000pop1$simulations[[g]])
  sim5000pop2$simulations[[g]] <- do.call(rbind, sim5000pop2$simulations[[g]])
  sim5000pop3$simulations[[g]] <- do.call(rbind, sim5000pop3$simulations[[g]])
}

# pop1_rast55 <- rasterize(do.call(rbind, sim5000pop1$simulations[1:55])[,c("x","y")], rast, fun = "count")
# pop2_rast55 <- rasterize(do.call(rbind, sim5000pop2$simulations[1:55])[,c("x","y")], rast, fun = "count")
# pop3_rast55 <- rasterize(do.call(rbind, sim5000pop3$simulations[1:55])[,c("x","y")], rast, fun = "count")
# writeRaster(pop1_rast55, filename = "../data/pop1_predictions.tif", format = "GTiff")
# writeRaster(pop2_rast55, filename = "../data/pop2_predictions.tif", format = "GTiff")
# writeRaster(pop3_rast55, filename = "../data/pop3_predictions.tif", format = "GTiff")

maps.pop1 <- maps.pop2 <- maps.pop3 <- list()
index = 1
for(w in c(1,3,7,20,55)){
  maps.pop1[[index]] <- rasterize(do.call(rbind, sim5000pop1$simulations[1:w])[,c("x","y")], rast, fun = "count")
  maps.pop2[[index]] <- rasterize(do.call(rbind, sim5000pop2$simulations[1:w])[,c("x","y")], rast, fun = "count")
  maps.pop3[[index]] <- rasterize(do.call(rbind, sim5000pop3$simulations[1:w])[,c("x","y")], rast, fun = "count")
  index = index + 1
}

# save(maps.pop1,maps.pop2,maps.pop3, file = "../data/densitymaps.RData")
# load("../data/densitymaps.RData")

#####----#----#----#----#----#----#----#----#----#----#----#----#####
# ggplot theme
mytheme<- theme(plot.margin = unit(c(0,0,0,0), "cm"),
                panel.spacing = unit(0,"null"),
                axis.title.x = element_text(size=26),
                axis.title.y = element_text(size=26, angle=90),
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                legend.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                legend.position = c(0.9, 0.78),
                legend.key.width = unit(0.2, "cm"),
                legend.key.height = unit(0.4, "cm"), 
                legend.title=element_text(size=15, angle = 90),
                legend.text=element_text(size=12),
                legend.direction = "vertical",
                strip.text.x = element_text(size = 30), 
                strip.background = element_blank(),
                axis.ticks = element_blank(),
                panel.background = element_rect())

#####----#----#----#----#----#----#----#----#----#----#----#----#####


# ggplot() +   
#   geom_tile(data = as.data.frame(coordinates(maps.pop1[[1]])), 
#             aes(x = x, y = y, fill = log(values(maps.pop1[[5]])))) +
#   geom_tile(data = as.data.frame(coordinates(maps.pop2[[1]])), 
#             aes(x = x, y = y, fill = log(values(maps.pop2[[5]])))) +
#   scale_fill_gradient(low=rgb(1,0,0,0.2), high="black",
#                       guide="colorbar",na.value=rgb(1,1,1,0)) +
#   coord_fixed() + 
#   xlab("")+
#   ylab("")+
#   mytheme 

for(i in seq_along(maps.pop1)){
  # logarithmise
  vals1 <- log(values(maps.pop1[[i]]))
  vals2 <- log(values(maps.pop2[[i]]))
  vals3 <- log(values(maps.pop3[[i]]))
  # range between 0 and 1
  values(maps.pop1[[i]]) <- (vals1 - min(c(vals1, vals2, vals3), na.rm = T)) / 
    (max(c(vals1, vals2, vals3), na.rm = T) - min(c(vals1, vals2, vals3), na.rm = T))
  values(maps.pop2[[i]]) <- (vals2 - min(c(vals1, vals2, vals3), na.rm = T)) / 
    (max(c(vals1, vals2, vals3), na.rm = T) - min(c(vals1, vals2, vals3), na.rm = T))
  values(maps.pop3[[i]]) <- (vals3 - min(c(vals1, vals2, vals3), na.rm = T)) / 
    (max(c(vals1, vals2, vals3), na.rm = T) - min(c(vals1, vals2, vals3), na.rm = T))
} 


# own colourramppalette with min max values
SHpalette <- function (colors, ...) 
{
  ramp <- colorRamp(colors, ...)
  function(n, minV, maxV) {
    x <- ramp(seq.int(minV, maxV, length.out = n))
    if (ncol(x) == 4L) 
      rgb(x[, 1L], x[, 2L], x[, 3L], x[, 4L], maxColorValue = 255)
    else rgb(x[, 1L], x[, 2L], x[, 3L], maxColorValue = 255)
  }
}

# population specific colour palettes
cp1 <- SHpalette(c(rgb(1,0,0,0),rgb(0,0,0,0.8)), alpha = TRUE)
cp1s <- SHpalette(c(rgb(1,0,0,0),rgb(0,0,0,1)), alpha = TRUE)
cp2 <- SHpalette(c(rgb(0,0,1,0),rgb(0,0,0,0.7)), alpha = TRUE)
cp2s <- SHpalette(c(rgb(0,0,1,0),rgb(0,0,0,1)), alpha = TRUE)
cp3 <- SHpalette(c(rgb(0,1,0,0),rgb(0,0,0,0.6)), alpha = TRUE)
cp3s <- SHpalette(c(rgb(0,1,0,0),rgb(0,0,0,1)), alpha = TRUE)

# geographic boundary of SA
bwch <- readOGR("../envVars/administrativeBoundary", layer = "BW_CH_adm0",
                encoding = "ESRI Shapefile", verbose = FALSE)
fr <- readOGR("../envVars/administrativeBoundary", layer = "FRA_adm0",
              encoding = "ESRI Shapefile", verbose = FALSE)
fr@data <- data.frame(ID = 3)
bwchfr <- rbind(bwch, fr)
# project to GK 3
bwchfr <- spTransform(bwchfr, CRS("+init=epsg:31467"))
# crop to extended SA
bwchfr <- crop(bwchfr, rast)

# overlay for specific generation time
for(g in seq_along(maps.pop1)){
  pdf(paste0("../figures/results/mapOverlayGen",round(exp(g-1)),".pdf"), width = 10 * NCOL(rast) / NROW(rast), height = 10)
  par(mar = c(0,0,0,0))
  plot(bwchfr, col = "grey80", xaxt = "n", yaxt = "n")
  plot(maps.pop1[[g]], col = cp1(100, minV = min(values(maps.pop1[[g]]), na.rm = T), 
                                 maxV = max(values(maps.pop1[[g]]), na.rm = T)), 
       add = T, xaxt = "n", yaxt = "n", legend=FALSE)
  plot(maps.pop2[[g]], col = cp2(100, minV = min(values(maps.pop2[[g]]), na.rm = T), 
                                 maxV = max(values(maps.pop2[[g]]), na.rm = T)), 
       add = T, xaxt = "n", yaxt = "n", legend=FALSE)
  plot(maps.pop3[[g]], col = cp3(100, minV = min(values(maps.pop3[[g]]), na.rm = T), 
                                 maxV = max(values(maps.pop3[[g]]), na.rm = T)), 
       add = T, xaxt = "n", yaxt = "n", legend=FALSE)
  dev.off()
  
  #only pop1
  pdf(paste0("../figures/results/mapOverlayGen",round(exp(g-1)),"_pop1.pdf"), width = 10 * NCOL(rast) / NROW(rast), height = 10)
  par(mar = c(0,0,0,0))
  plot(bwchfr, col = "grey80", xaxt = "n", yaxt = "n")
  plot(maps.pop1[[g]], col = cp1s(100, minV = min(values(maps.pop1[[g]]), na.rm = T), 
                                  maxV = max(values(maps.pop1[[g]]), na.rm = T)), 
       add = T, xaxt = "n", yaxt = "n", legend=FALSE)
  dev.off()
  
  #only pop2
  pdf(paste0("../figures/results/mapOverlayGen",round(exp(g-1)),"_pop2.pdf"), width = 10 * NCOL(rast) / NROW(rast), height = 10)
  par(mar = c(0,0,0,0))
  plot(bwchfr, col = "grey80", xaxt = "n", yaxt = "n")
  plot(maps.pop2[[g]], col = cp2s(100, minV = min(values(maps.pop2[[g]]), na.rm = T), 
                                  maxV = max(values(maps.pop2[[g]]), na.rm = T)), add = T, xaxt = "n", yaxt = "n", legend=FALSE)
  dev.off()
  
  #only pop3
  pdf(paste0("../figures/results/mapOverlayGen",round(exp(g-1)),"_pop3.pdf"), width = 10 * NCOL(rast) / NROW(rast), height = 10)
  par(mar = c(0,0,0,0))
  plot(bwchfr, col = "grey80", xaxt = "n", yaxt = "n")
  plot(maps.pop3[[g]], col = cp3s(100, minV = min(values(maps.pop3[[g]]), na.rm = T), 
                                  maxV = max(values(maps.pop3[[g]]), na.rm = T)), add = T, xaxt = "n", yaxt = "n", legend=FALSE)
  dev.off()
}



#--------##--------##--------##--------##--------##--------##--------##--------#
# count cells filled after 1 to 55 generations per population
disp.count <- data.frame(generation = seq_along(sim5000pop1$simulations),
                         pop1 = NA, pop2 = NA, pop3 = NA)

cores <- parallel::detectCores() - 1
cl <- parallel::makeCluster(cores, type = "SOCK")
doSNOW::registerDoSNOW(cl)

disp.count <- foreach(i = 1:3, .combine = cbind, .packages=c("raster")) %dopar% {
  counts <- numeric()
  for(w in seq_along(get(paste0("sim5000pop",i))$simulations)){
    ri <- rasterize(do.call(rbind, get(paste0("sim5000pop",i))$simulations[1:w])[,c("x","y")], rast, fun = "count")
    counts[w] <- sum(!is.na(values(ri)))
    print(w)
  } 
  counts
}

parallel::stopCluster(cl)


# save(disp.count, file = "../data/disp_counts.RData")
# load("../data/disp_counts.RData")

#--------##--------##--------##--------##--------##--------##--------##--------#
#--------##--------##--------##--------##--------##--------##--------##--------#

# plotting

# reduce simulations to one data.frame
# simulations <- vector("list", length(sim55$simulations))
# for(g in seq_along(simulations)){
#   simulations[[g]] <- do.call(rbind, sim55$simulations[[g]])
# }
# simulations <- do.call(rbind, simulations)
# simulations <- as.data.frame(simulations[,c("x","y")])
# simulations$value <- 1
# 
# # produce point shapefile
# simSpP <- SpatialPointsDataFrame(coords = simulations, data = simulations,
#                                  proj4string = CRS(proj4string(habitatSuitability)))
# 
# writeOGR(simSpP, dsn = "../data", layer="forwardSimulation50sp", driver="ESRI Shapefile")
# 
# #--------##--------##--------##--------##--------##--------##--------##--------#
# 
# # 5 year step-wise boundary shapes plot
# library(maptools)
# 
# # interval
# interval <- 5
# # buffer radius
# radius <- 200
# # output
# boundaries <- vector("list", floor(length(sim55$simulations)/interval))
# 
# # reduce to generations
# simulations <- vector("list", length(sim55$simulations))
# for(g in seq_along(sim55$simulations)){
#   simulations[[g]] <- do.call(rbind, sim55$simulations[[g]])
# }
# 
# # split into interval
# index = 1 # initiate loop counter
# for(g in seq(interval, length(sim55$simulations), interval)){
#   # chop up simulations
#   chopped <- as.data.frame(do.call(rbind, simulations[(g-interval+1):g])[,c("x","y")])
#   chopped$generation <- g
#   
#   
#   # to point shapefile
#   spdf <- SpatialPointsDataFrame(coords = chopped[,c("x","y")], data = chopped,
#                                  proj4string = CRS(proj4string(habitatSuitability)))
#   
#   # write to shapefile
#   writeOGR(spdf, dsn = "../data/generationIntervals", layer=paste0("loc_gen",g), driver="ESRI Shapefile")
#   
#   #test <- gBuffer(spdf, width=radius, quadsegs=1, capStyle="Round", id=g)
#   
#   # loop counter
#   index = index + 1
# }


# # ppp
# library(graphics)
# smoothScatter(simulations[,"x"], simulations[,"y"], nbin = 256, nrpoints = 100, 
#               colramp = colorRampPalette(c(rgb(1,1,1,0), blues9)), add = T)
