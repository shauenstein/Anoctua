#----------------------------------------------------#
#----------------------------------------------------#
#----------------------------------------------------#
# SH 26-01-2016
# Home Range (HR) analysis:  HR, short rast, dispersal
#----------------------------------------------------#
#----------------------------------------------------#
# set working directory
# setwd("C:/Users/SH/Dropbox/MSc_LittleOwl/Steinkauz_FR/littleowl")
# setwd("D:/InteamDocs/Uni/littleowl")

#----------------------------------------------------#

library(adehabitatHR)
library(sp)
library(raster)
# get rid of temporary files created while doiing raster calculations
showTmpFiles()
removeTmpFiles(h=0) # h being the number of hours before which files created should be removed

library(rgdal)
library(rgeos)
library(gdalUtils); gdal_setInstallation("C:/OSGEO4~2/bin/")
library(INLA)
library(spdep)
library(lme4)
library(nlme)
library(mgcv)
library(randomForest)
library(Hmisc)
library(lubridate)
library(RSAGA) 
#----------------------------------------------------#

#### Utitlisation Density ####

cleaned <- read.csv("data/natal_dispersal_cleaned.csv",
                    stringsAsFactors = FALSE)
addInfo <- read.csv("data/Additional_data_juveniles.csv")
#----------------------------------------------------#

# date as POSIXct
cleaned$date <- as.POSIXct(cleaned$date)

addInfo$Measurement.date <- as.POSIXct(strptime(addInfo$Measurement.date, "%d/%m/%Y", 
                                                tz = "Europe/Berlin"))
addInfo$DOB <- as.POSIXct(strptime(addInfo$DOB, "%d/%m/%Y", 
                                   tz = "Europe/Berlin"))
addInfo$LastKnownToBeAlive <- as.POSIXct(strptime(addInfo$LastKnownToBeAlive, "%d/%m/%Y", 
                                                  tz = "Europe/Berlin"))
addInfo$KnownToBeDead <- as.POSIXct(strptime(addInfo$KnownToBeDead, "%d/%m/%Y", 
                                             tz = "Europe/Berlin"))

# remove rows form addInfo
addInfo <- addInfo[addInfo$Name1 %in% unique(cleaned$id),]

# create an ltraj object to create steps from point data and to measure time difference
# ltraj needs an SPDF object as input
spdf <- SpatialPointsDataFrame(coords = cleaned[,c("posX","posY")], 
                               data = cleaned) 

# add correct projection: DHDN / Gauss-Kruger zone 3
proj4string(spdf) <- CRS("+init=epsg:31467")

ltdf <- as.ltraj(xy = spdf@coords, id = spdf$id, 
                 date = spdf$date, typeII = TRUE)

ltdfSummary <- summary(ltdf)
fixPeriod <- ltdfSummary$date.end - ltdfSummary$date.begin

#----------------------------------------------------#
# plot individual temporal scales
par(mar = c(3,3,0.5,1), cex.lab = 1.3, cex.axis = 1.3, las = 1, mgp = c(1.5,0.2,0), tcl = 0.2, pty = "s")
plot(NA, xlim = c(min(ltdfSummary$date.begin), max(ltdfSummary$date.end)), 
     ylim = c(1,126), type = "n", xlab = "Time", ylab = "Individuals", xaxt = "n", yaxt = "n")
axis.POSIXct(side = 1, at = seq.POSIXt(min(ltdfSummary$date.begin), 
                                       max(ltdfSummary$date.end),
                                       length.out = 5), 
             format = "%m.%Y", labels = TRUE)
for(i in seq_along(ltdf)) lines(c(ltdfSummary$date.begin[i], ltdfSummary$date.end[i]), 
                                rep(i,2), lend = 2, lwd = 3)
# plot fixes
for(i in seq_along(ltdf)){
  idFixes <- cleaned$date[cleaned$id == unique(cleaned$id)[i]]
  points(idFixes, rep(i, length(idFixes)), col = "white", pch = "|", cex = 0.3)
}
# time of death
for(i in seq_along(ltdf)) points(addInfo$KnownToBeDead[i], i, col = "red", pch = 134, cex = 0.8)
#----------------------------------------------------#
#----------------------------------------------------#
# movement characteristics

# starting in year x
spdf$YOB <- NA
for(i in seq_along(unique(spdf$id))){
  spdfId <- spdf[spdf$id == unique(spdf$id)[i],]
  spdf$YOB[spdf$id == unique(spdf$id)[i]] <- year(spdfId@data[1,"date"])
}
# convex hulls by ID
idLocs <- idCHull <- vector(mode = "list", length = length(unique(spdf$id)))
idArea <- 1:3
for(i in seq_along(unique(spdf$id))){
  idLocs[[i]] <- as(spdf[spdf$id == unique(spdf$id)[i],], "SpatialPoints")
  idCHull[[i]] <- gConvexHull(idLocs[[i]])
  idArea[i] <- sapply(slot(idCHull[[i]], "polygons"), slot, "area")/1e+06 #  km^2
}

# plot for sex
par(mar = c(3,3,0.5,1), cex.lab = 1, cex.axis = 1, las = 1, mgp = c(1.5,0.2,0), tcl = -0.2, pty = "s", xpd = FALSE)
hm <- hist(idArea[addInfo$Sex == "m"], breaks = seq(0,350,5), plot = FALSE)
hf <- hist(idArea[addInfo$Sex == "f"], breaks = seq(0,350,5), plot = FALSE)
plot(seq(2.5,347.5,5)-1.1, hm$counts, type = "h", lend = 1, ylim = c(0.4,17), lwd = 2, col = "royalblue3",
     xlab = expression(paste("Area of convex hull [k", m^2, "]")), ylab = "Frequency")
lines(seq(2.5,347.5,5)+1.1, hf$counts, type = "h", lend = 1, lwd = 2, col = "hotpink")
legend("topright", legend = c("",""), pch = c("\u2642", "\u2640"), col = c("royalblue3", "hotpink"), bty = "n", cex = 1.4)


# plot for feeding
par(mar = c(3,3,0.5,1), cex.lab = 1, cex.axis = 1, las = 1, mgp = c(1.5,0.2,0), tcl = -0.2, pty = "s", xpd = FALSE)
hfed <- hist(idArea[addInfo$fed == 1], breaks = seq(0,350,5), plot = FALSE)
hnfed <- hist(idArea[addInfo$fed == 0], breaks = seq(0,350,5), plot = FALSE)
plot(seq(2.5,347.5,5)-1.1, hfed$counts, type = "h", lend = 1, ylim = c(0.4,17), lwd = 2, col = "darkblue",
     xlab = expression(paste("Area of convex hull [k", m^2, "]")), ylab = "Frequency")
lines(seq(2.5,347.5,5)+1.1, hnfed$counts, type = "h", lend = 1, lwd = 2, col = "orange")
legend("topright", legend = c("fed","not fed"), lty = 1, col = c("darkblue", "orange"), bty = "n", cex = 1.4)

# t test for both sex and feeeding, fairly evenly distributed
t.test(idArea[addInfo$Sex == "m" & addInfo$fed == 1], idArea[addInfo$Sex == "f" & addInfo$fed == 1])
t.test(idArea[addInfo$Sex == "m" & addInfo$fed == 0], idArea[addInfo$Sex == "f" & addInfo$fed == 0])
t.test(idArea[addInfo$Sex == "m" & addInfo$fed == 1], idArea[addInfo$Sex == "m" & addInfo$fed == 0])
t.test(idArea[addInfo$Sex == "f" & addInfo$fed == 1], idArea[addInfo$Sex == "f" & addInfo$fed == 0])
t.test(idArea[addInfo$fed == 1], idArea[addInfo$fed == 0])

# avaerage values
mean(idArea[addInfo$Sex == "m" & addInfo$fed == 1], na.rm = T) # 12.96384
mean(idArea[addInfo$Sex == "f" & addInfo$fed == 1], na.rm = T) # 54.05542
mean(idArea[addInfo$Sex == "m" & addInfo$fed == 0], na.rm = T) # 20.52315
mean(idArea[addInfo$Sex == "f" & addInfo$fed == 0], na.rm = T) # 37.33235
#----------------------------------------------------#
# distance from start point
spdf$distToStart <- NA
for(i in seq_along(unique(spdf$id))){
  spdfId <- spdf[spdf$id == unique(spdf$id)[i],]
  startpoint <- as(spdfId[1,], "SpatialPoints")
  spdf$distToStart[spdf$id == unique(spdf$id)[i]] <- sapply(seq(NROW(spdfId)), 
                                                            function(x) gDistance(startpoint,
                                                                                  as(spdfId[x,], "SpatialPoints")))
  
}

par(mfrow = c(1,3), pty = "s", mar = c(2,1.5,1,1), oma = c(1,1.5,0,0), 
    mgp = c(1.4,0.2,0), tcl = 0.2, xpd = NA, cex.axis = 1.2, cex.lab = 1.2)
# starting 2009
plot(NA, 
     ylim = c(0,max(spdf$distToStart[spdf$YOB == 2009])), 
     xlim = range(spdf$date[spdf$YOB == 2009]), xaxt = "n",
     xlab = "Month/Year", ylab = "Distance to start location")
axis.POSIXct(1, at = seq.POSIXt(min(spdf$date[spdf$YOB == 2009]), max(spdf$date[spdf$YOB == 2009]), by = "month"),
             format = "%m/%Y")
for(i in seq_along(unique(spdf$id))) lines(spdf$date[spdf$id == unique(spdf$id)[i] & spdf$YOB == 2009 ], 
                                           spdf$distToStart[spdf$id == unique(spdf$id)[i] & spdf$YOB == 2009],
                                           col = rgb(0,0,0,0.3))

# starting 2010
plot(NA, 
     ylim = c(0,max(spdf$distToStart[spdf$YOB == 2010])), 
     xlim = range(spdf$date[spdf$YOB == 2010]), xaxt = "n",
     xlab = "Month/Year", ylab = "")
axis.POSIXct(1, at = seq.POSIXt(min(spdf$date[spdf$YOB == 2010]), max(spdf$date[spdf$YOB == 2010]), by = "month"),
             format = "%m/%Y")
for(i in seq_along(unique(spdf$id))) lines(spdf$date[spdf$id == unique(spdf$id)[i] & spdf$YOB == 2010 ], 
                                           spdf$distToStart[spdf$id == unique(spdf$id)[i] & spdf$YOB == 2010],
                                           col = rgb(0,0,0,0.3))

# starting 2011
plot(NA, 
     ylim = c(0,max(spdf$distToStart[spdf$YOB == 2011])), 
     xlim = range(spdf$date[spdf$YOB == 2011]), xaxt = "n",
     xlab = "Month/Year", ylab = "")
axis.POSIXct(1, at = seq.POSIXt(min(spdf$date[spdf$YOB == 2011]), max(spdf$date[spdf$YOB == 2011]), by = "month"),
             format = "%m/%Y")
for(i in seq_along(unique(spdf$id))) lines(spdf$date[spdf$id == unique(spdf$id)[i] & spdf$YOB == 2011 ], 
                                           spdf$distToStart[spdf$id == unique(spdf$id)[i] & spdf$YOB == 2011],
                                           col = rgb(0,0,0,0.3))

#----------------------------------------------------#
#----------------------------------------------------#
# estimate diffusion parameters
vv <- BRB.D(ltdf, Tmax= 50, Lmin=0, habitat=NULL)

# compute UD
ud <- BRB(ltdf, D = vv, Tmax = 50, Lmin = 0, hmin=100, grid = 200, b=F,
          same4all=FALSE, extent=0.2, type = "UD")
# ud <- BRB(ltdf, D = vv, Tmax = 5*3600, Lmin = 5, hmin=50, grid = 200, b=F,
#           same4all=FALSE, extent=0.01, type = "ID", maxt = 24*3600, radius = 150) 

vid <- getvolumeUD(ud, standardize = TRUE)

# plot standardised ud
# pdf("texting/ST1_appendix.pdf", width = 6.3, height = 50.4)
# par(mfrow = c(32,4), mar = c(0,0,0,0))
# for(i in 1:126){
#   image(vid[[i]])
#   plot(spdf[spdf$id == unique(spdf$id)[i],], add = T, pch = "+", cex = 0.6, col = rgb(0,0,0,0.2))
#   legend("topleft", legend = paste(unique(spdf$id)[i]), bty = "n")
# } 
# dev.off()

# png("texting/ST1_appendix.png", width = 6.3, height = 50.4, units = "in", res = 300)
# par(mfrow = c(32,4), mar = c(0,0,0,0))
# for(i in 1:126){
#   image(vid[[i]])
#   plot(spdf[spdf$id == unique(spdf$id)[i],], add = T, pch = "+", cex = 0.6, col = rgb(0,0,0,0.2))
#   legend("topleft", legend = paste(unique(spdf$id)[i]), bty = "n")
#   cat(i)
# } 
# dev.off()

# indiv <- 26
# # windows(10,10)
# image(vid[[indiv]])
# contour(vid[[indiv]], add=TRUE, nlevels=1, levels=30,
#         lwd=3, drawlabels=FALSE)
# contour(vid[[indiv]], add=TRUE, nlevels=1, levels=50,
#         lwd=3, drawlabels=FALSE)
# 
# plot(spdf[spdf$id == unique(spdf$id)[indiv],], add = T)

# # obtain values from BRB object, assign to spdf
# spdf$brb <- NA
# for(i in seq_along(unique(spdf$id))){
#   brb <- as(vid[[i]],"SpatialPixelsDataFrame")
#   proj4string(brb) <- proj4string(spdf)
#   spdf$brb[which(spdf$id == unique(spdf$id)[i])] <- over(as(spdf[spdf$id == unique(spdf$id)[i],], 
#                                                             "SpatialPoints"), brb)$n
# }

# check frequencies in classes
# freq <- sapply(1:126, function(x) hist(values(raster(vid[[x]])), breaks = seq(0,100,10), plot = FALSE)$counts)


# sample data
udDf <- data.frame(id = character(), 
                   x = numeric(), 
                   y = numeric(), 
                   ud = numeric(),
                   stringsAsFactors = FALSE)

set.seed(01022016)

for(i in seq_along(unique(spdf$id))){
  # transform to raster
  rast <- raster(vid[[i]])
  proj4string(rast) <- CRS("+init=epsg:31467")
  vals <- values(rast)
  coords <- coordinates(rast)
  
  # sample 100 (when possible) from every class
  for(c in seq.int(10)){
    crit <- which(vals >= ((c-1)*10) &
                    vals < (c*10))
    if(length(crit) > 0){
      if(length(crit) >= 100){
        index <- sample(crit, size = 100, replace = FALSE)
      }else{
        index <- crit
      }
      currentRow <- nrow(udDf) + 1
      # id
      udDf[currentRow : (currentRow + length(index) - 1), "id"] <- rep(unique(spdf$id)[i], length(index))
      # coordinates
      udDf[currentRow : (currentRow + length(index) - 1), c("x", "y")] <- coords[index, ] 
      # ud score
      udDf[currentRow : (currentRow + length(index) - 1), "ud"] <- vals[index] 
    }
  }
}



udSpDf <- SpatialPointsDataFrame(coords = udDf[,c("x","y")], 
                                 data = udDf) 
proj4string(udSpDf) <- CRS("+init=epsg:31467")

extUD <- extent(c(3475564, 3548127, 5395462, 5470803))
#----------------------------------------------------#
#----------------------------------------------------#
#----------------------------------------------------#
# predictors

#----------------------------------------------------#
# landcover
lc <- raster("envVars/landCover/lc_stk/w001000.adf")
lcc <- crop(lc, udSpDf)
#----------------------------------------------------#
# distance to nesting boxes
# nb <- readOGR(dsn = "raw", layer = "all_nestboxes")
# nb <- rasterize(coordinates(nb)[,1:2], raster(lcc), field = 1)
# nbDist <- distance(nb)
# writeRaster(nbDist, file = "envVars/nestingBoxes/distToNestingBoxes.tif", format = "GTiff")
nbDist <- raster("envVars/nestingBoxes/distToNestingBoxes.tif")
# set to 95% quantile
nbDist[values(nbDist) > 1464.261] <- 1464.261
#----------------------------------------------------#
# distance to forest
# f <- lcc
# f[values(f) != 6] <- NA
# do distance calculation with qgis: gdal_proximity.bat
fDist <- raster("envVars/forest/distToForest.tif")
# set to 95% quantile
fDist[values(fDist) > 768.4042] <- 768.4042
#----------------------------------------------------#
# elevation/terrain
elev <- raster("envVars/terrain/elevproj10m/w01000.adf")
proj4string(elev) <- CRS("+init=epsg:31467")
twiElev <- crop(elev, extUD+ c(-2000, 2000, -2000, 2000)) # buffer for  TWI
elev <- crop(elev, extUD) 
# twiElevFoc <- focal(twiElev, w = matrix(1, nrow = 11, ncol = 11), fun = mean, na.rm = TRUE)
# writeRaster(twiElevFoc, file = "envVars/terrain/twiElev10Foc11.tif", format = "GTiff")
twiElevFoc <- raster("envVars/terrain/twiElev10Foc11.tif")
# # compute TWI/CTI with rsaga
# # writeRaster(twiElev, file = "envVars/terrain/TWIelevation10.tif", format = "GTiff")
# rsaga.env(path="C:/Progra~1/SAGA-GIS")
# rsaga.get.modules("ta_preprocessor")
# rsaga.get.usage("ta_preprocessor", 4)
# rsaga.geoprocessor(lib="ta_preprocessor", module=4, param=list(ELEV="envVars/terrain/TWIelevation10.sgrd", FILLED="envVars/terrain/TWIelevation10Wl.sgrd"))
# rsaga.wetness.index(in.dem="envVars/terrain/TWIelevation10.sgrd", out.wetness.index="envVars/terrain/TWI.sgrd")
# rsaga.sgrd.to.esri(in.sgrds="envVars/terrain/TWI.sgrd", out.grids="envVars/terrain/TWI.asc", prec=1, out.path=getwd())
# twi <- raster("envVars/terrain/TWI.asc")
# 
# # derive SAGA Wetness Index:
# rsaga.get.modules("ta_hydrology")
# rsaga.get.usage("ta_hydrology", 15)
# rsaga.geoprocessor(lib="ta_hydrology", module=15, param=list(DEM="envVars/terrain/TWIelevation10Wl.sgrd", TWI="envVars/terrain/TWI.sgrd", SLOPE_TYPE=0)) #local slope
# rsaga.sgrd.to.esri(in.sgrds="envVars/terrain/TWI.sgrd", out.grids="envVars/terrain/TWI.asc", prec=1, out.path=getwd())

# compute TWI in ArcGIS
# cellsize=10
# fd = flowdirection(dem)
# sca = flowaccumulation(fd)
# slope = ( slope(dem) * 1.570796 ) / 90
# tan_slp = con( slope > 0, tan(slope), 0.001 )
# sca_scaled = ( sca + 1 ) * cellsize
# cti = ln ( sca_scaled / tan_slp )

# twi <- raster("envVars/terrain/twi10.tif")
# twi <- focal(twi, w = matrix(1, nrow = 11, ncol = 11), fun = mean, na.rm = TRUE)
# rangeTwi <- range(values(twi), na.rm = TRUE)
# twi <- (twi - rangeTwi[1]) / diff(rangeTwi)
# writeRaster(twi, file = "envVars/terrain/twi10standFoc11.tif")
twi <- raster("envVars/terrain/twi10standFoc11.tif")

#----------------------------------------------------#
# distance to major roads (on 20x20 m grid)
# major roads query from osm data:
#  "type"  =  'motorway' OR  "type"  =  'motorway_link' OR  
# "type"  =  'primary' OR  "type"  =  'primary_link'  OR  
# "type"  =  'secondary'  OR  "type"  =  'secondary_link'  OR  
# "type"  =  'tertiary'  OR  "type" = 'tertiary_link' 
# roads <- readOGR(dsn = "envVars/roads", layer = "majorRoads")
# gdal_rasterize("envVars/roads/majorRoads.shp", "envVars/roads/roads.tif", 
#                burn = 1, te = extUD[c(1,3,2,4)], tr = c(20,20), ignore.full_scan = TRUE, 
#                verbose = FALSE)
# do distance calculation with qgis: gdal_proximity.bat
# roadsDist <- raster("envVars/roads/distToRoads.tif")
#----------------------------------------------------#
# distance to highways (on 20x20 m grid)
# highways query from osm data:
#  "type"  =  'motorway' OR  "type"  =  'motorway_link' OR  
# roads <- readOGR(dsn = "envVars/roads", layer = "highways")
# gdal_rasterize("envVars/roads/highways.shp", "envVars/roads/highways.tif", 
#                burn = 1, te = extUD[c(1,3,2,4)], tr = c(20,20), ignore.full_scan = TRUE, 
#                verbose = FALSE)
# do distance calculation with qgis: gdal_proximity.bat
roadsDist <- raster("envVars/roads/distToHighways.tif")
# set to 95% quantile
roadsDist[values(roadsDist) > 13328.1] <- 13328.1


#----------------------------------------------------#
#----------------------------------------------------#
# extract values from predictor maps
location <- as(udSpDf, "SpatialPoints")

udSpDf$distToForest <- extract(fDist, location) # distance to forest
udSpDf$landuse <- as.factor(extract(lcc, location)) # landuse/landcover
udSpDf$distToNestingBox <- extract(nbDist, location) # distance to nesting box
udSpDf$altitude <- extract(elev, location) # altitude/elevation
udSpDf$tri <- extract(tri, location) # terrain ruggedness index
udSpDf$distToRoads <- extract(roadsDist, location) # distance to major roads

# reverse ud score + scale to 0,1
udRev <- udSpDf$ud/100
udRev <- sapply(seq_along(udRev), function(x) (1 + max(udRev)) - udRev[x] -1)
udSpDf$ud <- udRev

# translate landuse variable
# landclass <- udSpDf$landuse <- as.numeric(as.character(udSpDf$landuse))
# udSpDf$landuse[landclass == 1] <- "agriCult" # 1 Ackerflaeche/agriculture
# udSpDf$landuse[landclass == 6] <- "forest" # 6	Wald/forest
# udSpDf$landuse[landclass == 9] <- "scarceVeg" # 9	wenig oder keine Vegetation/scarce vegetation
# udSpDf$landuse[landclass == 10] <- "wetland" # 10	Feuchtgebiet/wetland
# udSpDf$landuse[landclass == 11] <- "water" # 11	Wasser/water
# udSpDf$landuse[landclass == 21] <- "orchard" # 21	Obstwiese/orchard
# udSpDf$landuse[landclass == 25] <- "other" # 25	sonstige Flaeche/other
# udSpDf$landuse[landclass == 30] <- "meadow" # 30	Wiese/meadow
# udSpDf$landuse[landclass == 41] <- "urban" # 41	Siedlung/settlement area
# udSpDf$landuse[landclass == 80] <- "hedge" # 80	Feldgehoelz oder Hecke/hedgerow


# omit NAs
udSpDf@data <- na.omit(udSpDf@data)

# unique ids
udSpDf$uniqueId <- seq.int(NROW(udSpDf))


write.csv(udSpDf@data, "data/sampledUD.csv", row.names = FALSE)


# # calculate minimum distance 
# udDist <- dnearneigh(coordinates(udSpDf), 0, 1000)
# distList <- nbdists(udDist, coordinates(udSpDf))
# d <- unlist(lapply(distList, FUN=function(x) min(x, na.rm = T)))
# hist(d, breaks = 1000, xlim=c(0,30))
#----------------------------------------------------#
#----------------------------------------------------#
#----------------------------------------------------#
#### Habitat selection #### 

udDf <- read.csv("data/sampledUD.csv")
addInfo <- read.csv("data/Additional_data_juveniles.csv")

# add sex to data
udDf$sex <- factor(NA); levels(udDf$sex) <- c(levels(addInfo$Sex), "unknown")
for(i in unique(addInfo$Name1)) udDf$sex[udDf$id == i] <- addInfo$Sex[addInfo$Name1 == i]
udDf$sex[is.na(udDf$sex)] <- factor("unknown")
#----------------------------------------------------#
# Predictor distribution
# hist(sqrt(udDf$distToForest)) # sqrt
# table(udDf$landuse) # unbalanced, but no solution
# hist(log(udDf$distToNestingBox + 5)) # log(x + 5)
# hist(udDf$altitude)# normal
# hist(sqrt(udDf$tri)) # sqrt
# hist(sqrt(udDf$distToRoads)) # sqrt

udDfs <- udDf
udDfs$distToForest <- sqrt(udDfs$distToForest)
udDfs$distToNestingBox <- log(udDfs$distToNestingBox + 5)
udDfs$tri <- sqrt(udDfs$tri)
udDfs$distToRoads <- sqrt(udDfs$distToRoads)

#----------------------------------------------------#
# scale variables
for(i in names(udDfs)[c(2,3,5,7,8,9,10)]) udDfs[,i] <- scale(udDfs[,i])

#----------------------------------------------------#
# predictor correlation
# plot(varclus(as.matrix(udDfs[,c(5,7,8,9,10)]))) # not an issue

#----------------------------------------------------#
# shrink ud 
udDfs$ud <- (udDfs$ud - 0.5) / 1.02 + 0.5

# categorise HR (2), Dispersal (1), Exclusion (0)
udDfs$disp <- ifelse(udDfs$ud >= 0.6, NA, 
                     ifelse(udDfs$ud >= 0.1, 1, 0))
udDfs$hr <- ifelse(udDfs$ud >= 0.6, 1, 0)

# dummy to predict with random effect turned off
dum <- rep(1, NROW(udDfs))

# random sample of data
# set.seed(04022016)
# udDfs <- udDfs[sample(seq.int(NROW(udDfs)), size = 10000, replace = FALSE), ]

# #----------------------------------------------------#
# INLA, use combined model BYM, which
## contains both the "iid" and the "besag" model

# define neighborhood # encompass diagonals in 20x20 grid
# neighs <- dnearneigh(as.matrix(udDf[,c("x", "y")]), d1 = 0.01, d2 = 100) 
# spId <- seq.int(NROW(udDf))
# nb2INLA("analysis/nc.txt", neighs)

# fmInla <- inla(ud ~ as.factor(landuse) + distToForest + I(distToForest^2) 
#                   + distToNestingBox + altitude + distToRoads
#                   + f(spId, model = "besag", graph = "analysis/nc.txt")
#                   + f(id ,model="iid"),
#                data = udDfs,
#                family = "gaussian", control.family = list(link = "logit"),
#                control.predictor = list(compute = TRUE),
#                control.compute = list(dic = TRUE, cpo = TRUE))

#----------------------------------------------------#
# GLMM
# fmGlmm <- glmer(ud ~ as.factor(landuse) + distToForest + I(distToForest^2) 
#                     + distToNestingBox + altitude + I(altitude^2)
#                     + distToRoads 
#                     + (1|id), 
#                 control = glmerControl(optCtrl = list(maxfun = 100000)),
#                 family = gaussian(link = "logit"), data = udDfs)

#----------------------------------------------------#
# continuous ud
system.time(
  fmGam <- gam(ud ~ as.factor(landuse) + s(distToForest, bs = "ts") 
               + s(distToNestingBox, bs = "ts") + s(altitude, bs = "ts") 
               + s(distToRoads, bs = "ts") 
               + s(id, bs = "re", by = dum),# + s(x, y), 
               family = gaussian(link = "logit"),
               data = udDfs)
) # 2335.943
summary(fmGam)
AIC(fmGam) # spat: -14845.59 nonspat: -14228.78

# # binary: HR selection
# system.time(
#   fmGamHr <- gam(hr ~ as.factor(landuse) + s(distToForest, bs = "cs") 
#                + s(distToNestingBox, bs = "cs") + s(altitude, bs = "cs") 
#                + s(distToRoads, bs = "cs") 
#                + s(id, bs = "re") + s(x, y), 
#                family = binomial,
#                data = udDfs)
# ) # 2971.032  
# summary(fmGamHr)
# AIC(fmGamHr) # 50599.35

# binary: Dispersal selection
system.time(
  fmGamDisp <- gam(disp ~ as.factor(landuse) + s(distToForest, bs = "ts") 
                   + s(distToNestingBox, bs = "ts") + s(altitude, bs = "ts") 
                   + s(distToRoads, bs = "ts") 
                   + s(id, bs = "re", by = dum),# + s(x, y), 
                   family = binomial,
                   data = udDfs)
) # 766.897
summary(fmGamDisp)
AIC(fmGamDisp) # 47032.96

#----------------------------------------------------#
# # random forest
# # multinomial response
# udDfs$cat <- ifelse(is.na(udDfs$disp), 2, udDfs$disp)
# udDfs$landuse <- as.factor(udDfs$landuse)
# system.time(
#   fRf <- randomForest(cat ~ landuse + distToForest 
#                       + distToNestingBox + altitude 
#                       + distToRoads,
#                       data = udDfs, ntree = 200,
#                       importance = TRUE, strata = id)
# )
#----------------------------------------------------#

load("analysis/fm.RData")


# windows(10,10)
# par(mfrow=c(3,2))
# plot(fmGam)
# windows(10,10)
# par(mfrow=c(3,2))
# plot(fmGamDisp)

#----------------------------------------------------#
# Effect plots, fix other variable on 0 (mean)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# landuse
newLc <- unique(udDfs$landuse)
predsLc <- predict(fmGam, type = "response", se.fit = TRUE,
                   newdata = data.frame("landuse" = as.factor(newLc), 
                                        "distToForest" = 0,
                                        "distToNestingBox" = 0,
                                        "altitude" = 0,
                                        "distToRoads" = 0,
                                        "id" = unique(udDf$id)[1], # any id is fine, since r.e. "turned off"
                                        "dum" = 0, # r.e. "turned off"
                                        "x" = 0,
                                        "y" = 0))

predsLcDisp <- predict(fmGamDisp, type = "response", se.fit = TRUE,
                       newdata = data.frame("landuse" = as.factor(newLc), 
                                            "distToForest" = 0,
                                            "distToNestingBox" = 0,
                                            "altitude" = 0,
                                            "distToRoads" = 0,
                                            "id" = unique(udDf$id)[1], # any id is fine, since r.e. "turned off"
                                            "dum" = 0, # r.e. "turned off"
                                            "x" = 0,
                                            "y" = 0))

par(mar = c(4,1.5,0.5,0.5), mgp = c(2, 0.2, 0), tcl = 0.2, las = 1, cex.axis = 1.2, cex.lab = 1.2, pty = "s")

plot(1:8-0.1, predsLc$fit-mean(predsLc$fit), ylim = c(-0.5,0.5), pch = "-", 
     cex = 2, col = "darkred", xlab = "", ylab = "Centered prediction", xaxt = "n")
axis(1, at = 1:8, labels = FALSE)
text(1:8, -0.555,labels = c("Agriculture", "Meadow", "Orchard", "Forest", 
                            "Urban", "Hedge", "Other", "Water"), 
     srt = 30, adj = 1, xpd = TRUE, cex = 1.2)

arrows(x0 = 1:8-0.1, y0 = predsLc$fit-mean(predsLc$fit), 
       y1 = predsLc$fit-mean(predsLc$fit) + 2 * predsLc$se.fit, 
       angle = 90, length = 0.1, col = "darkred")
arrows(x0 = 1:8-0.1, y0 = predsLc$fit-mean(predsLc$fit), 
       y1 = predsLc$fit-mean(predsLc$fit) - 2 * predsLc$se.fit, 
       angle = 90, length = 0.1, col = "darkred")

points(1:8+0.1, predsLcDisp$fit-mean(predsLcDisp$fit), pch = "-", 
       cex = 2, col = "darkblue")
arrows(x0 = 1:8+0.1, y0 = predsLcDisp$fit-mean(predsLcDisp$fit), 
       y1 = predsLcDisp$fit-mean(predsLcDisp$fit) + 2 * predsLcDisp$se.fit, 
       angle = 90, length = 0.1, col = "darkblue")
arrows(x0 = 1:8+0.1, y0 = predsLcDisp$fit-mean(predsLcDisp$fit), 
       y1 = predsLcDisp$fit-mean(predsLcDisp$fit) - 2 * predsLcDisp$se.fit, 
       angle = 90, length = 0.1, col = "darkblue")


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# distToForest

newDF <- seq(min(udDfs$distToForest), max(udDfs$distToForest), len = 100)
oldDF <- seq(min(udDf$distToForest), max(udDf$distToForest), len = 100)
predsDF <- predict(fmGam, type = "response", se.fit = TRUE,
                   newdata = data.frame("landuse" = factor(1), 
                                        "distToForest" = newDF,
                                        "distToNestingBox" = 0,
                                        "altitude" = 0,
                                        "distToRoads" = 0,
                                        "id" = unique(udDf$id)[1], # any id is fine, since r.e. "turned off"
                                        "dum" = 0, # r.e. "turned off"
                                        "x" = 0,
                                        "y" = 0))

predsDFDisp <- predict(fmGamDisp, type = "response", se.fit = TRUE,
                       newdata = data.frame("landuse" = factor(1), 
                                            "distToForest" = newDF,
                                            "distToNestingBox" = 0,
                                            "altitude" = 0,
                                            "distToRoads" = 0,
                                            "id" = unique(udDf$id)[1], # any id is fine, since r.e. "turned off"
                                            "dum" = 0, # r.e. "turned off"
                                            "x" = 0,
                                            "y" = 0))

par(mar = c(3.5,1.5,0.5,0.5), mgp = c(2, 0.2, 0), tcl = 0.2, las = 1, cex.axis = 1.2, cex.lab = 1.2, pty = "s")

plot(oldDF, predsDF$fit-mean(predsDF$fit), ylim = c(-0.5,0.5), type = "n",
     xlab = "Distance to forest [m]", ylab = "Centered prediction")

polygon(c(oldDF, rev(oldDF)), c(predsDF$fit-mean(predsDF$fit) + 2 * predsDF$se.fit,
                                rev(predsDF$fit-mean(predsDF$fit) - 2 * predsDF$se.fit)), 
        col=rgb(0.54,0,0,0.2),border=NA)

polygon(c(oldDF, rev(oldDF)), c(predsDFDisp$fit-mean(predsDFDisp$fit) + 2 * predsDFDisp$se.fit,
                                rev(predsDFDisp$fit-mean(predsDFDisp$fit) - 2 * predsDFDisp$se.fit)), 
        col=rgb(0,0,0.54,0.2),border=NA)

lines(oldDF, predsDF$fit-mean(predsDF$fit), col = "darkred", lwd = 2)

lines(oldDF, predsDFDisp$fit-mean(predsDFDisp$fit), col = "darkblue", lwd = 2)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# distToNestingBox
newDN <- seq(min(udDfs$distToNestingBox), max(udDfs$distToNestingBox), len = 100)
oldDN <- seq(min(udDf$distToNestingBox), max(udDf$distToNestingBox), len = 100)
predsDN <- predict(fmGam, type = "response", se.fit = TRUE,
                   newdata = data.frame("landuse" = factor(1), 
                                        "distToForest" = 0,
                                        "distToNestingBox" = newDN,
                                        "altitude" = 0,
                                        "distToRoads" = 0,
                                        "id" = unique(udDf$id)[1], # any id is fine, since r.e. "turned off"
                                        "dum" = 0, # r.e. "turned off"
                                        "x" = 0,
                                        "y" = 0))

predsDNDisp <- predict(fmGamDisp, type = "response", se.fit = TRUE,
                       newdata = data.frame("landuse" = factor(1), 
                                            "distToForest" = 0,
                                            "distToNestingBox" = newDN,
                                            "altitude" = 0,
                                            "distToRoads" = 0,
                                            "id" = unique(udDf$id)[1], # any id is fine, since r.e. "turned off"
                                            "dum" = 0, # r.e. "turned off"
                                            "x" = 0,
                                            "y" = 0))

par(mar = c(3.5,1.5,0.5,0.5), mgp = c(2, 0.2, 0), tcl = 0.2, las = 1, cex.axis = 1.2, cex.lab = 1.2, pty = "s")

plot(oldDN, predsDN$fit-mean(predsDN$fit), ylim = c(-0.5,0.5), type = "n",
     xlab = "Distance to nesting box [m]", ylab = "Centered prediction")

polygon(c(oldDN, rev(oldDN)), c(predsDN$fit-mean(predsDN$fit) + 2 * predsDN$se.fit,
                                rev(predsDN$fit-mean(predsDN$fit) - 2 * predsDN$se.fit)), 
        col=rgb(0.54,0,0,0.2),border=NA)

polygon(c(oldDN, rev(oldDN)), c(predsDNDisp$fit-mean(predsDNDisp$fit) + 2 * predsDNDisp$se.fit,
                                rev(predsDNDisp$fit-mean(predsDNDisp$fit) - 2 * predsDNDisp$se.fit)), 
        col=rgb(0,0,0.54,0.2),border=NA)

lines(oldDN, predsDN$fit-mean(predsDN$fit), col = "darkred", lwd = 2)

lines(oldDN, predsDNDisp$fit-mean(predsDNDisp$fit), col = "darkblue", lwd = 2)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# altitude

newAlt <- seq(min(udDfs$altitude), max(udDfs$altitude), len = 100)
oldAlt <- seq(min(udDf$altitude), max(udDf$altitude), len = 100)
predsAlt <- predict(fmGam, type = "response", se.fit = TRUE,
                    newdata = data.frame("landuse" = factor(1), 
                                         "distToForest" = 0,
                                         "distToNestingBox" = 0,
                                         "altitude" = newAlt,
                                         "distToRoads" = 0,
                                         "id" = unique(udDf$id)[1], # any id is fine, since r.e. "turned off"
                                         "dum" = 0, # r.e. "turned off"
                                         "x" = 0,
                                         "y" = 0))

predsAltDisp <- predict(fmGamDisp, type = "response", se.fit = TRUE,
                        newdata = data.frame("landuse" = factor(1), 
                                             "distToForest" = 0,
                                             "distToNestingBox" = 0,
                                             "altitude" = newAlt,
                                             "distToRoads" = 0,
                                             "id" = unique(udDf$id)[1], # any id is fine, since r.e. "turned off"
                                             "dum" = 0, # r.e. "turned off"
                                             "x" = 0,
                                             "y" = 0))

par(mar = c(3.5,1.5,0.5,0.5), mgp = c(2, 0.2, 0), tcl = 0.2, las = 1, cex.axis = 1.2, cex.lab = 1.2, pty = "s")

plot(oldAlt, predsAlt$fit-mean(predsAlt$fit), ylim = c(-0.5,0.5), type = "n",
     xlab = "Altitude [MASL]", ylab = "Centered prediction")

polygon(c(oldAlt, rev(oldAlt)), c(predsAlt$fit-mean(predsAlt$fit) + 2 * predsAlt$se.fit,
                                  rev(predsAlt$fit-mean(predsAlt$fit) - 2 * predsAlt$se.fit)), 
        col=rgb(0.54,0,0,0.2),border=NA)

polygon(c(oldAlt, rev(oldAlt)), c(predsAltDisp$fit-mean(predsAltDisp$fit) + 2 * predsAltDisp$se.fit,
                                  rev(predsAltDisp$fit-mean(predsAltDisp$fit) - 2 * predsAltDisp$se.fit)), 
        col=rgb(0,0,0.54,0.2),border=NA)

lines(oldAlt, predsAlt$fit-mean(predsAlt$fit), col = "darkred", lwd = 2)

lines(oldAlt, predsAltDisp$fit-mean(predsAltDisp$fit), col = "darkblue", lwd = 2)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# distToRoads

newDR <- seq(min(udDfs$distToRoads), max(udDfs$distToRoads), len = 100)
oldDR <- seq(min(udDf$distToRoads), max(udDf$distToRoads), len = 100)
predsDR <- predict(fmGam, type = "response", se.fit = TRUE,
                   newdata = data.frame("landuse" = factor(1), 
                                        "distToForest" = 0,
                                        "distToNestingBox" = 0,
                                        "altitude" = 0,
                                        "distToRoads" = newDR,
                                        "id" = unique(udDf$id)[1], # any id is fine, since r.e. "turned off"
                                        "dum" = 0, # r.e. "turned off"
                                        "x" = 0,
                                        "y" = 0))

predsDRDisp <- predict(fmGamDisp, type = "response", se.fit = TRUE,
                       newdata = data.frame("landuse" = factor(1), 
                                            "distToForest" = 0,
                                            "distToNestingBox" = 0,
                                            "altitude" = 0,
                                            "distToRoads" = newDR,
                                            "id" = unique(udDf$id)[1], # any id is fine, since r.e. "turned off"
                                            "dum" = 0, # r.e. "turned off"
                                            "x" = 0,
                                            "y" = 0))

par(mar = c(3.5,1.5,0.5,0.5), mgp = c(2, 0.2, 0), tcl = 0.2, las = 1, cex.axis = 1.2, cex.lab = 1.2, pty = "s")

plot(oldDR, predsDR$fit-mean(predsDR$fit), ylim = c(-0.5,0.5), type = "n",
     xlab = "Distance to major roads [m]", ylab = "Centered prediction")

polygon(c(oldDR, rev(oldDR)), c(predsDR$fit-mean(predsDR$fit) + 2 * predsDR$se.fit,
                                rev(predsDR$fit-mean(predsDR$fit) - 2 * predsDR$se.fit)), 
        col=rgb(0.54,0,0,0.2),border=NA)

polygon(c(oldDR, rev(oldDR)), c(predsDRDisp$fit-mean(predsDRDisp$fit) + 2 * predsDRDisp$se.fit,
                                rev(predsDRDisp$fit-mean(predsDRDisp$fit) - 2 * predsDRDisp$se.fit)), 
        col=rgb(0,0,0.54,0.2),border=NA)

lines(oldDR, predsDR$fit-mean(predsDR$fit), col = "darkred", lwd = 2)

lines(oldDR, predsDRDisp$fit-mean(predsDRDisp$fit), col = "darkblue", lwd = 2)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#########----------------------------##########
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# mapped predictions and residuals
# preds <- predict(fmGam, type = "response")
# predsDisp <- predict(fmGamDisp, type = "response")
# plot(udDf$x, udDf$y, pch=15, cex=0.5, col=terrain.colors(100)[preds*100], xlab="x", ylab="y")
# plot(udDf$x, udDf$y, pch=15, cex=0.5, col=terrain.colors(100)[predsDisp*100], xlab="x", ylab="y")
# points(udDfH$posX, udDfH$posY, pch = "+")

# scaled environmental variables on uniform grid

fDist <- raster("envVars/forest/distToForest.tif")
fDmean <- mean(sqrt(udDf$distToForest)); fDsd <- sd(sqrt(udDf$distToForest))

lc <- raster("envVars/landCover/lc_stk/w001000.adf")
lcc <- crop(lc, fDist)


nbDist <- raster("envVars/nestingBoxes/distToNestingBoxes.tif")
nbDmean <- mean(log(udDf$distToNestingBox + 5)); nbDsd <- sd(log(udDf$distToNestingBox + 5))

elev <- raster("envVars/terrain/elevproj10m/w01000.adf")
proj4string(elev) <- CRS("+init=epsg:31467")
elev <- projectRaster(elev, fDist)
Emean <- mean(udDf$altitude); Esd <- sd(udDf$altitude)

roadsDist <- raster("envVars/roads/distToRoads.tif")
rDmean <- mean(sqrt(udDf$distToRoads)); rDsd <- sd(sqrt(udDf$distToRoads))
roadsDist <- disaggregate(roadsDist, fact = 2)

# aggregate to 30x30 m grid
fDist <- aggregate(fDist, fact = 3)
lcc <- aggregate(lcc, fact = 3, fun = modal)
nbDist <- aggregate(nbDist, fact = 3)
elev <- aggregate(elev, fact = 3)
roadsDist <- aggregate(roadsDist, fact = 3)

xCol <- xFromCol(fDist)
yRow <- yFromRow(fDist)
xmean <- mean(udDf$x); xsd <- sd(udDf$x)
ymean <- mean(udDf$y); ysd <- sd(udDf$y)
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# predict and write predictions row-wise
con <- file("udPredictionsNonSpat.asc", "w")
con2 <- file("dispPredictionsNonSpat.asc", "w")

writeLines("NCOLS 2419 
           NROWS 2512 
           XLLCORNER 3475566.95105929 
           YLLCORNER 5395442.20915809 
           CELLSIZE 30 
           NODATA_value -9999", con =   con) 

writeLines("NCOLS 2419 
           NROWS 2512 
           XLLCORNER 3475566.95105929 
           YLLCORNER 5395442.20915809 
           CELLSIZE 30 
           NODATA_value -9999", con =   con2) 

for(i in seq(NROW(fDist))){
  newdat <- data.frame("landuse" = as.factor(lcc[i,]), 
                       "distToForest" = (sqrt(fDist[i,]) - fDmean) / fDsd,
                       "distToNestingBox" = (log(nbDist[i,] + 5) - nbDmean) / nbDsd,
                       "altitude" = (elev[i,] - Emean) / Esd,
                       "distToRoads" = (sqrt(roadsDist[i,]) - rDmean) / rDsd,
                       "id" = unique(udDf$id)[1], # any id is fine, since r.e. "turned off"
                       "dum" = 0, # r.e. "turned off"
                       "x" = (xCol - xmean) / xsd,
                       "y" = (yRow[i] - ymean) / ysd)
  
  
  preds <- predict(fmGam, newdata = newdat, 
                   type = "response")
  preds2 <- predict(fmGamDisp, newdata = newdat, 
                    type = "response")
  
  cat(preds, "\n", file = con)
  cat(preds2, "\n", file = con2)
  
  cat("row",i, " ")
}
close(con); close(con2)


# udPreds <- raster("analysis/udPredictions.tif")
# dispPreds <- raster("analysis/dispPredictions.tif")
# udPredsNonSpat <- raster("analysis/udPredictionsNonSpat.asc")
# dispPredsNonSpat <- raster("analysis/dispPredictionsNonSpat.tif")
# proj4string(udPredsNonSpat) <- CRS("+init=epsg:31467")
# writeRaster(udPredsNonSpat, file = "analysis/udPredictionsNonSpat.tif", overwrite = TRUE)

#------------------------------------------------------#
#------------------------------------------------------#
#------------------------------------------------------#
# plot prediction maps

library(ggplot2); library(grid); library(gridExtra)

#------------------------------------------------------#

udPredsNonSpat <- raster("analysis/udPredictionsNonSpat.tif")
dispPredsNonSpat <- raster("analysis/dispPredictionsNonSpat.tif")

mapsDf <- data.frame(ud = values(udPredsNonSpat), disp = values(dispPredsNonSpat), coordinates(udPredsNonSpat))
# center predictions
mapsDf$ud <- mapsDf$ud - mean(mapsDf$ud, na.rm = TRUE)
mapsDf$disp <- mapsDf$disp - mean(mapsDf$disp, na.rm = TRUE)
rm(udPredsNonSpat, dispPredsNonSpat)

p1 <- ggplot(data = mapsDf) + 
  geom_raster(aes(x, y, fill = ud)) + 
  scale_fill_gradient2("", limits = c(-0.3,0.8), low="grey94", mid = "#d7301f",high="#190000", 
                       midpoint = 0.25, guide = FALSE) + 
  scale_x_continuous(limits=range(mapsDf$x), expand = c(0, 0)) +
  scale_y_continuous(limits=range(mapsDf$y), expand = c(0, 0)) +
  coord_fixed() + 
  xlab("Easting")+
  ylab("Northing")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        panel.margin = unit(c(0,0,0,0),"cm"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20, angle=90),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 30), 
        strip.background = element_blank(),
        # axis.ticks = element_blank(),
        panel.background = element_blank()) 

p2 <- ggplot(data = mapsDf) + 
  geom_raster(aes(x, y, fill = disp)) + 
  scale_fill_gradient2("", limits = c(-0.3,0.8), low="grey94", mid = "#d7301f",high="#190000", 
                       midpoint = 0.25, guide = guide_colourbar(ticks = FALSE, title.position = "left", 
                                                                title.vjust = 0.4, title.hjust = 0.6,
                                                                scale_colour = "black")) +
  scale_x_continuous(limits=range(mapsDf$x), expand = c(0, 0)) +
  scale_y_continuous(limits=range(mapsDf$y), expand = c(0, 0)) +
  coord_fixed() + 
  xlab("Easting")+
  ylab("")+
  theme(plot.margin = unit(c(0,2.5,0,-0.75), "cm"),
        panel.margin = unit(c(0,0,0,0),"cm"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20, angle=90),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_blank(),
        legend.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(1.08, 0.7),
        legend.key.width = unit(0.4, "cm"),
        legend.key.height = unit(1.2, "cm"), 
        legend.title=element_blank(),
        legend.text=element_text(size=18),
        legend.direction = "vertical",
        strip.text.x = element_text(size = 30), 
        strip.background = element_blank(),
        # axis.ticks = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank())  

png("texting/Pmaps_centered.png", width = 12, height = 6, units = "in", res = 300)
grid.arrange(p1, p2, ncol=2)
dev.off()

#------------------------------------------------------#
#------------------------------------------------------#
#------------------------------------------------------#
# # predict randomFores
# predsRf <- predict(fRf, newdata = data.frame("landuse" = as.factor(values(lcc)), 
#                                              "distToForest" = (sqrt(values(fDist)) - fDmean) / fDsd,
#                                              "distToNestingBox" = (log(values(nbDist) + 5) - nbDmean) / nbDsd,
#                                              "altitude" = (values(elev) - Emean) / Esd,
#                                              "distToRoads" = (sqrt(values(roadsDist)) - rDmean) / rDsd),
#                    type = "response")
#                                              
# rfPredictions <- fDist
# values(rfPredictions) <- predsRf
# # writeRaster(rfPredictions, file = "rfPredictions.tif", format = "GTiff")
# 
# rfPreds <- raster("rfPredictions.tif")
