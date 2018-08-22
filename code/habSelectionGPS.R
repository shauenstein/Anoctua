#----------------------------------------------------#
#----------------------------------------------------#
#----------------------------------------------------#
# SH 26-01-2016
# Home Range (HR) analysis:  HR, short rast, dispersal
#----------------------------------------------------#
#----------------------------------------------------#
# set working directory
# setwd("C:/Users/SH/Dropbox/MSc_LittleOwl/Steinkauz_FR")
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
library(mgcv)
library(randomForest)
library(Hmisc)
#----------------------------------------------------#

#### Utitlisation Density ####

cleaned <- read.csv("data/natal_dispersal_cleaned.csv",
                    stringsAsFactors = FALSE)

#----------------------------------------------------#

# date as POSIXct
cleaned$date <- as.POSIXct(cleaned$date)


# create an ltraj object to create steps from point data and to measure time difference
# ltraj needs an SPDF object as input
spdf <- SpatialPointsDataFrame(coords = cleaned[,c("posX","posY")], 
                               data = cleaned) 

# add correct projection NAD83 UTM Zone 11
proj4string(spdf) <- CRS("+init=epsg:31467")

ltdf <- as.ltraj(xy = spdf@coords, id = spdf$id, 
                 date = spdf$date, typeII = TRUE)

# ltdf <- as.ltraj(xy = spdf@coords[spdf$id == unique(spdf$id)[78],], id = unique(spdf$id)[78], 
#                  date = spdf$date[spdf$id == unique(spdf$id)[78]], typeII = TRUE)
#----------------------------------------------------#

# estimate diffusion parameters
vv <- BRB.D(ltdf, Tmax= 50, Lmin=0, habitat=NULL)

# compute UD
ud <- BRB(ltdf, D = vv, Tmax = 50, Lmin = 0, hmin=100, grid = 200, b=F,
          same4all=FALSE, extent=0.2, type = "UD", maxt = 3*24*3600, 
          radius = 300)
# ud <- BRB(ltdf, D = vv, Tmax = 5*3600, Lmin = 5, hmin=50, grid = 200, b=F,
#           same4all=FALSE, extent=0.01, type = "ID", maxt = 24*3600, radius = 150) 

vid <- getvolumeUD(ud, standardize = TRUE)

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
spdf$ud <- NA
for(i in seq_along(unique(spdf$id))){
  ud <- as(vid[[i]],"SpatialPixelsDataFrame")
  proj4string(ud) <- CRS("+init=epsg:31467")
  spdf$ud[which(spdf$id == unique(spdf$id)[i])] <- over(as(spdf[spdf$id == unique(spdf$id)[i],], 
                                                            "SpatialPoints"), ud)$n
}


udSpDf <- spdf



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
#----------------------------------------------------#
# distance to forest
# f <- lcc
# f[values(f) != 6] <- NA
# do distance calculation with qgis: gdal_proximity.bat
fDist <- raster("envVars/forest/distToForest.tif")
#----------------------------------------------------#
# elevation/terrain
elev <- raster("envVars/terrain/elevproj10m/w01000.adf")
proj4string(elev) <- CRS("+init=epsg:31467")
elev <- crop(elev, udSpDf)
tri <- terrain(elev, opt = "TRI")
#----------------------------------------------------#
# distance to major roads (on 20x20 m grid)
# major roads query from osm data:
#  "type"  =  'motorway' OR  "type"  =  'motorway_link' OR  
# "type"  =  'primary' OR  "type"  =  'primary_link'  OR  
# "type"  =  'secondary'  OR  "type"  =  'secondary_link'  OR  
# "type"  =  'tertiary'  OR  "type" = 'tertiary_link' 
# roads <- readOGR(dsn = "envVars/roads", layer = "majorRoads")
# gdal_rasterize("envVars/roads/majorRoads.shp", "envVars/roads/roads.tif", 
#                burn = 1, te = extent(lcc)[c(1,3,2,4)], tr = c(20,20), ignore.full_scan = TRUE, 
#                verbose = FALSE)
# do distance calculation with qgis: gdal_proximity.bat
# roadsDist <- raster("envVars/roads/distToRoads.tif")
#----------------------------------------------------#
# distance to highways (on 20x20 m grid)
# highways query from osm data:
#  "type"  =  'motorway' OR  "type"  =  'motorway_link' OR  
# roads <- readOGR(dsn = "envVars/roads", layer = "highways")
# gdal_rasterize("envVars/roads/highways.shp", "envVars/roads/highways.tif", 
#                burn = 1, te = extent(udSpDf)[c(1,3,2,4)], tr = c(20,20), ignore.full_scan = TRUE, 
#                verbose = FALSE)
# do distance calculation with qgis: gdal_proximity.bat
roadsDist <- raster("envVars/roads/distToHighways.tif")

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

# check max values:
max(udSpDf$distToForest)
quantile(udSpDf$distToForest, 0.95) #768.4042 
max(udSpDf$distToNestingBox)
quantile(udSpDf$distToNestingBox, 0.95) # 1464.261
max(udSpDf$distToRoads)
quantile(udSpDf$distToRoads, 0.95)

# reverse ud score + scale to 0,1
udRev <- udSpDf$ud/100
udRev <- sapply(seq_along(udRev), function(x) (1 + max(udRev)) - udRev[x] -1)
udSpDf$ud <- udRev

# translate landuse variable
landclass <- udSpDf$landuse <- as.numeric(as.character(udSpDf$landuse))
udSpDf$landuse[landclass == 1] <- "agriCult" # 1 Ackerflaeche/agriculture
udSpDf$landuse[landclass == 6] <- "forest" # 6	Wald/forest
udSpDf$landuse[landclass == 9] <- "scarceVeg" # 9	wenig oder keine Vegetation/scarce vegetation
udSpDf$landuse[landclass == 10] <- "wetland" # 10	Feuchtgebiet/wetland
udSpDf$landuse[landclass == 11] <- "water" # 11	Wasser/water
udSpDf$landuse[landclass == 21] <- "orchard" # 21	Obstwiese/orchard
udSpDf$landuse[landclass == 25] <- "other" # 25	sonstige Flaeche/other
udSpDf$landuse[landclass == 30] <- "meadow" # 30	Wiese/meadow
udSpDf$landuse[landclass == 41] <- "urban" # 41	Siedlung/settlement area
udSpDf$landuse[landclass == 80] <- "hedge" # 80	Feldgehoelz oder Hecke/hedgerow


# omit NAs
udSpDf@data <- na.omit(udSpDf@data)

# unique ids
udSpDf$uniqueId <- seq.int(NROW(udSpDf))


write.csv(udSpDf@data, "data/gpsUD.csv", row.names = FALSE)


# # calculate minimum distance 
# udDist <- dnearneigh(coordinates(udSpDf), 0, 1000)
# distList <- nbdists(udDist, coordinates(udSpDf))
# d <- unlist(lapply(distList, FUN=function(x) min(x, na.rm = T)))
# hist(d, breaks = 1000, xlim=c(0,30))
#----------------------------------------------------#
#----------------------------------------------------#
#----------------------------------------------------#
#### Habitat selection on GPS locations ####

udDfH <- read.csv("data/gpsUD.csv")
udDfH$date <- udDfH$accuracy <- NULL
sex <- read.csv("data/Additional_data_juveniles.csv")

# add sex to data
udDfH$sex <- factor(NA); levels(udDfH$sex) <- c(levels(sex$Sex), "unknown")
for(i in unique(sex$Name1)) udDf$sex[udDfH$id == i] <- sex$Sex[sex$Name1 == i]
udDf$sex[is.na(udDfH$sex)] <- factor("unknown")
#----------------------------------------------------#
# Predictor distribution
# hist(sqrt(udDfH$distToForest)) # sqrt
# table(udDfH$landuse) # unbalanced, but no solution
# hist(log(udDfH$distToNestingBox + 5)) # log(x + 5)
# hist(udDfH$altitude)# normal
# hist(sqrt(udDfH$tri)) # sqrt
# hist(sqrt(udDfH$distToRoads)) # sqrt

udDfHs <- udDfH
udDfHs$distToForest <- sqrt(udDfHs$distToForest)
udDfHs$distToNestingBox <- log(udDfHs$distToNestingBox + 5)
udDfHs$tri <- sqrt(udDfHs$tri)
udDfHs$distToRoads <- sqrt(udDfHs$distToRoads)

#----------------------------------------------------#
# scale variables
for(i in names(udDfHs)[c(5,7,8,9,10)]) udDfHs[,i] <- scale(udDfHs[,i])

#----------------------------------------------------#
# predictor correlation
# plot(varclus(as.matrix(udDfHs[,c(5,7,8,9,10)]))) # not an issue

#----------------------------------------------------#
# shrink ud 
udDfHs$ud <- (udDfHs$ud - 0.5) / 1.02 + 0.5


#----------------------------------------------------#
# find breakpoints visually: HR, Dispersal, Exclusion

plot(density(udDfHs$ud))
abline(v=c(0.19, 0.67))
plot(udDfH$posX, udDfH$posY, col = ifelse(udDfH$ud > 0.67, 2, ifelse(udDfH$ud > 0.19, 3, 1)))

# #----------------------------------------------------#
# INLA, use combined model BYM, which
## contains both the "iid" and the "besag" model

# define neighborhood # encompass diagonals in 20x20 grid
# neighs <- dnearneigh(as.matrix(udDfH[,c("x", "y")]), d1 = 0.01, d2 = 100) 
# spId <- seq.int(NROW(udDfH))
# nb2INLA("analysis/nc.txt", neighs)

# fmInlaH <- inla(ud ~ landuse + distToForest + I(distToForest^2) 
#                   + distToNestingBox + altitude + distToRoads
#                   + f(spId, model = "besag", graph = "analysis/nc.txt")
#                   + f(id ,model="iid"),
#                data = udDfHs,
#                family = "gaussian", control.family = list(link = "logit"),
#                control.predictor = list(compute = TRUE),
#                control.compute = list(dic = TRUE, cpo = TRUE))

#----------------------------------------------------#
# # GLMM
# fmGlmmH <- glmer(ud ~ landuse + distToForest + I(distToForest^2) 
#                 + distToNestingBox + altitude + I(altitude^2)
#                 + distToRoads #
#                 + (1|id), 
#                 control = glmerControl(optCtrl = list(maxfun = 20000)),
#                 family = gaussian(link = "logit"), data = udDfHs)

#----------------------------------------------------#
# # GAM
system.time(
  fmGamH <- gam(ud ~ landuse + s(distToForest, bs = "cs") 
                  + s(distToNestingBox, bs = "cs") + s(altitude, bs = "cs") 
                  + s(distToRoads, bs = "cs") + s(id, bs = "re"), 
                family = gaussian(link = "logit"), select = TRUE,
                data = udDfHs)
) # 350 sec

summary(fmGamH)



system.time(
  fmGamH2 <- gam(ud ~ landuse + s(distToForest, bs = "ts") 
                + s(distToNestingBox, bs = "ts") + s(altitude, bs = "ts") 
                + s(distToRoads, bs = "ts") + s(id, bs = "re"), 
                family = gaussian(link = "logit"), select = TRUE,
                data = udDfHs)
) # 505.35 sec

summary(fmGamH2)


#----------------------------------------------------#
# randomForest
system.time(
  fmRfH <- randomForest(ud ~ landuse + distToForest
                        + distToNestingBox + altitude 
                        + distToRoads, 
                        data = udDfHs, importance = TRUE, 
                        ntree = 200)
) # 69 sec


# #----------------------------------------------------#
# #----------------------------------------------------#
# #----------------------------------------------------#
# # step length distribution
# stepspeed <- unlist(sapply(seq.int(length(ltdf)), 
#                            function(x) ltdf[[x]]$dist / ltdf[[x]]$dt)) # in m/s
# mean(stepspeed, na.rm = T)*3600 # mean 330 m/h; not really though
# median(stepspeed, na.rm = T)*3600 # median 22 m/h
