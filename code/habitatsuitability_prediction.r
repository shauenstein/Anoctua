# predict habitat suitability on template raster
# SH/11/04/2018


#####------##------##------##------##------##------##------##------##------##------#####

library(raster)
library(gdalUtils)
library(mgcv)
library(rgdal)
library(sp)

#####------##------##------##------##------##------##------##------##------##------#####
# load fitted gam object
load("../analysis/trainingGam.RData")
load("../analysis/model_fitting_data.RData")


#####------##------##------##------##------##------##------##------##------##------#####
# template raster
rast <- raster(xmn=3308355, xmx=3610195, ymn=5182816, ymx=5517416, 
               crs = CRS("+proj=tmerc +lat_0=0 +lon_0=9 +k=1 +x_0=3500000 +y_0=0 +ellps=bessel +towgs84=598.1,73.7,418.2,0.202,0.045,-2.455,6.7 +units=m +no_defs"), resolution = 20, vals=NA)

# prediction data frame
predictor_df <- data.frame(dum = rep(0, ncell(rast)), distToNestingBox = 0, id = unique(udDf$id)[1])

#####------##------##------##------##------##------##------##------##------##------#####
# prepare predictor data 
#####------##------##------##------##------#####
# landcover
lc <- raster("../envVars/envOnStudyExtent/landcover_20.tif")
predictor_df$landuse <- as.factor(values(lc))

#####------##------##------##------##------#####
# distance to forest
# forest <- lc
# forest[forest != 6] <- NA
# writeRaster(forest, filename = "../envVars/envOnStudyExtent/forest_20.tif", datatype = "INT1U")
# # do distance calculation with qgis: gdal_proximity.bat
# distToF <- raster("../envVars/envOnStudyExtent/distToF.tif")
# distToF[is.na(lc)] <- NA
# writeRaster(distToF, filename= "../envVars/envOnStudyExtent/forest_distance_20.tif", datatype = "INT2S")
dF <- raster("../envVars/envOnStudyExtent/forest_distance_20.tif")
dF[dF > 2000] <- 2000
predictor_df$distToForest <- (sqrt(values(dF)) - mean(sqrt(udDf$distToForest))) / sd(sqrt(udDf$distToForest))

#####------##------##------##------##------#####
# altitude
# gdalwarp(srcfile = "../envVars/elevproj10m/w001000.adf",
#          dstfile = "../envVars/elevation/elevation_20.tif",
#          s_srs = "+proj=tmerc +lat_0=0 +lon_0=9 +k=1 +x_0=3500000 +y_0=0 +ellps=bessel +units=m +no_defs", t_srs = proj4string(rast),
#          srcnodata = -32768, dstnodata = -999,
#          r = "bilinear", tr = c(20,20),
#          te = c(xmn=3308355, ymn=5182816, xmx=3610195, ymx=5517416),
#          overwrite = TRUE)
# gdalwarp(srcfile = "../envVars/elevation/srtm_severin.tif",
#          dstfile = "../envVars/elevation/elevation_rest_20.tif",
#          s_srs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", t_srs = proj4string(rast),
#          srcnodata = "n/a", dstnodata = -999,
#          r = "bilinear", tr = c(20,20),
#          te = c(xmn=3308355, ymn=5182816, xmx=3610195, ymx=5517416),
#          overwrite = TRUE)
# 
# mosaic_rasters(gdalfile = c("../envVars/elevation/elevation_20.tif",
#                             "../envVars/elevation/elevation_rest_20.tif"),
#                dst_dataset = "../envVars/envOnStudyExtent/elev_20.tif",
#                srcnodata = -999, r = "bilinear")
# intermediate <- raster("../envVars/envOnStudyExtent/elev_20.tif") # load gdalUtils output to save to INT1U
# intermediate[is.na(lc)] <- NA
# writeRaster(intermediate, filename = "../envVars/envOnStudyExtent/elevation_20.tif", datatype = "INT2S")
dem <- raster("../envVars/envOnStudyExtent/elevation_20.tif")
predictor_df$altitude <- (values(dem) - mean(udDf$altitude)) / sd(udDf$altitude)

#####------##------##------##------##------#####
# distance to settlements

# select settlements and merge all layers
# bw_places <-  readOGR(dsn = "../envVars/settlements/bw", layer = "places")
# bw_places <- bw_places[bw_places$type %in% c("hamlet", "suburb", "village"), c("osm_id", "type")]
# colnames(bw_places@data)[2] <- "fclass"
# ch_places <-  readOGR(dsn = "../envVars/settlements/ch", layer = "places")
# ch_places <- ch_places[ch_places$type %in% c("hamlet", "suburb", "village"), c("osm_id", "type")]
# colnames(ch_places@data)[2] <- "fclass"
# al_places <-  readOGR(dsn = "../envVars/settlements/alsace", layer = "gis.osm_places_free_1")
# al_places <- al_places[al_places$fclass %in% c("hamlet", "suburb", "village"), c("osm_id", "fclass")]
# lo_places <-  readOGR(dsn = "../envVars/settlements/lorraine", layer = "gis.osm_places_free_1")
# lo_places <- lo_places[lo_places$fclass %in% c("hamlet", "suburb", "village"), c("osm_id", "fclass")]
# fc_places <-  readOGR(dsn = "../envVars/settlements/franche_comte", layer = "gis.osm_places_free_1")
# fc_places <- fc_places[fc_places$fclass %in% c("hamlet", "suburb", "village"), c("osm_id", "fclass")]
# rp_places <-  readOGR(dsn = "../envVars/settlements/rheinland_pfalz", layer = "gis.osm_places_free_1")
# rp_places <- rp_places[rp_places$fclass %in% c("hamlet", "suburb", "village"), c("osm_id", "fclass")]
# sa_places <-  readOGR(dsn = "../envVars/settlements/saarland", layer = "gis.osm_places_free_1")
# sa_places <- sa_places[sa_places$fclass %in% c("hamlet", "suburb", "village"), c("osm_id", "fclass")]
# 
# settlements <- rbind(bw_places, ch_places, al_places, lo_places, fc_places, rp_places, sa_places)
# settlements <- spTransform(settlements, proj4string(rast))
# settlements <- crop(settlements, rast)
# writeOGR(settlements, dsn = "../envVars/settlements", layer = "settlements", driver = "ESRI Shapefile")

# gdal_rasterize(src_datasource = "../envVars/settlements/settlements.shp",
#                dst_filename = "../envVars/settlements/settlements.tif",
#                burn = 1, a_nodata = 255,
#                te = c(xmn=3308355, ymn=5182816, xmx=3610195, ymx=5517416),
#                tr = c(20,20), ot = "Byte)

# # do distance calculation with qgis: gdal_proximity.bat
# distToS <- raster("../envVars/envOnStudyExtent/distToS.tif")
# distToS[is.na(lc)] <- NA
# writeRaster(distToS, filename= "../envVars/envOnStudyExtent/settlements_distance_20.tif", datatype = "INT2S")
dS <- raster("../envVars/envOnStudyExtent/settlements_distance_20.tif")
dS[dS > 2000] <- 2000
predictor_df$distToSettlements <- (log(values(dS) + 15) - mean(log(udDf$distToSettlements + 15))) / 
  sd(log(udDf$distToSettlements + 15))

#####------##------##------##------##------#####
# distance to major roads

# select major road types and merge all layers
# bw_roads <-  readOGR(dsn = "../envVars/roads/bw", layer = "roads")
# bw_roads <- bw_roads[bw_roads$type %in% c("motorway", "motorway_link", "primary", "primary_link", "secondary", "secondary_link", "tertiary", "tertiary_link"), c("osm_id", "type")]
# colnames(bw_roads@data)[2] <- "fclass"
# ch_roads <-  readOGR(dsn = "../envVars/roads/ch", layer = "roads")
# ch_roads <- ch_roads[ch_roads$type %in% c("motorway", "motorway_link", "primary", "primary_link", "secondary", "secondary_link", "tertiary", "tertiary_link"), c("osm_id", "type")]
# colnames(ch_roads@data)[2] <- "fclass"
# al_roads <-  readOGR(dsn = "../envVars/roads/alsace", layer = "gis.osm_roads_free_1")
# al_roads <- al_roads[al_roads$fclass %in% c("motorway", "motorway_link", "primary", "primary_link", "secondary", "secondary_link", "tertiary", "tertiary_link"), c("osm_id", "fclass")]
# lo_roads <-  readOGR(dsn = "../envVars/roads/lorraine", layer = "gis.osm_roads_free_1")
# lo_roads <- lo_roads[lo_roads$fclass %in% c("motorway", "motorway_link", "primary", "primary_link", "secondary", "secondary_link", "tertiary", "tertiary_link"), c("osm_id", "fclass")]
# fc_roads <-  readOGR(dsn = "../envVars/roads/franche_comte", layer = "gis.osm_roads_free_1")
# fc_roads <- fc_roads[fc_roads$fclass %in% c("motorway", "motorway_link", "primary", "primary_link", "secondary", "secondary_link", "tertiary", "tertiary_link"), c("osm_id", "fclass")]
# rp_roads <-  readOGR(dsn = "../envVars/roads/rheinland_pfalz", layer = "gis.osm_roads_free_1")
# rp_roads <- rp_roads[rp_roads$fclass %in% c("motorway", "motorway_link", "primary", "primary_link", "secondary", "secondary_link", "tertiary", "tertiary_link"), c("osm_id", "fclass")]
# sa_roads <-  readOGR(dsn = "../envVars/roads/saarland", layer = "gis.osm_roads_free_1")
# sa_roads <- sa_roads[sa_roads$fclass %in% c("motorway", "motorway_link", "primary", "primary_link", "secondary", "secondary_link", "tertiary", "tertiary_link"), c("osm_id", "fclass")]
# 
# major_roads <- rbind(bw_roads, ch_roads, al_roads, lo_roads, fc_roads, rp_roads, sa_roads)
# major_roads <- spTransform(major_roads, proj4string(rast))
# major_roads <- crop(major_roads, rast)
# writeOGR(major_roads, dsn = "../envVars/roads", layer = "major_roads", driver = "ESRI Shapefile")

# gdal_rasterize(src_datasource = "../envVars/roads/major_roads.shp",
#                dst_filename = "../envVars/roads/major_roads.tif",
#                burn = 1, a_nodata = 255,
#                te = c(xmn=3308355, ymn=5182816, xmx=3610195, ymx=5517416),
#                tr = c(20,20), ot = "Byte")

# # do distance calculation with qgis: gdal_proximity.bat
# distToR <- raster("../envVars/envOnStudyExtent/distToR.tif")
# distToR[is.na(lc)] <- NA
# writeRaster(distToR, filename= "../envVars/envOnStudyExtent/roads_distance_20.tif", datatype = "INT2S")
dR <- raster("../envVars/envOnStudyExtent/roads_distance_20.tif")
dR[dR > 2000] <- 2000
predictor_df$distToRoads <- (sqrt(values(dR)) - mean(sqrt(udDf$distToRoads))) / sd(sqrt(udDf$distToRoads))

#####------##------##------##------##------##------##------##------##------##------#####

# save(predictor_df, file = "../analysis/prediction_data.RData")

load("../analysis/prediction_data.RData")

#####------##------##------##------##------##------##------##------##------##------#####

ppn <- 100

index_seq <- as.integer(seq(0, NROW(predictor_df), len = ppn+1))

hs <- 1:3
for(i in 1:ppn){
  hs[(index_seq[i]+1):index_seq[i+1]] <- predict(fmGam, 
                                                 newdata = predictor_df[(index_seq[i]+1):index_seq[i+1],], 
                                                 type = "response")
  cat(i)
}


hs_rast <- raster(xmn=3308355, xmx=3610195, ymn=5182816, ymx=5517416, 
                  crs = CRS("+proj=tmerc +lat_0=0 +lon_0=9 +k=1 +x_0=3500000 +y_0=0 +ellps=bessel +towgs84=598.1,73.7,418.2,0.202,0.045,-2.455,6.7 +units=m +no_defs"), resolution = 20, vals=hs)


writeRaster(hs_rast, filename = "../envVars/envOnStudyExtent/habitatSuitability_20.tif", format = "GTiff", datatype = "FLT4S")

hs_rast_scaled <- raster(xmn=3308355, xmx=3610195, ymn=5182816, ymx=5517416, 
                         crs = CRS("+proj=tmerc +lat_0=0 +lon_0=9 +k=1 +x_0=3500000 +y_0=0 +ellps=bessel +towgs84=598.1,73.7,418.2,0.202,0.045,-2.455,6.7 +units=m +no_defs"), resolution = 20, vals= (hs - min(hs, na.rm = T)) / (max(hs, na.rm = T) - min(hs, na.rm = T)))

writeRaster(hs_rast_scaled, filename = "../envVars/envOnStudyExtent/habitatSuitability_scaled_20.tif", format = "GTiff", datatype = "FLT4S")


#####------##------##------##------##------##------##------##------##------##------#####
#####------##------##------##------##------##------##------##------##------##------#####
# plot habitat suitability raster for publication
# load raster
habitatSuitability <- raster("../envVars/envOnStudyExtent/habitatSuitability_scaled_20.tif")

# showMethods("plot")
# getMethod("plot", c("Raster", "ANY"))
# getAnywhere(".plotraster2")
# getAnywhere(".rasterImagePlot")
args(raster:::.rasterImagePlot)

png(filename = "../figures/results/haitatSuitability.png", width = 1000, height = 1050)
par(mar = c(0,0.5,0,0))
plot(habitatSuitability, legend = F, axes = F, box = F)
dev.off()




