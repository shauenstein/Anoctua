# merging landcover maps: manually filled gaps, swiss, german and french data
# SH/09/04/2018


#####------##------##------##------##------##------##------##------##------##------#####

library(raster)
library(gdalUtils)

#####------##------##------##------##------##------##------##------##------##------#####
# template raster
rast <- raster(xmn=3308355, xmx=3610195, ymn=5182816, ymx=5517416, 
               crs = CRS("+proj=tmerc +lat_0=0 +lon_0=9 +k=1 +x_0=3500000 +y_0=0 +ellps=bessel +towgs84=598.1,73.7,418.2,0.202,0.045,-2.455,6.7 +units=m +no_defs"), resolution = 20, vals=NA)

#####------##------##------##------##------##------##------##------##------##------#####
# load individual files
# aggregate and crop to template properties
gdalwarp(srcfile = "../envVars/landcover/landcover2007.tif",
         dstfile = "../envVars/landcover/landcover2007_20.tif",
         s_srs = proj4string(rast), t_srs = proj4string(rast),
         srcnodata = 255, dstnodata = 255,
         r = "mode", tr = c(20,20),
         te = c(xmn=3308355, ymn=5182816, xmx=3610195, ymx=5517416),
         overwrite = TRUE)
gdalwarp(srcfile = "../envVars/landcover/OCS_2016_CESBIO.tif",
         dstfile = "../envVars/landcover/OCS_2016_CESBIO_20.tif",
         s_srs = CRS("+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"), t_srs = proj4string(rast),
         srcnodata = "n/a", dstnodata = 255,
         r = "mode", tr = c(20,20), 
         te = c(xmn=3308355, ymn=5182816, xmx=3610195, ymx=5517416),
         overwrite = TRUE)

# harmonise landcover values
lc_chger <- raster("../envVars/landcover/landcover2007_20.tif")
lc_chger[values(lc_chger) == 9] <- 30
lc_chger[values(lc_chger) == 10] <- 11
writeRaster(lc_chger, filename = "../envVars/landcover/lc_chger_20.tif", datatype = "INT1U")


lc_fr <- raster("../envVars/landcover/OCS_2016_CESBIO_20.tif")
lc_fr[values(lc_fr) %in% 11:12] <- 1
lc_fr[values(lc_fr) %in% 31:32] <- 6
lc_fr[values(lc_fr) %in% c(34,211)] <- 30
lc_fr[values(lc_fr) == 36] <- 80
lc_fr[values(lc_fr) %in% 41:44] <- 41
lc_fr[values(lc_fr) %in% c(45,46,53)] <- 25
lc_fr[values(lc_fr) == 51] <- 11
lc_fr[values(lc_fr) %in% 221:222] <- 21
writeRaster(lc_fr, filename = "../envVars/landcover/lc_fr_20.tif", datatype = "INT1U")

lc_chger_gaps <- raster("../envVars/envOnStudyExtent/landuseFilled.tif")
lc_chger_gaps[values(lc_chger_gaps) == 9] <- 30
lc_chger_gaps[values(lc_chger_gaps) == 10] <- 11
writeRaster(lc_chger_gaps, filename = "../envVars/landcover/lc_chger_gaps_20.tif", datatype = "INT1U")


# mosaic landcover maps 
mosaic_rasters(gdalfile = c("../envVars/landcover/lc_chger_20.tif",
                            "../envVars/landcover/lc_fr_20.tif",
                            "../envVars/landcover/lc_chger_gaps_20.tif"),
               dst_dataset = "../envVars/envOnStudyExtent/lc_20.tif",
               srcnodata = 255, r = "mode")
intermediate <- raster("../envVars/envOnStudyExtent/lc_20.tif") # load gdalUtils output to save to INT1U
writeRaster(intermediate, filename = "../envVars/envOnStudyExtent/landcover_20.tif", datatype = "INT1U")
