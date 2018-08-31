# conceptual figures little owl IBM #
# SH-18/01/17                       #
#-----------------------------------#

#--------------------------------------------------------------------------#
# packages
library(ggplot2)
library(raster)
library(sp)
library(rgdal)
library(rgeos)
library(ggsn) # scale bar + north arrow

#--------------------------------------------------------------------------#
# load data
#--------------------------------------------------------------------------#
ringRaw <- read.csv("../data/Steinkauz_Beringungsdaten_ab_2005.csv")
ringData <- data.frame(id = as.factor(ringRaw$STRRINGNR),
                       date = as.POSIXct(strptime(ringRaw$DTMDATE, "%d/%m/%Y %H:%M", 
                                                  tz = "Europe/Berlin")),
                       longitude = ringRaw$LNGLONG,
                       latitude = ringRaw$LNGLAT,
                       geoReg = ringRaw$STRTEXT.2) 
ringSDF <- SpatialPointsDataFrame(coords = ringData[,c("longitude","latitude")], 
                                  data = ringData) 
# add correct projection: wgs84
proj4string(ringSDF) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
# project to GK3
ringSDF <- spTransform(ringSDF, CRS("+init=epsg:31467"))

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

## Load administrative boundary BW (http://www.gadm.org)
bw <- readOGR("../envVars/administrativeBoundary", layer = "DEU_adm1",
              encoding = "ESRI Shapefile", verbose = FALSE)
bw <- bw[bw$NAME_1 == "Baden-Württemberg",]
# project to GK 3
bw <- spTransform(bw, CRS("+init=epsg:31467"))

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

## crop ringing data to bw
ringSDF <- crop(ringSDF, extent(bw))
ringBW <- gIntersects(ringSDF, bw, byid = TRUE)
ringBW <- ringSDF[c(ringBW),]
years <- format(ringBW$date, format="%Y")
# select last 5 years
ringBW <- ringBW[as.numeric(years) >= 2012,]

# rgdal::writeOGR(ringBW, dsn = "../data", layer = "ringBWa2012", driver = "ESRI Shapefile")

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# add Switzerland
ch <- readOGR("../envVars/administrativeBoundary", layer = "CHE_adm0",
              encoding = "ESRI Shapefile", verbose = FALSE)
# project to GK 3
ch <- spTransform(ch, CRS("+init=epsg:31467"))

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# add France
fr <- readOGR("../envVars/administrativeBoundary", layer = "FRA_adm0",
              encoding = "ESRI Shapefile", verbose = FALSE)
# project to GK 3
fr <- spTransform(fr, CRS("+init=epsg:31467"))

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# add study extent
rast <- raster(xmn=3308355, xmx=3610195, ymn=5182816, ymx=5517416, 
               crs = CRS("+proj=tmerc +lat_0=0 +lon_0=9 +k=1 +x_0=3500000 +y_0=0 +ellps=bessel +towgs84=598.1,73.7,418.2,0.202,0.045,-2.455,6.7 +units=m +no_defs"), resolution = 20, vals=NA)
SExt <- extent(rast)
chInS <- crop(ch, SExt)
frInS <- crop(fr, SExt)
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# add primary study extent
cleaned <- read.csv("../data/natal_dispersal_cleaned.csv",
                    stringsAsFactors = FALSE)
spdf <- SpatialPointsDataFrame(coords = cleaned[,c("posX","posY")], 
                               data = cleaned) 
# add correct projection: DHDN / Gauss-Kruger zone 3
proj4string(spdf) <- CRS("+init=epsg:31467")
# primary study extent
pSExt <- extent(spdf)  + rep(c(-10000,10000), 2)
pSBox <- data.frame(x = pSExt[c(1,2,2,1,1)], y = pSExt[c(3,3,4,4,3)])

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# get extent of bw and ch for plot extent
# add Switzerland
bwch <- readOGR("../envVars/administrativeBoundary", layer = "BW_CH_adm0",
                encoding = "ESRI Shapefile", verbose = FALSE)

fr <- readOGR("../envVars/administrativeBoundary", layer = "FRA_adm0",
              encoding = "ESRI Shapefile", verbose = FALSE)
fr@data <- data.frame(ID = 3)
bwchfr <- rbind(bwch, fr)
# project to GK 3
fr <- spTransform(fr, CRS("+init=epsg:31467"))
bwch <- spTransform(bwch, CRS("+init=epsg:31467"))
bwchfr <- spTransform(bwchfr, CRS("+init=epsg:31467"))
# crop to extended SA
bwchfr <- crop(bwchfr, rast)

extBwch <- extent(bwch)
boxBwch <- extent(bwch) + c(-10000,10000,-5000,25000)
SBox <- data.frame(x = boxBwch[c(1,2,2,1,1)], y = boxBwch[c(3,3,4,4,3)])
fr <- crop(fr, SBox)
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# load other topographic features

# cities 
cts <- readOGR("../envVars/topography/cities.kml", layer = "cities")
cts <- spTransform(cts, CRS("+init=epsg:31467"))

# Rhine and Neckar
rvs <- readOGR("../envVars/topography", layer = "Gewaessernetz (AWGN)_line",
               encoding = "ESRI Shapefile", verbose = FALSE)

# regions around Neckarbecken
region <-  readOGR("../envVars/topography", layer = "region",
                   encoding = "ESRI Shapefile", verbose = FALSE)

# Swabian Alb
schwAlb <- readOGR("../envVars/topography", layer = "schwaebische_alb_dissolve",
                   encoding = "ESRI Shapefile", verbose = FALSE)

# Black Forest
bf <- readOGR("../envVars/topography", layer = "schwarzwald_dissolve",
              encoding = "ESRI Shapefile", verbose = FALSE)

# sketched polygon for swiss mountains
mtsCH <- readOGR("../envVars/topography/mountainsCH.kml", layer = "mountainsCH")
mtsCH <- spTransform(mtsCH, CRS("+init=epsg:31467"))
mtsCH <- crop(as(mtsCH, "SpatialLines"), bwchfr)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# digital elevation model as background
dem <- raster("../envVars/envOnStudyExtent/elevation_rprj500.tif")
slo <- terrain(dem)
asp <- terrain(dem, opt = "aspect")
hills <- hillShade(slo, asp, angle = 40, direction = 270)

# dem <- projectRaster(dem, rast)

dem_spdf <- as(dem, "SpatialPixelsDataFrame")
dem_df <- as.data.frame(dem_spdf)
colnames(dem_df) <- c("demv", "x", "y")

hills_spdf <- as(hills, "SpatialPixelsDataFrame")
hills_df <- as.data.frame(hills_spdf)
colnames(hills_df) <- c("hillv", "x", "y")


# scalebar function 
# http://stackoverflow.com/questions/39067838/parsimonious-way-to-add-north-arrow-and-scale-bar-to-ggmap
scalebar = function(x,y,w,n,d, units="km"){
  # x,y = lower left coordinate of bar
  # w = width of bar
  # n = number of divisions on bar
  # d = distance along each division
  
  bar = data.frame( 
    xmin = seq(0.0, n*d, by=d) + x,
    xmax = seq(0.0, n*d, by=d) + x + d,
    ymin = y,
    ymax = y+w,
    z = rep(c(1,0),n)[1:(n+1)],
    fill.col = rep(c("black","white"),n)[1:(n+1)])
  
  labs = data.frame(
    xlab = c(seq(0.0, (n+1)*d, by=d) + x, x), 
    ylab = c(rep(y-w*1.5, n+2), y-3.5*w),
    text = c(as.character(seq(0.0, (n+1)*d/1000, by=d/1000)), units)
  )
  list(bar, labs)
}

#--------------------------------------------------------------------------#
#### wider study area ####
#--------------------------------------------------------------------------#

wsb = scalebar(x = extBwch[1]+0.63*diff(extBwch[1:2]), 
              y = extBwch[3]+0.08*diff(extBwch[3:4]), 
              w = 5000, n = 4, d = 25000)

# # transform to WGS84 epsg:4326
# bw <- spTransform(bw, CRS("+init=epsg:4326"))
# chInS <- spTransform(chInS, CRS("+init=epsg:4326"))
# schwAlb <- spTransform(schwAlb, CRS("+init=epsg:4326"))
# bf <- spTransform(bf, CRS("+init=epsg:4326"))
# mtsCH <- spTransform(mtsCH, CRS("+init=epsg:4326"))
# ch <- spTransform(ch, CRS("+init=epsg:4326"))
# cts <- spTransform(cts, CRS("+init=epsg:4326"))
# rvs <- spTransform(rvs, CRS("+init=epsg:4326"))
# ringBW <- spTransform(ringBW, CRS("+init=epsg:4326"))
# bwch <- spTransform(bwch, CRS("+init=epsg:4326"))
# extBwch <- extent(bwch)

# png("../figures/widerStudyArea.png", width = 1000, height = 1000*diff(extBwch[3:4]) / diff(extBwch[1:2]))
#x11(7,7*diff(extBwch[3:4]) / diff(extBwch[1:2]))
ggplot() +
  #xlim(3200000 ,  extent(bwch)[2]) +
  #ylim(extent(bw)[3],  extent(bw)[4]) + 
  xlab("") + 
  ylab("") +
  # geom_tile(data=hills_df, aes(x=x, y=y, col=hillv)) +
  # scale_colour_gradient(low = "white", high = "black", space = "Lab",
  #                       na.value = NA, guide = FALSE) +
  # geom_tile(data=dem_df, aes(x=x, y=y, fill=demv, alpha = 0.01)) + # need to make pixel smaller
  # scale_fill_gradientn(colours = terrain.colors(5), limits = range(dem_df$demv)) +
  # geom_tile(data=hills_df, aes(x,y,alpha=hillv), fill = "grey20") +
  # scale_alpha(range = c(0, 0.8)) +
  # pbox
  geom_path(data = pSBox, aes(x = x, y = y)) +
  # BW+CH+FR background fill
  geom_polygon(data = as.data.frame(bw@polygons[[1]]@Polygons[[3]]@coords),
               aes(x = V1, y = V2), fill = "grey90") +
  geom_polygon(data = as.data.frame(chInS@polygons[[1]]@Polygons[[2]]@coords), #v
               aes(x = x, y = y), fill = "grey90") +
  geom_polygon(data = as.data.frame(frInS@polygons[[1]]@Polygons[[1]]@coords), #v
               aes(x = x, y = y), fill = "grey90") +
  # regional landmarks
  geom_path(data = as.data.frame(schwAlb@polygons[[1]]@Polygons[[1]]@coords),
            aes(x = V1, y = V2), colour = "darkgrey") +
  geom_path(data = as.data.frame(bf@polygons[[1]]@Polygons[[1]]@coords),
            aes(x = V1, y = V2),colour = "darkgrey") +
  geom_path(data = as.data.frame(mtsCH@lines[[1]]@Lines[[1]]@coords), # v
            aes(x = x, y = y), colour = "darkgrey") +
  geom_path(data = as.data.frame(mtsCH@lines[[1]]@Lines[[2]]@coords), #v
            aes(x = x, y = y), colour = "darkgrey") +
  # BW+CH+FR outline
  geom_path(data = as.data.frame(ch@polygons[[1]]@Polygons[[4]]@coords),
            aes(x = V1, y = V2), colour ="grey80", size = 1) +
  geom_path(data =  as.data.frame(bw@polygons[[1]]@Polygons[[3]]@coords),
            aes(x = V1, y = V2), colour ="grey80", size = 1) +
  geom_path(data = as.data.frame(chInS@polygons[[1]]@Polygons[[2]]@coords), #v
            aes(x = x, y = y), colour ="grey80", size = 1) +
  geom_path(data = as.data.frame(fr@polygons[[1]]@Polygons[[1]]@coords), #v
            aes(x = x, y = y), colour ="grey80", size = 1) +
  geom_point(data = as.data.frame(cts@coords),
             aes(x = coords.x1, y = coords.x2), colour = "grey60", size = 6) +
  geom_path(data = as.data.frame(rvs@lines[[2]]@Lines[[1]]@coords),
            aes(x = V1, y = V2), colour = "black") +
  geom_path(data = as.data.frame(rvs@lines[[1]]@Lines[[1]]@coords),
            aes(x = V1, y = V2), colour = "black") +
  geom_point(data = as.data.frame(ringBW@coords),
             aes(x = longitude, y = latitude), size = 0.3) +
  geom_density_2d(data = as.data.frame(ringBW@coords),
                  aes(x = longitude, y = latitude),colour="grey50", size = 0.3) +
  # box
  geom_path(data = SBox, aes(x = x, y = y)) +
  # pbox
  geom_path(data = pSBox, aes(x = x, y = y)) +
  coord_fixed() +
  # annotate("text", x = 3460000, y = 5230000, label = "Northern") +
  # annotate("text", x = 3460000, y = 5220000, label = "Switzerland") +
  # annotate("text", x = 3490000, y = 5370000, label = "Baden-") +
  # annotate("text", x = 3490000, y = 5360000, label = "W\374rrtemberg") +
  # theme_bw()
  # scalebar(location = "bottomright", dist = 50, dd2km = TRUE,
  #          x.min = extBwch[1], x.max=extBwch[2], y.min=extBwch[3], y.max=extBwch[2],
  #          model = "WGS84") +
  # North arrow
  # geom_segment(arrow=arrow(length=unit(4,"mm")), aes(x=extBwch[1]+0.95*diff(extBwch[1:2]),
  #                                                    xend=extBwch[1]+0.95*diff(extBwch[1:2]),
  #                                                    y=extBwch[3]+0.95*diff(extBwch[3:4]),
  #                                                    yend=extBwch[3]+0.95*diff(extBwch[3:4])+30000),
  #              colour="black", size = 1.5) +
  # annotate(x=extBwch[1]+0.95*diff(extBwch[1:2]),
  #          y=extBwch[3]+0.95*diff(extBwch[3:4])-8000, label="N", colour="black", geom="text", size=7) +
  # scalebar
  geom_rect(data=wsb[[1]], aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=z), inherit.aes=F,
            show.legend = F,  color = "black", fill = wsb[[1]]$fill.col) +
  geom_text(data=wsb[[2]], aes(x=xlab, y=ylab, label=text), size = 10, inherit.aes=F, show.legend = F) +
  theme(plot.margin = unit(c(-2,-1.5,-2,-1.5), "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # legend.key.width = element_blank(),
        # legend.key.height = element_blank(),
        legend.title=element_blank(),
        legend.text=element_blank(),
        legend.direction = "vertical",
        strip.text.x = element_blank(),
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())
# dev.off()


#--------------------------------------------------------------------------#
#### data collection area ####
#--------------------------------------------------------------------------#

psb = scalebar(x = pSExt[1]+0.63*diff(pSExt[1:2]), 
              y = pSExt[3]+0.08*diff(pSExt[3:4]), 
              w = 1000, n = 4, d = 5000)


png("../figures/primaryStudyArea.png", width = 1000, height = 1000*diff(pSExt[3:4]) / diff(pSExt[1:2]))
# x11(7,7*diff(pSExt[3:4]) / diff(pSExt[1:2]))
ggplot() +
  xlim(range(pSBox$x)) +
  ylim(range(pSBox$y)) + 
  xlab("") + 
  ylab("") +
  coord_fixed() + 
  geom_path(data = pSBox, aes(x = x, y = y)) + 
  geom_path(data = as.data.frame(region@polygons[[2]]@Polygons[[1]]@coords), aes(x = V1, y = V2), colour = "darkgrey") +  
  geom_path(data = as.data.frame(region@polygons[[3]]@Polygons[[1]]@coords), aes(x = V1, y = V2), colour = "darkgrey") +  
  geom_path(data = as.data.frame(region@polygons[[4]]@Polygons[[1]]@coords), aes(x = V1, y = V2), colour = "darkgrey") +  
  geom_path(data = as.data.frame(region@polygons[[5]]@Polygons[[1]]@coords), aes(x = V1, y = V2), colour = "darkgrey") +  
  geom_path(data = as.data.frame(region@polygons[[6]]@Polygons[[1]]@coords), aes(x = V1, y = V2), colour = "darkgrey") +  
  geom_path(data = as.data.frame(region@polygons[[7]]@Polygons[[1]]@coords), aes(x = V1, y = V2), colour = "darkgrey") +  
  geom_path(data = as.data.frame(rvs@lines[[1]]@Lines[[1]]@coords), aes(x = V1, y = V2), colour = "black", size = 1) +  
  geom_point(data = as.data.frame(spdf), aes(x = posX, y = posY), size = 0.8, alpha = 0.1) +
  # North arrow
  # geom_segment(arrow=arrow(length=unit(4,"mm")), aes(x=pSExt[1]+0.95*diff(pSExt[1:2]),
  #                                                    xend=pSExt[1]+0.95*diff(pSExt[1:2]),
  #                                                    y=pSExt[3]+0.85*diff(pSExt[3:4]),
  #                                                    yend=pSExt[3]+0.95*diff(pSExt[3:4])+1000), 
  #              colour="black", size = 1.5) +
  # annotate(x=pSExt[1]+0.95*diff(pSExt[1:2]), 
  #          y=pSExt[3]+0.95*diff(pSExt[3:4])-8000, label="N", colour="black", geom="text", size=7) +
  # scalebar
  geom_rect(data=psb[[1]], aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=z), inherit.aes=F,
            show.legend = F,  color = "black", fill = psb[[1]]$fill.col) +
  geom_text(data=psb[[2]], aes(x=xlab, y=ylab, label=text), size = 10, inherit.aes=F, show.legend = F) +
  # theme
  theme(plot.margin = unit(c(-1.2,-1.2,-1.2,-1.2), "cm"),
        panel.margin = unit(c(0,0,0,0),"cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.text.x = element_blank(), 
        strip.background = element_blank(),
        # axis.ticks = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank())
dev.off()
