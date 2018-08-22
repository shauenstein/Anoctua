# set working directory
# setwd("C:/Users/SH/Dropbox/MSc_LittleOwl/Steinkauz_FR")
# setwd("D:/InteamDocs/Uni/littleowl")

# packages
library(sp)
library(adehabitatLT)
library(hab)
#----------------------------------------------------#

cleaned <- read.csv("data/natal_dispersal_cleaned.csv",
                    stringsAsFactors=F)

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

#----------------------------------------------------#

# obtain subsequent 5-min fixes: keep fixes with time lag == 5 min
ltdf <- subset(ltdf, dt == 300, rec = TRUE) # recompute 


















    


