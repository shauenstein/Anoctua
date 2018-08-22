# checking and cleaning the data

# set working directory
# setwd("C:/Users/SH/Dropbox/MSc_LittleOwl/Steinkauz_FR")
# setwd("D:/InteamDocs/Uni/littleowl")

# packages


#----------------------------------------------------#

# load raw data
raw <- read.csv("../raw/MovementData/natal_dispersal/natal_dispersal.csv")

#----------------------------------------------------#

# check individual patterns

# windows(10,10)
# plot(raw$X, raw$Y, type = "n")
# 
# for(i in unique(raw$individual)){
#   points(raw$X[raw$individual == i], raw$Y[raw$individual == i], col = 2)
#   Sys.sleep(0.5)
#   points(raw$X[raw$individual == i], raw$Y[raw$individual == i])
# }

#----------------------------------------------------#
# get timestamp from "datum" and "zeit"
# timestamp <- as.POSIXct(strptime(paste(raw$datum, raw$zeit), 
#                                  format="%Y/%m/%d %H:%M", 
#                                  tz = "Europe/Berlin"))


#----------------------------------------------------#
# start cleaning/ preparing
#----------------------------------------------------#

#Individuals to exclude
toExclude <- c("2010GB085.1") # Hit by train and then released again
raw <- raw[!(raw$individual %in% toExclude),]

# bring date and time to a readable format
date <- gsub(" 0:00:00", "", raw$datum)
date <- paste(date, raw$zeit)
date <- as.POSIXct(strptime(date, "%Y/%m/%d %H:%M", 
                            tz = "Europe/Berlin"))

out <- data.frame(id = raw$individual,
                  date= date,
                  posX = raw$X,
                  posY = raw$Y,
                  accuracy = raw$accuracy)

# Reorder by id and date
out <- out[order(out$id),]
ids <- unique(out$id)
for(i in ids){
  subset <- out[out$id == i,]
  subset <- subset[order(subset$date),]
  out[out$id == i,] <- subset
}

#Positions to exclude
#Fixes 18134 and 18135 are recorded at the same time:
which(diff(out$date) == 0)
out[18275,]; out[18276,]
# Delete row 18275, which has slightly worse accuracy
out <- out[-18275,]
# clean date: no doubles
which(diff(out$date) == 0)

#################################################################################

#Save relevant columns
write.csv(out, file="data/natal_dispersal_cleaned.csv",row.names=F)


