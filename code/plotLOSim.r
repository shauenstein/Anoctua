# SH 07-10-2016
#
# plot LOSim
#
#-#-#-#-#-#-#-#-#-#-#

plotLOSim <- function(LOSimObj, density = FALSE, plotBG = FALSE, plotHS = "used", CRS = NULL, ...){
  if(plotHS == "used"){
    library(raster)
    HSrast <- raster(LOSimObj$environment$values, 
                     xmn = LOSimObj$environment$lowerleft[1],
                     xmx = LOSimObj$environment$lowerleft[1] + 
                       NCOL(LOSimObj$environment$values)*
                       LOSimObj$environment$resolution,
                     ymn = LOSimObj$environment$lowerleft[3],
                     ymx = LOSimObj$environment$lowerleft[3] + 
                       NROW(LOSimObj$environment$values)*
                       LOSimObj$environment$resolution,
                     crs = CRS)
    plot(HSrast, ...)
  }
  if(plotHS == "empty"){
    # empty plot
    plot(1:3, type = "n", xlim = c(LOSimObj$environment$lowerleft[1],
                                   LOSimObj$environment$lowerleft[1] + 
                                     NCOL(LOSimObj$environment$values)*
                                     LOSimObj$environment$resolution),
         ylim = c(LOSimObj$environment$lowerleft[3],
                  LOSimObj$environment$lowerleft[3] + 
                    NROW(LOSimObj$environment$values)*
                    LOSimObj$environment$resolution), 
         ...)
  }
  
  # add lines
  for(g in seq_along(LOSimObj$simulations)){
    for(i in seq_along(LOSimObj$simulations[[g]])){
      lines(LOSimObj$simulations[[g]][[i]][,c("x","y")],...)
      if(plotBG) points(LOSimObj$simulations[[g]][[i]][1,"x"],
                        LOSimObj$simulations[[g]][[i]][1,"y"], 
                        col = 2, pch = 15, cex = 0.5)
    }
  }
}


