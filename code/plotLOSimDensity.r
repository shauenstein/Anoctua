plotLOSimDensity <- function(LOSimObj,
                             densityBins = 5,
                             basemap, 
                             bmrange = c(0,1), 
                             bmColorRange = c("grey94","#d7301f","#190000"),
                             maxPixels = 1000000,
                             legendTitlePosition = "left", 
                             legendTitleVJust = 0.4, 
                             legendTitleHJust = 0.6,
                             legendScaleCol = "black",
                             ggplotTheme = theme(plot.margin = unit(c(0,0,0,0), "cm"),
                                                 panel.margin = unit(0,"null"),
                                                 axis.title.x = element_text(size=26),
                                                 axis.title.y = element_text(size=26, angle=90),
                                                 axis.text.x = element_blank(),
                                                 axis.text.y = element_blank(),
                                                 legend.background = element_blank(),
                                                 panel.grid.major = element_blank(),
                                                 panel.grid.minor = element_blank(),
                                                 legend.position = c(0.87, 0.75),
                                                 legend.key.width = unit(0.4, "cm"),
                                                 legend.key.height = unit(0.8, "cm"), 
                                                 legend.title=element_text(size=10, angle = 90),
                                                 legend.text=element_text(size=12),
                                                 legend.direction = "vertical",
                                                 strip.text.x = element_text(size = 30), 
                                                 strip.background = element_blank(),
                                                 axis.ticks = element_blank(),
                                                 panel.background = element_blank())
                             ){
  
  library(raster)
  library(ggplot2)
  
  # reduce simulations to one data.frame
  simulations <- vector("list", length(LOSimObj$simulations))
  for(g in seq_along(LOSimObj$simulations)){
    simulations[[g]] <- do.call(rbind, LOSimObj$simulations[[g]])
  }
  simulations <- do.call(rbind, simulations)
  simulations <- as.data.frame(simulations[,c("x","y")])
  rm(LOSimObj)
  
  # must be raster
  if(class(basemap) != "RasterLayer") stop("Provide basemap of class raster.")
  # plot base map  
  bmExt <- extent(basemap)
  #   Convert raster to dataframes
  # bm <- rasterToPoints(basemap)
  # bm <- data.frame(bm)
  # colnames(bm) <- c("x","y","z")
  # 
  # plot with ggplot
  gplot(bm, maxpixels=maxPixels) + 
  geom_tile(aes(fill=value)) +
  scale_fill_gradient2("Predicted habitat suitability", limits = bmrange, 
                       low = bmColorRange[1], mid = bmColorRange[2],high = bmColorRange[3], 
                       midpoint = sum(bmrange)*0.5, guide = guide_colourbar(ticks = FALSE, 
                                                                            title.position = legendTitlePosition, 
                                                                            title.vjust = legendTitleVJust, 
                                                                            title.hjust = legendTitleHJust,
                                                                            scale_colour = legendScaleCol)) +
  coord_equal() + 
  xlab("") +
  ylab("") +
  stat_density2d(data = simulations, aes(x,y), bins = densityBins) +
  ggplotTheme 
}