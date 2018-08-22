# conceptual figures little owl IBM #
# SH-13/01/17                       #
#-----------------------------------#

#--------------------------------------------------------------------------#
# packages
library(RandomFields)
#library(dismo)
library(ggplot2)
library(raster)

#--------------------------------------------------------------------------#
# overview figure
#--------------------------------------------------------------------------#

# simulate environment 
# size in xcell x ycell
xcell <- 50; ycell <- 50
# lon lat coords
Xvec <- seq.int(xcell) ; Yvec <- seq.int(ycell) 
coords <- expand.grid(X = Xvec, Y = Yvec)
#--------------------------------------------------------------------------#
set.seed(4) # random seed
expCov <- RMexp(var = 1, scale = 40)
env <- as.vector(RFsimulate(expCov, x = Xvec, y = Yvec, spConform=FALSE))
env <- (env - min(env)) / diff(range(env))
par(mar=rep(0.5,4))
# plot(coords$X, coords$Y, pch=15, cex=3,col=terrain.colors(100)[env*100], xlab="", ylab="", yaxt="n")
pathCoords <- data.frame(x = c(29,33,35,33,28,21,13,9,17,26,34,36,39,41)-5,
                         y = c(30,35,42,46,41,40,36,33,35,34,31,25,20,26),
                         roost = c(1,rep(0,2),1,rep(0,3),1,rep(0,4),1,0))
# cut off 5 cells on each side and 10 from the bottom
env <- env[which(coords$X > 5 & coords$X <= 45 & coords$Y > 10)]
coords <- coords[which(coords$X > 5 & coords$X <= 45 & coords$Y > 10),]
coords$X <- coords$X - 5

# lines(pathCoords$x, pathCoords$y, type = "b")
# points(pathCoords$x[pathCoords$roost==1], pathCoords$y[pathCoords$roost==1], col = 2)

pdf("../figures/concept_overview.pdf", width = 10, height = 10)
ggplot(coords, aes(X, Y)) +
  geom_tile(aes(fill=env), colour = "grey70", size = 0.01) +
  scale_fill_gradientn(colours = terrain.colors(100), guide = FALSE) +
  coord_fixed() +
  xlab("") +
  ylab("") +
  geom_path(data = pathCoords, aes(x,y), size = 1) +
  geom_point(data = pathCoords, aes(x,y), colour = ifelse(pathCoords$roost ==1, "red","black"), size = 6) +
  theme(plot.margin = unit(rep(-.7,4), "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.width = element_blank(),
        legend.key.height = element_blank(), 
        legend.title=element_blank(),
        legend.text=element_blank(),
        legend.direction = "vertical",
        strip.text.x = element_blank(), 
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())
dev.off()
#--------------------------------------------------------------------------#
# concept figure
#--------------------------------------------------------------------------#
# grid
xcellc <- 21; ycellc <- 21
Xvecc <- seq(1, xcellc, len = xcellc); Yvecc <- seq(1, ycellc, len = ycellc) 
coordsc <- expand.grid(X = Xvecc, Y = Yvecc)

# location fixes
fixes <- data.frame(x = c(3,6,11),
                    y = c(0.5,4,11))
newfix1 <- data.frame(x = c(11,13),
                      y = c(11,17))

# perception range grid
t1 <- t2 <- rep(0, nrow(coordsc))
for(i in seq(nrow(coordsc))){
  d1 <- sqrt((fixes[2,1]-coordsc$X[i])^2 + (fixes[2,2]-coordsc$Y[i])^2)
  d2 <- sqrt((fixes[3,1]-coordsc$X[i])^2 + (fixes[3,2]-coordsc$Y[i])^2)
  if(d1 <= 9) t1[i] <- 1
  if(d2 <= 9) t2[i] <- 1
  }


# perception range circle
circleFun <- function(center = c(0,0),radius = 1, npoints = 100){
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + radius * cos(tt)
  yy <- center[2] + radius * sin(tt)
  return(data.frame(x = xx, y = yy))
}

pRange1 <- circleFun(c(6,4),9.5,npoints = 200)
# cut 
pRange1$x <- ifelse(pRange1$x < .5, .5, pRange1$x)
pRange1$y <- ifelse(pRange1$y < .5, .5, pRange1$y)
pRange2 <- circleFun(c(11,11),9.5,npoints = 200)

pdf("../figures/concept_concept.pdf", width = 10, height = 10)
ggplot(coordsc, aes(X, Y)) +
  geom_tile(aes(fill=0), colour = "grey70", size = 0.01, fill = "white") +
  geom_tile(aes(fill=t1), colour = "grey70", size = 0.01, fill = ifelse(t1 == 1,"grey90", rgb(0,0,0,0))) +
  geom_tile(aes(fill=t2), colour = "grey70", size = 0.01, fill = ifelse(t2 == 1,rgb(0.545,0,0,0.2), rgb(0,0,0,0))) +
  coord_fixed() +
  xlab("") +
  ylab("") +
  geom_path(data = fixes, aes(x,y)) +
  geom_point(data = fixes, aes(x,y), colour = c(rgb(0,0,0,0),"grey50","black"), 
             size = 5) +
  # potential fixes
  geom_path(data = newfix1, aes(x,y), linetype = 2) +
  geom_point(data = newfix1[2,], aes(x,y), colour = "black", shape = 4, 
             size = 6) +
  geom_path(data = pRange1, aes(x,y), colour = "grey70") +
  geom_path(data = pRange2, aes(x,y)) +
  theme(plot.margin = unit(rep(-.7,4), "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.width = element_blank(),
        legend.key.height = element_blank(), 
        legend.title=element_blank(),
        legend.text=element_blank(),
        legend.direction = "vertical",
        strip.text.x = element_blank(), 
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())
dev.off()
