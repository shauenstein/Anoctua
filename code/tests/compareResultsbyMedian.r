# corrleation of posterior medians
load("../data/posterior100_3.RData")
posterior100_2 <- posterior100
load("../data/posterior100.RData")

medians <- matrix(unlist(posterior100$median), ncol = 4, byrow = T)
medians_2 <- matrix(unlist(posterior100_2$median), ncol = 4, byrow = T)

cor(medians[,1],medians_2[,1])
cor(medians[,2],medians_2[,2])
cor(medians[,3],medians_2[,3])
cor(medians[,4],medians_2[,4])
