## Script to make disparity through time (dtt) plots from the data.

library(geomorph)
library(geiger)
source("./functions/prepare-data.R")
source("./functions/analysis.R")

## Get data. See 'Prepare_data.R':
load("./data/Geniculata_data.RData")

#####################################################
## Get phylogenetic distance and Procrustes distance data:

## Phylogenetic distance:
dist.tr <- cophenetic(tr)

## Pairwise procrustes distances between shapes:
pp.male <- make.proc.dist(cord.ind.male)
pp.female <- make.proc.dist(cord.ind.female)
pp.scu.male <- make.proc.dist(cord.scu.male)
pp.scu.female <- make.proc.dist(cord.scu.female)
pp.pro.male <- make.proc.dist(cord.pro.male)
pp.pro.female <- make.proc.dist(cord.pro.female)
pp.jug.male <- make.proc.dist(cord.jug.male)
pp.jug.female <- make.proc.dist(cord.jug.female)

## Correct order all results:
pp.male <- pp.male[colnames(dist.tr),colnames(dist.tr)]
pp.female <- pp.female[colnames(dist.tr),colnames(dist.tr)]
pp.scu.male <- pp.scu.male[colnames(dist.tr),colnames(dist.tr)]
pp.scu.female <- pp.scu.female[colnames(dist.tr),colnames(dist.tr)]
pp.pro.male <- pp.pro.male[colnames(dist.tr),colnames(dist.tr)]
pp.pro.female <- pp.pro.female[colnames(dist.tr),colnames(dist.tr)]
pp.jug.male <- pp.jug.male[colnames(dist.tr),colnames(dist.tr)]
pp.jug.female <- pp.jug.female[colnames(dist.tr),colnames(dist.tr)]

## Plotting:
#####################################################
par(mfrow = c(2,1))
par(mai = c(0.1,1,1,0.1))
plot(c(dist.tr), c(pp.male), ylim = c(0,0.4), xlim = c(0, 12), xlab = ""
   , ylab = "", axes = F)
abline(lm(c(pp.male) ~ c(dist.tr)), lwd = 2)
points(c(dist.tr), c(pp.scu.male), pch = 4, col = "grey")
abline(lm(c(pp.scu.male) ~ c(dist.tr)), lty = 1, lwd = 2, col = "grey")
points(c(dist.tr), c(pp.pro.male), pch = 4, col = "grey")
abline(lm(c(pp.pro.male) ~ c(dist.tr)), lty = 1, lwd = 2, col = "grey")
points(c(dist.tr), c(pp.jug.male), pch = 4, col = "grey")
abline(lm(c(pp.jug.male) ~ c(dist.tr)), lty = 1, lwd = 2, col = "grey")
axis(side = 2, at = c(0.0, 0.2, 0.4), cex.axis = 0.8)
points(x = 12, y = 0.35, pch = -0x2642L, cex = 3.0)
mtext("Procrustes distance", side = 2, line = 2.5)

par(mai = c(1,1,0.1,0.1))
plot(c(dist.tr), c(pp.female), ylim = c(0,0.4), xlim = c(0,12),xlab = ""
   , ylab = "", axes = F)
abline(lm(c(pp.female) ~ c(dist.tr)), lwd = 2)
points(c(dist.tr), c(pp.scu.female), pch = 4, col = "grey")
abline(lm(c(pp.scu.female) ~ c(dist.tr)), lty = 1, lwd = 2, col = "grey")
points(c(dist.tr), c(pp.pro.female), pch = 4, col = "grey")
abline(lm(c(pp.pro.female) ~ c(dist.tr)), lty = 1, lwd = 2, col = "grey")
points(c(dist.tr), c(pp.jug.female), pch = 4, col = "grey")
abline(lm(c(pp.jug.female) ~ c(dist.tr)), lty = 1, lwd = 2, col = "grey")
axis(side = 2, at = c(0.0, 0.2, 0.4), cex.axis = 0.8)
axis(side = 1, cex.axis = 0.8)
points(x = 12, y = 0.35, pch = -0x2640L, cex = 3.0)
mtext("Procrustes distance", side = 2, line = 2.5)
mtext("Phylogenetic distance", side = 1, line = 2.5)
