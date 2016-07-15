## This script reproduces all plots reported in the manuscript.
## Please refer to "Make_analysis.R" to perform the analysis and produce the objects that the
##       plots are based on. As an alternative you can just load the data with results of previous
##       simulations. We encourage the reader to run the analysis again for reproducibility.

## Uncomment to load previous simulations and results:
library(ape)
library(geomorph)
load("./data/all_analyses_results.RData")

## Figure 1: Not made in R.

######################################################################################
## Figure 2: Correlated evolutionary changes of shape.

## How many axes we need to explain the covariation of the genitalia?
## We need two PLS1 axes to explain 95% of the evolutionary covariation in the data.
int.tr$pls.singular.values[1]^2 / sum(int.tr$pls.singular.values^2) ## 79 % by the first.
int.tr$pls.singular.values[2]^2 / sum(int.tr$pls.singular.values^2) ## 16 % by the second, cumulative of 95%.

int.grafen$pls.singular.values[1]^2 / sum(int.grafen$pls.singular.values^2) ## 73 % by the first.
int.grafen$pls.singular.values[2]^2 / sum(int.grafen$pls.singular.values^2) ## 20 % by the second, cumulative of 93%.

## First plot the first and second axes of covariation in the data. Colors will be the same in the plots and in the phylogeny.
pdf("Figure_2_PLS_of_male_female_genitalia_unit.pdf", width = 7, height = 14)
par( mfrow = c(2,1) )
plot(x=int.tr$XScores[,1], y=int.tr$YScores[,1], xlab="Phylogenetic PLS1 of male genitalia", ylab="Phylogenetic PLS1 of female genitalia"
   , xlim=c(-0.15, +0.15),ylim=c(-0.1,0.1), type="n")
colfunc <- colorRampPalette(c("blue", "red"))
colp <- vector()
colp[order(int.tr$XScores[,1])] <- colfunc(13)
abline(h = 0, col="slategray2"); abline(v = 0, col="slategray2")
points(x=int.tr$XScores[,1], y=int.tr$YScores[,1], pch=19, col=colp)
adjust <- rep(0, times=13)
adjust[8] <- 0.002
text(x=int.tr$XScores[,1], y=int.tr$YScores[,1]+adjust, labels=gsub("C. ", "", rownames(int.tr$A1.matrix)), pos=4, )

plot(x=int.tr$XScores[,2], y=int.tr$YScores[,2], xlab="Phylogenetic PLS2 of male genitalia", ylab="Phylogenetic PLS2 of female genitalia"
     , type="n", xlim=c(-0.1, 0.1), ylim=c(-0.06,0.06))
colfunc <- colorRampPalette(c("darkgreen", "orange"))
colp <- vector()
colp[order(int.tr$XScores[,2])] <- colfunc(13)
adjust.x <- adjust.y <- rep(0, times=13)
adjust.x[11] <- - 0.033
adjust.x[8] <- - 0.037
adjust.y[1] <- -0.002
adjust.y[7] <- -0.001
adjust.y[12] <- +0.001
abline(h = 0, col="slategray2"); abline(v = 0, col="slategray2")
points(x=int.tr$XScores[,2], y=int.tr$YScores[,2], pch=19, col=colp)
text(x=int.tr$XScores[,2]+adjust.x, y=int.tr$YScores[,2]+adjust.y, labels=gsub("C. ", "", rownames(int.tr$A1.matrix)), pos=4, )

dev.off()

## Now we will plot the extreme shapes of the PLS axes. This will show which are the features of the shape that are related with the evolutionary correlation.

## PLS1 axis:
pdf("PLS1_extreme_shapes.pdf", width = 7, height = 7+(7/2))
par(mfrow=c(3,2))
m <- cord.gen.male[,,"C. rideri"]
lim <- range(c(m[,1], m[,2]))
plot(x=m[,1], y=m[,2], type="c", axes=FALSE, xlab = "", ylab = "", xlim=lim, ylim=lim)
points(x=m[,1], y=m[,2], type="p", pch=16, col="red")
m <- cord.gen.female[,,"C. rideri"]
lim <- range(c(m[,1], m[,2]))
plot(x=m[,1], y=m[,2], type="c", axes=FALSE, xlab = "", ylab = "", xlim=lim, ylim=lim)
points(x=m[,1], y=m[,2], type="p", pch=16, col="red")
m <- cord.gen.male[,,"C. nigritarsis"]
lim <- range(c(m[,1], m[,2]))
plot(x=m[,1], y=m[,2], type="c", axes=FALSE, xlab = "", ylab = "", xlim=lim, ylim=lim)
points(x=m[,1], y=m[,2], type="p", pch=16, col="#6A7F00")
m <- cord.gen.female[,,"C. nigritarsis"]
lim <- range(c(m[,1], m[,2]))
plot(x=m[,1], y=m[,2], type="c", axes=FALSE, xlab = "", ylab = "", xlim=lim, ylim=lim)
points(x=m[,1], y=m[,2], type="p", pch=16, col="#6A7F00")
m <- cord.gen.male[,,"C. immaculata"]
lim <- range(c(m[,1], m[,2]))
plot(x=m[,1], y=m[,2], type="c", axes=FALSE, xlab = "", ylab = "", xlim=lim, ylim=lim)
points(x=m[,1], y=m[,2], type="p", pch=16, col="blue")
m <- cord.gen.female[,,"C. immaculata"]
lim <- range(c(m[,1], m[,2]))
plot(x=m[,1], y=m[,2], type="c", axes=FALSE, xlab = "", ylab = "", xlim=lim, ylim=lim)
points(x=m[,1], y=m[,2], type="p", pch=16, col="blue")
dev.off()

## PLS2 axis:
pdf("PLS2_extreme_shapes.pdf", width = 7, height = 7+(7/2))
par(mfrow=c(3,2))
m <- cord.gen.male[,,"C. gravis"]
lim <- range(c(m[,1], m[,2]))
plot(x=m[,1], y=m[,2], type="c", axes=FALSE, xlab = "", ylab = "", xlim=lim, ylim=lim)
points(x=m[,1], y=m[,2], type="p", pch=16, col="darkgreen")
m <- cord.gen.female[,,"C. gravis"]
lim <- range(c(m[,1], m[,2]))
plot(x=m[,1], y=m[,2], type="c", axes=FALSE, xlab = "", ylab = "", xlim=lim, ylim=lim)
points(x=m[,1], y=m[,2], type="p", pch=16, col="darkgreen")
m <- cord.gen.male[,,"C. pengue"]
lim <- range(c(m[,1], m[,2]))
plot(x=m[,1], y=m[,2], type="c", axes=FALSE, xlab = "", ylab = "", xlim=lim, ylim=lim)
points(x=m[,1], y=m[,2], type="p", pch=16, col="#D49A00")
m <- cord.gen.female[,,"C. pengue"]
lim <- range(c(m[,1], m[,2]))
plot(x=m[,1], y=m[,2], type="c", axes=FALSE, xlab = "", ylab = "", xlim=lim, ylim=lim)
points(x=m[,1], y=m[,2], type="p", pch=16, col="#D49A00")
m <- cord.gen.male[,,"C. vanduzeei"]
lim <- range(c(m[,1], m[,2]))
plot(x=m[,1], y=m[,2], type="c", axes=FALSE, xlab = "", ylab = "", xlim=lim, ylim=lim)
points(x=m[,1], y=m[,2], type="p", pch=16, col="orange")
m <- cord.gen.female[,,"C. vanduzeei"]
lim <- range(c(m[,1], m[,2]))
plot(x=m[,1], y=m[,2], type="c", axes=FALSE, xlab = "", ylab = "", xlim=lim, ylim=lim)
points(x=m[,1], y=m[,2], type="p", pch=16, col="orange")
dev.off()

## To show if there is any difference with the set of grafen branch lengths:
pdf("Figure_2_PLS_of_male_female_genitalia_grafen.pdf", width = 7, height = 14)
par( mfrow = c(2,1) )
plot(x=int.grafen$XScores[,1], y=int.grafen$YScores[,1], xlab="Phylogenetic PLS1 of male genitalia", ylab="Phylogenetic PLS1 of female genitalia"
   , xlim=c(-0.4, +0.4),ylim=c(-0.4,0.4), type="n")
colfunc <- colorRampPalette(c("blue", "red"))
colp <- vector()
colp[order(int.grafen$XScores[,1])] <- colfunc(13)
abline(h = 0, col="slategray2"); abline(v = 0, col="slategray2")
points(x=int.grafen$XScores[,1], y=int.grafen$YScores[,1], pch=19, col=colp)
adjust.x <- rep(0, times=13)
adjust.y <- rep(0, times=13)
adjust.x[7] <- 0.18 -.03
adjust.x[13] <- 0.16 -.03
adjust.x[9] <- 0.14 -.03
adjust.x[1] <- 0.25 -.03
adjust.y[c(10)] <- -0.027
adjust.y[c(9)] <- -0.033
text(x=int.grafen$XScores[,1]+adjust.x, y=int.grafen$YScores[,1]+adjust.y
   , labels=gsub("C. ", "", rownames(int.grafen$A1.matrix)), pos=2, )

plot(x=int.grafen$XScores[,2], y=int.grafen$YScores[,2], xlab="Phylogenetic PLS2 of male genitalia", ylab="Phylogenetic PLS2 of female genitalia"
     , type="n", xlim=c(-0.3, 0.3), ylim=c(-0.2,0.2))
colfunc <- colorRampPalette(c("darkgreen", "orange"))
colp <- vector()
colp[order(int.grafen$XScores[,2])] <- colfunc(13)
adjust.x <- adjust.y <- rep(0, times=13)
adjust.x[3] <- 0.10 -.02
adjust.x[6] <- 0.10 -.02
adjust.x[1] <- 0.19 -.02
adjust.x[12] <- 0.12 -.02
abline(h = 0, col="slategray2"); abline(v = 0, col="slategray2")
points(x=int.grafen$XScores[,2], y=int.grafen$YScores[,2], pch=19, col=colp)
text(x=int.grafen$XScores[,2]+adjust.x, y=int.grafen$YScores[,2]+adjust.y
   , labels=gsub("C. ", "", rownames(int.grafen$A1.matrix)), pos=2, )

dev.off()

######################################################################################
## Figure 3: Results of the analysis of modularity using the CR coefficient.

pdf("CR_analysis_both_trees.pdf", width = 14, height = 7)
par(mfrow = c(1,2))
hist(mod.unit.test$random.CR, xlim = c(0,1.5), ylim=c(0,2000), col="slategray2", border = "white", xlab="CR values"
   , main="", breaks = 50)
hist(mod.unit.test$CR.boot, col="orangered", border = "white", xlab="CR values", breaks = 27, add=TRUE)
lines(x = c(mod.unit.test$CR,mod.unit.test$CR), y = c(0,1100), lwd = 2)
text(x = mod.unit.test$CR, y = 1200, labels = round(mod.unit.test$CR, digits = 2), cex=2)

hist(mod.grafen.test$random.CR, xlim = c(0,1.5), ylim=c(0,2000), col="slategray2", border = "white", xlab="CR values"
   , main="", breaks = 50)
hist(mod.grafen.test$CR.boot, col="orangered", border = "white", xlab="CR values", breaks = 27, add=TRUE)
lines(x = c(mod.grafen.test$CR, mod.grafen.test$CR), y = c(0,1100), lwd = 2)
text(x = mod.grafen.test$CR, y = 1200, labels = round(mod.grafen.test$CR, digits = 2), cex=2)
dev.off()

######################################################################################
## Figure 4: Pairwise Procrustes distances between the shape of male and female genitalia and
##       somatic characters in function of the phylogenetic distance.

## Phylogenetic distance:
dist.tr <- cophenetic(tr)

## Pairwise procrustes distances between shapes:
pp.male <- make.proc.dist(cord.gen.male)
pp.female <- make.proc.dist(cord.gen.female)
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

pdf("Figure_pairwise_procrustes_distance.pdf", height = 14, width = 7)
par(mfrow = c(2,1))
plot(c(dist.tr), c(pp.male), ylim = c(0,0.4), xlim = c(0, 12), xlab = ""
     , ylab = "", axes = F, col = "orangered", pch=20)
abline(lm(c(pp.male) ~ c(dist.tr)), lwd = 2, col = "black")
points(c(dist.tr), c(pp.scu.male), pch = 4, col = "slategray2")
abline(lm(c(pp.scu.male) ~ c(dist.tr)), lty = 1, lwd = 2, col = "slategray2")
points(c(dist.tr), c(pp.pro.male), pch = 4, col = "slategray2")
abline(lm(c(pp.pro.male) ~ c(dist.tr)), lty = 1, lwd = 2, col = "slategray2")
points(c(dist.tr), c(pp.jug.male), pch = 4, col = "slategray2")
abline(lm(c(pp.jug.male) ~ c(dist.tr)), lty = 1, lwd = 2, col = "slategray2")
axis(side = 2, at = c(0.0, 0.2, 0.4), cex.axis = 0.8)
## points(x = 12, y = 0.35, pch = -0x2642L, cex = 3.0)
mtext("Procrustes distance", side = 2, line = 2.5)

plot(c(dist.tr), c(pp.female), ylim = c(0,0.4), xlim = c(0,12),xlab = ""
     , ylab = "", axes = F, col = "orangered", pch = 20)
abline(lm(c(pp.female) ~ c(dist.tr)), lwd = 2, col = "black")
points(c(dist.tr), c(pp.scu.female), pch = 4, col = "slategray2")
abline(lm(c(pp.scu.female) ~ c(dist.tr)), lty = 1, lwd = 2, col = "slategray2")
points(c(dist.tr), c(pp.pro.female), pch = 4, col = "slategray2")
abline(lm(c(pp.pro.female) ~ c(dist.tr)), lty = 1, lwd = 2, col = "slategray2")
points(c(dist.tr), c(pp.jug.female), pch = 4, col = "slategray2")
abline(lm(c(pp.jug.female) ~ c(dist.tr)), lty = 1, lwd = 2, col = "slategray2")
axis(side = 2, at = c(0.0, 0.2, 0.4), cex.axis = 0.8)
axis(side = 1, cex.axis = 0.8)
## points(x = 12, y = 0.35, pch = -0x2640L, cex = 3.0)
mtext("Procrustes distance", side = 2, line = 2.5)
mtext("Phylogenetic distance", side = 1, line = 2.5)
dev.off()

######################################################################################
## Figure 5: Evolutionary rates of the male and female genitalia estimated using the ultrametric
##       phylogeny.

pdf("Figure_5_evolutionary_rates.pdf")
par(mar = c(4.0,4.0,2.0,2.0))
dd.null.uncor <- density(comp.gen.graphen$null.uncor, from = 0, to = 4)
dd.null.cor <- density(comp.gen.graphen$null.cor, from = 0, to = 4)
qq.uncor <- quantile(comp.gen.graphen$null.uncor, probs = 0.95)
qq.cor <- quantile(comp.gen.graphen$null.cor, probs = 0.95)
ii.uncor <- min(which(dd.null.uncor$x >= qq.uncor))
ii.cor <- min(which(dd.null.cor$x >= qq.cor))
plot(NA, main="", xlab = "", ylab = "", xlim = c(0,4), ylim = c(0,5), axes=FALSE, lty = 3, lwd = 3)
polygon(x = dd.null.uncor$x[c(ii.uncor, ii.uncor:length(dd.null.uncor$x), length(dd.null.uncor$x) )]
        , y = c(0, dd.null.uncor$y[ii.uncor:length(dd.null.uncor$x)], 0 )
        , col="grey", border = "white")
polygon(x = dd.null.cor$x[c(ii.cor, ii.cor:length(dd.null.cor$x), length(dd.null.cor$x) )]
        , y = c(0, dd.null.cor$y[ii.cor:length(dd.null.cor$x)], 0 )
        , col="grey", border = "white")
lines(dd.null.uncor, lty = 3, lwd = 3)
lines(dd.null.cor, lwd = 3)
obs.rate <- comp.gen.graphen$obs[1] / comp.gen.graphen$obs[2]
segments(x0 = obs.rate, y0 = 0, x1 = obs.rate, y1 = 1.5, lwd = 2, lty = 1)
text(x = obs.rate, y = 1.7, labels = round(obs.rate, 1))
legend(x = 2.5, y = 4, legend = c("Uncorrelated","Correlated"), bty = "n", lty = c(3,1), lwd = 3)
axis(side = 2, at = 0:5); axis(side = 1)
mtext(expression(paste(sigma["mult.M"]^2 / sigma["mult.F"]^2)), side = 1, line = 2.5)
mtext("Density", side = 2, line = 2)
dev.off()

######################################################################################
## Some additional plots not included in the manuscript:

## Plot the ultrametric tree with branch lengths by Graphen (1989) and the tree
##     where all branch lengths are equal to 1.
par(mfrow = c(1,2))
plot.phylo(tr.grafen, edge.width = 1.5, label.offset = 0.03, font = 4, direction = "upwards")
axisPhylo(side = 2)
plot.phylo(tr, edge.width = 1.5, label.offset = 0.03, font = 4, direction = "upwards")
axisPhylo(side = 2)

## Plot the phylomorphospace.
## With the speciational tree:
par(mfrow = c(1,2))
plotGMPhyloMorphoSpace(tr, cord.gen.male)
title("Branch lengths equal 1: Males")
plotGMPhyloMorphoSpace(tr, cord.gen.female)
title("Branch lengths equal 1: Females")
## With the ultrametric tree:
par(mfrow = c(1,2))
plotGMPhyloMorphoSpace(tr.grafen, cord.ind.male)
title("Ultrametric tree: Males")
plotGMPhyloMorphoSpace(tr.grafen, cord.ind.female)
title("Ultrametric tree: Females")
