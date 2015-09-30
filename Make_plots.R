## This script reproduces all plots reported in the manuscript.
## Please refer to "Make_analysis.R" to perform the analysis and produce the objects that the
##       plots are based on. As an alternative you can just load the data with results of previous
##       simulations. We encourage the reader to run the analysis again for reproducibility.

## Uncomment to load previous simulations and results:
## load("./data/all_analyses_results.RData")

## Figure 1: Not made in R.

######################################################################################
## Figure 2: Correlated evolutionary changes of shape.

par(mfrow = c(3,1))
cex <- 2

## Male and female genitalia.
par(mar=c(2.5,4,2.5,2))
two.b.pls.plot(pls.gen.grafen, "PLS1 of male genitalia contrasts", "PLS1 of female genitalia contrasts", cex)

## Male and somatic.
up1 <- round(max(pls.scu.m.grafen$x.scores), digits = 1) + 0.1
up2 <- round(max(pls.pro.m.grafen$x.scores), digits = 1) + 0.1
up <- max(up1,up2)
par(mar=c(2.5,4,2.5,2))
plot(pls.scu.m.grafen$x.scores, pls.scu.m.grafen$y.scores, pch = 21, col = "grey"
     , xlab = "", ylab = "", axes = FALSE, xlim = c(-up,up), ylim = c(-up,up), cex = cex)
points(pls.pro.m.grafen$x.scores, pls.pro.m.grafen$y.scores, cex = cex, pch = 24, col = "grey")
points(pls.jug.m.grafen$x.scores, pls.jug.m.grafen$y.scores, cex = cex, pch = 23, col = "grey")
abline(lm(pls.scu.m.grafen$y.scores ~ pls.scu.m.grafen$x.scores), lwd = cex, col = "black")
abline(lm(pls.pro.m.grafen$y.scores ~ pls.pro.m.grafen$x.scores), lwd = cex, lty = 1, col = "black")
abline(lm(pls.jug.m.grafen$y.scores ~ pls.jug.m.grafen$x.scores), lwd = cex, lty = 1, col = "black")
axis(side = 1); axis(side = 2)
mtext(side = 1, text = "PLS1 of male somatic contrasts", line = 2.5)
mtext(side = 2, text = "PLS1 of male genitalia contrasts", line = 2.5)

## Female and somatic.
up1 <- round(max(pls.scu.f.grafen$x.scores), digits = 1) + 0.1
up2 <- round(max(pls.pro.f.grafen$x.scores), digits = 1) + 0.1
up <- max(up1,up2)
par(mar=c(2.5,4,2.5,2))
plot(pls.scu.f.grafen$x.scores, pls.scu.f.grafen$y.scores, pch = 21, col = "grey"
     , xlab = "", ylab = "", axes = FALSE, xlim = c(-up,up), ylim = c(-up,up), cex = cex)
points(pls.pro.f.grafen$x.scores, pls.pro.f.grafen$y.scores, cex = cex, pch = 24, col = "grey")
points(pls.jug.f.grafen$x.scores, pls.jug.f.grafen$y.scores, cex = cex, pch = 23, col = "grey")
abline(lm(pls.scu.f.grafen$y.scores ~ pls.scu.f.grafen$x.scores), lwd = cex, col = "black")
abline(lm(pls.pro.f.grafen$y.scores ~ pls.pro.f.grafen$x.scores), lwd = cex, lty = 1, col = "black")
abline(lm(pls.jug.f.grafen$y.scores ~ pls.jug.f.grafen$x.scores), lwd = cex, lty = 1, col = "black")
axis(side = 1); axis(side = 2)
mtext(side = 1, text = "PLS1 of female somatic contrasts", line = 2.5)
mtext(side = 2, text = "PLS1 of female genitalia contrasts", line = 2.5)

######################################################################################
## Figure 3: Not made in R.

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

######################################################################################
## Figure 5: Evolutionary rates of the male and female genitalia estimated using the ultrametric
##       phylogeny.
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