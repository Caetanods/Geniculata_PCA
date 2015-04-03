## Analysis of integration between the shape of the male and the female genitalia.
## Calculate a PLS analysis of the contrasts data.

## Load packages and functions
library(geomorph)
library(geiger)
library(nsprcomp)
library(parallel)
library(Matrix)
library(grDevices)
source("./functions/analysis.R")
source("./functions/prepare-data.R")

## Get data. See 'Prepare_data.R':
load("./data/Geniculata_data.RData")

## Get phylogenetic contrasts:
mm.ind <- to.matrix(cord.ind.male)
mf.ind <- to.matrix(cord.ind.female)
mm.scu <- to.matrix(cord.scu.male)
mf.scu <- to.matrix(cord.scu.female)
mm.pro <- to.matrix(cord.pro.male)
mf.pro <- to.matrix(cord.pro.female)
mm.jug <- to.matrix(cord.jug.male)
mf.jug <- to.matrix(cord.jug.female)

## Use untrametric tree:
pic.grafen.ind.m <- sapply(1:dim(mm.ind)[2], function(x) pic(mm.ind[,x], tr.grafen) )
pic.grafen.ind.f <- sapply(1:dim(mf.ind)[2], function(x) pic(mf.ind[,x], tr.grafen) )
pic.grafen.scu.m <- sapply(1:dim(mm.scu)[2], function(x) pic(mm.scu[,x], tr.grafen) )
pic.grafen.scu.f <- sapply(1:dim(mf.scu)[2], function(x) pic(mf.scu[,x], tr.grafen) )
pic.grafen.pro.m <- sapply(1:dim(mm.pro)[2], function(x) pic(mm.pro[,x], tr.grafen) )
pic.grafen.pro.f <- sapply(1:dim(mf.pro)[2], function(x) pic(mf.pro[,x], tr.grafen) )
pic.grafen.jug.m <- sapply(1:dim(mm.jug)[2], function(x) pic(mm.jug[,x], tr.grafen) )
pic.grafen.jug.f <- sapply(1:dim(mf.jug)[2], function(x) pic(mf.jug[,x], tr.grafen) )

## Use unit tree:
pic.unit.ind.m <- sapply(1:dim(mm.ind)[2], function(x) pic(mm.ind[,x], tr) )
pic.unit.ind.f <- sapply(1:dim(mf.ind)[2], function(x) pic(mf.ind[,x], tr) )
pic.unit.scu.m <- sapply(1:dim(mm.scu)[2], function(x) pic(mm.scu[,x], tr) )
pic.unit.scu.f <- sapply(1:dim(mf.scu)[2], function(x) pic(mf.scu[,x], tr) )
pic.unit.pro.m <- sapply(1:dim(mm.pro)[2], function(x) pic(mm.pro[,x], tr) )
pic.unit.pro.f <- sapply(1:dim(mf.pro)[2], function(x) pic(mf.pro[,x], tr) )
pic.unit.jug.m <- sapply(1:dim(mm.jug)[2], function(x) pic(mm.jug[,x], tr) )
pic.unit.jug.f <- sapply(1:dim(mf.jug)[2], function(x) pic(mf.jug[,x], tr) )

## Two Block PLS analysis for ultrametric and speciational trees:

## Male genitalia vs. female genitalia:
pls.ind.grafen <- two.b.pls(pic.grafen.ind.m, pic.grafen.ind.f, verbose = TRUE)
pls.ind.unit <- two.b.pls(pic.unit.ind.m, pic.unit.ind.f, verbose = TRUE)

## Male genitalia vs. male scutellum and pronotum:
pls.scu.m.grafen <- two.b.pls(pic.grafen.ind.m, pic.grafen.scu.m, verbose = TRUE)
pls.scu.m.unit <- two.b.pls(pic.unit.ind.m, pic.unit.scu.m, verbose = TRUE)
pls.pro.m.grafen <- two.b.pls(pic.grafen.ind.m, pic.grafen.pro.m, verbose = TRUE)
pls.pro.m.unit <- two.b.pls(pic.unit.ind.m, pic.unit.pro.m, verbose = TRUE)
pls.jug.m.grafen <- two.b.pls(pic.grafen.ind.m, pic.grafen.jug.m, verbose = TRUE)
pls.jug.m.unit <- two.b.pls(pic.unit.ind.m, pic.unit.jug.m, verbose = TRUE)

## Female genitalia vs. female scutellum and pronotum:
pls.scu.f.grafen <- two.b.pls(pic.grafen.ind.f, pic.grafen.scu.f, verbose = TRUE)
pls.scu.f.unit <- two.b.pls(pic.unit.ind.f, pic.unit.scu.f, verbose = TRUE)
pls.pro.f.grafen <- two.b.pls(pic.grafen.ind.f, pic.grafen.pro.f, verbose = TRUE)
pls.pro.f.unit <- two.b.pls(pic.unit.ind.f, pic.unit.pro.f, verbose = TRUE)
pls.jug.f.grafen$pvalue <- two.b.pls(pic.grafen.ind.f, pic.grafen.jug.f, verbose = TRUE)
pls.jug.f.unit <- two.b.pls(pic.unit.ind.f, pic.unit.jug.f, verbose = TRUE)

## Making result plots:

## Plotting the first dimension of PLS results.
pdf("PLS_plots.pdf", width = 3, height = 9)
par(mfrow = c(3,1))
cex <- 2

## Male and female genitalia.
par(mar=c(2.5,4,2.5,2))
two.b.pls.plot(pls.ind.grafen, "PLS1 of male genitalia contrasts", "PLS1 of female genitalia contrasts", cex)

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

dev.off()
