## Prepare morphometric data for analysis.
## Landmarks and procrustes analysis exported from MorphoJ.

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

#################################################
## Some plotting:

## Plot the ultrametric tree with branch lengths by Graphen (1989) and the tree
##   where all branch lengths are equal to 1.
## pdf("Phylogeny_branch_lengths.pdf", width = 14)
## par(mfrow = c(1,2))
## plot.phylo(tr.grafen, edge.width = 1.5, label.offset = 0.03, font = 4, direction = "upwards")
## axisPhylo(side = 2)
## plot.phylo(tr, edge.width = 1.5, label.offset = 0.03, font = 4, direction = "upwards")
## axisPhylo(side = 2)
## dev.off()

## Plot the phylomorphospace.
## par(mfrow = c(1,2))
## plotGMPhyloMorphoSpace(tr, cord.ind.male)
## title("Branch lengths equal 1: Males")
## plotGMPhyloMorphoSpace(tr, cord.ind.female)
## title("Branch lengths equal 1: Females")

## par(mfrow = c(1,2))
## plotGMPhyloMorphoSpace(tr.grafen, cord.ind.male)
## title("Ultrametric tree: Males")
## plotGMPhyloMorphoSpace(tr.grafen, cord.ind.female)
## title("Ultrametric tree: Females")

###################################################
## Start analysis.

## Check if the male genitalia is the most modular structure.
## This is expected if they evolve faster and show more species level variation.

## From geomorph format back to matrix:
dat.male <- to.matrix(cord.ind.male)
dat.female <- to.matrix(cord.ind.female)

## Calculate Sparse Principal Component Analysis (SPCA):
## This analysis is NOT analogous to a phylogenetic PCA.
male.nspr <- nsprcomp(dat.male, scale. = TRUE)
female.nspr <- nsprcomp(dat.female, scale. = TRUE)
## Calculate percent of variation per SPC component:
male.var.prop <- round(male.nspr$sdev^2 / sum(male.nspr$sdev^2), digits = 3)
female.var.prop <- round(female.nspr$sdev^2 / sum(female.nspr$sdev^2), digits = 3)

## Plot the percent of variation explained by each axis for male and female
##   genitalia. Note that female genitalia are more autocorrelated.
## pdf("./manuscript_figures/SPCA_genitalia.pdf")
## par(mar = c(4.0,4.0,2.0,2.0))
## plot(1:length(male.var.prop), male.var.prop, type = "b", pch = 4,
##      , ylim = c(0,1.0), xlim = c(1,length(male.var.prop))
##      , xlab = "", ylab = ""
##      , axes = FALSE )
## lines(1:length(female.var.prop), female.var.prop, type = "b")
## axis(side = 2); axis(side = 1, at = seq(1, length(male.var.prop)), pos = -0.03)
## mtext("SPC axis", side = 1, line = 2)
## mtext("% variation explained", side = 2, line = 2.2)
## legend(x=9, y=0.8, legend = c("Male","Female"), bty = "n", pch = c(4,1))
## dev.off()

################################################################
## Analysis of rate of evolution between genitalia of male and female.

## Using the geo.comp.rates function:
## This function make a similar analysis to compare.evol.rates in package 'geomorph'.
## Here we compare two multivariate datasets under the same phylogeny instead of the same
##   dataset with two or more clades comparison.
## In addition, we implement a more realistic null model in which the datasets to be
##   compared are truly multivariate, since Adams (2014) generate a null model in which the
##   simulations draw from an uncorrelated multivariate-normal distribution.
## Note how the addition of the covariances to the null model changed the distribution:

## comp.graphen <- geo.comp.rates(tr.grafen, cord.ind.male, cord.ind.female, plot = FALSE)
## comp.unit <- geo.comp.rates(tr, cord.ind.male, cord.ind.female, plot = FALSE)
## save(comp.graphen, comp.unit, file = "./data/male_female_rates.RData")
load("./data/male_female_rates.RData")

## plot(dd.null.uncor, zero.line = FALSE, main="", xlab = "", ylab = ""
##    , xlim = c(0,4), ylim = c(0,5), axes=FALSE, lty = 3, lwd = 3)

## Make density plot for ultrametric tree:
pdf("manuscript_figures/Density_Genitalia_male_female_graphen.pdf")
par(mar = c(4.0,4.0,2.0,2.0))
dd.null.uncor <- density(comp.graphen$null.uncor, from = 0, to = 4)
dd.null.cor <- density(comp.graphen$null.cor, from = 0, to = 4)
qq.uncor <- quantile(comp.graphen$null.uncor, probs = 0.95)
qq.cor <- quantile(comp.graphen$null.cor, probs = 0.95)
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
obs.rate <- comp.graphen$obs[1] / comp.graphen$obs[2]
segments(x0 = obs.rate, y0 = 0, x1 = obs.rate, y1 = 1.5, lwd = 2, lty = 1)
text(x = obs.rate, y = 1.7, labels = round(obs.rate, 1))
legend(x = 2.5, y = 4, legend = c("Uncorrelated","Correlated"), bty = "n", lty = c(3,1), lwd = 3)
axis(side = 2, at = 0:5); axis(side = 1)
mtext(expression(paste(sigma["mult.M"]^2 / sigma["mult.F"]^2)), side = 1, line = 2.5)
mtext("Density", side = 2, line = 2)
dev.off()

################################################################
## Does the genitalia evolve faster under BM than somatic data?

## Analysis of rate of evolution between genitalia and somatic traits.
## Relative rates in function of the scutellum for males and females:
## Note that the rates of the genitalia are much faster.
## For the ultrametric tree:

## m.s.graphen <- geo.comp.rates(tr.grafen, cord.ind.male, cord.scu.male, plot = TRUE)
## f.s.graphen <- geo.comp.rates(tr.grafen, cord.ind.female, cord.scu.female, plot = TRUE)
## ## For the tree with branch lengths all equal to 1:
## m.s.unit <- geo.comp.rates(tr, cord.ind.male, cord.scu.male, plot = TRUE)
## f.s.unit <- geo.comp.rates(tr, cord.ind.female, cord.scu.female, plot = TRUE)
## save(m.s.graphen, f.s.graphen, m.s.unit, f.s.unit, file = "./data/genitalia_scutelum_rates.RData")
load("./data/genitalia_scutelum_rates.RData")

## Both male and female genitalia evolve faster independent of the branch lengths and correlation
##      of the data:
m.s.graphen$p.value ## Male in ultrametric tree.
f.s.graphen$p.value ## Female in ultrametric tree.
m.s.unit$p.value ## Male in speciational tree.
f.s.unit$p.value ## Female in speciational tree.

## Comparing rates of evolution between the pronotum and the genitalia:

## For the ultrametric tree:
m.p.graphen <- geo.comp.rates(tr.grafen, cord.ind.male, cord.pro.male, plot = TRUE)
f.p.graphen <- geo.comp.rates(tr.grafen, cord.ind.female, cord.pro.female, plot = TRUE)
## For the tree with branch lengths all equal to 1:
m.p.unit <- geo.comp.rates(tr, cord.ind.male, cord.pro.male, plot = TRUE)
f.p.unit <- geo.comp.rates(tr, cord.ind.female, cord.pro.female, plot = TRUE)
save(m.p.graphen, f.p.graphen, m.p.unit, f.p.unit, file = "./data/genitalia_pronotum_rates.RData")
