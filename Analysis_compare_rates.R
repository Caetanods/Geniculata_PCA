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


## Get data
tr <- read.tree("./data/tnt_geniculata.tre")
raw.male <- read.table("./data/Raw.coord_male_20lmk.txt", header = T, sep = "\t")[,-1]
raw.female <- read.table("./data/Raw.coord_female.txt", header = T, sep = "\t")[,-1]
raw.male <- to.geomorph(raw.male)
raw.female <- to.geomorph(raw.female)

##############################################################
## Define which landmarks will slide.
##slide <- define.sliders.2d(raw.male[,,1], nsliders = 18)
slide <- as.matrix(read.csv("./data/curveslide.csv"))

## Procrustes superposition step with sliding semilandmarks minimizing bending energy.
ind.male <- gpagen(raw.male, ProcD = FALSE, ShowPlot = FALSE, curves = slide)
ind.female <- gpagen(raw.female, ProcD = FALSE, ShowPlot = FALSE, curves = slide)

## Species mean value for each coordinate of the procrustes centralized data.
cord.male <- to.mean.shape(ind.male)
cord.female <- to.mean.shape(ind.female)
dimnames(cord.male)[[3]] <- dimnames(cord.female)[[3]]

## Set phylogenetic labels equal to data:
tr$tip.label <- gsub("_"," ", tr$tip.label)

## Get ultrametric trees by branch trasformations:
tr.grafen <- compute.brlen(tr)

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
## plotGMPhyloMorphoSpace(tr, cord.male)
## title("Branch lengths equal 1: Males")
## plotGMPhyloMorphoSpace(tr, cord.female)
## title("Branch lengths equal 1: Females")

## par(mfrow = c(1,2))
## plotGMPhyloMorphoSpace(tr.grafen, cord.male)
## title("Ultrametric tree: Males")
## plotGMPhyloMorphoSpace(tr.grafen, cord.female)
## title("Ultrametric tree: Females")

###################################################
## Start analysis.

## Check if the male genitalia is the most modular structure.
## This is expected if they evolve faster and show more species level variation.

## From geomorph format back to matrix:
dat.male <- to.matrix(cord.male)
dat.female <- to.matrix(cord.female)

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

## comp.graphen <- geo.comp.rates(tr.grafen, cord.male, cord.female, plot = FALSE)
## comp.unit <- geo.comp.rates(tr, cord.male, cord.female, plot = FALSE)
## save(comp.graphen, comp.unit, file = "./data/male_female_rates.RData")
load("./data/male_female_rates.RData")

## Make density plot for ultrametric tree:
## pdf("manuscript_figures/Density_Genitalia_male_female_graphen.pdf")
## par(mar = c(4.0,4.0,2.0,2.0))
## plot(density(comp.graphen$null.uncor), zero.line = FALSE, main=""
##    , xlab = "", ylab = ""
##    , xlim = c(0,4), ylim = c(0,5), axes=FALSE, lty = 3, lwd = 3)
## lines(density(comp.graphen$null.cor), lwd = 3)
## obs.rate <- comp.graphen$obs[1] / comp.graphen$obs[2]
## segments(x0 = obs.rate, y0 = 0, x1 = obs.rate, y1 = 1.5, lwd = 2, lty = 1)
## text(x = obs.rate, y = 1.7, labels = round(obs.rate, 1))
## legend(x = 2.5, y = 4, legend = c("Uncorrelated","Correlated"), bty = "n", lty = c(3,1), lwd = 3)
## axis(side = 2, at = 0:5); axis(side = 1)
## mtext(expression(paste(sigma["mult.M"]^2 / sigma["mult.F"]^2)), side = 1, line = 2.5)
## mtext("Density", side = 2, line = 2)
## dev.off()

## ## For the unit tree:
## pdf("manuscript_figures/Density_Genitalia_male_female_unit.pdf")
## par(mar = c(4.0,4.0,2.0,2.0))
## plot(density(comp.unit$null.uncor), zero.line = FALSE, main=""
##    , xlab = "", ylab = ""
##    , xlim = c(0,4), ylim = c(0,5), axes=FALSE, lty = 3, lwd = 3)
## lines(density(comp.unit$null.cor), lwd = 3)
## obs.rate <- comp.unit$obs[1] / comp.unit$obs[2]
## segments(x0 = obs.rate, y0 = 0, x1 = obs.rate, y1 = 1.5, lwd = 2, lty = 1)
## text(x = obs.rate, y = 1.7, labels = round(obs.rate, 1))
## legend(x = 2.5, y = 4, legend = c("Uncorrelated","Correlated"), bty = "n", lty = c(3,1), lwd = 3)
## axis(side = 2, at = 0:5); axis(side = 1)
## mtext(expression(paste(sigma["mult.M"]^2 / sigma["mult.F"]^2)), side = 1, line = 2.5)
## mtext("Density", side = 2, line = 2)
## dev.off()

## ## Plate version:
## pdf("manuscript_figures/Density_Genitalia_male_female_plate.pdf", width = 14)
## par(mfrow = c(1,2))
## par(mar = c(4.0,4.0,2.0,0.0))
## plot(density(comp.graphen$null.uncor), zero.line = FALSE, main=""
##    , xlab = "", ylab = ""
##    , xlim = c(0,4), ylim = c(0,5), axes=FALSE, lty = 3, lwd = 3)
## lines(density(comp.graphen$null.cor), lwd = 3)
## obs.rate <- comp.graphen$obs[1] / comp.graphen$obs[2]
## segments(x0 = obs.rate, y0 = 0, x1 = obs.rate, y1 = 1.5, lwd = 2, lty = 1)
## text(x = obs.rate, y = 1.7, labels = round(obs.rate, 1))
## legend(x = 2.5, y = 4, legend = c("Uncorrelated","Correlated"), bty = "n", lty = c(3,1), lwd = 3)
## axis(side = 2, at = 0:5); axis(side = 1)
## mtext(expression(paste(sigma["mult.M"]^2 / sigma["mult.F"]^2)), side = 1, line = 2.5)
## mtext("Density", side = 2, line = 2)
## par(mar = c(4.0,3.0,2.0,1.0))
## plot(density(comp.unit$null.uncor), zero.line = FALSE, main=""
##    , xlab = "", ylab = ""
##    , xlim = c(0,4), ylim = c(0,5), axes=FALSE, lty = 3, lwd = 3)
## lines(density(comp.unit$null.cor), lwd = 3)
## obs.rate <- comp.unit$obs[1] / comp.unit$obs[2]
## segments(x0 = obs.rate, y0 = 0, x1 = obs.rate, y1 = 1.5, lwd = 2, lty = 1)
## text(x = obs.rate, y = 1.7, labels = round(obs.rate, 1))
## legend(x = 2.5, y = 4, legend = c("Uncorrelated","Correlated"), bty = "n", lty = c(3,1), lwd = 3)
## axis(side = 2, at = 0:5); axis(side = 1)
## mtext(expression(paste(sigma["mult.M"]^2 / sigma["mult.F"]^2)), side = 1, line = 2.5)
## mtext("Density", side = 2, line = 2)
## dev.off()

## Using Monte Carlo simulation produces empirical p.values.
## Then we analysed multiple times to get a distribution of p.values:
## NOTE: We do not need to do this. We need to increase the number of simulations in the distribution
##       of the Monte Carlo for each analysis.
## p.graphen <- mclapply(1:20, FUN = function(x) geo.comp.rates(tr.grafen, cord.male, cord.female), mc.cores = 2)
## p.graphen <- lapply(p.graphen, FUN = function(x) x$p.value)
## p.graphen <- do.call(rbind, p.graphen)
## p.unit <- mclapply(1:20, FUN = function(x) geo.comp.rates(tr, cord.male, cord.female), mc.cores = 2)
## p.unit <- lapply(p.unit, FUN = function(x) x$p.value)
## p.unit <- do.call(rbind, p.unit)
## ## Summary of distributions of p.value using the ultrametric and the unit tree.
## summary(p.graphen)
## summary(p.unit)

################################################################
## Does the genitalia evolve faster under BM than somatic data?

## Analysis of rate of evolution between genitalia and scutelum.

## Get scutelum shape data:
raw.scu <- read.table("./data/Raw_scutelum_male_female.txt", header=T, sep="\t", stringsAsFactors=FALSE)[,-1]
raw.scu.male <- raw.scu[which(raw.scu$sexo == "m"),-2]
raw.scu.female <- raw.scu[which(raw.scu$sexo == "f"),-2]
raw.scu.male <- to.geomorph(raw.scu.male)
raw.scu.female <- to.geomorph(raw.scu.female)

## Procrustes superposition step with sliding semilandmarks minimizing bending energy.
## Again, we can use the same slide matrix from the previous datasets.
ind.scu.male <- gpagen(raw.scu.male, ProcD = FALSE, ShowPlot = FALSE, curves = slide)
ind.scu.female <- gpagen(raw.scu.female, ProcD = FALSE, ShowPlot = FALSE, curves = slide)

## Species mean value for each coordinate of the procrustes centralized data.
cord.scu.male <- to.mean.shape(ind.scu.male)
cord.scu.female <- to.mean.shape(ind.scu.female)

## Relative rates in function of the scutellum for males and females:
## Note that the rates of the genitalia are much faster.
## For the ultrametric tree:
## m.s.graphen <- geo.comp.rates(tr.grafen, cord.male, cord.scu.male, plot = TRUE)
## f.s.graphen <- geo.comp.rates(tr.grafen, cord.female, cord.scu.female, plot = TRUE)
## ## For the tree with branch lengths all equal to 1:
## m.s.unit <- geo.comp.rates(tr.grafen, cord.male, cord.scu.male, plot = TRUE)
## f.s.unit <- geo.comp.rates(tr.grafen, cord.female, cord.scu.female, plot = TRUE)
## save(m.s.graphen, f.s.graphen, m.s.unit, f.s.unit, file = "./data/genitalia_scutelum_rates.RData")
load("./data/genitalia_scutelum_rates.RData")

## Both male and female genitalia evolve faster independent of the branch lengths and correlation
##      of the data:
m.s.graphen$p.value ## Male in ultrametric tree.
f.s.graphen$p.value ## Female in ultrametric tree.
m.s.unit$p.value ## Male in speciational tree.
f.s.unit$p.value ## Female in speciational tree.

## Analysis of rate of evolution between genitalia and the metacoxa.
## Note that the metacoxa has linear measurements.
meta <- read.csv("./data/metacoxa_linear.csv")
agg <- aggregate(meta$right.metacoxae, by = list(meta$sex, meta$species), FUN = mean)
mm <- which(agg$Group.1 == "m")
ff <- which(agg$Group.1 == "f")
m.meta <- agg$x[mm] ## Mean spp values for males.
names(m.meta) <- agg$Group.2[mm]
f.meta <- agg$x[ff] ## Mean spp values for females.
names(f.meta) <- agg$Group.2[ff]

## IMPORTANT: The analysis below are biased. The 'geo.vec.comp.rates' function will return
##            a (near) significant difference between the rates even if the trait were simulated
##            under the null model. Maybe the is a result of testing datasets with very different
##            dimentionality.

## Testing rates differences with raw values:
m.f.graphen <- geo.vec.comp.rates(tr.grafen, m.meta, cord.male, plot = TRUE)
f.f.graphen <- geo.vec.comp.rates(tr.grafen, f.meta, cord.female, plot = TRUE)
m.f.unit <- geo.vec.comp.rates(tr, m.meta, cord.male, plot = TRUE)
f.f.unit <- geo.vec.comp.rates(tr, f.meta, cord.female, plot = TRUE)

## Testing rates differences with centred values:
scale.m.f.graphen <- geo.vec.comp.rates(tr.grafen, scale(m.meta)[,1], cord.male, plot = TRUE)
log.m.f.graphen <- geo.vec.comp.rates(tr.grafen, log(m.meta), cord.male, plot = TRUE)
