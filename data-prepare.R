## Prepare morphometric data for analysis.
## Landmarks and procrustes analysis exported from MorphoJ.

## Load packages
library(geomorph)
library(geiger)

## Get data and functions
source("./functions/prepare-data.R")
ind.male <- read.table("./data/Procrustes_male.txt", header = T, sep = "\t")[,-1]
ind.female <- read.table("./data/Procrustes_female.txt", header = T, sep = "\t")[,-1]
tr <- read.tree("./data/tnt_geniculata.tre")

## Mean value for each coordinate of the procrustes centralized data.
cord.male <- aggregate(ind.male[,-1], by = list(ind.male$Sps), FUN = mean)
cord.female <- aggregate(ind.female[,-1], by = list(ind.female$Sps), FUN = mean)

## Prepare data in correct format.
m.male <- to.geomorph(cord.male)
m.female <- to.geomorph(cord.female)

## Set phylogenetic labels equal to data:
tr$tip.label <- gsub("_"," ", tr$tip.label)

## Get ultrametric trees by branch trasformations:
tr.grafen <- compute.brlen(tr)

## Some plotting:

pdf("Transformed_tree.pdf")
plot.phylo(tr.grafen, edge.width = 1.5, label.offset = 0.03, font = 4)
edgelabels(round(tr.grafen$edge.length, 2), frame = "none", cex = 0.6, adj = c(0.5,-1))
axisPhylo(side = 1)
title(main = "Branch length computed by Grafen", cex = 0.8)
dev.off()

pdf("PC_phylomorphospace.pdf", width = 16, height = 8)
par(mfrow = c(1,2))
plotGMPhyloMorphoSpace(tr, m.male)
title("Males")
plotGMPhyloMorphoSpace(tr, m.female)
title("Females")
dev.off()

## Now using an ultrametric tree:
pdf("PC_ultrametric_phylomorphospace.pdf", width = 16, height = 8)
par(mfrow = c(1,2))
plotGMPhyloMorphoSpace(tr.grafen, m.male)
title("Males")
plotGMPhyloMorphoSpace(tr.grafen, m.female)
title("Females")
dev.off()

## Check if data have differ from a BM model:
male.sig.eq <- physignal(tr, m.male, method = "Kmult", iter = 999)
male.sig.gr <- physignal(tr.grafen, m.male, method = "Kmult", iter = 999)
female.sig.eq <- physignal(tr, m.female, method = "Kmult", iter = 999)
female.sig.gr <- physignal(tr.grafen, m.female, method = "Kmult", iter = 999)
print("Males:")
matrix(c(male.sig.eq$phy.signal, male.sig.gr$phy.signal, male.sig.eq$pvalue, male.sig.gr$pvalue)
     , byrow = TRUE, ncol = 2, dimnames = list(c("phy.signal","p.value"),c("all.equal","grafen")) )
print("Females:")
matrix(c(female.sig.eq$phy.signal, female.sig.gr$phy.signal, female.sig.eq$pvalue, female.sig.gr$pvalue)
     , byrow = TRUE, ncol = 2, dimnames = list(c("phy.signal","p.value"),c("all.equal","grafen")) )
## Females deviate more from BM than males.

## Calculating rates of evolution under a BM model (Adams, 2014, Sys Bio)
## Checking if male genitalia evolve faster than female genitalia under BM
gp <- factor(rep(1, times = length(tr$tip.label)))
names(gp) <- tr$tip.label

rate.male <- compare.evol.rates(tr, m.male, gp = gp, iter = 5)
rate.female <- compare.evol.rates(tr, m.female, gp = gp, iter = 5)
rate.male$sigma.d / rate.female$sigma.d

## The value are different but is this different significative?
## Simulate multivariate data in the phylogeny from a normal distribution.
## Maybe here two normal distributions, each equal to one of the extimated data.
## Did a modification of the Dean Adams function to do this simulation:
source("./functions/prepare-data.R")
## The function is working. However, the simulation is not working.
## Guess that the simulation is made to work with the ratios and not
## with the absolute value. Need to tweek a little bit to work then.
res <- sim.geomorpho(tr, m.male, iter = 99)

## Fitting different models to the data:
## Fitting an OU model by trasnforming the tree:
## Are the fit better under other models?

## Not sure if I can fit other models of evolution this way.
## Need to review the characteristics of the K statistics.
to.grafen.OU <- rescale(tr.grafen, model = "OU")

OU.it.k <- list()
alpha <- rexp(10, 2)
for(i in 1:10){
    male.grafen.OU <- to.grafen.OU(alpha[i], rate.male$sigma.d)
    OU.it.k[[i]] <- physignal(male.grafen.OU, m.male, method = "Kmult", iter = 999)
}
OU.it.k[[10]]

## Using the morphological integration analysis (Adams and Felice, 2014)
phylo.pls(m.male, m.female, tr, iter=5)
