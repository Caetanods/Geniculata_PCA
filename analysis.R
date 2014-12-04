## Prepare morphometric data for analysis.
## Landmarks and procrustes analysis exported from MorphoJ.

## Load packages and functions
library(geomorph)
library(geiger)
source("./functions/prepare-data.R")

## Get data
tr <- read.tree("./data/tnt_geniculata.tre")
raw.male <- read.table("./data/Raw.coord_male_20lmk.txt", header = T, sep = "\t")[,-1]
raw.female <- read.table("./data/Raw.coord_female.txt", header = T, sep = "\t")[,-1]
raw.male <- to.geomorph(raw.male)
raw.female <- to.geomorph(raw.female)

## Define which landmarks will slide.
slide <- define.sliders.2d(raw.male[,,1], nsliders = 18)

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

## Some plotting:

## pdf("Transformed_tree.pdf")
## plot.phylo(tr.grafen, edge.width = 1.5, label.offset = 0.03, font = 4)
## edgelabels(round(tr.grafen$edge.length, 2), frame = "none", cex = 0.6, adj = c(0.5,-1))
## axisPhylo(side = 1)
## title(main = "Branch length computed by Grafen", cex = 0.8)
## dev.off()

pdf("PC_phylomorphospace.pdf", width = 16, height = 8)
par(mfrow = c(1,2))
plotGMPhyloMorphoSpace(tr, cord.male)
title("Males")
plotGMPhyloMorphoSpace(tr, cord.female)
title("Females")
dev.off()

## Now using an ultrametric tree:
pdf("PC_ultrametric_phylomorphospace.pdf", width = 16, height = 8)
par(mfrow = c(1,2))
plotGMPhyloMorphoSpace(tr.grafen, cord.male)
title("Males")
plotGMPhyloMorphoSpace(tr.grafen, cord.female)
title("Females")
dev.off()

##############################
## Start analysis.

## Using the morphological integration analysis (Adams and Felice, 2014).
## Phylo with branch lengths equal to 1.
pdf("phylo.pls_unit.pdf")
int.unit <- phylo.pls(cord.male, cord.female, tr, iter=1000, verbose = TRUE)
dev.off()
## Phylo with ultrametric branch lengths.
pdf("phylo.pls_grafen.pdf")
int.grafen <- phylo.pls(cord.male, cord.female, tr.grafen, iter=1000, verbose = TRUE)
dev.off()

## Check pvalue.
int.unit$pvalue
int.grafen$pvalue
## Check PLS correlation.
int.unit[[1]]
int.grafen[[1]]

## Check if data differ from a single rate BM model:
male.sig.eq <- physignal(tr, cord.male, method = "Kmult", iter = 999)
male.sig.gr <- physignal(tr.grafen, cord.male, method = "Kmult", iter = 999)
female.sig.eq <- physignal(tr, cord.female, method = "Kmult", iter = 999)
female.sig.gr <- physignal(tr.grafen, cord.female, method = "Kmult", iter = 999)
print("Males:")
matrix(c(male.sig.eq$phy.signal, male.sig.gr$phy.signal, male.sig.eq$pvalue, male.sig.gr$pvalue)
     , byrow = TRUE, ncol = 2, dimnames = list(c("phy.signal","p.value"),c("all.equal","grafen")) )
print("Females:")
matrix(c(female.sig.eq$phy.signal, female.sig.gr$phy.signal, female.sig.eq$pvalue, female.sig.gr$pvalue)
     , byrow = TRUE, ncol = 2, dimnames = list(c("phy.signal","p.value"),c("all.equal","grafen")) )
## We can see that males and females do not show significant phylogenetic signal under BM.
## Values of K_mult ~0.5 indicates that the shape have more variance than expected under a BM model.

## Calculating rates of evolution under a BM model (Adams, 2014, Sys Bio).

## Using the tree with branch lengths equal to 1:
pdf(file = "Monte_Carlo_unit_plot.pdf")
res <- geo.comp.rates(tr, cord.male, cord.female, iter = 1999)
dev.off()

## Results
res$sigmaA / res$sigmaB
res$p.value

## Using the tree with branch lengths by Grafen:
pdf(file = "Monte_Carlo_grafen_plot.pdf")
res <- geo.comp.rates(tr.grafen, cord.male, cord.female, iter = 1999)
dev.off()

## Results
res$sigmaA / res$sigmaB
res$p.value

## Does the genitalia evolve faster under BM than somatic data?
## Calculate relative rates of genitalia:

## Get scutelum shape data:

## For males:


## Now I can use the metric of Sidlauskas (2008). Sidlauskas, B. 2008. Continuous and Arrested Morphological Diversification in Sister Clades of Characiform Fishes: A Phylomorphospace Approach. Evolution 62:3135–3156.
## The idea is to calculate the density index for males and for females. The calculus that Brian used in this paper was developed for two clades in the same tree, however, in our case we have more than one trait for each tip and a single group including all the tips.
## We might need to modify his calculations to do something more direct in this case.

## Remember also to analyse the data from the somatic measurements. The traits that were not derived from the genitalia will be used to understand the relative rate of evolution of the reproductive traits. We can thus compare female and male rates of evolution with the basal rate derived from the somatic traits.
