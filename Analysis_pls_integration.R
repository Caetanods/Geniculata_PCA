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

## Get phylogeny
tr <- read.tree("./data/tnt_geniculata.tre")

## Get genitalia data
raw.male <- read.table("./data/Raw.coord_male_20lmk.txt", header = T, sep = "\t")[,-1]
raw.female <- read.table("./data/Raw.coord_female.txt", header = T, sep = "\t")[,-1]
raw.male <- to.geomorph(raw.male)
raw.female <- to.geomorph(raw.female)

## Get scutelum data
raw.scu <- read.table("./data/Raw_scutelum_male_female.txt", header=T, sep="\t", stringsAsFactors=FALSE)[,-1]
raw.scu.male <- raw.scu[which(raw.scu$sexo == "m"),-2]
raw.scu.female <- raw.scu[which(raw.scu$sexo == "f"),-2]
raw.scu.male <- to.geomorph(raw.scu.male)
raw.scu.female <- to.geomorph(raw.scu.female)

## Decrease data for 20 landmarks
scu.male <- array(dim = c(20, 2, 63) , dimnames = dimnames(raw.scu.male))
scu.female <- array(dim = c(20, 2, 97), dimnames = dimnames(raw.scu.female) )

for(i in 1:dim(raw.scu.male)[3]){
    scu.male[,,i] <- raw.scu.male[-seq(2, 40, by = 2),,i]
}
for(i in 1:dim(raw.scu.female)[3]){
    scu.female[,,i] <- raw.scu.female[-seq(2, 40, by = 2),,i]
}


##############################################################
## Define which landmarks will slide.
## slide <- define.sliders.2d(raw.male[,,1], nsliders = 18)
slide <- as.matrix(read.csv("./data/curveslide_20.csv"))

## Procrustes superposition step with sliding semilandmarks minimizing bending energy.
ind.male <- gpagen(raw.male, ProcD = FALSE, ShowPlot = FALSE, curves = slide)
ind.female <- gpagen(raw.female, ProcD = FALSE, ShowPlot = FALSE, curves = slide)
scu.male <- gpagen(scu.male, ProcD = FALSE, ShowPlot = FALSE, curves = slide)
scu.female <- gpagen(scu.female, ProcD = FALSE, ShowPlot = FALSE, curves = slide)

## Species mean value for each coordinate of the procrustes centralized data.
cord.ind.male <- to.mean.shape(ind.male)
cord.ind.female <- to.mean.shape(ind.female)
cord.scu.male <- to.mean.shape(scu.male)
cord.scu.female <- to.mean.shape(scu.female)
dimnames(cord.ind.male)[[3]] <- dimnames(cord.ind.female)[[3]]

## Set phylogenetic labels equal to data:
tr$tip.label <- gsub("_"," ", tr$tip.label)

## Get ultrametric trees by branch trasformations:
tr.grafen <- compute.brlen(tr)

## Get phylogenetic contrasts:
mm.ind <- to.matrix(cord.ind.male)
mf.ind <- to.matrix(cord.ind.female)
mm.scu <- to.matrix(cord.scu.male)
mf.scu <- to.matrix(cord.scu.female)

## Use untrametric tree:
pic.grafen.ind.m <- sapply(1:dim(mm.ind)[2], function(x) pic(mm.ind[,x], tr.grafen) )
pic.grafen.ind.f <- sapply(1:dim(mf.ind)[2], function(x) pic(mf.ind[,x], tr.grafen) )
pic.grafen.scu.m <- sapply(1:dim(mm.scu)[2], function(x) pic(mm.scu[,x], tr.grafen) )
pic.grafen.scu.f <- sapply(1:dim(mf.scu)[2], function(x) pic(mf.scu[,x], tr.grafen) )

## Use unit tree:
pic.unit.ind.m <- sapply(1:dim(mm.ind)[2], function(x) pic(mm.ind[,x], tr) )
pic.unit.ind.f <- sapply(1:dim(mf.ind)[2], function(x) pic(mf.ind[,x], tr) )
pic.unit.scu.m <- sapply(1:dim(mm.scu)[2], function(x) pic(mm.scu[,x], tr) )
pic.unit.scu.f <- sapply(1:dim(mf.scu)[2], function(x) pic(mf.scu[,x], tr) )

## Two Block PLS analysis for ultrametric and speciational trees:

## Male genitalia vs. female genitalia:
pls.ind.grafen <- two.b.pls(pic.grafen.ind.m, pic.grafen.ind.f, verbose = TRUE)
pls.ind.unit <- two.b.pls(pic.unit.ind.m, pic.unit.ind.f, verbose = TRUE)

## Male genitalia vs. male scutellum:
pls.scu.m.grafen <- two.b.pls(pic.grafen.ind.m, pic.grafen.scu.m, verbose = TRUE)
pls.scu.m.unit <- two.b.pls(pic.unit.ind.m, pic.unit.scu.m, verbose = TRUE)

## Female genitalia vs. female scutellum:
pls.scu.f.grafen <- two.b.pls(pic.grafen.ind.f, pic.grafen.scu.f, verbose = TRUE)
pls.scu.f.unit <- two.b.pls(pic.unit.ind.f, pic.unit.scu.f, verbose = TRUE)

## Making result plots:

## Plotting the first dimension of PLS results.
par(mfrow = c(3,1))
cex <- 2.5
two.b.pls.plot(pls.ind.grafen, "PLS1 of male genitalia contrasts", "PLS1 of female genitalia contrasts", cex)
two.b.pls.plot(pls.scu.m.grafen, "PLS1 of male scutellum contrasts", "PLS1 of male genitalia contrasts", cex)
two.b.pls.plot(pls.scu.f.grafen, "PLS1 of female scutellum contrasts", "PLS1 of female genitalia contrasts", cex)
dev.copy2pdf()
