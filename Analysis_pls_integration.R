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


