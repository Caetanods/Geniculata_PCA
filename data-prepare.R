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

## Using the morphological integration analysis (Adams and Felice, 2014)

dimnames(m.male)[[3]]
phylo.pls(m.male[,,], m.female[,,], tr, iter=5)
