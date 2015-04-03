## This script prepare the data to use in all subsequent analysis.

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

## Get scutelum, pronotum and juga shape data:
raw.scu <- read.table("./data/Raw_scutelum_male_female.txt", header=T, sep="\t", stringsAsFactors=FALSE)[,-1]
raw.scu.male <- raw.scu[which(raw.scu$sexo == "m"),-2]
raw.scu.female <- raw.scu[which(raw.scu$sexo == "f"),-2]
raw.scu.male <- to.geomorph(raw.scu.male)
raw.scu.female <- to.geomorph(raw.scu.female)

raw.pro <- read.table("./data/Pronotum_male_female.txt", header=T, sep="\t", stringsAsFactors=FALSE)[,-1]
raw.pro.male <- raw.pro[which(raw.pro$sex == "m"),-2]
raw.pro.female <- raw.pro[which(raw.pro$sex == "f"),-2]
pro.male <- to.geomorph(raw.pro.male)
pro.female <- to.geomorph(raw.pro.female)

raw.jug <- read.table("./data/juga_male_female_raw.txt", header=T, sep="\t", stringsAsFactors=FALSE)[,-1]
raw.jug <- raw.jug[-which(raw.jug$Sex == "e"),]
raw.jug.male <- raw.jug[which(raw.jug$Sex == "m"),-2]
raw.jug.female <- raw.jug[which(raw.jug$Sex == "f"),-2]
jug.male <- to.geomorph(raw.jug.male)
jug.female <- to.geomorph(raw.jug.female)

## Decrease data to 20 landmarks
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
pro.male <- gpagen(pro.male, ProcD = FALSE, ShowPlot = FALSE, curves = slide)
pro.female <- gpagen(pro.female, ProcD = FALSE, ShowPlot = FALSE, curves = slide)
jug.male <- gpagen(jug.male, ProcD = FALSE, ShowPlot = FALSE, curves = slide)
jug.female <- gpagen(jug.female, ProcD = FALSE, ShowPlot = FALSE, curves = slide)

## Species mean value for each coordinate of the procrustes centralized data.
cord.ind.male <- to.mean.shape(ind.male)
cord.ind.female <- to.mean.shape(ind.female)
dimnames(cord.ind.male)[[3]] <- dimnames(cord.ind.female)[[3]]
cord.scu.male <- to.mean.shape(scu.male)
cord.scu.female <- to.mean.shape(scu.female)
cord.pro.male <- to.mean.shape(pro.male)
cord.pro.female <- to.mean.shape(pro.female)
cord.jug.male <- to.mean.shape(jug.male)
cord.jug.female <- to.mean.shape(jug.female)

## Set phylogenetic labels equal to data:
tr$tip.label <- gsub("_"," ", tr$tip.label)

## Get ultrametric trees by branch trasformations:
tr.grafen <- compute.brlen(tr)

## Save the important objects:
save(cord.ind.male, cord.ind.female, cord.scu.male, cord.scu.female, cord.pro.male
         , cord.pro.female, cord.jug.male, cord.jug.female, tr, tr.grafen
         , file = "./data/Geniculata_data.RData")
