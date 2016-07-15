## This script prepare the data to use in all subsequent analysis.

## Load packages and functions
library(geomorph) ## Version 2.1.7
library(geiger) ## Version 2.0.6
library(Matrix) ## Version 1.2
source("./functions/analysis.R")
source("./functions/prepare-data.R")

## Get phylogeny
tr <- read.tree("./data/tnt_geniculata.tre")

## Get genitalia data:
## Read from the original data files.
raw.gen.male <- read.table("./data/Raw.coord_male_20lmk.txt", header = T, sep = "\t")[,-1]
raw.gen.female <- read.table("./data/Raw.coord_female.txt", header = T, sep = "\t")[,-1]
raw.gen.male <- to.geomorph(raw.gen.male)
raw.gen.female <- to.geomorph(raw.gen.female)

## Get scutelum, pronotum and juga shape data:
## Read from the original data files.
raw.scu <- read.table("./data/Raw_scutelum_male_female.txt", header=T, sep="\t", stringsAsFactors=FALSE)[,-1]
raw.scu.male <- raw.scu[which(raw.scu$sexo == "m"),-2]
raw.scu.female <- raw.scu[which(raw.scu$sexo == "f"),-2]
raw.scu.male <- to.geomorph(raw.scu.male)
raw.scu.female <- to.geomorph(raw.scu.female)

raw.pro <- read.table("./data/Pronotum_male_female.txt", header=T, sep="\t", stringsAsFactors=FALSE)[,-1]
raw.pro.male <- raw.pro[which(raw.pro$sex == "m"),-2]
raw.pro.female <- raw.pro[which(raw.pro$sex == "f"),-2]
raw.pro.male <- to.geomorph(raw.pro.male)
raw.pro.female <- to.geomorph(raw.pro.female)

raw.jug <- read.table("./data/juga_male_female_raw.txt", header=T, sep="\t", stringsAsFactors=FALSE)[,-1]
raw.jug <- raw.jug[-which(raw.jug$Sex == "e"),]
raw.jug.male <- raw.jug[which(raw.jug$Sex == "m"),-2]
raw.jug.female <- raw.jug[which(raw.jug$Sex == "f"),-2]
raw.jug.male <- to.geomorph(raw.jug.male)
raw.jug.female <- to.geomorph(raw.jug.female)

## Decrease scutelum data to 20 landmarks.
## All somatic traits have equal number of landmarks.
reduce.male <- array(dim = c(20, 2, 63) , dimnames = dimnames(raw.scu.male))
reduce.female <- array(dim = c(20, 2, 97), dimnames = dimnames(raw.scu.female) )
for(i in 1:dim(raw.scu.male)[3]){
    reduce.male[,,i] <- raw.scu.male[-seq(2, 40, by = 2),,i]
}
for(i in 1:dim(raw.scu.female)[3]){
    reduce.female[,,i] <- raw.scu.female[-seq(2, 40, by = 2),,i]
}
raw.scu.male <- reduce.male
raw.scu.female <- reduce.female

##############################################################
## Define which landmarks will slide.
## slide.20 is for the 20 landmarks data.
slide.20 <- as.matrix(read.csv("./data/curveslide_20.csv"))

## Procrustes superposition step with sliding semilandmarks minimizing bending energy.
gen.male <- gpagen(raw.gen.male, ProcD = FALSE, curves = slide.20)
gen.female <- gpagen(raw.gen.female, ProcD = FALSE, curves = slide.20)
scu.male <- gpagen(raw.scu.male, ProcD = FALSE, curves = slide.20)
scu.female <- gpagen(raw.scu.female, ProcD = FALSE, curves = slide.20)
pro.male <- gpagen(raw.pro.male, ProcD = FALSE, curves = slide.20)
pro.female <- gpagen(raw.pro.female, ProcD = FALSE, curves = slide.20)
jug.male <- gpagen(raw.jug.male, ProcD = FALSE, curves = slide.20)
jug.female <- gpagen(raw.jug.female, ProcD = FALSE, curves = slide.20)

## Species mean value for each coordinate of the procrustes centralized data.
cord.gen.male <- to.mean.shape(gen.male)
cord.gen.female <- to.mean.shape(gen.female)
cord.scu.male <- to.mean.shape(scu.male)
cord.scu.female <- to.mean.shape(scu.female)
cord.pro.male <- to.mean.shape(pro.male)
cord.pro.female <- to.mean.shape(pro.female)
cord.jug.male <- to.mean.shape(jug.male)
cord.jug.female <- to.mean.shape(jug.female)

## Species mean value for centroid size
size.gen.male <- to.mean.centroid(gen.male)
size.gen.female <- to.mean.centroid(gen.female)
size.scu.male <- to.mean.centroid(scu.male)
size.scu.female <- to.mean.centroid(scu.female)
size.pro.male <- to.mean.centroid(pro.male)
size.pro.female <- to.mean.centroid(pro.female)
size.jug.male <- to.mean.centroid(jug.male)
size.jug.female <- to.mean.centroid(jug.female)

## Grand mean for the shape data:
mean.gen.male <- to.grand.mean(gen.male)
mean.gen.female <- to.grand.mean(gen.female)
mean.scu.male <- to.grand.mean(scu.male)
mean.scu.female <- to.grand.mean(scu.female)
mean.pro.male <- to.grand.mean(pro.male)
mean.pro.female <- to.grand.mean(pro.female)
mean.jug.male <- to.grand.mean(jug.male)
mean.jug.female <- to.grand.mean(jug.female)

## Format the phylogeny tip labels to be equal to the data:
tr$tip.label <- gsub("_"," ", tr$tip.label)

## Use the method of Graphen to get ultrametric tree:
tr.grafen <- compute.brlen(tr)

## Save data to downstream analyses:
save(cord.gen.male, cord.gen.female, cord.scu.male, cord.scu.female, cord.pro.male
   , cord.pro.female, cord.jug.male, cord.jug.female, tr, tr.grafen
   , size.gen.male, size.gen.female, size.scu.male, size.scu.female, size.pro.male
   , size.pro.female, size.jug.male, size.jug.female
   , mean.gen.male, mean.gen.female, mean.scu.male, mean.scu.female, mean.pro.male
   , mean.pro.female, mean.jug.male, mean.jug.female
   , file = "./data/Geniculata_data.RData")

## Now follow the 'Make_analysis.R' script for all the analysis made on the article.
