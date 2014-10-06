## Prepare morphometric data for analysis.
## Landmarks and procrustes analysis exported from MorphoJ.

library(geomorph)
library(geiger)

ind.male <- read.table("./data/Procrustes_male.txt", header = T, sep = "\t")[,-1]
ind.female <- read.table("./data/Procrustes_female.txt", header = T, sep = "\t")[,-1]

## Mean value for each coordinates of the procrustes centralized data.
cord.male <- aggregate(cord.male[,-1], by = list(cord.male$Sps), FUN = mean)
cord.female <- aggregate(cord.female[,-1], by = list(cord.female$Sps), FUN = mean)

## Phylogeny for comparative analysis.
tr <- read.tree("./data/tnt_geniculata.tre")

## Using the morphological integration analysis (Adams and Felice, 2014)
data(plethspecies)
Y.gpa <- gpagen(plethspecies$land)
phylo.pls(Y.gpa$coords[1:5,,],Y.gpa$coords[6:11,,],plethspecies$phy,iter=5)
