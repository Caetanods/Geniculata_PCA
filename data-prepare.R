## Prepare morphometric data for analysis.
## Landmarks and procrustes analysis exported from MorphoJ.

library(geomorph)

ind.male <- read.table("./data/Procrustes_male.txt", header = T, sep = "\t")[,-1]
ind.female <- read.table("./data/Procrustes_female.txt", header = T, sep = "\t")[,-1]

cord.male <- aggregate(cord.male[,-1], by = list(cord.male$Sps), FUN = mean)
cord.female <- aggregate(cord.female[,-1], by = list(cord.female$Sps), FUN = mean)


