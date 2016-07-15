## This script make all analyses in the manuscript.
## Please refer to "Prepare_data.R" for raw data, procrustes superposition and definition of the
##        sliding-landmarks.

## Load packages and functions
library(geiger) ## Version 2.0.6
library(geomorph) ## Version 3.0.0
library(parallel)
library(Matrix) ## Version 1.2
source("./functions/analysis.R")
source("./functions/prepare-data.R")

## Get data. See 'Prepare_data.R':
load("./data/Geniculata_data.RData")

########################################################################
## Test of allometry. Here we will check whether evolutionary allometry is an important factor in the evolution of Stink Bugs genitalia.
## Results of this analysis show that some of the traits have shapes which have an important effect of evolutionary allometry. Both the male and femanle genitalia traits do not show any effect from the correlation with the centroid size. Thus, we will do the analysis of modularity and rates of evolution with and without correcting in function of size.

reg.tr.gen.male <- to.procD.pgls(shape = cord.gen.male , size = log(size.gen.male), tree = tr)
reg.tr.gen.female <- to.procD.pgls(shape = cord.gen.female, size = log(size.gen.female), tree = tr)
reg.tr.scu.male <- to.procD.pgls(shape = cord.scu.male, size = log(size.scu.male), tree = tr)
reg.tr.scu.female <- to.procD.pgls(shape = cord.scu.female, size = log(size.scu.female), tree = tr)
reg.tr.pro.male <- to.procD.pgls(shape = cord.pro.male, size = log(size.pro.male), tree = tr)
reg.tr.pro.female <- to.procD.pgls(shape = cord.pro.female, size = log(size.pro.female), tree = tr)
reg.tr.jug.male <- to.procD.pgls(shape = cord.jug.male, size = log(size.jug.male), tree = tr)
reg.tr.jug.female <- to.procD.pgls(shape = cord.jug.female, size = log(size.jug.female), tree = tr)

reg.tr.gen.male$aov.table[1,7]
reg.tr.gen.female$aov.table[1,7]
reg.tr.scu.male$aov.table[1,7]
reg.tr.scu.female$aov.table[1,7] ## effect
reg.tr.pro.male$aov.table[1,7]
reg.tr.pro.female$aov.table[1,7] ## effect
reg.tr.jug.male$aov.table[1,7]
reg.tr.jug.female$aov.table[1,7]

reg.grafen.gen.male <- to.procD.pgls(shape = cord.gen.male , size = log(size.gen.male), tree = tr.grafen)
reg.grafen.gen.female <- to.procD.pgls(shape = cord.gen.female, size = log(size.gen.female), tree = tr.grafen)
reg.grafen.scu.male <- to.procD.pgls(shape = cord.scu.male, size = log(size.scu.male), tree = tr.grafen)
reg.grafen.scu.female <- to.procD.pgls(shape = cord.scu.female, size = log(size.scu.female), tree = tr.grafen)
reg.grafen.pro.male <- to.procD.pgls(shape = cord.pro.male, size = log(size.pro.male), tree = tr.grafen)
reg.grafen.pro.female <- to.procD.pgls(shape = cord.pro.female, size = log(size.pro.female), tree = tr.grafen)
reg.grafen.jug.male <- to.procD.pgls(shape = cord.jug.male, size = log(size.jug.male), tree = tr.grafen)
reg.grafen.jug.female <- to.procD.pgls(shape = cord.jug.female, size = log(size.jug.female), tree = tr.grafen)

reg.grafen.gen.male$aov.table[1,7]
reg.grafen.gen.female$aov.table[1,7]
reg.grafen.scu.male$aov.table[1,7]
reg.grafen.scu.female$aov.table[1,7]
reg.grafen.pro.male$aov.table[1,7]
reg.grafen.pro.female$aov.table[1,7]
reg.grafen.jug.male$aov.table[1,7] ## effect
reg.grafen.jug.female$aov.table[1,7]

########################################################################
## Tests of correlated evolution using phylogenetic PLS (PGLS) of the male and female
##    genitalia shape coordinates. This analysis generates a correlation of the shape of the male and female genitalia taking into accoun the phylogeny. Different from the PIC, the PGLS maintains the tip data that can be related with the phenotype values in the phylogeny.

int.tr <- phylo.integration(cord.gen.male[,,mm], cord.gen.female[,,mm], phy=tr)
int.grafen <- phylo.integration(cord.gen.male[,,mm], cord.gen.female[,,mm], phy=tr.grafen)

########################################################################
## Analysis of modularity.

## The above results show that the correlation between the female and male genitalia is stronger
##     than what is observed among the somatic traits. However, in order to further check if the genitalia
##     traits are more correlated than the somatic traits we will use a more directed test. This test
##     calculates the covariation among the landmarks within the groups against the landmarks between groups.
##     If the covariation within groups is stronger than among groups we say that there are modularity in the
##     data. For this test we are going to use the genitalia traits against the somatic traits as the two
##     test groups. This tests is somewhat related to PLS (see Adams, 2016). So we are also going to use the
##     contrast data.

## Get and prepare the data for the somatic traits:
## Change the data from 'geomorph' to a matrix format.
mm.gen <- to.matrix(cord.gen.male)
mf.gen <- to.matrix(cord.gen.female)
mm.scu <- to.matrix(cord.scu.male)
mf.scu <- to.matrix(cord.scu.female)
mm.pro <- to.matrix(cord.pro.male)
mf.pro <- to.matrix(cord.pro.female)
mm.jug <- to.matrix(cord.jug.male)
mf.jug <- to.matrix(cord.jug.female)

## Use ultrametric tree to calculate phylogenetic contrasts:
pic.grafen.gen.m <- sapply(1:dim(mm.gen)[2], function(x) pic(mm.gen[,x], tr.grafen) )
pic.grafen.gen.f <- sapply(1:dim(mf.gen)[2], function(x) pic(mf.gen[,x], tr.grafen) )
pic.grafen.scu.m <- sapply(1:dim(mm.scu)[2], function(x) pic(mm.scu[,x], tr.grafen) )
pic.grafen.scu.f <- sapply(1:dim(mf.scu)[2], function(x) pic(mf.scu[,x], tr.grafen) )
pic.grafen.pro.m <- sapply(1:dim(mm.pro)[2], function(x) pic(mm.pro[,x], tr.grafen) )
pic.grafen.pro.f <- sapply(1:dim(mf.pro)[2], function(x) pic(mf.pro[,x], tr.grafen) )
pic.grafen.jug.m <- sapply(1:dim(mm.jug)[2], function(x) pic(mm.jug[,x], tr.grafen) )
pic.grafen.jug.f <- sapply(1:dim(mf.jug)[2], function(x) pic(mf.jug[,x], tr.grafen) )

## Use speciational tree to calculate phylogenetic contrasts:
pic.unit.gen.m <- sapply(1:dim(mm.gen)[2], function(x) pic(mm.gen[,x], tr) )
pic.unit.gen.f <- sapply(1:dim(mf.gen)[2], function(x) pic(mf.gen[,x], tr) )
pic.unit.scu.m <- sapply(1:dim(mm.scu)[2], function(x) pic(mm.scu[,x], tr) )
pic.unit.scu.f <- sapply(1:dim(mf.scu)[2], function(x) pic(mf.scu[,x], tr) )
pic.unit.pro.m <- sapply(1:dim(mm.pro)[2], function(x) pic(mm.pro[,x], tr) )
pic.unit.pro.f <- sapply(1:dim(mf.pro)[2], function(x) pic(mf.pro[,x], tr) )
pic.unit.jug.m <- sapply(1:dim(mm.jug)[2], function(x) pic(mm.jug[,x], tr) )
pic.unit.jug.f <- sapply(1:dim(mf.jug)[2], function(x) pic(mf.jug[,x], tr) )

## Analysis for the ultrametric tree:
data.grafen <- cbind(rbind(pic.grafen.gen.m, pic.grafen.gen.f),
                     rbind(pic.grafen.scu.m, pic.grafen.scu.f),
                     rbind(pic.grafen.jug.m, pic.grafen.jug.f),
                     rbind(pic.grafen.pro.m, pic.grafen.pro.f)
                     )
## `data.graphen` is a matrix, with cols = landmarks and 24 rows, equal to the total number of individuals.
index <- rep("genitalia", times = 40) ## 40 is the number of traits. In our case, landmarks.
index <- c(index, rep("somatic", times = 40*3))
mod.grafen.test <- modularity.test(A=data.grafen, partition.gp = index, CI=TRUE, iter=10000)

## Now we repeat the analysis using the speciational tree instead of the graphen method ultrametric tree.
## Analysis for the ultrametric tree:
data.unit <- cbind(rbind(pic.unit.gen.m, pic.unit.gen.f),
                     rbind(pic.unit.scu.m, pic.unit.scu.f),
                     rbind(pic.unit.jug.m, pic.unit.jug.f),
                     rbind(pic.unit.pro.m, pic.unit.pro.f)
                     )
mod.unit.test <- modularity.test(A=data.unit, partition.gp = index, CI=TRUE, iter=10000)

## The analysis above will be performed again with the residuals of the PGLS linear regression on size. This analysis will check for the modularity between genitalia and somatic traits taking into account the evolutionary allometry in the shape.

## Use ultrametric tree to calculate phylogenetic contrasts with size corrected data:
ll <- 40
pic.size.grafen.gen.m <- sapply(1:ll, function(x) pic(reg.grafen.gen.male$residuals[,x], tr.grafen) )
pic.size.grafen.gen.f <- sapply(1:ll, function(x) pic(reg.grafen.gen.female$residuals[,x], tr.grafen) )
pic.size.grafen.scu.m <- sapply(1:ll, function(x) pic(reg.grafen.scu.male$residuals[,x], tr.grafen) )
pic.size.grafen.scu.f <- sapply(1:ll, function(x) pic(reg.grafen.scu.female$residuals[,x], tr.grafen) )
pic.size.grafen.pro.m <- sapply(1:ll, function(x) pic(reg.grafen.pro.male$residuals[,x], tr.grafen) )
pic.size.grafen.pro.f <- sapply(1:ll, function(x) pic(reg.grafen.pro.female$residuals[,x], tr.grafen) )
pic.size.grafen.jug.m <- sapply(1:ll, function(x) pic(reg.grafen.jug.male$residuals[,x], tr.grafen) )
pic.size.grafen.jug.f <- sapply(1:ll, function(x) pic(reg.grafen.jug.female$residuals[,x], tr.grafen) )

## Use speciational tree to calculate phylogenetic contrasts with size corrected data:
pic.size.tr.gen.m <- sapply(1:ll, function(x) pic(reg.tr.gen.male$residuals[,x], tr) )
pic.size.tr.gen.f <- sapply(1:ll, function(x) pic(reg.tr.gen.female$residuals[,x], tr) )
pic.size.tr.scu.m <- sapply(1:ll, function(x) pic(reg.tr.scu.male$residuals[,x], tr) )
pic.size.tr.scu.f <- sapply(1:ll, function(x) pic(reg.tr.scu.female$residuals[,x], tr) )
pic.size.tr.pro.m <- sapply(1:ll, function(x) pic(reg.tr.pro.male$residuals[,x], tr) )
pic.size.tr.pro.f <- sapply(1:ll, function(x) pic(reg.tr.pro.female$residuals[,x], tr) )
pic.size.tr.jug.m <- sapply(1:ll, function(x) pic(reg.tr.jug.male$residuals[,x], tr) )
pic.size.tr.jug.f <- sapply(1:ll, function(x) pic(reg.tr.jug.female$residuals[,x], tr) )

data.size.grafen <- cbind(rbind(pic.size.grafen.gen.m, pic.size.grafen.gen.f),
                          rbind(pic.size.grafen.scu.m, pic.size.grafen.scu.f),
                          rbind(pic.size.grafen.jug.m, pic.size.grafen.jug.f),
                          rbind(pic.size.grafen.pro.m, pic.size.grafen.pro.f)
                          )
mod.size.grafen.test <- modularity.test(A=data.size.grafen, partition.gp = index, CI=TRUE, iter=10000)

## Now we repeat the analysis using the speciational tree instead of the graphen method ultrametric tree.
## Analysis for the ultrametric tree:
data.size.tr <- cbind(rbind(pic.size.tr.gen.m, pic.size.tr.gen.f),
                     rbind(pic.size.tr.scu.m, pic.size.tr.scu.f),
                     rbind(pic.size.tr.jug.m, pic.size.tr.jug.f),
                     rbind(pic.size.tr.pro.m, pic.size.tr.pro.f)
                     )
mod.size.unit.test <- modularity.test(A=data.size.tr, partition.gp = index, CI=TRUE, iter=10000)

########################################################################
## Estimates of evolutionary rates

## Make estimates.
## Please uncomment and run the following lines to reproduce the analysis.
## Note that those are all Monte Carlo simulations. The p-values reperted in the manuscript can
##      vary slightly at each simulation. This variation should not change the results dicussed in
##      the article.

## Comparing rates between the male and female genitalia:
rate.index <- rep(c("one","two"), each=40)
comp.gen.tr.graphen <- compare.multi.evol.rates(A=cbind(to.matrix(cord.gen.male),to.matrix(cord.gen.female)),gp=rate.index,phy=tr.grafen,Subset=FALSE)
comp.gen.tr <- compare.multi.evol.rates(A=cbind(to.matrix(cord.gen.male),to.matrix(cord.gen.female)),gp=rate.index,phy=tr,Subset=FALSE)

## Rates between the male genitalia and the scutelum, pronotum and juga, under the grafen tree:
comp.scu.tr.m.grafen <- compare.multi.evol.rates(A=cbind(to.matrix(cord.gen.male),to.matrix(cord.scu.male)),gp=rate.index,phy=tr.grafen,Subset=FALSE)
comp.pro.tr.m.grafen <- compare.multi.evol.rates(A=cbind(to.matrix(cord.gen.male),to.matrix(cord.pro.male)),gp=rate.index,phy=tr.grafen,Subset=FALSE)
comp.jug.tr.m.grafen <- compare.multi.evol.rates(A=cbind(to.matrix(cord.gen.male),to.matrix(cord.jug.male)),gp=rate.index,phy=tr.grafen,Subset=FALSE)
## Rates between the male genitalia and the scutelum, pronotum and juga, under the unit tree:
comp.scu.tr.m <- compare.multi.evol.rates(A=cbind(to.matrix(cord.gen.male),to.matrix(cord.scu.male)),gp=rate.index,phy=tr,Subset=FALSE)
comp.pro.tr.m <- compare.multi.evol.rates(A=cbind(to.matrix(cord.gen.male),to.matrix(cord.pro.male)),gp=rate.index,phy=tr,Subset=FALSE)
comp.jug.tr.m <- compare.multi.evol.rates(A=cbind(to.matrix(cord.gen.male),to.matrix(cord.jug.male)),gp=rate.index,phy=tr,Subset=FALSE)

## Rates between the female genitalia and the scutelum, pronotum and juga, under the grafen tree:
comp.scu.tr.f.grafen <- compare.multi.evol.rates(A=cbind(to.matrix(cord.gen.female),to.matrix(cord.scu.female)),gp=rate.index,phy=tr.grafen,Subset=FALSE)
comp.pro.tr.f.grafen <- compare.multi.evol.rates(A=cbind(to.matrix(cord.gen.female),to.matrix(cord.pro.female)),gp=rate.index,phy=tr.grafen,Subset=FALSE)
comp.jug.tr.f.grafen <- compare.multi.evol.rates(A=cbind(to.matrix(cord.gen.female),to.matrix(cord.jug.female)),gp=rate.index,phy=tr.grafen,Subset=FALSE)
## Rates between the female genitalia and the scutelum, pronotum and juga, under the unit tree:
comp.scu.tr.f <- compare.multi.evol.rates(A=cbind(to.matrix(cord.gen.female),to.matrix(cord.scu.female)),gp=rate.index,phy=tr,Subset=FALSE)
comp.pro.tr.f <- compare.multi.evol.rates(A=cbind(to.matrix(cord.gen.female),to.matrix(cord.pro.female)),gp=rate.index,phy=tr,Subset=FALSE)
comp.jug.tr.f <- compare.multi.evol.rates(A=cbind(to.matrix(cord.gen.female),to.matrix(cord.jug.female)),gp=rate.index,phy=tr,Subset=FALSE)

## Now repeat the analyses with the residuals of the PGLS regression with size (centroid size).

## Compare rates of evolution of male and female controled by centroid size.
comp.size.gen.tr <- compare.multi.evol.rates(A=cbind(reg.tr.gen.male$residuals,reg.tr.gen.female$residuals),gp=rate.index,phy=tr,Subset=FALSE)
comp.size.gen.grafen <- compare.multi.evol.rates(A=cbind(reg.grafen.gen.male$residuals,reg.grafen.gen.female$residuals),gp=rate.index,phy=tr.grafen,Subset=FALSE)

## Rates between the male genitalia and the scutelum, pronotum and juga, under the unit tree:
comp.size.scu.m.tr <- compare.multi.evol.rates(A=cbind(reg.tr.gen.male$residuals,reg.tr.scu.female$residuals),gp=rate.index,phy=tr,Subset=FALSE)
comp.size.pro.m.tr <- compare.multi.evol.rates(A=cbind(reg.tr.gen.male$residuals,reg.tr.pro.female$residuals),gp=rate.index,phy=tr,Subset=FALSE)
comp.size.jug.m.tr <- compare.multi.evol.rates(A=cbind(reg.tr.gen.male$residuals,reg.tr.jug.female$residuals),gp=rate.index,phy=tr,Subset=FALSE)

## Rates between the male genitalia and the scutelum, pronotum and juga, under the grafen tree:
comp.size.scu.m.grafen <- compare.multi.evol.rates(A=cbind(reg.grafen.gen.male$residuals,reg.grafen.scu.female$residuals),gp=rate.index,phy=tr.grafen,Subset=FALSE)
comp.size.pro.m.grafen <- compare.multi.evol.rates(A=cbind(reg.grafen.gen.male$residuals,reg.grafen.pro.female$residuals),gp=rate.index,phy=tr.grafen,Subset=FALSE)
comp.size.jug.m.grafen <- compare.multi.evol.rates(A=cbind(reg.grafen.gen.male$residuals,reg.grafen.jug.female$residuals),gp=rate.index,phy=tr.grafen,Subset=FALSE)

## Rates between the female genitalia and the scutelum, pronotum and juga, under the unit tree:
comp.size.scu.f.tr <- compare.multi.evol.rates(A=cbind(reg.tr.gen.female$residuals,reg.tr.scu.female$residuals),gp=rate.index,phy=tr,Subset=FALSE)
comp.size.pro.f.tr <- compare.multi.evol.rates(A=cbind(reg.tr.gen.female$residuals,reg.tr.pro.female$residuals),gp=rate.index,phy=tr,Subset=FALSE)
comp.size.jug.f.tr <- compare.multi.evol.rates(A=cbind(reg.tr.gen.female$residuals,reg.tr.jug.female$residuals),gp=rate.index,phy=tr,Subset=FALSE)

## Rates between the female genitalia and the scutelum, pronotum and juga, under the unit tree:
comp.size.scu.f.grafen <- compare.multi.evol.rates(A=cbind(reg.grafen.gen.female$residuals,reg.grafen.scu.female$residuals),gp=rate.index,phy=tr.grafen,Subset=FALSE)
comp.size.pro.f.grafen <- compare.multi.evol.rates(A=cbind(reg.grafen.gen.female$residuals,reg.grafen.pro.female$residuals),gp=rate.index,phy=tr.grafen,Subset=FALSE)
comp.size.jug.f.grafen <- compare.multi.evol.rates(A=cbind(reg.grafen.gen.female$residuals,reg.grafen.jug.female$residuals),gp=rate.index,phy=tr.grafen,Subset=FALSE)

## ## Save all the workspace:
save.image("./data/all_analyses_results.RData")
