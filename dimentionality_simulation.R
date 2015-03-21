## This is to simulate the behaviour of sigma mult and the sigma mult null model simulated
##         with a diagonal matrix (i.e., no covariation) given two set of traits with
##         different dimentionality.

## We start with a shape. Evolve the shape. Create a new dataset only by changing
##         the dimentionality of the trait. Then estimate sigma mult, simulate the null
##         and check whether the method gives the expected result.

library(geomorph)
library(geiger)
library(MASS)
library(phytools)

source("./functions/prepare-data.R")
source("./functions/analysis.R")
source("./functions/simulations.R")

## Simulate the tree:
phy <- sim.bdtree(stop = "taxa", n = 50)
phy <- rescale(phy, model = "depth", 1)
plot(phy, direction = "upward"); axisPhylo(side = 2)

#############################################################################

## Generate R matrices. Different rates and dimensionality.

R1 <- Posdef(10, rexp(10, 2))

#####################################################
## Here we make simulations using diagonal matrices. This R matrices do not have a realistic
##     covariance structure, however, it is difficult to evaluate the expected difference
##     in rates of evolution under BM for more realistic covariance matrices since
##     the structure of covariance may result in differences of rates between matrices
##     that were simulated to be the same.
## The results of those simulations should hold if one include realistic R matrices that
##     are known to have a given evolutionary rate relationship.

## Doing simulations with diagonal matrices:
Rdiag1 <- diag(0.05, 10, 10)
Rdiag2 <- diag(0.10, 10, 10)
Rdiag3 <- diag(0.05, 20, 20)
Rdiag4 <- diag(0.10, 20, 20)

## Simulate the traits under a BM model:
trd1 <- sim.char(phy, Rdiag1, model = "BM")[,,1]
trd2 <- sim.char(phy, Rdiag1, model = "BM")[,,1]
trd3 <- sim.char(phy, Rdiag2, model = "BM")[,,1]
trd4 <- sim.char(phy, Rdiag3, model = "BM")[,,1]
trd5 <- sim.char(phy, Rdiag4, model = "BM")[,,1]

## Prepare the data to run the 'geo.comp.rates' function:
trd1 <- to.geomorph(cbind(rownames(trd1), trd1))
trd2 <- to.geomorph(cbind(rownames(trd2), trd2))
trd3 <- to.geomorph(cbind(rownames(trd3), trd3))
trd4 <- to.geomorph(cbind(rownames(trd4), trd4))
trd5 <- to.geomorph(cbind(rownames(trd5), trd5))

## Same rate and same dimensions.
diag.comp1 <- geo.comp.rates(phy, trd1, trd2, plot = TRUE)
diag.comp1$obs

## Same rate and different dimensions. (A is smaller)
diag.comp2 <- geo.comp.rates(phy, trd1, trd4, plot = TRUE)
diag.comp2$obs

## Different rates and same dimensions. (A is slower.)
diag.comp3 <- geo.comp.rates(phy, trd1, trd3, plot = TRUE)
diag.comp3$obs

## Different rates and different dimensions. (A is faster and smaller).
diag.comp4 <- geo.comp.rates(phy, trd3, trd4, plot = TRUE)
diag.comp4$obs

## Different rates and different dimensions. (A is faster and larger).
diag.comp5 <- geo.comp.rates(phy, trd5, trd1, plot = TRUE)
diag.comp5$obs
