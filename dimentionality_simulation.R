## This is to simulate the behaviour of sigma mult and the sigma mult null model simulated with a diagonal matrix (i.e., no covariation) given two set of traits with different dimentionality.

## Start with a shape. Evolve the shape. Create a new dataset only by changing the dimentionality of the trait. Then estimate sigma mult, simulate the null, check if it works.
## Add the information of the R matrix in the simulation step. How to get the R matrix from the Adams method? Maybe from the covariance matrix of the correlated evolution paper.

library(geomorph)
library(geiger)
library(MASS)
library(phytools)

source("./functions/prepare-data.R")
source("./functions/simulations.R")

## Simulate the tree:
phy <- sim.bdtree(stop = "taxa", n = 100)
phy <- rescale(phy, model = "depth", 1)
plot(phy, direction = "upward"); axisPhylo(side = 2)

#############################################################################

## Generate R matrices with the same parameters:
## No difference in dimensionality:
R1 <- Posdef(10, rexp(10, 5/100))
R2 <- Posdef(10, rexp(10, 1/100))

## Use R matrices to simulate data, estimate sigma_mult and simulate sigma_mult
##     distributions under BM.
tr1 <- sim.char(phy, R1, model = "BM")[,,1]
tr1 <- to.geomorph(cbind(rownames(tr1),tr1))
tr2 <- sim.char(phy, R2, model = "BM")[,,1]
tr2 <- to.geomorph(cbind(rownames(tr2),tr2))
comp <- geo.comp.corr.rates(phy, tr1, tr2)

cc1 <- comp.adams.harmon(phy, tr1, 100)
cc2 <- comp.adams.harmon(phy, tr2, 100)

## We expect no difference between rates under both null models (correlated or not).
dif.uncor <- cc1$sigma.uncor.sim / cc2$sigma.uncor.sim
dif.cor <- cc1$sigma.cor.sim / cc2$sigma.cor.sim
pdf("Sims_cor_uncor_no_difference.pdf", width = 14)
par(mfrow = c(1,2))
hist(dif.uncor, breaks = 50, border = "blue", main = "Uncorrelated", xlab = "sigma_mult")
abline(v=cc1$sigma.obs / cc2$sigma.obs)
hist(dif.cor, breaks = 50, border = "blue", main = "Correlated", xlab = "sigma_mult")
abline(v=cc1$sigma.obs / cc2$sigma.obs)
dev.off()

##########################################
## This part here came from another script.
## Simulating under known R matrices and phylogenies:

phy <- sim.bdtree(stop="taxa")
phy <- rescale(phy, "depth", 1)

## Generate R matrices with the same parameters:
R1 <- Posdef(10)
R2 <- Posdef(10)
## Use R matrices to simulate data, estimate sigma_mult and simulate sigma_mult
##     distributions under BM.
tr1 <- sim.to.geo(phy, R1)
cc1 <- comp.adams.harmon(phy, tr1, 1000)
tr2 <- sim.to.geo(phy, R2)
cc2 <- comp.adams.harmon(phy, tr2, 1000)

## We expect no difference between rates under both null models (correlated or not).
dif.uncor <- cc1$sigma.uncor.sim / cc2$sigma.uncor.sim
dif.cor <- cc1$sigma.cor.sim / cc2$sigma.cor.sim
pdf("Sims_cor_uncor_no_difference.pdf", width = 14)
par(mfrow = c(1,2))
hist(dif.uncor, breaks = 50, border = "blue", main = "Uncorrelated", xlab = "sigma_mult")
abline(v=cc1$sigma.obs / cc2$sigma.obs)
hist(dif.cor, breaks = 50, border = "blue", main = "Correlated", xlab = "sigma_mult")
abline(v=cc1$sigma.obs / cc2$sigma.obs)
dev.off()
