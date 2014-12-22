## Herein we are going to make simulations to test whether the morphology integration method of Adams and Felice (2014) have good behaviour under the conditions of our dataset.
## Our simulations replicate the procedures used by the authors.

## To do this we generated a evolutionary covariance matrix. Evolved traits under a multivariate Brownian motion model with equal rates (sigma^2 = 1.0). Then we divided the traits in two datasets with 40 traits each. For the uncorrelated simulations we generated 40 independent traits using a diagonal matrix.
## The number of traits in each block and the number of tips in the phylogeny is the same as our dataset.

## Load packages and functions
library(geomorph)
library(geiger)
library(parallel)
library(Matrix)
library(grDevices)
source("./functions/prepare-data.R")
source("./functions/simulations.R")

## Define parameters:
ntaxa <- 64
ntraits <- 30 ## Number of traits per block. Some even number.
nsim <- 200

## Simulate trees:
phy <- mclapply(1:nsim, FUN = function(x) sim.bdtree(stop = "taxa", n = ntaxa), mc.cores = 2)
phy <- mclapply(phy, FUN = function(x) rescale(x, model = "depth", 1), mc.cores = 2)
plot(phy[[1]], direction = "upward"); axisPhylo(side = 2)

##############################################
## No correlation simulation:
##############################################

## Prepare data:
uncorr.mat <- diag(1, nrow = ntraits)
tr1 <- mclapply(phy, FUN = function(x)
    sim.char(x, uncorr.mat, nsim = 1, model = "BM", root = runif(1,10,20)), mc.cores = 2)
tr2 <- mclapply(phy, FUN = function(x)
    sim.char(x, uncorr.mat, nsim = 1, model = "BM", root = runif(1,10,20)), mc.cores = 2)
geo.tr1 <- mclapply(1:nsim, FUN = function(x)
    to.geomorph(cbind(rownames(tr1[[x]][,,1]),tr1[[x]][,,1])), mc.cores = 2)
geo.tr2 <- mclapply(1:nsim, FUN = function(x)
    to.geomorph(cbind(rownames(tr2[[x]][,,1]),tr2[[x]][,,1])), mc.cores = 2)
pro.geo.tr1 <- mclapply(1:nsim, function(x) gpagen(geo.tr1[[x]], ShowPlot = FALSE)[[1]], mc.cores = 2)
pro.geo.tr2 <- mclapply(1:nsim, function(x) gpagen(geo.tr2[[x]], ShowPlot = FALSE)[[1]], mc.cores = 2)

## Make analysis:
res <- mclapply(1:nsim, function(x)
    phylo.pls.light(pro.geo.tr1[[x]], pro.geo.tr2[[x]], phy[[x]]), mc.cores = 2)

## Result analysis:
p.value <- sapply(1:nsim, function(x) res[[x]][[2]][1,1])
sum(p.value <= 0.05) / nsim

## Plot results:
par(mfrow = c(1,2))
plot(density(p.value))
abline(v = 0.05, col = "red")
qs <- quantile(p.value, probs = c(0.05,0.95))
mn <- mean(p.value)
abline(v = qs, lty = 3); abline(v = mn, col = "blue")
boxplot(p.value); abline(h = 0.05, col = "red")

##############################################
## With correlation simulation:
##############################################

R <- Posdef(80)
Rcor <- cov2cor(R)
