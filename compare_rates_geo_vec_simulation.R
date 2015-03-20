## This is a simulation to test the performance of the function geo.vec.comp.rates().

## Load packages and functions:
library(geomorph)
library(geiger)
library(nsprcomp)
library(parallel)
library(Matrix)
library(grDevices)
source("./functions/analysis.R")
source("./functions/prepare-data.R")

## Load tree:
tr <- read.tree("./data/tnt_geniculata.tre")
tr.grafen <- compute.brlen(tr)

## Sigma rate for the simulation
rate <- 0.05

## Simulate geometric data:
mt <- diag(rate, nrow = 40, ncol = 40)
geo.sim <- sim.char(tr.grafen, par = mt, model = "BM")[,,1]
geo.mt <- data.frame(species = rownames(geo.sim), geo.sim, stringsAsFactors = FALSE)
geo.run <- to.geomorph(geo.mt)

## Simulate scalar data:
vec.run <- sim.char(tr.grafen, par = rate, model = "BM")[,,1]

## Simulation under no difference between the rates:
## sim <- geo.vec.comp.rates(tr.grafen, vec.run, geo.run, plot = TRUE)
load("data/geo.vec.sim.RData")
