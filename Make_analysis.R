## This script make all analyses in the manuscript.
## Please refer to "Prepare_data.R" for raw data, procrustes superposition and definition of the
##        sliding-landmarks.

## Load packages and functions
library(geiger) ## Version 2.0.6
library(geomorph) ## Version 2.1.7
library(parallel)
library(Matrix) ## Version 1.2
source("./functions/analysis.R")
source("./functions/prepare-data.R")

########################################################################
## Tests of correlated evolution:

## Get data. See 'Prepare_data.R':
load("./data/Geniculata_data.RData")

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

## Two Block PLS analysis for ultrametric and speciational trees:
## Please note that significance values are calculated from simulations, thus
##      values of p will be different each time.
## Results produced bellow are the ones reported in Table 1.

## Male genitalia vs. female genitalia:
pls.gen.grafen <- two.b.pls(pic.grafen.gen.m, pic.grafen.gen.f, verbose = TRUE)
pls.gen.unit <- two.b.pls(pic.unit.gen.m, pic.unit.gen.f, verbose = TRUE)
## Male genitalia vs. male scutellum, pronotum and juga:
pls.scu.m.grafen <- two.b.pls(pic.grafen.gen.m, pic.grafen.scu.m, verbose = TRUE)
pls.scu.m.unit <- two.b.pls(pic.unit.gen.m, pic.unit.scu.m, verbose = TRUE)
pls.pro.m.grafen <- two.b.pls(pic.grafen.gen.m, pic.grafen.pro.m, verbose = TRUE)
pls.pro.m.unit <- two.b.pls(pic.unit.gen.m, pic.unit.pro.m, verbose = TRUE)
pls.jug.m.grafen <- two.b.pls(pic.grafen.gen.m, pic.grafen.jug.m, verbose = TRUE)
pls.jug.m.unit <- two.b.pls(pic.unit.gen.m, pic.unit.jug.m, verbose = TRUE)
## Female genitalia vs. female scutellum, pronotum and juga:
pls.scu.f.grafen <- two.b.pls(pic.grafen.gen.f, pic.grafen.scu.f, verbose = TRUE)
pls.scu.f.unit <- two.b.pls(pic.unit.gen.f, pic.unit.scu.f, verbose = TRUE)
pls.pro.f.grafen <- two.b.pls(pic.grafen.gen.f, pic.grafen.pro.f, verbose = TRUE)
pls.pro.f.unit <- two.b.pls(pic.unit.gen.f, pic.unit.pro.f, verbose = TRUE)
pls.jug.f.grafen <- two.b.pls(pic.grafen.gen.f, pic.grafen.jug.f, verbose = TRUE)
pls.jug.f.unit <- two.b.pls(pic.unit.gen.f, pic.unit.jug.f, verbose = TRUE)

########################################################################
## Estimates of evolutionary rates

## Make estimates.
## Please uncomment and run the following lines to reproduce the analysis.
## Note that those are all Monte Carlo simulations. The p-values reperted in the manuscript can
##      vary slightly at each simulation. This variation should not change the results dicussed in
##      the article.

# ## Comparing rates between the male and female genitalia:
# comp.gen.graphen <- geo.comp.rates(tr.grafen, cord.gen.male, cord.gen.female, plot = FALSE)
# comp.gen.unit <- geo.comp.rates(tr, cord.gen.male, cord.gen.female, plot = FALSE)
# ## Rates between the male genitalia and the scutelum, pronotum and juga:
# comp.scu.m.graphen <- geo.comp.rates(tr.grafen, cord.gen.male, cord.scu.male, plot = FALSE)
# comp.scu.m.unit <- geo.comp.rates(tr, cord.gen.male, cord.scu.male, plot = FALSE)
# comp.pro.m.graphen <- geo.comp.rates(tr.grafen, cord.gen.male, cord.pro.male, plot = FALSE)
# comp.pro.m.unit <- geo.comp.rates(tr, cord.gen.male, cord.pro.male, plot = FALSE)
# comp.jug.m.graphen <- geo.comp.rates(tr.grafen, cord.gen.male, cord.jug.male, plot = FALSE)
# comp.jug.m.unit <- geo.comp.rates(tr, cord.gen.male, cord.jug.male, plot = FALSE)
# ## Rates between the female genitalia and the scutelum, pronotum and juga:
# comp.scu.f.graphen <- geo.comp.rates(tr.grafen, cord.gen.female, cord.scu.female, plot = FALSE)
# comp.scu.f.unit <- geo.comp.rates(tr, cord.gen.female, cord.scu.female, plot = FALSE)
# comp.pro.f.graphen <- geo.comp.rates(tr.grafen, cord.gen.female, cord.pro.female, plot = FALSE)
# comp.pro.f.unit <- geo.comp.rates(tr, cord.gen.female, cord.pro.female, plot = FALSE)
# comp.jug.f.graphen <- geo.comp.rates(tr.grafen, cord.gen.female, cord.jug.female, plot = FALSE)
# comp.jug.f.unit <- geo.comp.rates(tr, cord.gen.female, cord.jug.female, plot = FALSE)

## Save results:
## save(list = ls(pattern = "^comp.*"), file = "./data/results_comp_rates.RData")

## Load results:
load("./data/results_comp_rates.RData")

## Checking the results of the analysis:

## Male vs. female genitalia:
comp.gen.graphen$obs[1]/comp.gen.graphen$obs[2]
comp.gen.graphen$p.value
comp.gen.unit$obs[1]/comp.gen.unit$obs[2]
comp.gen.unit$p.value

## Male genitalia vs. somatic traits:
comp.scu.m.graphen$obs[1]/comp.scu.m.graphen$obs[2]
comp.scu.m.graphen$p.value
comp.scu.m.unit$obs[1]/comp.scu.m.unit$obs[2]
comp.scu.m.unit$p.value
comp.pro.m.graphen$obs[1]/comp.pro.m.graphen$obs[2]
comp.pro.m.graphen$p.value
comp.pro.m.unit$obs[1]/comp.pro.m.unit$obs[2]
comp.pro.m.unit$p.value
comp.jug.m.graphen$obs[1]/comp.jug.m.graphen$obs[2]
comp.jug.m.graphen$p.value
comp.jug.m.unit$obs[1]/comp.jug.m.unit$obs[2]
comp.jug.m.unit$p.value

## Female genitalia vs. somatic traits:
comp.scu.f.graphen$obs[1]/comp.scu.f.graphen$obs[2]
comp.scu.f.graphen$p.value
comp.scu.f.unit$obs[1]/comp.scu.f.unit$obs[2]
comp.scu.f.unit$p.value
comp.pro.f.graphen$obs[1]/comp.pro.f.graphen$obs[2]
comp.pro.f.graphen$p.value
comp.pro.f.unit$obs[1]/comp.pro.f.unit$obs[2]
comp.pro.f.unit$p.value
comp.jug.f.graphen$obs[1]/comp.jug.f.graphen$obs[2]
comp.jug.f.graphen$p.value
comp.jug.f.unit$obs[1]/comp.jug.f.unit$obs[2]
comp.jug.f.unit$p.value

## Save all the workspace:
save.image("./data/all_analyses_results.RData")
