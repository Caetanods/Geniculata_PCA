## Example of simulation.

## Load packages
library(geomorph)
library(geiger)

## Define helping functions:
to.geomorph <- function(x){
    ## Function to put data in the correct format.
    arow <- ( dim(x)[2]-1 )/2
    amm <- nrow(x)
    mm <- array(NA, dim=c(arow,2,amm))
    for(i in 1:amm){
        sp <- x[i,][-1]
        cx <- sp[seq(1,arow*2,by=2)]
        cy <- sp[seq(2,arow*2,by=2)]
        names(cx) <- NULL
        names(cy) <- NULL
        cx <-as.numeric(as.matrix(cx))
        cy <- as.numeric(as.matrix(cy))
        cxy <- cbind(cx,cy)
        mm[,,i] <- cxy
    }
    dimnames(mm)[[3]] <- as.character(x[,1])
    return(mm)
}

## Define parameters:
ntaxa <- 64 ## Number of tips in the phylo.
ntraits <- 30 ## Number of traits per block. Need to be a even number.
nsim <- 200 ## Number of simulations.

## Simulate trees:
phy <- lapply(1:nsim, FUN = function(x) sim.bdtree(stop = "taxa", n = ntaxa) )
phy <- lapply(phy, FUN = function(x) rescale(x, model = "depth", 1) )

##############################################
## Simulation with no effect of correlation:
##############################################

## Generate R matrix as diagonal:
uncorr.mat <- diag(1, nrow = ntraits)
## Simulate two multidimensional traits under BM (tr1 and tr2):
tr1 <- lapply(phy, FUN = function(x) sim.char(x, uncorr.mat, nsim = 1, model = "BM") )
tr2 <- lapply(phy, FUN = function(x) sim.char(x, uncorr.mat, nsim = 1, model = "BM") )
## Format traits to simulate geometric morphometric data:
geo.tr1 <- lapply(1:nsim, FUN = function(x) to.geomorph(cbind(rownames(tr1[[x]][,,1]),tr1[[x]][,,1])) )
geo.tr2 <- lapply(1:nsim, FUN = function(x) to.geomorph(cbind(rownames(tr2[[x]][,,1]),tr2[[x]][,,1])) )
## Use 'gpagen' to simulate procrustes. Just to replicate the same procedure
##    applied to the empirical data.
pro.geo.tr1 <- lapply(1:nsim, function(x) gpagen(geo.tr1[[x]], ShowPlot = FALSE)[[1]] )
pro.geo.tr2 <- lapply(1:nsim, function(x) gpagen(geo.tr2[[x]], ShowPlot = FALSE)[[1]] )

## Make analysis:
res <- lapply(1:nsim, function(x) phylo.pls(pro.geo.tr1[[x]], pro.geo.tr2[[x]], phy[[x]]) )

## Get p-value and calculate percent of significance:
p.value <- sapply(1:nsim, function(x) res[[x]][[2]][1,1])
sum(p.value <= 0.05) / nsim
## ~ 0.5

## Plot results:
pdf("Example_sim.pdf", width = 14, height = 7)
par(mfrow = c(1,2))
plot(density(p.value), main = "Density of p.values")
abline(v = 0.05, col = "red")
qs <- quantile(p.value, probs = c(0.05,0.95,0.75))
mn <- mean(p.value)
abline(v = qs[-3], lty = 3); abline(v = mn, col = "blue")
legend(x = qs[3], y = 2.0, legend = c("0.05 level", "mean", "0.05; 0.95 quantiles"), lty = c(1,1,3), col = c("red","blue","black"), cex = 0.75)
boxplot(p.value, main = "Red line = 0.05"); abline(h = 0.05, col = "red")
dev.off()
