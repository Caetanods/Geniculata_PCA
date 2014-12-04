to.geomorph <- function(x){
    arow <- ( dim(x)[2]-1 )/2
    amm <- nrow(x)
    mm <- array(NA, dim=c(arow,2,amm))
    for(i in 1:amm){
        sp <- x[i,][-1]
        cx <- sp[seq(1,arow*2,by=2)]
        cy <- sp[seq(2,arow*2,by=2)]
        names(cx) <- NULL
        cx <- as.numeric(cx)
        names(cy) <- NULL
        cy <- as.numeric(cy)
        mm[,1,i] <- cx
        mm[,2,i] <- cy
    }
    dimnames(mm)[[3]] <- as.character(x[,1])
    return(mm)
}

to.mean.shape <- function(ind.coord){
    ## Function gets the output of 'gpagen' and calculates the mean shape per species.
    spp <- unique(dimnames(ind.coord$coords)[[3]])
    ind.per.spp <- lapply(spp, FUN = function(x) which(dimnames(ind.coord$coords)[[3]] == x) )
    land <- dim(ind.coord$coords[,,1])[1]
    coord <- array(dim = c(land,2,length(spp)), dimnames = list(NULL,NULL,spp))
    for(i in 1:length(spp)){
        coord[,,i] <- apply(ind.coord$coords[,,ind.per.spp[[i]]], c(1, 2), mean)
    }
    return(coord)
}

geo.comp.rates <- function (phy, A, B, iter = 999){
	## This is a modification of the function 'compare.evol.rates' from the package
	##		'geomorph' by Dean Adams. Please cite the original package and correspondent
	##		articles. See 'help(compare.evol.rates)' for more info.

	## Check objects block:
    A <- two.d.array(A) ## Make the data into "MorphoJ export" format.
	B <- two.d.array(B) ## Make the data into "MorphoJ export" format.
    ntaxa <- length(phy$tip.label)
    A.p <- ncol(A) ## This is the (number of landmarks * 2) + 1
    B.p <- ncol(B) ## This is the (number of landmarks * 2) + 1

	## Define the sigma.d function:
    sigma.d <- function(phy, x, N, p) {
        x <- prcomp(x)$x
        ones <- array(1, N)
        C <- vcv.phylo(phy)
        C <- C[rownames(x), rownames(x)]
        a.obs <- colSums(solve(C)) %*% x/sum(solve(C))
        eigC <- eigen(C)
        D.mat <- solve(eigC$vectors %*% diag(sqrt(eigC$values)) %*% 
            t(eigC$vectors))
        dist.adj <- as.matrix(dist(rbind((D.mat %*% (x - (ones %*% 
            a.obs))), 0)))
        vec.d2 <- dist.adj[N + 1, 1:N]^2
        sigma.d.all <- sum(vec.d2)/N/p
        return(sigma.d.all)
    }
    
    sigmad.obs.A <- sigma.d(phy, A, ntaxa, A.p)
	sigmad.obs.B <- sigma.d(phy, B, ntaxa, B.p)
	## Getting the faster rate:
	if(sigmad.obs.A >= sigmad.obs.B){
		sigmad.ratio <- sigmad.obs.A / sigmad.obs.B
		sigma.sim <- sigmad.obs.A
		who <- "A"
	} else {
		sigmad.ratio <- sigmad.obs.B / sigmad.obs.A
		sigma.sim <- sigmad.obs.B
		who <- "B"
	}
	
	## Simulate the null distribution of traits.
	## There is no difference between the set of traits.
    rate.A <- diag(sigma.sim, A.p)
	rate.B <- diag(sigma.sim, B.p)

	## Both rate.A and B are diag matrices. So the simulation here have covariance 0.
    A.sim <- sim.char(phy, rate.A, nsim = iter)
	B.sim <- sim.char(phy, rate.B, nsim = iter)

    sig.sim <- 1
    sigmad.sim.A <- rep(0, iter)
	sigmad.sim.B <- rep(0, iter)

	## Calculate distribution of sigma values.
	if(who == "A"){
		for (ii in 1:iter) {
			sigmad.sim.A[ii] <- sigma.d(phy, A.sim[, ,ii], ntaxa, A.p)
			sigmad.sim.B[ii] <- sigma.d(phy, B.sim[, ,ii], ntaxa, B.p)
			sig.sim <- ifelse(sigmad.sim.A[ii] / sigmad.sim.B[ii] >= sigmad.ratio, sig.sim + 1, sig.sim)
		}
	sim.ratio <- sigmad.sim.A / sigmad.sim.B
	} else {
		for (ii in 1:iter) {
			sigmad.sim.A[ii] <- sigma.d(phy, A.sim[, ,ii], ntaxa, A.p)
			sigmad.sim.B[ii] <- sigma.d(phy, B.sim[, ,ii], ntaxa, B.p)
			sig.sim <- ifelse(sigmad.sim.B[ii] / sigmad.sim.A[ii] >= sigmad.ratio, sig.sim + 1, sig.sim)
		}
	sim.ratio <- sigmad.sim.B / sigmad.sim.A
	}
	
	## Calculate the p value for the Monte Carlo:
	p.value <- sig.sim / (iter + 1)

	## Need to add the observed value to the distribution.
	sim.ratio <- append(sim.ratio, sigmad.ratio)
	
	## Make the plot:
    hist(sim.ratio, 30, freq = TRUE, col = "gray", xlab = "SigmaD")
    arrows(sigmad.ratio, 50, sigmad.ratio, 5, length = 0.1, lwd = 2, col = "red")
	ifelse(who == "A", st <- "A / B", st <- "B / A")
	legend("topright", paste("sigma ratio: ", st, "\n", "p value: ", round(p.value, 3), sep="")
					 , bty="n", text.col = "blue")

    return(list(sigmaA = sigmad.obs.A, sigmaB = sigmad.obs.B,
					null.sigmaA = sigmad.sim.A, null.sigmaB = sigmad.sim.B, p.value = p.value))
}
