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

sim.geomorpho <- function (phy, A, B, iter = 999){
	## This is a modification of the function 'compare.evol.rates' from the package
	##		'geomorph' by Dean Adams. Please cite the original package and correspondent
	##		articles. See 'help(compare.evol.rates)' for more info.

	## Check objects block:
    A <- two.d.array(A) ## Make the data into "MorphoJ export" format.
	B <- two.d.array(B)
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
		who <- "A"
	} else {
		sigmad.ratio <- sigmad.obs.B / sigmad.obs.A
		who <- "B"
	}
    rate.A <- diag(sigmad.obs.A, A.p)
	rate.B <- diag(sigmad.obs.B, B.p)
    A.sim <- sim.char(phy, rate.A, nsim = iter)
	B.sim <- sim.char(phy, rate.B, nsim = iter)
    sig.sim <- 1
    sigmad.sim.A <- rep(0, iter)
	sigmad.sim.B <- rep(0, iter)

	## Do the simulations
    for (ii in 1:iter) {
        sigmad.sim.A[ii] <- sigma.d(phy, A.sim[, ,ii], ntaxa, A.p)
		sigmad.sim.B[ii] <- sigma.d(phy, B.sim[, ,ii], ntaxa, B.p)
	}

	ifelse(who == "A", sim.ratio <- sigmad.sim.A / sigmad.sim.B, 
					sim.ratio <- sigmad.sim.B / sigmad.sim.A) 

    hist(sim.ratio, 30, freq = TRUE, col = "gray", xlab = "SigmaD")
    arrows(sigmad.ratio, 50, sigmad.ratio, 5, length = 0.1, lwd = 2)

    return(list(sigmaA = sigmad.obs.A, sigmaB = sigmad.obs.B,
					sim.sigmaA = sigmad.sim.A, sim.sigmaB = sigmad.sim.B))
}
