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

sim.geomorpho <- function (phy, A, iter = 999){
	## This is a modification of the function 'compare.evol.rates' from the package
	##		'geomorph' by Dean Adams. Please cite the original package and correspondent
	##		articles. See 'help(compare.evol.rates)' for more info.

    if (length(dim(A)) == 3) {
        if (is.null(dimnames(A)[[3]])) {
            stop("Data matrix does not include taxa names as dimnames for 3rd dimension.")
        }
        x <- two.d.array(A) ## Make the data into "MorphoJ export" format.
    }
    if (length(dim(A)) == 2) {
        if (is.null(rownames(A))) {
            stop("Data matrix does not include taxa names as dimnames for rows.")
        }
        x <- A
    }
    if (any(is.na(A)) == T) {
        stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")
    }
    if (class(phy) != "phylo") 
        stop("tree must be of class 'phylo.'")
    ntaxa <- length(phy$tip.label)
    N <- nrow(x) ## This is the number of species.
    if (N != dim(x)[1]) {
        stop("Number of taxa in data matrix and tree are not not equal.")
    }
    if (length(match(rownames(x), phy$tip.label)) != N) 
        stop("Data matrix missing some taxa present on the tree.")
    if (length(match(phy$tip.label, rownames(x))) != N) 
        stop("Tree missing some taxa in the data matrix.")
    p <- ncol(x) ## This is the (number of landmarks * 2) + 1

    sigma.d <- function(phy, x, N) {
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
    
    sigmad.obs <- sigma.d(phy, x, ntaxa)
    ## return(list(sigma.d = sigmad.obs))  ## The original function ends here and do not do the simulations.
    }
    if (nlevels(gp) > 1) {
        sigmad.obs <- sigma.d(phy, x, ntaxa, gp)
        rate.mat <- diag(sigmad.obs$sigma.all, p)
        x.sim <- sim.char(phy, rate.mat, nsim = iter)
        sig.sim <- 1
        if (nlevels(gp) > 2) {
            gp.sig.sim <- array(1, dim = c(dim(sigmad.obs$sigma.gp)[1], 
                dim(sigmad.obs$sigma.gp)[1]))
        }
        rate.val <- rep(0, iter)
        for (ii in 1:iter) {
            sigmad.sim <- sigma.d(phy, x.sim[, , ii], ntaxa, 
                gp)
            sig.sim <- ifelse(sigmad.sim$ratio >= sigmad.obs$ratio, 
                sig.sim + 1, sig.sim)
            rate.val[ii] <- sigmad.sim$ratio
            if (nlevels(gp) > 2) {
                gp.sig.sim <- ifelse(sigmad.sim$sigma.gp >= sigmad.obs$sigma.gp, 
                  gp.sig.sim + 1, gp.sig.sim)
            }
        }
        sig.sim <- sig.sim/(iter + 1)
        if (nlevels(gp) > 2) {
            gp.sig.sim <- gp.sig.sim/(iter + 1)
            rownames(gp.sig.sim) <- colnames(gp.sig.sim) <- levels(gp)
            gp.sig.sim[lower.tri(gp.sig.sim)] <- NA
        }
        rate.val[iter + 1] = sigmad.obs$ratio
        hist(rate.val, 30, freq = TRUE, col = "gray", xlab = "SigmaD ratio")
        arrows(sigmad.obs$ratio, 50, sigmad.obs$ratio, 5, length = 0.1, 
            lwd = 2)
        if (nlevels(gp) > 2) {
            return(list(sigma.d = sigmad.obs$sigma.all, sigmad.all = sigmad.obs$sigma.d.all, 
                sigmad.ratio = sigmad.obs$ratio, pvalue = sig.sim, 
                pairwise.pvalue = gp.sig.sim))
        }
        else if (nlevels(gp) == 2) {
            return(list(sigma.d = sigmad.obs$sigma.all, sigmad.all = sigmad.obs$sigma.d.all, 
                sigmad.ratio = sigmad.obs$ratio, pvalue = sig.sim))
        }
    }
}
