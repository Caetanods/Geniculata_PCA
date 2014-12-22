comp.adams.harmon <- function (phy, A, times){
	## Simulation function to compare estimate of sigma_mult from the data
	##		with a 'sim.char' simulated data using the same phylo and a R matrix
	## 		derived from the vcv R matrix estimate of 'ratematrix'.
	## Returns sigma_mult for A and a vector with sigma_mult values for data simulated
	## 		under R.

	## Check objects block:
    A <- two.d.array(A) ## Make the data into "MorphoJ export" format.
    ntaxa <- length(phy$tip.label)
    A.p <- ncol(A) ## This is the (number of landmarks * 2) + 1

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

	## Get the estimate of sigma_mult    
    sigmad.obs.A <- sigma.d(phy, A, ntaxa, A.p)
	
	## Get the estimate of the R matrix:
	rate.A <- ratematrix(phy, A)
	rate.A <- cov2cor(rate.A)

	## Make a diag matrix following Adams:
	rate.uncor <- diag(sigmad.obs.A, A.p)

	## Simulat chars under the 
    A.sim <- sim.char(phy, rate.A, nsim = times)
	A.uncor <- sim.char(phy, rate.uncor, nsim = times)

    sigmad.sim.A <- rep(0, times)
	sigmad.uncor <- rep(0, times)

	## Calculate sigma_mult from simulated data.
	for (ii in 1:times) {
		sigmad.sim.A[ii] <- sigma.d(phy, A.sim[, ,ii], ntaxa, A.p)
		sigmad.uncor[ii] <- sigma.d(phy, A.uncor[, ,ii], ntaxa, A.p)
	}

    return(list(sigma.obs = sigmad.obs.A, sigma.cor.sim = sigmad.sim.A, sigma.uncor.sim = sigmad.uncor))
}

## Simulate vcv matrix:
Posdef <- function (n, ev = rexp(n, 1/100)) {
  ## Function for simulating positive definite covariance matrices
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}

phylo.pls.light <- function(A1, A2, phy, iter = 999){
	## This is a light version of phylo.pls which does not make default plots and dramatically
	##   improves speed by dropping series of if tests.
	## Warning: This version of the function does not provide warning messages and only
	##   works with 2D morphometric data.
    x <- two.d.array(A1)
    y <- two.d.array(A2)
    num.taxa.X <- nrow(x)
    namesX <- rownames(x)
    num.taxa.Y <- nrow(y)
    namesY <- rownames(y)
    if (is.null(namesX) == FALSE && is.null(namesY) == FALSE) {
        mtch.A <- namesX[is.na(match(namesX, namesY))]
    }
    mtch.B <- namesX[is.na(match(namesX, phy$tip.label))]
    data.all <- cbind(x, y)
    Nspec <- nrow(x)
    C <- vcv.phylo(phy, anc.nodes = FALSE)
    C <- C[rownames(y), rownames(y)]
    x <- x[rownames(y), ]
    invC <- solve(C)
    one <- matrix(1, Nspec, 1)
    a <- t(t(one) %*% invC %*% data.all) * sum(sum(invC))^-1
    R <- t(data.all - one %*% t(a)) %*% invC %*% (data.all - 
        one %*% t(a)) * (Nspec - 1)^-1
    R12 <- R[1:dim(x)[2], (dim(x)[2] + 1):(dim(x)[2] + dim(y)[2])]
    pls <- svd(R12)
    U <- pls$u
    V <- pls$v
    eigC <- eigen(C)
    D.mat <- solve(eigC$vectors %*% diag(sqrt(eigC$values)) %*% 
        t(eigC$vectors))
    Phy.X <- D.mat %*% (data.all - one %*% t(a))
    x.phy <- Phy.X[, c(1:dim(x)[2])]
    y.phy <- Phy.X[, c((dim(x)[2] + 1):(dim(x)[2] + dim(y)[2]))]
    XScores <- x.phy %*% U
    YScores <- y.phy %*% V
    pls.obs <- cor(XScores[, 1], YScores[, 1])
    P.val <- 1
    pls.val <- rep(0, iter)
    for (ii in 1:iter) {
        y.r <- y[sample(nrow(y)), ]
        data.all.r <- cbind(x, y.r)
        a.r <- t(t(one) %*% invC %*% data.all.r) * sum(sum(invC))^-1
        R.r <- t(data.all.r - one %*% t(a.r)) %*% invC %*% (data.all.r - 
            one %*% t(a.r)) * (Nspec - 1)^-1
        R12.r <- R.r[1:dim(x)[2], (dim(x)[2] + 1):(dim(x)[2] + 
            dim(y.r)[2])]
        pls.r <- svd(R12.r)
        U.r <- pls.r$u
        V.r <- pls.r$v
        Phy.X.r <- D.mat %*% (data.all.r - one %*% t(a.r))
        x.phy.r <- Phy.X.r[, c(1:dim(x)[2])]
        y.phy.r <- Phy.X.r[, c((dim(x)[2] + 1):(dim(x)[2] + dim(y.r)[2]))]
        XScores.r <- x.phy.r %*% U.r[, 1]
        YScores.r <- y.phy.r %*% V.r[, 1]
        pls.r <- cor(XScores.r, YScores.r)
        pls.val[ii] <- pls.r
        P.val <- ifelse(pls.r >= pls.obs, P.val + 1, P.val)
    }
    pls.val[iter + 1] = pls.obs
    P.val <- P.val/(iter + 1)
    for (i in 1:iter) {
        y.r <- y[sample(nrow(y)), ]
        XY.vcv.r <- cov(cbind(x, y.r))
        S12.r <- XY.vcv.r[1:dim(x)[2], (dim(x)[2] + 1):(dim(x)[2] + 
            dim(y.r)[2])]
        S21.r <- t(S12.r)
        S11.r <- XY.vcv.r[1:dim(x)[2], 1:dim(x)[2]]
        S22.r <- XY.vcv.r[(dim(x)[2] + 1):(dim(x)[2] + dim(y.r)[2]), 
            (dim(x)[2] + 1):(dim(x)[2] + dim(y.r)[2])]
        pls.r <- svd(S12.r)
        U.r <- pls.r$u
        V.r <- pls.r$v
        XScores.r <- x %*% U.r[, 1]
        YScores.r <- y.r %*% V.r[, 1]
        PLS.r <- cor(XScores.r, YScores.r)
        pls.val[i] <- PLS.r
        P.val <- ifelse(PLS.r >= pls.obs, P.val + 1, P.val)
    }
    pls.val[iter + 1] = pls.obs
    P.val <- P.val/(iter + 1)
    return(list(`PLS Correlation` = pls.obs, pvalue = P.val, 
            `Block 1 PLS Scores` = XScores[, 1], `Block 2 PLS Scores` = YScores[, 
                1]))
}
