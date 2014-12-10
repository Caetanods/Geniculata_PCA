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
