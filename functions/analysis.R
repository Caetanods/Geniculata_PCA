physignal_extra <- function(A, phy, nsims=999, plot=TRUE){
    ## Function to calculate phylogenetic signal under a BM model using the distance methods by Adams.
    ## This function differ from 'physignal' from 'geomorph' because it calculates a distribution of
    ##    values for k given the observed data and tree. This distribution might deviate from 1 depending
    ##    on the correlation structure of the data.
    ## Function returns the vector with simulations results with tips shuffled and of simulations under
    ##    a mvBM model with the structure of evolutionary covariation estimated from the data.
    ## Function also returns the density of the observed value of K under both densities.
    ## A = geometric morphometics data.
    ## phy = phylogenetic tree.
    ## plot = Whether to plot the results.
    dt <- treedata(phy, to.matrix(A), sort=TRUE)
    out <- physignal(A=dt$data, phy=dt$phy, iter=nsims)
    get.Kstat(x=dt$data, phy=dt$phy)
    sigma <- ratematrix(dt$phy, dt$data)
    if( det(sigma) == 0 ) sigma <- nearPD( sigma )$mat
    sim <- sim.char(dt$phy, par=as.matrix(sigma), nsim=nsims)
    k.sim.bm <- sapply(1:999, function(x) get.Kstat(x=sim[,,x], phy=dt$phy) )
    k.sim.bm <- c(k.sim.bm, out$phy.signal)
    xlim <- range(c(out$random.K, k.sim.bm))
    d.k.sim <- ecdf(k.sim.bm)(out$phy.signal)
    d.k.shuffle <- ecdf(out$random.K)(out$phy.signal)
    den <- c(d.k.sim, d.k.shuffle)
    names(den) <- c("density_k_BM_sims","density_k_shuffle")
    if(plot){
        hist(out$random.K, xlim=xlim, freq=FALSE, border="white", col="gray", main = ""
             , xlab = "K values", )
        hist(k.sim.bm, xlim=xlim, add=TRUE, freq=FALSE, border="white", col="red")
        abline(v=out$phy.signal, col="blue", lwd = 2)
        legend(x="topright", legend=c("Shuffled tips","BM sims"), fill=c("gray","red")
             , border=c("white","white") )
    }
    return(list(density = den, k.sim.bm = k.sim.bm, k.shuffle = out$random.K))
}

get.Kstat <- function(x, phy){
    ## Function will calculate the K stat only, without performing simulations by shuffling the
    ##   tips of the phylogeny.
    ## x = A matrix, this can be generated with the function 'to.matrix'.
    ## phy = a 'phylo' object.
    x <- as.matrix(x)
    dt <- treedata(phy, x, sort=TRUE)
    N <- length(dt$phy$tip.label)
    phy.parts <- geomorph:::phylo.mat(dt$data, dt$phy)
    invC <- phy.parts$invC
    D.mat <- phy.parts$D.mat
    C <- phy.parts$C
    ones <- matrix(1, nrow(x), 1)
    a.obs <- colSums(invC) %*% x/sum(invC)
    distmat <- as.matrix(dist(rbind(as.matrix(dt$data), a.obs)))
    MSEobs.d <- sum(distmat[(1:N), (N + 1)]^2)
    dist.adj <- as.matrix(dist(rbind((D.mat %*% (x - (ones %*% a.obs))), 0)))
    MSE.d <- sum(dist.adj[(1:N), (N + 1)]^2)
    K.denom <- (sum(diag(C)) - N * solve(t(ones) %*% solve(C) %*% ones))/(N - 1)
    K.stat <- (MSEobs.d/MSE.d)/K.denom
    return( as.numeric(K.stat) )
}

sigma.d <- function(phy, x, N, p) {
    ## Function calculates the sigma^2_mult (Adams, 2014) for one group only.
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

pls.plot <- function(A1, A2, scores1, scores2, xlab, ylab, xlim=NULL, ylim=NULL){
    ## This function make a scatter.smooth plot of the PLS scores and add the
    ##   deformation grids. This is a custom function extracted from the
    ##   'pls.phylo' function in geomorph.
    A1.ref <- mshape(A1)
    pls1.min <- A1[, , which.min(scores1)]
    pls1.max <- A1[, , which.max(scores1)]
    A2.ref <- mshape(A2)
    pls2.min <- A2[, , which.min(scores2)]
    pls2.max <- A2[, , which.max(scores2)]
    par(mar = c(1, 1, 1, 1) + 0.1)
    split.screen(matrix(c(0.22, 1, 0.22, 1, 0.19, 0.39, 0, 
        0.19, 0.8, 1, 0, 0.19, 0, 0.19, 0.19, 0.39, 0, 0.19, 
        0.8, 1), byrow = T, ncol = 4))
    screen(1)
    scatter.smooth(scores1~scores2, axes = FALSE, main = ""
       , ylab = ylab, xlab = xlab, xlim = xlim, ylim = ylim, pch = 16)
    axis(side = 1); axis(side = 2)
    mtext(side = 1, xlab, line = 2.5)
    mtext(side = 2, ylab, line = 2.5)
    screen(2)
    geomorph:::tps(A1.ref, pls1.min, 20, sz = 0.7)
    screen(3)
    geomorph:::tps(A1.ref, pls1.max, 20, sz = 0.7)
    screen(4)
    geomorph:::tps(A2.ref, pls2.min, 20, sz = 0.7)
    screen(5)
    geomorph:::tps(A2.ref, pls2.max, 20, sz = 0.7)
    close.screen(all.screens = TRUE)
    par(mar = c(5.1, 4.1, 4.1, 2.1))
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

make.polygon <- function(times, mn.dispar, se.dispar, col){
    ## Function to plot the polygon of standard error of the mean for the dtt plot.
    ## Takes the times, mean disparity, standard error disparity, color.
    px <- c(times, times[length(times)]
          , times[(length(times) - 1):1]
          , times[1])
    py.down <- c(mn.dispar, mn.dispar[length(mn.dispar)] - se.dispar[length(se.dispar)]
               , mn.dispar[(length(mn.dispar) - 1):1] - se.dispar[(length(se.dispar) - 1):1]
               , mn.dispar[1])
    py.up <- c(mn.dispar, mn.dispar[length(mn.dispar)] + se.dispar[length(se.dispar)]
             , mn.dispar[(length(mn.dispar) - 1):1] + se.dispar[(length(se.dispar) - 1):1]
             , mn.dispar[1])
    polygon(px, py.down, col = col, border = col)
    polygon(px, py.up, col = col, border = col)
}

make.proc.dist <- function(x){
    ## Function to calculate pairwise partial Procrustes distances
    ##      of a geomorph object.
    ll <- dim(x)[3]
    res <- matrix(nrow = ll, ncol = ll)
    for(i in 1:ll ){
        for(j in 1:ll){
            res[i,j] <- proc.dist(x[,,i], x[,,j])
        }
    }
    colnames(res) <- rownames(res) <- dimnames(x)[[3]]
    return(res)
}

proc.dist <- function(x,y){
    ## Calculate the partial Procrustes distance between two shapes.
    ## x and y are two geomorph objects.
    sq <- (x - y)^2
    s.sq <- sum(apply(sq, 2, sum))
    rt.s.sq <- sqrt(s.sq)
    return(rt.s.sq)
}

to.procD.pgls <- function(shape, size, tree){
    ## Function perform the test of evolutionary allometry and calculates the P value for the multivariate regression using permutations.
    geo.dt <- geomorph.data.frame(coords=shape, size=size)
    res <- procD.pgls(f1=coords ~ size, phy=tree, data=geo.dt, iter=10000)
    return(res)
}

phylo.size.correct <- function(shape, mean.shape, size, tree){
    ## Function will perform the phylogenetic correction as described by Klingenberg and Marugan-Lobon, (2013) Systematic Biology.
    ## Calculate the coefficients of a regression of contrasts of centroid size and contrasts of shape. Calculate the deviance from the species mean shapes from the grand mean shape. Use the coefficients of the contrasts regression to calculate predicted values for the species mean shapes (using the raw centroid sizes). Returns the residuals of the predicted values subtracted from the deviance from the grand mean shape.
    m.shape <- to.matrix(shape)
    ll <- dim(shape)[3]
    dd <- dim(shape)[1]
    rownames(m.shape) <- dimnames(shape)[[3]]
    pic.shape <- sapply(1:ncol(m.shape), function(x) pic(x=m.shape[,x], phy=tree) )
    pic.size <- pic(x=size, phy=tree)
    regression.pic <- lm(pic.shape[,1:40] ~ pic.size -1)
    expected <- array( dim=c(dd, 2, ll ) )
    dev.mean <- array( dim=c(dd, 2, ll ) )
    residuals <- array( dim=c(dd, 2, ll ) )
    for(i in 1:ll){
        expected[,,i] <- matrix(size[i] * as.numeric(regression.pic$coefficients), nrow=20, ncol=2)
        dev.mean[,,i] <- shape[,,i] - mean.shape
        residuals[,,i] <- dev.mean[,,i] - expected[,,i]
    }
    dimnames( residuals ) <- list(NULL, NULL, names(size.gen.male))
    return(residuals)
}
