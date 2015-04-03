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

geo.comp.rates <- function (phy, A, B, plot = FALSE, iter = 999){
    ## Arguments:
    ## phy <- phylogeny of type 'phylo'
    ## A <- geomorph data. 2D morphometric data.
    ## B <- geomorph data. 2D morphometric data.
    
    ## This is a modification of the function 'compare.evol.rates' from the package
    ##	'geomorph' by Dean Adams. Please cite the original package and correspondent
    ##  articles. See 'help(compare.evol.rates)' for more info.

    ## The correlated null is generated using the R matrix estimated with 'ratematrix'
    ##  in geiger v.2 and scaled to the respective correlation matrix.
    ## When the number of traits is larger than the number of tips the R matrix
    ##  is converted into the nearest vcv matrix.
    ## The uncorrelated null is generated with the mean of the two sigma_mult values.
    ## Both for correlated and uncorrelated null models the function keep the relative
    ##  difference in the number of traits between matrices.
    ## See 'help(compare.evol.rates)' for more info.

    library(Matrix)	

    ## Check objects block:
    A <- two.d.array(A) ## Make the data into "MorphoJ export" format.
    B <- two.d.array(B) ## Make the data into "MorphoJ export" format.
    ntaxa <- length(phy$tip.label)
    A.p <- ncol(A) ## This is the (number of landmarks * 2) + 1
    B.p <- ncol(B) ## This is the (number of landmarks * 2) + 1

    ## Calculate sigma^2_mult for the matrix A and B.
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

    ## The uncorrelated case:
    rate.uncor.A <- diag(mean(sigmad.obs.A, sigmad.obs.B), A.p)
    rate.uncor.B <- diag(mean(sigmad.obs.A, sigmad.obs.B), B.p)
    A.sim.uncor <- sim.char(phy, rate.uncor.A, nsim = iter)
    B.sim.uncor <- sim.char(phy, rate.uncor.B, nsim = iter)

    ## The correlated case:
    rate.cor.A <- ratematrix(phy, A)
    rate.cor.B <- ratematrix(phy, B)
    ## Calculate nearest vcv matrix.
    if(A.p >= ntaxa) rate.cor.A <- nearPD(rate.cor.A)[[1]]
    if(B.p >= ntaxa) rate.cor.B <- nearPD(rate.cor.B)[[1]]
    ## Calculate the correlation matrix.
    rate.cor.A <- as.matrix(cov2cor(rate.cor.A))
    rate.cor.B <- as.matrix(cov2cor(rate.cor.B))
    A.sim.cor <- sim.char(phy, rate.cor.A, nsim = iter)
    B.sim.cor <- sim.char(phy, rate.cor.B, nsim = iter)

    sig.sim.uncor <- 1
    sigmad.sim.uncor.A <- rep(0, iter)
    sigmad.sim.uncor.B <- rep(0, iter)

    sig.sim.cor <- 1
    sigmad.sim.cor.A <- rep(0, iter)
    sigmad.sim.cor.B <- rep(0, iter)

    ## Calculate distribution of sigma values for correlated and uncorrelated.
    if(who == "A"){
        for (ii in 1:iter) {
            ## Null for uncorrelated:
            sigmad.sim.uncor.A[ii] <- sigma.d(phy, A.sim.uncor[, ,ii], ntaxa, A.p)
            sigmad.sim.uncor.B[ii] <- sigma.d(phy, B.sim.uncor[, ,ii], ntaxa, B.p)
            sig.sim.uncor <- ifelse(sigmad.sim.uncor.A[ii] / sigmad.sim.uncor.B[ii] >= sigmad.ratio, sig.sim.uncor + 1, sig.sim.uncor)
            
            ## Null for correlated:
            sigmad.sim.cor.A[ii] <- sigma.d(phy, A.sim.cor[, ,ii], ntaxa, A.p)
            sigmad.sim.cor.B[ii] <- sigma.d(phy, B.sim.cor[, ,ii], ntaxa, B.p)
            sig.sim.cor <- ifelse(sigmad.sim.cor.A[ii] / sigmad.sim.cor.B[ii] >= sigmad.ratio, sig.sim.cor + 1, sig.sim.cor)
        }
    sim.ratio.uncor <- sigmad.sim.uncor.A / sigmad.sim.uncor.B
    sim.ratio.cor <- sigmad.sim.cor.A / sigmad.sim.cor.B
    } else {
        for (ii in 1:iter) {
            ## Null for uncorrelated:
            sigmad.sim.uncor.A[ii] <- sigma.d(phy, A.sim.uncor[, ,ii], ntaxa, A.p)
            sigmad.sim.uncor.B[ii] <- sigma.d(phy, B.sim.uncor[, ,ii], ntaxa, B.p)
            sig.sim.uncor <- ifelse(sigmad.sim.uncor.B[ii] / sigmad.sim.uncor.A[ii] >= sigmad.ratio, sig.sim.uncor + 1, sig.sim.uncor)
            ## Null for correlated:
            sigmad.sim.cor.A[ii] <- sigma.d(phy, A.sim.cor[, ,ii], ntaxa, A.p)
            sigmad.sim.cor.B[ii] <- sigma.d(phy, B.sim.cor[, ,ii], ntaxa, B.p)
            sig.sim.cor <- ifelse(sigmad.sim.cor.B[ii] / sigmad.sim.cor.A[ii] >= sigmad.ratio, sig.sim.cor + 1, sig.sim.cor)
        }
        sim.ratio.uncor <- sigmad.sim.uncor.B / sigmad.sim.uncor.A
        sim.ratio.cor <- sigmad.sim.cor.B / sigmad.sim.cor.A
    }

    ## Calculate the p value for the Monte Carlo:
    p.value.cor <- sig.sim.cor / (iter + 1)
    p.value.uncor <- sig.sim.uncor / (iter + 1)

    ## Need to add the observed value to the distribution.
    sim.ratio.uncor <- append(sim.ratio.uncor, sigmad.ratio)
    sim.ratio.cor <- append(sim.ratio.cor, sigmad.ratio)

    if(plot == TRUE){        
        ## Make two histogram plots:
        par(mfrow = c(1,2))
        ## Uncorrelated data:
        hist(sim.ratio.uncor, 30, freq = TRUE, col = "gray", xlab = "SigmaD", main = "Uncorrelated null")
        arrows(sigmad.ratio, 50, sigmad.ratio, 5, length = 0.1, lwd = 2, col = "red")
        ifelse(who == "A", st <- "A / B", st <- "B / A")
        legend("topright", paste("sigma ratio: ", st, "\n", "p value: ", round(p.value.uncor, 3), sep="")
                                         , bty="n", text.col = "blue")
        ## Correlated data:
        hist(sim.ratio.cor, 30, freq = TRUE, col = "gray", xlab = "SigmaD", main = "Correlated null")
        arrows(sigmad.ratio, 50, sigmad.ratio, 5, length = 0.1, lwd = 2, col = "red")
        ifelse(who == "A", st <- "A / B", st <- "B / A")
        legend("topright", paste("sigma ratio: ", st, "\n", "p value: ", round(p.value.cor, 3), sep="")
                                         , bty="n", text.col = "blue")
    }

    ## Create return objects:
    obs <- c(sigmad.obs.A, sigmad.obs.B)
    names(obs) <- c("sigmaA","sigmaB")
    p.value <- c(p.value.uncor, p.value.cor)
    names(p.value) <- c("Uncorr","Corr")

    return(list(obs = obs, p.value = p.value, null.uncor = sim.ratio.uncor, null.cor = sim.ratio.cor))
}

geo.vec.comp.rates <- function (phy, A, B, plot = FALSE, iter = 999){
    ## Arguments:
    ## phy <- phylogeny of type 'phylo'
    ## A <- vector object. The scalar trait.
    ## B <- geomorph data. 2D morphometric data.
    ## This is a modification of the function 'compare.evol.rates' from the package
    ##	'geomorph' by Dean Adams. Please cite the original package and correspondent
    ##  articles. See 'help(compare.evol.rates)' for more info.

    ## This function compare rates of evolution of a scalar trait with a 2D morphometric data.
    ## The p.value for the difference is calculated using Monte Carlo simulations of the ratio
    ##  between the rates.

    ## The correlated null is generated using the R matrix estimated with 'ratematrix'
    ##  in geiger v.2 and scaled to the respective correlation matrix for the 2D morphometric trait.
    ## When the number of traits is larger than the number of tips the R matrix
    ##  is bent (converted) into the nearest positive-definite matrix.

    ## See 'help(compare.evol.rates)' for more info.

    library(Matrix)

    ## Check objects block:
    B <- two.d.array(B) ## Make the data into "MorphoJ export" format.
    ntaxa <- length(phy$tip.label)
    B.p <- ncol(B) ## This is the (number of landmarks * 2) + 1

    ## Calculate sigma^2_mult for the matrix B.
    sigmad.obs.B <- sigma.d(phy, B, ntaxa, B.p)
    ## Calculate sigma^2 for the vector A.
    fit <- fitContinuous(phy, A, model = "BM")
    sigmad.obs.A <- fit$opt$sigsq

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

    ## The scalar trait:
    A.sim <- sim.char(phy, sigmad.obs.A, model = "BM", nsim = iter)

    ## The uncorrelated case:
    rate.uncor.B <- diag(mean(sigmad.obs.A, sigmad.obs.B), B.p)
    B.sim.uncor <- sim.char(phy, rate.uncor.B, nsim = iter)

    ## The correlated case:
    rate.cor.B <- ratematrix(phy, B)
    ## Calculate nearest vcv matrix.
    if(B.p >= ntaxa) rate.cor.B <- nearPD(rate.cor.B)[[1]]
    ## Calculate the correlation matrix.
    rate.cor.B <- as.matrix(cov2cor(rate.cor.B))
    B.sim.cor <- sim.char(phy, rate.cor.B, nsim = iter)

    sig.sim.uncor <- 1
    sigmad.sim.uncor.B <- rep(0, iter)

    sig.sim.cor <- 1
    sigmad.sim.cor.B <- rep(0, iter)

    ## Calculate distribution of sigma values for correlated and uncorrelated.

    ## The vector estimate is constant:
    sigmad.sim.A <- sapply(seq(dim(A.sim)[3]), function(x) fitContinuous(phy, A.sim[,,x], model = "BM")$opt$sigsq )

    ## Now the correlated and uncorrelated. Guess I can make it faster somehow.
    if(who == "A"){
        for (ii in 1:iter) {
            ## Null for uncorrelated:
            sigmad.sim.uncor.B[ii] <- sigma.d(phy, B.sim.uncor[, ,ii], ntaxa, B.p)
            sig.sim.uncor <- ifelse(sigmad.sim.A[ii] / sigmad.sim.uncor.B[ii] >= sigmad.ratio, sig.sim.uncor + 1, sig.sim.uncor)
            
            ## Null for correlated:
            sigmad.sim.cor.B[ii] <- sigma.d(phy, B.sim.cor[, ,ii], ntaxa, B.p)
            sig.sim.cor <- ifelse(sigmad.sim.A[ii] / sigmad.sim.cor.B[ii] >= sigmad.ratio, sig.sim.cor + 1, sig.sim.cor)
        }
    sim.ratio.uncor <- sigmad.sim.A / sigmad.sim.uncor.B
    sim.ratio.cor <- sigmad.sim.A / sigmad.sim.cor.B
    } else {
        for (ii in 1:iter) {
            ## Null for uncorrelated:
            sigmad.sim.uncor.B[ii] <- sigma.d(phy, B.sim.uncor[, ,ii], ntaxa, B.p)
            sig.sim.uncor <- ifelse(sigmad.sim.uncor.B[ii] / sigmad.sim.A[ii] >= sigmad.ratio, sig.sim.uncor + 1, sig.sim.uncor)
            ## Null for correlated:
            sigmad.sim.cor.B[ii] <- sigma.d(phy, B.sim.cor[, ,ii], ntaxa, B.p)
            sig.sim.cor <- ifelse(sigmad.sim.cor.B[ii] / sigmad.sim.A[ii] >= sigmad.ratio, sig.sim.cor + 1, sig.sim.cor)
        }
        sim.ratio.uncor <- sigmad.sim.uncor.B / sigmad.sim.A
        sim.ratio.cor <- sigmad.sim.cor.B / sigmad.sim.A
    }

    ## Calculate the p value for the Monte Carlo:
    p.value.cor <- sig.sim.cor / (iter + 1)
    p.value.uncor <- sig.sim.uncor / (iter + 1)

    ## Need to add the observed value to the distribution.
    sim.ratio.uncor <- append(sim.ratio.uncor, sigmad.ratio)
    sim.ratio.cor <- append(sim.ratio.cor, sigmad.ratio)

    if(plot == TRUE){        
        ## Make two histogram plots:
        par(mfrow = c(1,2))
        ## Uncorrelated data:
        hist(sim.ratio.uncor, 30, freq = TRUE, col = "gray", xlab = "SigmaD", main = "Uncorrelated null")
        arrows(sigmad.ratio, 50, sigmad.ratio, 5, length = 0.1, lwd = 2, col = "red")
        ifelse(who == "A", st <- "A / B", st <- "B / A")
        legend("topright", paste("sigma ratio: ", st, "\n", "p value: ", round(p.value.uncor, 3), sep="")
                                         , bty="n", text.col = "blue")
        ## Correlated data:
        hist(sim.ratio.cor, 30, freq = TRUE, col = "gray", xlab = "SigmaD", main = "Correlated null")
        arrows(sigmad.ratio, 50, sigmad.ratio, 5, length = 0.1, lwd = 2, col = "red")
        ifelse(who == "A", st <- "A / B", st <- "B / A")
        legend("topright", paste("sigma ratio: ", st, "\n", "p value: ", round(p.value.cor, 3), sep="")
                                         , bty="n", text.col = "blue")
    }

    ## Create return objects:
    obs <- c(sigmad.obs.A, sigmad.obs.B)
    names(obs) <- c("sigmaA","sigmaB")
    p.value <- c(p.value.uncor, p.value.cor)
    names(p.value) <- c("Uncorr","Corr")

    return(list(obs = obs, p.value = p.value, null.uncor = sim.ratio.uncor, null.cor = sim.ratio.cor))
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

two.b.pls.plot <- function(x, xlab, ylab, cex = 1){
    ## Plot the results of two.b.pls.
    up <- round(max(x$x.scores), digits = 1) + 0.1
    plot(x$x.scores, x$y.scores, pch = 16, xlab = "", ylab = "", axes = FALSE, xlim = c(-up,up)
       , ylim = c(-up,up), cex = cex)
    axis(side = 1); axis(side = 2)
    mtext(side = 1, text = xlab, line = 2.5)
    mtext(side = 2, text = ylab, line = 2.5)
    abline(lm(x$y.scores~x$x.scores), lwd = cex)
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
