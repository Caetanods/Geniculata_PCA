to.geomorph <- function(x){
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

to.matrix <- function(x){
    ## Function to get data in geomorph format and get a matrix.
    dd <- dim(x)
    spp <- dimnames(x)[[3]]
    dat <- matrix(ncol = dd[1]*dd[2])
    for(i in 1:dd[3] ){
        dat <- rbind(dat, as.numeric(t(x[,1,i]),t(x[,2,i])))
    }
    dat <- dat[-1,]
    rownames(dat) <- spp
    return(dat)
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

    return(list(obs = obs, p.value = p.value, null.uncor = sim.ratio.uncor, null.cor = sim.ratio.uncor))
}
