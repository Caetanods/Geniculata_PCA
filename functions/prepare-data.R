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
