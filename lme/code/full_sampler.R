## I step
sampleRanEff <- function (ylist, fixefList, ranefList, fixSamp, errVar, tmatTilde, gammaMat) {
    nsample <- length(fixefList)
    nranef <- ncol(ranefList[[1]])

    ranSampList <- vector("list", nsample)
    for (ii in 1:nsample) {
        tmp1 <- ranefList[[ii]] %*% gammaMat %*% tmatTilde # z g t
        wmat <- chol2inv(chol(tcrossprod(ranefList[[ii]] %*% gammaMat, tmp1) + errVar * diag(1, nrow(ranefList[[ii]]))))
        postRanVar <- tmatTilde - crossprod(tmp1, wmat %*% tmp1)
        postRanMean <- drop(crossprod(tmp1, wmat %*% (ylist[[ii]] - fixefList[[ii]] %*% fixSamp)))
        ranSampList[[ii]] <- postRanMean + drop(crossprod(chol(postRanVar), rnorm(length(postRanMean))))
    }

    ranSampList
}

## P step
samplePars <- function (ylist, fixefList, ranefList, ranSampList) {
    nsample <- length(fixefList) # m
    nfix <- ncol(fixefList[[1]]) # p
    nranef <- ncol(ranefList[[1]]) # q
    ndim <- nranef * nranef # dim (l)
    nobsv <- sum(sapply(fixefList, nrow)) # n
    
    xztxz <- matrix(0.0, nfix + ndim, nfix + ndim)
    xzty <- matrix(0.0, nfix + ndim, 1)
    for (ii in 1:nsample) {
        idz <- rep(1:nranef, times = nranef)
        zz <- ranefList[[ii]][ , idz, drop = FALSE]
        idc <- rep(1:nranef, each = nranef)
        cc <- matrix(ranSampList[[ii]][idc], ncol = length(idc), nrow = nrow(zz), byrow = TRUE)
        ztilde <- zz * cc
        xtilde <- cbind(fixefList[[ii]], ztilde)
        xztxz <- xztxz + crossprod(xtilde, xtilde)
        xzty <- xzty + crossprod(xtilde, ylist[[ii]])        
    }

    xztxzInv <- chol2inv(chol(xztxz))
    alphaHat <- drop(xztxzInv %*% xzty)

    rss <- 0.0
    for (ii in 1:nsample) {
        idz <- rep(1:nranef, times = nranef)
        zz <- ranefList[[ii]][ , idz, drop = FALSE]
        idc <- rep(1:nranef, each = nranef)
        cc <- matrix(ranSampList[[ii]][idc], ncol = length(idc), nrow = nrow(zz), byrow = TRUE)        
        ztilde <- zz * cc
        xtilde <- cbind(fixefList[[ii]], ztilde)
        rss <- rss + sum((ylist[[ii]] - drop(xtilde %*% alphaHat))^2)
    }

    errVar <- rss / rchisq(1, nobsv - nfix - ndim)
    
    alphas <- drop(alphaHat + crossprod(chol(errVar * xztxzInv), rnorm(nfix + ndim)))

    fixSamp <- alphas[1:nfix]        
    gmat <- matrix(alphas[-(1:nfix)], nranef, nranef)

    vv <- diag(0.0, nranef)
    for (ii in 1:nsample) {
        vv <- vv + tcrossprod(ranSampList[[ii]], ranSampList[[ii]])
    }
    
    tmatTilde <- chol2inv(chol(rWishart(1, nsample - nranef, solve(vv))[ , , 1]))
    tmat <- tcrossprod(gmat %*% tmatTilde, gmat)
    
    list(errVar = errVar, fixSamp = fixSamp, tmatSamp = tmat, tmatTildeSamp = tmatTilde, gmatSamp = gmat)
}

## van Dyk-Meng Sampler
lmeSampler <- function (ylist, fixefList, ranefList, 
                        beta0, tmat0, errVar0,
                        niter = 1000) {
    tmat <- tmat0
    tmatTilde <- tmat0
    errVar <- errVar0
    fixs <- beta0
    its <- 0
    nsample <- length(ranefList)
    nranef <- ncol(ranefList[[1]])
    gmat <- diag(1, nranef)
    
    sampTmat <- list()
    sampErrVar <- list()
    sampBeta <- list()

    startTime <- proc.time()
    while (its < niter) {
        its <- its + 1

        if (its %% 100 == 0) cat("mixef iter: ", its, "\n")

        rans <- sampleRanEff(ylist, fixefList, ranefList, fixs, errVar, tmatTilde, gmat)
        pars <- samplePars(ylist, fixefList, ranefList, rans)
        
        fixs <- pars$fixSamp
        tmat <- pars$tmatSamp
        tmatTilde <- pars$tmatTildeSamp
        gmat <- pars$gmatSamp
        errVar <- pars$errVar
        
        sampErrVar[[its]] <- errVar
        sampBeta[[its]] <- fixs
        sampTmat[[its]] <- tmat
        
    }
    endTime <- proc.time()

    list(
        sigma = unlist(sampErrVar),
        tmat = sampTmat,
        betas = do.call(rbind, sampBeta),
        time = endTime[3] - startTime[3]
    )
}

