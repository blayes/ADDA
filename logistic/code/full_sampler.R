polyaGammaSampler <- function(yvec, xmat, ntrail = rep(1, length(yvec)),
                              mub0 = rep(0, ncol(xmat)), sigbInv0 = diag(0.01, ncol(xmat)),
                              niter = 1000) {
    
    library(pgdraw)
    ## Initialize.
    cts <- 0
    ndim <- ncol(xmat)
    nsample <- nrow(xmat)
    kappa <- yvec - ntrail / 2
    betas <- rep(0, ndim)
    psis <- xmat %*% betas
    w <- pgdraw(ntrail, psis)

    betasSamp <- matrix(0.0, niter, ndim)

    startTime <- proc.time()
    for (its in 1:niter) {
        if (its %% 500 == 0) cat("iter: ", its, "\n")
        ## draw w
        psis <- abs(drop(xmat %*% betas))
        w <- pgdraw(ntrail, psis)

        ## draw beta
        covBetasInv <- crossprod(xmat, w * xmat) + sigbInv0
        covBetas <- chol2inv(chol(covBetasInv))
        meanBetas <- covBetas %*% (crossprod(xmat, kappa) + sigbInv0 %*% mub0)
        betas <- meanBetas + crossprod(chol(covBetas), rnorm(ndim, 0, 1))

        cts <- cts + 1
        betasSamp[cts, ] <- betas

    }
    endTime = proc.time()

    list(
        'betas' = betasSamp,
        'time' = endTime[3] - startTime[3]
    )
}
