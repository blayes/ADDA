drawTau <- function(betas, sig, lambdaSq) {
    library(mgcv)
    ndim <- length(betas)
    taus <- numeric(ndim)

    for (pp in 1:ndim) {
        taus[pp] <- 1 / rig(1, mean = sqrt(lambdaSq) * sqrt(sig) / abs(betas[pp]),
                            scale = 1 / lambdaSq)
    }

    taus
}

drawBetaSig <- function(y, x, taus, a0, b0) {
    ndim <- ncol(x)
    nsamp <- nrow(x)

    A <- crossprod(x, x) + diag(1 / taus)
    Ainv <- chol2inv(chol(A))
    aa <- 0.5 * nsamp + a0
    bb <- 0.5 * (sum(y * y) - drop(crossprod(y, tcrossprod(x %*% Ainv, x) %*% y))) + b0

    sig <- 1 / rgamma(1, aa, bb)

    muBeta <- Ainv %*% crossprod(x, y)
    covBeta <- sig * Ainv

    list(sig = sig, beta = drop(muBeta + crossprod(chol(covBeta), rnorm(ndim))))
}

lassoSampler <- function(y, x, a0 = 0, b0 = 0,
                         niter = 10000) {

    y <- as.numeric(y)
    x <- as.matrix(x)

    sampBeta <- list()
    sampSigma <- list()

    ndim <- ncol(x)
    sigPrev <- runif(1)
    betaPrev <- runif(ndim)
    tauPrev <- runif(ndim)
    lambdaSq <- 1

    startTime <- proc.time()
    for (ii in 1:niter) {
        if ((ii %% 500) == 0) {
            cat("its: ", ii, "\n")
        }
        tauVec <- drawTau(betaPrev, sigPrev, lambdaSq)
        betaSig <- drawBetaSig(y, x, tauVec, a0, b0)
        betas <- betaSig$beta
        sig <- betaSig$sig

        sampBeta[[ii]] <- betas
        sampSigma[[ii]] <- sig

        sigPrev <- sig
        betaPrev <- betas
    }
    endTime <- proc.time()


    list('betas' = do.call(rbind, sampBeta),
         'sig' = unlist(sampSigma),
         'time' = endTime[3] - startTime[3]
         )
}
