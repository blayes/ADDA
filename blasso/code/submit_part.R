## n = 100, 1000, 10000
## k = 10, 25, 50
## f = 20, 40, 60, 80
cmdArgs <- commandArgs(trailingOnly = TRUE)

nobsid <- as.integer(cmdArgs[1]) ## n = 1, 2, 3
nsubid <- as.integer(cmdArgs[2]) ## k = 1, 2, 3
id <- as.integer(cmdArgs[3]) ## used to determine f, replication

nfrac <- c(20, 40, 60, 80)
nrep <- 1:10
ff <- rep(nfrac, times = length(nrep))
cc <- rep(nrep, each = length(nfrac))

fid <- ff[id]
cid <- cc[id]

frac <- fid / 100 ## fraction of communicating workers
nsubs <- c(10, 25, 50)[nsubid]

## create the data subsets to be sent to the workers
if (nobsid == 1) {
    cvtrain <- readRDS("/Shared/ssrivastva/dist_blasso/data/data_n_1.rds")
    partIds <- readRDS("/Shared/ssrivastva/dist_blasso/data/part_n_1.rds")

    train <- cvtrain[[cid]]
    y <- as.numeric(train$y)
    x <- as.matrix(train$x)

    splitIdx <- partIds[[cid]][[nsubid]]
} else if (nobsid == 2) {
    cvtrain <- readRDS("/Shared/ssrivastva/dist_blasso/data/data_n_2.rds")
    partIds <- readRDS("/Shared/ssrivastva/dist_blasso/data/part_n_2.rds")

    train <- cvtrain[[cid]]
    y <- as.numeric(train$y)
    x <- as.matrix(train$x)

    splitIdx <- partIds[[cid]][[nsubid]]
} else {
    cvtrain <- readRDS("/Shared/ssrivastva/dist_blasso/data/data_n_3.rds")
    partIds <- readRDS("/Shared/ssrivastva/dist_blasso/data/part_n_3.rds")

    train <- cvtrain[[cid]]
    y <- as.numeric(train$y)
    x <- as.matrix(train$x)

    splitIdx <- partIds[[cid]][[nsubid]]
}

library(parallel)
library(mgcv)

numCores <- detectCores()

if (numCores >= nsubs) {
    cat("**** working with ", numCores, " cores, but we only need ",  nsubs, " cores ****", "\n")
    cl <- makeCluster(nsubs)
} else {
    cat("**** PROBLEM: working with ", numCores, " cores, but need ",  nsubs, " cores :PROBLEM  ***", "\n")
}


## START: intitialize a few variables before running our algorithm
ndim <- ncol(x)
nsample <- nrow(x)

betas <- runif(ndim)
tauVec <- runif(ndim)
sig <- runif(1)

nwait <- round(frac * nsubs)
niter <- 20000
track <- rep(0, nsubs)

lambdaSq <- 2
a0 <- 0
b0 <- 0

sampBetas <- vector("list", niter)
sampSig <- vector("list", niter)
## END: intitialize a few variables before running our algorithm
clusterExport(cl, "lambdaSq")

subIdx <- split(1:ndim, splitIdx)

parList <- list()
for (ss in seq_along(subIdx)) {
    parList[[ss]] <- list(betas = betas[subIdx[[ss]]],
                          sig = sig
                          )
}


## Assign the variables appropriately in the workers' GlobalEnv. The 'recvd'
## variable is a consistency check to ensure we have succeeded.
sendData <- clusterApply(cl, parList,
                         function (pars) {
                             betas <<- pars$betas
                             sig <<- pars$sig
                             recvd <<- TRUE
                         }
                         )

workerBetas <- clusterEvalQ(cl, betas)
workerPars <- clusterEvalQ(cl, list(sig, lambdaSq))

drawPartTau <- function(partBetas, sig, lambdaSq) {
    library(mgcv)

    partDim <- length(partBetas)
    partTaus <- numeric(partDim)

    for (pp in 1:partDim) {
        partTaus[pp] <- 1 / rig(1, mean = sqrt(lambdaSq) * sqrt(sig) / abs(partBetas[pp]),
                                scale = 1 / lambdaSq)
    }

    partTaus
}

clusterExport(cl, "drawPartTau")

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

tauList <- clusterCall(cl, function(pars) drawPartTau(betas, sig, lambdaSq))

for (ss in seq_along(subIdx)) {
    tauVec[subIdx[[ss]]] <- tauList[[ss]]
}

betaSig <- drawBetaSig(y, x, tauVec, a0, b0)
betas <- betaSig$beta
sig <- betaSig$sig

startTime <- proc.time()
for (its in 1:niter) {
    if (its %% 500 == 0) {
        cat("iter: ", its, "\n")
    }

    parList <- list()
    for (ss in seq_along(subIdx)) {
        parList[[ss]] <- list(betas = betas[subIdx[[ss]]],
                              sig = sig
                              )
    }

    sendData <- clusterApply(cl, parList,
                             function (pars) {
                                 betas <<- pars$betas
                                 sig <<- pars$sig
                                 recvd <<- TRUE
                             }
                             )

    waitIdx <- sort(sample(1:nsubs, nwait))

    ## draw taus
    tauList <- clusterCall(cl, function(pars) drawPartTau(betas, sig, lambdaSq))
    for (ww in waitIdx) {
        tauVec[subIdx[[ww]]] <- tauList[[ww]]
    }

    betaSig <- drawBetaSig(y, x, tauVec, a0, b0)
    betas <- betaSig$beta
    sig <- betaSig$sig

    ## updated the books!
    sampBetas[[its]] <- betas
    sampSig[[its]] <- sig
    track[waitIdx] <- 1 + track[waitIdx]
}
endTime <- proc.time()

res <- list('sigma' = unlist(sampSig),
            'betas' = do.call(rbind, sampBetas),
            'track' = track,
            'time' = endTime[3] - startTime[3]
            )

fname <- paste0("/Shared/ssrivastva/dist_blasso/result/part/part_n_", nobsid, "_k_", nsubid, "_frac_", fid, "_rep_", cid, ".rds")
saveRDS(res, fname)

cat("Closing workers \n")

stopCluster(cl)

cat("Done! \n")
