## n = 1e6, 1e7
## k = 10, 25, 50
## f = 20, 40, 60, 80
cmdArgs <- commandArgs(trailingOnly = TRUE)

nobsid <- as.integer(cmdArgs[1]) ## n = 1, 2
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
    cvtrain <- readRDS("/Shared/ssrivastva/dist_mov_len/data/train_ml_1e6.rds")
    partIds <- readRDS("/Shared/ssrivastva/dist_mov_len/data/part_ml.rds")

    train <- cvtrain[[cid]]
    y0 <- as.numeric(train$y > 3)
    x0 <- as.matrix(train$x)
    ntrail0 <- as.numeric(rep(1, length(y0)))

    splitIdx <- partIds[[1]][[nsubid]]
    subIdx <- split(1:length(y0), splitIdx)

    dataList <- list()
    for (ss in seq_along(subIdx)) {
        dataList[[ss]] <- list(x      = x0[subIdx[[ss]], ],
                               ntrail = ntrail0[subIdx[[ss]]]
                               )
    }

} else if (nobsid == 2) {
    cvtrain <- readRDS("/Shared/ssrivastva/dist_mov_len/data/train_ml_1e7.rds")
    partIds <- readRDS("/Shared/ssrivastva/dist_mov_len/data/part_ml.rds")

    train <- cvtrain[[cid]]
    y0 <- as.numeric(train$y > 3)
    x0 <- as.matrix(train$x)
    ntrail0 <- as.numeric(rep(1, length(y0)))

    splitIdx <- partIds[[2]][[nsubid]]
    subIdx <- split(1:length(y0), splitIdx)

    dataList <- list()
    for (ss in seq_along(subIdx)) {
        dataList[[ss]] <- list(x      = x0[subIdx[[ss]], ],
                               ntrail = ntrail0[subIdx[[ss]]]
                               )
    }
} else {
    print("peace")
}

library(parallel)
library(pgdraw)

numCores <- detectCores()

if (numCores >= nsubs) {
    cat("**** working with ", numCores, " cores, but we only need ",  nsubs, " cores ****", "\n")
    cl <- makeCluster(nsubs)
} else {
    cat("**** PROBLEM: working with ", numCores, " cores, but need ",  nsubs, " cores :PROBLEM  ***", "\n")
}

## Assign the variables appropriately in the workers' GlobalEnv. The 'recvd'
## variable is a consistency check to ensure we have succeeded.
sendData <- clusterApply(cl, dataList,
                         function (vars) {
                             x <<- vars$x
                             ntrail <<- vars$ntrail
                             recvd <<- TRUE
                         }
                         )

workerX <- clusterEvalQ(cl, x)

## we see that the below loop outputs TRUE, so we are fine
for (ii in 1:nsubs) {
    cat(all.equal(workerX[[ii]], lapply(dataList, function(ll) ll$x)[[ii]]), "\n")
}

samplePG <- function(xmat, betas, ntrail) {
    library(pgdraw)

    psis <- as.numeric(xmat %*% betas)
    w <- pgdraw(ntrail, psis)

    w
}

clusterExport(cl, "samplePG")

## START: intitialize a few variables before running our algorithm
ndim <- ncol(x0)
kappa <- y0 - ntrail0 / 2
betas <- rep(0.0, ndim)
psis <- rep(0.0, nrow(x0))
w <- pgdraw(ntrail0, psis)

mub0 <- rep(0, ndim)
sigbInv0 <- diag(0.01, ndim)

nwait <- round(frac * nsubs)
niter <- 10000
betasSamp <- matrix(0.0, niter, ndim)
track <- rep(0, nsubs)

eps <- 0.1
## END: intitialize a few variables before running our algorithm

startTime <- proc.time()
for (its in 1:niter) {
    if (its %% 100 == 0) cat("iter: ", its, "\n")

    mup <- 1 - rbinom(1, 1, eps)


    waitIdx <- sort(sample(1:nsubs, nwait))
    wSampList <- clusterCall(cl, function(betaSamp) samplePG(x, betaSamp, ntrail), betaSamp = betas)

    if (mup == 1) {
        w[unlist(subIdx[waitIdx])] <- unlist(wSampList[waitIdx])
        ## draw beta
        covBetasInv <- crossprod(x0, w * x0) + sigbInv0
        covBetas <- chol2inv(chol(covBetasInv))
        meanBetas <- covBetas %*% (crossprod(x0, kappa) + sigbInv0 %*% mub0)
        betas <- meanBetas + crossprod(chol(covBetas), rnorm(ndim, 0, 1))

        ## updated the books!
        betasSamp[its, ] <- betas
        track[waitIdx] <- 1 + track[waitIdx]
    } else {
        w[unlist(subIdx)] <- unlist(wSampList)
        ## draw beta
        covBetasInv <- crossprod(x0, w * x0) + sigbInv0
        covBetas <- chol2inv(chol(covBetasInv))
        meanBetas <- covBetas %*% (crossprod(x0, kappa) + sigbInv0 %*% mub0)
        betas <- meanBetas + crossprod(chol(covBetas), rnorm(ndim, 0, 1))

        ## updated the books!
        betasSamp[its, ] <- betas
        track[1:nsubs] <- 1 + track[1:nsubs]
    }
}
endTime <- proc.time()

res <- list('betas' = betasSamp,
            'track' = track,
            'time' = endTime[3] - startTime[3]
            )

fname <- paste0("/Shared/ssrivastva/dist_mov_len/result/logistic/part_n_", nobsid, "_k_", nsubid, "_frac_", fid, "_rep_", cid, "_eps_0.1.rds")
saveRDS(res, fname)

cat("Closing workers \n")

stopCluster(cl)

cat("Done! \n")
