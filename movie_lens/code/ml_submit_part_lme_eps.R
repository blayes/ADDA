## NOTE: the original file was lost in a backup fudge by HPC!
## n = 1 (< 1 mil), 2 (> 1 mil)
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
    cvtrain <- readRDS("/Shared/ssrivastva/dist_mov_len/data/train_ml_100k.rds")
    partIds <- readRDS("/Shared/ssrivastva/dist_mov_len/data/part_ids_100k.rds")

    train <- cvtrain[[cid]]
    subIdx <- partIds[[cid]][[nsubid]]

    dataList <- list()
    for (ss in seq_along(subIdx)) {
        idx <- which(train$group %in% subIdx[[ss]])
        group <- train$group[idx]
        grpLbl <- sort(unique(group))
        ngroup <- length(grpLbl)
        ranefList0 <- list()
        fixefList0 <- list()
        ylist0 <- list()
        grpIdx0 <- list()
        for (gg in 1:ngroup) {
            grpIdx0[[gg]] <- which(train$group == grpLbl[gg])
            fixefList0[[gg]] <- train$x[grpIdx0[[gg]], , drop = FALSE]
            ranefList0[[gg]] <- train$z[grpIdx0[[gg]], , drop = FALSE]
            ylist0[[gg]] <- train$y[grpIdx0[[gg]]]
        }

        dataList[[ss]] <- list(x = fixefList0,
                               z = ranefList0,
                               y = ylist0,
                               users = subIdx[[ss]]
                               )
    }
    nranef0 <- ncol(dataList[[1]]$z[[1]])
    nfixef0 <- ncol(dataList[[1]]$x[[1]])
} else if (nobsid == 2) {
    cvtrain <- readRDS("/Shared/ssrivastva/dist_mov_len/data/train_ml_1000k.rds")
    partIds <- readRDS("/Shared/ssrivastva/dist_mov_len/data/part_ids_1000k.rds")

    train <- cvtrain[[cid]]
    subIdx <- partIds[[cid]][[nsubid]]

    dataList <- list()
    for (ss in seq_along(subIdx)) {
        idx <- which(train$group %in% subIdx[[ss]])
        group <- train$group[idx]
        grpLbl <- sort(unique(group))
        ngroup <- length(grpLbl)
        ranefList0 <- list()
        fixefList0 <- list()
        ylist0 <- list()
        grpIdx0 <- list()
        for (gg in 1:ngroup) {
            grpIdx0[[gg]] <- which(train$group == grpLbl[gg])
            fixefList0[[gg]] <- train$x[grpIdx0[[gg]], , drop = FALSE]
            ranefList0[[gg]] <- train$z[grpIdx0[[gg]], , drop = FALSE]
            ylist0[[gg]] <- train$y[grpIdx0[[gg]]]
        }

        dataList[[ss]] <- list(x = fixefList0,
                               z = ranefList0,
                               y = ylist0,
                               users = subIdx[[ss]]
                               )
    }
    nranef0 <- ncol(dataList[[1]]$z[[1]])
    nfixef0 <- ncol(dataList[[1]]$x[[1]])
} else {
    print("peace")
}

library(parallel)

numCores <- detectCores()

if (numCores >= nsubs) {
    cat("**** working with ", numCores, " cores, but we only need ",  nsubs, " cores ****", "\n")
    cl <- makeCluster(nsubs)
    ## set random seed
    ## RNGkind("L'Ecuyer-CMRG")
    ## clusterSetRNGStream(cl, iseed = 2020)
} else {
    cat("**** PROBLEM: working with ", numCores, " cores, but need ",  nsubs, " cores :PROBLEM  ***", "\n")
}

## Assign the variables appropriately in the workers' GlobalEnv. The 'recvd'
## variable is a consistency check to ensure we have succeeded.
sendData <- clusterApply(cl, dataList,
                         function (vars) {
                             x <<- vars$x
                             z <<- vars$z
                             y <<- vars$y
                             recvd <<- TRUE
                         }
                         )

workerX <- clusterEvalQ(cl, x)

## we see that the below loop outputs TRUE, so we are fine
for (ii in 1:nsubs) {
    cat(all.equal(workerX[[ii]], lapply(dataList, function(ll) ll$x)[[ii]]), "\n")
}

subsetSuffStat <- function (ylist, fixefList, ranefList, fixSamp, errVar, tmatTilde, gammaMat) {
    nsample <- length(fixefList)
    nranef <- ncol(ranefList[[1]])
    nfixef <- ncol(fixefList[[1]]) # p
    ndim <- nranef * nranef

    ranSampList <- vector("list", nsample)
    for (ii in 1:nsample) {
        tmp1 <- ranefList[[ii]] %*% gammaMat %*% tmatTilde # z g t
        wmat <- chol2inv(chol(tcrossprod(ranefList[[ii]] %*% gammaMat, tmp1) + errVar * diag(1, nrow(ranefList[[ii]]))))
        postRanVar <- tmatTilde - crossprod(tmp1, wmat %*% tmp1)
        postRanMean <- drop(crossprod(tmp1, wmat %*% (ylist[[ii]] - fixefList[[ii]] %*% fixSamp)))
        ranSampList[[ii]] <- postRanMean + drop(crossprod(chol(postRanVar), rnorm(length(postRanMean))))
    }
    rans <<- ranSampList ## need to have the random effects in the GlobalEnv

    vmat <- Reduce("+", lapply(rans, function(x) tcrossprod(x, x)))

    xztxz <- matrix(0.0, nfixef + ndim, nfixef + ndim)
    xzty <- matrix(0.0, nfixef + ndim, 1)
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


    list(xztxz = xztxz,
         xzty  = xzty,
         vmat = vmat
         )
}

clusterExport(cl, "subsetSuffStat")

calcRss <- function (ylist, fixefList, ranefList, ranSampList, alphaHat) {
    nsample <- length(fixefList)
    nranef <- ncol(ranefList[[1]]) # q
    ndim <- nranef * nranef

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

    rss
}

clusterExport(cl, "calcRss")

## START: intitialize a few variables before running our algorithm
nfixef <- nfixef0
nranef <- nranef0
ndim <- nranef * nranef
nsample <- length(unique(train$group))
nobsv <- length(train$y)

tmat <- diag(1, nranef)
gmat <- diag(1, nranef)
tmatTilde <- tmat
errVar <- 2
fixs <- rep(0, nfixef)

sampTmat <- list()
sampErrVar <- list()
sampBeta <- list()

nwait <- round(frac * nsubs)
niter <- 10000
track <- rep(0, nsubs)

eps <- 0.1
## END: intitialize a few variables before running our algorithm

## intitalize the sufficient statistics and update the parameter
## estimates accordingly
suffStat <- clusterCall(cl, function(pars) subsetSuffStat(y, x, z, pars$betas, pars$errVar, pars$tmatTilde, pars$gmat),
                        pars = list(betas = fixs, errVar = errVar, tmatTilde = tmatTilde, gmat = gmat))

xztxz <- Reduce("+", lapply(suffStat, function(x) x$xztxz))
xzty <- Reduce("+", lapply(suffStat, function(x) x$xzty))
vmat <- Reduce("+", lapply(suffStat, function(x) x$vmat))

xztxzInv <- chol2inv(chol(xztxz))
alphaHat <- xztxzInv %*% xzty

rss <- clusterCall(cl, function(pars) calcRss(y, x, z, rans, pars), pars = alphaHat)
sumRss <- Reduce("+", rss)

## draw tau, beta, Gamma
errVar <- sumRss / rchisq(1, nobsv - nfixef - ndim)

alphas <- drop(alphaHat + crossprod(chol(errVar * xztxzInv), rnorm(nfixef + ndim)))

fixs <- alphas[1:nfixef]
gmat <- matrix(alphas[-(1:nfixef)], nranef, nranef)
tmatTilde <- chol2inv(chol(rWishart(1, nsample - nranef, solve(vmat))[ , , 1]))
tmat <- tcrossprod(gmat %*% tmatTilde, gmat)
## done with the setup for running the algorithm

## start the sampling algorithm
startTime <- proc.time()
for (its in 1:niter) {
    if (its %% 10 == 0) {
        cat("iter: ", its, "\n")
        cat("errs: ", round(mean(unlist(sampErrVar)), 2), "\n")
        cat("beta: ", round(colMeans(do.call(rbind, sampBeta)), 2), "\n")
        cat("vars: ", round(colMeans(do.call(rbind, lapply(sampTmat, function(x) diag(x)))), 2), "\n")
    }

    mup <- 1 - rbinom(1, 1, eps)

    if (mup == 1) {
        waitIdx <- sort(sample(1:nsubs, nwait))
        ## draw random effects
        recvdSuffStat <- clusterCall(cl, function(pars) subsetSuffStat(y, x, z, pars$betas, pars$errVar, pars$tmatTilde, pars$gmat),
                                     pars = list(betas = fixs, errVar = errVar, tmatTilde = tmatTilde, gmat = gmat))


        for (ww in waitIdx) {
            suffStat[[ww]]$xztxz <- recvdSuffStat[[ww]]$xztxz
            suffStat[[ww]]$xzty <- recvdSuffStat[[ww]]$xzty
            suffStat[[ww]]$vmat <- recvdSuffStat[[ww]]$vmat
        }

        xztxz <- Reduce("+", lapply(suffStat, function(x) x$xztxz))
        xzty <- Reduce("+", lapply(suffStat, function(x) x$xzty))
        vmat <- Reduce("+", lapply(suffStat, function(x) x$vmat))

        xztxzInv <- chol2inv(chol(xztxz))
        alphaHat <- xztxzInv %*% xzty

        recvdRss <- clusterCall(cl, function(pars) calcRss(y, x, z, rans, pars), pars = alphaHat)


        ## update the random effects that have returned the results
        for (ww in waitIdx) {
            rss[[ww]] <- recvdRss[[ww]]
        }

        sumRss <- Reduce("+", rss)

        ## draw tau, beta, L
        errVar <- sumRss / rchisq(1, nobsv - nfixef - ndim)

        alphas <- drop(alphaHat + crossprod(chol(errVar * xztxzInv), rnorm(nfixef + ndim)))

        fixs <- alphas[1:nfixef]
        gmat <- matrix(alphas[-(1:nfixef)], nranef, nranef)
        tmatTilde <- chol2inv(chol(rWishart(1, nsample - nranef, solve(vmat))[ , , 1]))
        tmat <- tcrossprod(gmat %*% tmatTilde, gmat)

        ## updated the books!
        sampErrVar[[its]] <- errVar
        sampBeta[[its]] <- fixs
        sampTmat[[its]] <- tmat
        track[waitIdx] <- 1 + track[waitIdx]

    } else {
        suffStat <- clusterCall(cl, function(pars) subsetSuffStat(y, x, z, pars$betas, pars$errVar, pars$tmatTilde, pars$gmat),
                                pars = list(betas = fixs, errVar = errVar, tmatTilde = tmatTilde, gmat = gmat))

        xztxz <- Reduce("+", lapply(suffStat, function(x) x$xztxz))
        xzty <- Reduce("+", lapply(suffStat, function(x) x$xzty))
        vmat <- Reduce("+", lapply(suffStat, function(x) x$vmat))

        xztxzInv <- chol2inv(chol(xztxz))
        alphaHat <- xztxzInv %*% xzty

        rss <- clusterCall(cl, function(pars) calcRss(y, x, z, rans, pars), pars = alphaHat)
        sumRss <- Reduce("+", rss)

        ## draw tau, beta, Gamma
        errVar <- sumRss / rchisq(1, nobsv - nfixef - ndim)

        alphas <- drop(alphaHat + crossprod(chol(errVar * xztxzInv), rnorm(nfixef + ndim)))

        fixs <- alphas[1:nfixef]
        gmat <- matrix(alphas[-(1:nfixef)], nranef, nranef)
        tmatTilde <- chol2inv(chol(rWishart(1, nsample - nranef, solve(vmat))[ , , 1]))
        tmat <- tcrossprod(gmat %*% tmatTilde, gmat)


        ## updated the books!
        sampErrVar[[its]] <- errVar
        sampBeta[[its]] <- fixs
        sampTmat[[its]] <- tmat
        track <- 1 + track
    }
}
endTime <- proc.time()
res <- list('sigma' = unlist(sampErrVar),
            'tmat' = do.call(rbind, lapply(sampTmat, function(x) c(diag(x), x[lower.tri(x)]))),
            'betas' = do.call(rbind, sampBeta),
            'track' = track,
            'time' = endTime[3] - startTime[3]
            )

fname <- paste0("/Shared/ssrivastva/dist_mov_len/result/sub/ml_part_n_", nobsid, "_k_", nsubid, "_frac_", fid, "_rep_", cid, "_eps_0.1.rds")
saveRDS(res, fname)


cat("Closing workers \n")

stopCluster(cl)
