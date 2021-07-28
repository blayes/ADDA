cmdArgs <- commandArgs(trailingOnly = TRUE)

mtd <- as.numeric(cmdArgs[1])
id <- as.numeric(cmdArgs[2])

if (mtd == 1) {
    source("full_sampler_lme.R")

    cvtrain <- readRDS("/Shared/ssrivastva/dist_mov_len/data/train_ml_100k.rds")
    train <- cvtrain[[id]]
    rm(cvtrain)

    group <- train$group
    grpLbl <- sort(unique(group))
    ngroup <- length(grpLbl)
    ranefList0 <- list()
    fixefList0 <- list()
    ylist0 <- list()
    grpIdx0 <- list()
    for (gg in 1:ngroup) {
        grpIdx0[[gg]] <- which(group == grpLbl[gg])
        fixefList0[[gg]] <- train$x[grpIdx0[[gg]], , drop = FALSE]
        ranefList0[[gg]] <- train$z[grpIdx0[[gg]], , drop = FALSE]
        ylist0[[gg]] <- train$y[grpIdx0[[gg]]]
    }

    nranef0 <- ncol(ranefList0[[1]])
    nfixef0 <- ncol(fixefList0[[1]])

    res <- lmeSampler (ylist0, fixefList0, ranefList0,
                       rep(0, nfixef0), diag(1, nranef0), 2,
                       niter = 10000)

    fname <- paste0("/Shared/ssrivastva/dist_mov_len/result/full/ml_rep_", id, "_100k.rds")
    saveRDS(res, fname)
} else if (mtd == 2) {
    source("full_sampler_lme.R")

    cvtrain <- readRDS("/Shared/ssrivastva/dist_mov_len/data/train_ml_1000k.rds")
    train <- cvtrain[[id]]
    rm(cvtrain)

    group <- train$group
    grpLbl <- sort(unique(group))
    ngroup <- length(grpLbl)
    ranefList0 <- list()
    fixefList0 <- list()
    ylist0 <- list()
    grpIdx0 <- list()
    for (gg in 1:ngroup) {
        grpIdx0[[gg]] <- which(group == grpLbl[gg])
        fixefList0[[gg]] <- train$x[grpIdx0[[gg]], , drop = FALSE]
        ranefList0[[gg]] <- train$z[grpIdx0[[gg]], , drop = FALSE]
        ylist0[[gg]] <- train$y[grpIdx0[[gg]]]
    }

    nranef0 <- ncol(ranefList0[[1]])
    nfixef0 <- ncol(fixefList0[[1]])

    res <- lmeSampler (ylist0, fixefList0, ranefList0,
                       rep(0, nfixef0), diag(1, nranef0), 2,
                       niter = 10000)

    fname <- paste0("/Shared/ssrivastva/dist_mov_len/result/full/ml_rep_", id, "_1000k.rds")
    saveRDS(res, fname)
} else if (mtd == 3) {
    source("full_sampler_logistic.R")

    cvtrain <- readRDS("/Shared/ssrivastva/dist_mov_len/data/train_ml_1e6.rds")
    train <- cvtrain[[id]]
    rm(cvtrain)

    y0 <- as.numeric(train$y > 3)
    x0 <- as.matrix(train$x)
    ntrail0 <- as.numeric(rep(1, length(y0)))

    res <- polyaGammaSampler(yvec = y0, xmat = x0, ntrail = ntrail0, mub0 = rep(0, ncol(x0)), sigbInv0 = diag(0.01, ncol(x0)), niter = 10000) 

    fname <- paste0("/Shared/ssrivastva/dist_mov_len/result/full/ml_logistic_rep_", id, "_1e6.rds")
    saveRDS(res, fname)
} else if (mtd == 4) {
    source("full_sampler_logistic.R")

    cvtrain <- readRDS("/Shared/ssrivastva/dist_mov_len/data/train_ml_1e7.rds")
    train <- cvtrain[[id]]
    rm(cvtrain)

    y0 <- as.numeric(train$y > 3)
    x0 <- as.matrix(train$x)
    ntrail0 <- as.numeric(rep(1, length(y0)))

    res <- polyaGammaSampler(yvec = y0, xmat = x0, ntrail = ntrail0, mub0 = rep(0, ncol(x0)), sigbInv0 = diag(0.01, ncol(x0)), niter = 10000) 

    fname <- paste0("/Shared/ssrivastva/dist_mov_len/result/full/ml_logistic_rep_", id, "_1e7.rds")
    saveRDS(res, fname)
} else {
    print("peace")
}
