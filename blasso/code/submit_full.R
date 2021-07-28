cmdArgs <- commandArgs(trailingOnly = TRUE)

mtd <- as.numeric(cmdArgs[1])
id <- as.numeric(cmdArgs[2])

if (mtd == 1) {
    source("full_sampler.R")

    cvtrain <- readRDS("/Shared/ssrivastva/dist_blasso/data/data_n_1.rds")
    train <- cvtrain[[id]]
    rm(cvtrain)

    y0 <- as.numeric(train$y)
    x0 <- as.matrix(train$x)

    res <- lassoSampler(y0, x0, a0 = 0, b0 = 0, niter = 20000)
    fname <- paste0("/Shared/ssrivastva/dist_blasso/result/full/blasso_full_n_1_rep_", id, ".rds")
    saveRDS(res, fname)
} else if (mtd == 2) {
    source("full_sampler.R")

    cvtrain <- readRDS("/Shared/ssrivastva/dist_blasso/data/data_n_2.rds")
    train <- cvtrain[[id]]
    rm(cvtrain)

    y0 <- as.numeric(train$y)
    x0 <- as.matrix(train$x)

    res <- lassoSampler(y0, x0, a0 = 0, b0 = 0, niter = 20000)
    fname <- paste0("/Shared/ssrivastva/dist_blasso/result/full/blasso_full_n_2_rep_", id, ".rds")
    saveRDS(res, fname)
} else if (mtd == 3) {
    source("full_sampler.R")

    cvtrain <- readRDS("/Shared/ssrivastva/dist_blasso/data/data_n_3.rds")
    train <- cvtrain[[id]]
    rm(cvtrain)

    y0 <- as.numeric(train$y)
    x0 <- as.matrix(train$x)

    res <- lassoSampler(y0, x0, a0 = 0, b0 = 0, niter = 20000)
    fname <- paste0("/Shared/ssrivastva/dist_blasso/result/full/blasso_full_n_3_rep_", id, ".rds")
    saveRDS(res, fname)
} else {
    print("peace")
}
