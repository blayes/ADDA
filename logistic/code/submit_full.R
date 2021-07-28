cmdArgs <- commandArgs(trailingOnly = TRUE)

mtd <- as.numeric(cmdArgs[1])
id <- as.numeric(cmdArgs[2])

if (mtd == 1) {
    source("full_sampler.R")
    
    cvtrain <- readRDS("/Shared/ssrivastva/dist_logistic/data/data_n_1.rds")
    train <- cvtrain[[id]]
    rm(cvtrain)
    
    y0 <- as.numeric(train$y)
    x0 <- as.matrix(train$x)
    ntrail0 <- as.numeric(rep(10, length(y0)))

    
    res <- polyaGammaSampler(yvec = y0, xmat = x0, ntrail = ntrail0, mub0 = rep(0, ncol(x0)), sigbInv0 = diag(0.01, ncol(x0)), niter = 20000) 

    fname <- paste0("/Shared/ssrivastva/dist_logistic/result/sim/full/sim_full_n_1_rep_", id, ".rds")
    saveRDS(res, fname)
} else if (mtd == 2) {
    source("full_sampler.R")
    
    cvtrain <- readRDS("/Shared/ssrivastva/dist_logistic/data/data_n_2.rds")
    train <- cvtrain[[id]]
    rm(cvtrain)
    
    y0 <- as.numeric(train$y)
    x0 <- as.matrix(train$x)
    ntrail0 <- as.numeric(rep(10, length(y0)))

    res <- polyaGammaSampler(yvec = y0, xmat = x0, ntrail = ntrail0, mub0 = rep(0, ncol(x0)), sigbInv0 = diag(0.01, ncol(x0)), niter = 20000) 

    fname <- paste0("/Shared/ssrivastva/dist_logistic/result/sim/full/sim_full_n_2_rep_", id, ".rds")
    saveRDS(res, fname)
} else if (mtd == 3) {
    source("full_sampler.R")
    
    cvtrain <- readRDS("/Shared/ssrivastva/dist_logistic/data/data_n_3.rds")
    train <- cvtrain[[id]]
    rm(cvtrain)
    
    y0 <- as.numeric(train$y)
    x0 <- as.matrix(train$x)
    ntrail0 <- as.numeric(rep(10, length(y0)))

    res <- polyaGammaSampler(yvec = y0, xmat = x0, ntrail = ntrail0, mub0 = rep(0, ncol(x0)), sigbInv0 = diag(0.01, ncol(x0)), niter = 20000) 

    fname <- paste0("/Shared/ssrivastva/dist_logistic/result/sim/full/sim_full_n_3_rep_", id, ".rds")
    saveRDS(res, fname)
} else {
    print("peace")
}
