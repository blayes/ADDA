cmdArgs <- commandArgs(trailingOnly = TRUE)

mtd <- as.numeric(cmdArgs[1])
id <- as.numeric(cmdArgs[2])

if (mtd == 1) {
    source("full_sampler.R")
    
    cvtrain <- readRDS("/Shared/ssrivastva/dist_lme/data/data_n_1.rds")
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
                       niter = 20000) 

    fname <- paste0("/Shared/ssrivastva/dist_lme/result/sim/full/sim_full_n_1_rep_", id, ".rds")
    saveRDS(res, fname)
} else if (mtd == 2) {
    source("full_sampler.R")
    
    cvtrain <- readRDS("/Shared/ssrivastva/dist_lme/data/data_n_2.rds")
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
                       niter = 20000) 

    fname <- paste0("/Shared/ssrivastva/dist_lme/result/sim/full/sim_full_n_2_rep_", id, ".rds")
    saveRDS(res, fname)
} else if (mtd == 3) {
    source("full_sampler.R")
    
    cvtrain <- readRDS("/Shared/ssrivastva/dist_lme/data/data_n_3.rds")
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
                       niter = 20000) 
  
    fname <- paste0("/Shared/ssrivastva/dist_lme/result/sim/full/sim_full_n_3_rep_", id, ".rds")
    saveRDS(res, fname)
} else {
    print("peace")
}
