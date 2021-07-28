setwd("~/dist_lme/code")

rm(list=ls())

library(matrixStats)
library(Matrix)

## See mbest package at http://ptrckprry.com/code/ Perry (2017) in JRSS-B.
genData <- function (ngroup, nobs, nfixef, nranef) {
    library(matrixStats)
    library(Matrix)

    ## fixed effects coefficients
    fixef <- rep(c(-2, 2), length = nfixef)
    if (nranef == 3) {
        ranefCorr <- matrix(c(1, -0.4, 0.3,
                              -0.4, 1, 0.001,
                              0.3, 0.001, 1),
                            nranef, nranef)
    } else {
        ranefCorr <- as.matrix(bdiag(rep(list(matrix(c(1, -0.4, 0.3,
                                                       -0.4, 1, 0.001,
                                                       0.3, 0.001, 1),
                                                     3, 3)), 2)))
    }
    ranefCov <- outer(sqrt(1:nranef), sqrt(1:nranef)) * ranefCorr
    ranefCovSqrt <- chol(ranefCov)

    ## generate coefficients
    u <- matrix(rnorm(ngroup * nranef), ngroup, nranef)
    ranef <- u %*% ranefCovSqrt

    ## generate group
    group <- rep(1:ngroup, each = nobs / ngroup)

    ## generate feature  matrices with Pr(x[i,j] = +1) = P(x[i,j] = -1) = 1/2,
    x <- matrix(sample(c(-1, +1), nobs * nfixef, replace=TRUE), nobs, nfixef)
    z <- matrix(sample(c(-1, +1), nobs * nranef, replace=TRUE), nobs, nranef)

    ## compute linear predictors and generate observations
    mu <- drop(x %*% fixef) + rowSums(z * ranef[group, ])
    y <- rnorm(nobs, mean=mu, sd=1)

    list(ngroup = ngroup, nobs = nobs,
         fixef = fixef,
         ranefCov = ranefCov,
         group = group,
         ranef = ranef,
         x = x,
         z = z,
         y = y
         )
}

nobs <- 1e4
ngroup <- 1000
nfixef <- 4
nranef <- 3

set.seed(12345)

repData <- list()
for (cc in 1:10) {
    cat("cc: ", cc, "\n")
    repData[[cc]] <- genData(ngroup, nobs, nfixef, nranef)
}
    
saveRDS(repData, "/Shared/ssrivastva/dist_lme/data/data_n_1.rds")

rm(list = ls())

nobs <- 1e5
ngroup <- 10000
nfixef <- 4
nranef <- 3

set.seed(12345)

repData <- list()
for (cc in 1:10) {
    cat("cc: ", cc, "\n")
    repData[[cc]] <- genData(ngroup, nobs, nfixef, nranef)
}

    
saveRDS(repData, "/Shared/ssrivastva/dist_lme/data/data_n_2.rds")

rm(list = ls())

nobs <- 1e6
ngroup <- 100000
nfixef <- 4
nranef <- 3

set.seed(12345)

repData <- list()
for (cc in 1:10) {
    cat("cc: ", cc, "\n")
    repData[[cc]] <- genData(ngroup, nobs, nfixef, nranef)
}

saveRDS(repData, "/Shared/ssrivastva/dist_lme/data/data_n_3.rds")


### partition the full data into subsets ###
rm(list = ls())

set.seed(12345)

part1 <- part2 <- part3 <- list()
for (cc in 1:10) {
    part1[[cc]] <- list("k10" = rep(1:10, each = 100),
                        "k25" = rep(1:25, each = 40),
                        "k50" = rep(1:50, each = 20),
                        "k100" = rep(1:100, each = 10))
    
    part2[[cc]] <- list("k10" = rep(1:10, each = 1000),
                        "k25" = rep(1:25, each = 400),
                        "k50" = rep(1:50, each = 200),
                        "k100" = rep(1:100, each = 100))
    
    part3[[cc]] <- list("k10" = rep(1:10, each = 10000),
                        "k25" = rep(1:25, each = 4000),
                        "k50" = rep(1:50, each = 2000),
                        "k100" = rep(1:100, each = 1000))
}

saveRDS(part1, "/Shared/ssrivastva/dist_lme/data/part_n_1.rds")
saveRDS(part2, "/Shared/ssrivastva/dist_lme/data/part_n_2.rds")
saveRDS(part3, "/Shared/ssrivastva/dist_lme/data/part_n_3.rds")

