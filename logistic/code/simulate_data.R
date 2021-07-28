rm(list=ls())

setwd("~/dist_logistic/code")

library(matrixStats)
library(Matrix)

genData <- function (nobs, ndim, ntrail) {
    ## fixed effects coefficients
    fixef <- rep(c(-2, 2), length = ndim)

    ## generate feature  matrices with Pr(x[i,j] = +1) = P(x[i,j] = -1) = 1/2,
    x <- matrix(sample(c(-1, +1), nobs * ndim, replace=TRUE), nobs, ndim)

    ## compute linear predictors and generate observations
    mu <- drop(x %*% fixef)
    probs <- 1 / (1 + exp(- mu))
    y <- rbinom(n = nobs, size = ntrail, prob = probs)

    list(x = x, probs = probs, y = y)
}

nobs <- 1e4
ndim <- 10
ntrail <- 10

set.seed(12345)

repData <- list()
for (cc in 1:10) {
    repData[[cc]] <- genData(nobs, ndim, ntrail)
}
    
saveRDS(repData, "/Shared/ssrivastva/dist_logistic/data/data_n_1.rds")

rm(list = ls())

nobs <- 1e5
ndim <- 10
ntrail <- 10

set.seed(12345)

repData <- list()
for (cc in 1:10) {
    repData[[cc]] <- genData(nobs, ndim, ntrail)
}
    
saveRDS(repData, "/Shared/ssrivastva/dist_logistic/data/data_n_2.rds")


rm(list = ls())

nobs <- 1e6
ndim <- 10
ntrail <- 10

set.seed(12345)

repData <- list()
for (cc in 1:10) {
    repData[[cc]] <- genData(nobs, ndim, ntrail)
}
    
saveRDS(repData, "/Shared/ssrivastva/dist_logistic/data/data_n_3.rds")



### partition the full data into subsets ###
rm(list = ls())

set.seed(12345)

part1 <- part2 <- part3 <- list()
for (cc in 1:10) {
    part1[[cc]] <- list("k10" = sample(1:10, 1e4, replace = TRUE),
                        "k25" = sample(1:25, 1e4, replace = TRUE),
                        "k50" = sample(1:50, 1e4, replace = TRUE),
                        "k100" = sample(1:100, 1e4, replace = TRUE))
    
    part2[[cc]] <- list("k10" = sample(1:10, 1e5, replace = TRUE),
                        "k25" = sample(1:25, 1e5, replace = TRUE),
                        "k50" = sample(1:50, 1e5, replace = TRUE),
                        "k100" = sample(1:100, 1e5, replace = TRUE))
    
    part3[[cc]] <- list("k10" = sample(1:10, 1e6, replace = TRUE),
                        "k25" = sample(1:25, 1e6, replace = TRUE),
                        "k50" = sample(1:50, 1e6, replace = TRUE),
                        "k100" = sample(1:100, 1e6, replace = TRUE))        
}

saveRDS(part1, "/Shared/ssrivastva/dist_logistic/data/part_n_1.rds")
saveRDS(part2, "/Shared/ssrivastva/dist_logistic/data/part_n_2.rds")
saveRDS(part3, "/Shared/ssrivastva/dist_logistic/data/part_n_3.rds")


