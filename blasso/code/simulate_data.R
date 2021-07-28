rm(list=ls())

setwd("~/dist_blasso/code")

library(matrixStats)
library(Matrix)

genData <- function (nsamp, ndim) {
    fsparse <- 0.90
    nzero <- round(fsparse * ndim)
    sig <- 0.01

    ## beta
    beta0 <- c(rep(0, nzero), rep(c(-2, 2), length = ndim - nzero))
    ## x
    xmat <- matrix(rnorm(nsamp * ndim), nrow = nsamp, ncol = ndim)
    ## y
    yvec <- drop(xmat %*% beta0) + rnorm(nsamp, sd = sqrt(sig))

    list(x = xmat, beta = beta0, y = yvec)
}

nobs <- 50
ndim <- 50

set.seed(12345)

repData <- list()
for (cc in 1:10) {
    repData[[cc]] <- genData(nobs, ndim)
}

saveRDS(repData, "/Shared/ssrivastva/dist_blasso/data/data_n_1.rds")

rm(list = ls())

nobs <- 500
ndim <- 500

set.seed(12345)

repData <- list()
for (cc in 1:10) {
    repData[[cc]] <- genData(nobs, ndim)
}

saveRDS(repData, "/Shared/ssrivastva/dist_blasso/data/data_n_2.rds")

rm(list = ls())

nobs <- 5000
ndim <- 5000

set.seed(12345)

repData <- list()
for (cc in 1:10) {
    repData[[cc]] <- genData(nobs, ndim)
}

saveRDS(repData, "/Shared/ssrivastva/dist_blasso/data/data_n_3.rds")

### partition the full data into subsets ###
rm(list = ls())

set.seed(12345)

part1 <- part2 <- part3 <- list()
for (cc in 1:10) {
    part1[[cc]] <- list("k10" = rep(1:10, each = 5),
                        "k25" = rep(1:25, each = 2),
                        "k50" = rep(1:50, each = 1))

    part2[[cc]] <- list("k10" = rep(1:10, each = 50),
                        "k25" = rep(1:25, each = 20),
                        "k50" = rep(1:50, each = 10))

    part3[[cc]] <- list("k10" = rep(1:10, each = 500),
                        "k25" = rep(1:25, each = 200),
                        "k50" = rep(1:50, each = 100))
}

saveRDS(part1, "/Shared/ssrivastva/dist_blasso/data/part_n_1.rds")
saveRDS(part2, "/Shared/ssrivastva/dist_blasso/data/part_n_2.rds")
saveRDS(part3, "/Shared/ssrivastva/dist_blasso/data/part_n_3.rds")
