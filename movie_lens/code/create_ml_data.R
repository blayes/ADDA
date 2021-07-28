rm(list=ls())
setwd("~/dist_mov_len/code")
mlData <- readRDS("../data/ml_full.rds")

set.seed(12345)
group <- as.numeric(mlData$group)
grpLbl <- sort(unique(group))
ngroup <- length(grpLbl)

usrs <- split(1:nrow(mlData$x), group)
movRated <- sapply(usrs, length)
selUsrs <- usrs[which(movRated >= 200)]
selIdx20 <- lapply(selUsrs, function(x) sort(sample(x, 20)))

trainData <- list()
for (cc in 1:10) {
    cat("cc: ", cc, "\n")
    usrIdx <- sort(sample(1:length(selUsrs), 5000))
    ridx <- unlist(selIdx20[usrIdx])
    trainData[[cc]] <- list(x = mlData$x[ridx, ],
                            z = mlData$z[ridx, ],
                            y = mlData$y[ridx],
                            users = usrIdx,
                            group = group[ridx]
                            )
}

saveRDS(trainData, "/Shared/ssrivastva/dist_mov_len/data/train_ml_100k.rds")

rm(list=ls())
setwd("~/dist_mov_len/code")
mlData <- readRDS("../data/ml_full.rds")

set.seed(12345)
group <- as.numeric(mlData$group)
grpLbl <- sort(unique(group))
ngroup <- length(grpLbl)

usrs <- split(1:nrow(mlData$x), group)
movRated <- sapply(usrs, length)
selUsrs <- usrs[which(movRated >= 200)]
selIdx200 <- lapply(selUsrs, function(x) sort(sample(x, 200)))

trainData <- list()
for (cc in 1:10) {
    cat("cc: ", cc, "\n")
    usrIdx <- sort(sample(1:length(selUsrs), 5000))
    ridx <- unlist(selIdx200[usrIdx])
    trainData[[cc]] <- list(x = mlData$x[ridx, ],
                            z = mlData$z[ridx, ],
                            y = mlData$y[ridx],
                            users = usrIdx,
                            group = group[ridx]
                            )
}

saveRDS(trainData, "/Shared/ssrivastva/dist_mov_len/data/train_ml_1000k.rds")

rm(list = ls())

trainData <- readRDS("/Shared/ssrivastva/dist_mov_len/data/train_ml_100k.rds")

set.seed(12345)

parts <- list()
for (cc in 1:10) {
    users <- unique(trainData[[cc]]$group)
    idx10 <- sample(1:10, length(users), replace = TRUE)
    idx25 <- sample(1:25, length(users), replace = TRUE)
    idx50 <- sample(1:50, length(users), replace = TRUE)
    parts[[cc]] <- list("k10" = split(users, idx10),
                        "k25" = split(users, idx25),
                        "k50" = split(users, idx50))
}

saveRDS(parts, "/Shared/ssrivastva/dist_mov_len/data/part_ids_100k.rds")

rm(list = ls())

trainData <- readRDS("/Shared/ssrivastva/dist_mov_len/data/train_ml_1000k.rds")

set.seed(12345)

parts <- list()
for (cc in 1:10) {
    users <- unique(trainData[[cc]]$group)
    idx10 <- sample(1:10, length(users), replace = TRUE)
    idx25 <- sample(1:25, length(users), replace = TRUE)
    idx50 <- sample(1:50, length(users), replace = TRUE)
    parts[[cc]] <- list("k10" = split(users, idx10),
                        "k25" = split(users, idx25),
                        "k50" = split(users, idx50))
}

saveRDS(parts, "/Shared/ssrivastva/dist_mov_len/data/part_ids_1000k.rds")



######### LOGISTIC REGRESSION ########

rm(list=ls())
setwd("~/dist_mov_len/code")
mlData <- readRDS("../data/ml_full.rds")


trainData <- list()
for (cc in 1:10) {
    cat("cc: ", cc, "\n")
    ridx <- sample(1:nrow(mlData$x), 1e7)
    trainData[[cc]] <- list(x = mlData$x[ridx, ],
                            y = mlData$y[ridx]
                            )
}

saveRDS(trainData, "/Shared/ssrivastva/dist_mov_len/data/train_ml_1e7.rds")


trainData <- list()
for (cc in 1:10) {
    cat("cc: ", cc, "\n")
    ridx <- sample(1:nrow(mlData$x), 1e6)
    trainData[[cc]] <- list(x = mlData$x[ridx, ],
                            y = mlData$y[ridx]
                            )
}

saveRDS(trainData, "/Shared/ssrivastva/dist_mov_len/data/train_ml_1e6.rds")


### partition the full data into subsets ###
rm(list = ls())

part1 <- list("k10" = rep(1:10, each = 1e5),
              "k25" = rep(1:25, each = 40000),
              "k50" = rep(1:50, each = 20000))
part2 <- list("k10" = rep(1:10, each = 1e6),
              "k25" = rep(1:25, each = 400000),
              "k50" = rep(1:50, each = 200000))

saveRDS(list(part1, part2), "/Shared/ssrivastva/dist_mov_len/data/part_ml.rds")


