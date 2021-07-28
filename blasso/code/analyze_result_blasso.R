setwd("/Shared/ssrivastva/dist_blasso/result/")
rm(list = ls())

library(KernSmooth)
library(plyr)
library(ggplot2)

nfrac <- c(20, 40, 60, 80)

###### ACCURACY COMPUTATIONS
acc <- array(NA, dim = c(3, 3, 4, 10, 20),
             dimnames = list(paste0("n", 1:3),
                             paste0("k", 1:3),
                             paste0("f", 1:4),
                             paste0("cv", 1:10),
                             paste0("iter", 1:20)
                             ))

times <- array(NA, dim = c(3, 3, 4, 10, 2),
               dimnames = list(paste0("n", 1:3),
                               paste0("k", 1:3),
                               paste0("f", 1:4),
                               paste0("cv", 1:10),
                               c("full", "part")
                               ))

dims <- c(50, 500, 5000)
for (nn in 1:3) {
    for (kk in 1:3) {
        for (ff in 1:4) {
            for (cc in 1:10) {
                cat("kk: ", kk, " ff: ", ff, " cc: ", cc, "\n")
                fname <- paste0("part/part_n_", nn, "_k_", kk, "_frac_", nfrac[ff], "_rep_", cc, ".rds")
                if (file.exists(fname)) {
                    res <- readRDS(fname)
                    ffname <- paste0("full/blasso_full_n_", nn, "_rep_", cc, ".rds")
                    fres <- readRDS(ffname)
                    pbetas <- res$betas
                    fbetas <- fres$betas
                    times[nn, kk, ff, cc, "full"] <- fres$time
                    times[nn, kk, ff, cc, "part"] <- res$time

                    margAcc <- matrix(NA, dims[nn], 20)
                    for (dd in 1:dims[nn]) {
                        rr <- range(c(pbetas[100:20000, dd], fbetas[100:20000, dd]))
                        for (its in 1:20) {
                            pdens <- bkde(x         = pbetas[100:(its * 1000), dd],
                                          bandwidth = dpik(pbetas[1:(its * 1000), dd], gridsize = 1000, range.x = rr),
                                          gridsize = 1000,
                                          range.x   = rr
                                          )
                            fdens <- bkde(x         = fbetas[100:(its * 1000), dd],
                                          bandwidth = dpik(fbetas[1:(its * 1000), dd], gridsize = 1000, range.x = rr),
                                          gridsize = 1000,
                                          range.x   = rr
                                          )
                            margAcc[dd, its] <- max(1 - sum(abs(pdens$y - fdens$y) * diff(pdens$x)[1]) / 2, 0)
                        }
                        acc[nn, kk, ff, cc, ] <- colMeans(margAcc, na.rm = TRUE)
                    }
                }
            }
        }
    }
}

saveRDS(list(acc = acc, times = times), "res_blasso.rds")

dims <- c(50, 500, 5000)
for (nn in 1:3) {
    for (kk in 1:3) {
        for (ff in 1:4) {
            for (cc in 1:10) {
                cat("kk: ", kk, " ff: ", ff, " cc: ", cc, "\n")
                fname <- paste0("part/part_n_", nn, "_k_", kk, "_frac_", nfrac[ff], "_rep_", cc, ".rds")
                if (file.exists(fname)) {
                    res <- readRDS(fname)
                    ffname <- paste0("full/blasso_full_n_", nn, "_rep_", cc, ".rds")
                    fres <- readRDS(ffname)
                    pbetas <- res$betas
                    fbetas <- fres$betas
                    times[nn, kk, ff, cc, "full"] <- fres$time
                    times[nn, kk, ff, cc, "part"] <- res$time

                    margAcc <- matrix(NA, 50, 20)
                    di <- 0
                    for (dd in sample(1:dims[nn], 50)) {
                        di <- di + 1
                        rr <- range(c(pbetas[1:20000, dd], fbetas[1:20000, dd]))
                        for (its in 1:20) {
                            pdens <- bkde(x         = pbetas[1:(its * 1000), dd],
                                          bandwidth = dpik(pbetas[1:(its * 1000), dd], gridsize = 1000, range.x = rr),
                                          gridsize = 1000,
                                          range.x   = rr
                                          )
                            fdens <- bkde(x         = fbetas[1:(its * 1000), dd],
                                          bandwidth = dpik(fbetas[1:(its * 1000), dd], gridsize = 1000, range.x = rr),
                                          gridsize = 1000,
                                          range.x   = rr
                                          )
                            margAcc[di, its] <- max(1 - sum(abs(pdens$y - fdens$y) * diff(pdens$x)[1]) / 2, 0)
                        }
                        acc[nn, kk, ff, cc, ] <- colMeans(margAcc, na.rm = TRUE)
                    }
                }
            }
        }
    }
}

accBetaEst <- apply(acc, c(1:3, 5), mean, na.rm = TRUE)

accBetaEstDf <- adply(accBetaEst, 1:4)
colnames(accBetaEstDf) <- c("nsamp", "nsub", "frac", "niter", "est")

accBetaEstDf$niter <- as.numeric(gsub("iter", "", as.character(accBetaEstDf$niter)))

plt <- ggplot(accBetaEstDf, aes(x = niter, y = est, group = frac, colour = frac)) +
    xlab("Iterations (x 1000)") + ylab("Accuracy of " ~ beta) +
    scale_colour_brewer(name="",
                        breaks=c("f1", "f2", "f3", "f4"),
                        labels=c("r = 0.20", "r = 0.40", "r = 0.60", "r = 0.80"),
                        palette="Set1") +
    facet_grid(nsub  ~ nsamp,
               labeller = labeller(.rows = c(k1 = "k = 10", k2 = "k = 25", k3 = "k = 50"),
                                   .cols = c(n1 = "n = 50", n2 = "n = 500", n3 = "n = 5000")
                                   )
               ) +
    geom_line(stat = "identity", size = 2) + geom_point(size = 4) +
    theme_bw() +
    theme(strip.text.x = element_text(size=24, face="bold"),
          strip.text.y = element_text(size=24, face="bold"),
          axis.text.x = element_text(size=20, face="bold"),
          axis.text.y = element_text(size=20, face="bold"),
          axis.title.x = element_text(size=24, face="bold"),
          axis.title.y = element_text(size=24, face="bold"),
          legend.text = element_text(size = 24)
          ) +
    ylim(0.70, 1.00)

pdf("~/dist_blasso/blasso_beta.pdf", 12, 9)
print(plt)
dev.off()

####### MCSE COMPUTATIONS
setwd("/Shared/ssrivastva/dist_blasso/result/")
rm(list=ls())

library(mcmcse)

nfrac <- c(20, 40, 60, 80)

mcmcse <- array(NA, dim = c(3, 3, 4, 10, 20),
             dimnames = list(paste0("n", 1:3),
                             paste0("k", 1:3),
                             paste0("f", 1:4),
                             paste0("cv", 1:10),
                             paste0("iter", 1:20)
                             ))

dims <- c(50, 500, 5000)
for (nn in 1:3) {
    for (kk in 1:3) {
        for (ff in 1:4) {
            for (cc in 1:10) {
                cat("kk: ", kk, " ff: ", ff, " cc: ", cc, "\n")
                fname <- paste0("part/part_n_", nn, "_k_", kk, "_frac_", nfrac[ff], "_rep_", cc, ".rds")
                if (file.exists(fname)) {
                    res <- readRDS(fname)
                    ffname <- paste0("full/blasso_full_n_", nn, "_rep_", cc, ".rds")
                    fres <- readRDS(ffname)
                    pbetas <- res$betas
                    fbetas <- fres$betas

                    for (its in 1:20) {
                        mcmcse[nn, kk, ff, cc, its] <- mean(mcse.mat(pbetas[100:(its * 1000), ], method = "obm")[ , 2] - mcse.mat(fbetas[100:(its * 1000), ], method = "obm")[ , 2])
                    }
                }
            }
        }
    }
}

saveRDS(mcmcse, "~/dist_blasso/data/mcmcse_betas.rds")

betaEst <- apply(abs(mcmcse), c(1:3, 5), max, na.rm = TRUE)

library(plyr)
library(ggplot2)

betaEstDf <- adply(betaEst, 1:4)
colnames(betaEstDf) <- c("nsamp", "nsub", "frac", "niter", "est")

betaEstDf$niter <- as.numeric(gsub("iter", "", as.character(betaEstDf$niter)))

plt <- ggplot(betaEstDf, aes(x = niter, y = est, group = frac, colour = frac)) +
    xlab("Iterations (x 1000)") + ylab("SE Difference for" ~ beta) +
    scale_colour_brewer(name="",
                        breaks=c("f1", "f2", "f3", "f4"),
                        labels=c("r = 0.20", "r = 0.40", "r = 0.60", "r = 0.80"),
                        palette="Set1",
                        guide = guide_legend(override.aes = list(shape = NA))) +
    facet_grid(nsub  ~ nsamp,
               labeller = labeller(.rows = c(k1 = "k = 10", k2 = "k = 25", k3 = "k = 50"),
                                   .cols = c(n1 = "p = 50", n2 = "p = 500", n3 = "p = 5000")
                                   )
               ) +
    geom_line(stat = "identity", size = 4) + geom_point(size = 8) +
    theme_bw() +
    theme(strip.text.x = element_text(size=24, face=
"bold"),
          strip.text.y = element_text(size=24, face="bold"),
          axis.text.x = element_text(size=20, face="bold"),
          axis.text.y = element_text(size=20, face="bold"),
          axis.title.x = element_text(size=24, face="bold"),
          axis.title.y = element_text(size=24, face="bold"),
          legend.text = element_text(size = 24)
          ) + ylim(0, 0.002)

pdf("~/dist_blasso/mcse_blasso_beta.pdf", 12, 9)
print(plt)
dev.off()
