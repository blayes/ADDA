setwd("/Shared/ssrivastva/dist_mov_len/result/")
rm(list = ls())

library(KernSmooth)
library(mcmcse)
nfrac <- c(20, 40, 60, 80)

## ACCURACY and MCSE COMPUTATIONS
fullRes <- vector("list", 10)
for (cc in 1:10) {
    fullRes[[cc]] <- vector("list", 2)
    for (nn in 1:2) {
        if (nn == 1) {
            fname <- paste0("full/ml_logistic_rep_", cc, "_1e6.rds")
        } else {
            fname <- paste0("full/ml_logistic_rep_", cc, "_1e7.rds")
        }

        if (file.exists(fname)) {
            fres <- readRDS(fname)
            fbetas <- fres$betas
            fprob <- 1 / (1 + exp(-(drop(fbetas %*% c(rep(0, 4), 1, 0)))))

            fullRes[[cc]][[nn]] <- list(
                beta = fbetas,
                prob = fprob,
                time = fres$time)
        } else {
            fullRes[[cc]][[nn]] <- NA
        }
    }
}

##saveRDS(fullRes, "full/full_summ_logistic.rds")
fullRes <- readRDS("full/full_summ_logistic.rds")

accBeta <- array(NA, dim = c(2, 3, 4, 10, 6, 10, 3),
                 dimnames = list(paste0("n", 1:2),
                                 paste0("k", 1:3),
                                 paste0("f", 1:4),
                                 paste0("cv", 1:10),
                                 paste0("dim", 1:6),
                                 paste0("iter", 1:10),
                                 paste0("eps", 0:2)
                                 ))

accProb <- array(NA, dim = c(2, 3, 4, 10, 10, 3),
                 dimnames = list(paste0("n", 1:2),
                                 paste0("k", 1:3),
                                 paste0("f", 1:4),
                                 paste0("cv", 1:10),
                                 paste0("iter", 1:10),
                                 paste0("eps", 0:2)
                                 ))

mcmcseBeta <- array(NA, dim = c(2, 3, 4, 10, 6, 10, 3),
                 dimnames = list(paste0("n", 1:2),
                                 paste0("k", 1:3),
                                 paste0("f", 1:4),
                                 paste0("cv", 1:10),
                                 paste0("dim", 1:6),
                                 paste0("iter", 1:10),
                                 paste0("eps", 0:2)
                                 ))

mcmcseProb <- array(NA, dim = c(2, 3, 4, 10, 10, 3),
                 dimnames = list(paste0("n", 1:2),
                                 paste0("k", 1:3),
                                 paste0("f", 1:4),
                                 paste0("cv", 1:10),
                                 paste0("iter", 1:10),
                                 paste0("eps", 0:2)
                                 ))

times <- array(NA, dim = c(2, 3, 4, 10, 3, 2),
               dimnames = list(paste0("n", 1:2),
                               paste0("k", 1:3),
                               paste0("f", 1:4),
                               paste0("cv", 1:10),
                               paste0("eps", 0:2),
                               c("full", "part")
                               ))

for (nn in 1:2) {
    for (kk in 1:3) {
        for (ff in 1:4) {
            for (cc in 1:10) {
                fname <- paste0("/Shared/ssrivastva/dist_mov_len/result/logistic/part_n_", nn, "_k_", kk, "_frac_", nfrac[ff], "_rep_", cc, ".rds")
                if (file.exists(fname)) {
                    res <- readRDS(fname)
                    pbetas <- res$betas
                    pprob <- 1 / (1 + exp(-(drop(pbetas %*% c(rep(0, 4), 1, 0)))))
                    times[nn, kk, ff, cc, 1, "part"] <- res$time
                }

                if (all(!is.na(fullRes[[cc]][[nn]]))) {
                    times[nn, kk, ff, cc, 1, "full"] <- fullRes[[cc]][[nn]]$time
                }

                if (all(!is.na(fullRes[[cc]][[nn]])) & file.exists(fname)) {
                    cat("nn: ", nn, " kk: ", kk, " ff: ", ff, " cc: ", cc, "\n")

                    for (dd in 1:6) {
                        rr <- range(c(pbetas[100:10000, dd], fullRes[[cc]][[nn]]$beta[100:10000, dd]))
                        for (its in 1:10) {
                            pdensB <- bkde(x         = pbetas[100:(its * 1000), dd],
                                           bandwidth = dpik(pbetas[100:(its * 1000), dd], gridsize = 1000, range.x = rr),
                                           gridsize  = 1000,
                                           range.x   = rr
                                           )
                            fdensB <- bkde(x         = fullRes[[cc]][[nn]]$beta[100:(its * 1000), dd],
                                           bandwidth = dpik(fullRes[[cc]][[nn]]$beta[100:(its * 1000), dd],
                                                            gridsize = 1000, range.x = rr),
                                           gridsize = 1000,
                                           range.x   = rr
                                           )
                            accBeta[nn, kk, ff, cc, dd, its, 1] <- max(1 - sum(abs(pdensB$y - fdensB$y) * diff(pdensB$x)[1]) / 2, 0)
                            mcmcseBeta[nn, kk, ff, cc, dd, its, 1] <- mcse(pbetas[100:(its * 1000), dd], method = "obm")$se - mcse(fullRes[[cc]][[nn]]$beta[100:(its * 1000), dd], method = "obm")$se
                        }
                    }

                    rr <- range(c(pprob[100:10000], fullRes[[cc]][[nn]]$prob[100:10000]))
                    for (its in 1:10) {
                        pdensS <- bkde(x         = pprob[100:(its * 1000)],
                                       bandwidth = dpik(pprob[100:(its * 1000)], gridsize = 1000, range.x = rr),
                                       gridsize  = 1000,
                                       range.x   = rr
                                       )
                        fdensS <- bkde(x         = fullRes[[cc]][[nn]]$prob[100:(its * 1000)],
                                       bandwidth = dpik(fullRes[[cc]][[nn]]$prob[100:(its * 1000)], gridsize = 1000, range.x = rr),
                                       gridsize  = 1000,
                                       range.x   = rr
                                       )
                        accProb[nn, kk, ff, cc, its, 1] <- max(1 - sum(abs(pdensS$y - fdensS$y) * diff(pdensS$x)[1]) / 2, 0)
                        mcmcseProb[nn, kk, ff, cc, its, 1] <- mcse(pprob[100:(its * 1000)], method = "obm")$se - mcse(fullRes[[cc]][[nn]]$prob[100:(its * 1000)], method = "obm")$se
                    }
                }
            }
        }
    }
}

saveRDS(list(accBeta, mcmcseBeta, accProb, mcmcseProb), "~/dist_mov_len/data/log_eps0.rds")

for (nn in 1:2) {
    for (kk in 1:3) {
        for (ff in 1:4) {
            for (cc in 1:10) {
                fname <- paste0("/Shared/ssrivastva/dist_mov_len/result/logistic/part_n_", nn, "_k_", kk, "_frac_", nfrac[ff], "_rep_", cc, "_eps_0.1.rds")
                if (file.exists(fname)) {
                    res <- readRDS(fname)
                    pbetas <- res$betas
                    pprob <- 1 / (1 + exp(-(drop(pbetas %*% c(rep(0, 4), 1, 0)))))
                    times[nn, kk, ff, cc, 2, "part"] <- res$time
                }

                if (all(!is.na(fullRes[[cc]][[nn]]))) {
                    times[nn, kk, ff, cc, 2, "full"] <- fullRes[[cc]][[nn]]$time
                }

                if (all(!is.na(fullRes[[cc]][[nn]])) & file.exists(fname)) {
                    cat("nn: ", nn, " kk: ", kk, " ff: ", ff, " cc: ", cc, "\n")

                    for (dd in 1:6) {
                        rr <- range(c(pbetas[100:10000, dd], fullRes[[cc]][[nn]]$beta[100:10000, dd]))
                        for (its in 1:10) {
                            pdensB <- bkde(x         = pbetas[100:(its * 1000), dd],
                                           bandwidth = dpik(pbetas[100:(its * 1000), dd], gridsize = 1000, range.x = rr),
                                           gridsize  = 1000,
                                           range.x   = rr
                                           )
                            fdensB <- bkde(x         = fullRes[[cc]][[nn]]$beta[100:(its * 1000), dd],
                                           bandwidth = dpik(fullRes[[cc]][[nn]]$beta[100:(its * 1000), dd],
                                                            gridsize = 1000, range.x = rr),
                                           gridsize = 1000,
                                           range.x   = rr
                                           )
                            accBeta[nn, kk, ff, cc, dd, its, 2] <- max(1 - sum(abs(pdensB$y - fdensB$y) * diff(pdensB$x)[1]) / 2, 0)
                            mcmcseBeta[nn, kk, ff, cc, dd, its, 2] <- mcse(pbetas[100:(its * 1000), dd], method = "obm")$se - mcse(fullRes[[cc]][[nn]]$beta[100:(its * 1000), dd], method = "obm")$se
                        }
                    }

                    rr <- range(c(pprob[100:10000], fullRes[[cc]][[nn]]$prob[100:10000]))
                    for (its in 1:10) {
                        pdensS <- bkde(x         = pprob[100:(its * 1000)],
                                       bandwidth = dpik(pprob[100:(its * 1000)], gridsize = 1000, range.x = rr),
                                       gridsize  = 1000,
                                       range.x   = rr
                                       )
                        fdensS <- bkde(x         = fullRes[[cc]][[nn]]$prob[100:(its * 1000)],
                                       bandwidth = dpik(fullRes[[cc]][[nn]]$prob[100:(its * 1000)], gridsize = 1000, range.x = rr),
                                       gridsize  = 1000,
                                       range.x   = rr
                                       )
                        accProb[nn, kk, ff, cc, its, 2] <- max(1 - sum(abs(pdensS$y - fdensS$y) * diff(pdensS$x)[1]) / 2, 0)
                        mcmcseProb[nn, kk, ff, cc, its, 2] <- mcse(pprob[100:(its * 1000)], method = "obm")$se - mcse(fullRes[[cc]][[nn]]$prob[100:(its * 1000)], method = "obm")$se
                    }
                }
            }
        }
    }
}
saveRDS(list(accBeta, mcmcseBeta, accProb, mcmcseProb), "~/dist_mov_len/data/log_eps1.rds")

for (nn in 1:2) {
    for (kk in 1:3) {
        for (ff in 1:4) {
            for (cc in 1:10) {
                fname <- paste0("/Shared/ssrivastva/dist_mov_len/result/logistic/part_n_", nn, "_k_", kk, "_frac_", nfrac[ff], "_rep_", cc, "_eps_0.01.rds")
                if (file.exists(fname)) {
                    res <- readRDS(fname)
                    pbetas <- res$betas
                    pprob <- 1 / (1 + exp(-(drop(pbetas %*% c(rep(0, 4), 1, 0)))))
                    times[nn, kk, ff, cc, 3, "part"] <- res$time
                }

                if (all(!is.na(fullRes[[cc]][[nn]]))) {
                    times[nn, kk, ff, cc, 3, "full"] <- fullRes[[cc]][[nn]]$time
                }

                if (all(!is.na(fullRes[[cc]][[nn]])) & file.exists(fname)) {
                    cat("nn: ", nn, " kk: ", kk, " ff: ", ff, " cc: ", cc, "\n")

                    for (dd in 1:6) {
                        rr <- range(c(pbetas[100:10000, dd], fullRes[[cc]][[nn]]$beta[100:10000, dd]))
                        for (its in 1:10) {
                            pdensB <- bkde(x         = pbetas[100:(its * 1000), dd],
                                           bandwidth = dpik(pbetas[100:(its * 1000), dd], gridsize = 1000, range.x = rr),
                                           gridsize  = 1000,
                                           range.x   = rr
                                           )
                            fdensB <- bkde(x         = fullRes[[cc]][[nn]]$beta[100:(its * 1000), dd],
                                           bandwidth = dpik(fullRes[[cc]][[nn]]$beta[100:(its * 1000), dd],
                                                            gridsize = 1000, range.x = rr),
                                           gridsize = 1000,
                                           range.x   = rr
                                           )
                            accBeta[nn, kk, ff, cc, dd, its, 3] <- max(1 - sum(abs(pdensB$y - fdensB$y) * diff(pdensB$x)[1]) / 2, 0)
                            mcmcseBeta[nn, kk, ff, cc, dd, its, 3] <- mcse(pbetas[100:(its * 1000), dd], method = "obm")$se - mcse(fullRes[[cc]][[nn]]$beta[100:(its * 1000), dd], method = "obm")$se
                        }
                    }

                    rr <- range(c(pprob[100:10000], fullRes[[cc]][[nn]]$prob[100:10000]))
                    for (its in 1:10) {
                        pdensS <- bkde(x         = pprob[100:(its * 1000)],
                                       bandwidth = dpik(pprob[100:(its * 1000)], gridsize = 1000, range.x = rr),
                                       gridsize  = 1000,
                                       range.x   = rr
                                       )
                        fdensS <- bkde(x         = fullRes[[cc]][[nn]]$prob[100:(its * 1000)],
                                       bandwidth = dpik(fullRes[[cc]][[nn]]$prob[100:(its * 1000)], gridsize = 1000, range.x = rr),
                                       gridsize  = 1000,
                                       range.x   = rr
                                       )
                        accProb[nn, kk, ff, cc, its, 3] <- max(1 - sum(abs(pdensS$y - fdensS$y) * diff(pdensS$x)[1]) / 2, 0)
                        mcmcseProb[nn, kk, ff, cc, its, 3] <- mcse(pprob[100:(its * 1000)], method = "obm")$se - mcse(fullRes[[cc]][[nn]]$prob[100:(its * 1000)], method = "obm")$se
                    }
                }
            }
        }
    }
}
saveRDS(list(accBeta, mcmcseBeta, accProb, mcmcseProb), "~/dist_mov_len/data/log_eps2.rds")


res0 <- readRDS("~/dist_mov_len/data/log_eps0.rds")
res1 <- readRDS("~/dist_mov_len/data/log_eps1.rds")
res2 <- readRDS("~/dist_mov_len/data/log_eps2.rds")


accBeta[ , , , , , , 1] <- res0[[1]][ , , , , , , 1]
accBeta[ , , , , , , 2] <- res1[[1]][ , , , , , , 2]
accBeta[ , , , , , , 3] <- res2[[1]][ , , , , , , 3]
mcmcseBeta[ , , , , , , 1] <- abs(res0[[2]][ , , , , , , 1])
mcmcseBeta[ , , , , , , 2] <- abs(res1[[2]][ , , , , , , 2])
mcmcseBeta[ , , , , , , 3] <- abs(res2[[2]][ , , , , , , 3])

accProb[ , , , , , 1] <- res0[[3]][ , , , , , 1]
accProb[ , , , , , 2] <- res1[[3]][ , , , , , 2]
accProb[ , , , , , 3] <- res2[[3]][ , , , , , 3]
mcmcseProb[ , , , , , 1] <- abs(res0[[4]][ , , , , , 1])
mcmcseProb[ , , , , , 2] <- abs(res1[[4]][ , , , , , 2])
mcmcseProb[ , , , , , 3] <- abs(res2[[4]][ , , , , , 3])


accBetaEst <- apply(accBeta, c(1:3, 6, 7), mean, na.rm = TRUE)
accProbEst <- apply(accProb, c(1:3, 5, 6), mean, na.rm = TRUE)
mcmcseBetaEst <- apply(mcmcseBeta, c(1:3, 6, 7), mean, na.rm = TRUE)
mcmcseProbEst <- apply(mcmcseProb, c(1:3, 5, 6), mean, na.rm = TRUE)

library(plyr)
library(ggplot2)

accBetaEstDf <- adply(accBetaEst, 1:5)
colnames(accBetaEstDf) <- c("nsamp", "nsub", "frac", "niter", "eps", "est")
accProbEstDf <- adply(accProbEst, 1:5)
colnames(accProbEstDf) <- c("nsamp", "nsub", "frac", "niter", "eps", "est")

accBetaEstDf$niter <- as.numeric(gsub("iter", "", as.character(accBetaEstDf$niter)))
accProbEstDf$niter <- as.numeric(gsub("iter", "", as.character(accProbEstDf$niter)))


mcmcseBetaEstDf <- adply(mcmcseBetaEst, 1:5)
colnames(mcmcseBetaEstDf) <- c("nsamp", "nsub", "frac", "niter", "eps", "est")
mcmcseProbEstDf <- adply(mcmcseProbEst, 1:5)
colnames(mcmcseProbEstDf) <- c("nsamp", "nsub", "frac", "niter", "eps", "est")

mcmcseBetaEstDf$niter <- as.numeric(gsub("iter", "", as.character(mcmcseBetaEstDf$niter)))
mcmcseProbEstDf$niter <- as.numeric(gsub("iter", "", as.character(mcmcseProbEstDf$niter)))

plt <- ggplot(accBetaEstDf, aes(x = niter, y = est, group = interaction(frac, eps), colour = frac, shape = eps)) +
    xlab("Iterations (x 1000)") + ylab("Accuracy for" ~ bold(beta)) +
    scale_x_continuous(breaks = 1:10) +
    scale_colour_brewer(name="",
                        breaks=c("f1", "f2", "f3", "f4"),
                        labels=c("r = 0.20", "r = 0.40", "r = 0.60", "r = 0.80"),
                        palette="Set1",
                        guide = guide_legend(override.aes = list(shape = NA))) +
    scale_shape(name="",
                        breaks=c("eps0", "eps2", "eps1"),
                        labels=c("\u03B5 = 0.00", "\u03B5 = 0.01", "\u03B5 = 0.10")) +
    geom_line(stat = "identity", size = 4) + geom_point(size = 8) +
    facet_grid(nsamp ~ nsub,
               labeller = labeller(.cols = c(k1 = "k = 10", k2 = "k = 25", k3 = "k = 50"),
                                   .rows = c(n1 = "n ~ 1000,000", n2 = "n ~ 10,000,000")
                                   )
               ) +     theme_bw() +
    theme(strip.text.x = element_text(size=24, face="bold"),
          strip.text.y = element_text(size=24, face="bold"),
          axis.text.x = element_text(size=20, face="bold"),
          axis.text.y = element_text(size=20, face="bold"),
          axis.title.x = element_text(size=24, face="bold"),
          axis.title.y = element_text(size=24, face="bold"),
          legend.text = element_text(size = 24)
          ) +
    ylim(0.80, 1.00)


mplt <- ggplot(mcmcseBetaEstDf, aes(x = niter, y = est, group = interaction(frac, eps), colour = frac, shape = eps)) +
    xlab("Iterations (x 1000)") + ylab("SE Difference for" ~ bold(beta)) +
    scale_x_continuous(breaks = 1:10) +
    scale_colour_brewer(name="",
                        breaks=c("f1", "f2", "f3", "f4"),
                        labels=c("r = 0.20", "r = 0.40", "r = 0.60", "r = 0.80"),
                        palette="Set1",
                        guide = guide_legend(override.aes = list(shape = NA))) +
    scale_shape(name="",
                        breaks=c("eps0", "eps2", "eps1"),
                        labels=c("\u03B5 = 0.00", "\u03B5 = 0.01", "\u03B5 = 0.10")) +
    geom_line(stat = "identity", size = 4) + geom_point(size = 8) +
    facet_grid(nsamp ~ nsub,
               labeller = labeller(.cols = c(k1 = "k = 10", k2 = "k = 25", k3 = "k = 50"),
                                   .rows = c(n1 = "n ~ 1000,000", n2 = "n ~ 10,000,000")
                                   )
               ) +     theme_bw() +
    theme(strip.text.x = element_text(size=24, face="bold"),
          strip.text.y = element_text(size=24, face="bold"),
          axis.text.x = element_text(size=20, face="bold"),
          axis.text.y = element_text(size=20, face="bold"),
          axis.title.x = element_text(size=24, face="bold"),
          axis.title.y = element_text(size=24, face="bold"),
          legend.text = element_text(size = 24)
          ) +
    ylim(0, 0.0004)

plt1 <- ggplot(accProbEstDf, aes(x = niter, y = est, group = interaction(frac, eps), colour = frac, shape = eps)) +
    xlab("Iterations (x 1000)") + ylab("Accuracy for P(Y = 1 | X)") +
    scale_x_continuous(breaks = 1:10) +
    scale_colour_brewer(name="",
                        breaks=c("f1", "f2", "f3", "f4"),
                        labels=c("r = 0.20", "r = 0.40", "r = 0.60", "r = 0.80"),
                        palette="Set1",
                        guide = guide_legend(override.aes = list(shape = NA))) +
    scale_shape(name="",
                        breaks=c("eps0", "eps2", "eps1"),
                        labels=c("\u03B5 = 0.00", "\u03B5 = 0.01", "\u03B5 = 0.10")) +
    geom_line(stat = "identity", size = 4) + geom_point(size = 8) +
    facet_grid(nsamp ~ nsub,
               labeller = labeller(.cols = c(k1 = "k = 10", k2 = "k = 25", k3 = "k = 50"),
                                   .rows = c(n1 = "n ~ 1000,000", n2 = "n ~ 10,000,000")
                                   )
               ) +     theme_bw() +
    theme(strip.text.x = element_text(size=24, face="bold"),
          strip.text.y = element_text(size=24, face="bold"),
          axis.text.x = element_text(size=20, face="bold"),
          axis.text.y = element_text(size=20, face="bold"),
          axis.title.x = element_text(size=24, face="bold"),
          axis.title.y = element_text(size=24, face="bold"),
          legend.text = element_text(size = 24)
          ) +
    ylim(0.80, 1.00)

mplt1 <- ggplot(mcmcseProbEstDf, aes(x = niter, y = est, group = interaction(frac, eps), colour = frac, shape = eps)) +
    xlab("Iterations (x 1000)") + ylab("SE Difference for P(Y = 1 | X)") +
    scale_x_continuous(breaks = 1:10) +
    scale_colour_brewer(name="",
                        breaks=c("f1", "f2", "f3", "f4"),
                        labels=c("r = 0.20", "r = 0.40", "r = 0.60", "r = 0.80"),
                        palette="Set1",
                        guide = guide_legend(override.aes = list(shape = NA))) +
    scale_shape(name="",
                        breaks=c("eps0", "eps2", "eps1"),
                        labels=c("\u03B5 = 0.00", "\u03B5 = 0.01", "\u03B5 = 0.10")) +
    geom_line(stat = "identity", size = 4) + geom_point(size = 8) +
    facet_grid(nsamp ~ nsub,
               labeller = labeller(.cols = c(k1 = "k = 10", k2 = "k = 25", k3 = "k = 50"),
                                   .rows = c(n1 = "n ~ 1000,000", n2 = "n ~ 10,000,000")
                                   )
               ) +     theme_bw() +
    theme(strip.text.x = element_text(size=24, face="bold"),
          strip.text.y = element_text(size=24, face="bold"),
          axis.text.x = element_text(size=20, face="bold"),
          axis.text.y = element_text(size=20, face="bold"),
          axis.title.x = element_text(size=24, face="bold"),
          axis.title.y = element_text(size=24, face="bold"),
          legend.text = element_text(size = 24)
          ) +
    ylim(0, 0.0004)

cairo_pdf("~/dist_mov_len/data/log_mov_beta.pdf", family="DejaVu Sans", 14, 10)
print(plt)
dev.off()
cairo_pdf("~/dist_mov_len/data/log_mov_prob.pdf", family="DejaVu Sans", 14, 10)
print(plt1)
dev.off()

cairo_pdf("~/dist_mov_len/data/mcse_log_mov_beta.pdf", family="DejaVu Sans", 14, 10)
print(mplt)
dev.off()
cairo_pdf("~/dist_mov_len/data/mcse_log_mov_prob.pdf", family="DejaVu Sans", 14, 10)
print(mplt1)
dev.off()
