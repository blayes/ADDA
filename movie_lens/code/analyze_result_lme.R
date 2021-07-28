setwd("/Shared/ssrivastva/dist_mov_len/result/")
rm(list = ls())

library(KernSmooth)
library(plyr)
library(ggplot2)

nfrac <- c(20, 40, 60, 80)

### ACCURACY COMPUTATIONS
fullRes <- vector("list", 10)
for (cc in 1:10) {
    fullRes[[cc]] <- vector("list", 2)
    for (nn in 1:2) {
        if (nn == 1) {
            fname <- paste0("full/ml_rep_", cc, "_100k.rds")
        } else {
            fname <- paste0("full/ml_rep_", cc, "_1000k.rds")
        }

        if (file.exists(fname)) {
            fres <- readRDS(fname)
            fbetas <- fres$betas
            fdmat <- cbind(do.call(rbind, lapply(fres$tmat, function(x) diag(x))),
                           do.call(rbind, lapply(fres$tmat, function(x) x[lower.tri(x)]))
                           )
            fsig <- fres$sigma

            fullRes[[cc]][[nn]] <- list(
                beta = fbetas,
                dmat = fdmat,
                sig = fsig,
                time = fres$time)
        } else {
            fullRes[[cc]][[nn]] <- NA
        }
    }
}
## saveRDS(fullRes, "full/full_summ.rds")
## fullRes <- readRDS("full/full_summ.rds")

accBeta <- array(NA, dim = c(2, 3, 4, 10, 6, 10, 3),
                 dimnames = list(paste0("n", 1:2),
                                 paste0("k", 1:3),
                                 paste0("f", 1:4),
                                 paste0("cv", 1:10),
                                 paste0("dim", 1:6),
                                 paste0("iter", 1:10),
                                 paste0("eps", 0:2)
                                 ))

accDmat <- array(NA, dim = c(2, 3, 4, 10, 21, 10, 3),
                 dimnames = list(paste0("n", 1:2),
                                 paste0("k", 1:3),
                                 paste0("f", 1:4),
                                 paste0("cv", 1:10),
                                 paste0("dim", 1:21),
                                 paste0("iter", 1:10),
                                 paste0("eps", 0:2)
                                 ))

accSig <- array(NA, dim = c(2, 3, 4, 10, 10, 3),
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
                fname <- paste0("sub/ml_part_n_", nn, "_k_", kk, "_frac_", nfrac[ff], "_rep_", cc, ".rds")

                if (file.exists(fname)) {
                    res <- readRDS(fname)
                    pbetas <- res$betas
                    pdmat <- res$tmat
                    psig <- res$sigma
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
                        }
                    }

                    for (dd in 1:21) {
                        rr <- range(c(pdmat[100:10000, dd], fullRes[[cc]][[nn]]$dmat[100:10000, dd]))
                        for (its in 1:10) {
                            pdensD <- bkde(x         = pdmat[100:(its * 1000), dd],
                                           bandwidth = dpik(pdmat[100:(its * 1000), dd], gridsize = 1000, range.x = rr),
                                           gridsize  = 1000,
                                           range.x   = rr
                                           )
                            fdensD <- bkde(x         = fullRes[[cc]][[nn]]$dmat[100:(its * 1000), dd],
                                           bandwidth = dpik(fullRes[[cc]][[nn]]$dmat[100:(its * 1000), dd], gridsize = 1000, range.x = rr),
                                           gridsize  = 1000,
                                           range.x   = rr
                                           )
                            accDmat[nn, kk, ff, cc, dd, its, 1] <- max(1 - sum(abs(pdensD$y - fdensD$y) * diff(pdensD$x)[1]) / 2, 0)
                        }
                    }

                    rr <- range(c(psig[100:10000], fullRes[[cc]][[nn]]$sig[100:10000]))
                    for (its in 1:10) {
                        pdensS <- bkde(x         = psig[100:(its * 1000)],
                                       bandwidth = dpik(psig[100:(its * 1000)], gridsize = 1000, range.x = rr),
                                       gridsize  = 1000,
                                       range.x   = rr
                                       )
                        fdensS <- bkde(x         = fullRes[[cc]][[nn]]$sig[100:(its * 1000)],
                                       bandwidth = dpik(fullRes[[cc]][[nn]]$sig[100:(its * 1000)], gridsize = 1000, range.x = rr),
                                       gridsize  = 1000,
                                       range.x   = rr
                                       )
                        accSig[nn, kk, ff, cc, its, 1] <- max(1 - sum(abs(pdensS$y - fdensS$y) * diff(pdensS$x)[1]) / 2, 0)
                    }
                }
            }
        }
    }
}


for (nn in 1:2) {
    for (kk in 1:3) {
        for (ff in 1:4) {
            for (cc in 1:10) {
                fname <- paste0("sub/ml_part_n_", nn, "_k_", kk, "_frac_", nfrac[ff], "_rep_", cc, "_eps_0.1.rds")

                if (file.exists(fname)) {
                    res <- readRDS(fname)
                    pbetas <- res$betas
                    pdmat <- res$tmat
                    psig <- res$sigma
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
                        }
                    }

                    for (dd in 1:21) {
                        rr <- range(c(pdmat[100:10000, dd], fullRes[[cc]][[nn]]$dmat[100:10000, dd]))
                        for (its in 1:10) {
                            pdensD <- bkde(x         = pdmat[100:(its * 1000), dd],
                                           bandwidth = dpik(pdmat[100:(its * 1000), dd], gridsize = 1000, range.x = rr),
                                           gridsize  = 1000,
                                           range.x   = rr
                                           )
                            fdensD <- bkde(x         = fullRes[[cc]][[nn]]$dmat[100:(its * 1000), dd],
                                           bandwidth = dpik(fullRes[[cc]][[nn]]$dmat[100:(its * 1000), dd], gridsize = 1000, range.x = rr),
                                           gridsize  = 1000,
                                           range.x   = rr
                                           )
                            accDmat[nn, kk, ff, cc, dd, its, 2] <- max(1 - sum(abs(pdensD$y - fdensD$y) * diff(pdensD$x)[1]) / 2, 0)
                        }
                    }

                    rr <- range(c(psig[100:10000], fullRes[[cc]][[nn]]$sig[100:10000]))
                    for (its in 1:10) {
                        pdensS <- bkde(x         = psig[100:(its * 1000)],
                                       bandwidth = dpik(psig[100:(its * 1000)], gridsize = 1000, range.x = rr),
                                       gridsize  = 1000,
                                       range.x   = rr
                                       )
                        fdensS <- bkde(x         = fullRes[[cc]][[nn]]$sig[100:(its * 1000)],
                                       bandwidth = dpik(fullRes[[cc]][[nn]]$sig[100:(its * 1000)], gridsize = 1000, range.x = rr),
                                       gridsize  = 1000,
                                       range.x   = rr
                                       )
                        accSig[nn, kk, ff, cc, its, 2] <- max(1 - sum(abs(pdensS$y - fdensS$y) * diff(pdensS$x)[1]) / 2, 0)
                    }
                }
            }
        }
    }
}


for (nn in 1:2) {
    for (kk in 1:3) {
        for (ff in 1:4) {
            for (cc in 1:10) {
                fname <- paste0("sub/ml_part_n_", nn, "_k_", kk, "_frac_", nfrac[ff], "_rep_", cc, "_eps_0.01.rds")

                if (file.exists(fname)) {
                    res <- readRDS(fname)
                    pbetas <- res$betas
                    pdmat <- res$tmat
                    psig <- res$sigma
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
                        }
                    }

                    for (dd in 1:21) {
                        rr <- range(c(pdmat[100:10000, dd], fullRes[[cc]][[nn]]$dmat[100:10000, dd]))
                        for (its in 1:10) {
                            pdensD <- bkde(x         = pdmat[100:(its * 1000), dd],
                                           bandwidth = dpik(pdmat[100:(its * 1000), dd], gridsize = 1000, range.x = rr),
                                           gridsize  = 1000,
                                           range.x   = rr
                                           )
                            fdensD <- bkde(x         = fullRes[[cc]][[nn]]$dmat[100:(its * 1000), dd],
                                           bandwidth = dpik(fullRes[[cc]][[nn]]$dmat[100:(its * 1000), dd], gridsize = 1000, range.x = rr),
                                           gridsize  = 1000,
                                           range.x   = rr
                                           )
                            accDmat[nn, kk, ff, cc, dd, its, 3] <- max(1 - sum(abs(pdensD$y - fdensD$y) * diff(pdensD$x)[1]) / 2, 0)
                        }
                    }

                    rr <- range(c(psig[100:10000], fullRes[[cc]][[nn]]$sig[100:10000]))
                    for (its in 1:10) {
                        pdensS <- bkde(x         = psig[100:(its * 1000)],
                                       bandwidth = dpik(psig[100:(its * 1000)], gridsize = 1000, range.x = rr),
                                       gridsize  = 1000,
                                       range.x   = rr
                                       )
                        fdensS <- bkde(x         = fullRes[[cc]][[nn]]$sig[100:(its * 1000)],
                                       bandwidth = dpik(fullRes[[cc]][[nn]]$sig[100:(its * 1000)], gridsize = 1000, range.x = rr),
                                       gridsize  = 1000,
                                       range.x   = rr
                                       )
                        accSig[nn, kk, ff, cc, its, 3] <- max(1 - sum(abs(pdensS$y - fdensS$y) * diff(pdensS$x)[1]) / 2, 0)
                    }
                }
            }
        }
    }
}

accBetaEst <- apply(accBeta, c(1:3, 6, 7), mean, na.rm = TRUE)
accDmatEst <- apply(accDmat, c(1:3, 6, 7), mean, na.rm = TRUE)
accSigEst <- apply(accSig, c(1:3, 5, 6), mean, na.rm = TRUE)
accBetaEstDf <- adply(accBetaEst, 1:5)
accDmatEstDf <- adply(accDmatEst, 1:5)
accSigEstDf <- adply(accSigEst, 1:5)

colnames(accBetaEstDf) <- colnames(accDmatEstDf) <- colnames(accSigEstDf) <- c("nsamp", "nsub", "frac", "niter", "eps", "est")

accBetaEstDf$niter <- as.numeric(gsub("iter", "", as.character(accBetaEstDf$niter)))
accDmatEstDf$niter <- as.numeric(gsub("iter", "", as.character(accDmatEstDf$niter)))
accSigEstDf$niter <- as.numeric(gsub("iter", "", as.character(accSigEstDf$niter)))

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
    facet_grid(nsamp  ~ nsub,
               labeller = labeller(.cols = c(k1 = "k = 10", k2 = "k = 25", k3 = "k = 50"),
                                   .rows = c(n1 = "n ~ 100,000", n2 = "n ~ 1,000,000")
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
    ylim(0.40, 1.00)


plt1 <- ggplot(accDmatEstDf, aes(x = niter, y = est, group = interaction(frac, eps), colour = frac, shape = eps)) +
    xlab("Iterations (x 1000)") + ylab("Accuracy for" ~ bold(Sigma)) +
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
    facet_grid(nsamp  ~ nsub,
               labeller = labeller(.cols = c(k1 = "k = 10", k2 = "k = 25", k3 = "k = 50"),
                                   .rows = c(n1 = "n ~ 100,000", n2 = "n ~ 1,000,000")
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
    ylim(0.40, 1.00)


plt2 <- ggplot(accSigEstDf, aes(x = niter, y = est, group = interaction(frac, eps), colour = frac, shape = eps)) +
    xlab("Iterations (x 1000)") + ylab(expression("Accuracy for" ~ sigma^2)) +
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
    facet_grid(nsamp  ~ nsub,
               labeller = labeller(.cols = c(k1 = "k = 10", k2 = "k = 25", k3 = "k = 50"),
                                   .rows = c(n1 = "n ~ 100,000", n2 = "n ~ 1,000,000")
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
    ylim(0.40, 1.00)


cairo_pdf("~/dist_mov_len/data/mov_beta.pdf", family="DejaVu Sans", 15, 8)
print(plt)
dev.off()
cairo_pdf("~/dist_mov_len/data/mov_dmat.pdf", family="DejaVu Sans", 15, 8)
print(plt1)
dev.off()
cairo_pdf("~/dist_mov_len/data/mov_sig.pdf", family="DejaVu Sans", 15, 8)
print(plt2)
dev.off()

## TIME COMPARISON COMPUTATIONS
for (nn in 1:2) {
    for (kk in 1:3) {
        for (ff in 1:4) {
            for (cc in 1:10) {
                fname <- paste0("sub/ml_part_n_", nn, "_k_", kk, "_frac_", nfrac[ff], "_rep_", cc, ".rds")

                if (file.exists(fname)) {
                    res <- readRDS(fname)
                    pbetas <- res$betas
                    pdmat <- res$tmat
                    psig <- res$sigma
                    times[nn, kk, ff, cc, 1, "part"] <- res$time
                }

                if (all(!is.na(fullRes[[cc]][[nn]]))) {
                    times[nn, kk, ff, cc, 1, "full"] <- fullRes[[cc]][[nn]]$time
                }
            }
        }
    }
}

for (nn in 1:2) {
    for (kk in 1:3) {
        for (ff in 1:4) {
            for (cc in 1:10) {
                fname <- paste0("sub/ml_part_n_", nn, "_k_", kk, "_frac_", nfrac[ff], "_rep_", cc, "_eps_0.1.rds")

                if (file.exists(fname)) {
                    res <- readRDS(fname)
                    pbetas <- res$betas
                    pdmat <- res$tmat
                    psig <- res$sigma
                    times[nn, kk, ff, cc, 2, "part"] <- res$time
                }

                if (all(!is.na(fullRes[[cc]][[nn]]))) {
                    times[nn, kk, ff, cc, 2, "full"] <- fullRes[[cc]][[nn]]$time
                }
            }
        }
    }
}


for (nn in 1:2) {
    for (kk in 1:3) {
        for (ff in 1:4) {
            for (cc in 1:10) {
                fname <- paste0("sub/ml_part_n_", nn, "_k_", kk, "_frac_", nfrac[ff], "_rep_", cc, "_eps_0.01.rds")

                if (file.exists(fname)) {
                    res <- readRDS(fname)
                    pbetas <- res$betas
                    pdmat <- res$tmat
                    psig <- res$sigma
                    times[nn, kk, ff, cc, 3, "part"] <- res$time
                }

                if (all(!is.na(fullRes[[cc]][[nn]]))) {
                    times[nn, kk, ff, cc, 3, "full"] <- fullRes[[cc]][[nn]]$time
                }
            }
        }
    }
}

timeComp <- apply(times, c(1:3, 5, 6), mean, na.rm = TRUE)

fr1.0 <- as.numeric(timeComp[ , , "f1", "eps0", "full"] / timeComp[ , , "f1", "eps0", "part"])
fr2.0 <- as.numeric(timeComp[ , , "f2", "eps0", "full"] / timeComp[ , , "f2", "eps0", "part"])
fr3.0 <- as.numeric(timeComp[ , , "f3", "eps0", "full"] / timeComp[ , , "f3", "eps0", "part"])
fr4.0 <- as.numeric(timeComp[ , , "f4", "eps0", "full"] / timeComp[ , , "f4", "eps0", "part"])
fr1.1 <- as.numeric(timeComp[ , , "f1", "eps1", "full"] / timeComp[ , , "f1", "eps1", "part"])
fr2.1 <- as.numeric(timeComp[ , , "f2", "eps1", "full"] / timeComp[ , , "f2", "eps1", "part"])
fr3.1 <- as.numeric(timeComp[ , , "f3", "eps1", "full"] / timeComp[ , , "f3", "eps1", "part"])
fr4.1 <- as.numeric(timeComp[ , , "f4", "eps1", "full"] / timeComp[ , , "f4", "eps1", "part"])
fr1.2 <- as.numeric(timeComp[ , , "f1", "eps2", "full"] / timeComp[ , , "f1", "eps2", "part"])
fr2.2 <- as.numeric(timeComp[ , , "f2", "eps2", "full"] / timeComp[ , , "f2", "eps2", "part"])
fr3.2 <- as.numeric(timeComp[ , , "f3", "eps2", "full"] / timeComp[ , , "f3", "eps2", "part"])
fr4.2 <- as.numeric(timeComp[ , , "f4", "eps2", "full"] / timeComp[ , , "f4", "eps2", "part"])

ns <- rep(1:2, times = 3)
ks <- rep(c(10, 25, 50), each = 2)

gaindf <- rbind.data.frame(data.frame(nsamp = rep(paste0("n", ns), 4),
                                      nsub = rep(ks, 4),
                                      frac = c(rep("f1", 6), rep("f2", 6), rep("f3", 6), rep("f4", 6)),
                                      imp = c(fr1.0, fr2.0, fr3.0, fr4.0),
                                      eps = "eps0"),
                           data.frame(nsamp = rep(paste0("n", ns), 4),
                                      nsub = rep(ks, 4),
                                      frac = c(rep("f1", 6), rep("f2", 6), rep("f3", 6), rep("f4", 6)),
                                      imp = c(fr1.1, fr2.1, fr3.1, fr4.1),
                                      eps = "eps1"),
                           data.frame(nsamp = rep(paste0("n", ns), 4),
                                      nsub = rep(ks, 4),
                                      frac = c(rep("f1", 6), rep("f2", 6), rep("f3", 6), rep("f4", 6)),
                                      imp = c(fr1.2, fr2.2, fr3.2, fr4.2),
                                      eps = "eps2")
                           )

pltt <- ggplot(gaindf, aes(x = nsub, y = imp, group = interaction(frac, eps), colour = frac, shape = eps)) +
    xlab("Number of Subsets (k)") + ylab("Gain in Run Time") +
    scale_x_continuous(breaks = c(10, 25, 50), labels = c("10", "25", "50")) +
    scale_y_continuous(breaks = 1:6) +
    scale_colour_brewer(name="",
                        breaks=c("f1", "f2", "f3", "f4"),
                        labels=c("r = 0.20", "r = 0.40", "r = 0.60", "r = 0.80"),
                        palette="Set1",
                        guide = guide_legend(override.aes = list(shape = NA))) +
    scale_shape(name="",
                        breaks=c("eps0",  "eps2", "eps1"),
                        labels=c("\u03B5 = 0.00", "\u03B5 = 0.01", "\u03B5 = 0.10")) +
    geom_line(stat = "identity", size = 4) + geom_point(size = 8) +
    facet_grid(. ~ nsamp,
               labeller = labeller(
                   .cols = c(n1 = "n ~ 100,000", n2 = "n ~ 1,000,000")
               )
               ) +  theme_bw() +
    geom_hline(yintercept=1, linetype="dashed",
               color = "black", size=4) +
    theme(strip.text.x = element_text(size=24, face="bold"),
          strip.text.y = element_text(size=24, face="bold"),
          axis.text.x = element_text(size=20, face="bold"),
          axis.text.y = element_text(size=20, face="bold"),
          axis.title.x = element_text(size=24, face="bold"),
          axis.title.y = element_text(size=24, face="bold"),
          legend.text = element_text(size = 24)
          )

cairo_pdf("~/dist_mov_len/data/mov_time.pdf", family="DejaVu Sans", 20, 6)
print(pltt)
dev.off()

## MCSE COMPUTATIONS
setwd("/Shared/ssrivastva/dist_mov_len/result/")

library(mcmcse)

nfrac <- c(20, 40, 60, 80)

mcmcseBeta <- array(NA, dim = c(2, 3, 4, 10, 6, 10, 3),
                 dimnames = list(paste0("n", 1:2),
                                 paste0("k", 1:3),
                                 paste0("f", 1:4),
                                 paste0("cv", 1:10),
                                 paste0("dim", 1:6),
                                 paste0("iter", 1:10),
                                 paste0("eps", 0:2)
                                 ))

mcmcseDmat <- array(NA, dim = c(2, 3, 4, 10, 21, 10, 3),
                 dimnames = list(paste0("n", 1:2),
                                 paste0("k", 1:3),
                                 paste0("f", 1:4),
                                 paste0("cv", 1:10),
                                 paste0("dim", 1:21),
                                 paste0("iter", 1:10),
                                 paste0("eps", 0:2)
                                 ))

mcmcseSig <- array(NA, dim = c(2, 3, 4, 10, 10, 3),
                 dimnames = list(paste0("n", 1:2),
                                 paste0("k", 1:3),
                                 paste0("f", 1:4),
                                 paste0("cv", 1:10),
                                 paste0("iter", 1:10),
                                 paste0("eps", 0:2)
                                 ))

for (nn in 1:2) {
    for (kk in 1:3) {
        for (ff in 1:4) {
            for (cc in 1:10) {
                fname <- paste0("sub/ml_part_n_", nn, "_k_", kk, "_frac_", nfrac[ff], "_rep_", cc, ".rds")

                if (file.exists(fname)) {
                    res <- readRDS(fname)
                    pbetas <- res$betas
                    pdmat <- res$tmat
                    psig <- res$sigma
                }

                if (all(!is.na(fullRes[[cc]][[nn]])) & file.exists(fname)) {
                    cat("nn: ", nn, " kk: ", kk, " ff: ", ff, " cc: ", cc, "\n")

                    for (dd in 1:6) {
                        for (its in 1:10) {
                            mcmcseBeta[nn, kk, ff, cc, dd, its, 1] <- mcse(pbetas[100:(its * 1000), dd], method = "obm")$se - mcse(fullRes[[cc]][[nn]]$beta[100:(its * 1000), dd], method = "obm")$se
                        }
                    }

                    for (dd in 1:21) {
                        for (its in 1:10) {
                            mcmcseDmat[nn, kk, ff, cc, dd, its, 1] <- mcse(pdmat[100:(its * 1000), dd], method = "obm")$se - mcse(fullRes[[cc]][[nn]]$dmat[100:(its * 1000), dd], method = "obm")$se
                        }
                    }

                    for (its in 1:10) {
                        mcmcseSig[nn, kk, ff, cc, its, 1] <- mcse(psig[100:(its * 1000)], method = "obm")$se - mcse(fullRes[[cc]][[nn]]$sig[100:(its * 1000)], method = "obm")$se
                    }
                }
            }
        }
    }
}

saveRDS(list(mcmcseBeta, mcmcseDmat, mcmcseSig), "~/dist_mov_len/data/mcse_lme_eps0.rds")

for (nn in 1:2) {
    for (kk in 1:3) {
        for (ff in 1:4) {
            for (cc in 1:10) {
                fname <- paste0("sub/ml_part_n_", nn, "_k_", kk, "_frac_", nfrac[ff], "_rep_", cc, "_eps_0.1.rds")

                if (file.exists(fname)) {
                    res <- readRDS(fname)
                    pbetas <- res$betas
                    pdmat <- res$tmat
                    psig <- res$sigma
                }

                if (all(!is.na(fullRes[[cc]][[nn]])) & file.exists(fname)) {
                    cat("nn: ", nn, " kk: ", kk, " ff: ", ff, " cc: ", cc, "\n")

                    for (dd in 1:6) {
                        for (its in 1:10) {
                            mcmcseBeta[nn, kk, ff, cc, dd, its, 2] <- mcse(pbetas[100:(its * 1000), dd], method = "obm")$se - mcse(fullRes[[cc]][[nn]]$beta[100:(its * 1000), dd], method = "obm")$se
                        }
                    }

                    for (dd in 1:21) {
                        for (its in 1:10) {
                            mcmcseDmat[nn, kk, ff, cc, dd, its, 2] <- mcse(pdmat[100:(its * 1000), dd], method = "obm")$se - mcse(fullRes[[cc]][[nn]]$dmat[100:(its * 1000), dd], method = "obm")$se
                        }
                    }

                    for (its in 1:10) {
                        mcmcseSig[nn, kk, ff, cc, its, 2] <- mcse(psig[100:(its * 1000)], method = "obm")$se - mcse(fullRes[[cc]][[nn]]$sig[100:(its * 1000)], method = "obm")$se
                    }
                }
            }
        }
    }
}
saveRDS(list(mcmcseBeta, mcmcseDmat, mcmcseSig), "~/dist_mov_len/data/mcse_lme_eps1.rds")






for (nn in 1:2) {
    for (kk in 1:3) {
        for (ff in 1:4) {
            for (cc in 1:10) {
                fname <- paste0("sub/ml_part_n_", nn, "_k_", kk, "_frac_", nfrac[ff], "_rep_", cc, "_eps_0.01.rds")

                if (file.exists(fname)) {
                    res <- readRDS(fname)
                    pbetas <- res$betas
                    pdmat <- res$tmat
                    psig <- res$sigma
                }

                if (all(!is.na(fullRes[[cc]][[nn]])) & file.exists(fname)) {
                    cat("nn: ", nn, " kk: ", kk, " ff: ", ff, " cc: ", cc, "\n")

                    for (dd in 1:6) {
                        for (its in 1:10) {
                            mcmcseBeta[nn, kk, ff, cc, dd, its, 3] <- mcse(pbetas[100:(its * 1000), dd], method = "obm")$se - mcse(fullRes[[cc]][[nn]]$beta[100:(its * 1000), dd], method = "obm")$se
                        }
                    }

                    for (dd in 1:21) {
                        for (its in 1:10) {
                            mcmcseDmat[nn, kk, ff, cc, dd, its, 3] <- mcse(pdmat[100:(its * 1000), dd], method = "obm")$se - mcse(fullRes[[cc]][[nn]]$dmat[100:(its * 1000), dd], method = "obm")$se
                        }
                    }

                    for (its in 1:10) {
                        mcmcseSig[nn, kk, ff, cc, its, 3] <- mcse(psig[100:(its * 1000)], method = "obm")$se - mcse(fullRes[[cc]][[nn]]$sig[100:(its * 1000)], method = "obm")$se
                    }
                }
            }
        }
    }
}


saveRDS(list(mcmcseBeta, mcmcseDmat, mcmcseSig), "~/dist_mov_len/data/mcse_lme_eps2.rds")

res0 <- readRDS("~/dist_mov_len/data/mcse_lme_eps0.rds")
res1 <- readRDS("~/dist_mov_len/data/mcse_lme_eps1.rds")
res2 <- readRDS("~/dist_mov_len/data/mcse_lme_eps2.rds")

mcmcseBeta[ , , , , , , 1] <- abs(res0[[1]][ , , , , , , 1])
mcmcseBeta[ , , , , , , 2] <- abs(res1[[1]][ , , , , , , 2])
mcmcseBeta[ , , , , , , 3] <- abs(res2[[1]][ , , , , , , 3])

mcmcseDmat[ , , , , , , 1] <- abs(res0[[2]][ , , , , , , 1])
mcmcseDmat[ , , , , , , 2] <- abs(res1[[2]][ , , , , , , 2])
mcmcseDmat[ , , , , , , 3] <- abs(res2[[2]][ , , , , , , 3])

mcmcseSig[ , , , , , 1] <- abs(res0[[3]][ , , , , , 1])
mcmcseSig[ , , , , , 2] <- abs(res1[[3]][ , , , , , 2])
mcmcseSig[ , , , , , 3] <- abs(res2[[3]][ , , , , , 3])

betaEst <- apply(mcmcseBeta, c(1:3, 6, 7), mean, na.rm = TRUE)
dmatEst <- apply(mcmcseDmat, c(1:3, 6, 7), mean, na.rm = TRUE)
sigEst <- apply(mcmcseSig, c(1:3, 5, 6), mean, na.rm = TRUE)

betaEstDf <- adply(betaEst, 1:5)
dmatEstDf <- adply(dmatEst, 1:5)
sigEstDf <- adply(sigEst, 1:5)

colnames(betaEstDf) <- colnames(dmatEstDf) <- colnames(sigEstDf) <- c("nsamp", "nsub", "frac", "niter", "eps", "est")

betaEstDf$niter <- as.numeric(gsub("iter", "", as.character(betaEstDf$niter)))
dmatEstDf$niter <- as.numeric(gsub("iter", "", as.character(dmatEstDf$niter)))
sigEstDf$niter <- as.numeric(gsub("iter", "", as.character(sigEstDf$niter)))

plt <- ggplot(betaEstDf, aes(x = niter, y = est, group = interaction(frac, eps), colour = frac, shape = eps)) +
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
    facet_grid(nsamp  ~ nsub,
               labeller = labeller(.cols = c(k1 = "k = 10", k2 = "k = 25", k3 = "k = 50"),
                                   .rows = c(n1 = "n ~ 100,000", n2 = "n ~ 1,000,000")
                                   )
               ) +     theme_bw() +
    theme(strip.text.x = element_text(size=24, face="bold"),
          strip.text.y = element_text(size=24, face="bold"),
          axis.text.x = element_text(size=20, face="bold"),
          axis.text.y = element_text(size=20, face="bold"),
          axis.title.x = element_text(size=24, face="bold"),
          axis.title.y = element_text(size=24, face="bold"),
          legend.text = element_text(size = 24)
          ) + ylim(0, 0.03)


plt1 <- ggplot(dmatEstDf, aes(x = niter, y = est, group = interaction(frac, eps), colour = frac, shape = eps)) +
    xlab("Iterations (x 1000)") + ylab("SE Difference for" ~ bold(Sigma)) +
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
    facet_grid(nsamp  ~ nsub,
               labeller = labeller(.cols = c(k1 = "k = 10", k2 = "k = 25", k3 = "k = 50"),
                                   .rows = c(n1 = "n ~ 100,000", n2 = "n ~ 1,000,000")
                                   )
               ) +     theme_bw() +
    theme(strip.text.x = element_text(size=24, face="bold"),
          strip.text.y = element_text(size=24, face="bold"),
          axis.text.x = element_text(size=20, face="bold"),
          axis.text.y = element_text(size=20, face="bold"),
          axis.title.x = element_text(size=24, face="bold"),
          axis.title.y = element_text(size=24, face="bold"),
          legend.text = element_text(size = 24)
          ) + ylim(0, 0.03)



plt2 <- ggplot(sigEstDf, aes(x = niter, y = est, group = interaction(frac, eps), colour = frac, shape = eps)) +
    xlab("Iterations (x 1000)") + ylab(expression("SE Difference for" ~ sigma^2)) +
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
    facet_grid(nsamp  ~ nsub,
               labeller = labeller(.cols = c(k1 = "k = 10", k2 = "k = 25", k3 = "k = 50"),
                                   .rows = c(n1 = "n ~ 100,000", n2 = "n ~ 1,000,000")
                                   )
               ) +     theme_bw() +
    theme(strip.text.x = element_text(size=24, face="bold"),
          strip.text.y = element_text(size=24, face="bold"),
          axis.text.x = element_text(size=20, face="bold"),
          axis.text.y = element_text(size=20, face="bold"),
          axis.title.x = element_text(size=24, face="bold"),
          axis.title.y = element_text(size=24, face="bold"),
          legend.text = element_text(size = 24)
          ) + ylim(0, 0.03)


cairo_pdf("~/dist_mov_len/data/mcse_mov_beta.pdf", family="DejaVu Sans", 15, 8)
print(plt)
dev.off()
cairo_pdf("~/dist_mov_len/data/mcse_mov_dmat.pdf", family="DejaVu Sans", 15, 8)
print(plt1)
dev.off()
cairo_pdf("~/dist_mov_len/data/mcse_mov_sig.pdf", family="DejaVu Sans", 15, 8)
print(plt2)
dev.off()
