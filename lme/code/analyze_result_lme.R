setwd("/Shared/ssrivastva/dist_lme/result/sim/")
rm(list = ls())

library(KernSmooth)
library(mcmcse)
library(plyr)
library(ggplot2)

nfrac <- c(20, 40, 60, 80)
## ACCURACY COMPUTATIONS
accBeta <- array(NA, dim = c(3, 3, 4, 10, 4, 20),
                 dimnames = list(paste0("n", 1:3),
                                 paste0("k", 1:3),
                                 paste0("f", 1:4),
                                 paste0("cv", 1:10),
                                 paste0("dim", 1:4),
                                 paste0("iter", 1:20)
                                 ))

accDmat <- array(NA, dim = c(3, 3, 4, 10, 6, 20),
                 dimnames = list(paste0("n", 1:3),
                                 paste0("k", 1:3),
                                 paste0("f", 1:4),
                                 paste0("cv", 1:10),
                                 paste0("dim", 1:6),
                                 paste0("iter", 1:20)
                                 ))

accSig <- array(NA, dim = c(3, 3, 4, 10, 20),
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

for (nn in 1:3) {
    for (kk in 1:3) {
        for (ff in 1:4) {
            for (cc in 1:10) {
                cat("nn: ", nn, " kk: ", kk, " ff: ", ff, " cc: ", cc, "\n")
                fname <- paste0("part/px_sim_part_n_", nn, "_k_", kk, "_frac_", nfrac[ff], "_rep_", cc, ".rds")

                res <- readRDS(fname)
                ffname <- paste0("full/sim_full_n_", nn, "_rep_", cc, ".rds")

                fres <- readRDS(ffname)
                pbetas <- res$betas
                fbetas <- fres$betas
                pdmat <- res$tmat
                fdmat <- cbind(do.call(rbind, lapply(fres$tmat, function(x) diag(x))),
                               do.call(rbind, lapply(fres$tmat, function(x) x[lower.tri(x)]))
                               )
                psig <- res$sigma
                fsig <- fres$sigma
                times[nn, kk, ff, cc, "full"] <- fres$time
                times[nn, kk, ff, cc, "part"] <- res$time
                for (dd in 1:4) {
                    rr <- range(c(pbetas[100:20000, dd], fbetas[100:20000, dd]))
                    for (its in 1:20) {
                        pdensB <- bkde(x         = pbetas[100:(its * 1000), dd],
                                       bandwidth = dpik(pbetas[100:(its * 1000), dd], gridsize = 1000, range.x = rr),
                                       gridsize  = 1000,
                                       range.x   = rr
                                       )
                        fdensB <- bkde(x         = fbetas[100:(its * 1000), dd],
                                       bandwidth = dpik(fbetas[100:(its * 1000), dd], gridsize = 1000, range.x = rr),
                                       gridsize = 1000,
                                       range.x   = rr
                                       )
                        accBeta[nn, kk, ff, cc, dd, its] <- max(1 - sum(abs(pdensB$y - fdensB$y) * diff(pdensB$x)[1]) / 2, 0)
                    }
                }

                for (dd in 1:6) {
                    rr <- range(c(pdmat[100:20000, dd], fdmat[100:20000, dd]))
                    for (its in 1:20) {
                        pdensD <- bkde(x         = pdmat[100:(its * 1000), dd],
                                       bandwidth = dpik(pdmat[100:(its * 1000), dd], gridsize = 1000, range.x = rr),
                                       gridsize  = 1000,
                                       range.x   = rr
                                       )
                        fdensD <- bkde(x         = fdmat[100:(its * 1000), dd],
                                       bandwidth = dpik(fdmat[100:(its * 1000), dd], gridsize = 1000, range.x = rr),
                                       gridsize  = 1000,
                                       range.x   = rr
                                       )
                        accDmat[nn, kk, ff, cc, dd, its] <- max(1 - sum(abs(pdensD$y - fdensD$y) * diff(pdensD$x)[1]) / 2, 0)
                    }
                }

                rr <- range(c(psig[100:20000], fsig[100:20000]))
                for (its in 1:20) {
                    pdensS <- bkde(x         = psig[100:(its * 1000)],
                                   bandwidth = dpik(psig[100:(its * 1000)], gridsize = 1000, range.x = rr),
                                   gridsize  = 1000,
                                   range.x   = rr
                                   )
                    fdensS <- bkde(x         = fsig[100:(its * 1000)],
                                   bandwidth = dpik(fsig[100:(its * 1000)], gridsize = 1000, range.x = rr),
                                   gridsize  = 1000,
                                   range.x   = rr
                                   )
                    accSig[nn, kk, ff, cc, its] <- max(1 - sum(abs(pdensS$y - fdensS$y) * diff(pdensS$x)[1]) / 2, 0)

                }

            }
        }
    }
}

accBetaEst <- apply(accBeta, c(1:3, 6), mean, na.rm = TRUE)
accDmatEst <- apply(accDmat, c(1:3, 6), mean, na.rm = TRUE)
accSigEst <- apply(accSig, c(1:3, 5), mean, na.rm = TRUE)

accBetaEstDf <- adply(accBetaEst, 1:4)
colnames(accBetaEstDf) <- c("nsamp", "nsub", "frac", "niter", "est")
accDmatEstDf <- adply(accDmatEst, 1:4)
colnames(accDmatEstDf) <- c("nsamp", "nsub", "frac", "niter", "est")

accSigEstDf <- adply(accSigEst, 1:4)
colnames(accSigEstDf) <- c("nsamp", "nsub", "frac", "niter", "est")


accBetaEstDf$niter <- as.numeric(gsub("iter", "", as.character(accBetaEstDf$niter)))
accDmatEstDf$niter <- as.numeric(gsub("iter", "", as.character(accDmatEstDf$niter)))
accSigEstDf$niter <- as.numeric(gsub("iter", "", as.character(accSigEstDf$niter)))

plt <- ggplot(accBetaEstDf, aes(x = niter, y = est, group = frac, colour = frac)) +
    xlab("Iterations (x 1000)") + ylab("Accuracy for" ~ beta) +
    scale_colour_brewer(name="",
                        breaks=c("f1", "f2", "f3", "f4"),
                        labels=c("r = 0.20", "r = 0.40", "r = 0.60", "r = 0.80"),
                        palette="Set1",
                        guide = guide_legend(override.aes = list(shape = NA))) +
    facet_grid(nsub  ~ nsamp,
               labeller = labeller(.rows = c(k1 = "k = 10", k2 = "k = 25", k3 = "k = 50"),
                                   .cols = c(n1 = "n = 10,000", n2 = "n = 100,000", n3 = "n = 1,000,000")
                                   )
               ) +
    geom_line(stat = "identity", size = 2) + geom_point(size = 8) +
    theme_bw() +
    theme(strip.text.x = element_text(size=24, face="bold"),
          strip.text.y = element_text(size=24, face="bold"),
          axis.text.x = element_text(size=20, face="bold"),
          axis.text.y = element_text(size=20, face="bold"),
          axis.title.x = element_text(size=24, face="bold"),
          axis.title.y = element_text(size=24, face="bold"),
          legend.text = element_text(size = 24)
          ) +
    ylim(0.90, 1.00)


plt1 <- ggplot(accDmatEstDf, aes(x = niter, y = est, group = frac, colour = frac)) +
    xlab("Iterations (x 1000)") + ylab("Accuracy for" ~ Sigma) +
    scale_colour_brewer(name="",
                        breaks=c("f1", "f2", "f3", "f4"),
                        labels=c("r = 0.20", "r = 0.40", "r = 0.60", "r = 0.80"),
                        palette="Set1",
                        guide = guide_legend(override.aes = list(shape = NA))) +
    facet_grid(nsub  ~ nsamp,
               labeller = labeller(.rows = c(k1 = "k = 10", k2 = "k = 25", k3 = "k = 50"),
                                   .cols = c(n1 = "n = 10,000", n2 = "n = 100,000", n3 = "n = 1,000,000")
                                   )
               ) +
    geom_line(stat = "identity", size = 2) + geom_point(size = 8) +
    theme_bw() +
    theme(strip.text.x = element_text(size=24, face="bold"),
          strip.text.y = element_text(size=24, face="bold"),
          axis.text.x = element_text(size=20, face="bold"),
          axis.text.y = element_text(size=20, face="bold"),
          axis.title.x = element_text(size=24, face="bold"),
          axis.title.y = element_text(size=24, face="bold"),
          legend.text = element_text(size = 24)
          ) +
    ylim(0.90, 1.00)

plt2 <- ggplot(accSigEstDf, aes(x = niter, y = est, group = frac, colour = frac)) +
    xlab("Iterations (x 1000)") + ylab(expression("Accuracy for" ~ sigma^2)) +
    scale_colour_brewer(name="",
                        breaks=c("f1", "f2", "f3", "f4"),
                        labels=c("r = 0.20", "r = 0.40", "r = 0.60", "r = 0.80"),
                        palette="Set1",
                        guide = guide_legend(override.aes = list(shape = NA))) +
    geom_line(stat = "identity", size = 2) + geom_point(size = 8) +
    facet_grid(nsub  ~ nsamp,
               labeller = labeller(.rows = c(k1 = "k = 10", k2 = "k = 25", k3 = "k = 50"),
                                   .cols = c(n1 = "n = 10,000", n2 = "n = 100,000", n3 = "n = 1,000,000")
                                   )
               ) +
    theme_bw() +
    theme(strip.text.x = element_text(size=24, face="bold"),
          strip.text.y = element_text(size=24, face="bold"),
          axis.text.x = element_text(size=20, face="bold"),
          axis.text.y = element_text(size=20, face="bold"),
          axis.title.x = element_text(size=24, face="bold"),
          axis.title.y = element_text(size=24, face="bold"),
          legend.text = element_text(size = 24)
          ) +
    ylim(0.90, 1.00)


pdf("lme_beta.pdf", 12, 9)
print(plt)
dev.off()
pdf("lme_dmat.pdf", 12, 9)
print(plt1)
dev.off()
pdf("lme_sig.pdf", 12, 9)
print(plt2)
dev.off()

#### MCSE COMPUTATIONS

setwd("/Shared/ssrivastva/dist_lme/result/sim/")
rm(list = ls())

nfrac <- c(20, 40, 60, 80)

mcmcseBeta <- array(NA, dim = c(3, 3, 4, 10, 4, 20),
                 dimnames = list(paste0("n", 1:3),
                                 paste0("k", 1:3),
                                 paste0("f", 1:4),
                                 paste0("cv", 1:10),
                                 paste0("dim", 1:4),
                                 paste0("iter", 1:20)
                                 ))

mcmcseDmat <- array(NA, dim = c(3, 3, 4, 10, 6, 20),
                 dimnames = list(paste0("n", 1:3),
                                 paste0("k", 1:3),
                                 paste0("f", 1:4),
                                 paste0("cv", 1:10),
                                 paste0("dim", 1:6),
                                 paste0("iter", 1:20)
                                 ))

mcmcseSig <- array(NA, dim = c(3, 3, 4, 10, 20),
                 dimnames = list(paste0("n", 1:3),
                                 paste0("k", 1:3),
                                 paste0("f", 1:4),
                                 paste0("cv", 1:10),
                                 paste0("iter", 1:20)
                                 ))

for (nn in 1:3) {
    for (kk in 1:3) {
        for (ff in 1:4) {
            for (cc in 1:10) {
                cat("nn: ", nn, " kk: ", kk, " ff: ", ff, " cc: ", cc, "\n")
                fname <- paste0("part/px_sim_part_n_", nn, "_k_", kk, "_frac_", nfrac[ff], "_rep_", cc, ".rds")

                res <- readRDS(fname)
                ffname <- paste0("full/sim_full_n_", nn, "_rep_", cc, ".rds")

                fres <- readRDS(ffname)
                pbetas <- res$betas
                fbetas <- fres$betas
                pdmat <- res$tmat
                fdmat <- cbind(do.call(rbind, lapply(fres$tmat, function(x) diag(x))),
                               do.call(rbind, lapply(fres$tmat, function(x) x[lower.tri(x)]))
                               )
                psig <- res$sigma
                fsig <- fres$sigma

                for (its in 1:20) {
                    mcmcseBeta[nn, kk, ff, cc, , its] <- mcse.mat(pbetas[100:(its * 1000), ], method = "obm")[ , 2] - mcse.mat(fbetas[100:(its * 1000), ], method = "obm")[ , 2]
                }
                for (its in 1:20) {
                    mcmcseDmat[nn, kk, ff, cc, , its] <- mcse.mat(pdmat[100:(its * 1000), ], method = "obm")[ , 2] - mcse.mat(fdmat[100:(its * 1000), ], method = "obm")[ , 2]
                }
                for (its in 1:20) {
                    mcmcseSig[nn, kk, ff, cc, its] <- mcse(psig[100:(its * 1000)], method = "obm")$se - mcse(fsig[10:(its * 1000)], method = "obm")$se
                }
            }
        }
    }
}


saveRDS(mcmcseBeta, "~/dist_lme/data/mcse_betas.rds")
saveRDS(mcmcseDmat, "~/dist_lme/data/mcse_dmat.rds")
saveRDS(mcmcseSig, "~/dist_lme/data/mcse_sig.rds")

betaEst <- apply(abs(mcmcseBeta), c(1:3, 6), mean, na.rm = TRUE)
dmatEst <- apply(abs(mcmcseDmat), c(1:3, 6), mean, na.rm = TRUE)
sigEst <- apply(abs(mcmcseSig), c(1:3, 5), mean, na.rm = TRUE)

betaEstDf <- adply(betaEst, 1:4)
dmatEstDf <- adply(dmatEst, 1:4)
sigEstDf <- adply(sigEst, 1:4)

colnames(betaEstDf) <- colnames(dmatEstDf) <- colnames(sigEstDf) <- c("nsamp", "nsub", "frac", "niter", "est")

betaEstDf$niter <- as.numeric(gsub("iter", "", as.character(betaEstDf$niter)))
dmatEstDf$niter <- as.numeric(gsub("iter", "", as.character(dmatEstDf$niter)))
sigEstDf$niter <- as.numeric(gsub("iter", "", as.character(sigEstDf$niter)))

plt <-
ggplot(betaEstDf, aes(x = niter, y = est, group = frac, colour = frac)) +
    xlab("Iterations (x 1000)") + ylab("SE Difference for" ~ bold(beta)) +
    scale_colour_brewer(name="",
                        breaks=c("f1", "f2", "f3", "f4"),
                        labels=c("r = 0.20", "r = 0.40", "r = 0.60", "r = 0.80"),
                        palette="Set1",
                        guide = guide_legend(override.aes = list(shape = NA))) +
    facet_grid(nsub  ~ nsamp,
               labeller = labeller(.rows = c(k1 = "k = 10", k2 = "k = 25", k3 = "k = 50"),
                                   .cols = c(n1 = "n = 10,000", n2 = "n = 100,000", n3 = "n = 1,000,000")
                                   )
               ) +
    geom_line(stat = "identity", size = 4) + geom_point(size = 8) +
    theme_bw() +
    theme(strip.text.x = element_text(size=24, face="bold"),
          strip.text.y = element_text(size=24, face="bold"),
          axis.text.x = element_text(size=20, face="bold"),
          axis.text.y = element_text(size=20, face="bold"),
          axis.title.x = element_text(size=24, face="bold"),
          axis.title.y = element_text(size=24, face="bold"),
          legend.text = element_text(size = 24)
          ) + ylim(0, 0.001)


plt1 <-
ggplot(dmatEstDf, aes(x = niter, y = est, group = frac, colour = frac)) +
    xlab("Iterations (x 1000)") + ylab("SE Difference for" ~ bold(Sigma)) +
    scale_colour_brewer(name="",
                        breaks=c("f1", "f2", "f3", "f4"),
                        labels=c("r = 0.20", "r = 0.40", "r = 0.60", "r = 0.80"),
                        palette="Set1",
                        guide = guide_legend(override.aes = list(shape = NA))) +
    facet_grid(nsub  ~ nsamp,
               labeller = labeller(.rows = c(k1 = "k = 10", k2 = "k = 25", k3 = "k = 50"),
                                   .cols = c(n1 = "n = 10,000", n2 = "n = 100,000", n3 = "n = 1,000,000")
                                   )
               ) +
    geom_line(stat = "identity", size = 4) + geom_point(size = 8) +
    theme_bw() +
    theme(strip.text.x = element_text(size=24, face="bold"),
          strip.text.y = element_text(size=24, face="bold"),
          axis.text.x = element_text(size=20, face="bold"),
          axis.text.y = element_text(size=20, face="bold"),
          axis.title.x = element_text(size=24, face="bold"),
          axis.title.y = element_text(size=24, face="bold"),
          legend.text = element_text(size = 24)
          )+ ylim(0, 0.001)


plt2 <-
ggplot(sigEstDf, aes(x = niter, y = est, group = frac, colour = frac)) +
    xlab("Iterations (x 1000)") + ylab("SE Difference for" ~ sigma^2) +
    scale_colour_brewer(name="",
                        breaks=c("f1", "f2", "f3", "f4"),
                        labels=c("r = 0.20", "r = 0.40", "r = 0.60", "r = 0.80"),
                        palette="Set1",
                        guide = guide_legend(override.aes = list(shape = NA))) +
    facet_grid(nsub  ~ nsamp,
               labeller = labeller(.rows = c(k1 = "k = 10", k2 = "k = 25", k3 = "k = 50"),
                                   .cols = c(n1 = "n = 10,000", n2 = "n = 100,000", n3 = "n = 1,000,000")
                                   )
               ) +
    geom_line(stat = "identity", size = 4) + geom_point(size = 8) +
    theme_bw() +
    theme(strip.text.x = element_text(size=24, face="bold"),
          strip.text.y = element_text(size=24, face="bold"),
          axis.text.x = element_text(size=20, face="bold"),
          axis.text.y = element_text(size=20, face="bold"),
          axis.title.x = element_text(size=24, face="bold"),
          axis.title.y = element_text(size=24, face="bold"),
          legend.text = element_text(size = 24)
          )+ ylim(0, 0.001)

pdf("~/dist_lme/data/mcse_lme_beta.pdf", 12, 9)
print(plt)
dev.off()
pdf("~/dist_lme/data/mcse_lme_dmat.pdf", 12, 9)
print(plt1)
dev.off()
pdf("~/dist_lme/data/mcse_lme_sig.pdf", 12, 9)
print(plt2)
dev.off()
