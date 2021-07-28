setwd("/Shared/ssrivastva/dist_logistic/result/sim/")
rm(list = ls())


library(KernSmooth)
library(mcmcse)
library(plyr)
library(ggplot2)

nfrac <- c(20, 40, 60, 80)
## ACCURACY COMPUTATIONS

acc <- array(NA, dim = c(3, 3, 4, 10, 10, 20),
             dimnames = list(paste0("n", 1:3),
                             paste0("k", 1:3),
                             paste0("f", 1:4),
                             paste0("cv", 1:10),
                             paste0("dim", 1:10),
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
                cat("kk: ", kk, " ff: ", ff, " cc: ", cc, "\n")
                fname <- paste0("part/sim_part_n_", nn, "_k_", kk, "_frac_", nfrac[ff], "_rep_", cc, ".rds")
                res <- readRDS(fname)
                ffname <- paste0("full/sim_full_n_", nn, "_rep_", cc, ".rds")
                fres <- readRDS(ffname)
                pbetas <- res$betas
                fbetas <- fres$betas
                times[nn, kk, ff, cc, "full"] <- fres$time
                times[nn, kk, ff, cc, "part"] <- res$time
                for (dd in 1:10) {
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
                        acc[nn, kk, ff, cc, dd, its] <- max(1 - sum(abs(pdens$y - fdens$y) * diff(pdens$x)[1]) / 2, 0)
                    }
                }
            }
        }
    }
}

accEst <- apply(acc, c(1:3, 6), mean)

library(plyr)
library(ggplot2)

accEstDf <- adply(accEst, 1:4)
colnames(accEstDf) <- c("nsamp", "nsub", "frac", "niter", "est")

accEstDf$niter <- as.numeric(gsub("iter", "", as.character(accEstDf$niter)))

plt <- ggplot(accEstDf, aes(x = niter, y = est, group = frac, colour = frac)) +
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
    geom_line(stat = "identity", size = 4) + geom_point(size = 8) +
    theme_bw() +
    theme(strip.text.x = element_text(size=24, face="bold"),
          strip.text.y = element_text(size=24, face="bold"),
          axis.text.x = element_text(size=20, face="bold"),
          axis.text.y = element_text(size=20, face="bold"),
          axis.title.x = element_text(size=24, face="bold"),
          axis.title.y = element_text(size=24, face="bold"),
          legend.text = element_text(size = 24)
          ) + ylim(0.5, 1.00)

pdf("log_beta.pdf", 12, 9)
print(plt)
dev.off()

### FUNCTION of BETA

acc <- array(NA, dim = c(3, 3, 4, 10, 20),
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
                cat("kk: ", kk, " ff: ", ff, " cc: ", cc, "\n")
                fname <- paste0("part/sim_part_n_", nn, "_k_", kk, "_frac_", nfrac[ff], "_rep_", cc, ".rds")
                res <- readRDS(fname)
                ffname <- paste0("full/sim_full_n_", nn, "_rep_", cc, ".rds")
                fres <- readRDS(ffname)
                pbetas <- res$betas
                fbetas <- fres$betas

                pprobs <- 1 / (1 + exp(-(drop(pbetas %*% rep(1, 10)))))
                fprobs <- 1 / (1 + exp(-(drop(fbetas %*% rep(1, 10)))))

                rr <- range(c(pprobs[100:20000], fprobs[100:20000]))
                for (its in 1:20) {
                    pdens <- bkde(x         = pprobs[1:(its * 1000)],
                                  bandwidth = dpik(pprobs[1:(its * 1000)], gridsize = 1000, range.x = rr),
                                  gridsize = 1000,
                                  range.x   = rr
                                  )
                    fdens <- bkde(x         = fprobs[1:(its * 1000)],
                                  bandwidth = dpik(fprobs[1:(its * 1000)], gridsize = 1000, range.x = rr),
                                  gridsize = 1000,
                                  range.x   = rr
                                  )
                    acc[nn, kk, ff, cc, its] <- max(1 - sum(abs(pdens$y - fdens$y) * diff(pdens$x)[1]) / 2, 0)
                }
            }
        }
    }
}


accProbEst <- apply(acc, c(1:3, 5), mean, na.rm = TRUE)
accProbEstDf <- adply(accProbEst, 1:4)
colnames(accProbEstDf) <- c("nsamp", "nsub", "frac", "niter", "est")
accProbEstDf$niter <- as.numeric(gsub("iter", "", as.character(accProbEstDf$niter)))

plt1 <-
ggplot(accProbEstDf, aes(x = niter, y = est, group = frac, colour = frac)) +
    xlab("Iterations (x 1000)") + ylab("Accuracy for P(Y = 1 | X)") +
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
          ) + ylim(0.5, 1.00)


pdf("~/dist_logistic/data/log_prob.pdf", 12, 9)
print(plt1)
dev.off()


## MCSE COMPUTATION
setwd("/Shared/ssrivastva/dist_logistic/result/sim/")
rm(list = ls())

nfrac <- c(20, 40, 60, 80)

mcmcse <- array(NA, dim = c(3, 3, 4, 10, 10, 20),
             dimnames = list(paste0("n", 1:3),
                             paste0("k", 1:3),
                             paste0("f", 1:4),
                             paste0("cv", 1:10),
                             paste0("dim", 1:10),
                             paste0("iter", 1:20)
                             ))

for (nn in 1:3) {
    for (kk in 1:3) {
        for (ff in 1:4) {
            for (cc in 1:10) {
                cat("kk: ", kk, " ff: ", ff, " cc: ", cc, "\n")
                fname <- paste0("part/sim_part_n_", nn, "_k_", kk, "_frac_", nfrac[ff], "_rep_", cc, ".rds")
                res <- readRDS(fname)
                ffname <- paste0("full/sim_full_n_", nn, "_rep_", cc, ".rds")
                fres <- readRDS(ffname)
                pbetas <- res$betas
                fbetas <- fres$betas
                for (its in 1:20) {
                    mcmcse[nn, kk, ff, cc, , its] <- mcse.mat(pbetas[100:(its * 1000), ], method = "obm")[ , 2] - mcse.mat(fbetas[100:(its * 1000), ], method = "obm")[ , 2]
                }
            }
        }
    }
}

saveRDS(mcmcse, "~/dist_logistic/data/mcse_betas.rds")

mcmcse <- array(NA, dim = c(3, 3, 4, 10, 20),
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
                cat("kk: ", kk, " ff: ", ff, " cc: ", cc, "\n")
                fname <- paste0("part/sim_part_n_", nn, "_k_", kk, "_frac_", nfrac[ff], "_rep_", cc, ".rds")
                res <- readRDS(fname)
                ffname <- paste0("full/sim_full_n_", nn, "_rep_", cc, ".rds")
                fres <- readRDS(ffname)
                pbetas <- res$betas
                fbetas <- fres$betas

                pprobs <- 1 / (1 + exp(-(drop(pbetas %*% rep(1, 10)))))
                fprobs <- 1 / (1 + exp(-(drop(fbetas %*% rep(1, 10)))))

                for (its in 1:20) {
                    mcmcse[nn, kk, ff, cc, its] <- mcse(pprobs[100:(its * 1000)], method = "obm")$se - mcse(fbetas[100:(its * 1000)], method = "obm")$se
                }
            }
        }
    }
}

saveRDS(mcmcse, "~/dist_logistic/data/mcmcse_probs.rds")

#### MCMCSE
setwd("~/dist_logistic/code")
rm(list=ls())

mcmcseBeta <- readRDS("~/dist_logistic/data/mcse_betas.rds")
mcmcseProbs <- readRDS("~/dist_logistic/data/mcmcse_probs.rds")

betaEst <- apply(abs(mcmcseBeta), c(1:3, 6), mean, na.rm = TRUE)
probEst <- apply(abs(mcmcseProbs), c(1:3, 5), mean, na.rm = TRUE)

betaEstDf <- adply(betaEst, 1:4)
probEstDf <- adply(probEst, 1:4)

colnames(betaEstDf) <- colnames(probEstDf) <- c("nsamp", "nsub", "frac", "niter", "est")

betaEstDf$niter <- as.numeric(gsub("iter", "", as.character(betaEstDf$niter)))
probEstDf$niter <- as.numeric(gsub("iter", "", as.character(probEstDf$niter)))

plt <-
ggplot(betaEstDf, aes(x = niter, y = est, group = frac, colour = frac)) +
    xlab("Iterations (x 1000)") + ylab("SE Difference for" ~ beta) +
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
          ) + ylim (0.00, 0.015)


plt1 <-
ggplot(probEstDf, aes(x = niter, y = est, group = frac, colour = frac)) +
    xlab("Iterations (x 1000)") + ylab("SE Difference for P(Y = 1 | X)") +
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
          )+ ylim (0.00, 0.015)


pdf("~/dist_logistic/data/mcse_log_beta.pdf", 12, 9)
print(plt)
dev.off()
pdf("~/dist_logistic/data/mcse_log_prob.pdf", 12, 9)
print(plt1)
dev.off()
