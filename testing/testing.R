### Title:    Testing Missingness Simulation Routines
### Author:   Kyle M. Lang
### Created:  2021-03-01
### Modified: 2021-03-04

rm(list = ls(all = TRUE))

library(mvtnorm)
library(pROC)

source("../code/simMissingness.R")

n    <- 10000
p    <- 5
pm   <- 0.25
auc  <- 0.75
snr  <- NULL
beta <- runif(p, -1, 1)
type <- "low"
opt  <- "slopes"
std  <- TRUE

X <- rmvnorm(n, rep(10, p), diag(p))
                                        #X <- as.data.frame(matrix(runif(n * p), ncol = p))

tmp <- simLogisticMissingness(pm       = pm,
                              auc      = auc,
                              snr      = snr,
                              data     = X,
                              type     = type,
                              beta     = beta,
                              stdData  = TRUE,
                              optimize = opt,
                              smooth   = FALSE)

tmp <- simLinearMissingness(pm       = pm,
                            auc      = auc,
                            snr      = snr,
                            data     = X,
                            type     = type,
                            beta     = beta,
                            optimize = opt,
                            stdData  = std,
                            smooth   = FALSE)



r <- tmp$r
mean(r)

tmp$auc
tmp$snr
tmp$fit

plot(x = abs(tmp$eta), y = tmp$p)

fit <- glm(r ~ abs(tmp$eta), family = binomial)

n    <- 100000
p    <- 5
pm   <- 0.25
snr  <- 0.25
beta <- rep(5, p) #runif(p, -1, 1)

X <- rmvnorm(n, rep(0, p), diag(p))
     
eta  <- X %*% matrix(beta)

                                        #s2   <- (1 / snr) * var(eta)
                                        #eta2 <- eta + rnorm(length(eta), 0, sqrt(s2))

eta2 <- eta + (1 / sqrt(snr)) * rnorm(length(eta), 0, sd(eta))

r <- eta2 < quantile(eta2, pm)

mean(r)

rocOut <- roc(as.numeric(r), as.numeric(eta))
plot(rocOut)
auc(rocOut)

var(eta) / (var(eta2) - var(eta))

n   <- 1000000
s   <- 5
snr <- 0.5

sig   <- rnorm(n, 0, s)
noise <- (1 / snr) * rnorm(n, 0, s)

x <- sig + noise

var(sig) / var(noise)
sd(sig) / sd(noise)

var(sig) / (var(x) - var(sig))
sd(sig) / sqrt(var(x) - var(sig))
