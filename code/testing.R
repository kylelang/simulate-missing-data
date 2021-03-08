### Title:    Testing Missingness Simulation Routines
### Author:   Kyle M. Lang
### Created:  2021-03-01
### Modified: 2021-03-04

rm(list = ls(all = TRUE))

library(mvtnorm)
library(pROC)

source("simMissingness.R")

n    <- 1000
p    <- 5
pm   <- 0.25
auc  <- 0.75
snr  <- 0.5
beta <- rep(10, p) #runif(p, -1, 1)
type <- "tails"
opt  <- "slopes"

X <- as.data.frame(rmvnorm(n, rep(0, p), diag(p)))
                                        #X <- as.data.frame(matrix(runif(n * p), ncol = p))

tmp <- simMissingness(pm       = pm,
                      auc      = auc,
                      snr      = snr,
                      data     = X,
                      type     = type,
                      beta     = beta,
                      stdEta   = TRUE,
                      optimize = opt,
                      smooth   = FALSE)

r <- tmp$r
mean(r)

tmp$fit

plot(x = abs(tmp$eta), y = tmp$p)

fit <- glm(r ~ abs(tmp$eta), family = binomial)

rocOut <- roc(r, predict(fit), smooth = TRUE)
plot(rocOut)
auc(rocOut)
