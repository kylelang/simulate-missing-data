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
pm   <- 0.5
auc  <- 0.75
beta <- rep(10, p) #runif(p, -1, 1)
type <- "high"

                                        #X <- as.data.frame(rmvnorm(n, rep(0, p), diag(p)))
X <- as.data.frame(matrix(runif(n * p), ncol = p))

tmp <- simMissingness(pm     = pm,
                      auc    = auc,
                      data   = X,
                      type   = type,
                      beta   = beta,
                      stdEta = FALSE,
                      smooth = FALSE)

r <- tmp$r
mean(r)

tmp$auc

plot(x = tmp$eta, y = tmp$p)

fit <- glm(r ~ tmp$eta, family = binomial)

rocOut <- roc(r, predict(fit), smooth = TRUE)
plot(rocOut)
auc(rocOut)
