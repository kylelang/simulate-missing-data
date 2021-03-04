### Title:    Testing Missingness Simulation Routines
### Author:   Kyle M. Lang
### Created:  2021-03-01
### Modified: 2021-03-04

rm(list = ls(all = TRUE))

library(mvtnorm)
library(pROC)

source("simMissingness.R")

n    <- 1000
p    <- 1
pm   <- 0.5
auc  <- 0.65
beta <- 1 #c(-1, 2, 3, -4)
type <- "tails"

                                        #X <- as.data.frame(rmvnorm(n, rep(0, p), diag(p)))
X <- as.data.frame(matrix(runif(n * p), ncol = p))

tmp <- simMissingness(pm = pm, auc = auc, data = X, type = type, beta = beta)

r <- tmp$r
mean(r)

tmp$auc

plot(x = tmp$eta, y = tmp$p)

fit <- glm(r ~ abs(tmp$eta), family = binomial)

rocOut <- roc(r, abs(predict(fit)), smooth = TRUE)
plot(rocOut)
auc(rocOut)
