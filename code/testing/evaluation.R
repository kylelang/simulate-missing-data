### Title:    Simulation to Evaluate Missingness Generation
### Author:   Kyle M. Lang
### Created:  2021-03-04
### Modified: 2021-03-04

rm(list = ls(all = TRUE))

library(mvtnorm)
library(pROC)
library(parallel)

source("../simMissingness.R")
source("subroutines.R")

set.seed(235711)

nReps <- 250

n    <- c(500, 250, 100, 50)
p    <- c(25, 5, 1)
pm   <- c(0.1, 0.25, 0.5)
auc  <- seq(0.55, 0.95, 0.2)
type <- c("high", "low", "center", "tails")
dist <- c("norm", "unif", "gamma")

conds <- expand.grid(n = n, p = p, pm = pm, auc = auc, type = type, dist = dist,
                     stringsAsFactors = FALSE)

## Run simulation in parallel:
out <- mclapply(1 : nReps, doRep, conds = conds, mc.cores = 4)

## Time one full replication:
                                        #t1 <- system.time(
                                        #    out <- doRep(1, conds)
                                        #)
