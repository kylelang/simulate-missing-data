### Title:    Simulation to Evaluate Missingness Generation
### Author:   Kyle M. Lang
### Created:  2021-03-04
### Modified: 2021-03-08

rm(list = ls(all = TRUE))

library(mvtnorm)
library(pROC)
library(parallel)

source("../code/simMissingness.R")
source("subroutines.R")

outDir <- "../output/optSlopes/"

set.seed(235711)

## expNum:
## 1 = Logistic model & optimized slopes
## 2 = Logistic model & fixed slopes
## 3 = Logistic model & optimized noise
## 4 = Linear probability model & optimized noise
## 5 = Linear probability model & fixed SNR

expNum <- 5

nReps <- 250

                                        #auc  <- seq(0.55, 0.95, 0.2)
snr  <- c(0.5, 1.0, 1.5)
n    <- c(500, 250, 100, 50)
p    <- c(25, 5, 1)
pm   <- c(0.1, 0.25, 0.5)
type <- c("high", "low") #, "center", "tails")
dist <- c("norm", "unif", "gamma")

conds <- expand.grid(n    = n,
                     p    = p,
                     pm   = pm,
                     snr  = snr,
                                        #auc  = auc,
                     type = type,
                     dist = dist,
                     stringsAsFactors = FALSE)

saveRDS(conds, paste0(outDir, "conds.rds"))

## Run simulation in parallel:
out <- mclapply(1 : nReps,
                FUN      = doRep,
                conds    = conds,
                expNum   = expNum,
                smooth   = FALSE,
                mc.cores = 4)

saveRDS(out, paste0(outDir, "evaluation_output.rds"))

## Time one full replication:
t1 <- system.time(
    out <- doRep(1, conds, expNum = expNum)
)

