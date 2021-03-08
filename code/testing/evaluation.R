### Title:    Simulation to Evaluate Missingness Generation
### Author:   Kyle M. Lang
### Created:  2021-03-04
### Modified: 2021-03-08

rm(list = ls(all = TRUE))

library(mvtnorm)
library(pROC)
library(parallel)

source("../simMissingness.R")
source("subroutines.R")

outDir <- "../../output/noise/"

set.seed(235711)

## expNum:
## 1 = Optimized slopes
## 2 = Fixed slopes
## 3 = Optimized noise
expNum <- 3

nReps <- 250

auc  <- seq(0.55, 0.95, 0.2)
n    <- c(500, 250, 100, 50)
p    <- c(25, 5, 1)
pm   <- c(0.1, 0.25, 0.5)
type <- c("high", "low", "center", "tails")
dist <- c("norm", "unif", "gamma")

conds <- expand.grid(n    = n,
                     p    = p,
                     pm   = pm,
                     auc  = auc,
                     type = type,
                     dist = dist,
                     stringsAsFactors = FALSE)

saveRDS(conds, paste0(outDir, "conds.rds"))

## Run simulation in parallel:
out <- mclapply(1 : nReps,
                FUN         = doRep,
                conds       = conds,
                fixedSlopes = FALSE,
                smooth      = FALSE,
                mc.cores    = 4)

saveRDS(out, paste0(outDir, "evaluation_output.rds"))

## Time one full replication:
t1 <- system.time(
    out <- doRep(1, conds, expNum = expNum)
)

125 * t1 / 60^2
