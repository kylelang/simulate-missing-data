### Title:    Demonstrate the Missingness Simulation Routines
### Author:   Kyle M. Lang
### Created:  2021-03-01
### Modified: 2021-03-10

rm(list = ls(all = TRUE))

source("simMissingness.R")

n   <- 1000
p   <- 5
pm  <- 0.25
auc <- 0.8
snr <- 0.5

## Generate some synthetic data:
X <- as.data.frame(matrix(rnorm(n * p), ncol = p))
y <- rnorm(n)

###--------------------------------------------------------------------------###

## Simulate missingness via logistic regression
## - Proportion missing = pm
## - AUC for predicting the missingess = auc
## - Missing values in the positive tail of the linear predictor
## - AUC achieved via optimizing the slope values
out <- simLogisticMissingness(pm       = pm,
                              auc      = auc,
                              data     = X,
                              type     = "high",
                              optimize = "slopes")

## Print the achieved AUC:
out$fit

## Impose missing data on y:
y2        <- y
y2[out$r] <- NA

## Compute proportion missing:
mean(is.na(y2))

## WARNING:
## This approach is highly unstable at small sample sizes (i.e., N << 1000).
## In small samples, the optimization can fail or results (in terms of PM and
## AUC) can be negatively biased.

###--------------------------------------------------------------------------###

## Simulate missingness via a linear probability model
## - Proportion missing = pm
## - AUC for predicting the missingess = auc
## - Missing values in the negative tail of the linear predictor
## - AUC achieved via optimizing the SNR
out <- simLinearMissingness(pm       = pm,
                            auc      = auc,
                            data     = X,
                            type     = "low",
                            optimize = TRUE)


## Print the achieved AUC:
out$auc

## Print the achieved SNR:
out$snr

## Impose missing data on y:
y2        <- y
y2[out$r] <- NA

## Compute proportion missing:
mean(is.na(y2))

## WARNING:
## This approach is highly unstable at small sample sizes (i.e., N << 1000).
## In small samples, the optimization can fail or results (in terms of PM and
## AUC) can be negatively biased.

###--------------------------------------------------------------------------###

## Simulate missingness via a linear probability model with fixed SNR
## - Proportion missing = pm
## - SNR for the linear probability model = snr
## - Missing values in both tails of the linear predictor
out <- simLinearMissingness(pm       = pm,
                            snr      = snr,
                            data     = X,
                            type     = "tails",
                            optimize = FALSE)


## AUC is not calculated when simulating symmetric missingness via the linear
## probability model:
out$auc

## Print the achieved SNR:
out$snr

## Impose missing data on y:
y2        <- y
y2[out$r] <- NA

## Compute proportion missing:
mean(is.na(y2))

###--------------------------------------------------------------------------###

