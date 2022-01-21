### Title:    Demonstrate the Simulation Routines
### Author:   Kyle M. Lang
### Created:  2021-03-01
### Modified: 2022-01-21

rm(list = ls(all = TRUE))

source("simMissingness.R")
source("simData.R")

library(mvtnorm) # Need this to simluate data via hardEtAl2012DataSim()
library(pROC)    # Need this to calculate AUC in missingness simulations

n   <- 100
p   <- 5
pm  <- 0.25
auc <- 0.8
snr <- 0.5

## Generate some synthetic data:
X <- as.data.frame(matrix(rnorm(n * p), ncol = p))
y <- rnorm(n)

###--------------------------------------------------------------------------###

## Simulate simple missingness via logistic regression
## - Proportion missing = pm
## - Missing values in the positive tail of the linear predictor
out <- simLogisticMissingness0(pm      = pm,
                               data    = X,
                               type    = "high",
                               stdData = TRUE,
                               beta    = runif(5)
                               )

## Print the achieved AUC:
out$fit

## Impose missing data on y:
y2        <- y
y2[out$r] <- NA

## Compute proportion missing:
mean(is.na(y2))

## NOTE:
## This approach does not control the strength of the MAR assocation. Unlike
## the methods demonstrated below, this appraoch is computationally stable at
## small sample sizes and produces asymptotically accurate proportions of
## missingness, but the strength of the association varies as an unspecified
## function of the predictors and the beta weights.

###--------------------------------------------------------------------------###

## Simulate missingness via logistic regression
## - Proportion missing = pm
## - AUC for predicting the missingness = auc
## - Missing values in the positive tail of the linear predictor
## - AUC achieved via optimizing the slope values
out <- simLogisticMissingness(pm       = pm,
                              auc      = auc,
                              data     = X,
                              type     = "high",
                              optimize = "slopes",
                              stdData  = TRUE,
                              beta     = runif(5)
                              )

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
## - AUC for predicting the missingness = auc
## - Missing values in the negative tail of the linear predictor
## - AUC achieved via optimizing the SNR
out <- simLinearMissingness(pm       = pm,
                            auc      = auc,
                            data     = X,
                            preds    = paste0("V", 1 : 3),
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

## Simulate missingness via a linear probability model with fixed SNR
## - Proportion missing = pm
## - SNR for the linear probability model = snr
## - Missing values in the upper tail of the linear predictor
## - Only one MAR predictor
out <- simLinearMissingness(pm       = pm,
                            snr      = snr,
                            data     = X,
                            preds    = "V2",
                            type     = "high",
                            optimize = FALSE)

## Print the achieved AUC:
out$auc

## Print the achieved SNR:
out$snr

## Impose missing data on y:
y2        <- y
y2[out$r] <- NA

## Compute proportion missing:
mean(is.na(y2))

###--------------------------------------------------------------------------###

## Simulate data that matches Hardt et al (2012)
## - Sample size = n
## - Correlations involving auxiliary variables = rZ
## - Number of auxiliary variables = pZ

dat1 <- hardtEtAl2012DataSim(n = 10000, rZ = 0.5, pZ = 5)

## Check the results:
head(dat1)

summary(lm(y ~ x1 + x2, data = dat1))

cor(dat1)
