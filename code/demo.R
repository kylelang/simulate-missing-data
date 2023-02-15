### Title:    Demonstrate the Simulation Routines
### Author:   Kyle M. Lang
### Created:  2021-03-01
### Modified: 2022-01-21

rm(list = ls(all = TRUE))

source("simMissingness.R")
source("simData.R")

library(mvtnorm) # Need this to simluate data via hardEtAl2012DataSim()
library(pROC)    # Need this to calculate AUC in missingness simulations

n   <- 5000
p   <- 5
pm  <- 0.25
auc <- 0.8
snr <- 0.5
r2  <- 0.5

## Generate some synthetic data:
X <- as.data.frame(matrix(rnorm(n * p), ncol = p))

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

undebug(simLogisticMissingness)

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
                              stdData  = TRUE
                              #beta     = runif(5, 0.5, 2.0)
                              #tol      = c(0.1, 0.0001)
                              )


library(ggplot2)
library(dplyr)

beta <- rep(1, 5)
type <- "high"
model <- "logistic"

y   <- seq(0, 1000, 0.01)
tmp <- rep(NA, length(y))
i <- 0
for(x in y) {
    i <- i + 1
    tmp[i] <- .fSlopes(x,
                       target = auc,
                       weights = beta,
                       data = X,
                       type = type,
                       model = model)
}


?missing

s2 <- getSigma2(0.25, X, beta)
eta <- as.matrix(X) %*% matrix(beta)

c(0.1, 0.3, 0.5)^2

mean(eta)

y <- eta + rnorm(nrow(X), 0, sqrt(s2))
p <- pnorm(y, mean = quantile(y, 0.25))

m1 <- as.numeric(y < quantile(y, 0.5))
m2 <- as.numeric(y < runif(length(y)))
m3 <- rbinom(length(y), 1, p)

mean(m1)
mean(m2)
mean(m3)

auc(m1, as.numeric(eta))
auc(m2, as.numeric(eta))
auc(m3, as.numeric(eta))

dat <- data.frame(weight = y, objective = tmp)

mean(m)

dat %>%
    filter(weight < 1) %>%
    ggplot(aes(weight, objective)) + geom_line()

tmp <- fSlopes_dump

roc1 <- with(tmp, roc(m, eta))
with(tmp, auc(m, eta))

ls(roc1)

roc1$auc

tmp

dat <- with(tmp, data.frame(m = m, p = p, eta = eta))

ggplot(dat, aes(eta, p, color = m)) + geom_point()

?auc

tmp$p
tmp$m
?roc

## Print the achieved AUC:
out$fit

## Compute proportion missing:
mean(out$m)

## WARNING:
## This approach is highly unstable at small sample sizes (i.e., N << 1000).
## In small samples, the optimization can fail or results (in terms of PM and
## AUC) can be negatively biased.

###--------------------------------------------------------------------------###

## Simulate missingness via probit regression
## - Proportion missing = pm
## - AUC for predicting the missingness = auc
## - Missing values in the positive tail of the linear predictor
## - AUC achieved via optimizing the noise values
out <- simProbitMissingness(pm       = pm,
                            auc      = auc,
                            data     = X,
                            type     = "high",
                            optimize = "slopes",
                            stdData  = TRUE,
                            beta     = runif(5),
                            tol      = c(0.1, 0.0001)
                            )

dim(X)

undebug(simProbitMissingness)

## Print the achieved AUC:
out$fit

## Compute proportion missing:
mean(out$m)

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

###--------------------------------------------------------------------------###

## Simulate missingness via a thresholded latent variable
## - Proportion missing = pm
## - R^2 for the latent variable = r2
## - Missing values in the upper tail of the latent variable

debug(simLvtMissingness)

out <- simLvtMissingness(pm     = 0.3,
                         r2     = 0.5,
                         data   = X,
                         type   = "high",
                         retAuc = TRUE)

## Print the achieved AUC:
out$auc

## Compute proportion missing:
mean(out$m)

ls(out)

with(out, cor(eta, eta + noise))^2

dat <- data.frame(out$m, X)
colnames(dat) <- c("m", paste0("x", 1:ncol(X)))

fit <- glm(m ~ x1 + x2 + x3 + x4 + x5,
                                        #I(x1^2) + I(x2^2) + I(x3^2) + I(x4^2) + I(x5^2),
           data = dat,
           family = "binomial")

summary(fit)

library(DescTools)

PseudoR2(fit, "all")
