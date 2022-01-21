### Title:    Simulate MAR Missingness via Logistic Regression
### Author:   Kyle M. Lang
### Created:  2019-11-06
### Modified: 2022-01-21

###--------------------------------------------------------------------------###

## Define the missingess vector for the linear probability model:
.linProbMissingness <- function(eta, pm, type) 
    switch(type,
           high   = eta > quantile(eta, 1 - pm),
           low    = eta < quantile(eta, pm),
           center = eta > quantile(eta, (1 - pm) / 2) &
               eta < quantile(eta, pm + (1 - pm) / 2),
           tails  = eta < quantile(eta, pm / 2) |
               eta > quantile(eta, 1 - (pm / 2))
           )

###--------------------------------------------------------------------------###

## Objective function for the intercept:
.fIntercept <- function(intercept, eta, pm, type)
{
    p <- plogis(
        switch(type,
               center = intercept - abs(eta),
               tails  = abs(eta) - intercept,
               intercept + eta
               )
    )
    (mean(p) - pm)^2
}

###--------------------------------------------------------------------------###

## Objective function for the slopes:
.fSlopes <- function(slope, target, weights, data, type, ...)
{
    eta <- as.numeric(as.matrix(data) %*% matrix(slope * weights))

    if(type %in% c("center", "tails")) eta <- abs(eta)
    
    p <- plogis(eta)
    r <- rbinom(length(p), 1, p)
    
    (auc(roc(r, p, quiet = TRUE, ...)) - target)^2
}

###--------------------------------------------------------------------------###

## Objective function for the proportion of noise:
.fNoise <- function(weight, target, eta, noise, pm, type, logistic, ...)
{
    if(logistic & type %in% c("center", "tails")) eta <- abs(eta)

    ## Add noise to the linear predictor:
    eta2 <- eta + weight * noise
    
    if(logistic) {
        p   <- plogis(eta2)
        r   <- rbinom(length(p), 1, p)
        roc <- roc(r, p, quiet = TRUE, ...)
    }
    else {
        r   <- .linProbMissingness(eta = eta2, pm = pm, type = type)
        roc <- roc(r, eta, quiet = TRUE, ...)
    }
    
    (auc(roc) - target)^2
}

###--------------------------------------------------------------------------###

## Optimize the logistic regression intercept for a given value of the linear
## predictor (eta) to get a desired percent missing (pm):
.optIntercept <- function(pm,
                          eta,
                          type,
                          tol     = c(0.1, 0.001),
                          maxIter = 10,
                          ...)
{
    for(k in 1 : maxIter) {
        ## Define the search range:
        int <- k * range(eta)
        
        ## Optimize the objective over 'int':
        out <- optimize(f        = .fIntercept,
                        interval = int,
                        eta      = eta,
                        pm       = pm,
                        type     = type)
        
        ## Are we far enough from the boundary?
        dist   <- out$minimum - int
        check1 <- all(abs(dist) > tol[1] * diff(int))
        
        ## Are we within tolerance?
        check2 <- out$objective < tol[2]
        
        if(check1 & check2) break
    }
    ## Did we fail?
    if(!check1 | ! check2) stop("I could not optimize the intercept.")
    
    out
}

###--------------------------------------------------------------------------###

## Optimize the logistic regression slopes for a given set of predictors to get
## a desired area under the ROC (auc):
.optSlopes <- function(auc,
                       weights,
                       data,
                       type,
                       tol     = c(0.1, 0.001),
                       maxIter = 10,
                       ...)
{
    for(k in 1 : maxIter) {
        ## Define the search range:
        int <- k * range(data)
        
        ## Optimize the objective over 'int':
        out <- optimize(f        = .fSlopes,
                        interval = int,
                        target   = auc,
                        weights  = weights,
                        data     = data,
                        type     = type,
                        ...)
        
        ## Are we far enough from the boundary?
        dist   <- out$minimum - int
        check1 <- all(abs(dist) > tol[1] * diff(int))
        
        ## Are we within tolerance?
        check2 <- out$objective < tol[2]
        
        if(check1 & check2) break
    }
    ## Did we fail?
    if(!check1 | ! check2) stop("I could not optimize the slopes.")
    
    out
}

###--------------------------------------------------------------------------###

## Optimize the proportion of noise to get a desired area under the ROC (auc):
.optNoise <- function(auc,
                      eta,
                      noise,
                      pm,
                      type,
                      logistic,
                      tol     = c(0.1, 0.001),
                      maxIter = 10,
                      ...)
{
    for(k in 1 : maxIter) {
        ## Define the search range:
        int <- c(0, k * 2)
        
        ## Optimize the objective over 'int':
        out <- optimize(f        = .fNoise,
                        interval = int,
                        target   = auc,
                        eta      = eta,
                        noise    = noise,
                        pm       = pm,
                        type     = type,
                        logistic = logistic,
                        ...)
        
        ## Are we far enough from the boundary?
        dist   <- out$minimum - int
        check1 <- all(abs(dist) > tol[1] * diff(int))
        
        ## Are we within tolerance?
        check2 <- out$objective < tol[2]
        
        if(check1 & check2) break
    }
    ## Did we fail?
    if(!check1 | ! check2) stop("I could not optimize the noise.")
    
    out
}

###--------------------------------------------------------------------------###

## Simulate a nonresponse vector via logistic regression:
simLogisticMissingness <- function(pm,
                                   data,
                                   auc       = NULL,
                                   snr       = NULL,
                                   optimize  = "intercept",
                                   preds     = colnames(data),
                                   type      = "high",
                                   beta      = rep(1.0, length(preds)),
                                   stdEta    = FALSE,
                                   stdData   = !stdEta,
                                   ...)
{
                                        #if(is.null(snr) & is.null(auc))
                                        #    stop("You must define a value of either 'snr' or 'auc'")

    if(!is.data.frame(data))
        stop("'data' must be a data.frame")
    
    ## Extract the MAR predictors:
    data <- data[preds]
    
    ## Standardize the missing data predictors:
    if(stdData & optimize != "noise") data <- scale(data)
    
    ## Find optimal slope values:
    if(is.null(snr) & optimize == "slopes") {
        b1Fit <-
            .optSlopes(auc = auc, weights = beta, data = data, type = type, ...)
        beta    <- b1Fit$minimum * beta
    }
    
    ## Define the (centered) linear predictor:
    eta <- as.numeric(as.matrix(data) %*% matrix(beta))
    if(stdEta | optimize == "noise")
        eta <- as.numeric(scale(eta))

    ## Find the optimal proportion of noise:
    if(is.null(snr) & optimize == "noise") {
        noise    <- rnorm(length(eta), 0, sd(eta))
        noiseFit <- .optNoise(auc      = auc,
                              eta      = eta,
                              noise    = noise,
                              type     = type,
                              pm       = pm,
                              logistic = TRUE,
                              ...)

        ## Save the variance of the raw linear predictor:
        v0 <- var(eta)

        ## Add noise to the linear predictor:
        eta <- eta + noiseFit$minimum * noise
    }
    
    ## Define the noisy linear predictor in terms of the specified SNR:
    if(!is.null(snr)) {
        v0  <- var(eta)
        eta <- eta + (1 / snr) * rnorm(length(eta), 0, sd(eta))
    }
    
    ## Find an optimal intercept value:
    b0Fit <- .optIntercept(pm = pm, eta = eta, type = type, ...)
    
    ## Define the optimized intercept:
    intercept <- b0Fit$minimum
                                           
    ## Compute the probabilities of nonresponse:
    probs <- plogis(
        switch(type,
               high   = intercept + eta,
               low    = intercept - eta,
               center = intercept - abs(eta),
               tails  = abs(eta) - intercept
               )
    )

    r <- as.logical(rbinom(n = length(eta), size = 1, prob = probs))    
    
    if(is.null(snr))
        fit <- auc(roc(r, probs, smooth = TRUE))
                                        #    fit <- ifelse(optimize == "noise",
                                        #                  sqrt(noiseFit$objective) + auc,
                                        #                  sqrt(b1Fit$objective) + auc)
    else
        fit <- sqrt(v0) / sqrt(var(eta) - v0)
    
    list(r   = r,
         p   = probs,
         eta = eta,
         b0  = intercept,
         b1  = beta,
         fit = fit 
         )
}

###--------------------------------------------------------------------------###

## Simulate a nonresponse vector via logistic regression without controlling the
## strength of the MAR assocation:
simLogisticMissingness0 <- function(pm,
                                    data,
                                    preds   = colnames(data),
                                    type    = "high",
                                    beta    = rep(1.0, length(preds)),
                                    stdData = FALSE,
                                    ...)
    simLogisticMissingness(pm      = pm,
                           data    = data,
                           preds   = preds,
                           type    = type,
                           beta    = beta,
                           stdData = stdData,
                           stdEta  = FALSE,
                           ...)
                           
###--------------------------------------------------------------------------###

## Simulate a nonresponse vector via a linear probability model:
simLinearMissingness <- function(pm,
                                 data,
                                 auc      = NULL,
                                 snr      = NULL,
                                 optimize = TRUE,
                                 preds    = colnames(data),
                                 type     = "high",
                                 beta     = rep(1.0, length(preds)),
                                 stdData  = TRUE,
                                 ...)
{
    if(is.null(snr) & is.null(auc))
        stop("You must define a value of either 'snr' or 'auc'")

    if(optimize & !type %in% c("high", "low"))
        stop("'type' must be either 'high' or 'low'")
    
    if(!is.data.frame(data))
        stop("'data' must be a data.frame")
    
    ## Extract the MAR predictors:
    data <- data[preds]
    
    ## Standardize the missing data predictors:
    if(stdData) data <- scale(data)
    
    ## Define the (centered) linear predictor:
    eta <- as.numeric(as.matrix(data) %*% matrix(beta))
    
    ## Find the optimal proportion of noise:
    if(optimize) {
        noise    <- rnorm(length(eta), 0, sd(eta))
        noiseFit <- .optNoise(auc      = auc,
                              eta      = eta,
                              noise    = noise,
                              type     = type,
                              pm       = pm,
                              logistic = FALSE,
                              ...)
        
        ## Add noise to the linear predictor:
        eta2 <- eta + (noiseFit$minimum * noise)
    }
    else
        ## Define the noisy linear predictor in terms of the specified SNR:
        eta2 <- eta + (1 / snr) * rnorm(length(eta), 0, sd(eta))
    
    r <- .linProbMissingness(eta = eta2, pm = pm, type = type)

    ## Compute the achieved AUC:
    auc <- ifelse(type %in% c("high", "low"),
                  as.numeric(auc(roc(r, eta, quiet = TRUE, ...))),
                  NA)
    
    list(r   = r,
         eta = eta2,
         auc = auc,
         snr = sd(eta) / sqrt(var(eta2) - var(eta))
         )
}

###--------------------------------------------------------------------------###

## Generate deterministic missingness based on quantiles of the  MAR predictor
simSimpleMar <- function(data,
                         pm,
                         preds,
                         type,
                         beta    = rep(1, length(preds)),
                         stdData = FALSE)
{
    
    ## Extract the MAR predictors:
    data <- data[preds]
    
    ## Standardize the missing data predictors:
    if(stdData) data <- scale(data)
    
    ## Define the (centered) linear predictor:
    eta <- as.numeric(as.matrix(data) %*% matrix(beta))
    
    switch(type,
           high   = eta > quantile(eta, 1 - pm),
           low    = eta < quantile(eta, pm),
           center = eta > quantile(eta, 0.5 - pm / 2) &
               eta < quantile(eta, 0.5 + pm / 2),
           tails  = eta < quantile(eta, pm / 2) &
               eta > quantile(eta, 1 - pm / 2)
           )
}
