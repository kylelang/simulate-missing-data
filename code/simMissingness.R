### Title:    Simulate MAR Missingness via Logistic Regression
### Author:   Kyle M. Lang
### Created:  2019-11-06
### Modified: 2021-03-08

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
.fNoise <- function(weight, target, eta, noise, type, ...)
{
    if(type %in% c("center", "tails")) eta <- abs(eta)

    ## Add noise to the linear predictor:
    eta <- eta + weight * noise
    
    p <- plogis(eta)
    r <- rbinom(length(p), 1, p)

    (auc(roc(r, p, quiet = TRUE, ...)) - target)^2
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
                      type,
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
    if(!check1 | ! check2) stop("I could not optimize the noise.")
    
    out
}

###--------------------------------------------------------------------------###

#n    <- 1000
#p    <- 5
#pm   <- 0.25
#auc  <- 0.75
#beta <- rep(10, p)
#type <- "high"
#data <- as.data.frame(rmvnorm(n, rep(0, p), diag(p)))

#stdData <- FALSE
#stdEta  <- TRUE

## Simulate a nonresponse vector:
simMissingness <- function(pm,
                           data,
                           auc       = NULL,
                           snr       = NULL,
                           optimize  = "noise",
                           preds     = colnames(data),
                           type      = "high",
                           beta      = rep(1.0, length(preds)),
                           stdEta    = FALSE,
                           stdData   = !stdEta,
                           ...)
{
    if(is.null(snr) & is.null(auc))
        stop("You must define a value of either 'snr' or 'auc'")
    
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
        noiseFit <-
            .optNoise(auc = auc, eta = eta, noise = noise, type = type) #, ...)
        eta      <- eta + noiseFit$minimum * noise
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

    if(is.null(snr))
        fit <- ifelse(optimize == "noise",
                      sqrt(noiseFit$objective) + auc,
                      sqrt(b1Fit$objective) + auc)
    else
        fit <- v0 / (var(eta) - v0)
    
    list(r   = as.logical(rbinom(n = length(eta), size = 1, prob = probs)),
         p   = probs,
         eta = eta,
         b0  = intercept,
         b1  = beta,
         fit = fit 
         )
}
