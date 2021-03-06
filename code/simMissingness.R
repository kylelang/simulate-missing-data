### Title:    Simulate MAR Missingness via Logistic Regression
### Author:   Kyle M. Lang
### Created:  2019-11-06
### Modified: 2021-03-04

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
    eta <- as.numeric(data %*% matrix(slope * weights))

    if(type %in% c("center", "tails")) eta <- abs(eta)
    
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

## Simulate a nonresponse vector:
simMissingness <- function(pm,
                           data,
                           auc       = NULL,
                           preds     = colnames(data),
                           type      = "high",
                           beta      = rep(1.0, length(preds)),
                           stdEta    = FALSE,
                           stdData   = !stdEta,
                           ...)
{
    ## Standardize the missing data predictors:
    if(stdData) data <- scale(data)

    ## Find optimal slope values:
    if(!is.null(auc)) {
        optOut1 <-
            .optSlopes(auc = auc, weights = beta, data = data, type = type, ...)
        beta    <- optOut1$minimum * beta
    }

    ## Define the (centered) linear predictor:
    eta <- as.matrix(data) %*% matrix(beta)
    if(stdEta) eta <- scale(eta)
        
    ## Find an optimal intercept value:
    optOut2 <- .optIntercept(pm = pm, eta = eta, type = type, ...)

    ## Define the optimized intercept:
    intercept <- optOut2$minimum
                                           
    ## Compute the probabilities of nonresponse:
    probs <- plogis(
        switch(type,
               high   = intercept + eta,
               low    = intercept - eta,
               center = intercept - abs(eta),
               tails  = abs(eta) - intercept
               )
    )
    
    list(r   = as.logical(rbinom(n = length(eta), size = 1, prob = probs)),
         p   = probs,
         eta = eta,
         b0  = intercept,
         b1  = beta,
         auc = ifelse(is.null(auc), NA, sqrt(optOut1$objective) + auc)
         )
}
