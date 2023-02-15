### Title:    Simulate MAR Missingness via Logistic Regression
### Author:   Kyle M. Lang
### Created:  2019-11-06
### Modified: 2023-02-15

###--------------------------------------------------------------------------###

## Define the missingess vector by thresholding a latent variable:
.makeMissingness <- function(y, pm, type)
    switch(
        type,
        high   = y > quantile(y, 1 - pm),
        low    = y < quantile(y, pm),
        center = y > quantile(y, (1 - pm) / 2) & y < quantile(y, pm + (1 - pm) / 2),
        tails  = y < quantile(y, pm / 2) | y > quantile(y, 1 - (pm / 2))
    )

###--------------------------------------------------------------------------###

## Objective function for the intercept:
.fIntercept <- function(intercept, eta, pm, type, model) {
    f <- switch(model,
                logistic = plogis,
                probit   = pnorm,
                stop("I can only optimize the intercept for logistic or probit models. I don't know what to do with a '", model, "' model.")
    )
    p <- f(
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
.fSlopes <- function(slope, target, weights, data, type, model, ...) {
    eta <- as.numeric(as.matrix(data) %*% matrix(slope * weights))

    if(type %in% c("center", "tails")) eta <- abs(eta)

    p <- switch(model,
                logistic = plogis(eta),
                probit   = pnorm(eta),
                stop("I can only optimize the slopes for logistic or probit models. I don't know what to do with a '", model, "' model.")
                )
    m <- rbinom(length(p), 1, p)

    (auc(m, eta, quiet = TRUE, ...) - target)^2
}

###--------------------------------------------------------------------------###

## Objective function for the proportion of noise:
.fNoise <- function(weight, target, eta, noise, pm, type, model, ...) {
    check <- model %in% c("logistic", "probit") & type %in% c("center", "tails")
    if(check) eta <- abs(eta)

    ## Add noise to the linear predictor:
    eta2 <- eta + weight * noise

    if(model %in% c("logistic", "probit")) {
        p   <- switch(model,
                      logistic = plogis(eta2),
                      probit   = pnorm(eta2)
                      )
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
                          model,
                          tol     = c(0.001, 0.0001),
                          maxIter = 10,
                          ...) {
    for(k in 1 : maxIter) {
        ## Define the search range:
        int <- k * range(eta)

        ## Optimize the objective over 'int':
        out <- optimize(f        = .fIntercept,
                        interval = int,
                        eta      = eta,
                        pm       = pm,
                        type     = type,
                        model    = model)

        ## Are we far enough from the boundary?
        dist   <- out$minimum - int
        check1 <- all(abs(dist) > tol[1])

        ## Are we within tolerance?
        check2 <- out$objective < tol[2]

        if(check1 & check2) break
    }
    ## Did we fail?
    if(!check1 | !check2) stop("I could not optimize the intercept.")

    out
}

###--------------------------------------------------------------------------###

## Optimize the logistic regression slopes for a given set of predictors to get
## a desired area under the ROC (auc):
.optSlopes <- function(auc,
                       weights,
                       data,
                       type,
                       model,
                       tol     = c(0.001, 0.0001),
                       maxIter = 100,
                       ...) {
    for(k in 1:maxIter) {
        ## Define the search range:
        int <- c(0, max(data)) * k

        ## Optimize the objective over 'int':
        out <- optimize(f        = .fSlopes,
                        interval = int,
                        target   = auc,
                        weights  = weights,
                        data     = data,
                        type     = type,
                        model    = model,
                        ...)

        ## Are we far enough from the boundary?
        dist   <- out$minimum - int
        check1 <- all(abs(dist) > tol[1])

        ## Are we within tolerance?
        check2 <- out$objective < tol[2]

        print(paste("k:", k))
        print(paste("int:", paste(int, collapse = ":")))
        print(paste("max(data):", max(data)))
        print(paste("Objective:", out$objective))
        print(paste("Slope:", out$minimum))
        print(paste("check1:", check1))
        print(paste("check2:", check2))
        print("")

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
                      model,
                      tol     = c(0.1, 0.001),
                      maxIter = 10,
                      ...) {
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
                        model    = model,
                        ...)

        ## Are we far enough from the boundary?
        dist   <- out$minimum - int
        check1 <- all(abs(dist) > tol[1] * diff(int))

        ## Are we within tolerance?
        check2 <- out$objective < tol[2]

        if(check1 & check2) break
    }
    ## Did we fail?
    if(!check1 | !check2) stop("I could not optimize the noise.")

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
                                   ...) {
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
        b1Fit <- .optSlopes(auc     = auc,
                            weights = beta,
                            data    = data,
                            type    = type,
                            model   = "logistic",
                            ...)
        beta <- b1Fit$minimum * beta
    }

    ## Define the (centered) linear predictor:
    eta <- as.numeric(as.matrix(data) %*% matrix(beta))
    if(stdEta | optimize == "noise")
        eta <- as.numeric(scale(eta))

    ## Find the optimal proportion of noise:
    if(is.null(snr) & optimize == "noise") {
        noise    <- rnorm(length(eta), 0, sd(eta))
        noiseFit <- .optNoise(auc   = auc,
                              eta   = eta,
                              noise = noise,
                              type  = type,
                              pm    = pm,
                              model = "logistic",
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
    b0Fit <- .optIntercept(pm    = pm,
                           eta   = eta,
                           type  = type,
                           model = "logistic",
                           ...)

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

    m <- as.logical(rbinom(n = length(eta), size = 1, prob = probs))

    if(is.null(snr))
        fit <- auc(roc(m, probs, smooth = TRUE))
                                        #    fit <- ifelse(optimize == "noise",
                                        #                  sqrt(noiseFit$objective) + auc,
                                        #                  sqrt(b1Fit$objective) + auc)
    else
        fit <- sqrt(v0) / sqrt(var(eta) - v0)

    list(m   = m,
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

## Simulate a nonresponse vector via probit regression:
simProbitMissingness <- function(pm,
                                 data,
                                 auc       = NULL,
                                 snr       = NULL,
                                 optimize  = "intercept",
                                 preds     = colnames(data),
                                 type      = "high",
                                 beta      = rep(1.0, length(preds)),
                                 stdEta    = FALSE,
                                 stdData   = !stdEta,
                                 ...) {
                                        #if(is.null(snr) & is.null(auc))
                                        #    stop("You must define a value of either 'snr' or 'auc'")

    if(!is.data.frame(data))
        stop("'data' must be a data.frame")

    ## Extract the MAR predictors:
    data <- data[preds]

    ## Standardize the missing data predictors:
    if(stdData) data <- scale(data)

    ## Find optimal slope values:
    if(is.null(snr) & optimize == "slopes") {
        b1Fit <- .optSlopes(auc     = auc,
                            weights = beta,
                            data    = data,
                            type    = type,
                            model   = "probit",
                            ...)
        beta <- b1Fit$minimum * beta
    }

    ## Define the (centered) linear predictor:
    eta <- as.numeric(as.matrix(data) %*% matrix(beta))
                                        #if(stdEta | optimize == "noise")
                                        #    eta <- as.numeric(scale(eta))

    ## Find the optimal proportion of noise:
    if(is.null(snr) & optimize == "noise") {
        noise    <- rnorm(length(eta), 0, 1) #sd(eta))
        noiseFit <- .optNoise(auc   = auc,
                              eta   = eta,
                              noise = noise,
                              type  = type,
                              pm    = pm,
                              model = "probit",
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
    b0Fit <- .optIntercept(pm    = pm,
                           eta   = eta,
                           type  = type,
                           model = "probit",
                           ...)

    ## Define the optimized intercept:
    intercept <- b0Fit$minimum

    ## Compute the probabilities of nonresponse:
    probs <- pnorm(
        switch(type,
               high   = intercept + eta,
               low    = intercept - eta,
               center = intercept - abs(eta),
               tails  = abs(eta) - intercept
               )
    )

    m <- as.logical(rbinom(n = length(eta), size = 1, prob = probs))

    if(is.null(snr))
        fit <- auc(roc(m, probs, smooth = TRUE))
                                        #    fit <- ifelse(optimize == "noise",
                                        #                  sqrt(noiseFit$objective) + auc,
                                        #                  sqrt(b1Fit$objective) + auc)
    else
        fit <- sqrt(v0) / sqrt(var(eta) - v0)

    list(m   = m,
         p   = probs,
         eta = eta,
         b0  = intercept,
         b1  = beta,
         fit = fit
         )
}

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
                                 ...) {
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

    m <- .makeMissingness(y = eta2, pm = pm, type = type)

    ## Compute the achieved AUC:
    auc <- ifelse(type %in% c("high", "low"),
                  as.numeric(auc(roc(m, eta, quiet = TRUE, ...))),
                  NA)

    list(m   = m,
         eta = eta2,
         auc = auc,
         snr = sd(eta) / sqrt(var(eta2) - var(eta))
         )
}

###--------------------------------------------------------------------------###

## Compute the residual variance necessary to produce a given R^2:
.getSigma2 <- function(r2, eta, x, beta) {
    if(missing(eta)) {
        if(length(beta) > 1) {
            beta <- matrix(beta)
            vX <- t(beta) %*% cov(x) %*% beta
        }
        else
            vX <- beta^2 * var(x)
    }
    else
        vX <- var(eta)

    (vX / r2) - vX
}

###--------------------------------------------------------------------------###

## Simulate missingness based on a thresholded latent response variable:
simLvtMissingness <- function(pm,
                              data,
                              r2      = 1.0,
                              preds   = colnames(data),
                              type    = "high",
                              beta    = rep(1.0, length(preds)),
                              stdData = FALSE,
                              minimal = FALSE,
                              retAuc  = !minimal,
                              ...) {

    if(r2 <= 0 | r2 > 1)
        stop("'r2' must be in the interval (0, 1]")

    if(r2 == 1.0)
        message("You have set 'r2 = 1', so I will simulate a deterministic MAR relation.")

    if(!is.data.frame(data))
        stop("'data' must be a data.frame.")

    ## Extract the MAR predictors:
    data <- data[preds]

    ## Standardize the missing data predictors:
    if(stdData) data <- scale(data)

    ## Define the (deterministic) linear predictor:
    eta <- as.numeric(as.matrix(data) %*% matrix(beta))

    ## Add noise if r2 < 1.0:
    if(r2 < 1.0) {
        ## Define the latent residual variance:
        s2 <- .getSigma2(r2 = r2, eta = eta)

        ## Define the (stochastic) latent response vector:
        noise <- rnorm(length(eta), 0, sqrt(s2))
    }

    ## Threshold y to produce the missingness vector:
    m <- .makeMissingness(y = eta + noise, pm = pm, type = type)

    if(minimal) return(m)

    out <- list(m = m, eta = eta, s2 = s2, noise = noise)

    if(retAuc)
        out$auc <- pROC::auc(m, eta, quiet = TRUE, ...)

    out
}

###--------------------------------------------------------------------------###

## Generate deterministic missingness based on quantiles of the MAR predictor
simSimpleMar <- function(data,
                         pm,
                         preds,
                         type,
                         beta    = rep(1, length(preds)),
                         stdData = FALSE,
                         ...)
    simLvtMissingness(pm      = pm,
                      data    = data,
                      r2      = 1.0,
                      preds   = preds,
                      type    = type,
                      beta    = beta,
                      stdData = stdData,
                              ...)
