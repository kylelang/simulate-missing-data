### Title:    Subroutines for Simulation to Evaluate Missingness Generation
### Author:   Kyle M. Lang
### Created:  2021-03-04
### Modified: 2021-03-10

simData <- function(n, p, dist)
{
    if(dist == "unif")
        as.data.frame(
            matrix(runif(n * p), ncol = p)
        )
    else if(dist == "norm")
        as.data.frame(rmvnorm(n, rep(0, p), diag(p)))
    else if(dist == "gamma")
        as.data.frame(
            matrix(rgamma(n * p, 2.5, 0.5), ncol = p)
        )
    else
        stop("'dist' must be 'unif', 'norm', or 'gamma'")
}

###--------------------------------------------------------------------------###

## Run experiment within one design cell with optimized slopes:
runCell1 <- function(parms, data, ...)
{
    ## Generate some regression weights:
    beta <- runif(parms$p, -1, 1)

    ## Simulate the missingness:
    tmp <- try(
        with(parms,
             simMissingness(pm   = pm,
                            auc  = auc,
                            data = data,
                            type = type,
                            beta = beta,
                            ...)
             ),
        silent = TRUE
    )
    
    ## Return key summaries:
    if(class(tmp) != "try-error")
        list(pm = mean(tmp$r), auc = tmp$auc, e = NA)
    else
        list(pm = NA, auc = NA, e = tmp)
}

###--------------------------------------------------------------------------###

## Run experiment within one design cell with fixed slopes:
runCell2 <- function(parms, data, ...)
{
    ## Generate some regression weights:
    beta <- runif(parms$p, -1, 1)
    
    ## Simulate the missingness:
    tmp <- try(
        with(parms,
             simMissingness(pm     = pm,
                            auc    = NULL,
                            data   = data,
                            type   = type,
                            beta   = beta,
                            stdEta = TRUE,
                            ...)
             ),
        silent = TRUE
    )

    if(class(tmp) != "try-error") {
        if(parms$type %in% c("center", "tails"))
            eta <- abs(tmp$eta)
        else
            eta <- tmp$eta

        ## Compute the AUC:
        fit <- glm(tmp$r ~ eta, family = binomial)
        roc <- try(
            roc(tmp$r, predict(fit), smooth = TRUE, quiet = TRUE),
            silent = TRUE
        )
        
        if(class(roc) != "try-error")
            auc <- auc(roc)
        else
            auc <- NA
        
        ## Return key summaries:
        list(pm = mean(tmp$r), auc = auc, e = NA)
    }
    else
        list(pm = NA, auc = NA, e = tmp)
}

###--------------------------------------------------------------------------###

## Run experiment within one design cell with optimized noise:
runCell3 <- function(parms, data, ...)
{
    ## Generate some regression weights:
    beta <- runif(parms$p, -1, 1)

    ## Simulate the missingness:
    tmp <- try(
        with(parms,
             simMissingness(pm       = pm,
                            auc      = auc,
                            data     = data,
                            type     = type,
                            beta     = beta,
                            optimize = "noise",
                            ...)
             ),
        silent = TRUE
    )
    
    ## Return key summaries:
    if(class(tmp) != "try-error")
        list(pm = mean(tmp$r), auc = tmp$auc, e = NA)
    else
        list(pm = NA, auc = NA, e = tmp)
}

###--------------------------------------------------------------------------###

## Run experiment within one design cell with linear probability models and
## optimized noise:
runCell4 <- function(parms, data, ...)
{
    ## Generate some regression weights:
    beta <- runif(parms$p, -1, 1)

    ## Simulate the missingness:
    tmp <- try(
        with(parms,
             simLinearMissingness(pm       = pm,
                                  auc      = auc,
                                  data     = data,
                                  type     = type,
                                  beta     = beta,
                                  optimize = TRUE,
                                  stdData  = FALSE,
                                  ...)
             ),
        silent = TRUE
    )
    
    ## Return key summaries:
    if(class(tmp) != "try-error")
                                        #list(pm = mean(tmp$r), auc = tmp$auc, snr = tmp$snr, v1 = tmp$v1, v2 = tmp$v2, v3 = tmp$v3, w = tmp$w, e = NA)
        list(pm = mean(tmp$r), auc = tmp$auc, snr = tmp$snr, e = NA)
    else
                                        #list(pm = NA, auc = NA, snr = NA, v1 = NA, v2 = NA, v3 = NA, w = NA, e = tmp)
        list(pm = NA, auc = NA, snr = NA, e = tmp)
}

###--------------------------------------------------------------------------###

## Run experiment within one design cell with linear probability models and
## optimized noise:
runCell5 <- function(parms, data, ...)
{
    ## Generate some regression weights:
    beta <- runif(parms$p, -1, 1)

    ## Simulate the missingness:
    tmp <- try(
        with(parms,
             simLinearMissingness(pm       = pm,
                                  snr      = snr,
                                  data     = data,
                                  type     = type,
                                  beta     = beta,
                                  optimize = FALSE,
                                  stdData  = FALSE,
                                  ...)
             ),
        silent = TRUE
    )
    
    ## Return key summaries:
    if(class(tmp) != "try-error")
        list(pm = mean(tmp$r), auc = tmp$auc, snr = tmp$snr, e = NA)
    else
        list(pm = NA, auc = NA, snr = NA, e = tmp)
}

###--------------------------------------------------------------------------###

doRep <- function(rp, conds, expNum, ...)
{
    newData <- TRUE
    out <- err <- list()
    for(i in 1 : nrow(conds)) {
        n <- conds[i, "n"]
        p <- conds[i, "p"]
        
        ## Do we need to generate a new dataset?
        newData <- n == max(conds$n) & p == max(conds$p)
        
        ## Generate/process the predictor data:
        if(newData)
            dat0 <- simData(n = n, p = p, dist = conds[i, "dist"])
        
        dat1 <- dat0[1 : n, 1 : p]

        ## Choose the appropriate computational kernel:
        runCell <-
            switch(expNum, runCell2, runCell1, runCell3, runCell4, runCell5)
        
        ## Generate the appropriate missingness: 
        tmp      <- runCell(as.list(conds[i, ]), data = dat1, ...)
        out[[i]] <- with(tmp, c(pm, auc, snr))
        err[[i]] <- tmp$e
    }
    
    out           <- do.call(rbind, out)
    colnames(out) <- c("pm", "auc", "snr")
    
    err <- do.call(c, err)
    
    list(out = out, err = err)
}

