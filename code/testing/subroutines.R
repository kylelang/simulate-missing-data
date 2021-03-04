### Title:    Subroutines for Simulation to Evaluate Missingness Generation
### Author:   Kyle M. Lang
### Created:  2021-03-04
### Modified: 2021-03-04

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

runCell <- function(parms, data)
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
                            beta = beta)
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

doRep <- function(rp, conds)
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
        
        ## Generate the appropriate missingness: 
        tmp      <- runCell(as.list(conds[i, ]), data = dat1)
        out[[i]] <- with(tmp, c(pm, auc))
        err[[i]] <- tmp$e
    }
    
    out           <- do.call(rbind, out)
    colnames(out) <- c("pm", "auc")
    
    err <- do.call(c, err)
    
    list(conds = conds, out = out, err = err)
}

