### Title:    Data Simulation Routines
### Author:   Kyle M. Lang
### Created:  2021-03-21
### Modified: 2021-03-21

## Simulate data to match the setup described by Hardt et al (2012):
hardtEtAl2012DataSim <- function(n, rZ, pZ)
{
    ## Define the correlation matrix of all the variables. Start by filling the
    ## correlations involving the auxiliaries:
    sigma       <- matrix(rZ, 3 + pZ, 3 + pZ)
    diag(sigma) <- 1.0
    
    ## Add the correlations between X and Y. This value is implied by beta1 =
    ## beta2 = 1 and R^2 = 0.4:
    sigma[1, 2 : 3] <- sigma[2 : 3, 1] <- 0.2 * sqrt(7)
    
    ## Add the correlations between x1 and x2:
    sigma[2, 3] <- sigma[3, 2] <- 0.4
    
    ## Simulate MVN data assuming zero means (not specified in the paper):
    dat1 <- as.data.frame(rmvnorm(n, sigma = sigma))

    ## Add appropriate column names:
    cn <- c("y", "x1", "x2")
    if(pZ > 0) cn <- c(cn, paste0("z", 1 : pZ))
    
    colnames(dat1) <- cn 
    
    dat1
}

