### Title:    Process Simulation Results
### Author:   Kyle M. Lang
### Created:  2021-03-05
### Modified: 2021-03-10

rm(list = ls(all = TRUE))

outDir <- "../output/optSlopes/"

conds <- readRDS(paste0(outDir, "conds.rds"))
out   <- readRDS(paste0(outDir, "evaluation_output.rds"))

pool  <- data.frame(pm  = rep(0, nrow(conds)),
                    auc = rep(0, nrow(conds)),
                    snr = rep(0, nrow(conds))
                    )
error <- rep(0, nrow(conds))

## Aggregate the results:
for(x in out) {
    tmp             <- x$out
    tmp[is.na(tmp)] <- 0

    pool  <- pool + tmp
    error <- error + !is.na(x$e)
}

pool <- pool / length(out)

head(pool)

## Fixed Slopes:
                                        #prb <- 100 * (pool$pm - conds$pm) / conds$pm

## Optimized Slopes:
true <- conds[c("pm", "auc")]

prb           <- 100 * (pool[c("pm", "auc")] - true) / true
colnames(prb) <- paste0("prb", c("PM", "AUC"))

auc   <- data.frame(conds, auc = pool$auc)
prb   <- data.frame(conds, prb)
snr   <- data.frame(conds, snr = pool$snr)
error <- data.frame(conds, count = error)

## Fixed SNR:
                                        #true <- conds[c("pm", "snr")]

                                        #prb           <- 100 * (pool[c("pm", "snr")] - true) / true
                                        #colnames(prb) <- paste0("prb", c("PM", "SNR"))

                                        #auc   <- data.frame(conds, auc = pool$auc)
                                        #prb   <- data.frame(conds, prb)
                                        #error <- data.frame(conds, count = error)

