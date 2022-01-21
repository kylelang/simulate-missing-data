### Title:    Visualize Simulation Results
### Author:   Kyle M. Lang
### Created:  2021-03-05
### Modified: 2021-03-05

rm(list = ls(all = TRUE))

library(ggplot2)

source("processResults.R")

###--Optimized Slopes--------------------------------------------------------###

head(error)
error[1 : 6] <- lapply(error[1 : 6], as.factor)

dat <- error[error$dist == "gamma", ]

p <- ggplot(data = dat,
            mapping = aes(y = count, x = pm, color = n, shape = p)
            ) +
    geom_jitter()
p + facet_grid(rows = vars(type), cols = vars(auc))

dat <- error[error$n == 500, ]

p <- ggplot(data = dat,
            mapping = aes(y = count, x = pm, color = dist, shape = p)
            ) +
    geom_jitter()
p + facet_grid(rows = vars(type), cols = vars(auc))

head(prb)
prb[1 : 6] <- lapply(prb[1 : 6], as.factor)

dat <- prb[prb$dist == "norm" & prb$n == 500, ]

p <- ggplot(data = dat,
            mapping = aes(y = prbAUC, x = pm, color = n, shape = p)
            ) +
    geom_jitter()
p + facet_grid(rows = vars(type), cols = vars(auc))

dat <- prb[prb$n == 500, ]

p <- ggplot(data = dat,
            mapping = aes(y = prbAUC, x = pm, color = dist, shape = p)
            ) +
    geom_jitter()
p + facet_grid(rows = vars(type), cols = vars(auc))

###--Fixed Slopes------------------------------------------------------------###

head(prb)
prb[1 : 5] <- lapply(prb[1 : 5], as.factor)

p <- ggplot(data = prb,
            mapping = aes(y = prb, x = pm, color = n, shape = p)
            ) +
    geom_jitter()
p + facet_grid(rows = vars(type), cols = vars(dist))

head(auc)
auc[1 : 5] <- lapply(auc[1 : 5], as.factor)

p <- ggplot(data = auc,
            mapping = aes(y = auc, x = pm, color = n, shape = p)
            ) +
    geom_jitter()
p + facet_grid(rows = vars(type), cols = vars(dist))

###--Fixed SNR---------------------------------------------------------------###

sum(error$count)

head(prb)
prb[1 : 6] <- lapply(prb[1 : 6], as.factor)

dat <- prb[prb$dist == "norm", ]

p <- ggplot(data = dat,
            mapping = aes(y = prbSNR, x = pm, color = n, shape = p)
            ) +
    geom_jitter()
p + facet_grid(rows = vars(type), cols = vars(snr))
