source('functionsBonf.R')
library(pbapply)

## Simulation study for paper
## simulate signal for variance known case.
simulSignal <- function(mu, sd = 1, interval_res = NULL, even = FALSE) {
  signalLength <- length(mu)
  x <- rnorm(signalLength, mu, sd = sd)
  z <- c(0,cumsum(x))    
  z2 <- c(0, cumsum(x^2))
  
  if(is.null(interval_res)) {
    interval_res <- BonferroniIntervals(signalLength, var.known = TRUE, even = even) 
  }
  
  rejected_intervals <- getRejectedIntervals(z = z/sd, 
                                             bonferroniIntervals = interval_res$intervals,
                                             blocksize = interval_res$blockSize,
                                             alpha = alpha, 
                                             var.known = TRUE)
  conf_intervals <- rejectedToConf(rejected_intervals)
  ci_results <- getCIs(conf_intervals)
  return(c(lb = nrow(ci_results$disjoint), 
           cover = checkCI(ci_results, getCptIdx(mu))))
}



res_list = list()
set.seed(0919)

nsim <- 10000
alpha <- 0.1

## Null -----
n <- 1000
interval_res <- BonferroniIntervals(n, var.known = TRUE)
mu <- 0
res_list[[paste0('null', n)]] <- list(mu = rep(mu, n),
                                      result = pbreplicate(n = nsim, simulSignal(mu = rep(mu,n),
                                                                                 interval_res = interval_res)))

n <- 2000
interval_res <- BonferroniIntervals(n, var.known = TRUE)
mu <- 0
res_list[[paste0('null', n)]] <- list(mu = rep(mu, n),
                                      result = pbreplicate(n = nsim, simulSignal(mu = rep(mu,n),
                                                                                 interval_res = interval_res)))

n <- 3000
interval_res <- BonferroniIntervals(n, var.known = TRUE)
mu <- 0
res_list[[paste0('null', n, sep = "")]] <- list(mu = rep(mu, n),
                                                result = pbreplicate(n = nsim, simulSignal(mu = rep(mu,n),
                                                                                           interval_res = interval_res)))


## Blocks ----
# true num = 11
n <- 2048
mu=c(rep(0,204),rep(14.64,62),rep(-3.66,308-267),rep(7.32,472-308),rep(-7.32,512-472),rep(10.98,820-512))
mu=c(mu,rep(-4.39,902-820),rep(3.29,1332-902),rep(19.03,1557-1332),rep(7.68,1598-1557),rep(15.37,1659-1598))
mu=c(mu,rep(0,n+1-1659))
mu = mu/10
interval_res <- BonferroniIntervals(n, var.known = TRUE)
res_list[['blocks']] <- list(mu = mu,
                             result = pbreplicate(n = nsim, simulSignal(mu = mu,
                                                                        interval_res = interval_res)))

# fms -----
# true num = 6
n <- 497 
mu=c(rep(-0.18,138),rep(0.08,226-139),rep(1.07,243-226),rep(-0.53,300-243),rep(0.16,309-300),rep(-0.69,333-309))
mu=c(mu,rep(-0.16,n+1-333))
mu = mu/0.3
interval_res <- BonferroniIntervals(n, var.known = TRUE)
res_list[['fms']] <- list(mu = mu,
                          result = pbreplicate(n = nsim, simulSignal(mu = mu,
                                                                     interval_res = interval_res)))

# mix -------
# true num = 13
n <- 560
mu=c(rep(7,10),rep(-7,21-11),rep(6,41-21),rep(-6,61-41),rep(5,91-61),rep(-5,121-91))
mu=c(mu,c(rep(4,161-121),rep(-4,201-161),rep(3,251-201),rep(-3,301-251),rep(2,361-301),rep(-2,421-361)))
mu=c(mu,c(rep(1,491-421),rep(-1,n+1-491)))
mu=mu/4
interval_res <- BonferroniIntervals(n, var.known = TRUE)
res_list[['mix']] <- list(mu = mu,
                          result = pbreplicate(n = nsim, simulSignal(mu = mu,
                                                                     interval_res = interval_res)))

# teeth10 -----
# true num = 13
n <- 140
mu=(rep(c(rep(0,10),rep(1,10)),7))/0.4
interval_res <- BonferroniIntervals(n, var.known = TRUE)
res_list[['teeth10']] <- list(mu = mu,
                              result = pbreplicate(n = nsim, simulSignal(mu = mu,
                                                                         interval_res = interval_res)))

# stairs10 ------
# true num = 14
n <- 150
mu=rep(1:15, each = 10)/0.3
interval_res <- BonferroniIntervals(n, var.known = TRUE)
res_list[['stairs10']] <- list(mu = mu,
                               result = pbreplicate(n = nsim, simulSignal(mu = mu,
                                                                          interval_res = interval_res)))



lapply(res_list, function(x) length(x$mu))
lapply(res_list, function(x) getNumCpt(x$mu))
# N(\alpha) and \hat{p}_1
lapply(res_list, function(x) rowMeans(x$result))
lapply(res_list, function(x) apply(x$result, 1, sd))
# \hat{p}_2
lapply(res_list, function(x) mean(x$result['lb',] <= getNumCpt(x$mu)))
# lb - K
lapply(res_list, function(x) table(x$result['lb',] - getNumCpt(x$mu))/ncol(x$result))


