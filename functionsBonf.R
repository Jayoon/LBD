library(ggplot2)
### functions for Bonferroni Intervals

### Get grid value of J_l
dl <- function(n, l) {
  ml <- 2^l
  ceiling(ml/sqrt(2*log(exp(1)*n/ml)))
}

Jl <- function(n, l, even = FALSE) {
  first <- TRUE
  grid_val <- dl(n, l)
  lower <- ceiling(2^l/grid_val)
  upper <- ceiling(2^(l+1)/grid_val)-1
  for(size in lower:upper) {
    if(!even | (size * grid_val) %% 2 == 0) {
      start_pt <- 0:(floor(n/grid_val)-size)
      if(first) {
        result <- cbind(start_pt, start_pt + size)
        first <- FALSE
      } else {
        result <- rbind(result, cbind(start_pt, start_pt + size))
      }
    }
  }
  return(result*grid_val)
}

BonferroniIntervals <- function(n, 
                                var.known = FALSE,
                                three.split = FALSE,
                                even = FALSE) {
  largestwindow <- n/4
  sn <- ceiling(log2(log(n)))
  Bmax <- floor(log2(largestwindow)) - sn +1
  if(!var.known & !three.split & n < 55) {
    stop("Sample size should be greater than 55 when variance is unknown")
  } 
  if(!var.known & three.split & n < 2980) {
    stop("Sample size should be greater than 2980 when variance is unknown and using three split statistics")
  } 
  bonferroniIntervals <- NA
  blocksize <- numeric(length = Bmax)
  intervalCnt <- 0
  # If variance is known, start with ell = 1, 
  # if not start with ell = 2 when using one split, ell = 3 when using three splits
  block1start <- ifelse(var.known, 1, ifelse(three.split, 3, 2)) 
  for(ell in block1start:(sn-1)) { ### Block 1
    intervals <- cbind(Jl(n, ell, even = even), ell)
    intervalCnt = intervalCnt + nrow(intervals)
    if(ell == block1start) {
      bonferroniIntervals <- intervals
    } else {
      bonferroniIntervals <- rbind(bonferroniIntervals, intervals)
    }
  }
  bonferroniIntervals <- cbind(bonferroniIntervals, 1)
  blocksize[1] <- intervalCnt
  for(B in 2:Bmax) {
    intervals <- cbind(Jl(n, B-2+sn, even = even), B-2+sn)
    blocksize[B] <- nrow(intervals)
    bonferroniIntervals <- rbind(bonferroniIntervals, cbind(intervals, B))
  }
  colnames(bonferroniIntervals) <- c("lower", "upper", "level", "block")
  bonferroniIntervals <- data.frame(bonferroniIntervals)
  res <- list(intervals = bonferroniIntervals,
              blockSize = blocksize)
  res
}


### Compute the standardized mean difference in interval (j,k]
### w.r.t point m. 
### z = c(0, cusum(x)), z2 = c(0, cumsum(x^2))
### returns NaN if j+1 = k for var.known case (i.e. 1 obsv)
### returns error if m-j or k-m <= 1 for var.unknown case
testStat <- function(j, k, m, z, z2 = NULL, var.known = FALSE) {
  # m = ceiling((j+k)/2)
  na = m-j
  nb = k-m
  if(var.known) {
    sqrt(na*nb/(k-j))*abs((z[m+1]-z[j+1])/na - (z[k+1]-z[m+1])/nb)
  } else {
    # xx = diff(z[(j+1):(k+1)])    # gives x[(j+1):k]
    # sqrt(na*nb/(k-j))*abs((z[m+1]-z[j+1])/na - (z[k+1]-z[m+1])/nb)/sd(xx)
    if(na <= 1 | nb <= 1) {
      stop("Each interval should contain at least 2 observations to estimate the variance")
    }
    muA = (z[m+1]-z[j+1])/na
    muB = (z[k+1]-z[m+1])/nb
    pooled.var = (z2[m+1]-z2[j+1]-na*muA^2+z2[k+1]-z2[m+1]-nb*muB^2)/(na+nb-2)
    if(pooled.var < 0) {
      if(all.equal(pooled.var, 0)) {
        pooled.var = 0
      } else {
        stop("negative variance estimate")
      }
    }
    pooled.sd = sqrt(pooled.var)
    sqrt(na*nb/(k-j))*abs(muA-muB)/pooled.sd
  }
}


getCriticalVal <- function(bonferroniIntervals,
                           blocksize, alpha = 0.1, var.known = FALSE,
                           three.split = FALSE) {
  Bmax = length(blocksize)
  correction_factor <- sum(1/(1:Bmax))
  if(var.known) {
    qnorm(alpha/(2*correction_factor*(blocksize*(1:Bmax))),
          lower.tail = FALSE)
  } else if(!three.split) {
    qt(alpha/(2*correction_factor*(rep(blocksize,blocksize)*bonferroniIntervals$block)),
       lower.tail = FALSE, df = bonferroniIntervals$upper-bonferroniIntervals$lower-2)
  } else { ## Using maximum of the three splits
    qt(alpha/(2*3*correction_factor*(rep(blocksize,blocksize)*bonferroniIntervals$block)),
       lower.tail = FALSE, df = bonferroniIntervals$upper-bonferroniIntervals$lower-2)
  }
}


# z: c(0, cusum(x))
# bonferroniIntervals: data frame output of BonferroniIntervals function
# blocksize: vector output of BonferroniIntervals function
# alpha: type 1 error
getRejectedIntervals <- function(z, z2 = NULL, bonferroniIntervals,
                                 blocksize, alpha = 0.1, var.known = FALSE,
                                 three.split = FALSE) {
  if(var.known) {
    testvec <- mapply(function(j,k) testStat(j,k,ceiling((j+k)/2),z=z,var.known=TRUE), 
                      bonferroniIntervals[,1], bonferroniIntervals[,2])
  } else if(!three.split) { 
    testvec <- mapply(function(j,k) testStat(j,k,ceiling((j+k)/2),z=z,z2=z2, var.known = FALSE), 
                      bonferroniIntervals[,1], bonferroniIntervals[,2])
  } else { ## Using maximum of the three splits
    testvec <- mapply(function(j,k) max(testStat(j,k,ceiling((j+k)/2),z=z,z2=z2, var.known = FALSE),
                                        testStat(j,k,ceiling((3*j+k)/4),z=z,z2=z2, var.known = FALSE),
                                        testStat(j,k,ceiling((j+3*k)/4),z=z,z2=z2, var.known = FALSE)),
                      bonferroniIntervals[,1], bonferroniIntervals[,2])
  }
  critical_val = getCriticalVal(bonferroniIntervals,
                                blocksize, alpha = alpha, 
                                var.known = var.known, three.split = three.split)
  if(var.known) {
    select.idx <- (testvec > critical_val[bonferroniIntervals$block]) &
      !is.nan(testvec)
  } else {
    select.idx <- (testvec > critical_val) & !is.nan(testvec)
  }
  cbind(bonferroniIntervals[select.idx,], testvec[select.idx])
}

## Convert rejected intervals in the form of (s, e] to [s+1, e-1] 
rejectedToConf <- function(intervals) {
  intervals[,1] <- intervals[,1] + 1
  intervals[,2] <- intervals[,2] - 1
  return(intervals)
}

## Get the indices of changepoints from vector of means
getCptIdx <- function(mu) {
  which(diff(mu) != 0)
}

## Get the number of changepoints from vector of means
getNumCpt <- function(mu) {
  length(which(diff(mu) != 0))
}


## Given group of confidence intervals and true changepoints,
## Return TRUE if for each confidence interval, there is at least 
## one changepoint that is covered by the confidence interval and 
## FALSE otherwise
checkCI <- function(ci_results, cpts) {
  if(length(cpts) == 0) {
    return(TRUE)
  }
  tmp <- rep(FALSE, nrow(ci_results$minimal))
  for(i in (1:length(cpts))) {
    tmp <- (tmp | (ci_results$minimal[,1] <= cpts[i] & ci_results$minimal[,2] >= cpts[i]))
  }
  return(all(tmp))
}


lbd <- function(x, var.known = FALSE, alpha = 0.05, even = FALSE) {
  n <- length(x)
  interval_res <- BonferroniIntervals(n, var.known = var.known, even = even)
  z <- c(0,cumsum(x))
  z2 <- c(0, cumsum(x^2))
  rejected_intervals <- getRejectedIntervals(z = z, 
                                             z2 = z2, 
                                             bonferroniIntervals = interval_res$intervals,
                                             blocksize = interval_res$blockSize,
                                             alpha = alpha, 
                                             var.known = var.known)
  conf_intervals <- rejectedToConf(rejected_intervals)
  ci_results <- getCIs(conf_intervals)
  return(ci_results)
}


# Get a set of Maximum Disjoint Intervals given rejected intervals
# Input: 
# intervals: 1st column: lower, 2nd column: upper
getMaximumDisjoint <- function(intervals) {
  if(nrow(intervals) == 0) {
    return(intervals)
  }
  sorted <- intervals[order(intervals[,2], -intervals[,1]),]
  upper_bd <- -Inf
  selected <- rep(FALSE, nrow(sorted))
  for(i in 1:nrow(sorted)) {
    if(sorted[i,1] >= upper_bd) {
      selected[i] <- TRUE
      upper_bd <- sorted[i,2]
    }
  }
  sorted[selected,]
}

## intervals: m by 2 matrix where the first column is the lower bound and the second column is the upper bound 
##           The intervals are closed intervals. 
getCIs <- function(intervals) {
  if(nrow(intervals) == 0) {
    return(list(minimal = intervals,
                disjoint = intervals,
                rejected = intervals))
  }
  sorted <- intervals[order(intervals[,2], -intervals[,1]),]
  fb <- -Inf
  gb <- -Inf
  hb <- - Inf
  minimal_selected <- rep(FALSE, nrow(sorted))
  disjoint_selected <- rep(FALSE, nrow(sorted))
  for(i in 1:nrow(sorted)) {
    if(sorted[i,1] > fb) {
      disjoint_selected[i] <- TRUE
      fb <- sorted[i,2]
    }
    if((sorted[i,1] > gb) & (sorted[i,2] > hb)) {
      minimal_selected[i] <- TRUE
      gb <- sorted[i,1]
      hb <- sorted[i,2]
    }
  }
  return(list(minimal = sorted[minimal_selected,],
              disjoint = sorted[disjoint_selected,],
              rejected = sorted))
}

## plot intervals
plot_cis <- function(x, intervals, xlab = "", ylab = "", title = "") {
  if(nrow(intervals) == 0) return()
  df = data.frame(x = 1:length(x),
                  y = x)
  plt <- ggplot(data = df, aes(x = x, 
                                       y = y)) +
    geom_point(alpha = 0.4) +
    theme_classic() +
    theme(legend.position = "none") + 
    xlab(xlab) +
    ylab(ylab) +
    ggtitle(title) 
  if(nrow(intervals) >= 1) {
    plt <- plt + 
      geom_vline(dat = intervals[seq(1,nrow(intervals), 2),],
                 aes(xintercept = c(lower)),
                 col = "darkblue",
                 linetype = "longdash",
                 linewidth = 0.3) + 
      geom_vline(dat = intervals[seq(1,nrow(intervals), 2),],
                 aes(xintercept = c(upper)),
                 col = "darkblue",
                 linetype = "longdash",
                 linewidth = 0.3) 
  }
  if(nrow(intervals) >= 2) {
    plt <- plt + geom_vline(dat = intervals[seq(2,nrow(intervals), 2),],
                            aes(xintercept = c(lower)),
                            col = "darkorange",
                            linetype = "longdash",
                            linewidth = 0.3) + 
      geom_vline(dat = intervals[seq(2,nrow(intervals), 2),],
                 aes(xintercept = c(upper)),
                 col = "darkorange",
                 linetype = "longdash",
                 linewidth = 0.3)
  }
  plt
}


# z = (0, cumsum(x))
# disjoint.intervals output from function getMaximumDisjoint
# Returns a list of size 2
# result: matrix of estimated interval and mean value in that interval
# changepoints: vector of changepoint estimates
changepointEstimate <- function(z, disjoint.intervals) {
  n.cpt <- nrow(disjoint.intervals)
  if(n.cpt == 0) {
    result <- matrix(c(0, length(z)-1, z[length(z)]/(length(z)-1)), 
                     ncol = 3)
    colnames(result) <- c("lower", "upper", "mean")
    list(result = result, 
         changepoints = numeric(0))
  } else {
    cpt.list <- ceiling((disjoint.intervals$lower + disjoint.intervals$upper)/2)
    cpt.list <- unique(sort(cpt.list))
    result <- matrix(nrow = length(cpt.list) + 1,
                     ncol = 3)
    colnames(result) <- c("lower", "upper", "mean")
    result <- data.frame(result)
    result$lower <- c(0, cpt.list)
    result$upper <- c(cpt.list, n)
    result$mean <- (z[result$upper+1]-z[result$lower+1])/(result$upper-result$lower)
    list(result = result,
         changepoints = cpt.list)
  }
}


