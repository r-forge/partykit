## variable selection: given model scores, conduct
## all M-fluctuation tests for orderins in z
.mfluc_test <- function(model, trafo, data, subset, weights, j, SPLITONLY, ctrl)
{

  stopifnot(!SPLITONLY)
  stopifnot(is.null(.get_index(data, "yx")))
  z <- .get_var(data, j)[subset]
  estfun <- model$estfun[subset]
  cluster <- .get_var(data, "(cluster)")[subset]
  if(length(weights) == 0) {
      weights <- rep(1, NROW(estfun))
  } else {
      weights <- weights[subset]
  }
  obj <- model$obj

  if(length(unique(z)) < 2L) return(list(statistic = NA, p.value = NA))
  
  ## set up return values
  m <- NCOL(z)
  pval <- rep.int(NA_real_, m)
  stat <- rep.int(0, m)
  ifac <- rep.int(FALSE, m)
  
  ## estimating functions (dropping zero weight observations)
  process <- as.matrix(estfun)
  ww0 <- (weights > 0)
  ww0[!(seq_along(ww0) %in% subset)] <- FALSE
  process <- process[ww0, , drop = FALSE]
  if(!is.null(cluster)) cluster <- droplevels(cluster[ww0]) 
  weights <- weights[ww0]
  z <- z[ww0]
  k <- NCOL(process)
  n <- NROW(process)
  nobs <- if(ctrl$caseweights && any(weights != 1L)) sum(weights) else n

  
  ## stop if all collumn values in process are the same
  if(all(apply(process, 2, function(x) length(unique(x))) == 1)) 
    return(list(statistic = NA, p.value = NA))
  
  ## scale process
  process <- process/sqrt(nobs)
  vcov <- ctrl$vcov
  if(is.null(obj)) vcov <- "opg"
  if(vcov != "opg") {
    bread <- vcov(obj) * nobs
  }
  if(vcov != "info") {
    meat <- if(is.null(cluster)) {
      crossprod(process/sqrt(weights))
    } else {
      crossprod(as.matrix(apply(process/sqrt(weights), 2L, tapply, cluster, sum)))
    }
  }
  ## from strucchange
  root.matrix <- function(X) {
    if((ncol(X) == 1L)&&(nrow(X) == 1L)) return(sqrt(X)) else {
      X.eigen <- eigen(X, symmetric = TRUE)
      if(any(X.eigen$values < 0)) stop("Matrix is not positive semidefinite")
      sqomega <- sqrt(diag(X.eigen$values))
      V <- X.eigen$vectors
      return(V %*% sqomega %*% t(V))
    }
  }
  J12 <- root.matrix(switch(vcov,
         		    "opg" = chol2inv(chol(meat)),
         		    "info" = bread,
         		    "sandwich" = bread %*% meat %*% bread
  ))
  process <- t(J12 %*% t(process))  
  
  
  ## get critical values for supLM statistic
  from <- if(ctrl$trim > 1) ctrl$trim else ceiling(nobs * ctrl$trim)
  from <- max(from, ctrl$minbucket)
  to <- nobs - from
  lambda <- ((nobs - from) * to)/(from * (nobs - to))
  
  beta <- partykit:::mob_beta_suplm
  logp.supLM <- function(x, k, lambda)
  {
    if(k > 40L) {
      ## use Estrella (2003) asymptotic approximation
      logp_estrella2003 <- function(x, k, lambda)
        -lgamma(k/2) + k/2 * log(x/2) - x/2 + log(abs(log(lambda) * 
         (1 - k/x) + 2/x))
      ## FIXME: Estrella only works well for large enough x
      ## hence require x > 1.5 * k for Estrella approximation and
      ## use an ad hoc interpolation for larger p-values
      p <- ifelse(x <= 1.5 * k, 
        (x/(1.5 * k))^sqrt(k) * logp_estrella2003(1.5 * k, k, lambda), 
        logp_estrella2003(x, k, lambda))
    } else {
      ## use Hansen (1997) approximation
      nb <- ncol(beta) - 1L
      tau <- if(lambda < 1) lambda else 1/(1 + sqrt(lambda))
      beta <- beta[(((k - 1) * 25 + 1):(k * 25)),]
      dummy <- beta[,(1L:nb)] %*% x^(0:(nb-1))
      dummy <- dummy * (dummy > 0)
      pp <- pchisq(dummy, beta[,(nb+1)], lower.tail = FALSE, log.p = TRUE)
      if(tau == 0.5) {
        p <- pchisq(x, k, lower.tail = FALSE, log.p = TRUE)
      } else if(tau <= 0.01) {
        p <- pp[25L]
      } else if(tau >= 0.49) {
        p <- log((exp(log(0.5 - tau) + pp[1L]) + exp(log(tau - 0.49) + pchisq(x, k, lower.tail = FALSE, log.p = TRUE))) * 100)
	## if p becomes so small that 'correct' weighted averaging does not work, resort to 'naive' averaging
	if(!is.finite(p)) p <- mean(c(pp[1L], pchisq(x, k, lower.tail = FALSE, log.p = TRUE)))
      } else {
        taua <- (0.51 - tau) * 50
        tau1 <- floor(taua)
  	p <- log(exp(log(tau1 + 1 - taua) + pp[tau1]) + exp(log(taua-tau1) + pp[tau1 + 1L]))
	## if p becomes so small that 'correct' weighted averaging does not work, resort to 'naive' averaging
	if(!is.finite(p)) p <- mean(pp[tau1 + 0L:1L])
      }
    }
    return(as.vector(p))
  }
  
  ## compute statistic and p-value
  if(is.factor(z)) {
    oi <- order(z)
    proci <- process[oi, , drop = FALSE]
       
    ifac <- TRUE
    iord <- is.ordered(z) & (ctrl$ordinal != "chisq")
    
    ## order partitioning variable
    z <- z[oi]
    # re-apply factor() added to drop unused levels
    z <- factor(z, levels = unique(z))
    # compute segment weights
    segweights <- if(ctrl$caseweights) tapply(weights[oi], z, sum) else table(z)
    segweights <- as.vector(segweights)/nobs
    
    # compute statistic only if at least two levels are left
    if(length(segweights) < 2L) {
      stat <- 0
      pval <- NA_real_
    } else if(iord) {
      proci <- apply(proci, 2L, cumsum)
      tt0 <- head(cumsum(table(z)), -1L)
      tt <- head(cumsum(segweights), -1L)
      if(ctrl$ordinal == "max") {
        stat <- max(abs(proci[tt0, ] / sqrt(tt * (1-tt))))
        pval <- log(as.numeric(1 - mvtnorm::pmvnorm(
          lower = -stat, upper = stat,
          mean = rep(0, length(tt)),
          sigma = outer(tt, tt, function(x, y)
            sqrt(pmin(x, y) * (1 - pmax(x, y)) / 
            ((pmax(x, y) * (1 - pmin(x, y))))))
        )^k))
      } else {
        proci <- rowSums(proci^2)
        stat <- max(proci[tt0] / (tt * (1-tt)))
        pval <- log(strucchange::ordL2BB(segweights, nproc = k, nrep = ctrl$nrep)$computePval(stat, nproc = k))
      }
    } else {      
      stat <- sum(sapply(1L:k, function(j) (tapply(proci[,j], z, sum)^2)/segweights))
      pval <- pchisq(stat, k*(length(levels(z))-1), log.p = TRUE, lower.tail = FALSE)
    }
  } else {
    oi <- if(ctrl$breakties) {
      mm <- sort(unique(z))
      mm <- ifelse(length(mm) > 1L, min(diff(mm))/10, 1)
      order(z + runif(length(z), min = -mm, max = +mm))
    } else {
      order(z)
    }
    proci <- process[oi, , drop = FALSE]
    proci <- apply(proci, 2L, cumsum)
    tt0 <- if(ctrl$caseweights && any(weights != 1L)) cumsum(weights[oi]) else 1:n
    stat <- if(from < to) {
      xx <- rowSums(proci^2)
      xx <- xx[tt0 >= from & tt0 <= to]
      tt <- tt0[tt0 >= from & tt0 <= to]/nobs
      max(xx/(tt * (1 - tt)))	  
    } else {
      0
    }
    pval <- if(from < to) logp.supLM(stat, k, lambda) else NA
  }
  
  ## return version of pvalue that urp deals with
  rval <- list(statistic = log(stat), p.value = log1p(-exp(pval)))
  return(rval)
}

