.modeltrafo <- function(formula, data, ctrl, converged = NULL, fit, ...) {
  
  weights <- model.weights(data)
  if (is.null(weights)) weights <- integer(0)
  cluster <- data[["cluster"]]
  offset <- model.offset(data)
  
  ### <FIXME> handle offset and cluster </FIXME>
  if (!is.null(cluster)) stop("cluster not implemented")
  mf <- model.frame(formula, data, na.action = na.pass)
  cc <- complete.cases(mf)
  y <- model.response(mf)
  x <- model.matrix(formula, data = mf)
  
  ## function for model fitting
  modelfit <- function(subset, estfun = TRUE, object = FALSE, info = NULL, ...) {
    
    ## get subset of the data
    s <- subset[cc[subset]]
    ys <- y[s]
    xs <- x[s, , drop = FALSE]
    nobs <- NROW(xs)
    
    if (length(weights) > 0) {
      weights <- weights[cc[subset]]
      nobs <- sum(weights)
    } else {
      weights <- NULL
    }
    
    ## call the fit function
    args <- c(list(x = xs, y = ys, start = info$coef, weights = weights),
              list(object = TRUE, estfun = TRUE)[c("estfun", "object") %in% 
                                                   names(formals(fit))],
              ...)
    
    ret <- do.call("fit", args = args)
    
    
    ## if ret is not a list of object, estfun, ...
    if(class(ret)[1] != "list") {
      ret <- list(estfun = estfun(ret), coefficients = coef(ret), 
                  objfun = - logLik(ret), object = ret)
    }
    
    
    
    ## get convergence info
    if (is.null(converged)) {
      cv <- if (is.null(ret$object$converged)) TRUE else ret$object$converged
    } else {
      cv <- converged(ret$object, mf, subset)
    }
    
    ## correct dimension of estfun 
    ef <- matrix(0, nrow = NROW(x), ncol = NCOL(x))
    ef[subset,] <- ret$estfun
    
    
    ## return
    list(estfun = ef, coefficients = ret$coefficients, objfun = - ret$objfun,
      object = if (object) ret$object else NULL, nobs = nobs, 
      converged = cv)
    
  }
  
  return(modelfit)
}



mob2_control <- function(
  alpha = 0.05,
  mincriterion = 1 - alpha,
  minsize = 20L,
  minsplit = minsize, 
  minbucket = minsize, 
  minprob = 0.01, 
  stump = FALSE, 
  mtry = Inf, 
  maxdepth = Inf, 
  # nmax = Inf, # TODO: check if this works first
  splittry = 2L, 
  MIA = FALSE, 
  maxsurrogate = 0L, 
  numsurrogate = FALSE,
  majority = FALSE, 
  caseweights = TRUE, 
  applyfun = NULL, 
  cores = NULL, 
  testflavour = "mfluc", 
  splitflavour = "exhaustive",
  lookahead = FALSE,
  testtype = "Bonferroni",
  bonferroni = TRUE,
  nresample = 9999L,   # used for testtype = "MonteCarlo"
  breakties = FALSE,
  teststat = "quadratic",  # used for testflavour/splitflavour = "ctree"
  splitstat = "quadratic", # used for testflavour/splitflavour = "ctree"
  splittest = FALSE,        # used for testflavour/splitflavour = "ctree"
  numsplit = "left",
  catsplit = "binary",
  trim = 0.1, 
  parm = NULL
) {
  
  if(("Bonferroni" %in% testtype) != (bonferroni)) {
    warning("Arguments bonferroni and testtype must align. 
            Turning Bonferroni adjustment off.")
    bonferroni <- FALSE
    if("Bonferroni" %in% testtype) testtype <- testtype[testtype != "Bonferroni"]
    if(length(testtype) == 0) testtype <- "Univariate"
  }
  
  intersplit <- numsplit == "center"
  multiway <- catsplit == "multiway"
  
  if(testflavour == "exhaustive") { 
    alpha <- 1
    criterion <- "statistic"
  } else {
    criterion <- ifelse("Teststatistic" %in% testtype, 
                        "statistic", "p.value")
  }
  
  
  if("statistic" %in% criterion & mincriterion < 1 & mincriterion >= 0.9)
    warning("When criterion = 'statistic', mincriterion is the test statistic that must be exceeded. 
            Are you sure you chose the correct value for mincriterion?")
  
  
  c(.urp_control(criterion = criterion,
                 logmincriterion = log(mincriterion), minsplit = minsplit, 
                 minbucket = minbucket, minprob = minprob, nmax = Inf, 
                 stump = stump, lookahead = lookahead, mtry = mtry, 
                 maxdepth = maxdepth, multiway = multiway, splittry = splittry, 
                 MIA = MIA, maxsurrogate = maxsurrogate, numsurrogate = numsurrogate, 
                 majority = majority, caseweights = caseweights, 
                 applyfun = applyfun, cores = cores, testflavour = testflavour, 
                 splitflavour = splitflavour, bonferroni = bonferroni, trim = trim),
    list(breakties = breakties, testtype = testtype, nresample = nresample, 
         intersplit = intersplit, teststat = teststat, splitstat = splitstat, 
         splittest = splittest, parm = parm)
    )
}




## variable selection: given model scores, conduct
## all M-fluctuation tests for orderins in z
.fluct_test_split <- function(estfun, z, weights, obj = NULL, cluster = NULL, control)
{  
  
  if(length(weights) == 0) weights <- rep(1, NROW(estfun))
  
  ## set up return values
  m <- NCOL(z)
  pval <- rep.int(NA_real_, m)
  stat <- rep.int(0, m)
  ifac <- rep.int(FALSE, m)
  
  ## estimating functions (dropping zero weight observations)
  process <- as.matrix(estfun)
  ww0 <- (weights > 0)
  process <- process[ww0, , drop = FALSE]
  z <- z[ww0]
  k <- NCOL(process)
  n <- NROW(process)
  
  ## scale process
  process <- process/sqrt(n)
  vcov <- control$vcov
  if(is.null(obj)) vcov <- "opg"
  if(vcov != "opg") {
    bread <- vcov(obj) * n
  }
  if(vcov != "info") {
    meat <- if(is.null(cluster)) {
      crossprod(process)
    } else {
      crossprod(as.matrix(aggregate(process, by = list(cluster), FUN = sum)[, -1L, drop = FALSE]))
    }
  }
  J12 <- strucchange::root.matrix(switch(vcov,
                                         "opg" = chol2inv(chol(meat)),
                                         "info" = bread,
                                         "sandwich" = bread %*% meat %*% bread
  ))
  process <- t(J12 %*% t(process))  
  
  
  ## get critical values for supLM statistic
  from <- if(control$trim > 1) control$trim else ceiling(n * control$trim)
  from <- max(from, control$minbucket)
  to <- n - from
  lambda <- ((n - from) * to)/(from * (n - to))
  
  beta <- mob_beta_suplm
  logp.supLM <- function(x, k, lambda)
  {
    if(k > 40L) {
      ## use Estrella (2003) asymptotic approximation
      logp_estrella2003 <- function(x, k, lambda)
        -lgamma(k/2) + k/2 * log(x/2) - x/2 + log(abs(log(lambda) * (1 - k/x) + 2/x))
      ## FIXME: Estrella only works well for large enough x
      ## hence require x > 1.5 * k for Estrella approximation and
      ## use an ad hoc interpolation for larger p-values
      p <- ifelse(x <= 1.5 * k, (x/(1.5 * k))^sqrt(k) * logp_estrella2003(1.5 * k, k, lambda), logp_estrella2003(x, k, lambda))
    } else {
      ## use Hansen (1997) approximation
      nb <- ncol(beta) - 1L
      if(lambda<1) tau <- lambda
      else tau <- 1/(1 + sqrt(lambda))
      beta <- beta[(((k - 1) * 25 + 1):(k * 25)),]
      dummy <- beta[,(1L:nb)]%*%x^(0:(nb-1))
      dummy <- dummy*(dummy>0)
      pp <- pchisq(dummy, beta[,(nb+1)], lower.tail = FALSE, log.p = TRUE)
      if(tau == 0.5)
        p <- pchisq(x, k, lower.tail = FALSE, log.p = TRUE)
      else if(tau <= 0.01)
        p <- pp[25L]
      else if(tau >= 0.49)
        p <- log((exp(log(0.5 - tau) + pp[1L]) + exp(log(tau - 0.49) + pchisq(x, k, lower.tail = FALSE, log.p = TRUE))) * 100)
      else
      {
        taua <- (0.51 - tau) * 50
        tau1 <- floor(taua)
        p <- log(exp(log(tau1 + 1 - taua) + pp[tau1]) + exp(log(taua-tau1) + pp[tau1 + 1L]))
      }
    }
    return(as.vector(p))
  }
  
  ## compute statistic and p-value
  if(is.factor(z)) {
    proci <- process[order(z), , drop = FALSE]
    ifac <- TRUE
    iord <- is.ordered(z) & (control$ordinal != "chisq")
    
    ## order partitioning variable
    z <- z[order(z)]
    # re-apply factor() added to drop unused levels
    z <- factor(z, levels = unique(z))
    # compute segment weights
    segweights <- as.vector(table(z))/n
    
    # compute statistic only if at least two levels are left
    if(length(segweights) < 2L) {
      stat <- 0
      pval <- NA_real_
    } else if(iord) {
      proci <- apply(proci, 2L, cumsum)
      tt <- head(cumsum(segweights), -1L)
      if(control$ordinal == "max") {
        stat <- max(abs(proci[round(tt * n), ] / sqrt(tt * (1-tt))))
        pval <- log(as.numeric(1 - mvtnorm::pmvnorm(
          lower = -stat, upper = stat,
          mean = rep(0, length(tt)),
          sigma = outer(tt, tt, function(x, y)
            sqrt(pmin(x, y) * (1 - pmax(x, y)) / ((pmax(x, y) * (1 - pmin(x, y))))))
        )^k))
      } else {
        proci <- rowSums(proci^2)
        stat <- max(proci[round(tt * n)] / (tt * (1-tt)))
        pval <- log(strucchange::ordL2BB(segweights, nproc = k, nrep = control$nrep)$computePval(stat, nproc = k))
      }
    } else {      
      stat <- sum(sapply(1L:k, function(j) (tapply(proci[,j], z, sum)^2)/segweights))
      pval <- pchisq(stat, k*(length(levels(z))-1), log.p = TRUE, lower.tail = FALSE)
    }
  } else {
    oi <- if(control$breakties) {
      mm <- sort(unique(z))
      mm <- ifelse(length(mm) > 1L, min(diff(mm))/10, 1)
      order(z + runif(length(z), min = -mm, max = +mm))
    } else {
      order(z)
    }
    proci <- process[oi, , drop = FALSE]
    proci <- apply(proci, 2L, cumsum)
    stat <- if(from < to) {
      xx <- rowSums(proci^2)
      xx <- xx[from:to]
      tt <- (from:to)/n
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

.objfun_test_split <- function(trafo, info = info, x, bdr = NULL, j, ctrl, subset, 
                               weights, cluster, splitonly = TRUE, 
                               minbucket) {
  
  if (all(is.na(x[subset]))) { ### all x values are missing
    if (splitonly) return(NULL)
    return(list(statistic = NA, p.value = NA))
  }
  
  if (is.null(cluster)) cluster <- integer(0)
  if (splitonly) {
    if ((ctrl$multiway && ctrl$maxsurrogate == 0) &&
        is.factor(x) && nlevels(x[subset, drop = TRUE]) > 1) {
      index <- 1L:nlevels(x)
      if (length(weights) > 0) {
        xt <- xtabs(weights ~ x, subset = subset)
      } else {
        xt <- xtabs(~ x, subset = subset)
      }
      index[xt == 0] <- NA
      index[xt > 0 & xt < minbucket] <- nlevels(x) + 1L
      if (length(unique(index)) == 1) return(NULL)
      index <- unclass(factor(index))
      return(partysplit(as.integer(j),
                        index = as.integer(index)))
    }
  }
  
  ux <- NULL
  ORDERED <- is.ordered(x) || is.numeric(x)
  if (is.null(bdr)) {
    if (!is.factor(x)) {
      x[-subset] <- NA
      ux <- sort(unique(x))
    }
  } else {
    ix <- bdr[[j]]
    ux <- attr(ix, "levels")
  }
  
  linfo <- rinfo <- info
  maxlogLik <- nosplitll <- trafo(subset, info = info, estfun = FALSE)$objfun
  sp <- NULL
  
  if (ORDERED) {
    for (u in 1:length(ux)) {
      sleft <- subset[LEFT <- x[subset] <= ux[u]]
      sright <- subset[!LEFT]
      if (length(weights) > 0) {
        if (sum(weights[sleft]) < minbucket ||
            sum(weights[sright]) < minbucket)
          next();
      } else {
        if (length(sleft) < minbucket || 
            length(sright) < minbucket)
          next();
      }
      ltr <- trafo(sleft, info = linfo, estfun = FALSE)
      rtr <- trafo(sright, info = rinfo, estfun = FALSE)
      ll <- ltr$objfun + rtr$objfun
      linfo <- ltr$info
      rinfo <- rtr$info
      if (ll > maxlogLik) {
        sp <- u
        maxlogLik <- ll
      }
    }
  } else {
    splits <- mob_grow_getlevels(x)
    for (u in 1:nrow(splits)) {
      sleft <- subset[LEFT <- x[subset] %in% levels(x)[splits[u,]]]
      sright <- subset[!LEFT]
      if (length(weights) > 0) {
        if (sum(weights[sleft]) < minbucket ||
            sum(weights[sright]) < minbucket)
          next();
      } else {
        if (length(sleft) < minbucket || 
            length(sright) < minbucket)
          next();
      }
      ltr <- trafo(sleft, info = linfo, estfun = FALSE)
      rtr <- trafo(sright, info = rinfo, estfun = FALSE)
      ll <- ltr$objfun + rtr$objfun
      linfo <- ltr$info
      rinfo <- rtr$info
      if (ll > maxlogLik) {
        sp <- splits[u,] + 1L
        maxlogLik <- ll
      }
    }
  }
  
  if (!splitonly){
    ## split only if logLik improves due to splitting
    maxlogLik <- ifelse(maxlogLik == nosplitll, NA, maxlogLik)
    return(list(statistic = maxlogLik, p.value = NA))
  }
  if (is.null(sp) || all(is.na(sp))) return(NULL)
  if (ORDERED) {
    if (!is.ordered(x))
      ### interpolate split-points, see https://arxiv.org/abs/1611.04561
      if (ctrl$intersplit & sp < length(ux)) {
        sp <- (ux[sp] + ux[sp + 1]) / 2 
      } else {
        sp <- ux[sp]  ### x <= sp vs. x > sp
      }
    ret <- partysplit(as.integer(j), breaks = sp,
                      index = 1L:2L)
  } else {
    ret <- partysplit(as.integer(j),
                      index = as.integer(sp))
  }
  return(ret)
}

