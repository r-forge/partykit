mob2_control <- function(
  alpha = 0.05,
  minsplit = 20L, 
  minbucket = 20L, 
  minprob = 0.01, 
  stump = FALSE, 
  mtry = Inf, 
  maxdepth = Inf, 
  multiway = FALSE, 
  splittry = 2L, 
  MIA = FALSE, 
  maxsurrogate = 0L, 
  numsurrogate = FALSE,
  majority = TRUE, 
  caseweights = TRUE, 
  applyfun = NULL, 
  testflavour = "mfluc", 
  splitflavour = "exhaustive",
  nmax = Inf,
  breakties = FALSE,
  bonferroni = TRUE,
  testtype = "Bonferroni",
  intersplit = FALSE,
  lookahead = FALSE
) {
  
  if((testtype == "Bonferroni") != (bonferroni)) ## bonferroni actually not need
    stop("Arguments bonferroni and testtype must align.")
  
  c(.urp_control(criterion = "p.value",
                 logmincriterion = log(1-alpha), minsplit = minsplit, 
                 minbucket = minbucket, minprob = minprob, stump = stump, 
                 mtry = mtry, maxdepth = maxdepth, multiway = multiway, 
                 splittry = splittry, MIA = MIA, maxsurrogate = maxsurrogate, 
                 numsurrogate = numsurrogate,
                 majority = majority, caseweights = caseweights, 
                 applyfun = applyfun, testflavour = testflavour, 
                 splitflavour = splitflavour),
    list(nmax = nmax, breakties = breakties, testtype = testtype, 
         intersplit = intersplit, lookahead = lookahead)
    # list(teststat = teststat, splitstat = splitstat, splittest = splittest, pargs = pargs,
    #      testtype = testtype, nmax = nmax, nresample = nresample, lookahead = lookahead,
    #      intersplit = intersplit)
    )
}

# ## control splitting parameters
# mob2_control <- function(
#   alpha = 0.05, bonferroni = TRUE, minsize = NULL, maxdepth = Inf,
#   mtry = Inf, trim = 0.1, breakties = FALSE, parm = NULL, dfsplit = TRUE, prune = NULL, restart = TRUE,
#   verbose = FALSE, caseweights = TRUE, ytype = "vector", xtype = "matrix",
#   terminal = "object", inner = terminal, model = TRUE,
#   numsplit = "left", catsplit = "binary", vcov = "opg", ordinal = "chisq", nrep = 10000,
#   minsplit = minsize, minbucket = minsize,
#   applyfun = NULL, cores = NULL,
#   ## NEW
#   testflavour = "mfluc", 
#   splitflavour = "exhaustive",
#   pargs = GenzBretz(),
#   mincriterion = 1 - alpha,
#   logmincriterion = log(mincriterion),
#   minprob = 0.01,
#   stump = FALSE,
#   lookahead = FALSE,
#   MIA = FALSE,
#   maxsurrogate = 0L,
#   numsurrogate = FALSE,
#   splittry = 2L,
#   intersplit = FALSE,
#   majority = FALSE
# )
# {
#   
#   if (!caseweights)
#     stop("only caseweights currently implemented")
#   
#   ## no mtry if infinite or non-positive
#   if(is.finite(mtry)) {
#     mtry <- if(mtry < 1L) Inf else as.integer(mtry)
#   }
#   
#   ## data types for formula processing
#   ytype <- match.arg(ytype, c("vector", "data.frame", "matrix"))
#   xtype <- match.arg(xtype, c("data.frame", "matrix"))
#   
#   ## what to store in inner/terminal nodes
#   if(!is.null(terminal)) terminal <- as.vector(sapply(terminal, match.arg, c("estfun", "object")))
#   if(!is.null(inner))    inner    <- as.vector(sapply(inner,    match.arg, c("estfun", "object")))
#   
#   ## how to split and how to select splitting variables
#   numsplit <- match.arg(tolower(numsplit), c("left", "center", "centre"))
#   if(numsplit == "centre") numsplit <- "center"
#   catsplit <- match.arg(tolower(catsplit), c("binary", "multiway"))
#   vcov <- match.arg(tolower(vcov), c("opg", "info", "sandwich"))
#   ordinal <- match.arg(tolower(ordinal), c("l2", "max", "chisq"))
#   
#   ## apply infrastructure for determining split points
#   if(is.null(applyfun)) {
#     applyfun <- if(is.null(cores)) {
#       lapply
#     } else {
#       function(X, FUN, ...) parallel::mclapply(X, FUN, ..., mc.cores = cores)
#     }
#   }
#   
#   ## return list with all options
#   c(.urp_control(criterion = ifelse(testtype == "Teststatistic", 
#                                     "statistic", "p.value"),
#                  logmincriterion = logmincriterion, minsplit = minsplit, 
#                  minbucket = minbucket, minprob = minprob, stump = stump, 
#                  mtry = mtry, maxdepth = maxdepth, multiway = multiway, 
#                  splittry = splittry, MIA = MIA, maxsurrogate = maxsurrogate, 
#                  numsurrogate = numsurrogate,
#                  majority = majority, caseweights = caseweights, 
#                  applyfun = applyfun, testflavour = testflavour, 
#                  splitflavour = splitflavour),
#     list(
#       lookahead = lookahead, intersplit = intersplit))
#   return(rval)
# }


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
  
  if (ORDERED) {
    sp <- NULL
    maxlogLik <- -Inf
    linfo <- rinfo <- info
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
  
  if (!splitonly)
    return(statistic = maxlogLik, p.value = NA)
  
  if (all(is.na(sp))) return(NULL)
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
                      index = as.integer(sp) + 1L)
  }
  return(ret)
}

