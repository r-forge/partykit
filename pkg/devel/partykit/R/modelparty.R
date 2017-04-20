
.modeltrafo <- function(
  formula, 
  data, 
  ctrl, 
  converged = NULL, 
  fit, 
  ...
){
  
  weights <- model.weights(data)
  if (is.null(weights)) weights <- integer(0)
  cluster <- data[["cluster"]]
  offset <- model.offset(data)
  
  ### <FIXME> handle offset and cluster </FIXME>
  if (!is.null(cluster)) stop("cluster not implemented")
  mf <- model.frame(formula, data, na.action = na.pass)
  y <- model.response(mf)
  x <- model.matrix(formula, data = mf)
  cc <- complete.cases(x)
  
  ## function for model fitting
  modelfit <- function(subset, estfun = TRUE, object = "object" %in% ctrl$inner, 
                       info = NULL) {
    
    ## get subset of the data
    s <- subset[cc[subset]]
    if(is.matrix(y)) {
      ys <- y[s, , drop = FALSE]
    } else {
      ys <- y[s]
    }
    xs <- x[s, , drop = FALSE]
    attr(xs, "formula") <- formula
    attr(xs, "terms") <- terms(formula, data = data[s, ])
    attr(xs, "offset") <- offset
    
    nobs <- NROW(xs)
    
    if (length(weights) > 0) {
      weights <- weights[s]
      if (ctrl$caseweights) nobs <- sum(weights) else nobs <- sum(weights > 0)
    } else {
      weights <- NULL
    }
    
    if (length(offset) > 0) {
      offset <- offset[s]
    } else {
      offset <- NULL
    }
    
    
    ## call the fit function
    args <- c(list(x = if (ncol(xs) == 0) NULL else xs, y = ys, 
                   start = info$coef, weights = weights, offset = offset),
              list(object = TRUE, estfun = TRUE)[c("estfun", "object") %in% 
                                                   names(formals(fit))],
              list(...))

    ret <- do.call("fit", args = args)
    
    
    ## if ret is not a list of object, estfun, ...
    if(class(ret)[1] != "list") {
      ret <- list(estfun = sandwich::estfun(ret), coefficients = coef(ret), 
                  objfun = - logLik(ret), object = ret)
    }
    
    
    
    ## get convergence info
    if (is.null(converged)) {
      cv <- if (is.null(ret$object$converged)) TRUE else ret$object$converged
    } else {
      cv <- converged(ret$object, mf, subset)
    }
    
    ## correct dimension of estfun 
    ef <- NULL
    if(estfun) {
      ef <- matrix(0, nrow = NROW(x), ncol = NCOL(ret$estfun))
      ef[s,] <- ret$estfun
      if(!is.null(ctrl$parm)) ef <- ef[, ctrl$parm]
    }
    
    ## return
    list(estfun = ef, coefficients = ret$coefficients, objfun = - ret$objfun,
         object = if (object) ret$object else NULL, nobs = nobs, 
         converged = cv)
    
  }
  
  return(modelfit)
}



mob_control <- function(
  testflavour = "mfluc",
  testtype = "Bonferroni",
  alpha = 0.05,
  mincriterion = 1 - alpha,
  parm = NULL,
  bonferroni = TRUE,
  breakties = FALSE,
  nrep = 10000L,
  ordinal = "chisq",
  trim = 0.1,
  nresample = 9999L,
  teststat = "quadratic",
  splittry = 2,
  vcov = "opg",
  splitflavour = "exhaustive",
  splittest = FALSE,
  splitstat = "quadratic",
  catsplit = "binary",
  numsplit = "left",
  majority = FALSE,
  maxsurrogate = 0,
  MIA = FALSE,
  numsurrogate = FALSE,
  lookahead = FALSE,
  maxdepth = Inf,
  minbucket = minsize,
  minprob = 0.01,
  minsplit = minsize,
  minsize = NULL,
  stump = FALSE,
  caseweights = TRUE,
  dfsplit = TRUE,
  applyfun = NULL,
  cores = NULL,
  restart = TRUE, ### TODO: change default to FALSE
  inner = "object",
  model = TRUE,
  terminal = "object",
  mtry = Inf, 
  # deprecated
  verbose, xtype, ytype, prune
) {
  
  if (!missing("verbose"))
    warning("argument verbose deprecated")
  if (!missing("xtype"))
    warning("argument xtype deprecated")
  if (!missing("ytype"))
    warning("argument ytype deprecated")
  if (!missing("prune"))
    warning("argument prune deprecated")
    
  
  ## FIXME: Involve caseweights in stopping (minbucket, minsize, etc)
  ##        For now only used to compute correct n (e.g. for printing)
  
  if("estfun" %in% inner) {
    inner <- inner[inner != "estfun"]
    warning("estfun can no longer be stored in inner nodes")
    if(length(inner) == 0) inner <- NULL
  }
  
  if(("Bonferroni" %in% testtype) != (bonferroni)) {
    bonferroni <- FALSE
    if("Bonferroni" %in% testtype) 
      testtype <- testtype[testtype != "Bonferroni"]
    if(length(testtype) == 0) testtype <- "Univariate"
    if(testflavour == "ctree") 
      warning("Arguments bonferroni and testtype must align. 
            Turning Bonferroni adjustment off.")
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
  
  
  if ("statistic" %in% criterion) { 
    if (mincriterion < 1 & mincriterion >= 0.9)
      warning("When criterion = 'statistic', mincriterion is the test 
               statistic that must be exceeded. 
               Are you sure you chose the correct value for mincriterion?")
  } else {
    if (mincriterion <= 0 | mincriterion >= 1)
      stop("mincriterion = 1 - alpha must be between 0 and 1.")
  }
  
  ### make sure right criterion is used for exhaustive search
  if(testflavour == "exhaustive"){
    criterion <- "statistic"
    logmincriterion <- -Inf
  }
  
  c(.urp_control(criterion = criterion,
                 logmincriterion = log(mincriterion), minsplit = minsplit, 
                 minbucket = minbucket, minprob = minprob, nmax = Inf, 
                 stump = stump, lookahead = lookahead, mtry = mtry, 
                 maxdepth = maxdepth, multiway = multiway, splittry = splittry, 
                 MIA = MIA, maxsurrogate = maxsurrogate, 
                 numsurrogate = numsurrogate, majority = majority, 
                 caseweights = caseweights, applyfun = applyfun, cores = cores, 
                 saveinfo = TRUE, ### always
                 testflavour = testflavour, splitflavour = splitflavour, 
                 bonferroni = bonferroni, trim = trim),
    list(breakties = breakties, testtype = testtype, nresample = nresample, 
         intersplit = intersplit, teststat = teststat, splitstat = splitstat, 
         splittest = splittest, parm = parm, dfsplit = dfsplit, 
         restart = restart, model = model, vcov = vcov, ordinal = ordinal, 
         nrep = nrep, terminal = terminal, inner = inner)
  )
}

mob <- function
(
  formula, 
  data, 
  subset,
  na.action = na.pass, 
  weights = NULL, 
  offset,
  cluster, 
  fit,
  control = mob_control(...), 
  converged = NULL,
  ...
) {
  ### get the call and the calling environment for .urp_tree
  call <- match.call(expand.dots = FALSE)
  call$na.action <- na.action
  frame <- parent.frame()
  if (missing(data)) {
    data <- NULL
    data_asis <- FALSE
  } else {
    data_asis <- missing(weights) && missing(subset) && 
      missing(cluster) && missing(offset)
  }
  
  Formula <- as.Formula(formula)
  weights <- model.weights(data)
  mf <- model.frame(Formula, data, na.action = na.pass)
  y <- model.response(mf)
  if(length(Formula)[2] == 2) {
    modelfmla <- formula(Formula, lhs = 1L, rhs = 1L)
  } else {
    modelfmla <- formula(Formula, lhs = 1L, rhs = 0L)
  }
  x <- model.matrix(modelfmla, data = mf)
  
  ### if minsize is NULL, set to 10 * number of parameters
  if (is.null(control$minbucket) | is.null(control$minsplit)) {
    cluster <- data[["cluster"]]
    offset <- model.offset(data)
    cc <- complete.cases(mf)
    
    args <- c(list(x = if (ncol(x) == 0) NULL else x, y = y, 
                   weights = weights, offset = offset),
              list(object = FALSE, estfun = FALSE)[c("estfun", "object") %in% 
                                                     names(formals(fit))],
              list(...))
    mod0 <- do.call("fit", args = args)
    
    n_coef <- ifelse(is.list(mod0), length(mod0$coefficients), coef(mod0))
    if(n_coef == 0) n_coef <- 1
    n_y <- NCOL(y)
    minsize <- as.integer(ceiling(10L * n_coef/n_y))
    if (is.null(control$minbucket)) control$minbucket <- minsize
    if (is.null(control$minsplit)) control$minsplit <- minsize
  }
  

  
  # trafofun <- function(...) .modeltrafo(..., converged = converged, fit = fit)
  dots <- list(...)
  trafofun <- function(...) 
    do.call(".modeltrafo", c(list(...), dots, converged = converged, fit = fit))
  tree <- .urp_tree(call, frame, data = data, data_asis = data_asis, 
                    control = control, trafofun = trafofun, doFit = TRUE)
  
  ### prepare as modelparty
  mf <- tree$mf
  # weights <- model.weights(mf)
  if (length(weights) == 0) weights <- rep(1, nrow(mf))
  modelf <- as.Formula(tree$modelf)
  partf <- as.Formula(tree$partf)
  mtY <- terms(modelf, data = mf)
  mtZ <- delete.response(terms(partf, data = mf))
  
  fitted <- data.frame("(fitted)" = fitted_node(tree$nodes, mf),
                       "(weights)" = weights,
                       check.names = FALSE)
  mmf <- model.frame(modelf, data = mf)
  # y <- model.part(modelf, data = mmf, lhs = 1, rhs = 0)
  # if (length(y) == 1) y <- y[[1]]
  fitted[[3]] <- y
  names(fitted)[3] <- "(response)"
  
  control$ytype <- ifelse(is.vector(y), "vector", class(y))
  # x <- model.matrix(modelf, data = mmf)
  control$xtype <- "matrix" # TODO: find out when to use data.frame
  
  ## return party object
  rval <- party(tree$nodes, 
                data = if(control$model) mf else mf[0,],
                fitted = fitted,
                terms = tree$terms,
                info = list(
                  call = match.call(),
                  formula = formula,
                  Formula = as.Formula(formula),
                  terms = list(response = mtY, partitioning = mtZ),
                  fit = fit,
                  control = control,
                  dots = list(...),
                  nreg = NCOL(x)
                )
  )
  class(rval) <- c("modelparty", class(rval))
  
  ### add modelinfo (object) and estfun if not there yet, but wanted
  # TODO: check if this can be done prettier
  which_terminals <- nodeids(rval, terminal = TRUE)
  which_all <- nodeids(rval)
  
  idx <- lapply(which_all, .get_path, obj = tree$nodes)
  names(idx) <- which_all
  tree_ret <- unclass(rval)
  subset_term <- predict(rval, type = "node")
  
  for (i in which_all) {
    ichar <- as.character(i)
    iinfo <- tree_ret[[c(1, idx[[ichar]])]]$info
    
    if (i %in% which_terminals) winfo <- control$terminal else 
      winfo <- control$inner
    
    if (is.null(winfo)) {
      iinfo$object <- NULL
      iinfo$estfun <- NULL
    } else {
      if (is.null(iinfo) | any(is.null(iinfo[[winfo]])) | 
          any(! winfo %in% names(iinfo))) {
        iinfo <- tree$trafo(subset = which(subset_term == i),
                            estfun = ("estfun" %in% winfo),
                            object = ("object" %in% winfo))
      }
    }
    
    tree_ret[[c(1, idx[[ichar]])]]$info <- iinfo
  }
  
  class(tree_ret) <- class(rval)
  
  return(tree_ret)
}

## variable selection: given model scores, conduct
## all M-fluctuation tests for orderins in z
.fluct_test_split <- function(estfun, z, subset, weights, obj = NULL, 
                              cluster = NULL, control)
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
  ww0[!(seq_along(ww0) %in% subset)] <- FALSE
  process <- process[ww0, , drop = FALSE]
  cluster <- cluster[ww0]
  z <- z[ww0]
  k <- NCOL(process)
  n <- NROW(process)
  
  ## stop if all collumn values in process are the same
  if(all(apply(process, 2, function(x) length(unique(x))) == 1)) 
    return(list(statistic = NA, p.value = NA))
  
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
      crossprod(as.matrix(aggregate(process, 
        by = list(cluster), FUN = sum)[, -1L, drop = FALSE]))
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
        p <- log((exp(log(0.5 - tau) + pp[1L]) + exp(log(tau - 0.49) + 
                 pchisq(x, k, lower.tail = FALSE, log.p = TRUE))) * 100)
      else
      {
        taua <- (0.51 - tau) * 50
        tau1 <- floor(taua)
        p <- log(exp(log(tau1 + 1 - taua) + pp[tau1]) + 
                 exp(log(taua-tau1) + pp[tau1 + 1L]))
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
            sqrt(pmin(x, y) * (1 - pmax(x, y)) / 
            ((pmax(x, y) * (1 - pmin(x, y))))))
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

.objfun_test_split <- function(trafo, info = info, x, bdr = NULL, j, ctrl, 
                               subset, weights, cluster, splitonly = TRUE, 
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
    if (ORDERED) {
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
    ll <- ctrl$applyfun(1:length(ux), function(u) {
      sleft <- subset[LEFT <- x[subset] <= ux[u]]
      sright <- subset[!LEFT]
      if (length(weights) > 0) {
        if (sum(weights[sleft]) < minbucket ||
            sum(weights[sright]) < minbucket)
          return(-Inf);
      } else {
        if (length(sleft) < minbucket || 
            length(sright) < minbucket)
          return(-Inf);
      }
      if (ctrl$restart) {
        linfo <- NULL
        rinfo <- NULL
      }
      linfo <- trafo(sleft, info = linfo, estfun = FALSE)
      rinfo <- trafo(sright, info = rinfo, estfun = FALSE)
      ll <- linfo$objfun + rinfo$objfun
      return(ll)
    })
    maxlogLik <- max(unlist(ll))
    if(maxlogLik > nosplitll)
      sp <- which.max(unlist(ll))
    
  } else {
    xsubs <- factor(x[subset])
    ## stop if only one level left
    if(nlevels(xsubs) < 2) {
      if (splitonly) {
        return(NULL)
      } else {
        return(list(statistic = NA, p.value = NA))
      } 
    }
    splits <- .mob_grow_getlevels(xsubs)
    ll <- ctrl$applyfun(1:nrow(splits), function(u) {
      sleft <- subset[LEFT <- xsubs %in% levels(xsubs)[splits[u,]]]
      sright <- subset[!LEFT]
      if (length(weights) > 0) {
        if (sum(weights[sleft]) < minbucket ||
            sum(weights[sright]) < minbucket)
          return(-Inf);
      } else {
        if (length(sleft) < minbucket || 
            length(sright) < minbucket)
          return(-Inf);
      }
      if (ctrl$restart) {
        linfo <- NULL
        rinfo <- NULL
      }
      linfo <- trafo(sleft, info = linfo, estfun = FALSE)
      rinfo <- trafo(sright, info = rinfo, estfun = FALSE)
      ll <- linfo$objfun + rinfo$objfun
      return(ll)
    })
    maxlogLik <- max(unlist(ll))
    if(maxlogLik > nosplitll) {
      sp <- splits[which.max(unlist(ll)),] + 1L
      levs <- levels(x)
      if(length(sp) != length(levs)) {
        sp <- sp[levs]
        names(sp) <- levs
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
    ### interpolate split-points, see https://arxiv.org/abs/1611.04561
    if (!is.factor(x) & ctrl$intersplit & sp < length(ux)) {
      sp <- (ux[sp] + ux[sp + 1]) / 2 
    } else {
      sp <- ux[sp]  ### x <= sp vs. x > sp
    }
    if (is.factor(sp)) sp <- as.integer(sp)
    ret <- partysplit(as.integer(j), breaks = sp,
                      index = 1L:2L)
  } else {
    ret <- partysplit(as.integer(j),
                      index = as.integer(sp))
  }
  return(ret)
}

