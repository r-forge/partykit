## variable selection: given model scores, conduct
## all M-fluctuation tests for orderins in z
.mfluc_test <- function(model, trafo, data, subset, weights, j, SPLITONLY, ctrl)
{

  stopifnot(!SPLITONLY)
  stopifnot(is.null(data[["yx", type = "index"]]))
  z <- data[[j]][subset]
  estfun <- model$estfun[subset,,drop = FALSE]
  cluster <- data[["(cluster)"]][subset]
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
  
  ## return version of pvalue that .extree_node deals with
  rval <- list(statistic = log(stat), p.value = log1p(-exp(pval)))
  return(rval)
}

mob_control <- function(
  alpha = 0.05,
  mincriterion = 1 - alpha,
  parm = NULL,
  bonferroni = TRUE,
  breakties = FALSE,
  nrep = 10000L,
  ordinal = c("chisq", "max", "L2"), 
  trim = 0.1,
  splittry = 2,
  vcov = c("opg", "info", "sandwich"), 
  catsplit = "binary",
  numsplit = "left",
  majority = FALSE,
  maxsurrogate = 0,
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
  nmax = Inf,
  ytype = c("vector", "data.frame", "matrix"),
  # deprecated
  verbose, xtype, prune
) {
  
  if (!missing("verbose"))
    warning("argument verbose deprecated")
  if (!missing("xtype"))
    warning("argument xtype deprecated")
  if (!missing("prune"))
    warning("argument prune deprecated")
    
  
  ## FIXME: Involve caseweights in stopping (minbucket, minsize, etc)
  ##        For now only used to compute correct n (e.g. for printing)
  
  if("estfun" %in% inner) {
    inner <- inner[inner != "estfun"]
    warning("estfun can no longer be stored in inner nodes")
    if(length(inner) == 0) inner <- NULL
  }
  
  intersplit <- numsplit == "center"
  multiway <- catsplit == "multiway"
  
  c(extree_control(criterion = "p.value",
                  logmincriterion = log(mincriterion), minsplit = minsplit, 
                  minbucket = minbucket, minprob = minprob, nmax = nmax, 
                  stump = stump, lookahead = lookahead, mtry = mtry, 
                  maxdepth = maxdepth, multiway = multiway, splittry = splittry, 
                  MIA = FALSE, maxsurrogate = maxsurrogate, 
                  numsurrogate = numsurrogate, majority = majority, 
                  caseweights = caseweights, applyfun = applyfun, cores = cores, 
                  saveinfo = TRUE, ### always
                  testflavour = "mfluc", splitflavour = "exhaustive", 
                  bonferroni = bonferroni),
    list(breakties = breakties, 
         intersplit = intersplit, parm = parm, dfsplit = dfsplit, 
         restart = restart, model = model, vcov = match.arg(vcov), 
         ordinal = match.arg(ordinal), ytype = match.arg(ytype),
         nrep = nrep, terminal = terminal, inner = inner, trim = trim))
}


mob <- function
(
  formula, 
  data, 
  subset,
  na.action = na.omit, 
  weights = NULL, 
  offset,
  cluster, 
  fit,
  control = mob_control(...), 
  converged = NULL,
  doFit = TRUE,
  ...
) {

  ## check fitting function
  fitargs <- names(formals(fit))
  if(!all(c("y", "x", "start", "weights", "offset") %in% fitargs)) {
    stop("no suitable fitting function specified")
  }

  ## augment fitting function (if necessary)
  if(!all(c("estfun", "object") %in% fitargs)) {
    afit <- function(y,
      x = NULL, start = NULL, weights = NULL, offset = NULL, cluster = NULL, estfun = FALSE, object = FALSE)
    {
      obj <- if("cluster" %in% fitargs) {
        fit(y = y, x = x, start = start, weights = weights, offset = offset, cluster = cluster, ...)
      } else {
        fit(y = y, x = x, start = start, weights = weights, offset = offset, ...)
      }
      list(
        coefficients = coef(obj),
        objfun = -as.numeric(logLik(obj)),
        estfun = if(estfun) sandwich::estfun(obj) else NULL,
        object = if(object) obj else NULL
      )
    }
  } else {
    if("cluster" %in% fitargs) {
      afit <- function(y, x = NULL, start = NULL, weights = NULL, offset = NULL, cluster = NULL, estfun = FALSE, object = FALSE)
        fit(y = y, x = x, start = start, weights = weights, offset = offset, cluster = cluster, ..., estfun = estfun, object = object)
    } else {
      afit <- function(y, x = NULL, start = NULL, weights = NULL, offset = NULL, cluster = NULL, estfun = FALSE, object = FALSE)
        fit(y = y, x = x, start = start, weights = weights, offset = offset, ..., estfun = estfun, object = object)
    }
  }


    ## set up model.frame() call
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action", "weights", "offset", "cluster"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$yx <- "matrix"
    mf$nmax <- control$nmax
    mf$ytype <- control$ytype
    ## evaluate model.frame
    mf[[1L]] <- quote(partykit::extree_data)

    d <- eval(mf, parent.frame())
    subset <- .start_subset(d)

    weights <- model.weights(d)

    update <- function(subset, weights, control)
        extree_fit(data = d, trafo = afit, converged = converged, partyvars = d$variables$z, 
                   subset = subset, weights = weights, ctrl = control)
    if (!doFit) return(list(d = d, update = update))

    ### if minsize is NULL, set to 10 * number of parameters
    if (is.null(control$minbucket) | is.null(control$minsplit)) {
        ctrl <- control
        N <- sum(complete.cases(model.frame(d, yxonly = TRUE)))
        ctrl$minbucket <- ctrl$minsplit <- N 
        ctrl$logmincriterion <- Inf
        ctrl$stump <- TRUE
    
        tree <- update(subset = subset, weights = weights, control = ctrl)
        cf <- tree$trafo(subset = subset, weights = weights, info = NULL)$coefficient
        if (is.null(cf)) {
            n_coef <- 1
        } else {
            n_coef <- length(cf)
        }
        minsize <- as.integer(ceiling(10L * n_coef))
        if (is.null(control$minbucket)) control$minbucket <- minsize
        if (is.null(control$minsplit)) control$minsplit <- minsize
    }

    tree <- update(subset = subset, weights = weights, control = control)
    trafo <- tree$trafo

    ### prepare as modelparty
    mf <- model.frame(d)
    if (length(weights) == 0) weights <- rep(1, nrow(mf))

    fitted <- data.frame("(fitted)" = fitted_node(tree$nodes, mf),
                         "(weights)" = weights,
                         check.names = FALSE)

    fitted[[3]] <- y <- mf[, d$variables$y, drop = TRUE]
    names(fitted)[3] <- "(response)"

    control$ytype <- ifelse(is.vector(y), "vector", class(y))
    # x <- model.matrix(modelf, data = mmf)
    control$xtype <- "matrix" # TODO: find out when to use data.frame

    ## return party object
    rval <- party(tree$nodes, 
                  data = if(control$model) mf else mf[0,],
                  fitted = fitted,
                  terms = d$terms$all,
                  info = list(
                    call = match.call(),
                    formula = formula,
                    Formula = as.Formula(formula),
                    terms = list(response = d$terms$yx, partitioning = d$terms$z),
                    fit = fit,
                    control = control,
                    dots = list(...),
                    nreg = NCOL(d$yx$x)
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
        iinfo <- trafo(subset = which(subset_term == i), weights = weights, info = NULL,
                       estfun = ("estfun" %in% winfo),
                       object = ("object" %in% winfo))
      }
    }

    tree_ret[[c(1, idx[[ichar]])]]$info <- iinfo
  }

  class(tree_ret) <- class(rval)

  return(tree_ret)
}

model.frame.modelparty <- function(formula, ...)
{
  mf <- formula$data
  if(nrow(mf) > 0L) return(mf)

  dots <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0L)]
  mf <- formula$info$call
  mf <- mf[c(1L, match(c("formula", "data", "subset", "na.action"), names(mf), 0L))]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf[names(nargs)] <- nargs
  if(is.null(env <- environment(formula$info$terms))) env <- parent.frame()
  mf$formula <- Formula::Formula(as.formula(mf$formula))
  eval(mf, env)
}
