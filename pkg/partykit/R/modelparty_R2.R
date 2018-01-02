
.mfluc_select <- function(...)
    function(model, trafo, data, subset, weights, whichvar, ctrl) {
        args <- list(...)
        ctrl[names(args)] <- args
        .select(model, trafo, data, subset, weights, whichvar, ctrl, FUN = .mfluc_test)
    }

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
    ## correct scaling of estfun for variance estimate:
    ## - caseweights=FALSE: weights are integral part of the estfun -> squared in estimate
    ## - caseweights=TRUE: weights are just a factor in variance estimate -> require division by sqrt(weights)
    meat <- if(is.null(cluster)) {
      crossprod(if(ctrl$caseweights) process/sqrt(weights) else process)
    } else {
      crossprod(as.matrix(apply(if(ctrl$caseweights) process/sqrt(weights) else process, 2L, tapply, cluster, sum)))
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
  splittry = 1L,
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
                  maxsurrogate = maxsurrogate, 
                  numsurrogate = numsurrogate, majority = majority, 
                  caseweights = caseweights, applyfun = applyfun, cores = cores, 
                  saveinfo = TRUE, ### always
                  selectfun = .mfluc_select(),
                  splitfun = .objfun_split(),
                  svselectfun = .ctree_select(),
                  svsplitfun = .ctree_split(minbucket = 0),
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

  ### <FIXME> this is already in extree_fit </FIXME>
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
    mf[[1L]] <- quote(extree_data)

    d <- eval(mf, parent.frame())
    subset <- .start_subset(d)

    weights <- model.weights(model.frame(d))

    if (is.null(control$update)) control$update <- TRUE

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
        minsize <- as.integer(ceiling(10L * n_coef / NCOL(d$yx$y)))
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


## methods concerning call/formula/terms/etc.
## (default methods work for terms and update)

formula.modelparty <- function(x, extended = FALSE, ...)
  if(extended) x$info$Formula else x$info$formula

getCall.modelparty <- function(x, ...) x$info$call

weights.modelparty <- function(object, ...) {
  fit <- object$fitted
  ww <- if(!is.null(w <- fit[["(weights)"]])) w else rep.int(1L, NROW(fit))
  structure(ww, .Names = rownames(fit))
}

## methods concerning model/parameters/loglik/etc.

coef.modelparty <- function(object, node = NULL, drop = TRUE, ...) {
  if(is.null(node)) node <- nodeids(object, terminal = TRUE)
  cf <- do.call("rbind", nodeapply(object, ids = node, FUN = function(n) info_node(n)$coefficients))
  if(drop) drop(cf) else cf
}

refit.modelparty <- function(object, node = NULL, drop = TRUE, ...)
{
  ## by default use all ids
  if(is.null(node)) node <- nodeids(object)
  
  ## estimated coefficients
  cf <- nodeapply(object, ids = node, FUN = function(n) info_node(n)$coefficients)
  
  ## model.frame
  mf <- model.frame(object)
  weights <- weights(object)
  offset <- model.offset(mf)
  cluster <- mf[["(cluster)"]]
  
  ## fitted ids
  fitted <- object$fitted[["(fitted)"]]

  ## response variables
  Y <- switch(object$info$control$ytype,
    "vector" = Formula::model.part(object$info$Formula, mf, lhs = 1L)[[1L]],
    "matrix" = model.matrix(~ 0 + ., Formula::model.part(object$info$Formula, mf, lhs = 1L)),
    "data.frame" = Formula::model.part(object$info$Formula, mf, lhs = 1L)
  )
  hasx <- object$info$nreg >= 1L | attr(object$info$terms$response, "intercept") > 0L
  X <- if(!hasx) NULL else switch(object$info$control$xtype,
    "matrix" = model.matrix(object$info$terms$response, mf),
    "data.frame" = Formula::model.part(object$info$Formula, mf, rhs = 1L)
  )
  if(!is.null(X)) {
    attr(X, "formula") <- formula(object$info$Formula, rhs = 1L)
    attr(X, "terms") <- object$info$terms$response
    attr(X, "offset") <- object$info$call$offset
  }

  suby <- function(y, index) {
    if(object$info$control$ytype == "vector") y[index] else y[index, , drop = FALSE]
  }
  subx <- if(hasx) {
    function(x, index) {
      sx <- x[index, , drop = FALSE]
      attr(sx, "contrasts") <- attr(x, "contrasts")
      attr(sx, "xlevels")   <- attr(x, "xlevels")
      attr(sx, "formula")   <- attr(x, "formula")
      attr(sx, "terms")     <- attr(x, "terms")
      attr(sx, "offset")    <- attr(x, "offset")
      sx
    }
  } else {
    function(x, index) NULL
  }
  
  ## fitting function
  afit <- object$info$fit
  
  ## refit
  rval <- lapply(seq_along(node), function(i) {
    ix <- fitted %in% nodeids(object, from = node[i], terminal = TRUE)
    args <- list(y = suby(Y, ix), x = subx(X, ix), start = cf[[i]],
      weights = weights[ix], offset = offset[ix], cluster = cluster[ix], object = TRUE)
    args <- c(args, object$info$dots)
    do.call("afit", args)$object
  })
  names(rval) <- node

  ## drop?
  if(drop & length(rval) == 1L) rval <- rval[[1L]]
  
  ## return
  return(rval)
}

apply_to_models <- function(object, node = NULL, FUN = NULL, drop = FALSE, ...) {
  if(is.null(node)) node <- nodeids(object, terminal = FALSE)
  if(is.null(FUN)) FUN <- function(object, ...) object  
  rval <- if("object" %in% object$info$control$terminal) {
    nodeapply(object, node, function(n) FUN(info_node(n)$object))
  } else {
    lapply(refit.modelparty(object, node, drop = FALSE), FUN)
  }
  names(rval) <- node
  if(drop & length(node) == 1L) rval <- rval[[1L]]
  return(rval)
}

logLik.modelparty <- function(object, dfsplit = NULL, ...)
{
  if(is.null(dfsplit)) dfsplit <- object$info$control$dfsplit
  dfsplit <- as.integer(dfsplit)
  ids <- nodeids(object, terminal = TRUE)
  ll <- apply_to_models(object, node = ids, FUN = logLik)
  dfsplit <- dfsplit * (length(object) - length(ll))
  structure(
    sum(as.numeric(ll)),
    df = sum(sapply(ll, function(x) attr(x, "df"))) + dfsplit,
    nobs = nobs(object),
    class = "logLik"
  )
}

nobs.modelparty <- function(object, ...) {
  sum(unlist(nodeapply(object,
    nodeids(object, terminal = TRUE),
    function(n) info_node(n)$nobs
  )))
}

deviance.modelparty <- function(object, ...)
{
  ids <- nodeids(object, terminal = TRUE)
  dev <- apply_to_models(object, node = ids, FUN = deviance)
  sum(unlist(dev))
}

summary.modelparty <- function(object, node = NULL, ...)
{
  ids <- if(is.null(node)) nodeids(object, terminal = TRUE) else node
  rval <- apply_to_models(object, node = ids, FUN = summary)
  if(length(ids) == 1L) rval[[1L]] else rval
}

sctest.modelparty <- function(object, node = NULL, ...)
{
  ids <- if(is.null(node)) nodeids(object, terminal = FALSE) else node
  rval <- nodeapply(object, ids, function(n)
      info_node(n)$criterion[c("statistic", "p.value"),,drop = FALSE])
  names(rval) <- ids
  if(length(ids) == 1L) rval[[1L]] else rval  
}

print.modelparty <- function(x, node = NULL,
  FUN = NULL, digits = getOption("digits") - 4L,
  header = TRUE, footer = TRUE, title = NULL, objfun = "", ...)
{
  digits <- max(c(0, digits))
  if(objfun != "") objfun <- paste(" (", objfun, ")", sep = "")
  if(is.null(title)) title <- sprintf("Model-based recursive partitioning (%s)",
    deparse(x$info$call$fit))

  if(is.null(node)) {
    header_panel <- if(header) function(party) {      
      c(title, "", "Model formula:", deparse(party$info$formula), "", "Fitted party:", "")
    } else function(party) ""
  
    footer_panel <- if(footer) function(party) {
      n <- width(party)
      n <- format(c(length(party) - n, n))
      info <- nodeapply(x, ids = nodeids(x, terminal = TRUE),
        FUN = function(n) c(length(info_node(n)$coefficients), info_node(n)$objfun))
      k <- mean(sapply(info, "[", 1L))
      of <- format(sum(sapply(info, "[", 2L)), digits = getOption("digits"))

      c("", paste("Number of inner nodes:   ", n[1L]),
        paste("Number of terminal nodes:", n[2L]),
        paste("Number of parameters per node:", format(k, digits = getOption("digits"))),
        paste("Objective function", objfun, ": ", of, sep = ""), "")
    } else function (party) ""

    if(is.null(FUN)) {
      FUN <- function(x) c(sprintf(": n = %s", x$nobs), capture.output(print(x$coefficients)))
    }
    terminal_panel <- function(node) formatinfo_node(node,
      default = "*", prefix = NULL, FUN = FUN)

    print.party(x, terminal_panel = terminal_panel,
      header_panel = header_panel, footer_panel = footer_panel, ...)
  } else {
    node <- as.integer(node)
    info <- nodeapply(x, ids = node,
      FUN = function(n) info_node(n)[c("coefficients", "objfun", "criterion")])    
    for(i in seq_along(node)) {
      if(i == 1L) {
        cat(paste(title, "\n", collapse = ""))
      } else {
        cat("\n")
      }
      cat(sprintf("-- Node %i --\n", node[i]))
      cat("\nEstimated parameters:\n")
      print(info[[i]]$coefficients)
      cat(sprintf("\nObjective function:\n%s\n", format(info[[i]]$objfun)))
      cat("\nParameter instability tests:\n")
      print(info[[i]]$criterion[c("statistic", "p.value"),,drop = FALSE])
    }
  }
  invisible(x)
}

predict.modelparty <- function(object, newdata = NULL, type = "node", ...)
{
  ## predicted node ids
  node <- predict.party(object, newdata = newdata)
  if(identical(type, "node")) return(node)

  ## obtain fitted model objects
  ids <- sort(unique(as.integer(node)))
  mod <- apply_to_models(object, node = ids)

  ## obtain predictions
  pred <- if(is.character(type)) {
    function(object, newdata = NULL, ...) predict(object, newdata = newdata, type = type, ...)
  } else {
    type
  }
  if("newdata" %in% names(formals(pred))) {    
    ix <- lapply(seq_along(ids), function(i) which(node == ids[i]))
    preds <- lapply(seq_along(ids), function(i)
      pred(mod[[i]], newdata = newdata[ix[[i]], , drop = FALSE], ...))
    nc <- NCOL(preds[[1L]])
    rval <- if(nc > 1L) {
      matrix(0, nrow = length(node), ncol = nc, dimnames = list(names(node), colnames(preds[[1L]])))
    } else {
      rep(preds[[1L]], length.out = length(node))
    }      
    for(i in seq_along(ids)) {
      if(nc > 1L) {
        rval[ix[[i]], ] <- preds[[i]]
	rownames(rval) <- names(node)
      } else {
        rval[ix[[i]]] <- preds[[i]]
	names(rval) <- names(node)
      }
    }
  } else {
    rval <- lapply(mod, pred, ...)
    if(NCOL(rval[[1L]]) > 1L) {
      rval <- do.call("rbind", rval)
      rownames(rval) <- ids
      rval <- rval[as.character(node), , drop = FALSE]
      rownames(rval) <- names(node)
      rval <- drop(rval)
    } else {
      ## provide a c() method for factors locally
      c.factor <- function(...) {
        args <- list(...)
	lev <- levels(args[[1L]])
	args[[1L]] <- unclass(args[[1L]])
	rval <- do.call("c", args)
	factor(rval, levels = 1L:length(lev), labels = lev)
      }
      rval <- do.call("c", rval)
      names(rval) <- ids
      rval <- rval[as.character(node)]
      names(rval) <- names(node)
    }
  }
  return(rval)
}

fitted.modelparty <- function(object, ...)
{
  ## fitted nodes
  node <- predict.party(object, type = "node")

  ## obtain fitted model objects
  ids <- nodeids(object, terminal = TRUE)
  fit <- apply_to_models(object, node = ids, FUN = fitted)

  nc <- NCOL(fit[[1L]])
  rval <- if(nc > 1L) {
    matrix(0, nrow = length(ids), ncol = nc, dimnames = list(names(node), colnames(fit[[1L]])))
  } else {
    rep(fit[[1L]], length.out = length(node))
  }	 
  for(i in seq_along(ids)) {
    if(nc > 1L) {
      rval[node == ids[i], ] <- fit[[i]]
      rownames(rval) <- names(node)
    } else {
      rval[node == ids[i]] <- fit[[i]]
      names(rval) <- names(node)
    }
  }

  return(rval)
}

residuals.modelparty <- function(object, ...)
{
  ## fitted nodes
  node <- predict.party(object, type = "node")

  ## obtain fitted model objects
  ids <- nodeids(object, terminal = TRUE)
  res <- apply_to_models(object, node = ids, FUN = residuals)

  nc <- NCOL(res[[1L]])
  rval <- if(nc > 1L) {
    matrix(0, nrow = length(ids), ncol = nc, dimnames = list(names(node), colnames(res[[1L]])))
  } else {
    rep(res[[1L]], length.out = length(node))
  }	 
  for(i in seq_along(ids)) {
    if(nc > 1L) {
      rval[node == ids[i], ] <- res[[i]]
      rownames(rval) <- names(node)
    } else {
      rval[node == ids[i]] <- res[[i]]
      names(rval) <- names(node)
    }
  }

  return(rval)
}

plot.modelparty <- function(x, terminal_panel = NULL, FUN = NULL, tp_args = NULL, ...) {
  if(is.null(terminal_panel)) {
    if(is.null(FUN)) {
      FUN <- function(x) {
        cf <- x$coefficients
	cf <- matrix(cf, ncol = 1, dimnames = list(names(cf), ""))
        c(sprintf("n = %s", x$nobs), "Estimated parameters:",
          strwrap(capture.output(print(cf, digits = 4L))[-1L]))
      }
    }
    terminal_panel <- do.call("node_terminal", c(list(obj = x, FUN = FUN), tp_args))
    tp_args <- NULL
  }
  plot.party(x, terminal_panel = terminal_panel, tp_args = tp_args, ...)
}

### AIC-based pruning
prune.modelparty <- function(tree, type = "AIC", ...)
{
  
  ## prepare pruning function
  if(is.character(type)) {
    type <- tolower(type)
    type <- match.arg(type, c("aic", "bic", "none"))
    
    if("lmtree" %in% class(tree)) {
      type <- switch(type,
                     "aic" = {
                       function(objfun, df, nobs) (nobs[1L] * log(objfun[1L]) + 2 * df[1L]) < (nobs[1L] * log(objfun[2L]) + 2 * df[2L])
                     }, "bic" = {
                       function(objfun, df, nobs) (nobs[1L] * log(objfun[1L]) + log(nobs[2L]) * df[1L]) < (nobs[1L] * log(objfun[2L]) + log(nobs[2L]) * df[2L])
                     }, "none" = {
                       NULL
                     })
    } else {
      type <- switch(type,
                     "aic" = {
                       function(objfun, df, nobs) (2 * - objfun[1L] + 2 * df[1L]) < (2 * - objfun[2L] + 2 * df[2L])
                     }, "bic" = {
                       function(objfun, df, nobs) (2 * - objfun[1L] + log(n) * df[1L]) < (2 * - objfun[2L] + log(n) * df[2L])
                     }, "none" = {
                       NULL
                     })   
    }
  }
  if(!is.function(type)) {
    warning("Unknown specification of 'prune'")
    type <- NULL
  }
  
  ## degrees of freedom
  dfsplit <- tree$info$control$dfsplit
  
  ## turn node to list
  node <- tree$node
  nd <- as.list(node)
  
  ## if no pruning selected
  if(is.null(type)) return(nd)
  
  ## node information (IDs, kids, ...)
  id <- seq_along(nd)
  kids <- lapply(nd, "[[", "kids")
  tmnl <- sapply(kids, is.null)
  
  ## check nodes that only have terminal kids
  check <- sapply(id, function(i) !tmnl[i] && all(tmnl[kids[[i]]]))
  while(any(check)) {
    
    ## pruning node information
    pnode <- which(check)
    objfun <- sapply(nd, function(x) x$info$objfun)
    n <- nrow(tree$fitted)
    pok <- sapply(pnode, function(i) type(
      objfun = c(objfun[i], sum(objfun[kids[[i]]])),
      df = c(length(nd[[1]]$info$coefficients), length(kids[[i]]) * length(nd[[1]]$info$coefficients) + as.integer(dfsplit)),
      nobs = c(nd[[i]]$info$nobs, n)
    ))
    
    ## do any nodes need pruning?
    pnode <- pnode[pok]
    if(length(pnode) < 1L) break
    
    ## prune
    tree <- nodeprune.party(tree, ids = pnode)
    node <- tree$node
    nd <- as.list(node)
    
    ## node information
    kids <- lapply(nd, "[[", "kids")
    tmnl <- sapply(kids, is.null)
    id <- seq_along(nd)
    check <- sapply(id, function(i) !tmnl[i] && all(tmnl[kids[[i]]]))
  }
  
  ## return pruned tree
  return(tree)
}
