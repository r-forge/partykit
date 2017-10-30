## extensible tree (model) function
extree <- function(formula, data, subset, na.action = na.pass, weights, offset, cluster,
  ytype = "none", xtype = "none", ...)
{
  ## call
  cl <- match.call()

  noformula <- !inherits(formula, "formula")
  if(noformula) {

    ## formula needs to be a 'list' (if it is not a 'formula')
    if(!inherits(formula, "list")) stop("unsupported specification of 'formula'")    

    ## specified formula elements and overall call elements
    fonam <- names(formula)
    clnam <- names(cl)[-1L]
    vanam <- c("y", "x", "z", "weights", "offset", "cluster")
    
    ## y and z (and optionally x) need to be in formula
    if(!all(c("y", "z") %in% fonam)) stop("'formula' needs to specify at least a response 'y' and partitioning variables 'z'")
    if(!("x" %in% fonam)) formula$x <- NULL
    
    ## furthermore weights/offset/cluster may be in formula or call
    vars <- formula[vanam]
    names(vars) <- vanam
    if("weights" %in% clnam) {
      clvar <- try(weights, silent = TRUE)
      vars[["weights"]] <- c(vars[["weights"]], if(!inherits(clvar, "try-error")) clvar else deparse(cl$weights))
    }
    if("offset" %in% clnam) {
      clvar <- try(offset, silent = TRUE)
      vars[["offset"]] <- c(vars[["offset"]], if(!inherits(clvar, "try-error")) clvar else deparse(cl$offset))
    }
    if("cluster" %in% clnam) {
      clvar <- try(cluster, silent = TRUE)
      vars[["cluster"]] <- c(vars[["cluster"]], if(!inherits(clvar, "try-error")) clvar else deparse(cl$cluster))
    }
    
    ## sanity checking
    for(v in vanam) {
      if(!is.null(vars[[v]]) && !(is.numeric(vars[[v]]) | is.character(vars[[v]]) | is.logical(vars[[v]]))) {
        warning(sprintf("unknown specification of '%s', must be character, numeric, or logical", v))
        vars[v] <- list(NULL)
      }
    }
    if(!missing(subset)) warning("'subset' argument ignored in list specification of 'formula'")
    if(!missing(na.action)) warning("'na.action' argument ignored in list specification of 'formula'")    

  } else {

    mf <- match.call(expand.dots = FALSE)
    if(missing(data)) data <- environment(formula)
    m <- match(c("formula", "data", "subset", "na.action", "weights", "offset", "cluster"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE

    ## formula FIXME: y ~ . or y ~ x | .
    ## default: dot = "sequential"
    ## idea: check (x and) z vs. deparse(cl$weights), deparse(cl$offset), deparse(cl$cluster)
    oformula <- as.formula(formula)
    formula <- Formula::as.Formula(formula)
    if(length(formula)[2L] < 2L) {
      formula <- Formula::Formula(formula(Formula::as.Formula(formula(formula), ~ 0), rhs = 2L:1L))
      xreg <- FALSE
    } else {
      if(length(formula)[2L] > 2L) {
        formula <- Formula::Formula(formula(formula, rhs = 1L:2L))
        warning("Formula must not have more than two RHS parts")
      }
      xreg <- TRUE
    }
    mf$formula <- formula
    ## FIXME: allow m-part formulas with m > 3, e.g., for beta regression or heteroscedastic tobit etc.
    ## the m-th part is "z" and then there is a list of "x" matrices

    ## evaluate model.frame
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())

    ## extract terms in various combinations
    mt <- list(
      "all" = terms(formula, data = data),
      "y"   = terms(formula, data = data, rhs = 0L),
      "x"   = if(xreg) terms(formula, data = data, lhs = 0L, rhs = 1L) else NULL,
      "yx"  = if(xreg) terms(formula, data = data, rhs = 1L) else NULL,
      "z"   = terms(formula, data = data, lhs = 0L, rhs = 2L)
    )

    ## extract variable lists
    vars <- list(
      y = attr(mt$y, "term.labels"),
      x = attr(mt$x, "term.labels"),
      z = attr(mt$z, "term.labels"),
      weights = if("(weights)" %in% names(mf)) "(weights)" else NULL,
      offset  = if("(offset)"  %in% names(mf)) "(offset)"  else NULL,
      cluster = if("(cluster)" %in% names(mf)) "(cluster)" else NULL
    )
    ## FIXME check handling of univariate vs. multivariate responses
    ymult <- length(vars$y) >= 1L
    if(!ymult) vars$y <- names(mf)[1L]
    ## FIXME: store information which variable(s) went into (weights), (offset), (cluster)

    ## check wether offset was inside the formula
    if(!is.null(off <- attr(mt$x, "offset"))) {
      if(is.null(vars$offset)) mf[["(offset)"]] <- rep.int(0, nrow(mf))
      for(i in off) mf[["(offset)"]] <- mf[["(offset)"]] + mf[[i]]
      vars$offset <- "(offset)"
    }
  }
  
  ## canonicalize y/x/z term labels
  vanam <- if(noformula) names(data) else names(mf)
  ## z to numeric
  if(is.null(vars$z)) stop("at least one 'z' variable must be specified")
  if(is.integer(vars$z)) vars$z <- vanam[vars$z]
  if(is.character(vars$z)) vars$z <- vanam %in% vars$z
  if(is.logical(vars$z)) vars$z <- as.numeric(vars$z)
  if(any(is.na(vars$z))) vars$z[is.na(vars$z)] <- 0
  vars$z <- as.numeric(vars$z)
  ## all others to integer
  for(v in c("y", "x", "weights", "offset", "cluster")) {
    if(!is.null(vars[[v]])) {
      if(is.character(vars[[v]])) vars[[v]] <- match(vars[[v]], vanam)
      if(is.logical(vars[[v]])) vars[[v]] <- which(vars[[v]])
      if(any(is.na(vars[[v]]))) {
        vars[[v]] <- vars[[v]][!is.na(vars[[v]])]
        warning(sprintf("only found the '%s' variables: %s", v, paste(vanam[vars[[v]]], collapse = ", ")))
      }
    }
    vars[[v]] <- as.integer(vars[[v]])
  }
  if(is.null(vars$y)) stop("at least one 'y' variable must be specified")

  

  ## FIXME: separate object with options for: discretization, condensation, some NA handling

  if((ytype != "none") | (xtype != "none")) {

    ## pre-process y and x if desired
    if(ytype == "matrix") {
      mf[["(y)"]] <- model.matrix(~ 0 + ., Formula::model.part(formula, mf, lhs = TRUE))
      vars$y <- "(y)"
    }
    if(xreg && xtype == "matrix") {
      mf[["(x)"]] <- model.matrix(if(ymult) mt$x else mt$yx, mf)
      if(ncol(mf[["(x)"]]) < 1L) {
        mf[["(x)"]] <- NULL
        xreg <- FALSE
        vars$x <- NULL
      } else {
        attr(mf[["(x)"]], "formula") <- formula(formula, rhs = 1L)
        attr(mf[["(x)"]], "terms") <- if(ymult) mt$x else mt$yx
        attr(mf[["(x)"]], "offset") <- cl$offset
        vars$x <- "(x)"
      }
    }

  }

  ## FIXME: subsequently fitting, testing, splitting
  ## - fit: either pre-processed _and_ subsetted data --or-- full data object plus subset vector
  ## - test: additionally needs fit output --and-- fit function
  ## - split: additionally needs test output
  ## - tbd: control of all details

  return(list(data = if(noformula) data else mf, variables = vars)) ## add terms = mt
}
