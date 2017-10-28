extree <- function(formula, data, subset, na.action = na.pass, weights, offset, cluster,
  ytype = "none", xtype = "none", ...)
{
  ## FIXME: allow formula = vars
  ## plus optionally resolve weights/offset/cluster
  ## warn if na.action or subset are specified

  ## call
  cl <- match.call()
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

  ## FIXME: canonicalize y/x term labels to integer and z term labels to logical/numeric

  ## FIXME: separate object with options for: discretization, condensation, some NA handling

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

  ## FIXME: subsequently fitting, testing, splitting
  ## - fit: either pre-processed _and_ subsetted data --or-- full data object plus subset vector
  ## - test: additionally needs fit output --and-- fit function
  ## - split: additionally needs test output
  ## - tbd: control of all details

  return(list(data = mf, variables = vars)) ## add terms = mt
}
