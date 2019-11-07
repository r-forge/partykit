# -------------------------------------------------------------------
# Trafo returning the response (identity)
# -------------------------------------------------------------------
## Trafo with estfun = y
trafo_identity <- function(subset, data, weights = NULL, info = NULL, estfun = TRUE, object = TRUE) {
  
  if(! is.null(weights))  stop("weights must be null")

  ## Get subset 
  y <- data[[1, "origin"]][subset]  # FIXME: (ML, LS) data copy? no aggregation possible!


  ## Get weights for subset
  weights <- if(is.null(weights) || (length(weights)==0L)) rep.int(1, NROW(y))[subset] else weights[subset]

  ## Build estfun and set values not in subset to zero
  ef <- as.matrix(y[subset])
  ef[-subset, ] <- 0  # FIXME: (ML) zero or NA?

  ## Return list
  rval <- list(
    estfun = if(estfun) ef else NULL,
    unweighted = FALSE,  # FIXME: estfun is weighted, extree_fit reverts weighting
    coefficients = 1,  # FIXME: (ML) what is coef here ?
    objfun = 0,  # FIXME: (ML) what is the objfun here?
    object = NULL,
    nobs = NROW(y),  # FIXME: (ML, LS) needed?
    converged = TRUE  # FIXME: (ML, LS) always converged?
  )

  return(rval)

}


# -------------------------------------------------------------------
# Trafo for numerical response w/o regressor matrix
# -------------------------------------------------------------------
trafo_num <- function(subset, data, weights = NULL, offset = NULL, info = NULL, 
                      estfun = TRUE, object = FALSE) {

  ## Get data and apply offset
  y <- data[[1, "origin"]]  # FIXME: (ML, LS) data copy? no aggregation possible!
  if (!is.null(offset)) y <- y - offset

  ## Get weights for subset
  weights <- if(is.null(weights) || (length(weights)==0L)) rep.int(1, NROW(y))[subset] else weights[subset]

  ## Calculate res and rss for subset
  m <- mean(y[subset] * weights) / mean(weights)
  res <- y[subset] - m
  rss <- sum(res^2 * weights)

  ## Build estfun with original dimension and fill up subsetted indices
  ef <- matrix(0, ncol = NCOL(y), nrow = NROW(y))
  ef[subset, ] <- res * weights

  ## Return list
  rval <- list(
    estfun = if(estfun) ef else NULL,
    unweighted = FALSE,  # FIXME: estfun is weighted, extree_fit reverts weighting
    coefficients = c("mean" = m),
    objfun = -rss,
    object = if(object) list(nuisance = c("log(variance)" = log(rss/sum(weights)))) else NULL,
    nobs = NROW(d[[1, "origin"]]),  # FIXME: (ML, LS) needed?
    converged = TRUE  # FIXME: (ML, LS) always converged?
  )
  
  return(rval)

}


# -------------------------------------------------------------------
# Trafo for categorical response w/o regressor matrix
# -------------------------------------------------------------------
trafo_cat <- function(subset, data, weights = NULL, offset = NULL, info = NULL, 
                      estfun = TRUE, object = FALSE) {

  ## Get data and apply offset
  ys <- data[[1, "origin"]][subset]  # FIXME: (ML, LS) data copy? no aggregation possible!

  ## Get weights for subset
  weights <- if(is.null(weights) || (length(weights)==0L)) rep.int(1, NROW(ys)) else weights[subset]

  ## tables and probabilities
  tab <- tapply(weights, ys, sum)
  tab[is.na(tab)] <- 0L
  pr <- tab/sum(tab)
  alias <- tab == 0L
  ix1 <- which(!alias)[1L]
  if(estfun) ef <- matrix(0, nrow = length(ys), ncol = length(tab),
    dimnames = list(names(ys), names(tab)))
  
  ## Setup return list if alias < 2
  if(sum(!alias) < 2L) {
    return(list(
      estfun = NULL,
      unweighted = FALSE,  # FIXME: estfun is weighted, extree_fit reverts weighting
      coefficients = log(pr[-ix1]) - log(pr[ix1]),
      objfun = 0,
      object = NULL,
      nobs = NROW(d[[1, "origin"]]),  # FIXME: (ML, LS) needed?
      converged = TRUE  # FIXME: (ML, LS) always converged?
    ))
  }
  
  ## Setup return list if alias > 2
  rval <- list(
    estfun = NULL,
    unweighted = FALSE,  # FIXME: estfun is weighted, extree_fit reverts weighting
    coefficients = log(pr[-ix1]) - log(pr[ix1]),
    objfun = -sum(tab[tab > 0L] * log(pr[tab > 0L])),
    object = NULL,
    nobs = NROW(d[[1, "origin"]]),  # FIXME: (ML, LS) needed?
    converged = TRUE  # FIXME: (ML, LS) always converged?
  )
  
  ## Build estfun with original dimension and fill up subsetted indices
  if(estfun) {
    rval$estfun <- matrix(0, ncol = length(tab), nrow = NROW(data[[1, "origin"]]),
      dimnames = list(names(ys), names(tab)))[, -1L, drop = FALSE]
    
    cf <- log(pr) - log(pr[ix1])
    ef[] <- rep(-pr, each = nrow(ef))
    ef[cbind(1:nrow(ef), as.numeric(ys))] <- (1 - pr[ys])
    ef <- ef[, !alias, drop = FALSE]
    ef <- ef[, -1L, drop = FALSE]
    rval$estfun[subset, ] <- ef * weights
  }
  
  return(rval)

}


# -------------------------------------------------------------------
# Trafo for numerical response w/ regressor matrix (adapted of lmfit)
# -------------------------------------------------------------------
trafo_lm <- function(subset, data, weights = NULL, offset = NULL, info = NULL, 
                     estfun = TRUE, object = FALSE, ...) {

  ## Get data and apply offset
  ys <- data[[1, "origin"]][subset]  # FIXME: (ML, LS) data copy? no aggregation possible!
  xs <- data$yx$x[subset, ]  # FIXME: (ML) needs to be done nicer! data copy? no aggregation possible!

  ## Get weights for subset
  weights <- if(is.null(weights) || (length(weights)==0L)) rep.int(1, NROW(ys)) else weights[subset]

  ## add intercept-only regressor matrix (if missing)
  ## NOTE: does not have terms/formula
  if(is.null(xs)) xs <- matrix(1, nrow = NROW(ys), ncol = 1L,
    dimnames = list(NULL, "(Intercept)"))

  ## call lm fitting function
  if(is.null(weights) || identical(as.numeric(weights), rep.int(1, length(weights)))) {
    z <- lm.fit(xs, ys, offset = offset, ...)
    weights <- 1
  } else {
    z <- lm.wfit(xs, ys, w = weights, offset = offset, ...)
  }

  ## list structure
  rval <- list(
    estfun = NULL,
    unweighted = FALSE,  # FIXME: estfun is weighted, extree_fit reverts weighting
    coefficients = z$coefficients,
    objfun = -sum(weights * z$residuals^2),  # FIXME: (ML) changed to negative sum
    object = NULL,
    nobs = NROW(d[[1, "origin"]]),  # FIXME: (ML, LS) needed?
    converged = TRUE  # FIXME: (ML, LS) always converged?
  )

  ## add estimating functions (if desired)
  if(estfun) {
    rval$estfun <- as.vector(z$residuals) * weights * xs[, !is.na(z$coefficients), drop = FALSE]
  }

  ## add model (if desired)
  if(object) {
    class(z) <- c(if(is.matrix(z$fitted)) "mlm", "lm")
    z$offset <- if(is.null(offset)) 0 else offset
    z$contrasts <- attr(xs, "contrasts")
    z$xlevels <- attr(xs, "xlevels")

    cl <- as.call(expression(lm))
    cl$formula <- attr(xs, "formula")
    if(!is.null(offset)) cl$offset <- attr(xs, "offset")
    z$call <- cl
    z$terms <- attr(xs, "terms")

    rval$object <- z
  }

  return(rval)

}
