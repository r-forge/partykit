objfun <- function(x, newdata = NULL, ...)
{
  UseMethod("objfun")
}


objfun.lm <- function(x, newdata = NULL, weights = NULL, ...)
{
  if (!missing(...)) 
    warning("extra arguments discarded")
  if (inherits(x, "mlm")) 
    stop("'objfun.lm' does not support multiple responses")
  
  ## get residuals
  if(is.null(newdata)) {
    res <- x$residuals
  } else {
    yhat <- predict(x, newdata = newdata)
    modformula <- Formula::as.Formula(x$call$formula)
    yformula <- formula(modformula, lhs = 1, rhs = 0)
    y <- get_all_vars(yformula, data = newdata)
    res <- y - yhat
  }
  
  p <- x$rank
  N <- length(res)
  
  ## get weights
  if (is.null(w <- weights)) {
    w <- rep.int(1, N)
    if(is.null(newdata) & !is.null(x$weights)) w <- x$weights 
  } else {
    excl <- w == 0
    if (any(excl)) {
      res <- res[!excl]
      N <- length(res)
      w <- w[!excl]
    }
  }
  
  val <- w * res^2

  attr(val, "nobs") <- N
  attr(val, "df") <- p + 1
  class(val) <- "objfun"
  val
}


objfun.glm <- function(x, newdata = NULL, weights = NULL, log = TRUE, ...) 
  {
  if (!missing(...)) 
    warning("extra arguments discarded")
  
  ## yhat / mu
  yhat <- predict(x, newdata = newdata, type = "response")
  nobs <- sum(!is.na(yhat))
  
  ## weights
  if (is.null(weights)) {
    weights <- rep.int(1, nobs)
    if(is.null(newdata) & !is.null(x$weights)) weights <- x$prior.weights 
  }
  
  ## y
  if(is.null(newdata)) newdata <- x$data
  modformula <- Formula::as.Formula(x$call$formula)
  yformula <- formula(modformula, lhs = 1, rhs = 0)
  y <- get_all_vars(yformula, data = newdata)[,,drop = TRUE]

  
  ## info from family
  family <- family(x)
  fam <- family$family
  p <- x$rank
  if (fam %in% c("gaussian", "Gamma", "inverse.gaussian")) 
    p <- p + 1
  
  ## likelihood
  # dev <- family$dev.resids(y, yhat, w)
  # val <- sapply(1:n, function(i) 1 - family$aic(y[i], 1, yhat[i], w[i], dev[i])/2)
  
  # FIXME check what to do in binomial if tabular data is given
  # FIXME add options for all glm families
  if(!fam %in% c("gaussian", "poisson", "binomial")) stop("Haven't implemented objfun for family", fam, "yet. Let me know if you want me to do that!")
  val <- switch(fam,
                 gaussian = weights * dnorm(y, mean = yhat, sd = sigma(x), log = log),
                 poisson = weights * dpois(y, lambda = yhat, log = log),
                 binomial = {
                   eval(family$initialize)
                   m <- if (any(n > 1)) n else weights 
                   dbinom(m * y, size = n, prob = yhat, log = log)
                   })
  

  attr(val, "nobs") <- nobs
  attr(val, "df") <- p
  class(val) <- "logLik"
  val
}