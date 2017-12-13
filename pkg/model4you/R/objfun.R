
#' Residual sum of squares
#'
#' Returns the sum of the squared residuals for a given object.
#'
#' @param object model object.
#' @param ... passed on to specific methods.
#'
#' @return sum of the squared residuals.
#' 
#' @export
#' @examples
#' ## example from ?lm
#' ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
#' trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
#' group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
#' weight <- c(ctl, trt)
#' lm.D9 <- lm(weight ~ group)
#' rss(lm.D9)
rss <- function (object, ...) 
{
  UseMethod("rss")
}

#' @rdname rss
#' @export
rss.default <- function(object, ...) {
  sum(stats::residuals(object)^2)
}



#' Objective function
#'
#' Get the contributions of an objective function. For \code{\link[stats]{glm}}
#' these are the (weighted) log-likelihood contributions, for \code{\link[stats]{lm}} the
#' negative (weighted) squared error.
#'
#' @param x model object.
#' @param newdata optional. New data frame. Can be useful for model evaluation / benchmarking.
#' @param weights optional. Prior weights. See \code{\link[stats]{glm}} or \code{\link[stats]{lm}}.
#' @param log should the log-Likelihood contributions or the Likelhood contributions be returned?
#' @param ... further arguments passed on to \code{objfun} methods.
#' 
#' @importFrom stats dbinom dnorm dpois sigma
#'
#' @return vector of objective function contributions.
#' @examples
#' ## Example taken from ?stats::glm
#' ## Dobson (1990) Page 93: Randomized Controlled Trial :
#' counts <- c(18,17,15,20,10,20,25,13,12) 
#' outcome <- gl(3,1,9)
#' treatment <- gl(3,3)
#' print(d.AD <- data.frame(treatment, outcome, counts))
#' glm.D93 <- glm(counts ~ outcome + treatment, family = poisson())
#' logLik_contributions <- objfun(glm.D93)
#' sum(logLik_contributions)
#' logLik(glm.D93)
#' 
#' @export
objfun <- function(x, newdata = NULL, ...)
{
  UseMethod("objfun")
}

#' @rdname objfun
#' @export
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
  
  val <- - (w * res^2)

  attr(val, "nobs") <- N
  attr(val, "df") <- p + 1
  class(val) <- "objfun"
  val
}

#' @rdname objfun
#' @export
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
  # FIXME check what to do in binomial if tabular data is given
  # FIXME add options for all glm families
  if(!fam %in% c("gaussian", "poisson", "binomial")) 
    stop("Haven't implemented objfun for family", fam, 
         "yet. Let me know if you want me to do that!")
  val <- switch(fam,
                 gaussian = weights * dnorm(y, mean = yhat, sd = sigma(x), 
                                            log = log),
                 poisson = weights * dpois(y, lambda = yhat, log = log),
                 binomial = {
                   n <- 1L
                   eval(family$initialize)
                   m <- if (any(n > 1)) n else weights 
                   dbinom(m * y, size = n, prob = yhat, log = log)
                   })
  

  attr(val, "nobs") <- nobs
  attr(val, "df") <- p
  class(val) <- "logLik"
  val
}