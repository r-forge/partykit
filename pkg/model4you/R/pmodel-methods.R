
#' Objective function of personalised models
#'
#' Get the contributions of an objective function (e.g. likelihood contributions)
#' and the sum thereof (e.g. log-Likelihood).
#'
#' @param object object of class pmodel.identity (obtained by \code{pmodel(..., fun = identity)}). 
#' @param ... additional parameters passed on to \code{\link{objfun}}.
#' @param add_df it is not very clear what the degrees of freedom are in personalised models.
#'  With this argument you can add/substract degrees of freedom at your convenience. Default
#'  is \code{0} which means adding up the degrees of freedom of all individual models.
#' 
#' @details 
#' Note that \code{logLik.pmodel.identity} returns the sum of contributions and
#' thus not neccessarily the log-Likelihood, if the objective function is not
#' the log-Likelihood. For example if the base model is a linear model (\code{lm})
#' the objective function are the negative squared residuals and thus 
#' logLik.pmodel.identity will return the sum of squared residuals.
#'
#' For examples see \code{\link{pmodel}}.
#' 
objfun.pmodel.identity <- function(object, ...) {
  
  nobs <- length(object)
  data <- attr(object, "data")
  
  ofl <- lapply(seq_len(nobs), function(i) objfun(object[[i]], newdata = data[i, ], ...))
  of <- unlist(ofl)
  
  attr(of, "class") <- "objfun"
  attr(of, "nobs") <- nobs
  attr(of, "df") <- sapply(ofl, attr, "df")
  
  return(of)
}

#' @rdname objfun.pmodel.identity
#' 
#' @export
logLik.pmodel.identity <- function(object, add_df = 0, ...) {
  
  of <- objfun.pmodel.identity(object, ...)
  
  structure(
    sum(of),
    df = sum(attr(of, "df")) + add_df,
    nobs = object$nobs,
    class = "logLik"
  )
}