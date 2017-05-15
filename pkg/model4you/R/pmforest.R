
#' Fit function when model object is given
#'
#' Use update function to refit model and extract info such as coef, logLik and 
#' estfun.
#'
#' @param formula ignored but required by \code{.ctreetrafo}.
#' @param model model object.
#' @param data data.
#' @param weights weights.
#' @param cluster cluster.
#' @param ctrl control options from \code{ctree_control}.
#' @param parm which parameters should be used for instability test?
#'
#' @return A function returning a list of
#' \item{coefficients}{ \code{coef}. }
#' \item{objfun}{ \code{logLik}. }
#' \item{object}{ the model object. }
#' \item{converged}{ Did the model converge? }
#' \item{estfun}{ \code{estfun}. }
.modelfit <- function(formula, # ignored but needed for .ctreetrafo
                      model, data, weights, cluster, ctrl, parm = NULL) {
  
  fitfun <- function(subset, estfun = TRUE, object = FALSE, info = NULL) { 
    
    ## compute model on data
    if(length(weights) == 0) weights <- rep(1, NROW(data))
    mod <- tryCatch(update(object = model, data = data, subset = subset),
                    warning = function(w) list(converged = FALSE),
                    error = function(e) {
                      list(converged = FALSE)
                    })
    
    ## stop if error
    if(!(is.null(mod$converged)) && !(mod$converged)) 
      return(list(converged = FALSE))
    
    # mod <- try(update(object = model, data = dat),
    #            silent = TRUE)
    
    ## get convergence info
    if (is.null(ctrl$converged)) {
      conv <- if (is.null(mod$converged)) TRUE else mod$converged
    } else {
      conv <- ctrl$converged(mod, data, subset)
    }
    
    ## prepare return list
    ret <- list(coefficients = coef(mod), objfun = logLik(mod),
                object = if(object) mod else NULL,
                converged = conv)
    
    ## add estfun if wanted
    if(estfun) {
      ef <- estfun(mod)
      ret$estfun <- matrix(0, nrow = NROW(data), ncol = NCOL(ef))
      ret$estfun[subset,] <- ef
      if(!is.null(parm)) ret$estfun <- ret$estfun[, parm]
    }
    
    return(ret)
  }
  
  return(fitfun)
}



#' Compute model-based forest from model.
#'
#' Input a parametric model and get a forest.
#'
#' @param object a model object.
#' @param data data. If NULL (default) the data from the model object are used.
#' @param zformula formula describing which variable should be used for partitioning.
#' If \code{NULL} (default) all variables in data that are not in the model are used (i.e. \code{~ .}}.
#' @param ntree number of trees.
#' @param control control parameters, see \code{ctree_control}.
#' @param ... additional parameters passed on to model fit such as weights.
#'
#' @return cforest object
pmforest <- function(object, data = NULL, zformula = NULL, ntree = 2,
                     control = ctree_control(lookahead = TRUE, mincriterion = 0, saveinfo = FALSE, ...), 
                     ...) {
  
  args <- .prepare_args(object = object, data = data, zformula = zformula, 
                        control = control, ntree = ntree)
  
  ## call cforest
  args$ytrafo <- function(...) .modelfit(model = object, ...)
  ret <- do.call("cforest", args)
  ret$info$object <- object
  
  return(ret)
}