
#' Personalised model
#'
#' Compute personalised models from cforest object.
#'
#' @param x cforest object or matrix of weights.
#' @param object model object. If NULL the model in \code{x$info$object} is used.
#' @param newdata new data. If NULL cforest learning data is used. Ignored if \code{x} is a matrix.
#' @param OOB In case of using the learning data, should patient similarities be
#' computed out of bag?
#' @param fun function to apply on the personalised model before returning. The
#' default \code{coef} returns a matrix of personalised coefficients. For returning
#' the model objects use \code{identity}.
#'
#' @return depends on fun.
#' 
#' @example inst/examples/ex-pmodel.R
#' 
#' @export
pmodel <- function(x = NULL, object = NULL, newdata = NULL, OOB = TRUE, fun = coef) {
  
  ## compute similarity weights
  if(is.matrix(x)) {
    if(is.null(object)) stop("When x is a matrix, object must not be NULL. Please enter a model object.")
    pweights <- x
  } else {
    if(is.null(object)) object <- x$info$object
    pweights <- predict(x, type = "weights", newdata = newdata, OOB = OOB)
  }
  
  ## personalised model or model coefficients
  get_pmod <- function(w) {
    
    if(sum(w) == 0) stop("The weights for one observation are all 0. A solution may be increasing ntree.")
    
    ## compute the model
    pmod <- update(object, weights = w, subset = w > 0)
    
    ## return model or coefficients
    fun(pmod)
  }
  
  ret <- apply(pweights, 2, get_pmod)
  if(class(ret) == "matrix") ret <- t(ret)
  attr(ret, "modelcall") <- getCall(object)
  
  return(ret)
}