#' Compute model-based forest from model.
#'
#' Input a parametric model and get a forest.
#'
#' @param object a model object.
#' @param data data. If NULL (default) the data from the model object are used.
#' @param zformula formula describing which variable should be used for partitioning.
#' If \code{NULL} (default) all variables in data that are not in the model are used (i.e. \code{~ .}).
#' @param control control parameters, see \code{\link[partykit]{ctree_control}}.
#' @param ... additional parameters passed on to model fit such as weights.
#'
#' @return cforest object
#' 
#' @export
pmtree <- function(object, data = NULL, zformula = NULL, 
                     control = ctree_control(), 
                     ...) {
  
  args <- .prepare_args(object = object, data = data, zformula = zformula, 
                        control = control)
  
  ## call cforest
  args$ytrafo <- function(...) .modelfit(model = object, ...)
  ret <- do.call("ctree", args)
  ret$info$object <- object
  
  return(ret)
}