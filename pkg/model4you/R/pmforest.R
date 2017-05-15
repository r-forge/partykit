#' Compute model-based forest from model.
#'
#' Input a parametric model and get a forest.
#'
#' @param object a model object.
#' @param data data. If NULL (default) the data from the model object are used.
#' @param zformula formula describing which variable should be used for partitioning.
#' If \code{NULL} (default) all variables in data that are not in the model are used (i.e. \code{~ .}).
#' @param ntree number of trees.
#' @param control control parameters, see \code{\link[partykit]{ctree_control}}.
#' @param ... additional parameters passed on to model fit such as weights.
#'
#' @return cforest object
#' 
#' @export
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