#' Compute model-based forest from model.
#'
#' Input a parametric model and get a forest.
#'
#' @param model a model object.
#' @param data data. If NULL (default) the data from the model object are used.
#' @param zformula formula describing which variable should be used for partitioning.
#' Default is to use all variables in data that are not in the model (i.e. \code{~ .}).
#' @param ntree number of trees.
#' @param control control parameters, see \code{\link[partykit]{ctree_control}}.
#' @param ... additional parameters passed on to model fit such as weights.
#'
#' @return cforest object
#' 
#' @export
#' @importFrom partykit ctree_control
pmforest <- function(model, data = NULL, zformula = ~., ntree = 500L,
                     control = ctree_control(teststat = "quad", testtype = "Univ", 
                                             mincriterion = 0, saveinfo = FALSE, 
                                             lookahead = TRUE, ...), 
                     ...) {
  cl <- match.call()
  args <- .prepare_args(model = model, data = data, zformula = zformula, 
                        control = control, ntree = ntree)
  
  ## call cforest
  args$ytrafo <- function(...) .modelfit(model = model, ...)
  ret <- do.call("cforest", args)
  ret$info$model <- model
  ret$info$zformula <- zformula
  ret$call <- cl
  
  return(ret)
}