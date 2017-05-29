#' Compute model-based forest from model.
#'
#' Input a parametric model and get a forest.
#'
#' @param object a model object.
#' @param data data. If NULL (default) the data from the model object are used.
#' @param zformula formula describing which variable should be used for partitioning.
#' Default is to use all variables in data that are not in the model (i.e. \code{~ .}).
#' @param control control parameters, see \code{\link[partykit]{ctree_control}}.
#' @param coeffun function that takes the model object and returns the coefficients. 
#' Useful when \code{coef()} does not return all coefficients (e.g. \code{survreg}).
#' @param ... additional parameters passed on to model fit such as weights.
#'
#' @return ctree object
#' 
#' @example inst/examples/ex-pmtree.R
#' 
#' @export
#' @import partykit 
#' @importFrom partykit ctree_control nodeids nodeapply 
#' @importFrom stats predict
pmtree <- function(object, data = NULL, zformula = ~., 
                   control = ctree_control(), coeffun = coef,
                   ...) {
  
  args <- .prepare_args(object = object, data = data, zformula = zformula, 
                        control = control)
  
  ## call ctree
  args$ytrafo <- function(...) .modelfit(model = object, coeffun = coeffun, ...)
  ret <- do.call("ctree", args)
  
  ### add modelinfo to teminal nodes if not there yet, but wanted
  which_terminals <- nodeids(ret, terminal = TRUE)
  # which_all <- nodeids(ret)
  
  idx <- get_paths(nodeapply(ret)[[1]], which_terminals)
  names(idx) <- which_terminals
  tree_ret <- unclass(ret)
  subset_term <- predict(ret, type = "node")
  
  if(control$saveinfo) {
    for (i in which_terminals) {
      ichar <- as.character(i)
      iinfo <- tree_ret[[c(1, idx[[ichar]])]]$info
      subsi <- subset_term == i
      
      if (is.null(iinfo)) {
        umod <- update(object, subset = subsi)
        iinfo <- list(estfun = estfun(umod), coefficients = coeffun(umod),
                      objfun = logLik(umod), object = NULL)
        tree_ret[[c(1, idx[[ichar]])]]$info <- iinfo
      } 
      tree_ret[[c(1, idx[[ichar]])]]$info$nobs <- sum(subsi)
      
    }
  }
  
  ## prepare return object
  class(tree_ret) <- c("pmtree", class(ret))
  tree_ret$info$object <- object
  tree_ret$info$zformula <- if(is.null(zformula)) as.formula("~ .") else 
    as.formula(zformula)
  # tree_ret$data <- data
  tree_ret$nobs <- sum(unlist(
    nodeapply(tree_ret, ids = which_terminals, function(x) x$info$nobs)
    ))
  return(tree_ret)
}