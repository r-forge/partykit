#' Compute model-based tree from model.
#'
#' Input a parametric model and get a forest.
#'
#' @param model a model object.
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
#' @example inst/examples/ex-pmtree-methods.R
#' 
#' @export
#' @import partykit 
#' @importFrom partykit ctree_control nodeids nodeapply 
#' @importFrom stats predict
pmtree <- function(model, data = NULL, zformula = ~., 
                   control = ctree_control(), coeffun = coef,
                   ...) {
  
  ### nmax not possible because data come from model
  stopifnot(all(!is.finite(control$nmax)))

  args <- .prepare_args(model = model, data = data, zformula = zformula, 
                        control = control)
  
  ## call ctree
  args$ytrafo <- function(data, weights, control, ...) 
      .modelfit(data = data, weights = weights, control = control, 
                model = model, coeffun = coeffun, ...)
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
      idn <- idx[[ichar]]
      
      if(length(idn) > 1) idn <- c(1, idn)
      iinfo <- tree_ret[[idn]]$info
      subsi <- subset_term == i
      
      if (is.null(iinfo)) {
        di <- args$data[subsi, ]
        umod <- update(model, data = di)
        iinfo <- list(estfun = estfun(umod), coefficients = coeffun(umod),
                      objfun = ifelse(class(umod)[[1]] == "lm", 
                                      sum(objfun(umod)), 
                                      logLik(umod)), 
                      model = NULL)
        tree_ret[[idn]]$info <- iinfo
      } 
      tree_ret[[idn]]$info$nobs <- sum(subsi)
      
    }
  }
  
  ## prepare return object
  class(tree_ret) <- c("pmtree", class(ret))
  tree_ret$info$model <- model
  tree_ret$info$zformula <- if(is.null(zformula)) as.formula("~ .") else 
    as.formula(zformula)
  # tree_ret$data <- data
  tree_ret$nobs <- sum(unlist(
    nodeapply(tree_ret, ids = which_terminals, function(x) x$info$nobs)
    ))
  return(tree_ret)
}
