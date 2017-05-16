#' Compute model-based forest from model.
#'
#' Input a parametric model and get a forest.
#'
#' @param object a model object.
#' @param data data. If NULL (default) the data from the model object are used.
#' @param zformula formula describing which variable should be used for partitioning.
#' Default is to use all variables in data that are not in the model (i.e. \code{~ .}).
#' @param control control parameters, see \code{\link[partykit]{ctree_control}}.
#' @param ... additional parameters passed on to model fit such as weights.
#'
#' @return ctree object
#' 
#' @export
pmtree <- function(object, data = NULL, zformula = ~., 
                   control = ctree_control(),
                   ...) {
  
  args <- .prepare_args(object = object, data = data, zformula = zformula, 
                        control = control)
  
  ## call ctree
  args$ytrafo <- function(...) .modelfit(model = object, ...)
  ret <- do.call("ctree", args)
  
  ### add modelinfo to teminal nodes if not there yet, but wanted
  which_terminals <- nodeids(ret, terminal = TRUE)
  which_all <- nodeids(ret)
  
  idx <- lapply(which_all, partykit:::.get_path, obj = nodeapply(ret)[[1]])
  names(idx) <- which_all
  tree_ret <- unclass(ret)
  subset_term <- predict(ret, type = "node")
  
  if(control$saveinfo) {
    for (i in which_terminals) {
      ichar <- as.character(i)
      iinfo <- tree_ret[[c(1, idx[[ichar]])]]$info
      subsi <- subset_term == i
      
      if (is.null(iinfo)) {
        umod <- update(object, subset = subsi)
        iinfo <- list(estfun = estfun(umod), coefficients = coef(umod),
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
  return(tree_ret)
}



#' Print method for pmtree
#'
#' Print pmtree
#'
#' @param x pmtree object.
#' @param node node number, if any.
#' @param FUN formatinfo function.
#' @param digits number of digits.
#' @param footer should footer be included?
#' @param ... further arguments passed on to prin.party.
#'
#' @return print
#' @export
print.pmtree <- function(x, node = NULL,
                         FUN = NULL, digits = getOption("digits") - 4L,
                         footer = TRUE, ...)
{
  digits <- max(c(0, digits))
  title <- paste("Partitioned model:\n", paste(deparse(getCall(x$info$object)), 
                                               sep = "\n", collapse = "\n"),
                 "\nPartitioning variables:", deparse(x$info$zformula))
  
  if(is.null(node)) {
    header_panel <- function(party) ""
    
    footer_panel <- if(footer) function(party) {
      n <- width(party)
      n <- format(c(length(party) - n, n))
      info <- nodeapply(x, ids = nodeids(x, terminal = TRUE),
                        FUN = function(n) c(length(info_node(n)$coefficients), info_node(n)$objfun))
      k <- mean(sapply(info, "[", 1L))
      of <- format(sum(sapply(info, "[", 2L)), digits = getOption("digits"))
      
      c("", paste("Number of inner nodes:   ", n[1L]),
        paste("Number of terminal nodes:", n[2L]),
        paste("Number of parameters per node:", format(k, digits = getOption("digits"))),
        paste("Log-likelihood: ", of, sep = ""), "")
    } else function (party) ""
    
    if(is.null(FUN)) {
      FUN <- function(x) c(sprintf(": n = %s", x$nobs), capture.output(print(x$coefficients)))
    }
    terminal_panel <- function(node) formatinfo_node(node,
                                                     default = "*", prefix = NULL, FUN = FUN)
    
    print.party(x, terminal_panel = terminal_panel,
                header_panel = header_panel, footer_panel = footer_panel, ...)
  } else {
    node <- as.integer(node)
    info <- nodeapply(x, ids = node,
                      FUN = function(n) info_node(n)[c("coefficients", "objfun", "criterion")])    
    for(i in seq_along(node)) {
      if(i == 1L) {
        cat(paste(title, "\n", collapse = ""))
      } else {
        cat("\n")
      }
      cat(sprintf("-- Node %i --\n", node[i]))
      cat("\nEstimated parameters:\n")
      print(info[[i]]$coefficients)
      cat(sprintf("\nObjective function:\n%s\n", format(info[[i]]$objfun)))
      cat("\nParameter instability tests:\n")
      print(info[[i]]$criterion)
    }
  }
  invisible(x)
}
