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
  tree_ret
  tree_ret$info$zformula <- if(is.null(zformula)) as.formula("~ .") else 
    as.formula(zformula)
  # tree_ret$data <- data
  return(tree_ret)
}



#' Methods for pmtree
#'
#' Print and summary methods for pmtree objects.
#'
#' @param x object.
#' @param node node number, if any.
#' @param FUN formatinfo function.
#' @param digits number of digits.
#' @param footer should footer be included?
#' @param ... further arguments passed on to \code{\link[partykit]{print.party}}.
#'
#' @return print
#' @export
#' @importFrom partykit print.party info_node nodeapply width formatinfo_node
#' @importFrom utils capture.output
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



#' @rdname print.pmtree 
#' 
#' @param object object.
#' @export
summary.pmtree <- function(object, node = NULL) {
  
  ids <- if(is.null(node)) nodeids(object, terminal = TRUE) else node
  info <- nodeapply(object, ids = ids, function(x) x$info)
  
  ## coefficients
  coefs <- sapply(info, function(x) x$coefficients)
  colnames(coefs) <- paste("node", colnames(coefs))
  
  ## objective functions
  objfuns <- sapply(info, function(x) x$objfun)
  if(is.null(node)) {
    brobjf <- paste0("(", round(objfuns, 2), ")")
    sobjf <- paste(brobjf, collapse = " + ")
    objfs <- paste(sobjf, "=", round(sum(objfuns), 2))
  } else {
    objfs <- objfuns
    names(objfs) <- paste("node", ids)
  }
  
  ## call
  cl <- getCall(object$info$object)
  
  ## nobs
  nobs <- sapply(info, function(x) x$nobs)
  nullnobs <- sapply(nobs, is.null)
  if(1 %in% ids & any(nullnobs)) 
    nobs[nullnobs] <- nrow(object$data)
  names(nobs) <- paste("node", ids)
    
  
  
  ret <- list(ids = ids, call = cl, coefs = coefs, objfs = objfs, nobs = nobs)
  class(ret) <- "summary.pmtree"
  return(ret)
}


#' @rdname print.pmtree 
#' 
#' @param digits minimal number of significant digits, see \code{\link[base]{print.default}}
#' @export
print.summary.pmtree <- function(x, digits = 4) {
  cat("Stratified model for node(s)", paste(x$ids, collapse = ", "))
  cat("\n\nModel call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  
  cat("Coefficients:\n")
  print(x$coefs, digits = digits)
  
  cat("\nNumber of obervations:\n")
  print(unlist(x$nobs))
  
  cat("\nObjective function:\n")
  if(is.character(x$objfs)) cat(x$objfs) else {
    print(x$objfs)
  }
}
