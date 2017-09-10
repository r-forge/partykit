
#' pmtree predictions
#'
#' Compute predictions from pmtree object. 
#'
#' @param object pmtree object.
#' @param newdata an optional data frame in which to look for variables with
#'    which to predict, if omitted, \code{object$data} is used.
#' @param type character denoting the type of predicted value. The terminal node 
#' is returned for \code{"node"}. If \code{type = "pass"} the model 
#' predict method is used and arguments can be passed to it via \code{predict_args}.
#' @param predict_args If \code{type = "pass"} arguments can be passed on to the
#' model predict function.
#' @param ... ignored.
#'
#' @return predictions
#' 
#' @example inst/examples/ex-pmtree-methods.R
#' 
#' @importFrom partykit predict.party
#' @export
predict.pmtree <- function(object, newdata = NULL, type = "node", predict_args = list(), ...) {
  
  ## node
  # terminals <- nodeids(object, terminal = TRUE)
  node <- predict.party(object, newdata = newdata, type = "node")
  if(type == "node") return(node)
  
  if(is.null(newdata)) newdata <- object$data
  
  ## response
  if(type == "pass") {
    trdatnodes <- object$fitted["(fitted)"]
    
    # predict outcome using respective models
    unode <- sort(unique(node))
    newdata$.node <- node
    newdata$.id <- seq_len(NROW(newdata))
    pr <- lapply(unode, function(nd) {
      
      # model
      mod <- update(object$info$model, 
                    subset = (trdatnodes == nd))
      
      # prediction
      args <- c(list(object = mod,
                     newdata = newdata[newdata$.node == nd, ]),
                predict_args)
      pred <- do.call(predict, args = args)
      data.frame(pred, .id = newdata$.id[newdata$.node == nd])
    })
    pr <- do.call(rbind, pr)
    
    ## return sorted predictions
    pred <- pr[order(pr$.id), ]
    pred$.id <- NULL
    return(pred)

  }
}


#' Extract sum of objective functions
#'
#' Extract sum of objective functions of all terminal nodes. By default the degrees
#' of freedom from the models are used but optionally degrees of freedom for splits
#' can be incorporated.
#'
#' @param object pmtree object.
#' @param dfsplit degrees of freedom per selected split.
#' @param ... ignored.
#'
#' @return Returns an object of class \code{\link[stats]{logLik}}.
#' 
#' @export
logLik.pmtree <- function(object, dfsplit = 0, ...) {
  ids <- nodeids(object, terminal = TRUE) 
  info <- nodeapply(object, ids = ids, function(x) x$info)
  dfs <- (length(object) - width(object)) * dfsplit
  
  ll <- lapply(info, function(x) x$objfun)
  structure(
    sum(as.numeric(ll)),
    df = sum(sapply(ll, function(x) attr(x, "df"))) + dfs,
    nobs = object$nobs,
    class = "logLik"
  )
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
  title <- paste("Partitioned model:\n", paste(deparse(getCall(x$info$model)), 
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
summary.pmtree <- function(object, node = NULL, ...) {
  
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
  cl <- getCall(object$info$model)
  
  ## nobs
  nobs <- sapply(info, function(x) x$nobs)
  nullnobs <- sapply(nobs, is.null)
  if(1 %in% ids & any(nullnobs)) 
    nobs[nullnobs] <- object$nobs
  names(nobs) <- paste("node", ids)
  
  
  
  ret <- list(ids = ids, call = cl, coefs = coefs, objfs = objfs, nobs = nobs)
  class(ret) <- "summary.pmtree"
  return(ret)
}


#' @rdname print.pmtree 
#' 
#' @export
print.summary.pmtree <- function(x, digits = 4, ...) {
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
