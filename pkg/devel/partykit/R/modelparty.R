
## methods concerning call/formula/terms/etc.
## (default methods work for terms and update)

formula.modelparty <- function(x, extended = FALSE, ...)
  if(extended) x$info$Formula else x$info$formula

getCall.modelparty <- function(x, ...) x$info$call

model.frame.modelparty <- function(formula, ...)
{
  mf <- formula$data
  if(nrow(mf) > 0L) return(mf)

  dots <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0L)]
  mf <- formula$info$call
  mf <- mf[c(1L, match(c("formula", "data", "subset", "na.action"), names(mf), 0L))]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf[names(nargs)] <- nargs
  if(is.null(env <- environment(formula$info$terms))) env <- parent.frame()
  mf$formula <- Formula::Formula(as.formula(mf$formula))
  eval(mf, env)
}

weights.modelparty <- function(object, ...) {
  fit <- object$fitted
  ww <- if(!is.null(w <- fit[["(weights)"]])) w else rep.int(1L, NROW(fit))
  structure(ww, .Names = rownames(fit))
}

## methods concerning model/parameters/loglik/etc.

coef.modelparty <- function(object, node = NULL, drop = TRUE, ...) {
  if(is.null(node)) node <- nodeids(object, terminal = TRUE)
  cf <- do.call("rbind", nodeapply(object, ids = node, FUN = function(n) info_node(n)$coefficients))
  if(drop) drop(cf) else cf
}

refit.modelparty <- function(object, node = NULL, drop = TRUE, ...)
{
  ## by default use all ids
  if(is.null(node)) node <- nodeids(object)
  
  ## estimated coefficients
  cf <- nodeapply(object, ids = node, FUN = function(n) info_node(n)$coefficients)
  
  ## model.frame
  mf <- model.frame(object)
  weights <- weights(object)
  offset <- model.offset(mf)
  cluster <- mf[["(cluster)"]]
  
  ## fitted ids
  fitted <- object$fitted[["(fitted)"]]

  ## response variables
  Y <- switch(object$info$control$ytype,
    "vector" = Formula::model.part(object$info$Formula, mf, lhs = 1L)[[1L]],
    "matrix" = model.matrix(~ 0 + ., Formula::model.part(object$info$Formula, mf, lhs = 1L)),
    "data.frame" = Formula::model.part(object$info$Formula, mf, lhs = 1L)
  )
  hasx <- object$info$nreg >= 1L | attr(object$info$terms$response, "intercept") > 0L
  X <- if(!hasx) NULL else switch(object$info$control$xtype,
    "matrix" = model.matrix(object$info$terms$response, mf),
    "data.frame" = Formula::model.part(object$info$Formula, mf, rhs = 1L)
  )
  if(!is.null(X)) {
    attr(X, "formula") <- formula(object$info$Formula, rhs = 1L)
    attr(X, "terms") <- object$info$terms$response
    attr(X, "offset") <- object$info$call$offset
  }

  suby <- function(y, index) {
    if(object$info$control$ytype == "vector") y[index] else y[index, , drop = FALSE]
  }
  subx <- if(hasx) {
    function(x, index) {
      sx <- x[index, , drop = FALSE]
      attr(sx, "contrasts") <- attr(x, "contrasts")
      attr(sx, "xlevels")   <- attr(x, "xlevels")
      attr(sx, "formula")   <- attr(x, "formula")
      attr(sx, "terms")     <- attr(x, "terms")
      attr(sx, "offset")    <- attr(x, "offset")
      sx
    }
  } else {
    function(x, index) NULL
  }
  
  ## fitting function
  afit <- object$info$fit
  
  ## refit
  rval <- lapply(seq_along(node), function(i) {
    ix <- fitted %in% nodeids(object, from = node[i], terminal = TRUE)
    args <- list(y = suby(Y, ix), x = subx(X, ix), start = cf[[i]],
      weights = weights[ix], offset = offset[ix], cluster = cluster[ix], object = TRUE)
    args <- c(args, object$info$dots)
    do.call("afit", args)$object
  })
  names(rval) <- node

  ## drop?
  if(drop & length(rval) == 1L) rval <- rval[[1L]]
  
  ## return
  return(rval)
}

apply_to_models <- function(object, node = NULL, FUN = NULL, drop = FALSE, ...) {
  if(is.null(node)) node <- nodeids(object, terminal = FALSE)
  if(is.null(FUN)) FUN <- function(object, ...) object  
  rval <- if("object" %in% object$info$control$terminal) {
    nodeapply(object, node, function(n) FUN(info_node(n)$object))
  } else {
    lapply(refit.modelparty(object, node, drop = FALSE), FUN)
  }
  names(rval) <- node
  if(drop & length(node) == 1L) rval <- rval[[1L]]
  return(rval)
}

logLik.modelparty <- function(object, dfsplit = NULL, ...)
{
  if(is.null(dfsplit)) dfsplit <- object$info$control$dfsplit
  dfsplit <- as.integer(dfsplit)
  ids <- nodeids(object, terminal = TRUE)
  ll <- apply_to_models(object, node = ids, FUN = logLik)
  dfsplit <- dfsplit * (length(object) - length(ll))
  structure(
    sum(as.numeric(ll)),
    df = sum(sapply(ll, function(x) attr(x, "df"))) + dfsplit,
    nobs = nobs(object),
    class = "logLik"
  )
}

nobs.modelparty <- function(object, ...) {
  sum(unlist(nodeapply(object,
    nodeids(object, terminal = TRUE),
    function(n) info_node(n)$nobs
  )))
}

deviance.modelparty <- function(object, ...)
{
  ids <- nodeids(object, terminal = TRUE)
  dev <- apply_to_models(object, node = ids, FUN = deviance)
  sum(unlist(dev))
}

summary.modelparty <- function(object, node = NULL, ...)
{
  ids <- if(is.null(node)) nodeids(object, terminal = TRUE) else node
  rval <- apply_to_models(object, node = ids, FUN = summary)
  if(length(ids) == 1L) rval[[1L]] else rval
}

sctest.modelparty <- function(object, node = NULL, ...)
{
  ids <- if(is.null(node)) nodeids(object, terminal = FALSE) else node
  rval <- nodeapply(object, ids, function(n) info_node(n)$test)
  names(rval) <- ids
  if(length(ids) == 1L) rval[[1L]] else rval  
}

print.modelparty <- function(x, node = NULL,
  FUN = NULL, digits = getOption("digits") - 4L,
  header = TRUE, footer = TRUE, title = NULL, objfun = "", ...)
{
  digits <- max(c(0, digits))
  if(objfun != "") objfun <- paste(" (", objfun, ")", sep = "")
  if(is.null(title)) title <- sprintf("Model-based recursive partitioning (%s)",
    deparse(x$info$call$fit))

  if(is.null(node)) {
    header_panel <- if(header) function(party) {      
      c(title, "", "Model formula:", deparse(party$info$formula), "", "Fitted party:", "")
    } else function(party) ""
  
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
        paste("Objective function", objfun, ": ", of, sep = ""), "")
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
      FUN = function(n) info_node(n)[c("coefficients", "objfun", "test")])    
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
      print(info[[i]]$test)
    }
  }
  invisible(x)
}

predict.modelparty <- function(object, newdata = NULL, type = "node", ...)
{
  ## predicted node ids
  node <- predict.party(object, newdata = newdata)
  if(identical(type, "node")) return(node)

  ## obtain fitted model objects
  ids <- sort(unique(as.integer(node)))
  mod <- apply_to_models(object, node = ids)

  ## obtain predictions
  pred <- if(is.character(type)) {
    function(object, newdata = NULL, ...) predict(object, newdata = newdata, type = type, ...)
  } else {
    type
  }
  if("newdata" %in% names(formals(pred))) {    
    ix <- lapply(seq_along(ids), function(i) which(node == ids[i]))
    preds <- lapply(seq_along(ids), function(i)
      pred(mod[[i]], newdata = newdata[ix[[i]], , drop = FALSE], ...))
    nc <- NCOL(preds[[1L]])
    rval <- if(nc > 1L) {
      matrix(0, nrow = length(node), ncol = nc, dimnames = list(names(node), colnames(preds[[1L]])))
    } else {
      rep(preds[[1L]], length.out = length(node))
    }      
    for(i in seq_along(ids)) {
      if(nc > 1L) {
        rval[ix[[i]], ] <- preds[[i]]
	rownames(rval) <- names(node)
      } else {
        rval[ix[[i]]] <- preds[[i]]
	names(rval) <- names(node)
      }
    }
  } else {
    rval <- lapply(mod, pred, ...)
    if(NCOL(rval[[1L]]) > 1L) {
      rval <- do.call("rbind", rval)
      rownames(rval) <- ids
      rval <- rval[as.character(node), , drop = FALSE]
      rownames(rval) <- names(node)
      rval <- drop(rval)
    } else {
      ## provide a c() method for factors locally
      c.factor <- function(...) {
        args <- list(...)
	lev <- levels(args[[1L]])
	args[[1L]] <- unclass(args[[1L]])
	rval <- do.call("c", args)
	factor(rval, levels = 1L:length(lev), labels = lev)
      }
      rval <- do.call("c", rval)
      names(rval) <- ids
      rval <- rval[as.character(node)]
      names(rval) <- names(node)
    }
  }
  return(rval)
}

fitted.modelparty <- function(object, ...)
{
  ## fitted nodes
  node <- predict.party(object, type = "node")

  ## obtain fitted model objects
  ids <- nodeids(object, terminal = TRUE)
  fit <- apply_to_models(object, node = ids, FUN = fitted)

  nc <- NCOL(fit[[1L]])
  rval <- if(nc > 1L) {
    matrix(0, nrow = length(ids), ncol = nc, dimnames = list(names(node), colnames(fit[[1L]])))
  } else {
    rep(fit[[1L]], length.out = length(node))
  }	 
  for(i in seq_along(ids)) {
    if(nc > 1L) {
      rval[node == ids[i], ] <- fit[[i]]
      rownames(rval) <- names(node)
    } else {
      rval[node == ids[i]] <- fit[[i]]
      names(rval) <- names(node)
    }
  }

  return(rval)
}

residuals.modelparty <- function(object, ...)
{
  ## fitted nodes
  node <- predict.party(object, type = "node")

  ## obtain fitted model objects
  ids <- nodeids(object, terminal = TRUE)
  res <- apply_to_models(object, node = ids, FUN = residuals)

  nc <- NCOL(res[[1L]])
  rval <- if(nc > 1L) {
    matrix(0, nrow = length(ids), ncol = nc, dimnames = list(names(node), colnames(res[[1L]])))
  } else {
    rep(res[[1L]], length.out = length(node))
  }	 
  for(i in seq_along(ids)) {
    if(nc > 1L) {
      rval[node == ids[i], ] <- res[[i]]
      rownames(rval) <- names(node)
    } else {
      rval[node == ids[i]] <- res[[i]]
      names(rval) <- names(node)
    }
  }

  return(rval)
}

plot.modelparty <- function(x, terminal_panel = NULL, FUN = NULL, ...) {
  if(is.null(terminal_panel)) {
    if(is.null(FUN)) {
      FUN <- function(x) {
        cf <- x$coefficients
	cf <- matrix(cf, ncol = 1, dimnames = list(names(cf), ""))
        c(sprintf("n = %s", x$nobs), "Estimated parameters:",
          strwrap(capture.output(print(cf, digits = 4L))[-1L]))
      }
    }
    terminal_panel <- node_terminal(x, FUN = FUN)
  }
  plot.party(x, terminal_panel = terminal_panel, ...)
}
