## extensible tree (model) function
extree <- function(formula, data, subset, na.action = na.pass, weights, offset, cluster,
  yx = "none", ...)
{
  ## call
  cl <- match.call()

  ## 'formula' may either be a (multi-part) formula or a list
  noformula <- !inherits(formula, "formula")
  if(noformula) {

    ## formula needs to be a 'list' (if it is not a 'formula')
    if(!inherits(formula, "list")) stop("unsupported specification of 'formula'")    

    ## specified formula elements and overall call elements
    fonam <- names(formula)
    clnam <- names(cl)[-1L]
    vanam <- c("y", "x", "z", "weights", "offset", "cluster")
    
    ## y and z (and optionally x) need to be in formula
    if(!all(c("y", "z") %in% fonam)) stop("'formula' needs to specify at least a response 'y' and partitioning variables 'z'")
    if(!("x" %in% fonam)) formula$x <- NULL
    
    ## furthermore weights/offset/cluster may be in formula or call
    vars <- formula[vanam]
    names(vars) <- vanam
    if("weights" %in% clnam) {
      clvar <- try(weights, silent = TRUE)
      vars[["weights"]] <- c(vars[["weights"]], if(!inherits(clvar, "try-error")) clvar else deparse(cl$weights))
    }
    if("offset" %in% clnam) {
      clvar <- try(offset, silent = TRUE)
      vars[["offset"]] <- c(vars[["offset"]], if(!inherits(clvar, "try-error")) clvar else deparse(cl$offset))
    }
    if("cluster" %in% clnam) {
      clvar <- try(cluster, silent = TRUE)
      vars[["cluster"]] <- c(vars[["cluster"]], if(!inherits(clvar, "try-error")) clvar else deparse(cl$cluster))
    }
    
    ## sanity checking
    for(v in vanam) {
      if(!is.null(vars[[v]]) && !(is.numeric(vars[[v]]) | is.character(vars[[v]]) | is.logical(vars[[v]]))) {
        warning(sprintf("unknown specification of '%s', must be character, numeric, or logical", v))
        vars[v] <- list(NULL)
      }
    }
    if(!missing(subset)) warning("'subset' argument ignored in list specification of 'formula'")
    if(!missing(na.action)) warning("'na.action' argument ignored in list specification of 'formula'")    

    ## no terms (by default)
    mt <- NULL

  } else {

    ## set up model.frame() call
    mf <- match.call(expand.dots = FALSE)
    if(missing(data)) data <- environment(formula)
    m <- match(c("formula", "data", "subset", "na.action", "weights", "offset", "cluster"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf$dot <- "sequential"

    ## formula processing
    oformula <- as.formula(formula)
    formula <- Formula::as.Formula(formula)
    mf$formula <- formula
    npart <- length(formula)
    if(any(npart < 1L)) stop("'formula' must specify at least one left-hand and one right-hand side")
    npart <- npart[2L]

    ## evaluate model.frame
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())

    ## extract terms in various combinations
    mt <- list(
      "all" = terms(formula, data = data,                        dot = "sequential"),
      "y"   = terms(formula, data = data, rhs = 0L,              dot = "sequential"),
      "z"   = terms(formula, data = data, lhs = 0L, rhs = npart, dot = "sequential")
    )
    if(npart > 1L) {
      mt$yx <-terms(formula, data = data, rhs = 1L,              dot = "sequential")
      for(i in 1L:(npart-1L)) {
        mt[[paste("x", if(i == 1L) "" else i, sep = "")]] <- terms(
	            formula, data = data, lhs = 0L, rhs = i,     dot = "sequential")
      }
    }

    ## extract variable lists
    vars <- list(
      y = attr(mt$y, "term.labels"),
      x = unique(unlist(lapply(grep("^x", names(mt)), function(i) attr(mt[[i]], "term.labels")))),
      z = attr(mt$z, "term.labels"),
      weights = if("(weights)" %in% names(mf)) "(weights)" else NULL,
      offset  = if("(offset)"  %in% names(mf)) "(offset)"  else NULL,
      cluster = if("(cluster)" %in% names(mf)) "(cluster)" else NULL
    )
    ymult <- length(vars$y) >= 1L
    if(!ymult) vars$y <- names(mf)[1L]
    ## FIXME: store information which variable(s) went into (weights), (offset), (cluster)
    ## idea: check (x and) z vs. deparse(cl$weights), deparse(cl$offset), deparse(cl$cluster)

    ## check wether offset was inside the formula
    if(!is.null(off <- attr(mt$x, "offset"))) {
      if(is.null(vars$offset)) mf[["(offset)"]] <- rep.int(0, nrow(mf))
      for(i in off) mf[["(offset)"]] <- mf[["(offset)"]] + mf[[i]]
      vars$offset <- "(offset)"
    }
  }
  
  ## canonicalize y/x/z term labels
  vanam <- if(noformula) names(data) else names(mf)
  ## z to numeric
  if(is.null(vars$z)) stop("at least one 'z' variable must be specified")
  if(is.integer(vars$z)) vars$z <- vanam[vars$z]
  if(is.character(vars$z)) vars$z <- vanam %in% vars$z
  if(is.logical(vars$z)) vars$z <- as.numeric(vars$z)
  if(is.null(names(vars$z))) names(vars$z) <- vanam
  vars$z <- vars$z[vanam]
  if(any(is.na(vars$z))) vars$z[is.na(vars$z)] <- 0
  vars$z <- as.numeric(vars$z)
  ## all others to integer
  for(v in c("y", "x", "weights", "offset", "cluster")) {
    if(!is.null(vars[[v]])) {
      if(is.character(vars[[v]])) vars[[v]] <- match(vars[[v]], vanam)
      if(is.logical(vars[[v]])) vars[[v]] <- which(vars[[v]])
      if(any(is.na(vars[[v]]))) {
        vars[[v]] <- vars[[v]][!is.na(vars[[v]])]
        warning(sprintf("only found the '%s' variables: %s", v, paste(vanam[vars[[v]]], collapse = ", ")))
      }
    }
    vars[[v]] <- unique(as.integer(vars[[v]]))
  }
  if(is.null(vars$y)) stop("at least one 'y' variable must be specified")


  ## FIXME: separate object with options for: discretization, condensation, some NA handling
  ## below is just "proof-of-concept" implementation using plain model.matrix() which could
  ## be included as one option...
  if(yx == "matrix") {
  
    ## fake formula/terms if necessary
    if(noformula) {
      formula <- Formula::as.Formula(sprintf("%s ~ %s | %s",
        paste(vanam[vars$y], collapse = " + "),
        if(length(vars$x) > 0L) paste(vanam[vars$x], collapse = " + ") else "0",
        paste(vanam[vars$z > 0], collapse = " + ")
      ))
      mt <- list(
        "all" = terms(formula),
        "y"   = terms(formula, data = data, rhs = 0L),
        "z"   = terms(formula, data = data, lhs = 0L, rhs = 2L),
        "yx"  = terms(formula, data = data, rhs = 1L),
        "x"   = terms(formula, data = data, lhs = 0L, rhs = 1L)
      )
      ymult <- length(vars$y) > 1L
      npart <- 2L
    }

    yx <- list("y" = model.matrix(~ 0 + ., Formula::model.part(formula, if(noformula) data else mf, lhs = TRUE)))         
    for(i in (1L:npart)[-npart]) {
      ni <- paste("x", if(i == 1L) "" else i, sep = "")
      ti <- if(!ymult & npart == 2L) mt$yx else mt[[ni]]
      yx[[ni]] <- model.matrix(ti, if(noformula) data else mf)
      if(ncol(yx[[ni]]) < 1L) {
        yx[[ni]] <- NULL
      } else {
        attr(yx[[ni]], "formula") <- formula(formula, rhs = i)
        attr(yx[[ni]], "terms") <- ti
        attr(yx[[ni]], "offset") <- cl$offset
      }    
    }

  }

  ## FIXME: subsequently fitting, testing, splitting
  ## - fit: either pre-processed _and_ subsetted data --or-- full data object plus subset vector
  ## - test: additionally needs fit output --and-- fit function
  ## - split: additionally needs test output
  ## - tbd: control of all details

  return(list(
    data = if(noformula) data else mf,
    variables = vars,
    terms = mt,
    yx = yx
  ))
}
