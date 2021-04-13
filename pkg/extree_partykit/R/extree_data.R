
## dim method
dim.extree_data <- function(x) {
    dim(x$data)
}

## dimnames method
dimnames.extree_data <- function(x) {
    dimnames(x$data)
}

## str method
str.extree_data <- function(object, max.level = 1, give.attr = FALSE, ...) {
    cat("'extree_data':\n")
    str(unclass(object), max.level = max.level, give.attr = give.attr, ...)
}


## print method
print.extree_data <- function(x, maxvars = 5, ...) { 
  
    cat("'extree_data' with data of dimension", paste(dim(x$data), collapse = "x"), "\n")
    yx_inum <- extree_variable(x, variable = "yx", type = "inum")
    if(class(yx_inum) == "inumtotal") 
      cat("  dimension of yx reduced to", paste(dim(attr(yx_inum, "levels")), collapse = "x"), "\n")
    cat("\n")
    
    ## paste first 5 variables
    pastevars <- function(vars){
        
        lengthvars <- ifelse(is.logical(vars), sum(vars), length(vars))
        
        if(lengthvars > maxvars) {
            return(paste0(paste(names(x$data)[vars][1:maxvars], collapse = ", "), ", ..."))
        } else {
            return(paste(names(x$data)[vars], collapse = ", "))
        } 
    }
        
    cat("Variables: \n")
    cat(length(x$variables$y), "response variable(s) y:", pastevars(x$variables$y), "\n")
    if(length(x$variables$x) > 0) 
        cat(length(x$variables$x), "model variable(s) x:", pastevars(x$variables$x), "\n")
    cat(sum(x$variables$z), "split variable(s) z:", pastevars(x$variables$z == 1), "\n")
}

## FIXME (HS) summary method: how many variables of which type, which binning, scores?



## extensible tree (model) function
extree_data <- function(formula, data, subset, na.action = na.pass, weights, offset, cluster,
    strata, scores = NULL, yx = c("none", "matrix"), ytype = c("vector", "data.frame", "matrix"), 
    nmax = c("yx" = Inf, "z" = Inf), ...)
{
    ## call
    cl <- match.call()
    yx <- match.arg(yx, choices = c("none", "matrix"))
    ytype <- match.arg(ytype, choices = c("vector", "data.frame", "matrix"))
    
    if (!all(names(nmax) %in% c("yx", "z"))) {
      stop("names of 'nmax' should be equal to 'yx' and 'z'")
    }
    
    ## FIXME: (SD) Currently only NULL, numeric vector but not a matrix is 
    ## allowed as offset.
    if (!missing(offset) && !is.vector(offset)) {
      stop("unsupported specification of 'offset', should be a numeric vector")
    }
    
    ## 'formula' may either be a (multi-part) formula or a list
    noformula <- !inherits(formula, "formula")
    if(noformula) {
        
        ## formula needs to be a 'list' (if it is not a 'formula')
        if(!inherits(formula, "list")) stop("unsupported specification of 'formula'")   
        
        ## elements of formula needs to be variable names of data 
        ischar <- sapply(formula, is.character)
        indata <- sapply(formula, function(nam) {all(nam %in% names(data))})
        if (!all(ischar)) stop("elements in 'formula' need to be of type character")
        if (!all(indata)) stop(sprintf("variable '%s' not found in 'data'", names(indata[!indata])))
        
        ## specified formula elements and overall call elements
        fonam <- names(formula)
        clnam <- names(cl)[-1L]
        vanam <- c("y", "x", "z", "weights", "offset", "cluster", "strata")
        
        ## y and z (and optionally x) need to be in formula
        if(!all(c("y", "z") %in% fonam)) stop("'formula' needs to specify at least a response 'y' and partitioning variables 'z'")
        if(!("x" %in% fonam)) formula$x <- NULL 
        
        ## furthermore weights/offset/cluster/strata may be in formula or call
        vars <- formula[vanam]
        names(vars) <- vanam
        if("weights" %in% clnam) {
            clvar <- try(weights, silent = TRUE)
            vars[["weights"]] <- c(vars[["weights"]], if(!inherits(clvar, "try-error")) clvar else deparse(cl$weights))
        }
        if("offset" %in% clnam) {
            clvar <- try(offset, silent = TRUE)
            if(!inherits(clvar, "try-error")) {
            ## Added this to handle numeric offsets properly 
              if (is.numeric(offset)) {
                data[,"(offset)"] <- offset
                clvar <- "(offset)"
              } 
            } else {
              clvar <- deparse(cl$offset)
            }
            vars[["offset"]] <- c(vars[["offset"]], clvar)
        }
        if("cluster" %in% clnam) {
            clvar <- try(cluster, silent = TRUE)
            vars[["cluster"]] <- c(vars[["cluster"]], if(!inherits(clvar, "try-error")) clvar else deparse(cl$cluster))
        }
        if("strata" %in% clnam) {
            clvar <- try(strata, silent = TRUE)
            vars[["strata"]] <- c(vars[["strata"]], if(!inherits(clvar, "try-error")) clvar else deparse(cl$strata))
        }
        
        if ("offset" %in% c(fonam, clnam)) {
          ## sum offset 
          data[, "(offset)"] <- rowSums(data[vars[["offset"]]])
          vars[["offset"]] <- "(offset)" 
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
        mf$na.action <- na.action ### evaluate na.action
        if(missing(data)) data <- environment(formula)
        m <- match(c("formula", "data", "subset", "na.action", "weights", "offset", "cluster", "strata"), names(mf), 0L)
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
        offset <- model.offset(mf)
        
        mf$"(offset)" <- offset
        
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
            y = .get_term_labels(mt$y),
            x = unique(unlist(lapply(grep("^x", names(mt)), function(i) .get_term_labels(mt[[i]])))),
            z = .get_term_labels(mt$z),
            weights = if("(weights)" %in% names(mf)) "(weights)" else NULL,
            offset  = if("(offset)"  %in% names(mf)) "(offset)"  else NULL,
            cluster = if("(cluster)" %in% names(mf)) "(cluster)" else NULL,
            strata = if("(strata)" %in% names(mf)) "(strata)" else NULL
        )
        ymult <- length(vars$y) >= 1L
        if(!ymult) vars$y <- names(mf)[1L]
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
    for(v in c("y", "x", "weights", "offset", "cluster", "strata")) {
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
    
    
    ## FIXME: subsequently fitting, testing, splitting
    ## - fit: either pre-processed _and_ subsetted data --or-- full data object plus subset vector
    ## - test: additionally needs fit output --and-- fit function
    ## - split: additionally needs test output
    ## - tbd: control of all details
    
    ret <- list(
        data = if(noformula) data else mf,
        model.frame = if(noformula) FALSE else TRUE,
        variables = vars,
        terms = mt,
        scores = NULL,
        zindex = NULL,
        missings = NULL,
        yxmissings = NULL,
        yx = NULL
    )
    
    mf <- ret$data
    yxvars <- c(vars$y, vars$x, vars$offset, vars$cluster)
    zerozvars <- which(vars$z == 0)
    
    # In future, allow vector of scores in case of a *single* ordered factor
    ret$scores <- vector(mode = "list", length = length(ret$variables$z))
    names(ret$scores) <- names(mf)
    if (!is.null(scores)) {
      if(!inherits(scores, "list")) stop("unsupported specification of 'scores', should be a list")  
      if(!(names(scores) %in% names(mf))) stop("names of 'scores' must match names of variables") 
      levs <- sapply(mf, nlevels)
      len <- sapply(scores, length)
      if(!all(levs[names(scores)] == len[names(scores)])) stop("number of scores in 'scores' must match number of levels of ordered factors")
      ret$scores[names(scores)] <- scores
    }
    
    if (length(nmax) == 1L) {
      if(is.null(names(nmax))) {
        nmax <- c("yx" = nmax, "z" = nmax)
      } else {
        temp <- c("yx" = Inf, "z" = Inf)
        temp[names(nmax)] <- nmax
        nmax <- temp
      }
    }
    
    # create zindex 
    ### <FIXME> make meanlevels an argument and make sure intersplit is TRUE </FIXME>
    ret$zindex <- inum::inum(mf, ignore = names(mf)[zerozvars], total = FALSE, 
        nmax = nmax["z"], meanlevels = FALSE)
    
    if (is.finite(nmax["yx"])) {
        ret$yxindex <- inum::inum(mf[, yxvars, drop = FALSE], total = TRUE, 
            as.interval = names(mf)[vars$y], complete.cases.only = TRUE, 
            nmax = nmax["yx"], meanlevels = FALSE)
        yxmf <- attr(ret$yxindex, "levels")
        yxmf[["(weights)"]] <- NULL
        attr(ret$yxindex, "levels") <- yxmf
    } else {
        ret$yxindex <- NULL
        yxmf <- mf
    }
    
    ## Get missings (only for relevant variables: yxvars and vars$z)
    get_miss <- function(i) {
        if (i %in% yxvars || vars$z[i] == 1)
            which(is.na(ret$data[, i])) else NA
    }
    ret$missings <- lapply(seq_len(ncol(ret$data)), get_miss)
    names(ret$missings) <- names(ret$data)
    ret$yxmissings <- sort(unique(do.call("c", ret$missings[yxvars])))
    
    ## FIXME: separate object with options for: discretization, condensation, some NA handling
    ## below is just "proof-of-concept" implementation using plain model.matrix() which could
    ## be included as one option...
    if (yx == "matrix") {
        ## fake formula/terms if necessary
        formula <- Formula::as.Formula(sprintf("%s ~ %s | %s",
            paste(vanam[vars$y], collapse = " + "),
            if(length(vars$x) > 0L) paste(vanam[vars$x], collapse = " + ") else "0",
            paste(vanam[vars$z > 0], collapse = " + ")
        ))
        ## SD: "all", "z" are commented out for now since they are not used afterwards
        mt <- list(
            # "all" = terms(formula),
            "y"   = terms(formula, data = data, rhs = 0L),
            # "z"   = terms(formula, data = data, lhs = 0L, rhs = 2L),
            "x"   = terms(formula, data = data, lhs = 0L, rhs = 1L),
            "yx"  = terms(formula, data = data, rhs = 1L)
        )
        ymult <- length(vars$y) > 1L
        npart <- 2L
        
        if (ytype == "vector" && !ymult) {
            yx <- list("y" = yxmf[, vanam[vars$y], drop = TRUE])
        } else if (ytype == "data.frame") {
            yx <- list("y" = yxmf[vanam[vars$y]])
        } else { ### ytype = "matrix"
            Ytmp <- model.matrix(~ 0 + ., Formula::model.part(formula, yxmf, lhs = TRUE))
            ### <FIXME> are there cases where Ytmp already has missings? </FIXME>
            if (is.finite(nmax["yx"])) {
              Ymat <- Ytmp
            } else {
                if (length(ret$yxmissings) == 0) {
                    Ymat <- Ytmp
                } else {
                  ### FIXME: (SD) why are NAs set to 0? why whole line set to 0 in case of multiple outcomes?
                    Ymat <- matrix(0, nrow = NROW(yxmf), ncol = NCOL(Ytmp))
                    Ymat[-ret$yxmissings,] <- Ytmp
                }
            }
            yx <- list("y" = Ymat)
        }
        for(i in (1L:npart)[-npart]) {
            ni <- paste("x", if(i == 1L) "" else i, sep = "")
            ti <- if(!ymult & npart == 2L) mt$yx else mt[[ni]]
            Xtmp <- model.matrix(ti, yxmf)
            if (is.finite(nmax["yx"])) {
                Xmat <- Xtmp
            } else {
                if (length(ret$yxmissings) == 0) {
                    Xmat <- Xtmp
                } else {
                    Xmat <- matrix(0, nrow = NROW(yxmf), ncol = NCOL(Xtmp))
                    Xmat[-ret$yxmissings,] <- Xtmp
                }
            }
            yx[[ni]] <- Xmat
            if(ncol(yx[[ni]]) < 1L) {
                yx[[ni]] <- NULL
            } else {
                attr(yx[[ni]], "formula") <- formula(formula, rhs = i)
                attr(yx[[ni]], "terms") <- ti
                attr(yx[[ni]], "offset") <- yxmf[["(offset)"]]
            }    
        }
        ret$yx <- yx
    }
    
    class(ret) <- "extree_data"
    ret
}

model.frame.extree_data <- function(formula, yxonly = FALSE, ...) {
    if (!yxonly) 
        return(formula$data)
    if (!is.null(formula$yxindex))
        return(attr(formula$yxindex, "levels"))
    vars <- formula$variables
    return(formula$data[, c(vars$y, vars$x, vars$offset, vars$cluster),drop = FALSE])
}    

### <FIXME> document how to extract slots fast </FIXME>
extree_variable <- function(x, index = NULL, variable = NULL,  
  type = c("original", "inum", "scores", "missings")) {
  type <- match.arg(type, choices = c("original", "inum", "scores", "missings"))
  
  if (!is.null(variable)) {
    match.arg(variable, c("yx", "y")) # FIXME: (SD) Allow x and z? 
    if (!is.null(index)) stop("Specify either 'index' or 'variable' - not both")
  }
  
  ## Check
  if(length(index) > 1 & type == "scores") 
    stop("For length(index) > 1 type = 'scores' is currently not implemented.")
  
  switch(type, 
    "original" = {
      if (!is.null(variable)) {
        switch(variable, 
          # model.frame(x$terms[[variable]], model.frame(x)) # Shortcut for original data not binned one
          "yx" = {
            return(model.frame(x, yxonly = TRUE)) # FIXME: SD is attr(formula$yxindex, "levels") in model.frame really what we want?
          }, 
          "y" = {
            # FIXME: (SD) binning: return yx? 
            # or returned binned y: mf <- model.frame(x$terms$y, attr(x$yxindex, "levels"))
            # or return original? (currently, see also above)
            mf <- model.frame(x)
            yid <- x$variables$y
            if (length(yid) == 1) return(mf[[yid]]) else return(mf[yid])
          })
      }
      mf <- model.frame(x)
      ### [[.data.frame needs lots of memory
      class(mf) <- "list"
      if (length(index) == 1) return(mf[[index]]) else return(mf[index])
    },
    "inum" = {
      if (!is.null(variable)) {
        var <- x$yxindex
        switch(variable, 
          "yx" = {
            return(var)
          }, 
          "y" = {
            mf <- model.frame(x)
            varnam <- names(mf)[x$variables$y]
            attr(var, "levels") <- attr(var, "levels")[, varnam, drop = FALSE]
            return(var)
          }
        )
      }
      if (length(index) == 1) return(x$zindex[[index]]) else return(x$zindex[index])
    },
    "scores" = {
      # f <- x[[index]]
      f <- extree_variable(x, index, variable, type = "original")
      if (is.ordered(f)) {
        sc <- x$scores[[index]]
        if (is.null(sc)) sc <- 1:nlevels(f)
        return(sc)
      }
      return(NULL)
    },
    "missings" = {
      if (!is.null(variable)) {
        switch(variable, 
          "yx" = {
            return(x$yxmissings)
          }, 
          "y" = {
            yid <- x$variables$y
            if (length(yid) == 1) return(x$missings[[yid]]) else return(x$missings[yid])
          })
      } 
      if (length(index) == 1) return(x$missings[[index]]) else return(x$missings[index])
    }
  )
}

## for handling of non-standard variable names within extree_data() by mimicking 
## handling of such variables in model.frame() etc. in "stats" prompted by
## https://stackoverflow.com/questions/64660889/ctree-ignores-variables-with-non-syntactic-names

.deparse_variables <- function(x) paste(deparse(x, width.cutoff = 500L,
  backtick = !is.symbol(x) && is.language(x)), collapse = " ")

.get_term_labels <- function(terms, delete_response = FALSE) {
  ## ## with just standard variable names one could use:
  ## attr(terms, "term.labels")

  ## delete response from terms (if needed)
  if(delete_response) terms <- delete.response(terms)

  ## with non-standard variable names -> deparse to handle `...` correctly
  vapply(attr(terms, "variables"), .deparse_variables, " ")[-1L]
}

### <FIXME> (HS) delete if possible </FIXME>
"[[.extree_data" <- function(x, i, 
    type = c("original", "inum", "scores", "missings")) { 
    
    warning("[[.extree_data is deprecated. Please use extree_variable().")
    extree_variable(x = x, i = i, type = type)
}
