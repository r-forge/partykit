
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
print.extree_data <- function(x, ...) { 
    # FIXME (HS) add argument maxvars to let user choose the number of printed vars
    cat("'extree_data' with data of dimension", paste(dim(x$data), collapse = "x"), "\n\n")
    
    ## paste first 5 variables
    pastevars <- function(vars){
        
        lengthvars <- ifelse(is.logical(vars), sum(vars), length(vars))
        
        if(lengthvars > 5) {
            return(paste0(paste(names(x$data)[vars][1:5], collapse = ", "), ", ..."))
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
    
    ## 'formula' may either be a (multi-part) formula or a list
    noformula <- !inherits(formula, "formula")
    if(noformula) {
        
        ## formula needs to be a 'list' (if it is not a 'formula')
        if(!inherits(formula, "list")) stop("unsupported specification of 'formula'")    
        
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
            vars[["offset"]] <- c(vars[["offset"]], if(!inherits(clvar, "try-error")) clvar else deparse(cl$offset))
        }
        if("cluster" %in% clnam) {
            clvar <- try(cluster, silent = TRUE)
            vars[["cluster"]] <- c(vars[["cluster"]], if(!inherits(clvar, "try-error")) clvar else deparse(cl$cluster))
        }
        if("strata" %in% clnam) {
            clvar <- try(strata, silent = TRUE)
            vars[["strata"]] <- c(vars[["strata"]], if(!inherits(clvar, "try-error")) clvar else deparse(cl$strata))
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
            cluster = if("(cluster)" %in% names(mf)) "(cluster)" else NULL,
            strata = if("(strata)" %in% names(mf)) "(strata)" else NULL
        )
        ymult <- length(vars$y) >= 1L
        if(!ymult) vars$y <- names(mf)[1L]
        ## FIXME: store information which variable(s) went into (weights), (offset), (cluster)
        ## (strata)
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
        variables = vars,
        terms = mt
    )
    
    mf <- ret$data
    yxvars <- c(vars$y, vars$x, vars$offset, vars$cluster)
    zerozvars <- which(vars$z == 0)
    
    ret$scores <- vector(mode = "list", length = length(ret$variables$z))
    names(ret$scores) <- names(mf)
    if (!is.null(scores))
        ret$scores[names(scores)] <- scores
    
    if (length(nmax) == 1L) nmax <- c("yx" = nmax, "z" = nmax)
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
          ## FIXME: (SD) model.part not practical for high-dim data  
            Ytmp <- model.matrix(~ 0 + ., Formula::model.part(formula, yxmf, lhs = TRUE))
            ### <FIXME> are there cases where Ytmp already has missings? </FIXME>
            if (length(ret$yxmissings) == 0) {
                Ymat <- Ytmp
            } else {
              ### FIXME: (SD) why are NAs set to 0? why whole line set to 0 in case of multiple outcomes?
                Ymat <- matrix(0, nrow = NROW(yxmf), ncol = NCOL(Ytmp))
                Ymat[-ret$yxmissings,] <- Ytmp
            }
            yx <- list("y" = Ymat)
        }
        for(i in (1L:npart)[-npart]) {
            ni <- paste("x", if(i == 1L) "" else i, sep = "")
            ti <- if(!ymult & npart == 2L) mt$yx else mt[[ni]]
            yx[[ni]] <- model.matrix(ti, yxmf)
            ## FIXME: (SD)  if model.matrix not right dim consider code below. 
            # if (length(ret$yxmissings) == 0) {
            #     Xmat <- Xtmp
            # } else {
            #     Xmat <- matrix(0, nrow = NROW(yxmf), ncol = NCOL(Xtmp))
            #     ## FIXME: (SD) Error if no missings in x, but in y! 
            #     ## Either do not replace NAs by 0 or do it with is.NA() - not ret$xymissing! 
            #     Xmat[-ret$yxmissings,] <- Xtmp
            # }
            # yx[[ni]] <- Xmat
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
extree_variable <- function(x, i, 
    type = c("original", "index", "scores", "missings")) {
    type <- match.arg(type, choices = c("original", "index", "scores", "missings"))
    
    ## Check
    if(length(i) > 1 & type == "scores") 
        stop("For length(i) > 1 type = 'scores' is currently not implemented.")
    
    switch(type, 
        "original" = {
            if (length(i) == 1 && i == "yx") return(model.frame(x, yxonly = TRUE))
            mf <- model.frame(x)
            ### [[.data.frame needs lots of memory
            class(mf) <- "list"
            if (length(i) == 1) return(mf[[i]]) else return(mf[i])
        },
        "index" = {
            if (i == "yx" || (length(i) == 1 & i %in% c(x$variables$y, x$variables$x)))
            # if (i == "yx" || i %in% c(x$variables$y, x$variables$x))
                return(x$yxindex) ### may be NULL
            if (length(i) == 1) return(x$zindex[[i]]) else return(x$zindex[i])
        },
        "scores" = {
            # f <- x[[i]]
            f <- extree_variable(x, i, type = "original")
            if (is.ordered(f)) {
                sc <- x$scores[[i]]
                if (is.null(sc)) sc <- 1:nlevels(f)
                return(sc)
            }
            return(NULL)
        },
        "missings" = {
            if (i == "yx" || i %in% c(x$variables$y, x$variables$x))
                return(x$yxmissings)
            x$missings[[i]]
        }
    )
}


### <FIXME> (HS) delete if possible </FIXME>
"[[.extree_data" <- function(x, i, 
    type = c("original", "index", "scores", "missings")) { 
    
    warning("[[.extree_data is deprecated. Please use extree_variable().")
    extree_variable(x = x, i = i, type = type)
}
