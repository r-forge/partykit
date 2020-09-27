extree <- function(data, 
    trafo,
    control = extree_control(#TODO: rename selectfun -> var_select, splitfun -> split_select, 
        ...), 
    converged = NULL,
    ...) {
    
    ## check / preprocess extree data
    subset <- .start_subset(data = data)
    weights <- model.weights(model.frame(data))
    
    
    
    
    ## trafo preprocessing
    mytrafo <- function(subset, weights, info = NULL, estfun = TRUE, object = TRUE) {
        trafo(subset, data = data, weights, info = NULL, estfun = TRUE, object = TRUE)
    }
    
    ## TODO: converged preprocessing (if needed)
    
    ## set up trafo
    update <- function(subset, weights, control, doFit = TRUE) {
        extree_fit(data = data, trafo = mytrafo, converged = converged,
            partyvars = data$variables$z, subset = subset,
            weights = weights, ctrl = control, doFit = doFit)  
    }
    
    ## fit
    tree <- update(subset = subset, weights = weights, control = control)
    
    ## TODO: prepare extree object
    
}


## helper function for .preprocess_var_select
.get_strategy_function <- function(strategy, select_type = "var") {
    
    
    ## find matching objects and set up list
    ## FIXME: better solution than via search()?
    snam <- sprintf(paste0("^", select_type, "_select_%s"), strategy)
    onam <- unlist(lapply(search(), objects))
    onam <- unique(onam[grep(snam, onam)])
    strategy <- lapply(onam, get)
    
    ## supply variable type a name
    nam <- gsub(pattern = paste0(snam, "_"), replacement = "", onam)
    if (all(nam %in% onam)) nam <- regmatches(onam, regexpr("[a-z]+$", onam))
    names(strategy) <- nam
    
    
    ## drop non-functions
    strategy <- strategy[sapply(strategy, is.function)]
    
    return(strategy)
}




.preprocess_select <- function(select, select_type = "var") {
    
    ## function: return as is
    if (is.function(select)) return(select)
    
    ## character: return appropriate function
    if (is.character(select)) {
        return(
            .get_strategy_function(select, 
                select_type = select_type)
            )
    }
    
    ## list: go through all elements and return accordingly
    if (is.list(select)) {
        
        get_strategy <- function(select_nam) {
            
            ## return function if function
            if (is.function(select[[select_nam]])) {
                sfun <- list(select[[select_nam]])
                names(sfun) <- select_nam
                return(sfun)
            } 
            
            ## get appropriat function if character
            if (is.character(select[[select_nam]])) {
                sfun <- .get_strategy_function(select[[select_nam]], 
                    select_type = select_type)
                return(sfun)
            } 
            
            ## if none of the above -> ERROR
            stop("select can only be functions or characters.")
            
        }
        
        ## go through all list elements and choose approriate function
        select_list <- sapply(names(select), FUN = get_strategy, simplify = TRUE, USE.NAMES = FALSE)
        return(select_list)
    }
}



# <FIXME> (HS) better name of function
.get_varclass <- function(select_list, data, j) {
    ### which class is variable?
    varclass_orig <- class(extree_variable(x = data, i = j))
    
    ### Use most appropriate class (1st), if more than one is available
    ### Remove varclass if no var_select function is available
    varclass <- varclass_orig[varclass_orig %in% names(select_list)][1]
    
    ### if no function for this class is provided use default function
    if(length(varclass) == 0 | is.na(varclass)) {
        if (!("default" %in% names(select_list))) 
            stop("The is no specific or default select function for split variables of class ",
                varclass_orig, ". Please provide one.")
        varclass <- "default"
    } 
    
    return(varclass)
}

# <FIXME> (HS) Better name for function
selector <- function(select, model, trafo, data, subset, weights, whichvar, control, j) {
    
    
    # <FIXME> (HS) add check if function(s) return what we want
    
    if(is.function(select)) {
        ## if var_select is function, apply function
        ret <- select(model = model, trafo = trafo, data = data, 
            subset = subset, weights = weights, j = j, 
            split_only = FALSE, control = control)
        
    } else if (is.list(select) && all(sapply(select, is.function))) {
        
        ## if var_select is list of functions, apply appropriate function
        varclass <- .get_varclass(select_list = select, data = data, 
            j = j)
        
        ### Run appropriate var_select function
        ret <- select[[varclass]](model = model, trafo = trafo, data = data, 
            subset = subset, weights = weights, j = j, 
            split_only = FALSE, control = control)
        
    } else {
        ## a future option would be a character which hints to functions 
        ## to be used.
        stop("Selection strategy must currently be a function or a named list of functions.")
    }
    
    return(ret)
}


## new selectfunction
var_select_loop <- function(model, trafo, data, subset, weights, whichvar, 
    control, var_select) {

    ## set up return list + criterion matrix
    ret <- list(criterion = matrix(NA, nrow = 2L, ncol = ncol(model.frame(data))))
    rownames(ret$criterion) <- c("statistic", "p.value")
    colnames(ret$criterion) <- names(model.frame(data))
    if (length(whichvar) == 0) return(ret)

    ## loop over all relevant variables and use var_select function supplied
    for (j in whichvar) {
        
        tst <- selector(select = var_select, model = model, trafo = trafo, 
            data = data, subset = subset, weights = weights, j = j, 
            control = control)
        
        ret$criterion["statistic",j] <- tst$statistic
        ret$criterion["p.value",j] <- tst$p.value
    }
    ret
}

## Split function new
split_select_loop <- function(model, trafo, data, subset, weights, whichvar, 
    control, split_select) {
    
    if (length(whichvar) == 0) return(NULL)
    
    ## loop over all vars in whichvar (stop if split is found)
    for (j in whichvar) {
        
        ret <- selector(select = split_select, model = model, trafo = trafo, 
            data = data, subset = subset, weights = weights, j = j, 
            control = control)
        
        ### check if trafo can be successfully applied to all daugther nodes 
        ### (converged = TRUE)
        if (control$lookahead & !is.null(ret)) {
            sp <- kidids_split(ret, model.frame(data), obs = subset)
            conv <- sapply(unique(na.omit(sp)), function(i)
                isTRUE(trafo(subset[sp == i & !is.na(sp)], weights = weights)$converged))
            if (!all(conv)) ret <- NULL
            ## FIXME: allow option to keep estfun and move on --> update = FALSE
        }
        
        ## stop if a split was found, otherwise continue with next possible var
        if (!is.null(ret)) break()
    }
    return(ret)
}

## Select function old
.select <- function(model, trafo, data, subset, weights, whichvar, ctrl, FUN) {
    ret <- list(criterion = matrix(NA, nrow = 2L, ncol = ncol(model.frame(data))))
    rownames(ret$criterion) <- c("statistic", "p.value")
    colnames(ret$criterion) <- names(model.frame(data))
    if (length(whichvar) == 0) return(ret)
    ### <FIXME> allow joint MC in the absense of missings; fix seeds
    ### write ctree_test / ... with whichvar and loop over variables there
    ### </FIXME>
    for (j in whichvar) {
        tst <- FUN(model = model, trafo = trafo, data = data, 
                   subset = subset, weights = weights, j = j, 
                   SPLITONLY = FALSE, ctrl = ctrl)
        ret$criterion["statistic",j] <- tst$statistic
        ret$criterion["p.value",j] <- tst$p.value
    }
    ret
}

## Split function old
.split <- function(model, trafo, data, subset, weights, whichvar, ctrl, FUN) {
    if (length(whichvar) == 0) return(NULL)
    for (j in whichvar) {
        x <- model.frame(data)[[j]]
        if (ctrl$multiway && is.factor(x) && !is.ordered(x) &&
            (ctrl$maxsurrogate == 0) && nlevels(x[subset, drop = TRUE]) > 1) 
        {
            index <- 1L:nlevels(x)
            xt <- libcoin::ctabs(ix = unclass(x), weights = weights, subset = subset)[-1]
            index[xt == 0] <- NA
            ### maybe multiway is not so smart here as
            ### nodes with nobs < minbucket could result
            index[xt > 0 & xt < ctrl$minbucket] <- nlevels(x) + 1L
            if (length(unique(index)) == 1) {
                ret <- NULL
            } else {
                index <- unclass(factor(index))
                ret <- partysplit(as.integer(j), index = as.integer(index))
            }
        } else {
            ret <- FUN(model = model, trafo = trafo, data = data, 
                       subset = subset, weights = weights, j = j, 
                       SPLITONLY = TRUE, ctrl = ctrl)
        }
        ### check if trafo can be successfully applied to all daugther nodes 
        ### (converged = TRUE)
        if (ctrl$lookahead & !is.null(ret)) {
            sp <- kidids_split(ret, model.frame(data), obs = subset)
            conv <- sapply(unique(na.omit(sp)), function(i)
                    isTRUE(trafo(subset[sp == i & !is.na(sp)], weights = weights)$converged))
            if (!all(conv)) ret <- NULL
        }
        if (!is.null(ret)) break()
    }
    return(ret)
}

.objfun_select <- function(...)
    function(model, trafo, data, subset, weights, whichvar, ctrl) {
        args <- list(...)
        ctrl[names(args)] <- args
        .select(model, trafo, data, subset, weights, whichvar, ctrl, FUN = .objfun_test)
    }

.objfun_split <- function(...)
    function(model, trafo, data, subset, weights, whichvar, ctrl) {
        args <- list(...)
        ctrl[names(args)] <- args
        .split(model, trafo, data, subset, weights, whichvar, ctrl, FUN = .objfun_test) 
    }

### unbiased recursive partitioning: set up new node
.extree_node <- function
(
    id = 1L, 			### id of this node
    data, 			### full data, readonly
    trafo,
    selectfun, 			### variable selection
    splitfun,                   ### split selection
    svselectfun,                ### same for surrogate splits
    svsplitfun,                 ### same for surrogate splits
    partyvars, 			### partytioning variables
                                ### a subset of 1:ncol(model.frame(data))
    weights = integer(0L),	### optional case weights
    subset, 			### subset of 1:nrow(data)
                                ### for identifying obs for this node
    ctrl, 			### extree_control()
    info = NULL,
    extree_env = NULL           ### environment for depth and maxid
) {
  ## initialize tracker: is this a terminal node?
  terminal <- FALSE

  ## is there any info to be preserved?

  ## depth keeps track of the depth of the tree: has to be < than maxdepth
  ## maxit is the largest id in the left subtree
  if (is.null(extree_env)) {
    extree_env <- new.env()
    assign("depth", 0L, envir = extree_env)
  }
  depth <- get("depth", envir = extree_env)
  assign("maxid", id, envir = extree_env)
  if (depth >= ctrl$maxdepth) terminal <- TRUE
  ## FIXME: previously the code stopped here with:
  ##   return(partynode(as.integer(id)))
  ## However, we should collect "info" if desired!
  ## Analogously for all subsequent checks for "terminal".

  ## check for stumps
  if (!terminal && id > 1L && ctrl$stump) terminal <- TRUE

  ## sw is basically the number of observations
  ## which has to be > minsplit in order to consider
  ## the node for splitting
  ## FIXME: currently computed even if we know terminal=TRUE already
  sw <- if (length(weights) > 0L) {
    if (ctrl$caseweights) sum(weights[subset]) else sum(weights[subset] > 0L)
  } else {
    length(subset)
  }
  if (!terminal && sw < ctrl$minsplit) terminal <- TRUE

  ## FIXME: split variable selection structure should be set up
  ## with NAs if terminal=TRUE but "criterion" selected for "info"
  svars <- which(partyvars > 0)
  if (ctrl$mtry < Inf) {
    mtry <- min(length(svars), ctrl$mtry)
    svars <- .resample(svars, mtry, prob = partyvars[svars])
  } 

  ## FIXME: include a condition like
  ##  if(!terminal || "trafo" %in% ctrl$save)
  thismodel <- trafo(subset = subset, weights = weights, info = info,
    estfun = TRUE, object = TRUE)
  if (is.null(thismodel)) terminal <- TRUE

  ## update sample size constraints on possible splits
  ## need to do this here because selectfun might consider splits
  mb <- ctrl$minbucket
  mp <- ctrl$minprob
  swp <- ceiling(sw * mp)
  if (mb < swp) mb <- as.integer(swp)
  thisctrl <- ctrl
  thisctrl$minbucket <- mb

  ## split variable selection:
  ## - either already returns a finished "partysplit" object
  ## - or "criterion" matrix (typically statistics & p-value) for all variables
  ## FIXME: need to get "empty" criterion matrix with all NAs
  varsel <- selectfun(model = thismodel, trafo = trafo, data = data,
    subset = subset, weights = weights, whichvar = svars, ctrl = thisctrl)

  if (terminal || inherits(varsel, "partysplit")) {
    thissplit <- if(terminal) NULL else varsel
    info <- nodeinfo <- thismodel[!(names(thismodel) %in% c("estfun"))] ## FIXME: needs updating
    info$nobs <- sw
    if (!ctrl$saveinfo) info <- NULL
  } else {
    ## criterion matrix and some sanity checking
    p <- varsel$criterion
    if(!is.matrix(p) || NCOL(p) != length(svars)) {
      stop("variable selection function does not provide a valid criterion matrix")
    }

    ## adjust p-values (for non-NA p-values), if any
    if("p.value" %in% rownames(p)) p["p.value", ] <- ctrl$padjust(p["p.value", ])
    if("log.p.value" %in% rownames(p)) p["log.p.value", ] <- ctrl$padjust(p["log.p.value", ], log = TRUE)

    ## determine criterion
    criterion <- ctrl$criterion
    critvalue <- ctrl$critvalue
    ## use log p-values if available (and p-values wanted)
    if(criterion == "p.value" && "log.p.value" %in% rownames(p)) {
      criterion <- "log.p.value"
      critvalue <- log(critvalue)
    }
    ## minimize p-values? (or maximize other criterion)
    minp <- criterion %in% c("p.value", "log.p.value")
    ## extract criterion values
    crit <- p[criterion, , drop = TRUE]

    ## no admissible splits in any variable, e.g., all
    ## split variables constant
    if (all(is.na(crit)) || (minp && all(crit >= critvalue)) || (!minp && all(crit <= critvalue))) {
      terminal <- TRUE
    }

    ## placeholder for missing criterion values
    crit[is.na(crit)] <- if(minp) Inf else -Inf

    if(!terminal) {
      ## optimal criterion, trying to break ties if any
      critopt <- if(minp) min(crit, na.rm = TRUE) else max(crit, na.rm = TRUE)
      ties <- which(abs(crit - critopt) < sqrt(.Machine$double.xmin))
      if(length(ties) > 1L) {
        ## break ties in p-values (if this is the "criterion") by using statistics (if available)
        if(minp && any(c("statistic", "log.statistic") %in% rownames(p))) {
          stat <- if("log.statistic" %in% rownames(p)) {
	    p["log.statistic", ties, drop = TRUE]
          } else {
	    p["statistic", ties, drop = TRUE]
          }
          ## subtract small value from (log-)p-value based on rank of (log-)statistic
          crit[ties] <- crit[ties] - rank(stat)/length(ties) * 1e-8      
        } else {
          ## FIXME: Can we do anything here? Always choose the first? Randomize?
        }
      }

      ## store up to "splittry" best variables
      jsel <- order(crit, decreasing = !minp)
      splittry <- pmin(sum(is.finite(crit)), ctrl$splittry)
      jsel <- jsel[1L:splittry]      
      jsel <- if(minp) jsel[crit[jsel] < critvalue] else jsel[crit[jsel] > critvalue]
    } else {
      jsel <- integer()
    }

    if(!terminal) { ## FIXME: if(!terminal || "criterion" %in% save)
      ## always include the criterion used as "criterion"
      p <- rbind(p, criterion = crit)

      ## switch to "statistic" and "p.value" omitting logs
      if(("log.statistic" %in% rownames(p)) && !("statistic" %in% rownames(p))) {
        p <- rbind(p, statistic = exp(p["statistic", ]))
      }
      if(("log.p.value" %in% rownames(p)) && !("p.value" %in% rownames(p))) {
        p <- rbind(p, p.value = exp(p["log.p.value", ]))
      }
      p <- p[-which(rownames(p) %in% c("log.statistic", "log.p.value")), , drop = FALSE]

      ## store optimal p-value (FIXME: Always p-value - or optimal criterion?)
      jopt <- if(minp) which.min(crit) else which.max(crit)
      if("p.value" %in% rownames(p)) {
	popt <- p["p.value", jopt]
      } else {
        popt <- NA_real_
      }
      names(popt) <- colnames(p)[jopt]
    }

    ## FIXME: Would we want this? I guess it's clearer/easier to keep these in.
    ## ## report on tests actually performed only
    ## p <- p[,!is.na(p["statistic",]) & is.finite(p["statistic",]),
    ##      drop = FALSE]
    info <- nodeinfo <- c(list(criterion = p, p.value = popt),  ## FIXME: update according to "save"
                varsel[!(names(varsel) %in% c("criterion", "converged"))],
                thismodel[!(names(thismodel) %in% c("estfun"))])
    info$nobs <- sw
    if (!ctrl$saveinfo) info <- NULL

    if (!is.null(varsel$splits)) {
      ### selectfun may return of a list of partysplit objects; use these for
      ### splitting; selectfun is responsible for making sure lookahead is implemented
      thissplit <- varsel$splits[[jsel[1L]]]
    } else {
      ### try to find an admissible split in data[, jsel]
      thissplit <- splitfun(model = thismodel, trafo = trafo, data = data, subset = subset, 
                  weights = weights, whichvar = jsel, ctrl = thisctrl)
    }
  }

  ### failed split search:
  if (is.null(thissplit))
    return(partynode(as.integer(id), info = info))

  ### successful split search: set-up node
  ret <- partynode(as.integer(id))
  ret$split <- thissplit
  ret$info <- info

  ### determine observations for splitting (only non-missings)
  snotNA <- subset[!subset %in% data[[varid_split(thissplit), type = "missings"]]]
  if (length(snotNA) == 0)
    return(partynode(as.integer(id), info = info))
  ### and split observations
  kidids <- kidids_node(ret, model.frame(data), obs = snotNA)

  ### compute probability of going left / right
  prob <- tabulate(kidids) / length(kidids) 
  # names(dimnames(prob)) <- NULL
  if (ctrl$majority)  ### go with majority
    prob <- as.double((1L:length(prob)) %in% which.max(prob))
  if (is.null(ret$split$prob))
    ret$split$prob <- prob

  ### compute surrogate splits
  if (ctrl$maxsurrogate > 0L) {
    pv <- partyvars
    pv[varid_split(thissplit)] <- 0
    pv <- which(pv > 0)
    if (ctrl$numsurrogate)
      pv <- pv[sapply(model.frame(data)[, pv], function(x) is.numeric(x) || is.ordered(x))]
    ret$surrogates <- .extree_surrogates(kidids, data = data, 
      weights = weights, subset = snotNA, 
      whichvar = pv,
      selectfun = svselectfun, splitfun = svsplitfun, ctrl = ctrl)
  }
  kidids <- kidids_node(ret, model.frame(data), obs = subset)

  ### proceed recursively
  kids <- vector(mode = "list", length = max(kidids)) 
  nextid <- id + 1L
  for (k in 1L:max(kidids)) {
    nextsubset <- subset[kidids == k]
    assign("depth", depth + 1L, envir = extree_env)
    kids[[k]] <- .extree_node(id = nextid, data = data, 
      trafo = trafo,
      selectfun = selectfun, splitfun = splitfun,
      svselectfun = svselectfun, svsplitfun = svsplitfun, 
      partyvars = partyvars, 
      weights = weights, subset = nextsubset, 
      ctrl = ctrl, info = nodeinfo, extree_env = extree_env)
    ### was: nextid <- max(nodeids(kids[[k]])) + 1L
    nextid <- get("maxid", envir = extree_env) + 1L
  }
  ret$kids <- kids

  return(ret)
}

### unbiased recursive partitioning: surrogate splits
.extree_surrogates <- function
(
    split, 			### integer vector with primary kidids
    data, 			### full data, readonly
    weights,
    subset, 			### subset of 1:nrow(data) with
				### non-missings in primary split
    whichvar, 			### partytioning variables
    selectfun, 			### variable selection and split
				### function
    splitfun,
    ctrl			### ctree_control()
) {

    if (length(whichvar) == 0) return(NULL)
    ms <- max(split)
    if (ms != 2) return(NULL) ### ie no multiway splits!
    dm <- matrix(0, nrow = nrow(model.frame(data)), ncol = ms)
    dm[cbind(subset, split)] <- 1
    thismodel <- list(estfun = dm)
    varsel <- selectfun(model = thismodel, trafo = NULL, data = data, subset = subset, 
                    weights = weights, whichvar = whichvar, ctrl = ctrl)
    p <- varsel$criterion
    ### partykit always used p-values, so expect some differences
    crit <- p[ctrl$criterion,,drop = TRUE]
    ### crit is maximised, but there might be ties
    ties <- which(abs(crit - max(crit, na.rm = TRUE)) < .Machine$double.eps)
    if (length(ties) > 1) {
        ### add a small value (< 1/1000) to crit derived from order of 
        ### teststat
        crit[ties] <- crit[ties] + 
            order(p["statistic", ties]) / (sum(ties) * 1000)
    }

    ret <- vector(mode = "list", length = min(c(length(whichvar), 
                                                ctrl$maxsurrogate)))

    for (i in 1L:length(ret)) {
        jsel <- which.max(crit)
        thisctrl <- ctrl
        thisctrl$minbucket <- 0L
        sp <- splitfun(model = thismodel, trafo = NULL, data = data, subset = subset, 
                       weights = weights, whichvar = jsel, ctrl = ctrl)
        if (is.null(sp)) next
        ret[[i]] <- sp
        tmp <- kidids_split(ret[[i]], model.frame(data), obs = subset)

        ### <FIXME> this needs fixing for multiway "split"
        tab <- table(tmp, split)
        if (tab[1, 1] < tab[1, 2]) {
            indx <- ret[[i]]$index
            ret[[i]]$index[indx == 1] <- 2L
            ret[[i]]$index[indx == 2] <- 1L
        }
        ### </FIXME>
        crit[which.max(crit)] <- -Inf
    }
    ret <- ret[!sapply(ret, is.null)]
    if (length(ret) == 0L) ret <- NULL
    return(ret)
}

extree_fit <- function(data, trafo, converged, selectfun = ctrl$selectfun, 
                       splitfun = ctrl$splitfun, svselectfun = ctrl$svselectfun, 
                       svsplitfun = ctrl$svsplitfun, partyvars, subset, weights, ctrl, doFit = TRUE) {
    ret <- list()

    ### <FIXME> use data$vars$z as default for partyvars </FIXME>
    ### <FIXME> try to avoid doFit </FIXME>

    nf <- names(formals(trafo))
    if (all(c("subset", "weights", "info", "estfun", "object") %in% nf)) {
        mytrafo <- trafo
    } else {
        stopifnot(all(c("y", "x", "offset", "weights", "start") %in% nf))
        stopifnot(!is.null(yx <- data$yx))
        mytrafo <- function(subset, weights, info, estfun = FALSE, object = FALSE, ...) {
            iy <- data[["yx", type = "index"]]
            if (is.null(iy)) {
                NAyx <- data[["yx", type = "missing"]]
                y <- yx$y
                x <- yx$x
                offset <- attr(yx$x, "offset")
                ### <FIXME> other ways of handling NAs necessary? </FIXME>
                subset <- subset[!(subset %in% NAyx)]
                if (NCOL(y) > 1) {
                    y <- y[subset,,drop = FALSE]
                } else {
                    y <- y[subset]
                }
                if (!is.null(x)) {
                    ax <- attributes(x)
                    ax$dim <- NULL
                    ax$dimnames <- NULL
                    x <- x[subset,,drop = FALSE]
                    for (a in names(ax)) attr(x, a) <- ax[[a]] ### terms, formula, ... for predict
                }
                w <- weights[subset]
                offset <- offset[subset]
                cluster <- data[["(cluster)"]][subset]
                if (all(c("estfun", "object") %in% nf)) { 
                    m <- trafo(y = y, x = x, offset = offset, weights = w, start = info$coef, 
                               cluster = cluster, estfun = estfun, object = object, ...)
                } else {
                    obj <- trafo(y = y, x = x, offset = offset, weights = w, start = info$coef, 
                                 cluster = cluster, ...)
                    m <- list(coefficients = coef(obj),
                              objfun = -as.numeric(logLik(obj)),
                              estfun = NULL, object = NULL)
                    if (estfun) m$estfun <- sandwich::estfun(obj)
                    if (object) m$object <- obj
                }
                if (!is.null(ef <- m$estfun)) {
                    ### ctree expects unweighted scores
                    if (!isTRUE(m$unweighted) && is.null(selectfun) && ctrl$testflavour == "ctree") 
                        m$estfun <- m$estfun / w
                    Y <- matrix(0, nrow = nrow(model.frame(data)), ncol = ncol(ef))
                    Y[subset,] <- m$estfun
                    m$estfun <- Y
                }
            } else {
                w <- libcoin::ctabs(ix = iy, subset = subset, weights = weights)[-1]
                offset <- attr(yx$x, "offset")
                cluster <- model.frame(data, yxonly = TRUE)[["(cluster)"]]
                if (all(c("estfun", "object") %in% nf)) { 
                    m <- trafo(y = yx$y, x = yx$x, offset = offset, weights = w, start = info$coef, 
                               cluster = cluster,
                               estfun = estfun, object = object, ...)
                } else {
                    obj <- trafo(y = yx$y, x = yx$x, offset = offset, weights = w, start = info$coef, 
                                 cluster = cluster, ...)
                    m <- list(coefficients = coef(obj),
                              objfun = -as.numeric(logLik(obj)),
                              estfun = NULL, object = NULL)
                    if (estfun) m$estfun <- sandwich::estfun(obj)
                    if (object) m$object <- obj
                    if (!is.null(obj$unweighted)) 
                        m$unweighted <- obj$unweighted
                    m$converged <- obj$converged ### may or may not exist
                }
                ### <FIXME> unweight scores in ctree or weight scores in
                ### mfluc (means: for each variable again) </FIXME>
                ### ctree expects unweighted scores
                if (!is.null(m$estfun))  {
                    if (!isTRUE(m$unweighted) && is.null(selectfun) && ctrl$testflavour == "ctree") 
                        m$estfun <- m$estfun / w
                }
                if (!is.null(ef <- m$estfun))
                    m$estfun <- rbind(0, ef)
            }
            return(m)
        }
    }
                 
    if (!ctrl$update) {
        rootestfun <- mytrafo(subset = subset, weights = weights)
        updatetrafo <- function(subset, weights, info, ...)
            return(rootestfun)
    } else {
        updatetrafo <- function(subset, weights, info, ...) {
            ret <- mytrafo(subset = subset, weights = weights, info = info, ...)
            if (is.null(ret$converged)) ret$converged <- TRUE
            conv <- TRUE
            if (is.function(converged)) conv <- converged(subset, weights)
            ret$converged <- ret$converged && conv
            if (!ret$converged) return(NULL)
            ret
        }
    }

    nm <- c("model", "trafo", "data", "subset", "weights", "whichvar", "ctrl")
    stopifnot(all(nm == names(formals(selectfun))))
    stopifnot(all(nm == names(formals(splitfun))))
    stopifnot(all(nm == names(formals(svselectfun))))
    stopifnot(all(nm == names(formals(svsplitfun))))

    if (!doFit) return(mytrafo)

    list(nodes = .extree_node(id = 1, data = data, trafo = updatetrafo, selectfun = selectfun, 
                              splitfun = splitfun, svselectfun = svselectfun, svsplitfun = svsplitfun, 
                              partyvars = partyvars, weights = weights, subset = subset, ctrl = ctrl),
         trafo = mytrafo)
}



### control arguments needed in this file
extree_control <- function(
  ## unbiased inference-based variable selection  
  criterion = "p.value",
  critvalue = 0.05,
  padjust = TRUE,

  ## technical
  applyfun = NULL, 
  cores = NULL,

  minsplit = 20L,
  minbucket = 7L, 
  minprob = 0.01, 
  nmax = Inf,
  stump = FALSE,
  lookahead = FALSE, ### try trafo() for daugther nodes before implementing the split
  maxsurrogate = 0L, 
  numsurrogate = FALSE,
  mtry = Inf,
  maxdepth = Inf, 
  multiway = FALSE, 
  splittry = 2L,
  majority = FALSE, 
  caseweights = TRUE, 
  saveinfo = TRUE,
  update = NULL,
  
  selectfun, ## FIXME: add default (ctree?)
  splitfun, ## FIXME: add default
  svselectfun, ## FIXME: add default
  svsplitfun, ## FIXME: add default

  ## legacy
  bonferroni = FALSE
) {

  ## p-value adjustment
  if(is.logical(padjust)) {
    padjust <- if(padjust) "sidak" else "none"
  }
  if(is.character(padjust)) {
    padjust <- switch(tolower(padjust),
      "bonferroni" = p_adjust_bonferroni,
      "sidak" = p_adjust_sidak,
      "none" = function(p, ...) p,
      stop(paste("unknown padjust method:", padjust))
    )
  }
  if(!is.function(padjust)) stop("padjust must be a logical, character, or function")

  ## optionally employ parallel computing facilities
  if(is.null(applyfun)) {
    if(is.null(cores)) cores <- 1L
    applyfun <- if(cores == 1L) {
      lapply
    } else if(.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(cores) 
      on.exit(parallel::stopCluster(cl))
      function(X, FUN, ...) parallel::parLapply(cl, X, FUN, ...)
    } else {
      function(X, FUN, ...) parallel::mclapply(X, FUN, ..., mc.cores = cores)    
    }
  }

  ### well, it is implemented but not correctly so (FIXME)
  if (multiway & maxsurrogate > 0L)
      stop("surrogate splits currently not implemented for multiway splits")
  
  ## var_select preprocessing
  if(is.list(selectfun) || is.character(selectfun) || "j" %in% names(formals(selectfun))) {
      
      var_sel <- .preprocess_select(selectfun, select_type = "var")
      
      selectfun <- function(model, trafo, data, subset, weights,
  	  whichvar, ctrl) {
  	  var_select_loop(model, trafo, data, subset, weights, whichvar, ctrl,
  	      var_select = var_sel)
      }
  }
  
  ## split_select preprocessing
  if(is.list(splitfun) || is.character(splitfun) || "j" %in% names(formals(splitfun))) {
      
      split_sel <- .preprocess_select(splitfun, select_type = "split")
      
      splitfun <- function(model, trafo, data, subset, weights,
  	  whichvar, ctrl) {
  	  split_select_loop(model = model, trafo = trafo, data = data,
  	      subset = subset, weights = weights, whichvar = whichvar,
  	      control = ctrl, split_select = split_sel)
      }
  }

  list(
    criterion = criterion,
    critvalue = critvalue,
    padjust = padjust,
    
    applyfun = applyfun,
    
       minsplit = minsplit, minbucket = minbucket, 
       minprob = minprob, stump = stump, nmax = nmax,
       lookahead = lookahead, mtry = mtry,
       maxdepth = maxdepth, multiway = multiway, splittry = splittry,
       maxsurrogate = maxsurrogate, 
       numsurrogate = numsurrogate, majority = majority,
       caseweights = caseweights, 
       saveinfo = saveinfo, update = update,
       selectfun = selectfun, splitfun = splitfun, svselectfun =
       svselectfun, svsplitfun = svsplitfun
  )
}


.objfun_test <- function(model, trafo, data, subset, weights, j, SPLITONLY, ctrl)
{

  x <- data[[j]]
  NAs <- data[[j, type = "missing"]]
  if (all(subset %in% NAs)) { 
    if (SPLITONLY) return(NULL)
    return(list(statistic = NA, p.value = NA))
  }

  ix <- data[[j, type = "index"]]
  ux <- attr(ix, "levels")
  ixtab <- libcoin::ctabs(ix = ix, weights = weights, subset = subset)[-1]
  ORDERED <- is.ordered(x) || is.numeric(x)
  
  linfo <- rinfo <- model
  minlogLik <- nosplitll <- trafo(subset = subset, weights = weights, info = model, estfun = FALSE)$objfun
  sp <- NULL
  
  if (ORDERED) {
    ll <- ctrl$applyfun(which(ixtab > 0), function(u) {
      sleft <- subset[LEFT <- (ix[subset] <= u)]
      sright <- subset[!LEFT]
      if (length(weights) > 0 && ctrl$caseweights) {
        if (sum(weights[sleft]) < ctrl$minbucket ||
            sum(weights[sright]) < ctrl$minbucket)
          return(Inf);
      } else {
        if (length(sleft) < ctrl$minbucket || 
            length(sright) < ctrl$minbucket)
          return(Inf);
      }
      if (ctrl$restart) {
        linfo <- NULL
        rinfo <- NULL
      }
      linfo <- trafo(subset = sleft, weights = weights, info = linfo, estfun = FALSE)
      rinfo <- trafo(subset = sright, weights = weights, info = rinfo, estfun = FALSE)
      ll <- linfo$objfun + rinfo$objfun
      return(ll)
    })
    minlogLik <- min(unlist(ll))
    if(minlogLik < nosplitll)
      sp <- which(ixtab > 0)[which.min(unlist(ll))]
    
  } else {
    xsubs <- factor(x[subset])
    ## stop if only one level left
    if(nlevels(xsubs) < 2) {
      if (SPLITONLY) {
        return(NULL)
      } else {
        return(list(statistic = NA, p.value = NA))
      } 
    }
    splits <- .mob_grow_getlevels(xsubs)
    ll <- ctrl$applyfun(1:nrow(splits), function(u) {
      sleft <- subset[LEFT <- xsubs %in% levels(xsubs)[splits[u,]]]
      sright <- subset[!LEFT]
      if (length(weights) > 0 && ctrl$caseweights) {
        if (sum(weights[sleft]) < ctrl$minbucket ||
            sum(weights[sright]) < ctrl$minbucket)
          return(Inf);
      } else {
        if (length(sleft) < ctrl$minbucket || 
            length(sright) < ctrl$minbucket)
          return(Inf);
      }
      if (ctrl$restart) {
        linfo <- NULL
        rinfo <- NULL
      }
      linfo <- trafo(subset = sleft, weights = weights, info = linfo, estfun = FALSE)
      rinfo <- trafo(subset = sright, weights = weights, info = rinfo, estfun = FALSE)
      ll <- linfo$objfun + rinfo$objfun
      return(ll)
    })
    minlogLik <- min(unlist(ll))
    if(minlogLik < nosplitll) {
      sp <- splits[which.min(unlist(ll)),] + 1L
      levs <- levels(x)
      if(length(sp) != length(levs)) {
        sp <- sp[levs]
        names(sp) <- levs
      }
    }
  }
  
  if (!SPLITONLY){
    ## split only if logLik improves due to splitting
    minlogLik <- ifelse(minlogLik == nosplitll, NA, minlogLik)
    return(list(statistic = -minlogLik, p.value = NA)) ### .extree_node maximises
  }
  if (is.null(sp) || all(is.na(sp))) return(NULL)
  if (ORDERED) {
    ### interpolate split-points, see https://arxiv.org/abs/1611.04561
    if (!is.factor(x) & ctrl$intersplit & sp < length(ux)) {
      sp <- (ux[sp] + ux[sp + 1]) / 2 
    } else {
      sp <- ux[sp]  ### x <= sp vs. x > sp
    }
    if (is.factor(sp)) sp <- as.integer(sp)
    ret <- partysplit(as.integer(j), breaks = sp,
                      index = 1L:2L)
  } else {
    ret <- partysplit(as.integer(j),
                      index = as.integer(sp))
  }
  return(ret)
}

.start_subset <- function(data) {
    ret <- 1:NROW(model.frame(data))
    if (length(data$yxmissings) > 0)
        ret <- ret[!(ret %in% data$yxmissings)]
    ret
}

