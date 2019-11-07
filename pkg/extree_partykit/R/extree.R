extree <- function(data, 
    trafo,
    control = extree_control(#TODO: rename selectfun -> var_select, splitfun -> split_select, 
        ...), 
    converged = NULL,
    ...) {
    
    ## check / preprocess extree data
    subset <- partykit:::.start_subset(data = data)
    weights <- model.weights(model.frame(data))
    
    
    
    
    ## trafo preprocessing
    mytrafo <- function(subset, weights, info = NULL, estfun = TRUE, object = TRUE) {
        trafo(subset, data = data, weights, info = NULL, estfun = TRUE, object = TRUE)
    }
    
    ## TODO: converged preprocessing (if needed)
    
    ## set up trafo
    update <- function(subset, weights, control, doFit = TRUE) {
        partykit::extree_fit(data = data, trafo = mytrafo, converged = converged,
            partyvars = data$variables$z, subset = subset,
            weights = weights, ctrl = control, doFit = doFit)  
    }
    
    ## fit
    tree <- update(subset = subset, weights = weights, control = control)
    
    ## TODO: prepare extree object
    
}


## helper function for .preprocess_var_select
.get_strategy_function <- function(strategy, var_type = NULL, 
    select_type = "var") {
    
    if(is.null(var_type)) var_type <- ""
    if(length(var_type) != 1) stop("This function allows only one var_type")
    
    ## find matching objects and set up list
    ## FIXME: better solution than via search()?
    snam <- sprintf(paste0("^", select_type, "_select_%s_%s"), strategy, var_type)
    onam <- unlist(lapply(search(), objects))
    onam <- unique(onam[grep(snam, onam)])
    strategy <- lapply(onam, get)
    names(strategy) <- ifelse(var_type == "", gsub(snam, "", onam), var_type)
    
    ## drop non-functions
    strategy <- strategy[sapply(strategy, is.function)]
    
    return(strategy)
}




.preprocess_select <- function(select, select_type = "var") {
    
    ## function: return as is
    if(is.function(select)) return(select)
    
    ## character: return appropriate function
    if(is.character(select)) {
        return(
            .get_strategy_function(select, strategy = "", 
                select_type = select_type)
            )
    }
    
    ## list: go through all elements and return accordingly
    if(is.list(select)) {
        
        get_strategy <- function(select_nam) {
            
            ## return function if function
            if(is.function(select[[select_nam]])) {
                return(select[[select_nam]])
            } 
            
            ## get appropriat function if character
            if(is.character(select[[select_nam]])) {
                return(.get_strategy_function(select[[select_nam]], 
                    var_type = select_nam, select_type = select_type))
            } 
            
            ## if none of the above -> ERROR
            stop("select can only be functions or characters.")
            
        }
        
        ## go through all list elements and choose approriate function
        select_list <- lapply(names(select), get_strategy)
        return(select_list)
    }
}



# <FIXME> (HS) better name of function
.get_varclass <- function(select_list, data, j) {
    ### which class is variable?
    varclass <- class(extree_variable(x = data, i = j))
    
    ### Use most appropriate class (1st), if more than one is available
    ### Remove varclass if no var_select function is available
    varclass <- varclass[varclass %in% names(select_list)][1]
    
    ### if no function for this class is provided use default function
    if(length(varclass) == 0 | is.na(varclass)) {
        stopifnot("default" %in% names(select_list)) # <FIXME> (HS) improve error message
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

    ## set up return list + criteria matrix
    ret <- list(criteria = matrix(NA, nrow = 2L, ncol = ncol(model.frame(data))))
    rownames(ret$criteria) <- c("statistic", "p.value")
    colnames(ret$criteria) <- names(model.frame(data))
    if (length(whichvar) == 0) return(ret)

    ## loop over all relevant variables and use var_select function supplied
    for (j in whichvar) {
        
        tst <- selector(select = var_select, model = model, trafo = trafo, 
            data = data, subset = subset, weights = weights, j = j, 
            control = control)
        
        ret$criteria["statistic",j] <- tst$statistic
        ret$criteria["p.value",j] <- tst$p.value
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
        }
        
        ## stop if a split was found, otherwise continue with next possible var
        if (!is.null(ret)) break()
    }
    return(ret)
}

## Select function old
.select <- function(model, trafo, data, subset, weights, whichvar, ctrl, FUN) {
    ret <- list(criteria = matrix(NA, nrow = 2L, ncol = ncol(model.frame(data))))
    rownames(ret$criteria) <- c("statistic", "p.value")
    colnames(ret$criteria) <- names(model.frame(data))
    if (length(whichvar) == 0) return(ret)
    ### <FIXME> allow joint MC in the absense of missings; fix seeds
    ### write ctree_test / ... with whichvar and loop over variables there
    ### </FIXME>
    for (j in whichvar) {
        tst <- FUN(model = model, trafo = trafo, data = data, 
                   subset = subset, weights = weights, j = j, 
                   SPLITONLY = FALSE, ctrl = ctrl)
        ret$criteria["statistic",j] <- tst$statistic
        ret$criteria["p.value",j] <- tst$p.value
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
    cenv = NULL			### environment for depth and maxid
) {

    ### depth keeps track of the depth of the tree
    ### which has to be < than maxdepth
    ### maxit is the largest id in the left subtree
    if (is.null(cenv)) {
        cenv <- new.env()
        assign("depth", 0L, envir = cenv)
    }
    depth <- get("depth", envir = cenv)
    assign("maxid", id, envir = cenv)
    if (depth >= ctrl$maxdepth)
        return(partynode(as.integer(id)))

    ### check for stumps
    if (id > 1L && ctrl$stump) 
        return(partynode(as.integer(id)))

    ### sw is basically the number of observations
    ### which has to be > minsplit in order to consider
    ### the node for splitting
    if (length(weights) > 0L) {
        if (ctrl$caseweights) {
            sw <- sum(weights[subset]) 
        } else {
            sw <- sum(weights[subset] > 0L)
        }
    } else {
        sw <- length(subset)
    }
    if (sw < ctrl$minsplit) 
        return(partynode(as.integer(id)))

    svars <- which(partyvars > 0)
    if (ctrl$mtry < Inf) {
        mtry <- min(sum(partyvars > 0), ctrl$mtry)
        svars <- .resample(svars, mtry, prob = partyvars[partyvars > 0])
    } 

    thismodel <- trafo(subset = subset, weights = weights, info = info, 
                       estfun = TRUE, object = TRUE)
    if (is.null(thismodel))
        return(partynode(as.integer(id)))

    ### update sample size constraints on possible splits
    ### need to do this here because selectfun might consider splits
    mb <- ctrl$minbucket
    mp <- ctrl$minprob
    swp <- ceiling(sw * mp)
    if (mb < swp) mb <- as.integer(swp)
    thisctrl <- ctrl
    thisctrl$minbucket <- mb

    ### compute test statistics and p-values
    ### for _unbiased_ variable selection
    sf <- selectfun(model = thismodel, trafo = trafo, data = data, subset = subset, weights = weights, 
                    whichvar = svars, ctrl = thisctrl)

    if (inherits(sf, "partysplit")) {
        thissplit <- sf
        info <- nodeinfo <- thismodel[!(names(thismodel) %in% c("estfun"))]
        info$nobs <- sw
        if (!ctrl$saveinfo) info <- NULL
    } else {
        if (ctrl$bonferroni) 
            ### make sure to correct for _non-constant_ variables only
            sf$criteria["p.value",] <- sf$criteria["p.value",] * 
                                       sum(!is.na(sf$criteria["p.value",]))
        ### selectfun might return other things later to be used for info
        p <- sf$criteria

        crit <- p[ctrl$criterion,,drop = TRUE]
        if (all(is.na(crit))) 
            return(partynode(as.integer(id)))

        crit[is.na(crit)] <- -Inf
        ### crit is maximised, but there might be ties
        ties <- which(abs(crit - max(crit, na.rm = TRUE)) < sqrt(.Machine$double.xmin))
        if (length(ties) > 1) {
            ### add a small value (< 1/1000) to crit derived from rank of 
            ### teststat
            crit[ties] <- crit[ties] + 
                rank(p["statistic", ties]) / (sum(ties) * 1000)
        }

        ### switch from log(1 - pval) to pval for info slots
        ### switch from log(statistic) to statistic
        ### criterion stays on log scale to replicate variable selection
        p <- rbind(p, criterion = crit)
        p["statistic",] <- exp(p["statistic",])
        p["p.value",] <- -expm1(p["p.value",])
        pmin <- p["p.value", which.max(crit)]
        names(pmin) <- colnames(model.frame(data))[which.max(crit)]

        ### report on tests actually performed only
        p <- p[,!is.na(p["statistic",]) & is.finite(p["statistic",]),
               drop = FALSE]
        info <- nodeinfo <- c(list(criterion = p, p.value = pmin), 
                              sf[!(names(sf) %in% c("criteria", "converged"))],
                              thismodel[!(names(thismodel) %in% c("estfun"))])
        info$nobs <- sw
        if (!ctrl$saveinfo) info <- NULL

        ### nothing "significant"
        if (all(crit <= ctrl$logmincriterion))
            return(partynode(as.integer(id), info = info))

        ### at most ctrl$splittry variables with meaningful criterion
        st <- pmin(sum(is.finite(crit)), ctrl$splittry)
        jsel <- rev(order(crit))[1:st]
        jsel <- jsel[crit[jsel] > ctrl$logmincriterion]
        if (!is.null(sf$splits)) {
            ### selectfun may return of a list of partysplit objects; use these for
            ### splitting; selectfun is responsible for making sure lookahead is implemented
            thissplit <- sf$splits[[jsel[1]]]
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
        assign("depth", depth + 1L, envir = cenv)
        kids[[k]] <- .extree_node(id = nextid, data = data, 
            trafo = trafo,
            selectfun = selectfun, splitfun = splitfun,
            svselectfun = svselectfun, svsplitfun = svsplitfun, 
            partyvars = partyvars, 
            weights = weights, subset = nextsubset, 
            ctrl = ctrl, info = nodeinfo, cenv = cenv)
        ### was: nextid <- max(nodeids(kids[[k]])) + 1L
        nextid <- get("maxid", envir = cenv) + 1L
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
    sf <- selectfun(model = thismodel, trafo = NULL, data = data, subset = subset, 
                    weights = weights, whichvar = whichvar, ctrl = ctrl)
    p <- sf$criteria
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
extree_control <- function
(
    criterion, 
    logmincriterion, 
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
    applyfun = NULL, 
    cores = NULL,
    saveinfo = TRUE,
    bonferroni = FALSE,
    update = NULL,
    selectfun, 
    splitfun, 
    svselectfun, 
    svsplitfun
) {

    ## apply infrastructure for determining split points
    if (is.null(applyfun)) {
        applyfun <- if(is.null(cores)) {
            lapply
        } else {
            function(X, FUN, ...)
                parallel::mclapply(X, FUN, ..., mc.cores = cores)
        }
    }

    ### well, it is implemented but not correctly so
    if (multiway & maxsurrogate > 0L)
        stop("surrogate splits currently not implemented for multiway splits")
    
    ## var_select preprocessing
    if("j" %in% names(formals(selectfun)) | is.list(selectfun)) {
        
        selectfun <- function(model, trafo, data, subset, weights, 
            whichvar, ctrl) {
            var_sel <- .preprocess_select(selectfun, select_type = "var")
            var_select_loop(model, trafo, data, subset, weights, whichvar, ctrl, 
                var_select = var_sel)
        }
    }
    
    ## split_select preprocessing
    if("j" %in% names(formals(splitfun)) | is.list(splitfun)) {
        
        splitfun <- function(model, trafo, data, subset, weights,
            whichvar, ctrl) {
            split_sel <- .preprocess_select(splitfun, select_type = "split")
            split_select_loop(model = model, trafo = trafo, data = data, 
                subset = subset, weights = weights, whichvar = whichvar, 
                control = ctrl, split_select = split_sel)
        }
    }

    list(criterion = criterion, logmincriterion = logmincriterion,
         minsplit = minsplit, minbucket = minbucket, 
         minprob = minprob, stump = stump, nmax = nmax,
         lookahead = lookahead, mtry = mtry,
         maxdepth = maxdepth, multiway = multiway, splittry = splittry,
         maxsurrogate = maxsurrogate, 
         numsurrogate = numsurrogate, majority = majority,
         caseweights = caseweights, applyfun = applyfun,
         saveinfo = saveinfo, bonferroni = bonferroni, update = update,
         selectfun = selectfun, splitfun = splitfun, svselectfun =
         svselectfun, svsplitfun = svsplitfun)
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
