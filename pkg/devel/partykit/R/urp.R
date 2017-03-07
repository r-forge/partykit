
.urp_fit_1d <- function
(
    data, 				### full data, readonly
    partyvars, 				### partytioning variables,
					### a subset of 1:ncol(data)
    cluster = integer(0), 		### a blocking factor w/o NA
    ctrl				### control
) {

    ### transform partytioning variables
    X <- vector(mode = "list", length = NCOL(data))
    X[partyvars] <- lapply(partyvars, function(j) {
        x <- data[[j]]
        if (is.numeric(x)) {
            ret <- x
        } else if (is.ordered(x)) {
            sc <- 1:nlevels(x)
            if (!is.null(attr(x, "scores")))
                sc <- attr(x, "scores")
            ret <- matrix(sc[x], ncol = 1)
        } else if (is.factor(x)) {
            ret <- matrix(0, nrow = nrow(data), ncol = nlevels(x))
            ret[cbind(1:length(x), unclass(x))] <- 1 ### model.matrix(~ x - 1)
            ret[is.na(x),] <- NA
        } else {
            stop("cannot handle class", class(x))
        }
        storage.mode(ret) <- "double"
        ret
    })

    ### this is the update function
    updatefun <- function(trafo, subset, weights) {

        ### compute statistics and (optionally) p-values
        ### for a subset of observations and variables
        ### y is used for node ids when computing surrogate splits
        selectfun <- function(y = NULL, trafo, subset = integer(0), 
                              weights = integer(0), whichvar, info = NULL) 
        {

            ret <- list(criteria = matrix(NA, nrow = 2L, ncol = ncol(data)))
            colnames(ret$criteria) <- names(data)
            rownames(ret$criteria) <- c("statistic", "p.value")

            if (is.null(y)) {
                ### nrow(Y) = nrow(data)!!!
                tr <- trafo(subset, info = info)
                if (!tr$converged) return(ret)
            } else {
                ### y is kidids in .csurr and nothing else
                stopifnot(length(y) == length(subset))
                Y <- matrix(0, nrow = NROW(data), ncol = max(y))
                Y[cbind(subset, y)] <- 1 ### model.matrix(~ as.factor(y) - 1)
                storage.mode(Y) <- "double"
                tr <- list(estfun = Y)
            }
            Y <- tr$estfun
            if (!is.null(tr$index)) {
                if (length(tr$index) != nrow(data))
                    stop("incorrect index")
                ### index == 0 means NA
                subset <- subset[tr$index[subset] > 0]
                Y <- Y[tr$index + 1L,,drop = FALSE]
            }

            for (j in whichvar) {
                tst <- switch(ctrl$testflavour, 
                    "ctree" = .ctree_test_split(x = data[[j]], bdr = NULL, j = j, ctrl = ctrl, 
                                                X = X, Y = Y, iy = NULL, subset = subset, 
                                                weights = weights, cluster = cluster,
                                                splitonly = FALSE, minbucket =
                                                ctrl$minbucket),
                    "exhaustive" = .objfun_test_split(trafo = trafo, info = info, x = data[[j]], 
                                                      bdr = NULL, j = j, ctrl = ctrl, 
                                                      subset = subset,  weights = weights, 
                                                      cluster = cluster, splitonly = FALSE, 
                                                      minbucket = ctrl$minbucket),
                    "mfluc" = .fluct_test_split(estfun = Y, z = data[[j]], weights = weights,
                                                cluster = cluster, control = ctrl),
                    stop(ctrl$testflavour, "not yet implemented")
                )
                ### <FIXME> minbucket is updated in .urp_node but only after testing... </FIXME>
                ret$criteria["statistic", j] <- tst$statistic
                ret$criteria["p.value", j] <- tst$p.value
            }
            if (ctrl$bonferroni)
                ret$criteria["p.value",] <- ret$criteria["p.value",] * length(whichvar)

            ret <- c(ret, tr[names(tr) != "estfun"])

            ### compute best fitting cutpoint (according to minbucket)
            ### for a subset of observations and variables
            ### y is used for node ids when computing surrogate splits
            ### splitfun as part of the return object of selectfun allows
            ### returning splits found during selections (ie, exhaustive 
            ### searchs)
            ### splitfun already knows Y by lexical scoping, no need to 
            ### compute it twice
            ret$splitfun <- function(whichvar, minbucket) 
            {
                for (j in whichvar) {
                    ret <- switch(ctrl$splitflavour, 
                        "ctree" = .ctree_test_split(x = data[[j]], bdr = NULL, j = j, ctrl = ctrl,
                                                    X = X, Y = Y, iy = NULL, subset = subset, 
                                                    weights = weights, cluster = cluster,
                                                    splitonly = TRUE, minbucket =
                                                    minbucket),
                        "exhaustive" = .objfun_test_split(trafo = trafo, info = info, x = data[[j]], 
                                                          bdr = NULL, j = j, ctrl = ctrl, 
                                                          subset = subset,  weights = weights, 
                                                          cluster = cluster, splitonly = TRUE, 
                                                          minbucket = ctrl$minbucket),
                        stop(ctrl$splitflavour, "not yet implemented")
                    )
                    ### check if trafo can be successfully applied to all daugther nodes 
                    ### (converged = TRUE)
                    if (ctrl$lookahead & !is.null(ret)) {
                        sp <- kidids_split(ret, data, obs = subset)
                        conv <- sapply(unique(sp), function(i)
                            trafo(subset[sp == i], info = info, estfun = FALSE)$converged)
                        if (!all(conv)) ret <- NULL
                    }
                    if (!is.null(ret)) break()
                }
                ret
            }
            ret
        }

        tree <- .urp_node(id = 1L, data = data, 
                          selectfun = function(...) 
                              selectfun(..., trafo = trafo),
                          partyvars = partyvars, weights = weights, 
                          subset = subset, ctrl = ctrl)

        return(tree)
    }
    
    return(updatefun)
}

### faster but approximate version
.urp_fit_2d <- function
(
    data,                               ### full data, readonly
    partyvars,                          ### partytioning variables,
                                        ### a subset of 1:ncol(data)
    cluster = integer(0),                 ### a blocking factor w/o NA
    ctrl                                ### ctree_control()
) {

    bdr <- inum::inum(data, nmax = ctrl$nmax)
    X <- vector(mode = "list", length = NCOL(data))
    names(X) <- colnames(data)
    X[partyvars] <- lapply(partyvars, function(j) {
        x <- attr(bdr[[j]], "levels")
        if (is.logical(x)) {
            X <- rbind(0, diag(2))
        } else if (is.numeric(x)) {
            X <- rbind(0, matrix(x, ncol = 1L))
        } else if (is.factor(x) && !is.ordered(x)) {
            X <- rbind(0, diag(nlevels(x)))
        } else if (is.ordered(x)) {
            sc <- attr(data[[j]], "scores")
            if (is.null(sc)) sc <- 1:nlevels(x)
            X <- rbind(0, matrix(sc, ncol = 1L))
        } else {
            stop("cannot handle predictors of class", " ", sQuote(class(x)))
        }
        storage.mode(X) <- "double"
        return(X)
    })

    return(function(trafo, subset, weights) {

        ### compute statistics and (optionally) p-values
        ### for a subset of observations and variables
        ### y is used for node ids when computing surrogate splits
        selectfun <- function(y = NULL, trafo, subset = integer(0), 
                              weights = integer(0), whichvar, info = NULL) 
        {
    
            ret <- list(criteria = matrix(NA, nrow = 2L, ncol = ncol(data)))
            colnames(ret$criteria) <- names(data)
            rownames(ret$criteria) <- c("statistic", "p.value")

            if (is.null(y)) {
                tr <- trafo(subset = subset, info = info)
                if (!tr$converged) return(ret)
            } else {
                ### y is kidids in .csurr and nothing else
                stopifnot(length(y) == length(subset))
                Y <- rbind(0, max(y))
                iy <- numeric(NROW(data))
                iy[subset] <- as.integer(y)
                tr <- list(estfun = Y, index = iy)
            }
            Y <- tr$estfun
            iy <- tr$index
            if (is.null(iy)) stop("trafo did not return index")

            for (j in whichvar) {
                tst <- switch(ctrl$testflavour, 
                    "ctree" = .ctree_test_split(x = data[[j]], bdr = bdr, j = j, ctrl = ctrl,
                                                X = X, Y = Y, iy = iy, subset = subset,
                                                weights = weights, cluster = cluster,
                                                splitonly = FALSE, minbucket =
                                                ctrl$minbucket),
                    stop(ctrl$testflavour, "not yet implemented")
                )
                ret$criteria["statistic", j] <- tst$statistic
                ret$criteria["p.value", j] <- tst$p.value
            }
            if (ctrl$bonferroni)
                ret$criteria["p.value",] <- ret$criteria["p.value",] * 
                    length(whichvar)

            ret <- c(ret, tr[!(names(tr) %in% c("estfun", "index"))])

            ### compute best fitting cutpoint (according to minbucket)
            ### for a subset of observations and variables
            ### y is used for node ids when computing surrogate splits
            ret$splitfun <- function(whichvar,  minbucket)
            {
                for (j in whichvar) {
                    ret <- switch(ctrl$splitflavour, 
                        "ctree" = .ctree_test_split(x = data[[j]], bdr = bdr, j = j, ctrl = ctrl,
                                                    X = X, Y = Y, iy = iy, subset = subset,
                                                    weights = weights, cluster = cluster,
                                                    splitonly = TRUE, minbucket =
                                                    minbucket),
                        stop(ctrl$splitflavour, "not yet implemented")
                    )
                    if (!is.null(ret)) break()
                }
                return(ret)
           }
           ret
        }

        tree <- .urp_node(id = 1L, data = data, 
                          selectfun = function(...) 
                              selectfun(..., trafo = trafo),
                          partyvars = partyvars, weights = weights,
                          subset = subset, ctrl = ctrl)
        return(tree)
    })
}

.urp_fit <- function
(
    data, 
    partyvars, 
    cluster, 
    ctrl
) {

    if (ctrl$nmax < Inf)
        return(.urp_fit_2d(data = data, partyvars = partyvars, 
                           cluster = cluster, ctrl = ctrl))
    return(.urp_fit_1d(data = data, partyvars = partyvars,
                       cluster = cluster, ctrl = ctrl))
}



### unbiased recursive partitioning: set up new node
.urp_node <- function
(
    id = 1L, 			### id of this node
    data, 			### full data, readonly
    selectfun, 			### variable selection
                                ### and split function
    svselectfun = selectfun,    ### same for surrogate splits
    partyvars, 			### partytioning variables
                                ### a subset of 1:ncol(data)
    weights = integer(0L),	### optional case weights
    subset, 			### subset of 1:nrow(data)
                                ### for identifying obs for this node
    ctrl, 			### .urp_control()
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

    svars <- partyvars
    if (ctrl$mtry < Inf) {
        mtry <- min(length(partyvars), ctrl$mtry)
        svars <- .resample(partyvars, mtry)
    } 

    ### compute test statistics and p-values
    ### for _unbiased_ variable selection
    sf <- selectfun(subset = subset, weights = weights, whichvar = svars, info = info)
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
    names(pmin) <- colnames(data)[which.max(crit)]

    info <- c(list(criterion = p, p.value = pmin), 
                   sf[!(names(sf) %in% c("criteria", "splitfun"))])

    ### nothing "significant"
    if (all(crit <= ctrl$logmincriterion))
        return(partynode(as.integer(id), info = info))

    ### update sample size constraints on possible splits
    mb <- ctrl$minbucket
    mp <- ctrl$minprob
    swp <- ceiling(sw * mp)
    if (mb < swp) mb <- as.integer(swp)
    ### at most ctrl$splittry variables with meaningful criterion
    st <- pmin(sum(is.finite(crit)), ctrl$splittry)
    jsel <- rev(order(crit))[1:st]
    jsel <- jsel[crit[jsel] > ctrl$logmincriterion]
    ### try to find an admissible split in data[, jsel]
    thissplit <- sf$splitfun(whichvar = jsel, minbucket = mb)

    ### failed split search:
    if (is.null(thissplit))
        return(partynode(as.integer(id), info = info))

    ### successful split search: set-up node
    ret <- partynode(as.integer(id))
    ret$split <- thissplit
    ret$info <- info

    ### determine observations for splitting (only non-missings)
    snotNA <- subset[!is.na(data[subset, varid_split(thissplit)])]
    ### and split observations
    kidids <- kidids_node(ret, data, obs = snotNA)

    ### compute probability of going left / right
    prob <- tabulate(kidids) / length(kidids) 
    # names(dimnames(prob)) <- NULL
    if (ctrl$majority)  ### go with majority
        prob <- as.double((1L:length(prob)) %in% which.max(prob))
    if (is.null(ret$split$prob))
        ret$split$prob <- prob

    ### compute surrogate splits
    if (ctrl$maxsurrogate > 0L) {
        pv <- svars[svars != varid_split(thissplit)]
        if (ctrl$numsurrogate)
            pv <- pv[sapply(data[, pv], is.numeric)]
        ret$surrogates <- .urp_surrogates(kidids, data = data, 
            weights = weights, subset = snotNA, 
            partyvars = pv,
            selectfun = svselectfun, ctrl = ctrl)
    }
    kidids <- kidids_node(ret, data, obs = subset)

    ### proceed recursively
    kids <- vector(mode = "list", length = max(kidids)) 
    nextid <- id + 1L
    for (k in 1L:max(kidids)) {
        nextsubset <- subset[kidids == k]
        assign("depth", depth + 1L, envir = cenv)
        kids[[k]] <- .urp_node(id = nextid, data = data, 
            selectfun = selectfun, partyvars = partyvars, 
            weights = weights, subset = nextsubset, 
            ctrl = ctrl, info = info, cenv = cenv)
        ### was: nextid <- max(nodeids(kids[[k]])) + 1L
        nextid <- get("maxid", envir = cenv) + 1L
    }
    ret$kids <- kids

    return(ret)
}

### unbiased recursive partitioning: surrogate splits
.urp_surrogates <- function
(
    split, 			### integer vector with primary kidids
    data, 			### full data, readonly
    weights,
    subset, 			### subset of 1:nrow(data) with
				### non-missings in primary split
    partyvars, 			### partytioning variables
    selectfun, 			### variable selection and split
				### function
    ctrl			### ctree_control()
) {

    sf <- selectfun(y = split, subset = subset, weights = weights, whichvar = partyvars)
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

    ret <- vector(mode = "list", length = min(c(length(partyvars), 
                                                ctrl$maxsurrogate)))

    for (i in 1L:length(ret)) {
        jsel <- which.max(crit)
        sp <- sf$splitfun(whichvar = jsel, minbucket = 0L)
        if (is.null(sp)) next
        ret[[i]] <- sp
        tmp <- kidids_split(ret[[i]], data, obs = subset)

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


### parse formula and grow unbiased tree in a generic way
.urp_tree <- function
(
    call, 		### match.call of user-visible function
    frame, 		### parent.frame of user-visible function
    data = NULL, 
    data_asis = FALSE, 
    control, 		### .urp_control() or more
    trafofun, 		### function for transformations
    doFit = TRUE	### grow tree or set-up only
) {

    ### call and frame come from user-visible functions, like ctree()
    ### strata for cforest only
    m <- match(c("formula", "data", "subset", "weights", "offset", "cluster", "strata"),
               names(call), 0L)
    mf <- call[c(1, m)]
    formula <- eval(mf$formula, frame)
    na.action <- eval(call$na.action, frame)

    f <- if (inherits(formula, "Formula")) formula else Formula(formula)
    ### formula must feature one lhs and one rhs
    if (length(length(f)) != 2L)
        stop("incorrect formula") 
    ### three-part formula allowed in rhs 
    if (!(length(f)[2] %in% 1L:2L))
        stop("incorrect formula")
    ### only simple formula allowed in lhs 
    if (length(f)[1] != 1)
        stop("incorrect formula")
    mf$formula <- f

    if (!is.null(data) && data_asis) {
        mf <- data
    } else {
        mf[[1]] <- quote(stats::get_all_vars)
        mf <- eval(mf, frame) 
    }
    ### Note: get terms from data and not mf as the latter
    ### might include variables "weights" or "subset"
    mfterms <- terms(f, data = data) 
    ### there might be dots in formula, fdot
    ### is formula with dots replaced
    fdot <- attr(mfterms, "Formula_without_dot")
    if (!is.null(fdot)) fdot <- Formula(formula(fdot, collapse = TRUE))
    av <- if (!is.null(fdot)) all.vars(fdot) else all.vars(f)
    nm <- names(mf) 
    nm[!(nm %in% av)] <- paste("(", nm[!(nm %in% av)], ")", sep = "")
    names(mf) <- nm 

    ### <FIXME> do we want this? </FIXME>
    if (!data_asis)
        mf <- na.action(mf)

    weights <- model.weights(mf) 
    if (is.null(weights)) weights <- integer(0)
    cluster <- mf[["(cluster)"]] 
    if (!is.null(cluster)) {
        if (!is.factor(cluster)) stop("cluster is not a factor")
    }
    subset <- mf[["(subset)"]] 
    if (!is.null(subset)) {
        subset <- which(subset)
    } else {
        subset <- 1L:nrow(mf)
    }

    if (length(f)[2] == 1L) { ### y ~ z or y ~ .
        if (is.null(fdot)) fdot <- f
        modelf <- formula(fdot, lhs = 1L, rhs = 0L)
        partf <- formula(fdot, lhs = 0L, rhs = 1L)
    } else if (length(f)[2] == 2L) { ### y ~ x | z
        if (!is.null(fdot))
            stop("dots are not allowed in multipart formulas")
        modelf <- formula(f, lhs = 1L, rhs = 1L)
        partf <- formula(f, lhs = 0L, rhs = 2L)
    } 
    zvars <- rownames(attr(terms(partf, data = mf), "factors"))

    ### returns a _function_ (trafo, subset, weights)
    ### for growing the tree, weights = integer(0) must work
    treefun <- .urp_fit(mf, partyvars = match(zvars, colnames(mf)), 
                        cluster = cluster, ctrl = control)
    if (!isTRUE(all.equal(names(formals(treefun)), 
                          c("trafo", "subset", "weights"))))
        stop(".urp_fit returned incorrect results")
    ### returns a _function_ (subset) for computing estfun (essentially)
    trafo <- trafofun(modelf, data = mf, ctrl = control)
    if (!isTRUE(all.equal(names(formals(trafo))[1L], "subset")))
        stop("trafofun return incorrect")
    
    ret <- list(treefun = treefun, trafo = trafo, mf = mf, terms = mfterms,
                partyvars = match(zvars, colnames(mf)))
    if (!doFit)
        return(ret)

    ### grow the tree
    tree <- treefun(trafo, subset = subset, weights)
    ret$nodes <- tree
    return(ret)
}

### control arguments needed in this file
.urp_control <- function
(
    criterion, 
    logmincriterion, 
    minsplit = 20L,
    minbucket = 7L, 
    minprob = 0.01, 
    nmax = Inf,
    stump = FALSE,
    lookahead = FALSE, ### try trafo() for daugther nodes before implementing the split
    MIA = FALSE,
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
    testflavour = c("ctree", "exhaustive", "mfluc"),
    bonferroni = FALSE,
    splitflavour = c("ctree", "exhaustive"),
    vcov = c("opg", "info", "sandwich"), ### mob_control
    ordinal = c("chisq", "max", "L2"),   ### mob_control
    nrep = 10000,                        ### mob_control
    trim = 0.1                           ### mob_control, TODO: minprob maxprob
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

    if (MIA && maxsurrogate > 0)
        warning("Mixing MIA splits with surrogate splits does not make sense")

    if (MIA && majority)
        warning("Mixing MIA splits with majority does not make sense")

    list(criterion = criterion, logmincriterion = logmincriterion,
         minsplit = minsplit, minbucket = minbucket, 
         minprob = minprob, stump = stump, nmax = nmax,
         lookahead = lookahead, mtry = mtry,
         maxdepth = maxdepth, multiway = multiway, splittry = splittry,
         MIA = MIA, maxsurrogate = maxsurrogate, 
         numsurrogate = numsurrogate, majority = majority,
         caseweights = caseweights, applyfun = applyfun,
         testflavour = match.arg(testflavour), 
         bonferroni = bonferroni,
         splitflavour = match.arg(splitflavour),
         vcov = match.arg(vcov),
         ordinal = match.arg(ordinal),
         nrep = nrep, trim = trim)
}
