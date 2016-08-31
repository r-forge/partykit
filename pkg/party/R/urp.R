
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
        ### length(partyvars) == 1 will lead to sample.int instead of sample; 
        ### see example(sample)
        resample <- function(x, ...) x[sample.int(length(x), ...)]
        svars <- resample(partyvars, mtry)
    } 

    ### compute test statistics and p-values
    ### for _unbiased_ variable selection
    sf <- selectfun(subset = subset, whichvar = svars)
    ### selectfun might return other things later to be used for info
    p <- sf$criteria
    crit <- p[ctrl$criterion,,drop = TRUE]
    if (all(is.na(crit))) 
        return(partynode(as.integer(id)))

    crit[is.na(crit)] <- -Inf
    ### crit is maximised, but there might be ties
    ties <- which(abs(crit - max(crit, na.rm = TRUE)) < .Machine$double.eps)
    if (length(ties) > 1) {
        ### add a small value (< 1/1000) to crit derived from order of 
        ### teststat
        crit[ties] <- crit[ties] + 
            order(p["statistic", ties]) / (sum(ties) * 1000)
    }

    ### switch from log(1 - pval) to pval for info slots
    p <- rbind(p, criterion = crit)
    p["p.value",] <- -expm1(p["p.value",])
    pmin <- p["p.value", which.max(crit)]
    names(pmin) <- colnames(data)[which.max(crit)]

    info <- c(list(criterion = p, p.value = pmin), 
                   sf[!(names(sf) %in% c("criteria", "splitfun"))])

    ### nothing "significant"
    if (all(crit < ctrl$logmincriterion))
        return(partynode(as.integer(id), info = info))

    ### update sample size constraints on possible splits
    mb <- ctrl$minbucket
    mp <- ctrl$minprob
    swp <- ceiling(sw * mp)
    if (mb < swp) mb <- as.integer(swp)
    jsel <- rev(order(crit))[1:ctrl$splittry]
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
    ret$split$prob <- prob

    ### compute surrogate splits
    if (ctrl$maxsurrogate > 0L)
        ret$surrogates <- .urp_surrogates(kidids, data = data, 
            subset = snotNA, 
            partyvars = svars[svars != varid_split(thissplit)],
            selectfun = svselectfun, ctrl = ctrl)
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
            ctrl = ctrl, cenv = cenv)
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
    subset, 			### subset of 1:nrow(data) with
				### non-missings in primary split
    partyvars, 			### partytioning variables
    selectfun, 			### variable selection and split
				### function
    ctrl			### ctree_control()
) {

    sf <- selectfun(y = split, subset = subset, whichvar = partyvars)
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
    growfun, 		### function for growing trees
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
    mfterms <- terms(f, data = mf) 
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
    treefun <- growfun(mf, partyvars = match(zvars, colnames(mf)), 
                       cluster = cluster, ctrl = control)
    if (!isTRUE(all.equal(names(formals(treefun)), 
                          c("trafo", "subset", "weights"))))
        stop("growfun return incorrect")
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
    stump = FALSE,
    maxsurrogate = 0L, 
    mtry = Inf,
    maxdepth = Inf, 
    multiway = FALSE, 
    splittry = 2L,
    majority = FALSE, 
    caseweights = TRUE, 
    applyfun = NULL, 
    cores = NULL
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

    list(criterion = criterion, logmincriterion = logmincriterion,
         minsplit = minsplit, minbucket = minbucket, 
         minprob = minprob, stump = stump, mtry = mtry,
         maxdepth = maxdepth, multiway = multiway, splittry = splittry,
         maxsurrogate = maxsurrogate, majority = majority,
         caseweights = caseweights, applyfun = applyfun)
}
