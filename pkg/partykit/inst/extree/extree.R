
### unbiased recursive partitioning: set up new node
.extree_node <- function
(
    id = 1L, 			### id of this node
    data, 			### full data, readonly
    updatefun,
    selectfun, 			### variable selection
    splitfun,                   ### split selection
    svselectfun = selectfun,    ### same for surrogate splits
    svsplitfun = splitfun,      ### same for surrogate splits
    partyvars, 			### partytioning variables
                                ### a subset of 1:ncol(model.frame(data))
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

    thismodel <- updatefun(subset = subset, weights = weights, info = info)
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
    sf <- selectfun(thismodel, subset = subset, weights = weights, 
                    whichvar = svars, ctrl = thisctrl)

    if (inherits(sf, "partysplit")) {
        thissplit <- sf
        info <- nodeinfo <- thismodel[!(names(thismodel) %in% c("estfun"))]
        if (!ctrl$saveinfo) info <- NULL
    } else {
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
                              sf[!(names(sf) %in% c("criterion"))],
                              thismodel[!(names(thismodel) %in% c("estfun"))])
        if (!ctrl$saveinfo) info <- NULL

        ### nothing "significant"
        if (all(crit <= ctrl$logmincriterion))
            return(partynode(as.integer(id), info = info))

        ### at most ctrl$splittry variables with meaningful criterion
        st <- pmin(sum(is.finite(crit)), ctrl$splittry)
        jsel <- rev(order(crit))[1:st]
        jsel <- jsel[crit[jsel] > ctrl$logmincriterion]
        ### try to find an admissible split in data[, jsel]
        thissplit <- splitfun(thismodel, subset = subset, weights = weights, 
                              whichvar = jsel, ctrl = thisctrl)
    }

    ### failed split search:
    if (is.null(thissplit))
        return(partynode(as.integer(id), info = info))

    ### successful split search: set-up node
    ret <- partynode(as.integer(id))
    ret$split <- thissplit
    ret$info <- info

    ### determine observations for splitting (only non-missings)
    snotNA <- subset[!subset %in% missings(data, varid_split(thissplit))]
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
            pv <- pv[sapply(model.frame(data)[, pv], function(x) is.numeric(x) || is.ordered(x))]
        ret$surrogates <- .extree_surrogates(kidids, data = data, 
            weights = weights, subset = snotNA, 
            partyvars = pv,
            selectfun = svselectfun, splitfun = svsplitfun, ctrl = ctrl)
    }
    kidids <- kidids_node(ret, model.frame(data), obs = subset)

    ### proceed recursively
    kids <- vector(mode = "list", length = max(kidids)) 
    nextid <- id + 1L
    for (k in 1L:max(kidids)) {
        nextsubset <- subset[kidids == k]
        assign("depth", depth + 1L, envir = cenv)
        kids[[k]] <- .extree_node(id = nextid, data = data, updatefun =
        updatefun,
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
    partyvars, 			### partytioning variables
    selectfun, 			### variable selection and split
				### function
    splitfun,
    ctrl			### ctree_control()
) {

    dm <- matrix(0, nrow = nrow(model.frame(data)), ncol = max(split))
    dm[cbind(subset, split)] <- 1
    thismodel <- list(estfun = dm)
    sf <- selectfun(thismodel, subset = subset, weights = weights, whichvar = partyvars,
    ctrl = ctrl)
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
        thisctrl <- ctrl
        thisctrl$minbucket <- 0L
        sp <- splitfun(thismodel, subset = subset, weights = weights, whichvar = jsel,
        ctrl = ctrl)
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

.extree_fit <- function(data, trafo, partyvars, subset, weights, ctrl) {
    ret <- list()
    if (!ctrl$update) {
        rootestfun <- trafo(subset = subset, weights = weights)
        updatefun <- function(subset, weights, info)
            return(rootestfun)
    } else {
        updatefun <- function(subset, weights, info) {
            ret <- trafo(subset = subset, weights = weights, info = info)
            if (!ret$converged) return(NULL)
            ret
        }
    }

    selectfun <- function(model, subset, weights, whichvar, ctrl) {
        ret <- list(criteria = matrix(NA, nrow = 2L, ncol = ncol(model.frame(data))))
        rownames(ret$criteria) <- c("statistic", "p.value")
        colnames(ret$criteria) <- names(model.frame(data))
        ### seed <- .Random.seed
        for (j in whichvar) {
            ### set.seed(seed)
            tst <- switch(ctrl$testflavour,
                "ctree" = .ctree_test(model = model, trafo = trafo, data = data, 
                                      subset = subset, weights = weights, j = j, 
                                      SPLITONLY = FALSE, ctrl = ctrl),
                "exhaustive" = .objfun_test(model = model, trafo = trafo, data = data, 
                                            subset = subset, weights = weights, j = j, 
                                            SPLITONLY = FALSE, ctrl = ctrl),
                "mfluc" = .mfluc_test(model = model, trafo = trafo, data = data, 
                                      subset = subset, weights = weights, j = j, 
                                      SPLITONLY = FALSE, ctrl = ctrl),
                stop("not implemented")
            )
            ret$criteria["statistic",j] <- tst$statistic
            ret$criteria["p.value",j] <- tst$p.value
        }
        if (ctrl$bonferroni) 
            ### make sure to correct for _non-constant_ variables only
            ret$criteria["p.value",] <- ret$criteria["p.value",] * 
                                        sum(!is.na(ret$criteria["p.value",]))
        ret
    }
    
    splitfun <- function(model, subset, weights, whichvar, ctrl) {
        for (j in whichvar) {
            x <- model.frame(data)[[j]]
            if (ctrl$multiway && is.factor(x) && !is.ordered(x) &&
                (ctrl$maxsurrogate == 0) && nlevels(x[subset, drop = TRUE]) > 1) 
            {
                index <- 1L:nlevels(x)
                if (length(weights) > 0) {
                    xt <- xtabs(weights ~ x, subset = subset)
                } else {
                    xt <- xtabs(~ x, subset = subset)
                }
                index[xt == 0] <- NA
                ### maybe multiway is not so smart here as
                ### nodes with nobs < minbucket could result
                index[xt > 0 & xt < minbucket] <- nlevels(x) + 1L
                if (length(unique(index)) == 1) {
                    ret <- NULL
                } else {
                    index <- unclass(factor(index))
                    ret <- partysplit(as.integer(j), index = as.integer(index))
                }
            } else {
                ret <- switch(ctrl$splitflavour, 
               "ctree" = .ctree_test(model = model, trafo = trafo, data = data, 
                                      subset = subset, weights = weights, j = j, 
                                      SPLITONLY = TRUE, ctrl = ctrl),
                "exhaustive" = .objfun_test(model = model, trafo = trafo, data = data, 
                                            subset = subset, weights = weights, j = j, 
                                            SPLITONLY = TRUE, ctrl = ctrl),
                stop("not implemented"))
            }
            ### check if trafo can be successfully applied to all daugther nodes 
            ### (converged = TRUE)
            if (ctrl$lookahead & !is.null(ret)) {
                sp <- kidids_split(ret, model.frame(data), obs = subset)
                conv <- sapply(unique(sp), function(i)
                        trafo(subset[sp == i], weights = weights)$converged)
                if (!all(conv)) ret <- NULL
            }
            if (!is.null(ret)) break()
        }
        return(ret)
    }
    .extree_node(id = 1, data = data, updatefun = updatefun, selectfun =
    selectfun, splitfun = splitfun, partyvars = partyvars, weights = weights, subset = subset, ctrl = ctrl)
}
