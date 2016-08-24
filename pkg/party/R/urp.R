
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
    weights = integer(0),	### optional case weights
    subset, 			### subset of 1:nrow(data)
                                ### for identifying obs for this node
    ctrl, 			### ctree_control()
    ### needs: stump, maxdepth, caseweights, minsplit, mtry, criterion,
    ###        logmincriterion, minbucket, minprob, splittry,
    ###        majority, maxsurrogate
    cenv = NULL			### environment for depth and maxid
) {

    ### check for stumps
    if (id > 1 && ctrl$stump) 
        return(partynode(as.integer(id)))

    ### depth keeps track of the depth of the tree
    ### which has to be < than maxdepth
    if (is.null(cenv)) {
        cenv <- new.env()
        depth <- 0L
    } else {
        depth <- get("depth", envir = cenv)
        if (depth >= ctrl$maxdepth)
            return(partynode(as.integer(id)))
    }
    assign("maxid", id, envir = cenv)

    ### sw is basically the number of observations
    ### which has to be > minsplit in order to consider
    ### the node for splitting
    ### <FIXME> check ctrl$caseweights
    if (length(weights) > 0) {
        sw <- sum(weights[subset])
    } else {
        sw <- length(subset)
    }
    ### </FIXME>
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
    p <- sf$p
    crit <- p[ctrl$criterion,,drop = TRUE]
    crit[is.na(crit)] <- -Inf
    ### crit is maximised, but there might be ties
    ties <- which(abs(crit - max(crit, na.rm = TRUE)) < .Machine$double.eps)
    if (length(ties) > 1) {
        ### add a small value (< 1/1000) to crit derived from order of teststat
        crit[ties] <- crit[ties] + 
            order(p["statistic", ties]) / (sum(ties) * 1000)
    }

    ### switch from log(1 - pval) to pval
    p <- rbind(p, criterion = crit)
    p["p.value",] <- -expm1(p["p.value",])
    pmin <- p["p.value", which.max(crit)]
    names(pmin) <- colnames(data)[which.max(crit)]

    ### nothing "significant"
    if (all(crit < ctrl$logmincriterion)) {
        return(partynode(as.integer(id), 
                         info = list(criterion = p,
                                     p.value = pmin)))
    }

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
        return(partynode(as.integer(id), 
                         info = list(criterion = p,
                                     p.value = pmin)))           

    ### successful split search: set-up node
    ret <- partynode(as.integer(id))
    ret$split <- thissplit
    ret$info <- list(criterion = p, p.value = pmin)
    thissurr <- NULL

    ### determine observations for splitting (only non-missings)
    snotNA <- subset[!is.na(data[subset, varid_split(thissplit)])]
    ### and split observations
    kidids <- kidids_node(ret, data, obs = snotNA)

    ### compute probability of going left / right
    prob <- tabulate(kidids) / length(kidids) ### was: prop.table(table(kidids))
    # names(dimnames(prob)) <- NULL
    if (ctrl$majority)  ### go with majority
        prob <- as.double((1L:length(prob)) %in% which.max(prob))
    ret$split$prob <- prob

    ### compute surrogate splits
    if (ctrl$maxsurrogate > 0)
        ret$surrogates <- .urp_surrogates(kidids, data = data, 
            subset = snotNA, 
            partyvars = svars[svars != varid_split(thissplit)],
            selectfun = svselectfun, ctrl = ctrl)
    kidids <- kidids_node(ret, data, obs = subset)

    ### proceed recursively
    kids <- vector(mode = "list", length = max(kidids)) 
    nextid <- id + 1
    for (k in 1L:max(kidids)) {
        nextsubset <- subset[kidids == k]
        assign("depth", depth + 1, envir = cenv)
        kids[[k]] <- .urp_node(id = nextid, data = data, 
            selectfun = selectfun, partyvars = partyvars, 
            weights = weights, subset = nextsubset, 
            ctrl = ctrl, cenv = cenv)
        ### was: nextid <- max(nodeids(kids[[k]])) + 1
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
    p <- sf$p
    ### partykit always used p-values, so expect some differences
    crit <- p[ctrl$criterion,,drop = TRUE]
    ### crit is maximised, but there might be ties
    ties <- which(abs(crit - max(crit, na.rm = TRUE)) < .Machine$double.eps)
    if (length(ties) > 1) {
        ### add a small value (< 1/1000) to crit derived from order of teststat
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
    if (length(ret) == 0) ret <- NULL
    return(ret)
}

