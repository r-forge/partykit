
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


    ### depth keeps track of the depth of the tree
    ### which has to be < than maxdepth
    if (is.null(cenv)) {
        cenv <- new.env()
        depth <- 0L
        assign("maxid", id, envir = cenv)
    } else {
        depth <- get("depth", envir = cenv)
        assign("maxid", id, envir = cenv)
        if (depth >= ctrl$maxdepth)
            return(partynode(as.integer(id)))
    }

    ### check for stumps
    if (id > 1 && ctrl$stump) 
        return(partynode(as.integer(id)))

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

    info <- c(list(criterion = p, p.value = pmin), sf[!(names(sf) %in% c("p", "splitfun"))])

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

.urp_tree <- function(call, frame, na.action, control, 
                      growfun, trafofun, scores = NULL, doFit = TRUE)
{

    m <- match(c("formula", "data", "subset", "weights", "na.action"),
               names(call), 0)
    mf <- call[c(1, m)]
    formula <- eval(mf$formula)

    f <- if (inherits(formula, "Formula")) formula else Formula(formula)
    if (length(length(f)) != 2)
        stop("incorrect formula")  
    if (!(length(f)[2] %in% 1:3))
        stop("incorrect formula")
    mf$formula <- f
    mf$drop.unused.levels <- FALSE      
    mf$na.action <- na.pass
    mf[[1]] <- quote(stats::model.frame)
    mf1 <- eval(mf, frame)
    mfterms <- terms(mf1)
    fdot <- attr(mfterms, "Formula_without_dot")

    weights <- model.weights(mf1)
    if (is.null(weights)) weights <- integer(0)
    mf1[["(weights)"]] <- NULL
    if (!is.null(o <- attr(attr(mf1, "terms"), "offset")))
        mf1[[attr(attr(mf1, "terms"), "offset")]] <- NULL

    av <- all.vars(f)
    av <- av[av != "."]
    if (!all(av %in% colnames(mf1))) {
         mf[[1]] <- quote(stats::get_all_vars)
         mf$drop.unused.levels <- NULL
         mf$na.action <- NULL
         mf$weights <- NULL
         mf2 <- eval(mf, frame)
         mf2 <- mf2[, !(colnames(mf2) %in% colnames(mf1)), drop = FALSE]
         mf1 <- cbind(mf1, mf2)
    }
    mf <- na.action(mf1)

    if (length(f)[2] == 1) { ### y ~ z _or_ y ~ .
        if (is.null(fdot)) fdot <- f
        modelf <- formula(fdot, lhs = 1, rhs = 0)
        partf <- formula(fdot, lhs = 0, rhs = 1)
        blockf <- NULL
    } else if (length(f)[2] == 2) { ### y ~ x | z
        if (!is.null(fdot))
            stop("dots are not allowed in multipart formulas")
        modelf <- formula(f, lhs = 1, rhs = 1)
        partf <- formula(f, lhs = 0, rhs = 2)
        blockf <- NULL
    } else if (length(f)[2] == 3) { ### y ~ x | z | block
        if (!is.null(fdot))
            stop("dots are not allowed in multipart formulas")
        modelf <- formula(f, lhs = 1, rhs = 1)
        partf <- formula(f, lhs = 0, rhs = 2)
        blockf <- formula(f, lhs = 0, rhs = 3)
    }
    zvars <- rownames(attr(terms(partf, data = mf), "factors"))

    ### <FIXME> should be xtrafo
    if (!is.null(scores)) {
        for (n in names(scores)) {
            sc <- scores[[n]]
            if (is.ordered(mf[[n]]) && 
                nlevels(mf[[n]]) == length(sc)) {
                attr(mf[[n]], "scores") <- as.numeric(sc)
            } else {
                warning("scores for variable ", sQuote(n), " ignored")
            }
        }
    }
    #### </FIXME>

    block <- NULL
    if (!is.null(blockf)) {
        block <- model.frame(blockf, data = mf)
        if (length(block) != 1 || !is.factor(block[[1]]))
            stop("block is not a single factor")
    }

    ### returns a _function_ (trafo, subset, weights)
    ### for growing the tree
    treefun <- growfun(mf, partyvars = match(zvars, colnames(mf)), 
                       block = block, ctrl = control)
    ### returns a _function_ (subset) for computing estfun (essentially)
    trafo <- trafofun(modelf, data = mf, weights = weights, block = block, 
                      ctrl = control)
    
    if (length(weights) > 0)    
        mf[["(weights)"]] <- weights    
    ret <- list(treefun = treefun, trafo = trafo, mf = mf, terms = mfterms)
    if (!doFit)
        return(ret)

    ### grow the tree
    tree <- treefun(trafo, subset = 1:nrow(mf), weights)
    ret$nodes <- tree
    return(ret)
}

.urp_control <- function(criterion, logmincriterion, minsplit = 20L,
                         minbucket = 7L, minprob = 0.01, stump = FALSE,
                         maxsurrogate = 0L, mtry = Inf,
                         maxdepth = Inf, multiway = FALSE, splittry = 2L,
                         majority = FALSE, applyfun = NULL, cores = NULL) {

    ## apply infrastructure for determining split points
    if (is.null(applyfun)) {
        applyfun <- if(is.null(cores)) {
            lapply
        } else {
            function(X, FUN, ...)
                parallel::mclapply(X, FUN, ..., mc.cores = cores)
        }
    }

    if (multiway & maxsurrogate > 0)
        stop("surrogate splits currently not implemented for multiway splits")

    list(criterion = criterion, logmincriterion = logmincriterion,
         minsplit = minsplit, minbucket = minbucket, 
         minprob = minprob, stump = stump, mtry = mtry,
         maxdepth = maxdepth, multiway = multiway, splittry = splittry,
         maxsurrogate = maxsurrogate, majority = majority,
         applyfun = applyfun)
}
