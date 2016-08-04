
### unbiased recursive partitioning: surrogate splits
.urp_surrogates <- function(split, data, subset, whichvar, selectfun, ctrl) {

    sf <- selectfun(y = split, subset = subset, whichvar = whichvar)
    p <- sf$p
    ### partykit always used p-values, so expect some differences
    crit <- p[ctrl$criterion,,drop = TRUE]
    ### crit is maximised, but there might be ties
    ties <- which(abs(crit - max(crit, na.rm = TRUE)) < .Machine$double.eps)
    if (length(ties) > 1) {
        ### add a small value (< 1/1000) to crit derived from order of teststats
        crit[ties] <- crit[ties] + order(p["statistic", ties]) / (sum(ties) * 1000)
    }

    ret <- vector(mode = "list", length = min(c(length(whichvar), ctrl$maxsurrogate)))

    for (i in 1L:length(ret)) {
        jsel <- which.max(crit)
        sp <- sf$splitfun(whichvar = jsel, 0L)
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


### unbiased recursive partitioning: set up new node
.urp_node <- function(id = 1, data, selectfun, inputs, 
                      weights = integer(0), subset, ctrl, cenv = NULL) {

    if (id > 1 && ctrl$stump) return(partynode(as.integer(id)))

    if (is.null(cenv)) {
        cenv <- new.env()
        depth <- 0L
    } else {
        depth <- get("depth", envir = cenv)
        if (depth >= ctrl$maxdepth)
            return(partynode(as.integer(id)))
    }

    if (length(weights) > 0) {
        sw <- sum(weights[subset])
    } else {
        sw <- length(subset)
    }
    if (sw < ctrl$minsplit) return(partynode(as.integer(id)))

    inp <- inputs
    if (ctrl$mtry < Inf) {
        mtry <- min(sum(inp), ctrl$mtry)
        ### sum(inp) == 1 will lead to sample.int instead of sample; see example(sample)
        resample <- function(x, ...) x[sample.int(length(x), ...)]
        s <- resample(which(inp), mtry)
        inp <- logical(length(inp))
        inp[s] <- TRUE
    } 

    sf <- selectfun(subset = subset, whichvar = which(inp))
    p <- sf$p
    crit <- p[ctrl$criterion,,drop = TRUE]
    crit[is.na(crit)] <- -Inf
    ### crit is maximised, but there might be ties
    ties <- which(abs(crit - max(crit, na.rm = TRUE)) < .Machine$double.eps)
    if (length(ties) > 1) {
        ### add a small value (< 1/1000) to crit derived from order of teststats
        crit[ties] <- crit[ties] + order(p["statistic", ties]) / (sum(ties) * 1000)
    }

    ### format p values
    fmP <- function(p) {
        if (all(is.na(p["p.value",]))) return(NA)
        ret <- 1 - exp(p["p.value", which.max(crit)])
        names(ret) <- colnames(data)[which.max(crit)]
        ret
    }

    if (all(crit < ctrl$mincriterion)) {
        return(partynode(as.integer(id), 
                         info = list(criterion = p,
                         p.value = fmP(p))))
    }

    mb <- ctrl$minbucket
    mp <- ctrl$minprob
    swp <- ceiling(sw * mp)
    if (mb < swp) mb <- as.integer(swp)
    jsel <- rev(order(crit))[1:ctrl$splittry]
    jsel <- jsel[crit[jsel] > ctrl$mincriterion]
    thissplit <- sf$splitfun(whichvar = jsel, minbucket = mb)

    if (is.null(thissplit))
        return(partynode(as.integer(id), 
                         info = list(criterion = p,
                                     p.value = fmP(p))))           

    ret <- partynode(as.integer(id))
    ret$split <- thissplit
    ret$info <- list(criterion = p[,inp], p.value = fmP(p))
    thissurr <- NULL
    s <- subset[!is.na(data[subset, varid_split(thissplit)])]
    kidids <- kidids_node(ret, data, obs = s)
    prob <- prop.table(table(kidids))
    names(dimnames(prob)) <- NULL
    if (ctrl$majority)  ### go with majority
        prob <- as.double((1L:length(prob)) %in% which.max(prob))
    ret$split$prob <- prob

    if (ctrl$maxsurrogate > 0) {
        inp <- inputs
        inp[thissplit$varid] <- FALSE
        ret$surrogates <- .urp_surrogates(kidids, data = data, subset = s, 
                                          whichvar = which(inp), 
                                          selectfun = selectfun, 
                                          ctrl = ctrl)
    }
    ### <FIXME> we do this twice, not really needed </FIXME>
    kidids <- kidids_node(ret, data, obs = subset)

    kids <- vector(mode = "list", length = max(kidids)) ## Z: was 1:max(kidids)
    nextid <- id + 1
    for (k in 1L:max(kidids)) {
        nextsubset <- subset[kidids == k]
        assign("depth", depth + 1, envir = cenv)
        kids[[k]] <- .urp_node(id = nextid, data = data, selectfun = selectfun, 
                               inputs = inputs, weights = weights, 
                               subset = nextsubset, ctrl = ctrl, cenv = cenv)
        nextid <- max(nodeids(kids[[k]])) + 1
    }
    ret$kids <- kids

    return(ret)
}
