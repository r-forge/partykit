
### surrogate splits
.csurr <- function(split, data, subset, whichvar, ctrl) {

    if (ctrl$nmax < Inf) {
        Y <- rbind(0, max(split))
        iy <- split
    } else {
        Y <- matrix(NA, nrow = length(split), ncol = max(split))
        split[split == 0] <- NA
        Y[!is.na(split),] <- model.matrix(~ factor(split) - 1)
        storage.mode(Y) <- "double"
        iy <- NULL
    }

    p <- ctrl$var_select(Y, iy = iy, subset = subset, whichvar = whichvar)
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
        sp <- ctrl$var_split(Y, iy = iy, subset = subset, whichvar = jsel, 0L)
        if (is.null(sp)) next
        ret[[i]] <- sp
        tmp <- kidids_split(ret[[i]], data, obs = subset)
        tmps <- split[subset]
        tab <- table(tmp, tmps)
        if (tab[1, 1] < tab[1, 2]) {
            indx <- ret[[i]]$index
            ret[[i]]$index[indx == 1] <- 2L
            ret[[i]]$index[indx == 2] <- 1L
        }
        crit[which.max(crit)] <- -Inf
    }
    ret <- ret[!sapply(ret, is.null)]
    if (length(ret) == 0) ret <- NULL
    return(ret)
}


### set up new node for conditional inference tree
.cnode <- function(id = 1, data, trafo, inputs, weights = integer(0), 
                   subset, ctrl, cenv = NULL) {

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

    Y <- trafo(subset)
    p <- ctrl$var_select(Y, subset = subset, whichvar = which(inp))
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
    thissplit <- ctrl$var_split(Y, subset = subset, whichvar = jsel , minbucket = mb)

    if (is.null(thissplit))
        return(partynode(as.integer(id), 
                         info = list(criterion = p,
                                     p.value = fmP(p))))           

    ret <- partynode(as.integer(id))
    ret$split <- thissplit
    ret$info <- list(criterion = p[,inp], p.value = fmP(p))
    thissurr <- NULL
    kidids <- kidids_node(ret, data, obs = subset)
    prob <- prop.table(table(kidids))
    names(dimnames(prob)) <- NULL
    if (ctrl$majority)  ### go with majority
        prob <- (1L:length(prob)) %in% which.max(prob)
    ret$split$prob <- prob

    if (ctrl$maxsurrogate > 0) {
        inp <- inputs
        inp[thissplit$varid] <- FALSE
        tmp <- integer(nrow(data))
        tmp[subset] <- kidids
        ret$surrogates <- .csurr(tmp, data = data, subset = subset, whichvar = which(inp), ctrl)
        kidids <- kidids_node(ret, data, obs = subset)
    }

    kids <- vector(mode = "list", length = max(kidids)) ## Z: was 1:max(kidids)
    nextid <- id + 1
    for (k in 1L:max(kidids)) {
        nextsubset <- subset[kidids == k]
        assign("depth", depth + 1, envir = cenv)
        kids[[k]] <- .cnode(nextid, data, trafo, inputs, weights = weights, 
                            subset = nextsubset, ctrl, cenv)
        nextid <- max(nodeids(kids[[k]])) + 1
    }
    ret$kids <- kids

    return(ret)
}

.ctree_fit <- function(data, response, trafo = NULL, 
                       weights = integer(0), subset, block = integer(0), ctrl)
{

    inputs <- !(colnames(data) %in% response)


    if (missing(subset)) subset <- 1:NROW(data)

    if (ctrl$nmax < Inf) {
        bdr <- libcoin:::BDR(data, ctrl$nmax)
        X <- lapply(bdr, function(x) attr(x, "X"))
    } else {
        Y <- partykit:::.y2infl(data, response, ytrafo = trafo)
        X <- lapply(data, function(x) {
            if (is.numeric(x)) {
                ret <- matrix(x, ncol = 1)
            } else if (is.ordered(x)) {
                sc <- 1:nlevels(x)
                if (!is.null(attr(x, "scores")))
                    sc <- attr(x, "scores")
                ret <- matrix(sc[x], ncol = 1)
            } else if (is.factor(x)) {
                ret <- model.matrix( ~ x - 1)
            } else {
                stop("cannot handle class", class(x))
            }
            storage.mode(ret) <- "double"
            ret
        })
        bdr <- vector(mode = "list", length = ncol(data))
        names(bdr) <- colnames(data)
    }
    if (is.null(trafo)) {
        Y <- X[[response]]
        trafo <- function(...) Y
    }

    ctrl$var_select <- function(Y, iy =  bdr[[response]], subset, whichvar) {
    
        ret <- matrix(NA, nrow = 2L, ncol = ncol(data))
        colnames(ret) <- names(data)
        rownames(ret) <- c("statistic", "p.value")
        
        for (j in whichvar) {
            lev <- LinStatExpCov(X = X[[j]], ix = bdr[[j]], 
                                 Y = Y, iy = iy, subset = subset,
                                 weights = weights, block = block, 
                                 B = ifelse(ctrl$testtype == "MonteCarlo", 
                                            ctrl$nresample, 0L))
            tst <- doTest(lev, type = ifelse(ctrl$teststat == "quad", "quadform", "maxstat"), 
                          lower = TRUE, 
                          log = TRUE)
            ret["statistic", j] <- log(tst$TestStatistic)
            ret["p.value", j] <- tst$p.value
        }
        if (ctrl$testtype == "Bonferroni")
            ret["p.value",] <- ret["p.value",] * length(whichvar)
        ret
     }

     ctrl$var_split <- function(Y, iy = bdr[[response]], subset, whichvar, 
                                minbucket)
     {
         ret <- NULL
         for (j in whichvar) {
             x <- data[[j]]
             ORDERED <- is.ordered(x) || is.numeric(x)
             if ((ctrl$multiway && ctrl$maxsurrogate == 0) && is.factor(x)) {
                 index <- 1L:nlevels(x)
                 if (length(weights) > 0) {
                     xt <- xtabs(~ x, subset = subset)
                 } else {
                     xt <- xtabs(weights ~ x, subset = subset)
                 }
                 index[xt == 0] <- NA
                 index[xt < mb] <- nlevels(x) + 1L
                 index <- unclass(factor(index))
                 ret <- partysplit(as.integer(j), index = as.integer(index))
                 break()
            } else {
                 ix <- bdr[[j]]
                 if (!is.null(ix)) {
                     X <- numeric(0) 
                     ux <- attr(bdr[[j]], "levels")
                 } else {
                     if (is.factor(x)) {
                         X <- unclass(x)
                     } else {
                         X <- match(x, ux <- sort(unique(x)))
                         X[is.na(X)] <- 0L
                         storage.mode(X) <- "integer"
                     }
                 }
                 lev <- LinStatExpCov(X = X, ix = bdr[[j]],
                                 Y = Y, iy = iy, subset = subset,
                                 weights = weights, block = block, B = 0L)
                 sp <- doTest(lev, type = ifelse(ctrl$teststat == "quad", "quadform", "maxstat"), 
                              minbucket = minbucket, 
                              ordered = ORDERED)$index
                 if (!all(is.na(sp))) {
                     if (length(sp) == 1) {
                         if (!is.ordered(x))
                             sp <- ux[sp]
                         ret <- partysplit(as.integer(j), breaks = sp, index = 1L:2L)
                     } else {
                         ret <- partysplit(as.integer(j), index = as.integer(sp) + 1L)
                     }
                     break
                 }
            }
        }
        ret
    }

    tree <- .cnode(1L, data, trafo, inputs, weights, subset, ctrl)

    return(tree)
}


ctree_control <- function(teststat = c("quad", "max"),
    testtype = c("Bonferroni", "MonteCarlo", "Univariate", "Teststatistic"),
    nmax = 20L,
    mincriterion = 0.95, minsplit = 50L, minbucket = 7L, minprob = 0.01,
    stump = FALSE, nresample = 9999L, maxsurrogate = 0L, mtry = Inf, maxdepth = Inf, 
    multiway = FALSE, splittry = 2L, majority = FALSE,
    applyfun = NULL, cores = NULL) {

    teststat <- match.arg(teststat)
    testtype <- match.arg(testtype)

    ## apply infrastructure for determining split points
    if (is.null(applyfun)) {
        applyfun <- if(is.null(cores)) {
            lapply
        } else {
            function(X, FUN, ...) 
                parallel::mclapply(X, FUN, ..., mc.cores = cores)
        }
    }

    list(teststat = teststat, criterion = ifelse(testtype == "Teststatistic", "statistic", "p.value"),
         testtype = testtype, nmax = nmax, mincriterion = log(mincriterion),
         minsplit = minsplit, minbucket = minbucket, 
         minprob = minprob, stump = stump, nresample = nresample, mtry = mtry, 
         maxdepth = maxdepth, multiway = multiway, splittry = splittry,
         maxsurrogate = maxsurrogate, majority = majority, 
         applyfun = applyfun)
}

ctree <- function(formula, data, weights, subset, na.action = na.pass, 
                  control = ctree_control(...), ytrafo = NULL, 
                  scores = NULL, ...) {

    if (missing(data))
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action"),
               names(mf), 0)
    mf <- mf[c(1, m)]
    
    ### only necessary for extended model formulae 
    ### e.g. multivariate responses
    formula <- Formula::Formula(formula)
    mf$formula <- formula
    mf$drop.unused.levels <- FALSE
    mf$na.action <- na.action
    mf[[1]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())

    response <- names(Formula::model.part(formula, mf, lhs = 1))
    weights <- model.weights(mf)
    dat <- mf[, colnames(mf) != "(weights)"]
    if (!is.null(scores)) {
        for (n in names(scores)) {
            sc <- scores[[n]]
            if (is.ordered(dat[[n]]) && 
                nlevels(dat[[n]]) == length(sc)) {
                attr(dat[[n]], "scores") <- as.numeric(sc)
            } else {
                warning("scores for variable ", sQuote(n), " ignored")
            }
        }
    }

    ### <FIXME> implement y ~ x | block or y ~ 1 | x | block ? </FIXME>
    block <- integer(0)
    tree <- .ctree_fit(dat, response = response, trafo = ytrafo, 
                       weights = weights, block = block, ctrl = control)

    if (is.null(weights)) weights <- rep(1, nrow(dat))
    fitted <- data.frame("(fitted)" = fitted_node(tree, dat), 
                         "(weights)" = weights,
                         check.names = FALSE)
    fitted[[3]] <- dat[, response, drop = length(response) == 1]
    names(fitted)[3] <- "(response)"
    ret <- party(tree, data = dat, fitted = fitted, 
                 info = list(call = match.call(), control = control))
    class(ret) <- c("constparty", class(ret))

    ### doesn't work for Surv objects
    # ret$terms <- terms(formula, data = mf)
    ret$terms <- terms(mf)
    ### need to adjust print and plot methods
    ### for multivariate responses
    ### if (length(response) > 1) class(ret) <- "party"
    return(ret)
}
