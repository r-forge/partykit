
### surrogate splits
.csurr <- function(split, data, subset, whichvar, selectfun, splitfun, ctrl) {

    p <- selectfun(y = split, subset = subset, whichvar = whichvar)
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
        sp <- splitfun(y = split, subset = subset, whichvar = jsel, 0L)
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


### set up new node for conditional inference tree
.cnode <- function(id = 1, data, selectfun, splitfun, inputs, 
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

    p <- selectfun(subset = subset, whichvar = which(inp))
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
    thissplit <- splitfun(subset = subset, whichvar = jsel, minbucket = mb)

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
        ret$surrogates <- .csurr(kidids, data = data, subset = s, whichvar = which(inp), 
                                 selectfun = selectfun, splitfun = splitfun, ctrl)
    }
    ### <FIXME> we do this twice, not really needed </FIXME>
    kidids <- kidids_node(ret, data, obs = subset)

    kids <- vector(mode = "list", length = max(kidids)) ## Z: was 1:max(kidids)
    nextid <- id + 1
    for (k in 1L:max(kidids)) {
        nextsubset <- subset[kidids == k]
        assign("depth", depth + 1, envir = cenv)
        kids[[k]] <- .cnode(id = nextid, data = data, selectfun = selectfun, 
                            splitfun = splitfun, inputs = inputs, weights = weights, 
                            subset = nextsubset, ctrl = ctrl, cenv = cenv)
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
    }

    selectfun <- function(y = NULL, subset, whichvar) 
    {
    
        ret <- matrix(NA, nrow = 2L, ncol = ncol(data))
        colnames(ret) <- names(data)
        rownames(ret) <- c("statistic", "p.value")

        if (is.null(y)) {
            if (!is.null(trafo))
                Y[subset,] <- trafo(data = data, subset = subset, weights = weights)
            iy <- bdr[[response]]
        } else {
            ### y is kidids in .csurr and nothing else
            stopifnot(length(y) == length(subset))
            if (ctrl$nmax < Inf) {
                Y <- rbind(0, max(y))
                iy <- numeric(NROW(data))
                iy[subset] <- as.integer(y)
            } else {
                Y <- matrix(0, nrow = NROW(data), ncol = max(y))
                Y[cbind(subset, y)] <- 1 ### model.matrix(~ as.factor(y) - 1)
                storage.mode(Y) <- "double"
                iy <- NULL
            }
        }

        for (j in whichvar) {
            lev <- LinStatExpCov(X = X[[j]], ix = bdr[[j]], 
                                 Y = Y, iy = iy, subset = subset,
                                 weights = weights, block = block, 
                                 B = ifelse(ctrl$testtype == "MonteCarlo", 
                                            ctrl$nresample, 0L))
            tst <- doTest(lev, teststat = ctrl$teststat, 
                          pvalue = ctrl$testtype != "Teststatistic",
                          lower = TRUE, 
                          log = TRUE)
            ret["statistic", j] <- log(tst$TestStatistic)
            ret["p.value", j] <- tst$p.value
        }
        if (ctrl$testtype == "Bonferroni")
            ret["p.value",] <- ret["p.value",] * length(whichvar)
        ret
    }

    splitfun <- function(y = NULL, subset, whichvar,  minbucket)
    {
        ret <- NULL
        if (is.null(y)) {
            if (!is.null(trafo))
                Y[subset,] <- trafo(data = data, subset = subset, weights = weights)
            iy <- bdr[[response]]
        } else {
            ### y is kidids in .csurr and nothing else
            stopifnot(length(y) == length(subset))
            if (ctrl$nmax < Inf) {
                Y <- rbind(0, max(y))
                iy <- numeric(NROW(data))
                iy[subset] <- as.integer(y)
            } else {
                Y <- matrix(0, nrow = NROW(data), ncol = max(y))
                Y[cbind(subset, y)] <- 1 ### model.matrix(~ as.factor(y) - 1)
                storage.mode(Y) <- "double"
                iy <- NULL
            }
        }

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
                        x[-subset] <- NA
                        X <- match(x, ux <- sort(unique(x)))
                        X[is.na(X)] <- 0L
                        storage.mode(X) <- "integer"
                    }
                }
                lev <- LinStatExpCov(X = X, ix = bdr[[j]],
                                Y = Y, iy = iy, subset = subset,
                                weights = weights, block = block, B = 0L, varonly = TRUE)
                sp <- doTest(lev, teststat = ctrl$splitstat,
                             minbucket = minbucket, pvalue = FALSE,
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

    tree <- .cnode(id = 1L, data = data, selectfun = selectfun, splitfun = splitfun, 
                   inputs = inputs, weights = weights, subset = subset, ctrl = ctrl)

    return(tree)
}


ctree_control <- function(teststat = c("quadratic", "maximum"), splitstat = c("maximum", "quadratic"),
    testtype = c("Bonferroni", "MonteCarlo", "Univariate", "Teststatistic"),
    nmax = 20L,
    mincriterion = 0.95, minsplit = 20L, minbucket = 7L, minprob = 0.01,
    stump = FALSE, nresample = 9999L, maxsurrogate = 0L, mtry = Inf, maxdepth = Inf, 
    multiway = FALSE, splittry = 2L, majority = FALSE,
    applyfun = NULL, cores = NULL) {

    teststat <- match.arg(teststat)
    splitstat <- match.arg(splitstat)
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

    list(teststat = teststat, splitstat = splitstat, 
         criterion = ifelse(testtype == "Teststatistic", "statistic", "p.value"),
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
