
.ctree_test_split <- function(x, bdr = NULL, j, ctrl, X, Y, iy = NULL, subset, 
                              weights, cluster, splitonly = TRUE, minbucket) {

    ### <FIXME> check if x and/or Y have only unique values and return immediately? </FIXME>

    ### <FIXME> MIA splits only estimated in the presence of missings
    ###         new missings in predict(<>, newdata) will cause trouble
    ### </FIXME>
    MIA <- ctrl$MIA && any(is.na(x[subset]))

    if (all(is.na(x[subset]))) { ### all x values are missing
        if (splitonly) return(NULL)
        return(list(statistic = NA, p.value = NA))
    }

    if (is.null(cluster)) cluster <- integer(0)
    if (splitonly) {
        if ((ctrl$multiway && ctrl$maxsurrogate == 0) &&
            is.factor(x) && nlevels(x[subset, drop = TRUE]) > 1) {
            index <- 1L:nlevels(x)
            if (length(weights) > 0) {
                xt <- xtabs(weights ~ x, subset = subset)
            } else {
                xt <- xtabs(~ x, subset = subset)
            }
            index[xt == 0] <- NA
            index[xt < minbucket] <- nlevels(x) + 1L
            index <- unclass(factor(index))
            return(partysplit(as.integer(j),
                              index = as.integer(index)))
        }
    }

    X <- X[[j]]
    ux <- NULL
    ORDERED <- is.ordered(x) || is.numeric(x)
    if (is.null(bdr)) {
        ix <- NULL
        if (ctrl$splittest || splitonly) {
            ### integer X trigger maximally selected stats
            if (is.factor(x)) {
                X <- unclass(x)
                attr(X, "levels") <- levels(x)
            } else {
                x[-subset] <- NA
                ux <- sort(unique(x))
                X <- cut.default(x, breaks = c(-Inf, ux, Inf),
                                 labels = FALSE, right = TRUE)
                # X[is.na(X)] <- 0L (NAs are handled by LinStatExpCov)
                attr(X, "levels") <- ux 
                storage.mode(X) <- "integer"
            }
            if (MIA) {
                Xlev <- attr(X, "levels")
                Xleft <- X + 1L
                Xleft[is.na(Xleft)] <- 1L
                Xright <- X
                Xright[is.na(Xright)] <- as.integer(length(Xlev) + 1L)
                attr(Xleft, "levels") <- c(NA, Xlev)
                attr(Xright, "levels") <- c(Xlev, NA)
                ixleft <- ixright <- ix
            }
        } else {
            MIA <- FALSE
        } 
    } else {
        ix <- bdr[[j]]
        if (ctrl$splittest || splitonly) {
            X <- numeric(0) 
            ux <- attr(ix, "levels")
            if (MIA) {
                Xlev <- attr(ix, "levels")
                ixleft <- ix + 1L
                ixright <- ix
                ixright[ixright == 0L] <- as.integer(length(Xlev) + 1L)
                attr(ixleft, "levels") <- c(NA, Xlev)
                attr(ixright, "levels") <- c(Xlev, NA)
                Xleft <- Xright <- X
           }
        } else {
            MIA <- FALSE
        }
    }

    if (splitonly) {
        B <- 0L
        varonly <- TRUE
        pvalue <- FALSE
        teststat <- ctrl$splitstat
    } else {
        B <- ifelse(ctrl$testtype == "MonteCarlo",
                    ctrl$nresample, 0L)
        if (ctrl$splittest) {
            if (ctrl$teststat != ctrl$splitstat)
                warning("Using different test statistics for testing and splitting")
            teststat <- ctrl$splitstat
            if (B == 0) 
                stop("MonteCarlo approximation mandatory for splittest = TRUE")
        } else {
            teststat <- ctrl$teststat
        }
        varonly <- ctrl$testtype == "MonteCarlo" && 
                   teststat == "maxtype"
        pvalue <- ctrl$testtype != "Teststatistic"
    }
    ### see libcoin/src/C_ordered_Xfactor_block
    if (length(cluster) > 0) varonly <- FALSE 

    ### if (MIA) use tst as fallback
    ### compute linear statistic + expecation and covariance
    lev <- LinStatExpCov(X = X, Y = Y, ix = ix, iy = iy, subset = subset,
                         weights = weights, block = cluster,
                         B = B, varonly = varonly)
    ### compute test statistic and log(1 - p-value)
    tst <- doTest(lev, teststat = teststat, pvalue = pvalue,
                  lower = TRUE, log = TRUE, ordered = ORDERED, 
                  minbucket = minbucket, pargs = ctrl$pargs)

    if (MIA) {
        ### compute linear statistic + expecation and covariance
        lev <- LinStatExpCov(X = Xleft, Y = Y, ix = ixleft, iy = iy, subset = subset,
                             weights = weights, block = cluster,
                             B = B, varonly = varonly)
        ### compute test statistic and log(1 - p-value)
        tstleft <- doTest(lev, teststat = teststat, pvalue = pvalue,
                          lower = TRUE, log = TRUE, ordered = ORDERED, 
                          minbucket = minbucket, pargs = ctrl$pargs)
        ### compute linear statistic + expecation and covariance
        lev <- LinStatExpCov(X = Xright, Y = Y, ix = ixright, iy = iy, subset = subset,
                             weights = weights, block = cluster,
                             B = B, varonly = varonly)
        ### compute test statistic and log(1 - p-value)
        tstright <- doTest(lev, teststat = teststat, pvalue = pvalue,
                           lower = TRUE, log = TRUE, ordered = ORDERED, 
                           minbucket = minbucket, pargs = ctrl$pargs)
    }

    if (!splitonly) {
        if (MIA) {
            tst <- tstleft
            if (tst$TestStatistic < tstright$TestStatistic)
                tst <- tstright
        }
        return(list(statistic = log(pmax(tst$TestStatistic, .Machine$double.eps)), 
                    p.value = tst$p.value))
    }

    ret <- NULL
    if (MIA) {
        if (ORDERED) {
            if (tstleft$TestStatistic >= tstright$TestStatistic) {
                if (tst$index == 1) { ### case C
                    ret <- partysplit(as.integer(j), breaks = Inf, 
                                      index = 1L:2L, prob = as.double(0:1))
                } else {
                    sp <- tstleft$index - 1L ### case A
                    if (!is.ordered(x)) {
                        ### interpolate split-points, see https://arxiv.org/abs/1611.04561
                        if (ctrl$intersplit & sp < length(ux)) {     
                            sp <- (ux[sp] + ux[sp + 1]) / 2
                        } else {
                            sp <- ux[sp]  ### X <= sp vs. X > sp
                        }
                    }
                    ret <- partysplit(as.integer(j), breaks = sp,
                                      index = 1L:2L, prob = as.double(rev(0:1)))
                }
            } else {
                ### case C was handled above (tstleft = tstright in this case)
                sp <- tstright$index ### case B
                if (!is.ordered(x)) {
                    ### interpolate split-points, see https://arxiv.org/abs/1611.04561
                    if (ctrl$intersplit & sp < length(ux)) {     
                        sp <- (ux[sp] + ux[sp + 1]) / 2
                    } else {
                        sp <- ux[sp]  ### X <= sp vs. X > sp
                    }
                }
                ret <- partysplit(as.integer(j), breaks = sp,
                                  index = 1L:2L, prob = as.double(0:1))
            }
        } else {
            sp <- tstleft$index[-1L] ### tstleft = tstright for unordered factors
            if (all(sp == 0)) { ### case C
                ret <- partysplit(as.integer(j), index = as.integer(tst$index) + 1L)
            } else { ### always case A
                ret <- partysplit(as.integer(j),
                                  index = as.integer(sp) + 1L, 
                                  prob = as.double(rev(0:1)))
            }
        }
    } else {
        sp <- tst$index
        if (all(is.na(sp))) return(NULL)
        if (ORDERED) {
            if (!is.ordered(x))
                ### interpolate split-points, see https://arxiv.org/abs/1611.04561
                if (ctrl$intersplit & sp < length(ux)) {
                    sp <- (ux[sp] + ux[sp + 1]) / 2 
                } else {
                    sp <- ux[sp]  ### X <= sp vs. X > sp
                }
                ret <- partysplit(as.integer(j), breaks = sp,
                                  index = 1L:2L)
        } else {
            ret <- partysplit(as.integer(j),
                              index = as.integer(sp) + 1L)
        }
    }
    return(ret)
}

.ctreetrafo <- function
(
    formula, 
    data,
    ctrl, 
    ytrafo,
    converged = NULL
) {

    weights <- model.weights(data)
    if (is.null(weights)) weights <- integer(0)
    cluster <- data[["cluster"]]
    offset <- model.offset(data)
    if (!is.null(offset)) warning("offset ignored by trafo")

    if (ctrl$nmax < Inf) {
        if (is.function(ytrafo)) 
            return(ytrafo(formula, data = data, weights = weights, 
                          cluster = cluster, ctrl = ctrl))
        f <- Formula(formula)
        mf <- model.frame(formula = f, data = data)
        for (nm in names(mf)[names(mf) %in% names(data)]) {
            sc <- attr(data[[nm]], "scores")
            if (!is.null(sc)) attr(mf[[nm]], "scores") <- sc
        }
        y <- model.part(f, data = mf, lhs = 1, rhs = 0)
        bdr <- inum::inum(y, nmax = ctrl$nmax, total = TRUE, 
                          complete.cases.only = TRUE)
        y <- attr(bdr, "levels")
        index <- c(bdr)
        attr(index, "levels") <- 1:NROW(y)
        cn <- colnames(y)
        ### <FIXME> y lost its scores attribute, readd! </FIXME>
        Y <- .y2infl(y, cn[cn != "(weights)"], ytrafo = ytrafo)
        ### first row corresponds to missings
        Y <- rbind(0, Y)  
        return(function(subset, ...)
            list(estfun = Y, index = index, 
                 converged =  if (is.null(converged)) 
                                  TRUE 
                              else 
                                  converged(Y, mf, subset)))
    } else {
        if (is.function(ytrafo))
            return(ytrafo(formula, data = data, weights = weights, 
                          cluster = cluster, ctrl = ctrl))
        f <- Formula(formula)
        mf <- model.frame(formula = f, data = data, na.action = na.pass)
        for (nm in names(mf)[names(mf) %in% names(data)]) {
            sc <- attr(data[[nm]], "scores")
            if (!is.null(sc)) attr(mf[[nm]], "scores") <- sc
        }
        y <- model.part(f, data = mf, lhs = 1, rhs = 0)
        cc <- complete.cases(y)
        ### do not subset y before calling .y2infl as the scores attribute
        ### would get lost...
        Yi <- .y2infl(y, colnames(y), ytrafo = ytrafo)[cc,,drop = FALSE]
        Y <- matrix(NA, nrow = nrow(mf), ncol = NCOL(Yi))
        Y[cc,] <- Yi
        #    colnames(Y) <- colnames(Yi)
        storage.mode(Y) <- "double"
        return(function(subset, ...)
            list(estfun = Y, converged = if (is.null(converged)) 
                                             TRUE 
                                         else 
                                             converged(Y, mf, subset)))
    }
}

.ctreegrow <- function
(
    data, 
    partyvars, 
    cluster, 
    ctrl
) {

    if (ctrl$nmax < Inf)
        return(.ctree_fit_2d(data = data, partyvars = partyvars, 
                             cluster = cluster, ctrl = ctrl))
    return(.ctree_fit_1d(data = data, partyvars = partyvars,
                         cluster = cluster, ctrl = ctrl))
}

### conditional inference trees
.ctree_fit_1d <- function
(
    data, 				### full data, readonly
    partyvars, 				### partytioning variables,
					### a subset of 1:ncol(data)
    cluster = integer(0), 		### a blocking factor w/o NA
    ctrl				### ctree_control()
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
        } else {
            stop("cannot handle class", class(x))
        }
        storage.mode(ret) <- "double"
        ret
    })

    ### this is the update function
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
                tst <- .ctree_test_split(x = data[[j]], bdr = NULL, j = j, ctrl = ctrl, 
                                         X = X, Y = Y, iy = NULL, subset = subset, 
                                         weights = weights, cluster = cluster,
                                         splitonly = FALSE, minbucket = ctrl$minbucket)
                ### <FIXME> minbucket is updated in .urp_node but only after testing... </FIXME>
                ret$criteria["statistic", j] <- tst$statistic
                ret$criteria["p.value", j] <- tst$p.value
            }
            if (ctrl$testtype == "Bonferroni")
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
                    ret <- .ctree_test_split(x = data[[j]], bdr = NULL, j = j, ctrl = ctrl,
                                             X = X, Y = Y, iy = NULL, subset = subset, 
                                             weights = weights, cluster = cluster,
                                             splitonly = TRUE, minbucket = minbucket)
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
    })
}


### conditional inference trees, faster but approximate version
.ctree_fit_2d <- function
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
                tst <- .ctree_test_split(x = data[[j]], bdr = bdr, j = j, ctrl = ctrl,
                                         X = X, Y = Y, iy = iy, subset = subset,
                                         weights = weights, cluster = cluster,
                                         splitonly = FALSE, minbucket = ctrl$minbucket)
                ret$criteria["statistic", j] <- tst$statistic
                ret$criteria["p.value", j] <- tst$p.value
            }
            if (ctrl$testtype == "Bonferroni")
                ret$criteria["p.value",] <- ret$criteria["p.value",] * 
                    length(whichvar)

            ret <- c(ret, tr[!(names(tr) %in% c("estfun", "index"))])

            ### compute best fitting cutpoint (according to minbucket)
            ### for a subset of observations and variables
            ### y is used for node ids when computing surrogate splits
            ret$splitfun <- function(whichvar,  minbucket)
            {
                for (j in whichvar) {
                    ret <- .ctree_test_split(x = data[[j]], bdr = bdr, j = j, ctrl = ctrl,
                                             X = X, Y = Y, iy = iy, subset = subset,
                                             weights = weights, cluster = cluster,
                                             splitonly = TRUE, minbucket = minbucket)
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

ctree_control <- function
(
    teststat = c("quadratic", "maximum"), 
    splitstat = c("quadratic", "maximum"), ### much better for q > 1, max was default
    splittest = FALSE,
    testtype = c("Bonferroni", "MonteCarlo", 
                 "Univariate", "Teststatistic"),
    pargs = GenzBretz(),
    nmax = Inf, 
    alpha = 0.05, 
    mincriterion = 1 - alpha, 
    logmincriterion = log(mincriterion), 
    minsplit = 20L, 
    minbucket = 7L, 
    minprob = 0.01, 
    stump = FALSE, 
    lookahead = FALSE,	### try trafo() for daugther nodes before implementing the split
    nresample = 9999L, 
    MIA = FALSE,	### DOI: 10.1016/j.patrec.2008.01.010
    maxsurrogate = 0L, 
    mtry = Inf, 
    maxdepth = Inf, 
    multiway = FALSE, 
    splittry = 2L, 
    intersplit = FALSE,
    majority = FALSE, 
    caseweights = TRUE, 
    applyfun = NULL, 
    cores = NULL
) {

    teststat <- match.arg(teststat)
    splitstat <- match.arg(splitstat)
    testtype <- match.arg(testtype)

    if (!caseweights)
        stop("only caseweights currently implemented in ctree")

    c(.urp_control(criterion = ifelse(testtype == "Teststatistic", 
                                      "statistic", "p.value"),
                   logmincriterion = logmincriterion, minsplit = minsplit, 
                   minbucket = minbucket, minprob = minprob, stump = stump, 
                   mtry = mtry, maxdepth = maxdepth, multiway = multiway, 
                   splittry = splittry, MIA = MIA, maxsurrogate = maxsurrogate, 
                   majority = majority, caseweights = caseweights, 
                   applyfun = applyfun),
      list(teststat = teststat, splitstat = splitstat, splittest = splittest, pargs = pargs,
           testtype = testtype, nmax = nmax, nresample = nresample, lookahead = lookahead,
           intersplit = intersplit))
}

ctree <- function
(
    formula, 
    data, 
    weights, 
    subset,
    offset,
    cluster, 
    na.action = na.pass, 
    control = ctree_control(...), 
    ytrafo = NULL, 
    converged = NULL,
    scores = NULL,
    ...
) {

    ### get the call and the calling environment for .urp_tree
    call <- match.call(expand.dots = FALSE)
    call$na.action <- na.action
    frame <- parent.frame()
    if (missing(data)) {
        data <- NULL
        data_asis <- FALSE
    } else {
        data_asis <- missing(weights) && missing(subset) && 
                     missing(cluster) && missing(offset)
    }

    ### <FIXME> should be xtrafo
    if (!is.null(scores)) {
        if (missing(data))
            stop("can deal with scores with data being missing")
        for (n in names(scores)) {
            sc <- scores[[n]]
            if (is.ordered(data[[n]]) &&
                nlevels(data[[n]]) == length(sc)) {
                attr(data[[n]], "scores") <- as.numeric(sc)
            } else {
                warning("scores for variable ", sQuote(n), " ignored")
            }
        }
    }
    #### </FIXME>

    trafofun <- function(...) .ctreetrafo(..., ytrafo = ytrafo, converged = converged)
    tree <- .urp_tree(call, frame, data = data, data_asis = data_asis, control = control,
                      growfun = .ctreegrow, trafofun = trafofun,
                      doFit = TRUE)
    mf <- tree$mf
    weights <- model.weights(mf)
    if (is.null(weights)) weights <- rep(1, nrow(mf))

    fitted <- data.frame("(fitted)" = fitted_node(tree$node, mf), 
                         "(weights)" = weights,
                         check.names = FALSE)
    mf2 <- model.frame(Formula(formula), data = mf, na.action = na.pass)
    y <- model.part(Formula(formula), data = mf2, 
                    lhs = 1, rhs = 0)
    if (length(y) == 1) y <- y[[1]]
    fitted[[3]] <- y
    names(fitted)[3] <- "(response)"
    ret <- party(tree$node, data = mf, fitted = fitted, 
                 info = list(call = match.call(), control = control))
    ret$update <- tree$treefun
    ret$trafo <- tree$trafo
    class(ret) <- c("constparty", class(ret))

    ### doesn't work for Surv objects
    # ret$terms <- terms(formula, data = mf)
    ret$terms <- tree$terms
    ### need to adjust print and plot methods
    ### for multivariate responses
    ### if (length(response) > 1) class(ret) <- "party"
    return(ret)
}

#.logrank_trafo <- function(x, ties.method = c("logrank", "HL")) {
#    ties.method <- match.arg(ties.method)
#    time <- x[,1]
#    event <- x[,2]
#    n <- length(time)
#    ot <- order(time, event)
#    rt <- rank(time, ties.method = "max")
#    mt <- rank(time, ties.method = "min") - 1
#    fact <- switch(ties.method, "logrank" = event / (n - mt),
#                                "HL" = event/(n - rt + 1)
#                  )   
#    event - cumsum(fact[ot])[rt]
#}

.logrank_trafo <- function(...)
    return(coin::logrank_trafo(...))

### convert response y to influence function h(y)
.y2infl <- function(data, response, ytrafo = NULL) {

    if (length(response) == 1) {
        if (!is.null(ytrafo[[response]])) {
            yfun <- ytrafo[[response]]
            rtype <- "user-defined"
        } else {
            rtype <- class(data[[response]])[1]
            if (rtype == "integer") rtype <- "numeric"
        }
        response <- data[[response]]

        infl <- switch(rtype,
            "user-defined" = yfun(response),
            "factor" = { 
                X <- model.matrix(~ response - 1)
                if (nlevels(response) > 2) return(X)
                return(X[,-1, drop = FALSE])
            },
            "ordered" = {
                sc <- attr(response, "scores")
                if (is.null(sc)) sc <- 1L:nlevels(response)
                sc <- as.numeric(sc)
                return(matrix(sc[as.integer(response)], ncol = 1))
            },
            "numeric" = response,
            "Surv" = .logrank_trafo(response)
        )
    } else {
        ### multivariate response
        infl <- lapply(response, .y2infl, data = data)
        tmp <- do.call("cbind", infl)
        attr(tmp, "assign") <- rep(1L:length(infl), sapply(infl, NCOL))
        infl <- tmp
    }
    if (!is.matrix(infl)) infl <- matrix(infl, ncol = 1)
    storage.mode(infl) <- "double"
    return(infl)
}

sctest.constparty <- function(object, node = NULL, ...)
{

    if(is.null(node)) {
        ids <- nodeids(object, terminal = FALSE) ### all nodes
    } else {
        ids <- node
    }

    rval <- nodeapply(object, ids, function(n) {
        crit <- info_node(n)$criterion
        if (is.null(crit)) return(NULL)
        ret <- crit[c("statistic", "p.value"),,drop = FALSE]
        ret["statistic",] <- exp(ret["statistic",])
        ret <- ret[,-1,drop = FALSE] ### remove NA for response
        ret
    })
    names(rval) <- ids
    if(length(ids) == 1L)
        return(rval[[1L]])
    return(rval)
}
