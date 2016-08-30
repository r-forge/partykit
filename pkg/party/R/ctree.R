
.ctree_test_split <- function(x, bdr = NULL, j, ctrl, X, Y, iy = NULL, subset, 
                              weights, cluster, splitonly = TRUE, minbucket) {

    if (splitonly) {
        if ((ctrl$multiway && ctrl$maxsurrogate == 0) &&
            is.factor(x)) {
            index <- 1L:nlevels(x)
            if (length(weights) > 0) {
                xt <- xtabs(~ x, subset = subset)
            } else {
                xt <- xtabs(weights ~ x, subset = subset)
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
        } 
    } else {
         ix <- bdr[[j]]
         if (ctrl$splittest || splitonly) {
             X <- numeric(0) 
             ux <- attr(ix, "levels")
         } 
    }

    if (splitonly) {
        B <- 0L
        varonly <- TRUE
        pvalue <- FALSE
        teststat <- ctrl$splitstat
    } else {
        if (ctrl$splittest) {
            if (ctrl$teststat != ctrl$splitstat)
                warning("Using different test statistics for testing and splitting")
            teststat <- ctrl$splitstat
        } else {
            teststat <- ctrl$teststat
        }
        B <- ifelse(ctrl$testtype == "MonteCarlo",
                    ctrl$nresample, 0L)
        varonly <- ctrl$testtype == "MonteCarlo" && 
                   teststat == "maxtype"
        pvalue <- ctrl$testtype != "Teststatistic"
    }

    ### compute linear statistic + expecation and covariance
    lev <- LinStatExpCov(X = X, Y = Y, ix = ix, iy = iy, subset = subset,
                         weights = weights, block = cluster,
                         B = B, varonly = varonly)
    ### compute test statistic and log(1 - p-value)
    tst <- doTest(lev, teststat = teststat, pvalue = pvalue,
                  lower = TRUE, log = TRUE, ordered = ORDERED, minbucket = minbucket)

    if (splitonly) {
         ret <- NULL
         sp <- tst$index
         if (!all(is.na(sp))) {
             if (length(sp) == 1) {
                 if (!is.ordered(x))
                     sp <- ux[sp]
                 ret <- partysplit(as.integer(j), breaks = sp,
                                   index = 1L:2L)
             } else {
                  ret <- partysplit(as.integer(j),
                                    index = as.integer(sp) + 1L)
             }
         }
         return(ret)
    }  

    return(list(statistic = log(tst$TestStatistic), p.value = tst$p.value))
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
        y <- model.part(f, data = mf, lhs = 1, rhs = 0)
        bdr <- BDR::BDR(y, nmax = ctrl$nmax, total = TRUE, 
                        complete.cases.only = TRUE)
        y <- attr(bdr, "levels")
        index <- c(bdr)
        attr(index, "levels") <- 1:NROW(y)
        cn <- colnames(y)
        Y <- partykit:::.y2infl(y, cn[cn != "(weights)"], ytrafo = ytrafo)
        ### first row corresponds to missings
        Y <- rbind(0, Y)  
        return(function(subset)
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
        y <- model.part(f, data = mf, lhs = 1, rhs = 0)
        cc <- complete.cases(y)
        Yi <- partykit:::.y2infl(y[cc,,drop = FALSE], colnames(y), ytrafo = ytrafo)
        Y <- matrix(NA, nrow = nrow(mf), ncol = NCOL(Yi))
        Y[cc,] <- Yi
        #    colnames(Y) <- colnames(Yi)
        storage.mode(Y) <- "double"
        return(function(subset)
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
                              weights = integer(0), whichvar) 
        {

            ret <- list(criteria = matrix(NA, nrow = 2L, ncol = ncol(data)))
            colnames(ret$criteria) <- names(data)
            rownames(ret$criteria) <- c("statistic", "p.value")

            if (is.null(y)) {
                ### nrow(Y) = nrow(data)!!!
                tr <- trafo(subset)
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

    bdr <- BDR::BDR(data, nmax = ctrl$nmax)
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
                              weights = integer(0), whichvar) 
        {
    
            ret <- list(criteria = matrix(NA, nrow = 2L, ncol = ncol(data)))
            colnames(ret$criteria) <- names(data)
            rownames(ret$criteria) <- c("statistic", "p.value")

            if (is.null(y)) {
                tr <- trafo(subset = subset)
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
    splitstat = c("quadratic", "maximum"), ### much better for q > 1
    splittest = FALSE,
    testtype = c("Bonferroni", "MonteCarlo", 
                 "Univariate", "Teststatistic"),
    nmax = Inf, 
    alpha = 0.05, 
    mincriterion = 1 - alpha, 
    logmincriterion = log(mincriterion), 
    minsplit = 20L, 
    minbucket = 7L, 
    minprob = 0.01, 
    stump = FALSE, 
    nresample = 9999L, 
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
                   splittry = splittry, maxsurrogate = maxsurrogate, 
                   majority = majority, caseweights = caseweights, 
                   applyfun = applyfun),
      list(teststat = teststat, splitstat = splitstat, splittest = splittest,
           testtype = testtype, nmax = nmax, nresample = nresample))
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
