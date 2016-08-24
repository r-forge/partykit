
ctreegrow <- function(data, partyvars, block, ctrl) {
    if (ctrl$nmax < Inf)
        return(.ctree_fit_2d(data = data, partyvars = partyvars, 
                             block = block, ctrl = ctrl))
    return(.ctree_fit_1d(data = data, partyvars = partyvars,
                         block = block, ctrl = ctrl))
}

ctreetrafo <- function(formula, data, weights, block, ctrl, ytrafo) {
    if (ctrl$nmax < Inf) {
        if (!is.function(ytrafo))
            return(.cfit_2d(formula, data = data, weights = weights, block = block,
                            nmax = ctrl$nmax, ytrafo = ytrafo))
        return(ytrafo(formula, data = data, weights = weights, block = block,
                      nmax = ctrl$nmax))
    } else {
        if (!is.function(ytrafo))
            return(.cfit(formula, data = data, weights = weights, block = block,
                         ytrafo = ytrafo))
        return(ytrafo(formula, data = data, weights = weights, block = block))
    }
}
    

.cfit <- function(formula, data, weights = NULL, block = NULL, ytrafo = NULL) {
    f <- Formula(formula)
    mf <- model.frame(formula = f, data = data)
    y <- model.part(f, data = mf, lhs = 1, rhs = 0)
    Y <- partykit:::.y2infl(y, colnames(y), ytrafo = ytrafo)
    function(subset) {
        list(estfun = Y)
    }
}

.cfit_2d <- function(formula, data, weights = NULL, block = NULL,
                   nmax = 25, ytrafo = ytrafo) {
    f <- Formula(formula)
    mf <- model.frame(formula = f, data = data)
    y <- model.part(f, data = mf, lhs = 1, rhs = 0)
    bdr <- BDR::BDR(y, nmax = nmax, total = TRUE, complete.cases.only = TRUE)
    y <- attr(bdr, "levels")
    index <- c(bdr)
    attr(index, "levels") <- 1:NROW(y)
    cn <- colnames(y)
    Y <- partykit:::.y2infl(y, cn[cn != "(weights)"], ytrafo = ytrafo)
    ### first row corresponds to missings
    Y <- rbind(0, Y)
    function(subset) {
        list(estfun = Y, index = index)
    }
}

### conditional inference trees
.ctree_fit_1d <- function
(
    data, 				### full data, readonly
    partyvars, 				### partytioning variables,
					### a subset of 1:ncol(data)
    block = integer(0), 		### a blocking factor w/o NA
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

            ret <- list(p = matrix(NA, nrow = 2L, ncol = ncol(data)))
            colnames(ret$p) <- names(data)
            rownames(ret$p) <- c("statistic", "p.value")

            if (is.null(y)) {
                ### nrow(Y) = nrow(data)!!!
                Y <- trafo(subset)$estfun
            } else {
                ### y is kidids in .csurr and nothing else
                stopifnot(length(y) == length(subset))
                Y <- matrix(0, nrow = NROW(data), ncol = max(y))
                Y[cbind(subset, y)] <- 1 ### model.matrix(~ as.factor(y) - 1)
                storage.mode(Y) <- "double"
            }

            for (j in whichvar) {
                ### compute linear statistic + expecation and covariance
                lev <- LinStatExpCov(X = X[[j]], Y = Y, subset = subset,
                                     weights = weights, block = block, 
                                     B = ifelse(ctrl$testtype == "MonteCarlo", 
                                                ctrl$nresample, 0L))
                ### compute test statistic and log(1 - p-value)
                tst <- doTest(lev, teststat = ctrl$teststat, 
                              pvalue = ctrl$testtype != "Teststatistic",
                              lower = TRUE, log = TRUE)
                ret$p["statistic", j] <- log(tst$TestStatistic)
                ret$p["p.value", j] <- tst$p.value
            }
            if (ctrl$testtype == "Bonferroni")
                ret$p["p.value",] <- ret$p["p.value",] * length(whichvar)

            ### compute best fitting cutpoint (according to minbucket)
            ### for a subset of observations and variables
            ### y is used for node ids when computing surrogate splits
            ### splitfun as part of the return object of selectfun allows
            ### returning splits found during selections (ie, exhaustive searchs)
            ### splitfun already knows Y by lexical scoping, no need to compute it twice
            ret$splitfun <- function(whichvar, minbucket) 
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
                        index[xt < minbucket] <- nlevels(x) + 1L
                        index <- unclass(factor(index))
                        ret <- partysplit(as.integer(j), index = as.integer(index))
                        break()
                    } else {
                        ### X being an integer triggers maximally selected stats
                        if (is.factor(x)) {
                            X <- unclass(x)
                        } else {
                            x[-subset] <- NA
                            ux <- sort(unique(x))
                            X <- cut.default(x, breaks = c(-Inf, ux, Inf),
                                             labels = FALSE, right = TRUE)
                            # X[is.na(X)] <- 0L (NAs are handled by LinStatExpCov)
                            attr(X, "levels") <- ux
                            storage.mode(X) <- "integer"
                        }
                        lev <- LinStatExpCov(X = X, Y = Y, subset = subset,
                                             weights = weights, block = block, 
                                             B = 0L, varonly = TRUE)
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
           ret
        }

        tree <- .urp_node(id = 1L, data = data, 
                          selectfun = function(...) selectfun(..., trafo = trafo),
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
    block = integer(0),                 ### a blocking factor w/o NA
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
    
            ret <- list(p = matrix(NA, nrow = 2L, ncol = ncol(data)))
            colnames(ret$p) <- names(data)
            rownames(ret$p) <- c("statistic", "p.value")

            if (is.null(y)) {
                tr <- trafo(subset = subset)
                Y <- tr$estfun
                iy <- tr$index
            } else {
                ### y is kidids in .csurr and nothing else
                stopifnot(length(y) == length(subset))
                Y <- rbind(0, max(y))
                iy <- numeric(NROW(data))
                iy[subset] <- as.integer(y)
            }

            for (j in whichvar) {
                ix <- bdr[[j]]
                lev <- LinStatExpCov(X = X[[j]], ix = ix, 
                                     Y = Y, iy = iy, subset = subset,
                                     weights = weights, block = block, 
                                     B = ifelse(ctrl$testtype == "MonteCarlo", 
                                                ctrl$nresample, 0L))
                tst <- doTest(lev, teststat = ctrl$teststat, 
                              pvalue = ctrl$testtype != "Teststatistic",
                              lower = TRUE, log = TRUE)
                ret$p["statistic", j] <- log(tst$TestStatistic)
                ret$p["p.value", j] <- tst$p.value
            }
            if (ctrl$testtype == "Bonferroni")
                ret$p["p.value",] <- ret$p["p.value",] * length(whichvar)


            ### compute best fitting cutpoint (according to minbucket)
            ### for a subset of observations and variables
            ### y is used for node ids when computing surrogate splits
            ret$splitfun <- function(whichvar,  minbucket)
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
                        index[xt < minbucket] <- nlevels(x) + 1L
                        index <- unclass(factor(index))
                        ret <- partysplit(as.integer(j), index = as.integer(index))
                        break()
                    } else {
                        ix <- bdr[[j]]
                        X <- numeric(0) 
                        ux <- attr(ix, "levels")
                        lev <- LinStatExpCov(X = X, ix = ix,
                                             Y = Y, iy = iy, subset = subset,
                                             weights = weights, block = block, 
                                             B = 0L, varonly = TRUE)
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
               return(ret)
           }
           ret
        }

        tree <- .urp_node(id = 1L, data = data, 
                          selectfun = function(...) selectfun(..., trafo = trafo),
                          partyvars = partyvars, weights = weights,
                          subset = subset, ctrl = ctrl)
        return(tree)
    })
}

ctree_control <- function(teststat = c("quadratic", "maximum"), 
                          splitstat = c("maximum", "quadratic"),
                          testtype = c("Bonferroni", "MonteCarlo", 
                                       "Univariate", "Teststatistic"),
                          nmax = Inf, 
                          alpha = 0.05, mincriterion = 1 - alpha, 
                          logmincriterion = log(mincriterion), minsplit = 20L, 
                          minbucket = 7L, minprob = 0.01, stump = FALSE, 
                          nresample = 9999L, maxsurrogate = 0L, mtry = Inf, 
                          maxdepth = Inf, multiway = FALSE, splittry = 2L, 
                          majority = FALSE, applyfun = NULL, cores = NULL) {

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

    if (multiway & maxsurrogate > 0)
        stop("surrogate splits currently not implemented for multiway splits")

    list(teststat = teststat, splitstat = splitstat, 
         criterion = ifelse(testtype == "Teststatistic", "statistic", "p.value"),
         testtype = testtype, nmax = nmax, logmincriterion = logmincriterion,
         minsplit = minsplit, minbucket = minbucket, 
         minprob = minprob, stump = stump, nresample = nresample, mtry = mtry, 
         maxdepth = maxdepth, multiway = multiway, splittry = splittry,
         maxsurrogate = maxsurrogate, majority = majority, 
         applyfun = applyfun)
}

ctree <- function(formula, data, weights, subset, na.action = na.pass, 
                  control = ctree_control(...), ytrafo = NULL, 
                  scores = NULL, ...) {

    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action"),
               names(mf), 0)
    mf <- mf[c(1, m)]

    f <- if (inherits(formula, "Formula")) formula else Formula(formula)
    if (length(length(f)) != 2)
        stop("incorrect formula")  
    if (!(length(f)[2] %in% 1:3))
        stop("incorrect formula")
    mf$formula <- f
    mf$drop.unused.levels <- FALSE      
    mf$na.action <- na.pass
    mf[[1]] <- quote(stats::model.frame)
    mf1 <- eval(mf, parent.frame())
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
         mf2 <- eval(mf, parent.frame())
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

    block <- NULL
    if (!is.null(blockf)) {
        block <- model.frame(blockf, data = mf)
        if (length(block) != 1 || !is.factor(block[[1]]))
            stop("block is not a single factor")
    }

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

    treefun <- ctreegrow(mf, partyvars = match(zvars, colnames(mf)), 
                         block = block, ctrl = control)
    ytrafo <- ctreetrafo(modelf, data = mf, weights = weights, block = block, 
                         ctrl = control, ytrafo = ytrafo)
    tree <- treefun(ytrafo, subset = 1:nrow(mf), weights)

    if (length(weights) == 0)
        weights <- rep.int(1, nrow(mf))
    fitted <- data.frame("(fitted)" = fitted_node(tree, mf), 
                         "(weights)" = weights,
                         check.names = FALSE)
    y <- model.part(Formula(modelf), data = mf, lhs = 1, rhs = 0)
    if (length(y) == 1) y <- y[[1]]
    fitted[[3]] <- y
    names(fitted)[3] <- "(response)"
    ret <- party(tree, data = mf, fitted = fitted, 
                 info = list(call = match.call(), control = control))
    ret$update <- treefun
    class(ret) <- c("constparty", class(ret))

    ### doesn't work for Surv objects
    # ret$terms <- terms(formula, data = mf)
    ret$terms <- mfterms
    ### need to adjust print and plot methods
    ### for multivariate responses
    ### if (length(response) > 1) class(ret) <- "party"
    return(ret)
}
