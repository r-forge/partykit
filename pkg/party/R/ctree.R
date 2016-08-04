
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
                Y <- trafo(data = data, subset = subset, 
                           weights = weights)
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

    bdr <- libcoin:::BDR(data, ctrl$nmax)
    X <- lapply(bdr, function(x) attr(x, "X"))

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
                tr <- trafo(data = data, subset = subset, 
                           weights = weights)
                Y <- attr(tr[[1]], "X")
                iy <- tr[[1]]
            } else {
                ### y is kidids in .csurr and nothing else
                stopifnot(length(y) == length(subset))
                Y <- rbind(0, max(y))
                iy <- numeric(NROW(data))
                iy[subset] <- as.integer(y)
            }

            for (j in whichvar) {
                lev <- LinStatExpCov(X = X[[j]], ix = bdr[[j]], 
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
                        ux <- attr(bdr[[j]], "levels")
                        lev <- LinStatExpCov(X = X, ix = bdr[[j]],
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
                          nmax = Inf, mincriterion = 0.95, minsplit = 20L, 
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
    if (is.null(weights)) weights <- integer(0)
    fitdat <- mf[, colnames(mf) != "(weights)"]

    ### <FIXME> should be xtrafo
    if (!is.null(scores)) {
        for (n in names(scores)) {
            sc <- scores[[n]]
            if (is.ordered(fitdat[[n]]) && 
                nlevels(fitdat[[n]]) == length(sc)) {
                attr(fitdat[[n]], "scores") <- as.numeric(sc)
            } else {
                warning("scores for variable ", sQuote(n), " ignored")
            }
        }
    }
    #### </FIXME>

    ### <FIXME> implement y ~ x | block or y ~ 1 | x | block ? </FIXME>
    block <- integer(0)

    if (!is.function(ytrafo)) {
        Y <- partykit:::.y2infl(fitdat, response, ytrafo = ytrafo)
        ytrafo <- function(...) Y
    }

    if (control$nmax < Inf) {
        treefun <- .ctree_fit_2d(fitdat, partyvars = which(!(colnames(fitdat) %in% response)), 
                          block = block, ctrl = control)
        Y <- libcoin:::BDR(fitdat[response], nmax = control$nmax)
        tree <- treefun(function(...) Y, subset = 1:nrow(fitdat), weights)
    } else {
        treefun <- .ctree_fit_1d(fitdat, partyvars = which(!(colnames(fitdat) %in% response)), 
                                 block = block, ctrl = control)
        tree <- treefun(ytrafo, subset = 1:nrow(fitdat), weights)
    }

    if (length(weights) == 0)
        weights <- rep.int(1, nrow(fitdat))
    fitted <- data.frame("(fitted)" = fitted_node(tree, fitdat), 
                         "(weights)" = weights,
                         check.names = FALSE)
    fitted[[3]] <- fitdat[, response, drop = length(response) == 1]
    names(fitted)[3] <- "(response)"
    ret <- party(tree, data = fitdat, fitted = fitted, 
                 info = list(call = match.call(), control = control))
    ret$update <- treefun
    class(ret) <- c("constparty", class(ret))

    ### doesn't work for Surv objects
    # ret$terms <- terms(formula, data = mf)
    ret$terms <- terms(mf)
    ### need to adjust print and plot methods
    ### for multivariate responses
    ### if (length(response) > 1) class(ret) <- "party"
    return(ret)
}
