
.ctree_test <- function(model, trafo, data, subset, weights, j, SPLITONLY = FALSE, ctrl) {

    ix <- index(data, j)
    iy <- model$index
    Y <- estfun(model)

    if (!is.null(ix) && !is.null(iy))
        return(.ctree_test_2d(data = data, j = j, Y = Y, iy = iy, 
                              subset = subset, weights = weights, 
                              SPLITONLY = SPLITONLY, ctrl = ctrl))

    NAy <- missings(data, "response")
    NAx <- missings(data, j)
    if (ctrl$MIA) {
        subsetNArm <- subset[!(subset %in% NAy)]
    } else {
        subsetNArm <- subset[!(subset %in% c(NAy, NAx))]
    }

    if (is.null(ix)) {
        if (!is.null(iy)) {
            iy <- NULL
            Y <- Y[index + 1L,, drop = FALSE]
        }
        return(.ctree_test_1d(data = data, j = j, Y = Y, subset = subsetNArm, 
                              weights = weights, SPLITONLY = SPLITONLY, ctrl = ctrl))
    } else {
        if (is.null(iy))
            return(.ctree_test_1d(data = data, j = j, Y = Y, subset = subsetNArm, 
                                  weights = weights, SPLITONLY = SPLITONLY, ctrl = ctrl))
    }
}

.ctree_split <- function(model, trafo, data, subset, weights, j, SPLITONLY = TRUE, ctrl)
    .ctree_test(model, trafo, data, subset, weights, j, SPLITONLY, ctrl)

.partysplit <- function(varid, breaks = NULL, index = NULL, right = TRUE, 
                        prob = NULL, info = NULL) {
    ret <- list(varid = varid, breaks = breaks, index = index, right = right, 
                prob = prob, info = info)
    class(ret) <- "partysplit"
    ret
}

.ctree_test_1d <- function(data, j, Y, subset, weights, SPLITONLY = FALSE, ctrl) {

    x <- model.frame(data)[[j]]
    MIA <- FALSE
    if (ctrl$MIA) MIA <- any(is.na(x[subset]))

    ### X for (ordered) factors is always dummy matrix
    if (is.factor(x) || is.ordered(x))
        X <- index(data, j)

    scores <- scores(data, j)
    ORDERED <- is.ordered(x) || is.numeric(x)

    ux <- Xleft <- Xright <- NULL
    
    if (ctrl$splittest || SPLITONLY) {
        MAXSELECT <- TRUE
        if (is.numeric(x)) {
            X <- index(data, j)
            ux <- levels(X)
        }
        if (MIA) {
            Xlev <- attr(X, "levels")
            Xleft <- X + 1L
            Xleft[is.na(Xleft)] <- 1L
            Xright <- X
            Xright[is.na(Xright)] <- as.integer(length(Xlev) + 1L)
            attr(Xleft, "levels") <- c(NA, Xlev)
            attr(Xright, "levels") <- c(Xlev, NA)
        } 
    } else {
        MAXSELECT <- FALSE
        if (is.numeric(x))
            X <- matrix(as.double(x), ncol = 1)
        MIA <- FALSE
    }
    cluster <- model.frame(data)[["(cluster)"]]

    .ctree_test_internal(x = x, X = X, ix = NULL, Xleft = Xleft, Xright = Xright, 
                         ixleft = NULL, ixright = NULL, ux = ux, scores = scores, 
                         j = j, Y = Y, iy = NULL, subset = subset, weights = weights, 
                         cluster = cluster, MIA = MIA, SPLITONLY = SPLITONLY, 
                         MAXSELECT = MAXSELECT, ORDERED = ORDERED, ctrl = ctrl)
}


.ctree_test_2d <- function(data, Y, iy, j, subset, weights, SPLITONLY = FALSE, ctrl) {

    x <- model.frame(data)[[j]]
    ix <- index(data, j)
    X <- ux <- attr(ix, "levels")
    MIA <- FALSE
    if (ctrl$MIA) MIA <- any(ix[subset] == 0)

    ### X for (ordered) factors is always dummy matrix
    if (is.factor(x) || is.ordered(x))
        X <- integer(0)

    scores <- scores(data, j)
    ORDERED <- is.ordered(x) || is.numeric(x)

    if (ctrl$splittest || SPLITONLY) {
        MAXSELECT <- TRUE
        X <- integer(0)

        if (is.numeric(ux)) {
            X <- index(data, j)
            ux <- levels(X)
        }
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
        MAXSELECT <- FALSE
        MIA <- FALSE
    }
    cluster <- model.frame(data)[["(cluster)"]]

    .ctree_test_internal(x = x, X = X, ix = ix, Xleft = Xleft, Xright = Xright, 
                         ixleft = ixleft, ixright = ixright, ux = ux, scores = scores, 
                         j = j, Y = Y, iy = iy, subset = subset, weights = weights, 
                         cluster = cluster, MIA = MIA, SPLITONLY = SPLITONLY, 
                         MAXSELECT = MAXSELECT, ORDERED = ORDERED, ctrl = ctrl)
}



.ctree_test_internal <- function(x, X, ix, Xleft, Xright, ixleft, ixright, ux, scores, j, Y , iy, 
                                 subset, weights, cluster, MIA, SPLITONLY, MAXSELECT, ORDERED, ctrl) {

    if (SPLITONLY) {
        nperm <- 0L
        varonly <- TRUE
        pvalue <- FALSE
        teststat <- ctrl$splitstat
    } else {
        nperm <- ifelse("MonteCarlo" %in% ctrl$testtype,
                        ctrl$nresample, 0L)
        if (ctrl$splittest) {
            if (ctrl$teststat != ctrl$splitstat)
                warning("Using different test statistics for testing and splitting")
            teststat <- ctrl$splitstat
            if (nperm == 0) 
               stop("MonteCarlo approximation mandatory for splittest = TRUE")
        } else {
           teststat <- ctrl$teststat
        }
        varonly <- "MonteCarlo" %in% ctrl$testtype && 
                   teststat == "maxtype"
        pvalue <- !("Teststatistic" %in% ctrl$testtype) 
    }

    ### see libcoin
    if (!is.null(cluster)) varonly <- FALSE 

    ### if (MIA) use tst as fallback
    ### compute linear statistic + expecation and covariance
    lev <- LinStatExpCov(X = X, Y = Y, subset = subset,
                         weights = weights, block = cluster,
                         nperm = nperm, varonly = varonly)
    if (is.ordered(x) && !ctrl$splittest) 
        lev <- matrix(scores, nrow = 1) %*% lev

    if (varonly) {
        vars <- lev$Variance
    } else {
        vars <- diag(vcov(lev))
    }
    if (all(vars < ctrl$tol)) return(list(statistic = NA, p.value = NA))

    ### compute test statistic and log(1 - p-value)
    tst <- doTest(lev, teststat = teststat, pvalue = pvalue,
                  lower = TRUE, log = TRUE, ordered = ORDERED, 
                  maxselect = MAXSELECT,
                  minbucket = ctrl$minbucket, pargs = ctrl$pargs)

    if (MIA) {
        ### compute linear statistic + expecation and covariance
        lev <- LinStatExpCov(X = Xleft, Y = Y, subset = subset,
                             weights = weights, block = cluster,
                             nperm = nperm, varonly = varonly)
        ### compute test statistic and log(1 - p-value)
        tstleft <- doTest(lev, teststat = teststat, pvalue = pvalue,
                          lower = TRUE, log = TRUE, ordered = ORDERED, 
                          minbucket = minbucket, pargs = ctrl$pargs)
        ### compute linear statistic + expecation and covariance
        lev <- LinStatExpCov(X = Xright, Y = Y, subset = subset,
                             weights = weights, block = cluster,
                             nperm = nperm, varonly = varonly)
        ### compute test statistic and log(1 - p-value)
        tstright <- doTest(lev, teststat = teststat, pvalue = pvalue,
                           lower = TRUE, log = TRUE, ordered = ORDERED, 
                           minbucket = minbucket, pargs = ctrl$pargs)
    }

    if (!SPLITONLY) {
        if (MIA) {
            tst <- tstleft
            if (tst$TestStatistic < tstright$TestStatistic)
                tst <- tstright
        }        
        return(list(statistic = log(pmax(tst$TestStatistic, .Machine$double.eps)), 
                    p.value = tst$p.value))
    }

    ret <- NULL
    if (MIA && !any(is.na(tst$index))) {
        if (ORDERED) {
            if (tstleft$TestStatistic >= tstright$TestStatistic) {
                if (all(tst$index == 1)) { ### case C
                    ret <- .partysplit(as.integer(j), breaks = Inf, 
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
                    ret <- .partysplit(as.integer(j), breaks = sp,
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
                ret <- .partysplit(as.integer(j), breaks = sp,
                                  index = 1L:2L, prob = as.double(0:1))
            }
        } else {
            sp <- tstleft$index[-1L] ### tstleft = tstright for unordered factors
            if (length(unique(sp)) == 1L) { ### case C
                ret <- .partysplit(as.integer(j), index = as.integer(tst$index) + 1L)
            } else { ### always case A
                ret <- .partysplit(as.integer(j),
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
                ret <- .partysplit(as.integer(j), breaks = sp,
                                  index = 1L:2L)
        } else {
            ret <- .partysplit(as.integer(j),
                              index = as.integer(sp) + 1L)
        }
    }
    return(ret)

}


library("partykit")

model.frame.partydata <- function(object)
    object

scores <- function(object, varid) {
    x <- model.frame(object)[[varid]]
    if (is.ordered(x)) return(1:nlevels(x))
    return(NULL)
}

index <- function(object, varid) {
    x <- model.frame(object)[[varid]]
    if (is.factor(x) || is.ordered(x))
        return(x)
    ux <- sort(unique(x))
    X <- .bincode(x, breaks = c(-Inf, ux, Inf),
                  right = TRUE)
    attr(X, "levels") <- ux 
    storage.mode(X) <- "integer"
    X
}

missings <- function(object, varid) {
    if (is.character(varid)) varid <- 1L
    which(is.na(model.frame(object)[[varid]]))
}

d <- iris[,c(5, 1:4)]
class(d) <- c("partydata", "data.frame")
Y <- model.matrix(~ Species - 1, data = iris)
subset <- integer(0)
weights <- integer(0)
control <- ctree_control()

.ctree_test_1d(data = d, j = 1, Y = Y, subset = subset, weights = weights, SPLITONLY = FALSE, ctrl = control)

.ctree_test_1d(data = d, j = 1, Y = Y, subset = subset, weights = weights, SPLITONLY = TRUE, ctrl = control)

(ct <- ctree(Species ~ Sepal.Length, data = iris, stump = TRUE))
info_node(node_party(ct))

source("extree.R")

control$update <- FALSE

estfun <- function(object) object$estfun

.extree_fit(data = d, trafo = function(subset, weights, ...) list(estfun = Y), partyvars = 2:5,
subset = 1:nrow(d), integer(0), ctrl = control)

ctree(Species ~ ., data = d)
