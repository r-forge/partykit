
.ctree_test <- function(model, trafo, data, subset, weights, j, SPLITONLY = FALSE, ctrl) {

    ix <- .get_index(data, j)
    iy <- .get_index(data, "yx")
    Y <- model$estfun

    if (!is.null(iy))
        return(.ctree_test_2d(data = data, j = j, Y = Y, iy = iy, 
                              subset = subset, weights = weights, 
                              SPLITONLY = SPLITONLY, ctrl = ctrl))

    NAyx <- .get_NAs(data, "yx")
    NAz <- .get_NAs(data, j)
    if (ctrl$MIA) {
        subsetNArm <- subset[!(subset %in% NAyx)]
    } else {
        subsetNArm <- subset[!(subset %in% c(NAyx, NAz))]
    }

    return(.ctree_test_1d(data = data, j = j, Y = Y, subset = subsetNArm, 
                          weights = weights, SPLITONLY = SPLITONLY, ctrl = ctrl))
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

    x <- .get_var(data, j)
    MIA <- FALSE
    if (ctrl$MIA) MIA <- (length(NAs <- .get_NAs(data, j)) > 0)

    ### X for (ordered) factors is always dummy matrix
    if (is.factor(x) || is.ordered(x))
        X <- .get_index(data, j)

    scores <- .get_scores(data, j)
    ORDERED <- is.ordered(x) || is.numeric(x)

    ux <- Xleft <- Xright <- NULL
    
    if (ctrl$splittest || SPLITONLY) {
        MAXSELECT <- TRUE
        if (is.numeric(x)) {
            X <- .get_index(data, j)
            ux <- levels(X)
        }
        if (MIA) {
            Xlev <- attr(X, "levels")
            Xleft <- X + 1L
            Xleft[NAs] <- 1L
            Xright <- X
            Xright[NAs] <- as.integer(length(Xlev) + 1L)
            attr(Xleft, "levels") <- c(NA, Xlev)
            attr(Xright, "levels") <- c(Xlev, NA)
        } 
    } else {
        MAXSELECT <- FALSE
        if (is.numeric(x)) {
            if (storage.mode(x) == "double") {
                X <- x
            } else {
                X <- as.double(x) ### copy when necessary
            }
        }
        MIA <- FALSE
    }
    cluster <- .get_var(data, "(cluster)")

    .ctree_test_internal(x = x, X = X, ix = NULL, Xleft = Xleft, Xright = Xright, 
                         ixleft = NULL, ixright = NULL, ux = ux, scores = scores, 
                         j = j, Y = Y, iy = NULL, subset = subset, weights = weights, 
                         cluster = cluster, MIA = MIA, SPLITONLY = SPLITONLY, 
                         MAXSELECT = MAXSELECT, ORDERED = ORDERED, ctrl = ctrl)
}


.ctree_test_2d <- function(data, Y, iy, j, subset, weights, SPLITONLY = FALSE, ctrl) {

    x <- .get_var(data, j)
    ix <- .get_index(data, j)
    X <- ux <- attr(ix, "levels")
    MIA <- FALSE
    if (ctrl$MIA) MIA <- any(ix[subset] == 0)

    ### X for (ordered) factors is always dummy matrix
    if (is.factor(x) || is.ordered(x))
        X <- integer(0)

    scores <- .get_scores(data, j)
    ORDERED <- is.ordered(x) || is.numeric(x)

    if (ctrl$splittest || SPLITONLY) {
        MAXSELECT <- TRUE
        X <- integer(0)

        if (is.numeric(ux)) {
            X <- .get_index(data, j)
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
    cluster <- .get_var(data, "(cluster)")

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
    lev <- LinStatExpCov(X = X, Y = Y, ix = ix, iy = iy, subset = subset,
                         weights = weights, block = cluster,
                         nperm = nperm, varonly = varonly, checkNAs = FALSE)
    if (is.ordered(x) && !ctrl$splittest) 
        lev <- matrix(scores, nrow = 1) %*% lev

    ### check if either X or Y were unique
    if (all(lev$Variance < ctrl$tol)) 
        return(list(statistic = NA, p.value = NA))

    ### compute test statistic and log(1 - p-value)
    tst <- doTest(lev, teststat = teststat, pvalue = pvalue,
                  lower = TRUE, log = TRUE, ordered = ORDERED, 
                  maxselect = MAXSELECT,
                  minbucket = ctrl$minbucket, pargs = ctrl$pargs)

    if (MIA) {
        ### compute linear statistic + expecation and covariance
        lev <- LinStatExpCov(X = Xleft, Y = Y, ix = ixleft, iy = iy, subset = subset,
                             weights = weights, block = cluster,
                             nperm = nperm, varonly = varonly, checkNAs = FALSE)
        ### compute test statistic and log(1 - p-value)
        tstleft <- doTest(lev, teststat = teststat, pvalue = pvalue,
                          lower = TRUE, log = TRUE, ordered = ORDERED, 
                          minbucket = minbucket, pargs = ctrl$pargs)
        ### compute linear statistic + expecation and covariance
        lev <- LinStatExpCov(X = Xright, Y = Y, ix = ixright, iy = iy, subset = subset,
                             weights = weights, block = cluster,
                             nperm = nperm, varonly = varonly, checkNAs = FALSE)
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

model.frame.extree_data <- function(object, yxonly = FALSE, ...) {
    if (!yxonly) 
        return(object$data)
    if (!is.null(object$yxindex))
        return(attr(object$yxindex, "levels"))
    vars <- object$variables
    return(object$data[, c(vars$y, vars$x, vars$offset),drop = FALSE])
}    

.get_var <- function(object, varid) {
    mf <- model.frame(object)
    ### [[.data.frame needs lots of memory
    class(mf) <- "list"
    mf[[varid]]
}

.get_scores <- function(object, varid) {
    x <- .get_var(object, varid)
    if (is.ordered(x)) {
        sc <- object$scores[[varid]]
        if (is.null(sc)) sc <- 1:nlevels(x)
        return(sc)
    }
    return(NULL)
}

.get_index <- function(object, varid) {
    if (varid %in% c(object$variables$y, object$variables$x))
        return(object$yxindex) ### may be NULL

    return(object$zindex[[varid]])
}

.get_NAs <- function(object, varid) {
    if (is.character(varid)) {
        if (varid == "yx")
            return(object$yxmissings)
        stop("varid not found")
    }
    object$missings[[varid]]
}

Ctree <- function(formula, data, subset, na.action = na.pass, weights, offset, cluster, control, ...) {

    ## set up model.frame() call
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action", "weights", "offset", "cluster"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$yx <- "matrix"
    mf$nmax <- control$nmax
    ## evaluate model.frame
    mf[[1L]] <- quote(extree)

    d <- eval(mf, parent.frame())

    Y <- d$yx$y
    subset <- 1:nrow(model.frame(d))
    control$update <- FALSE

    .extree_fit(data = d, trafo = function(subset, weights, ...) list(estfun = Y), 
                partyvars = d$variables$z, 
                subset = subset, .get_var(d, "(weights)"), ctrl = control)

}

library("partykit")

source("extree.R")

iris$block <- sample(gl(3, 1), nrow(iris), replace = TRUE)
iris$x <- sample(gl(3, 1, ordered = TRUE), nrow(iris), replace = TRUE)
iris$Species[1] <- NA
iris$Petal.Length[2] <- NA

a1 <- extree(Species ~ Petal.Length | Sepal.Length + Sepal.Width + x, data = iris, na.action = na.pass, 
       offset = Petal.Width, cluster = block, nmax = Inf, yx = "none")

object.size(iris)
object.size(a1)

head(model.frame(a1))
head(model.frame(a1, yxonly = TRUE))
.get_var(a1, varid = 1)
.get_scores(a1, varid = 5)
.get_index(a1, varid = 2)
.get_NAs(a1, varid = 2)
.get_NAs(a1, varid = "yx")

a2 <- extree(Species ~ Petal.Length | Sepal.Length + Sepal.Width + x, data = iris, na.action = na.pass, 
       offset = Petal.Width, cluster = block, nmax = Inf, yx = "matrix")

head(model.frame(a2))
head(model.frame(a2, yxonly = TRUE))
.get_var(a2, varid = 1)
.get_scores(a2, varid = 5)
.get_index(a2, varid = 2)
.get_NAs(a2, varid = 2)
.get_NAs(a2, varid = "yx")


a3 <- extree(Species ~ Petal.Length | Sepal.Length + Sepal.Width + x, data = iris, na.action = na.pass, 
       offset = Petal.Width, cluster = block, nmax = 5, yx = "none")

head(model.frame(a3))
head(model.frame(a3, yxonly = TRUE))
.get_var(a3, varid = 1)
.get_scores(a3, varid = 5)
.get_index(a3, varid = 2)
.get_NAs(a3, varid = 2)
.get_NAs(a3, varid = "yx")


a4 <- extree(Species ~ Petal.Length | Sepal.Length + Sepal.Width + x, data = iris, na.action = na.pass, 
       offset = Petal.Width, cluster = block, nmax = 5, yx = "matrix")

head(model.frame(a1))
head(model.frame(a4, yxonly = TRUE))
.get_var(a4, varid = 1)
.get_scores(a4, varid = 5)
.get_index(a4, varid = 2)
.get_NAs(a4, varid = 2)
.get_NAs(a4, varid = "yx")

data(iris)
iris <- iris[sample(1:nrow(iris), 10000, replace = TRUE),,drop = FALSE]
dim(iris)

system.time(ct <- ctree(Species ~ ., data = iris))
info_node(node_party(ct))

system.time(ct2 <- Ctree(Species ~ ., data = iris, control = ctree_control()))
info_node(ct2)

