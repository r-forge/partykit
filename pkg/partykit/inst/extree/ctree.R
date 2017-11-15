
.ctree_test <- function(model, trafo, data, subset, weights, j, SPLITONLY = FALSE, ctrl) {

    ix <- data[[j, type = "index"]]
    iy <- data[["yx", type = "index"]]
    Y <- model$estfun

    if (!is.null(iy)) {
        stopifnot(NROW(levels(iy)) == (NROW(Y) - 1))
        return(.ctree_test_2d(data = data, j = j, Y = Y, iy = iy, 
                              subset = subset, weights = weights, 
                              SPLITONLY = SPLITONLY, ctrl = ctrl))
    }

    stopifnot(NROW(Y) == length(ix))

    NAyx <- data[["yx", type = "missing"]]
    NAz <- data[[j, type = "missing"]]
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

    x <- data[[j]]
    MIA <- FALSE
    if (ctrl$MIA) MIA <- (length(NAs <- data[[j, type = "missing"]]) > 0)

    ### X for (ordered) factors is always dummy matrix
    if (is.factor(x) || is.ordered(x))
        X <- data[[j, type = "index"]]

    scores <- data[[j, type = "scores"]]
    ORDERED <- is.ordered(x) || is.numeric(x)

    ux <- Xleft <- Xright <- NULL
    
    if (ctrl$splittest || SPLITONLY) {
        MAXSELECT <- TRUE
        if (is.numeric(x)) {
            X <- data[[j, type = "index"]]
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
    cluster <- data[["(cluster)"]]

    .ctree_test_internal(x = x, X = X, ix = NULL, Xleft = Xleft, Xright = Xright, 
                         ixleft = NULL, ixright = NULL, ux = ux, scores = scores, 
                         j = j, Y = Y, iy = NULL, subset = subset, weights = weights, 
                         cluster = cluster, MIA = MIA, SPLITONLY = SPLITONLY, 
                         MAXSELECT = MAXSELECT, ORDERED = ORDERED, ctrl = ctrl)
}


.ctree_test_2d <- function(data, Y, iy, j, subset, weights, SPLITONLY = FALSE, ctrl) {

    x <- data[[j]]
    ix <- data[[j, type = "index"]]
    ux <- attr(ix, "levels")

    MIA <- FALSE
    if (ctrl$MIA) MIA <- any(ix[subset] == 0)

    ### X for (ordered) factors is always dummy matrix
    if (is.factor(x) || is.ordered(x))
        X <- integer(0)

    scores <- data[[j, type = "scores"]]
    ORDERED <- is.ordered(x) || is.numeric(x)

    if (ctrl$splittest || SPLITONLY) {
        MAXSELECT <- TRUE
        X <- integer(0)

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
        if (is.numeric(x))
            X <- matrix(c(0, as.double(attr(ix, "levels"))), ncol = 1)
    }
    cluster <- data[["(cluster)"]]

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

.y2infl <- partykit:::.y2infl

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
    MIA = FALSE,	### DOI: 10.1016/j.patrec.2008.01.010
    nresample = 9999L, 
    tol = sqrt(.Machine$double.eps),
    maxsurrogate = 0L, 
    numsurrogate = FALSE,
    mtry = Inf, 
    maxdepth = Inf, 
    multiway = FALSE, 
    splittry = 2L, 
    intersplit = FALSE,
    majority = FALSE, 
    caseweights = TRUE, 
    applyfun = NULL, 
    cores = NULL,
    saveinfo = TRUE,
    update = FALSE
) {

    testtype <- match.arg(testtype, several.ok = TRUE)
    if (length(testtype) == 4) testtype <- testtype[1]
    ttesttype <- testtype
    if (length(testtype) > 1) {
        stopifnot(all(testtype %in% c("Bonferroni", "MonteCarlo")))
        ttesttype <- "MonteCarlo"
    }

    splitstat <- match.arg(splitstat)
    teststat <- match.arg(teststat)

    if (!caseweights)
        stop("only caseweights currently implemented in ctree")

    c(extree_control(criterion = ifelse("Teststatistic" %in% testtype, 
                                      "statistic", "p.value"),
                     logmincriterion = logmincriterion, minsplit = minsplit, 
                     minbucket = minbucket, minprob = minprob, 
                     nmax = nmax, stump = stump, lookahead = lookahead,
                     mtry = mtry, maxdepth = maxdepth, multiway = multiway, 
                     splittry = splittry, MIA = MIA, maxsurrogate = maxsurrogate, 
                     numsurrogate = numsurrogate,
                     majority = majority, caseweights = caseweights, 
                     applyfun = applyfun, saveinfo = saveinfo,  ### always
                     testflavour = "ctree", 
                     bonferroni = "Bonferroni" %in% testtype, 
                     splitflavour = "ctree", update = update),
      list(teststat = teststat, splitstat = splitstat, splittest = splittest, pargs = pargs,
           testtype = ttesttype, nresample = nresample, tol = tol,
           intersplit = intersplit))
}


Ctree <- function(formula, data, subset, na.action = na.pass, weights, offset, cluster,
                  control = ctree_control(), ytrafo = NULL, converged = NULL, scores = NULL, 
                  doFit = TRUE, ...) {

    ## set up model.frame() call
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action", "weights", "offset", "cluster", "scores"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$yx <- "none"
    if (is.function(ytrafo)) {
        if (all(c("y", "x") %in% names(formals(ytrafo))))
            mf$yx <- "matrix"
    }
    mf$nmax <- control$nmax
    ## evaluate model.frame
    mf[[1L]] <- quote(extree_data)

    d <- eval(mf, parent.frame())
    subset <- 1:nrow(model.frame(d))

    weights <- model.weights(d)

    if (is.function(ytrafo)) {
        control$update <- TRUE
        nf <- names(formals(ytrafo))
        if (all(c("data", "weights", "control") %in% nf))
            ytrafo <- ytrafo(data = d, weights = weights, control = control)
        stopifnot(all(c("subset", "weights", "info", "estfun", "object") %in% nf) ||
                  all(c("y", "x", "weights", "offset", "start") %in% nf))
    } else {
        control$update <- FALSE
        stopifnot(length(d$variables$x) == 0)
        mfyx <- model.frame(d, yxonly = TRUE)
        mfyx[["(weights)"]] <- mfyx[["(offset)"]] <- NULL
        yvars <- names(mfyx)
        for (yvar in yvars) {
            sc <- d[[yvar, "scores"]]
            if (!is.null(sc))
                attr(mfyx[[yvar]], "scores") <- sc
        }
        Y <- .y2infl(mfyx, response = d$variables$y, ytrafo = ytrafo)
        if (!is.null(iy <- d[["yx", type = "index"]])) {
            Y <- rbind(0, Y)
        } 
        ytrafo <- function(subset, weights, info, estfun, object, ...) list(estfun = Y)
    }
    if (is.function(converged)) {
        stopifnot(all(c("data", "weights", "control") %in% names(formals(converged))))
        converged <- converged(d, weights, control = control)
    } else {
        converged <- TRUE
    }            

    update <- function(subset, weights, control)
        extree_fit(data = d, trafo = ytrafo, converged = converged, partyvars = d$variables$z, 
                   subset = subset, weights = weights, ctrl = control)
    if (!doFit) return(list(d = d, update = update))
    tree <- update(subset = subset, weights = weights, control = control)
    trafo <- tree$trafo
    tree <- tree$nodes
    
    mf <- model.frame(d)
    if (is.null(weights)) weights <- rep(1, nrow(mf))

    fitted <- data.frame("(fitted)" = fitted_node(tree, mf), 
                         "(weights)" = weights,
                         check.names = FALSE)
    fitted[[3]] <- mf[, d$variables$y, drop = TRUE]
    names(fitted)[3] <- "(response)"
    ret <- party(tree, data = mf, fitted = fitted, 
                 info = list(call = match.call(), control = control))
    ret$update <- update
    ret$trafo <- trafo
    class(ret) <- c("constparty", class(ret))

    ### doesn't work for Surv objects
    # ret$terms <- terms(formula, data = mf)
    ret$terms <- d$terms$all
    ### need to adjust print and plot methods
    ### for multivariate responses
    ### if (length(response) > 1) class(ret) <- "party"
    return(ret)
}

library("partykit")

source("extree.R")

iris$block <- sample(gl(3, 1), nrow(iris), replace = TRUE)
iris$x <- sample(gl(3, 1, ordered = TRUE), nrow(iris), replace = TRUE)
iris$Species[1] <- NA
iris$Petal.Length[2] <- NA

a1 <- extree_data(Species ~ Petal.Length | Sepal.Length + Sepal.Width + x, data = iris, na.action = na.pass, 
       offset = Petal.Width, cluster = block, nmax = Inf, yx = "none")

object.size(iris)
object.size(a1)

head(model.frame(a1))
head(model.frame(a1, yxonly = TRUE))
a1[[1]]
a1[[5, "scores"]]
a1[[2, "index"]]
a1[[2, "missing"]]
a1[["yx", "missing"]]

a2 <- extree_data(Species ~ Petal.Length | Sepal.Length + Sepal.Width + x, data = iris, na.action = na.pass, 
       offset = Petal.Width, cluster = block, nmax = Inf, yx = "matrix")

head(model.frame(a2))
head(model.frame(a2, yxonly = TRUE))
a2[[1]]
a2[[5, "scores"]]
a2[[2, "index"]]
a2[[2, "missing"]]
a2[["yx", "missing"]]

a3 <- extree_data(Species ~ Petal.Length | Sepal.Length + Sepal.Width + x, data = iris, na.action = na.pass, 
       offset = Petal.Width, cluster = block, nmax = 5, yx = "none")

head(model.frame(a3))
head(model.frame(a3, yxonly = TRUE))
a3[[1]]
a3[[5, "scores"]]
a3[[2, "index"]]
a3[[2, "missing"]]
a3[["yx", "missing"]]


a4 <- extree_data(Species ~ Petal.Length | Sepal.Length + Sepal.Width + x, data = iris, na.action = na.pass, 
       offset = Petal.Width, cluster = block, nmax = 5, yx = "matrix")

head(model.frame(a1))
head(model.frame(a4, yxonly = TRUE))
a4[[1]]
a4[[5, "scores"]]
a4[[2, "index"]]
a4[[2, "missing"]]
a4[["yx", "missing"]]

data(iris)
iris <- iris[sample(1:nrow(iris), 10000, replace = TRUE),,drop = FALSE]
dim(iris)

system.time(ct <- ctree(Species ~ ., data = iris))
info_node(node_party(ct))

ct

system.time(ct2 <- Ctree(Species ~ ., data = iris))
info_node(node_party(ct2))

ct2

ctrl <- ctree_control(nmax = 50)
system.time(ct3 <- Ctree(Species ~ ., data = iris, control = ctrl))
info_node(node_party(ct3))

ct3

ct4 <- Ctree(Species ~ ., data = iris, control = ctrl, doFit = FALSE)
i <- 1:nrow(iris)
ctrl$update <- FALSE
system.time(replicate(100, ct4$update(sort(sample(i, floor(length(i)/2))), NULL, ctrl)))

library("randomForest")
system.time(randomForest(Species ~ ., data = iris, ntree = 100))

### diabetes

logit <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...)
    glm(y ~ 0 + x, family = binomial, start = start, weights = weights, ...)

data("PimaIndiansDiabetes", package = "mlbench")

PID <- PimaIndiansDiabetes
# PID <- PID[sample(1:nrow(PID), 10000, replace = TRUE),]

library("sandwich")

system.time(d1 <- Ctree(diabetes ~ glucose | pregnant + pressure + triceps + insulin + mass + pedigree + age,
      data = PID, ytrafo = logit))

d1

ctrl <- ctree_control()
ctrl$splitflavour <- "exhaustive"
ctrl$restart <- FALSE
ctrl$breakties <- TRUE

system.time(d2 <- Ctree(diabetes ~ glucose | pregnant + pressure + triceps + insulin + mass + pedigree + age,
      data = PID, ytrafo = logit, control = ctrl))

d2

glmtree(diabetes ~ glucose | pregnant + pressure + triceps + insulin + mass + pedigree + age,
      data = PID, family = binomial())

ctrl <- ctree_control(stump = TRUE)
ctrl$splitflavour <- "exhaustive"
ctrl$testflavour <- "mfluc"
ctrl$restart <- FALSE
ctrl$breakties <- TRUE
ctrl$trim <- .1

source("modelparty.R")

system.time(d2 <- Ctree(diabetes ~ glucose | pregnant + pressure + triceps + insulin + mass + pedigree + age,
      data = PID, ytrafo = logit, control = ctrl))

d2

ctrl <- ctree_control(nmax = 25)

system.time(d3 <- Ctree(diabetes ~ glucose | pregnant + pressure + triceps + insulin + mass + pedigree + age,
      data = PID, ytrafo = logit, control = ctrl))

d3


warnings()

     ### regression
     airq <- subset(airquality, !is.na(Ozone))

ctrl <- ctree_control()#maxsurr = 3)

     ct1 <- ctree(Ozone ~ ., data = airq, control = ctrl)
     ct2 <- Ctree(Ozone ~ ., data = airq, control = ctrl)

max(abs(predict(ct1) - predict(ct2)))

data("GBSG2", package = "TH.data")
library("survival")
GBSG2ct1 <- ctree(Surv(time, cens) ~ ., data = GBSG2)
GBSG2ct2 <- Ctree(Surv(time, cens) ~ ., data = GBSG2)

### multivariate responses
airct1 <- ctree(Ozone + Temp ~ ., data = airq)
airct2 <- Ctree(Ozone + Temp ~ ., data = airq)
all.equal(predict(airct1), predict(airct2))


ct <- Ctree(dist + I(dist^2) ~ speed, data = cars)

