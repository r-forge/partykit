# -------------------------------------------------------------------
# - NAME:   extree_examples_2019_11_08.R
# - AUTHOR: Achim, Heidi, Lisa, Moritz
# - DATE:   2019-11-07
# -------------------------------------------------------------------
# - PURPOSE: Some examples illustrating the proposed extree structure
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2019-11-08 on thinkmoritz
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Prelimanaries
# -------------------------------------------------------------------
library("partykit")
source("selection_modules.R")
source("trafo_functions.R")

# -------------------------------------------------------------------
# EXAMPLE 1: Use trafo_identity() and guide_test() for numerical response
# -------------------------------------------------------------------
## Prepare data
airq <- subset(airquality, !is.na(Ozone))
airq_dat <- extree_data(Ozone ~ Wind + Temp,
    data = airq, yx = "matrix")

## Set up variable selection list
var_select_guide_numeric <- function(model, trafo, data, subset, weights, j, 
    split_only = FALSE, control) {
    
    estfun <- model$estfun[subset]
    
    # categorize estfun if not already a factor
    est_cat <- if (is.factor(estfun)) estfun else cut(estfun, 
        breaks = c(min(estfun), mean(estfun), max(estfun)), 
        include.lowest = TRUE, right = TRUE)
    
    ## get possible split variable
    # sv <- data$zindex[[j]][subset] ## FIXME: always the same as below?
    sv <- extree_variable(data, i = j, type = "index")[subset]
    ## FIXME: can copying be avoided here^?
    
    ## categorize possible split variable
    breaks <- unique(quantile(sv, c(0, 0.25, 0.5, 0.75, 1)))
    if (length(breaks) < 5) breaks <- c(min(sv), mean(sv), max(sv))
    sv_cat <- cut(sv, breaks = breaks, 
        include.lowest = TRUE, right = TRUE)
    
    ## independence test 
    test <- chisq.test(x = est_cat, y = sv_cat) ## FIXME: use libcoin?
    res <- list(statistic = test$statistic, p.value = test$p.value) 
    ## FIXME: (ML, LS) return log(1 - p-value) instead?
    
    return(res)
}

var_select_guide_factor <- function(model, trafo, data, subset, weights, j, 
    split_only = FALSE, control) {
    
    estfun <- model$estfun[subset]
    
    ## categorize estfun if not already a factor
    if(is.factor(estfun)) est_cat <- estfun else {
        breaks <- unique(quantile(estfun, c(0, 0.5, 1)))
        if(length(breaks) < 3) breaks <- c(min(estfun), mean(estfun), max(estfun))
        est_cat <- cut(estfun, breaks = breaks, 
            include.lowest = TRUE, right = TRUE)
    }
    
    ## get possible split variable
    # sv_cat <- data$zindex[[j]][subset] ## FIXME: always the same as below?
    sv_cat <- extree_variable(data, i = j, type = "index")[subset]
    ## FIXME: can copying be avoided here^?
    
    ## independence test
    test <- chisq.test(x = est_cat, y = sv_cat)
    res <- list(statistic = test$statistic, p.value = test$p.value) 
    ## FIXME: (ML, LS) return log(1 - p-value) instead?
    
    return(res)
}


## Set up split selection: split_select with median without index j
split_select_median_numeric <- function(model, trafo, data, subset, weights, 
    whichvar, ctrl) {

    if (length(whichvar) == 0) return(NULL)

    ## split FIRST variable at median 
    ## FIXME: loop if necessary
    j <- whichvar[1]
    x <- model.frame(data)[[j]][subset] ## FIXME: better via extree_variable or model.frame?
    ret <- partysplit(as.integer(j), breaks = median(x))

    return(ret)
}

## Call extree
ctrl1 <- extree_control(criterion = "p.value",
    logmincriterion = log(1 - 0.05),
    update = TRUE,
    selectfun = list(
        numeric = var_select_guide_numeric,
        default = var_select_guide_factor
    ),
    splitfun = split_select_median_numeric,
    svselectfun = NULL, ## FIXME: we need a better default here!!!
    svsplitfun = NULL,
    minsplit = 50)


## trafo
trafo_identity <- function(subset, data, weights = NULL, info = NULL, 
    estfun = TRUE, object = TRUE) {
    
    ## Extract response and "subset"
    y <- extree_variable(data, i = 1, type = "original")  # FIXME: (ML, LS) data copy? no aggregation possible!
    y[-subset] <- NA  
    
    ## Return list
    rval <- list(
        estfun = if (estfun) y else NULL,
        unweighted = TRUE,  
        converged = TRUE 
    )
    
    return(rval)
}


tr1 <- extree(data = airq_dat, trafo = trafo_identity, 
    control = c(ctrl1, restart = TRUE))

tr1


# -------------------------------------------------------------------
# EXAMPLE 2: Use trafo_identity() and guide_test() for categorial response
# -------------------------------------------------------------------
## Prepare data
iris_dat <-  extree_data(Species ~ Petal.Width + Petal.Length,
    data = iris, yx = "matrix")

## Call extree
ctrl2 <- extree_control(criterion = "p.value",
    logmincriterion = log(1 - 0.05),
    update = TRUE,
    selectfun = "guide",
    splitfun = split_select_median_numeric,
    svselectfun = NULL, ## FIXME: we need a better default here!!!
    svsplitfun = NULL,
    minsplit = 50)

tr2 <- extree(data = iris_dat, trafo = trafo_identity,
    control = c(ctrl2, restart = TRUE))

tr2

# -------------------------------------------------------------------
# EXAMPLE 3: Define multiway for anorexia data
# -------------------------------------------------------------------
library("MASS")
data("anorexia")
ano_dat <- extree_data(Postwt ~ Prewt + Treat, data = anorexia, yx = "matrix")

## Split_select multiway
split_select_multiway_factor <- function(model, trafo, data, subset, weights, j,
    split_only = TRUE, control) {

    ## copied from .split
    x <- model.frame(data)[[j]] ## FIXME: use extree_variable instead?

    index <- 1L:nlevels(x)
    xt <- libcoin::ctabs(ix = unclass(x), weights = weights, subset = subset)[-1]
    index[xt == 0] <- NA
    index[xt > 0 & xt < control$minbucket] <- nlevels(x) + 1L
    if (length(unique(index)) == 1) {
        ret <- NULL
    } else {
        index <- unclass(factor(index))
        ret <- partysplit(as.integer(j), index = as.integer(index))
    }
    ## ---

    return(ret)
}

ctrl3 <- extree_control(criterion = "p.value",
    logmincriterion = log(1 - 0.05),
    update = TRUE,
    selectfun = list(
        numeric = var_select_guide_numeric,
        default = var_select_guide_factor
    ),
    splitfun = list(
        numeric = split_select_median_numeric,
        factor = split_select_multiway_factor
    ),
    svselectfun = NULL,
    svsplitfun = NULL,
    minsplit = 50)

tr3 <- extree(data = ano_dat, trafo = trafo_identity,
    control = c(ctrl3, restart = TRUE))

tr3


# # -------------------------------------------------------------------
# # EXAMPLE 4: As example 3, but with character arguments for split/select fun
# # -------------------------------------------------------------------
# ctrl4 <- extree_control(criterion = "p.value",
#     logmincriterion = log(1 - 0.05),
#     update = TRUE,
#     selectfun = "guide",
#     splitfun = list(
#         numeric = "median_numeric",
#         default = "multiway_factor"
#     ),
#     svselectfun = NULL,
#     svsplitfun = NULL,
#     minsplit = 50)
# 
# 
# tr4 <- extree(data = ano_dat, trafo = trafo_identity,
#     control = c(ctrl4, restart = TRUE))
# 
# tr4


# -------------------------------------------------------------------
# EXAMPLE 5: Compare extree with ctree
# -------------------------------------------------------------------
airq <- subset(airquality, !is.na(Ozone))
airq_dat <- extree_data(Ozone ~ Wind + Temp,
    data = airq, yx = "matrix")

## Call original ctree
tr5_ctree <- ctree(Ozone ~ Wind + Temp, data = airq)

## Set up control
ctrl <- extree_control(criterion = "p.value",
  logmincriterion = log(1 - 0.05),
  update = TRUE,
  selectfun = partykit:::.ctree_select(),
  splitfun = partykit:::.ctree_split(),
  svselectfun = NULL,
  svsplitfun = NULL,
  minsplit = 2)

## Add ctree specific control arguments
ctrl5 <- c(ctrl, ctree_control()[!names(ctree_control()) %in% names(ctrl)])

## Call extree 
tr5 <- extree(data = airq_dat, trafo = tr5_ctree$trafo,
    control = ctrl5)


tr5_ctree
tr5
