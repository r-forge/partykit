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
var_select_guide_numeric <- function(model, trafo, data, subset, weights, j, split_only = FALSE, control) {
    res <- var_select_num(estfun = model$estfun, data = data, subset = subset, j = j)
    return(as.list(res))
}

var_select_guide_factor <- function(model, trafo, data, subset, weights, j, split_only = FALSE, control) {
    res <- var_select_cat(estfun = model$estfun, data = data, subset = subset, j = j)
    return(as.list(res))
}

var_select_guide <- list(
    numeric = var_select_guide_numeric,
    default = var_select_guide_factor
)

## Set up split selection: split_select with median without index j
split_select_median_numeric <- function(model, trafo, data, subset, weights, whichvar, ctrl) {

    if (length(whichvar) == 0) return(NULL)

    ## split FIRST variable at median
    j <- whichvar[1]
    x <- model.frame(data)[[j]][subset]
    ret <- partysplit(as.integer(j), breaks = median(x))

    return(ret)
}

## Call extree
tr1 <- extree(data = airq_dat, trafo = trafo_num,
    control = c(extree_control(criterion = "p.value",
        logmincriterion = log(1 - 0.05),
        update = TRUE,
        selectfun = var_select_guide,
        splitfun = split_select_median_numeric,
        svselectfun = .ctree_select(),
        svsplitfun = .ctree_split(),
        minsplit = 20),
        restart = TRUE))

tr1


# -------------------------------------------------------------------
# EXAMPLE 2: Use trafo_identity() and guide_test() for categorial response
# -------------------------------------------------------------------
## Prepare data
iris_dat <-  extree_data(Species ~ Petal.Width + Petal.Length,
    data = iris, yx = "matrix")

## Call extree
tr2 <- extree(data = iris_dat, trafo = trafo_identity,
    control = c(extree_control(criterion = "p.value",
        logmincriterion = log(1 - 0.05),
        update = TRUE,
        selectfun = var_select_guide,
        splitfun = split_select_median_numeric,
        svselectfun = .ctree_select(),
        svsplitfun = .ctree_split(),
        minsplit = 70),
        restart = TRUE))

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
    x <- model.frame(data)[[j]]

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
    svselectfun = .ctree_select(),
    svsplitfun = .ctree_split(),
    minsplit = 70)

tr3 <- extree(data = ano_dat, trafo = trafo_identity,
    control = c(ctrl3, restart = TRUE))

tr3


# -------------------------------------------------------------------
# EXAMPLE 4: As example 3, but with character arguments for split/select fun
# -------------------------------------------------------------------
ctrl4 <- extree_control(criterion = "p.value",
    logmincriterion = log(1 - 0.05),
    update = TRUE,
    selectfun = "guide",
    splitfun = list(
        numeric = "median_numeric",
        default = "multiway_factor"
    ),
    svselectfun = .ctree_select(),
    svsplitfun = .ctree_split(),
    minsplit = 70)


tr4 <- extree(data = ano_dat, trafo = trafo_identity,
    control = c(ctrl4, restart = TRUE))

tr4


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
  selectfun = .ctree_select(),
  splitfun = .ctree_split(),
  svselectfun = .ctree_select(),
  svsplitfun = .ctree_split(),
  minsplit = 2)

## Add ctree specific control arguments
ctrl5 <- c(ctrl, ctree_control()[!names(ctree_control()) %in% names(ctrl)])

## Call extree 
tr5 <- extree(data = airq_dat, trafo = tr5_ctree$trafo,
    control = ctrl5)
