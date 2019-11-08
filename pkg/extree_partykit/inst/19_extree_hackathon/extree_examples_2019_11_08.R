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
# EXAMPLE 1: Use trafo_num() and guide_test() for airquality data
# -------------------------------------------------------------------
## Prepare data
airq <- subset(airquality, !is.na(Ozone))
airq_dat <- extree_data(Ozone ~ Wind + Temp,
    data = airq, yx = "matrix")

## Set up variable selection list
var_select_guide_num <- function(model, trafo, data, subset, weights, j, split_only = FALSE, control) {
    res <- var_select_num(estfun = model$estfun, data = data, subset = subset, j = j)
    return(as.list(res))
}

var_select_guide_cat <- function(model, trafo, data, subset, weights, j, split_only = FALSE, control) {
    res <- var_select_cat(estfun = model$estfun, data = data, subset = subset, j = j)
    return(as.list(res))
}

var_select_guide <- list(
    numeric = var_select_guide_num,
    default = var_select_guide_cat
)


## Set up split selection: split_select with median without index j
split_select_median <- function(model, trafo, data, subset, weights, whichvar, ctrl) {
    # args <- list(...)

    print(whichvar)

    if (length(whichvar) == 0) return(NULL)

    ## split FIRST variable at median
    j <- whichvar[1]
    x <- model.frame(data)[[j]][subset]
    ret <- partysplit(as.integer(j), breaks = median(x))

    return(ret)
}

## Call extree
tr1_guide <- extree(data = airq_dat, trafo = trafo_num,
    control = c(extree_control(criterion = "p.value",
        logmincriterion = log(1 - 0.05),
        update = TRUE,
        selectfun = var_select_guide,
        splitfun = split_select_median,
        svselectfun = var_select_guide,
        svsplitfun = split_select_median,
        minsplit = 70),
        restart = TRUE))

tr1_guide


# -------------------------------------------------------------------
# EXAMPLE 2: Use trafo_identity() and guide_test() for iris data
# -------------------------------------------------------------------
## Prepare data
iris_dat <-  extree_data(Species ~ Petal.Width + Petal.Length,
    data = iris, yx = "matrix")

## Call extree
tr2_guide <- extree(data = iris_dat, trafo = trafo_identity,
    control = c(extree_control(criterion = "p.value",
        logmincriterion = log(1 - 0.05),
        update = TRUE,
        selectfun = var_select_guide,
        splitfun = split_select_median,
        svselectfun = var_select_guide,
        svsplitfun = split_select_median,
        minsplit = 70),
        restart = TRUE))

tr2_guide

# -------------------------------------------------------------------
# EXAMPLE 3: Define split_select with median and multiway for anorexia data
# -------------------------------------------------------------------
library("MASS")
data("anorexia")
ano_dat <- extree_data(Postwt ~ Prewt + Treat, data = anorexia, yx = "matrix")

## Split_select with median
split_selectmedian_num <- function(model, trafo, data, subset, weights, j,
    split_only = TRUE, control) {

    ## split variable at median
    x <- model.frame(data)[[j]][subset]
    ret <- partysplit(as.integer(j), breaks = median(x))

    return(ret)
}


## Split_select multiway
split_selectmultiway_cat <- function(model, trafo, data, subset, weights, j,
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

## Split_select combined without index j
split_select <- function(model, trafo, data, subset, weights, whichvar, ctrl) {

    split_select_loop(model = model, trafo = trafo, data = data,
        subset = subset, weights = weights, whichvar = whichvar,
        control = ctrl, split_select = list(numeric = split_selectmedian_num,
            factor = split_selectmultiway_cat))
}


ctrl_ano <- extree_control(criterion = "p.value",
    logmincriterion = log(1 - 0.05),
    update = TRUE,
    selectfun = list(
        numeric = var_select_guide_num,
        default = var_select_guide_cat
    ),
    splitfun = list(
        numeric = split_selectmedian_num,
        factor = split_selectmultiway_cat
    ),
    svselectfun = var_select_guide,
    svsplitfun = split_select,
    minsplit = 70)

tr_ano <- extree(data = ano_dat, trafo = trafo_identity,
    control = c(ctrl_ano, restart = TRUE))

tr_ano

# -------------------------------------------------------------------
# EXAMPLE 4: As example 3, but with character arguments for split/select fun
# -------------------------------------------------------------------
var_select_awesome_numeric <- var_select_guide_num
var_select_awesome_default <- var_select_guide_cat

split_select_awesome_numeric <- split_selectmedian_num
split_select_awesome_default <- split_selectmultiway_cat

ctrl2_ano <- extree_control(criterion = "p.value",
    logmincriterion = log(1 - 0.05),
    update = TRUE,
    selectfun = "awesome",
    #selectfun = list(
    #    numeric = "awesome",
    #    default = "awesome"
    #),
    splitfun = "awesome",
    svselectfun = var_select_guide,
    svsplitfun = split_select,
    minsplit = 70)


tr2_ano <- extree(data = ano_dat, trafo = trafo_identity,
    control = c(ctrl2_ano, restart = TRUE))

tr2_ano


# -------------------------------------------------------------------
# EXAMPLE 5: Compare extree with ctree
# -------------------------------------------------------------------
airq <- subset(airquality, !is.na(Ozone))
airq_dat <- extree_data(Ozone ~ Wind + Temp,
    data = airq, yx = "matrix")

## Call original ctree
tr1_ctree <- ctree(Ozone ~ Wind + Temp, data = airq)

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
ctrl <- c(ctrl, ctree_control()[!names(ctree_control()) %in% names(ctrl)])

## Call extree 
tr2_ctree <- extree(data = airq_dat, trafo = tr1_ctree$trafo,
    control = ctrl)
