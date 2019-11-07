# -------------------------------------------------------------------
# - NAME:   extree_examples_2019_11_08.R
# - AUTHOR: Achim, Heidi, Lisa, Moritz
# - DATE:   2019-11-07
# -------------------------------------------------------------------
# - PURPOSE: Some examples illustrating the proposed extree structure
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2019-11-07 on thinkmoritz
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Prelimanaries
# -------------------------------------------------------------------
library(partykit)
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


## Set up split selection: split_select with median
split_select1 <- function(model, trafo, data, subset, weights, whichvar, ctrl) {
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
        splitfun = split_select1,
        svselectfun = var_select_guide,
        svsplitfun = split_select1,
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
        splitfun = split_select1,
        svselectfun = var_select_guide,
        svsplitfun = split_select1,
        minsplit = 70),
        restart = TRUE))

tr2_guide


# -------------------------------------------------------------------
# EXAMPLE 3: Compare extree with ctree
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
