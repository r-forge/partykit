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
# EXAMPLE 1: Use trafo_num() and guide_test()
# -------------------------------------------------------------------
# - airquality

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

# FIXME: allow that I can give list directly to selectfun
# or even just selectfun = "guide" and the functions need to be called
# var_select_guide_factor, var_select_guide_numeric

var_select_guide_call <- function(model, trafo, data, subset, weights, whichvar, ctrl) {
    var_select_loop(model, trafo, data, subset, weights, whichvar, ctrl,
        var_select = var_select_guide)
}

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

## Trafo with estfun = y, objfun = -MSE
trafo1 <- function(subset, data, weights, info = NULL, estfun = TRUE, object = TRUE) {
    estfun <- matrix(0, ncol = NCOL(data$yx$y), nrow = nrow(data$data))
    estfun[subset,] <- as.matrix(data$yx$y)[subset, ]
    list(estfun = estfun, objfun = sum((data$yx$y[subset] - mean(data$yx$y[subset]))^2), converged = TRUE)
}




tr_guide <- extree(data = airq_dat, trafo = trafo1,
    control = c(extree_control(criterion = "p.value",
        logmincriterion = log(1 - 0.05),
        update = TRUE,
        selectfun = var_select_guide_call,
        splitfun = split_select1,
        svselectfun = var_select_guide_call,
        svsplitfun = split_select1,
        minsplit = 70),
        restart = TRUE))

tr_guide



## airquality data
airq <- subset(airquality, !is.na(Ozone))
data <- extree_data(Ozone ~ Wind + Temp,
    data = airq, yx = "matrix")
subset <- partykit:::.start_subset(data = data)





