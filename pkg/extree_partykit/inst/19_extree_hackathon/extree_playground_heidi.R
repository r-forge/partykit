library("partykitx")

### --- Example 1 --- ###
### - airquality with only numeric variables
### - Trafo with estfun = y, objfun = -MSE
### - varselect with exhaustive search
### - split_select with median

## airquality data
airq <- subset(airquality, !is.na(Ozone))
my_data <- extree_data(Ozone ~ Wind + Temp, 
    data = airq, yx = "matrix")
## README <Z> to get the numeric vector d$yx$y we need to indicate yx = "matrix"?

## Trafo with estfun = y, objfun = -MSE
## README <Z>
## - +MSE ?
## - also include coefficients = mean(data$yx$y[subset]) ?
## - always rely on yx being available?
trafo1 <- function(subset, data, weights, info = NULL, estfun = TRUE, object = TRUE) {
    estfun <- matrix(0, ncol = NCOL(data$yx$y), nrow = nrow(data$data))
    estfun[subset,] <- as.matrix(data$yx$y)[subset, ]
    list(estfun = estfun, objfun = sum((data$yx$y[subset] - mean(data$yx$y[subset]))^2), converged = TRUE)
}

## split_select with median
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

## run tree
tr1 <- extree(data = my_data, trafo = trafo1, 
    control = c(extree_control(
        criterion = "statistic",
        critvalue = -65000,
        update = TRUE,
        selectfun = partykitx:::.objfun_select(),
        splitfun = split_select1,
        svselectfun = partykitx:::.objfun_select(),
        svsplitfun = split_select1),
        restart = TRUE))

tr1
## README <Z> $trafo does not seem to pass on arguments, deliberately?
## plot(tr1)


### --- Example 2 --- ###
### Based on selection_modules.R from Lisa
### - iris data
### - Trafo with estfun = y
### - varselect with GUIDE test
### - split_select with median


## Iris data
d <-  extree_data(Species ~ Petal.Width + Petal.Length,
    data = iris, yx = "matrix")
## README <Z> to get the factor d$yx$y we need to indicate yx = "matrix"?

## Trafo with estfun = y
## README <Z> nevertheless include coefficients and/or objfun?
trafo2 <- function(subset, data, weights, info = NULL, estfun = TRUE, object = TRUE) {
    estfun <- data$yx$y  ## data[[1, "original"]]
    estfun[-subset] <- NA
     
    list(estfun = estfun, converged = TRUE)
}

## varselect via GUIDE
source("selection_modules.R")

varselect2 <- function(model, trafo, data, subset, weights, j, split_only = FALSE, control) {
    res <- varselect(estfun = model$estfun, data = data, subset = subset, j = j)
    return(as.list(res))
}

varselect2_call <- function(model, trafo, data, subset, weights, whichvar, ctrl) {
    # args <- list(...)
    # ctrl[names(args)] <- args
    # partykit:::.select(model, trafo, data, subset, weights, whichvar, ctrl, 
    #     FUN = varselect2)
    varselect_loop(model, trafo, data, subset, weights, whichvar, ctrl, 
        varselect = varselect2)
}


## split select with median
## see above


## run tree
tr2 <- extree(data = d, trafo = trafo2, 
    control = c(extree_control(criterion = "p.value",
        update = TRUE,
        selectfun = varselect2_call,
        splitfun = split_select1,
        svselectfun = varselect2_call,
        svsplitfun = split_select1,
        minsplit = 90),
        restart = TRUE))

tr2



### --- Example 3 --- ###
### Based on selection_modules.R from Lisa
### - iris data
### - Trafo with estfun = y
### - varselect with GUIDE test, but separate functions
### - split_select with median

### how could we allow the following?
## FUN = list(numeric = ..., factor = ..., ordered = ..., default = ...)

varselect3_num <- function(model, trafo, data, subset, weights, j, split_only = FALSE, control) {
    res <- varselect_num(estfun = model$estfun, data = data, subset = subset, j = j)
    return(as.list(res))
}

varselect3_cat <- function(model, trafo, data, subset, weights, j, split_only = FALSE, control) {
    res <- varselect_cat(estfun = model$estfun, data = data, subset = subset, j = j)
    return(as.list(res))
}

varselect3 <- list(
    numeric = varselect3_num,
    default = varselect3_cat
)

## allow that I can give list directly to selectfun
## or even just selectfun = "guide" and the functions need to be called
## varselect_guide_factor, varselect_guide_numeric

varselect3_call <- function(model, trafo, data, subset, weights, whichvar, ctrl) {
    varselect_loop(model, trafo, data, subset, weights, whichvar, ctrl, 
        varselect = varselect3)
}

tr3 <- extree(data = d, trafo = trafo2, 
    control = c(extree_control(criterion = "p.value",
        logmincriterion = log(1 - 0.05),
        update = TRUE,
        selectfun = varselect3_call,
        splitfun = split_select1,
        svselectfun = varselect3_call,
        svsplitfun = split_select1,
        minsplit = 70),
        restart = TRUE))

tr3

all.equal(tr2, tr3)


### --- Example 4 --- ###
### Based on selection_modules.R from Lisa
### - anorexia data
### - Trafo with estfun = y
### - varselect with GUIDE test, but separate functions
### - split_select with median

## data
library("MASS")
data("anorexia")
d4 <- extree_data(Postwt ~ Prewt + Treat, data = anorexia, yx = "matrix")


## split_select with median
split_select4_num <- function(model, trafo, data, subset, weights, j, 
    split_only = TRUE, control) {
    
    ## split variable at median
    x <- model.frame(data)[[j]][subset]
    ret <- partysplit(as.integer(j), breaks = median(x))
    
    return(ret)
}


## split_select multiway
split_select4_cat <- function(model, trafo, data, subset, weights, j, 
    split_only = TRUE, control) {
    
    
    ## --- copied from .split
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

split_select4 <- function(model, trafo, data, subset, weights, whichvar, ctrl) {
    
    split_select_loop(model = model, trafo = trafo, data = data, 
        subset = subset, weights = weights, whichvar = whichvar, 
        control = ctrl, split_select = list(numeric = split_select4_num,
            factor = split_select4_cat))
}


## tree
tr4 <- extree(data = d4, trafo = trafo2, 
    control = c(extree_control(criterion = "p.value",
        logmincriterion = log(1 - 0.05),
        update = TRUE,
        selectfun = varselect3_call,
        splitfun = split_select4,
        svselectfun = varselect3_call,
        svsplitfun = split_select4,
        minsplit = 70),
        restart = TRUE))

tr4


### --- Example 5 --- ###
### Based on selection_modules.R from Lisa
### - anorexia data
### - Trafo with estfun = y
### - varselect with GUIDE test, but separate functions
### - split_select with median

ctrl5 <- extree_control(criterion = "p.value",
    logmincriterion = log(1 - 0.05),
    update = TRUE,
    selectfun = list(
        numeric = varselect3_num,
        default = varselect3_cat
    ),
    splitfun = list(
        numeric = split_select4_num,
        factor = split_select4_cat
    ),
    svselectfun = varselect3_call,
    svsplitfun = split_select4,
    minsplit = 70)

tr5 <- extree(data = d4, trafo = trafo2, 
    control = c(ctrl5, restart = TRUE))

tr5

all.equal(tr4, tr5)

### --- Example 6 --- ###
### Based on selection_modules.R from Lisa
### - anorexia data
### - Trafo with estfun = y
### - varselect with GUIDE test, but separate functions
### - split_select with median

varselect_awesome_numeric <- varselect3_num
varselect_awesome_default <- varselect3_cat

split_select_awesome_numeric <- split_select4_num
split_select_awesome_default <- split_select4_cat

ctrl6 <- extree_control(criterion = "p.value",
    logmincriterion = log(1 - 0.05),
    update = TRUE,
    selectfun = list(
        numeric = varselect_awesome_numeric, 
        default = "awesome_default"
    ),
    splitfun = list(
        numeric = "awesome_numeric",
        default = "awesome_default"
    ),
    svselectfun = varselect3_call,
    svsplitfun = split_select4,
    minsplit = 70)

tr6 <- extree(data = d4, trafo = trafo2, 
    control = c(ctrl6, restart = TRUE))

all.equal(tr6, tr4)


### --- Example 7 --- ###
### Based on selection_modules.R from Lisa
### - anorexia data
### - Trafo with estfun = y
### - varselect with GUIDE test, but separate functions
### - split_select with median

ctrl7 <- extree_control(criterion = "p.value",
    logmincriterion = log(1 - 0.05),
    update = TRUE,
    selectfun = "awesome",
    splitfun = "awesome",
    svselectfun = varselect3_call,
    svsplitfun = split_select4,
    minsplit = 70)


tr7 <- extree(data = d4, trafo = trafo2, 
    control = c(ctrl7, restart = TRUE))

all.equal(tr7, tr4)

