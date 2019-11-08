library("partykit")

### --- Example 1 --- ###
### - airquality with only numeric variables
### - Trafo with estfun = y, objfun = -MSE
### - var_select with exhaustive search
### - split_select with median

## airquality data
airq <- subset(airquality, !is.na(Ozone))
my_data <- extree_data(Ozone ~ Wind + Temp, 
    data = airq, yx = "matrix")

## Trafo with estfun = y, objfun = -MSE
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
    control = c(extree_control(criterion = "statistic",
        logmincriterion = log(1 - 0.04),
        update = TRUE,
        selectfun = partykit:::.objfun_select(),
        splitfun = split_select1,
        svselectfun = partykit:::.objfun_select(),
        svsplitfun = split_select1),
        restart = TRUE))

tr1
# plot(tr1)


### --- Example 2 --- ###
### Based on selection_modules.R from Lisa
### - iris data
### - Trafo with estfun = y
### - var_select with GUIDE test
### - split_select with median


## Iris data
d <-  extree_data(Species ~ Petal.Width + Petal.Length,
    data = iris, yx = "matrix")

## Trafo with estfun = y
trafo2 <- function(subset, data, weights, info = NULL, estfun = TRUE, object = TRUE) {
    estfun <- data$yx$y  ## data[[1, "original"]]
    estfun[-subset] <- NA
     
    list(estfun = estfun, converged = TRUE)
}

## var_select via GUIDE
source("selection_modules.R")

var_select2 <- function(model, trafo, data, subset, weights, j, split_only = FALSE, control) {
    res <- var_select(estfun = model$estfun, data = data, subset = subset, j = j)
    return(as.list(res))
}

var_select2_call <- function(model, trafo, data, subset, weights, whichvar, ctrl) {
    # args <- list(...)
    # ctrl[names(args)] <- args
    # partykit:::.select(model, trafo, data, subset, weights, whichvar, ctrl, 
    #     FUN = var_select2)
    var_select_loop(model, trafo, data, subset, weights, whichvar, ctrl, 
        var_select = var_select2)
}


## split select with median
## see above


## run tree
tr2 <- extree(data = d, trafo = trafo2, 
    control = c(extree_control(criterion = "p.value",
        logmincriterion = log(1 - 0.05),
        update = TRUE,
        selectfun = var_select2_call,
        splitfun = split_select1,
        svselectfun = var_select2_call,
        svsplitfun = split_select1,
        minsplit = 70),
        restart = TRUE))

tr2



### --- Example 3 --- ###
### Based on selection_modules.R from Lisa
### - iris data
### - Trafo with estfun = y
### - var_select with GUIDE test, but separate functions
### - split_select with median

### how could we allow the following?
## FUN = list(numeric = ..., factor = ..., ordered = ..., default = ...)

var_select3_num <- function(model, trafo, data, subset, weights, j, split_only = FALSE, control) {
    res <- var_select_num(estfun = model$estfun, data = data, subset = subset, j = j)
    return(as.list(res))
}

var_select3_cat <- function(model, trafo, data, subset, weights, j, split_only = FALSE, control) {
    res <- var_select_cat(estfun = model$estfun, data = data, subset = subset, j = j)
    return(as.list(res))
}

var_select3 <- list(
    numeric = var_select3_num,
    default = var_select3_cat
)

## allow that I can give list directly to selectfun
## or even just selectfun = "guide" and the functions need to be called
## var_select_guide_factor, var_select_guide_numeric

var_select3_call <- function(model, trafo, data, subset, weights, whichvar, ctrl) {
    var_select_loop(model, trafo, data, subset, weights, whichvar, ctrl, 
        var_select = var_select3)
}

tr3 <- extree(data = d, trafo = trafo2, 
    control = c(extree_control(criterion = "p.value",
        logmincriterion = log(1 - 0.05),
        update = TRUE,
        selectfun = var_select3_call,
        splitfun = split_select1,
        svselectfun = var_select3_call,
        svsplitfun = split_select1,
        minsplit = 70),
        restart = TRUE))

tr3

all.equal(tr2, tr3)


### --- Example 4 --- ###
### Based on selection_modules.R from Lisa
### - anorexia data
### - Trafo with estfun = y
### - var_select with GUIDE test, but separate functions
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
        selectfun = var_select3_call,
        splitfun = split_select4,
        svselectfun = var_select3_call,
        svsplitfun = split_select4,
        minsplit = 70),
        restart = TRUE))

tr4


### --- Example 5 --- ###
### Based on selection_modules.R from Lisa
### - anorexia data
### - Trafo with estfun = y
### - var_select with GUIDE test, but separate functions
### - split_select with median

ctrl5 <- extree_control(criterion = "p.value",
    logmincriterion = log(1 - 0.05),
    update = TRUE,
    selectfun = list(
        numeric = var_select3_num,
        default = var_select3_cat
    ),
    splitfun = list(
        numeric = split_select4_num,
        factor = split_select4_cat
    ),
    svselectfun = var_select3_call,
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
### - var_select with GUIDE test, but separate functions
### - split_select with median

var_select_awesome_numeric <- var_select3_num
var_select_awesome_default <- var_select3_cat

split_select_awesome_numeric <- split_select4_num
split_select_awesome_default <- split_select4_cat

ctrl6 <- extree_control(criterion = "p.value",
    logmincriterion = log(1 - 0.05),
    update = TRUE,
    selectfun = list(
        numeric = "awesome_numeric", 
        default = "awesome_default"
    ),
    splitfun = list(
        numeric = "awesome_numeric",
        default = "awesome_default"
    ),
    svselectfun = var_select3_call,
    svsplitfun = split_select4,
    minsplit = 70)

tr6 <- extree(data = d4, trafo = trafo2, 
    control = c(ctrl6, restart = TRUE))

all.equal(tr6, tr4)


### --- Example 7 --- ###
### Based on selection_modules.R from Lisa
### - anorexia data
### - Trafo with estfun = y
### - var_select with GUIDE test, but separate functions
### - split_select with median

ctrl7 <- extree_control(criterion = "p.value",
    logmincriterion = log(1 - 0.05),
    update = TRUE,
    selectfun = "awesome",
    splitfun = "awesome",
    svselectfun = var_select3_call,
    svsplitfun = split_select4,
    minsplit = 70)


tr7 <- extree(data = d4, trafo = trafo2, 
    control = c(ctrl7, restart = TRUE))

all.equal(tr7, tr4)

