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
    list(estfun = estfun, objfun = -sum((data$yx$y - mean(data$yx$y))^2), converged = TRUE)
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
    
    # estfun <- matrix(NA, ncol = NCOL(data$yx$y), nrow = nrow(data$data))
    # estfun[subset,] <- as.matrix(data$yx$y)[subset, ]
    # estfun[subset,] <- as.data.frame(data$yx$y)[subset, ]
    list(estfun = estfun, converged = TRUE)
}

## var_select via GUIDE
source("selection_modules.R")

var_select2 <- function(model, trafo, data, subset, weights, j, SPLITONLY = FALSE, ctrl) {
    ## TODO: what about subset???
    res <- var_select(estfun = model$estfun, data = data, j = j)
    return(as.list(res))
}

var_select2_call <- function(...) {
    var_sel_fun <- function(model, trafo, data, subset, weights, whichvar, ctrl) {
        args <- list(...)
        ctrl[names(args)] <- args
        partykit:::.select(model, trafo, data, subset, weights, whichvar, ctrl, 
            FUN = var_select2)
    }
    return(var_sel_fun)
}

### how could we allow the following?
## FUN = list(numeric = ..., factor = ..., ordered = ..., default = ...)

### or create the list automatically
## FUN = "foo"
## class(z)
## try(do.call(sprintf("var_select_%s_%s"), FUN, class(z), ...))


## split select with median
## see above


## run tree
tr2 <- extree(data = d, trafo = trafo2, 
    control = c(extree_control(criterion = "p.value",
        logmincriterion = log(1 - 0.05),
        update = TRUE,
        selectfun = var_select2_call(),
        splitfun = split_select1,
        svselectfun = var_select2_call(),
        svsplitfun = split_select1,
        minsplit = 70),
        restart = TRUE))

tr2

