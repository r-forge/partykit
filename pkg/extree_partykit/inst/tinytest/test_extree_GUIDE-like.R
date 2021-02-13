library("partykitx")


### ---- GUIDE variable selection (adapted version of selection_modules.R) ------ ##

# variable selection for numeric splitvariable with index j
var_select_guide_numeric <- function(model, trafo, data, subset, weights, j, 
  split_only = FALSE, control) {
  
  estfun <- model$estfun[subset]
  
  # categorize estfun if not already a factor
  if(is.factor(estfun)) {
    est_cat <- estfun 
  } else {
    breaks <- unique(quantile(estfun, c(0, 0.5, 1)))
    if(length(breaks) < 3) breaks <- c(min(estfun), mean(estfun), max(estfun))
    est_cat <- cut(estfun, breaks = breaks, 
      include.lowest = TRUE, right = TRUE)
  }
  
  # get possible split variable
  sv <- data$zindex[[j]][subset]
  
  # categorize possible split variable
  breaks <- unique(quantile(sv, c(0,0.25, 0.5, 0.75,1)))
  if(length(breaks) < 5) breaks <- c(min(sv), mean(sv), max(sv))
  sv_cat <- cut(sv, breaks = breaks, 
    include.lowest = TRUE, right = TRUE)
  
  # independence test 
  if(length(unique(est_cat)) < 2) return(NULL)
  test <- chisq.test(x = est_cat, y = sv_cat)
  res <- c(test$p.value, test$statistic)
  names(res) <- c("p.value", "statistic")  ## FIXME: (ML, LS) return log(p-value) instead?
  
  return(as.list(res))
}


# variable selection for categorical splitvariable with index j
var_select_guide_factor <- function(model, trafo, data, subset, weights, j, split_only = FALSE, control) {
  
  estfun <- model$estfun[subset]
  # categorize estfun if not already a factor
  if(is.factor(estfun)) {
    est_cat <- estfun 
  } else {
    breaks <- unique(quantile(estfun, c(0, 0.5, 1)))
    if(length(breaks) < 3) breaks <- c(min(estfun), mean(estfun), max(estfun))
    est_cat <- cut(estfun, breaks = breaks, 
      include.lowest = TRUE, right = TRUE)
  }
  
  # get possible split variable
  sv_cat <- data$zindex[[j]][subset]
  
  # independence test
  if(length(unique(est_cat)) < 2) return(NULL)
  test <- chisq.test(x = est_cat, y = sv_cat)
  res <- c(test$p.value, test$statistic)
  names(res) <- c("p.value", "statistic")   ## FIXME: (ML, LS) return log(p-value) instead?
  
  return(as.list(res))
}


## split_select with median
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

### --- Example --- ###
### - iris data
### - Trafo with estfun = y
### - var_select with GUIDE test, but separate functions
### - split_select with median

## Iris data
d <-  extree_data(Species ~ Petal.Width + Petal.Length,
  data = iris, yx = "matrix")


## Trafo with estfun = y
trafo_y <- function(subset, data, weights, info = NULL, estfun = TRUE, object = TRUE) {
  estfun <- data$yx$y  ## data[[1, "original"]]
  estfun[-subset] <- NA
  
  list(estfun = estfun, converged = TRUE)
}





### (1) separate functions for different types of data with argument j
tr1 <- extree(data = d, trafo = trafo_y, 
  control = c(extree_control(criterion = "p.value",
    logmincriterion = log(1 - 0.05),
    update = TRUE,
    varselect = list(
      numeric = var_select_guide_numeric,
      default = var_select_guide_factor
    ),
    splitselect = split_select_median,
    svarselect = list(
      numeric = var_select_guide_numeric,
      default = var_select_guide_factor
    ),
    ssplitselect = split_select_median,
    minbucket = 70,
    lookahead = TRUE),
    restart = TRUE))

# Warnings due to too small tables for chisquare test
tr1

### (2) one function
var_select_guide <- list(
  numeric = var_select_guide_numeric,
  default = var_select_guide_factor
)

var_select_guide_call <- function(model, trafo, data, subset, weights, whichvar, ctrl) {
  var_select_loop(model, trafo, data, subset, weights, whichvar, ctrl, 
    var_select = var_select_guide)
}

tr2 <- extree(data = d, trafo = trafo_y, 
  control = c(extree_control(criterion = "p.value",
    logmincriterion = log(1 - 0.05),
    update = TRUE,
    varselect = var_select_guide_call,
    splitselect = split_select_median,
    svarselect = var_select_guide_call,
    ssplitselect = split_select_median,
    minbucket = 70,
    lookahead = TRUE),
    restart = TRUE))

tr2


all.equal(tr1, tr2)


### (3) character --> function checks if there is a function called
### paste0("^", select_type, "_select_%s"), strategy
### see .get_strategy_function and .get_varclass
### -> naming strategy has to be var_select_type_varclass (i.e. var_select_guide_numeric)
tr3 <- extree(data = d, trafo = trafo_y, 
  control = c(extree_control(criterion = "p.value",
    logmincriterion = log(1 - 0.05),
    update = TRUE,
    varselect = "guide",
    splitselect = split_select_median,
    svarselect = "guide",
    ssplitselect = split_select_median,
    minbucket = 70,
    lookahead = TRUE),
    restart = TRUE))

all.equal(tr1, tr3)


### (4) separate functions for different types of data with argument j

# !!! does not work (yet?) !!!

# var_select_guide_numeric_call <- function(model, trafo, data, subset, weights, whichvar, ctrl) {
#   var_select_loop(model, trafo, data, subset, weights, whichvar, ctrl, 
#     var_select = var_select_guide_numeric)
# }
# 
# 
# tr4 <- extree(data = d, trafo = trafo_y, 
#   control = c(extree_control(criterion = "p.value",
#     logmincriterion = log(1 - 0.05),
#     update = TRUE,
#     varselect = list(
#       numeric = var_select_guide_numeric_call,
#       default = var_select_guide_factor
#     ),
#     splitselect = split_select_median,
#     svarselect = list(
#       numeric = var_select_guide_numeric_call,
#       default = var_select_guide_factor
#     ),
#     ssplitselect = split_select_median,
#     minbucket = 70,
#     lookahead = TRUE),
#     restart = TRUE))
# 
# all.equal(tr1, tr4)