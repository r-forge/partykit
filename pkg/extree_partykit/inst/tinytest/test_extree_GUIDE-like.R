### ---- GUIDE variable selection (adapted version of selection_modules.R) ------ ##

## TODO: (HS) support weights, also ask Lisa if we need to implement
## any other functionality.

# variable selection for numeric splitvariable with index j
varselect_guide_numeric <- function(model, trafo, data, subset, weights, j, 
  split_only = FALSE, control) {
  
  ## TODO: (HS) allow matrix estfun
  estfun <- model$estfun[subset]
  
  # categorize estfun if not already a factor
  if(is.factor(estfun)) {
    est_cat <- estfun 
  } else {
    ## TODO (HS): GUIDE uses mean, maybe we should use that in general? -> ask Lisa
    breaks <- unique(quantile(estfun, c(0, 0.5, 1)))
    if(length(breaks) < 3) breaks <- c(min(estfun), mean(estfun), max(estfun))
    est_cat <- cut(estfun, breaks = breaks, 
      include.lowest = TRUE, right = TRUE)
  }
  
  # get possible split variable
  sv <- data$zindex[[j]][subset]
  ## TODO: (HS) use sv <- extree_variable(data, i = j, type = "original")
  
  # categorize possible split variable
  breaks <- unique(quantile(sv, c(0,0.25, 0.5, 0.75,1)))
  if(length(breaks) < 5) breaks <- c(min(sv), mean(sv), max(sv))
  sv_cat <- cut(sv, breaks = breaks, 
    include.lowest = TRUE, right = TRUE)
  
  # independence test 
  ## TODO: (HS) for internal implementation use coin testing instead
  ## of chisq.test
  if(length(unique(est_cat)) < 2) return(NULL)
  test <- chisq.test(x = est_cat, y = sv_cat)
  res <- c(test$p.value, test$statistic)
  names(res) <- c("p.value", "statistic")  
  ## TODO: (ML, LS) return log(p-value) and log-statistic instead
  ## -> c("log.p.value", "log.statistic")
  
  ## TODO: (HS) return matrix
  return(res)
}


# variable selection for categorical splitvariable with index j
varselect_guide_factor <- function(model, trafo, data, subset, weights, j, split_only = FALSE, control) {
  
  ## TODO: (HS) same TODOs as above
  
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
  names(res) <- c("p.value", "statistic") 
  
  return(res)
}


## TODO: (HS) implement a split_select_objfun, based on 
## trafo(...)$objfun on subsets for all possible splits (like MOB or GUIDE)

## split_select with median
split_select_median <- function(model, trafo, data, subset, weights, whichvar, ctrl) {
  # args <- list(...)
  
  # print(whichvar)
  
  if (length(whichvar) == 0) return(NULL)
  
  ## split FIRST variable at median
  ## TODO: (HS) consider what to do if first does not work (use control info)
  j <- whichvar[1]
  x <- model.frame(data)[[j]][subset] ## TODO: (HS) use extree_variable here
  ret <- partysplit(as.integer(j), breaks = median(x))
  
  return(ret)
}

### --- Example --- ###
### - iris data
### - Trafo with estfun = y
### - varselect with GUIDE test, but separate functions
### - split_select with median

## Iris data
d <-  extree_data(Species ~ Petal.Width + Petal.Length,
  data = iris, yx = "matrix") ## TODO: (HS) we shouldn't need "matrix" here

## Trafo with estfun = y
trafo_y <- function(subset, data, weights, info = NULL, estfun = TRUE, object = TRUE) {
  estfun <- extree_variable(data, variable = "y")
  estfun[-subset] <- NA
  
  list(estfun = estfun, converged = TRUE)
}





### (1) separate functions for different types of data with argument j
suppressWarnings(  # FIXME: (ML) Resolve warning in Chi-squared test
  tr1 <- extree(data = d, trafo = trafo_y, 
    control = c(extree_control(criterion = "p.value",
      critvalue = 0.05,
      update = TRUE,
      varselect = list(
        numeric = varselect_guide_numeric,
        default = varselect_guide_factor
      ),
      splitselect = split_select_median,
      svarselect = list(
        numeric = varselect_guide_numeric,
        default = varselect_guide_factor
      ),
      ssplitselect = split_select_median,
      minbucket = 70,
      lookahead = TRUE),
      restart = TRUE))
)
# Warnings due to too small tables for chisquare test
## TODO: (HS) implement without warnings
# tr1

### (2) one function
varselect_guide <- list(
  numeric = varselect_guide_numeric,
  default = varselect_guide_factor
)

varselect_guide_call <- function(model, trafo, data, subset, weights, whichvar, ctrl) {
  varselect_loop(model, trafo, data, subset, weights, whichvar, ctrl, 
    varselect = varselect_guide)
}
suppressWarnings(  # FIXME: (ML) Resolve warning in Chi-squared test
  tr2 <- extree(data = d, trafo = trafo_y, 
    control = c(extree_control(criterion = "p.value",
      logmincriterion = log(1 - 0.05),
      update = TRUE,
      varselect = varselect_guide_call,
      splitselect = split_select_median,
      svarselect = varselect_guide_call,
      ssplitselect = split_select_median,
      minbucket = 70,
      lookahead = TRUE),
      restart = TRUE))
)
# tr2


expect_equal(tr1, tr2)


### (3) character --> function checks if there is a function called
### paste0("^", select_type, "_select_%s"), strategy
### see .get_strategy_function and .get_varclass
### -> naming strategy has to be varselect_type_varclass (i.e. varselect_guide_numeric)

## FIXME: (ML) Resolve tree, which currently leads to an error:
##               * `varselect = "guide"`, `svarselect = "guide"` are the problem
##               * these lead to an empty `strategy` in .get_strategy_function()
##               * maybe due to execution env?
##suppressWarnings(  # FIXME: (ML) Resolve warning in Chi-squared test
##  tr3 <- extree(data = d, trafo = trafo_y, 
##    control = c(extree_control(criterion = "p.value",
##      logmincriterion = log(1 - 0.05),
##      update = TRUE,
##      varselect = "guide",
##      splitselect = split_select_median,
##      svarselect = "guide",
##      ssplitselect = split_select_median,
##      minbucket = 70,
##      lookahead = TRUE),
##      restart = TRUE))
##)
##expect_equal(tr1, tr3)


### (4) separate functions for different types of data with argument j

# !!! does not work (yet?) !!!

# varselect_guide_numeric_call <- function(model, trafo, data, subset, weights, whichvar, ctrl) {
#   varselect_loop(model, trafo, data, subset, weights, whichvar, ctrl, 
#     varselect = varselect_guide_numeric)
# }
# 
# 
# tr4 <- extree(data = d, trafo = trafo_y, 
#   control = c(extree_control(criterion = "p.value",
#     logmincriterion = log(1 - 0.05),
#     update = TRUE,
#     varselect = list(
#       numeric = varselect_guide_numeric_call,
#       default = varselect_guide_factor
#     ),
#     splitselect = split_select_median,
#     svarselect = list(
#       numeric = varselect_guide_numeric_call,
#       default = varselect_guide_factor
#     ),
#     ssplitselect = split_select_median,
#     minbucket = 70,
#     lookahead = TRUE),
#     restart = TRUE))
# 
# all.equal(tr1, tr4)
