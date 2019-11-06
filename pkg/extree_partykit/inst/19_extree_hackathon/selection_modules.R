trafo <- function(data){
  estfun <- data$yx$y
  return(estfun)
}



# variable selection for numeric splitvariable with index j
var_select_num <- function(estfun, data, subset, j){
  
  estfun <- estfun[subset]
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
  test <- chisq.test(x = est_cat, y = sv_cat)
  res <- c(test$p.value, test$statistic)
  names(res) <- c("p.value", "statistic")  ## FIXME: (ML, LS) return log(p-value) instead?

  return(res)
}


# variable selection for categorical splitvariable with index j
var_select_cat <- function(estfun, data, subset, j){
  
  estfun <- estfun[subset]
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
  test <- chisq.test(x = est_cat, y = sv_cat)
  res <- c(test$p.value, test$statistic)
  names(res) <- c("p.value", "statistic")   ## FIXME: (ML, LS) return log(p-value) instead?
  
  return(res)
}



# unifying function for variable selection using var_select_cat or var_select_num
var_select <- function(estfun, data, subset, j){
  if(class(data[[j]]) == "factor"){
    res <- var_select_cat(estfun = estfun, data = data, subset = subset, j = j)
  }
  if(class(data[[j]]) == "numeric"){
    res <- var_select_num(estfun = estfun, data = data, subset = subset, j = j)
  }
  return(res)
}



# general function for complete data set and all types of splitvariables
var_select_all <- function(estfun, data, subset){
  
  estfun <- estfun[subset]
  # categorize estfun if not already a factor
  if(is.factor(estfun)) {
    est_cat <- estfun 
  } else {
    breaks <- unique(quantile(estfun, c(0, 0.5, 1)))
    if(length(breaks) < 3) breaks <- c(min(estfun), mean(estfun), max(estfun))
    est_cat <- cut(estfun, breaks = breaks, 
                   include.lowest = TRUE, right = TRUE)
  }
  
  # select possible splitvariables
  is.splitvar <- !sapply(data$zindex, FUN = is.null)
  
  # store p-value and test statistic in matrix 'res'
  res <- matrix(nrow = sum(is.splitvar), ncol = 2)
  colnames(res) <- c("p.value", "statistic")         ## FIXME: (ML, LS) return log(p-value) instead?
  rownames(res) <- names(is.splitvar[is.splitvar])
    
  # independence test over all possible splitvariables
  for(j in 1:length(is.splitvar)){
    if(is.splitvar[j]){
      sv <- data$zindex[[j]][subset]
      # categorize possible splitvariable if not already a factor
      if(is.factor(sv)) {
        sv_cat <- sv
      } else {
        breaks <- unique(quantile(sv, c(0,0.25, 0.5, 0.75,1)))
        if(length(breaks) < 5) breaks <- c(min(sv), mean(sv), max(sv))
        sv_cat <- if(is.factor(sv)) sv else cut(sv, breaks = breaks, 
                                                include.lowest = TRUE, right = TRUE)
      }
      test <- chisq.test(x = est_cat, y = sv_cat)
      res[names(is.splitvar[j]),"p.value"] <- test$p.value
      res[names(is.splitvar[j]),"statistic"] <- test$statistic
    }
  }
  
  return(res)
}




if(FALSE){
  
  library("partykit")
  
  d <-  extree_data(Species ~ .,
                    data = iris, yx = "matrix")
  
  ef <- trafo(d)
  
  varsel4 <- var_select(ef, d, subset = c(1:NROW(d$data)), 4)
  varsel5 <- var_select(ef, d, subset = c(1:NROW(d$data)), 5)
  varsel <- var_select_all(ef, d, subset = c(1:NROW(d$data)))
  
  varsel4
  varsel5
  varsel
  
}