trafo <- function(data){
  estfun <- data$yx$y
  return(estfun)
}



# variable selection for numeric splitvariable with index j
var_select_num <- function(estfun, data, j){
  
  # categorize estfun if not already a factor
  est_cat <- if(is.factor(estfun)) estfun else cut(estfun, breaks = quantile(estfun, c(0, 0.5, 1)), 
                                                   include.lowest = TRUE, right = TRUE)
  
  # get possible split variable
  sv <- data$zindex[[j]]
  
  # categorize possible split variable
  sv_cat <- cut(sv, breaks = quantile(sv, c(0,0.25, 0.5, 0.75,1)), 
                include.lowest = TRUE, right = TRUE)
  
  # independence test 
  test <- chisq.test(x = est_cat, y = sv_cat)
  res <- c(test$p.value, test$statistic)
  names(res) <- c("p.value", "statistic")

  return(res)
}


# variable selection for categorical splitvariable with index j
var_select_cat <- function(estfun, data, j){
  
  # categorize estfun if not already a factor
  est_cat <- if(is.factor(estfun)) estfun else cut(estfun, breaks = quantile(estfun, c(0, 0.5, 1)), 
                                                   include.lowest = TRUE, right = TRUE)
  
  # get possible split variable
  sv_cat <- data$zindex[[j]]
  
  # independence test
  test <- chisq.test(x = est_cat, y = sv_cat)
  res <- c(test$p.value, test$statistic)
  names(res) <- c("p.value", "statistic")
  
  return(res)
}



# unifying function for variable selection using var_select_cat or var_select_num
var_select <- function(estfun, data, j){
  if(class(data[[j]]) == "factor"){
    res <- var_select_cat(estfun, data, j)
  }
  if(class(data[[j]]) == "numeric"){
    res <- var_select_num(estfun, data, j)
  }
  return(res)
}

# general function for complete data set and all types of splitvariables
var_select_all <- function(estfun, data){
  
  # categorize estfun if not already a factor
  est_cat <- if(is.factor(estfun)) estfun else cut(estfun, breaks = quantile(estfun, c(0, 0.5, 1)), 
                                                    include.lowest = TRUE, right = TRUE)
  
  # select possible splitvariables
  is.splitvar <- !sapply(data$zindex, FUN = is.null)
  
  # store p-value and test statistic in matrix 'res'
  res <- matrix(nrow = sum(is.splitvar), ncol = 2)
  colnames(res) <- c("p.value", "statistic")
  rownames(res) <- names(is.splitvar[is.splitvar])
    
  # independence test over all possible splitvariables
  for(j in 1:length(is.splitvar)){
    if(is.splitvar[j]){
      sv <- data$zindex[[j]]
      # categorize possible splitvariables
      sv_cat <- if(is.factor(sv)) sv else cut(sv, breaks = quantile(sv, c(0,0.25, 0.5, 0.75,1)), 
                                              include.lowest = TRUE, right = TRUE)
      test <- chisq.test(x = est_cat, y = sv_cat)
      res[names(is.splitvar[j]),"p.value"] <- test$p.value
      res[names(is.splitvar[j]),"statistic"] <- test$statistic
    }
  }
  
  return(res)
}




if(FALSE){
  
  library("partykit")
  
  d <-  extree_data(Species ~ Sepal.Width + Sepal.Length | Petal.Width + Petal.Length,
                    data = iris, yx = "matrix")
  
  ef <- trafo(d)
  
  varsel4 <- var_select(ef, d, 4)
  varsel5 <- var_select(ef, d, 5)
  
  varsel <- var_select_all(ef, d)
  
  

}