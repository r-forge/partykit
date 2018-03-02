ccprune <- function(object, costfunction = NULL) {
  
  is.null(costfunction) costfunction <- "err"
  
  if(costfunction == "err"){
     cf <- function(object){
      if(inherits(object, "partynode"))
        err <- sum((object$info$object$model$y - object$info$object$fitted.values)^2)
      if(inherits(object, "modelparty"))
        err <- sum((object$data$y - fitted(object))^2)
      return(err)
    }
    
  }
  
  if(is.function(costfunction)) cf <- costfunction
  
  ## FIX ME: get T1 out of T_max = object
  T1 <- object
  Tree <- T1
  T_all <- list(T1)
  
  while(width(Tree) > 1){
    
    # ids of inner nodes
    inids <- nodeids(Tree)[!(nodeids(Tree) %in% nodeids(Tree, terminal = TRUE))]
    
    alpha_min <- Inf
    t_min_id <- integer()
    
    for(i in inids){
      
      # costfunction evaluated at the inner node i
      cfi <- partykit::nodeapply(Tree, ids = i, FUN = cf)[[1]]
      
      # costfunction evaluated at the branch with node i as root node 
      tnids <- partykit::nodeapply(Tree, ids = i, FUN = function(node) nodeids(node, terminal = TRUE))[[1]]
      cfb <- sum(unlist(partykit::nodeapply(Tree, ids = tnids, FUN = cfbranch)))
      
      #replace current alpha min if the new alpha is smaller
      alpha <- (cfi-cfb) / (length(tnids) - 1)
      if(alpha <= alpha_min) {
        if(alpha < alpha_min) t_min_id <- i 
        if(alpha == alpha_min) t_min_id <- cbind(t_min_id, i)
        alpha_min <- alpha
      }
    }
    
    ## cut off branch at node t_min_id
    Tree <- nodeprune(Tree, ids = t_min_id)
    T_all[[length(T_all)+1]] <- Tree
    
  }
  
  ## FIX ME: select best tree out of T_all by cross validation

}
