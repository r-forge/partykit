ccprune <- function(object, costfunction = NULL, nrfolds = 10) {
  
  # per default use residual sum of squares
  if(is.null(costfunction)) costfunction <- "RSS"   
  
  if(costfunction == "RSS"){
     cf <- function(object){
      if(inherits(object, "partynode"))
        RSS <- sum((object$info$object$model$y - object$info$object$fitted.values)^2)
      if(inherits(object, "modelparty"))
        RSS <- sum((object$data$y - fitted(object))^2)
      return(RSS)
    }
  }
  
  if(is.function(costfunction)) cf <- costfunction
  
  
  
  ###############################################################
  ## 1st step: select finite step of subtrees which are candidates for the optimal tree
  
  ## FIX ME: get first Tree out of T_max = object
  Tree <- object
  T_all <- list(Tree)
  alpha_all <- numeric()
  
  while(width(Tree) > 1){
    
    # ids of inner nodes
    inids <- nodeids(Tree)[!(nodeids(Tree) %in% nodeids(Tree, terminal = TRUE))]
    
    alpha_min <- Inf
    t_min_id <- integer()
    
    for(i in inids){
      
      # costfunction evaluated at tree pruned at inner node i
      Tp <- nodeprune(Tree[[i]], ids = 1)
      cfi <- cf(Tp)
      
      # costfunction evaluated at the branch with node i as root node 
      cfb <- cf(Tree[[i]])
      
      #replace current alpha min if the new alpha is smaller
      alpha <- (cfi-cfb) / (width(Tree[[i]]) - 1)
      
      if(alpha <= alpha_min) {
        if(alpha < alpha_min) t_min_id <- i 
        if(alpha == alpha_min) t_min_id <- cbind(t_min_id, i)
        alpha_min <- alpha
      }
    }
    
    ## cut off branch at node t_min_id
    Tree <- nodeprune(Tree, ids = t_min_id)
    T_all[[length(T_all)+1]] <- Tree
    alpha_all <- c(alpha_all, alpha_min)
  }
  
  
  
  
  ###############################################################
  ## 2nd step: select best tree out of T_all by cross validation

  nobs <- NROW(object$data)
  alpha_mid <- sqrt(alpha_all * c(alpha_all[-1], alpha_all[length(alpha_all)]+0.1))
  buildcl <- object$info$call
  testids <- list()
  allids <- c(1:nobs)
  
  ## store results from cross-validation in a list of matrices (one matrix for each fold)
  # for each tree in the newly built sequence of minimal cost-complex tree store 
  # the corresponding alpha and prediction error on test data
  cvres <- list()
  
  set.seed(7)
  
  for(j in 1:nrfolds){
    testids[[j]] <- sample(allids, round(nobs/nrfolds), replace = FALSE)
    allids <- allids[!(allids %in% testids[[j]])]
  
    # split in to learning data and testdata
    learndata <- object$data[!(c(1:nobs) %in% testids[[j]]),]
    testdata <- object$data[(c(1:nobs) %in% testids[[j]]),]
    
    # build large tree on learning data
    buildcl$data <- learndata
    Tree <- eval(buildcl, parent.frame())
    
    resmat <- matrix(ncol = 2, nrow = 0)
    colnames(resmat) <- c("alpha", "pred.err")
    
    while(width(Tree) > 1){
      
      # ids of inner nodes
      inids <- nodeids(Tree)[!(nodeids(Tree) %in% nodeids(Tree, terminal = TRUE))]
      
      alpha_min <- Inf
      t_min_id <- integer()
      
      for(i in inids){
        
        # costfunction evaluated at tree pruned at inner node i
        Tp <- nodeprune(Tree[[i]], ids = 1)
        cfi <- cf(Tp)
        
        # costfunction evaluated at the branch with node i as root node 
        cfb <- cf(Tree[[i]])
        
        #replace current alpha min if the new alpha is smaller
        alpha <- (cfi-cfb) / (width(Tree[[i]]) - 1)
        
        if(alpha <= alpha_min) {
          if(alpha < alpha_min) t_min_id <- i 
          if(alpha == alpha_min) t_min_id <- cbind(t_min_id, i)
          alpha_min <- alpha
        }
      }
      
      ## cut off branch at node t_min_id
      Tree <- nodeprune(Tree, ids = t_min_id)
      
      resmat <- rbind(resmat, c(alpha_min,
                                mean((testdata$y - predict(Tree, newdata = testdata, 
                                                           type = "response"))^2)))
      
    }
    
    # stores results for fold j
    cvres[[j]] <- resmat
  }
  
  ## estimate error for each tree of T_all by averaging
  # over the errors of the cross validation trees with alpha
  # closest to the alpha of the tree of T_all
  averr <- numeric(length = length(T_all))
  
  ## FIX ME: add 0 as first element to alpha_mid
  # as the maximal Tree corresponds to the minimizing tree for alpha = 0
  alpha_mid <- c(0, alpha_mid)
  
  for(i in 1:length(T_all)){
    
    Tree <- T_all[[i]]
    alpha <- alpha_mid[[i]]
    
    err <- numeric(length = length(nrfolds))
    
    for(j in 1:nrfolds){
      err[j] <- cvres[[j]][which.min(abs(cvres[[j]][,"alpha"] - alpha)), "pred.err"]
    }
    
    averr[i] <- mean(err)
  }
  
  Tree <- T_all[[which.min(averr)]]
  
  return(Tree)
}
