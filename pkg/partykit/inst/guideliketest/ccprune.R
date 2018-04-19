ccprune <- function(object, costfunction = NULL, nrfolds = 10, SE1rule = FALSE, seed = 7) {
  
  cl <- match.call()
  
  # per default use residual sum of squares
  if(is.null(costfunction)) costfunction <- "RSS"   
  
  if(costfunction == "RSS"){
     cf <- function(object, newdata = NULL, mean = FALSE){
       if(is.null(newdata)) {
         if(mean) return(mean((object$data$y - fitted(object))^2)) 
         return(sum((object$data$y - fitted(object))^2))
       } else {
         if(mean) return(mean((newdata$y - predict(object, newdata = newdata, type = "response"))^2))
         return(sum((newdata$y - predict(object, newdata = newdata, type = "response"))^2))
       }
    }
  }
  
  
  ## FIX ME: 
  # - degrees of freedom for AIC?
  # - evaluation on newdata for cross validation
  if(costfunction == "AIC"){
    cf <- function(object, newdata = NULL, mean = FALSE){
      ll <- logLik_modelparty_lm_nd(object, newdata = newdata)
      aic <- (-2 * as.numeric(ll) + 2 * attr(ll, "df"))
      return(aic)
    }
  }
  
  
  if(costfunction == "loglikelihood"){
    cf <- function(object, newdata = NULL, mean = FALSE){
      ll <- as.numeric(logLik_modelparty_lm_nd(object, newdata = newdata))
      return(ll)
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
      cfi <- cf(Tp, mean = FALSE)  
      ## sum instead of mean because the nr of obs in this subtree is of importance
      # (compare R(T): includes probability of being in node t which is nobs(t)/nobs)
      
      # costfunction evaluated at the branch with node i as root node 
      cfb <- cf(Tree[[i]], mean = FALSE)  
      ## sum instead of mean because the nr of obs in this subtree is of importance
      # (compare R(T): includes probability of being in node t which is nobs(t)/nobs)
      
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
  
  set.seed(seed)
  
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
        cfi <- cf(Tp, mean = FALSE)
        
        # costfunction evaluated at the branch with node i as root node 
        cfb <- cf(Tree[[i]], mean = FALSE)
        
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
                                cf(Tree, newdata = testdata, mean = TRUE)))
      
    }
    
    # stores results for fold j
    cvres[[j]] <- resmat
  }
  
  ## estimate error for each tree of T_all by averaging
  # over the errors of the cross validation trees with alpha
  # closest to the alpha of the tree of T_all
  averr <- numeric(length = length(T_all))
  if(SE1rule) serr <- numeric(length = length(T_all))
  
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
    if(SE1rule) serr[i] <- sd(err)/sqrt(length(err))
  }
  
  # select optimal tree
  mintreeid <- which.min(averr)
  if(SE1rule) {
    # select the smallest tree of those beneath the upper error limit (mean + standard error)
    errlimit <- averr[mintreeid] + serr[mintreeid]
    candidatesid <- c(1:length(T_all))[averr <= errlimit]
    Tree <- T_all[[max(candidatesid)]]
  } else {
    Tree <- T_all[[mintreeid]]
  }
  
  return(list(Tree = Tree,
              averr = averr,
              call = cl))
}
