ccprune <- function(object, costfunction = NULL, nrfolds = 10, SE1rule = FALSE, seed = 7) {
  
  cl <- match.call()
  
  # per default use residual sum of squares
  if(is.null(costfunction)) costfunction <- "RSS"   
  
  if(costfunction == "RSS"){
    cf <- function(predcoef, data, mean_error = FALSE){
      if(mean_error) return(mean((data$y - (predcoef[1] + predcoef[2] * data$x))^2)) 
      return(sum((data$y - (predcoef[1] + predcoef[2] * data$x))^2))
    }
  }

  
  ## other costfunctions have to be adapted to new structure of cf
  if(FALSE){
    ## FIX ME: 
    # - degrees of freedom for AIC?
    # - evaluation on newdata for cross validation
    if(costfunction == "AIC"){
      cf <- function(object, newdata = NULL, mean_error = FALSE){
        ll <- logLik_modelparty_lm_nd(object, newdata = newdata)
        aic <- (-2 * as.numeric(ll) + 2 * attr(ll, "df"))
        return(aic)
      }
    }
    
    
    if(costfunction == "loglikelihood"){
      cf <- function(object, newdata = NULL, mean_error = FALSE){
        ll <- as.numeric(logLik_modelparty_lm_nd(object, newdata = newdata))
        return(ll)
      }
    }
  }  
   
   
  if(is.function(costfunction)) cf <- costfunction
  
  ## function to evaluate cf for given node IDs (in the original tree 'object')
  eval_cf <- function(object, nids, mean_error = FALSE){
    cfsum <- 0
    for(i in nids) {
      node_cf <- cf(predcoef = coef(nodeprune(object[[i]], ids = 1)),
                    data = object[[i]]$data, mean_error = mean_error)
      if(mean_error) node_cf <- node_cf * NROW(object[[i]]$data) / NROW(object$data)
      cfsum <- cfsum + node_cf
    }
    return(cfsum)
  }
  
  
  ###############################################################
  ## 1st step: select finite step of subtrees which are candidates for the optimal tree
  
  ## FIX ME: get first Tree out of T_max = object
  alpha_all <- c(0)
  tnids_all <- inids_all <- list()
  
  # ids of terminal nodes, all nodes and inner nodes
  tnids <- nodeids(object, terminal = TRUE)
  anids <- nodeids(object)
  inids <- nodeids(object)[!(anids %in% tnids)]
  tnids_all[[1]] <- tnids
  inids_all[[1]] <- inids
  
  
  while(length(tnids) > 1){
    
    alpha_min <- Inf
    t_min_id <- integer()
    
    for(i in inids){
      
      # ids of terminal nodes of branch with root node i
      bnids <- anids[anids %in% names(object[[i]]) & anids %in% tnids]
      
      # costfunction evaluated at tree pruned at inner node i
      cfp <- eval_cf(object, i, mean_error = FALSE)  
      ## sum instead of mean because the nr of obs in this subtree is of importance
      # (compare R(T): includes probability of being in node t which is nobs(t)/nobs)
      
      # costfunction evaluated at the branch with node i as root node 
      cfb <- eval_cf(object, bnids, mean_error = FALSE)  
      ## sum instead of mean because the nr of obs in this subtree is of importance
      # (compare R(T): includes probability of being in node t which is nobs(t)/nobs)
      
      #replace current alpha min if the new alpha is smaller
      alpha <- (cfp-cfb) / (length(bnids) - 1)
      
      if(alpha <= alpha_min) {
        if(alpha < alpha_min) t_min_id <- i 
        if(alpha == alpha_min) t_min_id <- cbind(t_min_id, i)
        alpha_min <- alpha
      }
    }
    
    # in case of more than one split minimizing alpha
    # select the minimal tree (hence the largest branch to be cut off)
    if(length(t_min_id)>1) {
      t_min_id_smallest <- t_min_id[1]
      for(i in 2:length(t_min_id)) {
        if(width(object[[t_min_id[i]]]) > width(object[[t_min_id_smallest]])) t_min_id_smallest <- t_min_id[i]
      }
      t_min_id <- t_min_id_smallest
    }
      
    
    ## cut off branch at node t_min_id
    # by turning inner node t_min_id into a terminal node
    # and deleting its child nodes from the list of terminal nodes
    inids <- inids[!(inids %in% names(object[[t_min_id]]))]
    tnids <- sort(c(t_min_id, tnids[!(tnids %in% names(object[[t_min_id]]))]))
    
    inids_all[[length(inids_all)+1]] <- inids
    tnids_all[[length(tnids_all)+1]] <- tnids
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
    
    ## randomly select ids for learning and testing data from allids
    # last fold: testing data is all the data left after all earlier selections
    if(j != nrfolds) {
      testids[[j]] <- sample(allids, round(nobs/nrfolds), replace = FALSE)
      allids <- allids[!(allids %in% testids[[j]])]
    } else {
      testids[[j]] <- allids
    }  
    # split in to learning data and testdata
    learndata <- object$data[!(c(1:nobs) %in% testids[[j]]),]
    testdata <- object$data[(c(1:nobs) %in% testids[[j]]),]
    
    # build large tree on learning data
    buildcl$data <- learndata
    Tree <- eval(buildcl, parent.frame())
    
    tnids <- nodeids(Tree, terminal = TRUE)
    anids <- nodeids(Tree)
    inids <- nodeids(Tree)[!(anids %in% tnids)]
    
    resmat <- matrix(ncol = 2, nrow = 0)
    colnames(resmat) <- c("alpha", "pred.err")
    
    while(length(tnids) > 1){
      
      alpha_min <- Inf
      t_min_id <- integer()
      
      for(i in inids){
        
        # ids of terminal nodes of branch with root node i
        bnids <- anids[anids %in% names(Tree[[i]]) & anids %in% tnids]
        
        # costfunction evaluated at tree pruned at inner node i
        cfp <- eval_cf(Tree, i, mean_error = FALSE)  
        ## sum instead of mean because the nr of obs in this subtree is of importance
        # (compare R(T): includes probability of being in node t which is nobs(t)/nobs)
        
        # costfunction evaluated at the branch with node i as root node 
        cfb <- eval_cf(Tree, bnids, mean_error = FALSE)  
        ## sum instead of mean because the nr of obs in this subtree is of importance
        # (compare R(T): includes probability of being in node t which is nobs(t)/nobs)
        
        #replace current alpha min if the new alpha is smaller
        alpha <- (cfp-cfb) / (length(bnids) - 1)
        
        if(alpha <= alpha_min) {
          if(alpha < alpha_min) t_min_id <- i 
          if(alpha == alpha_min) t_min_id <- cbind(t_min_id, i)
          alpha_min <- alpha
        }
      }
      
      
      # in case of more than one split minimizing alpha
      # select the minimal tree (hence the largest branch to be cut off)
      if(length(t_min_id)>1) {
        t_min_id_smallest <- t_min_id[1]
        for(k in 2:length(t_min_id)) {
          if(width(Tree[[t_min_id[k]]]) > width(Tree[[t_min_id_smallest]])) t_min_id_smallest <- t_min_id[k]
        }
        t_min_id <- t_min_id_smallest
      }
      
      
      ## cut off branch at node t_min_id
      # by turning inner node t_min_id into a terminal node
      # and deleting its child nodes from the list of terminal nodes
      inids <- inids[!(inids %in% names(Tree[[t_min_id]]))]
      tnids <- sort(c(t_min_id, tnids[!(tnids %in% names(Tree[[t_min_id]]))]))
      
      
      ## get predicted node ids for testdata in original tree
      prednodes_test <- predict(Tree, newdata = testdata)
      # if the predicted node has been pruned: replace by parent node which is also terminal
      for(m in 1:length(prednodes_test)){
        if(!prednodes_test[m] %in% tnids) {
          for(n in tnids) {
            if(prednodes_test[m] %in% as.numeric(names(Tree[[n]]))) prednodes_test[m] <- n
          }
        }
      }
      
      testerror <- 0
      for(l in c(1:NROW(testdata))){
        testerror <- testerror + cf(predcoef = coef(nodeprune(Tree[[prednodes_test[l]]], ids = 1)),
                                    data = testdata[l,], mean_error = FALSE)
      }
      testerror <- testerror/NROW(testdata)
        
      resmat <- rbind(resmat, c(alpha_min, testerror))
      
    }
    
    # stores results for fold j
    cvres[[j]] <- resmat
  }
  
  ## estimate error for each tree of T_all by averaging
  # over the errors of the cross validation trees with alpha
  # closest to the alpha of the tree of T_all
  averr <- numeric(length = length(tnids_all))
  if(SE1rule) serr <- numeric(length = length(tnids_all))
  
  for(i in 1:length(tnids_all)){
    
    tnids <- tnids_all[[i]]
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
    candidatesid <- c(1:length(tnids_all))[averr <= errlimit]
    tnids <- tnids_all[[max(candidatesid)]]
  } else {
    tnids <- tnids_all[[mintreeid]]
  }
  
  ## build optimal tree from list of its terminal node ids tnids
  ## FIX ME: check that data is taken along correctly
  Tree <- nodeprune(object, ids = tnids)
  
  return(list(Tree = Tree,
              averr = averr,
              call = cl))
}
