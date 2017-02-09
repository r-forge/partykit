distforest <- function(formula, data, family = NO(), decorrelate = "none", ntree = 500L,
                       perturb = list(replace = FALSE, fraction = 0.632), fitted.OOB = TRUE,
                       control = ctree_control(...), ocontrol = list(), ...)
{
  ## keep call
  cl <- match.call(expand.dots = TRUE)
  np <- if(inherits(family, "gamlss.family")) family$nopar else length(family$link)
  
  ## wrapper function to apply distfit in cforest
  # input: data, family, weights
  # output: scores (estfun)
  modelscores_decor <- function(data, weights = NULL) {
    
    y <- data[,1]
    #if(survival::is.Surv(y)) y <- data[,1] else y <- as.vector(data[,"y"])
    
    model <- distfit(y, family = family, weights = weights, start = NULL,
                     vcov = (decorrelate == "vcov"), type.hessian = "analytic", estfun = TRUE)
    
    ef <- as.matrix(sandwich::estfun(model))
    #n <- NROW(ef)
    #ef <- ef/sqrt(n)
    
    if(decorrelate != "none") {
      n <- NROW(ef)
      ef <- ef/sqrt(n)
      
      vcov <- if(decorrelate == "vcov") {
        vcov(model) * n
      } else {
        solve(crossprod(ef))
      }
      
      root.matrix <- function(X) {
        if((ncol(X) == 1L)&&(nrow(X) == 1L)) return(sqrt(X)) else {
          X.eigen <- eigen(X, symmetric = TRUE)
          if(any(X.eigen$values < 0)) stop("Matrix is not positive semidefinite")
          sqomega <- sqrt(diag(X.eigen$values))
          V <- X.eigen$vectors
          return(V %*% sqomega %*% t(V))
        }
      }
      ef <- as.matrix(t(root.matrix(vcov) %*% t(ef)))
    }
    return(ef)
  }
  
  
  # rval <- cforest(formula, data, ytrafo = modelscores_decor, ntree = ntree)
  
  ## call cforest
  m <- match.call(expand.dots = FALSE)
  m$ytrafo <- modelscores_decor
  # m$data <- data
  m$family <- m$decorrelate <- NULL
  # for(n in names(ocontrol)) m[[n]] <- ocontrol[[n]]
  # if("..." %in% names(m)) m[["..."]] <- NULL
  # if("type.tree" %in% names(m)) m[["type.tree"]] <- NULL
  m[[1L]] <- as.name("cforest")
  rval <- eval(m, parent.frame())
  
  ### calculate fitted value, fitted distribution parameters, loglikelihood and log scores for every observation
  fitted <- data.frame(idx = 1:nrow(data))
  fitted.par <- data.frame(matrix(0, nrow = nrow(data), ncol = np))
  loglik <- data.frame(idx = 1:nrow(data))
  logscore <- data.frame(matrix(0, nrow = nrow(data), ncol = np))
  
  # extract weights
  w <- predict.cforest(rval, type = "weights", OOB = fitted.OOB)
  for(i in 1:nrow(data)){
    wi <- w[,i]
    # personalized model for observation data[i,]
    pm <-  distfit(data$y, family = family, weights = wi)
    fitted[i,] <- predict(pm)
    fitted.par[i,] <- coef(pm, type = "parameter")
    loglik[i,] <- pm$ddist(data[i,1], log = TRUE)
    logscore[i,] <- pm$familylist$sdist(data[i,1], eta = coef(pm, type = "link"), sum = FALSE)
  }
  
  names(fitted) <- "(response)"
  rval$fitted <- fitted
  
  names(fitted.par) <- names(logscore) <- names(coef(pm, type = "parameter"))
  rval$fitted.par <- fitted.par
  
  rval$loglik <- loglik
  rval$logscore <- logscore
  
  rval$family <- family
  rval$npar <- np
  
  class(rval) <- c("distforest", class(rval))
  return(rval)
}






###################
# methods for class 'distforest'

predict.distforest <- function (object, newdata = NULL, type = c("response", "prob", 
                                                              "weights", "node"), OOB = FALSE, ...) 
{
  responses <- object$fitted[["(response)"]]
  forest <- object$nodes
  nd <- object$data
  vmatch <- 1:ncol(nd)
  if (!is.null(newdata)) {
    nd <- model.frame(delete.response(object$terms), data = newdata, 
                      na.action = na.pass)
    OOB <- FALSE
    vmatch <- match(names(object$data), names(nd))
  }
  nam <- rownames(nd)
  type <- match.arg(type)
  if (type == "node") 
    return(lapply(forest, fitted_node, data = nd, vmatch = vmatch, 
                  ...))
  rw <- object$weights
  applyfun <- lapply
  if (!is.null(object$info)) 
    applyfun <- object$info$control$applyfun
  bw <- applyfun(1:length(forest), function(b) {
    ids <- nodeids(forest[[b]], terminal = TRUE)
    fnewdata <- fitted_node(forest[[b]], nd, vmatch = vmatch, 
                            ...)
    fdata <- fitted_node(forest[[b]], object$data, ...)
    tw <- rw[[b]]
    if (OOB) 
      tw <- as.integer(tw == 0)
    pw <- sapply(ids, function(i) tw * (fdata == i))
    return(pw[, match(fnewdata, ids), drop = FALSE])
  })
  w <- Reduce("+", bw)
  if (!is.matrix(w)) 
    w <- matrix(w, ncol = 1)
  if (type == "weights") {
    ret <- w
    colnames(ret) <- nam
    rownames(ret) <- rownames(responses)
    return(ret)
  }
  
  
  if(type == "response"){
    if(!is.null(newdata)) {
      # get weights for new data
      nw <- predict.cforest(object, newdata = nd, type = "weights", OOB = FALSE)
      
      # calculate prediction for the first observation before the loop in order to get the number of parameters
      pm <-  distfit(object$data$y, family = object$family, weights = nw[,1])
      pred.val1 <- predict(pm)
      pred.par1 <- coef(pm, type = "parameter")
        
      pred.val <- data.frame(idx = 1:nrow(nd))
      pred.val[1,] <- pred.val1
      pred.par <- data.frame(matrix(0, nrow = nrow(nd), ncol = length(pred.par1)))
      pred.par[1,] <- pred.par1
      
      if(nrow(nd)>=2){
        for(i in 2:nrow(nd)){
          nwi <- nw[,i]
          # personalized model
          pm <-  distfit(object$data$y, family = object$family, weights = nwi)
          pred.val[i,] <- predict(pm)
          pred.par[i,] <- coef(pm, type = "parameter")
        }
      }
      pred <- cbind(pred.val, pred.par)
      colnames(pred) <- c("(response)", names(coef(pm, type = "parameter")))
      return(pred)
    } else return(cbind(object$fitted, object$fitted.par))
  }
}



logLik.distforest <- function(object) {
  structure(sum(object$loglik), df = object$npar, class = "logLik")
}










