distforest <- function(formula, data, family = NO(), decorrelate = "none", ntree = 500L, fit = TRUE,
                       perturb = list(replace = FALSE, fraction = 0.632), fitted.OOB = TRUE,
                       cens = "none", censpoint = NULL,
                       control = ctree_control(teststat = "quad", testtype = "Univ", mincriterion = 0, ...), 
                       ocontrol = list(), ...)
{
  ## keep call
  cl <- match.call(expand.dots = TRUE)
  if(is.function(family)) family <- family()
  np <- if(inherits(family, "gamlss.family")) family$nopar else length(family$link)
  
  resp.name <- as.character(formula[2])
  
  ## wrapper function to apply distfit in cforest
  ytrafo <- function(formula, data, weights = NULL, cluster = cluster, ctrl = control) {
    
    if(!(is.null(cluster))) stop("FIX: cluster ignored by trafo-function")
    if(!(is.numeric(formula[[3]]))) {
      #print(formula)
      stop("covariates can only be used as splitting variables (formula has to be of type y~1|x or y~0|x)")
    }
    
    decorrelate <- if(is.null(ctrl$decorrelate)) "none" else ctrl$decorrelate  # FIX ME: include in ctrl?
    
    modelscores_decor <- function(subset, estfun = TRUE, object = TRUE, info = NULL) {
      
      ys <- data[subset,resp.name]
      subweights <- if(is.null(weights) || (length(weights)==0L)) weights else weights[subset] ## FIX ME: scores with or without weights?
      # start <- if(!(is.null(info$coefficients))) info$coefficients else NULL
      start <- info$coefficients
      
      model <- distfit(ys, family = family, weights = subweights, start = start,
                       vcov = (decorrelate == "vcov"), type.hessian = "analytic", 
                       estfun = estfun, cens = cens, censpoint = censpoint, ...)
      
      
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
      
      if(estfun) {
        estfun <- matrix(0, ncol = ncol(ef), nrow = nrow(data)) 
        estfun[subset,] <- ef
      } else estfun <- NULL
      
      object <-  if(object) model else NULL
      
      ret <- list(estfun = estfun,
                  coefficients = coef(model, type = "parameter"),
                  objfun = logLik(model),  # optional function to be maximized (FIX: negative?/minimize?)
                  object = object,
                  converged = TRUE  # FIX ME: warnings is distfit does not converge
      )
      return(ret)
    }
    return(modelscores_decor)
  }    
  
  
  # rval <- cforest(formula, data, ytrafo = modelscores_decor, ntree = ntree)
  
  ## call cforest
  m <- match.call(expand.dots = FALSE)
  m$ytrafo <- ytrafo
  m$control <- control
  # m$data <- data
  m$family <- m$decorrelate <- m$cens <- m$censpoint <- NULL
  # for(n in names(ocontrol)) m[[n]] <- ocontrol[[n]]
  # if("..." %in% names(m)) m[["..."]] <- NULL
  # if("type.tree" %in% names(m)) m[["type.tree"]] <- NULL
  m[[1L]] <- as.name("cforest")
  rval <- eval(m, parent.frame())
  
  if(fit) {
    ### calculate fitted value, fitted distribution parameters, loglikelihood (and log scores) for every observation
    fitted <- data.frame(idx = 1:nrow(data))
    fitted.par <- data.frame(matrix(0, nrow = nrow(data), ncol = np))
    loglik <- data.frame(idx = 1:nrow(data))
    # logscore <- data.frame(matrix(0, nrow = nrow(data), ncol = np))
    
    # extract weights
    w <- predict.cforest(rval, type = "weights", OOB = fitted.OOB)
    
    for(i in 1:nrow(data)){
      wi <- w[,i]
      # personalized model for observation data[i,]
      pm <-  distfit(data[,resp.name], family = family, weights = wi, vcov = FALSE, cens = cens, censpoint = censpoint)
      fitted[i,] <- predict(pm, type = "response")
      fitted.par[i,] <- coef(pm, type = "parameter")
      if(inherits(family, "gamlss.family")) {
        if(("censored" %in% strsplit(family$family[2], " ")[[1]]) && (!survival::is.Surv(data[i,resp.name]))) {
          if(cens == "left") loglik[i,] <- pm$ddist(survival::Surv(data[i,resp.name], data[i,resp.name] > censpoint, type = "left"), log = TRUE)
          if(cens == "right") loglik[i,] <- pm$ddist(survival::Surv(data[i,resp.name], data[i,resp.name] < censpoint, type = "right"), log = TRUE)
          ## FIX ME: interval censored
          #if(cen == "interval") y <- survival::Surv(y, ((y > censpoint[1]) * (y < censpoint[2])), type = "interval")
        } else loglik[i,] <- pm$ddist(data[i,resp.name], log = TRUE)
      } else loglik[i,] <- pm$ddist(data[i,resp.name], log = TRUE)
      
      # logscore[i,] <- pm$familylist$sdist(data[i,resp.name], eta = coef(pm, type = "link"), sum = FALSE)
    }
    
    if(is.null(weights) || (length(weights)==0L || is.function(weights))) weights <- numeric(nrow(data)) + 1
    rval$fitted <- as.data.frame(weights, nrow = length(weights))
    colnames(rval$fitted) <- "(weights)"
    colnames(fitted) <- "(fitted.response)"
    rval$fitted$`(response)` <- data[,resp.name]
    rval$fitted$`(fitted.response)` <- fitted
    #colnames(rval$fitted) <- c("(weights)", "(response)", "(fitted.response)")
    
    
    names(fitted.par) <- names(coef(pm, type = "parameter"))
    # names(logscore) <- names(coef(pm, type = "parameter"))
    rval$fitted.par <- fitted.par
    
    rval$loglik <- sum(loglik)
    # rval$logscore <- logscore
  }
  
  rval$info$family <- family
  rval$info$npar <- np
  rval$info$formula <- formula
  
  class(rval) <- c("distforest", class(rval))
  return(rval)
}






###################
# methods for class 'distforest'

predict.distforest <- function (object, newdata = NULL, type = c("response", "parameter", "prob", 
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
      resp.name <- as.character(object$info$call$formula[2])
      pm <-  distfit(object$data[,resp.name], family = object$info$family, weights = nw[,1], vcov = FALSE)
      pred.val1 <- predict(pm, type = "response")
        
      pred.val <- data.frame(idx = 1:nrow(nd))
      pred.val[1,] <- pred.val1
      
      if(nrow(nd)>=2){
        for(i in 2:nrow(nd)){
          nwi <- nw[,i]
          # personalized model
          pm <-  distfit(object$data[,resp.name], family = object$info$family, weights = nwi, vcov = FALSE)
          pred.val[i,] <- predict(pm, type = "response")
        }
      }
      pred <- pred.val
      colnames(pred) <- c("(fitted.response)")
      return(pred)
    } else return(object$fitted$`(fitted.response)`)
  }
  
  if(type == "parameter"){
    if(!is.null(newdata)) {
      
      # get weights for new data
      nw <- predict.cforest(object, newdata = nd, type = "weights", OOB = FALSE)
      
      # calculate prediction for the first observation before the loop in order to get the number of parameters
      resp.name <- as.character(object$info$call$formula[2])
      pm <-  distfit(object$data[,resp.name], family = object$info$family, weights = nw[,1], vcov = FALSE)
      pred.par1 <- coef(pm, type = "parameter")

      pred.par <- data.frame(matrix(0, nrow = nrow(nd), ncol = length(pred.par1)))
      pred.par[1,] <- pred.par1
      
      if(nrow(nd)>=2){
        for(i in 2:nrow(nd)){
          nwi <- nw[,i]
          # personalized model
          pm <-  distfit(object$data[,resp.name], family = object$info$family, weights = nwi, vcov = FALSE)
          pred.par[i,] <- coef(pm, type = "parameter")
        }
      }
      pred <- pred.par
      colnames(pred) <- c(names(coef(pm, type = "parameter")))
      return(pred)
    } else return(object$fitted.par)
  }
}



logLik.distforest <- function(object, newdata = NULL) {
  if(is.null(newdata)) {
    return(structure(object$loglik, df = object$npar, class = "logLik"))
  } else {
    ll <- 0
    pred.par <- predict(object, newdata = newdata, type = "parameter")
    distlist <- if(inherits(object$info$family, "gamlss.family")) make_dist_list(object$info$family) else object$info$family
    for(i in 1:(nrow(newdata))){
      par <- pred.par[i,]
      eta <-  as.numeric(distlist$linkfun(par))
      ll <- ll + distlist$ddist(newdata[i,paste(object$info$formula[[2]])], eta = eta,  log=TRUE)
    }
      
    return(ll)
  }
}










