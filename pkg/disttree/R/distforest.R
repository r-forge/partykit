distforest <- function(formula, data, na.action = na.pass, cluster, family = NO(), 
                       type.tree = "ctree", decorrelate = "none", offset,
                       cens = "none", censpoint = NULL, weights = NULL,
                       control = ctree_control(teststat = "quad", testtype = "Univ", mincriterion = 0, ...), 
                       ocontrol = list(),
                       ntree = 500L, fit = TRUE, perturb = list(replace = FALSE, fraction = 0.632), fitted.OOB = TRUE,
                       cores = NULL, applyfun = NULL,
                       mtry = ceiling(sqrt(nvar)),
                       ...)
{
  ## keep call
  cl <- match.call(expand.dots = TRUE)
  if(missing(data)) data <- environment(formula)
  ## for gamlss.family function: turn into gamlss. family object
  if(is.function(family)) family <- family()
  
  resp.name <- as.character(formula[2])
  np <- if(inherits(family, "gamlss.family")) family$nopar else length(family$link)
  nvar <- (length(formula[[3]]) + 1)/2    # number of variables (RHS minus nr of arithmetic symbols)
  
  # check input arguments
  type.tree <- match.arg(type.tree, c("mob", "ctree"))
  # if(!(type.tree %in% c("mob", "ctree"))) stop("unknown argument for type.tree (can only be mob or ctree)")
  if(!(decorrelate) %in% c("none", "opg", "vcov")) stop("unknown argument for decorrelate (can only be none, opg or vcov)")
  
  m <- match.call(expand.dots = FALSE)
  
  if(type.tree == "ctree") {
    
    # select arguments for cforest and put them in the right order
    mo <- match(c("formula", "data", "weights", "subset", "offset", "cluster", 
                  "strata", "na.action", "control", "ntree", "perturb", "mtry", "applyfun", "cores"), names(m), 0L)
    m <- m[c(1L, mo)]
    
    ## wrapper function to apply distfit in cforest
    ytrafo <- function(formula, data, weights = NULL, cluster = cluster, ctrl = control) {
      
      if(!(is.null(cluster))) stop("FIX: cluster ignored by trafo-function")
      if(!(is.numeric(formula[[3]]))) {
        #print(formula)
        stop("covariates can only be used as splitting variables (formula has to be of type y~1|x or y~0|x)")
      }
      
      # decorrelate <- if(is.null(ctrl$decorrelate)) "none" else ctrl$decorrelate  # FIX ME: include in ctrl?
      
      modelscores_decor <- function(subset, estfun = TRUE, object = TRUE, info = NULL) {
        
        ys <- data[subset,resp.name]
        subweights <- if(is.null(weights) || (length(weights)==0L)) weights else weights[subset] ## FIX ME: scores with or without weights?
        # start <- if(!(is.null(info$coefficients))) info$coefficients else NULL
        start <- info$coefficients
        
        model <- distfit(ys, family = family, weights = subweights, start = start,
                         vcov = (decorrelate == "vcov"), type.hessian = "analytic", 
                         estfun = estfun, cens = cens, censpoint = censpoint, ...)
        
        
        ef <- as.matrix(model$estfun)
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
                    converged = model$converged  # FIX ME: warnings if distfit does not converge
        )
        return(ret)
      }
      return(modelscores_decor)
    }    
    
    ## call cforest
    m$ytrafo <- ytrafo
    #for(n in names(ocontrol)) m[[n]] <- ocontrol[[n]]
    if("..." %in% names(m)) m[["..."]] <- NULL
    #if("type.tree" %in% names(m)) m[["type.tree"]] <- NULL
    m[[1L]] <- as.name("cforest")
    rval <- eval(m, parent.frame())
  }
  
  
  # first approach to building a forest using mob
  ## FIX ME: applies mob and extracts nodes slot
  if(type.tree == "mob") {
    
    if(!("control" %in% names(cl))) {
      control <- mob_control()
      warning("for type.tree = 'mob' the control argument is by default set to mob_control()") 
    }
    
    cl$na.action <- na.action
    frame <- parent.frame()

    ## FIX ME: weights?
    # weights <- rep.int(1, nrow(data))
    
    cl$weights <- NULL ### NOTE: trees are unweighted, weights enter sampling!
    
    ## glue code for calling distfit() with given family in mob()
    dist_family_fit <- function(y, x = NULL, start = NULL, weights = NULL, offset = NULL,
                                cluster = NULL, vcov = FALSE, estfun = TRUE, 
                                object = FALSE, type.hessian = "analytic", ...)
    {
      if(!(is.null(x) || NCOL(x) == 0L)) warning("x not used")
      if(!is.null(offset)) warning("offset not used")
      
      mod <- distfit(y, family = family, weights = weights, start = start,
                     vcov = vcov, estfun = estfun, type.hessian = type.hessian,
                     cens = cens, censpoint = censpoint, ...)
      
      rval <- list(
        coefficients = mod$par,
        objfun = - mod$loglik,      #FIX ME: before: return minimized function, now: maximized?
        estfun = if(estfun) mod$estfun else NULL,   ## rval$estfun contains the scores of the positive loglik 
        object = if(object) mod else NULL
      )
      return(rval)
    }
    
    tree <- mob(formula, data = data,
                fit= dist_family_fit, control = control)
    
    tree$mf <- model.frame(formula, data)
    strata <- tree$mf[["(strata)"]]
    if (!is.null(strata)) {
      if (!is.factor(strata)) stop("strata is not a single factor")
    }
    
    probw <- NULL
    weights <- model.weights(tree$mf)
    if (!is.null(weights)) {
      probw <- weights / sum(weights)
    } else {
      weights <- integer(0)
    }
    # nvar <- length(tree$partyvars) # doesn't work for modelparty tree
    control$mtry <- mtry
    #control$applyfun <- applyfun
    
    idx <- 1L:nrow(tree$mf)
    if (is.null(strata)) {
      size <- nrow(tree$mf)
      if (!perturb$replace) size <- floor(size * perturb$fraction)
      rw <- replicate(ntree, sample(idx, size = size, replace = perturb$replace, prob = probw),
                      simplify = FALSE)
    } else {
      frac <- if (!perturb$replace) perturb$fraction else 1
      rw <- replicate(ntree, function() 
        do.call("c", tapply(idx, strata, function(i) sample(i, size = length(i) * frac, 
                                                            replace = perturb$replace, prob = probw[i]))))
    }
    
    ## FIX ME: hand over applyfun (via control)?
    # if(is.null(applyfun)) applyfun <- control$applyfun
    # cores <- NULL
    # applyfun <- tree$info$control$applyfun()
    ## apply infrastructure for determining split points
    if (is.null(applyfun)) {
      applyfun <- if(is.null(cores)) {
        lapply  
      } else {
        function(X, FUN, ...)
          parallel::mclapply(X, FUN, ..., mc.cores = cores)
      }
    }
    
    forest <- applyfun(1:ntree, function(b) {
      sdata <- data[rw[[b]],]
      mob(formula, data = sdata, fit= dist_family_fit, control = control)$node
    })
    
    #forest <- list()
    #for(b in 1:ntree){
    #  forest[[b]] <- mob(formula, data = data, subset = rw[[b]], fit= dist_family_fit, control = mob_control())$node
    #}
    
    
    fitted <- data.frame(idx = idx)  
    mf <- model.frame(Formula::Formula(formula), data = tree$mf, na.action = na.pass)
    y <- Formula::model.part(Formula::Formula(formula), data = mf, lhs = 1, rhs = 0)
    if (length(y) == 1) y <- y[[1]]
    fitted[[2]] <- y
    names(fitted)[2] <- "(response)"
    fitted <- fitted[2]
    if (length(weights) > 0)
      fitted[["(weights)"]] <- weights
    
    ### turn subsets in weights (maybe we can avoid this?)
    rw <- lapply(rw, function(x) tabulate(x, nbins = length(idx)))
    
    control$applyfun <- applyfun
    
    ## FIXME: export constparties() from partykit
    rval <- constparties(nodes = forest, data = tree$mf, weights = rw,
                                   fitted = fitted, terms = terms(mf), 
                                   info = list(call = match.call(), control = control))
    
    #rval$fit <- tree$fit
    class(rval) <- c("distforest", class(rval))
  }

  
  
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
      loglik[i,] <- pm$ddist(data[i,resp.name], log = TRUE)
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
  
  rval$info$call <- cl
  rval$info$family <-  family
  #rval$info$family <-  family$family.name
  #rval$info$familylist <- family
  rval$info$npar <- np
  rval$info$formula <- formula
  rval$info$cens <- cens
  rval$info$censpoint
  
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
      #nw <- predict.cforest(object, newdata = nd, type = "weights", OOB = OOB)
      nw <- w
      
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



logLik.distforest <- function(object, newdata = NULL, ...) {
  if(is.null(newdata)) {
    if(!is.null(object$loglik)) return(structure(object$loglik, df = NA, class = "logLik"))
    newdata <- data
  } 
  ll <- 0
  pred.par <- predict(object, newdata = newdata, type = "parameter")
  distlist <- if(inherits(object$info$family, "gamlss.family")) make_dist_list(object$info$family) else object$info$family
  
  if(inherits(object$info$family, "gamlss.family") && ("censored" %in% strsplit(object$info$family[[1]], " ")[[2]])) {
    cens <- object$info$cens
    censpoint <- object$info$censpoint
    for(i in 1:(nrow(newdata))){
      par <- pred.par[i,]
      eta <-  as.numeric(distlist$linkfun(par))
      ydata <- newdata[i,paste(object$info$formula[[2]])]
      if(!survival::is.Surv(ydata)) {
        if(cens == "left") ll <- ll + distlist$ddist(survival::Surv(ydata, ydata > censpoint, type = "left"), eta = eta, log = TRUE)
        if(cens == "right") ll <- ll + distlist$ddist(survival::Surv(ydata, ydata < censpoint, type = "right"), eta = eta, log = TRUE)
        ## FIX ME: interval censored
      } else ll <- ll + distlist$ddist(ydata, eta = eta,  log=TRUE)
    }
  } else {
    for(i in 1:(nrow(newdata))){
      par <- pred.par[i,]
      eta <-  as.numeric(distlist$linkfun(par))
      ll <- ll + distlist$ddist(newdata[i,paste(object$info$formula[[2]])], eta = eta,  log=TRUE)
    }
  }
  return(structure(ll, df = NA, class = "logLik"))
}

## -----------------------------------------------------------------------------

## FIXME: currently copied from partykit:::constparties
## -> check with TH whether this can be documented/exported from partykit
constparties <- function(nodes, data, weights, fitted = NULL, terms = NULL, info = NULL) {

    stopifnot(all(sapply(nodes, function(x) inherits(x, "partynode"))))
    stopifnot(inherits(data, "data.frame"))
    stopifnot(inherits(weights, "list"))

    if(!is.null(fitted)) {
        stopifnot(inherits(fitted, "data.frame"))
        stopifnot(nrow(data) == 0L | nrow(data) == nrow(fitted))
        if (nrow(data) == 0L)
            stopifnot("(response)" %in% names(fitted))
    } else {
        stopifnot(nrow(data) > 0L)
        stopifnot(!is.null(terms))
        fitted <- data.frame("(response)" = model.response(model.frame(terms, data = data, 
                                                                       na.action = na.pass)),
                             check.names = FALSE)
    }

    ret <- list(nodes = nodes, data = data, weights = weights, fitted = fitted)
    class(ret) <- c("constparties", "parties")

    if(!is.null(terms)) {
        stopifnot(inherits(terms, "terms"))
        ret$terms <- terms
    }

    if (!is.null(info))
        ret$info <- info

    ret
}
