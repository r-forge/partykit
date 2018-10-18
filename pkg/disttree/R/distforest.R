distforest <- function(formula, data, na.action = na.pass, cluster, family = NO(), bd = NULL,
                       type.tree = "ctree", decorrelate = "none", offset,
                       censtype = "none", censpoint = NULL, weights = NULL,
                       control = partykit::ctree_control(teststat = "quad", testtype = "Univ", 
                                                         mincriterion = 0, ...), 
                       ocontrol = list(), type.hessian = c("checklist", "analytic", "numeric"),
                       ntree = 500L, fit = TRUE, 
                       perturb = list(replace = FALSE, fraction = 0.632), 
                       fitted.OOB = TRUE, cores = NULL, applyfun = NULL,
                       mtry = ceiling(sqrt(nvar)),
                       ...)
{
  ## keep call
  cl <- match.call(expand.dots = TRUE)
  if(missing(data)) data <- environment(formula)
  
  
  ## prepare family:
  # check format of the family input and if necessary transform it to the required familiy list
  # family input can be of one of the following formats:
  # - gamlss.family object
  # - gamlss.family function
  # - character string with the name of a gamlss.family object
  # - function generating a list with the required information about the distribution
  # - character string with the name of a function generating a list with the required information about the distribution
  # - list with the required information about the distribution
  # - character string with the name of a distribution for which a list generating function is provided in disttree
  {
    if(is.character(family)) {
      getfamily <- try(getAnywhere(paste("dist", family, sep = "_")), silent = TRUE)
      if(length(getfamily$objs) == 0L) getfamily <- try(getAnywhere(family), silent = TRUE)
      if(length(getfamily$objs) == 0L) {
        stop("unknown 'family' specification")
      } else {
        gamlssobj <- ("gamlss.dist" %in% unlist(strsplit(getfamily$where[1], split = ":")))
        family <- getfamily[[2]][[1]]() #first found is chosen 
        family$gamlssobj <- gamlssobj
      }
      #if(!(inherits(family, "try-error")))family <- family[[2]]$`package:disttree`()    
      # FIX ME: better selection of dist function
    }
    
    # if family is a gamlss family object or gamlss family function
    if(is.function(family)) family <- family()
    if(inherits(family, "gamlss.family")) family <- disttree::make_dist_list(family, bd = bd)
    
    if(!is.list(family)) stop ("unknown family specification")
    if(!(all(c("ddist", "sdist", "link", "linkfun", "linkinv", "mle", "startfun") %in% names(family)))) stop("family needs to specify a list with ...")
    # linkinvdr only used in the method vcov for type = "parameter"
    
  }
  
  np <- length(family$link)
  nvar <- (length(formula[[3]]) + 1)/2    # number of variables (RHS minus nr of arithmetic symbols)
  
  # check input arguments
  type.tree <- match.arg(type.tree, c("mob", "ctree"))
  # if(!(type.tree %in% c("mob", "ctree"))) stop("unknown argument for type.tree (can only be mob or ctree)")
  if(!(decorrelate) %in% c("none", "opg", "vcov")) stop("unknown argument for decorrelate (can only be none, opg or vcov)")
  
  m <- match.call(expand.dots = FALSE)
  #m$drop.unused.levels <- TRUE
  
  ## formula
  oformula <- as.formula(formula)
  formula <- Formula::as.Formula(formula)
  if(length(formula)[2L] > 1L) {
    formula <- Formula::Formula(formula(formula, rhs = 2L))  
    ## FIX ME: if rhs has more than 1 element it is here assumed that partitioning variables are handed over on 2nd slot
    warning("formula must not have more than one RHS parts (only partitioning variables allowed)")
  }
  
  m$formula <- formula
  
  
  if(type.tree == "ctree") {
    
    # select arguments for cforest and put them in the right order
    mnames <- match(c("formula", "data", "weights", "subset", "offset", "cluster", 
                      "strata", "na.action", "control", "ntree", "perturb", "mtry", "applyfun", "cores"), names(m), 0L)
    m <- m[c(1L, mnames)]
    
    
    ## wrapper function to apply distfit in cforest
    ytrafo <- function(data, weights = NULL, control) {
      
      Y <- model.frame(data, yxonly = TRUE)
      if(dim(Y)[2] > 1) stop("response variable has to be univariate") 
      Y <- Y[,1]
      
      modelscores_decor <- function(subset, weights, estfun = TRUE, object = TRUE, info = NULL) {
        
        ys <- Y[subset]
        subweights <- if(is.null(weights) || (length(weights)==0L)) weights else weights[subset]
        # start <- if(!(is.null(info$coefficients))) info$coefficients else NULL
        start <- info$coefficients
        
        model <- disttree::distfit(ys, family = family, weights = subweights, start = start,
                                   vcov = (decorrelate == "vcov"), type.hessian = type.hessian, 
                                   estfun = estfun, censtype = censtype, censpoint = censpoint,
                                   ocontrol = ocontrol, ...)
        
        if(estfun) {
          ef <- as.matrix(model$estfun)
          
          if(decorrelate != "none") {
            n <- NROW(ef)
            ef <- ef/sqrt(n)
            
            vcov <- if(decorrelate == "vcov") {
              vcov(model, type = "link") * n
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
          
          estfun <- matrix(0, ncol = ncol(ef), nrow = nrow(data$data)) 
          estfun[subset,] <- ef
          ### now automatically if(!(is.null(weights) || (length(weights)==0L))) estfun <- estfun / weights # estfun has to be unweighted for ctree
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
    m[[1L]] <- quote(partykit::cforest)
    rval <- eval(m, parent.frame())
  }
  
  
  # first approach to building a forest using mob
  ## FIX ME: applies mob and extracts nodes slot
  if(type.tree == "mob") {
    
    if(!("control" %in% names(cl))) {
      control <- mob_control()
      ## FIX ME: change control arguments such as mincriterion etc for forest (compare cforest)
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
                                object = FALSE, ...)
    {
      if(!(is.null(x) || NCOL(x) == 0L)) warning("x not used")
      if(!is.null(offset)) warning("offset not used")
      
      model <- disttree::distfit(y, family = family, weights = weights, start = start,
                                 vcov = vcov, estfun = estfun, type.hessian = type.hessian,
                                 censtype = censtype, censpoint = censpoint,
                                 ocontrol = ocontrol, ...)
      
      ef <- NULL
      if(estfun) {
        ef <- as.matrix(model$estfun)
        
        if(decorrelate != "none") {
          n <- NROW(ef)
          ef <- ef/sqrt(n)
          
          vcov <- if(decorrelate == "vcov") {
            vcov(model, type = "link") * n
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
      }
      estfun <- ef
      
      rval <- list(
        coefficients = model$par,
        objfun = - model$loglik,      
        estfun = estfun,   ## rval$estfun contains the scores of the positive loglik 
        object = if(object) model else NULL
      )
      return(rval)
    }
    
    tree <- partykit::mob(formula, data = data,
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
      partykit::mob(formula, data = sdata, fit= dist_family_fit, control = control)$node
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
    rval <- partykit:::constparties(nodes = forest, data = tree$mf, weights = rw,
                                    fitted = fitted, terms = terms(mf), 
                                    info = list(call = match.call(), control = control))
    
    #rval$fit <- tree$fit
    class(rval) <- c("distforest", class(rval))
  }
  
  
  
  if(fit) {
    ### calculate (fitted value,) fitted distribution parameters, loglikelihood (and log scores) for every observation
    # fitted <- data.frame(idx = 1:nrow(data))
    fitted.par <- data.frame(matrix(0, nrow = nrow(data), ncol = np))
    loglik <- data.frame(idx = 1:nrow(data))
    # logscore <- data.frame(matrix(0, nrow = nrow(data), ncol = np))
    
    # extract weights
    w <- partykit::predict.cforest(rval, type = "weights", OOB = fitted.OOB)
    
    Y <- rval$fitted$`(response)`
    
    for(i in 1:nrow(data)){
      wi <- w[,i]
      # personalized model for observation data[i,]
      pm <-  disttree::distfit(Y, family = family, weights = wi, vcov = FALSE, censtype = censtype, censpoint = censpoint, ocontrol = ocontrol, ...)
      # fitted[i,] <- predict(pm, type = "response")
      fitted.par[i,] <- coef(pm, type = "parameter")
      loglik[i,] <- if(is.function(pm$ddist)) pm$ddist(Y[i], log = TRUE) else NA
      # logscore[i,] <- pm$family$sdist(Y[i], eta = coef(pm, type = "link"), sum = FALSE)
    }
    
    if(is.null(weights) || (length(weights)==0L || is.function(weights))) weights <- numeric(nrow(data)) + 1
    rval$fitted$`(weights)` <- weights
    # rval$fitted$`(fitted.response)` <- fitted
    
    
    names(fitted.par) <- names(coef(pm, type = "parameter"))
    # names(logscore) <- names(coef(pm, type = "parameter"))
    rval$fitted.par <- fitted.par
    
    rval$loglik <- sum(loglik)
    # rval$logscore <- logscore
  }
  
  rval$info$call <- cl
  rval$info$family <-  family
  rval$info$npar <- np
  rval$info$formula <- formula
  rval$info$censtype <- censtype
  rval$info$censpoint <- censpoint
  #rval$data <- mf
  
  class(rval) <- c("distforest", class(rval))
  return(rval)
}




###################
# methods for class 'distforest'

predict.distforest <- function (object, newdata = NULL, 
                                type = c("parameter", "response", "prob",  "weights", "node"),  ## FIX ME: check type = "prob" ?
                                OOB = FALSE, FUN = NULL, simplify = TRUE, scale = TRUE, ...) 
{
  
  # per default 'type' is set to 'parameter'
  if(length(type)>1) type <- type[1]
  
  # for type = "prob" the argument scale has to set to FALSE
  if(type == "prob" & scale) {
    warning("for type 'prob' the weights have to be integers, i.e. the argument 'scale' has to be set to FALSE")
    scale <- FALSE
  }
  
  if(type %in% c("prob",  "weights", "node"))
    return(partykit::predict.cforest(object = object, newdata = newdata,
                                     type = type, OOB = OOB, FUN = FUN, 
                                     simplify = simplify, scale = scale, ...))
  
  if(type == "parameter" & is.null(newdata)) return(object$fitted.par)
  
  # get weights
  w <- partykit::predict.cforest(object = object, newdata = newdata, 
                                 type = "weights", OOB = OOB, FUN = FUN, 
                                 simplify = simplify, scale = scale, ...)
  
  responses <- object$fitted[["(response)"]]
  
  # number of (possibly new) observations to make predictions for
  n <- if(is.null(newdata)) nrow(object$data) else nrow(newdata)   # n <- ncol(w)
  
  # for type "parameter" the number of columns is set to the number of parameters
  # for type "response" it is set to 1
  pred <- if(type == "parameter") {
    data.frame(matrix(0, nrow = n, ncol = ncol(object$fitted.par))) 
  } else {
    data.frame(matrix(0, nrow = n, ncol = 1))
  }      
  
  for(i in 1:ncol(w)){
    wi <- w[,i]       ## FIX ME: what to do if all weights for observation i are 0?
    # personalized model
    pm <-  disttree::distfit(responses, family = object$info$family, weights = wi, vcov = FALSE, 
                             ocontrol = object$call$ocontrol,
                             censtype = object$info$censtype, censpoint = object$info$censpoint)
    pred[i,] <- predict(pm, type = type)
  }
  
  colnames(pred) <- if(type == "parameter") c(names(coef(pm, type = "parameter"))) else c("(fitted.response)")
  return(pred)
}




logLik.distforest <- function(object, newdata = NULL, weights = NULL, ...) {
  if(is.null(newdata)) {
    if(!is.null(weights)) stop("for weighted loglikelihood hand over data as newdata")
    if(!is.null(object$loglik)) return(structure(object$loglik, df = NA, class = "logLik"))
    newdata <- object$data
  } 
  
  #store modelframe for newdata
  mf <- model.frame(formula = object$info$formula, data = newdata)
  
  # extract response
  Y <- model.response(mf, "numeric")
  
  # check weights
  if(is.null(weights)) weights <- rep.int(1, NROW(newdata)) else stopifnot(NROW(newdata) == nrow(weights))
  
  ll <- 0
  pred.par <- predict(object, newdata = newdata, type = "parameter")
  distlist <- if(inherits(object$info$family, "gamlss.family")) disttree::make_dist_list(object$info$family) else object$info$family
  ## FIXME: get family for all types of input for 'family'
  
  if(object$info$family$gamlssobj && object$info$family$censored) {
    censtype <- object$info$censtype
    censpoint <- object$info$censpoint
    for(i in 1:(nrow(newdata))){
      par <- pred.par[i,]
      eta <-  as.numeric(distlist$linkfun(par))
      weight <- weights[i]
      ydata <- Y
      if(!survival::is.Surv(ydata)) {
        if(censtype == "left") ll <- ll + weight * distlist$ddist(survival::Surv(ydata, ydata > censpoint, type = "left"), eta = eta, log = TRUE)
        if(censtype == "right") ll <- ll + weight * distlist$ddist(survival::Surv(ydata, ydata < censpoint, type = "right"), eta = eta, log = TRUE)
        ## FIX ME: interval censored
      } else ll <- ll + distlist$ddist(ydata, eta = eta,  log=TRUE)
    }
  } else {
    for(i in 1:(nrow(newdata))){
      par <- pred.par[i,]
      eta <-  as.numeric(distlist$linkfun(par))
      weight <- weights[i]
      ll <- ll + weight * distlist$ddist(Y[i], eta = eta,  log=TRUE)
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



## FIXME: error handling: if no split found -> start for next ieteration includes NA
# weights are all zero
# simlist <- sim_fun(fun = fun, kappa.start = 1, nsteps = 4, stepsize = 1,
#+                 ntree = 10, nobs = 100, nrep = 2, type.gfun = "nonparametric",  
#+                 covariate = c("x1","x2"), plot = "none", g.interact = FALSE, type.tree = "ctree",
#+                 ctree_minbucket = 20, ctree_mincrit = 0.2,
#+                 cforest_minbucket = 20, cforest_mincrit = 0)
