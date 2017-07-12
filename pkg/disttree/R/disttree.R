## high-level convenience interface to mob() and ctree()
# FIX ME: default settings for family, decorrelate only necesary for type.tree == "ctree"
# FIX ME: 'starting weights' in trees?
disttree <- function(formula, data, na.action, cluster, family = NO(),
                     type.tree = "mob", decorrelate = "none", offset,
                     cens = "none", censpoint = NULL, weights = NULL,
                     control = mob_control(), ocontrol = list(), ...)
{
  ## keep call
  cl <- match.call(expand.dots = TRUE)
  if(missing(data)) data <- environment(formula)
  ## for gamlss.family function: turn into gamlss. family object
  if(is.function(family)) family <- family()
  
  resp.name <- as.character(formula[2])
  np <- if(inherits(family, "gamlss.family")) family$nopar else length(family$link)
  
  # check input arguments
  type.tree <- match.arg(type.tree, c("mob", "ctree"))
  # if(!(type.tree %in% c("mob", "ctree"))) stop("unknown argument for type.tree (can only be mob or ctree)")
  if(!(decorrelate) %in% c("none", "opg", "vcov")) stop("unknown argument for decorrelate (can only be none, opg or vcov)")
  
  m <- match.call(expand.dots = FALSE)
  # m$drop.unused.levels <- TRUE
  
  ## FIX ME: doesn't work for ctree
  ## formula
  #oformula <- as.formula(formula)
  #formula <- as.Formula(formula)
  #if(length(formula)[2L]  >= 2L) {
  #  stop("formula can have only one RHS consisting of the partitioning variables") 
  #}
  #if(length(formula)[1L]  >= 2L) {
  #  stop("formula can only have one LHS consisting of the response variable")
  #}
  #m$formula <- formula
  
  if(type.tree == "mob") {
    
    # select arguments for mob and put them in the right order
    mo <- match(c("formula", "data", "subset", "na.action", "weights", "offset", "cluster", "control"), names(m), 0L)
    m <- m[c(1L, mo)]
    
    ## glue code for calling distfit() with given family in mob()
    dist_family_fit <- function(y, x = NULL, start = NULL, weights = NULL, offset = NULL,
                                cluster = NULL, vcov = FALSE, estfun = TRUE, 
                                object = FALSE, type.hessian = "analytic", ...)
    {
      if(!(is.null(x) || NCOL(x) == 0L)) warning("x not used")
      if(!is.null(offset)) warning("offset not used")
      
      model <- distfit(y, family = family, weights = weights, start = start,
                      vcov = vcov, estfun = estfun, type.hessian = type.hessian,
                      cens = cens, censpoint = censpoint, ocontrol = ocontrol, ...)
      
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
        estfun = estfun, 
        object = if(object) model else NULL
      )
      return(rval)
    }
    
    ## call mob
    m$fit <- dist_family_fit
    # m$family <- m$censpoint <- m$cens <- NULL
    # m$family <- m$ocontrol <- NULL
    # for(n in names(ocontrol)) m[[n]] <- ocontrol[[n]]
    if("..." %in% names(m)) m[["..."]] <- NULL
    # if("type.tree" %in% names(m)) m[["type.tree"]] <- NULL
    m[[1L]] <- as.name("mob")
    # m[[1L]] <- as.name("partykitR1::mob")
    rval <- eval(m, parent.frame())
    
    rval$fitted$`(weights)` <- if(length(weights)>0) weights else rep.int(1, nrow(data)) 
    rval$fitted$`(response)` <- data[,paste(formula[[2]])]
    rval$fitted$`(fitted.response)` <- predict(rval, type = "response")
    rval$coefficients <- coef(rval)    # rval is returned from mob -> no type argument needed
    rval$loglik <- logLik(rval)
  }
  
  
  if(type.tree == "ctree") {
    
    # select arguments for ctree and put them in the right order
    mo <- match(c("formula", "data", "weights", "subset", "offset", "cluster", "na.action", "control"), names(m), 0L)
    m <- m[c(1L, mo)]
    
    ## wrapper function to apply distfit in ctree
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
                         estfun = estfun, cens = cens, censpoint = censpoint, ocontrol = ocontrol, ...)

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
          
          estfun <- matrix(0, ncol = ncol(ef), nrow = nrow(data)) 
          estfun[subset,] <- ef
          if(!(is.null(weights) || (length(weights)==0L))) estfun <- estfun / weights # estfun has to be unweighted for ctree
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
    
    ## call ctree
    m$ytrafo <- ytrafo
    # for(n in names(ocontrol)) m[[n]] <- ocontrol[[n]]
    if("..." %in% names(m)) m[["..."]] <- NULL
    #if("type.tree" %in% names(m)) m[["type.tree"]] <- NULL
    m[[1L]] <- as.name("ctree")
    rval <- eval(m, parent.frame())

    # number of terminal nodes
    n_tn <- width(rval)
    # predicted terminal nodes for the given data
    pred_tn <- predict(rval, type = "node")
    # ids of terminal nodes
    id_tn <- as.vector(unique(pred_tn))
    
    if(is.null(weights) || (length(weights)==0L)) weights <- numeric(nrow(data)) + 1
    
    # get coefficients for terminal nodes:
    # first iteration out of loop:
    model1 <- distfit(y = data[(id_tn[1]==pred_tn),resp.name], family = family, weights = weights[(id_tn[1]==pred_tn)], start = NULL,
                     vcov = FALSE, type.hessian = "analytic", 
                     estfun = FALSE, cens = cens, censpoint = censpoint, ocontrol = ocontrol, ...)
    coefficients_par <- matrix(nrow = n_tn, ncol = length(model1$par))
    # coefficients_eta <- matrix(nrow = n_tn, ncol = length(model1$eta)) 
    colnames(coefficients_par) <- names(model1$par)
    # colnames(coefficients_eta) <- names(model1$eta)
    rownames(coefficients_par) <- as.character(id_tn)
    # rownames(coefficients_eta) <- as.character(id_tn)
    
    coefficients_par[1,] <- model1$par
    # coefficients_eta[1,] <- model1$eta
    
    loglik <- sum(model1$ddist(data[(id_tn[1]==pred_tn),resp.name], log = TRUE))
    
    if(n_tn>1){
      for(i in (2:n_tn)){
        model <- distfit(y = data[(id_tn[i]==pred_tn),resp.name], family = family, weights = weights[(id_tn[i]==pred_tn)], start = NULL,
                         vcov = FALSE, type.hessian = "analytic", 
                         estfun = FALSE, cens = cens, censpoint = censpoint, ocontrol = ocontrol, ...)
        coefficients_par[i,] <- model$par
        # coefficients_eta[i,] <- model$eta
        loglik <- loglik + sum(model$ddist(data[(id_tn[i]==pred_tn),resp.name], log = TRUE))
      }
    }
    
    rval$coefficients <- coefficients_par
    rval$fitted$`(fitted.response)` <- predict(rval, type = "response")
    rval$loglik <- loglik
  }
  
  
  
  ## extend class and keep original call/family/control
  rval$info$call <- cl
  rval$info$family <- family   # FIX ME: family list is only generated within distfit -> not always returned
  #rval$info$family <-  family$family.name
  #rval$info$familylist <- family
  rval$info$ocontrol <- ocontrol
  rval$info$formula <- rval$info$call$formula
  rval$info$censpoint <- censpoint
  rval$info$cens <- cens
  
  groupcoef <- rval$coefficients
  if(!(is.null(groupcoef))){
    if(is.vector(groupcoef)) {
      groupcoef <- t(as.matrix(groupcoef))
      rownames(groupcoef) <- 1
    }
    rval$fitted.par <- groupcoef[paste(rval$fitted[,1]),]
    rownames(rval$fitted.par) <- c(1: (length(rval$fitted.par[,1])))
    rval$fitted.par <- as.data.frame(rval$fitted.par)
  }
  
  
  class(rval) <- c("disttree", class(rval))
  return(rval)
}


## methods
print.disttree <- function(x, title = NULL, objfun = "negative log-likelihood", ...)
{
  familyname <- if(inherits(x$info$family, "gamlss.family")) {
    paste(x$info$family[[1]][2], "Distribution")
  } else {x$info$family$family.name}
  if(is.null(title)) title <- sprintf("Distributional regression tree (%s)", familyname)
  partykit::print.modelparty(x, title = title, objfun = objfun, ...)
}



predict.disttree <- function (object, newdata = NULL, type = c("parameter", "node", "response"), OOB = FALSE, ...) 
{
  if((type == "node") || (type == "response")) {
    # if mob was applied
    if(inherits(object, "modelparty")){
      return(predict.modelparty(object = object, newdata = newdata, type = type, OOB = OOB, ...))
    }
    # if ctree was applied
    if(inherits(object, "constparty")){
      return(partykit::predict.party(object = object, newdata = newdata, type = type, OOB = OOB, ...))
    }
  }
  if(type == "parameter") {
    if(inherits(object, "constparty")) pred.subgroup <- partykit::predict.party(object, newdata =  newdata, type = "node")
    if(inherits(object, "modelparty")) pred.subgroup <- partykit::predict.modelparty(object, newdata =  newdata, type = "node")
    groupcoef <- coef(object)
    if(is.vector(groupcoef)) {
      groupcoef <- t(as.data.frame(groupcoef))
      rownames(groupcoef) <- 1
    }
    pred.par <- groupcoef[paste(pred.subgroup),]
    rownames(pred.par) <- c(1: (length(pred.par[,1])))
    pred.par <- as.data.frame(pred.par)
    return(pred.par)
  }
}

  


coef.disttree <- function(object, ...){
  object$coefficients
}


logLik.disttree <- function(object, newdata = NULL, ...) {
  if(is.null(newdata)) {
    return(structure(object$loglik, df = ncol(coef(object))*width(object) + width(object)-1 , class = "logLik"))
  } else {
    ll <- 0
    # predicted nodes for the new dataset
    pred.node <- predict(object, newdata = newdata, type = "node")
    # coefficients in the terminal nodes
    coef_tn <- coef(object)
    # number of terminal nodes
    n_tn <- width(object) # <- nrow(coef_tn)
    # id of terminal nodes
    id_tn <- rownames(coef_tn)
    # get link fun and ddist from distribution list
    distlist <- if(inherits(object$info$family, "gamlss.family")) make_dist_list(object$info$family) else object$info$family
    linkfun <- distlist$linkfun
    ddist <- distlist$ddist
    
    
    if(inherits(object$info$family, "gamlss.family") && ("censored" %in% strsplit(object$info$family[[1]], " ")[[2]])) {
      cens <- object$info$cens
      censpoint <- object$info$censpoint
      for(i in 1:n_tn){
        par <- coef_tn[i,]
        eta <-  as.numeric(linkfun(par))
        # response variable of the observations that end up in this terminal node
        nobs_tn <- newdata[pred.node == id_tn[i], paste(object$info$formula[[2]])]
        if(!survival::is.Surv(nobs_tn)) {
          if(cens == "left") ll <- ll + ddist(survival::Surv(nobs_tn, nobs_tn > censpoint, type = "left"), eta = eta, log = TRUE, sum = TRUE)
          if(cens == "right") ll <- ll + ddist(survival::Surv(nobs_tn, nobs_tn < censpoint, type = "right"), eta = eta, log = TRUE, sum = TRUE)
          ## FIX ME: interval censored
        } else ll <- ll + ddist(nobs_tn, eta = eta,  log=TRUE, sum = TRUE)
      }
    } else {
      for(i in 1:n_tn){
        par <- coef_tn[i,]
        eta <-  as.numeric(linkfun(par))
        # response variable of the observations that end up in this terminal node
        nobs_tn <- newdata[pred.node == id_tn[i], paste(object$info$formula[[2]])]
        ll <- ll + ddist(nobs_tn, eta = eta,  log=TRUE, sum = TRUE)
      }
    }
    return(structure(ll, df = ncol(coef(object))*width(object) + width(object)-1 , class = "logLik"))
  }
}



## predict.disttree <- function(object, newdata = NULL,
##   type = c("worth", "rank", "best", "node"), ...)
## {
##   ## type of prediction
##   type <- match.arg(type)
##   
##   ## nodes can be handled directly
##   if(type == "node") return(partykit::predict.modelparty(object, newdata = newdata, type = "node", ...))
##   
##   ## get default newdata otherwise
##   if(is.null(newdata)) newdata <- model.frame(object)
##   
##   pred <- switch(type,
##     "worth" = worth,
##     "rank" = function(obj, ...) rank(-worth(obj)),
##     "best" = function(obj, ...) {
##       wrth <- worth(obj)
##       factor(names(wrth)[which.max(wrth)], levels = names(wrth))
##     }
##   )
##   partykit::predict.modelparty(object, newdata = newdata, type = pred, ...)
## }


## FIX: adapt for disttree_mob
## plot.disttree <- function(x, terminal_panel = node_histogram,
##   tp_args = list(...), tnex = NULL, drop_terminal = NULL, ...)
## {
##   if(is.null(tnex)) tnex <- if(is.null(terminal_panel)) 1L else 2L
##   if(is.null(drop_terminal)) drop_terminal <- !is.null(terminal_panel)
##   partykit::plot.modelparty(x, terminal_panel = terminal_panel,
##     tp_args = tp_args, tnex = tnex, drop_terminal = drop_terminal, ...)
## }


if(FALSE){
  tr <- disttree(dist ~ speed, data = cars)
  # trc <- disttree(dist ~ speed, data = cars, type.tree = "ctree")
  print(tr)
  
  plot(tr)
  plot(as.constparty(tr))
}


if(FALSE){
  y1 <- rNO(200,1,0.5)
  y2 <- rNO(200,5,2)
  y3 <- rNO(200,50,4)
  y4 <- rNO(200,100,7)
  x1 <- vector(mode = "numeric",length = length(y1)) + 1
  x2 <- vector(mode = "numeric",length = length(y2)) + 2
  x3 <- vector(mode = "numeric",length = length(y3)) + 3
  x4 <- vector(mode = "numeric",length = length(y4)) + 4
  d <- as.data.frame(cbind(c(y1,y2,y3,y4),c(x1,x2,x3,x4)))
  colnames(d) <- c("y","x")
  
  test_mob <- disttree(y~x, data = d, family = dist_list_normal, type.tree = "mob")
  test_ctree <- disttree(y~x, data = d, family = dist_list_normal, type.tree = "ctree", control = ctree_control())
  print(test_mob)
  plot(test_mob)
  plot(as.constparty(test_mob))
  print(test_ctree)
  plot(test_ctree)
}
