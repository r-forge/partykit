## high-level convenience interface to mob() and ctree()
# FIX ME: default settings for family, decorrelate only necesary for type.tree == "ctree"
# FIX ME: 'starting weights' in trees?
disttree <- function(formula, data, na.action, cluster, family = NO(), bd = NULL,
                     type.tree = "mob", decorrelate = "none", offset,
                     censtype = "none", censpoint = NULL, weights = NULL,
                     control = mob_control(), ocontrol = list(), ...)
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
    if(inherits(family, "gamlss.family")) family <- make_dist_list(family, bd = bd)
    
    if(!is.list(family)) stop ("unknown family specification")
    if(!(all(c("ddist", "sdist", "link", "linkfun", "linkinv", "mle", "startfun") %in% names(family)))) stop("family needs to specify a list with ...")
    # linkinvdr only used in the method vcov for type = "parameter"
  }
  
  np <- length(family$link)
  
  # check input arguments
  type.tree <- match.arg(type.tree, c("mob", "ctree"))
  if(!(type.tree %in% c("mob", "ctree"))) stop("unknown argument for type.tree (can only be mob or ctree)")
  if(!(decorrelate) %in% c("none", "opg", "vcov")) stop("unknown argument for decorrelate (can only be none, opg or vcov)")
  
  m <- match.call(expand.dots = FALSE)
  #m$drop.unused.levels <- TRUE
  
  ## formula
  oformula <- as.formula(formula)
  formula <- as.Formula(formula)
  if(length(formula)[2L] > 1L) {
    formula <- Formula(formula(formula, rhs = 2L))  
    ## FIX ME: if rhs has more than 1 element it is here assumed that partitioning variables are handed over on 2nd slot
    warning("formula must not have more than one RHS parts (only partitioning variables allowed)")
  }
  
  m$formula <- formula
  

  
  if(type.tree == "mob") {
    
    # select arguments for mob and put them in the right order
    mnames <- match(c("formula", "data", "subset", "na.action", "weights", "offset", "cluster", "control"), names(m), 0L)
    m <- m[c(1L, mnames)]
    
    ## glue code for calling distfit() with given family in mob()
    dist_family_fit <- function(y, x = NULL, start = NULL, weights = NULL, offset = NULL,
                                cluster = NULL, vcov = FALSE, estfun = TRUE, 
                                object = FALSE, type.hessian = "analytic", ...)
    {
      if(!(is.null(x) || NCOL(x) == 0L)) warning("x not used")
      if(!is.null(offset)) warning("offset not used")
      
      model <- distfit(y, family = family, weights = weights, start = start,
                      vcov = vcov, estfun = estfun, type.hessian = type.hessian,
                      censtype= censtype, censpoint = censpoint, ocontrol = ocontrol, ...)
      
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
    # m$family <- m$censpoint <- m$censtype <- NULL
    # m$family <- m$ocontrol <- NULL
    # for(n in names(ocontrol)) m[[n]] <- ocontrol[[n]]
    if("..." %in% names(m)) m[["..."]] <- NULL
    # if("type.tree" %in% names(m)) m[["type.tree"]] <- NULL
    m[[1L]] <- as.name("mob")
    # m[[1L]] <- as.name("partykitR1::mob")
    rval <- eval(m, parent.frame())
    
    rval$fitted$`(weights)` <- if(length(weights)>0) weights else rep.int(1, nrow(data)) 
    rval$fitted$`(response)` <- model.response(rval$data)
    rval$fitted$`(fitted.response)` <- predict(rval, type = "response")
    rval$coefficients <- coef(rval)    # rval is returned from mob -> no type argument needed
    rval$loglik <- logLik(rval)
  }
  
  
  if(type.tree == "ctree") {
    
    # select arguments for ctree and put them in the right order
    mnames <- match(c("formula", "data", "weights", "subset", "offset", "cluster", "na.action", "control"), names(m), 0L)
    m <- m[c(1L, mnames)]
    
    ## wrapper function to apply distfit in ctree
    ytrafo <- function(formula, data, weights = NULL, cluster = cluster, ctrl = control) {
      
      if(!(is.null(cluster))) stop("FIX: cluster ignored by trafo-function")
      
      cl <- match.call()
      if(missing(data)) data <- environment(formula)
      mf <- match.call(expand.dots = FALSE)
      mfnames <- match(c("formula", "data"), names(mf), 0L)
      mf <- mf[c(1L, mfnames)]
      mf$drop.unused.levels <- TRUE
      
      ## formula
      oformula <- as.formula(formula)
      formula <- as.Formula(formula)
      if(length(formula)[2L] > 2L) {
        formula <- Formula(formula(formula, rhs = 2L))  
        ## FIX ME: if rhs has more than 1 element it is here assumed that partitioning variables are handed over on 2nd slot
        warning("formula must not have more than one RHS parts (only partitioning variables allowed)")
      }
      
      mf$formula <- formula
      
      ## evaluate model.frame
      mf[[1L]] <- as.name("model.frame")
      mf <- eval(mf, parent.frame())
      
      ## extract terms, model matrix, response
      #mt <- terms(formula, data = data)
      #mtZ <- terms(formula, data = data, rhs = 1L)
      #attributes(mtZ)$intercept <- 0
      Y <- model.response(mf, "numeric")
      #Z <- model.matrix(mtZ, mf)
      
      # decorrelate <- if(is.null(ctrl$decorrelate)) "none" else ctrl$decorrelate  # FIX ME: include in ctrl?
      
      modelscores_decor <- function(subset, estfun = TRUE, object = TRUE, info = NULL) {

        ys <- Y[subset]
        subweights <- if(is.null(weights) || (length(weights)==0L)) weights else weights[subset] ## FIX ME: scores with or without weights?
        # start <- if(!(is.null(info$coefficients))) info$coefficients else NULL
        start <- info$coefficients
        
        model <- distfit(ys, family = family, weights = subweights, start = start,
                         vcov = (decorrelate == "vcov"), type.hessian = "analytic", 
                         estfun = estfun, censtype = censtype, censpoint = censpoint, ocontrol = ocontrol, ...)

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
    
    ## get coefficients for terminal nodes:
    Y <- rval$fitted$`(response)`
    # first iteration out of loop:
    model1 <- distfit(y = Y[(id_tn[1]==pred_tn)], family = family, weights = weights[(id_tn[1]==pred_tn)], start = NULL,
                     vcov = FALSE, type.hessian = "analytic", 
                     estfun = FALSE, censtype = censtype, censpoint = censpoint, ocontrol = ocontrol, ...)
    coefficients_par <- matrix(nrow = n_tn, ncol = length(model1$par))
    # coefficients_eta <- matrix(nrow = n_tn, ncol = length(model1$eta)) 
    colnames(coefficients_par) <- names(model1$par)
    # colnames(coefficients_eta) <- names(model1$eta)
    rownames(coefficients_par) <- as.character(id_tn)
    # rownames(coefficients_eta) <- as.character(id_tn)
    
    coefficients_par[1,] <- model1$par
    # coefficients_eta[1,] <- model1$eta
    
    loglik <- sum(model1$ddist(Y[(id_tn[1]==pred_tn)], log = TRUE))
    
    if(n_tn>1){
      for(i in (2:n_tn)){
        model <- distfit(y = Y[(id_tn[i]==pred_tn)], family = family, weights = weights[(id_tn[i]==pred_tn)], start = NULL,
                         vcov = FALSE, type.hessian = "analytic", 
                         estfun = FALSE, censtype = censtype, censpoint = censpoint, ocontrol = ocontrol, ...)
        coefficients_par[i,] <- model$par
        # coefficients_eta[i,] <- model$eta
        loglik <- loglik + sum(model$ddist(Y[(id_tn[i]==pred_tn)], log = TRUE))
      }
    }
    
    rval$coefficients <- coefficients_par
    rval$fitted$`(fitted.response)` <- predict(rval, type = "response")
    rval$loglik <- loglik
  }
  
  
  
  ## extend class and keep original call/family/control
  rval$info$call <- cl
  rval$info$family <- family   
  rval$info$ocontrol <- ocontrol
  rval$info$formula <- m$formula
  rval$info$censpoint <- censpoint
  rval$info$censtype <- censtype
  
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
    rownames(pred.par) <- c(1: (NROW(pred.par)))
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
    linkfun <- object$info$family$linkfun
    ddist <- object$info$family$ddist
    
    
    if(object$info$family$gamlssobj && object$info$family$censored) {
      censtype <- object$info$censtype
      censpoint <- object$info$censpoint
      for(i in 1:n_tn){
        par <- coef_tn[i,]
        eta <-  as.numeric(linkfun(par))
        # response variable of the observations that end up in this terminal node
        nobs_tn <- newdata[pred.node == id_tn[i], paste(object$info$formula[[2]])]
        if(length(nobs_tn) > 0L){
          if(!survival::is.Surv(nobs_tn)) {
            if(censtype == "left") ll <- ll + ddist(survival::Surv(nobs_tn, nobs_tn > censpoint, type = "left"), eta = eta, log = TRUE, sum = TRUE)
            if(censtype == "right") ll <- ll + ddist(survival::Surv(nobs_tn, nobs_tn < censpoint, type = "right"), eta = eta, log = TRUE, sum = TRUE)
            ## FIX ME: interval censored
          } else ll <- ll + ddist(nobs_tn, eta = eta,  log=TRUE, sum = TRUE)
        }
      }
    } else {
      for(i in 1:n_tn){
        par <- coef_tn[i,]
        eta <-  as.numeric(linkfun(par))
        # response variable of the observations that end up in this terminal node
        nobs_tn <- newdata[pred.node == id_tn[i], paste(object$info$formula[[2]])]
        if(length(nobs_tn) > 0L) ll <- ll + ddist(nobs_tn, eta = eta,  log=TRUE, sum = TRUE)
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
# argument 'type' can be "density" (FIX), "coef" and "hist" (FIX)
#plot.disttree <- function(x, type = "coef",
#   tp_args = list(...), tnex = NULL, drop_terminal = NULL, ...)
# {
#  if(type == "density"){
#    node_density <- function (tree, xscale = NULL, yscale = NULL, horizontal = FALSE,
#                              main = "", xlab = "", ylab = "Density", id = TRUE, rug = TRUE,
#                              fill = "lightgrey", col = "black", lwd = 0.5, ...) {
#      yobs <- tree$data[,as.character(tree$info$formula[[2]])]
#      ylines <- 1.5
#      if (is.null(xscale)) xscale <- c(-5.1,50)
#      if (is.null(yscale)) yscale <- c(-0.05,0.25)
#      xr <- xscale
#      yr <- yscale
#      
#      if (horizontal) {
#        yyy <- xscale
#        xscale <- yscale
#        yscale <- yyy
#      }
#      
#      rval <- function(node) {
#        yrange <- seq(from = -20, to = 200)/4
#        ydens <- node$info$object$ddist(yrange)
#        
#        top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3, 
#                                                widths = unit(c(ylines, 1, 1), c("lines", "null", "lines")), 
#                                                heights = unit(c(1, 1), c("lines", "null"))), 
#                           width = unit(1, "npc"), 
#                           height = unit(1, "npc") - unit(2, "lines"), 
#                           name = paste("node_density",node$id, sep = ""))
#        pushViewport(top_vp)
#        grid.rect(gp = gpar(fill = "white", col = 0))
#        top <- viewport(layout.pos.col = 2, layout.pos.row = 1)
#        pushViewport(top)
#        mainlab <- paste(ifelse(id, paste("Node", node$id, "(n = "), "n = "), node$info$nobs, ifelse(id, ")", ""), sep = "")
#        
#        grid.text(mainlab)
#        popViewport()
#        plot <- viewport(layout.pos.col = 2, layout.pos.row = 2, 
#                         xscale = xscale, yscale = yscale, 
#                         name = paste("node_density",  node$id, "plot", sep = ""))
#        pushViewport(plot)
#        yd <- ydens
#        xd <- yrange
#        if (horizontal) {
#          yyy <- xd
#          xd <- yd
#          yd <- yyy
#          yyy <- xr
#          xr <- yr
#          yr <- yyy
#          rxd <- rep(0, length(xd))
#          ryd <- rev(yd)
#        } else {
#          rxd <- rev(xd)
#          ryd <- rep(0, length(yd))
#        }
#        
#        if (rug) {
#          nodeobs <- node$info$object$y
#          if (horizontal) {
#            grid.rect(x = xscale[1], y = nodeobs , height = 0, width = xscale[1], 
#                      default.units = "native", just = c("right", "bottom"),
#                      gp = gpar(lwd = 2, col = gray(0, alpha = 0.18)))
#          } else {
#            grid.rect(x = nodeobs, y = yscale[1], 
#                      width = 0, height = abs(yscale[1]), default.units = "native", 
#                      just = c("center", "bottom"),
#                      gp = gpar(lwd = 2, col = gray(0, alpha = 0.18)))
#          }
#        }
#        
#        
#        grid.polygon(x = c(xd, rxd), y = c(yd, ryd), default.units = "native",
#                     gp = gpar(col = "black", fill = fill, lwd = lwd))
#        grid.xaxis()
#        grid.yaxis()
#        grid.rect(gp = gpar(fill = "transparent"))
#        upViewport(2)
#      }
#      return(rval)
#    }
#    class(node_density) <- "grapcon_generator"
#    
#    plot.modelparty(x, tnex = 1.7, drop = TRUE,
#         terminal_panel = node_density)
#  }
#  
#  if(type == "coef") plot.modelparty(x)
#  
#  #if(type == "hist"){
#  #  terminal_panel = node_histogram
#  #  if(is.null(tnex)) tnex <- if(is.null(terminal_panel)) 1L else 2L
#  #  if(is.null(drop_terminal)) drop_terminal <- !is.null(terminal_panel)
#  #  partykit::plot.modelparty(x, terminal_panel = terminal_panel,
#  #                            tp_args = tp_args, tnex = tnex, drop_terminal = drop_terminal, ...)
#  #}
#   
#}


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
