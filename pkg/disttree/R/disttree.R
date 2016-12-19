## high-level convenience interface to mob() and ctree()
# FIX ME: default settings for family, decorrelate only necesary for type.tree == "ctree"
disttree <- function(formula, data, na.action, cluster, family = NO(),
                     type.tree = "mob", decorrelate = "none",
                     control = mob_control(...), ocontrol = list(), ...)
{
  ## keep call
  cl <- match.call(expand.dots = TRUE)
  
  # check input arguments
  if(!(type.tree %in% c("mob", "ctree"))) stop("unknown argument for type.tree (can only be mob or ctree)")
  if(!(decorrelate) %in% c("none", "opg", "vcov")) stop("unknown argument for decorrelate (can only be none, opg or vcov)")
  
  if(type.tree == "mob") {
    
    ## glue code for calling distfit() with given family in mob()
    dist_family_fit <- function(y, x = NULL, start = NULL, weights = NULL, offset = NULL,
                                cluster = NULL, vcov = FALSE, estfun = TRUE, 
                                object = FALSE, type.hessian = "analytic", ...)
    {
      if(!(is.null(x) || NCOL(x) == 0L)) warning("x not used")
      if(!is.null(offset)) warning("offset not used")
      rval <- distfit(y, family = family, weights = weights, start = start,
                      vcov = vcov, estfun = estfun, type.hessian = type.hessian, ...)
      rval <- list(
        coefficients = rval$par,
        objfun = -rval$loglik,
        estfun = if(estfun) rval$estfun else NULL,   ## rval$estfun contains the scores of the positive loglik 
        object = if(object) rval else NULL
      )
      return(rval)
    }
    
    ## call mob
    m <- match.call(expand.dots = FALSE)
    m$fit <- dist_family_fit
    m$family <- NULL
    # m$family <- m$ocontrol <- NULL
    for(n in names(ocontrol)) m[[n]] <- ocontrol[[n]]
    if("..." %in% names(m)) m[["..."]] <- NULL
    if("type.tree" %in% names(m)) m[["type.tree"]] <- NULL
    m[[1L]] <- as.name("mob")
    rval <- eval(m, parent.frame())
  }
  
  
  if(type.tree == "ctree") {
    
    ## wrapper function to apply distfit in ctree
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
    
    
    ## call ctree
    m <- match.call(expand.dots = FALSE)
    m$ytrafo <- modelscores_decor
    # m$ocontrol <- NULL
    # m$family <- m$ocontrol <- NULL
    for(n in names(ocontrol)) m[[n]] <- ocontrol[[n]]
    if("..." %in% names(m)) m[["..."]] <- NULL
    if("type.tree" %in% names(m)) m[["type.tree"]] <- NULL
    m[[1L]] <- as.name("ctree")
    rval <- eval(m, parent.frame())
  
  }
  
  
  
  ## extend class and keep original call/family/control
  rval$info$call <- cl
  rval$info$family <- family
  rval$info$ocontrol <- ocontrol
  class(rval) <- c("disttree", class(rval))
  return(rval)
}


## methods
print.disttree <- function(x,
  title = NULL, objfun = "negative log-likelihood", ...)
{
  if(is.null(title)) title <- sprintf("Distributional regression tree (%s)", x$info$family$family[2L])
  partykit::print.modelparty(x, title = title, objfun = objfun, ...)
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
  
  test <- disttree(y~x, data = d, family = dist_list_normal)
  print(test)
  plot(test)
  plot(as.constparty(test))
}
