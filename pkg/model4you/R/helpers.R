
#' Fit function when model object is given
#'
#' Use update function to refit model and extract info such as coef, logLik and 
#' estfun.
#'
#' @param formula ignored but required by \code{.ctreetrafo}.
#' @param model model object.
#' @param data data.
#' @param weights weights.
#' @param cluster cluster.
#' @param ctrl control options from \code{ctree_control}.
#' @param parm which parameters should be used for instability test?
#'
#' @return A function returning a list of
#' \item{coefficients}{ \code{coef}. }
#' \item{objfun}{ \code{logLik}. }
#' \item{object}{ the model object. }
#' \item{converged}{ Did the model converge? }
#' \item{estfun}{ \code{estfun}. }
.modelfit <- function(formula, # ignored but needed for .ctreetrafo
                      model, data, weights, cluster, ctrl, parm = NULL) {
  
  fitfun <- function(subset, estfun = TRUE, object = FALSE, info = NULL) { 
    
    ## compute model on data
    if(length(weights) == 0) weights <- rep(1, NROW(data))
    mod <- tryCatch(update(object = model, data = data, subset = subset),
                    warning = function(w) list(converged = FALSE),
                    error = function(e) {
                      list(converged = FALSE)
                    })
    
    ## stop if error
    if(!(is.null(mod$converged)) && !(mod$converged)) 
      return(list(converged = FALSE))
    
    # mod <- try(update(object = model, data = dat),
    #            silent = TRUE)
    
    ## get convergence info
    if (is.null(ctrl$converged)) {
      conv <- if (is.null(mod$converged)) TRUE else mod$converged
    } else {
      conv <- ctrl$converged(mod, data, subset)
    }
    
    ## prepare return list
    ret <- list(coefficients = coef(mod), objfun = logLik(mod),
                object = if(object) mod else NULL,
                converged = conv)
    
    ## add estfun if wanted
    if(estfun) {
      ef <- estfun(mod)
      ret$estfun <- matrix(0, nrow = NROW(data), ncol = NCOL(ef))
      ret$estfun[subset,] <- ef
      if(!is.null(parm)) ret$estfun <- ret$estfun[, parm]
    }
    
    return(ret)
  }
  
  return(fitfun)
}

.prepare_args <- function(object, data, zformula, control, ...) {
  
  if (is.null(modcall <- getCall(object))) 
    stop("Need an object with call component, see getCall")
  
  ## get arguments for cforest call
  args <- list(...)
  # args$ntree <- ntree
  args$control <- control
  
  # ## arguments used in model
  # modargs <- as.list(modcall)[-1]
  
  ## formula and data
  if(is.null(data)) data <- eval(modcall$data)
  args$data <- data
  modformula <- eval(modcall$formula)
  
  ## in case I switch to mob
  # if(is.null(zformula)) zformula <- formula(~ .)
  # mobformula <- as.Formula(modformula, zformula)
  
  ## cforest formula
  if(is.null(zformula)) zformula <- "~ ."
  if(!is.character(zformula)) zformula <- paste(as.character(zformula), collapse = " ")
  modvars <- all.vars(modformula)
  args$formula <- as.formula(
    paste(paste(modvars, collapse = " + "), zformula)
  )
  
  return(args)
}