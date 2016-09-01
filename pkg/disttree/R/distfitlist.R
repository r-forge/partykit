distfitlist <- function(y, family, weights = NULL, start = NULL, vcov = TRUE, type.hessian = "analytic", estfun = TRUE, bd = 1, fixed = NULL, fixed.values = NULL, ...)
{
  ## match call
  cl <- match.call()
  

  
  ## check whether the input is a gamlss.family object (or function) or a list of the required type
  if(is.function(family)) family <- family()
  if(inherits(family, "gamlss.family")) family <- make_dist_list(family)
  # if(is.character(family)) family <- ...
  if(!is.list(family)) stop ("unknown family specification")
  if(!(all(c("ddist", "sdist", "hdist", "link", "use.optim") %in% names(family)))) stop("family needs to specify a list with ...")

  
  
  # FIXME:
  # the input argument fixed can be a character string or a list of character strings with the name(s) of the parameter(s) which are fixed
  # in this case their values must be set in the argument fixed.values
  
  # 2 families include fixed parameter(s): LNO() ... nu is fixed
  #                                        NET() ... nu and tau are fixed
  # any(family$parameters == FALSE)
  
  
  ## number of observations
  ny <- NROW(y)
  
  ## weights
  if(is.null(weights)) weights <- 1
  if(length(weights) == 1L) weights <- rep.int(weights, ny)
  weights <- as.vector(weights)
  
  ## store y and select observations with weight > 0 
  #FIXME# y.store <- y          
  #FIXME# y <- y[weights > 0]
  
  ## number of observations = sum of weights (i.e., case weights)
  ## FIXME ## also need proportionality weights, i.e., weights = sum(weights > 0) ?
  nobs <- sum(weights)
  
  ## notation:
  # par ... distribution parameters (mu, sigma, nu, tau)
  # eta ... coefficients of the linear predictor, here: intercept (g(mu)=eta[1], g(sigma)=eta[2], g(nu)=eta[3], g(tau)=eta[4])
  
  # if(np > 0L) m <- family$mu.linkinv(eta[1L])          # m ... mu           eta[1] ... g(mu)        g ... link function
  # if(np > 1L) s <- family$sigma.linkinv(eta[2L])       # s ... sigma        eta[2] ... g(sigma)     g ... link function
  # if(np > 2L) v <- family$nu.linkinv(eta[3L])          # v ... nu           eta[3] ... g(nu)        g ... link function
  # if(np > 3L) t <- family$tau.linkinv(eta[4L])         # t ... tau          eta[4] ... g(tau)       g ... link function
  

  
  ## set up negative log-likelihood
  nll <- function(eta) {
    nloglik <- -family$ddist(y, eta, log = TRUE, type = "link")
    nloglik <- sum(weights * nloglik)
    return(nloglik)
  }
  
  
  ## set up gradient
  grad <- function(eta, sum = TRUE) {
    gr <- -weights * family$sdist(y, eta, type = "link")
    if(sum) gr <- colSums(gr)     # *1/nobs ? scale, doesn't influence optimization
    return(gr)
  }
  
  
  ## calculate initial values if necessary or otherwise transform initial values for the distribution parameters to initial values for the intercepts
  if(is.null(start)){
    starteta <- family$start.eta(y = rep(y, round(weights)))
    startpar <- family$link.inv(starteta)
    # startpar <- family$start(y)
  } else {
    startpar <- start
    starteta <- family$link.fun(startpar)
    names(startpar) <- names(family$link.inv(starteta))   ## FIX ME
  }

  
  if(family$use.optim) {
    ## optimize log-likelihood
    opt <- optim(par = starteta, fn = nll, gr = grad, method = "BFGS",
                 hessian = (type.hessian == "numeric"), control = list(...))
    
    ## extract parameters
    eta <- opt$par
    par <- family$link.inv(eta)
    
    ## loglikelihood value
    loglik = -opt$value
    
    
  } else {
    par <- family$closed.mle(y, weights)
    eta <- family$link.fun(par)
    loglik <- sum(weights * family$ddist(y, par, log = TRUE, type = "parameter"))
    
    # use optim if numerically calculated hessian is required
    if(type.hessian == "numeric") {
      starteta <- eta
      opt <- optim(par = starteta, fn = nll, gr = grad, method = "BFGS",
                   hessian = (type.hessian == "numeric"), control = list(...))
    } else {opt <- NULL}
  }
  
  
  names(eta) <- names(starteta)
  names(par) <- names(startpar)
  
  
  ## hess matrix for distribution parameter  (FIX ME: until now only analytic, even if type.hessian = "numeric")
  if(type.hessian == "numeric") {
    hess.eta <- -opt$hessian
    hess.eta <- as.matrix(hess.eta)
    hess.par <- NULL    ## FIX
    hess.par <- as.matrix(hess.par)
  } else {
    hess.eta <- family$hdist(y, eta, type = "link", weights = weights)
    hess.eta <- as.matrix(hess.eta)
    hess.par <- family$hdist(y, par, type = "parameter", weights = weights)
    hess.par <- as.matrix(hess.par)
  }
      
  
  
  ## variance-covariance matrix estimate 
  if(vcov){
    
    # vcov for distribution parameter
    vc.par <- solve(-hess.par)
    vc.par <- as.matrix(vc.par)
    colnames(vc.par) <- rownames(vc.par) <- colnames(hess.par)
    
    #vcov for link coefficients eta
    vc.eta <- solve(-hess.eta)
    vc.eta <- as.matrix(vc.eta)
    colnames(vc.eta) <- rownames(vc.eta) <- colnames(hess.eta)
    
  } else {
    vc.par <- NULL
    vc.eta <- NULL
  }
  
  
  ## estfun
  # each column represents one distribution parameter (1.col -> dldm * dmdpar = "dldeta.mu", 2.col -> dldd * dddpar = "dldeta.sigma", ...)
  if(estfun) {
    
    # estfun for distribution parameter
    ef.par <- weights * family$sdist(y, par, type = "parameter")  
    
    # estfun for link coefficients eta
    ef.eta <- weights * family$sdist(y, eta, type = "link")
    
  } else {
    ef.par <- NULL
    ef.eta <- NULL                    
  }
  
  
  ##### FIX ME: additional functions with set parameters
  
  ## density function
  ddist.par <- function(x, log = FALSE) family$ddist(x, par = par, log = log, type = "parameter")
  
  ## density function
  #ddist <- get(paste0("d",family$family[1]))
  
  ## cumulative distribution function
  #pdist <- get(paste0("p",family$family[1]))
  
  ## quantile function
  #qdist <- get(paste0("q",family$family[1]))
  
  ## random function
  #rdist <- get(paste0("r",family$family[1]))
  

  
  
  ## return value 
  # FIX: return bd? df = np?
  rval <- list(
    np = length(par),
    df = length(par),
    y = y,
    weights = weights,
    family = family$family.name,
    startpar = startpar,
    starteta = starteta,
    opt = opt,
    par = par,
    eta = eta,
    hess = hess.par,
    hess.eta = hess.eta,
    vcov = vc.par,
    vcov.eta = vc.eta,
    loglik = loglik,
    call = cl,
    ny = ny,        
    nobs = nobs,
    estfun = ef.par,
    estfun.eta = ef.eta
    #ddist = ddist,
    #pdist = pdist,
    #qdist = qdist,
    #rdist = rdist
  )
  
  #if(any(family$family%in%.distfit.bi.list)) rval <- c(rval, bd = bd)
  
  class(rval) <- "distfit"
  return(rval)
}









## print, summary?, predict?
nobs.distfit <- function(object, ...) {
  object$nobs
}

coef.distfit <- function(object, type = "parameter" , ...) {
  if(type == "link") return(object$eta)
  if(type == "parameter") return(object$par)
  ## FIXME: else, warning
}

vcov.distfit <- function(object, type = "parameter", ...) {
  if(type == "link") return(object$vcov.eta)
  if(type == "parameter") return(object$vcov)
  ## delta method
  #delta.m <- diag(object$npar)
  #delta.m[1,1] <- object$family$mu.dr(object$eta[1])
  #if(object$npar > 1L) delta.m[2,2] <- object$family$sigma.dr(object$eta[2])
  #if(object$npar > 2L) delta.m[3,3] <- object$family$nu.dr(object$eta[3])
  #if(object$npar > 3L) delta.m[4,4] <- object$family$tau.dr(object$eta[4])
  #colnames(delta.m) <- rownames(delta.m) <- names(object$par)
  #return(delta.m %*% object$vcov %*% delta.m)
}

estfun.distfit <- function(object, type = "parameter", ...) {                         
  if(type == "link") return(object$estfun.eta)
  if(type == "parameter") return(object$estfun)
}

logLik.distfit <- function(object, ...) {
  structure(object$loglik, df = object$df, class = "logLik")
}

bread.distfit <- function(object, type = "parameter", ...) {
  if(type == "link") return(object$vcov.eta * object$nobs)
  if(type == "parameter") return(object$vcov * object$nobs)
}

confint.distfit <- function(object, parm, level = 0.95, type = "parameter", ...) {
  np <- object$npar
  if(type == "link"){ 
    vcov <- object$vcov.eta
    coef <- object$eta
  }
  if(type == "parameter"){ 
    vcov <- object$vcov
    coef <- object$par
  }
  
  left <- (1-level)/2
  right <- 1-left
  
  if(missing(parm)){
    use.parm <- rep(TRUE,length = np)
  } else {
    use.parm <- logical(length = np)
    if(("mu" %in% parm)    || (paste0(object$family$mu.link,"(mu)") %in% parm)       || 1 %in% parm) use.parm[1] <- TRUE
    if(("sigma" %in% parm) || (paste0(object$family$sigma.link,"(sigma)") %in% parm) || 2 %in% parm) use.parm[2] <- TRUE
    if(("nu" %in% parm)    || (paste0(object$family$nu.link,"(nu)") %in% parm)       || 3 %in% parm) use.parm[3] <- TRUE
    if(("tau" %in% parm)   || (paste0(object$family$tau.link,"(tau)") %in% parm)     || 4 %in% parm) use.parm[4] <- TRUE
  }
  
  confint <- NULL
  if((np > 0L) && use.parm[1]){
    confint1 <- c(coef[1] + qnorm(left) * sqrt(vcov[1,1]), coef[1] + qnorm(right) * sqrt(vcov[1,1]))
    confint <- rbind(confint, confint1)
  } 
  if((np > 1L) && use.parm[2]){
    confint2 <- c(coef[2] + qnorm(left) * sqrt(vcov[2,2]), coef[2] + qnorm(right) * sqrt(vcov[2,2]))
    confint <- rbind(confint, confint2)
  }
  if((np > 2L) && use.parm[3]){
    confint3 <- c(coef[3] + qnorm(left) * sqrt(vcov[3,3]), coef[3] + qnorm(right) * sqrt(vcov[3,3]))
    confint <- rbind(confint, confint3)
  }
  if((np > 3L) && use.parm[4]){ 
    confint4 <- c(coef[4] + qnorm(left) * sqrt(vcov[4,4]), coef[4] + qnorm(right) * sqrt(vcov[4,4]))
    confint <- rbind(confint, confint4)
  }
  
  confint <- as.matrix(confint)
  colnames(confint) <- c(paste0(left," %"), paste0(right," %"))
  rownames(confint) <- names(coef[use.parm])
  
  confint
}


## FIXME: summary (further information)
summary.distfit <- function (object, type = "parameter", ...){
  if(type == "link"){ 
    coef <- object$eta
    vcov <- object$vcov.eta
  }
  if(type == "parameter"){ 
    coef <- object$par
    vcov <- object$vcov
  }
  
  se <- sqrt(diag(vcov))
  TAB <- cbind(Estimate = coef,
               StdErr = se)
  
  
  
  sumlist <- list(Call = object$call,
                  Family = object$family,
                  Parameter = TAB,
                  Estimated_covariance_matrix = vcov,
                  LogLikelihood = object$loglik
                  # LogLikelihood = logLik(object)
  )
  class(sumlist) <- "summary.distfit"
  sumlist 
}
