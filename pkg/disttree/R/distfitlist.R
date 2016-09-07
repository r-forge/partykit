distfitlist <- function(y, family, weights = NULL, start = NULL, vcov = TRUE, type.hessian = "analytic", estfun = TRUE, bd = 1, fixed = NULL, fixed.values = NULL, ...)
{
  ## match call
  cl <- match.call()

  
  ## check whether the input is a gamlss.family object (or function) or a list of the required type
  if(is.function(family)) family <- family()
  if(inherits(family, "gamlss.family")) family <- make_dist_list(family)
  # if(is.character(family)) family <- ...
  if(!is.list(family)) stop ("unknown family specification")
  if(!(all(c("ddist", "sdist", "link", "linkfun", "linkinv", "mle", "startfun") %in% names(family)))) stop("family needs to specify a list with ...")

  
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
  # eta ... coefficients of the linear predictor, here: intercept (g1(mu)=eta[1], g2(sigma)=eta[2], g3(nu)=eta[3], g4(tau)=eta[4])
  
  # if(np > 0L) m <- family$mu.linkinv(eta[1L])          # m ... mu           eta[1] ... g1(mu)        g1 ... link function
  # if(np > 1L) s <- family$sigma.linkinv(eta[2L])       # s ... sigma        eta[2] ... g2(sigma)     g2 ... link function
  # if(np > 2L) v <- family$nu.linkinv(eta[3L])          # v ... nu           eta[3] ... g3(nu)        g3 ... link function
  # if(np > 3L) t <- family$tau.linkinv(eta[4L])         # t ... tau          eta[4] ... g4(tau)       g4 ... link function
  

  
  ## set up negative log-likelihood
  nll <- function(eta) {
    nloglik <- - family$ddist(y, eta, log = TRUE, weights = weights, sum = TRUE)
    return(nloglik)
  }
  
  
  ## set up gradient
  grad <- function(eta) {
    gr <- - family$sdist(y, eta, weights = weights, sum = TRUE)
    return(gr)
  }
  
  
  ## calculate initial values if necessary or otherwise transform initial values for the distribution parameters to initial values for the intercepts
  if(is.null(start)){
    starteta <- family$startfun(y, weights = weights)
    # starteta <- family$startfun(y = rep(y, round(weights)))
  } else {
    starteta <- family$linkfun(start)
  }

  
  if(!family$mle) {
    ## optimize negative log-likelihood
    opt <- optim(par = starteta, fn = nll, gr = grad, method = "BFGS",
                 hessian = (type.hessian == "numeric"), control = list(...))
    
    ## extract parameters
    eta <- opt$par
    names(eta) <- names(starteta)
    par <- family$linkinv(eta)
    
    ## loglikelihood value
    loglik = -opt$value
    
    
  } else {
    eta <- family$startfun(y, weights)
    par <- family$linkinv(eta)
    loglik <- family$ddist(y, eta, log = TRUE, weights = weights, sum = TRUE)
    
    # use optim if numerically calculated hessian is required
    if(type.hessian == "numeric") {
      opt <- optim(par = eta, fn = nll, gr = grad, method = "BFGS",
                   hessian = (type.hessian == "numeric"), control = list(...))
    } else {opt <- NULL}
  }

  
  ## hess matrix and variance-covariance matrix estimate
  if(vcov) {
    if(type.hessian == "numeric") {
      hess <- -opt$hessian
      hess <- as.matrix(hess)
      colnames(hess) <- rownames(hess) <- names(eta)
    } else {
      hess <- family$hdist(y, eta, weights = weights)
      hess <- as.matrix(hess)
    }      
    
    
    #vcov for link coefficients eta
    vc <- solve(-hess)
    vc <- as.matrix(vc)
    colnames(vc) <- rownames(vc) <- colnames(hess)
    
  } else {
    vc <- NULL
  }
  
  
  ## estfun
  # each column represents one distribution parameter (1.col -> dldm * dmdpar = "dldeta.mu", 2.col -> dldd * dddpar = "dldeta.sigma", ...)
  if(estfun) {
    # estfun for link coefficients eta
    ef <- weights * family$sdist(y, eta, sum = FALSE)   ## FIX ME: cut out rows with weight = 0?
  } else {
    ef <- NULL                    
  }
  
  
  ##### FIX ME: additional functions with set parameters
  
  ## density function
  ddist.par <- function(x, log = FALSE) family$ddist(x, eta = eta, log = log)
  
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
    npar = length(par),
    df = length(par),
    y = y,
    weights = weights,
    family = family$family.name,
    starteta = starteta,
    opt = opt,
    par = par,
    eta = eta,
    hess = hess,
    vcov = vc,
    loglik = loglik,
    call = cl,
    ny = ny,        
    nobs = nobs,
    estfun = ef,
    familylist = family
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

coef.distfit <- function(object, type = "link" , ...) {
  if(type == "link") return(object$eta)
  if(type == "parameter") return(object$par)
  ## FIXME: else, warning
}

vcov.distfit <- function(object, type = "link", ...) {
  if(type == "link") return(object$vcov)
  if(type == "parameter"){
    ## delta method
    delta.m <- diag(object$familylist$linkinvdr(object$eta))
    colnames(delta.m) <- rownames(delta.m) <- names(object$par)
    return(delta.m %*% object$vcov %*% delta.m)
  }
}
  
estfun.distfit <- function(object, ...) {                         
  return(object$estfun)
}

logLik.distfit <- function(object, ...) {
  structure(object$loglik, df = object$df, class = "logLik")
}

bread.distfit <- function(object, ...) {
  return(object$vcov * object$nobs)
}

confint.distfit <- function(object, parm, type = "link", level = 0.95, ...) {
  np <- object$npar
  
  if(type == "link"){ 
    coef <- object$eta
    vcov <- object$vcov
  }
  if(type == "parameter"){ 
    coef <- object$par
    vcov <- vcov(object, type = "parameter")
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
    vcov <- object$vcov
  }
  if(type == "parameter"){ 
    coef <- object$par
    vcov <- vcov(object, type = "parameter")
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
