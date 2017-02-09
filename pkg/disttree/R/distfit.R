distfit <- function(y, family, weights = NULL, start = NULL, vcov = TRUE, type.hessian = "analytic", estfun = TRUE, bd = NULL, fixed = NULL, fixed.values = NULL, ...)
{
  ## match call
  cl <- match.call()

  
  ## check whether the input is a gamlss.family object (or function) or a list of the required type
  if(is.function(family)) family <- family()
  if(inherits(family, "gamlss.family")) family <- make_dist_list(family, bd = bd)
  # if(is.character(family)) family <- ...     # for biniomial distributions: bd should be handed over once, but not appear in the list from here on
  if(!is.list(family)) stop ("unknown family specification")
  if(!(all(c("ddist", "sdist", "link", "linkfun", "linkinv", "mle", "startfun") %in% names(family)))) stop("family needs to specify a list with ...")
  # linkinvdr only used in the method vcov for type = "parameter"
  
  if(type.hessian == "numeric") family$hdist <- NULL
  if(type.hessian == "analytic" && is.null(family$hdist)) stop("analytic calculation of hessian matrix not possible without list element hdist ...")

  
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
  # nobs <- sum(weights)
  
  ## notation:
  # par ... distribution parameters (mu, sigma, nu, tau)
  # eta ... coefficients of the linear predictor, here: intercept (g1(mu)=eta[1], g2(sigma)=eta[2], g3(nu)=eta[3], g4(tau)=eta[4])
  
  # if(np > 0L) m <- family$mu.linkinv(eta[1L])          # m ... mu           eta[1] ... g1(mu)        g1 ... link function
  # if(np > 1L) s <- family$sigma.linkinv(eta[2L])       # s ... sigma        eta[2] ... g2(sigma)     g2 ... link function
  # if(np > 2L) v <- family$nu.linkinv(eta[3L])          # v ... nu           eta[3] ... g3(nu)        g3 ... link function
  # if(np > 3L) t <- family$tau.linkinv(eta[4L])         # t ... tau          eta[4] ... g4(tau)       g4 ... link function
  

  
  ## set up negative log-likelihood
  nll <- function(eta) {
    
    #rval <- family$ddist(y, eta, log = TRUE, weights = weights, sum = FALSE)
    
    rval <- suppressWarnings(try(family$ddist(y, eta, log = TRUE, weights = weights, sum = FALSE), silent = TRUE))
    if(inherits(rval, "try-error")) {
      #print("try-error in ddist")
      return(1.7e+308)
    } else {
      if(any(is.na(rval))) {
        #print("NAs in ddist")
        return(1.7e+308)
      } else {
        if(any(abs(rval) == Inf)) {
          #print(c("Infinity in ddist"))
          return(1.7e+308)
        } else {
          nloglik <- - sum(weights * rval)
          if(abs(nloglik) == Inf) return(1.7e+308)
          return(nloglik)
        }
      }
    }
  }  
  
  
  ## set up gradient
  grad <- function(eta) {
    gr <- - family$sdist(y, eta, weights = weights, sum = TRUE)
    if(any(gr == -Inf)) {
      gr[gr==-Inf] = -1
      #print(c(eta,"gr = -Inf"))
    }
    if(any(gr ==  Inf)) {
      gr[gr== Inf] =  1
      #print(c(eta,"gr = Inf"))
    }
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
    opt <- try(optim(par = starteta, fn = nll, gr = grad, method = "L-BFGS-B",
                 hessian = (type.hessian == "numeric"), control = list(...)))
    if(inherits(opt, "try-error")) {
      #print("opt try-error with gr=grad, L-BFGS-B")
      opt <- try(optim(par = starteta, fn = nll, gr = grad, method = "BFGS",
                       hessian = (type.hessian == "numeric"), control = list(...)))
      if(inherits(opt, "try-error")) {
        #print("opt try-error with gr = grad, BFGS")
        opt <- optim(par = starteta, fn = nll, method = "BFGS",
                         hessian = (type.hessian == "numeric"), control = list(...))
      }
    }
    
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
    opt <- NULL
  }

  
  ## hess matrix (second-order partial derivatives of the (positive) log-likelihood function) and Hessian estimate for the variance-covariance matrix
  if(vcov) {
    if(is.null(family$hdist)) {
      if(family$mle) {
        nhess <- try(optim(par = eta, fn = nll, gr = grad, method = "L-BFGS-B",
                         hessian = TRUE, control = list(...))$hessian)
        if(inherits(nhess, "try-error")) {
          #print("opt try-error in hess with gr=grad, L-BFGS-B")
          nhess <- try(optim(par = eta, fn = nll, gr = grad, method = "BFGS",
                       hessian = TRUE, control = list(...))$hessian)
          if(inherits(nhess, "try-error")) {
            #print("opt try-error in hess with gr=grad, BFGS")
            nhess <- optim(par = eta, fn = nll, method = "BFGS",
                               hessian = TRUE, control = list(...))$hessian
          }
        }
        hess <- -nhess
      } else {
        hess <- -opt$hessian
      }
    } else {
      hess <- family$hdist(y, eta, weights = weights)
    }
    hess <- as.matrix(hess)
    colnames(hess) <- rownames(hess) <- names(eta)
    
    
    #vcov for link coefficients eta
    vc <- solve(-hess)
    vc <- as.matrix(vc)
    colnames(vc) <- rownames(vc) <- colnames(hess)
    
  } else {
    hess <- NULL
    vc <- NULL
  }
  
  
  ## estfun
  # each column represents one distribution parameter (1.col -> dldm * dmdpar = "dldeta.mu", 2.col -> dldd * dddpar = "dldeta.sigma", ...)
  # (first-order partial derivatives of the (positive) log-likelihood function)
  if(estfun) {
    # estfun for link coefficients eta
    ef <- weights * family$sdist(y, eta, sum = FALSE)   ## FIX ME: cut out rows with weight = 0? -> No! index is of importance for independence tests (relation to covariates)
  } else {
    ef <- NULL                    
  }
  
  
  ##### FIX ME: additional functions with set parameters
  
  ## density function
  ddist <- function(x, log = FALSE) family$ddist(x, eta = eta, log = log)
  
  ## density function
  # ddist <- family$ddist
  
  ## cumulative distribution function
  # pdist <- family$pdist
  
  ## quantile function
  # qdist <- family$qdist
  
  ## random function
  # rdist <- family$rdist
  

  
  
  ## return value 
  rval <- list(
    npar = length(par),
    y = y,
    ny = ny,
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
    estfun = ef,
    familylist = family,
    ddist = ddist
    #pdist = pdist,
    #qdist = qdist,
    #rdist = rdist
  )
  
  class(rval) <- "distfit"
  return(rval)
}





## print, summary?, predict?
nobs.distfit <- function(object, ...) {
  object$ny
}

coef.distfit <- function(object, type = "link" , ...) {
  if(type == "link") return(object$eta)
  if(type == "parameter") return(object$par)
  ## FIXME: else, warning
}

# FIX: complete with other types ?
predict.distfit <- function(object, type = "response", OOB = FALSE, ...){
  # calculation of the expected value 
  # of the given distribution with the calculated parameters
  if(type == "response"){
    f <- function(x){x * object$ddist(x, log = FALSE)}
    expv <- integrate(f,-Inf,Inf )
    return(expv[[1]])
  }
  
  ## FIX: if censored distribution -> change integration boundaries
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
  structure(object$loglik, df = object$npar, class = "logLik")
}

bread.distfit <- function(object, ...) {
  return(object$vcov * object$ny)
}

confint.distfit <- function(object, parm, level = 0.95, type = "link", ...) {
  np <- object$npar
  
  if(type == "link"){ 
    coef <- object$eta
    vcov <- object$vcov
    # FIX ME: vcov on link scale: values around zero might be negative => error using sqrt
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
    if(("mu" %in% parm) || (paste0(object$familylist$link[1],"(mu)") %in% parm) || 1 %in% parm) use.parm[1] <- TRUE
    if(np > 1L) {
      if(("sigma" %in% parm) || (paste0(object$familylist$link[2],"(sigma)") %in% parm) || 2 %in% parm) use.parm[2] <- TRUE
      if(np > 2L) {
        if(("nu" %in% parm) || (paste0(object$familylist$link[3],"(nu)") %in% parm) || 3 %in% parm) use.parm[3] <- TRUE
        if(np > 3L) if(("tau" %in% parm) || (paste0(object$familylist$link[4],"(tau)") %in% parm) || 4 %in% parm) use.parm[4] <- TRUE
      }
    }
  }
  
  confint <- NULL
  if(use.parm[1]) {
    confint1 <- c(coef[1] + qnorm(left) * sqrt(vcov[1,1]), coef[1] + qnorm(right) * sqrt(vcov[1,1]))
    confint <- rbind(confint, confint1)
  } 
  if(np > 1L) {
    if(use.parm[2]) {
      confint2 <- c(coef[2] + qnorm(left) * sqrt(vcov[2,2]), coef[2] + qnorm(right) * sqrt(vcov[2,2]))
      confint <- rbind(confint, confint2)
    }
    if(np > 2L) {
      if(use.parm[3]) {
        confint3 <- c(coef[3] + qnorm(left) * sqrt(vcov[3,3]), coef[3] + qnorm(right) * sqrt(vcov[3,3]))
        confint <- rbind(confint, confint3)
      }
      if(np > 3L) {
        if(use.parm[4]) { 
          confint4 <- c(coef[4] + qnorm(left) * sqrt(vcov[4,4]), coef[4] + qnorm(right) * sqrt(vcov[4,4]))
          confint <- rbind(confint, confint4)
        }
      }
    }
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


plot.distfit <- function(object,
                         main = "", xlab = "", ylab = "Density",
                         fill = "lightgray", col = "darkred", lwd = 1.5,
                         ...)
{
  ## FIX ME: barplot instead of hist for discrete distributions
  if(isTRUE(all.equal(object$y, round(object$y)))) {
    
    ## barplot instead of hist:
    #ytab <- table(object$y)
    #values <- as.numeric(names(ytab))
    #absfreq <- as.numeric(ytab)
    #relfreq <- absfreq / sum(absfreq)
    #relytab <- ytab / sum(absfreq)
    #bars <- barplot(relytab, xlim = c(min(values), max(values)), ylim = c(-0.1, 1.1),
    #                xlab = xlab , ylab = ylab)
    #abline(h = 0)
    #lines(values*1.2 + 0.7, object$ddist(values), type = "b", lwd = 2, pch = 19, col = 2)
    
    # using hist:
    histogram <- c(
      list(x = object$y, main = main, xlab = xlab, ylab = ylab, col = fill),
      list(...)
    )
    histogram$freq <- FALSE
    histogram$probability <- TRUE
    histogram$breaks <- seq(from = min(object$y), to = max(object$y) + 1) - 0.5
    # histogram$breaks <- seq(from = min(object$y), to = 2*max(object$y) + 1)/2 - 0.25
    histogram <- do.call("hist", histogram)
    yrange <- seq(from = min(object$y), to = max(object$y))
    lines(yrange, object$ddist(yrange), col = col, lwd = lwd)
    
  } else {
    
    histogram <- c(
      list(x = object$y, main = main, xlab = xlab, ylab = ylab, col = fill),
      list(...)
    )
    histogram$freq <- FALSE
    histogram$probability <- TRUE
  
    histogram <- do.call("hist", histogram)
    yrange <- seq(from = histogram$breaks[1L],
                  to = histogram$breaks[length(histogram$breaks)],
                  length.out = 100L)
    lines(yrange, object$ddist(yrange), col = col, lwd = lwd)
    ## FIXME: ddist arguments/par?
  }
}




if(FALSE) {
  
  family <- ZABI()
  y <- rZABI(1000, bd = 10, mu = 0.5, sigma = 0.2)
  ny <- length(y)
  start <- c(0.8, 0.1)
  weights <- rbinom(ny, 1, 0.75)
  
  df <- distfit(y, family, weights = weights, start = start, bd = 10)
  df2 <- distfit(y, family, bd = 10)
  
}


if(FALSE) {
  ####### example normal distribution
  family <- NO()
  y <- rNO(1000, mu = 5, sigma = 2)
  ny <- length(y)
  start <- c(2, 1)
  weights <- rbinom(ny, 1, 0.75)
  
  df <- distfit(y, family, start = start)
  df2 <- distfit(y, family = dist_list_normal, start = start)
  df3 <- distfit(y, family, start = start, type.hessian = "numeric")
  df4 <- distfit(y, family = dist_list_normal, start = start, type.hessian = "numeric")
  
  coef(df)
  coef(df2)
  coef(df3)
  coef(df4)
  vcov(df)
  vcov(df2)
  vcov(df3)
  vcov(df4)
  solve(vcov(df))
  solve(vcov(df2))
  solve(vcov(df3))
  solve(vcov(df4))
}





if(FALSE){
  ######### examples
  ## simulate artifical negative binomial data
  set.seed(0)
  y <- rnbinom(1000, size = 1, mu = 2)
  
  ## simple distfit
  df <- distfit(y, family = NBI)
  coef(df)
  confint(df)
  logLik(df)
  
  ## using tabulated data
  ytab <- table(y)
  df2 <- distfit(as.numeric(names(ytab)), family = NBI, weights = ytab)
  coef(df2)
  confint(df2)
  logLik(df2)
  
  ## coefficients tests
  if(require("lmtest")) {
    coeftest(df)
    coeftest(df2)
  }
  
  ## censored logistic example
  if(require("crch") & require("gamlss.cens")) {
    library("crch")
    data("RainIbk", package = "crch")
    m1 <- crch(rain ~ 1, data = RainIbk, left = 0, dist = "logistic")
    
    library("gamlss.cens")
    gen.cens(LO, type = "left")
    RainIbk$rains <- Surv(RainIbk$rain, RainIbk$rain > 0, type = "left")
    m2 <- distfit(RainIbk$rains, family = LOlc)
    # FIX ME: calculation of starting values for censored distributions
    # m2 <- distfit(RainIbk$rains, family = LOlc, start = c(1,1))
    
    coef(m1)
    coef(m2)
    logLik(m1)
    logLik(m2)
    vcov(m1)
    vcov(m2)
  }
  
  
  
  ## censored normal example
  y <- rnorm(500,3,2)
  y[y<0] <- 0
  m1 <- crch(y ~ 1, left = 0, dist = "gaussian")
  m2 <- distfit(Surv(y, y>0, type="left"), family = NOlc)
  m3 <- distfit(y, family = dist_list_cens_normal)
  
  coef(m1)
  coef(m2)
  coef(m3)
  logLik(m1)
  logLik(m2)
  logLik(m3)
  vcov(m1)
  vcov(m2)
  vcov(m3)
  solve(vcov(m1))
  solve(vcov(m2))
  solve(vcov(m3))
  
}



if(FALSE){
  library("disttree")
  family <- dist_list_weibull
  family$family.name
  y <- survival:::rsurvreg(1000, mean = 5, scale = 2) 
  # y <- rweibull(1000, shape = 1/2, scale = exp(5))
  weights <- rbinom(length(y), 4, 0.8)
  t <- distfit(y, family = family, weights = weights, vcov = TRUE, type.hessian = "numeric")
  t$par
  start <- c(3, 1)
  ts <- distfit(y, family = family, weights = weights, start = start, vcov = FALSE)
  ts$par

  # different parametrizations in WEI, WEI2 and WEI3 (calculated based on mu, sigma)
}




if(FALSE){
  family <- BE()
  family$nopar
  family$family
  y <- rBE(10000, 0.5, 0.3)
  weights <- rbinom(length(y), 1, 0.8)
  t <- distfit(y, family = BE(), weights = weights)
  t$par
}

if(FALSE){
  family <- BCCG()
  family$nopar
  family$family
  y <- rBCCG(10000, 1, 0.1, 0.1)
  weights <- rbinom(length(y), 1, 0.8)
  t <- distfit(y, family = BCCG(), weights = weights)
  t$par
}

if(FALSE){
  family <- BCPE()
  family$nopar
  family$family
  y <- rBCPE(10000, 5, 0.1, 1, 2)
  weights <- rbinom(length(y), 1, 0.8)
  t <- distfit(y, family = BCPE(), weights = weights)
  t$par
}

if(FALSE){
  family <- BCT()
  family$nopar
  family$family
  y <- rBCT(10000, 5, 0.1, 1, 2)
  weights <- rbinom(length(y), 1, 0.8)
  t <- distfit(y, family = BCT(), weights = weights)
  t$par
}

if(FALSE){
  family <- EGB2()
  family$nopar
  family$family
  y <- rEGB2(10000, 0, 1, 1, 0.5)
  weights <- rbinom(length(y), 1, 0.8)
  t <- distfit(y, family = EGB2(), weights = weights)
  t$par
}

if(FALSE){
  family <- LO()
  family$nopar
  family$family
  y <- rLO(1000, 0, 1)
  weights <- rbinom(length(y), 1, 0.8)
  t <- distfit(y, family = LO(), weights = weights)
  t$par
}

if(FALSE){
  family <- JSUo()
  family$nopar
  family$family
  y <- rJSUo(1000, 0, 1, 0, 0.5)
  weights <- rbinom(length(y), 1, 0.8)
  t <- distfit(y, family = JSUo(), weights = weights)
  t$par
}


if(FALSE){
  family <- LNO()
  family$nopar
  family$family
  y <- rLNO(10000, 1, 0.8, 0.3)
  weights <- rbinom(length(y), 1, 0.8)
  t <- distfit(y, family = LNO(), weights = weights)
  t$par
}

