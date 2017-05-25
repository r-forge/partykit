distfit <- function(y, family, weights = NULL, start = NULL, start.eta = NULL, 
                    vcov = TRUE, type.hessian = "analytic", estfun = TRUE, 
                    bd = NULL, fixed = NULL, fixed.values = NULL, 
                    cens = "none", censpoint = NULL, ocontrol = list(), ...)
{
  ## start on par scale
  ## start.eta on link scale
  
  
  ## FIX ME: what to do if weights consists of zeros only
  ## FIX ME: hand over ocontrol arguments to optim in the right form/order
  
  ## match call
  cl <- match.call()

  ## number of observations
  ny <- NROW(y)
  
  ## check weights
  if(is.null(weights) || (length(weights)==0L)) weights <- as.vector(rep.int(1, ny))
  if(length(weights) != ny) stop("number of observations and length of weights are not equal")
    
  
  ## check whether the input is a gamlss.family object (or function) or a list of the required type
  if(is.function(family)) family <- family()
  
  # if family is a gamlss.dist family object and the distribution is censored:
  # data has to be converted to a Survival object (necessary input arguments: cens and censpoint)
  if(inherits(family, "gamlss.family")) {
    if(("censored" %in% strsplit(family$family[2], " ")[[1]]) && (!survival::is.Surv(y))) {
      if(cens == "none" || is.null(censpoint)) stop("for censored gamlss.dist family objects the censoring type and point(s) have to be set (input arguments cens and censpoint)")
      if(cens == "left") y <- survival::Surv(y, y > censpoint, type = "left")
      if(cens == "right") y <- survival::Surv(y, y < censpoint, type = "right")
      ## FIX ME: interval censored
      #if(cen == "interval") y <- survival::Surv(y, ((y > censpoint[1]) * (y < censpoint[2])), type = "interval")
    }
  } else {
    if(survival::is.Surv(y)) {
      yc <- y
      y <- y[,1]
      # if(!("censored" %in% (strsplit(family$family.name, " ")[[1]]))) {
      # warning("response is Surv object but given distribution is not censored")
      # }
    }
  }
  
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
    
    rval <- suppressWarnings(try(family$ddist(y, eta, log = TRUE, weights = weights, sum = TRUE), silent = TRUE))
    if(inherits(rval, "try-error")) {
      print("try-error in ddist")
      return(1.7e+308)
    } else {
      if(any(is.na(rval))) {
        print("NAs in ddist")
        return(1.7e+308)
      } else {
        if(any(abs(rval) == Inf)) {
          print(c("Infinity in ddist"))
          return(1.7e+308)
        } else {
          nloglik <- - rval
          if(abs(nloglik) == Inf) return(1.7e+308)
          return(nloglik)
        }
      }
    }
  }  
  
  
  ## set up gradient
  grad <- function(eta) {
    gr <- - family$sdist(y, eta, weights = weights, sum = TRUE)
    if(any(is.na(gr))) {
      gr[is.na(gr)] =  1.7e+308
      print(c(eta,"gradient = NaN"))
    } else {
      if(any(gr == -Inf)) {
        gr[gr==-Inf] = -1.7e+308
        print(c(eta,"gradient = -Inf"))
      }
      if(any(gr ==  Inf)) {
        gr[gr== Inf] =  1.7e+308
        print(c(eta,"gradient = Inf"))
      }
    }
    return(gr)
  }
  
  
  ## calculate initial values if necessary or otherwise transform initial values for the distribution parameters to initial values on the link scale
  if(is.null(start) && is.null(start.eta)){
    starteta <- family$startfun(y, weights = weights)
    if(any(is.na(starteta))) {
      print(y)
      print(weights)
      print("start = NaN")
    # starteta <- family$startfun(y = rep(y, round(weights)))
    }
  } else {
    if(!(is.null(start))){
      # check that those parameters with a log-function in the linkfun are not negative
      if(start[(family$link == "log" | family$link == "logit")]<0) 
        start[(family$link == "log" | family$link == "logit")] <- 1e-10
      starteta <- family$linkfun(start)
    }
    if(!(is.null(start.eta))){
      starteta <- start.eta
    }
    if(any(is.na(starteta))) {
      warning("NAs in starteta")
      print(y)
      print(weights)
      print(start)
    }
  }
  
  if(!family$mle) {
    ## optimize negative log-likelihood
    opt <- try(optim(par = starteta, fn = nll, gr = grad, method = "L-BFGS-B",
                 hessian = (type.hessian == "numeric"), control = ocontrol))
    if(inherits(opt, "try-error")) {
      #print("opt try-error with gr=grad, L-BFGS-B")
      opt <- try(optim(par = starteta, fn = nll, gr = grad, method = "BFGS",
                       hessian = (type.hessian == "numeric"), control = ocontrol))
      if(inherits(opt, "try-error")) {
        #print("opt try-error with gr = grad, BFGS")
        print(starteta)
        opt <- optim(par = starteta, fn = nll, gr = grad, method = "Nelder-Mead",
                         hessian = (type.hessian == "numeric"), control = ocontrol)
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
                         hessian = TRUE, control = ocontrol)$hessian)
        if(inherits(nhess, "try-error")) {
          #print("opt try-error in hess with gr=grad, L-BFGS-B")
          nhess <- try(optim(par = eta, fn = nll, gr = grad, method = "BFGS",
                       hessian = TRUE, control = ocontrol)$hessian)
          if(inherits(nhess, "try-error")) {
            #print("opt try-error in hess with gr=grad, BFGS")
            nhess <- optim(par = eta, fn = nll, method = "BFGS",
                               hessian = TRUE, control = ocontrol)$hessian
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
    vc <- try(solve(-hess), silent = TRUE)
    if(inherits(vc, "try-error")) {
      vc <- try(qr.solve(-hess), silent = TRUE)
      if(inherits(vc, "try-error")) {
        vc <- try(chol2inv(chol(-hess)))
        if(inherits(vc, "try-error")) {
          print(-hess)
          print("hessian matrix is 'numerically' singular")
        }
      }
    }
    

    # if(inherits(vc, "try-error")) {
    #   vc <- MASS::ginv(-hess)
    # }
    
    vc <- as.matrix(vc)
    colnames(vc) <- rownames(vc) <- colnames(hess)
    
  } else {
    hess <- NULL
    vc <- NULL
  }
  
  
  ## estfun
  # each column represents one distribution parameter (1.col -> dldm * dmdpar = "dldeta.mu", 2.col -> dldd * dddpar = "dldeta.sigma", ...)
  # (first-order partial derivatives of the distribution function)
  if(estfun) {
    # estfun for link coefficients eta
    ef <- weights * family$sdist(y, eta, sum = FALSE)   ## FIX ME: cut out rows with weight = 0? -> No! index is of importance for independence tests (relation to covariates)
  } else {
    ef <- NULL                    
  }
  
  
  ##### additional functions with set parameters
  
  ## density function
  if(inherits(cl$family, "gamlss.family") && ("censored" %in% strsplit(family$family.name, " ")[[1]])){
    ddist <- function(x, log = FALSE) {
      if(!survival::is.Surv(x)){
        if(cens == "left") family$ddist(survival::Surv(x, x > censpoint, type = "left"), eta = eta, log = log)
        if(cens == "right") family$ddist(survival::Surv(x, x < censpoint, type = "right"), eta = eta, log = log)
        ## FIX ME: interval censored
      } else family$ddist(x, eta = eta,  log=log)
    }
  } else ddist <- function(x, log = FALSE) family$ddist(x, eta = eta, log = log)
  
  ## cumulative distribution function
  pdist <- function(q, lower.tail = TRUE, log.p = FALSE) family$pdist(q, eta = eta, lower.tail = lower.tail, log.p = log.p)
  
  ## quantile function
  qdist <- function(p, lower.tail = TRUE, log.p = FALSE) family$qdist(p, eta = eta, lower.tail = lower.tail, log.p = log.p)
  
  ## random function
  rdist <- if(is.null(family$rdist)) NULL else function(n) family$rdist(n, eta = eta)

  
  

  
  
  ## return value 
  rval <- list(
    npar = length(par),
    y = y,
    ny = ny,
    weights = weights,
    family = family$family.name,
    starteta = starteta,
    start = start,
    opt = opt,
    par = par,
    eta = eta,
    hess = hess,
    vcov = vc,
    loglik = loglik,
    call = cl,
    estfun = ef,
    familylist = family,
    ddist = ddist,
    pdist = pdist,
    qdist = qdist,
    rdist = rdist
  )
  
  class(rval) <- "distfit"
  return(rval)
}





## print, summary?, predict?
nobs.distfit <- function(object, ...) {
  object$ny
}

coef.distfit <- function(object, type = "parameter" , ...) {
  if(type == "link") return(object$eta)
  if(type == "parameter") return(object$par)
  ## FIXME: else, warning
}

## FIX: censored logistic distribution
predict.distfit <- function(object, type = "response", OOB = FALSE, ...){
  # calculation of the expected value 
  # of the given distribution with the calculated parameters
  if(type == "response") {
    
    if("censored" %in% strsplit(object$family, " ")[[1]])
    {
      par <- coef(object, type = "parameter")
      expv <- par[1]
      #lat.expv <- par[1]
      #object$ddist() / object$pdist()
      return(expv)
    } 
    
    
    f <- function(x){x * object$ddist(x, log = FALSE)}
    expv <- try(integrate(f,-Inf, Inf), silent = TRUE)
    if(inherits(expv, "try-error")) {
      expv <- try(integrate(f,-Inf, Inf, rel.tol = 1e-03))
      if(inherits(expv, "try-error")) {
        expv <- try(integrate(f,-Inf, Inf, rel.tol = 1e-02))
        if(inherits(expv, "try-error")) {
          print("rel.tol had to be set to 0.1 to calculated expected values for predictions")
          #print(coef(object))
          expv <- integrate(f,-Inf, Inf, rel.tol = 1e-01)
        }
      }
    }
    return(expv[[1]])
  }
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
    system.time(m1 <- crch(rain ~ 1, data = RainIbk, left = 0, dist = "logistic"))
    
    library("gamlss.cens")
    gen.cens(LO, type = "left")
    #RainIbk$rains <- Surv(RainIbk$rain, RainIbk$rain > 0, type = "left")
    system.time(m2 <- distfit(RainIbk$rain, family = LOlc, cens = "left", censpoint = 0))
    # FIX ME: calculation of starting values for censored distributions
    # m2 <- distfit(RainIbk$rains, family = LOlc, start = c(1,1))
    
    dist_list_cens_log <- make_censored_dist_list(dist = "logistic", type = "left", censpoint = 0)
    system.time(m3 <- distfit(RainIbk$rain, family = dist_list_cens_log))
    
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
  
  
  
  ## censored normal example
  gen.cens(NO, type = "left")
  y <- rnorm(500,-1,2)
  y <- pmax(y,0)
  system.time(m1 <- crch(y ~ 1, left = 0, dist = "gaussian"))
  system.time(m2 <- distfit(y, censpoint = 0, cens="left", family = NOlc))
  system.time(m3 <- distfit(y, family = dist_list_cens_normal))
  
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

