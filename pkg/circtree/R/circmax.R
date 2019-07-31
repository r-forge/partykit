circmax <- function(formula, data, subset, na.action,
  model = TRUE, y = TRUE, x = FALSE,
  control = circmax_control(...), ...) {    # '...' arguments go also into circmax_control()

  ## Call
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)    # reorder arguments
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE   # dropping unused factor variables in model.frame

  ## Formula
  oformula <- as.formula(formula)
  formula <- Formula::as.Formula(formula)
  if(length(formula)[2L] < 2L) {    # formula consits of lhs (=repsonse) and rhs (=linear predictors)
    formula <- Formula::as.Formula(formula(formula), ~ 1)
  } else {
    if(length(formula)[2L] > 2L) {
      formula <- Formula::Formula(formula(formula, rhs = 1L:2L))
      warning("formula must not have more than two RHS parts")
    }
  }
  mf$formula <- formula

  ## Evaluate model.frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## Extract terms, model matrix, response
  mt <- terms(formula, data = data)
  mtX <- terms(formula, data = data, rhs = 1L)
  mtZ <- delete.response(terms(formula, data = data, rhs = 2L))
  Y <- model.response(mf, "numeric")
  X <- model.matrix(mtX, mf)
  Z <- model.matrix(mtZ, mf)

  ## Convert response to values between 0 and 2pi
  Y <- Y %% (2 * pi)

  ## Sanity check
  if(length(Y) < 1) stop("empty model")
  n <- length(Y)

  ## Call the actual workhorse: circmax_fit()
  rval <- circmax_fit(X, Y, Z, control)

  ## Further model information
  rval$call <- cl
  rval$formula <- oformula
  rval$terms <- list(location = mtX, concentration = mtZ, full = mt)
  rval$levels <- list(location = .getXlevels(mtX, mf), 
    concentration = .getXlevels(mtZ, mf), full = .getXlevels(mt, mf))
  rval$contrasts <- list(location = attr(X, "contrasts"), concentration = attr(Z, "contrasts"))
  if(model) rval$model <- mf
  if(y) rval$y <- Y
  if(x) rval$x <- list(location = X, concentration = Z)
  class(rval) <- "circmax"
  return(rval)
}

circmax_control <- function(maxit = 5000, start = NULL, method = "Nelder-Mead", 
  solve_kappa = "Newton-Fourier", gradient = FALSE, hessian = TRUE, ...) {

  ## Extract control arguments
  ctrl <- c(
    list(maxit = maxit, start = start, method = method, gradient = gradient, 
      solve_kappa = solve_kappa, hessian = hessian), 
    list(...)
  )
  if(!is.null(ctrl$fnscale)) warning("fnscale must not be modified")
  ctrl$fnscale <- 1
  if(is.null(ctrl$reltol)) ctrl$reltol <- .Machine$double.eps^(1/1.2)
  ctrl
}

circmax_fit <- function(x, y, z = NULL, control) {

  ## Dimensions
  n <- length(y)
  if(is.null(z)) z <- matrix(1, n, 1, dimnames = list(rownames(x), "(Intercept)"))
  m <- ncol(x)  
  p <- ncol(z)
  stopifnot(n == nrow(x), n == nrow(z))

  ## clean up control arguments
  method <- control$method
  gradient <- control$gradient
  solve_kappa <- control$solve_kappa
  hessian <- control$hessian
  control$method <- control$gradient <- control$solve_kappa <- control$hessian <- NULL

  ## Negative log-likelihood and gradients   
  nll <- function(par) {
    beta <- par[1:m]
    gamma <- par[m + (1:p)]
    mu <- 2 * atan(x[, 1, drop = FALSE] %*% beta[1]) + 2 * atan(x[, -1, drop = FALSE] %*% beta[-1])
    kappa <- exp(z %*% gamma)
    ll <- dvonmises(y, mu = mu, kappa = kappa, log = TRUE)
    -sum(ll)
  }

  gr <- function(par) {
    beta <- par[1:m]
    gamma <- par[m + (1:p)]
    mu0 <- 2 * atan(x[, 1, drop = FALSE] %*% beta[1])
    if(m > 1){
      mu <- 2 * atan(x[, -1, drop = FALSE] %*% beta[-1])
    } else { 
      mu <- 0
    }
    kappa <- exp(z %*% gamma)
    
    gr_mu0 <- 2 * kappa * sin(y - mu0 - mu) / (tan(mu0 / 2)^2 + 1) 
    if(m > 1){
      gr_mu <- 2 * kappa * sin(y - mu0 - mu) / (tan(mu / 2)^2 + 1) 
    }
    gr_kappa <- kappa * (cos(y - mu0 - mu) - 
      besselI(kappa, nu = 0, expon.scaled = TRUE) /besselI(kappa, nu = 1, expon.scaled = TRUE))
    gr <- numeric(length(par))
    gr[1] <- t(x[, 1, drop = FALSE]) %*% gr_mu0
    if(m > 1){
      gr[2:m] <- t(x[, -1, drop = FALSE]) %*% gr_mu
    }
    gr[m + (1:p)] <- t(z) %*% gr_kappa
    return(-gr)
  }

  ## Starting values (by default zeros)
  if(is.null(control$start)) {

    ## MLE according to Bettina Gruen
    circfam <- dist_vonmises()
    eta <- circfam$startfun(y, weights = NULL, solve_kappa = solve_kappa)
    
    start <- rep(0, m + p)
    start[1] <- eta[1] 
    start[m + 1] <- eta[2] 
  } else {
    start <- control$start
    stopifnot(length(start) == m + p)
  }
  control$start <- NULL

  ## Optimization
  if(gradient){ 
    if(method != "L-BFGS") warning("switched to method 'L-BFGS-B' to take into 'gradient = TRUE'")
    opt <- optim(par = start, fn = nll, control = control, method = "L-BFGS-B", gr = gr, hessian = hessian)
  } else { 
    opt <- optim(par = start, fn = nll, method = method, control = control, hessian = hessian)
  }

  ## Collect information
  names(opt)[1:2] <- c("coefficients", "loglik")
  opt$coefficients <- list(
    location = opt$coefficients[1:m],
    concentration = opt$coefficients[m + 1:p]
  )

  names(opt$coefficients$location) <- colnames(x)
  names(opt$coefficients$concentration) <- colnames(z)
  opt$loglik <- -opt$loglik
  opt$nobs <- n
  opt$df <- m + p
  
  return(opt)
}

logLik.circmax <- function(object, ...) {
  structure(object$loglik, df = object$df, class = "logLik")
}

coef.circmax <- function(object, model = c("full", "location", "concentration"), ...) {
  model <- match.arg(model)
  cf <- object$coefficients
  switch(model,
    "location" = {
      cf$location
    },
    "concentration" = {
      cf$concentration
    },
    "full" = {
      structure(c(cf$location, cf$concentration),
        .Names = c(names(cf$location), paste("(concentration)", names(cf$concentration), sep = "_")))
    }
  )
}

print.circmax <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Maximum likelihood estimation for the von Mises distribution\n\n")
  if(x$convergence > 0) {
    cat("Model did not converge\n")
  } else {
    if(length(x$coefficients$location)) {
      cat("Coefficients (location model with tanhalf link):\n")
      print.default(format(x$coefficients$location, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    } else {
      cat("No coefficients (in location model)\n\n")
    }
    if(length(x$coefficients$concentration)) {
      cat("Coefficients (concentration model (density kappa) with log link):\n")
      print.default(format(x$coefficients$concentration, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    } else {
      cat("No coefficients (in concentration model)\n\n")
    }
    cat(paste("Log-likelihood: ", format(x$loglik, digits = digits), "\n", sep = ""))
    if(length(x$df)) {
      cat(paste("Df: ", format(x$df, digits = digits), "\n", sep = ""))
    }
    cat("\n")
  }

  invisible(x)
}

terms.circmax <- function(x, model = c("location", "concentration", "full"), ...) x$terms[[match.arg(model)]]

model.frame.circmax <- function(formula, ...) {
  if(!is.null(formula$model)) return(formula$model)
  formula$terms <- formula$terms$full
  formula$call$formula <- formula$formula <- formula(formula$terms)
  NextMethod()
} 

model.matrix.circmax <- function(object, model = c("location", "concentration"), ...) {
  model <- match.arg(model)
  rval <- if(!is.null(object$x[[model]])) object$x[[model]]
    else model.matrix(object$terms[[model]], model.frame(object), contrasts = object$contrasts[[model]])
  return(rval)
}

predict.circmax <- function(object, newdata = NULL,
  type = c("location", "concentration", "parameter"),
  na.action = na.pass, ...) {

  ## Types of prediction
  type <- match.arg(type)

  ## Obtain model.frame/model.matrix
  tnam <- switch(type,
    "location" = "location",
    "concentration" = "concentration",
    "full")  
  if(is.null(newdata)) {
    X <- model.matrix(object, model = "location")
    Z <- model.matrix(object, model = "concentration")
  } else {
    mf <- model.frame(delete.response(object$terms[[tnam]]), newdata, na.action = na.action, xlev = object$levels[[tnam]])
    if(type != "concentration") X <- model.matrix(delete.response(object$terms$location), mf, contrasts = object$contrasts$location)
    if(type != "location") Z <- model.matrix(object$terms$concentration, mf, contrasts = object$contrasts$concentration)
  }

  ## Predicted parameters
  if(type != "concentration") location <- drop(2 * atan(X[, 1, drop = FALSE] %*% object$coefficients$location[1]) 
                                       + 2 * atan(X[, -1, drop = FALSE] %*% object$coefficients$location[-1]))
  if(type != "location") concentration <- exp(drop(Z %*% object$coefficients$concentration))

  ## Convert location to values between 0 and 2pi
  if (type != "concentration") {
    location <- location %% (2 * pi)
  }

  ## Compute result
  rval <- switch(type,
    "location" = location,
    "concentration" = concentration,
    "parameter" = data.frame(location, concentration)
  )
  return(rval)
}


estfun.circmax <- function(x, ...){
  # FIXME: Check if that is correct.

  ## Observed data and fit
  if(is.null(x$y) || is.null(x$x)) {
    mf <- model.frame(x)
    x$y <- model.response(mf)
    x$x <- list(
      "location" = model.matrix(x$terms$location, mf),
      "concentration" = model.matrix(x$terms$concentration, mf)
    )
  }

  ## Calculate distribution parameters
  mu <- drop(2 * atan(x$x$location[, 1, drop = FALSE] %*% x$coefficients$location[1]) 
    + 2 * atan(x$x$location[, -1, drop = FALSE] %*% x$coefficients$location[-1]))
  kappa <- exp(drop(x$x$concentration %*% x$coefficients$concentration))

  ## Calculate scores
  rval <- cbind(
    drop(2 * kappa * sin(x$y - mu) / ((tan(mu/2))^2 + 1) ) * x$x$location[, , drop = FALSE],
    drop(kappa * (cos(x$y - mu) -
      besselI(kappa, nu = 1, expon.scaled = TRUE) / besselI(kappa, nu = 0, expon.scaled = TRUE))) *
      x$x$concentration[, , drop = FALSE]
  )

  ## Convert to matrix with nice column names
  rval <- as.matrix(rval)
  colnames(rval) <- c(colnames(x$x$location), 
    paste("(concentration)", colnames(x$x$concentration), sep = "_"))

  return(rval)
}


vcov.circmax <- function(object, ...){
  # FIXME: Check if that is correct, and include analytical hessian.

   ## Observed data and fit
  if(is.null(object$y) || is.null(object$x)) {
    mf <- model.frame(object)
    object$y <- model.response(mf)
    object$x <- list(
      "location" = model.matrix(object$terms$location, mf),
      "concentration" = model.matrix(object$terms$concentration, mf)
    )
  }

  if (!is.null(object$hessian)){
    rval <- solve(as.matrix(object$hessian))
  } else {
    stop("Restart optimization with 'hessian = TRUE'")
    ### Calculate distribution parameters
    #mu <- drop(2 * atan(object$x$location[, 1, drop = FALSE] %*% object$coefficients$location[1]) 
    #  + 2 * atan(object$x$location[, -1, drop = FALSE] %*% object$coefficients$location[-1]))
    #kappa <- exp(drop(object$x$concentration %*% object$coefficients$concentration))
  
    ### Calculate hessian
    #hessian <- cbind(
    #)
    #rval <- solve(as.matrix(object$hessian))
  }

  ## Convert to matrix with nice column names
  rval <- as.matrix(rval)
  colnames(rval) <- c(colnames(object$x$location), 
    paste("(concentration)", colnames(object$x$concentration), sep = "_"))

  return(rval)
}

circmax_simulate <- function(n = 1000, beta = c(3, 5, 2), gamma = c(3, 3), seed = 111) {
  set.seed(seed)
  m <- length(beta) - 1  # here: number of betas minus intercept
  p <- length(gamma) -1  # here: number of gammas minus intercept

  #d <- sapply(1:(m + p), function(x) rnorm(n, runif(1, 0, 2 * pi), 0.2))
  #d <- sapply(1:(m + p), function(x) runif(n, 0, 1))
  d <- sapply(1:(m + p), function(x) rnorm(n, 0, 0.2))
  colnames(d) <- paste0("x", 1:(m + p))

  mu <- 2 * atan(beta[1]) + 2 * atan(crossprod(t(d[, 1:m, drop = FALSE]), beta[-1]))
  kappa <- exp(gamma[1] + crossprod(t(d[, m + 1:p, drop = FALSE]), gamma[-1]))
  d <- data.frame(d)
  d$y <- NULL
  for(i in 1:n) {
    d[i, "y"] <- circular::rvonmises(1, mu = circular::circular(mu[i]), kappa = kappa[i])
  }
  return(d)
}


