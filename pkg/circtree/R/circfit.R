circfit <- function(y, x = NULL, start = NULL, subset = NULL, na.action = NULL, 
                    weights = NULL, offset = NULL, ...,
                    vcov = TRUE, estfun = TRUE, object = FALSE, fit_control = circfit_control()) {

  # TODO: Why is circfit quite often called with object TRUE and what are the start values used for?!
    # if(object) cat("now object = TRUE\n")

  # Check unsupported arguments
  if(!(is.null(x) || NCOL(x) == 0L)) warning("no regression coefficients 
    are currently taken into account..")
  if(!is.null(subset)) warning("'subset' currently not supported and therefore not used")
  if(!is.null(na.action)) warning("'na.action' currently not supported and therefore not used")
  if(!is.null(offset)) warning("'offset' currently not supported and therefore not used")

  ## Convenience variables
  n <- NROW(y)
  allequy <- (length(unique(y)) == 1)
  cl <- match.call()
  
  ## Weights
  if(is.null(weights) || (length(weights)==0L)) weights <- as.vector(rep.int(1, n))
  if(unique(weights) == 0L) stop("weights are not allowed to be all zero")
  if(length(weights) != n) stop("number of observations and length of weights are not equal")

  ## Control parameters
  solve_kappa <- fit_control$solve_kappa
  useC <- fit_control$useC
  ncores <- fit_control$ncores
  fit_control$solve_kappa <- fit_control$useC <- fit_control$ncores <- NULL

  ## Get Von Mises family functions (for interchangeability the same distlist is used)
  circfam <- dist_vonmises(useC = useC, ncores = ncores) 

  ## MLE according to Bettina Gruen
  eta <- circfam$startfun(y, weights = weights, solve_kappa = solve_kappa)
  par <- circfam$linkinv(eta)

  ## Compute negative loglik
  nll <- -circfam$ddist(y, eta, log = TRUE, weights = weights, sum = TRUE)

  # FIXME: Check if estfun and vcov are correct (at least they are identical to Lisa's).
  if(estfun) {
    if(allequy) {
      ef <- matrix(0, ncol = length(eta), nrow = n)
    } else {
      ef <- as.matrix(weights * circfam$sdist(y, eta, sum = FALSE))
    }
  } else {
    ef <- NULL
  }

  if(vcov) {
    hess <- circfam$hdist(y, eta, weights = weights)

  ## Convert to matrix with nice column names
    hess <- as.matrix(hess)
    colnames(hess) <- rownames(hess) <- names(eta)

    vc <- try(solve(-hess), silent = TRUE)
    if(inherits(vc, "try-error")) {
      vc <- try(qr.solve(-hess), silent = TRUE)
      if(inherits(vc, "try-error")) {
        vc <- try(chol2inv(chol(-hess)))
        if(inherits(vc, "try-error")) {
          warning("hessian matrix is not invertible ('numerically' singular)")
        }
      }
    }
  
    ## Convert to matrix with nice column names
    vc <- as.matrix(vc)
    colnames(vc) <- rownames(vc) <- colnames(hess)

  } else {
    hess <- NULL
    vc <- NULL
  }

  # FIXME: Which model should I return, is there a mistake in the mobster?!
  if(vcov) object <- TRUE
  model <- list(
    npar = length(par),
    y = y,
    ny = n,
    weights = weights,
    family = circfam,
    start = start,
    par = par,
    eta = eta,
    hess = hess,
    vcov = vc,
    loglik = -nll,
    call = cl,
    estfun = ef
  )
  class(model) <- "circfit"
  
  list(coefficients = par,
       objfun = nll,
       estfun = ef,
       object = if(object) model else NULL)
}


## Different methods for circfit
# FIXME: All models go for the fit$object: good, bad?!

nobs.circfit <- function(object, ...) {
  return(object$ny)
}


coef.circfit <- function(object, type = "parameter" , ...) {
  if(type == "link") return(object$eta)
  if(type == "parameter") return(object$par)
}


predict.circfit <- function(object, type = "parameter", ...){
  if(type == "parameter") {
    return(object$par)
  } else {
    # FIXME: Implement type = 'response'
    stop("currently just type = 'parameter' implemented")
  }
}


vcov.circfit <- function(object, type = "link", ...) {
  if(type == "link") return(object$vcov)
  if(type == "parameter"){
    ## delta method
    delta.m <- diag(object$family$linkinvdr(object$eta))
    colnames(delta.m) <- rownames(delta.m) <- names(object$par)
    return(delta.m %*% object$vcov %*% delta.m)
  }
}


estfun.circfit <- function(object, ...) {
  return(object$estfun)
}


logLik.circfit <- function(object, ...) {
  structure(object$loglik, df = object$npar, class = "logLik")
}


bread.circfit <- function(object, type = c("parameter", "link"), ...) {
  type <- match.arg(type)
  if(type == "parameter") return(vcov(object, type = "parameter") * object$ny)
  if(type == "link") return(object$vcov * object$ny)
}


print.circfit <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Fitted ", x$family$family.name, "\n\n")
  if(length(x$par)) {
    cat("Distribution parameter(s):\n")
    print.default(format(x$par, digits = digits), print.gap = 2, quote = FALSE)
    cat("\n")
  } else {
    cat("No parameters \n\n")
  }
  cat(paste("Log-likelihood: ", format(x$loglik, digits = digits), "\n", sep = ""))
  if(length(x$npar)) {
    cat(paste("Df: ", format(x$npar, digits = digits), "\n", sep = ""))
  }
  cat("\n")
  
  invisible(x)
}
