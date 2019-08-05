
## Density Von Mises 
## CAUTION: Different to dist_vonmises$dd(), as the latter takes only 2-values (no matrices!!)
dvonmises <- function(y, mu, kappa, log = FALSE) {
   if (any(kappa < 0)) {
      stop("kappa must be non-negative")
   }
   be <- besselI(kappa, nu = 0, expon.scaled = TRUE)
   out <- - log(2 * pi * be) + kappa * (cos(y - mu) - 1)
   if (!log) {
      out <- exp(out)
   }
   out
}


## Family for bamlss
# scaled besselI ( times exp(-kappa) in density, not in score)  and with hessian
vonmises_bamlss <- function(...) {
   f <- list(
         "family" = "vonmises",
         "names" = c("mu", "kappa"),
         "links" = c("tanhalf", "log"),
         "d" = function(y, par, log = FALSE) {
            if (any(par$kappa < 0)) {
               stop("kappa must be non-negative")
            }
            be <- besselI(par$kappa, nu = 0, expon.scaled = TRUE)
            out <- - log(2 * pi * be) + par$kappa * (cos(y - par$mu) - 1)
            if (!log) {
               out <- exp(out)
            }
            out
         },
         "score" = list(
            "mu" = function(y, par, ...) {
               drop(2 * par$kappa * sin(y - par$mu) / ((tan(par$mu/2))^2 + 1) )
            },
            "kappa" = function(y, par, ...) {
               drop(par$kappa * (cos(y - par$mu)  
                    - besselI(par$kappa, nu = 1, expon.scaled = TRUE) / besselI(par$kappa, nu = 0, expon.scaled = TRUE)))
            }
         ),
         "hess" = list(
            "mu" = function(y, par, ...) {
               ta <- tan(par$mu/2)
               drop(4 * par$kappa / (ta^2 + 1)^2 * (sin(y - par$mu) * ta + cos(y - par$mu)))
            },
            "kappa" = function(y, par, ...) { 
               be0 <- besselI(par$kappa, nu = 0, expon.scaled = TRUE)
               be1 <- besselI(par$kappa, nu = 1, expon.scaled = TRUE)
               be2 <- besselI(par$kappa, nu = 2, expon.scaled = TRUE)
               drop(-par$kappa * (cos(y - par$mu) - be1 / be0) - par$kappa^2 * ((2 * be1^2 - be0^2 - 
                 be2 * be0) / (2 * be0^2)))
            }
         )
   )
   class(f) <- "family.bamlss"
   return(f)
}


## Distlist Von Mises Distribution

dist_vonmises <- function(useC = FALSE, ncores = 1) {
  parnames <- c("mu", "kappa")
  etanames <- c("tan(mu/2)", "log(kappa)")

  # Setting the logical flag here, can be evaluated
  # inside the different methods using prent frame if
  # object is not ripped appart.
  useC   <- useC
  stopifnot(is.logical(useC))
  stopifnot(is.numeric(ncores))
  ncores <- floor(ncores)
  if ( ncores < 1 ) stop("ncores argument has to be a positive integer.")

  # Von Mises distribution function. Can be used to
  # retrieve densities, log-densities, and log-likelihood sum.
  ddist <- function(y, eta, log = TRUE, weights = NULL, sum = FALSE) {

    # Use R-version of the density function.
    if ( ! useC ) {
       # par <- linkinv(eta)  # CAUTION: Lisa is evaluating eta, as this simplifies the equation and, therefore, fastens the evaluation 
       par <- cbind(2 * atan(eta[1]), exp(eta[2]))  ## use cbind such that par is a data.frame and not a list

       if (any(par[2] < 0, na.rm = TRUE)) {
         stop("kappa must be non-negative")
       }
       be <- besselI(par[2][[1]], nu = 0, expon.scaled = TRUE)  ## NOTE: kappa needs to be numeric, but par[2] is still a data.frame
       val <- - log(2 * pi * be) + par[2] * (cos(y - par[1]) - 1)
       if (!log) {
         val <- exp(val)
       }

       if(sum) {
         if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y))  # CAUTION: All the weights are just copied from Lisa
         val <- sum(weights * val, na.rm = TRUE)
       }

       return(val)

    # Else use the C function.
    } else {
       if(is.null(weights) || (length(weights)==0L)) weights <- rep(1., length(y))
       .Call("circ_density", as.numeric(y), as.numeric(eta), as.numeric(weights),
             as.integer(log[1L]), as.integer(sum[1L]), as.integer(ncores[1L]), PACKAGE = "circtree")
    }

  }

  # Von Mises score function.
  sdist <- function(y, eta, weights = NULL, sum = FALSE) {

    # Use the score function
    if ( ! useC ) {

      #par <- linkinv(eta)
      par <- c(2 * atan(eta[1]), exp(eta[2]))
      # CAUTION: this is the score as in bamlss: d(l)/d(eta) = d(l)/d(mu) * d(mu)/d(eta) and d(l)/d(eta) = d(l)/d(kappa) * d(mu)/d(kappa)
      score <- cbind(drop(2 * par[2] * sin(y - par[1]) / ((tan(par[1]/2))^2 + 1) ),
                     drop(par[2] * (cos(y - par[1])
                      - besselI(par[2], nu = 1, expon.scaled = TRUE) / besselI(par[2], nu = 0, expon.scaled = TRUE))))

      score <- as.matrix(score)
      colnames(score) <- etanames
      if(sum) {
        if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y))
        # if score == Inf replace score with 1.7e308 because Inf*0 would lead to NaN -> gradient is NaN
        score[score==Inf] = 1.7e308
        score <- colSums(weights * score, na.rm = TRUE)
      }

    # Else use the C function
    } else {
       score <- .Call("circ_score", as.numeric(y), as.numeric(eta), as.numeric(weights),
             as.integer(sum[1L]), as.integer(ncores[1L]), PACKAGE = "circtree")
       colnames(score) <- etanames
       if ( nrow(score) == 1 ) score <- score[1,]
    }

    return(score)
  }

  hdist <- function(y, eta, weights = NULL) {

    #par <- linkinv(eta)
    par <- c(2 * atan(eta[1]), exp(eta[2]))

    if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y))

    # CAUTION: Lisa's approach is different to Niki's. The whole hessian matrix is constructed. 
    #          For bamlss we just need d^2(l)/d(eta.mu)^2 and d^2(l)/d(eta.sigma)^2.
    #          Also not completely clear over which dimension and why is summed?
    #calc d^2(l)/d(eta.mu)^2
    ta <- tan(par[1]/2)
    d2ld.etamu2 <- sum(weights * drop(-4 * par[2] / (ta^2 + 1)^2 * (sin(y - par[1]) * ta + cos(y - par[1]))), na.rm = TRUE)

    #calc d^2(l)/d(eta.mu eta.sigma)

    d2ld.etamu.d.etasigma <- sum(weights * drop(2 * par[2] * sin(y - par[1]) / ((tan(par[1]/2))^2 + 1) ), na.rm = TRUE)
          # should be 0 for exact parameters (here: observed hess)

    # calc d^2(l)/d(eta.sigma)^2
    be0 <- besselI(par[2], nu = 0, expon.scaled = TRUE)
    be1 <- besselI(par[2], nu = 1, expon.scaled = TRUE)
    be2 <- besselI(par[2], nu = 2, expon.scaled = TRUE)
    d2ld.etasigma2 <- sum(weights * drop(par[2] *
      (cos(y - par[1]) - be1 / be0) + par[2]^2 * ((2 * be1^2 - be0^2 - be2 * be0) / (2 * be0^2))), na.rm = TRUE)

    # par <- c(eta[1], exp(eta[2]))                           
    # d2ld.etamu2 <- sum(weights * rep.int(-1/par[2]^2, ny))
    # d2ld.etamu.d.etasigma <- sum(weights * (-2)*(y-par[1])/par[2]^2), na.rm = TRUE)          # should be 0 for exact parameters (here: observed hess)
    # d2ld.etasigma2 <- sum(weights * (-2)*(y-par[1])^2/par[2]^2, na.rm = TRUE)         

    hess <- matrix(c(d2ld.etamu2, d2ld.etamu.d.etasigma, d2ld.etamu.d.etasigma, d2ld.etasigma2), nrow = 2)
    colnames(hess) <- rownames(hess) <-  etanames

    return(hess)
  }

  # CAUTION: Should the inverse links be on the parameters? And do the circular calls properly work?!

  pdist <- function(q, eta) circular::pvonmises(q, mu = 2 * atan(eta[1]), kappa = exp(eta[2]))
  qdist <- function(p, eta) circular::qvonmises(q, mu = 2 * atan(eta[1]), kappa = exp(eta[2]))
  rdist <- function(n, eta) circular::rvonmises(n, mu = 2 * atan(eta[1]), kappa = exp(eta[2]))

  # CAUTION: is the name important
  link <- c("tanhalf", "log")

  linkfun <- function(par) {
    eta <- c(tan(par[1] / 2), log(par[2]))
    names(eta) <- etanames
    return(eta)
  }

  linkinv <- function(eta) {
    par <- c(2 * atan(eta[1]), exp(eta[2]))
    names(par) <- parnames
    return(par)
  }

  # CAUTION: seems to be the same as mu.eta = the derivative of the inverse link (check!)
  linkinvdr <- function(eta) {
    dpardeta <- c(2 / (eta[1]^2 + 1), exp(eta[2]))
    names(dpardeta) <- parnames
    return(dpardeta)
  }

  #startfun <- function(y, weights = NULL){
  #  mu <- 0
  #  kappa <- 1

  #  starteta <- c(tan(mu/2), log(kappa))
  #  names(starteta) <- etanames
  #  return(starteta)
  #}

  startfun <- function(y, weights = NULL, solve_kappa = "Newton-Fourier") {
    x <- cbind(cos(y), sin(y))
    if (is.null(weights) || (length(weights)==0L)) {
        xbar <- colMeans(x)
    } else {
        xbar <- colMeans(weights * x) / mean(weights)
    }
    mu <- atan(xbar[2] /xbar[1]) + (xbar[1] < 0) * sign(xbar[2]) * pi
    rbar <- sqrt(sum(xbar^2))

    # Calling solver function (iteratively estimate kappa).
    if(solve_kappa == "Newton-Fourier"){
      kappa <- do.call(solve_kappa_Newton_Fourier, list(r = rbar, useC = useC, ncores = ncores))
    } else if(solve_kappa == "Uniroot"){
      kappa <- do.call(solve_kappa_uniroot, list(r = rbar))
    } else if(solve_kappa == "Banerjee_et_al_2005"){
      kappa <- do.call(solve_kappa_Banerjee_et_al_2005, list(r = rbar))
    }

    starteta <- c(tan(mu / 2), log(kappa))
    names(starteta) <- etanames
    return(starteta)
  }

  mle <- TRUE

  list(family.name = "Von Mises Distribution",
       ddist = ddist,
       sdist = sdist,
       hdist = hdist,
       pdist = pdist,
       qdist = qdist,
       rdist = rdist,
       link = link,
       linkfun = linkfun,
       linkinv = linkinv,
       linkinvdr = linkinvdr,
       startfun = startfun,
       mle = mle,
       gamlssobj = FALSE,
       censored = FALSE
  )

}

## MLE Von Mises Distribution

## Package movMF implements different kappa solvers:
## o Newton Fourier is used by default.
## o Uniroot seems to provide a safe option.
## o Banerjee_et_al_2005 provides a quick approximation.
solve_kappa_Newton_Fourier <- function (r, tol = 1e-06, maxiter = 100L, useC = FALSE, ncores = 1) {

    # Using C?
    if ( ! useC ) {

       lower <- movMF:::Rinv_lower_Amos_bound(r, 0)
       upper <- movMF:::Rinv_upper_Amos_bound(r, 0)
       iter <- 1L
       while (iter <= maxiter) {
           A <- movMF:::A(lower, 2)
           Aprime <- movMF:::Aprime(lower, 2, A = A)
           lower <- lower - (A - r)/Aprime
           A <- movMF:::A(upper, 2)
           upper <- upper - (A - r)/Aprime
           if ((upper - lower) < tol * (lower + upper)) {
               if ((upper - lower) < -tol * (lower + upper))
                   stop("no convergence")
               break
           }
           iter <- iter + 1L
       }
       return((lower + upper)/2)

    # Else using the C function
    } else {
        .Call("solve_kappa_Newton_Fourier", as.numeric(r),
              as.integer(maxiter), as.integer(ncores[1L]), PACKAGE = "circtree")
    }
}

solve_kappa_uniroot <- function(r, tol = 1e-06) {
    interval <- c(movMF:::Rinv_lower_Amos_bound(r, 0),
                  movMF:::Rinv_upper_Amos_bound(r, 0))
    if (abs(diff(interval)) < tol)
        mean(interval)
    else uniroot(function(kappa) movMF:::A(kappa, 2) - r,
                 interval = interval, tol = tol)$root
}

solve_kappa_Banerjee_et_al_2005 <- function(r) {
   r * (2 - r^2)/(1 - r^2)
}

#startfun <- function(y, weights = NULL, solve_kappa =
#solve_kappa_Newton_Fourier) {
#    x <- cbind(cos(y), sin(y))
#    if (is.null(weights) || (length(weights)==0L)) {
#        xbar <- colMeans(x)
#    } else {
#        xbar <- colMeans(weights * x) / mean(weights)
#    }
#    mu <- atan(xbar[2] /xbar[1]) + (xbar[1] < 0) * sign(xbar[2]) * pi
#    rbar <- sqrt(sum(xbar^2))
#    kappa <- solve_kappa(rbar)
#
#    starteta <- c(tan(mu / 2), log(kappa))
#    names(starteta) <- etanames
#    return(starteta)
#}

