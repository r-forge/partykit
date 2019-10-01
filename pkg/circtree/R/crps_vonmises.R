# -------------------------------------------------------------------
# - NAME:   crps_vonmises.R
# - AUTHOR: Moritz N. Lang, Lisa Schlosser
# - DATE:   2019-05-21
# -------------------------------------------------------------------
# - PURPOSE: Circular CRPS (von Mises) based on numeric integration using
#            the charististic equation
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2019-10-01 on thinkmoritz
# -------------------------------------------------------------------

##YEAR: 2017
##COPYRIGHT HOLDER: Ludmila Simkova
## Function cfC_vonMises() from Package 'CharFun'
cfC_vonMises <- function(t, mu = 0, kappa = 1){
  szt <- dim(t)
  t <- c(t)
  cf <- unlist(lapply(t, function(t) tryCatch((Bessel::BesselI(kappa, 
      abs(t), TRUE)/Bessel::BesselI(kappa, 0, TRUE)) * exp((0+1i) * 
      t * mu), error = function(e) 0)))
  cf[t == 0] <- 1
  dim(cf) <- szt
  return(cf)
}



## CRPS von Mises
crps_vonmises <- function(y, mu, kappa, sum = FALSE) {

  ## Perform some input checks
  if(any(y < -pi) || any(y > pi) || any(mu < -pi) || any(mu > pi) || any(kappa < 0 ))
    stop("y and mu must be in the interval of [-pi, pi], and kappa must be non negative!") 

  if(!inherits(y, c("numeric", "integer")) || !inherits(mu, c("numeric", "integer")) ||
    !inherits(kappa, c("numeric", "integer"))) {
    stop("Input 'y', 'mu', and 'kappa' must be numeric vectors...")
  }

  ## Create data.frame (fails if any has not the same length or length equal one)
  dat <- data.frame("y" = y, "mu" = mu, "kappa" = kappa)

  ## Get mu with smallest absolute distance
  idx <- apply(abs(matrix(dat$mu, ncol = 3, nrow = nrow(dat), byrow = FALSE) +
    matrix(c(-2 * pi, 0, 2 * pi), ncol = 3, nrow = nrow(dat), byrow = TRUE) - dat$y), 1, which.min)
  dat$mu <- dat$mu + c(-2 * pi, 0, 2 * pi)[idx]

  ## Calculate CRPS based on characteristic function  
  rval <- sapply(1:nrow(dat), function(i){

    int_fun <- function(x) {
      1 / pi * abs(cfC_vonMises(x, dat[i, "mu"], dat[i, "kappa"]) - 
        exp((0+1i) * x * dat[i, "y"]))^2 / x^2
    }

    crps <- integrate(int_fun, 0, Inf)$value
    return(crps)
  })

  if(sum) rval <- sum(rval)
  return(rval)
}
