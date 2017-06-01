dist_gaussian <- function() {

  parnames <- c("mu", "sigma")
  etanames <- c("mu", "log(sigma)")
  
  
  ddist <-  function(y, eta, log = TRUE, weights = NULL, sum = FALSE) {     
    
    val <- -1/2 * (log(2*pi) + 2*eta[2] + exp(log((y-eta[1])^2) - 2*eta[2]))
    if(!log) val <- exp(val)
    
    # par <- c(eta[1], exp(eta[2]))
    # val <- dnorm(y, mean = par[1], sd = par[2], log = log)
    
    if(sum) {
      if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y))
      val <- sum(weights * val, na.rm = TRUE)
    }
    return(val)
  }
  
  
  sdist <- function(y, eta, weights = NULL, sum = FALSE) {   
      
    score <- cbind(exp(-2*eta[2]) * (y-eta[1]), 
                   -1 + exp(-2*eta[2] + log((y-eta[1])^2)))
    
    # par <- c(eta[1], exp(eta[2])) 
    # score <- cbind(1/par[2]^2 * (y-par[1]), 
    #                (-1/par[2] + ((y - par[1])^2)/(par[2]^3)) * exp(eta[2]))
    
    score <- as.matrix(score)
    colnames(score) <- etanames
    if(sum) {
      if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y))
      # if score == Inf replace score with 1.7e308 because Inf*0 would lead to NaN -> gradient is NaN
      score[score==Inf] = 1.7e308
      score <- colSums(weights * score, na.rm = TRUE)
    }
    return(score)
  }
  
  
  hdist <- function(y, eta, weights = NULL) {    
    ny <- length(y)
    if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, ny)
    
    d2ld.etamu2 <- sum(weights * rep.int(-exp(-2*eta[2]), ny))
    d2ld.etamu.d.etasigma <- sum(weights * (-2)*(y-eta[1]) * exp(-2*eta[2]), na.rm = TRUE)          # should be 0 for exact parameters (here: observed hess)
    d2ld.etasigma2 <- sum(weights * (-2)*exp(log((y-eta[1])^2) - 2*eta[2]), na.rm = TRUE)    
    
    # par <- c(eta[1], exp(eta[2]))                           
    # d2ld.etamu2 <- sum(weights * rep.int(-1/par[2]^2, ny))
    # d2ld.etamu.d.etasigma <- sum(weights * (-2)*(y-par[1])/par[2]^2), na.rm = TRUE)          # should be 0 for exact parameters (here: observed hess)
    # d2ld.etasigma2 <- sum(weights * (-2)*(y-par[1])^2/par[2]^2, na.rm = TRUE)         
    
    hess <- matrix(c(d2ld.etamu2, d2ld.etamu.d.etasigma, d2ld.etamu.d.etasigma, d2ld.etasigma2), nrow = 2)
    colnames(hess) <- rownames(hess) <-  etanames
    
    return(hess)
  }
  
  
  ## additional functions pdist, qdist, rdist
  pdist <- function(q, eta, lower.tail = TRUE, log.p = FALSE) pnorm(q, mean = eta[1], sd = eta[2], 
                                                                            lower.tail = lower.tail, log.p = log.p)
  qdist <- function(p, eta, lower.tail = TRUE, log.p = FALSE) qnorm(p, mean = eta[1], sd = eta[2], 
                                                                            lower.tail = lower.tail, log.p = log.p)
  rdist <- function(n, eta) rnorm(n, mean = eta[1], sd = eta[2])
  
  
  link <- c("identity", "log")
  
  linkfun <- function(par) {
    eta <- c(par[1], log(par[2]))
    names(eta) <- etanames
    return(eta)
  }
  
  
  linkinv <- function(eta) {
    par <- c(eta[1], exp(eta[2]))
    names(par) <- parnames
    return(par)
  }
  
  
  linkinvdr <- function(eta) {
    dpardeta <- c(1, exp(eta[2]))
    names(dpardeta) <- parnames
    return(dpardeta)
  }
  

  startfun <- function(y, weights = NULL){
    if(is.null(weights) || (length(weights)==0L)) {
      mu <- mean(y)
      sigma <- sqrt(1/length(y) * sum((y - mu)^2))
    } else {
      mu <- weighted.mean(y, weights)
      sigma <- sqrt(1/sum(weights) * sum(weights * (y - mu)^2))
    }
    starteta <- c(mu, log(sigma))
    names(starteta) <- etanames
    return(starteta)
  }
  
  mle <- TRUE

  list(family.name = "Normal Distribution",
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
       mle = mle
  )

}



##------------------------------------------------------------------------------

## FIX ME: adapt initial values for Weibull distribution
###### dist_list for Weibull distribution
if(FALSE) {
  
  dist_list_weibull <- list()
  
  parnames <- c("mean", "scale")
  etanames <- c("mean", "log(scale)")
  
  
  ddist <-  function(y, eta, log = TRUE, weights = NULL, sum = FALSE) {     
    par <- c(eta[1], exp(eta[2]))
    
    if(any(par==Inf)) {
      print(c(eta,"initial values are Inf"))
      #par[2] <- 1
      # par <- c(2,1)
      return(1.7e308)   ### FIX
    }
    
    # val <- 1/(par[2]*exp(par[1])) * (y/exp(par[1]))^(1/par[2]-1) * exp(-(y/exp(par[1]))^(1/par[2]))
    # if(log) val <- log(val)
    
    val <- -eta[2] + expm1(-eta[2]) * log(y) - eta[1]*exp(-eta[2]) - exp(exp(-eta[2]) * (log(y) - eta[1]))
    # (y/exp(eta[1]))^exp(-eta[2]) = exp(exp(-eta[2]) * (log(y) - eta[1]))
    if(!log) val <- exp(val)
    
    # val <- survival::dsurvreg(y, mean = par[1], scale = par[2])
    # val <- dweibull(y, shape = 1/par[2], scale = exp(par[1]), log = log)
    
    # if(val<=0) print(c("negative val", eta))    ### FIX
    #if(any(val == -Inf)) print(c("infinite values in ddist",par, y))
    #if(any(is.na(val))) print(c("NAs in ddist", par, y))
    #if(any(val == -Inf)) print(c("infinite values in ddist",par))
    #if(any(is.na(val))) print(c("NAs in ddist", par))
    
    if(sum) {
      if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y))
      val <- sum(weights * val, na.rm = TRUE)
    }
    if(any(abs(val) == Inf) || any(is.na(val))) return(1.7e308)
    return(val)
  }
  
  
  sdist <- function(y, eta, weights = NULL, sum = FALSE) {   
       
    score_m <- exp(-eta[2]) * (exp(exp(-eta[2]) * (log(y) - eta[1])) -1) 
    score_s <- -1 - (log(y) - eta[1]) * (exp(-eta[2])) * (1 - exp(exp(-eta[2]) * (log(y) - eta[1])))
      
    # par <- c(eta[1], exp(eta[2]))  
    # score_m <- (1/par[2]) * (-1 + (y/exp(par[1]))^(1/par[2]))
    # score_s <- ((-1/par[2]) - (log(y)-par[1])/(par[2]^2) * (1 - (y/exp(par[1]))^(1/par[2]))) * par[2]
    
    score <- cbind(score_m, score_s)
    score <- as.matrix(score)
    colnames(score) <- etanames
    if(sum) {
      if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y))
      # if score == Inf replace score with 1.7e308 because Inf*0 would lead to NaN -> gradient is NaN
      score[score==Inf] = 1.7e308
      score <- colSums(weights * score, na.rm = TRUE)
    }
    # if(any(is.null(score)) || any(is.na(score))) print(c("gr is NULL or NA", score, "y", y, eta, weights))
    return(score)
  }
  
  
  hdist <- function(y, eta, weights = NULL) {    
    ny <- length(y)
    if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, ny)
    
    d2ld.etamu2 <- sum(weights * 
                         (-exp(-2*eta[2]) * exp(exp(-eta[2]) * (log(y) - eta[1]))), 
                       na.rm = TRUE) 
    d2ld.etamu.d.etasigma <- sum(weights * 
                                   exp(-eta[2]) * (1 - exp(exp(-eta[2]) * (log(y) - eta[1])) * (1 + exp(-eta[2]) * (log(y) - eta[1]))), 
                                 na.rm = TRUE)
    d2ld.etasigma2 <- sum(weights * 
                            (log(y) - eta[1]) * exp(-eta[2]) * (1 - exp(exp(-eta[2]) * (log(y) - eta[1])) * (1 + (log(y) - eta[1]) * exp(-eta[2]))),
                          na.rm = TRUE)
    
    # par <- c(eta[1], exp(eta[2]))                           
    # d2ld.etamu2 <- sum(weights * (-1/par[2]^2) * ((y/exp(par[1]))^(1/par[2])))
    # d2ld.etamu.d.etasigma <- sum(weights * (1/par[2]) * (1 - (y/exp(par[1]))^(1/par[2]) * (1 + (log(y)-par[1])/par[2])))
    # d2ld.etasigma2 <- sum(weights * 
    #                         ((log(y)-par[1])/par[2]) * 
    #                         (1 - (y/exp(par[1]))^(1/par[2]) * (1 + ((log(y)-par[1])/par[2]))))
    # d2ld.etasigma2 <- sum(weights * 
    #                         (log(y)-par[1])/par[2] - ((log(y)-par[1])/par[2]) * (y/exp(par[1]))^(1/par[2])
    #                       - ((log(y)-par[1])/par[2])^2 * (y/exp(par[1]))^(1/par[2]))
    
    hess <- matrix(c(d2ld.etamu2, d2ld.etamu.d.etasigma, d2ld.etamu.d.etasigma, d2ld.etasigma2), nrow = 2)
    colnames(hess) <- rownames(hess) <-  etanames
    
    return(hess)
  }
  
  
  

  ## additional functions pdist, qdist, rdist
  
    pdist <- function(q, eta, lower.tail = TRUE, log.p = FALSE) pweibull(q, shape = 1/exp(eta[2]), scale = exp(eta[1]), lower.tail = lower.tail, log.p = log.p)
    qdist <- function(p, eta, lower.tail = TRUE, log.p = FALSE) qweibull(p, shape = 1/exp(eta[2]), scale = exp(eta[1]), lower.tail = lower.tail, log.p = log.p)
    rdist <- function(n, eta) rweibull(n, shape = 1/exp(eta[2]), scale = exp(eta[1]))
  
  
  
  link <- c("identity", "log")
  
  linkfun <- function(par) {
    eta <- c(par[1], log(par[2]))
    names(eta) <- etanames
    return(eta)
  }
  
  
  linkinv <- function(eta) {
    par <- c(eta[1], exp(eta[2]))
    names(par) <- parnames
    return(par)
  }
  
  
  linkinvdr <- function(eta) {
    dpardeta <- c(1, exp(eta[2]))
    names(dpardeta) <- parnames
    return(dpardeta)
  }
  
  
  ## FIX ME: adapt starting values for Weibull distribution
  startfun <- function(y, weights = NULL){
    
    # the 0.632 quantile of the distribution is an estimator for lambda
    # with lambda being the scale parameter in the rweibull parametrization
    # -> mean = log(lambda)
    # lambda <- if(is.null(weights) || (length(weights)==0L)) quantile(y, p=0.632) else quantile(rep(y, round(weights)), p=0.632)
    # mean <- log(lambda)
    
    # using the initial values for the exponential distribution (scale_survreg = 1)
    # mean = log(1/lambda) with lambda being the rate parameter of the Exp. distribution
    # (lambda = 1/scale_dweibull)
    lambda <- if(is.null(weights) || (length(weights)==0L)) length(y)/sum(y) else sum(weights)/sum(weights * y)
    mean <- log(1/lambda)
    
    scale <- 1
    starteta <- c(mean, log(scale))
    names(starteta) <- etanames
    return(starteta)
  }
  
  mle <- FALSE
  
  dist_list_weibull <- list(family.name = "Weibull Distribution",
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
                           mle = mle
  )
}






###### dist_list for Poisson distribution
if(FALSE) {
  
  dist_list_poisson <- list()
  
  parnames <- c("mu")
  etanames <- c("log(mu)")
  
  
  ddist <-  function(y, eta, log = TRUE, weights = NULL, sum = FALSE) {     
    par <- exp(eta)
    #val <- par^y / gamma(y+1) * exp(-par)
    #if(log) val <- log(val)
    val <- dpois(x = y, lambda = par, log = log)
    if(sum) {
      if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y))
      val <- sum(weights * val)
    }
    return(val)
  }
  
  
  sdist <- function(y, eta, weights = NULL, sum = FALSE) {   
    par <- exp(eta)                           
    score <- (y - par)
    score <- as.matrix(score)
    colnames(score) <- etanames
    if(sum) {
      if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y))
      # if score == Inf replace score with 1.7e308 because Inf*0 would lead to NaN -> gradient is NaN
      score[score==Inf] = 1.7e308
      score <- colSums(weights * score)
    }
    return(score)
  }
  
  
  hdist <- function(y, eta, weights = NULL) {
    if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y))
    
    par <- exp(eta)                           
    hess <- rep(-par, length(y))
    hess <- as.matrix(sum(weights * hess))
    colnames(hess) <- rownames(hess) <-  etanames
    
    return(hess)
  }
  
  
  ## additional functions pdist, qdist, rdist
  {
    pdist <- function(q, eta, lower.tail = TRUE, log.p = FALSE) ppois(q, lambda = exp(eta[1]), lower.tail = lower.tail, log.p = log.p)
    qdist <- function(p, eta, lower.tail = TRUE, log.p = FALSE) qpois(p, lambda = exp(eta[1]), lower.tail = lower.tail, log.p = log.p)
    rdist <- function(n, eta) rpois(n, lambda = exp(eta[1]))
  }

  
  link <- c("log")
  
  linkfun <- function(par) {
    eta <- log(par)
    names(eta) <- etanames
    return(eta)
  }
  
  
  linkinv <- function(eta) {
    par <- exp(eta)
    names(par) <- parnames
    return(par)
  }
  
  
  linkinvdr <- function(eta) {
    dpardeta <- exp(eta)
    names(dpardeta) <- parnames
    return(dpardeta)
  }
  
  
  startfun <- function(y, weights = NULL){
    mu <- if(is.null(weights) || (length(weights)==0L)) mean(y) else weighted.mean(y, weights)
    starteta <- log(mu)
    names(starteta) <- etanames
    return(starteta)
  }
  
  
  mle <- TRUE
  
  
  dist_list_poisson <- list(family.name = "Poisson Distribution",
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
                            mle = mle
  )
}





###### dist_list for exponential distribution
if(FALSE) {
  
  dist_list_exp <- list()
  
  parnames <- c("lambda")
  etanames <- c("log(lambda)")
  
  
  ddist <-  function(y, eta, log = TRUE, weights = NULL, sum = FALSE) {     
    par <- exp(eta)
    #val <- par * exp(-par * y)
    #if(log) val <- log(val)
    val <- dexp(x = y, rate = par, log = log)
    if(sum) {
      if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y))
      val <- sum(weights * val)
    }
    return(val)
  }
  
  
  sdist <- function(y, eta, weights = NULL, sum = FALSE) {
    par <- exp(eta)                           
    score <- 1 - y * par
    score <- as.matrix(score)
    colnames(score) <- etanames
    if(sum) {
      if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y))
      # if score == Inf replace score with 1.7e308 because Inf*0 would lead to NaN -> gradient is NaN
      score[score==Inf] = 1.7e308
      score <- colSums(weights * score)
    }
    return(score)
  }
  
  
  hdist <- function(y, eta, weights = NULL) {    
    if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y))
    
    par <- exp(eta)                           
    hess <- -y * par
    hess <- as.matrix(sum(weights * hess))
    colnames(hess) <- rownames(hess) <-  etanames
    
    return(hess)
  }
  

  ## additional functions pdist, qdist, rdist
  {
    pdist <- function(q, eta, lower.tail = TRUE, log.p = FALSE) pexp(q, rate = exp(eta), lower.tail = lower.tail, log.p = log.p)
    qdist <- function(p, eta, lower.tail = TRUE, log.p = FALSE) qexp(p, rate = exp(eta), lower.tail = lower.tail, log.p = log.p)
    rdist <- function(n, eta) rexp(n, rate = exp(eta))
  }
  
  
  link <- c("log")
  
  linkfun <- function(par) {
    eta <- log(par)
    names(eta) <- etanames
    return(eta)
  }
  
  
  linkinv <- function(eta) {
    par <- exp(eta)
    names(par) <- parnames
    return(par)
  }
  
  
  linkinvdr <- function(eta) {
    dpardeta <- exp(eta)
    names(dpardeta) <- parnames
    return(dpardeta)
  }
  
  
  startfun <- function(y, weights = NULL){
    lambda <- if(is.null(weights) || (length(weights)==0L)) length(y)/sum(y) else sum(weights)/sum(weights * y)
    starteta <- log(lambda)
    names(starteta) <- etanames
    return(starteta)
  }
  
  mle <- TRUE
  
  
  dist_list_exp <- list(family.name = "Exponential Distribution",
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
                        mle = mle
  )
}





###### dist_list for gamma distribution
if(FALSE) {
  
  dist_list_gamma <- list()
  
  parnames <- c("shape", "scale")
  etanames <- c("log(shape)", "log(scale)")
  
  
  ddist <-  function(y, eta, log = TRUE, weights = NULL, sum = FALSE) {     
    par <- c(exp(eta[1]), exp(eta[2]))
    # val <- -par[1] * log(par[2]) + (par[1]-1) * log(y) - y / par[2] - lgamma(par[1])
    # if(!log) val <- exp(val)
    val <- dgamma(x = y, shape = par[1], scale = par[2], log = log)
    if(sum) {
      if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y))
      val <- sum(weights * val)
    }
    return(val)
  }
  
  
  sdist <- function(y, eta, weights = NULL, sum = FALSE) {   
    par <- c(exp(eta[1]), exp(eta[2]))                           
    score <- cbind((-log(par[2]) + log(y) - 1/lgamma(par[1]) * digamma(par[1])) * par[1], 
                   (-par[1]/par[2] + y/par[2]^2) * par[2])
    score <- as.matrix(score)
    colnames(score) <- etanames
    if(sum) {
      if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y))
      # if score == Inf replace score with 1.7e308 because Inf*0 would lead to NaN -> gradient is NaN
      score[score==Inf] = 1.7e308
      score <- colSums(weights * score)
    }
    return(score)
  }
  
  
  hdist <- function(y, eta, weights = NULL) {    
    ny <- length(y)
    if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, ny)
    
    par <- c(exp(eta[1]), exp(eta[2]))                          
    
    d2ld.etamu2 <- sum(weights * ((1/lgamma(par[1])^2 * digamma(par[1])^2 - 1/lgamma(par[1]) * trigamma(par[1])) * par[1]^2 + (-log(par[2]) + log(y) - 1/lgamma(par[1]) * digamma(par[1])) * par[1]))
    d2ld.etamu.d.etasigma <- sum(weights * (-par[1]))          # FIX ME: should be ... for exact parameters (here ~ e-17 due to calculations)
    d2ld.etasigma2 <- sum(weights * (par[1]/par[2]^2 - 2*y/par[2]^3) * par[2]^2 + (-par[1]/par[2] + y/par[2]^2) * par[2])         
    
    hess <- matrix(c(d2ld.etamu2, d2ld.etamu.d.etasigma, d2ld.etamu.d.etasigma, d2ld.etasigma2), nrow = 2)
    colnames(hess) <- rownames(hess) <-  etanames
    
    return(hess)
  }
  
  
  ## additional functions pdist, qdist, rdist
  pdist <- function(q, eta, lower.tail = TRUE, log.p = FALSE) pgamma(q, shape = exp(eta[1]), scale = exp(eta[2]), 
                                                                    lower.tail = lower.tail, log.p = log.p)
  qdist <- function(p, eta, lower.tail = TRUE, log.p = FALSE) qgamma(p, shape = exp(eta[1]), scale = exp(eta[2]),
                                                                    lower.tail = lower.tail, log.p = log.p)
  rdist <- function(n, eta) rgamma(n, shape = exp(eta[1]), scale = exp(eta[2]))

  
  
  link <- c("log", "log")
  
  linkfun <- function(par) {
    eta <- c(log(par[1]), log(par[2]))
    names(eta) <- etanames
    return(eta)
  }
  
  
  linkinv <- function(eta) {
    par <- c(exp(eta[1]), exp(eta[2]))
    names(par) <- parnames
    return(par)
  }
  
  
  linkinvdr <- function(eta) {
    dpardeta <- c(exp(eta[1]), exp(eta[2]))
    names(dpardeta) <- parnames
    return(dpardeta)
  }
  
  
  startfun <- function(y, weights = NULL){
    y.m <- if(is.null(weights) || (length(weights)==0L)) mean(y) else weighted.mean(y, weights)
    y.sd <- if(is.null(weights) || (length(weights)==0L)) sd(y) else sqrt(Hmisc::wtd.var(y, weights))
    shape <- (y.m/y.sd)^2
    scale <- y.m/shape     # <- y.sd^2/y.m
    starteta <- c(log(shape), log(scale))
    names(starteta) <- etanames
    return(starteta)
  }
  
  mle <- FALSE
  
  dist_list_gamma <- list(family.name = "Gamma Distribution",
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
                          mle = mle
  )
}




