###### Gaussian distribution
dist_gaussian <- function() {

  # parnames <- c("mu", "sigma")
  # etanames <- c("mu", "log(sigma)")
  
  
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
    colnames(score) <- c("mu", "log(sigma)")
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
    colnames(hess) <- rownames(hess) <-  c("mu", "log(sigma)")
    
    return(hess)
  }
  
  
  ## additional functions pdist, qdist, rdist on link scale
  pdist <- function(q, eta, lower.tail = TRUE, log.p = FALSE) pnorm(q, mean = eta[1], sd = exp(eta[2]), 
                                                                            lower.tail = lower.tail, log.p = log.p)
  qdist <- function(p, eta, lower.tail = TRUE, log.p = FALSE) qnorm(p, mean = eta[1], sd = exp(eta[2]), 
                                                                            lower.tail = lower.tail, log.p = log.p)
  rdist <- function(n, eta) rnorm(n, mean = eta[1], sd = exp(eta[2]))
  
  
  link <- c("identity", "log")
  
  linkfun <- function(par) {
    eta <- c(par[1], log(par[2]))
    names(eta) <- c("mu", "log(sigma)")
    return(eta)
  }
  
  
  linkinv <- function(eta) {
    par <- c(eta[1], exp(eta[2]))
    names(par) <- c("mu", "sigma")
    return(par)
  }
  
  
  linkinvdr <- function(eta) {
    dpardeta <- c(1, exp(eta[2]))
    names(dpardeta) <- c("mu", "sigma")
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
    names(starteta) <- c("mu", "log(sigma)")
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
       mle = mle,
       gamlssobj = FALSE,
       censored = FALSE
  )

}





###### censored gaussian or logistic distributions
dist_crch <- function(dist = c("gaussian","logistic"), 
                      truncated = FALSE,
                      type = c("left", "right", "interval"),
                      censpoint = 0)
{
  
  if(dist == "normal") dist <- "gaussian"
  if(length(dist) > 1) dist <- dist[1]
  if(length(type) > 1) type <- type[1]
  
  if((type == "interval") && !(length(censpoint) == 2)) stop("for type 'intervall' two censoring points have to be set")
  if((type %in% c("left", "right")) && !(length(censpoint) ==1)) stop("for type 'left' or 'right' one censoring point has to be set")
  
  if(truncated) family.name <- if(dist == "gaussian") "truncated Normal Distribution" else "truncated Logistic Distribution"
  if(!truncated) family.name <- if(dist == "gaussian") "censored Normal Distribution" else "censored Logistic Distribution"
  
  if(type == "interval") {
    if(truncated) list.name <- paste("dist_list_trunc", type, censpoint[1], censpoint[2], sep = "_")
    if(!truncated) list.name <- paste("dist_list_cens", type, censpoint[1], censpoint[2], sep = "_")
    left <- censpoint[1]
    right <- censpoint[2]
  } else {
    if(truncated) list.name <- paste("dist_list_trunc", type, censpoint, sep = "_")
    if(!truncated) list.name <- paste("dist_list_cens", type, censpoint, sep = "_")
    if(type == "left") {
      left <- censpoint
      right <- Inf
    }
    if(type == "right"){
      left <- -Inf
      right <- censpoint
    }
  }
  
  
  dist_list <- list()
  
  # parnames <- c("mu", "sigma")
  # etanames <- c("mu", "log(sigma)")
  
  
  
  
  

  if(truncated) {
    if(dist == "gaussian") {
      
      ddist <-  function(y, eta, log = TRUE, weights = NULL, sum = FALSE) {     
        par <- c(eta[1], exp(eta[2]))
        val <- crch::dtnorm(x = y, mean = par[1], sd = par[2], left = left, right = right, log = log)
        if(sum) {
          if(is.null(weights) || (length(weights)==0L)) weights <- if(is.matrix(y)) rep.int(1, dim(y)[1]) else rep.int(1, length(y))
          val <- sum(weights * val, na.rm = TRUE)
        }
        return(val)
      }
      
      
      sdist <- function(y, eta, weights = NULL, sum = FALSE) {   
        par <- c(eta[1], exp(eta[2]))
        # y[y==0] <- 1e-323
        
        score_m <- crch:::stnorm(x = y, mean = par[1], sd = par[2], which = "mu", left = left, right = right)
        score_s <- crch:::stnorm(x = y, mean = par[1], sd = par[2], which = "sigma", left = left, right = right) * exp(eta[2]) # inner derivation exp(eta[2])
        score <- cbind(score_m, score_s)
        score <- as.matrix(score)
        colnames(score) <- c("mu", "log(sigma)")
        if(sum) {
          if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y)[1])
          # if score == Inf replace score with 1.7e308 because Inf*0 would lead to NaN (0 in weights)
          score[score==Inf] = 1.7e308
          score <- colSums(weights * score, na.rm = TRUE)
          #if(any(is.nan(score))) print(c(eta, "y", y))
        }
        return(score)
      }
      
      
      hdist <- function(y, eta, weights = NULL) {    
        ny <- length(y)
        if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, ny)
        
        par <- c(eta[1], exp(eta[2]))                           
        # y[y==0] <- 1e-323
        
        d2mu <- crch:::htnorm(x = y, mean = par[1], sd = par[2], which = "mu", left = left, right = right)
        d2sigma <- crch:::htnorm(x = y, mean = par[1], sd = par[2], which = "sigma", left = left, right = right)
        dmudsigma <- crch:::htnorm(x = y, mean = par[1], sd = par[2], which = "mu.sigma", left = left, right = right) # FIX: order?
        dsigmadmu <- crch:::htnorm(x = y, mean = par[1], sd = par[2], which = "sigma.mu", left = left, right = right) # FIX: order?
        dsigma <- crch:::stnorm(x = y, mean = par[1], sd = par[2], which = "sigma", left = left, right = right)
        
        d2ld.etamu2 <- sum(weights * d2mu, na.rm = TRUE)
        d2ld.etamu.d.etasigma <- sum(weights * dmudsigma * par[2], na.rm = TRUE)
        d2ld.etasigma.d.etamu <- sum(weights * dsigmadmu * par[2], na.rm = TRUE)
        d2ld.etasigma2 <- sum(weights * (d2sigma * exp(2*eta[2]) + dsigma * par[2]), na.rm = TRUE)         
        
        hess <- matrix(c(d2ld.etamu2, d2ld.etamu.d.etasigma, d2ld.etasigma.d.etamu, d2ld.etasigma2), nrow = 2)
        colnames(hess) <- rownames(hess) <-  c("mu", "log(sigma)")
        
        return(hess)
      }
      
      ## additional functions pdist, qdist, rdist on link scale
      pdist <- function(q, eta, lower.tail = TRUE, log.p = FALSE) crch::ptnorm(q, mean = eta[1], sd = exp(eta[2]), 
                                                                               lower.tail = lower.tail, log.p = log.p, 
                                                                               left = left, right = right)
      qdist <- function(p, eta, lower.tail = TRUE, log.p = FALSE) crch::qtnorm(p, mean = eta[1], sd = exp(eta[2]), 
                                                                               lower.tail = lower.tail, log.p = log.p, 
                                                                               left = left, right = right)
      rdist <- function(n, eta) crch::rtnorm(n, mean = eta[1], sd = exp(eta[2]), left = left, right = right)
    }
    
    
    
    
    if(dist == "logistic") {
      
      ddist <-  function(y, eta, log = TRUE, weights = NULL, sum = FALSE) {     
        par <- c(eta[1], exp(eta[2]))
        val <- crch::dtlogis(x = y, location = par[1], scale = par[2], left = left, right = right, log = log)
        if(sum) {
          if(is.null(weights) || (length(weights)==0L)) weights <- if(is.matrix(y)) rep.int(1, dim(y)[1]) else rep.int(1, length(y))
          val <- sum(weights * val, na.rm = TRUE)
        }
        return(val)
      }
      
      
      sdist <- function(y, eta, weights = NULL, sum = FALSE) {   
        par <- c(eta[1], exp(eta[2]))
        # y[y==0] <- 1e-323
        
        score_m <- crch:::stlogis(x = y, location = par[1], scale = par[2], which = "mu", left = left, right = right)
        score_s <- crch:::stlogis(x = y, location = par[1], scale = par[2], which = "sigma", left = left, right = right) * exp(eta[2]) # inner derivation exp(eta[2])
        score <- cbind(score_m, score_s)
        score <- as.matrix(score)
        colnames(score) <- c("mu", "log(sigma)")
        if(sum) {
          if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y)[1])
          # if score == Inf replace score with 1.7e308 because Inf*0 would lead to NaN (0 in weights)
          score[score==Inf] = 1.7e308
          score <- colSums(weights * score, na.rm = TRUE)
          #if(any(is.nan(score))) print(c(eta, "y", y))
        }
        return(score)
      }
      
      
      hdist <- function(y, eta, weights = NULL) {    
        ny <- length(y)
        if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, ny)
        
        par <- c(eta[1], exp(eta[2]))                           
        # y[y==0] <- 1e-323
        
        d2mu <- crch:::htlogis(x = y, location = par[1], scale = par[2], which = "mu", left = left, right = right)
        d2sigma <- crch:::htlogis(x = y, location = par[1], scale = par[2], which = "sigma", left = left, right = right)
        dmudsigma <- crch:::htlogis(x = y, location = par[1], scale = par[2], which = "mu.sigma", left = left, right = right) # FIX: order?
        dsigmadmu <- crch:::htlogis(x = y, location = par[1], scale = par[2], which = "sigma.mu", left = left, right = right) # FIX: order?
        dsigma <- crch:::stlogis(x = y, location = par[1], scale = par[2], which = "sigma", left = left, right = right)
        
        d2ld.etamu2 <- sum(weights * d2mu, na.rm = TRUE)
        d2ld.etamu.d.etasigma <- sum(weights * dmudsigma * par[2], na.rm = TRUE)
        d2ld.etasigma.d.etamu <- sum(weights * dsigmadmu * par[2], na.rm = TRUE)
        d2ld.etasigma2 <- sum(weights * (d2sigma * exp(2*eta[2]) + dsigma * par[2]), na.rm = TRUE)         
        
        hess <- matrix(c(d2ld.etamu2, d2ld.etamu.d.etasigma, d2ld.etasigma.d.etamu, d2ld.etasigma2), nrow = 2)
        colnames(hess) <- rownames(hess) <-  c("mu", "log(sigma)")
        
        return(hess)
      }
      
      ## additional functions pdist, qdist, rdist on link scale
      # FIX ME: better par instead of eta?
      pdist <- function(q, eta, lower.tail = TRUE, log.p = FALSE) crch::ptlogis(q, location = eta[1], scale = exp(eta[2]), 
                                                                                lower.tail = lower.tail, log.p = log.p, 
                                                                                left = left, right = right)
      qdist <- function(p, eta, lower.tail = TRUE, log.p = FALSE) crch::qtlogis(p, location = eta[1], scale = exp(eta[2]), 
                                                                                lower.tail = lower.tail, log.p = log.p, 
                                                                                left = left, right = right)
      rdist <- function(n, eta) crch::rtlogis(n, location = eta[1], scale = exp(eta[2]), left = left, right = right)
    }
    
  }
  
  
  
  
  
  
  if(!truncated) {
    if(dist == "gaussian") {
      
      ddist <-  function(y, eta, log = TRUE, weights = NULL, sum = FALSE) {     
        par <- c(eta[1], exp(eta[2]))
        val <- crch::dcnorm(x = y, mean = par[1], sd = par[2], left = left, right = right, log = log)
        if(sum) {
          if(is.null(weights) || (length(weights)==0L)) weights <- if(is.matrix(y)) rep.int(1, dim(y)[1]) else rep.int(1, length(y))
          val <- sum(weights * val, na.rm = TRUE)
        }
        return(val)
      }
      
      
      sdist <- function(y, eta, weights = NULL, sum = FALSE) {   
        par <- c(eta[1], exp(eta[2]))
        # y[y==0] <- 1e-323
        
        score_m <- crch:::scnorm(x = y, mean = par[1], sd = par[2], which = "mu", left = left, right = right)
        score_s <- crch:::scnorm(x = y, mean = par[1], sd = par[2], which = "sigma", left = left, right = right) * exp(eta[2]) # inner derivation exp(eta[2])
        score <- cbind(score_m, score_s)
        score <- as.matrix(score)
        colnames(score) <- c("mu", "log(sigma)")
        if(sum) {
          if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y)[1])
          # if score == Inf replace score with 1.7e308 because Inf*0 would lead to NaN (0 in weights)
          score[score==Inf] = 1.7e308
          score <- colSums(weights * score, na.rm = TRUE)
          #if(any(is.nan(score))) print(c(eta, "y", y))
        }
        return(score)
      }
      
      
      hdist <- function(y, eta, weights = NULL) {    
        ny <- length(y)
        if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, ny)
        
        par <- c(eta[1], exp(eta[2]))                           
        # y[y==0] <- 1e-323
        
        d2mu <- crch:::hcnorm(x = y, mean = par[1], sd = par[2], which = "mu", left = left, right = right)
        d2sigma <- crch:::hcnorm(x = y, mean = par[1], sd = par[2], which = "sigma", left = left, right = right)
        dmudsigma <- crch:::hcnorm(x = y, mean = par[1], sd = par[2], which = "mu.sigma", left = left, right = right) # FIX: order?
        dsigmadmu <- crch:::hcnorm(x = y, mean = par[1], sd = par[2], which = "sigma.mu", left = left, right = right) # FIX: order?
        dsigma <- crch:::scnorm(x = y, mean = par[1], sd = par[2], which = "sigma", left = left, right = right)
        
        d2ld.etamu2 <- sum(weights * d2mu, na.rm = TRUE)
        d2ld.etamu.d.etasigma <- sum(weights * dmudsigma * par[2], na.rm = TRUE)
        d2ld.etasigma.d.etamu <- sum(weights * dsigmadmu * par[2], na.rm = TRUE)
        d2ld.etasigma2 <- sum(weights * (d2sigma * exp(2*eta[2]) + dsigma * par[2]), na.rm = TRUE)         
        
        hess <- matrix(c(d2ld.etamu2, d2ld.etamu.d.etasigma, d2ld.etasigma.d.etamu, d2ld.etasigma2), nrow = 2)
        colnames(hess) <- rownames(hess) <-  c("mu", "log(sigma)")
        
        return(hess)
      }
      
      ## additional functions pdist, qdist, rdist on link scale
      pdist <- function(q, eta, lower.tail = TRUE, log.p = FALSE) crch::pcnorm(q, mean = eta[1], sd = exp(eta[2]), 
                                                                               lower.tail = lower.tail, log.p = log.p, 
                                                                               left = left, right = right)
      qdist <- function(p, eta, lower.tail = TRUE, log.p = FALSE) crch::qcnorm(p, mean = eta[1], sd = exp(eta[2]), 
                                                                               lower.tail = lower.tail, log.p = log.p, 
                                                                               left = left, right = right)
      rdist <- function(n, eta) crch::rcnorm(n, mean = eta[1], sd = exp(eta[2]), left = left, right = right)
    }
    
    
    
  
    if(dist == "logistic") {
      
      ddist <-  function(y, eta, log = TRUE, weights = NULL, sum = FALSE) {     
        par <- c(eta[1], exp(eta[2]))
        val <- crch::dclogis(x = y, location = par[1], scale = par[2], left = left, right = right, log = log)
        if(sum) {
          if(is.null(weights) || (length(weights)==0L)) weights <- if(is.matrix(y)) rep.int(1, dim(y)[1]) else rep.int(1, length(y))
          val <- sum(weights * val, na.rm = TRUE)
        }
        return(val)
      }
      
      
      sdist <- function(y, eta, weights = NULL, sum = FALSE) {   
        par <- c(eta[1], exp(eta[2]))
        # y[y==0] <- 1e-323
        
        score_m <- crch:::sclogis(x = y, location = par[1], scale = par[2], which = "mu", left = left, right = right)
        score_s <- crch:::sclogis(x = y, location = par[1], scale = par[2], which = "sigma", left = left, right = right) * exp(eta[2]) # inner derivation exp(eta[2])
        score <- cbind(score_m, score_s)
        score <- as.matrix(score)
        colnames(score) <- c("mu", "log(sigma)")
        if(sum) {
          if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y)[1])
          # if score == Inf replace score with 1.7e308 because Inf*0 would lead to NaN (0 in weights)
          score[score==Inf] = 1.7e308
          score <- colSums(weights * score, na.rm = TRUE)
          #if(any(is.nan(score))) print(c(eta, "y", y))
        }
        return(score)
      }
      
      
      hdist <- function(y, eta, weights = NULL) {    
        ny <- length(y)
        if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, ny)
        
        par <- c(eta[1], exp(eta[2]))                           
        # y[y==0] <- 1e-323
        
        d2mu <- crch:::hclogis(x = y, location = par[1], scale = par[2], which = "mu", left = left, right = right)
        d2sigma <- crch:::hclogis(x = y, location = par[1], scale = par[2], which = "sigma", left = left, right = right)
        dmudsigma <- crch:::hclogis(x = y, location = par[1], scale = par[2], which = "mu.sigma", left = left, right = right) # FIX: order?
        dsigmadmu <- crch:::hclogis(x = y, location = par[1], scale = par[2], which = "sigma.mu", left = left, right = right) # FIX: order?
        dsigma <- crch:::sclogis(x = y, location = par[1], scale = par[2], which = "sigma", left = left, right = right)
        
        d2ld.etamu2 <- sum(weights * d2mu, na.rm = TRUE)
        d2ld.etamu.d.etasigma <- sum(weights * dmudsigma * par[2], na.rm = TRUE)
        d2ld.etasigma.d.etamu <- sum(weights * dsigmadmu * par[2], na.rm = TRUE)
        d2ld.etasigma2 <- sum(weights * (d2sigma * exp(2*eta[2]) + dsigma * par[2]), na.rm = TRUE)         
        
        hess <- matrix(c(d2ld.etamu2, d2ld.etamu.d.etasigma, d2ld.etasigma.d.etamu, d2ld.etasigma2), nrow = 2)
        colnames(hess) <- rownames(hess) <-  c("mu", "log(sigma)")
        
        return(hess)
      }
      
      ## additional functions pdist, qdist, rdist on link scale
      # FIX ME: better par instead of eta?
      pdist <- function(q, eta, lower.tail = TRUE, log.p = FALSE) crch::pclogis(q, location = eta[1], scale = exp(eta[2]), 
                                                                                lower.tail = lower.tail, log.p = log.p, 
                                                                                left = left, right = right)
      qdist <- function(p, eta, lower.tail = TRUE, log.p = FALSE) crch::qclogis(p, location = eta[1], scale = exp(eta[2]), 
                                                                                lower.tail = lower.tail, log.p = log.p, 
                                                                                left = left, right = right)
      rdist <- function(n, eta) crch::rclogis(n, location = eta[1], scale = exp(eta[2]), left = left, right = right)
    }
    
  }
    
    
  link <- c("identity", "log")
  
  linkfun <- function(par) {
    eta <- c(par[1], log(par[2]))
    names(eta) <- c("mu", "log(sigma)")
    return(eta)
  }
  
  
  linkinv <- function(eta) {
    par <- c(eta[1], exp(eta[2]))
    names(par) <- c("mu", "sigma")
    return(par)
  }
  
  
  linkinvdr <- function(eta) {
    dpardeta <- c(1, exp(eta[2]))
    names(dpardeta) <- c("mu", "sigma")
    return(dpardeta)
  }
  
  
  startfun <- function(y, weights = NULL){
    
    # alternative version using for mle = TRUE
    if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y))
    
    # observations with positive weights
    y_posweights <- y[weights>0]
    # if only one observation or only equal observations (with positive weights)
    # set mu to the value of this one observation 
    # and sigma to 1e-10
    if(length(unique(y_posweights)) == 1) {
      starteta <- c(unique(y_posweights), log(1e-10))
    } else {
      
      # if only one positive observation error in optim (non-finite values)
      if(sum(y_posweights>0) < 2){
        starteta <- c(0, log(1e-10))
      } else {
        
        x <- cbind(rep(1, length(y)))
        colnames(x) <- "(Intercept)"
        
        starteta <- try(unlist(crch::crch.fit(x = x, z = x, y = y, left = 0, right = Inf,
                                        truncated = truncated, weights = weights, dist = dist, 
                                        control = crch::crch.control(method = "L-BFGS-B", 
                                                               reltol = 1e-8, factr = 1e7,
                                                               maxit = 100,
                                                               hessian = FALSE))$coefficients), 
                        silent = TRUE)
        if(inherits(starteta, "try-error")){
          warning("Error for method L-BFGS-B in optim, applied method BFGS instead")
          starteta <- try(unlist(crch::crch.fit(x = x, z = x, y = y, left = 0, right = Inf,
                                          truncated = truncated, weights = weights, dist = dist, 
                                          control = crch::crch.control(method = "BFGS", 
                                                                 reltol = 1e-8, factr = 1e7,
                                                                 maxit = 100,
                                                                 hessian = FALSE))$coefficients), 
                          silent = TRUE)
        }
        if(inherits(starteta, "try-error")){
          warning("Error for method BFGS in optim, applied method Nelder-Mead instead")
          starteta <- try(unlist(crch::crch.fit(x = x, z = x, y = y, left = 0, right = Inf,
                                          truncated = truncated, weights = weights, dist = dist, 
                                          control = crch::crch.control(method = "Nelder-Mead", 
                                                                 reltol = 1e-8, factr = 1e7,
                                                                 #maxit = 100,
                                                                 hessian = FALSE))$coefficients), 
                          silent = TRUE)
        }
        if(inherits(starteta, "try-error")){
          warning("Error for method Nelder-Mead in optim, applied method SANN instead")
          starteta <- try(unlist(crch::crch.fit(x = x, z = x, y = y, left = 0, right = Inf,
                                          truncated = truncated, weights = weights, dist = dist, 
                                          control = crch::crch.control(method = "SANN", 
                                                                 reltol = 1e-8, factr = 1e7,
                                                                 #maxit = 100,
                                                                 hessian = FALSE))$coefficients), 
                          silent = TRUE)
        }
        if(inherits(starteta, "try-error")) {
          print("ERROR for all optimization methods")
          print(c("summary y:", summary(y)))
          print(c("summary weights:", summary(weights)))
          save(y, file = "~/svn/partykit/pkg/disttree/demo/error_examples/y.rda")
          save(weights, file = "~/svn/partykit/pkg/disttree/demo/error_examples/weights.rda")
        }
      }
    }    
    names(starteta) <- c("mu", "log(sigma)")
    return(starteta)
    
    
    ## for mle = FALSE (to use optim())
    ## FIX ME: in case the data is censored with other censoring points, optional ?
    #if(type == "interval"){
    #  yc[yc<censpoint[1]] <- censpoint[1]
    #  yc[yc>censpoint[2]] <- censpoint[2]
    #}
    #if(type == "left") yc <- pmax(censpoint,y) 
    #if(type == "right") yc <- pmin(censpoint,y) 
    
    #if(is.null(weights) || (length(weights)==0L)) {
    #  mu <- mean(yc)
    #  sigma <- sqrt(1/length(yc) * sum((yc - mu)^2))
    #} else {
    #  mu <- weighted.mean(yc, weights)
    #  sigma <- sqrt(1/sum(weights) * sum(weights * (yc - mu)^2))
    #}
    
    #starteta <- c(mu, log(sigma))
    #names(starteta) <- c("mu", "log(sigma)")
    #return(starteta)
  }
  
  #mle <- FALSE
  mle <- TRUE
  
  dist_list <- list(family.name = family.name,
                    list.name = list.name,
                    type = type,
                    dist = dist,
                    censpoint = censpoint,
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
                    censored = !truncated,
                    truncated = truncated
  )
  
  return(dist_list)
}





###### Weibull distribution
## FIX ME: adapt initial values for Weibull distribution
dist_weibull <- function() {
  
  # parnames <- c("mean", "scale")
  # etanames <- c("mean", "log(scale)")
  
  
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
    colnames(score) <- c("mean", "log(scale)")
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
    colnames(hess) <- rownames(hess) <-  c("mean", "log(scale)")
    
    return(hess)
  }
  
  
  
  
  ## additional functions pdist, qdist, rdist on link scale
  # FIX ME: par instead of eta better?
  pdist <- function(q, eta, lower.tail = TRUE, log.p = FALSE) pweibull(q, shape = 1/exp(eta[2]), scale = exp(eta[1]), lower.tail = lower.tail, log.p = log.p)
  qdist <- function(p, eta, lower.tail = TRUE, log.p = FALSE) qweibull(p, shape = 1/exp(eta[2]), scale = exp(eta[1]), lower.tail = lower.tail, log.p = log.p)
  rdist <- function(n, eta) rweibull(n, shape = 1/exp(eta[2]), scale = exp(eta[1]))
  
  
  
  link <- c("identity", "log")
  
  linkfun <- function(par) {
    eta <- c(par[1], log(par[2]))
    names(eta) <- c("mean", "log(scale)")
    return(eta)
  }
  
  
  linkinv <- function(eta) {
    par <- c(eta[1], exp(eta[2]))
    names(par) <- c("mean", "scale")
    return(par)
  }
  
  
  linkinvdr <- function(eta) {
    dpardeta <- c(1, exp(eta[2]))
    names(dpardeta) <- c("mean", "scale")
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
    names(starteta) <- c("mean", "log(scale)")
    return(starteta)
  }
  
  mle <- FALSE
  
  list(family.name = "Weibull Distribution",
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




###### Poisson distribution
dist_poisson <- function() {
  
  # parnames <- c("mu")
  # etanames <- c("log(mu)")
  
  
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
    colnames(score) <- c("log(mu)")
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
    colnames(hess) <- rownames(hess) <-  c("log(mu)")
    
    return(hess)
  }
  
  
  ## additional functions pdist, qdist, rdist on link scale
  # FIX ME: par instead of eta better?
  pdist <- function(q, eta, lower.tail = TRUE, log.p = FALSE) ppois(q, lambda = exp(eta), lower.tail = lower.tail, log.p = log.p)
  qdist <- function(p, eta, lower.tail = TRUE, log.p = FALSE) qpois(p, lambda = exp(eta), lower.tail = lower.tail, log.p = log.p)
  rdist <- function(n, eta) rpois(n, lambda = exp(eta))
  
  
  
  link <- c("log")
  
  linkfun <- function(par) {
    eta <- log(par)
    names(eta) <- c("log(mu)")
    return(eta)
  }
  
  
  linkinv <- function(eta) {
    par <- exp(eta)
    names(par) <- c("mu")
    return(par)
  }
  
  
  linkinvdr <- function(eta) {
    dpardeta <- exp(eta)
    names(dpardeta) <- c("mu")
    return(dpardeta)
  }
  
  
  startfun <- function(y, weights = NULL){
    mu <- if(is.null(weights) || (length(weights)==0L)) mean(y) else weighted.mean(y, weights)
    starteta <- log(mu)
    names(starteta) <- c("log(mu)")
    return(starteta)
  }
  
  
  mle <- TRUE
  
  
  list(family.name = "Poisson Distribution",
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





###### Exponential distribution
dist_exponential <- function() {
  
  # parnames <- c("lambda")
  # etanames <- c("log(lambda)")
  
  
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
    colnames(score) <- c("log(lambda)")
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
    colnames(hess) <- rownames(hess) <-  c("log(lambda)")
    
    return(hess)
  }
  
  
  ## additional functions pdist, qdist, rdist on link-scale
  # FIX ME: par instead of eta better?
  pdist <- function(q, eta, lower.tail = TRUE, log.p = FALSE) pexp(q, rate = exp(eta), lower.tail = lower.tail, log.p = log.p)
  qdist <- function(p, eta, lower.tail = TRUE, log.p = FALSE) qexp(p, rate = exp(eta), lower.tail = lower.tail, log.p = log.p)
  rdist <- function(n, eta) rexp(n, rate = exp(eta))
  
  
  link <- c("log")
  
  linkfun <- function(par) {
    eta <- log(par)
    names(eta) <- c("log(lambda)")
    return(eta)
  }
  
  
  linkinv <- function(eta) {
    par <- exp(eta)
    names(par) <- c("lambda")
    return(par)
  }
  
  
  linkinvdr <- function(eta) {
    dpardeta <- exp(eta)
    names(dpardeta) <- c("lambda")
    return(dpardeta)
  }
  
  
  startfun <- function(y, weights = NULL){
    lambda <- if(is.null(weights) || (length(weights)==0L)) length(y)/sum(y) else sum(weights)/sum(weights * y)
    starteta <- log(lambda)
    names(starteta) <- c("log(lambda)")
    return(starteta)
  }
  
  mle <- TRUE
  
  
  list(family.name = "Exponential Distribution",
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





###### Gamma distribution
dist_gamma <- function() {
  
  # parnames <- c("shape", "scale")
  # etanames <- c("log(shape)", "log(scale)")
  
  
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
    colnames(score) <- c("log(shape)", "log(scale)")
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
    colnames(hess) <- rownames(hess) <-  c("log(shape)", "log(scale)")
    
    return(hess)
  }
  
  
  ## additional functions pdist, qdist, rdist on link-scale
  # FIX ME: par instead of eta better?
  pdist <- function(q, eta, lower.tail = TRUE, log.p = FALSE) pgamma(q, shape = exp(eta[1]), scale = exp(eta[2]), 
                                                                     lower.tail = lower.tail, log.p = log.p)
  qdist <- function(p, eta, lower.tail = TRUE, log.p = FALSE) qgamma(p, shape = exp(eta[1]), scale = exp(eta[2]),
                                                                     lower.tail = lower.tail, log.p = log.p)
  rdist <- function(n, eta) rgamma(n, shape = exp(eta[1]), scale = exp(eta[2]))
  
  
  
  link <- c("log", "log")
  
  linkfun <- function(par) {
    eta <- c(log(par[1]), log(par[2]))
    names(eta) <- c("log(shape)", "log(scale)")
    return(eta)
  }
  
  
  linkinv <- function(eta) {
    par <- c(exp(eta[1]), exp(eta[2]))
    names(par) <- c("shape", "scale")
    return(par)
  }
  
  
  linkinvdr <- function(eta) {
    dpardeta <- c(exp(eta[1]), exp(eta[2]))
    names(dpardeta) <- c("shape", "scale")
    return(dpardeta)
  }
  
  
  startfun <- function(y, weights = NULL){
    y.m <- if(is.null(weights) || (length(weights)==0L)) mean(y) else weighted.mean(y, weights)
    y.sd <- if(is.null(weights) || (length(weights)==0L)) sd(y) else sqrt(Hmisc::wtd.var(y, weights))
    shape <- (y.m/y.sd)^2
    scale <- y.m/shape     # <- y.sd^2/y.m
    starteta <- c(log(shape), log(scale))
    names(starteta) <- c("log(shape)", "log(scale)")
    return(starteta)
  }
  
  mle <- FALSE
  
  list(family.name = "Gamma Distribution",
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





###### Zero truncated negative binomial
dist_ztnbinom <- function() {
  
  ## parnames <- c("mu", "sigma")
  ## etanames <- c("log(mu)", "log(sigma)")
  
  ddist <- function(y, eta, log = TRUE, weights = NULL, sum = FALSE) {
    
    rval <- countreg::dztnbinom(y, mu = exp(eta[1]), sigma = exp(eta[2]), log = log)
    
    if(sum) {
      if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y))
      rval <- sum(weights * rval, na.rm = TRUE)
    }
    return(rval)
  }
  
  sdist <- function(y, eta, weights = NULL, sum = FALSE) {
    
    par  <- c(exp(eta[1]), exp(eta[2]))
    rval <- countreg::sztnbinom(y, mu = par[1], sigma = par[2], drop = FALSE)
    rval <- cbind(rval[,1] * par[1], rval[,2] * par[2])
    colnames(rval) <- c("log(mu)", "log(sigma)")
    
    if(sum) {
      if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y))
      rval[rval==Inf] = 1.7e308  ## Taken from dist_gaussian
      rval <- colSums(weights * rval, na.rm = TRUE)
    }
    return(rval)
  }
  
  hdist <- function(y, eta, weights = NULL) {
    ny <- length(y)
    if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, ny)
    
    par <- c(exp(eta[1]), exp(eta[2]))
    s   <- countreg::sztnbinom(y, mu = weights * par[1], sigma = weights * par[2],
                               drop = FALSE)
    h   <- countreg::hztnbinom(y, mu = weights * par[1], sigma = weights * par[2],
                               drop = FALSE, parameter = c("mu", "sigma", "mu.sigma"))
    
    par <- as.matrix(par)
    
    hess <- -s * par + h[c("mu", "sigma"), ] * par^2
    hess <- cbind(hess, h["mu.sigma", ] * par[1, ] * par[2, ])
    hess <- colSums(hess)
    
    rval <- matrix(c(hess[1], rep(hess[3], 2), hess[2]), ncol = 2)
    colnames(rval) <- rownames(rval) <- c("log(mu)", "log(sigma)")
    
    return(rval)
  }
  
  pdist <- function(q, eta, lower.tail = TRUE, log.p = FALSE) {
    countreg::pztnbinom(q, mu = exp(eta[1]), sigma = exp(eta[2]),
                        lower.tail = lower.tail, log.p = log.p)
  }
  
  qdist <- function(p, eta, lower.tail = TRUE, log.p = FALSE) {
    countreg::qztnbinom(p, mu = exp(eta[1]), sigma = exp(eta[2]),
                        lower.tail = lower.tail, log.p = log.p)
  }
  
  rdist <- function(n, eta) {
    countreg::rztnbinom(n, mu = exp(eta[1]), sigma = exp(eta[2]))
  }
  
  link <- c("log", "log")
  
  linkfun <- function(par) {
    eta <- c(log(par[1]), log(par[2]))
    names(eta) <- c("log(mu)", "log(sigma)")
    return(eta)
  }
  
  linkinv <- function(eta) {
    par <- c(exp(eta[1]), exp(eta[2]))
    names(par) <- c("mu", "sigma")
    return(par)
  }
  
  linkinvdr <- function(eta) {
    dpardeta <- c(exp(eta[1]), exp(eta[2]))
    names(dpardeta) <- c("mu", "sigma")
    return(dpardeta)
  }
  
  startfun <- function(y, weights = NULL) {
    ## I guess there are better starting values than those.
    ## However, MoMs are hart du calculate analytically.
    if(is.null(weights) || (length(weights)==0L)) {
      mu    <- mean(y) - 1 + 0.00001
      sigma <- 1
    } else {
      mu    <- weighted.mean(y, weights) - 1 + 0.00001
      sigma <- 1
    }
    starteta <- c(log(mu), log(sigma))
    names(starteta) <- c("log(mu)", "log(sigma)")
    return(starteta)
  }
  
  list(family.name = "ztnbinom",
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
       mle = FALSE,
       gamlssobj = FALSE,
       censored = FALSE
  )
}





###### Binomial
dist_binomial <- function() {
  
  # parnames <- "mu"
  # etanames <- "logit(mu)"
  
  ddist <- function(y, eta, log = TRUE, weights = NULL, sum = FALSE) {
    
    par  <- 1 / (1 + exp(-eta))
    rval <- dbinom(y, size = 1, prob = par, log = log)
    
    if(sum) {
      if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y))
      rval <- sum(weights * rval, na.rm = TRUE)
    }
    return(rval)
  }
  
  sdist <- function(y, eta, weights = NULL, sum = FALSE) {
    
    par  <- 1 / (1 + exp(-eta))
    rval <- matrix(y - par, ncol = 1)
    colnames(rval) <- "logit(mu)"
    
    if(sum) {
      if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y))
      rval[rval==Inf] = 1.7e308
      rval <- colSums(weights * rval, na.rm = TRUE)
    }
    return(rval)
  } 
  
  hdist <- function(y, eta, weights = NULL) {
    
    par  <- 1 / (1 + exp(-eta))
    rval <- matrix(par * (1 - par), ncol = 1, nrow = 1,
                   dimnames = list(colnames = "logit(mu)",
                                   rownames = "logit(mu)"))
    
    return(rval)
  }
  
  pdist <- function(q, eta, lower.tail = TRUE, log.p = FALSE) {
    par <- 1 / (1 + exp(-eta))
    pbinom(q, size = 1, prob = par, lower.tail = lower.tail, log.p = log.p)
  }
  
  qdist <- function(p, eta, lower.tail = TRUE, log.p = FALSE) {
    par <- 1 / (1 + exp(-eta))
    qbinom(p, size = 1, prob = par, lower.tail = lower.tail, log.p = log.p)
  }
  
  rdist <- function(n, eta) {
    par <- 1 / (1 + exp(-eta))
    rbinom(n, size = 1, prob = par)
  }
  
  link <- "logit"
  
  linkfun <- function(par) {
    eta <- log(par) - log(1 - par)
    names(eta) <- "logit(mu)"
    return(eta)
  }
  
  linkinv <- function(eta) {
    par <- 1 / (1 + exp(-eta))
    names(par) <- "mu"
    return(par)
  }
  
  linkinvdr <- function(eta) {
    expeta   <- exp(-eta)
    dpardeta <- expeta / (1 + expeta)^2
    names(dpardeta) <- "mu"
    return(dpardeta)
  }
  
  startfun <- function(y, weights = NULL) {
    if(is.logical(y)) y <- as.numeric(y)
    if(is.factor(y)) {
      if(length(levels(y)) > 2) stop("factor variable can only have 2 levels for binomial distribution")
      y <- as.numeric(y)-1
    }
    if(is.null(weights) || (length(weights)==0L)) {
      par <- mean(y)
    } else {
      par <- weighted.mean(y, weights)
    }
    starteta <- log(par) - log(1 - par)
    names(starteta) <- "logit(mu)"
    return(starteta)
  }
  
  list(family.name = "Binomial",
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
       mle = TRUE,
       gamlssobj = FALSE,
       censored = FALSE
  )
}






##########################################
#complete distribution lists



#### dist_list_normal
{
  dist_list_normal <- list()
  
  # parnames <- c("mu", "sigma")
  # etanames <- c("mu", "log(sigma)")
  
  
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
    colnames(score) <- c("mu", "log(sigma)")
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
    colnames(hess) <- rownames(hess) <-  c("mu", "log(sigma)")
    
    return(hess)
  }
  
  
  
  ## additional functions pdist, qdist, rdist on link scale
  pdist <- function(q, eta, lower.tail = TRUE, log.p = FALSE) pnorm(q, mean = eta[1], sd = exp(eta[2]), 
                                                                    lower.tail = lower.tail, log.p = log.p)
  qdist <- function(p, eta, lower.tail = TRUE, log.p = FALSE) qnorm(p, mean = eta[1], sd = exp(eta[2]), 
                                                                    lower.tail = lower.tail, log.p = log.p)
  rdist <- function(n, eta) rnorm(n, mean = eta[1], sd = exp(eta[2]))
  
  
  link <- c("identity", "log")
  
  linkfun <- function(par) {
    eta <- c(par[1], log(par[2]))
    names(eta) <- c("mu", "log(sigma)")
    return(eta)
  }
  
  
  linkinv <- function(eta) {
    par <- c(eta[1], exp(eta[2]))
    names(par) <- c("mu", "sigma")
    return(par)
  }
  
  
  linkinvdr <- function(eta) {
    dpardeta <- c(1, exp(eta[2]))
    names(dpardeta) <- c("mu", "sigma")
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
    names(starteta) <- c("mu", "log(sigma)")
    return(starteta)
  }
  
  mle <- TRUE
  
  dist_list_normal <- list(family.name = "Normal Distribution",
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





#### dist_list_cens_normal
{
  dist_list_cens_normal <- list()
  
  # parnames <- c("mu", "sigma")
  # etanames <- c("mu", "log(sigma)")
  
  ddist <-  function(y, eta, log = TRUE, weights = NULL, sum = FALSE, left = 0, right = Inf) {     
    par <- c(eta[1], exp(eta[2]))
    val <- crch::dcnorm(x = y, mean = par[1], sd = par[2], left = left, right = right, log = log)
    if(sum) {
      if(is.null(weights) || (length(weights)==0L)) weights <- if(is.matrix(y)) rep.int(1, dim(y)[1]) else rep.int(1, length(y))
      val <- sum(weights * val, na.rm = TRUE)
    }
    return(val)
  }
  
  
  sdist <- function(y, eta, weights = NULL, sum = FALSE, left = 0, right = Inf) {   
    par <- c(eta[1], exp(eta[2]))
    # y[y==0] <- 1e-323
    
    score_m <- crch:::scnorm(x = y, mean = par[1], sd = par[2], which = "mu", left = left, right = right)
    score_s <- crch:::scnorm(x = y, mean = par[1], sd = par[2], which = "sigma", left = left, right = right) * exp(eta[2]) # inner derivation exp(eta[2])
    score <- cbind(score_m, score_s)
    score <- as.matrix(score)
    colnames(score) <- c("mu", "log(sigma)")
    if(sum) {
      if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y)[1])
      # if score == Inf replace score with 1.7e308 because Inf*0 would lead to NaN (0 in weights)
      score[score==Inf] = 1.7e308
      score <- colSums(weights * score, na.rm = TRUE)
      #if(any(is.nan(score))) print(c(eta, "y", y))
    }
    return(score)
  }
  
  
  hdist <- function(y, eta, weights = NULL, left = 0, right = Inf) {    
    ny <- length(y)
    if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, ny)
    
    par <- c(eta[1], exp(eta[2]))                           
    # y[y==0] <- 1e-323
    
    d2mu <- crch:::hcnorm(x = y, mean = par[1], sd = par[2], which = "mu", left = left, right = right)
    d2sigma <- crch:::hcnorm(x = y, mean = par[1], sd = par[2], which = "sigma", left = left, right = right)
    dmudsigma <- crch:::hcnorm(x = y, mean = par[1], sd = par[2], which = "mu.sigma", left = left, right = right) # FIX: order?
    dsigmadmu <- crch:::hcnorm(x = y, mean = par[1], sd = par[2], which = "sigma.mu", left = left, right = right) # FIX: order?
    dsigma <- crch:::scnorm(x = y, mean = par[1], sd = par[2], which = "sigma", left = left, right = right)
    
    d2ld.etamu2 <- sum(weights * d2mu, na.rm = TRUE)
    d2ld.etamu.d.etasigma <- sum(weights * dmudsigma * par[2], na.rm = TRUE)
    d2ld.etasigma.d.etamu <- sum(weights * dsigmadmu * par[2], na.rm = TRUE)
    d2ld.etasigma2 <- sum(weights * (d2sigma * exp(2*eta[2]) + dsigma * par[2]), na.rm = TRUE)         
    
    hess <- matrix(c(d2ld.etamu2, d2ld.etamu.d.etasigma, d2ld.etasigma.d.etamu, d2ld.etasigma2), nrow = 2)
    colnames(hess) <- rownames(hess) <-  c("mu", "log(sigma)")
    
    return(hess)
  }
  
  
  ## additional functions pdist, qdist, rdist on link-scale
  # FIX ME: par instead of eta better?
  pdist <- function(q, eta, lower.tail = TRUE, log.p = FALSE) crch::pcnorm(q, mean = eta[1], sd = exp(eta[2]), 
                                                                            lower.tail = lower.tail, log.p = log.p, 
                                                                            left = left, right = right)
  qdist <- function(p, eta, lower.tail = TRUE, log.p = FALSE) crch::qcnorm(p, mean = eta[1], sd = exp(eta[2]), 
                                                                            lower.tail = lower.tail, log.p = log.p, 
                                                                            left = left, right = right)
  rdist <- function(n, eta) crch::rcnorm(n, mean = eta[1], sd = exp(eta[2]), left = left, right = right)
  
  
  link <- c("identity", "log")
  
  linkfun <- function(par) {
    eta <- c(par[1], log(par[2]))
    names(eta) <- c("mu", "log(sigma)")
    return(eta)
  }
  
  
  linkinv <- function(eta) {
    par <- c(eta[1], exp(eta[2]))
    names(par) <- c("mu", "sigma")
    return(par)
  }
  
  
  linkinvdr <- function(eta) {
    dpardeta <- c(1, exp(eta[2]))
    names(dpardeta) <- c("mu", "sigma")
    return(dpardeta)
  }
  
  startfun <- function(y, weights = NULL){
    
    ## alternative version using for mle = TRUE
    if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y))
    
    # observations with positive weights
    y_posweights <- y[weights>0]
    # if only one observation or only equal observations (with positive weights)
    # set mu to the value of this one observation 
    # and sigma to 1e-10
    if(length(unique(y_posweights)) == 1) {
      starteta <- c(unique(y_posweights), log(1e-10))
    } else {
      
      # if only one positive observation error in optim (non-finite values)
      if(sum(y_posweights>0) < 2){
        starteta <- c(0, log(1e-10))
      } else {
        
        x <- cbind(rep(1, length(y)))
        colnames(x) <- "(Intercept)"
        
        starteta <- try(unlist(crch::crch.fit(x = x, z = x, y = y, left = 0, right = Inf,
                                        truncated = FALSE, weights = weights, 
                                        control = crch::crch.control(method = "L-BFGS-B", 
                                                               reltol = 1e-8, factr = 1e7,
                                                               maxit = 100,
                                                               hessian = FALSE))$coefficients), 
                        silent = TRUE)
        if(inherits(starteta, "try-error")){
          warning("Error for method L-BFGS-B in optim, applied method BFGS instead")
          starteta <- try(unlist(crch::crch.fit(x = x, z = x, y = y, left = 0, right = Inf,
                                          truncated = FALSE, weights = weights, 
                                          control = crch::crch.control(method = "BFGS", 
                                                                 reltol = 1e-8, factr = 1e7,
                                                                 maxit = 100,
                                                                 hessian = FALSE))$coefficients), 
                          silent = TRUE)
        }
        if(inherits(starteta, "try-error")){
          warning("Error for method BFGS in optim, applied method Nelder-Mead instead")
          starteta <- try(unlist(crch::crch.fit(x = x, z = x, y = y, left = 0, right = Inf,
                                          truncated = FALSE, weights = weights, 
                                          control = crch::crch.control(method = "Nelder-Mead", 
                                                                 reltol = 1e-8, factr = 1e7,
                                                                 #maxit = 100,
                                                                 hessian = FALSE))$coefficients), 
                          silent = TRUE)
        }
        if(inherits(starteta, "try-error")){
          warning("Error for method Nelder-Mead in optim, applied method SANN instead")
          starteta <- try(unlist(crch::crch.fit(x = x, z = x, y = y, left = 0, right = Inf,
                                          truncated = FALSE, weights = weights, 
                                          control = crch::crch.control(method = "SANN", 
                                                                 reltol = 1e-8, factr = 1e7,
                                                                 #maxit = 100,
                                                                 hessian = FALSE))$coefficients), 
                          silent = TRUE)
        }
        if(inherits(starteta, "try-error")) {
          print("ERROR for all optimization methods")
          print(c("summary y:", summary(y)))
          print(c("summary weights:", summary(weights)))
          save(y, file = "~/svn/partykit/pkg/disttree/demo/error_examples/y.rda")
          save(weights, file = "~/svn/partykit/pkg/disttree/demo/error_examples/weights.rda")
        }
      }
    }
    names(starteta) <- c("mu", "log(sigma)")
    return(starteta)
    
    ## for mle = FALSE (to use optim())
    #yc <- pmax(0,y)  # optional ?
    #if(is.null(weights)) {
    #  mu <- mean(yc)
    #  sigma <- sqrt(1/length(yc) * sum((yc - mu)^2))
    #} else {
    #  mu <- weighted.mean(yc, weights)
    #  sigma <- sqrt(1/sum(weights) * sum(weights * (yc - mu)^2))
    #}
    #starteta <- c(mu, log(sigma))
    #names(starteta) <- c("mu", "log(sigma)")
    #return(starteta)
  }
  
  
  mle <- TRUE
  #mle <- FALSE
  
  dist_list_cens_normal <- list(family.name = "censored Normal Distribution",
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
                                censored = TRUE
  )
}






#### dist_list_trunc_normal
{
  dist_list_trunc_normal <- list()
  
  # parnames <- c("mu", "sigma")
  # etanames <- c("mu", "log(sigma)")
  
  ddist <-  function(y, eta, log = TRUE, weights = NULL, sum = FALSE, left = 0, right = Inf) {     
    par <- c(eta[1], exp(eta[2]))
    val <- crch::dtnorm(x = y, mean = par[1], sd = par[2], left = left, right = right, log = log)
    if(sum) {
      if(is.null(weights) || (length(weights)==0L)) weights <- if(is.matrix(y)) rep.int(1, dim(y)[1]) else rep.int(1, length(y))
      val <- sum(weights * val, na.rm = TRUE)
    }
    return(val)
  }
  
  
  sdist <- function(y, eta, weights = NULL, sum = FALSE, left = 0, right = Inf) {   
    par <- c(eta[1], exp(eta[2]))
    # y[y==0] <- 1e-323
    
    score_m <- crch:::stnorm(x = y, mean = par[1], sd = par[2], which = "mu", left = left, right = right)
    score_s <- crch:::stnorm(x = y, mean = par[1], sd = par[2], which = "sigma", left = left, right = right) * exp(eta[2]) # inner derivation exp(eta[2])
    score <- cbind(score_m, score_s)
    score <- as.matrix(score)
    colnames(score) <- c("mu", "log(sigma)")
    if(sum) {
      if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y)[1])
      # if score == Inf replace score with 1.7e308 because Inf*0 would lead to NaN (0 in weights)
      score[score==Inf] = 1.7e308
      score <- colSums(weights * score, na.rm = TRUE)
      #if(any(is.nan(score))) print(c(eta, "y", y))
    }
    return(score)
  }
  
  
  hdist <- function(y, eta, weights = NULL, left = 0, right = Inf) {    
    ny <- length(y)
    if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, ny)
    
    par <- c(eta[1], exp(eta[2]))                           
    # y[y==0] <- 1e-323
    
    d2mu <- crch:::htnorm(x = y, mean = par[1], sd = par[2], which = "mu", left = left, right = right)
    d2sigma <- crch:::htnorm(x = y, mean = par[1], sd = par[2], which = "sigma", left = left, right = right)
    dmudsigma <- crch:::htnorm(x = y, mean = par[1], sd = par[2], which = "mu.sigma", left = left, right = right) # FIX: order?
    dsigmadmu <- crch:::htnorm(x = y, mean = par[1], sd = par[2], which = "sigma.mu", left = left, right = right) # FIX: order?
    dsigma <- crch:::stnorm(x = y, mean = par[1], sd = par[2], which = "sigma", left = left, right = right)
    
    d2ld.etamu2 <- sum(weights * d2mu, na.rm = TRUE)
    d2ld.etamu.d.etasigma <- sum(weights * dmudsigma * par[2], na.rm = TRUE)
    d2ld.etasigma.d.etamu <- sum(weights * dsigmadmu * par[2], na.rm = TRUE)
    d2ld.etasigma2 <- sum(weights * (d2sigma * exp(2*eta[2]) + dsigma * par[2]), na.rm = TRUE)         
    
    hess <- matrix(c(d2ld.etamu2, d2ld.etamu.d.etasigma, d2ld.etasigma.d.etamu, d2ld.etasigma2), nrow = 2)
    colnames(hess) <- rownames(hess) <-  c("mu", "log(sigma)")
    
    return(hess)
  }
  
  
  ## additional functions pdist, qdist, rdist on link-scale
  # FIX ME: par instead of eta better?
  pdist <- function(q, eta, lower.tail = TRUE, log.p = FALSE) crch::ptnorm(q, mean = eta[1], sd = exp(eta[2]), 
                                                                           lower.tail = lower.tail, log.p = log.p, 
                                                                           left = left, right = right)
  qdist <- function(p, eta, lower.tail = TRUE, log.p = FALSE) crch::qtnorm(p, mean = eta[1], sd = exp(eta[2]), 
                                                                           lower.tail = lower.tail, log.p = log.p, 
                                                                           left = left, right = right)
  rdist <- function(n, eta) crch::rtnorm(n, mean = eta[1], sd = exp(eta[2]), left = left, right = right)
  
  
  link <- c("identity", "log")
  
  linkfun <- function(par) {
    eta <- c(par[1], log(par[2]))
    names(eta) <- c("mu", "log(sigma)")
    return(eta)
  }
  
  
  linkinv <- function(eta) {
    par <- c(eta[1], exp(eta[2]))
    names(par) <- c("mu", "sigma")
    return(par)
  }
  
  
  linkinvdr <- function(eta) {
    dpardeta <- c(1, exp(eta[2]))
    names(dpardeta) <- c("mu", "sigma")
    return(dpardeta)
  }
  
  startfun <- function(y, weights = NULL){
    
    if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y))
    ypos <- y[y > 0]
    y_posweights <- y[y>0 & weights>0]
    
    # if only one or none positive observation or only equal positive observations
    # set mu to the value of this one observation or 0 for no positive observations
    # and sigma to 1e-10
    if(length(unique(y_posweights)) < 2){
      if(length(unique(y_posweights)) == 1) {
        starteta <- c(unique(y_posweights), log(1e-10))
      } else {
        starteta <- c(0, log(1e-10))}
    } else {
      
      xpos <- cbind(rep(1, length(ypos)))
      colnames(xpos) <- "(Intercept)"
      wpos <- weights[y>0]
      
      # calculate starteta using crch::crch.fit
      starteta <- try(unlist(crch::crch.fit(x = xpos, z = xpos, y = ypos, left = 0, right = Inf,
                                      truncated = TRUE, weights = wpos, 
                                      control = crch::crch.control(method = "L-BFGS-B", 
                                                             reltol = 1e-8, factr = 1e7,
                                                             maxit = 100,
                                                             hessian = FALSE))$coefficients), 
                      silent = TRUE)
      if(inherits(starteta, "try-error")){
        warning("Error for method L-BFGS-B in optim, applied method BFGS instead")
        starteta <- try(unlist(crch::crch.fit(x = xpos, z = xpos, y = ypos, left = 0, right = Inf,
                                        truncated = TRUE, weights = wpos, 
                                        control = crch::crch.control(method = "BFGS", 
                                                               reltol = 1e-8, factr = 1e7,
                                                               maxit = 100,
                                                               hessian = FALSE))$coefficients), 
                        silent = TRUE)
      }
      if(inherits(starteta, "try-error")){
        warning("Error for method BFGS in optim, applied method Nelder-Mead instead")
        starteta <- try(unlist(crch::crch.fit(x = xpos, z = xpos, y = ypos, left = 0, right = Inf,
                                        truncated = TRUE, weights = wpos, 
                                        control = crch::crch.control(method = "Nelder-Mead", 
                                                               reltol = 1e-8, factr = 1e7,
                                                               #maxit = 100,
                                                               hessian = FALSE))$coefficients), 
                        silent = TRUE)
      }
      if(inherits(starteta, "try-error")){
        warning("Error for method Nelder-Mead in optim, applied method SANN instead")
        starteta <- try(unlist(crch::crch.fit(x = xpos, z = xpos, y = ypos, left = 0, right = Inf,
                                        truncated = TRUE, weights = wpos, 
                                        control = crch::crch.control(method = "SANN", 
                                                               reltol = 1e-8, factr = 1e7,
                                                               #maxit = 100,
                                                               hessian = FALSE))$coefficients), 
                        silent = TRUE)
      }
      if(inherits(starteta, "try-error")) {
        print("ERROR for all optimization methods")
        print(c("summary y:", summary(y)))
        print(c("summary ypos:", summary(ypos)))
        print(c("summary weights:", summary(weights)))
        print(c("summary wpos:", summary(wpos)))
        save(y, file = "~/svn/partykit/pkg/disttree/demo/error_examples/y.rda")
        save(weights, file = "~/svn/partykit/pkg/disttree/demo/error_examples/weights.rda")
      }
    }
    names(starteta) <- c("mu", "log(sigma)")
    return(starteta)
  }
  
  
  mle <- TRUE
  
  dist_list_trunc_normal <- list(family.name = "truncated Normal Distribution",
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
                                 censored = FALSE   ## FIX ME: new argument 'truncated'?
  )
}






#### dist_list_hurdle_normal
{
  dist_list_hurdle_normal <- list()
  
  # parnames <- c("mu", "sigma", "nu")
  # etanames <- c("mu", "log(sigma)", "log(nu/(1-nu))")
  
  ddist <-  function(y, eta, log = TRUE, weights = NULL, sum = FALSE, left = 0, right = Inf) {     
    par <- c(eta[1], exp(eta[2]), exp(eta[3])/(1 + exp(eta[3])))
    
    if(log) {
      val <- (y>0) * (crch::dtnorm(x = y, mean = par[1], sd = par[2], left = left, right = right, log = TRUE) + log(par[3])) +
        (y==0) * log(1-par[3])
    } else {
      val <- (y>0) * (crch::dtnorm(x = y, mean = par[1], sd = par[2], left = left, right = right, log = FALSE) * par[3]) +
        (y==0) * (1-par[3])
    }
    if(sum) {
      if(is.null(weights) || (length(weights)==0L)) weights <- if(is.matrix(y)) rep.int(1, dim(y)[1]) else rep.int(1, length(y))
      val <- sum(weights * val, na.rm = TRUE)
    }
    return(val)
  }
  
  
  sdist <- function(y, eta, weights = NULL, sum = FALSE, left = 0, right = Inf) {   
    par <- c(eta[1], exp(eta[2]), exp(eta[3])/(1 + exp(eta[3])))
    # y[y==0] <- 1e-323
    
    ## score / estfun (first-order partial derivatives of the (positive) log-likelihood function by eta)
    score_m <- crch:::stnorm(x = y, mean = par[1], sd = par[2], which = "mu", left = left, right = right)
    score_s <- crch:::stnorm(x = y, mean = par[1], sd = par[2], which = "sigma", left = left, right = right) * exp(eta[2])       # inner derivation exp(eta[2])
    score_p <- (y>0) * 1/(1+exp(eta[3])) + (y==0) * (-exp(eta[3]/(1+exp(eta[3]))))
    score <- cbind(score_m, score_s, score_p)
    score <- as.matrix(score)
    colnames(score) <- c("mu", "log(sigma)", "log(nu/(1-nu))")
    if(sum) {
      if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y)[1])
      # if score == Inf replace score with 1.7e308 because Inf*0 would lead to NaN (0 in weights)
      score[score==Inf] = 1.7e308
      score <- colSums(weights * score, na.rm = TRUE)
      #if(any(is.nan(score))) print(c(eta, "y", y))
    }
    return(score)
  }
  
  
  hdist <- function(y, eta, weights = NULL, left = 0, right = Inf) {    
    ny <- length(y)
    if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, ny)
    
    par <- c(eta[1], exp(eta[2]), exp(eta[3])/(1 + exp(eta[3])))                           
    # y[y==0] <- 1e-323
    
    d2mu <- crch:::htnorm(x = y, mean = par[1], sd = par[2], which = "mu", left = left, right = right)
    d2sigma <- crch:::htnorm(x = y, mean = par[1], sd = par[2], which = "sigma", left = left, right = right)
    d2nu <- (y>0) * (-1/(par[3]^2)) + (y==0) * (-1/((1-par[3])^2))
    dmudsigma <- crch:::htnorm(x = y, mean = par[1], sd = par[2], which = "mu.sigma", left = left, right = right) # FIX: order?
    dmudnu <- 0
    dsigmadmu <- crch:::htnorm(x = y, mean = par[1], sd = par[2], which = "sigma.mu", left = left, right = right) # FIX: order?
    dsigmadnu <- 0
    dnudmu <- 0
    dnudsigma <- 0
    dsigma <- crch:::stnorm(x = y, mean = par[1], sd = par[2], which = "sigma", left = left, right = right)
    
    d2ld.etamu2 <- sum(weights * d2mu, na.rm = TRUE)
    d2ld.etamu.d.etasigma <- sum(weights * dmudsigma * par[2], na.rm = TRUE)
    d2ld.etasigma.d.etamu <- sum(weights * dsigmadmu * par[2], na.rm = TRUE)
    d2ld.etasigma2 <- sum(weights * (d2sigma * exp(2*eta[2]) + dsigma * par[2]), na.rm = TRUE)      
    d2ld.etamu.d.etanu <- 0
    d2ld.etasigma.d.etanu <- 0 
    d2ld.etanu.d.etamu <- 0
    d2ld.etanu.d.etasigma <- 0 
    d2ld.etanu2 <- sum(weights * (-exp(eta[3])/((1+exp(eta[3]))^2)), na.rm = TRUE)
    
    hess <- matrix(c(d2ld.etamu2, d2ld.etamu.d.etasigma, d2ld.etamu.d.etanu,
                     d2ld.etasigma.d.etamu, d2ld.etasigma2, d2ld.etasigma.d.etanu,
                     d2ld.etanu.d.etamu, d2ld.etanu.d.etasigma, d2ld.etanu2), 
                   nrow = 3)
    colnames(hess) <- rownames(hess) <-  c("mu", "log(sigma)", "log(nu/(1-nu))")
    
    return(hess)
  }
  
  
  ## additional functions pdist, qdist, rdist on link-scale
  # FIX ME: par instead of eta better?
  # FIX ME: adaption to distribution with third parameter correct?
  #pdist <- function(q, eta, lower.tail = TRUE, log.p = FALSE) crch::ptnorm(q, mean = eta[1], sd = exp(eta[2]), 
  #                                                                         lower.tail = lower.tail, log.p = log.p, 
  #                                                                         left = left, right = right) * pbinom(q, 1, prob = exp(eta[3])/(1+exp(eta[3])))
  #qdist <- function(p, eta, lower.tail = TRUE, log.p = FALSE) crch::qtnorm(p, mean = eta[1], sd = exp(eta[2]), 
  #                                                                         lower.tail = lower.tail, log.p = log.p, 
  #                                                                         left = left, right = right) * qbinom(p, 1, prob = exp(eta[3])/(1+exp(eta[3])))
  #rdist <- function(n, eta) crch::rtnorm(n, mean = eta[1], sd = exp(eta[2]), left = left, right = right) * rbinom(n, 1, prob = exp(eta[3])/(1+exp(eta[3])))
  
  
  link <- c("identity", "log", "logit")
  
  linkfun <- function(par) {
    eta <- c(par[1], log(par[2]), log(par[3]/(1-par[3])))
    names(eta) <- c("mu", "log(sigma)", "log(nu/(1-nu))")
    return(eta)
  }
  
  
  linkinv <- function(eta) {
    par <- c(eta[1], exp(eta[2]), exp(eta[3])/(1 + exp(eta[3]))) 
    names(par) <- c("mu", "sigma", "nu")
    return(par)
  }
  
  
  linkinvdr <- function(eta) {
    dpardeta <- c(1, exp(eta[2]), exp(eta[3])/(1 + exp(eta[3]))^2)
    names(dpardeta) <- c("mu", "sigma", "nu")
    return(dpardeta)
  }
  
  
  startfun <- function(y, weights = NULL){
    
    if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y))
    ypos <- y[y > 0]
    y_posweights <- y[y>0 & weights>0]
    
    # if only one or none positive observation or only equal positive observations
    # set mu to the value of this one observation or 0 for no positive observations
    # and sigma to 1e-10
    if(length(unique(y_posweights)) < 2){
      if(length(unique(y_posweights)) == 1) {
        starteta <- c(unique(y_posweights), log(1e-10), qlogis(weighted.mean(y > 0, w = weights)))
      } else {
        starteta <- c(0, log(1e-10), qlogis(weighted.mean(y > 0, w = weights)))}
    } else {
      
      xpos <- cbind(rep(1, length(ypos)))
      colnames(xpos) <- "(Intercept)"
      wpos <- weights[y>0]

      # calculate starteta using crch::crch.fit
      starteta <- try(c(unlist(crch::crch.fit(x = xpos, z = xpos, y = ypos, left = 0, right = Inf,
                                        truncated = TRUE, weights = wpos, 
                                        control = crch::crch.control(method = "L-BFGS-B", 
                                                               reltol = 1e-8, factr = 1e7,
                                                               maxit = 100,
                                                               hessian = FALSE))$coefficients), 
                        if(any(y[weights>0]==0)) qlogis(weighted.mean(y > 0, w = weights)) else qlogis(1-1e-16)), 
                      silent = TRUE)
      if(inherits(starteta, "try-error")){
        warning("Error for method L-BFGS-B in optim, applied method BFGS instead")
        starteta <- try(c(unlist(crch::crch.fit(x = xpos, z = xpos, y = ypos, left = 0, right = Inf,
                                          truncated = TRUE, weights = wpos, 
                                          control = crch::crch.control(method = "BFGS", 
                                                                 reltol = 1e-8, factr = 1e7,
                                                                 maxit = 100,
                                                                 hessian = FALSE))$coefficients), 
                          if(any(y[weights>0]==0)) qlogis(weighted.mean(y > 0, w = weights)) else qlogis(1-1e-16)), 
                        silent = TRUE)
      }
      if(inherits(starteta, "try-error")){
        warning("Error for method BFGS in optim, applied method Nelder-Mead instead")
        starteta <- try(c(unlist(crch::crch.fit(x = xpos, z = xpos, y = ypos, left = 0, right = Inf,
                                          truncated = TRUE, weights = wpos, 
                                          control = crch::crch.control(method = "Nelder-Mead", 
                                                                 reltol = 1e-8, factr = 1e7,
                                                                 #maxit = 100,
                                                                 hessian = FALSE))$coefficients), 
                          if(any(y[weights>0]==0)) qlogis(weighted.mean(y > 0, w = weights)) else qlogis(1-1e-16)),
                        silent = TRUE)
      }
      if(inherits(starteta, "try-error")) {
        warning("Error for method Nelder-Mead in optim, applied method SANN instead")
        starteta <- try(c(unlist(crch::crch.fit(x = xpos, z = xpos, y = ypos, left = 0, right = Inf,
                                          truncated = TRUE, weights = wpos, 
                                          control = crch::crch.control(method = "SANN", 
                                                                 reltol = 1e-8, factr = 1e7,
                                                                 #maxit = 100,
                                                                 hessian = FALSE))$coefficients), 
                          if(any(y[weights>0]==0)) qlogis(weighted.mean(y > 0, w = weights)) else qlogis(1-1e-16)),
                        silent = TRUE)
      }
      if(inherits(starteta, "try-error")) {
        print("ERROR for all optimization methods")
        print(c("summary y:", summary(y)))
        print(c("summary ypos:", summary(ypos)))
        print(c("summary weights:", summary(weights)))
        print(c("summary wpos:", summary(wpos)))
        save(y, file = "~/svn/partykit/pkg/disttree/demo/error_examples/y.rda")
        save(weights, file = "~/svn/partykit/pkg/disttree/demo/error_examples/weights.rda")
      }
    }
    names(starteta) <- c("mu", "log(sigma)", "log(nu/(1-nu))")
    return(starteta)
  }
  
  mle <- TRUE
    
  dist_list_hurdle_normal <- list(family.name = "hurdle Normal Distribution",
                                  ddist = ddist, 
                                  sdist = sdist, 
                                  hdist = hdist,
                                  #pdist = pdist,
                                  #qdist = qdist,
                                  #rdist = rdist,
                                  link = link, 
                                  linkfun = linkfun, 
                                  linkinv = linkinv, 
                                  linkinvdr = linkinvdr,
                                  startfun = startfun,
                                  mle = mle,
                                  gamlssobj = FALSE,
                                  censored = FALSE   ## FIX ME: new argument 'truncated'?
  )
}









