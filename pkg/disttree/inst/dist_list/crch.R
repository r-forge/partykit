###### function that returns dist_list for censored distributions
dist_crch <- function(dist = c("normal","logistic"), 
                                     type = c("left", "right", "interval"),
                                     censpoint = 0)
{
  if((type == "interval") && !(length(censpoint) == 2)) stop("for type 'intervall' two censoring points have to be set")
  if((type %in% c("left", "right")) && !(length(censpoint) ==1)) stop("for type 'left' or 'right' one censoring point has to be set")
  
  family.name <- if(dist == "normal") "censored Normal Distribution" else "censored Logistic Distribution"
  
  if(type == "interval") {
    list.name <- paste("dist_list_cens", type, censpoint[1], censpoint[2], sep = "_")
    left <- censpoint[1]
    right <- censpoint[2]
  } else {
    list.name <- paste("dist_list_cens", type, censpoint, sep = "_")
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
  
  parnames <- c("mu", "sigma")
  etanames <- c("mu", "log(sigma)")
  
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
    colnames(score) <- etanames
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
    colnames(hess) <- rownames(hess) <-  etanames
    
    return(hess)
  }
  
  
  
  if(dist == "normal") {
    
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
      colnames(score) <- etanames
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
      colnames(hess) <- rownames(hess) <-  etanames
      
      return(hess)
    }
    
    ## additional functions pdist, qdist, rdist
    pdist <- function(q, eta, lower.tail = TRUE, log.p = FALSE) crch:::pcnorm(q, mean = eta[1], sd = eta[2], 
                                                                              lower.tail = lower.tail, log.p = log.p, 
                                                                              left = left, right = right)
    qdist <- function(p, eta, lower.tail = TRUE, log.p = FALSE) crch:::qcnorm(p, mean = eta[1], sd = eta[2], 
                                                                              lower.tail = lower.tail, log.p = log.p, 
                                                                              left = left, right = right)
    rdist <- function(n, eta) crch:::rcnorm(n, mean = eta[1], sd = eta[2], left = left, right = right)
  }
  
  
  
  
  
  
  
  if(dist == "logistic") {
    
    ddist <-  function(y, eta, log = TRUE, weights = NULL, sum = FALSE) {     
      par <- c(eta[1], exp(eta[2]))
      val <- crch::dclogis(x = y, mean = par[1], sd = par[2], left = left, right = right, log = log)
      if(sum) {
        if(is.null(weights) || (length(weights)==0L)) weights <- if(is.matrix(y)) rep.int(1, dim(y)[1]) else rep.int(1, length(y))
        val <- sum(weights * val, na.rm = TRUE)
      }
      return(val)
    }
    
    
    sdist <- function(y, eta, weights = NULL, sum = FALSE) {   
      par <- c(eta[1], exp(eta[2]))
      # y[y==0] <- 1e-323
      
      score_m <- crch:::sclogis(x = y, mean = par[1], sd = par[2], which = "mu", left = left, right = right)
      score_s <- crch:::sclogis(x = y, mean = par[1], sd = par[2], which = "sigma", left = left, right = right) * exp(eta[2]) # inner derivation exp(eta[2])
      score <- cbind(score_m, score_s)
      score <- as.matrix(score)
      colnames(score) <- etanames
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
      
      d2mu <- crch:::hclogis(x = y, mean = par[1], sd = par[2], which = "mu", left = left, right = right)
      d2sigma <- crch:::hclogis(x = y, mean = par[1], sd = par[2], which = "sigma", left = left, right = right)
      dmudsigma <- crch:::hclogis(x = y, mean = par[1], sd = par[2], which = "mu.sigma", left = left, right = right) # FIX: order?
      dsigmadmu <- crch:::hclogis(x = y, mean = par[1], sd = par[2], which = "sigma.mu", left = left, right = right) # FIX: order?
      dsigma <- crch:::sclogis(x = y, mean = par[1], sd = par[2], which = "sigma", left = left, right = right)
      
      d2ld.etamu2 <- sum(weights * d2mu, na.rm = TRUE)
      d2ld.etamu.d.etasigma <- sum(weights * dmudsigma * par[2], na.rm = TRUE)
      d2ld.etasigma.d.etamu <- sum(weights * dsigmadmu * par[2], na.rm = TRUE)
      d2ld.etasigma2 <- sum(weights * (d2sigma * exp(2*eta[2]) + dsigma * par[2]), na.rm = TRUE)         
      
      hess <- matrix(c(d2ld.etamu2, d2ld.etamu.d.etasigma, d2ld.etasigma.d.etamu, d2ld.etasigma2), nrow = 2)
      colnames(hess) <- rownames(hess) <-  etanames
      
      return(hess)
    }
    
    ## additional functions pdist, qdist, rdist
    pdist <- function(q, eta, lower.tail = TRUE, log.p = FALSE) crch:::pclogis(q, mean = eta[1], sd = eta[2], 
                                                                              lower.tail = lower.tail, log.p = log.p, 
                                                                              left = left, right = right)
    qdist <- function(p, eta, lower.tail = TRUE, log.p = FALSE) crch:::qclogis(p, mean = eta[1], sd = eta[2], 
                                                                              lower.tail = lower.tail, log.p = log.p, 
                                                                              left = left, right = right)
    rdist <- function(n, eta) crch:::rclogis(n, mean = eta[1], sd = eta[2], left = left, right = right)
  }

  
  
  
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
    # FIX ME: in case the data is censored with other censoring points, optional ?
    if(type == "interval"){
      yc[yc<censpoint[1]] <- censpoint[1]
      yc[yc>censpoint[2]] <- censpoint[2]
    }
    if(type == "left") yc <- pmax(censpoint,y) 
    if(type == "right") yc <- pmin(censpoint,y) 
    
    if(is.null(weights) || (length(weights)==0L)) {
      mu <- mean(yc)
      sigma <- sqrt(1/length(yc) * sum((yc - mu)^2))
    } else {
      mu <- weighted.mean(yc, weights)
      sigma <- sqrt(1/sum(weights) * sum(weights * (yc - mu)^2))
    }
    starteta <- c(mu, log(sigma))
    names(starteta) <- etanames
    return(starteta)
  }
  
  mle <- FALSE
  
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
                    mle = mle
  )
  
  return(dist_list)
}
