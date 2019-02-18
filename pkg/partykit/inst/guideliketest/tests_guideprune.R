library("partykit")
library("Formula")
library("parallel")
# source("~/svn/partykit/pkg/partykit/inst/guideliketest/guidelike.R")
# source("~/partykit/inst/guideliketest/guidelike.R")
#source("guidelike.R")
# source("~/svn/partykit/pkg/partykit/inst/guideliketest/ccprune.R")
# source("~/partykit/inst/guideliketest/ccprune.R")
#source("ccprune.R")

# function to compute adjusted Rand Index
adj_rand_index <- function(x, y) {
  
  tab <- table(x, y)
  a <- rowSums(tab)
  b <- colSums(tab)
  
  M <- sum(choose(tab, 2))
  N <- choose(length(x), 2)
  A <- sum(choose(a, 2))
  B <- sum(choose(b, 2))
  
  c(ARI = (M - (A * B) / N) / (0.5 * (A + B) - (A * B) / N))
}



###########
# DGP

## FIX ME: 
# for vary_beta = "all":
# beta1 includes the effect of 2*delta, because otherwise (1*delta)
# the values for z1>0 would all be smaller than for z1>0 and the other way round
# -> similar scenario as for vary_beta = "beta0"
# -> guide1 could easily find split

dgp_stump <- function(nobs = 100, delta = 1, xi = 0.3,  
                      # xi can also be a vector with 2 elements for nrsteps = 2
                      sigma = 1, seed = 7, only_intercept = FALSE, 
                      binary_regressor = FALSE, binary_beta = TRUE,
                      vary_beta = c("all", "beta0", "beta1"),
                      beta0 = NULL, beta1 = NULL,
                      nrsteps = 1, nrlevels = 2, shift = 0,
                      z1dist = c("unif", "norm"),
                      nrsplitvar = 1 ){ # nrsplitvar not needed but to match form of dgp_tree
  
  if(nrsteps == 1 & nrlevels > 2) stop("for 1 step only two levels are possible")
  if(nrsteps == 1 & length(xi)>1) stop("for 1 step only one split point xi can be handed over")
  if(shift>0 & nrsteps != 2) stop("shift is only defined for 2 steps (can be used to shift xi1 and xi2 in order to avoid symmetry)")
  if(shift>0 & nrlevels != 2) stop("shift can only be applied for nrlevels = 2 (to avoid symmetry)")
  if(length(vary_beta) > 1) vary_beta <- vary_beta[1]
  if(!vary_beta %in% c("all", "beta0", "beta1")) stop("vary_beta has to be one of the following options: 'all', 'beta0', 'beta1'")
  if(vary_beta == "beta0" & is.null(beta1)) stop("if only beta0 varies, a fixed value has to set for beta1")
  if(vary_beta == "beta1" & is.null(beta0)) stop("if only beta1 varies, a fixed value has to set for beta0")
  if(length(z1dist) > 1) z1dist <- z1dist[1]
  if(!z1dist %in% c("norm", "unif")) stop("z1dist has to be one of the following options: 'norm', 'unif'")
  if(!nrsteps %in% c(1,2)) stop("nrsteps can only be 1 or 2")
  if(length(xi)>1 & xi[1]>xi[2]) stop("for more than 1 split point (xi) the points have to be sorted (increasing)")
  
  
  set.seed(seed)
  
  # beta0 and/or beta1 symmetric around 0? (otherwise two split points are handed over)
  sym <- length(xi)==1 & nrsteps>1
  
  if(z1dist == "norm") {
    z1 <- rnorm(nobs, 0, 1) 
    z2 <- runif(nobs,-1,1)
  } else {
    z1 <- runif(nobs,-1,1)
    z2 <- rnorm(nobs, 0, 1) 
  }
  z3 <- rnorm(nobs, 0, 1)
  z4 <- rnorm(nobs, 0, 1)
  z5 <- rnorm(nobs, 0, 1)
  z6 <- rnorm(nobs, 0, 1)
  z7 <- runif(nobs, -1, 1)
  z8 <- runif(nobs, -1, 1)
  z9 <- runif(nobs, -1, 1)
  z10 <- runif(nobs, -1, 1)
  
  id <- numeric(length(z1))
  
  x <- if(binary_regressor) (-1)^rbinom(nobs, 1, 0.5) else runif(nobs, min = -1, max = 1)
  
  # for binary beta: one step at break point xi, 
  # for continuous beta: linear function
  if(nrsteps == 1) {
    if(vary_beta == "all"){
      if(binary_beta){
        beta0 <- delta * (-1)^(z1<xi)
        beta1 <- delta * (-1)^(z1<xi) * (-1)   # opposite signs  for beta0 and beta1   
      } else {
        beta0 <- delta * z1
        beta1 <- delta * z1 * (-1)   # opposite signs for beta0 and beta1
      }
    }
    
    if(vary_beta == "beta0"){
      beta0 <- if(binary_beta) delta * (-1)^(z1<xi) else delta * z1
      beta1 <- beta1
    }
    
    if(vary_beta == "beta1"){
      beta0 <- beta0
      beta1 <- if(binary_beta) delta * (-1)^(z1<xi) * (-1) else delta * z1 * (-1) 
      # opposite signs for beta0 and beta1
    }
    
    if(binary_beta) id <- 1+(z1>=xi)
  }  
  
  
  # for binary beta: two steps (one upwards and one downwards) at +/- xi, 
  # for continuous beta: quadratic function
  if(nrsteps == 2) {
    
    # set split points
    if(sym) {
      xi1 <- -xi 
      xi2 <- xi 
    } else {
      xi1 <- xi[1]
      xi2 <- xi[2]
    }
    if(shift>0) {
      xi1 <- xi1 + shift
      xi2 <- xi2 + shift
    }
    
    if(nrlevels == 2){
      if(vary_beta == "all"){
        if(binary_beta){
          beta0 <- delta * (-1)^(xi1<z1 & z1<xi2)
          beta1 <- delta * (-1)^(xi1<z1 & z1<xi2) * (-1)   # opposite signs for beta0 and beta1
        } else {
          if(shift>0){
            beta0 <- delta * ((z1-shift)^2 * 1/(1+shift)^2 * (z1<shift) + (z1-shift)^2 * 1/(1-shift)^2 * (z1>=shift)) * 2 -1
            beta1 <- (delta * ((z1-shift)^2 * 1/(1+shift)^2 * (z1<shift) + (z1-shift)^2 * 1/(1-shift)^2 * (z1>=shift)) * 2 -1) * (-1)   # opposite signs if both betas vary
          } else {
            beta0 <- delta * 2 * z1^2 - 1
            beta1 <- (delta * 2 * z1^2 - 1) * (-1)   # opposite signs for beta0 and beta1
          }
        }
      }
      
      if(vary_beta == "beta0"){
        beta0 <- if(binary_beta) delta * (-1)^(xi1<z1 & z1<xi2) else {
          if(shift>0){
            delta * ((z1-shift)^2 * 1/(1+shift)^2 * (z1<shift) + (z1-shift)^2 * 1/(1-shift)^2 * (z1>=shift)) * 2 -1
          } else {
            delta * 2 * z1^2 - 1
          }
        }
        beta1 <- beta1
      }
      
      if(vary_beta == "beta1"){
        beta0 <- beta0
        beta1 <- if(binary_beta) delta * (-1)^(xi1<z1 & z1<xi2) else {
          if(shift>0){
            delta * ((z1-shift)^2 * 1/(1+shift)^2 * (z1<shift) + (z1-shift)^2 * 1/(1-shift)^2 * (z1>=shift)) * 2 -1
          } else {
            delta * 2 * z1^2 - 1
          }
        }
      }
      
      if(binary_beta) id <- 1 + (z1 >= xi1) + (z1 >= xi2)  ## 2 levels but 3 subgroups in a tree (attention for stump = TRUE)
    }
    
    if(nrlevels == 3){
      if(vary_beta == "all"){
        if(binary_beta){
          beta0 <- delta * (-1)^(xi1<z1) * 0^(z1>=xi2)
          beta1 <- delta * (-1)^(xi1<z1) * 0^(z1>=xi2) * (-1)   # opposite signs for beta0 and beta1
        } else {
          beta0 <- (delta * z1^2 * 0.5^(z1>0)) * 2 - 1
          beta1 <- ((delta * z1^2 * 0.5^(z1>0)) * 2 - 1) * (-1)   # opposite signs for beta0 and beta1
        }
      }
      
      if(vary_beta == "beta0"){
        beta0 <- if(binary_beta) delta * (-1)^(xi1<z1) * 0^(z1>=xi2) else (delta * z1^2 * 0.5^(z1>0)) * 2 - 1
        beta1 <- beta1
      }
      
      if(vary_beta == "beta1"){
        beta0 <- beta0
        beta1 <- if(binary_beta) delta * (-1)^(xi1<z1) * 0^(z1>=xi2) else (delta * z1^2 * 0.5^(z1>0)) * 2 - 1
      }
      
      if(binary_beta) id <- 1 + (z1 >= xi1) + (z1 >= xi2)
    }
  }  
  
  
  mu <- if(only_intercept) beta0 else beta0 + beta1 * x
  
  y <- rnorm(nobs, mu, sigma)
  
  d <- data.frame(y = y, x = x, 
                  z1 = z1, z2 = z2, z3 = z3, z4 = z4, z5 = z5, z6 = z6, z7 = z7, z8 = z8, z9 = z9, z10 = z10,
                  beta0 = beta0, beta1 = beta1, mu = mu, sigma = rep.int(sigma, times = length(y)), id = id)
  
  return(d)
}



dgp_tree <- function(nobs = 100, delta = 1, xi = c(-0.3, 0.3),  
                     sigma = 1, seed = 7, only_intercept = FALSE, 
                     binary_regressor = FALSE, binary_beta = TRUE,
                     vary_beta = "all",
                     beta0 = NULL, beta1 = NULL,
                     nrsplitvar = 2,
                     z1dist = c("unif", "norm"),
                     nrsteps = 1, nrlevels = 2, shift = NULL # not needed but to match form of dgp_stump
)
{
  
  # check input values
  if(vary_beta != "all") stop("vary_beta can only be set to 'all' in dgp_tree")
  if(!binary_beta) stop("coefficients can only be binary in dgp_tree")
  if(!is.null(beta0) | !is.null(beta1)) warning("values for beta0 or beta1 are ignored since vary_beta='all' for dgp_tree")
  if(length(z1dist) > 1) z1dist <- z1dist[1]
  if(!z1dist %in% c("norm", "unif")) stop("z1dist has to be one of the following options: 'norm', 'unif'")
  if(nrsplitvar > 2) stop("not more than two splitvariables can be selected")
  #if(nrsplitvar == 2 & length(xi) < 2) stop("for two splitvariables two splitpoints need to be set in argument 'xi'")
  if(nrsplitvar == 1 & length(xi) == 2) if(xi[2]<xi[1]) stop("for one splitvariable with two splitpoints the first splitpoint
                                                             has to be smaller than or equal to the second")
  
  set.seed(seed)
  
  # beta0 and/or beta1 symmetric around 0? (otherwise two split points are handed over)
  if(nrsplitvar == 1 & nrsteps == 2){
    sym <- (length(xi)==1)
    if(sym) {
      xi1 <- -xi
      xi2 <- xi
    } else {
      xi1 <- xi[1]
      xi2 <- xi[2]
    }
  }  
  
  if(nrsplitvar == 2){
    if(length(xi)==1){
      xi1 <- xi2 <- xi
    } else {
      xi1 <- xi[1]
      xi2 <- xi[2]
    }
  }
  
  if(z1dist == "norm") {
    z1 <- rnorm(nobs, 0, 1) 
    z2 <- rnorm(nobs, 0, 1)
    z3 <- runif(nobs,-1,1)
    z4 <- runif(nobs,-1,1)
  } else {
    z1 <- runif(nobs,-1,1)
    z2 <- runif(nobs,-1,1) 
    z3 <- rnorm(nobs, 0, 1)
    z4 <- rnorm(nobs, 0, 1)
  }
  z5 <- rnorm(nobs, 0, 1)
  z6 <- rnorm(nobs, 0, 1)
  z7 <- rnorm(nobs, 0, 1)
  z8 <- runif(nobs, -1, 1)
  z9 <- runif(nobs, -1, 1)
  z10 <- runif(nobs, -1, 1)
  
  id <- numeric(length(z1))
  
  x <- if(binary_regressor) (-1)^rbinom(nobs, 1, 0.5) else runif(nobs, min = -1, max = 1)
  
  
  if(nrsplitvar == 1) {

    beta0 <- delta * (-1)^(z1<xi1)
    beta1 <- delta * (-1)^(z1<xi2) * (-1)   # opposite signs for beta0 and beta1 

    id <- 1+(z1>=xi1)+(z1>=xi2) 
  }
  
  
  
  if(nrsplitvar == 2) {
    
    ## first a split in z2 at split point xi2 => beta1 changes
    ## then a split in z1 at split point xi1 => beta0 changes
    # resulting subgroups:
    # group1 (z2<xi2):           beta0 = 0,      beta1 = delta
    # group2 (z2<xi2 & z1<xi1):  beta0 = -delta, beta1 = -delta
    # group3 (z2<xi2 & z1=>xi1): beta0 = +delta, beta1 = -delta
    
    #if(z2<xi2){
    #  beta0 <- 0 
    #  beta1 <- +delta
    #} else {
    #  if(z1<xi1){
    #    beta0 <- -delta 
    #    beta1 <- -delta
    #  } else {
    #    beta0 <- +delta 
    #    beta1 <- +delta
    #  }
    #}
    
    beta0 <- delta * (-1)^(z1<xi1) * 0^(z2<xi2)
    #beta0 <- delta * (-1)^(z1<xi1) * (z2>=xi2) + delta * (z2<xi2)
    beta1 <- delta * (-1)^(z2>=xi2)
    
    id <- 1 + (z2>=xi2) + (z2>=xi2)*(z1>=xi1)
  }
  
  
  mu <- if(only_intercept) beta0 else beta0 + beta1 * x
  
  y <- rnorm(nobs, mu, sigma)
  
  d <- data.frame(y = y, x = x, 
                  z1 = z1, z2 = z2, z3 = z3, z4 = z4, z5 = z5, z6 = z6, z7 = z7, z8 = z8, z9 = z9, z10 = z10,
                  beta0 = beta0, beta1 = beta1, mu = mu, sigma = rep.int(sigma, times = length(y)), id = id)
  
  return(d)
}

## illustration of data set
if(FALSE){
  d <- dgp_tree(nobs = 400, delta = 1, xi = c(-0.1, 0.2),  
                sigma = 1, seed = 7, only_intercept = FALSE, 
                binary_regressor = FALSE, binary_beta = TRUE,
                vary_beta = "all",
                beta0 = NULL, beta1 = NULL,
                nrsplitvar = 2,
                z1dist = "unif")
  
  plot(d[d$z2>=0.2 & d$z1 >= -0.1,]$x, d[d$z2>=0.2 & d$z1 >= -0.1,]$y, col = "red")
  points(d[d$z2>=0.2 & d$z1< -0.1,]$x, d[d$z2>=0.2 & d$z1< -0.1,]$y, col = "blue")
  points(d[d$z2<0.2,]$x, d[d$z2<0.2,]$y, col = "green")
  l <-  lm(d[d$z2<0.2,]$y ~ d[d$z2<0.2,]$x)
  abline(l)
  l2 <- lm(d[d$z2>=0.2 & d$z1< -0.1,]$y ~ d[d$z2>=0.2 & d$z1< -0.1,]$x)
  abline(l2)
  l3 <- lm(d[d$z2>=0.2 & d$z1 >= -0.1,]$y ~ d[d$z2>=0.2 & d$z1 >= -0.1,]$x)
  abline(l3)
}



# more complex version of dgp_tree allowing for more variations
if(FALSE){
  dgp_tree <- function(nobs = 100, delta = 1, xi = c(-0.3, 0.3),  
                       sigma = 1, seed = 7, only_intercept = FALSE, 
                       binary_regressor = TRUE, binary_beta = TRUE,
                       vary_beta = c("all", "beta0", "beta1"),
                       beta0 = NULL, beta1 = NULL,
                       nrsplitvar = 1,
                       z1dist = c("unif", "norm"),
                       nrsteps = NULL, nrlevels = NULL, shift = NULL # not needed but to match form of dgp_stump
  )
  {
    
    # check input values
    if(length(vary_beta) > 1) vary_beta <- vary_beta[1]
    if(!vary_beta %in% c("all", "beta0", "beta1")) stop("vary_beta has to be one of the following options: 'all', 'beta0', 'beta1'")
    if(vary_beta == "beta0" & is.null(beta1)) stop("if only beta0 varies, a fixed value has to set for beta1")
    if(vary_beta == "beta1" & is.null(beta0)) stop("if only beta1 varies, a fixed value has to set for beta0")
    if(length(z1dist) > 1) z1dist <- z1dist[1]
    if(!z1dist %in% c("norm", "unif")) stop("z1dist has to be one of the following options: 'norm', 'unif'")
    if(nrsplitvar > 2) stop("not more than two splitvariables can be selected")
    #if(nrsplitvar == 2 & length(xi) < 2) stop("for two splitvariables two splitpoints need to be set in argument 'xi'")
    if(nrsplitvar == 1 & length(xi) == 2) if(xi[2]<xi[1]) stop("for one splitvariable with two splitpoints the first splitpoint
                                                             has to be smaller than or equal to the second")
    
    
    
    set.seed(seed)
    
    # beta0 and/or beta1 symmetric around 0? (otherwise two split points are handed over)
    # (only for one splitvariable, because for two splitvariables xi has to have two elements)
    sym <- (nrsplitvar == 1 & length(xi)==1 & nrsteps>1)
    if(sym) {
      xi1 <- -xi
      xi2 <- xi
    } else {
      xi1 <- xi[1]
      xi2 <- xi[2]
    }
    
    equalxi <- (nrsplitvar == 2 & length(xi)==1)
    if(equalxi) xi1 <- xi2 <- xi
    
    if(z1dist == "norm") {
      z1 <- rnorm(nobs, 0, 1) 
      z2 <- rnorm(nobs, 0, 1)
      z3 <- runif(nobs,-1,1)
      z4 <- runif(nobs,-1,1)
    } else {
      z1 <- runif(nobs,-1,1)
      z2 <- runif(nobs,-1,1) 
      z3 <- rnorm(nobs, 0, 1)
      z4 <- rnorm(nobs, 0, 1)
    }
    z5 <- rnorm(nobs, 0, 1)
    z6 <- rnorm(nobs, 0, 1)
    z7 <- rnorm(nobs, 0, 1)
    z8 <- runif(nobs, -1, 1)
    z9 <- runif(nobs, -1, 1)
    z10 <- runif(nobs, -1, 1)
    
    id <- numeric(length(z1))
    
    x <- if(binary_regressor) (-1)^rbinom(nobs, 1, 0.5) else runif(nobs, min = -1, max = 1)
    
    
    if(nrsplitvar == 1) {
      # for binary beta: one step in each parameter at break point xi1 or xi2, 
      # for continuous beta: plogis function
      if(vary_beta == "all"){
        if(binary_beta){
          beta0 <- delta * (-1)^(z1<xi1)
          beta1 <- delta * (-1)^(z1<xi2) * (-1)   # opposite signs for beta0 and beta1 
        } else {
          beta0 <- delta * plogis(z1, xi1, 1/5) * 2 - 1
          beta1 <- delta * (plogis(z1, xi2, 1/5) * 2 - 1) * (-1)   # opposite signs for beta0 and beta1
        }
      }
      
      
      if(vary_beta == "beta0"){
        beta0 <- if(binary_beta) {
          delta * (-1)^(z1<xi1) * 0^(z1>=xi2) 
        } else {
          delta * ((1-plogis(z1, xi1, abs(xi2-xi1)/8))^(z1<(xi1+xi2)/2) * (1/2 * plogis(z1, xi2, abs(xi2-xi1)/8))^(z1>(xi1+xi2)/2) * 2 - 1)
        }
        beta1 <- beta1
      }
      
      if(vary_beta == "beta1"){
        beta0 <- beta0
        beta1 <- if(binary_beta) {
          delta * (-1)^(z1<xi1) * 0^(z1>=xi2) 
        } else {
          delta * ((1-plogis(z1, xi1, abs(xi2-xi1)/8))^(z1<(xi1+xi2)/2) * (1/2 * plogis(z1, xi2, abs(xi2-xi1)/8))^(z1>(xi1+xi2)/2) * 2 - 1)
        }
      }
      
      if(binary_beta) id <- 1+(z1>=xi1)+(z1>=xi2) 
    }
    
    
    
    if(nrsplitvar == 2) {
      
      # for binary beta: one step in each parameter at break point xi1 or xi2, 
      # for continuous beta: plogis function
      if(vary_beta == "all"){
        if(binary_beta){
          #if(z2<xi2){
          #  beta0 <- -delta 
          #  beta1 <- -delta
          #} else {
          #  if(z1<xi1){
          #    beta0 <- -delta 
          #    beta1 <- delta
          #  } else {
          #    beta0 <- delta 
          #    beta1 <- delta
          #  }
          #}
          
          beta0 <- (z2<xi2) * (-delta) + (z2>=xi2) * ((-delta) * (z1<xi1) + delta * (z1>=xi1))
          beta1 <- (z2<xi2) * (-delta) + (z2>=xi2) * delta
          
        } else {
          beta0 <- delta * ((plogis(z2, xi2, 1/5) * (plogis(z1, xi1, 1/5))) * 2 - 1)
          beta1 <- delta * (plogis(z2, xi2, 1/5) * 2 - 1)
        }
      }
      
      
      if(vary_beta == "beta0"){
        beta0 <- if(binary_beta) {
          (z2<xi2) * (-delta) + (z2>=xi2) * (delta) * (z1<xi1)
        } else {
          ## TO DO: check for 2 splitvariables
          delta * (plogis(z2, xi2, 1/5) * (2-plogis(z1, xi1, 1/5)) - 1)
        }
        beta1 <- beta1
      }
      
      
      if(vary_beta == "beta1"){
        beta0 <- beta0
        beta1 <- if(binary_beta) {
          (z2<xi2) * (-delta) + (z2>=xi2) * (delta) * (z1<xi1) 
        } else {
          ## TO DO: check for 2 splitvariables
          delta * (plogis(z2, xi2, 1/5) * (2-plogis(z1, xi1, 1/5)) - 1)
        }
      }
      
      if(binary_beta) id <- 1 + (z2>=xi2) + (z2>=xi2)*(z1>=xi1)
    }
    
    
    mu <- if(only_intercept) beta0 else beta0 + beta1 * x
    
    y <- rnorm(nobs, mu, sigma)
    
    d <- data.frame(y = y, x = x, 
                    z1 = z1, z2 = z2, z3 = z3, z4 = z4, z5 = z5, z6 = z6, z7 = z7, z8 = z8, z9 = z9, z10 = z10,
                    beta0 = beta0, beta1 = beta1, mu = mu, sigma = rep.int(sigma, times = length(y)), id = id)
    
    return(d)
  }
}


if(FALSE){
  xi1 <- -0.4
  xi2 <- 0.3
  d <- dgp_tree(400, sigma = 0.2, binary_regressor = FALSE, vary_beta = "all", xi = c(xi1, xi2))
  plot(y = d$y[d$z1 < xi2 & d$z1 > xi1], x = d$x[d$z1 < xi2 & d$z1 > xi1], ylim = c(-2,2))
  points(y = d$y[d$z1 >= xi2], x = d$x[d$z1 >= xi2], col = "blue")
  points(y = d$y[d$z1 < xi1], x = d$x[d$z1 < xi1], col = "red")
  l <- lm(y~x, data = d)
  abline(l)
  
  dt <- partykit::lmtree(y~x|z1 + z2 + z3, data = d)
  plot(dt)
}

## FIX ME: no decorrelation for length(guide_parm) = 1
evaltests <- function(formula, data, 
                      testfun = c("guide", "ctree", "mfluc", "ctree_cat", "mfluc_cat", 
                                  "ctree_bin", "mfluc_bin", "ctree_cat_bin", "mfluc_cat_bin"), 
                      subset, weights, stump = FALSE, alpha = 0.05,
                      whichscores = c("all", "residuals"),
                      ctree_max = c(FALSE, TRUE),
                      guide_testtype = c("sum", "max", "coin"), 
                      decorrelate = "vcov",
                      guide_parm = NULL,
                      xgroups = NULL, ygroups = NULL){
  
  na.action = na.pass
  converged = NULL
  scores = NULL
  doFit = TRUE
  offset = NULL 
  cluster = NULL
  if(testfun == "guide" & length(whichscores) == 1) stop("for GUIDE use argument 'guide_parm' to decide which scores should be used")
  if(length(whichscores) > 1) whichscores <- whichscores[1]
  if(length(ctree_max) > 1) ctree_max <- ctree_max[1]
  splittest <- ctree_max
  
  if(length(testfun) > 1) testfun <- testfun[1]
  if(!testfun %in% c("guide", "ctree", "mfluc", "ctree_cat", "mfluc_cat", 
                     "ctree_bin", "mfluc_bin", "ctree_cat_bin", "mfluc_cat_bin")) 
    stop("testfun hast to be one of the following options: 'guide', 'ctree', 'mfluc', 'ctree_cat', 'mfluc_cat',
         'ctree_bin', 'mfluc_bin', 'ctree_cat_bin', 'mfluc_cat_bin'")
  
  
  if(testfun == "guide"){
    
    if(length(guide_testtype) > 1) guide_testtype <- guide_testtype[1]
    if(!guide_testtype %in% c("max", "sum", "coin")) stop("testfun hast to be one of the following options: 'max', 'sum', 'coin'")
    
    control <- partykit:::extree_control(criterion = "p.value",
                                         selectfun = .guide_select(guide_decorrelate = "none",    # if required, decorrelation already within ytrafo 
                                                                   guide_testtype = guide_testtype,
                                                                   guide_parm = guide_parm,
                                                                   xgroups = xgroups,
                                                                   ygroups = ygroups,
                                                                   interaction = FALSE),  
                                         splitfun = partykit:::.objfun_split(restart = TRUE,
                                                                             intersplit = TRUE),
                                         svselectfun = .guide_select(guide_decorrelate = "none",    # if required, decorrelation already within ytrafo 
                                                                     guide_testtype = guide_testtype,
                                                                     guide_parm = guide_parm,
                                                                     xgroups = xgroups,
                                                                     ygroups = ygroups,
                                                                     interaction = FALSE), 
                                         svsplitfun = partykit:::.objfun_split(restart = TRUE,
                                                                               intersplit = TRUE),
                                         bonferroni = TRUE,
                                         logmincriterion = log(1-alpha),
                                         update = TRUE,
                                         stump = stump)
  }
  
  if(testfun == "ctree"){
    control <- partykit:::extree_control(criterion = "p.value",
                                         selectfun = partykit:::.ctree_select(teststat = "quadratic", splitstat = "quadratic", 
                                                                              splittest = splittest, pargs = GenzBretz(),
                                                                              testtype = "MonteCarlo", nresample = 9999L, 
                                                                              tol = sqrt(.Machine$double.eps),
                                                                              intersplit = TRUE, MIA = FALSE),
                                         splitfun = partykit:::.objfun_split(restart = TRUE,
                                                                             intersplit = TRUE),
                                         svselectfun = partykit:::.ctree_select(teststat = "quadratic", splitstat = "quadratic", 
                                                                                splittest = splittest, pargs = GenzBretz(),
                                                                                testtype = "MonteCarlo", nresample = 9999L, 
                                                                                tol = sqrt(.Machine$double.eps),
                                                                                intersplit = TRUE, MIA = FALSE), 
                                         svsplitfun = partykit:::.objfun_split(restart = TRUE,
                                                                               intersplit = TRUE),
                                         logmincriterion = log(1-alpha),
                                         update = TRUE,
                                         bonferroni = TRUE,
                                         stump = stump)
  }
  
  if(testfun == "ctree_cat"){
    control <- partykit:::extree_control(criterion = "p.value",
                                         selectfun = .ctree_cat_select(teststat = "quadratic", splitstat = "quadratic", 
                                                                       splittest = splittest, pargs = GenzBretz(),
                                                                       testtype = "MonteCarlo", nresample = 9999L, 
                                                                       tol = sqrt(.Machine$double.eps),
                                                                       intersplit = TRUE, MIA = FALSE,
                                                                       xgroups = xgroups),
                                         splitfun = partykit:::.objfun_split(restart = TRUE,
                                                                             intersplit = TRUE),
                                         svselectfun = .ctree_cat_select(teststat = "quadratic", splitstat = "quadratic", 
                                                                         splittest = splittest, pargs = GenzBretz(),
                                                                         testtype = "MonteCarlo", nresample = 9999L, 
                                                                         tol = sqrt(.Machine$double.eps),
                                                                         intersplit = TRUE, MIA = FALSE,
                                                                         xgroups = xgroups), 
                                         svsplitfun = partykit:::.objfun_split(restart = TRUE,
                                                                               intersplit = TRUE),
                                         logmincriterion = log(1-alpha),
                                         update = TRUE,
                                         bonferroni = TRUE,
                                         stump = stump)
  }
  
  if(testfun == "ctree_bin"){
    control <- partykit:::extree_control(criterion = "p.value",
                                         selectfun = .ctree_bin_select(teststat = "quadratic", splitstat = "quadratic", 
                                                                       splittest = splittest, pargs = GenzBretz(),
                                                                       testtype = "MonteCarlo", nresample = 9999L, 
                                                                       tol = sqrt(.Machine$double.eps),
                                                                       intersplit = TRUE, MIA = FALSE,
                                                                       xgroups = xgroups),
                                         splitfun = partykit:::.objfun_split(restart = TRUE,
                                                                             intersplit = TRUE),
                                         svselectfun = .ctree_bin_select(teststat = "quadratic", splitstat = "quadratic", 
                                                                         splittest = splittest, pargs = GenzBretz(),
                                                                         testtype = "MonteCarlo", nresample = 9999L, 
                                                                         tol = sqrt(.Machine$double.eps),
                                                                         intersplit = TRUE, MIA = FALSE,
                                                                         xgroups = xgroups), 
                                         svsplitfun = partykit:::.objfun_split(restart = TRUE,
                                                                               intersplit = TRUE),
                                         logmincriterion = log(1-alpha),
                                         update = TRUE,
                                         bonferroni = TRUE,
                                         stump = stump)
  }
  
  if(testfun == "ctree_cat_bin"){
    control <- partykit:::extree_control(criterion = "p.value",
                                         selectfun = .ctree_cat_bin_select(teststat = "quadratic", splitstat = "quadratic", 
                                                                           splittest = splittest, pargs = GenzBretz(),
                                                                           testtype = "MonteCarlo", nresample = 9999L, 
                                                                           tol = sqrt(.Machine$double.eps),
                                                                           intersplit = TRUE, MIA = FALSE,
                                                                           xgroups = xgroups),
                                         splitfun = partykit:::.objfun_split(restart = TRUE,
                                                                             intersplit = TRUE),
                                         svselectfun = .ctree_cat_bin_select(teststat = "quadratic", splitstat = "quadratic", 
                                                                             splittest = splittest, pargs = GenzBretz(),
                                                                             testtype = "MonteCarlo", nresample = 9999L, 
                                                                             tol = sqrt(.Machine$double.eps),
                                                                             intersplit = TRUE, MIA = FALSE,
                                                                             xgroups = xgroups), 
                                         svsplitfun = partykit:::.objfun_split(restart = TRUE,
                                                                               intersplit = TRUE),
                                         logmincriterion = log(1-alpha),
                                         update = TRUE,
                                         bonferroni = TRUE,
                                         stump = stump)
  }
  
  if(testfun == "mfluc"){
    control <- partykit:::extree_control(criterion = "p.value",
                                         selectfun = partykit:::.mfluc_select(breakties = FALSE, 
                                                                              intersplit = TRUE, parm = NULL, 
                                                                              dfsplit = TRUE, 
                                                                              restart = TRUE, model = TRUE, 
                                                                              vcov = "sandwich", 
                                                                              ordinal = "chisq", 
                                                                              ytype = "vector",
                                                                              nrep = 10000L, terminal = "object", 
                                                                              inner = "object", trim = 0.1),  
                                         splitfun = partykit:::.objfun_split(restart = TRUE,
                                                                             intersplit = TRUE),
                                         svselectfun = partykit:::.mfluc_select(breakties = FALSE, 
                                                                                intersplit = TRUE, parm = NULL, 
                                                                                dfsplit = TRUE, 
                                                                                restart = TRUE, model = TRUE, 
                                                                                vcov = "sandwich", 
                                                                                ordinal = "chisq", 
                                                                                ytype = "vector",
                                                                                nrep = 10000L, terminal = "object", 
                                                                                inner = "object", trim = 0.1), 
                                         svsplitfun = partykit:::.objfun_split(restart = TRUE,
                                                                               intersplit = TRUE),
                                         logmincriterion = log(1-alpha),
                                         update = TRUE,
                                         bonferroni = TRUE,
                                         stump = stump)
  }
  
  if(testfun == "mfluc_cat"){
    control <- partykit:::extree_control(criterion = "p.value",
                                         selectfun = .mfluc_cat_select(breakties = FALSE, 
                                                                       intersplit = TRUE, parm = NULL, 
                                                                       dfsplit = TRUE, 
                                                                       restart = TRUE, model = TRUE, 
                                                                       vcov = "sandwich", 
                                                                       ordinal = "chisq", 
                                                                       ytype = "vector",
                                                                       nrep = 10000L, terminal = "object", 
                                                                       inner = "object", trim = 0.1,
                                                                       xgroups = xgroups),  
                                         splitfun = partykit:::.objfun_split(restart = TRUE,
                                                                             intersplit = TRUE),
                                         svselectfun = .mfluc_cat_select(breakties = FALSE, 
                                                                         intersplit = TRUE, parm = NULL, 
                                                                         dfsplit = TRUE, 
                                                                         restart = TRUE, model = TRUE, 
                                                                         vcov = "sandwich", 
                                                                         ordinal = "chisq", 
                                                                         ytype = "vector",
                                                                         nrep = 10000L, terminal = "object", 
                                                                         inner = "object", trim = 0.1,
                                                                         xgroups = xgroups), 
                                         svsplitfun = partykit:::.objfun_split(restart = TRUE,
                                                                               intersplit = TRUE),
                                         logmincriterion = log(1-alpha),
                                         update = TRUE,
                                         bonferroni = TRUE,
                                         stump = stump)
  }
  
  if(testfun == "mfluc_bin"){
    control <- partykit:::extree_control(criterion = "p.value",
                                         selectfun = .mfluc_bin_select(breakties = FALSE, 
                                                                       intersplit = TRUE, parm = NULL, 
                                                                       dfsplit = TRUE, 
                                                                       restart = TRUE, model = TRUE, 
                                                                       vcov = "sandwich", 
                                                                       ordinal = "chisq", 
                                                                       ytype = "vector",
                                                                       nrep = 10000L, terminal = "object", 
                                                                       inner = "object", trim = 0.1,
                                                                       xgroups = xgroups),  
                                         splitfun = partykit:::.objfun_split(restart = TRUE,
                                                                             intersplit = TRUE),
                                         svselectfun = .mfluc_bin_select(breakties = FALSE, 
                                                                         intersplit = TRUE, parm = NULL, 
                                                                         dfsplit = TRUE, 
                                                                         restart = TRUE, model = TRUE, 
                                                                         vcov = "sandwich", 
                                                                         ordinal = "chisq", 
                                                                         ytype = "vector",
                                                                         nrep = 10000L, terminal = "object", 
                                                                         inner = "object", trim = 0.1,
                                                                         xgroups = xgroups), 
                                         svsplitfun = partykit:::.objfun_split(restart = TRUE,
                                                                               intersplit = TRUE),
                                         logmincriterion = log(1-alpha),
                                         update = TRUE,
                                         bonferroni = TRUE,
                                         stump = stump)
  }
  
  if(testfun == "mfluc_cat_bin"){
    control <- partykit:::extree_control(criterion = "p.value",
                                         selectfun = .mfluc_cat_bin_select(breakties = FALSE, 
                                                                           intersplit = TRUE, parm = NULL, 
                                                                           dfsplit = TRUE, 
                                                                           restart = TRUE, model = TRUE, 
                                                                           vcov = "sandwich", 
                                                                           ordinal = "chisq", 
                                                                           ytype = "vector",
                                                                           nrep = 10000L, terminal = "object", 
                                                                           inner = "object", trim = 0.1,
                                                                           xgroups = xgroups),  
                                         splitfun = partykit:::.objfun_split(restart = TRUE,
                                                                             intersplit = TRUE),
                                         svselectfun = .mfluc_cat_bin_select(breakties = FALSE, 
                                                                             intersplit = TRUE, parm = NULL, 
                                                                             dfsplit = TRUE, 
                                                                             restart = TRUE, model = TRUE, 
                                                                             vcov = "sandwich", 
                                                                             ordinal = "chisq", 
                                                                             ytype = "vector",
                                                                             nrep = 10000L, terminal = "object", 
                                                                             inner = "object", trim = 0.1,
                                                                             xgroups = xgroups), 
                                         svsplitfun = partykit:::.objfun_split(restart = TRUE,
                                                                               intersplit = TRUE),
                                         logmincriterion = log(1-alpha),
                                         update = TRUE,
                                         bonferroni = TRUE,
                                         stump = stump)
  }
  
  control$inner <- "object"
  control$terminal <- "object"
  
  ## set up model.frame() call
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$yx <- "matrix"
  mf$nmax <- control$nmax
  mf$ytype <- control$ytype
  ## evaluate model.frame
  mf[[1L]] <- quote(partykit::extree_data)
  
  d <- eval(mf, parent.frame())
  subset <- partykit:::.start_subset(d)
  
  weights <- model.weights(model.frame(d))
  
  if (is.null(control$update)) control$update <- TRUE
  
  # wrapper function to use lm() as trafo
  ytrafo <- function(subset, weights, estfun = TRUE, object = TRUE, info = NULL) {
    
    lmformula <- if(is.null(d$terms$yx)) Formula(d$terms$all) else Formula(d$terms$yx)   ## FIX ME: only regressors for lm, drop splitting variables
    sdata <- d$data[subset,]
    
    ## FIX ME: error in lm if weights are handed over
    ## FIX ME: scores with or without weights?
    #subweights <- if(is.null(weights) || (length(weights)==0L)) weights else weights[subset]
    
    model <- lm(formula = lmformula, data = sdata) #, weights = subweights) error in lm if weights are handed over
    
    ## FIX ME: add argument 'decorrelate' to control
    #decorrelate <- "vcov"
    #decorrelate <- "opg"
    #decorrelate <- "none"
    
    if(estfun) {
      ef <- as.matrix(sandwich::estfun(model)) 
      
      
      #if(ctrl$decorrelate != "none") {
      if(decorrelate != "none") {
        n <- NROW(ef)
        ef <- ef/sqrt(n)
        
        vcov <- if(decorrelate == "vcov") {
          vcov(model) * n
        } else {
          solve(crossprod(ef))
        }
        
        root.matrix <- function(X) {
          if((ncol(X) == 1L)&&(nrow(X) == 1L)) return(sqrt(X)) else {
            X.eigen <- eigen(X, symmetric = TRUE)
            if(any(X.eigen$values < 0)) stop("Matrix is not positive semidefinite")
            sqomega <- sqrt(diag(X.eigen$values))
            V <- X.eigen$vectors
            return(V %*% sqomega %*% t(V))
          }
        }
        ef <- as.matrix(t(root.matrix(vcov) %*% t(ef)))
      }
      
      estfun <- matrix(0, ncol = ncol(ef), nrow = nrow(d$data)) 
      estfun[subset,] <- ef
      ## FIX ME: return matrix of size nobs times par or only the scores from those observations selected by subset?
      ## FIX ME: estfun of lm-object returns unwheighted scores?
    
      if(whichscores == "residuals") estfun <- as.matrix(estfun[,1], ncol = 1)
      
    } else estfun <- NULL
    
    object <-  if(object) model else NULL
    
    ret <- list(estfun = estfun,
                unweighted = TRUE, # unweighted = TRUE would prevent estfun / w in extree_fit
                coefficients = coef(model),
                objfun = -logLik(model),  # optional function to be minimized 
                object = object,
                decorrelated = (decorrelate != "none"),    
                converged = !is.null(model)  # FIX ME: better check whether coefficients are returned ?
    )
    return(ret)
  }
  
  
  converged <- TRUE
  control$model <- TRUE
  
  update <- function(subset, weights, control, doFit = TRUE)
    extree_fit(data = d, trafo = ytrafo, converged = converged, partyvars = d$variables$z, 
               subset = subset, weights = weights, ctrl = control, doFit = doFit)
  if (!doFit) return(list(d = d, update = update))
  
  tree <- update(subset = subset, weights = weights, control = control)
  trafo <- tree$trafo
  
  ### prepare as modelparty
  mf <- model.frame(d)
  if (length(weights) == 0) weights <- rep(1, nrow(mf))
  
  fitted <- data.frame("(fitted)" = fitted_node(tree$nodes, mf),
                       "(weights)" = weights,
                       check.names = FALSE)
  
  fitted[[3]] <- y <- mf[, d$variables$y, drop = TRUE]
  names(fitted)[3] <- "(response)"
  
  control$ytype <- ifelse(is.vector(y), "vector", class(y))
  # x <- model.matrix(modelf, data = mmf)
  control$xtype <- "matrix" # TODO: find out when to use data.frame
  
  ## return party object
  rval <- party(tree$nodes, 
                data = if(control$model) mf else mf[0,],
                fitted = fitted,
                terms = d$terms$all,
                info = list(
                  call = match.call(),
                  formula = formula,
                  Formula = as.Formula(formula),
                  terms = list(response = d$terms$yx, partitioning = d$terms$z),
                  fit = ytrafo,
                  control = control,
                  dots = list(),
                  nreg = NCOL(d$yx$x)
                )
  )
  class(rval) <- c("modelparty", class(rval))
  
  ### add modelinfo (object) and estfun if not there yet, but wanted
  # TODO: check if this can be done prettier
  which_terminals <- nodeids(rval, terminal = TRUE)
  which_all <- nodeids(rval)
  
  idx <- lapply(which_all, partykit:::.get_path, obj = tree$nodes)
  names(idx) <- which_all
  tree_ret <- unclass(rval)
  subset_term <- predict(rval, type = "node")
  
  for (i in which_all) {
    ichar <- as.character(i)
    iinfo <- tree_ret[[c(1, idx[[ichar]])]]$info
    
    if (i %in% which_terminals) winfo <- "object" else winfo <- control$inner
    #if (i %in% which_terminals) winfo <- control$terminal else winfo <- control$inner
    
    if (is.null(winfo)) {
      iinfo$object <- NULL
      iinfo$estfun <- NULL
    } else {
      if (is.null(iinfo) | any(is.null(iinfo[[winfo]])) | 
          any(! winfo %in% names(iinfo))) {
        iinfo <- trafo(subset = which(subset_term == i), weights = weights, info = NULL,
                       estfun = ("estfun" %in% winfo),
                       object = ("object" %in% winfo))
      }
    }
    
    tree_ret[[c(1, idx[[ichar]])]]$info <- iinfo
  }
  
  class(tree_ret) <- class(rval)
  
  return(tree_ret)
  
  if(FALSE){
    tree <- tree$nodes
    
    mf <- model.frame(d)
    if (is.null(weights)) weights <- rep(1, nrow(mf))
    
    fitted <- data.frame("(fitted)" = fitted_node(tree, mf), 
                         "(weights)" = weights,
                         check.names = FALSE)
    fitted[[3]] <- mf[, d$variables$y, drop = TRUE]
    names(fitted)[3] <- "(response)"
    ret <- party(tree, data = mf, fitted = fitted, 
                 info = list(call = match.call(), control = control))
    ret$update <- update
    ret$trafo <- trafo
    class(ret) <- c("constparty", class(ret))
    
    ### doesn't work for Surv objects
    # ret$terms <- terms(formula, data = mf)
    ret$terms <- d$terms$all
    ### need to adjust print and plot methods
    ### for multivariate responses
    ### if (length(response) > 1) class(ret) <- "party"
    return(ret)
  }
  
}


############
# compare on one data set
if(FALSE){
  d <- dgp_stump(400, vary_beta = "beta1", beta0 = 3)
  d <- dgp_stump(400, vary_beta = "all", xi = -0.5)
  d <- dgp_stump(1000, vary_beta = "all", xi = -0.5, binary_regressor = FALSE)
  d <- dgp_stump(1000, vary_beta = "all", xi = -0.0, binary_regressor = FALSE)
  
  d <- dgp_stump(400, vary_beta = "all", xi = -0.4, delta = 3, binary_regressor = TRUE)
  d <- dgp_stump(200, vary_beta = "beta1", beta0 = 0, xi = 0.0, binary_regressor = FALSE)
  d <- dgp_stump(200, vary_beta = "beta1", beta0 = 0, xi = -0.5, binary_regressor = FALSE)
  d <- dgp_stump(200, vary_beta = "beta1", beta0 = 0, xi = 0.0, binary_regressor = FALSE)
  
  d <- dgp_stump(200, vary_beta = "beta1", beta0 = 0, xi = 0.0, binary_regressor = FALSE, binary_beta = TRUE, delta = 5)
  d <- dgp_stump(200, vary_beta = "beta1", beta0 = 0, xi = 0.0, binary_regressor = FALSE, binary_beta = FALSE, delta = 5)
  
  d <- dgp_stump(200, vary_beta = "all", xi = 0.7, binary_regressor = FALSE, binary_beta = TRUE, delta = 1, z1dist = "unif")
  d <- dgp_stump(nobs = 1000, xi = 0.3, binary_regressor = FALSE, vary_beta = "all", nrsteps = 2, nrlevels = 3, binary_beta = TRUE)
  
  d <- dgp_stump(nobs = 250, xi = 0.3, binary_regressor = FALSE, vary_beta = "all", nrsteps = 2, nrlevels = 2, binary_beta = TRUE)
  
  d <- dgp_tree(nobs = 400, xi = c(0,0.5), binary_regressor = TRUE, vary_beta = "all", 
                 binary_beta = TRUE, nrsplitvar = 2)
  d <- dgp_tree(nobs = 400, xi = c(0,-0.4), binary_regressor = FALSE, vary_beta = "all", 
                binary_beta = TRUE, nrsplitvar = 2, sigma = 0.1)
  d <- dgp_tree(nobs = 400, xi = c(0), binary_regressor = FALSE, vary_beta = "all", 
                binary_beta = TRUE, nrsplitvar = 2, sigma = 0.1)
  
  plot(d$x[d$id==1], d$y[d$id==1], ylim = c(min(d$y), max(d$y)))
  points(d$x[d$id==2], d$y[d$id==2], col = "red")
  points(d$x[d$id==3], d$y[d$id==3], col = "green")
  l <- lm(y~x, data = d)
  abline(l)
  
  ctest <- evaltests(y~x|z1+z2+z3+z4+z5+z6+z7+z8+z9+z10, data = d, testfun = "ctree", stump = FALSE)
  mtest <- evaltests(y~x|z1+z2+z3+z4+z5+z6+z7+z8+z9+z10, data = d, testfun = "mfluc", stump = FALSE)
  cctest <- evaltests(y~x|z1+z2+z3+z4+z5+z6+z7+z8+z9+z10, data = d, testfun = "ctree_cat", stump = FALSE)
  mctest <- evaltests(y~x|z1+z2+z3+z4+z5+z6+z7+z8+z9+z10, data = d, testfun = "mfluc_cat", stump = FALSE)
  gstest12 <- evaltests(y~x|z1+z2+z3+z4+z5+z6+z7+z8+z9+z10, data = d, testfun = "guide", 
                        guide_testtype = "sum", guide_parm = c(1,2),
                        xgroups = NULL, ygroups = NULL, stump = FALSE)
  gmtest12 <- evaltests(y~x|z1+z2+z3+z4+z5+z6+z7+z8+z9+z10, data = d, testfun = "guide", 
                        guide_testtype = "max", guide_parm = c(1,2),
                        xgroups = NULL, ygroups = NULL, stump = FALSE)
  gctest12 <- evaltests(y~x|z1+z2+z3+z4+z5+z6+z7+z8+z9+z10, data = d, testfun = "guide", 
                        guide_testtype = "coin", guide_parm = c(1,2),
                        xgroups = NULL, ygroups = NULL, stump = FALSE)
  gstest1 <- evaltests(y~x|z1+z2+z3+z4+z5+z6+z7+z8+z9+z10, data = d, testfun = "guide", 
                       guide_testtype = "sum", guide_parm = c(1),
                       xgroups = NULL, ygroups = NULL, stump = FALSE)
  gmtest1 <- evaltests(y~x|z1+z2+z3+z4+z5+z6+z7+z8+z9+z10, data = d, testfun = "guide", 
                       guide_testtype = "max", guide_parm = c(1),
                       xgroups = NULL, ygroups = NULL, stump = FALSE)
  gctest1 <- evaltests(y~x|z1+z2+z3+z4+z5+z6+z7+z8+z9+z10, data = d, testfun = "guide", 
                       guide_testtype = "coin", guide_parm = c(1),
                       xgroups = NULL, ygroups = NULL, stump = FALSE)
  gstest2 <- evaltests(y~x|z1+z2+z3+z4+z5+z6+z7+z8+z9+z10, data = d, testfun = "guide", 
                       guide_testtype = "sum", guide_parm = c(2),
                       xgroups = NULL, ygroups = NULL, stump = FALSE)
  gmtest2 <- evaltests(y~x|z1+z2+z3+z4+z5+z6+z7+z8+z9+z10, data = d, testfun = "guide", 
                       guide_testtype = "max", guide_parm = c(2),
                       xgroups = NULL, ygroups = NULL, stump = FALSE)
  gctest2 <- evaltests(y~x|z1+z2+z3+z4+z5+z6+z7+z8+z9+z10, data = d, testfun = "guide", 
                       guide_testtype = "coin", guide_parm = c(2),
                       xgroups = NULL, ygroups = NULL, stump = FALSE)
  
  
  plot(ctest)
  plot(mtest)
  plot(cctest)
  plot(mctest)
  plot(gmtest12)
  plot(gstest12)
  plot(gctest12)
  plot(gmtest1)
  plot(gstest1)
  plot(gctest1)
  plot(gmtest2)
  plot(gstest2)
  plot(gctest2)
  
  plot(ctest, terminal_panel = node_bivplot)
  plot(mtest, terminal_panel = node_bivplot)
  plot(cctest, terminal_panel = node_bivplot)
  plot(mctest, terminal_panel = node_bivplot)
  plot(gmtest12, terminal_panel = node_bivplot)
  plot(gstest12, terminal_panel = node_bivplot)
  plot(gctest12, terminal_panel = node_bivplot)
  plot(gmtest1, terminal_panel = node_bivplot)
  plot(gstest1, terminal_panel = node_bivplot)
  plot(gctest1, terminal_panel = node_bivplot)
  plot(gmtest2, terminal_panel = node_bivplot)
  plot(gstest2, terminal_panel = node_bivplot)
  plot(gctest2, terminal_panel = node_bivplot)
}




###################################################
# various data sets


sim <- function(nobs = 100, nrep = 100, seed = 7, stump = TRUE, nrsteps = 1,
                nrlevels = 2, shift = 0, nrsplitvar = 1,
                formula = y~x|z1+z2+z3+z4+z5+z6+z7+z8+z9+z10,
                vary_beta = "all", beta0 = 0, beta1 = 1, xi = 0, delta = 1, 
                binary_regressor = FALSE, binary_beta = TRUE, z1dist = "unif",
                sigma = 1, only_intercept = FALSE, alpha = 0.05,
                return_matrices = FALSE,
                compare_pruning = FALSE,
                test = c("ctree", "mfluc", "ctree_max", 
                         "ctree_cat", "ctree_max_cat", "mfluc_cat",
                         "ctree_bin", "ctree_max_bin", "mfluc_bin",
                         "ctree_cat_bin", "ctree_max_cat_bin", "mfluc_cat_bin",
                         "ctree_resid", "mfluc_resid", "ctree_resid_max", 
                         "ctree_resid_cat", "ctree_resid_max_cat", "mfluc_resid_cat",
                         "ctree_resid_bin", "ctree_resid_max_bin", "mfluc_resid_bin",
                         "ctree_resid_cat_bin", "ctree_resid_max_cat_bin", "mfluc_resid_cat_bin",
                         "guide_sum_12", "guide_max_12", "guide_coin_12",
                         "guide_sum_1", "guide_max_1", "guide_coin_1",
                         "guide_sum_2", "guide_max_2", "guide_coin_2",
                         "guide_sum_1_cor", "guide_max_1_cor", "guide_coin_1_cor",
                         "guide_sum_2_cor", "guide_max_2_cor", "guide_coin_2_cor"))
{
  set.seed(seed)
  
  ## call
  cl <- match.call()
  
  ## check input arguments for comparing different pruning strategies
  if(compare_pruning){
    if(stump) stop("different pruning strategies can only be compared for stump = FALSE")
  }
  
  ## FIX ME: better way of handing over two split points (vector xi) for stump=FALSE
  if(!stump & length(xi)<2 & nrsplitvar == 1){
    if(xi == -0.5) xi <- c(-0.5, 0)
    if(xi == -0.4) xi <- c(-0.4, 0.3)
    if(xi == 0) xi <- c(0, 0.5)
  }

    
  if(stump){
    pval_z1 <- matrix(rep(NA, length(test)*nrep), ncol = length(test))
    pval <- matrix(rep(NA, length(test)*nrep), ncol = length(test))
    sv <- matrix(rep(NA, length(test)*nrep), ncol = length(test))
    colnames(pval) <- colnames(pval_z1) <- colnames(sv) <- test
  } else {
    nrsubgr <- matrix(rep(NA, length(test)*nrep), ncol = length(test))
    colnames(nrsubgr) <- test
  }
  err_coef <- matrix(rep(NA, length(test)*nrep), ncol = length(test))
  err_total <- matrix(rep(NA, length(test)*nrep), ncol = length(test))
  ari <- matrix(rep(NA, length(test)*nrep), ncol = length(test))
  colnames(err_coef) <- colnames(err_total) <- colnames(ari) <- test
  
  if(compare_pruning) {
    if(stump){
      pval_z1_p <- matrix(rep(NA, length(test)*nrep), ncol = length(test))
      pval_p <- matrix(rep(NA, length(test)*nrep), ncol = length(test))
      sv_p <- matrix(rep(NA, length(test)*nrep), ncol = length(test))
      colnames(pval_p) <- colnames(pval_z1_p) <- colnames(sv_p) <- test
    } else {
      nrsubgr_p <- matrix(rep(NA, length(test)*nrep), ncol = length(test))
      colnames(nrsubgr_p) <- test
    }
    err_coef_p <- matrix(rep(NA, length(test)*nrep), ncol = length(test))
    err_total_p <- matrix(rep(NA, length(test)*nrep), ncol = length(test))
    ari_p <- matrix(rep(NA, length(test)*nrep), ncol = length(test))
    colnames(err_coef_p) <- colnames(err_total_p) <- colnames(ari_p) <- test
  }  
  
  
  for(i in 1:nrep){
    
    dgp <- if(stump) dgp_stump else dgp_tree
    d <- dgp(nobs = nobs, vary_beta = vary_beta, beta0 = beta0, beta1 = beta1, xi = xi, 
             binary_regressor = binary_regressor, delta = delta,
             binary_beta = binary_beta, z1dist = z1dist, nrsteps = nrsteps,
             nrlevels = nrlevels, shift = shift, nrsplitvar = nrsplitvar,
             seed = seed + i, sigma = sigma, only_intercept = only_intercept)
    
    # FIXME: different name for function (get_resval ?)
    compute_pval <- function(test) {
      test <- match.arg(test, c("ctree", "mfluc", "ctree_max", 
                                "ctree_cat", "ctree_max_cat", "mfluc_cat",
                                "ctree_bin", "ctree_max_bin", "mfluc_bin",
                                "ctree_cat_bin", "ctree_max_cat_bin", "mfluc_cat_bin",
                                "ctree_resid", "mfluc_resid", "ctree_resid_max", 
                                "ctree_resid_cat", "ctree_resid_max_cat", "mfluc_resid_cat",
                                "ctree_resid_bin", "ctree_resid_max_bin", "mfluc_resid_bin",
                                "ctree_resid_cat_bin", "ctree_resid_max_cat_bin", "mfluc_resid_cat_bin",
                                "guide_sum_12", "guide_max_12", "guide_coin_12",
                                "guide_sum_1", "guide_max_1", "guide_coin_1",
                                "guide_sum_2", "guide_max_2", "guide_coin_2",
                                "guide_sum_1_cor", "guide_max_1_cor", "guide_coin_1_cor",
                                "guide_sum_2_cor", "guide_max_2_cor", "guide_coin_2_cor"))
      testres <- switch(test,
                        "ctree" = evaltests(formula, data = d, testfun = "ctree", stump = stump, decorrelate = "vcov"),
                        "ctree_max" = evaltests(formula, data = d, testfun = "ctree", stump = stump, decorrelate = "vcov", ctree_max = TRUE),
                        "mfluc" = evaltests(formula, data = d, testfun = "mfluc", stump = stump, decorrelate = "vcov"),
                        
                        "ctree_cat" = evaltests(formula, data = d, testfun = "ctree_cat", stump = stump, decorrelate = "vcov"),
                        "ctree_max_cat" = evaltests(formula, data = d, testfun = "ctree_cat", stump = stump, decorrelate = "vcov", ctree_max = TRUE),
                        "mfluc_cat" = evaltests(formula, data = d, testfun = "mfluc_cat", stump = stump, decorrelate = "vcov"),
                        
                        "ctree_bin" = evaltests(formula, data = d, testfun = "ctree_bin", stump = stump, decorrelate = "vcov"),
                        "ctree_max_bin" = evaltests(formula, data = d, testfun = "ctree_bin", stump = stump, decorrelate = "vcov", ctree_max = TRUE),
                        "mfluc_bin" = evaltests(formula, data = d, testfun = "mfluc_bin", stump = stump, decorrelate = "vcov"),
                        
                        "ctree_cat_bin" = evaltests(formula, data = d, testfun = "ctree_cat_bin", stump = stump, decorrelate = "vcov"),
                        "ctree_max_cat_bin" = evaltests(formula, data = d, testfun = "ctree_cat_bin", stump = stump, decorrelate = "vcov", ctree_max = TRUE),
                        "mfluc_cat_bin" = evaltests(formula, data = d, testfun = "mfluc_cat_bin", stump = stump, decorrelate = "vcov"),
                        
                        "ctree_resid" = evaltests(formula, data = d, testfun = "ctree", stump = stump, decorrelate = "vcov", whichscores = "residuals"),
                        "ctree_resid_max" = evaltests(formula, data = d, testfun = "ctree", stump = stump, decorrelate = "vcov", ctree_max = TRUE, whichscores = "residuals"),
                        "mfluc_resid" = evaltests(formula, data = d, testfun = "mfluc", stump = stump, decorrelate = "vcov", whichscores = "residuals"),
                        
                        "ctree_resid_cat" = evaltests(formula, data = d, testfun = "ctree_cat", stump = stump, decorrelate = "vcov", whichscores = "residuals"),
                        "ctree_resid_max_cat" = evaltests(formula, data = d, testfun = "ctree_cat", stump = stump, decorrelate = "vcov", ctree_max = TRUE, whichscores = "residuals"),
                        "mfluc_resid_cat" = evaltests(formula, data = d, testfun = "mfluc_cat", stump = stump, decorrelate = "vcov", whichscores = "residuals"),
                        
                        "ctree_resid_bin" = evaltests(formula, data = d, testfun = "ctree_bin", stump = stump, decorrelate = "vcov", whichscores = "residuals"),
                        "ctree_resid_max_bin" = evaltests(formula, data = d, testfun = "ctree_bin", stump = stump, decorrelate = "vcov", ctree_max = TRUE, whichscores = "residuals"),
                        "mfluc_resid_bin" = evaltests(formula, data = d, testfun = "mfluc_bin", stump = stump, decorrelate = "vcov", whichscores = "residuals"),
                        
                        "ctree_resid_cat_bin" = evaltests(formula, data = d, testfun = "ctree_cat_bin", stump = stump, decorrelate = "vcov", whichscores = "residuals"),
                        "ctree_resid_max_cat_bin" = evaltests(formula, data = d, testfun = "ctree_cat_bin", stump = stump, decorrelate = "vcov", ctree_max = TRUE, whichscores = "residuals"),
                        "mfluc_resid_cat_bin" = evaltests(formula, data = d, testfun = "mfluc_cat_bin", stump = stump, decorrelate = "vcov", whichscores = "residuals"),
                        
                        "guide_sum_12" = evaltests(formula, data = d, testfun = "guide", 
                                                   guide_testtype = "sum", guide_parm = c(1,2),
                                                   xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "vcov"),
                        "guide_max_12" = evaltests(formula, data = d, testfun = "guide", 
                                                   guide_testtype = "max", guide_parm = c(1,2),
                                                   xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "vcov"),
                        "guide_coin_12" = evaltests(formula, data = d, testfun = "guide", 
                                                    guide_testtype = "coin", guide_parm = c(1,2),
                                                    xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "vcov"),
                        "guide_sum_1" = evaltests(formula, data = d, testfun = "guide", 
                                                  guide_testtype = "sum", guide_parm = c(1),
                                                  xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "vcov"),
                        "guide_max_1" = evaltests(formula, data = d, testfun = "guide", 
                                                  guide_testtype = "max", guide_parm = c(1),
                                                  xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "vcov"),
                        "guide_coin_1" = evaltests(formula, data = d, testfun = "guide", 
                                                   guide_testtype = "coin", guide_parm = c(1),
                                                   xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "vcov"),
                        "guide_sum_2" = evaltests(formula, data = d, testfun = "guide", 
                                                  guide_testtype = "sum", guide_parm = c(2),
                                                  xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "vcov"),
                        "guide_max_2" = evaltests(formula, data = d, testfun = "guide", 
                                                  guide_testtype = "max", guide_parm = c(2),
                                                  xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "vcov"),
                        "guide_coin_2" = evaltests(formula, data = d, testfun = "guide", 
                                                   guide_testtype = "coin", guide_parm = c(2),
                                                   xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "vcov"),
                        "guide_sum_1_cor" = evaltests(formula, data = d, testfun = "guide", 
                                                      guide_testtype = "sum", guide_parm = c(1),
                                                      xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "none"),
                        "guide_max_1_cor" = evaltests(formula, data = d, testfun = "guide", 
                                                      guide_testtype = "max", guide_parm = c(1),
                                                      xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "none"),
                        "guide_coin_1_cor" = evaltests(formula, data = d, testfun = "guide", 
                                                       guide_testtype = "coin", guide_parm = c(1),
                                                       xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "none"),
                        "guide_sum_2_cor" = evaltests(formula, data = d, testfun = "guide", 
                                                      guide_testtype = "sum", guide_parm = c(2),
                                                      xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "none"),
                        "guide_max_2_cor" = evaltests(formula, data = d, testfun = "guide", 
                                                      guide_testtype = "max", guide_parm = c(2),
                                                      xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "none"),
                        "guide_coin_2_cor" = evaltests(formula, data = d, testfun = "guide", 
                                                       guide_testtype = "coin", guide_parm = c(2),
                                                       xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "none")
      )
      
      predresp <- predict(testres, type = "response")
      prednode <- predict(testres, type = "node")
      predcoef <- if(is.matrix(coef(testres))) coef(testres)[as.character(prednode),] else coef(testres)
      
      # calculate error of predicted coefficients by taking the euclidean norm of the differences of predicted and true coefficients,
      # averaged over all observations
      err_coef <- mean(sqrt(rowSums((predcoef - d[,c("beta0", "beta1")])^2)))
      # mean squared error of predicted responses
      err_total <- mean(sqrt((predresp - d[,"y"])^2))
      
      ari <- as.numeric(adj_rand_index(prednode, d$id))
      
      if(stump){
        returnlist <- list(pval = as.numeric(info_node(testres$node)$p.value),
                           pval_z1 = as.numeric(info_node(testres$node)$criterion["p.value","z1"]),
                           sv = names(info_node(testres$node)$p.value),
                           err_coef = err_coef,
                           err_total = err_total,
                           ari = ari)
      } else {
        returnlist <- list(err_coef = err_coef,
                           err_total = err_total,
                           nrsubgr = width(testres),
                           ari = ari)
      }
      
      
      ## if the same tree should be built using post-pruning
      if(compare_pruning){
        testres_p <- switch(test,
                          "ctree" = evaltests(formula, data = d, testfun = "ctree", stump = stump, decorrelate = "vcov", alpha = 1),
                          "ctree_max" = evaltests(formula, data = d, testfun = "ctree", stump = stump, decorrelate = "vcov", ctree_max = TRUE, alpha = 1),
                          "mfluc" = evaltests(formula, data = d, testfun = "mfluc", stump = stump, decorrelate = "vcov", alpha = 1),
                          
                          "ctree_cat" = evaltests(formula, data = d, testfun = "ctree_cat", stump = stump, decorrelate = "vcov", alpha = 1),
                          "ctree_max_cat" = evaltests(formula, data = d, testfun = "ctree_cat", stump = stump, decorrelate = "vcov", ctree_max = TRUE, alpha = 1),
                          "mfluc_cat" = evaltests(formula, data = d, testfun = "mfluc_cat", stump = stump, decorrelate = "vcov", alpha = 1),
                          
                          "ctree_bin" = evaltests(formula, data = d, testfun = "ctree_bin", stump = stump, decorrelate = "vcov", alpha = 1),
                          "ctree_max_bin" = evaltests(formula, data = d, testfun = "ctree_bin", stump = stump, decorrelate = "vcov", ctree_max = TRUE, alpha = 1),
                          "mfluc_bin" = evaltests(formula, data = d, testfun = "mfluc_bin", stump = stump, decorrelate = "vcov", alpha = 1),
                          
                          "ctree_cat_bin" = evaltests(formula, data = d, testfun = "ctree_cat_bin", stump = stump, decorrelate = "vcov", alpha = 1),
                          "ctree_max_cat_bin" = evaltests(formula, data = d, testfun = "ctree_cat_bin", stump = stump, decorrelate = "vcov", ctree_max = TRUE, alpha = 1),
                          "mfluc_cat_bin" = evaltests(formula, data = d, testfun = "mfluc_cat_bin", stump = stump, decorrelate = "vcov", alpha = 1),
                          
                          "ctree_resid" = evaltests(formula, data = d, testfun = "ctree", stump = stump, decorrelate = "vcov", whichscores = "residuals", alpha = 1),
                          "ctree_resid_max" = evaltests(formula, data = d, testfun = "ctree", stump = stump, decorrelate = "vcov", ctree_max = TRUE, whichscores = "residuals", alpha = 1),
                          "mfluc_resid" = evaltests(formula, data = d, testfun = "mfluc", stump = stump, decorrelate = "vcov", whichscores = "residuals", alpha = 1),
                          
                          "ctree_resid_cat" = evaltests(formula, data = d, testfun = "ctree_cat", stump = stump, decorrelate = "vcov", whichscores = "residuals", alpha = 1),
                          "ctree_resid_max_cat" = evaltests(formula, data = d, testfun = "ctree_cat", stump = stump, decorrelate = "vcov", ctree_max = TRUE, whichscores = "residuals", alpha = 1),
                          "mfluc_resid_cat" = evaltests(formula, data = d, testfun = "mfluc_cat", stump = stump, decorrelate = "vcov", whichscores = "residuals", alpha = 1),
                          
                          "ctree_resid_bin" = evaltests(formula, data = d, testfun = "ctree_bin", stump = stump, decorrelate = "vcov", whichscores = "residuals", alpha = 1),
                          "ctree_resid_max_bin" = evaltests(formula, data = d, testfun = "ctree_bin", stump = stump, decorrelate = "vcov", ctree_max = TRUE, whichscores = "residuals", alpha = 1),
                          "mfluc_resid_bin" = evaltests(formula, data = d, testfun = "mfluc_bin", stump = stump, decorrelate = "vcov", whichscores = "residuals", alpha = 1),
                          
                          "ctree_resid_cat_bin" = evaltests(formula, data = d, testfun = "ctree_cat_bin", stump = stump, decorrelate = "vcov", whichscores = "residuals", alpha = 1),
                          "ctree_resid_max_cat_bin" = evaltests(formula, data = d, testfun = "ctree_cat_bin", stump = stump, decorrelate = "vcov", ctree_max = TRUE, whichscores = "residuals", alpha = 1),
                          "mfluc_resid_cat_bin" = evaltests(formula, data = d, testfun = "mfluc_cat_bin", stump = stump, decorrelate = "vcov", whichscores = "residuals", alpha = 1),
                          
                          "guide_sum_12" = evaltests(formula, data = d, testfun = "guide", 
                                                     guide_testtype = "sum", guide_parm = c(1,2),
                                                     xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "vcov", alpha = 1),
                          "guide_max_12" = evaltests(formula, data = d, testfun = "guide", 
                                                     guide_testtype = "max", guide_parm = c(1,2),
                                                     xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "vcov", alpha = 1),
                          "guide_coin_12" = evaltests(formula, data = d, testfun = "guide", 
                                                      guide_testtype = "coin", guide_parm = c(1,2),
                                                      xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "vcov", alpha = 1),
                          "guide_sum_1" = evaltests(formula, data = d, testfun = "guide", 
                                                    guide_testtype = "sum", guide_parm = c(1),
                                                    xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "vcov", alpha = 1),
                          "guide_max_1" = evaltests(formula, data = d, testfun = "guide", 
                                                    guide_testtype = "max", guide_parm = c(1),
                                                    xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "vcov", alpha = 1),
                          "guide_coin_1" = evaltests(formula, data = d, testfun = "guide", 
                                                     guide_testtype = "coin", guide_parm = c(1),
                                                     xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "vcov", alpha = 1),
                          "guide_sum_2" = evaltests(formula, data = d, testfun = "guide", 
                                                    guide_testtype = "sum", guide_parm = c(2),
                                                    xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "vcov", alpha = 1),
                          "guide_max_2" = evaltests(formula, data = d, testfun = "guide", 
                                                    guide_testtype = "max", guide_parm = c(2),
                                                    xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "vcov", alpha = 1),
                          "guide_coin_2" = evaltests(formula, data = d, testfun = "guide", 
                                                     guide_testtype = "coin", guide_parm = c(2),
                                                     xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "vcov", alpha = 1),
                          "guide_sum_1_cor" = evaltests(formula, data = d, testfun = "guide", 
                                                        guide_testtype = "sum", guide_parm = c(1),
                                                        xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "none", alpha = 1),
                          "guide_max_1_cor" = evaltests(formula, data = d, testfun = "guide", 
                                                        guide_testtype = "max", guide_parm = c(1),
                                                        xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "none", alpha = 1),
                          "guide_coin_1_cor" = evaltests(formula, data = d, testfun = "guide", 
                                                         guide_testtype = "coin", guide_parm = c(1),
                                                         xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "none", alpha = 1),
                          "guide_sum_2_cor" = evaltests(formula, data = d, testfun = "guide", 
                                                        guide_testtype = "sum", guide_parm = c(2),
                                                        xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "none", alpha = 1),
                          "guide_max_2_cor" = evaltests(formula, data = d, testfun = "guide", 
                                                        guide_testtype = "max", guide_parm = c(2),
                                                        xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "none", alpha = 1),
                          "guide_coin_2_cor" = evaltests(formula, data = d, testfun = "guide", 
                                                         guide_testtype = "coin", guide_parm = c(2),
                                                         xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "none", alpha = 1)
        )
        
        testres_p <- ccprune(testres_p, nrfolds = 5)$Tree
        
        #predresp_p <- predict(testres_p, type = "response")
        ## FIX ME: newdata has to be set in order to get correct node ids
        prednode_p <- predict(testres_p, type = "node", newdata = testres_p$data)
        predcoef_p <- if(is.matrix(coef(testres_p))) coef(testres_p)[as.character(prednode_p),] else coef(testres_p)
        if(is.vector(predcoef_p)) predcoef_p <- matrix(predcoef_p, nrow = 1, ncol = 2)
        predresp_p <- predcoef_p[,1]  + predcoef_p[,2] * testres_p$data$x
        
        # calculate error of predicted coefficients by taking the euclidean norm of the differences of predicted and true coefficients,
        # averaged over all observations
        err_coef_p <- mean(sqrt(rowSums((predcoef_p - d[,c("beta0", "beta1")])^2)))
        # mean squared error of predicted responses
        err_total_p <- mean(sqrt((predresp_p - d[,"y"])^2))
        
        ari_p <- as.numeric(adj_rand_index(prednode_p, d$id))
        
        
        returnlist <- c(returnlist,
                        list(err_coef_p = err_coef_p,
                             err_total_p = err_total_p,
                             nrsubgr_p = width(testres_p),
                             ari_p = ari_p))
      }
      
      return(returnlist)
    }
    
    resmat <- sapply(test, compute_pval)
    
    if(stump){
      pval[i,] <- unlist(resmat["pval",])
      pval_z1[i,] <- unlist(resmat["pval_z1",])
      sv[i,] <- unlist(resmat["sv",])
    } else {
      nrsubgr[i,] <- unlist(resmat["nrsubgr",])
    }
    err_coef[i,] <- unlist(resmat["err_coef",])
    err_total[i,] <- unlist(resmat["err_total",])
    ari[i,] <- unlist(resmat["ari",])
    
    if(compare_pruning){
      
      nrsubgr_p[i,] <- unlist(resmat["nrsubgr_p",])
      
      err_coef_p[i,] <- unlist(resmat["err_coef_p",])
      err_total_p[i,] <- unlist(resmat["err_total_p",])
      ari_p[i,] <- unlist(resmat["ari_p",])
    }
    
    
  }
  
  
  if(stump){
    # proportion of cases where a split is / is not found
    prop_nosplit <- colMeans(pval > alpha)
    prop_split <- colMeans(pval <= alpha)
    
    # proportion of cases where the correct variable (z1) has the smallest p-value (but not necessarily smaller than alpha)
    prop_z1 <- colMeans(sv =="z1")
    
    # proportion of cases where a split is found in the correct variable (z1)
    pval_T <- pval
    pval_T[which(sv !="z1")] <- 1
    prop_Tsplit <- colMeans(pval_T < alpha)
    
    # proportion of cases where a split is found in a noise variable (z2, ..., z10)
    pval_F <- pval
    pval_F[which(sv =="z1")] <- 1
    prop_Fsplit <- colMeans(pval_F < alpha)
    
    # average of smallest p-value (regardless of the corresponding splitting variable)
    pval_min <- colMeans(pval)
    
    # average of p-value for the correct variable (z1)
    pval_true <- colMeans(pval_z1)
    
  } else {
    # average of number of subgroups
    nrsubgr <- colMeans(nrsubgr)
  }
  
  # average of errors and adjusted rand index
  err_coef <- colMeans(err_coef)
  err_total <- colMeans(err_total)
  ari <- colMeans(ari)
  
  ## FIX ME: which list to return for stump = FALSE
  if(stump){
    returnlist <- list(prop_nosplit = prop_nosplit,
                       prop_split = prop_split,
                       prop_Fsplit = prop_Fsplit,
                       prop_Tsplit = prop_Tsplit,
                       prop_z1 = prop_z1,
                       pval_min = pval_min,
                       pval_true = pval_true, 
                       err_coef = err_coef,
                       err_total = err_total,
                       ari = ari)
  } else {
    returnlist <- list(nrsubgr = nrsubgr,
                       err_coef = err_coef,
                       err_total = err_total,
                       ari = ari)
  }
  
  if(return_matrices) returnlist <- c(returnlist, list(sv = sv,
                                                       pval = pval,
                                                       pval_z1 = pval_z1))
  
  
  if(compare_pruning){
    
    # average of number of subgroups
    nrsubgr_p <- colMeans(nrsubgr_p)
    
    # average of errors and adjusted rand index
    err_coef_p <- colMeans(err_coef_p)
    err_total_p <- colMeans(err_total_p)
    ari_p <- colMeans(ari_p)
    
    ## FIX ME: which list to return for stump = FALSE
    returnlist <- c(returnlist, 
                    list(nrsubgr_p = nrsubgr_p,
                         err_coef_p = err_coef_p,
                         err_total_p = err_total_p,
                         ari_p = ari_p))
  }
  
  return(returnlist)
}





if(FALSE){
  simp <- sim(250, nrep = 2, stump = FALSE, nrsplitvar = 2, 
              delta = 1, binary_regressor = FALSE, 
              compare_pruning = TRUE, test = c("ctree", "mfluc", "guide_sum_12"))
  
  
  
  simtest <- sim(nobs = 200, nrep = 5, seed = 7,
                 formula = y~x|z1+z2+z3+z4+z5+z6+z7+z8+z9+z10,
                 vary_beta = "all", xi = c(-0.4,0.3), delta = 1,
                 stump = FALSE, binary_regressor = FALSE,
                 binary_beta = TRUE, z1dist = "unif",
                 sigma = 1, only_intercept = FALSE,
                 test = c("ctree","mfluc", "guide_max_1_cor", "guide_coin_12", "guide_max_2_cor"))
}


############# wrapper function applying sim over varying variables of interest

simwrapper <- function(nobs = 200, nrep = 100, seed = 7,
                       delta = seq(from = 1, to = 5, by = 2),
                       xi = c(0, 0.8), vary_beta = c("all", "beta0", "beta1"),
                       binary_regressor = c(TRUE, FALSE),
                       binary_beta = c(TRUE, FALSE),
                       nrsteps = 1, shift = 0, nrlevels = 2,
                       nrsplitvar = 1,
                       only_intercept = c(TRUE, FALSE),
                       compare_pruning = FALSE,
                       test = c("ctree", "mfluc", "ctree_max", 
                                "ctree_cat", "ctree_max_cat", "mfluc_cat",
                                "ctree_bin", "ctree_max_bin", "mfluc_bin",
                                "ctree_cat_bin", "ctree_max_cat_bin", "mfluc_cat_bin",
                                "ctree_resid", "mfluc_resid", "ctree_resid_max", 
                                "ctree_resid_cat", "ctree_resid_max_cat", "mfluc_resid_cat",
                                "ctree_resid_bin", "ctree_resid_max_bin", "mfluc_resid_bin",
                                "ctree_resid_cat_bin", "ctree_resid_max_cat_bin", "mfluc_resid_cat_bin",
                                "guide_sum_12", "guide_max_12", "guide_coin_12",
                                "guide_sum_1", "guide_max_1", "guide_coin_1",
                                "guide_sum_2", "guide_max_2", "guide_coin_2",
                                "guide_sum_1_cor", "guide_max_1_cor", "guide_coin_1_cor",
                                "guide_sum_2_cor", "guide_max_2_cor", "guide_coin_2_cor"),
                       beta0=0, beta1 = 1,
                       stump = TRUE, z1dist = "unif", sigma = 1, alpha = 0.05)  # only needed for dgp_stump
{
  
  prs <- expand.grid(delta = delta, xi = xi, vary_beta = vary_beta,
                     binary_regressor = binary_regressor,
                     binary_beta = binary_beta,
                     only_intercept = only_intercept)
  
  rmid <- which((prs$only_intercept == TRUE & prs$vary_beta != "beta0"))
  if(length(rmid) > 0) prs <- prs[-rmid,]
  rmid <- which(prs$binary_beta == FALSE & prs$xi != 0)
  if(length(rmid) > 0) prs <- prs[-rmid,]
  rmid <- which(prs$only_intercept == TRUE & prs$binary_regressor == FALSE)
  if(length(rmid) > 0) prs <- prs[-rmid,]
  
  rownames(prs) <- c(1:NROW(prs))
  
  nprs <- nrow(prs)
  ntest <- length(test)
  
  if(stump){
    
    prop_nosplit <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    prop_split <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    prop_Fsplit <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    prop_Tsplit <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    prop_z1 <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    pval_min <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    pval_true <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    
  } else {
    nrsubgr <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    if(compare_pruning) nrsubgr_p <- matrix(rep(NA, ntest * nprs), ncol = ntest)
  }
  err_coef <- matrix(rep(NA, ntest * nprs), ncol = ntest)
  err_total <- matrix(rep(NA, ntest * nprs), ncol = ntest)
  ari <- matrix(rep(NA, ntest * nprs), ncol = ntest)
  
  if(compare_pruning){
    err_coef_p <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    err_total_p <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    ari_p <- matrix(rep(NA, ntest * nprs), ncol = ntest)
  }  
  
  
  for(i in 1:nprs) {
    reslist <- sim(nobs = nobs, nrep = nrep, stump = stump, seed = seed,
                   formula = y~x|z1+z2+z3+z4+z5+z6+z7+z8+z9+z10,
                   beta0 = beta0, beta1 = beta1, z1dist = z1dist,
                   sigma = sigma, alpha = alpha,
                   test = test, nrsteps = nrsteps,
                   nrsplitvar = nrsplitvar,
                   nrlevels = nrlevels, shift = shift,
                   vary_beta = prs$vary_beta[i], 
                   binary_regressor = prs$binary_regressor[i],
                   binary_beta = prs$binary_beta[i],
                   xi = prs$xi[i],
                   only_intercept = prs$only_intercept[i],
                   delta = prs$delta[i],
                   compare_pruning = compare_pruning)
    
    if(stump){
      prop_nosplit[i,] <- reslist$prop_nosplit
      prop_split[i,] <- reslist$prop_split
      prop_Fsplit[i,] <- reslist$prop_Fsplit
      prop_Tsplit[i,] <- reslist$prop_Tsplit
      prop_z1[i,] <- reslist$prop_z1
      pval_min[i,] <- reslist$pval_min
      pval_true[i,] <- reslist$pval_true
    }
    else {
      nrsubgr[i,] <- reslist$nrsubgr
      if(compare_pruning) nrsubgr_p[i,] <- reslist$nrsubgr_p
    }
    
    err_coef[i,] <- reslist$err_coef
    err_total[i,] <- reslist$err_total
    ari[i,] <- reslist$ari
    
    if(compare_pruning){
      err_coef_p[i,] <- reslist$err_coef_p
      err_total_p[i,] <- reslist$err_total_p
      ari_p[i,] <- reslist$ari_p
    }
  }
  
  rval <- data.frame()
  for(i in 1:ntest) rval <- rbind(rval, prs)
  rval$test <- gl(ntest, nprs, labels = test)
  if(stump){
    rval$prop_nosplit <- as.vector(prop_nosplit)
    rval$prop_split <- as.vector(prop_split)
    rval$prop_Fsplit <- as.vector(prop_Fsplit)
    rval$prop_Tsplit <- as.vector(prop_Tsplit)
    rval$prop_z1 <- as.vector(prop_z1)
    rval$pval_min <- as.vector(pval_min)
    rval$pval_true <- as.vector(pval_true)
  } else {
    rval$nrsubgr <- as.vector(nrsubgr)
    if(compare_pruning) rval$nrsubgr_p <- as.vector(nrsubgr_p)
  }
  rval$err_coef <- as.vector(err_coef)
  rval$err_total <- as.vector(err_total)
  rval$ari <- as.vector(ari)
  
  if(compare_pruning){
    rval$err_coef_p <- as.vector(err_coef_p)
    rval$err_total_p <- as.vector(err_total_p)
    rval$ari_p <- as.vector(ari_p)
  }
  
  rval$delta <- factor(rval$delta)
  rval$vary_beta <- factor(rval$vary_beta)
  rval$binary_regressor <- factor(rval$binary_regressor)
  rval$binary_beta <- factor(rval$binary_beta)
  rval$xi <- factor(rval$xi)
  rval$only_intercept <- factor(rval$only_intercept)
  
  return(rval)
}





simwrapper_p <- function(nobs = 200, nrep = 100, seed = 7, nrsteps = 1, nrlevels = 2, shift = 0,
                         delta = seq(from = 1, to = 5, by = 2),
                         xi = c(0, 0.8), vary_beta = c("all", "beta0", "beta1"),
                         binary_regressor = c(TRUE, FALSE),
                         binary_beta = c(TRUE, FALSE),
                         only_intercept = c(TRUE, FALSE),
                         compare_pruning = FALSE,
                         test = c("ctree", "mfluc", "ctree_max", 
                                  "ctree_cat", "ctree_max_cat", "mfluc_cat",
                                  "ctree_bin", "ctree_max_bin", "mfluc_bin",
                                  "ctree_cat_bin", "ctree_max_cat_bin", "mfluc_cat_bin",
                                  "ctree_resid", "mfluc_resid", "ctree_resid_max", 
                                  "ctree_resid_cat", "ctree_resid_max_cat", "mfluc_resid_cat",
                                  "ctree_resid_bin", "ctree_resid_max_bin", "mfluc_resid_bin",
                                  "ctree_resid_cat_bin", "ctree_resid_max_cat_bin", "mfluc_resid_cat_bin",
                                  "guide_sum_12", "guide_max_12", "guide_coin_12",
                                  "guide_sum_1", "guide_max_1", "guide_coin_1",
                                  "guide_sum_2", "guide_max_2", "guide_coin_2",
                                  "guide_sum_1_cor", "guide_max_1_cor", "guide_coin_1_cor",
                                  "guide_sum_2_cor", "guide_max_2_cor", "guide_coin_2_cor"),
                         beta0 = 0, beta1 = 1,
                         stump = TRUE, z1dist = "unif", sigma = 1, alpha = 0.05,
                         nrsplitvar = 1)
{
  
  cl <- match.call()
  
  prs <- expand.grid(delta = delta, xi = xi, vary_beta = vary_beta,
                     binary_regressor = binary_regressor,
                     binary_beta = binary_beta,
                     only_intercept = only_intercept)
  
  rmid <- which(prs$only_intercept == TRUE & prs$vary_beta != "beta0")
  if(length(rmid) > 0) prs <- prs[-rmid,]
  rmid <- which(prs$binary_beta == FALSE & prs$xi != 0)
  if(length(rmid) > 0) prs <- prs[-rmid,]
  rmid <- which(prs$only_intercept == TRUE & prs$binary_regressor == FALSE)
  if(length(rmid) > 0) prs <- prs[-rmid,]
  
  rownames(prs) <- c(1:NROW(prs))
  
  nprs <- nrow(prs)
  ntest <- length(test)
  
  simres <- mclapply(1:nprs, 
                  function(i) {
                    reslist <- sim(nobs = nobs, nrep = nrep, stump = stump, seed = seed,
                                   formula = y~x|z1+z2+z3+z4+z5+z6+z7+z8+z9+z10,
                                   beta0 = beta0, beta1 = beta1, z1dist = z1dist,
                                   sigma = sigma, alpha = alpha,
                                   test = test, nrsteps = nrsteps,
                                   nrlevels = nrlevels, shift = shift,
                                   nrsplitvar = nrsplitvar,
                                   vary_beta = prs$vary_beta[i], 
                                   binary_regressor = prs$binary_regressor[i],
                                   binary_beta = prs$binary_beta[i],
                                   xi = prs$xi[i],
                                   only_intercept = prs$only_intercept[i],
                                   delta = prs$delta[i],
                                   compare_pruning = compare_pruning)
                    
                    
                    #prop_nosplit[i,] <- reslist$prop_nosplit
                    #prop_split[i,] <- reslist$prop_split
                    #prop_Fsplit[i,] <- reslist$prop_Fsplit
                    #prop_Tsplit[i,] <- reslist$prop_Tsplit
                    
                    return(reslist)
                  },
                  mc.cores = detectCores() - 1
  )
  
  
  if(stump){
    
    prop_nosplit <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    prop_split <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    prop_Fsplit <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    prop_Tsplit <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    prop_z1 <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    pval_min <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    pval_true <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    
    for(i in 1:nprs){
      prop_nosplit[i,] <- simres[[i]]$prop_nosplit
      prop_split[i,] <- simres[[i]]$prop_split
      prop_Fsplit[i,] <- simres[[i]]$prop_Fsplit
      prop_Tsplit[i,] <- simres[[i]]$prop_Tsplit
      prop_z1[i,] <- simres[[i]]$prop_z1
      pval_min[i,] <- simres[[i]]$pval_min
      pval_true[i,] <- simres[[i]]$pval_true
    }
    
  } else {
    nrsubgr <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    for(i in 1:nprs) nrsubgr[i,] <- simres[[i]]$nrsubgr
    
    if(compare_pruning){
      nrsubgr_p <- matrix(rep(NA, ntest * nprs), ncol = ntest)
      for(i in 1:nprs) nrsubgr_p[i,] <- simres[[i]]$nrsubgr_p
    }
  }
  
  err_coef <- matrix(rep(NA, ntest * nprs), ncol = ntest)
  err_total <- matrix(rep(NA, ntest * nprs), ncol = ntest)
  ari <- matrix(rep(NA, ntest * nprs), ncol = ntest)
  
  for(i in 1:nprs){
    err_coef[i,] <- simres[[i]]$err_coef
    err_total[i,] <- simres[[i]]$err_total
    ari[i,] <- simres[[i]]$ari
  }
  
  if(compare_pruning){
    err_coef_p <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    err_total_p <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    ari_p <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    
    for(i in 1:nprs){
      err_coef_p[i,] <- simres[[i]]$err_coef_p
      err_total_p[i,] <- simres[[i]]$err_total_p
      ari_p[i,] <- simres[[i]]$ari_p
    }
  }
  
  res <- data.frame()
  for(i in 1:ntest) res <- rbind(res, prs)
  res$test <- gl(ntest, nprs, labels = test)
  
  
  if(stump){
    res$prop_nosplit <- as.vector(prop_nosplit)
    res$prop_split <- as.vector(prop_split)
    res$prop_Fsplit <- as.vector(prop_Fsplit)
    res$prop_Tsplit <- as.vector(prop_Tsplit)
    res$prop_z1 <- as.vector(prop_z1)
    res$pval_min <- as.vector(pval_min)
    res$pval_true <- as.vector(pval_true)
  } else {
    res$nrsubgr <- as.vector(nrsubgr)
    if(compare_pruning) res$nrsubgr_p <- as.vector(nrsubgr_p)
  }
  
  res$err_coef <- as.vector(err_coef)
  res$err_total <- as.vector(err_total)
  res$ari <- as.vector(ari)
  
  if(compare_pruning){
    res$err_coef_p <- as.vector(err_coef_p)
    res$err_total_p <- as.vector(err_total_p)
    res$ari_p <- as.vector(ari_p)
  }
  
  res$delta <- factor(res$delta)
  res$vary_beta <- factor(res$vary_beta)
  res$binary_regressor <- factor(res$binary_regressor)
  res$binary_beta <- factor(res$binary_beta)
  res$xi <- factor(res$xi)
  res$only_intercept <- factor(res$only_intercept)
  
  return(list(res = res,
              call = cl))
}



### function to prepare data set for 3-way factorial analysis
prep_3way <- function(simres)
{
  pval<- simres$pval
  # prepare data.frame
  {
    pval_T <- pval
    pval_T[which(simres$sv !="z1")] <- 1
    
    head(pval_T)
    
    d <- rep(colnames(pval)[1], NROW(pval))                 
    for(i in 2:NCOL(pval)){
      d <- c(d, rep(colnames(pval)[i], NROW(pval)))
    }
    
    d <- cbind(d, 
               rep(0, length(d)), 
               rep(0, length(d)), 
               rep(0, length(d)))    
    
    colnames(d) <- c("strategy", "res_scores", "bin", "cat")
    
    ## set values for factor variables representing the three switches
    {
      d[d[,"strategy"] == "ctree", "res_scores"] <- "scores"
      d[d[,"strategy"] == "ctree", "bin"] <- "lin"
      d[d[,"strategy"] == "ctree", "cat"] <- "lin"
      
      d[d[,"strategy"] == "ctree_max", "res_scores"] <- "scores"
      d[d[,"strategy"] == "ctree_max", "bin"] <- "lin"
      d[d[,"strategy"] == "ctree_max", "cat"] <- "max"
      
      d[d[,"strategy"] == "mfluc", "res_scores"] <- "scores"
      d[d[,"strategy"] == "mfluc", "bin"] <- "lin"
      d[d[,"strategy"] == "mfluc", "cat"] <- "max"
      
      d[d[,"strategy"] == "ctree_cat", "res_scores"] <- "scores"
      d[d[,"strategy"] == "ctree_cat", "bin"] <- "lin"
      d[d[,"strategy"] == "ctree_cat", "cat"] <- "cat"
      
      d[d[,"strategy"] == "ctree_max_cat", "res_scores"] <- "scores"
      d[d[,"strategy"] == "ctree_max_cat", "bin"] <- "lin"
      d[d[,"strategy"] == "ctree_max_cat", "cat"] <- "cat"
      
      d[d[,"strategy"] == "mfluc_cat", "res_scores"] <- "scores"
      d[d[,"strategy"] == "mfluc_cat", "bin"] <- "lin"
      d[d[,"strategy"] == "mfluc_cat", "cat"] <- "cat"
      
      d[d[,"strategy"] == "ctree_bin", "res_scores"] <- "scores"
      d[d[,"strategy"] == "ctree_bin", "bin"] <- "bin"
      d[d[,"strategy"] == "ctree_bin", "cat"] <- "lin"
      
      d[d[,"strategy"] == "ctree_max_bin", "res_scores"] <- "scores"
      d[d[,"strategy"] == "ctree_max_bin", "bin"] <- "bin"
      d[d[,"strategy"] == "ctree_max_bin", "cat"] <- "max"
      
      d[d[,"strategy"] == "mfluc_bin", "res_scores"] <- "scores"
      d[d[,"strategy"] == "mfluc_bin", "bin"] <- "bin"
      d[d[,"strategy"] == "mfluc_bin", "cat"] <- "max"
      
      d[d[,"strategy"] == "ctree_cat_bin", "res_scores"] <- "scores"
      d[d[,"strategy"] == "ctree_cat_bin", "bin"] <- "bin"
      d[d[,"strategy"] == "ctree_cat_bin", "cat"] <- "cat"
      
      d[d[,"strategy"] == "ctree_max_cat_bin", "res_scores"] <- "scores"
      d[d[,"strategy"] == "ctree_max_cat_bin", "bin"] <- "bin"
      d[d[,"strategy"] == "ctree_max_cat_bin", "cat"] <- "cat"
      
      d[d[,"strategy"] == "mfluc_cat_bin", "res_scores"] <- "scores"
      d[d[,"strategy"] == "mfluc_cat_bin", "bin"] <- "bin"
      d[d[,"strategy"] == "mfluc_cat_bin", "cat"] <- "cat"
      
      
      
      d[d[,"strategy"] == "ctree_resid", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "ctree_resid", "bin"] <- "lin"
      d[d[,"strategy"] == "ctree_resid", "cat"] <- "lin"
      
      d[d[,"strategy"] == "ctree_resid_max", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "ctree_resid_max", "bin"] <- "lin"
      d[d[,"strategy"] == "ctree_resid_max", "cat"] <- "max"
      
      d[d[,"strategy"] == "mfluc_resid", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "mfluc_resid", "bin"] <- "lin"
      d[d[,"strategy"] == "mfluc_resid", "cat"] <- "max"
      
      d[d[,"strategy"] == "ctree_resid_cat", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "ctree_resid_cat", "bin"] <- "lin"
      d[d[,"strategy"] == "ctree_resid_cat", "cat"] <- "cat"
      
      d[d[,"strategy"] == "ctree_resid_max_cat", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "ctree_resid_max_cat", "bin"] <- "lin"
      d[d[,"strategy"] == "ctree_resid_max_cat", "cat"] <- "cat"
      
      d[d[,"strategy"] == "mfluc_resid_cat", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "mfluc_resid_cat", "bin"] <- "lin"
      d[d[,"strategy"] == "mfluc_resid_cat", "cat"] <- "cat"
      
      d[d[,"strategy"] == "ctree_resid_bin", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "ctree_resid_bin", "bin"] <- "bin"
      d[d[,"strategy"] == "ctree_resid_bin", "cat"] <- "lin"
      
      d[d[,"strategy"] == "ctree_resid_max_bin", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "ctree_resid_max_bin", "bin"] <- "bin"
      d[d[,"strategy"] == "ctree_resid_max_bin", "cat"] <- "max"
      
      d[d[,"strategy"] == "mfluc_resid_bin", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "mfluc_resid_bin", "bin"] <- "bin"
      d[d[,"strategy"] == "mfluc_resid_bin", "cat"] <- "max"
      
      d[d[,"strategy"] == "ctree_resid_cat_bin", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "ctree_resid_cat_bin", "bin"] <- "bin"
      d[d[,"strategy"] == "ctree_resid_cat_bin", "cat"] <- "cat"
      
      d[d[,"strategy"] == "ctree_resid_max_cat_bin", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "ctree_resid_max_cat_bin", "bin"] <- "bin"
      d[d[,"strategy"] == "ctree_resid_max_cat_bin", "cat"] <- "cat"
      
      d[d[,"strategy"] == "mfluc_resid_cat_bin", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "mfluc_resid_cat_bin", "bin"] <- "bin"
      d[d[,"strategy"] == "mfluc_resid_cat_bin", "cat"] <- "cat"
      
      
      ## GUIDE
      
      d[d[,"strategy"] == "guide_sum_1", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "guide_sum_1", "bin"] <- "bin"
      d[d[,"strategy"] == "guide_sum_1", "cat"] <- "cat"
      
      d[d[,"strategy"] == "guide_coin_1", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "guide_coin_1", "bin"] <- "bin"
      d[d[,"strategy"] == "guide_coin_1", "cat"] <- "cat"
      
      d[d[,"strategy"] == "guide_max_1", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "guide_max_1", "bin"] <- "bin"
      d[d[,"strategy"] == "guide_max_1", "cat"] <- "cat"
      
      d[d[,"strategy"] == "guide_sum_1_cor", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "guide_sum_1_cor", "bin"] <- "bin"
      d[d[,"strategy"] == "guide_sum_1_cor", "cat"] <- "cat"
      
      d[d[,"strategy"] == "guide_coin_1_cor", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "guide_coin_1_cor", "bin"] <- "bin"
      d[d[,"strategy"] == "guide_coin_1_cor", "cat"] <- "cat"
      
      d[d[,"strategy"] == "guide_max_1_cor", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "guide_max_1_cor", "bin"] <- "bin"
      d[d[,"strategy"] == "guide_max_1_cor", "cat"] <- "cat"
      
      
      d[d[,"strategy"] == "guide_sum_2", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "guide_sum_2", "bin"] <- "bin"
      d[d[,"strategy"] == "guide_sum_2", "cat"] <- "cat"
      
      d[d[,"strategy"] == "guide_coin_2", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "guide_coin_2", "bin"] <- "bin"
      d[d[,"strategy"] == "guide_coin_2", "cat"] <- "cat"
      
      d[d[,"strategy"] == "guide_max_2", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "guide_max_2", "bin"] <- "bin"
      d[d[,"strategy"] == "guide_max_2", "cat"] <- "cat"
      
      d[d[,"strategy"] == "guide_sum_2_cor", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "guide_sum_2_cor", "bin"] <- "bin"
      d[d[,"strategy"] == "guide_sum_2_cor", "cat"] <- "cat"
      
      d[d[,"strategy"] == "guide_coin_2_cor", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "guide_coin_2_cor", "bin"] <- "bin"
      d[d[,"strategy"] == "guide_coin_2_cor", "cat"] <- "cat"
      
      d[d[,"strategy"] == "guide_max_2_cor", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "guide_max_2_cor", "bin"] <- "bin"
      d[d[,"strategy"] == "guide_max_2_cor", "cat"] <- "cat"
      
      d[d[,"strategy"] == "guide_sum_12", "res_scores"] <- "scores"
      d[d[,"strategy"] == "guide_sum_12", "bin"] <- "bin"
      d[d[,"strategy"] == "guide_sum_12", "cat"] <- "cat"
      
      d[d[,"strategy"] == "guide_max_12", "res_scores"] <- "scores"
      d[d[,"strategy"] == "guide_max_12", "bin"] <- "bin"
      d[d[,"strategy"] == "guide_max_12", "cat"] <- "cat"
      
      d[d[,"strategy"] == "guide_coin_12", "res_scores"] <- "scores"
      d[d[,"strategy"] == "guide_coin_12", "bin"] <- "bin"
      d[d[,"strategy"] == "guide_coin_12", "cat"] <- "cat"
    }
    
    
    d <- data.frame(strategy = factor(d[,"strategy"]),
                    res_scores = factor(d[,"res_scores"]),
                    bin = factor(d[,"bin"]),
                    cat = factor(d[,"cat"]))
    
    d$pval <- as.vector(pval)
    d$pvalT <- as.vector(pval_T)
    d$pval_z1 <- as.vector(simres$pval_z1)
  }
  
  return(d)
}




### simulate stump
if(FALSE){
  simres <- simwrapper_p(nobs = 250, nrep = 100, nrsteps = 1,
                         nrsplitvar = 1,
                         delta = seq(from = 0, to = 1, by = 0.25),
                         xi = c(0, 0.2, 0.5, 0.8), 
                         vary_beta = c("all", "beta0", "beta1"),
                         binary_regressor = c(TRUE, FALSE),
                         #binary_regressor = FALSE,
                         binary_beta = c(TRUE, FALSE),
                         only_intercept = c(TRUE, FALSE),
                         #only_intercept = FALSE,
                         test = c("ctree", "mfluc", "ctree_max",
                                  "ctree_cat", "mfluc_cat", "ctree_max_cat",
                                  "ctree_bin", "mfluc_bin", "ctree_max_bin",
                                  "ctree_cat_bin", "mfluc_cat_bin", "ctree_max_cat_bin",
                                  "guide_sum_12", 
                                  "guide_coin_12",
                                  "guide_sum_1_cor"),
                         beta0 = 0, beta1 = 1,
                         stump = TRUE, z1dist = "unif", sigma = 1, alpha = 0.05)
  
  
  save(simres, file = "~/svn/partykit/pkg/partykit/inst/guideliketest/sim/simres20181120_stump.rda")
  # save(simres, file = "~/partykit/inst/guideliketest/sim/simres20180501_tree.rda")
}

### simulate tree
if(FALSE){
  
  simres <- simwrapper_p(nobs = 250, nrep = 100, seed = 7, 
                         nrsteps = 1, nrlevels = 2,
                         delta = seq(from = 0, to = 1, by = 0.25),
                         xi = c(0, 0.2, 0.5, 0.8), 
                         vary_beta = "all",
                         #vary_beta = c("all", "beta0", "beta1"),
                         binary_regressor = FALSE,
                         binary_beta = TRUE,
                         only_intercept = FALSE,
                         test = c("ctree", "mfluc", "ctree_max",
                                  "ctree_cat", "mfluc_cat", "ctree_max_cat",
                                  "ctree_bin", "mfluc_bin", "ctree_max_bin",
                                  "ctree_cat_bin", "mfluc_cat_bin", "ctree_max_cat_bin",
                                  "guide_sum_12", 
                                  "guide_coin_12",
                                  "guide_sum_1_cor"),
                         beta0 = NULL, beta1 = NULL,
                         stump = FALSE, z1dist = "unif", sigma = 1, alpha = 0.05,
                         nrsplitvar = 2)
  
  save(simres, file = "~/svn/partykit/pkg/partykit/inst/guideliketest/sim/simres20181120_tree.rda")
  # save(simres, file = "~/partykit/inst/guideliketest/sim/simres20180501_tree.rda")
}



if(FALSE){
  subdata <- subset(simres, test %in% c("ctree", "ctree_cat", "ctree_bin", 
                                        "ctree_cat_bin", "guide_coin_12"))
  subdata <- subset(simres, test %in% c("mfluc", "mfluc_cat", "mfluc_bin", 
                                        "mfluc_cat_bin", "guide_coin_12"))
  subdata <- subset(simres, test %in% c("guide_coin_12", "guide_sum_12",
                                        "guide_sum_1_cor", "mfluc"))
  subdata$test <- factor(subdata$test)
  
  xyplot(nrsubgr ~ delta | vary_beta + xi, groups = ~ test, 
         data = subdata,
         subset = binary_beta == TRUE & binary_regressor == FALSE & 
           only_intercept == FALSE, 
         type = "b", auto.key = TRUE, ylim = c(0,5), panel=function(...) {
           panel.xyplot(...)
           panel.abline(h=3)
         })
  
  xyplot(err_coef ~ delta | vary_beta + xi, groups = ~ test, 
         data = subdata,
         subset = binary_beta == TRUE & binary_regressor == FALSE & 
           only_intercept == FALSE, 
         type = "b", auto.key = TRUE, ylim = c(0,1))
  
  xyplot(err_total ~ delta | vary_beta + xi, groups = ~ test, 
         data = subdata,
         subset = binary_beta == TRUE & binary_regressor == FALSE & 
           only_intercept == FALSE, 
         type = "b", auto.key = TRUE, ylim = c(0,1))
  
  xyplot(ari ~ delta | vary_beta + xi, groups = ~ test, 
         data = subdata,
         subset = binary_beta == TRUE & binary_regressor == FALSE & 
           only_intercept == FALSE, 
         type = "b", auto.key = TRUE, ylim = c(0,1))
}





if(FALSE){
  library("lattice")
  load("~/svn/partykit/pkg/partykit/inst/guideliketest/sim/simres20180224_2step.rda")
  
  
  subdata <- subset(simres, test %in% c("mfluc", "mfluc_cat", "mfluc_bin", "mfluc_cat_bin", "guide_coin_12"))
  subdata <- subset(simres, test %in% c("ctree", "ctree_cat", "ctree_bin", "ctree_cat_bin", "guide_coin_12"))
  subdata <- subset(simres, test %in% c("ctree", "ctree_cat", "ctree_bin", "ctree_cat_bin", "guide_coin_12", 
                                        "mfluc_cat", "mfluc_bin", "mfluc_cat_bin", "guide_coin_12"))
  subdata <- subset(simres, test %in% c("ctree_cat_bin", "mfluc_cat_bin", "guide_coin_12"))
  subdata <- subset(simres, test %in% c("ctree_bin", "mfluc_bin", "guide_coin_12"))
  subdata <- subset(simres, test %in% c("ctree_cat", "mfluc_cat", "guide_coin_12"))
  subdata$test <- factor(subdata$test)
  xyplot(prop_Tsplit ~ delta | xi + vary_beta + binary_beta, groups = ~ test, 
         data = subdata,
         type = "b", auto.key = TRUE, 
         subset = (binary_regressor == TRUE & only_intercept == FALSE))
  
  
  subdata <- subset(simres, test %in% c("ctree", "mfluc",
                                        "ctree_max", "ctree_max_cat",
                                        "guide_sum_12", 
                                        "guide_coin_12",
                                        "guide_sum_1_cor"))
  subdata$test <- factor(subdata$test)
  xyplot(prop_Tsplit ~ delta | xi + vary_beta + binary_beta, groups = ~ test, 
         data = subdata,
         type = "b", auto.key = TRUE, 
         subset = (binary_regressor == TRUE & only_intercept == FALSE) 
  )
  
  subdata <- subset(simres, test %in% c("ctree", "ctree_cat", "ctree_max",
                                        "guide_sum_12", 
                                        "guide_coin_12",
                                        "guide_sum_1_cor"))
  subdata$test <- factor(subdata$test)
  xyplot(prop_Tsplit ~ delta | xi + vary_beta + binary_beta, groups = ~ test, 
         data = subdata,
         type = "b", auto.key = TRUE, 
         subset = (binary_regressor == TRUE & only_intercept == FALSE) 
  )
  
  
  
  xyplot(prop_split ~ vary_beta | xi + binary_beta + delta, groups = ~ test, data = simres, type = "b", auto.key = TRUE)
  xyplot(prop_split ~ delta | xi + vary_beta + binary_beta + binary_regressor + only_intercept, groups = ~ test, data = simres, type = "b", auto.key = TRUE)
  
  
  # models considering all scores
  s_allscores <- subset(simres, test %in% c("ctree", "mfluc", "ctree_cat", "mfluc_cat", 
                                            "guide_sum_12", "guide_max_12", "guide_coin_12"))
  s_allscores$test <- factor(s_allscores$test)
  xyplot(prop_split ~ xi | vary_beta + binary_beta + delta, groups = ~ test, data = s_allscores, type = "b", auto.key = TRUE)
  xyplot(prop_split ~ delta | xi + vary_beta + binary_beta, groups = ~ test, data = s_allscores, type = "b", auto.key = TRUE)
  
  # ctree vs ctree_max
  s_ctree <- subset(simres, test %in% c("ctree", "ctree_max", "ctree_cat", "ctree_max_cat"))
  s_ctree$test <- factor(s_ctree$test)
  xyplot(prop_split ~ xi | vary_beta + binary_beta + delta, groups = ~ test, data = s_ctree, type = "b", auto.key = TRUE)
  xyplot(prop_split ~ delta | xi + vary_beta + binary_beta, groups = ~ test, data = s_ctree, type = "b", auto.key = TRUE)
  xyplot(prop_split ~ delta | xi + vary_beta + binary_beta, groups = ~ test, 
         data = s_ctree, type = "b", auto.key = TRUE, 
         subset = binary_regressor == TRUE & only_intercept == FALSE)
  xyplot(prop_Tsplit ~ delta | xi + vary_beta + binary_beta, groups = ~ test, 
         data = s_ctree, type = "b", auto.key = TRUE, 
         subset = binary_regressor == TRUE & only_intercept == FALSE)
  
  # all GUIDE models
  s_guide <- subset(simres, test %in% c("guide_sum_12", "guide_max_12", "guide_coin_12",
                                        "guide_sum_1", "guide_max_1", "guide_coin_1",
                                        "guide_sum_2", "guide_max_2", "guide_coin_2",
                                        "guide_sum_1_cor", "guide_max_1_cor", "guide_coin_1_cor",
                                        "guide_sum_2_cor", "guide_max_2_cor", "guide_coin_2_cor"))
  s_guide$test <- factor(s_guide$test)
  xyplot(prop_split ~ xi | vary_beta + binary_beta + delta, groups = ~ test, data = s_guide, type = "b", auto.key = TRUE)
  xyplot(prop_split ~ delta | xi + vary_beta + binary_beta, groups = ~ test, data = s_guide, type = "b", auto.key = TRUE)
  
  # all GUIDE models sum
  s_guide <- subset(simres, test %in% c("guide_sum_12",
                                        "guide_sum_1", 
                                        "guide_sum_2", 
                                        "guide_sum_1_cor", 
                                        "guide_sum_2_cor"))
  s_guide$test <- factor(s_guide$test)
  xyplot(prop_split ~ xi | vary_beta + binary_beta + delta, groups = ~ test, data = s_guide, type = "b", auto.key = TRUE)
  xyplot(prop_split ~ delta | xi + vary_beta + binary_beta, groups = ~ test, data = s_guide, type = "b", auto.key = TRUE)
  
  
  # GUIDE_12 models
  s_guide12 <- subset(simres, test %in% c("guide_sum_12", "guide_max_12", "guide_coin_12"))
  s_guide12$test <- factor(s_guide12$test)
  xyplot(prop_split ~ xi | vary_beta + binary_beta + delta, groups = ~ test, data = s_guide12, type = "b", auto.key = TRUE)
  xyplot(prop_split ~ delta | xi + vary_beta + binary_beta, groups = ~ test, data = s_guide12, type = "b", auto.key = TRUE)
  
  
  # GUIDE_1 models
  s_guide1 <- subset(simres, test %in% c("guide_sum_1", "guide_max_1", "guide_coin_1",
                                         "guide_sum_1_cor", "guide_max_1_cor", "guide_coin_1_cor"))
  s_guide1$test <- factor(s_guide1$test)
  xyplot(prop_split ~ xi | vary_beta + binary_beta + delta, groups = ~ test, data = s_guide1, type = "b", auto.key = TRUE)
  xyplot(prop_split ~ delta | xi + vary_beta + binary_beta, groups = ~ test, data = s_guide1, type = "b", auto.key = TRUE)
  
  
  # GUIDE_2 models
  s_guide2 <- subset(simres, test %in% c("guide_sum_2", "guide_max_2", "guide_coin_2",
                                         "guide_sum_2_cor", "guide_max_2_cor", "guide_coin_2_cor"))
  s_guide2$test <- factor(s_guide2$test)
  xyplot(prop_split ~ xi | vary_beta + binary_beta + delta, groups = ~ test, data = s_guide2, type = "b", auto.key = TRUE)
  xyplot(prop_split ~ delta | xi + vary_beta + binary_beta, groups = ~ test, data = s_guide2, type = "b", auto.key = TRUE)
  
  
  # ctree, mfluc, GUIDE_12
  s_scores2 <- subset(simres, test %in% c("ctree", "mfluc", "ctree_max", 
                                          "ctree_cat", "mfluc_cat", 
                                          "guide_sum_12", "guide_coin_12",
                                          "guide_sum_1_cor"))
  s_scores2$test <- factor(s_scores2$test)
  xyplot(prop_split ~ xi | vary_beta + binary_beta + delta, groups = ~ test, data = s_scores2, type = "b", auto.key = TRUE)
  xyplot(prop_split ~ delta | xi + vary_beta + binary_beta, groups = ~ test, data = s_scores2, type = "b", auto.key = TRUE)
  xyplot(prop_split ~ delta | xi + vary_beta + binary_beta, groups = ~ test, data = s_scores2, 
         subset = binary_regressor == TRUE & only_intercept == FALSE,
         type = "b", auto.key = TRUE)
  
  
  
  
  tab <- xtabs(prop_Tsplit ~ delta + xi + binary_beta + binary_regressor + test,
               data = simres)
  ftable(tab, row.vars = c("xi", "binary_beta", "binary_regressor", "delta"), 
         col.vars = "test")
}





if(FALSE){
  b1fix_contbeta_contreg <- simcomp(nobs = 100, nrep = 100, seed = 7,
                                    formula = y~x|z1+z2+z3+z4+z5+z6+z7+z8+z9+z10,
                                    z1dist = "unif", sigma = 1, only_intercept = FALSE,
                                    vary_beta = "beta0", beta1 = 1, xi = 0.0, delta = 2,
                                    stump = TRUE, binary_regressor = FALSE,
                                    binary_beta = FALSE)
  
  save(b1fix_contbeta_contreg, file = "~/svn/partykit/pkg/partykit/inst/guideliketest/sim/b1fix_contbeta_contreg.rda")
  
  
  
  load("~/svn/partykit/pkg/partykit/inst/guideliketest/sim/b1fix_contbeta_contreg.rda")
  simtest <- b1fix_contbeta_contreg
  
  # mean over all p-values (incl. p-values > 0.05)
  mean_all <- cbind(mean(simtest$pval_ctest),
                    mean(simtest$pval_mtest),
                    mean(simtest$pval_cctest),
                    mean(simtest$pval_mctest),
                    mean(simtest$pval_gstest12),
                    mean(simtest$pval_gmtest12),
                    mean(simtest$pval_gctest12),
                    mean(simtest$pval_gstest1),
                    mean(simtest$pval_gmtest1),
                    mean(simtest$pval_gctest1),
                    mean(simtest$pval_gstest2),
                    mean(simtest$pval_gmtest2),
                    mean(simtest$pval_gctest2))
  colnames(mean_all) <-  c("c","m", "cc", "mc", "gs12","gm12","gc12", "gs1","gm1","gc1", "gs2","gm2","gc2")
  mean_all
  
  # percentage of trees with correct split variable (incl. p-values > 0.05)
  prop_var <- cbind(sum(simtest$sv_ctest == "z1")/length(simtest$sv_ctest),
                    sum(simtest$sv_mtest == "z1")/length(simtest$sv_mtest),
                    sum(simtest$sv_cctest == "z1")/length(simtest$sv_cctest),
                    sum(simtest$sv_mctest == "z1")/length(simtest$sv_mctest),
                    sum(simtest$sv_gstest12 == "z1")/length(simtest$sv_gstest12),
                    sum(simtest$sv_gmtest12 == "z1")/length(simtest$sv_gmtest12),
                    sum(simtest$sv_gctest12 == "z1")/length(simtest$sv_gctest12),
                    sum(simtest$sv_gstest1 == "z1")/length(simtest$sv_gstest1),
                    sum(simtest$sv_gmtest1 == "z1")/length(simtest$sv_gmtest1),
                    sum(simtest$sv_gctest1 == "z1")/length(simtest$sv_gctest1),
                    sum(simtest$sv_gstest2 == "z1")/length(simtest$sv_gstest2),
                    sum(simtest$sv_gmtest2 == "z1")/length(simtest$sv_gmtest2),
                    sum(simtest$sv_gctest2 == "z1")/length(simtest$sv_gctest2))
  colnames(prop_var) <-  c("c","m", "cc", "mc", "gs12","gm12","gc12", "gs1","gm1","gc1", "gs2","gm2","gc2")
  barplot(prop_var, ylim = c(0,1))
  
  # percentage of trees with p-values <= 0.05
  prop_005 <- cbind(sum(simtest$pval_ctest <= 0.05)/length(simtest$sv_ctest),
                    sum(simtest$pval_mtest <= 0.05)/length(simtest$sv_mtest),
                    sum(simtest$pval_cctest <= 0.05)/length(simtest$sv_cctest),
                    sum(simtest$pval_mctest <= 0.05)/length(simtest$sv_mctest),
                    sum(simtest$pval_gstest12 <= 0.05)/length(simtest$sv_gstest12),
                    sum(simtest$pval_gmtest12 <= 0.05)/length(simtest$sv_gmtest12),
                    sum(simtest$pval_gctest12 <= 0.05)/length(simtest$sv_gctest12),
                    sum(simtest$pval_gstest1 <= 0.05)/length(simtest$sv_gstest1),
                    sum(simtest$pval_gmtest1 <= 0.05)/length(simtest$sv_gmtest1),
                    sum(simtest$pval_gctest1 <= 0.05)/length(simtest$sv_gctest1),
                    sum(simtest$pval_gstest2 <= 0.05)/length(simtest$sv_gstest2),
                    sum(simtest$pval_gmtest2 <= 0.05)/length(simtest$sv_gmtest2),
                    sum(simtest$pval_gctest2 <= 0.05)/length(simtest$sv_gctest2))
  colnames(prop_005) <-  c("c","m", "cc", "mc", "gs12","gm12","gc12", "gs1","gm1","gc1", "gs2","gm2","gc2")
  barplot(prop_005, ylim = c(0,1))
  
  # percentage of trees with correct split variable and p-values <= 0.05
  prop_var005 <- cbind(sum(simtest$sv_ctest == "z1" & simtest$pval_ctest <= 0.05)/length(simtest$sv_ctest),
                       sum(simtest$sv_mtest == "z1" & simtest$pval_mtest <= 0.05)/length(simtest$sv_mtest),
                       sum(simtest$sv_cctest == "z1" & simtest$pval_cctest <= 0.05)/length(simtest$sv_cctest),
                       sum(simtest$sv_mctest == "z1" & simtest$pval_mctest <= 0.05)/length(simtest$sv_mctest),
                       sum(simtest$sv_gstest12 == "z1" & simtest$pval_gstest12 <= 0.05)/length(simtest$sv_gstest12),
                       sum(simtest$sv_gmtest12 == "z1" & simtest$pval_gmtest12 <= 0.05)/length(simtest$sv_gmtest12),
                       sum(simtest$sv_gctest12 == "z1" & simtest$pval_gctest12 <= 0.05)/length(simtest$sv_gctest12),
                       sum(simtest$sv_gstest1 == "z1" & simtest$pval_gstest1 <= 0.05)/length(simtest$sv_gstest1),
                       sum(simtest$sv_gmtest1 == "z1" & simtest$pval_gmtest1 <= 0.05)/length(simtest$sv_gmtest1),
                       sum(simtest$sv_gctest1 == "z1" & simtest$pval_gctest1 <= 0.05)/length(simtest$sv_gctest1),
                       sum(simtest$sv_gstest2 == "z1" & simtest$pval_gstest2 <= 0.05)/length(simtest$sv_gstest2),
                       sum(simtest$sv_gmtest2 == "z1" & simtest$pval_gmtest2 <= 0.05)/length(simtest$sv_gmtest2),
                       sum(simtest$sv_gctest2 == "z1" & simtest$pval_gctest2 <= 0.05)/length(simtest$sv_gctest2))
  colnames(prop_var005) <-  c("c","m", "cc", "mc", "gs12","gm12","gc12", "gs1","gm1","gc1", "gs2","gm2","gc2")
  barplot(prop_var005, ylim = c(0,1))
  
  
  # boxplot of all p-values
  boxplot(simtest$pval_ctest, simtest$pval_mtest, 
          simtest$pval_cctest, simtest$pval_mctest, 
          simtest$pval_gstest12, simtest$pval_gmtest12, simtest$pval_gctest12,
          simtest$pval_gstest1, simtest$pval_gmtest1, simtest$pval_gctest1,
          simtest$pval_gstest2, simtest$pval_gmtest2, simtest$pval_gctest2,
          outline = FALSE, 
          names = c("c","m", "cc", "mc",
                    "gs12","gm12","gc12",
                    "gs1","gm1","gc1",
                    "gs2","gm2","gc2"))
  
  
  # boxplot of all p-values if correct split variable was found
  boxplot(simtest$pval_ctest[simtest$sv_ctest == "z1"], simtest$pval_mtest[simtest$sv_mtest == "z1"], 
          simtest$pval_cctest[simtest$sv_cctest == "z1"], simtest$pval_mctest[simtest$sv_mctest == "z1"], 
          simtest$pval_gstest12[simtest$sv_gstest12 == "z1"], simtest$pval_gmtest12[simtest$sv_gmtest12 == "z1"], 
          simtest$pval_gctest12[simtest$sv_gctest12 == "z1"],
          simtest$pval_gstest1[simtest$sv_gstest1 == "z1"], simtest$pval_gmtest1[simtest$sv_gmtest1 == "z1"], 
          simtest$pval_gctest1[simtest$sv_gctest1 == "z1"],
          simtest$pval_gstest2[simtest$sv_gstest2 == "z1"], simtest$pval_gmtest2[simtest$sv_gmtest2 == "z1"], 
          simtest$pval_gctest2[simtest$sv_gctest2 == "z1"],
          outline = FALSE, 
          names = c("c","m", "cc", "mc",
                    "gs12","gm12","gc12",
                    "gs1","gm1","gc1",
                    "gs2","gm2","gc2"))
  
  
  boxplot(simtest$pval_ctest, simtest$pval_mtest, 
          simtest$pval_cctest, simtest$pval_mctest, 
          simtest$pval_gstest12, simtest$pval_gmtest12, simtest$pval_gctest12,
          simtest$pval_gstest1, simtest$pval_gmtest1, simtest$pval_gctest1,
          #simtest$pval_gstest2, simtest$pval_gmtest2, simtest$pval_gctest2,
          outline = FALSE, 
          names = c("c", "m", "cc", "mc",
                    "gs12","gm12", "gc12",
                    "gs1","gm1","gc1")) #,
  #"gs2","gm2","gc2"))
  
  
}



