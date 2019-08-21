# -------------------------------------------------------------------
# - NAME:   circfit.R
# - AUTHOR: Moritz N. Lang, Lisa Schlosser
# - DATE:   2019-08-21
# -------------------------------------------------------------------
# - PURPOSE: Wrapper function for distfit plus S3 methods 
# -------------------------------------------------------------------

## Wrapper function for distfit

circfit <- function(y, weights = NULL, start = NULL, start.eta = NULL,
                    response_range = NULL,
                    vcov = TRUE, type.hessian =  c("checklist", "analytic", "numeric"), 
                    method = "L-BFGS-B", estfun = TRUE, optim.control = list(), ...){
  
  ## Get and modify call
  cl <- match.call()
  cl2 <- cl
  cl2[[1]] <- quote(disttree::distfit)
  cl2$family <- dist_vonmises()
  cl2$response_range <- NULL
  
  ## Transform data to parameter range (-pi, pi)
  y <- angle_trans(y, 
                   start = response_range[1], 
                   end = response_range[2])
  cl2$y <- y
  
  ## Evaluate call
  fit <- eval(cl2)
  fit$call <- cl
  fit$response_range <- attr(y, "response_range") 
  
  class(fit) <- c("circfit", class(fit))
  fit
}

