# -------------------------------------------------------------------
# - NAME:   circtree.R
# - AUTHOR: Moritz N. Lang, Lisa Schlosser
# - DATE:   2019-07-31
# -------------------------------------------------------------------
# - PURPOSE: Wrapper function for disttree plus S3 methods 
# -------------------------------------------------------------------

## Wrapper function for disttree
circtree <- function(formula,
                     data,
                     response_range = NULL, ## TODO: or default c(0,2*pi) with check of values and stop function
                     subset,
                     na.action = na.pass,
                     weights,
                     offset,
                     cluster,
                     control = disttree_control(...),
                     converged = NULL,
                     scores = NULL,
                     doFit = TRUE,
                     ...) {
  
  ## Get and modify call
  cl <- match.call()
  cl2 <- cl
  cl2[[1]] <- quote(disttree::disttree)
  cl2$family <- dist_vonmises()
  cl2$response_range <- NULL

  formula <- Formula::as.Formula(formula)
  if(length(formula)[2L] > 1L) {
    formula <- Formula::Formula(formula(formula, rhs = 2L))  
    warning("formula must not have more than one RHS parts (only partitioning variables allowed)")
  }
  response.name <- as.character(formula[[2]])

  ## Transform data to parameter range (-pi, pi)
  data[, response.name] <- angle_trans(data[, response.name], 
                                       start = response_range[1], 
                                       end = response_range[2])
  cl2$data <- data

  ## Evaluate call
  tree <- eval(cl2)
  tree$info$call <- cl
  tree$info$response_range <- attr(data[, response.name], "response_range") 
  
  class(tree) <- c("circtree", class(tree))
  tree
}


## Print method
print.circtree <- function(x, title = NULL, objfun = "negative log-likelihood", ...){
  if(is.null(title)) title <- sprintf("Circular regression tree (von Mises Distribution)")
  partykit::print.modelparty(x, title = title, objfun = objfun, ...)
}


## TODO: Does kappa stay the same for response scale??!
## Coef method
coef.circtree <- function(object, response_range = FALSE, ...){
  
  ## Check response range
  if(!(is.logical(response_range) | (is.numeric(response_range) & all(is.finite(response_range)) & 
    length(response_range) == 2))) 
    stop("argument 'range' has to be logical or a numeric vector of length 2")

  ## If TRUE, use response_range of object
  if(identical(response_range, TRUE)) response_range <- object$info$response_range 
  
  ## Call 
  cf <- partykit:::coef.modelparty(object)
 
  ## Transform if necessary 
  if(is.numeric(response_range)) {
    cf[, "mu"] <- angle_retrans(cf[, "mu"], response_range[1], response_range[2])
    cf[,"kappa"] <- cf[, "kappa"] * (2 * pi)^2 / diff(response_range)^2
  }
  
  return(cf)
}



## logLik method
logLik.circtree <- function(object, newdata = NULL, weights = NULL, ...) {

  ## Get call
  cl <- match.call()

  ## Transform newdata to same range as fitted parameters
  if(!is.null(newdata)){
    formula <- Formula::as.Formula(object$info$formula)
    if(length(formula)[2L] > 1L) {
      formula <- Formula::Formula(formula(formula, rhs = 2L))
      warning("formula must not have more than one RHS parts (only partitioning variables allowed)")
    }
    response.name <- as.character(formula[[2]])

    newdata[, response.name] <- angle_trans(newdata[, response.name],
                                            start = object$info$response_range[1],
                                            end = object$info$response_range[2])
    cl$newdata
  }

  ## Evaluate call
  cl[[1]] <- quote(disttree:::logLik.disttree)
  eval(cl)
}


## Predict method
predict.circtree <- function (object, newdata = NULL, type = c("parameter", "response", "node"), 
                              response_range = FALSE, OOB = FALSE, ...) {

  type <- match.arg(type) 

  ## Check response range
  if(!(is.logical(response_range) | (is.numeric(response_range) & all(is.finite(response_range)) &
    length(response_range) == 2)))
    stop("argument 'range' has to be logical or a numeric vector of length 2")

  ## If TRUE, use response_range of object
  if(identical(response_range, TRUE)) response_range <- object$info$response_range

  ## Call
  cl <- match.call()
  cl[[1]] <- quote(disttree:::predict.disttree)
  rval <- eval(cl)

  ## Transform if necessary
  if(type == "response" & is.numeric(response_range)){
    rval <- angle_retrans(rval,
                          start = response_range[1],
                          end = response_range[2])
  } else if(type == "parameter" & is.numeric(response_range)){
    rval[, "mu"] <- angle_retrans(rval[, "mu"], response_range[1], response_range[2])
    rval[,"kappa"] <- rval[, "kappa"] * (2 * pi)^2 / diff(response_range)^2
  }
  return(rval)
}


## Simulate data
circtree_simulate <- function(n = 1000, mu = c(0, 2, 5), kappa = c(3, 3, 1), 
  response_range = c(0, 2 * pi), seed = 111){
  # FIXME: Extend to more general cases
  set.seed(seed)
  d <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
  d$group <- ifelse(d$x1 < 0, 1, ifelse(d$x2 < 0, 2, 3))
  d$mu <- mu[d$group]
  d$kappa <- kappa[d$group]

  for(i in 1:n){
    d[i, "y"] <- circular::rvonmises(1, mu = circular::circular(d$mu[i]), kappa = d$kappa[i])
  }
  d$y <- angle_trans(d$y, 0, 2*pi)
  d$y <- angle_retrans(d$y, response_range[1], response_range[2])
  return(d)
}
