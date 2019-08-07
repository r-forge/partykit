# -------------------------------------------------------------------
# - NAME:   circtree.R
# - AUTHOR: Moritz N. Lang, Lisa Schlosser
# - DATE:   2019-07-31
# -------------------------------------------------------------------
# - PURPOSE: Wrapper function for distextree plus S3 methods 
# -------------------------------------------------------------------

## Wrapper function for distextree
circtree <- function(formula,
                     data,
                     response_range = NULL, ## TODO: or default c(0,2*pi) with check of values and stop function
                     subset,
                     na.action = na.pass,
                     weights,
                     offset,
                     cluster,
                     control = disttree::distextree_control(...),
                     converged = NULL,
                     scores = NULL,
                     doFit = TRUE,
                     ...) {
  
  ## Get and modify call
  cl <- match.call()
  cl2 <- cl
  cl2[[1]] <- quote(disttree::distextree)
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
coef.circtree <- function(object, type = c("parameter", "response"), ...){
  type <- match.arg(type)

  cf <- partykit:::coef.modelparty(object)
  
  if(type == "response"){
    response_range <- attr(object$fitted[["(response)"]], "response_range")
    cf[, "mu"] <- angle_retrans(cf[, "mu"], response_range[1], response_range[2])
    cf[,"kappa"] <- cf[, "kappa"] * (2 * pi)^2 / diff(response_range)^2
  }

  return(cf)
}


## TODO: This S3-method should be defined for distextree, so no for circtree necessary anymorge
fitted.circtree <- function(object, ...){

  rval <- predict.circtree(object, newdata = NULL, type = "response", OOB = FALSE, ...)
  return(rval)

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
                                  start = attr(object$fitted[["(response)"]], "response_range")[1],
                                  end = attr(object$fitted[["(response)"]], "response_range")[2])
    cl$newdata
  }

  ## Evaluate call
  cl[[1]] <- quote(disttree:::logLik.distextree)
  eval(cl)
}


## Predict method
predict.circtree <- function (object, newdata = NULL, type = c("parameter", "response", "node"), 
                              OOB = FALSE, ...) {

  type <- match.arg(type) 

  ## For 'type=response' transform to response_range
  if(type != "response"){
    cl <- match.call()
    cl[[1]] <- quote(disttree:::predict.distextree)
    eval(cl)
  } else {
    cl <- match.call()
    cl[[1]] <- quote(disttree:::predict.distextree)
    response <- eval(cl)
    
    rval <- angle_retrans(response,
                          start = attr(object$fitted[["(response)"]], "response_range")[1],
                          end = attr(object$fitted[["(response)"]], "response_range")[2])
    return(rval)
  }
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
