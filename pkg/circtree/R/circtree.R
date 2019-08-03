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

  ## Transform data to parameter range (-pi, pi)
  if(!is.null(response_range) && length(response_range) != 2)
    stop("argument 'response_range' hast to be defined by 2 values 
      (start and end value of the circular interval)")
  
  formula <- Formula::as.Formula(formula)
  if(length(formula)[2L] > 1L) {
    formula <- Formula::Formula(formula(formula, rhs = 2L))  
    warning("formula must not have more than one RHS parts (only partitioning variables allowed)")
  }
  response.name <- as.character(formula[[2]])

  data[, response.name] <- angle_trans(data[, response.name], 
                                       start = response_range[1], 
                                       end = response_range[2])
  cl2$data <- data

  ## Evaluate call
  tree <- eval(cl2)
  tree$info$call <- cl

  ## Retransform response 
  tree$fit$`(response)` <- angle_retrans(tree$fit$`(response)`, 
                                         start = attr(cl2$data[,response.name], "response_range")[1],
                                         end = attr(cl2$data[,response.name], "response_range")[2])
  # TODO: to be fixed in distextree, use only fit or fitted
  tree$fitted$`(response)` <- angle_retrans(tree$fitted$`(response)`, 
                                         start = attr(cl2$data[,response.name], "response_range")[1],
                                         end = attr(cl2$data[,response.name], "response_range")[2])
  
  ## TODO: should we change the coefficients stored in tree structure accoring to response_range?
  #for(i in 1:length(tree)){
  #  tree[[i]]$node$info$coefficients["mu"] <- angle_retrans(tree[[i]]$node$info$coefficients["mu"],
  #                                                          start = attr(cl2$data[,response.name], "response_range")[1],
  #                                                          end = attr(cl2$data[,response.name], "response_range")[2])
  #}
  
  class(tree) <- c("circtree", class(tree))
  tree
}


## Print method
print.circtree <- function(x, title = NULL, objfun = "negative log-likelihood", ...){
  if(is.null(title)) title <- sprintf("Circular regression tree (von Mises Distribution)")
  partykit::print.modelparty(x, title = title, objfun = objfun, ...)
}


## Coef method
coef.circtree <- function(object, ...){
  #response_range <- attr(object$fitted$`(response)`, "response_range")
  cf <- partykit:::coef.modelparty(object)
  return(cf)
}


## logLik method
logLik.circtree <- function(object, newdata = NULL, weights = NULL, ...) {

  ## Get call
  cl <- match.call()

  ## Transform newdata to same range as fit
  if(!is.null(newdata)){
    formula <- Formula::as.Formula(object$info$call$formula)
    if(length(formula)[2L] > 1L) {
      formula <- Formula::Formula(formula(formula, rhs = 2L))
      warning("formula must not have more than one RHS parts (only partitioning variables allowed)")
    }
    response.name <- as.character(formula[[2]])

    newdata[, response.name] <- angle_trans(newdata[, response.name],
                                  start = attr(object$fitted$`(response)`, "response_range")[1],
                                  end = attr(object$fitted$`(response)`, "response_range")[2])
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

  ## Retransform 'type=response' to response range
  if(type != "response"){
    cl <- match.call()
    cl[[1]] <- quote(disttree:::predict.distextree)
    eval(cl)
  } else {
    cl <- match.call()
    cl[[1]] <- quote(disttree:::predict.distextree)
    cl$type <- "parameter"
    parameters <- eval(cl)
    formula <- Formula::as.Formula(object$info$call$formula)
    if(length(formula)[2L] > 1L) {
      formula <- Formula::Formula(formula(formula, rhs = 2L))
      warning("formula must not have more than one RHS parts (only partitioning variables allowed)")
    }
    response.name <- as.character(formula[[2]])

    rval <- angle_retrans(parameters$mu,
                        start = attr(object$fitted$`(response)`, "response_range")[1],
                        end = attr(object$fitted$`(response)`, "response_range")[2])
    return(rval)
  }
}


## Simulate data
circtree_simulate <- function(n = 1000, mu = c(0, 2, 5), kappa = c(3, 3, 1), seed = 111){
  # FIXME: Extend to more general cases
  set.seed(seed)
  d <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
  d$group <- ifelse(d$x1 < 0, 1, ifelse(d$x2 < 0, 2, 3))
  d$mu <- mu[d$group]
  d$kappa <- kappa[d$group]

  for(i in 1:n){
    d[i, "y"] <- circular::rvonmises(1, mu = circular::circular(d$mu[i]), kappa = d$kappa[i])
  }
  return(d)
}
