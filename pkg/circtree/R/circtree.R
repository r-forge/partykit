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
                     circ_range = (-pi,pi), ## TODO: or default c(0,2*pi) with check of values and stop function
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
  cl <- match.call()
  cl2 <- cl
  cl2[[1]] <- quote(disttree::distextree)
  cl2$family <- dist_vonmises()
  cl2$circ_range <- NULL

  #if(!is.null(circ_range) && length(circ_range) != 2)
  #  stop("argument 'circ_range' hast to be defined by 2 values (start and end value of the circular interval)")
  
  formula <- Formula::as.Formula(formula)
  if(length(formula)[2L] > 1L) {
    formula <- Formula::Formula(formula(formula, rhs = 2L))  
    # NOTE: (LS) if rhs has more than 1 element it is here assumed that partitioning variables are handed over on 2nd slot
    warning("formula must not have more than one RHS parts (only partitioning variables allowed)")
  }
  response.name <- as.character(formula[[2]])
  
  if(any(data[,response.name]<circ_range[1] || data[,response.name]>circ_range[2]))
    stop("response variable has to be defined on the interval (-pi, pi]")
  if(circ_range != c(-pi,pi))
    stop("other response intervals than (-pi,pi] are currently not provided")
  
  #data[,response.name] <- angle_trans(data[,response.name], 
  #                                    start = circ_range[1], 
  #                                    end = circ_range[2])
  #cl2$data <- data

  tree <- eval(cl2)
  tree$info$call <- cl

  ## TODO: to be fixed in distextree, use only fit or fitted
  #tree$fit$`(response)` <- angle_retrans(tree$fit$`(response)`, 
  #                                       start = attr(cl2$data[,response.name], "circ_range")[1],
  #                                       end = attr(cl2$data[,response.name], "circ_range")[2])
  #tree$fitted$`(response)` <- angle_retrans(tree$fitted$`(response)`, 
  #                                       start = attr(cl2$data[,response.name], "circ_range")[1],
  #                                       end = attr(cl2$data[,response.name], "circ_range")[2])
  
  ## TODO: change coefficients stored in tree structure accoring to circ_range
  #for(i in 1:length(tree)){
  #  tree[[i]]$node$info$coefficients["mu"] <- angle_retrans(tree[[i]]$node$info$coefficients["mu"],
  #                                                          start = attr(cl2$data[,response.name], "circ_range")[1],
  #                                                          end = attr(cl2$data[,response.name], "circ_range")[2])
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
  #circ_range <- attr(object$fitted$`(response)`, "circ_range")
  cf <- partykit:::coef.modelparty(object)
  #cf[,"mu"] <- circtree:::angle_retrans(cf[,"mu"], start = circ_range[1], end = circ_range[2])
  return(cf)
}


## logLik method
logLik.circtree <- function(object, newdata = NULL, weights = NULL, ...) {
  cl <- match.call()
  cl[[1]] <- quote(disttree:::logLik.distextree)
  eval(cl)
}


## Predict method
predict.circtree <- function (object, newdata = NULL, type = c("parameter", "response", "node"), 
                              OOB = FALSE, ...) {
  cl <- match.call()
  cl[[1]] <- quote(disttree:::predict.distextree)
  eval(cl)
}


## Simulate data
circtree_simulate <- function(n = 1000, mu = c(0, 2, 1), kappa = c(3, 3, 1), seed = 111){
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
