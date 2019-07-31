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


## Coef method
coef.circtree <- function(object, ...){
  partykit:::coef.modelparty(object)
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
