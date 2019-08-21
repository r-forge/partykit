# -------------------------------------------------------------------
# - NAME:   circforest.R
# - AUTHOR: Moritz N. Lang, Lisa Schlosser
# - DATE:   2019-08-04
# -------------------------------------------------------------------
# - PURPOSE: Wrapper function for distforest plus S3 methods 
# -------------------------------------------------------------------


## Wrapper function for distforest
circforest <- function(formula,
                       data,       
                       response_range = NULL, ## TODO: or default c(0,2*pi) with check of values and stop function
                       subset, 
                       na.action = na.pass,
                       weights,
                       offset, 
                       cluster,
                       strata,
                       control = disttree_control(
                         teststat = "quad", testtype = "Univ", mincriterion = 0,
                         saveinfo = FALSE, minsplit = 20, minbucket = 7, splittry = 2, ...),
                       ntree = 500L, 
                       fit.par = FALSE,
                       perturb = list(replace = FALSE, fraction = 0.632),
                       mtry = ceiling(sqrt(nvar)), 
                       applyfun = NULL,
                       cores = NULL, 
                       trace = FALSE,
                       ...) {

  ## To pass R CMD check: "no visible binding for global variable"
  nvar <- NULL

  ## Get original formula
  oformula <- as.formula(formula)
  
  ## Get and modify call
  cl <- match.call()
  cl2 <- cl
  cl2[[1]] <- quote(disttree::distforest)
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
  forest <- eval(cl2)
  forest$info$call <- cl
  forest$info$response_range <- attr(data[, response.name], "response_range") 
  
  class(forest) <- c("circforest", class(forest))
  forest
}



predict.circforest <- function(object, newdata = NULL, 
                               type = c("parameter", "response", "weights", "node"), 
                               OOB = TRUE, 
                               response_range = FALSE, ...) {
  
  ## Check response range
  if(!(is.logical(response_range) | (is.numeric(response_range) & all(is.finite(response_range)) & 
                                     length(response_range) == 2))) 
    stop("argument 'range' has to be logical or a numeric vector of length 2")
  
  ## If TRUE, use response_range of object
  if(identical(response_range, TRUE)) response_range <- object$info$response_range 
  
  type <- match.arg(type) 
  
  cl <- match.call()
  cl[[1]] <- quote(disttree:::predict.distforest)
  
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


logLik.circforest <- function(object, newdata = NULL, weights = NULL, ...){
  
  ## Get call
  cl <- match.call()
  
  # Get response name
  formula <- Formula::as.Formula(object$info$formula)
  if(length(formula)[2L] > 1L) {
    formula <- Formula::Formula(formula(formula, rhs = 2L))
    warning("formula must not have more than one RHS parts (only partitioning variables allowed)")
  }
  response.name <- as.character(formula[[2]])

  ## Transform response from newdata to same range as fitted parameters: 
  # for newdata: we expect newdata to be on the same range as defined by response_range
  if(!is.null(newdata)) {
    newdata[, response.name] <- angle_trans(newdata[, response.name],
                                            start = object$info$response_range[1],
                                            end = object$info$response_range[2])
                                            
    cl$newdata <- newdata
  }


  ## Evaluate call
  cl[[1]] <- quote(disttree:::logLik.distforest)
  eval(cl)
}


## varimp
varimp.circforest <- function(object, nperm = 1L, ...){

  cl <- match.call()

  # Get response name
  formula <- Formula::as.Formula(object$info$formula)
  if(length(formula)[2L] > 1L) {
    formula <- Formula::Formula(formula(formula, rhs = 2L))
    warning("formula must not have more than one RHS parts (only partitioning variables allowed)")
  }
  response.name <- as.character(formula[[2]])

  object$data[, response.name] <- angle_retrans(object$data[, response.name],
                                                start = object$info$response_range[1],
                                                end = object$info$response_range[2])

  ## Evaluate call
  cl$object <- object
  cl[[1]] <- quote(disttree:::varimp.distforest)
  eval(cl)
}


