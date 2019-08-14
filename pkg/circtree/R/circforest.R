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
  
  class(forest) <- c("circforest", class(forest))
  forest
}



predict.circforest <- function(object, newdata = NULL, 
                               type = c("parameter", "response", "weights", "node"), 
                               OOB = TRUE, 
                               scale = TRUE, ...) {
  
  type <- match.arg(type) 
  
  cl <- match.call()
  cl[[1]] <- quote(disttree:::predict.distforest)
  
  ## NOTE: predict.distforest uses object$fitted[["(response)"]] as response values for parameter estimation,
  # therefore object$fitted[["(response)"]] has to be transformed to (-pi,pi] within predict.circforest
  if(type == "parameter" | type == "response"){
    object$fitted[["(response)"]] <- angle_trans(object$fitted[["(response)"]],
                                                 start = attr(object$fitted[["(response)"]], "response_range")[1],
                                                 end = attr(object$fitted[["(response)"]], "response_range")[2])
    cl$object <- object
  }  
  
  ## For 'type=response' transform to response_range
  if(type != "response"){
    eval(cl)
  } else {
    response <- eval(cl)
    rval <- angle_retrans(response,
                          start = attr(object$fitted[["(response)"]], "response_range")[1],
                          end = attr(object$fitted[["(response)"]], "response_range")[2])
    return(rval)
  }
}


logLik.circforest <- function(object, newdata = NULL, weights = NULL, ...){
  
  ## Get call
  cl <- match.call()
  
  ## NOTE: logLik.distforest calls predict.distforest as a first step which 
  # uses object$fitted[["(response)"]] as response values for parameter estimation,
  # therefore object$fitted[["(response)"]] has to be transformed to (-pi,pi] within predict.circforest
  object$fitted[["(response)"]] <- angle_trans(object$fitted[["(response)"]],
                                               start = attr(object$fitted[["(response)"]], "response_range")[1],
                                               end = attr(object$fitted[["(response)"]], "response_range")[2])
  cl$object <- object
  
  
  
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
                                            start = attr(object$fitted$`(response)`, "response_range")[1],
                                            end = attr(object$fitted$`(response)`, "response_range")[2])
                                            
    cl$newdata <- newdata
  }


  ## Evaluate call
  cl[[1]] <- quote(disttree:::logLik.distforest)
  eval(cl)
}
