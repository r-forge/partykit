# -------------------------------------------------------------------
# - NAME:   circforest.R
# - AUTHOR: Moritz N. Lang, Lisa Schlosser
# - DATE:   2019-08-04
# -------------------------------------------------------------------
# - PURPOSE: Wrapper function for distexforest plus S3 methods 
# -------------------------------------------------------------------

## Wrapper function for distexforest
circforest <- function(formula,
                       data,       
                       response_range = NULL, ## TODO: or default c(0,2*pi) with check of values and stop function
                       subset, 
                       na.action = na.pass,
                       weights,
                       offset, 
                       cluster,
                       strata,
                       control = distextree_control(
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
  
  ## Get and modify call
  cl <- match.call()
  cl2 <- cl
  cl2[[1]] <- quote(disttree::distexforest)
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
  forest <- eval(cl2)
  forest$info$call <- cl
  
  ## Retransform response 
  forest$fit$`(response)` <- angle_retrans(forest$fit$`(response)`, 
                                           start = attr(cl2$data[,response.name], "response_range")[1],
                                           end = attr(cl2$data[,response.name], "response_range")[2])
  # TODO: to be fixed in distextree, use only fit or fitted
  forest$fitted$`(response)` <- angle_retrans(forest$fitted$`(response)`, 
                                              start = attr(cl2$data[,response.name], "response_range")[1],
                                              end = attr(cl2$data[,response.name], "response_range")[2])
  
  class(forest) <- c("circforest", class(forest))
  forest
}



predict.circforest <- function(object, newdata = NULL, 
                               type = c("parameter", "response", "weights", "node"), 
                               OOB = TRUE, 
                               scale = TRUE, ...) {
  
  type <- match.arg(type) 
  
  cl <- match.call()
  cl[[1]] <- quote(disttree:::predict.distexforest)
  
  ## NOTE: predict.distexforest uses object$fitted[["(response)"]] as response values for parameter estimation,
  # therefore object$fitted[["(response)"]] has to be transformed to (-pi,pi] within predict.circforest
  if(type == "parameter" | type == "response"){
    object$fitted[["(response)"]] <- angle_trans(object$fitted[["(response)"]],
                                                 start = attr(object$fitted$`(response)`, "response_range")[1],
                                                 end = attr(object$fitted$`(response)`, "response_range")[2])
    cl$object <- object
  }  
  
  ## For 'type=response' use 'type=parameter' and transform parameter mu
  if(type != "response"){
    eval(cl)
  } else {
    cl$type <- "parameter"
    parameters <- eval(cl)
    formula <- if(is.name(object$info$call$formula)) eval(object$info$call$formula) else object$info$call$formula
    formula <- Formula::as.Formula(formula)
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


logLik.circforest <- function(object, newdata = NULL, weights = NULL, ...){
  
  ## Get call
  cl <- match.call()
  
  ## NOTE: logLik.distexforest calls predict.distexforest as a first step which 
  # uses object$fitted[["(response)"]] as response values for parameter estimation,
  # therefore object$fitted[["(response)"]] has to be transformed to (-pi,pi] within predict.circforest
  object$fitted[["(response)"]] <- angle_trans(object$fitted[["(response)"]],
                                               start = attr(object$fitted$`(response)`, "response_range")[1],
                                               end = attr(object$fitted$`(response)`, "response_range")[2])
  cl$object <- object
  
  
  
  # Get response name
  formula <- if(is.name(object$info$call$formula)) eval(object$info$call$formula) else object$info$call$formula
  formula <- Formula::as.Formula(formula)
  object$info$call$formula <- formula  ## TODO: necessary? or should eval(formula) be performed in logLikd.distexforest as well?
  cl$object <- object
  if(length(formula)[2L] > 1L) {
    formula <- Formula::Formula(formula(formula, rhs = 2L))
    warning("formula must not have more than one RHS parts (only partitioning variables allowed)")
  }
  response.name <- as.character(formula[[2]])
    
  
  ## Transform response (possibly from newdata) to same range as fitted parameters: 
  # for newdata: should we expect newdata to be on the same range as defined by response_range?
  # otherwise: logLik.distexforest uses object$fitted[["(response)"]] which has already 
  #            been transformed do (-pi,pi] within logLik.circforest
  if(!is.null(newdata)) {
    newdata[, response.name] <- angle_trans(newdata[, response.name]) ## TODO: guess start and end?
    cl$newdata <- newdata
  }


  ## Evaluate call
  cl[[1]] <- quote(disttree:::logLik.distexforest)
  eval(cl)
}