
#' Variable Importance for pmforest
#'
#' See \code{\link[partykit]{varimp.cforest}}.
#'
#' @param object DESCRIPTION.
#' @param nperm the number of permutations performed.
#' @param OOB a logical determining whether the importance is computed from the 
#' out-of-bag sample or the learning sample (not suggested).
#' @param risk the risk to be evaluated. By default the objective function 
#' (e.g. log-Likelihood) is used.
#' @param ... passed on to \code{\link[partykit]{varimp.cforest}}.
#'
#' @return A vector of ‘mean decrease in accuracy’ importance scores.
varimp.pmforest <- function(object, nperm = 1L, OOB = TRUE,
                            risk = function(x, ...) - objfun(x, sum = TRUE, ...), 
                             ...) {
  
  varimp.cforest(object = object, nperm = nperm, OOB = OOB, risk = risk, 
                 ...)
  
}
