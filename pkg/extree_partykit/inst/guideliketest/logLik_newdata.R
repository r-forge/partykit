# logLik method for lm with optional argument newdata
logLik_lm_nd <- function (object, REML = FALSE, newdata = NULL, ...) 
{
  if (inherits(object, "mlm")) 
    stop("'logLik.lm' does not support multiple responses")
  if(!is.null(newdata)){
    ## FIX ME: better way of choosing response variable in newdata
    respname <- as.character(object$terms[[2]])
    res <- predict(object = object, newdata = newdata) - newdata[,respname]
  } 
    else res <- object$residuals
  p <- object$rank
  N <- length(res)
  if (is.null(w <- object$weights)) {
    w <- rep.int(1, N)
  }
  else {
    excl <- w == 0
    if (any(excl)) {
      res <- res[!excl]
      N <- length(res)
      w <- w[!excl]
    }
  }
  N0 <- N
  if (REML) 
    N <- N - p
  val <- 0.5 * (sum(log(w)) - N * (log(2 * pi) + 1 - log(N) + 
                                     log(sum(w * res^2))))
  if (REML) 
    val <- val - sum(log(abs(diag(object$qr$qr)[1L:p])))
  attr(val, "nall") <- N0
  attr(val, "nobs") <- N
  attr(val, "df") <- p + 1
  class(val) <- "logLik"
  val
}

# logLik method for modelparty and lm in nodes with optional argument newdata
logLik_modelparty_lm_nd <- function (object, dfsplit = NULL, newdata = NULL, ...) 
{
  if (is.null(dfsplit)) 
    dfsplit <- object$info$control$dfsplit
  dfsplit <- as.integer(dfsplit)
  ids <- partykit::nodeids(object, terminal = TRUE)
  
  if(!is.null(newdata)){
    ids_newdata <- predict(object, newdata = newdata, type = "node")
    newdata$ids <- ids_newdata
    ll <- list()
    for(i in ids){
      if(i %in% ids_newdata){
        ll <- c(ll, partykit:::apply_to_models(object, node = i, 
                                               FUN = function(object, REML = FALSE, ...) 
                                                 logLik_lm_nd(object, REML, 
                                                              newdata = newdata[newdata$ids == i,], ...)))
      }
    }
  } else {
    ll <- apply_to_models(object, node = ids, FUN = logLik)
  }
  
  
  ## FIX ME: Should every split be counted as degree of freedom? (in ccprune this is already penalized by nr of terminal nodes)
  ##         Should df of the models in the terminal nodes be multiplied with nr of terminal nodes? (same as above)
  
  #dfmod <- attr(ll[[1]], "df")
  #if(length(dfsplit) == 0)
  #  dfsplit <- (dfmod * length(ids))  + (length(object) - length(ll))
  
  if(length(dfsplit) == 0) dfsplit <- 1
  dfsplit <- dfsplit * (length(object) - length(ll))
  structure(sum(as.numeric(ll)), df = sum(sapply(ll, function(x) attr(x, 
                                                                      "df"))) + dfsplit, nobs = nobs(object), class = "logLik")
}
