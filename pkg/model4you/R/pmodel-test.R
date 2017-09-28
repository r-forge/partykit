#' Test if personalised models improve upon 
#' base model.
#'
#' This is a rudimentary test if there is 
#' heterogeneity in the model parameters.
#'
#' @param forest pmforest object.
#' @param pmodels pmodel.identity object (pmodel(..., fun = identity)).
#' @param data data.
#' @param B number of bootstrap samples.
#'
#' @return list where the first element is the p-value und the 
#' second element is a data.frame with all neccessary infos to
#' compute the p-value.
#' 
#' @examples
#' \dontrun{
#' set.seed(123)
#' n <- 160
#' trt <- factor(rep(0:1, each = n/2))
#' y <- 4 + (trt == 1) + rnorm(n)
#' z <- matrix(rnorm(n * 2), ncol = 2)
#' 
#' dat <- data.frame(y, trt, z)
#' 
#' mod <- lm(y ~ trt, data = dat)
#' 
#' ## Note that ntree should usually be higher
#' frst <- pmforest(mod, ntree = 20) 
#' pmods <- pmodel(frst, fun = identity)
#' 
#' ## Note that B should be at least 100
#' ## The low B is just for demonstration 
#' ## purposes.
#' tst <- test_heterogeneity(forest = frst, 
#'                           pmodels = pmods, 
#'                           model = mod,
#'                           B = 10) 
#' tst$pvalue
#' plot(tst)
#' }
#' 
#' @importFrom stats simulate
#' 
#' @export
test_heterogeneity <- function(forest, pmodels = NULL, data = NULL, B = 100) {
  
  ## simulate new y
  model <- forest$info$model
  ys_new <- simulate(model, nsim = B)
  
  ## get info from model
  modcall <- getCall(model)
  modformula <- Formula::as.Formula(eval(modcall$formula))
  yformula <- formula(modformula, lhs = 1, rhs = 0)
  ynam <- all.vars(yformula)
  if(is.null(data)) {
    dat <- .get_model_data(model)
  }
  
  ## compute log-Likelihoods under the Null-Hypothesis
  get_loglik_yn <- function(yn) {
    dat[ , ynam] <- yn
    m_new <- update(model, data = dat)
    frst_new <- update(object = forest, model = m_new, data = dat)
    pmods_new <- pmodel(frst_new, fun = identity)
    c(base_model = logLik(m_new), 
      personalised_models = logLik(pmods_new))
  }
  logliks_null <- sapply(ys_new, get_loglik_yn)
  
  ## compute log-Likelihoods for orignial data
  logliks <- c(base_model = logLik(model), 
               personalised_models = logLik(pmodels))
  
  ## put all likelihoods in a data.frame with the needed info
  lls <- rbind.data.frame(logliks, t(logliks_null))
  lls$sim <- grepl("sim", row.names(lls))
  lls$hypothesis <- c("H[1]", "H[0]")[lls$sim + 1]
  lls$diff <- lls$base_model - lls$personalised_models
  
  ## create the return object
  ret <- list(pvalue = sum(lls$diff[!lls$sim] > lls$diff[lls$sim]) / B,
              logLiks = lls)
  class(ret) <- "heterogeneity_test"
  return(ret)
}


#' @rdname test_heterogeneity
#'
#' @param x object of class heterogeneity_test.
#' @param ... ignored.
#'
#' @importFrom ggplot2 geom_rug geom_vline scale_linetype element_blank labs aes_string
#' 
#' @export
plot.heterogeneity_test <- function(x, ...) {
  lls <- x$logLiks
  sim <- lls[lls$sim, ]
  orig <- lls[!lls$sim, ]
  ggplot() + 
    geom_line(data = sim, aes(x = diff), stat = "density") +
    geom_rug(data = lls, aes_string(x = "diff", linetype = "hypothesis")) +
    geom_vline(data = orig, aes(xintercept = diff),
               linetype = 2) +
    scale_linetype(labels = expression(H[0], data)) + 
    theme(legend.title=element_blank()) +
    labs(x = "Difference in log-Likelihoods\n between base model and personalised models")
}











