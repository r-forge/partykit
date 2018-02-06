library("partykit")
library("Formula")


###########
# DGP

dgp <- function(nobs = 100, delta = 1, xi = 0,
                sigma = 1, seed = 7, only_intercept = FALSE, 
                binary_regressor = TRUE, binary_beta = TRUE,
                vary_beta = c("all", "beta0", "beta1"),
                beta0 = NULL, beta1 = NULL,
                z1dist = c("norm", "unif")){
  
  if(length(vary_beta) > 1) vary_beta <- vary_beta[1]
  if(!vary_beta %in% c("all", "beta0", "beta1")) stop("vary_beta has to be one of the following options: 'all', 'beta0', 'beta1'")
  if(length(z1dist) > 1) z1dist <- z1dist[1]
  if(!z1dist %in% c("norm", "unif")) stop("z1dist has to be one of the following options: 'norm', 'unif'")
  
  
  set.seed(seed)
  
  if(z1dist == "norm") {
    z1 <- rnorm(nobs, 0, 1) 
    z2 <- runif(nobs,-1,1)
  } else {
    z1 <- runif(nobs,-1,1)
    z2 <- rnorm(nobs, 0, 1) 
  }
  z3 <- rnorm(nobs, 0, 1)
  z4 <- rnorm(nobs, 0, 1)
  z5 <- rnorm(nobs, 0, 1)
  z6 <- rnorm(nobs, 0, 1)
  z7 <- runif(nobs, -1, 1)
  z8 <- runif(nobs, -1, 1)
  z9 <- runif(nobs, -1, 1)
  z10 <- runif(nobs, -1, 1)
  
  #if(!only_intercept)  {
  x <- if(binary_regressor) (-1)^rbinom(nobs, 1, 0.5) else runif(nobs, min = -1, max = 1)
  #}
  
  e <- rnorm(nobs, 0, sigma)
  
  if(vary_beta == "all"){
    if(binary_beta){
      beta0 <- delta * (-1)^(z1<=xi)
      beta1 <- delta * (-1)^(z1<=xi) * (-1)   # opposite signs if both betas vary
    } else {
      beta0 <- delta * z1
      beta1 <- delta * z1 * (-1)   # opposite signs if both betas vary
    }
  }
  
  if(vary_beta == "beta0"){
    beta0 <- if(binary_beta) delta * (-1)^(z1<=xi) else delta * z1
    beta1 <- beta1
  }
  
  if(vary_beta == "beta1"){
    beta0 <- beta0
    beta1 <- if(binary_beta) delta * (-1)^(z1<=xi) else delta * z1
  }
  
  
  y <- if(only_intercept) beta0 + e else beta0 + beta1 * x + e
  
  d <- data.frame(y = y, x = x, 
                  z1 = z1, z2 = z2, z3 = z3, z4 = z4, z5 = z5, z6 = z6, z7 = z7, z8 = z8, z9 = z9, z10 = z10,
                  beta0 = beta0, beta1 = beta1)
  
  return(d)
}



## FIX ME: no decorrelation for length(guide_parm) = 1
evaltests <- function(formula, data, testfun = c("guide", "ctree", "mfluc", "ctree_cat", "mfluc_cat"), 
                      subset, weights, stump = TRUE,
                      guide_testtype = c("sum", "max", "coin"), 
                      decorrelate = "vcov",
                      guide_parm = NULL,
                      xgroups = NULL, ygroups = NULL){
  
  na.action = na.pass
  converged = NULL
  scores = NULL
  doFit = TRUE
  offset = NULL 
  cluster = NULL
  
  if(length(testfun) > 1) testfun <- testfun[1]
  if(!testfun %in% c("guide", "ctree", "mfluc", "ctree_cat", "mfluc_cat")) 
    stop("testfun hast to be one of the following options: 'guide', 'ctree', 'mfluc', 'ctree_cat', 'mfluc_cat'")
  
  
  if(testfun == "guide"){
    
    if(length(guide_testtype) > 1) guide_testtype <- guide_testtype[1]
    if(!guide_testtype %in% c("max", "sum", "coin")) stop("testfun hast to be one of the following options: 'max', 'sum', 'coin'")
    
    control <- partykit:::extree_control(criterion = "p.value",
                                         selectfun = .guide_select(guide_decorrelate = "none",    # if required, decorrelation already within ytrafo 
                                                                   guide_testtype = guide_testtype,
                                                                   guide_parm = guide_parm,
                                                                   xgroups = xgroups,
                                                                   ygroups = ygroups,
                                                                   interaction = FALSE),  
                                         splitfun = partykit:::.objfun_split(restart = TRUE,
                                                                             intersplit = TRUE),
                                         svselectfun = .guide_select(guide_decorrelate = "none",    # if required, decorrelation already within ytrafo 
                                                                     guide_testtype = guide_testtype,
                                                                     guide_parm = guide_parm,
                                                                     xgroups = xgroups,
                                                                     ygroups = ygroups,
                                                                     interaction = FALSE), 
                                         svsplitfun = partykit:::.objfun_split(restart = TRUE,
                                                                               intersplit = TRUE),
                                         logmincriterion = log(1-0.05),
                                         update = TRUE,
                                         stump = stump)
  }
  
  if(testfun == "ctree"){
    control <- partykit:::extree_control(criterion = "p.value",
                                         selectfun = partykit:::.ctree_select(teststat = "quadratic", splitstat = "quadratic", 
                                                                              splittest = FALSE, pargs = GenzBretz(),
                                                                              testtype = "MonteCarlo", nresample = 9999L, 
                                                                              tol = sqrt(.Machine$double.eps),
                                                                              intersplit = TRUE, MIA = FALSE),
                                         splitfun = partykit:::.objfun_split(restart = TRUE,
                                                                             intersplit = TRUE),
                                         svselectfun = partykit:::.ctree_select(teststat = "quadratic", splitstat = "quadratic", 
                                                                                splittest = FALSE, pargs = GenzBretz(),
                                                                                testtype = "MonteCarlo", nresample = 9999L, 
                                                                                tol = sqrt(.Machine$double.eps),
                                                                                intersplit = TRUE, MIA = FALSE), 
                                         svsplitfun = partykit:::.objfun_split(restart = TRUE,
                                                                               intersplit = TRUE),
                                         logmincriterion = log(1-0.05),
                                         update = TRUE,
                                         stump = stump)
  }
  
  if(testfun == "ctree_cat"){
    control <- partykit:::extree_control(criterion = "p.value",
                                         selectfun = .ctree_cat_select(teststat = "quadratic", splitstat = "quadratic", 
                                                                              splittest = FALSE, pargs = GenzBretz(),
                                                                              testtype = "MonteCarlo", nresample = 9999L, 
                                                                              tol = sqrt(.Machine$double.eps),
                                                                              intersplit = TRUE, MIA = FALSE,
                                                                              xgroups = xgroups),
                                         splitfun = partykit:::.objfun_split(restart = TRUE,
                                                                             intersplit = TRUE),
                                         svselectfun = .ctree_cat_select(teststat = "quadratic", splitstat = "quadratic", 
                                                                                splittest = FALSE, pargs = GenzBretz(),
                                                                                testtype = "MonteCarlo", nresample = 9999L, 
                                                                                tol = sqrt(.Machine$double.eps),
                                                                                intersplit = TRUE, MIA = FALSE,
                                                                                xgroups = xgroups), 
                                         svsplitfun = partykit:::.objfun_split(restart = TRUE,
                                                                               intersplit = TRUE),
                                         logmincriterion = log(1-0.05),
                                         update = TRUE,
                                         stump = stump)
  }
  
  if(testfun == "mfluc"){
    control <- partykit:::extree_control(criterion = "p.value",
                              selectfun = partykit:::.mfluc_select(breakties = FALSE, 
                                                                   intersplit = TRUE, parm = NULL, 
                                                                   dfsplit = TRUE, 
                                                                   restart = TRUE, model = TRUE, 
                                                                   vcov = "sandwich", 
                                                                   ordinal = "chisq", 
                                                                   ytype = "vector",
                                                                   nrep = 10000L, terminal = "object", 
                                                                   inner = "object", trim = 0.1),  
                              splitfun = partykit:::.objfun_split(restart = TRUE,
                                                                  intersplit = TRUE),
                              svselectfun = partykit:::.mfluc_select(breakties = FALSE, 
                                                                     intersplit = TRUE, parm = NULL, 
                                                                     dfsplit = TRUE, 
                                                                     restart = TRUE, model = TRUE, 
                                                                     vcov = "sandwich", 
                                                                     ordinal = "chisq", 
                                                                     ytype = "vector",
                                                                     nrep = 10000L, terminal = "object", 
                                                                     inner = "object", trim = 0.1), 
                              svsplitfun = partykit:::.objfun_split(restart = TRUE,
                                                                    intersplit = TRUE),
                              logmincriterion = log(1-0.05),
                              update = TRUE,
                              stump = stump)
  }
  
  if(testfun == "mfluc_cat"){
    control <- partykit:::extree_control(criterion = "p.value",
                                         selectfun = .mfluc_cat_select(breakties = FALSE, 
                                                                              intersplit = TRUE, parm = NULL, 
                                                                              dfsplit = TRUE, 
                                                                              restart = TRUE, model = TRUE, 
                                                                              vcov = "sandwich", 
                                                                              ordinal = "chisq", 
                                                                              ytype = "vector",
                                                                              nrep = 10000L, terminal = "object", 
                                                                              inner = "object", trim = 0.1,
                                                                              xgroups = xgroups),  
                                         splitfun = partykit:::.objfun_split(restart = TRUE,
                                                                             intersplit = TRUE),
                                         svselectfun = .mfluc_cat_select(breakties = FALSE, 
                                                                                intersplit = TRUE, parm = NULL, 
                                                                                dfsplit = TRUE, 
                                                                                restart = TRUE, model = TRUE, 
                                                                                vcov = "sandwich", 
                                                                                ordinal = "chisq", 
                                                                                ytype = "vector",
                                                                                nrep = 10000L, terminal = "object", 
                                                                                inner = "object", trim = 0.1,
                                                                                xgroups = xgroups), 
                                         svsplitfun = partykit:::.objfun_split(restart = TRUE,
                                                                               intersplit = TRUE),
                                         logmincriterion = log(1-0.05),
                                         update = TRUE,
                                         stump = stump)
  }
  
  
  control$inner <- "object"
  control$terminal <- "object"
  
  ## set up model.frame() call
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$yx <- "matrix"
  mf$nmax <- control$nmax
  mf$ytype <- control$ytype
  ## evaluate model.frame
  mf[[1L]] <- quote(partykit::extree_data)
  
  d <- eval(mf, parent.frame())
  subset <- partykit:::.start_subset(d)
  
  weights <- model.weights(model.frame(d))
  
  if (is.null(control$update)) control$update <- TRUE
  
  # wrapper function to use lm() as trafo
  ytrafo <- function(subset, weights, estfun = TRUE, object = TRUE, info = NULL) {
    
    lmformula <- if(is.null(d$terms$yx)) Formula(d$terms$all) else Formula(d$terms$yx)   ## FIX ME: only regressors for lm, drop splitting variables
    sdata <- d$data[subset,]
    
    ## FIX ME: error in lm if weights are handed over
    ## FIX ME: scores with or without weights?
    #subweights <- if(is.null(weights) || (length(weights)==0L)) weights else weights[subset]
    
    model <- lm(formula = lmformula, data = sdata) #, weights = subweights) error in lm if weights are handed over
    
    ## FIX ME: add argument 'decorrelate' to control
    #decorrelate <- "vcov"
    #decorrelate <- "opg"
    #decorrelate <- "none"
    
    if(estfun) {
      ef <- as.matrix(sandwich::estfun(model)) 
      
      
      #if(ctrl$decorrelate != "none") {
      if(decorrelate != "none") {
        n <- NROW(ef)
        ef <- ef/sqrt(n)
        
        vcov <- if(decorrelate == "vcov") {
          vcov(model) * n
        } else {
          solve(crossprod(ef))
        }
        
        root.matrix <- function(X) {
          if((ncol(X) == 1L)&&(nrow(X) == 1L)) return(sqrt(X)) else {
            X.eigen <- eigen(X, symmetric = TRUE)
            if(any(X.eigen$values < 0)) stop("Matrix is not positive semidefinite")
            sqomega <- sqrt(diag(X.eigen$values))
            V <- X.eigen$vectors
            return(V %*% sqomega %*% t(V))
          }
        }
        ef <- as.matrix(t(root.matrix(vcov) %*% t(ef)))
      }
      
      estfun <- matrix(0, ncol = ncol(ef), nrow = nrow(d$data)) 
      estfun[subset,] <- ef
      ## FIX ME: return matrix of size nobs times par or only the scores from those observations selected by subset?
      ## FIX ME: estfun of lm-object returns unwheighted scores?
      
    } else estfun <- NULL
    
    
    
    object <-  if(object) model else NULL
    
    ret <- list(estfun = estfun,
                unweighted = TRUE, # unweighted = TRUE would prevent estfun / w in extree_fit
                coefficients = coef(model),
                objfun = -logLik(model),  # optional function to be minimized 
                object = object,
                decorrelated = (decorrelate != "none"),    
                converged = !is.null(model)  # FIX ME: better check whether coefficients are returned ?
    )
    return(ret)
  }
  
  
  converged <- TRUE
  control$model <- TRUE
  
  update <- function(subset, weights, control, doFit = TRUE)
    extree_fit(data = d, trafo = ytrafo, converged = converged, partyvars = d$variables$z, 
               subset = subset, weights = weights, ctrl = control, doFit = doFit)
  if (!doFit) return(list(d = d, update = update))
  
  tree <- update(subset = subset, weights = weights, control = control)
  trafo <- tree$trafo
  
  ### prepare as modelparty
  mf <- model.frame(d)
  if (length(weights) == 0) weights <- rep(1, nrow(mf))
  
  fitted <- data.frame("(fitted)" = fitted_node(tree$nodes, mf),
                       "(weights)" = weights,
                       check.names = FALSE)
  
  fitted[[3]] <- y <- mf[, d$variables$y, drop = TRUE]
  names(fitted)[3] <- "(response)"
  
  control$ytype <- ifelse(is.vector(y), "vector", class(y))
  # x <- model.matrix(modelf, data = mmf)
  control$xtype <- "matrix" # TODO: find out when to use data.frame
  
  ## return party object
  rval <- party(tree$nodes, 
                data = if(control$model) mf else mf[0,],
                fitted = fitted,
                terms = d$terms$all,
                info = list(
                  call = match.call(),
                  formula = formula,
                  Formula = as.Formula(formula),
                  terms = list(response = d$terms$yx, partitioning = d$terms$z),
                  fit = ytrafo,
                  control = control,
                  dots = list(),
                  nreg = NCOL(d$yx$x)
                )
  )
  class(rval) <- c("modelparty", class(rval))
  
  ### add modelinfo (object) and estfun if not there yet, but wanted
  # TODO: check if this can be done prettier
  which_terminals <- nodeids(rval, terminal = TRUE)
  which_all <- nodeids(rval)
  
  idx <- lapply(which_all, partykit:::.get_path, obj = tree$nodes)
  names(idx) <- which_all
  tree_ret <- unclass(rval)
  subset_term <- predict(rval, type = "node")
  
  for (i in which_all) {
    ichar <- as.character(i)
    iinfo <- tree_ret[[c(1, idx[[ichar]])]]$info
    
    if (i %in% which_terminals) winfo <- "object" else winfo <- control$inner
    #if (i %in% which_terminals) winfo <- control$terminal else winfo <- control$inner
    
    if (is.null(winfo)) {
      iinfo$object <- NULL
      iinfo$estfun <- NULL
    } else {
      if (is.null(iinfo) | any(is.null(iinfo[[winfo]])) | 
          any(! winfo %in% names(iinfo))) {
        iinfo <- trafo(subset = which(subset_term == i), weights = weights, info = NULL,
                       estfun = ("estfun" %in% winfo),
                       object = ("object" %in% winfo))
      }
    }
    
    tree_ret[[c(1, idx[[ichar]])]]$info <- iinfo
  }
  
  class(tree_ret) <- class(rval)
  
  return(tree_ret)
  
  if(FALSE){
    tree <- tree$nodes
    
    mf <- model.frame(d)
    if (is.null(weights)) weights <- rep(1, nrow(mf))
    
    fitted <- data.frame("(fitted)" = fitted_node(tree, mf), 
                         "(weights)" = weights,
                         check.names = FALSE)
    fitted[[3]] <- mf[, d$variables$y, drop = TRUE]
    names(fitted)[3] <- "(response)"
    ret <- party(tree, data = mf, fitted = fitted, 
                 info = list(call = match.call(), control = control))
    ret$update <- update
    ret$trafo <- trafo
    class(ret) <- c("constparty", class(ret))
    
    ### doesn't work for Surv objects
    # ret$terms <- terms(formula, data = mf)
    ret$terms <- d$terms$all
    ### need to adjust print and plot methods
    ### for multivariate responses
    ### if (length(response) > 1) class(ret) <- "party"
    return(ret)
  }
  
}


############
# compare on one data set
if(FALSE){
d <- dgp(400, vary_beta = "beta1", beta0 = 3)
d <- dgp(400, vary_beta = "all", xi = -0.5)
d <- dgp(1000, vary_beta = "all", xi = -0.5, binary_regressor = FALSE)
d <- dgp(1000, vary_beta = "all", xi = -0.0, binary_regressor = FALSE)

d <- dgp(400, vary_beta = "all", xi = -0.4, delta = 3, binary_regressor = TRUE)
d <- dgp(200, vary_beta = "beta1", beta0 = 0, xi = 0.0, binary_regressor = FALSE)
d <- dgp(200, vary_beta = "beta1", beta0 = 0, xi = -0.5, binary_regressor = FALSE)
d <- dgp(200, vary_beta = "beta1", beta0 = 0, xi = 0.0, binary_regressor = FALSE)

d <- dgp(200, vary_beta = "beta1", beta0 = 0, xi = 0.0, binary_regressor = FALSE, binary_beta = TRUE, delta = 5)
d <- dgp(200, vary_beta = "beta1", beta0 = 0, xi = 0.0, binary_regressor = FALSE, binary_beta = FALSE, delta = 5)

d <- dgp(200, vary_beta = "beta1", beta0 = 0, xi = 0.5, binary_regressor = FALSE, binary_beta = TRUE, delta = 5, z1dist = "unif")


ctest <- evaltests(y~x|z1+z2+z3, data = d, testfun = "ctree")
mtest <- evaltests(y~x|z1+z2+z3, data = d, testfun = "mfluc")
cctest <- evaltests(y~x|z1+z2+z3, data = d, testfun = "ctree_cat")
mctest <- evaltests(y~x|z1+z2+z3, data = d, testfun = "mfluc_cat")
gstest12 <- evaltests(y~x|z1+z2+z3, data = d, testfun = "guide", 
                    guide_testtype = "sum", guide_parm = c(1,2),
                    xgroups = NULL, ygroups = NULL)
gmtest12 <- evaltests(y~x|z1+z2+z3, data = d, testfun = "guide", 
                    guide_testtype = "max", guide_parm = c(1,2),
                    xgroups = NULL, ygroups = NULL)
gctest12 <- evaltests(y~x|z1+z2+z3, data = d, testfun = "guide", 
                    guide_testtype = "coin", guide_parm = c(1,2),
                    xgroups = NULL, ygroups = NULL)
gstest1 <- evaltests(y~x|z1+z2+z3, data = d, testfun = "guide", 
                      guide_testtype = "sum", guide_parm = c(1),
                      xgroups = NULL, ygroups = NULL)
gmtest1 <- evaltests(y~x|z1+z2+z3, data = d, testfun = "guide", 
                      guide_testtype = "max", guide_parm = c(1),
                      xgroups = NULL, ygroups = NULL)
gctest1 <- evaltests(y~x|z1+z2+z3, data = d, testfun = "guide", 
                      guide_testtype = "coin", guide_parm = c(1),
                      xgroups = NULL, ygroups = NULL)
gstest2 <- evaltests(y~x|z1+z2+z3, data = d, testfun = "guide", 
                      guide_testtype = "sum", guide_parm = c(2),
                      xgroups = NULL, ygroups = NULL)
gmtest2 <- evaltests(y~x|z1+z2+z3, data = d, testfun = "guide", 
                      guide_testtype = "max", guide_parm = c(2),
                      xgroups = NULL, ygroups = NULL)
gctest2 <- evaltests(y~x|z1+z2+z3, data = d, testfun = "guide", 
                      guide_testtype = "coin", guide_parm = c(2),
                      xgroups = NULL, ygroups = NULL)


plot(ctest)
plot(mtest)
plot(cctest)
plot(mctest)
plot(gmtest12)
plot(gstest12)
plot(gctest12)
plot(gmtest1)
plot(gstest1)
plot(gctest1)
plot(gmtest2)
plot(gstest2)
plot(gctest2)

plot(ctest, terminal_panel = node_bivplot)
plot(mtest, terminal_panel = node_bivplot)
plot(cctest, terminal_panel = node_bivplot)
plot(mctest, terminal_panel = node_bivplot)
plot(gmtest12, terminal_panel = node_bivplot)
plot(gstest12, terminal_panel = node_bivplot)
plot(gctest12, terminal_panel = node_bivplot)
plot(gmtest1, terminal_panel = node_bivplot)
plot(gstest1, terminal_panel = node_bivplot)
plot(gctest1, terminal_panel = node_bivplot)
plot(gmtest2, terminal_panel = node_bivplot)
plot(gstest2, terminal_panel = node_bivplot)
plot(gctest2, terminal_panel = node_bivplot)
}




###################################################
# various data sets


sim <- function(nobs = 100, nrep = 100, seed = 7, stump = TRUE,
                    formula = y~x|z1+z2+z3+z4+z5+z6+z7+z8+z9+z10,
                    vary_beta = "all", beta0 = 0, beta1 = 1, xi = 0, delta = 1, 
                    binary_regressor = TRUE, binary_beta = TRUE, z1dist = "norm",
                    sigma = 1, only_intercept = FALSE, alpha = 0.05,
                    test = c("ctree", "mfluc", "ctree_cat", "mfluc_cat",
                             "guide_sum_12", "guide_max_12", "guide_coin_12",
                             "guide_sum_1", "guide_max_1", "guide_coin_1",
                             "guide_sum_2", "guide_max_2", "guide_coin_2",
                             "guide_sum_1_cor", "guide_max_1_cor", "guide_coin_1_cor",
                             "guide_sum_2_cor", "guide_max_2_cor", "guide_coin_2_cor"))
{
  set.seed(seed)
  
  ## call
  cl <- match.call()
  
  pval <- matrix(rep(NA, length(test)*nrep), ncol = length(test))
  sv <- matrix(rep(NA, length(test)*nrep), ncol = length(test))
  colnames(pval) <- colnames(sv) <- test
  
  
  #pval_ctest <- pval_mtest <- pval_cctest <- pval_mctest <- pval_gstest12 <- pval_gmtest12 <- pval_gctest12 <- pval_gstest1 <- 
  #  pval_gmtest1 <-   pval_gctest1 <- pval_gstest2 <- pval_gmtest2 <- pval_gctest2 <- numeric(length = nrep)
  #
  #sv_ctest <- sv_mtest <- sv_cctest <- sv_mctest <- sv_gstest12 <- sv_gmtest12 <- sv_gctest12 <- sv_gstest1 <- 
  #  sv_gmtest1 <-   sv_gctest1 <- sv_gstest2 <- sv_gmtest2 <- sv_gctest2 <- character(length = nrep)
  
  for(i in 1:nrep){
    
    d <- dgp(nobs = nobs, vary_beta = vary_beta, beta0 = beta0, beta1 = beta1, xi = xi, 
             binary_regressor = binary_regressor, delta = delta,
             binary_beta = binary_beta, z1dist = z1dist,
             seed = seed + i, sigma = sigma, only_intercept = only_intercept)
    

    compute_pval <- function(test) {
      test <- match.arg(test, c("ctree", "mfluc", "ctree_cat", "mfluc_cat",
                                "guide_sum_12", "guide_max_12", "guide_coin_12",
                                "guide_sum_1", "guide_max_1", "guide_coin_1",
                                "guide_sum_2", "guide_max_2", "guide_coin_2",
                                "guide_sum_1_cor", "guide_max_1_cor", "guide_coin_1_cor",
                                "guide_sum_2_cor", "guide_max_2_cor", "guide_coin_2_cor"
                                ))
      testres <- switch(test,
                        "ctree" = evaltests(formula, data = d, testfun = "ctree", stump = stump, decorrelate = "vcov"),
                        "mfluc" = evaltests(formula, data = d, testfun = "mfluc", stump = stump, decorrelate = "vcov"),
                        "ctree_cat" = evaltests(formula, data = d, testfun = "ctree_cat", stump = stump, decorrelate = "vcov"),
                        "mfluc_cat" = evaltests(formula, data = d, testfun = "mfluc_cat", stump = stump, decorrelate = "vcov"),
                        "guide_sum_12" = evaltests(formula, data = d, testfun = "guide", 
                                              guide_testtype = "sum", guide_parm = c(1,2),
                                              xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "vcov"),
                        "guide_max_12" = evaltests(formula, data = d, testfun = "guide", 
                                              guide_testtype = "max", guide_parm = c(1,2),
                                              xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "vcov"),
                        "guide_coin_12" = evaltests(formula, data = d, testfun = "guide", 
                                              guide_testtype = "coin", guide_parm = c(1,2),
                                              xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "vcov"),
                        "guide_sum_1" = evaltests(formula, data = d, testfun = "guide", 
                                             guide_testtype = "sum", guide_parm = c(1),
                                             xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "vcov"),
                        "guide_max_1" = evaltests(formula, data = d, testfun = "guide", 
                                             guide_testtype = "max", guide_parm = c(1),
                                             xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "vcov"),
                        "guide_coin_1" = evaltests(formula, data = d, testfun = "guide", 
                                             guide_testtype = "coin", guide_parm = c(1),
                                             xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "vcov"),
                        "guide_sum_2" = evaltests(formula, data = d, testfun = "guide", 
                                             guide_testtype = "sum", guide_parm = c(2),
                                             xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "vcov"),
                        "guide_max_2" = evaltests(formula, data = d, testfun = "guide", 
                                             guide_testtype = "max", guide_parm = c(2),
                                             xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "vcov"),
                        "guide_coin_2" = evaltests(formula, data = d, testfun = "guide", 
                                             guide_testtype = "coin", guide_parm = c(2),
                                             xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "vcov"),
                        "guide_sum_1_cor" = evaltests(formula, data = d, testfun = "guide", 
                                                  guide_testtype = "sum", guide_parm = c(1),
                                                  xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "none"),
                        "guide_max_1_cor" = evaltests(formula, data = d, testfun = "guide", 
                                                  guide_testtype = "max", guide_parm = c(1),
                                                  xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "none"),
                        "guide_coin_1_cor" = evaltests(formula, data = d, testfun = "guide", 
                                                   guide_testtype = "coin", guide_parm = c(1),
                                                   xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "none"),
                        "guide_sum_2_cor" = evaltests(formula, data = d, testfun = "guide", 
                                                  guide_testtype = "sum", guide_parm = c(2),
                                                  xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "none"),
                        "guide_max_2_cor" = evaltests(formula, data = d, testfun = "guide", 
                                                  guide_testtype = "max", guide_parm = c(2),
                                                  xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "none"),
                        "guide_coin_2_cor" = evaltests(formula, data = d, testfun = "guide", 
                                                   guide_testtype = "coin", guide_parm = c(2),
                                                   xgroups = NULL, ygroups = NULL, stump = stump, decorrelate = "none")
                        )
      list(pval = as.numeric(info_node(testres$node)$p.value),
           sv = names(info_node(testres$node)$p.value))
    }
    
    resmat <- sapply(test, compute_pval)
    pval[i,] <- unlist(resmat["pval",])
    sv[i,] <- unlist(resmat["sv",])
    
  }
  
  
  prop_nosplit <- colMeans(pval > alpha)
  
  prop_split <- colMeans(pval <= alpha)
  
  pval_T <- pval
  pval_T[which(sv !="z1")] <- 1
  prop_T <- colMeans(pval_T < alpha)
  
  pval_F <- pval
  pval_F[which(sv =="z1")] <- 1
  prop_F <- colMeans(pval_F < alpha)
  
  return(list(prop_nosplit = prop_nosplit,
              prop_split = prop_split,
              prop_F = prop_F,
              prop_T = prop_T))
}





if(FALSE){
  simtest <- sim(nobs = 100, nrep = 100, seed = 7,
                 formula = y~x|z1+z2+z3+z4+z5+z6+z7+z8+z9+z10,
                 vary_beta = "beta1", beta0=1, xi = 0.8, delta = 2,
                 stump = TRUE, binary_regressor = FALSE,
                 binary_beta = FALSE, z1dist = "unif",
                 sigma = 1, only_intercept = FALSE,
                 test = c("ctree","mfluc"))
}


############# wrapper function applying simcomp over varying variables of interest

simwrapper <- function(nobs = 200, nrep = 100,
                       delta = seq(from = 1, to = 5, by = 2),
                       xi = c(0, 0.8), vary_beta = c("all", "beta0", "beta1"),
                       binary_regressor = c(TRUE, FALSE),
                       binary_beta = c(TRUE, FALSE),
                       only_intercept = c(TRUE, FALSE),
                       test = c("ctree", "mfluc", "ctree_cat", "mfluc_cat",
                                "guide_sum_12", "guide_max_12", "guide_coin_12",
                                "guide_sum_1", "guide_max_1", "guide_coin_1",
                                "guide_sum_2", "guide_max_2", "guide_coin_2",
                                "guide_sum_1_cor", "guide_max_1_cor", "guide_coin_1_cor",
                                "guide_sum_2_cor", "guide_max_2_cor", "guide_coin_2_cor"),
                       beta0=0, beta1 = 1,
                       stump = TRUE, z1dist = "unif", sigma = 1, alpha = 0.05)
{
  
  prs <- expand.grid(delta = delta, xi = xi, vary_beta = vary_beta,
                     binary_regressor = binary_regressor,
                     binary_beta = binary_beta,
                     only_intercept = only_intercept)
  
  rmid <- which((prs$only_intercept == TRUE & prs$vary_beta != "beta0"))
  if(length(rmid) > 0) prs <- prs[-rmid,]
  rmid <- which(prs$binary_beta == FALSE & prs$xi != 0)
  if(length(rmid) > 0) prs <- prs[-rmid,]
  rmid <- which(prs$only_intercept == TRUE & prs$binary_regressor == FALSE)
  if(length(rmid) > 0) prs <- prs[-rmid,]
  
  rownames(prs) <- c(1:NROW(prs))

  nprs <- nrow(prs)
  ntest <- length(test)
  prop_nosplit <- matrix(rep(NA, ntest * nprs), ncol = ntest)
  prop_split <- matrix(rep(NA, ntest * nprs), ncol = ntest)
  prop_F <- matrix(rep(NA, ntest * nprs), ncol = ntest)
  prop_T <- matrix(rep(NA, ntest * nprs), ncol = ntest)
  
  for(i in 1:nprs) {
    reslist <- sim(nobs = nobs, nrep = nrep, stump = stump,
                   formula = y~x|z1+z2+z3+z4+z5+z6+z7+z8+z9+z10,
                   beta0 = beta0, beta1 = beta1, z1dist = z1dist,
                   sigma = sigma, alpha = alpha,
                   test = test,
                   vary_beta = prs$vary_beta[i], 
                   binary_regressor = prs$binary_regressor[i],
                   binary_beta = prs$binary_beta[i],
                   xi = prs$xi[i],
                   only_intercept = prs$only_intercept[i],
                   delta = prs$delta[i])
    
    prop_nosplit[i,] <- reslist$prop_nosplit
    prop_split[i,] <- reslist$prop_split
    prop_F[i,] <- reslist$prop_F
    prop_T[i,] <- reslist$prop_T
  }
  
  rval <- data.frame()
  for(i in 1:ntest) rval <- rbind(rval, prs)
  rval$test <- gl(ntest, nprs, labels = test)
  rval$prop_nosplit <- as.vector(prop_nosplit)
  rval$prop_split <- as.vector(prop_split)
  rval$prop_F <- as.vector(prop_F)
  rval$prop_T <- as.vector(prop_T)
  rval$delta <- factor(rval$delta)
  rval$vary_beta <- factor(rval$vary_beta)
  rval$binary_regressor <- factor(rval$binary_regressor)
  rval$binary_beta <- factor(rval$binary_beta)
  rval$xi <- factor(rval$xi)
  rval$only_intercept <- factor(rval$only_intercept)
  
  return(rval)
}





simres <- simwrapper(nobs = 250, nrep = 100,
                     delta = seq(from = 0.5, to = 1.5, by = 0.5),
                     xi = c(0, 0.5, 0.8), vary_beta = c("all", "beta0", "beta1"),
                     #binary_regressor = c(TRUE, FALSE),
                     binary_regressor = FALSE,
                     binary_beta = c(TRUE, FALSE),
                     #only_intercept = c(TRUE, FALSE),
                     only_intercept = FALSE,
                     test = c("ctree", "mfluc", "ctree_cat", "mfluc_cat",
                              "guide_sum_12", "guide_max_12", "guide_coin_12",
                              "guide_sum_1", "guide_max_1", "guide_coin_1",
                              "guide_sum_2", "guide_max_2", "guide_coin_2",
                              "guide_sum_1_cor", "guide_max_1_cor", "guide_coin_1_cor",
                              "guide_sum_2_cor", "guide_max_2_cor", "guide_coin_2_cor"),
                     beta0 = 0, beta1 = 1,
                     stump = TRUE, z1dist = "unif", sigma = 1, alpha = 0.05)


save(simres, file = "~/svn/partykit/pkg/partykit/inst/guideliketest/sim/simres20180207.rda")

library("lattice")
load("~/svn/partykit/pkg/partykit/inst/guideliketest/sim/simres20180206.rda")

xyplot(prop_split ~ vary_beta | xi + binary_beta, groups = ~ test, data = simres, type = "b", auto.key = TRUE)
xyplot(prop_split ~ delta | xi + vary_beta + binary_beta, groups = ~ test, data = simres, type = "b", auto.key = TRUE)


# models considering all scores
s_allscores <- subset(simres, test %in% c("ctree", "mfluc", "ctree_cat", "mfluc_cat", 
                                          "guide_sum_12", "guide_max_12", "guide_coin_12"))
s_allscores$test <- factor(s_allscores$test)
xyplot(prop_split ~ xi | vary_beta + binary_beta, groups = ~ test, data = s_allscores, type = "b", auto.key = TRUE)
xyplot(prop_split ~ delta | xi + vary_beta + binary_beta, groups = ~ test, data = s_allscores, type = "b", auto.key = TRUE)


# all GUIDE models
s_guide <- subset(simres, test %in% c("guide_sum_12", "guide_max_12", "guide_coin_12",
                                      "guide_sum_1", "guide_max_1", "guide_coin_1",
                                      "guide_sum_2", "guide_max_2", "guide_coin_2",
                                      "guide_sum_1_cor", "guide_max_1_cor", "guide_coin_1_cor",
                                      "guide_sum_2_cor", "guide_max_2_cor", "guide_coin_2_cor"))
s_guide$test <- factor(s_guide$test)
xyplot(prop_split ~ xi | vary_beta + binary_beta, groups = ~ test, data = s_guide, type = "b", auto.key = TRUE)
xyplot(prop_split ~ delta | xi + vary_beta + binary_beta, groups = ~ test, data = s_guide, type = "b", auto.key = TRUE)


# GUIDE_12 models
s_guide12 <- subset(simres, test %in% c("guide_sum_12", "guide_max_12", "guide_coin_12"))
s_guide12$test <- factor(s_guide12$test)
xyplot(prop_split ~ xi | vary_beta + binary_beta, groups = ~ test, data = s_guide12, type = "b", auto.key = TRUE)
xyplot(prop_split ~ delta | xi + vary_beta + binary_beta, groups = ~ test, data = s_guide12, type = "b", auto.key = TRUE)


# GUIDE_1 models
s_guide1 <- subset(simres, test %in% c("guide_sum_1", "guide_max_1", "guide_coin_1",
                                       "guide_sum_1_cor", "guide_max_1_cor", "guide_coin_1_cor"))
s_guide1$test <- factor(s_guide1$test)
xyplot(prop_split ~ xi | vary_beta + binary_beta, groups = ~ test, data = s_guide1, type = "b", auto.key = TRUE)
xyplot(prop_split ~ delta | xi + vary_beta + binary_beta, groups = ~ test, data = s_guide1, type = "b", auto.key = TRUE)


# GUIDE_2 models
s_guide2 <- subset(simres, test %in% c("guide_sum_2", "guide_max_2", "guide_coin_2",
                                       "guide_sum_2_cor", "guide_max_2_cor", "guide_coin_2_cor"))
s_guide2$test <- factor(s_guide2$test)
xyplot(prop_split ~ xi | vary_beta + binary_beta, groups = ~ test, data = s_guide2, type = "b", auto.key = TRUE)
xyplot(prop_split ~ delta | xi + vary_beta + binary_beta, groups = ~ test, data = s_guide2, type = "b", auto.key = TRUE)


# ctree, mfluc, GUIDE_12_max
s_scores2 <- subset(simres, test %in% c("ctree", "mfluc", "ctree_cat", "mfluc_cat", "guide_max_12"))
s_scores2$test <- factor(s_scores2$test)
xyplot(prop_split ~ xi | vary_beta + binary_beta, groups = ~ test, data = s_scores2, type = "b", auto.key = TRUE)
xyplot(prop_split ~ delta | xi + vary_beta + binary_beta, groups = ~ test, data = s_scores2, type = "b", auto.key = TRUE)




tab <- xtabs(prop_T ~ delta + xi + binary_beta + binary_regressor + test,
             data = simres)
ftable(tab, row.vars = c("xi", "binary_beta", "binary_regressor", "delta"), 
       col.vars = "test")






if(FALSE){
  b1fix_contbeta_contreg <- simcomp(nobs = 100, nrep = 100, seed = 7,
                                    formula = y~x|z1+z2+z3+z4+z5+z6+z7+z8+z9+z10,
                                    z1dist = "unif", sigma = 1, only_intercept = FALSE,
                                    vary_beta = "beta0", beta1 = 1, xi = 0.0, delta = 2,
                                    stump = TRUE, binary_regressor = FALSE,
                                    binary_beta = FALSE)
  
  save(b1fix_contbeta_contreg, file = "~/svn/partykit/pkg/partykit/inst/guideliketest/sim/b1fix_contbeta_contreg.rda")
  
  
  
  load("~/svn/partykit/pkg/partykit/inst/guideliketest/sim/b1fix_contbeta_contreg.rda")
  simtest <- b1fix_contbeta_contreg
  
  # mean over all p-values (incl. p-values > 0.05)
  mean_all <- cbind(mean(simtest$pval_ctest),
                    mean(simtest$pval_mtest),
                    mean(simtest$pval_cctest),
                    mean(simtest$pval_mctest),
                    mean(simtest$pval_gstest12),
                    mean(simtest$pval_gmtest12),
                    mean(simtest$pval_gctest12),
                    mean(simtest$pval_gstest1),
                    mean(simtest$pval_gmtest1),
                    mean(simtest$pval_gctest1),
                    mean(simtest$pval_gstest2),
                    mean(simtest$pval_gmtest2),
                    mean(simtest$pval_gctest2))
  colnames(mean_all) <-  c("c","m", "cc", "mc", "gs12","gm12","gc12", "gs1","gm1","gc1", "gs2","gm2","gc2")
  mean_all
  
  # percentage of trees with correct split variable (incl. p-values > 0.05)
  prop_var <- cbind(sum(simtest$sv_ctest == "z1")/length(simtest$sv_ctest),
                    sum(simtest$sv_mtest == "z1")/length(simtest$sv_mtest),
                    sum(simtest$sv_cctest == "z1")/length(simtest$sv_cctest),
                    sum(simtest$sv_mctest == "z1")/length(simtest$sv_mctest),
                    sum(simtest$sv_gstest12 == "z1")/length(simtest$sv_gstest12),
                    sum(simtest$sv_gmtest12 == "z1")/length(simtest$sv_gmtest12),
                    sum(simtest$sv_gctest12 == "z1")/length(simtest$sv_gctest12),
                    sum(simtest$sv_gstest1 == "z1")/length(simtest$sv_gstest1),
                    sum(simtest$sv_gmtest1 == "z1")/length(simtest$sv_gmtest1),
                    sum(simtest$sv_gctest1 == "z1")/length(simtest$sv_gctest1),
                    sum(simtest$sv_gstest2 == "z1")/length(simtest$sv_gstest2),
                    sum(simtest$sv_gmtest2 == "z1")/length(simtest$sv_gmtest2),
                    sum(simtest$sv_gctest2 == "z1")/length(simtest$sv_gctest2))
  colnames(prop_var) <-  c("c","m", "cc", "mc", "gs12","gm12","gc12", "gs1","gm1","gc1", "gs2","gm2","gc2")
  barplot(prop_var, ylim = c(0,1))
  
  # percentage of trees with p-values <= 0.05
  prop_005 <- cbind(sum(simtest$pval_ctest <= 0.05)/length(simtest$sv_ctest),
                    sum(simtest$pval_mtest <= 0.05)/length(simtest$sv_mtest),
                    sum(simtest$pval_cctest <= 0.05)/length(simtest$sv_cctest),
                    sum(simtest$pval_mctest <= 0.05)/length(simtest$sv_mctest),
                    sum(simtest$pval_gstest12 <= 0.05)/length(simtest$sv_gstest12),
                    sum(simtest$pval_gmtest12 <= 0.05)/length(simtest$sv_gmtest12),
                    sum(simtest$pval_gctest12 <= 0.05)/length(simtest$sv_gctest12),
                    sum(simtest$pval_gstest1 <= 0.05)/length(simtest$sv_gstest1),
                    sum(simtest$pval_gmtest1 <= 0.05)/length(simtest$sv_gmtest1),
                    sum(simtest$pval_gctest1 <= 0.05)/length(simtest$sv_gctest1),
                    sum(simtest$pval_gstest2 <= 0.05)/length(simtest$sv_gstest2),
                    sum(simtest$pval_gmtest2 <= 0.05)/length(simtest$sv_gmtest2),
                    sum(simtest$pval_gctest2 <= 0.05)/length(simtest$sv_gctest2))
  colnames(prop_005) <-  c("c","m", "cc", "mc", "gs12","gm12","gc12", "gs1","gm1","gc1", "gs2","gm2","gc2")
  barplot(prop_005, ylim = c(0,1))
  
  # percentage of trees with correct split variable and p-values <= 0.05
  prop_var005 <- cbind(sum(simtest$sv_ctest == "z1" & simtest$pval_ctest <= 0.05)/length(simtest$sv_ctest),
                       sum(simtest$sv_mtest == "z1" & simtest$pval_mtest <= 0.05)/length(simtest$sv_mtest),
                       sum(simtest$sv_cctest == "z1" & simtest$pval_cctest <= 0.05)/length(simtest$sv_cctest),
                       sum(simtest$sv_mctest == "z1" & simtest$pval_mctest <= 0.05)/length(simtest$sv_mctest),
                       sum(simtest$sv_gstest12 == "z1" & simtest$pval_gstest12 <= 0.05)/length(simtest$sv_gstest12),
                       sum(simtest$sv_gmtest12 == "z1" & simtest$pval_gmtest12 <= 0.05)/length(simtest$sv_gmtest12),
                       sum(simtest$sv_gctest12 == "z1" & simtest$pval_gctest12 <= 0.05)/length(simtest$sv_gctest12),
                       sum(simtest$sv_gstest1 == "z1" & simtest$pval_gstest1 <= 0.05)/length(simtest$sv_gstest1),
                       sum(simtest$sv_gmtest1 == "z1" & simtest$pval_gmtest1 <= 0.05)/length(simtest$sv_gmtest1),
                       sum(simtest$sv_gctest1 == "z1" & simtest$pval_gctest1 <= 0.05)/length(simtest$sv_gctest1),
                       sum(simtest$sv_gstest2 == "z1" & simtest$pval_gstest2 <= 0.05)/length(simtest$sv_gstest2),
                       sum(simtest$sv_gmtest2 == "z1" & simtest$pval_gmtest2 <= 0.05)/length(simtest$sv_gmtest2),
                       sum(simtest$sv_gctest2 == "z1" & simtest$pval_gctest2 <= 0.05)/length(simtest$sv_gctest2))
  colnames(prop_var005) <-  c("c","m", "cc", "mc", "gs12","gm12","gc12", "gs1","gm1","gc1", "gs2","gm2","gc2")
  barplot(prop_var005, ylim = c(0,1))
  
  
  # boxplot of all p-values
  boxplot(simtest$pval_ctest, simtest$pval_mtest, 
          simtest$pval_cctest, simtest$pval_mctest, 
          simtest$pval_gstest12, simtest$pval_gmtest12, simtest$pval_gctest12,
          simtest$pval_gstest1, simtest$pval_gmtest1, simtest$pval_gctest1,
          simtest$pval_gstest2, simtest$pval_gmtest2, simtest$pval_gctest2,
          outline = FALSE, 
          names = c("c","m", "cc", "mc",
                    "gs12","gm12","gc12",
                    "gs1","gm1","gc1",
                    "gs2","gm2","gc2"))
  
  
  # boxplot of all p-values if correct split variable was found
  boxplot(simtest$pval_ctest[simtest$sv_ctest == "z1"], simtest$pval_mtest[simtest$sv_mtest == "z1"], 
          simtest$pval_cctest[simtest$sv_cctest == "z1"], simtest$pval_mctest[simtest$sv_mctest == "z1"], 
          simtest$pval_gstest12[simtest$sv_gstest12 == "z1"], simtest$pval_gmtest12[simtest$sv_gmtest12 == "z1"], 
          simtest$pval_gctest12[simtest$sv_gctest12 == "z1"],
          simtest$pval_gstest1[simtest$sv_gstest1 == "z1"], simtest$pval_gmtest1[simtest$sv_gmtest1 == "z1"], 
          simtest$pval_gctest1[simtest$sv_gctest1 == "z1"],
          simtest$pval_gstest2[simtest$sv_gstest2 == "z1"], simtest$pval_gmtest2[simtest$sv_gmtest2 == "z1"], 
          simtest$pval_gctest2[simtest$sv_gctest2 == "z1"],
          outline = FALSE, 
          names = c("c","m", "cc", "mc",
                    "gs12","gm12","gc12",
                    "gs1","gm1","gc1",
                    "gs2","gm2","gc2"))
  
  
  boxplot(simtest$pval_ctest, simtest$pval_mtest, 
          simtest$pval_cctest, simtest$pval_mctest, 
          simtest$pval_gstest12, simtest$pval_gmtest12, simtest$pval_gctest12,
          simtest$pval_gstest1, simtest$pval_gmtest1, simtest$pval_gctest1,
          #simtest$pval_gstest2, simtest$pval_gmtest2, simtest$pval_gctest2,
          outline = FALSE, 
          names = c("c", "m", "cc", "mc",
                    "gs12","gm12", "gc12",
                    "gs1","gm1","gc1")) #,
  #"gs2","gm2","gc2"))
  
  
}