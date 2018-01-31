library("partykit")
library("Formula")

comptests <- function(formula, data, testfun = c("guide", "ctree", "mfluc"), 
                      subset, weights, stump = TRUE,
                      guide_testtype = c("sum", "max", "coin"), 
                      decorrelate = "vcov",
                      guide_parm = NULL,
                      xgroups = NULL,
                      ygroups = NULL){
  
  na.action = na.pass
  converged = NULL
  scores = NULL
  doFit = TRUE
  offset = NULL 
  cluster = NULL
  
  if(length(testfun) > 1) testfun <- testfun[1]
  if(!testfun %in% c("guide", "ctree", "mfluc")) stop("testfun hast to be one of the following options: 'guide', 'ctree', 'mfluc'")
  
  
  if(testfun == "guide"){
    
    if(length(guide_testtype) > 1) guide_testtype <- guide_testtype[1]
    if(!guide_testtype %in% c("max", "sum", "coin")) stop("testfun hast to be one of the following options: 'max', 'sum', 'coin'")
    
    control <- partykit:::extree_control(criterion = "p.value",
                                         selectfun = .guide_select(guide_decorrelate = "none",    # because of decorrelation within ytrafo 
                                                                   guide_testtype = guide_testtype,
                                                                   guide_parm = guide_parm,
                                                                   xgroups = xgroups,
                                                                   ygroups = ygroups,
                                                                   interaction = FALSE),  
                                         splitfun = partykit:::.objfun_split(restart = TRUE,
                                                                             intersplit = TRUE),
                                         svselectfun = .guide_select(guide_decorrelate = "none",    # because of decorrelation within ytrafo 
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
  
  
  control$inner <- "object"
  control$temrinal <- "object"
  
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
    decorrelate <- "vcov"
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



###########
# DGP

dgp <- function(n = 100, delta = 1, xi = 0,
                sigma = 1, seed = 7,
                only_intercept = FALSE, 
                binary_regressor = TRUE, 
                vary_beta = c("all", "beta0", "beta1"),
                binary_beta = TRUE,
                beta0 = NULL,
                beta1 = NULL){
  
  if(length(vary_beta) > 1) vary_beta <- vary_beta[1]
  if(!vary_beta %in% c("all", "beta0", "beta1")) stop("vary_beta has to be one of the following options: 'all', 'beta0', 'beta1'")
  
  set.seed(seed)
  
  z1 <- rnorm(n, 0, 1)
  z2 <- rnorm(n, 0, 1)
  z3 <- rnorm(n, 0, 1)
  z4 <- rnorm(n, 0, 1)
  z5 <- rnorm(n, 0, 1)
  z6 <- runif(n, -1, 1)
  z7 <- runif(n, -1, 1)
  z8 <- runif(n, -1, 1)
  z9 <- runif(n, -1, 1)
  z10 <- runif(n, -1, 1)
  
  if(!only_intercept)  x <- if(binary_regressor) (-1)^rbinom(n, 1, 0.5) else runif(100, min = -1, max = 1)
  
  e <- rnorm(n, 0, sigma)
  
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




############
# compare

d <- dgp(400, vary_beta = "beta1", beta0 = 3)
d <- dgp(400, vary_beta = "all", xi = -0.5)
d <- dgp(1000, vary_beta = "all", xi = -0.5, binary_regressor = FALSE)
d <- dgp(1000, vary_beta = "all", xi = -0.0, binary_regressor = FALSE)

d <- dgp(400, vary_beta = "all", xi = -0.4, delta = 3, binary_regressor = TRUE)
d <- dgp(200, vary_beta = "beta0", beta1 = 2, xi = 0.0, binary_regressor = FALSE)
d <- dgp(200, vary_beta = "beta1", beta0 = 0, xi = -0.5, binary_regressor = FALSE)
d <- dgp(200, vary_beta = "beta1", beta0 = 0, xi = 0.0, binary_regressor = FALSE)


ctest <- comptests(y~x|z1+z2+z3, data = d, testfun = "ctree")
mtest <- comptests(y~x|z1+z2+z3, data = d, testfun = "mfluc")
gstest12 <- comptests(y~x|z1+z2+z3, data = d, testfun = "guide", 
                    guide_testtype = "sum", guide_parm = c(1,2),
                    xgroups = NULL, ygroups = NULL)
gmtest12 <- comptests(y~x|z1+z2+z3, data = d, testfun = "guide", 
                    guide_testtype = "max", guide_parm = c(1,2),
                    xgroups = NULL, ygroups = NULL)
gctest12 <- comptests(y~x|z1+z2+z3, data = d, testfun = "guide", 
                    guide_testtype = "coin", guide_parm = c(1,2),
                    xgroups = NULL, ygroups = NULL)
gstest1 <- comptests(y~x|z1+z2+z3, data = d, testfun = "guide", 
                      guide_testtype = "sum", guide_parm = c(1),
                      xgroups = NULL, ygroups = NULL)
gmtest1 <- comptests(y~x|z1+z2+z3, data = d, testfun = "guide", 
                      guide_testtype = "max", guide_parm = c(1),
                      xgroups = NULL, ygroups = NULL)
gctest1 <- comptests(y~x|z1+z2+z3, data = d, testfun = "guide", 
                      guide_testtype = "coin", guide_parm = c(1),
                      xgroups = NULL, ygroups = NULL)
gstest2 <- comptests(y~x|z1+z2+z3, data = d, testfun = "guide", 
                      guide_testtype = "sum", guide_parm = c(2),
                      xgroups = NULL, ygroups = NULL)
gmtest2 <- comptests(y~x|z1+z2+z3, data = d, testfun = "guide", 
                      guide_testtype = "max", guide_parm = c(2),
                      xgroups = NULL, ygroups = NULL)
gctest2 <- comptests(y~x|z1+z2+z3, data = d, testfun = "guide", 
                      guide_testtype = "coin", guide_parm = c(2),
                      xgroups = NULL, ygroups = NULL)


plot(ctest)
plot(mtest)
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
plot(gmtest12, terminal_panel = node_bivplot)
plot(gstest12, terminal_panel = node_bivplot)
plot(gctest12, terminal_panel = node_bivplot)
plot(gmtest1, terminal_panel = node_bivplot)
plot(gstest1, terminal_panel = node_bivplot)
plot(gctest1, terminal_panel = node_bivplot)
plot(gmtest2, terminal_panel = node_bivplot)
plot(gstest2, terminal_panel = node_bivplot)
plot(gctest2, terminal_panel = node_bivplot)





###################################################
# various data sets


simcomp <- function(nobs = 100, nrep = 100, seed = 7, stump = TRUE,
                    vary_beta = "all", beta0 = 0, beta1 = 1, xi = 0, delta = 1, binary_regressor = TRUE)
{
  set.seed(seed)
  
  ## call
  cl <- match.call()
  
  pval_ctest <- pval_mtest <- pval_gstest12 <- pval_gmtest12 <- pval_gctest12 <- pval_gstest1 <- 
    pval_gmtest1 <-   pval_gctest1 <- pval_gstest2 <- pval_gmtest2 <- pval_gctest2 <- numeric(length = nrep)
  
  sv_ctest <- sv_mtest <- sv_gstest12 <- sv_gmtest12 <- sv_gctest12 <- sv_gstest1 <- 
    sv_gmtest1 <-   sv_gctest1 <- sv_gstest2 <- sv_gmtest2 <- sv_gctest2 <- character(length = nrep)
  
  for(i in 1:nrep){
    
    d <- dgp(n = nobs, vary_beta = vary_beta, beta0 = beta0, beta1 = beta1, xi = xi, 
             binary_regressor = binary_regressor, seed = seed + i)
    
    ctest <- comptests(y~x|z1+z2+z3, data = d, testfun = "ctree", stump = stump)
    mtest <- comptests(y~x|z1+z2+z3, data = d, testfun = "mfluc", stump = stump)
    gstest12 <- comptests(y~x|z1+z2+z3, data = d, testfun = "guide", 
                          guide_testtype = "sum", guide_parm = c(1,2),
                          xgroups = NULL, ygroups = NULL, stump = stump)
    gmtest12 <- comptests(y~x|z1+z2+z3, data = d, testfun = "guide", 
                          guide_testtype = "max", guide_parm = c(1,2),
                          xgroups = NULL, ygroups = NULL, stump = stump)
    gctest12 <- comptests(y~x|z1+z2+z3, data = d, testfun = "guide", 
                          guide_testtype = "coin", guide_parm = c(1,2),
                          xgroups = NULL, ygroups = NULL, stump = stump)
    gstest1 <- comptests(y~x|z1+z2+z3, data = d, testfun = "guide", 
                         guide_testtype = "sum", guide_parm = c(1),
                         xgroups = NULL, ygroups = NULL, stump = stump)
    gmtest1 <- comptests(y~x|z1+z2+z3, data = d, testfun = "guide", 
                         guide_testtype = "max", guide_parm = c(1),
                         xgroups = NULL, ygroups = NULL, stump = stump)
    gctest1 <- comptests(y~x|z1+z2+z3, data = d, testfun = "guide", 
                         guide_testtype = "coin", guide_parm = c(1),
                         xgroups = NULL, ygroups = NULL, stump = stump)
    gstest2 <- comptests(y~x|z1+z2+z3, data = d, testfun = "guide", 
                         guide_testtype = "sum", guide_parm = c(2),
                         xgroups = NULL, ygroups = NULL, stump = stump)
    gmtest2 <- comptests(y~x|z1+z2+z3, data = d, testfun = "guide", 
                         guide_testtype = "max", guide_parm = c(2),
                         xgroups = NULL, ygroups = NULL, stump = stump)
    gctest2 <- comptests(y~x|z1+z2+z3, data = d, testfun = "guide", 
                         guide_testtype = "coin", guide_parm = c(2),
                         xgroups = NULL, ygroups = NULL, stump = stump)
    
    
    pval_ctest[i] <- info_node(ctest$node)$p.value
    pval_mtest[i] <- info_node(mtest$node)$p.value
    pval_gstest12[i] <- info_node(gstest12$node)$p.value
    pval_gmtest12[i] <- info_node(gmtest12$node)$p.value
    pval_gctest12[i] <- info_node(gctest12$node)$p.value
    pval_gstest1[i] <- info_node(gstest1$node)$p.value
    pval_gmtest1[i] <- info_node(gmtest1$node)$p.value
    pval_gctest1[i] <- info_node(gctest1$node)$p.value
    pval_gstest2[i] <- info_node(gstest2$node)$p.value
    pval_gmtest2[i] <- info_node(gmtest2$node)$p.value
    pval_gctest2[i] <- info_node(gctest2$node)$p.value
    
    sv_ctest[i] <- names(info_node(ctest$node)$p.value)
    sv_mtest[i] <- names(info_node(mtest$node)$p.value)
    sv_gstest12[i] <- names(info_node(gstest12$node)$p.value)
    sv_gmtest12[i] <- names(info_node(gmtest12$node)$p.value)
    sv_gctest12[i] <- names(info_node(gctest12$node)$p.value)
    sv_gstest1[i] <- names(info_node(gstest1$node)$p.value)
    sv_gmtest1[i] <- names(info_node(gmtest1$node)$p.value)
    sv_gctest1[i] <- names(info_node(gctest1$node)$p.value)
    sv_gstest2[i] <- names(info_node(gstest2$node)$p.value)
    sv_gmtest2[i] <- names(info_node(gmtest2$node)$p.value)
    sv_gctest2[i] <- names(info_node(gctest2$node)$p.value)

  }
  
 return(list(pval_ctest = pval_ctest,
             pval_mtest = pval_mtest,
             pval_gstest12 = pval_gstest12,
             pval_gmtest12 = pval_gmtest12,
             pval_gctest12 = pval_gctest12,
             pval_gstest1 = pval_gstest1,
             pval_gmtest1 = pval_gmtest1,
             pval_gctest1 = pval_gctest1,
             pval_gstest2 = pval_gstest2,
             pval_gmtest2 = pval_gmtest2,
             pval_gctest2 = pval_gctest2,
             sv_ctest = sv_ctest,
             sv_mtest = sv_mtest,
             sv_gstest12 = sv_gstest12,
             sv_gmtest12 = sv_gmtest12,
             sv_gctest12 = sv_gctest12,
             sv_gstest1 = sv_gstest1,
             sv_gmtest1 = sv_gmtest1,
             sv_gctest1 = sv_gctest1,
             sv_gstest2 = sv_gstest2,
             sv_gmtest2 = sv_gmtest2,
             sv_gctest2 = sv_gctest2,
             call = cl))
}



simtest <- simcomp(nobs = 100, nrep = 100, seed = 7,
                   vary_beta = "beta1", beta0=0, xi = 0.5,
                   stump = TRUE, binary_regressor = FALSE)





# mean over all p-values (incl. p-values > 0.05)
mean_all <- cbind(mean(simtest$pval_ctest),
                  mean(simtest$pval_mtest),
                  mean(simtest$pval_gstest12),
                  mean(simtest$pval_gmtest12),
                  mean(simtest$pval_gctest12),
                  mean(simtest$pval_gstest1),
                  mean(simtest$pval_gmtest1),
                  mean(simtest$pval_gctest1),
                  mean(simtest$pval_gstest2),
                  mean(simtest$pval_gmtest2),
                  mean(simtest$pval_gctest2))
colnames(mean_all) <-  c("c","m", "gs12","gm12","gc12", "gs1","gm1","gc1", "gs2","gm2","gc2")
mean_all

# percentage of trees with correct split variable (incl. p-values > 0.05)
prop_var <- cbind(sum(simtest$sv_ctest == "z1")/length(simtest$sv_ctest),
                    sum(simtest$sv_mtest == "z1")/length(simtest$sv_mtest),
                    sum(simtest$sv_gstest12 == "z1")/length(simtest$sv_gstest12),
                    sum(simtest$sv_gmtest12 == "z1")/length(simtest$sv_gmtest12),
                    sum(simtest$sv_gctest12 == "z1")/length(simtest$sv_gctest12),
                    sum(simtest$sv_gstest1 == "z1")/length(simtest$sv_gstest1),
                    sum(simtest$sv_gmtest1 == "z1")/length(simtest$sv_gmtest1),
                    sum(simtest$sv_gctest1 == "z1")/length(simtest$sv_gctest1),
                    sum(simtest$sv_gstest2 == "z1")/length(simtest$sv_gstest2),
                    sum(simtest$sv_gmtest2 == "z1")/length(simtest$sv_gmtest2),
                    sum(simtest$sv_gctest2 == "z1")/length(simtest$sv_gctest2))
colnames(prop_var) <-  c("c","m", "gs12","gm12","gc12", "gs1","gm1","gc1", "gs2","gm2","gc2")
barplot(prop_var, ylim = c(0,1))

# percentage of trees with p-values <= 0.05
prop_005 <- cbind(sum(simtest$pval_ctest <= 0.05)/length(simtest$sv_ctest),
                  sum(simtest$pval_mtest <= 0.05)/length(simtest$sv_mtest),
                  sum(simtest$pval_gstest12 <= 0.05)/length(simtest$sv_gstest12),
                  sum(simtest$pval_gmtest12 <= 0.05)/length(simtest$sv_gmtest12),
                  sum(simtest$pval_gctest12 <= 0.05)/length(simtest$sv_gctest12),
                  sum(simtest$pval_gstest1 <= 0.05)/length(simtest$sv_gstest1),
                  sum(simtest$pval_gmtest1 <= 0.05)/length(simtest$sv_gmtest1),
                  sum(simtest$pval_gctest1 <= 0.05)/length(simtest$sv_gctest1),
                  sum(simtest$pval_gstest2 <= 0.05)/length(simtest$sv_gstest2),
                  sum(simtest$pval_gmtest2 <= 0.05)/length(simtest$sv_gmtest2),
                  sum(simtest$pval_gctest2 <= 0.05)/length(simtest$sv_gctest2))
colnames(prop_005) <-  c("c","m", "gs12","gm12","gc12", "gs1","gm1","gc1", "gs2","gm2","gc2")
barplot(prop_005, ylim = c(0,1))

# percentage of trees with correct split variable and p-values <= 0.05
prop_var005 <- cbind(sum(simtest$sv_ctest == "z1" & simtest$pval_ctest <= 0.05)/length(simtest$sv_ctest),
                       sum(simtest$sv_mtest == "z1" & simtest$pval_mtest <= 0.05)/length(simtest$sv_mtest),
                       sum(simtest$sv_gstest12 == "z1" & simtest$pval_gstest12 <= 0.05)/length(simtest$sv_gstest12),
                       sum(simtest$sv_gmtest12 == "z1" & simtest$pval_gmtest12 <= 0.05)/length(simtest$sv_gmtest12),
                       sum(simtest$sv_gctest12 == "z1" & simtest$pval_gctest12 <= 0.05)/length(simtest$sv_gctest12),
                       sum(simtest$sv_gstest1 == "z1" & simtest$pval_gstest1 <= 0.05)/length(simtest$sv_gstest1),
                       sum(simtest$sv_gmtest1 == "z1" & simtest$pval_gmtest1 <= 0.05)/length(simtest$sv_gmtest1),
                       sum(simtest$sv_gctest1 == "z1" & simtest$pval_gctest1 <= 0.05)/length(simtest$sv_gctest1),
                       sum(simtest$sv_gstest2 == "z1" & simtest$pval_gstest2 <= 0.05)/length(simtest$sv_gstest2),
                       sum(simtest$sv_gmtest2 == "z1" & simtest$pval_gmtest2 <= 0.05)/length(simtest$sv_gmtest2),
                       sum(simtest$sv_gctest2 == "z1" & simtest$pval_gctest2 <= 0.05)/length(simtest$sv_gctest2))
colnames(prop_var005) <-  c("c","m", "gs12","gm12","gc12", "gs1","gm1","gc1", "gs2","gm2","gc2")
barplot(prop_var005, ylim = c(0,1))



boxplot(simtest$pval_ctest, simtest$pval_mtest, 
        simtest$pval_gstest12, simtest$pval_gmtest12, simtest$pval_gctest12,
        simtest$pval_gstest1, simtest$pval_gmtest1, simtest$pval_gctest1,
        simtest$pval_gstest2, simtest$pval_gmtest2, simtest$pval_gctest2,
        outline = FALSE, 
        names = c("c","m",
                  "gs12","gm12","gc12",
                  "gs1","gm1","gc1",
                  "gs2","gm2","gc2"))


boxplot(simtest$pval_ctest, simtest$pval_mtest, 
        simtest$pval_gstest12, simtest$pval_gmtest12, #simtest$pval_gctest12,
        #simtest$pval_gstest1, simtest$pval_gmtest1, simtest$pval_gctest1,
        simtest$pval_gstest2, simtest$pval_gmtest2, #simtest$pval_gctest2,
        outline = FALSE, 
        names = c("c","m",
                  "gs12","gm12", #"gc12",
                  #"gs1","gm1","gc1",
                  "gs2","gm2")) #,"gc2"))


