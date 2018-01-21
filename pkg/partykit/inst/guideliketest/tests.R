library("partykit")
library("Formula")

comptests <- function(formula, data, testfun = c("guide", "ctree", "mfluc"), 
                      subset, weights, 
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
                                                                   ygroups = ygroups),  
                                         splitfun = partykit:::.objfun_split(restart = TRUE,
                                                                             intersplit = TRUE),
                                         svselectfun = .guide_select(guide_decorrelate = "none",    # because of decorrelation within ytrafo 
                                                                     guide_testtype = guide_testtype,
                                                                     guide_parm = guide_parm,
                                                                     xgroups = xgroups,
                                                                     ygroups = ygroups), 
                                         svsplitfun = partykit:::.objfun_split(restart = TRUE,
                                                                               intersplit = TRUE),
                                         logmincriterion = log(1-0.05),
                                         update = TRUE)
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
                                         update = TRUE)
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
                              update = TRUE)
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



ctest <- comptests(dist~speed, data = cars, testfun = "ctree")
mtest <- comptests(dist~speed, data = cars, testfun = "mfluc")
gstest <- comptests(dist~speed, data = cars, testfun = "guide", guide_testtype = "sum")
gmtest <- comptests(dist~speed, data = cars, testfun = "guide", guide_testtype = "max")

plot(ctest)
plot(mtest)
plot(gmtest)
plot(gstest)




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


ctest <- comptests(y~x|z1+z2+z3, data = d, testfun = "ctree")
mtest <- comptests(y~x|z1+z2+z3, data = d, testfun = "mfluc")
gstest <- comptests(y~x|z1+z2+z3, data = d, testfun = "guide", 
                    guide_testtype = "sum", guide_parm = c(1,2),
                    xgroups = NULL, ygroups = NULL)
gmtest <- comptests(y~x|z1+z2+z3, data = d, testfun = "guide", 
                    guide_testtype = "max", guide_parm = c(1,2),
                    xgroups = NULL, ygroups = NULL)

ctest

plot(ctest)
plot(mtest)
plot(gmtest)
plot(gstest)

plot(ctest, terminal_panel = node_bivplot)
plot(mtest, terminal_panel = node_bivplot)
plot(gmtest, terminal_panel = node_bivplot)
plot(gstest, terminal_panel = node_bivplot)





#####################################################
# data set

d <- extree_data(dist~speed, data = cars)


dg <- extree_fit(d, trafo = ytrafo, 
                 weights = rep.int(1, NROW(d$data)),
                 subset = c(1 : NROW(d$data)),
                 partyvars = d$variables$z,
                 converged = NULL,
                 ctrl = partykit:::extree_control(criterion = "p.value",
                                                  selectfun = .guide_select(guide_decorrelate = "none",    # because of decorrelation within ytrafo 
                                                                            guide_testtype = "sum"),  
                                                  splitfun = partykit:::.objfun_split(restart = TRUE,
                                                                                      intersplit = TRUE),
                                                  svselectfun = .guide_select(guide_decorrelate = "none",    # because of decorrelation within ytrafo 
                                                                              guide_testtype = "sum"), 
                                                  svsplitfun = partykit:::.objfun_split(restart = TRUE,
                                                                                        intersplit = TRUE),
                                                  logmincriterion = log(1-0.05),
                                                  update = TRUE))

dg$nodes
ctree(dist~speed, data = cars, ytrafo = ytrafo, control = ctree_control(update = TRUE,
                                                                        logmincriterion = log(1-0.05)))

party(dg$nodes, data = cars)
plot(party(dg$nodes, data =cars))
plot(ctree(dist~speed, data = cars))




#####################################################################################

# control arguments for ctree as testfunction:
ctest_control <- c(partykit:::extree_control(criterion = "p.value",
                                             logmincriterion = logmincriterion, minsplit = minsplit, 
                                             minbucket = minbucket, minprob = minprob, 
                                             nmax = nmax, stump = stump, lookahead = lookahead,
                                             mtry = mtry, maxdepth = maxdepth, multiway = multiway, 
                                             splittry = splittry, maxsurrogate = maxsurrogate, 
                                             numsurrogate = numsurrogate,
                                             majority = majority, caseweights = caseweights, 
                                             applyfun = applyfun, saveinfo = saveinfo,  ### always
                                             selectfun = .ctree_select(),
                                             splitfun = if (splitflavour == "ctree") .ctree_split() else .objfun_test(),  ## FIX ME .objfun_split()  (also fix in ctree)
                                             svselectfun = .ctree_select(),
                                             svsplitfun =.ctree_split(minbucket = 0),
                                             bonferroni = "Bonferroni" %in% testtype, 
                                             update = update),
                   list(teststat = teststat, splitstat = splitstat, splittest = splittest, pargs = pargs,
                        testtype = ttesttype, nresample = nresample, tol = tol,
                        intersplit = intersplit, MIA = MIA))

testctrl <- partykit:::extree_control(criterion = "p.value", 
                                      logmincriterion = log(0.05), 
                                      minsplit = 20L,
                                      minbucket = 7L, 
                                      minprob = 0.01, 
                                      nmax = Inf,
                                      stump = FALSE,
                                      lookahead = FALSE, ### try trafo() for daugther nodes before implementing the split
                                      maxsurrogate = 0L, 
                                      numsurrogate = FALSE,
                                      mtry = Inf,
                                      maxdepth = Inf, 
                                      multiway = FALSE, 
                                      splittry = 2L,
                                      majority = FALSE, 
                                      caseweights = TRUE, 
                                      applyfun = NULL, 
                                      cores = NULL,
                                      saveinfo = TRUE,
                                      bonferroni = FALSE,
                                      update = TRUE,   # FIX ME: per default set to NULL ?
                                      selectfun = "ctree", 
                                      splitfun = "ctree", 
                                      svselectfun = "ctree", 
                                      svsplitfun = "ctree"
)



ctest <- ctree(dist~speed, data = cars, ytrafo = ytrafo)
ctest <- extree_fit(data = d, trafo = ytrafo, ctrl = ctree_control(update = TRUE, splitflavour = "exhaustive"),
                    weights = rep.int(1, NROW(d)))
