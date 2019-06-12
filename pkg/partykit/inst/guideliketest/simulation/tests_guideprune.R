#################################
## Source file for replication of the simulation results presented in:
## The Power of Unbiased Recursive Partitioning: A Unifying View of CTree, MOB, and GUIDE (2019)
## by Lisa Schlosser and Torsten Hothorn and and Achim Zeileis
##
## Content: Data generating proceses and evaluation functions
#################################



# function to compute adjusted Rand Index
adj_rand_index <- function(x, y) {
  
  tab <- table(x, y)
  a <- rowSums(tab)
  b <- colSums(tab)
  
  M <- sum(choose(tab, 2))
  N <- choose(length(x), 2)
  A <- sum(choose(a, 2))
  B <- sum(choose(b, 2))
  
  c(ARI = (M - (A * B) / N) / (0.5 * (A + B) - (A * B) / N))
}



###########
# DGPs

dgp_stump <- function(nobs = 100, delta = 1, xi = 0.3,  
                      sigma = 1, seed = 7, 
                      changetype = c("abrupt", "continuous"),
                      variation = c("all", "beta0", "beta1"),
                      beta0 = NULL, beta1 = NULL)
{
  
  if(length(xi)>1) stop("for stump scenario only one split point xi can be handed over")
  if(length(variation) > 1) variation <- variation[1]
  if(!variation %in% c("all", "beta0", "beta1")) stop("variation has to be one of the following options: 'all', 'beta0', 'beta1'")
  if(variation == "beta0" & is.null(beta1)) stop("if only beta0 varies, a fixed value has to set for beta1")
  if(variation == "beta1" & is.null(beta0)) stop("if only beta1 varies, a fixed value has to set for beta0")
  if(length(changetype) > 1) changetype <- changetype[1]
  
  set.seed(seed)
  
  z1 <- runif(nobs,-1,1)
  z2 <- rnorm(nobs, 0, 1) 
  z3 <- rnorm(nobs, 0, 1)
  z4 <- rnorm(nobs, 0, 1)
  z5 <- rnorm(nobs, 0, 1)
  z6 <- rnorm(nobs, 0, 1)
  z7 <- runif(nobs, -1, 1)
  z8 <- runif(nobs, -1, 1)
  z9 <- runif(nobs, -1, 1)
  z10 <- runif(nobs, -1, 1)
  
  id <- numeric(length(z1))
  
  x <- runif(nobs, min = -1, max = 1)
  
  
  # for changetype abrupt: one step at break point xi, 
  # for changetype continuous: linear function
  if(variation == "all"){
    if(changetype == "abrupt"){
      beta0 <- delta * (-1)^(z1<xi)
      beta1 <- delta * (-1)^(z1<xi) * (-1)   # opposite signs  for beta0 and beta1   
    } else {
      beta0 <- delta * z1
      beta1 <- delta * z1 * (-1)   # opposite signs for beta0 and beta1
    }
  }
  
  if(variation == "beta0"){
    beta0 <- if(changetype == "abrupt") delta * (-1)^(z1<xi) else delta * z1
    beta1 <- beta1
  }
  
  if(variation == "beta1"){
    beta0 <- beta0
    beta1 <- if(changetype == "abrupt") delta * (-1)^(z1<xi) * (-1) else delta * z1 * (-1) 
    # opposite signs for beta0 and beta1
  }
  
  if(changetype == "abrupt") id <- 1+(z1>=xi)
  
  mu <- beta0 + beta1 * x
  
  y <- rnorm(nobs, mu, sigma)
  
  d <- data.frame(y = y, x = x, 
                  z1 = z1, z2 = z2, z3 = z3, z4 = z4, z5 = z5, z6 = z6, z7 = z7, z8 = z8, z9 = z9, z10 = z10,
                  beta0 = beta0, beta1 = beta1, mu = mu, sigma = rep.int(sigma, times = length(y)), id = id)
  
  return(d)
}



dgp_tree <- function(nobs = 100, delta = 1, xi = c(-0.3, 0.3),  
                     sigma = 1, seed = 7,  
                     changetype = "abrupt",
                     variation = "all",
                     beta0 = NULL, beta1 = NULL)
{
  # check input values
  if(variation != "all") stop("variation can only be set to 'all' in dgp_tree")
  if(changetype != "abrupt") stop("changetype can only be abrupt in dgp_tree")
  if(!is.null(beta0) | !is.null(beta1)) warning("values for beta0 or beta1 are ignored since variation='all' for dgp_tree")
  
  set.seed(seed)
  
  if(length(xi)==1){
    xi1 <- xi2 <- xi
  } else {
    xi1 <- xi[1]
    xi2 <- xi[2]
  }
  
  z1 <- runif(nobs,-1,1)
  z2 <- runif(nobs,-1,1) 
  z3 <- rnorm(nobs, 0, 1)
  z4 <- rnorm(nobs, 0, 1)
  z5 <- rnorm(nobs, 0, 1)
  z6 <- rnorm(nobs, 0, 1)
  z7 <- rnorm(nobs, 0, 1)
  z8 <- runif(nobs, -1, 1)
  z9 <- runif(nobs, -1, 1)
  z10 <- runif(nobs, -1, 1)
  
  id <- numeric(length(z1))
  
  x <- runif(nobs, min = -1, max = 1)
  
  ## first a split in z2 at split point xi2 => beta1 changes
  ## then a split in z1 at split point xi1 => beta0 changes
  # resulting subgroups:
  # group1 (z2<xi2):           beta0 = 0,      beta1 = delta
  # group2 (z2<xi2 & z1<xi1):  beta0 = -delta, beta1 = -delta
  # group3 (z2<xi2 & z1=>xi1): beta0 = +delta, beta1 = -delta
  
  #if(z2<xi2){
  #  beta0 <- 0 
  #  beta1 <- +delta
  #} else {
  #  if(z1<xi1){
  #    beta0 <- -delta 
  #    beta1 <- -delta
  #  } else {
  #    beta0 <- +delta 
  #    beta1 <- +delta
  #  }
  #}
  
  beta0 <- delta * (-1)^(z1<xi1) * 0^(z2<xi2)
  #beta0 <- delta * (-1)^(z1<xi1) * (z2>=xi2) + delta * (z2<xi2)
  beta1 <- delta * (-1)^(z2>=xi2)
  
  id <- 1 + (z2>=xi2) + (z2>=xi2)*(z1>=xi1)
  
  mu <- beta0 + beta1 * x
  
  y <- rnorm(nobs, mu, sigma)
  
  d <- data.frame(y = y, x = x, 
                  z1 = z1, z2 = z2, z3 = z3, z4 = z4, z5 = z5, z6 = z6, z7 = z7, z8 = z8, z9 = z9, z10 = z10,
                  beta0 = beta0, beta1 = beta1, mu = mu, sigma = rep.int(sigma, times = length(y)), id = id)
  
  return(d)
}



###########
# evaluation of one test on one given dataset
evaltests <- function(formula, data, 
                      scenario = c("stump", "tree"), 
                      testfun = c("guide", "ctree", "mfluc", "ctree_cat", "mfluc_cat", 
                                  "ctree_bin", "mfluc_bin", "ctree_cat_bin", "mfluc_cat_bin"), 
                      subset, weights, alpha = 0.05,
                      whichscores = c("all", "residuals"),
                      ctree_max = c(FALSE, TRUE),
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
  if(testfun == "guide" & length(whichscores) == 1) stop("for GUIDE use argument 'guide_parm' to decide which scores should be used")
  if(length(whichscores) > 1) whichscores <- whichscores[1]
  if(length(ctree_max) > 1) ctree_max <- ctree_max[1]
  splittest <- ctree_max
  if(length(scenario)>1) scenario <- scenario[1]
  if(length(testfun) > 1) testfun <- testfun[1]
  if(!testfun %in% c("guide", "ctree", "mfluc", "ctree_cat", "mfluc_cat", 
                     "ctree_bin", "mfluc_bin", "ctree_cat_bin", "mfluc_cat_bin")) 
    stop("testfun hast to be one of the following options: 'guide', 'ctree', 'mfluc', 'ctree_cat', 'mfluc_cat',
         'ctree_bin', 'mfluc_bin', 'ctree_cat_bin', 'mfluc_cat_bin'")
  
  
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
                                         bonferroni = TRUE,
                                         logmincriterion = log(1-alpha),
                                         update = TRUE,
                                         stump = scenario == "stump")
  }
  
  if(testfun == "ctree"){
    control <- partykit:::extree_control(criterion = "p.value",
                                         selectfun = partykit:::.ctree_select(teststat = "quadratic", splitstat = "quadratic", 
                                                                              splittest = splittest, pargs = GenzBretz(),
                                                                              testtype = "MonteCarlo", nresample = 9999L, 
                                                                              tol = sqrt(.Machine$double.eps),
                                                                              intersplit = TRUE, MIA = FALSE),
                                         splitfun = partykit:::.objfun_split(restart = TRUE,
                                                                             intersplit = TRUE),
                                         svselectfun = partykit:::.ctree_select(teststat = "quadratic", splitstat = "quadratic", 
                                                                                splittest = splittest, pargs = GenzBretz(),
                                                                                testtype = "MonteCarlo", nresample = 9999L, 
                                                                                tol = sqrt(.Machine$double.eps),
                                                                                intersplit = TRUE, MIA = FALSE), 
                                         svsplitfun = partykit:::.objfun_split(restart = TRUE,
                                                                               intersplit = TRUE),
                                         logmincriterion = log(1-alpha),
                                         update = TRUE,
                                         bonferroni = TRUE,
                                         stump = scenario == "stump")
  }
  
  if(testfun == "ctree_cat"){
    control <- partykit:::extree_control(criterion = "p.value",
                                         selectfun = .ctree_cat_select(teststat = "quadratic", splitstat = "quadratic", 
                                                                       splittest = splittest, pargs = GenzBretz(),
                                                                       testtype = "MonteCarlo", nresample = 9999L, 
                                                                       tol = sqrt(.Machine$double.eps),
                                                                       intersplit = TRUE, MIA = FALSE,
                                                                       xgroups = xgroups),
                                         splitfun = partykit:::.objfun_split(restart = TRUE,
                                                                             intersplit = TRUE),
                                         svselectfun = .ctree_cat_select(teststat = "quadratic", splitstat = "quadratic", 
                                                                         splittest = splittest, pargs = GenzBretz(),
                                                                         testtype = "MonteCarlo", nresample = 9999L, 
                                                                         tol = sqrt(.Machine$double.eps),
                                                                         intersplit = TRUE, MIA = FALSE,
                                                                         xgroups = xgroups), 
                                         svsplitfun = partykit:::.objfun_split(restart = TRUE,
                                                                               intersplit = TRUE),
                                         logmincriterion = log(1-alpha),
                                         update = TRUE,
                                         bonferroni = TRUE,
                                         stump = scenario == "stump")
  }
  
  if(testfun == "ctree_bin"){
    control <- partykit:::extree_control(criterion = "p.value",
                                         selectfun = .ctree_bin_select(teststat = "quadratic", splitstat = "quadratic", 
                                                                       splittest = splittest, pargs = GenzBretz(),
                                                                       testtype = "MonteCarlo", nresample = 9999L, 
                                                                       tol = sqrt(.Machine$double.eps),
                                                                       intersplit = TRUE, MIA = FALSE,
                                                                       xgroups = xgroups),
                                         splitfun = partykit:::.objfun_split(restart = TRUE,
                                                                             intersplit = TRUE),
                                         svselectfun = .ctree_bin_select(teststat = "quadratic", splitstat = "quadratic", 
                                                                         splittest = splittest, pargs = GenzBretz(),
                                                                         testtype = "MonteCarlo", nresample = 9999L, 
                                                                         tol = sqrt(.Machine$double.eps),
                                                                         intersplit = TRUE, MIA = FALSE,
                                                                         xgroups = xgroups), 
                                         svsplitfun = partykit:::.objfun_split(restart = TRUE,
                                                                               intersplit = TRUE),
                                         logmincriterion = log(1-alpha),
                                         update = TRUE,
                                         bonferroni = TRUE,
                                         stump = scenario == "stump")
  }
  
  if(testfun == "ctree_cat_bin"){
    control <- partykit:::extree_control(criterion = "p.value",
                                         selectfun = .ctree_cat_bin_select(teststat = "quadratic", splitstat = "quadratic", 
                                                                           splittest = splittest, pargs = GenzBretz(),
                                                                           testtype = "MonteCarlo", nresample = 9999L, 
                                                                           tol = sqrt(.Machine$double.eps),
                                                                           intersplit = TRUE, MIA = FALSE,
                                                                           xgroups = xgroups),
                                         splitfun = partykit:::.objfun_split(restart = TRUE,
                                                                             intersplit = TRUE),
                                         svselectfun = .ctree_cat_bin_select(teststat = "quadratic", splitstat = "quadratic", 
                                                                             splittest = splittest, pargs = GenzBretz(),
                                                                             testtype = "MonteCarlo", nresample = 9999L, 
                                                                             tol = sqrt(.Machine$double.eps),
                                                                             intersplit = TRUE, MIA = FALSE,
                                                                             xgroups = xgroups), 
                                         svsplitfun = partykit:::.objfun_split(restart = TRUE,
                                                                               intersplit = TRUE),
                                         logmincriterion = log(1-alpha),
                                         update = TRUE,
                                         bonferroni = TRUE,
                                         stump = scenario == "stump")
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
                                         logmincriterion = log(1-alpha),
                                         update = TRUE,
                                         bonferroni = TRUE,
                                         stump = scenario == "stump")
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
                                         logmincriterion = log(1-alpha),
                                         update = TRUE,
                                         bonferroni = TRUE,
                                         stump = scenario == "stump")
  }
  
  if(testfun == "mfluc_bin"){
    control <- partykit:::extree_control(criterion = "p.value",
                                         selectfun = .mfluc_bin_select(breakties = FALSE, 
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
                                         svselectfun = .mfluc_bin_select(breakties = FALSE, 
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
                                         logmincriterion = log(1-alpha),
                                         update = TRUE,
                                         bonferroni = TRUE,
                                         stump = scenario == "stump")
  }
  
  if(testfun == "mfluc_cat_bin"){
    control <- partykit:::extree_control(criterion = "p.value",
                                         selectfun = .mfluc_cat_bin_select(breakties = FALSE, 
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
                                         svselectfun = .mfluc_cat_bin_select(breakties = FALSE, 
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
                                         logmincriterion = log(1-alpha),
                                         update = TRUE,
                                         bonferroni = TRUE,
                                         stump = scenario == "stump")
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
    
    lmformula <- if(is.null(d$terms$yx)) Formula(d$terms$all) else Formula(d$terms$yx)   
    sdata <- d$data[subset,]
    
    model <- lm(formula = lmformula, data = sdata)
    
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
      
      if(whichscores == "residuals") estfun <- as.matrix(estfun[,1], ncol = 1)
      
    } else estfun <- NULL
    
    object <-  if(object) model else NULL
    
    ret <- list(estfun = estfun,
                unweighted = TRUE, # unweighted = TRUE would prevent estfun / w in extree_fit
                coefficients = coef(model),
                objfun = -logLik(model),  # optional function to be minimized 
                object = object,
                decorrelated = (decorrelate != "none"),    
                converged = !is.null(model)  
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
  control$xtype <- "matrix" 
  
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
  
}



###########
# evaluation of one test on various data sets (generated based on input arguments)
sim <- function(scenario = c("stump", "tree"),
                nobs = 100, nrep = 100, seed = 7,  
                formula = y~x|z1+z2+z3+z4+z5+z6+z7+z8+z9+z10,
                variation = "all", xi = 0, delta = 1, 
                changetype = c("abrupt", "continuous"),
                return_matrices = FALSE,
                compare_pruning = FALSE,
                test = c("ctree", "mfluc", "ctree_max", 
                         "ctree_cat", "ctree_max_cat", "mfluc_cat",
                         "ctree_bin", "ctree_max_bin", "mfluc_bin",
                         "ctree_cat_bin", "ctree_max_cat_bin", "mfluc_cat_bin",
                         "ctree_resid", "mfluc_resid", "ctree_resid_max", 
                         "ctree_resid_cat", "ctree_resid_max_cat", "mfluc_resid_cat",
                         "ctree_resid_bin", "ctree_resid_max_bin", "mfluc_resid_bin",
                         "ctree_resid_cat_bin", "ctree_resid_max_cat_bin", "mfluc_resid_cat_bin",
                         "guide_sum_12", "guide_max_12", "guide_coin_12",
                         "guide_sum_1", "guide_max_1", "guide_coin_1",
                         "guide_sum_2", "guide_max_2", "guide_coin_2",
                         "guide_sum_1_cor", "guide_max_1_cor", "guide_coin_1_cor",
                         "guide_sum_2_cor", "guide_max_2_cor", "guide_coin_2_cor"),
                beta0 = 0, beta1 = 1, 
                sigma = 1, alpha = 0.05)
{
  set.seed(seed)
  
  ## call
  cl <- match.call()
  
  ## check input arguments for comparing different pruning strategies
  if(compare_pruning){
    if(scenario == "stump") stop("different pruning strategies can only be compared for stump = FALSE")
  }
  
  if(length(scenario)>1) scenario <- scenario[1]
  if(scenario == "tree") beta0 <- beta1 <- NULL
  if(length(changetype) > 1) changetype <- changetype[1]

  
  if(scenario == "stump"){
    pval_z1 <- matrix(rep(NA, length(test)*nrep), ncol = length(test))
    pval <- matrix(rep(NA, length(test)*nrep), ncol = length(test))
    sv <- matrix(rep(NA, length(test)*nrep), ncol = length(test))
    colnames(pval) <- colnames(pval_z1) <- colnames(sv) <- test
  } else {
    nrsubgr <- matrix(rep(NA, length(test)*nrep), ncol = length(test))
    colnames(nrsubgr) <- test
  }
  err_coef <- matrix(rep(NA, length(test)*nrep), ncol = length(test))
  err_total <- matrix(rep(NA, length(test)*nrep), ncol = length(test))
  ari <- matrix(rep(NA, length(test)*nrep), ncol = length(test))
  colnames(err_coef) <- colnames(err_total) <- colnames(ari) <- test
  
  if(compare_pruning) {
    if(scenario == "stump"){
      pval_z1_p <- matrix(rep(NA, length(test)*nrep), ncol = length(test))
      pval_p <- matrix(rep(NA, length(test)*nrep), ncol = length(test))
      sv_p <- matrix(rep(NA, length(test)*nrep), ncol = length(test))
      colnames(pval_p) <- colnames(pval_z1_p) <- colnames(sv_p) <- test
    } else {
      nrsubgr_p <- matrix(rep(NA, length(test)*nrep), ncol = length(test))
      colnames(nrsubgr_p) <- test
    }
    err_coef_p <- matrix(rep(NA, length(test)*nrep), ncol = length(test))
    err_total_p <- matrix(rep(NA, length(test)*nrep), ncol = length(test))
    ari_p <- matrix(rep(NA, length(test)*nrep), ncol = length(test))
    colnames(err_coef_p) <- colnames(err_total_p) <- colnames(ari_p) <- test
  }  
  
  
  for(i in 1:nrep){
    
    dgp <- if(scenario == "stump") dgp_stump else dgp_tree
    d <- dgp(nobs = nobs, variation = variation, 
             beta0 = beta0, beta1 = beta1, xi = xi, 
             delta = delta,
             changetype = changetype, 
             seed = seed + i, sigma = sigma)
    
    
    compute_pval <- function(test) {
      test <- match.arg(test, c("ctree", "mfluc", "ctree_max", 
                                "ctree_cat", "ctree_max_cat", "mfluc_cat",
                                "ctree_bin", "ctree_max_bin", "mfluc_bin",
                                "ctree_cat_bin", "ctree_max_cat_bin", "mfluc_cat_bin",
                                "ctree_resid", "mfluc_resid", "ctree_resid_max", 
                                "ctree_resid_cat", "ctree_resid_max_cat", "mfluc_resid_cat",
                                "ctree_resid_bin", "ctree_resid_max_bin", "mfluc_resid_bin",
                                "ctree_resid_cat_bin", "ctree_resid_max_cat_bin", "mfluc_resid_cat_bin",
                                "guide_sum_12", "guide_max_12", "guide_coin_12",
                                "guide_sum_1", "guide_max_1", "guide_coin_1",
                                "guide_sum_2", "guide_max_2", "guide_coin_2",
                                "guide_sum_1_cor", "guide_max_1_cor", "guide_coin_1_cor",
                                "guide_sum_2_cor", "guide_max_2_cor", "guide_coin_2_cor"))
      testres <- switch(test,
                        "ctree" = evaltests(formula, data = d, testfun = "ctree", scenario = scenario, decorrelate = "vcov"),
                        "ctree_max" = evaltests(formula, data = d, testfun = "ctree", scenario = scenario, decorrelate = "vcov", ctree_max = TRUE),
                        "mfluc" = evaltests(formula, data = d, testfun = "mfluc", scenario = scenario, decorrelate = "vcov"),
                        
                        "ctree_cat" = evaltests(formula, data = d, testfun = "ctree_cat", scenario = scenario, decorrelate = "vcov"),
                        "ctree_max_cat" = evaltests(formula, data = d, testfun = "ctree_cat", scenario = scenario, decorrelate = "vcov", ctree_max = TRUE),
                        "mfluc_cat" = evaltests(formula, data = d, testfun = "mfluc_cat", scenario = scenario, decorrelate = "vcov"),
                        
                        "ctree_bin" = evaltests(formula, data = d, testfun = "ctree_bin", scenario = scenario, decorrelate = "vcov"),
                        "ctree_max_bin" = evaltests(formula, data = d, testfun = "ctree_bin", scenario = scenario, decorrelate = "vcov", ctree_max = TRUE),
                        "mfluc_bin" = evaltests(formula, data = d, testfun = "mfluc_bin", scenario = scenario, decorrelate = "vcov"),
                        
                        "ctree_cat_bin" = evaltests(formula, data = d, testfun = "ctree_cat_bin", scenario = scenario, decorrelate = "vcov"),
                        "ctree_max_cat_bin" = evaltests(formula, data = d, testfun = "ctree_cat_bin", scenario = scenario, decorrelate = "vcov", ctree_max = TRUE),
                        "mfluc_cat_bin" = evaltests(formula, data = d, testfun = "mfluc_cat_bin", scenario = scenario, decorrelate = "vcov"),
                        
                        "ctree_resid" = evaltests(formula, data = d, testfun = "ctree", scenario = scenario, decorrelate = "vcov", whichscores = "residuals"),
                        "ctree_resid_max" = evaltests(formula, data = d, testfun = "ctree", scenario = scenario, decorrelate = "vcov", ctree_max = TRUE, whichscores = "residuals"),
                        "mfluc_resid" = evaltests(formula, data = d, testfun = "mfluc", scenario = scenario, decorrelate = "vcov", whichscores = "residuals"),
                        
                        "ctree_resid_cat" = evaltests(formula, data = d, testfun = "ctree_cat", scenario = scenario, decorrelate = "vcov", whichscores = "residuals"),
                        "ctree_resid_max_cat" = evaltests(formula, data = d, testfun = "ctree_cat", scenario = scenario, decorrelate = "vcov", ctree_max = TRUE, whichscores = "residuals"),
                        "mfluc_resid_cat" = evaltests(formula, data = d, testfun = "mfluc_cat", scenario = scenario, decorrelate = "vcov", whichscores = "residuals"),
                        
                        "ctree_resid_bin" = evaltests(formula, data = d, testfun = "ctree_bin", scenario = scenario, decorrelate = "vcov", whichscores = "residuals"),
                        "ctree_resid_max_bin" = evaltests(formula, data = d, testfun = "ctree_bin", scenario = scenario, decorrelate = "vcov", ctree_max = TRUE, whichscores = "residuals"),
                        "mfluc_resid_bin" = evaltests(formula, data = d, testfun = "mfluc_bin", scenario = scenario, decorrelate = "vcov", whichscores = "residuals"),
                        
                        "ctree_resid_cat_bin" = evaltests(formula, data = d, testfun = "ctree_cat_bin", scenario = scenario, decorrelate = "vcov", whichscores = "residuals"),
                        "ctree_resid_max_cat_bin" = evaltests(formula, data = d, testfun = "ctree_cat_bin", scenario = scenario, decorrelate = "vcov", ctree_max = TRUE, whichscores = "residuals"),
                        "mfluc_resid_cat_bin" = evaltests(formula, data = d, testfun = "mfluc_cat_bin", scenario = scenario, decorrelate = "vcov", whichscores = "residuals"),
                        
                        "guide_sum_12" = evaltests(formula, data = d, testfun = "guide", 
                                                   guide_testtype = "sum", guide_parm = c(1,2),
                                                   xgroups = NULL, ygroups = NULL, scenario = scenario, decorrelate = "vcov"),
                        "guide_max_12" = evaltests(formula, data = d, testfun = "guide", 
                                                   guide_testtype = "max", guide_parm = c(1,2),
                                                   xgroups = NULL, ygroups = NULL, scenario = scenario, decorrelate = "vcov"),
                        "guide_coin_12" = evaltests(formula, data = d, testfun = "guide", 
                                                    guide_testtype = "coin", guide_parm = c(1,2),
                                                    xgroups = NULL, ygroups = NULL, scenario = scenario, decorrelate = "vcov"),
                        "guide_sum_1" = evaltests(formula, data = d, testfun = "guide", 
                                                  guide_testtype = "sum", guide_parm = c(1),
                                                  xgroups = NULL, ygroups = NULL, scenario = scenario, decorrelate = "vcov"),
                        "guide_max_1" = evaltests(formula, data = d, testfun = "guide", 
                                                  guide_testtype = "max", guide_parm = c(1),
                                                  xgroups = NULL, ygroups = NULL, scenario = scenario, decorrelate = "vcov"),
                        "guide_coin_1" = evaltests(formula, data = d, testfun = "guide", 
                                                   guide_testtype = "coin", guide_parm = c(1),
                                                   xgroups = NULL, ygroups = NULL, scenario = scenario, decorrelate = "vcov"),
                        "guide_sum_2" = evaltests(formula, data = d, testfun = "guide", 
                                                  guide_testtype = "sum", guide_parm = c(2),
                                                  xgroups = NULL, ygroups = NULL, scenario = scenario, decorrelate = "vcov"),
                        "guide_max_2" = evaltests(formula, data = d, testfun = "guide", 
                                                  guide_testtype = "max", guide_parm = c(2),
                                                  xgroups = NULL, ygroups = NULL, scenario = scenario, decorrelate = "vcov"),
                        "guide_coin_2" = evaltests(formula, data = d, testfun = "guide", 
                                                   guide_testtype = "coin", guide_parm = c(2),
                                                   xgroups = NULL, ygroups = NULL, scenario = scenario, decorrelate = "vcov"),
                        "guide_sum_1_cor" = evaltests(formula, data = d, testfun = "guide", 
                                                      guide_testtype = "sum", guide_parm = c(1),
                                                      xgroups = NULL, ygroups = NULL, scenario = scenario, decorrelate = "none"),
                        "guide_max_1_cor" = evaltests(formula, data = d, testfun = "guide", 
                                                      guide_testtype = "max", guide_parm = c(1),
                                                      xgroups = NULL, ygroups = NULL, scenario = scenario, decorrelate = "none"),
                        "guide_coin_1_cor" = evaltests(formula, data = d, testfun = "guide", 
                                                       guide_testtype = "coin", guide_parm = c(1),
                                                       xgroups = NULL, ygroups = NULL, scenario = scenario, decorrelate = "none"),
                        "guide_sum_2_cor" = evaltests(formula, data = d, testfun = "guide", 
                                                      guide_testtype = "sum", guide_parm = c(2),
                                                      xgroups = NULL, ygroups = NULL, scenario = scenario, decorrelate = "none"),
                        "guide_max_2_cor" = evaltests(formula, data = d, testfun = "guide", 
                                                      guide_testtype = "max", guide_parm = c(2),
                                                      xgroups = NULL, ygroups = NULL, scenario = scenario, decorrelate = "none"),
                        "guide_coin_2_cor" = evaltests(formula, data = d, testfun = "guide", 
                                                       guide_testtype = "coin", guide_parm = c(2),
                                                       xgroups = NULL, ygroups = NULL, scenario = scenario, decorrelate = "none")
      )
      
      predresp <- predict(testres, type = "response")
      prednode <- predict(testres, type = "node")
      predcoef <- if(is.matrix(coef(testres))) coef(testres)[as.character(prednode),] else coef(testres)
      
      # calculate error of predicted coefficients by taking the euclidean norm of the differences of predicted and true coefficients,
      # averaged over all observations
      err_coef <- mean(sqrt(rowSums((predcoef - d[,c("beta0", "beta1")])^2)))
      # mean squared error of predicted responses
      err_total <- mean(sqrt((predresp - d[,"y"])^2))
      
      ari <- as.numeric(adj_rand_index(prednode, d$id))
      
      if(scenario == "stump"){
        returnlist <- list(pval = as.numeric(info_node(testres$node)$p.value),
                           pval_z1 = as.numeric(info_node(testres$node)$criterion["p.value","z1"]),
                           sv = names(info_node(testres$node)$p.value),
                           err_coef = err_coef,
                           err_total = err_total,
                           ari = ari)
      } else {
        returnlist <- list(err_coef = err_coef,
                           err_total = err_total,
                           nrsubgr = width(testres),
                           ari = ari)
      }
      
      
      ## if the same tree should be built using post-pruning
      if(compare_pruning){
        testres_p <- switch(test,
                            "ctree" = evaltests(formula, data = d, testfun = "ctree", scenario = scenario, decorrelate = "vcov", alpha = 1),
                            "ctree_max" = evaltests(formula, data = d, testfun = "ctree", scenario = scenario, decorrelate = "vcov", ctree_max = TRUE, alpha = 1),
                            "mfluc" = evaltests(formula, data = d, testfun = "mfluc", scenario = scenario, decorrelate = "vcov", alpha = 1),
                            
                            "ctree_cat" = evaltests(formula, data = d, testfun = "ctree_cat", scenario = scenario, decorrelate = "vcov", alpha = 1),
                            "ctree_max_cat" = evaltests(formula, data = d, testfun = "ctree_cat", scenario = scenario, decorrelate = "vcov", ctree_max = TRUE, alpha = 1),
                            "mfluc_cat" = evaltests(formula, data = d, testfun = "mfluc_cat", scenario = scenario, decorrelate = "vcov", alpha = 1),
                            
                            "ctree_bin" = evaltests(formula, data = d, testfun = "ctree_bin", scenario = scenario, decorrelate = "vcov", alpha = 1),
                            "ctree_max_bin" = evaltests(formula, data = d, testfun = "ctree_bin", scenario = scenario, decorrelate = "vcov", ctree_max = TRUE, alpha = 1),
                            "mfluc_bin" = evaltests(formula, data = d, testfun = "mfluc_bin", scenario = scenario, decorrelate = "vcov", alpha = 1),
                            
                            "ctree_cat_bin" = evaltests(formula, data = d, testfun = "ctree_cat_bin", scenario = scenario, decorrelate = "vcov", alpha = 1),
                            "ctree_max_cat_bin" = evaltests(formula, data = d, testfun = "ctree_cat_bin", scenario = scenario, decorrelate = "vcov", ctree_max = TRUE, alpha = 1),
                            "mfluc_cat_bin" = evaltests(formula, data = d, testfun = "mfluc_cat_bin", scenario = scenario, decorrelate = "vcov", alpha = 1),
                            
                            "ctree_resid" = evaltests(formula, data = d, testfun = "ctree", scenario = scenario, decorrelate = "vcov", whichscores = "residuals", alpha = 1),
                            "ctree_resid_max" = evaltests(formula, data = d, testfun = "ctree", scenario = scenario, decorrelate = "vcov", ctree_max = TRUE, whichscores = "residuals", alpha = 1),
                            "mfluc_resid" = evaltests(formula, data = d, testfun = "mfluc", scenario = scenario, decorrelate = "vcov", whichscores = "residuals", alpha = 1),
                            
                            "ctree_resid_cat" = evaltests(formula, data = d, testfun = "ctree_cat", scenario = scenario, decorrelate = "vcov", whichscores = "residuals", alpha = 1),
                            "ctree_resid_max_cat" = evaltests(formula, data = d, testfun = "ctree_cat", scenario = scenario, decorrelate = "vcov", ctree_max = TRUE, whichscores = "residuals", alpha = 1),
                            "mfluc_resid_cat" = evaltests(formula, data = d, testfun = "mfluc_cat", scenario = scenario, decorrelate = "vcov", whichscores = "residuals", alpha = 1),
                            
                            "ctree_resid_bin" = evaltests(formula, data = d, testfun = "ctree_bin", scenario = scenario, decorrelate = "vcov", whichscores = "residuals", alpha = 1),
                            "ctree_resid_max_bin" = evaltests(formula, data = d, testfun = "ctree_bin", scenario = scenario, decorrelate = "vcov", ctree_max = TRUE, whichscores = "residuals", alpha = 1),
                            "mfluc_resid_bin" = evaltests(formula, data = d, testfun = "mfluc_bin", scenario = scenario, decorrelate = "vcov", whichscores = "residuals", alpha = 1),
                            
                            "ctree_resid_cat_bin" = evaltests(formula, data = d, testfun = "ctree_cat_bin", scenario = scenario, decorrelate = "vcov", whichscores = "residuals", alpha = 1),
                            "ctree_resid_max_cat_bin" = evaltests(formula, data = d, testfun = "ctree_cat_bin", scenario = scenario, decorrelate = "vcov", ctree_max = TRUE, whichscores = "residuals", alpha = 1),
                            "mfluc_resid_cat_bin" = evaltests(formula, data = d, testfun = "mfluc_cat_bin", scenario = scenario, decorrelate = "vcov", whichscores = "residuals", alpha = 1),
                            
                            "guide_sum_12" = evaltests(formula, data = d, testfun = "guide", 
                                                       guide_testtype = "sum", guide_parm = c(1,2),
                                                       xgroups = NULL, ygroups = NULL, scenario = scenario, decorrelate = "vcov", alpha = 1),
                            "guide_max_12" = evaltests(formula, data = d, testfun = "guide", 
                                                       guide_testtype = "max", guide_parm = c(1,2),
                                                       xgroups = NULL, ygroups = NULL, scenario = scenario, decorrelate = "vcov", alpha = 1),
                            "guide_coin_12" = evaltests(formula, data = d, testfun = "guide", 
                                                        guide_testtype = "coin", guide_parm = c(1,2),
                                                        xgroups = NULL, ygroups = NULL, scenario = scenario, decorrelate = "vcov", alpha = 1),
                            "guide_sum_1" = evaltests(formula, data = d, testfun = "guide", 
                                                      guide_testtype = "sum", guide_parm = c(1),
                                                      xgroups = NULL, ygroups = NULL, scenario = scenario, decorrelate = "vcov", alpha = 1),
                            "guide_max_1" = evaltests(formula, data = d, testfun = "guide", 
                                                      guide_testtype = "max", guide_parm = c(1),
                                                      xgroups = NULL, ygroups = NULL, scenario = scenario, decorrelate = "vcov", alpha = 1),
                            "guide_coin_1" = evaltests(formula, data = d, testfun = "guide", 
                                                       guide_testtype = "coin", guide_parm = c(1),
                                                       xgroups = NULL, ygroups = NULL, scenario = scenario, decorrelate = "vcov", alpha = 1),
                            "guide_sum_2" = evaltests(formula, data = d, testfun = "guide", 
                                                      guide_testtype = "sum", guide_parm = c(2),
                                                      xgroups = NULL, ygroups = NULL, scenario = scenario, decorrelate = "vcov", alpha = 1),
                            "guide_max_2" = evaltests(formula, data = d, testfun = "guide", 
                                                      guide_testtype = "max", guide_parm = c(2),
                                                      xgroups = NULL, ygroups = NULL, scenario = scenario, decorrelate = "vcov", alpha = 1),
                            "guide_coin_2" = evaltests(formula, data = d, testfun = "guide", 
                                                       guide_testtype = "coin", guide_parm = c(2),
                                                       xgroups = NULL, ygroups = NULL, scenario = scenario, decorrelate = "vcov", alpha = 1),
                            "guide_sum_1_cor" = evaltests(formula, data = d, testfun = "guide", 
                                                          guide_testtype = "sum", guide_parm = c(1),
                                                          xgroups = NULL, ygroups = NULL, scenario = scenario, decorrelate = "none", alpha = 1),
                            "guide_max_1_cor" = evaltests(formula, data = d, testfun = "guide", 
                                                          guide_testtype = "max", guide_parm = c(1),
                                                          xgroups = NULL, ygroups = NULL, scenario = scenario, decorrelate = "none", alpha = 1),
                            "guide_coin_1_cor" = evaltests(formula, data = d, testfun = "guide", 
                                                           guide_testtype = "coin", guide_parm = c(1),
                                                           xgroups = NULL, ygroups = NULL, scenario = scenario, decorrelate = "none", alpha = 1),
                            "guide_sum_2_cor" = evaltests(formula, data = d, testfun = "guide", 
                                                          guide_testtype = "sum", guide_parm = c(2),
                                                          xgroups = NULL, ygroups = NULL, scenario = scenario, decorrelate = "none", alpha = 1),
                            "guide_max_2_cor" = evaltests(formula, data = d, testfun = "guide", 
                                                          guide_testtype = "max", guide_parm = c(2),
                                                          xgroups = NULL, ygroups = NULL, scenario = scenario, decorrelate = "none", alpha = 1),
                            "guide_coin_2_cor" = evaltests(formula, data = d, testfun = "guide", 
                                                           guide_testtype = "coin", guide_parm = c(2),
                                                           xgroups = NULL, ygroups = NULL, scenario = scenario, decorrelate = "none", alpha = 1)
        )
        
        testres_p <- ccprune(testres_p, nrfolds = 5)$Tree
        
        prednode_p <- predict(testres_p, type = "node", newdata = testres_p$data)
        predcoef_p <- if(is.matrix(coef(testres_p))) coef(testres_p)[as.character(prednode_p),] else coef(testres_p)
        if(is.vector(predcoef_p)) predcoef_p <- matrix(predcoef_p, nrow = 1, ncol = 2)
        predresp_p <- predcoef_p[,1]  + predcoef_p[,2] * testres_p$data$x
        
        # calculate error of predicted coefficients by taking the euclidean norm of the differences of predicted and true coefficients,
        # averaged over all observations
        err_coef_p <- mean(sqrt(rowSums((predcoef_p - d[,c("beta0", "beta1")])^2)))
        # mean squared error of predicted responses
        err_total_p <- mean(sqrt((predresp_p - d[,"y"])^2))
        
        ari_p <- as.numeric(adj_rand_index(prednode_p, d$id))
        
        
        returnlist <- c(returnlist,
                        list(err_coef_p = err_coef_p,
                             err_total_p = err_total_p,
                             nrsubgr_p = width(testres_p),
                             ari_p = ari_p))
      }
      
      return(returnlist)
    }
    
    resmat <- sapply(test, compute_pval)
    
    if(scenario == "stump"){
      pval[i,] <- unlist(resmat["pval",])
      pval_z1[i,] <- unlist(resmat["pval_z1",])
      sv[i,] <- unlist(resmat["sv",])
    } else {
      nrsubgr[i,] <- unlist(resmat["nrsubgr",])
    }
    err_coef[i,] <- unlist(resmat["err_coef",])
    err_total[i,] <- unlist(resmat["err_total",])
    ari[i,] <- unlist(resmat["ari",])
    
    if(compare_pruning){
      
      nrsubgr_p[i,] <- unlist(resmat["nrsubgr_p",])
      
      err_coef_p[i,] <- unlist(resmat["err_coef_p",])
      err_total_p[i,] <- unlist(resmat["err_total_p",])
      ari_p[i,] <- unlist(resmat["ari_p",])
    }
    
    
  }
  
  
  if(scenario == "stump"){
    # proportion of cases where a split is / is not found
    prop_nosplit <- colMeans(pval > alpha)
    prop_split <- colMeans(pval <= alpha)
    
    # proportion of cases where the correct variable (z1) has the smallest p-value (but not necessarily smaller than alpha)
    prop_z1 <- colMeans(sv =="z1")
    
    # proportion of cases where a split is found in the correct variable (z1)
    pval_T <- pval
    pval_T[which(sv !="z1")] <- 1
    prop_Tsplit <- colMeans(pval_T < alpha)
    
    # proportion of cases where a split is found in a noise variable (z2, ..., z10)
    pval_F <- pval
    pval_F[which(sv =="z1")] <- 1
    prop_Fsplit <- colMeans(pval_F < alpha)
    
    # average of smallest p-value (regardless of the corresponding splitting variable)
    pval_min <- colMeans(pval)
    
    # average of p-value for the correct variable (z1)
    pval_true <- colMeans(pval_z1)
    
  } else {
    # average of number of subgroups
    nrsubgr <- colMeans(nrsubgr)
  }
  
  # average of errors and adjusted rand index
  err_coef <- colMeans(err_coef)
  err_total <- colMeans(err_total)
  ari <- colMeans(ari)
  
  if(scenario == "stump"){
    returnlist <- list(prop_nosplit = prop_nosplit,
                       prop_split = prop_split,
                       prop_Fsplit = prop_Fsplit,
                       prop_Tsplit = prop_Tsplit,
                       prop_z1 = prop_z1,
                       pval_min = pval_min,
                       pval_true = pval_true, 
                       err_coef = err_coef,
                       err_total = err_total,
                       ari = ari)
  } else {
    returnlist <- list(nrsubgr = nrsubgr,
                       err_coef = err_coef,
                       err_total = err_total,
                       ari = ari)
  }
  
  if(return_matrices) returnlist <- c(returnlist, list(sv = sv,
                                                       pval = pval,
                                                       pval_z1 = pval_z1))
  
  
  if(compare_pruning){
    
    # average of number of subgroups
    nrsubgr_p <- colMeans(nrsubgr_p)
    
    # average of errors and adjusted rand index
    err_coef_p <- colMeans(err_coef_p)
    err_total_p <- colMeans(err_total_p)
    ari_p <- colMeans(ari_p)
    
    returnlist <- c(returnlist, 
                    list(nrsubgr_p = nrsubgr_p,
                         err_coef_p = err_coef_p,
                         err_total_p = err_total_p,
                         ari_p = ari_p))
  }
  
  return(returnlist)
}




###########
## wrapper function applying sim with different tests and over varying variables of interest
simwrapper <- function(nobs = 200, nrep = 100, seed = 7, 
                       scenario = c("stump", "tree"),
                       delta = seq(from = 1, to = 5, by = 2),
                       xi = c(0, 0.8), 
                       variation = c("all", "beta0", "beta1"),
                       changetype = c("abrupt", "continuous"),
                       compare_pruning = FALSE,
                       test = c("ctree", "mfluc", "ctree_max", 
                                "ctree_cat", "ctree_max_cat", "mfluc_cat",
                                "ctree_bin", "ctree_max_bin", "mfluc_bin",
                                "ctree_cat_bin", "ctree_max_cat_bin", "mfluc_cat_bin",
                                "ctree_resid", "mfluc_resid", "ctree_resid_max", 
                                "ctree_resid_cat", "ctree_resid_max_cat", "mfluc_resid_cat",
                                "ctree_resid_bin", "ctree_resid_max_bin", "mfluc_resid_bin",
                                "ctree_resid_cat_bin", "ctree_resid_max_cat_bin", "mfluc_resid_cat_bin",
                                "guide_sum_12", "guide_max_12", "guide_coin_12",
                                "guide_sum_1", "guide_max_1", "guide_coin_1",
                                "guide_sum_2", "guide_max_2", "guide_coin_2",
                                "guide_sum_1_cor", "guide_max_1_cor", "guide_coin_1_cor",
                                "guide_sum_2_cor", "guide_max_2_cor", "guide_coin_2_cor"),
                       beta0 = 0, beta1 = 1, 
                       sigma = 1, alpha = 0.05)
{
  
  cl <- match.call()
  if(length(scenario)>1) scenario <- scenario[1]
  if(scenario == "tree") beta0 <- beta1 <- NULL
  if(scenario == "tree" & "continuous" %in% changetype) stop("changetype can only be abrupt for tree scenario")
  
  prs <- expand.grid(delta = delta, xi = xi, 
                     variation = variation,
                     changetype = changetype)
  
  rmid <- which(prs$changetype == "continuous" & prs$xi != 0)
  if(length(rmid) > 0) prs <- prs[-rmid,]
  
  rownames(prs) <- c(1:NROW(prs))
  
  nprs <- nrow(prs)
  ntest <- length(test)
  
  simres <- mclapply(1:nprs, 
                     function(i) {
                       reslist <- sim(nobs = nobs, nrep = nrep, scenario = scenario, 
                                      seed = seed,
                                      formula = y~x|z1+z2+z3+z4+z5+z6+z7+z8+z9+z10,
                                      beta0 = beta0, beta1 = beta1,
                                      sigma = sigma, alpha = alpha,
                                      test = test, 
                                      variation = prs$variation[i], 
                                      changetype = prs$changetype[i],
                                      xi = prs$xi[i],
                                      delta = prs$delta[i],
                                      compare_pruning = compare_pruning)
                       
                       
                       #prop_nosplit[i,] <- reslist$prop_nosplit
                       #prop_split[i,] <- reslist$prop_split
                       #prop_Fsplit[i,] <- reslist$prop_Fsplit
                       #prop_Tsplit[i,] <- reslist$prop_Tsplit
                       
                       return(reslist)
                     },
                     mc.cores = detectCores() - 1
  )
  
  
  if(scenario == "stump"){
    
    prop_nosplit <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    prop_split <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    prop_Fsplit <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    prop_Tsplit <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    prop_z1 <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    pval_min <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    pval_true <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    
    for(i in 1:nprs){
      prop_nosplit[i,] <- simres[[i]]$prop_nosplit
      prop_split[i,] <- simres[[i]]$prop_split
      prop_Fsplit[i,] <- simres[[i]]$prop_Fsplit
      prop_Tsplit[i,] <- simres[[i]]$prop_Tsplit
      prop_z1[i,] <- simres[[i]]$prop_z1
      pval_min[i,] <- simres[[i]]$pval_min
      pval_true[i,] <- simres[[i]]$pval_true
    }
    
  } else {
    nrsubgr <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    for(i in 1:nprs) nrsubgr[i,] <- simres[[i]]$nrsubgr
    
    if(compare_pruning){
      nrsubgr_p <- matrix(rep(NA, ntest * nprs), ncol = ntest)
      for(i in 1:nprs) nrsubgr_p[i,] <- simres[[i]]$nrsubgr_p
    }
  }
  
  err_coef <- matrix(rep(NA, ntest * nprs), ncol = ntest)
  err_total <- matrix(rep(NA, ntest * nprs), ncol = ntest)
  ari <- matrix(rep(NA, ntest * nprs), ncol = ntest)
  
  for(i in 1:nprs){
    err_coef[i,] <- simres[[i]]$err_coef
    err_total[i,] <- simres[[i]]$err_total
    ari[i,] <- simres[[i]]$ari
  }
  
  if(compare_pruning){
    err_coef_p <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    err_total_p <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    ari_p <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    
    for(i in 1:nprs){
      err_coef_p[i,] <- simres[[i]]$err_coef_p
      err_total_p[i,] <- simres[[i]]$err_total_p
      ari_p[i,] <- simres[[i]]$ari_p
    }
  }
  
  res <- data.frame()
  for(i in 1:ntest){res <- rbind(res, prs)}
  res$test <- gl(ntest, nprs, labels = test)
  
  
  if(scenario == "stump"){
    res$prop_nosplit <- as.vector(prop_nosplit)
    res$prop_split <- as.vector(prop_split)
    res$prop_Fsplit <- as.vector(prop_Fsplit)
    res$prop_Tsplit <- as.vector(prop_Tsplit)
    res$prop_z1 <- as.vector(prop_z1)
    res$pval_min <- as.vector(pval_min)
    res$pval_true <- as.vector(pval_true)
  } else {
    res$nrsubgr <- as.vector(nrsubgr)
    if(compare_pruning) res$nrsubgr_p <- as.vector(nrsubgr_p)
  }
  
  res$err_coef <- as.vector(err_coef)
  res$err_total <- as.vector(err_total)
  res$ari <- as.vector(ari)
  
  if(compare_pruning){
    res$err_coef_p <- as.vector(err_coef_p)
    res$err_total_p <- as.vector(err_total_p)
    res$ari_p <- as.vector(ari_p)
  }
  
  res$delta <- factor(res$delta)
  res$variation <- factor(res$variation)
  res$changetype <- factor(res$changetype)
  res$xi <- factor(res$xi)
  
  return(list(res = res,
              call = cl))
}



###########
## function to prepare data set for full factorial analysis 
## (restructured based on building blocks)
prep_3way <- function(simres)
{
  pval<- simres$pval
  # prepare data.frame
  {
    pval_T <- pval
    pval_T[which(simres$sv !="z1")] <- 1
    
    head(pval_T)
    
    d <- rep(colnames(pval)[1], NROW(pval))                 
    for(i in 2:NCOL(pval)){
      d <- c(d, rep(colnames(pval)[i], NROW(pval)))
    }
    
    d <- cbind(d, 
               rep(0, length(d)), 
               rep(0, length(d)), 
               rep(0, length(d)))    
    
    colnames(d) <- c("strategy", "res_scores", "bin", "cat")
    
    ## set values for factor variables representing the three building blocks
    {
      d[d[,"strategy"] == "ctree", "res_scores"] <- "scores"
      d[d[,"strategy"] == "ctree", "bin"] <- "lin"
      d[d[,"strategy"] == "ctree", "cat"] <- "lin"
      
      d[d[,"strategy"] == "ctree_max", "res_scores"] <- "scores"
      d[d[,"strategy"] == "ctree_max", "bin"] <- "lin"
      d[d[,"strategy"] == "ctree_max", "cat"] <- "max"
      
      d[d[,"strategy"] == "mfluc", "res_scores"] <- "scores"
      d[d[,"strategy"] == "mfluc", "bin"] <- "lin"
      d[d[,"strategy"] == "mfluc", "cat"] <- "max"
      
      d[d[,"strategy"] == "ctree_cat", "res_scores"] <- "scores"
      d[d[,"strategy"] == "ctree_cat", "bin"] <- "lin"
      d[d[,"strategy"] == "ctree_cat", "cat"] <- "cat"
      
      d[d[,"strategy"] == "ctree_max_cat", "res_scores"] <- "scores"
      d[d[,"strategy"] == "ctree_max_cat", "bin"] <- "lin"
      d[d[,"strategy"] == "ctree_max_cat", "cat"] <- "cat"
      
      d[d[,"strategy"] == "mfluc_cat", "res_scores"] <- "scores"
      d[d[,"strategy"] == "mfluc_cat", "bin"] <- "lin"
      d[d[,"strategy"] == "mfluc_cat", "cat"] <- "cat"
      
      d[d[,"strategy"] == "ctree_bin", "res_scores"] <- "scores"
      d[d[,"strategy"] == "ctree_bin", "bin"] <- "bin"
      d[d[,"strategy"] == "ctree_bin", "cat"] <- "lin"
      
      d[d[,"strategy"] == "ctree_max_bin", "res_scores"] <- "scores"
      d[d[,"strategy"] == "ctree_max_bin", "bin"] <- "bin"
      d[d[,"strategy"] == "ctree_max_bin", "cat"] <- "max"
      
      d[d[,"strategy"] == "mfluc_bin", "res_scores"] <- "scores"
      d[d[,"strategy"] == "mfluc_bin", "bin"] <- "bin"
      d[d[,"strategy"] == "mfluc_bin", "cat"] <- "max"
      
      d[d[,"strategy"] == "ctree_cat_bin", "res_scores"] <- "scores"
      d[d[,"strategy"] == "ctree_cat_bin", "bin"] <- "bin"
      d[d[,"strategy"] == "ctree_cat_bin", "cat"] <- "cat"
      
      d[d[,"strategy"] == "ctree_max_cat_bin", "res_scores"] <- "scores"
      d[d[,"strategy"] == "ctree_max_cat_bin", "bin"] <- "bin"
      d[d[,"strategy"] == "ctree_max_cat_bin", "cat"] <- "cat"
      
      d[d[,"strategy"] == "mfluc_cat_bin", "res_scores"] <- "scores"
      d[d[,"strategy"] == "mfluc_cat_bin", "bin"] <- "bin"
      d[d[,"strategy"] == "mfluc_cat_bin", "cat"] <- "cat"
      
      
      
      d[d[,"strategy"] == "ctree_resid", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "ctree_resid", "bin"] <- "lin"
      d[d[,"strategy"] == "ctree_resid", "cat"] <- "lin"
      
      d[d[,"strategy"] == "ctree_resid_max", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "ctree_resid_max", "bin"] <- "lin"
      d[d[,"strategy"] == "ctree_resid_max", "cat"] <- "max"
      
      d[d[,"strategy"] == "mfluc_resid", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "mfluc_resid", "bin"] <- "lin"
      d[d[,"strategy"] == "mfluc_resid", "cat"] <- "max"
      
      d[d[,"strategy"] == "ctree_resid_cat", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "ctree_resid_cat", "bin"] <- "lin"
      d[d[,"strategy"] == "ctree_resid_cat", "cat"] <- "cat"
      
      d[d[,"strategy"] == "ctree_resid_max_cat", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "ctree_resid_max_cat", "bin"] <- "lin"
      d[d[,"strategy"] == "ctree_resid_max_cat", "cat"] <- "cat"
      
      d[d[,"strategy"] == "mfluc_resid_cat", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "mfluc_resid_cat", "bin"] <- "lin"
      d[d[,"strategy"] == "mfluc_resid_cat", "cat"] <- "cat"
      
      d[d[,"strategy"] == "ctree_resid_bin", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "ctree_resid_bin", "bin"] <- "bin"
      d[d[,"strategy"] == "ctree_resid_bin", "cat"] <- "lin"
      
      d[d[,"strategy"] == "ctree_resid_max_bin", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "ctree_resid_max_bin", "bin"] <- "bin"
      d[d[,"strategy"] == "ctree_resid_max_bin", "cat"] <- "max"
      
      d[d[,"strategy"] == "mfluc_resid_bin", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "mfluc_resid_bin", "bin"] <- "bin"
      d[d[,"strategy"] == "mfluc_resid_bin", "cat"] <- "max"
      
      d[d[,"strategy"] == "ctree_resid_cat_bin", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "ctree_resid_cat_bin", "bin"] <- "bin"
      d[d[,"strategy"] == "ctree_resid_cat_bin", "cat"] <- "cat"
      
      d[d[,"strategy"] == "ctree_resid_max_cat_bin", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "ctree_resid_max_cat_bin", "bin"] <- "bin"
      d[d[,"strategy"] == "ctree_resid_max_cat_bin", "cat"] <- "cat"
      
      d[d[,"strategy"] == "mfluc_resid_cat_bin", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "mfluc_resid_cat_bin", "bin"] <- "bin"
      d[d[,"strategy"] == "mfluc_resid_cat_bin", "cat"] <- "cat"
      
      
      ## GUIDE
      
      d[d[,"strategy"] == "guide_sum_1", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "guide_sum_1", "bin"] <- "bin"
      d[d[,"strategy"] == "guide_sum_1", "cat"] <- "cat"
      
      d[d[,"strategy"] == "guide_coin_1", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "guide_coin_1", "bin"] <- "bin"
      d[d[,"strategy"] == "guide_coin_1", "cat"] <- "cat"
      
      d[d[,"strategy"] == "guide_max_1", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "guide_max_1", "bin"] <- "bin"
      d[d[,"strategy"] == "guide_max_1", "cat"] <- "cat"
      
      d[d[,"strategy"] == "guide_sum_1_cor", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "guide_sum_1_cor", "bin"] <- "bin"
      d[d[,"strategy"] == "guide_sum_1_cor", "cat"] <- "cat"
      
      d[d[,"strategy"] == "guide_coin_1_cor", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "guide_coin_1_cor", "bin"] <- "bin"
      d[d[,"strategy"] == "guide_coin_1_cor", "cat"] <- "cat"
      
      d[d[,"strategy"] == "guide_max_1_cor", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "guide_max_1_cor", "bin"] <- "bin"
      d[d[,"strategy"] == "guide_max_1_cor", "cat"] <- "cat"
      
      
      d[d[,"strategy"] == "guide_sum_2", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "guide_sum_2", "bin"] <- "bin"
      d[d[,"strategy"] == "guide_sum_2", "cat"] <- "cat"
      
      d[d[,"strategy"] == "guide_coin_2", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "guide_coin_2", "bin"] <- "bin"
      d[d[,"strategy"] == "guide_coin_2", "cat"] <- "cat"
      
      d[d[,"strategy"] == "guide_max_2", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "guide_max_2", "bin"] <- "bin"
      d[d[,"strategy"] == "guide_max_2", "cat"] <- "cat"
      
      d[d[,"strategy"] == "guide_sum_2_cor", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "guide_sum_2_cor", "bin"] <- "bin"
      d[d[,"strategy"] == "guide_sum_2_cor", "cat"] <- "cat"
      
      d[d[,"strategy"] == "guide_coin_2_cor", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "guide_coin_2_cor", "bin"] <- "bin"
      d[d[,"strategy"] == "guide_coin_2_cor", "cat"] <- "cat"
      
      d[d[,"strategy"] == "guide_max_2_cor", "res_scores"] <- "residuals"
      d[d[,"strategy"] == "guide_max_2_cor", "bin"] <- "bin"
      d[d[,"strategy"] == "guide_max_2_cor", "cat"] <- "cat"
      
      d[d[,"strategy"] == "guide_sum_12", "res_scores"] <- "scores"
      d[d[,"strategy"] == "guide_sum_12", "bin"] <- "bin"
      d[d[,"strategy"] == "guide_sum_12", "cat"] <- "cat"
      
      d[d[,"strategy"] == "guide_max_12", "res_scores"] <- "scores"
      d[d[,"strategy"] == "guide_max_12", "bin"] <- "bin"
      d[d[,"strategy"] == "guide_max_12", "cat"] <- "cat"
      
      d[d[,"strategy"] == "guide_coin_12", "res_scores"] <- "scores"
      d[d[,"strategy"] == "guide_coin_12", "bin"] <- "bin"
      d[d[,"strategy"] == "guide_coin_12", "cat"] <- "cat"
    }
    
    
    d <- data.frame(strategy = factor(d[,"strategy"]),
                    res_scores = factor(d[,"res_scores"]),
                    bin = factor(d[,"bin"]),
                    cat = factor(d[,"cat"]))
    
    d$pval <- as.vector(pval)
    d$pvalT <- as.vector(pval_T)
    d$pval_z1 <- as.vector(simres$pval_z1)
  }
  
  return(d)
}
