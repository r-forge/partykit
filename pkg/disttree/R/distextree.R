############
# new version of disttree using extree directly without applying ctree

distextree <- function(formula, data, subset, weights, family = NO(), na.action = na.pass, offset, cluster,
                       control = distextree_control(...), converged = NULL, scores = NULL, 
                       doFit = TRUE, bd = NULL,
                       #type.tree = "ctree",   ## should be in distextree_control 
                       decorrelate = "none",
                       censtype = "none", censpoint = NULL,
                       ocontrol = list(), ...) {
  
  ## keep call
  cl <- match.call(expand.dots = TRUE)
  if(missing(data)) data <- environment(formula)
  
  # check input arguments
  ## FIX ME: not necessary to set type of tree, but tests etc.
  # type.tree <- match.arg(type.tree, c("mob", "ctree"))
  # if(!(type.tree %in% c("mob", "ctree"))) stop("unknown argument for type.tree (can only be mob or ctree)")
  if(!(decorrelate) %in% c("none", "opg", "vcov")) stop("unknown argument for decorrelate (can only be none, opg or vcov)")
  # check formula
  oformula <- as.formula(formula)
  formula <- as.Formula(formula)
  if(length(formula)[2L] > 1L) {
    formula <- Formula(formula(formula, rhs = 2L))  
    ## FIX ME: if rhs has more than 1 element it is here assumed that partitioning variables are handed over on 2nd slot
    warning("formula must not have more than one RHS parts (only partitioning variables allowed)")
  }
  
  ## set up model.frame() call
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset", "cluster", "scores"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$yx <- "none"
  
  # function ytrafo is fixed within distextree (does not depend on x)
  #if (is.function(ytrafo)) {
  #  if (all(c("y", "x") %in% names(formals(ytrafo))))
  #    mf$yx <- "matrix"
  #}
  
  mf$nmax <- control$nmax
  ## evaluate model.frame
  mf[[1L]] <- quote(extree_data)
  
  d <- eval(mf, parent.frame())
  #subset <- partykit:::.start_subset(d)  ## FIX ME: export function?
  # for now: function .start_subset copied directly into this code
  subset <- 1:NROW(model.frame(d))
  if (length(d$yxmissings) > 0) subset <- subset[!(subset %in% d$yxmissings)]
  
  weights <- model.weights(model.frame(d))
  
  
  #############################################
  # INSERTED HERE:
  ## prepare family:
  # check format of the family input and if necessary transform it to the required familiy list
  # family input can be of one of the following formats:
  # - gamlss.family object
  # - gamlss.family function
  # - character string with the name of a gamlss.family object
  # - function generating a list with the required information about the distribution
  # - character string with the name of a function generating a list with the required information about the distribution
  # - list with the required information about the distribution
  # - character string with the name of a distribution for which a list generating function is provided in disttree
  {
    if(is.character(family)) {
      getfamily <- try(getAnywhere(paste("dist", family, sep = "_")), silent = TRUE)
      if(length(getfamily$objs) == 0L) getfamily <- try(getAnywhere(family), silent = TRUE)
      if(length(getfamily$objs) == 0L) {
        stop("unknown 'family' specification")
      } else {
        gamlssobj <- ("gamlss.dist" %in% unlist(strsplit(getfamily$where[1], split = ":")))
        family <- getfamily[[2]][[1]]() #first found is chosen 
        family$gamlssobj <- gamlssobj
      }
      #if(!(inherits(family, "try-error")))family <- family[[2]]$`package:disttree`()    
      # FIX ME: better selection of dist function
    }
    
    # if family is a gamlss family object or gamlss family function
    if(is.function(family)) family <- family()
    if(inherits(family, "gamlss.family")) family <- make_dist_list(family, bd = bd)
    
    if(!is.list(family)) stop ("unknown family specification")
    if(!(all(c("ddist", "sdist", "link", "linkfun", "linkinv", "mle", "startfun") %in% names(family)))) stop("family needs to specify a list with ...")
    # linkinvdr only used in the method vcov for type = "parameter"
  }
  
  np <- length(family$link)
  #############################################
  
  
  #############################################
  # INSERT HERE: ytrafo
  ## wrapper function to apply distfit in ctree
  ytrafo <- function(data, weights = NULL, control) {
    
    Y <- model.frame(data, yxonly = TRUE)
    if(dim(Y)[2] > 1) stop("response variable has to be univariate") 
    Y <- Y[,1]
    
    modelscores_decor <- function(subset, weights, estfun = TRUE, object = TRUE, info = NULL) {
      
      ys <- Y[subset]
      subweights <- if(is.null(weights) || (length(weights)==0L)) weights else weights[subset] ## FIX ME: scores with or without weights?
      # start <- if(!(is.null(info$coefficients))) info$coefficients else NULL
      start <- info$coefficients
      
      model <- distfit(ys, family = family, weights = subweights, start = start,
                       vcov = (decorrelate == "vcov"), type.hessian = "analytic", 
                       estfun = estfun, censtype = censtype, censpoint = censpoint, ocontrol = ocontrol, ...)
      
      if(estfun) {
        ef <- as.matrix(model$estfun)
        
        if(decorrelate != "none") {
          n <- NROW(ef)
          ef <- ef/sqrt(n)
          
          vcov <- if(decorrelate == "vcov") {
            vcov(model, type = "link") * n
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
        
        estfun <- matrix(0, ncol = ncol(ef), nrow = nrow(data$data)) 
        estfun[subset,] <- ef
        ### now automatically if(!(is.null(weights) || (length(weights)==0L))) estfun <- estfun / weights # estfun has to be unweighted for ctree
      } else estfun <- NULL
      
      
      
      object <-  if(object) model else NULL
      
      ret <- list(estfun = estfun,
                  coefficients = coef(model, type = "parameter"),
                  objfun = logLik(model),  # optional function to be maximized (FIX: negative?/minimize?)
                  object = object,
                  converged = model$converged  # FIX ME: warnings if distfit does not converge
      )
      return(ret)
    }
    
    return(modelscores_decor)
  }    
  
  #############################################
  
  
  
  if (is.null(control$update))
    control$update <- TRUE
  nf <- names(formals(ytrafo))
  if (all(c("data", "weights", "control") %in% nf))
    ytrafo <- ytrafo(data = d, weights = weights, control = control)
  nf <- names(formals(ytrafo))
  stopifnot(all(c("subset", "weights", "info", "estfun", "object") %in% nf) ||
              all(c("y", "x", "weights", "offset", "start") %in% nf))
  
  if (is.function(converged)) {
    stopifnot(all(c("data", "weights", "control") %in% names(formals(converged))))
    converged <- converged(d, weights, control = control)
  } else {
    converged <- TRUE
  }            
  
  update <- function(subset, weights, control, doFit = TRUE)
    extree_fit(data = d, trafo = ytrafo, converged = converged, partyvars = d$variables$z, 
                          subset = subset, weights = weights, ctrl = control, doFit = doFit)
  if (!doFit) return(list(d = d, update = update))
  tree <- update(subset = subset, weights = weights, control = control)
  trafo <- tree$trafo
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
  
  ################################################
  ## INSERTED HERE:
  # distributional fit: calculate coefficients for terminal nodes using distfit()
  
  # number of terminal nodes
  n_tn <- width(ret)
  # predicted terminal nodes for the given data
  pred_tn <- predict(ret, type = "node")
  # ids of terminal nodes
  id_tn <- as.vector(unique(pred_tn))
  
  if(is.null(weights) || (length(weights)==0L)) weights <- numeric(nrow(data)) + 1
  
  ## get coefficients for terminal nodes:
  Y <- ret$fitted$`(response)`
  # first iteration out of loop:
  model1 <- distfit(y = Y[(id_tn[1]==pred_tn)], family = family, weights = weights[(id_tn[1]==pred_tn)], start = NULL,
                    vcov = FALSE, type.hessian = "analytic", 
                    estfun = FALSE, censtype = censtype, censpoint = censpoint, ocontrol = ocontrol, ...)
  coefficients_par <- matrix(nrow = n_tn, ncol = length(model1$par))
  # coefficients_eta <- matrix(nrow = n_tn, ncol = length(model1$eta)) 
  colnames(coefficients_par) <- names(model1$par)
  # colnames(coefficients_eta) <- names(model1$eta)
  rownames(coefficients_par) <- as.character(id_tn)
  # rownames(coefficients_eta) <- as.character(id_tn)
  
  coefficients_par[1,] <- model1$par
  # coefficients_eta[1,] <- model1$eta
  
  loglik <- sum(model1$ddist(Y[(id_tn[1]==pred_tn)], log = TRUE))
  
  if(n_tn>1){
    for(i in (2:n_tn)){
      model <- distfit(y = Y[(id_tn[i]==pred_tn)], family = family, weights = weights[(id_tn[i]==pred_tn)], start = NULL,
                       vcov = FALSE, type.hessian = "analytic", 
                       estfun = FALSE, censtype = censtype, censpoint = censpoint, ocontrol = ocontrol, ...)
      coefficients_par[i,] <- model$par
      # coefficients_eta[i,] <- model$eta
      loglik <- loglik + sum(model$ddist(Y[(id_tn[i]==pred_tn)], log = TRUE))
    }
  }
  
  ret$coefficients <- coefficients_par
  ret$fitted$`(fitted.response)` <- predict(ret, type = "response")
  ret$loglik <- loglik

  ## extend class and keep original call/family/control
  ret$info$call <- cl
  ret$info$family <- family   
  ret$info$ocontrol <- ocontrol
  ret$info$formula <- formula
  ret$info$censpoint <- censpoint
  ret$info$censtype <- censtype
  
  groupcoef <- ret$coefficients
  if(!(is.null(groupcoef))){
    if(is.vector(groupcoef)) {
      groupcoef <- t(as.matrix(groupcoef))
      rownames(groupcoef) <- 1
    }
    ret$fitted.par <- groupcoef[paste(ret$fitted[,1]),]
    rownames(ret$fitted.par) <- c(1: (length(ret$fitted.par[,1])))
    ret$fitted.par <- as.data.frame(ret$fitted.par)
  }
  
###################################################
  
  
  class(ret) <- c("disttree", class(ret))
  
  ### doesn't work for Surv objects
  # ret$terms <- terms(formula, data = mf)
  ret$terms <- d$terms$all
  ### need to adjust print and plot methods
  ### for multivariate responses
  ### if (length(response) > 1) class(ret) <- "party"
  return(ret)
}





## FIX ME: distextree_control()


distextree_control <- function(type.tree = NULL,
                               criterion, 
                               logmincriterion, 
                               minsplit = 20L,
                               minbucket = 7L, 
                               minprob = 0.01, 
                               nmax = Inf,
                               stump = FALSE,
                               lookahead = FALSE, ### try trafo() for daugther nodes before implementing the split
                               MIA = FALSE,
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
                               testflavour = c("ctree", "exhaustive", "mfluc"),
                               bonferroni = FALSE,
                               splitflavour = c("ctree", "exhaustive"),
                               update = NULL) 
{
  
  if(!is.null(type.tree)) {
    
    if(!type.tree %in% c("ctree", "mob")) {
      stop("type.tree can only be set to 'ctree' or 'mob'")
    } else {
      
      if(type.tree == "ctree"){
        if(!("ctree" %in% splitflavour & "ctree" %in% testflavour)){
          stop("for type.tree = 'ctree' testflavour and splitflavour can not be set to other than 'ctree'")
        } else {
          
          testflavour <- "ctree"
          splitflavour <- "ctree"
          
          if(FALSE){
            control <- ctree_control(teststat = c("quadratic", "maximum"), 
                                     splitstat = c("quadratic", "maximum"), ### much better for q > 1, max was default
                                     splittest = FALSE,
                                     testtype = c("Bonferroni", "MonteCarlo", 
                                                  "Univariate", "Teststatistic"),
                                     pargs = GenzBretz(),
                                     nmax = c("yx" = Inf, "z" = Inf),
                                     alpha = 0.05, 
                                     mincriterion = 1 - alpha, 
                                     logmincriterion = log(mincriterion), 
                                     minsplit = minsplit, 
                                     minbucket = minbucket, 
                                     minprob = minprob, 
                                     stump = FALSE, 
                                     lookahead = FALSE,	### try trafo() for daugther nodes before implementing the split
                                     MIA = FALSE,	### DOI: 10.1016/j.patrec.2008.01.010
                                     nresample = 9999L, 
                                     tol = sqrt(.Machine$double.eps),
                                     maxsurrogate = 0L, 
                                     numsurrogate = FALSE,
                                     mtry = Inf, 
                                     maxdepth = Inf, 
                                     multiway = FALSE, 
                                     splittry = 2L, 
                                     intersplit = FALSE,
                                     majority = FALSE, 
                                     caseweights = TRUE, 
                                     applyfun = NULL, 
                                     cores = NULL,
                                     saveinfo = TRUE,
                                     update = NULL)
          }
          
        }
      }
      
      if(type.tree == "mob"){        
        if(!("exhaustive" %in% splitflavour & "mfluc" %in% testflavour)){
          stop("for type.tree = 'ctree' testflavour can not be set to other than 'mfluc' and splitflavour can not be set to other than 'exhaustive'")
        } else {
          
          testflavour <- "mfluc"
          splitflavour <- "exhaustive"
          
          if(FALSE){
            control <- mob_control(alpha = 0.05, 
                                   bonferroni = TRUE, 
                                   minsize = minbucket,   ## FIX ME: minsize minimal number of obs in node after splitting? 
                                   maxdepth = Inf,
                                   mtry = Inf, 
                                   trim = 0.1, 
                                   breakties = FALSE, 
                                   parm = NULL, 
                                   dfsplit = TRUE, 
                                   prune = NULL, 
                                   restart = TRUE,
                                   verbose = FALSE, 
                                   caseweights = TRUE, 
                                   ytype = "vector", 
                                   xtype = "matrix",
                                   terminal = "object", 
                                   inner = terminal, 
                                   model = TRUE,
                                   numsplit = "left", 
                                   catsplit = "binary", 
                                   vcov = "opg", 
                                   ordinal = "chisq", 
                                   nrep = 10000,
                                   minsplit = minsize,   ## FIX ME: per default minsplit and minbucket are set to the same value? both minsize?
                                   minbucket = minsize,
                                   applyfun = NULL, 
                                   cores = NULL)
          }
          
        }
      }
    }
  }
  
  
  
  
  
  
  
  control <- extree_control(criterion = criterion, 
                            logmincriterion = logmincriterion, 
                            minsplit = minsplit,
                            minbucket = minbucket, 
                            minprob = minprob, 
                            nmax = nmax,
                            stump = stump,
                            lookahead = lookahead, ### try trafo() for daugther nodes before implementing the split
                            MIA = MIA,
                            maxsurrogate = maxsurrogate, 
                            numsurrogate = numsurrogate,
                            mtry = mtry,
                            maxdepth = maxdepth, 
                            multiway = multiway, 
                            splittry = splittry,
                            majority = majority, 
                            caseweights = caseweights, 
                            applyfun = applyfun, 
                            cores = cores,
                            saveinfo = saveinfo,
                            testflavour = testflavour,
                            bonferroni = bonferroni,
                            splitflavour = splitflavour,
                            update = update)
}






if(FALSE){
  
ctree_control <- function
  (
    
    teststat = c("quadratic", "maximum"), 
    splitstat = c("quadratic", "maximum"), ### much better for q > 1, max was default
    splittest = FALSE,
    testtype = c("Bonferroni", "MonteCarlo", 
                 "Univariate", "Teststatistic"),
    pargs = GenzBretz(),
    nmax = c("yx" = Inf, "z" = Inf),
    alpha = 0.05, 
    mincriterion = 1 - alpha, 
    logmincriterion = log(mincriterion), 
    minsplit = 20L, 
    minbucket = 7L, 
    minprob = 0.01, 
    stump = FALSE, 
    lookahead = FALSE,	### try trafo() for daugther nodes before implementing the split
    MIA = FALSE,	### DOI: 10.1016/j.patrec.2008.01.010
    nresample = 9999L, 
    tol = sqrt(.Machine$double.eps),
    maxsurrogate = 0L, 
    numsurrogate = FALSE,
    mtry = Inf, 
    maxdepth = Inf, 
    multiway = FALSE, 
    splittry = 2L, 
    intersplit = FALSE,
    majority = FALSE, 
    caseweights = TRUE, 
    applyfun = NULL, 
    cores = NULL,
    saveinfo = TRUE,
    update = NULL
  ) {
    
    testtype <- match.arg(testtype, several.ok = TRUE)
    if (length(testtype) == 4) testtype <- testtype[1]
    ttesttype <- testtype
    if (length(testtype) > 1) {
      stopifnot(all(testtype %in% c("Bonferroni", "MonteCarlo")))
      ttesttype <- "MonteCarlo"
    }
    
    splitstat <- match.arg(splitstat)
    teststat <- match.arg(teststat)
    
    if (!caseweights)
      stop("only caseweights currently implemented in ctree")
    
    c(extree_control(criterion = ifelse("Teststatistic" %in% testtype, 
                                        "statistic", "p.value"),
                     logmincriterion = logmincriterion, minsplit = minsplit, 
                     minbucket = minbucket, minprob = minprob, 
                     nmax = nmax, stump = stump, lookahead = lookahead,
                     mtry = mtry, maxdepth = maxdepth, multiway = multiway, 
                     splittry = splittry, MIA = MIA, maxsurrogate = maxsurrogate, 
                     numsurrogate = numsurrogate,
                     majority = majority, caseweights = caseweights, 
                     applyfun = applyfun, saveinfo = saveinfo,  ### always
                     testflavour = "ctree", 
                     bonferroni = "Bonferroni" %in% testtype, 
                     splitflavour = "ctree", update = update),
      list(teststat = teststat, splitstat = splitstat, splittest = splittest, pargs = pargs,
           testtype = ttesttype, nresample = nresample, tol = tol,
           intersplit = intersplit))
  }
  
  


  
extree_control <- function
  (
    criterion, 
    logmincriterion, 
    minsplit = 20L,
    minbucket = 7L, 
    minprob = 0.01, 
    nmax = Inf,
    stump = FALSE,
    lookahead = FALSE, ### try trafo() for daugther nodes before implementing the split
    MIA = FALSE,
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
    testflavour = c("ctree", "exhaustive", "mfluc"),
    bonferroni = FALSE,
    splitflavour = c("ctree", "exhaustive"),
    update = NULL
  ) {
    
    ## apply infrastructure for determining split points
    if (is.null(applyfun)) {
      applyfun <- if(is.null(cores)) {
        lapply
      } else {
        function(X, FUN, ...)
          parallel::mclapply(X, FUN, ..., mc.cores = cores)
      }
    }
    
    ### well, it is implemented but not correctly so
    if (multiway & maxsurrogate > 0L)
      stop("surrogate splits currently not implemented for multiway splits")
    
    if (MIA && maxsurrogate > 0)
      warning("Mixing MIA splits with surrogate splits does not make sense")
    
    if (MIA && majority)
      warning("Mixing MIA splits with majority does not make sense")
    
    list(criterion = criterion, logmincriterion = logmincriterion,
         minsplit = minsplit, minbucket = minbucket, 
         minprob = minprob, stump = stump, nmax = nmax,
         lookahead = lookahead, mtry = mtry,
         maxdepth = maxdepth, multiway = multiway, splittry = splittry,
         MIA = MIA, maxsurrogate = maxsurrogate, 
         numsurrogate = numsurrogate, majority = majority,
         caseweights = caseweights, applyfun = applyfun,
         saveinfo = saveinfo, testflavour = match.arg(testflavour), 
         bonferroni = bonferroni,
         splitflavour = match.arg(splitflavour), update = update)
  }
}