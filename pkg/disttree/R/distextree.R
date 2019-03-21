########################################################################
# new version of disttree using extree directly without applying ctree #
########################################################################

distextree <- function(formula, 
                       data, 
                       subset, 
                       weights, 
                       family = NO(), 
                       na.action = na.pass, 
                       offset, cluster,
                       control = distextree_control(...), 
                       type.hessian = c("checklist", "analytic", "numeric"),
                       converged = NULL, 
                       scores = NULL, 
                       doFit = TRUE, 
                       bd = NULL, # terminal_objects = FALSE,
                       decorrelate = "none", 
                       censtype = "none", 
                       censpoint = NULL,
                       ocontrol = list(), 
                       ...) {
  
  decorrelate <- match.arg(decorrelate)

  ## Keep call
  cl <- match.call(expand.dots = TRUE)
  if(missing(data)) data <- environment(formula)
  
  # Check formula
  oformula <- as.formula(formula)
  formula <- Formula::as.Formula(formula)
  if(length(formula)[2L] > 1L) {
    formula <- Formula::Formula(formula(formula, rhs = 2L))  
    ## FIX ME: if rhs has more than 1 element it is here assumed that partitioning variables are handed over on 2nd slot
    warning("formula must not have more than one RHS parts (only partitioning variables allowed)")
  }
  
  ## Set up model.frame() call
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", 
               "offset", "cluster", "scores"), names(mf), 0L)

  mf <- mf[c(1L, m)]
  mf$yx <- "matrix"      # FIX ME: in ctree "none" ?
  mf$nmax <- control$nmax
  mf$ytype <- "vector"
  
  ## Evaluate model.frame
  mf[[1L]] <- quote(partykit::extree_data)
  
  d <- eval(mf, parent.frame())

  subset <- partykit:::.start_subset(d)
  weights <- model.weights(model.frame(d))
  
  if(is.null(control$partyvars)) control$partyvars <- d$variables$z
  
  # Set up family 
  family <- distfamily(family)
  np <- length(family$link)
  
  #Y <- d$yx[[1]]
  # check whether d$yx really contains only the response (and no covariates)
  if(length(d$yx) > 1) stop("covariates can only be used as split variables")
  
  if(NCOL(d$yx[[1]]) > 1) stop("response variable has to be univariate") 
  if(inherits(d$yx[[1]], "interval")) stop("can not deal with binned intervals yet") 
  
  ## Set up wrapper function for distexfit
  ytrafo <- function(subset, weights, estfun = FALSE, object = FALSE, info = NULL) {
    
    ys <- d$yx[[1]][subset]  # necessary to get response data into the function
    subweights <- if(is.null(weights) || (length(weights)==0L)) weights else weights[subset]
    ## FIX ME: scores with or without weights?
    # start <- if(!(is.null(info$coefficients))) info$coefficients else NULL
    start <- info$coefficients
    
    model <- disttree::distexfit(ys, family = family, weights = subweights, start = start,
                                 vcov = (decorrelate == "vcov"), type.hessian = type.hessian, 
                                 estfun = estfun, censtype = censtype, censpoint = censpoint, 
                                 ocontrol = ocontrol)
    
    if(estfun) {
      ef <- as.matrix(model$estfun) # distexfit returns weighted scores!
      
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
      
      estfun <- matrix(0, ncol = ncol(ef), nrow = nrow(d$data)) 
      estfun[subset,] <- ef
      ### now automatically if(!(is.null(weights) || (length(weights)==0L))) estfun <- estfun / weights 
    } else {
      estfun <- NULL
    }
    
    rval <- list(estfun = estfun,
                unweighted = FALSE, # unweighted = TRUE would prevent estfun / w in extree_fit
                coefficients = coef(model, type = "parameter"),
                objfun = -logLik(model),  # optional function to be minimized 
                object = if(object) model else NULL,
                converged = model$converged  # FIX ME: warnings if distexfit does not converge
    )
    return(rval)
  }

  # adaption to new version of extree (different structure of ytrafo, compare to ctree)
  if (is.null(control$update))
    control$update <- TRUE
  nf <- names(formals(ytrafo))
  #if (all(c("data", "weights", "control") %in% nf))     
  #  ytrafo <- ytrafo(data = d, weights = weights, control = control)
  #nf <- names(formals(ytrafo))
  stopifnot(all(c("subset", "weights", "info", "estfun", "object") %in% nf) ||
              all(c("y", "x", "weights", "offset", "start") %in% nf))
  #############################################
   
  ## FIX ME: implement function checking whether all response values are equal in one node
  # if so, return FALSE
  # still to do: which control arguments should be used here? (now not used)
  # allow to hand over additional conditions in further 'converged' functions?
  converged <- function(data, weights, control){
    #if(is.null(weights)) weights <- rep.int(1, NROW(data$yx[[1]]))
    convfun <- function(subset, weights){
      ys <- data$yx[[1]][subset]
      ws <- if(is.null(weights)) rep.int(1, NROW(ys)) else weights[subset]
      conv <- length(unique(ys[ws > 0]))>1
      return(conv)
    }
    return(convfun)
  }
  
  # not necessary since converged is fixed here (defined above), 
  # but has to be checked if further 'converged' functions are handed over
  # FIX ME: various 'converged' functions?
  # FIX ME: how is this function applied on single nodes? 
  if (is.function(converged)) {
    stopifnot(all(c("data", "weights", "control") %in% names(formals(converged))))
    converged <- converged(d, weights, control = control)
  } else {
    converged <- TRUE
  }            

  ## Set up wrapper for extree_fit with predefined fit function
  update <- function(subset, weights, control, doFit = TRUE) {
    partykit::extree_fit(data = d, trafo = ytrafo, converged = converged, partyvars = control$partyvars, 
                         subset = subset, weights = weights, ctrl = control, doFit = doFit)
  }

  if (!doFit) return(list(d = d, update = update))

  tree <- update(subset = subset, weights = weights, control = control)
  trafo <- tree$trafo
  
  mf <- model.frame(d)
  if (is.null(weights)) weights <- rep(1, nrow(mf))
  
  fitted <- data.frame("(fitted)" = fitted_node(tree$nodes, mf), 
                       "(weights)" = weights,
                       check.names = FALSE)

  fitted[[3]] <- mf[, d$variables$y, drop = TRUE]
  names(fitted)[3] <- "(response)"

  rval <- partykit::party(tree$nodes, 
                data = if(control$model) mf else mf[0,],
                fitted = fitted,
                terms = d$terms$all,
                info = list(
                  call = cl,
                  formula = formula,
                  family = family,
                  #terms = list(response = d$terms$yx, partitioning = d$terms$z), #FIXME: Braucht man?
                  fit = distexfit,
                  control = control#,
                  #dots = list(...) #FIXME: Braucht man?!
              )
  )

  ## COMMENT LISA's VERSION (begin)
  #rval <- partykit::party(tree$nodes, data = mf, fitted = fitted, 
  #             info = list(call = match.call(), control = control))

  #rval$update <- update #FIXME: Auch hier, braucht man das?!
  #rval$trafo <- trafo #FIXME: Will man das noch haben?
  #class(rval) <- c("constparty", class(rval))
  #
  ### INSERTED HERE:
  ## distributional fit: calculate coefficients for terminal nodes using distexfit()
  #
  ### FIX ME: first check whether there is already a fitted model in each of the nodes, 
  ## if so, extract coefficients instead of calculating them again in the following lines
  #
  ## number of terminal nodes
  #n_tn <- width(rval)
  ## predicted terminal nodes for the given data
  #pred_tn <- predict(rval, type = "node")
  ## ids of terminal nodes
  #id_tn <- as.vector(unique(pred_tn))
  #
  #if(is.null(weights) || (length(weights)==0L)) weights <- numeric(nrow(data)) + 1
  #
  ### get coefficients for terminal nodes:
  #Y <- rval$fitted$`(response)`
  ## first iteration out of loop:
  #model1 <- disttree::distexfit(y = Y[(id_tn[1]==pred_tn)], family = family, weights = weights[(id_tn[1]==pred_tn)], start = NULL,
  #                            vcov = FALSE, type.hessian = type.hessian, 
  #                            estfun = FALSE, censtype = censtype, censpoint = censpoint, ocontrol = ocontrol)
  #coefficients_par <- matrix(nrow = n_tn, ncol = length(model1$par))
  ## coefficients_eta <- matrix(nrow = n_tn, ncol = length(model1$eta)) 
  #colnames(coefficients_par) <- names(model1$par)
  ## colnames(coefficients_eta) <- names(model1$eta)
  #rownames(coefficients_par) <- as.character(id_tn)
  ## rownames(coefficients_eta) <- as.character(id_tn)
  #
  #coefficients_par[1,] <- model1$par
  ## coefficients_eta[1,] <- model1$eta
  #
  #loglik <- sum(model1$ddist(Y[(id_tn[1]==pred_tn)], log = TRUE))
  #
  #if(n_tn>1){
  #  for(i in (2:n_tn)){
  #    model <- disttree::distexfit(y = Y[(id_tn[i]==pred_tn)], family = family, weights = weights[(id_tn[i]==pred_tn)], start = NULL,
  #                               vcov = FALSE, type.hessian = type.hessian, 
  #                               estfun = FALSE, censtype = censtype, censpoint = censpoint, ocontrol = ocontrol)
  #    coefficients_par[i,] <- model$par
  #    # coefficients_eta[i,] <- model$eta
  #    loglik <- loglik + sum(model$ddist(Y[(id_tn[i]==pred_tn)], log = TRUE))
  #  }
  #}
  #
  #rval$coefficients <- coefficients_par
  ## rval$fitted$`(fitted.response)` <- predict(rval, type = "response")
  #rval$loglik <- loglik

  ### extend class and keep original call/family/control
  #rval$info$call <- cl
  #rval$info$family <- family  
  #rval$info$ocontrol <- ocontrol #TODO: Muss man anders loesen
  #rval$info$formula <- formula
  #rval$info$censpoint <- censpoint #TODO: Muss man anders loesen
  #rval$info$censtype <- censtype #TODO: Muss man anders loesen
  #
  #groupcoef <- rval$coefficients #TODO: Was sind groupcoef?!
  #if(!(is.null(groupcoef))){
  #  if(is.vector(groupcoef)) {
  #    groupcoef <- t(as.matrix(groupcoef))
  #    rownames(groupcoef) <- 1
  #  }
  #  rval$fitted.par <- groupcoef[paste(rval$fitted[,1]),]
  #  rownames(rval$fitted.par) <- c(1: (length(rval$fitted.par[,1])))
  #  rval$fitted.par <- as.data.frame(rval$fitted.par)
  #}
  # 
  # class(rval) <- c("disttree", class(rval))
  # 
  # ### doesn't work for Surv objects #TODO: Das heisst?!
  # # rval$terms <- terms(formula, data = mf)
  # rval$terms <- d$terms$all  
  # ### need to adjust print and plot methods
  # ### for multivariate responses   
  # ### if (length(response) > 1) class(rval) <- "party"
  # return(rval)

  ## COMMENT LISA's VERSION (end)

  class(rval) <- c("modelparty", class(rval))  # TODO: Either model or constparty object!

  ### Add modelinfo (object) and estfun if not there yet, but wanted
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

    if (i %in% which_terminals) winfo <- control$terminal else
      winfo <- control$inner

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




distextree_control <- function(type.tree = NULL,
                               criterion = c("p.value", "statistic"),
                               logmincriterion = log(0.95), 
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
                               bonferroni = TRUE,       # FIX ME: different than extree_control default
                               update = NULL,
                               splitflavour = c("ctree", "exhaustive"),  
                               testflavour = c("ctree", "mfluc", "guide"),   # FIX ME: "exhaustive"
                               guide_interaction = FALSE,
                               #guide_unweighted = FALSE,
                               nresample = 9999L,
                               intersplit = FALSE,
                               tol = sqrt(.Machine$double.eps),
                               partyvars = NULL,   # per default it will be set to partyvars = data$variables$z after applying extree_data
                               
                               ## FIX ME: additional arguments for ctree
                               teststat = c("quadratic", "maximum"), 
                               splitstat = c("quadratic", "maximum"), ### much better for q > 1, max was default
                               splittest = FALSE,
                               testtype = c("Bonferroni", "MonteCarlo", 
                                            "Univariate", "Teststatistic"),
                               pargs = GenzBretz(),
                               ## FIX ME: additional arguments for exhaustive
                               restart = TRUE,
                               breakties = FALSE,
                               parm = NULL,  ## FIX ME: match with partyvars above
                               dfsplit = TRUE,
                               vcov = c("opg", "info", "sandwich"),
                               ordinal = c("chisq", "max", "L2"),
                               ytype = c("vector", "data.frame", "matrix"),
                               #nrep = 10000L,
                               terminal = "object",
                               model = TRUE,
                               inner = "object",
                               trim = 0.1,
                               ## FIX ME: additional arguments for guide
                               guide_parm = NULL,  # a vector of indices of the parameters (incl. intercept) for which estfun should be considered
                               guide_testtype = c("max", "sum", "coin"),
                               interaction = FALSE,
                               guide_decorrelate = "vcov",   # needs to be set to other than "none" for testtype max and sum 
                               # unless ytrafo returns decorrelated scores
                               # FIX ME: c("none","vcov","opg")
                               xgroups = NULL,  # number of categories for split variables (optionally breaks can be handed over)
                               ygroups = NULL,  # number of categories for scores (optionally breaks can be handed over)
                               weighted.scores = FALSE   # logical, should scores be weighted in GUIDE 
                               ) 
{
  
  add_control <- NULL
  
  if (length(criterion) == 2) criterion <- criterion[1]
  if (length(testflavour) == 3 & is.null(type.tree)) testflavour <- testflavour[1]
  if (length(splitflavour) == 2 & is.null(type.tree)) splitflavour <- splitflavour[1]
  
  
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

        }
      }
      
      if(type.tree == "mob"){        
        if(!("exhaustive" %in% splitflavour & "mfluc" %in% testflavour)){
          stop("for type.tree = 'mob' testflavour can not be set to other than 'mfluc' and splitflavour can not be set to other than 'exhaustive'")
        } else {
          
          testflavour <- "mfluc"
          splitflavour <- "exhaustive"  
          
        }
      }
    }
  }
  
  
  
  
  
  if(testflavour == "ctree" | splitflavour == "ctree") {
    
    testtype <- match.arg(testtype, several.ok = TRUE)
    if (length(testtype) == 4) testtype <- testtype[1]
    ttesttype <- testtype
    if (length(testtype) > 1) {
      stopifnot(all(testtype %in% c("Bonferroni", "MonteCarlo")))
      ttesttype <- "MonteCarlo"
    }
    
    if (MIA && maxsurrogate > 0)
      warning("Mixing MIA splits with surrogate splits does not make sense")
    
    if (MIA && majority)
      warning("Mixing MIA splits with majority does not make sense")
    
    splitstat <- match.arg(splitstat)
    teststat <- match.arg(teststat)
    
    if (!caseweights)
      stop("only caseweights currently implemented in ctree")
    
    add_control <- list(teststat = teststat, 
                        splitstat = splitstat, 
                        splittest = splittest, 
                        pargs = pargs,
                        testtype = testtype, 
                        MIA = MIA,
                        nresample = nresample,     # FIX ME
                        tol = tol,                 # FIX ME
                        intersplit = intersplit,                  # FIX ME
                        terminal = terminal,
                        model = model,
                        inner = inner)
    
    criterion2 = ifelse("Teststatistic" %in% testtype, "statistic", "p.value")
    if(criterion != criterion2) stop("'criterion' does not match 'testtype'")
    
    bonferroni2 = "Bonferroni" %in% testtype
    if(bonferroni != bonferroni2) stop("'bonferroni' does not match 'testtype'")
    
  }
  
  if(testflavour == "mfluc" | splitflavour == "exhaustive") {
    add_control <- c(add_control,
                     list(restart = restart,
                          breakties = breakties,
                          parm = parm,
                          dfsplit = dfsplit,
                          vcov = vcov,
                          ordinal = ordinal,
                          ytype = ytype,
                          #nrep = nrep,
                          terminal = terminal,
                          model = model,
                          inner = inner,
                          trim = trim,
                          intersplit = intersplit))                  # FIX ME 
  }
  
  if(testflavour == "guide") {
    
    if(!criterion == "p.value" & guide_interaction) stop("For testflavour GUIDE with interaction tests only 'p.value' can be selected as criterion")
    add_control <- c(add_control,
                     list(guide_parm = guide_parm,  # a vector of indices of parameters for which estfun should be considered
                          guide_testtype = guide_testtype,
                          interaction = interaction,
                          guide_decorrelate = guide_decorrelate,   # needs to be set to other than "none" for testtype max and sum 
                          # unless ytrafo returns decorrelated scores
                          # FIX ME: c("none","vcov","opg")
                          xgroups = xgroups,  # number of categories for split variables (optionally breaks can be handed over)
                          ygroups = ygroups,  # number of categories for scores (optionally breaks can be handed over)
                          weighted.scores = weighted.scores)   # logical, should scores be weighted
                     )
  }
  
  
  if (splitflavour == "ctree") splitfun = partykit:::.ctree_split() 
  if (splitflavour == "exhaustive") splitfun = partykit:::.objfun_split()   
  
  if (testflavour == "ctree") selectfun = partykit:::.ctree_select()
  if (testflavour == "mfluc") selectfun = partykit:::.mfluc_select()
  if (testflavour == "guide") selectfun = partykit:::.guide_select()
  
  svselectfun = partykit:::.ctree_select()
  svsplitfun = partykit:::.ctree_split(minbucket = 0)
  
  ## FIX ME: argumets numsplit in mob and intersplit in ctree/extree
  # intersplit <- numsplit == "center"
  
  control <- c(partykit:::extree_control(criterion = criterion, 
                                         logmincriterion = logmincriterion, 
                                         minsplit = minsplit,
                                         minbucket = minbucket, 
                                         minprob = minprob, 
                                         nmax = nmax,
                                         stump = stump,
                                         lookahead = lookahead, ### try trafo() for daugther nodes before implementing the split
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
                                         bonferroni = bonferroni,  
                                         update = update,
                                         selectfun = selectfun,
                                         splitfun = splitfun,
                                         svselectfun = svselectfun,
                                         svsplitfun = svsplitfun),   
               add_control)
  
  return(control)
}



if(FALSE) {
# use different version of .select than partykit:::select for GUIDE as selectfun 
# (returns p.values from curvature and interaction tests instead of p.value and teststatistic)
.select_g <- function(model, trafo, data, subset, weights, whichvar, ctrl, FUN) {
  ret <- list(criteria = matrix(NA, nrow = 2L, ncol = ncol(model.frame(data))))
  rownames(ret$criteria) <- c("statistic", "p.value")
  colnames(ret$criteria) <- names(model.frame(data))
  if (length(whichvar) == 0) return(ret)
  ### <FIXME> allow joint MC in the absense of missings; fix seeds
  ### write ctree_test / ... with whichvar and loop over variables there
  ### </FIXME>
  for (j in whichvar) {
    tst <- FUN(model = model, trafo = trafo, data = data, 
               subset = subset, weights = weights, j = j, 
               SPLITONLY = FALSE, ctrl = ctrl)
    ret$criteria["statistic",j] <- - as.numeric(tst["p.curv"])
    ret$criteria["p.value",j] <- - as.numeric(tst["p.min"])
    
    # for testflavour = "guide" only "p.value" can be chosen as criterion
    # because 2 p.values need to be returned (from curvature test and the min from curvature and interaction test)
    # the min is stored as 'p.value' in the returned matrix
    # if this min is from an interaction test, two covariates will have the same p.value
    # in this case the one with the lower p.value from the curvature test should be chosen
    # in extree: among the two covariates with the lowest p.value (highest negativ p.value) the one with the higher test statistic is chosen (ranked)
    # -> the negative p.value from the curvature test is stored as "statistic"
  }
  ret
}


.guide_select <- function(...)
  function(model, trafo, data, subset, weights, whichvar, ctrl) {
    args <- list(...)
    ctrl[names(args)] <- args
    .select_g(model, trafo, data, subset, weights, whichvar, ctrl, FUN = .guide_test)   # optional: use partykit:::.select
  }

.guide_test <- function(model, trafo, data, subset, weights, j, SPLITONLY = FALSE, ctrl) {
  
  ## TO DO: include SPLITONLY, MIA, ... ?
  ix <- data$zindex[[j]] ### data[[j, type = "index"]]
  iy <- data$yxindex ### data[["yx", type = "index"]]
  Y <- model$estfun  ## model from distexfit, returns wheigthed scores
  if(ctrl$guide_unweighted) Y <- Y/weights  ## FIX ME: influence of weights only on categorization
  x <- data[[j]]
  if(!is.null(subset)) {
    Y <- if(is.vector(Y)) Y[subset] else Y[subset,]
    x <- x[subset]
  }
  
  # if all values of the selected covariate are equal return highest possible p.value
  if(length(unique(x))<2) return(c(p.min = 1, p.curv = 1))
  
  # only select those other covariates which are also in partyvars
  ix_others <- c(1:NCOL(data$data))[data$variables$z + ctrl$partyvars == 2]
  ix_others <- ix_others[!ix_others == j]
  
  # split Y into 2 parts based on whether residuals (here: scores) are positive or negative
  # separately for each parameter
  for(k in 1:model$object$npar){
    respos <- (Y[,k]>0)
    #respos <- factor((Y[,k]>0), levels = c(FALSE,TRUE), labels = c(0,1))
    Y <- cbind(Y, respos)
    colnames(Y)[(model$object$npar + k)] <- paste0("rp",k)
  }
  
  if(is.numeric(x)){
    x_cat <- rep.int(1, length(x))
    q1 <- quantile(x, 0.25)
    q2 <- quantile(x, 0.50)
    q3 <- quantile(x, 0.75)
    for(l in 1: length(x)){
      if(x[l] > q1 & x[l] <= q2) x_cat[l] <- 2
      if(x[l] > q2 & x[l] <= q3) x_cat[l] <- 3
      if(x[l] > q3) x_cat[l] <- 4
    }
    x_cat <- factor(x_cat, levels = c(1:4))
  } else {
    x_cat <- x
  }
  
  ## compute curvature test (for each parameter separately)
  p.curv <- chisq.test(x = x_cat, y = Y[,(model$object$npar+1)])$p.value
  if(model$object$npar > 1){
    for(k in 2:model$object$npar){
      p <- chisq.test(x = x_cat, y = Y[,(model$object$npar+k)])$p.value
      if(p < p.curv) p.curv <- p
    }
  }  
  
  p.min <- p.curv
  
  ## compute interaction test (for each parameter and for each of the other covariates separately)
  # only keep test if p.value is smaller than the one resulting from the curvature test
  if(ctrl$guide_interaction & length(ix_others)>0){
    
    for(v in ix_others){
      
      xo <- data[[v]]
      if(!is.null(subset)) xo <- xo[subset]
      
      # only consider other covariate for interaction test if not all of its values are equal
      if(length(unique(xo))>1){
        
        x_cat_2d <- rep.int(1, length(x))
        
        if(is.factor(x) & is.factor(xo)){
          c1 <- length(levels(x))
          c2 <- length(levels(xo))
          for(l in 1:length(x)){
            for(m in 1:c1){
              for(n in 1:c2){
                if(x[l] == levels(x)[m] & xo[l] == levels(xo)[n]) x_cat_2d[l] <- (m-1)*c2+n
              }
            }
          }
        }
        
        if(is.factor(x) & is.numeric(xo)){
          c1 <- length(levels(x))
          med_xo <- median(xo)
          for(l in 1:length(x)){
            for(m in 1:c1){
              if(x[l] == levels(x)[m]){
                x_cat_2d[l] <- if(xo[l] > med_xo) c1+m else m
              }
            }
          }
        }
        
        if(is.numeric(x) & is.factor(xo)){
          c2 <- length(levels(xo))
          med_x <- median(x)
          for(l in 1:length(x)){
            for(n in 1:c2){
              if(xo[l] == levels(xo)[n]){
                x_cat_2d[l] <- if(x[l] > med_x) c2+n else n
              }
            }
          }
        }
        
        if(is.numeric(x) & is.numeric(xo)){
          med_x <- median(x)
          med_xo <- median(xo)
          for(l in 1:length(x)){
            if(x[l] <= med_x) {
              x_cat_2d[l] <- if(xo[l] <= med_x) 1 else 2
            } else {
              x_cat_2d[l] <- if(xo[l] <= med_x) 3 else 4
            }
          }
        }
        
        p.int <- chisq.test(x = x_cat_2d, y = Y[,(model$object$npar+1)])$p.value
        if(model$object$npar > 1){
          for(k in 2:model$object$npar){
            p <- chisq.test(x = x_cat_2d, y = Y[,(model$object$npar+k)])$p.value
            if(p < p.int) p.int <- p
          }
        }  
        
        if(p.int < p.min) p.min <- p.int
      }
    }
  }
  
  ret <- c(p.min = p.min, p.curv = p.curv)
  
  return(ret)

}

}
