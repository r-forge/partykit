#### TO DO:
# - add control arguments for .guide_test:
#       guide_interaction = FALSE 
#       parm (indices of residuals/coefficients to be tested)
#       guide_unwheighted = FALSE 
# - if(guide_interaction) choose .select_g within .guide_select, else .select as usual


## arguments for guide_select -> guide_test: 
#guide_parm = NULL,  # a vector of indices of parameters for which estfun should be considered
#guide_testtype = c("max", "sum", "coin"),
#interaction = FALSE,
#guide_decorrelate = "vcov",   # needs to be set to other than "none" for testtype max and sum 
## unless ytrafo returns decorrelated scores
## FIX ME: c("none","vcov","opg")
#xgroups = NULL,  # number of categories for split variables (optionally breaks can be handed over)
#ygroups = NULL,  # number of categories for scores (optionally breaks can be handed over)
#weighted.scores = FALSE   # logical, should scores be weighted


.guide_select <- function(...) {
  function(model, trafo, data, subset, weights, whichvar, ctrl) {
    args <- list(...)
    ctrl[names(args)] <- args
    
    if(!"guide_parm" %in% names(ctrl)) ctrl$guide_parm <- NULL
    if(!"guide_testtype" %in% names(ctrl)) ctrl$guide_testtype <- c("max", "sum", "coin")
    if(!"interaction" %in% names(ctrl)) ctrl$interaction <- FALSE
    if(!"weighted.scores" %in% names(ctrl)) ctrl$weighted.scores <- FALSE
    if(!"xgroups" %in% names(ctrl)) ctrl$xgroups <- NULL
    if(!"ygroups" %in% names(ctrl)) ctrl$ygroups <- NULL
    if(!"guide_decorrelate" %in% names(ctrl)) ctrl$guide_decorrelate <- "vcov"
    
    if(ctrl$interaction){
      partykit:::.select_int(model, trafo, data, subset, weights, whichvar, ctrl, FUN = .guide_test_int)
    } else {
      partykit:::.select(model, trafo, data, subset, weights, whichvar, ctrl, FUN = .guide_test)
    }
  }
}


# already in extree.R
if(FALSE){
  .select <- function(model, trafo, data, subset, weights, whichvar, ctrl, FUN) {
    ret <- list(criteria = matrix(NA, nrow = 2L, ncol = ncol(model.frame(data))))
    rownames(ret$criteria) <- c("statistic", "p.value")
    colnames(ret$criteria) <- names(model.frame(data))
    if (length(whichvar) == 0) return(ret)
    for (j in whichvar) {
      tst <- FUN(model = model, trafo = trafo, data = data, 
                 subset = subset, weights = weights, j = j, 
                 SPLITONLY = FALSE, ctrl = ctrl)
      ret$criteria["statistic",j] <- tst$statistic
      ret$criteria["p.value",j] <- tst$p.value
    }
    ret
  }
}

.guide_test <- function(model, trafo, data, subset, weights, j, SPLITONLY = FALSE, ctrl) {
  
  ## test function returning 'p.value' = log(1-pval) and 'statistic' = log(stat) of independence tests
  
  
  ## TO DO:  - fix option coin
  #          - include SPLITONLY, MIA, ... ?
  
  #ix <- data$zindex[[j]] ### data[[j, type = "index"]]
  #iy <- data$yxindex ### data[["yx", type = "index"]]
  
  Y <- model$estfun  
  
  ## check control arguments
  
  if(length(ctrl$guide_testtype) > 1) ctrl$guide_testtype <- ctrl$guide_testtype[1]
  if(!ctrl$guide_testtype %in% c("max", "sum", "coin")) stop("guide_testtype has to be one of the following options: max, sum, coin")
  
  if(is.null(ctrl$guide_parm)) {
    # indices of residuals to be considered (for intercept and regressors)
    ctrl$guide_parm <- c(1:NCOL(Y))
    #if(length(data$variables$x) > 0) ctrl$guide_parm <- c(1, data$variables$x - length(data$variables$y) + 1)
    #if(length(data$variables$x)== 0) ctrl$guide_parm <- c(1, c(1:length(data$variables$z))[data$variables$z > 0] - length(data$variables$y) + 1)
  }
  
  if(!ctrl$guide_decorrelate %in% c("vcov", "opg", "none")) stop("guide_decorrelate has to be set to one of the following options: none, vcov, opg")
      
  # check whether model$estfun is already decorrelated
  if(!is.null(model$decorrelated)) {
    
    if(!model$decorrelated) {
      # check if argument decorrelate is set to TRUE for testtype max or sum
      if((ctrl$guide_testtype == "sum" | ctrl$guide_testtype == "max") & ctrl$guide_decorrelate == "none") {
        stop("guide_decorrelate has to be set to TRUE for guide_testtype 'sum' or 'max'")
      }
    }
    
    if(model$decorrelated) {
      if(ctrl$guide_decorrelate != "none") stop("scores have already been decorrelated within the trafo function")
      ctrl$guide_decorrelate <- "none"
    }
  } else {
    # check if argument decorrelate is set to TRUE for testtype max or sum
    if((ctrl$guide_testtype == "sum" | ctrl$guide_testtype == "max") & ctrl$guide_decorrelate == "none") {
      stop("guide_decorrelate has to be set to TRUE for guide_testtype 'sum' or 'max'")
    }
  }
  
  if(ctrl$guide_decorrelate != "none") {
    ef <- as.matrix(Y)
    n <- NROW(ef)
    ef <- ef/sqrt(n)
    
    if(is.null(vcov(model$object))) warning("vcov is not stored in model object and therefore not used for decorrelation")
    vcov <- if(ctrl$guide_decorrelate == "vcov" & !is.null(vcov(model$object))) {
      vcov(model$object) * n
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
  
  
  # should scores be weighted ?
  if(is.null(weights) | length(weights) == 0) weights <- rep.int(1, NCOL(Y)) 
  if(!ctrl$weighted.scores && !model$unweighted) Y <- Y/weights  ## FIX ME: influence of weights only on categorization
  x <- data[[j]]
  if(!is.null(subset)) {
    Y <- if(is.vector(Y)) Y[subset] else Y[subset,]
    npar <- NCOL(Y)
    x <- x[subset]
  }
  
  # if all values of the selected covariate are equal return highest possible p.value 
  # and Teststatistic = 0
  if(length(unique(x))<2) return(list(p.value = log(1-1), statistic = log(0)))
  
  
  # categorize residuals
  if(is.null(ctrl$ygroups)){
    # split Y into 2 parts based on whether residuals (here: scores) are positive or negative
    # separately for each parameter
    Ybin <- data.frame(factor(Y[,ctrl$guide_parm[1]]>0))
    colnames(Ybin)[NCOL(Ybin)] <- paste0("rp", ctrl$guide_parm[1])
    
    if(length(ctrl$guide_parm)>1){
      for(k in ctrl$guide_parm[-1]){
        respos <- factor(Y[,k]>0)
        Ybin <- cbind(Ybin, respos)
        colnames(Ybin)[NCOL(Ybin)] <- paste0("rp",k)
      }
    }  
    
  } else {
    
    if(length(ctrl$ygroups)>1) ybreaks <- ctrl$ygroups
    if(length(ctrl$ygroups)==1) ybreaks <- quantile(Y, c(0:ctrl$ygroups)/ctrl$ygroups)
    
    # split Y according to ybreaks
    # separately for each parameter
    Ybin <- data.frame(cut(Y[,1], breaks = ybreaks, labels = c(1:ctrl$ygroups), include.lowest = TRUE))
    colnames(Ybin)[NCOL(Ybin)] <- paste0("res", ctrl$guide_parm[1])
    
    if(length(ctrl$guide_parm)>1){
      for(k in ctrl$guide_parm[-1]){
        respos <- cut(Y[,k], breaks = ybreaks, labels = c(1:ctrl$ygroups), include.lowest = TRUE)
        Ybin <- cbind(Ybin, respos)
        colnames(Ybin)[NCOL(Ybin)] <- paste0("res",k)
      }
    }
  }    
    
    
  # categorize split variable
  if(is.null(ctrl$xgroups)) ctrl$xgroups <- 4
  if(length(ctrl$xgroups)>1) {
    xbreaks <- ctrl$xgroups
  } else {
    xbreaks <- quantile(x, c(0:ctrl$xgroups)/ctrl$xgroups)
  }
  
  if(is.numeric(x)){
    x_cat <- cut(x, breaks = xbreaks, labels = c(1:ctrl$xgroups), include.lowest = TRUE)
  } else {
    x_cat <- x
  }

  
  ## compute curvature test (for each parameter separately)
  ## TO DO: depending on ctrl argument use chisq.test or coin::independence_test
  
  if(ctrl$guide_testtype == "coin"){
    require("coin")
    ip <- new("IndependenceProblem", x=data.frame(x_cat = x_cat), y=Ybin)
    tst_curv <- independence_test(ip, teststat = "quadratic")    # coin:::independence_test(ip, teststat = "quadratic")
    ret <- list(p.value = log(1 - as.numeric(pvalue(tst_curv))), statistic = log(as.numeric(statistic(tst_curv)))) 
    # ret <- list(p.value = log(1 - as.numeric(coin::pvalue(tst_curv))), statistic = log(as.numeric(coin::statistic(tst_curv)))) 
  }
  
  if(ctrl$guide_testtype == "sum" | ctrl$guide_testtype == "max") {
    tst_curv <- chisq.test(x = x_cat, y = Ybin[,1])
    ret <- list(p.value = log(1 - as.numeric(tst_curv$p.value)), statistic = log(as.numeric(tst_curv$statistic)))
    
    if(length(ctrl$guide_parm)>1){
      
      if(ctrl$guide_testtype == "sum"){
        sumstat <- tst_curv$statistic
        for(k in 2:length(ctrl$guide_parm)){
          sumstat <- sumstat + chisq.test(x = x_cat, y = Ybin[,k])$statistic
        }
        p.val <- 1-pchisq(sumstat, df = 3*length(ctrl$guide_parm))
        stat <- sumstat
      }
      
      if(ctrl$guide_testtype == "max"){
        allstat <- numeric(length = length(ctrl$guide_parm))
        allstat[1] <- tst_curv$statistic
        for(k in 2:length(ctrl$guide_parm)){
          allstat[k] <- chisq.test(x = x_cat, y = Ybin[,k])$statistic
        }
        stat <- max(allstat)
        p.val <- 1 - pchisq(stat, df = 3)^length(ctrl$guide_parm)
      }  
      
      ret <- list(p.value = log(1 - as.numeric(p.val)), statistic = log(as.numeric(stat)))
    }  
  }
  
  return(ret)
  
}






##### ctree and mfluc tests with categorization of explanatory variables

## ctree_cat
.ctree_cat_select <- function(...)
  function(model, trafo, data, subset, weights, whichvar, ctrl) {
    args <- list(...)
    ctrl[names(args)] <- args
    partykit:::.select(model, trafo, data, subset, weights, whichvar, ctrl, FUN = .ctree_cat_test)
  }

.ctree_cat_test <- function(model, trafo, data, subset, weights, j, SPLITONLY = FALSE, ctrl) {
  
  ix <- data$zindex[[j]] ### data[[j, type = "index"]]
  iy <- data$yxindex ### data[["yx", type = "index"]]
  Y <- model$estfun
  
  if (!is.null(iy)) {
    stopifnot(NROW(levels(iy)) == (NROW(Y) - 1))
    return(.ctree_cat_test_2d(data = data, j = j, Y = Y, iy = iy, 
                          subset = subset, weights = weights, 
                          SPLITONLY = SPLITONLY, ctrl = ctrl))
  }
  
  stopifnot(NROW(Y) == length(ix))
  
  NAyx <- data$yxmissings ### data[["yx", type = "missings"]]
  NAz <- data$missings[[j]] ### data[[j, type = "missings"]]
  if (ctrl$MIA && (ctrl$splittest || SPLITONLY)) {
    subsetNArm <- subset[!(subset %in% NAyx)]
  } else {
    subsetNArm <- subset[!(subset %in% c(NAyx, NAz))]
  }
  
  return(.ctree_cat_test_1d(data = data, j = j, Y = Y, subset = subsetNArm, 
                        weights = weights, SPLITONLY = SPLITONLY, ctrl = ctrl))
}



.ctree_cat_test_1d <- function(data, j, Y, subset, weights, SPLITONLY = FALSE, ctrl) {
  
  x <- data[[j]]
  
  # categorize split variable
  if(!is.numeric(x) & !is.null(ctrl$xgroups)) warning("xgroups is only used for categorization of numerical split variables")
  if(is.numeric(x)){
    if(is.null(ctrl$xgroups)) ctrl$xgroups <- 4
    if(length(ctrl$xgroups)>1) {
      xbreaks <- ctrl$xgroups
    } else {
      xbreaks <- quantile(x, c(0:ctrl$xgroups)/ctrl$xgroups)
    }
    x <- cut(x, breaks = xbreaks, labels = c(1:ctrl$xgroups), include.lowest = TRUE)
  }
  
  MIA <- FALSE
  if (ctrl$MIA) {
    NAs <- data$missings[[j]] ### data[[j, type = "missings"]]
    MIA <- (length(NAs) > 0)
  }
  
  ### X for (ordered) factors is always dummy matrix
  if (is.factor(x) || is.ordered(x))
    X <- inum::enum(x)     ## FIX ME: split variable has been categorized before, apply enum here (in extree_data inum is applied)
    #X <- data$zindex[[j]] ### data[[j, type = "index"]]
  
  scores <- data[[j, type = "scores"]]
  ORDERED <- is.ordered(x) || is.numeric(x)
  
  ux <- Xleft <- Xright <- NULL
  
  if (ctrl$splittest || SPLITONLY) {
    MAXSELECT <- TRUE
    if (is.numeric(x)) {
      X <- data$zindex[[j]] ###data[[j, type = "index"]]
      ux <- levels(X)
    }
    if (MIA) {
      Xlev <- attr(X, "levels")
      Xleft <- X + 1L
      Xleft[NAs] <- 1L
      Xright <- X
      Xright[NAs] <- as.integer(length(Xlev) + 1L)
      attr(Xleft, "levels") <- c(NA, Xlev)
      attr(Xright, "levels") <- c(Xlev, NA)
    } 
  } else {
    MAXSELECT <- FALSE
    if (is.numeric(x)) {
      if (storage.mode(x) == "double") {
        X <- x
      } else {
        X <- as.double(x) ### copy when necessary
      }
    }
    MIA <- FALSE
  }
  cluster <- data[["(cluster)"]]
  
  partykit:::.ctree_test_internal(x = x, X = X, ix = NULL, Xleft = Xleft, Xright = Xright, 
                       ixleft = NULL, ixright = NULL, ux = ux, scores = scores, 
                       j = j, Y = Y, iy = NULL, subset = subset, weights = weights, 
                       cluster = cluster, MIA = MIA, SPLITONLY = SPLITONLY, 
                       MAXSELECT = MAXSELECT, ORDERED = ORDERED, ctrl = ctrl)
}


.ctree_cat_test_2d <- function(data, Y, iy, j, subset, weights, SPLITONLY = FALSE, ctrl) {
  
  x <- data[[j]]
  
  # categorize split variable
  if(!is.numeric(x) & !is.null(ctrl$xgroups)) warning("xgroups is only used for categorization of numerical split variables")
  if(is.numeric(x)){
    if(is.null(ctrl$xgroups)) ctrl$xgroups <- 4
    if(length(ctrl$xgroups)>1) {
      xbreaks <- ctrl$xgroups
    } else {
      xbreaks <- quantile(x, c(0:ctrl$xgroups)/ctrl$xgroups)
    }
    x <- cut(x, breaks = xbreaks, labels = c(1:ctrl$xgroups), include.lowest = TRUE)
  }
  
  ix <- inum::enum(x)     ## FIX ME: split variable has been categorized before, apply enum here (in extree_data inum is applied)
  #ix <- data$zindex[[j]] ### data[[j, type = "index"]]
  ux <- attr(ix, "levels")
  
  MIA <- FALSE
  if (ctrl$MIA) MIA <- any(ix[subset] == 0)
  
  ### X for (ordered) factors is always dummy matrix
  if (is.factor(x) || is.ordered(x))
    X <- integer(0)
  
  scores <- data[[j, type = "scores"]]
  ORDERED <- is.ordered(x) || is.numeric(x)
  
  if (ctrl$splittest || SPLITONLY) {
    MAXSELECT <- TRUE
    X <- integer(0)
    
    if (MIA) {
      Xlev <- attr(ix, "levels")
      ixleft <- ix + 1L
      ixright <- ix
      ixright[ixright == 0L] <- as.integer(length(Xlev) + 1L)
      attr(ixleft, "levels") <- c(NA, Xlev)
      attr(ixright, "levels") <- c(Xlev, NA)
      Xleft <- Xright <- X
    } 
  } else {
    MAXSELECT <- FALSE
    MIA <- FALSE
    if (is.numeric(x))
      X <- matrix(c(0, as.double(attr(ix, "levels"))), ncol = 1)
  }
  cluster <- data[["(cluster)"]]
  
  partykit:::.ctree_test_internal(x = x, X = X, ix = ix, Xleft = Xleft, Xright = Xright, 
                       ixleft = ixleft, ixright = ixright, ux = ux, scores = scores, 
                       j = j, Y = Y, iy = iy, subset = subset, weights = weights, 
                       cluster = cluster, MIA = MIA, SPLITONLY = SPLITONLY, 
                       MAXSELECT = MAXSELECT, ORDERED = ORDERED, ctrl = ctrl)
}



## mfluc_cat
.mfluc_cat_select <- function(...)
  function(model, trafo, data, subset, weights, whichvar, ctrl) {
    args <- list(...)
    ctrl[names(args)] <- args
    partykit:::.select(model, trafo, data, subset, weights, whichvar, ctrl, FUN = .mfluc_cat_test)
  }



.mfluc_cat_test <- function(model, trafo, data, subset, weights, j, SPLITONLY, ctrl)
{
  
  stopifnot(!SPLITONLY)
  stopifnot(is.null(data[["yx", type = "index"]]))
  z <- data[[j]][subset]
  
  # categorize split variable
  if(!is.numeric(z) & !is.null(ctrl$xgroups)) warning("xgroups is only used for categorization of numerical split variables")
  if(is.numeric(z)){
    if(is.null(ctrl$xgroups)) ctrl$xgroups <- 4
    if(length(ctrl$xgroups)>1) {
      xbreaks <- ctrl$xgroups
    } else {
      xbreaks <- quantile(z, c(0:ctrl$xgroups)/ctrl$xgroups)
    }
    z <- cut(z, breaks = xbreaks, labels = c(1:ctrl$xgroups), include.lowest = TRUE)
  }
  
  estfun <- model$estfun[subset,,drop = FALSE]
  cluster <- data[["(cluster)"]][subset]
  if(length(weights) == 0) {
    weights <- rep(1, NROW(estfun))
  } else {
    weights <- weights[subset]
  }
  obj <- model$obj
  
  if(length(unique(z)) < 2L) return(list(statistic = NA, p.value = NA))
  
  ## set up return values
  m <- NCOL(z)
  pval <- rep.int(NA_real_, m)
  stat <- rep.int(0, m)
  ifac <- rep.int(FALSE, m)
  
  ## estimating functions (dropping zero weight observations)
  process <- as.matrix(estfun)
  ww0 <- (weights > 0)
  process <- process[ww0, , drop = FALSE]
  if(!is.null(cluster)) cluster <- droplevels(cluster[ww0]) 
  weights <- weights[ww0]
  z <- z[ww0]
  k <- NCOL(process)
  n <- NROW(process)
  nobs <- if(ctrl$caseweights && any(weights != 1L)) sum(weights) else n
  
  
  ## stop if all collumn values in process are the same
  if(all(apply(process, 2, function(x) length(unique(x))) == 1)) 
    return(list(statistic = NA, p.value = NA))
  
  ## scale process
  process <- process/sqrt(nobs)
  vcov <- ctrl$vcov
  if(is.null(obj)) vcov <- "opg"
  if(vcov != "opg") {
    bread <- vcov(obj) * nobs
  }
  if(vcov != "info") {
    ## correct scaling of estfun for variance estimate:
    ## - caseweights=FALSE: weights are integral part of the estfun -> squared in estimate
    ## - caseweights=TRUE: weights are just a factor in variance estimate -> require division by sqrt(weights)
    meat <- if(is.null(cluster)) {
      crossprod(if(ctrl$caseweights) process/sqrt(weights) else process)
    } else {
      crossprod(as.matrix(apply(if(ctrl$caseweights) process/sqrt(weights) else process, 2L, tapply, cluster, sum)))
    }
  }
  ## from strucchange
  root.matrix <- function(X) {
    if((ncol(X) == 1L)&&(nrow(X) == 1L)) return(sqrt(X)) else {
      X.eigen <- eigen(X, symmetric = TRUE)
      if(any(X.eigen$values < 0)) stop("Matrix is not positive semidefinite")
      sqomega <- sqrt(diag(X.eigen$values))
      V <- X.eigen$vectors
      return(V %*% sqomega %*% t(V))
    }
  }
  J12 <- root.matrix(switch(vcov,
                            "opg" = chol2inv(chol(meat)),
                            "info" = bread,
                            "sandwich" = bread %*% meat %*% bread
  ))
  process <- t(J12 %*% t(process))  
  
  
  ## get critical values for supLM statistic
  from <- if(ctrl$trim > 1) ctrl$trim else ceiling(nobs * ctrl$trim)
  from <- max(from, ctrl$minbucket)
  to <- nobs - from
  lambda <- ((nobs - from) * to)/(from * (nobs - to))
  
  beta <- partykit:::mob_beta_suplm
  logp.supLM <- function(x, k, lambda)
  {
    if(k > 40L) {
      ## use Estrella (2003) asymptotic approximation
      logp_estrella2003 <- function(x, k, lambda)
        -lgamma(k/2) + k/2 * log(x/2) - x/2 + log(abs(log(lambda) * 
                                                        (1 - k/x) + 2/x))
      ## FIXME: Estrella only works well for large enough x
      ## hence require x > 1.5 * k for Estrella approximation and
      ## use an ad hoc interpolation for larger p-values
      p <- ifelse(x <= 1.5 * k, 
                  (x/(1.5 * k))^sqrt(k) * logp_estrella2003(1.5 * k, k, lambda), 
                  logp_estrella2003(x, k, lambda))
    } else {
      ## use Hansen (1997) approximation
      nb <- ncol(beta) - 1L
      tau <- if(lambda < 1) lambda else 1/(1 + sqrt(lambda))
      beta <- beta[(((k - 1) * 25 + 1):(k * 25)),]
      dummy <- beta[,(1L:nb)] %*% x^(0:(nb-1))
      dummy <- dummy * (dummy > 0)
      pp <- pchisq(dummy, beta[,(nb+1)], lower.tail = FALSE, log.p = TRUE)
      if(tau == 0.5) {
        p <- pchisq(x, k, lower.tail = FALSE, log.p = TRUE)
      } else if(tau <= 0.01) {
        p <- pp[25L]
      } else if(tau >= 0.49) {
        p <- log((exp(log(0.5 - tau) + pp[1L]) + exp(log(tau - 0.49) + pchisq(x, k, lower.tail = FALSE, log.p = TRUE))) * 100)
        ## if p becomes so small that 'correct' weighted averaging does not work, resort to 'naive' averaging
        if(!is.finite(p)) p <- mean(c(pp[1L], pchisq(x, k, lower.tail = FALSE, log.p = TRUE)))
      } else {
        taua <- (0.51 - tau) * 50
        tau1 <- floor(taua)
        p <- log(exp(log(tau1 + 1 - taua) + pp[tau1]) + exp(log(taua-tau1) + pp[tau1 + 1L]))
        ## if p becomes so small that 'correct' weighted averaging does not work, resort to 'naive' averaging
        if(!is.finite(p)) p <- mean(pp[tau1 + 0L:1L])
      }
    }
    return(as.vector(p))
  }
  
  ## compute statistic and p-value
  if(is.factor(z)) {
    oi <- order(z)
    proci <- process[oi, , drop = FALSE]
    
    ifac <- TRUE
    iord <- is.ordered(z) & (ctrl$ordinal != "chisq")
    
    ## order partitioning variable
    z <- z[oi]
    # re-apply factor() added to drop unused levels
    z <- factor(z, levels = unique(z))
    # compute segment weights
    segweights <- if(ctrl$caseweights) tapply(weights[oi], z, sum) else table(z)
    segweights <- as.vector(segweights)/nobs
    
    # compute statistic only if at least two levels are left
    if(length(segweights) < 2L) {
      stat <- 0
      pval <- NA_real_
    } else if(iord) {
      proci <- apply(proci, 2L, cumsum)
      tt0 <- head(cumsum(table(z)), -1L)
      tt <- head(cumsum(segweights), -1L)
      if(ctrl$ordinal == "max") {
        stat <- max(abs(proci[tt0, ] / sqrt(tt * (1-tt))))
        pval <- log(as.numeric(1 - mvtnorm::pmvnorm(
          lower = -stat, upper = stat,
          mean = rep(0, length(tt)),
          sigma = outer(tt, tt, function(x, y)
            sqrt(pmin(x, y) * (1 - pmax(x, y)) / 
                   ((pmax(x, y) * (1 - pmin(x, y))))))
        )^k))
      } else {
        proci <- rowSums(proci^2)
        stat <- max(proci[tt0] / (tt * (1-tt)))
        pval <- log(strucchange::ordL2BB(segweights, nproc = k, nrep = ctrl$nrep)$computePval(stat, nproc = k))
      }
    } else {      
      stat <- sum(sapply(1L:k, function(j) (tapply(proci[,j], z, sum)^2)/segweights))
      pval <- pchisq(stat, k*(length(levels(z))-1), log.p = TRUE, lower.tail = FALSE)
    }
  } else {
    oi <- if(ctrl$breakties) {
      mm <- sort(unique(z))
      mm <- ifelse(length(mm) > 1L, min(diff(mm))/10, 1)
      order(z + runif(length(z), min = -mm, max = +mm))
    } else {
      order(z)
    }
    proci <- process[oi, , drop = FALSE]
    proci <- apply(proci, 2L, cumsum)
    tt0 <- if(ctrl$caseweights && any(weights != 1L)) cumsum(weights[oi]) else 1:n
    stat <- if(from < to) {
      xx <- rowSums(proci^2)
      xx <- xx[tt0 >= from & tt0 <= to]
      tt <- tt0[tt0 >= from & tt0 <= to]/nobs
      max(xx/(tt * (1 - tt)))	  
    } else {
      0
    }
    pval <- if(from < to) logp.supLM(stat, k, lambda) else NA
  }
  
  ## return version of pvalue that .extree_node deals with
  rval <- list(statistic = log(stat), p.value = log1p(-exp(pval)))
  return(rval)
}



if(FALSE){
  # use different version of .select than partykit:::select for GUIDE as selectfun depending on 
  # whether or not interaction tests should be made
  # (returns p.values from curvature and interaction tests instead of p.value and teststatistic)
  .select_int <- function(model, trafo, data, subset, weights, whichvar, ctrl, FUN) {
    ret <- list(criteria = matrix(NA, nrow = 2L, ncol = ncol(model.frame(data))))
    rownames(ret$criteria) <- c("statistic", "p.value")
    colnames(ret$criteria) <- names(model.frame(data))
    if (length(whichvar) == 0) return(ret)
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
  
  
  .guide_test_int <- function(model, trafo, data, subset, weights, j, SPLITONLY = FALSE, ctrl) {
    
    ## TO DO: include SPLITONLY, MIA, ... ?
    ##        add ctrl argument 'guide_parm', a vector of indices of the parameters that should be considered (only within the regressors, starting at 1)
    
    #ix <- data$zindex[[j]] ### data[[j, type = "index"]]
    #iy <- data$yxindex ### data[["yx", type = "index"]]
    
    Y <- model$estfun  ## model from distfit returns wheigthed scores
    if(ctrl$guide_unweighted && !model$unweighted) Y <- Y/weights  ## FIX ME: influence of weights only on categorization
    x <- data[[j]]
    if(!is.null(subset)) {
      Y <- if(is.vector(Y)) Y[subset] else Y[subset,]
      npar <- NCOL(Y)
      x <- x[subset]
    }
    
    # if all values of the selected covariate are equal return highest possible p.value
    if(length(unique(x))<2) return(c(p.min = 1, p.curv = 1))
    
    # split Y into 2 parts based on whether residuals (here: scores) are positive or negative
    # separately for each parameter
    Ybin <- data.frame(factor(Y[,ctrl$guide_parm[1]]>0))
    colnames(Ybin)[NCOL(Ybin)] <- paste0("rp", ctrl$guide_parm[1])
    if(length(ctrl$guide_parm)>1){
      for(k in ctrl$guide_parm[-1]){
        respos <- (Y[,k]>0)
        Ybin <- cbind(Ybin, respos)
        colnames(Ybin)[NCOL(Ybin)] <- paste0("rp",k)
      }
    }  
    
    if(is.numeric(x)){
      x_cat <- cut(x, breaks = quantile(x, c(0:4/4)), labels = c(1:4), include.lowest = TRUE)
    } else {
      x_cat <- x
    }
    
    ## compute curvature test (for each parameter separately)
    p.curv <- chisq.test(x = x_cat, y = Ybin[,1])$p.value
    if(length(ctrl$guide_parm)>1){
      for(k in 2:length(ctrl$guide_parm)){
        p <- chisq.test(x = x_cat, y = Ybin[,k])$p.value
        if(p < p.curv) p.curv <- p
      }
    }  
    
    p.min <- p.curv
    
    
    ## compute interaction test (for each parameter and for each of the other covariates separately)
    # only keep test if p.value is smaller than the one resulting from the curvature test
    
    ## FIX ME: select correct covariates (data$variables$z, data$variables$x, ctrl$partyvars, ctrl_guide_parm)
    # only select those other covariates which are also in partyvars
    ix_others <- c(1:NCOL(data$data))[data$variables$z + ctrl$partyvars == 2]
    ix_others <- ix_others[!ix_others == j]
    
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
  
  
  
  # wrapper function to use lm() as trafo
  ## FIX ME: get the right data set (d, data, ... ?)
  ytrafo <- function(subset, weights, estfun = TRUE, object = TRUE, info = NULL) {
    
    formula <- formula(d$terms$all)
    sdata <- d$data[subset,]
    
    ## FIX ME: error in lm if weights are handed over
    ## FIX ME: scores with or without weights?
    #subweights <- if(is.null(weights) || (length(weights)==0L)) weights else weights[subset]
    
    model <- lm(formula = formula, data = sdata) #, weights = subweights) error in lm if weights are handed over
    
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
}  