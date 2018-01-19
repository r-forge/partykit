#### TO DO:
# - add control arguments for .guide_test:
#       guide_interaction = FALSE 
#       parm (indices of residuals/coefficients to be tested)
#       guide_unwheighted = FALSE 
# - if(guide_interaction) choose .select_g within .guide_select, else .select as usual


.guide_select <- function(guide_parm = NULL,  # a vector of indices of parameters for which estfun should be considered
                          guide_testtype = c("max", "sum", "coin"),
                          interaction = FALSE,
                          guide_decorrelate = "vcov",   # needs to be set to other than "none" for testtype max and sum 
                                                        # unless ytrafo returns decorrelated scores
                          # FIX ME: c("none","vcov","opg")
                          xgroups = NULL,  # number of categories for split variables (optionally breaks can be handed over)
                          ygroups = NULL,  # number of categories for scores (optionally breaks can be handed over)
                          weighted.scores = FALSE   # logical, should scores be weighted
                          ) {
  function(model, trafo, data, subset, weights, whichvar, ctrl) {
    args <- list(guide_parm = guide_parm, guide_testtype = guide_testtype, guide_decorrelate = guide_decorrelate, 
                 xgroups = xgroups, ygroups = ygroups, weighted.scores = weighted.scores)
    ctrl[names(args)] <- args
    if(interaction){
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
  
  
  ## TO DO:  - include input arguments xgroups and ygroups
  #          - include SPLITONLY, MIA, ... ?
  
  #ix <- data$zindex[[j]] ### data[[j, type = "index"]]
  #iy <- data$yxindex ### data[["yx", type = "index"]]
  
  Y <- model$estfun  
  
  ## check control arguments
  
  if(length(ctrl$guide_testtype) > 1) ctrl$guide_testtype <- ctrl$guide_testtype[1]
  if(!ctrl$guide_testtype %in% c("max", "sum", "coin")) stop("guide_testtype has to be one of the following options: max, sum, coin")
  
  if(is.null(ctrl$guide_parm)) ctrl$guide_parm <- c(1:length(data$variables$z))[data$variables$z > 0]
  
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
    n <- NROW(ef)
    ef <- ef/sqrt(n)
    
    vcov <- if(ctrl$guide_decorrelate == "vcov") {
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
  
  
  # should scores be weighted ?
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
  
  if(is.numeric(x)){
    x_cat <- cut(x, breaks = quantile(x, c(0:4/4)), labels = c(1:4), include.lowest = TRUE)
  } else {
    x_cat <- x
  }
  
  
  ## compute curvature test (for each parameter separately)
  ## TO DO: depending on ctrl argument use chisq.test or coin::independence_test
  
  if(ctrl$guide_testtype == "coin"){
    tst_curv <- coin:::independence_test(x_cat ~ Ybin)
    ret <- list(p.value = log(1 - as.numeric(coin::pvalue(tst_curv))), statistic = log(as.numeric(coin::statistic(tst_curv)))) 
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
    ##        add ctrl argument 'guide_parm', a vector of indices of the parameters that should be considered
    
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