distexforest <- function
(
    formula,
    data,   
    subset, 
    na.action = na.pass,
    weights,
    offset, 
    cluster,
    family = NO(),
    strata,
    control = distextree_control(
        teststat = "quad", testtype = "Univ", mincriterion = 0,
        saveinfo = FALSE, minsplit = 20, minbucket = 7, splittry = 2, ...),
    ntree = 500L, 
    fit.par = FALSE,
    perturb = list(replace = FALSE, fraction = 0.632),
    mtry = ceiling(sqrt(nvar)), 
    applyfun = NULL,
    cores = NULL, 
    trace = FALSE,
    ...
) {
   
    ## Get original formula
    oformula <- as.formula(formula)

    ## Get family
    if(!inherits(family, "disttree.family"))  #FIXME: (LS) Better way to extract prepared family?
        family <- distfamily(family)

    ## Set minsize to 10 * number of parameters, if NULL
    if (is.null(control$minbucket) | is.null(control$minsplit)) {
      n_coef <- length(family$link)
      minsize <- as.integer(ceiling(10L * n_coef)) #FIXME: (ML) Adapt for multivariate repsone.
      if (is.null(control$minbucket)) control$minbucket <- minsize
      if (is.null(control$minsplit)) control$minsplit <- minsize
    }
   
    ### get the call and the calling environment for .urp_tree
    call <- match.call(expand.dots = FALSE)
    oweights <- NULL
    if (!missing(weights))
        oweights <- weights
    m <- match(c("formula", "data", "subset", "na.action", "offset", "cluster", 
                 "family", "control"), names(call), 0L)
    distextreecall <- call[c(1L, m)]
    distextreecall$doFit <- FALSE
    if (!is.null(oweights))
        distextreecall$weights <- 1:NROW(oweights)
    distextreecall$control <- control ### put ... into ctree_control()
    distextreecall[[1L]] <- quote(disttree::distextree)
    tree <- eval(distextreecall, parent.frame())

    control$update <- TRUE

    d <- tree$d
    updatefun <- tree$update

    nvar <- sum(d$variables$z > 0)
    control$mtry <- mtry
    control$applyfun <- lapply
 
    strata <- d[["(strata)"]]
    if (!is.null(strata)) {
        if (!is.factor(strata)) stop("strata is not a single factor")
    }
    
    probw <- NULL
    iweights <- model.weights(model.frame(d))
    if (!is.null(oweights)) {
        if (is.matrix(oweights)) {
            weights <- oweights[iweights,,drop = FALSE]
        } else {
            weights <- oweights[iweights]
        }
    } else {
        weights <- NULL
    }
    rm(oweights)
    rm(iweights)
    N <- nrow(model.frame(d))
    rw <- NULL
    if (!is.null(weights)) {
        if (is.matrix(weights)) {
            if (ncol(weights) == ntree && nrow(weights) == N) {
                rw <- unclass(as.data.frame(weights))
                rw <- lapply(rw, function(w) 
                    rep(1:length(w), w))
                weights <- integer(0)
            } else {
                stop(sQuote("weights"), "argument incorrect")
            }
        } else {
            probw <- weights / sum(weights)
        }
    } else {
        weights <- integer(0)
    }

    idx <- partykit:::.start_subset(d)
    if (is.null(rw)) {
        if (is.null(strata)) {
            size <- N
            if (!perturb$replace) size <- floor(size * perturb$fraction)
            rw <- replicate(ntree, 
                            sample(idx, size = size, 
                                   replace = perturb$replace, prob = probw),
                            simplify = FALSE)
        } else {
            frac <- if (!perturb$replace) perturb$fraction else 1
            rw <- replicate(ntree, function() 
                  do.call("c", tapply(idx, strata, 
                          function(i) 
                              sample(i, size = length(i) * frac, 
                                     replace = perturb$replace, prob = probw[i]))))
        }
    }

    ## apply infrastructure for determining split points
    ## use RNGkind("L'Ecuyer-CMRG") to make this reproducible
    if (is.null(applyfun)) {
        applyfun <- if(is.null(cores)) {
            lapply  
        } else {
            function(X, FUN, ...)
                parallel::mclapply(X, FUN, ..., mc.set.seed = TRUE, mc.cores = cores)
        }
    }

    trafo <- updatefun(sort(rw[[1]]), integer(0), control, doFit = FALSE)
    if (trace) pb <- txtProgressBar(style = 3) 
    forest <- applyfun(1:ntree, function(b) {
        if (trace) setTxtProgressBar(pb, b/ntree)
        ret <- updatefun(sort(rw[[b]]), integer(0), control)
        # trafo <<- ret$trafo
        ret$nodes  #FIXME: (ML) If control$saveinfo = TRUE, should we really return all in node$info? 
    })
    if (trace) close(pb)

    fitted <- data.frame(idx = 1:N)  
    mf <- model.frame(d)
    fitted[[2]] <- mf[, d$variables$y, drop = TRUE]
    names(fitted)[2] <- "(response)"
    if (length(weights) > 0)
        fitted[["(weights)"]] <- weights

    ### turn subsets in weights (maybe we can avoid this?)
    rw <- lapply(rw, function(x) as.integer(tabulate(x, nbins = length(idx))))

    control$applyfun <- applyfun

    ret <- partykit:::constparties(nodes = forest, data = mf, weights = rw,
                        fitted = fitted, terms = d$terms$all,
                        info = list(call = match.call(), 
                                    control = control,
                                    family = family))

    ret$info$call$formula <- oformula   # FIXME: (ML) Tweak to get real formula for logLik in varimp()

    ret$trafo <- trafo
    ret$predictf <- d$terms$z
    class(ret) <- c("distexforest", class(ret))

    if(fit.par) {
      ## adaptive local log-likelihood
      ## calculate fitted values, fitted distribution parameters, 
      ##   loglikelihood for every observation

      np <- length(family$link)

      fitted.par <- data.frame(matrix(0, nrow = nrow(data), ncol = np))
      loglik <- data.frame(idx = 1:nrow(data))

      # extract weights
      w <- partykit::predict.cforest(ret, type = "weights", OOB = FALSE)  #FIXME: (LS) allow for additional argument OOB in distexforest

      Y <- ret$fitted$`(response)`

      for(i in 1:nrow(data)){
        wi <- w[,i]
        # personalized model for observation data[i,]
        pm <-  disttree::distexfit(Y, family = family, weights = wi, vcov = FALSE, 
          optim.control = control$optim.control, ...)
        fitted.par[i,] <- coef(pm, type = "parameter")
        loglik[i,] <- if(is.function(pm$ddist)) pm$ddist(Y[i], log = TRUE) else NA
      }

      if(is.null(weights) || (length(weights)==0L || is.function(weights))) 
        weights <- numeric(nrow(data)) + 1

      ret$fitted$`(weights)` <- weights

      names(fitted.par) <- names(coef(pm, type = "parameter"))
      ret$fitted.par <- fitted.par

      ret$loglik <- sum(loglik)
    }
    return(ret)
}



predict.distexforest <- function(object, newdata = NULL, 
                                 type = c("parameter", "response", "weights", "node"), 
                                 OOB = FALSE, 
                                 scale = TRUE, ...) {

    responses <- object$fitted[["(response)"]]
    forest <- object$nodes
    nd <- object$data
    vmatch <- 1:ncol(nd)
    NOnewdata <- TRUE
    if (!is.null(newdata)) {
        factors <- which(sapply(nd, is.factor))
        xlev <- lapply(factors, function(x) levels(nd[[x]]))
        names(xlev) <- names(nd)[factors]
        nd <- model.frame(object$predictf, ### all variables W/O response
                          data = newdata, na.action = na.pass, xlev = xlev)
        OOB <- FALSE
        vmatch <- match(names(object$data), names(nd))
        NOnewdata <- FALSE
    }
    nam <- rownames(nd)

    type <- match.arg(type)

    ### return terminal node ids for data or newdata
    if (type == "node")
        return(lapply(forest, fitted_node, data = nd, vmatch = vmatch, ...))

    ### extract weights
    rw <- object$weights

    w <- 0L

    applyfun <- lapply
    if (!is.null(object$info))
        applyfun <- object$info$control$applyfun

    fdata <- lapply(forest, fitted_node, data = object$data, ...)
    if (NOnewdata && OOB) {
        fnewdata <- list()
    } else {
        fnewdata <- lapply(forest, fitted_node, data = nd, vmatch = vmatch, ...)
    }

    w <- partykit:::.rfweights(fdata, fnewdata, rw, scale)

#    for (b in 1:length(forest)) {
#        ids <- nodeids(forest[[b]], terminal = TRUE)
#        fnewdata <- fitted_node(forest[[b]], nd, vmatch = vmatch, ...)
#        fdata <- fitted_node(forest[[b]], object$data, ...)
#        tw <- rw[[b]]
#        pw <- sapply(ids, function(i) tw * (fdata == i))
#        ret <- pw[, match(fnewdata, ids), drop = FALSE]
#        ### obs which are in-bag for this tree don't contribute
#        if (OOB) ret[,tw > 0] <- 0
#        w <- w + ret
#    }
#
#    #w <- Reduce("+", bw)
#    if (!is.matrix(w)) w <- matrix(w, ncol = 1)

    colnames(w) <- nam
    rownames(w) <- rownames(responses)
    
    if (type == "weights") {
        ret <- w
        return(ret)
    }
    
    if(type == "parameter" | type == "response") {
      
      if(NOnewdata & type == "parameter" & OOB == FALSE & !is.null(object$fitted.par)) return(object$fitted.par)
        
      family <- object$info$family
      # if(!inherits(family, "disttree.family"))  
      #   family <- distfamily(family)
      np <- length(family$link)
      
      if(type == "parameter") pred <- matrix(0, nrow = nrow(nd), ncol = np)
      if(type == "response") pred <- matrix(0, nrow = nrow(nd), ncol = 1)
      # loglik <- data.frame(idx = 1:nrow(data))  ## FIXME: should we provide loglik for newdata?
      
      for(i in 1:nrow(nd)){
        wi <- w[,i]
        # personalized model for observation nd[i,]
        pm <-  disttree::distexfit(responses, family = family, weights = wi, vcov = FALSE, 
                                   optim.control = object$info$control$optim.control, ...)
        pred[i,] <- predict.distexfit(pm, type = type)
        # loglik[i,] <- if(is.function(pm$ddist)) pm$ddist(responses[i], log = TRUE) else NA
      }
      
      if(type == "parameter") {
        pred <- data.frame(pred)
        names(pred) <- names(coef.distexfit(pm, type = "parameter"))
      }
      if(type == "response") {
        pred <- data.frame(pred)
        names(pred) <- "(fitted.response)"
      }
      return(pred)
    }
}


## TODO: how to deal with NAs in fitted parameters (possibly due to all zero weights)
logLik.distexforest <- function(object, newdata = NULL, weights = NULL, ...){
  
  pred.par <- predict.distexforest(object = object, newdata = newdata, type = "parameter")
  
  if(is.null(newdata)) {
    responses <- object$fitted[["(response)"]]
  } else {
    ## FIXME: can we extract the response without building the whole data.frame?
    nd <- object$data
    vmatch <- 1:ncol(nd)
    factors <- which(sapply(nd, is.factor))
    xlev <- lapply(factors, function(x) levels(nd[[x]]))
    names(xlev) <- names(nd)[factors]
    formula <- if(is.name(object$info$call$formula)) eval(object$info$call$formula) else object$info$call$formula
    nd <- model.frame(formula, ## FIXME: use formula with response only
                      data = newdata, na.action = na.pass, xlev = xlev)
    responses <- nd[,as.character(formula[[2]])]
  }

  if(is.null(weights)) weights <- rep.int(1, NROW(responses))
  
  ## FIXME: check format returned by linkfun / linkinv (in case of multidimensional data)
  ll <- sapply(1:NROW(responses),  function(i) object$info$family$ddist(responses[i], 
           eta = as.numeric(object$info$family$linkfun(pred.par[i,])),
           log = TRUE, weights = weights[i]))
  ll <- sum(ll, na.rm = TRUE)
  return(structure(ll, df = NA, class = "logLik"))
}



## varimp
varimp.distexforest <- function(object, nperm = 1L, ...){

  # define function for parallelization
  applyfun <- function(X, FUN, ...) parallel::mclapply(X, FUN, ..., mc.cores = pmax(1, parallel::detectCores() - 1))
  
  # function to permutate chosen variable and then calculate mean loglikelihood
  riskfun <- function(permute = NULL, newdata = object$data) {
    if(!is.null(permute)) newdata[[permute]] <- sample(newdata[[permute]])
    logLik(object, newdata = newdata)
  }

  # apply for all covariates except for dswrf_mean_min and 
  # dswrf_sprd_min (columns 30 and 33) as they are always 0
  
  # using only one core
  # risk_all <- replicate(nperm, sapply(c(5:29, 31, 32, 34: ncol(testdata)), meancrps))
  # risk <- rowMeans(risk_all)

  ## FIXME: better way to get dataset without response
  splitvar <- names(object$data) %in% attr(object$predictf, "term.labels")
  splitid <- c(1:NCOL(object$data))[splitvar]
  
  # or parallel
  risklist <- applyfun(1:nperm, 
                       function(i){
                         set.seed(i)
                         sapply(splitid, riskfun)
                       })
  risk <- Reduce("+", risklist) / length(risklist)
  
  names(risk) <- names(object$data)[splitid]
  vimp <- risk - riskfun(newdata = object$data)
  vimp <- sort(vimp, decreasing = TRUE)

  return(vimp)
}





## copied methods from cforest
model.frame.distexforest <- function(formula, ...) {
    class(formula) <- "party"
    model.frame(formula, ...)
}

### FIXME: (ML) does not get a full distextree object
gettree.distexforest <- function(object, tree = 1L, ...) {
    consttree <- partykit:::gettree.cforest(object, tree = tree, ...)
    #d.response <- consttree$fitted$`(response)`[consttree$fitted$`(weights)` == 1]
    d <- object$data[consttree$fitted$`(weights)` == 1,]
    ctrl <- object$info$control
    ctrl$saveinfo <- TRUE
    newtree <- distextree(object$info$call$formula, 
                          data = d, control = ctrl)
    return(newtree)
}
    

## alternative version of gettree.distexforest
#gettree.distexforest <- function(object, tree = 1L, ...) {
#  consttree <- partykit:::gettree.cforest(object, tree = tree, ...)
#  coefs <- matrix(ncol = length(object$info$family$link), nrow = width(consttree))
#  terminal <- 1
#  terminal.ids <- NULL
#  for(i in 1:length(consttree)){
#    if(is.terminal(consttree[[i]]$node)){
#      d <- consttree$fitted$`(response)`[consttree$fitted$`(weights)` == 1 & consttree$fitted$`(fitted)` == i]
#      cf <- distexfit(d, family = object$info$family)
#      # t[[i]]$node$info$coefficients <- cf$par
#      coefs[terminal, ] <- cf$par
#      terminal.ids <- c(terminal.ids, i)
#      terminal <- terminal+1
#    }
#  }
#  rownames(coefs) <- as.character(terminal.ids)
#  colnames(coefs) <- names(cf$par)
# return(list(treestructure = consttree,
#              coefficients = coefs))
#}
  