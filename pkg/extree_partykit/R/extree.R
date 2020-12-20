extree <- function(data, trafo, control = extree_control(...),  converged = NULL, ...) {
    
    ## check / preprocess extree data
    subset <- .start_subset(data = data)
    weights <- model.weights(model.frame(data))
    
    
    ## trafo preprocessing
    mytrafo <- function(subset, weights, info = NULL, estfun = TRUE, object = TRUE) {
        trafo(subset, data = data, weights, info = NULL, estfun = TRUE, object = TRUE)
    }
    
    ## TODO: converged preprocessing (if needed)
    
    ## set up trafo
    update <- function(subset, weights, control, doFit = TRUE) {
        extree_fit(data = data, trafo = mytrafo, converged = converged,
            partyvars = data$variables$z, subset = subset,
            weights = weights, ctrl = control, doFit = doFit)  
    }
    
    ## fit
    tree <- update(subset = subset, weights = weights, control = control)
    
    ## TODO: prepare extree object
    
}

extree_fit <- function(data, trafo, converged, varselect = ctrl$varselect, 
                       splitselect = ctrl$splitselect, svarselect = ctrl$svarselect, 
                       ssplitselect = ctrl$ssplitselect, partyvars, subset, weights, ctrl, doFit = TRUE) {
  ret <- list()
  
  ### <FIXME> use data$vars$z as default for partyvars </FIXME>
  ### <FIXME> try to avoid doFit </FIXME>
  
  nf <- names(formals(trafo))
  if (all(c("subset", "weights", "info", "estfun", "object") %in% nf)) {
    mytrafo <- trafo
  } else {
    stopifnot(all(c("y", "x", "offset", "weights", "start") %in% nf))
    stopifnot(!is.null(yx <- data$yx))
    mytrafo <- function(subset, weights, info, estfun = FALSE, object = FALSE, ...) {
      iy <- extree_variable(data, i = "yx", type = "index")
      if (is.null(iy)) {
        NAyx <- extree_variable(data, i ="yx", type = "missings")
        y <- yx$y
        x <- yx$x
        offset <- attr(yx$x, "offset")
        ### <FIXME> other ways of handling NAs necessary? </FIXME>
        subset <- subset[!(subset %in% NAyx)]
        if (NCOL(y) > 1) {
          y <- y[subset,,drop = FALSE]
        } else {
          y <- y[subset]
        }
        if (!is.null(x)) {
          ax <- attributes(x)
          ax$dim <- NULL
          ax$dimnames <- NULL
          x <- x[subset,,drop = FALSE]
          for (a in names(ax)) attr(x, a) <- ax[[a]] ### terms, formula, ... for predict
        }
        w <- weights[subset]
        offset <- offset[subset]
        cluster <- data[["(cluster)"]][subset]
        if (all(c("estfun", "object") %in% nf)) { 
          m <- trafo(y = y, x = x, offset = offset, weights = w, start = info$coef, 
                     cluster = cluster, estfun = estfun, object = object, ...)
        } else {
          obj <- trafo(y = y, x = x, offset = offset, weights = w, start = info$coef, 
                       cluster = cluster, ...)
          m <- list(coefficients = coef(obj),
                    objfun = -as.numeric(logLik(obj)),
                    estfun = NULL, object = NULL)
          if (estfun) m$estfun <- sandwich::estfun(obj)
          if (object) m$object <- obj
        }
        if (!is.null(ef <- m$estfun)) {
          ### ctree expects unweighted scores
          if (!isTRUE(m$unweighted) && is.null(varselect) && ctrl$unweighted) 
            m$estfun <- m$estfun / w
          Y <- matrix(0, nrow = nrow(model.frame(data)), ncol = ncol(ef))
          Y[subset,] <- m$estfun
          m$estfun <- Y
        }
      } else {
        w <- libcoin::ctabs(ix = iy, subset = subset, weights = weights)[-1]
        offset <- attr(yx$x, "offset")
        cluster <- model.frame(data, yxonly = TRUE)[["(cluster)"]]
        if (all(c("estfun", "object") %in% nf)) { 
          m <- trafo(y = yx$y, x = yx$x, offset = offset, weights = w, start = info$coef, 
                     cluster = cluster,
                     estfun = estfun, object = object, ...)
        } else {
          obj <- trafo(y = yx$y, x = yx$x, offset = offset, weights = w, start = info$coef, 
                       cluster = cluster, ...)
          m <- list(coefficients = coef(obj),
                    objfun = -as.numeric(logLik(obj)),
                    estfun = NULL, object = NULL)
          if (estfun) m$estfun <- sandwich::estfun(obj)
          if (object) m$object <- obj
          if (!is.null(obj$unweighted)) 
            m$unweighted <- obj$unweighted
          m$converged <- obj$converged ### may or may not exist
        }
        ### <FIXME> unweight scores in ctree or weight scores in
        ### mfluc (means: for each variable again) </FIXME>
        ### ctree expects unweighted scores
        if (!is.null(m$estfun))  {
          if (!isTRUE(m$unweighted) && is.null(varselect) && ctrl$unweighted) 
            m$estfun <- m$estfun / w
        }
        if (!is.null(ef <- m$estfun))
          m$estfun <- rbind(0, ef)
      }
      return(m)
    }
  }
  
  if (!ctrl$update) {
    rootestfun <- mytrafo(subset = subset, weights = weights)
    updatetrafo <- function(subset, weights, info, ...)
      return(rootestfun)
  } else {
    updatetrafo <- function(subset, weights, info, ...) {
      ret <- mytrafo(subset = subset, weights = weights, info = info, ...)
      if (is.null(ret$converged)) ret$converged <- TRUE
      conv <- TRUE
      if (is.function(converged)) conv <- converged(subset, weights)
      ret$converged <- ret$converged && conv
      if (!ret$converged) return(NULL)
      ret
    }
  }
  
  nm <- c("model", "trafo", "data", "subset", "weights", "whichvar", "ctrl")
  stopifnot(all(nm == names(formals(varselect))))
  stopifnot(all(nm == names(formals(splitselect))))
  stopifnot(all(nm == names(formals(svarselect))))
  stopifnot(all(nm == names(formals(ssplitselect))))
  
  if (!doFit) return(mytrafo)
  
  list(nodes = .extree_node(id = 1, data = data, trafo = updatetrafo, varselect = varselect, 
                            splitselect = splitselect, svarselect = svarselect, ssplitselect = ssplitselect, 
                            partyvars = partyvars, weights = weights, subset = subset, ctrl = ctrl),
       trafo = mytrafo)
}



### unbiased recursive partitioning: set up new node
.extree_node <- function
(
    id = 1L, 			### id of this node
    data, 			### full data, readonly
    trafo,
    varselect, 			### variable selection
    splitselect,                ### split selection
    svarselect,                 ### same for surrogate splits
    ssplitselect,               ### same for surrogate splits
    partyvars, 			### partytioning variables
                                ### a subset of 1:ncol(model.frame(data))
    weights = integer(0L),	### optional case weights
    subset, 			### subset of 1:nrow(data)
                                ### for identifying obs for this node
    ctrl, 			### extree_control()
    info = NULL,
    extree_env = NULL           ### environment for depth and maxid
) {
  ## initialize tracker: is this a terminal node?
  terminal <- FALSE

  ## is there any info to be preserved?

  ## depth keeps track of the depth of the tree: has to be < than maxdepth
  ## maxit is the largest id in the left subtree
  if (is.null(extree_env)) {
    extree_env <- new.env()
    assign("depth", 0L, envir = extree_env)
  }
  depth <- get("depth", envir = extree_env)
  assign("maxid", id, envir = extree_env)

  ## FIXME: (Z) initialize extree_env$criterion with NULL

  if (depth >= ctrl$maxdepth) terminal <- TRUE
  ## FIXME: (Z) previously the code stopped here with:
  ##   return(partynode(as.integer(id)))
  ## However, we should collect "info" if desired!
  ## Analogously for all subsequent checks for "terminal".

  ## check for stumps
  if (!terminal && id > 1L && ctrl$stump) terminal <- TRUE

  ## sw is basically the number of observations
  ## which has to be > minsplit in order to consider
  ## the node for splitting
  ## FIXME: (Z) currently computed even if we know terminal=TRUE already
  sw <- if (length(weights) > 0L) {
    if (ctrl$caseweights) sum(weights[subset]) else sum(weights[subset] > 0L)
  } else {
    length(subset)
  }
  if (!terminal && sw < ctrl$minsplit) terminal <- TRUE

  ## FIXME: (Z) split variable selection structure should be set up
  ## with NAs if terminal=TRUE but "criterion" selected for "info"
  svars <- which(partyvars > 0)
  if (ctrl$mtry < Inf) {
    mtry <- min(length(svars), ctrl$mtry)
    svars <- .resample(svars, mtry, prob = partyvars[svars])
  } 

  ## FIXME: (Z) include a condition like
  ##  if(!terminal || "trafo" %in% ctrl$save)
  thismodel <- trafo(subset = subset, weights = weights, info = info,
    estfun = TRUE, object = TRUE)
  if (is.null(thismodel)) terminal <- TRUE

  ## FIXME: (Z) canonicalize saveinfo to list(inner, terminal) with character vectors

  ## update sample size constraints on possible splits
  ## need to do this here because varselect might consider splits
  mb <- ctrl$minbucket
  mp <- ctrl$minprob
  swp <- ceiling(sw * mp)
  if (mb < swp) mb <- as.integer(swp)
  thisctrl <- ctrl
  thisctrl$minbucket <- mb
  ## FIXME: (Z) add canonicalized saveinfo argument

  ## split variable selection:
  ## - either already returns a finished "partysplit" object
  ## - or "criterion" matrix (typically statistics & p-value) for all variables
  ## FIXME: (Z) need to get "empty" criterion matrix with all NAs
  ## -> store in extree_env after computing varsel for the first time in the root node
  ## except when root node is already terminal
  varsel <- varselect(model = thismodel, trafo = trafo, data = data,
    subset = subset, weights = weights, whichvar = svars, ctrl = thisctrl)

  if (terminal || inherits(varsel, "partysplit")) {
    thissplit <- if(terminal) NULL else varsel
    info <- nodeinfo <- thismodel[!(names(thismodel) %in% c("estfun"))] ## FIXME: (Z) needs updating
    info$nobs <- sw
    if (!ctrl$saveinfo) info <- NULL
  } else {
    ## criterion matrix and some sanity checking
    p <- varsel$criterion
    if(!is.matrix(p) || NCOL(p) != length(svars)) {
      stop("variable selection function does not provide a valid criterion matrix")
    }

    ## adjust p-values (for non-NA p-values), if any
    if("p.value" %in% rownames(p)) p["p.value", ] <- ctrl$padjust(p["p.value", ])
    if("log.p.value" %in% rownames(p)) p["log.p.value", ] <- ctrl$padjust(p["log.p.value", ], log = TRUE)

    ## determine criterion
    criterion <- ctrl$criterion
    critvalue <- ctrl$critvalue
    ## use log p-values if available (and p-values wanted)
    if(criterion == "p.value" && "log.p.value" %in% rownames(p)) {
      criterion <- "log.p.value"
      critvalue <- log(critvalue)
    }
    ## minimize p-values? (or maximize other criterion)
    minp <- criterion %in% c("p.value", "log.p.value")
    ## extract criterion values
    crit <- p[criterion, , drop = TRUE]

    ## no admissible splits in any variable, e.g., all
    ## split variables constant
    if (all(is.na(crit)) || (minp && all(crit >= critvalue)) || (!minp && all(crit <= critvalue))) {
      terminal <- TRUE
      popt <- NA_real_
    }

    ## placeholder for missing criterion values
    crit[is.na(crit)] <- if(minp) Inf else -Inf

    if(!terminal) {
      ## optimal criterion, trying to break ties if any
      critopt <- if(minp) min(crit, na.rm = TRUE) else max(crit, na.rm = TRUE)
      ties <- which(abs(crit - critopt) < sqrt(.Machine$double.xmin))
      if(length(ties) > 1L) {
        ## break ties in p-values (if this is the "criterion") by using statistics (if available)
        if(minp && any(c("statistic", "log.statistic") %in% rownames(p))) {
          stat <- if("log.statistic" %in% rownames(p)) {
	    p["log.statistic", ties, drop = TRUE]
          } else {
	    p["statistic", ties, drop = TRUE]
          }
          ## subtract small value from (log-)p-value based on rank of (log-)statistic
          crit[ties] <- crit[ties] - rank(stat)/length(ties) * 1e-8      
        } else {
          ## FIXME: (Z) Can we do anything here? Always choose the first? Randomize?
        }
      }

      ## store up to "splittry" best variables
      jsel <- order(crit, decreasing = !minp)
      splittry <- pmin(sum(is.finite(crit)), ctrl$splittry)
      jsel <- jsel[1L:splittry]      
      jsel <- if(minp) jsel[crit[jsel] < critvalue] else jsel[crit[jsel] > critvalue]
    } else {
      jsel <- integer()
    }

    if(!terminal) { ## FIXME: (Z) if(!terminal || "criterion" %in% save)
      ## always include the criterion used as "criterion"
      p <- rbind(p, criterion = crit)

      ## switch to "statistic" and "p.value" omitting logs
      if(("log.statistic" %in% rownames(p)) && !("statistic" %in% rownames(p))) {
        p <- rbind(p, statistic = exp(p["statistic", ]))
      }
      if(("log.p.value" %in% rownames(p)) && !("p.value" %in% rownames(p))) {
        p <- rbind(p, p.value = exp(p["log.p.value", ]))
      }
      if (any(rownames(p) %in% c("log.statistic", "log.p.value"))) {
        p <- p[-which(rownames(p) %in% c("log.statistic", "log.p.value")), , drop = FALSE]  
      }
    
      ## store optimal p-value (FIXME: (Z) Always p-value - or optimal criterion?)
      jopt <- if(minp) which.min(crit) else which.max(crit)
      if("p.value" %in% rownames(p)) {
	      popt <- p["p.value", jopt]
      } else {
        popt <- NA_real_
      }
      names(popt) <- colnames(p)[jopt]
    }

    ## FIXME: (Z) Would we want this? I guess it's clearer/easier to keep these in.
    ## ## report on tests actually performed only
    ## p <- p[,!is.na(p["statistic",]) & is.finite(p["statistic",]),
    ##      drop = FALSE]
    info <- nodeinfo <- c(list(criterion = p, p.value = popt),  ## FIXME: (Z) update according to "save"
                varsel[!(names(varsel) %in% c("criterion", "converged"))],
                thismodel[!(names(thismodel) %in% c("estfun"))])
    info$nobs <- sw
    if (!ctrl$saveinfo) info <- NULL

    if (!is.null(varsel$splits)) {
      ### varselect may return of a list of partysplit objects; use these for
      ### splitting; varselect is responsible for making sure lookahead is implemented
      thissplit <- varsel$splits[[jsel[1L]]]
    } else {
      ### try to find an admissible split in data[, jsel]
      thissplit <- splitselect(model = thismodel, trafo = trafo, data = data, subset = subset, 
                  weights = weights, whichvar = jsel, ctrl = thisctrl)
    }
  }

  ### failed split search:
  if (is.null(thissplit))
    return(partynode(as.integer(id), info = info))

  ### successful split search: set-up node
  ret <- partynode(as.integer(id))
  ret$split <- thissplit
  ret$info <- info

  ### determine observations for splitting (only non-missings)
  snotNA <- subset[!subset %in% extree_variable(data, i = varid_split(thissplit), type = "missings")]
  if (length(snotNA) == 0)
    return(partynode(as.integer(id), info = info))
  ### and split observations
  kidids <- kidids_node(ret, model.frame(data), obs = snotNA)

  ### compute probability of going left / right
  prob <- tabulate(kidids) / length(kidids) 
  # names(dimnames(prob)) <- NULL
  if (ctrl$majority)  ### go with majority
    prob <- as.double((1L:length(prob)) %in% which.max(prob))
  if (is.null(ret$split$prob))
    ret$split$prob <- prob

  ### compute surrogate splits
  if (ctrl$maxsurrogate > 0L) {
    pv <- partyvars
    pv[varid_split(thissplit)] <- 0
    pv <- which(pv > 0)
    if (ctrl$numsurrogate)
      pv <- pv[sapply(model.frame(data)[, pv], function(x) is.numeric(x) || is.ordered(x))]
    ret$surrogates <- .extree_surrogates(kidids, data = data, 
      weights = weights, subset = snotNA, 
      whichvar = pv,
      varselect = svarselect, splitselect = ssplitselect, ctrl = ctrl)
  }
  kidids <- kidids_node(ret, model.frame(data), obs = subset)

  ### proceed recursively
  kids <- vector(mode = "list", length = max(kidids)) 
  nextid <- id + 1L
  for (k in 1L:max(kidids)) {
    nextsubset <- subset[kidids == k]
    assign("depth", depth + 1L, envir = extree_env)
    kids[[k]] <- .extree_node(id = nextid, data = data, 
      trafo = trafo,
      varselect = varselect, splitselect = splitselect,
      svarselect = svarselect, ssplitselect = ssplitselect, 
      partyvars = partyvars, 
      weights = weights, subset = nextsubset, 
      ctrl = ctrl, info = nodeinfo, extree_env = extree_env) ## FIXME: (Z) update saveinfo in ctrl
    ### was: nextid <- max(nodeids(kids[[k]])) + 1L
    nextid <- get("maxid", envir = extree_env) + 1L
  }
  ret$kids <- kids

  return(ret)
}



### control arguments needed in this file
extree_control <- function(
  ## unbiased inference-based variable selection  
  criterion = "p.value",
  critvalue = 0.05,
  padjust = TRUE,

  ## tree construction
  minsplit = 20L,
  minbucket = 7L, 
  minprob = 0.01, 
  stump = FALSE,
  maxsurrogate = 0L, 
  numsurrogate = FALSE,
  mtry = Inf,
  maxdepth = Inf, 
  splittry = 2L,
  majority = FALSE, 
  caseweights = TRUE, 
  update = NULL,

  ## currently not used/implemented
  nmax = Inf,        ## better name/settings?
  lookahead = FALSE, ## try trafo() for daugther nodes before implementing the split
  multiway = FALSE,  ## passed to split_select?

  ## TODO
  save = NULL,         ## FIXME: (Z) -> save instead of saveinfo
  varselect = NULL,    ## FIXME: (Z) add default (ctree?)
  splitselect = NULL,  ## FIXME: (Z) add default
  svarselect = NULL,   ## FIXME: (Z) add default
  ssplitselect = NULL, ## FIXME: (Z) add default

  ## technical (FIXME: (Z) not yet used, pass to select functions?)
  applyfun = NULL, 
  cores = NULL,

  ## allow extensions
  ...,

  ## legacy
  bonferroni = NULL,
  saveinfo = TRUE,
  selectfun = NULL,
  splitfun = NULL,
  svselectfun = NULL,
  svsplitfun = NULL,
  unweighted = TRUE ## FIXME: (Z) get rid of this option? needed for: m$estfun <- m$estfun / w
) {

  ## p-value adjustment (legacy option: bonferroni)
  if(!is.null(bonferroni)) padjust <- bonferroni
  if(is.logical(padjust)) {
    padjust <- if(padjust) "sidak" else "none"
  }
  if(is.character(padjust)) {
    padjust <- switch(tolower(padjust),
      "bonferroni" = p_adjust_bonferroni,
      "sidak" = p_adjust_sidak,
      "none" = function(p, ...) p,
      stop(paste("unknown padjust method:", padjust))
    )
  }
  if(!is.function(padjust)) stop("padjust must be a logical, character, or function")

  ## optionally employ parallel computing facilities
  if(is.null(applyfun)) {
    if(is.null(cores)) cores <- 1L
    applyfun <- if(cores == 1L) {
      lapply
    } else if(.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(cores) 
      on.exit(parallel::stopCluster(cl))
      function(X, FUN, ...) parallel::parLapply(cl, X, FUN, ...)
    } else {
      function(X, FUN, ...) parallel::mclapply(X, FUN, ..., mc.cores = cores)    
    }
  }

  ### well, it is implemented but not correctly so (FIXME)
  if (multiway & maxsurrogate > 0L) stop("surrogate splits currently not implemented for multiway splits")
  
  ## variable selection preprocessing
  if(!is.null(selectfun)) varselect <- selectfun
  if(is.list(varselect) || is.character(varselect) || "j" %in% names(formals(varselect))) {
    var_sel <- .preprocess_select(varselect, select_type = "var")
    varselect <- function(model, trafo, data, subset, weights, whichvar, ctrl) {
      var_select_loop(model, trafo, data, subset, weights, whichvar, ctrl, var_select = var_sel)
    }
  }
  
  ## split selection preprocessing
  if(!is.null(splitfun)) splitselect <- splitfun
  if(is.list(splitselect) || is.character(splitselect) || "j" %in% names(formals(splitselect))) {
    split_sel <- .preprocess_select(splitselect, select_type = "split")
    splitselect <- function(model, trafo, data, subset, weights, whichvar, ctrl) {
      split_select_loop(model = model, trafo = trafo, data = data, subset = subset,
        weights = weights, whichvar = whichvar, control = ctrl, split_select = split_sel)
    }
  }

  list(
    criterion = criterion,
    critvalue = critvalue,
    padjust = padjust,
    
    minsplit = minsplit,
    minbucket = minbucket, 
    minprob = minprob,
    stump = stump,
    maxsurrogate = maxsurrogate, 
    numsurrogate = numsurrogate,
    mtry = mtry,
    maxdepth = maxdepth,
    splittry = splittry,
    majority = majority,
    caseweights = caseweights, 
    update = update,

    varselect = varselect,
    splitselect = splitselect,
    svarselect = svarselect,
    ssplitselect = ssplitselect,

    applyfun = applyfun,
    save = save,

    nmax = nmax,
    lookahead = lookahead,
    multiway = multiway,
    saveinfo = saveinfo,
    unweighted = unweighted,
    ...
  )
}


### data preprocessing -> initial subset
.start_subset <- function(data) {
    ret <- 1:NROW(model.frame(data))
    if (length(data$yxmissings) > 0)
        ret <- ret[!(ret %in% data$yxmissings)]
    ret
}



### unbiased recursive partitioning: surrogate splits
.extree_surrogates <- function
(
  split, 			### integer vector with primary kidids
  data, 			### full data, readonly
  weights,
  subset, 			### subset of 1:nrow(data) with
  ### non-missings in primary split
  whichvar, 			### partytioning variables
  varselect, 			### variable selection and split
  ### function
  splitselect,
  ctrl			### ctree_control()
) {
  
  if (length(whichvar) == 0) return(NULL)
  ms <- max(split)
  if (ms != 2) return(NULL) ### ie no multiway splits!
  dm <- matrix(0, nrow = nrow(model.frame(data)), ncol = ms)
  dm[cbind(subset, split)] <- 1
  thismodel <- list(estfun = dm)
  varsel <- varselect(model = thismodel, trafo = NULL, data = data, subset = subset, 
                      weights = weights, whichvar = whichvar, ctrl = ctrl)
  p <- varsel$criterion
  ### partykit always used p-values, so expect some differences
  crit <- p[ctrl$criterion,,drop = TRUE]
  ### crit is maximised, but there might be ties
  ties <- which(abs(crit - max(crit, na.rm = TRUE)) < .Machine$double.eps)
  if (length(ties) > 1) {
    ### add a small value (< 1/1000) to crit derived from order of 
    ### teststat
    crit[ties] <- crit[ties] + 
      order(p["statistic", ties]) / (sum(ties) * 1000)
  }
  
  ret <- vector(mode = "list", length = min(c(length(whichvar), 
                                              ctrl$maxsurrogate)))
  
  for (i in 1L:length(ret)) {
    jsel <- which.max(crit)
    thisctrl <- ctrl
    thisctrl$minbucket <- 0L
    sp <- splitselect(model = thismodel, trafo = NULL, data = data, subset = subset, 
                      weights = weights, whichvar = jsel, ctrl = ctrl)
    if (is.null(sp)) next
    ret[[i]] <- sp
    tmp <- kidids_split(ret[[i]], model.frame(data), obs = subset)
    
    ### <FIXME> this needs fixing for multiway "split"
    tab <- table(tmp, split)
    if (tab[1, 1] < tab[1, 2]) {
      indx <- ret[[i]]$index
      ret[[i]]$index[indx == 1] <- 2L
      ret[[i]]$index[indx == 2] <- 1L
    }
    ### </FIXME>
    crit[which.max(crit)] <- -Inf
  }
  ret <- ret[!sapply(ret, is.null)]
  if (length(ret) == 0L) ret <- NULL
  return(ret)
}
