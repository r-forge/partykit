
### unbiased recursive partitioning: set up new node
.extree_node <- function
(
    id = 1L, 			### id of this node
    data, 			### full data, readonly
    trafo,
    selectfun, 			### variable selection
    splitfun,                   ### split selection
    svselectfun,                ### same for surrogate splits
    svsplitfun,                 ### same for surrogate splits
    partyvars, 			### partytioning variables
                                ### a subset of 1:ncol(model.frame(data))
    weights = integer(0L),	### optional case weights
    subset, 			### subset of 1:nrow(data)
                                ### for identifying obs for this node
    ctrl, 			### extree_control()
    info = NULL,
    cenv = NULL			### environment for depth and maxid
) {

    ### depth keeps track of the depth of the tree
    ### which has to be < than maxdepth
    ### maxit is the largest id in the left subtree
    if (is.null(cenv)) {
        cenv <- new.env()
        assign("depth", 0L, envir = cenv)
    }
    depth <- get("depth", envir = cenv)
    assign("maxid", id, envir = cenv)
    if (depth >= ctrl$maxdepth)
        return(partynode(as.integer(id)))

    ### check for stumps
    if (id > 1L && ctrl$stump) 
        return(partynode(as.integer(id)))

    ### sw is basically the number of observations
    ### which has to be > minsplit in order to consider
    ### the node for splitting
    if (length(weights) > 0L) {
        if (ctrl$caseweights) {
            sw <- sum(weights[subset]) 
        } else {
            sw <- sum(weights[subset] > 0L)
        }
    } else {
        sw <- length(subset)
    }
    if (sw < ctrl$minsplit) 
        return(partynode(as.integer(id)))

    svars <- which(partyvars > 0)
    if (ctrl$mtry < Inf) {
        mtry <- min(sum(partyvars > 0), ctrl$mtry)
        svars <- .resample(partyvars, mtry, prob = partyvars)
    } 

    thismodel <- trafo(subset = subset, weights = weights, info = info, 
                       estfun = TRUE, object = TRUE)
    if (is.null(thismodel))
        return(partynode(as.integer(id)))

    ### update sample size constraints on possible splits
    ### need to do this here because selectfun might consider splits
    mb <- ctrl$minbucket
    mp <- ctrl$minprob
    swp <- ceiling(sw * mp)
    if (mb < swp) mb <- as.integer(swp)
    thisctrl <- ctrl
    thisctrl$minbucket <- mb

    ### compute test statistics and p-values
    ### for _unbiased_ variable selection
    sf <- selectfun(thismodel, subset = subset, weights = weights, 
                    whichvar = svars, ctrl = thisctrl)

    if (inherits(sf, "partysplit")) {
        thissplit <- sf
        info <- nodeinfo <- thismodel[!(names(thismodel) %in% c("estfun"))]
        if (!ctrl$saveinfo) info <- NULL
    } else {
        ### selectfun might return other things later to be used for info
        p <- sf$criteria
        crit <- p[ctrl$criterion,,drop = TRUE]
        if (all(is.na(crit))) 
            return(partynode(as.integer(id)))

        crit[is.na(crit)] <- -Inf
        ### crit is maximised, but there might be ties
        ties <- which(abs(crit - max(crit, na.rm = TRUE)) < sqrt(.Machine$double.xmin))
        if (length(ties) > 1) {
            ### add a small value (< 1/1000) to crit derived from rank of 
            ### teststat
            crit[ties] <- crit[ties] + 
                rank(p["statistic", ties]) / (sum(ties) * 1000)
        }

        ### switch from log(1 - pval) to pval for info slots
        ### switch from log(statistic) to statistic
        ### criterion stays on log scale to replicate variable selection
        p <- rbind(p, criterion = crit)
        p["statistic",] <- exp(p["statistic",])
        p["p.value",] <- -expm1(p["p.value",])
        pmin <- p["p.value", which.max(crit)]
        names(pmin) <- colnames(model.frame(data))[which.max(crit)]

        ### report on tests actually performed only
        p <- p[,!is.na(p["statistic",]) & is.finite(p["statistic",]),
               drop = FALSE]
        info <- nodeinfo <- c(list(criterion = p, p.value = pmin), 
                              sf[!(names(sf) %in% c("criteria", "converged"))],
                              thismodel[!(names(thismodel) %in% c("estfun"))])
        if (!ctrl$saveinfo) info <- NULL

        ### nothing "significant"
        if (all(crit <= ctrl$logmincriterion))
            return(partynode(as.integer(id), info = info))

        ### at most ctrl$splittry variables with meaningful criterion
        st <- pmin(sum(is.finite(crit)), ctrl$splittry)
        jsel <- rev(order(crit))[1:st]
        jsel <- jsel[crit[jsel] > ctrl$logmincriterion]
        ### try to find an admissible split in data[, jsel]
        thissplit <- splitfun(thismodel, subset = subset, weights = weights, 
                              whichvar = jsel, ctrl = thisctrl)
    }

    ### failed split search:
    if (is.null(thissplit))
        return(partynode(as.integer(id), info = info))

    ### successful split search: set-up node
    ret <- partynode(as.integer(id))
    ret$split <- thissplit
    ret$info <- info

    ### determine observations for splitting (only non-missings)
    snotNA <- subset[!subset %in% data[[varid_split(thissplit), type = "missings"]]]
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
        if (ctrl$numsurrogate)
            pv <- pv[sapply(model.frame(data)[, pv], function(x) is.numeric(x) || is.ordered(x))]
        ret$surrogates <- .extree_surrogates(kidids, data = data, 
            weights = weights, subset = snotNA, 
            partyvars = pv,
            selectfun = svselectfun, splitfun = svsplitfun, ctrl = ctrl)
    }
    kidids <- kidids_node(ret, model.frame(data), obs = subset)

    ### proceed recursively
    kids <- vector(mode = "list", length = max(kidids)) 
    nextid <- id + 1L
    for (k in 1L:max(kidids)) {
        nextsubset <- subset[kidids == k]
        assign("depth", depth + 1L, envir = cenv)
        kids[[k]] <- .extree_node(id = nextid, data = data, 
            trafo = trafo,
            selectfun = selectfun, splitfun = splitfun,
            svselectfun = svselectfun, svsplitfun = svsplitfun, 
            partyvars = partyvars, 
            weights = weights, subset = nextsubset, 
            ctrl = ctrl, info = nodeinfo, cenv = cenv)
        ### was: nextid <- max(nodeids(kids[[k]])) + 1L
        nextid <- get("maxid", envir = cenv) + 1L
    }
    ret$kids <- kids

    return(ret)
}

### unbiased recursive partitioning: surrogate splits
.extree_surrogates <- function
(
    split, 			### integer vector with primary kidids
    data, 			### full data, readonly
    weights,
    subset, 			### subset of 1:nrow(data) with
				### non-missings in primary split
    partyvars, 			### partytioning variables
    selectfun, 			### variable selection and split
				### function
    splitfun,
    ctrl			### ctree_control()
) {

    ms <- max(split)
    if (ms == 2) return(NULL) ### ie no multiway splits!
    dm <- matrix(0, nrow = nrow(model.frame(data)), ncol = ms)
    dm[cbind(subset, split)] <- 1
    thismodel <- list(estfun = dm)
    sf <- selectfun(thismodel, subset = subset, weights = weights, whichvar = which(partyvars > 0),
    ctrl = ctrl)
    p <- sf$criteria
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

    ret <- vector(mode = "list", length = min(c(sum(partyvars > 0), 
                                                ctrl$maxsurrogate)))

    for (i in 1L:length(ret)) {
        jsel <- which.max(crit)
        thisctrl <- ctrl
        thisctrl$minbucket <- 0L
        sp <- splitfun(thismodel, subset = subset, weights = weights, whichvar = jsel,
        ctrl = ctrl)
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

extree_fit <- function(data, trafo, converged, selectfun = NULL, 
                       splitfun = NULL, svselectfun = NULL, 
                       svsplitfun = NULL, partyvars, subset, weights, ctrl) {
    ret <- list()

    nf <- names(formals(trafo))
    if (all(c("subset", "weights", "info", "estfun", "object") %in% nf)) {
        mytrafo <- trafo
    } else {
        stopifnot(all(c("y", "x", "offset", "weights", "start") %in% nf))
        stopifnot(!is.null(yx <- data$yx))
        mytrafo <- function(subset, weights, info, estfun = FALSE, object = FALSE, ...) {
            iy <- data[["yx", type = "index"]]
            if (is.null(iy)) {
                NAyx <- data[["yx", type = "missing"]]
                if (length(NAyx) > 0) {
                    y <- matrix(0, nrow = nrow(model.frame(data)), ncol = ncol(yx$y))
                    y[-NAyx, ,drop = FALSE] <- yx$y
                    x <- matrix(0, nrow = nrow(model.frame(data)), ncol = ncol(yx$x))
                    x[-NAyx, ,drop = FALSE] <- yx$x
                    offset <- numeric(nrow(model.frame(data)))
                    offset[-NAyx] <- attr(yx$x, "offset")
                } else {
                    y <- yx$y
                    x <- yx$x
                    offset <- attr(yx$x, "offset")
                }
                subset <- subset[!(subset %in% NAyx)]
                y <- y[subset,,drop = FALSE]
                if (!is.null(yx$x)) {
                    ax <- attributes(yx$x)
                    ax$dim <- NULL
                    ax$dimnames <- NULL
                    x <- x[subset,,drop = FALSE]
                    for (a in names(ax)) attr(x, a) <- ax[[a]] ### terms, formula, ... for predict
                }
                w <- weights[subset]
                offset <- offset[subset]
                if (all(c("estfun", "object") %in% nf)) { 
                    m <- trafo(y = y, x = x, offset = offset, weights = w, start = info$coef, 
                               estfun = estfun, object = object, ...)
                } else {
                    obj <- trafo(y = y, x = x, offset = offset, weights = w, start = info$coef, ...)
                    m <- list(coefficients = coef(obj),
                              objfun = logLik(obj),
                              estfun = NULL, object = NULL)
                    if (estfun) m$estfun <- sandwich::estfun(obj)
                    if (object) m$object <- obj
                }
                if (!is.null(ef <- m$estfun)) {
                    ### ctree expects unweighted scores
                    if (is.null(selectfun) && ctrl$testflavour == "ctree")
                        m$estfun <- m$estfun / w
                    Y <- matrix(0, nrow = nrow(model.frame(data)), ncol = ncol(ef))
                    Y[subset,] <- m$estfun
                    m$estfun <- Y
                }
            } else {
                w <- libcoin::ctabs(ix = iy, subset = subset, weights = weights)[-1]
                offset <- attr(yx$x, "offset")
                if (all(c("estfun", "object") %in% nf)) { 
                    m <- trafo(y = yx$y, x = yx$x, offset = offset, weights = w, start = info$coef, 
                               estfun = estfun, object = object, ...)
                } else {
                    obj <- trafo(y = yx$y, x = yx$x, offset = offset, weights = w, start = info$coef, ...)
                    m <- list(coefficients = coef(obj),
                              objfun = logLik(obj),
                              estfun = NULL, object = NULL)
                    if (estfun) m$estfun <- sandwich::estfun(obj)
                    if (object) m$object <- obj
                }
                m$converged <- obj$converged ### may or may not exist
                ### ctree expects unweighted scores
                if (is.null(selectfun) && ctrl$testflavour == "ctree") {
                    if (!is.null(m$estfun)) m$estfun <- m$estfun / w
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

    .selectfun <- function(model, subset, weights, whichvar, ctrl) {
        ret <- list(criteria = matrix(NA, nrow = 2L, ncol = ncol(model.frame(data))))
        rownames(ret$criteria) <- c("statistic", "p.value")
        colnames(ret$criteria) <- names(model.frame(data))
        ### <FIXME> allow joint MC in the absense of missings; fix seeds
        ### </FIXME>
        for (j in whichvar) {
            tst <- switch(ctrl$testflavour,
                "ctree" = .ctree_test(model = model, trafo = mytrafo, data = data, 
                                      subset = subset, weights = weights, j = j, 
                                      SPLITONLY = FALSE, ctrl = ctrl),
                "exhaustive" = .objfun_test(model = model, trafo = mytrafo, data = data, 
                                            subset = subset, weights = weights, j = j, 
                                            SPLITONLY = FALSE, ctrl = ctrl),
                "mfluc" = .mfluc_test(model = model, trafo = mytrafo, data = data, 
                                      subset = subset, weights = weights, j = j, 
                                      SPLITONLY = FALSE, ctrl = ctrl),
                stop("not implemented")
            )
            ret$criteria["statistic",j] <- tst$statistic
            ret$criteria["p.value",j] <- tst$p.value
        }
        if (ctrl$bonferroni) 
            ### make sure to correct for _non-constant_ variables only
            ret$criteria["p.value",] <- ret$criteria["p.value",] * 
                                        sum(!is.na(ret$criteria["p.value",]))
        ret
    }
    
    .splitfun <- function(model, subset, weights, whichvar, ctrl) {
        for (j in whichvar) {
            x <- model.frame(data)[[j]]
            if (ctrl$multiway && is.factor(x) && !is.ordered(x) &&
                (ctrl$maxsurrogate == 0) && nlevels(x[subset, drop = TRUE]) > 1) 
            {
                index <- 1L:nlevels(x)
                xt <- libcoin::ctabs(ix = unclass(x), weights = weights, subset = subset)[-1]
                index[xt == 0] <- NA
                ### maybe multiway is not so smart here as
                ### nodes with nobs < minbucket could result
                index[xt > 0 & xt < ctrl$minbucket] <- nlevels(x) + 1L
                if (length(unique(index)) == 1) {
                    ret <- NULL
                } else {
                    index <- unclass(factor(index))
                    ret <- partysplit(as.integer(j), index = as.integer(index))
                }
            } else {
                ret <- switch(ctrl$splitflavour, 
                "ctree" = .ctree_test(model = model, trafo = mytrafo, data = data, 
                                      subset = subset, weights = weights, j = j, 
                                      SPLITONLY = TRUE, ctrl = ctrl),
                "exhaustive" = .objfun_test(model = model, trafo = mytrafo, data = data, 
                                            subset = subset, weights = weights, j = j, 
                                            SPLITONLY = TRUE, ctrl = ctrl),
                stop("not implemented"))
            }
            ### check if trafo can be successfully applied to all daugther nodes 
            ### (converged = TRUE)
            if (ctrl$lookahead & !is.null(ret)) {
                sp <- kidids_split(ret, model.frame(data), obs = subset)
                conv <- sapply(unique(sp), function(i)
                        trafo(subset[sp == i], weights = weights)$converged)
                if (!all(conv)) ret <- NULL
            }
            if (!is.null(ret)) break()
        }
        return(ret)
    }

    .svselectfun <- function(..., ctrl) {
        ctrl$testflavour <- "ctree"
        ctrl$splitflavour <- "ctree"
        .selectfun(..., ctrl = ctrl)
    }

    .svsplitfun <- function(..., ctrl) {
        ctrl$testflavour <- "ctree"
        ctrl$splitflavour <- "ctree"
        .splitfun(..., ctrl = ctrl)
    }

    if (is.null(selectfun)) {
        selectfun <- .selectfun
    } else {
        stopifnot(all(c("model", "subset", "weights", "whichvar", "ctrl") %in% names(formals(selectfun))))
    }

    if (is.null(splitfun)) {
        splitfun <- .splitfun
    } else {
        stopifnot(all(c("model", "subset", "weights", "whichvar", "ctrl") %in% names(formals(splitfun))))
    }

    if (is.null(svselectfun)) {
        svselectfun <- .svselectfun
    } else {
        stopifnot(all(c("model", "subset", "weights", "whichvar", "ctrl") %in% names(formals(svselectfun))))
    }

    if (is.null(svsplitfun)) {
        svsplitfun <- .svsplitfun
    } else {
        stopifnot(all(c("model", "subset", "weights", "whichvar", "ctrl") %in% names(formals(svsplitfun))))
    }

    list(nodes = .extree_node(id = 1, data = data, trafo = updatetrafo, selectfun = selectfun, 
                              splitfun = splitfun, svselectfun = svselectfun, svsplitfun = svsplitfun, 
                              partyvars = partyvars, weights = weights, subset = subset, ctrl = ctrl),
         trafo = mytrafo)
}

## extensible tree (model) function
extree_data <- function(formula, data, subset, na.action = na.pass, weights, offset, cluster,
  scores = NULL, yx = "none", nmax = c("yx" = Inf, "z" = Inf), ...)
{
  ## call
  cl <- match.call()

  ## 'formula' may either be a (multi-part) formula or a list
  noformula <- !inherits(formula, "formula")
  if(noformula) {

    ## formula needs to be a 'list' (if it is not a 'formula')
    if(!inherits(formula, "list")) stop("unsupported specification of 'formula'")    

    ## specified formula elements and overall call elements
    fonam <- names(formula)
    clnam <- names(cl)[-1L]
    vanam <- c("y", "x", "z", "weights", "offset", "cluster")
    
    ## y and z (and optionally x) need to be in formula
    if(!all(c("y", "z") %in% fonam)) stop("'formula' needs to specify at least a response 'y' and partitioning variables 'z'")
    if(!("x" %in% fonam)) formula$x <- NULL
    
    ## furthermore weights/offset/cluster may be in formula or call
    vars <- formula[vanam]
    names(vars) <- vanam
    if("weights" %in% clnam) {
      clvar <- try(weights, silent = TRUE)
      vars[["weights"]] <- c(vars[["weights"]], if(!inherits(clvar, "try-error")) clvar else deparse(cl$weights))
    }
    if("offset" %in% clnam) {
      clvar <- try(offset, silent = TRUE)
      vars[["offset"]] <- c(vars[["offset"]], if(!inherits(clvar, "try-error")) clvar else deparse(cl$offset))
    }
    if("cluster" %in% clnam) {
      clvar <- try(cluster, silent = TRUE)
      vars[["cluster"]] <- c(vars[["cluster"]], if(!inherits(clvar, "try-error")) clvar else deparse(cl$cluster))
    }
    
    ## sanity checking
    for(v in vanam) {
      if(!is.null(vars[[v]]) && !(is.numeric(vars[[v]]) | is.character(vars[[v]]) | is.logical(vars[[v]]))) {
        warning(sprintf("unknown specification of '%s', must be character, numeric, or logical", v))
        vars[v] <- list(NULL)
      }
    }
    if(!missing(subset)) warning("'subset' argument ignored in list specification of 'formula'")
    if(!missing(na.action)) warning("'na.action' argument ignored in list specification of 'formula'")    

    ## no terms (by default)
    mt <- NULL

  } else {

    ## set up model.frame() call
    mf <- match.call(expand.dots = FALSE)
    mf$na.action <- na.action ### evaluate na.action
    if(missing(data)) data <- environment(formula)
    m <- match(c("formula", "data", "subset", "na.action", "weights", "offset", "cluster"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf$dot <- "sequential"

    ## formula processing
    oformula <- as.formula(formula)
    formula <- Formula::as.Formula(formula)
    mf$formula <- formula
    npart <- length(formula)
    if(any(npart < 1L)) stop("'formula' must specify at least one left-hand and one right-hand side")
    npart <- npart[2L]

    ## evaluate model.frame
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())

    ## extract terms in various combinations
    mt <- list(
      "all" = terms(formula, data = data,                        dot = "sequential"),
      "y"   = terms(formula, data = data, rhs = 0L,              dot = "sequential"),
      "z"   = terms(formula, data = data, lhs = 0L, rhs = npart, dot = "sequential")
    )
    if(npart > 1L) {
      mt$yx <-terms(formula, data = data, rhs = 1L,              dot = "sequential")
      for(i in 1L:(npart-1L)) {
        mt[[paste("x", if(i == 1L) "" else i, sep = "")]] <- terms(
	            formula, data = data, lhs = 0L, rhs = i,     dot = "sequential")
      }
    }

    ## extract variable lists
    vars <- list(
      y = attr(mt$y, "term.labels"),
      x = unique(unlist(lapply(grep("^x", names(mt)), function(i) attr(mt[[i]], "term.labels")))),
      z = attr(mt$z, "term.labels"),
      weights = if("(weights)" %in% names(mf)) "(weights)" else NULL,
      offset  = if("(offset)"  %in% names(mf)) "(offset)"  else NULL,
      cluster = if("(cluster)" %in% names(mf)) "(cluster)" else NULL
    )
    ymult <- length(vars$y) >= 1L
    if(!ymult) vars$y <- names(mf)[1L]
    ## FIXME: store information which variable(s) went into (weights), (offset), (cluster)
    ## idea: check (x and) z vs. deparse(cl$weights), deparse(cl$offset), deparse(cl$cluster)

    ## check wether offset was inside the formula
    if(!is.null(off <- attr(mt$x, "offset"))) {
      if(is.null(vars$offset)) mf[["(offset)"]] <- rep.int(0, nrow(mf))
      for(i in off) mf[["(offset)"]] <- mf[["(offset)"]] + mf[[i]]
      vars$offset <- "(offset)"
    }
  }
  
  ## canonicalize y/x/z term labels
  vanam <- if(noformula) names(data) else names(mf)
  ## z to numeric
  if(is.null(vars$z)) stop("at least one 'z' variable must be specified")
  if(is.integer(vars$z)) vars$z <- vanam[vars$z]
  if(is.character(vars$z)) vars$z <- vanam %in% vars$z
  if(is.logical(vars$z)) vars$z <- as.numeric(vars$z)
  if(is.null(names(vars$z))) names(vars$z) <- vanam
  vars$z <- vars$z[vanam]
  if(any(is.na(vars$z))) vars$z[is.na(vars$z)] <- 0
  vars$z <- as.numeric(vars$z)
  ## all others to integer
  for(v in c("y", "x", "weights", "offset", "cluster")) {
    if(!is.null(vars[[v]])) {
      if(is.character(vars[[v]])) vars[[v]] <- match(vars[[v]], vanam)
      if(is.logical(vars[[v]])) vars[[v]] <- which(vars[[v]])
      if(any(is.na(vars[[v]]))) {
        vars[[v]] <- vars[[v]][!is.na(vars[[v]])]
        warning(sprintf("only found the '%s' variables: %s", v, paste(vanam[vars[[v]]], collapse = ", ")))
      }
    }
    vars[[v]] <- unique(as.integer(vars[[v]]))
  }
  if(is.null(vars$y)) stop("at least one 'y' variable must be specified")


  ## FIXME: subsequently fitting, testing, splitting
  ## - fit: either pre-processed _and_ subsetted data --or-- full data object plus subset vector
  ## - test: additionally needs fit output --and-- fit function
  ## - split: additionally needs test output
  ## - tbd: control of all details

  ret <- list(
    data = if(noformula) data else mf,
    variables = vars,
    terms = mt
  )

  mf <- ret$data
  yxvars <- c(vars$y, vars$x, vars$offset)
  zerozvars <- which(vars$z == 0)

  ret$scores <- vector(mode = "list", length = length(ret$variables$z))
  names(ret$scores) <- names(mf)
  if (!is.null(scores))
      ret$scores[names(scores)] <- scores

  if (length(nmax) == 1) nmax <- c("yx" = nmax, "z" = nmax)
  ret$zindex <- inum::inum(mf, ignore = names(mf)[zerozvars], total = FALSE, 
                           nmax = nmax["z"], meanlevels = FALSE)
  ret$yxindex <- NULL
  yxmf <- mf
  if (is.finite(nmax["yx"])) {
      ret$yxindex <- inum::inum(mf[, yxvars, drop = FALSE], total = TRUE, 
                                as.interval = names(mf)[vars$y], complete.cases.only = TRUE, 
                                nmax = nmax["yx"], meanlevels = FALSE)
      yxmf <- attr(ret$yxindex, "levels")
      yxmf[["(weights)"]] <- NULL
      attr(ret$yxindex, "levels") <- yxmf
  }

  ret$missings <- lapply(ret$data, function(x) which(is.na(x)))
  ret$yxmissings <- sort(unique(do.call("c", ret$missings[yxvars])))

  ## FIXME: separate object with options for: discretization, condensation, some NA handling
  ## below is just "proof-of-concept" implementation using plain model.matrix() which could
  ## be included as one option...
  if (yx == "matrix") {

    ## fake formula/terms if necessary
    formula <- Formula::as.Formula(sprintf("%s ~ %s | %s",
      paste(vanam[vars$y], collapse = " + "),
      if(length(vars$x) > 0L) paste(vanam[vars$x], collapse = " + ") else "0",
      paste(vanam[vars$z > 0], collapse = " + ")
    ))
    mt <- list(
      "all" = terms(formula),
      "y"   = terms(formula, data = data, rhs = 0L),
      "z"   = terms(formula, data = data, lhs = 0L, rhs = 2L),
      "yx"  = terms(formula, data = data, rhs = 1L),
      "x"   = terms(formula, data = data, lhs = 0L, rhs = 1L)
    )
    ymult <- length(vars$y) > 1L
    npart <- 2L

    yx <- list("y" = model.matrix(~ 0 + ., Formula::model.part(formula, yxmf, lhs = TRUE)))         
    for(i in (1L:npart)[-npart]) {
      ni <- paste("x", if(i == 1L) "" else i, sep = "")
      ti <- if(!ymult & npart == 2L) mt$yx else mt[[ni]]
      yx[[ni]] <- model.matrix(ti, yxmf)
      if(ncol(yx[[ni]]) < 1L) {
        yx[[ni]] <- NULL
      } else {
        attr(yx[[ni]], "formula") <- formula(formula, rhs = i)
        attr(yx[[ni]], "terms") <- ti
        attr(yx[[ni]], "offset") <- yxmf[["(offset)"]]
      }    
    }
    ret$yx <- yx
  }

  class(ret) <- "extree_data"
  ret
}

model.frame.extree_data <- function(object, yxonly = FALSE, ...) {
    if (!yxonly) 
        return(object$data)
    if (!is.null(object$yxindex))
        return(attr(object$yxindex, "levels"))
    vars <- object$variables
    return(object$data[, c(vars$y, vars$x, vars$offset),drop = FALSE])
}    

"[[.extree_data" <- function(x, i, type = c("original", "index", "scores", "missings")) {
    type <- match.arg(type)
    switch(type, 
        "original" = {
            mf <- model.frame(x)
            ### [[.data.frame needs lots of memory
            class(mf) <- "list"
            return(mf[[i]])
        },
        "index" = {
            if (i == "yx") return(x$yxindex)
            if (i %in% c(x$variables$y, x$variables$x))
                return(x$yxindex) ### may be NULL
            return(x$zindex[[i]])
        },
        "scores" = {
            f <- x[[i]]
            if (is.ordered(f)) {
                sc <- x$scores[[i]]
                if (is.null(sc)) sc <- 1:nlevels(f)
                return(sc)
            }
            return(NULL)
        },
        "missings" = {
            if (i == "yx" || i %in% c(x$variables$y, x$variables$x))
                return(x$yxmissings)
            x$missings[[i]]
        }
    )
}

model.weights.extree_data <- function(x)
    x[["(weights)"]]

model.offset.extree_data <- function(x)
    model.offset(model.frame(x, yxonly = TRUE))


### control arguments needed in this file
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
    update = TRUE
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


.objfun_test <- function(model, trafo, data, subset, weights, j, SPLITONLY, ctrl)
{

  x <- data[[j]]
  NAs <- data[[j, type = "missing"]]
  if (all(subset %in% NAs)) { 
    if (SPLITONLY) return(NULL)
    return(list(statistic = NA, p.value = NA))
  }

  ix <- data[[j, type = "index"]]
  ux <- attr(ix, "levels")
  ixtab <- libcoin::ctabs(ix = ix, weights = weights, subset = subset)[-1]
  ORDERED <- is.ordered(x) || is.numeric(x)
  
  linfo <- rinfo <- model
  maxlogLik <- nosplitll <- trafo(subset = subset, weights = weights, info = model, estfun = FALSE)$objfun
  sp <- NULL
  
  if (ORDERED) {
    ll <- ctrl$applyfun(which(ixtab > 0), function(u) {
      sleft <- subset[LEFT <- (ix[subset] <= u)]
      sright <- subset[!LEFT]
      if (length(weights) > 0) {
        if (sum(weights[sleft]) < ctrl$minbucket ||
            sum(weights[sright]) < ctrl$minbucket)
          return(-Inf);
      } else {
        if (length(sleft) < ctrl$minbucket || 
            length(sright) < ctrl$minbucket)
          return(-Inf);
      }
      if (ctrl$restart) {
        linfo <- NULL
        rinfo <- NULL
      }
      linfo <- trafo(subset = sleft, weights = weights, info = linfo, estfun = FALSE)
      rinfo <- trafo(subset = sright, weights = weights, info = rinfo, estfun = FALSE)
      ll <- linfo$objfun + rinfo$objfun
      return(ll)
    })
    maxlogLik <- max(unlist(ll))
    if(maxlogLik > nosplitll)
      sp <- which(ixtab > 0)[which.max(unlist(ll))]
    
  } else {
    xsubs <- factor(x[subset])
    ## stop if only one level left
    if(nlevels(xsubs) < 2) {
      if (SPLITONLY) {
        return(NULL)
      } else {
        return(list(statistic = NA, p.value = NA))
      } 
    }
    splits <- .mob_grow_getlevels(xsubs)
    ll <- ctrl$applyfun(1:nrow(splits), function(u) {
      sleft <- subset[LEFT <- xsubs %in% levels(xsubs)[splits[u,]]]
      sright <- subset[!LEFT]
      if (length(weights) > 0) {
        if (sum(weights[sleft]) < ctrl$minbucket ||
            sum(weights[sright]) < ctrl$minbucket)
          return(-Inf);
      } else {
        if (length(sleft) < ctrl$minbucket || 
            length(sright) < ctrl$minbucket)
          return(-Inf);
      }
      if (ctrl$restart) {
        linfo <- NULL
        rinfo <- NULL
      }
      linfo <- trafo(subset = sleft, weights = weights, info = linfo, estfun = FALSE)
      rinfo <- trafo(subset = sright, weights = weights, info = rinfo, estfun = FALSE)
      ll <- linfo$objfun + rinfo$objfun
      return(ll)
    })
    maxlogLik <- max(unlist(ll))
    if(maxlogLik > nosplitll) {
      sp <- splits[which.max(unlist(ll)),] + 1L
      levs <- levels(x)
      if(length(sp) != length(levs)) {
        sp <- sp[levs]
        names(sp) <- levs
      }
    }
  }
  
  if (!SPLITONLY){
    ## split only if logLik improves due to splitting
    maxlogLik <- ifelse(maxlogLik == nosplitll, NA, maxlogLik)
    return(list(statistic = maxlogLik, p.value = NA))
  }
  if (is.null(sp) || all(is.na(sp))) return(NULL)
  if (ORDERED) {
    ### interpolate split-points, see https://arxiv.org/abs/1611.04561
    if (!is.factor(x) & ctrl$intersplit & sp < length(ux)) {
      sp <- (ux[sp] + ux[sp + 1]) / 2 
    } else {
      sp <- ux[sp]  ### x <= sp vs. x > sp
    }
    if (is.factor(sp)) sp <- as.integer(sp)
    ret <- partysplit(as.integer(j), breaks = sp,
                      index = 1L:2L)
  } else {
    ret <- partysplit(as.integer(j),
                      index = as.integer(sp))
  }
  return(ret)
}
