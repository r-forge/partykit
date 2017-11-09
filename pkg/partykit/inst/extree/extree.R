
### unbiased recursive partitioning: set up new node
.extree_node <- function
(
    id = 1L, 			### id of this node
    data, 			### full data, readonly
    trafo,
    selectfun, 			### variable selection
    splitfun,                   ### split selection
    svselectfun = selectfun,    ### same for surrogate splits
    svsplitfun = splitfun,      ### same for surrogate splits
    partyvars, 			### partytioning variables
                                ### a subset of 1:ncol(model.frame(data))
    weights = integer(0L),	### optional case weights
    subset, 			### subset of 1:nrow(data)
                                ### for identifying obs for this node
    ctrl, 			### .urp_control()
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

    thismodel <- trafo(subset = subset, weights = weights, info = info)
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
                              sf[!(names(sf) %in% c("criteria"))],
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
    snotNA <- subset[!subset %in% .get_NAs(data, varid_split(thissplit))]
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

    dm <- matrix(0, nrow = nrow(model.frame(data)), ncol = max(split))
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

.extree_fit <- function(data, trafo, partyvars, subset, weights, ctrl) {
    ret <- list()

    nf <- names(formals(trafo))
    if (all(c("subset", "weights", "info") %in% nf)) {
        mytrafo <- trafo
    } else {
        stopifnot(all(c("y", "x", "offset", "weights", "start") %in% nf))
        stopifnot(!is.null(yx <- data$yx))
        mytrafo <- function(subset, weights, info, ...) {
            iy <- .get_index(data, "yx")
            if (is.null(iy)) {
                NAyx <- .get_missings(data, "yx")
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
                y <- y[-subset,,drop = FALSE]
		x <- x[-subset,,drop = FALSE]
                weights <- weights[-subset]
                offset <- offset[-subset]
                m <- trafo(y = y, x = x, offset = offset, weights = weights, start = info$coef, ...)
                if (!is.null(ef <- m$estfun)) {
                    Y <- matrix(0, nrow = nrow(model.frame(data)), ncol = ncol(ef))
                    Y[subset,] <- ef
                    m$estfun <- Y
                }
            } else {
                w <- inum::ctabs(ix = iy, subset = subset, weights = weights)[-1]
                offset <- attr(yx$x, "offset")
                m <- trafo(y = yx$y, x = yx$x, offset = offset, weights = w, start = info$coef, ...)
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
            ret <- trafo(subset = subset, weights = weights, info = info, ...)
            if (!ret$converged) return(NULL)
            ret
        }
    }

    selectfun <- function(model, subset, weights, whichvar, ctrl) {
        ret <- list(criteria = matrix(NA, nrow = 2L, ncol = ncol(model.frame(data))))
        rownames(ret$criteria) <- c("statistic", "p.value")
        colnames(ret$criteria) <- names(model.frame(data))
        ### <FIXME> allow joint MC in the absense of missings; fix seeds
        ### </FIXME>
        for (j in whichvar) {
            tst <- switch(ctrl$testflavour,
                "ctree" = .ctree_test(model = model, trafo = trafo, data = data, 
                                      subset = subset, weights = weights, j = j, 
                                      SPLITONLY = FALSE, ctrl = ctrl),
                "exhaustive" = .objfun_test(model = model, trafo = trafo, data = data, 
                                            subset = subset, weights = weights, j = j, 
                                            SPLITONLY = FALSE, ctrl = ctrl),
                "mfluc" = .mfluc_test(model = model, trafo = trafo, data = data, 
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
    
    splitfun <- function(model, subset, weights, whichvar, ctrl) {
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
                index[xt > 0 & xt < minbucket] <- nlevels(x) + 1L
                if (length(unique(index)) == 1) {
                    ret <- NULL
                } else {
                    index <- unclass(factor(index))
                    ret <- partysplit(as.integer(j), index = as.integer(index))
                }
            } else {
                ret <- switch(ctrl$splitflavour, 
                "ctree" = .ctree_test(model = model, trafo = trafo, data = data, 
                                      subset = subset, weights = weights, j = j, 
                                      SPLITONLY = TRUE, ctrl = ctrl),
                "exhaustive" = .objfun_test(model = model, trafo = trafo, data = data, 
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

    .extree_node(id = 1, data = data, trafo = updatetrafo, selectfun = selectfun, 
                 splitfun = splitfun, partyvars = partyvars, weights = weights, 
                 subset = subset, ctrl = ctrl)
}

## extensible tree (model) function
extree <- function(formula, data, subset, na.action = na.pass, weights, offset, cluster,
  yx = "none", nmax = c("yx" = Inf, "z" = Inf), ...)
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


  ## FIXME: separate object with options for: discretization, condensation, some NA handling
  ## below is just "proof-of-concept" implementation using plain model.matrix() which could
  ## be included as one option...
  if(yx == "matrix") {
  
    ## fake formula/terms if necessary
    if(noformula) {
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
    }
  }

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

  if (length(nmax) == 1) nmax <- c("yx" = nmax, "z" = nmax)
  ret$zindex <- inum::inum(mf, ignore = names(mf)[zerozvars], total = FALSE, 
                           nmax = nmax["z"])
  ret$yxindex <- NULL
  yxmf <- mf
  if (is.finite(nmax["yx"])) {
      ret$yxindex <- inum::inum(mf[, yxvars, drop = FALSE], total = TRUE, 
                                as.interval = names(mf)[vars$y], complete.cases.only = TRUE, 
                                nmax = nmax["yx"])
      yxmf <- attr(ret$yxindex, "levels")
      yxmf[["(weights)"]] <- NULL
      attr(ret$yxindex, "levels") <- yxmf
  }

  ret$missings <- lapply(ret$data, function(x) which(is.na(x)))
  ret$yxmissings <- sort(unique(do.call("c", ret$missings[yxvars])))

  if (yx == "matrix") {
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

