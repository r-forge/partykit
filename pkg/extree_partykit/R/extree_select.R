## helper function for .preprocess_varselect
.get_strategy_function <- function(strategy, select_type = "var") {
  
  
  ## find matching objects and set up list
  ## FIXME: (Z) better solution than via search()?
  snam <- sprintf(paste0("^", select_type, "select_%s"), strategy)
  onam <- unlist(lapply(search(), objects))
  onam <- unique(onam[grep(snam, onam)])
  strategy <- lapply(onam, get)
  
  ## supply variable type a name
  nam <- gsub(pattern = paste0(snam, "_"), replacement = "", onam)
  if (all(nam %in% onam)) nam <- regmatches(onam, regexpr("[a-z]+$", onam))
  names(strategy) <- nam
  
  
  ## drop non-functions
  strategy <- strategy[sapply(strategy, is.function)]
  
  return(strategy)
}




.preprocess_select <- function(select, select_type = "var") {
  
  ## function: return as is
  if (is.function(select)) return(select)
  
  ## character: return appropriate function
  if (is.character(select)) {
    return(
      .get_strategy_function(select, 
                             select_type = select_type)
    )
  }
  
  ## list: go through all elements and return accordingly
  if (is.list(select)) {
    
    get_strategy <- function(select_nam) {
      
      ## return function if function
      if (is.function(select[[select_nam]])) {
        sfun <- list(select[[select_nam]])
        names(sfun) <- select_nam
        return(sfun)
      } 
      
      ## get appropriat function if character
      if (is.character(select[[select_nam]])) {
        sfun <- .get_strategy_function(select[[select_nam]], 
                                       select_type = select_type)
        return(sfun)
      } 
      
      ## if none of the above -> ERROR
      stop("select can only be functions or characters.")
      
    }
    
    ## go through all list elements and choose approriate function
    select_list <- sapply(names(select), FUN = get_strategy, simplify = TRUE, USE.NAMES = FALSE)
    return(select_list)
  }
}



# <FIXME> (HS) better name of function
.get_varclass <- function(select_list, data, j) {
  ### which class is variable?
  varclass_orig <- class(extree_variable(x = data, i = j))
  
  ### Use most appropriate class (1st), if more than one is available
  ### Remove varclass if no varselect function is available
  varclass <- varclass_orig[varclass_orig %in% names(select_list)][1]
  
  ### if no function for this class is provided use default function
  if(length(varclass) == 0 | is.na(varclass)) {
    if (!("default" %in% names(select_list))) 
      stop("The is no specific or default select function for split variables of class ",
           varclass_orig, ". Please provide one.")
    varclass <- "default"
  } 
  
  return(varclass)
}

# <FIXME> (HS) Better name for function
selector <- function(select, model, trafo, data, subset, weights, whichvar, 
  control, j, split_only = FALSE) {
  
  
  # <FIXME> (HS) add check if function(s) return what we want
  
  if(is.function(select)) {
    ## if varselect is function, apply function
    ret <- select(model = model, trafo = trafo, data = data, 
                  subset = subset, weights = weights, j = j, 
                  split_only = split_only, control = control)
    
  } else if (is.list(select) && all(sapply(select, is.function))) {
    
    ## if varselect is list of functions, apply appropriate function
    varclass <- .get_varclass(select_list = select, data = data, 
                              j = j)
    
    ### Run appropriate varselect function
    ret <- select[[varclass]](model = model, trafo = trafo, data = data, 
                              subset = subset, weights = weights, j = j, 
                              split_only = split_only, control = control)
    
  } else {
    ## a future option would be a character which hints to functions 
    ## to be used.
    stop("Selection strategy must currently be a function or a named list of functions.")
  }
  
  return(ret)
}


## new selectfunction
varselect_loop <- function(model, trafo, data, subset, weights, whichvar, 
                            control, varselect) {
  
  ## if whichvar is empty return NULL
  if (length(whichvar) == 0L) return(NULL)
  
  ## set up return list + criterion matrix
  svarnams <- names(data$data)[whichvar]
  ret <- list(criterion = matrix(NA_real_, nrow = 2L, ncol = length(whichvar), 
    dimnames = list(c("statistic", "p.value"), svarnams)))
  ## TODO: (HS) it should be possible to have other things than statistic and p.value
  
  # znams <- names(data$data)[data$variables$z > 0] #attr(data$variables$z, "variable_names")
  # ret <- list(criterion = matrix(NA_real_, nrow = 2L, ncol = length(znams)))
  # ## FIXME: it is not ideal to use these names first and fix them later
  # rownames(ret$criterion) <- c("statistic", "p.value")
  # colnames(ret$criterion) <- znams
  
  ## loop over all relevant variables and use varselect function supplied
  for (i in seq_along(whichvar)) {
    
    tst <- selector(select = varselect, model = model, trafo = trafo, 
                    data = data, subset = subset, weights = weights, j = whichvar[i], 
                    control = control)
    
    logs <- "log.statistic" %in% names(tst)
    logp <- "log.p.value" %in% names(tst)
    
    if (is.null(tst)) {
      ret$criterion["statistic", i] <- NA
      ret$criterion["p.value", i] <- NA
    } else {
      ret$criterion["statistic", i] <- if(logs) tst$log.statistic else tst$statistic
      ret$criterion["p.value", i] <- if(logp) tst$log.p.value else tst$p.value
    }
  }
  if(logs) rownames(ret$criterion)[1L] <- "log.statistic"
  if(logp) rownames(ret$criterion)[2L] <- "log.p.value"
  ret
}

## Split function new
splitselect_loop <- function(model, trafo, data, subset, weights, whichvar, 
                              control, splitselect) {
  
  if (length(whichvar) == 0) return(NULL)
  
  ## loop over all vars in whichvar (stop if split is found)
  for (j in whichvar) {
    
    ret <- selector(select = splitselect, model = model, trafo = trafo, 
                    data = data, subset = subset, weights = weights, j = j, 
                    control = control, split_only = TRUE)
    
    ### check if trafo can be successfully applied to all daugther nodes 
    ### (converged = TRUE)
    if (control$lookahead & !is.null(ret)) {
      sp <- kidids_split(ret, model.frame(data), obs = subset)
      conv <- sapply(unique(na.omit(sp)), function(i)
        isTRUE(trafo(subset[sp == i & !is.na(sp)], weights = weights)$converged))
      if (!all(conv)) ret <- NULL
      ## FIXME: (Z) allow option to keep estfun and move on --> update = FALSE
    }
    
    ## stop if a split was found, otherwise continue with next possible var
    if (!is.null(ret)) break()
  }
  return(ret)
}

## Select function old
.select <- function(model, trafo, data, subset, weights, whichvar, ctrl, FUN) {
  znams <- names(data$data)[data$variables$z > 0] 
  ret <- list(criterion = matrix(NA, nrow = 2L, ncol = length(znams)))
  rownames(ret$criterion) <- c("statistic", "p.value")
  colnames(ret$criterion) <- znams
  if (length(whichvar) == 0) return(ret)
  ### <FIXME> allow joint MC in the absense of missings; fix seeds
  ### write ctree_test / ... with whichvar and loop over variables there
  ### </FIXME>
  for (i in seq_along(whichvar)) {
    tst <- FUN(model = model, trafo = trafo, data = data, 
               subset = subset, weights = weights, j = whichvar[i], 
               SPLITONLY = FALSE, ctrl = ctrl)
    ret$criterion["statistic", i] <- tst$statistic
    ret$criterion["p.value", i] <- tst$p.value
  }
  ret
}

## Split function old
.split <- function(model, trafo, data, subset, weights, whichvar, ctrl, FUN) {
  if (length(whichvar) == 0) return(NULL)
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
      ret <- FUN(model = model, trafo = trafo, data = data, 
                 subset = subset, weights = weights, j = j, 
                 SPLITONLY = TRUE, ctrl = ctrl)
    }
    ### check if trafo can be successfully applied to all daugther nodes 
    ### (converged = TRUE)
    if (ctrl$lookahead & !is.null(ret)) {
      sp <- kidids_split(ret, model.frame(data), obs = subset)
      conv <- sapply(unique(na.omit(sp)), function(i)
        isTRUE(trafo(subset[sp == i & !is.na(sp)], weights = weights)$converged))
      if (!all(conv)) ret <- NULL
    }
    if (!is.null(ret)) break()
  }
  return(ret)
}

.objfun_select <- function(...)
  function(model, trafo, data, subset, weights, whichvar, ctrl) {
    args <- list(...)
    ctrl[names(args)] <- args
    .select(model, trafo, data, subset, weights, whichvar, ctrl, FUN = .objfun_test)
  }

.objfun_split <- function(...)
  function(model, trafo, data, subset, weights, whichvar, ctrl) {
    args <- list(...)
    ctrl[names(args)] <- args
    .split(model, trafo, data, subset, weights, whichvar, ctrl, FUN = .objfun_test) 
  }

.objfun_test <- function(model, trafo, data, subset, weights, j, SPLITONLY, ctrl)
{
  
  x <- extree_variable(data, i = j)
  NAs <- extree_variable(data, i= j, type = "missing")
  if (all(subset %in% NAs)) { 
    if (SPLITONLY) return(NULL)
    return(list(statistic = NA, p.value = NA))
  }
  
  ix <- extree_variable(data, i = j, type = "index")
  ux <- attr(ix, "levels")
  ixtab <- libcoin::ctabs(ix = ix, weights = weights, subset = subset)[-1]
  ORDERED <- is.ordered(x) || is.numeric(x)
  
  linfo <- rinfo <- model
  minlogLik <- nosplitll <- trafo(subset = subset, weights = weights, info = model, estfun = FALSE)$objfun
  sp <- NULL
  
  if (ORDERED) {
    ll <- ctrl$applyfun(which(ixtab > 0), function(u) {
      sleft <- subset[LEFT <- (ix[subset] <= u)]
      sright <- subset[!LEFT]
      if (length(weights) > 0 && ctrl$caseweights) {
        if (sum(weights[sleft]) < ctrl$minbucket ||
            sum(weights[sright]) < ctrl$minbucket)
          return(Inf);
      } else {
        if (length(sleft) < ctrl$minbucket || 
            length(sright) < ctrl$minbucket)
          return(Inf);
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
    minlogLik <- min(unlist(ll))
    if(minlogLik < nosplitll)
      sp <- which(ixtab > 0)[which.min(unlist(ll))]
    
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
      if (length(weights) > 0 && ctrl$caseweights) {
        if (sum(weights[sleft]) < ctrl$minbucket ||
            sum(weights[sright]) < ctrl$minbucket)
          return(Inf);
      } else {
        if (length(sleft) < ctrl$minbucket || 
            length(sright) < ctrl$minbucket)
          return(Inf);
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
    minlogLik <- min(unlist(ll))
    if(minlogLik < nosplitll) {
      sp <- splits[which.min(unlist(ll)),] + 1L
      levs <- levels(x)
      if(length(sp) != length(levs)) {
        sp <- sp[levs]
        names(sp) <- levs
      }
    }
  }
  
  if (!SPLITONLY){
    ## split only if logLik improves due to splitting
    minlogLik <- ifelse(minlogLik == nosplitll, NA, minlogLik)
    return(list(statistic = -minlogLik, p.value = NA)) ### .extree_node maximises
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
