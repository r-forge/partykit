nodeprune2.party <- function(x, ids, ...) {
  
  ### map names to nodeids
  if (!is.numeric(ids))
    ids <- match(ids, names(x))
  stopifnot(ids %in% nodeids(x))
  
  ### compute indices path to each node
  ### to be pruned off
  idxs <- lapply(ids, .get_path, obj = node_party(x))
  
  ### [[.party is NOT [[.list
  cls <- class(x)
  x <- unclass(x)
  ni <- which(names(x) == "node")
  
  for (i in 1:length(idxs)) {
    
    idx <- c(ni, idxs[[i]])
    ### check if we already pruned-off this node
    tmp <- try(x[[idx]], silent = TRUE)
    if (inherits(tmp, "try-error"))
      next()
    
    ### node ids of off-pruned daugther nodes
    idrm <- nodeids(x[[idx]])[-1]
    
    ### prune node by introducing a "new" terminal node
    x[[idx]] <- partynode(id = id_node(x[[idx]]),
                          info = info_node(x[[idx]]))
    
    ### constparty only: make sure the node ids in
    ### fitted are corrected
    if (length(idrm) > 0) {
      if(!is.null(x$fitted) && 
         "(fitted)" %in% names(x$fitted)) {
        j <- x$fitted[["(fitted)"]] %in% idrm
        x$fitted[["(fitted)"]][j] <- ids[i]
      }
    }
  }
  
  ### reindex to 1:max(nodeid)
  class(x) <- cls
  oldids <- nodeids(x)
  newids <- 1:length(nodeids(x))
  nodeids(x) <- newids
  
  for (i in seq_along(oldids)) {
    if (oldids[i] != newids[i]) {
      x$fitted[["(fitted)"]][x$fitted[["(fitted)"]] == oldids[i]] <- newids[i]
    }
  }
  
  return(x)
}



prune.modelparty <- function(object, type = "AIC")
{
  
  ## prepare pruning function
  if(is.character(type)) {
    type <- tolower(type)
    type <- match.arg(type, c("aic", "bic", "none"))
    
    if("lmtree" %in% class(object)) {
      type <- switch(type,
                     "aic" = {
                       function(objfun, df, nobs) (nobs[1L] * log(- objfun[1L]) + 2 * df[1L]) < (nobs[1L] * log(- objfun[2L]) + 2 * df[2L])
                     }, "bic" = {
                       function(objfun, df, nobs) (nobs[1L] * log(- objfun[1L]) + log(nobs[2L]) * df[1L]) < (nobs[1L] * log(- objfun[2L]) + log(nobs[2L]) * df[2L])
                     }, "none" = {
                       NULL
                     })
    } else {
      type <- switch(type,
                     "aic" = {
                       function(objfun, df, nobs) (2 * - objfun[1L] + 2 * df[1L]) < (2 * - objfun[2L] + 2 * df[2L])
                     }, "bic" = {
                       function(objfun, df, nobs) (2 * - objfun[1L] + log(n) * df[1L]) < (2 * - objfun[2L] + log(n) * df[2L])
                     }, "none" = {
                       NULL
                     })   
    }
  }
  if(!is.function(type)) {
    warning("Unknown specification of 'prune'")
    type <- NULL
  }
  
  ## degrees of freedom
  dfsplit <- object$info$control$dfsplit
  
  ## turn node to list
  node <- object$node
  nd <- as.list(node)
  
  ## if no pruning selected
  if(is.null(type)) return(nd)
  
  ## node information (IDs, kids, ...)
  id <- seq_along(nd)
  kids <- lapply(nd, "[[", "kids")
  tmnl <- sapply(kids, is.null)
  
  ## check nodes that only have terminal kids
  check <- sapply(id, function(i) !tmnl[i] && all(tmnl[kids[[i]]]))
  while(any(check)) {
    
    ## pruning node information
    pnode <- which(check)
    objfun <- sapply(nd, function(x) x$info$objfun)
    n <- nrow(object$fitted)
    pok <- sapply(pnode, function(i) type(
      objfun = c(objfun[i], sum(objfun[kids[[i]]])),
      df = c(length(nd[[1]]$info$coefficients), length(kids[[i]]) * length(nd[[1]]$info$coefficients) + as.integer(dfsplit)),
      nobs = c(nd[[i]]$info$nobs, n)
    ))
    
    ## do any nodes need pruning?
    pnode <- pnode[pok]
    if(length(pnode) < 1L) break
    
    ## prune
    object <- nodeprune2.party(object, ids = pnode)
    node <- object$node
    nd <- as.list(node)
    
    ## node information
    kids <- lapply(nd, "[[", "kids")
    tmnl <- sapply(kids, is.null)
    id <- seq_along(nd)
    check <- sapply(id, function(i) !tmnl[i] && all(tmnl[kids[[i]]]))
  }
  
  ## return pruned tree
  return(object)
}