
cforest <- function
(
    formula,
    data,   
    weights,
    subset, 
    offset, 
    cluster,
    na.action = na.pass,
    control = ctree_control(...),
    ytrafo = NULL, 
    scores = NULL, 
    ntree = 500L, 
    perturb = list(replace = FALSE, fraction = 0.632),
    mtry = ceiling(sqrt(nvar)), 
    applyfun = NULL,
    cores = NULL, 
    trace = FALSE,
    ...
) {
   
    ### get the call and the calling environment for .urp_tree
    call <- match.call(expand.dots = FALSE)
    call$na.action <- na.action
    frame <- parent.frame()
    if (missing(data)) {   
        data <- NULL
        data_asis <- FALSE
    } else {
        data_asis <- missing(weights) && missing(subset) &&
                     missing(cluster) && missing(offset)   
    }

    ### <FIXME> should be xtrafo
    if (!is.null(scores)) {
        if (missing(data)) 
            stop("can deal with scores with data being missing")
        for (n in names(scores)) {
            sc <- scores[[n]]
            if (is.ordered(data[[n]]) &&
                nlevels(data[[n]]) == length(sc)) {
                attr(data[[n]], "scores") <- as.numeric(sc)
            } else {
                warning("scores for variable ", sQuote(n), " ignored")
            }
        }
    }
    #### </FIXME>

    trafofun <- function(...) .ctreetrafo(..., ytrafo = ytrafo)
    tree <- .urp_tree(call, frame, data = data, data_asis = data_asis, control = control,
                      growfun = .ctreegrow, trafofun = trafofun,
                      doFit = FALSE)
    
    probw <- NULL
    weights <- model.weights(tree$mf)
    if (!is.null(weights)) {
        probw <- weights / sum(weights)
    } else {
        weights <- integer(0)
    }
    nvar <- length(tree$partyvars)
    control$mtry <- mtry
    applyfun <- control$applyfun

    idx <- 1L:nrow(tree$mf)
    size <- nrow(tree$mf)
    if (!perturb$replace) size <- floor(size * perturb$fraction)
    ### this is different now as we sample subsets and not weights
    rw <- replicate(ntree, sample(idx, size = size, replace = perturb$replace, prob = probw),
                    simplify = FALSE)
#    } else {
#        stopifnot(nrow(weights) == nrow(dat) && ncol(weights) == ntree)
#        rw <- as.data.frame(weights)
#        class(rw) <- "list"
#    }

        ## apply infrastructure for determining split points
    if (is.null(applyfun)) {
        applyfun <- if(is.null(cores)) {
            lapply  
        } else {
            function(X, FUN, ...)
                parallel::mclapply(X, FUN, ..., mc.cores = cores)
        }
    }
    ### no parallelization in ctree, only in cforest
    control$applyfun <- lapply

    if (trace) pb <- txtProgressBar(style = 3) 
    forest <- applyfun(1:ntree, function(b) {
        if (trace) setTxtProgressBar(pb, b/ntree)
        tree$treefun(tree$trafo, rw[[b]], integer(0))
    })
    if (trace) close(pb)

    fitted <- data.frame(idx = idx)  
    y <- model.part(Formula(formula), data = mf <- model.frame(Formula(formula), data = tree$mf),
                    lhs = 1, rhs = 0)
    if (length(y) == 1) y <- y[[1]]
    fitted[[2]] <- y
    names(fitted)[2] <- "(response)"
    fitted <- fitted[2]

    ### turn subsets in weights (maybe we can avoid this?)
    rw <- lapply(rw, function(x) tabulate(x, nbins = length(idx)))

    ret <- partykit:::constparties(nodes = forest, data = tree$mf, weights = rw,
                        fitted = fitted, terms = terms(mf), 
                        info = list(call = match.call(), control = control))
    class(ret) <- c("cforest", class(ret))

    return(ret)
}
