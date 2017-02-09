
### constructor for forest objects
constparties <- function(nodes, data, weights, fitted = NULL, terms = NULL, info = NULL) {

    stopifnot(all(sapply(nodes, function(x) inherits(x, "partynode"))))
    stopifnot(inherits(data, "data.frame"))
    stopifnot(inherits(weights, "list"))

    if(!is.null(fitted)) {
        stopifnot(inherits(fitted, "data.frame"))
        stopifnot(nrow(data) == 0L | nrow(data) == nrow(fitted))
        if (nrow(data) == 0L)
            stopifnot("(response)" %in% names(fitted))
    } else {
        stopifnot(nrow(data) > 0L)
        stopifnot(!is.null(terms))
        fitted <- data.frame("(response)" = model.response(model.frame(terms, data = data, 
                                                                       na.action = na.pass)),
                             check.names = FALSE)
    }

    ret <- list(nodes = nodes, data = data, weights = weights, fitted = fitted)
    class(ret) <- c("constparties", "parties")

    if(!is.null(terms)) {
        stopifnot(inherits(terms, "terms"))
        ret$terms <- terms
    }

    if (!is.null(info))
        ret$info <- info

    ret
}

.perturb <- function(replace = FALSE, fraction = .632) {
    ret <- function(prob) {
        if (replace) {
            rw <- rmultinom(1, size = length(prob), prob = prob)
        } else {
            rw <- integer(length(prob))
            i <- sample(1:length(prob), ceiling(fraction * length(prob)), prob = prob)
            rw[i] <- 1L
        }
        as.integer(rw)
    }
    ret
}

cforest <- function
(
    formula,
    data,   
    weights,
    subset, 
    offset, 
    cluster,
    strata,
    na.action = na.pass,
    control = ctree_control(mincriterion = 0, ...),
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

    call$weights <- NULL ### NOTE: trees are unweighted, weights enter sampling!
    trafofun <- function(...) .ctreetrafo(..., ytrafo = ytrafo)
    tree <- .urp_tree(call, frame, data = data, data_asis = data_asis, control = control,
                      trafofun = trafofun, doFit = FALSE)

    strata <- tree$mf[["(strata)"]]
    if (!is.null(strata)) {
        if (!is.factor(strata)) stop("strata is not a single factor")
    }
    
    probw <- NULL
    weights <- model.weights(tree$mf)
    if (!is.null(weights)) {
        probw <- weights / sum(weights)
    } else {
        weights <- integer(0)
    }
    nvar <- length(tree$partyvars)
    control$mtry <- mtry
    control$applyfun <- NULL

    idx <- 1L:nrow(tree$mf)
    if (is.null(strata)) {
        size <- nrow(tree$mf)
        if (!perturb$replace) size <- floor(size * perturb$fraction)
        rw <- replicate(ntree, sample(idx, size = size, replace = perturb$replace, prob = probw),
                        simplify = FALSE)
    } else {
        frac <- if (!perturb$replace) perturb$fraction else 1
        rw <- replicate(ntree, function() 
            do.call("c", tapply(idx, strata, function(i) sample(i, size = length(i) * frac, 
                   replace = perturb$replace, prob = probw[i]))))
    }

    ## apply infrastructure for determining split points
    if (is.null(applyfun)) {
        applyfun <- if(is.null(cores)) {
            lapply  
        } else {
            function(X, FUN, ...)
                parallel::mclapply(X, FUN, ..., mc.cores = cores)
        }
    }

    if (trace) pb <- txtProgressBar(style = 3) 
    forest <- applyfun(1:ntree, function(b) {
        if (trace) setTxtProgressBar(pb, b/ntree)
        tree$treefun(tree$trafo, rw[[b]], integer(0))
    })
    if (trace) close(pb)

    fitted <- data.frame(idx = idx)  
    mf <- model.frame(Formula(formula), data = tree$mf, na.action = na.pass)
    y <- model.part(Formula(formula), data = mf, lhs = 1, rhs = 0)
    if (length(y) == 1) y <- y[[1]]
    fitted[[2]] <- y
    names(fitted)[2] <- "(response)"
    fitted <- fitted[2]
    if (length(weights) > 0)
        fitted[["(weights)"]] <- weights

    ### turn subsets in weights (maybe we can avoid this?)
    rw <- lapply(rw, function(x) tabulate(x, nbins = length(idx)))

    control$applyfun <- applyfun

    ret <- constparties(nodes = forest, data = tree$mf, weights = rw,
                        fitted = fitted, terms = terms(mf), 
                        info = list(call = match.call(), control = control))
    ### ret$update <- tree$treefun # not useful
    ret$trafo <- tree$trafo
    class(ret) <- c("cforest", class(ret))

    return(ret)
}

predict.cforest <- function(object, newdata = NULL, type = c("response", "prob", "weights", "node"), 
                            OOB = FALSE, FUN = NULL, simplify = TRUE, ...) {

    responses <- object$fitted[["(response)"]]
    forest <- object$nodes
    nd <- object$data
    vmatch <- 1:ncol(nd)
    if (!is.null(newdata)) {
        nd <- model.frame(delete.response(object$terms), 
                          data = newdata, na.action = na.pass)
        OOB <- FALSE
        vmatch <- match(names(object$data), names(nd))
    }
    nam <- rownames(nd)

    type <- match.arg(type)

    ### return terminal node ids for data or newdata
    if (type == "node")
        return(lapply(forest, fitted_node, data = nd, vmatch = vmatch, ...))

    ### extract weights
    rw <- object$weights

    # w <- matrix(0L, nrow = NROW(responses), ncol = length(nam))

    applyfun <- lapply
    if (!is.null(object$info))
        applyfun <- object$info$control$applyfun

    bw <- applyfun(1:length(forest), function(b) {
        ids <- nodeids(forest[[b]], terminal = TRUE)
        fnewdata <- fitted_node(forest[[b]], nd, vmatch = vmatch, ...)
        fdata <- fitted_node(forest[[b]], object$data, ...)
        tw <- rw[[b]]
        if (OOB) tw <- as.integer(tw == 0)
        pw <- sapply(ids, function(i) tw * (fdata == i))
        return(pw[, match(fnewdata, ids), drop = FALSE])
    })

    w <- Reduce("+", bw)
    if (!is.matrix(w)) w <- matrix(w, ncol = 1)

    if (type == "weights") {
        ret <- w
        colnames(ret) <- nam
        rownames(ret) <- rownames(responses)
        return(ret)
    }
    
    pfun <- function(response) {

        if (is.null(FUN)) {

            rtype <- class(response)[1]
            if (rtype == "ordered") rtype <- "factor"
            if (rtype == "integer") rtype <- "numeric"

            FUN <- switch(rtype,
                "Surv" = if (type == "response") .pred_Surv_response else .pred_Surv,
                "factor" = if (type == "response") .pred_factor_response else .pred_factor,
                "numeric" = if (type == "response") .pred_numeric_response else .pred_ecdf)
        }

        ret <- vector(mode = "list", length = ncol(w))
        for (j in 1:ncol(w))
            ret[[j]] <- FUN(response, w[,j])
        ret <- as.array(ret)
        dim(ret) <- NULL
        names(ret) <- nam
         
        if (simplify)
            ret <- .simplify_pred(ret, names(ret), names(ret))
        ret
    }
    if (!is.data.frame(responses)) {
        ret <- pfun(responses)
    } else {
        ret <- lapply(responses, pfun)
        if (all(sapply(ret, is.atomic)))
            ret <- as.data.frame(ret)
        names(ret) <- colnames(responses)
    }
    ret
}
