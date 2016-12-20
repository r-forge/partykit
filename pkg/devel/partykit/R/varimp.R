
logLik.constparty <- function(object, newdata, weights, ...) {

    y <- object$fitted[["(response)"]]
    nd <- factor(object$fitted[["(fitted)"]])
    if (missing(newdata)) {
        pnd <- nd
    } else {
        pnd <- factor(predict(object, newdata = newdata, type = "node", ...))
    }
    ll <- switch(class(y)[1], 
           "integer" = {
               -(y - tapply(y, nd, mean)[pnd])^2
           },
           "numeric" = {
               -(y - tapply(y, nd, mean)[pnd])^2
           },
          "factor" = {
              probs <- do.call("rbind", 
                  tapply(y, nd, function(x) prop.table(table(x))))[pnd,]
              log(pmax(probs[cbind(1:length(y), unclass(y))], 
                        sqrt(.Machine$double.eps)))
          },
          "ordered" = {
              probs <- do.call("rbind", 
                  tapply(y, nd, function(x) prop.table(table(x))))[pnd,]
              log(pmax(probs[cbind(1:length(y), unclass(y))], 
                       sqrt(.Machine$double.eps)))
          },
          "Surv" = stop("not yet implemented"),
          stop("not yet implemented")   
    )
    if (missing(weights)) return(sum(ll) / length(y))
    return(sum(weights * ll) / sum(weights))
}

miscls <- function(object, newdata, weights, ...) {

    y <- object$fitted[["(response)"]]
    stopifnot(is.factor(y))
    nd <- factor(object$fitted[["(fitted)"]])
    if (missing(newdata)) {
        pnd <- nd
    } else {
        pnd <- factor(predict(object, newdata = newdata, type = "node", ...))
    }
    ll <- unclass(y) != c(tapply(y, nd, function(x) which.max(table(x))))[pnd]

    if (missing(weights)) return(sum(ll) / length(y))
    return(sum(weights * ll) / sum(weights))
}

varimp <- function(object, nperm = 1L, ...)
    UseMethod("varimp")

varimp.constparty <- function(object, nperm = 1L, risk = c("loglik", "misclassification"), conditions = NULL, 
                              mincriterion = 0, ...) {

    if (!is.function(risk)) {
        risk <- match.arg(risk)
        ### risk is _NEGATIVE_ log-likelihood
        risk <- switch(risk, "loglik" = function(...) -logLik(...),
                             "misclassification" = miscls)
    }

    if (mincriterion > 0) 
        stop("mincriterion not yet implemented") ### use nodeprune

    psplitids <- unique(do.call("c", 
        nodeapply(node_party(object), 
                  ids = nodeids(node_party(object)),
                  FUN = function(x) split_node(x)$varid)))
    vnames <- names(object$data)
    psplitvars <- vnames[psplitids]
    ret <- numeric(length(psplitvars))
    names(ret) <- psplitvars

    for (vn in psplitvars) {
        cvn <- conditions[[vn]]
        if (is.null(cvn)) {
            perm <- vn
        } else {
            blocks <- .get_psplits(object, cvn) 
            if (length(blocks) == 0) blocks <- factor(rep(1, nrow(object$data)))
            perm <- vector(mode = "list", length = 1)
            names(perm) <- vn
            perm[[vn]] <- blocks
           }
        for (p in 1:nperm)
            ret[vn] <- ret[vn] + risk(object, newdata = object$data, perm = perm, ...)
    }
    ret <- (ret - risk(object, newdata = object$data)) / nperm

    ret
}

gettree <- function(object, tree = 1L, ...)
    UseMethod("gettree")

gettree.cforest <- function(object, tree = 1L, ...) {
    ret <- party(object$nodes[[tree]], data = object$data, fitted = object$fitted)
    ret$terms <- object$terms
    class(ret) <- c("constparty", class(ret))
    ret
}

.create_cond_list <- function(object, threshold) {

    d <- object$data
    response <- names(d)[attr(object$terms, "response")]
    xnames <- all.vars(object$terms)
    xnames <- xnames[xnames != response]

    ret <- lapply(xnames, function(x) {
        tmp <- ctree(as.formula(paste(x, "~", paste(xnames[xnames != x], collapse = "+"))),
                     data = d, control = ctree_control(teststat = "quad", testtype = "Univariate",
                                                       stump = TRUE))
        pval <- info_node(node_party(tmp))$criterion["p.value",]
        pval[is.na(pval)] <- 1
        ret <- names(pval)[pval < threshold]
        if (length(ret) == 0) return(NULL)
        return(ret)
    })
    names(ret) <- xnames
    return(ret)
}

.get_psplits <- function(object, xnames) {

    d <- object$data
    ret <- lapply(xnames, function(xn) {
        id <- which(colnames(d) == xn)
        psplits <- nodeapply(node_party(object), 
            ids = nodeids(node_party(object)),
            FUN = function(x) {
                if (is.null(x)) return(NULL)
                if (is.terminal(x)) return(NULL)
                if (split_node(x)$varid == id)
                    return(split_node(x))
                return(NULL)
            })
        psplits <- psplits[!sapply(psplits, is.null)]
        if (length(psplits) > 0)
            return(do.call("interaction", lapply(psplits, kidids_split, data = d))[, drop = TRUE])
        return(NULL)
    })
    ret <- ret[!sapply(ret, is.null)]
    if (length(ret) > 0)
        return(factor(do.call("interaction", ret)[, drop = TRUE], exclude = NULL))
    return(NULL)
}

varimp.cforest <- function(object, nperm = 1L, OOB = TRUE, risk = c("loglik", "misclassification"), 
                           conditional = FALSE, threshold = .2, ...) {

    ret <- matrix(0, nrow = length(object$nodes), ncol = ncol(object$data))
    colnames(ret) <- names(object$data)

    if (conditional) {
        conditions <- .create_cond_list(object, threshold)
    } else {
        conditions <- NULL
    }

    for (b in 1:length(object$nodes)) {
        tree <- gettree(object, b)
        if (OOB) {
            oobw <- as.integer(object$weights[[b]] == 0)
            vi <- varimp(tree, nperm = nperm, risk = risk, conditions = conditions, 
                         weights = oobw, ...)
        } else {
            vi <- varimp(tree, nperm = nperm, risk = risk, conditions = conditions, 
                         ...)
        }
        ret[b, match(names(vi), colnames(ret))] <- vi
    }
    colMeans(ret)
}

