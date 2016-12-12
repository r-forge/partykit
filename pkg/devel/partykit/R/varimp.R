
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
               (y - tapply(y, nd, mean)[pnd])^2
           },
           "numeric" = {
               (y - tapply(y, nd, mean)[pnd])^2
           },
          "factor" = {
              probs <- do.call("rbind", tapply(y, nd, function(x) prop.table(table(x))))[pnd,]
              -log(pmax(probs[cbind(1:length(y), unclass(y))], sqrt(.Machine$double.eps)))
          },
          "ordered" = {
              probs <- do.call("rbind", tapply(y, nd, function(x) prop.table(table(x))))[pnd,]
              -log(pmax(probs[cbind(1:length(y), unclass(y))], sqrt(.Machine$double.eps)))
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

varimp.constparty <- function(object, nperm = 1L, weights, risk = logLik, conditions = NULL, 
                              mincriterion = 0, ...) {

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

    for (p in 1:nperm) {
        for (vn in psplitvars) {
            if (is.null(conditions[[vn]])) {
                ret[vn] <- ret[vn] + risk(object, newdata = object$data, 
                                          perm = vn, ...)
            } else {
                perm <- do.call("c", tapply(1:nrow(object$data), conditions[[vn]], sample))
                perm <- list(perm)
                names(perm) <- vn
                ret[vn] <- ret[vn] + risk(object, newdata = object$data, 
                                          perm = perm, ...)
           }
        }
    }
    ret <- (ret - risk(object)) / nperm

    ret
}

varimp.cforest <- function(object, nperm = 1L, OOB = TRUE, risk = logLik, ...) {

    ret <- matrix(0, nrow = length(object$nodes), ncol = ncol(object$data))
    colnames(ret) <- names(object$data)

    for (b in 1:length(object$nodes)) {
        tree <- party(object$nodes[[b]], data = object$data, fitted = object$fitted)
        class(tree) <- c("constparty", class(tree))
        if (OOB) {
            oobw <- as.integer(object$weights[[b]] == 0)
            vi <- varimp(tree, nperm = nperm, risk = risk, weights = oobw, ...)
        } else {
            vi <- varimp(tree, nperm = nperm, risk = risk, ...)
        }
        ret[b, match(names(vi), colnames(ret))] <- vi
    }
    colMeans(ret)
}

library("partykit")

### regression
airq <- subset(airquality, !is.na(Ozone))
airct <- ctree(Ozone ~ ., data = airq)

mean((airq$Ozone - predict(airct))^2)
logLik(airct)
logLik(airct, airq, perm = "Temp")

varimp(airct)

aircf <- cforest(Ozone ~ ., data = airq, ntree = 100)

varimp(aircf)

#library("party")
#
#aircf2 <- party::cforest(Ozone ~ ., data = airq)
#party::varimp(aircf2)

ict <- cforest(Species ~ ., data = iris)
varimp(ict)
varimp(ict, risk = miscls)

library("party")

ict2 <- party::cforest(Species ~ ., data = iris)
party::varimp(ict2)
