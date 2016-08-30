
.glmtrafo <- function(formula, data, ctrl, converged = NULL) {

    weights <- model.weights(data)
    if (is.null(weights)) weights <- integer(0)
    cluster <- data[["cluster"]]
    offset <- model.offset(data)

    ### <FIXME> handle offset and cluster </FIXME>

    if (ctrl$nmax < Inf) {
        if (!is.null(cluster)) stop("cluster not implemented")   
        mf <- model.frame(formula, data, na.action = na.pass)
        bdr <- BDR::BDR(mf, complete.cases.only = TRUE, total = TRUE)
        mf2 <- as.data.frame(bdr)
        iy <- c(bdr)
        attr(iy, "levels") <- 1:nrow(mf2)
        mfs <- model.frame(formula, data = mf2)
        y <- model.response(mfs)
        x <- model.matrix(formula, data = mf2)
        return(function(subset) {
            w <- c(libcoin::ctabs(iy, weights = weights, subset = subset)[-1L])
            mod <- glm(y ~ x + 0, family = ctrl$family, weights = w)
            Y <- sandwich::estfun(mod)
            Y <- Y / w
            Y[w == 0,] <- 0
            ret <- rbind(0, Y)
            list(estfun = ret, index = iy, coef = coef(mod), logLik = logLik(mod),
                 converged = if (is.null(converged)) 
                     mod$converged else converged(mod, mf, subset))
        })
    }
    if (!is.null(cluster)) stop("cluster not implemented")
    mf <- model.frame(formula, data, na.action = na.pass)
    cc <- complete.cases(mf)
    y <- model.response(mf)
    x <- model.matrix(formula, data = mf)
    return(function(subset) {
        s <- subset[cc[subset]]
        ys <- y[s]
        xs <- x[s, , drop = FALSE]
        if (length(weights) > 0) {
            w <- weights[cc[subset]]
            mod <- glm(ys ~ xs + 0, family = ctrl$family, weights = w)
        } else {
             mod <- glm(ys ~ xs + 0, family = ctrl$family)
        }
        ret <- matrix(0, nrow = NROW(x), ncol = NCOL(x))
        Y <- sandwich::estfun(mod)
        if (length(weights) > 0) {
            Y <- Y / w
            Y[w == 0,] <- 0
        }
        ret[subset,] <- Y
        storage.mode(ret) <- "double"
        list(estfun = ret, coef = coef(mod), logLik = logLik(mod),
             converged = if (is.null(converged)) 
                 mod$converged else converged(mod, mf, subset))
    })
}

glmtree <- function(formula, data, weights, subset, na.action = na.pass, 
                    family = gaussian(), control = ctree_control(), ...) {

    call <- match.call(expand.dots = FALSE)
    call$na.action <- na.action
    frame <- parent.frame()
    control$family <- family

    tree <- .urp_tree(call, frame, control = control,
                      growfun = .ctreegrow, trafofun = .glmtrafo,
                      doFit = TRUE)

    ### <FIXME> change this to modelparty 
    mf <- tree$mf
    weights <- model.weights(mf)
    if (is.null(weights)) weights <- rep(1, nrow(mf))

    fitted <- data.frame("(fitted)" = fitted_node(tree$node, mf),
                         "(weights)" = weights,
                         check.names = FALSE)
    y <- model.part(Formula(formula), data = mf, lhs = 1, rhs = 0)
    if (length(y) == 1) y <- y[[1]]
    fitted[[3]] <- y
    names(fitted)[3] <- "(response)"
    ret <- party(tree$node, data = mf, fitted = fitted,
                 info = list(call = match.call(), control = control))
    ret$update <- tree$treefun
    ret$trafo <- tree$trafo
    class(ret) <- c("constparty", class(ret))

    ### doesn't work for Surv objects
    # ret$terms <- terms(formula, data = mf)
    ret$terms <- tree$terms
    ### need to adjust print and plot methods
    ### for multivariate responses
    ### if (length(response) > 1) class(ret) <- "party"
    ### </FIXME>
    return(ret)
}

