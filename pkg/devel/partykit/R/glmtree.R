## simple wrapper function to specify fitter and return class
glmtree <- function(formula, data, subset, na.action, weights, offset, cluster,
  family = gaussian, epsilon = 1e-8, maxit = 25, ...)
{
  ## use dots for setting up mob_control
  control <- mob_control(...)

  ## keep call
  cl <- match.call(expand.dots = TRUE)

  ## extend formula if necessary
  f <- Formula::Formula(formula)
  if(length(f)[2L] == 1L) {
    attr(f, "rhs") <- c(list(1), attr(f, "rhs"))
    formula[[3L]] <- formula(f)[[3L]]
  } else {
    f <- NULL
  }

  ## process family
  if(inherits(family, "family")) {
    fam <- TRUE
  } else {
    fam <- FALSE
    if(is.character(family)) family <- get(family)
    if(is.function(family)) family <- family()
  }

  ## call mob
  m <- match.call(expand.dots = FALSE)
  if(!is.null(f)) m$formula <- formula
  m$fit <- glmfit
  m$control <- control
  m$epsilon <- epsilon
  m$maxit <- maxit
  if("..." %in% names(m)) m[["..."]] <- NULL
  if(!fam) m$family <- family
  m[[1L]] <- as.call(quote(partykit::mob))
  rval <- eval(m, parent.frame())
  
  ## extend class and keep original call
  rval$info$call <- cl
  rval$info$family <- family$family
  class(rval) <- c("glmtree", class(rval))
  return(rval)
}

## actual fitting function for mob()
glmfit <- function(y, x, start = NULL, weights = NULL, offset = NULL, cluster = NULL, ...,
  estfun = FALSE, object = FALSE)
{
  ## catch control arguments
  args <- list(...)
  control <- list()
  for(n in c("epsilon", "maxit")) {
    if(n %in% names(args)) {
      control[[n]] <- args[[n]]
      args[[n]] <- NULL
    }
  }
  args$control <- do.call("glm.control", control)
  
  ## add intercept-only regressor matrix (if missing)
  ## NOTE: does not have terms/formula
  if(is.null(x)) x <- matrix(1, nrow = NROW(y), ncol = 1L,
    dimnames = list(NULL, "(Intercept)"))
  
  ## call glm fitting function
  args <- c(list(x = x, y = y, start = start, weights = weights, offset = offset), args)
  z <- do.call("glm.fit", args)

  ## degrees of freedom
  df <- z$rank
  if(z$family$family %in% c("gaussian", "Gamma", "inverse.gaussian")) df <- df + 1
  if(substr(z$family$family, 1L, 5L) != "quasi") objfun <- z$aic/2 - df else objfun <- z$deviance

  ## list structure
  rval <- list(
    coefficients = z$coefficients,
    objfun = objfun,
    estfun = NULL,
    object = NULL
  )

  ## add estimating functions (if desired)
  if(estfun) {
    wres <- as.vector(z$residuals) * z$weights
    dispersion <- if(substr(z$family$family, 1L, 17L) %in% c("poisson", "binomial", "Negative Binomial")) {
      1
    } else {
      sum(wres^2, na.rm = TRUE)/sum(z$weights, na.rm = TRUE)
    }
    rval$estfun <- wres * x/dispersion
  }

  ## add model (if desired)
  if(object) {
    class(z) <- c("glm", "lm")
    z$offset <- if(is.null(offset)) 0 else offset
    z$contrasts <- attr(x, "contrasts")
    z$xlevels <- attr(x, "xlevels")    

    cl <- as.call(expression(glm))
    cl$formula <- attr(x, "formula")	
    if(!is.null(offset)) cl$offset <- attr(x, "offset")
    z$call <- cl
    z$terms <- attr(x, "terms")

    rval$object <- z
  }

  return(rval)
}

## methods
print.glmtree <- function(x,
  title = NULL, objfun = NULL, ...)
{
  if(is.null(title)) title <- sprintf("Generalized linear model tree (family: %s)", x$info$family)
  if(is.null(objfun)) objfun <- if(substr(x$info$family, 1L, 5L) != "quasi") "negative log-likelihood" else "deviance"
  print.modelparty(x, title = title, objfun = objfun, ...)
}

predict.glmtree <- function(object, newdata = NULL, type = "response", ...)
{
  ## FIXME: possible to get default?
  if(is.null(newdata) & !identical(type, "node")) stop("newdata has to be provided")
  predict.modelparty(object, newdata = newdata, type = type, ...)
}

plot.glmtree <- function(x, terminal_panel = node_bivplot,
  tp_args = list(), tnex = NULL, drop_terminal = NULL, ...)
{
  nreg <- if(is.null(tp_args$which)) x$info$nreg else length(tp_args$which)
  if(nreg < 1L & missing(terminal_panel)) {
    plot.constparty(as.constparty(x),
      tp_args = tp_args, tnex = tnex, drop_terminal = drop_terminal, ...)
  } else {
    if(is.null(tnex)) tnex <- if(is.null(terminal_panel)) 1L else 2L * nreg
    if(is.null(drop_terminal)) drop_terminal <- !is.null(terminal_panel)
    plot.modelparty(x, terminal_panel = terminal_panel,
      tp_args = tp_args, tnex = tnex, drop_terminal = drop_terminal, ...)
  }
}

### proof of concept for reimplementation of glmtree using urp_tree
.glmtrafo <- function(formula, data, ctrl, converged = NULL) {

    weights <- model.weights(data)
    if (is.null(weights)) weights <- integer(0)
    cluster <- data[["cluster"]]
    offset <- model.offset(data)

    ### <FIXME> handle offset and cluster </FIXME>

    if (ctrl$nmax < Inf) {
        if (!is.null(cluster)) stop("cluster not implemented")   
        mf <- model.frame(formula, data, na.action = na.pass)
        bdr <- inum::inum(mf, complete.cases.only = TRUE, total = TRUE)
        mf2 <- as.data.frame(bdr)
        iy <- c(bdr)
        attr(iy, "levels") <- 1:nrow(mf2)
        mfs <- model.frame(formula, data = mf2)
        y <- model.response(mfs)
        x <- model.matrix(formula, data = mf2)
        
        glmfit <- function(subset, estfun = TRUE, object = FALSE, info = NULL, ...) {
          nobs <- length(subset)
          w <- c(libcoin::ctabs(iy, weights = weights, subset = subset)[-1L])
          mod <- glm.fit(x = x, y = y, weights = w, start = info$coef, family = ctrl$family)
          
          ret <- NULL
          if(estfun) {
            wres <- as.vector(mod$residuals) * mod$weights
            dispersion <- if(substr(mod$family$family, 1L, 17L) %in% c("poisson", "binomial", "Negative Binomial")) {
              1
            } else {
              sum(wres^2, na.rm = TRUE)/sum(mod$weights, na.rm = TRUE)
            }
            Y <- wres * x/dispersion
            Y <- Y / w
            Y[w == 0,] <- 0
            ret <- rbind(0, Y)
          }
          
          # mod <- glm(y ~ x + 0, family = ctrl$family, weights = w, start = info$coef)
          # ret <- NULL
          # if (estfun) {
          #   Y <- sandwich::estfun(mod)
          #   Y <- Y / w
          #   Y[w == 0,] <- 0
          #   ret <- rbind(0, Y)
          # }
          class(mod) <- c("glm", "lm")
          list(estfun = ret, index = iy, coefficients = coef(mod), objfun = logLik(mod),
               object = if (object) mod else NULL, nobs = nobs,
               converged = if (is.null(converged)) 
                 mod$converged else converged(mod, mf, subset))
        }
        return(glmfit)
    }
    if (!is.null(cluster)) stop("cluster not implemented")
    mf <- model.frame(formula, data, na.action = na.pass)
    cc <- complete.cases(mf)
    y <- model.response(mf)
    x <- model.matrix(formula, data = mf)
    glmfit <- function(subset, estfun = TRUE, object = FALSE, info = NULL, ...) {
        s <- subset[cc[subset]]
        ys <- y[s]
        xs <- x[s, , drop = FALSE]
        if (length(weights) > 0) {
            w <- weights[cc[subset]]
            nobs <- sum(w)
            mod <- glm.fit(x = xs, y = ys, weights = w, start = info$coef, family = ctrl$family)
            # mod <- glm(ys ~ xs + 0, family = ctrl$family, weights = w, start = info$coef)
        } else {
             nobs <- NROW(ys)
             # mod <- glm(ys ~ xs + 0, family = ctrl$family, start = info$coef)
             mod <- glm.fit(x = xs, y = ys, start = info$coef, family = ctrl$family)
        }
        

        
        ## degrees of freedom
        df <- mod$rank
        if(mod$family$family %in% c("gaussian", "Gamma", "inverse.gaussian")) df <- df + 1
        if(substr(mod$family$family, 1L, 5L) != "quasi") objfun <- mod$aic/2 - df else objfun <- mod$deviance

        
        ## add estimating functions (if desired)
        ret <- NULL
        if(estfun) {
          ret <- matrix(0, nrow = NROW(x), ncol = NCOL(x))
          wres <- as.vector(mod$residuals) * mod$weights
          dispersion <- if(substr(mod$family$family, 1L, 17L) %in% c("poisson", "binomial", "Negative Binomial")) {
            1
          } else {
            sum(wres^2, na.rm = TRUE)/sum(mod$weights, na.rm = TRUE)
          }
          Y <- wres * xs/dispersion
          ret[subset,] <- Y
          storage.mode(ret) <- "double"
        }
        # ret <- NULL
        # if (estfun) {
        #     ret <- matrix(0, nrow = NROW(x), ncol = NCOL(x))
        #     Y <- sandwich::estfun(mod)
        #     if (length(weights) > 0) {
        #         Y <- Y / w
        #         Y[w == 0,] <- 0
        #     }
        #     ret[subset,] <- Y
        #     storage.mode(ret) <- "double"
        # }
        
        class(mod) <- c("glm", "lm")
        ## add model (if desired)
        if(object) {
          mod$offset <- if(is.null(offset)) 0 else offset
          mod$contrasts <- attr(xs, "contrasts")
          mod$xlevels <- attr(xs, "xlevels")    
          
          cl <- as.call(expression(glm))
          cl$formula <- attr(xs, "formula")	
          if(!is.null(offset)) cl$offset <- attr(x, "offset")
          mod$call <- cl
          mod$terms <- attr(xs, "terms")
        }
        
        

        list(estfun = ret, coefficients = coef(mod), objfun = logLik(mod),
             object = if (object) mod else NULL, nobs = nobs,
             converged = if (is.null(converged)) 
                 mod$converged else converged(mod, mf, subset))
    }
    return(glmfit)
}




glmtree2 <- function
(
    formula, 
    data, 
    weights, 
    subset,
    offset,
    cluster, 
    na.action = na.pass, 
    family = gaussian(), 
    epsilon = 1e-8, ## TODO: make use of this
    maxit = 25, ## TODO: make use of this
    converged = NULL,
    scores = NULL,
    ...
) {
  ## use dots for setting up mob_control
  control <- mob2_control(...)
  control$family <- family

  mob2(fit = .glmtrafo, formula = formula, data = data, weights = weights,
       subset = subset, offset = offset, cluster = cluster, na.action = na.action,
       control = control, converged = converged, scores = scores, ...)
}


mob2 <- function
(
  fit,
  formula, 
  data, 
  weights, 
  subset,
  offset,
  cluster, 
  na.action = na.pass, 
  control = mob2_control(...), 
  converged = NULL,
  scores = NULL,
  ...
) {
  
  
  ### make sure right criterion is used for exhaustive search
  if(control$testflavour == "exhaustive"){
    control$criterion <- "statistic"
    control$logmincriterion <- -Inf
  }
  
  
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
  
  trafofun <- function(...) fit(..., converged = converged)
  tree <- .urp_tree(call, frame, data = data, data_asis = data_asis, control = control,
                    trafofun = trafofun, doFit = TRUE)
  
  ### prepare as modelparty
  mf <- tree$mf
  weights <- model.weights(mf)
  if (is.null(weights)) weights <- rep(1, nrow(mf))
  mtY <- terms(tree$modelf, data = mf)
  mtZ <- delete.response(terms(tree$partf, data = mf))
  
  fitted <- data.frame("(fitted)" = fitted_node(tree$nodes, mf),
                       "(weights)" = weights,
                       check.names = FALSE)
  y <- model.part(as.Formula(tree$modelf), data = mf, lhs = 1, rhs = 0)
  if (length(y) == 1) y <- y[[1]]
  fitted[[3]] <- y
  names(fitted)[3] <- "(response)"
  
  ## return party object
  rval <- party(tree$nodes, 
                data = mf,
                fitted = fitted,
                terms = tree$terms,
                info = list(
                  call = match.call(),
                  formula = formula,
                  Formula = as.Formula(formula),
                  terms = list(response = mtY, partitioning = mtZ),
                  fit = fit,
                  control = control,
                  dots = list(...)
                  # nreg = max(0L, as.integer(xreg) * (nyx - NCOL(Y))))
                )
  )
  class(rval) <- c("modelparty", class(rval))
  
  ### add modelinfo if not there yet 
  # TODO: check if this can be done prettier
  terminals <- nodeids(rval, terminal = TRUE)
  idx <- lapply(terminals, .get_path, obj = tree$nodes)
  tree_ret <- unclass(rval)
  subset_term <- predict(rval, type = "node")
  
  for (i in 1:length(idx)) {
    
    if(is.null(tree_ret[[c(1, idx[[i]])]]$info)) {
      tree_ret[[c(1, idx[[i]])]]$info <- tree$trafo(subset = which(subset_term == terminals[i]), 
                                                    estfun = FALSE)
    }
  }
  
  class(tree_ret) <- class(rval)

  return(tree_ret)

}
