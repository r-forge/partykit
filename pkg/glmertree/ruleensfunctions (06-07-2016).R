#' Derive initial rule ensemble
#'
#' \code{initial} derives a large initial ensemble of prediction rules.
#' 
#' @param formula regression formula; a symbolic description of the model to be fit.
#' @param data a data frame containing the variables in the model.
#' @param type character. type of base learners to be included in ensemble. Defaults to "both" (intial ensemble included both rules and linear functions). Other option may be "rules" (for prediction rules only) or "linear" (for linear functions only).
#' @param weights an optional vercotr of weights to be used in deriving the trees.
#' @param sampfrac numeric. fraction of randomly selected training observations to produce each tree. Defaults to .5.
#' @param ntrees numeric. total number of trees to be generated. Defaults to 100.
#' @param seed numeric. random seed to be used in deriving the final ensemble (for reproducability). Defaults to 42.
#' @param maxdepth numeric. maximal depth of trees to be grown. Defaults to 3, resulting in trees with max 8 terminal nodes.
#' @param memory.par numeric. learning rate for sequentially induced trees. Defaults to .01.
#' @return a list with \code{data} (the original dataset), \code{rulevars} dataset with 0-1 coded rule memberships for each rule and observation, \code{rules} dataframe with rule descriptions
#' @examples airq.init.ens <- initial(Ozone ~ ., data=airquality[complete.cases(airquality),])
initial <- function(formula, data, type = "both", weights = rep(1, times = nrow(data)), sampfrac = .5, 
                    ntrees = 100, seed = 42, maxdepth = 3, memory.par = 0.01)
{
  ###################
  ## Preliminaries ##
  ###################
  
  if(missing(formula) || missing(data)) {
    stop("Arguments formula and dataset are required")
  }
  if(!is.formula(formula) & !Formula::is.Formula(formula)) {
    stop("Argument formula should be of class formula")
  }
  if(length(sampfrac) != 1 || sampfrac < 0.01 || sampfrac > 0.99) {
    stop("Bad value for 'sampfrac'")
  }
  if(length(type) != 1 || (type != "rules" & type != "both" & type != "linear")) {
    stop("Argument type should equal 'both', 'rules' or 'linear'")
  }
  set.seed(seed)
  y_name <- paste(attr(Formula::as.Formula(formula), "lhs"))
  y_orig <- data[,y_name]
  y_learn <- data[,y_name]
  x_names <- attr(terms(formula, data = data), "term.labels")
  n <- nrow(data)
  
  #####################################
  ## Step 1: Derive prediction rules ##
  #####################################
  
  if(type != "linear") {
    trees <- list()
    rules <- list()
    subsamples <- list()
    for(i in 1:ntrees) {
      # Take subsample of dataset
      subsamples[[i]] <- sample(1:n, round(sampfrac * n))
      subsampledata <- data[subsamples[[i]],]
      subsampledata[,y_name] <- y_learn[subsamples[[i]]]
      # Make sure ctree() can find object specified by weights argument 
      # (ctree will look for objects in specified dataset and environment(formula)):
      environment(formula) <- environment()
      # Grow ctree (should later be glmertrees) on subsample:
      trees[[i]] <- ctree(formula, data = subsampledata, weights = weights[subsamples[[i]]], maxdepth = maxdepth)
      # Collect rules from tree:
      rules[[i]] <- partykit:::.list.rules.party(trees[[i]])
      # Substract predictions from current y:
      y_learn <- y_learn - memory.par * predict(trees[[i]], newdata = data)
    }
    # Keep unique rules only:
    rules <- unique(unlist(rules))
    # Keep non-empty rules only:
    rules <- rules[!rules==""]
    # To do: Should also test for rules that have exact same support!
    
    # Create dataframe with 0-1 coded rules:
    rulevars <- data.frame(rule1 = as.numeric(with(data, eval(parse(text = rules[[1]])))))
    for(i in 2:length(rules)) {
      rulevars[,paste("rule", i, sep="")] <- as.numeric(with(data, eval(parse(text = rules[[i]]))))
    }
  }
  
  ####################
  ## Return results ##
  ####################
  
  result <- list(data = data, rulevars = rulevars, call = match.call(), weights = weights,
                 x_names = x_names, y_name = y_name, weights = weights, type = type,
                 rules = data.frame(rule = paste("rule", 1:length(rules), sep = ""), 
                                    descriptions = rules))
  return(result)
}


#' Derive final rule ensemble
#'
#' \code{final} derives a final prediction rule ensemble with penalized regression
#' 
#' @param object an object resulting from application of \code{intial}.
#' @param seed numeric. set random seed for replicability. Defaults to 42.
#' @param alpha numeric. elastic net mixing parameter. 
#' @param dfmax numeric. maximal number of terms in the final model.
#' @param thres numeric. threshold for 
#' @param standardize logical. Standardize rules and predictor variables before estimating the regression model?
#' @param winsorize logical. Winsorize linear terms?
#' @param wins.frac numeric. numeric. Only used when winsorize = TRUE, quantiles of data distribution to be used to winsorize variables 
#' @param normalize logical. Normalize linear variables before entering estimating the regression model? Normalizing gives lienar terms the same a priori influence as a typical rule.
#' @param nfolds numeric. number of folds to be used in performing cross validation for determining penalty parameter
#' @param mod.sel.crit character. 
#' @param remove.duplicates logical. Remove rules from the ensemble which have the exact same support in training data?
#' @return a list with 
#' @examples airq.init.ens <- initial(Ozone ~ ., data=airquality[complete.cases(airquality),])
#' airq.final.ens1 <- final(airq.init.ens, normalize = FALSE, winsorize = FALSE, standardize = TRUE)
#' coef(airq.final.ens1$glmnet.fit, s="lambda.min")
#' plot(airq.final.ens1$glmnet.fit)

final <- function(object, seed = 42, alpha = 1, dfmax = ncol(object$data) + ncol(object$rulevars) + 1, 
                  thres = 1e-07, standardize = FALSE, winsorize = TRUE, wins.frac = .025,
                  normalize = TRUE, nfolds = 10, mod.sel.crit = "deviance", remove.duplicates = TRUE) 
{ 
  if(length(wins.frac) != 1 || wins.frac < 0 || wins.frac > 0.5) {
    stop("Bad value for 'wins.fraction'.")
  }
  set.seed(seed)
  type <- object$type
  
  if(remove.duplicates) {
    rulevars <- object$rulevars[,!duplicated(t(object$rulevars))]
  } else {
    rulevars <- object$rulevars
  }
  x <- object$data[,object$x_names]
  
  # Derive linear terms (if type is 'linear' or 'both'):
  if(type != "rules") {
    # Winsorize variables (section 5 of F&P(2008)):
    # To do: only do this for numeric variables!
    if(winsorize) {
      for(i in 1:ncol(x)) {
        lim <- quantile(x[,i], probs = c(wins.frac, 1-wins.frac))
        x[x[,i] < lim[1],i] <- lim[1]
        x[x[,i] > lim[2],i] <- lim[2]
      }
    }
    if(normalize) { 
      # Normalize linear terms (section 5 of F&P(2008)):
      for(i in ncol(x)) {
        x[,i] <- 0.4 * x[,i] / sd(x[,i], na.rm = TRUE)
      }
    } # To do: For predictions, coefficients of linear terms should be unnormalized!
  }
  
  if(type == "both") {
    x <- as.matrix(cbind(rulevars, x))
  }
  if(type == "linear") {
    x <- as.matrix(x)
  }
  if(type == "rules") {
    x <- as.matrix(rulevars)
  }
  y <- object$data[,object$y_name]
  glmnet.fit <- cv.glmnet(as.matrix(x), y, alpha = alpha, standardize = standardize, nfolds = nfolds,
                          type.measure = mod.sel.crit, thres = thres, dfmax = dfmax, 
                          weights = object$weights)
  
  lmin_ind <- which(glmnet.fit$lambda == glmnet.fit$lambda.min)
  l1se_ind <- which(glmnet.fit$lambda == glmnet.fit$lambda.1se)
  cat("Model with minimum cv error: \n  lambda = ", glmnet.fit$lambda[lmin_ind],  
      "\n  number of terms = ", glmnet.fit$nzero[lmin_ind], 
      "\n  mean cv error (se) = ", glmnet.fit$cvm[lmin_ind], 
      " (", glmnet.fit$cvsd[lmin_ind], ")", 
      "\n\nModel with cv error within 1se of minimum: \n  lambda = ", glmnet.fit$lambda[l1se_ind],  
      "\n  number of terms = ", glmnet.fit$nzero[l1se_ind], 
      "\n  mean cv error (se) = ", glmnet.fit$cvm[l1se_ind], 
      " (", glmnet.fit$cvsd[l1se_ind], ")", sep="")
  
  result <- list(glmnet.fit = glmnet.fit, call = match.call(), weights = object$weights, 
                 data = object$data, winsorize = winsorize, normalize = normalize, 
                 rules = object$rules, x_names = colnames(x))
  return(result)
}





coef.final <- function(object, penalty.par.val = "lambda.1se")
{
  if(object$normalize){} # to put coefficients of linear terms in original metric, use std(x)*coef/.4
  coefs <- as(coef(object$glmnet.fit, s = penalty.par.val), Class = "matrix")
  coefs <- data.frame(coefs = round(coefs[,1], digits = 4), rule = rownames(coefs))
  coefs <- merge(coefs, object$rules, all.x=T)
  coefs <- coefs[coefs$coefs != 0,]
  coefs[order(abs(coefs$coefs), decreasing = TRUE),]
}




predict.final <- function(object, newdata = NULL, penalty.par.val = "lambda.1se")
{
  if(is.null(newdata)){newdata <- object$data}
  rules <- as.character(object$rules$descriptions)
  for(i in 1:length(rules)) {
    newdata[,paste("rule", i, sep="")] <- as.numeric(with(newdata, eval(parse(text = rules[[i]]))))
  }
  return(predict.cv.glmnet(object$glmnet.fit, newx = as.matrix(newdata[object$x_names]), 
                           s = penalty.par.val))
}