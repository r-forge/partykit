




#' Derive initial rule ensemble
#'
#' \code{initial} derives a large initial ensemble of prediction rules.
#' 
#' @param formula regression formula; a symbolic description of the model to be fit.
#' @param data a data frame containing the variables in the model.
#' @param type character. type of base learners to be included in ensemble. Defaults to "both" (intial 
#' ensemble included both rules and linear functions). Other option may be "rules" (for prediction rules 
#' only) or "linear" (for linear functions only).
#' @param weights an optional vercotr of weights to be used in deriving the trees.
#' @param sampfrac numeric. fraction of randomly selected training observations to produce each tree. 
#' Defaults to .5.
#' @param ntrees numeric. total number of trees to be generated. Defaults to 100.
#' @param seed numeric. random seed to be used in deriving the final ensemble (for reproducability). 
#' Defaults to 42.
#' @param maxdepth numeric. maximal depth of trees to be grown. Defaults to 3, resulting in trees with max 
#' 15 nodes (8 terminal and 7 inner nodes), and therefore 15 max rules.
#' @param memory.par numeric. learning rate for sequentially induced trees. Defaults to .01.
#' @param remove.duplicates logical. Remove rules from the ensemble which have the exact same support in 
#' training data?
#' @param max.rules numeric. Approximate maximum number of rules to be generated. The number of rules in 
#' the final ensemble will be smaller, due to the omission of rules with identical conditions or support.
#' @return a list with \code{data} (the original dataset), \code{rulevars} dataset with 0-1 coded rule 
#' memberships for each rule and observation, \code{rules} dataframe with rule descriptions
#' @examples airq.init.ens <- initial(Ozone ~ ., data=airquality[complete.cases(airquality),])
#' @export
#' @import partykit Formula
initial <- function(formula, data, type = "both", weights = rep(1, times = nrow(data)), sampfrac = .5, 
                    ntrees = 100, seed = 42, maxdepth = 3, memory.par = 0.01, remove.duplicates = TRUE,
                    max.rules = 2500)
{
  ###################
  ## Preliminaries ##
  ###################
  
  if(missing(formula) || missing(data)) {
    stop("Arguments formula and dataset are required")
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
    rules <- vector()
    counter <- 0
    for(i in 1:ntrees) {
      while(length(rules) <= max.rules) {
        # Take subsample of dataset
        counter <- counter + 1
        subsample <- sample(1:n, round(sampfrac * n))
        subsampledata <- data[subsample,]
        subsampledata[,y_name] <- y_learn[subsample]
        # Make sure ctree() can find object specified by weights argument: 
        environment(formula) <- environment()
        # Grow ctree (should later be glmertrees) on subsample:
        tree <- ctree(formula, data = subsampledata, weights = weights[subsample], maxdepth = maxdepth)
        # Collect rules from tree:
        rules <- append(rules, unlist(partykit:::.list.rules.party(tree)))
        # Substract predictions from current y:
        y_learn <- y_learn - memory.par * predict(tree, newdata = data)
      }
    }
    if(counter < ntrees){
      cat("\nA total of", counter, "trees were grown on subsumples when the maximum number of rules was reached. \n")
    }
    # Keep unique, non-empty rules only:
    rules <- unique(rules[!rules==""])
    
    # Create dataframe with 0-1 coded rules:
    rulevars <- data.frame(rule1 = as.numeric(with(data, eval(parse(text = rules[[1]])))))
    for(i in 2:length(rules)) {
      rulevars[,paste("rule", i, sep="")] <- as.numeric(with(data, eval(parse(text = rules[[i]]))))
    }
  }
  
  if(remove.duplicates) {
    rulevars <- rulevars[,!duplicated(t(rulevars))]
    cat("\nA total of", length(rulevars[,duplicated(t(rulevars))]), 
        "generated rules had support identical to earlier rules and were removed from the initial ensemble. \n")
    duplicates.removed <- colnames(rulevars[,duplicated(t(rulevars))])
  } else {
    duplicates.removed <- NULL
  }
  
  ####################
  ## Return results ##
  ####################
  
  cat("\nA total of", ncol(rulevars), "rules were included in the initial ensemble.")  
  result <- list(data = data, rulevars = rulevars, call = match.call(), weights = weights,
                 x_names = x_names, y_name = y_name, type = type,
                 remove.duplicates = remove.duplicates, duplicates.removed = duplicates.removed, 
                 trees = trees, rules = data.frame(rule = paste("rule", 1:length(rules), sep = ""), 
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
#' @param thres numeric. threshold for convergence. 
#' @param standardize logical. Standardize rules and predictor variables before estimating the regression
#' model?
#' @param wins.frac numeric. Quantiles of data distribution to be used for winsorizing linear predictors.
#' When set to 0, no winsorizing is performed.
#' @param normalize logical. Normalize linear variables before estimating the regression model? Normalizing gives lienar terms the same a priori influence as a typical rule.
#' @param nfolds numeric. Number of folds to be used in performing cross validation for determining penalty parameter
#' @param mod.sel.crit character.  
#' @return a list with many elements 
#' @examples airq.init.ens <- initial(Ozone ~ ., data=airquality[complete.cases(airquality),])
#' airq.final.ens <- final(airq.init.ens)
#' @export
#' @import glmnet
final <- function(object, seed = 42, alpha = 1, dfmax = ncol(object$data) + ncol(object$rulevars) + 1, 
                  thres = 1e-07, standardize = FALSE, wins.frac = .025,
                  normalize = TRUE, nfolds = 10, mod.sel.crit = "deviance") 
{ 
  
  ###################
  ## Preliminaries ##
  ###################
  
  if(length(wins.frac) != 1 || wins.frac < 0 || wins.frac > 0.5) {
    stop("Bad value for 'wins.fraction'.")
  }
  set.seed(seed)
  type <- object$type
  
  ######################################################
  ## Prepare rules, linear terms and outcome variable ##
  ######################################################
  
  x <- object$data[,object$x_names]
  rulevars <- object$rulevars
  
  if(type != "rules") {
    # Winsorize variables (section 5 of F&P(2008)):
    if(wins.frac > 0) {
      for(i in 1:ncol(x)) {
        if(is.numeric(x[,i])) { 
          lim <- quantile(x[,i], probs = c(wins.frac, 1-wins.frac))
          x[x[,i] < lim[1],i] <- lim[1]
          x[x[,i] > lim[2],i] <- lim[2]
        }
      }
    }
    x_sds <- apply(x[,apply(x, 2, is.numeric)], 2, sd, na.rm = TRUE)
    if(normalize) { 
      # Normalize linear terms (section 5 of F&P(2008)):
      x_scales <- apply(x[,apply(x, 2, is.numeric)], 2, sd, na.rm = TRUE)/0.4
      x <- scale(x[,apply(x, 2, is.numeric)], center = FALSE, scale = x_scales)
    } else {
      x_scales <- NULL
    } 
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
  
  ##################################################
  ## Perform penalized regression on the ensemble ##
  ##################################################
  
  glmnet.fit <- cv.glmnet(as.matrix(x), y, alpha = alpha, standardize = standardize, nfolds = nfolds,
                          type.measure = mod.sel.crit, thres = thres, dfmax = dfmax, 
                          weights = object$weights)
  
  ####################
  ## Return results ##
  ####################
  
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
                 data = object$data, normalize = normalize, x_scales = x_scales, type = type,
                 rules = object$rules, x_names = colnames(x), initial.ens = object, x_sds = x_sds,
                 rulevars = rulevars, varnames = object$x_names)
  class(result) <- "final"
  return(result)
}
