#' Derive an unbiased prediction rule ensemble
#'
#' \code{upre} derives a sparse ensemble of rules and/or linear functions for prediction
#' 
#' @param formula regression formula; a symbolic description of the model to be fit.
#' @param data matrix or data frame containing the variables in the model.
#' @param type character. Type of base learners to be included in ensemble. Defaults to "both" (intial 
#' ensemble included both rules and linear functions). Other option may be "rules" (for prediction rules 
#' only) or "linear" (for linear functions only).
#' @param weights an optional vector of observation weights to be used for deriving the ensemble.
#' @param sampfrac numeric value greater than 0, and smaller than or equal to 1. Fraction of randomly
#' selected training observations used to produce each tree. Setting this to values < 1 will result in 
#' subsamples being drawn without replacement (i.e., subsampling). Setting this equal to 1 will result in 
#' bootstrap sampling.
#' @param ntrees numeric. Total number of trees to be generated.
#' @param seed numeric. Random seed to be used in deriving the final ensemble (for reproducability). 
#' Defaults to 42.
#' @param maxdepth numeric. Maximal depth of trees to be grown. Defaults to 3, resulting in trees with max 
#' 15 nodes (8 terminal and 7 inner nodes), and therefore max 15 rules.
#' @param learnrate numeric. Learning rate for sequentially induced trees.
#' @param removeduplicates logical. Remove rules from the ensemble which have the exact same support in 
#' training data?
#' @param maxrules numeric. Approximate maximum number of rules to be generated. The number of rules in 
#' the final ensemble will be smaller, due to the omission of rules with identical conditions or support.
#' @param alpha numeric. Elastic net mixing parameter. 
#' @param dfmax numeric. Maximal number of terms in the final model.
#' @param mtry numeric. Number of randomly selected predictor variables for creating each split in each 
#' tree.
#' @param thres numeric. Threshold for convergence. 
#' @param standardize logical. Standardize rules and predictor variables before estimating the regression
#' model?
#' @param winsfrac numeric. Quantiles of data distribution to be used for winsorizing linear predictors.
#' When set to 0, no winsorizing is performed.
#' @param normalize logical. Normalize linear variables before estimating the regression model? 
#' Normalizing gives linear terms the same a priori influence as a typical rule.
#' @param nfolds numeric. Number of folds to be used in performing cross validation for determining 
#' penalty parameter
#' @param mod.sel.crit character. Model selection criterion to be used for deriving the final ensemble. 
#' The default is type.measure="deviance", which uses squared-error for gaussian models (a.k.a 
#' type.measure="mse"). type.measure="mse" or type.measure="mae" (mean absolute error) measure the 
#' deviation from the fitted mean to the response.
#' @param verbose logical. Should information on the initial and final ensemble be printed to the command line?
#' @return a list with many elements 
#' @examples 
#' \dontrun{
#'  airq.ens <- upre(Ozone ~ ., data=airquality[complete.cases(airquality),])
#' }
#' @import glmnet partykit Formula datasets
#' @details Note that variable names supplied in the formula may not start with the word 'rule' 
#' @export
upre <- function(formula, data, type = "both", weights = rep(1, times = nrow(data)), sampfrac = .5, 
                    ntrees = 500, seed = 42, maxdepth = 3, learnrate = 0.01, removeduplicates = TRUE,
                    maxrules = 2000, alpha = 1, dfmax = 5*ncol(data), mtry = Inf, 
                    thres = 1e-07, standardize = FALSE, winsfrac = .025, normalize = TRUE, nfolds = 10, 
                    mod.sel.crit = "deviance", verbose = TRUE)   
{
  ###################
  ## Preliminaries ##
  ###################
  
  if(missing(formula) || missing(data)) {
    stop("Arguments formula and dataset are required")
  }
  if(length(sampfrac) != 1 || sampfrac < 0.01 || sampfrac > 1) {
    stop("Bad value for 'sampfrac'")
  }
  if(length(type) != 1 || (type != "rules" & type != "both" & type != "linear")) {
    stop("Argument type should equal 'both', 'rules' or 'linear'")
  }
  if(length(winsfrac) != 1 || winsfrac < 0 || winsfrac > 0.5) {
    stop("Bad value for 'winsfraction'.")
  }
  if(!is.logical(verbose)) {
    stop("Bad value for 'verbose'.")
  }
  set.seed(seed)
  y_name <- paste(attr(Formula::as.Formula(formula), "lhs"))
  y_orig <- data[,y_name]
  y_learn <- data[,y_name]
  x_names <- attr(terms(formula, data = data), "term.labels")
  n <- nrow(data)
  
  #############################
  ## Derive prediction rules ##
  #############################
  
  if(type != "linear") {
    rules <- vector()
    treecount <- 0
    for(i in 1:ntrees) {
      while(length(rules) <= maxrules) {
        # Take subsample of dataset
        treecount <- treecount + 1
        if(sampfrac == 1) { # then bootstrap
          subsample <- sample(1:n, size = n, replace = TRUE)
        } else { # else subsample
          subsample <- sample(1:n, size = round(sampfrac * n), replace = FALSE)
        }
        subsampledata <- data[subsample,]
        subsampledata[,y_name] <- y_learn[subsample]
        # Make sure ctree() can find object specified by weights argument: 
        environment(formula) <- environment()
        # Grow ctree (should later be glmertrees) on subsample:
        tree <- ctree(formula, data = subsampledata, weights = weights[subsample], maxdepth = maxdepth,
                      mtry = mtry)
        # Collect rules from tree:
        rules <- append(rules, unlist(partykit:::.list.rules.party(tree)))
        # Substract predictions from current y:
        y_learn <- y_learn - learnrate * predict(tree, newdata = data)
      }
    }
    nrules <- length(rules)
    if(verbose){
      cat("\nA total of", treecount, "trees were grown, and a total of", nrules, "rules were generated initally.")
    }
    # Keep unique, non-empty rules only:
    rules <- unique(rules[!rules==""])
    if(verbose) {
      cat("\n\nA total of", nrules - length(rules), "rules were duplicates and removed from the initial ensemble.")
    }
    # Create dataframe with 0-1 coded rules:
    rulevars <- data.frame(rule1 = as.numeric(with(data, eval(parse(text = rules[[1]])))))
    for(i in 2:length(rules)) {
      rulevars[,paste("rule", i, sep="")] <- as.numeric(with(data, eval(parse(text = rules[[i]]))))
    }
  }
  if(removeduplicates) {
    duplicates <- duplicated(t(rulevars))
    duplicates.removed <- colnames(rulevars[,duplicates])
    rulevars <- rulevars[,!duplicates]
    rules <- rules[!duplicates]
    if(verbose) {
      cat("\n\nA total of", length(duplicates.removed), "generated rules had support identical to earlier rules and were removed from the initial ensemble ($duplicates.removed shows which, if any).")
    }
  } else {
    duplicates.removed <- NULL
  }
  if(verbose) {
    cat("\n\nAn initial ensemble was succesfully created, consisting of", ncol(rulevars), "rules.")  
  }

  ######################################################
  ## Prepare rules, linear terms and outcome variable ##
  ######################################################
  
  x <- data[,x_names]
  if(type != "rules") {
    # Winsorize variables (section 5 of F&P(2008)):
    if(winsfrac > 0) {
      for(i in 1:ncol(x)) {
        if(is.numeric(x[,i])) { 
          lim <- quantile(x[,i], probs = c(winsfrac, 1-winsfrac))
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
  y <- data[,y_name]
  
  ##################################################
  ## Perform penalized regression on the ensemble ##
  ##################################################
  
  glmnet.fit <- cv.glmnet(as.matrix(x), y, alpha = alpha, standardize = standardize, nfolds = nfolds,
                          type.measure = mod.sel.crit, thres = thres, dfmax = dfmax, weights = weights)
  
  ####################
  ## Return results ##
  ####################
  
  lmin_ind <- which(glmnet.fit$lambda == glmnet.fit$lambda.min)
  l1se_ind <- which(glmnet.fit$lambda == glmnet.fit$lambda.1se)
  if(verbose) {
    cat("\n\nFinal ensemble with minimum cv error: \n  lambda = ", glmnet.fit$lambda[lmin_ind],  
        "\n  number of terms = ", glmnet.fit$nzero[lmin_ind], 
        "\n  mean cv error (se) = ", glmnet.fit$cvm[lmin_ind], 
        " (", glmnet.fit$cvsd[lmin_ind], ")", 
        "\n\nEnsemble with cv error within 1se of minimum: \n  lambda = ", glmnet.fit$lambda[l1se_ind],  
        "\n  number of terms = ", glmnet.fit$nzero[l1se_ind], 
        "\n  mean cv error (se) = ", glmnet.fit$cvm[l1se_ind], 
        " (", glmnet.fit$cvsd[l1se_ind], ")", sep="")
  }
  result <- list(glmnet.fit = glmnet.fit, call = match.call(), weights = weights, 
                 data = data, normalize = normalize, x_scales = x_scales, type = type,
                 rules = data.frame(rule = names(rulevars), 
                                    descriptions = rules), 
                 x_names = x_names, y_name = y_name, baselearn_names = colnames(x), x_sds = x_sds,
                 rulevars = rulevars, removeduplicates = removeduplicates, 
                 duplicates.removed = duplicates.removed)
  class(result) <- "upre"
  return(result)
}





#' Coefficients for the final prediction rule ensemble
#'
#' \code{coef.upre} gets coefficients for each of the prediction rules and linear terms in the ensemble
#' 
#' @param object an object resulting from application of \code{upre()} 
#' @param penalty.par.val character. Should model be selected with lambda giving minimum cv error ("lambda.min"), 
#' or lambda giving cv error that is within 1 standard error of the minimum cv error ("lambda.1se")?
#' @param print logical. Should the coefficients of the base learners with non-zero coeffcients in the
#' final ensemble be printed to the command line?
#' @param ... additional arguments to be passed to \code{\link[glmnet]{coef.glmnet}}.
#' @return returns a dataframe with 3 columns: coefs (coefficients), rule (rule or variable name) and
#' descriptions (<NA> for linear terms, conditions for rules). In the command line, the non zero coefficients
#' are printed
#' @examples
#' \dontrun{
#'  airq.ens <- upre(Ozone ~ ., data=airquality[complete.cases(airquality),])
#'  coef(airq.ens)
#' }
#' @export
coef.upre <- function(object, penalty.par.val = "lambda.min", print = TRUE, ...)
{
  coefs <- as(coef.glmnet(object$glmnet.fit, s = penalty.par.val, ...), Class = "matrix")
  # coefficients for normalized variables should be unnormalized 
  if(object$normalize) {
    coefs[names(object$x_scales),] <- coefs[names(object$x_scales),] * (1/object$x_scales)
  }
  coefs <- data.frame(coefs = round(coefs[,1], digits = 4), rule = rownames(coefs))
  coefs <- merge(coefs, object$rules, all.x=T)
  nonzerocoefs <- coefs[coefs$coefs != 0,]
  if(print) {
    print(nonzerocoefs[order(abs(nonzerocoefs$coefs), decreasing = TRUE),])
  }
  return(coefs)
}





#' Predicted values based on the final prediction rule ensemble
#'
#' \code{predict.upre} generates predictions based on the final ensemble for training, or for new 
#' test observations
#' 
#' @param object an object resulting from application of final()
#' @param newdata optional dataframe of new observations, including all predictor variables used to
#' generate the initial ensemble
#' @param penalty.par.val character. Should model be selected with lambda giving minimum cv error ("lambda.min"), 
#' @param ... additional arguments to be passed (currently not used). 
#' @details When newdata is not provided, training data included in the specified object is used.
#' @examples
#' \dontrun{
#'  airq.ens <- upre(Ozone ~ ., data = airquality[complete.cases(airquality),])
#'  predict(airq.ens, newdata = airquality[complete.cases(airquality),])
#'  predict(airq.ens)
#' }
#' @import Matrix
#' @export
#' 
predict.upre <- function(object, newdata = NULL, penalty.par.val = "lambda.min", ...)
{
  if(is.null(newdata)) {
    newdata <- data.frame(object$data)
  } else {
    newdata <- data.frame(newdata)
  }
  # evaluate all rules with non-zero coefficients for the new dataset:
  coefs <- as(coef.glmnet(object$glmnet.fit, s = penalty.par.val), Class = "matrix")
  nonzerorulenames <- names(coefs[coefs!=0,])[-1] # -1 is for removing the intercept
  nonzerorules <- as.character(object$rules$descriptions[object$rules$rule %in% nonzerorulenames])
  for(i in 1:length(nonzerorules)) {
    newdata[,nonzerorulenames[i]] <- as.numeric(with(data.frame(newdata), eval(parse(text = nonzerorules[i]))))
  }
  # set all rules with zero coefficients to 0:
  zerorulenames <- names(coefs[coefs==0,])
  zerorules <- as.character(object$rules$descriptions[object$rules$rule %in% zerorulenames])
  newdata[,zerorulenames] <- 0
  # linear terms that were normalized before application of glmnet should be normalized prior to
  # applying predict.glmnet:
  if(object$normalize) {
    newdata[,names(object$x_scales)] <- scale(newdata[,names(object$x_scales)], center = FALSE, 
                                              scale = object$x_scales)
  }
  # Make it a sparse matrix, to avoid memory problems:
  newdata <- Matrix::Matrix(as.matrix(newdata[object$baselearn_names]), sparse = TRUE)
  return(predict.cv.glmnet(object$glmnet.fit, newx = newdata, s = penalty.par.val))
}





#' Create partial dependence plot for a single variable
#' 
#' \code{singleplot} generates a partial dependence plot to assess the effect of a single predictor 
#' variable, on the predictions of the ensemble
#' 
#' @param object an object of class \code{upre}
#' @param varname character vector of length one, specifying the variable for which the partial dependence plot
#' should be created.
#' @param penalty.par.val character. Should model be selected with lambda giving minimum cv error ("lambda.min"), 
#' @param nvals optional numeric vector or length one. For how many values of x should the partial dependence
#' plot be created?
#' @details By default, a partial dependence plot will be created for each unique observed value of the 
#' specified predictor variable. When the number of unique observed values is large, this may take a long time
#' to compute. Specifying the nvals argument can substantially reduce computing time. When the nvals argument
#' is supplied, values for the minimum, maximum, and nvals - 2 intermediate values of the predictor 
#' variable will be plotted. Providing the name of a variable that does not appear in the final prediction rule ensemble
#' will result in an error. 
#' @examples 
#' \dontrun{
#'  airq.ens <- upre(Ozone ~ ., data = airquality[complete.cases(airquality),])
#'  singleplot(airq.ens, "Temp")
#' }
#' @export
#' @import graphics akima
singleplot <- function(object, varname, penalty.par.val = "lambda.min", nvals = NULL) 
{
  # preliminaries:
  if(length(varname) != 1) {
    stop("A partial dependence plot should be requested for 1 variable")
  }
  if(!is.character(varname)) {
    stop("Specified varname should be of mode character")
  }
  
  # Generate expanded dataset:
  if(is.null(nvals)) {
    newx <- unique(object$data[,varname])
  } else {
    newx <- seq(min(object$data[,varname]), max(object$data[,varname]), length = nvals)
  }
  exp_dataset <- object$data[rep(row.names(object$data), times = length(newx)),]
  exp_dataset[,varname] <- rep(newx, each = nrow(object$data))

  # get predictions:
  exp_dataset$predy <- predict.upre(object, newdata = exp_dataset)

  # create 2D plot:  
  plot(aggregate(exp_dataset$predy, by = exp_dataset[varname], data = exp_dataset, FUN = mean), 
       type = "l", ylab = "predicted y", xlab = varname, main = paste("partial dependence on", varname))
  
  # To be implemented:
  # qntl = trimming factor for plotting numeric variables. Plots are shown for variable values in the range [quantile (qntl) - quantile(1-qntl)]. (Ignored for categorical variables (factors).)
  # nval = maximum number of abscissa evaluation points for numeric variables. (Ignored for categorical variables (factors).)
  # nav = maximum number of observations used for averaging calculations. (larger values provide higher accuracy with a diminishing return; computation grows linearly with nav)
  # catvals = vector of names for values (levels) of categorical variable (factor). (Ignored for numeric variables or length(vars) > 1)
  # samescale = plot vertical scaling flag .
  # samescale = TRUE / FALSE => do/don't require same vertical scale for all plots.
  # horiz = plot orientation flag for categorical variable barplots
  # horiz = T/F => do/don't plot bars horizontally
  # las = label orientation flag for categorical variable plots (horiz = F, only)
  # las = 1 => horizontal orientation of value (level) names stored in catvals (if present)
  # las = 2 => vertical orientation of value (level) names stored in catvals (if present)
  # cex.names = expansion factor for axis names (bar labels) for categorical variable barplots
  # col = color of barplot for categorical variables
  # denqnt = quantile for data density tick marks along upper plot boundary  for numeric variables ( < 1)
  # denqnt <= 0 => no data density tick marks displayed 
}





#' Create partial dependence plot for a pair of predictor variables
#' 
#' \code{pairplot} generates a partial dependence plot to assess the effects of a pair of predictor 
#' variables, on the predictions of the ensemble
#' 
#' @param object an object of class \code{upre}
#' @param varnames character vector of length two
#' @param penalty.par.val character. Should model be selected with lambda giving minimum cv error ("lambda.min"), 
#' @param phi numeric. See \code{persp()} documentation.
#' @param theta numeric. See \code{persp()} documentation.
#' @param col character. Optional color to be used for surface in 3D plot.
#' @param nvals optional numeric vector or length  two. For how many values of x1 and x2 should 
#' partial dependence be plotted?
#' @details By default, partial dependence will be plotted for each combination of unique observed values 
#' of the specified predictor variables. When the number of unique observed values is large, this may 
#' take a long time to compute. Specifying the nvals argument can substantially reduce computing time. 
#' When the nvals argument is supplied, values for the minimum, maximum, and nvals - 2 intermediate values 
#' of the predictor variable will be plotted.
#' @details Providing the names of two variables that do not appear in the final prediction rule ensemble
#' will result in an error. 
#' @examples 
#' \dontrun{
#'  airq.ens <- upre(Ozone ~ ., data = airquality[complete.cases(airquality),])
#'  pairplot(airq.ens, c("Temp", "Wind"))
#' }
#' @export
#' @import graphics akima
pairplot <- function(object, varnames, penalty.par.val = "lambda.min", phi = 45, theta = 315, 
                       col = "cyan", nvals = NULL) 
{
  # preliminaries:
  if(length(varnames) != 2) {
    stop("Partial dependence should be requested for 2 variables")
   }
  if(!is.character(varnames)) {
    stop("Specified varname should be of mode character")
  }

  # generate expanded dataset: 
  if(is.null(nvals)){
    newx1 <- unique(object$data[,varnames[1]])
    newx2 <- unique(object$data[,varnames[2]])
  } else {
    newx1 <- seq(min(object$data[,varnames[1]]), max(object$data[,varnames[1]]), length = nvals[1])
    newx2 <- seq(min(object$data[,varnames[2]]), max(object$data[,varnames[2]]), length = nvals[2])
  }
  nobs1 <- length(newx1)
  nobs2 <- length(newx2)
  nobs <- nobs1*nobs2
  exp_dataset <- object$data[rep(row.names(object$data), times = nobs),]
  exp_dataset[,varnames[1]] <- rep(newx1, each = nrow(object$data)*nobs2)
  exp_dataset[,varnames[2]] <- rep(rep(newx2, each = nrow(object$data)), times = nobs1)
  
  # get predictions:  
  pred_vals <- predict.upre(object, newdata = exp_dataset, penalty.par.val = penalty.par.val)

  # create 3D plot:
  if(is.null(nvals)) nvals <- 3
  xyz <- akima::interp(exp_dataset[,varnames[1]], exp_dataset[,varnames[2]], pred_vals, duplicate = "mean")
  persp(xyz, xlab = varnames[1], ylab = varnames[2], zlab = "predicted y", phi = phi, theta = theta,
        col = col, ticktype = "detailed", nticks = max(nvals))
  
  # to be implemented:
  # type = flag for type of plot when both var1 and var2 are numeric
  # type = "image" => heat map plot
  # type = "persp" => perspective mesh plot
  # type = "contour" => contour plot
  # chgvars = flag for changing plotting relationship when both var1 and var2 are categorical (factors)
  # chgvars = FALSE => plot the partial dependence on the variable (factor) with the most values (levels), for each of the  respective values (levels) of the other variable (factor)
  # chgvars = TRUE => reverse this relationship
  # qntl = trimming factor for plotting numeric variables. Plots are shown for variable values in the range [quantile (qntl) - quantile(1-qntl)]. (Ignored for categorical variables (factors).)
  # nval = maximum number of evaluation points for numeric variables. (Ignored for categorical variables).
  # nav = maximum number of observations used for averaging calculations. (larger values provide higher accuracy with a diminishing return; computation grows linearly with nav)
  # vals1 = vector of names for values (levels) of var1 if it is categorical (factor). (Ignored if var1 is numeric)
  # vals2 = vector of names for values (levels) of var2 if it is categorical (factor). (Ignored if var2 is numeric) 
  # horiz = plot orientation for categorical variable barplots
  # horiz = T/F => do/don't plot bars horizontally
  # las = label orientation flag for categorical variable plots (horiz = F, only)
  # las =1 => horizontal orientation of value (level) names stored in vals1 and/or vals2 (if present).
  # las =2 => vertical orientation of value (level) names stored in vals1 and/or vals2 (if present).
  # cex.names = expansion factor for axis names (bar labels)  for categorical variable barplots 
}





#' Calculate importances of base learners (rules and linear terms) and input variables
#' 
#' \code{importance} calculates importances for rules, linear terms and input variables in the ensemble,
#' and provides a bar plot of variable importances
#' 
#' @param object an object of class \code{upre}
#' @param plot logical. Should variable importances be plotted?
#' @param ylab character. Only used when \code{plot = TRUE}. PLotting label for y-axis.
#' @param ... further arguments to be passed to \code{barplot}
#' @return A list with two dataframes: $baseimps, giving the importance for each baselearner (not) in the
#' ensemble, and $varimps, giving the importance for each predictor variable (not) in the ensemble
#' @examples 
#' \dontrun{
#'  airq.ens <- upre(Ozone ~ ., data=airquality[complete.cases(airquality),])
#'  importance(airq.ens)
#' }
#' @export
importance <- function(object, plot = TRUE, ylab = "Importance", main = "Variable importances", ...) 
{
  # Get sds for every linear function and/or rule:
  if(object$type == "both") {
    sds <- c(0, apply(object$rulevars, 2, sd, na.rm = TRUE), object$x_sds)    
  }
  if(object$type == "linear") { 
    sds <- c(0, object$x_sds)
  }
  if(object$type == "rules") { 
    sds <- c(0, apply(object$rulevars, 2, sd, na.rm = TRUE))
  }
  names(sds)[1] <- "(Intercept)"
  sds <- sds[order(names(sds))]
  coefs <- coef.upre(object, print = FALSE)
  coefs$descriptions <- as.character(coefs$descriptions)
  coefs$descriptions[is.na(coefs$descriptions)] <- 
    paste(as.character(coefs$rule)[is.na(coefs$descriptions)], " ", sep = "")
  coefs <- coefs[order(coefs$rule),]
  if(all(names(sds) != coefs$rule)) {
    stop("There seems to be a problem with the ordering or size of the coefficient and sd vectors. Importances cannot be calculated.")
  }
  # Calculate importances (absolute value of the coefficient * stand dev, see F&P section 6):
  baseimps <- data.frame(coefs, sds = sds, imps = abs(coefs$coefs)*sds)
  # Calculate a variable counting of the number of terms in each rule:
  baseimps$nterms <- NA
  for(i in 1:nrow(baseimps)) {
    # If there is no "&" in the rule description, there is only 1 term/variable in the base learner: 
    if(gregexpr("&", baseimps$descriptions)[[i]][1] == -1) {
      baseimps$nterms[i] <- 1 
    } else {
      baseimps$nterms[i] <- length(gregexpr("&", baseimps$descriptions)[[i]]) + 1
    }
  }
  varimps <- data.frame(varnames = object$x_names, imps = 0)
  # For every input variable:
  for(i in 1:nrow(varimps)) {
    # For every baselearner:
    for(j in 1:nrow(baseimps)) {
      # if the variable name appears in the rule:
      if(gregexpr(paste(varimps$varnames[i], " ", sep = ""), baseimps$descriptions[j])[[1]][1] != -1) {
        # then count the number of times it appears in the rule:
        n_occ <- length(
          gregexpr(paste(varimps$varnames[i], " ", sep = ""), baseimps$descriptions[j])[[1]]
        )
        # and add it to the importance of the variable:
        varimps$imps[i] <- varimps$imps[i] + (n_occ * baseimps$imps[j] / baseimps$nterms[j])
      }
    }
  }
  varimps <- varimps[order(varimps$imps, decreasing = TRUE),]
  baseimps <- baseimps[order(baseimps$imps, decreasing = TRUE),]
  barplot(height = varimps$imps, names.arg = varimps$varnames, ylab = ylab, main = main)
  return(list(varimps = varimps, baseimps = baseimps[, c("imps", "coefs", "descriptions")]))
}





#' Compute boostrapped null interaction models for reference distributions of interaction test statistics
#' 
#' \code{bsnullinteract} calculates null interaction models on bootstrapped datasets, for deriving a 
#' reference distribution of the test statistic calculated with \code{interact}
#' 
#' @param object an object of class \code{upre}
#' @param nsamp the number of bootstrapped null interaction models to be derived
#' @param seed numeric. Random seed to be used in deriving the final ensemble (for reproducability). 
#' Defaults to 42.
#' @return A list of null interactrion models
#' @examples 
#' \dontrun{
#'  airq.ens <- upre(Ozone ~ ., data=airquality[complete.cases(airquality),])
#'  nullmods <- bsnullinteract(airq.ens)
#' }
#' @details Can be computationally intensive.
#' @export
bsnullinteract <- function(object, nsamp = 10, seed = 42) {
  # preliminaries:
  set.seed(seed)
  # create a call for creating the bootstrapped null model (only dataset and maxdepth should be changed):
  bsnullmodcall <- object$call
  bsnullmodcall$maxdepth <- 1
  bsnullmodcall$verbose <- FALSE  
  # create a call for a model allowing for interactions, grown on the bootstrapped datasets without 
  # interactions
  bsintmodcall <- bsnullmodcall
  if(is.null(object$call$maxdepth)) {
    bsintmodcall$maxdepth <- 3
  } else {
    bsintmodcall$maxdepth <- object$call$maxdepth
  }
  # compute boostrapped null dataset (i.e., dataset with no interactions):
  bs.ens <- list()
  cat("This may take a while. Computing null model ")
  for(i in 1:nsamp) {
    cat(i, "of", nsamp, "... ")
    # step 1: Take bootstrap sample {x_ip, y_ip}:
    bsdataset <- object$data[sample(1:nrow(object$data), nrow(object$data), replace = TRUE),]
    # step 2: Build F_null, a model involving main effects only using {x_ip, y_ip}:
    bsnullmodcall$data <- bsdataset
    bs.ens.null <- eval(match.call(upre, call = bsnullmodcall))
    # step 3: Calculate residuals from predictions y^hat_ip using x_ip and F_null:
    yhatip <- bsdataset[object$y_name] - predict.upre(bs.ens.null)
    # step 4: Calculate predictions for original x, using F_null:
    fipx <- predict.upre(bs.ens.null, newdata = object$data)
    # step 5: Calculate ybar, by adding residuals from step 3 to predictions from step 4:
    bsdataset[,object$y_name] <- yhatip + fipx
    # step 6: Build a model using (x,ybar), using the same procedure as was originally applied to (x,y):
    bsintmodcall$data <- bsdataset
    bs.ens[[i]] <- eval(match.call(upre, call = bsintmodcall))
  }
  cat("Done!")
  return(bs.ens)
}





# Internal function for calculating H statistic (section 8.1, equation 45):
interact <- function(object, varname, k = 10) {
  # Eq. 39 states that the partial dependence of a function F on a subset of variables x_s is the 
  # expected value of F(x_s), averaged over all values of x_/s 
  # Eq. 40 states that this value can be estimated from the data by summing, for each value of x_s,
  # the predicted values for each of the N data values of x_/s, and dividing this sum by N.
  # Eq. 45 shows how the test statistic H can be computed for this:
  #
  # Calculate the predicted value F(x) of the full model for each observation:
  predsx <- predict(object, newdata = object$data)
  # Calculate the expected value of F_j(x_j), over all observed values x_/j,
  # and the expected value of F_/j(x_/j), over all observed values x_j:
  exp_dataset <- object$data[rep(row.names(object$data), times = nrow(object$data)),]
  exp_dataset[,varname] <- rep(object$data[,varname], each = nrow(object$data))
  ids <- caret::createFolds(1:nrow(exp_dataset), k = k)
  for(i in 1:k) {  
    cat(".")
    exp_dataset[ids[[i]], "yhat"] <- predict.upre(object, newdata = exp_dataset[ids[[i]],])
  }
  # expected value of F_j(x_j), over all observed values x_/j:
  exp_dataset$i_xj <- rep(1:nrow(object$data), each = nrow(object$data))
  preds_xj <- aggregate(yhat ~ i_xj, data = exp_dataset, FUN = mean)$V1
  # expected value of F_/j(x_/j), over all observed values x_j:
  exp_dataset$i_xnotj <-  rep(1:nrow(object$data), times = nrow(object$data))
  preds_xnotj <- aggregate(yhat ~ i_xnotj, data = exp_dataset, FUN = mean)$V1
  return(sum((predsx - preds_xj - preds_xnotj)^2) / sum(scale(predsx, scale = FALSE)^2))
}





#' Calculate interaction test statistic for user-specified variable 
#' 
#' \code{interaction} calculates a statistic for testing whether a user-supplied variable interacts with
#' any other variable in the ensemble.
#' 
#' @param object an object of class \code{upre}
#' @param varname character. Variable for which interation test statistic should be calculated.
#' @param k integer. Calculating interaction test statistics is a computationally intensive, so  
#' calculations are split up in several parts to prevent memory allocation errors. If a memory allocation
#' error still occurs, increase k.
#' @param nullmods object with bootstrapped null interaction models, resulting from application of 
#' \code{bsnullinteract}.
#' @examples 
#' \dontrun{
#'  airq.ens <- upre(Ozone ~ ., data=airquality[complete.cases(airquality),])
#'  interaction(airq.ens, "Temp")
#' }
#' @details Can be computationally intensive, especially when nullmods is specified.
#' @return If nullmods is not specified: the test statistic of the interaction strength
#' If nullmods is specified:
#'  $H = the test statistic of the interaction strength
#'  $nullH = a vector of test statistics of the interaction strength for each of the bootstrapped null
#'  interaction models
#' @import caret
#' @export
interaction <- function(object, varname, nullmods = NULL, k = 10)
{
  cat("Please be patient, this may take a while (", k*(length(nullmods)+1), "dots ). ")
  # Calculate H_j for the original dataset:
  H <- interact(object = object, varname = varname, k = k)
  # Calculate mean and sd of H_j for the bootstrapped null models  
  if(!is.null(nullmods)) {
    nullint <- list()
    for(i in 1:length(nullmods)) {
      nullint[[i]] <- interact(object = nullmods[[i]], varname = varname, k = k)  
    }
    nullH <- unlist(nullint)
  }
  cat("\n")
  if(is.null(nullmods)) {return(H)}
  if(!is.null(nullmods)) {return(list(H, nullH))}  
  # Add check whether specified x_j appears in the final ensemble at all. If not, function should stop.
} 