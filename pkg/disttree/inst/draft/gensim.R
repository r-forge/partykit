# FIX ME:
# error in bamlss (seedconst7): "Error in smooth.construct.tp.smooth.spec(object, dk$data, dk$knots) : \n  object 'C_construct_tprs' not found\n"
# error in gamboostLSS (seedconst723, fix.mu): "Error in names(pr) <- nm : \n  'names' attribute [300] must be the same length as the vector [10]\n"
# error in gamlss (seedconst723, covsep, 3 out of 10): ""Error in while (RATIO > tol & nit < maxit) { : \n  missing value where TRUE/FALSE needed\n" 
# error in gamboostLSS: "Error in as.families(fname = \"NOlc\") : \n  ‘fname’ specifies no valid gamlss family\n"
# error in gamlss: "Error in eval(parse(text = object$family[[1]])) : object 'NOlc' not found\n"
# error in bamlss: "Error in smooth.construct.tp.smooth.spec(object, dk$data, dk$knots) : \n  object 'C_construct_tprs' not found\n"

gensim <- function(seedconst = 7, nrep = 100, ntree = 100,
                   nsteps = 9, stepsize = 1, kappa.start = 1,
                   formula = y~x1+x2, 
                   nobs = 400, testnobs = 200L,
                   tree_minsplit = 14, tree_mincrit = 0.95,
                   forest_minsplit = 7, forest_mincrit = 0, 
                   type.tree = "ctree",
                   onecov = TRUE,
                   censNO = TRUE,
                   cov.sep = FALSE,
                   fix.mu = FALSE,
                   fix.sigma = FALSE,
                   mu.sigma.interaction = FALSE,
                   gamboost_cvr = FALSE,
                   eval_disttree = TRUE,
                   eval_distforest = TRUE,
                   eval_gamlss = TRUE,
                   eval_bamlss = FALSE,
                   eval_gamboostLSS = FALSE,
                   eval_randomForest = FALSE,
                   eval_cforest = FALSE)
{
  
  cl <- match.call()
  
  ### preliminaries
  library("disttree")
  library("gamlss")
  library("lattice")
  library("crch")
  library("latex2exp")
  library("parallel")
  library("gamlss.cens")
  library("gamboostLSS")
  library("bamlss")
  library("randomForest")
  library("gamboostLSS")
  library("Formula")  # FIX ME: should be imported in disttree
  library("survival")
  
  ## define distribution list:
  # dist_list_normal
  {
    
    dist_list_normal <- list()
    
    parnames <- c("mu", "sigma")
    etanames <- c("mu", "log(sigma)")
    
    
    ddist <-  function(y, eta, log = TRUE, weights = NULL, sum = FALSE) {     
      
      val <- -1/2 * (log(2*pi) + 2*eta[2] + exp(log((y-eta[1])^2) - 2*eta[2]))
      if(!log) val <- exp(val)
      
      # par <- c(eta[1], exp(eta[2]))
      # val <- dnorm(y, mean = par[1], sd = par[2], log = log)
      
      if(sum) {
        if(is.null(weights)) weights <- rep.int(1, length(y))
        val <- sum(weights * val, na.rm = TRUE)
      }
      return(val)
    }
    
    
    sdist <- function(y, eta, weights = NULL, sum = FALSE) {   
      
      score <- cbind(exp(-2*eta[2]) * (y-eta[1]), 
                     -1 + exp(-2*eta[2] + log((y-eta[1])^2)))
      
      # par <- c(eta[1], exp(eta[2])) 
      # score <- cbind(1/par[2]^2 * (y-par[1]), 
      #                (-1/par[2] + ((y - par[1])^2)/(par[2]^3)) * exp(eta[2]))
      
      score <- as.matrix(score)
      colnames(score) <- etanames
      if(sum) {
        if(is.null(weights)) weights <- rep.int(1, length(y))
        # if score == Inf replace score with 1.7e308 because Inf*0 would lead to NaN -> gradient is NaN
        score[score==Inf] = 1.7e308
        score <- colSums(weights * score, na.rm = TRUE)
      }
      return(score)
    }
    
    
    hdist <- function(y, eta, weights = NULL) {    
      ny <- length(y)
      if(is.null(weights)) weights <- rep.int(1, ny)
      
      d2ld.etamu2 <- sum(weights * rep.int(-exp(-2*eta[2]), ny))
      d2ld.etamu.d.etasigma <- sum(weights * (-2)*(y-eta[1]) * exp(-2*eta[2]), na.rm = TRUE)          # should be 0 for exact parameters (here: observed hess)
      d2ld.etasigma2 <- sum(weights * (-2)*exp(log((y-eta[1])^2) - 2*eta[2]), na.rm = TRUE)    
      
      # par <- c(eta[1], exp(eta[2]))                           
      # d2ld.etamu2 <- sum(weights * rep.int(-1/par[2]^2, ny))
      # d2ld.etamu.d.etasigma <- sum(weights * (-2)*(y-par[1])/par[2]^2), na.rm = TRUE)          # should be 0 for exact parameters (here: observed hess)
      # d2ld.etasigma2 <- sum(weights * (-2)*(y-par[1])^2/par[2]^2, na.rm = TRUE)         
      
      hess <- matrix(c(d2ld.etamu2, d2ld.etamu.d.etasigma, d2ld.etamu.d.etasigma, d2ld.etasigma2), nrow = 2)
      colnames(hess) <- rownames(hess) <-  etanames
      
      return(hess)
    }
    
    
    ## additional functions pdist, qdist, rdist
    pdist <- pnorm
    qdist <- qnorm
    rdist <- rnorm  
    
    
    link <- c("identity", "log")
    
    linkfun <- function(par) {
      eta <- c(par[1], log(par[2]))
      names(eta) <- etanames
      return(eta)
    }
    
    
    linkinv <- function(eta) {
      par <- c(eta[1], exp(eta[2]))
      names(par) <- parnames
      return(par)
    }
    
    
    linkinvdr <- function(eta) {
      dpardeta <- c(1, exp(eta[2]))
      names(dpardeta) <- parnames
      return(dpardeta)
    }
    
    
    startfun <- function(y, weights = NULL){
      if(is.null(weights)) {
        mu <- mean(y)
        sigma <- sqrt(1/length(y) * sum((y - mu)^2))
      } else {
        mu <- weighted.mean(y, weights)
        sigma <- sqrt(1/sum(weights) * sum(weights * (y - mu)^2))
      }
      starteta <- c(mu, log(sigma))
      names(starteta) <- etanames
      return(starteta)
    }
    
    mle <- TRUE
    
    dist_list_normal <- list(family.name = "Normal Distribution",
                             ddist = ddist, 
                             sdist = sdist, 
                             hdist = hdist,
                             pdist = pdist,
                             qdist = qdist,
                             rdist = rdist,
                             link = link, 
                             linkfun = linkfun, 
                             linkinv = linkinv, 
                             linkinvdr = linkinvdr,
                             startfun = startfun,
                             mle = mle,
                             gamlssobj = FALSE,
                             censored = FALSE
    )
  }
  
  
  # dist_list_cens_normal
  {
    
    dist_list_cens_normal <- list()
    
    parnames <- c("mu", "sigma")
    etanames <- c("mu", "log(sigma)")
    
    ddist <-  function(y, eta, log = TRUE, weights = NULL, sum = FALSE, left = 0, right = Inf) {     
      par <- c(eta[1], exp(eta[2]))
      val <- crch::dcnorm(x = y, mean = par[1], sd = par[2], left = left, right = right, log = log)
      if(sum) {
        if(is.null(weights)) weights <- if(is.matrix(y)) rep.int(1, dim(y)[1]) else rep.int(1, length(y))
        val <- sum(weights * val, na.rm = TRUE)
      }
      return(val)
    }
    
    
    sdist <- function(y, eta, weights = NULL, sum = FALSE, left = 0, right = Inf) {   
      par <- c(eta[1], exp(eta[2]))
      # y[y==0] <- 1e-323
      
      score_m <- crch:::scnorm(x = y, mean = par[1], sd = par[2], which = "mu", left = left, right = right)
      score_s <- crch:::scnorm(x = y, mean = par[1], sd = par[2], which = "sigma", left = left, right = right) * exp(eta[2]) # inner derivation exp(eta[2])
      score <- cbind(score_m, score_s)
      score <- as.matrix(score)
      colnames(score) <- etanames
      if(sum) {
        if(is.null(weights)) weights <- rep.int(1, length(y)[1])
        # if score == Inf replace score with 1.7e308 because Inf*0 would lead to NaN (0 in weights)
        score[score==Inf] = 1.7e308
        score <- colSums(weights * score, na.rm = TRUE)
        #if(any(is.nan(score))) print(c(eta, "y", y))
      }
      return(score)
    }
    
    
    hdist <- function(y, eta, weights = NULL, left = 0, right = Inf) {    
      ny <- length(y)
      if(is.null(weights)) weights <- rep.int(1, ny)
      
      par <- c(eta[1], exp(eta[2]))                           
      # y[y==0] <- 1e-323
      
      d2mu <- crch:::hcnorm(x = y, mean = par[1], sd = par[2], which = "mu", left = left, right = right)
      d2sigma <- crch:::hcnorm(x = y, mean = par[1], sd = par[2], which = "sigma", left = left, right = right)
      dmudsigma <- crch:::hcnorm(x = y, mean = par[1], sd = par[2], which = "mu.sigma", left = left, right = right) # FIX: order?
      dsigmadmu <- crch:::hcnorm(x = y, mean = par[1], sd = par[2], which = "sigma.mu", left = left, right = right) # FIX: order?
      dsigma <- crch:::scnorm(x = y, mean = par[1], sd = par[2], which = "sigma", left = left, right = right)
      
      d2ld.etamu2 <- sum(weights * d2mu, na.rm = TRUE)
      d2ld.etamu.d.etasigma <- sum(weights * dmudsigma * par[2], na.rm = TRUE)
      d2ld.etasigma.d.etamu <- sum(weights * dsigmadmu * par[2], na.rm = TRUE)
      d2ld.etasigma2 <- sum(weights * (d2sigma * exp(2*eta[2]) + dsigma * par[2]), na.rm = TRUE)         
      
      hess <- matrix(c(d2ld.etamu2, d2ld.etamu.d.etasigma, d2ld.etasigma.d.etamu, d2ld.etasigma2), nrow = 2)
      colnames(hess) <- rownames(hess) <-  etanames
      
      return(hess)
    }
    
    
    ## additional functions pdist, qdist, rdist
    pdist <- function(q, eta, lower.tail = TRUE, log.p = FALSE) crch:::pcnorm(q, mean = eta[1], sd = eta[2], 
                                                                              lower.tail = lower.tail, log.p = log.p, 
                                                                              left = left, right = right)
    qdist <- function(p, eta, lower.tail = TRUE, log.p = FALSE) crch:::qcnorm(p, mean = eta[1], sd = eta[2], 
                                                                              lower.tail = lower.tail, log.p = log.p, 
                                                                              left = left, right = right)
    rdist <- function(n, eta) crch:::rcnorm(n, mean = eta[1], sd = eta[2], left = left, right = right)
    
    
    link <- c("identity", "log")
    
    linkfun <- function(par) {
      eta <- c(par[1], log(par[2]))
      names(eta) <- etanames
      return(eta)
    }
    
    
    linkinv <- function(eta) {
      par <- c(eta[1], exp(eta[2]))
      names(par) <- parnames
      return(par)
    }
    
    
    linkinvdr <- function(eta) {
      dpardeta <- c(1, exp(eta[2]))
      names(dpardeta) <- parnames
      return(dpardeta)
    }
    
    
    startfun <- function(y, weights = NULL){
      yc <- pmax(0,y)  # optional ?
      if(is.null(weights)) {
        mu <- mean(yc)
        sigma <- sqrt(1/length(yc) * sum((yc - mu)^2))
      } else {
        mu <- weighted.mean(yc, weights)
        sigma <- sqrt(1/sum(weights) * sum(weights * (yc - mu)^2))
      }
      starteta <- c(mu, log(sigma))
      names(starteta) <- etanames
      return(starteta)
    }
    
    mle <- FALSE
    
    dist_list_cens_normal <- list(family.name = "censored Normal Distribution",
                                  ddist = ddist, 
                                  sdist = sdist, 
                                  hdist = hdist,
                                  pdist = pdist,
                                  qdist = qdist,
                                  rdist = rdist,
                                  link = link, 
                                  linkfun = linkfun, 
                                  linkinv = linkinv, 
                                  linkinvdr = linkinvdr,
                                  startfun = startfun,
                                  mle = mle,
                                  gamlssobj = FALSE,
                                  censored = TRUE
    )
  }
  
  
  
  #if(censNO != family$censored) stop("if censNO = TRUE the family has to be censored and vice versa")
  family <- if(censNO) dist_list_cens_normal else dist_list_normal
  
  ######## data generating process based on a predefined tree or a parameter function
  # if tree: 2 splits => 3 terminal nodes
  # 10 split variables (x1, ..., x10) are availabe (5 numeric, 2 binary, 3 categorical)
  #
  # input: n            ..... nr of observations
  #        family       ..... distribution of the generated observations
  #        split.matrix ..... indices of the variables used for the splits together with the corresponding split points
  #        parm.matrix  ..... set of distribution parameters for each subgroup
  #        fun          ..... parameter function
  #
  # output: data.frame with generated observations y (generated seperatly in each subgroups with given distributions parameters),
  #                         the given split variables x1, ..., x10 for each observation and
  #                         the index of the subgroup for each observation or
  #                         the distribution parameter for each observation
  
  # FIX: until now only complete lists can be handed over (dist_list...)
  dgp <- function(n, family = dist_list_normal, 
                  fun = NULL,
                  split.matrix = matrix(nrow = 2, ncol = 2), 
                  par.matrix = matrix(nrow = 3, ncol = 2), 
                  round.sp = 3)
  {
    
    # generating the possible split variables
    x1 <- runif(n,-0.4,1)
    x2 <- runif(n,-10,10)
    x3 <- runif(n,0,100)
    x4 <- x1 + rnorm(n, sd = 0.1)
    x5 <- x1 + rnorm(n, sd = 0.3)
    x6 <- rbinom(n,1,0.5)
    x7 <- rbinom(n,1,0.5)
    x8 <- sample(1:4, n, replace = TRUE) 
    x9 <- sample(1:5, n, replace = TRUE) 
    x10 <- sample(1:7, n, replace = TRUE)
    #x11 <- runif(n,-0.5,1)
    x <- cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)
    # reduce nr of possible split points by rounding values of split variables
    x <- round(x, digits = round.sp)
    
    y <- vector(mode = "numeric", length = n)
    
    # getting the random function
    if(is.function(family)) family <- family()
    if(inherits(family, "gamlss.family")) {
      rfun <- get(paste0("r",family$family[[1]]))
    } else {
      if(family$family.name == "Normal Distribution") rfun <- rNO
      if(family$family.name == "censored Normal Distribution") rfun <- rNO
      if(family$family.name == "Poisson Distribution") rfun <- rPO
      if(family$family.name == "Weibull Distribution") rfun <- survival:::rsurvreg 
    }
    
    
    if(is.null(fun)) {
      index <- vector(mode = "numeric", length = n)
      sv = split.matrix[,1]    # index of split variables
      sp = split.matrix[,2]    # split points
      par = par.matrix
      
      # splitting based on the given split variables and split points and
      # generating the observations seperatly for each subgroup (with the corresponding given distribution parameters)
      for(k in 1:n){
        if(x[k, sv[1]] <= sp[1]){
          y[k] <- do.call(rfun, as.list(c(1,par[1,])))
          index[k] <- 2L
        } else {
          if(x[k, sv[2]] <= sp[2]){
            y[k] <- do.call(rfun, as.list(c(1,par[2,]))) 
            index[k] <- 4L
          } else {
            y[k] <- do.call(rfun, as.list(c(1,par[3,])))
            index[k] <- 5L
          }
        }
      }
      d <- as.data.frame(cbind(y, x, index))
      colnames(d) <- c("y", paste0("x", c(1:10)), "index")
      
    } else {
      
      dpar <- fun(x)
      
      if(!(family$family.name == "Poisson Distribution")){
        y <- rfun(n, dpar[,1], dpar[,2])
        d <- as.data.frame(cbind(y, x, dpar))
        colnames(d) <- c("y", paste0("x", c(1:10)), "mu", "sigma")
      } else {
        y <- rfun(n, dpar)
        d <- as.data.frame(cbind(y, x, dpar))
        colnames(d) <- c("y", paste0("x", c(1:10)), "lambda")
      }
    }
    
    
    if(family$family.name == "censored Normal Distribution") {
      #d$ystar <- d$y
      d$y <- pmax(d$y, 0)
    }
    
    return(d)
  }
  
  
  
  
  
  ## function to estimate standard deviation of cforest for a new observation
  cf_getsd <- function(cf, newdata = NULL){
    
    # get IDs of predicted nodes for the learning data cfdata (does not have to be handed over as cf was learned on cfdata)
    pred.node.learn <- predict(cf, type = "node")
    
    # get IDs of predicted nodes for the new observations
    if(is.null(newdata)) {
      pred.node.new <- pred.node.learn
      newdata <- cf$data
    } else {
      pred.node.new <- predict(cf, newdata = newdata, type = "node")
    }
    
    sdnew <- numeric(length = NROW(newdata))
    
    for(i in 1:NROW(newdata)){
      for(t in 1:ntree){
        nodedata <- cf$data[(pred.node.learn[[t]] == pred.node.new[[t]][i]),]
        sdnew[i] <- sd(nodedata[,paste(cf$terms[[2]])])
      }
    }
    
    return(sdnew)   
  }
  
  
  
  ## function to estimate standard deviation of randomForest for a new observation
  # (in randomForest the argument 'keep.inbag' must be set to TRUE)
  rf_getsd <- function(rf, newobs, rfdata){
    
    # get predictions for the new observations from all trees
    pred.newobs <- predict(rf, predict.all = TRUE, newdata = newobs)
    
    # vector where the standard deviations from all trees are stored
    sd_trees <- numeric(length = rf$ntree)
    
    # loop over all trees of the forest
    for(i in 1:rf$ntree){
      
      # get data used to build this tree
      obsid <- rep.int(c(1:NROW(rfdata)), as.vector(rf$inbag[,i]))
      obs_tree <- rfdata[obsid,]
      rownames(obs_tree) <- c(1:NROW(obs_tree))
      # get predictions for this data from this tree
      pred.obs_tree <- predict(rf, newdata = obs_tree, predict.all = TRUE)$individual[,i]
      
      # get prediction for the new observation from this tree
      pred.newobs_tree <- pred.newobs$individual[,i]
      
      # get part of the data that ends up in the same terminal node (has the same prediction)
      obs_node <- obs_tree[pred.obs_tree == pred.newobs_tree,]
      
      sd_trees[i] <- sd(obs_node$y)
    }
    
    # average of sd over all trees
    sd_newobs <- mean(sd_trees, na.rm = TRUE)
    
    return(sd_newobs)   
  }
  
  
  
  
  
  ### set up function to get logLik from disttree object for new data 
  #dtll.newdata <- function(object, newdata){
  #  ll <- 0
  #  nd <- newdata[,-c(1,ncol(newdata)-1, ncol(newdata))]
  #  if(is.null(coef(object))){
  #    readline(prompt="Press [enter] to continue")
  #    print(object$info$call)
  #    print(object$data)
  #    return(NULL)
  #  } else {
  #    pred.par <- predict(object, newdata = nd, type = "parameter")
  #    for(i in 1:(nrow(newdata))){
  #      eta <- as.numeric(object$info$family$linkfun(pred.par[i,]))
  #      ll <- ll + object$info$family$ddist(newdata[i,paste(object$info$formula[[2]])], eta = eta, log=TRUE)
  #    }
  #    if(is.na(ll)) print("dt.ll = NA")
  #    return(ll)
  #  }
  #}
  
  
  
  ### set up function to get logLik from gamlss object for new data and normal distribution
  gll.newdata <- function(object, newdata, data){
    if(is.null(object)) {
      warning("error in gamlss")
      return(NA)
    }
    ll <- 0
    nd <- newdata[,-c(1,ncol(newdata)-1, ncol(newdata))]
    pred.par <- cbind(predict(object, newdata = nd, what = "mu", type = "response", data = data),
                      predict(object, newdata = nd, what = "sigma", type = "response", data = data))
    #np <- ncol(pred.par)
    distfun <- if(censNO) {
      function(y, mean, sd) crch::dcnorm(x = y, mean = mean, sd = sd, left = 0, right = Inf, log = TRUE)
    } else {
      function(y, mean, sd) dnorm(x = y, mean = mean, sd = sd, log = TRUE)
    }
    for(i in 1:(nrow(newdata))){
      par <- pred.par[i,]
      ll <- ll + distfun(newdata[i,1], mean = par[1], sd = par[2])
    }
    if(is.na(ll)) print("g.ll = NA")
    if(is.null(ll)) print("g.ll = NULL")
    return(ll)
  }
  
  
  
  ### set up function to get logLik from bamlss object for new data and normal distribution
  bll.newdata <- function(object, newdata){
    if(is.null(object)) {
      warning("error in bamlss")
      return(NA)
    }
    ll <- 0
    nd <- newdata[,-c(1,ncol(newdata)-1, ncol(newdata))]
    pred.par <- predict(object, newdata = nd, what = "parameters", type = "parameter")
    pred.par <- cbind(pred.par$mu, pred.par$sigma)
    #np <- ncol(pred.par)
    distfun <- if(censNO) {
      function(y, mean, sd) crch::dcnorm(x = y, mean = mean, sd = sd, left = 0, right = Inf, log = TRUE)
    } else {
      function(y, mean, sd) dnorm(x = y, mean = mean, sd = sd, log = TRUE)
    }
    for(i in 1:(nrow(newdata))){
      par <- pred.par[i,]
      ll <- ll + distfun(newdata[i,1], mean = par[1], sd = par[2])
    }
    if(is.na(ll)) print("b.ll = NA")
    if(is.null(ll)) print("b.ll = NULL")
    return(ll)
  }
  
  
  
  ### set up function to get logLik from gamboostLSS object for new data and normal distribution
  gbll.newdata <- function(object, newdata){
    if(is.null(object)) {
      warning("error in gamboostLSS")
      return(NA)
    }
    ll <- 0
    nd <- newdata[,-c(1,ncol(newdata)-1, ncol(newdata))]
    pred.par <- predict(object, newdata = nd, parameter = list("mu", "sigma"), type = "response")
    pred_mu <- pred.par[[1]]
    pred_sigma <- pred.par[[2]]
    distfun <- if(censNO) {
      function(y, mean, sd) crch::dcnorm(x = y, mean = mean, sd = sd, left = 0, right = Inf, log = TRUE)
    } else {
      function(y, mean, sd) dnorm(x = y, mean = mean, sd = sd, log = TRUE)
    }
    for(i in 1:(nrow(newdata))){
      ll <- ll + distfun(newdata[i,1], mean = pred_mu[i], sd = pred_sigma[i], log = TRUE)
    }
    if(is.na(ll)) print("gb.ll = NA")
    if(is.null(ll)) print("gb.ll = NULL")
    return(ll)
  }
  
  
  
  ### set up function to get logLik from randomForest object for new data and normal distribution
  rfll.newdata <- function(object, newdata, learndata){
    if(is.null(object)) {
      warning("error in randomForest")
      return(NA)
    }
    if(!("y" %in% colnames(newdata))) stop("to get the loglikelihood for testdata argument 'newdata' has to include the response")
    ll <- 0
    pred_mu <- predict(object, newdata = newdata, type = "response")
    pred_sd <- numeric(length(NROW(newdata)))
    for(i in 1:NROW(newdata)){
      pred_sd[i] <- rf_getsd(object, newobs = newdata[i,], rfdata = learndata)
    }
    
    distfun <- if(censNO) {
      function(y, mean, sd) crch::dcnorm(x = y, mean = mean, sd = sd, left = 0, right = Inf, log = TRUE)
    } else {
      function(y, mean, sd) dnorm(x = y, mean = mean, sd = sd, log = TRUE)
    }
    for(i in 1:(nrow(newdata))){
      ll <- ll + distfun(newdata[i,1], mean = pred_mu[i], sd = pred_sd[i], log = TRUE)
    }
    
    if(is.na(ll)) print("rf.ll = NA")
    if(is.null(ll)) print("rf.ll = NULL")
    return(ll)
  }
  
  
  
  ### set up function to get logLik from cforest object for new data and normal distribution
  cfll.newdata <- function(object, newdata){
    if(is.null(object)) {
      warning("error in cforest")
      return(NA)
    }
    ll <- 0
    nd <- newdata[,-c(1,ncol(newdata)-1, ncol(newdata))]
    pred_mu <- predict(object, newdata = nd, type = "response")
    pred_sd <- cf_getsd(object, newdata = newdata)
    
    distfun <- if(censNO) {
      function(y, mean, sd) crch::dcnorm(x = y, mean = mean, sd = sd, left = 0, right = Inf, log = TRUE)
    } else {
      function(y, mean, sd) dnorm(x = y, mean = mean, sd = sd, log = TRUE)
    }
    for(i in 1:(nrow(newdata))){
      ll <- ll + distfun(newdata[i,1], mean = pred_mu[i], sd = pred_sd[i], log = TRUE)
    }

    if(is.na(ll)) print("cf.ll = NA")
    if(is.null(ll)) print("cf.ll = NULL")
    return(ll)
  }
  
  
  
  
  
  ######################################
  # smooth parameter function
  
  # parameter functions
  #fun <- function(x, kappa){
  #  mu <- 10 * (1-plogis(round(1.6^kappa,0) * x[,2]/10) + (2*plogis(round(1.6^kappa,0) * x[,2]/10) -1) * rbinom(nrow(x),1,exp(-(3*x[,1]-1.5)^(2*round(1.49^kappa,0)))))
  #  sigma <- 0.01 + mu/3
  #  par <- cbind(mu, sigma)
  #  return(par)
  #}
  
  
  if(onecov){
    if(fix.mu){
      if(fix.sigma){
        fun <- function(x, kappa){
          mu <- 0
          sigma <- 3
          par <- cbind(mu, sigma)
          return(par)
        }
      } else {
        fun <- function(x, kappa){
          mu <- 0
          sigma <- 1+3*abs(x[,1])
          par <- cbind(mu, sigma)
          return(par)
        }
      }
    } else {
      if(fix.sigma){
        fun <- function(x, kappa){
          mu <- 4 + 8 * (exp(-(4*x[,1]-0.7)^(2*kappa)))
          #mu <- 2 + 10 * (exp(-(4*x[,1]-2)^(2*kappa)))
          if(censNO) mu[x[,1]<0] <- 0
          sigma <- 3
          par <- cbind(mu, sigma)
          return(par)
        }
      } else {
        if(mu.sigma.interaction){
          fun <- function(x, kappa){
            mu <- 4 + 8 * (exp(-(4*x[,1]-0.7)^(2*kappa)))
            #mu <- 2 + 10 * (exp(-(4*x[,1]-2)^(2*kappa)))
            if(censNO) mu[x[,1]<0] <- 0
            sigma <- 2 + mu/4
            par <- cbind(mu, sigma)
            return(par)
          }
        } else {
          fun <- function(x, kappa){
            mu <- 4 + 8 * (exp(-(4*x[,1]-0.7)^(2*kappa)))
            #mu <- 2 + 10 * (exp(-(4*x[,1]-2)^(2*kappa)))
            if(censNO) mu[x[,1]<0] <- 0
            sigma <- 1+3*abs(x[,1])
            par <- cbind(mu, sigma)
            return(par)
          }
        }
      }
    }
    
  } else {
    # should the covariates have separated effects, e.g. mu depends on x1 only and sigma on x2 only (for fix.mu = fix.sigma = FALSE)
    if(cov.sep){
      if(mu.sigma.interaction) stop("if cov.sep is TRUE no interaction between mu and sigma can be included (mu.sigma.interaction has to be FALSE)")
      if(fix.mu) {
        if(fix.sigma) {
          fun <- function(x, kappa){
            mu <- 0
            sigma <- 3
            par <- cbind(mu, sigma)
            return(par)
          }
        } else {
          fun <- function(x, kappa){
            mu <- 0
            sigma <- 0.1 + 5 * (1-plogis((kappa^(1.8)) * (x[,2]-3)/10))
            par <- cbind(mu, sigma)
            return(par)
          }
        }
      } else {
        if(fix.sigma) {
          fun <- function(x, kappa){
            mu <- 10 * (exp(-(4*x[,1]-2)^(2*kappa)))
            sigma <- 3
            par <- cbind(mu, sigma)
            return(par)
          }
        } else {
          fun <- function(x, kappa){
            mu <- 10 * (exp(-(4*x[,1]-2)^(2*kappa)))
            sigma <- 0.1 + 5 * (1-plogis((kappa^(1.8)) * (x[,2]-3)/10))
            par <- cbind(mu, sigma)
            return(par)
          }
        }
      }
      # or should mu depend on x1 and x2 and sigma on x2 (for fix.mu = fix.sigma = mu.sigma.interaction = FALSE)
    } else {
      if(fix.mu) {
        if(fix.sigma) {
          fun <- function(x, kappa){
            mu <- 0
            sigma <- 3
            par <- cbind(mu, sigma)
            return(par)
          }
        } else {
          fun <- function(x, kappa){
            mu <- 0
            sigma <- 0.1 + abs(x[,2])
            par <- cbind(mu, sigma)
            return(par)
          }
        }
      } else {
        if(fix.sigma) {
          fun <- function(x, kappa){
            mu <- 10 * (1-plogis((kappa^(1.8)) * (x[,2]-3)/10) 
                        + (2*plogis((kappa^(1.8)) * (x[,2]-3)/10) -1) 
                        * rbinom(NROW(x),1,exp(-(4*x[,1]-2)^(2*kappa))))
            sigma <- 3
            par <- cbind(mu, sigma)
            return(par)
          }
        } else {
          if(mu.sigma.interaction){
            fun <- function(x, kappa){
              mu <- 10 * (1-plogis((kappa^(1.8)) * (x[,2]-3)/10) 
                          + (2*plogis((kappa^(1.8)) * (x[,2]-3)/10) -1) 
                          * rbinom(NROW(x),1,exp(-(4*x[,1]-2)^(2*kappa))))
              sigma <- 2 + mu/2.5
              par <- cbind(mu, sigma)
              return(par)
            }
          } else {
            fun <- function(x, kappa){
              mu <- 10 * (1-plogis((kappa^(1.8)) * (x[,2]-3)/10) 
                          + (2*plogis((kappa^(1.8)) * (x[,2]-3)/10) -1) 
                          * rbinom(NROW(x),1,exp(-(4*x[,1]-2)^(2*kappa))))
              sigma <- 1 + abs(x[,2])
              par <- cbind(mu, sigma)
              return(par)
            }
          }
        }
      }
    }
  }
  

  ### simulation
  {
    
    # matrix of results
    resmat <- NULL
    
    if(eval_disttree){
      res_dt <- mclapply(1:(nsteps+1),
                         function(i){
                           
                           # for each step the parameter function is defined
                           f <- function(x) {fun(x, kappa.start + (i-1)*stepsize)} 
                           
                           rmse.exp.true.dt <- numeric(length = nrep)
                           rmse.exp.obs.dt <- numeric(length = nrep)
                           rmse.mu.dt <- numeric(length = nrep)
                           rmse.sigma.dt <- numeric(length = nrep)
                           loglik.dt <- numeric(length = nrep)
                           
                           for(j in 1:nrep){
                             set.seed(seedconst+((i-1)*nrep)+j)
                             learndata <- dgp(nobs, family = family, round.sp = 4, fun = f)
                             newdata <- dgp(testnobs, family = family, round.sp = 4, fun = f)
                             nd <- newdata[,-c(1,ncol(newdata)-1,ncol(newdata))]
                             
                             if(type.tree == "mob"){
                               dt <- disttree(formula, data=learndata, family=family, type.tree = "mob", 
                                              control = mob_control(restart = FALSE, numsplit = "center", 
                                                                    alpha = 1-tree_mincrit, minsplit = tree_minsplit))
                             }
                             if(type.tree == "ctree"){
                               dt <- disttree(formula, data=learndata, family=family, type.tree = "ctree", 
                                              control = ctree_control(teststat = "quad", testtype = "Bonferroni", intersplit = TRUE,
                                                                      mincriterion = tree_mincrit, minsplit = tree_minsplit))
                             }
                             
                             # get true and predicted parameters for newdata
                             true_mu <- newdata[,"mu"]
                             true_sigma <- newdata[,"sigma"]
                             
                             dt_coef <- if(is.vector(coef(dt))){
                               t(as.data.frame(coef(dt)))
                             } else {
                               coef(dt)[paste(as.vector(predict(dt, newdata = nd, type = "node"))),]
                             }
                             dt_mu <- as.vector(dt_coef[,"mu"]) 
                             dt_sigma <- as.vector(dt_coef[,"sigma"]) 
                             
                             if(family$censored){
                               #calculate expected value for censored data
                               true_exp <- pnorm(true_mu/true_sigma) * (true_mu + true_sigma * (dnorm(true_mu/true_sigma) / pnorm(true_mu/true_sigma)))
                               dt_exp <- pnorm(dt_mu/dt_sigma) * (dt_mu + dt_sigma * (dnorm(dt_mu/dt_sigma) / pnorm(dt_mu/dt_sigma)))
                             } else {
                               # for non censored normal data the expected value is the location paramter
                               true_exp <- true_mu
                               dt_exp <- dt_mu
                             }
                             
                             # calculate RMSEs
                             rmse.exp.true.dt[j] <- sqrt(mean((dt_exp - true_exp)^2))
                             rmse.exp.obs.dt[j] <- sqrt(mean((dt_exp - newdata[,"y"])^2))
                             rmse.mu.dt[j] <- sqrt(mean((dt_mu - true_mu)^2))
                             rmse.sigma.dt[j] <- sqrt(mean((dt_sigma - true_sigma)^2))
                             
                             # calculate loglikelihood for newdata
                             loglik.dt[j] <- logLik(dt, newdata = newdata)
                           }
                           
                           
                           return(as.data.frame(c(mean(rmse.exp.true.dt),
                                                  mean(rmse.exp.obs.dt),
                                                  mean(rmse.mu.dt),
                                                  mean(rmse.sigma.dt),
                                                  mean(loglik.dt)),
                                                row.names = c( "av.rmse.exp.true.dt",
                                                               "av.rmse.exp.obs.dt",
                                                               "av.rmse.mu.dt",
                                                               "av.rmse.sigma.dt",
                                                               "av.loglik.dt")))
                           
                         },
                         mc.cores = detectCores() - 1
      )
      
      resmat_dt <- matrix(ncol = 5, nrow = (nsteps+1))
      colnames(resmat_dt) <- rownames(res_dt[[1]])
      for(k in 1:(nsteps+1)) resmat_dt[k,] <- t(as.vector(res_dt[[k]]))
      
      resmat <- cbind(resmat, resmat_dt)
    }
    
    
    
    if(eval_distforest){
      res_df <- mclapply(1:(nsteps+1),
                         function(i){
                           
                           # for each step the parameter function is defined
                           f <- function(x) {fun(x, kappa.start + (i-1)*stepsize)} 
                           
                           rmse.exp.true.df <- numeric(length = nrep)
                           rmse.exp.obs.df <- numeric(length = nrep)
                           rmse.mu.df <- numeric(length = nrep)
                           rmse.sigma.df <- numeric(length = nrep)
                           loglik.df <- numeric(length = nrep)
                           
                           for(j in 1:nrep){
                             set.seed(seedconst+((i-1)*nrep)+j)
                             learndata <- dgp(nobs, family = family, round.sp = 4, fun = f)
                             newdata <- dgp(testnobs, family = family, round.sp = 4, fun = f)
                             nd <- newdata[,-c(1,ncol(newdata)-1,ncol(newdata))]
                             
                             if(type.tree == "mob"){
                               df <- distforest(formula, data=learndata, family=family, type.tree = "mob", ntree = ntree,
                                                control = mob_control(restart = FALSE, numsplit = "center", 
                                                                      alpha = 1-forest_mincrit, minsplit = forest_minsplit))
                             }
                             if(type.tree == "ctree"){
                               df <- distforest(formula, data=learndata, family=family, type.tree = "ctree", ntree = ntree,
                                                control = ctree_control(teststat = "quad", testtype = "Univ", intersplit = TRUE,
                                                                        mincriterion = forest_mincrit, minsplit = forest_minsplit))
                             }
                             
                             # get true and predicted parameters for newdata
                             true_mu <- newdata[,"mu"]
                             true_sigma <- newdata[,"sigma"]
                             
                             df_coef <- predict(df, newdata = nd, type = "parameter") # returns fitted values and fitted parameters
                             df_mu <- df_coef$mu
                             df_sigma <- df_coef$sigma
                             
                             if(family$censored){
                               #calculate expected value for censored data
                               true_exp <- pnorm(true_mu/true_sigma) * (true_mu + true_sigma * (dnorm(true_mu/true_sigma) / pnorm(true_mu/true_sigma)))
                               df_exp <- pnorm(df_mu/df_sigma) * (df_mu + df_sigma * (dnorm(df_mu/df_sigma) / pnorm(df_mu/df_sigma)))
                             } else {
                               # for non censored normal data the expected value is the location paramter
                               true_exp <- true_mu
                               df_exp <- df_mu
                             }
                             
                             # calculate RMSEs
                             rmse.exp.true.df[j] <- sqrt(mean((df_exp - true_exp)^2))
                             rmse.exp.obs.df[j] <- sqrt(mean((df_exp - newdata[,"y"])^2))
                             rmse.mu.df[j] <- sqrt(mean((df_mu - true_mu)^2))
                             rmse.sigma.df[j] <- sqrt(mean((df_sigma - true_sigma)^2))
                             
                             # calculate loglikelihood for newdata
                             loglik.df[j] <- as.numeric(logLik(df, newdata = newdata))
                           }
                           
                           
                           return(as.data.frame(c(mean(rmse.exp.true.df),
                                                  mean(rmse.exp.obs.df),
                                                  mean(rmse.mu.df),
                                                  mean(rmse.sigma.df),
                                                  mean(loglik.df)),
                                                row.names = c( "av.rmse.exp.true.df",
                                                               "av.rmse.exp.obs.df",
                                                               "av.rmse.mu.df",
                                                               "av.rmse.sigma.df",
                                                               "av.loglik.df")))
                           
                         },
                         mc.cores = detectCores() - 1
      )
      
      resmat_df <- matrix(ncol = 5, nrow = (nsteps+1))
      colnames(resmat_df) <- rownames(res_df[[1]])
      for(k in 1:(nsteps+1)) resmat_df[k,] <- t(as.vector(res_df[[k]]))
      
      resmat <- cbind(resmat, resmat_df)
    }
    
    
    
    if(eval_gamlss){
      res_g <- mclapply(1:(nsteps+1),
                        function(i){
                          
                          if(censNO) gen.cens(NO, type = "left")
                          
                          # for each step the parameter function is defined
                          f <- function(x) {fun(x, kappa.start + (i-1)*stepsize)} 
                          
                          rmse.exp.true.g <- numeric(length = nrep)
                          rmse.exp.obs.g <- numeric(length = nrep)
                          rmse.mu.g <- numeric(length = nrep)
                          rmse.sigma.g <- numeric(length = nrep)
                          loglik.g <- numeric(length = nrep)
                          
                          # FIX ME: notation only works with x1 and x2 as regressors
                          formula_rh <- formula[[3]]
                          if("x1" %in% as.character(formula_rh)){
                            mu.formula <- y~pb(x1)
                            sigma.formula <- ~pb(x1)
                            if("x2" %in% as.character(formula_rh)){
                              mu.formula <- y~pb(x1)+pb(x2)
                              sigma.formula <- ~pb(x1)+pb(x2)
                            }
                          } 
                          
                          
                          for(j in 1:nrep){
                            set.seed(seedconst+((i-1)*nrep)+j)
                            learndata <- dgp(nobs, family = family, round.sp = 4, fun = f)
                            newdata <- dgp(testnobs, family = family, round.sp = 4, fun = f)
                            nd <- newdata[,-c(1,ncol(newdata)-1,ncol(newdata))]
                            
                            if(censNO){
                              learndata$y <- Surv(learndata$y, learndata$y>0, type="left")
                              #newdata$y <- Surv(newdata$y, newdata$y>0, type="left")
                              g <- gamlss(formula = mu.formula, sigma.formula = sigma.formula, data = learndata, family = NOlc())
                            } else {
                              g <- gamlss(formula = mu.formula, sigma.formula = sigma.formula, data = learndata, family = NO())
                            }
                            
                            # get true and predicted parameters for newdata
                            true_mu <- newdata[,"mu"]
                            true_sigma <- newdata[,"sigma"]
                            
                            g_mu <- predict(g, newdata = nd, what = "mu", type = "response", data = learndata)
                            g_sigma <- predict(g, newdata = nd, what = "sigma", type = "response", data = learndata)
                            
                            if(family$censored){
                              #calculate expected value for censored data
                              true_exp <- pnorm(true_mu/true_sigma) * (true_mu + true_sigma * (dnorm(true_mu/true_sigma) / pnorm(true_mu/true_sigma)))
                              g_exp <- pnorm(g_mu/g_sigma) * (g_mu + g_sigma * (dnorm(g_mu/g_sigma) / pnorm(g_mu/g_sigma)))
                            } else {  
                              # for non censored normal data the expected value is the location paramter
                              true_exp <- true_mu
                              g_exp <- g_mu
                            }
                            
                            
                            # calculate RMSEs
                            rmse.exp.true.g[j] <- sqrt(mean((g_exp - true_exp)^2))
                            rmse.exp.obs.g[j] <- sqrt(mean((g_exp - newdata[,"y"])^2))
                            rmse.mu.g[j] <- sqrt(mean((g_mu  - true_mu)^2))
                            rmse.sigma.g[j] <- sqrt(mean((g_sigma - true_sigma)^2))
                            
                            # calculate loglikelihood for newdata
                            loglik.g[j] <- gll.newdata(g, newdata = newdata, data = learndata)
                          }
                          
                          
                          return(as.data.frame(c(mean(rmse.exp.true.g, na.rm = TRUE),
                                                 mean(rmse.exp.obs.g, na.rm = TRUE),
                                                 mean(rmse.mu.g, na.rm = TRUE),
                                                 mean(rmse.sigma.g, na.rm = TRUE),
                                                 mean(loglik.g, na.rm = TRUE)),
                                               row.names = c( "av.rmse.exp.true.g",
                                                              "av.rmse.exp.obs.g",
                                                              "av.rmse.mu.g",
                                                              "av.rmse.sigma.g",
                                                              "av.loglik.g")))
                          
                        },
                        mc.cores = detectCores() - 1
      )
      
      resmat_g <- matrix(ncol = 5, nrow = (nsteps+1))
      colnames(resmat_g) <- rownames(res_g[[1]])
      for(k in 1:(nsteps+1)) resmat_g[k,] <- t(as.vector(res_g[[k]]))
      
      resmat <- cbind(resmat, resmat_g)
    }
    
    
    ## FIX ME: set up bamlss for censored data (as.surv(y))
    if(eval_bamlss){
      res_b <- mclapply(1:(nsteps+1),
                        function(i){
                          
                          # for each step the parameter function is defined
                          f <- function(x) {fun(x, kappa.start + (i-1)*stepsize)} 
                          
                          rmse.exp.true.b <- numeric(length = nrep)
                          rmse.exp.obs.b <- numeric(length = nrep)
                          rmse.mu.b <- numeric(length = nrep)
                          rmse.sigma.b <- numeric(length = nrep)
                          loglik.b <- numeric(length = nrep)
                          
                          # FIX ME: notation only works with x1 and x2 as regressors
                          formula_rh <- formula[[3]]
                          if("x1" %in% as.character(formula_rh)){
                            mu.formula <- y~s(x1)
                            sigma.formula <- ~s(x1)
                            if("x2" %in% as.character(formula_rh)){
                              mu.formula <- y~s(x1)+s(x2)
                              sigma.formula <- ~s(x1)+s(x2)
                            }
                          }
                          
                          for(j in 1:nrep){
                            set.seed(seedconst+((i-1)*nrep)+j)
                            learndata <- dgp(nobs, family = family, round.sp = 4, fun = f)
                            newdata <- dgp(testnobs, family = family, round.sp = 4, fun = f)
                            nd <- newdata[,-c(1,ncol(newdata)-1,ncol(newdata))]
                            
                            if(family$censore) {
                              ## FIX ME: set up for bamlss
                              #learndata$y <- Surv(learndata$y, learndata$y>0, type="left")
                              b <- bamlss(bamlss.formula(list(mu.formula, sigma.formula), family = cnorm_bamlss), data = learndata)
                            } else {
                              b <- bamlss(bamlss.formula(list(mu.formula, sigma.formula), family = "gaussian"), data = learndata)
                            }
                            
                            # get true and predicted parameters for newdata
                            true_mu <- newdata[,"mu"]
                            true_sigma <- newdata[,"sigma"]
                            
                            b_coef <- predict(b, newdata = nd, what = "parameters", type = "parameter", data = learndata)
                            b_mu <- b_coef$mu
                            b_sigma <- b_coef$sigma
                            
                            if(family$censored){
                              
                              #calculate expected value for censored data
                              true_exp <- pnorm(true_mu/true_sigma) * (true_mu + true_sigma * (dnorm(true_mu/true_sigma) / pnorm(true_mu/true_sigma)))
                              b_exp <- pnorm(b_mu/b_sigma) * (b_mu + b_sigma * (dnorm(b_mu/b_sigma) / pnorm(b_mu/b_sigma)))
                              
                            } else {
                              
                              # for non censored normal data the expected value is the location paramter
                              true_exp <- true_mu
                              b_exp <- b_mu
                            }
                            
                            
                            # calculate RMSEs
                            rmse.exp.true.b[j] <- sqrt(mean((b_exp - true_exp)^2))
                            rmse.exp.obs.b[j] <- sqrt(mean((b_exp - newdata[,"y"])^2))
                            rmse.mu.b[j] <- sqrt(mean((b_mu  - true_mu)^2))
                            rmse.sigma.b[j] <- sqrt(mean((b_sigma - true_sigma)^2))
                            
                            # calculate loglikelihood for newdata
                            loglik.b[j] <- bll.newdata(b, newdata = newdata)
                          }
                          
                          
                          return(as.data.frame(c(mean(rmse.exp.true.b),
                                                 mean(rmse.exp.obs.b),
                                                 mean(rmse.mu.b),
                                                 mean(rmse.sigma.b),
                                                 mean(loglik.b)),
                                               row.names = c( "av.rmse.exp.true.b",
                                                              "av.rmse.exp.obs.b",
                                                              "av.rmse.mu.b",
                                                              "av.rmse.sigma.b",
                                                              "av.loglik.b")))
                          
                        },
                        mc.cores = detectCores() - 1
      )
      
      resmat_b <- matrix(ncol = 5, nrow = (nsteps+1))
      colnames(resmat_b) <- rownames(res_b[[1]])
      for(k in 1:(nsteps+1)) resmat_b[k,] <- t(as.vector(res_b[[k]]))
      
      resmat <- cbind(resmat, resmat_b)
    }
    
    
    ## FIX ME: censored?
    if(eval_gamboostLSS){
      
      ##grid for cvrisk to find optimal mstop:
      
      #for method="cyclic"
      # grid <- make.grid(max = c(mu = 400, sigma = 200), min = 10, 
      #                       length.out = 10, dense_mu_grid = TRUE)
      
      #for method="noncyclic"
      grid <- make.grid(max = 300, min = 10, length.out = 10)
      
      res_gb <- mclapply(1:(nsteps+1),
                         function(i){
                           
                           if(censNO) gen.cens(NO, type = "left")
                           
                           # for each step the parameter function is defined
                           f <- function(x) {fun(x, kappa.start + (i-1)*stepsize)} 
                           
                           rmse.exp.true.gb <- numeric(length = nrep)
                           rmse.exp.obs.gb <- numeric(length = nrep)
                           rmse.mu.gb <- numeric(length = nrep)
                           rmse.sigma.gb <- numeric(length = nrep)
                           loglik.gb <- numeric(length = nrep)
                           
                           # FIX ME: notation only works with x1 and x2 as regressors
                           formula_rh <- formula[[3]]
                           if("x1" %in% as.character(formula_rh)){
                             mu.formula <- y~bbs(x1)
                             sigma.formula <- y~bbs(x1)
                             if("x2" %in% as.character(formula_rh)){
                               mu.formula <- y~bbs(x1)+bbs(x2)
                               sigma.formula <- y~bbs(x1)+bbs(x2)
                             }
                           } 
                           
                           for(j in 1:nrep){
                             set.seed(seedconst+((i-1)*nrep)+j)
                             learndata <- dgp(nobs, family = family, round.sp = 4, fun = f)
                             newdata <- dgp(testnobs, family = family, round.sp = 4, fun = f)
                             nd <- newdata[,-c(1,ncol(newdata)-1,ncol(newdata))]
                             
                             # FIX ME: censored data
                             if(family$censored) {
                               learndata$y <- Surv(learndata$y, learndata$y>0, type="left")
                               #newdata$y <- Surv(newdata$y, newdata$y>0, type="left")
                               gb <- gamboostLSS(formula = list(mu = mu.formula, sigma =sigma.formula), data = learndata, 
                                                 families = as.families(fname = NOlc()), method = "noncyclic",
                                                 control = boost_control(mstop = 400L))
                               if(gamboost_cvr){
                                 cvr <- cvrisk(gb, grid = grid)
                                 mstop(gb) <- mstop(cvr)
                               }
                             } else {
                               #gb <- gamboostLSS(formula = list(mu = mu.formula, sigma =sigma.formula), data = learndata, families = GaussianLSS())
                               #gb <- gamboostLSS(formula = list(mu = mu.formula, sigma =sigma.formula), data = learndata, 
                               #                  families = GaussianLSS(), control = boost_control(mstop = list(mu = 400, sigma = 200)))
                               gb <- gamboostLSS(formula = list(mu = mu.formula, sigma =sigma.formula), data = learndata, 
                                                 families = GaussianLSS(), method = "noncyclic",
                                                 control = boost_control(mstop = 400L))
                               if(gamboost_cvr){
                                 cvr <- cvrisk(gb, grid = grid)
                                 mstop(gb) <- mstop(cvr)
                               }
                             }                            
                             
                             # get true and predicted parameters for newdata
                             true_mu <- newdata[,"mu"]
                             true_sigma <- newdata[,"sigma"]
                             
                             gb.pred.par <- predict(gb, newdata = nd, parameter = list("mu","sigma"), type = "response")
                             gb_mu <- gb.pred.par[[1]]
                             gb_sigma <- gb.pred.par[[2]]
                             
                             if(family$censored){
                               
                               #calculate expected value for censored data
                               true_exp <- pnorm(true_mu/true_sigma) * (true_mu + true_sigma * (dnorm(true_mu/true_sigma) / pnorm(true_mu/true_sigma)))
                               gb_exp <- pnorm(gb_mu/gb_sigma) * (gb_mu + gb_sigma * (dnorm(gb_mu/gb_sigma) / pnorm(gb_mu/gb_sigma)))
                             } else {
                               
                               # for non censored normal data the expected value is the location paramter
                               true_exp <- true_mu
                               gb_exp <- gb_mu
                             }
                             
                             
                             # calculate RMSEs
                             rmse.exp.true.gb[j] <- sqrt(mean((gb_exp - true_exp)^2))
                             rmse.exp.obs.gb[j] <- sqrt(mean((gb_exp - newdata[,"y"])^2))
                             rmse.mu.gb[j] <- sqrt(mean((gb_mu - true_mu)^2))
                             rmse.sigma.gb[j] <- sqrt(mean((gb_sigma - true_sigma)^2))
                             
                             # calculate loglikelihood for newdata
                             loglik.gb[j] <- gbll.newdata(gb, newdata = newdata)
                           }
                           
                           
                           return(as.data.frame(c(mean(rmse.exp.true.gb),
                                                  mean(rmse.exp.obs.gb),
                                                  mean(rmse.mu.gb),
                                                  mean(rmse.sigma.gb),
                                                  mean(loglik.gb)),
                                                row.names = c( "av.rmse.exp.true.gb",
                                                               "av.rmse.exp.obs.gb",
                                                               "av.rmse.mu.gb",
                                                               "av.rmse.sigma.gb",
                                                               "av.loglik.gb")))
                           
                         },
                         mc.cores = detectCores() - 1
      )
      
      resmat_gb <- matrix(ncol = 5, nrow = (nsteps+1))
      colnames(resmat_gb) <- rownames(res_gb[[1]])
      for(k in 1:(nsteps+1)) resmat_gb[k,] <- t(as.vector(res_gb[[k]]))
      
      resmat <- cbind(resmat, resmat_gb)
    }
    
    
    
    if(eval_randomForest){
      res_rf <- mclapply(1:(nsteps+1),
                         function(i){
                           
                           # for each step the parameter function is defined
                           f <- function(x) {fun(x, kappa.start + (i-1)*stepsize)} 
                           
                           rmse.exp.true.rf <- numeric(length = nrep)
                           rmse.exp.obs.rf <- numeric(length = nrep)
                           rmse.mu.rf <- numeric(length = nrep)
                           rmse.sigma.rf <- numeric(length = nrep)
                           loglik.rf <- numeric(length = nrep)
                           
                           for(j in 1:nrep){  
                             set.seed(seedconst+((i-1)*nrep)+j)
                             learndata <- dgp(nobs, family = family, round.sp = 4, fun = f)
                             newdata <- dgp(testnobs, family = family, round.sp = 4, fun = f)
                             nd <- newdata[,-c(1,ncol(newdata)-1,ncol(newdata))]
                             
                             rf <- randomForest(formula, data = learndata, ntree = ntree, 
                                                nodesize = forest_minsplit, 
                                                keep.inbag = TRUE, replace = FALSE)
                             
                             # get true and predicted parameters for newdata
                             true_mu <- newdata[,"mu"]
                             true_sigma <- newdata[,"sigma"]
                             
                             rf_mu <- predict(rf, newdata = newdata)
                             
                             ## get 'fitted.sigma.rf'
                             rf_sigma <- numeric(length(NROW(newdata)))
                             for(l in 1:NROW(newdata)){
                               rf_sigma[l] <- rf_getsd(rf, newobs = newdata[l,], rfdata = learndata)
                             }
                             
                             if(family$censored){
                               
                               #calculate expected value for censored data
                               true_exp <- pnorm(true_mu/true_sigma) * (true_mu + true_sigma * (dnorm(true_mu/true_sigma) / pnorm(true_mu/true_sigma)))
                               rf_exp <- pnorm(rf_mu/rf_sigma) * (rf_mu + rf_sigma * (dnorm(rf_mu/rf_sigma) / pnorm(rf_mu/rf_sigma)))
                               
                             } else {
                               
                               # for non censored normal data the expected value is the location paramter
                               true_exp <- true_mu
                               rf_exp <- rf_mu
                             }                            
                             
                             
                             # calculate RMSEs
                             rmse.exp.true.rf[j] <- sqrt(mean((rf_exp - true_exp)^2))
                             rmse.exp.obs.rf[j] <- sqrt(mean((rf_exp - newdata[,"y"])^2))
                             rmse.mu.rf[j] <- sqrt(mean((rf_mu - true_mu)^2))
                             rmse.sigma.rf[j] <- sqrt(mean((rf_sigma - true_sigma)^2))
                             
                             # calculate loglikelihood for newdata
                             loglik.rf[j] <- rfll.newdata(rf, newdata = newdata, learndata = learndata)
                           }
                           
                           
                           return(as.data.frame(c(mean(rmse.exp.true.rf),
                                                  mean(rmse.exp.obs.rf),
                                                  mean(rmse.mu.rf),
                                                  mean(rmse.sigma.rf),
                                                  mean(loglik.rf)),
                                                row.names = c( "av.rmse.exp.true.rf",
                                                               "av.rmse.exp.obs.rf",
                                                               "av.rmse.mu.rf",
                                                               "av.rmse.sigma.rf",
                                                               "av.loglik.rf")))
                         },
                         mc.cores = detectCores() - 1
      )
      
      resmat_rf <- matrix(ncol = 5, nrow = (nsteps+1))
      colnames(resmat_rf) <- rownames(res_rf[[1]])
      for(k in 1:(nsteps+1)) resmat_rf[k,] <- t(as.vector(res_rf[[k]]))
      
      resmat <- cbind(resmat, resmat_rf)
    }
    
    
    
    if(eval_cforest){
      res_cf <- mclapply(1:(nsteps+1),
                         function(i){
                           
                           # for each step the parameter function is defined
                           f <- function(x) {fun(x, kappa.start + (i-1)*stepsize)} 
                           
                           rmse.exp.true.cf <- numeric(length = nrep)
                           rmse.exp.obs.cf <- numeric(length = nrep)
                           rmse.mu.cf <- numeric(length = nrep)
                           rmse.sigma.cf <- numeric(length = nrep)
                           loglik.cf <- numeric(length = nrep)
                           
                           for(j in 1:nrep){
                             set.seed(seedconst+((i-1)*nrep)+j)
                             learndata <- dgp(nobs, family = family, round.sp = 4, fun = f)
                             newdata <- dgp(testnobs, family = family, round.sp = 4, fun = f)
                             nd <- newdata[,-c(1,ncol(newdata)-1,ncol(newdata))]
                             
                             cf <- cforest(formula, data = learndata, ntree = ntree,
                                           control = ctree_control(teststat = "quad", testtype = "Univ", 
                                                                   mincriterion = forest_mincrit, minsplit = forest_minsplit, 
                                                                   intersplit = TRUE))
                             
                             
                             # get true and predicted parameters for newdata
                             true_mu <- newdata[,"mu"]
                             true_sigma <- newdata[,"sigma"]
                             
                             cf_mu <- predict(cf, newdata = newdata, type = "response")
                             
                             ## get 'fitted.sigma.cf'
                             cf_sigma <- cf_getsd(cf, newdata = newdata)
                             
                             
                             if(family$censored){
                               
                               #calculate expected value for censored data
                               true_exp <- pnorm(true_mu/true_sigma) * (true_mu + true_sigma * (dnorm(true_mu/true_sigma) / pnorm(true_mu/true_sigma)))
                               cf_exp <- pnorm(cf_mu/cf_sigma) * (cf_mu + cf_sigma * (dnorm(cf_mu/cf_sigma) / pnorm(cf_mu/cf_sigma)))
                               
                             } else {
                               
                               # for non censored normal data the expected value is the location paramter
                               true_exp <- true_mu
                               cf_exp <- cf_mu
                             }                            
                             
                             
                             # calculate RMSEs
                             rmse.exp.true.cf[j] <- sqrt(mean((cf_exp - true_exp)^2))
                             rmse.exp.obs.cf[j] <- sqrt(mean((cf_exp - newdata[,"y"])^2))
                             rmse.mu.cf[j] <- sqrt(mean((cf_mu - true_mu)^2))
                             rmse.sigma.cf[j] <- sqrt(mean((cf_sigma - true_sigma)^2))
                             
                             # calculate loglikelihood for newdata
                             loglik.cf[j] <- cfll.newdata(cf, newdata = newdata)
                           }
                           
                           
                           return(as.data.frame(c(mean(rmse.exp.true.cf),
                                                  mean(rmse.exp.obs.cf),
                                                  mean(rmse.mu.cf),
                                                  mean(rmse.sigma.cf),
                                                  mean(loglik.cf)),
                                                row.names = c( "av.rmse.exp.true.cf",
                                                               "av.rmse.exp.obs.cf",
                                                               "av.rmse.mu.cf",
                                                               "av.rmse.sigma.cf",
                                                               "av.loglik.cf")))
                         },
                         mc.cores = detectCores() - 1
      )
      
      resmat_cf <- matrix(ncol = 5, nrow = (nsteps+1))
      colnames(resmat_cf) <- rownames(res_cf[[1]])
      for(k in 1:(nsteps+1)) resmat_cf[k,] <- t(as.vector(res_cf[[k]]))
      
      resmat <- cbind(resmat, resmat_cf)
    }
    
    
    
    
    
    
    x.axis <- seq(kappa.start, (kappa.start + (nsteps)*stepsize), stepsize)
    
    simres <- list(resmat = resmat, x.axis = x.axis)
    
  }
  
  
  
  simres$fun <- fun
  simres$call <- cl
  return(simres)
}





#simres <- gensim(seedconst = 74, nrep = 10, ntree = 300L,
#                 eval_disttree = TRUE,
#                 eval_distforest = FALSE,
#                 eval_gamlss = FALSE,
#                 eval_bamlss = FALSE,
#                 eval_gamboostLSS = FALSE,
#                 eval_randomForest = TRUE,
#                 eval_cforest = TRUE)


## plot functions
# plot RMSE
#plot_rmse <- function(simres, type = c("exp", "par")){
#  
#  if(type == "exp"){
#    rmse <- cbind(simres$resmat[,"av.rmse.exp.true.dt"], simres$resmat[,"av.rmse.exp.true.df"], 
#                  simres$resmat[,"av.rmse.exp.true.g"], simres$resmat[,"av.rmse.exp.true.b"],
#                  #simres$resmat[,"av.rmse.exp.true.rf"],
#                  simres$resmat[,"av.rmse.exp.obs.dt"], simres$resmat[,"av.rmse.exp.obs.df"], 
#                  simres$resmat[,"av.rmse.exp.obs.g"], simres$resmat[,"av.rmse.exp.obs.b"]#,
#                  #simres$resmat[,"av.rmse.exp.obs.rf"]
#                  )
#    colnames(rmse) <- c("dt.true", "df.true", "g.true", "b.true", #"rf.true", 
#                        "dt.obs", "df.obs", "g.obs", "b.obs" #, "rf.obs"
#                        )
#    ylim <- c(min(rmse), max(rmse))
#    
#    plot(x = simres$x.axis, y = rmse[,"dt.true"], type = "l", col = "forestgreen", ylim = ylim,
#         xlab = "kappa", ylab = "RMSE", main = "disttree vs distforest vs gamlss vs randomForest")
#    lines(x = simres$x.axis, y = rmse[,"dt.obs"], type = "l", lty = 2, col = 'forestgreen')
#    lines(x = simres$x.axis, y = rmse[,"df.true"], type = "l", col = 'red')
#    lines(x = simres$x.axis, y = rmse[,"df.obs"], type = "l", lty = 2, col = 'red')
#    lines(x = simres$x.axis, y = rmse[,"g.true"], type = "l", col = 'blue')
#    lines(x = simres$x.axis, y = rmse[,"g.obs"], type = "l", lty = 2, col = 'blue')
#    lines(x = simres$x.axis, y = rmse[,"b.true"], type = "l", col = 'purple')
#    lines(x = simres$x.axis, y = rmse[,"b.obs"], type = "l", lty = 2, col = 'purple')
#    #lines(x = simres$x.axis, y = rmse[,"rf.true"], type = "l", col = 'yellow')
#    #lines(x = simres$x.axis, y = rmse[,"rf.obs"], type = "l", lty = 2, col = 'yellow')
#    
#    legend('topleft', c("disttree", "distforest", "gamlss", "bamlss"), 
#           col = c('forestgreen', 'red', 'blue', 'purple'), lty = 1, cex = 0.7)
#  } 
#  
#  if(type == "par"){
#    rmse <- cbind(simres$resmat[,"av.rmse.mu.dt"], simres$resmat[,"av.rmse.mu.df"], 
#                  simres$resmat[,"av.rmse.mu.g"], simres$resmat[,"av.rmse.mu.b"], 
#                  simres$resmat[,"av.rmse.sigma.dt"], simres$resmat[,"av.rmse.sigma.df"], 
#                  simres$resmat[,"av.rmse.sigma.g"], simres$resmat[,"av.rmse.sigma.b"])
#    colnames(rmse) <- c("mu.dt", "mu.df", "mu.g", "mu.b", "sigma.dt", "sigma.df", "sigma.g", "sigma.b")
#    ylim <- c(min(rmse), max(rmse))
#    plot(x = simres$x.axis, y = rmse[,"mu.dt"], type = "l", col = "forestgreen", ylim = ylim,
#         xlab = "kappa", ylab = "RMSE", main = "disttree vs distforest vs gamlss")
#    lines(x = simres$x.axis, y = rmse[,"sigma.dt"], type = "l", lty = 2, col = 'forestgreen')
#    lines(x = simres$x.axis, y = rmse[,"mu.df"], type = "l", col = 'red')
#    lines(x = simres$x.axis, y = rmse[,"sigma.df"], type = "l", lty = 2, col = 'red')
#    lines(x = simres$x.axis, y = rmse[,"mu.g"], type = "l", col = 'blue')
#    lines(x = simres$x.axis, y = rmse[,"sigma.g"], type = "l", lty = 2, col = 'blue')
#    lines(x = simres$x.axis, y = rmse[,"mu.b"], type = "l", col = 'purple')
#    lines(x = simres$x.axis, y = rmse[,"sigma.b"], type = "l", lty = 2, col = 'purple')
#    
#    legend('topleft', c("disttree", "distforest", "gamlss", "bamlss"), 
#           col = c('forestgreen', 'red', 'blue', 'purple'), lty = 1, cex = 0.7)
#  }
#}

# plot loglikelihood
#plot_ll <- function(simlist){
#  ll <- cbind(simres$resmat[,"av.loglik.dt"], simres$resmat[,"av.loglik.df"], 
#              simres$resmat[,"av.loglik.g"], simres$resmat[,"av.loglik.b"])
#  colnames(ll) <- c("dt", "df","g","b")
#  ylim <- c(min(ll), max(ll))
#  plot(x = simres$x.axis, y = ll[,"dt"], type = "l", col = "forestgreen", ylim = ylim,
#       xlab = "kappa", ylab = "log-likelihood", main = "disttree vs distforest vs gamlss")
#  lines(x = simres$x.axis, y = ll[,"df"], type = "l", col = 'red')
#  lines(x = simres$x.axis, y = ll[,"g"], type = "l", col = 'blue')
#  lines(x = simres$x.axis, y = ll[,"b"], type = "l", col = 'purple')
#  legend('topleft', c("disttree", "distforest", "gamlss", "bamlss"), 
#         col = c('forestgreen', 'red', 'blue', 'purple'), lty = 1, cex = 0.7)
#}



#plot_rmse(simres, type = "exp")
#plot_rmse(simres, type = "par")
#plot_ll(simres)







