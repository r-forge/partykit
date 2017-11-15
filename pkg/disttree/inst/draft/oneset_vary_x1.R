################################
### simulation function for one data set
sim_oneset <- function(parfun,
                       nobs = 500,
                       seedconst = 7, ntree = 100,
                       formula = y~x1+x2+x3+x4+x5+x6, 
                       tree_minsplit = 25, tree_minbucket = 10, tree_mincrit = 0.95, 
                       forest_minsplit = 7, forest_minbucket = 4, forest_mincrit = 0, 
                       forest_mtry = 3,
                       type.tree = "ctree",
                       censNO = TRUE,
                       pred_fix = FALSE,
                       pred_fix_x1 = -0.5,
                       pred_fix_x2 = -0.5,
                       pred_fix_x4 = 0,
                       nrunif = 3,
                       nrcovform = 6)
{
  ### preliminaries
  library("disttree")
  library("gamlss")
  library("lattice")
  library("crch")
  library("latex2exp")
  library("parallel")
  library("gamlss.cens")
  library("Formula")  # FIX ME: should be imported in disttree
  library("survival")
  gen.cens(NO, type = "left")
  
  cl <- match.call()
  
  
  
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
  #        noise_sd     ..... optional parameter to set the standard deviation of the noise variable
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
    x1 <- runif(n,-1,1)
    x2 <- runif(n,-1,1)
    x3 <- runif(n,-1,1)
    if(nrunif == 6){
      x4 <- runif(n,-1,1)
      x5 <- runif(n,-1,1)
      x6 <- runif(n,-1,1)
    }
    if(nrunif == 3){
      x4 <- rbinom(n,1,0.5)
      x5 <- rbinom(n,1,0.5)
      x6 <- rbinom(n,1,0.5)
    }
    x <- cbind(x1,x2,x3,x4,x5,x6)
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
      sp = split.matrix[,4]    # split points
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
        colnames(d) <- c("y", paste0("x", c(1:NCOL(x))), "mu", "sigma")
      } else {
        y <- rfun(n, dpar)
        d <- as.data.frame(cbind(y, x, dpar))
        colnames(d) <- c("y", paste0("x", c(1:NCOL(x))), "lambda")
      }
    }
    
    
    if(family$family.name == "censored Normal Distribution") {
      #d$ystar <- d$y
      d$y <- pmax(d$y, 0)
    }
    
    return(d)
  }
  
  
  
  
  # generate data
  set.seed(seedconst)
  
  data <- dgp(nobs, family = family, round.sp = 4, fun = parfun)
  learndata <- data[,-c(NCOL(data)-1, NCOL(data))]
  
  if(pred_fix){
    # data sets with fixed covariates except for one of them (x1, x2 or x4)
    data_x1 <- data_x2 <- data_x4 <- learndata[,-1]
    data_x1[,c("x3","x5","x6")] <-
      data_x2[,c("x3","x5","x6")] <-
      data_x4[,c("x3","x5","x6")] <- 0
    data_x1[,c("x2", "x4")] <- cbind(rep.int(pred_fix_x2, NROW(learndata)),
                                     rep.int(pred_fix_x4, NROW(learndata)))
    data_x2[,c("x1", "x4")] <- cbind(rep.int(pred_fix_x1, NROW(learndata)),
                                     rep.int(pred_fix_x4, NROW(learndata)))
    data_x4[,c("x1", "x2")] <- cbind(rep.int(pred_fix_x1, NROW(learndata)),
                                     rep.int(pred_fix_x2, NROW(learndata)))
  }  
  
  
  pred <- list()
  pred_x1 <- list()
  pred_x2 <- list()
  pred_x4 <- list()
  
  
  # eval disttree
  {
    if(type.tree == "mob"){
      dt <- disttree(formula, data=learndata, family=family, type.tree = "mob", 
                     control = mob_control(restart = FALSE, numsplit = "center", minbucket = tree_minbucket,
                                           alpha = 1-tree_mincrit, minsplit = tree_minsplit))
    }
    if(type.tree == "ctree"){
      dt <- disttree(formula, data=learndata, family=family, type.tree = "ctree", 
                     control = ctree_control(teststat = "quad", testtype = "Bonferroni", intersplit = TRUE,
                                             minbucket = tree_minbucket, mincriterion = tree_mincrit, minsplit = tree_minsplit))
    }
    
    
    pdt <- predict(dt, type = "parameter")
    dt_mu <- pdt$mu
    dt_sigma <- pdt$sigma
    
    if(censNO){
      dt_exp <- pnorm(dt_mu/dt_sigma) * (dt_mu + dt_sigma * (dnorm(dt_mu/dt_sigma) / pnorm(dt_mu/dt_sigma)))
      ## FIX ME: calculation of dt_exp for sigma set to 0.0001
      # idea: 
      if(any(is.na(dt_exp))){
        dt_exp[dt_sigma <= 0.0002] <- pmax(0, dt_mu[dt_sigma <= 0.0002])    
      }
    } else dt_exp <- dt_mu
    
    dt_pred <- cbind(data[,c("y","mu","sigma")], 
                     dt_exp, dt_mu, dt_sigma)
    colnames(dt_pred) <- c("y","true.mu","true.sigma", 
                           "fitted.exp.dt", "fitted.mu.dt","fitted.sigma.dt")
    
    pred <- c(pred, list(dt_pred))
    
    
    if(pred_fix){
      dt_pred_x1 <- cbind(data[,c("y","x1")], parfun(data_x1),
                        predict(dt, newdata = data_x1, type = "parameter"))
      colnames(dt_pred_x1) <- c("y","x1","true.mu","true.sigma", 
                              "fitted.mu.dt","fitted.sigma.dt")
      
      dt_pred_x2 <- cbind(data[,c("y","x2")], parfun(data_x2),
                              predict(dt, newdata = data_x2, type = "parameter"))
      colnames(dt_pred_x2) <- c("y","x2","true.mu","true.sigma", 
                                    "fitted.mu.dt","fitted.sigma.dt")
      
      dt_pred_x4 <- cbind(data[,c("y","x4")], parfun(data_x4),
                          predict(dt, newdata = data_x4, type = "parameter"))
      colnames(dt_pred_x4) <- c("y","x4","true.mu","true.sigma", 
                                "fitted.mu.dt","fitted.sigma.dt")
      
      pred_x1 <- c(pred_x1, list(dt_pred_x1))
      pred_x2 <- c(pred_x2, list(dt_pred_x2))
      pred_x4 <- c(pred_x4, list(dt_pred_x4))
    }
  }
  
  
  # eval distforest
  {
    if(type.tree == "mob"){
      df <- distforest(formula, data=learndata, family=family, type.tree = "mob", ntree = ntree,
                       control = mob_control(restart = FALSE, numsplit = "center", 
                                             alpha = 1-forest_mincrit, minsplit = forest_minsplit,
                                             minbucket = forest_minbucket, mtry = forest_mtry))
    }
    if(type.tree == "ctree"){
      df <- distforest(formula, data=learndata, family=family, type.tree = "ctree", ntree = ntree,
                       control = ctree_control(teststat = "quad", testtype = "Univ", intersplit = TRUE,
                                               mincriterion = forest_mincrit, minsplit = forest_minsplit,
                                               minbucket = forest_minbucket, mtry = forest_mtry))
    }
    
    pdf <- predict(df, type = "parameter")
    df_mu <- pdf$mu
    df_sigma <- pdf$sigma
    
    if(censNO){
      df_exp <- pnorm(df_mu/df_sigma) * (df_mu + df_sigma * (dnorm(df_mu/df_sigma) / pnorm(df_mu/df_sigma)))
      ## FIX ME: calculation of df_exp for sigma set to 0.0001
      # idea: 
      if(any(is.na(df_exp))){
        df_exp[df_sigma <= 0.0002] <- pmax(0, df_mu[df_sigma <= 0.0002])    
      }
    } else df_exp <- df_mu
    
    df_pred <- cbind(data[,c("y","mu","sigma")], 
                     df_exp, df_mu, df_sigma)
    colnames(df_pred) <- c("y","true.mu","true.sigma", 
                           "fitted.exp.df", "fitted.mu.df","fitted.sigma.df")
    
    pred <- c(pred, list(df_pred))
    
    
    if(pred_fix){
      df_pred_x1 <- cbind(data[,c("y","x1")], parfun(data_x1),
                              predict(df, newdata = data_x1, type = "parameter"))
      colnames(df_pred_x1) <- c("y","x1","true.mu","true.sigma", 
                                    "fitted.mu.df","fitted.sigma.df")
      
      df_pred_x2 <- cbind(data[,c("y","x2")], parfun(data_x2),
                              predict(df, newdata = data_x2, type = "parameter"))
      colnames(df_pred_x2) <- c("y","x2","true.mu","true.sigma", 
                                    "fitted.mu.df","fitted.sigma.df")
      
      df_pred_x4 <- cbind(data[,c("y","x4")], parfun(data_x4),
                          predict(df, newdata = data_x4, type = "parameter"))
      colnames(df_pred_x4) <- c("y","x4","true.mu","true.sigma", 
                                "fitted.mu.df","fitted.sigma.df")
      
      pred_x1 <- c(pred_x1, list(df_pred_x1))
      pred_x2 <- c(pred_x2, list(df_pred_x2))
      pred_x4 <- c(pred_x4, list(df_pred_x4))
    }
  }
  
  
  # eval gamlss
  {
    
    if(nrcovform == 6){    
      mu.formula <- y~pb(x1)+pb(x2)+pb(x3)+pb(x4)+pb(x5)+pb(x6)
      sigma.formula <- ~pb(x1)+pb(x2)+pb(x3)+pb(x4)+pb(x5)+pb(x6)
    }
    if(nrcovform == 3){    
      mu.formula <- y~pb(x1)+pb(x2)+pb(x3)
      sigma.formula <- ~pb(x1)+pb(x2)+pb(x3)
    }
    if(nrcovform == 4){    
      mu.formula <- y~pb(x1)+pb(x2)+pb(x4)+pb(x5)
      sigma.formula <- ~pb(x1)+pb(x2)+pb(x4)+pb(x5)
    }
    if(censNO){
      g_learndata <- learndata
      g_learndata$y <- Surv(g_learndata$y, g_learndata$y>0, type="left")
      g <- gamlss(formula = mu.formula, sigma.formula = sigma.formula, data = g_learndata, family = cens("NO", type = "left")())
    } else {
      g <- gamlss(formula = mu.formula, sigma.formula = sigma.formula, data = learndata, family = NO())
    }
    
    g_mu <- predict(g, what = "mu", type = "response", data = g_learndata)
    g_sigma <- predict(g, what = "sigma", type = "response", data = g_learndata)
    
    if(censNO){
      #calculate expected value for censored data
      g_exp <- pnorm(g_mu/g_sigma) * (g_mu + g_sigma * (dnorm(g_mu/g_sigma) / pnorm(g_mu/g_sigma)))
    } else {  
      # for non censored normal data the expected value is the location paramter
      g_exp <- g_mu
    }
    
    g_pred <- cbind(data[,c("y","mu","sigma")], 
                     g_exp, g_mu, g_sigma)
    colnames(g_pred) <- c("y","true.mu","true.sigma", 
                           "fitted.exp.g", "fitted.mu.g","fitted.sigma.g")
    
    pred <- c(pred, list(g_pred))
    
    
    if(pred_fix){
      g_pred_x1 <- cbind(data[,c("y","x1")], parfun(data_x1),
                             predict(g, what = "mu", type = "response", newdata = data_x1, data = g_learndata),
                             predict(g, what = "sigma", type = "response", newdata = data_x1, data = g_learndata))
      colnames(g_pred_x1) <- c("y","x1","true.mu","true.sigma", 
                                    "fitted.mu.g","fitted.sigma.g")
      
      g_pred_x2 <- cbind(data[,c("y","x2")], parfun(data_x2),
                         predict(g, what = "mu", type = "response", newdata = data_x2, data = g_learndata),
                         predict(g, what = "sigma", type = "response", newdata = data_x2, data = g_learndata))
      colnames(g_pred_x2) <- c("y","x2","true.mu","true.sigma", 
                               "fitted.mu.g","fitted.sigma.g")
      
      g_pred_x4 <- cbind(data[,c("y","x4")], parfun(data_x4),
                             predict(g, what = "mu", type = "response", newdata = data_x4, data = g_learndata),
                             predict(g, what = "sigma", type = "response", newdata = data_x4, data = g_learndata))
      colnames(g_pred_x4) <- c("y","x4","true.mu","true.sigma", 
                                    "fitted.mu.g","fitted.sigma.g")
      
      pred_x1 <- c(pred_x1, list(g_pred_x1))
      pred_x2 <- c(pred_x2, list(g_pred_x2))
      pred_x4 <- c(pred_x4, list(g_pred_x4))
    }
  }
  
  
  res <- list()
  res$call <- cl
  res$data <- data
  res$fun <- parfun
  res$pred <- pred
  res$pred_x1 <- pred_x1
  res$pred_x2 <- pred_x2
  res$pred_x4 <- pred_x4
  names(res$pred) <- names(res$pred_x1) <- names(res$pred_x2) <- 
    names(res$pred_x4) <- c("dt", "df", "g")
  res$data_x1 <- data_x1
  res$data_x2 <- data_x2
  res$data_x4 <- data_x4
  res$family <- family
  
  #names(res) <- c("call", "data", "fun", "pred", "pred_x1", "pred_x2", "pred_x4")
  
  return(res)
}



################################
### plotting function for one data set
plot_oneset <- function(oneset,
                  only_dt = FALSE,
                  only_df = FALSE,
                  only_g = FALSE,
                  compare_mu = FALSE, 
                  w_sigma = FALSE,
                  compare_sigma_area = FALSE,
                  compare_sigma_line = FALSE,
                  add_dt = FALSE,
                  add_df = FALSE,
                  add_g = FALSE,
                  nomodel = FALSE,
                  ylim = NULL,
                  plot_legend = FALSE) 
{
  ## HCL palette
  pal <- hcl(c(10, 128, 260, 290, 30, 90, 180), 100, 50)
  names(pal) <- c("distforest", "disttree", "gamlss", "randomForest", "bamlss", "gamboostLSS", "cforest")
  
  pallight <- hcl(c(10, 128, 260, 290, 30, 90, 180), 100, 50, alpha = 0.25)
  names(pallight) <- c("distforest", "disttree", "gamlss", "randomForest", "bamlss", "gamboostLSS", "cforest")
  
  transpgrey <- rgb(0.190,0.190,0.190, alpha = 0.2)
  
  #learndata <- oneset$data[,-c(NCOL(oneset$data)-1, NCOL(oneset$data))]
  
  if(nomodel){
    plotdata <- oneset$data[,c("y","x1","mu","sigma")]
    colnames(plotdata) <- c("y","x1","true.mu","true.sigma")
    sp <- plotdata[order(plotdata["x1"]),]
    
    if(is.null(ylim)) ylim <- c(min(sp$true.mu - sp$true.sigma), max(sp$true.mu + sp$true.sigma))
    par(mar=c(5.1,4.1,4.1,3.1))
    plot(y = sp$y, x = sp$x, type = "p", col="grey", main = "True parameters", cex.main = 1.2, #xaxt="n", yaxt="n", 
         xlab = expression(x[1]), ylab = expression(mu(x[1])), ylim = ylim)
    lines(x = sp$x, y = sp$true.mu, type = "l", col = 'black')
    polygon(c(sp$x, rev(sp$x)), c(sp$true.mu + sp$true.sigma, rev(sp$true.mu - sp$true.sigma)),
            col = transpgrey, border = "transparent")
    legend("topleft", legend = expression(mu %+-% sigma), bty = "n")
  }
  
  # disttree
  if(only_dt){
    
    plotdata <- oneset$pred_x1$dt
    sp <- plotdata[order(plotdata["x1"]),]
    
    par(mar=c(5.1,4.1,4.1,3.1))
    if(is.null(ylim)) ylim <- c(min(plotdata[,c(1,3,5)] - sp$true.sigma), max(plotdata[,c(1,3,5)] + sp$true.sigma))
    plot(y = sp$y, x = sp$x1, type = "p", col="grey", main = "disttree", col.main = pal["disttree"], cex.main = 1.2, #xaxt="n", yaxt="n", 
         xlab = expression(x[1]), ylab = expression(mu(x[1])), ylim = ylim)
    lines(x = sp$x1, y = sp$true.mu, type = "l", col = 'black')
    lines(x = sp$x1, y = sp$true.mu + sp$true.sigma, type = "l", col = 'grey')
    lines(x = sp$x1, y = sp$true.mu - sp$true.sigma, type = "l", col = 'grey')
    lines(x = sp$x1, y = sp$fitted.mu.dt, type = "l", col = pal["disttree"], lwd = 2.5)
    polygon(c(sp$x1, rev(sp$x1)), c(sp$fitted.mu.dt + sp$fitted.sigma.dt, rev(sp$fitted.mu.dt - sp$fitted.sigma.dt)),
            col = pallight["disttree"], border = "transparent")
    legend("topleft", legend = expression(mu %+-% sigma), bty = "n")
  }
  
  # distforest
  if(only_df){
    
    plotdata <- oneset$pred_x1$df
    sp <- plotdata[order(plotdata["x1"]),]
    if(is.null(ylim)) ylim <- c(min(plotdata[,c(1,3,5)] - sp$true.sigma), max(plotdata[,c(1,3,5)] + sp$true.sigma))
    par(mar=c(5.1,4.1,4.1,3.1))
    plot(y = sp$y, x = sp$x1, type = "p", col="grey", main = "distforest", col.main = pal["distforest"], cex.main = 1.2, #xaxt="n", yaxt="n", 
         xlab = expression(x[1]), ylab = expression(mu(x[1])), ylim = ylim)
    lines(x = sp$x1, y = sp$true.mu, type = "l", col = 'black')
    lines(x = sp$x1, y = sp$true.mu + sp$true.sigma, type = "l", col = 'grey')
    lines(x = sp$x1, y = sp$true.mu - sp$true.sigma, type = "l", col = 'grey')
    lines(x = sp$x1, y = sp$fitted.mu.df, type = "l", col = pal["distforest"], lwd = 2.5)
    polygon(c(sp$x1, rev(sp$x1)), c(sp$fitted.mu.df + sp$fitted.sigma.df, rev(sp$fitted.mu.df - sp$fitted.sigma.df)),
            col = pallight["distforest"], border = "transparent")
    legend("topleft", legend = expression(mu %+-% sigma), bty = "n")
  }
  
  # gamlss
  if(only_g){

    plotdata <- oneset$pred_x1$g
    sp <- plotdata[order(plotdata["x1"]),]
    if(is.null(ylim)) ylim <- c(min(plotdata[,c(1,3,5)] - sp$true.sigma), max(plotdata[,c(1,3,5)] + sp$true.sigma))
    par(mar=c(5.1,4.1,4.1,3.1))
    plot(y = sp$y, x = sp$x1, type = "p", col="grey",  main = "gamlss", col.main = pal["gamlss"], cex.main = 1.2, #xaxt="n", yaxt="n", 
         xlab = expression(x[1]), ylab = expression(mu(x[1])), ylim = ylim)
    lines(x = sp$x1, y = sp$true.mu, type = "l", col = 'black')
    lines(x = sp$x1, y = sp$true.mu + sp$true.sigma, type = "l", col = 'grey')
    lines(x = sp$x1, y = sp$true.mu - sp$true.sigma, type = "l", col = 'grey')
    lines(x = sp$x1, y = sp$fitted.mu.g, type = "l", col = pal["gamlss"], lwd = 2.5)
    polygon(c(sp$x1, rev(sp$x1)), c(sp$fitted.mu.g + sp$fitted.sigma.g, rev(sp$fitted.mu.g - sp$fitted.sigma.g)),
            col = pallight["gamlss"], border = "transparent")
    legend("topleft", legend = expression(mu %+-% sigma), bty = "n")
    
  }

  
  
  # compare location parameter mu
  if(compare_mu){
    
    plotdata <- cbind(oneset$data[,c("y","x1")], oneset$fun(oneset$data_x1))
    coln <- c("y","x1","true.mu","true.sigma")
    methodnames <- NULL
    
    if(add_dt) {
      plotdata <- cbind(plotdata, oneset$pred_x1$dt[,c("fitted.mu.dt","fitted.sigma.dt")])
      coln <- c(coln, "fitted.mu.dt","fitted.sigma.dt")
      methodnames <- c(methodnames, "disttree")
    }
    if(add_df) {
      plotdata <- cbind(plotdata, oneset$pred_x1$df[,c("fitted.mu.df","fitted.sigma.df")])
      coln <- c(coln, "fitted.mu.df","fitted.sigma.df")
      methodnames <- c(methodnames, "distforest")
    }
    if(add_g) {
      plotdata <- cbind(plotdata, oneset$pred_x1$g[,c("fitted.mu.g","fitted.sigma.g")])
      coln <- c(coln, "fitted.mu.g","fitted.sigma.g")
      methodnames <- c(methodnames, "gamlss")
    }
    
    colnames(plotdata) <- coln  
    sp <- plotdata[order(plotdata["x1"]),]
    
    ymin <- min(sp[,(c(2:(length(coln)/2))*2)-1]) - max(sp$true.sigma)-1
    ymax <- max(sp[,(c(2:(length(coln)/2))*2)-1]) + max(sp$true.sigma)+1
    if(is.null(ylim)) ylim <- c(ymin, ymax)
    
    
    par(mar=c(5.1,4.1,4.1,3.1))
    plot(y = sp$y, x = sp$x1, type = "p", col="grey", #xaxt="n", yaxt="n", 
         xlab = expression(x[1]), ylab = expression(mu(x[1])), ylim = ylim)
    if(w_sigma){
      lines(x = sp$x1, y = sp$true.mu + sp$true.sigma, type = "l", col = 'black', lty = 2)
      lines(x = sp$x1, y = sp$true.mu - sp$true.sigma, type = "l", col = 'black', lty = 2)
    }
      
    
    lines(x = sp$x1, y = sp$true.mu, type = "l", col = 'black')
    if(add_dt) {
      lines(x = sp$x1, y = sp$fitted.mu.dt, type = "l", col = pal["disttree"], lwd = 2)
      if(w_sigma){
        polygon(c(sp$x1, rev(sp$x1)), c(sp$fitted.mu.dt + sp$fitted.sigma.dt, rev(sp$fitted.mu.dt - sp$fitted.sigma.dt)),
                          col = pallight["disttree"], border = "transparent")
      }
    }
    if(add_df) {
      lines(x = sp$x1, y = sp$fitted.mu.df, type = "l", col = pal["distforest"], lwd = 2)
      if(w_sigma){
        polygon(c(sp$x1, rev(sp$x1)), c(sp$fitted.mu.df + sp$fitted.sigma.df, rev(sp$fitted.mu.df - sp$fitted.sigma.df)),
                col = pallight["distforest"], border = "transparent")
      }}
    if(add_g) {
      lines(x = sp$x1, y = sp$fitted.mu.g, type = "l", col = pal["gamlss"], lwd = 2)
      if(w_sigma){
        polygon(c(sp$x1, rev(sp$x1)), c(sp$fitted.mu.g + sp$fitted.sigma.g, rev(sp$fitted.mu.g - sp$fitted.sigma.g)),
                col = pallight["gamlss"], border = "transparent")
      }
    }
    if(plot_legend) legend("topleft", legend =  methodnames, col = pal[methodnames], lty = 1, lwd = 2, bty = "n")
  }
  
  # variance
  if(compare_sigma_area){
    
    plotdata <- cbind(oneset$data[,c("y","x1")], oneset$fun(oneset$data_x1))
    coln <- c("y","x1","true.mu","true.sigma")
    methodnames <- NULL
    
    if(add_dt) {
      plotdata <- cbind(plotdata, oneset$pred_x1$dt[,c("fitted.mu.dt","fitted.sigma.dt")])
      coln <- c(coln, "fitted.mu.dt","fitted.sigma.dt")
      methodnames <- c(methodnames, "disttree")
    }
    if(add_df) {
      plotdata <- cbind(plotdata, oneset$pred_x1$df[,c("fitted.mu.df","fitted.sigma.df")])
      coln <- c(coln, "fitted.mu.df","fitted.sigma.df")
      methodnames <- c(methodnames, "distforest")
    }
    if(add_g) {
      plotdata <- cbind(plotdata, oneset$pred_x1$g[,c("fitted.mu.g","fitted.sigma.g")])
      coln <- c(coln, "fitted.mu.g","fitted.sigma.g")
      methodnames <- c(methodnames, "gamlss")
    }

    colnames(plotdata) <- coln  
    sp <- plotdata[order(plotdata["x1"]),]
    
    ymin <- min(sp[,(c(2:(length(coln)/2))*2)-1]) - max(sp$true.sigma) - 1
    ymax <- max(sp[,(c(2:(length(coln)/2))*2)-1]) + max(sp$true.sigma) + 1
    if(is.null(ylim)) ylim <- c(ymin, ymax)
    
    par(mar=c(5.1,4.1,4.1,3.1))
    plot(y = sp$y, x = sp$x1, type = "p", col="grey", #xaxt="n", yaxt="n", 
         xlab = expression(x[1]), ylab = expression(mu(x[1])), ylim = ylim)
    
    lines(x = sp$x1, y = sp$true.mu + sp$true.sigma, type = "l", col = 'black')
    lines(x = sp$x1, y = sp$true.mu - sp$true.sigma, type = "l", col = 'black')
    
    if(add_dt) polygon(c(sp$x1, rev(sp$x1)), c(sp$fitted.mu.dt + sp$fitted.sigma.dt, rev(sp$fitted.mu.dt - sp$fitted.sigma.dt)),
            col = pallight["disttree"], border = "transparent")
    if(add_df) polygon(c(sp$x1, rev(sp$x1)), c(sp$fitted.mu.df + sp$fitted.sigma.df, rev(sp$fitted.mu.df - sp$fitted.sigma.df)),
            col = pallight["distforest"], border = "transparent")
    if(add_g) polygon(c(sp$x1, rev(sp$x1)), c(sp$fitted.mu.g + sp$fitted.sigma.g, rev(sp$fitted.mu.g - sp$fitted.sigma.g)),
            col = pallight["gamlss"], border = "transparent")

    if(plot_legend) legend("topleft", legend =  methodnames, col = pal[methodnames], lty = 1, lwd = 2, bty = "n")
  }

  # variance
  if(compare_sigma_line){
    
    plotdata <- cbind(oneset$data[,c("y","x1")], oneset$fun(oneset$data_x1))
    coln <- c("y","x1","true.mu","true.sigma")
    methodnames <- NULL
    
    if(add_dt) {
      plotdata <- cbind(plotdata, oneset$pred_x1$dt[,c("fitted.mu.dt","fitted.sigma.dt")])
      coln <- c(coln, "fitted.mu.dt","fitted.sigma.dt")
      methodnames <- c(methodnames, "disttree")
    }
    if(add_df) {
      plotdata <- cbind(plotdata, oneset$pred_x1$df[,c("fitted.mu.df","fitted.sigma.df")])
      coln <- c(coln, "fitted.mu.df","fitted.sigma.df")
      methodnames <- c(methodnames, "distforest")
    }
    if(add_g) {
      plotdata <- cbind(plotdata, oneset$pred_x1$g[,c("fitted.mu.g","fitted.sigma.g")])
      coln <- c(coln, "fitted.mu.g","fitted.sigma.g")
      methodnames <- c(methodnames, "gamlss")
    }

    
    colnames(plotdata) <- coln  
    sp <- plotdata[order(plotdata["x1"]),]
    
    par(mar=c(5.1,4.1,4.1,3.1))
    ymin <- min(sp[,c(2:(length(coln)/2))*2])
    ymax <- max(sp[,c(2:(length(coln)/2))*2])
    if(is.null(ylim)) ylim <- c(ymin, ymax)
    
    plot(x = sp$x1, y = sp$true.sigma, type = "l", col = 'black', #xaxt="n", yaxt="n", 
         xlab = expression(x[1]), ylab = expression(mu(x[1])), ylim = c(ymin, ymax))
    
    if(add_dt) lines(x = sp$x1, y = sp$fitted.sigma.dt, type = "l", col = pal["disttree"], lwd = 2)
    if(add_df) lines(x = sp$x1, y = sp$fitted.sigma.df, type = "l", col = pal["distforest"], lwd = 2)
    if(add_g) lines(x = sp$x1, y = sp$fitted.sigma.g, type = "l", col = pal["gamlss"], lwd = 2)
    
    if(plot_legend) legend("topleft", legend =  methodnames, col = pal[methodnames], lty = 1, lwd = 2, bty = "n")
  }
  
}

















##################################################################
# produce and save plots
if(FALSE){
  kappa <- 1
  fun <- function(x){
    mu <- 0.2 + 
      (7 * (exp(-(3*x[,1]-1.2)^(2*kappa)))) * as.numeric(x[,1] < 0.4) + 
      (7 * (exp(-(3*x[,1]-1.2)^(2*kappa))) * (1/2) + 7/2) * as.numeric(x[,1]>= 0.4) #+
    #(3 * (exp(-(3*x[,1]-1.2)^(2*kappa)))) +
    (3 * (exp(-(3*x[,1]-1.2)^(2*kappa)))) * x[,4]
    sigma <- 0.5 + 
      mu/5 + 
      #(x[,1]>0.4)  #*x[,4]#+
      2*(x[,1] < -0.7) 
    #plogis((kappa^(1.8)) * 4 * (x[,2] + 0.7)) #* x[,4]
    #4*(1-plogis((kappa^(1.8)) * 2 * (x[,1] + 0.5))) * x[,4]
    #5 * abs(x[,1]) * x[,4]
    par <- cbind(mu, sigma)
    return(par)
  }
  
  
  
  library("gamlss.cens")
  gen.cens("NO", type = "left")
  
  oneset1 <- sim_oneset(parfun = fun, nobs = 500, seedconst = 74, ntree = 100,
                        formula = y~x1+x2+x3+x4+x5+x6, 
                        #formula = y~x1+x2+x4+x5,
                        #formula = y~x1+x2+x3,
                        tree_minsplit = 20, tree_minbucket = 7, tree_mincrit = 0.95, 
                        forest_minsplit = 20, forest_minbucket = 7, forest_mincrit = 0, 
                        forest_mtry = 3, type.tree = "ctree", 
                        censNO = TRUE, 
                        pred_fix = TRUE,
                        pred_fix_x1 = 0.9,
                        pred_fix_x2 = 0.9,
                        pred_fix_x4 = 1,
                        gamboost_cvr = FALSE,
                        eval_disttree = TRUE,
                        eval_distforest = TRUE,
                        eval_gamlss = TRUE,
                        nrunif = 3,
                        nrcovform = 6)
  
  
  kappa <- 10
  fun <- function(x){
    mu <- 0.2 + 
      (7 * (exp(-(3*x[,1]-1.2)^(2*kappa)))) * as.numeric(x[,1] < 0.4) + 
      (7 * (exp(-(3*x[,1]-1.2)^(2*kappa))) * (1/2) + 7/2) * as.numeric(x[,1]>= 0.4) #+
    #(3 * (exp(-(3*x[,1]-1.2)^(2*kappa)))) +
    (3 * (exp(-(3*x[,1]-1.2)^(2*kappa)))) * x[,4]
    sigma <- 0.5 + 
      mu/5 + 
      #(x[,1]>0.4)  #*x[,4]#+
      2*(x[,1] < -0.7) 
    #plogis((kappa^(1.8)) * 4 * (x[,2] + 0.7)) #* x[,4]
    #4*(1-plogis((kappa^(1.8)) * 2 * (x[,1] + 0.5))) * x[,4]
    #5 * abs(x[,1]) * x[,4]
    par <- cbind(mu, sigma)
    return(par)
  }
  
  
  
  library("gamlss.cens")
  gen.cens("NO", type = "left")
  
  oneset10 <- sim_oneset(parfun = fun, nobs = 500, seedconst = 74, ntree = 100,
                         formula = y~x1+x2+x3+x4+x5+x6, 
                         #formula = y~x1+x2+x4+x5,
                         #formula = y~x1+x2+x3,
                         tree_minsplit = 20, tree_minbucket = 7, tree_mincrit = 0.95, 
                         forest_minsplit = 20, forest_minbucket = 7, forest_mincrit = 0, 
                         forest_mtry = 3, type.tree = "ctree", 
                         censNO = TRUE, 
                         pred_fix = TRUE,
                         pred_fix_x1 = 0.9,
                         pred_fix_x2 = 0.9,
                         pred_fix_x4 = 1,
                         gamboost_cvr = FALSE,
                         eval_disttree = TRUE,
                         eval_distforest = TRUE,
                         eval_gamlss = TRUE,
                         nrunif = 3,
                         nrcovform = 6)
  
  
  
  
  jpeg(file = "plot_oneset1_mu.jpeg")
  plot_oneset(oneset1, only_dt = FALSE, only_df = FALSE, only_g = FALSE,
              compare_mu = TRUE, 
              add_dt = TRUE,
              add_df = TRUE,
              add_g = TRUE,
              nomodel = FALSE)
  dev.off()
  
  jpeg(file = "plot_oneset1_sigma.jpeg")
  plot1sigma <- plot_oneset(oneset1, only_dt = FALSE, only_df = FALSE, only_g = FALSE,
                            compare_mu = FALSE,
                            compare_sigma_line = TRUE, 
                            add_dt = TRUE,
                            add_df = TRUE,
                            add_g = TRUE,
                            nomodel = FALSE)
  dev.off()
  
  jpeg(file = "plot_oneset10_mu.jpeg")
  plot_oneset(oneset10, only_dt = FALSE, only_df = FALSE, only_g = FALSE,
              compare_mu = TRUE, 
              add_dt = TRUE,
              add_df = TRUE,
              add_g = TRUE,
              nomodel = FALSE)
  dev.off()
  
  jpeg(file = "plot_oneset10_sigma.jpeg")
  plot1sigma <- plot_oneset(oneset10, only_dt = FALSE, only_df = FALSE, only_g = FALSE,
                            compare_mu = FALSE,
                            compare_sigma_line = TRUE, 
                            add_dt = TRUE,
                            add_df = TRUE,
                            add_g = TRUE,
                            nomodel = FALSE)
  dev.off()
}