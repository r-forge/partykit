##########
# compare to formula=y~x1+x2: due to large scale of x2 all models get very wiggly on 1-dim plot
# no structure recognizeable


sim_onecov <- function(kappa = 1, nobs = 400,
                       seedconst = 7, ntree = 100,
                       formula = y~x1, 
                       tree_minsplit = 14, tree_mincrit = 0.95,
                       forest_minsplit = 7, forest_mincrit = 0, 
                       type.tree = "ctree",
                       noise_sd = 0,
                       censNO = TRUE,
                       censpar = 0.7,
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
                  round.sp = 3,
                  noise_sd = 0)
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
    
    if(noise_sd > 0) d$y <- d$y + rnorm(n, mean = 0, sd = noise_sd)
    
    if(family$family.name == "censored Normal Distribution") {
      #d$ystar <- d$y
      d$y <- pmax(d$y, 0)
    }
    
    return(d)
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
  
  
  if(fix.mu){
    if(fix.sigma){
      fun <- function(x){
        mu <- 0
        sigma <- 3
        if(censNO) sigma[x[,1]<0] <- 0.00001
        par <- cbind(mu, sigma)
        return(par)
      }
    } else {
      fun <- function(x){
        mu <- 0
        sigma <- 1+3*abs(x[,1])
        if(censNO) sigma[x[,1]<0] <- 0.00001
        par <- cbind(mu, sigma)
        return(par)
      }
    }
  } else {
    if(fix.sigma){
      fun <- function(x){
        mu <- 4 + 8 * (exp(-(4*x[,1]-censpar)^(2*kappa)))
        #mu <- 2 + 10 * (exp(-(4*x[,1]-2)^(2*kappa)))
        if(censNO) mu[x[,1]<0] <- 0
        sigma <- 3
        if(censNO) sigma[x[,1]<0] <- 0.00001
        par <- cbind(mu, sigma)
        return(par)
      }
    } else {
      if(mu.sigma.interaction){
        fun <- function(x){
          mu <- 4 + 8 * (exp(-(4*x[,1]-censpar)^(2*kappa)))
          #mu <- 2 + 10 * (exp(-(4*x[,1]-2)^(2*kappa)))
          if(censNO) mu[x[,1]<0] <- 0
          sigma <- 2 + mu/4
          if(censNO) sigma[x[,1]<0] <- 0.00001
          par <- cbind(mu, sigma)
          return(par)
        }
      } else {
        fun <- function(x){
          mu <- 4 + 8 * (exp(-(4*x[,1]-censpar)^(2*kappa)))
          #mu <- 2 + 10 * (exp(-(4*x[,1]-2)^(2*kappa)))
          if(censNO) mu[x[,1]<0] <- 0
          sigma <- 1+3*abs(x[,1])
          if(censNO) sigma[x[,1]<0] <- 0.00001
          par <- cbind(mu, sigma)
          return(par)
        }
      }
    }
  }
  
  
  
  
  
  
  # generate data
  set.seed(seedconst)
  
  learndata <- dgp(nobs, family = family, round.sp = 4, fun = fun, noise_sd = noise_sd)
  
  res_onecov <- list()
  res_onecov$call <- cl
  res_onecov$learndata <- learndata
  res_onecov$fun <- fun
  
  listnames <- c("call", "learndata", "fun")
  
  if(eval_disttree){
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
    res_onecov <- c(res_onecov, list(dt))
    listnames <- c(listnames, "dt")
  }
  
  if(eval_distforest){
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
    res_onecov <- c(res_onecov, list(df))
    listnames <- c(listnames, "df")
  }
  
  if(eval_gamlss){
    
    formula_rh <- formula[[3]]
    if("x1" %in% as.character(formula_rh)){
      mu.formula <- y~pb(x1)
      sigma.formula <- ~pb(x1)
      if("x2" %in% as.character(formula_rh)){
        mu.formula <- y~pb(x1)+pb(x2)
        sigma.formula <- ~pb(x1)+pb(x2)
      }
    } 
    
    if(censNO){
      g_learndata <- learndata
      g_learndata$y <- Surv(g_learndata$y, g_learndata$y>0, type="left")
      g <- gamlss(formula = mu.formula, sigma.formula = sigma.formula, data = g_learndata, family = NOlc())
    } else {
      g <- gamlss(formula = mu.formula, sigma.formula = sigma.formula, data = learndata, family = NO())
    }
    res_onecov <- c(res_onecov, list(g))
    listnames <- c(listnames, "g")
  }
  
  if(eval_bamlss){
    
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
    
    if(censNO) {
      ## FIX ME: set up for bamlss
      #learndata$y <- Surv(learndata$y, learndata$y>0, type="left")
      b <- bamlss(bamlss.formula(list(mu.formula, sigma.formula), family = cnorm_bamlss), data = learndata)
    } else {
      b <- bamlss(bamlss.formula(list(mu.formula, sigma.formula), family = "gaussian"), data = learndata)
    }
    res_onecov <- c(res_onecov, list(b))
    listnames <- c(listnames, "b")
  }
  
  if(eval_gamboostLSS){
    #for method="noncyclic"
    grid <- make.grid(max = 300, min = 10, length.out = 10)
    
    
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
    
    
    if(censNO) {
      g_learndata <- learndata
      g_learndata$y <- Surv(g_learndata$y, g_learndata$y>0, type="left")
      
      gb <- gamboostLSS(formula = list(mu = mu.formula, sigma =sigma.formula), data = g_learndata, 
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
    res_onecov <- c(res_onecov, list(gb))
    listnames <- c(listnames, "gb")
  }
  
  if(eval_randomForest){
    rf <- randomForest(formula, data = learndata, ntree = ntree, 
                       nodesize = forest_minsplit, 
                       keep.inbag = TRUE, replace = FALSE)
    res_onecov <- c(res_onecov, list(rf))
    listnames <- c(listnames, "rf")
  }
  
  if(eval_cforest){
    cf <- cforest(formula, data = learndata, ntree = ntree,
                  control = ctree_control(teststat = "quad", testtype = "Univ", 
                                          mincriterion = forest_mincrit, minsplit = forest_minsplit, 
                                          intersplit = TRUE))
    res_onecov <- c(res_onecov, list(cf))
    listnames <- c(listnames, "cf")
  }
  
  names(res_onecov) <- listnames
  return(res_onecov)
}







##########################3
# test
if(FALSE){
  sim_onecov_test <- sim_onecov(kappa = 1, nobs = 400,
                         seedconst = 7, ntree = 100,
                         formula = y~x1, 
                         tree_minsplit = 14, tree_mincrit = 0.95,
                         forest_minsplit = 7, forest_mincrit = 0, 
                         type.tree = "ctree",
                         noise_sd = 1,
                         censNO = TRUE,
                         censpar = 2, # = 0.7
                         fix.mu = FALSE,
                         fix.sigma = FALSE,
                         mu.sigma.interaction = FALSE,
                         gamboost_cvr = FALSE,
                         eval_disttree = TRUE,
                         eval_distforest = TRUE,
                         eval_gamlss = TRUE,
                         eval_bamlss = TRUE,
                         eval_gamboostLSS = TRUE,
                         eval_randomForest = TRUE,
                         eval_cforest = TRUE)
}


### plotting function for one data set
plot_onecov <- function(learndata, 
                  dt=NULL, df=NULL, g=NULL, gb=NULL, b=NULL, rf=NULL, cf=NULL, 
                  only_dt = FALSE,
                  only_df = FALSE,
                  only_g = FALSE,
                  only_gb = FALSE,
                  only_b = FALSE,
                  only_rf = FALSE,
                  only_cf = FALSE,
                  compare_mu=FALSE, 
                  compare_sigma_area=FALSE,
                  compare_sigma_line=FALSE,
                  add_dt = FALSE,
                  add_df = FALSE,
                  add_g = FALSE,
                  add_gb = FALSE,
                  add_b = FALSE,
                  add_rf = FALSE,
                  add_cf = FALSE,
                  nomodel = FALSE) 
{
  ## HCL palette
  pal <- hcl(c(10, 128, 260, 290, 30, 90, 180), 100, 50)
  names(pal) <- c("forest", "tree", "gamlss", "randomForest", "bamlss", "gamboostLSS", "cforest")
  
  pallight <- hcl(c(10, 128, 260, 290, 30, 90, 180), 100, 50, alpha = 0.25)
  names(pallight) <- c("forest", "tree", "gamlss", "randomForest", "bamlss", "gamboostLSS", "cforest")
  
  transpgrey <- rgb(0.190,0.190,0.190, alpha = 0.2)
  
  
  
  ## function to estimate standard deviation of cforest for a new observation
  cf_getsd <- function(cf, newdata = NULL){
    
    ntree <- cf$info$call$ntree
    
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
  
  
  # if(only_rf | (add_rf & (compare_sigma_area | compare_sigma_line)))
  if(only_rf | add_rf) {
    rf_sd <- numeric(length(NROW(learndata)))
    for(l in 1:NROW(learndata)){
      rf_sd[l] <- rf_getsd(rf, newobs = learndata[l,], rfdata = learndata)
    }
  }
  
  # if(only_cf | (add_cf & (compare_sigma_area | compare_sigma_line)))
  if(only_cf | add_cf) {
    cf_mu <- predict(cf, type = "response")
    cf_sd <- cf_getsd(cf)
  }
  
  plotdata <- learndata[,c("y","x1","mu","sigma")]
  colnames(plotdata) <- c("y","x1","true.mu","true.sigma")
  sp <- plotdata[order(plotdata["x1"]),]

  if(nomodel){
    par(mar=c(5.1,4.1,4.1,3.1))
    plot(y = sp$y, x = sp$x, type = "p", col="grey", main = "True parameters", cex.main = 1.2, xaxt="n", yaxt="n", xlab = "", ylab = "")
    lines(x = sp$x, y = sp$true.mu, type = "l", col = 'black')
    polygon(c(sp$x, rev(sp$x)), c(sp$true.mu + sp$true.sigma, rev(sp$true.mu - sp$true.sigma)),
            col = transpgrey, border = "transparent")
    #lines(x = sp$x, y = sp$true.mu + sp$true.sigma, type = "l", lty = 1, col = 'grey')
    #lines(x = sp$x, y = sp$true.mu - sp$true.sigma, type = "l", lty = 1, col = 'grey')
    #legend('topleft', c(TeX('$\\mu$'), TeX('$\\mu \\pm \\sigma$')), 
    #       col = c('black', transpgrey), lty = 1, bty = "n", lwd = 2.5)
    legend(x = -0.5, y = 14, expression(mu %+-% sigma), bty = "n")
  }
  
  # disttree
  if(only_dt){
    
    plotdata <- cbind(learndata[,c("y","x1","mu","sigma")],
                      dt$fitted.par)
    colnames(plotdata) <- c("y","x1","true.mu","true.sigma", 
                            "fitted.mu.dt","fitted.sigma.dt")
    sp <- plotdata[order(plotdata["x1"]),]
    
    par(mar=c(5.1,4.1,4.1,3.1))
    plot(y = sp$y, x = sp$x, type = "p", col="grey", main = "disttree", col.main = pal["tree"], cex.main = 1.2, xaxt="n", yaxt="n", xlab = "", ylab = "")
    lines(x = sp$x, y = sp$true.mu, type = "l", col = 'black')
    lines(x = sp$x, y = sp$true.mu + sp$true.sigma, type = "l", col = 'grey')
    lines(x = sp$x, y = sp$true.mu - sp$true.sigma, type = "l", col = 'grey')
    lines(x = sp$x, y = sp$fitted.mu.dt, type = "l", col = pal["tree"], lwd = 2.5)
    polygon(c(sp$x, rev(sp$x)), c(sp$fitted.mu.dt + sp$fitted.sigma.dt, rev(sp$fitted.mu.dt - sp$fitted.sigma.dt)),
            col = pallight["tree"], border = "transparent")
    legend(x = -0.5, y = 14, expression(mu %+-% sigma), bty = "n")
  }
  
  # distforest
  if(only_df){
    
    plotdata <- cbind(learndata[,c("y","x1","mu","sigma")],
                      df$fitted.par)
    colnames(plotdata) <- c("y","x1","true.mu","true.sigma", 
                            "fitted.mu.df","fitted.sigma.df")
    sp <- plotdata[order(plotdata["x1"]),]
    par(mar=c(5.1,4.1,4.1,3.1))
    plot(y = sp$y, x = sp$x, type = "p", col="grey", main = "distforest", col.main = pal["forest"], cex.main = 1.2, xaxt="n", yaxt="n", xlab = "", ylab = "")
    lines(x = sp$x, y = sp$true.mu, type = "l", col = 'black')
    #polygon(c(sp$x, rev(sp$x)), c(sp$true.mu + sp$true.sigma, rev(sp$true.mu - sp$true.sigma)),
    #col = transpgrey, border = "transparent")
    lines(x = sp$x, y = sp$true.mu + sp$true.sigma, type = "l", col = 'grey')
    lines(x = sp$x, y = sp$true.mu - sp$true.sigma, type = "l", col = 'grey')
    lines(x = sp$x, y = sp$fitted.mu.df, type = "l", col = pal["forest"], lwd = 2.5)
    polygon(c(sp$x, rev(sp$x)), c(sp$fitted.mu.df + sp$fitted.sigma.df, rev(sp$fitted.mu.df - sp$fitted.sigma.df)),
            col = pallight["forest"], border = "transparent")
    #lines(x = sp$x, y = sp$fitted.mu.df + sp$fitted.sigma.df, type = "l", col = 'orange')
    #lines(x = sp$x, y = sp$fitted.mu.df - sp$fitted.sigma.df, type = "l", col = 'orange')
    #legend('topleft', c(TeX('$\\mu$'), TeX('$\\mu \\pm \\sigma$'), 
    #                    TeX('$\\hat{\\mu}$'), TeX('$\\hat{\\mu} \\pm \\hat{\\sigma}$')), 
    #       col = c('black','grey',pal["forest"], pallight["forest"]), lty = 1, bty = "n", lwd = 2.5)
    legend(x = -0.5, y = 14, expression(mu %+-% sigma), bty = "n")
  }
  
  # gamlss
  if(only_g){
    plotdata <- cbind(learndata[,c("y","x1","mu","sigma")],
                      g$mu.fv, g$sigma.fv)
    
    colnames(plotdata) <- c("y","x1","true.mu","true.sigma", 
                            "fitted.mu.g","fitted.sigma.g")
    sp <- plotdata[order(plotdata["x1"]),]
    
    par(mar=c(5.1,4.1,4.1,3.1))
    plot(y = sp$y, x = sp$x, type = "p", col="grey",  main = "gamlss", col.main = pal["gamlss"], cex.main = 1.2, xaxt="n", yaxt="n", xlab = "", ylab = "")
    lines(x = sp$x, y = sp$true.mu, type = "l", col = 'black')
    #polygon(c(sp$x, rev(sp$x)), c(sp$true.mu + sp$true.sigma, rev(sp$true.mu - sp$true.sigma)),
    #col = transpgrey, border = "transparent")
    lines(x = sp$x, y = sp$true.mu + sp$true.sigma, type = "l", col = 'grey')
    lines(x = sp$x, y = sp$true.mu - sp$true.sigma, type = "l", col = 'grey')
    lines(x = sp$x, y = sp$fitted.mu.g, type = "l", col = pal["gamlss"], lwd = 2.5)
    polygon(c(sp$x, rev(sp$x)), c(sp$fitted.mu.g + sp$fitted.sigma.g, rev(sp$fitted.mu.g - sp$fitted.sigma.g)),
            col = pallight["gamlss"], border = "transparent")
    #lines(x = sp$x, y = sp$fitted.mu.g + sp$fitted.sigma.g, type = "l", col = 'lightblue')
    #lines(x = sp$x, y = sp$fitted.mu.g - sp$fitted.sigma.g, type = "l", col = 'lightblue')
    #legend('topleft', c(TeX('$\\mu$'), TeX('$\\mu \\pm \\sigma$'), 
    #                    TeX('$\\hat{\\mu}$'), TeX('$\\hat{\\mu} \\pm \\hat{\\sigma}$')),  
    #       col = c('black','grey',pal["gamlss"], pallight["gamlss"]), lty = 1, bty = "n", lwd = 2.5)
    legend(x = -0.5, y = 14, expression(mu %+-% sigma), bty = "n")
  }
  
  # gamboostLSS
  if(only_gb){
    plotdata <- cbind(learndata[,c("y","x1","mu","sigma")],
                      fitted(gb)$mu, fitted(gb)$sigma)
    
    colnames(plotdata) <- c("y","x1","true.mu","true.sigma",
                            "fitted.mu.gb","fitted.sigma.gb")
    sp <- plotdata[order(plotdata["x1"]),]
    
    
    par(mar=c(5.1,4.1,4.1,3.1))
    plot(y = sp$y, x = sp$x, type = "p", col="grey",  main = "gamboostLSS", col.main = pal["gamboostLSS"], cex.main = 1.2, xaxt="n", yaxt="n", xlab = "", ylab = "")
    lines(x = sp$x, y = sp$true.mu, type = "l", col = 'black')
    #polygon(c(sp$x, rev(sp$x)), c(sp$true.mu + sp$true.sigma, rev(sp$true.mu - sp$true.sigma)),
    #col = transpgrey, border = "transparent")
    lines(x = sp$x, y = sp$true.mu + sp$true.sigma, type = "l", col = 'grey')
    lines(x = sp$x, y = sp$true.mu - sp$true.sigma, type = "l", col = 'grey')
    lines(x = sp$x, y = sp$fitted.mu.gb, type = "l", col = pal["gamboostLSS"], lwd = 2.5)
    polygon(c(sp$x, rev(sp$x)), c(sp$fitted.mu.gb + sp$fitted.sigma.gb, rev(sp$fitted.mu.gb - sp$fitted.sigma.gb)),
            col = pallight["gamboostLSS"], border = "transparent")
    legend(x = -0.5, y = 14, expression(mu %+-% sigma), bty = "n")
  }
  
  # bamlss
  if(only_b){
    plotdata <- cbind(learndata[,c("y","x1","mu","sigma")],
                      fitted(b)$mu[,"Mean"], fitted(b)$sigma[,"Mean"])
    
    colnames(plotdata) <- c("y","x1","true.mu","true.sigma",
                            "fitted.mu.b","fitted.sigma.b")
    sp <- plotdata[order(plotdata["x1"]),]
    
    
    par(mar=c(5.1,4.1,4.1,3.1))
    plot(y = sp$y, x = sp$x, type = "p", col="grey",  main = "bamlss", col.main = pal["bamlss"], cex.main = 1.2, xaxt="n", yaxt="n", xlab = "", ylab = "")
    lines(x = sp$x, y = sp$true.mu, type = "l", col = 'black')
    #polygon(c(sp$x, rev(sp$x)), c(sp$true.mu + sp$true.sigma, rev(sp$true.mu - sp$true.sigma)),
    #col = transpgrey, border = "transparent")
    lines(x = sp$x, y = sp$true.mu + sp$true.sigma, type = "l", col = 'grey')
    lines(x = sp$x, y = sp$true.mu - sp$true.sigma, type = "l", col = 'grey')
    lines(x = sp$x, y = sp$fitted.mu.b, type = "l", col = pal["bamlss"], lwd = 2.5)
    polygon(c(sp$x, rev(sp$x)), c(sp$fitted.mu.b + sp$fitted.sigma.b, rev(sp$fitted.mu.b - sp$fitted.sigma.b)),
            col = pallight["bamlss"], border = "transparent")
    legend(x = -0.5, y = 14, expression(mu %+-% sigma), bty = "n")
  }
  
  # randomForest 
  if(only_rf){
    plotdata <- cbind(learndata[,c("y","x1","mu","sigma")],
                      rf$predicted, rf_sd)
    
    colnames(plotdata) <- c("y","x1","true.mu","true.sigma",
                            "fitted.mu.rf","fitted.sigma.rf")
    sp <- plotdata[order(plotdata["x1"]),]
    
    
    par(mar=c(5.1,4.1,4.1,3.1))
    plot(y = sp$y, x = sp$x, type = "p", col="grey",  main = "randomForest", col.main = pal["randomForest"], cex.main = 1.2, xaxt="n", yaxt="n", xlab = "", ylab = "")
    lines(x = sp$x, y = sp$true.mu, type = "l", col = 'black')
    #polygon(c(sp$x, rev(sp$x)), c(sp$true.mu + sp$true.sigma, rev(sp$true.mu - sp$true.sigma)),
    #col = transpgrey, border = "transparent")
    lines(x = sp$x, y = sp$true.mu + sp$true.sigma, type = "l", col = 'grey')
    lines(x = sp$x, y = sp$true.mu - sp$true.sigma, type = "l", col = 'grey')
    lines(x = sp$x, y = sp$fitted.mu.rf, type = "l", col = pal["randomForest"], lwd = 2.5)
    polygon(c(sp$x, rev(sp$x)), c(sp$fitted.mu.rf + sp$fitted.sigma.rf, rev(sp$fitted.mu.rf - sp$fitted.sigma.rf)),
            col = pallight["randomForest"], border = "transparent")
    legend(x = -0.5, y = 14, expression(mu %+-% sigma), bty = "n")
  }
  
  # cforest
  if(only_cf){
    plotdata <- cbind(learndata[,c("y","x1","mu","sigma")],
                      cf_mu, cf_sd)
    
    colnames(plotdata) <- c("y","x1","true.mu","true.sigma",
                            "fitted.mu.cf","fitted.sigma.cf")
    sp <- plotdata[order(plotdata["x1"]),]
    
    
    par(mar=c(5.1,4.1,4.1,3.1))
    plot(y = sp$y, x = sp$x, type = "p", col="grey",  main = "cforest", col.main = pal["cforest"], cex.main = 1.2, xaxt="n", yaxt="n", xlab = "", ylab = "")
    lines(x = sp$x, y = sp$true.mu, type = "l", col = 'black')
    #polygon(c(sp$x, rev(sp$x)), c(sp$true.mu + sp$true.sigma, rev(sp$true.mu - sp$true.sigma)),
    #col = transpgrey, border = "transparent")
    lines(x = sp$x, y = sp$true.mu + sp$true.sigma, type = "l", col = 'grey')
    lines(x = sp$x, y = sp$true.mu - sp$true.sigma, type = "l", col = 'grey')
    lines(x = sp$x, y = sp$fitted.mu.cf, type = "l", col = pal["cforest"], lwd = 2.5)
    polygon(c(sp$x, rev(sp$x)), c(sp$fitted.mu.cf + sp$fitted.sigma.cf, rev(sp$fitted.mu.cf - sp$fitted.sigma.cf)),
            col = pallight["cforest"], border = "transparent")
    legend(x = -0.5, y = 14, expression(mu %+-% sigma), bty = "n")
  }
  
  
  
  # compare location parameter mu
  if(compare_mu){
    plotdata <- learndata[,c("y","x1","mu","sigma")]
    coln <- c("y","x1","true.mu","true.sigma")
    
    if(add_dt) {
      plotdata <- cbind(plotdata, dt$fitted.par)
      coln <- c(coln, "fitted.mu.dt","fitted.sigma.dt")
    }
    if(add_df) {
      plotdata <- cbind(plotdata, df$fitted.par)
      coln <- c(coln, "fitted.mu.df","fitted.sigma.df")
    }
    if(add_g) {
      plotdata <- cbind(plotdata, g$mu.fv, g$sigma.fv)
      coln <- c(coln, "fitted.mu.g","fitted.sigma.g")
    }
    if(add_gb) {
      plotdata <- cbind(plotdata, fitted(gb)$mu, fitted(gb)$sigma)
      coln <- c(coln, "fitted.mu.gb","fitted.sigma.gb")
    }
    if(add_b) {
      plotdata <- cbind(plotdata, fitted(b)$mu[,"Mean"], fitted(b)$sigma[,"Mean"])
      coln <- c(coln, "fitted.mu.b","fitted.sigma.b")
    }
    if(add_rf) {
      plotdata <- cbind(plotdata, rf$predicted, rf_sd)
      coln <- c(coln, "fitted.mu.rf","fitted.sigma.rf")
    }
    if(add_cf) {
      plotdata <- cbind(plotdata, cf_mu, cf_sd)
      coln <- c(coln, "fitted.mu.cf","fitted.sigma.cf")
    }
    
    colnames(plotdata) <- coln  
    sp <- plotdata[order(plotdata["x1"]),]
    
    
    par(mar=c(5.1,4.1,4.1,3.1))
    plot(y = sp$y, x = sp$x, type = "p", col="grey", xaxt="n", yaxt="n", xlab = "", ylab = "")
    #title(main = "disttree", col.main = pal["tree"], cex.main = 1.2, font.main = 2, adj = 0.19)
    #title(main = "distforest", col.main = pal["forest"], cex.main = 1.2, font.main = 2, adj = 0.5)
    #title(main = "gamlss", col.main = pal["gamlss"], cex.main = 1.2, font.main = 2, adj = 0.8)
    #title(main = "vs.", col.main = "black", cex.main = 1.2, font.main = 2, adj = 0.35)
    #title(main = "vs.", col.main = "black", cex.main = 1.2, font.main = 2, adj = 0.65)
    
    lines(x = sp$x, y = sp$true.mu, type = "l", col = 'black')
    if(add_dt) lines(x = sp$x, y = sp$fitted.mu.dt, type = "l", col = pal["tree"], lwd = 2)
    if(add_df) lines(x = sp$x, y = sp$fitted.mu.df, type = "l", col = pal["forest"], lwd = 2)
    if(add_g) lines(x = sp$x, y = sp$fitted.mu.g, type = "l", col = pal["gamlss"], lwd = 2)
    if(add_gb) lines(x = sp$x, y = sp$fitted.mu.gb, type = "l", col = pal["gamboostLSS"], lwd = 2)
    if(add_b) lines(x = sp$x, y = sp$fitted.mu.b, type = "l", col = pal["bamlss"], lwd = 2)
    if(add_rf) lines(x = sp$x, y = sp$fitted.mu.rf, type = "l", col = pal["randomForest"], lwd = 2)
    if(add_cf) lines(x = sp$x, y = sp$fitted.mu.cf, type = "l", col = pal["cforest"], lwd = 2)
    #legend('topleft', c(TeX('$\\mu$'), TeX('$\\hat{\\mu}$ disttree'),
    #                    TeX('$\\hat{\\mu}$ distforest'), TeX('$\\hat{\\mu}$ gamlss')), 
    #       col = c('black',pal["tree"],pal["forest"],pal["gamlss"]), lty = 1, bty = "n", lwd = 2.5)
    legend(x = -0.5, y = 13.95, expression(mu), bty = "n")
  }
  
  # variance
  if(compare_sigma_area){
    coln <- c("y","x1","true.mu","true.sigma")
    
    if(add_dt) {
      plotdata <- cbind(plotdata, dt$fitted.par)
      coln <- c(coln, "fitted.mu.dt","fitted.sigma.dt")
    }
    if(add_df) {
      plotdata <- cbind(plotdata, df$fitted.par)
      coln <- c(coln, "fitted.mu.df","fitted.sigma.df")
    }
    if(add_g) {
      plotdata <- cbind(plotdata, g$mu.fv, g$sigma.fv)
      coln <- c(coln, "fitted.mu.g","fitted.sigma.g")
    }
    if(add_gb) {
      plotdata <- cbind(plotdata, fitted(gb)$mu, fitted(gb)$sigma)
      coln <- c(coln, "fitted.mu.gb","fitted.sigma.gb")
    }
    if(add_b) {
      plotdata <- cbind(plotdata, fitted(b)$mu[,"Mean"], fitted(b)$sigma[,"Mean"])
      coln <- c(coln, "fitted.mu.b","fitted.sigma.b")
    }
    if(add_rf) {
      plotdata <- cbind(plotdata, rf$predicted, rf_sd)
      coln <- c(coln, "fitted.mu.rf","fitted.sigma.rf")
    }
    if(add_cf) {
      plotdata <- cbind(plotdata, cf_mu, cf_sd)
      coln <- c(coln, "fitted.mu.cf","fitted.sigma.cf")
    }
    
    colnames(plotdata) <- coln  
    sp <- plotdata[order(plotdata["x1"]),]
    
    par(mar=c(5.1,4.1,4.1,3.1))
    plot(y = sp$y, x = sp$x, type = "p", col="grey", xaxt="n", yaxt="n", xlab = "", ylab = "")
    #title(main = "disttree", col.main = pal["tree"], cex.main = 1.2, font.main = 2, adj = 0.19)
    #title(main = "distforest", col.main = pal["forest"], cex.main = 1.2, font.main = 2, adj = 0.5)
    #title(main = "gamlss", col.main = pal["gamlss"], cex.main = 1.2, font.main = 2, adj = 0.8)
    #title(main = "vs.", col.main = "black", cex.main = 1.2, font.main = 2, adj = 0.35)
    #title(main = "vs.", col.main = "black", cex.main = 1.2, font.main = 2, adj = 0.65)
    
    lines(x = sp$x, y = sp$true.mu + sp$true.sigma, type = "l", col = 'black')
    lines(x = sp$x, y = sp$true.mu - sp$true.sigma, type = "l", col = 'black')
    
    if(add_dt) polygon(c(sp$x, rev(sp$x)), c(sp$fitted.mu.dt + sp$fitted.sigma.dt, rev(sp$fitted.mu.dt - sp$fitted.sigma.dt)),
            col = pallight["tree"], border = "transparent")
    if(add_df) polygon(c(sp$x, rev(sp$x)), c(sp$fitted.mu.df + sp$fitted.sigma.df, rev(sp$fitted.mu.df - sp$fitted.sigma.df)),
            col = pallight["forest"], border = "transparent")
    if(add_g) polygon(c(sp$x, rev(sp$x)), c(sp$fitted.mu.g + sp$fitted.sigma.g, rev(sp$fitted.mu.g - sp$fitted.sigma.g)),
            col = pallight["gamlss"], border = "transparent")
    if(add_gb) polygon(c(sp$x, rev(sp$x)), c(sp$fitted.mu.gb + sp$fitted.sigma.gb, rev(sp$fitted.mu.gb - sp$fitted.sigma.gb)),
                      col = pallight["gambootLSS"], border = "transparent")
    if(add_b) polygon(c(sp$x, rev(sp$x)), c(sp$fitted.mu.b + sp$fitted.sigma.b, rev(sp$fitted.mu.b - sp$fitted.sigma.b)),
                      col = pallight["bamlss"], border = "transparent")
    if(add_rf) polygon(c(sp$x, rev(sp$x)), c(sp$fitted.mu.rf + sp$fitted.sigma.rf, rev(sp$fitted.mu.rf - sp$fitted.sigma.rf)),
                      col = pallight["randomForest"], border = "transparent")
    if(add_cf) polygon(c(sp$x, rev(sp$x)), c(sp$fitted.mu.cf + sp$fitted.sigma.cf, rev(sp$fitted.mu.cf - sp$fitted.sigma.cf)),
                      col = pallight["cforest"], border = "transparent")
    legend(x = -0.5, y = 14, expression(mu %+-% sigma), bty = "n")
  }

  
  # variance
  if(compare_sigma_line){
    
    plotdata <- learndata[,c("y","x1","mu","sigma")]
    coln <- c("y","x1","true.mu","true.sigma")
    
    if(add_dt) {
      plotdata <- cbind(plotdata, dt$fitted.par)
      coln <- c(coln, "fitted.mu.dt","fitted.sigma.dt")
    }
    if(add_df) {
      plotdata <- cbind(plotdata, df$fitted.par)
      coln <- c(coln, "fitted.mu.df","fitted.sigma.df")
    }
    if(add_g) {
      plotdata <- cbind(plotdata, g$mu.fv, g$sigma.fv)
      coln <- c(coln, "fitted.mu.g","fitted.sigma.g")
    }
    if(add_gb) {
      plotdata <- cbind(plotdata, fitted(gb)$mu, fitted(gb)$sigma)
      coln <- c(coln, "fitted.mu.gb","fitted.sigma.gb")
    }
    if(add_b) {
      plotdata <- cbind(plotdata, fitted(b)$mu[,"Mean"], fitted(b)$sigma[,"Mean"])
      coln <- c(coln, "fitted.mu.b","fitted.sigma.b")
    }
    if(add_rf) {
      plotdata <- cbind(plotdata, rf$predicted, rf_sd)
      coln <- c(coln, "fitted.mu.rf","fitted.sigma.rf")
    }
    if(add_cf) {
      plotdata <- cbind(plotdata, cf_mu, cf_sd)
      coln <- c(coln, "fitted.mu.cf","fitted.sigma.cf")
    }
    
    colnames(plotdata) <- coln  
    sp <- plotdata[order(plotdata["x1"]),]
    
    par(mar=c(5.1,4.1,4.1,3.1))
    ymin <- min(sp[,c(2:(length(coln)/2))*2])
    ymax <- max(sp[,c(2:(length(coln)/2))*2])
    #title(main = "disttree", col.main = pal["tree"], cex.main = 1.2, font.main = 2, adj = 0.19)
    #title(main = "distforest", col.main = pal["forest"], cex.main = 1.2, font.main = 2, adj = 0.5)
    #title(main = "gamlss", col.main = pal["gamlss"], cex.main = 1.2, font.main = 2, adj = 0.8)
    #title(main = "vs.", col.main = "black", cex.main = 1.2, font.main = 2, adj = 0.35)
    #title(main = "vs.", col.main = "black", cex.main = 1.2, font.main = 2, adj = 0.65)
    
    plot(x = sp$x, y = sp$true.sigma, type = "l", col = 'black', 
         xaxt="n", yaxt="n", xlab = "", ylab = "", ylim = c(ymin, ymax))
    
    if(add_dt) lines(x = sp$x, y = sp$fitted.sigma.dt, type = "l", col = pal["tree"], lwd = 2)
    if(add_df) lines(x = sp$x, y = sp$fitted.sigma.df, type = "l", col = pal["forest"], lwd = 2)
    if(add_g) lines(x = sp$x, y = sp$fitted.sigma.g, type = "l", col = pal["gamlss"], lwd = 2)
    if(add_gb) lines(x = sp$x, y = sp$fitted.sigma.gb, type = "l", col = pal["gamboostLSS"], lwd = 2)
    if(add_b) lines(x = sp$x, y = sp$fitted.sigma.b, type = "l", col = pal["bamlss"], lwd = 2)
    if(add_rf) lines(x = sp$x, y = sp$fitted.sigma.rf, type = "l", col = pal["randomForest"], lwd = 2)
    if(add_cf) lines(x = sp$x, y = sp$fitted.sigma.cf, type = "l", col = pal["cforest"], lwd = 2)
    
    legend(x = -0.5, y = 14, expression(sigma), bty = "n")
  }
  
}



###############################################
# test
if(FALSE){
l <- sim_onecov_test

plot_onecov(learndata = l$learndata, dt = l$dt, df = l$df, g = l$g, b = l$b, rf = l$rf, cf = l$cf, 
            compare_mu = TRUE, 
            add_dt = TRUE, add_df = TRUE, add_g = TRUE, add_b = TRUE, add_rf = FALSE, add_cf = FALSE)

plot_onecov(learndata = l$learndata, dt = l$dt, df = l$df, g = l$g, b = l$b, rf = l$rf, cf = l$cf, 
            compare_sigma_area = TRUE, 
            add_dt = TRUE, add_df = TRUE, add_g = TRUE, add_b = TRUE, add_rf = TRUE, add_cf = TRUE)

plot_onecov(learndata = l$learndata, dt = l$dt, g = l$g, b = l$b, 
            only_b = TRUE)
}



#################################################
if(FALSE){
plot1(learndata, dt, df, g, gb, b, rf, cf, compare_mu = TRUE, add_dt = TRUE,
    add_df = TRUE, add_g = TRUE, add_gb = TRUE, add_b = TRUE)#, add_rf = TRUE, add_cf = TRUE)


plot1(learndata, dt, df, g, gb, b, rf, cf, compare_sigma_area = TRUE, add_dt = TRUE,
      add_df = TRUE, add_g = TRUE, add_gb = TRUE, add_b = TRUE)

plot1(learndata, dt, df, g, gb, b, rf, cf, compare_sigma_line = TRUE, add_dt = TRUE,
      add_df = TRUE, add_g = TRUE, add_gb = TRUE, add_b = FALSE)


plot1(learndata, dt, df, g, gb, b, rf, cf, only_dt = TRUE)
plot1(learndata, dt, df, g, gb, b, rf, cf, only_df = TRUE)
plot1(learndata, dt, df, g, gb, b, rf, cf, only_g = TRUE)
plot1(learndata, dt, df, g, gb, b, rf, cf, only_gb = TRUE)
plot1(learndata, dt, df, g, gb, b, rf, cf, only_b = TRUE)
plot1(learndata, dt, df, g, gb, b, rf, cf, only_rf = TRUE)
plot1(learndata, dt, df, g, gb, b, rf, cf, only_cf = TRUE)
}