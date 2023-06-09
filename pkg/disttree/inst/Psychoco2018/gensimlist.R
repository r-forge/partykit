### generates simlist which will then be plotted in slides.Rnw
gensimlist <- function(seed = 7, nrep = 100, ntree = 100, covariate = c("x1"))
{
  library("disttree")
  library("gamlss")
  library("parallel")
  
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
  
  
  
  
  ## data generating process
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
  
  
  ### get logLik from disttree object for new data 
  dtll.newdata <- function(object, newdata){
    ll <- 0
    nd <- newdata[,-c(1,ncol(newdata)-1, ncol(newdata))]
    if(is.null(coef(object))){
      readline(prompt="Press [enter] to continue")
      print(object$info$call)
      print(object$data)
      return(NULL)
    } else {
      pred.par <- predict(object, newdata = nd, type = "parameter")
      for(i in 1:(nrow(newdata))){
        eta <- as.numeric(object$info$family$linkfun(pred.par[i,]))
        ll <- ll + object$info$family$ddist(newdata[i,paste(object$info$formula[[2]])], eta = eta, log=TRUE)
      }
      if(is.na(ll)) print("dt.ll = NA")
      return(ll)
    }
  }
  
  
  
  ### get logLik from gamlss object for new data and normal distribution
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
    distfun <- if(object$family[1] == "NOlc") {
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
  
  
  
  ### simulation wrapper function
  # (optionally with plot)
  sim_fun <- function(fun, nsteps = 10, ntree = 10, nobs = 200, nrep = 10, 
                      family = dist_list_normal, kappa.start = 0, stepsize = 1,
                      type.gfun = "linear", covariate = covariate,
                      g.interact = FALSE, var.nobs = FALSE, kappa.fix = NULL,
                      plot = c("none", "RMSE", "loglikelihood"), type.tree = "mob",
                      ctree_minbucket = 7, ctree_mincrit = 0, 
                      cforest_minbucket = 7, cforest_mincrit = 0)
  {
    
    if(var.nobs && (is.null(kappa.fix))) stop("for varying number of observations a value for kappa.fix has to be set")
    
    
    if(var.nobs) {
      nobs.start <- 100
      f <- function(x) fun(x, kappa.fix)
    }
    
    
    #### parallelization
    
    res <- mclapply(1:(nsteps+1),
                    function(i){
                      
                      # for each step the parameter function is defined
                      if(var.nobs) {
                        nobs <- nobs.start + (i-1)*stepsize
                      } else {
                        f <- function(x) {fun(x, kappa.start + (i-1)*stepsize)} 
                      }
                      
                      rmse.exp.true.dt <- numeric(length = nrep)
                      rmse.exp.true.df <- numeric(length = nrep)
                      rmse.exp.true.g <- numeric(length = nrep)
                      
                      rmse.exp.obs.dt <- numeric(length = nrep)
                      rmse.exp.obs.df <- numeric(length = nrep)
                      rmse.exp.obs.g <- numeric(length = nrep)
                      
                      rmse.mu.dt <- numeric(length = nrep)
                      rmse.mu.df <- numeric(length = nrep)
                      rmse.mu.g <- numeric(length = nrep)
                      
                      rmse.sigma.dt <- numeric(length = nrep)
                      rmse.sigma.df <- numeric(length = nrep)
                      rmse.sigma.g <- numeric(length = nrep)
                      
                      loglik.dt <- numeric(length = nrep)
                      loglik.df <- numeric(length = nrep)
                      loglik.g <- numeric(length = nrep)
                      
                      #data <- list()
                      
                      for(j in 1:nrep){
                        learndata <- dgp(nobs, family = family, round.sp = 7, fun = f)
                        newdata <- dgp(500, family = family, round.sp = 7, fun = f)
                        nd <- newdata[,-c(1,ncol(newdata)-1,ncol(newdata))]
                        #data[[j]] <- list(learndata = learndata, newdata = newdata, nd = nd)
                        
                        if(("x1" %in% covariate) && !("x2" %in%covariate)){
                          if(type.tree == "mob"){
                            dt <- disttree(y~x1, data=learndata, family=family, 
                                           type.tree = "mob", control = mob_control(restart = FALSE))
                            df <- distforest(y~x1, data=learndata, family=family, ntree=ntree, type.tree = "mob", control = mob_control(restart = FALSE))
                          }
                          if(type.tree == "ctree"){
                            control <- ctree_control(teststat = "quad", testtype = "Bonferroni", 
                                                     mincriterion = ctree_mincrit, minbucket = ctree_minbucket)
                            dt <- disttree(y~x1, data=learndata, family=family, 
                                           type.tree = "ctree", control = control)
                            control <- ctree_control(teststat = "quad", testtype = "Univ", 
                                                     mincriterion = cforest_mincrit, minbucket = cforest_minbucket)
                            df <- distforest(y~x1, data=learndata, family=family, ntree=ntree, type.tree = type.tree, control = control)
                          }
                          
                          if(family$family.name == "censored Normal Distribution") {
                            if(type.gfun =="linear") g <- gamlss(Surv(y, y>0, type="left") ~ x1, sigma.formula = ~ x1, data = learndata, trace = FALSE, family = NOlc())
                            if(type.gfun =="smooth") g <- gamlss(Surv(y, y>0, type="left") ~ cs(x1), sigma.formula = ~ cs(x1), data = learndata, trace = FALSE, family = NOlc())
                            if(type.gfun =="nonparametric")  g <- gamlss(Surv(y, y>0, type="left") ~ pb(x1), sigma.formula = ~ pb(x1), data = learndata, trace = FALSE, family = NOlc())
                          } else {
                            if(type.gfun =="linear") g <- gamlss(y ~ x1, sigma.formula = ~ x1, data = learndata, trace = FALSE)
                            if(type.gfun =="smooth") g <- gamlss(y ~ cs(x1), sigma.formula = ~ cs(x1), data = learndata, trace = FALSE)
                            if(type.gfun =="nonparametric")  g <- gamlss(y ~ pb(x1), sigma.formula = ~ pb(x1), data = learndata, trace = FALSE)
                          }
                        }
                        if(("x1" %in% covariate) && ("x2" %in% covariate)){
                          if(type.tree == "mob"){
                            dt <- disttree(y~x1+x2, data=learndata, family=family, 
                                           type.tree = "mob", control = mob_control(restart = FALSE))
                            df <- distforest(y~x1+x2, data=learndata, family=family, ntree=ntree, type.tree = "mob", control = mob_control(restart = FALSE))
                          }
                          if(type.tree == "ctree"){
                            control <- ctree_control(teststat = "quad", testtype = "Bonferroni", 
                                                     mincriterion = ctree_mincrit, minbucket = ctree_minbucket)
                            dt <- disttree(y~x1+x2, data=learndata, family=family, 
                                           type.tree = type.tree, control = control)
                            control <- ctree_control(teststat = "quad", testtype = "Univ", 
                                                     mincriterion = cforest_mincrit, minbucket = cforest_minbucket)
                            df <- distforest(y~x1+x2, data=learndata, family=family, ntree=ntree, type.tree = type.tree, control = control)
                          }
                          
                          #df <- distforest(y~x1+x2, data=learndata, family=family, ntree=ntree, fitted.OOB = FALSE)
                          #df <- distforest(y~x1+x2, data=learndata, family=family, ntree=ntree, control = ctree_control(teststat = "quad", testtype = "Univ", mincriterion = 0, minsplit = 10))
                          
                          if(family$family.name == "censored Normal Distribution") {
                            if(type.gfun =="linear") {
                              if(g.interact) {
                                g <- gamlss(Surv(y, y>0, type="left") ~ x1*x2, sigma.formula =~x1*x2, data = learndata, trace = FALSE, family = NOlc())
                              } else {
                                g <-  gamlss(Surv(y, y>0, type="left") ~ x1+x2, sigma.formula =~x1+x2, data = learndata, trace = FALSE, family = NOlc())
                              } 
                            }
                            if(type.gfun =="smooth") {
                              if(g.interact) {
                                g <- gamlss(Surv(y, y>0, type="left") ~ cs(x1)+cs(x2)+x1:x2, sigma.formula =~cs(x1)+cs(x2)+x1:x2, data = learndata, trace = FALSE, family = NOlc())
                              } else {
                                g <- gamlss(Surv(y, y>0, type="left") ~ cs(x1)+cs(x2), sigma.formula =~cs(x1)+cs(x2), data = learndata, trace = FALSE, family = NOlc())
                              }
                            } 
                            if(type.gfun =="nonparametric") {
                              if(g.interact) {
                                g <- gamlss(Surv(y, y>0, type="left") ~ pb(x1)+pb(x2)+x1:x2, sigma.formula =~pb(x1)+pb(x2)+x1:x2, data = learndata, trace = FALSE, family = NOlc())
                              } else {
                                g <- try(gamlss(Surv(y, y>0, type="left") ~ pb(x1)+pb(x2), sigma.formula =~pb(x1)+pb(x2), data = learndata, trace = FALSE, family = NOlc()))
                                if(inherits(g, "try-error")) {
                                  print("try-error in gamlss")
                                  g <- NULL
                                }
                              }
                            }
                            
                          } else {
                            
                            {
                              if(type.gfun =="linear") {
                                if(g.interact) {
                                  g <- gamlss(y ~ x1*x2, sigma.formula =~x1*x2, data = learndata, trace = FALSE)
                                } else {
                                  g <-  gamlss(y ~ x1+x2, sigma.formula =~x1+x2, data = learndata, trace = FALSE)
                                } 
                              }
                              if(type.gfun =="smooth") {
                                if(g.interact) {
                                  g <- gamlss(y ~ cs(x1)+cs(x2)+x1:x2, sigma.formula =~cs(x1)+cs(x2)+x1:x2, data = learndata, trace = FALSE)
                                } else {
                                  g <- gamlss(y ~ cs(x1)+cs(x2), sigma.formula =~cs(x1)+cs(x2), data = learndata, trace = FALSE)
                                }
                              } 
                              if(type.gfun =="nonparametric") {
                                if(g.interact) {
                                  g <- gamlss(y ~ pb(x1)+pb(x2)+x1:x2, sigma.formula =~pb(x1)+pb(x2)+x1:x2, data = learndata, trace = FALSE)
                                } else {
                                  g <- try(gamlss(y ~ pb(x1)+pb(x2), sigma.formula =~pb(x1)+pb(x2), data = learndata, trace = FALSE))
                                  if(inherits(g, "try-error")) {
                                    print("try-error in gamlss")
                                    g <- NULL
                                  }
                                }
                              }
                            }
                          }
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
                        
                        df_coef <- predict(df, newdata = nd, type = "parameter") # returns fitted values and fitted parameters
                        df_mu <- df_coef$mu
                        df_sigma <- df_coef$sigma
                        
                        g_mu <- predict(g, newdata = nd, what = "mu", type = "response", data = learndata)
                        g_sigma <- predict(g, newdata = nd, what = "sigma", type = "response", data = learndata)
                        
                        if("censored" %in% strsplit(family$family.name, " ")[[1]]){
                          
                          #calculate expected value for censored data
                          true_exp <- pnorm(true_mu/true_sigma) * (true_mu + true_sigma * (dnorm(true_mu/true_sigma) / pnorm(true_mu/true_sigma)))
                          dt_exp <- pnorm(dt_mu/dt_sigma) * (dt_mu + dt_sigma * (dnorm(dt_mu/dt_sigma) / pnorm(dt_mu/dt_sigma)))
                          df_exp <- pnorm(df_mu/df_sigma) * (df_mu + df_sigma * (dnorm(df_mu/df_sigma) / pnorm(df_mu/df_sigma)))
                          g_exp <- pnorm(g_mu/g_sigma) * (g_mu + g_sigma * (dnorm(g_mu/g_sigma) / pnorm(g_mu/g_sigma)))
                          
                        } else {
                          
                          # for non censored normal data the expected value is the location paramter
                          true_exp <- true_mu
                          dt_exp <- dt_mu
                          df_exp <- df_mu
                          g_exp <- g_mu
                        }
                        
                        
                        # calculate RMSEs
                        rmse.exp.true.dt[j] <- sqrt(mean((dt_exp - true_exp)^2))
                        rmse.exp.true.df[j] <- sqrt(mean((df_exp - true_exp)^2))
                        rmse.exp.true.g[j] <- sqrt(mean((g_exp - true_exp)^2))
                        
                        rmse.exp.obs.dt[j] <- sqrt(mean((dt_exp - newdata[,"y"])^2))
                        rmse.exp.obs.df[j] <- sqrt(mean((df_exp - newdata[,"y"])^2))
                        rmse.exp.obs.g[j] <- sqrt(mean((g_exp - newdata[,"y"])^2))
                        
                        rmse.mu.dt[j] <- sqrt(mean((dt_mu - true_mu)^2))
                        rmse.mu.df[j] <- sqrt(mean((df_mu - true_mu)^2))
                        rmse.mu.g[j] <- sqrt(mean((g_mu  - true_mu)^2))
                        
                        rmse.sigma.dt[j] <- sqrt(mean((dt_sigma - true_sigma)^2))
                        rmse.sigma.df[j] <- sqrt(mean((df_sigma - true_sigma)^2))
                        rmse.sigma.g[j] <- sqrt(mean((g_sigma - true_sigma)^2))
                        
                        # calculate loglikelihood for newdata
                        loglik.dt[j] <- dtll.newdata(dt, newdata = newdata)
                        loglik.df[j] <- as.numeric(logLik(df, newdata = newdata))
                        loglik.g[j] <- gll.newdata(g, newdata = newdata, data = learndata)
                      }
                      
                      
                      return(as.data.frame(c(mean(rmse.exp.true.dt),
                                             mean(rmse.exp.true.df),
                                             mean(rmse.exp.true.g, na.rm = TRUE),
                                             mean(rmse.exp.obs.dt),
                                             mean(rmse.exp.obs.df),
                                             mean(rmse.exp.obs.g, na.rm = TRUE),
                                             mean(rmse.mu.dt),
                                             mean(rmse.mu.df),
                                             mean(rmse.mu.g, na.rm = TRUE),
                                             mean(rmse.sigma.dt),
                                             mean(rmse.sigma.df),
                                             mean(rmse.sigma.g, na.rm = TRUE),
                                             mean(loglik.dt),
                                             mean(loglik.df),
                                             mean(loglik.g, na.rm = TRUE)), 
                                           row.names = c( "av.rmse.exp.true.dt",
                                                          "av.rmse.exp.true.df",
                                                          "av.rmse.exp.true.g",
                                                          "av.rmse.exp.obs.dt",
                                                          "av.rmse.exp.obs.df",
                                                          "av.rmse.exp.obs.g",
                                                          "av.rmse.mu.dt",
                                                          "av.rmse.mu.df",
                                                          "av.rmse.mu.g",
                                                          "av.rmse.sigma.dt",
                                                          "av.rmse.sigma.df",
                                                          "av.rmse.sigma.g",
                                                          "av.loglik.dt",
                                                          "av.loglik.df",
                                                          "av.loglik.g")))
                      
                    },
                    mc.cores = detectCores() - 1
    )
    #stopCluster(cl)
    
    # res is a list with nsteps+1 elements, each element is a data frame with the average values of RMSE and loglikelihood
    # transform to one data frame
    resframe <- res[[1]]
    for(k in 2: length(res)){
      resframe <- cbind(resframe, res[[k]])
    }  
    colnames(resframe) <- paste(c(1:length(res)))
    
    x.axis <- if(var.nobs) {
      seq(100, (100 + (nsteps)*stepsize), stepsize) 
    } else {
      seq(kappa.start, (kappa.start + (nsteps)*stepsize), stepsize)
    }
    
    
    simlist <- list(
      av.rmse.exp.true.dt = as.numeric(resframe["av.rmse.exp.true.dt",]),
      av.rmse.exp.true.df = as.numeric(resframe["av.rmse.exp.true.df",]),
      av.rmse.exp.true.g = as.numeric(resframe["av.rmse.exp.true.g",]),
      
      av.rmse.exp.obs.dt = as.numeric(resframe["av.rmse.exp.obs.dt",]),
      av.rmse.exp.obs.df = as.numeric(resframe["av.rmse.exp.obs.df",]),
      av.rmse.exp.obs.g = as.numeric(resframe["av.rmse.exp.obs.g",]),
      
      av.rmse.mu.dt = as.numeric(resframe["av.rmse.mu.dt",]),
      av.rmse.mu.df = as.numeric(resframe["av.rmse.mu.df",]),
      av.rmse.mu.g = as.numeric(resframe["av.rmse.mu.g",]),
      
      av.rmse.sigma.dt = as.numeric(resframe["av.rmse.sigma.dt",]),
      av.rmse.sigma.df = as.numeric(resframe["av.rmse.sigma.df",]),
      av.rmse.sigma.g = as.numeric(resframe["av.rmse.sigma.g",]),
      
      av.loglik.dt = as.numeric(resframe["av.loglik.dt",]),
      av.loglik.df = as.numeric(resframe["av.loglik.df",]),
      av.loglik.g
      = as.numeric(resframe["av.loglik.g",]),
      
      x.axis = x.axis
    )
    
    
    return(simlist)
  }
  
  
  
  
  
  # parameter function
  #fun <- function(x, kappa){cbind(10 * exp(-((4*x[,1] - 2)^(2*kappa))), rep.int(2, nrow(x)))}
  fun <- function(x, kappa){cbind(10 * exp(-((4*x[,1] - 2)^(2*kappa))), 0.5 + 2*abs(x[,1]))}
  
  set.seed(seed)
  
  
  simlist <- sim_fun(fun = fun, kappa.start = 1, nsteps = 10, stepsize = 7, 
                     ntree = ntree, nobs = 1000, nrep = nrep, type.gfun = "nonparametric",
                     covariate = c("x1"), plot = "none", type.tree = "ctree",
                     ctree_minbucket = 20, ctree_mincrit = 0.9,
                     cforest_minbucket = 20, cforest_mincrit = 0)
  
  
  return(simlist)
}
  



## HCL palette
pal <- hcl(c(10, 128, 260), 100, 50)
names(pal) <- c("forest", "tree", "gamlss")

pallight <- hcl(c(10, 128, 260), 100, 50, alpha = 0.25)
names(pallight) <- c("forest", "tree", "gamlss")

transpgrey <- rgb(0.190,0.190,0.190, alpha = 0.2)


# plots for various data sets
# results given in simlist
# plot RMSE
plot_rmse <- function(simlist, type = c("exp", "par")){
  
  if(type == "exp"){
    #dev.off() 
    rmse <- cbind(simlist$av.rmse.exp.true.dt, simlist$av.rmse.exp.true.df, simlist$av.rmse.exp.true.g,
                  simlist$av.rmse.exp.obs.dt, simlist$av.rmse.exp.obs.df, simlist$av.rmse.exp.obs.g)
    colnames(rmse) <- c("dt.true", "df.true", "g.true", "dt.obs", "df.obs", "g.obs")
    ylim <- c(min(rmse), max(rmse))
    par(mar=c(5.1,4.1,4.1,3.1))
    plot(x = simlist$x.axis, y = rmse[,"dt.true"], type = "l", col = pal["tree"], ylim = ylim,
         xlab = "", ylab = "", xaxt="n", yaxt="n", lwd = 2)
    title(main = "disttree", col.main = pal["tree"], cex.main = 1.2, font.main = 2, adj = 0.19)
    title(main = "distforest", col.main = pal["forest"], cex.main = 1.2, font.main = 2, adj = 0.5)
    title(main = "gamlss", col.main = pal["gamlss"], cex.main = 1.2, font.main = 2, adj = 0.8)
    title(main = "vs.", col.main = "black", cex.main = 1.2, font.main = 2, adj = 0.35)
    title(main = "vs.", col.main = "black", cex.main = 1.2, font.main = 2, adj = 0.65)
    lines(x = simlist$x.axis, y = rmse[,"dt.obs"], type = "l", lty = 2, col = pal["tree"], lwd = 2)
    lines(x = simlist$x.axis, y = rmse[,"df.true"], type = "l", col = pal["forest"], lwd = 2)
    lines(x = simlist$x.axis, y = rmse[,"df.obs"], type = "l", lty = 2, col = pal["forest"], lwd = 2)
    lines(x = simlist$x.axis, y = rmse[,"g.true"], type = "l", col = pal["gamlss"], lwd = 2)
    lines(x = simlist$x.axis, y = rmse[,"g.obs"], type = "l", lty = 2, col = pal["gamlss"], lwd = 2)
    #legend('left', c("disttree", "distforest", "gamlss"), col = c(pal["tree"], pal["forest"], pal["gamlss"]), lty = 1)
    mtext(text= "smooth", side = 1, line = 1, adj = 0)
    mtext(text= "steep", side = 1, line = 1, adj = 1)
    mtext(text= TeX('$\\kappa$'), side = 1, line = 1)
    mtext(text= "RMSE", side = 2, line = 1)
    mtext(text= "disttree", side = 2, col = pal["tree"], las = 1, line = 0.2, padj = 5)
    mtext(text= "distforest", side = 2, col = pal["forest"], las = 1, line = 0.2, padj = 13)
    mtext(text= "gamlss", side = 2, col = pal["gamlss"], las = 1, line = 0.3, padj = 16.5)
    
    
    # add arrows
    parsave <- par(new = TRUE,  mar = c(0,0,0,0))
    plot(0,0,xlim=c(0,1),ylim=c(0,1),type="n",xlab='', ylab='', col='white', axes = FALSE) 
    
    #add arrow on the right side (downwards, for RMSE)
    arrows(0.99,0.47, 0.99,0.3,lwd=1) ## add arrow
    text(0.99,0.5, 'better')
    segments(0.99,0.53,0.99,0.7,lwd=1)
    
    #add arrows below
    arrows(0.45,0.085,0.3,0.085,lwd=1) ## add arrow
    arrows(0.6,0.085,0.75,0.085,lwd=1) ## add arrow
    
    # plot parameters  reset to prior values
    par(parsave)
  } 
  
  if(type == "par"){
    rmse <- cbind(simlist$av.rmse.mu.dt, simlist$av.rmse.mu.df, simlist$av.rmse.mu.g, 
                  simlist$av.rmse.sigma.dt, simlist$av.rmse.sigma.df, simlist$av.rmse.sigma.g)
    colnames(rmse) <- c("mu.dt", "mu.df", "mu.g", "sigma.dt", "sigma.df", "sigma.g")
    ylim <- c(min(rmse), max(rmse))
    par(mar=c(5.1,4.1,4.1,3.1))
    plot(x = simlist$x.axis, y = rmse[,"mu.dt"], type = "l", col = pal["tree"], ylim = ylim,
         xlab = "", ylab = "", xaxt="n", yaxt="n", lwd = 2)
    title(main = "disttree", col.main = pal["tree"], cex.main = 1.2, font.main = 2, adj = 0.19)
    title(main = "distforest", col.main = pal["forest"], cex.main = 1.2, font.main = 2, adj = 0.5)
    title(main = "gamlss", col.main = pal["gamlss"], cex.main = 1.2, font.main = 2, adj = 0.8)
    title(main = "vs.", col.main = "black", cex.main = 1.2, font.main = 2, adj = 0.35)
    title(main = "vs.", col.main = "black", cex.main = 1.2, font.main = 2, adj = 0.65)
    lines(x = simlist$x.axis, y = rmse[,"sigma.dt"], type = "l", lty = 2, col = pal["tree"], lwd = 2)
    lines(x = simlist$x.axis, y = rmse[,"mu.df"], type = "l", col = pal["forest"], lwd = 2)
    lines(x = simlist$x.axis, y = rmse[,"sigma.df"], type = "l", lty = 2, col = pal["forest"], lwd = 2)
    lines(x = simlist$x.axis, y = rmse[,"mu.g"], type = "l", col = pal["gamlss"], lwd = 2)
    lines(x = simlist$x.axis, y = rmse[,"sigma.g"], type = "l", lty = 2, col = pal["gamlss"], lwd = 2)
    #legend('topleft', c("disttree", "distforest", "gamlss"), col = c(pal["tree"], pal["forest"],pal["gamlss"]), lty = 1)
    mtext(text= "smooth", side = 1, line = 1, adj = 0)
    mtext(text= "steep", side = 1, line = 1, adj = 1)
    mtext(text= TeX('$\\kappa$'), side = 1, line = 1)
    mtext(text= "RMSE", side = 2, line = 1)
    mtext(text= "disttree", side = 2, col = pal["tree"], las = 1, line = 0.2, padj = -7)
    mtext(text= "distforest", side = 2, col = pal["forest"], las = 1, line = 0.2, padj = 8)
    mtext(text= "gamlss", side = 2, col = pal["gamlss"], las = 1, line = 0.3, padj = 14)
    
    # add arrows
    parsave <- par(new = TRUE,  mar = c(0,0,0,0))
    plot(0,0,xlim=c(0,1),ylim=c(0,1),type="n",xlab='',ylab='',col='white', axes = FALSE) 
    
    #add arrow on the right side (downwards, for RMSE)
    arrows(0.99,0.47,0.99,0.3,lwd=1) ## add arrow
    text(0.99,0.5, 'better')
    segments(0.99,0.53,0.99,0.7,lwd=1)
    
    #add arrows below
    arrows(0.45,0.085,0.3,0.085,lwd=1) ## add arrow
    arrows(0.6,0.085,0.75,0.085,lwd=1) ## add arrow
    
    # plot parameters  reset to prior values
    par(parsave)
  }
}

# plot loglikelihood
plot_ll <- function(simlist){
  ll <- cbind(simlist$av.loglik.dt, simlist$av.loglik.df, simlist$av.loglik.g)
  colnames(ll) <- c("dt", "df","g")
  ylim <- c(min(ll), max(ll))
  par(mar=c(5.1,4.1,4.1,3.1))
  plot(x = simlist$x.axis, y = ll[,"dt"], type = "l", col = pal["tree"], ylim = ylim,
       xlab = "", ylab = "", xaxt="n", yaxt="n", lwd = 2)
  title(main = "disttree", col.main = pal["tree"], cex.main = 1.2, font.main = 2, adj = 0.19)
  title(main = "distforest", col.main = pal["forest"], cex.main = 1.2, font.main = 2, adj = 0.5)
  title(main = "gamlss", col.main = pal["gamlss"], cex.main = 1.2, font.main = 2, adj = 0.8)
  title(main = "vs.", col.main = "black", cex.main = 1.2, font.main = 2, adj = 0.35)
  title(main = "vs.", col.main = "black", cex.main = 1.2, font.main = 2, adj = 0.65)
  
  lines(x = simlist$x.axis, y = ll[,"df"], type = "l", col = pal["forest"], lwd = 2)
  lines(x = simlist$x.axis, y = ll[,"g"], type = "l", col = pal["gamlss"], lwd = 2)
  #legend('topright', c("disttree", "distforest", "gamlss"), col = c(pal["tree"], pal["forest"],pal["gamlss"]), lty = 1)
  #mtext(text= "disttree", side = 3, col = pal["tree"], line = 1, adj = 0.3)
  #mtext(text= "distforest", side = 3, col = pal["forest"], line = 1, adj = 0.5)
  #mtext(text= "gamlss", side = 3, col = pal["gamlss"], line = 1, adj = 0.7)
  #mtext(text= "vs.", side = 3, col = "black", line = 1, adj = 0.4)
  #mtext(text= "vs.", side = 3, col = "black", line = 1, adj = 0.6)
  mtext(text= "smooth", side = 1, line = 1, adj = 0)
  mtext(text= "steep", side = 1, line = 1, adj = 1)
  mtext(text= TeX('$\\kappa$'), side = 1, line = 1)
  mtext(text= "Log-Likelihood", side = 2, line = 1)
  mtext(text= "disttree", side = 2, col = pal["tree"], las = 1, line = 0.2, padj = 16.7)
  mtext(text= "distforest", side = 2, col = pal["forest"], las = 1, line = 0.2, padj = -8)
  mtext(text= "gamlss", side = 2, col = pal["gamlss"], las = 1, line = 0.3, padj = -16)
  
  # add arrows
  parsave <- par(new = TRUE,  mar = c(0,0,0,0))
  plot(0,0,xlim=c(0,1),ylim=c(0,1),type="n",xlab='', ylab='', col='white', axes = FALSE) 
  
  #add arrow on the right side (upwards, for LL)
  arrows(0.99,0.53,0.99,0.7,lwd=1) ## add arrow
  text(0.99,0.5, 'better')
  segments(0.99,0.3,0.99,0.47,lwd=1)
  
  #add arrows below
  arrows(0.45,0.085,0.3,0.085,lwd=1) ## add arrow
  arrows(0.6,0.085,0.75,0.085,lwd=1) ## add arrow
  
  # plot parameters  reset to prior values
  par(parsave)
}






if(FALSE){
  
  simlist <- gensimlist(seed = 9, nrep = 70, ntree = 100, covariate = c("x1", "x2"))
  
  
  
  
  
  
  
  
  plot_ll(simlist)
  plot_rmse(simlist, type = "par")
  plot_rmse(simlist, type = "exp")
  
  
}