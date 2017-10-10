## IDEA: learn models on the first 29 years and evaluate it on the 30th year
# average the resulting CRPS and loglikelihood over the 30 days within one station
# and repeat this procedure over all stations (with a sufficient number of years of observations)

rain_axams_pred <- function(seedconst = 7, ntree = 100,
                      tree_minsplit = 50, tree_minbucket = 20, tree_mincrit = 0.95,
                      forest_minsplit = 50, forest_minbucket = 20, forest_mincrit = 0,
                      forest_mtry = 27,   # if frac == TRUE: nvar/3 = 30
                      type.tree = "ctree",
                      gamboost_cvr = FALSE,
                      frac = FALSE)
{
  
  cl <- match.call()
  
  library("disttree")
  library("gamlss")
  library("gamlss.dist")
  library("gamlss.cens")
  gen.cens(NO, type = "left")
  library("gamboostLSS")
  library("crch")
  library("scoringRules")
  library("parallel")
  
  
  #setwd("~/svn/partykit/pkg/disttree/inst/draft/")
  #setwd("~/disttree/inst/draft/")
  
  
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
  
  
  
  
  # load observations
  #load("~/svn/partykit/pkg/disttree/inst/draft/rainData/rain.rda")
  #load("~/disttree/inst/draft/rainData/rain.rda")
  load("rainData/rain.rda")
  #head(data)
  Sys.setenv("TZ" = "UTC")
  rain  <- data
  #dim(rain)
  rain$robs <- (rain$obs)^(1/1.6)
  rain$day <- as.POSIXlt(rain$date)$yday
  rain$year <- as.POSIXlt(rain$date)$year
  rain$hour <- as.POSIXlt(rain$date)$hour
  #head(rain)
  ## save complete data
  rain_all <- rain

  
  
  ############
  # formula
  if(frac){
    
    # tree and forest formula
    dt.formula <- df.formula <- robs ~ tppow_mean + tppow_sprd + tppow_min + tppow_max + 
      #tp_frac +                                                       # include? (mainly zeros), or just the 6h values?
      tp_frac0612 + tp_frac1218 + tp_frac1824 + tp_frac2430 + 
      tppow_mean0612 + tppow_mean1218 + tppow_mean1824 + tppow_mean2430 + 
      tppow_sprd0612 + tppow_sprd1218 + tppow_sprd1824 + tppow_sprd2430 + 
      capepow_mean + capepow_sprd + capepow_min + capepow_max + 
      #cape_frac +                                                    # include? (mainly zeros), or just the 6h values?
      cape_frac0612 + cape_frac1218 + cape_frac1824 + cape_frac2430 + 
      capepow_mean0612 + capepow_mean1218 + capepow_mean1224 + capepow_mean1230 +
      capepow_sprd0612 + capepow_sprd1218 + capepow_sprd1224 + capepow_sprd1230 +
      dswrf_mean_mean + dswrf_mean_max + #dswrf_mean_min + 
      dswrf_sprd_mean + dswrf_sprd_max + #dswrf_sprd_min +
      msl_mean_mean + msl_mean_min + msl_mean_max + 
      msl_sprd_mean + msl_sprd_min + msl_sprd_max +
      pwat_mean_mean + pwat_mean_min + pwat_mean_max + 
      pwat_sprd_mean + pwat_sprd_min + pwat_sprd_max +
      tmax_mean_mean + tmax_mean_min + tmax_mean_max +
      tmax_sprd_mean + tmax_sprd_min + tmax_sprd_max +
      tcolc_mean_mean + tcolc_mean_min + tcolc_mean_max +
      tcolc_sprd_mean + tcolc_sprd_min + tcolc_sprd_max +
      t500_mean_mean + t500_mean_min + t500_mean_max +
      t700_mean_mean + t700_mean_min + t700_mean_max +
      t850_mean_mean + t850_mean_min + t850_mean_max +
      t500_sprd_mean + t500_sprd_min + t500_sprd_max +
      t700_sprd_mean + t700_sprd_min + t700_sprd_max +
      t850_sprd_mean + t850_sprd_min + t850_sprd_max +
      tdiff500850_mean + tdiff500850_min + tdiff500850_max +
      tdiff700850_mean + tdiff700850_min + tdiff700850_max +
      tdiff500700_mean + tdiff500700_min + tdiff500700_max + 
      msl_diff
    
    # gamlss formula
    g.mu.formula <- robs ~ pb(tppow_mean) + 
      pb(tppow_mean1218 * capepow_mean1218) + 
      pb(tppow_max) + 
      pb(dswrf_mean_mean) +
      pb(tcolc_mean_mean) + 
      pb(msl_diff) + 
      pb(pwat_mean_mean) + 
      pb(tdiff500850_mean) #+  pb(tdiff700850_mean) 
    
    g.sigma.formula <- ~ pb(tppow_sprd) + 
      pb(tppow_sprd1218 * capepow_mean1218) + 
      pb(dswrf_sprd_mean) +
      pb(tcolc_sprd_mean) + 
      pb(tdiff500850_mean) #+  pb(tdiff700850_mean) 
    
    
    
    # gamboostLSS formula
    gb.mu.formula <- gb.sigma.formula <- 
      robs ~ bbs(tppow_mean) + bbs(tppow_sprd) + bbs(tppow_min) + bbs(tppow_max) + 
      #bbs(tp_frac) +                                                   # include? (mainly zeros), or just the 6h values?
      bbs(tp_frac0612) + bbs(tp_frac1218) + bbs(tp_frac1824) + bbs(tp_frac2430) +
      bbs(tppow_mean0612) + bbs(tppow_mean1218) + bbs(tppow_mean1824) + bbs(tppow_mean2430) + 
      bbs(tppow_sprd0612) + bbs(tppow_sprd1218) + bbs(tppow_sprd1824) + bbs(tppow_sprd2430) +
      bbs(capepow_mean) + bbs(capepow_sprd) + bbs(capepow_min) + bbs(capepow_max) + 
      #bbs(cape_frac) +                                                 # include? (mainly zeros), or just the 6h values?
      bbs(cape_frac0612) + bbs(cape_frac1218) + bbs(cape_frac1824) + bbs(cape_frac2430) +
      bbs(capepow_mean0612) + bbs(capepow_mean1218) + bbs(capepow_mean1224) + bbs(capepow_mean1230) +
      bbs(capepow_sprd0612) + bbs(capepow_sprd1218) + bbs(capepow_sprd1224) + bbs(capepow_sprd1230) +
      bbs(dswrf_mean_mean) + bbs(dswrf_mean_max) +  #bbs(dswrf_mean_min) +
      bbs(dswrf_sprd_mean) + bbs(dswrf_sprd_max) + #bbs(dswrf_sprd_min) +
      bbs(msl_mean_mean) + bbs(msl_mean_min) + bbs(msl_mean_max) + 
      bbs(msl_sprd_mean) + bbs(msl_sprd_min) + bbs(msl_sprd_max) +
      bbs(pwat_mean_mean) + bbs(pwat_mean_min) + bbs(pwat_mean_max) + 
      bbs(pwat_sprd_mean) + bbs(pwat_sprd_min) + bbs(pwat_sprd_max) +
      bbs(tmax_mean_mean) + bbs(tmax_mean_min) + bbs(tmax_mean_max) +
      bbs(tmax_sprd_mean) + bbs(tmax_sprd_min) + bbs(tmax_sprd_max) +
      bbs(tcolc_mean_mean) + bbs(tcolc_mean_min) + bbs(tcolc_mean_max) +
      bbs(tcolc_sprd_mean) + bbs(tcolc_sprd_min) + bbs(tcolc_sprd_max) +
      bbs(t500_mean_mean) + bbs(t500_mean_min) + bbs(t500_mean_max) +
      bbs(t700_mean_mean) + bbs(t700_mean_min) + bbs(t700_mean_max) +
      bbs(t850_mean_mean) + bbs(t850_mean_min) + bbs(t850_mean_max) +
      bbs(t500_sprd_mean) + bbs(t500_sprd_min) + bbs(t500_sprd_max) +
      bbs(t700_sprd_mean) + bbs(t700_sprd_min) + bbs(t700_sprd_max) +
      bbs(t850_sprd_mean) + bbs(t850_sprd_min) + bbs(t850_sprd_max) +
      bbs(tdiff500850_mean) + bbs(tdiff500850_min) + bbs(tdiff500850_max) +
      bbs(tdiff700850_mean) + bbs(tdiff700850_min) + bbs(tdiff700850_max) +
      bbs(tdiff500700_mean) + bbs(tdiff500700_min) + bbs(tdiff500700_max) + 
      bbs(msl_diff)
    
    
  } else {
    
    # tree and forest formula
    dt.formula <- df.formula <- 
      robs ~ tppow_mean + tppow_sprd + tppow_min + tppow_max + 
      #tp_frac +                                                       # include? (mainly zeros), or just the 6h values?
      #tp_frac0612 + tp_frac1218 + tp_frac1824 + tp_frac2430 + 
      tppow_mean0612 + tppow_mean1218 + tppow_mean1824 + tppow_mean2430 + 
      tppow_sprd0612 + tppow_sprd1218 + tppow_sprd1824 + tppow_sprd2430 + 
      capepow_mean + capepow_sprd + capepow_min + capepow_max + 
      #cape_frac +                                                    # include? (mainly zeros), or just the 6h values?
      #cape_frac0612 + cape_frac1218 + cape_frac1824 + cape_frac2430 + 
      capepow_mean0612 + capepow_mean1218 + capepow_mean1224 + capepow_mean1230 +
      capepow_sprd0612 + capepow_sprd1218 + capepow_sprd1224 + capepow_sprd1230 +
      dswrf_mean_mean + dswrf_mean_max + #dswrf_mean_min + 
      dswrf_sprd_mean + dswrf_sprd_max + #dswrf_sprd_min +
      msl_mean_mean + msl_mean_min + msl_mean_max + 
      msl_sprd_mean + msl_sprd_min + msl_sprd_max +
      pwat_mean_mean + pwat_mean_min + pwat_mean_max + 
      pwat_sprd_mean + pwat_sprd_min + pwat_sprd_max +
      tmax_mean_mean + tmax_mean_min + tmax_mean_max +
      tmax_sprd_mean + tmax_sprd_min + tmax_sprd_max +
      tcolc_mean_mean + tcolc_mean_min + tcolc_mean_max +
      tcolc_sprd_mean + tcolc_sprd_min + tcolc_sprd_max +
      t500_mean_mean + t500_mean_min + t500_mean_max +
      t700_mean_mean + t700_mean_min + t700_mean_max +
      t850_mean_mean + t850_mean_min + t850_mean_max +
      t500_sprd_mean + t500_sprd_min + t500_sprd_max +
      t700_sprd_mean + t700_sprd_min + t700_sprd_max +
      t850_sprd_mean + t850_sprd_min + t850_sprd_max +
      tdiff500850_mean + tdiff500850_min + tdiff500850_max +
      tdiff700850_mean + tdiff700850_min + tdiff700850_max +
      tdiff500700_mean + tdiff500700_min + tdiff500700_max +
      msl_diff
    
    # gamlss formula
    g.mu.formula <- robs ~ pb(tppow_mean) + 
      pb(tppow_mean1218 * capepow_mean1218) + 
      pb(tppow_max) + 
      pb(dswrf_mean_mean) +
      pb(tcolc_mean_mean) + 
      pb(msl_diff) + 
      pb(pwat_mean_mean) + 
      pb(tdiff500850_mean) #+  pb(tdiff700850_mean) 
    
    g.sigma.formula <- ~ pb(tppow_sprd) + 
      pb(tppow_sprd1218 * capepow_mean1218) + 
      pb(dswrf_sprd_mean) +
      pb(tcolc_sprd_mean) + 
      pb(tdiff500850_mean) #+  pb(tdiff700850_mean) 
    
    
    
    # gamboostLSS formula
    gb.mu.formula <- gb.sigma.formula <- 
      robs ~ bbs(tppow_mean) + bbs(tppow_sprd) + bbs(tppow_min) + bbs(tppow_max) + 
      #bbs(tp_frac) +                                                   # include? (mainly zeros), or just the 6h values?
      #bbs(tp_frac0612) + bbs(tp_frac1218) + bbs(tp_frac1824) + bbs(tp_frac2430) +
      bbs(tppow_mean0612) + bbs(tppow_mean1218) + bbs(tppow_mean1824) + bbs(tppow_mean2430) + 
      bbs(tppow_sprd0612) + bbs(tppow_sprd1218) + bbs(tppow_sprd1824) + bbs(tppow_sprd2430) +
      bbs(capepow_mean) + bbs(capepow_sprd) + bbs(capepow_min) + bbs(capepow_max) + 
      #bbs(cape_frac) +                                                 # include? (mainly zeros), or just the 6h values?
      #bbs(cape_frac0612) + bbs(cape_frac1218) + bbs(cape_frac1824) + bbs(cape_frac2430) +
      bbs(capepow_mean0612) + bbs(capepow_mean1218) + bbs(capepow_mean1224) + bbs(capepow_mean1230) +
      bbs(capepow_sprd0612) + bbs(capepow_sprd1218) + bbs(capepow_sprd1224) + bbs(capepow_sprd1230) +
      bbs(dswrf_mean_mean) + bbs(dswrf_mean_max) +  #bbs(dswrf_mean_min) +
      bbs(dswrf_sprd_mean) + bbs(dswrf_sprd_max) + #bbs(dswrf_sprd_min) +
      bbs(msl_mean_mean) + bbs(msl_mean_min) + bbs(msl_mean_max) + 
      bbs(msl_sprd_mean) + bbs(msl_sprd_min) + bbs(msl_sprd_max) +
      bbs(pwat_mean_mean) + bbs(pwat_mean_min) + bbs(pwat_mean_max) + 
      bbs(pwat_sprd_mean) + bbs(pwat_sprd_min) + bbs(pwat_sprd_max) +
      bbs(tmax_mean_mean) + bbs(tmax_mean_min) + bbs(tmax_mean_max) +
      bbs(tmax_sprd_mean) + bbs(tmax_sprd_min) + bbs(tmax_sprd_max) +
      bbs(tcolc_mean_mean) + bbs(tcolc_mean_min) + bbs(tcolc_mean_max) +
      bbs(tcolc_sprd_mean) + bbs(tcolc_sprd_min) + bbs(tcolc_sprd_max) +
      bbs(t500_mean_mean) + bbs(t500_mean_min) + bbs(t500_mean_max) +
      bbs(t700_mean_mean) + bbs(t700_mean_min) + bbs(t700_mean_max) +
      bbs(t850_mean_mean) + bbs(t850_mean_min) + bbs(t850_mean_max) +
      bbs(t500_sprd_mean) + bbs(t500_sprd_min) + bbs(t500_sprd_max) +
      bbs(t700_sprd_mean) + bbs(t700_sprd_min) + bbs(t700_sprd_max) +
      bbs(t850_sprd_mean) + bbs(t850_sprd_min) + bbs(t850_sprd_max) +
      bbs(tdiff500850_mean) + bbs(tdiff500850_min) + bbs(tdiff500850_max) +
      bbs(tdiff700850_mean) + bbs(tdiff700850_min) + bbs(tdiff700850_max) +
      bbs(tdiff500700_mean) + bbs(tdiff500700_min) + bbs(tdiff500700_max) +
      bbs(msl_diff)
    
  }
  
  
  
  ## choose one month:
  # July
  rain <- rain[(((rain$year %% 4 == 0) & rain$day>=183) | (!(rain$year %% 4 == 0) & rain$day>=182)),]
  rain <- rain[(((rain$year %% 4 == 0) & rain$day<=213) | (!(rain$year %% 4 == 0) & rain$day<=212)),]
  
  # predictions start in 1985
  rain <- rain[rain$year>=85,]
  
  #rain_all_July85 <- rain
  
  
  ########################################
  # only one station: Axams (EHYD)  9103309
               
  stationname <- "Axams (EHYD)"
  stationnr <- 9103309
                        
  # load predictions
  #load(paste0("~/svn/partykit/pkg/disttree/inst/draft/rainData/prepared/GEFSV2_prepared_", stationnr, ".rda"))
  #load(paste0("~/disttree/inst/draft/rainData/prepared/GEFSV2_prepared_", stationnr, ".rda"))
  load(paste0("rainData/prepared/GEFSV2_prepared_", stationnr, ".rda", sep = ""))
  prediction <- prepared
  Sys.setenv("TZ" = "UTC")
  prediction$day <-as.POSIXlt(prediction$init)$yday
  prediction$year <- as.POSIXlt(prediction$init)$year
  #head(prediction)
                        
  rain <- rain[rain$station == stationname,]
                        
  # observations end in 2012
  prediction <- prediction[prediction$year<=112,]
  
  # combine data frames
  
  # predictions with init 07-01 are for day 07-02 
  # observations are the sum of precipitation of the last 24h
  # match prediction with init 07-01 with observation on 07-02 
  # (pred starts at 00UTC and predicts from 06UTC until 30UTC = 06UTC of next day, 
  # observations are meassured from 06UTC of previous day to 06UTC of current day)
  
  raindata <- cbind(rain$date, rain$obs, rain$robs, prediction)
  #head(raindata[,c(1:10)])
  colnames(raindata)[c(1:3)] <- c("date", "obs", "robs")
                        
  
  #remove row 825
  raindata <- raindata[-825,]
                        

                        
  # table(raindata[, "year"])
  
  ######################################################
  # only keep variables with sufficient values
  raindata <- raindata[, c("robs", "year",
                           "tppow_mean", "tppow_sprd", "tppow_min", "tppow_max", 
                           "tppow_mean0612", "tppow_mean1218", "tppow_mean1824", "tppow_mean2430", 
                           "tppow_sprd0612", "tppow_sprd1218", "tppow_sprd1824", "tppow_sprd2430",
                           "tp_frac", 
                           "tp_frac0612", "tp_frac1218", "tp_frac1824", "tp_frac2430", 
                           "capepow_mean", "capepow_sprd", "capepow_min", "capepow_max", 
                           "capepow_mean0612", "capepow_mean1218", "capepow_mean1224", "capepow_mean1230",
                           "capepow_sprd0612", "capepow_sprd1218", "capepow_sprd1224", "capepow_sprd1230",
                           "cape_frac",
                           "cape_frac0612", "cape_frac1218", "cape_frac1824", "cape_frac2430", 
                           "dswrf_mean_mean", "dswrf_mean_min", "dswrf_mean_max",
                           "dswrf_sprd_mean", "dswrf_sprd_min", "dswrf_sprd_max",
                           "msl_mean_mean", "msl_mean_min", "msl_mean_max",
                           "msl_sprd_mean", "msl_sprd_min", "msl_sprd_max",
                           "pwat_mean_mean", "pwat_mean_min", "pwat_mean_max",
                           "pwat_sprd_mean", "pwat_sprd_min", "pwat_sprd_max",
                           "tcolc_mean_mean", "tcolc_mean_min", "tcolc_mean_max",
                           "tcolc_sprd_mean", "tcolc_sprd_min", "tcolc_sprd_max",
                           "tmax_mean_mean", "tmax_mean_min", "tmax_mean_max",
                           "tmax_sprd_mean", "tmax_sprd_min", "tmax_sprd_max",
                           "t500_mean_mean", "t500_mean_min", "t500_mean_max",
                           "t700_mean_mean", "t700_mean_min", "t700_mean_max",
                           "t850_mean_mean", "t850_mean_min", "t850_mean_max",
                           "t500_sprd_mean", "t500_sprd_min", "t500_sprd_max",
                           "t700_sprd_mean", "t700_sprd_min", "t700_sprd_max",
                           "t850_sprd_mean", "t850_sprd_min", "t850_sprd_max",
                           "tdiff500850_mean", "tdiff500850_min", "tdiff500850_max",
                           "tdiff700850_mean", "tdiff700850_min", "tdiff700850_max",
                           "tdiff500700_mean", "tdiff500700_min", "tdiff500700_max")]
  
  
  #new vairable msl_diff
  raindata$msl_diff <- raindata$msl_mean_max - raindata$msl_mean_min

  # table(raindata$year)
  
  
  # learning the models on 29 years and evaluating predictions on the 30th year
  learndata <- raindata[raindata$year < 110,]
  testdata <- raindata[raindata$year %in% c(110, 111, 112),]
  
  
                        
  ##############################################################
  # fitting the models
  set.seed(seedconst)
  
  if(type.tree == "mob"){
    dt <- disttree(dt.formula, 
                   data = learndata, family = dist_list_cens_normal, 
                   censtype = "left", censpoint = 0, type.tree = "mob", 
                   control = mob_control(restart = FALSE, numsplit = "center", 
                                         alpha = 1-tree_mincrit, minsplit = tree_minsplit,
                                         minbucket = tree_minbucket))
  }
  if(type.tree == "ctree"){
    dt <- disttree(dt.formula, 
                   data = learndata, family = dist_list_cens_normal, 
                   censtype = "left", censpoint = 0, type.tree = "ctree", 
                   control = ctree_control(teststat = "quad", testtype = "Bonferroni", intersplit = TRUE,
                                           mincriterion = tree_mincrit, minsplit = tree_minsplit,
                                           minbucket = tree_minbucket))
  }
                        
                        
  if(type.tree == "mob"){
    df <- distforest(df.formula, 
                     data = learndata, family = dist_list_cens_normal, type.tree = "mob", 
                     ntree = ntree, censtype = "left", censpoint = 0,
                     control = mob_control(restart = FALSE, numsplit = "center", 
                                           alpha = 1-forest_mincrit, minsplit = forest_minsplit,
                                           minbucket = forest_minbucket), mtry = forest_mtry)
  }
  if(type.tree == "ctree"){
    df <- distforest(df.formula, 
                     data = learndata, family = dist_list_cens_normal, type.tree = "ctree", 
                     ntree = ntree, censtype = "left", censpoint = 0,
                     control = ctree_control(teststat = "quad", testtype = "Univ", intersplit = TRUE,
                                             mincriterion = forest_mincrit, minsplit = forest_minsplit,
                                             minbucket = forest_minbucket), mtry = forest_mtry)
  }
                        
                        
  g_learndata <- learndata
  g_learndata$robs <- Surv(g_learndata$robs, g_learndata$robs>0, type="left")
  
  g <- try(gamlss(formula = g.mu.formula, sigma.formula = g.sigma.formula, data = g_learndata, 
                  family = cens("NO", type = "left"),
                  control = gamlss.control(n.cyc = 100),
                  i.control = glim.control(cyc = 100, bf.cyc = 100)))
  if(inherits(g, "try-error")) {
    g <- NA
    g_error <- paste0(stationname,"_model")
  } else g_error <- NA
  
  gb <- gamboostLSS(formula = list(mu = gb.mu.formula, sigma = gb.sigma.formula), data = g_learndata, 
                    families = as.families(fname = cens("NO", type = "left")), method = "noncyclic",
                    control = boost_control(mstop = 1000L))
  #gbOOB <- gamboostLSS(formula = list(mu = gb.mu.formula, sigma = gb.sigma.formula), data = g_learndata, 
  #                  families = as.families(fname = cens("NO", type = "left")), method = "noncyclic",
  #                  control = boost_control(mstop = 1000L, risk = "oobag"))
  if(gamboost_cvr){
    grid <- seq(50,1000, by = 25)
    cvr <- cvrisk(gb, grid = grid)
    mstop(gb) <- mstop(cvr)
    cvr_opt <- mstop(cvr) 
  } else cvr_opt <- NA
  
                        
  mi <- try(crch(formula = robs ~ tppow_mean | tppow_sprd, 
                 data = learndata, dist = "gaussian", left = 0, link.scale = "identity"))
  if(inherits(mi, "try-error")) {
    mi <- NA
    mi_error <- paste0(stationname,"_model")
  } else mi_error <- NA
  
  ml <- try(crch(formula = robs ~ tppow_mean | log(tppow_sprd + 0.001), 
                 data = learndata, dist = "gaussian", left = 0, link.scale = "log"))
  if(inherits(ml, "try-error")) {
    ml <- NA
    ml_error <- paste0(stationname,"_model")
  } else ml_error <- NA
  
  mq <- try(crch(formula = robs ~ tppow_mean | I(tppow_sprd^2), 
                 data = learndata, dist = "gaussian", left = 0, link.scale = "quadratic"))
  if(inherits(mq, "try-error")) {
    mq <- NA
    mq_error <- paste0(stationname,"_model")
  } else mq_error <- NA
  
  
                        
                        
  ## get predicted parameter for testdata
  
  # disttree
  pdt <- predict(dt, newdata = testdata, type = "parameter")
  dt_mu <- pdt$mu
  dt_sigma <- pdt$sigma
  dt_exp <- pnorm(dt_mu/dt_sigma) * (dt_mu + dt_sigma * (dnorm(dt_mu/dt_sigma) / pnorm(dt_mu/dt_sigma)))
  ## FIX ME: calculation of dt_exp for sigma set to 0.0001
  # idea: 
  if(any(is.na(dt_exp))){
    dt_exp[dt_sigma <= 0.0002] <- pmax(0, dt_mu[dt_sigma <= 0.0002])    
  }
  
  # distforest
  pdf <- predict(df, newdata = testdata, type = "parameter")
  df_mu <- pdf$mu
  df_sigma <- pdf$sigma
  df_exp <- pnorm(df_mu/df_sigma) * (df_mu + df_sigma * (dnorm(df_mu/df_sigma) / pnorm(df_mu/df_sigma)))
  ## FIX ME: calculation of df_exp for sigma set to 0.0001
  # idea: 
  if(any(is.na(df_exp))){
    df_exp[df_sigma <= 0.0002] <- pmax(0, df_mu[df_sigma <= 0.0002])    
  }
  
  # gamlss
  if(!(all(is.na(g)))){
    g_mu <- try(predict(g, newdata = testdata, what = "mu", type = "response", data = g_learndata))
    g_sigma <- try(predict(g, newdata = testdata, what = "sigma", type = "response", data = g_learndata))
    if(inherits(g_mu, "try-error") | inherits(g_sigma, "try-error")) {
      g_mu <- g_sigma <- g_exp <- NA
      g_error <- paste0(stationname,"_pred")
    } else g_exp <- pnorm(g_mu/g_sigma) * (g_mu + g_sigma * (dnorm(g_mu/g_sigma) / pnorm(g_mu/g_sigma)))
  } else g_mu <- g_sigma <- g_exp <- NA
  g_na <- any(c(all(is.na(g)), all(is.na(g_mu)), all(is.na(g_sigma))))
  
  # gamboostLSS
  pgb <- predict(gb, newdata = testdata, parameter = c("mu","sigma"), type = "response")
  gb_mu <- pgb$mu
  gb_sigma <- pgb$sigma
  gb_exp <- pnorm(gb_mu/gb_sigma) * (gb_mu + gb_sigma * (dnorm(gb_mu/gb_sigma) / pnorm(gb_mu/gb_sigma)))
  # pgbOOB <- predict(gbOOB, newdata = testdata, parameter = c("mu","sigma"), type = "response")
  # gbOOB_mu <- pgbOOB$mu
  # gbOOB_sigma <- pgbOOB$sigma
  # gbOOB_exp <- pnorm(gbOOB_mu/gbOOB_sigma) * (gbOOB_mu + gbOOB_sigma * (dnorm(gbOOB_mu/gbOOB_sigma) / pnorm(gbOOB_mu/gbOOB_sigma)))
  
  
  # EMOS
  if(!(all(is.na(mi)))){
    mi_mu <- try(predict(mi, type = "location", newdata = testdata))     # returns parameter on response scale
    mi_sigma <- try(predict(mi, type = "scale", newdata = testdata))
    if(inherits(mi_mu, "try-error") | inherits(mi_sigma, "try-error")) {
      mi_mu <- mi_sigma <- mi_exp <- NA
      mi_error <- paste0(stationname,"_pred")
    } else mi_exp <- pnorm(mi_mu/mi_sigma) * (mi_mu + mi_sigma * (dnorm(mi_mu/mi_sigma) / pnorm(mi_mu/mi_sigma)))
  } else mi_mu <- mi_sigma <- mi_exp <- NA
  mi_na <- any(c(all(is.na(mi)), all(is.na(mi_mu)), all(is.na(mi_sigma))))
  
  if(!(all(is.na(ml)))){
    ml_mu <- try(predict(ml, type = "location", newdata = testdata))     # returns parameter on response scale
    ml_sigma <- try(predict(ml, type = "scale", newdata = testdata))
    if(inherits(ml_mu, "try-error") | inherits(ml_sigma, "try-error")) {
      ml_mu <- ml_sigma <- ml_exp <- NA
      ml_error <- paste0(stationname,"_pred")
    } else ml_exp <- pnorm(ml_mu/ml_sigma) * (ml_mu + ml_sigma * (dnorm(ml_mu/ml_sigma) / pnorm(ml_mu/ml_sigma)))
  } else ml_mu <- ml_sigma <- ml_exp <- NA
  ml_na <- any(c(all(is.na(ml)), all(is.na(ml_mu)), all(is.na(ml_sigma))))
  
  if(!(all(is.na(mq)))){
    mq_mu <- try(predict(mq, type = "location", newdata = testdata))     # returns parameter on response scale
    mq_sigma <- try(predict(mq, type = "scale", newdata = testdata))
    if(inherits(mq_mu, "try-error") | inherits(mq_sigma, "try-error")) {
      mq_mu <- mq_sigma <- mq_exp <- NA
      mq_error <- paste0(stationname,"_pred")
    } else mq_exp <- pnorm(mq_mu/mq_sigma) * (mq_mu + mq_sigma * (dnorm(mq_mu/mq_sigma) / pnorm(mq_mu/mq_sigma)))
  } else mq_mu <- mq_sigma <- mq_exp <- NA
  mq_na <- any(c(all(is.na(mq)), all(is.na(mq_mu)), all(is.na(mq_sigma))))
  
                        
  # CPRS
  crps_dt <- mean(crps_cnorm(testdata$robs, location = dt_mu, scale = dt_sigma, lower = 0, upper = Inf), na.rm = TRUE)
  crps_df <- mean(crps_cnorm(testdata$robs, location = df_mu, scale = df_sigma, lower = 0, upper = Inf), na.rm = TRUE)
  crps_g <- if(!g_na) mean(crps_cnorm(testdata$robs, location = g_mu, scale = g_sigma, lower = 0, upper = Inf), na.rm = TRUE) else NA
  crps_gb <- mean(crps_cnorm(testdata$robs, location = gb_mu, scale = gb_sigma, lower = 0, upper = Inf), na.rm = TRUE)
  #crps_gbOOB <- mean(crps_cnorm(testdata$robs, location = gbOOB_mu, scale = gbOOB_sigma, lower = 0, upper = Inf), na.rm = TRUE)
  crps_mi <- if(!mi_na) mean(crps_cnorm(testdata$robs, location = mi_mu, scale = mi_sigma, lower = 0, upper = Inf), na.rm = TRUE) else NA
  crps_ml <- if(!ml_na) mean(crps_cnorm(testdata$robs, location = ml_mu, scale = ml_sigma, lower = 0, upper = Inf), na.rm = TRUE) else NA
  crps_mq <- if(!mq_na) mean(crps_cnorm(testdata$robs, location = mq_mu, scale = mq_sigma, lower = 0, upper = Inf), na.rm = TRUE) else NA
  
  # RMSE
  rmse_dt <- sqrt(mean((dt_exp - testdata[,"robs"])^2, na.rm = TRUE))
  rmse_df <- sqrt(mean((df_exp - testdata[,"robs"])^2, na.rm = TRUE))
  rmse_g <- if(!g_na) sqrt(mean((g_exp - testdata[,"robs"])^2, na.rm = TRUE)) else NA
  rmse_gb <- sqrt(mean((gb_exp - testdata[,"robs"])^2, na.rm = TRUE))
  #rmse_gbOOB <- sqrt(mean((gbOOB_exp - testdata[,"robs"])^2, na.rm = TRUE))
  rmse_mi <- if(!mi_na) sqrt(mean((mi_exp - testdata[,"robs"])^2, na.rm = TRUE)) else NA
  rmse_ml <- if(!ml_na) sqrt(mean((ml_exp - testdata[,"robs"])^2, na.rm = TRUE)) else NA
  rmse_mq <- if(!mq_na) sqrt(mean((mq_exp - testdata[,"robs"])^2, na.rm = TRUE)) else NA
  
  # loglikelihood
  dtll <- dfll <- gll <-  gbll <- mill <- mlll <- mqll <- numeric(length = NROW(testdata))
  for(j in 1:(NROW(testdata))){
    
    eta_dt <- as.numeric(dist_list_cens_normal$linkfun(cbind(dt_mu, dt_sigma)[j,]))
    eta_df <- as.numeric(dist_list_cens_normal$linkfun(cbind(df_mu, df_sigma)[j,]))
    eta_g <- if(!g_na) as.numeric(dist_list_cens_normal$linkfun(cbind(g_mu, g_sigma)[j,])) else NA
    eta_gb <- as.numeric(dist_list_cens_normal$linkfun(cbind(gb_mu, gb_sigma)[j,]))
    eta_mi <- if(!mi_na) as.numeric(dist_list_cens_normal$linkfun(cbind(mi_mu, mi_sigma)[j,])) else NA
    eta_ml <- if(!ml_na) as.numeric(dist_list_cens_normal$linkfun(cbind(ml_mu, ml_sigma)[j,])) else NA
    eta_mq <- if(!mq_na) as.numeric(dist_list_cens_normal$linkfun(cbind(mq_mu, mq_sigma)[j,])) else NA
    
    dtll[j] <- dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_dt, log=TRUE)
    dfll[j] <- dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_df, log=TRUE)
    gll[j] <- if(!g_na) dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_g, log=TRUE) else NA
    gbll[j] <- dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_gb, log=TRUE)
    mill[j] <- if(!mi_na) dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_mi, log=TRUE) else NA
    mlll[j] <- if(!mi_na) dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_ml, log=TRUE) else NA
    mqll[j] <- if(!mi_na) dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_mq, log=TRUE) else NA
    
  }
  
  dtll <- mean(dtll, na.rm = TRUE)
  dfll <- mean(dfll, na.rm = TRUE)
  gll <- if(!g_na) mean(gll, na.rm = TRUE) else NA
  gbll <- mean(gbll, na.rm = TRUE)
  mill <- if(!mi_na) mean(mill, na.rm = TRUE) else NA      
  mlll <- if(!ml_na) mean(mlll, na.rm = TRUE) else NA 
  mqll <- if(!mq_na) mean(mqll, na.rm = TRUE) else NA 
  ll <- c(dtll, dfll, gll, gbll, mill, mlll, mqll) 
  rmse <- c(rmse_dt, rmse_df, rmse_g, rmse_gb, rmse_mi, rmse_ml, rmse_mq)
  crps <- c(crps_dt, crps_df, crps_g, crps_gb, crps_mi, crps_ml, crps_mq)
                        
  rainres <- list()
  rainres$results <- rbind(ll, rmse, crps)
  colnames(rainres$results) <-
    c("disttree", "distforest", "gamlss", "gamboostLSS", "EMOS id", "EMOS log", "EMOS quad")
  rownames(rainres$results) <- c("ll", "rmse", "crps")
  
  rainres$cvr_opt <- cvr_opt
  rainres$g_error <- g_error
  rainres$mi_error <- mi_error
  rainres$ml_error <- ml_error
  rainres$mq_error <- mq_error
  rainres$call <- cl
  rainres$dt <- dt
  rainres$df <- df
  rainres$g <- g
  rainres$gb <- gb
  rainres$ml <- ml
  rainres$learndata <- learndata
  rainres$testdata <- testdata
  
  return(rainres)
}                    
  





if(FALSE){
  setwd("~/svn/partykit/pkg/disttree/inst/draft")
  #source("rain_pred.R")
  library("gamlss.cens")
  gen.cens("NO", type = "left")
  
  #save(res, file = "~/svn/partykit/pkg/disttree/inst/draft/rain_pred.rda")
  #save(res, file = "~/disttree/inst/draft/rain_pred.rda")

  res <- rain_axams_pred(seedconst = 7, ntree = 100,
                         forest_mtry = 27,   # if frac == TRUE: nvar/3 = 30
                         gamboost_cvr = TRUE,
                         frac = FALSE)
  
  save(res, file = paste0("rain_Axams_pred_", res$call$frac, "_", res$call$seedconst, ".rda"))
  
  
  load("rain_Axams_pred_FALSE_7.rda")
  
  ###########
  # predictions for one day (in each of the three years) 
  # (19th of July 2011 is missing)
  pday <- 24  # 2 (hohe Beobachtung zu niedrig geschaetzt), 4, 15, evtl. auch 7, 8, 23 (eine 0-Beobachtung und 2 sehr aehnliche), 
  {
    pdays <- if(pday<19) c(pday, pday + 31, pday + 61) else c(pday, pday + 30, pday + 61)
    pdf <- predict(res$df, newdata = res$testdata[pdays,], type = "parameter")
    df_mu <- pdf$mu
    df_sigma <- pdf$sigma
    df_exp <- pnorm(df_mu/df_sigma) * (df_mu + df_sigma * (dnorm(df_mu/df_sigma) / pnorm(df_mu/df_sigma)))
    cbind(df_exp, res$testdata[pdays,"robs"])
    
    
    # plot predicted distributions together with observations
    #set.seed(res$call$seedconst)
    set.seed(7)
    x <- c(0.01, sort(runif(500,0.01,8)))
    y1 <- crch::dcnorm(x, mean = df_mu[1], sd = df_sigma[1], left = 0)
    y2 <- crch::dcnorm(x, mean = df_mu[2], sd = df_sigma[2], left = 0)
    y3 <- crch::dcnorm(x, mean = df_mu[3], sd = df_sigma[3], left = 0)
    dayending <- if(pday > 3) "th" else switch(pday, "1" = {dayending <- "st"}, "2" = {dayending <- "nd"}, "3" = {dayending <- "rd"})
    
    # point mass (slightly shifted)
    pm1 <- c(0.04, crch::dcnorm(-1, mean = df_mu[1], sd = df_sigma[1], left = 0))
    pm2 <- c(-0.03, crch::dcnorm(-1, mean = df_mu[2], sd = df_sigma[2], left = 0))
    pm3 <- c(-0.1, crch::dcnorm(-1, mean = df_mu[3], sd = df_sigma[3], left = 0))
    
    # predictions
    pred1 <- c(res$testdata[pdays,"robs"][1], crch::dcnorm(res$testdata[pdays,"robs"][1], mean = df_mu[1], sd = df_sigma[1], left = 0))
    pred2 <- c(res$testdata[pdays,"robs"][2], crch::dcnorm(res$testdata[pdays,"robs"][2], mean = df_mu[2], sd = df_sigma[2], left = 0))
    pred3 <- c(res$testdata[pdays,"robs"][3], crch::dcnorm(res$testdata[pdays,"robs"][3], mean = df_mu[3], sd = df_sigma[3], left = 0))
    
    #legendheight
    lh1 <- crch::dcnorm(0.01, mean = df_mu[1], sd = df_sigma[1], left = 0)
    lh2 <- crch::dcnorm(0.01, mean = df_mu[2], sd = df_sigma[2], left = 0)
    lh3 <- crch::dcnorm(0.01, mean = df_mu[3], sd = df_sigma[3], left = 0)
    
    plot(x = x, y = y1, type = "l", col = "red", 
         main = paste0("July ", pday, dayending), ylab = "density", 
         ylim = c(0,max(y1, y2, y3, pm1, pm2, pm3) + 0.01),
         xlim = c(-1.5,8))
    
    lines(x = x, y = y2, type = "l", col = "blue")
    lines(x = x, y = y3, type = "l", col = "darkgreen")
    #legend('topright', c("2010", "2011", "2012"), col = c("red", "blue", "darkgreen"), lty = 1, cex = 1)
    
    # plot point mass
    lines(x = c(pm1[1], pm1[1]), y = c(pm1[2], 0), col = "red", type = "l", lwd = 1)
    lines(x = c(pm2[1], pm2[1]), y = c(pm2[2], 0), col = "blue", type = "l", lwd = 1)
    lines(x = c(pm3[1], pm3[1]), y = c(pm3[2], 0), col = "darkgreen", type = "l", lwd = 1)
    points(x = pm1[1], y = pm1[2], col = "red", pch = 19)
    points(x = pm2[1], y = pm2[2], col = "blue", pch = 19)
    points(x = pm3[1], y = pm3[2], col = "darkgreen", pch = 19)
    
    
    # plot predictions
    points(x = pred1[1], y = pred1[2], col = "red", pch = 4)
    points(x = pred2[1], y = pred2[2], col = "blue", pch = 4)
    points(x = pred3[1], y = pred3[2], col = "darkgreen", pch = 4)
    
    lines(x = c(pred1[1], pred1[1]), y = c(pred1[2], 0), col = "darkgrey", type = "l", lty = 2)
    lines(x = c(pred2[1], pred2[1]), y = c(pred2[2], 0), col = "darkgrey", type = "l", lty = 2)
    lines(x = c(pred3[1], pred3[1]), y = c(pred3[2], 0), col = "darkgrey", type = "l", lty = 2)
    
    # add labels
    text(x = -0.8, y = lh1, labels = "2010", col = "red", cex = 0.8)
    text(x = -0.8, y = lh2, labels = "2011", col = "blue", cex = 0.8)
    text(x = -0.8, y = lh3, labels = "2012", col = "darkgreen", cex = 0.8)
    
  }
  
  
  
  ###################
  ## pit histograms
  if(FALSE){
    
    ## get predicted parameter
    # for learndata and testdata
    
    # disttree
    pdt <- predict(res$dt, type = "parameter")
    dt_mu_l <- pdt$mu 
    dt_sigma_l <- pdt$sigma 
    pdt <- predict(res$dt, type = "parameter", newdata = res$testdata)
    dt_mu_t <- pdt$mu 
    dt_sigma_t <- pdt$sigma 
    
    #distforest
    pdf <- predict(res$df, type = "parameter")
    df_mu_l <- pdf$mu
    df_sigma_l <- pdf$sigma
    pdf <- predict(res$df, type = "parameter", newdata = res$testdata)
    df_mu_t <- pdf$mu
    df_sigma_t <- pdf$sigma
    
    #gamlss
    g_mu_l <- predict(res$g, what = "mu", type = "response", data = g_raindata)
    g_sigma_l <- predict(res$g, what = "sigma", type = "response", data = g_raindata)
    g_mu_t <- predict(res$g, what = "mu", type = "response", data = g_raindata, newdata = res$testdata)
    g_sigma_t <- predict(res$g, what = "sigma", type = "response", data = g_raindata, newdata = res$testdata)
    
    #gamboostLSS
    pgb <- predict(res$gb, parameter = c("mu","sigma"), type = "response")
    gb_mu_l <- pgb$mu
    gb_sigma_l <- pgb$sigma
    pgb <- predict(res$gb, parameter = c("mu","sigma"), type = "response", newdata = res$testdata)
    gb_mu_t <- pgb$mu
    gb_sigma_t <- pgb$sigma
    
    #EMOS
    ml_mu_l <- predict(res$ml, type = "location")     # returns parameter on response scale
    ml_sigma_l <- predict(res$ml, type = "scale")
    ml_mu_t <- predict(res$ml, type = "location", newdata = res$testdata)     # returns parameter on response scale
    ml_sigma_t <- predict(res$ml, type = "scale", newdata = res$testdata)
    
    
  
    # disttree
    set.seed(7)
    #hist(pnorm(res$learndata[,"robs"], dt_mu_l, dt_sigma_l))
    pit_dt_l <- pnorm(res$learndata[,"robs"], dt_mu_l, dt_sigma_l)
    pit_dt_l[which(res$learndata[,"robs"]==0)] <- pit_dt_l[which(res$learndata[,"robs"]==0)]*runif(length(pit_dt_l[which(res$learndata[,"robs"]==0)]),0,1)
    hist(pit_dt_l)
    #hist(pnorm(res$testdata[,"robs"], dt_mu_t, dt_sigma_t))
    set.seed(7)
    pit_dt_t <- pnorm(res$testdata[,"robs"], dt_mu_t, dt_sigma_t)
    pit_dt_t[which(res$testdata[,"robs"]==0)] <- pit_dt_t[which(res$testdata[,"robs"]==0)]*runif(length(pit_dt_t[which(res$testdata[,"robs"]==0)]),0,1)
    hist(pit_dt_t)
    
    # distforest
    set.seed(7)
    #hist(pnorm(res$learndata[,"robs"], df_mu_l, df_sigma_l))
    pit_df_l <- pnorm(res$learndata[,"robs"], df_mu_l, df_sigma_l)
    pit_df_l[which(res$learndata[,"robs"]==0)] <- pit_df_l[which(res$learndata[,"robs"]==0)]*runif(length(pit_df_l[which(res$learndata[,"robs"]==0)]),0,1)
    hist(pit_df_l)
    set.seed(7)
    #hist(pnorm(res$testdata[,"robs"], df_mu_t, df_sigma_t))
    pit_df_t <- pnorm(res$testdata[,"robs"], df_mu_t, df_sigma_t)
    pit_df_t[which(res$testdata[,"robs"]==0)] <- pit_df_t[which(res$testdata[,"robs"]==0)]*runif(length(pit_df_t[which(res$testdata[,"robs"]==0)]),0,1)
    hist(pit_df_t)
    
    # gamlss
    set.seed(7)
    #hist(pnorm(res$learndata[,"robs"], g_mu_l, g_sigma_l))
    pit_g_l <- pnorm(res$learndata[,"robs"], g_mu_l, g_sigma_l)
    pit_g_l[which(res$learndata[,"robs"]==0)] <- pit_g_l[which(res$learndata[,"robs"]==0)]*runif(length(pit_g_l[which(res$learndata[,"robs"]==0)]),0,1)
    hist(pit_g_l)
    set.seed(7)
    #hist(pnorm(raindata[,"robs"], g_mu_t, g_sigma_t))
    pit_g_t <- pnorm(raindata[,"robs"], g_mu_t, g_sigma_t)
    pit_g_t[which(raindata[,"robs"]==0)] <- pit_g_t[which(raindata[,"robs"]==0)]*runif(length(pit_g_t[which(raindata[,"robs"]==0)]),0,1)
    hist(pit_g_t)
  
    # gamboostLSS
    set.seed(7)
    #hist(pnorm(res$learndata[,"robs"], gb_mu_l, gb_sigma_l))
    pit_gb_l <- pnorm(res$learndata[,"robs"], gb_mu_l, gb_sigma_l)
    pit_gb_l[which(res$learndata[,"robs"]==0)] <- pit_gb_l[which(res$learndata[,"robs"]==0)]*runif(length(pit_gb_l[which(res$learndata[,"robs"]==0)]),0,1)
    hist(pit_gb_l)
    set.seed(7)
    #hist(pnorm(res$testdata[,"robs"], gb_mu_t, gb_sigma_t))
    pit_gb_t <- pnorm(res$testdata[,"robs"], gb_mu_t, gb_sigma_t)
    pit_gb_t[which(res$testdata[,"robs"]==0)] <- pit_gb_t[which(res$testdata[,"robs"]==0)]*runif(length(pit_gb_t[which(res$testdata[,"robs"]==0)]),0,1)
    hist(pit_gb_t)
    
    # EMOS
    set.seed(7)
    #hist(pnorm(res$learndata[,"robs"], ml_mu_l, ml_sigma_l))
    pit_ml_l <- pnorm(res$learndata[,"robs"], ml_mu_l, ml_sigma_l)
    pit_ml_l[which(res$learndata[,"robs"]==0)] <- pit_ml_l[which(res$learndata[,"robs"]==0)]*runif(length(pit_ml_l[which(res$learndata[,"robs"]==0)]),0,1)
    hist(pit_ml_l)
    set.seed(7)
    #hist(pnorm(res$testdata[,"robs"], ml_mu_t, ml_sigma_t))
    pit_ml_t <- pnorm(res$testdata[,"robs"], ml_mu_t, ml_sigma_t)
    pit_ml_t[which(res$testdata[,"robs"]==0)] <- pit_ml_t[which(res$testdata[,"robs"]==0)]*runif(length(pit_ml_t[which(res$testdata[,"robs"]==0)]),0,1)
    hist(pit_ml_t)
    
  }
  
  
  
  ##############################
  # variable importance
  vimp_ll <- varimp(res$df, nperm = 1L)
  vimp10_ll <- varimp(res$df, nperm = 10L)
  vimp100_ll <- varimp(res$df, nperm = 100L)
  vimp1000_ll <- varimp(res$df, nperm = 1000L)
  
  save(vimp_ll, file = "~/svn/partykit/pkg/disttree/inst/draft/vimp_ll.rda")
  save(vimp10_ll, file = "~/svn/partykit/pkg/disttree/inst/draft/vimp10_ll.rda")
  save(vimp100_ll, file = "~/svn/partykit/pkg/disttree/inst/draft/vimp100_ll.rda")
  save(vimp1000_ll, file = "~/svn/partykit/pkg/disttree/inst/draft/vimp1000_ll.rda")
  
  
  # locally redefine logLik.distforest to use crps in varimp()
  logLik.distforest <- function(object, newdata = NULL, ...) {
    if(is.null(newdata)) {
      newdata <- object$data
    } 
    
    # predict parameter
    pdf <- predict(object, newdata = newdata, type = "parameter")
    df_mu <- pdf$mu
    df_sigma <- pdf$sigma
    
    # calculate CRPS
    crps <- mean(crps_cnorm(newdata$robs, location = df_mu, scale = df_sigma, lower = 0, upper = Inf), na.rm = TRUE)
  
    return(structure(crps, df = NA, class = "logLik"))
  }
  
  vimp_crps <- varimp(res$df, nperm = 1L)
  vimp10_crps <- varimp(res$df, nperm = 10L)
  vimp100_crps <- varimp(res$df, nperm = 100L)
  vimp1000_crps <- varimp(res$df, nperm = 1000L)
  
  save(vimp_crps, file = "~/svn/partykit/pkg/disttree/inst/draft/vimp_crps.rda")
  save(vimp10_crps, file = "~/svn/partykit/pkg/disttree/inst/draft/vimp10_crps.rda")
  save(vimp100_crps, file = "~/svn/partykit/pkg/disttree/inst/draft/vimp100_crps.rda")
  save(vimp1000_crps, file = "~/svn/partykit/pkg/disttree/inst/draft/vimp1000_crps.rda")
  
  rm(logLik.distforest)  
  
  
  load("vimp1000_crps.rda")
  
  most_impvar <- sort(vimp1000_crps, decreasing = TRUE)[c(1:10)]
  #most_impvar <- most_impvar[most_impvar>2.7]
  names_mi <- names(most_impvar)
  covnames <- colnames(res$learndata)
  id_mi <- match(names_mi, covnames)
  
  datafix <- res$learndata
  datafix1 <- datafix2 <- datafix3 <- datafix4 <- 
    datafix5 <- datafix5 <- datafix6 <- datafix7 <- datafix8 <- datafix9 <- datafix10 <- datafix
  
  # fix all variables except for the most important ones to their mean
  for(i in 3:NCOL(datafix))  if(i!=id_mi[1]) datafix1[,i] <- rep(mean(datafix[,i], na.rm = TRUE), NROW(datafix))
  for(i in 3:NCOL(datafix))  if(i!=id_mi[2]) datafix2[,i] <- rep(mean(datafix[,i], na.rm = TRUE), NROW(datafix))
  for(i in 3:NCOL(datafix))  if(i!=id_mi[3]) datafix3[,i] <- rep(mean(datafix[,i], na.rm = TRUE), NROW(datafix))
  for(i in 3:NCOL(datafix))  if(i!=id_mi[4]) datafix4[,i] <- rep(mean(datafix[,i], na.rm = TRUE), NROW(datafix))
  for(i in 3:NCOL(datafix))  if(i!=id_mi[5]) datafix5[,i] <- rep(mean(datafix[,i], na.rm = TRUE), NROW(datafix))
  for(i in 3:NCOL(datafix))  if(i!=id_mi[6]) datafix6[,i] <- rep(mean(datafix[,i], na.rm = TRUE), NROW(datafix))
  for(i in 3:NCOL(datafix))  if(i!=id_mi[7]) datafix7[,i] <- rep(mean(datafix[,i], na.rm = TRUE), NROW(datafix))
  for(i in 3:NCOL(datafix))  if(i!=id_mi[8]) datafix8[,i] <- rep(mean(datafix[,i], na.rm = TRUE), NROW(datafix))
  for(i in 3:NCOL(datafix))  if(i!=id_mi[9]) datafix9[,i] <- rep(mean(datafix[,i], na.rm = TRUE), NROW(datafix))
  for(i in 3:NCOL(datafix))  if(i!=id_mi[10]) datafix10[,i] <- rep(mean(datafix[,i], na.rm = TRUE), NROW(datafix))
  
  
  pred1 <- predict(res$df, newdata  = datafix1, type = "parameter")
  plotdata1 <- cbind(datafix1, pred1)
  sp1 <- plotdata1[order(datafix1[,id_mi[1]]),]
  plot(x = sp1[,id_mi[1]], y = sp1[, "mu"], type = "l", xlab = names_mi[1], ylab = "mu")
  plot(x = sp1[,id_mi[1]], y = sp1[, "sigma"], type = "l", xlab = names_mi[1], ylab = "sigma")
  
  pred2 <- predict(res$df, newdata  = datafix2, type = "parameter")
  plotdata2 <- cbind(datafix2, pred2)
  sp2 <- plotdata2[order(datafix2[,id_mi[2]]),]
  plot(x = sp2[,id_mi[2]], y = sp2[, "mu"], type = "l", xlab = names_mi[2], ylab = "mu")
  plot(x = sp2[,id_mi[2]], y = sp2[, "sigma"], type = "l", xlab = names_mi[2], ylab = "sigma")
  
  pred3 <- predict(res$df, newdata  = datafix3, type = "parameter")
  plotdata3 <- cbind(datafix3, pred3)
  sp3 <- plotdata3[order(datafix3[,id_mi[3]]),]
  plot(x = sp3[,id_mi[3]], y = sp3[, "mu"], type = "l", xlab = names_mi[3], ylab = "mu")
  plot(x = sp3[,id_mi[3]], y = sp3[, "sigma"], type = "l", xlab = names_mi[3], ylab = "sigma")
  
  pred4 <- predict(res$df, newdata  = datafix4, type = "parameter")
  plotdata4 <- cbind(datafix4, pred4)
  sp4 <- plotdata4[order(datafix4[,id_mi[4]]),]
  plot(x = sp4[,id_mi[4]], y = sp4[, "mu"], type = "l", xlab = names_mi[4], ylab = "mu")
  plot(x = sp4[,id_mi[4]], y = sp4[, "sigma"], type = "l", xlab = names_mi[4], ylab = "sigma")
  
  pred5 <- predict(res$df, newdata  = datafix5, type = "parameter")
  plotdata5 <- cbind(datafix5, pred5)
  sp5 <- plotdata5[order(datafix5[,id_mi[5]]),]
  plot(x = sp5[,id_mi[5]], y = sp5[, "mu"], type = "l", xlab = names_mi[5], ylab = "mu")
  plot(x = sp5[,id_mi[5]], y = sp5[, "sigma"], type = "l", xlab = names_mi[5], ylab = "sigma")
  
  pred6 <- predict(res$df, newdata  = datafix6, type = "parameter")
  plotdata6 <- cbind(datafix6, pred6)
  sp6 <- plotdata6[order(datafix6[,id_mi[6]]),]
  plot(x = sp6[,id_mi[6]], y = sp6[, "mu"], type = "l", xlab = names_mi[6], ylab = "mu")
  plot(x = sp6[,id_mi[6]], y = sp6[, "sigma"], type = "l", xlab = names_mi[6], ylab = "sigma")
  
  pred7 <- predict(res$df, newdata  = datafix7, type = "parameter")
  plotdata7 <- cbind(datafix7, pred7)
  sp7 <- plotdata7[order(datafix7[,id_mi[7]]),]
  plot(x = sp7[,id_mi[7]], y = sp7[, "mu"], type = "l", xlab = names_mi[7], ylab = "mu")
  plot(x = sp7[,id_mi[7]], y = sp7[, "sigma"], type = "l", xlab = names_mi[7], ylab = "sigma")
  
  pred8 <- predict(res$df, newdata  = datafix8, type = "parameter")
  plotdata8 <- cbind(datafix8, pred8)
  sp8 <- plotdata8[order(datafix8[,id_mi[8]]),]
  plot(x = sp8[,id_mi[8]], y = sp8[, "mu"], type = "l", xlab = names_mi[8], ylab = "mu")
  plot(x = sp8[,id_mi[8]], y = sp8[, "sigma"], type = "l", xlab = names_mi[8], ylab = "sigma")
  
  pred9 <- predict(res$df, newdata  = datafix9, type = "parameter")
  plotdata9 <- cbind(datafix9, pred9)
  sp9 <- plotdata9[order(datafix9[,id_mi[9]]),]
  plot(x = sp9[,id_mi[9]], y = sp9[, "mu"], type = "l", xlab = names_mi[9], ylab = "mu")
  plot(x = sp9[,id_mi[9]], y = sp9[, "sigma"], type = "l", xlab = names_mi[9], ylab = "sigma")
  
  pred10 <- predict(res$df, newdata  = datafix10, type = "parameter")
  plotdata10 <- cbind(datafix10, pred10)
  sp10 <- plotdata10[order(datafix10[,id_mi[10]]),]
  plot(x = sp10[,id_mi[10]], y = sp10[, "mu"], type = "l", xlab = names_mi[10], ylab = "mu")
  plot(x = sp10[,id_mi[10]], y = sp10[, "sigma"], type = "l", xlab = names_mi[10], ylab = "sigma")
  
}


### FIX ME:
# in disttree (and distforest): for nodes with all equal observations (all 0), sigma is set to 0.0001
# this leads to extremely low values in the loglikelihood if the observed value is not 0
# replace? warning?



