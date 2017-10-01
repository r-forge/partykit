## IDEA: learn models on the first 29 years and evaluate it on the 30th year
# average the resulting CRPS and loglikelihood over the 30 days within one station
# and repeat this procedure over all stations (with a sufficient number of years of observations)

rain_pred <- function(seedconst = 7, ntree = 100,
                      tree_minsplit = 50, tree_minbucket = 20, tree_mincrit = 0.95,
                      forest_minsplit = 50, forest_minbucket = 20, forest_mincrit = 0,
                      forest_mtry = 27,
                      type.tree = "ctree",
                      gamboost_cvr = FALSE)
{
  
  cl <- match.call()
  
  library("disttree")
  library("gamlss")
  library("gamlss.dist")
  library("gamlss.cens")
  gen.cens(NO, type = "left")
  library("gamboostLSS")
  library("bamlss")
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
  
  
  # load station list
  #load("~/svn/partykit/pkg/disttree/inst/draft/rainData/ehyd.statlist.rda")
  #load("~/disttree/inst/draft/rainData/ehyd.statlist.rda")
  load("rainData/ehyd.statlist.rda")
  stations <- cbind(ehyd.statlist$station, ehyd.statlist$stationname)
  colnames(stations) <- c("stationnr", "stationname")
  
  
  
  ### check which stations provide observations for the years 1985 - 2012 (both included)
  {
    complete <- list()
    nryears <- numeric(length = NROW(stations))
    all85_112 <- logical(length = NROW(stations))
    all31 <- logical(length = NROW(stations))
    
    for(i in 1:NROW(stations)){
      
      rain <- rain_all
      stationname <- stations[i, "stationname"]
      stationnr <- stations[i, "stationnr"]
      
      # load predictions
      load(paste0("rainData/prepared/GEFSV2_prepared_", stationnr, ".rda", sep = ""))
      prediction <- prepared
      Sys.setenv("TZ" = "UTC")
      prediction$day <-as.POSIXlt(prediction$init)$yday
      prediction$year <- as.POSIXlt(prediction$init)$year
      
      rain <- rain[rain$station == stationname,]
      
      ## choose one month:
      # July
      rain <- rain[(((rain$year %% 4 == 0) & rain$day>=183) | (!(rain$year %% 4 == 0) & rain$day>=182)),]
      rain <- rain[(((rain$year %% 4 == 0) & rain$day<=213) | (!(rain$year %% 4 == 0) & rain$day<=212)),]
      
      all85_112[i] <- all(c(85:112) %in% as.numeric(row.names(table(rain$year)))) 
      all31[i] <- (min(table(rain$year)) == 31) & (max(table(rain$year)) == 31)
      
      complete[[i]] <- table(rain$year)
      #unique(rain$year)
      nryears[i] <- length(unique(rain$year))
    }
    
    complete_stations <- which(all85_112 & all31)
  }
  
  rain <- rain_all
  
  
  ############
  # formula
  {
  
  # tree and forest formula
  dt.formula <- df.formula <- robs ~ tppow_mean + tppow_sprd + tppow_min + tppow_max + #tp_frac + 
    tppow_mean0612 + tppow_mean1218 + tppow_mean1824 + tppow_mean2430 + 
    tppow_sprd0612 + tppow_sprd1218 + tppow_sprd1824 + tppow_sprd2430 + 
    capepow_mean + capepow_sprd + capepow_min + capepow_max + #cape_frac +
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
    tdiff500700_mean + tdiff500700_min + tdiff500700_max
  
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
  gb.mu.formula <- robs ~ bbs(tppow_mean) + bbs(tppow_sprd) + bbs(tppow_min) + bbs(tppow_max) + #bbs(tp_frac) +  
    bbs(tppow_mean0612) + bbs(tppow_mean1218) + bbs(tppow_mean1824) + bbs(tppow_mean2430) + 
    bbs(tppow_sprd0612) + bbs(tppow_sprd1218) + bbs(tppow_sprd1824) + bbs(tppow_sprd2430) +
    bbs(capepow_mean) + bbs(capepow_sprd) + bbs(capepow_min) + bbs(capepow_max) + #bbs(cape_frac) +
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
    bbs(tdiff500700_mean) + bbs(tdiff500700_min) + bbs(tdiff500700_max)
  
  gb.sigma.formula <- robs ~ bbs(tppow_mean) + bbs(tppow_sprd) + bbs(tppow_min) + bbs(tppow_max) + #bbs(tp_frac) +  
    bbs(tppow_mean0612) + bbs(tppow_mean1218) + bbs(tppow_mean1824) + bbs(tppow_mean2430) + 
    bbs(tppow_sprd0612) + bbs(tppow_sprd1218) + bbs(tppow_sprd1824) + bbs(tppow_sprd2430) +
    bbs(capepow_mean) + bbs(capepow_sprd) + bbs(capepow_min) + bbs(capepow_max) + #bbs(cape_frac) +
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
    bbs(tdiff500700_mean) + bbs(tdiff500700_min) + bbs(tdiff500700_max)
  
  
  # bamlss formula
  b.mu.formula <- robs ~ s(tppow_mean) + s(tppow_sprd) + s(tppow_min) + s(tppow_max) + #s(tp_frac) +  
    s(tppow_mean0612) + s(tppow_mean1218) + s(tppow_mean1824) + s(tppow_mean2430) + 
    s(tppow_sprd0612) + s(tppow_sprd1218) + s(tppow_sprd1824) + s(tppow_sprd2430) +
    s(capepow_mean) + s(capepow_sprd) + s(capepow_min) + s(capepow_max) + #s(cape_frac) +
    s(capepow_mean0612) + s(capepow_mean1218) + s(capepow_mean1224) + s(capepow_mean1230) +
    s(capepow_sprd0612) + s(capepow_sprd1218) + s(capepow_sprd1224) + s(capepow_sprd1230) +
    s(dswrf_mean_mean) + s(dswrf_mean_max) +  #s(dswrf_mean_min) +
    s(dswrf_sprd_mean) + s(dswrf_sprd_max) + #s(dswrf_sprd_min) +
    s(msl_mean_mean) + s(msl_mean_min) + s(msl_mean_max) + 
    s(msl_sprd_mean) + s(msl_sprd_min) + s(msl_sprd_max) +
    s(pwat_mean_mean) + s(pwat_mean_min) + s(pwat_mean_max) + 
    s(pwat_sprd_mean) + s(pwat_sprd_min) + s(pwat_sprd_max) +
    s(tmax_mean_mean) + s(tmax_mean_min) + s(tmax_mean_max) +
    s(tmax_sprd_mean) + s(tmax_sprd_min) + s(tmax_sprd_max) +
    s(tcolc_mean_mean) + s(tcolc_mean_min) + s(tcolc_mean_max) +
    s(tcolc_sprd_mean) + s(tcolc_sprd_min) + s(tcolc_sprd_max) +
    s(t500_mean_mean) + s(t500_mean_min) + s(t500_mean_max) +
    s(t700_mean_mean) + s(t700_mean_min) + s(t700_mean_max) +
    s(t850_mean_mean) + s(t850_mean_min) + s(t850_mean_max) +
    s(t500_sprd_mean) + s(t500_sprd_min) + s(t500_sprd_max) +
    s(t700_sprd_mean) + s(t700_sprd_min) + s(t700_sprd_max) +
    s(t850_sprd_mean) + s(t850_sprd_min) + s(t850_sprd_max) +
    s(tdiff500850_mean) + s(tdiff500850_min) + s(tdiff500850_max) +
    s(tdiff700850_mean) + s(tdiff700850_min) + s(tdiff700850_max) +
    s(tdiff500700_mean) + s(tdiff500700_min) + s(tdiff500700_max)
  
  b.sigma.formula <- ~ s(tppow_mean) + s(tppow_sprd) + s(tppow_min) + s(tppow_max) + #s(tp_frac) +  
    s(tppow_mean0612) + s(tppow_mean1218) + s(tppow_mean1824) + s(tppow_mean2430) + 
    s(tppow_sprd0612) + s(tppow_sprd1218) + s(tppow_sprd1824) + s(tppow_sprd2430) +
    s(capepow_mean) + s(capepow_sprd) + s(capepow_min) + s(capepow_max) + #s(cape_frac) +
    s(capepow_mean0612) + s(capepow_mean1218) + s(capepow_mean1224) + s(capepow_mean1230) +
    s(capepow_sprd0612) + s(capepow_sprd1218) + s(capepow_sprd1224) + s(capepow_sprd1230) +
    s(dswrf_mean_mean) + s(dswrf_mean_max) +  #s(dswrf_mean_min) +
    s(dswrf_sprd_mean) + s(dswrf_sprd_max) + #s(dswrf_sprd_min) +
    s(msl_mean_mean) + s(msl_mean_min) + s(msl_mean_max) + 
    s(msl_sprd_mean) + s(msl_sprd_min) + s(msl_sprd_max) +
    s(pwat_mean_mean) + s(pwat_mean_min) + s(pwat_mean_max) + 
    s(pwat_sprd_mean) + s(pwat_sprd_min) + s(pwat_sprd_max) +
    s(tmax_mean_mean) + s(tmax_mean_min) + s(tmax_mean_max) +
    s(tmax_sprd_mean) + s(tmax_sprd_min) + s(tmax_sprd_max) +
    s(tcolc_mean_mean) + s(tcolc_mean_min) + s(tcolc_mean_max) +
    s(tcolc_sprd_mean) + s(tcolc_sprd_min) + s(tcolc_sprd_max) +
    s(t500_mean_mean) + s(t500_mean_min) + s(t500_mean_max) +
    s(t700_mean_mean) + s(t700_mean_min) + s(t700_mean_max) +
    s(t850_mean_mean) + s(t850_mean_min) + s(t850_mean_max) +
    s(t500_sprd_mean) + s(t500_sprd_min) + s(t500_sprd_max) +
    s(t700_sprd_mean) + s(t700_sprd_min) + s(t700_sprd_max) +
    s(t850_sprd_mean) + s(t850_sprd_min) + s(t850_sprd_max) +
    s(tdiff500850_mean) + s(tdiff500850_min) + s(tdiff500850_max) +
    s(tdiff700850_mean) + s(tdiff700850_min) + s(tdiff700850_max) +
    s(tdiff500700_mean) + s(tdiff500700_min) + s(tdiff500700_max)
  }
  
  
  
  ## prepare rain data of all stations but only July and from 1985
  rain <- rain_all
  ## choose one month:
  # July
  rain <- rain[(((rain$year %% 4 == 0) & rain$day>=183) | (!(rain$year %% 4 == 0) & rain$day>=182)),]
  rain <- rain[(((rain$year %% 4 == 0) & rain$day<=213) | (!(rain$year %% 4 == 0) & rain$day<=212)),]
  
  # predictions start in 1985
  rain <- rain[rain$year>=85,]
  
  rain_all_July85 <- rain
  
  
  ########################################
  # loop over all stations (which provide sufficient data)
  
  rainres <- mclapply(1:length(complete_stations),
                      function(i){
                        
                        set.seed(seedconst)
                        
                        rain <- rain_all_July85
                        
                        stationname <- stations[complete_stations[i], "stationname"]
                        stationnr <- stations[complete_stations[i], "stationnr"]
                        
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
                        
                        
                        ########################################
                        # remove rows with missing values or NAs
                        # (as a result the years 93 and 94 are dropped, for year 92 only one observation is left,
                        #  in 9 other years 1-3 observations are dropped)
                        
                        #remove row 825
                        raindata <- raindata[-825,]
                        #remove rows 223-310
                        raindata <- raindata[-c(223:310),]
                        #remove rows 111 176 178 220 221 241 320 364 367
                        raindata <- raindata[-c(111, 176, 178, 220, 221, 241, 320, 364, 367),]
                        #remove rows 5 215 217 301 314
                        raindata <- raindata[-c(5, 215, 217, 301, 314),]
                        #remove rows 251 264 
                        raindata <- raindata[-c(251, 264),]
                        
                        # table(raindata[, "year"])
                        
                        ######################################################
                        # only keep variables with sufficient values
                        raindata <- raindata[, c("robs", "year",
                                                 "tppow_mean", "tppow_sprd", "tppow_min", "tppow_max", "tp_frac", 
                                                 "tppow_mean0612", "tppow_mean1218", "tppow_mean1824", "tppow_mean2430", 
                                                 "tppow_sprd0612", "tppow_sprd1218", "tppow_sprd1824", "tppow_sprd2430",
                                                 "capepow_mean", "capepow_sprd", "capepow_min", "capepow_max", "cape_frac",
                                                 "capepow_mean0612", "capepow_mean1218", "capepow_mean1224", "capepow_mean1230",
                                                 "capepow_sprd0612", "capepow_sprd1218", "capepow_sprd1224", "capepow_sprd1230",
                                                 "dswrf_mean_mean", "dswrf_mean_min", "dswrf_mean_max",
                                                 "dswrf_sprd_mean", "dswrf_sprd_min", "dswrf_sprd_max",
                                                 "msl_mean_mean", "msl_mean_min", "msl_mean_max",
                                                 "msl_sprd_mean", "msl_sprd_min", "msl_sprd_max",
                                                 "pwat_mean_mean", "pwat_mean_min", "pwat_mean_max",
                                                 "pwat_sprd_mean", "pwat_sprd_min", "pwat_sprd_max",
                                                 "tmax_mean_mean", "tmax_mean_min", "tmax_mean_max",
                                                 "tmax_sprd_mean", "tmax_sprd_min", "tmax_sprd_max",
                                                 "tcolc_mean_mean", "tcolc_mean_min", "tcolc_mean_max",
                                                 "tcolc_sprd_mean", "tcolc_sprd_min", "tcolc_sprd_max",
                                                 "t500_mean_mean", "t500_mean_min", "t500_mean_max",
                                                 "t700_mean_mean", "t700_mean_min", "t700_mean_max",
                                                 "t850_mean_mean", "t850_mean_min", "t850_mean_max",
                                                 "t500_sprd_mean", "t500_sprd_min", "t500_sprd_max",
                                                 "t700_sprd_mean", "t700_sprd_min", "t700_sprd_max",
                                                 "t850_sprd_mean", "t850_sprd_min", "t850_sprd_max",
                                                 "tdiff500850_mean", "tdiff500850_min", "tdiff500850_max",
                                                 "tdiff700850_mean", "tdiff700850_min", "tdiff700850_max",
                                                 "tdiff500700_mean", "tdiff500700_min", "tdiff500700_max")]
                        
                        ######################################################
                        #new vairable msl_diff
                        raindata$msl_diff <- raindata$msl_mean_max - raindata$msl_mean_min
                        
                        ## FIX ME: check years (85-112)
                        # table(raindata$year)
                        
                        # learning the models on 29 years and evaluating predictions on the 30th year
                        learndata <- raindata[raindata$year < 112,]
                        testdata <- raindata[raindata$year == 112,]
                        
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
                                        family = cens("NO", type = "left")))
                        if(inherits(g, "try-error")) {
                          g <- NA
                          g_error <- paste0(stationname,"_model")
                        } else g_error <- NA
                        
                        gb <- gamboostLSS(formula = list(mu = gb.mu.formula, sigma = gb.sigma.formula), data = g_learndata, 
                                          families = as.families(fname = cens("NO", type = "left")), method = "noncyclic",
                                          control = boost_control(mstop = 1000L))
                        if(gamboost_cvr){
                          grid <- c(seq(50,500, by = 25), seq(510, 1000, by = 10))
                          cvr <- cvrisk(gb, grid = grid)
                          mstop(gb) <- mstop(cvr)
                          cvr_opt <- mstop(cvr) 
                        } else cvr_opt <- NA
                        
                        #b <- bamlss(list(b.mu.formula, b.sigma.formula), family = "cnorm", 
                        #            data = learndata, sampler = FALSE, optimizer = boost, 
                        #            stop.criterion = "BIC", plot = FALSE) #, nu = 0.1)
                        
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
                        
                        # gamboostLSS
                        pgb <- predict(gb, newdata = testdata, parameter = list("mu","sigma"), type = "response")
                        gb_mu <- pgb[[1]]
                        gb_sigma <- pgb[[2]]
                        gb_exp <- pnorm(gb_mu/gb_sigma) * (gb_mu + gb_sigma * (dnorm(gb_mu/gb_sigma) / pnorm(gb_mu/gb_sigma)))
                        
                        
                        # bamlss
                        #pb <- predict(b, newdata = testdata, type = "parameter")
                        #b_mu <- pb$mu
                        #b_sigma <- pb$sigma
                        #b_exp <- pnorm(b_mu/b_sigma) * (b_mu + b_sigma * (dnorm(b_mu/b_sigma) / pnorm(b_mu/b_sigma)))
                        
                        
                        # gamlss
                        if(!(all(is.na(g)))){
                          g_mu <- try(predict(g, newdata = testdata, what = "mu", type = "response", data = g_learndata))
                          g_sigma <- try(predict(g, newdata = testdata, what = "sigma", type = "response", data = g_learndata))
                          if(inherits(g_mu, "try-error") | inherits(g_sigma, "try-error")) {
                            g_mu <- g_sigma <- g_exp <- NA
                            g_error <- paste0(stationname,"_pred")
                          } else g_exp <- pnorm(g_mu/g_sigma) * (g_mu + g_sigma * (dnorm(g_mu/g_sigma) / pnorm(g_mu/g_sigma)))
                        }
                        g_na <- any(c(all(is.na(g)), is.na(g_mu), is.na(g_sigma)))
                        
                        # EMOS
                        if(!(all(is.na(mi)))){
                          mi_mu <- try(predict(mi, type = "location", newdata = testdata))     # returns parameter on response scale
                          mi_sigma <- try(predict(mi, type = "scale", newdata = testdata))
                          if(inherits(mi_mu, "try-error") | inherits(mi_sigma, "try-error")) {
                            mi_mu <- mi_sigma <- mi_exp <- NA
                            mi_error <- paste0(stationname,"_pred")
                          } else mi_exp <- pnorm(mi_mu/mi_sigma) * (mi_mu + mi_sigma * (dnorm(mi_mu/mi_sigma) / pnorm(mi_mu/mi_sigma)))
                        }
                        mi_na <- any(c(all(is.na(mi)), is.na(mi_mu), is.na(mi_sigma)))
                        
                        if(!(all(is.na(ml)))){
                          ml_mu <- try(predict(ml, type = "location", newdata = testdata))     # returns parameter on response scale
                          ml_sigma <- try(predict(ml, type = "scale", newdata = testdata))
                          if(inherits(ml_mu, "try-error") | inherits(ml_sigma, "try-error")) {
                            ml_mu <- ml_sigma <- ml_exp <- NA
                            ml_error <- paste0(stationname,"_pred")
                          } else ml_exp <- pnorm(ml_mu/ml_sigma) * (ml_mu + ml_sigma * (dnorm(ml_mu/ml_sigma) / pnorm(ml_mu/ml_sigma)))
                        }
                        ml_na <- any(c(all(is.na(ml)), is.na(ml_mu), is.na(ml_sigma)))
                        
                        if(!(all(is.na(mq)))){
                          mq_mu <- try(predict(mq, type = "location", newdata = testdata))     # returns parameter on response scale
                          mq_sigma <- try(predict(mq, type = "scale", newdata = testdata))
                          if(inherits(mq_mu, "try-error") | inherits(mq_sigma, "try-error")) {
                            mq_mu <- mq_sigma <- mq_exp <- NA
                            mq_error <- paste0(stationname,"_pred")
                          } else mq_exp <- pnorm(mq_mu/mq_sigma) * (mq_mu + mq_sigma * (dnorm(mq_mu/mq_sigma) / pnorm(mq_mu/mq_sigma)))
                        }
                        mq_na <- any(c(all(is.na(mq)), is.na(mq_mu), is.na(mq_sigma)))
                        
                        
                        # CPRS
                        crps_dt <- sum(crps_cnorm(testdata$robs, location = dt_mu, scale = dt_sigma, lower = 0, upper = Inf))
                        crps_df <- sum(crps_cnorm(testdata$robs, location = df_mu, scale = df_sigma, lower = 0, upper = Inf))
                        crps_g <- if(!g_na) sum(crps_cnorm(testdata$robs, location = g_mu, scale = g_sigma, lower = 0, upper = Inf)) else NA
                        #crps_b <- sum(crps_cnorm(testdata$robs, location = b_mu, scale = b_sigma, lower = 0, upper = Inf))
                        crps_gb <- sum(crps_cnorm(testdata$robs, location = gb_mu, scale = gb_sigma, lower = 0, upper = Inf))
                        crps_mi <- if(!mi_na) sum(crps_cnorm(testdata$robs, location = mi_mu, scale = mi_sigma, lower = 0, upper = Inf)) else NA
                        crps_ml <- if(!ml_na) sum(crps_cnorm(testdata$robs, location = ml_mu, scale = ml_sigma, lower = 0, upper = Inf)) else NA
                        crps_mq <- if(!mq_na) sum(crps_cnorm(testdata$robs, location = mq_mu, scale = mq_sigma, lower = 0, upper = Inf)) else NA
                        
                        # RMSE
                        rmse_dt <- sqrt(mean((dt_exp - testdata[,"robs"])^2))
                        rmse_df <- sqrt(mean((df_exp - testdata[,"robs"])^2))
                        rmse_g <- if(!g_na) sqrt(mean((g_exp - testdata[,"robs"])^2)) else NA
                        rmse_gb <- sqrt(mean((gb_exp - testdata[,"robs"])^2))
                        #rmse_b <- sqrt(mean((b_exp - testdata[,"robs"])^2))
                        rmse_mi <- if(!mi_na) sqrt(mean((mi_exp - testdata[,"robs"])^2)) else NA
                        rmse_ml <- if(!ml_na) sqrt(mean((ml_exp - testdata[,"robs"])^2)) else NA
                        rmse_mq <- if(!mq_na) sqrt(mean((mq_exp - testdata[,"robs"])^2)) else NA
                        
                        # loglikelihood
                        dtll <- dfll <- gll <-  gbll <- mill <- mlll <- mqll <- 0
                        # bll <- 0
                        for(j in 1:(nrow(testdata))){
                          
                          eta_dt <- as.numeric(dist_list_cens_normal$linkfun(cbind(dt_mu, dt_sigma)[j,]))
                          eta_df <- as.numeric(dist_list_cens_normal$linkfun(cbind(df_mu, df_sigma)[j,]))
                          eta_g <- if(!g_na) as.numeric(dist_list_cens_normal$linkfun(cbind(g_mu, g_sigma)[j,])) else NA
                          eta_gb <- as.numeric(dist_list_cens_normal$linkfun(cbind(gb_mu, gb_sigma)[j,]))
                          #eta_b <- as.numeric(dist_list_cens_normal$linkfun(cbind(b_mu, b_sigma)[j,]))
                          eta_mi <- if(!mi_na) as.numeric(dist_list_cens_normal$linkfun(cbind(mi_mu, mi_sigma)[j,])) else NA
                          eta_ml <- if(!ml_na) as.numeric(dist_list_cens_normal$linkfun(cbind(ml_mu, ml_sigma)[j,])) else NA
                          eta_mq <- if(!mq_na) as.numeric(dist_list_cens_normal$linkfun(cbind(mq_mu, mq_sigma)[j,])) else NA
                          
                          dtll_j <- dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_dt, log=TRUE)
                          dfll_j <- dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_df, log=TRUE)
                          gll_j <- if(!g_na) dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_g, log=TRUE) else NA
                          gbll_j <- dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_gb, log=TRUE)
                          #bll_j <- dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_b, log=TRUE)
                          mill_j <- if(!mi_na) dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_mi, log=TRUE) else NA
                          mlll_j <- if(!mi_na) dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_ml, log=TRUE) else NA
                          mqll_j <- if(!mi_na) dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_mq, log=TRUE) else NA
                          
                          dtll <- if(is.na(dtll_j)) {
                            print(eta_dt, testdata[j,"robs"]) 
                            dtll + (-5)
                          } else {dtll + dtll_j}
                          
                          dfll <- if(is.na(dfll_j)) {
                            print(eta_df, testdata[j,"robs"]) 
                            dfll + (-5)
                          } else {dfll + dfll_j}    ## FIX ME: NAs from distforest
                          
                          gll <- if(!g_na) {
                            if(is.na(gll_j)) {
                              print(eta_g, testdata[j,"robs"]) 
                              gll + (-5)
                            } else {gll + gll_j}
                          } else NA
                          
                          gbll <- if(is.na(gbll_j)) {
                            print(eta_gb, testdata[j,"robs"]) 
                            gbll + (-5)
                          } else {gbll + gbll_j}
                          
                          #bll <- if(is.na(bll_j)) {
                          #  print(eta_b, testdata[j,"robs"]) 
                          #  bll + (-5)
                          #} else {bll + bll_j}   
                          
                          mill <- if(!mi_na) {
                            if(is.na(mill_j)) {
                              print(eta_mi, testdata[j,"robs"]) 
                              mill + (-5)
                            } else {mill + mill_j} 
                          } else NA      
                          
                          mlll <- if(!ml_na) {
                            if(is.na(mlll_j)) {
                              print(eta_ml, testdata[j,"robs"]) 
                              mlll + (-5)
                            } else {mlll + mlll_j} 
                          } else NA
                          
                          mqll <- if(!mq_na) {
                            if(is.na(mqll_j)) {
                              print(eta_mq, testdata[j,"robs"]) 
                              mqll + (-5)
                            } else {mqll + mqll_j} 
                          } else NA
                        }
                        
                        ll <- c(dtll, dfll, gll, gbll, mill, mlll, mqll) 
                        rmse <- c(rmse_dt, rmse_df, rmse_g, rmse_gb, rmse_mi, rmse_ml, rmse_mq)
                        crps <- c(crps_dt, crps_df, crps_g, crps_gb, crps_mi, crps_ml, crps_mq)
                        
                        res <- list()
                        res$results <- rbind(ll, rmse, crps)
                        colnames(res$results) <-
                          c("disttree", "distforest", "gamlss", "gamboostLSS", "EMOS id", "EMOS log", "EMOS quad")
                        rownames(res$results) <- c("ll", "rmse", "crps")
                        
                        res$cvr_opt <- cvr_opt
                        res$g_error <- g_error
                        res$mi_error <- mi_error
                        res$ml_error <- ml_error
                        res$mq_error <- mq_error
                        
                        return(res)
                      },
                      mc.cores = detectCores() - 1
  )
                      
  rainres$call <- cl
  rainres$complete_stations <- complete_stations
  
  return(rainres)
}                    
  



setwd("~/svn/partykit/pkg/disttree/inst/draft")
#source("rain_pred.R")
library("gamlss.cens")
gen.cens("NO", type = "left")

#save(res, file = "~/svn/partykit/pkg/disttree/inst/draft/rain_pred.rda")
#save(res, file = "~/disttree/inst/draft/rain_pred.rda")
res <- rain_pred(seedconst = 7, ntree = 200,
                 tree_minsplit = 50, tree_minbucket = 25, tree_mincrit = 0.95,
                 forest_minsplit = 50, forest_minbucket = 25, forest_mincrit = 0,
                 forest_mtry = 27,
                 type.tree = "ctree",
                 gamboost_cvr = FALSE)

if(FALSE){
  save(res, file = "rain_pred.rda")
  
  ll <- res[[1]]$results["ll",]
  for(i in 2:length(res)) ll <- rbind(ll, res[[i]]$results["ll",])
  
  rmse <- res[[1]]$results["rmse",]
  for(i in 2:length(res)) rmse <- rbind(rmse, res[[i]]$results["rmse",])
  
  crps <- res[[1]]$results["crps",]
  for(i in 2:length(res))  crps <- rbind(crps, res[[i]]$results["crps",])
  
  colnames(ll) <- colnames(rmse) <- colnames(crps) <- colnames(res[[1]]$results)
  
  boxplot(rmse)
  boxplot(crps)
  boxplot(ll)
  
}


### FIX ME:
# in disttree (and distforest): for nodes with all equal observations (all 0), sigma is set to 0.0001
# this leads to extremely low values in the loglikelihood if the observed value is not 0
# replace? warning?



