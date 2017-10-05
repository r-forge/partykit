## choose one station: 
# 9101121 Lech (EHYD)
# 9103135 Fieberbrunn (EHYD)
# 9102954 Nassereith (EHYD)
# 9102756 Steinberg am Rofan (EHYD)  
# 9101113 Zuers (EHYD)

# 9102913 Walchsee (EHYD)
# 9102905 Koessen (EHYD)
# 9103309 Axams (EHYD)
# 9103176 Niederndorferberg (EHYD)
# 9119479 Ammerwald (EHYD)


rain_cross <- function(stationname, seedconst = 7, ntree = 100,
   tree_minsplit = 50, tree_minbucket = 20, tree_mincrit = 0.95,
   forest_minsplit = 50, forest_minbucket = 20, forest_mincrit = 0,
   forest_mtry = 27,   # if frac == TRUE: nvar/3 = 30
   frac = FALSE,
   type.tree = "ctree",
   gamboost_cvr = FALSE)
{
  
  if(!is.character(stationname)) stop("argument 'stationname' has to be a character")
  if(!("(EHYD)" %in% strsplit(stationname, split = " ")[[1]])) stationname <- paste(stationname, "(EHYD)", sep = " ")
  
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
  #cbind(ehyd.statlist$station, ehyd.statlist$stationname, ehyd.statlist$height)
  #stations <- cbind(ehyd.statlist$station, ehyd.statlist$stationname)
  stationnr <- ehyd.statlist$station[match(stationname, ehyd.statlist$stationname)]
  
  
  
  
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
  #head(rain)
  
  ## choose one month:
  # July
  rain <- rain[(((rain$year %% 4 == 0) & rain$day>=183) | (!(rain$year %% 4 == 0) & rain$day>=182)),]
  rain <- rain[(((rain$year %% 4 == 0) & rain$day<=213) | (!(rain$year %% 4 == 0) & rain$day<=212)),]
  
  # predictions start in 1985
  rain <- rain[rain$year>=85,]
  
  # observations end in 2012
  prediction <- prediction[prediction$year<=112,]
  
  
  # combine data frames
  
  # predictions with init 07-01 are for day 07-02 
  # observations are the sum of precipitation of the last 24h
  # match prediction with init 07-01 with observation on 07-02 
  # (pred starts at 00UTC and predicts from 06UTC until 30UTC = 06UTC of next day, 
  # observations are meassured from 06UTC of previous day to 06UTC of current day)
  
  raindata <- cbind(rain$date, rain$obs, rain$robs, prediction)
  colnames(raindata)[c(1:3)] <- c("date", "obs", "robs")
  
  
  
  
  
  
  
  
  ########################################
  # remove rows with missing values or NAs
  # (as a result the years 93 and 94 are dropped, for year 92 only one observation is left,
  #  in 9 other years 1-3 observations are dropped)
  
  #remove row 825
  raindata <- raindata[-825,]
  #remove rows 223-310
  #raindata <- raindata[-c(223:310),]
  #remove rows 111 176 178 220 221 241 320 364 367
  #raindata <- raindata[-c(111, 176, 178, 220, 221, 241, 320, 364, 367),]
  #remove rows 5 215 217 301 314
  #raindata <- raindata[-c(5, 215, 217, 301, 314),]
  #remove rows 251 264 
  #raindata <- raindata[-c(251, 264),]
  
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
  
  ######################################################
  #new vairable msl_diff
  raindata$msl_diff <- raindata$msl_mean_max - raindata$msl_mean_min
  
  
  #dim(raindata)
  #any(is.na(raindata))
  

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
  
  
  
  
  
  ############################################################
  # cross validation (10x10)
  
  
  nrep_cross <- 10
  #grid to find optimal mstop for method="noncyclic" in gamboostLSS
  grid <- seq(50, 1000, by = 25)
  #grid <- c(seq(100,400, by = 100), seq(420, 1000, by = 20))
  #grid <- make.grid(max = 500, min = 10, length.out = 10)
  
  rainres <- mclapply(1:nrep_cross,
                      function(k){
                        
                        set.seed(seedconst * k)
                        
                        crps_dt <- crps_df <- crps_g <- crps_gb <- crps_mi <- crps_ml <- crps_mq <- numeric(length = 10L)
                        rmse_dt <- rmse_df <- rmse_g <- rmse_gb <- rmse_mi <- rmse_ml <- rmse_mq <- numeric(length = 10L)
                        ll_dt <- ll_df <- ll_g <- ll_gb <- ll_mi <- ll_ml <- ll_mq <- numeric(length = 10L)
                        cvr_opt <- numeric(length = 10L)
                        g_error_seed <- numeric(length = 10L)
                        mi_error_seed <- numeric(length = 10L)
                        ml_error_seed <- numeric(length = 10L)
                        mq_error_seed <- numeric(length = 10L)
                        evaltime_user <- matrix(0, ncol = 8, nrow = 10L)
                        evaltime_system <- matrix(0, ncol = 8, nrow = 10L)
                        evaltime_elapsed <- matrix(0, ncol = 8, nrow = 10L)
                        
                        # randomly split data in 10 roughly equal parts
                        # (here: 867 observations -> 7 groups of 87 and 3 groups of 86 obs)
                        id <- c(1:NROW(raindata))
                        testid <- list()
                        for(i in 1:7){
                          testid[[i]] <- sample(id, 87, replace = FALSE)
                          id <- id[!(id %in% testid[[i]])]
                        }
                        for(i in 8:10){
                          testid[[i]] <- sample(id, 86, replace = FALSE)
                          id <- id[!(id %in% testid[[i]])]
                        }
                        
                        
                        
                        for(i in 1:10) {
                          
                          testdata <- raindata[testid[[i]], ]
                          learndata <- raindata[-testid[[i]], ]
                          
                          if(type.tree == "mob"){
                            dt_time <- system.time(dt <- disttree(dt.formula, 
                                                                  data = learndata, family = dist_list_cens_normal, 
                                                                  censtype = "left", censpoint = 0, type.tree = "mob", 
                                                                  control = mob_control(restart = FALSE, numsplit = "center", 
                                                                                        alpha = 1-tree_mincrit, minsplit = tree_minsplit,
                                                                                        minbucket = tree_minbucket)))
                          }
                          if(type.tree == "ctree"){
                            dt_time <- system.time(dt <- disttree(dt.formula, 
                                                                  data = learndata, family = dist_list_cens_normal, 
                                                                  censtype = "left", censpoint = 0, type.tree = "ctree", 
                                                                  control = ctree_control(teststat = "quad", testtype = "Bonferroni", intersplit = TRUE,
                                                                                          mincriterion = tree_mincrit, minsplit = tree_minsplit,
                                                                                          minbucket = tree_minbucket)))
                          }
                          
                          
                          if(type.tree == "mob"){
                            df_time <- system.time(df <- distforest(df.formula, 
                                                                    data = learndata, family = dist_list_cens_normal, type.tree = "mob", 
                                                                    ntree = ntree, censtype = "left", censpoint = 0,
                                                                    control = mob_control(restart = FALSE, numsplit = "center", 
                                                                                          alpha = 1-forest_mincrit, minsplit = forest_minsplit,
                                                                                          minbucket = forest_minbucket), mtry = forest_mtry))
                          }
                          if(type.tree == "ctree"){
                            df_time <- system.time(df <- distforest(df.formula, 
                                                                    data = learndata, family = dist_list_cens_normal, type.tree = "ctree", 
                                                                    ntree = ntree, censtype = "left", censpoint = 0, #fitted.OOB = FALSE,
                                                                    control = ctree_control(teststat = "quad", testtype = "Univ", intersplit = TRUE,
                                                                                            mincriterion = forest_mincrit, minsplit = forest_minsplit,
                                                                                            minbucket = forest_minbucket), mtry = forest_mtry))
                          }
                          
                          
                          
                          
                          g_learndata <- learndata
                          g_learndata$robs <- Surv(g_learndata$robs, g_learndata$robs>0, type="left")
                          
                          g_time <- system.time(g <- try(gamlss(formula = g.mu.formula, sigma.formula = g.sigma.formula, data = g_learndata, 
                                                                family = cens("NO", type = "left"),
                                                                control = gamlss.control(n.cyc = 100),
                                                                i.control = glim.control(cyc = 100, bf.cyc = 100))))
                          if(inherits(g, "try-error")) {
                            g_time <- NA
                            g <- NA
                            g_error_seed[i] <- seedconst * k
                          }
                          
                          
                          #gb <- gb
                          gb_time <- system.time(gb <- gamboostLSS(formula = list(mu = gb.mu.formula, sigma = gb.sigma.formula), data = g_learndata,
                                                                   families = as.families(fname = cens("NO", type = "left")), method = "noncyclic",
                                                                   control = boost_control(mstop = 1000L)))
                          if(gamboost_cvr){
                            gb_cvr_time <- system.time(cvr <- cvrisk(gb, grid = grid))
                            mstop(gb) <- mstop(cvr)
                            cvr_opt[i] <- mstop(cvr) 
                          } else cvr_opt[i] <- gb_cvr_time <- NA
                          

                          mi_time <- system.time(mi <- try(crch(formula = robs ~ tppow_mean | tppow_sprd, 
                                                                data = learndata, dist = "gaussian", left = 0, link.scale = "identity")))
                          if(inherits(mi, "try-error")) {
                            mi_time <- NA
                            mi <- NA
                            mi_error_seed[i] <- seedconst * k
                          }
                          
                          ml_time <- system.time(ml <- try(crch(formula = robs ~ tppow_mean | tppow_sprd, 
                                                                data = learndata, dist = "gaussian", left = 0, link.scale = "identity")))
                          if(inherits(ml, "try-error")) {
                            ml_time <- NA
                            ml <- NA
                            ml_error_seed[i] <- seedconst * k
                          }
                          
                          mq_time <- system.time(mq <- try(crch(formula = robs ~ tppow_mean | tppow_sprd, 
                                                                data = learndata, dist = "gaussian", left = 0, link.scale = "identity")))
                          if(inherits(mq, "try-error")) {
                            mq_time <- NA
                            mq <- NA
                            mq_error_seed[i] <- seedconst * k
                          }
                          
                          
                          ## get predicted parameter
                          
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
                              g_error_seed[i] <- seedconst * k
                            } else g_exp <- pnorm(g_mu/g_sigma) * (g_mu + g_sigma * (dnorm(g_mu/g_sigma) / pnorm(g_mu/g_sigma)))
                          } else g_mu <- g_sigma <- g_exp <- NA
                          g_na <- any(c(all(is.na(g)), all(is.na(g_mu)), all(is.na(g_sigma))))
                          
                          # gamboostLSS
                          pgb <- predict(gb, newdata = testdata, parameter = list("mu","sigma"), type = "response")
                          gb_mu <- pgb[[1]]
                          gb_sigma <- pgb[[2]]
                          gb_exp <- pnorm(gb_mu/gb_sigma) * (gb_mu + gb_sigma * (dnorm(gb_mu/gb_sigma) / pnorm(gb_mu/gb_sigma)))
                          
                          # EMOS
                          if(!(all(is.na(mi)))){
                            mi_mu <- try(predict(mi, type = "location", newdata = testdata))     # returns parameter on response scale
                            mi_sigma <- try(predict(mi, type = "scale", newdata = testdata))
                            if(inherits(mi_mu, "try-error") | inherits(mi_sigma, "try-error")) {
                              mi_mu <- mi_sigma <- mi_exp <- NA
                              mi_error_seed[i] <- seedconst * k
                            } else mi_exp <- pnorm(mi_mu/mi_sigma) * (mi_mu + mi_sigma * (dnorm(mi_mu/mi_sigma) / pnorm(mi_mu/mi_sigma)))
                          } else mi_mu <- mi_sigma <- mi_exp <- NA
                          mi_na <- any(c(all(is.na(mi)), all(is.na(mi_mu)), all(is.na(mi_sigma))))
                          
                          if(!(all(is.na(ml)))){
                            ml_mu <- try(predict(ml, type = "location", newdata = testdata))     # returns parameter on response scale
                            ml_sigma <- try(predict(ml, type = "scale", newdata = testdata))
                            if(inherits(ml_mu, "try-error") | inherits(ml_sigma, "try-error")) {
                              ml_mu <- ml_sigma <- ml_exp <- NA
                              ml_error_seed[i] <- seedconst * k
                            } else ml_exp <- pnorm(ml_mu/ml_sigma) * (ml_mu + ml_sigma * (dnorm(ml_mu/ml_sigma) / pnorm(ml_mu/ml_sigma)))
                          } else ml_mu <- ml_sigma <- ml_exp <- NA
                          ml_na <- any(c(all(is.na(ml)), all(is.na(ml_mu)), all(is.na(ml_sigma))))
                          
                          if(!(all(is.na(mq)))){
                            mq_mu <- try(predict(mq, type = "location", newdata = testdata))     # returns parameter on response scale
                            mq_sigma <- try(predict(mq, type = "scale", newdata = testdata))
                            if(inherits(mq_mu, "try-error") | inherits(mq_sigma, "try-error")) {
                              mq_mu <- mq_sigma <- mq_exp <- NA
                              mq_error_seed[i] <- seedconst * k
                            } else mq_exp <- pnorm(mq_mu/mq_sigma) * (mq_mu + mq_sigma * (dnorm(mq_mu/mq_sigma) / pnorm(mq_mu/mq_sigma)))
                          } else mq_mu <- mq_sigma <- mq_exp <- NA
                          mq_na <- any(c(all(is.na(mq)), all(is.na(mq_mu)), all(is.na(mq_sigma))))
                          
                          
                          # CPRS
                          crps_dt[i] <- mean(crps_cnorm(testdata$robs, location = dt_mu, scale = dt_sigma, lower = 0, upper = Inf), na.rm = TRUE)
                          crps_df[i] <- mean(crps_cnorm(testdata$robs, location = df_mu, scale = df_sigma, lower = 0, upper = Inf), na.rm = TRUE)
                          crps_g[i] <- if(!g_na) mean(crps_cnorm(testdata$robs, location = g_mu, scale = g_sigma, lower = 0, upper = Inf), na.rm = TRUE) else NA
                          crps_gb[i] <- mean(crps_cnorm(testdata$robs, location = gb_mu, scale = gb_sigma, lower = 0, upper = Inf), na.rm = TRUE)
                          crps_mi[i] <- if(!mi_na) mean(crps_cnorm(testdata$robs, location = mi_mu, scale = mi_sigma, lower = 0, upper = Inf), na.rm = TRUE) else NA
                          crps_ml[i] <- if(!ml_na) mean(crps_cnorm(testdata$robs, location = ml_mu, scale = ml_sigma, lower = 0, upper = Inf), na.rm = TRUE) else NA
                          crps_mq[i] <- if(!mq_na) mean(crps_cnorm(testdata$robs, location = mq_mu, scale = mq_sigma, lower = 0, upper = Inf), na.rm = TRUE) else NA
                          
                          # RMSE
                          rmse_dt[i] <- sqrt(mean((dt_exp - testdata[,"robs"])^2, na.rm = TRUE))
                          rmse_df[i] <- sqrt(mean((df_exp - testdata[,"robs"])^2, na.rm = TRUE))
                          rmse_g[i] <- if(!g_na) sqrt(mean((g_exp - testdata[,"robs"])^2, na.rm = TRUE)) else NA
                          rmse_gb[i] <- sqrt(mean((gb_exp - testdata[,"robs"])^2, na.rm = TRUE))
                          rmse_mi[i] <- if(!mi_na) sqrt(mean((mi_exp - testdata[,"robs"])^2, na.rm = TRUE)) else NA
                          rmse_ml[i] <- if(!ml_na) sqrt(mean((ml_exp - testdata[,"robs"])^2, na.rm = TRUE)) else NA
                          rmse_mq[i] <- if(!mq_na) sqrt(mean((mq_exp - testdata[,"robs"])^2, na.rm = TRUE)) else NA
                          
                          # loglikelihood
                          dtll <- dfll <- gll <-  gbll <- mill <- mlll <- mqll <- numeric(length = NROW(testdata))
                          for(j in 1:(NROW(testdata))){
                            
                            eta_dt <- as.numeric(dist_list_cens_normal$linkfun(cbind(dt_mu, dt_sigma)[j,]))
                            eta_df <- as.numeric(dist_list_cens_normal$linkfun(cbind(df_mu, df_sigma)[j,]))
                            eta_g <- if(!g_na) as.numeric(dist_list_cens_normal$linkfun(cbind(g_mu, g_sigma)[j,])) else NA
                            eta_gb <- as.numeric(dist_list_cens_normal$linkfun(cbind(gb_mu, gb_sigma)[j,]))
                            #eta_b <- as.numeric(dist_list_cens_normal$linkfun(cbind(b_mu, b_sigma)[j,]))
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
                          
                          ll_dt[i] <- mean(dtll, na.rm = TRUE)
                          ll_df[i] <- mean(dfll, na.rm = TRUE)
                          ll_g[i] <- if(!g_na) mean(gll, na.rm = TRUE) else NA
                          ll_gb[i] <- mean(gbll, na.rm = TRUE)
                          ll_mi[i] <- if(!mi_na) mean(mill, na.rm = TRUE) else NA      
                          ll_ml[i] <- if(!ml_na) mean(mlll, na.rm = TRUE) else NA 
                          ll_mq[i] <- if(!mq_na) mean(mqll, na.rm = TRUE) else NA 
                          
                          evaltime_user[i,] <- c(dt_time[1], df_time[1], g_time[1], gb_time[1], gb_cvr_time[1], mi_time[1], ml_time[1], mq_time[1])
                          evaltime_system[i,] <- c(dt_time[2], df_time[2], g_time[2], gb_time[2], gb_cvr_time[2], mi_time[2], ml_time[2], mq_time[2])
                          evaltime_elapsed[i,] <- c(dt_time[3], df_time[3], g_time[3], gb_time[3], gb_cvr_time[3], mi_time[3], ml_time[3], mq_time[3])
                          
                        }
                        
                        res <- list()
                        
                        res$ll <- cbind(ll_dt, ll_df, ll_g, ll_gb, ll_mi, ll_ml, ll_mq)
                        res$rmse <- cbind(rmse_dt, rmse_df, rmse_g, rmse_gb, rmse_mi, rmse_ml, rmse_mq)
                        res$crps <- cbind(crps_dt, crps_df, crps_g, crps_gb, crps_mi, crps_ml, crps_mq)
                        colnames(res$rmse) <- colnames(res$ll) <- colnames(res$crps) <- 
                          c("disttree", "distforest", "gamlss", "gamboostLSS", "EMOS id", "EMOS log", "EMOS quad")
                        
                        res$cvr_opt <- cvr_opt
                        res$evaltime_user <- evaltime_user
                        res$evaltime_system <- evaltime_system
                        res$evaltime_elapsed <- evaltime_elapsed
                        colnames(res$evaltime_user) <- colnames(res$evaltime_system) <- colnames(res$evaltime_elapsed) <- 
                          c("disttree", "distforest", "gamlss", "gamboostLSS", "gamboostLSS_cvrisk", "EMOS id", "EMOS log", "EMOS quad")
                        res$g_error_seed <- g_error_seed
                        res$mi_error_seed <- mi_error_seed
                        res$ml_error_seed <- ml_error_seed
                        res$mq_error_seed <- mq_error_seed
                        
                        return(res)
                      },
                      mc.cores = detectCores() - 1
  )
  
  rainres$call <- cl
  
  return(rainres)
}




if(FALSE){
  library("gamlss.cens")
  gen.cens("NO", type = "left")
  setwd("~/svn/partykit/pkg/disttree/inst/draft/")
  rainres <- rain_cross(stationname = "Axams", 
                        seedconst = 7, ntree = 100,
                        forest_mtry = 27,
                        gamboost_cvr = TRUE,
                        frac = FALSE)
  #save(rainres, file = paste0("~/svn/partykit/pkg/disttree/inst/draft/rain_", rainres$call$stationname, rainres$call$seedconst, "frac_", rainres$call$frac, ".rda"))
  #save(rainres, file = paste0("~/disttree/inst/draft/rain_", rainres$call$stationname, rainres$call$seedconst, "frac_", rainres$call$frac, ".rda"))
  save(rainres, file = paste0("rain_", rainres$call$stationname, rainres$call$seedconst, "frac_", rainres$call$frac, ".rda"))
}





  
### FIX ME:
# in disttree (and distforest): for nodes with all equal observations (all 0), sigma is set to 0.0001
# this leads to extremely low values in the loglikelihood if the observed value is not 0
# replace? warning?

### TO DO:
# variable importance

if(FALSE){
  rain_rmse <- matrix(0, ncol = ncol(rainres[[1]]$rmse), length(rainres)-1)
  rain_ll <- matrix(0, ncol = ncol(rainres[[1]]$ll), length(rainres)-1)
  rain_crps <- matrix(0, ncol = ncol(rainres[[1]]$crps), length(rainres)-1)
  
  for(i in 1:(length(rainres)-1)){
    rain_rmse[i,] <- colMeans(rainres[[i]]$rmse, na.rm = TRUE)
    rain_ll[i,] <- colMeans(rainres[[i]]$ll, na.rm = TRUE)
    rain_crps[i,] <- colMeans(rainres[[i]]$crps, na.rm = TRUE)
  }
  
  colnames(rain_rmse) <- colnames(rain_ll) <- colnames(rain_crps) <-
    colnames(rainres[[1]]$rmse)
  
  boxplot(rain_rmse)
  boxplot(rain_ll)
  boxplot(rain_crps)
  
  unlist(lapply(1:10, function(i)rainres[[i]]$g_error_seed))
  unlist(lapply(1:10, function(i)rainres[[i]]$mi_error_seed))
  unlist(lapply(1:10, function(i)rainres[[i]]$ml_error_seed))
  unlist(lapply(1:10, function(i)rainres[[i]]$mq_error_seed))
  unlist(lapply(1:10, function(i)rainres[[i]]$cvr_opt))
  
  # scoring rules
  
  # skills score mit Referenz-Modell:
  # rain_crps / rain_crps[,7]
  # boxplot(rain_crps / rain_crps[,7])
  # boxplot(rain_crps[,-1] / rain_crps[,7])
  # boxplot(1 - rain_crps / rain_crps[,7])
  # abline(h = 0, col = 2, lwd = 2)
}




