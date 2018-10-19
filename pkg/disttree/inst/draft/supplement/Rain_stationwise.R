#######################################################
### Probabilistic Forecasting on Precipitation Data ###
#######################################################

## Replication material for Supplement 2 (Stationwise Evaluation) of the paper 
## Distributional Regression Forests for Probabilistic Precipitation Forecasting in Complex Terrain (2018)
## by Lisa Schlosser and Torsten Hothorn and Reto Stauffer and Achim Zeileis
## URL: http://arxiv.org/abs/1804.02921

## This demo includes the application on selected observation stations
## Full replication for station Axams can be obtained with
## demo("RainAxams", package = "disttree")
## Full replication of all other results presented in the paper can be obtained with
## demo("RainTyrol", package = "disttree")

## Computation time: approximately .. minutes (on our machines, using .. kernels)


library("disttree")

#####
# further packages
library("gamlss")
library("gamlss.dist")
library("gamboostLSS")
library("crch")
library("scoringRules")
library("RainTyrol")

# if gamlss.cens family object should be used as family
library("gamlss.cens")
gen.cens(NO, type = "left")

assign("NO",  gamlss.dist::NO,  pos = ".GlobalEnv")
assign("dNO", gamlss.dist::dNO, pos = ".GlobalEnv")
assign("pNO", gamlss.dist::pNO, pos = ".GlobalEnv")
assign("qNO", gamlss.dist::qNO, pos = ".GlobalEnv")
assign("rNO", gamlss.dist::rNO, pos = ".GlobalEnv")
gamlss.cens::gen.cens(NO, type = "left")
assign("NOlc",  NOlc,  pos = ".GlobalEnv")
assign("dNOlc", dNOlc, pos = ".GlobalEnv")
assign("pNOlc", pNOlc, pos = ".GlobalEnv")
assign("qNOlc", qNOlc, pos = ".GlobalEnv")


# set function for parallelization
applyfun <- function(X, FUN, ...) parallel::mclapply(X, FUN, ..., mc.cores = pmax(1, parallel::detectCores() - 1))


# HCL palette
pal <- hcl(c(10, 128, 260, 290, 50), 100, 50)


#####

# formula 
{
  
  # tree and forest formula
  dt.formula <- df.formula <- 
    robs ~ tppow_mean + tppow_sprd + tppow_min + tppow_max + 
    tppow_mean0612 + tppow_mean1218 + tppow_mean1824 + tppow_mean2430 + 
    tppow_sprd0612 + tppow_sprd1218 + tppow_sprd1824 + tppow_sprd2430 + 
    capepow_mean + capepow_sprd + capepow_min + capepow_max + 
    capepow_mean0612 + capepow_mean1218 + capepow_mean1224 + capepow_mean1230 +
    capepow_sprd0612 + capepow_sprd1218 + capepow_sprd1224 + capepow_sprd1230 +
    dswrf_mean_mean + dswrf_mean_max + 
    dswrf_sprd_mean + dswrf_sprd_max +
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
  
  # formula for prespecified GAM
  g.mu.formula <- robs ~ pb(tppow_mean) + 
    pb(tppow_mean1218 * capepow_mean1218) + 
    pb(tppow_max) + 
    pb(dswrf_mean_mean) +
    pb(tcolc_mean_mean) + 
    pb(msl_diff) + 
    pb(pwat_mean_mean) + 
    pb(tdiff500850_mean)  
  
  g.sigma.formula <- ~ pb(tppow_sprd) + 
    pb(tppow_sprd1218 * capepow_mean1218) + 
    pb(dswrf_sprd_mean) +
    pb(tcolc_sprd_mean) + 
    pb(tdiff500850_mean)  
  
  # formula for boosted GAM
  gb.mu.formula <- gb.sigma.formula <- 
    robs ~ bbs(tppow_mean) + bbs(tppow_sprd) + bbs(tppow_min) + bbs(tppow_max) + 
    bbs(tppow_mean0612) + bbs(tppow_mean1218) + bbs(tppow_mean1824) + bbs(tppow_mean2430) + 
    bbs(tppow_sprd0612) + bbs(tppow_sprd1218) + bbs(tppow_sprd1824) + bbs(tppow_sprd2430) +
    bbs(capepow_mean) + bbs(capepow_sprd) + bbs(capepow_min) + bbs(capepow_max) + 
    bbs(capepow_mean0612) + bbs(capepow_mean1218) + bbs(capepow_mean1224) + bbs(capepow_mean1230) +
    bbs(capepow_sprd0612) + bbs(capepow_sprd1218) + bbs(capepow_sprd1224) + bbs(capepow_sprd1230) +
    bbs(dswrf_mean_mean) + bbs(dswrf_mean_max) +
    bbs(dswrf_sprd_mean) + bbs(dswrf_sprd_max) +
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


data("RainTyrol")

stationeval <- function(station) {
  
  #####
  # get observations and covariates for selected station
  RainData <- RainTyrol[RainTyrol$station == as.character(station), ]
  rownames(RainData) <- c(1:NROW(RainData))
  
  # learning data: 24 years (1985 - 2008, both inlcuded)
  # testing data: 4 successive years (2009, 2010, 2011, 2012)
  learndata <- RainData[RainData$year < 2009,]
  testdata <- RainData[RainData$year %in% c(2009, 2010, 2011, 2012),]
  
  # define matrix to store computation times
  fit_time <- matrix(ncol = 5, nrow = 6)
  colnames(fit_time) <- c("user.self", "sys.self", "elapsed", "user.child", "sys.child")
  rownames(fit_time) <- c("disttree", "distforest", "gamlss", "gamboostLSS", "gamboostLSS_cvr", "EMOS")

  #####
  # fitting the model
  set.seed(7)
  
  # fit distributional tree
  fit_time["disttree",] <- system.time(dt <- disttree(dt.formula, 
                                                      data = learndata, family = dist_list_cens_normal, 
                                                      censtype = "left", censpoint = 0, type.tree = "ctree", 
                                                      control = ctree_control(teststat = "quad", testtype = "Bonferroni", intersplit = TRUE,
                                                                              mincriterion = 0.95, minsplit = 50,
                                                                              minbucket = 20)))
  
  # fit distributional forest
  fit_time["distforest",] <- system.time(df <- distforest(df.formula, 
                                                          data = learndata, family = dist_list_cens_normal, type.tree = "ctree", 
                                                          ntree = 100, censtype = "left", censpoint = 0, mtry = 27,
                                                          control = ctree_control(teststat = "quad", testtype = "Univ", intersplit = TRUE,
                                                                                  mincriterion = 0, minsplit = 50,
                                                                                  minbucket = 20)))
  
  # fit prespecified GAMLSS (covariates selected based on meteorological expert knowledge)
  g_learndata <- learndata
  g_learndata$robs <- Surv(g_learndata$robs, g_learndata$robs>0, type="left")
  fit_time["gamlss",] <- system.time(g <- try(gamlss(formula = g.mu.formula, sigma.formula = g.sigma.formula, data = g_learndata, 
                                                     family = cens("NO", type = "left"),
                                                     control = gamlss.control(n.cyc = 100),
                                                     i.control = glim.control(cyc = 100, bf.cyc = 100)),
                                              silent = TRUE))
  
  if(inherits(g, "try-error")){
    fit_time["gamlss",] <- NA
    warning("Error in gamlss, repeated with 50 itertations only (instead of 100)")
    g <- try(gamlss(formula = g.mu.formula, sigma.formula = g.sigma.formula, data = g_learndata, 
                    family = cens("NO", type = "left"),
                    control = gamlss.control(n.cyc = 50),
                    i.control = glim.control(cyc = 50, bf.cyc = 50)),
             silent = TRUE)
  }
  if(inherits(g, "try-error")){
    warning("Error in gamlss, returns NA instead of a model")
    g <- NA
  }

  # fit boosted GAMLSS
  fit_time["gamboostLSS",] <- system.time(gb <- gamboostLSS(formula = list(mu = gb.mu.formula, sigma = gb.sigma.formula), 
                                                            data = g_learndata, 
                                                            families =  as.families(fname = cens("NO", type = "left")), 
                                                            method = "noncyclic",
                                                            control = boost_control(mstop = 1000L)))
  
  # find optimal value for mstop (evalution is very time-consuming)
  grid <- seq(50,1000, by = 25)  
  fit_time["gamboostLSS_cvr",] <- system.time(cvr <- cvrisk(gb, grid = grid))
  mstop(gb) <- mstop(cvr)
  
  # fit linear model with only total precipitation as covariate (Ensemble Model Output Statistics, EMOS)
  fit_time["EMOS",] <- system.time(ml <- crch(formula = robs ~ tppow_mean | log(tppow_sprd + 0.001), 
                                              data = learndata, dist = "gaussian", 
                                              left = 0, link.scale = "log"))
    

  
  ## get predicted parameter of all models for testdata
  
  pred_time <- matrix(ncol = 5, nrow = NROW(testdata))
  colnames(pred_time) <- c("disttree", "distforest", "gamlss", "gamboostLSS", "EMOS")
  
  # distributional tree
  pdt <- data.frame()
  for(i in 1:NROW(testdata)){
    pred_time[i, "disttree"] <- system.time(predpar <- predict(dt, newdata = testdata[i,], type = "parameter"))["elapsed"]
    pdt <- rbind(pdt, predpar)
  }
  rownames(pdt) <- c(1:NROW(pdt))
  dt_mu <- pdt$mu
  dt_sigma <- pdt$sigma
  
  # distributional forest
  pdf <- data.frame()
  for(i in 1:NROW(testdata)){
    pred_time[i, "distforest"] <- system.time(predpar <- predict(df, newdata = testdata[i,], type = "parameter"))["elapsed"]
    pdf <- rbind(pdf, predpar)
  }
  df_mu <- pdf$mu
  df_sigma <- pdf$sigma
  
  # prespecified GAMLSS
  if(!is.na(g)){
    g_mu <- g_sigma <- numeric()
    for(i in 1:NROW(testdata)){
      pred_time[i, "gamlss"] <- system.time(predmu <- predict(g, newdata = testdata[i,], what = "mu", type = "response", data = g_learndata))["elapsed"] +
        system.time(predsigma <- predict(g, newdata = testdata[i,], what = "sigma", type = "response", data = g_learndata))["elapsed"]
      g_mu <- c(g_mu, predmu)
      g_sigma <- c(g_sigma, predsigma)
    }
  } else pred_time[, "gamlss"] <- g_mu <- g_sigma <- as.numeric(rep.int(NA, times = NROW(testdata)))
  pg <- data.frame(g_mu, g_sigma)
  colnames(pg) <- c("mu", "sigma")
  
  # boosted GAMLSS
  pgb <- data.frame()
  for(i in 1:NROW(testdata)){
    pred_time[i, "gamboostLSS"] <- system.time(predpar <- predict(gb, newdata = testdata[i,], parameter = c("mu","sigma"), type = "response"))["elapsed"]
    pgb <- rbind(pgb, predpar)
  }
  gb_mu <- pgb$mu
  gb_sigma <- pgb$sigma
  
  # EMOS
  ml_mu <- ml_sigma <- numeric()
  for(i in 1:NROW(testdata)){
    pred_time[i, "EMOS"] <- system.time(predmu <- predict(ml, type = "location", newdata = testdata[i,]))["elapsed"] +
      system.time(predsigma <- predict(ml, type = "scale", newdata = testdata[i,]))["elapsed"]
    ml_mu <- c(ml_mu, predmu)
    ml_sigma <- c(ml_sigma, predsigma)
  }
  pml <- data.frame(ml_mu, ml_sigma)
  colnames(pml) <- c("mu", "sigma")
  rownames(pml) <- c(1:NROW(pml))
  
  # store parameter
  par <- list(pdt = pdt,
              pdf = pdf,
              pg = pg,
              pgb = pgb,
              pml = pml,
              observations = testdata$robs)
  
  # CPRS
  crps_dt <- crps_cnorm(testdata$robs, location = dt_mu, scale = dt_sigma, lower = 0, upper = Inf)
  crps_df <- crps_cnorm(testdata$robs, location = df_mu, scale = df_sigma, lower = 0, upper = Inf)
  crps_g  <- if(!is.na(g)) crps_cnorm(testdata$robs, location = g_mu, scale = g_sigma, lower = 0, upper = Inf) else as.numeric(rep.int(NA, times = NROW(testdata)))
  crps_gb <- crps_cnorm(testdata$robs, location = gb_mu, scale = gb_sigma, lower = 0, upper = Inf)
  crps_ml <- crps_cnorm(testdata$robs, location = ml_mu, scale = ml_sigma, lower = 0, upper = Inf) 
  
  crps <- cbind(crps_dt, crps_df, crps_g, crps_gb, crps_ml)
  
  # loglikelihood
  dtll <- dfll <- gll <-  gbll <- mlll <- numeric(length = NROW(testdata))
  for(j in 1:(NROW(testdata))){
    
    eta_dt <- as.numeric(dist_list_cens_normal$linkfun(cbind(dt_mu, dt_sigma)[j,]))
    eta_df <- as.numeric(dist_list_cens_normal$linkfun(cbind(df_mu, df_sigma)[j,]))
    eta_g  <- if(!is.na(g)) as.numeric(dist_list_cens_normal$linkfun(cbind(g_mu, g_sigma)[j,])) else NA 
    eta_gb <- as.numeric(dist_list_cens_normal$linkfun(cbind(gb_mu, gb_sigma)[j,]))  
    eta_ml <- as.numeric(dist_list_cens_normal$linkfun(cbind(ml_mu, ml_sigma)[j,]))
    
    dtll[j] <- dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_dt, log=TRUE)
    dfll[j] <- dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_df, log=TRUE)
    gll[j]  <- if(!is.na(g)) dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_g, log=TRUE) else NA
    gbll[j] <- dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_gb, log=TRUE)
    mlll[j] <- dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_ml, log=TRUE)
    
  }

  ll <- cbind(dtll, dfll, gll, gbll, mlll) 

  colnames(ll) <- colnames(crps) <- c("disttree", "distforest", "gamlss", "gamboostLSS", "EMOS")

  
  
  ##  calculations for PIT histograms
  set.seed(7)
  
  # distributional tree
  pit_dt <- cbind(0, pnorm(testdata[,"robs"], mean = dt_mu, sd = dt_sigma))
  pit_dt[testdata[,"robs"]>0, 1] <- pit_dt[testdata[,"robs"]>0, 2]
  
  # distributional forest
  pit_df <- cbind(0, pnorm(testdata[,"robs"], mean = df_mu, sd = df_sigma))
  pit_df[testdata[,"robs"]>0, 1] <- pit_df[testdata[,"robs"]>0, 2]
  
  # prespecified GAMLSS
  if(!is.na(g)) {
    pit_g <- cbind(0, pnorm(testdata[,"robs"], mean = g_mu, sd = g_sigma))
    pit_g[testdata[,"robs"]>0, 1] <- pit_g[testdata[,"robs"]>0, 2]
  } else pit_g <- NA
  
  # boosted GAMLSS
  pit_gb <- cbind(0, pnorm(testdata[,"robs"], mean = gb_mu, sd = gb_sigma))
  pit_gb[testdata[,"robs"]>0, 1] <- pit_gb[testdata[,"robs"]>0, 2]
  
  # EMOS
  pit_ml <- cbind(0, pnorm(testdata[,"robs"], mean = ml_mu, sd = ml_sigma))
  pit_ml[testdata[,"robs"]>0, 1] <- pit_ml[testdata[,"robs"]>0, 2]


  
  ## Variable importance
  set.seed(7)
  nperm <- 50
  
  
  # function to permutate chosen variable and then calculate mean crps
  meancrps <- function(permute = NULL, newdata = testdata) {
    if(!is.null(permute)) newdata[[permute]] <- sample(newdata[[permute]])
    p <- predict(df, newdata = newdata, type = "parameter")
    mean(crps_cnorm(newdata$robs, location = p$mu, scale = p$sigma, lower = 0))
  }
  
  # apply for all covariates except for dswrf_mean_min and 
  # dswrf_sprd_min (columns 30 and 33) as they are always 0
  
  # using only one core
  # risk_all <- replicate(nperm, sapply(c(5:29, 31, 32, 34: ncol(testdata)), meancrps))
  # risk <- rowMeans(risk_all)
  
  # or parallel
  risklist <- applyfun(1:nperm, 
                       function(i){
                         set.seed(i)
                         sapply(c(5:29, 31, 32, 34: ncol(testdata)), meancrps)
                       })
  risk <- Reduce("+", risklist) / length(risklist)
  
  names(risk) <- names(testdata)[c(5:29, 31, 32, 34: ncol(testdata))]
  vimp_crps <- risk - meancrps(newdata = testdata)
  vimp_crps <- sort(vimp_crps, decreasing = TRUE)
  
  
  ## Cross validation
  nrep_cross <- 10
  seed <- 7
  
  res_cross <- applyfun(1:nrep_cross,
                        function(i){
                          
                          set.seed(seed*i)
                          
                          # randomly split data in 7 parts each including 4 years
                          years <- 1985:2012
                          testyears <- list()
                          for(j in 1:7){
                            testyears[[j]] <- sample(years, 4, replace = FALSE)
                            years <- years[!(years %in% testyears[[j]])]
                          }
                          
                          crps <- matrix(nrow = 7, ncol = 7)
                          reslist <- list()
                          for(k in 1:7){
                            test <- testyears[[k]]
                            train <- c(1985:2012)[!c(1985:2012) %in% test]
                            
                            res <- evalmodels(station = station,
                                              train = train,
                                              test = test,
                                              gamboost_cvr = TRUE,
                                              distfamily = "gaussian")
                            
                            crps[k,] <- res$crps
                            reslist[[k]] <- res
                          }
                          
                          colnames(crps) <- names(res$crps)
                          return(reslist)
                        }
  )
  
  
  results <- list(station = station,
                  par = par,
                  ll = ll,
                  crps = crps,
                  fit_time = fit_time,
                  pred_time = pred_time,
                  pit_dt = pit_dt,
                  pit_dt = pit_dt,
                  pit_df = pit_df,
                  pit_g = pit_g,
                  pit_gb = pit_gb,
                  pit_ml = pit_ml,
                  vimp_crps = vimp_crps,
                  res_cross = res_cross)
  
  return(results)
}
  




## evaluate for 15 selected stations and collect results in one .rda file ('stationwise.rda')
if(FALSE){
  
  stationlist <- c("Axams", "Lech", "Zuers", "See im Paznaun", "Jungholz", 
                   "Ladis-Neuegg", "Oetz", "Ochsengarten-Obergut",
                   "Ginzling", "Rotholz", "Walchsee", "Koessen", 
                   "Innervillgraten", "Matrei in Osttirol", 
                   "St.Johann im Walde")
  
  
  for(station in stationlist) {
    results <- stationeval(station = station)
    save(results, file = paste0("res_", 
                                gsub("-", "", gsub(".", "", gsub(" ", "", results$station, fixed = T), fixed = T), fixed = T),
                                ".rda"))
  }
  
  ## collect results
  stationwise <- list()
  for (i in Sys.glob("res_*.rda")) {
    load(i)
    stationwise[[results$station]] <- results
    # stationwise[[substr(i, 5, nchar(i) - 4)]] <- results
  }
  save(stationwise, file = "stationwise.rda")
  
  # remove individual .rda-files
  file.remove(Sys.glob("res_*.rda"))
}


## plot results for one selected station
if(FALSE){
  
  station <- "Lech"
  load("stationwise.rda")
  results <- stationwise[[station]]
  
  ## (i) map
  data("StationsTyrol", package = "RainTyrol")
  data("MapTyrol", package = "RainTyrol")
  library("sp")
  sp <- SpatialPointsDataFrame(subset(StationsTyrol, select = c(lon, lat)),
                               data = subset(StationsTyrol, select = -c(lon, lat)),
                               proj4string = raster::crs(MapTyrol$RasterLayer))
  
  
  ## convenience functions to unify results
  cross <- Reduce("c", results$res_cross)
  
  get_cross <- function(what = "crps", which = NULL) {
    ## methods
    if(is.null(which)) which <- c("distforest", "gamlss", "gamboostLSS", "gamboostLSS_cvr", "emos_log")
    if(what != "evaltime") which <- which[which != "gamboostLSS_cvr"]
    
    ## extract desired quantity (crps, evaltime, predtime)
    rval <- lapply(cross, "[[", what)
    rval <- if(what == "crps") {
      sapply(rval, function(x) x[which])
    } else {
      sapply(rval, function(x) x[which, "elapsed"])
    }
    
    ## reshape
    rval <- t(rval)
    rval <- aggregate(rval, list(rep(1:10, each = 7)), mean, na.rm = TRUE)[, -1]
    rval <- as.matrix(rval)
    which[which == "emos_log"] <- "EMOS"
    colnames(rval) <- which
    
    return(rval)
  }
  
  fmt <- function(x) {
    ## select and aggregate (if necessary)
    if(!is.null(dim(x))) {
      x <- if("elapsed" %in% colnames(x)) x[, "elapsed"] else colMeans(x, na.rm = TRUE)
    }
    
    ## round and format to three digits
    x <- format(round(x, digits = 3), nsmall = 3)
    
    ## collect results across methods (if necessary)
    if(length(x) == 1L) return(x)
    nam <- c("distforest", "gamlss", "gamboostLSS", "gamboostLSS_cvr", "EMOS")
    x <- x[nam]
    x[is.na(x)] <- ""
    names(x) <- nam
    return(x)
  }
  
  
  ## (ii) summary table
  tab <- cbind(
    "crps_24to4" =     fmt(results$crps),
    "fittime_24to4" =  fmt(results$fit_time),
    "predtime_24to4" = fmt(results$pred_time),
    "crps_cv" =        fmt(get_cross("crps")),
    "fittime_cv" =     fmt(get_cross("evaltime")),
    "predtime_cv" =    fmt(get_cross("predtime"))
  )
  
  
  ## (iii) crps skill score
  crps_cross <- get_cross("crps")
  crps_cross_mean <- colMeans(crps_cross, na.rm = TRUE)
  crpss_cross <- 1 - crps_cross[, c("distforest", "gamlss", "gamboostLSS")] / crps_cross[, "EMOS"]
  crpss_24to4 <- 1 - colMeans(results$crps[, c("distforest", "gamlss", "gamboostLSS")], na.rm = TRUE) / mean(results$crps[, "EMOS"], na.rm = TRUE)
  
  ## (iv) variable importance
  varimp10 <- rev(head(sort(results$vimp_crps, decreasing = TRUE), 10))
  names(varimp10) <- gsub("pow", "", names(varimp10), fixed = TRUE)
  
  ## Plot CRPS skill score
  par(mar = c(2.5, 4, 1, 2))
  boxplot(crpss_cross, 
          ylim = c(min(-0.01, min(crpss_cross, na.rm = TRUE) - 0.01), 
                   max(-0.01, max(crpss_cross, na.rm = TRUE) + 0.01)),
          main = "",
          names = c("Distributional\nforest", "Prespecified\nGAMLSS", "Boosted\nGAMLSS"),
          ylab = "CRPS skill score", col = "lightgray", las = 1) 
  #axis(1, 0:4, c("", "Distributional\nforest", "Prespecified\nGAMLSS", "Boosted\nGAMLSS", ""))
  abline(h = 0, col = hcl(50, 100, 50), lwd = 2)
  # lines(crpss_24to4, col = hcl(128, 100, 50), lwd = 2, type = "o", pch = 19)
  
  ## Plot Variable Importance
  par(mar = c(5, 10, 1, 2))
  barplot(varimp10, horiz = TRUE, las = 1, axes = FALSE,
          xlab = "Mean decrease in CRPS")
  axis(1, at = seq(0, ceiling(max(varimp10) * 5)/5, 0.02), las = 1, mgp = c(0, 1, 0))
  
    
  ## Plot residual QQ plots
  library("countreg")
  set.seed(7)
  par(mfrow = c(2, 2))
  qqrplot(results$pit_df, nsim = 100, main = "Distributional forest", ylim = c(-5, 5), col = gray(0.04, alpha = 0.01), pch = 19)
  qqrplot(results$pit_ml, nsim = 100, main = "EMOS",		    ylim = c(-5, 5), col = gray(0.04, alpha = 0.01), pch = 19)
  if(!all(is.na(results$pit_g))) qqrplot(results$pit_g,  nsim = 100, main = "Prespecified GAMLSS",   ylim = c(-5, 5), col = gray(0.04, alpha = 0.01), pch = 19)
  qqrplot(results$pit_gb, nsim = 100, main = "Boosted GAMLSS",	    ylim = c(-5, 5), col = gray(0.04, alpha = 0.01), pch = 19)

  ## Plot PIT histograms  
  set.seed(4)
  par(mfrow = c(2, 2))
  pithist(results$pit_df, nsim = 1000, breaks = seq(0, 1, length.out = 9), main = "Distributional forest", ylim = c(0, 1.9))
  pithist(results$pit_ml, nsim = 1000, breaks = seq(0, 1, length.out = 9), main = "EMOS", ylim = c(0, 1.9))
  if(!all(is.na(results$pit_g))) pithist(results$pit_g,  nsim = 1000, breaks = seq(0, 1, length.out = 9), main = "Prespecified GAMLSS", ylim = c(0, 1.9))
  pithist(results$pit_gb, nsim = 1000, breaks = seq(0, 1, length.out = 9), main = "boosted GAMLSS", ylim = c(0, 1.9))
  
}
