evalmodels <- function(station, train, test, 
                       ntree = 100,
                       tree_minsplit = 50, tree_minbucket = 20, tree_mincrit = 0.95,
                       forest_minsplit = 50, forest_minbucket = 20, forest_mincrit = 0,
                       forest_mtry = 27,
                       type.tree = "ctree",
                       gamboost_cvr = FALSE)
{
  cl <- match.call()
  
  ## check input arguments train and test
  if(any(!(c(train, test) %in% c(1985:2012)))) stop("training and testing years can only be chosen out of 1985-2012")
  ## in the 
  ## check station
  if(!is.character(station)) stop("argument 'station' has to be a character")
  
  ## load data
  data("RainTyrol", package = "RainTyrol")
  data("StationsTyrol", package = "RainTyrol")
  if(!(station %in% StationsTyrol$name)) stop("selected station is not among the 95 availble observation stations")
  rain <- RainTyrol[RainTyrol$station == station, ]
  
  #library("disttree")
  library("gamlss")
  library("gamlss.dist")
  library("gamlss.cens")
  gen.cens(NO, type = "left")
  library("gamboostLSS")
  
  
  ############
  # formula
  
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
  
  # gamlss formula
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
  
  
  # gamboostLSS formula
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
 
  
   
  ###################################################3
  # evaluation
  
  testdata <- rain[rain$year %in% test, ]
  learndata <- rain[rain$year %in% train, ]
  
  if(type.tree == "mob"){
    dt_time <- system.time(dt <- disttree::disttree(dt.formula, 
                                                    data = learndata, family = disttree::dist_list_cens_normal, 
                                                    censtype = "left", censpoint = 0, type.tree = "mob", 
                                                    control = partykit::mob_control(restart = FALSE, numsplit = "center", 
                                                                                    alpha = 1-tree_mincrit, minsplit = tree_minsplit,
                                                                                    minbucket = tree_minbucket)))
  }
  if(type.tree == "ctree"){
    dt_time <- system.time(dt <- disttree::disttree(dt.formula, 
                                                    data = learndata, family = disttree::dist_list_cens_normal, 
                                                    censtype = "left", censpoint = 0, type.tree = "ctree", 
                                                    control = partykit::ctree_control(teststat = "quad", testtype = "Bonferroni", intersplit = TRUE,
                                                                                      mincriterion = tree_mincrit, minsplit = tree_minsplit,
                                                                                      minbucket = tree_minbucket)))
  }
                          
                          
  if(type.tree == "mob"){
    df_time <- system.time(df <- disttree::distforest(df.formula, 
                                            data = learndata, family = disttree::dist_list_cens_normal, type.tree = "mob", 
                                            ntree = ntree, censtype = "left", censpoint = 0,
                                            control = partykit::mob_control(restart = FALSE, numsplit = "center", 
                                                                  alpha = 1-forest_mincrit, minsplit = forest_minsplit,
                                                                  minbucket = forest_minbucket), mtry = forest_mtry))
  }
  if(type.tree == "ctree"){
    df_time <- system.time(df <- disttree::distforest(df.formula, 
                                            data = learndata, family = disttree::dist_list_cens_normal, type.tree = "ctree", 
                                            ntree = ntree, censtype = "left", censpoint = 0, #fitted.OOB = FALSE,
                                            control = partykit::ctree_control(teststat = "quad", testtype = "Univ", intersplit = TRUE,
                                                                    mincriterion = forest_mincrit, minsplit = forest_minsplit,
                                                                    minbucket = forest_minbucket), mtry = forest_mtry))
  }
                          
                          
                          
                          
  g_learndata <- learndata
  g_learndata$robs <- survival::Surv(g_learndata$robs, g_learndata$robs>0, type="left")
  
  g_time <- system.time(g <- try(gamlss(formula = g.mu.formula, sigma.formula = g.sigma.formula, data = g_learndata, 
                                        family = cens("NO", type = "left"),
                                        control = gamlss.control(n.cyc = 100),
                                        i.control = glim.control(cyc = 100, bf.cyc = 100))))
  
  if(inherits(g, "try-error")) {
    g_time <- NA
    g <- NA
    g_error <- TRUE
  } else g_error <- FALSE
  
  
  #gb <- gb
  gb_time <- system.time(gb <- gamboostLSS(formula = list(mu = gb.mu.formula, sigma = gb.sigma.formula), data = g_learndata,
                                           families = as.families(fname = cens("NO", type = "left")), 
                                           method = "noncyclic",
                                           control = boost_control(mstop = 1000L)))
  if(gamboost_cvr){
    grid <- seq(50, 1000, by = 25)
    gb_cvr_time <- system.time(cvr <- cvrisk(gb, grid = grid))
    mstop(gb) <- mstop(cvr)
    cvr_opt <- mstop(cvr) 
  } else cvr_opt <- gb_cvr_time <- NA
  
                          
  mi_time <- system.time(mi <- try(crch::crch(formula = robs ~ tppow_mean | tppow_sprd, 
                                        data = learndata, dist = "gaussian", left = 0, link.scale = "identity")))
  if(inherits(mi, "try-error")) {
    mi_time <- NA
    mi <- NA
    mi_error <- TRUE
  } else mi_error <- FALSE
  
  ml_time <- system.time(ml <- try(crch::crch(formula = robs ~ tppow_mean | log(tppow_sprd + 0.001), 
                                        data = learndata, dist = "gaussian", left = 0, link.scale = "log")))
  if(inherits(ml, "try-error")) {
    ml_time <- NA
    ml <- NA
    ml_error <- TRUE
  } else ml_error <- FALSE
  
  mq_time <- system.time(mq <- try(crch::crch(formula = robs ~ tppow_mean | I(tppow_sprd^2), 
                                        data = learndata, dist = "gaussian", left = 0, link.scale = "quadratic")))
  if(inherits(mq, "try-error")) {
    mq_time <- NA
    mq <- NA
    mq_error <- TRUE
  } else mq_error <- FALSE
  
  
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
      g_error <- TRUE
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
  crps_dt <- mean(scoringRules::crps_cnorm(testdata$robs, location = dt_mu, scale = dt_sigma, lower = 0, upper = Inf), na.rm = TRUE)
  crps_df <- mean(scoringRules::crps_cnorm(testdata$robs, location = df_mu, scale = df_sigma, lower = 0, upper = Inf), na.rm = TRUE)
  crps_g <- if(!g_na) mean(scoringRules::crps_cnorm(testdata$robs, location = g_mu, scale = g_sigma, lower = 0, upper = Inf), na.rm = TRUE) else NA
  crps_gb <- mean(scoringRules::crps_cnorm(testdata$robs, location = gb_mu, scale = gb_sigma, lower = 0, upper = Inf), na.rm = TRUE)
  crps_mi <- if(!mi_na) mean(scoringRules::crps_cnorm(testdata$robs, location = mi_mu, scale = mi_sigma, lower = 0, upper = Inf), na.rm = TRUE) else NA
  crps_ml <- if(!ml_na) mean(scoringRules::crps_cnorm(testdata$robs, location = ml_mu, scale = ml_sigma, lower = 0, upper = Inf), na.rm = TRUE) else NA
  crps_mq <- if(!mq_na) mean(scoringRules::crps_cnorm(testdata$robs, location = mq_mu, scale = mq_sigma, lower = 0, upper = Inf), na.rm = TRUE) else NA
  
  # RMSE
  rmse_dt <- sqrt(mean((dt_exp - testdata[,"robs"])^2, na.rm = TRUE))
  rmse_df <- sqrt(mean((df_exp - testdata[,"robs"])^2, na.rm = TRUE))
  rmse_g <- if(!g_na) sqrt(mean((g_exp - testdata[,"robs"])^2, na.rm = TRUE)) else NA
  rmse_gb <- sqrt(mean((gb_exp - testdata[,"robs"])^2, na.rm = TRUE))
  rmse_mi <- if(!mi_na) sqrt(mean((mi_exp - testdata[,"robs"])^2, na.rm = TRUE)) else NA
  rmse_ml <- if(!ml_na) sqrt(mean((ml_exp - testdata[,"robs"])^2, na.rm = TRUE)) else NA
  rmse_mq <- if(!mq_na) sqrt(mean((mq_exp - testdata[,"robs"])^2, na.rm = TRUE)) else NA
  
  # loglikelihood
  dtll <- dfll <- gll <-  gbll <- mill <- mlll <- mqll <- numeric(length = NROW(testdata))
  for(j in 1:(NROW(testdata))){
    
    eta_dt <- as.numeric(disttree::dist_list_cens_normal$linkfun(cbind(dt_mu, dt_sigma)[j,]))
    eta_df <- as.numeric(disttree::dist_list_cens_normal$linkfun(cbind(df_mu, df_sigma)[j,]))
    eta_g <- if(!g_na) as.numeric(disttree::dist_list_cens_normal$linkfun(cbind(g_mu, g_sigma)[j,])) else NA
    eta_gb <- as.numeric(disttree::dist_list_cens_normal$linkfun(cbind(gb_mu, gb_sigma)[j,]))
    eta_mi <- if(!mi_na) as.numeric(disttree::dist_list_cens_normal$linkfun(cbind(mi_mu, mi_sigma)[j,])) else NA
    eta_ml <- if(!ml_na) as.numeric(disttree::dist_list_cens_normal$linkfun(cbind(ml_mu, ml_sigma)[j,])) else NA
    eta_mq <- if(!mq_na) as.numeric(disttree::dist_list_cens_normal$linkfun(cbind(mq_mu, mq_sigma)[j,])) else NA
    
    dtll[j] <- disttree::dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_dt, log=TRUE)
    dfll[j] <- disttree::dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_df, log=TRUE)
    gll[j] <- if(!g_na) disttree::dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_g, log=TRUE) else NA
    gbll[j] <- disttree::dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_gb, log=TRUE)
    mill[j] <- if(!mi_na) disttree::dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_mi, log=TRUE) else NA
    mlll[j] <- if(!mi_na) disttree::dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_ml, log=TRUE) else NA
    mqll[j] <- if(!mi_na) disttree::dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_mq, log=TRUE) else NA
    
  }
  
  ll_dt <- mean(dtll, na.rm = TRUE)
  ll_df <- mean(dfll, na.rm = TRUE)
  ll_g <- if(!g_na) mean(gll, na.rm = TRUE) else NA
  ll_gb <- mean(gbll, na.rm = TRUE)
  ll_mi <- if(!mi_na) mean(mill, na.rm = TRUE) else NA      
  ll_ml <- if(!ml_na) mean(mlll, na.rm = TRUE) else NA 
  ll_mq <- if(!mq_na) mean(mqll, na.rm = TRUE) else NA 
  
  evaltime <- data.frame(
    user <- c(dt_time[1], df_time[1], g_time[1], gb_time[1], gb_cvr_time[1], mi_time[1], ml_time[1], mq_time[1]),
    system <- c(dt_time[2], df_time[2], g_time[2], gb_time[2], gb_cvr_time[2], mi_time[2], ml_time[2], mq_time[2]),             
    elapsed <- c(dt_time[3], df_time[3], g_time[3], gb_time[3], gb_cvr_time[3], mi_time[3], ml_time[3], mq_time[3]))
  colnames(evaltime) <- c("user", "system", "elapsed")
  rownames(evaltime) <- c("disttree", "distforest", "gamlss", "gamboostLSS", "gamboostLSS_cvr", "emos_id", "emos_log", "emos_quad")
  
  
  ## store results in a list
  res <- list()
  res$call <- cl
  
  res$ls <- c(ll_dt, ll_df, ll_g, ll_gb, ll_mi, ll_ml, ll_mq)
  res$rmse <- c(rmse_dt, rmse_df, rmse_g, rmse_gb, rmse_mi, rmse_ml, rmse_mq)
  res$crps <- c(crps_dt, crps_df, crps_g, crps_gb, crps_mi, crps_ml, crps_mq)
  names(res$rmse) <- names(res$ls) <- names(res$crps) <- 
    c("disttree", "distforest", "gamlss", "gamboostLSS", "emos_id", "emos_log", "emos_quad")
  
  res$cvr_opt <- cvr_opt
  res$evaltime <- evaltime
  
  # store errors if any occured
  if(g_error) res$g_error <- g_error
  if(mi_error) res$mi_error <- mi_error
  if(ml_error) res$ml_error <- ml_error
  if(mq_error) res$mq_error <- mq_error
  
  return(res)
}
 