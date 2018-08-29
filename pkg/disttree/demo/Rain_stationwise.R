#######################################################
### Probabilistic Forecasting on Precipitation Data ###
#######################################################

## Replication material for: 
## Distributional Regression Forests for Probabilistic Precipitation Forecasting in Complex Terrain (2018)
## by Lisa Schlosser and Torsten Hothorn and Reto Stauffer and Achim Zeileis
## URL: http://arxiv.org/abs/1804.02921

## This demo includes the application on selected stations
## (models learned on 24 years and evaluated on 4 years)
## Full replication of all other results can be obtained with
## demo("RainTyrol", package = "disttree")

## Computation time: approximately 18 minutes (on our machines, using 1 kernel)


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
gen.cens(LO, type = "left")

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
applyfun <- function(X, FUN, ...) parallel::mclapply(X, FUN, ..., mc.cores = pmax(1, detectCores() - 1))


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
  
  # formula for boosted GAM if used as replacement of gamlss in case of an error
  g_gb.mu.formula <- robs ~ bbs(tppow_mean) + 
    bspatial(tppow_mean1218, capepow_mean1218) + 
    #bbs(tppow_mean1218) + bbs(capepow_mean1218) +
    bbs(tppow_max) + 
    bbs(dswrf_mean_mean) +
    bbs(tcolc_mean_mean) + 
    bbs(msl_diff) + 
    bbs(pwat_mean_mean) + 
    bbs(tdiff500850_mean)  
  
  g_gb.sigma.formula <- robs ~ bbs(tppow_sprd) + 
    bspatial(tppow_sprd1218, capepow_mean1218) + 
    #bbs(tppow_sprd1218) + bbs(capepow_mean1218) + 
    bbs(dswrf_sprd_mean) +
    bbs(tcolc_sprd_mean) + 
    bbs(tdiff500850_mean)
  
  
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




stationeval <- function(station) {
  
  #####
  # load observations and covariates 
  load(paste("~/svn/partykit/pkg/disttree/data/Rain", station, ".rda", sep = ""))
  RainData <- get(paste("Rain", station, sep = ""))
  
  # learning data: 24 years (1985 - 2008, both inlcuded)
  # testing data: 4 successive years (2009, 2010, 2011, 2012)
  learndata <- RainData[RainData$year < 2009,]
  testdata <- RainData[RainData$year %in% c(2009, 2010, 2011, 2012),]
  
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
  
  # fit prespecified GAM (covariates selected based on meteorological expert knowledge)
  g_learndata <- learndata
  g_learndata$robs <- Surv(g_learndata$robs, g_learndata$robs>0, type="left")
  fit_time["gamlss",] <- system.time(g <- try(gamlss(formula = g.mu.formula, sigma.formula = g.sigma.formula, data = g_learndata, 
                                            family = cens("NO", type = "left"),
                                            control = gamlss.control(n.cyc = 100),
                                            i.control = glim.control(cyc = 100, bf.cyc = 100)),
                                     silent = TRUE))
  
  # if(inherits(g, "try-error")) g <- gamboostLSS(formula = list(mu = g_gb.mu.formula, sigma = g_gb.sigma.formula), 
  #                    data = g_learndata, 
  #                    families =  as.families(fname = cens("NO", type = "left")), 
  #                    method = "noncyclic",
  #                    control = boost_control(mstop = 2000L))
  
  # compare
  # predg <- cbind(predict(g, newdata = testdata, what = "mu", type = "response"), predict(g, newdata = testdata, what = "sigma", type = "response")
  # predg2 <- cbind(predict(g2, newdata = testdata, parameter = "mu", type = "response"), predict(g2, newdata = testdata, parameter = "sigma", type = "response"))
  
  
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
    warning("Error in gamlss, repeated with 10 itertations only (instead of 100)")
    g <- try(gamlss(formula = g.mu.formula, sigma.formula = g.sigma.formula, data = g_learndata, 
                    family = cens("NO", type = "left"),
                    control = gamlss.control(n.cyc = 10),
                    i.control = glim.control(cyc = 10, bf.cyc = 10)),
             silent = TRUE)
  }
  if(inherits(g, "try-error")){
    warning("Error in gamlss, repeated with 5 itertations only (instead of 100)")
    g <- gamlss(formula = g.mu.formula, sigma.formula = g.sigma.formula, data = g_learndata, 
                family = cens("NO", type = "left"),
                control = gamlss.control(n.cyc = 5),
                i.control = glim.control(cyc = 5, bf.cyc = 5))
  }

  # fit boosted GAM
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
  ## FIX ME: for loop very time consuming!!!
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
  
  # prespecified GAM
  g_mu <- g_sigma <- numeric()
  for(i in 1:NROW(testdata)){
    pred_time[i, "gamlss"] <- system.time(predmu <- predict(g, newdata = testdata[i,], what = "mu", type = "response", data = g_learndata))["elapsed"] +
      system.time(predsigma <- predict(g, newdata = testdata[i,], what = "sigma", type = "response", data = g_learndata))["elapsed"]
    g_mu <- c(g_mu, predmu)
    g_sigma <- c(g_sigma, predsigma)
  }
  pg <- data.frame(g_mu, g_sigma)
  colnames(pg) <- c("mu", "sigma")

  # boosted GAM
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
  crps_g  <- crps_cnorm(testdata$robs, location = g_mu, scale = g_sigma, lower = 0, upper = Inf)
  crps_gb <- crps_cnorm(testdata$robs, location = gb_mu, scale = gb_sigma, lower = 0, upper = Inf)
  crps_ml <- crps_cnorm(testdata$robs, location = ml_mu, scale = ml_sigma, lower = 0, upper = Inf) 
  
  crps <- cbind(crps_dt, crps_df, crps_g, crps_gb, crps_ml)
  
  # loglikelihood
  dtll <- dfll <- gll <-  gbll <- mlll <- numeric(length = NROW(testdata))
  for(j in 1:(NROW(testdata))){
    
    eta_dt <- as.numeric(dist_list_cens_normal$linkfun(cbind(dt_mu, dt_sigma)[j,]))
    eta_df <- as.numeric(dist_list_cens_normal$linkfun(cbind(df_mu, df_sigma)[j,]))
    eta_g  <- as.numeric(dist_list_cens_normal$linkfun(cbind(g_mu, g_sigma)[j,])) 
    eta_gb <- as.numeric(dist_list_cens_normal$linkfun(cbind(gb_mu, gb_sigma)[j,]))  
    eta_ml <- as.numeric(dist_list_cens_normal$linkfun(cbind(ml_mu, ml_sigma)[j,]))
    
    dtll[j] <- dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_dt, log=TRUE)
    dfll[j] <- dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_df, log=TRUE)
    gll[j]  <- dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_g, log=TRUE)
    gbll[j] <- dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_gb, log=TRUE)
    mlll[j] <- dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_ml, log=TRUE)
    
  }

  ll <- cbind(dtll, dfll, gll, gbll, mlll) 

  colnames(ll) <- colnames(crps) <- c("disttree", "distforest", "gamlss", "gamboostLSS", "EMOS")
  
  #results <- data.frame(ll, crps)
  
  
  
  ##  PIT histograms
  
  set.seed(7)
  
  # distributional tree
  pit_dt <- cbind(0, pnorm(testdata[,"robs"], mean = dt_mu, sd = dt_sigma))
  pit_dt[testdata[,"robs"]>0, 1] <- pit_dt[testdata[,"robs"]>0, 2]
  
  # distributional forest
  pit_df <- cbind(0, pnorm(testdata[,"robs"], mean = df_mu, sd = df_sigma))
  pit_df[testdata[,"robs"]>0, 1] <- pit_df[testdata[,"robs"]>0, 2]
  
  # prespecified GAM
  pit_g <- cbind(0, pnorm(testdata[,"robs"], mean = g_mu, sd = g_sigma))
  pit_g[testdata[,"robs"]>0, 1] <- pit_g[testdata[,"robs"]>0, 2]
  
  # boosted GAM
  pit_gb <- cbind(0, pnorm(testdata[,"robs"], mean = gb_mu, sd = gb_sigma))
  pit_gb[testdata[,"robs"]>0, 1] <- pit_gb[testdata[,"robs"]>0, 2]
  
  # EMOS
  pit_ml <- cbind(0, pnorm(testdata[,"robs"], mean = ml_mu, sd = ml_sigma))
  pit_ml[testdata[,"robs"]>0, 1] <- pit_ml[testdata[,"robs"]>0, 2]


  
  ## Variable importance
  
  set.seed(7)
  
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
  
  vimp_crps <- varimp(df, nperm = 1L)
  
  rm(logLik.distforest) 
  
  
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
  

station <- "Axams"

results <- stationeval(station = station)
save(results, file = paste0("~/svn/partykit/pkg/disttree/demo/results_stationwise/res_", results$station, ".rda"))



## plot results
if(FALSE){
  
  # PIT histograms
  set.seed(4)
  par(mfrow = c(2,2))
  library("countreg")
  
  pithist(results$pit_df, nsim = 1000, breaks = seq(0, 1, length.out = 9), main = "Distributional forest", ylim = c(0,1.5))
  pithist(results$pit_ml, nsim = 1000, breaks = seq(0, 1, length.out = 9), main = "EMOS", ylim = c(0,1.5))
  pithist(results$pit_g,  nsim = 1000, breaks = seq(0, 1, length.out = 9), main = "Prespecified GAMLSS", ylim = c(0,1.5))
  pithist(results$pit_gb, nsim = 1000, breaks = seq(0, 1, length.out = 9), main = "boosted GAMLSS", ylim = c(0,1.5))
  
  
  
  # QQR plots
  set.seed(7)
  par(mfrow = c(2, 2))
  qqrplot(results$pit_df, nsim = 100, main = "Distributional forest", ylim = c(-5, 5), col = gray(0.04, alpha = 0.01), pch = 19)
  qqrplot(results$pit_ml, nsim = 100, main = "EMOS",                  ylim = c(-5, 5), col = gray(0.04, alpha = 0.01), pch = 19)
  qqrplot(results$pit_g,  nsim = 100, main = "Prespecified GAMLSS",   ylim = c(-5, 5), col = gray(0.04, alpha = 0.01), pch = 19)
  qqrplot(results$pit_gb, nsim = 100, main = "Boosted GAMLSS",        ylim = c(-5, 5), col = gray(0.04, alpha = 0.01), pch = 19)
  
  
  # variable importance: plot top 10
  par(mfrow = c(1, 1), mar = c(5, 10, 2, 4))
  barplot(sort(results$vimp_crps, decreasing = FALSE)[(length(results$vimp_crps)-9):length(results$vimp_crps)], 
          horiz = TRUE, las = 1, axes = FALSE,
          xlab = "Variable importance: mean decrease in CRPS",
          font.axis = 3, #list(family="HersheySerif", face=3),
          names.arg = gsub("pow", "", names(sort(results$vimp_crps, decreasing = FALSE)[(length(results$vimp_crps)-9):length(results$vimp_crps)])))
  axis(1, at = seq(0,1.6,0.2), las = 1, mgp=c(0,1,0))
  
  
  
  ## cross validation
  # extract CRPS 
  crps_cross <- matrix(nrow = 10, ncol = 7)
  # loop over all repetitions
  for(i in 1:length(results$res_cross)){
    #loop over all 7 folds (for 7 methods)
    crps_cross_int <- matrix(nrow = length(results$res_cross[[1]]), ncol = 7)
    for(j in 1:length(results$res_cross[[1]])){
      crps_cross_int[j,] <- results$res_cross[[i]][[j]]$crps
    }
    crps_cross[i,] <- colMeans(crps_cross_int, na.rm = TRUE)
  }
  colnames(crps_cross) <- names(results$res_cross[[1]][[1]]$crps) 
  
  #save(crps_cross, file = "crps_cross.rda")
  
  ## boxplot of crps_cross for station Axams
  boxplot(1 - crps_cross[,c(2,3,4)] / crps_cross[,6], ylim = c(-0.005, 0.065),
          names = c("Distributional forest", "Prespecified GAMLSS", "Boosted GAMLSS"),
          ylab = "CRPS skill score", col = "lightgray") 
  abline(h = 0, col = pal[5], lwd = 2)
  
  
}