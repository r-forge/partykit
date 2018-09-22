#######################################################
### Probabilistic Forecasting on Precipitation Data ###
#######################################################

## Replication material for supplement of: 
## Distributional Regression Forests for Probabilistic Precipitation Forecasting in Complex Terrain (2018)
## by Lisa Schlosser and Torsten Hothorn and Reto Stauffer and Achim Zeileis
## URL: http://arxiv.org/abs/1804.02921

## This demo includes the application on selected stations
## (models learned on 24 years and evaluated on 4 years)
## Full replication of all other results can be obtained with
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
gen.cens(LO, type = "left")

# if gamlss.tr family object should be used as family
library("gamlss.tr")
gen.trun(0, "NO", type = "left")

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

assign("LO",  gamlss.dist::NO,  pos = ".GlobalEnv")
assign("dLO", gamlss.dist::dNO, pos = ".GlobalEnv")
assign("pLO", gamlss.dist::pNO, pos = ".GlobalEnv")
assign("qLO", gamlss.dist::qNO, pos = ".GlobalEnv")
assign("rLO", gamlss.dist::rNO, pos = ".GlobalEnv")
gamlss.cens::gen.cens(NO, type = "left")
assign("LOlc",  NOlc,  pos = ".GlobalEnv")
assign("dLOlc", dNOlc, pos = ".GlobalEnv")
assign("pLOlc", pNOlc, pos = ".GlobalEnv")
assign("qLOlc", qNOlc, pos = ".GlobalEnv")

# distribution list for left-censored at zero logistic distribution
dist_list_cens_log <- dist_crch(dist = "logistic", type = "left", censpoint = 0)

##### 
# HCL palette
pal <- hcl(c(10, 128, 260, 290, 50), 100, 50)


#####

# formula for gaussian and logistic model and for truncated part of the hgaussian model
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



# formula for binomial model
{
  
  # tree and forest formula
  dt.formula_bin <- df.formula_bin <- 
    robs>0 ~ tppow_mean + tppow_sprd + tppow_min + tppow_max + 
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
  g.mu.formula_bin <- as.factor(robs>0) ~ pb(tppow_mean) + 
    pb(tppow_mean1218 * capepow_mean1218) + 
    pb(tppow_max) + 
    pb(dswrf_mean_mean) +
    pb(tcolc_mean_mean) + 
    pb(msl_diff) + 
    pb(pwat_mean_mean) + 
    pb(tdiff500850_mean)  
  
  
  # formula for boosted GAM
  gb.mu.formula_bin <-  
    as.factor(robs>0) ~ bbs(tppow_mean) + bbs(tppow_sprd) + bbs(tppow_min) + bbs(tppow_max) + 
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


stationeval <- function(station, method, distribution)
  {
  
  #####
  # get observations and covariates for selected station
  RainData <- RainTyrol[RainTyrol$station == as.character(station), ]
  rownames(RainData) <- c(1:NROW(RainData))
  
  # learning data: 24 years (1985 - 2008, both inlcuded)
  # testing data: 4 successive years (2009, 2010, 2011, 2012)
  learndata <- RainData[RainData$year < 2009,]
  testdata <- RainData[RainData$year %in% c(2009, 2010, 2011, 2012),]
  
  if(distribution == "hgaussian"){
    learndata_pos <- learndata[learndata$robs>0,]
    testdata_pos <- testdata[testdata$robs>0,]
  }

  
  #####
  # fitting the model
  set.seed(7)
  
  if(distribution %in% c("gaussian", "logistic")){
    
    if(method == "disttree"){
      dt <- disttree(dt.formula, 
                     data = learndata, 
                     family = if(distribution == "gaussian") dist_list_cens_normal else dist_list_cens_log, 
                     censtype = "left", censpoint = 0, type.tree = "ctree", 
                     control = ctree_control(teststat = "quad", testtype = "Bonferroni", intersplit = TRUE,
                                             mincriterion = 0.95, minsplit = 50,
                                             minbucket = 20))
      predpar <- predict(dt, newdata = testdata, type = "parameter")

    }
    
    if(method == "distforest"){
      df <- distforest(df.formula, 
                       data = learndata, 
                       family = if(distribution == "gaussian") dist_list_cens_normal else dist_list_cens_log, 
                       type.tree = "ctree", 
                       ntree = 100, censtype = "left", censpoint = 0, mtry = 27,
                       control = ctree_control(teststat = "quad", testtype = "Univ", intersplit = TRUE,
                                               mincriterion = 0, minsplit = 50,
                                               minbucket = 20))
      predpar <- predict(df, newdata = testdata, type = "parameter")
    }
    
    if(method == "gamlss"){
      g_learndata <- learndata
      g_learndata$robs <- Surv(g_learndata$robs, g_learndata$robs>0, type="left")
      g <- try(gamlss(formula = g.mu.formula, sigma.formula = g.sigma.formula, data = g_learndata, 
                  family = if(distribution == "gaussian") cens("NO", type = "left") else cens("LO", type = "left"),
                  control = gamlss.control(n.cyc = 100),
                  i.control = glim.control(cyc = 100, bf.cyc = 100)),
               silent = TRUE)
      if(inherits(g, "try-error")){
        warning("Error in gamlss, repeated with 50 itertations only (instead of 100)")
        g <- try(gamlss(formula = g.mu.formula, sigma.formula = g.sigma.formula, data = g_learndata, 
                        family = if(distribution == "gaussian") cens("NO", type = "left") else cens("LO", type = "left"),
                        control = gamlss.control(n.cyc = 50),
                        i.control = glim.control(cyc = 50, bf.cyc = 50)),
                 silent = TRUE)
      }
      if(inherits(g, "try-error")){
        g <- NA
      }
      
      if(!is.na(g)){
        predpar <- data.frame(predict(g, newdata = testdata, what = "mu", type = "response", data = g_learndata),
                              predict(g, newdata = testdata, what = "sigma", type = "response", data = g_learndata))
      } else {
        predpar <- data.frame(as.numeric(rep.int(NA, times = NROW(testdata))),
                              as.numeric(rep.int(NA, times = NROW(testdata))))
      }
      colnames(predpar) <- c("mu", "sigma")
    }
    
    if(method == "gamboostLSS"){
      g_learndata <- learndata
      g_learndata$robs <- Surv(g_learndata$robs, g_learndata$robs>0, type="left")
      if(distribution == "gaussian"){
        gb <- gamboostLSS(formula = list(mu = gb.mu.formula, sigma = gb.sigma.formula), 
                          data = g_learndata, 
                          families =  as.families(fname = cens("NO", type = "left")), 
                          method = "noncyclic",
                          control = boost_control(mstop = 1000L))
      }
      if(distribution == "logistic"){
        gb <- gamboostLSS(formula = list(mu = gb.mu.formula, sigma = gb.sigma.formula), 
                          data = g_learndata, 
                          families =  as.families(fname = cens("LO", type = "left")), 
                          method = "noncyclic",
                          control = boost_control(mstop = 1000L))
      }
      # find optimal value for mstop (evalution is very time-consuming)
      grid <- seq(50,1000, by = 25)  ## FIX ME: changed stepsize from 25 to 100 due to memory problems (Error in mcfork(detached) : unable to fork, possible reason: Cannot allocate memory)"
      cvr <- cvrisk(gb, grid = grid)
      mstop(gb) <- mstop(cvr)
      predpar <- data.frame(predict(gb, newdata = testdata, parameter = "mu", type = "response"),
                            predict(gb, newdata = testdata, parameter = "sigma", type = "response"))
      colnames(predpar) <- c("mu", "sigma")
    }
    
    if(method == "EMOS"){
      ml <- crch(formula = robs ~ tppow_mean | log(tppow_sprd + 0.001), 
                 data = learndata, dist = as.character(distribution), 
                 left = 0, link.scale = "log")
      predpar <- data.frame(predict(ml, type = "location", newdata = testdata),
                            predict(ml, type = "scale", newdata = testdata))
      colnames(predpar) <- c("mu", "sigma")
      
    }
  }
  
  
  
  if(distribution == "hgaussian"){
    
    if(method == "disttree"){
      # first a binomial tree to decide whether there is any precipitation at all
      dt_bin <- disttree(dt.formula_bin, 
                         data = learndata,
                         family = dist_binomial, 
                         type.tree = "ctree", 
                         control = ctree_control(teststat = "quad", testtype = "Bonferroni", intersplit = TRUE,
                                                 mincriterion = 0.95, minsplit = 50,
                                                 minbucket = 20))
      
      # fit a tree on the positive observations only (with truncated normal distribution)
      dt_tr <- disttree(dt.formula, 
                        data = learndata_pos,
                        family = dist_list_trunc_normal, 
                        type.tree = "ctree", 
                        control = ctree_control(teststat = "quad", testtype = "Bonferroni", intersplit = TRUE,
                                                mincriterion = 0.95, minsplit = 50,
                                                minbucket = 20))
      
      predpar <- predict(dt_tr, newdata = testdata, type = "parameter")
      predpar$nu <- predict(dt_bin, newdata = testdata, type = "parameter")$mu
      
    }
    
    if(method == "distforest"){
      # first a binomial tree to decide whether there is any precipitation at all
      df_bin <- distforest(df.formula_bin, 
                           data = learndata, 
                           family = dist_binomial, type.tree = "ctree", 
                           ntree = 100, mtry = 27,
                           control = ctree_control(teststat = "quad", testtype = "Univ", intersplit = TRUE,
                                                   mincriterion = 0, minsplit = 50,
                                                   minbucket = 20))
      
      
      # fit a forest on the positive observations only (with truncated normal distribution)
      df_tr <- distforest(df.formula, 
                          data = learndata_pos, 
                          family = dist_list_trunc_normal, type.tree = "ctree", 
                          ntree = 100, mtry = 27,
                          control = ctree_control(teststat = "quad", testtype = "Univ", intersplit = TRUE,
                                                  mincriterion = 0, minsplit = 50,
                                                  minbucket = 20))
      
      predpar <- predict(df_tr, newdata = testdata, type = "parameter")
      predpar$nu <- predict(df_bin, newdata = testdata, type = "parameter")$mu
    }
    
    if(method == "gamlss"){
      g_bin <- try(gamlss(formula = g.mu.formula_bin, 
                          data = learndata, 
                          family = BI(),
                          control = gamlss.control(n.cyc = 100),
                          i.control = glim.control(cyc = 100, bf.cyc = 100)), silent = TRUE)
      if(inherits(g_bin, "try-error")){
        warning("Error in gamlss for binomial model, repeated with 50 itertations only (instead of 100)")
        g_bin <- try(gamlss(formula = g.mu.formula_bin, 
                            data = learndata, 
                            family = BI(),
                            control = gamlss.control(n.cyc = 50),
                            i.control = glim.control(cyc = 50, bf.cyc = 50)), silent = TRUE)
      }
      if(inherits(g_bin, "try-error")){
        g_bin <- NA
      }
      
      
      g_tr <- try(gamlss(formula = g.mu.formula, sigma.formula = g.sigma.formula, 
                         data = learndata_pos, 
                         family = trun(0, "NO", type = "left"),
                         control = gamlss.control(n.cyc = 100),
                         i.control = glim.control(cyc = 100, bf.cyc = 100)), silent = TRUE)
      if(inherits(g_tr, "try-error")){
        warning("Error in gamlss for truncated model, repeated with 50 itertations only (instead of 100)")
        g_tr <- try(gamlss(formula = g.mu.formula, sigma.formula = g.sigma.formula, 
                           data = learndata_pos, 
                           family = trun(0, "NO", type = "left"),
                           control = gamlss.control(n.cyc = 50),
                           i.control = glim.control(cyc = 50, bf.cyc = 50)), silent = TRUE)
      }
      if(inherits(g_tr, "try-error")){
        g_tr <- NA
      }
      
      if(!(all(is.na(g_bin)) | all(is.na(g_tr)))){
        predpar <- data.frame(predict(g_tr, newdata = testdata, what = "mu", type = "response", data = learndata_pos),
                              predict(g_tr, newdata = testdata, what = "sigma", type = "response", data = learndata_pos))
        predpar$nu <- as.numeric(predict(g_bin, newdata = testdata, what = "mu", type = "response", data = learndata))
      } else {
        predpar <- data.frame(as.numeric(rep.int(NA, times = NROW(testdata))),
                              as.numeric(rep.int(NA, times = NROW(testdata))),
                              as.numeric(rep.int(NA, times = NROW(testdata))))
      }
    }
    
    if(method == "gamboostLSS"){
      gb_bin <- mboost::mboost(formula = gb.mu.formula_bin, 
                               data = learndata, 
                               family = Binomial(type = "adaboost",
                                                 link = "logit"), 
                               control = boost_control(mstop = 1000L))
      # find optimal value for mstop (evalution is very time-consuming)
      grid <- seq(50,1000, by = 25)
      cvr_bin <- cvrisk(gb_bin, grid = grid)
      mstop(gb_bin) <- mstop(cvr_bin)
      
      gb_tr <- gamboostLSS(formula = list(mu = gb.mu.formula, sigma = gb.sigma.formula), 
                           data = learndata_pos, 
                           families = as.families(fname = trun(0, "NO", type = "left")), method = "noncyclic",
                           control = boost_control(mstop = 1000L))
      # find optimal value for mstop (evalution is very time-consuming)
      grid <- seq(50,1000, by = 25)
      cvr_tr <- cvrisk(gb_tr, grid = grid)
      mstop(gb_tr) <- mstop(cvr_tr)
      predpar <- data.frame(predict(gb_tr, newdata = testdata, parameter = "mu", type = "response"),
                            predict(gb_tr, newdata = testdata, parameter = "sigma", type = "response"))
      predpar$nu <- as.numeric(predict(gb_bin, newdata = testdata, parameter = c("mu"), type = "response"))
    }
    
    if(method == "EMOS"){
      ml_bin <- gamlss(formula = (robs>0) ~ tppow_mean,
                       data = learndata, family = BI())
      
      ml_tr <- crch(formula = robs ~ tppow_mean | log(tppow_sprd + 0.001), 
                    data = learndata_pos, truncated = TRUE, 
                    dist = "gaussian", left = 0, link.scale = "log")
      predpar <- data.frame(predict(ml_tr, type = "location", newdata = testdata),
                            predict(ml_tr, type = "scale", newdata = testdata))
      predpar$nu <- predict(ml_bin, newdata = testdata, what = "mu", type = "response", data = learndata)
    }
    
    colnames(predpar) <- c("mu", "sigma", "nu")
  }
  
  
  
  # loglikelihood
  {
    linkfun <- switch(as.character(distribution), "gaussian" = dist_list_cens_normal$linkfun,
                      "logistic" = dist_list_cens_log$linkfun,
                      "hgaussian" = dist_list_hurdle_normal$linkfun)
    ddist <- switch(as.character(distribution), "gaussian" = dist_list_cens_normal$ddist,
                    "logistic" = dist_list_cens_log$ddist,
                    "hgaussian" = dist_list_hurdle_normal$ddist)
    ll <- numeric(length = NROW(testdata))
    for(j in 1:(NROW(testdata))){
      eta <- as.numeric(linkfun(predpar[j,]))
      if(distribution == "hgaussian"){
        if(!is.na(predpar[j,3]) & predpar[j,3] == 1) {
          nu <- 1-1e-10
          eta[3] <- log(nu/(1-nu))
        }
        if(!is.na(predpar[j,3]) & predpar[j,3] == 0) {
          nu <- 1e-10
          eta[3] <- log(nu/(1-nu))
        }
      }
      ll[j] <- ddist(testdata[j,"robs"], eta = eta, log=TRUE)
    }
  }    
  
  # CPRS
  {
    crps <- numeric(length = NROW(testdata))
    
    if(distribution == "gaussian") crps <- crps_cnorm(testdata$robs, location = predpar$mu, scale = predpar$sigma, lower = 0, upper = Inf)
    if(distribution == "logistic") crps <- crps_clogis(testdata$robs, location = predpar$mu, scale = predpar$sigma, lower = 0, upper = Inf)
    
    if(distribution == "hgaussian"){
      int_upper <- ceiling(max(RainData$robs))  # alternative: quantile(RainData$robs, probs = 0.999)
      for(j in 1:(NROW(testdata))){
        crps[j] <- if(!any(is.na(predpar[j,]))) { 
          integrate(function(z){(1 - predpar$nu[j] + predpar$nu[j] * ptnorm(z, mean = predpar$mu[j], sd = predpar$sigma[j], left = 0) - (testdata$robs[j] <= z))^2},
                    lower = 0, upper = int_upper)$value
        } else NA
      }
    }
  }    
  
  
  results <- data.frame(ll, crps)
  rm(RainData)
  return(results)
  
}






wrapper <- function(stationlist = c("Axams", "Lech", "Zuers", "See im Paznaun", "Jungholz", 
                                    "Ladis-Neuegg", "Oetz", "Ochsengarten-Obergut",
                                    "Ginzling", "Rotholz", "Walchsee", "Koessen", 
                                    "Innervillgraten", "Matrei in Osttirol", 
                                    "St.Johann im Walde"),  
                    methodlist = c("distforest", "gamlss", "gamboostLSS", "EMOS"),
                    distributionlist = c("gaussian", "hgaussian", "logistic"))
{
  set <- expand.grid(station = stationlist, 
                     method = methodlist, 
                     distribution = distributionlist)
  
  nset <- NROW(set)
  ll <- crps <- matrix(nrow = nset, ncol = 123) # 4 months (July) with one missing day
  
  for(i in 1:nset){
    res <- stationeval(station = set$station[i], 
                       method = set$method[i], 
                       distribution = set$distribution[i])
    
    crps[i,] <- res$crps
    ll[i,] <- res$ll
    
  }
  
  means <- cbind(set, rowMeans(crps), rowMeans(ll))
  medians <- cbind(set, apply(crps, 1, median), apply(ll, 1, median))
  colnames(medians) <- colnames(means) <- c(colnames(set), "crps", "ll")
  
  results <- list(crps = crps, 
                  ll = ll, 
                  set = set,
                  medians = medians,
                  means = means)
  
  return(results)
}


results <- wrapper()

save(results, file = "Rain_distributions.rda")


## plots
if(FALSE){
  
  
  ##  map
  data("StationsTyrol", package = "RainTyrol")
  data("MapTyrol", package = "RainTyrol")
  library("sp")
  sp <- SpatialPointsDataFrame(subset(StationsTyrol, select = c(lon, lat)),
                               data = subset(StationsTyrol, select = -c(lon, lat)),
                               proj4string = raster::crs(MapTyrol$RasterLayer))
  
  plot(MapTyrol$SpatialPolygons)
  points(sp, pch = 21, col = "darkgrey", las = 1, cex = 0.6)
  points(sp[c(5, 6, 20, 23, 25, 32, 33, 46, 47, 56, 57, 70, 79, 82, 83),], 
         pch = 21, bg = hcl(325, 100, 70), cex = 1.5)
  # stationname beneath
  text(sp[c(25, 47, 57, 70),], 
       labels = StationsTyrol$name[c(25, 47, 57, 70)],
       pos=1, cex=0.75)
  # stationname left
  text(sp[c(20, 32, 46, 79, 83),], 
       labels = StationsTyrol$name[c(20, 32, 46, 79, 83)],
       pos=2, cex=0.75)
  # stationname above
  text(sp[c(6, 33, 82),], 
       labels = StationsTyrol$name[c(6, 33, 82)],
       pos=3, cex=0.75)
  # stationname right
  text(sp[c(5, 23, 56),], 
       labels = StationsTyrol$name[c(5, 23, 56)],
       pos=4, cex=0.75)
  
  
  ## xy-plot
  
  # HCL palette
  pal <- hcl(c(10, 128, 260, 290, 50), 100, 50)
  
  load("Rain_distributions.rda")
  stations <- c("Axams", "Lech", "Zuers", "See im Paznaun", "Jungholz", 
                "Ladis-Neuegg", "Oetz", "Ochsengarten-Obergut",
                "Ginzling", "Rotholz", "Walchsee", "Koessen", 
                "Innervillgraten", "Matrei in Osttirol", 
                "St.Johann im Walde")
  
  
  means <- results$means[(results$means$method != "disttree") & 
                           results$means$station %in% stations, ]
  rownames(means) <- c(1:NROW(means))
  means$station <- factor(means$station, levels = stations)
  means$method <- factor(means$method, 
                         levels = c("distforest",
                                    "gamlss",
                                    "gamboostLSS",
                                    "EMOS"),
                         labels = c("Distributional forest",
                                    "Prespecified GAMLSS",
                                    "Boosted GAMLSS",
                                    "EMOS"))
  
  strip.background.settings <- trellis.par.get("strip.background")
  strip.background.settings$col <- "gray"
  trellis.par.set("strip.background", strip.background.settings)
  
  strip.border.settings <- trellis.par.get("strip.border")
  strip.border.settings$col <- "black"
  trellis.par.set("strip.border", strip.border.settings)
  
  superpose.line.settings <- trellis.par.get("superpose.line")
  superpose.line.settings$col <- hcl(c(128, 260, 290, 50), 100, 50)
  trellis.par.set("superpose.line", superpose.line.settings)
  
  superpose.symbol.settings <- trellis.par.get("superpose.symbol")
  superpose.symbol.settings$pch <- c(16,24,25,15)
  superpose.symbol.settings$col <- hcl(c(128, 260, 290, 50), 100, 50)
  superpose.symbol.settings$fill <- hcl(c(128, 260, 290, 50), 100, 50)
  trellis.par.set("superpose.symbol", superpose.symbol.settings)
  
  xyplot(crps ~ distribution | station, groups = ~ method, 
         data = means, auto.key = TRUE, 
         type = "o", lwd = 2, lty = 1,  
         drop.unused.levels = TRUE, layout = c(3,5),
         as.table = TRUE, ylab = "CRPS", xlab = "Distribution")
  
  
  #### boxplots
  
  ## CRPS skill score by distribution (reference: gaussian)
  means$crps_ss_dist <- means$crps
  means$crps_ss_dist[means$distribution == "gaussian"] <- 
    1 - means$crps[means$distribution == "gaussian"] / means$crps[means$distribution == "gaussian"]
  means$crps_ss_dist[means$distribution == "logistic"] <- 
    1 - means$crps[means$distribution == "logistic"] / means$crps[means$distribution == "gaussian"]
  means$crps_ss_dist[means$distribution == "hgaussian"] <- 
    1 - means$crps[means$distribution == "hgaussian"] / means$crps[means$distribution == "gaussian"]
  
  means_sel <- means[means$distribution %in% c("logistic", "hgaussian"),]
  means_sel$distribution <- factor(means_sel$distribution, 
                                   levels(means_sel$distribution)[c(2:3)])
  
  # CRPS skill score for distributional forest
  s_df <- matrix(ncol = 2, nrow = 15)
  s_df[,1] <- means_sel[means_sel$method == "Distributional forest" & means_sel$distribution == "hgaussian", "crps_ss_dist"]
  s_df[,2] <- means_sel[means_sel$method == "Distributional forest" & means_sel$distribution == "logistic", "crps_ss_dist"]
  
  # CRPS skill score for prespecified GAMLSS
  s_g <- matrix(ncol = 2, nrow = 15)
  s_g[,1] <- means_sel[means_sel$method == "Prespecified GAMLSS" & means_sel$distribution == "hgaussian", "crps_ss_dist"]
  s_g[,2] <- means_sel[means_sel$method == "Prespecified GAMLSS" & means_sel$distribution == "logistic", "crps_ss_dist"]
  
  # CRPS skill score for boosted GAMLSS
  s_gb <- matrix(ncol = 2, nrow = 15)
  s_gb[,1] <- means_sel[means_sel$method == "Boosted GAMLSS" & means_sel$distribution == "hgaussian", "crps_ss_dist"]
  s_gb[,2] <- means_sel[means_sel$method == "Boosted GAMLSS" & means_sel$distribution == "logistic", "crps_ss_dist"]
  
  # CRPS skill score for EMOS
  s_emos <- matrix(ncol = 2, nrow = 15)
  s_emos[,1] <- means_sel[means_sel$method == "EMOS" & means_sel$distribution == "hgaussian", "crps_ss_dist"]
  s_emos[,2] <- means_sel[means_sel$method == "EMOS" & means_sel$distribution == "logistic", "crps_ss_dist"]
  
  colnames(s_df) <- colnames(s_g) <- colnames(s_gb) <- colnames(s_emos) <- c("hgaussian", "logistic")
  
  # boxplots with matplot
  par(mfrow = c(1,4), mar = c(4,4,3,0.7))
  matplot(t(s_df[,]), type = "l", lwd = 2, 
          col = gray(0.5, alpha = 0.2),
          lty = 1, axes = FALSE,  main = "Distributional forest",
          xlab = "", ylab = "CRPS skill score", xlim = c(0.5, 2.5),
          ylim = c(-0.10, 0.05))
  boxplot(s_df, add = TRUE, col = "transparent")
  abline(h = 0, col = pal[5], lwd = 2)
  
  par(mar = c(4,2.9,3,1.8))
  matplot(t(s_g[,]), type = "l", lwd = 2, 
          col = gray(0.5, alpha = 0.2),
          lty = 1, axes = FALSE, main = "Prespecified GAMLSS",
          xlab = "", ylab = "", xlim = c(0.5, 2.5),
          ylim = c(-0.10, 0.05))
  boxplot(s_g, add = TRUE, col = "transparent", yaxt = 'n')
  abline(h = 0, col = pal[5], lwd = 2)
  
  par(mar = c(4,1.8,3,2.9))
  matplot(t(s_gb[,]), type = "l", lwd = 2, 
          col = gray(0.5, alpha = 0.2),
          lty = 1, axes = FALSE, main = "Boosted GAMLSS",
          xlab = "", ylab = "", xlim = c(0.5, 2.5),
          ylim = c(-0.10, 0.05))
  boxplot(s_gb, add = TRUE, col = "transparent", yaxt = 'n')
  abline(h = 0, col = pal[5], lwd = 2)
  
  par(mar = c(4,0.7,3,4))
  matplot(t(s_emos[,]), type = "l", lwd = 2, 
          col = gray(0.5, alpha = 0.2),
          lty = 1, axes = FALSE,  main = "EMOS",
          xlab = "", ylab = "", xlim = c(0.5, 2.5),
          ylim = c(-0.10, 0.05))
  boxplot(s_emos, add = TRUE, col = "transparent", yaxt = 'n')
  abline(h = 0, col = pal[5], lwd = 2)
  
  
  
  ## CRPS skill score by method (reference: EMOS)
  means$crps_ss_method <- means$crps
  means$crps_ss_method[means$method == "Distributional forest"] <- 
    1 - means$crps[means$method == "Distributional forest"] / means$crps[means$method == "EMOS"]
  means$crps_ss_method[means$method == "Prespecified GAMLSS"] <- 
    1 - means$crps[means$method == "Prespecified GAMLSS"] / means$crps[means$method == "EMOS"]
  means$crps_ss_method[means$method == "Boosted GAMLSS"] <- 
    1 - means$crps[means$method == "Boosted GAMLSS"] / means$crps[means$method == "EMOS"]
  means$crps_ss_method[means$method == "EMOS"] <- 
    1 - means$crps[means$method == "EMOS"] / means$crps[means$method == "EMOS"]
  
  means_sel <- means[means$method %in% c("Distributional forest", "Prespecified GAMLSS", "Boosted GAMLSS"),]
  means_sel$method <- factor(means_sel$method, levels(means_sel$method)[c(1:3)])
  
  # CRPS skill score for gaussian
  s_gaussian <- matrix(ncol = 3, nrow = 15)
  s_gaussian[,1] <- means_sel[means_sel$method == "Distributional forest" & means_sel$distribution == "gaussian", "crps_ss_method"]
  s_gaussian[,2] <- means_sel[means_sel$method == "Prespecified GAMLSS" & means_sel$distribution == "gaussian", "crps_ss_method"]
  s_gaussian[,3] <- means_sel[means_sel$method == "Boosted GAMLSS" & means_sel$distribution == "gaussian", "crps_ss_method"]
  
  # CRPS skill score for boosted hgaussian
  s_hgaussian <- matrix(ncol = 3, nrow = 15)
  s_hgaussian[,1] <- means_sel[means_sel$method == "Distributional forest" & means_sel$distribution == "hgaussian", "crps_ss_method"]
  s_hgaussian[,2] <- means_sel[means_sel$method == "Prespecified GAMLSS" & means_sel$distribution == "hgaussian", "crps_ss_method"]
  s_hgaussian[,3] <- means_sel[means_sel$method == "Boosted GAMLSS" & means_sel$distribution == "hgaussian", "crps_ss_method"]
  
  # CRPS skill score for prespecified logistic
  s_logistic <- matrix(ncol = 3, nrow = 15)
  s_logistic[,1] <- means_sel[means_sel$method == "Distributional forest" & means_sel$distribution == "logistic", "crps_ss_method"]
  s_logistic[,2] <- means_sel[means_sel$method == "Prespecified GAMLSS" & means_sel$distribution == "logistic", "crps_ss_method"]
  s_logistic[,3] <- means_sel[means_sel$method == "Boosted GAMLSS" & means_sel$distribution == "logistic", "crps_ss_method"]
  
  colnames(s_gaussian) <- colnames(s_hgaussian) <- colnames(s_logistic) <- c("Distributional forest", "Prespecified GAMLSS", "Boosted GAMLSS")
  
  
  par(mfrow = c(1,3), mar = c(5,4,3,0))
  matplot(t(s_gaussian[,]), type = "l", lwd = 2, 
          col = gray(0.5, alpha = 0.2),
          lty = 1, axes = FALSE,  main = "gaussian", ylim = c(-0.12,0.23), 
          xlab = "", ylab = "CRPS skill score", xlim = c(0.5, 3.5))
  boxplot(s_gaussian, add = TRUE, col = "transparent", axes = FALSE)
  abline(h = 0, col = pal[5], lwd = 2)
  axis(1, 0:4, c("",
                 "Distr.
                 forest",
                 "Presp.
                 GAMLSS",
                 "Boosted
                 GAMLSS", ""), las=2)
  axis(2, seq(-0.15, 0.25, 0.05), c(seq(-0.15, 0, 0.05), seq(0.05, 0.25, 0.05)), las = 0)
  axis(3, 0:4, lwd.ticks = 0, labels = FALSE)
  axis(4, seq(-0.15, 0.25, 0.05), lwd.ticks = 0, labels = FALSE)
  
  par(mar = c(5,2,3,2))
  matplot(t(s_hgaussian[,]), type = "l", lwd = 2, 
          col = gray(0.5, alpha = 0.2),
          lty = 1, axes = FALSE,  main = "hgaussian", ylim = c(-0.12,0.23), 
          xlab = "", ylab = "", xlim = c(0.5, 3.5))
  boxplot(s_hgaussian, add = TRUE, col = "transparent", axes = FALSE)
  abline(h = 0, col = pal[5], lwd = 2)
  axis(2, seq(-0.15, 0.25, 0.05), lwd.ticks = 0, labels = FALSE)
  axis(1, 0:4, c("",
                 "Distr.
                 forest",
                 "Presp.
                 GAMLSS",
                 "Boosted
                 GAMLSS", ""), las=2)
  axis(3, 0:4, lwd.ticks = 0, labels = FALSE)
  axis(4, seq(-0.15, 0.25, 0.05), lwd.ticks = 0, labels = FALSE)
  
  par(mar = c(5,0,3,4))
  matplot(t(s_logistic[,]), type = "l", lwd = 2, 
          col = gray(0.5, alpha = 0.2),
          lty = 1, axes = FALSE,  main = "logistic", ylim = c(-0.12,0.23), 
          xlab = "", ylab = "", xlim = c(0.5, 3.5))
  boxplot(s_logistic, add = TRUE, col = "transparent", axes = FALSE)
  abline(h = 0, col = pal[5], lwd = 2)
  axis(2, seq(-0.15, 0.25, 0.05), lwd.ticks = 0, labels = FALSE)
  axis(1, 0:4, c("",
                 "Distr.
                 forest",
                 "Presp.
                 GAMLSS",
                 "Boosted
                 GAMLSS", ""), las=2)
  axis(3, 0:4, lwd.ticks = 0, labels = FALSE)
  axis(4, seq(-0.15, 0.25, 0.05), lwd.ticks = 0, labels = FALSE)
  
  
  
  
  ## CRPS 
  
  # CRPS for gaussian
  c_gaussian <- matrix(ncol = 4, nrow = 15)
  c_gaussian[,1] <- means[means$method == "Distributional forest" & means$distribution == "gaussian", "crps"]
  c_gaussian[,2] <- means[means$method == "Prespecified GAMLSS" & means$distribution == "gaussian", "crps"]
  c_gaussian[,3] <- means[means$method == "Boosted GAMLSS" & means$distribution == "gaussian", "crps"]
  c_gaussian[,4] <- means[means$method == "EMOS" & means$distribution == "gaussian", "crps"]
  
  # CRPS for hgaussian
  c_hgaussian <- matrix(ncol = 4, nrow = 15)
  c_hgaussian[,1] <- means[means$method == "Distributional forest" & means$distribution == "hgaussian", "crps"]
  c_hgaussian[,2] <- means[means$method == "Prespecified GAMLSS" & means$distribution == "hgaussian", "crps"]
  c_hgaussian[,3] <- means[means$method == "Boosted GAMLSS" & means$distribution == "hgaussian", "crps"]
  c_hgaussian[,4] <- means[means$method == "EMOS" & means$distribution == "hgaussian", "crps"]
  
  # CRPS for logistic
  c_logistic <- matrix(ncol = 4, nrow = 15)
  c_logistic[,1] <- means[means$method == "Distributional forest" & means$distribution == "logistic", "crps"]
  c_logistic[,2] <- means[means$method == "Prespecified GAMLSS" & means$distribution == "logistic", "crps"]
  c_logistic[,3] <- means[means$method == "Boosted GAMLSS" & means$distribution == "logistic", "crps"]
  c_logistic[,4] <- means[means$method == "EMOS" & means$distribution == "logistic", "crps"]
  
  colnames(c_gaussian) <- colnames(c_hgaussian) <- colnames(c_logistic) <- c("Distributional forest", "Prespecified GAMLSS", "Boosted GAMLSS", "EMOS")
  
  
  par(mfrow = c(1,3), mar = c(5,4,3,0))
  matplot(t(c_gaussian[,]), type = "l", lwd = 2, 
          col = gray(0.5, alpha = 0.2),
          lty = 1, axes = FALSE,  main = "gaussian", 
          xlab = "", ylab = "CRPS", xlim = c(0.5, 4.5),
          ylim = c(0.65,1.15))
  boxplot(c_gaussian, add = TRUE, col = "transparent", axes = FALSE)
  axis(2, seq(0.6, 1.2, 0.1), seq(0.6, 1.2, 0.1), las = 0)
  axis(1, 0:5, c("",
                 "Distr.
                 forest",
                 "Presp.
                 GAMLSS",
                 "Boosted
                 GAMLSS",
                 "EMOS", ""), las=2)
  axis(3, 0:5, lwd.ticks = 0, labels = FALSE)
  axis(4, seq(0.6, 1.2, 0.1), lwd.ticks = 0, labels = FALSE)
  
  par(mar = c(5,2,3,2))
  matplot(t(c_hgaussian[,]), type = "l", lwd = 2, 
          col = gray(0.5, alpha = 0.2),
          lty = 1, axes = FALSE,  main = "hgaussian", 
          xlab = "", ylab = "", xlim = c(0.5, 4.5),
          ylim = c(0.65,1.15))
  boxplot(c_hgaussian, add = TRUE, col = "transparent", axes = FALSE)
  axis(2, seq(0.6, 1.2, 0.1), lwd.ticks = 0, labels = FALSE)
  axis(1, 0:5, c("",
                 "Distr.
                 forest",
                 "Presp.
                 GAMLSS",
                 "Boosted
                 GAMLSS",
                 "EMOS", ""), las=2)
  axis(3, 0:5, lwd.ticks = 0, labels = FALSE)
  axis(4, seq(0.6, 1.2, 0.1), lwd.ticks = 0, labels = FALSE)
  
  par(mar = c(5,0,3,4))
  matplot(t(c_logistic[,]), type = "l", lwd = 2, 
          col = gray(0.5, alpha = 0.2),
          lty = 1, axes = FALSE,  main = "logistic", 
          xlab = "", ylab = "CRPS", xlim = c(0.5, 4.5),
          ylim = c(0.65,1.15))
  boxplot(c_logistic, add = TRUE, col = "transparent", axes = FALSE)
  axis(2, seq(0.6, 1.2, 0.1), lwd.ticks = 0, labels = FALSE)
  axis(1, 0:5, c("",
                 "Distr.
                 forest",
                 "Presp.
                 GAMLSS",
                 "Boosted
                 GAMLSS",
                 "EMOS", ""), las=2)
  axis(3, 0:5, lwd.ticks = 0, labels = FALSE)
  axis(4, seq(0.6, 1.2, 0.1), lwd.ticks = 0, labels = FALSE)

}

