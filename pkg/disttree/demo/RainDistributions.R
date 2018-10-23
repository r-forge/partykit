#######################################################
### Probabilistic Forecasting on Precipitation Data ###
#######################################################

## Replication material for Supplement 1 (Different Response Distributions) of the paper 
## Distributional Regression Forests for Probabilistic Precipitation Forecasting in Complex Terrain (2018)
## by Lisa Schlosser and Torsten Hothorn and Reto Stauffer and Achim Zeileis
## URL: http://arxiv.org/abs/1804.02921

## This demo includes the application on selected stations emplying different distributions
## (models learned on 24 years and evaluated on 4 years)
## Full replication for station Axams can be obtained with
## demo("RainAxams", package = "disttree")
## Full replication of all other results presented in the paper can be obtained with
## demo("RainTyrol", package = "disttree")
## Full replication of Supplement 2 (Stationwise Evaluation) can be obtained with
## demo("RainStationwise", package = "disttree")



## Computation time: approximately 22 hours  
# (on our machines, using 15 kernels)


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

assign("LO",  gamlss.dist::LO,  pos = ".GlobalEnv")
assign("dLO", gamlss.dist::dLO, pos = ".GlobalEnv")
assign("pLO", gamlss.dist::pLO, pos = ".GlobalEnv")
assign("qLO", gamlss.dist::qLO, pos = ".GlobalEnv")
assign("rLO", gamlss.dist::rLO, pos = ".GlobalEnv")
gamlss.cens::gen.cens(LO, type = "left")
assign("LOlc",  LOlc,  pos = ".GlobalEnv")
assign("dLOlc", dLOlc, pos = ".GlobalEnv")
assign("pLOlc", pLOlc, pos = ".GlobalEnv")
assign("qLOlc", qLOlc, pos = ".GlobalEnv")

gamlss.tr::gen.trun(0, "NO", type = "left")
assign("NOtr",  NOtr,  pos = ".GlobalEnv")
assign("dNOtr", dNOtr, pos = ".GlobalEnv")
assign("pNOtr", pNOtr, pos = ".GlobalEnv")
assign("qNOtr", qNOtr, pos = ".GlobalEnv")

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


# evaluation function for one specific combination of station, method and distribution
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





# wrapper function applying stationeval for each a list of stations, methods and distributions
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

library("lattice")
##### 
# HCL palette
pal <- hcl(c(10, 128, 260, 290, 50), 100, 50)

load("Rain_distributions.rda")
levels(results$means$distribution) <- c("cgaussian", "hgaussian", "clogistic")
levels(results$means$station)[15] <- "St. Johann im Walde"


stations <- c("Axams", "Lech", "Zuers", "See im Paznaun", "Jungholz", 
              "Ladis-Neuegg", "Oetz", "Ochsengarten-Obergut",
              "Ginzling", "Rotholz", "Walchsee", "Koessen", 
              "Innervillgraten", "Matrei in Osttirol", 
              "St. Johann im Walde")

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



# crps skill score by distribution (reference: gaussian)
means$crps_ss_dist <- means$crps
means$crps_ss_dist[means$distribution == "cgaussian"] <- 
  1 - means$crps[means$distribution == "cgaussian"] / means$crps[means$distribution == "cgaussian"]
means$crps_ss_dist[means$distribution == "clogistic"] <- 
  1 - means$crps[means$distribution == "clogistic"] / means$crps[means$distribution == "cgaussian"]
means$crps_ss_dist[means$distribution == "hgaussian"] <- 
  1 - means$crps[means$distribution == "hgaussian"] / means$crps[means$distribution == "cgaussian"]

means_sel <- means[means$distribution %in% c("clogistic", "hgaussian"),]
means_sel$distribution <- factor(means_sel$distribution, 
                                 levels(means_sel$distribution)[c(2:3)])

# CRPS skill score for distributional forest
s_df <- matrix(ncol = 2, nrow = 15)
s_df[,1] <- means_sel[means_sel$method == "Distributional forest" & means_sel$distribution == "hgaussian", "crps_ss_dist"]
s_df[,2] <- means_sel[means_sel$method == "Distributional forest" & means_sel$distribution == "clogistic", "crps_ss_dist"]

# CRPS skill score for prespecified GAMLSS
s_g <- matrix(ncol = 2, nrow = 15)
s_g[,1] <- means_sel[means_sel$method == "Prespecified GAMLSS" & means_sel$distribution == "hgaussian", "crps_ss_dist"]
s_g[,2] <- means_sel[means_sel$method == "Prespecified GAMLSS" & means_sel$distribution == "clogistic", "crps_ss_dist"]

# CRPS skill score for boosted GAMLSS
s_gb <- matrix(ncol = 2, nrow = 15)
s_gb[,1] <- means_sel[means_sel$method == "Boosted GAMLSS" & means_sel$distribution == "hgaussian", "crps_ss_dist"]
s_gb[,2] <- means_sel[means_sel$method == "Boosted GAMLSS" & means_sel$distribution == "clogistic", "crps_ss_dist"]

# CRPS skill score for EMOS
s_emos <- matrix(ncol = 2, nrow = 15)
s_emos[,1] <- means_sel[means_sel$method == "EMOS" & means_sel$distribution == "hgaussian", "crps_ss_dist"]
s_emos[,2] <- means_sel[means_sel$method == "EMOS" & means_sel$distribution == "clogistic", "crps_ss_dist"]

colnames(s_df) <- colnames(s_g) <- colnames(s_gb) <- colnames(s_emos) <- c("hgaussian", "clogistic")


# boxplots with matplot
par(mfrow = c(1,4), mar = c(4,4,3,0.7))
matplot(t(s_df[,]), type = "l", lwd = 2, 
        col = gray(0.5, alpha = 0.2),
        lty = 1, axes = FALSE,  main = "Distributional forest",
        xlab = "", ylab = "CRPS skill score", xlim = c(0.5, 2.5),
        ylim = c(-0.05, 0.05))
boxplot(s_df, add = TRUE, col = "transparent")
abline(h = 0, col = pal[5], lwd = 2)
#boxplot(crps_ss_dist~distribution, data = means_sel, 
#        subset = (method == "Distributional forest"), main = "Distributional forest",
#        ylim = c(-0.12,0.07), las = 2, ylab = "CRPS skill score")
#abline(h = 0, col = pal[5], lwd = 2)
par(mar = c(4,2.9,3,1.8))
matplot(t(s_g[,]), type = "l", lwd = 2, 
        col = gray(0.5, alpha = 0.2),
        lty = 1, axes = FALSE, main = "Prespecified GAMLSS",
        xlab = "", ylab = "", xlim = c(0.5, 2.5),
        ylim = c(-0.05, 0.05))
boxplot(s_g, add = TRUE, col = "transparent", yaxt = 'n')
abline(h = 0, col = pal[5], lwd = 2)
#boxplot(crps_ss_dist~distribution, data = means_sel, 
#        subset = (method == "Prespecified GAMLSS"), main = "Prespecified GAMLSS",
#        ylim = c(-0.12,0.07), las = 2)
#abline(h = 0, col = pal[5], lwd = 2)
par(mar = c(4,1.8,3,2.9))
matplot(t(s_gb[,]), type = "l", lwd = 2, 
        col = gray(0.5, alpha = 0.2),
        lty = 1, axes = FALSE, main = "Boosted GAMLSS",
        xlab = "", ylab = "", xlim = c(0.5, 2.5),
        ylim = c(-0.05, 0.05))
boxplot(s_gb, add = TRUE, col = "transparent", yaxt = 'n')
abline(h = 0, col = pal[5], lwd = 2)
#boxplot(crps_ss_dist~distribution, data = means_sel, 
#        subset = (method == "Boosted GAMLSS"), main = "Boosted GAMLSS",
#        ylim = c(-0.12,0.07), las = 2)
#abline(h = 0, col = pal[5], lwd = 2)
par(mar = c(4,0.7,3,4))
matplot(t(s_emos[,]), type = "l", lwd = 2, 
        col = gray(0.5, alpha = 0.2),
        lty = 1, axes = FALSE,  main = "EMOS",
        xlab = "", ylab = "", xlim = c(0.5, 2.5),
        ylim = c(-0.05, 0.05))
boxplot(s_emos, add = TRUE, col = "transparent", yaxt = 'n')
abline(h = 0, col = pal[5], lwd = 2)
#boxplot(crps_ss_dist~distribution, data = means_sel, 
#        subset = (method == "EMOS"), main = "EMOS",
#        ylim = c(-0.12,0.07), las = 2)
#abline(h = 0, col = pal[5], lwd = 2)



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
s_gaussian[,1] <- means_sel[means_sel$method == "Distributional forest" & means_sel$distribution == "cgaussian", "crps_ss_method"]
s_gaussian[,2] <- means_sel[means_sel$method == "Prespecified GAMLSS" & means_sel$distribution == "cgaussian", "crps_ss_method"]
s_gaussian[,3] <- means_sel[means_sel$method == "Boosted GAMLSS" & means_sel$distribution == "cgaussian", "crps_ss_method"]

# CRPS skill score for boosted hgaussian
s_hgaussian <- matrix(ncol = 3, nrow = 15)
s_hgaussian[,1] <- means_sel[means_sel$method == "Distributional forest" & means_sel$distribution == "hgaussian", "crps_ss_method"]
s_hgaussian[,2] <- means_sel[means_sel$method == "Prespecified GAMLSS" & means_sel$distribution == "hgaussian", "crps_ss_method"]
s_hgaussian[,3] <- means_sel[means_sel$method == "Boosted GAMLSS" & means_sel$distribution == "hgaussian", "crps_ss_method"]

# CRPS skill score for prespecified logistic
s_logistic <- matrix(ncol = 3, nrow = 15)
s_logistic[,1] <- means_sel[means_sel$method == "Distributional forest" & means_sel$distribution == "clogistic", "crps_ss_method"]
s_logistic[,2] <- means_sel[means_sel$method == "Prespecified GAMLSS" & means_sel$distribution == "clogistic", "crps_ss_method"]
s_logistic[,3] <- means_sel[means_sel$method == "Boosted GAMLSS" & means_sel$distribution == "clogistic", "crps_ss_method"]

colnames(s_gaussian) <- colnames(s_hgaussian) <- colnames(s_logistic) <- c("Distributional forest", "Prespecified GAMLSS", "Boosted GAMLSS")


par(mfrow = c(1,3), mar = c(5,4,3,0))
matplot(t(s_gaussian[,]), type = "l", lwd = 2, 
        col = gray(0.5, alpha = 0.2),
        lty = 1, axes = FALSE,  main = "cgaussian", ylim = c(-0.12,0.23), 
        xlab = "", ylab = "CRPS skill score", xlim = c(0.5, 3.5))
boxplot(s_gaussian, add = TRUE, col = "transparent", axes = FALSE)
abline(h = 0, col = pal[5], lwd = 2)
#boxplot(crps_ss_method~method, data = means_sel, 
#        subset = (distribution == "cgaussian"), main = "cgaussian",
#        ylim = c(-0.12,0.23), axes = FALSE, ylab = "CRPS skill score")
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
#boxplot(crps_ss_method~method, data = means_sel, 
#        subset = (distribution == "hgaussian"), main = "hgaussian",
#        ylim = c(-0.12,0.23), axes = FALSE)
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
        lty = 1, axes = FALSE,  main = "clogistic", ylim = c(-0.12,0.23), 
        xlab = "", ylab = "", xlim = c(0.5, 3.5))
boxplot(s_logistic, add = TRUE, col = "transparent", axes = FALSE)
abline(h = 0, col = pal[5], lwd = 2)
#boxplot(crps_ss_method~method, data = means_sel, 
#        subset = (distribution == "clogistic"), main = "clogistic",
#        ylim = c(-0.12,0.23), axes = FALSE)
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
c_gaussian[,1] <- means[means$method == "Distributional forest" & means$distribution == "cgaussian", "crps"]
c_gaussian[,2] <- means[means$method == "Prespecified GAMLSS" & means$distribution == "cgaussian", "crps"]
c_gaussian[,3] <- means[means$method == "Boosted GAMLSS" & means$distribution == "cgaussian", "crps"]
c_gaussian[,4] <- means[means$method == "EMOS" & means$distribution == "cgaussian", "crps"]

# CRPS for hgaussian
c_hgaussian <- matrix(ncol = 4, nrow = 15)
c_hgaussian[,1] <- means[means$method == "Distributional forest" & means$distribution == "hgaussian", "crps"]
c_hgaussian[,2] <- means[means$method == "Prespecified GAMLSS" & means$distribution == "hgaussian", "crps"]
c_hgaussian[,3] <- means[means$method == "Boosted GAMLSS" & means$distribution == "hgaussian", "crps"]
c_hgaussian[,4] <- means[means$method == "EMOS" & means$distribution == "hgaussian", "crps"]

# CRPS for logistic
c_logistic <- matrix(ncol = 4, nrow = 15)
c_logistic[,1] <- means[means$method == "Distributional forest" & means$distribution == "clogistic", "crps"]
c_logistic[,2] <- means[means$method == "Prespecified GAMLSS" & means$distribution == "clogistic", "crps"]
c_logistic[,3] <- means[means$method == "Boosted GAMLSS" & means$distribution == "clogistic", "crps"]
c_logistic[,4] <- means[means$method == "EMOS" & means$distribution == "clogistic", "crps"]

colnames(c_gaussian) <- colnames(c_hgaussian) <- colnames(c_logistic) <- c("Distributional forest", "Prespecified GAMLSS", "Boosted GAMLSS", "EMOS")


par(mfrow = c(1,3), mar = c(5,4,3,0))
matplot(t(c_gaussian[,]), type = "l", lwd = 2, 
        col = gray(0.5, alpha = 0.2),
        lty = 1, axes = FALSE,  main = "cgaussian", 
        xlab = "", ylab = "CRPS", xlim = c(0.5, 4.5),
        ylim = c(0.65,1.15))
boxplot(c_gaussian, add = TRUE, col = "transparent", axes = FALSE)
#boxplot(crps~method, data = means, subset = distribution == "cgaussian", main = "cgaussian",
#        ylim = c(0.65,1.15), axes = FALSE, ylab = "CRPS")
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
#boxplot(crps~method, data = means, subset = distribution == "hgaussian", main = "hgaussian",
#        ylim = c(0.65,1.15), las = 2, axes = FALSE)
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
        lty = 1, axes = FALSE,  main = "clogistic", 
        xlab = "", ylab = "CRPS", xlim = c(0.5, 4.5),
        ylim = c(0.65,1.15))
boxplot(c_logistic, add = TRUE, col = "transparent", axes = FALSE)
#boxplot(crps~method, data = means, subset = distribution == "clogistic", main = "clogistic",
#        ylim = c(0.65,1.15), las = 2, axes = FALSE)
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






#######################################
# comparison of two distributional tree/forest hurdle models:
# one-part hurdle vs two part hurdle



#####
# formula for truncated model
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
  
  
}


# select stations
stationlist <- c("Axams", "Lech", "Zuers", "See im Paznaun", "Jungholz", 
                 "Ladis-Neuegg", "Oetz", "Ochsengarten-Obergut",
                 "Ginzling", "Rotholz", "Walchsee", "Koessen", 
                 "Innervillgraten", "Matrei in Osttirol", 
                 "St.Johann im Walde")

results_hurdle <- list()


for(station in stationlist){
  
  #####
  # get observations and covariates for selected station
  RainData <- RainTyrol[RainTyrol$station == as.character(station), ]
  rownames(RainData) <- c(1:NROW(RainData))
  
  # learning data: 24 years (1985 - 2008, both inlcuded)
  # testing data: 4 successive years (2009, 2010, 2011, 2012)
  learndata <- RainData[RainData$year < 2009,]
  learndata_pos <- learndata[learndata$robs>0,]
  testdata <- RainData[RainData$year %in% c(2009, 2010, 2011, 2012),]
  testdata_pos <- testdata[testdata$robs>0,]
  
  
  
  #####
  # fitting the models
  
  set.seed(7)
  
  ## fit distributional tree
  # first a binomial tree to decide whether there is any precipitation at all
  dt_bin <- disttree(dt.formula_bin, 
                     data = learndata,
                     family = dist_binomial, 
                     #censtype = "left", censpoint = 0, 
                     type.tree = "ctree", 
                     control = ctree_control(teststat = "quad", testtype = "Bonferroni", intersplit = TRUE,
                                             mincriterion = 0.95, minsplit = 50,
                                             minbucket = 20))
  
  # fit a tree on the positive observations only (with truncated normal distribution)
  dt_tr <- disttree(dt.formula, 
                    data = learndata_pos,
                    family = dist_list_trunc_normal, 
                    #censtype = "left", censpoint = 0, 
                    type.tree = "ctree", 
                    control = ctree_control(teststat = "quad", testtype = "Bonferroni", intersplit = TRUE,
                                            mincriterion = 0.95, minsplit = 50,
                                            minbucket = 20))
  
  
  ## fit a hurdle tree using a 3-parametric truncated normal distribution with additional
  # third parameter nu modelling robs>0
  dt_h <- disttree(dt.formula, 
                   data = learndata,
                   family = dist_list_hurdle_normal, 
                   #censtype = "left", censpoint = 0, 
                   type.tree = "ctree", 
                   control = ctree_control(teststat = "quad", testtype = "Bonferroni", intersplit = TRUE,
                                           mincriterion = 0.95, minsplit = 50,
                                           minbucket = 20))
  
  
  ## fit distributional forest
  # first a binomial tree to decide whether there is any precipitation at all
  df_bin <- distforest(df.formula_bin, 
                       data = learndata, 
                       family = dist_binomial, type.tree = "ctree", 
                       # censtype = "left", censpoint = 0,
                       ntree = 100, mtry = 27,
                       control = ctree_control(teststat = "quad", testtype = "Univ", intersplit = TRUE,
                                               mincriterion = 0, minsplit = 50,
                                               minbucket = 20))
  
  
  # fit a forest on the positive observations only (with truncated normal distribution)
  df_tr <- distforest(df.formula, 
                      data = learndata_pos, 
                      family = dist_list_trunc_normal, type.tree = "ctree", 
                      # censtype = "left", censpoint = 0,
                      ntree = 100, mtry = 27,
                      control = ctree_control(teststat = "quad", testtype = "Univ", intersplit = TRUE,
                                              mincriterion = 0, minsplit = 50,
                                              minbucket = 20))
  
  ## fit a hurdle forest using a 3-parametric truncated normal distribution with additional
  # third parameter nu modelling robs>0
  df_h <- distforest(df.formula, 
                     data = learndata, 
                     family = dist_list_hurdle_normal, type.tree = "ctree", 
                     # censtype = "left", censpoint = 0,
                     ntree = 100, mtry = 27,
                     control = ctree_control(teststat = "quad", testtype = "Univ", intersplit = TRUE,
                                             mincriterion = 0, minsplit = 50,
                                             minbucket = 20))
  
  
  #####
  # get predicted parameter of all models for testdata
  # distributional tree
  # for hurdle model:
  pdt_h <- predict(dt_h, newdata = testdata, type = "parameter")
  dt_h_mu <- pdt_h$mu
  dt_h_sigma <- pdt_h$sigma
  dt_h_nu <- pdt_h$nu
  dt_h_exp <- (dt_h_mu + dt_h_sigma * (dnorm(dt_h_mu/dt_h_sigma) / pnorm(dt_h_mu/dt_h_sigma))) * dt_h_nu
  
  # for 2 part model
  pdt_2p <- predict(dt_tr, newdata = testdata, type = "parameter")
  dt_2p_mu <- pdt_2p$mu
  dt_2p_sigma <- pdt_2p$sigma
  dt_2p_nu <- predict(dt_bin, newdata = testdata, type = "parameter")[,1]
  dt_2p_exp <- (dt_2p_mu + dt_2p_sigma * (dnorm(dt_2p_mu/dt_2p_sigma) / pnorm(dt_2p_mu/dt_2p_sigma))) * dt_2p_nu
  
  # distributional forest
  # for hurdle model:
  pdf_h <- predict(df_h, newdata = testdata, type = "parameter")
  df_h_mu <- pdf_h$mu
  df_h_sigma <- pdf_h$sigma
  df_h_nu <- pdf_h$nu
  df_h_exp <- (df_h_mu + df_h_sigma * (dnorm(df_h_mu/df_h_sigma) / pnorm(df_h_mu/df_h_sigma))) * df_h_nu
  
  # for 2 part model
  pdf_2p <- predict(df_tr, newdata = testdata, type = "parameter")
  df_2p_mu <- pdf_2p$mu
  df_2p_sigma <- pdf_2p$sigma
  df_2p_nu <- predict(df_bin, newdata = testdata, type = "parameter")[,1]
  df_2p_exp <- (df_2p_mu + df_2p_sigma * (dnorm(df_2p_mu/df_2p_sigma) / pnorm(df_2p_mu/df_2p_sigma))) * df_2p_nu
  
  
  
  
  # CPRS
  {
    crps_dt_h <- crps_dt_2p <- crps_df_h <- crps_df_2p <- numeric(length = NROW(testdata))
    crps_dt_h_error <- logical(length = NROW(testdata))
    
    int_upper <- ceiling(max(RainData$robs))  # alternative: quantile(RainData$robs, probs = 0.999)
    
    for(j in 1:(NROW(testdata))){
      crps_dt_h_try <- try(integrate(function(z){(1 - dt_h_nu[j] + dt_h_nu[j] * ptnorm(z, mean = dt_h_mu[j], sd = dt_h_sigma[j], left = 0) - (testdata$robs[j] <= z))^2},
                                 lower = 0, upper = int_upper)$value, silent = TRUE)
      if(inherits(crps_dt_h_try, "try-error")) {
        crps_dt_h[j] <- NA 
        crps_dt_h_error[j] <- TRUE
      } else {
        crps_dt_h[j] <- crps_dt_h_try
        crps_dt_h_error[j] <- FALSE
      }
      crps_dt_2p[j] <- integrate(function(z){(1 - dt_2p_nu[j] + dt_2p_nu[j] * ptnorm(z, mean = dt_2p_mu[j], sd = dt_2p_sigma[j], left = 0) - (testdata$robs[j] <= z))^2},
                                 lower = 0, upper = int_upper)$value
      crps_df_h[j] <- integrate(function(z){(1 - df_h_nu[j] + df_h_nu[j] * ptnorm(z, mean = df_h_mu[j], sd = df_h_sigma[j], left = 0) - (testdata$robs[j] <= z))^2},
                                lower = 0, upper = int_upper)$value
      crps_df_2p[j] <- integrate(function(z){(1 - df_2p_nu[j] + df_2p_nu[j] * ptnorm(z, mean = df_2p_mu[j], sd = df_2p_sigma[j], left = 0) - (testdata$robs[j] <= z))^2},
                                 lower = 0, upper = int_upper)$value
    }
    
    crps_dt_h <- if(all(crps_dt_h_error)) NA else mean(crps_dt_h, na.rm = TRUE)
    crps_dt_2p <- mean(crps_dt_2p, na.rm = TRUE)
    crps_df_h <- mean(crps_df_h, na.rm = TRUE)
    crps_df_2p <- mean(crps_df_2p, na.rm = TRUE)
  }
  
  
  ## error: expected value - observation
  expobs_dt_h <- abs(dt_h_exp - testdata[,"robs"])
  expobs_dt_2p <- abs(dt_2p_exp - testdata[,"robs"])
  expobs_df_h <- abs(df_h_exp - testdata[,"robs"])
  expobs_df_2p <- abs(df_2p_exp - testdata[,"robs"])
  
  # RMSE
  rmse_dt_h <- sqrt(mean(expobs_dt_h^2))
  rmse_dt_2p <- sqrt(mean(expobs_dt_2p^2))
  rmse_df_h <- sqrt(mean(expobs_df_h^2))
  rmse_df_2p <- sqrt(mean(expobs_df_2p^2))
  
  # loglikelihood
  {
    dthll <- dt2pll <- dfhll <- df2pll <- numeric(length = NROW(testdata))
    for(j in 1:(NROW(testdata))){
      
      eta_dt_h <- as.numeric(dist_list_hurdle_normal$linkfun(cbind(dt_h_mu, dt_h_sigma, dt_h_nu)[j,]))
      eta_dt_2p <- as.numeric(dist_list_hurdle_normal$linkfun(cbind(dt_2p_mu, dt_2p_sigma, dt_2p_nu)[j,]))
      eta_df_h <- as.numeric(dist_list_hurdle_normal$linkfun(cbind(df_h_mu, df_h_sigma, df_h_nu)[j,]))
      eta_df_2p <- as.numeric(dist_list_hurdle_normal$linkfun(cbind(df_2p_mu, df_2p_sigma, df_2p_nu)[j,]))
      
      dthll[j] <- dist_list_hurdle_normal$ddist(testdata[j,"robs"], eta = eta_dt_h, log=TRUE)
      dt2pll[j] <- dist_list_hurdle_normal$ddist(testdata[j,"robs"], eta = eta_dt_2p, log=TRUE)
      dfhll[j] <- dist_list_hurdle_normal$ddist(testdata[j,"robs"], eta = eta_df_h, log=TRUE)
      df2pll[j] <- dist_list_hurdle_normal$ddist(testdata[j,"robs"], eta = eta_df_2p, log=TRUE)
    }
    
    dthll <- mean(dthll, na.rm = TRUE)
    dt2pll <- mean(dt2pll, na.rm = TRUE)
    dfhll <- mean(dfhll, na.rm = TRUE)
    df2pll <- mean(df2pll, na.rm = TRUE)
  }
  
  
  # or 
  # loglikelihood separated
  #{
  #  dthll <- dt2pll <- dfhll <- df2pll <- gll <-  gbll <- mlll <- numeric(length = NROW(testdata))
  #  for(j in 1:(NROW(testdata))){
  #    
  #    eta_dt_h <- as.numeric(dist_list_trunc_normal$linkfun(cbind(dt_h_mu, dt_h_sigma)[j,]))
  #    eta_dt_2p <- as.numeric(dist_list_trunc_normal$linkfun(cbind(dt_2p_mu, dt_2p_sigma)[j,]))
  #    eta_df_h <- as.numeric(dist_list_trunc_normal$linkfun(cbind(df_h_mu, df_h_sigma)[j,]))
  #    eta_df_2p <- as.numeric(dist_list_trunc_normal$linkfun(cbind(df_2p_mu, df_2p_sigma)[j,]))
  #    eta_g  <- as.numeric(dist_list_trunc_normal$linkfun(cbind(g_mu, g_sigma)[j,])) 
  #    eta_gb <- as.numeric(dist_list_trunc_normal$linkfun(cbind(gb_mu, gb_sigma)[j,]))  
  #    eta_ml <- as.numeric(dist_list_trunc_normal$linkfun(cbind(ml_mu, ml_sigma)[j,]))
  #    
  #    dist_list_binomial <- dist_binomial()
  #    etanu_dt_h <- as.numeric(dist_list_binomial$linkfun(cbind(dt_h_nu)[j,]))
  #    etanu_dt_2p <- as.numeric(dist_list_binomial$linkfun(cbind(dt_2p_nu)[j,]))
  #    etanu_df_h <- as.numeric(dist_list_binomial$linkfun(cbind(df_h_nu)[j,]))
  #    etanu_df_2p <- as.numeric(dist_list_binomial$linkfun(cbind(df_2p_nu)[j,]))
  #    etanu_g  <- as.numeric(dist_list_binomial$linkfun(cbind(g_nu)[j,])) 
  #    etanu_gb <- as.numeric(dist_list_binomial$linkfun(cbind(gb_nu)[j,]))  
  #    etanu_ml <- as.numeric(dist_list_binomial$linkfun(cbind(ml_nu)[j,]))
  #    
  #    dthll[j] <- dist_list_trunc_normal$ddist(testdata[j,"robs"], eta = eta_dt_h, log=TRUE) * (testdata[j,"robs"]>0) + 
  #      dist_list_binomial$ddist(testdata[j,"robs"]>0, eta = etanu_dt_h, log=TRUE)
  #    dt2pll[j] <- dist_list_trunc_normal$ddist(testdata[j,"robs"], eta = eta_dt_2p, log=TRUE) * (testdata[j,"robs"]>0) + 
  #      dist_list_binomial$ddist(testdata[j,"robs"]>0, eta = etanu_dt_2p, log=TRUE) 
  #    dfhll[j] <- dist_list_trunc_normal$ddist(testdata[j,"robs"], eta = eta_df_h, log=TRUE) * (testdata[j,"robs"]>0) + 
  #      dist_list_binomial$ddist(testdata[j,"robs"]>0, eta = etanu_df_h, log=TRUE)
  #    df2pll[j] <- dist_list_trunc_normal$ddist(testdata[j,"robs"], eta = eta_df_2p, log=TRUE) * (testdata[j,"robs"]>0) + 
  #      dist_list_binomial$ddist(testdata[j,"robs"]>0, eta = etanu_df_2p, log=TRUE)
  #    gll[j]  <- dist_list_trunc_normal$ddist(testdata[j,"robs"], eta = eta_g, log=TRUE) * (testdata[j,"robs"]>0) + 
  #      dist_list_binomial$ddist(testdata[j,"robs"]>0, eta = etanu_g, log=TRUE)
  #    gbll[j] <- dist_list_trunc_normal$ddist(testdata[j,"robs"], eta = eta_gb, log=TRUE) * (testdata[j,"robs"]>0) + 
  #      dist_list_binomial$ddist(testdata[j,"robs"]>0, eta = etanu_gb, log=TRUE)
  #    mlll[j] <- dist_list_trunc_normal$ddist(testdata[j,"robs"], eta = eta_ml, log=TRUE) * (testdata[j,"robs"]>0) + 
  #      dist_list_binomial$ddist(testdata[j,"robs"]>0, eta = etanu_ml, log=TRUE)
  
  #  }
  
  #  dthll <- mean(dthll, na.rm = TRUE)
  #  dt2pll <- mean(dt2pll, na.rm = TRUE)
  #  dfhll <- mean(dfhll, na.rm = TRUE)
  #  df2pll <- mean(df2pll, na.rm = TRUE)
  #  gll  <- mean(gll, na.rm = TRUE) 
  #  gbll <- mean(gbll, na.rm = TRUE)
  #  mlll <- mean(mlll, na.rm = TRUE)
  #}
  
  ll <- c(dthll, dt2pll, dfhll, df2pll) 
  rmse <- c(rmse_dt_h, rmse_dt_2p, rmse_df_h, rmse_df_2p)
  crps <- c(crps_dt_h, crps_dt_2p, crps_df_h, crps_df_2p)
  
  # compare results
  results_hurdle[[station]] <- rbind(ll, rmse, crps)
  colnames(results_hurdle[[station]]) <- c("tree_h", "tree_2p", "forest_h", "forest_2p")

  
}


save(results_hurdle, file = "Rain_hurdleversions_forest.rda")

load("Rain_hurdleversions_forest.rda")
crps_forests <- matrix(ncol = 2, nrow = length(names(results_hurdle)))
rownames(crps_forests) <- names(results_hurdle)
colnames(crps_forests) <- c("forest_h", "forest_2p")
for(i in 1:length(names(results_hurdle))){
  crps_forests[i,] <- results_hurdle[[i]]["crps", c("forest_h", "forest_2p")]
}

save(crps_forests, file = "crps_forests.rda")



## plot
load("crps_forests.rda")
crps_ss_forest_2p <- data.frame(1 - crps_forests[,"forest_2p"] / crps_forests[,"forest_h"])
names(crps_ss_forest_2p) <- "Three-parametric forest"
par(mfrow = c(1,1), mar = c(5,4,2,2))
boxplot(crps_ss_forest_2p, ylab = "CRPS skill score", names = c("Three-parametric forest"))
abline(h = 0, col = pal[5], lwd = 2)
axis(1, 0:2, c("","One-part vs. two-part hurdle forest", ""), las=1)

