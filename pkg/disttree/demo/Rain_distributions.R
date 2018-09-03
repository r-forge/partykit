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


load("~/svn/partykit/pkg/RainTyrol/data/RainTyrol.rda")


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


  # specify distribution family
  
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
        warning("Error in gamlss, repeated with 10 itertations only (instead of 100)")
        g <- try(gamlss(formula = g.mu.formula, sigma.formula = g.sigma.formula, data = g_learndata, 
                        family = if(distribution == "gaussian") cens("NO", type = "left") else cens("LO", type = "left"),
                        control = gamlss.control(n.cyc = 10),
                        i.control = glim.control(cyc = 10, bf.cyc = 10)),
                 silent = TRUE)
      }
      if(inherits(g, "try-error")){
        warning("Error in gamlss, repeated with 5 itertations only (instead of 100)")
        g <- gamlss(formula = g.mu.formula, sigma.formula = g.sigma.formula, data = g_learndata, 
                    family = if(distribution == "gaussian") cens("NO", type = "left") else cens("LO", type = "left"),
                    control = gamlss.control(n.cyc = 5),
                    i.control = glim.control(cyc = 5, bf.cyc = 5))
      }
      predpar <- data.frame(predict(g, newdata = testdata, what = "mu", type = "response", data = g_learndata),
                            predict(g, newdata = testdata, what = "sigma", type = "response", data = g_learndata))
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
        warning("Error in gamlss for binomial model, repeated with 10 itertations only (instead of 100)")
        g_bin <- try(gamlss(formula = g.mu.formula_bin, 
                            data = learndata, 
                            family = BI(),
                            control = gamlss.control(n.cyc = 10),
                            i.control = glim.control(cyc = 10, bf.cyc = 10)), silent = TRUE)
      }
      if(inherits(g_bin, "try-error")){
        warning("Error in gamlss for binomial model, repeated with 5 itertations only (instead of 100)")
        g_bin <- gamlss(formula = g.mu.formula_bin, 
                            data = learndata, 
                            family = BI(),
                            control = gamlss.control(n.cyc = 5),
                            i.control = glim.control(cyc = 5, bf.cyc = 5))
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
        warning("Error in gamlss for truncated model, repeated with 10 itertations only (instead of 100)")
        g_tr <- try(gamlss(formula = g.mu.formula, sigma.formula = g.sigma.formula, 
                           data = learndata_pos, 
                           family = trun(0, "NO", type = "left"),
                           control = gamlss.control(n.cyc = 10),
                           i.control = glim.control(cyc = 10, bf.cyc = 10)), silent = TRUE)
      }
      if(inherits(g_tr, "try-error")){
        warning("Error in gamlss for truncated model, repeated with 10 itertations only (instead of 100)")
        g_tr <- gamlss(formula = g.mu.formula, sigma.formula = g.sigma.formula, 
                           data = learndata_pos, 
                           family = trun(0, "NO", type = "left"),
                           control = gamlss.control(n.cyc = 5),
                           i.control = glim.control(cyc = 5, bf.cyc = 5))
      }
      
      predpar <- data.frame(predict(g_tr, newdata = testdata, what = "mu", type = "response", data = learndata_pos),
                       predict(g_tr, newdata = testdata, what = "sigma", type = "response", data = learndata_pos))
      predpar$nu <- as.numeric(predict(g_bin, newdata = testdata, what = "mu", type = "response", data = learndata))
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
                      "hgaussian" = dist_list_hgaussian_normal$linkfun)
    ddist <- switch(as.character(distribution), "gaussian" = dist_list_cens_normal$ddist,
                    "logistic" = dist_list_cens_log$ddist,
                    "hgaussian" = dist_list_hgaussian_normal$ddist)
    ll <- numeric(length = NROW(testdata))
    for(j in 1:(NROW(testdata))){
      eta <- as.numeric(linkfun(predpar[j,]))
      if(distribution == "hgaussian"){
        if(predpar[j,3] == 1) {
          nu <- 1-1e-10
          eta[3] <- log(nu/(1-nu))
        }
        if(predpar[j,3] == 0) {
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
        crps[j] <- integrate(function(z){(1 - predpar$nu[j] + predpar$nu[j] * ptnorm(z, mean = predpar$mu[j], sd = predpar$sigma[j], left = 0) - (testdata$robs[j] <= z))^2},
                             lower = 0, upper = int_upper)$value
      }
    }
  }    
  
  
  results <- data.frame(ll, crps)
  rm(RainData)
  return(results)

}
  





wrapper <- function(stationlist = c("Axams", "Lech", "Walchsee", "Vils", "Oetz", "Zuers", "Innervillgraten", "Pass Thurn",
                                    "Koessen", "See im Paznaun", "Matrei in Osttirol", "Jungholz", "Rotholz",
                                    "Ochsengarten-Obergut", "Ladis-Neuegg", "St.Johann im Walde"),        
                                    # "Laengenfeld"
                    methodlist = c("disttree", "distforest", "gamlss", "gamboostLSS", "EMOS"),
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

save(results, file = "~/svn/partykit/pkg/disttree/demo/Rain_distributions.rda")


## plots
if(FALSE){
  library("lattice")
  #xyplot(crps ~ method | station, groups = ~ distribution, data = results$means, auto.key = TRUE)
  
  strip.background.settings <- trellis.par.get("strip.background")
  str(strip.background.settings)
  strip.background.settings$col <- "gray"
  str(strip.background.settings)
  trellis.par.set("strip.background", strip.background.settings)
  
  strip.border.settings <- trellis.par.get("strip.border")
  str(strip.border.settings)
  strip.border.settings$col <- "black"
  str(strip.border.settings)
  trellis.par.set("strip.border", strip.border.settings)
  
  superpose.line.settings <- trellis.par.get("superpose.line")
  str(superpose.line.settings)
  superpose.line.settings$col <- hcl(c(10, 128, 260, 290, 50), 100, 50)
  str(superpose.line.settings)
  trellis.par.set("superpose.line", superpose.line.settings)
  
  superpose.symbol.settings <- trellis.par.get("superpose.symbol")
  str(superpose.symbol.settings)
  superpose.symbol.settings$pch <- c(3,16,24,25,15)
  superpose.symbol.settings$col <- hcl(c(10, 128, 260, 290, 50), 100, 50)
  superpose.symbol.settings$fill <- hcl(c(10, 128, 260, 290, 50), 100, 50)
  str(superpose.symbol.settings)
  trellis.par.set("superpose.symbol", superpose.symbol.settings)
  
  xyplot(crps ~ distribution | station, groups = ~ method, data = results$means, auto.key = TRUE, type = "o", lwd = 2, lty = 1)
  
  xyplot(crps ~ distribution | station, groups = ~ method, data = results$means, 
         subset = station %in% c("Rotholz", "Innervillgraten", "Pass Thurn",
                                 "See im Paznaun", "Axams", "Lech", "Walchsee", "Vils"),
         auto.key = TRUE, type = "o", lwd = 2, lty = 1)
  
  
  ##########
  # boxplots
  
  # crps
  par(mfrow = c(1,3))
  boxplot(crps~method, data = results$means, subset = distribution == "gaussian", main = "gaussian",
          ylim = c(0.65,1.15))
  boxplot(crps~method, data = results$means, subset = distribution == "logistic", main = "logistic",
          ylim = c(0.65,1.15))
  boxplot(crps~method, data = results$means, subset = distribution == "hgaussian", main = "hgaussian",
          ylim = c(0.65,1.15))
  
  
  # crps skill score by method (reference: EMOS)
  means <- results$means
  means$crps_ss_method <- means$crps
  means$crps_ss_method[means$method == "disttree"] <- 
    1 - means$crps[means$method == "disttree"] / means$crps[means$method == "EMOS"]
  means$crps_ss_method[means$method == "distforest"] <- 
    1 - means$crps[means$method == "distforest"] / means$crps[means$method == "EMOS"]
  means$crps_ss_method[means$method == "gamlss"] <- 
    1 - means$crps[means$method == "gamlss"] / means$crps[means$method == "EMOS"]
  means$crps_ss_method[means$method == "gamboostLSS"] <- 
    1 - means$crps[means$method == "gamboostLSS"] / means$crps[means$method == "EMOS"]
  means$crps_ss_method[means$method == "EMOS"] <- 
    1 - means$crps[means$method == "EMOS"] / means$crps[means$method == "EMOS"]
  
  means_sel <- means[means$method %in% c("distforest", "gamlss", "gamboostLSS"),]
  means_sel$method <- factor(means_sel$method, levels(means_sel$method)[c(2:4)])
  
  par(mfrow = c(1,3))
  boxplot(crps_ss_method~method, data = means_sel, 
          subset = (distribution == "gaussian"), main = "gaussian",
          ylim = c(-0.12,0.23))
  abline(h = 0, col = pal[5], lwd = 2)
  boxplot(crps_ss_method~method, data = means_sel, 
          subset = (distribution == "logistic"), main = "logistic",
          ylim = c(-0.12,0.23))
  abline(h = 0, col = pal[5], lwd = 2)
  boxplot(crps_ss_method~method, data = means_sel, 
          subset = (distribution == "hgaussian"), main = "hgaussian",
          ylim = c(-0.12,0.23))
  abline(h = 0, col = pal[5], lwd = 2)
  
  
  
  # crps skill score by distribution (reference: gaussian)
  means$crps_ss_dist <- means$crps
  means$crps_ss_dist[means$distribution == "gaussian"] <- 
    1 - means$crps[means$distribution == "gaussian"] / means$crps[means$distribution == "gaussian"]
  means$crps_ss_dist[means$distribution == "logistic"] <- 
    1 - means$crps[means$distribution == "logistic"] / means$crps[means$distribution == "gaussian"]
  means$crps_ss_dist[means$distribution == "hgaussian"] <- 
    1 - means$crps[means$distribution == "hgaussian"] / means$crps[means$distribution == "gaussian"]
  
  means_sel <- means[means$distribution %in% c("logistic", "hgaussian"),]
  means_sel$distribution <- factor(means_sel$distribution, levels(means_sel$distribution)[c(2:3)])
  
  par(mfrow = c(1,5))
  boxplot(crps_ss_dist~distribution, data = means_sel, 
          subset = (method == "disttree"), main = "disttree",
          ylim = c(-0.12,0.12))
  abline(h = 0, col = pal[5], lwd = 2)
  boxplot(crps_ss_dist~distribution, data = means_sel, 
          subset = (method == "distforest"), main = "distforest",
          ylim = c(-0.12,0.12))
  abline(h = 0, col = pal[5], lwd = 2)
  boxplot(crps_ss_dist~distribution, data = means_sel, 
          subset = (method == "gamlss"), main = "gamlss",
          ylim = c(-0.12,0.12))
  abline(h = 0, col = pal[5], lwd = 2)
  boxplot(crps_ss_dist~distribution, data = means_sel, 
          subset = (method == "gamboostLSS"), main = "gamboostLSS",
          ylim = c(-0.12,0.12))
  abline(h = 0, col = pal[5], lwd = 2)
  boxplot(crps_ss_dist~distribution, data = means_sel, 
          subset = (method == "EMOS"), main = "EMOS",
          ylim = c(-0.12,0.12))
  abline(h = 0, col = pal[5], lwd = 2)
}



{
  
  summary(results$crps[results$set$method == "gamboostLSS"])
  summary(results$crps[results$set$method == "gamlss"])
  summary(results$crps[results$set$method == "distforest"])
  
  
  library("lattice")
  #xyplot(crps ~ method | station, groups = ~ distribution, data = results$medians, auto.key = TRUE)
  
  strip.background.settings <- trellis.par.get("strip.background")
  str(strip.background.settings)
  strip.background.settings$col <- "gray"
  str(strip.background.settings)
  trellis.par.set("strip.background", strip.background.settings)
  
  strip.border.settings <- trellis.par.get("strip.border")
  str(strip.border.settings)
  strip.border.settings$col <- "black"
  str(strip.border.settings)
  trellis.par.set("strip.border", strip.border.settings)
  
  superpose.line.settings <- trellis.par.get("superpose.line")
  str(superpose.line.settings)
  superpose.line.settings$col <- hcl(c(10, 128, 260, 290, 50), 100, 50)
  str(superpose.line.settings)
  trellis.par.set("superpose.line", superpose.line.settings)
  
  superpose.symbol.settings <- trellis.par.get("superpose.symbol")
  str(superpose.symbol.settings)
  superpose.symbol.settings$pch <- c(3,16,24,25,15)
  superpose.symbol.settings$col <- hcl(c(10, 128, 260, 290, 50), 100, 50)
  superpose.symbol.settings$fill <- hcl(c(10, 128, 260, 290, 50), 100, 50)
  str(superpose.symbol.settings)
  trellis.par.set("superpose.symbol", superpose.symbol.settings)
  
  xyplot(crps ~ distribution | station, groups = ~ method, data = results$medians, auto.key = TRUE, type = "o", lwd = 2, lty = 1)
  
  
  
  ##########
  # boxplots
  
  # crps
  par(mfrow = c(1,3))
  boxplot(crps~method, data = results$medians, subset = distribution == "gaussian", main = "gaussian",
          ylim = c(0.4,1))
  boxplot(crps~method, data = results$medians, subset = distribution == "logistic", main = "logistic",
          ylim = c(0.4,1))
  boxplot(crps~method, data = results$medians, subset = distribution == "hgaussian", main = "hgaussian",
          ylim = c(0.4,1))
  
  
  # crps skill score by method (reference: EMOS)
  medians <- results$medians
  medians$crps_ss_method <- medians$crps
  medians$crps_ss_method[medians$method == "disttree"] <- 
    1 - medians$crps[medians$method == "disttree"] / medians$crps[medians$method == "EMOS"]
  medians$crps_ss_method[medians$method == "distforest"] <- 
    1 - medians$crps[medians$method == "distforest"] / medians$crps[medians$method == "EMOS"]
  medians$crps_ss_method[medians$method == "gamlss"] <- 
    1 - medians$crps[medians$method == "gamlss"] / medians$crps[medians$method == "EMOS"]
  medians$crps_ss_method[medians$method == "gamboostLSS"] <- 
    1 - medians$crps[medians$method == "gamboostLSS"] / medians$crps[medians$method == "EMOS"]
  medians$crps_ss_method[medians$method == "EMOS"] <- 
    1 - medians$crps[medians$method == "EMOS"] / medians$crps[medians$method == "EMOS"]
  
  medians_sel <- medians[medians$method %in% c("distforest", "gamlss", "gamboostLSS"),]
  medians_sel$method <- factor(medians_sel$method, levels(medians_sel$method)[c(2:4)])
  
  par(mfrow = c(1,3))
  boxplot(crps_ss_method~method, data = medians_sel, 
          subset = (distribution == "gaussian"), main = "gaussian",
          ylim = c(-0.3,0.3))
  abline(h = 0, col = pal[5], lwd = 2)
  boxplot(crps_ss_method~method, data = medians_sel, 
          subset = (distribution == "logistic"), main = "logistic",
          ylim = c(-0.3,0.3))
  abline(h = 0, col = pal[5], lwd = 2)
  boxplot(crps_ss_method~method, data = medians_sel, 
          subset = (distribution == "hgaussian"), main = "hgaussian",
          ylim = c(-0.3,0.3))
  abline(h = 0, col = pal[5], lwd = 2)
  
  
  
  # crps skill score by distribution (reference: gaussian)
  medians$crps_ss_dist <- medians$crps
  medians$crps_ss_dist[medians$distribution == "gaussian"] <- 
    1 - medians$crps[medians$distribution == "gaussian"] / medians$crps[medians$distribution == "gaussian"]
  medians$crps_ss_dist[medians$distribution == "logistic"] <- 
    1 - medians$crps[medians$distribution == "logistic"] / medians$crps[medians$distribution == "gaussian"]
  medians$crps_ss_dist[medians$distribution == "hgaussian"] <- 
    1 - medians$crps[medians$distribution == "hgaussian"] / medians$crps[medians$distribution == "gaussian"]
  
  medians_sel <- medians[medians$distribution %in% c("logistic", "hgaussian"),]
  medians_sel$distribution <- factor(medians_sel$distribution, levels(medians_sel$distribution)[c(2:3)])
  
  par(mfrow = c(1,4))
  boxplot(crps_ss_dist~distribution, data = medians_sel, 
          subset = (method == "disttree"), main = "distforest",
          ylim = c(-0.4,0.2))
  abline(h = 0, col = pal[5], lwd = 2)
  boxplot(crps_ss_dist~distribution, data = medians_sel, 
          subset = (method == "distforest"), main = "distforest",
          ylim = c(-0.4,0.2))
  abline(h = 0, col = pal[5], lwd = 2)
  boxplot(crps_ss_dist~distribution, data = medians_sel, 
          subset = (method == "gamlss"), main = "gamlss",
          ylim = c(-0.4,0.2))
  abline(h = 0, col = pal[5], lwd = 2)
  boxplot(crps_ss_dist~distribution, data = medians_sel, 
          subset = (method == "gamboostLSS"), main = "gamboostLSS",
          ylim = c(-0.4,0.2))
  abline(h = 0, col = pal[5], lwd = 2)
}
