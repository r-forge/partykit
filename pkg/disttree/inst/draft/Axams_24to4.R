#######################################################
### Probabilistic Forecasting on Precipitation Data ###
#######################################################

Axams_24to4 <- function(){
  
  library("disttree")
  
  #####
  # load observations and covariates 
  data("RainAxams")
  
  
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
  
  
  
  #####
  # further packages
  library("gamlss")
  library("gamlss.dist")
  library("gamboostLSS")
  library("crch")
  library("scoringRules")
  
  # if gamlss.cens family object should be used as family
  library("gamlss.cens")
  gen.cens(NO, type = "left")
  
  
  # learning data: 24 years (1985 - 2008, both inlcuded)
  # testing data: 4 successive years (2009, 2010, 2011, 2012)
  learndata <- RainAxams[RainAxams$year < 2009,]
  testdata <- RainAxams[RainAxams$year %in% c(2009, 2010, 2011, 2012),]
  
  
  
  ############
  # fitting the models
  
  set.seed(7)
  
  # fit distributional tree
  dt <- disttree(dt.formula, 
                 data = learndata, family = dist_list_cens_normal, 
                 censtype = "left", censpoint = 0, type.tree = "ctree", 
                 control = ctree_control(teststat = "quad", testtype = "Bonferroni", intersplit = TRUE,
                                         mincriterion = 0.95, minsplit = 50,
                                         minbucket = 20, saveinfo = FALSE))
  
  
  # fit distributional forest
  df <- distforest(df.formula, 
                   data = learndata, family = dist_list_cens_normal, type.tree = "ctree", 
                   ntree = 100, censtype = "left", censpoint = 0, mtry = 27,
                   control = ctree_control(teststat = "quad", testtype = "Univ", intersplit = TRUE,
                                           mincriterion = 0, minsplit = 50,
                                           minbucket = 20, saveinfo = FALSE))
  
  ############
  # fit other heteroscedastic censored gaussian models
  
  # fit prespecified GAM (covariates selected based on meteorological expert knowledge)
  g_learndata <- learndata
  g_learndata$robs <- Surv(g_learndata$robs, g_learndata$robs>0, type="left")
  g <- gamlss(formula = g.mu.formula, sigma.formula = g.sigma.formula, data = g_learndata, 
              family = cens("NO", type = "left"),
              control = gamlss.control(n.cyc = 100),
              i.control = glim.control(cyc = 100, bf.cyc = 100))
  
  # fit boosted GAM
  gb <- gamboostLSS(formula = list(mu = gb.mu.formula, sigma = gb.sigma.formula), data = g_learndata, 
                    families = as.families(fname = cens("NO", type = "left")), method = "noncyclic",
                    control = boost_control(mstop = 1000L))
  # find optimal value for mstop (evalution is very time-consuming)
  grid <- seq(50,1000, by = 25)
  cvr <- cvrisk(gb, grid = grid)
  mstop(gb) <- mstop(cvr)
  
  # fit linear model with only total precipitation as covariate (Ensemble Model Output Statistics, EMOS)
  ml <- crch(formula = robs ~ tppow_mean | log(tppow_sprd + 0.001), 
             data = learndata, dist = "gaussian", left = 0, link.scale = "log")
  
 
  res <- list(dt = dt,
              df = df,
              g = g,
              gb = gb,
              ml = ml,
              learndata = learndata,
              g_learndata = g_learndata,
              testdata = testdata)
   
  return(res)
}