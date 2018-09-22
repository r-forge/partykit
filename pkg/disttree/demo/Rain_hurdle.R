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
# HCL palette
pal <- hcl(c(10, 128, 260, 290, 50), 100, 50)


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


#####
# further packages
library("gamlss")
library("gamlss.dist")
library("gamboostLSS")
library("crch")
library("scoringRules")

# if gamlss.tr family object should be used as family
library("gamlss.tr")
gen.trun(0, "NO", type = "left")



# select stations
stationlist <- c("Axams", "Lech", "Walchsee", "Oetz", "Koessen", "Vils", 
                 "See", "Jungholz", "Zuers")
results <- list()

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
  # fit other heteroscedastiCc censored gaussian models
  
  # fit prespecified GAM (covariates selected based on meteorological expert knowledge)
  #g_learndata <- learndata
  #g_learndata$robs <- Surv(g_learndata$robs, g_learndata$robs>0, type="left")
  g_bin <- gamlss(formula = g.mu.formula_bin, 
                  data = learndata, 
                  family = BI(),
                  control = gamlss.control(n.cyc = 100),
                  i.control = glim.control(cyc = 100, bf.cyc = 100))
  
  g_tr <- try(gamlss(formula = g.mu.formula, sigma.formula = g.sigma.formula, 
                     data = learndata_pos, 
                     family = trun(0, "NO", type = "left"),
                     control = gamlss.control(n.cyc = 100),
                     i.control = glim.control(cyc = 100, bf.cyc = 100)), silent = TRUE)
  if(inherits(g_tr, "try-error")){
    warning("Error in gamlss, nr of iterations reduced to 50 (from 100)")
    g_tr <- try(gamlss(formula = g.mu.formula, sigma.formula = g.sigma.formula, 
                       data = learndata_pos, 
                       family = trun(0, "NO", type = "left"),
                       control = gamlss.control(n.cyc = 50),
                       i.control = glim.control(cyc = 50, bf.cyc = 50)), silent = TRUE)
  }  
  if(inherits(g_tr, "try-error")){
    warning("Error in gamlss, nr of iterations reduced to 10 (from 50)")
    g_tr <- gamlss(formula = g.mu.formula, sigma.formula = g.sigma.formula, 
                   data = learndata_pos, 
                   family = trun(0, "NO", type = "left"),
                   control = gamlss.control(n.cyc = 10),
                   i.control = glim.control(cyc = 10, bf.cyc = 10))
  }  
  # fit boosted GAM
  library("mboost")
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
  
  
  
  # fit linear model with only total precipitation as covariate (Ensemble Model Output Statistics, EMOS)
  ml_bin <- gamlss(formula = (robs>0) ~ tppow_mean,
                   data = learndata, family = BI())
  
  ml_tr <- crch(formula = robs ~ tppow_mean | log(tppow_sprd + 0.001), 
                data = learndata_pos, truncated = TRUE, 
                dist = "gaussian", left = 0, link.scale = "log")
  
  
  
  
  
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
  
  
  # prespecified GAM
  g_mu <- predict(g_tr, newdata = testdata, what = "mu", type = "response", data = learndata_pos)
  g_sigma <- predict(g_tr, newdata = testdata, what = "sigma", type = "response", data = learndata_pos)
  g_nu <- predict(g_bin, newdata = testdata, what = "mu", type = "response", data = learndata)
  g_exp <- (g_mu + g_sigma * (dnorm(g_mu/g_sigma) / pnorm(g_mu/g_sigma))) * g_nu
  
  # boosted GAM
  pgb <- predict(gb_tr, newdata = testdata, parameter = c("mu","sigma"), type = "response")
  gb_mu <- pgb$mu
  gb_sigma <- pgb$sigma
  gb_nu <- predict(gb_bin, newdata = testdata, parameter = c("mu"), type = "response")
  gb_exp <- (gb_mu + gb_sigma * (dnorm(gb_mu/gb_sigma) / pnorm(gb_mu/gb_sigma))) * gb_nu
  
  # EMOS
  ml_mu <- predict(ml_tr, type = "location", newdata = testdata)
  ml_sigma <- predict(ml_tr, type = "scale", newdata = testdata)
  ml_nu <- predict(ml_bin, newdata = testdata, what = "mu", type = "response", data = learndata)
  ml_exp <- (ml_mu + ml_sigma * (dnorm(ml_mu/ml_sigma) / pnorm(ml_mu/ml_sigma))) * ml_nu
  
  
  # CPRS
  {
    crps_dt_h <- crps_dt_2p <- crps_df_h <- crps_df_2p <- 
      crps_g  <- crps_gb <- crps_ml <- numeric(length = NROW(testdata))
    
    int_upper <- ceiling(max(RainData$robs))  # alternative: quantile(RainData$robs, probs = 0.999)
    
    for(j in 1:(NROW(testdata))){
      crps_dt_h[j] <- integrate(function(z){(1 - dt_h_nu[j] + dt_h_nu[j] * ptnorm(z, mean = dt_h_mu[j], sd = dt_h_sigma[j], left = 0) - (testdata$robs[j] <= z))^2},
                                lower = 0, upper = int_upper)$value
      crps_dt_2p[j] <- integrate(function(z){(1 - dt_2p_nu[j] + dt_2p_nu[j] * ptnorm(z, mean = dt_2p_mu[j], sd = dt_2p_sigma[j], left = 0) - (testdata$robs[j] <= z))^2},
                                 lower = 0, upper = int_upper)$value
      crps_df_h[j] <- integrate(function(z){(1 - df_h_nu[j] + df_h_nu[j] * ptnorm(z, mean = df_h_mu[j], sd = df_h_sigma[j], left = 0) - (testdata$robs[j] <= z))^2},
                                lower = 0, upper = int_upper)$value
      crps_df_2p[j] <- integrate(function(z){(1 - df_2p_nu[j] + df_2p_nu[j] * ptnorm(z, mean = df_2p_mu[j], sd = df_2p_sigma[j], left = 0) - (testdata$robs[j] <= z))^2},
                                 lower = 0, upper = int_upper)$value
      crps_g[j] <- integrate(function(z){(1 - g_nu[j] + g_nu[j] * ptnorm(z, mean = g_mu[j], sd = g_sigma[j], left = 0) - (testdata$robs[j] <= z))^2},
                             lower = 0, upper = int_upper)$value
      crps_gb[j] <- integrate(function(z){(1 - gb_nu[j] + gb_nu[j] * ptnorm(z, mean = gb_mu[j], sd = gb_sigma[j], left = 0) - (testdata$robs[j] <= z))^2},
                              lower = 0, upper = int_upper)$value
      crps_ml[j] <- integrate(function(z){(1 - ml_nu[j] + ml_nu[j] * ptnorm(z, mean = ml_mu[j], sd = ml_sigma[j], left = 0) - (testdata$robs[j] <= z))^2},
                              lower = 0, upper = int_upper)$value
    }
    
    crps_dt_h <- mean(crps_dt_h, na.rm = TRUE)
    crps_dt_2p <- mean(crps_dt_2p, na.rm = TRUE)
    crps_df_h <- mean(crps_df_h, na.rm = TRUE)
    crps_df_2p <- mean(crps_df_2p, na.rm = TRUE)
    crps_g  <- mean(crps_g, na.rm = TRUE) 
    crps_gb <- mean(crps_gb, na.rm = TRUE)
    crps_ml <- mean(crps_ml, na.rm = TRUE)
  }
  
  
  ## error: expected value - observation
  expobs_dt_h <- abs(dt_h_exp - testdata[,"robs"])
  expobs_dt_2p <- abs(dt_2p_exp - testdata[,"robs"])
  expobs_df_h <- abs(df_h_exp - testdata[,"robs"])
  expobs_df_2p <- abs(df_2p_exp - testdata[,"robs"])
  expobs_g  <- abs(g_exp - testdata[,"robs"])
  expobs_gb <- abs(gb_exp - testdata[,"robs"])
  expobs_ml <- abs(ml_exp - testdata[,"robs"])
  
  # RMSE
  rmse_dt_h <- sqrt(mean(expobs_dt_h^2))
  rmse_dt_2p <- sqrt(mean(expobs_dt_2p^2))
  rmse_df_h <- sqrt(mean(expobs_df_h^2))
  rmse_df_2p <- sqrt(mean(expobs_df_2p^2))
  rmse_g <- sqrt(mean(expobs_g^2))
  rmse_gb <- sqrt(mean(expobs_gb^2))
  rmse_ml <- sqrt(mean(expobs_ml^2))
  
  # loglikelihood
  {
    dthll <- dt2pll <- dfhll <- df2pll <- gll <-  gbll <- mlll <- numeric(length = NROW(testdata))
    for(j in 1:(NROW(testdata))){
      
      eta_dt_h <- as.numeric(dist_list_hurdle_normal$linkfun(cbind(dt_h_mu, dt_h_sigma, dt_h_nu)[j,]))
      eta_dt_2p <- as.numeric(dist_list_hurdle_normal$linkfun(cbind(dt_2p_mu, dt_2p_sigma, dt_2p_nu)[j,]))
      eta_df_h <- as.numeric(dist_list_hurdle_normal$linkfun(cbind(df_h_mu, df_h_sigma, df_h_nu)[j,]))
      eta_df_2p <- as.numeric(dist_list_hurdle_normal$linkfun(cbind(df_2p_mu, df_2p_sigma, df_2p_nu)[j,]))
      eta_g  <- as.numeric(dist_list_hurdle_normal$linkfun(cbind(g_mu, g_sigma, g_nu)[j,])) 
      eta_gb <- as.numeric(dist_list_hurdle_normal$linkfun(cbind(gb_mu, gb_sigma, gb_nu)[j,]))  
      eta_ml <- as.numeric(dist_list_hurdle_normal$linkfun(cbind(ml_mu, ml_sigma, ml_nu)[j,]))
      
      dthll[j] <- dist_list_hurdle_normal$ddist(testdata[j,"robs"], eta = eta_dt_h, log=TRUE)
      dt2pll[j] <- dist_list_hurdle_normal$ddist(testdata[j,"robs"], eta = eta_dt_2p, log=TRUE)
      dfhll[j] <- dist_list_hurdle_normal$ddist(testdata[j,"robs"], eta = eta_df_h, log=TRUE)
      df2pll[j] <- dist_list_hurdle_normal$ddist(testdata[j,"robs"], eta = eta_df_2p, log=TRUE)
      gll[j]  <- dist_list_hurdle_normal$ddist(testdata[j,"robs"], eta = eta_g, log=TRUE)
      gbll[j] <- dist_list_hurdle_normal$ddist(testdata[j,"robs"], eta = eta_gb, log=TRUE)
      mlll[j] <- dist_list_hurdle_normal$ddist(testdata[j,"robs"], eta = eta_ml, log=TRUE)
      
    }
    
    dthll <- mean(dthll, na.rm = TRUE)
    dt2pll <- mean(dt2pll, na.rm = TRUE)
    dfhll <- mean(dfhll, na.rm = TRUE)
    df2pll <- mean(df2pll, na.rm = TRUE)
    gll  <- mean(gll, na.rm = TRUE) 
    gbll <- mean(gbll, na.rm = TRUE)
    mlll <- mean(mlll, na.rm = TRUE)
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
  
  ll <- c(dthll, dt2pll, dfhll, df2pll, gll, gbll, mlll) 
  rmse <- c(rmse_dt_h, rmse_dt_2p, rmse_df_h, rmse_df_2p, rmse_g, rmse_gb, rmse_ml)
  crps <- c(crps_dt_h, crps_dt_2p, crps_df_h, crps_df_2p, crps_g, crps_gb, crps_ml)
  
  # compare results
  results[[station]] <- rbind(ll, rmse, crps)
  colnames(results[[station]]) <- c("tree_h", "tree_2p", "forest_h", "forest_2p", "prespGAM", "boostGAM", "EMOS")
  #results
  
}

save(results, file = "~/svn/partykit/pkg/disttree/demo/Rain_hurdle.rda")



if(FALSE) {
  load("~/svn/partykit/pkg/disttree/demo/Rain_hurdle.rda")
  
  results$Axams["crps", c("forest_h", "forest_2p")]
  results$Lech["crps", c("forest_h", "forest_2p")]
  results$Walchsee["crps", c("forest_h", "forest_2p")]
  results$Oetz["crps", c("forest_h", "forest_2p")]
  results$Koessen["crps", c("forest_h", "forest_2p")]
  results$Vils["crps", c("forest_h", "forest_2p")]
  results$See["crps", c("forest_h", "forest_2p")]
  results$Jungholz["crps", c("forest_h", "forest_2p")]
  
}