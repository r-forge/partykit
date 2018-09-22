#######################################################
### Probabilistic Forecasting on Precipitation Data ###
#######################################################

## Replication material for supplement of: 
## Distributional Regression Forests for Probabilistic Precipitation Forecasting in Complex Terrain (2018)
## by Lisa Schlosser and Torsten Hothorn and Reto Stauffer and Achim Zeileis
## URL: http://arxiv.org/abs/1804.02921

## This demo includes the application of two different 
## versions of the forest hurdle model on selected stations
## (models learned on 24 years and evaluated on 4 years)
## Full replication of all other results can be obtained with
## demo("RainTyrol", package = "disttree")

## Computation time: approximately ... (on our machines, using .. kernel)


library("disttree")
library("RainTyrol")
library("crch")

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


#####
# further packages
library("scoringRules")



# select stations
stationlist <- c("Axams", "Lech", "Zuers", "See im Paznaun", "Jungholz", 
                 "Ladis-Neuegg", "Oetz", "Ochsengarten-Obergut",
                 "Ginzling", "Rotholz", "Walchsee", "Koessen", 
                 "Innervillgraten", "Matrei in Osttirol", 
                 "St.Johann im Walde")

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
    
    int_upper <- ceiling(max(RainData$robs))  # alternative: quantile(RainData$robs, probs = 0.999)
    
    for(j in 1:(NROW(testdata))){
      crps_dt_h <- try(integrate(function(z){(1 - dt_h_nu[j] + dt_h_nu[j] * ptnorm(z, mean = dt_h_mu[j], sd = dt_h_sigma[j], left = 0) - (testdata$robs[j] <= z))^2},
                                lower = 0, upper = int_upper)$value, silent = TRUE)
      if(inherits(crps_dt_h, "try-error")) {
        crps_dt_h[j] <- NA 
        crps_dt_h_error <- TRUE
      } else crps_dt_h[j] <- crps_dt_h
      crps_dt_2p[j] <- integrate(function(z){(1 - dt_2p_nu[j] + dt_2p_nu[j] * ptnorm(z, mean = dt_2p_mu[j], sd = dt_2p_sigma[j], left = 0) - (testdata$robs[j] <= z))^2},
                                 lower = 0, upper = int_upper)$value
      crps_df_h[j] <- integrate(function(z){(1 - df_h_nu[j] + df_h_nu[j] * ptnorm(z, mean = df_h_mu[j], sd = df_h_sigma[j], left = 0) - (testdata$robs[j] <= z))^2},
                                lower = 0, upper = int_upper)$value
      crps_df_2p[j] <- integrate(function(z){(1 - df_2p_nu[j] + df_2p_nu[j] * ptnorm(z, mean = df_2p_mu[j], sd = df_2p_sigma[j], left = 0) - (testdata$robs[j] <= z))^2},
                                 lower = 0, upper = int_upper)$value
    }
    
    crps_dt_h <- if(crps_dt_h_error) NA else mean(crps_dt_h, na.rm = TRUE)
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
  results[[station]] <- rbind(ll, rmse, crps)
  colnames(results[[station]]) <- c("tree_h", "tree_2p", "forest_h", "forest_2p")
  #results
  
}


save(results, file = "Rain_hurdleversions_forest.rda")

crps_forests <- matrix(ncol = 2, nrow = length(names(results)))
rownames(crps_forests) <- names(results)
colnames(crps_forests) <- c("forest_h", "forest_2p")
for(i in 1:length(names(results))){
  crps_forests[i,] <- results[[i]]["crps", c("forest_h", "forest_2p")]
}

save(crps_forests, file = "crps_forests.rda")

if(FALSE){
  
  # HCL palette
  pal <- hcl(c(10, 128, 260, 290, 50), 100, 50)
  
  load("crps_forests.rda")
  ## CRPS skill score boxplot with reference two-part hurdle model
  crps_ss_forest_2p <- data.frame(1 - crps_forests[,"forest_2p"] / crps_forests[,"forest_h"])
  names(crps_ss_forest_2p) <- "Three-parametric forest"
  boxplot(crps_ss_forest_2p, ylab = "CRPS skill score", names = c("Three-parametric forest"))
  abline(h = 0, col = pal[5], lwd = 2)
  axis(1, 0:2, c("","Three-parametric forest", ""), las=1)
  

  sum(crps_forests[,"forest_h"] > crps_forests[,"forest_2p"])
  summary(crps_forests[,"forest_h"] - crps_forests[,"forest_2p"])
  mean(crps_forests[,"forest_h"] - crps_forests[,"forest_2p"])
  colMeans(crps_forests)

}