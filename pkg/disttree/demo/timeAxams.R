#######################################################
### Probabilistic Forecasting on Precipitation Data ###
#######################################################

## Replication material for: 
## Distributional Regression Forests for Probabilistic Precipitation Forecasting in Complex Terrain (2018)
## by Lisa Schlosser and Torsten Hothorn and Reto Stauffer and Achim Zeileis
## URL: http://arxiv.org/abs/1804.02921

## Timing computation time

## This demo includes the application on one station (Axams) 
## (models learned on 24 years and evaluated on 4 years)
## Full replication of all other results can be obtained with
## demo("RainTyrol", package = "disttree")

## Computation time: approximately 18 minutes (on our machines, using 1 kernel)


library("disttree")

#####
# load observations and covariates 
data("RainAxams", package = "disttree")


##### 
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


fit_time <- matrix(ncol = 5, nrow = 6)

colnames(fit_time) <- c("user.self", "sys.self", "elapsed", "user.child", "sys.child")
rownames(fit_time) <- c("disttree", "distforest", "gamlss", "gamboostLSS", "gamboostLSS_cvr", "EMOS")

#####
# fitting the models

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


#####
# fit other heteroscedastic censored gaussian models

# fit prespecified GAM (covariates selected based on meteorological expert knowledge)
g_learndata <- learndata
g_learndata$robs <- Surv(g_learndata$robs, g_learndata$robs>0, type="left")
fit_time["gamlss",] <- system.time(g <- gamlss(formula = g.mu.formula, sigma.formula = g.sigma.formula, data = g_learndata, 
                                               family = cens("NO", type = "left"),
                                               control = gamlss.control(n.cyc = 100),
                                               i.control = glim.control(cyc = 100, bf.cyc = 100)))

# fit boosted GAM
fit_time["gamboostLSS",] <- system.time(gb <- gamboostLSS(formula = list(mu = gb.mu.formula, sigma = gb.sigma.formula), data = g_learndata, 
                                                          families = as.families(fname = cens("NO", type = "left")), method = "noncyclic",
                                                          control = boost_control(mstop = 1000L)))
# find optimal value for mstop (evalution is very time-consuming)
grid <- seq(50,1000, by = 25)
fit_time["gamboostLSS_cvr",] <- system.time(cvr <- cvrisk(gb, grid = grid))
mstop(gb) <- mstop(cvr)

# fit linear model with only total precipitation as covariate (Ensemble Model Output Statistics, EMOS)
fit_time["EMOS",] <- system.time(ml <- crch(formula = robs ~ tppow_mean | log(tppow_sprd + 0.001), 
                                            data = learndata, dist = "gaussian", left = 0, link.scale = "log"))



###################################################
# measure computation time for one observation
# in a loop

dt_time <- df_time <- g_time <- gb_time <- ml_time <- matrix(ncol = 5, nrow = NROW(testdata))

for(i in 1:NROW(testdata)){
  dt_time[i,] <- system.time(predict(dt, newdata = testdata[i,], type = "parameter"))
  df_time[i,] <- system.time(predict(df, newdata = testdata[i,], type = "parameter"))
  g_time[i,] <- system.time(predict(g, newdata = testdata[i,], what = "mu", type = "response", data = g_learndata)) +
    system.time(predict(g, newdata = testdata[i,], what = "sigma", type = "response", data = g_learndata))
  gb_time[i,] <- system.time(predict(gb, newdata = testdata[i,], parameter = c("mu","sigma"), type = "response"))
  ml_time[i,] <- system.time(predict(ml, type = "location", newdata = testdata[i,])) + 
    system.time(predict(ml, type = "scale", newdata = testdata[i,]))
}

colnames(dt_time) <- colnames(df_time) <- colnames(g_time) <- colnames(gb_time) <- 
  colnames(ml_time) <- names(system.time(predict(dt, newdata = testdata[i,], type = "parameter")))  

time_user.self <- colMeans(cbind(dt_time[,"user.self"], 
                                 df_time[,"user.self"], 
                                 g_time[,"user.self"], 
                                 gb_time[,"user.self"], 
                                 ml_time[,"user.self"]))

time_sys.self <- colMeans(cbind(dt_time[,"sys.self"], 
                                 df_time[,"sys.self"], 
                                 g_time[,"sys.self"], 
                                 gb_time[,"sys.self"], 
                                 ml_time[,"sys.self"]))

time_elapsed <- colMeans(cbind(dt_time[,"elapsed"], 
                                 df_time[,"elapsed"], 
                                 g_time[,"elapsed"], 
                                 gb_time[,"elapsed"], 
                                 ml_time[,"elapsed"]))

time_user.child <- colMeans(cbind(dt_time[,"user.child"], 
                                 df_time[,"user.child"], 
                                 g_time[,"user.child"], 
                                 gb_time[,"user.child"], 
                                 ml_time[,"user.child"]))

time_sys.child <- colMeans(cbind(dt_time[,"sys.child"], 
                                df_time[,"sys.child"], 
                                g_time[,"sys.child"], 
                                gb_time[,"sys.child"], 
                                ml_time[,"sys.child"]))

names(time_user.self) <- names(time_sys.self) <- names(time_elapsed) <- 
  names(time_user.child) <- names(time_sys.child) <- c("disttree", "distforest", "gamlss", "gamboostLSS", "EMOS")


fit_time
time_user.self

time <- list(fit_time = fit_time,
             time_user.self = time_user.self,
             time_sys.self = time_sys.self,
             time_elapsed = time_elapsed,
             time_user.child = time_user.child,
             time_sys.child = time_sys.child)

save(time, file = "~/svn/partykit/pkg/disttree/demo/timeAxams.rda")
