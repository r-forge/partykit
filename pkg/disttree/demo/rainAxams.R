#######################################################
### Probabilistic Forecasting on Precipitation Data ###
#######################################################

library("disttree")
set.seed(7)


#####
# load observations and covariates 
data("rainAxams")


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
learndata <- rainAxams[rainAxams$year < 109,]
testdata <- rainAxams[rainAxams$year %in% c(109, 110, 111, 112),]



############
# fitting the models

# fit distributional tree
dt <- disttree(dt.formula, 
               data = learndata, family = dist_list_cens_normal, 
               censtype = "left", censpoint = 0, type.tree = "ctree", 
               control = ctree_control(teststat = "quad", testtype = "Bonferroni", intersplit = TRUE,
                                       mincriterion = 0.95, minsplit = 50,
                                       minbucket = 20))


# visualization
plot(dt)

# estimated coefficients
coef(dt)


# fit distributional forest
df <- distforest(df.formula, 
                 data = learndata, family = dist_list_cens_normal, type.tree = "ctree", 
                 ntree = 100, censtype = "left", censpoint = 0, mtry = 27,
                 control = ctree_control(teststat = "quad", testtype = "Univ", intersplit = TRUE,
                                         mincriterion = 0, minsplit = 50,
                                         minbucket = 20))



############
# predict expected total precipitation for one day in July of the 4 test years and
# plot corresponding predicted density functions for 2009, 2010, 2011 and 2012

# predictions for one day (in each of the four years) 
# here July 24
pday <- 24  
# (19th of July 2011 is missing)
pdays <- if(pday<19) c(pday, pday + 31, pday + 62, pday + 92) else c(pday, pday + 31, pday + 61, pday + 92)
# get predicted parameters
pdf <- predict(df, newdata = testdata[pdays,], type = "parameter")
df_mu <- pdf$mu
df_sigma <- pdf$sigma
# caclulate predicted expected value
df_exp <- pnorm(df_mu/df_sigma) * (df_mu + df_sigma * (dnorm(df_mu/df_sigma) / pnorm(df_mu/df_sigma)))
# predicted expected values together with the corresponding observations
cbind(df_exp, testdata[pdays,"robs"])


# prepare data for needed for the plot
# distribution functions
x <- c(0.01, sort(runif(500,0.01,8)))
y1 <- crch::dcnorm(x, mean = df_mu[1], sd = df_sigma[1], left = 0)
y2 <- crch::dcnorm(x, mean = df_mu[2], sd = df_sigma[2], left = 0)
y3 <- crch::dcnorm(x, mean = df_mu[3], sd = df_sigma[3], left = 0)
y4 <- crch::dcnorm(x, mean = df_mu[4], sd = df_sigma[4], left = 0)
dayending <- if(pday > 3) "th" else switch(pday, "1" = {dayending <- "st"}, "2" = {dayending <- "nd"}, "3" = {dayending <- "rd"})

# point mass (slightly shifted)
pm1 <- c(0.05, crch::dcnorm(-1, mean = df_mu[1], sd = df_sigma[1], left = 0))
pm2 <- c(0.01, crch::dcnorm(-1, mean = df_mu[2], sd = df_sigma[2], left = 0))
pm3 <- c(-0.03, crch::dcnorm(-1, mean = df_mu[3], sd = df_sigma[3], left = 0))
pm4 <- c(-0.07, crch::dcnorm(-1, mean = df_mu[4], sd = df_sigma[4], left = 0))

# predictions
pred1 <- c(testdata[pdays,"robs"][1], crch::dcnorm(testdata[pdays,"robs"][1], mean = df_mu[1], sd = df_sigma[1], left = 0))
pred2 <- c(testdata[pdays,"robs"][2], crch::dcnorm(testdata[pdays,"robs"][2], mean = df_mu[2], sd = df_sigma[2], left = 0))
pred3 <- c(testdata[pdays,"robs"][3], crch::dcnorm(testdata[pdays,"robs"][3], mean = df_mu[3], sd = df_sigma[3], left = 0))
pred4 <- c(testdata[pdays,"robs"][4], crch::dcnorm(testdata[pdays,"robs"][4], mean = df_mu[4], sd = df_sigma[4], left = 0))

#legendheight
lh1 <- crch::dcnorm(0.01, mean = df_mu[1], sd = df_sigma[1], left = 0)
lh2 <- crch::dcnorm(0.01, mean = df_mu[2], sd = df_sigma[2], left = 0)
lh3 <- crch::dcnorm(0.01, mean = df_mu[3], sd = df_sigma[3], left = 0)
lh4 <- crch::dcnorm(0.01, mean = df_mu[4], sd = df_sigma[4], left = 0)


## PLOT
plot(x = x, y = y1, type = "l", col = "red", 
     main = paste0("July ", pday), ylab = "Density", 
     xlab = expression(Total~precipitation~"["~mm^(1/1.6)~"/"~"24h"~"]"),
     ylim = c(0,max(y1, y2, y3, y4, pm1, pm2, pm3, pm4) + 0.01),
     xlim = c(-1.5,8))
lines(x = x, y = y2, type = "l", col = "blue")
lines(x = x, y = y3, type = "l", col = "darkgreen")
lines(x = x, y = y4, type = "l", col = "purple")
legend('topright', c("predicted distribution", "point mass at censoring point", "observation"),
       bty = "n", col = "black", lty = c(1, NA, NA), pch = c(NA, 19, 4), cex = 0.8)

# plot point mass
lines(x = c(pm1[1], pm1[1]), y = c(pm1[2], 0), col = "red", type = "l", lwd = 1)
lines(x = c(pm2[1], pm2[1]), y = c(pm2[2], 0), col = "blue", type = "l", lwd = 1)
lines(x = c(pm3[1], pm3[1]), y = c(pm3[2], 0), col = "darkgreen", type = "l", lwd = 1)
lines(x = c(pm4[1], pm4[1]), y = c(pm4[2], 0), col = "purple", type = "l", lwd = 1)
points(x = pm1[1], y = pm1[2], col = "red", pch = 19)
points(x = pm2[1], y = pm2[2], col = "blue", pch = 19)
points(x = pm3[1], y = pm3[2], col = "darkgreen", pch = 19)
points(x = pm4[1], y = pm4[2], col = "purple", pch = 19)

# plot predictions
points(x = pred1[1], y = pred1[2], col = "red", pch = 4)
points(x = pred2[1], y = pred2[2], col = "blue", pch = 4)
points(x = pred3[1], y = pred3[2], col = "darkgreen", pch = 4)
points(x = pred4[1], y = pred4[2], col = "purple", pch = 4)
lines(x = c(pred1[1], pred1[1]), y = c(pred1[2], 0), col = "darkgrey", type = "l", lty = 2)
lines(x = c(pred2[1], pred2[1]), y = c(pred2[2], 0), col = "darkgrey", type = "l", lty = 2)
lines(x = c(pred3[1], pred3[1]), y = c(pred3[2], 0), col = "darkgrey", type = "l", lty = 2)
lines(x = c(pred4[1], pred4[1]), y = c(pred4[2], 0), col = "darkgrey", type = "l", lty = 2)

# add labels
text(x = -0.8, y = lh1, labels = "2009", col = "red", cex = 0.8)
text(x = -0.8, y = lh2, labels = "2010", col = "blue", cex = 0.8)
text(x = -0.8, y = lh3, labels = "2011", col = "darkgreen", cex = 0.8)
text(x = -0.8, y = lh4, labels = "2012", col = "purple", cex = 0.8)









############
# Comparison to other heteroscedastic censored gaussian models

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
                  control = boost_control(mstop = 500L))
# find optimal value for mstop (evalution is very time-consuming)
grid <- seq(50,1000, by = 25)
cvr <- cvrisk(gb, grid = grid)
mstop(gb) <- mstop(cvr)

# fit linear model with only total precipitation as covariate (Ensemble Model Output Statistics, EMOS)
ml <- crch(formula = robs ~ tppow_mean | log(tppow_sprd + 0.001), 
           data = learndata, dist = "gaussian", left = 0, link.scale = "log")



## get predicted parameter for testdata

# distributional tree
pdt <- predict(dt, newdata = testdata, type = "parameter")
dt_mu <- pdt$mu
dt_sigma <- pdt$sigma
dt_exp <- pnorm(dt_mu/dt_sigma) * (dt_mu + dt_sigma * (dnorm(dt_mu/dt_sigma) / pnorm(dt_mu/dt_sigma)))

# distributional forest
pdf <- predict(df, newdata = testdata, type = "parameter")
df_mu <- pdf$mu
df_sigma <- pdf$sigma
df_exp <- pnorm(df_mu/df_sigma) * (df_mu + df_sigma * (dnorm(df_mu/df_sigma) / pnorm(df_mu/df_sigma)))

# prespecified GAM
g_mu <- predict(g, newdata = testdata, what = "mu", type = "response", data = g_learndata)
g_sigma <- predict(g, newdata = testdata, what = "sigma", type = "response", data = g_learndata)
g_exp <- pnorm(g_mu/g_sigma) * (g_mu + g_sigma * (dnorm(g_mu/g_sigma) / pnorm(g_mu/g_sigma)))

# boosted GAM
pgb <- predict(gb, newdata = testdata, parameter = c("mu","sigma"), type = "response")
gb_mu <- pgb$mu
gb_sigma <- pgb$sigma
gb_exp <- pnorm(gb_mu/gb_sigma) * (gb_mu + gb_sigma * (dnorm(gb_mu/gb_sigma) / pnorm(gb_mu/gb_sigma)))

# EMOS
ml_mu <- predict(ml, type = "location", newdata = testdata)
ml_sigma <- predict(ml, type = "scale", newdata = testdata)
ml_exp <- pnorm(ml_mu/ml_sigma) * (ml_mu + ml_sigma * (dnorm(ml_mu/ml_sigma) / pnorm(ml_mu/ml_sigma)))


# CPRS
crps_dt <- mean(crps_cnorm(testdata$robs, location = dt_mu, scale = dt_sigma, lower = 0, upper = Inf), na.rm = TRUE)
crps_df <- mean(crps_cnorm(testdata$robs, location = df_mu, scale = df_sigma, lower = 0, upper = Inf), na.rm = TRUE)
crps_g  <- mean(crps_cnorm(testdata$robs, location = g_mu, scale = g_sigma, lower = 0, upper = Inf), na.rm = TRUE)
crps_gb <- mean(crps_cnorm(testdata$robs, location = gb_mu, scale = gb_sigma, lower = 0, upper = Inf), na.rm = TRUE)
crps_ml <- mean(crps_cnorm(testdata$robs, location = ml_mu, scale = ml_sigma, lower = 0, upper = Inf), na.rm = TRUE) 

# RMSE
rmse_dt <- sqrt(mean((dt_exp - testdata[,"robs"])^2, na.rm = TRUE))
rmse_df <- sqrt(mean((df_exp - testdata[,"robs"])^2, na.rm = TRUE))
rmse_g  <- sqrt(mean((g_exp - testdata[,"robs"])^2, na.rm = TRUE))
rmse_gb <- sqrt(mean((gb_exp - testdata[,"robs"])^2, na.rm = TRUE))
rmse_ml <- sqrt(mean((ml_exp - testdata[,"robs"])^2, na.rm = TRUE))

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

dtll <- mean(dtll, na.rm = TRUE)
dfll <- mean(dfll, na.rm = TRUE)
gll  <- mean(gll, na.rm = TRUE) 
gbll <- mean(gbll, na.rm = TRUE)
mlll <- mean(mlll, na.rm = TRUE)

ll <- c(dtll, dfll, gll, gbll, mlll) 
rmse <- c(rmse_dt, rmse_df, rmse_g, rmse_gb, rmse_ml)
crps <- c(crps_dt, crps_df, crps_g, crps_gb, crps_ml)

# compare results
results <- rbind(ll, rmse, crps)
colnames(results) <- c("tree", "forest", "prespGAM", "boostGAM", "EMOS")
results





## PIT histograms

# prepare data
{
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
  
}




##
# PIT histograms
set.seed(4)

par(mfrow = c(2,3))
library("countreg")

pithist(pit_dt, nsim = 1000, breaks = seq(0, 1, length.out = 9), main = "Tree", ylim = c(0,1.5))
pithist(pit_df, nsim = 1000, breaks = seq(0, 1, length.out = 9), main = "Forest", ylim = c(0,1.5))
pithist(pit_g,  nsim = 1000, breaks = seq(0, 1, length.out = 9), main = "prespGAM", ylim = c(0,1.5))
pithist(pit_gb, nsim = 1000, breaks = seq(0, 1, length.out = 9), main = "boostGAM", ylim = c(0,1.5))
pithist(pit_ml, nsim = 1000, breaks = seq(0, 1, length.out = 9), main = "EMOS", ylim = c(0,1.5))


