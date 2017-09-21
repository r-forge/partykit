library("disttree")
library("gamlss")
library("gamlss.dist")
library("gamlss.cens")
gen.cens(NO, type = "left")
library("gamboostLSS")
library("bamlss")
library("crch")
library("scoringRules")
library("parallel")



# dist_list_cens_normal
{
  
  dist_list_cens_normal <- list()
  
  parnames <- c("mu", "sigma")
  etanames <- c("mu", "log(sigma)")
  
  ddist <-  function(y, eta, log = TRUE, weights = NULL, sum = FALSE, left = 0, right = Inf) {     
    par <- c(eta[1], exp(eta[2]))
    val <- crch::dcnorm(x = y, mean = par[1], sd = par[2], left = left, right = right, log = log)
    if(sum) {
      if(is.null(weights)) weights <- if(is.matrix(y)) rep.int(1, dim(y)[1]) else rep.int(1, length(y))
      val <- sum(weights * val, na.rm = TRUE)
    }
    return(val)
  }
  
  
  sdist <- function(y, eta, weights = NULL, sum = FALSE, left = 0, right = Inf) {   
    par <- c(eta[1], exp(eta[2]))
    # y[y==0] <- 1e-323
    
    score_m <- crch:::scnorm(x = y, mean = par[1], sd = par[2], which = "mu", left = left, right = right)
    score_s <- crch:::scnorm(x = y, mean = par[1], sd = par[2], which = "sigma", left = left, right = right) * exp(eta[2]) # inner derivation exp(eta[2])
    score <- cbind(score_m, score_s)
    score <- as.matrix(score)
    colnames(score) <- etanames
    if(sum) {
      if(is.null(weights)) weights <- rep.int(1, length(y)[1])
      # if score == Inf replace score with 1.7e308 because Inf*0 would lead to NaN (0 in weights)
      score[score==Inf] = 1.7e308
      score <- colSums(weights * score, na.rm = TRUE)
      #if(any(is.nan(score))) print(c(eta, "y", y))
    }
    return(score)
  }
  
  
  hdist <- function(y, eta, weights = NULL, left = 0, right = Inf) {    
    ny <- length(y)
    if(is.null(weights)) weights <- rep.int(1, ny)
    
    par <- c(eta[1], exp(eta[2]))                           
    # y[y==0] <- 1e-323
    
    d2mu <- crch:::hcnorm(x = y, mean = par[1], sd = par[2], which = "mu", left = left, right = right)
    d2sigma <- crch:::hcnorm(x = y, mean = par[1], sd = par[2], which = "sigma", left = left, right = right)
    dmudsigma <- crch:::hcnorm(x = y, mean = par[1], sd = par[2], which = "mu.sigma", left = left, right = right) # FIX: order?
    dsigmadmu <- crch:::hcnorm(x = y, mean = par[1], sd = par[2], which = "sigma.mu", left = left, right = right) # FIX: order?
    dsigma <- crch:::scnorm(x = y, mean = par[1], sd = par[2], which = "sigma", left = left, right = right)
    
    d2ld.etamu2 <- sum(weights * d2mu, na.rm = TRUE)
    d2ld.etamu.d.etasigma <- sum(weights * dmudsigma * par[2], na.rm = TRUE)
    d2ld.etasigma.d.etamu <- sum(weights * dsigmadmu * par[2], na.rm = TRUE)
    d2ld.etasigma2 <- sum(weights * (d2sigma * exp(2*eta[2]) + dsigma * par[2]), na.rm = TRUE)         
    
    hess <- matrix(c(d2ld.etamu2, d2ld.etamu.d.etasigma, d2ld.etasigma.d.etamu, d2ld.etasigma2), nrow = 2)
    colnames(hess) <- rownames(hess) <-  etanames
    
    return(hess)
  }
  
  
  ## additional functions pdist, qdist, rdist
  pdist <- function(q, eta, lower.tail = TRUE, log.p = FALSE) crch:::pcnorm(q, mean = eta[1], sd = eta[2], 
                                                                            lower.tail = lower.tail, log.p = log.p, 
                                                                            left = left, right = right)
  qdist <- function(p, eta, lower.tail = TRUE, log.p = FALSE) crch:::qcnorm(p, mean = eta[1], sd = eta[2], 
                                                                            lower.tail = lower.tail, log.p = log.p, 
                                                                            left = left, right = right)
  rdist <- function(n, eta) crch:::rcnorm(n, mean = eta[1], sd = eta[2], left = left, right = right)
  
  
  link <- c("identity", "log")
  
  linkfun <- function(par) {
    eta <- c(par[1], log(par[2]))
    names(eta) <- etanames
    return(eta)
  }
  
  
  linkinv <- function(eta) {
    par <- c(eta[1], exp(eta[2]))
    names(par) <- parnames
    return(par)
  }
  
  
  linkinvdr <- function(eta) {
    dpardeta <- c(1, exp(eta[2]))
    names(dpardeta) <- parnames
    return(dpardeta)
  }
  
  
  startfun <- function(y, weights = NULL){
    yc <- pmax(0,y)  # optional ?
    if(is.null(weights)) {
      mu <- mean(yc)
      sigma <- sqrt(1/length(yc) * sum((yc - mu)^2))
    } else {
      mu <- weighted.mean(yc, weights)
      sigma <- sqrt(1/sum(weights) * sum(weights * (yc - mu)^2))
    }
    starteta <- c(mu, log(sigma))
    names(starteta) <- etanames
    return(starteta)
  }
  
  mle <- FALSE
  
  dist_list_cens_normal <- list(family.name = "censored Normal Distribution",
                                ddist = ddist, 
                                sdist = sdist, 
                                hdist = hdist,
                                pdist = pdist,
                                qdist = qdist,
                                rdist = rdist,
                                link = link, 
                                linkfun = linkfun, 
                                linkinv = linkinv, 
                                linkinvdr = linkinvdr,
                                startfun = startfun,
                                mle = mle,
                                gamlssobj = FALSE,
                                censored = TRUE
  )
}



# load station list
#load("~/svn/partykit/pkg/disttree/inst/draft/rainData/ehyd.statlist.rda")
#load("~/disttree/inst/draft/rainData/ehyd.statlist.rda")
#cbind(ehyd.statlist$station, ehyd.statlist$stationname, ehyd.statlist$height)



# load observations
#load("~/svn/partykit/pkg/disttree/inst/draft/rainData/rain.rda")
load("~/disttree/inst/draft/rainData/rain.rda")
#head(data)
Sys.setenv("TZ" = "UTC")
rain  <- data
#dim(rain)
rain$robs <- (rain$obs)^(1/1.6)
rain$day <- as.POSIXlt(rain$date)$yday
rain$year <- as.POSIXlt(rain$date)$year
rain$hour <- as.POSIXlt(rain$date)$hour
#head(rain)
## save complete data
rain_all <- rain

## choose one station: 
# 9101121 Lech
# 9103135 Fieberbrunn
# 9102954 Nassereith
# 9102756 Steinberg am Rofan (EHYD)  
# 9101113 Zuers (EHYD)

# 9102913 Walchsee (EHYD)
# 9102905 Koessen (EHYD)
# 9103309 Axams (EHYD)
# 9103176 Niederndorferberg (EHYD)
# 9119479 Ammerwald (EHYD)



# load predictions
#load("~/svn/partykit/pkg/disttree/inst/draft/rainData/prepared/GEFSV2_prepared_9103309.rda")
load("~/disttree/inst/draft/rainData/prepared/GEFSV2_prepared_9103309.rda")
prediction <- prepared
Sys.setenv("TZ" = "UTC")
prediction$day <-as.POSIXlt(prediction$init)$yday
prediction$year <- as.POSIXlt(prediction$init)$year
#head(prediction)

rain <- rain[rain$station == "Axams (EHYD)",]
#head(rain)

## choose one month:
# July
rain <- rain[(((rain$year %% 4 == 0) & rain$day>=183) | (!(rain$year %% 4 == 0) & rain$day>=182)),]
rain <- rain[(((rain$year %% 4 == 0) & rain$day<=213) | (!(rain$year %% 4 == 0) & rain$day<=212)),]

# predictions start in 1985
rain <- rain[rain$year>=85,]

# observations end in 2012
prediction <- prediction[prediction$year<=112,]

#head(rain)
#dim(rain)





# combine data frames

# predictions with init 07-01 are for day 07-02 
# observations are the sum of precipitation of the last 24h
# match prediction with init 07-01 with observation on 07-02 
# (pred starts at 00UTC and predicts from 06UTC until 30UTC = 06UTC of next day, 
# observations are meassured from 06UTC of previous day to 06UTC of current day)

#dim(rain)
#dim(prediction)
#head(rain)
#head(prediction[,c(1:10)])

raindata <- cbind(rain$date, rain$obs, rain$robs, prediction)
#head(raindata[,c(1:10)])
colnames(raindata)[c(1:3)] <- c("date", "obs", "robs")








########################################
# select variables with sufficient data

#dim(raindata)
#colnames(raindata)


if(FALSE)
{
  #total precipitation mean
  which(is.na(raindata$tp_mean))
  which(is.na(raindata$tp_mean0612))
  which(is.na(raindata$tp_mean1218))
  which(is.na(raindata$tp_mean1824))
  which(is.na(raindata$tp_mean2430))
  
  # total precipitation spread
  which(is.na(raindata$tp_sprd))
  which(is.na(raindata$tp_sprd0612))
  which(is.na(raindata$tp_sprd1218))
  which(is.na(raindata$tp_sprd1824))
  which(is.na(raindata$tp_sprd2430))
  
  #total precipitation min and max
  which(is.na(raindata$tp_min))
  which(is.na(raindata$tp_max))
  
  #total precipitation frac
  which(is.na(raindata$tp_frac))
  which(is.na(raindata$tp_frac0612))
  which(is.na(raindata$tp_frac1218))
  which(is.na(raindata$tp_frac1824))
  which(is.na(raindata$tp_frac2430))
  
  #total cape mean
  which(is.na(raindata$cape_mean))
  which(is.na(raindata$cape_mean0612))
  which(is.na(raindata$cape_mean1218))
  which(is.na(raindata$cape_mean1224))
  which(is.na(raindata$cape_mean1230))
  
  # total cape spread
  which(is.na(raindata$cape_sprd))
  which(is.na(raindata$cape_sprd0612))
  which(is.na(raindata$cape_sprd1218))
  which(is.na(raindata$cape_sprd1224))
  which(is.na(raindata$cape_sprd1230))
  
  #total cape min and max
  which(is.na(raindata$cape_min))
  which(is.na(raindata$cape_max))
  
  #total cape frac
  which(is.na(raindata$cape_frac))
  which(is.na(raindata$cape_frac0612))
  which(is.na(raindata$cape_frac1218))
  which(is.na(raindata$cape_frac1824))
  which(is.na(raindata$cape_frac2430))
  
  # sunshine
  which(is.na(raindata$dswrf_mean_mean))
  which(is.na(raindata$dswrf_mean_min))
  which(is.na(raindata$dswrf_mean_max))
  which(is.na(raindata$dswrf_sprd_mean))
  which(is.na(raindata$dswrf_sprd_min))
  which(is.na(raindata$dswrf_sprd_max))
  
  #mean sealevel pressure
  which(is.na(raindata$msl_mean_mean))
  which(is.na(raindata$msl_mean_min))
  which(is.na(raindata$msl_mean_max))
  which(is.na(raindata$msl_sprd_mean))
  which(is.na(raindata$msl_sprd_min))
  which(is.na(raindata$msl_sprd_max))
  
  # preciptal water
  which(is.na(raindata$pwat_mean_mean))
  which(is.na(raindata$pwat_mean_min))
  which(is.na(raindata$pwat_mean_max))
  which(is.na(raindata$pwat_sprd_mean))
  which(is.na(raindata$pwat_sprd_min))
  which(is.na(raindata$pwat_sprd_max))
  
  # total column-integrated condensate
  which(is.na(raindata$tcolc_mean_mean))
  which(is.na(raindata$tcolc_mean_min))
  which(is.na(raindata$tcolc_mean_max))
  which(is.na(raindata$tcolc_sprd_mean))
  which(is.na(raindata$tcolc_sprd_min))
  which(is.na(raindata$tcolc_sprd_max))
  
  # 2m maximum temperature
  which(is.na(raindata$tmax_mean_mean))
  which(is.na(raindata$tmax_mean_min))
  which(is.na(raindata$tmax_mean_max))
  which(is.na(raindata$tmax_sprd_mean))
  which(is.na(raindata$tmax_sprd_min))
  which(is.na(raindata$tmax_sprd_max))

  # temperature on 850 hPA
  which(is.na(raindata$t850_mean_mean))
  which(is.na(raindata$t850_mean_min))
  which(is.na(raindata$t850_mean_max))
  which(is.na(raindata$t850_sprd_mean))
  which(is.na(raindata$t850_sprd_min))
  which(is.na(raindata$t850_sprd_max))
  
  # temperature on 700 hPA
  which(is.na(raindata$t700_mean_mean))
  which(is.na(raindata$t700_mean_min))
  which(is.na(raindata$t700_mean_max))
  which(is.na(raindata$t700_sprd_mean))
  which(is.na(raindata$t700_sprd_min))
  which(is.na(raindata$t700_sprd_max))
  
  # temperature on 500 hPA
  which(is.na(raindata$t500_mean_mean))
  which(is.na(raindata$t500_mean_min))
  which(is.na(raindata$t500_mean_max))
  which(is.na(raindata$t500_sprd_mean))
  which(is.na(raindata$t500_sprd_min))
  which(is.na(raindata$t500_sprd_max))
  
  # temperature differences
  which(is.na(raindata$tdiff700850_mean))
  which(is.na(raindata$tdiff700850_min))
  which(is.na(raindata$tdiff700850_max))
  which(is.na(raindata$tdiff500850_mean))
  which(is.na(raindata$tdiff500850_min))
  which(is.na(raindata$tdiff500850_max))
  which(is.na(raindata$tdiff500700_mean))
  which(is.na(raindata$tdiff500700_min))
  which(is.na(raindata$tdiff500700_max))
}


######################################################
#remove row 825
raindata <- raindata[-825,]
######################################################
#remove rows 223-310
raindata <- raindata[-c(223:310),]
######################################################
#remove rows 111 176 178 220 221 241 320 364 367
raindata <- raindata[-c(111, 176, 178, 220, 221, 241, 320, 364, 367),]
######################################################
#remove rows 5 215 217 301 314
raindata <- raindata[-c(5, 215, 217, 301, 314),]
######################################################
#remove rows 251 264 
raindata <- raindata[-c(251, 264),]


######################################################
#only keep variables with sufficient values
raindata <- raindata[, c("robs",
                         "tppow_mean", "tppow_sprd", "tppow_min", "tppow_max", "tp_frac", 
                         "tppow_mean0612", "tppow_mean1218", "tppow_mean1824", "tppow_mean2430", 
                         "tppow_sprd0612", "tppow_sprd1218", "tppow_sprd1824", "tppow_sprd2430",
                         "capepow_mean", "capepow_sprd", "capepow_min", "capepow_max", "cape_frac",
                         "capepow_mean0612", "capepow_mean1218", "capepow_mean1224", "capepow_mean1230",
                         "capepow_sprd0612", "capepow_sprd1218", "capepow_sprd1224", "capepow_sprd1230",
                         "dswrf_mean_mean", "dswrf_mean_min", "dswrf_mean_max",
                         "dswrf_sprd_mean", "dswrf_sprd_min", "dswrf_sprd_max",
                         "msl_mean_mean", "msl_mean_min", "msl_mean_max",
                         "msl_sprd_mean", "msl_sprd_min", "msl_sprd_max",
                         "pwat_mean_mean", "pwat_mean_min", "pwat_mean_max",
                         "pwat_sprd_mean", "pwat_sprd_min", "pwat_sprd_max",
                         "tmax_mean_mean", "tmax_mean_min", "tmax_mean_max",
                         "tmax_sprd_mean", "tmax_sprd_min", "tmax_sprd_max",
                         "tcolc_mean_mean", "tcolc_mean_min", "tcolc_mean_max",
                         "tcolc_sprd_mean", "tcolc_sprd_min", "tcolc_sprd_max",
                         "t500_mean_mean", "t500_mean_min", "t500_mean_max",
                         "t700_mean_mean", "t700_mean_min", "t700_mean_max",
                         "t850_mean_mean", "t850_mean_min", "t850_mean_max",
                         "t500_sprd_mean", "t500_sprd_min", "t500_sprd_max",
                         "t700_sprd_mean", "t700_sprd_min", "t700_sprd_max",
                         "t850_sprd_mean", "t850_sprd_min", "t850_sprd_max",
                         "tdiff500850_mean", "tdiff500850_min", "tdiff500850_max",
                         "tdiff700850_mean", "tdiff700850_min", "tdiff700850_max",
                         "tdiff500700_mean", "tdiff500700_min", "tdiff500700_max")]

######################################################
#new vairable msl_diff
raindata$msl_diff <- raindata$msl_mean_max - raindata$msl_mean_min


#dim(raindata)
#any(is.na(raindata))

# tree and forest formula
dt.formula <- df.formula <- robs ~ tppow_mean + tppow_sprd + tppow_min + tppow_max + #tp_frac + 
  tppow_mean0612 + tppow_mean1218 + tppow_mean1824 + tppow_mean2430 + 
  tppow_sprd0612 + tppow_sprd1218 + tppow_sprd1824 + tppow_sprd2430 + 
  capepow_mean + capepow_sprd + capepow_min + capepow_max + #cape_frac +
  capepow_mean0612 + capepow_mean1218 + capepow_mean1224 + capepow_mean1230 +
  capepow_sprd0612 + capepow_sprd1218 + capepow_sprd1224 + capepow_sprd1230 +
  dswrf_mean_mean + dswrf_mean_max + #dswrf_mean_min + 
  dswrf_sprd_mean + dswrf_sprd_max + #dswrf_sprd_min +
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
  tdiff500700_mean + tdiff500700_min + tdiff500700_max

# gamlss formula
g.mu.formula <- robs ~ pb(tppow_mean) + 
  pb(tppow_mean1218 * capepow_mean1218) + 
  pb(tppow_max) + 
  pb(dswrf_mean_mean) +
  pb(tcolc_mean_mean) + 
  pb(msl_diff) + 
  pb(pwat_mean_mean) + 
  pb(tdiff500850_mean) #+  pb(tdiff700850_mean) 

g.sigma.formula <- ~ pb(tppow_sprd) + 
  pb(tppow_sprd1218 * capepow_mean1218) + 
  pb(dswrf_sprd_mean) +
  pb(tcolc_sprd_mean) + 
  pb(tdiff500850_mean) #+  pb(tdiff700850_mean) 

  

# gamboostLSS formula
gb.mu.formula <- robs ~ bbs(tppow_mean) + bbs(tppow_sprd) + bbs(tppow_min) + bbs(tppow_max) + #bbs(tp_frac) +  
  bbs(tppow_mean0612) + bbs(tppow_mean1218) + bbs(tppow_mean1824) + bbs(tppow_mean2430) + 
  bbs(tppow_sprd0612) + bbs(tppow_sprd1218) + bbs(tppow_sprd1824) + bbs(tppow_sprd2430) +
  bbs(capepow_mean) + bbs(capepow_sprd) + bbs(capepow_min) + bbs(capepow_max) + #bbs(cape_frac) +
  bbs(capepow_mean0612) + bbs(capepow_mean1218) + bbs(capepow_mean1224) + bbs(capepow_mean1230) +
  bbs(capepow_sprd0612) + bbs(capepow_sprd1218) + bbs(capepow_sprd1224) + bbs(capepow_sprd1230) +
  bbs(dswrf_mean_mean) + bbs(dswrf_mean_max) +  #bbs(dswrf_mean_min) +
  bbs(dswrf_sprd_mean) + bbs(dswrf_sprd_max) + #bbs(dswrf_sprd_min) +
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
  bbs(tdiff500700_mean) + bbs(tdiff500700_min) + bbs(tdiff500700_max)

gb.sigma.formula <- robs ~ bbs(tppow_mean) + bbs(tppow_sprd) + bbs(tppow_min) + bbs(tppow_max) + #bbs(tp_frac) +  
  bbs(tppow_mean0612) + bbs(tppow_mean1218) + bbs(tppow_mean1824) + bbs(tppow_mean2430) + 
  bbs(tppow_sprd0612) + bbs(tppow_sprd1218) + bbs(tppow_sprd1824) + bbs(tppow_sprd2430) +
  bbs(capepow_mean) + bbs(capepow_sprd) + bbs(capepow_min) + bbs(capepow_max) + #bbs(cape_frac) +
  bbs(capepow_mean0612) + bbs(capepow_mean1218) + bbs(capepow_mean1224) + bbs(capepow_mean1230) +
  bbs(capepow_sprd0612) + bbs(capepow_sprd1218) + bbs(capepow_sprd1224) + bbs(capepow_sprd1230) +
  bbs(dswrf_mean_mean) + bbs(dswrf_mean_max) +  #bbs(dswrf_mean_min) +
  bbs(dswrf_sprd_mean) + bbs(dswrf_sprd_max) + #bbs(dswrf_sprd_min) +
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
  bbs(tdiff500700_mean) + bbs(tdiff500700_min) + bbs(tdiff500700_max)


# bamlss formula
b.mu.formula <- robs ~ s(tppow_mean) + s(tppow_sprd) + s(tppow_min) + s(tppow_max) + #s(tp_frac) +  
  s(tppow_mean0612) + s(tppow_mean1218) + s(tppow_mean1824) + s(tppow_mean2430) + 
  s(tppow_sprd0612) + s(tppow_sprd1218) + s(tppow_sprd1824) + s(tppow_sprd2430) +
  s(capepow_mean) + s(capepow_sprd) + s(capepow_min) + s(capepow_max) + #s(cape_frac) +
  s(capepow_mean0612) + s(capepow_mean1218) + s(capepow_mean1224) + s(capepow_mean1230) +
  s(capepow_sprd0612) + s(capepow_sprd1218) + s(capepow_sprd1224) + s(capepow_sprd1230) +
  s(dswrf_mean_mean) + s(dswrf_mean_max) +  #s(dswrf_mean_min) +
  s(dswrf_sprd_mean) + s(dswrf_sprd_max) + #s(dswrf_sprd_min) +
  s(msl_mean_mean) + s(msl_mean_min) + s(msl_mean_max) + 
  s(msl_sprd_mean) + s(msl_sprd_min) + s(msl_sprd_max) +
  s(pwat_mean_mean) + s(pwat_mean_min) + s(pwat_mean_max) + 
  s(pwat_sprd_mean) + s(pwat_sprd_min) + s(pwat_sprd_max) +
  s(tmax_mean_mean) + s(tmax_mean_min) + s(tmax_mean_max) +
  s(tmax_sprd_mean) + s(tmax_sprd_min) + s(tmax_sprd_max) +
  s(tcolc_mean_mean) + s(tcolc_mean_min) + s(tcolc_mean_max) +
  s(tcolc_sprd_mean) + s(tcolc_sprd_min) + s(tcolc_sprd_max) +
  s(t500_mean_mean) + s(t500_mean_min) + s(t500_mean_max) +
  s(t700_mean_mean) + s(t700_mean_min) + s(t700_mean_max) +
  s(t850_mean_mean) + s(t850_mean_min) + s(t850_mean_max) +
  s(t500_sprd_mean) + s(t500_sprd_min) + s(t500_sprd_max) +
  s(t700_sprd_mean) + s(t700_sprd_min) + s(t700_sprd_max) +
  s(t850_sprd_mean) + s(t850_sprd_min) + s(t850_sprd_max) +
  s(tdiff500850_mean) + s(tdiff500850_min) + s(tdiff500850_max) +
  s(tdiff700850_mean) + s(tdiff700850_min) + s(tdiff700850_max) +
  s(tdiff500700_mean) + s(tdiff500700_min) + s(tdiff500700_max)

b.sigma.formula <- ~ s(tppow_mean) + s(tppow_sprd) + s(tppow_min) + s(tppow_max) + #s(tp_frac) +  
  s(tppow_mean0612) + s(tppow_mean1218) + s(tppow_mean1824) + s(tppow_mean2430) + 
  s(tppow_sprd0612) + s(tppow_sprd1218) + s(tppow_sprd1824) + s(tppow_sprd2430) +
  s(capepow_mean) + s(capepow_sprd) + s(capepow_min) + s(capepow_max) + #s(cape_frac) +
  s(capepow_mean0612) + s(capepow_mean1218) + s(capepow_mean1224) + s(capepow_mean1230) +
  s(capepow_sprd0612) + s(capepow_sprd1218) + s(capepow_sprd1224) + s(capepow_sprd1230) +
  s(dswrf_mean_mean) + s(dswrf_mean_max) +  #s(dswrf_mean_min) +
  s(dswrf_sprd_mean) + s(dswrf_sprd_max) + #s(dswrf_sprd_min) +
  s(msl_mean_mean) + s(msl_mean_min) + s(msl_mean_max) + 
  s(msl_sprd_mean) + s(msl_sprd_min) + s(msl_sprd_max) +
  s(pwat_mean_mean) + s(pwat_mean_min) + s(pwat_mean_max) + 
  s(pwat_sprd_mean) + s(pwat_sprd_min) + s(pwat_sprd_max) +
  s(tmax_mean_mean) + s(tmax_mean_min) + s(tmax_mean_max) +
  s(tmax_sprd_mean) + s(tmax_sprd_min) + s(tmax_sprd_max) +
  s(tcolc_mean_mean) + s(tcolc_mean_min) + s(tcolc_mean_max) +
  s(tcolc_sprd_mean) + s(tcolc_sprd_min) + s(tcolc_sprd_max) +
  s(t500_mean_mean) + s(t500_mean_min) + s(t500_mean_max) +
  s(t700_mean_mean) + s(t700_mean_min) + s(t700_mean_max) +
  s(t850_mean_mean) + s(t850_mean_min) + s(t850_mean_max) +
  s(t500_sprd_mean) + s(t500_sprd_min) + s(t500_sprd_max) +
  s(t700_sprd_mean) + s(t700_sprd_min) + s(t700_sprd_max) +
  s(t850_sprd_mean) + s(t850_sprd_min) + s(t850_sprd_max) +
  s(tdiff500850_mean) + s(tdiff500850_min) + s(tdiff500850_max) +
  s(tdiff700850_mean) + s(tdiff700850_min) + s(tdiff700850_max) +
  s(tdiff500700_mean) + s(tdiff500700_min) + s(tdiff500700_max)




# evaluate models on the whole data set
if(FALSE){
  
  dt <- disttree(dt.formula, 
                 data = raindata, family = dist_list_cens_normal, 
                 censtype = "left", censpoint = 0, type.tree = "ctree")
  ## FIX ME: error if type.tree = "mob"
  
  df <- distforest(df.formula, 
                   data = raindata, family = dist_list_cens_normal, ntree = 200, #mtry = 40,
                   censtype = "left", censpoint = 0, type.tree = "ctree")
  
  g_raindata <- raindata
  g_raindata$robs <- Surv(g_raindata$robs, g_raindata$robs>0, type="left")
  g <- gamlss(formula = g.mu.formula, sigma.formula = g.sigma.formula, data = g_raindata, family = cens("NO", type = "left"))
  
  gb <- gamboostLSS(formula = list(mu = gb.mu.formula, sigma = gb.sigma.formula), data = g_raindata, 
                    families = as.families(fname = cens("NO", type = "left")), method = "noncyclic",
                    control = boost_control(mstop = 400L))
  
  b <- bamlss(list(b.mu.formula, b.sigma.formula), family = "cnorm", 
              data = raindata, sampler = FALSE, optimizer = boost, 
              stop.criterion = "BIC", plot = FALSE) #, nu = 0.1)
  
  mi <- crch(formula = robs ~ tppow_mean | tppow_sprd, 
             data = raindata, dist = "gaussian", left = 0, link.scale = "identity")
  
  ml <- crch(formula = robs ~ tppow_mean | log(tppow_sprd + 0.001), 
             data = raindata, dist = "gaussian", left = 0, link.scale = "log")
  
  mq <- crch(formula = robs ~ tppow_mean | I(tppow_sprd^2), 
             data = raindata, dist = "gaussian", left = 0, link.scale = "quadratic")
  
  
  
  ## get predicted parameter
  
  # disttree
  if(is.vector(coef(dt))){
    pdt <- t(as.data.frame(coef(dt)))[as.vector(predict(dt, type = "node")),]
    rownames(pdt) <- paste(c(1 : NROW(pdt)))
  } else {
    pdt <- coef(dt)[paste(as.vector(predict(dt, type = "node"))),]
  }
  dt_mu <- as.vector(pdt[,"mu"]) 
  dt_sigma <- as.vector(pdt[,"sigma"]) 
  
  
  #distforest
  pdf <- predict(df, type = "parameter")
  df_mu <- pdf$mu
  df_sigma <- pdf$sigma
  
  #gamlss
  g_mu <- predict(g, what = "mu", type = "response", data = g_raindata)
  g_sigma <- predict(g, what = "sigma", type = "response", data = g_raindata)
  
  #bamlss
  pb <- predict(b, type = "parameter")
  b_mu <- pb$mu
  b_sigma <- pb$sigma
  
  #gamboostLSS
  pgb <- predict(gb, parameter = list("mu","sigma"), type = "response")
  gb_mu <- pgb[[1]]
  gb_sigma <- pgb[[2]]
  
  #EMOS
  mi_mu <- predict(mi, type = "location")     # returns parameter on response scale
  mi_sigma <- predict(mi, type = "scale")
  
  ml_mu <- predict(ml, type = "location")     # returns parameter on response scale
  ml_sigma <- predict(ml, type = "scale")
  
  mq_mu <- predict(mq, type = "location")     # returns parameter on response scale
  mq_sigma <- predict(mq, type = "scale")
  
  
  
  ### CRPS
  dt_crps <- sum(crps_cnorm(raindata$robs, location = dt_mu, scale = dt_sigma, lower = 0, upper = Inf))
  df_crps <- sum(crps_cnorm(raindata$robs, location = df_mu, scale = df_sigma, lower = 0, upper = Inf))
  g_crps <- sum(crps_cnorm(raindata$robs, location = g_mu, scale = g_sigma, lower = 0, upper = Inf))
  b_crps <- sum(crps_cnorm(raindata$robs, location = b_mu, scale = b_sigma, lower = 0, upper = Inf))
  gb_crps <- sum(crps_cnorm(raindata$robs, location = gb_mu, scale = gb_sigma, lower = 0, upper = Inf))
  mi_crps <- sum(crps_cnorm(raindata$robs, location = mi_mu, scale = mi_sigma, lower = 0, upper = Inf))
  ml_crps <- sum(crps_cnorm(raindata$robs, location = ml_mu, scale = ml_sigma, lower = 0, upper = Inf))
  mq_crps <- sum(crps_cnorm(raindata$robs, location = mq_mu, scale = mq_sigma, lower = 0, upper = Inf))
  
  #c(dt_crps, df_crps, g_crps, b_crps, gb_crps, mi_crps, ml_crps, mq_crps)
  c(dt_crps, df_crps, g_crps, gb_crps, mi_crps, ml_crps, mq_crps)
  
  
  # pit histogram for one station
  if(FALSE){
    # disttree
    #hist(pnorm(raindata[,"robs"], dt_mu, dt_sigma))
    pit_dt <- pnorm(raindata[,"robs"], dt_mu, dt_sigma)
    pit_dt[which(raindata[,"robs"]==0)] <- pit_dt[which(raindata[,"robs"]==0)]*runif(length(pit_dt[which(raindata[,"robs"]==0)]),0,1)
    hist(pit_dt)
    
    # distforest
    #hist(pnorm(testdata[,"robs"], df_mu, df_sigma))
    pit_df <- pnorm(raindata[,"robs"], df_mu, df_sigma)
    pit_df[which(raindata[,"robs"]==0)] <- pit_df[which(raindata[,"robs"]==0)]*runif(length(pit_df[which(raindata[,"robs"]==0)]),0,1)
    hist(pit_df)
    
    # gamlss
    #hist(pnorm(raindata[,"robs"], g_mu, g_sigma))
    pit_g <- pnorm(raindata[,"robs"], g_mu, g_sigma)
    pit_g[which(raindata[,"robs"]==0)] <- pit_g[which(raindata[,"robs"]==0)]*runif(length(pit_g[which(raindata[,"robs"]==0)]),0,1)
    hist(pit_g)
    
    # bamlss
    #hist(pnorm(raindata[,"robs"], b_mu, b_sigma))
    pit_b <- pnorm(raindata[,"robs"], b_mu, b_sigma)
    pit_b[which(raindata[,"robs"]==0)] <- pit_b[which(raindata[,"robs"]==0)]*runif(length(pit_b[which(raindata[,"robs"]==0)]),0,1)
    hist(pit_b)
    
    # gamboostLSS
    #hist(pnorm(raindata[,"robs"], gb_mu, gb_sigma))
    pit_gb <- pnorm(raindata[,"robs"], gb_mu, gb_sigma)
    pit_gb[which(raindata[,"robs"]==0)] <- pit_gb[which(raindata[,"robs"]==0)]*runif(length(pit_gb[which(raindata[,"robs"]==0)]),0,1)
    hist(pit_gb)
    
    # EMOS
    #hist(pnorm(raindata[,"robs"], mi_mu, mi_sigma))
    pit_mi <- pnorm(raindata[,"robs"], mi_mu, mi_sigma)
    pit_mi[which(raindata[,"robs"]==0)] <- pit_mi[which(raindata[,"robs"]==0)]*runif(length(pit_mi[which(raindata[,"robs"]==0)]),0,1)
    hist(pit_mi)
  }
  
}





###################
# too many missing values in the following variables
if(FALSE){
  # meridional windspeed on 700 hPA
  which(is.na(raindata$v700_mean_mean))
  which(is.na(raindata$v700_mean_min))
  which(is.na(raindata$v700_mean_max))
  which(is.na(raindata$v700_sprd_mean))
  which(is.na(raindata$v700_sprd_min))
  which(is.na(raindata$v700_sprd_max))
  
  # meridional windspeed on 500 hPA
  which(is.na(raindata$v500_mean_mean))
  which(is.na(raindata$v500_mean_min))
  which(is.na(raindata$v500_mean_max))
  which(is.na(raindata$v500_sprd_mean))
  which(is.na(raindata$v500_sprd_min))
  which(is.na(raindata$v500_sprd_max))
  
  # vertical wind speed (upwards) on 850 hPA
  which(is.na(raindata$w850_mean_mean))
  which(is.na(raindata$w850_mean_min))
  which(is.na(raindata$w850_mean_max))
  which(is.na(raindata$w850_sprd_mean))
  which(is.na(raindata$w850_sprd_min))
  which(is.na(raindata$w850_sprd_max))
}










############################################################
# cross validation (10x10)


seedconst <- 7
nrep_cross <- 10
#grid to find optimal mstop for method="noncyclic" in gamboostLSS
grid <- make.grid(max = 500, min = 10, length.out = 10)

rainres <- mclapply(1:nrep_cross,
                    function(k){
                      
                      set.seed(seedconst * k)
                    
                      crps_dt <- crps_df <- crps_g <- crps_gb <- crps_mi <- crps_ml <- crps_mq <- numeric(length = 10L)
                      # crps_b <- numeric(length = 10L)
                      rmse_dt <- rmse_df <- rmse_g <- rmse_gb <- rmse_mi <- rmse_ml <- rmse_mq <- numeric(length = 10L)
                      # rmse_b <- numeric(length = 10L)
                      ll_dt <- ll_df <- ll_g <- ll_gb <- ll_mi <- ll_ml <- ll_mq <- numeric(length = 10L)
                      # ll_b <- numeric(length = 10L)
                      cvr_opt <- numeric(length = 10L)
                      g_error_seed <- numeric(length = 10L)
                      evaltime_user <- matrix(0, ncol = 8, nrow = 10L)
                      evaltime_system <- matrix(0, ncol = 8, nrow = 10L)
                      evaltime_elapsed <- matrix(0, ncol = 8, nrow = 10L)
                      
                      # randomly split data in 10 roughly equal parts
                      # (here: 765 observations -> 7 groups of 76 and 3 groups of 77 obs)
                      id <- c(1:NROW(raindata))
                      testid <- list()
                      for(i in 1:7){
                        testid[[i]] <- sample(id, 76, replace = FALSE)
                        id <- id[!(id %in% testid[[i]])]
                      }
                      for(i in 8:10){
                        testid[[i]] <- sample(id, 77, replace = FALSE)
                        id <- id[!(id %in% testid[[i]])]
                      }
                      
                      
                      
                      for(i in 1:10) {
                        
                        testdata <- raindata[testid[[i]], ]
                        learndata <- raindata[-testid[[i]], ]
                        
                        
                        dt_time <- system.time(dt <- disttree(dt.formula, 
                                                              data = learndata, family = dist_list_cens_normal, 
                                                              censtype = "left", censpoint = 0, type.tree = "ctree"))
                        
                        ## FIX ME: error if type.tree = "mob"
                        
                        #df <- df
                        df_time <- system.time(df <- distforest(df.formula, 
                                                                data = learndata, family = dist_list_cens_normal, ntree = 200, #mtry = 40, fitted.OOB = FALSE,
                                                                censtype = "left", censpoint = 0, type.tree = "ctree"))
                        
                        g_learndata <- learndata
                        g_learndata$robs <- Surv(g_learndata$robs, g_learndata$robs>0, type="left")
                        
                        #gb <- gb
                        gb_time <- system.time(gb <- gamboostLSS(formula = list(mu = gb.mu.formula, sigma = gb.sigma.formula), data = g_learndata,
                                                                 families = as.families(fname = cens("NO", type = "left")), method = "noncyclic",
                                                                 control = boost_control(mstop = 500L)))
                        gb_cvr_time <- system.time(cvr <- cvrisk(gb, grid = grid))
                        mstop(gb) <- mstop(cvr)
                        cvr_opt[i] <- mstop(cvr) 
                        
                        #b <- b
                        #b <- bamlss(list(b.mu.formula, b.sigma.formula), family = "cnorm", 
                        #            data = learndata, sampler = FALSE, optimizer = boost, 
                        #            stop.criterion = "BIC", plot = FALSE)
                        
                        g_time <- system.time(g <- try(gamlss(formula = g.mu.formula, sigma.formula = g.sigma.formula, data = g_learndata, 
                                                          family = cens("NO", type = "left"))))
                        if(inherits(g, "try-error")) {
                          g_time <- NA
                          g <- NA
                          g_error_seed[i] <- seedconst * k
                        }
                        
                        mi_time <- system.time(mi <- crch(formula = robs ~ tppow_mean | tppow_sprd, 
                                                          data = raindata, dist = "gaussian", left = 0, link.scale = "identity"))
                        
                        ml_time <- system.time(ml <- crch(formula = robs ~ tppow_mean | log(tppow_sprd + 0.001), 
                                                          data = raindata, dist = "gaussian", left = 0, link.scale = "log"))
                        
                        mq_time <- system.time(mq <- crch(formula = robs ~ tppow_mean | I(tppow_sprd^2), 
                                                          data = raindata, dist = "gaussian", left = 0, link.scale = "quadratic"))
                        
                        
                        
                        ## get predicted parameter
                        
                        # disttree
                        if(is.vector(coef(dt))){
                          pdt <- t(as.data.frame(coef(dt)))[as.vector(predict(dt, newdata = testdata, type = "node")),]
                          rownames(pdt) <- paste(c(1 : NROW(pdt)))
                        } else {
                          pdt <- coef(dt)[paste(as.vector(predict(dt, newdata = testdata, type = "node"))),]
                        }
                        dt_mu <- as.vector(pdt[,"mu"]) 
                        dt_sigma <- as.vector(pdt[,"sigma"]) 
                        dt_exp <- pnorm(dt_mu/dt_sigma) * (dt_mu + dt_sigma * (dnorm(dt_mu/dt_sigma) / pnorm(dt_mu/dt_sigma)))
                        
                        
                        # distforest
                        pdf <- predict(df, newdata = testdata, type = "parameter")
                        df_mu <- pdf$mu
                        df_sigma <- pdf$sigma
                        df_exp <- pnorm(df_mu/df_sigma) * (df_mu + df_sigma * (dnorm(df_mu/df_sigma) / pnorm(df_mu/df_sigma)))
                        
                        
                        # gamboostLSS
                        pgb <- predict(gb, newdata = testdata, parameter = list("mu","sigma"), type = "response")
                        gb_mu <- pgb[[1]]
                        gb_sigma <- pgb[[2]]
                        gb_exp <- pnorm(gb_mu/gb_sigma) * (gb_mu + gb_sigma * (dnorm(gb_mu/gb_sigma) / pnorm(gb_mu/gb_sigma)))
                        
                        
                        # bamlss
                        #pb <- predict(b, newdata = testdata, type = "parameter")
                        #b_mu <- pb$mu
                        #b_sigma <- pb$sigma
                        #b_exp <- pnorm(b_mu/b_sigma) * (b_mu + b_sigma * (dnorm(b_mu/b_sigma) / pnorm(b_mu/b_sigma)))
                        
                        
                        # gamlss
                        if(!(all(is.na(g)))){
                          g_mu <- try(predict(g, newdata = testdata, what = "mu", type = "response", data = g_learndata))
                          g_sigma <- try(predict(g, newdata = testdata, what = "sigma", type = "response", data = g_learndata))
                          if(inherits(g_mu, "try-error") | inherits(g_sigma, "try-error")) {
                            g_mu <- g_sigma <- g_exp <- NA
                            g_error_seed[i] <- seedconst * k
                          } else g_exp <- pnorm(g_mu/g_sigma) * (g_mu + g_sigma * (dnorm(g_mu/g_sigma) / pnorm(g_mu/g_sigma)))
                        }
                        
                        
                        # EMOS
                        mi_mu <- predict(mi, type = "location", newdata = testdata)     # returns parameter on response scale
                        mi_sigma <- predict(mi, type = "scale", newdata = testdata)
                        mi_exp <- pnorm(mi_mu/mi_sigma) * (mi_mu + mi_sigma * (dnorm(mi_mu/mi_sigma) / pnorm(mi_mu/mi_sigma)))
                        
                        ml_mu <- predict(ml, type = "location", newdata = testdata)     # returns parameter on response scale
                        ml_sigma <- predict(ml, type = "scale", newdata = testdata)
                        ml_exp <- pnorm(ml_mu/ml_sigma) * (ml_mu + ml_sigma * (dnorm(ml_mu/ml_sigma) / pnorm(ml_mu/ml_sigma)))
                        
                        mq_mu <- predict(mq, type = "location", newdata = testdata)     # returns parameter on response scale
                        mq_sigma <- predict(mq, type = "scale", newdata = testdata)
                        mq_exp <- pnorm(mq_mu/mq_sigma) * (mq_mu + mq_sigma * (dnorm(mq_mu/mq_sigma) / pnorm(mq_mu/mq_sigma)))
                        
                        g_error <- any(c(all(is.na(g)), is.na(g_mu), is.na(g_sigma)))
                        
                        # CPRS
                        crps_dt[i] <- sum(crps_cnorm(testdata$robs, location = dt_mu, scale = dt_sigma, lower = 0, upper = Inf))
                        crps_df[i] <- sum(crps_cnorm(testdata$robs, location = df_mu, scale = df_sigma, lower = 0, upper = Inf))
                        crps_g[i] <- if(!g_error) sum(crps_cnorm(testdata$robs, location = g_mu, scale = g_sigma, lower = 0, upper = Inf)) else NA
                        #crps_b[i] <- sum(crps_cnorm(testdata$robs, location = b_mu, scale = b_sigma, lower = 0, upper = Inf))
                        crps_gb[i] <- sum(crps_cnorm(testdata$robs, location = gb_mu, scale = gb_sigma, lower = 0, upper = Inf))
                        crps_mi[i] <- sum(crps_cnorm(testdata$robs, location = mi_mu, scale = mi_sigma, lower = 0, upper = Inf))
                        crps_ml[i] <- sum(crps_cnorm(testdata$robs, location = ml_mu, scale = ml_sigma, lower = 0, upper = Inf))
                        crps_mq[i] <- sum(crps_cnorm(testdata$robs, location = mq_mu, scale = mq_sigma, lower = 0, upper = Inf))
                        
                        # RMSE
                        rmse_dt[i] <- sqrt(mean((dt_exp - testdata[,"robs"])^2))
                        rmse_df[i] <- sqrt(mean((df_exp - testdata[,"robs"])^2))
                        rmse_g[i] <- if(!g_error) sqrt(mean((g_exp - testdata[,"robs"])^2)) else NA
                        rmse_gb[i] <- sqrt(mean((gb_exp - testdata[,"robs"])^2))
                        #rmse_b[i] <- sqrt(mean((b_exp - testdata[,"robs"])^2))
                        rmse_mi[i] <- sqrt(mean((mi_exp - testdata[,"robs"])^2))
                        rmse_ml[i] <- sqrt(mean((ml_exp - testdata[,"robs"])^2))
                        rmse_mq[i] <- sqrt(mean((mq_exp - testdata[,"robs"])^2))
                        
                        # loglikelihood
                        dtll <- dfll <- gll <-  gbll <- mill <- mlll <- mqll <- 0
                        # bll <- 0
                        for(j in 1:(nrow(testdata))){
                          
                          eta_dt <- as.numeric(dist_list_cens_normal$linkfun(cbind(dt_mu, dt_sigma)[j,]))
                          eta_df <- as.numeric(dist_list_cens_normal$linkfun(cbind(df_mu, df_sigma)[j,]))
                          eta_g <- if(!g_error) as.numeric(dist_list_cens_normal$linkfun(cbind(g_mu, g_sigma)[j,])) else NA
                          eta_gb <- as.numeric(dist_list_cens_normal$linkfun(cbind(gb_mu, gb_sigma)[j,]))
                          #eta_b <- as.numeric(dist_list_cens_normal$linkfun(cbind(b_mu, b_sigma)[j,]))
                          eta_mi <- as.numeric(dist_list_cens_normal$linkfun(cbind(mi_mu, mi_sigma)[j,]))
                          eta_ml <- as.numeric(dist_list_cens_normal$linkfun(cbind(ml_mu, ml_sigma)[j,]))
                          eta_mq <- as.numeric(dist_list_cens_normal$linkfun(cbind(mq_mu, mq_sigma)[j,]))
                          
                          dtll_j <- dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_dt, log=TRUE)
                          dfll_j <- dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_df, log=TRUE)
                          gll_j <- if(!g_error) dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_g, log=TRUE) else NA
                          gbll_j <- dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_gb, log=TRUE)
                          #bll_j <- dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_b, log=TRUE)
                          mill_j <- dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_mi, log=TRUE)
                          mlll_j <- dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_ml, log=TRUE)
                          mqll_j <- dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_mq, log=TRUE)
                          
                          dtll <- if(is.na(dtll_j)) {
                            print(eta_dt, testdata[j,"robs"]) 
                            dtll + (-5)
                          } else {dtll + dtll_j}
                          
                          dfll <- if(is.na(dfll_j)) {
                            print(eta_df, testdata[j,"robs"]) 
                            dfll + (-5)
                          } else {dfll + dfll_j}    ## FIX ME: NAs from distforest
                          
                          gll <- if(!g_error) {
                            if(is.na(gll_j)) {
                              print(eta_g, testdata[j,"robs"]) 
                              gll + (-5)
                            } else {gll + gll_j}
                          } else NA
                          
                          gbll <- if(is.na(gbll_j)) {
                            print(eta_gb, testdata[j,"robs"]) 
                            gbll + (-5)
                          } else {gbll + gbll_j}
                          
                          #bll <- if(is.na(bll_j)) {
                          #  print(eta_b, testdata[j,"robs"]) 
                          #  bll + (-5)
                          #} else {bll + bll_j}   
                          
                          mill <- if(is.na(mill_j)) {
                            print(eta_mi, testdata[j,"robs"]) 
                            mill + (-5)
                          } else {mill + mill_j} 
                          
                          mlll <- if(is.na(mlll_j)) {
                            print(eta_ml, testdata[j,"robs"]) 
                            mlll + (-5)
                          } else {mlll + mlll_j} 
                          
                          mqll <- if(is.na(mqll_j)) {
                            print(eta_mq, testdata[j,"robs"]) 
                            mqll + (-5)
                          } else {mqll + mqll_j} 
                        }
                        
                        ll_dt[i] <- dtll
                        ll_df[i] <- dfll
                        ll_g[i] <- if(!g_error) gll else NA
                        ll_gb[i] <- gbll
                        #ll_b[i] <- bll
                        ll_mi[i] <- mill
                        ll_ml[i] <- mlll
                        ll_mq[i] <- mqll
                        
                        evaltime_user[i,] <- c(dt_time[1], df_time[1], g_time[1], gb_time[1], gb_cvr_time[1], mi_time[1], ml_time[1], mq_time[1])
                        evaltime_system[i,] <- c(dt_time[2], df_time[2], g_time[2], gb_time[2], gb_cvr_time[2], mi_time[2], ml_time[2], mq_time[2])
                        evaltime_elapsed[i,] <- c(dt_time[3], df_time[3], g_time[3], gb_time[3], gb_cvr_time[3], mi_time[3], ml_time[3], mq_time[3])
                        
                      }
                      
                      res <- list()
                      res$ll <- cbind(ll_dt, ll_df, ll_g, ll_gb, ll_mi, ll_ml, ll_mq)
                      res$rmse <- cbind(rmse_dt, rmse_df, rmse_g, rmse_gb, rmse_mi, rmse_ml, rmse_mq)
                      res$crps <- cbind(crps_dt, crps_df, crps_g, crps_gb, crps_mi, crps_ml, crps_mq)
                      #res$ll <- cbind(ll_dt, ll_df, ll_g, ll_gb, ll_b, ll_mi, ll_ml, ll_mq)
                      #res$rmse <- cbind(rmse_dt, rmse_df, rmse_g, rmse_gb, rmse_b, rmse_mi, rmse_ml, rmse_mq)
                      #res$crps <- cbind(crps_dt, crps_df, crps_g, crps_gb, crps_b, crps_mi, crps_ml, crps_mq)
                      colnames(res$rmse) <- colnames(res$ll) <- colnames(res$crps) <- 
                        c("disttree", "distforest", "gamlss", "gamboostLSS", "EMOS id", "EMOS log", "EMOS quad")
                       #c("disttree", "distforest", "gamlss", "gamboostLSS", "bamlss", "EMOS id", "EMOS log", "EMOS quad")
                      res$cvr_opt <- cvr_opt
                      res$evaltime_user <- evaltime_user
                      res$evaltime_system <- evaltime_system
                      res$evaltime_elapsed <- evaltime_elapsed
                      colnames(res$evaltime_user) <- colnames(res$evaltime_system) <- colnames(res$evaltime_elapsed) <- 
                        c("disttree", "distforest", "gamlss", "gamboostLSS", "gamboostLSS_cvrisk", "EMOS id", "EMOS log", "EMOS quad")
                      res$g_error_seed <- g_error_seed
                      
                      return(res)
                    },
                    mc.cores = detectCores() - 1
)


#save(rainres, file = "~/svn/partykit/pkg/disttree/inst/draft/rain_axams7.rda")
save(rainres, file = "~/disttree/inst/draft/rain_axams7.rda")
  



if(FALSE){
  rain_rmse <- matrix(0, ncol = ncol(rainres[[1]]$rmse), length(rainres))
  rain_ll <- matrix(0, ncol = ncol(rainres[[1]]$ll), length(rainres))
  rain_crps <- matrix(0, ncol = ncol(rainres[[1]]$crps), length(rainres))
  
  for(i in 1:length(rainres)){
    rain_rmse[i,] <- colMeans(rainres[[i]]$rmse)
    rain_ll[i,] <- colMeans(rainres[[i]]$ll)
    rain_crps[i,] <- colMeans(rainres[[i]]$crps)
  }
  
  colnames(rain_rmse) <- colnames(rain_ll) <- colnames(rain_crps) <-
    colnames(rainres[[1]]$rmse)
  
  boxplot(rain_rmse)
  boxplot(rain_ll)
  boxplot(rain_crps)
  
  # scoring rules
  
  # skills score mit Referenz-Modell:
  # rain_crps / rain_crps[,8]
  # boxplot(rain_crps / rain_crps[,8])
  # boxplot(rain_crps[,-5] / rain_crps[,8])
  # boxplot(1 - rain_crps[,-5] / rain_crps[,8])
  # abline(h = 0, col = 2, lwd = 2)
}




