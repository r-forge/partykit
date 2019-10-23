# -------------------------------------------------------------------
# - NAME:   circforest_fit_models.R
# - AUTHOR: Moritz N. Lang, Lisa Schlosser
# - DATE:   2019-09-13
# -------------------------------------------------------------------
# - PURPOSE:
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2019-10-23 on thinkmoritz
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Preliminaries
# -------------------------------------------------------------------
## Set time to utc (just in case)
Sys.setenv('TZ'='UTC')

## Load packages
library("optparse")
library("zoo")
library("circtree")
library("disttree")
library("verification") # Careful, version 1.35 needed!!
library("circglmbayes")

## Create folder for outputs
if (! dir.exists("results")) dir.create("results")

## Small helper functions

## Exponential weights
#calc_expweights <- function(theta, n){rev(sapply(1:n, function(x) theta^(x -1) * (1 - theta)))}
calc_expweights <- function(theta, n){rev(sapply(1:n, function(x) theta^(x -1) / 
  ((theta^n - 1) /(theta - 1))))}

## Exponential centered weights
calc_expweights_centered <- function(theta, n){sapply(-n:n, function(x) theta^(abs(x) -1) / 
  sum(theta^(abs(-n:n) -1)))}

ddff2uv <- function (dd, ff = NULL) {
    if (class(dd) %in% c("zoo", "data.frame")) {
        if (sum(!c("ff", "dd") %in% names(dd)) > 0) 
            stop("necessary colums \"ff\" and/or \"dd\" missing")
        ff = as.numeric(dd$ff)
        dd = as.numeric(dd$dd)
        print("x")
    }
    else if (NCOL(dd) == 2) {
        ff <- dd[, 1]
        dd <- dd[, 2]
    }
    metrad <- dd * pi/180
    u <- ff * (-sin(metrad))
    v <- ff * (-cos(metrad))
    rad <- 2 * pi - metrad - pi/2
    rad <- ifelse(rad >= 2 * pi, rad - 2 * pi, rad)
    data.frame(u, v, rad)
}


# -------------------------------------------------------------------
# NAMELIST PARAMETERS
# -------------------------------------------------------------------
option_list <- list(
  make_option(c("-v", "--verbose"), action = "store_true", default = TRUE,
    help = "Print extra output [default]"),
  make_option(c("-q", "--quietly"), action = "store_false",
    dest = "verbose", help = "Print little output"),
  make_option("--run_name", type = "character", default = "v14",
    help = "Run name or version of script used for output name [default \"%default\"]"),
  make_option("--station", type = "character", default = "ibk",
    help = "Weather Station used for fitting (e.g., 'ibk', 'vie') [default \"%default\"]"),
  make_option("--lag", type = "integer", default = 1,
    help = "Lag in 10min time steps [default %default]"),
  make_option("--lowff", action = "store_true", default = FALSE,
    help = "Use wind speeds (ff) below 1 m/s (Warning: TRUE not meaningful for Vie) [default]"),
  make_option("--plot", action = "store_true", default = FALSE,
    help = "Plot validation [default]"),
  make_option(c("--seed"), type = "integer", default = 123,
    help="Set seed for bootstrapping [default %default]")
  )

opt <- parse_args(OptionParser(option_list = option_list))

# -------------------------------------------------------------------
# Pre-process data
# -------------------------------------------------------------------
## Load data
tmp <- readRDS(sprintf("data/circforest_prepared_data_%s_lag%sh_2014-01-01_2018-12-31_update20191023.rds", opt$station, opt$lag))

if (opt$station == "ibk") {
  ## Subset data to five years and to full hours 
  tmp <- window(tmp, start = "2014-01-01", end = "2018-12-31")
  d <- tmp[as.POSIXlt(index(tmp))$min == 0L, ]; rm(tmp); gc()

  ## Remove very low values of ff 
  if(!opt$lowff){
    d <- d[!d[, "response.ff"] < 1,] ## ff below 1 not very meaningful, but at least correctly measured
  }

  ## Remove covariates with more than 10 nans
  idx_names <- names(which(apply(d, 2, function(x) sum(is.na(x))) > nrow(d) / 100 * 5)) ## remove sattelberg and ellboegen
  d <- d[, ! (names(d) %in% idx_names)]

  ## Omit nans
  d <- na.omit(d)

  ## Calculate u and v for lagged response (ibk)
  tmp <- ddff2uv(dd = as.numeric(d$innsbruck.dd), ff = as.numeric(d$innsbruck.ff))
  d$innsbruck.u <- tmp$u
  d$innsbruck.v <- tmp$v

  ## Transform response.dd from 0-360 degree to [-pi, pi]
  d$response.dd <- d$response.dd / 360 * 2*pi
  d$response.dd[d$response.dd > pi] <- d$response.dd[d$response.dd > pi] - 2*pi

} else if (opt$station == "vie") {
  ## Subset data to five years and to full hours
  tmp <- window(tmp, start = "2014-01-01", end = "2018-12-31")
  d <- tmp[as.POSIXlt(index(tmp))$min == 0L, ]; rm(tmp); gc()

  ## Remove very low values of ff
  if(!opt$lowff){
    d <- d[!d[, "response.ff"] < 1,] ## ff == 0.5 -> dd = NA; ff == 0 -> dd = 0
  }

  ## Remove covariates with more than 10 nans
  idx_names <- names(which(apply(d, 2, function(x) sum(is.na(x))) > nrow(d) / 100 * 5))
  d <- d[, ! (names(d) %in% idx_names)]

  ## Omit nans
  d <- na.omit(d)

  ## Calculate u and v for lagged response (windmessanlage 34)
  tmp <- ddff2uv(dd = as.numeric(d$wien_schwechat_flughafen.dd), ff = as.numeric(d$wien_schwechat_flughafen.ff))
  d$wien_schwechat_flughafen.u <- tmp$u
  d$wien_schwechat_flughafen.v <- tmp$v

  ## Transform response.dd from 0-360 degree to [-pi, pi]
  d$response.dd <- d$response.dd / 360 * 2*pi
  d$response.dd[d$response.dd > pi] <- d$response.dd[d$response.dd > pi] - 2*pi

} else {
  stop("Station not supported, currently station must be 'vie' or 'ibk'...")
} 


# -------------------------------------------------------------------
# Fit persistence model with previous time points (no CV)
# -------------------------------------------------------------------
## Calculate exponential weights
exp_weights <- calc_expweights(0.5, 6)

pred_pers_hour <- data.frame(mu = rep(NA, nrow(d)), kappa = rep(NA, nrow(d)))
for (i in seq(1:nrow(d))) {
  if (opt$verbose) cat(sprintf("Fitting persistence (previous time points) %2d%%\r", (i/nrow(d) * 100)%/%5 * 5))

  i_plt <- as.POSIXlt(index(d)[i], origin = "1970-01-01")
  
  ## Skip 29th of February 
  if ((i_plt$mon + 1 == 2) & (i_plt$mday == 29)) {
    pred_pers_hour[i, ] <- c("mu" = NA, "kappa" = NA)
    next
  }

  ## Subset training data set
  idx <- as.POSIXct(sprintf("%04d-%02d-%02d %02d:%02d:00", i_plt$year + 1900, i_plt$mon + 1, i_plt$mday, 
    i_plt$hour, i_plt$min), origin = "1970-01-01") + 60 * 60 * (seq(-5, -0) - opt$lag)

  train <- subset(d, index(d) %in% idx)

  ## Go to next iteration if train has no values
  if(nrow(train) == 0){
    pred_pers_hour[i, ] <- c("mu" = NA, "kappa" = NA)
    next
  }

  ## Get weights and rescale to sum up to one
  train_weights <- exp_weights[idx %in% index(train)]
  train_weights <- train_weights / sum(train_weights) * 1

  ## Perform fit only if values vary
  train_response <- as.numeric(train$response.dd)

  if (!all(train_response == train_response[1])) {
    ## Fit persistency
    pers_hour_fit <- try(distfit(train_response,
                           family = dist_vonmises(), weights = train_weights))

    ## Predict parameters
    if (class(pers_hour_fit) == "try-error") {
      pred_pers_hour[i, ] <- c("mu" = NA, "kappa" = NA)
    } else {
      pred_pers_hour[i, ] <- coef(pers_hour_fit, type = "parameter")
    }

  } else {
    ## Set coefficients to single value and kappa to ~Infinity
    pred_pers_hour[i, ] <- c("mu" = unique(train_response), "kappa" = 1e+16)
  }
  
}

if (opt$verbose) cat("\n")   

## Save predictions
saveRDS(pred_pers_hour, file = sprintf("results/circforest_pred_pers_hour_%s_lag%s_%s%s.rds", 
  opt$station, opt$lag, ifelse(opt$lowff, "with_lowff_", ""), opt$run_name))
rm(pred_pers_hour, pers_hour_fit, train); gc()


# -------------------------------------------------------------------
# Fit persistence model with previous days (no CV)
# -------------------------------------------------------------------
## Calculate exponential weights
exp_weights_centered <- calc_expweights_centered(0.5, 3)

pred_pers_day <- data.frame(mu = rep(NA, nrow(d)), kappa = rep(NA, nrow(d)))
for (i in seq(1:nrow(d))) {
  if (opt$verbose) cat(sprintf("Fitting persistence (previous days) %2d%%\r", (i/nrow(d) * 100)%/%5 * 5))

  i_plt <- as.POSIXlt(index(d)[i], origin = "1970-01-01")
  
  ## Skip 29th of February 
  if ((i_plt$mon + 1 == 2) & (i_plt$mday == 29)) {
    pred_pers_day[i, ] <- c("mu" = NA, "kappa" = NA)
    next
  }

  ## Subset training data set
  idx <- as.POSIXct(sprintf("%04d-%02d-%02d %02d:%02d:00", i_plt$year + 1900, i_plt$mon + 1, i_plt$mday, 
    i_plt$hour, i_plt$min), origin = "1970-01-01") - 60 * 60 * 24 + 60 * 60 * seq(-3, 3)

  train <- subset(d, index(d) %in% idx)

  ## Go to next iteration if train has no values
  if(nrow(train) == 0){
    pred_pers_day[i, ] <- c("mu" = NA, "kappa" = NA)
    next
  }
  
  ## Get weights and rescale to sum up to one
  train_weights <- exp_weights_centered[idx %in% index(train)]
  train_weights <- train_weights / sum(train_weights) * 1

  ## Perform fit only if values vary
  train_response <- as.numeric(train$response.dd)

  if (!all(train_response == train_response[1])) {
    ## Fit persistency
    pers_day_fit <- try(distfit(train_response,
                           family = dist_vonmises(), weights = train_weights))

    ## Predict parameters
    if (class(pers_day_fit) == "try-error") {
      pred_pers_day[i, ] <- c("mu" = NA, "kappa" = NA)
    } else {
      pred_pers_day[i, ] <- coef(pers_day_fit, type = "parameter")
    }

  } else {
    ## Set coefficients to single value and kappa to ~Infinity
    pred_pers_day[i, ] <- c("mu" = unique(train_response), "kappa" = 1e+16)
  }

}

if (opt$verbose) cat("\n")   
  

## Save predictions
saveRDS(pred_pers_day, file = sprintf("results/circforest_pred_pers_day_%s_lag%s_%s%s.rds", 
  opt$station, opt$lag, ifelse(opt$lowff, "with_lowff_", ""), opt$run_name))
rm(pred_pers_day, pers_day_fit, train); gc()


# -------------------------------------------------------------------
# Fit climatology (with CV)
# -------------------------------------------------------------------
## Fit models with cross-validation
#cvID <- sort(rep(1:5, ceiling(nrow(d) / 5)))[1:nrow(d)]
cvID <- as.numeric(factor(as.POSIXlt(index(d))$year))

pred_clim <- lapply(unique(cvID), function(x) data.frame(mu = rep(NA, sum(cvID == x)), 
  kappa = rep(NA, sum(cvID == x))))

for (cv in unique(cvID)) { 
  if (opt$verbose) cat(sprintf("Fitting models cv %s/%s\n", cv, max(cvID)))

  train <- d[cv != cvID, ]
  test <- d[cv == cvID, ]
 
  d_range <- as.POSIXlt(range(index(d)))

  for (i in seq(1:nrow(test))) {

    if (opt$verbose) cat(sprintf("Fitting climatology %2d%%\n", (i/nrow(test) * 100)%/%5 * 5))
  
    i_plt <- as.POSIXlt(index(test)[i], origin = "1970-01-01")
    
    ## Skip 29th of February 
    if ((i_plt$mon + 1 == 2) & (i_plt$mday == 29)) {
      pred_clim[[cv]][i, ] <- c("mu" = NA, "kappa" = NA)
      next
    }
  
    ## Subset training data set
    i_years <- seq.int(min(d_range$year), max(d_range$year)) + 1900

    idx <- lapply(i_years, function(x) 
      as.POSIXct(sprintf("%04d-%02d-%02d %02d:%02d:00", x, i_plt$mon + 1, i_plt$mday, 
      i_plt$hour, i_plt$min), origin = "1970-01-01") + 60 * 60 * 24 * seq(-15, 15))
  
    idx <- do.call(c, idx)
    train_subset <- subset(train, index(train) %in% idx)

    ## Go to next iteration if train has no values
    if(nrow(train_subset) == 0){
      pred_clim[[cv]][i, ] <- c("mu" = NA, "kappa" = NA)
      next
    }
    
    train_response <- as.numeric(train_subset$response.dd)

    ## Perform fit only if values vary
    if (!all(train_response == train_response[1])) {

      ## Fit climatology
      climfit <- try(distfit(train_response,
                             family = dist_vonmises()))

      ## Predict parameters
      if (class(climfit) == "try-error") {
        pred_clim[[cv]][i, ] <- c("mu" = NA, "kappa" = NA)
      } else {
        pred_clim[[cv]][i, ] <- coef(climfit, type = "parameter")
      }

    } else {
      ## Set coefficients to single value and kappa to ~Infinity
      pred_clim[[cv]][i, ] <- c("mu" = unique(train_response), "kappa" = 1e+16)
    }

  }
}

if (opt$verbose) cat("\n")  

## Save predictions
saveRDS(pred_clim, file = sprintf("results/circforest_pred_clim_%s_lag%s_%s%s.rds", 
  opt$station, opt$lag, ifelse(opt$lowff, "with_lowff_", ""), opt$run_name))
rm(pred_clim, climfit, train, test, train_subset); gc()


# -------------------------------------------------------------------
# Fit circular linear model (with CV)
# -------------------------------------------------------------------
## Fit models with cross-validation
#cvID <- sort(rep(1:5, ceiling(nrow(d) / 5)))[1:nrow(d)]
cvID <- as.numeric(factor(as.POSIXlt(index(d))$year))

pred_lm <- lapply(unique(cvID), function(x) data.frame(mu = rep(NA, sum(cvID == x)), 
  kappa = rep(NA, sum(cvID == x))))

d.lm <- d
d.lm[d.lm$response.dd < 0] <- d.lm$response.dd[d.lm$response.dd < 0] + 2 * pi

for (cv in unique(cvID)) { 
  if (opt$verbose) cat(sprintf("Fitting models cv %s/%s\n", cv, max(cvID)))

  train <- d.lm[cv != cvID, ]
  test <- d.lm[cv == cvID, ]
 
  d_range <- as.POSIXlt(range(index(d)))

  for (i in seq(1:nrow(test))) {

    if (opt$verbose) cat(sprintf("Fitting circular linear model %2d%%\n", (i/nrow(test) * 100)%/%5 * 5))
  
    i_plt <- as.POSIXlt(index(test)[i], origin = "1970-01-01")
    
    ## Skip 29th of February 
    if ((i_plt$mon + 1 == 2) & (i_plt$mday == 29)) {
      pred_lm[[cv]][i, ] <- c("mu" = NA, "kappa" = NA)
      next
    }
  
    ## Subset training data set
    i_years <- seq.int(min(d_range$year), max(d_range$year)) + 1900

    idx <- lapply(i_years, function(x) 
      as.POSIXct(sprintf("%04d-%02d-%02d %02d:%02d:00", x, i_plt$mon + 1, i_plt$mday, 
      i_plt$hour, i_plt$min), origin = "1970-01-01") + 60 * 60 * 24 * seq(-15, 15))
  
    idx <- do.call(c, idx)
    train_subset <- subset(train, index(train) %in% idx)

    ## Go to next iteration if train has no values
    if(nrow(train_subset) == 0){
      pred_lm[[cv]][i, ] <- c("mu" = NA, "kappa" = NA)
      next
    }
  
    ## Fit lm model
    if(opt$station == "ibk"){
      f.lm <- as.formula(response.dd ~ innsbruck.u + innsbruck.v + innsbruck.ff)
    } else if (opt$station == "vie"){
      f.lm <- as.formula(response.dd ~ wien_schwechat_flughafen.u + wien_schwechat_flughafen.v + wien_schwechat_flughafen.ff)
    } else {
      stop("Station not supported, currently station must be 'vie' or 'ibk'...")
    }
    lmfit <- try(circGLM(f.lm, data = train_subset, skipDichSplit = FALSE))

    ## Predict parameters
    if (any(class(lmfit) %in% "try-error")) {
      pred_lm[[cv]][i, ] <- c("mu" = NA, "kappa" = NA)
    } else {
      pred_tmp <- predict(lmfit, newdata = subset(d.lm, index(d.lm) %in% index(test)[i]))
      pred_tmp <- pred_tmp %% (2*pi)
      pred_tmp[pred_tmp > pi] <- 
        pred_tmp[pred_tmp > pi] - 2 * pi
      pred_lm[[cv]][i, "mu"] <- pred_tmp

      pred_lm[[cv]][i, "kappa"] <- coef(lmfit)["Kappa", "Estimate"]
    }
  }
}

if (opt$verbose) cat("\n")  

## Save predictions
saveRDS(pred_lm, file = sprintf("results/circforest_pred_lm_%s_lag%s_%s%s.rds", 
  opt$station, opt$lag, ifelse(opt$lowff, "with_lowff_", ""), opt$run_name))
rm(pred_lm, lmfit, train, test, train_subset); gc()


# -------------------------------------------------------------------
# Fit tree and forest
# -------------------------------------------------------------------
## Transform zoo object to data.frame
timepoints <- index(d)
tmp_time <- as.POSIXlt(index(d))
d <- as.data.frame(d)

## Add day of the year and time
d$daytime <- cut(tmp_time$hour, seq(0,24,by=3), include.lowest = TRUE)
d$doy <- tmp_time$yday + 1

rm(tmp_time); gc()

## Set up formula
f <- as.formula(paste("response.dd ~ ", paste(names(d)[-grep("response", names(d))], collapse= "+")))

## Fit single tree for visualization
m_ct.plot <- circtree(formula = f,
                 data = d,
                 control = disttree_control(mincriterion = (1 - .Machine$double.eps),
                                            minbucket = 2000))

saveRDS(m_ct.plot, file = sprintf("results/circforest_model_tree_%s_lag%s_%s%s_4plotting.rds", 
  opt$station, opt$lag, ifelse(opt$lowff, "with_lowff_", ""), opt$run_name))


## Fit models with cross-validation
#cvID <- sort(rep(1:5, ceiling(nrow(d) / 5)))[1:nrow(d)]
vip_cf <- list()

for(cv in unique(cvID)) { 
  if (opt$verbose) cat(sprintf("Fitting models cv %s/%s\n", cv, max(cvID)))

  train <- d[cv != cvID, ]
  test <- d[cv == cvID, ]
  
  ## Fit tree
  m_ct <- circtree(formula = f,
                   data = train,
                   control = disttree_control(mincriterion = (1 - .Machine$double.eps),
                                              minbucket = 2000))

  # Fit forest
  m_cf <- circforest(formula = f,
                     data = train,
                     ntree = 100,
                     mtry = ceiling(1. * length(all.vars(f[[3]]))),
                     perturb = list(replace = FALSE, fraction = 0.3),
                     control = disttree_control(nmax = c("yx" = Inf, "z" = 50)))

  #vip_cf[[cv]] <- disttree:::varimp.distforest(m_cf, nperm = 10)

  ## Predict models
  pred_ct.tmp <- predict(m_ct, newdata = test, type = "parameter")
  pred_cf.tmp <- predict(m_cf, newdata = test, type = "parameter")

  ## Save predictions (and models)
  saveRDS(m_ct, 
    file = sprintf("results/circforest_model_tree_%s_lag%s_%s%s_cv%s.rds", 
      opt$station, opt$lag, ifelse(opt$lowff, "with_lowff_", ""), opt$run_name, cv))
  saveRDS(pred_ct.tmp, 
    file = sprintf("results/circforest_pred_tree_%s_lag%s_%s%s_cv%s.rds",
       opt$station, opt$lag, ifelse(opt$lowff, "with_lowff_", ""), opt$run_name, cv))
  saveRDS(pred_cf.tmp, 
    file = sprintf("results/circforest_pred_forest_%s_lag%s_%s%s_cv%s.rds",
      opt$station, opt$lag, ifelse(opt$lowff, "with_lowff_", ""), opt$run_name, cv))
  #saveRDS(vip_cf, 
  #  file = sprintf("results/circforest_vip_forest_%s_lag%s_%s_cv%s.rds",
  #    opt$station, opt$lag, opt$run_name, cv))
  rm(pred_ct.tmp, pred_cf.tmp, m_cf, m_ct, train, test); gc()
}


# -------------------------------------------------------------------
# Combine predictions
# -------------------------------------------------------------------
## Load predictions
pred_ct <- pred_cf <- list()
for(cv in unique(cvID)) {
  pred_ct.tmp <- readRDS(file = sprintf("results/circforest_pred_tree_%s_lag%s_%s%s_cv%s.rds",
    opt$station, opt$lag, ifelse(opt$lowff, "with_lowff_", ""), opt$run_name, cv))
  pred_cf.tmp <- readRDS(file = sprintf("results/circforest_pred_forest_%s_lag%s_%s%s_cv%s.rds", 
    opt$station, opt$lag, ifelse(opt$lowff, "with_lowff_", ""), opt$run_name, cv))
  pred_ct[[cv]] <- pred_ct.tmp
  pred_cf[[cv]] <- pred_cf.tmp
  rm(pred_cf.tmp, pred_ct.tmp); gc()
}

pred_pers_hour <- readRDS(file = sprintf("results/circforest_pred_pers_hour_%s_lag%s_%s%s.rds", 
  opt$station, opt$lag, ifelse(opt$lowff, "with_lowff_", ""), opt$run_name))
pred_pers_day <- readRDS(file = sprintf("results/circforest_pred_pers_day_%s_lag%s_%s%s.rds", 
  opt$station, opt$lag, ifelse(opt$lowff, "with_lowff_", ""), opt$run_name))
pred_clim <- readRDS(file = sprintf("results/circforest_pred_clim_%s_lag%s_%s%s.rds", 
  opt$station, opt$lag, ifelse(opt$lowff, "with_lowff_", ""), opt$run_name))
pred_lm <- readRDS(file = sprintf("results/circforest_pred_lm_%s_lag%s_%s%s.rds", 
  opt$station, opt$lag, ifelse(opt$lowff, "with_lowff_", ""), opt$run_name))

## Combine predictions
pred <- list(tree = do.call("rbind", pred_ct), forest = do.call("rbind", pred_cf), 
  climatology = do.call("rbind", pred_clim), linear_model = do.call("rbind", pred_lm),
  persistence_hour = pred_pers_hour, persistence_day = pred_pers_day)
obs  <- d$response.dd

## Remove nans (if any)
idx <- do.call(rbind, sapply(pred, function(x) which(is.na(x) | x$mu < -pi | x$mu > pi, arr.ind=TRUE)))
idx <- unique(idx[,1]) 
if (length(idx) > 0) {
  pred_naomit <- lapply(pred, function(x) x[-idx, ])
  obs_naomit  <- obs[-idx]
  timepoints_naomit <- timepoints[-idx]
} else {
  pred_naomit <- pred
  obs_naomit  <- obs
  timepoints <- timepoints
}

## Save results
save(pred, pred_naomit, obs, obs_naomit, timepoints, timepoints_naomit,
  file = sprintf("results/circforest_results_%s_lag%s_%s%s.rda", opt$station, opt$lag, ifelse(opt$lowff, "with_lowff_", ""), 
    opt$run_name))


# -------------------------------------------------------------------
# Validate predictions based on grimit's crps
# -------------------------------------------------------------------
load(file = sprintf("results/circforest_results_%s_lag%s_%s%s.rda", opt$station, opt$lag, ifelse(opt$lowff, "with_lowff_", ""), 
  opt$run_name))

## Validate models based on grimit's crps
crps <- lapply(pred_naomit, function(x) sapply(1:nrow(x), function(i)
  as.numeric(crps.circ(x = obs_naomit[i], mu = x[i, "mu"], kappa = x[i, "kappa"]))))
crps <- data.frame(do.call(cbind, crps))

## Calculate mean scores per hour over single month
crps.agg <- aggregate(crps, list(month = as.POSIXlt(timepoints_naomit)$mon + 1, 
  hour = as.POSIXlt(timepoints_naomit)$hour + 1), FUN = mean)

crps.agg <- within(crps.agg, {
  month = factor(month)
  hour = factor(hour)
  })

## Calculate boot-strapped mean values
set.seed(opt$seed)
kboot <- 500
crps.boot <- data.frame(matrix(NA, ncol = ncol(crps), nrow = kboot, dimnames = list(NULL, names(crps))))
for (i in 1:kboot) {
   s <- sample(1 : nrow(crps), nrow(crps), replace = TRUE)
   crps.boot[i,] <- apply(coredata(crps)[s, ], 2, mean)
}

## Calculate skill scores
crps_skill <- (1 - crps[, c("climatology", "persistence_hour", "persistence_day", 
  "linear_model", "tree", "forest")] / crps[, "climatology"]) * 100

crps.agg_skill <- cbind(crps.agg[, c("month", "hour")], 
  (1 - crps.agg[, c("climatology", "persistence_hour", "persistence_day", "linear_model", "tree", 
  "forest")] / crps.agg[, "climatology"]) * 100)

crps.boot_skill <- (1 - crps.boot[, c("climatology", "persistence_hour", "persistence_day", 
  "linear_model", "tree", "forest")] / crps.boot[, "climatology"]) * 100

## Save validation
save(crps, crps.boot, crps_skill, crps.boot_skill, crps.agg, crps.agg_skill,
  file = sprintf("results/circforest_validation_%s_lag%s_%s%s.rda", opt$station, opt$lag, ifelse(opt$lowff, "with_lowff_", ""), 
  opt$run_name))


## -------------------------------------------------------------------
## Plotting
## -------------------------------------------------------------------
if (opt$plot) {

  ## Plot crps skill scores
  load(file = sprintf("results/circforest_validation_%s_lag%s_%s%s.rda", opt$station, opt$lag, ifelse(opt$lowff, "with_lowff_", ""), 
    opt$run_name))

  ## Plot crps scores
  X11(width = 8, height = 4.5)
  par(mar = c(3.1, 4.1, 2.1, 2.1))
  boxplot(crps[, c("climatology", "persistence", "tree", "forest")], ylim = c(0, 1), col = gray(0.6),
    ylab = "CRPS [rad]", main = "Raw values")
  text(c(mean(unlist(crps["climatology"]), na.rm = TRUE),
         mean(unlist(crps["persistence"]), na.rm = TRUE),
         mean(unlist(crps["tree"]), na.rm = TRUE),
         mean(unlist(crps["forest"]), na.rm = TRUE)), "*", cex=2.5 , col = "red")
  dev.print(pdf, sprintf("results/_plot_circforest_validation_crpsraw_%s_lag%s_%s%s.pdf", 
    opt$station, opt$lag, ifelse(opt$lowff, "with_lowff_", ""), opt$run_name))

  X11(width = 8, height = 4.5)
  par(mar = c(3.1, 4.1, 2.1, 2.1))
  boxplot(crps.boot[, c("climatology", "persistence", "tree", "forest")],
    ylim = c(0, 1), col = gray(0.6), ylab = "CRPS [rad]", main = "Raw values")
  text(c(mean(unlist(crps.boot["climatology"]), na.rm = TRUE),
         mean(unlist(crps.boot["persistence"]), na.rm = TRUE),
         mean(unlist(crps.boot["tree"]), na.rm = TRUE),
         mean(unlist(crps.boot["forest"]), na.rm = TRUE)), "*", cex=2.5 , col = "red")
  dev.print(pdf, sprintf("results/_plot_circforest_validation_crpsraw_boot_%s_lag%s_%s%s.pdf", 
    opt$station, opt$lag, ifelse(opt$lowff, "with_lowff_", ""), opt$run_name))

  X11(width = 8, height = 4.5)
  par(mar = c(3.1, 4.1, 2.1, 2.1))
  boxplot(crps_skill, ylim = c(-60, 110), col = gray(0.6),
    ylab = "CRPS skill ccore [%]", main = "Raw values")
  text(c(mean(unlist(crps_skill["climatology"]), na.rm = TRUE),
         mean(unlist(crps_skill["persistence"]), na.rm = TRUE),
         mean(unlist(crps_skill["tree"]), na.rm = TRUE),
         mean(unlist(crps_skill["forest"]), na.rm = TRUE)), "*", cex=2.5 , col = "red")
  dev.print(pdf, sprintf("results/_plot_circforest_validation_crpsskill_%s_lag%s_%s%s.pdf", 
    opt$station, opt$lag, ifelse(opt$lowff, "with_lowff_", ""), opt$run_name))

  X11(width = 8, height = 4.5)
  par(mar = c(3.1, 4.1, 2.1, 2.1))
  boxplot(crps.boot_skill, ylim = c(-60, 110), col = gray(0.6),
    ylab = "CRPS skill ccore [%]", main = "Boot-strapped mean values")
  dev.print(pdf, sprintf("results/_plot_circforest_validation_crpsskill_boot_%s_lag%s_%s%s.pdf", 
    opt$station, opt$lag, ifelse(opt$lowff, "with_lowff_", ""), opt$run_name))

  library(ggplot2)
  library(reshape2)
  library(ggpubr)
  
  crps.agg_skill.m <- melt(crps.agg_skill, variable.name = "model")
  crps.agg_skill.m$model <- plyr::revalue(crps.agg_skill.m$model,
    c("climatology" = "Climatology",
    "persistence" = "Persistence",
    "tree" = "Tree",
    "forest" = "Forest"))
  
  theme_set(theme_bw(base_size = 12) +
     theme(panel.grid.major = element_line(linetype = "dotted", colour = "grey80"),
           panel.grid.minor = element_blank(),
           plot.title = element_text(hjust = 0.5)))
  
  p1 <- ggplot(crps.agg_skill.m, aes(x = model, y = value)) +
        geom_hline(yintercept = 0, linetype ="solid", colour = "gray80") +
        stat_boxplot(geom = "errorbar", width = 0.2) +
        geom_boxplot(fill = "gray60")
  p1 <- p1 + labs(x = "", y = "CRPS skill score [%]") 

  dev.new(width=8, height=4.5)
  print(ggarrange(p1, legend = "none"))
  pdf_file <- sprintf("results/_plot_circforest_validation_crpsskill_agg_%s_lag%s_%s%s.pdf",
    opt$station, opt$lag, ifelse(opt$lowff, "with_lowff_", ""), opt$run_name)
  ggsave(pdf_file)

  crps.m <- melt(crps[, c("climatology", "persistence", "tree", "forest")], variable.name = "model")
  crps.m$model <- plyr::revalue(crps.m$model,
    c("climatology" = "Climatology",
    "persistence" = "Persistence",
    "tree" = "Tree",
    "forest" = "Forest"))
  
  theme_set(theme_bw(base_size = 12) +
     theme(panel.grid.major = element_line(linetype = "dotted", colour = "grey80"),
           panel.grid.minor = element_blank(),
           plot.title = element_text(hjust = 0.5)))
  
  p2 <- ggplot(crps.m, aes(x = model, y = value)) +
        geom_hline(yintercept = 0, linetype ="solid", colour = "gray80") +
        stat_boxplot(geom = "errorbar", width = 0.2) +
        geom_boxplot(fill = "gray60")
  p2 <- p2 + labs(x = "", y = "CRPS [rad]") 

  dev.new(width=8, height=4.5)
  print(ggarrange(p2, legend = "none"))
  pdf_file <- sprintf("results/_plot_circforest_validation_crpsraw2_%s_lag%s_%s%s.pdf",
    opt$station, opt$lag, ifelse(opt$lowff, "with_lowff_", ""), opt$run_name)
  ggsave(pdf_file)


  ## Plot single tree
  m_ct.plot <- readRDS(m_ct.plot, file = sprintf("results/circforest_model_tree_%s_lag%s_%s%s_4plotting.rds",    
    opt$station, opt$lag, ifelse(opt$lowff, "with_lowff_", ""), opt$run_name))

  X11(width = 18, height = 10)
  par(mar = c(3.1, 4.1, 2.1, 2.1))
  plot(m_ct.plot, ep_args = list(justmin = 10), tp_args = list(type = "response", plot_type = "geographics"), 
    ip_args = list(pval = FALSE))
  dev.print(pdf, sprintf("results/_plot_circforest_exampletree_%s_lag%s_%s%s.pdf", 
    opt$station, opt$lag, ifelse(opt$lowff, "with_lowff_", ""), opt$run_name))
}

