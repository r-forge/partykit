# -------------------------------------------------------------------
# - NAME:   circforest_fit_models_ibk.R
# - AUTHOR: Moritz N. Lang, Lisa Schlosser
# - DATE:   2019-06-03
# -------------------------------------------------------------------
# - PURPOSE:
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2019-06-06 on thinkmoritz
# -------------------------------------------------------------------

## Set time to utc (just in case)
Sys.setenv('TZ'='UTC')

## Load packages
library("zoo")
library("disttree")
library("circmax")
library("verification") # Careful, version 1.35 needed!!

if (! dir.exists("results")) dir.create("results")

# -------------------------------------------------------------------
# Pre-process data
# -------------------------------------------------------------------
## Load data / Remove rows with NAs
tmp <- readRDS("data/circforest_prepared_data_ibk_lag6.rds")
tmp <- tmp[, -grep("obertauern|montana|mariazell|guetsch|altdorf", names(tmp))]

tmp <- na.omit(tmp)

## Subset data to full hours
d <- tmp[as.POSIXlt(index(tmp))$min == 0L, ]; rm(tmp); gc()

## Remove very low values of ff and zero wind direction
d <- d[!d[, "dd.response"] == 0,]
d <- d[!d[, "ff.response"] < 1,]

## Transform response.dd from 0-360 degree to [-pi, pi]
d$dd.response <- d$dd.response / 360 * 2*pi
d$dd.response[d$dd.response > pi] <- d$dd.response[d$dd.response > pi] - 2*pi

## -------------------------------------------------------------------
## Fit climatology
## -------------------------------------------------------------------
#d_range <- as.POSIXlt(range(index(d)))
#
#pred_clim <- data.frame(mu = rep(NA, nrow(d)), kappa = rep(NA, nrow(d)))
#for(i in seq(1:nrow(d))){
#  cat(sprintf("Fitting climatology %2d%%\n", (i/nrow(d) * 100)%/%5 * 5))
#
#  i_plt <- as.POSIXlt(index(d)[i], origin = "1970-01-01")
#  
#  if((i_plt$mon + 1 == 2) & (i_plt$mday == 29)){
#    pred_clim[i, ] <- c("mu" = NA, "kappa" = NA)
#    next
#  }
#
#  i_years <- seq.int(min(d_range$year), max(d_range$year))[!(seq.int(min(d_range$year), 
#    max(d_range$year)) %in% i_plt$year)] + 1900
#
#  idx <- lapply(i_years, function(x) 
#    as.POSIXct(sprintf("%04d-%02d-%02d %02d:%02d:00", x, i_plt$mon + 1, i_plt$mday, 
#    i_plt$hour, i_plt$min), origin = "1970-01-01") + 60 * 60 * 24 * seq(-3, 3))
#
#  idx <- do.call(c, idx)
#  
#  train <- subset(d, index(d) %in% idx)
#
#  if(nrow(train) > 2){
#    climfit <- distexfit(train$dd.response,
#                         family = dist_vonmises())
#    pred_clim[i, ] <- coef(climfit, type = "parameter")
#  } else {
#    pred_clim[i, ] <- c("mu" = NA, "kappa" = NA)
#  }
#}
#
#save(pred_clim, file = "results/circforest_pred_clim_ibk_lag6_v3.rda")

# -------------------------------------------------------------------
# Fit persistence model (no CV)
# -------------------------------------------------------------------

pred_pers <- data.frame(mu = rep(NA, nrow(d)), kappa = rep(NA, nrow(d)))
for(i in seq(1:nrow(d))){
  cat(sprintf("Fitting persistence %2d%%\n", (i/nrow(d) * 100)%/%5 * 5))

  i_plt <- as.POSIXlt(index(d)[i], origin = "1970-01-01")
  
  if((i_plt$mon + 1 == 2) & (i_plt$mday == 29)){
    pred_pers[i, ] <- c("mu" = NA, "kappa" = NA)
    next
  }

  idx <- as.POSIXct(sprintf("%04d-%02d-%02d %02d:%02d:00", i_plt$year + 1900, i_plt$mon + 1, i_plt$mday, 
    i_plt$hour, i_plt$min), origin = "1970-01-01") + 60 * 60 * seq(-3, -1)

  train <- subset(d, index(d) %in% idx)

  if(nrow(train) > 2){
    climfit <- distexfit(train$dd.response,
                         family = dist_vonmises())
    pred_pers[i, ] <- coef(climfit, type = "parameter")
  } else {
    pred_pers[i, ] <- c("mu" = NA, "kappa" = NA)
  }
}

save(pred_pers, file = "results/circforest_pred_pers_ibk_lag6_v3.rda")

# -------------------------------------------------------------------
# Fit climatologyi (with CV)
# -------------------------------------------------------------------

## Fit models with cross-validation
cvID <- sort(rep(1:5, ceiling(nrow(d) / 5)))[1:nrow(d)]

pred_clim <- rep(list(data.frame(mu = rep(NA, nrow(d)), kappa = rep(NA, nrow(d)))), 5)

for(cv in unique(cvID)) { 
  cat(sprintf("Fitting models cv %s/%s\n", cv, max(cvID)))

  train <- d[cv != cvID, ]
  test <- d[cv == cvID, ]
 
  d_range <- as.POSIXlt(range(index(train)))

  for(i in seq(1:nrow(test))){
    cat(sprintf("Fitting climatology %2d%%\n", (i/nrow(test) * 100)%/%5 * 5))
  
    i_plt <- as.POSIXlt(index(test)[i], origin = "1970-01-01")
    
    if((i_plt$mon + 1 == 2) & (i_plt$mday == 29)){
      pred_clim[i, ] <- c("mu" = NA, "kappa" = NA)
      next
    }
  
    i_years <- seq.int(min(d_range$year), max(d_range$year))[!(seq.int(min(d_range$year), 
    max(d_range$year)) %in% i_plt$year)] + 1900

  
    idx <- lapply(i_years, function(x) 
      as.POSIXct(sprintf("%04d-%02d-%02d %02d:%02d:00", x, i_plt$mon + 1, i_plt$mday, 
      i_plt$hour, i_plt$min), origin = "1970-01-01") + 60 * 60 * 24 * seq(-3, 3))
  
    idx <- do.call(c, idx)
    
    train_subset <- subset(train, index(train) %in% idx)
  
    if(nrow(train) > 2){
      climfit <- distexfit(train_subset$dd.response,
                           family = dist_vonmises())
      pred_clim[[cv]][i, ] <- coef(climfit, type = "parameter")
    } else {
      pred_clim[[cv]][i, ] <- c("mu" = NA, "kappa" = NA)
    }
  }
}
  
save(pred_clim, file = "results/circforest_pred_clim_ibk_lag6_v3.rda")

# -------------------------------------------------------------------
# Fit Tree and Forest
# -------------------------------------------------------------------

## Transform zoo object to data.frame
tmp_time <- as.POSIXlt(index(d))
d <- as.data.frame(d)

## Add day of the year and time
d$daytime <- cut(tmp_time$hour, seq(0,24,by=3), include.lowest = TRUE)
d$doy <- tmp_time$yday + 1

## Set up formula
f <- as.formula(paste("dd.response ~ ", paste(names(d)[-grep("response", names(d))], collapse= "+")))

## Fit models with cross-validation
cvID <- sort(rep(1:5, ceiling(nrow(d) / 5)))[1:nrow(d)]

for(cv in unique(cvID)) { 
  cat(sprintf("Fitting models cv %s/%s\n", cv, max(cvID)))

  train <- d[cv != cvID, ]
  test <- d[cv == cvID, ]
  
  ## Fit tree
  m_dt <- distextree(formula = f,
                           data = train,
                           family = dist_vonmises(),
                           control = distextree_control(maxdepth = 4))
  
  ## Fit forest
  m_df <- distexforest(formula = f,
                       data = train,
                       family = dist_vonmises(),
                       ntree = 100,
                       mtry = 0.9 * ceiling(length(all.vars(f[[3]]))),
                       perturb = list(replace = FALSE, fraction = 0.3),
                       control = distextree_control(nmax = c("yx" = Inf, "z" = 50)))

  ## Predict models
  pred_dt.tmp <- predict(m_dt, newdata = test, type = "parameter")
  pred_df.tmp <- predict(m_df, newdata = test, type = "parameter")


  save(pred_dt.tmp, pred_df.tmp, 
    file = sprintf("results/circforest_pred_forest_ibk_lag6_v3_cv%s.rda", cv))
}

pred_dt <- pred_df <- list()
for(cv in unique(cvID)) {
  load(file = sprintf("results/circforest_pred_forest_ibk_lag6_v3_cv%s.rda", cv))
  pred_dt[[cv]] <- pred_dt.tmp
  pred_df[[cv]] <- pred_df.tmp
}


load(file = "results/circforest_pred_clim_ibk_lag6_v3.rda")

pred <- list(tree = do.call("rbind", pred_dt), forest = do.call("rbind", pred_df), 
  climatology = do.call("rbind", pred_clim), persistence = pred_pers)
obs  <- d$dd.response

## Remove nans (if any)
idx <- do.call(rbind, sapply(pred, function(x) which(is.na(x), arr.ind=TRUE)))
idx <- unique(idx[,1]) 
if(length(idx) > 0){
  pred_naomit <- lapply(pred, function(x) x[-idx, ])
  obs_naomit  <- d$dd.response[-idx]
} else{
  pred_naomit <- pred
  obs_naomit  <- d$dd.response
}

## Save results
#save(m_dt, m_df, file = "results/circforest_models_ibk_lag6.rds")
save(pred, pred_naomit, obs, obs_naomit, file = "results/circforest_results_ibk_lag6_v3.rds")

# -------------------------------------------------------------------
# Validate 
# -------------------------------------------------------------------

load(file = "results/circforest_results_ibk_lag6_v3.rds")

## Validate models
crps <- lapply(pred_naomit, function(x) crps_vonmises(mu = x$mu, kappa = x$kappa, 
  y = obs_naomit, sum = FALSE)) 

crps_grimit <- lapply(pred_naomit, function(x) sapply(1:nrow(x), function(i) 
  as.numeric(crps.circ(x = obs_naomit[i], mu = x[i, "mu"], kappa = x[i, "kappa"]))))

## Save validation
save(crps, crps_grimit, file = "results/circforest_validation_ibk_lag6_v3.rds")

## Plot tree
#circmax:::plot.circtree(dt[[1]], ep_args = list(justmin = 10), tp_args = list(type = "geographics"), 
#  ip_args = list(pval = FALSE))
