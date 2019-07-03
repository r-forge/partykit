# -------------------------------------------------------------------
# - NAME:   circforest_fit_models_ibk.R
# - AUTHOR: Moritz N. Lang, Lisa Schlosser
# - DATE:   2019-06-03
# -------------------------------------------------------------------
# - PURPOSE:
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2019-06-06 on thinkmoritz
# -------------------------------------------------------------------

## Model version for output files
my_vers <- "v4"

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
#
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
#
## Calculate exponential weights
calc_expweights <- function(theta, n){rev(sapply(1:n, function(x) theta^(x -1) * (1 - theta)))}
exp_weights <- calc_expweights(0.4, 6)

pred_pers <- data.frame(mu = rep(NA, nrow(d)), kappa = rep(NA, nrow(d)))
for(i in seq(1:nrow(d))){
  cat(sprintf("Fitting persistence %2d%%\n", (i/nrow(d) * 100)%/%5 * 5))

  i_plt <- as.POSIXlt(index(d)[i], origin = "1970-01-01")
  
  ## Skip 29th of February 
  if((i_plt$mon + 1 == 2) & (i_plt$mday == 29)){
    pred_pers[i, ] <- c("mu" = NA, "kappa" = NA)
    next
  }

  ## Subset training data set
  idx <- as.POSIXct(sprintf("%04d-%02d-%02d %02d:%02d:00", i_plt$year + 1900, i_plt$mon + 1, i_plt$mday, 
    i_plt$hour, i_plt$min), origin = "1970-01-01") + 60 * 60 * seq(-6, -1)

  train <- subset(d, index(d) %in% idx)
  train_weights <- exp_weights[idx %in% index(train)]

  ## Fit persistency
  persfit <- try(distexfit(as.numeric(train$dd.response),
                           family = dist_vonmises(), weights = train_weights))

  ## Predict parameters
  if(class(persfit) == "try-error") {
    pred_pers[i, ] <- c("mu" = NA, "kappa" = NA)
  } else {
    pred_pers[i, ] <- coef(persfit, type = "parameter")
  }
}

## Save predictions
saveRDS(pred_pers, file = sprintf("results/circforest_pred_pers_ibk_lag6_%s.rds", my_vers))
rm(pred_pers, persfit, train); gc()

# -------------------------------------------------------------------
# Fit climatology (with CV)
# -------------------------------------------------------------------

## Fit models with cross-validation
cvID <- sort(rep(1:5, ceiling(nrow(d) / 5)))[1:nrow(d)]

pred_clim <- lapply(unique(cvID), function(x) data.frame(mu = rep(NA, sum(cvID == x)), 
  kappa = rep(NA, sum(cvID == x))))

for(cv in unique(cvID)) { 
  cat(sprintf("Fitting models cv %s/%s\n", cv, max(cvID)))

  train <- d[cv != cvID, ]
  test <- d[cv == cvID, ]
 
  d_range <- as.POSIXlt(range(index(train)))

  for(i in seq(1:nrow(test))){
    cat(sprintf("Fitting climatology %2d%%\n", (i/nrow(test) * 100)%/%5 * 5))
  
    i_plt <- as.POSIXlt(index(test)[i], origin = "1970-01-01")
    
    ## Skip 29th of February 
    if((i_plt$mon + 1 == 2) & (i_plt$mday == 29)){
      pred_clim[[cv]][i, ] <- c("mu" = NA, "kappa" = NA)
      next
    }
  
    ## Subset training data set
    i_years <- seq.int(min(d_range$year), max(d_range$year))[!(seq.int(min(d_range$year), 
    max(d_range$year)) %in% i_plt$year)] + 1900

    idx <- lapply(i_years, function(x) 
      as.POSIXct(sprintf("%04d-%02d-%02d %02d:%02d:00", x, i_plt$mon + 1, i_plt$mday, 
      i_plt$hour, i_plt$min), origin = "1970-01-01") + 60 * 60 * 24 * seq(-3, 3))
  
    idx <- do.call(c, idx)
    
    train_subset <- subset(train, index(train) %in% idx)
  
    ## Fit climatology
    climfit <- try(distexfit(train_subset$dd.response,
                             family = dist_vonmises()))

    ## Predict parameters
    if(class(climfit) == "try-error") {
      pred_clim[[cv]][i, ] <- c("mu" = NA, "kappa" = NA)
    } else {
      pred_clim[[cv]][i, ] <- coef(climfit, type = "parameter")
    }
  }
}
  
## Save predictions
saveRDS(pred_clim, file = sprintf("results/circforest_pred_clim_ibk_lag6_%s.rds", my_vers))
rm(pred_clim, climfit, train, test, train_subset); gc()

# -------------------------------------------------------------------
# Fit tree and forest
# -------------------------------------------------------------------

## Transform zoo object to data.frame
tmp_time <- as.POSIXlt(index(d))
d <- as.data.frame(d)

## Add day of the year and time
d$daytime <- cut(tmp_time$hour, seq(0,24,by=3), include.lowest = TRUE)
d$doy <- tmp_time$yday + 1

rm(tmp_time); gc()

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
                       mtry = ceiling(1. * length(all.vars(f[[3]]))),
                       perturb = list(replace = FALSE, fraction = 0.3),
                       control = distextree_control(nmax = c("yx" = Inf, "z" = 50)))

  ## Predict models
  pred_dt.tmp <- predict(m_dt, newdata = test, type = "parameter")
  pred_df.tmp <- predict(m_df, newdata = test, type = "parameter")

  ## Save predictions (and models)
  saveRDS(m_dt, 
    file = sprintf("results/circforest_model_tree_ibk_lag6_%s_cv%s.rds", my_vers, cv))
  saveRDS(pred_dt.tmp, 
    file = sprintf("results/circforest_pred_tree_ibk_lag6_%s_cv%s.rds", my_vers, cv))
  saveRDS(pred_df.tmp, 
    file = sprintf("results/circforest_pred_forest_ibk_lag6_%s_cv%s.rds", my_vers, cv))
  rm(pred_dt.tmp, pred_df.tmp, m_df, m_dt, train, test); gc()
}

# -------------------------------------------------------------------
# Combine predictions
# -------------------------------------------------------------------

## Load predictions
pred_dt <- pred_df <- list()
for(cv in unique(cvID)) {
  pred_dt.tmp <- readRDS(file = sprintf("results/circforest_pred_tree_ibk_lag6_%s_cv%s.rds", my_vers, cv))
  pred_df.tmp <- readRDS(file = sprintf("results/circforest_pred_forest_ibk_lag6_%s_cv%s.rds", my_vers, cv))
  pred_dt[[cv]] <- pred_dt.tmp
  pred_df[[cv]] <- pred_df.tmp
  rm(pred_df.tmp, pred_dt.tmp); gc()
}

pred_pers <- readRDS(file = sprintf("results/circforest_pred_pers_ibk_lag6_%s.rds", my_vers))
pred_clim <- readRDS(file = sprintf("results/circforest_pred_clim_ibk_lag6_%s.rds", my_vers))

## Combine predictions
pred <- list(tree = do.call("rbind", pred_dt), forest = do.call("rbind", pred_df), 
  climatology = do.call("rbind", pred_clim), persistence = pred_pers)
obs  <- d$dd.response

## Remove nans (if any)
idx <- do.call(rbind, sapply(pred, function(x) which(is.na(x), arr.ind=TRUE)))
idx <- unique(idx[,1]) 
if(length(idx) > 0){
  pred_naomit <- lapply(pred, function(x) x[-idx, ])
  obs_naomit  <- obs[-idx]
} else{
  pred_naomit <- pred
  obs_naomit  <- obs
}

## Save results
save(pred, pred_naomit, obs, obs_naomit, 
  file = sprintf("results/circforest_results_ibk_lag6_%s.rda", my_vers))

# -------------------------------------------------------------------
# Validate predictions
# -------------------------------------------------------------------

load(file = sprintf("results/circforest_results_ibk_lag6_%s.rda", my_vers))

### Validate models
#crps <- lapply(pred_naomit, function(x) crps_vonmises(mu = x$mu, kappa = x$kappa, 
#  y = obs_naomit, sum = FALSE)) 

crps_grimit <- lapply(pred_naomit, function(x) sapply(1:nrow(x), function(i) 
  as.numeric(crps.circ(x = obs_naomit[i], mu = x[i, "mu"], kappa = x[i, "kappa"]))))

## Save validation
save(crps_grimit, file = sprintf("results/circforest_validation_ibk_lag6_%s.rda", my_vers))

## -------------------------------------------------------------------
## Plotting
## -------------------------------------------------------------------
if(FALSE){
  my_vers <- "v4"

  ## Plot CRPS
  load(file = sprintf("results/circforest_validation_ibk_lag6_%s.rda", my_vers))

  X11(width = 8, height = 4.5)
  par(mar = c(3.1, 4.1, 2.1, 2.1))
  boxplot(sapply(crps_grimit[c("climatology", "persistence", "tree", "forest")], 
    function(x) (1 - x/crps_grimit[["climatology"]]) * 100), ylim = c(-60, 110), col = gray(0.6), 
    ylab = "CRPS skill ccore [%]")
  dev.print(pdf, sprintf("results/_plot_circforest_validation_crpsskill_ibk_lag6_%s.pdf", my_vers))
  
  X11(width = 8, height = 4.5)
  par(mar = c(3.1, 4.1, 2.1, 2.1))
  boxplot(crps_grimit[c("climatology", "persistence", "tree", "forest")], ylim = c(0, 1), col = gray(0.6), 
    ylab = "CRPS [m/s]")
  dev.print(pdf, sprintf("results/_plot_circforest_validation_crpsraw_ibk_lag6_%s.pdf", my_vers))

  ## Plot single tree
  cv <- 1
  m_dt <- readRDS(file = sprintf("results/circforest_model_tree_ibk_lag6_%s_cv%s.rds", my_vers, cv))
  circmax:::plot.circtree(m_dt, ep_args = list(justmin = 10), tp_args = list(type = "geographics"), 
    ip_args = list(pval = FALSE))
}
