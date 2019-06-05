# -------------------------------------------------------------------
# - NAME:   circforest_fit_models_ibk.R
# - AUTHOR: Moritz N. Lang, Lisa Schlosser
# - DATE:   2019-06-03
# -------------------------------------------------------------------
# - PURPOSE:
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2019-06-05 on thinkmoritz
# -------------------------------------------------------------------

## Set time to utc (just in case)
Sys.setenv('TZ'='UTC')

## Load packages
library("zoo")
library("disttree")
library("circmax")

# -------------------------------------------------------------------
# Pre-process data
# -------------------------------------------------------------------
## Load data / Remove rows with NAs
tmp <- readRDS("data/circforest_prepared_data_ibk_lag3.rds")
tmp <- na.omit(tmp)

## Subset data to full hours
d <- tmp[as.POSIXlt(index(tmp))$min == 0L, ]; rm(tmp); gc()

## Transform response.dd from 0-360 degree to [-pi, pi]
d$dd.response <- d$dd.response / 360 * 2*pi
d$dd.response[d$dd.response > pi] <- d$dd.response[d$dd.response > pi] - 2*pi

# -------------------------------------------------------------------
# Fit climatology
# -------------------------------------------------------------------
d_range <- as.POSIXlt(range(index(d)))

pred_clim <- data.frame(mu = rep(NA, nrow(d)), kappa = rep(NA, nrow(d)))
for(i in seq(1:nrow(d))){
  cat(sprintf("Fitting climatology %2d%%\n", (i/nrow(d) * 100)%/%5 * 5))

  i_plt <- as.POSIXlt(index(d)[i], origin = "1970-01-01")
  i_years <- seq.int(min(d_range$year), max(d_range$year))[!(seq.int(min(d_range$year), 
    max(d_range$year)) %in% i_plt$year)] + 1900

  idx <- lapply(i_years, function(x) 
    as.POSIXct(sprintf("%04d-%02d-%02d %02d:%02d:00", x, i_plt$mon + 1, i_plt$mday, 
    i_plt$hour, i_plt$min), origin = "1970-01-01") + 60 * 60 * seq(-3, 3))

  idx <- do.call(c, idx)
  
  train <- subset(d, index(d) %in% idx)

  if(nrow(train) > 2){
    climfit <- distexfit(train$dd.response,
                         family = dist_vonmises())
    pred_clim[i, ] <- coef(climfit, type = "parameter")
  } else {
    pred_clim[i, ] <- c("mu" = NA, "kappa" = NA)
  }
}

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

m_dt <- m_df <- pred_dt <- pred_df <- list()
for(cv in unique(cvID)) { 
  cat(sprintf("Fitting models cv %s/%s\n", cv, max(cvID)))

  train <- d[cv != cvID, ]
  test <- d[cv == cvID, ]
  
  ## Fit tree
  m_dt[[cv]] <- distextree(formula = f,
                   data = train,
                   family = dist_vonmises(),
                   control = distextree_control(maxdepth = 4))
  
  ## Fit forest
  m_df[[cv]] <- distexforest(formula = f,
                     data = train,
                     family = dist_vonmises(),
                     ntree = 100)

  ## Predict models
  pred_dt[[cv]] <- predict(m_dt[[cv]], newdata = test, type = "parameter")
  pred_df[[cv]] <- predict(m_df[[cv]], newdata = test, type = "parameter")
}

# -------------------------------------------------------------------
# Validate and save models
# -------------------------------------------------------------------
pred <- list(tree = do.call("rbind", pred_dt), forest = do.call("rbind", pred_df), 
  climatology = pred_clim)

## Validate models
idx <- do.call(rbind, sapply(pred, function(x) which(is.na(x), arr.ind=TRUE)))
idx <- unique(idx[,1]) 
pred_naomit <- lapply(pred, function(x) x[-idx, ])

crps <- lapply(pred_naomit, function(x) crps_vonmises(mu = x$mu, kappa = x$kappa, 
  y = d$dd.response[-idx], sum = TRUE)) 

## Save output
if (! dir.exists("results")) dir.create("results")
#save(m_dt, m_df, file = "results/circforest_models_ibk_lag3.rds")
save(pred, pred_naomit, crps, file = "results/circforest_validation_ibk_lag3.rds")
