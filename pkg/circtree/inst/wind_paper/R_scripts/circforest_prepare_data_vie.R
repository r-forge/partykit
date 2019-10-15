# -------------------------------------------------------------------
# - NAME:   circforest_prepare_data_vie.R
# - AUTHOR: Moritz N. Lang
# - DATE:   2019-05-27
# -------------------------------------------------------------------
# - PURPOSE: Create data file for VIE with temporal/spatial differences
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2019-10-15 on thinkmoritz
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Preliminary
# -------------------------------------------------------------------
## Load packages
library(zoo)
library(foehnix)  # Reto's package
library(STAGE)  # Reto's package

data_in <- "/home/moritz/Projects/profcast_project/scripts/Rdevelopment/data_out"
data_out <- "data"

mylag <- 6


# -------------------------------------------------------------------
# Small helper functions
# -------------------------------------------------------------------
## Helper function to make a "strict" time series
make_strict <- function(x, dates = NULL) {
  if(is.null(dates)){
    delta <- as.numeric(min(diff(index(x))), unit = "secs")
    all <- seq(as.numeric(min(index(x))), as.numeric(max(index(x))), by = delta)
    x <- merge(x, zoo(, as.POSIXct(all, origin = "1970-01-01")), all = c(FALSE, TRUE))
  } else {
    x <- merge(x, zoo(, as.POSIXct(dates, origin = "1970-01-01")), all = c(FALSE, TRUE))
  }
  return(x)
}

calc_anglediff <- function(ff1, dd1, ff2, dd2, name) {
  z <- merge(ff1, dd1, ff2, dd2)
  #z <- na.omit(z)
  out <- data.frame(ff = as.numeric(z$ff1 - z$ff2))
  a <- as.numeric(z$dd1 - z$dd2)
  a <- matrix(a, ncol = 3, nrow = length(a), byrow = FALSE)
  b <- matrix(c(-360, 0, 360), byrow = TRUE, ncol = 3, nrow = nrow(a))
  c <- (a - b)
  out$dd <- c[cbind(1:nrow(c), apply(c, 1, function(x) {ret <- which.min(abs(x)); if(length(ret) == 0) ret <- NA; return(ret)}))]
  #out$dd <- sapply(1:length(dd1), function(x) {tmp <- as.numeric(dd1[x] - dd2[x]) + c(-360, 0, 360); tmp[which.min(abs(tmp))]})

  out <- zoo(out, index(z))
  if(!is.null(name)){
    names(out) <- paste0("diff", name, "_", names(out))
  } else {
    names(out) <- paste0("diff", "_", names(out))
  }
  return(out)
}


# -------------------------------------------------------------------
# Load observations from databases using LOWVIS mysql database for the observations (STAGE::lowvis)
# -------------------------------------------------------------------

# Getting lowvisDescription file
obsfile <- sprintf("data/_pure_observations_lowvisdb_update20191015.rds")
if (! file.exists(obsfile)) {

  ## Loaded distinct sensorids from database directly
  sensors <- c(2, 12, 22, 32, 161, 200, 201, 205, 206, 208, 420, 502, 507, 522,
               527, 542, 547, 562, 567, 582, 587, 590, 594, 602, 607, 610, 614, 622, 627,
               630, 634, 662, 666, 667, 670, 674, 682, 687, 690, 694, 702, 707, 710, 714)

  ## Loading year by year
  obs_pure <- list()
  for (yr in 2013:2018) {
    cat(sprintf("   Fetching observations for year %04d\n", yr))
    bgn <- as.POSIXct(sprintf("%04d-01-01 00:00:00", yr))
    end <- as.POSIXct(sprintf("%04d-12-31 12:59:59", yr))
    tmp <- lowvis("obs", sensorids = sensors, begin = bgn, end = end)
    obs_pure[[sprintf("year_%04d",yr)]] <- tmp; rm(tmp)
  }
  hold <- obs_pure
  obs_pure <- do.call("rbind", obs_pure)
  attr(obs_pure, "sensors") <- sensors
  attr(obs_pure, "created") <- Sys.time()
  attr(obs_pure, "createdon") <- Sys.info()["nodename"]

  ## Save both objects
  saveRDS(file = obsfile, obs_pure)
} else {
  cat(sprintf("   File \"%s\" exists, load and continue\n", obsfile))
  obs_pure <- readRDS(obsfile)
}

# -------------------------------------------------------------------
# Prepare lowvis data
# -------------------------------------------------------------------
## Create lagged observations (obs_pure)   
obs_pure <- make_strict(obs_pure)
if(unique(diff(index(data))) != 10) stop("wrong temporal resolution, suspected to 10min intervals")
data <- lag(obs_pure, -mylag)
rm(obs_pure); gc()

## Calculate lagged measurements
## 1 hour change
for(i_name in names(data)[-grep("vv|sws|lvp|cei|rvr|vis|dd|ff|ch30min", names(data))]){
  cmd <- sprintf("data$%1$s_ch1h <- data$%1$s - lag(data$%1$s, -6)",
    i_name)
  eval(parse(text = cmd))
}

tmp_names <- names(data)[grep("^dd", names(data))]
tmp_names <- gsub('dd', '', tmp_names)
for(idx in tmp_names){
  diff <- calc_anglediff(data[, paste0("ff", idx)], data[, paste0("dd", idx)], 
    lag(data[, paste0("ff", idx)], -6), lag(data[, paste0("dd", idx)], -6), paste0(idx, "_ch1h"))
  data <- merge(data, diff, all = c(TRUE, FALSE))
}

## 3 hour change
for(i_name in names(data)[-grep("vv|sws|lvp|cei|rvr|vis|dd|ff|ch30min|ch1h", names(data))]){
  cmd <- sprintf("data$%1$s_ch3h <- data$%1$s - lag(data$%1$s, -18)",
    i_name)
  eval(parse(text = cmd))
}

tmp_names <- names(data)[grep("^dd", names(data))]
tmp_names <- gsub('dd', '', tmp_names)
for(idx in tmp_names){
  diff <- calc_anglediff(data[, paste0("ff", idx)], data[, paste0("dd", idx)], 
    lag(data[, paste0("ff", idx)], -18), lag(data[, paste0("dd", idx)], -18), paste0(idx, "_ch3h"))
  data <- merge(data, diff, all = c(TRUE, FALSE))
}

# -------------------------------------------------------------------
# Merge 'data' with 'response' and calculate angle differences
# -------------------------------------------------------------------
#obs_pure_nolag <- readRDS(obsfile)
#dd.response <- obs_pure_nolag$dd34
#ff.response <- obs_pure_nolag$ff34
#rm(obs_pure_nolag); gc()
tawes_around_nolag <- readRDS("data/tawes_around_vie.rds")
dd.response <- tawes_around_nolag$dd.station11036
ff.response <- tawes_around_nolag$ff.station11036
rm(tawes_around_nolag); gc()

data <- merge(dd.response, data, all = c(TRUE, FALSE))
data <- merge(ff.response, data, all = c(TRUE, FALSE))

diffEXB <- calc_anglediff(data$ffEXB, data$ddEXB, ff.response, dd.response, "EXB")
diffTOW <- calc_anglediff(data$ffTOW, data$ddTOW, ff.response, dd.response, "TOW")
diffOM29 <- calc_anglediff(data$ffOM29, data$ddOM29, ff.response, dd.response, "OM29")
diffARS <- calc_anglediff(data$ffARS, data$ddARS, ff.response, dd.response, "ARS")

data <- merge(data, diffEXB, diffTOW, diffOM29, diffARS)


# -------------------------------------------------------------------
# Prepare tawes data
# -------------------------------------------------------------------
tawes_around <- readRDS("data/tawes_around_vie.rds")
tawes_around <- lag(tawes_around, -mylag)

## Calculate lagged tawes measurements
## 1 hour change
for(i_name in names(tawes_around)[-grep("dd|ff|ch30min", names(tawes_around))]){
  cmd <- sprintf("tawes_around$%1$s_ch1h <- tawes_around$%1$s - lag(tawes_around$%1$s, -6)",
    i_name)
  eval(parse(text = cmd))
}

tmp_names <- names(tawes_around)[grep("^dd", names(tawes_around))]
tmp_names <- gsub('dd', '', tmp_names)
for(idx in tmp_names){
  diff <- calc_anglediff(tawes_around[, paste0("ff", idx)], tawes_around[, paste0("dd", idx)],
    lag(tawes_around[, paste0("ff", idx)], -6), lag(tawes_around[, paste0("dd", idx)], -6), paste0(idx, "_ch1h"))
  tawes_around <- merge(tawes_around, diff, all = c(TRUE, FALSE))
}

## 3 hour change
for(i_name in names(tawes_around)[-grep("dd|ff|ch30min|ch1h", names(tawes_around))]){
  cmd <- sprintf("tawes_around$%1$s_ch3h <- tawes_around$%1$s - lag(tawes_around$%1$s, -18)",
    i_name)
  eval(parse(text = cmd))
}

tmp_names <- names(tawes_around)[grep("^dd", names(tawes_around))]
tmp_names <- gsub('dd', '', tmp_names)
for(idx in tmp_names){
  diff <- calc_anglediff(tawes_around[, paste0("ff", idx)], tawes_around[, paste0("dd", idx)],
    lag(tawes_around[, paste0("ff", idx)], -18), lag(tawes_around[, paste0("dd", idx)], -18), paste0(idx, "_ch3h"))
  tawes_around <- merge(tawes_around, diff, all = c(TRUE, FALSE))
}


# -------------------------------------------------------------------
# Merge 'data'  with 'tawes' and calculate angle differences
# -------------------------------------------------------------------
data <- merge(data, tawes_around, all = c(TRUE, FALSE))

for(idx in unique(regmatches(names(tawes_around), regexpr("station[0-9]{5}", names(tawes_around))))){
  diff <- calc_anglediff(data[, paste0("ff.", idx)], data[, paste0("dd.", idx)], ff.response, dd.response, idx)
  data <- merge(data, diff, all = c(TRUE, FALSE))
}

## Save data
if (! dir.exists(data_out)) dir.create(data_out)
saveRDS(data, paste0(data_out, "/circforest_prepared_data_vie_lag", mylag, "_update20191015.rds"))


#data <- as.data.frame(data)
#
## Remove all columns with too little data
#tmp <- which(colSums(is.na(data)) / nrow(data) > 0.1)
#if (length(tmp) > 0) data <- data[, -tmp]
#
## Remove NA's
#data <- na.omit(data)
#cat(sprintf("   We ended up with a matrix of size: %d x %d\n", nrow(data), ncol(data)))

