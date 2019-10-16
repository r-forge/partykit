# -------------------------------------------------------------------
# - NAME:   circforest_prepare_data_vie.R
# - AUTHOR: Moritz N. Lang
# - DATE:   2019-10-16
# -------------------------------------------------------------------
# - PURPOSE: Create data file for IBK with temporal/spatial differences
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2019-10-16 on thinkmoritz
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Preliminary
# -------------------------------------------------------------------
## Load packages
library(zoo)
library(foehnix)  # Reto's package
library(STAGE)  # Reto's package

mylag <- 3
response_station <- "11036"

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

calc_anglediff <- function(ff1, dd1, ff2, dd2, name = NULL) {
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

make_nicename <- function(name){
  name <- as.character(name)
  name <- gsub('^ | $', '', name)
  name <- gsub(' ', '_', name)
  name <- gsub('-', '_', name)
  name <- gsub('\\/', '_', name)
  name <- gsub('(\\(|\\))', "", name)
  name <- paste0(name, collapse = "")
  name <- tolower(name)
  name
}

## Create data directory if necessary
if (!dir.exists("data")) dir.create("data")

# -------------------------------------------------------------------
# Read observations from databases using LOWVIS mysql database for the observations (STAGE::lowvis)
# -------------------------------------------------------------------
lowvis_file <- sprintf("data/_observations_lowvisdb_update20191015.rds")
if (! file.exists(lowvis_file)) {

  ## Loaded distinct sensorids from database directly
  sensors <- c(2, 12, 22, 32, 161, 200, 201, 205, 206, 208, 420, 502, 507, 522,
               527, 542, 547, 562, 567, 582, 587, 590, 594, 602, 607, 610, 614, 622, 627,
               630, 634, 662, 666, 667, 670, 674, 682, 687, 690, 694, 702, 707, 710, 714)

  ## Loading year by year
  lowvis_data <- list()
  for (yr in 2013:2018) {
    cat(sprintf("   Fetching observations for year %04d\n", yr))
    bgn <- as.POSIXct(sprintf("%04d-01-01 00:00:00", yr))
    end <- as.POSIXct(sprintf("%04d-12-31 12:59:59", yr))
    tmp <- lowvis("obs", sensorids = sensors, begin = bgn, end = end)
    lowvis_data[[sprintf("year_%04d",yr)]] <- tmp; rm(tmp)
  }

  lowvis_data <- do.call("rbind", lowvis_data)
  attr(lowvis_data, "sensors") <- sensors
  attr(lowvis_data, "created") <- Sys.time()
  attr(lowvis_data, "createdon") <- Sys.info()["nodename"]

  ## Save both objects
  saveRDS(file = lowvis_file, lowvis_data)
  rm(lowvis_data); gc()
} 

# -------------------------------------------------------------------
# Read TAWES around airport from databeses at ZAMG (STAGE::sybase)
# -------------------------------------------------------------------

loww <- list(lat = 48.110833, lon = 16.570833)
stations <- sybasestations("*")

stations$d <- sqrt((stations$lon - loww$lon)^2 + (stations$lat-loww$lat)^2)
stations <- subset(stations,d < .4 & tawes == 1)

params <- c("dd", "ff", "ffx", "tl", "rf", "p", "pred")

sybase_stationlist <- list()
for (i in 1:nrow(stations)) {
  outfile <- sprintf("data/_observations_sybase_tawes_around_vie_%d_update20191015.rds", stations$statnr[i])
  sybase_stationlist[[i]] <- outfile
  if (file.exists(outfile)) next
  
  tmp <- sybase("tawes", stations$statnr[i], parameter = params,
    begin = "2013-01-01", end = "2018-12-31", archive = TRUE)
  
  saveRDS(file = outfile, tmp)
  rm(tmp); gc()
}
sybase_stationlist <- do.call("c", sybase_stationlist)

# -------------------------------------------------------------------
# Read data sets
# -------------------------------------------------------------------
cat("\nPrepare data ...\n")

datafile <- sprintf("data/circforest_prepared_data_vie_lag%sh_update20191015.rds", mylag)

## Load or create (and save) data set
if(file.exists(datafile)) {

  cat(sprintf("File \"%s\" exists, loading data ...\n", datafile))
  data <- readRDS(datafile)

} else {

  file_names <- sybase_stationlist

  ## Read response and make strict zoo object
  response <- eval(parse(text = sprintf("response <- make_strict(readRDS('%s'))", file_names[grep(response_station, file_names)])))
  if(unique(diff(index(response))) != 10) stop("wrong temporal resolution, suspected to 10min intervals")

  ## Subset to full hours
  response <- response[as.POSIXlt(index(response))$min == 0L, ]

  ## Calculate lagged response for spatial differences
  eval(parse(text = sprintf("response_lag <- lag(response, %s)", -mylag)))

  # -------------------------------------------------------------------
  # Prepare tawes data
  # -------------------------------------------------------------------
  stat_names <- NULL
  data <- NULL
  for(file in file_names){
    tmp <- sub('data/_observations_sybase_tawes_around_vie_', '', file)
    tmp <- sub('_update20191015.rds', '', tmp)
    tmp <- make_nicename(subset(stations, statnr == tmp)$name)

    cat(sprintf("\nPrepare data for station '%s'\n", tmp))
    
    eval(parse(text = sprintf("%s <- make_strict(readRDS(file))", tmp)))

    ## Check if 10min interval 
    eval(parse(text = sprintf("if(unique(diff(index(%s))) != 10) stop('wrong temporal resolution, suspected to 10min intervals')", tmp)))

    ## Subset to full hours
    eval(parse(text = sprintf("%1$s <- %1$s[as.POSIXlt(index(%1$s))$min == 0L, ]", tmp)))

    ## Make lags
    eval(parse(text = sprintf("%1$s <- lag(%1$s, %2$s)", tmp, -mylag)))
    eval(parse(text = sprintf("%1$s <- make_strict(%1$s, index(response_lag))", tmp)))

    ## Create differences to Vienna airport
    if (tmp != make_nicename(subset(stations, statnr == response_station)$name)) {
      for(i_param in names(eval(parse(text = tmp)))){
        if(grepl("dd", i_param) & "ff" %in% names(eval(parse(text = tmp)))){
          eval(parse(text = sprintf("%1$s$dd_diff <- calc_anglediff(response_lag$ff, response_lag$dd, 
            %1$s$ff, %1$s$dd)$diff_dd", tmp)))
        } else if (i_param %in% c("p", "psta")){
          eval(parse(text = sprintf("%1$s$%2$s_diff <- response_lag$p - %1$s$%2$s", tmp, i_param)))
        } else if (i_param %in% c("pred", "predsta")){
          eval(parse(text = sprintf("%1$s$%2$s_diff <- response_lag$pred - %1$s$%2$s", tmp, i_param)))
        } else if (i_param %in% c("t", "tl")){
          eval(parse(text = sprintf("%1$s$%2$s_diff <- response_lag$t - %1$s$%2$s", tmp, i_param)))
        } else if (i_param %in% c("rh", "rf")){
          eval(parse(text = sprintf("%1$s$%2$s_diff <- response_lag$rf - %1$s$%2$s", tmp, i_param)))
        } else if (i_param %in% c("ffx", "ffx1s")){
          eval(parse(text = sprintf("%1$s$%2$s_diff <- response_lag$ffx - %1$s$%2$s", tmp, i_param)))
        } else if (i_param %in% c("ff")){
          eval(parse(text = sprintf("%1$s$%2$s_diff <- response_lag$ff - %1$s$%2$s", tmp, i_param)))
        }
      }
    }
 
    ## Get ff maximum, minimum and mean over last three hours
    for(i_param in names(eval(parse(text = tmp)))){
      if (i_param %in% c("ff", "ffx", "ffx1s")){
        eval(parse(text = sprintf("%1$s$%2$s_max <- rollapply(%1$s$%2$s, width = 3, FUN = max, fill = NA, align = 'right')", 
          tmp, i_param)))
        eval(parse(text = sprintf("%1$s$%2$s_min <- rollapply(%1$s$%2$s, width = 3, FUN = min, fill = NA, align = 'right')", 
          tmp, i_param)))
        eval(parse(text = sprintf("%1$s$%2$s_mean <- rollapply(%1$s$%2$s, width = 3, FUN = mean, fill = NA, align = 'right')", 
          tmp, i_param)))
      }
    }

    ## Calculate temporal changes
    for(i_param in names(eval(parse(text = tmp)))[!grepl("(min|max|mean|diff)", names(eval(parse(text = tmp))))]){
      if(grepl("dd", i_param)){
        eval(parse(text = sprintf("%1$s$%2$s_ch1h <- calc_anglediff(%1$s$ff, %1$s$dd, 
          lag(%1$s$ff, -1), lag(%1$s$dd, -1))$diff_dd", tmp, i_param)))
        eval(parse(text = sprintf("%1$s$%2$s_ch3h <- calc_anglediff(%1$s$ff, %1$s$dd, 
          lag(%1$s$ff, -3), lag(%1$s$dd, -3))$diff_dd", tmp, i_param)))
      } else {
        eval(parse(text = sprintf("%1$s$%2$s_ch1h <- %1$s$%2$s - lag(%1$s$%2$s, -1)", tmp, i_param)))
        eval(parse(text = sprintf("%1$s$%2$s_ch3h <- %1$s$%2$s - lag(%1$s$%2$s, -3)", tmp, i_param)))
      }
    }

    if(is.null(data)) {
      names(response) <- paste0(names(response), ".response")
      eval(parse(text = sprintf("names(%1$s) <- paste0(names(%1$s), '.%1$s')", tmp)))
      eval(parse(text = sprintf("data <- merge(response, %s)", tmp)))
    } else {
      eval(parse(text = sprintf("names(%1$s) <- paste0(names(%1$s), '.%1$s')", tmp)))
      eval(parse(text = sprintf("data <- merge(data, %s)", tmp)))
    }  

    stat_names <- c(stat_names, tmp)
    eval(parse(text = sprintf("rm(%s)", tmp)))

  }


  # -------------------------------------------------------------------
  # Prepare lowvis data
  # -------------------------------------------------------------------
  lowvis_data <- make_strict(readRDS(lowvis_file))

  ## Check if 10min interval 
  if(unique(diff(index(lowvis_data))) != 10) stop("wrong temporal resolution, suspected to 10min intervals")

  ## Subset to full hours
  lowvis_data <- lowvis_data[as.POSIXlt(index(lowvis_data))$min == 0L, ]

  # Make lags
  lowvis_data <- lag(lowvis_data, -mylag)
  lowvis_data <- make_strict(lowvis_data, index(response_lag))

  ## Convert airport measurements in kt to m/s
  lowvis_data[, grep("^ff", names(lowvis_data))] <- lowvis_data[, grep("^ff", names(lowvis_data))] / 1.9438444924574
  
  ## Calculate lagged measurements
  ## 1 hour change
  for(i_name in names(lowvis_data)[-grep("vv|sws|lvp|cei|rvr|vis|dd|ff|ch30min", names(lowvis_data))]){
    cmd <- sprintf("lowvis_data$%1$s_ch1h <- lowvis_data$%1$s - lag(lowvis_data$%1$s, -1)",
      i_name)
    eval(parse(text = cmd))
  }
  
  tmp_names <- names(lowvis_data)[grep("^dd", names(lowvis_data))]
  tmp_names <- gsub('dd', '', tmp_names)
  for(idx in tmp_names){
    diff <- calc_anglediff(lowvis_data[, paste0("ff", idx)], lowvis_data[, paste0("dd", idx)],
      lag(lowvis_data[, paste0("ff", idx)], -1), lag(lowvis_data[, paste0("dd", idx)], -1), paste0(idx, "_ch1h"))
    lowvis_data <- merge(lowvis_data, diff, all = c(TRUE, FALSE))
  }
  
  ## 3 hour change
  for(i_name in names(lowvis_data)[-grep("vv|sws|lvp|cei|rvr|vis|dd|ff|ch30min|ch1h", names(lowvis_data))]){
    cmd <- sprintf("lowvis_data$%1$s_ch3h <- lowvis_data$%1$s - lag(lowvis_data$%1$s, -3)",
      i_name)
    eval(parse(text = cmd))
  }
  
  tmp_names <- names(lowvis_data)[grep("^dd", names(lowvis_data))]
  tmp_names <- gsub('dd', '', tmp_names)
  for(idx in tmp_names){
    diff <- calc_anglediff(lowvis_data[, paste0("ff", idx)], lowvis_data[, paste0("dd", idx)],
      lag(lowvis_data[, paste0("ff", idx)], -3), lag(lowvis_data[, paste0("dd", idx)], -3), paste0(idx, "_ch3h"))
    lowvis_data <- merge(lowvis_data, diff, all = c(TRUE, FALSE))
  }

  ## Calculate spatial differences 
  diffEXB <- calc_anglediff(response_lag$ff, response_lag$dd, lowvis_data$ffEXB, lowvis_data$ddEXB, "EXB")
  diffTOW <- calc_anglediff(response_lag$ff, response_lag$dd, lowvis_data$ffTOW, lowvis_data$ddTOW, "TOW")
  diffOM29 <- calc_anglediff(response_lag$ff, response_lag$dd, lowvis_data$ffOM29, lowvis_data$ddOM29, "OM29")
  diffARS <- calc_anglediff(response_lag$ff, response_lag$dd, lowvis_data$ffARS, lowvis_data$ddARS, "ARS")
  
  lowvis_data <- merge(lowvis_data, diffEXB, diffTOW, diffOM29, diffARS)

  data <- merge(data, lowvis_data)

  cat(sprintf("File \"%s\" does not exists, saving data ...\n", datafile))
  saveRDS(data, datafile)

}
