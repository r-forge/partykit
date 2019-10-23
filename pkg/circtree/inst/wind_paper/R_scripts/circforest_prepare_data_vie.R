# -------------------------------------------------------------------
# - NAME:   circforest_prepare_data_vie.R
# - AUTHOR: Moritz N. Lang
# - DATE:   2019-10-23
# -------------------------------------------------------------------
# - PURPOSE: Create data file for VIE with temporal/spatial differences
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2019-10-23 on thinkmoritz
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Preliminary
# -------------------------------------------------------------------
## Load packages
library(zoo)
library(foehnix)  # Reto's package
library(STAGE)  # Reto's package

version <- "update20191023"
mylag <- 1 ## in hour
start_date <- "2014-01-01"
end_date <- "2018-12-31"
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


## Helper function to find smallest not-absolute value
calc_angle_dist <- function(a, b, direction = FALSE, unit = c("rad", "deg")) {
  unit <- match.arg(unit)
  if(!is.logical(direction)) stop("arg 'absolute' needs to be logical")

  d <- abs(b - a)

  ## Find period smallest difference
  if (unit == "rad") {
    rval <- pmin(d, 2*pi - d)
  } else {
    rval <- pmin(d, 360 - d)
  }

  ## If direction should be maintained
  if (direction) {
    rval <- sign(b - a) * rval
  }

  return(rval)
}

## Helper function to check duplicates
rmduplicated <- function(x) {
  idx <- which(duplicated(index(x)))
  if (length(idx) > 0) x <- x[-idx, ]
  x
}

## Helper function to make nice station names
make_nice_stationname <- function(name){
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
lowvis_file <- sprintf("data/_observations_lowvisdb_vie_%s_%s_%s.rds", start_date, end_date, version)
if (!file.exists(lowvis_file)) {

  cat(sprintf("File \"%s\" does not exists, prepare data ...\n", lowvis_file))

  ## Loaded distinct sensorids from database directly
  sensors <- c(2, 12, 22, 32, 161, 200, 201, 205, 206, 208, 420, 502, 507, 522,
               527, 542, 547, 562, 567, 582, 587, 590, 594, 602, 607, 610, 614, 622, 627,
               630, 634, 662, 666, 667, 670, 674, 682, 687, 690, 694, 702, 707, 710, 714)

  ## Loading year by year
  lowvis_data <- list()
  for (yr in (as.POSIXlt(start_date)$year:as.POSIXlt(end_date)$year) + 1900) {
    cat(sprintf("   Fetching observations for year %04d\n", yr))
    bgn <- as.POSIXct(sprintf("%04d-01-01 00:00:00", yr))
    end <- as.POSIXct(sprintf("%04d-12-31 12:59:59", yr))
    tmp <- lowvis("obs", sensorids = sensors, begin = bgn, end = end)
    lowvis_data[[sprintf("year_%04d",yr)]] <- tmp; rm(tmp)
  }

  lowvis_data <- do.call("rbind", lowvis_data)

  ## Transform kt to m/s
  lowvis_data[, grep("^ff", names(lowvis_data))] <- lowvis_data[, grep("^ff", names(lowvis_data))] / 1.9438444924574

  attr(lowvis_data, "sensors") <- sensors
  attr(lowvis_data, "created") <- Sys.time()
  attr(lowvis_data, "createdon") <- Sys.info()["nodename"]

  ## Modify prefix 
  pattern <- "(11|16|29|34|TOW|FMA|KLI|ARS|EXB|OM34|666|OM11|OM16|OM29|OM34)"
  name_postfix <- gsub(pattern, "", names(lowvis_data), perl = TRUE)
  name_prefix <- regmatches(names(lowvis_data), regexpr(pattern, names(lowvis_data)))
  names(lowvis_data) <- paste0(name_prefix, ".", name_postfix)
  names(lowvis_data)[grep("^[0-9]", names(lowvis_data))] <- paste0("st", names(lowvis_data)[grep("^[0-9]", names(lowvis_data))])

  ## Save lowvis data
  saveRDS(lowvis_data, file = lowvis_file)
  rm(lowvis_data); gc()
}


# -------------------------------------------------------------------
# Read TAWES around airport from databeses at ZAMG (STAGE::sybase)
# -------------------------------------------------------------------

## Find stations around airport
loww <- list(lat = 48.110833, lon = 16.570833)
stations <- sybasestations("*")

stations$d <- sqrt((stations$lon - loww$lon)^2 + (stations$lat-loww$lat)^2)
stations <- subset(stations,d < .4 & tawes == 1)

params <- c("dd", "ff", "ffx", "tl", "rf", "p", "pred")

sybase_stationlist <- list()
for (i in 1:nrow(stations)) {
  sybase_file <- sprintf("data/_observations_sybase_tawes_around_vie_%d_%s_%s_%s.rds", 
    stations$statnr[i], start_date, end_date, version)
  tmp_info <- data.frame(statnr = NA, name = NA, file = NA)
  tmp_info$file <- sybase_file
  tmp_info$statnr <- stations$statnr[i]
  tmp_info$name <- make_nice_stationname(subset(stations, statnr == statnr[i])$name)
  sybase_stationlist[[i]] <- tmp_info
  if (file.exists(sybase_file)){
    next
  } else{ 
    cat(sprintf("File \"%s\" does not exists, prepare data ...\n", sybase_file))
 
    tmp <- try(sybase("tawes", stations$statnr[i], parameter = params,
      begin = as.POSIXct(start_date), end = as.POSIXct(end_date), archive = TRUE))
  
    if(class(try) == "try-error") next 

    saveRDS(file = sybase_file, tmp)
  }
}
sybase_stationlist <- do.call("rbind", sybase_stationlist)


# -------------------------------------------------------------------
# Read and combine data sets with temporal and spatial differences
# -------------------------------------------------------------------
cat("\nPrepare data ...\n")

## Load or create (and save) data set
datafile <- sprintf("data/circforest_prepared_data_vie_lag%sh_%s_%s_%s.rds", mylag, start_date, end_date, version)

if(file.exists(datafile)) {

  cat(sprintf("File \"%s\" exists, loading data ...\n", datafile))
  data <- readRDS(datafile)

} else {

  file_names <- sybase_stationlist[, "file"]
  
  ## Reduce number of stations
  file_names <- file_names[!grepl("bruckneudorf|brunn_am_gebirge|eisenstadt|neusiedl_am_see|wien_stammersdorf|
    wien_mariabrunn|zwerndorf_marchegg|altenburg|baden|mistelbach|podersdorf", sybase_stationlist[, "name"])]

  ## Read response and make strict zoo object, plus lagged response for spatial differences
  response <- eval(parse(text = sprintf("response <- make_strict(readRDS('%s'))", file_names[grep(response_station, file_names)])))
  if(unique(diff(index(response))) != 10) stop("wrong temporal resolution, suspected to 10min intervals")

  ## Subset to full hours
  response <- response[as.POSIXlt(index(response))$min == 0L, ]

  ## Subset to period between start_date and end_date
  response <- window(response, start = as.POSIXct(start_date), end = as.POSIXct(end_date))

  ## Calculate lagged response for spatial differences
  eval(parse(text = sprintf("response_lag <- lag(response, %s)", -mylag)))


  # -------------------------------------------------------------------
  # Prepare tawes data
  # -------------------------------------------------------------------
  stat_names <- NULL
  data <- NULL
  for(file in file_names){
    tmp <- sub('data/_observations_sybase_tawes_around_vie_', '', file)
    tmp <- sub("data/_observations_sybase_tawes_around_vie_", "", file)
    tmp <- sub(sprintf("_%s_%s_%s.rds", start_date, end_date, version), "", tmp)
    tmp <- make_nice_stationname(subset(stations, statnr == tmp)$name)

    cat(sprintf("\nPrepare data for station '%s'\n", tmp))
    
    eval(parse(text = sprintf("%s <- make_strict(readRDS(file))", tmp)))

    ## Check if 10min interval 
    eval(parse(text = sprintf("if(unique(diff(index(%s))) != 10) stop('wrong temporal resolution, suspected to 10min intervals')", tmp)))

    ## Subset to full hours
    eval(parse(text = sprintf("%1$s <- %1$s[as.POSIXlt(index(%1$s))$min == 0L, ]", tmp)))

    ## Subset to period between start_date and end_date
    eval(parse(text = sprintf("%1$s <- window(%s, start = as.POSIXct('%s'), end = as.POSIXct('%s'))", tmp, start_date, end_date)))

    ## Make lags
    eval(parse(text = sprintf("%1$s <- lag(%1$s, %2$s)", tmp, -mylag)))
    eval(parse(text = sprintf("%1$s <- make_strict(%1$s, index(response_lag))", tmp)))

    ## Calculate spatial differences
    if (!grepl(response_station, file)){
      for(i_param in names(eval(parse(text = tmp)))){
        if(grepl("dd", i_param)){
          eval(parse(text = sprintf("%1$s$dd.diff_resp <- calc_angle_dist(response_lag$dd, %1$s$dd, direction = TRUE, unit = 'deg')", 
          tmp)))
        } else if (i_param %in% c("p", "psta")){
          eval(parse(text = sprintf("%1$s$%2$s.diff_resp <- response_lag$p - %1$s$%2$s", tmp, i_param)))
        } else if (i_param %in% c("pred", "predsta")){
          eval(parse(text = sprintf("%1$s$%2$s.diff_resp <- response_lag$pred - %1$s$%2$s", tmp, i_param)))
        } else if (i_param %in% c("t", "tl")){
          eval(parse(text = sprintf("%1$s$%2$s.diff_resp <- response_lag$t - %1$s$%2$s", tmp, i_param)))
        } else if (i_param %in% c("rh", "rf")){
          eval(parse(text = sprintf("%1$s$%2$s.diff_resp <- response_lag$rf - %1$s$%2$s", tmp, i_param)))
        } else if (i_param %in% c("ffx", "ffx1s")){
          eval(parse(text = sprintf("%1$s$%2$s.diff_resp <- response_lag$ffx - %1$s$%2$s", tmp, i_param)))
        } else if (i_param %in% c("ff")){
          eval(parse(text = sprintf("%1$s$%2$s.diff_resp <- response_lag$ff - %1$s$%2$s", tmp, i_param)))
        }
      }
    }
 
    ## Get ff maximum, minimum and mean over last three hours
    for(i_param in names(eval(parse(text = tmp)))){
      if (i_param %in% c("ff", "ffx", "ffx1s")){
        eval(parse(text = sprintf("%1$s$%2$s.max3h <- rollapply(%1$s$%2$s, width = 3, FUN = max, fill = NA, align = 'right')", 
          tmp, i_param)))
        eval(parse(text = sprintf("%1$s$%2$s.min3h <- rollapply(%1$s$%2$s, width = 3, FUN = min, fill = NA, align = 'right')", 
          tmp, i_param)))
        eval(parse(text = sprintf("%1$s$%2$s.mean3h <- rollapply(%1$s$%2$s, width = 3, FUN = mean, fill = NA, align = 'right')", 
          tmp, i_param)))
      }
    }

    ## Calculate temporal changes
    for(i_param in names(eval(parse(text = tmp)))[!grepl("(min|max|mean|diff)", names(eval(parse(text = tmp))))]){
      if(grepl("dd", i_param)){
        eval(parse(text = sprintf("%1$s$%2$s.ch1h <- calc_angle_dist(%1$s$dd, lag(%1$s$dd, -1), direction = TRUE, unit = 'deg')", 
          tmp, i_param)))
        eval(parse(text = sprintf("%1$s$%2$s.ch3h <- calc_angle_dist(%1$s$dd, lag(%1$s$dd, -3), direction = TRUE, unit = 'deg')", 
          tmp, i_param)))
      } else {
        eval(parse(text = sprintf("%1$s$%2$s.ch1h <- %1$s$%2$s - lag(%1$s$%2$s, -1)", tmp, i_param)))
        eval(parse(text = sprintf("%1$s$%2$s.ch3h <- %1$s$%2$s - lag(%1$s$%2$s, -3)", tmp, i_param)))
      }
    }

    if(is.null(data)) {
      names(response) <- paste0("response.", names(response))
      eval(parse(text = sprintf("names(%1$s) <- paste0('%1$s.', names(%1$s))", tmp)))
      eval(parse(text = sprintf("data <- merge(response, %s)", tmp)))
    } else {
      eval(parse(text = sprintf("names(%1$s) <- paste0('%1$s.', names(%1$s))", tmp)))
      eval(parse(text = sprintf("data <- merge(data, %s)", tmp)))
    }  

    stat_names <- c(stat_names, tmp)
    eval(parse(text = sprintf("rm(%s)", tmp)))

  }

  # -------------------------------------------------------------------
  # Prepare lowvis data
  # -------------------------------------------------------------------
  cat(sprintf("\nPrepare data for station lowvis data'\n"))

  lowvis_data <- make_strict(readRDS(lowvis_file))

  ## Check if 10min interval 
  if(unique(diff(index(lowvis_data))) != 10) stop("wrong temporal resolution, suspected to 10min intervals")

  ## Subset to full hours
  lowvis_data <- lowvis_data[as.POSIXlt(index(lowvis_data))$min == 0L, ]

  # Make lags
  lowvis_data <- lag(lowvis_data, -mylag)
  lowvis_data <- make_strict(lowvis_data, index(response_lag))

  ## Calculate spatial differences (just for wind speed and direction)
  for(i_param in names(lowvis_data)[!grepl("(2min|_min|_max|ch|diff)", names(lowvis_data))]){
    if(grepl("dd", i_param)){
      eval(parse(text = sprintf("%1$s$%2$s.diff_resp <- calc_angle_dist(response_lag$dd, %1$s$%2$s, direction = TRUE, unit = 'deg')",
      "lowvis_data", i_param)))
    } else if (grepl("ff", i_param)) {
      eval(parse(text = sprintf("%1$s$%2$s.diff_resp <- response_lag$ff - %1$s$%2$s",
      "lowvis_data", i_param)))
    }
  }

  ## Get ff maximum, minimum and mean over last three hours
  for(i_param in names(lowvis_data)[!grepl("(2min|_min|_max|ch|diff)", names(lowvis_data))]){
    if (grepl("ff", i_param)){
      eval(parse(text = sprintf("%1$s$%2$s.max3h <- rollapply(%1$s$%2$s, width = 3, FUN = max, fill = NA, align = 'right')",
        "lowvis_data", i_param)))
      eval(parse(text = sprintf("%1$s$%2$s.min3h <- rollapply(%1$s$%2$s, width = 3, FUN = min, fill = NA, align = 'right')",
        "lowvis_data", i_param)))
      eval(parse(text = sprintf("%1$s$%2$s.mean3h <- rollapply(%1$s$%2$s, width = 3, FUN = mean, fill = NA, align = 'right')",
        "lowvis_data", i_param)))
    }
  }

  ## Calculate temporal changes (just for 10min_mean wind variables and other variables)
  for(i_param in names(lowvis_data)[!grepl("(2min|_min|_max|ch|diff|mean3h|max3h|min3h)", names(lowvis_data))]){
    if(grepl("dd", i_param)){
      eval(parse(text = sprintf("%1$s$%2$s.ch1h <- calc_angle_dist(%1$s$%2$s, lag(%1$s$%2$s, -1), direction = TRUE, unit = 'deg')",
        "lowvis_data", i_param)))
      eval(parse(text = sprintf("%1$s$%2$s.ch3h <- calc_angle_dist(%1$s$%2$s, lag(%1$s$%2$s, -3), direction = TRUE, unit = 'deg')",
        "lowvis_data", i_param)))
    } else {
      eval(parse(text = sprintf("%1$s$%2$s.ch1h <- %1$s$%2$s - lag(%1$s$%2$s, -1)", "lowvis_data", i_param)))
      eval(parse(text = sprintf("%1$s$%2$s.ch3h <- %1$s$%2$s - lag(%1$s$%2$s, -3)", "lowvis_data", i_param)))
    }
  }

  data <- merge(data, lowvis_data)

  cat(sprintf("File \"%s\" does not exists, saving data ...\n", datafile))
  saveRDS(data, datafile)

}
