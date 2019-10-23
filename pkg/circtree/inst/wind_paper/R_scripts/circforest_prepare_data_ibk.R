# -------------------------------------------------------------------
# - NAME:   circforest_prepare_data_ibk.R
# - AUTHOR: Moritz N. Lang
# - DATE:   2019-05-27
# -------------------------------------------------------------------
# - PURPOSE: Create data file for IBK with temporal/spatial differences
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2019-10-23 on thinkmoritz
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Preliminary
# -------------------------------------------------------------------
## Load packages
library(zoo)
devtools::load_all("/home/moritz/Projects/profcast_project/scripts/RProfcast")

data_in <- "/home/moritz/Data/circtree/Daten_Deborah/qc_agg10min"

version <- "update20191023"
mylag <- 1 ## in hour
start_date <- "2014-01-01"
end_date <- "2018-12-31"
response_station <- "innsbruck"

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

## Helper function to make nice sensor names
make_nice_sensorname <- function(x) {
  x <- gsub(" ", "_", x)  
  x <- gsub("-", "", x)  
  x <- gsub("Windgeschwindigkeit", "ff", x)
  x <- gsub("Windrichtung", "dd", x)
  x <- gsub("Temperatur", "Temp", x)
  x <- gsub("Virtuell", "Virt", x)
  x <- gsub("relative_Feuchte", "rh", x)
  x <- gsub("Mittel", "mean", x)
  x <- gsub("Maximum", "max", x)
  x <- gsub("Minimum", "min", x)
  x <- gsub("10_min", "10min", x)
  x <- gsub("2_min", "2min", x)
  x <- tolower(x)
  return(x)
}

## Create data directory if necessary
if (!dir.exists("data")) dir.create("data")

# -------------------------------------------------------------------
# Read observations from databases using LOWVIS sql3 database for the observations 
#   (Profcast::observations)
# -------------------------------------------------------------------
lowvis_file <- sprintf("data/_observations_lowvisdb_ibk_%s_%s_%s.rds", start_date, end_date, version)
if (!file.exists(lowvis_file)) {

  cat(sprintf("File \"%s\" does not exists, prepare data ...\n", lowvis_file))

  ## Looping over sensors
  sensors <- readRDS("innsbruck_sensor_600.rds")
  
  lowvis_data <- list()
  for (i in 1:nrow(sensors)) {
    cat(sprintf("Loading observations for \"%s\"\n", sensors$ID[i]))
    tmp <- Profcast::observations(sensors$ID[i], station = "innsbruck", interval = "600", 
          date = c(as.Date(start_date), as.Date(end_date)), verbose = TRUE)
    tmp <- rmduplicated(tmp)

    ## Transform kt to m/s
    if (grepl("Windgeschwindigkeit", sensors$name[i])) {
    tmp <- tmp / 1.9438444924574
    }

    ## Merge and rename
    names(tmp) <- make_nice_sensorname(sensors$name[i])
    lowvis_data[[as.character(sensors$ID[i])]] <- tmp
  }

  lowvis_data <- do.call("merge", lowvis_data)
  attr(lowvis_data, "sensors") <- sensors$ID
  attr(lowvis_data, "created") <- Sys.time()
  attr(lowvis_data, "createdon") <- Sys.info()["nodename"]

  ## Modify prefix 
  pattern <- "^wma_[^_.]+_"
  name_postfix <- gsub(pattern, "", names(lowvis_data), perl = TRUE)
  name_prefix <- regmatches(names(lowvis_data), regexpr(pattern, names(lowvis_data)))
  name_prefix <- gsub("_$", "", name_prefix)
  names(lowvis_data) <- paste0(name_prefix, ".", name_postfix)

  ## Save lowvis data
  saveRDS(lowvis_data, file = lowvis_file)
}

# -------------------------------------------------------------------
# Read TAWES around airport from databeses at ZAMG (STAGE::sybase)
# -------------------------------------------------------------------
## COMMENT: Not necessary, files prepared and quality checked by Deborah


# -------------------------------------------------------------------
# Read and combine data sets with temporal and spatial differences
# -------------------------------------------------------------------
cat("\nPrepare data ...\n")

## Load or create (and save) data set
datafile <- sprintf("data/circforest_prepared_data_ibk_lag%sh_%s_%s_%s.rds", mylag, start_date, end_date, version)

if(file.exists(datafile)) {

  cat(sprintf("File \"%s\" exists, loading data ...\n", datafile))
  data <- readRDS(datafile)

} else {

  file_names <- list.files(data_in, pattern = "agg_10min_.*.rds", full.names = TRUE)
  file_names <- file_names[!grepl("guetsch|obertauern|montana|mariazell|altdorf", file_names)]

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
    tmp <- sub('.*agg_10min_', '', file)
    tmp <- sub('.rds', '', tmp)

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

  ## Calculate spatial differences (just for wind speed, direction and temperature)
  for(i_param in names(lowvis_data)[!grepl("(2min|_min|_max|ch|diff)", names(lowvis_data))]){
    if(grepl("dd", i_param)){
      eval(parse(text = sprintf("%1$s$%2$s.diff_resp <- calc_angle_dist(response_lag$dd, %1$s$%2$s, direction = TRUE, unit = 'deg')",
      "lowvis_data", i_param)))
    } else if (grepl("ff", i_param)) {
      eval(parse(text = sprintf("%1$s$%2$s.diff_resp <- response_lag$ff - %1$s$%2$s",
      "lowvis_data", i_param)))
    } else if (grepl("temp", i_param)) {
      eval(parse(text = sprintf("%1$s$%2$s.diff_resp <- response_lag$tl - %1$s$%2$s",
      "lowvis_data", i_param)))
    }
  }

  ## Get ff and temp maximum, minimum and mean over last three hours
  for(i_param in names(lowvis_data)[!grepl("(2min|_min|_max|ch|diff)", names(lowvis_data))]){
    if (grepl("(ff|temp)", i_param)){
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
