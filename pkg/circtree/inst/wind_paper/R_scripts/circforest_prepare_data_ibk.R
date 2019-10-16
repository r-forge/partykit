# -------------------------------------------------------------------
# - NAME:   circforest_prepare_data_ibk.R
# - AUTHOR: Moritz N. Lang
# - DATE:   2019-05-27
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

data_in <- "/home/moritz/Data/circtree/Daten_Deborah/qc_agg10min"

mylag <- 3
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

calc_anglediff <- function(ff1, dd1, ff2, dd2) {
  z <- merge(ff1, dd1, ff2, dd2)
  #z <- na.omit(z)
  out <- data.frame(ff = as.numeric(z$ff1 - z$ff2))
  a <- as.numeric(z$dd1 - z$dd2)
  a <- matrix(a, ncol = 3, nrow = length(a), byrow = FALSE)
  b <- matrix(c(-360, 0, 360), byrow = TRUE, ncol =3, nrow = nrow(a))
  c <- (a - b)
  out$dd <- c[cbind(1:nrow(c), apply(c, 1, function(x) {ret <- which.min(abs(x)); if(length(ret) == 0) ret <- NA; return(ret)}))]
  #out$dd <- sapply(1:length(dd1), function(x) {tmp <- as.numeric(dd1[x] - dd2[x]) + c(-360, 0, 360); tmp[which.min(abs(tmp))]})

  out <- zoo(out, index(z))
  names(out) <- paste0("diff", "_", names(out))
  return(out)
}

# -------------------------------------------------------------------
# Read data sets
# -------------------------------------------------------------------
cat("\nPrepare data ...\n")

if (! dir.exists("data")) dir.create("data")
datafile <- sprintf("data/circforest_prepared_data_ibk_lag%sh_update20191015.rds", mylag)

## Load or create (and save) data set
if(file.exists(datafile)) {

  cat(sprintf("File \"%s\" exists, loading data ...\n", datafile))
  data <- readRDS(datafile)

} else {

  file_names <- list.files(data_in, pattern = "agg_10min_.*.rds", full.names = TRUE)
  file_names <- file_names[-grep("guetsch|obertauern|montana|mariazell|altdorf", file_names)]

  ## Read response and make strict zoo object, plus lagged response for spatial differences
  response <- eval(parse(text = sprintf("response <- make_strict(readRDS('%s'))", file_names[grep(response_station, file_names)])))
  if(unique(diff(index(response))) != 10) stop("wrong temporal resolution, suspected to 10min intervals")

  ## Subset to full hours
  response <- response[as.POSIXlt(index(response))$min == 0L, ]

  ## Calculate lagged response for spatial differences
  eval(parse(text = sprintf("response_lag <- lag(response, %s)", -mylag)))
 
  ## Read other data and make strict zoo object
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

    ## Make lags
    eval(parse(text = sprintf("%1$s <- lag(%1$s, %2$s)", tmp, -mylag)))
    eval(parse(text = sprintf("%1$s <- make_strict(%1$s, index(response_lag))", tmp)))

    ## Create differences to Innsbruck
    if (!grepl(response_station, file)){
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
    for(i_param in names(eval(parse(text = tmp)))){
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

  cat(sprintf("File \"%s\" does not exists, saving data ...\n", datafile))
  saveRDS(data, datafile)

}
