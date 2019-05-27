# -------------------------------------------------------------------
# - NAME:   circforest_prepare_data_ibk.R
# - AUTHOR: Moritz N. Lang
# - DATE:   2019-05-27
# -------------------------------------------------------------------
# - PURPOSE: Create data file for IBK with temporal/spatial differences
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2019-05-27 on thinkmoritz
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Preliminary
# -------------------------------------------------------------------
## Load packages
library(zoo)

data_in <- "/home/moritz/Data/circtree/Daten_Deborah/qc_agg10min"
data_out <- "/home/moritz/Projects/profcast_project/scripts/Rdevelopment/data_out"

mylag <- 3

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

datafile <- sprintf("%s/circforest_prepared_data_ibk_lag%s.rda", data_out, mylag)

## Load or create (and save) data set
if(file.exists(datafile)) {

  cat(sprintf("File \"%s\" exists, loading data ...\n", datafile))
  dat_ibk <- readRDS(datafile)

} else {

  file_names <- list.files(data_in, pattern = "agg_10min_.*.rds", full.names = TRUE)

  ## Read innsbruck and make strict zoo object
  response <- eval(parse(text = sprintf("response <- make_strict(readRDS('%s'))", file_names[grep("innsbruck", file_names)])))
  eval(parse(text = sprintf("innsbruck_tmp <- lag(response, %s)", -mylag)))
 
  ## Read other data and make strict zoo object
  stat_names <- NULL
  dat_ibk <- NULL
  for(file in file_names){
    tmp <- sub('.*agg_10min_', '', file)
    tmp <- sub('.rds', '', tmp)

    cat(sprintf("\nPrepare data for station '%s'\n", tmp))
    
    eval(parse(text = sprintf("%s <- make_strict(readRDS(file))", tmp)))

    eval(parse(text = sprintf("%1$s <- lag(%1$s, %2$s)", tmp, -mylag)))
    eval(parse(text = sprintf("%1$s <- make_strict(%1$s, index(innsbruck_tmp))", tmp)))

    ## Create differences to Innsbruck
    for(i_param in names(eval(parse(text = tmp)))){
      if(grepl("dd", i_param) & "ff" %in% names(eval(parse(text = tmp)))){
        eval(parse(text = sprintf("%1$s$dd_diff <- calc_anglediff(innsbruck_tmp$ff, innsbruck_tmp$dd, 
          %1$s$ff, %1$s$dd)$diff_dd", tmp)))
      } else if (i_param %in% c("p", "psta")){
        eval(parse(text = sprintf("%1$s$%2$s_diff <- innsbruck_tmp$p - %1$s$%2$s", tmp, i_param)))
      } else if (i_param %in% c("pred", "predsta")){
        eval(parse(text = sprintf("%1$s$%2$s_diff <- innsbruck_tmp$pred - %1$s$%2$s", tmp, i_param)))
      } else if (i_param %in% c("t", "tl")){
        eval(parse(text = sprintf("%1$s$%2$s_diff <- innsbruck_tmp$t - %1$s$%2$s", tmp, i_param)))
      } else if (i_param %in% c("rh", "rf")){
        eval(parse(text = sprintf("%1$s$%2$s_diff <- innsbruck_tmp$rf - %1$s$%2$s", tmp, i_param)))
      } else if (i_param %in% c("ffx", "ffx1s")){
        eval(parse(text = sprintf("%1$s$%2$s_diff <- innsbruck_tmp$ffx - %1$s$%2$s", tmp, i_param)))
      } else if (i_param %in% c("ff")){
        eval(parse(text = sprintf("%1$s$%2$s_diff <- innsbruck_tmp$ff - %1$s$%2$s", tmp, i_param)))
      }
    }
 
    ## Get ffx maximum and mean over last hour
    for(i_param in names(eval(parse(text = tmp)))){
      if (i_param %in% c("ffx", "ffx1s")){
        eval(parse(text = sprintf("%1$s$%2$s_mean <- rollapply(%1$s$%2$s, width = 6, FUN = mean, fill = NA, align = 'right')", 
          tmp, i_param)))
        eval(parse(text = sprintf("%1$s$%2$s_max <- rollapply(%1$s$%2$s, width = 6, FUN = max, fill = NA, align = 'right')", 
          tmp, i_param)))
      }
    }

    ## Calculate temporal changes
    for(i_param in names(eval(parse(text = tmp)))){
      if(grepl("dd", i_param)){
        eval(parse(text = sprintf("%1$s$%2$s_ch30min <- calc_anglediff(%1$s$ff, %1$s$dd, 
          lag(%1$s$ff, -3), lag(%1$s$dd, -3))$diff_dd", tmp, i_param)))
        eval(parse(text = sprintf("%1$s$%2$s_ch60min <- calc_anglediff(%1$s$ff, %1$s$dd, 
          lag(%1$s$ff, -6), lag(%1$s$dd, -6))$diff_dd", tmp, i_param)))
        eval(parse(text = sprintf("%1$s$%2$s_ch90min <- calc_anglediff(%1$s$ff, %1$s$dd, 
          lag(%1$s$ff, -9), lag(%1$s$dd, -9))$diff_dd", tmp, i_param)))
      } else {
        eval(parse(text = sprintf("%1$s$%2$s_ch30min <- %1$s$%2$s - lag(%1$s$%2$s, -3)", tmp, i_param)))
        eval(parse(text = sprintf("%1$s$%2$s_ch60min <- %1$s$%2$s - lag(%1$s$%2$s, -6)", tmp, i_param)))
        eval(parse(text = sprintf("%1$s$%2$s_ch90min <- %1$s$%2$s - lag(%1$s$%2$s, -9)", tmp, i_param)))
      }
    }

    if(is.null(dat_ibk)) {
      names(response) <- paste0(names(response), ".response")
      eval(parse(text = sprintf("names(%1$s) <- paste0(names(%1$s), '.%1$s')", tmp)))
      eval(parse(text = sprintf("dat_ibk <- merge(response, %s)", tmp)))
    } else {
      eval(parse(text = sprintf("names(%1$s) <- paste0(names(%1$s), '.%1$s')", tmp)))
      eval(parse(text = sprintf("dat_ibk <- merge(dat_ibk, %s)", tmp)))
    }  

    stat_names <- c(stat_names, tmp)
    eval(parse(text = sprintf("rm(%s)", tmp)))

  }

  cat(sprintf("File \"%s\" does not exists, saving data ...\n", datafile))
  saveRDS(dat_ibk, datafile)

}
