# -------------------------------------------------------------------
# - NAME:   circforest_prepare_data_vie.R
# - AUTHOR: Moritz N. Lang
# - DATE:   2019-05-27
# -------------------------------------------------------------------
# - PURPOSE: Create data file for VIE with temporal/spatial differences
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2019-09-11 on thinkmoritz
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
# Load observations from databases
# Using LOWVIS mysql database for the observations (STAGE::lowvis)
# and the sqlite3 for the wind-observations (Profcast::observations).
# -------------------------------------------------------------------

# Getting lowvisDescription file
obsfile <- sprintf("%s/_pure_observations_lowvisdb.rda", data_in)
if (! file.exists(obsfile)) {

  # -----------------------------------------
  # Getting LOWVIS observations
  # -----------------------------------------
  # Loaded distinct sensorids from database directly
  sensors <- c(200,201,205,206,208,502,507,522,527,542,547,562,567,
               582,587,590,10,20,30,0,12,22,32,2,121,131,151,9002,420,
               9000,141,161,9006,9008,9009,9010,9007,602,607,610,614,622,
               627,630,634,9102,165,166,9104,164,314,334,344,364,901,911,
               921,931,941,951,9001,9111,9116,9129,9134,9211,9216,9229,
               9234,9011,9016,9029,9034,907,917,927,937,947,957,662,667,
               670,674,682,687,690,694,702,707,710,714,594,9103,9311,9312,
               9313,9314,9105,666,608,9043)

  # Loading year by year
  obs_pure <- list()
  for (yr in 2013:2017) {
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

  # -----------------------------------------
  # Getting ffmax_10min and ffmean_2min for WMA's
  # -----------------------------------------
  stations <- stations("LOWW", relevant = TRUE)
  # Have some duplicated entries in the DB which are not in the DB which is
  # a bit freaky but ok.
  rmduplicated <- function(x) {
    idx <- which(duplicated(index(x)))
    if (length(idx) > 0) x <- x[-idx, ]
    x
  }
  # Looping over stations
  ffobs <- list()
  for (i in 1:nrow(stations)) {
    if (is.na(stations$sensors[i])) next
    cat(sprintf("   Loading observations for \"%s\"\n", stations$statnr[i]))
    tmp <- Profcast::observations(unlist(stations$sensors[i]), type = "10min",
          date = c(as.Date("2013-01-01"), Sys.Date()), verbose = opt$verbose)
    tmp <- lapply(tmp, rmduplicated)
    # Remove all flagged values here
    tmp <- lapply(tmp, function(x) subset(x, flag == 0)$value)
    # merge and rename
    tmp <- do.call("merge", tmp)
    names(tmp) <- sprintf("%s.%s", as.character(stations$statnr[i]), names(tmp))
    ffobs[[as.character(stations$statnr[i])]] <- tmp
  }
  ffobs <- do.call("merge", ffobs)
  attr(ffobs, "created") <- Sys.time()
  attr(ffobs, "createdon") <- Sys.info()["nodename"]

  # Save both objects
  save(file = obsfile, obs_pure, ffobs)
} else {
  cat(sprintf("   File \"%s\" exists, load and continue\n", obsfile))
  tmp <- load(obsfile)
  if("res" %in% tmp) obs_pure <- res; rm(res)  # workaround as vector might be called "res"
}

# -------------------------------------------------------------------
# Prepate data
# -------------------------------------------------------------------
obs_pure <- make_strict(obs_pure)

# Create lagged observations (obs_pure)   
data <- lag(obs_pure, -mylag)

# Calculate lagged measurements
for(i_name in names(data)[-grep("vv|sws|lvp|cei|rvr|vis|dd|ff", names(data))]){
  cmd <- sprintf("data$%1$s_ch30min <- data$%1$s - lag(data$%1$s, -3)",
    i_name)
  eval(parse(text = cmd))
}

tmp_names <- names(data)[grep("^dd", names(data))]
tmp_names <- gsub('dd', '', tmp_names)
for(idx in tmp_names){
  diff <- calc_anglediff(data[, paste0("ff", idx)], data[, paste0("dd", idx)], 
    lag(data[, paste0("ff", idx)], -3), lag(data[, paste0("dd", idx)], -3), paste0(idx, "_ch30min"))
  data <- merge(data, diff, all = c(TRUE, FALSE))
}

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

data <- merge(obs_pure$dd16, data, all = c(TRUE, FALSE))
data <- merge(obs_pure$ff16, data, all = c(TRUE, FALSE))

diffEXB <- calc_anglediff(data$ffEXB, data$ddEXB, data$ff16, data$dd16, "EXB")
diffTOW <- calc_anglediff(data$ffTOW, data$ddTOW, data$ff16, data$dd16, "TOW")
diffOM29 <- calc_anglediff(data$ffOM29, data$ddOM29, data$ff16, data$dd16, "OM29")
diffARS <- calc_anglediff(data$ffARS, data$ddARS, data$ff16, data$dd16, "ARS")

data <- merge(data, diffEXB, diffTOW, diffOM29, diffARS)

tawes_around <- readRDS(paste0(data_in, "/tawes_stations_around_airport.rds"))  ## generated in Profcast git-repo
tawes_around <- lag(tawes_around, -mylag)

data <- merge(data, tawes_around, all = c(TRUE, FALSE))

for(idx in unique(regmatches(names(tawes_around), regexpr("station[0-9]{5}", names(tawes_around))))){
  diff <- calc_anglediff(data[, paste0("ff.", idx)], data[, paste0("dd.", idx)], data$ff16, data$dd16, idx)
  data <- merge(data, diff, all = c(TRUE, FALSE))
}

names(data)[names(data) %in% "obs_pure$dd16"] <- "dd.response"
names(data)[names(data) %in% "obs_pure$ff16"] <- "ff.response"

if (! dir.exists(data_out)) dir.create(data_out)
saveRDS(data, paste0(data_out, "/circforest_prepared_data_vie_lag", mylag, ".rds"))


#data <- as.data.frame(data)
#
## Remove all columns with too little data
#tmp <- which(colSums(is.na(data)) / nrow(data) > 0.1)
#if (length(tmp) > 0) data <- data[, -tmp]
#
## Remove NA's
#data <- na.omit(data)
#cat(sprintf("   We ended up with a matrix of size: %d x %d\n", nrow(data), ncol(data)))

