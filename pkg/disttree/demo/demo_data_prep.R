setwd("~/svn/partykit/pkg/disttree/inst/draft")
#setwd("~/disttree/demo")

########################################
# prepare dataset for demo

{
  # load observations
  #load("~/svn/partykit/pkg/disttree/inst/draft/rainData/rain.rda")
  #load("~/disttree/inst/draft/rainData/rain.rda")
  load("rainData/rain.rda")
  #head(data)
  Sys.setenv("TZ" = "UTC")
  rain  <- data
  rm(data)
  #dim(rain)
  rain$robs <- (rain$obs)^(1/1.6)
  rain$day <- as.POSIXlt(rain$date)$yday
  rain$year <- as.POSIXlt(rain$date)$year
  rain$hour <- as.POSIXlt(rain$date)$hour
  #head(rain)
  ## save complete data
  rain_all <- rain
  
  
  # load station list
  #load("~/svn/partykit/pkg/disttree/inst/draft/rainData/ehyd.statlist.rda")
  #load("~/disttree/inst/draft/rainData/ehyd.statlist.rda")
  load("rainData/ehyd.statlist.rda")
  stations <- data.frame(ehyd.statlist$station, ehyd.statlist$stationname, 
                         ehyd.statlist$lon, ehyd.statlist$lat, ehyd.statlist$height)
  colnames(stations) <- c("stationnr", "stationname",
                          "lon", "lat", "height")
  
  
  
  ### check which stations provide observations for the years 1985 - 2012 (both included)
  {
    complete <- list()
    nryears <- numeric(length = NROW(stations))
    all85_112 <- logical(length = NROW(stations))
    all31 <- logical(length = NROW(stations))
    
    for(i in 1:NROW(stations)){
      
      rain <- rain_all
      stationname <- as.character(stations[i, "stationname"])
      stationnr <- stations[i, "stationnr"]
      
      # load predictions
      load(paste0("rainData/prepared/GEFSV2_prepared_", stationnr, ".rda", sep = ""))
      prediction <- prepared
      Sys.setenv("TZ" = "UTC")
      prediction$day <-as.POSIXlt(prediction$init)$yday
      prediction$year <- as.POSIXlt(prediction$init)$year
      
      rain <- rain[rain$station == stationname,]
      
      ## choose one month:
      # July
      rain <- rain[(((rain$year %% 4 == 0) & rain$day>=183) | (!(rain$year %% 4 == 0) & rain$day>=182)),]
      rain <- rain[(((rain$year %% 4 == 0) & rain$day<=213) | (!(rain$year %% 4 == 0) & rain$day<=212)),]
      
      all85_112[i] <- all(c(85:112) %in% as.numeric(row.names(table(rain$year)))) 
      all31[i] <- (min(table(rain$year)) == 31) & (max(table(rain$year)) == 31)
      
      complete[[i]] <- table(rain$year)
      #unique(rain$year)
      nryears[i] <- length(unique(rain$year))
    }
    
    complete_stations <- which(all85_112 & all31)
  }
  
  rain <- rain_all
  
  stations_comp <- stations[complete_stations,]
  
  
  
  
  ## prepare rain data of all stations but only July and from 1985
  rain <- rain_all
  ## choose one month:
  # July
  rain <- rain[(((rain$year %% 4 == 0) & rain$day>=183) | (!(rain$year %% 4 == 0) & rain$day>=182)),]
  rain <- rain[(((rain$year %% 4 == 0) & rain$day<=213) | (!(rain$year %% 4 == 0) & rain$day<=212)),]
  
  # predictions start in 1985
  rain <- rain[rain$year>=85,]
  
  rain_all_July85 <- rain
  
  
  #######################
  # get covariates (numerical ensemble predictions) for each station
  rainlist <- list()
  
  for(i in 1:nrow(stations_comp)){
    rain <- rain_all_July85
    stationname <- as.character(stations[complete_stations[i], "stationname"])
    stationnr <- stations[complete_stations[i], "stationnr"]
    
    # load predictions
    #load(paste0("~/svn/partykit/pkg/disttree/inst/draft/rainData/prepared/GEFSV2_prepared_", stationnr, ".rda"))
    #load(paste0("~/disttree/inst/draft/rainData/prepared/GEFSV2_prepared_", stationnr, ".rda"))
    load(paste0("rainData/prepared/GEFSV2_prepared_", stationnr, ".rda", sep = ""))
    prediction <- prepared
    Sys.setenv("TZ" = "UTC")
    prediction$day <-as.POSIXlt(prediction$init)$yday
    prediction$year <- as.POSIXlt(prediction$init)$year
    #head(prediction)
    
    rain <- rain[rain$station == stationname,]
    
    # observations end in 2012
    prediction <- prediction[prediction$year<=112,]
    
    # combine data frames
    
    # predictions with init 07-01 are for day 07-02 
    # observations are the sum of precipitation of the last 24h
    # match prediction with init 07-01 with observation on 07-02 
    # (pred starts at 00UTC and predicts from 06UTC until 30UTC = 06UTC of next day, 
    # observations are meassured from 06UTC of previous day to 06UTC of current day)
    
    raindata <- cbind(rain$date, rain$obs, rain$robs, prediction)
    #head(raindata[,c(1:10)])
    colnames(raindata)[c(1:3)] <- c("date", "obs", "robs")
    
    
    
    ########################################
    # select variables with sufficient data
    
    #dim(raindata)
    #colnames(raindata)
    
    
    #remove row 825
    raindata <- raindata[-825,]
    
    #remove rows 223-310
    #raindata <- raindata[-c(223:310),]
    #remove rows 111 176 178 220 221 241 320 364 367
    #raindata <- raindata[-c(111, 176, 178, 220, 221, 241, 320, 364, 367),]
    #remove rows 5 215 217 301 314
    #raindata <- raindata[-c(5, 215, 217, 301, 314),]
    #remove rows 251 264 
    #raindata <- raindata[-c(251, 264),]
    
    # table(raindata[, "year"])
    
    ######################################################
    # only keep variables with sufficient values
    raindata <- raindata[, c("robs", "year",
                             "tppow_mean", "tppow_sprd", "tppow_min", "tppow_max", 
                             "tppow_mean0612", "tppow_mean1218", "tppow_mean1824", "tppow_mean2430", 
                             "tppow_sprd0612", "tppow_sprd1218", "tppow_sprd1824", "tppow_sprd2430",
                             "tp_frac", 
                             "tp_frac0612", "tp_frac1218", "tp_frac1824", "tp_frac2430", 
                             "capepow_mean", "capepow_sprd", "capepow_min", "capepow_max", 
                             "capepow_mean0612", "capepow_mean1218", "capepow_mean1224", "capepow_mean1230",
                             "capepow_sprd0612", "capepow_sprd1218", "capepow_sprd1224", "capepow_sprd1230",
                             "cape_frac",
                             "cape_frac0612", "cape_frac1218", "cape_frac1824", "cape_frac2430", 
                             "dswrf_mean_mean", "dswrf_mean_min", "dswrf_mean_max",
                             "dswrf_sprd_mean", "dswrf_sprd_min", "dswrf_sprd_max",
                             "msl_mean_mean", "msl_mean_min", "msl_mean_max",
                             "msl_sprd_mean", "msl_sprd_min", "msl_sprd_max",
                             "pwat_mean_mean", "pwat_mean_min", "pwat_mean_max",
                             "pwat_sprd_mean", "pwat_sprd_min", "pwat_sprd_max",
                             "tcolc_mean_mean", "tcolc_mean_min", "tcolc_mean_max",
                             "tcolc_sprd_mean", "tcolc_sprd_min", "tcolc_sprd_max",
                             "tmax_mean_mean", "tmax_mean_min", "tmax_mean_max",
                             "tmax_sprd_mean", "tmax_sprd_min", "tmax_sprd_max",
                             "t500_mean_mean", "t500_mean_min", "t500_mean_max",
                             "t700_mean_mean", "t700_mean_min", "t700_mean_max",
                             "t850_mean_mean", "t850_mean_min", "t850_mean_max",
                             "t500_sprd_mean", "t500_sprd_min", "t500_sprd_max",
                             "t700_sprd_mean", "t700_sprd_min", "t700_sprd_max",
                             "t850_sprd_mean", "t850_sprd_min", "t850_sprd_max",
                             "tdiff500850_mean", "tdiff500850_min", "tdiff500850_max",
                             "tdiff700850_mean", "tdiff700850_min", "tdiff700850_max",
                             "tdiff500700_mean", "tdiff500700_min", "tdiff500700_max")]
    
    ######################################################
    #new vairable msl_diff
    raindata$msl_diff <- raindata$msl_mean_max - raindata$msl_mean_min
    
    ## FIX ME: check years (85-112)
    # table(raindata$year)
    rainlist$raindata <- raindata
    names(rainlist)[length(rainlist)] <- strsplit(stationname, split = " ")[[1]][1]
  }
}

# save data
{
  stations <- stations_comp
  save(stations, file = "~/svn/partykit/pkg/RainTyrol/data/stations.rda")
  rainTyrol <- rainlist
  save(rainTyrol, file = "~/svn/partykit/pkg/RainTyrol/data/rainTyrol.rda")
  rainAxams <- rainlist$Axams
  save(rainAxams, file = "~/svn/partykit/pkg/disttree/data/rainAxams.rda")
}



