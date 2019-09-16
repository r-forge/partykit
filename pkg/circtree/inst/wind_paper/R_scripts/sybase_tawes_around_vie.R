# -------------------------------------------------------------------
# - NAME:        sybase_tawes_around_vie.R
# - AUTHOR:      Reto Stauffer
# - DATE:        2018-01-04
# -------------------------------------------------------------------
# - DESCRIPTION:
# -------------------------------------------------------------------
# - EDITORIAL:   2018-01-04, RS: Created file on pc24-c707.
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2019-09-16 on thinkmoritz
# -------------------------------------------------------------------

library("STAGE")
library("zoo")

if (!dir.exists("data")) dir.create("data")

loww <- list(lat = 48.110833, lon = 16.570833)
stations <- sybasestations("*")

stations$d <- sqrt((stations$lon - loww$lon)^2 + (stations$lat-loww$lat)^2)

stations <- subset(stations,d < .5 & tawes == 1)

params <- c("dd", "ff", "ffx", "tl", "rf", "p", "pred")

for (i in 1:nrow(stations)) {
  ##
  outfile <- sprintf("data/STAGEobs_tawes_wind_%d.rds", stations$statnr[i])
  if (file.exists(outfile)) next
  ##
  obs <- sybase("tawes", stations$statnr[i], parameter = params,
        begin = "2013-01-01", end = "2019-09-01", archive = TRUE)
  ##
  saveRDS(file = outfile, obs)
}

saveRDS(file = "data/sybase_tawes_around_vie_stationlist.rds", stations)

data_list <- list()
for (i in 1:nrow(stations)) {

  outfile <- sprintf("data/STAGEobs_tawes_wind_%d.rds", stations$statnr[i])
  data_list[[i]] <- readRDS(outfile)
  names(data_list[[i]]) <- sprintf("%s.station%s", names(data_list[[i]]), stations$statnr[i])

}

data <-  do.call(merge, data_list)

saveRDS(file = "data/tawes_around_vie.rds", data)

