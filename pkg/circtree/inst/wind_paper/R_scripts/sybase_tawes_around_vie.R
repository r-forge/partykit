# -------------------------------------------------------------------
# - NAME:        sybase_tawes_around_vie.R
# - AUTHOR:      Reto Stauffer
# - DATE:        2018-01-04
# -------------------------------------------------------------------
# - DESCRIPTION:
# -------------------------------------------------------------------
# - EDITORIAL:   2018-01-04, RS: Created file on pc24-c707.
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2019-09-13 on thinkmoritz
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
  outfile <- sprintf("STAGEobs_tawes_wind_%d.rda", stations$statnr[i])
  if (file.exists(outfile)) next
  ##
  obs <- sybase("tawes", stations$statnr[i], parameter = params,
        begin = "2013-01-01", end = "2019-09-01", archive = TRUE)
  ##
  save(file = outfile, obs)
}

save(file = "data/sybase_tawes_around_vie_stationlist.rda", stations)

for (i in 1:nrow(stations)) {

  outfile <- sprintf("STAGEobs_tawes_wind_%d.rda", stations$statnr[i])
  load(outfile)

  stations$statnr[i] <- obs
  eval(parse(text = sprintf("station%s <- obs", attr(obs, "info")$station)))

  if (i == 1){
    name_old <- attr(obs, "info")$station
  } else if (i == 2){
    eval(parse(text = sprintf("tawes_around_vie <- merge(%s, %s)", name_old, attr(obs, "info")$station)))
  } else if (i > 2){
    eval(parse(text = sprintf("tawes_around_vie <- merge(tawes_around, %s)", attr(obs, "info")$station)))
  }

}
