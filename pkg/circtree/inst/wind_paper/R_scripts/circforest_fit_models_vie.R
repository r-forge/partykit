# -------------------------------------------------------------------
# - NAME:   circforest_fit_models_vie.R
# - AUTHOR: Moritz N. Lang, Lisa Schlosser
# - DATE:   2019-06-03
# -------------------------------------------------------------------
# - PURPOSE:
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2019-09-11 on thinkmoritz
# -------------------------------------------------------------------

## Set time to utc (just in case)
Sys.setenv('TZ'='UTC')

## Load packages
library("zoo")
library("circtree")
library("disttree")
library("verification") # Careful, version 1.35 needed!!

if (! dir.exists("results")) dir.create("results")

# -------------------------------------------------------------------
# Pre-process data
# -------------------------------------------------------------------

## Load data / Remove rows with NAs
tmp <- readRDS("data/circforest_prepared_data_vie_lag6.rds")

## Remove rows with more than 10 percent nans
idx_names <- names(which(apply(tmp, 2, function(x) sum(is.na(x))) > nrow(tmp) / 10))
tmp <- tmp[, ! (names(tmp) %in% idx_names)]

tmp <- na.omit(tmp)

## Subset data to full hours
d <- tmp[as.POSIXlt(index(tmp))$min == 0L, ]; rm(tmp); gc()

## Remove very low values of ff and zero wind direction
d <- d[!d[, "dd.response"] == 0,]
d <- d[!d[, "ff.response"] < 1,]

## Transform response.dd from 0-360 degree to [-pi, pi]
d$dd.response <- d$dd.response / 360 * 2*pi
d$dd.response[d$dd.response > pi] <- d$dd.response[d$dd.response > pi] - 2*pi


## Transform zoo object to data.frame
tmp_time <- as.POSIXlt(index(d))
d <- as.data.frame(d)

## Add day of the year and time
d$daytime <- cut(tmp_time$hour, seq(0,24,by=3), include.lowest = TRUE)
d$doy <- tmp_time$yday + 1

rm(tmp_time); gc()

# -------------------------------------------------------------------
# Tryouts
# -------------------------------------------------------------------
## Set up formula
f <- as.formula(paste("dd.response ~ ", paste(names(d)[-grep("response", names(d))], collapse= "+")))

## Fit tree
m_ct <- circtree(formula = f,
                 data = d,
                 control = disttree_control(mincriterion = (1 - .Machine$double.eps),
                                            minbucket = 2000))
