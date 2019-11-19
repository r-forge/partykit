# -------------------------------------------------------------------
# - NAME:   circforest_plot_exampletree.R
# - AUTHOR: Moritz N. Lang
# - DATE:   2019-11-12
# -------------------------------------------------------------------
# - PURPOSE:
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2019-11-19 on thinkmoritz
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Preliminaries
# -------------------------------------------------------------------
## Set time to utc (just in case)
Sys.setenv('TZ'='UTC')

## Load packages
library("optparse")
library("zoo")
library("disttree")
library("verification") # Careful, version 1.35 needed!!
library("circglmbayes")
devtools::load_all(".")

## Create folder for outputs
if (! dir.exists("results")) dir.create("results")

## Small helper functions

## Exponential weights
#calc_expweights <- function(theta, n){rev(sapply(1:n, function(x) theta^(x -1) * (1 - theta)))}
calc_expweights <- function(theta, n){rev(sapply(1:n, function(x) theta^(x -1) /
  ((theta^n - 1) /(theta - 1))))}

## Exponential centered weights
calc_expweights_centered <- function(theta, n){sapply(-n:n, function(x) theta^(abs(x) -1) /
  sum(theta^(abs(-n:n) -1)))}

ddff2uv <- function (dd, ff = NULL) {
    if (class(dd) %in% c("zoo", "data.frame")) {
        if (sum(!c("ff", "dd") %in% names(dd)) > 0)
            stop("necessary colums \"ff\" and/or \"dd\" missing")
        ff = as.numeric(dd$ff)
        dd = as.numeric(dd$dd)
        print("x")
    }
    else if (NCOL(dd) == 2) {
        ff <- dd[, 1]
        dd <- dd[, 2]
    }
    metrad <- dd * pi/180
    u <- ff * (-sin(metrad))
    v <- ff * (-cos(metrad))
    rad <- 2 * pi - metrad - pi/2
    rad <- ifelse(rad >= 2 * pi, rad - 2 * pi, rad)
    data.frame(u, v, rad)
}


# -------------------------------------------------------------------
# NAMELIST PARAMETERS
# -------------------------------------------------------------------
option_list <- list(
  make_option(c("-v", "--verbose"), action = "store_true", default = TRUE,
    help = "Print extra output [default]"),
  make_option(c("-q", "--quietly"), action = "store_false",
    dest = "verbose", help = "Print little output"),
  make_option("--run_name", type = "character", default = "v14",
    help = "Run name or version of script used for output name [default \"%default\"]"),
  make_option("--station", type = "character", default = "ibk",
    help = "Weather Station used for fitting (e.g., 'ibk', 'vie') [default \"%default\"]"),
  make_option("--lag", type = "integer", default = 1,
    help = "Lag in 10min time steps [default %default]"),
  make_option("--lowff", action = "store_true", default = FALSE,
    help = "Use wind speeds (ff) below 1 m/s (Warning: TRUE not meaningful for Vie) [default]"),
  make_option("--plot", action = "store_true", default = FALSE,
    help = "Plot validation [default]"),
  make_option(c("--seed"), type = "integer", default = 123,
    help="Set seed for bootstrapping [default %default]")
  )

opt <- parse_args(OptionParser(option_list = option_list))

if(opt$lag != 1 || opt$station != "ibk") stop("check namelist parameters...")


# -------------------------------------------------------------------
# Pre-process data
# -------------------------------------------------------------------
## Load data
tmp <- readRDS(sprintf("data/circforest_prepared_data_%s_lag%sh_2014-01-01_2018-12-31_update20191023.rds", opt$station, opt$lag))

if (opt$station == "ibk") {
  ## Subset data to five years and to full hours 
  tmp <- window(tmp, start = "2014-01-01", end = "2018-12-31")
  d <- tmp[as.POSIXlt(index(tmp))$min == 0L, ]; rm(tmp); gc()

  ## Remove very low values of ff 
  if(!opt$lowff){
    d <- d[!d[, "response.ff"] < 1,] ## ff below 1 not very meaningful, but at least correctly measured
  }

  ## Remove covariates with more than 10 nans
  idx_names <- names(which(apply(d, 2, function(x) sum(is.na(x))) > nrow(d) / 100 * 5)) ## remove sattelberg and ellboegen
  d <- d[, ! (names(d) %in% idx_names)]

  ## Omit nans
  d <- na.omit(d)

  ## Calculate u and v for lagged response (ibk)
  tmp <- ddff2uv(dd = as.numeric(d$innsbruck.dd), ff = as.numeric(d$innsbruck.ff))
  d$innsbruck.u <- tmp$u
  d$innsbruck.v <- tmp$v

  ## Transform response.dd from 0-360 degree to [-pi, pi]
  d$response.dd <- d$response.dd / 360 * 2*pi
  d$response.dd[d$response.dd > pi] <- d$response.dd[d$response.dd > pi] - 2*pi

} else if (opt$station == "vie") {
  ## Subset data to five years and to full hours
  tmp <- window(tmp, start = "2014-01-01", end = "2018-12-31")
  d <- tmp[as.POSIXlt(index(tmp))$min == 0L, ]; rm(tmp); gc()

  ## Remove very low values of ff
  if(!opt$lowff){
    d <- d[!d[, "response.ff"] < 1,] ## ff == 0.5 -> dd = NA; ff == 0 -> dd = 0
  }

  ## Remove covariates with more than 10 nans
  idx_names <- names(which(apply(d, 2, function(x) sum(is.na(x))) > nrow(d) / 100 * 5))
  d <- d[, ! (names(d) %in% idx_names)]

  ## Omit nans
  d <- na.omit(d)

  ## Calculate u and v for lagged response (windmessanlage 34)
  tmp <- ddff2uv(dd = as.numeric(d$wien_schwechat_flughafen.dd), ff = as.numeric(d$wien_schwechat_flughafen.ff))
  d$wien_schwechat_flughafen.u <- tmp$u
  d$wien_schwechat_flughafen.v <- tmp$v

  ## Transform response.dd from 0-360 degree to [-pi, pi]
  d$response.dd <- d$response.dd / 360 * 2*pi
  d$response.dd[d$response.dd > pi] <- d$response.dd[d$response.dd > pi] - 2*pi

} else {
  stop("Station not supported, currently station must be 'vie' or 'ibk'...")
}

# -------------------------------------------------------------------
# Plot example tree
# -------------------------------------------------------------------

## Transform zoo object to data.frame
timepoints <- index(d)
tmp_time <- as.POSIXlt(index(d))
d <- as.data.frame(d)

## Add day of the year and time
d$daytime <- cut(tmp_time$hour, seq(0,24,by=3), include.lowest = TRUE)
d$doy <- tmp_time$yday + 1

rm(tmp_time); gc()

d$steinach.pred.diff_resp <- round((d$innsbruck.pred - mean(d$innsbruck.pred)) - 
  (d$steinach.pred - mean(d$steinach.pred)), 1)
d_example <- subset(d, select = c(response.dd, innsbruck.dd, innsbruck.ff, landeck.pred.diff_resp,
    kufstein.pred.diff_resp, steinach.pred.diff_resp))

#par(mfrow = c(2, 1))
#plot(d$steinach.p[1:100], type = "l", col = 2, ylim = 
#  range(c(d$steinach.p[1:100], d$steinach.pred[1:100])), main = "Steinach")
#lines(d$steinach.pred[1:100], col = 3)
#legend("topleft", lty = 1, col = c(2, 3), legend = c(paste0("mean(p) = ", 
#  round(mean(d$steinach.p), 2)), 
#  paste0("mean(pred) = ", round(mean(d$steinach.pred), 2))))
#
#plot(d$innsbruck.p[1:100], type = "l", col = 2, ylim = 
#  range(c(d$innsbruck.p[1:100], d$innsbruck.pred[1:100])), main = "Innsbruck")
#lines(d$innsbruck.pred[1:100], col = 3)
#legend("topleft", lty = 1, col = c(2, 3), legend = c(paste0("mean(p) = ", 
#  round(mean(d$innsbruck.p), 2)), 
#  paste0("mean(pred) = ", round(mean(d$innsbruck.pred), 2))))

names(d_example)[names(d_example) == "innsbruck.dd"] <- "direction"
names(d_example)[names(d_example) == "innsbruck.ff"] <- "speed"
names(d_example)[names(d_example) == "landeck.pred.diff_resp"] <- "dpressure_west"
names(d_example)[names(d_example) == "kufstein.pred.diff_resp"] <- "dpressure_east"
names(d_example)[names(d_example) == "steinach.pred.diff_resp"] <- "dpressure_south"

## Set up formula
f <- as.formula(paste("response.dd ~ ", paste(names(d_example)[-grep("response", names(d_example))], collapse= "+")))


## Fit single tree for visualization
m_ct.plot <- circtree(formula = f,
                 data = d_example,
                 control = disttree_control(mincriterion = (1 - .Machine$double.eps),
                                            minbucket = 2000, maxdepth = 3))

## Plot single tree
pdf(file = sprintf("results/_plot_circforest_finalexampletree_%s_lag%s_%s%s.pdf",
  opt$station, opt$lag, ifelse(opt$lowff, "with_lowff_", ""), opt$run_name), width = 18, height = 10)
par(mar = c(3.1, 4.1, 2.1, 2.1))
plot(m_ct.plot, ep_args = list(justmin = 10), tp_args = list(template = "geographics"),
  ip_args = list(pval = FALSE), gp = gpar(cex = 1.4))
dev.off()

