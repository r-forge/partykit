# -------------------------------------------------------------------
# - NAME:   circular_crps_simulation.R
# - AUTHOR: Moritz N. Lang
# - DATE:   2019-10-01
# -------------------------------------------------------------------
# - PURPOSE:
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2019-10-09 on thinkmoritz
# -------------------------------------------------------------------

## Load libraries
library("disttree")
library("circtree")
library("verification") # Careful, version 1.35 needed!!
library("raster")
library('sp')
library("colorspace")
library("latex2exp")

## Create folder for outputs
if (! dir.exists("results")) dir.create("results")

## Helper function to plot grid search
toraster <- function(x, which, brks = NULL, brks.col = NULL, type = "absolute", mtext  = NULL, ...) {
  xylab <- names(x)[grep("^mu|^kappa", names(x))]
  ra <- rasterFromXYZ(as.matrix(x[,c(xylab[1], xylab[2], which)]))
  if(any(is.null(c(brks, brks.col))) & (type == "absolute")) {
    brks <- pretty(x[, which], n = 8)
    if(all(brks >= 0) & brks[1] == 0) brks <- c(-brks[2], brks)
    brks.col <- rev(heat_hcl(length(brks), h = c(300, 200), c = c(60, 0), 
      l = c(25, 95), power = c(0.7, 1.3)))
  } else if(any(is.null(c(brks, brks.col))) & (type == "relative")){
    brks <- pretty(x[, which], n = 8)
    tmp <- seq(-max(abs(brks)), max(abs(brks)), diff(brks)[1])
    brks.col <-  diverge_hcl(length(tmp) - 1, h = c(260, 0), c = 80, 
      l = c(30, 90), power = 1.)[round(tmp,1) %in% round(brks,1)]
  }

  xylab[grep("^kappa", xylab)] <- paste0("log(", xylab[grep("^kappa", xylab)], ")")
  ##xylab[grep("^kappa", xylab)] <- paste0("", xylab[grep("^kappa", xylab)], "")

  plot(ra, breaks = brks, legend = FALSE, col = brks.col, xlab = xylab[1], ylab = xylab[2], ...)
  plot(ra, breaks = brks, legend.only = TRUE, 
    smallplot = c(0.8, 0.85, 0.3, 0.75), col = brks.col, 
    legend.args = list(text = ifelse(type == "absolute", "CRPS [rad]", "CRPS [%]"), side = 4, cex = 0.8,
    line = 3))
  if(!is.null(mtext)) mtext(side=3, line=0.25, cex=0.7, mtext)
}

for(i_obs in c(-pi, -pi/2, 0, pi/2, pi)){
  ## Create grid
  grid <- expand.grid(mu = seq(-pi, pi, by = 5 * 2 * pi / 360), kappa = seq(0.1, 6.1, by = 0.4))
  ##grid <- expand.grid(mu = seq(-pi, pi, by = 5 * 2 * pi / 360), kappa = seq(0.1, 2, by = 0.12))

  ## Calculate crps  
  crps <- NULL
  for(i in 1:nrow(grid)){
    grid$crps_grimit[i] <- as.numeric(crps.circ(x = i_obs, mu = grid[i, "mu"], kappa = exp(grid[i, "kappa"])))
    grid$crps_own[i] <- crps_vonmises(y = i_obs, mu = grid[i, "mu"], kappa = exp(grid[i, "kappa"]), sum = FALSE)
    ##grid$crps_grimit[i] <- as.numeric(crps.circ(x = i_obs, mu = grid[i, "mu"], kappa = grid[i, "kappa"]))
    ##grid$crps_own[i] <- crps_vonmises(y = i_obs, mu = grid[i, "mu"], kappa = grid[i, "kappa"], sum = FALSE)
  }
  
  ## Calculate crps percentage change
  grid$crps_pchange <- (grid$crps_own - grid$crps_grimit) / abs(grid$crps_grimit) * 100
  grid$crps_diff <- grid$crps_own - grid$crps_grimit
  
  pdf(file = sprintf("results/circular_crps_simulation_obs_%s%2.1f_grid.pdf", ifelse(sign(i_obs) == -1, "neg", "pos"), 
    abs(i_obs)), width = 8, height = 6)
  ##pdf(file = sprintf("results/circular_crps_simulation_obs_%s%2.1f_grid_smallkappa.pdf", ifelse(sign(i_obs) == -1, "neg", "pos"), 
  ##  abs(i_obs)), width = 8, height = 6)
  par(mfrow = c(2,2), mar = c(5, 5, 4, 6) + 0.1, oma = c(0, 1, 3, 2.5))
  toraster(grid, which = "crps_grimit", main = "CRPS 'Grimit et al. (2006)'", asp = 0)
  toraster(grid, which = "crps_own", main = "CRPS 'Characteristic Function'", asp = 0)
  toraster(grid, which = "crps_diff", main = "CRPS Difference", mtext = TeX('$CRPS_{CharFun} - CRPS_{Grimit}$'), asp = 0)
  toraster(grid, which = "crps_pchange", main = "CRPS Percentage Change", 
    mtext = TeX('$(CRPS_{CharFun} - CRPS_{Grimit}) / |CRPS_{Grimit}| \\times 100$'), type = "relative", asp = 0)
  #mtext(sprintf("Circular CRPS Comparison (observation = %2.2f)", i_obs), side = 3, line = 1, outer = TRUE, cex = 1.5)
  dev.off()
}


# -------------------------------------------------------------------
# Comparison of different charFuns
# -------------------------------------------------------------------

library("CharFun")
library("verification")
library("scoringRules")
library("CharFun")
library("latex2exp")

angle_dist <- function(a, b) {
  d <- abs(b - a)
  return(pmin(d, 2*pi - d))
}

## CRPS von Mises
crps_vonmises_new <- function(y, mu, kappa, sum = FALSE, na.rm = FALSE) {

  ## Perform some input checks
  if(any(y < -pi) || any(y > pi) || any(mu < -pi) || any(mu > pi) || any(kappa < 0 ))
    stop("y and mu must be in the interval of [-pi, pi], and kappa must be non negative!")

  if(!inherits(y, c("numeric", "integer")) || !inherits(mu, c("numeric", "integer")) ||
    !inherits(kappa, c("numeric", "integer"))) {
    stop("Input 'y', 'mu', and 'kappa' must be numeric vectors...")
  }

  ## Create data.frame (fails if any has not the same length or length equal one)
  dat <- data.frame("y" = y, "mu" = mu, "kappa" = kappa)

  ## Get mu with smallest absolute distance
  idx <- apply(abs(matrix(dat$mu, ncol = 3, nrow = nrow(dat), byrow = FALSE) +
    matrix(c(-2 * pi, 0, 2 * pi), ncol = 3, nrow = nrow(dat), byrow = TRUE) - dat$y), 1, which.min)
  dat$mu <- dat$mu + c(-2 * pi, 0, 2 * pi)[idx]

  ## Calculate CRPS based on characteristic function  
  rval <- sapply(1:nrow(dat), function(i){
    if(dat[i, "kappa"] > 1500){
      crps <- scoringRules::crps_norm(y = dat[i, "y"], mean = dat[i, "mu"], sd = sqrt(1 / dat[i, "kappa"]))
      return(crps)

    } else {
      n <- 100000
      crps <- mean(angle_dist(rvm(n, dat[i, "mu"], dat[i, "kappa"]), dat[i, "y"])) - 
        0.5 * mean(angle_dist(rvm(n, dat[i, "mu"], dat[i, "kappa"]), rvm(n, dat[i, "mu"], dat[i, "kappa"])))
      return(crps)
    }
  })

  if(sum) rval <- sum(rval, na.rm = na.rm)
  return(rval)
}


## Define characteristic functions
charFun_normal <- function(t, mean, sd){
  szt <- dim(t)
  t <- c(t)
  cf <- exp((0+1i) * mean * t - 1/2 * sd^2 * t^2)
  cf[t == 0] <- 1
  dim(cf) <- szt
  return(cf)
}

charFun_vonmises <- function(t, mu, kappa) {
  besselI(kappa, abs(t), expon.scaled = TRUE) / besselI(kappa, 0, expon.scaled = TRUE) * 
    exp(complex(real = 0, imaginary = (t * mu)))
}

charFun_unif <- function(t, a, b) {
  (exp((0+1i) * t * b) - exp((0+1i) * t * a)) / (complex(real = 0, imaginary = 1) * t * (b - a))
}


## Define integrals for crps 
int_vonmises <- function(x) {
  1 / pi * abs(charFun_vonmises(x, mu, kappa) -
    exp((0+1i) * x * y))^2 / x^2
}

int_normal <- function(x) {
  1 / pi * abs(charFun_normal(x, mu, sqrt(1 / kappa)) -
    exp((0+1i) * x * y))^2 / x^2
}

int_unif <- function(x) {
  1 / pi * abs(charFun_unif(x, 0, 2 * pi) -
    exp((0+1i) * x * y))^2 / x^2
}


## Calculate crps comparison for different kappas
y <- 0
mu <- pi
kappa_list <- c(.Machine$double.eps, 0.05, 0.1, 0.5, 1, 2, 10, 200, 1000, 1500)

crps_list <- list()
for(kappa in kappa_list){

  crps_list[[as.character(kappa)]] <- list("N_charfun" = try(integrate(int_normal, 0, Inf)$value),
                       "N_scoringRules" = try(scoringRules::crps_norm(y, mu, sqrt(1 / kappa))),
                       "vM_charfun" = try(integrate(int_vonmises, 0, Inf)$value),
                       "vM_grimit" = try(verification::crps.circ(y, mu, kappa)),
                       "Unif_charfun" = try(integrate(int_unif, 0, Inf)$value),
                       "Unif_scoringRules" = try(scoringRules::crps_unif(y, 0, 2 *pi)))
  
  crps_list[[as.character(kappa)]] <- sapply(crps_list[[as.character(kappa)]], 
    function(x) if(any(class(x) %in% "try-error")){NA} else {as.numeric(x)})
}

crps_scores <- do.call(rbind, crps_list)

# -------------------------------------------------------------------
# Plot comparison of curves using charFuns
# -------------------------------------------------------------------
my_col <-  colorspace::qualitative_hcl(4)

## Fix observations and mu
y <- 0
mu <- pi
  
for(kappa in kappa_list){
  crps_tmp <- list("N_charfun" = try(integrate(int_normal, 0, Inf)$value),
                   "N_scoringRules" = try(scoringRules::crps_norm(y, mu, sqrt(1 / kappa))),
                   "vM_charfun" = try(integrate(int_vonmises, 0, Inf)$value),
                   "vM_grimit" = try(verification::crps.circ(y, mu, kappa)),
                   "Unif_charfun" = try(integrate(int_unif, 0, Inf)$value),
                   "Unif_scoringRules" = try(scoringRules::crps_unif(y, 0, 2 *pi)))
  
  crps_tmp <- sapply(crps_tmp, 
    function(x) if(any(class(x) %in% "try-error")){NA} else {round(as.numeric(x), 2)})
  
  pdf(file = sprintf("results/circular_crps_simulation_charfuns_kappa%.2f.pdf", kappa), 
    width = 10, height = 7)
  par(mfrow = c(1, 4))

  curve(int_vonmises, from = 0, to = 5, col = my_col[1], lwd = 1.5, main = "vonMises Distribution")
  legend("topright", legend = c(TeX(sprintf('$CRPS_{CharFun} = %.2f$', crps_tmp["vM_charfun"])), 
    TeX(sprintf('$CRPS_{Grimit} = %.2f$', crps_tmp["vM_grimit"]))))
  
  curve(int_normal, from = 0, to = 5, col = my_col[2], lwd = 1.5, main = "Normal Distribution")
  legend("topright", legend = c(TeX(sprintf('$CRPS_{CharFun} = %.2f$', crps_tmp["N_charfun"])), 
    TeX(sprintf('$CRPS_{ScoringRules} = %.2f$', crps_tmp["N_scoringRules"]))))
  
  curve(int_unif, from = 0, to = 5, col = my_col[3], lwd = 1.5, main = "Uniform Distribution")
  legend("topright", legend = c(TeX(sprintf('$CRPS_{CharFun} = %.2f$', crps_tmp["Unif_charfun"])), 
    TeX(sprintf('$CRPS_{ScoringRules} = %.2f$', crps_tmp["Unif_scoringRules"]))))
  
  curve(int_vonmises(x) - int_unif(x), from = 0, to = 5, col = my_col[3], lwd = 1.5, main = "Differences")
  curve(int_vonmises(x) - int_normal(x), from = 0, to = 5, col = my_col[2], lwd = 1.5, add = TRUE)
  legend("topright", legend = c(TeX("$Curve_{vM} - Curve_{N}$"), TeX("$Curve_{vM} - Curve_{Unif}$")), 
    col = c(my_col[2], my_col[3]), pch = c(20, 20))
  
  dev.off()
}
