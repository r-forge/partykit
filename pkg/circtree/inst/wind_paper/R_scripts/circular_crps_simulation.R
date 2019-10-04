# -------------------------------------------------------------------
# - NAME:   circular_crps_simulation.R
# - AUTHOR: Moritz N. Lang
# - DATE:   2019-10-01
# -------------------------------------------------------------------
# - PURPOSE:
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2019-10-02 on thinkmoritz
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

  xylab[grep("^kappa", xylab)] <- paste0("", xylab[grep("^kappa", xylab)], "")
  plot(ra, breaks = brks, legend = FALSE, col = brks.col, xlab = xylab[1], ylab = xylab[2], ...)
  plot(ra, breaks = brks, legend.only = TRUE, 
    smallplot = c(0.8, 0.85, 0.3, 0.75), col = brks.col, 
    legend.args = list(text = ifelse(type == "absolute", "CRPS [rad]", "CRPS [%]"), side = 4, cex = 0.8,
    line = 3))
  if(!is.null(mtext)) mtext(side=3, line=0.25, cex=0.7, mtext)
}

for(i_obs in c(-pi, -pi/2, 0, pi/2, pi)){
  ## Create grid
  grid <- expand.grid(mu = seq(-pi, pi, by = 5 * 2 * pi / 360), kappa = seq(0.1, 2, by = 0.12))

  ## Calculate crps  
  crps <- NULL
  for(i in 1:nrow(grid)){
    grid$crps_grimit[i] <- as.numeric(crps.circ(x = i_obs, mu = grid[i, "mu"], kappa = grid[i, "kappa"]))
    grid$crps_own[i] <- crps_vonmises(y = i_obs, mu = grid[i, "mu"], kappa = grid[i, "kappa"], sum = FALSE)
  }
  
  ## Calculate crps percentage change
  grid$crps_pchange <- (grid$crps_own - grid$crps_grimit) / abs(grid$crps_grimit) * 100
  grid$crps_diff <- grid$crps_own - grid$crps_grimit
  
  pdf(file = sprintf("results/circular_crps_simulation_obs_%s%2.1f_grid_smallkappa.pdf", ifelse(sign(i_obs) == -1, "neg", "pos"), 
    abs(i_obs)), width = 8, height = 6)
  par(mfrow = c(2,2), mar = c(5, 5, 4, 6) + 0.1, oma = c(0, 1, 3, 2.5))
  toraster(grid, which = "crps_grimit", main = "CRPS 'Grimit et al. (2006)'", asp = 0)
  toraster(grid, which = "crps_own", main = "CRPS 'Characteristic Function'", asp = 0)
  toraster(grid, which = "crps_diff", main = "CRPS Difference", mtext = TeX('$CRPS_{CharFun} - CRPS_{Grimit}$'), asp = 0)
  toraster(grid, which = "crps_pchange", main = "CRPS Percentage Change", 
    mtext = TeX('$(CRPS_{CharFun} - CRPS_{Grimit}) / |CRPS_{Grimit}| \\times 100$'), type = "relative", asp = 0)
  mtext(sprintf("Circular CRPS Comparison (observation = %2.2f)", i_obs), side = 3, line = 1, outer = TRUE, cex = 1.5)
  dev.off()
}
