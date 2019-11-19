# -------------------------------------------------------------------
# - NAME:   circforest_plot_density.R
# - AUTHOR: Moritz N. Lang, Lisa Schlosser
# - DATE:   2019-11-19
# -------------------------------------------------------------------
# - PURPOSE:
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2019-11-19 on thinkmoritz
# -------------------------------------------------------------------

library("circular")
library("latex2exp")

set.seed(123)

coefs = c(pi, 3)
stack = 10
cex = 1
label = TRUE
circlab = NULL
polygons = TRUE
rug = TRUE
template = "mathematical"
zero = "0"
rotation = "counter"
response_range = c(0, 2 * pi)
circlab <- c('$0$', '$\\pi/2$', '$\\pi$','$3\\pi/2$')

# Create data
X <- rvonmises(500, coefs[1], coefs[2])

#par(mfrow = c(1,2))

## circular density
pdf(file = "density_circular.pdf", width = 6, height = 3)
par(mar = c(3.1, 2.1, 2.1, 2.1))

# Empty Plot
plot(NA, xlim = c(-2,1), ylim = c(-1,1),
  type = "l", axes = FALSE, xlab = "",ylab = "", asp = 1)

# Histogram
Xhist <- hist(X, plot = FALSE, breaks = seq(0, 360, stack) * pi / 180)
idx <- (Xhist$density != 0)
Xscaled <- Xhist$density / max(Xhist$density)


# Draw densities (lines or polygons)
if (! polygons) {
  segments(cos(Xhist$mids[idx]), sin(Xhist$mids[idx]),
           (Xscaled[idx] + 1) * cos(Xhist$mids[idx]),
           (Xscaled[idx] + 1) * sin(Xhist$mids[idx]))
  #points(1 * cos(Xhist$mids), 1 * sin(Xhist$mids), pch = 20, col = c(NA, gray(0.2, alpha = 0.7))[idx + 1]) 
} else {
  segwidth <- diff(Xhist$mids[1:2])     # width of the segment
  segn     <- ceiling(segwidth / 0.005) # level of detail for the curvature of the segment
  for (i in which(idx)) {
    x <- seq(Xhist$mids[i] - segwidth/2, Xhist$mids[i] + segwidth/2, length = segn)
    polygon(c(cos(x), (Xscaled[i] + 1) * cos(rev(x))),
            c(sin(x), (Xscaled[i] + 1) * sin(rev(x))), border = NA, col = "gray80")
  }
}

# Plot rugs
if(rug){
  segments(cos(X), sin(X),
           (0.9) * cos(X),
           (0.9) * sin(X), gray(0.2, alpha = 0.7))
}

# Plot circle
circ_ln <- seq(0, 360, 0.5) * pi / 180
lines(cos(circ_ln), sin(circ_ln))

# Plot ticks and labels
circ_sgm <- seq(0, 360, 45) * pi / 180
segments(.8 * cos(circ_sgm), .8 * sin(circ_sgm), 1 * cos(circ_sgm), 1 * sin(circ_sgm))
if(label & requireNamespace("latex2exp", quietly = TRUE)){
  text(.75, 0, latex2exp::TeX(circlab[1]), adj = c(1, 0.5))
  text(0, .75, latex2exp::TeX(circlab[2]), adj = c(0.5, 1))
  text(-.75, 0, latex2exp::TeX(circlab[3]), adj = c(0, 0.5))
  text(0, -.75, latex2exp::TeX(circlab[4]), adj = c(0.5, 0))
}

# Plot fitted density
if(!is.null(coefs)){
  arrows(0, 0, 0.5 * cos(coefs[1]), 0.5 * sin(coefs[1]), col = 2, length = 0, lwd = 2)
  #arrows(0, 0, 0.5 * cos(coefs[1]), 0.5 * sin(coefs[1]), col = 2, length = 0.1, code = 2, lwd = 2)
  circular::plot.function.circular(function(x)
    dvonmises(x, circular::circular(coefs[1]), coefs[2]), add = TRUE, col = 2, lty = 1, lwd = 1.5)
}
dev.off()



## linear plot
breaks <- seq(from = 0, to = 2*pi, by = segwidth)
pdf(file = "density_linear.pdf", width = 6, height = 3)
par(mar = c(3.1, 2.1, 2.1, 2.1))
hist(as.numeric(X), border = NA, col = "gray80", axes = F, xlab = "", probability = T,
     ylab = "", main = "", breaks = breaks)
axis(1, at=seq(0,2*pi, by=pi/2), line = -0.38, 
     labels=c(latex2exp::TeX(circlab[1]),
              latex2exp::TeX(circlab[2]),
              latex2exp::TeX(circlab[3]),
              latex2exp::TeX(circlab[4]),
              latex2exp::TeX("$2\\pi$")))
plot.function(function(x)
  dvonmises(x, circular::circular(coefs[1]), coefs[2]), col = 2, lty = 1, lwd = 1.5, 
  from = 0, to = 2 *pi, add = TRUE)
rug(as.numeric(X), ticksize=0.036, side=1, lwd=0.5)
rug(seq(0,2*pi, by=pi/2), ticksize=0.08, side=1, lwd=1, line = 0.4)
rug(pi, ticksize=0.08, side=1, lwd=1, line = -0.4, col = "red")
dev.off()


