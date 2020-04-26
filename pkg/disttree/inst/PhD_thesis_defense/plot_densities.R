set.seed(7)
nobs <- 1000
x <- c(1:nobs)/nobs
mu1 <- 0.7
sigma1 <- 0.05
mu2 <- 0.4
sigma2 <- 0.1
y1 <- rnorm(nobs, mu1, sigma1)
y2 <- rnorm(nobs, mu2, sigma2)
mu <- mean(c(y1,y2))
sigma <- sd(c(y1,y2))
yd <- rnorm(10000*nobs, mean = mu, sd = sigma)
yd1 <- rnorm(10000*nobs, mean = mu1, sd = sigma1)
yd2 <- rnorm(10000*nobs, mean = mu2, sd = sigma2)

setwd("~/svn/partykit/pkg/disttree/inst/defensio/")
par(mar=c(0,0,0,0))
min(c(y1, y2))
max(c(y1, y2))
width <- 11
lwd <- 14

## HCL palette
pal <- hcl(c(10, 128, 260, 290, 50), 100, 50)
names(pal) <- c("forest", "tree", "gamlss", "gamboostLSS", "EMOS")
pallight <- hcl(c(10, 128, 260, 290, 70), 100, 50, alpha = 0.25)
names(pallight) <- c("forest", "tree", "gamlss", "gamboostLSS", "EMOS")
transpgray <- rgb(0.190,0.190,0.190, alpha = 0.2)

pdf(file = "density_all.pdf", width = width)
hist(c(y1,y2), breaks = c(0:50)/50, freq = FALSE, 
     xaxt="n", yaxt="n", ann=FALSE, border = "white", col = "white",
     ylim = c(0,4.5))
lines(density(yd), col = pal["forest"], lwd = lwd)
box(lwd=1)
dev.off()

pdf(file = "density_all_hist.pdf", width = width)
hist(c(y1,y2), breaks = c(0:50)/50, freq = FALSE, 
     xaxt="n", yaxt="n", ann=FALSE, border = transpgray, col = transpgray,
     ylim = c(0,4.5))
lines(density(yd), col = pal["forest"], lwd = lwd)
box(lwd=1)
dev.off()

pdf(file = "density_1.pdf", width = width)
hist(y1, breaks = c(0:50)/50, freq = FALSE, 
     xaxt="n", yaxt="n", ann=FALSE, border = "white", col = "white",
     ylim = c(0,9))
lines(density(yd1), col = pal["forest"], lwd = lwd)
lines(x = c(-1:5)/10, y = rep.int(0, times = 7), col = pal["forest"], lwd = lwd)
lines(x = c(9:11)/10, y = rep.int(0, times = 3), col = pal["forest"], lwd = lwd)
box(lwd=1)
dev.off()

pdf(file = "density_1_hist.pdf", width = width)
hist(y1, breaks = c(0:50)/50, freq = FALSE, 
     xaxt="n", yaxt="n", ann=FALSE, border = transpgray[1], col = transpgray[1],
     ylim = c(0,9))
lines(density(yd1), col = pal["forest"], lwd = lwd)
lines(x = c(-1:5)/10, y = rep.int(0, times = 7), col = pal["forest"], lwd = lwd)
lines(x = c(9:11)/10, y = rep.int(0, times = 3), col = pal["forest"], lwd = lwd)
box(lwd=1)
dev.off()

pdf(file = "density_2.pdf", width = width)
hist(y2, breaks = c(0:50)/50, freq = FALSE, 
     xaxt="n", yaxt="n", ann=FALSE, border = "white", col = "white",
     ylim = c(0,9))
lines(density(yd2), col = pal["forest"], lwd = lwd)
lines(x = c(8:11)/10, y = rep.int(0, times = 4), col = pal["forest"], lwd = lwd)
box(lwd=1)
dev.off()

pdf(file = "density_2_hist.pdf", width = width)
hist(y2, breaks = c(0:50)/50, freq = FALSE, 
     xaxt="n", yaxt="n", ann=FALSE, border = transpgray, col = transpgray,
     ylim = c(0,9))
lines(density(yd2), col = pal["forest"], lwd = lwd)
lines(x = c(8:11)/10, y = rep.int(0, times = 4), col = pal["forest"], lwd = lwd)
box(lwd=1)
dev.off()

