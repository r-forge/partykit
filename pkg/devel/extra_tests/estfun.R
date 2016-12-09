
library("mlt")
library("sandwich")
library("libcoin")
set.seed(29)

n <- 100
x <- runif(n)
ct <- .5
teststat <- c("maximum", "quadratic")

for (ts in teststat) {

print(ts)

ret <- matrix(NA, ncol = 4, nrow = 100) 

for (i in 1:NROW(ret)) {
mydf <- data.frame(y = rnorm(n, mean = c(0, 1)[(x > ct) + 1], 
                                sd = c(1, 2)[(x > ct) + 1]),
                   x = x)

vy <- numeric_var("y", support = quantile(mydf$y, prob = c(.1, .9)))
by <- as.basis(~ y, data = vy)
by <- Bernstein_basis(vy, ui = "incre", order = 1)

m <- ctm(by)
fm <- mlt(m, data = mydf)

Y <- estfun(fm)
xct <-0:50/50
ix <- findInterval(x, vec = xct)

lstat <- LinStatExpCov(ix = ix, Y = Y)
tst <- doTest(lstat, teststat = ts)
ret[i, 1] <- xct[tst$index]

cov(Y)
cor(Y)

rt <- t(chol(n * vcov(fm)))
Ys <- t(rt %*% t(Y))
cov(Ys)
cor(Ys)

lstat <- LinStatExpCov(ix = ix, Y = Ys)
tst <- doTest(lstat, teststat = ts)
ret[i, 2] <- xct[tst$index]

lstat <- LinStatExpCov(ix = ix, Y = mydf$y)
tst <- doTest(lstat, teststat = ts)
ret[i, 3] <- xct[tst$index]

lstat <- LinStatExpCov(ix = ix, Y = cbind(Y, mydf$y))
tst <- doTest(lstat, teststat = ts)
ret[i, 4] <- xct[tst$index]



}

print(colSums((ret - ct)^2))
print(summary(ret[,1] - ret[,2]))
print(summary(ret[,1] - ret[,3]))
print(summary(ret[,1] - ret[,4]))

}

