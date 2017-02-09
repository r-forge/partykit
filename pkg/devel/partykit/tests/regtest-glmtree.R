
library("partykit")
library("Formula")

set.seed(29)
n <- 1000
x <- runif(n)
z <- runif(n)
y <- rnorm(n, mean = x * c(-1, 1)[(z > 0.5) + 1], sd = 3)
d <- data.frame(y = y, x = x, z = z)

ctrl <- partykit:::ctree_control()
ctrl$family <- gaussian()
ctrl$splitflavour <- "exhaustive"

system.time(m1 <- glmtree2(y ~ x | z, data = d, control = ctrl))
m1

ctrl$splitflavour <- "ctree"

system.time(m2 <- glmtree2(y ~ x | z, data = d, control = ctrl))
m2

system.time(m3 <- glmtree(y ~ x | z, data = d))
m3
