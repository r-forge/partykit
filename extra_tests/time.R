
library("partykit")

data("iris")

a <- ctree(Species ~ ., data = iris, control = ctree_control(nmax = Inf))


library("proftools")

s <- 1:nrow(iris)
Y <- model.matrix(~ Species - 1, data = iris)
pd <- profileExpr(for (i in 1:500) a$update(function(...) list(estfun = Y, converged = TRUE), weights = integer(0), 
                  subset = s))

funSummary(pd)
srcSummary(pd)
hotPaths(pd, total.pct = 10.0)
plotProfileCallGraph(pd)
calleeTreeMap(pd)
flameGraph(pd)


