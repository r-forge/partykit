
set.seed(29)

library("partykit")
library("inum")
library("randomForest")

ntree <- 100
mtry <- 3
nmax <- 25

data("iris")
oiris <- iris
i <- sample(1:nrow(iris), 10000, replace = TRUE)
iris <- iris[i,]
dim(iris)

system.time(rf <- randomForest(Species ~ ., data = iris, mtry = mtry))

system.time(cf1 <- partykit:::cforest(Species ~ ., data = iris, mtry = mtry, ntree = ntree))

system.time(cf2 <- partykit:::cforest(Species ~ ., data = iris, mtry = mtry, 
                                      nmax = nmax, ntree = ntree))

system.time(cf3 <- party:::cforest(Species ~ ., data = iris, 
                                   control = party:::cforest_unbiased(ntree = ntree, 
                                                                      mtry = mtry)))

system.time(p <- predict(rf, newdata = oiris, type = "prob"))
system.time(p1 <- predict(cf1, newdata = oiris, type = "prob"))
system.time(p2 <- predict(cf2, newdata = oiris, type = "prob"))
system.time(p3 <- predict(cf3, newdata = oiris, type = "prob"))
p3 <- do.call("rbind", p3)

max(abs(p - p1))
max(abs(p - p2))
max(abs(p - p3))

