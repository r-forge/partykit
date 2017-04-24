
library("partykit")
set.seed(29)

n <- 500
p <- 10
x <- as.data.frame(matrix(runif(n * p), ncol = p))
x$y <- rnorm(n)

ns <- 10

cf <- cforest(y ~ ., data = x, ntree = 100, mtry = 3, trace = TRUE,
    perturb = list(replace = TRUE, fraction = 0.632),
    control = ctree_control(minsplit = 2,
                            mincriterion = 0, minbucket = ns))

cfy <- predict(cf, OOB = FALSE)
cfyOOB <- predict(cf, OOB = TRUE)

mean((x$y - cfy)^2)
mean((x$y - cfyOOB)^2) ### this should be 1 
### gets closer to 1 as ns increases

nt <- table(predict(gettree(cf, 10), type = "node"))
length(nt) ### number of terminal nodes 
table(nt)

library("randomForest")

(rf <- randomForest(y ~ ., data = x, ntree = 100, mtry = 3, nodesize = ns))

rfy <- predict(rf, newdata = x)
rfyOOB <- predict(rf)

mean((x$y - rfy)^2) ### more overfitting than in cforest
mean((x$y - rfyOOB)^2)

plot(cfy, rfy)
abline(a = 0, b = 1)

plot(cfyOOB, rfyOOB)
abline(a = 0, b = 1)

tr <- randomForest::getTree(rf, 10)
table(tr[, "status"]) ### -1 is terminal
### number of terminal nodes is larger than in cforest


