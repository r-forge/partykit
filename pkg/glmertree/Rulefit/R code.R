platform <- "windows"
rfhome <- getwd()
source("rulefit.r")
library(akima) 

library(MASS)
names(Boston)
rfbosthouse <- rulefit(x=Boston, y="medv", cat.vars="chas", sparse=3, model.type='both')
rulefit(x=Boston, y="medv", cat.vars="chas", sparse=3, model.type='both', test.reps = 0)
rulefit(x=Boston, y="medv", cat.vars="chas", sparse=3, model.type='both', test.reps = 0, test.frac=.2)
rulefit(x=Boston, y="medv", cat.vars="chas", sparse=3, model.type='both', test.reps = 0, test.frac=.5)
rulefit(x=Boston, y="medv", cat.vars="chas", sparse=3, model.type='both', test.reps = 1)
rulefit(x=Boston, y="medv", cat.vars="chas", sparse=3, model.type='both', test.reps = 1, test.frac=.2)
rulefit(x=Boston, y="medv", cat.vars="chas", sparse=3, model.type='both', test.reps = 1, test.frac=.5)
rulefit(x=Boston, y="medv", cat.vars="chas", sparse=3, model.type='both', test.reps = 2)
rulefit(x=Boston, y="medv", cat.vars="chas", sparse=3, model.type='both', test.reps = 2, test.frac=.2)
rulefit(x=Boston, y="medv", cat.vars="chas", sparse=3, model.type='both', test.reps = 2, test.frac=.5)
rulefit(x=Boston, y="medv", cat.vars="chas", sparse=3, model.type='both', test.reps = 10)
rulefit(x=Boston, y="medv", cat.vars="chas", sparse=3, model.type='both', test.reps = 10, test.frac=.2)
rulefit(x=Boston, y="medv", cat.vars="chas", sparse=3, model.type='both', test.reps = 10, test.frac=.5)




library(glmnet)
glmnet(x = as.matrix(Boston[,-14]), y = Boston$medv)

library(lars)
lars(x = as.matrix(Boston[,-14]), y = Boston$medv, type = "forward.stagewise", normalize = FALSE)
