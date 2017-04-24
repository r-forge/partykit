

library("mvtnorm")
library("party")
library("partykit")

set.seed(29)

#settings
N <- 253 
p <- 3
ntree<-50 #is higher in the original analysis but effect stays the same

#generate noise data
x <- rmvnorm(n=N,mean=rep(0,p),diag(rep(1,p)))
y <- rnorm(n=N)#, mean = c(0, 4)[(x[,1] > 0) + 1])
myData <- data.frame(y,x)

layout(matrix(1:4, nrow = 2))

partykit <- c()
partykit$model <- partykit::cforest(y ~ ., data=myData,trace=TRUE,ntree=ntree, mtry = 3,
    control = partykit::ctree_control(teststat = "quad",  
                            testtype = "Univ", mincriterion = 0, minsplit = 10,minbucket = 10))
partykit$y_hat <- predict(partykit$model, OOB=TRUE)
plot(y,partykit$y_hat) #does not look reansonable
c2 <- mean((y-partykit$y_hat)^2)
title(paste0('partykit, MSE=',round(c2, digits = 2)))

#maybe OOB argument is just ignored?
partykit$NoOOB <- predict(partykit$model, OOB=FALSE)
plot(y,partykit$NoOOB) 
c3 <- mean((y-partykit$NoOOB)^2)
title(paste0('partykit, OOB=FALSE, MSE=',round(c3, digits = 2))) 

w <- do.call("cbind", partykit$model$weights)


mycontrols <- party::cforest_unbiased(trace=TRUE,ntree=ntree,minsplit = 10,minbucket =
10, mtry = 3)
party <- c()
party$model <- party::cforest(y ~ .,controls=mycontrols,data=myData, weights = w)
pT <- do.call("cbind", party$model@prediction_weights(OOB = TRUE))



party$y_hat <- predict(party$model, OOB=TRUE)
plot(y,party$y_hat) #looks reasonable to me. high bias towards 0, no correlation 
c1 <- mean((y - party$y_hat)^2)
title(paste0('party, MSE=',round(c1, digits = 2)))
party$NoOOB <- predict(party$model, OOB = FALSE)
plot(y, party$NoOOB)
c1 <- mean((y - party$NoOOB)^2)
title(paste0('party OOB = FALSE, MSE=',round(c1, digits = 2)))

party::prettytree(party$model@ensemble[[10]])

partykit$model$nodes[[10]]

pT <- do.call("cbind", party$model@prediction_weights(OOB = TRUE))
pF <- do.call("cbind", party$model@prediction_weights(OOB = FALSE))

pkT <- predict(partykit$model, OOB = TRUE, type = "weights")
pkF <- predict(partykit$model, OOB = FALSE, type = "weights")

max(abs(pF - pkF))

max(abs(pT - pkT))

max(abs(party$y_hat - partykit$y_hat))
max(abs(party$NoOOB - partykit$NoOOB))

