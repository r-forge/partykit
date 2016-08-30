library(partykit)
library(glmnet)
library(MatrixModels)
source("C:/Users/tobii/Desktop/swReg/upre/R/upre.R")

## test on air quality data:
airq <- airquality[complete.cases(airquality),]
dim(airq)
# some tests to see what ctree does with factors:
plot(ctree(Ozone ~ ., data=airq))
airq.f <- airq
airq.f$Ozone <- factor(airq$Ozone > median(airq$Ozone))
plot(ctree(Ozone ~ ., data=airq.f))
airq.f2 <- airq
airq.f2$Ozone <- cut(airq$Ozone, breaks = 4)
plot(ctree(Ozone ~ ., data=airq.f2))

plot(ctree(Ozone ~ offset(Temp) + Temp, data=airq.f))
plot(ctree(Ozone ~ Temp, data=airq.f))
tmp <- airq$Temp
# ctree does not require a family. Automatically uses continous or (multi)nominal methods. 
# if learning rate > 0, ctree cannot be used for factors, as it does not take offset argument, 
# but glmtree does. Note: this gives a different tree than ctree!
plot(glmtree(Ozone ~ 1 | Temp, data=airq.f, family = "binomial", offset = tmp))

predict(glmtree(Ozone ~ 1 | Temp, data=airq.f, family = "binomial"), newdata = airq.f)





airq.ens <- upre(Ozone ~ ., data=airq)
airq.ens.f <- upre(Ozone ~ ., data=airq.f, learnrate = 0)
airq.ens.f2 <- upre(Ozone ~ ., data=airq.f, learnrate = .000001)
predict(airq.ens.f)
predict(airq.ens.f, type = "response")
predict(airq.ens.f, type = "class")
predict(airq.ens.f2)
predict(airq.ens.f2, type = "response")
predict(airq.ens.f2, type = "class")

airq.coefs <- coef(airq.ens)
airq.coefs[airq.coefs$coefficient>0,]
importance(airq.ens)
plot(airq.ens$glmnet.fit)

airq.f.coefs <- coef(airq.ens.f)
airq.f.coefs[airq.f.coefs$coefficient!=0,]
importance(airq.ens.f)
plot(airq.ens.f$glmnet.fit)

singleplot(airq.ens, "Temp")
airq.nullmods <- bsnullinteract(airq.ens, nsamp = 5, seed = 42)
interact(airq.ens, c("Temp","Wind"), nullmods = airq.nullmods)
singleplot(airq.ens.f, "Temp")
singleplot(airq.ens.f2, "Temp")



## test on psychological data:
library(foreign)
cardata <- read.spss("data Carillo et al.sav", to.data.frame = TRUE)
summary(cardata)
car.ens <- upre(bdi ~ ., data=cardata[1:100,], lambda = seq(0, 10, by = .1))
plot(car.ens$glmnet.fit)
coefs <- coef(car.ens)
# something going wrong here:
preds <- predict(car.ens, newdata = cardata[101:111,])
cor(preds, cardata$bdi[101:111])
cor(preds, cardata$bdi[101:111])^2 * var(cardata$bdi[101:111])
coefs
singleplot(car.ens, "ntot")
pairplot(car.ens, c("n3","ntot"))
importance(car.ens)

# convert a var to factor and see what happens:
cardata.fac <- cardata
cardata.fac[,"open4"] <- factor(cut(cardata.fac[,"open4"], breaks = 8, labels = FALSE))
car.ens.fac <- upre(bdi ~ ., data=cardata.fac[1:100,], ctreecontrol = ctree_control(alpha = .10))
plot(car.ens.fac$glmnet.fit)
coefs <- coef(car.ens.fac)
predict(car.ens.fac, newdata = cardata.fac[101:111,])
singleplot(car.ens.fac, "open4")
pairplot(car.ens.fac, c("n1","open4"))
importance(car.ens.fac)
importance(car.ens.fac, global = FALSE)


## Test on simulated data:
source("simulated data.R")
# Derive rule ensemble:
sim.ens <- upre(y ~ ., data=simdata2[1:400,])
plot(sim.ens$glmnet.fit)
# Get coefficients:
coefs <- coef.upre(sim.ens, penalty.par.val = "lambda.min")
# Get predictions and assess accuracy:
preds <- predict.upre(sim.ens, newdata = simdata[401:500,-11], penalty.par.val = "lambda.min")
preds <- predict.upre(sim.ens, newdata = simdata[401:500,], penalty.par.val = "lambda.min")
cor(preds, simdata$y[401:500])
varnames <- paste("x", 1:10, sep="")
for(i in 1:10) {
  singleplot(sim.ens, varnames[i], nvals = 10)
}
# X1 and 2, x3 and 4, x5 and 6 should show interaction:
pairplot(sim.ens, c("x1","x2"), nvals = c(5,5), phi = 45, theta = 315)
pairplot(sim.ens, c("x3","x4"), nvals = c(10,10))
pairplot(sim.ens, c("x5","x6"), nvals = c(10,10))
# Calculate importances:
finalimps <- importance(sim.ens)
sum(finalimps$varimps$imp)
sum(finalimps$baseimps$imp)
# calculate local importances:
importance(sim.ens, global = FALSE)
# Calculate interactions:
nullmods <- bsnullinteract(sim.ens, nsamp = 2)
int.bs <- interact(sim.ens, paste("x", 1:10, sep=""), nullmods)
int.obs <- interact(sim.ens, paste("x", 1:10, sep=""))




# To implement:
#
# Non-continuous output variables:
# ctree also allows for:
# censored (see examples in ?ctree, requires survival package) 
# ordered (single ordered factor output var), 
# nominal (single unordered factor output var) and 
# multivariate (multiple output vars, specified with + in left hand of formula statement) 
# for non-negative counts, could use glmtree() and specify familiy = "poisson"

# glmnet allows for: 
# Quantitative for family="gaussian", or family="poisson" (non-negative counts)
# For family="binomial" should be either a factor with two levels, or a two-column matrix of 
# counts or proportions (the second column is treated as the target class; for a factor, 
# the last level in alphabetical order is the target class). For family="multinomial", can be 
# a nc>=2 level factor, or a matrix with nc columns of counts or proportions. For either "binomial"
# or "multinomial", if y is presented as a vector, it will be coerced into a factor. For family="cox", 
# y should be a two-column matrix with columns named 'time' and 'status'. The latter is a binary 
# variable, with '1' indicating death, and '0' indicating right censored. The function Surv() in 
# package survival produces such a matrix. For family="mgaussian", y is a matrix of quantitative 
# responses.
# censored:
data("GBSG2", package = "TH.data")
library("survival")
(stree <- ctree(Surv(time, cens) ~ ., data = GBSG2))
GBSG2$time
GBSG2$cens






## Compare with rulefit:
sim.ens <- upre(y ~ ., data=simdata2[1:400,])
uprepreds <- predict(sim.ens, newdata = simdata2[401:500,])
importance(sim.ens)

platform <- "windows"
rfhome <- "C:/Users/tobii/Desktop/swReg/Rulefit"
setwd(rfhome)
source("rulefit.r")
library(akima) 
rfmod <- rulefit(simdata2[1:400,paste("x", 1:10, sep = "")], simdata2[1:400,"y"])
varimps <- varimp(rfmod, plot = F)
varimps
singleplot("x1")
pairplot("x3", "x4", type = "persp")
rulefitpreds <- rfpred(simdata2[401:500,paste("x", 1:10, sep = "")])

cor(uprepreds, simdata2[401:500, "y"])
cor(rulefitpreds, simdata2[401:500, "y"])
cor(rulefitpreds, uprepreds)

null.models <- intnull(ntimes=10, quiet=F)
int <- interact(paste("x", 1:10, sep = ""), null.mods = null.models)
int








##########
## LARS ##
##########

library(lars)
bdi_ind <- which(names(cardata) == "bdi")
lars(as.matrix(cbind(cardata[1:100,-bdi_ind],  car.init.ens$rulevars)), cardata[1:100,"bdi"], type = "forward.stagewise", use.Gram = F, trace = T)