library("model4you")
library("survival")

### survreg
set.seed(1)
data(GBSG2, package = "TH.data")

## base model
bmod <- survreg(Surv(time, cens) ~ horTh, data = GBSG2, model = TRUE)
survreg_plot(bmod)

## partitioned model
tr <- pmtree(bmod)
plot(tr, terminal_panel = node_pmterminal(tr, plotfun = survreg_plot, 
                                          confint = TRUE))
summary(tr)
summary(tr, node = 1:2)

logLik(bmod)
logLik(tr)


### glm binomial 
set.seed(2)
n <- 1000
trt <- factor(rep(1:2, each = n/2))
age <- sample(40:60, size = n, replace = TRUE)
eff <- -1 + I(trt == 2) + 1 * I(trt == 2) * I(age > 50)
expit <- function(x) 1/(1 + exp(-x))

success <- rbinom(n = n, size = 1, prob = expit(eff))

dat <- data.frame(success, trt, age)
library("plyr")
dattab <- ddply(.data = dat, .variables = .(trt, age), 
                function(x) data.frame(nsuccess = sum(x$success),
                                       nfailure = NROW(x) - sum(x$success)))

bmod1 <- glm(success ~ trt, family = binomial)
bmod2 <- glm(success ~ trt, family = "binomial")
bmod3 <- glm(success ~ trt, data = dat, family = binomial)
bmod4 <- glm(cbind(nsuccess, nfailure) ~ trt, data = dattab, family = binomial)

(tr1 <- pmtree(bmod1, data = dat))
(tr2 <- pmtree(bmod2, data = dat))
(tr3 <- pmtree(bmod3))
(tr4 <- pmtree(bmod4))

(mtr1 <- glmtree(success ~ trt | age, data = dat, family = binomial))
# (mtr2 <- glmtree(cbind(nsuccess, nfailure) ~ trt | age, data = dattab, family = binomial))

library("strucchange")
sctest(tr3)
sctest(tr4)

## check logLik
logLik(mtr1)
logLik(tr1)

sum(objfun(tr1, newdata = dat))
objfun(tr1, newdata = dat, sum = TRUE)
sum(objfun(tr1))
objfun(tr1, sum = TRUE)


### linear model and missings
data("MathExam14W", package = "psychotools")

## scale points achieved to [0, 100] percent
MathExam14W$tests <- 100 * MathExam14W$tests/26
MathExam14W$pcorrect <- 100 * MathExam14W$nsolved/13

## select variables to be used
MathExam <- MathExam14W[ , c("pcorrect", "group", "tests", "study",
                             "attempt", "semester", "gender")]

## compute base model
bmod_math <- lm(pcorrect ~ group, data = MathExam)

## compute tree
(tr_math <- pmtree(bmod_math, control = ctree_control(maxdepth = 2)))

## create data with NAs
Math_mx <- Math_mz <- MathExam
Math_mx$group[1:2] <- NA
Math_mz$tests[1:20] <- NA

bmod_math_mx <- lm(pcorrect ~ group, data = Math_mx)

(tr_math_mx1 <- pmtree(bmod_math, control = ctree_control(maxdepth = 2), data = Math_mx))
(tr_math_mx2 <- pmtree(bmod_math_mx, control = ctree_control(maxdepth = 2)))
(tr_math_mz <- pmtree(bmod_math, control = ctree_control(maxdepth = 2), data = Math_mz))


## check logLik
(tr_math_mob <- lmtree(pcorrect ~  group | ., data = MathExam, maxdepth = 2))

logLik(bmod_math)
logLik(tr_math)
logLik(tr_math_mob)

sum(bmod_math$residuals^2)