library("model4you")
library("survival")

## function to simulate the data
sim_data <- function(n = 500, p = 10, beta = 3, sd = 1){
  
  ## treatment
  lev <- c("C", "A")
  a <- rep(factor(lev, labels = lev, levels = lev), length = n)
  
  ## correlated z variables
  sigma <- diag(p) 
  sigma[sigma == 0] <- 0.2
  ztemp <- mvtnorm::rmvnorm(n, sigma = sigma)
  z <- (pnorm(ztemp) * 2 * pi) - pi  
  colnames(z) <- paste0("z", 1:ncol(z))
  z1 <- z[,1]
  
  ## outcome
  y <- 7 + 0.2 * (a %in% "A") + beta * cos(z1) * (a %in% "A") + rnorm(n, 0, sd)
  
  data.frame(y = y, a = a, z)
}


## simulate data
set.seed(123)
beta <- 3
ntrain <- 500
ntest <- 100
simdata <- sim_data(p = 5, beta = beta, n = ntrain)
tsimdata <- sim_data(p = 5, beta = beta, n = ntest)
cens <- rep(1, ntrain)
tcens <- rep(1, ntest)

## base model
basemodel_lm <- lm(y ~ a, data = simdata)
basemodel_wb <- survreg(Surv(y, cens) ~ a, data = simdata)

## forest
frst_lm <- pmforest(basemodel_lm, ntree = 50, 
                    perturb = list(replace = FALSE, fraction = 0.632),
                    control = ctree_control(mincriterion = 0))
# vi_lm <- varimp(frst_lm)
frst_wb <- pmforest(basemodel_wb, ntree = 50, 
                    perturb = list(replace = FALSE, fraction = 0.632),
                    control = ctree_control(mincriterion = 0))

## personalised models
coefs_lm <- pmodel(x = frst_lm, newdata = tsimdata)
summary(coefs_lm)


coeffun <- function(model) {
  ## model coefficients
  coefs <- c(coef(model), scale = model$scale)
  
  ## difference in median survival 
  p = 0.5
  coefs["median_s0"] <- qweibull(p = p, shape = 1/coefs["scale"], 
                       scale = exp(coefs["(Intercept)"]))
  coefs["median_s1"] <- qweibull(p = p, shape = 1/coefs["scale"], 
                       scale = exp(coefs["(Intercept)"] + coefs["aA"]))
  coefs["median_sdiff"] <- coefs["median_s1"] - coefs["median_s0"]
  
  return(coefs)
}
coefs_wb <- pmodel(x = frst_wb, newdata = cbind(tsimdata, cens = tcens),
                   fun = coeffun)
summary(coefs_wb)



## dependence plot
# library("ggplot2")
# dp_lm <- cbind(coefs_lm, tsimdata)
# dp_wb <- cbind(coefs_wb, tsimdata)
# ggplot(tsimdata) + 
#   stat_function(fun = function(z1) 0.2 + beta * cos(z1), aes(color = "true treatment\neffect")) + 
#   geom_point(data = dp_lm, aes(y = aA, x = z1, color = "estimates lm"), alpha = 0.5) + 
#   geom_point(data = dp_wb, aes(y = median_sdiff, x = z1, color = "estimates wb"), alpha = 0.5)

