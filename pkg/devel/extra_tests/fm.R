
library("partykit")
library("survival")

#set.seed(29)
n <- 100
d <- data.frame(y = runif(n), y1 = runif(n), e = TRUE, x = runif(n), z = runif(n), x1 = gl(2, n/2),
                b = sample(gl(2, n/2)), o2 = rep(2, n), o3 = rep(3, n))
d
w <- as.integer(round(runif(n) * 100))
w
ctree(y ~  z, data = d, weights = w, subset = x < .5)

ctree(y ~ x | z , data = d, weights = w)

ctree(y ~ 1 | z , data = d, weights = w)

ctree(y ~ 1 | z , data = d, cluster = b)

try(ctree(y ~ 1 | z, data = d, weights = w, cluster = b))

try(ctree(y ~ x | z, data = d, weights = w, cluster = b))

try(ctree(y ~ x + offset(o2 + o3) | z, data = d, weights = w, cluster = b))

ctree(Surv(y, e) ~ x:x1 | z , data = d, weights = w)

try(ctree(Surv(y, e) ~ x:x1 | z , data = d, weights = w, cluster = b))

try(ctree(y + y1 ~ x:x1 | z , data = d, weights = w , cluster = b))

try(ctree(y + y1 ~ . , data = d, weights = w))

try(ctree(y ~ . | z , data = d, weights = w, cluster = b))

