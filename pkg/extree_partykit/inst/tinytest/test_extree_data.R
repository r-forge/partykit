# -------------------------------------------------------------------
# - NAME:   test_extree_data.R 
# - AUTHOR: Heidi Seibold
# - DATE:   2021-01-08
# -------------------------------------------------------------------
# - PURPOSE: Test extree_data(): Does it contain the expected entries.
# -------------------------------------------------------------------



# -------------------------------------------------------------------
# AIRQUALITY EXAMPLE
# -------------------------------------------------------------------

airq <- subset(airquality, !is.na(Ozone))

aq1 <- extree_data(formula = Ozone ~ Wind + Temp,
  data = airq, yx = "matrix")

aq2 <- extree_data(formula = list(y = "Ozone", z = c("Wind", "Temp")),
  data = airq, yx = "matrix")

aq3 <- extree_data(formula = Ozone ~ Wind + Temp,
  data = airq)

aq4 <- extree_data(formula = list(y = "Ozone", z = c("Wind", "Temp")),
  data = airq, yx = "matrix")

## Check variable names of z
expect_equal(attributes(aq1$variables$z), attributes(aq2$variables$z))
expect_equal(attributes(aq1$variables$z), attributes(aq3$variables$z))
expect_equal(attributes(aq1$variables$z), attributes(aq4$variables$z))

## Check terms
expect_equal(aq1$terms, aq3$terms)
expect_null(aq2$terms)
expect_null(aq4$terms)

## Check zindex
#expect_equal(aq1$zindex, aq2$zindex) FIXME: (ML) More variables in aq2 than in aq1


# -------------------------------------------------------------------
# Pima Indians diabetes EXAMPLE
# -------------------------------------------------------------------
data("PimaIndiansDiabetes", package = "mlbench")

# with formula
ed1 <- extree_data(diabetes ~ glucose |  mass + pedigree + age,
  data = PimaIndiansDiabetes)

# without formula (useful for big data sets)
data_list <- list(
  y = c("diabetes"),
  x = c("glucose"),
  z = c("mass", "pedigree", "age")
)
ed2 <- extree_data(formula = data_list, data = PimaIndiansDiabetes)

## Check print
expect_identical(print(ed1), print(ed2))




# -------------------------------------------------------------------
# Working with missings EXAMPLE
# -------------------------------------------------------------------

# generate example data
set.seed(2020)
n <- 15
X <- rnorm(n = n)
Z <- rnorm(n = n)
Y <- rnorm(n = n, mean = X + I(Z > 0))
Zmiss <- Z
Zmiss[3:5] <- NA
O <- sample(1:2, size = n, replace = TRUE)
Omiss <- O
Omiss[1:2] <- NA
dat <- data.frame(Y, X, Z, Zmiss, O, Omiss)

exd <- extree_data(formula = list(y = c("Y"), x = c("X"), z = c("Zmiss"), 
  offset = c("Omiss")), data = dat)


## check missings in variable Zmiss
expect_equal(extree_variable(x = exd, index = "Zmiss", type = "missings"), c(3, 4, 5))
expect_equal(exd$missings$Zmiss, c(3, 4, 5))

## check missings in yx (this includes the offset, of course)
expect_equal(extree_variable(x = exd, variable = "yx", type = "missings"), c(1, 2)) # missings in offset Omiss
expect_equal(exd$missings$Y, integer())
expect_equal(exd$missings$X, integer())
expect_equal(exd$missings$`(offset)`, c(1, 2))

## missings should not be computed in unused variables
expect_equal(extree_variable(x = exd, index = "Z", type = "missings"), NA)
expect_equal(extree_variable(x = exd, index = "O", type = "missings"), NA)

## check that extree
exd$missings$Zmiss


# -------------------------------------------------------------------
# Tests for standard and rare inputs 
# -------------------------------------------------------------------
d <- data.frame(y = rep(1:5, each = 2), z = 1:10)
d$y[c(1, 10)] <- NA
w <- seq(0, 1, length.out = nrow(d))

# formula
exf1 <- extree_data(y ~ z, data = d)
exf2 <- extree_data(list(y = "y", z = "z"), data = d)
expect_equivalent(exf1$data, exf2$data)

## Subset
sub1 <- c(2, 7)
exs1 <- extree_data(y ~ z, data = d, subset = sub1)
expect_equivalent(exs1$data, d[sub1,])
# double entries allowed
sub2 <- c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0)
exs2 <- extree_data(y ~ z, data = d, subset = sub2)
expect_equivalent(exs2$data, d[sub2,])

## na.action
exna1 <- extree_data(y ~ z, data = d, na.action = na.exclude)
expect_false(any(is.na(exna1)))
exna2 <- extree_data(y ~ z, data = d, na.action = na.pass)
expect_equivalent(exna2$data, d)
expect_error(extree_data(y ~ z, data = d, na.action = na.fail))

## weights 
exw1 <- extree_data(y ~ z, data = d, weights = w)
expect_equal(exw1$data$`(weights)`, w)

expect_error(extree_data(y ~ z, data = d, weights = w[-1]), 
    pattern = "variable lengths differ")

expect_silent(extree_data(y ~ z, data = d, weights = sub2))

## offset
# Formula
# w part of d
exof1 <- extree_data(y ~ z, data = d, offset = w)
expect_equal(exof1$data$`(offset)`, w)
expect_equal(exof1$data[,exof1$variables$offset], w) 
d$w <- w
exof2 <- extree_data(y ~ offset(w) | z, data = d)
# w not part of d
d$w <- NULL
exof3 <- extree_data(y ~ offset(w) | z, data = d)
expect_equal(exof2$data$`(offset)`, w, exof3$data$`(offset)`) 
expect_equal(exof2$data[,exof2$variables$offset], w, exof3$data[,exof3$variables$offset])
# Multiple offsets 
exof4 <- extree_data(y ~ offset(w) + offset(log(w))| z, data = d)
exof5 <- extree_data(y ~ offset(log(w))| z, offset = w, data = d)
expect_equal(exof4$data[, exof4$variables$offset], exof5$data[, exof5$variables$offset],  w + log(w))
expect_error(extree_data(y ~ z, offset = cbind(w, log(w)), data = d), "unsupported specification of 'offset'")
# expect_error(extree_data(y ~ z, offset = cbind(w, log(w)), data = d), "number of offsets is 20, should equal 10")
exof6 <- extree_data(y ~ offset(log(w)) + offset(-w)| z, offset = w, data = d)
expect_equal(exof6$data[, exof6$variables$offset],  log(w))
# FIXME: (SD) should matrix of offsets work as well for multiple y's? --> lm()
# d$w <- NULL
# w2 <- log(w)
# exof7 <- extree_data(cbind(y, log(y)) ~ z, offset = as.matrix(cbind(w, w2)), data = d)
# exof8 <- extree_data(y + log(y) ~ z, offset = as.matrix(cbind(w, w2)), data = d)
## List instead of formula
d$w <- w
d$logw <- log(w)
expect_error(extree_data(list(y = "y", z = "z", offset = w), data = d))
exol1 <- extree_data(list(y = "y", z = "z"), offset = "w", data = d)
expect_equal(exol1$data$w, w) ## FIXME: (SD) should data differ to data above?
expect_equal( exol1$data[, exol1$variables$offset], w)
d$w <- NULL
exol2 <- extree_data(list(y = "y", z = "z"), offset = w, data = d)
expect_equal(exol2$data[, exol2$variables$offset], w)
# multiple offsets: 
d$w <- w
exol3 <- extree_data(list(y = "y", z = "z", offset = c("w", "logw")), offset = w, data = d)
expect_equal(exol3$data[, exol3$variables$offset], 2*w+log(w))

## cluster
clc <- sample(c(1, 2, 3), size = nrow(d), replace = TRUE)
cl <- as.factor(clc)
exc1 <- extree_data(y ~ z, data = d, cluster = cl)
expect_equivalent(exc1$data[,exc1$variables$cluster], cl)
exc2 <- extree_data(y ~ z, data = d, cluster = clc)
exc3 <- extree_data(list(y = "y", z = "z"), cluster = clc, data = d)
expect_error(extree_data(list(y = "y", z = "z", cluster = cl), data = d), "need to be of type character") #(SD) Do we want that factors are not allowed for cluster in formula?
expect_error(extree_data(list(y = "y", z = "z", cluster = "cl"), data = d), "not found in 'data'")
d$clc <- clc
d$cl <- cl
exc2 <- extree_data(list(y = "y", z = "z", cluster = "clc"), data = d)
exc3 <- extree_data(list(y = "y", z = "z", cluster = "cl"), data = d)
expect_equal(exc2$data, exc3$data) 

## scores
# sc <-  seq(0.1, 1, length.out = 10)
# d$zf <- ordered(d$z)
# exsc1 <- extree_data(y ~ zf, data = d)
# expect_error(extree_data(y ~ zf, data = d, scores = sc), "unsupported specification") 
# exsc3 <- extree_data(y ~ zf, data = d, scores = list(zf = sc))
# expect_equal(exsc3$scores$zf, sc)
# exsc4 <- extree_data(zf ~ y, data = d, scores = list(zf = sc))
# # FIXME (SD): Expect error if nr of scores != nr of levels? 
# expect_error(extree_data(y ~ zf, data = d, scores = list(zm = sc)), 
#   "names of 'scores' must match names")
# expect_error(extree_data(y ~ zf, data = d, scores = list(zf = sc[c(1:8)])), 
#   "number of scores in 'scores' must match number of levels of ordered factors")

## yx 
exyx1 <- extree_data(y ~ z, data = d, yx = "matrix")
expect_equal(exyx1$data$y, exyx1$yx$y)
expect_inherits(exyx1$yx$y, "integer")

## ytype 
exyt1 <- extree_data(y ~ z, data = d, yx = "matrix", ytype = "data.frame")
expect_inherits(exyt1$yx$y, "data.frame")
exyt2 <- extree_data(y ~ z, data = d, yx = "matrix", ytype = "matrix")
expect_inherits(exyt2$yx$y, "matrix")
exyt3 <- extree_data(y ~ z, data = d, yx = "none", ytype = "data.frame")
expect_equal(exyt3, exf1) # ytype is ignored if yx = "none" 

## nmax
# # convert to numeric since no binning for integer vars
dn <- data.frame(apply(d, MARGIN = 2, as.numeric))
dn$x <- as.numeric(sample(1:20, size = 10))
exn1 <- extree_data(y ~ z, data = dn, yx = "matrix", ytype = "matrix",
  nmax = c("yx" = 3, "z" = 3))
expect_equivalent(exn1$yxindex, inum::inum(dn[, "y", drop = FALSE], total = TRUE, 
  complete.cases.only = TRUE, nmax = 3, as.interval = "y"))
exn1$zindex #FIXME: (SD) why is there an element y?  
exn2 <- extree_data(y + x ~ z, data = dn, yx = "matrix", ytype = "matrix",
  nmax = c("yx" = 3, "z" = 3))
exn3 <- extree_data(y + x ~ z, data = dn, yx = "matrix", ytype = "vector", 
  nmax = c("yx" = 3, "z" = 3))
exn4 <- extree_data(y + x ~ z, data = dn, nmax = c("yx" = 3))
expect_true(length(levels(exn4$zindex$z)) != 3)
expect_error(extree_data(y ~ z, data = dn, nmax = c("heo" = 20)), "names of 'nmax'")
exn5 <- extree_data(y + x ~ z, data = dn, nmax = c("z" = 3))
expect_true(is.null(exn5$yxindex))
expect_true(length(levels(exn5$zindex$z))== 3)

