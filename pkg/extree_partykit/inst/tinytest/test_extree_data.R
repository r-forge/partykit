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
expect_equal(extree_variable(x = exd, i = "Zmiss", type = "missings"), c(3, 4, 5))
expect_equal(exd$missings$Zmiss, c(3, 4, 5))

## check missings in yx (this includes the offset, of course)
expect_equal(extree_variable(x = exd, i = "yx", type = "missings"), c(1, 2)) # missings in offset Omiss
expect_equal(exd$missings$Y, integer())
expect_equal(exd$missings$X, integer())
expect_equal(exd$missings$Omiss, c(1, 2))

## missings should not be computed in unused variables
expect_equal(extree_variable(x = exd, i = "Z", type = "missings"), NA)
expect_equal(extree_variable(x = exd, i = "O", type = "missings"), NA)

## check that extree
exd$missings$Zmiss


# -------------------------------------------------------------------
# Tests for standard and rare inputs 
# -------------------------------------------------------------------
# formula
ex0 <- extree_data(y ~ z, data = d)
ex1 <- extree_data(list(y = "y", z = "z"), data = d)
expect_equivalent(ex0$data, ex1$data)

## Subset
sub1 <- c(2, 7)
ex1 <- extree_data(y ~ z, data = d, subset = sub1)
expect_equivalent(ex1$data, d[sub1,])
# double entries allowed
sub2 <- c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0)
ex2 <- extree_data(y ~ z, data = d, subset = sub2)
expect_equivalent(ex2$data, d[sub2,])

## na.action
ex3 <- extree_data(y ~ z, data = d, na.action = na.exclude)
expect_false(any(is.na(ex3)))

ex4 <- extree_data(y ~ z, data = d, na.action = na.pass)
expect_equivalent(ex4$data, d)

expect_error(extree_data(y ~ z, data = d, na.action = na.fail))

## weights 
w <- seq(0, 1, length.out = 10)
ex5 <- extree_data(y ~ z, data = d, weights = w)
expect_equal(ex5$data$`(weights)`, w)

expect_error(extree_data(y ~ z, data = d, weights = w[-1]), 
    pattern = "variable lengths differ")

expect_silent(extree_data(y ~ z, data = d, weights = sub2))

## offset
# Formula
ex6 <- extree_data(y ~ z, data = d, offset = w)
expect_equal(ex6$data$`(offset)`, w)
d$w <- w
ex7 <- extree_data(y ~ offset(w) | z, data = d)
d$w <- NULL
ex8 <- extree_data(y ~ offset(w) | z, data = d)
ex9 <- extree_data(y ~ z, offset = w, data = d)
# FIXME: (SD) following currently failes due to wrongly specified off <- attr(mt$x, "offset")
# expect_equal(ex7$data$`(offset)`, w, ex7$data$`(offset)`) 
# expect_equal(ex7$data[,ex7$variables$offset], w) 
# expect_equal(ex8$data[,ex8$variables$offset], w)
expect_equal(ex9$data[,ex9$variables$offset], w) 
expect_equivalent(ex7, ex8)
# List instead of formula
d$w <- w
ex8 <- extree_data(list(y = "y", z = "z", offset = "w"), data = d)
expect_equal(ex8$data$w, w) ## FIXME: (SD) should data differ to data above?
expect_equal( ex8$data[, ex8$variables$offset], w)
ex9 <- extree_data(list(y = "y", z = "z"), offset = "w", data = d)
d$w <- NULL
w1 <- w
ex10 <- extree_data(list(y = "y", z = "z"), offset = w1, data = d)
ex11 <- extree_data(list(y = "y", z = "z", offset = w1), data = d)
## FIXME: (SD) variables$offset gives wrong index due to vars[[v]] <- unique(as.integer(vars[[v]])) & 
## proper handling of numerics as offset --> add w1 to data & use correct index! 
# expect_equal(ex10$data[, ex10$variables$offset], w1) 
# expect_equal(ex11$data[, ex11$variables$offset], w1) 
expect_equal(ex10, ex11)

## cluster
cl <- factor(sample(c(1, 2, 3), size = nrow(d), replace = TRUE))
clc <- sample(c(1, 2, 3), size = nrow(d), replace = TRUE)
expect_silent(extree_data(y ~ z, data = d, cluster = cl))
ex7 <- extree_data(y ~ z, data = d, cluster = clc)
expect_warning(extree_data(list(y = "y", z = "z", cluster = cl), data = d), 
    pattern = "'cluster', must be character") #FIXME: (SD) Do we want this? 
ex8 <- extree_data(list(y = "y", z = "z", cluster = clc), data = d)
# expect_equal(ex7$data, ex8$data) #FIXME: (SD) Currently fails, ex8: cluster not in data + wrong index in variables


## scores
d$zf <- ordered(d$z)
ex8 <- extree_data(y ~ zf, data = d)
ex9 <- extree_data(y ~ zf, data = d, scores = seq(0.1, 1, length.out = 10)) 
expect_equal(ex8, ex9) #SD: Do we want this? 
ex10 <- extree_data(y ~ zf, data = d, scores = seq(1, 10, length.out = 8))
expect_equal(ex10, ex9)

## yx 
ex11 <- extree_data(y ~ z, data = d, yx = "matrix")
expect_equal(ex11$data$y, ex11$yx$y)
expect_inherits(ex11$yx$y, "integer")

## ytype 
ex12 <- extree_data(y ~ z, data = d, yx = "matrix", ytype = "data.frame")
expect_inherits(ex12$yx$y, "data.frame")
ex13 <- extree_data(y ~ z, data = d, yx = "matrix", ytype = "matrix")
expect_inherits(ex13$yx$y, "matrix")
ex14 <- extree_data(y ~ z, data = d, yx = "none", ytype = "data.frame")
expect_equal(ex14, ex0) # ytype is ignored if yx = "none" 

## nmax 
# ex <- extree_data(y ~ z, data = d, yx = "matrix", ytype = "matrix",
#     nmax = c("yx" = 3, "z" = 3))