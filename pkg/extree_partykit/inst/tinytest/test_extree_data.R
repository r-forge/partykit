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
