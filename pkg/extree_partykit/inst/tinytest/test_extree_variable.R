# -------------------------------------------------------------------
# - NAME:   test_extree_variable.R 
# - AUTHOR: Susanne Dandl
# - DATE:   2021-03-15
# -------------------------------------------------------------------
# - PURPOSE: Test extree_variable(): Does it contain the expected entries.
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# IRIS Example 
# -------------------------------------------------------------------
data("iris", package = "datasets")
# one y
ed <- extree_data(Species ~ Sepal.Width + Sepal.Length | Petal.Width + Petal.Length, 
  data = iris, nmax = c("yx" = 25, "z" = 10), yx = "none")
# multiple y
edm <- extree_data(Species + Sepal.Length ~ Sepal.Width | Petal.Width + Petal.Length, 
  data = iris, nmax = c("yx" = 25, "z" = 10), yx = "matrix")
# no binning
edb <- extree_data(Species ~ Sepal.Width + Sepal.Length | Petal.Width + Petal.Length, 
  data = iris, yx = "matrix") # Question: (SD) why is edb$zindex$Petal.Width & Pegal.Length specified?  
# only y - no x
ednox <- extree_data(Species ~ Petal.Width + Petal.Length, data = iris, 
  nmax = c("yx" = 25, "z" = 10))

# Original 
## numeric index
ev1 <- extree_variable(ed, 4, type = "original")
expect_equal(ev1, ed$data$Petal.Width)

ev2 <- extree_variable(ed, c(2, 4), type = "original")
expect_equal(ev2, as.list(ed$data[, c("Sepal.Width", "Petal.Width")])) 

ev3 <- extree_variable(ed, index = 1, type = "original") 

## character index
expect_null(extree_variable(ed, index = "y", type = "original")) # returns NULL because no variable named "y"
ec1 <- extree_variable(ed, index = "Sepal.Length", type = "original")
expect_equal(ec1, ed$data$Sepal.Length)
ec2 <- extree_variable(ed, index = c("Sepal.Length", "Sepal.Width"), type = "original")
expect_equal(ec2, as.list(ed$data[, c("Sepal.Length", "Sepal.Width")]))

## y and x variables
ey1 <- extree_variable(edb, variable = "y", type = "original")
expect_equal(ey1, iris$Species)
expect_equal(extree_variable(ed, variable = "y", type = "original"), ed$data$Species) 
exy <- extree_variable(ed, variable = "yx", type = "original")
expect_equal(data.frame(exy), ed$data[, c("Species", "Sepal.Width", "Sepal.Length")]) 
expect_error(extree_variable(ed, variable = "x"), "'arg' should be one of")
expect_equal(names(extree_variable(edm, variable = c("yx"))), c("Species", "Sepal.Length", "Sepal.Width"))
expect_equal(extree_variable(ed, variable = "y"), extree_variable(ed, i = 1)) 
# only y -  no x
expect_equal(ey1, extree_variable(ednox, variable = "y"))

# inum
expect_equal(attr(extree_variable(ed, variable = "y", type = "inum"), "levels"), 
  attr(ed$yxindex, "levels")[, 1, drop = FALSE])
ediy <- extree_variable(ed, variable = "y", type = "inum")
expect_equal(as.numeric(ediy), as.numeric(ed$yxindex))
expect_inherits(ediy, "inumtotal")
expect_equal(dim(attr(ediy, "levels")), c(nrow(attr(ed$yxindex, "levels")), 1))
extree_variable(ed, index = c(3, 4), type = "inum") # FIXME: (SD) how to handle this if binning? --> currently only z variables
extree_variable(edb, index = c(3, 4), type = "inum") 
expect_equivalent(class(extree_variable(edb, index = 5, type = "inum")), c("enum", "integer")) # FIXME: (SD) is this fine? 
expect_null(extree_variable(edb, variable = "yx", type = "inum"))
# only y no x
expect_equal(unique(extree_variable(ednox, variable = "y", type = "inum")), c(1,2,3))

# scores
# d <- data.frame(y = rep(1:5, each = 2), z = 1:10)
# d$y[c(1, 10)] <- NA
# sc <-  seq(0.1, 1, length.out = 10)
# d$zf <- ordered(d$z)
# exsc3 <- extree_data(y ~ zf, data = d, scores = list(zf = sc))
# expect_equal(extree_variable(exsc3, index = 2, type = "scores"), sc)
# expect_null(extree_variable(exsc3, index = 1, type = "scores"))

# missing
## data 
misid <- c(1, 10, 30)
misid2 <- misid + 5
iris2 <- iris
iris2$Sepal.Width[misid] <- NA 
iris2$Species[misid2] <- NA
edmis <- extree_data(Species ~ Sepal.Width + Sepal.Length | Petal.Width + Petal.Length, 
  data = iris2, nmax = c("yx" = 25, "z" = 10), yx = "matrix")
## tests
expect_equal(extree_variable(edmis, variable = "yx", type = "missings"), 
  sort(c(misid, misid2)))
expect_equal(extree_variable(edmis, variable = "y", type = "missings"), misid2)
em1 <- extree_variable(edmis, index = c("Sepal.Width", "Sepal.Length"), type = "missings")
expect_equal(em1$Sepal.Width, misid)
expect_equal(em1$Sepal.Length, integer(0))
## multiple y 
expect_equal(extree_variable(edm, variable = "y", type = "missings"), 
  list(Species = integer(0), Sepal.Length = integer(0)))

# index and variable 
expect_error(extree_variable(ed, index = 2, variable = "yx"), 
  "Specify either 'index' or 'variable' - not both")
expect_error(extree_variable(ed), "Specify either 'index' or 'variable'")
