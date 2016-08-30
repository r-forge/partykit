
library("partyNG")

data("iris")

ctree(Species ~ ., data = iris, control = ctree_control(nmax = Inf))

ctree(Species ~ ., data = iris, control = ctree_control(nmax = Inf,
      splittest = TRUE, nresample = 1000, testtype = "MonteCarlo"))



partykit::ctree(Species ~ ., data = iris)

ctree(Species ~ ., data = iris, control = ctree_control(nmax = 50))

ctree(Species ~ ., data = iris, control = ctree_control(nmax = 50, splittest = TRUE,
                     nresample = 1000, testtype = "MonteCarlo"))

