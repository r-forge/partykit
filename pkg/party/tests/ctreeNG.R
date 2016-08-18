
library("partyNG")

data("iris")

ctree(Species ~ ., data = iris, control = ctree_control(nmax = Inf))

partykit::ctree(Species ~ ., data = iris)

ctree(Species ~ ., data = iris, control = ctree_control(nmax = 50))

