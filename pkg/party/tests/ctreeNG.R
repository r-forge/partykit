
library("partyNG")

data("iris")

ctree(Species ~ ., data = iris, control = ctree_control(nmax = 100))