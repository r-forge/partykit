
library("profmem")

d <- iris
profmem(class(d) <- "list")

profmem(x <- d[["Sepal.Width"]])

profmem(y <- as.double(x))



