
library("ATR")
library("partykit")

ctrl <- ctree_control(maxsurrogate = 3)

### AirQuality: Standard example
data("airquality", package = "datasets")
tree <- rotate(ctree(Wind ~ . , data = airquality, 
                                 control = ctrl))
#pdf("tmp.pdf", width = 12, height = 10)
plot(tree, main = "Standard example")
#dev.off()

### C02: example with factors
data("CO2")
tree <- rotate(ctree(uptake ~ . , data = as.data.frame(CO2), control = ctrl))
plot(tree, main = " Example with factors", tnex = 1)

### AirQuality: Example of a  very large tree
data("airquality", package = "datasets")
ctrl <- ctree_control(mincriterion = 1e-11, minsplit = 1, maxsurrogate = 3)
tree <- rotate(ctree(Wind ~ . , data = airquality, control = ctrl))
plot(tree, main = "Example of a large tree", tnex = 1.5, cex=.5, nobs.loc='top')

### AirQuality: Example of removed nobs
data("airquality", package = "datasets")
ctrl <- ctree_control(mincriterion = 1e-11, minsplit = 1, maxsurrogate = 3)
tree <- rotate(ctree(Wind ~ . , data = airquality, control = ctrl))
plot(tree, main = "Example with removed 'number of observations'", tnex = 1.5, cex=.5, tree.offset=-0.1, remove.nobs=TRUE)


### AirQuality: Example of a one-leave tree
data("airquality", package = "datasets")
ctrl <- ctree_control(mincriterion = 1, maxdepth = 1, maxsurrogate = 3)
tree <- rotate(ctree(Wind ~ . , data = airquality, control = ctrl))
plot(tree, main = "Example of a one-leave tree")

### AirQuality: keyword argument "debug" shows viewport-hitboxes
data("airquality", package = "datasets")
ctrl <- ctree_control(maxsurrogate = 3)
tree <- rotate(ctree(Wind ~ . , data = airquality, control = ctrl))
plot(tree, main = "Debug modus", tnex = 1.5, debug=TRUE)

# options for ctree_control:
# c <- ctree_control(mincriterion = 0.215, minsplit = 1, minbucket = 7, 
# stump = FALSE, nresample = 9999, maxsurrogate = 0, 
# mtry = 0, savesplitstats = TRUE, maxdepth = 0, remove_weights = FALSE)

