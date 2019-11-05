library("partykit")

airq <- subset(airquality, !is.na(Ozone))
my_data <- extree_data(Ozone ~ Wind + Temp, 
    data = airq, yx = "matrix")

## This currently needs to be done within extree
trafo <- function(subset, data, weights, info = NULL, estfun = TRUE, object = TRUE) {
    estfun <- matrix(0, ncol = NCOL(data$yx$y), nrow = nrow(data$data))
    estfun[subset,] <- as.matrix(data$yx$y)[subset, ]
    list(estfun = estfun, objfun = -sum((data$yx$y - mean(data$yx$y))^2), converged = TRUE)
}

my_split <- function(...) {
    msfn <- function(model, trafo, data, subset, weights, whichvar, ctrl) {
        # args <- list(...)
        
        if (length(whichvar) == 0) return(NULL)
        
        ## split LAST variable at median
        for (j in whichvar) {
            x <- model.frame(data)[[j]][subset]
            ret <- partysplit(as.integer(j), breaks = median(x))
        }
        return(ret)
    }
    
    return(msfn)
}


# debug(my_split)
# undebug(partykit:::.extree_node)

tr <- extree(data = my_data, trafo = trafo, 
    control = c(extree_control(criterion = "statistic",
        logmincriterion = log(1 - 0.04),
        update = TRUE,
        selectfun = partykit:::.objfun_select(),
        splitfun = my_split(),
        svselectfun = partykit:::.objfun_select(),
        svsplitfun = my_split()),
        restart = TRUE))
