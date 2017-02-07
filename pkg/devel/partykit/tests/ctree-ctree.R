
### compare partykit::ctree and party::ctree wrt predictions
### and computing time. But not during R CMD check

library("partykit")
library("party")

set.seed(29)

n <- 200
dgp <- function(n, response = c("numeric", "ordered", "binary", "factor"), na = FALSE) {

    response <- match.arg(response)
    x <- runif(n * 10)
    x <- matrix(x, ncol = 10)
    x <- cbind(x, x[,1:2] + cbind(rnorm(n, sd = 0.1), rnorm(n, sd = 0.1)))
    x1 <- sample(gl(4, n / 4))
    x2 <- sample(ordered(gl(4, n / 4)))
    y <- as.integer(interaction(x[,1] > 0.5, x[,2] > 0.25))
    y <- rnorm(n, mean = y)
    if (na) x[sample(1:(n * 10), n)] <- NA
    if (response == "numeric")
        return(data.frame(x, x1, x2, y = y))
    if (response == "binary") {
        y <- factor(y < median(y))
        return(data.frame(x, x1, x2, y = y))
    } else {
        y <- cut(y, 
            breaks = c(-Inf, quantile(y, c(0.25, 0.5, 0.75)), Inf))
    }
    if (response == "ordered")
        return(data.frame(x, x1, x2, y = ordered(y)))
    if (response == "factor")
        return(data.frame(x, x1, x2, y = y))
}

args <- list(response = c("numeric", "ordered", "binary", "factor"),
             na = c(FALSE, TRUE),
             teststat = c("quad", "max"),
             testtype = c("Bonferroni", "Univariate", "Teststatistic"),
             mincriterion = c(0.8, 0.9, 0.95),
             maxsurrogate = 0:2, maxdepth = 3)

gr <- do.call("expand.grid", args)
tgr <- gr
gr[["response"]] <- as.character(gr[["response"]])
gr[["teststat"]] <- as.character(gr[["teststat"]])
gr[["testtype"]] <- as.character(gr[["testtype"]])

### random splits are done differently in party::ctree
gr <- subset(gr, !na | (na & maxsurrogate > 0))


fun <- function(args) {

    learn <- dgp(n, response = args$response, na = args$na)
    test <- dgp(n, response = args$response, na = args$na)

    ct <- party:::ctree_control
    ctrl <- do.call("ct", args[-(1:2)])
    ot <- system.time(oldmod <- party:::ctree(y ~ ., data = learn, control = ctrl))[1]
    old <- predict(oldmod, newdata = test, type = "node")
    ct <- partykit:::ctree_control
    ctrl <- do.call("ct", args[-(1:2)])
    ctrl$majority <- TRUE
    ctrl$nmax <- Inf
    ctrl$splitstat <- "maximum"
    ctrl$numsurrogate <- TRUE
    nt <- system.time(newmod <- try(partykit:::ctree(y ~ ., data = learn, control = ctrl)))[1]
    if (inherits(newmod, "try-error")) {
        new <- as.list(rep(1, length(old)))
    } else {
        new <- predict(newmod, newdata = test, type = "node")
    }
    tb <- table(old, new) > 0
    err <- !(unique(rowSums(tb)) == 1 && unique(colSums(tb)) == 1)
    print(err)
    cat("---\n")
    list(error = err, time = c(ot, nt))
}

ret <- lapply(1:nrow(gr), function(i) { print(i); fun(gr[i, ,drop = FALSE]); })

tm <- t(sapply(ret, function(x) x$time))
err <- sapply(ret, function(x) x$error)

save(ret, tm, tgr, fun, gr, dgp, err, file = "results-regtest.Rda")

gr <- cbind(gr, err = err)
subset(gr, !na & err > 0.01)

### most of the differences come from the fact that in the absence of
### a surrogate splits observations are randomly assigned to daugther nodes
### slightly differently in the two implementations
