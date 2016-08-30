
### compare partyNG::ctree and party::ctree wrt predictions
### and computing time. But not during R CMD check

library("partyNG")
library("partykit")

set.seed(29)

n <- 200
dgp <- function(n, response = c("numeric", "ordered", "binary", "factor"), na = FALSE) {

    response <- match.arg(response)
    x <- round(runif(n * 10), 1)
    x <- matrix(x, ncol = 10)
    x <- cbind(x, x[,1:2] + cbind(rnorm(n, sd = 0.1), rnorm(n, sd = 0.1)))
    x1 <- sample(gl(4, n / 4))
    x2 <- sample(ordered(gl(4, n / 4)))
    y <- as.integer(interaction(x[,1] > 0.5, x[,2] > 0.25))
    y <- round(rnorm(n, mean = y), 1)
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

    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))

    ct <- partykit::ctree_control
    ctrl <- do.call("ct", args[-(1:2)])
    set.seed(29)
    ot <- system.time(oldmod <- partykit::ctree(y ~ ., data = learn, control = ctrl))[1]
    if (inherits(oldmod, "try-error")) {
        old <- as.list(rep(1, NROW(test)))
    } else {
        old <- predict(oldmod, newdata = test, simp = FALSE, 
            type = "node")
    }
    ct <- partyNG:::ctree_control
    ctrl <- do.call("ct", args[-(1:2)])
    ctrl$nmax <- Inf
    ctrl$splitstat <- "maximum"
    set.seed(29)
    nt <- system.time(newmod <- try(partyNG:::ctree(y ~ ., data = learn, control = ctrl)))[1]
    if (inherits(newmod, "try-error")) {
        new <- as.list(rep(1, nrow(test)))
    } else {
        new <- predict(newmod, newdata = test, simp = FALSE, 
                       type = "node")
    }
    ret <- list(error = max(abs(unlist(old) - unlist(new))), time = c(ot, nt))
    cat("Error: ", ret$error, "\n")
    ret
}

set.seed(29)

ret <- lapply(1:nrow(gr), function(i) { print(i); fun(gr[i, ,drop = FALSE]); })

tm <- t(sapply(ret, function(x) x$time))
err <- sapply(ret, function(x) x$error)

save(ret, tm, tgr, fun, gr, dgp, err, file = "results-regtest_new.Rda")

### scores
y <- gl(3, 10, ordered = TRUE)
x <- rnorm(length(y))
x <- ordered(cut(x, 3))
d <- data.frame(y = y, x = x)

### partyNG with scores
ct11 <- partyNG::ctree(y ~ x, data = d)
ct12 <- partyNG::ctree(y ~ x, data = d, 
                        scores = list(y = c(1, 4, 5)))
ct13 <- partyNG::ctree(y ~ x, data = d, 
                        scores = list(y = c(1, 4, 5), x = c(1, 5, 6)))

### party with scores
ct21 <- partykit::ctree(y ~ x, data = d)
ct22 <- partykit::ctree(y ~ x, data = d, 
                     scores = list(y = c(1, 4, 5)))
ct23 <- partykit::ctree(y ~ x, data = d, 
                     scores = list(y = c(1, 4, 5), x = c(1, 5, 6)))

all.equal(fitted(ct11), fitted(ct21))
all.equal(fitted(ct12), fitted(ct22))
all.equal(fitted(ct13), fitted(ct23))

### ytrafo
y <- runif(100, max = 3)
x <- rnorm(length(y))
d <- data.frame(y = y, x = x)

### partyNG with scores
ct11 <- partyNG::ctree(y ~ x, data = d)
ct12 <- partyNG::ctree(y ~ x, data = d,
                        ytrafo = list(y = sqrt))


### spotted by Peter Philip Stephensen (DREAM) <PSP@dreammodel.dk>
### splits x >= max(x) where possible in partyNG::ctree
library("partyNG")
nAge = 30
d = data.frame(Age=rep(1:nAge,2),y=c(rep(1,nAge),rep(0,nAge)),
n=rep(0,2*nAge))
ntot = 100
alpha = .5
d[d$y==1,]$n = floor(ntot * alpha * d[d$y==1,]$Age / nAge)
d[d$y==0,]$n = ntot - d[d$y==1,]$n
ctrl = partyNG::ctree_control(maxdepth=3, minbucket = min(d$n) + 1)
tree = partyNG::ctree(y ~ Age, weights=as.integer(n), data=d, control=ctrl)
tree

(w1 <- predict(tree, type = "node"))

ctrl = partykit::ctree_control(maxdepth=3, minbucket = min(d$n) + 1)
tree2 = partykit::ctree(y ~ Age, weights=d$n, data=d, control=ctrl)
tree2

(w2 <- predict(tree2, type = "node"))

stopifnot(max(abs(w1 - w2)) == 0)

gr <- cbind(gr, err = err)
subset(gr, !na & err > 0)

### most of the differences come from the fact that in the absence of
### a surrogate splits observations are randomly assigned to daugther nodes
### slightly differently in the two implementations
