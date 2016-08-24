
library("partyNG")
library("sandwich")
library("libcoin")
data("PimaIndiansDiabetes", package = "mlbench")

glmtree <- function(formula, data, family = gaussian(), 
                    control = ctree_control(), ...) {

    if (control$nmax < Inf) {
        glmtrafo <- function(formula, data, weights, block, ctrl) {
            if (!is.null(block)) stop("block not implemented")   
            mf <- model.frame(formula, data, na.action = na.pass)
            bdr <- BDR::BDR(mf, complete.cases.only = TRUE, total = TRUE)
            mf2 <- as.data.frame(bdr)
            iy <- c(bdr)
            attr(iy, "levels") <- 1:nrow(mf2)
            mfs <- model.frame(formula, data = mf2)
            y <- model.response(mfs)
            x <- model.matrix(formula, data = mf2)
            function(subset) {
                w <- c(libcoin::ctabs(iy, weights = weights, subset = subset)[-1L])
                mod <- glm(y ~ x + 0, family = family, weights = w)
                Y <- estfun(mod)
                Y <- Y / w
                Y[w == 0,] <- 0
                ret <- rbind(0, Y)
                list(estfun = ret, index = iy, coef = coef(mod), logLik = logLik(mod))
            }
        }
    } else {
        glmtrafo <- function(formula, data, weights, block, ctrl) {
            if (!is.null(block)) stop("block not implemented")
            mf <- model.frame(formula, data, na.action = na.pass)
            cc <- complete.cases(mf)
            y <- model.response(mf)
            x <- model.matrix(formula, data = mf)
            function(subset) {
                s <- subset[cc[subset]]
                ys <- y[s]
                xs <- x[s, , drop = FALSE]
                if (length(weights) > 0) {
                    w <- weights[cc[subset]]
                    mod <- glm(ys ~ xs + 0, family = family, weights = w)
                } else {
                    mod <- glm(ys ~ xs + 0, family = family)
                }
                ret <- matrix(0, nrow = NROW(x), ncol = NCOL(x))
                ret[subset,] <- estfun(mod)
                storage.mode(ret) <- "double"
                list(estfun = ret, coef = coef(mod), logLik = logLik(mod))
            }
        }
    }
    ctree(formula, data, ytrafo = glmtrafo, control = control, ...)
}


     ## recursive partitioning of a logistic regression model
system.time(     pid_tree <- partykit::glmtree(diabetes ~ glucose | pregnant +
       pressure + triceps + insulin + mass + pedigree + age,
       data = PimaIndiansDiabetes, family = binomial))

pid_tree

nd <- factor(predict(pid_tree, type = "node"))
m <- glm(diabetes ~ glucose:nd, data = PimaIndiansDiabetes, family = binomial)
logLik(m)
AIC(m)


     ## recursive partitioning of a logistic regression model
system.time(     pid_tree1 <- glmtree(diabetes ~ glucose | pregnant +
       pressure + triceps + insulin + mass + pedigree + age,
       data = PimaIndiansDiabetes, family = binomial))

pid_tree1

info_node(node_party(pid_tree1))

nd1 <- factor(predict(pid_tree1, type = "node"))
m1 <- glm(diabetes ~ glucose:nd1, data = PimaIndiansDiabetes, family = binomial)
logLik(m1)
AIC(m1)


table(nd, nd1)



     ## recursive partitioning of a logistic regression model
system.time(     pid_tree2 <- glmtree(diabetes ~ glucose | pregnant +
       pressure + triceps + insulin + mass + pedigree + age,
       data = PimaIndiansDiabetes, family = binomial, 
       control = ctree_control(nmax = 20)))

pid_tree2

info_node(node_party(pid_tree2))

nd2 <- factor(predict(pid_tree2, type = "node"))
m2 <- glm(diabetes ~ glucose:nd2, data = PimaIndiansDiabetes, family = binomial)
logLik(m2)
AIC(m2)


table(nd, nd2)
table(nd1, nd2)


