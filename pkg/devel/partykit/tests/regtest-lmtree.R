library("partykit")

set.seed(29)
n <- 1000
x <- runif(n)
z <- runif(n)
y <- rnorm(n, mean = x * c(-1, 1)[(z > 0.7) + 1], sd = 3)
z_noise <- factor(sample(1:3, size = n, replace = TRUE))
d <- data.frame(y = y, x = x, z = z, z_noise = z_noise)


fmla <- as.formula("y ~ x | z + z_noise")

(m_mob <- mob2(formula = fmla, data = d, fit = partykit:::lmfit))
(m_lm2 <- lmtree2(formula = fmla, data = d))
(m_lm1 <- lmtree(formula = fmla, data = d))


mods <- nodeapply(m_lm2, ids = nodeids(m_lm2, terminal = TRUE), function(x) x$info$object)
sum(sapply(mods, function(x) sum(x$residuals^2)))


## pruning
lm_big <- lmtree2(formula = fmla, data = d, testflavour = "exhaustive", maxdepth = 2)
lm_aic <- prune.modelparty(lm_big, type = "AIC")
lm_bic <- prune.modelparty(lm_big, type = "BIC")

width(lm_big)
width(lm_aic)
width(lm_bic)

glm_big <- glmtree2(formula = fmla, data = d, testflavour = "exhaustive", maxdepth = 2)
glm_aic <- prune.modelparty(glm_big, type = "AIC")
glm_bic <- prune.modelparty(glm_big, type = "BIC")

width(glm_big)
width(glm_aic)
width(glm_bic)

get_test <- function(x, terminal = TRUE) 
  nodeapply(x, ids = nodeids(x, terminal = terminal), function(n) n$info$criterion["criterion", ])

sapply(get_test(lm_big), function(x) x)
sapply(get_test(glm_big), function(x) x)

sapply(nodeapply(lm_big, ids = nodeids(lm_big, terminal = TRUE), function(n) n$info$nobs), function(x) x)
sapply(nodeapply(glm_big, ids = nodeids(glm_big, terminal = TRUE), function(n) n$info$nobs), function(x) x)


lm_nodes <- nodeapply(lm_big, ids = nodeids(lm_big))
glm_nodes <- nodeapply(glm_big, ids = nodeids(glm_big))

unlist(sapply(lm_nodes, function(x) x$split$varid))
unlist(sapply(glm_nodes, function(x) x$split$varid))

lm_nodes[[3]]$info$coefficients
glm_nodes[[3]]$info$coefficients

lm_nodes[[3]]$info$nobs
glm_nodes[[3]]$info$nobs

lm_nodes[[2]]$split
glm_nodes[[2]]$split

all.equal(lm_nodes[[3]]$split, glm_nodes[[3]]$split)


glmtree(formula = fmla, data = d, alpha = 0.43)
lmtree(formula = fmla, data = d, alpha = 0.43)


