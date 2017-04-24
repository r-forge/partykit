
library("partykit")
library("sandwich")

set.seed(29)
n <- 1000
x <- runif(n)
z <- runif(n)
y <- rnorm(n, mean = x * c(-1, 1)[(z > 0.7) + 1], sd = 3)
z_noise <- factor(sample(1:3, size = n, replace = TRUE))
d <- data.frame(y = y, x = x, z = z, z_noise = z_noise)


fmla <- as.formula("y ~ x | z + z_noise")
fmly <- gaussian()
fit <- partykit:::glmfit

# versions of the data
d1 <- d
d1$z <- signif(d1$z, digits = 1)

k <- 20
zs_noise <- matrix(rnorm(n*k), nrow = n)
colnames(zs_noise) <- paste0("z_noise_", 1:k)
d2 <- cbind(d, zs_noise)
fmla2 <- as.formula(paste("y ~ x | z + z_noise +",
                          paste0("z_noise_", 1:k, collapse = " + ")))


d3 <- d2
d3$z <- factor(sample(1:3, size = n, replace = TRUE, prob = c(0.1, 0.5, 0.4)))
d3$y <- rnorm(n, mean = x * c(-1, 1)[(d3$z == 2) + 1], sd = 3)

## check weights
w <- rep(1, n)
w[1:10] <- 2
(mw1 <- glmtree(formula = fmla, data = d, weights = w))
(mw2 <- glmtree(formula = fmla, data = d, weights = w, caseweights = FALSE))



## check dfsplit
(mmfluc2 <- mob(formula = fmla, data = d, fit = partykit:::glmfit))
(mmfluc3 <- glmtree(formula = fmla, data = d))
(mmfluc3_dfsplit <- glmtree(formula = fmla, data = d, dfsplit = 10))


## check tests
library("strucchange")
sctest(mmfluc3, node = 1) # does not yet work

x <- mmfluc3
(tst3 <- nodeapply(x, ids = nodeids(x), function(n) n$info$criterion))




## check logLik and AIC
logLik(mmfluc2)
logLik(mmfluc3)
logLik(mmfluc3_dfsplit)
logLik(glm(y ~ x, data = d))

AIC(mmfluc3)
AIC(mmfluc3_dfsplit)



## check inner and terminal
options <- list(NULL,
                "object",
                "estfun",
                c("object", "estfun"))

arguments <- list("inner",
                  "terminal",
                  c("inner", "terminal"))


for (o in options) {
  print(o)
  x <- glmtree(formula = fmla, data = d, inner = o)
  str(nodeapply(x, ids = nodeids(x), function(n) n$info[c("object", "estfun")]), 2)
}

for (o in options) {
  print(o)
  x <- glmtree(formula = fmla, data = d, terminal = o)
  str(nodeapply(x, ids = nodeids(x), function(n) n$info[c("object", "estfun")]), 2)
}


## check model
m_mt <- glmtree(formula = fmla, data = d, model = TRUE)
m_mf <- glmtree(formula = fmla, data = d, model = FALSE)

dim(m_mt$data)
dim(m_mf$data)


## Check restart (restart = TRUE should take longer, but results should be the same)
t_rt <- system.time(m_rt <- glmtree(formula = fmla, data = d, restart = TRUE, testflavour = "exhaustive"))
t_rf <- system.time(m_rf <- glmtree(formula = fmla, data = d, restart = FALSE, testflavour = "exhaustive"))
t_rt["elapsed"] > t_rf["elapsed"] # good if TRUE

## Check if Bonferroni correction leads to a smaller tree
(m_mc <- glmtree(formula = fmla2, data = d2, family = fmly, testflavour = "ctree",
                  testtype = "MonteCarlo", bonferroni = FALSE))

(m_bmc <- glmtree(formula = fmla2, data = d2, family = fmly, testflavour = "ctree",
                   testtype = c("Bonferroni", "MonteCarlo"), bonferroni = TRUE))

width(m_mc)
width(m_bmc)

logLik(m_mc)
logLik(m_bmc)


## check multiway
(m_mult <- glmtree(formula = fmla2, data = d3, catsplit = "multiway", minsize = 80))


## check parm
fmla_p <- as.formula("y ~ x + z_noise + z_noise_1 | z + z_noise_2")
(m_interc <- glmtree(formula = fmla_p, data = d2, parm = 1))

(m_p3 <- glmtree(formula = fmla_p, data = d2, parm = 3))


## check trim
(m_tt <- glmtree(formula = fmla, data = d, trim = 0.2))

(m_tf <- glmtree(formula = fmla, data = d, trim = 300, minsize = 300))



## check breakties
m_bt <- glmtree(formula = fmla, data = d1, breakties = TRUE)
m_df <- glmtree(formula = fmla, data = d1, breakties = FALSE)

all.equal(m_bt, m_df)

unclass(m_bt)$node$info$criterion
unclass(m_df)$node$info$criterion


## check testtype
tts <- list("Bonferroni", "MonteCarlo", "Univariate", "Teststatistic",
            c("Bonferroni", "MonteCarlo"))

ms <- list()
for(tt in tts) {
  nam <- paste(tt, collapse = "")
  ms[[nam]] <- glmtree(formula = fmla, data = d, family = fmly, testflavour = "ctree",
                        testtype = tt, bonferroni = "Bonferroni" %in% tt)
}

ms





## check splittry
(m_s2 <- glmtree(formula = fmla2, data = d2, splittry = 1, bonferroni = FALSE,
                  testflavour = "ctree", testtype = "Univariate"))
(m_s5 <- glmtree(formula = fmla2, data = d2, splittry = 5, bonferroni = FALSE,
                  testflavour = "ctree", testtype = "Univariate"))



## check lookahead
smpl <- sample(1:NROW(d), size = 15)
(m_l <- glmtree(formula = fmla, data = d[smpl, ], family = fmly,
                 lookahead = TRUE,
                 testflavour = "exhaustive",
                 minbucket = 2, minsplit = 2))



### example from mob vignette
data("PimaIndiansDiabetes", package = "mlbench")

logit <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
  glm(y ~ 0 + x, family = binomial, start = start, ...)
}

pid_formula <- diabetes ~ glucose | pregnant + pressure + triceps +
  insulin + mass + pedigree + age

pid_tree <- mob(pid_formula, data = PimaIndiansDiabetes, fit = logit)
pid_tree
nodeapply(pid_tree, ids = nodeids(pid_tree), function(n) n$info$criterion)



