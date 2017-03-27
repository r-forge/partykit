
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
(mw1 <- glmtree2(formula = fmla, data = d, weights = w))
(mw2 <- glmtree2(formula = fmla, data = d, weights = w, caseweights = FALSE))



## check dfsplit
(mmfluc2 <- mob2(formula = fmla, data = d, fit = partykit:::glmfit))
(mmfluc3 <- glmtree2(formula = fmla, data = d))
(mmfluc4 <- glmtree(formula = fmla, data = d))
(mmfluc3_dfsplit <- glmtree2(formula = fmla, data = d, dfsplit = 10))
(mmfluc4_dfsplit <- glmtree(formula = fmla, data = d, dfsplit = 10))


## check tests
library("strucchange")
sctest(mmfluc3, node = 1) # does not yet work
sctest(mmfluc4, node = 1)

x <- mmfluc3
tst3 <- nodeapply(x, ids = nodeids(x), function(n) n$info$criterion)
x <- mmfluc4
tst4 <- nodeapply(x, ids = nodeids(x), function(n) n$info$test)

# should be the same, is not -> TODO: figure out why
lapply(nodeids(x), function(i) cbind(tst3[[i]][c("statistic", "p.value"), 
                                               c("z", "z_noise")],
                                     tst4[[i]]))


## check logLik and AIC
logLik(mmfluc2)
logLik(mmfluc3)
logLik(mmfluc4)
logLik(mmfluc3_dfsplit)
logLik(mmfluc4_dfsplit)

AIC(mmfluc3)
AIC(mmfluc3_dfsplit)
AIC(mmfluc4)
AIC(mmfluc4_dfsplit)



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
  x <- glmtree2(formula = fmla, data = d, inner = o)
  str(nodeapply(x, ids = nodeids(x), function(n) n$info[c("object", "estfun")]), 2)
}

for (o in options) {
  print(o)
  x <- glmtree2(formula = fmla, data = d, terminal = o)
  str(nodeapply(x, ids = nodeids(x), function(n) n$info[c("object", "estfun")]), 2)
}


## check model
m_mt <- glmtree2(formula = fmla, data = d, model = TRUE)
m_mf <- glmtree2(formula = fmla, data = d, model = FALSE)

dim(m_mt$data)
dim(m_mf$data)


## Check restart (restart = TRUE should take longer, but results should be the same)
t_rt <- system.time(m_rt <- glmtree2(formula = fmla, data = d, restart = TRUE, testflavour = "exhaustive"))
t_rf <- system.time(m_rf <- glmtree2(formula = fmla, data = d, restart = FALSE, testflavour = "exhaustive"))
t_rt["elapsed"] > t_rf["elapsed"] # good if TRUE

## Check if Bonferroni correction leads to a smaller tree
(m_mc <- glmtree2(formula = fmla2, data = d2, family = fmly, testflavour = "ctree",
                  testtype = "MonteCarlo", bonferroni = FALSE))

(m_bmc <- glmtree2(formula = fmla2, data = d2, family = fmly, testflavour = "ctree",
                   testtype = c("Bonferroni", "MonteCarlo"), bonferroni = TRUE))

width(m_mc)
width(m_bmc)

logLik(m_mc)
logLik(m_bmc)


## check multiway
(m_mult <- glmtree2(formula = fmla2, data = d3, catsplit = "multiway", minsize = 80))
(mo_mult <- glmtree(formula = fmla2, data = d3, catsplit = "multiway", minsize = 80))


## check parm
fmla_p <- as.formula("y ~ x + z_noise + z_noise_1 | z + z_noise_2")
(m_interc <- glmtree2(formula = fmla_p, data = d2, parm = 1))
(mo_interc <- glmtree(formula = fmla_p, data = d2, parm = 1))

(m_p3 <- glmtree2(formula = fmla_p, data = d2, parm = 3))
(mo_p3 <- glmtree(formula = fmla_p, data = d2, parm = 3))


## check trim 
(m_tt <- glmtree2(formula = fmla, data = d, trim = 0.2))
(mo_tt <- glmtree(formula = fmla, data = d, trim = 0.2))

(m_tf <- glmtree2(formula = fmla, data = d, trim = 300, minsize = 300))
(mo_tf <- glmtree(formula = fmla, data = d, trim = 300, minsize = 300))



## check breakties
m_bt <- glmtree2(formula = fmla, data = d1, breakties = TRUE)
m_df <- glmtree2(formula = fmla, data = d1, breakties = FALSE)

all.equal(m_bt, m_df)

unclass(m_bt)$node$info$criterion
unclass(m_df)$node$info$criterion


## check testtype
tts <- list("Bonferroni", "MonteCarlo", "Univariate", "Teststatistic",
            c("Bonferroni", "MonteCarlo"))

ms <- list()
for(tt in tts) {
  nam <- paste(tt, collapse = "")
  ms[[nam]] <- glmtree2(formula = fmla, data = d, family = fmly, testflavour = "ctree",
                        testtype = tt, bonferroni = "Bonferroni" %in% tt)
}

ms





## check splittry 
(m_s2 <- glmtree2(formula = fmla2, data = d2, splittry = 1, bonferroni = FALSE,
                  testflavour = "ctree", testtype = "Univariate"))
(m_s5 <- glmtree2(formula = fmla2, data = d2, splittry = 5, bonferroni = FALSE,
                  testflavour = "ctree", testtype = "Univariate"))



## check lookahead
smpl <- sample(1:NROW(d), size = 15)
(m_l <- glmtree2(formula = fmla, data = d[smpl, ], family = fmly,
                 lookahead = TRUE, 
                 testflavour = "exhaustive",
                 minbucket = 2, minsplit = 2))





## ctree like testing
glmtree_args_c <- c(glmtree_args, testflavour = "ctree")

mob_ctrl_c <- partykit:::mob2_control(testflavour = "ctree")
mob_ctrl_c$family <- fmly
mob_args_c <- list(fit = fit,
                   formula = fmla,
                   data = d,
                   control = mob_ctrl_c)

(m_glmtree2_c <- do.call(glmtree2, args = glmtree_args_c))
# (m_mob2_c <- do.call(mob2, args = mob_args_c))


## exhaustive search
glmtree_args_e <- c(glmtree_args, testflavour = "exhaustive")

mob_ctrl_e <- partykit:::mob2_control(testflavour = "exhaustive")
mob_ctrl_e$family <- fmly
mob_args_e <- list(fit = fit,
                   formula = fmla,
                   data = d,
                   control = mob_ctrl_e)

(m_glmtree2_e <- do.call(glmtree2, args = glmtree_args_e))
m_mob2_e <- do.call(mob2, args = mob_args_e)


width(m_glmtree2_e)
width(m_mob2_e)



### example from mob vignette
data("PimaIndiansDiabetes", package = "mlbench")

logit <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
  glm(y ~ 0 + x, family = binomial, start = start, ...)
}

pid_formula <- diabetes ~ glucose | pregnant + pressure + triceps +
  insulin + mass + pedigree + age

pid_tree <- mob(pid_formula, data = PimaIndiansDiabetes, fit = logit)
pid_tree
nodeapply(pid_tree, ids = nodeids(pid_tree), function(n) n$info$test)

pid_tree2 <- mob2(pid_formula, data = PimaIndiansDiabetes, fit = logit)
pid_tree2
nodeapply(pid_tree2, ids = nodeids(pid_tree2), function(n) n$info$criterion)


