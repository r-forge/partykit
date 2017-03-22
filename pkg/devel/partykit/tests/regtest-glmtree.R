
library("partykit")
library("sandwich")

set.seed(29)
n <- 1000
x <- runif(n)
z <- runif(n)
y <- rnorm(n, mean = x * c(-1, 1)[(z > 0.7) + 1], sd = 3)
z_noise <- factor(sample(1:3, size = n, replace = TRUE))
d <- data.frame(y = y, x = x, z = z, z_noise = z_noise)

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



fmla <- as.formula("y ~ x | z + z_noise")
fmly <- gaussian()
fit <- partykit:::glmfit


(m2 <- mob2(formula = fmla, data = d, fit = partykit:::glmfit, 
           control = partykit:::mob2_control(testflavour = "ctree",
                                             splitflavour = "ctree",
                                             bonferroni = FALSE,
                                             testtype = "Univariate")))

(m2 <- mob2(formula = fmla, data = d, fit = partykit:::glmfit, 
           control = partykit:::mob2_control(testflavour = "ctree",
                                             splitflavour = "ctree")))


(mmfluc2 <- mob2(formula = fmla, data = d, fit = partykit:::glmfit))
(mmfluc3 <- glmtree2(formula = fmla, data = d))
(mmfluc4 <- glmtree(formula = fmla, data = d))

x <- mmfluc3
info <- nodeapply(x, ids = nodeids(x, terminal = TRUE),
                  FUN = function(n) info_node(n)$objfun)

## Check if Bonferroni correction leads to a smaller tree
(m_mc <- glmtree2(formula = fmla2, data = d2, family = fmly, testflavour = "ctree",
                  testtype = "MonteCarlo", bonferroni = FALSE))

(m_bmc <- glmtree2(formula = fmla2, data = d2, family = fmly, testflavour = "ctree",
                   testtype = c("Bonferroni", "MonteCarlo"), bonferroni = TRUE))

width(m_mc)
width(m_bmc)

# logLik(m_mc)
# logLik(m_bmc)




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
                  testflavour = "ctree"))
(m_s5 <- glmtree2(formula = fmla2, data = d2, splittry = 5, bonferroni = FALSE,
                  testflavour = "ctree"))



## check lookahead
smpl <- sample(1:NROW(d), size = 15)
(m_l <- glmtree2(formula = fmla, data = d[smpl, ], family = fmly,
                 lookahead = TRUE, 
                 testflavour = "exhaustive",
                 minbucket = 2, minsplit = 2))


## default settings
glmtree_args <- list(formula = fmla, 
                     data = d,
                     family = fmly)

mob_ctrl <- partykit:::mob2_control()
mob_ctrl$family <- fmly
mob_args <- list(fit = fit,
                 formula = fmla,
                 data = d,
                 control = mob_ctrl)


(m_glmtree <- do.call(glmtree, args = glmtree_args))
(m_glmtree2 <- do.call(glmtree2, args = glmtree_args))
(m_mob2 <- do.call(mob2, args = mob_args))


## ctree like testing
glmtree_args_c <- c(glmtree_args, testflavour = "ctree")

mob_ctrl_c <- partykit:::mob2_control(testflavour = "ctree")
mob_ctrl_c$family <- fmly
mob_args_c <- list(fit = fit,
                   formula = fmla,
                   data = d,
                   control = mob_ctrl_c)

(m_glmtree2_c <- do.call(glmtree2, args = glmtree_args_c))
(m_mob2_c <- do.call(mob2, args = mob_args_c))


## exhaustive search
glmtree_args_e <- c(glmtree_args, testflavour = "exhaustive")

mob_ctrl_e <- partykit:::mob2_control(testflavour = "exhaustive")
mob_ctrl_e$family <- fmly
mob_args_e <- list(fit = fit,
                   formula = fmla,
                   data = d,
                   control = mob_ctrl_e)

(m_glmtree2_e <- do.call(glmtree2, args = glmtree_args_e))
(m_mob2_e <- do.call(mob2, args = mob_args_e))


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

pid_tree2 <- mob2(pid_formula, data = PimaIndiansDiabetes, fit = logit)
pid_tree2
