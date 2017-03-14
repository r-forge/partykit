
library("partykit")
library("Formula")

set.seed(29)
n <- 1000
x <- runif(n)
z <- runif(n)
y <- rnorm(n, mean = x * c(-1, 1)[(z > 0.5) + 1], sd = 3)
z_noise <- factor(sample(1:3, size = n, replace = TRUE))
d <- data.frame(y = y, x = x, z = z, z_noise = z_noise)

fmla <- as.formula("y ~ x | z + z_noise")
fmly <- gaussian()
fit <- partykit:::.glmtrafo



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


## Check if Bonferroni correction leads to a smaller tree
k <- 20
zs_noise <- matrix(rnorm(n*k), nrow = n)
colnames(zs_noise) <- paste0("z_noise", 1:k)
d2 <- cbind(d, zs_noise)
fmla2 <- as.formula(paste("y ~ x | z + z_noise +", 
                          paste0("z_noise", 1:k, collapse = " + ")))

(m_mc <- glmtree2(formula = fmla2, data = d2, family = fmly, testflavour = "ctree",
                 testtype = "MonteCarlo", bonferroni = FALSE))

(m_bmc <- glmtree2(formula = fmla2, data = d2, family = fmly, testflavour = "ctree",
                  testtype = c("Bonferroni", "MonteCarlo"), bonferroni = TRUE))

width(m_mc)
width(m_bmc)

## check lookahead
smpl <- sample(1:NROW(d), size = 15)
(m_l <- glmtree2(formula = fmla, data = d[smpl, ], family = fmly,
                 lookahead = TRUE, 
                 testflavour = "exhaustive",
                 minbucket = 2, minsplit = 2))


## check nmax 
(m1 <- glmtree2(formula = fmla, data = d, nmax = 1000, testflavour = "ctree", splitflavour = "ctree"))
(m2 <- glmtree2(formula = fmla, data = d, nmax = Inf, testflavour = "ctree", splitflavour = "ctree"))



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






