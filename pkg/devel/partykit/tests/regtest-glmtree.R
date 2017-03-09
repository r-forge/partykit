
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


system.time(m_glmtree <- do.call(glmtree, args = glmtree_args))
system.time(m_glmtree2 <- do.call(glmtree2, args = glmtree_args))
system.time(m_mob2 <- do.call(mob2, args = mob_args))


## ctree like testing
glmtree_args_c <- c(glmtree_args, testflavour = "ctree")

mob_ctrl_c <- partykit:::mob2_control(testflavour = "ctree")
mob_ctrl_c$family <- fmly
mob_args_c <- list(fit = fit,
                 formula = fmla,
                 data = d,
                 control = mob_ctrl_c)

system.time(m_glmtree2_c <- do.call(glmtree2, args = glmtree_args_c))
system.time(m_mob2_c <- do.call(mob2, args = mob_args_c))


## exhaustive search
glmtree_args_e <- c(glmtree_args, testflavour = "exhaustive")

mob_ctrl_e <- partykit:::mob2_control(testflavour = "exhaustive")
mob_ctrl_e$family <- fmly
mob_args_e <- list(fit = fit,
                   formula = fmla,
                   data = d,
                   control = mob_ctrl_e)

system.time(m_glmtree2_e <- do.call(glmtree2, args = glmtree_args_e))
system.time(m_mob2_e <- do.call(mob2, args = mob_args_e))


m_glmtree2_e_small <- glmtree2(fmla, data = d, testflavour = "exhaustive", maxdepth = 2)
m_glmtree2_c_small <- glmtree2(fmla, data = d, testflavour = "ctree", maxdepth = 2)
m_glmtree2_e_small
m_glmtree2_c_small

## compare
m_glmtree
m_glmtree2
m_mob2

m_glmtree2_c
m_mob2_c

m_glmtree2_e
m_mob2_e

width(m_glmtree2_e)
width(m_mob2_e)

## check lookahead
# smpl <- sample(1:NROW(d), size = 15)
# m_l <- glmtree2(formula = fmla, data = d[smpl, ], family = fmly, 
#                 lookahead = TRUE, testflavour = "exhaustive", 
#                 minbucket = 6, minsplit = 2)
# 
# m_d <- glmtree2(formula = fmla, data = d, family = fmly, 
#                  testflavour = "exhaustive", 
#                 minbucket = 2, minsplit = 2)
