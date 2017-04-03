# library("partykit")
# 
# set.seed(29)
# n <- 1000
# x <- runif(n)
# z <- runif(n)
# y <- rnorm(n, mean = x * c(-1, 1)[(z > 0.7) + 1], sd = 3)
# z_noise <- factor(sample(1:3, size = n, replace = TRUE))
# d <- data.frame(y = y, x = x, z = z, z_noise = z_noise)
# 
# 
# fmla <- as.formula("y ~ x | z + z_noise")
# 
# (m_mob <- mob2(formula = fmla, data = d, fit = partykit:::lmfit))
# (m_lm2 <- lmtree2(formula = fmla, data = d))
# (m_lm1 <- lmtree(formula = fmla, data = d))
# 
# 
# mods <- nodeapply(m_lm2, ids = nodeids(m_lm2, terminal = TRUE), function(x) x$info$object)
# sum(sapply(mods, function(x) sum(x$residuals^2)))
# 
# ## glm versus lm / logLik versus sum of squared residuals
# lm_big <- lmtree2(formula = fmla, data = d, testflavour = "exhaustive")
# glm_big <- glmtree2(formula = fmla, data = d, testflavour = "exhaustive")
# 
# width(lm_big)
# width(glm_big)
# 
# AIC(lm_big)
# AIC(glm_big)
# # AIC of lmtree must be worse (greater) or equal than of glmtree because AIC assesses the logLik, 
# # which is optimized in glm but noth in lm (there it is the sum of squared residuals).
# # In splitting the different objective functions may lead to differences, because the logLik
# # contains the term log(sum of squared residuals), i.e. we optimize
# # log(sum of squared residuals in left node) + log(sum of squared residuals in right node);
# # In lmtree we optimize
# # sum of squared residuals in left node + sum of squared residuals in right node
# 
# 
# ## pruning
# lm_aic <- prune.modelparty(lm_big, type = "AIC")
# lm_bic <- prune.modelparty(lm_big, type = "BIC")
# 
# width(lm_big)
# width(lm_aic)
# width(lm_bic)
# 
# glm_aic <- prune.modelparty(glm_big, type = "AIC")
# glm_bic <- prune.modelparty(glm_big, type = "BIC")
# 
# width(glm_big)
# width(glm_aic)
# width(glm_bic)
# 
# 
# 
# 
