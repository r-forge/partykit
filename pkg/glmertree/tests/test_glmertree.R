library("glmertree")

## load artificial example data
data("DepressionDemo", package = "glmertree")

## fit normal linear regression LMM tree for continuous outcome
lt <- lmertree(depression ~ treatment | cluster | age + anxiety + duration,
               data = DepressionDemo)
#print(lt)
#plot(lt, which = "all") # default behavior, may also be "tree" or "ranef" 
coef(lt)
#ranef(lt)
#predict(lt, type = "response") # default behavior, may also be "node"
predict(lt, re.form = NA) # excludes random effects, see ?lme4::predict.merMod
#residuals(lt)
VarCorr(lt) # see lme4::VarCorr


## fit logistic regression GLMM tree for binary outcome
gt <- glmertree(depression_bin ~ treatment | cluster | age + anxiety + duration,
                data = DepressionDemo)
#print(gt)
#plot(gt, which = "all") # default behavior, may also be "tree" or "ranef" 
coef(gt)
#ranef(gt)
#predict(gt, type = "response") # default behavior, may also be "node" or "link"
predict(gt, re.form = NA) # excludes random effects, see ?lme4::predict.merMod
#residuals(gt)
VarCorr(gt) # see lme4::VarCorr