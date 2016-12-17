library(foreign)
library(glmertree)
library(lmerTest)
source("C:\\Users\\tobii\\Desktop\\swReg\\glmertree\\R\\glmertree.R")

# DATA PREPARATION
metadata <- read.dta("Database IPDMA CBT PHA Version 11.dta")
metadata[metadata == 999] <- NA
metadata[metadata == 888] <- NA
vars <- c("studyid", "Tx_group", "Age", "Gender", "education", "ComorbidAnxietyDisorder", "HRSDt0", "HRSDt1")
factors <- c("studyid", "Tx_group", "Gender", "education", "ComorbidAnxietyDisorder")
metadata$education <- factor(metadata$education, ordered = T)
for (i in 1:length(factors)) {metadata[,factors[i]] <- factor(metadata[,factors[i]])}
metadata <- metadata[vars] # select only relevant variables
metadata <- metadata[complete.cases(metadata[,vars]),] # select only complete data
metadata <- metadata[!metadata$Tx_group == "placebo",] # remove placebo observations
metadata$Tx_group <- factor(metadata$Tx_group)
summary(metadata)

## calculate trees
hrsd <- lm(HRSDt1 ~ HRSDt0, data = metadata) 
metadata$HRSDfit <- fitted(hrsd) 
lm_f <- lmtree(HRSDt1 ~ Tx_group + offset(HRSDfit) | Age + Gender +
                 education + ComorbidAnxietyDisorder + HRSDt0, data = metadata)
lm_o <- lmtree(HRSDt1 ~ Tx_group | Age + Gender +
                 education + ComorbidAnxietyDisorder + HRSDt0, offset = HRSDfit, data = metadata)
## when joint=TRUE, prediction for random effects becomes cumbersome
#lmer_f <- lmertree(HRSDt1 ~ Tx_group | (1 | studyid) + offset(HRSDfit) |
#                     Age + Gender + education + ComorbidAnxietyDisorder + HRSDt0,
#                   data = metadata, ranefstart = metadata$HRSDfit)
lmer_f <- lmertree(HRSDt1 ~ Tx_group | (1 | studyid) + offset(HRSDfit) |
                     Age + Gender + education + ComorbidAnxietyDisorder + HRSDt0,
                   data = metadata, ranefstart = metadata$HRSDfit)
plot(lm_f); print(lm_f)
plot(lm_o); print(lm_o)
plot(lmer_f); print(lmer_f)

## Fit GLMM wit h pre-specified interactions:
metadata2 <- metadata
metadata2$education <- as.numeric(metadata2$education)
GLMM <- lmer(HRSDt1 - HRSDfit ~ (1 | studyid) +  Tx_group*(Age + Gender + education + 
           ComorbidAnxietyDisorder + HRSDt0), data = metadata2)
summary(GLMM)



## calculate predictions
metadata$lmtreepred <- predict(lm_f, newdata = metadata)
metadata$lmertreepred <- predict(lmer_f, newdata = metadata) 
mean(metadata$HRSDt1); hist(metadata$HRSDt1)
mean(metadata$lmtreepred); hist(metadata$lmtreepred)
mean(metadata$lmertreepred); hist(metadata$lmertreepred)
cor(metadata$lmertreepred, metadata$HRSDt1)
cor(metadata$lmtreepred, metadata$HRSDt1)
cor(metadata$lmertreepred, metadata$lmtreepred)

## calculate amount of variance explained
lmer_f$varcorr[[1]][[1]] # variance component for random effects
var(metadata$HRSDt1) # total variance
lmer_f$varcor[[1]][[1]] / var(metadata$HRSDt1) # variance explained by random effects
var(metadata$lmertreepred) / var(metadata$HRSDt1) # Total variance explained by lmertree
var(metadata$lmtreepred) / var(metadata$HRSDt1) # variance explained by the lmtree
var(metadata$HRSDt0) / var(metadata$HRSDt1) # variance explained by HRSDt0
cor(metadata$HRSDt1, metadata$lmertreepred)^2 # lmertree explains 11% of variance in outcome variable
cor(metadata$HRSDt1, metadata$lmtreepred)^2 # lmtree explains 9% of variance outcome variable


## Calculate effect sizes: Cohen's d = x1 - x2 / pooled_sd
mu_mu <- predict(lm(HRSDt1 ~ HRSDt0, data = metadata), newdata = list(HRSDt0 = mean(metadata$HRSDt0)))
mean(metadata$HRSDt0)

lmer_es <- cbind(aggregate(metadata$HRSDt1, by = list(predict(lmer_f, type = "node")), FUN = sd), 
      m1 = coef(lmer_f)[,1] + mu_mu, m2 = coef(lmer_f)[,1] + coef(lmer_f)[,2] + mu_mu, 
      m1minm2 = coef(lmer_f)[,2])
lmer_es$cohensd <- (lmer_es$m1minm2)/lmer_es$x
lmer_es
aggregate(metadata$HRSDt1, by = list(metadata$Tx_group, predict(lmer_f, type = "node")), FUN = length)

lm_es <- cbind(aggregate(metadata$HRSDt1, by = list(predict(lm_o, type = "node")), FUN = sd), 
               m1 = coef(lm_o)[,1] + mu_mu, m2 = coef(lm_o)[,1] + coef(lm_o)[,2] + mu_mu, 
               m1minm2 = coef(lm_o)[,2])
lm_es$cohensd <- lm_es$m1minm2/lm_es$x
lm_es
aggregate(metadata$HRSDt1, by = list(metadata$Tx_group, predict(lm_o, type = "node")), FUN = length)


## 50-fold CV:
metadata <- metadata[,vars]
library(peperr)
set.seed(32)
samplething <- resample.indices(n = 694, sample.n = 50, method = "cv")
testdata <- list()
traindata <- list()
lmertrees <- list()
lmtrees <- list()
GLMMs <- list()
GLMMs2 <- list()
traindata2 <- list()
testdata2 <- list()
for (i in 1:50) {
  print(i)
  trainids <- samplething$sample.index[[i]]
  testids <- samplething$not.in.sample[[i]]
  traindata[[i]] <- metadata[trainids,]
  testdata[[i]] <- metadata[testids,]
  # fit lm(HRSDt1~HRSDt0) and predictions, based on training data 
  traindata[[i]]$HRSDfit <- fitted(lm(HRSDt1 ~ HRSDt0, data = traindata[[i]]))
  testdata[[i]]$HRSDfit <- predict(lm(HRSDt1 ~ HRSDt0, data = traindata[[i]]), newdata = testdata[[i]])

  # fit GLMM:
  traindata2[[i]] <- traindata[[i]]
  traindata2[[i]]$education <- as.numeric(traindata2[[i]]$education)
  testdata2[[i]] <- testdata[[i]]
  testdata2[[i]]$education <- as.numeric(testdata2[[i]]$education)
  
  traindata2[[i]] <- data.frame(studyid = traindata2[[i]]$studyid, 
                                HRSDt1 = traindata2[[i]]$HRSDt1,
                                HRSDfit = traindata2[[i]]$HRSDfit,
                                model.matrix(~ Tx_group*(Age + Gender + education + 
                                                           ComorbidAnxietyDisorder + HRSDt0), 
                                             data  = traindata2[[i]]))
  
  testdata2[[i]] <- data.frame(studyid = testdata2[[i]]$studyid, 
                               HRSDt1 = testdata2[[i]]$HRSDt1,
                               HRSDfit = testdata2[[i]]$HRSDfit,
                               model.matrix(~Tx_group*(Age + Gender + education + 
                                                         ComorbidAnxietyDisorder + HRSDt0), 
                                            data  = testdata2[[i]]))
  
  newform <- formula(paste("(HRSDt1 - HRSDfit) ~ (1|studyid) + ", 
                           paste(names(traindata2[[i]])[-(1:4)], collapse = " + ")))
  
  GLMMs[[i]] <- lmer(newform, data = traindata2[[i]])
  
  newform2 <- formula(paste("(HRSDt1 - HRSDfit) ~ (1|studyid) + ", 
                            paste(rownames(summary(GLMMs[[1]])$coefficients)[-1][summary(GLMMs[[1]])$coefficients[-1,"Pr(>|t|)"] < .05], 
                                  collapse = " + ")))
  
  GLMMs2[[i]] <- lmer(newform2, data = traindata2[[i]])
  # calculate GLMM predictions:
  testdata[[i]]$GLMMpreds <- predict(GLMMs2[[i]], newdata = testdata2[[i]])
  
  # grow lmtree
  lmtrees[[i]] <- lmtree(HRSDt1 ~ Tx_group + offset(HRSDfit) | Age + Gender + 
                           education + ComorbidAnxietyDisorder + HRSDt0, data = traindata[[i]])  
  # calculate lmtree predictions
  testdata[[i]]$lmtreepreds <- predict(lmtrees[[i]], newdata = testdata[[i]])

  # grow lmertree
  lmertrees[[i]] <- lmertree(HRSDt1 ~ Tx_group | (1 | studyid) + offset(HRSDfit) | 
                              Age + Gender + education + ComorbidAnxietyDisorder + HRSDt0,
                            data = traindata[[i]], ranefstart = traindata[[i]]$HRSDfit)
 
  # calculate lmertree predictions
  testdata[[i]]$lmertreepreds <- predict(lmertrees[[i]], newdata = testdata[[i]])
}  



## Assess 50-fold CV results
lmtreecors <- list()
lmertreecors <- list()
GLMMcors <- list()
for (i in 1:50){
  lmertreecors[[i]] <- cor(testdata[[i]]$lmertreepreds, testdata[[i]]$HRSDt1)
  lmtreecors[[i]] <- cor(testdata[[i]]$lmtreepreds, testdata[[i]]$HRSDt1)  
  GLMMcors[[i]] <- cor(testdata[[i]]$GLMMpreds, testdata[[i]]$HRSDt1)
}
mean(unlist(lmtreecors))
mean(unlist(lmertreecors))
mean(unlist(GLMMcors))
var(unlist(lmtreecors))
var(unlist(lmertreecors))
var(unlist(GLMMcors))



## Assess tree stability:
library(stablelearner)
# create trees without offset (otherwise, stabletree gives errors):
lm_tree <- lmtree(HRSDt1 - HRSDfit ~ Tx_group | Age + Gender +
                    education + ComorbidAnxietyDisorder + HRSDt0, data = metadata)
lmer_tree <- lmertree(HRSDt1 - HRSDfit ~ Tx_group | (1 | studyid) |
                        Age + Gender + education + ComorbidAnxietyDisorder + HRSDt0,
                      data = metadata)
set.seed(15)
subsample <- function (S = 500) {
  sampfun <- function(n) replicate(S, sample(1L:n, round(.9*n), replace = FALSE))
  list(method = "Subsampling", sampler = sampfun)
}

s_lm_tree <- stabletree(lm_tree, sampler = subsample)
s_lmer_tree <- stabletree(lmer_tree, sampler = subsample)
summary(s_lm_tree)
plot(s_lm_tree)
summary(s_lmer_tree)
plot(s_lmer_tree)




