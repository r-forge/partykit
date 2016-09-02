#install.packages("glmertree", repos="http://R-Forge.R-project.org")
#library(glmertree)
library(partykit)
library(lme4)

# source most recent glmertree code:
source("C:/Users/tobii/Desktop/swReg/glmertree/R/glmertree.r")

## Fit GLMM trees on training data:
for (c in 1:50) {
  print(c)
  load(paste("datasets", c, sep=""))
  load(paste("descriptions", c, sep=""))
  GLMMtrees <- list()
  for (i in 1:length(datasets)) {
    print(i)
    if (descriptions[[i]][[3]] == "number of covariates = 5") {
      GLMMtreeformula <- Y ~ T | cluster | X1 + X2 + X3 + X4 + X5
    }
    if (descriptions[[i]][[3]] == "number of covariates = 10") {
      GLMMtreeformula <- Y ~ T | cluster | X1 + X2 + X3 + X4 + X5 + X6 + X7 + 
        X8 + X9 + X10
    }
    GLMMtrees[[i]] <- lmertree(GLMMtreeformula, data=datasets[[i]], maxdepth=4)
  }
  save(GLMMtrees, file=paste("GLMMtrees", c, sep=""))
  rm(list=ls())
  source("C:/Users/tobii/Desktop/swReg/glmertree/R/glmertree.r")
}


## Fit GLMMs on training data:
for (c in 1:50) {
  print(c)
  load(paste("datasets", c, sep=""))
  load(paste("descriptions", c, sep=""))
  GLMMs <- list()
  for (i in 1:length(datasets)) {
    if (descriptions[[i]][[3]] == "number of covariates = 5") {
    GLMMformula <- Y ~ (1 | cluster) + T + X1 + X2 + X3 + X4 + X5 + 
      T:(X1 + X2 + X3 + X4 + X5) + (X1 + X2 + X3 + X4 + X5)^2 + 
      T:(X1 + X2 + X3 + X4 + X5)^2
    }
    if (descriptions[[i]][[3]] == "number of covariates = 10") {
      if(descriptions[[i]][[1]] == "N = 200") {
        GLMMformula <- Y ~ (1 | cluster) + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + T +
          T:(X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10) + 
          (X1 + X2 + X3 + X4 + X5)^2 +        
          T:(X1 + X2 + X3 + X4 + X5)^2
      } else {
        GLMMformula <- Y ~ (1 | cluster) + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + T +
          T:(X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10) + 
          (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10)^2 +        
          T:(X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10)^2
      }
    }
    dataset <- datasets[[i]]
    dataset$T <- as.numeric(dataset$T)-1  
    GLMMs[[i]] <- lmer(GLMMformula, data = dataset)
  }
  save(GLMMs, file = paste("GLMMs", c, sep = ""))
  rm(list=ls())
}

## Get tree and model characteristics:
for(c in 1:50){
  print(c)
  load(paste("GLMMtrees", c, sep=""))
  load(paste("GLMMs", c, sep=""))
  load(paste("datasets", c, sep=""))
  load(paste("descriptions", c, sep=""))
  # Characteristics of GLMMs: (get the names of the significant fixed-effect predictor variables):
  names.GLMMs <- list()
  for(i in 1:length(datasets)) {
    names.GLMMs[[i]] <- names(summary(GLMMs[[i]])$coefficients[-1,"t value"][
      summary(GLMMs[[i]])$coefficients[-1,"t value"] > 1.96])
  }
  save(names.GLMMs, file = paste("names.GLMMs", c, sep = ""))
  
  # Characteristics of GLMM trees:
  treesize.GLMM <- list()
  for (i in 1:length(datasets)) {
    treesize.GLMM[[i]] <- length(GLMMtrees[[i]]$tree)
  }
  save(treesize.GLMM, file=paste("treesize.GLMM", c, sep=""))

  varids.GLMM <- list()
  for (i in 1:length(datasets)) {
    varids.GLMM[[i]] <- list()
    for (j in 1:length(GLMMtrees[[i]]$tree)) {
      varids.GLMM[[i]][[j]] <- GLMMtrees[[i]]$tree[[j]]$node$split$varid
    }    
  }    
  splits.GLMM <- data.frame(varids=unlist(varids.GLMM))
  save(splits.GLMM, file=paste("splits.GLMM", c, sep=""))
}


## Evaluate predictive accuracy with test data:
source("C:/Users/tobii/Desktop/swReg/glmertree/R/glmertree.r")

for(c in 1:50){
  print(c)
  load(paste("GLMMs",c,sep=""))
  load(paste("GLMMtrees",c,sep=""))
  load(paste("testdatasets", c, sep=""))
  load(paste("testdescriptions", c, sep=""))

  ## Compare treatment difference estimates of GMMM, GLM tree and GLMM tree, and true difference:
  treatdiffs <- list()
  for (i in 1:length(testdata)) {
    print(i)
    diff <- as.numeric(substr(testdescriptions[[i]][[4]], 31, 33))
    tmp <- testdata[[i]][[1]]
    tmp$true_d_hat <- NA
    if(testdescriptions[[i]][[4]] == "type = linear") {
      tmp$true_d_hat <- (testdata[[i]][[1]]$X2-30)*(testdata[[i]][[1]]$X1-10)*-.25 + 
        (testdata[[i]][[1]]$X2-30)*(testdata[[i]][[1]]$X5-70)*.25
    }
    if(testdescriptions[[i]][[4]] == "type = piecewise") {
      r <- testdata[[i]][[1]]
      node3 = r$X2<=30 & r$X1<=17
      node4 = r$X2<=30 & r$X1>17
      node6 = r$X2>30 & r$X5<=63
      node7 = r$X2>30 & r$X5>63
      tmp$true_d_hat <- NA
      tmp$true_d_hat[node3] <- -28.39720+20.17224
      tmp$true_d_hat[node4] <- 39.59950-13.75540
      tmp$true_d_hat[node6] <- -39.53344+13.75933
      tmp$true_d_hat[node7] <- 28.37829-20.14495
    }
    if(testdescriptions[[i]][[4]] == "type = both") {
      tmp$true_d_hat <- (testdata[[i]][[1]]$X2-30)*(testdata[[i]][[1]]$X1-10)*-.25 + 
        (testdata[[i]][[1]]$X2-30)*(testdata[[i]][[1]]$X5-70)*.25
      tmp$true_d_hat <- tmp$true_d_hat*.5
      r <- testdata[[i]][[1]]
      node3 = r$X2<=30 & r$X1<=17
      node4 = r$X2<=30 & r$X1>17
      node6 = r$X2>30 & r$X5<=63
      node7 = r$X2>30 & r$X5>63
      tmp$true_d_hat[node3] <- tmp$true_d_hat[node3]+.5*(-28.39720+20.17224)
      tmp$true_d_hat[node4] <- tmp$true_d_hat[node4]+.5*(39.59950-13.75540)
      tmp$true_d_hat[node6] <- tmp$true_d_hat[node6]+.5*(-39.53344+13.75933)
      tmp$true_d_hat[node7] <- tmp$true_d_hat[node7]+.5*(28.37829-20.14495)
    }
    
    treatdiffs[[i]] <- data.frame(true_d_hat=tmp$true_d_hat)
    treatdiffs[[i]]$GLMMtree_d_hat <-
      predict(GLMMtrees[[i]], newdata = testdata[[i]][[2]], type = "response", allow.new.levels = TRUE) -
      predict(GLMMtrees[[i]], newdata = testdata[[i]][[1]], type = "response", allow.new.levels = TRUE)
    dataset <- testdata[[i]]
    dataset[[1]]$T <- 0
    dataset[[2]]$T <- 1
    treatdiffs[[i]]$GLMM_d_hat <- predict(GLMMs[[i]], newdata=dataset[[2]], allow.new.levels = TRUE) - 
      predict(GLMMs[[i]], newdata=dataset[[1]], allow.new.levels = TRUE)
  }
  save(treatdiffs, file=paste("treatdiffs", c, sep=""))
}