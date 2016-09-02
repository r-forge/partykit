#install.packages("glmertree", repos="http://R-Forge.R-project.org")
library(glmertree)

## Fit GLMM trees on training data:
for (c in 1:50) {
  load(paste("datasets", c, sep=""))
  load(paste("descriptions", c, sep=""))
  GLMMtrees <- list()
  for (i in 1:length(datasets)) {
    print(i)
    if (descriptions[[i]][[3]] == "number of covariates = 5") {
      GLMMtreeformula <- Y ~ T | cluster | X1 + X2 + X3 + X4 + X5
    }
    if (descriptions[[i]][[3]] == "number of covariates = 15") {
      GLMMtreeformula <- Y ~ T | cluster | X1 + X2 + X3 + X4 + X5 + X6 + X7 + 
        X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15
    }
    GLMMtrees[[i]] <- lmertree(GLMMtreeformula, data=datasets[[i]], maxdepth=4)
  }
  save(GLMMtrees, file=paste("GLMMtrees", c, sep=""))
  rm(list=ls())
}

## Fit GLM trees on training data:
for (c in 1:50) {
  load(paste("datasets", c, sep=""))
  load(paste("descriptions", c, sep=""))
  GLMtrees <- list()
  for (i in 1:length(datasets)) {
    print(i)
    if (descriptions[[i]][[3]] == "number of covariates = 5") {
      GLMtreeformula <- Y ~ T | X1 + X2 + X3 + X4 + X5
    }
    if (descriptions[[i]][[3]] == "number of covariates = 15") {
      GLMtreeformula <- Y ~ T | X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + 
        X10 + X11 + X12 + X13 + X14 + X15
    }
    GLMtrees[[i]] <- lmtree(GLMtreeformula, data=datasets[[i]], maxdepth=4)
  }
  save(GLMtrees, file=paste("GLMtrees", c, sep=""))
  rm(list=ls())
}

## Fit GLMMs on training data:
for (c in 1:50) {
  load(paste("datasets", c, sep=""))
  load(paste("descriptions", c, sep=""))
  GLMMs <- list()
  for (i in 1:length(datasets)) {
    print(i)
    #if (descriptions[[i]][[3]] == "number of covariates = 5") {
      GLMMformula <- Y ~ (1 | cluster) + T + X1 + X2 + X3 + X4 + X5 + 
        T:(X1 + X2 + X3 + X4 + X5) + (X1 + X2 + X3 + X4 + X5)^2 + 
        T:(X1 + X2 + X3 + X4 + X5)^2
    #}
    #if (descriptions[[i]][[3]] == "number of covariates = 15") {
    #  GLMMformula <- Y ~ (1 | cluster) + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + T +
    #    T:(X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15) + 
    #    (X1 + X2 + X3 + X4 + X5)^2 +        
    #    T:(X1 + X2 + X3 + X4 + X5)^2
    #}
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
  load(paste("GLMtrees", c, sep=""))
  load(paste("GLMMs", c, sep=""))
  
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

  breaks.GLMM <- list()
  varids.GLMM <- list()
  for (i in 1:length(datasets)) {
    breaks.GLMM[[i]] <- list()
    varids.GLMM[[i]] <- list()
    for (j in 1:length(GLMMtrees[[i]]$tree)) {
      varids.GLMM[[i]][[j]] <- GLMMtrees[[i]]$tree[[j]]$node$split$varid
      breaks.GLMM[[i]][[j]] <- GLMMtrees[[i]]$tree[[j]]$node$split$breaks
    }    
  }    
  splits.GLMM <- list()
  for (i in 1:length(datasets)) {
    splits.GLMM[[i]] <- list()
    for (j in 1:length(varids.GLMM[[i]])) {
      splits.GLMM[[i]][[j]] <- c(varids.GLMM[[i]][[j]], breaks.GLMM[[i]][[j]])
    }
  }
  splits.GLMM <- data.frame(varids=unlist(splits.GLMM)[1:length(unlist(splits.GLMM))%%2==1],
                               breaks=unlist(splits.GLMM)[1:length(unlist(splits.GLMM))%%2==0])
  save(splits.GLMM, file=paste("splits.GLMM", c, sep=""))

  # Characteristics of GLM trees:
  treesize.GLM <- list()
  for (i in 1:length(datasets)) {
    treesize.GLM[[i]] <- length(GLMtrees[[i]])
  }
  save(treesize.GLM, file=paste("treesize.GLM", c, sep=""))

  breaks.GLM <- list()
  varids.GLM <- list()
  for (i in 1:length(datasets)) {
    breaks.GLM[[i]] <- list()
    varids.GLM[[i]] <- list()
    for (j in 1:length(GLMtrees[[i]])) {
      varids.GLM[[i]][[j]] <- GLMtrees[[i]][[j]]$node$split$varid
      breaks.GLM[[i]][[j]] <- GLMtrees[[i]][[j]]$node$split$breaks
    }    
  }    
  splits.GLM <- list()
  for (i in 1:length(datasets)) {
    if(length(varids.GLM[[i]])>0) {
      splits.GLM[[i]] <- list()
      for (j in 1:length(varids.GLM[[i]])) {  
        splits.GLM[[i]][[j]] <- c(varids.GLM[[i]][[j]], breaks.GLM[[i]][[j]])
      }
    }
  }
  splits.GLM <- data.frame(varids=unlist(splits.GLM)[1:length(unlist(splits.GLM))%%2==1],
                            breaks=unlist(splits.GLM)[1:length(unlist(splits.GLM))%%2==0])
  save(splits.GLM, file=paste("splits.GLM", c, sep=""))
}


## Evaluate predictive accuracy with test data:
for(c in 1:50){
  print(c)
  load(paste("GLMMs",c,sep=""))
  load(paste("GLMtrees",c,sep=""))
  load(paste("GLMtrees",c,sep=""))
  load(paste("testdatasets", c, sep=""))
  load(paste("testdescriptions", c, sep=""))

  ## Compare treatment difference estimates of GMMM, GLM tree and GLMM tree, and true difference:
  treatdiffs <- list()
  for (i in 1:length(testdata)) {
    diff <- as.numeric(substr(testdescriptions[[i]][[4]], 31, 33))
    tmp <- testdata[[i]][[1]]
    tmp$true_d_hat <- NA
    tmp$true_d_hat[tmp$X2<=30 & tmp$X1<=17] <- diff
    tmp$true_d_hat[tmp$X2<=30 & tmp$X1>17] <- 0
    tmp$true_d_hat[tmp$X2>30 & tmp$X1<=63] <- 0 
    tmp$true_d_hat[tmp$X2>30 & tmp$X5>63] <- -diff
    treatdiffs[[i]] <- data.frame(true_d_hat=tmp$true_d_hat)
    treatdiffs[[i]]$GLMMtree_d_hat <- NA
    treatdiffs[[i]]$GLMMtree_d_hat <- predict(GLMMtrees[[i]], newdata=testdata[[i]][[1]], 
              type="response") - predict(GLMMtrees[[i]], newdata=testdata[[i]][[2]], type="response")
    treatdiffs[[i]]$GLMtree_d_hat <- predict(GLMtrees[[i]], newdata=testdata[[i]][[1]], 
              type="response") - predict(GLMtrees[[i]], 
              newdata=testdata[[i]][[2]], type="response")
    dataset <- testdata[[i]]
    dataset[[1]]$T <- as.numeric(dataset[[1]]$T)-1
    dataset[[2]]$T <- as.numeric(dataset[[2]]$T)-1
    treatdiffs[[i]]$GLMM_d_hat <- predict(GLMMs[[i]], newdata=dataset[[1]]) - 
      predict(GLMMs[[i]], newdata=dataset[[2]])
  }
  save(treatdiffs, file=paste("treatdiffs", c, sep=""))
}