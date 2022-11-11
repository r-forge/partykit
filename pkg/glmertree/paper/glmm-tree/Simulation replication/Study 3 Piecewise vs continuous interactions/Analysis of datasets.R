#install.packages("glmertree", repos="http://R-Forge.R-project.org")
library(glmertree)
library(partykit)
library(lme4)
library(lmerTest)

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
    if (descriptions[[i]][[3]] == "number of covariates = 15") {
      GLMMtreeformula <- Y ~ T | cluster | X1 + X2 + X3 + X4 + X5 + X6 + X7 + 
        X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15
    }
    GLMMtrees[[i]] <- lmertree(GLMMtreeformula, data=datasets[[i]], maxdepth=4)
  }
  save(GLMMtrees, file=paste("GLMMtrees", c, sep=""))
  rm(list=ls())
}


## Fit GLMMs on training data:
for (c in 1:50) {
  print(c)
  load(paste("datasets", c, sep=""))
  load(paste("descriptions", c, sep=""))
  GLMMs <- list()
  GLMMcoefs <- list()
  GLMMs2 <- list()
  varmeans <-  list()
  for (i in 1:length(datasets)) {
    print(i)
    dataset <- datasets[[i]]
    if (descriptions[[i]][[3]] == "number of covariates = 5") {
    GLMMformula <- Y ~ (1 | cluster) + T + X1 + X2 + X3 + X4 + X5 + 
      T:(X1 + X2 + X3 + X4 + X5) + (X1 + X2 + X3 + X4 + X5)^2 + 
      T:(X1 + X2 + X3 + X4 + X5)^2
      varmeans[[i]] <- attr(scale(dataset[,paste("X", 1:5, sep="")], 
                                  scale = FALSE, center = TRUE), "scaled:center")
      dataset[,paste("X", 1:5, sep = "")] <- scale(dataset[,paste("X", 1:5, sep = "")], 
                                                    scale = FALSE, center = TRUE)

    }
    if (descriptions[[i]][[3]] == "number of covariates = 15") {
      if(descriptions[[i]][[1]] == "N = 200") {
        GLMMformula <- Y ~ (1 | cluster) + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + 
          X11 + X12 + T +
          T:(X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12) + 
          (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12)^2 +        
          T:(X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12)^2
        varmeans[[i]] <- attr(scale(dataset[,paste("X", 1:12, sep="")], 
                                    scale = FALSE, center = TRUE), "scaled:center")
        dataset[,paste("X", 1:12, sep = "")] <- scale(dataset[,paste("X", 1:12, sep = "")], 
                                                      scale = FALSE, center = TRUE)
      } else {
        GLMMformula <- Y ~ (1 | cluster) + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + 
          X11 + X12 + X13 + X14 + X15 + T +
          T:(X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15) + 
          (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15)^2 +        
          T:(X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15)^2
        varmeans[[i]] <- attr(scale(dataset[,paste("X", 1:15, sep="")], 
                                    scale = FALSE, center = TRUE), "scaled:center")
        dataset[,paste("X", 1:15, sep = "")] <- scale(dataset[,paste("X", 1:15, sep = "")], 
                                                      scale = FALSE, center = TRUE)
      }
    }
    dataset$T <- as.numeric(dataset$T)-1  
    # Fit full GLMM model:
    GLMMs[[i]] <- lmer(GLMMformula, data = dataset)
    GLMMcoefs[[i]] <- summary(GLMMs[[i]])$coefficients
    # select only significant fixed-effect predictors:
    sig_pred_names <- rownames(GLMMcoefs[[i]][-1,])[which(GLMMcoefs[[i]][-1,][,"Pr(>|t|)"] < .05)]  
    # Fit new GLMM with only significant fixed predictor variables:
    if(length(sig_pred_names) > 0) {
      GLMMs2[[i]] <- lmer(paste("Y ~ (1|cluster) + ",
                                paste(sig_pred_names, collapse = " + ")),
                          data = dataset)
    } else {
      GLMMs2[[i]] <- lmer(Y ~ (1|cluster), data = dataset)
    }
  }
  save(GLMMcoefs, file = paste("GLMMcoefs", c, sep = ""))
  save(GLMMs2, file = paste("GLMMsalt", c, sep = ""))
  save(varmeans, file = paste("varmeans", c, sep = ""))
  rm(list=ls())
}


## Get tree and model characteristics:
for(c in 1:50){
  print(c)
  load(paste("GLMMtrees", c, sep=""))
  load(paste("GLMMcoefs", c, sep=""))
  load(paste("datasets", c, sep=""))
  load(paste("descriptions", c, sep=""))
  # Characteristics of GLMMs: (get the names of the significant fixed-effect predictor variables):
  names.GLMMs <- list()
  for(i in 1:length(datasets)) {
    names.GLMMs[[i]] <- rownames(GLMMcoefs[[i]][GLMMcoefs[[i]][,"Pr(>|t|)"]<.05,])
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

for(c in 1:50){
  print(c)
  load(paste("GLMMsalt",c,sep=""))
  load(paste("GLMMtrees",c,sep=""))
  load(paste("testdatasets", c, sep=""))
  load(paste("testdescriptions", c, sep=""))
  load(paste("varmeans", c, sep = ""))
  ## Compare treatment difference estimates of GMMM, GLM tree and GLMM tree, and true difference:
  treatdiffs <- list()
  for (i in 1:length(testdata)) {
    #diff <- as.numeric(substr(testdescriptions[[i]][[4]], 31, 33))
    # Calculate true treatment effect differences:
    tmp <- testdata[[i]][[1]]
    tmp$true_d_hat <- NA
    if(testdescriptions[[i]][[4]] == "type = linear") {
      X1X2 <- (tmp$X1-10)*(tmp$X2-30)
      X2X5 <- (tmp$X2-30)*(tmp$X5-70)
      tmp$true_d_hat <- -0.151*X1X2 + 0.151*X2X5
    }
    if(testdescriptions[[i]][[4]] == "type = both") {
      tmp <- testdata[[i]][[1]]
      X1X2 <- (tmp$X1-10)*(tmp$X2-30)
      X2X5 <- (tmp$X2-30)*(tmp$X5-70)
      X1X2d <- ifelse(tmp$X2<=30 & tmp$X1<=17, 1, 0)
      X2X5d <- ifelse(tmp$X2>30 & tmp$X5>63, 1, 0)
      X1X2.2 <- X1X2*X1X2d
      X2X5.2 <- X2X5*X2X5d
      tmp$true_d_hat <- -0.151*X1X2.2 + 0.151*X2X5.2
    }
    if(testdescriptions[[i]][[4]] == "type = piecewise") {
      r <- testdata[[i]][[1]]
      node3 = tmp$X2<=30 & tmp$X1<=17
      node4 = tmp$X2<=30 & tmp$X1>17
      node6 = tmp$X2>30 & tmp$X5<=63
      node7 = tmp$X2>30 & tmp$X5>63
      tmp$true_d_hat <- -5*node3 + 5*node7
    }
    treatdiffs[[i]] <- data.frame(true_d_hat=tmp$true_d_hat)
    # Calculate predicted treatment effect differences:
    treatdiffs[[i]]$GLMMtree_d_hat <-
      predict(GLMMtrees[[i]], newdata = testdata[[i]][[2]], type = "response", allow.new.levels = TRUE) -
      predict(GLMMtrees[[i]], newdata = testdata[[i]][[1]], type = "response", allow.new.levels = TRUE)
    
    dataset <- testdata[[i]]
    dataset[[1]]$T <- 0
    dataset[[2]]$T <- 1
    if (testdescriptions[[i]][[3]] == "number of covariates = 5") {
      dataset[[1]][,paste("X", 1:5, sep = "")] <- scale(dataset[[1]][,paste("X", 1:5, sep="")], 
                                  scale = FALSE, center = varmeans[[i]])
      dataset[[2]][,paste("X", 1:5, sep = "")] <- scale(dataset[[2]][,paste("X", 1:5, sep="")], 
                                                        scale = FALSE, center = varmeans[[i]])
    }
    if (testdescriptions[[i]][[3]] == "number of covariates = 15") {
      if(testdescriptions[[i]][[1]] == "N = 200") {
        dataset[[1]][,paste("X", 1:12, sep = "")] <- scale(dataset[[1]][,paste("X", 1:12, sep="")], 
                                    scale = FALSE,  center = varmeans[[i]])
        dataset[[2]][,paste("X", 1:12, sep = "")] <- scale(dataset[[2]][,paste("X", 1:12, sep="")], 
                                                      scale = FALSE,  center = varmeans[[i]])
      } else {
        dataset[[1]][,paste("X", 1:15, sep = "")] <- scale(dataset[[1]][,paste("X", 1:15, sep = "")], 
                                                      scale = FALSE, center = varmeans[[i]])
        dataset[[2]][,paste("X", 1:15, sep = "")] <- scale(dataset[[2]][,paste("X", 1:15, sep = "")], 
                                                           scale = FALSE, center = varmeans[[i]])
      }
    }
    treatdiffs[[i]]$GLMM_y_hat0 <- predict(GLMMs2[[i]], newdata=dataset[[1]], allow.new.levels = TRUE)
    treatdiffs[[i]]$GLMM_y_hat1 <- predict(GLMMs2[[i]], newdata=dataset[[2]], allow.new.levels = TRUE)
    treatdiffs[[i]]$GLMM_d_hat <- treatdiffs[[i]]$GLMM_y_hat1 - treatdiffs[[i]]$GLMM_y_hat0
  }
  save(treatdiffs, file=paste("treatdiffs", c, sep=""))
}