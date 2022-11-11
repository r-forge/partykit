#install.packages("glmertree", repos="http://R-Forge.R-project.org")
library(glmertree)
library(REEMtree)

## Fit GLMM trees on training data:
for (c in 1:50) {
  load(paste("datasets", c, sep = ""))
  load(paste("descriptions", c, sep = ""))
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
    GLMMtrees[[i]] <- lmertree(GLMMtreeformula, data = datasets[[i]], 
                               maxdepth = 4)
  }
  save(GLMMtrees, file = paste("GLMMtrees", c, sep = ""))
  rm(list = ls())
}

## Fit GLM trees on training data:
for (c in 1:50) {
  load(paste("datasets", c, sep = ""))
  load(paste("descriptions", c, sep = ""))
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
    GLMtrees[[i]] <- lmtree(GLMtreeformula, data = datasets[[i]], 
                            maxdepth = 4)
  }
  save(GLMtrees, file = paste("GLMtrees", c, sep = ""))
  rm(list = ls())
}

## Fit MERTs on training data:
for(c in 1:50) {
  load(paste("datasets", c, sep = ""))
  load(paste("descriptions", c, sep = ""))
  MERTs <- list()
  for (i in 1:length(datasets)) {
    print(i)
    if (descriptions[[i]][[3]] == "number of covariates = 5") {
      MERTformula <- Y ~ X1 + X2 + X3 + X4 + X5
    }
    if (descriptions[[i]][[3]] == "number of covariates = 15") {
      MERTformula <- Y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + 
        X10 + X11 + X12 + X13 + X14 + X15
    }
    MERTs[[i]] <- list()
    MERTs[[i]][[1]] <- REEMtree(MERTformula, random = ~1 | cluster, data = datasets[[i]][datasets[[i]]$T == "1",])
    MERTs[[i]][[2]] <- REEMtree(MERTformula, random = ~1 | cluster, data = datasets[[i]][datasets[[i]]$T == "2",])
  }
  save(MERTs, file = paste("MERTs", c, sep = ""))
  rm(list = ls())
}


## Get tree and model characteristics:
for(c in 1:50){
  print(c)
  load(paste("GLMMtrees", c, sep = ""))
  load(paste("GLMtrees", c, sep = ""))
  
  ## Get GLMM tree characteristics:
  treesize.GLMM <- list()
  splits.GLMM <- list()
  for (i in 1:length(GLMMtrees)) {
    treesize.GLMM[[i]] <- length(GLMMtrees[[i]]$tree)
    splits.GLMM[[i]] <- list()
    for (j in 1:length(GLMMtrees[[i]]$tree)) {
      splits.GLMM[[i]][[j]] <- c(GLMMtrees[[i]]$tree[[j]]$node$split$varid, 
                                 GLMMtrees[[i]]$tree[[j]]$node$split$breaks)
    }
  }
  save(splits.GLMM, file = paste("splits.GLMM", c, sep = ""))
  save(treesize.GLMM, file = paste("treesize.GLMM", c, sep = ""))
  
  ## Get GLM tree characteristics:
  treesize.GLM <- list()
  splits.GLM <- list()
  for (i in 1:length(GLMtrees)) {
    treesize.GLM[[i]] <- length(GLMtrees[[i]])
    splits.GLM[[i]] <- list()
    for (j in 1:length(GLMtrees[[i]])) {
      splits.GLM[[i]][[j]] <- c(GLMtrees[[i]][[j]]$node$split$varid,
                                GLMtrees[[i]][[j]]$node$split$breaks)
    }
  }
  save(splits.GLM, file = paste("splits.GLM", c, sep = ""))
  save(treesize.GLM, file = paste("treesize.GLM", c, sep = ""))
}


## Assess predictive accuracy with test data:
for(c in 1:50){
  print(c)
  load(paste("GLMtrees", c, sep = ""))
  load(paste("GLMMtrees", c, sep = ""))
  load(paste("MERTs", c, sep = ""))
  load(paste("testdatasets", c, sep = ""))
  load(paste("testdescriptions", c, sep = ""))

  ## Get treatment difference estimates:
  treatdiffs <- list()
  for (i in 1:length(testdata)) {
    diff <- as.numeric(substr(testdescriptions[[i]][[4]], 31, 33))
    tmp <- testdata[[i]][[1]]
    tmp$true_d_hat <- NA
    tmp$true_d_hat[tmp$X2 <= 30 & tmp$X1 <= 17] <- diff
    tmp$true_d_hat[tmp$X2 <= 30 & tmp$X1 > 17] <- 0
    tmp$true_d_hat[tmp$X2 > 30 & tmp$X1 <= 63] <- 0 
    tmp$true_d_hat[tmp$X2 > 30 & tmp$X5 > 63] <- -diff
    treatdiffs[[i]] <- data.frame(true_d_hat = tmp$true_d_hat)
    treatdiffs[[i]]$GLMMtree_d_hat <- 
      predict(GLMMtrees[[i]], newdata = testdata[[i]][[1]], type = "response", re.form = NA) - 
      predict(GLMMtrees[[i]], newdata = testdata[[i]][[2]], type = "response", re.form = NA)
    treatdiffs[[i]]$GLMtree_d_hat <- 
      predict(GLMtrees[[i]], newdata = testdata[[i]][[1]], type = "response", re.form = NA) - 
      predict(GLMtrees[[i]], newdata = testdata[[i]][[2]], type = "response", re.form = NA)
    treatdiffs[[i]]$MERT_d_hat <- 
      predict(MERTs[[i]][[1]], newdata = testdata[[i]][[1]], EstimateRandomEffects = FALSE) - 
      predict(MERTs[[i]][[2]], newdata = testdata[[i]][[2]], EstimateRandomEffects = FALSE)
  }
  save(treatdiffs, file = paste("treatdiffs", c, sep = ""))
}