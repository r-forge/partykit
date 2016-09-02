# install.packages("glmertree", repos="http://R-Forge.R-project.org") 
library(glmertree)
names(datasets[[1]])

# Fit trees on training data
for (c in 1:1) {
  load(paste("datasets", c, sep=""))
  load(paste("descriptions", c, sep=""))
  GLMMtrees <- list()
  GLMtrees <- list()
  for (i in 1:length(datasets)) {
    print(i)
    if (descriptions[[i]][[3]] == "number of covariates = 5") {
      GLMformula <- Y ~ T | X1 + X2 + X3 + X4 + X5
      GLMMformula <- Y ~ T | cluster | X1 + X2 + X3 + X4 + X5
    }
    if (descriptions[[i]][[3]] == "number of covariates = 15") {
      GLMformula <- Y ~ T | X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15
      GLMMformula <- Y ~ T | cluster | X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15
    }
  GLMMtrees[[i]] <- lmertree(GLMMformula, data=datasets[[i]], maxdepth=4)
  GLMtrees[[i]] <- lmtree(GLMformula, data=datasets[[i]], maxdepth=4)
  }
  save(GLMMtrees, file=paste("GLMMtrees", c, sep=""))
  save(GLMtrees, file=paste("GLMtrees", c, sep=""))
}

# Get tree characteristics
for(c in 1:50){
  load(paste("GLMMtrees", c, sep=""))
  load(paste("GLMtrees", c, sep=""))

  # characteristics of GLMM trees
  treesize.GLMM <- list()
  for (i in 1:length(datasets)) {
    treesize.GLMM[[i]] <- length(GLMMtrees[[i]]$tree)
  }
  save(treesize.GLMM, file=paste("treesize.GLMM", c, sep=""))
  # characteristics of GLM trees
  treesize.GLM <- list()
  for (i in 1:length(datasets)) {
    treesize.GLM[[i]] <- length(GLMtrees[[i]])
  }
  save(treesize.GLM, file=paste("treesize.GLM", c, sep=""))
}

# Evaluate tree sizes
GLMMtreesize <- vector()
GLMtreesize <- vector()
for (c in 1:50){
  load(paste("treesize.GLMM", c, sep=""))
  GLMMtreesize <- c(GLMMtreesize, unlist(treesize.GLMM))
  load(paste("treesize.GLM", c, sep=""))
  GLMtreesize <- c(GLMtreesize, unlist(treesize.GLM))  
}
treesizes <- data.frame(GLMMtreesize, GLMtreesize)

load("descriptions1")
for (i in 1:length(descriptions)) {
  treesizes[i,"N"] <- substr(descriptions[[i]][[1]], 5, 8) # N
  treesizes[i,"rho"] <- substr(descriptions[[i]][[2]], 7, 9) # is rho (correlation between covariates)
  treesizes[i,"np"] <- substr(descriptions[[i]][[3]], 24, 25) # is number of covariates
  treesizes[i,"treatdiff"] <- substr(descriptions[[i]][[4]], 31, 33) # is treatment effect difference
  treesizes[i,"corUbi"] <- substr(descriptions[[i]][[5]], 32, 64) # is correlation between U and bi
  treesizes[i,"numbclus"] <- substr(descriptions[[i]][[6]], 34, 35) # is number of random intercept values
  treesizes[i,"sigmabi"] <- substr(descriptions[[i]][[7]], 12, 14) # max and -min random intercept value    
}
save(treesizes,file="treesizes")

treesizes$N <- factor(treesizes$N)
treesizes$rho <-  factor(treesizes$rho)            
treesizes$np <- factor(treesizes$np)              
treesizes$treatdiff <- factor(treesizes$treatdiff)          
treesizes$corUbi <- factor(treesizes$corUbi)           
treesizes$numbclus <- factor(treesizes$numbclus)        
treesizes$sigmabi <- factor(treesizes$sigmabi)

treesizes$N <- rep(treesizes$N[1:length(descriptions)], times=1)
treesizes$rho <- rep(treesizes$rho[1:length(descriptions)], times=1)
treesizes$np <- rep(treesizes$np[1:length(descriptions)], times=1)
treesizes$treatdiff <- rep(treesizes$treatdiff[1:length(descriptions)], times=1)
treesizes$corUbi <- rep(treesizes$corUbi[1:length(descriptions)], times=1)
treesizes$numbclus <- rep(treesizes$numbclus[1:length(descriptions)], times=1)
treesizes$sigmabi <- rep(treesizes$sigmabi[1:length(descriptions)], times=1)

treesizes
table(treesizes[,1]) # GLMM trees
prop.table(table(treesizes[,1])) # GLMM trees
table(treesizes[,2]) # GLM trees
prop.table(table(treesizes[,2])) # GLM trees

treesizes.long <- data.frame(treesize=stack(treesizes[,1:2]), 
                             N=rep(treesizes$N, 2), rho=rep(treesizes$rho, 2), np=rep(treesizes$np, 2), 
                             treatdiff=rep(treesizes$treatdiff, 2),  corUbm=rep(treesizes$corUbi, 2), 
                             numbclus=rep(treesizes$numbclus, 2), sigmabm=rep(treesizes$sigmabi, 2),
                             datasetID=factor(rep(1:nrow(treesizes), 2)))
tmp <- as.character(treesizes.long$treesize.ind)

tmp[treesizes.long$treesize.ind=="GLMMtreesize"] <- "GLMM tree"
tmp[treesizes.long$treesize.ind=="GLMtreesize"] <- "GLM tree"
treesizes.long$treesize.ind <- factor(tmp)

mean(treesizes[,1]);sd(treesizes[,1]) # GLMM tree
mean(treesizes[,2]);sd(treesizes[,2]) # GLM tree

# create lattice xyplots
library(lattice)
treesize.anova <- aov(treesize.values ~ treesize.ind + N + rho + np + treatdiff + numbclus + sigmabm + corUbm +
                        treesize.ind*(N + rho + np + treatdiff + numbclus + sigmabm + corUbm), 
                      data=treesizes.long)
summary(treesize.anova)[[1]]["Sum Sq"] / sum(summary(treesize.anova)[[1]]["Sum Sq"] )
# N, sigmabm & corUbi have main and interaction effects with eta squared > .01
treesizes.long$sigmabm <- factor(as.numeric(as.character(treesizes.long$sigmabm)), ordered=T)
treesizes.long$N <- factor(as.numeric(as.character(treesizes.long$N)), ordered=T)
aggdata <- aggregate(formula=treesize.values ~ treesize.ind + N + sigmabm + corUbm, FUN=mean, data=treesizes.long)

xyplot(treesize.values ~ sigmabm | N + corUbm, data = aggdata, groups=treesize.ind, type="b",
       ylab="tree size", xlab="sigma_b", par.settings=standard.theme("pdf",color=F), abline=c(1,0), 
       auto.key=list(space="top", columns=2, title="Algorithm type", cex.title=1,lines=T, points=T))

