# install.packages("glmertree", repos="http://R-Forge.R-project.org")
library(glmertree)


## Evaluate tree size:

GLMMtreesize <- vector()
for (c in 1:50){
  load(paste("treesize.GLMM", c, sep=""))
  GLMMtreesize <- c(GLMMtreesize, unlist(treesize.GLMM))
}
treespecs <- data.frame(GLMMtreesize)

load("descriptions1")
for (i in 1:length(descriptions)) {
  treespecs[i,"N"] <- substr(descriptions[[i]][[1]], 5, 8)
  treespecs[i,"rho"] <- substr(descriptions[[i]][[2]], 7, 9)
  treespecs[i,"np"] <- substr(descriptions[[i]][[3]], 24, 25)
  treespecs[i,"type"] <- substr(descriptions[[i]][[4]], 8, 18)
  treespecs[i,"corUbi"] <- substr(descriptions[[i]][[5]], 32, 64)
  treespecs[i,"numbclus"] <- substr(descriptions[[i]][[6]], 34, 35)
  treespecs[i,"sigmabi"] <- substr(descriptions[[i]][[7]], 12, 14)  
}

treespecs$N <- rep(factor(treespecs$N[1:length(descriptions)]), times=50)
treespecs$rho <- rep(factor(treespecs$rho[1:length(descriptions)]), times=50)
treespecs$np <- rep(factor(treespecs$np[1:length(descriptions)]), times=50)
treespecs$type <- rep(factor(treespecs$type[1:length(descriptions)]), times=50)
treespecs$corUbi <- rep(factor(treespecs$corUbi[1:length(descriptions)]), times=50)
treespecs$numbclus <- rep(factor(treespecs$numbclus[1:length(descriptions)]), times=50)
treespecs$sigmabi <- rep(factor(treespecs$sigmabi[1:length(descriptions)]), times=50)

treespecs.long <- data.frame(treesize=treespecs["GLMMtreesize"], N=treespecs$N, rho=treespecs$rho, 
                             np=treespecs$np, type=treespecs$type, corUb=treespecs$corUbi, 
                             numbclus=treespecs$numbclus, sigmab=treespecs$sigmabi, datasetID=
                               factor(1:nrow(treespecs)))
treespecs.long <- rbind(treespecs.long, treespecs.long)

s# For GLMMs:
treespecs$noGLMMpreds <- NA
for(c in 1:50) {
  load(paste("names.GLMMs", c, sep=""))
       for(i in 1:length(names.GLMMs)){
         treespecs$noGLMMpreds[(c-1)*length(descriptions) + i] <- 
           length(names.GLMMs[[i]])
  }
}
mean(treespecs$noGLMMpreds)
tapply(treespecs$noGLMMpreds, treespecs$type, mean)
tapply(treespecs$noGLMMpreds, treespecs$type, sd)

mean(treespecs$GLMMtreesize)
(mean(treespecs$GLMMtreesize) - 1) / 2
tapply(treespecs$GLMMtreesize, treespecs$type, mean)
tapply(treespecs$GLMMtreesize, treespecs$type, sd)

## Evaluate predictive accuracy:

treespecs$cor_GLMMtree_true <- NA
treespecs$cor_GLMM_true <- NA
treespecs$mean_true_diffs <- NA
treespecs$mean_GLMMtree_diffs <- NA
treespecs$mean_GLMM_diffs <- NA
treespecs$sd_true_diffs <- NA
treespecs$sd_GLMMtree_diffs <- NA
treespecs$sd_GLMM_diffs <- NA

for(c in 1:50) {
  print(c)
  load(paste("treatdiffs", c, sep=""))
  for(i in 1:length(descriptions)) {
    treespecs[(c-1)*length(descriptions)+i, "cor_GLMMtree_true"] <- 
      with(treatdiffs[[i]], cor(GLMMtree_d_hat, true_d_hat))
    treespecs[(c-1)*length(descriptions)+i, "cor_GLMM_true"] <- 
      with(treatdiffs[[i]], cor(GLMM_d_hat, true_d_hat))
    treespecs[(c-1)*length(descriptions)+i, "mean_true_diffs"] <-
      with(treatdiffs[[i]], mean(true_d_hat))
    treespecs[(c-1)*length(descriptions)+i, "mean_GLMMtree_diffs"] <- 
      with(treatdiffs[[i]], mean(GLMMtree_d_hat))
    treespecs[(c-1)*length(descriptions)+i, "mean_GLMM_diffs"] <- 
      with(treatdiffs[[i]], mean(GLMM_d_hat))
    treespecs[(c-1)*length(descriptions)+i, "sd_true_diffs"] <-
      with(treatdiffs[[i]], sd(true_d_hat))
    treespecs[(c-1)*length(descriptions)+i, "sd_GLMMtree_diffs"] <- 
      with(treatdiffs[[i]], sd(GLMMtree_d_hat))
    treespecs[(c-1)*length(descriptions)+i, "sd_GLMM_diffs"] <- 
      with(treatdiffs[[i]], sd(GLMM_d_hat))
  }
}

# if a correlation cannot be calculated, the tree has 1 node, and cor should be 0:
table(is.na(treespecs$cor_GLMMtree_true))
treespecs$cor_GLMMtree_true[is.na(treespecs$cor_GLMMtree_true)] <- 0
table(is.na(treespecs$cor_GLMM_true))
treespecs$cor_GLMM_true[is.na(treespecs$cor_GLMM_true)] <- 0


# Get descriptives:
apply(treespecs[,9:16], 2, mean)
apply(treespecs[,9:16], 2, sd)
treespecs$type <- factor(treespecs$type, ordered = TRUE, levels = 
                           c("linear","both","piecewise"))
tapply(treespecs[,9], treespecs$type, mean, na.rm = TRUE) # GLMM tree
tapply(treespecs[,10], treespecs$type, mean, na.rm = TRUE) # GLMM

treespecs.long <- data.frame(correlation = stack(treespecs[c("cor_GLMMtree_true", "cor_GLMM_true")]), 
                             treespecs.long)
tmp <- as.character(treespecs.long$correlation.ind)
tmp[treespecs.long$correlation.ind=="cor_GLMMtree_true"] <- "GLMM tree"
tmp[treespecs.long$correlation.ind=="cor_GLMM_true"] <- "GLMM"
treespecs.long$correlation.ind <- factor(tmp)



## Evaluate tree accuracy:

# Recover splitting variables and split points:

tmp <- data.frame()
for (c in 1:50){
  load(paste("splits.GLMM",c,sep=""))
  tmp <- data.frame(rbind(tmp, splits.GLMM))
}
GLMMsplits <- data.frame(GLMMsplitvar=tmp[,1])

tmp <- vector() 
tmpnodeno <- vector()
for (j in 1:nrow(treespecs)) {
  print(j)
  tmp <- c(tmp, rep(treespecs$GLMMtreesize[j], times=(treespecs$GLMMtreesize[j]-1)/2))
  if(treespecs$GLMMtreesize[j] > 1) {
    tmpnodeno <- c(tmpnodeno, seq(1,(treespecs$GLMMtreesize[j]-1)/2))
  }
} 
GLMMsplits$treesize <- tmp
GLMMsplits$nodenumber <- tmpnodeno
GLMMsplits$GLMMsplitvar <- factor(GLMMsplits$GLMMsplitvar)
save(GLMMsplits, file="GLMMsplits")
load("GLMMsplits")



# Assess accuracy of splits:

# Accuracy of first split:
table(GLMMsplits[GLMMsplits$nodenumber==1, "GLMMsplitvar"])

# identify which split belongs to which tree:
GLMMtreeno <- rep(1:7200, times=(treespecs$GLMMtreesize-1)/2)
GLMMsplits$treeno <- GLMMtreeno



# Check whether the right predictor variables were recovered:
treespecs$treeno <- 1:7200
treespecs$trueGLMMtree <- NA
treespecs$trueGLMM <- NA

## Check resulting trees: 
# For piecewise interactions, true tree should be recovered:
GLMMtruefirstsplit <- GLMMsplits[GLMMsplits$GLMMsplitvar==4 & GLMMsplits$nodenumber == 1 &
                                 GLMMsplits$treesize==7, "treeno"]
GLMMtruesecsplit <- GLMMsplits[GLMMsplits$GLMMsplitvar==3 & GLMMsplits$nodenumber == 2 &
                               GLMMsplits$treesize==7, "treeno"]
GLMMtruethirdsplit <- GLMMsplits[GLMMsplits$GLMMsplitvar==7 & GLMMsplits$nodenumber == 3 &
                                 GLMMsplits$treesize==7, "treeno"]
GLMMtruetreenos <- GLMMtruefirstsplit[GLMMtruefirstsplit %in% GLMMtruesecsplit][
      GLMMtruefirstsplit[GLMMtruefirstsplit %in% GLMMtruesecsplit] %in% GLMMtruethirdsplit]
treespecs$trueGLMMtree[GLMMtruetreenos] <- TRUE 

tmp <- treespecs$trueGLMMtree
# In linear and combined datasets, check whether X2, X1 and X5 are all involved in tree:
for(i in c(which(treespecs$type == "linear"), which(treespecs$type == "both"))) {
  treespecs$trueGLMMtree[i] <- 
    all(c(3,4,7) %in% GLMMsplits[GLMMsplits$treeno==i,"GLMMsplitvar"]) & 
      !any(c(5,6,8:17) %in% GLMMsplits[GLMMsplits$treeno==i,"GLMMsplitvar"])
}

# For GLMMs:
for(c in 1:50) {
  load(paste("names.GLMMs", c, sep=""))
  print(c)
  for(i in 1:length(names.GLMMs)){
    treespecs$trueGLMM[(c-1)*length(descriptions) + i] <-
    sum(c("X1:X2", "X2:X5", "X2:X1", "X5:X2", "X1:X2:T", "X2:X5:T", "T:X1:X2", "T:X2:X5",
          "T:X2:X1", "T:X5:X2", "X1:T:X2", "X2:T:X5", "X2:T:X1", "X5:T:X2") 
        %in% names.GLMMs[[i]]) == 4
  }
}

# Put results together in dataframe and save it:
treespecs.long <- data.frame(truevars = c(treespecs$trueGLMMtree, treespecs$trueGLMM),
                             truevars.ind = rep(c("GLMMtree","GLMM"), each = 7200), treespecs.long)
mean(treespecs$trueGLMMtree)
mean(treespecs$trueGLMM)
save("treespecs.long", file="treespecs_studyIII.dat")
save("treespecs", file="treespecs.dat")










## Create plots of outcomes:
load("treespecs_studyIII.dat")
library(lattice)

tapply(treespecs.long$truevars, treespecs.long$truevars.ind, mean)
tapply(treespecs.long$correlation.values, treespecs.long$correlation.ind, mean, na.rm = TRUE)
tapply(treespecs.long$correlation.values, treespecs.long$correlation.ind, sd, na.rm = TRUE)

treespecs.long$N <- factor(as.numeric(as.character(treespecs.long$N)), ordered=T)
treespecs.long$sigmab <- factor(as.numeric(as.character(treespecs.long$sigmab)), ordered=T)
treespecs.long$type <- factor(treespecs.long$type, ordered = TRUE, levels = 
                                c("linear", "both", "piecewise"))
# xyplot for predictive accuracy on test data:
correlation.anova <- aov(correlation.values ~ correlation.ind + N + rho + np + type +  
                          sigmab + corUb + correlation.ind*(N + rho + np + type +  
                            sigmab + corUb), data=treespecs.long)
summary(correlation.anova)[[1]]["Sum Sq"] / sum(summary(correlation.anova)[[1]]["Sum Sq"] )
#N, type & np have main and/or interaction effects with eta squared > .01
aggdata.cor <- aggregate(formula=correlation.values ~ correlation.ind + np + type + N, FUN=mean, 
                     data=treespecs.long)
xyplot(correlation.values ~ type | N + np, data = aggdata.cor, groups=correlation.ind, type="b",
       ylab="correlation", par.settings=standard.theme("pdf",color=F), 
       auto.key=list(space="top", columns=2, title=" ", cex.title=1,lines=T, points=T))


# xyplot for accuracy of variable recovery:
treespecs.long$N <- as.numeric(as.character(treespecs.long$N))
treespecs.long$sigmab <- as.numeric(as.character(treespecs.long$sigmab))
treespecs.long$type <- factor(treespecs.long$type, ordered = F, levels = 
                                c("linear", "both", "piecewise"))
treeacc.glm <- glm(truevars ~ truevars.ind + N + rho + np + type + sigmab + corUb +
                     truevars.ind*(N + rho + np + type + sigmab + corUb), 
                   data=treespecs.long, family="binomial")
summary(treeacc.glm)
#N, np & type have main and/or interaction effects
treespecs.long$N <- factor(as.numeric(as.character(treespecs.long$N)), ordered=T)
treespecs.long$type <- factor(treespecs.long$type, ordered = TRUE, levels = 
                                c("linear", "both", "piecewise"))
aggdata.acc <- aggregate(formula=truevars ~ truevars.ind + N + type + np, FUN=mean, 
                     data=treespecs.long)
xyplot(truevars ~ type | N + np, data = aggdata.acc, groups=truevars.ind, type="b",
       ylab="tree accuracy", par.settings=standard.theme("pdf",color=F), 
       auto.key=list(space="top", columns=2, title=" ", cex.title=1,lines=T, points=T))

