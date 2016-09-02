# install.packages("glmertree", repos="http://R-Forge.R-project.org")
library(glmertree)


## Evaluate tree size:

GLMMtreesize <- vector()
GLMtreesize <- vector()
for (c in 1:50){
  load(paste("treesize.GLMM", c, sep=""))
  GLMMtreesize <- c(GLMMtreesize, unlist(treesize.GLMM))
  load(paste("treesize.GLM", c, sep=""))
  GLMtreesize <- c(GLMtreesize, unlist(treesize.GLM))  
}
treespecs <- data.frame(GLMMtreesize, GLMtreesize)

load("descriptions1")
for (i in 1:length(descriptions)) {
  treespecs[i,"N"] <- substr(descriptions[[i]][[1]], 5, 8)
  treespecs[i,"rho"] <- substr(descriptions[[i]][[2]], 7, 9)
  treespecs[i,"np"] <- substr(descriptions[[i]][[3]], 24, 25)
  treespecs[i,"treatdiff"] <- substr(descriptions[[i]][[4]], 31, 33)
  treespecs[i,"corUbi"] <- substr(descriptions[[i]][[5]], 32, 64)
  treespecs[i,"numbclus"] <- substr(descriptions[[i]][[6]], 34, 35)
  treespecs[i,"sigmabi"] <- substr(descriptions[[i]][[7]], 12, 14)  
}

treespecs$N <- rep(factor(treespecs$N[1:length(descriptions)]), times=50)
treespecs$rho <- rep(factor(treespecs$rho[1:length(descriptions)]), times=50)
treespecs$np <- rep(factor(treespecs$np[1:length(descriptions)]), times=50)
treespecs$treatdiff <- rep(factor(treespecs$treatdiff[1:length(descriptions)]), times=50)
treespecs$corUbi <- rep(factor(treespecs$corUbi[1:length(descriptions)]), times=50)
treespecs$numbclus <- rep(factor(treespecs$numbclus[1:length(descriptions)]), times=50)
treespecs$sigmabi <- rep(factor(treespecs$sigmabi[1:length(descriptions)]), times=50)

table(treespecs$GLMMtreesize)
prop.table(table(treespecs$GLMMtreesize))
mean(treespecs$GLMMtreesize)
sd(treespecs$GLMMtreesize)

table(treespecs$GLMtreesize)
prop.table(table(treespecs$GLMtreesize))
mean(treespecs$GLMtreesize)
sd(treespecs$GLMtreesize)

treespecs.long <- data.frame(treesize=rbind(stack(treespecs[c("GLMMtreesize","GLMtreesize")]),
                                            data.frame(values = rep(NA, times = nrow(treespecs)), 
                                                       ind = rep("GLMMsize", times = nrow(treespecs)))), 
  N=rep(treespecs$N, 3), rho=rep(treespecs$rho, 3), np=rep(treespecs$np, 3), 
  treatdiff=rep(treespecs$treatdiff, 3),  corUb=rep(treespecs$corUbi, 3), 
  numbclus=rep(treespecs$numbclus, 3), sigmab=rep(treespecs$sigmabi, 3),
  datasetID=factor(rep(1:nrow(treespecs), 3)))
tmp <- as.character(treespecs.long$treesize.ind)
tmp[treespecs.long$treesize.ind=="GLMMtreesize"] <- "GLMM tree"
tmp[treespecs.long$treesize.ind=="GLMtreesize"] <- "GLM tree"
treespecs.long$treesize.ind <- factor(tmp)



## Evaluate predictive accuracy:

treespecs$cor_GLMMtree_true <- NA
treespecs$cor_GLMtree_true <- NA
treespecs$cor_GLMM_true <- NA
treespecs$mean_true_diffs <- NA
treespecs$mean_GLMMtree_diffs <- NA
treespecs$mean_GLMtree_diffs <- NA
treespecs$mean_GLMM_diffs <- NA
treespecs$sd_true_diffs <- NA
treespecs$sd_GLMMtree_diffs <- NA
treespecs$sd_GLMtree_diffs <- NA
treespecs$sd_GLMM_diffs <- NA

for(c in 1:50) {
  print(c)
  load(paste("treatdiffs", c, sep=""))
  #treespecs[(c-1)*length(descriptions)+(1:length(descriptions)),10:12] <- matrix(unlist(lapply(treatdiffs, cor)), ncol=9, byrow=T)[,c(2,3,6)]
  #treespecs[(c-1)*length(descriptions)+(1:length(descriptions)),13:15] <- t(sapply(treatdiffs, apply, 2, mean))
  #treespecs[(c-1)*length(descriptions)+(1:length(descriptions)),16:18] <- t(sapply(treatdiffs, apply, 2, sd))  
  for(i in 1:length(descriptions)) {
    treespecs[(c-1)*length(descriptions)+i, "cor_GLMMtree_true"] <- 
      with(treatdiffs[[i]], cor(GLMMtree_d_hat, true_d_hat))
    treespecs[(c-1)*length(descriptions)+i, "cor_GLMtree_true"] <-   
      with(treatdiffs[[i]], cor(GLMtree_d_hat, true_d_hat))
    treespecs[(c-1)*length(descriptions)+i, "cor_GLMM_true"] <- 
      with(treatdiffs[[i]], cor(GLMM_d_hat, true_d_hat))
  
    treespecs[(c-1)*length(descriptions)+i, "mean_true_diffs"] <-
      with(treatdiffs[[i]], mean(true_d_hat))
    treespecs[(c-1)*length(descriptions)+i, "mean_GLMMtree_diffs"] <- 
      with(treatdiffs[[i]], mean(GLMMtree_d_hat))
    treespecs[(c-1)*length(descriptions)+i, "mean_GLMtree_diffs"] <- 
      with(treatdiffs[[i]], mean(GLMtree_d_hat))
    treespecs[(c-1)*length(descriptions)+i, "mean_GLMM_diffs"] <- 
      with(treatdiffs[[i]], mean(GLMM_d_hat))

    treespecs[(c-1)*length(descriptions)+i, "sd_true_diffs"] <-
      with(treatdiffs[[i]], sd(true_d_hat))
    treespecs[(c-1)*length(descriptions)+i, "sd_GLMMtree_diffs"] <- 
      with(treatdiffs[[i]], sd(GLMMtree_d_hat))
    treespecs[(c-1)*length(descriptions)+i, "sd_GLMtree_diffs"] <- 
      with(treatdiffs[[i]], sd(GLMtree_d_hat))
    treespecs[(c-1)*length(descriptions)+i, "sd_GLMM_diffs"] <- 
      with(treatdiffs[[i]], sd(GLMM_d_hat))
  }
}

apply(treespecs[,10:20], 2, mean)
apply(treespecs[,10:20], 2, sd)

treespecs.long <- data.frame(correlation = stack(treespecs[
  c("cor_GLMMtree_true","cor_GLMtree_true", "cor_GLMM_true")]), treespecs.long)
tmp <- as.character(treespecs.long$correlation.ind)
tmp[treespecs.long$correlation.ind=="cor_GLMMtree_true"] <- "GLMM tree"
tmp[treespecs.long$correlation.ind=="cor_GLMtree_true"] <- "GLM tree"
tmp[treespecs.long$correlation.ind=="cor_GLMM_true"] <- "GLMM"
treespecs.long$correlation.ind <- factor(tmp)



## Evaluate tree accuracy:

# Recover splitting variables and split points:

tmp <- data.frame()
for (c in 1:50){
  load(paste("splits.GLM", c, sep=""))
  tmp <- data.frame(rbind(tmp, splits.GLM))
}
GLMsplits <- data.frame(GLMsplitvar=tmp[,1], GLMsplitval=tmp[,2])

tmp <- data.frame()
for (c in 1:50){
  load(paste("splits.GLMM",c,sep=""))
  tmp <- data.frame(rbind(tmp, splits.GLMM))
}
GLMMsplits <- data.frame(GLMMsplitvar=tmp[,1], GLMMsplitval=tmp[,2])

tmp <- vector() 
tmpnodeno <- vector()
for (j in 1:nrow(treespecs)) {
  print(j)
  tmp <- c(tmp, rep(treespecs$GLMtreesize[j], times=(treespecs$GLMtreesize[j]-1)/2))
  tmpnodeno <- c(tmpnodeno, seq(1,(treespecs$GLMtreesize[j]-1)/2))
} 
GLMsplits$treesize <- tmp 
GLMsplits$nodenumber <- tmpnodeno

tmp <- vector() 
tmpnodeno <- vector()
for (j in 1:nrow(treespecs)) {
  print(j)
  tmp <- c(tmp, rep(treespecs$GLMMtreesize[j], times=(treespecs$GLMMtreesize[j]-1)/2))
  tmpnodeno <- c(tmpnodeno, seq(1,(treespecs$GLMMtreesize[j]-1)/2))
} 
GLMMsplits$treesize <- tmp
GLMMsplits$nodenumber <- tmpnodeno

GLMMsplits$GLMMsplitvar <- factor(GLMMsplits$GLMMsplitvar)
GLMsplits$mobsplitvar <- factor(GLMsplits$GLMsplitvar)
save(GLMMsplits, file="GLMMsplits")
save(GLMMsplits, file="GLMMsplits")
load("GLMMsplits"); load("GLMMsplits")



# Assess accuracy of splits:

# Accuracy of first split:
table(GLMMsplits[GLMMsplits$nodenumber==1, "GLMMsplitvar"])
mean(GLMMsplits[GLMMsplits$nodenumber==1, "GLMMsplitval"])
sd(GLMMsplits[GLMMsplits$nodenumber==1, "GLMMsplitval"])
table(GLMsplits[GLMsplits$nodenumber==1, "GLMsplitvar"])
mean(GLMsplits[GLMsplits$nodenumber==1 & GLMsplits$GLMsplitvar==4, "GLMsplitval"])
sd(GLMsplits[GLMsplits$nodenumber==1 & GLMsplits$GLMsplitvar==4, "GLMsplitval"])

# identify which split belongs to which tree:
length(GLMMtrees)
648*50
GLMtreeno <- rep(1:32400, times=(treespecs$GLMtreesize-1)/2)
GLMMtreeno <- rep(1:32400, times=(treespecs$GLMMtreesize-1)/2)
GLMMsplits$treeno <- GLMMtreeno
GLMsplits$treeno <- GLMtreeno



# Check which trees were accurately recovered: 
# 1) true number of splits shoulde be recovered
# 2) true splitting variables should be recovered 
# 3) true splitting value plus or minus 5 = .5*\sigma should be recovered

# Find correct GLMM trees:
GLMMtruefirstsplit <- GLMMsplits[GLMMsplits$treesize==7 & GLMMsplits$GLMMsplitvar==4 & 
                                   GLMMsplits$nodenumber == 1 &
                                   GLMMsplits$GLMMsplitval>25 & GLMMsplits$GLMMsplitval<35,"treeno"]
GLMMtruesecsplit <- GLMMsplits[GLMMsplits$treesize==7 & GLMMsplits$GLMMsplitvar==3 & 
                                 GLMMsplits$nodenumber == 2 &
                                     GLMMsplits$GLMMsplitval>12 & GLMMsplits$GLMMsplitval<22,"treeno"]
GLMMtruethirdsplit <- GLMMsplits[GLMMsplits$treesize==7 & GLMMsplits$GLMMsplitvar==7 & 
                                   GLMMsplits$nodenumber == 3 &
                                       GLMMsplits$GLMMsplitval>58 & GLMMsplits$GLMMsplitval<68,"treeno"]
GLMMtruetreenos <- GLMMtruefirstsplit[GLMMtruefirstsplit %in% GLMMtruesecsplit][
      GLMMtruefirstsplit[GLMMtruefirstsplit %in% GLMMtruesecsplit] %in% GLMMtruethirdsplit]
GLMMsplits$truetree <- GLMMsplits$treeno %in% GLMMtruetreenos

# Find correct GLM trees:
GLMtruefirstsplit <- GLMsplits[GLMsplits$treesize==7 & GLMsplits$GLMsplitvar==4 & 
                                 GLMsplits$nodenumber == 1 &
                                 GLMsplits$GLMsplitval>25 & GLMsplits$GLMsplitval<35,"treeno"]
GLMtruesecsplit <- GLMsplits[GLMsplits$treesize==7 & GLMsplits$GLMsplitvar==3 & 
                               GLMsplits$nodenumber == 2 &
                               GLMsplits$GLMsplitval>12 & GLMsplits$GLMsplitval<22,"treeno"]
GLMtruethirdsplit <- GLMsplits[GLMsplits$treesize==7 & GLMsplits$GLMsplitvar==7 & 
                                 GLMsplits$nodenumber == 3 &
                                 GLMsplits$GLMsplitval>58 & GLMsplits$GLMsplitval<68,"treeno"]
GLMtruetreenos <- GLMtruefirstsplit[GLMtruefirstsplit %in% GLMtruesecsplit][
  GLMtruefirstsplit[GLMtruefirstsplit %in% GLMtruesecsplit] %in% GLMtruethirdsplit]
GLMsplits$truetree <- GLMsplits$treeno %in% GLMtruetreenos



## Assess GLMM accuracy:

# there should be a significant effect of X2:X1, X2:X5, T:X2:X1, T:X2:X5 
treespecs$trueGLMM <- NA
for(c in 1:50) {
  load(paste("GLMMs", c, sep = ""))
  print(c)
  for(i in 1:length(descriptions)) {
    treespecs$trueGLMM[(c-1)*length(descriptions) + i] <- 
      sum(c("X2", "X1:X2", "X2:X5", "T:X1:X2", "T:X2:X5") %in% 
            rownames(summary(GLMMs[[i]])$coefficients[abs(summary(
              GLMMs[[i]])$coefficients[,"t value"])>=1.96,])) == 5
  }
}

# Put results together in dataframe and save it:
treespecs$treeno <- 1:32400
treespecs$trueGLMMtree <- treespecs$treeno %in% unique(GLMMsplits$treeno[GLMMsplits$truetree])
treespecs$trueGLMtree <- treespecs$treeno %in% unique(GLMsplits$treeno[GLMsplits$truetree])
truetree <- stack(treespecs[c("trueGLMMtree","trueGLMtree", "trueGLMM")])

treespecs.long <- data.frame(truetree, treespecs.long)
treespecs.long$ind[treespecs.long$ind=="trueGLMMtree"] <- "GLMM tree"
treespecs.long$ind[treespecs.long$ind=="trueGLMtree"] <- "GLM tree"
treespecs.long$ind[treespecs.long$ind=="trueGLMM"] <- "GLMM"
treespecs.long$truetree.values <- factor(treespecs.long$values)
prop.table(table(treespecs$trueGLMtree))
prop.table(table(treespecs$trueGLMMtree))
prop.table(table(treespecs$trueGLMM))
save("treespecs.long", file="treespecs_long.dat")
save("treespecs", file="treespecs.dat")




## Create plots of outcomes:
library(lattice)
load("treespecs_long.dat")

# xyplot for treesize:
treespecs.long2 <- treespecs.long[treespecs.long$treesize.ind!="GLMMsize",]
treespecs.long2$treesize.ind <- factor(treespecs.long2$treesize.ind)

treesize.anova <- aov(treesize.values ~ treesize.ind + N + rho + np + treatdiff + numbclus + sigmab + 
                        corUb + treesize.ind*(N + rho + np + treatdiff + numbclus + sigmab + corUb), 
                      data=treespecs.long2)
summary(treesize.anova)[[1]]["Sum Sq"] / sum(summary(treesize.anova)[[1]]["Sum Sq"] )
#N, sigmab & corUbi have main and interaction effects with eta squared > .01
treespecs.long2$N <- factor(as.numeric(as.character(treespecs.long2$N)), ordered=T)
treespecs.long2$sigmab <- factor(as.numeric(as.character(treespecs.long2$sigmab)), ordered=T)
aggdata.size <- aggregate(formula=treesize.values ~ treesize.ind + N + sigmab + corUb, FUN=mean, 
                     data=treespecs.long2)
levels(aggdata.size$corUb)[levels(aggdata.size$corUb)=="bi and splitting U correlated"] <- "b correlated with splitting U"
levels(aggdata.size$corUb)[levels(aggdata.size$corUb)=="bi and non-splitting U correlated"] <- "b correlated with non-splitting U"
levels(aggdata.size$corUb)[levels(aggdata.size$corUb)=="uncorrelated"] <- "b and U uncorrelated"
xyplot(treesize.values ~ sigmab | N + corUb, data = aggdata.size, groups=treesize.ind, type="b", 
       ylab="tree size", xlab="sigma_b", par.settings=standard.theme("pdf",color=F), abline=c(7,0),
       auto.key=list(space="top", columns=2, title=" ", cex.title=1,lines=T, points=T))

# xyplot for predictive accuracy on test data:
correlation.anova <- aov(correlation.values ~ correlation.ind + N + rho + np + treatdiff + numbclus + 
                          sigmab + corUb + correlation.ind*(N + rho + np + treatdiff + numbclus + 
                            sigmab + corUb), data=treespecs.long2)
summary(correlation.anova)[[1]]["Sum Sq"] / sum(summary(correlation.anova)[[1]]["Sum Sq"] )
#N, treatdiffs & sigmab have main and/or interaction effects with eta squared > .01
treespecs.long2$treatdiff <- factor(as.numeric(as.character(treespecs.long2$treatdiff)), ordered=T)
treespecs.long2$correlation.ind <- factor(treespecs.long2$correlation.ind)
aggdata.cor <- aggregate(formula=correlation.values ~ correlation.ind + N + sigmab + treatdiff, FUN=mean, 
                     data=treespecs.long2)
xyplot(correlation.values ~ sigmab | N + treatdiff, data = aggdata.cor, groups=correlation.ind, type="b",
       ylab="correlation", xlab="sigma_b", par.settings=standard.theme("pdf",color=F), 
       auto.key=list(space="top", columns=2, title=" ", cex.title=1,lines=T, points=T))

# xyplot for accuracy of tree recovery:
treespecs.long2$N <- factor(treespecs.long2$N, ordered = FALSE)
treespecs.long2$sigmab <- factor(treespecs.long2$sigmab, ordered = FALSE)
treeacc.glm <- glm(truetree.values+1 ~ correlation.ind + N + rho + np + treatdiff + numbclus + sigmab + corUb +
                     correlation.ind*(N + rho + np + treatdiff + numbclus + sigmab + corUb), 
                   data=treespecs.long2, family="binomial")
summary(treeacc.glm)
#N, corUb & sigmab have strongest main and/or interaction effects
treespecs.long2$truetree.values <- as.numeric(treespecs.long2$truetree.values)-1
aggdata.acc <- aggregate(formula=truetree.values ~ correlation.ind + N + sigmab + corUb, FUN=mean, 
                     data=treespecs.long2)
levels(aggdata.acc$corUb)[levels(aggdata.acc$corUb)=="bi and splitting U correlated"] <- "b correlated with splitting U"
levels(aggdata.acc$corUb)[levels(aggdata.acc$corUb)=="bi and non-splitting U correlated"] <- "b correlated with non-splitting U"
levels(aggdata.acc$corUb)[levels(aggdata.acc$corUb)=="uncorrelated"] <- "b and U uncorrelated"
xyplot(truetree.values ~ sigmab | N + corUb, data = aggdata.acc, groups=correlation.ind, type="b",
       ylab="tree accuracy", xlab="sigma_b", par.settings=standard.theme("pdf",color=F), 
       auto.key=list(space="top", columns=2, title=" ", cex.title=1,lines=T, points=T))

