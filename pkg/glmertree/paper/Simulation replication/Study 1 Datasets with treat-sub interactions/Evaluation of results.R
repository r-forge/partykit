# install.packages("glmertree", repos="http://R-Forge.R-project.org")
library(glmertree)

## Evaluate tree size:
GLMMtreesize <- vector()
GLMtreesize <- vector()
for (c in 1:50){
  load(paste("treesize.GLMM", c, sep = ""))
  GLMMtreesize <- c(GLMMtreesize, unlist(treesize.GLMM))
  load(paste("treesize.GLM", c, sep = ""))
  GLMtreesize <- c(GLMtreesize, unlist(treesize.GLM))
}
treespecs <- data.frame(GLMMtreesize, GLMtreesize)

## Create a wide dataframe for collecting all simulation results:
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
treespecs$N <- rep(factor(treespecs$N[1:length(descriptions)]), times = 50)
treespecs$rho <- rep(factor(treespecs$rho[1:length(descriptions)]), times = 50)
treespecs$np <- rep(factor(treespecs$np[1:length(descriptions)]), times = 50)
treespecs$treatdiff <- rep(factor(treespecs$treatdiff[1:length(descriptions)]), times = 50)
treespecs$corUbi <- rep(factor(treespecs$corUbi[1:length(descriptions)]), times = 50)
treespecs$numbclus <- rep(factor(treespecs$numbclus[1:length(descriptions)]), times = 50)
treespecs$sigmabi <- rep(factor(treespecs$sigmabi[1:length(descriptions)]), times = 50)

## Get some summary statistics about tree size:
table(treespecs$GLMMtreesize)
prop.table(table(treespecs$GLMMtreesize))
mean(treespecs$GLMMtreesize)
sd(treespecs$GLMMtreesize)
table(treespecs$GLMtreesize)
prop.table(table(treespecs$GLMtreesize))
mean(treespecs$GLMtreesize)
sd(treespecs$GLMtreesize)

## Make a long dataframe for collecting all simulation results:
treespecs.long <- data.frame(
  treesize = rbind(stack(treespecs[c("GLMMtreesize","GLMtreesize")]),
                   data.frame(values = rep(NA, times = nrow(treespecs)), 
                              ind = rep("MERTsize", times = nrow(treespecs)))), 
  N = rep(treespecs$N, 3), rho = rep(treespecs$rho, 3), 
  np = rep(treespecs$np, 3), treatdiff = rep(treespecs$treatdiff, 3), 
  corUb = rep(treespecs$corUbi, 3), numbclus = rep(treespecs$numbclus, 3), 
  sigmab = rep(treespecs$sigmabi, 3), 
  datasetID = factor(rep(1:nrow(treespecs), 3)))
tmp <- as.character(treespecs.long$treesize.ind)
tmp[treespecs.long$treesize.ind == "GLMMtreesize"] <- "GLMM tree"
tmp[treespecs.long$treesize.ind == "GLMtreesize"] <- "GLM tree"
tmp[treespecs.long$treesize.ind == "MERTsize"] <- "MERT"
treespecs.long$treesize.ind <- factor(tmp)



## Evaluate predictive accuracy:
treespecs$cor_GLMMtree_true <- NA
treespecs$cor_GLMtree_true <- NA
treespecs$cor_MERT_true <- NA
treespecs$mean_true_diffs <- NA
treespecs$mean_GLMMtree_diffs <- NA
treespecs$mean_GLMtree_diffs <- NA
treespecs$mean_MERT_diffs <- NA
treespecs$sd_true_diffs <- NA
treespecs$sd_GLMMtree_diffs <- NA
treespecs$sd_GLMtree_diffs <- NA
treespecs$sd_MERT_diffs <- NA

for(c in 1:50) {
  print(c)
  load(paste("treatdiffs", c, sep = ""))
  for(i in 1:length(descriptions)) {
    treespecs[(c - 1) * length(descriptions) + i, "cor_GLMMtree_true"] <- 
      with(treatdiffs[[i]], cor(GLMMtree_d_hat, true_d_hat))
    treespecs[(c - 1) * length(descriptions) + i, "cor_GLMtree_true"] <-   
      with(treatdiffs[[i]], cor(GLMtree_d_hat, true_d_hat))
    treespecs[(c - 1) * length(descriptions) + i, "cor_MERT_true"] <- 
      with(treatdiffs[[i]], cor(MERT_d_hat, true_d_hat))
    treespecs[(c - 1) * length(descriptions) + i, "mean_true_diffs"] <-
      with(treatdiffs[[i]], mean(true_d_hat))
    treespecs[(c - 1) * length(descriptions) + i, "mean_GLMMtree_diffs"] <- 
      with(treatdiffs[[i]], mean(GLMMtree_d_hat))
    treespecs[(c - 1) * length(descriptions) + i, "mean_GLMtree_diffs"] <- 
      with(treatdiffs[[i]], mean(GLMtree_d_hat))
    treespecs[(c - 1) * length(descriptions) + i, "mean_MERT_diffs"] <- 
      with(treatdiffs[[i]], mean(MERT_d_hat))
    treespecs[(c - 1) * length(descriptions) + i, "sd_true_diffs"] <-
      with(treatdiffs[[i]], sd(true_d_hat))
    treespecs[(c - 1) * length(descriptions) + i, "sd_GLMMtree_diffs"] <- 
      with(treatdiffs[[i]], sd(GLMMtree_d_hat))
    treespecs[(c - 1) * length(descriptions) + i, "sd_GLMtree_diffs"] <- 
      with(treatdiffs[[i]], sd(GLMtree_d_hat))
    treespecs[(c - 1) * length(descriptions) + i, "sd_MERT_diffs"] <- 
      with(treatdiffs[[i]], sd(MERT_d_hat))
  }
}

## Get some summary statistics about predictive accuracy:
apply(treespecs[,10:20], 2, mean)
apply(treespecs[,10:20], 2, sd)

## Put predictive accuracy findings in long dataframe:
treespecs.long <- data.frame(treespecs.long, correlation = stack(
    treespecs[c("cor_GLMMtree_true","cor_GLMtree_true", "cor_MERT_true")]))
tmp <- as.character(treespecs.long$correlation.ind)
tmp[treespecs.long$correlation.ind=="cor_GLMMtree_true"] <- "GLMM tree"
tmp[treespecs.long$correlation.ind=="cor_GLMtree_true"] <- "GLM tree"
tmp[treespecs.long$correlation.ind=="cor_MERT_true"] <- "MERT"
treespecs.long$correlation.ind <- factor(tmp)



## Evaluate tree accuracy:

# Get splitting variables and split points:

GLMsplits <- matrix(ncol = 2, dimnames = list(NULL, c("GLMsplitvar", "GLMsplitval")))
for (c in 1:50){
  load(paste("splits.GLM", c, sep= ""))
  GLMsplits <- rbind(GLMsplits, matrix(unlist(splits.GLM), ncol = 2, byrow = TRUE))
}
GLMsplits <- data.frame(GLMsplits[-1,])

GLMMsplits <- matrix(ncol = 2, dimnames = list(NULL, c("GLMMsplitvar", "GLMMsplitval")))
for (c in 1:50){
  load(paste("splits.GLMM", c, sep= ""))
  GLMMsplits <- rbind(GLMMsplits, matrix(unlist(splits.GLMM), ncol = 2, byrow = TRUE))
}
GLMMsplits <- data.frame(GLMMsplits[-1,])

size <- vector() 
nodeno <- vector()
for (j in 1:nrow(treespecs)) {
  print(j)
  size <- c(size, rep(treespecs$GLMtreesize[j], times = (treespecs$GLMtreesize[j]-1) / 2))
  nodeno <- c(nodeno, seq(1, (treespecs$GLMtreesize[j]-1) / 2))
} 
GLMsplits$treesize <- size
GLMsplits$nodenumber <- nodeno

size <- vector() 
nodeno <- vector()
for (j in 1:nrow(treespecs)) {
  print(j)
  size <- c(size, rep(treespecs$GLMMtreesize[j], times = (treespecs$GLMMtreesize[j]-1) / 2))
  nodeno <- c(nodeno, seq(1, (treespecs$GLMMtreesize[j]-1) / 2))
} 
GLMMsplits$treesize <- size
GLMMsplits$nodenumber <- nodeno

GLMMsplits$GLMMsplitvar <- factor(GLMMsplits$GLMMsplitvar)
GLMsplits$GLMsplitvar <- factor(GLMsplits$GLMsplitvar)
save(GLMMsplits, file="GLMMsplits")
save(GLMMsplits, file="GLMMsplits")
load("GLMMsplits"); load("GLMMsplits")



## Assess split accuracy:

## Accuracy of first split:
table(GLMMsplits[GLMMsplits$nodenumber == 1, "GLMMsplitvar"])
mean(GLMMsplits[GLMMsplits$nodenumber == 1, "GLMMsplitval"])
sd(GLMMsplits[GLMMsplits$nodenumber == 1, "GLMMsplitval"])
table(GLMsplits[GLMsplits$nodenumber == 1, "GLMsplitvar"])
mean(GLMsplits[GLMsplits$nodenumber == 1 & GLMsplits$GLMsplitvar == 4, "GLMsplitval"])
sd(GLMsplits[GLMsplits$nodenumber == 1 & GLMsplits$GLMsplitvar == 4, "GLMsplitval"])

## identify which split belongs to which tree:
GLMsplits$treeno <- rep(1:32400, times = (treespecs$GLMtreesize-1) / 2)
GLMMsplits$treeno <- rep(1:32400, times = (treespecs$GLMMtreesize-1) / 2)




## Check which trees were accurately recovered: 
## 1) true number of splits shoulde be recovered
## 2) true splitting variables should be recovered 
## 3) true splitting value plus or minus 5 = .5*\sigma should be recovered

## Get accuracy of GLMM trees:
GLMMtruefirstsplit <- GLMMsplits[GLMMsplits$treesize == 7 & 
                                   GLMMsplits$GLMMsplitvar == 4 & GLMMsplits$nodenumber == 1 &
                                   GLMMsplits$GLMMsplitval > 25 & GLMMsplits$GLMMsplitval < 35, "treeno"]
GLMMtruesecsplit <- GLMMsplits[GLMMsplits$treesize == 7 & 
                                 GLMMsplits$GLMMsplitvar == 3 & GLMMsplits$nodenumber == 2 &
                                 GLMMsplits$GLMMsplitval > 12 & GLMMsplits$GLMMsplitval < 22, "treeno"]
GLMMtruethirdsplit <- GLMMsplits[GLMMsplits$treesize == 7 & 
                                   GLMMsplits$GLMMsplitvar == 7 & GLMMsplits$nodenumber == 3 &
                                   GLMMsplits$GLMMsplitval > 58 & GLMMsplits$GLMMsplitval < 68, "treeno"]
GLMMtruetreenos <- GLMMtruefirstsplit[GLMMtruefirstsplit %in% GLMMtruesecsplit][
  GLMMtruefirstsplit[GLMMtruefirstsplit %in% GLMMtruesecsplit] %in% GLMMtruethirdsplit]
GLMMsplits$truetree <- GLMMsplits$treeno %in% GLMMtruetreenos

## Get accuracy of GLM trees:
GLMtruefirstsplit <- GLMsplits[GLMsplits$treesize == 7 & 
                                 GLMsplits$GLMsplitvar == 4 & GLMsplits$nodenumber == 1 &
                                 GLMsplits$GLMsplitval > 25 & GLMsplits$GLMsplitval < 35, "treeno"]
GLMtruesecsplit <- GLMsplits[GLMsplits$treesize == 7 & 
                               GLMsplits$GLMsplitvar == 3 & GLMsplits$nodenumber == 2 &
                               GLMsplits$GLMsplitval > 12 & GLMsplits$GLMsplitval < 22,"treeno"]
GLMtruethirdsplit <- GLMsplits[GLMsplits$treesize == 7 & 
                                 GLMsplits$GLMsplitvar == 7 & GLMsplits$nodenumber == 3 &
                                 GLMsplits$GLMsplitval > 58 & GLMsplits$GLMsplitval < 68,"treeno"]
GLMtruetreenos <- GLMtruefirstsplit[GLMtruefirstsplit %in% GLMtruesecsplit][
  GLMtruefirstsplit[GLMtruefirstsplit %in% GLMtruesecsplit] %in% GLMtruethirdsplit]
GLMsplits$truetree <- GLMsplits$treeno %in% GLMtruetreenos


## Add tree accuracy to results:
treespecs$treeno <- 1:32400
treespecs$trueGLMMtree <- treespecs$treeno %in% unique(GLMMsplits$treeno[GLMMsplits$truetree])
treespecs$trueGLMtree <- treespecs$treeno %in% unique(GLMsplits$treeno[GLMsplits$truetree])
truetree <- stack(cbind(treespecs[c("trueGLMMtree","trueGLMtree")], trueMERT = rep(NA, times = nrow(treespecs))))
truetree$ind <- as.character(truetree$ind)
truetree$ind[truetree$ind == "trueGLMMtree"] <- "GLMM tree"
truetree$ind[truetree$ind=="trueGLMtree"] <- "GLM tree"
truetree$ind[truetree$ind=="trueMERT"] <- "MERT"
names(truetree) <- c("truetree.values", "truetree.ind")

treespecs.long <- data.frame(truetree, treespecs.long)
treespecs.long$truetree.values <- factor(treespecs.long$truetree.values)
prop.table(table(treespecs$trueGLMtree))
prop.table(table(treespecs$trueGLMMtree))
save("treespecs.long", file = "treespecs_long.dat")
save("treespecs", file = "treespecs.dat")
save("treespecs.long", file = "treespecs_studyI.dat")



## Create plots of outcomes:
library(lattice)
load("treespecs_long.dat")

## Treesize:
treespecs.long2 <- treespecs.long[treespecs.long$treesize.ind != "MERT",]
treespecs.long2$treesize.ind <- factor(treespecs.long2$treesize.ind)
treesize.anova <- aov(treesize.values ~ treesize.ind + N + rho + np + treatdiff + numbclus + sigmab + 
                        corUb + treesize.ind*(N + rho + np + treatdiff + numbclus + sigmab + corUb), 
                      data=treespecs.long2)
summary(treesize.anova)[[1]]["Sum Sq"] / sum(summary(treesize.anova)[[1]]["Sum Sq"] )
#N, sigmab & corUbi have main and interaction effects with eta squared > .01
treespecs.long2$N <- factor(as.numeric(as.character(treespecs.long2$N)), ordered = TRUE)
treespecs.long2$sigmab <- factor(as.numeric(as.character(treespecs.long2$sigmab)), ordered = TRUE)
aggdata.size <- aggregate(formula=treesize.values ~ treesize.ind + N + sigmab + corUb, FUN = mean, 
                     data = treespecs.long2)
levels(aggdata.size$corUb)[levels(aggdata.size$corUb) == "bi and splitting U correlated"] <- "b correlated with splitting U"
levels(aggdata.size$corUb)[levels(aggdata.size$corUb) == "bi and non-splitting U correlated"] <- "b correlated with non-splitting U"
levels(aggdata.size$corUb)[levels(aggdata.size$corUb) == "uncorrelated"] <- "b and U uncorrelated"
xyplot(treesize.values ~ sigmab | N + corUb, data = aggdata.size, groups = treesize.ind, type = "b", 
       ylab = "tree size", xlab = "sigma_b", par.settings = standard.theme("pdf", color = FALSE), abline = c(7,0),
       auto.key = list(space = "top", columns = 2, title = " ", cex.title = 1, lines = TRUE, points = TRUE))



## Predictive accuracy on test data:
correlation.anova <- aov(correlation.values ~ 
                           correlation.ind + N + rho + np + treatdiff + 
                           numbclus + sigmab + corUb + 
                           correlation.ind*(N + rho + np + treatdiff + 
                                              numbclus + sigmab + corUb), 
                         data = treespecs.long)
summary(correlation.anova)[[1]]["Sum Sq"] / sum(summary(correlation.anova)[[1]]["Sum Sq"] )
#N, treatdiffs & sigmab have main and/or interaction effects with eta squared > .01
treespecs.long$treatdiff <- factor(as.numeric(as.character(treespecs.long$treatdiff)), ordered = TRUE)
treespecs.long$N <- factor(as.numeric(as.character(treespecs.long$N)), ordered = TRUE)
treespecs.long$sigmab <- factor(as.numeric(as.character(treespecs.long$sigmab)), ordered = TRUE)
treespecs.long$correlation.ind <- factor(treespecs.long$correlation.ind)
aggdata.cor <- aggregate(formula = correlation.values ~ correlation.ind + N + sigmab + treatdiff, FUN = mean, 
                     data=treespecs.long)
xyplot(correlation.values ~ sigmab | N + treatdiff, data = aggdata.cor, groups = correlation.ind, type = "b",
       ylab = "correlation", xlab = "sigma_b", par.settings = standard.theme("pdf", color = FALSE), 
       auto.key = list(space = "top", columns = 2, title = " ", cex.title = 1, lines = TRUE, points = TRUE))


## Tree accuracy:
treespecs.long2$N <- factor(treespecs.long2$N, ordered = FALSE)
treespecs.long2$sigmab <- factor(treespecs.long2$sigmab, ordered = FALSE)
treeacc.glm <- glm(truetree.values ~ truetree.ind + N + rho + np + treatdiff + numbclus + sigmab + corUb +
                     truetree.ind*(N + rho + np + treatdiff + numbclus + sigmab + corUb), 
                   data=treespecs.long2, family = "binomial")
summary(treeacc.glm)
#N, corUb & sigmab have strongest main and/or interaction effects
treespecs.long2$truetree.values <- as.numeric(treespecs.long2$truetree.values)-1
aggdata.acc <- aggregate(formula=truetree.values ~ truetree.ind + N + sigmab + corUb, FUN=mean, 
                     data=treespecs.long2)
levels(aggdata.acc$corUb)[levels(aggdata.acc$corUb) == "bi and splitting U correlated"] <- "b correlated with splitting U"
levels(aggdata.acc$corUb)[levels(aggdata.acc$corUb) == "bi and non-splitting U correlated"] <- "b correlated with non-splitting U"
levels(aggdata.acc$corUb)[levels(aggdata.acc$corUb) == "uncorrelated"] <- "b and U uncorrelated"
xyplot(truetree.values ~ sigmab | N + corUb, data = aggdata.acc, groups = truetree.ind, type = "b",
       ylab = "tree accuracy", xlab = "sigma_b", par.settings = standard.theme("pdf", color = FALSE), 
       auto.key = list(space = "top", columns = 2, title = " ", cex.title = 1, lines = TRUE, points = TRUE))