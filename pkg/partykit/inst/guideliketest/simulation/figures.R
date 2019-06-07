#################################
## Replication code for the figures presented in:
## The Power of Unbiased Recursive Partitioning: A Unifying View of CTree, MOB, and GUIDE (2019)
## by Lisa Schlosser and Torsten Hothorn and and Achim Zeileis
#################################

library("lattice")
library("latticeExtra")
library("scales")
library("colorspace")

## HCL palette
pal <- qualitative_hcl(5, "Dark 3")
names(pal) <- c("CTree", "MOB", "GUIDE", "GUIDE+scores", "CTree+max")


# If existing, load .rda files containing the simulation results,
# otherwise call simulation.R to calculate the results.
# (computation time for full replication of simulation results by simulation.R: 
# approximately 70 hours, on our machines, using 3 kernels)

if(file.exists("sim_stump.rda") & file.exists("sim_3way.rda") & file.exists("sim_tree.rda")){
  load("sim_stump.rda")
  load("sim_3way.rda")
  load("sim_tree.rda")
} else {
  source("simulation.R")
  save(sim_stump, file = "sim_stump.rda")
  save(sim_3way, file = "sim_3way.rda")
  save(sim_tree, file = "sim_tree.rda")
}



#############
###  reproduction of figures


### Figure 3: "stump" scenario
subdata <- subset(sim_stump, test %in% c("ctree","mfluc",
                                         "guide_sum_12", "guide_sum_1_cor"))
subdata$test <- factor(subdata$test, 
                       levels = c("ctree","mfluc",
                                  "guide_sum_12", "guide_sum_1_cor"),
                       labels = c("CTree","MOB",
                                  "GUIDE+scores", "GUIDE"))

subdataxi <- subset(subdata, xi %in% c(0,0.8))

par.settings <- list(superpose.symbol = list(col = pal[c("CTree", "MOB", "GUIDE+scores", "GUIDE")], 
                                             fill = pal[c("CTree", "MOB", "GUIDE+scores", "GUIDE")]), 
                     superpose.line = list(col = pal[c("CTree", "MOB", "GUIDE+scores", "GUIDE")]))

useOuterStrips(xyplot(prop_Tsplit ~ delta | variation + xi, groups = ~ test, 
                      data = subdataxi, 
                      subset = changetype == "abrupt", 
                      type = "l", lty = c(1,1,1,1),
                      lwd = 2,
                      scales=list(x=list(at=seq(1,11,2))),
                      key =  list(text = list(levels(subdataxi$test)),
                                  lines = list(lty = c(1,1,1,1), type = "l", pch = 1, cex = 0.8,
                                               col = pal[c("CTree", "MOB", "GUIDE+scores", "GUIDE")], lwd = 2)),
                      layout= c(3,2),
                      index.cond = list(c(2,3,1),c(2,1)),
                      xlab = expression(delta),
                      ylab = expression("Selection probability of"~Z[1]),
                      par.settings = par.settings),
               
               strip =strip.custom(factor.levels = c(expression("varying"~beta[0]~"and"~beta[1]), 
                                                     expression("varying"~beta[0]), 
                                                     expression("varying"~beta[1])),
                                   bg = "gray90"),
               strip.left=strip.custom(factor.levels = c(expression(xi == 0~" (50%)"),                                                              expression(xi == 0.8~" (90%)" )),
                                       bg = "gray90"))



### Figure 4: "stump" scenario, full factorial analysis
xyplot(pval_z1 ~ cat | xi, 
       groups = interaction(bin, res_scores), 
       data = sim_3way, 
       ylim = c(0,0.5),
       type = c("a"), 
       col = c(pal[1], pal[3], pal[1], pal[3]),
       lty = c(1,1,2,2),
       lwd = 2,
       scales = list(alternating=FALSE),
       xlab = "Categorization",
       ylab = expression("p-value of"~Z[1]),
       strip = strip.custom(factor.levels = c(expression(xi==0~"(50%)"~" and"~delta==0.3), 
                                              expression(xi==0.8~"(90%)"~" and"~delta==1)),
                            bg = "gray90"),
       key = list(text = list(c("Goodness of fit:", "  residuals", "  scores", "Dichotomization:", "  yes", "  no")),
                  columns = 2,
                  lines = list(lty = c(0,1,2,0,1,1), 
                               type = "l", lwd = 2,
                               col = c("black","black","black","black", pal[1], pal[3]))))



### Figure 5: "tree" scenario
par.settings <- list(superpose.symbol = list(col = pal[c("CTree", "MOB", "GUIDE+scores", "GUIDE")], 
                                             fill = pal[c("CTree", "MOB", "GUIDE+scores", "GUIDE")]), 
                     superpose.line = list(col = pal[c("CTree", "MOB", "GUIDE+scores", "GUIDE")]))
useOuterStrips(xyplot(ari ~ delta | xi + pruning, groups = ~ test, 
                      data = sim_tree, 
                      type = "l", lty = c(1,1,1,1),
                      lwd = 2,
                      key =  list(text = list(levels(sim_tree$test)[1:4]),
                                  lines = list(lty = c(1,1,1,1), type = "l", pch = 1, cex = 0.8,
                                               col = pal[c("CTree", "MOB", "GUIDE+scores", "GUIDE")], lwd = 2)),
                      layout= c(4,2),
                      scales=list(x=list(at=seq(1,11,2))),
                      xlab = expression(delta),
                      ylab = "Adjusted Rand Index",
                      #index.cond = list(c(2,3,1),c(2,1)),
                      par.settings = par.settings),
               strip = strip.custom(factor.levels = c(expression(xi==0), 
                                                      expression(xi==0.2), 
                                                      expression(xi==0.5), 
                                                      expression(xi==0.8)),
                                    bg = "gray90"),
               strip.left=strip.custom(factor.levels = c("Post-pruning", "Pre-pruning"),
                                       bg = "gray90"))



### Appendix, Figure 6: "stump" scenario without significance level 
subdata <- subset(sim_stump, test %in% c("ctree", "mfluc", "guide_sum_12", "guide_sum_1_cor"))
subdata$test <- factor(subdata$test, 
                       levels = c("ctree", "mfluc",
                                  "guide_sum_12", "guide_sum_1_cor"),
                       labels = c("CTree", "MOB",
                                  "GUIDE+scores", "GUIDE"))

subdataxi <- subset(subdata, xi %in% c(0,0.8))

par.settings <- list(superpose.symbol = list(col = pal[c("CTree", "MOB", "GUIDE+scores", "GUIDE")], 
                                             fill = pal[c("CTree", "MOB", "GUIDE+scores", "GUIDE")]), 
                     superpose.line = list(col = pal[c("CTree", "MOB", "GUIDE+scores", "GUIDE")]))

library("latticeExtra")
useOuterStrips(xyplot(prop_z1 ~ delta | vary_beta + xi, groups = ~ test, 
                      data = subdataxi, 
                      subset = binary_beta == TRUE & binary_regressor == FALSE & 
                        only_intercept == FALSE,# & xi == 0, 
                      type = "l", lty = c(1,1,1,1),
                      lwd = 2,
                      scales=list(x=list(at=seq(1,11,2))),
                      key =  list(text = list(levels(subdataxi$test)),
                                  lines = list(lty = c(1,1,1,1), type = "l", pch = 1, cex = 0.8,
                                               col = pal[c("CTree", "MOB", "GUIDE+scores", "GUIDE")], lwd = 2)),
                      layout= c(3,2),
                      index.cond = list(c(2,3,1),c(2,1)),
                      xlab = expression(delta),
                      ylab = expression("Prop. of"~Z[1]~"having lowest p-value"),
                      par.settings = par.settings),
               
               strip =strip.custom(factor.levels = c(expression("varying"~beta[0]~"and"~beta[1]), 
                                                     expression("varying"~beta[0]), 
                                                     expression("varying"~beta[1])),
                                   bg = "gray90"),
               strip.left=strip.custom(factor.levels = c(expression(xi == 0~" (50%)"),                                                              
                                                         expression(xi == 0.8~" (90%)")),
                                       bg = "gray90"))



### Appendix, Figure 7: "stump" scenario for continuous parameter functions
subdata <- subset(sim_stump, test %in% c("ctree", "ctree_max","mfluc", "guide_sum_1_cor","guide_sum_12"))
subdata$test <- factor(subdata$test, 
                       levels = c("ctree", "ctree_max","mfluc", "guide_sum_1_cor","guide_sum_12"),
                       labels = c("CTree","CTree+max","MOB", "GUIDE","GUIDE+scores"))
par.settings <- list(superpose.symbol = list(col = pal[c("CTree","CTree+max","MOB", "GUIDE","GUIDE+scores")], 
                                             fill = pal[c("CTree","CTree+max","MOB", "GUIDE","GUIDE+scores")]), 
                     superpose.line = list(col = pal[c("CTree","CTree+max","MOB", "GUIDE","GUIDE+scores")]))

xyplot(prop_Tsplit ~ delta | vary_beta, groups = ~ test, 
       data = subdata, 
       subset = binary_beta == FALSE & binary_regressor == FALSE & 
         only_intercept == FALSE, 
       scales=list(x=list(at=seq(1,11,2))),
       type = "l", 
       layout= c(3,1),
       lwd = 2,
       key =  list(text = list(levels(subdata$test)),
                   lines = list(lty = c(1,1,1,1,1), 
                                type = "l", 
                                col = pal[c("CTree","CTree+max","MOB", "GUIDE","GUIDE+scores")], 
                                lwd = 2)), 
       xlab = expression(delta),
       ylab = expression("Selection probability of"~Z[1]),
       index.cond = list(c(2,3,1)),
       strip = strip.custom(factor.levels = c(expression("varying"~beta[0]~" and "~beta[1]),
                                              expression("varying"~beta[0]), 
                                              expression("varying"~beta[1])),
                            bg = "gray90"),
       par.settings = par.settings)



### Appendix, Figure 8: "stump" scenario with modified versions of CTree and MOB 
# (CTree+dich and MOB+dich applying a dichotomization of residuals/scores)
subdata <- subset(sim_stump, test %in% c("ctree", "ctree_bin", 
                                         "mfluc", "mfluc_bin",
                                         "guide_sum_12"))

subdata <- subset(subdata, delta %in% c(0, 0.1, 0.2, 0.3, 0.4, 0.5))
subdata$test <- factor(subdata$test, 
                       levels = c("ctree", "ctree_bin", 
                                  "mfluc", "mfluc_bin",
                                  "guide_sum_12"),
                       labels = c("CTree","CTree+dich",
                                  "MOB", "MOB+dich",
                                  "GUIDE+scores"))

par.settings <- list(superpose.symbol = list(col = pal[c("CTree", "CTree", 
                                                         "MOB", "MOB",
                                                         "GUIDE+scores")], 
                                             fill = pal[c("CTree", "CTree", 
                                                          "MOB", "MOB",
                                                          "GUIDE+scores")]), 
                     superpose.line = list(col = pal[c("CTree", "CTree", 
                                                       "MOB", "MOB",
                                                       "GUIDE+scores")]))

xyplot(prop_Tsplit ~ delta | xi, groups = ~ test, 
       data = subdata, 
       subset = binary_beta == TRUE & binary_regressor == FALSE & 
         only_intercept == FALSE & vary_beta == "all", 
       type = "l",
       lty = c(1,2,#1,2,
               1,2,1),
       lwd = 2,
       layout= c(4,1),
       key =  list(text = list(levels(subdata$test)),
                   lines = list(lty = c(1,2,#1,2,
                                        1,2,1), 
                                type = "l", pch = 1, cex = 0.8,
                                col = pal[c("CTree", "CTree", 
                                            "MOB", "MOB",
                                            "GUIDE+scores")], lwd = 2)), 
       xlab = expression(delta),
       ylab = expression("Selection probability of"~Z[1]),
       strip = strip.custom(factor.levels = c(expression(xi==0), 
                                              expression(xi==0.2), 
                                              expression(xi==0.5), 
                                              expression(xi==0.8)),
                            bg = "gray90"),
       par.settings = par.settings)



### Appendix, Figure 9: "stump" scenario with modified versions of CTree and MOB 
# (CTree+cat and MOB+cat applying a categorization of split variables)
subdata <- subset(sim_stump, test %in% c("ctree", "ctree_cat",
                                         "mfluc", "mfluc_cat",
                                         "guide_sum_12"))
subdata <- subset(subdata, delta %in% c(0, 0.1, 0.2, 0.3, 0.4, 0.5))
subdata$test <- factor(subdata$test, 
                       levels = c("ctree", "ctree_cat", 
                                  "mfluc", "mfluc_cat",
                                  "guide_sum_12"),
                       labels = c("CTree","CTree+cat",
                                  "MOB", "MOB+cat",
                                  "GUIDE+scores"))

par.settings <- list(superpose.symbol = list(col = pal[c("CTree", "CTree", 
                                                         "MOB", "MOB",
                                                         "GUIDE+scores")], 
                                             fill = pal[c("CTree", "CTree", 
                                                          "MOB", "MOB",
                                                          "GUIDE+scores")]), 
                     superpose.line = list(col = pal[c("CTree", "CTree", 
                                                       "MOB", "MOB",
                                                       "GUIDE+scores")]))

xyplot(prop_Tsplit ~ delta | xi, groups = ~ test, 
       data = subdata, 
       subset = binary_beta == TRUE & binary_regressor == FALSE & 
         only_intercept == FALSE & vary_beta == "all", 
       type = "l", 
       lty = c(1,2,#1,2,
               1,2,1),
       lwd = 2,
       layout= c(4,1),
       key =  list(text = list(levels(subdata$test)),
                   lines = list(lty = c(1,2,#1,2,
                                        1,2,1), 
                                type = "l", pch = 1, cex = 0.8,
                                col = pal[c("CTree", "CTree", 
                                            "MOB", "MOB",
                                            "GUIDE+scores")], lwd = 2)), 
       xlab = expression(delta),
       ylab = expression("Selection probability of"~Z[1]),
       strip = strip.custom(factor.levels = c(expression(xi==0), 
                                              expression(xi==0.2), 
                                              expression(xi==0.5), 
                                              expression(xi==0.8)),
                            bg = "gray90"),
       par.settings = par.settings)
