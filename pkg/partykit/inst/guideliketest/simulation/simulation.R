#################################
## Replication code for the simulation results presented in:
## The Power of Unbiased Recursive Partitioning: A Unifying View of CTree, MOB, and GUIDE (2019)
## by Lisa Schlosser and Torsten Hothorn and and Achim Zeileis
## (corresponding figures can be reproduced by function calls provided in figures.R)
#################################

library("partykit")
library("Formula")
library("parallel")

# load source files
source("guidelike.R")
source("tests_guideprune.R")
source("ccprune.R")



#############
### reproduction of simulation results



### "stump" scenario as discussed in Section 5.1. and in the Appendix
sim_stump <- simwrapper(scenario = "stump",
                        nobs = 250, nrep = 100, seed = 7,
                        delta = seq(from = 0, to = 1, by = 0.1),
                        xi = c(0, 0.2, 0.5, 0.8), 
                        variation = c("all", "beta0", "beta1"),
                        changetype = c("abrupt", "continuous"),
                        test = c("ctree", "mfluc", "ctree_max",
                                 "ctree_cat", "mfluc_cat", "ctree_max_cat",
                                 "ctree_bin", "mfluc_bin", "ctree_max_bin",
                                 "ctree_cat_bin", "mfluc_cat_bin", "ctree_max_cat_bin",
                                 "guide_sum_12", "guide_coin_12", "guide_sum_1_cor"))

sim_stump <- sim_stump$res




### "stump" scenario, full factorial analysis as discussed in Section 5.2.

# for delta = 0.3 and xi = 0
sim3way1 <- sim(scenario = "stump",
                nobs = 250, nrep = 100, seed = 7,
                delta = 0.3,
                xi = 0, 
                changetype = "abrupt",
                variation = "all", 
                return_matrices = TRUE,
                test = c("guide_sum_1_cor",
                         "mfluc_resid_bin",
                         "ctree_resid_bin",
                         "mfluc_resid_cat",
                         "mfluc_resid",
                         "ctree_resid",
                         "guide_sum_12",
                         "mfluc_bin",
                         "ctree_bin",
                         "mfluc_cat",
                         "mfluc",
                         "ctree"))
# prepare results for illustrations
sim_3way1 <- prep_3way(sim_3way1)


  
# for delta = 1 and xi = 0.8
sim_3way2 <- sim(scenario = "stump",
                 nobs = 250, nrep = 100, seed = 7,
                 delta = 1,
                 xi = 0.8, 
                 changetype = "abrupt",
                 variation = "all",
                 return_matrices = TRUE,
                 test = c("guide_sum_1_cor",
                          "mfluc_resid_bin",
                          "ctree_resid_bin",
                          "mfluc_resid_cat",
                          "mfluc_resid",
                          "ctree_resid",
                          "guide_sum_12",
                          "mfluc_bin",
                          "ctree_bin",
                          "mfluc_cat",
                          "mfluc",
                          "ctree"))
# prepare results for illustrations
sim_3way2 <- prep_3way(sim_3way2)  

sim_3way1$xi <- 0
sim_3way2$xi <- 0.8

# combine results for illustration in two panels as presented in Section 5.2. (Figure 4)
sim_3way <- rbind(sim_3way1, sim_3way2)
sim_3way$xi <- factor(sim_3way$xi, levels = c(0,0.8), labels = c(0,0.8))



### "tree" scenario as discussed in Section 5.3.
sim_tree <- simwrapper(scenario = "tree",
                       nobs = 250, nrep = 100, seed = 7, 
                       delta = seq(from = 0, to = 1, by = 0.1),
                       xi = c(0, 0.2, 0.5, 0.8), 
                       variation = "all", 
                       changetype = "abrupt",
                       compare_pruning = TRUE,
                       test = c("ctree", "mfluc",
                                "guide_sum_12", 
                                "guide_sum_1_cor"))

# prepare data for illustration as presented in Section 5.3. (Figure 5)
sim_tree <- sim_tree$res
# extract data for xyplot comparing ARI and ARI_p
sim_tree <- sim_tree[,c(1,2,5,6,7,10,13)]
# extended data set with post-pruning results being represented as separated tests
sim_tree_ex <- sim_tree[,c(1,2,3,4,6)]
sim_tree_ex <- rbind(sim_tree_ex, sim_tree_ex)

sim_tree_ex$test <- factor(sim_tree_ex$test,
                           levels = c("ctree","mfluc",
                                      "guide_sum_12", "guide_sum_1_cor"),
                           labels = c("CTree","MOB",
                                      "GUIDE+scores", "GUIDE"))

sim_tree_ex$nrsubgr[(length(sim_tree$nrsubgr)+1) : (2*length(sim_tree$nrsubgr))] <- sim_tree$nrsubgr_p
sim_tree_ex$ari[(length(sim_tree$ari)+1) : (2*length(sim_tree$ari))] <- sim_tree$ari_p
sim_tree_ex$pruning <- c(rep("pre", NROW(sim_tree)), rep("post", NROW(sim_tree)))
sim_tree <- sim_tree_ex
rm(sim_tree_ex)