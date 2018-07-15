library("partykit")
library("Formula")
library("parallel")
library("lattice")
#setwd("~/svn/partykit/pkg/partykit/inst/guideliketest/draft/")
source("../guidelike.R")
source("../tests.R")


simres <- simwrapper_p(nobs = 250, nrep = 100, nrsteps = 1,
                       delta = seq(from = 0, to = 1, by = 0.25),
                       xi = c(0, 0.2, 0.5, 0.8), 
                       vary_beta = c("all", "beta0", "beta1"),
                       binary_regressor = c(TRUE, FALSE),
                       #binary_regressor = FALSE,
                       binary_beta = c(TRUE, FALSE),
                       only_intercept = c(TRUE, FALSE),
                       #only_intercept = FALSE,
                       test = c("ctree", "mfluc", "ctree_max",
                                "ctree_cat", "mfluc_cat","ctree_max_cat",
                                "ctree_bin", "mfluc_bin","ctree_max_bin",
                                "ctree_cat_bin", "mfluc_cat_bin","ctree_max_cat_bin",
                                "guide_sum_12", 
                                "guide_coin_12",
                                "guide_sum_1_cor",
                                "guide_sum_2_cor"),
                       beta0 = 0, beta1 = 1,
                       stump = TRUE, z1dist = "unif", sigma = 1, alpha = 0.05)


save(simres, file = "~/svn/partykit/pkg/partykit/inst/guideliketest/sim/simres20180709_1step.rda")



