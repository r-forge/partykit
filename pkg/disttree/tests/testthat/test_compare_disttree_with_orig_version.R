# -------------------------------------------------------------------
# - NAME:   test_compare_disttree_with_orig_version.R 
# - AUTHOR: Moritz N. Lang, Lisa Schlosser
# - DATE:   2019-08-12
# -------------------------------------------------------------------
# - PURPOSE: Test new implementation (build on extree()) with original implemention 
#            build on mob() and ctree()
# -------------------------------------------------------------------

context("comparison of new implementation 'disttree()' with original version")


# -------------------------------------------------------------------
# Cars example
# -------------------------------------------------------------------
suppressWarnings(RNGversion("3.5.0"))
set.seed(123)
### original disttree implementation with old control arguments (
#tr_cars <- disttree(dist ~ speed, data = cars, type.tree = "ctree")

## new disttree implementation with new control arguments (
tr_cars <- distextree(dist ~ speed, data = cars, type.tree = "ctree", 
  control = distextree_control(minsplit = 20L, minbucket = 7L, splittry = 2L))

expect_known_value(apply(coef(tr_cars), 2, sort), file = "test_compare_disttree_with_orig_version_tr_cars_coef.rds", update = FALSE, tolerance = 1e-6)
expect_known_value(logLik(tr_cars), file = "test_compare_disttree_with_orig_version_tr_cars_logLik.rds", update = FALSE)
expect_known_value(predict(tr_cars), file = "test_compare_disttree_with_orig_version_tr_cars_pred.rds", update = FALSE, tolerance = 1e-6)
expect_known_output(print(as.constparty(tr_cars)), file = "test_compare_disttree_with_orig_version_tr_cars_constparty.txt", update = FALSE, tolerance = 1e-6)


# -------------------------------------------------------------------
# Demo: application to precipitation data at Axams
# -------------------------------------------------------------------

# further packages
foo_package <- function(x, repos = getOption("repos")) {
  if (!require(x, character.only = TRUE)) {
     install.packages(x, dependencies = TRUE, repos = repos)
     library(x, character.only = TRUE)
  }
}

invisible(lapply(c("crch", "scoringRules", "gamlss.cens"), function(x) foo_package(x)))
invisible(foo_package("RainTyrol", repos = "http://r-forge.r-project.org"))

# if gamlss.cens family object should be used as family
gen.cens(NO, type = "left")

assign("NO",  gamlss.dist::NO,  pos = ".GlobalEnv")
assign("dNO", gamlss.dist::dNO, pos = ".GlobalEnv")
assign("pNO", gamlss.dist::pNO, pos = ".GlobalEnv")
assign("qNO", gamlss.dist::qNO, pos = ".GlobalEnv")
assign("rNO", gamlss.dist::rNO, pos = ".GlobalEnv")
gamlss.cens::gen.cens(NO, type = "left")
assign("NOlc",  NOlc,  pos = ".GlobalEnv")
assign("dNOlc", dNOlc, pos = ".GlobalEnv")
assign("pNOlc", pNOlc, pos = ".GlobalEnv")
assign("qNOlc", qNOlc, pos = ".GlobalEnv")

## formula 
{  # tree and forest formula
  dt.formula <- df.formula <- 
    robs ~ tppow_mean + tppow_sprd + tppow_min + tppow_max + 
    tppow_mean0612 + tppow_mean1218 + tppow_mean1824 + tppow_mean2430 + 
    tppow_sprd0612 + tppow_sprd1218 + tppow_sprd1824 + tppow_sprd2430 + 
    capepow_mean + capepow_sprd + capepow_min + capepow_max + 
    capepow_mean0612 + capepow_mean1218 + capepow_mean1224 + capepow_mean1230 +
    capepow_sprd0612 + capepow_sprd1218 + capepow_sprd1224 + capepow_sprd1230 +
    dswrf_mean_mean + dswrf_mean_max + 
    dswrf_sprd_mean + dswrf_sprd_max +
    msl_mean_mean + msl_mean_min + msl_mean_max + 
    msl_sprd_mean + msl_sprd_min + msl_sprd_max +
    pwat_mean_mean + pwat_mean_min + pwat_mean_max + 
    pwat_sprd_mean + pwat_sprd_min + pwat_sprd_max +
    tmax_mean_mean + tmax_mean_min + tmax_mean_max +
    tmax_sprd_mean + tmax_sprd_min + tmax_sprd_max +
    tcolc_mean_mean + tcolc_mean_min + tcolc_mean_max +
    tcolc_sprd_mean + tcolc_sprd_min + tcolc_sprd_max +
    t500_mean_mean + t500_mean_min + t500_mean_max +
    t700_mean_mean + t700_mean_min + t700_mean_max +
    t850_mean_mean + t850_mean_min + t850_mean_max +
    t500_sprd_mean + t500_sprd_min + t500_sprd_max +
    t700_sprd_mean + t700_sprd_min + t700_sprd_max +
    t850_sprd_mean + t850_sprd_min + t850_sprd_max +
    tdiff500850_mean + tdiff500850_min + tdiff500850_max +
    tdiff700850_mean + tdiff700850_min + tdiff700850_max +
    tdiff500700_mean + tdiff500700_min + tdiff500700_max +
    msl_diff
}


data("RainTyrol")
station <- "Axams"

# get observations and covariates for selected station
RainData <- RainTyrol[RainTyrol$station == as.character(station), ]
rownames(RainData) <- c(1:NROW(RainData))

# learning data: 24 years (1985 - 2008, both inlcuded)
# testing data: 4 successive years (2009, 2010, 2011, 2012)
learndata <- RainData[RainData$year < 2009,]
testdata <- RainData[RainData$year %in% c(2009, 2010, 2011, 2012),]

# fitting the model
set.seed(123)

### original disttree implementation with old control arguments (
#fit_time <- system.time(tr_axams <- disttree(dt.formula, 
#                                             data = learndata, family = dist_list_cens_normal, 
#                                             censtype = "left", censpoint = 0, type.tree = "ctree", 
#                                             control = ctree_control(teststat = "quad", testtype = "Bonferroni", 
#                                             intersplit = TRUE, mincriterion = 0.95, 
#                                             minsplit = 50, minbucket = 20)))

## new disttree implementation with new control arguments (
fit_time <- system.time(tr_axams <- distextree(dt.formula, 
                                               data = learndata, family = dist_list_cens_normal, 
                                               control = distextree_control(type.tree = "ctree", splittry = 2L, 
                                                                            teststat = "quad", testtype = "Bonferroni", 
                                                                            intersplit = TRUE, mincriterion = 0.95, 
                                                                            minsplit = 50, minbucket = 20)))

expect_known_value(apply(coef(tr_axams), 2, sort), file = "test_compare_disttree_with_orig_version_tr_axams_coef.rds", update = FALSE)
expect_known_value(logLik(tr_axams), file = "test_compare_disttree_with_orig_version_tr_axams_logLik.rds", update = FALSE)
expect_known_value(predict(tr_axams, newdata = testdata), file = "test_compare_disttree_with_orig_version_tr_axams_pred.rds", update = FALSE)
expect_known_value(fit_time[3L], file = "test_compare_disttree_with_orig_version_tr_axams_fittime.rds", update = FALSE, tolerance = 0.2)
expect_known_output(print(as.constparty(tr_axams)), file = "test_compare_disttree_with_orig_version_tr_axams_constparty.txt", update = FALSE, tolerance = 1e-6)

