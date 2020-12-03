# -------------------------------------------------------------------
# - NAME:   test_compare_distforest_with_orig_version.R 
# - AUTHOR: Moritz N. Lang, Lisa Schlosser
# - DATE:   2020-11-12
# -------------------------------------------------------------------
# - PURPOSE: Test distforest(): Compare new implementation (build on extree())
#            with original implementation build on ctree() using 'tinytest'.
# -------------------------------------------------------------------
# - REMARK: Tests are only run at home, where environmentent variable TT_AT_HOME=TRUE.
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# PRELIMINARIES
# -------------------------------------------------------------------
## To get the same seeds as in the saved versions
suppressWarnings(RNGversion("3.5.0"))

## Load necessary packages  ## Must be in Depends, Imports or Suggests
library("crch")
library("scoringRules")
library("gamlss.cens")
library("RainTyrol")

## Helper function to produce print statement
get_print_statement <- function(print_statement){
  sink(tmp_file <- tempfile())
  on.exit(sink()) 
  invisible(force(print(print_statement))) 

  con <- file(tmp_file, "r", encoding = "Latin1")
  output <- readLines(con); close(con)

  return(output)
} 

## Helper function to suppress cat() [original written by Hadley Wickham]
suppress_cat <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 


# -------------------------------------------------------------------
# CARS EXAMPLE
# -------------------------------------------------------------------
set.seed(123)

# -------------------------------------------------------------------
# Fit forests
# -------------------------------------------------------------------
## Original distforest implementation with old control arguments 
#f_cars_old <- distforest(dist ~ speed, data = cars, type.tree = "ctree")

## New disttree implementation with new control arguments
f_cars_new <- distforest(dist ~ speed, data = cars, type.tree = "ctree")

# -------------------------------------------------------------------
# Run tests
# -------------------------------------------------------------------
## Compare logLik
con <- file("files/test_compare_distforest_with_orig_version_f_cars_logLik.txt", "r", 
  encoding = "Latin1")
logLik_print_cars_old <- readLines(con); close(con)

expect_identical(get_print_statement(logLik(f_cars_new)), logLik_print_cars_old)

## Compare predicted parameter
expect_equal(predict(f_cars_new, type = "parameter", newdata = cars), 
             readRDS(file = "files/test_compare_distforest_with_orig_version_f_cars_par.rds"))

## Compare predicted response
expect_equal(predict(f_cars_new, type = "response", newdata = cars),
             readRDS(file = "files/test_compare_distforest_with_orig_version_f_cars_response.rds"))

## Compare predicted node
expect_equal(predict(f_cars_new, type = "node", newdata = cars), 
             readRDS(file = "files/test_compare_distforest_with_orig_version_f_cars_nodes.rds"))

## Compare predicted weights
expect_equal(predict(f_cars_new, type = "weights", newdata = cars), 
             readRDS(file = "files/test_compare_distforest_with_orig_version_f_cars_weights.rds"))


# -------------------------------------------------------------------
# APPLICATION FOR PRECIPITATION FORECASTING AT AXAMS
# -------------------------------------------------------------------
set.seed(123)

## Assign families
assign("NO",  gamlss.dist::NO,  pos = ".GlobalEnv")
assign("dNO", gamlss.dist::dNO, pos = ".GlobalEnv")
assign("pNO", gamlss.dist::pNO, pos = ".GlobalEnv")
assign("qNO", gamlss.dist::qNO, pos = ".GlobalEnv")
assign("rNO", gamlss.dist::rNO, pos = ".GlobalEnv")

suppress_cat(gamlss.cens::gen.cens(NO, type = "left"))
assign("NOlc",  NOlc,  pos = ".GlobalEnv")
assign("dNOlc", dNOlc, pos = ".GlobalEnv")
assign("pNOlc", pNOlc, pos = ".GlobalEnv")
assign("qNOlc", qNOlc, pos = ".GlobalEnv")

## Create formula 
df_formula <- 
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

## Load data
data("RainTyrol")
station <- "Axams"

## Get observations and covariates for selected station
RainData <- RainTyrol[RainTyrol$station == as.character(station), ]
rownames(RainData) <- c(1:NROW(RainData))

## Subset to learning and testing data
learndata <- RainData[RainData$year < 2009,]  # 24 years (1985 - 2008)
testdata <- RainData[RainData$year %in% c(2009, 2010, 2011, 2012),]  # 4 successive years

# -------------------------------------------------------------------
# Fit forests
# -------------------------------------------------------------------
## Original distforest implementation (w/ old control arguments)
#fit_time_old <- system.time(
#  f_axams_old <- distforest(df_formula, 
#                            data = learndata, family = dist_list_cens_normal, 
#                            censtype = "left", censpoint = 0, type.tree = "ctree", 
#                            control = ctree_control(teststat = "quad", testtype = "Bonferroni", 
#                                                    intersplit = TRUE, mincriterion = 0.95, 
#                                                    minsplit = 70, minbucket = 40), 
#                            ntree = 50)
#)

## New distforest implementation (w/ new control arguments)
fit_time_new <- system.time(
  f_axams_new <- distforest(df_formula, 
                            data = learndata, family = dist_list_cens_normal, 
                            control = disttree_control(type.tree = "ctree", splittry = 2L, 
                                                       teststat = "quad", testtype = "Bonferroni", 
                                                       intersplit = TRUE, mincriterion = 0.95, 
                                                       minsplit = 70, minbucket = 40), 
                            ntree = 50)
)

# -------------------------------------------------------------------
# Run tests
# -------------------------------------------------------------------
## Compare logLik
con <- file("files/test_compare_distforest_with_orig_version_f_axams_logLik.txt", "r", encoding = "Latin1")
logLik_print_axams_old <- readLines(con); close(con)

expect_identical(get_print_statement(logLik(f_axams_new, newdata = learndata)), 
  logLik_print_axams_old)

## Compare predicted parameter
expect_equal(predict(f_axams_new, type = "parameter", newdata = testdata), 
             readRDS(file = "files/test_compare_distforest_with_orig_version_f_axams_par.rds"))

## Compare predicted response
expect_equal(predict(f_axams_new, type = "response", newdata = testdata), 
             readRDS(file = "files/test_compare_distforest_with_orig_version_f_axams_response.rds"))

## Compare predicted node
expect_equal(predict(f_axams_new, type = "node", newdata = learndata), 
             readRDS(file = "files/test_compare_distforest_with_orig_version_f_axams_nodes.rds"))

## Compare predicted weights
expect_equal(predict(f_axams_new, type = "weights", newdata = learndata), 
             readRDS(file = "files/test_compare_distforest_with_orig_version_f_axams_weights.rds"))

## Check for estimation time
expect_true(fit_time_new[3L] < 40)
