context("comparison of new implementation 'distexforest()' with older version 'distforest()'")

library("disttree")

## load and prepare data (rain data for observation station Axams, see demo distforest)

# load observations and covariates 
data("RainAxams", package = "disttree")

# define function for parallelization
applyfun <- function(X, FUN, ...) parallel::mclapply(X, FUN, ..., mc.cores = pmax(1, detectCores() - 1))

# formula
{
  df.formula <- 
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



#####
# further packages
library("gamlss")
library("gamlss.dist")

# if gamlss.cens family object should be used as family
library("gamlss.cens")
gen.cens(NO, type = "left")


# learning data: 24 years (1985 - 2008, both inlcuded)
# testing data: 4 successive years (2009, 2010, 2011, 2012)
learndata <- RainAxams[RainAxams$year < 2009,]
testdata <- RainAxams[RainAxams$year %in% c(2009, 2010, 2011, 2012),]

## fit four distributional forests:
# distforest vs. distexforest
# prepared family list vs. gamlss.family object

set.seed(7)
df <- distforest(df.formula, 
                 data = learndata, family = dist_list_cens_normal, type.tree = "ctree", 
                 ntree = 100, censtype = "left", censpoint = 0, mtry = 27,
                 control = ctree_control(teststat = "quad", testtype = "Univ", intersplit = TRUE,
                                         mincriterion = 0, minsplit = 50,
                                         minbucket = 20))
set.seed(7)
dfg <- distforest(df.formula, 
                  data = learndata, family = NOlc(), type.tree = "ctree", 
                  ntree = 100, censtype = "left", censpoint = 0, mtry = 27,
                  control = ctree_control(teststat = "quad", testtype = "Univ", intersplit = TRUE,
                                          mincriterion = 0, minsplit = 50,
                                          minbucket = 20))

set.seed(7)
def <- distexforest(df.formula, 
                    data = learndata, family = dist_list_cens_normal, 
                    ntree = 100, mtry = 27,
                    control = distextree_control(teststat = "quad", testtype = "Univ", intersplit = TRUE,
                                                 mincriterion = 0, minsplit = 50,
                                                 minbucket = 20))

set.seed(7)
defg <- distexforest(df.formula, 
                     data = learndata, family = distfamily(family = NOlc(), censpoint = 0), 
                     ntree = 100, mtry = 27,
                     control = distextree_control(teststat = "quad", testtype = "Univ", intersplit = TRUE,
                                                  mincriterion = 0, minsplit = 50,
                                                  minbucket = 20))

pdf <- df$fitted.par
pdef <- def$fitted.par
pdfg <- dfg$fitted.par
pdefg <- defg$fitted.par
head(cbind(pdf,pdef,pdfg,pdefg))

# differences between distforest and distexforest
expect_equal(pdf,pdef)
expect_equal(pdfg,pdefg)

# differences between family objects
expect_equal(pdf,pdfg)
expect_equal(pdef,pdefg)

sum(abs(pdf-pdef))
sum(abs(pdfg-pdefg))
sum(abs(pdf-pdfg))
sum(abs(pdef-pdefg))

# using newdata
pdf2 <- predict(df, type = "parameter", newdata = testdata)
pdef2 <- predict(def, type = "parameter", newdata = testdata)
head(cbind(pdf2,pdef2))
sum(abs(pdf2-pdef2))
expect_equal(pdf2,pdef2)



## Check classes
expect_true(all(class(df) %in% c("distforest", "cforest", "constparties", "parties")))
expect_true(all(class(dfg) %in% c("distforest", "cforest", "constparties", "parties")))
expect_true(all(class(def) %in% c("distexforest", "constparties", "parties")))
expect_true(all(class(defg) %in% c("distexforest", "constparties", "parties")))

