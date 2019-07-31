context("comparison of new implementation 'distextree()' with older version 'disttree()'")

## Check classes
tr <- distextree(dist ~ speed, data = cars)
expect_true(all(class(tr) %in% c("distextree", "modelparty", "party")))

## Compare new with old disttree for tree type equal to 'mob' (requires .mfluc_select which is not available in partykit_1.2-4)
# m.old <- disttree(dist ~ speed, data = cars, type.tree = "mob")
# m.new <- distextree(dist ~ speed, data = cars, type.tree = "mob")

# expect_equal(coef(m.old), coef(m.new))
# expect_equal(as.numeric(logLik(m.old)), as.numeric(logLik(m.new))) # attr(logLik(c.old), 'nobs') missing

## Compare new with old disttree for tree type equal to 'ctree'
## with old control arguments
c.old <- disttree(dist ~ speed, data = cars, type.tree = "ctree")
c.new <- distextree(dist ~ speed, data = cars, type.tree = "ctree", 
  control = distextree_control(minsplit = 20L, minbucket = 7L, splittry = 2L))

#expect_equal(coef(c.old), coef(c.new)) # FIXME: (ML) mismatches, average diff: 1.89e-06
expect_equal(coef(c.old), coef(c.new), tolerance = 1e-6)
expect_equal(as.numeric(logLik(c.old)), as.numeric(logLik(c.new))) # attr(logLik(c.old), 'nobs') missing

#expect_equal(attr(logLik(m.new), "df"), attr(logLik(m.old), "df"))
expect_equal(attr(logLik(c.new), "df"), attr(logLik(c.new), "df"))
#expect_equal(as.numeric(logLik(m.new)), as.numeric(logLik(m.old)))
expect_equal(as.numeric(logLik(c.new)), as.numeric(logLik(c.old)))
#expect_equal(nrow(cars), nobs(m.new))
expect_equal(nrow(cars), nobs(c.new))

if(FALSE){ 
  ## application to precipitation data
  ## FIXME: (LS) different order of terminal nodes returned by coef()
  
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
  
  ## evaluation function for one selected station
  # stationeval <- function(station) {
    
  station <- "Axams"
  
    #####
    # get observations and covariates for selected station
    RainData <- RainTyrol[RainTyrol$station == as.character(station), ]
    rownames(RainData) <- c(1:NROW(RainData))
    
    # learning data: 24 years (1985 - 2008, both inlcuded)
    # testing data: 4 successive years (2009, 2010, 2011, 2012)
    learndata <- RainData[RainData$year < 2009,]
    testdata <- RainData[RainData$year %in% c(2009, 2010, 2011, 2012),]
    
    # define matrix to store computation times
    fit_time <- matrix(ncol = 5, nrow = 2)
    colnames(fit_time) <- c("user.self", "sys.self", "elapsed", "user.child", "sys.child")
    rownames(fit_time) <- c("disttree", "distextree")
    
    #####
    # fitting the model
    set.seed(7)
    
    # fit distributional tree
    fit_time["disttree",] <- system.time(dt <- disttree(dt.formula, 
                                                        data = learndata, family = dist_list_cens_normal, 
                                                        censtype = "left", censpoint = 0, type.tree = "ctree", 
                                                        control = ctree_control(teststat = "quad", testtype = "Bonferroni", intersplit = TRUE,
                                                                                mincriterion = 0.95, minsplit = 50,
                                                                                minbucket = 20)))
    
    # fit distributional forest
    fit_time["distextree",] <- system.time(det <- distextree(dt.formula, 
                                                             data = learndata, family = dist_list_cens_normal, 
                                                             control = distextree_control(type.tree = "ctree", splittry = 2L, 
                                                                                          teststat = "quad", testtype = "Bonferroni", intersplit = TRUE,
                                                                                          mincriterion = 0.95, minsplit = 50,
                                                                                          minbucket = 20)))
    
    
    # expect_equal(coef(dt), coef(det)) ## FIXME: (LS) different order returned by coef()
    expect_equal(coef(dt)[rownames(coef(det)),], coef(det))
    expect_equal(as.numeric(logLik(dt)), as.numeric(logLik(det)))
    
    
    ## further evaluations
    
    ## get predicted parameter of all models for testdata
    
    pred_time <- matrix(ncol = 2, nrow = NROW(testdata))
    colnames(pred_time) <- c("disttree", "distextree")

    expect_equal(fit_time["disttree","elapsed"], fit_time["distextree","elapsed"], tolerance = 0.1)
    
    # disttree
    pdt <- data.frame()
    for(i in 1:NROW(testdata)){
      pred_time[i, "disttree"] <- system.time(predpar <- predict(dt, newdata = testdata[i,], type = "parameter"))["elapsed"]
      pdt <- rbind(pdt, predpar)
    }
    rownames(pdt) <- c(1:NROW(pdt))
    dt_mu <- pdt$mu
    dt_sigma <- pdt$sigma
    
    # distextree
    pdet <- data.frame()
    for(i in 1:NROW(testdata)){
      pred_time[i, "distextree"] <- system.time(predpar <- disttree:::predict.disttree(det, newdata = testdata[i,], type = "parameter"))["elapsed"]
      pdet <- rbind(pdet, predpar)
    }
    rownames(pdet) <- c(1:NROW(pdet))
    det_mu <- pdet$mu
    det_sigma <- pdet$sigma

    expect_equal(sum(pred_time[, "disttree"]), sum(pred_time[, "distextree"]), tolerance = 0.3)
    
    
    # store parameter
    par <- list(pdt = pdt,
                pdet = pdet)
    expect_equal(pdt, pdet)
    
    # CPRS
    crps_dt <- crps_cnorm(testdata$robs, location = dt_mu, scale = dt_sigma, lower = 0, upper = Inf)
    crps_det <- crps_cnorm(testdata$robs, location = det_mu, scale = det_sigma, lower = 0, upper = Inf)

    crps <- cbind(crps_dt, crps_det)
    expect_equal(crps_dt, crps_det)
    
    # loglikelihood
    dtll <- detll <- numeric(length = NROW(testdata))
    for(j in 1:(NROW(testdata))){
      
      eta_dt <- as.numeric(dist_list_cens_normal$linkfun(cbind(dt_mu, dt_sigma)[j,]))
      eta_det <- as.numeric(dist_list_cens_normal$linkfun(cbind(det_mu, det_sigma)[j,]))
      
      dtll[j] <- dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_dt, log=TRUE)
      detll[j] <- dist_list_cens_normal$ddist(testdata[j,"robs"], eta = eta_det, log=TRUE)
    }
    
    ll <- cbind(dtll, detll) 
    expect_equal(dtll, detll)
    
    colnames(ll) <- colnames(crps) <- c("disttree", "distextree")
    
  
    results <- list(station = station,
                    par = par,
                    ll = ll,
                    crps = crps,
                    fit_time = fit_time,
                    pred_time = pred_time)

}
