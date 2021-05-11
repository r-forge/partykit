# library("partykitx")
# library("tinytest")

## example as in ?lmtree
data("BostonHousing", package = "mlbench")
BostonHousing <- transform(BostonHousing,
  chas = factor(chas, levels = 0:1, labels = c("no", "yes")),
  rad = factor(rad, ordered = TRUE))
bhdat <- extree_data(medv ~ log(lstat) + I(rm^2) | zn + indus + chas + nox + 
    age + dis + rad + tax + crim + b + ptratio, data = BostonHousing,
  yx = "matrix")


## simple trafo for lmtree/glmtree
mob_trafo <- function(fit) {
  function(subset, data, weights, info, estfun, object) {
    
    ## FIXME: (AZ) Improve formula handling
    x <- data$yx$x[subset, ]
    attr(x, "formula") <- formula(data$terms$yx)
    y <- extree_variable(data, variable = "y")[subset]
    
    ## FIXME: (AZ) add start, offset, and cluster 
    trafo <- fit(y = y, 
      x = x, ## FIXME: (AZ) can we use extree_variable here?
      weights = weights[subset], estfun = estfun, object = object)
    
    if(estfun) {
      ef <- matrix(nrow = nrow(data$data), ncol = ncol(trafo$estfun))
      ef[subset, ] <- trafo$estfun
      trafo$estfun <- ef
    }
    
    return(trafo)
  }
}


## varselect m-fluctuation test
varselect_mfluc <- function(model, trafo, data, subset, weights, j, 
  split_only = FALSE, control) {
  partykitx:::.mfluc_test(model = model, trafo = trafo, data = data, subset = subset, 
    weights = weights, j = j, SPLITONLY = split_only, ctrl = control)
}


## splitselect maximising objective function
splitselect_objfun <- function(model, trafo, data, subset, weights, j, 
  split_only = TRUE, control) {
  partykitx:::.objfun_test(model = model, trafo = trafo, data = data, subset = subset, 
    weights = weights, j = j, SPLITONLY = TRUE, ctrl = control)
}




## lmtree via extree
mobtr <- extree(data = bhdat, trafo = mob_trafo(partykit:::lmfit),  
  control = extree_control(criterion = "p.value",
    critvalue = 0.05,
    update = TRUE,
    varselect = varselect_mfluc,
    splitselect = splitselect_objfun,
    svarselect = varselect_mfluc,
    ssplitselect = splitselect_objfun,
    minbucket = 40, 
    # mob specific
    trim = 0.1,
    breakties = FALSE,
    ordinal = "chisq",
    restart = TRUE,
    intersplit = FALSE))


## glmtree via extree
mobtr2 <- extree(data = bhdat, trafo = mob_trafo(partykit:::glmfit), 
  control = extree_control(criterion = "p.value",
    critvalue = 0.05,
    update = TRUE,
    varselect = varselect_mfluc,
    splitselect = splitselect_objfun,
    svarselect = varselect_mfluc,
    ssplitselect = splitselect_objfun,
    minbucket = 40, 
    # mob specific
    trim = 0.1,
    breakties = FALSE,
    ordinal = "chisq",
    restart = TRUE,
    intersplit = FALSE))


## compare to partykit::lmtree
lmtr <- partykit::lmtree(medv ~ log(lstat) + I(rm^2) | zn + indus + chas + nox +
    age + dis + rad + tax + crim + b + ptratio, data = BostonHousing, minsize = 40)

# p-values look good
test1_mobtr <- nodeapply(mobtr$nodes, ids = 1, FUN = info_node)[[1]]$criterion
(sorted_pvalues_mobtr <- sort(test1_mobtr["p.value", ]))

test1_lmtr <- nodeapply(lmtr, ids = 1, FUN = info_node)[[1]]$test
(sorted_pvalues_lmtr <- sort(test1_lmtr["p.value", ]))

expect_equal(sorted_pvalues_mobtr, sorted_pvalues_lmtr, tolerance = 0.001)
expect_equal(names(sorted_pvalues_mobtr), names(sorted_pvalues_lmtr))

mobtr
lmtr

print(party(mobtr$nodes, data = bhdat$data), 
  terminal_panel = function(node) c(sprintf(": n = %s", node$info$nobs), 
    capture.output(print(node$info$coefficients))))

print(party(mobtr2$nodes, data = bhdat$data), 
  terminal_panel = function(node) c(sprintf(": n = %s", node$info$nobs), 
    capture.output(print(node$info$coefficients))))

expect_equal(width(mobtr$node), width(lmtr))
## Note: width can differ if split cannot be found in top selected variable
## because extree will also look at the second best (and so on), but old mob
## does not.



#### simple example with using strucchange::gefp #####

bhdat2 <- extree_data(medv ~ log(lstat) + I(rm^2) | tax + ptratio, data = BostonHousing,
  yx = "matrix")

## extree glmfit (for comparison)
tr1 <- extree(data = bhdat2, trafo = mob_trafo(partykit:::glmfit), 
  control = extree_control(criterion = "p.value",
    critvalue = 0.05,
    update = TRUE,
    varselect = varselect_mfluc,
    splitselect = splitselect_objfun,
    svarselect = varselect_mfluc,
    ssplitselect = splitselect_objfun,
    minbucket = 40, 
    # mob specific
    trim = 0.1,
    breakties = FALSE,
    ordinal = "chisq",
    restart = TRUE,
    intersplit = FALSE))


## varselect with strucchange::gefp
varselect_strucchange_numeric <- function(model, trafo, data, subset, weights, j, 
  split_only = FALSE, control) {
  
  z <- extree_variable(data, index = j)[subset]
  
  fluc <- strucchange::gefp(model$object, order.by = z, fit = NULL, 
    scores = function(...) model$estfun[subset, ])
  
  strucchange::sctest(fluc, 
    functional = strucchange::supLM(control$trim))[c("statistic", "p.value")]
  
}


## tree that uses strucchange::gefp
tr2 <- extree(data = bhdat2, trafo = mob_trafo(partykit:::glmfit), 
  control = extree_control(criterion = "p.value",
    critvalue = 0.05,
    update = TRUE,
    varselect = varselect_strucchange_numeric,
    splitselect = splitselect_objfun,
    minbucket = 40, 
    # mob specific
    trim = 0.1,
    breakties = FALSE,
    ordinal = "chisq",
    restart = TRUE,
    intersplit = FALSE))



## compare the two trees
print(party(tr1$nodes, data = bhdat2$data), 
      terminal_panel = function(node) c(sprintf(": n = %s", node$info$nobs), 
                                        capture.output(print(node$info$coefficients))))

print(party(tr2$nodes, data = bhdat2$data), 
  terminal_panel = function(node) c(sprintf(": n = %s", node$info$nobs), 
    capture.output(print(node$info$coefficients))))


## test that width is the same
expect_equal(width(tr1$node), width(tr2$node))


## check p-values in first split
# p-values look good
test1_tr1 <- nodeapply(tr1$nodes, ids = 1, FUN = info_node)[[1]]$criterion
(sorted_pvalues_tr1 <- sort(test1_tr1["p.value", ]))

test1_tr2 <- nodeapply(tr2$nodes, ids = 1, FUN = info_node)[[1]]$criterion
(sorted_pvalues_tr2 <- sort(test1_tr2["p.value", ]))

expect_equal(sorted_pvalues_tr1, sorted_pvalues_tr2, tolerance = 0.001)
