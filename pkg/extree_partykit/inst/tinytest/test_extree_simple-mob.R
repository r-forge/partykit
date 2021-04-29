library("partykitx")

data("BostonHousing", package = "mlbench")
BostonHousing <- transform(BostonHousing,
  chas = factor(chas, levels = 0:1, labels = c("no", "yes")),
  rad = factor(rad, ordered = TRUE))
bhdat <- extree_data(medv ~ lstat + rm | zn + indus + chas + nox + 
    age + dis + rad + tax + crim + b + ptratio, data = BostonHousing)
## FIXME: (HS) original example in (?lmtree) does not work here because data
## are already transformed (log(lstat) instead lstat in bhdat$data)


simple_lm_trafo <- function(subset, data, weights, info, estfun, object) {
  
  if(is.null(weights)) {
    obj <- lm(data$terms$yx, data = data$data[subset, ]) ## ? why does subset = subset not work?
  } else {
    obj <- lm(data$terms$yx, data = data$data[subset, ], weights = weights)
  }
  
  if(estfun) {
    ef_subset <- sandwich::estfun(obj) 
    ef <- matrix(nrow = nrow(data$data), ncol = ncol(ef_subset))
    ef[subset, ] <- ef_subset
  } 
  
  list(
    coefficients = coef(obj),
    objfun = -as.numeric(logLik(obj)),
    estfun = if(estfun) ef else NULL,
    object = if(object) obj else NULL
  )
}

varselect_mfluc <- function(model, trafo, data, subset, weights, j, 
  split_only = FALSE, control) {
  partykitx:::.mfluc_test(model = model, trafo = trafo, data = data, subset = subset, 
    weights = weights, j = j, SPLITONLY = split_only, ctrl = control)
}

splitselect_objfun <- function(model, trafo, data, subset, weights, j, 
  split_only = TRUE, control) {
  partykitx:::.objfun_test(model = model, trafo = trafo, data = data, subset = subset, 
    weights = weights, j = j, SPLITONLY = TRUE, ctrl = control)
}





mobtr <- extree(data = bhdat, trafo = simple_lm_trafo, 
  control = extree_control(criterion = "p.value",
    logmincriterion = log(1 - 0.05),
    update = TRUE,
    varselect = varselect_mfluc,
    splitselect = splitselect_objfun,
    svarselect = varselect_mfluc,
    ssplitselect = splitselect_objfun,
    # mob specific
    trim = 0.1,
    breakties = FALSE,
    ordinal = "chisq",
    restart = TRUE,
    intersplit = FALSE))

lmtr <- partykit::lmtree(medv ~ lstat + rm | zn + indus + chas + nox + 
    age + dis + rad + tax + crim + b + ptratio, data = BostonHousing)

# We still have some issues with the p-values, but the ordering seems fine! :)
test1_mobtr <- nodeapply(mobtr$nodes, ids = 1, FUN = info_node)[[1]]$criterion
(sorted_pvalues_mobtr <- sort(test1_mobtr["p.value", ]))

test1_lmtr <- nodeapply(lmtr, ids = 1, FUN = info_node)[[1]]$test
(sorted_pvalues_lmtr <- sort(test1_lmtr["p.value", ]))

expect_equal(sorted_pvalues_mobtr, sorted_pvalues_lmtr)
expect_equal(names(sorted_pvalues_mobtr), names(sorted_pvalues_lmtr))
