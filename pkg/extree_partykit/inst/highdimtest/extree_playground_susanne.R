# -------------------------------------------------------------------
# - NAME:   extree_playground_susanne.R
# - AUTHOR: Susanne
# - DATE:   2020-09-23
# -------------------------------------------------------------------
# - PURPOSE: Some examples illustrating high-dimension performance
# - BASIS: ../19_extree_hackathon/extree_examples_2019_11_08.R
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Prelimanaries
# -------------------------------------------------------------------
library(devtools)
devtools::load_all()
# library(partykit)

#--------------------------------------------------------------------
# Data 
#--------------------------------------------------------------------
data("airquality", package = "datasets")
nsim <- nrow(airquality)

add_noise_var <- function(data, dim) {
  nsim <- nrow(data)
  x <- matrix(runif(nsim * dim), nrow = nsim, ncol = dim)
  colnames(x) <- paste("X", 1:ncol(x), sep = "")
  cbind(data, x)
}

# -------------------------------------------------------------------
# EXAMPLE 1: Use trafo_identity() and guide_test() for numerical response
# -------------------------------------------------------------------
## Set up variable selection list
var_select_guide_numeric <- function(model, trafo, data, subset, weights, j, 
                                     split_only = FALSE, control) {
  
  estfun <- model$estfun[subset]
  
  # categorize estfun if not already a factor
  est_cat <- if (is.factor(estfun)) estfun else cut(estfun, 
                                                    breaks = c(min(estfun), mean(estfun), max(estfun)), 
                                                    include.lowest = TRUE, right = TRUE)
  
  ## get possible split variable
  # sv <- data$zindex[[j]][subset] ## FIXME: always the same as below?
  sv <- extree_variable(data, i = j, type = "index")[subset]
  ## FIXME: can copying be avoided here^?
  
  ## categorize possible split variable
  breaks <- unique(quantile(sv, c(0, 0.25, 0.5, 0.75, 1)))
  if (length(breaks) < 5) breaks <- c(min(sv), mean(sv), max(sv))
  sv_cat <- cut(sv, breaks = breaks, 
                include.lowest = TRUE, right = TRUE)
  
  ## independence test 
  test <- chisq.test(x = est_cat, y = sv_cat) ## FIXME: use libcoin?
  res <- list(statistic = test$statistic, p.value = test$p.value) 
  ## FIXME: (ML, LS) return log(1 - p-value) instead?
  
  return(res)
}

var_select_guide_factor <- function(model, trafo, data, subset, weights, j, 
                                    split_only = FALSE, control) {
  
  estfun <- model$estfun[subset]
  
  ## categorize estfun if not already a factor
  if(is.factor(estfun)) est_cat <- estfun else {
    breaks <- unique(quantile(estfun, c(0, 0.5, 1)))
    if(length(breaks) < 3) breaks <- c(min(estfun), mean(estfun), max(estfun))
    est_cat <- cut(estfun, breaks = breaks, 
                   include.lowest = TRUE, right = TRUE)
  }
  
  ## get possible split variable
  # sv_cat <- data$zindex[[j]][subset] ## FIXME: always the same as below?
  sv_cat <- extree_variable(data, i = j, type = "index")[subset]
  ## FIXME: can copying be avoided here^?
  
  ## independence test
  test <- chisq.test(x = est_cat, y = sv_cat)
  res <- list(statistic = test$statistic, p.value = test$p.value) 
  ## FIXME: (ML, LS) return log(1 - p-value) instead?
  
  return(res)
}


## Set up split selection: split_select with median without index j
split_select_median_numeric <- function(model, trafo, data, subset, weights, 
                                        whichvar, ctrl) {
  
  if (length(whichvar) == 0) return(NULL)
  
  ## split FIRST variable at median 
  ## FIXME: loop if necessary
  j <- whichvar[1]
  x <- model.frame(data)[[j]][subset] ## FIXME: better via extree_variable or model.frame?
  ret <- partysplit(as.integer(j), breaks = median(x))
  
  return(ret)
}


## trafo
trafo_identity <- function(subset, data, weights = NULL, info = NULL, 
                           estfun = TRUE, object = TRUE) {
  
  ## Extract response and "subset"
  y <- extree_variable(data, i = 1, type = "original")  # FIXME: (ML, LS) data copy? no aggregation possible!
  y[-subset] <- NA  
  
  ## Return list
  rval <- list(
    estfun = if (estfun) y else NULL,
    unweighted = TRUE,  
    converged = TRUE 
  )
  
  return(rval)
}

#-----------------------------------------------------------------------------
# Try out high-dim datasets
#-----------------------------------------------------------------------------

ctrl1 <- extree_control(criterion = "p.value",
  critvalue = 0.05,
  update = TRUE,
  selectfun = list(
    numeric = var_select_guide_numeric,
    default = var_select_guide_factor
  ),
  svselectfun = NULL, 
  svsplitfun = NULL,
  splitfun = split_select_median_numeric,
  minsplit = 50)

dim <- 100
yx <- "matrix"
bins <- Inf
datahd <- add_noise_var(airquality, dim) 

## One integer target: 
# Binning results in errors & wrong bins in case of integer outcomes
dl <- list(y = "Ozone", z = c("Wind", "Temp", paste("X", 1:dim, sep = "")))
extree_data(dl, data = datahd, yx = "matrix", 
  nmax = c("yx" = 10, "z" = Inf), ytype = "matrix") 

airq_dat2 <- extree_data(dl, data = datahd, yx = "none", # no error but also no binning
  nmax = c("yx" = 10, "z" = Inf), ytype = "matrix") 
airq_dat3 <- extree_data(dl, data = datahd, yx = "matrix", # no error but wrong output
  nmax = c("yx" = 10, "z" = Inf), ytype = "vector") 
airq_dat4 <- extree_data(dl, data = datahd, yx = "matrix", # works because no binning
  ytype = "matrix") 

# One numeric target: 
datahd$Ozonen = as.numeric(datahd$Ozone)
dln <- list(y = "Ozonen", z = c("Wind", "Temp"))
extree_data(dln, data = datahd, yx = "matrix", 
  nmax = c("yx" = 10), ytype = "matrix") 
# Error in if (length(ux) > nmax) ux <- unique(quantile(x, prob = 1:(nmax -  : 
#     missing value where TRUE/FALSE needed
extree_data(dln, data = datahd, yx = "matrix", #numeric outcome 
  nmax = c("yx" = 10), ytype = "vector") # error due to binning
# Error in if (length(ux) > nmax) ux <- unique(quantile(x, prob = 1:(nmax -  : 
#     missing value where TRUE/FALSE needed

## Multiple targets: 
datahd$Ozone2 <- sample(datahd$Ozone)
dlmulti <- list(y = c("Ozonen", "Ozone2"), z = c("Wind", "Temp", paste("X", 1:dim, sep = "")))
airq_dat <- extree_data(dlmulti, data = datahd, yx = "matrix",
  nmax = c("yx" = c(10, 10, 20), "z" = Inf), ytype = "matrix")

############################################################################

# Test predict: 
dl <- list(y = "Ozone", z = c("Wind", "Temp"))
airq_dat <- extree_data(Ozone ~ Wind + Temp, data = airquality)

xtr <- extree(data = airq_dat, trafo = trafo_identity, 
  control = ctrl1)
# class list 

ptr <- party(xtr$nodes, data = airq_dat$data)

pred <- predict(ptr)
predict(ptr, type = "response")
predict(ptr, type = "node")

############################################################################

#--- Runtime Study ---- 
library(ggplot2)
library(data.table) 
library(tictoc)
library(gridExtra)

fit_tree_airquality <- function(dim, bins = Inf, yx = "none") {
  datahd <- add_noise_var(airquality, dim) 
  dl <- list(y = "Ozone", z = c("Wind", "Temp", paste("X", 1:dim, sep = "")))
  
  
  tic()
  airq_dat <- extree_data(dl, data = datahd, yx = yx,
    nmax = c(yx = bins, z = bins))
  
  a <- toc()
  tic()
  xtr <- extree(data = airq_dat, trafo = trafo_identity,
    control = ctrl1)
  b <- toc()
  ptr <- party(xtr$nodes, data = airq_dat$data)
  tic()
  pred <- predict(ptr)
  c <- toc()
  list(
    time_data = as.numeric(a$toc - a$tic), 
    time_extree = as.numeric(b$toc - b$tic), 
    time_pred = as.numeric(c$toc - c$tic))
}

dim <- c(1000, 5000, 10000, 20000)
# bins <- c(10, 100, Inf)
yx <- c("matrix", "none")

# grid <- expand.grid(dim, bins, yx)
# grid.time <- as.data.table(grid)
# names(grid.time) = c("dims", "bins", "yx")

grid <- expand.grid(dim, yx)
grid.time <- as.data.table(grid)
names(grid.time) = c("dims", "yx")
grid.time[,  c("time_data", "time_extree", "time_predict") := data.table::rbindlist(lapply(seq_len(nrow(grid.time)), function(row) {
  setup <- grid.time[row, ]
  xtree <- fit_tree_airquality(dim = setup$dims, yx = as.character(setup$yx)) #, bins = setup$bins
  list(time_data = xtree$time_data,
    time_extree = xtree$time_extree, 
    time_predict = xtree$time_pred)
}))
]


# grid.time$bins = as.factor(grid.time$bins)

a <- ggplot(grid.time, aes(x = dims, y = time_data, color = yx)) + 
  # geom_point(aes(shape = bins)) +
  geom_point() + 
  ylab("extree_data time elapsed (sec)")

b <- ggplot(grid.time, aes(x = dims, y = time_extree, color = yx)) +
  # geom_point(aes(shape = bins)) +
  geom_point() + 
  ylab("extree time elapsed (sec)")

c <- ggplot(grid.time, aes(x = dims, y = time_predict, color = yx)) +
  geom_point() + 
  ylab("predict time elapsed (sec)")

d <- grid.arrange(a, b, c, ncol = 3)

# ggsave(filename = "extree_data_time.pdf", plot = a, height = 3.5, width = 4)
# ggsave(filename = "extree_time.pdf", plot = b, height = 3.5, width = 4)
# ggsave(filename = "extree_pedict_time.pdf", plot = c, height = 3.5, width = 4)
ggsave(filename = "time.pdf", plot = d, height = 3.5, width = 12)
