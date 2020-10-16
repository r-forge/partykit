# -------------------------------------------------------------------
# - NAME:   extree_playground_susanne.R
# - AUTHOR: Susanne
# - DATE:   2020-09-23
# -------------------------------------------------------------------
# - PURPOSE: Some examples illustrating high-dimension performance
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

fit_tree_airquality <- function(dim, bins = Inf, yx = "None") {
  
  datahd <- add_noise_var(airquality, dim) 
  dl <- list(y = "Ozone", z = c("Wind", "Temp", paste("X", 1:dim, sep = "")))
  
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
  
  tic()
  airq_dat <- extree_data(dl, data = datahd, yx = yx,
                          nmax = c(yx = bins, z = bins))
  
  a <- toc()
  time_data <- as.numeric(a$toc - a$tic)
  tic()
  xtr <- extree(data = airq_dat, trafo = trafo_identity, 
                control = c(ctrl1, restart = TRUE))
  b <- toc()
  list(time_data = time_data, time_extree = as.numeric(b$toc - b$tic))
}

# Study: 

dim <- c(1000, 5000, 10000)
bins <- c(10, 100, Inf)
yx <- c("none", "matrix")
grid <- expand.grid(dim, bins, yx)
grid.time <- as.data.table(grid)
names(grid.time) = c("dims", "bins", "yx")
grid.time[,  c("time_data", "time_extree") := data.table::rbindlist(lapply(seq_len(nrow(grid.time)), function(row) {
  setup <- grid.time[row, ]
  xtree <- fit_tree_airquality(dim = setup$dims, bins = setup$bins, yx = as.character(setup$yx))  
  xtree$time_data
  list(
    time_data = xtree$time_data, 
    time_extree = xtree$time_extree
  )
}))
]

grid.time$bins = as.factor(grid.time$bins)

a <- ggplot(grid.time, aes(x = dims, y = time_data, color = yx)) + 
  geom_point(aes(shape = bins)) +
  ylab("extree_data time elapsed (sec)")

b <- ggplot(grid.time, aes(x = dims, y = time_extree, color = yx)) +
  geom_point(aes(shape = bins)) +
  ylab("extree time elapsed (sec)")

ggsave(filename = "extree_data_time.pdf", plot = a, height = 3.5, width = 4)
ggsave(filename = "extree_time.pdf", plot = b, height = 3.5, width = 4)


# Profiling: 
dim <- 100
yx <- "matrix"
bins <- Inf
datahd <- add_noise_var(airquality, dim) 
datahd$Ozone2 <- sample(datahd$Ozone)
dl <- list(y = c("Ozone", "Ozone2"), z = c("Wind", "Temp", paste("X", 1:dim, sep = "")))

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

system.time(airq_dat <- extree_data(dl, data = datahd, yx = yx,
  nmax = c(yx = bins, z = bins), ytype = "matrix"))

# ERROR occurs: 
dl2 <- list(y = c("Ozone", "Ozone2"), x = "Day", 
  z = c("Wind", "Temp", paste("X", 1:dim, sep = "")))

airq_dat2 <- extree_data(dl2, data = datahd, yx = yx,
  nmax = c(yx = bins, z = bins), ytype = "matrix") 


# Improvement? 
# user  system elapsed 
# 22.128   0.006  22.138 

# user  system elapsed 
# 1.586   0.000   1.586 


# Test: 
head(datahd$Ozone)
head(datahd$Ozone2)
head(airq_dat$yx$y)
