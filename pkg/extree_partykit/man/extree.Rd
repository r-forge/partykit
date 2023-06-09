\name{extree}
\alias{extree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Extensible Trees
}
\description{
  Extensible trees provide the basic infrastructure to define
  tree algorithms via transformation, variable selection, and
  split point selection functions.
}
\usage{
extree(data, trafo, control = extree_control(...), 
    converged = NULL, doFit = TRUE, ...)
}
\arguments{
  \item{data}{an object of class \code{\link{extree_data}}, see
\code{\link{extree_data}}.
}
  \item{trafo}{a function with arguments \code{subset}, \code{data},
	  \code{weights}, \code{info}, \code{estfun} and \code{object}.
}
  \item{control}{list of control arguments generated by
	  \code{\link{extree_control}}.  Among others the variable and split
	  point selection (\code{selectfun}, \code{splitfun}) are specified
	  here.
}
  \item{converged}{an optional function with arguments \code{subset},
	  \code{weights}.
}
  \item{...}{Additional arguments passed on to \code{extree_fit}, see \code{\link{extree_fit}}.
}
}
\details{
This basic tree algorithm can be used to define your own tree algorithm 
variants. 

\code{trafo} defines how you want to preprocess your data for variable and split
point selection.  As an example, \code{mob} computes a model and returns
information such as \code{estfun} (the empirical estimating functions / score
contribution matrix, see also \code{\link[sandwich]{estfun}}), \code{objfun}
(value of the minimized objective function, usually negative log-Likelihood),
\code{coef} (estimated model coefficients), and \code{converged} (logical, has
the model converged?).

\code{selectfun} defines how to select the split variable.
\code{splitfun} defines how to select the split point. 
Details in \code{\link{extree_control}}.
}
\value{Currently: A list of \code{nodes} (an object of class \code{partynode})
	and \code{trafo} (the encapsulated transformation function).
	This will likely change soon.
}
\examples{
data(airquality, package = "datasets")
airq <- subset(airquality, !is.na(Ozone))
airq_dat <- extree_data(Ozone ~ Wind + Temp,
  data = airq, yx = "matrix")

### Set up trafo function to preprocess data for variable and split point selection 
trafo_identity <- function(subset, data, weights = NULL, info = NULL, 
  estfun = TRUE, object = TRUE) {
  ### Extract response and "subset"
  y <- extree_variable(data, index = 1, type = "original")  
  y[-subset] <- NA  
  
  ### Return list
  return(list(
    estfun = if (estfun) y else NULL,
    unweighted = TRUE,  
    converged = TRUE 
  )
  )
}

### Set up function to guide variable selection
### To find a split variable an independence test is conducted. 
### The function returns a list with values of test statistics and p-values 
var_select_guide <- function(model, trafo, data, subset, weights, j,
  split_only = FALSE, control) {
  
  estfun <- model$estfun[subset]
  
  ### categorize estfun if not already a factor
  if(is.factor(estfun)) est_cat <- estfun else {
    breaks <- unique(quantile(estfun, c(0, 0.25, 0.5, 0.75, 1)))
    if(length(breaks) < 5) breaks <- c(min(estfun), mean(estfun), max(estfun))
    est_cat <- cut(estfun, breaks = breaks,
      include.lowest = TRUE, right = TRUE)
  }
  
  ### get possible split variable
  sv_cat <- extree_variable(data, index = j, type = "inum")[subset]
  
  ### independence test
  test <- chisq.test(x = est_cat, y = sv_cat)
  res <- list(statistic = test$statistic, p.value = test$p.value)
  
  return(res)
}


### Set up split selection
### As a split point the median is used of the split variable
split_select_median_numeric <- function(model, trafo, data, subset, weights, 
  whichvar, ctrl) {
  
  if (length(whichvar) == 0) return(NULL)
  
  ### define split point as median of variable with lowest p-value 
  j <- whichvar[1] # select variable with lowest p-value
  x <- extree_variable(data, index = j, type = "original")[subset]
  ret <- partysplit(as.integer(j), breaks = median(x))
  
  return(ret)
}

### Set extree control 
ctrl <- extree_control(criterion = "p.value", # split variable selection criterion 
  logmincriterion = log(1 - 0.05),
  update = TRUE,
  selectfun = var_select_guide,
  splitfun = split_select_median_numeric,
  svselectfun = NULL, 
  svsplitfun = NULL,
  minsplit = 50)


### Call extree 
tr <- extree(data = airq_dat, trafo = trafo_identity, 
  control = c(ctrl, restart = TRUE))
print(tr$nodes)
ptr <- party(tr$nodes, data = airq_dat$data)
print(ptr)
plot(ptr)
}
\keyword{tree}


