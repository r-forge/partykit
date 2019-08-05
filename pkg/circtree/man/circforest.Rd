\name{circforest}
\alias{circforest}
\alias{predict.circforest}
\alias{logLik.circforest}

\title{Distributional Regression Forests for a Circular Response}

\description{
  Distributional forests based on maximum-likelihood estimation of parameters for
  a circular response employing the von Mises distribution.
}

\usage{
circforest(formula, data, response_range = NULL, subset, 
           na.action = na.pass, weights, offset, cluster, strata, 
           control = distextree_control(teststat = "quad", testtype = "Univ", mincriterion = 0,
           saveinfo = FALSE, minsplit = 20, minbucket = 7, splittry = 2, ...),
           ntree = 500L, fit.par = FALSE, perturb = list(replace = FALSE, fraction = 0.632),
           mtry = ceiling(sqrt(nvar)), applyfun = NULL, cores = NULL, trace = FALSE, ...)
\method{predict}{circforest}(object, newdata = NULL,
        type = c("parameter", "response", "weights", "node"),
        OOB = TRUE, scale = TRUE, \dots)
}

\arguments{
  \item{formula}{A symbolic description of the model to be fit. This
    should be of type \code{y ~ x1 + x2}
    where \code{y} should be the response variable
    and \code{x1} and \code{x2} are used as partitioning variables.}
  \item{data}{An optional data frame containing the variables in the model.}
  \item{response_range}{Circular range of response variable.}
  \item{subset}{An optional vector specifying a subset of observations to be
    used for fitting.}
  \item{na.action}{A function which indicates what should happen when the data
    contain \code{NA}s.}
  \item{weights}{Optional numeric vector of case weights.}
  \item{offset}{Optional numeric vector with \emph{a priori} known component to
    be included in the linear predictor for the location. \strong{Currently, not supported.}}
  \item{cluster}{An optional vector (typically numeric or factor) with a
    cluster ID to be employed for clustered covariances in the parameter
    stability tests.}
  \item{strata}{ an optional factor for stratified sampling.}
  \item{control}{Control arguments passed to \code{\link[partykit]{extree_fit}} via \code{\link[disttree]{distextree_control}}.}
  \item{ntree}{
    Number of trees to grow for the forest.}
  \item{fit.par}{logical. if TRUE, fitted and predicted values and predicted parameters are calculated 
    for the learning data (together with loglikelihood)}
  \item{perturb}{
    a list with arguments \code{replace} and \code{fraction} determining which type of
    resampling with \code{replace = TRUE} referring to the n-out-of-n bootstrap and
    \code{replace = FALSE} to sample splitting. \code{fraction} is the number of observations to 
    draw without replacement. }
  \item{mtry}{
    number of input variables randomly sampled as candidates
    at each node for random forest like algorithms. Bagging, as special case
    of a random forest without random input variable sampling, can
    be performed by setting \code{mtry} either equal to \code{Inf} or
    manually equal to the number of input variables.}
  \item{applyfun}{an optional \code{\link[base]{lapply}}-style function with arguments
                  \code{function(X, FUN, \dots)}. It is used for computing the variable selection criterion.
                  The default is to use the basic \code{lapply}
                  function unless the \code{cores} argument is specified (see below).}
  \item{cores}{numeric. If set to an integer the \code{applyfun} is set to
               \code{\link[parallel]{mclapply}} with the desired number of \code{cores}.}

  \item{trace}{a logical indicating if a progress bar shall be printed while
          the forest grows.}
  \item{object}{
    An object as returned by \code{distexforest}}
  \item{newdata}{
    An optional data frame containing test data.}
  \item{type}{
    a character string denoting the type of predicted value
    returned. For \code{"parameter"} the predicted distributional parameters
    are returned and for \code{"response"} the expectation is returned. \code{"weights"} returns an 
    integer vector of prediction weights. For \code{type = "node"}, a list of terminal node ids for 
    each of the trees in the forest ist returned.}
  \item{OOB}{
    a logical defining out-of-bag predictions (only if \code{newdata = NULL}).}
  \item{scale}{a logical indicating scaling of the nearest neighbor weights
               by the sum of weights in the corresponding terminal node of
               each tree. In the simple regression forest, predicting
               the conditional mean by nearest neighbor weights will be
               equivalent to (but slower!) the aggregation of means.}
  \item{\dots}{ additional arguments. }
}

\details{
  Distributional regression forests for a circular response is an application of model-based recursive 
  partitioning and unbiased recursive partitioning based on the implementation in 
  \code{\link[disttree]{distexforest}} using the infrastructure of \code{\link[partykit]{extree_fit}}.
}

\value{
  An object of S3 class \code{circforest} inheriting from class \code{distexforest}.
}

\examples{
sdat <- circtree_simulate()
m1.circforest <- circforest(y ~ x1 + x2, data = sdat, ntree = 50)
}


\seealso{
  \code{\link[disttree]{distexforest}}, \code{\link[partykit]{mob}}, \code{\link[partykit]{ctree}}
}

\keyword{forest, tree, parametric modelling, circular response}
