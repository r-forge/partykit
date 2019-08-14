\name{distforest}
\alias{distforest}
\alias{predict.distforest}
\alias{logLik.distforest}
\alias{gettree.distforest}
\alias{varimp.distforest}
\encoding{latin1}

\title{Distributional Regression Forests}

\description{
  Forests based on maximum-likelihood estimation of parameters for
  specified distribution families, for example from the GAMLSS family 
  (for generalized additive models for location, scale, and shape).
}

\usage{
distforest(formula, data, subset, na.action = na.pass, weights,
             offset, cluster, family = NO(), strata, 
             control = disttree_control(teststat = "quad", testtype = "Univ", 
             mincriterion = 0, saveinfo = FALSE, minsplit = 20, minbucket = 7, 
             splittry = 2, ...), 
             ntree = 500L, fit.par = FALSE, 
             perturb = list(replace = FALSE, fraction = 0.632), 
             mtry = ceiling(sqrt(nvar)), applyfun = NULL, cores = NULL, 
             trace = FALSE, ...)
\method{predict}{distforest}(object, newdata = NULL,
        type = c("parameter", "response", "weights", "node"),
        OOB = FALSE, scale = TRUE, \dots)
}

\arguments{
  \item{formula}{a symbolic description of the model to be fit. This
    should be of type \code{y ~ x1 + x2}
    where \code{y} should be the response variable
    and \code{x1} and \code{x2} are used as partitioning variables.}
  \item{data}{a data frame containing the variables in the model.}
  \item{subset}{an optional vector specifying a subset of observations to be
                 used in the fitting process.}
  \item{na.action}{a function which indicates what should happen when the data
                   contain missing value.}
  \item{weights}{
    an optional vector of weights to be used in the fitting
    process. Non-negative integer valued weights are
    allowed as well as non-negative real weights.
    Observations are sampled (with or without replacement)
    according to probabilities \code{weights / sum(weights)}.
    The fraction of observations to be sampled (without replacement)
    is computed based on the sum of the weights if all weights
    are integer-valued and based on the number of weights greater zero
    else. Alternatively, \code{weights} can be a double matrix defining
    case weights for all \code{ncol(weights)} trees in the forest directly.
    This requires more storage but gives the user more control.}
  \item{offset}{an optional vector of offset values.}
  \item{cluster}{an optional factor indicating independent clusters.
                  Highly experimental, use at your own risk.}
  \item{family}{specification of the response distribution.
    Either a \code{\link[gamlss.dist]{gamlss.family}} object, a list generating function or a family list.}
  \item{strata}{an optional factor for stratified sampling.}
  \item{control}{a list with control parameters, see
                 \code{\link{disttree_control}}. The default values that are not set within
                 the call of \code{\link[disttree]{distforest}} correspond to those
                 of the default values used by \code{\link[disttree]{disttree}} from the
                 \code{disttree} package. \code{saveinfo = FALSE} leads to less
                 memory hungry representations of trees. Note that arguments
                 \code{mtry}, \code{cores} and \code{applyfun} in
                 \code{\link{disttree_control}} are ignored for \code{\link{distforest}},
                 because they are already set.}
  \item{ntree}{
    number of trees to grow for the forest.}
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
    an object as returned by \code{distforest}}
  \item{newdata}{
    an optional data frame containing test data.}
  \item{type}{
    a character string denoting the type of predicted value
    returned. For \code{"parameter"} the predicted distributional parameters
    are returned and for \code{"response"} the expectation is returned. \code{"weights"} returns an 
    integer vector of prediction weights. For \code{type = "node"}, a list of terminal node ids for 
    each of the trees in the forest is returned.}
  \item{OOB}{
    a logical defining out-of-bag predictions (only if \code{newdata = NULL}).}
  \item{scale}{a logical indicating scaling of the nearest neighbor weights
               by the sum of weights in the corresponding terminal node of
               each tree. In the simple regression forest, predicting
               the conditional mean by nearest neighbor weights will be
               equivalent to (but slower!) the aggregation of means.}
  \item{\dots}{arguments to be used to form the default \code{control} argument
    if it is not supplied directly.}
}

\details{
  Distributional regression forests are an application of model-based recursive partitioning
  (implemented in \code{\link[partykit]{mob}}, \code{\link[partykit]{ctree}} and \code{\link[partykit]{cforest}}) to parametric model fits based on the GAMLSS family of distributions.
  
  Distributional regression trees, see \code{\link{disttree}}, are fitted to each
  of the \code{ntree} perturbed samples of the learning sample. Most of the hyper parameters in
  \code{\link{disttree_control}} regulate the construction of the distributional regression trees.

  Hyper parameters you might want to change are:

  1. The number of randomly preselected variables \code{mtry}, which is fixed
     to the square root of the number of input variables.

  2. The number of trees \code{ntree}. Use more trees if you have more variables.

  3. The depth of the trees, regulated by \code{mincriterion}. Usually unstopped and unpruned
     trees are used in random forests. To grow large trees, set \code{mincriterion} to a small value.


  The aggregation scheme works by averaging observation weights extracted
  from each of the \code{ntree} trees and NOT by averaging predictions directly
  as in \code{\link[randomForest]{randomForest}}.
  See Schlosser et al. (2019), Hothorn et al. (2004), and Meinshausen (2006) for a description.

  Predictions can be computed using \code{\link{predict}}. For observations
  with zero weights, predictions are computed from the fitted tree
  when \code{newdata = NULL}.
}

\value{
  An object of class \code{distforest}.
}

\references{
Breiman L (2001).
  Random Forests.
  \emph{Machine Learning}, \bold{45}(1), 5--32.

Hothorn T, Lausen B, Benner A, Radespiel-Troeger M (2004).
  Bagging Survival Trees.
  \emph{Statistics in Medicine}, \bold{23}(1), 77--91.

Hothorn T, B{\"u}hlmann P, Dudoit S, Molinaro A, Van der Laan MJ (2006a).
  Survival Ensembles.
  \emph{Biostatistics}, \bold{7}(3), 355--373.

Hothorn T, Hornik K, Zeileis A (2006b).
  Unbiased Recursive Partitioning: A Conditional Inference Framework.
  \emph{Journal of Computational and Graphical Statistics}, \bold{15}(3), 651--674.

Hothorn T, Zeileis A (2015).
  partykit: A Modular Toolkit for Recursive Partytioning in R.
  \emph{Journal of Machine Learning Research}, \bold{16}, 3905--3909.

Meinshausen N (2006).
  Quantile Regression Forests.
  \emph{Journal of Machine Learning Research}, \bold{7}, 983--999.

Schlosser L, Hothorn T, Stauffer R, Zeileis A (2019).
  Distributional Regression Forests for Probabilistic Precipitation Forecasting in Complex Terrain.
  \emph{arXiv 1804.02921}, arXiv.org E-Print Archive.
  \url{http://arxiv.org/abs/1804.02921v3}

Strobl C, Boulesteix AL, Zeileis A, Hothorn T (2007).
  Bias in Random Forest Variable Importance Measures: Illustrations, Sources and a Solution.
  \emph{BMC Bioinformatics}, \bold{8}, 25.
  \url{http://www.biomedcentral.com/1471-2105/8/25}

Strobl C, Malley J, Tutz G (2009).
  An Introduction to Recursive Partitioning: Rationale, Application, and Characteristics of
  Classification and Regression Trees, Bagging, and Random Forests.
  \emph{Psychological Methods}, \bold{14}(4), 323--348.
}

\examples{
## basic example: distributional regression forest for cars data
df <- distforest(dist ~ speed, data = cars)

## prediction of fitted mean and visualization
nd <- data.frame(speed = 4:25)
nd$mean  <- predict(df, newdata = nd, type = "response")[["(fitted.response)"]]
plot(dist ~ speed, data = cars)
lines(mean ~ speed, data = nd)

}

\keyword{random forests, distributional regression trees, parametric modeling}
