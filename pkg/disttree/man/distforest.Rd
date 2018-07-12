\name{distforest}

\alias{distforest}
\alias{predict.distforest}
\alias{logLik.distforest}

\title{Distributional Regression Forests}

\description{
  Forests based on maximum-likelihood estimation of parameters for
  distributions from the GAMLSS family (for generalized additive
  models for location, scale, and shape).
}

\usage{
distforest(formula, data, na.action = na.pass, cluster, family = NO(), bd = NULL,
           type.tree = "ctree", decorrelate = "none", offset,
           censtype = "none", censpoint = NULL, weights = NULL,
           control = partykit::ctree_control(teststat = "quad", testtype = "Univ", mincriterion = 0, ...),
           ocontrol = list(), type.hessian = c("checklist", "analytic", "numeric"),
           ntree = 500L, fit = TRUE, perturb = list(replace = FALSE, fraction = 0.632), fitted.OOB = TRUE,
           cores = NULL, applyfun = NULL,
           mtry = ceiling(sqrt(nvar)),
           ...)        
}

\arguments{
  \item{formula}{A symbolic description of the model to be fit. This
    should be of type \code{y ~ x1 + x2}
    where \code{y} should be the response variable
    and \code{x1} and \code{x2} are used as partitioning variables.}
  \item{data}{An optional data frame containing the variables in the model.}
  \item{na.action}{A function which indicates what should happen when the data
    contain \code{NA}s.}
  \item{cluster}{An optional vector (typically numeric or factor) with a
    cluster ID to be employed for clustered covariances in the parameter
    stability tests.}
  \item{family}{specification of the response distribution.
    Either a \code{\link[gamlss.dist]{gamlss.family}} object, a list generating function or a family list.}
  \item{bd}{binomial denominator: additional parameter needed for binomial gamlss.families}
  \item{type.tree}{Specification of the type of tree learner, either
    \code{"mob"} or \code{"ctree"}.}
  \item{decorrelate}{Specification of the type of decorrelation for the
    empirical estimating functions (or scores) either \code{"none"} or
    \code{"opg"} (for the outer product of gradients) or \code{"vcov"}
    (for the variance-covariance matrix, assuming this is an estimate
    of the Fisher information).}
  \item{offset}{FIX ME.}
  \item{censtype}{Can either be 'none', 'left' or 'right' to set the type of censoring for censored response.}
  \item{censpoint}{numeric value. Censoring point can be set for censored response.}
  \item{weights}{optional numeric vector of case weights.}
  \item{control}{Control arguments passed to \code{\link[partykit]{mob}} or \code{\link[partykit]{ctree}}.}
  \item{ocontrol}{List with control parameters passed to
    \code{\link[stats]{optim}} in \code{distfit}.}
  \item{type.hessian}{Can either be 'checklist', 'analytic' or 'numeric' to decide how the hessian matrix should be calculated in the fitting process in \code{distfit}. For 'checklist' it is checked whether a function 'hdist' is given in the family list. If so, 'type.hessian' is set to 'analytic', otherwise to 'numeric'.}
  \item{ntree}{Number of trees to grow for the forest.}
  \item{fit}{logical. if TRUE, fitted and predicted values and  predicted parameters are calculated for the learning data (together with loglikelihood)}
  \item{perturb}{a list with arguments \code{replace} and \code{fraction} determining which type of
    resampling with \code{replace = TRUE} referring to the n-out-of-n bootstrap and
    \code{replace = FALSE} to sample splitting. \code{fraction} is the number of observations to draw without
    replacement. }
  \item{fitted.OOB}{logical. if fitted.OOB=TRUE the weights for each observation of the learning data are predicted by \code{\link[partykit]{predict.cforest}} with the argument OOB=TRUE (only relevant if fit=TRUE)}
  \item{cores}{numeric. If set to an integer the \code{applyfun} is set to
               \code{\link[parallel]{mclapply}} with the desired number of \code{cores}.}
  \item{applyfun}{an optional \code{\link[base]{lapply}}-style function with arguments
                  \code{function(X, FUN, \dots)}. It is used for computing the variable selection criterion.
                  The default is to use the basic \code{lapply}
                  function unless the \code{cores} argument is specified (see below).}
  \item{mtry}{number of input variables randomly sampled as candidates
    at each node for random forest like algorithms. Bagging, as special case
    of a random forest without random input variable sampling, can
    be performed by setting \code{mtry} either equal to \code{Inf} or
    manually equal to the number of input variables.}
  \item{\dots}{further arguments passed to \code{\link[stats]{optim}} in \code{\link[disttree]{distfit}}.}
}

\details{
  Distributional regression forests are an application of model-based recursive partitioning
  (implemented in \code{\link[partykit]{mob}}, \code{\link[partykit]{ctree}} and \code{\link[partykit]{cforest}}) to parametric model fits based on the GAMLSS family of distribtuions.
}

\value{
  An object of S3 class \code{distforest} inheriting from class \code{cforest}.
}

\seealso{\code{\link[partykit]{mob}}, \code{\link[partykit]{ctree}}, \code{\link[partykit]{cforest}}, \code{\link{distfit}}}

\examples{
df <- distforest(dist ~ speed, data = cars)
predict(df)
}

\keyword{regression trees, random forests}
