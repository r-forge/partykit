\name{disttree}

\alias{disttree}
\alias{print.disttree}
\alias{predict.disttree}
\alias{coef.disttree}
\alias{logLik.disttree}

\title{Distributional Regression Tree}

\description{
  Trees based on maximum-likelihood estimation of parameters for
  distributions from the GAMLSS family (for generalized additive
  models for location, scale, and shape).
}

\usage{
disttree(formula, data, na.action, cluster, family = NO(),
         type.tree = "mob", decorrelate = "none", offset,
         cens = "none", censpoint = NULL, weights = NULL,
         control = mob_control(...), ocontrol = list(), ...)
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
    Either a \code{\link[gamlss.dist]{gamlss.family}} object or
    family list (see below for details, FIX ME).}
  \item{type.tree}{Specification of the type of tree learner, either
    \code{"mob"} or \code{"ctree"}.}
  \item{decorrelate}{Specification of the type of decorrelation for the
    empirical estimating functions (or scores) either \code{"none"} or
    \code{"opg"} (for the outer product of gradients) or \code{"vcov"}
    (for the variance-covariance matrix, assuming this is an estimate
    of the Fisher information).}
  \item{offset}{FIX ME.}
  \item{cens}{Can either be 'none', 'left' or 'right' to set the type of censoring for censored response.}
  \item{censpoint}{numeric value. Censoring point can be set for censored response.}
  \item{weights}{optional numeric vector of case weights.}
  \item{control}{Arguments passed to \code{\link[partykit]{mob_control}}.}
  \item{ocontrol}{List with control parameters passed to
    \code{\link[stats]{optim}}.}
  \item{\dots}{further arguments passed to \code{control} for
    \code{\link[stats]{optim}}.}
}

\details{
  Distributional regression trees are an application of model-based recursive partitioning and unbiased recursive partitioning
  (implemented in \code{\link[partykit]{mob}} and \code{\link[partykit]{ctree}}) to parametric model fits
  based on the GAMLSS family of distribtuions.
}

\value{
  An object of S3 class \code{"disttree"} inheriting from class \code{"modelparty"} or \code{"constparty"}.
}

\seealso{\code{\link[partykit]{mob}}, \code{\link[partykit]{ctree}}, \code{\link{distfit}}}

\examples{
tr <- disttree(dist ~ speed, data = cars)
print(tr)

plot(tr)
plot(as.constparty(tr))
}

\keyword{tree, parametric modelling}
