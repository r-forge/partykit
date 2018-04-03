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
disttree(formula, data, na.action, cluster, family = NO(), bd = NULL,
         type.tree = "mob", decorrelate = "none", offset,
         censtype = "none", censpoint = NULL, weights = NULL,
         terminal_objects = FALSE,
         control = partykit::mob_control(), ocontrol = list(), ...)
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
  \item{censpoint}{Numeric value. Censoring point can be set for censored response.}
  \item{weights}{Optional numeric vector of case weights.}
  \item{terminal_objects}{Logical argument deciding whether a model object should be fit to the terminal node.}
  \item{control}{Control arguments passed to \code{\link[partykit]{mob}} or \code{\link[partykit]{ctree}}.}
  \item{ocontrol}{List with control parameters passed to
    \code{\link[stats]{optim}} in \code{\link[disttree]{distfit}}.}
  \item{\dots}{Further arguments passed to \code{\link[stats]{optim}} in \code{\link[disttree]{distfit}}.}
}

\details{
  Distributional regression trees are an application of model-based recursive partitioning and unbiased recursive partitioning
  (implemented in \code{\link[partykit]{mob}} and \code{\link[partykit]{ctree}}) to parametric model fits
  based on the GAMLSS family of distribtuions.
}

\value{
  An object of S3 class \code{disttree} inheriting from class \code{modelparty} or \code{constparty}.
}

\seealso{\code{\link[partykit]{mob}}, \code{\link[partykit]{ctree}}, \code{\link[disttree]{distfit}}}

\examples{
tr <- disttree(dist ~ speed, data = cars)
print(tr)

plot(tr)
plot(as.constparty(tr))
}

\keyword{tree, parametric modelling}
