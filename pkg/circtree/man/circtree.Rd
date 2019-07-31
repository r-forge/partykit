\name{circtree}
\alias{circtree}
\alias{print.circtree}
\alias{predict.circtree}
\alias{coef.circtree}
\alias{logLik.circtree}
\alias{vcov.circtree}
\alias{estfun.circtree}
\alias{bread.circtree}
\alias{nobs.circtree}

\title{Regression Tree for a Circular Response}

\description{
  Trees based maximum-likelihood estimation of parameters for
  a circular response employing a von Mises distribution.
}

\usage{
circtree(formula, data, start, subset, na.action, weights, offset,
  control = partykit::mob_control(), fit_control = circfit_control(...), ...)
}

\arguments{
  \item{formula}{A symbolic description of the model to be fit. This
    should be of type \code{y ~ x1 + x2}
    where \code{y} should be the response variable
    and \code{x1} and \code{x2} are used as partitioning variables.}
  \item{data}{An optional data frame containing the variables in the model.}
  \item{start}{Starting values of distribution parameters used in the optimization. \strong{Currently, not supported.}}
  \item{subset}{An optional vector specifying a subset of observations to be
    used for fitting. \strong{Currently, not supported.}}
  \item{na.action}{A function which indicates what should happen when the data
    contain \code{NA}s. \strong{Currently, not supported.}}
  \item{weights}{Optional numeric vector of case weights.}
  \item{offset}{Optional numeric vector with \emph{a priori} known component to
    be included in the linear predictor for the location. \strong{Currently, not supported.}}
  \item{control}{Control arguments passed to \code{\link[partykit]{mob}}.}
  \item{fit_control}{A list of control parameters passed to
    \code{\link[stats]{optim}} in \code{\link{circfit}}. Default is \code{\link{circfit_control}()}}
  \item{\dots}{Arguments to be used to form the default \code{control} argument
    if it is not supplied directly.}
}

\details{
  Regression trees employing a von Mises distribution is an application of model-based recursive partitioning 
    and unbiased recursive partitioning based on the implementation in \code{\link[partykit]{mob}}.
}

\value{
  An object of S3 class \code{circtree} inheriting from class \code{modelparty} or \code{constparty}.
}

\examples{
sdat <- circtree_simulate()

m1.circtree <- circtree(y ~ x1 + x2, data = sdat)
#plot(m1.circtree, tp_args = list(kernel_density = TRUE)) 
#plot(m1.circtree) ## fails with cran version
}


\seealso{\code{\link[partykit]{mob}}, \code{\link{circfit}}}

\keyword{tree, parametric modelling}
