\name{circfit}
\alias{circfit}
\alias{bread.circfit}
\alias{coef.circfit}
\alias{estfun.circfit}
\alias{logLik.circfit}
\alias{nobs.circfit}
\alias{vcov.circfit}
\alias{predict.circfit}
\alias{print.circfit}

\title{Maximum-Likelihood Fitting of the von Mises Distribution.}

\description{
  The function \code{circfit} carries out maximum-likelihood
  estimation of parameters for the von Mises distribution.
}

\usage{
  circfit(y, x = NULL, start = NULL, subset = NULL, na.action = NULL,
    weights = NULL, offset = NULL, ...,
    vcov = TRUE, estfun = TRUE, object = FALSE, fit_control = circfit_control())
}

\arguments{
  \item{y}{Numeric vector of the response.}
  \item{x}{Covariates not supported. Use \code{\link{circmax}} instead.}
  \item{start}{Starting values of distribution parameters used in the optimization. \strong{Currently, not supported.}}
  \item{subset}{An optional vector specifying a subset of observations to be
    used for fitting. \strong{Currently, not supported.}}
  \item{na.action}{A function which indicates what should happen when the data
    contain \code{NA}s. \strong{Currently, not supported.}}
  \item{weights}{Optional numeric vector of case weights.}
  \item{offset}{Optional numeric vector with \emph{a priori} known component to
    be included in the linear predictor for the location. \strong{Currently, not supported.}}
  \item{vcov}{Logical. Specifies whether or not a variance-covariance matrix should be calculated and returned.}
  \item{estfun}{Logical. Should the matrix of observation-wise
    score contributions be returned?}
  \item{object}{Logical. If true, the fitted model is returned.}
  \item{fit_control}{A list of control parameters passed to
    \code{\link[stats]{optim}} in \code{\link{circfit}}. Default is \code{\link{circfit_control}()}}
  \item{\dots}{Further arguments passed to \code{\link[partykit]{mob}}. \strong{Really?}}
}

\details{
  The function \code{circfit} fits a von Mises distribution to the data.
  The resulting object of class \code{circfit} comes with a set of 
  standard methods to generic functions including \code{coef}, \code{estfun}, \code{vcov}, \code{predict} and \code{logLik}.
}

\seealso{\code{\link{circtree}}, \code{\link[stats]{optim}}}

\keyword{distribution}
