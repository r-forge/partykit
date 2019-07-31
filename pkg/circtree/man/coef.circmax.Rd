\name{coef.circmax}
\alias{coef.circmax}
\alias{print.circmax}
\alias{logLik.circmax}
\alias{terms.circmax}
\alias{model.frame.circmax}
\alias{model.matrix.circmax}
\alias{estfun.circmax}
\alias{vcov.circmax}

\title{Methods for CIRCMAX Objects}

\description{
Methods for extracting information from fitted \code{circmax} objects.
}

\usage{
\method{coef}{circmax}(object, model = c("full", "location", "concentration"), \dots)
\method{terms}{circmax}(x, model = c("location", "concentration", "full"), \dots)
}

\arguments{
  \item{object, x}{an object of class \code{"circmax"}.}
  \item{model}{model for which coefficients shall be returned.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\details{
  In addition to the methods above, a set of standard extractor functions for
  \code{"circmax"} objects is available, including methods to the generic
  functions \code{\link[base]{print}}, \code{\link[stats]{logLik}}, 
  and \code{\link[stats]{model.frame}}. Additionally, \code{\link[sandwich]{estfun}},
  \code{\link[stats]{vcov}} provide methods for 'robust' inference.
}

\seealso{\code{\link{circmax}}}

\keyword{regression}

