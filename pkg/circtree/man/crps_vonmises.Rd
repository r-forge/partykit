\name{crps_vonmises}
\alias{crps_vonmises}

\title{Circular CRPS}

\description{
Continuous Ranked Probability Score (CRPS) for a circular response
following the von Mises distribution.
}

\usage{
crps_vonmises(y, mu, kappa, sum = FALSE)
}

\arguments{
  \item{y}{numeric. Circular response.}
  \item{mu}{numeric. Location parameter of the von Mises distribution.}
  \item{kappa}{numeric. Concentration parameter of the von Mises distribution.}
  \item{sum}{logical. Should the sum of the CRPS-values over all response values be returned.}
}

\value{
For \code{sum=TRUE} the sum of the CRPS-values of all response values is returned. Otherwise
a vector of the same length as \code{y} with the observation-wise CRPS-values is returned.
}

\examples{
set.seed(20190731)
y <- rnorm(100, 0, 1)
y[abs(y)>pi] <- 0
testfit <- circfit(y)
testcrps <- crps_vonmises(y, mu = as.numeric(testfit$coefficients[1]), 
kappa = as.numeric(testfit$coefficients[2]))
}


\seealso{\code{\link{circtree}}}

\keyword{CPRS, circular}
