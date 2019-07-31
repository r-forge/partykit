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
  set.seed(123)
  sdat <- circtree_simulate(n = 100)
  m1.circtree <- circtree(y ~ x1 + x2, data = sdat)
  foo <- function(x, deg = FALSE){
    if(deg) x <- x *  pi / 180
    tmp <- x %% (2 * pi)
    tmp <- ifelse(tmp > pi, -(pi - (tmp - pi)), tmp)
    if(deg) tmp <- tmp *  180 / pi
    tmp
  }
   
  testcrps <- crps_vonmises(foo(sdat$y), mu = predict(m1.circtree)$mu, predict(m1.circtree)$kappa) 
}


\seealso{\code{\link{circtree}}}

\keyword{CPRS, circular response}
