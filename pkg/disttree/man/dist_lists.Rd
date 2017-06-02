\name{dist_lists}
\alias{dist_gaussian}
\alias{dist_exponential}
\alias{dist_weibull}
\alias{dist_gamma}
\alias{dist_poisson}

\title{Family List Generating Functions}

\description{
 The functions \code{dist_gaussian},  \code{exponential}, \code{dist_weibull}, \code{dist_gamma} and \code{dist_poisson} g   generate a distribution family list with all the required elements to fit a distribution in \code{distfit}. 
}

\usage{
  dist_gaussian()
  dist_exponential()
  dist_weibull()
  dist_gamma()
  dist_poisson()
}

\arguments{}

\details{
  The functions \code{dist_gaussian},  \code{exponential}, \code{dist_weibull}, \code{dist_gamma} and \code{dist_poisson} g   generate a distribution family list with all the required elements to fit a distribution in \code{distfit}. These lists include a density function, a score function, a hessian function, starting values, link functions and inverse link functions.
}

\value{
These functions return a list with functions of the corresponding distribution family as required by \code{distfit}.
}


\examples{
## get the family list for a Gaussian distribution family
dist_gaussian()
}

\keyword{distribution family list}
