\name{circfit}
\alias{circfit}

\title{Maximum-Likelihood Fitting for a Circular Response}

\description{
  The function \code{circfit} carries out maximum-likelihood estimation of 
  parameters for a circular response employing the von Mises distribution. 
  The parameters can be transformed through link functions but do not depend 
  on further covariates (i.e., are constant across observations).
}

\usage{
circfit(y, weights = NULL, start = NULL, start.eta = NULL,
        response_range = NULL,
        vcov = TRUE, type.hessian =  c("checklist", "analytic", "numeric"), 
        method = "L-BFGS-B", estfun = TRUE, optim.control = list(), ...)
}

\arguments{
  \item{y}{numeric vector of the response}
  \item{weights}{optional numeric vector of case weights.}
  \item{start}{starting values for the distribution parameters handed over to \code{\link[stats]{optim}}}
  \item{start.eta}{starting values for the distribution parameters on the link scale handed over to 
    \code{\link[stats]{optim}}.}
  \item{response_range}{either a logical value indicating whether the response should be transformed
  to its original range (TRUE) or kept on the interval (-pi,pi] or a two-dimensional vector specifying 
  a range of the circular response.}
  \item{vcov}{logical. Specifies whether or not a variance-covariance matrix should be calculated and returned.}
  \item{type.hessian}{Can either be 'checklist', 'analytic' or 'numeric' to decide how the hessian matrix should be 
    calculated in the fitting process in \code{distfit}. For 'checklist' it is checked whether a function 'hdist' 
    is given in the family list. If so, 'type.hessian' is set to 'analytic', otherwise to 'numeric'.}
  \item{method}{Optimization which should be applied in \code{\link[stats]{optim}}}
  \item{estfun}{logical. Should the matrix of observation-wise score contributions (or empirical estimating functions) be returned?}
  \item{optim.control}{A list with \code{\link[stats]{optim}} control parameters.}
  \item{\dots}{further arguments passed to \code{\link[stats]{optim}}.}

  
}

\details{
  The function \code{circfit} fits the parameter of the von Mises distribution to a
  circular response variable by applying \code{\link[disttree]{distfit}}.
}

\value{
  An object of S3 class \code{circfit} inheriting from class \code{distfit}.
}

\examples{

## example on parameter range:
sdat.par <- circtree_simulate(response_range = c(-pi, pi))
cf.par <- circfit(sdat.par$y)


## example on response range (0, 2pi):
sdat.rad <- circtree_simulate(response_range = c(0, 2*pi))
cf.rad <- circfit(sdat.rad$y)

## example on response range (0, 360):
sdat.deg <- circtree_simulate(response_range = c(0, 360))
cf.deg <- circfit(sdat.deg$y)
}


\seealso{
  \code{\link[disttree]{distfit}}
}

\keyword{maximum likelihood, von Mises, circular response}
