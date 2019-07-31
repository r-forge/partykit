\name{circmax}
\alias{circmax}
\alias{circmax_fit}
\alias{circmax_control}

\title{Circular Regression with Maximum Likelihood Estimation}

\description{
  Fit a regression model for a circular response
    by maximum likelihood estimation employing the von Mises distribution.
}
\usage{
circmax(formula, data, subset, na.action,
  model = TRUE, y = TRUE, x = FALSE,
  control = circmax_control(\dots), \dots)

circmax_fit(x, y, z = NULL, control)

circmax_control(maxit = 5000, start = NULL, method = "Nelder-Mead",
  solve_kappa = "Newton-Fourier", 
  gradient = FALSE, hessian = TRUE, \dots)
}
\arguments{
  \item{formula}{a formula expression of the form \code{y ~ x | z} where
    \code{y} is the response and \code{x} and \code{z} are regressor variables
    for the location and the concentration of the von Mises distribution.}
  \item{data}{an optional data frame containing the variables occurring in the
    formulas; y has to be given in radians.}
  \item{subset}{an optional vector specifying a subset of observations to be
    used for fitting.}
  \item{na.action}{a function which indicates what should happen when the data
    contain \code{NA}s.}
  \item{model}{logical. If \code{TRUE} \emph{model frame} is
    included as a component of the returned value.}
  \item{x, y}{for \code{circmax}: logical. If \code{TRUE} the model matrix and
    response vector used for fitting are returned as components of the returned value.
    For \code{circmax_fit}: \code{x} is a design matrix with regressors for the
    location and \code{y} is a vector of observations given in radians.}
  \item{z}{a design matrix with regressors for the concentration.}
  \item{\dots}{arguments to be used to form the default \code{control} argument
    if it is not supplied directly.}
  \item{control, maxit, start}{a list of control parameters passed to \code{\link[stats]{optim}}.}
  \item{method}{The \code{method} to be used for optimization.}
  \item{solve_kappa}{Which kappa solver should be used for the starting values for kappa.
    By default a Newton Fourier is used (\code{"Newton-Fourier"}). Alternatively,
    a uniroot provides a safe option (\code{"Uniroot"}) or code{"Banerjee_et_al_2005"}
    provides a quick approximation).}
  \item{gradient}{logical. Should gradients be used for optimization? If \code{TRUE},
    the default \code{method} is \code{"BFGS"}. Otherwise \code{method = "Nelder-Mead"}
    is used.}
  \item{hessian}{logical or character. Should a numeric approximation of the
    (negative) Hessian matrix by \code{\link[stats]{optim}} be computed?}
}

\details{
\code{circmax} fits a regression model for a circular response assuming a von Mises distribution.


\code{circmax_fit} is the lower level function where the parameters of the von Mises distribution 
  are fitted by maximum likelihood estimation.
}

\value{
  An object of class \code{"circmax"}.
}

\examples{
## Example 1: Simulated Data:

sdat <- circmax_simulate(n = 1000, beta = c(3, 5, 2), gamma = c(3, 3))

(m1.circmax <- circmax(y ~ x1 + x2 | x3, data = sdat))

## Example 2: Periwinkle Dataset of Fisher and Lee, 1992:
require("circular")
distance <- c(107, 46, 33, 67, 122, 69, 43, 30, 12, 25, 37, 69, 5, 83, 
  68, 38, 21, 1, 71, 60, 71, 71, 57, 53, 38, 70, 7, 48, 7, 21, 27)
directdeg <- c(67, 66, 74, 61, 58, 60, 100, 89, 171, 166, 98, 60, 197, 
  98, 86, 123, 165, 133, 101, 105, 71, 84, 75, 98, 83, 71, 74, 91, 38, 200, 56)
cdirect <- circular(directdeg * 2 * pi/360)
plot(as.numeric(cdirect) ~ distance, ylim = c(0, 4*pi), pch = 20)
points(as.numeric(cdirect) + 2*pi ~ distance, pch = 20)

(m2.circ <- lm.circular(type = "c-l", y = cdirect, x = distance, init = 0.0))
(m2.circmax <- circmax(cdirect ~ distance, data = data.frame(cbind(distance, cdirect))))
}
\keyword{circular, regression, von Mises}
