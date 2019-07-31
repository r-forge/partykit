\name{circfit_control}
\alias{circfit_control}

\title{Auxiliary Function for Controlling Circular Regression Tree Fitting,}

\description{
  Auxiliary function for \code{circfit} fitting. Specifies a list of values passed
  to \code{\link{optim}}.
}
\usage{
circfit_control(solve_kappa = "Newton-Fourier", useC = FALSE, ncores = 1, ...)
}
\arguments{
  \item{solve_kappa}{Which kappa solver should be used for the starting values for kappa.
    By default a "Newton-Fourier" is used. Alternatively,
    a "Uniroot" provides a safe option and "Banerjee_et_al_2005"
    provides a quick approximation.}
  \item{useC}{Should score function and solver be calculated in C?}
  \item{ncores}{If useC = TRUE, number of cores for parallelization with openMP.}
  \item{\dots}{additional parameters passed to \code{\link{optim}}. }
}
\value{
  A list with components named as the arguments.
}

\seealso{\code{\link{circfit}}, \code{\link{optim}}}

\keyword{distribution}
