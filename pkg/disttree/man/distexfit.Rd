\name{distexfit}
\alias{distexfit}
\alias{bread.distexfit}
\alias{coef.distexfit}
\alias{confint.distexfit}
\alias{estfun.distexfit}
\alias{logLik.distexfit}
\alias{nobs.distexfit}
\alias{vcov.distexfit}
\alias{predict.distexfit}
\alias{print.distexfit}
\alias{summary.distexfit}
\alias{print.summary.distexfit}
\alias{getSummary.distexfit}
\alias{residuals.distexfit}
\alias{plot.distexfit}

\title{Maximum-Likelihood Fitting of Parametric Distributions}

\description{
  The function \code{distexfit} carries out maximum-likelihood
  estimation of parameters for a specified distribution family, 
  for example from the GAMLSS family (for generalized additive 
  models for location, scale, and shape). The parameters can be 
  transformed through link functions but do not depend on further 
  covariates (i.e., are constant across observations).
}

\usage{
distexfit(y, family = NO(), weights = NULL, start = NULL, start.eta = NULL,
          vcov = TRUE, type.hessian =  c("checklist", "analytic", "numeric"),
          method = "L-BFGS-B", estfun = TRUE, optim.control = list(), ...)
}

\arguments{
  \item{y}{numeric vector of the response}
  \item{family}{specification of the response distribution.
    Either a \code{\link[gamlss.dist]{gamlss.family}} object, a list generating function or a family list.}
  \item{weights}{optional numeric vector of case weights.}
  \item{start}{starting values for the distribution parameters handed over to \code{\link[stats]{optim}}}
  \item{start.eta}{starting values for the distribution parameters on the link scale handed over to \code{\link[stats]{optim}}.}
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
  The function \code{distexfit} fits distributions,
  similar to \code{\link[MASS]{fitdistr}} from \pkg{MASS} (Venables and Ripley 2002)
  but based on GAMLSS families (Stasinopoulos and Rigby 2007).
  
  Provides analytical gradients and hessian, can be plugged into
  \code{\link[partykit]{mob}}.

  The resulting object of class \code{distexfit} comes with a set of 
  standard methods to generic functions including \code{coef}, \code{estfun}, \code{vcov}, \code{predict} and \code{logLik}.
}

\value{
  \code{distexfit} returns an object of class \code{distexfit} which is a list with
  the following components:
  \item{npar}{number of parameter}
  \item{y}{numeric vector of the response}
  \item{ny}{number of observations}
  \item{weights}{numeric vector of case weights handed over as input argument}
  \item{family}{employed distribution family list of class \code{disttree.family}}
  \item{start}{used starting values in \code{\link[stats]{optim}} that were handed over as input argument}
  \item{starteta}{starting values on the link scale used in \code{\link[stats]{optim}}}
  \item{opt}{list returned by \code{\link[stats]{optim}}}
  \item{converged}{logical. TRUE if \code{\link[stats]{optim}} returns convergence = 0 and FALSE else.}
  \item{par}{fitted distribution parameters (on parameter scale)}
  \item{eta}{fitted distribution parameters (on link scale)}
  \item{hess}{hessian matrix}
  \item{vcov}{variance-covariance matrix}
  \item{loglik}{value of the maximized log-likelihood function}
  \item{call}{function call}
  \item{estfun}{matrix with the scores for the estimated parameters. Each line represents an observation and each column a parameter.}
  \item{ddist}{density function with the estimated distribution parameters already plugged in}
  \item{pdist}{probability function with the estimated distribution parameters already plugged in}
  \item{qdist}{quantile function with the estimated distribution parameters already plugged in}
  \item{rdist}{random number generating function with the estimated distribution parameters already plugged in}
  \item{method}{optimization method applied in \code{\link[stats]{optim}}}
}

\references{
Stasinopoulos DM, Rigby RA (2007). 
  Generalized Additive Models for Location Scale and Shape (GAMLSS) in R,
  \emph{Journal of Statistical Software}, \bold{23}(7), 1-46.
  \doi{10.18637/jss.v023.i07}

Venables WN, Ripley BD (2002).
  \emph{Modern Applied Statistics with S}.
  4th Edition. Springer-Verlag, New York.
}

\seealso{\code{\link[gamlss.dist]{gamlss.family}}, \code{\link[stats]{optim}}}

\examples{
## simulate artifical negative binomial data
set.seed(0)
y <- rnbinom(1000, size = 1, mu = 2)
  
## simple distexfit
df <- distexfit(y, family = NBI)

}

\keyword{distribution}
