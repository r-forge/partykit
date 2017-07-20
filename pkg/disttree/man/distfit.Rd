\name{distfit}
\alias{distfit}
\alias{bread.distfit}
\alias{coef.distfit}
\alias{confint.distfit}
\alias{estfun.distfit}
\alias{logLik.distfit}
\alias{nobs.distfit}
\alias{vcov.distfit}
\alias{predict.distfit}
\alias{print.distfit}
\alias{summary.distfit}
\alias{print.summary.distfit}
\alias{getSummary.distfit}
\alias{residuals.distfit}
\alias{plot.distfit}

\title{Maximum-Likelihood Fitting of Parametric Distributions}

\description{
  The function \code{distfit} carries out maximum-likelihood
  estimation of parameters for distributions from the GAMLSS
  family (for generalized additive models for location, scale,
  and shape). The parameters can be transformed through
  link functions but do not depend on further covariates
  (i.e., are constant across observations).
}

\usage{
distfit(y, family, weights = NULL, start = NULL, start.eta = NULL, 
        vcov = TRUE, type.hessian = "analytic", estfun = TRUE, 
        bd = NULL, fixed = NULL, fixed.values = NULL, 
        cens = "none", censpoint = NULL, ocontrol = list(), ...)
}

\arguments{
  \item{y}{numeric vector of the response}
  \item{family}{specification of the response distribution.
    Either a \code{\link[gamlss.dist]{gamlss.family}} object, a list generating function or a family list.}
  \item{weights}{optional numeric vector of case weights.}
  \item{start}{starting values for the distribution parameters handed over to \code{\link[stats]{optim}}}
  \item{start.eta}{starting values for the distribution parameters on the link scale handed over to \code{\link[stats]{optim}}.}
  \item{vcov}{logical. Specifies whether or not a variance-covariance matrix should be calculated and returned.}
  \item{type.hessian}{character. Can be "analytic" or "numeric" and specifies how the hessian matrix should be calculated.}
  \item{estfun}{logical. Should the matrix of observation-wise
    score contributions (or empirical estimating functions) be
    returned?}
  \item{bd}{binomial denominator}
  \item{fixed}{FIXME}
  \item{fixed.values}{FIXME}
  \item{cens}{Can either be 'none', 'left' or 'right' to set the type of censoring for censored response.}
  \item{censpoint}{numeric value. Censoring point can be set for censored response.}
  \item{ocontrol}{List with control parameters passed to
    \code{\link[stats]{optim}}.}
  \item{\dots}{further arguments passed to \code{\link[stats]{optim}}.}
}

\details{
  The function \code{distfit} fits distributions,
  similar to \code{\link[MASS]{fitdistr}} from \pkg{MASS} (Venables and Ripley 2002)
  but based on GAMLSS families (Stasinopoulos and Rigby 2007).
  
  Provides analytical gradients and hessian, can be plugged into
  \code{\link[partykit]{mob}}.

  The resulting object of class \code{distfit} comes with a set of 
  standard methods to generic functions including \code{coef}, \code{estfun}, \code{vcov}, \code{predict} and \code{logLik}.
}

\value{
  \code{distfit} returns an object of class \code{distfit} which is a list with
  the following components:
  \item{npar}{number of parameter}
  \item{y}{numeric vector of the response}
  \item{ny}{number of observations}
  \item{weights}{numeric vector of case weights handed over as input argument}
  \item{family}{name of the distribution family}
  \item{familylist}{if the input argument 'family' was already a list, this list is returned. Otherwise the generated family list is returned.}
  \item{start}{used starting values in \code{\link[stats]{optim}} that were handed over as input argument}
  \item{starteta}{starting value on the link scale used in \code{\link[stats]{optim}}}
  \item{opt}{list returned by \code{\link[stats]{optim}}}
  \item{converged}{logical. TRUE if \code{\link[stats]{optim}} returnes convergence = 0 and FALSE else.}
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
  
## simple distfit
df <- distfit(y, family = NBI)
coef(df)
confint(df)
logLik(df)

## using tabulated data
ytab <- table(y)
df2 <- distfit(as.numeric(names(ytab)), family = NBI, weights = ytab)
coef(df2)
confint(df2)
logLik(df2)

## coefficients tests
if(require("lmtest")) {
  coeftest(df)
  coeftest(df2)
}

## censored logistic example
if(require("crch") & require("gamlss.cens")) {
    library("crch")
    data("RainIbk", package = "crch")
    m1 <- crch(rain ~ 1, data = RainIbk, left = 0, dist = "logistic")
  
    library("gamlss.cens")
    gen.cens(LO, type = "left")
    m2 <- distfit(RainIbk$rain, family = LOlc, cens = "left", censpoint = 0)
    
    #dist_list_cens_log <- make_censored_dist_list(dist = "logistic", type = "left", censpoint = 0)
    #m3 <- distfit(RainIbk$rain, family = dist_list_cens_log)
    
    coef(m1)
    coef(m2, type = "link")
    #coef(m3)
    logLik(m1)
    logLik(m2)
    #logLik(m3)
}
}

\keyword{distribution}
