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
\alias{summary.distfit}
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
  \item{y}{numeric vector.}
  \item{family}{specification of the response distribution.
    Either a \code{\link[gamlss.dist]{gamlss.family}} object or
    family list (see below for details, FIX ME).}
  \item{weights}{optional numeric vector of case weights.}
  \item{start}{starting values for the distribution parameters handed over to \code{optim}}
  \item{start.eta}{{starting values for the distribution parameters on the link scale handed over to \code{optim}}}
  \item{vcov}{logical or character. FIXME}
  \item{type.hessian}{FIXME}
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
  \item{\dots}{further arguments passed to \code{control} for
    \code{\link[stats]{optim}}.}
}

\details{
  The function \code{distfit} fits distributions,
  similar to \code{fitdistr} from MASS (Venables and Ripley 2002)
  but based on GAMLSS families (Stasinopoulos and Rigby 2007).
  
  Provides analytical gradients and hessian, can be plugged into
  \code{\link[partykit]{mob}}.

  The resulting object of class \code{"distfit"} comes with a set of 
  standard methods to generic functions including \code{coef}, \code{estfun}, \code{vcov}, \code{predict} and \code{logLik}. FIX ME
}

\value{
  \code{distfit} returns an object of class \code{"distfit"} which is a list with
  the following components:
  \item{npar}{npar,}
  \item{y}{y,}
  \item{ny}{ny,}
  \item{weights}{weights,}
  \item{family}{family,}
  \item{familylist}{familylist,}
  \item{start}{start,}
  \item{starteta}{starteta,}
  \item{opt}{opt,}
  \item{converged}{converged,}
  \item{par}{par,}
  \item{eta}{eta,}
  \item{hess}{hessian,}
  \item{vcov}{vc,}
  \item{loglik}{ll,}
  \item{call}{cl,}
  \item{estfun}{ef,}
  \item{ddist}{ddist,}
  \item{pdist}{pdist,}
  \item{qdist}{qdist,}
  \item{rdist}{rdist.}
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
    coef(m2)
    #coef(m3)
    logLik(m1)
    logLik(m2)
    #logLik(m3)
}
}

\keyword{distribution}
