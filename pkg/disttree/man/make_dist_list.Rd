\name{make_dist_list}
\alias{make_dist_list}


\title{Family list generating function applied in distfit}

\description{
  The function \code{make_dist_list} generates a list with all the necessary elements to fit a parametric model in \code{distfit}
}

\usage{make_dist_list <- function(family, bd = NULL) }

\arguments{
  \item{family}{a \code{\link[gamlss.dist]{gamlss.family}} object}
  \item{bd}{binomial denominator}
}

\details{
The function \code{make_dist_list} is applied in distfit in case the family object is a \code{\link[gamlss.dist]{gamlss.family}} object.
  The required information and functions are extracted, built and structured such that the returned list can be used to fit a parametric model in \code{distfit}.
}

\value{
  \code{make_dist_list} returns a list with
  the following components:
  \item{family.name }{family.name}
  \item{ddist}{ddist,}
  \item{sdist}{sdist,}
  \item{hdist}{hdist,}
  \item{pdist}{pdist,}
  \item{qdist}{qdist,}
  \item{rdist}{rdist,}
  \item{link}{link,}
  \item{linkfun}{linkfun,}
  \item{linkinv}{linkinv,}
  \item{linkinvdr}{linkinvdr,}
  \item{startfun}{startfun,}	
  \item{mle}{mle,}	
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

\seealso{\code{\link[gamlss.dist]{gamlss.family}}}

\examples{
library(gamlss.dist)
familylist <- make_dist_list(family = NO())
}

\keyword{distribution, family}
