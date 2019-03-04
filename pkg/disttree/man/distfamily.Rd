\name{distfamily}
\alias{distfamily}


\title{Preparation of family list as required in distfit}

\description{
  The function \code{distfamily} prepares the required family list in \code{distfit} from
  the handed over family object (which can already be a list as well).
}

\usage{
  distfamily(family) 
} 

\arguments{
  \item{family}{can be one of the following: 
  \code{\link[gamlss.dist]{gamlss.family}} object,
  \code{\link[gamlss.dist]{gamlss.family}} function,
  character string with the name of a \code{\link[gamlss.dist]{gamlss.family}} object,
  function generating a list with the required information about the distribution,
  character string with the name of a function generating a list with the required information about the distribution,
  list with the required information about the distribution,
  character string with the name of a distribution for which a list generating function is provided in \code{\link[disttree]{disttree}}
}
}

\details{
The function \code{distfamily} is applied in \code{distfit}.
Depending on the type (character, list, or a \code{\link[gamlss.dist]{gamlss.family}} object).
}

\value{
  \code{distfamily} returns a list with
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
  \item{gamlssobj}{gamlssobj,}	
  \item{censored}{censored}	
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
library(disttree)
family <- distfamily(family = NO())
}

\keyword{distribution, family}
