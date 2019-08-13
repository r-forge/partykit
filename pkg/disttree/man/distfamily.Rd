\name{distfamily}
\alias{distfamily}


\title{Preparation of family object of class \code{disttree.family} as employed in \code{distfit}, \code{disttree}, and \code{distforest}}

\description{
  The function \code{distfamily} prepares the required family object that is employed within
  \code{distfit} to estimate the parameters of the specified distribution family.
}

\usage{
  distfamily(family, bd = NULL, censpoint = NULL) 
} 

\arguments{
  \item{family}{can be one of the following: 
  \code{\link[gamlss.dist]{gamlss.family}} object,
  \code{\link[gamlss.dist]{gamlss.family}} function,
  character string with the name of a \code{\link[gamlss.dist]{gamlss.family}} object,
  function generating a family object with the required information about the distribution,
  character string with the name of a function generating a family object with the required information about the distribution, list with the required information about the distribution,
  character string with the name of a distribution for which a family  generating function is provided in \code{\link[disttree]{disttree}}
}
  \item{bd}{optional argument for binomial distributions specifying the binomial denominator}
  \item{censpoint}{censoring point for a censored \code{\link[gamlss.dist]{gamlss.family}} object}
}

\details{
The function \code{distfamily} is applied within \code{distfit}, \code{disttree}, and \code{distforest}. It generates a family object of class \code{disttree.family}.
If \code{family} is a \code{\link[gamlss.dist]{gamlss.family}} object the function \code{distfamily_gamlss} is called within \code{distfamily}.
}

\value{
  \code{distfamily} returns a family object of class \code{disttree.family} in form of a list with
  the following components:
  \item{family.name }{character string with the name of the specified distribution family}
  \item{ddist}{density function of the specified distribution family.}
  \item{sdist}{score function (1st partial derivatives) of the specified distribution family.}
  \item{hdist}{hessian function (2nd partial derivatives) of the specified distribution family.}
  \item{pdist}{distribution function of the specified distribution family.}
  \item{qdist}{quantile function of the specified distribution family.}
  \item{rdist}{random generation function of the specified distribution family.}
  \item{link}{character strings of the applied link functions.}
  \item{linkfun}{link functions.}
  \item{linkinv}{inverse link functions.}
  \item{linkinvdr}{derivative of the inverse link functions.}
  \item{startfun}{function generating the starting values for the employed optimization.}	
  \item{mle}{logical. Indicates whether a closed form solution exists (TRUE) for the maximum-likelihood optimization or whether a numerical optimization should be employed to estimate parameters (FALSE).}	
  \item{gamlssobj}{logical. Indicates whether the family has been obtained from a \code{\link[gamlss.dist]{gamlss.family}} object.}	
  \item{censored}{logical. Indicates whether the specified distribution family is censored.}	  
  \item{censpoint}{numeric. Censoring point (only if censored and gamlssobj),}
  \item{censtype}{character. Type of censoring ("left", "right") (only if censored and gamlssobj).}
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

\seealso{\code{disttree.family}, \code{\link[gamlss.dist]{gamlss.family}}}

\examples{
library(disttree)
family <- distfamily(family = NO())
}

\keyword{distribution, family}
