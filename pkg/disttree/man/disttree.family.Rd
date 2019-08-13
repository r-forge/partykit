\name{disttree.family}
\alias{family.disttree}
\alias{dist_gaussian}
\alias{dist_exponential}
\alias{dist_weibull}
\alias{dist_gamma}
\alias{dist_poisson}
\alias{dist_crch}
\alias{dist_list_normal}
\alias{dist_list_cens_normal}
\alias{dist_binomial}
\alias{dist_list_hurdle_normal}
\alias{dist_list_trunc_normal}
\alias{dist_ztnbinom}

\title{Family List Generating Functions}

\description{
 The functions \code{dist_gaussian}, \code{dist_crch},  \code{dist_exponential}, \code{dist_weibull}, \code{dist_gamma} and \code{dist_poisson} generate a distribution family object of class \code{disttree.family} with all the required elements to fit a distribution in \code{distfit}. 
 
 Complete distribution family lists are provided for example by \code{dist_list_normal} and \code{dist_list_cens_normal}.
}

\usage{
  dist_gaussian()
  dist_crch(dist = c("gaussian", "logistic"), truncated = FALSE,
            type = c("left", "right", "interval"), censpoint = 0)
  dist_exponential()
  dist_weibull()
  dist_gamma()
  dist_poisson()
}

\arguments{
\item{dist}{\code{character}. Either a gaussian ('gaussian') or a logistic ('logistic') distribution can be selected.}
\item{truncated}{\code{logical}. If TRUE truncated family list is generated with 'censpoint' interpreted as truncation points, If FALSE censored family list is generated. Default is FALSE
}
\item{type}{\code{character}. Type of censoring can be selectes ('left', 'right' or 'interval')}
\item{censpoint}{\code{numeric}. Censoring point can be set (per default set to 0).}
}

\details{
  The functions \code{dist_gaussian}, \code{dist_crch},  \code{dist_exponential}, \code{dist_weibull}, \code{dist_gamma} and \code{dist_poisson} generate a distribution family list with all the required elements to fit a distribution in \code{distfit}. These lists include a density function, a score function, a hessian function, starting values, link functions and inverse link functions.
  
  Complete distribution family lists are provided for example by \code{dist_list_normal} and \code{dist_list_cens_normal} for the normal and censored normal distribution respectively.
}

\value{
These functions return a family of class \code{disttree.family} with functions of the corresponding distribution family as required by \code{distfit}, \code{disttree}, and \code{distforest}.
}

\seealso{\code{\link[disttree]{distfamily}}}



\examples{
## get the family list for a Gaussian distribution family
dist_gaussian()
}

\keyword{distribution family}
