\name{circtree}
\alias{circtree}
\alias{print.circtree}
\alias{predict.circtree}
\alias{coef.circtree}
\alias{logLik.circtree}

\title{Distributional Regression Tree for a Circular Response}

\description{
  Distributional trees based on maximum-likelihood estimation of parameters for
  a circular response employing the von Mises distribution.
}

\usage{
circtree(formula, data, response_range = NULL, subset, na.action = na.pass,
         weights, offset, cluster, control = disttree::distextree_control(...),
         converged = NULL, scores = NULL, doFit = TRUE, ...)
}

\arguments{
  \item{formula}{A symbolic description of the model to be fit. This
    should be of type \code{y ~ x1 + x2}
    where \code{y} should be the response variable
    and \code{x1} and \code{x2} are used as partitioning variables.}
  \item{data}{An optional data frame containing the variables in the model.}
  \item{response_range}{Circular range of response variable.}
  \item{subset}{An optional vector specifying a subset of observations to be
    used for fitting.}
  \item{na.action}{A function which indicates what should happen when the data
    contain \code{NA}s.}
  \item{weights}{Optional numeric vector of case weights.}
  \item{offset}{Optional numeric vector with \emph{a priori} known component to
    be included in the linear predictor for the location. \strong{Currently, not supported.}}
  \item{cluster}{An optional vector (typically numeric or factor) with a
    cluster ID to be employed for clustered covariances in the parameter
    stability tests.}
  \item{control}{Control arguments passed to \code{\link[partykit]{extree_fit}} via \code{\link[disttree]{distextree_control}}.}
  \item{converged}{An optional function for checking user-defined criteria before splits are implemented.}
  \item{scores}{An optional named list of scores to be attached to ordered factors.}
  \item{doFit}{A logical indicating if the tree shall be grown (TRUE) or not FALSE}
  \item{\dots}{Arguments to be used to form the default \code{control} argument
    if it is not supplied directly.}
}

\details{
  Distributional regression trees for a circular response is an application of model-based recursive 
  partitioning and unbiased recursive partitioning based on the implementation in 
  \code{\link[disttree]{distextree}} using the infrastructure of \code{\link[partykit]{extree_fit}}.
}

\value{
  An object of S3 class \code{circtree} inheriting from class \code{distextree}.
}

\examples{
sdat <- circtree_simulate()

m1.circtree <- circtree(y ~ x1 + x2, data = sdat)
#plot(m1.circtree, tp_args = list(kernel_density = TRUE)) 
#plot(m1.circtree) ## fails with cran version
}


\seealso{
  \code{\link[disttree]{distextree}}, \code{\link[partykit]{mob}}, \code{\link[partykit]{ctree}}
}

\keyword{tree, parametric modelling, circular response}
