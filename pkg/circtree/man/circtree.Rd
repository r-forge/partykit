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
         weights, offset, cluster, control = disttree_control(...),
         converged = NULL, scores = NULL, doFit = TRUE, ...)
}

\arguments{
  \item{formula}{a symbolic description of the model to be fit. This
    should be of type \code{y ~ x1 + x2}
    where \code{y} should be the response variable
    and \code{x1} and \code{x2} are used as partitioning variables.}
  \item{data}{an optional data frame containing the variables in the model.}
  \item{response_range}{either a logical value indicating whether the response should be transformed
  to its original range (TRUE) or kept on the interval (-pi,pi] or a two-dimensional vector specifying 
  a range of the circular response.}
  \item{subset}{an optional vector specifying a subset of observations to be
    used in the fitting process.}
  \item{na.action}{a function which indicates what should happen when the data
    contain missing value.}
  \item{weights}{optional numeric vector of case weights.}
  \item{offset}{an optional vector of offset values.}
  \item{cluster}{an optional factor indicating independent clusters.
    Highly experimental, use at your own risk.}
  \item{control}{control arguments passed to \code{\link[partykit]{extree_fit}} via \code{\link[disttree]{disttree_control}}.}
  \item{converged}{an optional function for checking user-defined criteria before splits are implemented.}
  \item{scores}{an optional named list of scores to be attached to ordered factors.}
  \item{doFit}{a logical indicating if the tree shall be grown (TRUE) or not FALSE}
  \item{\dots}{arguments to be used to form the default \code{control} argument
    if it is not supplied directly.}
}

\details{
  Distributional regression trees for a circular response are an application of model-based recursive 
  partitioning and unbiased recursive partitioning based on the implementation in 
  \code{\link[disttree]{disttree}} using the infrastructure of \code{\link[partykit]{extree_fit}}.
}

\value{
  An object of S3 class \code{circtree} inheriting from class \code{disttree}.
}

\examples{

## example on parameter range:
sdat.par <- circtree_simulate(response_range = c(-pi, pi))
ct.par <- circtree(y ~ x1 + x2, data = sdat.par)
plot(ct.par)

## example on response range (0, 2pi):
sdat.rad <- circtree_simulate(response_range = c(0, 2*pi))
ct.rad <- circtree(y ~ x1 + x2, data = sdat.rad)
## default: type = "response"
plot(ct.rad, tp_args = list(type = "parameter"))
plot(ct.rad, tp_args = list(type = "response"))

## example on response range (0, 360):
sdat.deg <- circtree_simulate(response_range = c(0, 360))
ct.deg <- circtree(y ~ x1 + x2, data = sdat.deg)
plot(ct.deg, tp_args = list(type = "parameter"))
plot(ct.deg, tp_args = list(type = "response"))
plot(ct.deg, tp_args = list(plot_type = "geographics"))

## example on response range (0, 24):
sdat.hour <- circtree_simulate(response_range = c(0, 24))
ct.hour <- circtree(y ~ x1 + x2, data = sdat.hour, response_range = c(0, 24))
plot(ct.hour, tp_args = list(type = "parameter"))
plot(ct.hour, tp_args = list(type = "response", plot_type = "time"))
plot(ct.hour, tp_args = list(type = "response", plot_type = "time", 
  circlab = c("no", "mo", "mi", "ev")))
}


\seealso{
  \code{\link[disttree]{disttree}}, \code{\link[disttree]{distfit}}, 
  \code{\link[partykit]{extree_fit}}
}

\keyword{regression tree, parametric modelling, circular response}
