\name{disttree}

\alias{disttree}
\alias{print.disttree}
\alias{coef.disttree}
\alias{predict.disttree}
\alias{logLik.disttree}
\alias{fitted.disttree}


\title{Distributional Regression Tree}

\description{
  Trees based on maximum-likelihood estimation of parameters for
  specified distribution families, for example from the GAMLSS family 
  (for generalized additive models for location, scale, and shape).
}

\usage{
disttree(formula, data, subset, na.action = na.pass, weights, offset,
           cluster, family = NO(), control = disttree_control(...), 
           converged = NULL, scores = NULL, doFit = TRUE, ...)
}

\arguments{
  \item{formula}{a symbolic description of the model to be fit. This
    should be of type \code{y ~ x1 + x2}
    where \code{y} should be the response variable
    and \code{x1} and \code{x2} are used as partitioning variables.}
  \item{data}{an optional data frame containing the variables in the model.}
  \item{subset}{an optional vector specifying a subset of observations to be
                 used in the fitting process.}
  \item{na.action}{a function which indicates what should happen when the data
                   contain missing value.}
  \item{weights}{optional numeric vector of case weights.}
  \item{offset}{an optional vector of offset values.}
  \item{cluster}{an optional factor indicating independent clusters.
                  Highly experimental, use at your own risk.}
  \item{family}{specification of the response distribution.
    Either a \code{\link[gamlss.dist]{gamlss.family}} object, a list generating function or a family list.}
  \item{control}{control arguments passed to \code{\link[partykit]{extree_fit}} via \code{\link[disttree]{disttree_control}}.}
  %\item{type.hessian}{Can either be 'checklist', 'analytic' or 'numeric' to decide how the hessian matrix should be calculated in the fitting process in \code{distfit}. For 'checklist' it is checked whether a function 'hdist' is given in the family list. If so, 'type.hessian' is set to 'analytic', otherwise to 'numeric'.}
  \item{converged}{an optional function for checking user-defined criteria before splits are implemented.}
  \item{scores}{an optional named list of scores to be attached to ordered factors.}
  \item{doFit}{a logical indicating if the tree shall be grown (TRUE) or not (FALSE).}
  %\item{bd}{binomial denominator: additional parameter needed for binomial gamlss.families}
  %\item{decorrelate}{Specification of the type of decorrelation for the
  %  empirical estimating functions (or scores) either \code{"none"} or
  %  \code{"opg"} (for the outer product of gradients) or \code{"vcov"}
  %  (for the variance-covariance matrix, assuming this is an estimate
  %  of the Fisher information).}
  %\item{censtype}{Can either be 'none', 'left' or 'right' to set the type of censoring for censored response.}
  %\item{censpoint}{numeric value. Censoring point can be set for censored response.}
  %\item{ocontrol}{List with control parameters passed to
  %  \code{\link[stats]{optim}} in \code{\link[disttree]{distfit}}.}
  \item{\dots}{arguments to be used to form the default \code{control} argument
    if it is not supplied directly.}
}

\details{
  Distributional regression trees are an application of model-based recursive partitioning
  and unbiased recursive partitioning (implemented in \code{\link[partykit]{extree_fit}})
  to parametric model fits based on the GAMLSS family of distributions.
}

\value{
  An object of S3 class \code{disttree} inheriting from class \code{modelparty}.
}

\seealso{\code{\link[partykit]{mob}}, \code{\link[partykit]{ctree}}, \code{\link[partykit]{extree_fit}}, \code{\link[disttree]{distfit}}}

\examples{
tr <- disttree(dist ~ speed, data = cars)
print(tr)

plot(tr)
plot(as.constparty(tr))
}

\keyword{regression tree, parametric modeling}
