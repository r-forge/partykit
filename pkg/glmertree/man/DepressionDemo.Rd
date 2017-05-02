\name{DepressionDemo}
\alias{DepressionDemo}

\title{Artificial Motivating Dataset}

\description{
Artificial dataset to illustrate fitting (G)LMM trees.
}

\usage{
data("DepressionDemo")
}

\format{
A data frame containing 150 observations on 6 variables:
\describe{
  \item{depression}{numeric. Continuous treatment outcome variable.}
  \item{treatment}{factor. Binary treatment variable.}
  \item{cluster}{factor. Indicator for cluster with 10 levels.}
  \item{age}{numeric. Continuous partitioning variable.}
  \item{anxiety}{numeric. Continuous partitioning variable.}
  \item{duration}{numeric. Continuous partitioning variable.}
}}

\details{
The data contains three treatment effect subgroups with respect to the 
continuous partitioning variables anxiety and duration. The \code{cluster}
variable has an additional random intercept that should be accounted for.
The outcome is assesse by the continuous variable depression.
}

\seealso{\code{\link{lmertree}}}

\examples{
data("DepressionDemo", package = "glmertree")
summary(DepressionDemo)
}

\keyword{datasets}