\name{a_m_data}
\alias{a_m_data}

\title{Artificial Motivating Dataset}

\description{
Artificial dataset to illustrate fitting (G)LMM trees.
}

\usage{
data("exampledata")
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
data("a_m_data", package = "glmertree")
summary(a_m_data)
}

\keyword{datasets}