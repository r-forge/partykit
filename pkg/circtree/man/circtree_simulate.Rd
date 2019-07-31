\name{circtree_simulate}
\alias{circtree_simulate}
\title{Simulated Data Set for \code{\link{circtree}}}
\usage{
  circtree_simulate(n = 1000, mu = c(0, 2, 1), kappa = c(3, 3, 1), seed = 111)
}
\arguments{
\item{n}{The number of Observations.}

\item{mu}{The distribution parameters for the location part. \strong{Currently exactly three parameters necessary.}}

\item{kappa}{The distribution parameters for the concentration part. \strong{Currently exactly three parameters necessary.}}

\item{seed}{Sets the 'seed' to a numeric value.}
}
\value{
  Data frame with simualated covariates and respective response. 
}
\description{
  This function creates artifical data set for testing the regression trees employing a 
    von Mises distribution.
}
