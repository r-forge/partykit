\name{circmax_simulate}
\alias{circmax_simulate}
\title{Simulated Data Set for \code{\link[circtree]{circmax}}}
\usage{
  circmax_simulate(n = 1000, beta = c(3, 5, 2), gamma = c(3, 3), seed = 111)
}
\arguments{
\item{n}{The number of Observations.}

\item{beta}{The coefficients for the intercept and the covariates of the location part.}

\item{gamma}{The coefficients for the intercept and the covariates of the concentration part.}

\item{seed}{Sets the 'seed' to a numeric value.}
}
\value{
  Data frame with simualated covariates and respective response. 
}
\description{
  This function creates artifical data set for testing the regression models for a circular 
    response by maximum likelihood estimation.
}
