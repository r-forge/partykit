\name{demo_obs+cov}
\alias{demo_obs+cov}

\title{Observations and covariates for all 95 stations}

\description{Observations of precipitation sums and numerical ensemble predictions of meteorological quantities.}

\usage{data("demo_obs+cov")}

\format{
A list of data frames, one data frame for each of the 95 observation stations. Each data frame consists of power transformed observations and numberical ensemble predictions of meteorological quantities.
\describe{
  \item{robs}{numeric. Power transformed observations.}
}
}

\references{
   Hamill T M,  Bates G T, Whitaker J S, Murray D R, Fiorino M, Galarneau Jr. T J, Zhu Y, Lapenta W (2013).
   {NOAA}'s Second-Generation Global Medium-Range Ensemble Reforecast Dataset.
   \emph{Bulletin of the American Meteorological Society}, \bold{94}(10), 1553--1565.
    \doi{10.1175/BAMS-D-12-00014.1}
}

\examples{
data("demo_obs+cov")
head(datalist$Axams)
}

\keyword{datasets}
