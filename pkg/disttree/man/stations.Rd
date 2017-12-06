\name{stations}
\alias{stations}

\title{Observation stations}

\description{All 95 observations stations including all necessary information about each station.}

\usage{data("stations")}

\format{
A data frame containing 95 observations stations and 5 variables.
\describe{
  \item{stationnr}{numeric. Stationnumber}
  \item{stationname}{character. Stationname.}
  \item{lon}{numeric. Longitutde.}
  \item{lat}{numeric. Latitude.}
  \item{height}{numeric. Height.}
  }
}

\references{
{B}undesministerium fuer {L}and und {F}orstwirtschaft, {U}mwelt und {W}asserwirtschaft ({BMLFUW}), {A}bteilung {IV/4} - {W}asserhaushalt (2016).
{Available at \url{http://ehyd.gv.at}},
{Accessed: {F}ebruary 29 2016}
}


\examples{
data("stations")
head(stations)
}

\keyword{datasets}
