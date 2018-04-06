\name{stationsTyrol}
\alias{stationsTyrol}

\title{Observation stations}

\description{All 95 observations stations including all necessary information about each station.}

\usage{data("stationsTyrol")}

\format{
A data frame containing 95 observation stations and 5 variables.
\describe{
  \item{name}{character. Stationname.}
  \item{id}{numeric. Stationnumber}
  \item{lon}{numeric. Longitutde.}
  \item{lat}{numeric. Latitude.}
  \item{alt}{numeric. Altitude.}
  }
}

\references{
{B}undesministerium fuer {L}and und {F}orstwirtschaft, {U}mwelt und {W}asserwirtschaft ({BMLFUW}), {A}bteilung {IV/4} - {W}asserhaushalt (2016).
{Available at \url{http://ehyd.gv.at}},
{Accessed: {F}ebruary 29 2016}
}


\examples{
data("stationsTyrol")
head(stationsTyrol)
}

\keyword{datasets}
