\name{rainTyrol}
\alias{rainTyrol}

\title{Observations and covariates for all 95 stations}

\description{Observations of precipitation sums and numerical ensemble predictions of meteorological quantities.}

\usage{data("rainTyrol")}

\format{
A list of data frames, one data frame for each of the 95 observation stations. Each data frame consists of power transformed observations and numerical ensemble predictions of meteorological quantities. The base variables of the numerical ensemble predictions are listed below. For each of them variations such as ensemble mean/standard deviation/minimum/maximum are provided.
\describe{
  \item{robs}{numeric. Observed total precipitation (power transformed).}
  \item{year}{integer. Year in which the observation was taken.}
  \item{tppow}{numeric. Predicted total precipitation (power transformed).}
  \item{capepow}{numeric. Predicted convective available potential energy (power transformed).}
  \item{dswrf}{numeric. Predicted downwards short wave radiation flux (``sunshine'').}
  \item{msl}{numeric. Predicted mean sea level pressure.}
  \item{pwat}{numeric. Predicted preciptable water.}
  \item{tmax}{numeric. Predicted 2m maximum temperature.}
  \item{tcolc}{numeric. Predicted total column-integrated condensate.}
  \item{t500}{numeric. Predicted temperature on 500 hPa.}
  \item{t700}{numeric. Predicted temperature on 700 hPa.}
  \item{t8500}{numeric. Predicted temperature on 850 hPa.}
  \item{tdiff500700}{numeric. Predicted temperature difference 500 to 850 hPa.}
  \item{tdiff500850}{numeric. Predicted temperature difference 500 to 700 hPa.}
  \item{tdiff700850}{numeric. Predicted temperature difference 700 to 850 hPa.}
}
}

\references{
   Hamill T M,  Bates G T, Whitaker J S, Murray D R, Fiorino M, Galarneau Jr. T J, Zhu Y, Lapenta W (2013).
   {NOAA}'s Second-Generation Global Medium-Range Ensemble Reforecast Dataset.
   \emph{Bulletin of the American Meteorological Society}, \bold{94}(10), 1553--1565.
    \doi{10.1175/BAMS-D-12-00014.1}
  

  
{BMLFUW} (2016).
{B}undesministerium f\"ur {L}and und {F}orstwirtschaft, {U}mwelt und {W}asserwirtschaft ({BMLFUW}), {A}bteilung {IV/4} -- {W}asserhaushalt.
Available at \url{http://ehyd.gv.at}.
Accessed: 2016--02--29.
}


\examples{
data("rainTyrol")
head(rainTyrol$Axams)
}

\keyword{datasets}
