\name{rainAxams}
\alias{rainAxams}

\title{Observations and covariates for station Axams}

\description{
Observations of precipitation sums and weather forecasts of a set of
meteorological quantities from an ensemble prediction system for one
specific site.  This site is Axams located in the Eastern European Alps
(11.28E 47.23N, 890 meters a.m.s.l.).
}

\usage{data("rainAxams")}

\details{
The site is maintained by the hydrographical service Tyrol and provides daily
precipitation sums reported at 06~UTC. Before published, the observations have
been quality-controlled by the maintainer.

The forecast data is based on the second-generation global ensemble reforecast
dataset and consists of range of different meteorological quantities for day
one (forecast horizon +6 to +30 hours ahead). The forecasts have been
bi-linearly interpolated to the station location.
}

\format{
A \code{data.frame} consisting of the station's name, observation day and year,
power transformed observations of daily precipitation sums and the corresponding 
meteorological ensemble predictions for station Axams. The base variables of 
the numerical ensemble predictions are listed below. For each of them variations 
such as ensemble mean/standard deviation/minimum/maximum are included in the dataset.
All ``power transformed'' values use the same power parameter \code{p=1/1.6}.

\describe{
  \item{station}{\code{character}. Name of the observation station.}
  \item{robs}{\code{numeric}. Observed total precipitation (power transformed).}
  \item{year}{\code{integer}. Year in which the observation was taken.}
  \item{day}{\code{integer}. Day for which the observation was taken.}
  \item{tppow_mean, tppow_sprd, tppow_min, tppow_max, 
        tppow_mean0612, tppow_mean1218, tppow_mean1824, tppow_mean2430,
        ppow_sprd0612, tppow_sprd1218, tppow_sprd1824, tppow_sprd2430}{\code{numeric}.
        Predicted total precipitation (power transformed).}
  \item{capepow_mean, capepow_sprd, capepow_min, capepow_max, 
        capepow_mean0612, capepow_mean1218, capepow_mean1224, capepow_mean1230,
        capepow_sprd0612, capepow_sprd1218, capepow_sprd1224, capepow_sprd1230}{\code{numeric}.
        Predicted convective available potential energy (power transformed).}
  \item{dswrf_mean_mean, dswrf_mean_min, dswrf_mean_max, 
        dswrf_sprd_mean, dswrf_sprd_min, dswrf_sprd_max}{\code{numeric}.
        Predicted downwards shortwave radiation flux (``sunshine'').}
  \item{msl_diff, msl_mean_mean, msl_mean_min, msl_mean_max,
        msl_sprd_mean, msl_sprd_min, msl_sprd_max}{\code{numeric}.
        Predicted mean sea level pressure.}
  \item{pwat_mean_mean, pwat_mean_min, pwat_mean_max,
        pwat_sprd_mean, pwat_sprd_min, pwat_sprd_max}{\code{numeric}.
        Predicted precipitable water.}
  \item{tcolc_mean_mean, tcolc_mean_min, tcolc_mean_max,
        tcolc_sprd_mean, tcolc_sprd_min, tcolc_sprd_max}{\code{numeric}.
        Predicted total column-integrated condensate.}
  \item{tmax_mean_mean, tmax_mean_min, tmax_mean_max,
        tmax_sprd_mean, tmax_sprd_min, tmax_sprd_max}{\code{numeric}.
        Predicted 2m maximum temperature.}
  \item{t500_mean_mean, t500_mean_min, t500_mean_max,
        t500_sprd_mean, t500_sprd_min, t500_sprd_max}{\code{numeric}.
        Predicted temperature on 500 hPa.}
  \item{t700_mean_mean, t700_mean_min, t700_mean_max,
        t700_sprd_mean, t700_sprd_min, t700_sprd_max}{\code{numeric}.
        Predicted temperature on 700 hPa.}
  \item{t850_mean_mean, t850_mean_min, t850_mean_max,
        t850_sprd_mean, t850_sprd_min, t850_sprd_max}{\code{numeric}.
        Predicted temperature on 850 hPa.}
  \item{tdiff500850_mean, tdiff500850_min, tdiff500850_max}{\code{numeric}.
        Predicted temperature difference 500 hPa to 850 hPa.}
  \item{tdiff700850_mean, tdiff700850_min, tdiff700850_max}{\code{numeric}.
        Predicted temperature difference 700 hPa to 850 hPa.}
  \item{tdiff500700_mean, tdiff500700_min, tdiff500700_max}{\code{numeric}.
        Predicted temperature difference 500 hPa to 700 hPa.}
}
}

\references{
   Hamill T M,  Bates G T, Whitaker J S, Murray D R, Fiorino M, Galarneau Jr. T J, Zhu Y, Lapenta W (2013).
   {NOAA}'s Second-Generation Global Medium-Range Ensemble Reforecast Dataset.
   \emph{Bulletin of the American Meteorological Society}, \bold{94}(10), 1553--1565.
    \doi{10.1175/BAMS-D-12-00014.1}
      
   {BMLFUW} (2016).
   {B}undesministerium f\"ur {L}and und {F}orstwirtschaft, {U}mwelt und {W}asserwirtschaft ({BMLFUW}),
   {A}bteilung {IV/4} -- {W}asserhaushalt. Available at \url{http://ehyd.gv.at}. Accessed: 2016--02--29.
}

\examples{
data("rainAxams")
head(rainAxams)
colnames(rainAxams)
}

\keyword{datasets}
