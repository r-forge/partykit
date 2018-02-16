\name{rainAxams}
\alias{rainAxams}

\title{Observations and covariates for station Axams}

\description{Observations of precipitation sums and numerical ensemble predictions of meteorological quantities.}

\usage{data("rainAxams")}

\format{
A data frames consisting of power transformed observations of precipitation sums and numerical ensemble predictions of meteorological quantities for the station Axams. The base variables of the numerical ensemble predictions are listed below. For each of them variations such as ensemble mean/standard deviation/minimum/maximum are provided.
\describe{
  \item{robs}{numeric. Observed total precipitation (power transformed).}
  \item{year}{integer. Year in which the observation was taken.}
  \item{tppow_max, tppow_mean, tppow_mean0612, tppow_mean1218,
        tppow_mean1824, tppow_mean2430, tppow_min, tppow_sprd,
        tppow_sprd0612, tppow_sprd1218, tppow_sprd1824, tppow_sprd2430}{numeric. Predicted total precipitation (power transformed).}
  \item{capepow_max, capepow_mean, capepow_mean0612,
        capepow_mean1218, capepow_mean1224, capepow_mean1230, capepow_min,
        capepow_sprd, capepow_sprd0612, capepow_sprd1218, capepow_sprd1224,
        capepow_sprd1230}{numeric. Predicted convective available potential energy (power transformed).}
  \item{dswrf_mean_max, dswrf_mean_mean, dswrf_mean_min,
        dswrf_sprd_max, dswrf_sprd_mean, dswrf_sprd_min}{numeric. Predicted downwards short wave radiation flux (``sunshine'').}
  \item{msl_diff, msl_mean_max, msl_mean_mean, msl_mean_min,
        msl_sprd_max, msl_sprd_mean, msl_sprd_min}{numeric. Predicted mean sea level pressure.}
  \item{pwat_mean_max, pwat_mean_mean, pwat_mean_min,
        pwat_sprd_max, pwat_sprd_mean, pwat_sprd_min}{numeric. Predicted preciptable water.}
  \item{tmax_mean_max, tmax_mean_mean, tmax_mean_min, 
        tmax_sprd_max, tmax_sprd_mean, tmax_sprd_min}{numeric. Predicted 2m maximum temperature.}
  \item{tcolc_mean_max, tcolc_mean_mean, tcolc_mean_min, 
        tcolc_sprd_max, tcolc_sprd_mean, tcolc_sprd_min}{numeric. Predicted total column-integrated condensate.}
  \item{t500_mean_max, t500_mean_mean, t500_mean_min,
        t500_sprd_max, t500_sprd_mean, t500_sprd_min}{numeric. Predicted temperature on 500 hPa.}
  \item{t700_mean_max, t700_mean_mean, t700_mean_min,
        t700_sprd_max, t700_sprd_mean, t700_sprd_min}{numeric. Predicted temperature on 700 hPa.}
  \item{t850_mean_max, t850_mean_mean, t850_mean_min,
        t850_sprd_max, t850_sprd_mean, t850_sprd_min}{numeric. Predicted temperature on 850 hPa.}
  \item{tdiff500700_max, tdiff500700_mean, tdiff500700_min}{numeric. Predicted temperature difference 500 to 700 hPa.}
  \item{tdiff500850_max, tdiff500850_mean, tdiff500850_min}{numeric. Predicted temperature difference 500 to 850 hPa.}
  \item{tdiff700850_max, tdiff700850_mean, tdiff700850_min}{numeric. Predicted temperature difference 700 to 850 hPa.}
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
data("rainAxams")
head(rainAxams)
colnames(rainAxams)
}

\keyword{datasets}
