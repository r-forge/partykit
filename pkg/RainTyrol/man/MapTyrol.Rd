\name{MapTyrol}
\alias{MapTyrol}
\encoding{UTF-8}

\title{Topographic data for Tyrol}

\description{Topographic data to plot a map of Tyrol and surrounding areas.}

\source{\url{https://www.data.gv.at/katalog/dataset/vgd-stichtagsdaten-1-250-000}, \url{http://www.earthenv.org/DEM}}

\usage{data("MapTyrol")}

\format{
A \code{list} of two objects: a \code{RasterLayer} containing topographic data of Tyrol and surrounding areas and a \code{SpatialPolygons} representing the border of Tyrol.
}

\references{
Robinson N, Regetz J, Guralnick R P (2014).
{{E}arth{E}nv-{DEM90}: A Nearly-Global, Void-Free, Multi-Scale Smoothed, 90m Digital Elevation Model From Fused {ASTER} and {SRTM} Data},
\emph{ISPRS Journal of Photogrammetry and Remote Sensing}, \bold{87}, 57--67.
    \doi{10.1016/j.isprsjprs.2013.11.002}

EarthEnv-DEM90e website: \url{http://www.earthenv.org/DEM.html}

{B}undesamt für {E}ich- und {V}ermessungswesen

\url{https://www.data.gv.at/katalog/dataset/vgd-stichtagsdaten-1-250-000}
}

\examples{
data("MapTyrol", package = "RainTyrol")
}

\keyword{datasets}
