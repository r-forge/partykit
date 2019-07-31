\name{dvonmises}
\alias{dvonmises}
\title{Von Mises Density}
\usage{
dvonmises(y, mu, kappa, log = FALSE)
}
\arguments{
\item{y}{vector of observations.}

\item{mu}{vector of location parameters.}

\item{kappa}{vector of concentration parameters.}

\item{log}{logical; if TRUE, probabilities p are given as log(p)}
}
\value{
Von Mises Density
}
\description{
Density function for the von Mises distribution with location parameter mu 
and concentration parameter kappa.
}
