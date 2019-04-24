\name{GrowthCurveDemo}
\alias{GrowthCurveDemo}

\title{Artificial dataset for partitioning of linear growth curve models}

\description{
Artificial dataset to illustrate fitting of LMM trees with growth curve models  
the terminal nodes. 
}

\usage{
data("GrowthCurveDemo")
}

\format{
A data frame containing 1250 repeated observations on 250 persons:
\describe{
  \item{person}{numeric. Indicator linking repeated measurements to persons.}
  \item{time}{factor. Indicator for timepoint.}
  \item{y}{numeric. Response variable.}
  \item{x1}{numeric. Potential partitioning variable.}
  \item{x2}{numeric. Potential partitioning variable.}
  \item{x3}{numeric. Potential partitioning variable.}
  \item{x4}{numeric. Potential partitioning variable.}
  \item{x5}{numeric. Potential partitioning variable.}
  \item{x6}{numeric. Potential partitioning variable.}
  \item{x7}{numeric. Potential partitioning variable.}
  \item{x8}{numeric. Potential partitioning variable.}
  \item{x9}{numeric. Potential partitioning variable.}
  \item{x10}{numeric. Potential partitioning variable.}
  \item{x11}{numeric. Potential partitioning variable.}
  \item{x12}{numeric. Potential partitioning variable.}
  \item{x13}{numeric. Potential partitioning variable.}
  \item{x14}{numeric. Potential partitioning variable.}
  \item{x15}{numeric. Potential partitioning variable.}
  \item{x16}{numeric. Potential partitioning variable.}
  \item{x17}{numeric. Potential partitioning variable.}
  \item{x18}{numeric. Potential partitioning variable.}
  \item{x19}{numeric. Potential partitioning variable.}
  \item{x20}{numeric. Potential partitioning variable.}
  \item{x21}{numeric. Potential partitioning variable.}
  \item{x22}{numeric. Potential partitioning variable.}
  \item{x23}{numeric. Potential partitioning variable.}
  \item{x24}{numeric. Potential partitioning variable.}
  \item{x25}{numeric. Potential partitioning variable.}
  \item{x26}{numeric. Potential partitioning variable.}
  \item{x27}{numeric. Potential partitioning variable.}
  \item{x28}{numeric. Potential partitioning variable.}
}}

\details{Data were generated so that \code{x1}, \code{x2} and \code{x3} are
true partitioning variables, while \code{x4} through \code{x28} are noise
variables.}

\seealso{\code{\link{lmertree}}, \code{\link{glmertree}}}

\examples{
\donttest{
data("GrowthCurveDemo", package = "glmertree")
dim(GrowthCurveDemo)
names(GrowthCurveDemo)

## Fit LMM tree, while estimating random intercept for time
form <- formula(paste0("y ~ time | person | ", 
                paste0("x", 1:28, collapse = " + ")))
form
lt.default <- lmertree(form, data = GrowthCurveDemo)
plot(lt.default, which = "tree") ## yields too large tree
VarCorr(lt.default)

## Account for measurement level of the partitioning variables:
lt.cluster <- lmertree(form, cluster = person, data = GrowthCurveDemo)
plot(lt.cluster, which = "tree") ## yields correct tree
VarCorr(lt.cluster) ## yields slightly larger ranef variance

## Also estimate random slopes of time:
form.s <- formula(paste0("y ~ time | ((1+time)|person) | ", 
                  paste0("x", 1:28, collapse = " + ")))
form.s
lt.s.cluster <- lmertree(form.s, cluster = person, data = GrowthCurveDemo)
plot(lt.s.cluster, which = "tree") ## same tree as before
VarCorr(lt.s.cluster)

}
}

\keyword{datasets}