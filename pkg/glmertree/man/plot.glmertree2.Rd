\name{plot.glmertree2}
\alias{plot.lmertree2}
\alias{plot.glmertree2}

\encoding{latin1}

\title{Plotting (Generalized) Linear Mixed Model Trees}

\description{
  \code{plot} method for \code{(g)lmertree} objects.
}

\usage{
%\method{plot}{lmertree2}
plot.lmertree2(x, which = "all", ask = TRUE, type = "extended", 
    observed = TRUE, fitmean = "marginal", tp_args = list(), 
    drop_terminal = TRUE, terminal_panel = NULL, \dots)
%\method{plot}{glmertree2}
plot.glmertree2(x, which = "all", ask = TRUE, type = "extended", 
    observed = TRUE, fitmean = "marginal", tp_args = list(), 
    drop_terminal = TRUE, terminal_panel = NULL, \dots)
}

\arguments{
  \item{x}{an object of class \code{lmertree} or \code{glmertree}.}
  \item{which}{character; \code{"all"} (default), \code{"tree"} or 
      \code{"random"}. Specifies whether, tree, random effects or both should
      be plotted.}
  \item{ask}{logical. Should user be asked for input, before a new figure
      is drawn?}
  \item{type}{character; \code{"extended"} (default) or \code{"simple"}.
      \code{"extended"} yields a plotted tree with observed data and/or
      fitted means in the terminal nodes; \code{"simple"} yields a plotted
      tree with fixed and/or random effects coefficients presented in the
      terminal nodes.}
  \item{observed}{logical. Should observed datapoint be plotted?}
  \item{fitmean}{character. \code{"marginal"} (default), \code{"combined"} 
      or \code{"none"}. Specifies whether and how fitted means should be 
      visualized. Only used when predictor variables for the node-specific 
      (G)LMs are specified. If \code{"marginal"}, fitted values 
      will be calculated, with all remaining (random and fixed effects) fixed 
      at their (population and sample) means. If \code{"combined"}, fitted 
      values will calculated, based on observed values of the remaining 
      (random and fixed-effects) predictor variables.} 
  \item{tp_args}{list of arguments to be passed to panel generating function
      \code{node_glmertree}. See arguments \code{node_bivplot} in 
      \code{\link[partykit]{panelfunctions}}.}
  \item{drop_terminal}{logical. Should all terminal nodes be plotted at the 
      bottom?}
  \item{terminal_panel}{an optional panel function to be passed to 
      \code{plot.party()}. See \code{\link[partykit]{party-plot}} 
      documentation for details.}
  \item{\dots}{Additional arguments to be passed to \code{plot.party()}.
      See \code{\link[partykit]{party-plot}} documentation for details.}
}

  
\details{
  The code is still under development and might change in future versions.
}

\value{
The function returns a list with the following objects:
  \item{tree}{The final \code{lmtree}/\code{glmtree}.}
  \item{lmer}{The final \code{lmer} random-effects model.}
  \item{ranef}{The corresponding random effects of \code{lmer}.} 
  \item{varcorr}{The corresponding \code{VarCorr(lmer)}.}
  \item{variance}{The corresponding \code{attr(VarCorr(lmer), "sc")^2}.}
  \item{data}{The dataset specified with the \code{data} argument
    including added auxiliary variables \code{.ranef} and \code{.tree}
    from the last iteration.}
  \item{loglik}{The log-likelihood value of the last iteration.}
  \item{iterations}{The number of iterations used to estimate the \code{lmertree}.} 
  \item{maxit}{The maximum number of iterations specified with the \code{maxit} argument.}
  \item{ranefstart}{The random effects used as an offset, as specified with
    the \code{ranefstart} argument.}
  \item{formula}{The formula as specified with the \code{formula} argument.}
  \item{randomformula}{The formula as specified with the \code{randomformula} argument.}
  \item{abstol}{The prespecified value for the change in log-likelihood to evaluate
    convergence, as specified with the \code{abstol} argument.} 
  \item{mob.control}{A list containing control parameters passed to
    \code{lmtree()}, as specified with \dots.}
  \item{lmer.control}{A list containing control parameters passed to
    \code{lmer()}, as specified in the \code{lmer.control} argument.}
  \item{joint}{Whether the fixed effects from the tree were (re-)estimated jointly along 
    with the random effects, specified with the \code{joint} argument.}  
}

\references{
Fokkema M, Smits N, Zeileis A, Hothorn T, Kelderman H (2018).
  \dQuote{Detecting Treatment-Subgroup Interactions in Clustered Data
  with Generalized Linear Mixed-Effects Model Trees}. Behavior Research 
  Methods, 50(5), 2016-2034. \url{https://doi.org/10.3758/s13428-017-0971-x}
}


\examples{
\donttest{
## load artificial example data
data("DepressionDemo", package = "glmertree")

## fit normal linear regression LMM tree for continuous outcome
lt <- lmertree(depression ~ treatment + age | cluster | anxiety + duration,
  data = DepressionDemo)
plot.lmertree2(lt)
plot.lmertree2(lt, type = "simple")
plot.lmertree2(lt, which = "tree", fitmean = "combined")
plot.lmertree2(lt, which = "tree", fitmean = "none")
plot.lmertree2(lt, which = "tree", observed = FALSE)

## fit logistic regression GLMM tree for binary outcome
gt <- glmertree(depression_bin ~ treatment | cluster | 
  age + anxiety + duration, data = DepressionDemo)
plot.glmertree2(gt)
plot.glmertree2(gt, type = "simple")
plot.glmertree2(gt, which = "tree", fitmean = "combined")
plot.glmertree2(gt, which = "tree", fitmean = "none")
plot.glmertree2(gt, which = "tree", observed = FALSE)
}
}

\seealso{\code{\link{lmertree}}, \code{\link{glmertree}}, 
    \code{\link[partykit]{party-plot}}.}

\keyword{hplot}