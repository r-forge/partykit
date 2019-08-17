\name{disttree_control}
\alias{disttree_control}

\title{Auxiliary Function for Controlling disttree Fitting}

\description{
  Auxiliary function for \code{disttree} fitting. Specifies a list of control values for fitting
  a distributional regression tree or forest.
}
\usage{
disttree_control(type.tree = NULL, type.hessian = c("checklist",
                 "analytic", "numeric"), decorrelate = c("none", "opg",
                 "vcov"), method = "L-BFGS-B", optim.control = list(),
                 lower = -Inf, upper = Inf, minsplit = NULL, minbucket =
                 NULL, splittry = 1L, splitflavour = c("ctree",
                 "exhaustive"), testflavour = c("ctree", "mfluc",
                 "guide"), terminal = "object", model = TRUE, inner = "object",
                 restart = TRUE, breakties = FALSE, parm = NULL, dfsplit = TRUE,
                 vcov = c("opg", "info", "sandwich"), ordinal = c("chisq", "max", "L2"),
                 ytype = c("vector", "data.frame", "matrix"), trim = 0.1, 
                 guide_interaction = FALSE, interaction = FALSE, guide_parm = NULL, 
                 guide_testtype = c("max", "sum", "coin"), guide_decorrelate = "vcov", 
                 xgroups = NULL, ygroups = NULL, weighted.scores = FALSE, ...)
}
\arguments{
  \item{type.tree}{\code{NULL} or character specifying which type of tree should be fitted: 
    Either based on model-based recursive partitioning \code{type.tree="mob"} or 
    unbiased recursive partitioning \code{type.tree="ctree"}.}
  \item{type.hessian}{Add description.}
  \item{decorrelate}{Add description.}
  \item{method}{optimization method passed to \code{\link{optim}}.}
  \item{optim.control}{a list with further arguments to be passed to 'fn and 'gr' 
    in  \code{\link{optim}.}}
  \item{lower, upper}{bounds on the variables for the \code{"L-BFGS-B"} method,
    or bounds in which to search for method \code{"Brent"} passed to \code{\link{optim}}.}
  \item{minsplit, minbucket}{integer. The minimum number of observations in a node.
    If \code{NULL}, the default is to use 10 times the number of parameters
    to be estimated (divided by the number of responses per observation
    if that is greater than 1).}
  \item{splittry}{number of variables that are inspected for admissible splits
                  if the best split doesn`t meet the sample size constraints. 
                  FIXME: (ML) set to 1L, mob default.}
  \item{splitflavour}{use exhaustive search (\code{mob}) over splits instead of maximally
                      selected statistics (\code{ctree}). This feature may change.}
  \item{testflavour}{Add description}
  \item{terminal}{Add description}
  \item{model}{Add description}
  \item{inner}{Add description}
  \item{restart}{Add description}
  \item{breakties}{Add description}
  \item{parm}{Add description}
  \item{dfsplit}{Add description}
  \item{vcov}{Add description}
  \item{ordinal}{Add description}
  \item{ytype}{Add description}
  \item{trim}{Add description}
  \item{guide_interaction}{Add description}
  \item{interaction}{Add description}
  \item{guide_parm}{Add description}
  \item{guide_testtype}{Add description}
  \item{guide_decorrelate}{Add description}
  \item{xgroups}{Add description}
  \item{ygroups}{Add description}
  \item{weighted.scores}{Add description}
  \item{...}{Add description}
}

\value{
  A list with components named as the arguments.
}

\seealso{\code{\link{disttree}}, \code{\link[partykit]{extree_fit}}}
