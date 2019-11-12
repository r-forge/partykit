\name{extree_control}
\alias{extree_control}
\title{ Control for extensible trees }
\description{

  Various parameters that control aspects of the `extree' fit.

}
\usage{
extree_control(
    criterion, 
    logmincriterion, 
    minsplit = 20L,
    minbucket = 7L, 
    minprob = 0.01, 
    nmax = Inf,
    stump = FALSE,
    lookahead = FALSE, 
    maxsurrogate = 0L, 
    numsurrogate = FALSE,
    mtry = Inf,
    maxdepth = Inf, 
    multiway = FALSE, 
    splittry = 2L,
    majority = FALSE, 
    caseweights = TRUE, 
    applyfun = NULL, 
    cores = NULL,
    saveinfo = TRUE,
    bonferroni = FALSE,
    update = NULL,
    selectfun, 
    splitfun, 
    svselectfun,
    svsplitfun 
)
}
\arguments{
  \item{criterion}{ a character specifying the criterion used for variable and split point selection. Either "p.value" or "statistic". }
  \item{logmincriterion}{ the value of the test statistic or 1 - p-value that
                       must be exceeded in order to implement a split on
                       the log-scale. }
  \item{minsplit}{ the minimum sum of weights in a node in order to be considered
                   for splitting. }
  \item{minbucket}{ the minimum sum of weights in a terminal node. }
  \item{minprob}{ proportion of observations needed to establish a terminal node.}
  \item{nmax}{}
  \item{stump}{ a logical determining whether a stump (a tree with a maximum of three
                nodes only) is to be computed. }
  \item{lookahead}{ a logical determining whether a split is implemented only
                    after checking if \code{trafo} in both daughter nodes converges.}
  \item{maxsurrogate}{ number of surrogate splits to evaluate.}
  \item{numsurrogate}{ a logical for backward-compatibility with party. If
                       \code{TRUE}, only at least ordered variables are considered for surrogate splits.}
  \item{mtry}{ number of input variables randomly sampled as candidates
               at each node for random forest like algorithms. The default
               \code{mtry = Inf} means that no random selection takes place.}
  \item{maxdepth}{ maximum depth of the tree. The default \code{maxdepth = Inf}
                   means that no restrictions are applied to tree sizes.}
  \item{multiway}{ CURRENTLY IGNORED. A logical indicating if multiway splits for all factor levels
                   are implemented for unordered factors.}
  \item{splittry}{ number of variables that are inspected for admissible splits
                   if the best split doesn't meet the sample size constraints.}
  \item{majority}{ if \code{FALSE}, observations which can't be classified to a
                   daughter node because of missing information are randomly
                   assigned (following the node distribution). If \code{TRUE},
                   they go with the majority.}
  \item{caseweights}{ a logical interpreting \code{weights} as case weights.}
  \item{applyfun}{an optional \code{\link[base]{lapply}}-style function with arguments
                  \code{function(X, FUN, \dots)}. It is used for computing the variable selection criterion.
                  The default is to use the basic \code{lapply}
                  function unless the \code{cores} argument is specified (see below).}
  \item{cores}{numeric. If set to an integer the \code{applyfun} is set to
               \code{\link[parallel]{mclapply}} with the desired number of \code{cores}.}
  \item{saveinfo}{logical. Store information about variable selection
                  procedure in \code{info} slot of each \code{partynode}.}
  \item{update}{logical. If \code{TRUE}, the data transformation is updated
                in every node.}
  \item{bonferroni}{}
  \item{selectfun}{}
  \item{splitfun}{}
  \item{svselectfun}{}
  \item{svsplitfun}{}
}
\details{
	TBD
}
\value{
  A list.
}
\keyword{misc}
