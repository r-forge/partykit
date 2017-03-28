\name{mob2_control}
\alias{mob2_control}

\title{Control Parameters for Model-Based Partitioning}

\description{
  Various parameters that control aspects the fitting algorithm
  for recursively partitioned \code{\link{mob2}} models.
}

\usage{
mob2_control(testflavour = "mfluc", testtype = "Bonferroni", alpha = 0.05,
 parm = NULL, bonferroni = TRUE, breakties = FALSE, nrep = 10000L,
 ordinal = "chisq", trim = 0.1, nresample = 9999L, teststat = "quadratic",
 splittry = 2, vcov = "opg", splitflavour = "exhaustive",
 splittest = FALSE, splitstat = "quadratic", catsplit = "binary",
 numsplit = "left", majority = FALSE, maxsurrogate = 0, MIA = FALSE,
 numsurrogate = FALSE, lookahead = FALSE, maxdepth = Inf,
 minbucket = minsize, minprob = 0.01, minsplit = minsize, minsize = 20L,
 stump = FALSE, caseweights = TRUE, dfsplit = TRUE, applyfun = NULL,
 cores = NULL, restart = FALSE, inner = "object", model = TRUE,
 terminal = "object", mtry = Inf)
}

\arguments{

Variable selection:

\item{testflavour}{split variable selection strategy. The default,  "mfluc",
uses M-fluctuation tests; "ctree" uses permutation tests as in ctree and
"exhaustive" goes through all variables and splits and selects the one
corresponding to the greatest sum of objective functions in the leaf node
models.}
 
\item{alpha}{a double, the significance level for variable selection. If
testflavour = "exhaustive" this is ignored.}
 
\item{parm}{numeric or character. Number or name of model parameters included
in the parameter instability tests (by default all parameters are included).}
 
\item{splittry}{number of variables that are inspected for admissible splits if
the best split doesn't meet the sample size constraints}
 
\item{testtype}{Only relevant if testflavour =  "ctree". A character specifying
how to compute the distribution of the test statistic ("Bonferroni",
"MonteCarlo", "Univariate", "Teststatistic")}
 
\item{mincriterion}{the value of the test statistic or 1 - p-value that must be
exceeded in order to implement a split.}
 
\item{nresample}{number of permutations for testtype = "MonteCarlo".}
 
\item{teststat}{Only relevant if testflavour =  "ctree". A character specfying
the type of the test statistic ("quadratic" or "maximum") to be applied for
variable selection}
 
\item{bonferroni}{logical. Should p vlaues be Bonferroni corrected?}
 
\item{breakties}{logical. Only relevant iif testflavour =  "mfluc".  Should
ties in numeric variables be broken randomly for computing the associated
parameter instability test?}
 
\item{nrep}{numeric. Only relevant if testflavour =  "mfluc". Number of
replications in the simulation of p-values for the ordinal "L2" statistic (if
used).}
 
\item{ordinal}{character indicating which type of parameter instability test
should be employed for ordinal partitioning variables (i.e., ordered factors).
Only relevant if testflavour =  "mfluc". This can be "chisq", "max", or "L2".
If "chisq" then the variable is treated as unordered and a chi-squared test is
performed. If "L2", then a maxLM-type test as for numeric variables is carried
out but correcting for ties. This requires simulation of p-values via catL2BB
and requires some computation time. For "max" a weighted double maximum test is
used that computes p-values via pmvnorm.}
 
\item{trim}{numeric. Only relevant if testflavour =  "mfluc". This specifies
the trimming in the parameter instability test for the numerical variables. If
smaller than 1, it is interpreted as the fraction relative to the current node
size.}
 
\item{vcov}{character indicating which type of covariance matrix estimator
should be employed in the parameter instability tests. Only relevant if
testflavour =  "mfluc". The default is the outer product of gradients ("opg").
Alternatively, vcov = "info" employs the information matrix and vcov =
"sandwich" the sandwich matrix (both of which are only sensible for maximum
likelihood estimation).}



Split point selection:
 
\item{splitflavour}{split point selection strategy. The default, "exhaustive",
goes through all possible split points and chooses the one corresponding to the
greatest sum of objective functions in the leaf node models; "ctree" usese
permutation tests}
 
\item{catsplit}{character indicating how (unordered) categorical variables
should be splitted. By default the best "binary"split is searched (by
minimizing the objective function). Alternatively, if set to "multiway", the
node is simply splitted into all levels of the categorical variable (equal to
multiway = TRUE in ctree_control). }
 
\item{numsplit}{character indicating how splits for numeric variables should be
justified. Because any splitpoint in the interval between the last observation
from the left child segment and the first observation from the right child
segment leads to the same observed split, two options are available in
mob_control: Either, the split is"left"-justified (the default for backward
compatibility; equal to intersplit = FALSE in ctree_control) or
"center"-justified using the midpoint of the possible interval (equal to
intersplit = TRUE in ctree_control). }
 
\item{splittest}{Only relevant if splitflavour =  "ctree". A logical changing
linear (the default FALSE) to maximally selected statistics for variable
selection. Currently needs testtype = "MonteCarlo".}
 
\item{splitstat}{Only relevant if splitflavour =  "ctree". A character
specfying the type of the test statistic ("quadratic" or "maximum") to be
applied for split point selection}





Splitting in case of missings:
 
\item{majority}{if FALSE, observations wich can't be classified to a daughter
node because of missing information are randomly assigned (following the node
distribution). If TRUE (default), they go with the majority.}
 
\item{maxsurrogate}{number of surrogate splits to evaluate. Note that currently
only surrogate splits in ordered variables are implemented}
 
\item{MIA}{a logical determining the treatment of NA as a category in split,
see Twala et al. (2008)}
 
\item{numsurrogate}{a logical for backward-compatibility with party. If TRUE,
only numeric variables are considered for surrogate splits}



Stopping:
 
\item{lookahead}{a logical determining whether a split is implemented only
after checking if tests in both daughter nodes can be performed.}
 
\item{maxdepth}{maximum depth of the tree. The default maxdepth = Inf means
that no restrictions are applied to tree sizes}
 
\item{minbucket}{the minimum sum of weights in a terminal node}
 
\item{minprob}{proportion of observations needed to establish a terminal node}
 
\item{minsplit}{the minimum sum of weights in a node in order to be considered
for splitting}
 
\item{minsize}{minsplit}
 
\item{stump}{a logical determining whether a stump (a tree with a maximum three
nodes) is to be computed}



Observation number and degrees of freedom:
 
\item{caseweights}{a logical interpreting weights as case weights}
 
\item{dfsplit}{logical or numeric. as.integer(dfsplit) is the degrees of
freedom per selected split employed when computing information criteria etc.}



Computation:
 
\item{applyfun}{an optional lapply-style function with arguments function(X,
FUN, …). It is used for computing the variable selection criterion. The default
is to use the basic lapply function unless the cores argument is specified (see
below).}
 
\item{cores}{numeric. If set to an integer the applyfun is set to mclapply with
the desired number of cores.}
 
\item{restart}{logical. When determining the optimal split point in a numerical
variable: Should model estimation be restarted with NULL starting values for
each split? If FALSE (default), then the parameter estimates from the previous
split point are used as starting values for the next split point (because in
practice the difference are often not huge).}



Return options:
 
\item{inner}{character. Specification of whether the model object should be
stored in each node. If NULL, no additional information is stored.}
 
\item{model}{logical. Should the full model frame be stored in the resulting
object?}
 
\item{terminal}{character. Specification of which additional information
("estfun", "object", or both) should be stored in each node. If NULL, no
additional information is stored.}




For extensions:
 
\item{mtry}{number of input variables randomly sampled as candidates at each
node for random forest like algorithms. The default mtry = Inf means that no
random selection takes place}

}



\details{
The most important decisions to be made when computing a model-based tree is on
how to do variable and split point selection.

  \subsection{Variable selection (\code{testflavour})}{
  \itemize{
    \item \bold{M-Fluctuation tests (\code{"mfluc"}):} This test type is the default
    for model-based recursive partitioning. It searches for break points in the 
    score contributions (see \code{\link[sandwich]{estfun}}). By default p-values are Bonferroni corrected. To turn this 
    off, set \code{bonferroni = FALSE}.
    
    \item \bold{Conditional inference tests (\code{"ctree"}):} This test type has 
    in the past been used only for conditional inference trees (\code{\link{ctree}}). In
    \code{\link{mob}} it is used to find variables that correlate with the model
    score contributions (\code{\link[sandwich]{estfun}}). There are several options on
    how to do conditional inference tests:
    \itemize{
    \item \code{testtype = "Bonferroni"} Bonferroni adjusted
    p-values based on the asymptotic distribution are used (default).
    \item \code{testtype = "Univariate"} Non-adjusted
    p-values based on the asymptotic distribution are used.
    \item \code{testtype = "MonteCarlo"} Non-adjusted
    p-values permutation based are used.
    \item \code{testtype = c("MonteCarlo", "Bonferroni")} Bonferroni adjusted
    p-values permutation based are used.
    \item \code{testtype = "Teststatistic"} Instead of the p-value the test 
    statistic is used. In this case \code{mincriterion} is the value of the 
    test statistic that must be exceeded.
    }
    Using the p-value to select the variable means choosing the variable corresponding to most significant
    deviation from independence between scores and partitioning variable. Using the test
    statistic means choosing the variable with the strongest correlation to the scores.
    
    \item \bold{Exhaustive search (\code{"exhaustive"}):} In exhaustive search all 
    all possible two-way splits of all variables are implemented and the objective
    function (\code{\link{logLik}}) is computed. The variable is selected that leads
    to the greatest sum in objective functions of two child nodes.
  }
  }
  
  \subsection{Split point selection (\code{splitflavour})}{
  \itemize{
    \item \bold{Exhaustive search (\code{"exhaustive"}):} Exhaustive search goes through
    all possible split points of the selected variable and chooses the one leading to
    the maximum sum of objective function (\code{\link{logLik}}) in the child nodes.
    
    \item \bold{Conditional inference tests (\code{"ctree"}):} This version of split
    point selection goes through all possible split points of the selected variable,
    computes the conditional inference test with Null hypothesis of no difference in
    score contributions in the two child nodes (or in other words: independence between
    child nodes and scores) and chooses the split point corresponding to the greatest 
    test statistic.
  }
  }

  See \code{\link{mob}} for more details and references.
}

\seealso{\code{\link{mob}}}

\value{
  A list of class \code{mob_control} containing the control parameters.
}

\keyword{misc}
