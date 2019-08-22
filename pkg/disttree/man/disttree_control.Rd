\name{disttree_control}
\alias{disttree_control}

\title{Auxiliary Function for Controlling \code{disttree} Fitting}

\description{
  Auxiliary function for \code{disttree} fitting. Specifies a list of control values for fitting
  a distributional regression tree or forest. These \code{disttree} specific control values are set 
  in addition to the control values of \code{ctree_control} and can vary from its default values.
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
  \item{type.hessian}{Can either be "checklist", "analytic" or "numeric" to decide how the hessian matrix should be calculated in the fitting process in \code{distfit}. For "checklist" it is checked whether a function "hdist" is given in the family list. If so, "type.hessian" is set to "analytic", otherwise to "numeric".}
  \item{decorrelate}{specification of the type of decorrelation for the
    empirical estimating functions (or scores) either \code{"none"} or
    \code{"opg"} (for the outer product of gradients) or \code{"vcov"}
    (for the variance-covariance matrix, assuming this is an estimate
    of the Fisher information).}
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
  \item{testflavour}{employ permutation tests (\code{ctree}) or M-fluctuation tests (\code{mfluc}).}
  \item{terminal}{character. Specification of which additional information ("estfun", "object", or both) should be stored in each terminal node. If NULL, no additional information is stored. Note that the information slot 'object' contains a slot 'estfun' as well. 
  FIXME: (LS) Should estfun always be returned within object?}
  \item{model}{logical. Should the full model frame be stored in the resulting object?}
  \item{inner}{character. Specification of which additional information ("estfun", "object", or both) should be stored in each inner node. If NULL, no additional information is stored. Note that the information slot 'object' contains a slot 'estfun' as well. 
  FIXME: (LS) Should estfun always be returned within object?}
  \item{restart}{logical. When determining the optimal split point in a numerical variable: Should model estimation be restarted with NULL starting values for each split? The default is TRUE. If FALSE, then the parameter estimates from the previous split point are used as starting values for the next split point (because in practice the difference are often not huge). (Note that in that case a for loop is used instead of the applyfun for fitting models across sample splits.)}
  \item{breakties}{logical. If M-fluctuation tests are applied, should ties in numeric variables be broken randomly for computing the associated parameter instability test?}
  \item{parm}{numeric or character. Number or name of model parameters included in the parameter instability tests if M-fluctuation tests are applied (by default all parameters are included).
  FIXME: (LS) is it really applied?}
  \item{dfsplit}{logical or numeric. as.integer(dfsplit) is the degrees of freedom per selected split employed when computing information criteria etc.
  FIXME: (LS) is it really applied?}
  \item{vcov}{character indicating which type of covariance matrix estimator should be employed in the parameter instability tests if M-fluctuation tests are applied. The default is the outer product of gradients ("opg"). Alternatively, vcov = "info" employs the information matrix and vcov = "sandwich" the sandwich matrix (both of which are only sensible for maximum likelihood estimation).}
  \item{ordinal}{character indicating which type of parameter instability test should be employed for ordinal partitioning variables (i.e., ordered factors) if M-fluctuation tests are applied. This can be "chisq", "max", or "L2". If "chisq" then the variable is treated as unordered and a chi-squared test is performed. If "L2", then a maxLM-type test as for numeric variables is carried out but correcting for ties. This requires simulation of p-values via catL2BB and requires some computation time. For "max" a weighted double maximum test is used that computes p-values via pmvnorm.}
  \item{ytype}{character. For type.tree "mob": Specification of how mob should preprocess y variable. Possible choice are: "vector", i.e., only one variable; "matrix", i.e., the model matrix of all variables; "data.frame", i.e., a data frame of all variables.,
  FIXME: (LS) handle multidim. response?}
  \item{trim}{numeric. This specifies the trimming in the parameter instability test for the numerical variables if M-fluctuation tests are applied. If smaller than 1, it is interpreted as the fraction relative to the current node size.}
  \item{guide_interaction}{logical. Should interaction tests be evaluated as well?}
  \item{interaction}{Add description}
  \item{guide_parm}{a vector of indices of the parameters (incl. intercept) for which estfun should be considered in chi-squared tests.}
  \item{guide_testtype}{character specifying whether a maximal selection ("max"), the summed up test statistic ("sum"), or COIN ("coin") should be employed.}
  \item{guide_decorrelate}{Add description}
  \item{xgroups}{integer. Number of categories for split variables to be employed in chi-squared tests (optionally breaks can be handed over).}
  \item{ygroups}{integer. Number of categories for scores to be employed in chi-squared tests (optionally breaks can be handed over).}
  \item{weighted.scores}{logical. Should scores be weighted in GUIDE}
  \item{...}{additional \code{ctree_control} arguments.}
}

\value{
  A list with components named as the arguments.
}

\seealso{\code{\link{ctree_control}}, \code{\link{disttree}}, \code{\link[partykit]{extree_fit}}}
