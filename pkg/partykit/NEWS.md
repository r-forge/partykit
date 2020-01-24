# partykit 1.2-6

* Fix warning in predict.cforest
* Don't rely on suggested packages in man and tests

# partykit 1.2-5

### Bugfixes

* Trying to split in a variable where all observations were missing
          nevertheless produced a split, as reported by Kevin Ummel.


# partykit 1.2-4

### Bugfixes

* Update reference output, fix `RNGversion`.


# partykit 1.2-3

### New Features

* `varimp()` runs in parallel mode, optionally.
* `weights` in `cforest()` can now be used
	to specify a matrix of weights (number of observations
	times number of trees) to be used for tree induction (this was
	always possible in `party::cforest`. This was advertised in
	the documentation but actually not implemented so far.

### Bugfixes

* `predict()` did not pay attention to `xlev`; this caused
	problems when empty factor levels were removed prior to tree fitting. See
       <https://stackoverflow.com/questions/51043870/lmtree-suspect-behavior-with-factors>.
* `nodeprune()` may have got fitted terminal node numbers wrong,
	spotted by Jason Parker.


# partykit 1.2-2

### Bugfixes

* In `mob()` using the `cluster` argument with a `factor`
    variable sometimes lead to `NA`s in the covariance matrix estimate
    if empty categories occured in subgroups. The problem had been introduced
    in version 1.2-0 and has been fixed now.
* Methods for the `sctest()` generic from the `strucchange`
    package are now dynamically registered if `strucchange` is attached.
    Alternatively, the methods can be called directly using their full names
    `sctest.constparty()` and `sctest.modelparty()`.
* The `prune.modelparty()` function is now fully exported but
    it is also registered with the `prune` generic from `rpart`.


# partykit 1.2-1

### New Features

* New `scale` argument for `predict()` in `cforest()`.
	For simple regression forests, predicting
	the conditional mean by nearest neighbor weights with
	`scale = TRUE` is now
	equivalent to the aggregation of means. The unscaled version
	proposed in <doi:10.1002/sim.1593> can be obtained with
	`scale = FALSE`.

### Bugfixes

* Bug fix for case weights in `mob()` in previous version (1.2-0)
    introduced a bug in the handling of proportionality weights. Both cases
    are handled correctly now.
* `glmtree()` can now handle `caseweights = TRUE` correctly
    for `vcov` other than the default `"opg"`. Internally,
    the `glm` objects are adjusted by correcting the dispersion
    estimate and the degrees of freedom.
* `lookahead` did not work in the presence of missing values.
* Calling `partykit::ctree` did not work when partykit was not
    attached.
* `node_inner()` now allows to set a different `gpar(fontsize = ...)`
    in the inner nodes compared to the overall tree.
* `splittest` asked for Monte-Carlo p-values, even
	when the test statistic was used as criterion.


# partykit 1.2-0

### New Features

* We welcome Heidi Seibold as co-author!
* Internal re-organisation for `ctree()` by means of new extensible
	tree infrastructure (available in  `extree_data()` and 
	`extree_fit()`). Certain parts of the new infrastructure
	are still experimental. `ctree()` is fully backward
	compatible.
* Use `libcoin` for computing linear test statistics 
	and p-values for `ctree()`.
* Use `inum` for binning (the new `nmax` argument).
* Quadratic test statistics for splitpoint selection are now available for `ctree()` via 
	`ctree_control(splitstat = "quadratic")`.
* Maximally selected test statistics for variable selection are now available for `ctree()` via 
	`ctree_control(splittest = TRUE)`.
* Missing values can be treated as a separate category, also
	for splits in numeric variables in `ctree()` via
	`ctree_control(MIA = TRUE)`.
* Permutation variable importance, including conditional variable
	importance, was added to `partykit`.
* New `offset` argument in `ctree()`.
* New `get_paths` for computing paths to nodes.
* `node_barplot` gained a `text` argument
	that can be used to draw text labels for the
	percentages displayed.
* The `margins` used in `plot.party()` can now
	also be set by the user.

### Bugfixes

* Bug fix in `mob()` if `weights` are used and `caseweights = TRUE` (the default).
	The statistics for the parameter instability tests were computed incorrectly and
        consequently the selection of splitting variables and also the stopping criterion
        were affected/incorrect.` 
* Avoid log(p) values of `-Inf` inside `mob()` by replacing weighted averaging
	with naive averaging in the response surface regression output in case
        the p values are below machine precision.
* The `as.party()` method for `rpart` objects
	without any splits only returned a naked `partynode`
        rather than a full `party`. This has been corrected
        now.
* `nodeapply` did not produce the same results for permutations of `ids`.
	Spotted by Heidi Seibold.
* Out-of-bag predictions in `predict.cforest()` were incorrect.
* `perm` in `predict` was only considered when
	`newdata` was given. Spotted by Heidi Seibold.
* Don't try to search for binary splits in unordered factors with more than 31
	levels. This potentially caused an integer overrun in previous versions.
	`party::ctree()` uses an approximation for binary split
	searches in unordered factors; thus, using `party` might be an
	alternative.


# partykit 1.1-1

* Proper support of quasi-families in `glmtree()`
	and hence `palmtree()`.
* NA handling by following the majority was 
	potentially incorrect in `ctree()`.
* Minor speed improvements.
* Breaking ties before variable selection was
	suboptimal for very small log-p-values.


# partykit 1.1-0

* Added a new function `palmtree()` that fits partially
    additive (generalized) linear model trees. These employ
    model-based recursive partitioning (`mob()`) based on (generalized)
    linear models with some local (i.e., leaf-specific) and
    some global (i.e., constant throughout the tree) regression
    coefficients.
* Splits in ordinal variables are now represented correctly 
    in the (still internal) `list.rule` method.
* Kaplan-Meier curves in `"constparty"` trees were
    plotted incorrectly due to use a wrong scaling of the
    x-axis. Spotted by Peter Calhoun <calhoun.peter@gmail.com>.
* Use `quote(stats::model.frame)` instead of
    `as.name("model.frame")`.
* The `as.party()` methods for `rpart` and
    `Weka_tree` now have a `data = TRUE` argument
    so that by default the data is preserved in the `party`
    object (instead of an empty model frame).
* The `predict` method for `cforest` objects
    did not work for one-row data frames, fixed now.
* Added `rot` and `just` arguments to `node_barplot`
    for more fine control of x-axis labeling (e.g., with 45 degree
    rotation).


# partykit 1.0-5

* The `partykit` package has now been published in
    _Journal of Machine Learning Research_, 16, 3905-3909.
    <http://jmlr.org/papers/v16/hothorn15a.html>
* Added support for setting background in panel functions.
* The `as.list()` method for `partynode` objects
    erroneously created an object `thisnode` in the calling
    environment which is avoided now.


# partykit 1.0-4

* Bug fix in `plot()` method for `constparty` objects.
	In the previous partykit version clipping was accidentally 
        also applied to the axes labels.
* For `constparty` objects `plot(..., type = "simple")`
	did not work correctly whereas `plot(as.simpleparty(...))`
        yielded the desired visualization. Now internally `as.simpleparty()`
        is called also in the former case.
* The `as.simpleparty()` method now preserves p-values from
	`constparty` objects (if any).
* Added a `getCall()` method for `"party"` objects.
* In the `predict()` method for `"lmtree"` and `"glmtree"`
	objects the `offset` (if any) was sometimes ignored. It is
        now always used in the prediction.


# partykit 1.0-3

* Import `logrank_trafo` from `coin`.


# partykit 1.0-2

* `nodeprune(..., ids = 1)` did not prune the tree to the root
    node. Fixed now.
* `predict.cforest()` used `na.omit` instead of `na.pass`.
* `predict.party()` now features new `perm` argument for
        permuting splits in specific variables (useful for computing
        permutation variable importances).
* `NAMESPACE` updates.


# partykit 1.0-1

* The support for (generalized) linear model trees with just
    a constant regressor has been improved. Now `lmtree(y ~ x1 + x2)`
    is short for `lmtree(y ~ 1 | x1 + x2)`, analogously for `glmtree()`.
    Plotting now also works properly in this case.
* The `as.party()` method for `"rpart"` objects did not
    work if one of the partitioning variables was a `"character"`
    variable rather than a `"factor"`. A suitable work-around has
    been added.
* The `node_barplot()` panel function can now also be used
    for multivariate responses, e.g., when all responses are numeric and
    on the same scale.
* The package now also includes a new data set `HuntingSpiders`
    which is essentially a copy of the `spider` data from the package
    `mvpart` that is currently archived on CRAN. The documentation has
    been improved somewhat and is likely to change further to explain how
    the data has been transformed in De'ath (2002).
* The survival tree example for the GBSG2 data was broken due to
    the response being (incorrectly) also part of the explanatory variables.
    Fixed by using the latest `Formula` package (at least version 1.2-1).


# partykit 1.0-0

* Version 1.0-0 published. This version is described in the MLOSS paper 
	accepted for publication by the Journal of Machine Learning Research today.
* The unbiased version of `cforest()` (with `replace = FALSE`) is 
	now the default (as in `party()`).
* Register all S3 methods in `NAMESPACE`.


# partykit 0.8-4

* Extended `mob()` interface by a `cluster` argument. This can be
  a vector (numeric, integer, factor) with cluster IDs that are
  then passed on to the 'fit' function (if supported) and used
  for clustering the covariance matrix in the parameter stability
  tests. `lmtree()` and `glmtree()` hence both gained a `cluster`
  argument which is used only for cluster covariances but not
  for the model estimation (i.e., corresponding to a working
  independence model).
* Optionally, the parameters' variance-covariance matrix in `mob()`
  can now be estimated by the sandwich matrix instead of the default
  outer-product-of-gradients (OPG) matrix or the information matrix.
* Reimplementation of `cforest()` available with extended 
  prediction facilities. Both the internal representation and the user interface
  are still under development are likely to change in future versions.
* Added multicore support to `mob()`, `ctree()`, and `cforest()`.
  If control argument `cores` is specified (e.g., `cores = 4`) then the
  search for the best variable or split point (often involving numerous model fits in
  `mob()` or resampling in `ctree()`) is carried out using `parallel::mclapply()`
  rathern than sequential `for()` or `sapply()`. Additionally, other
  `applyfun`s can be provided, e.g., using networks of workstations etc.
* Bug fix in `mob()` that occurred when regressor variables and
  partitioning variables overlapped and were not sorted in the
  underlying model frame.


# partykit 0.8-3

* `mvpart` was archived 2014-12-15.


# partykit 0.8-2

* Fixed an uninitialized memory issue reported by valgrind.


# partykit 0.8-1

* partykit now depends on R version >= 3.1.0 in order to import the
        `depth()` generic from the `grid` package.
* The print methods for `party`/`partynode` objects with only a root node
    was modified. Now, the terminal panel function is also applied
    if there is only a root node (while previously it was not).
*  `ctree()` now catches `sum(weights) <= 1` situations before they 
    lead to an error.
*  Code from suggested packages is included by using `::` syntax as
    required by recent R versions.
*  Argument `ytrafo` of `ctree()` can now be a function which will be
    updated in every node.
*  A small demo briefly illustrating some memory and speed properties
    has been added. It can be run interactively via
    `demo("memory-speed", package = "partykit").`
*  Section 3 of the "constparty" vignette now shows how properties of
    a new tree algorithm can be assessed by using `partykit` building
    blocks.


# partykit 0.8-0

* Major improved version of `partykit`. The previously existing functions
    in the package were tested and enhanced, new functions and
    extensive vignettes added.
* Extended and improved introductory documentation. The basic classes
    and class constructors `partynode()`/`partysplit()`/`party()` are introduced in 
    much more detail now in `vignette("partykit", package = "partykit")`.
* The class `constparty` (inheriting from `party`) for representing `party`
    objects with constant fits in the nodes (along with coercion methods
    for `rpart`, `J48`, etc.) is now described in more detail in the new
    `vignette("constparty", package = "partykit")`.
* The package now includes a reimplementation of the model-based
    recursive partitioning algorithm (MOB) using `partykit` infrastructure.
    The generic algorithm, the corresponding convenience interfaces
    `lmtree()` and `glmtree()` as well as various illustrations and possible
    extensions are described in detail in the new
    `vignette("mob", package = "partykit")`.
* Improved implementation of conditional inference trees (CTree), see
    the new `vignette("ctree", package = "partykit")` for details.
* New `nodeprune()` generic for pruning nodes in all `party` trees and
    `partynode` objects.
* Deal with empty levels in `ctree()` for `teststat = "quad"`
    (bug reported by Wei-Yin Loh <loh_at_stat.wisc.edu>).
* In `predict()` method for `constparty` objects, `type = "prob"` now returns
    ECDF for numeric responses and `type = "response"` the (weighted) mean.
* New panel function `node_ecdf()` for plotting empirical cumulative
    distribution functions in the terminal leaves of `constparty` trees.

# partykit 0.1-6

* Bug fix in `as.party()` method for J48 trees with ordered factors.


# partykit 0.1-5

* Fix C code problems reported by clang under OS X.

# partykit 0.1-4

* Added `node_surv()` for plotting survival ctrees. Accompanying
    infrastructure for survival trees was enhanced.
* `ctree()` now checks for (and does not allow) `x >= max(x)` splits.

# partykit 0.1-3

* Added `ipred` to the list of suggested packages due to usage of
        GlaucomaM and GBSG2 data in tests/examples.


# partykit 0.1-2

* The `node_terminal()` panel-generating function is now customizable
    by a `FUN` argument that is passed to `formatinfo()`.
* The `plot()` method for `simpleparty` object now sets up a formatting
    function passed to `formatinfo()`, both in `print()` and `plot()`.
* Fixed bug in `pmmlTreeModel()` for processing label IDS in splits when
    not all levels are present.
* Cleaned up unused variables in C code and partial argument matching
    in R code.

# partykit 0.1-1

* First CRAN release.
* See `vignette("partykit", package = "partykit")` for a (somewhat rough)
    introduction to the package and its classes/methods.
