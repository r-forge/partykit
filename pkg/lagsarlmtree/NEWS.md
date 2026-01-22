# lagsarlmtree 1.0-2

* As the demo in the package leverages the `strucchange` package, this is
  declared now as a "Suggests" dependency in the DESCRIPTION.
  
* Further small improvements in documentation and links.


# lagsarlmtree 1.0-1

* The German Economic Review paper was finally published:
  Wagner and Zeileis (2019). "Heterogeneity and Spatial Dependence of
  Regional Growth in the EU: A Recursive Partitioning Approach."
  _German Economic Review_, **20**(1), 67-82.
  [doi:10.1111/geer.12146](https://doi.org/10.1111/geer.12146).
  The preprint version is still at:
  <https://www.zeileis.org/papers/Wagner+Zeileis-2019.pdf>

* Roger Bivand joined as a contributor to the package and adapted the
  code to transition from the old `spdep` to the new `spatialreg` package.
  Furthermore, impacts were added and extra arguments passed through to
  `lagsarlm()` from the `lagsarlmtree()` function.


# lagsarlmtree 1.0-0

* First CRAN release of the `lagsarlmtree` package accompanying the acceptance
  of the corresponding manuscript "Heterogeneity and Spatial Dependence of
  Regional Growth in the EU: A Recursive Partitioning Approach" in the
  _German Economic Review_ ([doi:10.1111/geer.12146](https://doi.org/10.1111/geer.12146)).
  A preprint version is available online at:
  <https://www.zeileis.org/papers/Wagner+Zeileis-2019.pdf>

* Replication materials for the manuscript are available in the package as:
  `demo("GrowthNUTS2", package = "lagsarlmtree")`

