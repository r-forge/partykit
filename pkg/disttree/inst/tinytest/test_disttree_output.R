# -------------------------------------------------------------------
# - NAME:   test_disttree_output.R
# - AUTHOR: Moritz N. Lang
# - DATE:   2020-11-16
# -------------------------------------------------------------------
# - PURPOSE: Test disttree output for plausability using tinytest.
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# PRELIMINARIES
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# RUN TESTS
# -------------------------------------------------------------------

## Check class of disttree output
tr <- disttree(dist ~ speed, data = cars)
expect_true(all(class(tr) %in% c("disttree", "modelparty", "party")))

##TODO: (ML) Implement further tests for family etc.

