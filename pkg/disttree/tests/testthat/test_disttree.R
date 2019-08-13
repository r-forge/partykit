# -------------------------------------------------------------------
# - NAME:   test_disttree.R
# - AUTHOR: Moritz N. Lang, Lisa Schlosser
# - DATE:   2019-08-12
# -------------------------------------------------------------------
# - PURPOSE: Tests for plausibility of tree
# -------------------------------------------------------------------

context("test of plausibility of disttree()'")

## Check classes
tr <- disttree(dist ~ speed, data = cars)
expect_true(all(class(tr) %in% c("disttree", "modelparty", "party")))

##TODO: (ML) Implement further tests for family etc.
