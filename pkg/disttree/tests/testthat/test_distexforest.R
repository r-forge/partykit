# -------------------------------------------------------------------
# - NAME:   test_distexforest.R
# - AUTHOR: Moritz N. Lang, Lisa Schlosser
# - DATE:   2019-08-12
# -------------------------------------------------------------------
# - PURPOSE: Tests for plausibility of forest
# -------------------------------------------------------------------

context("test of plausibility of distexforest()'")

## Check classes
tr <- distexforest(dist ~ speed, data = cars)
expect_true(all(class(tr) %in% c("distexforest", "constparties", "parties")))

##TODO: (ML) Implement further tests for varimp() etc.


