context("comparison of new implementation 'distextree()' with older version 'disttree()'")

tr <- distextree(dist ~ speed, data = cars)
expect_true(all(class(tr) %in% c("modelparty", "party")))

