
# -------------------------------------------------------------------
# EXAMPLES FOR TESTING WITH TINYTEST
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Run some sample tests
# -------------------------------------------------------------------
expect_equal(1 + 1, 2)          # TRUE
expect_equivalent(2, c(x=2))    # TRUE


# -------------------------------------------------------------------
# Run some sample tests only at home, where env variable TT_AT_HOME=TRUE
# -------------------------------------------------------------------
if(at_home()){


  expect_silent(1+1)            # TRUE
  expect_silent(print("hihi"))  # TRUE, nothing goes to screen
}
