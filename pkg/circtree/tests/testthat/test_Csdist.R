context("C routines for sdist")

## Small helper function
compare_times <- function(t_R, t_C, label = "R vs C", output = c("text", "numeric")) {
    output <- match.arg(output)
    t_R <- as.numeric(t_R[3L], units = "secs")
    t_C <- as.numeric(t_C[3L], units = "secs")
    inc <- 100 * (t_R/t_C - 1);
    if(output == "text"){
      cat(sprintf("\nTiming for \"%s\"\n", label))
      cat(sprintf("   R:  %10.5f sek   C: %10.5f sek   (%4.1f percent %s)\n",
                  t_R, t_C, inc, ifelse(inc < 0, "slower in C", "faster in C")))
    } else {
      inc
    }
}

## Set up families functions
dR  <- dist_vonmises(useC = FALSE, ncores = 1)
dC  <- dist_vonmises(useC = TRUE, ncores = 1)
dC2 <- dist_vonmises(useC = TRUE, ncores = 2)
dC5 <- dist_vonmises(useC = TRUE, ncores = 5)

## Test results
set.seed(111) 
eta <- c(tan(0.8/2), log(0.8))
n <- 1000000
y <- rnorm(n)

dR_out <- dR$sdist(y, eta)
dC_out <- dC$sdist(y, eta)
C2_out <- dC2$sdist(y, eta)
C5_out <- dC5$sdist(y, eta)

test_that("scores (sum = TRUE) are the same", {
  expect_equal(
    dR_out,
    dC_out
  )
  expect_equal(
    dC_out,
    C2_out
  )
  expect_equal(
    C2_out,
    C5_out
  )
})

dR_out <- dR$sdist(y, eta, sum = TRUE)
dC_out <- dC$sdist(y, eta, sum = TRUE)
C2_out <- dC2$sdist(y, eta, sum = TRUE)
C5_out <- dC5$sdist(y, eta, sum = TRUE)

test_that("scores (sum = FALSE) are the same", {
  expect_equal(
    dR_out,
    dC_out
  )
  expect_equal(
    dC_out,
    C2_out
  )
  expect_equal(
    C2_out,
    C5_out
  )
})


## Test performance 
t_R  <- system.time(dR$sdist(y, eta))
t_C  <- system.time(dC$sdist(y, eta))
t_C2 <- system.time(dC2$sdist(y, eta))
t_C5 <- system.time(dC5$sdist(y, eta))

test_that("time for c-routines are faster", {
  expect_true(compare_times(t_R, t_C, output = "numeric") > -20)
  expect_true(compare_times(t_R, t_C2, output = "numeric") > -20)
  expect_true(compare_times(t_R, t_C5, output = "numeric") > -20)
})

#compare_times(t_R, t_C, "sdist one core")
#compare_times(t_R, t_C2, "sdist 2 cores")
#compare_times(t_R, t_C5, "sdist 5 cores")

