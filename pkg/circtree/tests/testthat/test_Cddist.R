context("C routines for ddist")

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

dR <-  dist_vonmises(useC = FALSE, ncores = 1)
dC <-  dist_vonmises(useC = TRUE, ncores = 1)
dC2 <- dist_vonmises(useC = TRUE, ncores = 2)
dC5 <- dist_vonmises(useC = TRUE, ncores = 5)


## Test results
set.seed(111)
eta <- c(tan(0.8/2), log(0.8))
n <- 1000000
y <- rnorm(n)

dR_out <- dR$ddist(y, eta)
dC_out <- dC$ddist(y, eta)
C2_out <- dC2$ddist(y, eta)
C5_out <- dC5$ddist(y, eta)

test_that("scores (sum = FALSE) are the same", {
  expect_equal(
    dR_out,
    dC_out,
    tolerance = 1e-6
  )
  expect_equal(
    dC_out,
    C2_out,
    tolerance = 1e-6
  )
  expect_equal(
    C2_out,
    C5_out,
    tolerance = 1e-6
  )
})

dR_out <- dR$ddist(y, eta, sum = TRUE)
dC_out <- dC$ddist(y, eta, sum = TRUE)
C2_out <- dC2$ddist(y, eta, sum = TRUE)
C5_out <- dC5$ddist(y, eta, sum = TRUE)

test_that("scores (sum = TRUE) are the same", {
  expect_equal(
    dR_out,
    dC_out,
    tolerance = 1e-6
  )
  expect_equal(
    dC_out,
    C2_out,
    tolerance = 1e-6
  )
  expect_equal(
    C2_out,
    C5_out,
    tolerance = 1e-6
  )
})

dR_out <- dR$ddist(y, eta, log = FALSE)
dC_out <- dC$ddist(y, eta, log = FALSE)
C2_out <- dC2$ddist(y, eta, log = FALSE)
C5_out <- dC5$ddist(y, eta, log = FALSE)

test_that("scores (sum = TRUE) are the same", {
  expect_equal(
    dR_out,
    dC_out,
    tolerance = 1e-6
  )
  expect_equal(
    dC_out,
    C2_out,
    tolerance = 1e-6
  )
  expect_equal(
    C2_out,
    C5_out,
    tolerance = 1e-6
  )
})

## Test performance
t_R  <- system.time(dR$ddist(y, eta)) 
t_C  <- system.time(dC$ddist(y, eta)) 
t_C2 <- system.time(dC2$ddist(y, eta))
t_C5 <- system.time(dC5$ddist(y, eta))

test_that("time for c-routines are faster", {
  expect_true(compare_times(t_R, t_C, output = "numeric") > -200)  ## C slower than R using windows
  expect_true(compare_times(t_R, t_C2, output = "numeric") > -200)  ## C slower than R using windows
  expect_true(compare_times(t_R, t_C5, output = "numeric") > -200)  ## C slower than R using windows
})

#compare_times(t_R, t_C, "ddist one core")
#compare_times(t_R, t_C2, "ddist 2 cores")
#compare_times(t_R, t_C5, "ddist 5 cores")
