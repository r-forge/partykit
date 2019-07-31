context("C routines for kappa estimation")

## Test results
test_that("computed values are correct", {
  expect_equal(
    solve_kappa_Newton_Fourier(0.1,100,useC = FALSE),
    solve_kappa_Newton_Fourier(0.1,100,useC = TRUE)
  )
  expect_equal(
    .Call("Ginv", 0.1, 0.2, 0.3),
    movMF:::Ginv(0.1,0.2,0.3)
  )
  expect_equal(
    .Call("Rinv_lower_Amos_bound", 0.3, 0.3),
    movMF:::Rinv_lower_Amos_bound(0.3, 0.3)
  )
  expect_equal(
    .Call("Rinv_upper_Amos_bound", 0.3, 0.3),
    movMF:::Rinv_upper_Amos_bound(0.3, 0.3)
  )
  expect_equal(
    .Call("Rinv_upper_Amos_bound", 0.3, 2),
    movMF:::Rinv_upper_Amos_bound(0.3, 2)
  )
  expect_equal(
    .Call("A_PCF", 0.3, 1.1),
    movMF:::A(0.3, 1.1)
  )
  expect_equal(
    .Call("A_PCF", 0.3, 2),
    movMF:::A(0.3, 2)
  )
  expect_equal(
    .Call("Aprime_PCF", 0.1, 1.1, 0.3),
    movMF:::Aprime(0.1, 1.1, A = 0.3)
  )
  expect_equal(
    .Call("Aprime_PCF", 0.1, 0.1, 2),
    movMF:::Aprime(0.1, 0.1, A = 2)
  )
})

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

## Test performance
r <- rnorm(5)
fun <- function(r, ...) solve_kappa_Newton_Fourier(r, 100, ...)

n <- 100000
t_R <- system.time(sapply(rep(0.1,n), fun, useC = FALSE))
t_C <- system.time(sapply(rep(0.1,n), fun, useC = TRUE))


test_that("time for c-routines are faster", {
  expect_true(compare_times(t_R, t_C, output = "numeric") > -10)

})

#compare_times(t_R, t_C, "solve kappa")





