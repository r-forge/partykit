# -------------------------------------------------------------------
# - NAME:   test_input_variablenames.R
# - AUTHOR: Moritz N. Lang
# - DATE:   2020-11-23
# -------------------------------------------------------------------
# - PURPOSE: Test non-standard variable names employing ctree() and mob().
# -------------------------------------------------------------------
# - REMARK: - Actually written for `partykit', needs to be adapted or moved.
#           - Tests are only run at home, where env variable TT_AT_HOME=TRUE.
#           - Perform comparison based on nodeapply(), tree$node, or print statement?
# -------------------------------------------------------------------

if(at_home()){

  # -------------------------------------------------------------------
  # PRELIMINARIES
  # -------------------------------------------------------------------
  ## Check if partykit is installed
  if(!require("partykit")) stop("Please install the package 'partykit' ...")

  ## Helper function to produce print statement
  get_print_statement <- function(print_statement){
    sink(tmp_file <- tempfile())
    on.exit(sink())
    invisible(force(print(print_statement)))

    con <- file(tmp_file, "r", encoding = "Latin1")
    output <- readLines(con); close(con)

    return(output)
  }

  
  # -------------------------------------------------------------------
  # RUN TEST FOR SPLIT VARIABLE NAMES WITH BACKTICKS
  # -------------------------------------------------------------------
  
  # -------------------------------------------------------------------
  # Check ctree for single space with cars dataset
  # -------------------------------------------------------------------
  ct_cars1 <- ctree(dist ~ speed, data = cars) 
  
  mycars <- cars
  names(mycars)[1] <- "my speed"
  ct_cars2 <- ctree(dist ~ `my speed`, data = mycars) 
  
  expect_equivalent(nodeapply(ct_cars1), nodeapply(ct_cars2))
  
  
  # -------------------------------------------------------------------
  # Check ctree for single space and log() with cars dataset
  # -------------------------------------------------------------------
  ct_cars3 <- ctree(dist ~ log(speed), data = cars) 
  
  mycars <- cars
  names(mycars)[1] <- "my speed"
  ct_cars4 <- ctree(dist ~ log(`my speed`), data = mycars) 
  
  expect_equivalent(nodeapply(ct_cars3), nodeapply(ct_cars4))
  
  
  # -------------------------------------------------------------------
  # Check ctree for several spaces with simulated data
  # -------------------------------------------------------------------
  my_data <- data.frame(
    Y = factor(rep(LETTERS[1:3], each = 10)),
    x1 = 1:30,
    x2 = c(1:10, 2:11, 3:12)
  )
  
  ct_sim1 <- ctree(Y ~ ., data = my_data)
  
  my_data2 <- my_data
  names(my_data2)[2] <- "x 1 2 3"
  
  ct_sim2 <- ctree(Y ~ ., data = my_data2)
  
  expect_equivalent(nodeapply(ct_sim1), nodeapply(ct_sim2))
  
  
  # -------------------------------------------------------------------
  # Check lmtree single space with cars dataset  
  # -------------------------------------------------------------------
  lmt_cars1 <- lmtree(dist ~ speed, data = cars) 
  
  mycars <- cars
  names(mycars)[1] <- "my speed"
  lmt_cars2 <- lmtree(dist ~ `my speed`, data = mycars) 
  
  expect_equivalent(nodeapply(lmt_cars1), nodeapply(lmt_cars2))
  
  
  # -------------------------------------------------------------------
  # Check mob for several spaces with Pima Indians diabetes data (official example)
  # -------------------------------------------------------------------
  if(require("mlbench")) {
  
    data("PimaIndiansDiabetes", package = "mlbench")
    
    logit <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
      glm(y ~ 0 + x, family = binomial, start = start, ...)
    }
    
    mob_pid1 <- mob(diabetes ~ glucose | pregnant + pressure + triceps + insulin +
      mass + pedigree + age, data = PimaIndiansDiabetes, fit = logit)
  
    PimaIndiansDiabetes2 <- PimaIndiansDiabetes
    names(PimaIndiansDiabetes2)[2] <- "glu co se"
    names(PimaIndiansDiabetes2)[1] <- "pre gnant"
    names(PimaIndiansDiabetes2)[7] <- "pedi gree"
    names(PimaIndiansDiabetes2)[6] <- "m a s s"
  
    mob_pid2 <- mob(diabetes ~ `glu co se` | `pre gnant` + pressure + triceps + 
      insulin + `m a s s` + `pedi gree` + age, data = PimaIndiansDiabetes2, 
      fit = logit)
  
    expect_equivalent(nodeapply(mob_pid1), nodeapply(mob_pid2))
  }
  
  # -------------------------------------------------------------------
  # RUN TEST FOR RESPONSE VARIABLE NAMES WITH BACKTICKS
  # -------------------------------------------------------------------
  
  # -------------------------------------------------------------------
  # Check ctree for several spaces with cars dataset
  # -------------------------------------------------------------------
  ct_cars1 <- ctree(dist ~ speed, data = cars) 
  
  my_cars2 <- cars
  names(my_cars2)[2] <- "d i s t"
  ct_cars5 <- ctree(`d i s t` ~ speed, data = my_cars2) 
  
  expect_equivalent(nodeapply(ct_cars1), nodeapply(ct_cars5))
  
  
  # -------------------------------------------------------------------
  # Check lmtree for several spaces with cars dataset
  # -------------------------------------------------------------------
  lmt_cars1 <- lmtree(dist ~ speed, data = cars) 
  lmt_cars3 <- lmtree(`d i s t` ~ speed, data = my_cars2) 
  
  #TODO: (ML) Compare nodeapply(), tree$node, or print statement?
  expect_identical(get_print_statement(lmt_cars1$node), 
                   get_print_statement(lmt_cars3$node))
} 
