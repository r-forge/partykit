## Trigger tiny tests within /inst/tinytest/

if(requireNamespace("tinytest", quietly=TRUE)){
##  home <- identical(Sys.getenv("TT_AT_HOME"), "TRUE")  ## To trigger all tests only at home.
  tinytest::test_package("disttree")
}
