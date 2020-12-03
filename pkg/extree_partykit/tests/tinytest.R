## Trigger tiny tests within /inst/tinytest/
  
if(requireNamespace("tinytest", quietly=TRUE)){
  tinytest::test_package("partykitx")
}

