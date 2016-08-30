setwd("C:/Users/tobii/Desktop/swReg")
library(devtools)
library(roxygen2)


# process documentation:
setwd("./swReg")
document()

setwd("..")
check("swReg")
uninstall("swReg")
build("swReg")
install("swReg")





## create package:
#create("upre")

# process documentation:
setwd("./upre")
document()

showPDFmanual <- function(package, lib.loc=NULL)
{
  path <- find.package(package, lib.loc)
  system(paste(shQuote(file.path(R.home("bin"), "R")),
               "CMD", "Rd2pdf",
               shQuote(path)))
} 
showPDFmanual("upre")

# check and build
setwd("..")
check("upre")
uninstall("upre")
build("upre")
install("upre")



# It is typically a design mistake to use ::: in your code since the corresponding object has probably 
# been kept internal for a good reason. Consider contacting the package maintainer if you feel the need 
# to access the object for anything but mere inspection. 