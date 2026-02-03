
if (file.exists("packages.bib")) file.remove("packages.bib")
pkgversion <- function(pkg) {
    pkgbib(pkg)
    packageDescription(pkg)$Version
}
pkgbib <- function(pkg) {
    x <- citation(package = pkg, auto = TRUE)[[1]]
    b <- toBibtex(x)
    b <- gsub("Buehlmann", "B{\\\\\"u}hlmann", b)
    b[1] <- paste("@Manual{pkg:", pkg, ",", sep = "")
#    if (is.na(b["url"])) {
#        b[length(b)] <- paste("   URL = {http://CRAN.R-project.org/package=",
#                              pkg, "},", sep = "")
#    }
    b <- b[names(b) != "url"]
    if (is.na(b["doi"])) {
        b[length(b)] <- paste("   DOI = {10.32614/CRAN.package.",
                              pkg, "}", sep = "")
        b <- c(b, "}")
    }
    if (length(grep("commit", b["note"])) > 0)
        b["note"] <- paste0(gsub(", commit.*", "", b["note"]), "},")
    cat(b, sep = "\n", file = "packages.bib", append = TRUE)
}
pkg <- function(pkg)
    paste("\\\\pkg{", pkg, "} \\\\citep[version~",
          pkgversion(pkg), ",][]{pkg:", pkg, "}", sep = "")


pkgs <- c("mvtnorm", "partykit", "party", "psychotree", "C50", # "mvpart",
          "survival", "rpart", "modeltools", "mlbench", "pmml", "RWeka", "TH.data")
p <- sapply(pkgs, require, character.only = TRUE)
if (!all(p))
    try(install.packages(pkgs[!p], repos = "https://CRAN.R-project.org"))

sapply(pkgs, pkg)
