# examplepkg

A tiny example R package set up for CRAN submission workflow.

## Install (development)
```r
# install.packages("devtools")
devtools::install()
```

## Example
```r
library(examplepkg)
hello("you")
```

## Development / Checks
```r
# install.packages(c("devtools","roxygen2","testthat","rhub"))
devtools::document()        # generates NAMESPACE and man/ .Rd files from roxygen
devtools::test()
devtools::check()           # runs R CMD check locally
rhub::check_for_cran()      # optional: remote checks on multiple platforms
devtools::build()           # creates source tarball for CRAN
```
