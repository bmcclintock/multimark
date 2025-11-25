# multimark  [![R-CMD-check](https://github.com/bmcclintock/multimark/workflows/R-CMD-check/badge.svg)](https://github.com/bmcclintock/multimark/actions) [![CRAN_Downloads](https://cranlogs.r-pkg.org/badges/multimark)](https://cran.r-project.org/package=multimark) [![CRAN_Downloads](https://cranlogs.r-pkg.org/badges/grand-total/multimark)](https://cran.r-project.org/package=multimark)

R package for capture-recapture analysis with multiple non-invasive marks

## Installation instructions

### CRAN release
The package is available at [![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/multimark)](https://cran.r-project.org/package=multimark). To install it:
``` R
install.packages("multimark")
```

### Install from Github
To install the latest (stable) version of the package from Github: [![R-CMD-check](https://github.com/bmcclintock/multimark/workflows/R-CMD-check/badge.svg)](https://github.com/bmcclintock/multimark/actions)
``` R
library(remotes)
install_github("bmcclintock/multimark")
```

To install the latest (**unstable**) version of the package from Github: [![R-CMD-check](https://github.com/bmcclintock/multimark/actions/workflows/multimark.yml/badge.svg?branch=develop)](https://github.com/bmcclintock/multimark/actions/workflows/multimark.yml)
``` R
library(remotes)
install_github("bmcclintock/multimark@develop")
```

## References
McClintock, B.T. (2015) [multimark: an R package for analysis of capture–recapture data consisting of multiple “noninvasive” marks](https://doi.org/10.1002/ece3.1676). *Ecology and Evolution*, 5(21), 4920-4931.

Maronde, L., McClintock, B.T., Breitenmoser, U., Zimmermann, F. (2020) [Spatial capture–recapture with multiple noninvasive marks: An application to camera-trapping data of the European wildcat (*Felis silvestris*) using R package multimark](https://doi.org/10.1002/ece3.6990). *Ecology and Evolution*, 10(24), 13968–13979.
