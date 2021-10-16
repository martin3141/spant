
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Spectroscopy Analysis Tools (spant) <img src="man/figures/logo.png" align="right" width=130/>

[![R build
status](https://github.com/martin3141/spant/workflows/R-CMD-check/badge.svg)](https://github.com/martin3141/spant/actions)
[![](http://cranlogs.r-pkg.org/badges/spant)](http://cran.rstudio.com/web/packages/spant/index.html)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/spant)](https://cran.r-project.org/package=spant)
[![Coverage
Status](https://coveralls.io/repos/github/martin3141/spant/badge.svg?branch=master)](https://coveralls.io/github/martin3141/spant?branch=master)

## Overview

spant provides a full suite of tools to build automated analysis
pipelines for Magnetic Resonance Spectroscopy (MRS) data. The following
features and algorithms are included:

-   Advanced fully-automated metabolite fitting algorithm - ABfit
    <https://onlinelibrary.wiley.com/doi/10.1002/mrm.28385>.
-   Robust retrospective frequency and phase correction - RATS
    <https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.27605>.
-   Flexible data types to support single voxel, dynamic and
    spectroscopic imaging data types.
-   Raw data import from individual coils and dynamic measurements, eg
    support for importing individual FIDs from Siemens TWIX formatted
    data.
-   Publication quality plotting.
-   Extensive set of pre-processing steps (phasing, coil-combination,
    zero-filling, HSVD filtering…)
-   Quantum mechanical based simulation for experimental design and
    basis-set generation.
-   Set of metabolite, macromolecule and lipid parameters for typical
    brain analyses.
-   Voxel registration to anatomical images for partial volume
    concentration corrections.

## Installation

Download and install the latest version of R
(<https://cloud.r-project.org/>), or with your package manager if using
a recent Linux distribution, eg `sudo apt install r-base`.

It is also strongly recommended to install RStudio Desktop
(<https://rstudio.com/products/rstudio/download>) to provide a modern
environment for interactive data analysis.

Once R and RStudio have been installed, open the RStudio application and
type the following in the Console (lower left panel) to install the
latest stable version of spant:

``` r
install.packages("spant", dependencies = TRUE)
```

Or the the development version from GitHub (requires the `devtools`
package):

``` r
install.packages("devtools")
devtools::install_github("martin3141/spant", ref = "devel")
```

### Ubuntu 20.04 installation

CRAN packages need to be compiled on Linux, and therefore you may need
to ensure some additional system libraries are installed. spant may be
installed from a clean installation of Ubuntu 20.04 with the following
commands pasted into the terminal:

``` ubuntu
sudo apt install -y r-base libcurl4-openssl-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libtiff5-dev
mkdir -p ~/R/x86_64-pc-linux-gnu-library/3.6
Rscript -e 'install.packages("spant", dependencies = TRUE)'
```

### Ubuntu 21.10 installation

``` ubuntu
sudo apt install -y r-base libcurl4-openssl-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libtiff5-dev
mkdir -p ~/R/x86_64-pc-linux-gnu-library/4.0
Rscript -e 'install.packages("spant", dependencies = TRUE)'
```

## Documentation

Quick introduction to the basic analysis workflow :
<https://martin3141.github.io/spant/articles/spant-intro.html>

Short tutorials : <https://martin3141.github.io/spant/articles/>

Function reference : <https://martin3141.github.io/spant/reference/>

Once the spant library has been loaded with `library(spant)`, type
`?spant` on the console for instructions on how to access the offline
documentation. Note that offline help on the available functions can be
quickly shown in RStudio using `?function_name`, eg `?read_mrs`.
