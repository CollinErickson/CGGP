
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SGGP

[![Travis-CI Build
Status](https://travis-ci.org/CollinErickson/SGGP.svg?branch=master)](https://travis-ci.org/CollinErickson/SGGP)
[![Coverage
Status](https://img.shields.io/codecov/c/github/CollinErickson/SGGP/master.svg)](https://codecov.io/github/CollinErickson/SGGP?branch=master)

The goal of SGGP is to provide a sequential design of experiment
algorithm that can efficiently use many points and interpolate exactly.

## Installation

You can install SGGP from github with:

``` r
# install.packages("devtools")
devtools::install_github("CollinErickson/SGGP")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
## basic example code
library(SGGP)
d <- 8
SG = SGcreate(rep(0, d), rep(1, d),201)
```
