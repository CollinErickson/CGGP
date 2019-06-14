This was accepted to CRAN for the first time yesterday, but I found a
significant mistake in one of the default parameters for one of the
functions. I made the one change and am resubmitting since we found
this to have a significant impact on our model quality. In the future,
releases will be weeks/months apart as suggested by CRAN.

## Test environments
* local Windows 7 install, R 3.6.0
* ubuntu 14.04 (on travis-ci), R 3.6.0
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes


## Downstream dependencies

There are currently no downstream dependencies for this package.
