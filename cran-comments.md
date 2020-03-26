I received an email on 3/19/20 saying I had a broken donttest example.
I fixed the error shown in the email, then I fixed a second
error I found running R CMD with --run-donttest.
I also fixed another test and some documentation.

I reran tests on Travis, win-builder (devel and release),
and R-hub.

## Test environments
* local Ubuntu 18.04 install, R 3.6.3
* ubuntu 16.04 (on travis-ci), R 3.6.2
* win-builder (devel and release)
* R-hub (Windows, Ubuntu Linux, Fedora Linux)

## R CMD check results

0 errors | 0 warnings | 0 notes


## Downstream dependencies

There are currently no downstream dependencies for this package.
