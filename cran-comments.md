I received an email on 8/19/23 from Kurt Hornik saying that I had to 
update the documentation for my package as specified in this bug:
https://github.com/r-lib/roxygen2/issues/1491.

I ran CMD CHECK on all the following and had no issues.

## Test environments
* local Ubuntu 20.04.2 LTS install, R 4.0.3
* local Windows 10, R 4.3.2
* ubuntu 16.04.6 LTS (on travis-ci), R 4.0.2
* win-builder (devel and release)
* R-hub (Windows, Ubuntu Linux, Fedora Linux, Debian Linux)

## R CMD check results

On personal Ubuntu and Windows there are no ERRORs, WARNINGs, or NOTEs.

On R-hub all 4 have "Status: OK".
Debian Linux (1/14/24 OK)
Ubuntu Linux and Fedora Linux have a NOTE for a slow example.

On winbuilder (check_win_release and check_win_devel) it has "Status: OK".

## Downstream dependencies

There are currently no downstream dependencies for this package.
