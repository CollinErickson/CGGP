I received an email on 4/16/20 from Kurt Hornik saying that I had to add 
rmarkdown to Suggests in the DESCRIPTION because of a change to knitr.
I have made this change and a few other small changes to the package.

I ran CMD CHECK on all the following and had no issues.

## Test environments
* local Ubuntu 20.04.2 LTS install, R 4.0.3
* local Windows 10, R 4.0.5
* ubuntu 16.04.6 LTS (on travis-ci), R 4.0.2
* win-builder (devel and release)
* R-hub (Windows, Ubuntu Linux, Fedora Linux, Debian Linux)

## R CMD check results

On personal Ubuntu and Windows there are no ERRORs, WARNINGs, or NOTEs.

On R-hub all 4 have "Status: OK".

On winbuilder (check_win_release and check_win_devel) it has "Status: OK".

## Downstream dependencies

There are currently no downstream dependencies for this package.
