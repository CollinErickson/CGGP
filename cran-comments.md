I received an email on 8/19/23 from Kurt Hornik saying that I had to 
update the documentation for my package as specified in this bug:
https://github.com/r-lib/roxygen2/issues/1491. I made the requested change.
Apologies for taking so long.

I ran CMD CHECK on all the following and had no issues except for R-hub
Debian Linux. On that one there is a PREPERROR. Uwe said that this is okay
since it works on their tests. The issue only arises when I remove the CXX11
requirement.


## Test environments
* local Windows 11, R 4.3.2
* ubuntu 22.04 (on GitHub Actions), R 4.3.2
* win-builder (devel and release)
* R-hub (Windows, Ubuntu Linux, Fedora Linux, Debian Linux)


## R CMD check results

On personal Windows 11 (1/19/24) there are 0 notes/errors/warnings

On R-hubWindows Server, Fedora Linux, and Ubuntu Linux (1/20/24) there are a few
NOTEs, but no real issues.

On R-hub Devian Linux (1/20/24) there is a PREPERROR. Uwe said that this is okay.

On winbuilder (check_win_release, 1/20/24) it has "Status: OK".

On winbuilder (check_win_devel, 1/20/24) it has "Status: OK".

On GitHub Actions (1/20/24) it has "Status: OK".


## Downstream dependencies

There are currently no downstream dependencies for this package.
