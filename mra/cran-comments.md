## Test environments
* local Windows 7 Pro, R 3.1.3
* passes devtools::build_win()

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Trent McDonald <tmcdonald@west-inc.com>'
New submission
Package was archived on CRAN
CRAN repository db overrides:
  X-CRAN-Comment: Archived on 2015-03-03 as misuse of \donttest was not
corrected.

All \donttest have been removed. Only reason for past use of \donttest was 
long examples; but, examples run in ~ 5 seconds on my laptop. 

## Downstream dependencies
None known, but devtools:revdep_check() does not work because the previous 
version has been archived. 

## SRC 
Contains Fortran code.  Upon install, several *.mod (Fortran module) files 
are generated.  I have not included these module files. 