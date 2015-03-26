## Resubmission comments
* Changed title in DESCRIPTION to title-case, as required.
* I have read the current version of CRAN Repository Policy
(Version $Revision: 3298$).  I agree to them and believe I 
have complied with all.  

## Current Version
2.16.3
Previous submission: 2.16.1

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