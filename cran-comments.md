## Resubmission comments
* Fixed Fortran bound error from previous version. 

## Current Version
2.16.4
Previous submissions: 2.16.1, 2.16.3

## Test environments
* local Windows 7 Pro, R 3.1.3
* Linix - Ubuntu (via Travis)
* passes devtools::build_win()

## R CMD check results
No ERRORs or WARNINGs. 

One NOTE: 
  Maintainer: 'Trent McDonald <tmcdonald@west-inc.com>'
  New submission
  Package was archived on CRAN
  CRAN repository db overrides:
    X-CRAN-Comment: Archived on 2015-03-03 as misuse of \donttest was not
      corrected.
      
  Response: \donttest has been removed from all examples 

Bound checking was on during R CMD check.  FFLAGS in Makeconf were
   -O3 -mtune=core2 -fbounds-check

Examples run in ~ 5 seconds on my laptop. 


## Downstream dependencies
None known
