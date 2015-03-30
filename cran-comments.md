## Resubmission comments
* Fixed Fortran bound error from last previous version. 
* Examples now pass bound checking. 

## Current Version
2.16.4
Previous submissions: 2.16.1, 2.16.3

## Test environments
* local Windows 7 Pro, R 3.1.3
* Linix - Ubuntu (via Travis)
* passes devtools::build_win()

## R CMD check results
No ERRORs or WARNINGs. 

Bound checking turned on this time.  FFLAGS in Makeconf were
   -O3 -mtune=core2 -fbounds-check

Examples run in ~ 5 seconds on my laptop. 

## Downstream dependencies
None known
