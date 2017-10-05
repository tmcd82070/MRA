## Resubmission comments

*Comment: -> So apparently Bryan Manly is at least authors if not also copyright holder here? Why is his name not mention in the Authors@R field?*

* Added Bryan Manly to contributors list 
* Bryan Manly does not hold copyright on any routines used here.  Little of Bryan's work remains.  In addition, Bryan worked for the same company as the maintainer (Trent McDonald) when developing the routines, and this company released the code in approximately 2005.  

*Comment: Is anything from Numerical recipes still included? If so you probably must not use GPL and must not diestribute the package given the Numerical Recipes license?*

* This comment was very old.  
* No Numerical Recipies code is included in the current package. 


## Current Version
2.16.10
Previous submissions: 2.16.9

## Test environments
* local Windows 10 Pro (Version	10.0.15063 Build 15063), R 3.4.2
* Linix - Ubuntu (via Travis)
* passes devtools::build_win() with the following note:

*checking CRAN incoming feasibility ... NOTE
Maintainer: 'Trent McDonald <tmcdonald@west-inc.com>'

New submission

Package was archived on CRAN

Possibly mis-spelled words in DESCRIPTION:
  Horvitz (18:153)
  Huggin's (17:248)
  Seber (17:111, 17:208)
  covariates (17:56)
  logit (18:25)

CRAN repository db overrides:
  X-CRAN-Comment: Archived on 2017-04-24 as check errors were not
    corrected despite reminders.*

## R CMD check --as-cran results
No ERRORs or WARNINGs.  Same note as above. 

## Downstream dependencies
None known
