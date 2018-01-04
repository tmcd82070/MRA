## Resubmission comments

Version 2.16.11 fixes the following warnings encountered 
by version 2.16.10 on Linus-like machines

`
* installing *source* package ‘mra’ ...
** package ‘mra’ successfully unpacked and MD5 sums checked
** libs
gfortran  -fPIC -Wall -g -O2  -c  Mrawin.f90 -o Mrawin.o
Mrawin.f90:226:4:

     ptr_np => np
    1
Warning: Pointer at (1) in pointer assignment might outlive the pointer target [-Wtarget-lifetime]
Mrawin.f90:230:4:

     ptr_first => first
    1
Warning: Pointer at (1) in pointer assignment might outlive the pointer target [-Wtarget-lifetime]
Mrawin.f90:231:4:

     ptr_last => last
    1
Warning: Pointer at (1) in pointer assignment might outlive the pointer target [-Wtarget-lifetime]
Mrawin.f90:245:4:

     ptr_dead => idead
    1
Warning: Pointer at (1) in pointer assignment might outlive the pointer target [-Wtarget-lifetime]
Mrawin.f90:397:4:

     ptr_np => np
    1
Warning: Pointer at (1) in pointer assignment might outlive the pointer target [-Wtarget-lifetime]
Mrawin.f90:401:4:

     ptr_first => first
    1
Warning: Pointer at (1) in pointer assignment might outlive the pointer target [-Wtarget-lifetime]
`

## Current Version
2.16.11
Previous submissions: 2.16.10

## Test environments
* local Windows 10 Pro (Version	10.0.15063 Build 15063), R 3.4.2
* Linux - Ubuntu (via Travis)

## devtools::build_win()

* passes devtools::build_win() with the following note:

`
checking CRAN incoming feasibility ... NOTE
Maintainer: 'Trent McDonald <tmcdonald@west-inc.com>'

New submission

Package was archived on CRAN

Possibly mis-spelled words in DESCRIPTION:
  Horvitz (18:153)
  Huggin's (17:248)
  Seber (17:111, 17:208)

CRAN repository db overrides:
  X-CRAN-Comment: Archived on 2017-12-19 as check problems were not
    corrected despite reminders.
`

## R CMD check --as-cran results
No ERRORs or WARNINGs.  Same note as above. 

## Downstream dependencies
None known
