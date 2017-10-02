!	
!	MRAWIN
!
!	This file contains routines to build a DLL which is callable from
!	S and R to do Capture-recapture analysis.  This was based on the
!	previous versions of MRAWIN.
!
!	This is basically Version 5 of MRAWIN.
!
!	Trent McDonald created this version, from MRAWIN4, by updating
!	Bryan Manly's Fortran 77 code to Fortran 95. This version was
!	motivated by convergence problems that Eric Reghr found in MRAWIN4.
!
!	Compile this with > lf95 mrawin.f95 -dll -ml msvb
!
!	tlm - 6dec04 - named version 5
!
!   Trent McDonald modified this so that it could be called from both R and S-plus.
!   This required removing all characters from argument list.  Character strings
!   are passed differently by R and S .Fortran.
!   At that time Trent also converted to MS Visual Studio and Fortran .NET.
!   Now, to compile the DLL, open the 'solution' mrawin in Visual Studio,
!   change the 'solution configuration' to 'Release', and hit the 'start' icon.
!   What you are doing is building the solution.  Then copy the .dll from the bin
!   directory to where you want.
!
!   The compile > 'lf95 mrawin.f95 -dll -ml msvb' works for S-Plus, but not R.
!
!   tlm - 1July05 - named version 5.1
!
!   I switched to using the free G95 compiler to produce the DLL.  This compiler does
!   not give the "Warning: 0000F changed..." message that LF95 DLL's do.
!
!   To compile with G95> g95 -shared -o mrawin.dll mrawin.f95
!
!   tlm - 31 Jan 07
!
!   But, now that I am using Numerical Recipies to do the SVD, compile with
!   g95 -shared -o mrawin.dll mrawin.f95 nr.f90 nrutil.f90 nrtype.f90
!
!   tlm - 21 Feb 07
!
!   Now I am converted over to GFORTRAN so I could make this a package and post on
!   CRAN.  To make the package, I had to include NR code in this file.
!   i.e., we could not get the multiple files to compile on the CRAN side.
!   So, the GFORTRAN statement that makes the
!   DLL is :
!   gfortran -shared -o mra.dll mrawin.f95
!
!   The above does not give any warnings.  Use the following to see warnings:
!   gfortran -Wall -Wextra -shared -o mra.dll mrawin.f90
!
!   tlm - 18 oct 07
!
!   Added the Huggins model routines.
!   This is version 5.4 of MRAWIN
!
!   tlm - 29 nov 07
!
!   Version 2.3 of MRA: Added 'trace' parameter to keep MRA from creating a log file unless it is needed. 
!
!   For R, this code must compile on multiple platforms.  gfortran and 
!   the real(kind=x) syntax works on Windows, but does not work on other
!   machines.  Bryan Ripley suggested using the old fashioned 'double precision' 
!   rather than real(kind).  So I am.  However, for Numerical Recipies, 
!   I am not sure this will work, so I am leaving in some of the original NR code, but 
!   commented out.  In the remainder of this program, I changed 
!   the following: 
!        'real(sp)' to 'real' 
!        'real(dp)' and 'real(kind=dbl)' to 'double precison'.  
!        'integer(i4b)' to 'integer'
!   To convert back, global search and replace in reverse.
!
!   tlm - 6-nov-09 and 30-jan-10
!
!   Eric Reghr was doing simulations with p's and phi's close to 1.  The logit link was 
!   having issues with the asymptote there.  He asked that I implement the sine link.
!   so I did.  But, I called it the cosine link. 
!
!   tlm 21-feb-2010 (while in Kruger NP, South Africa)



!   ----- These are all the Numerical Recipie routines I need.
!	It seems that I need to include these routines in one file for CRAN
!	to compile.  We tried many different combinations of including and
!	compilation order, but finally punted.
!
!	For local construction of a binary zip package, these 3 include statements worked fine. Must be in this order.
!	But, now we don't need them.
!include 'nrtype.f90'
!include 'nrutil.f90'
!include 'nr.f90'

! ------------------------------
MODULE nrtype
    !INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
    !INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
    !INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
    !INTEGER, PARAMETER :: SP = KIND(1.0)
    !INTEGER, PARAMETER :: DP = KIND(1.0D0)
    !INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
    !INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
    !INTEGER, PARAMETER :: LGT = KIND(.true.)

    !REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_sp
    !REAL(SP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_sp
    !REAL(SP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_sp
    !REAL(SP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_sp
    !REAL(SP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_sp
    !REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_dp
    !REAL(DP), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858_dp
    !REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp

    real, PARAMETER :: PI=3.141592653589793238462643383279502884197
    real, PARAMETER :: PIO2=1.57079632679489661923132169163975144209858
    real, PARAMETER :: TWOPI=6.283185307179586476925286766559005768394
    real, PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967
    real, PARAMETER :: EULER=0.5772156649015328606065120900824024310422
    double precision, PARAMETER :: PI_D=3.141592653589793238462643383279502884197
    double precision, PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858
    double precision, PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394

    TYPE sprs2_sp
        integer :: n,len
        real, DIMENSION(:), POINTER :: val
        integer, DIMENSION(:), POINTER :: irow
        integer, DIMENSION(:), POINTER :: jcol
    END TYPE sprs2_sp
    TYPE sprs2_dp
        integer :: n,len
        double precision, DIMENSION(:), POINTER :: val
        integer, DIMENSION(:), POINTER :: irow
        integer, DIMENSION(:), POINTER :: jcol
    END TYPE sprs2_dp
END MODULE nrtype

!   --------------------------
!
!    This is not all of module NRUTIL, just the parts I need.
!
MODULE nrutil
    USE nrtype
    IMPLICIT NONE
    INTERFACE outerprod
        MODULE PROCEDURE outerprod_r,outerprod_d
    END INTERFACE

     CONTAINS

    FUNCTION outerprod_r(a,b)
    real, DIMENSION(:), INTENT(IN) :: a,b
    real, DIMENSION(size(a),size(b)) :: outerprod_r
    outerprod_r = spread(a,dim=2,ncopies=size(b)) * &
        spread(b,dim=1,ncopies=size(a))
    END FUNCTION outerprod_r

    FUNCTION outerprod_d(a,b)
    double precision, DIMENSION(:), INTENT(IN) :: a,b
    double precision, DIMENSION(size(a),size(b)) :: outerprod_d
    outerprod_d = spread(a,dim=2,ncopies=size(b)) * &
        spread(b,dim=1,ncopies=size(a))
    END FUNCTION outerprod_d

END MODULE nrutil

!   --------------------------

MODULE nr
!
!    Only interface blocks here.  Code is in separate files in
!    various places within Numerical Recipies.  I copied the
!    appropriate code into this file below (search on pythag_dp)
!
    INTERFACE pythag
        FUNCTION pythag_dp(a,b)
        USE nrtype
        double precision, INTENT(IN) :: a,b
        double precision :: pythag_dp
        END FUNCTION pythag_dp

        FUNCTION pythag_sp(a,b)
        USE nrtype
        real, INTENT(IN) :: a,b
        real :: pythag_sp
        END FUNCTION pythag_sp
    END INTERFACE
END MODULE nr


!   --------------------------
module constants
!
!    Constants used in the program
!
    use nrtype
    implicit none

    ! Numerical recipies DP set in module NRTYPE
    !integer, parameter :: dbl = DP              ! this needed for my routines
                                                ! DP and dbl MUST be equal here
                                                !   NOT NEEDED WHEN SWITCHED TO 'double precision'


    double precision :: SVD_ZERO = 0.5E-6  ! This is the value that determines when a singular value is zero.  I.e., when a
                                        ! singular value is less than this, it is considered zero.  This is key when
                                        ! computing rank of the variance-covariance matrix. Mark uses 0.5e-6.


!    integer, parameter :: logfile = 10 ! file handle of the log file, when tracing.
    integer :: logfile = 10 ! used this assignable handle for debugging

!   These are set at run time by calling subroutine set_constants
    double precision :: max_e_able     ! largest number that won't overflow when exp'ed
    double precision :: min_e_able     ! smallest negative number that won't underflow when exp'ed
    double precision :: max_log_able   ! largest number we can take log of
    double precision :: min_log_able   ! smallest positive number we can take log of

    ! These are used when taking numeric second derivatives.  Both relate to computing bit to add to params
    double precision :: eps    = 10e-16       ! My method from SAS manual. I set value of eps by changing
                                            !until standard errors agreed as much as possible with those of Mark.
    double precision :: deltax = 1.0D-8       ! Gary's method.   Got this value from mark.f



    logical, parameter :: central_diffs = .true.  ! whether to take central differences to approximate derivative/gradient.
                              !.False. = one-sided difference which are less accurate, but less costly
                              ! because we don't have to compute log like twice

!    double precision, parameter :: delta = 0.0001   ! Delta muliplier for computing derivatives

    double precision, parameter :: missing = -999999.0D0  ! value to signify missing value to S or R. Used, e.g., for residuals.

!   Parameters needed for GOF tests, but I changed to computing GOF tests in pure R code
!    integer, parameter :: chi_tab_rows = 4 ! # of rows in the chi square tables of mrawin_gof
!    integer, parameter :: orow = 1  ! row number for observed values in fit_table
!    integer, parameter :: erow = 2  ! row number for expected values in fit_table
!    integer, parameter :: oerow = 3 ! row number for Chisq contribution in fit_table
!    integer, parameter :: userow = 4 ! indicator for whether to use the cell in fit_table
!    double precision, parameter :: chi_ruleofthumb = 2.0D0  ! Rule of thumb value for whether to use a cell in a GOF chi square table
!    integer, parameter :: HL_nbins = 10 !number of bins for Hos-Lem statistic

    integer, parameter :: chat_rot = 5 ! C-hat rule of thumb. If all cells in a Test 2 or Test 3 Chi-square table are greater or equal this, use it in computation of c-hat


    double precision, parameter :: pi_mult = 4.0D0   ! Applies to sine link only.  Change this to give sine link more range. 
                                                     ! Sine link is 0 when x less than this number.  Sine link is 1 when x greater than this number.

end module

! ---------------------------------------------------------------------------

module globevars
!
!    This module contains global variables, which include pointers to the huge covariate arrays.
!    This is used so that the calls to the log likelihood procedure can be simple.

    use constants
    implicit none

!    The following pointers to global variables are defined and used in calculating the log likelihood:
!    ptr_nan = pointer to number of animals
!    ptr_nx = pointer to number of capture covariates
!    ptr_ny = pointer to number of survival covariates
!    ptr_np = pointer to nx+ny = total number of parameters
!    ptr_ns = pointer to number of capture occasions
!    ptr_first = pointer to vector containing first capture occasion, nan X 1
!    ptr_last = pointer to vector containing last capture occasion, nan X 1
!    ptr_hist = pointer to capture histories matrix, which is nan X ns
!    ptr_dead = pointer to vector containing death on capture indicators, nan X ns
!    ptr_capX = pointer to capture covariate 3-D array, nan X ns X nx
!    ptr_survX = pointer to survival covariate 3-D array, nan X ns X ny

    integer, pointer :: ptr_nan, ptr_nx, ptr_ny, ptr_ns, ptr_np
    double precision, dimension(:,:,:), pointer :: ptr_capX ! Capture covariates for CJS and Huggins
    double precision, dimension(:,:,:), pointer :: ptr_survX ! Survival covariates for CJS
    double precision, dimension(:,:,:), pointer :: ptr_capY ! RE-capture covariates for Huggins
    integer, dimension(:,:), pointer :: ptr_hist, ptr_dead
    integer, dimension(:), pointer :: ptr_first, ptr_last, ptr_remove
    double precision, dimension(:), pointer :: ptr_intervals

    integer :: trace  ! whether to output intermediate results. Input from user, but stored here.
                      ! values for trace parameter, 0 = no output, 1 = output of end results, 
                      ! 2 = results plus likelihood output, 3 = results plus likelihood plus gradient output
                      ! If trace > 0 and if the VA09AD routine is called, printing occurs every |trace| 
                      ! interations and also on exit.  Output is of the form: Function value, x(1), x(2),...x(n), G(1), ...G(n).
                      ! Intermediate printing from within VA09AD is surpressed if trace > maxit+1.
                      ! Values of X and G are surpressed from within VA09AD if trace < 0 (only final results printed).
                      
    integer :: link   ! The link function to use. The link is specified in R, before calling this. 
                      ! Valid values are:
                      !     1 = logistic
                      !     2 = cosine
                      !     3 = hazard
                      ! Error checking for valid values should be done in R before calling this routine.
                      ! If you ever get around to coding the inverse normal distribution, you could add probit. 
                                            

end module

! ----------------------------------------------------------------------------
subroutine set_constants()
!
!   This routine MUST be called upon entry from R in order t set certain
!   constants that live in module constants.  I would not need this if I could
!   figure out how to set these in module constants. I.e., with an initialization
!   routine in the "contains" section.  But, I can't, so this is the work around.
!
    use constants
    implicit none

    max_e_able = log( huge(max_e_able) - 1.0D0)
    min_e_able = log( tiny(min_e_able) )
    max_log_able = huge( max_log_able ) - 1.0D0
    min_log_able = tiny(min_log_able)


end subroutine set_constants

! ----------------------------------------------------------------------------
subroutine cjsmod( nan, &
        ns, &
        nx, &
        ny, &
        ng, &
        hist, &
        group, &
        algorithm, &
        cov_meth, &
        input_trace, &
        input_link, &
        nhat_v_meth, &
        capX, &
        survX, &
        cap_init, &
        sur_init, &
        max_obs_fn, &
        beta_tol_vec, &
        loglik, &
        deviance, &
        aic, &
        qaic, &
        c_hat, &
        chisq_vif, &
        df_vif, &
        parameters, &
        se_param, &
        covariance, &
        p_hat, &
        se_p_hat, &
        s_hat, &
        se_s_hat, &
        n_hat, &
        se_n_hat, &
        exit_code, &
        pos_def_code, &
        df, &
        intervals )
!
!    Main routine to do CJS Mark-Recapture Analysis.  This is called from
!    R with input parameters.
!
!    Note that previous versions were called MRAWIN.
!
!    One previous version, MRAWIN4, used direct assess files on
!    the hard disk to store covariates and capture indicators in order to
!    cut down on memory requirements.  This version will not work that way.  ]
!    All vectors will be loaded into memory.  This will reduce the size of
!    the largest problem that this program can tackle.
!
!    All this has now been made into an official R package.
!    See the R help files for F.cjs.estim for up-to-date documentations
!    on this routine.  The comments below may be out of date.
!    tlm - 31oct07
!
!    Inputs:
!    nan = number of animals = number of rows in matricies
!    ns = number of samples = number of columns in matricies
!    nx = number of capture covariates, including intercept (i.e., ~1+x would be nx=2)
!    ny = number of survival covariates
!    hist = nan X ns matrix of 0,1,2 indicators for capture
!    capX = nan X (ns*nx) matrix of covariates
!    survX = nan X (ns*ny) matrix of covariates
!    cap_init and sur.init = initial values for parameters
!    algorithm = integer indicator of the optimization algorithm to use (see subroutine estim)
!    cov_meth = integer indicating the method to use to compute covaraince of parameters
!    nhat_v_meth = integer indicating estimator to use to compute v(nhat)
!
!    output:
!    loglik = log likelihood value at the final parameters
!    deviance = deviance
!    aic = aic
!    qaic = qaic
!    c.hat = variance inflation factor
!    chisq_vif = chisquare statistic used to compute the vif
!    input_link = integer value specifying the link function, from the user, "input" to this routine
!    df_vif = degrees of freedom for vif
!    parameters = vector of parameters
!    se_param = vector of standard errors
!    covariance = covariance matrix for parameters estimates
!    p_hat = nan X ns matrix of estimated probability of capture.
!    se_p_hat = standard errors for p_hat
!    s_hat = nan X ns matrix of estiamted survivals
!    se_s_hat = standard errors for s_hat
!    n_hat = ns X 1 vector of population size estimates
!    se_n_hat = nx X 1 vector of standard errors for n_hat
!    exit_code = integer code returned by minimizer
!    pos_def_code = code for positive definiteness of estimated covariance matrix
!   df = number of degrees of freedom = number of real parameters = rank of VC matrix
!
    use constants
    use globevars
    implicit none

!   Can't have this if compiling with g95
!    dll_export mrawin

!    integer, parameter :: real_dbl=selected_real_kind(p=13,r=200)
!    integer, parameter :: dbl=selected_real_kind(p=13,r=200)

!    Input parameters
    integer, intent(inout), target :: nan, ns, nx, ny
    integer, intent(inout) :: ng, algorithm, cov_meth, input_trace, max_obs_fn, input_link
    integer, intent(inout), dimension(nan,ns), target :: hist
    integer, intent(inout), dimension(nan) :: group
    double precision, intent(inout), dimension(nan,ns,nx), target :: capX
    double precision, intent(inout), dimension(nan,ns,ny), target :: survX
    double precision, intent(inout), dimension(nx) :: cap_init
    double precision, intent(inout), dimension(ny) :: sur_init
    double precision, intent(inout), dimension(nx+ny) :: beta_tol_vec
    double precision, intent(inout), dimension(ns), target :: intervals  ! only 1:(ns-1) elements are used

!    Output parameters
    double precision, intent(inout) :: loglik, deviance, aic, qaic, c_hat, &
                                     chisq_vif, df_vif
    double precision, intent(inout), dimension(nx+ny) :: parameters, se_param
    double precision, intent(inout), dimension(nx+ny, nx+ny) :: covariance
    double precision, intent(inout), dimension(nan,ns) :: p_hat, s_hat, se_p_hat, se_s_hat
    double precision, intent(inout), dimension(ns) :: n_hat, se_n_hat
    integer, intent(inout) :: exit_code, pos_def_code, nhat_v_meth, df
!    character(len=*), intent(inout) :: message
!   double precision, intent(inout) :: fit_chisq, fit_chidf
!    double precision, intent(inout), dimension(nan,ns) :: residuals

!   GOF variables, if wanted to compute them here
!   double precision :: t4_chisq, t4_chidf, t5_chisq, t5_chidf, HL_chi, HL_df, roc
!   double precision, dimension(nan,ns) :: residuals
!   double precision, dimension(chi_tab_rows, ns)  :: t4_table
!    double precision, dimension(chi_tab_rows, nan) :: t5_table
!    double precision, dimension(chi_tab_rows, HL_nbins) :: HL_table


!    Local parameters
    integer :: i,j, ioerr
    integer, target :: np
    integer, dimension(nan, ns), target :: idead
    character(len=10) :: label
    integer, dimension(nan), target :: first, last
    character(len=10) :: date
    character(len=12) :: time

!    ---- Set some program constants at run time
    call set_constants()

!    ---- Associate pointers in globevar module
    ptr_nan => nan
    ptr_nx => nx
    ptr_ny => ny
    ptr_ns => ns
    ptr_capX => capX
    ptr_survX => survX
    ptr_hist => hist
    ptr_intervals => intervals

!    ---- Set maximization parameters in globevars
    trace = input_trace
    link = input_link
    
!    ---- Total number of parameters
    np = nx + ny
    ptr_np => np

!    ---- Compute first and last occasions of capture
    call location( nan, ns, hist, first, last)
    ptr_first => first
    ptr_last => last

!    ---- Find the dead ones, change their 2's to 1's, save indicator of deadness
    idead = 0
    do i = 1, nan
            do j= 1, ns
              if (hist(i,j) >= 2 ) then
                 idead(i,j) = 1
                 hist(i,j) = 1
          else if (hist(i,j) < 0) then
             hist(i,j) = 0   ! Just to be sure we only have 0's and 1's in history matrix
              end if
        end do
    end do
    ptr_dead => idead

!    ---- Open a log file, if tracing, primarily for debugging
!    if( trace /= 0 ) then 
!        OPEN(logfile,FILE="mra.log",status="replace",iostat=ioerr)
!        if ( ioerr /= 0 ) then
!            ! Cannot open log file, can't do trace
!            trace = 0
!        end if
!    end if
    
!    if( trace /= 0 ) then
!        !  Write header to log file.
!        call date_and_time( date, time )
!        date = date(1:4) // "-" // date(5:6) // "-" // date(7:8)
!        time = time(1:2) // ":" // time(3:4) // ":" // time(5:10)
!        write(logfile,9000) date, time
!        9000 FORMAT(/" CORMACK-JOLLY-SEBER OPEN POPULATION MODEL."/ &
!                         " Date and time of run: ",a,1x,a/)
!    end if



!    ---- If input c_hat is < 0, do Test 2 and 3.
!        Arrive at a VIF = c_hat for the problem. Output from this routine is
!        c_hat, chisq_vif, and df_vif. These tests are not adjusted for dead animals
!        All variance are adjusted for this c_hat
    if (c_hat <= 0.0D0) then
        call tests(nan, ns, hist, ng, group, c_hat, chisq_vif, i)  ! On return, i is degree of freedom for tests
    else
        call tests(nan, ns, hist, ng, group, aic, chisq_vif, i)  ! aic is just a dummy placeholder, i is df
!        if( trace /= 0 ) write(logfile,"(1x,a,f7.4,a)") "User specified c-hat = ", c_hat, " used."
    end if
    df_vif = i   ! convert to real


!    ---- Do the estimation using supplied initial values
    do i = 1,nx
        parameters(i) = cap_init(i)
    end do
    do i = 1,ny
        parameters(i+nx) = sur_init(i)
    end do

    ! output is parameters, loglik, and covariance and codes. previous values are destroyed
    call CJS_estim(np, algorithm, cov_meth, parameters, loglik, covariance, exit_code, &
        pos_def_code, df, max_obs_fn, beta_tol_vec)

    if( exit_code == 1 ) then

!        ---- Apply the VIF to the covariance matrix
        covariance = covariance * c_hat

!         ---- Check for negative variances, and compute standard errors
        se_param = -1.0D0
        forall( i=1:(nx+ny), covariance(i,i) > 0.0D0 )
            se_param(i) = sqrt( covariance(i,i) )
        end forall


!        ---- Compute model results (aic, dev, etc.)
        deviance=-2.0D0 * loglik
        aic=deviance + 2*df

        if (c_hat > 1.0D0) then
            qaic = deviance/c_hat + 2*df
        else
            qaic = aic
        end if

!        ---- Calculate probability of capture and survival, and SE for both estimates
        call CJS_probs_and_vars(nan,ns,np,parameters,covariance,p_hat,s_hat, se_p_hat,se_s_hat)

!        ---- Calculate estimate of N-hat and variance.
        call est_n_hat(nan,ns,np,covariance,p_hat,se_p_hat,nhat_v_meth, n_hat,se_n_hat)

    else
!        ---- Maximization failed
        se_param = -1.0D0
        deviance = -1.0D0
        aic = -1.0D0
        qaic = -1.0D0
        p_hat = -1.0D0
        s_hat = -1.0D0
        se_p_hat = -1.0D0
        se_s_hat = -1.0D0
    end if


!   ---- Log the results
!    if( trace /= 0 ) then 
!        write(logfile,9030)
!        9030 FORMAT(/" FINAL ESTIMATES OF PARAMETERS WITH STANDARD ERRORS"/)
!        write(logfile,9031) c_hat
!        9031 format(" Standard errors adjusted for c_hat =",f7.4)
!        if( cov_meth == 1 ) then
!            write(logfile,9041)
!            9041 format(" SE = approximation from numerical differentiation"/)
!        else
!            write(logfile,9050)
!            9050 format(" SE = approximation from maximization"/)
!        end if
!    
!        write(logfile,9060)
!        9060 format("     Param   Estimate         SE"/ &
!                        " ==============================="/ )
!        do i = 1, nx+ny
!            write(label,"(i5)") i
!            if( i <= nx ) then
!                label = "Cap" // adjustl(label)
!            else
!                label = "Sur" // adjustl(label)
!            end if
!    
!            WRITE(logfile,9040) label,parameters(i),se_param(i)
!            9040 format (1X,A9,1x,2(F10.6,1x))
!    
!        end do
!    
!        WRITE(logfile,9035) df, nx+ny, loglik, deviance, aic, qaic
!         9035 FORMAT(" ======================================="// &
!                    " Number of parameters (df) =",I8/ &
!                    "    Number of coefficients =",I8/ &
!                    "            Log-likelihood =",F14.6/ &
!                    "                  Deviance =",F14.6/ &
!                    "                       AIC =",F14.6/ &
!                    "                      QAIC =",F14.6)
!    
!    !    ---- Clean up
!        close(logfile)
!    end if

end subroutine cjsmod


! ----------------------------------------------------------------------------

subroutine hugginsmodel( &
        nan, &
        ns, &
        nx, &
        ny, &
        hist, &
        remove, &
        algorithm, &
        cov_meth, &
        nhat_v_meth, &
        capX,&
        capY, &
        p_init, &
        c_init, &
        input_trace, &
        input_link, &
        maxfn, &
        beta_tol_vec, &
        loglik, &
        deviance, &
        aic, &
        parameters, &
        se_param, &
        covariance, &
        p_hat, &
        se_p_hat, &
        c_hat, &
        se_c_hat, &
        n_hat, &
        se_n_hat, &
        n_ci_low, &
        n_ci_high, &
        exit_code, &
        pos_def_code, &
        df )

!
!    Main routine to estimate Huggins closed population Mark-Recapture model.
!    This is the entry point for R.  I.e., R calls this and supplies the input parameters.
!
!    This routine is very similar to 'cjs' above.  I initially copied 'cjs', then
!    modified it.
!
!    Some input parameters:
!        capX = X matrix for initial captures
!        capY = X matrix for recaptures.
!        p_hat = returned probabilities of initial capture
!        c_hat = returned probabilities of recapture
!
!    This is part of the "MRA" package, and as such is documented in
!    the file 'F.huggins.estim.Rd'.
!    tlm - 31oct07
!
    use constants
    use globevars
    implicit none

!       Can't have this if compiling with g95 or gfortran, but use with Lahey
!    dll_export huggins


!    Input parameters
    integer, intent(inout), target :: nan, ns, nx, ny
    integer, intent(inout) :: algorithm, cov_meth, input_trace, maxfn, input_link
    integer, intent(inout), dimension(nan,ns), target :: hist
    integer, intent(inout), dimension(nx), target :: remove   ! = 1 to remove a capture covar from recapture equation
    double precision, intent(inout), dimension(nan,ns,nx), target :: capX  ! initial capture covars
    double precision, intent(inout), dimension(nan,ns,ny), target :: capY  ! recapture covars
    double precision, intent(inout), dimension(nx) :: p_init
    double precision, intent(inout), dimension(ny) :: c_init
    double precision, intent(inout), dimension(nx+ny) :: beta_tol_vec
    

!    Output parameters
    double precision, intent(inout) :: loglik, deviance, aic
    double precision, intent(inout), dimension(nx+ny) :: parameters, se_param
    double precision, intent(inout), dimension(nx+ny, nx+ny) :: covariance
    double precision, intent(inout), dimension(nan,ns) :: p_hat, se_p_hat, c_hat, se_c_hat
    double precision, intent(inout) :: n_hat, se_n_hat, n_ci_low, n_ci_high
    integer, intent(inout) :: exit_code, pos_def_code, nhat_v_meth, df

!    Local parameters
    integer :: i, ioerr
    !integer :: j, k   ! debugging
    integer, target :: np
    character(len=10) :: label
    integer, dimension(nan), target :: first
    character(len=10)  :: date
    character(len=12) :: time

!    ---- Set maximization parameters in globevars
    trace = input_trace
    link = input_link

!    ---- Open a log file, if called for by the user
!    if (trace /= 0) then 
!        OPEN(logfile,FILE="mra.log",status="replace",iostat=ioerr)
!        if ( ioerr /= 0 ) then
!            ! Cannot open log file, so can't do trace
!            trace = 0
!        end if
!    end if

!    if (trace /= 0) then
!        ! ---- Write header to log file.
!        call date_and_time( date, time )
!        date = date(1:4) // "-" // date(5:6) // "-" // date(7:8)
!        time = time(1:2) // ":" // time(3:4) // ":" // time(5:10)
!        write(logfile,9000) date, time
!        9000 FORMAT(/" HUGGINS CLOSED POPULATION MODEL."/ &
!                         " Date and time of run: ",a,1x,a/)
!    end if

    
!    ---- Set some program constants at run time
    call set_constants()

!    ---- Associate pointers in globevar module
    ptr_nan => nan
    ptr_nx => nx
    ptr_ny => ny
    ptr_ns => ns

    ptr_capX => capX  ! initial capture covariates
    ptr_capY => capY  ! recapture covariates
    ptr_hist => hist
    ptr_remove => remove

!    ---- Total number of parameters
    np = nx + ny
    ptr_np => np

!    ---- Compute first and last occasions of capture, but, for Huggins model, we only need first
    call first_capture( nan, ns, hist, first)
    ptr_first => first

!!    ---- Do the estimation using supplied initial values
    do i = 1,nx
        parameters(i) = p_init(i)
    end do
    if (ny >= 1) then
        do i = nx+1, ptr_np
            parameters(i) = c_init(i-nx)
        end do
    end if


   
    call Huggins_estim(ptr_np, algorithm, cov_meth, parameters, loglik, covariance, exit_code, &
        pos_def_code, df, maxfn, beta_tol_vec)

    
    
    
    ! output from estim is parameters, loglik, covariance, df, and codes.
    ! previous values are destroyed
    if( exit_code == 1 ) then

!         ---- Maximization okay

!         ---- Check for negative variances, and compute standard errors
        se_param = -1.0D0
        forall( i=1:ptr_np, covariance(i,i) > 0.0D0 )
            se_param(i) = sqrt( covariance(i,i) )
        end forall

!        ---- Compute model fits (aic, dev, etc.)
        deviance=-2.0D0 * loglik
        aic=deviance + 2.0D0*df

!        ---- Calculate probability of initial and subsequent captures, and SE for both
        call huggins_pc_hat(nan,ns,nx,ny,ptr_np,parameters,covariance,p_hat,se_p_hat,c_hat,se_c_hat)

!        ---- Calculate estimate of N-hat and variance.
        call huggins_n_hat(nan,ns,np,nx,parameters,covariance,p_hat,nhat_v_meth, n_hat,se_n_hat,n_ci_low,n_ci_high)
    else
!        ---- Maximization failed
        se_param = -1.0D0
        deviance = -1.0D0
        aic = -1.0D0
        p_hat = -1.0D0
        se_p_hat = -1.0D0
        c_hat = -1.0D0
        se_p_hat = -1.0D0
        n_hat = -1.0D0
        se_n_hat = -1.0D0
    end if


!   ---- Log the results
!    if ( trace /= 0 ) then 
!        write(logfile,9030)
!        9030 FORMAT(/" FINAL ESTIMATES OF PARAMETERS WITH STANDARD ERRORS"/)
!        if( cov_meth == 1 ) then
!            write(logfile,9041)
!            9041 format(" SE = approximation from numerical differentiation"/)
!        else
!            write(logfile,9050)
!            9050 format(" SE = approximation from maximization"/)
!        end if
!    
!        write(logfile,9060)
!        9060 format("     Param   Estimate         SE"/ &
!                        " ==============================="/ )
!        do i = 1, nx
!            write(label,"(i5)") i
!            label = "  Cap" // adjustl(label)
!    
!            WRITE(logfile,9040) label,parameters(i),se_param(i)
!            9040 format (1X,A9,1x,2(F10.6,1x))
!        end do
!        if (ny >= 1) then
!            do i = nx+1, nx+ny
!
!         write(label,"(i5)") (i - nx)
!                label = "Recap" // adjustl(label)
!    
!                WRITE(logfile,9040) label,parameters(i),se_param(i)
!            end do
!        end if
!    
!        WRITE(logfile,9035) df, nx+ny, loglik, deviance, aic
!         9035 FORMAT(" ======================================="// &
!                    " Number of parameters (df) =",I8/ &
!                    "    Number of coefficients =",I8/ &
!                    "            Log-likelihood =",F14.6/ &
!                    "                  Deviance =",F14.6/ &
!                    "                       AIC =",F14.6)
!    
!        WRITE(logfile,9100) n_hat, se_n_hat, n_ci_low, n_ci_high
!         9100 FORMAT(" ======================================="// &
!                    "  Population size estimate =",F14.6/ &
!                    "       SE(Population size) =",F14.6/ &
!                    "     Lower 95% CI endpoint =",F14.6/ &
!                    "     Upper 95% CI endpoint =",F14.6)
!    
!    !    ---- Clean up
!        close(logfile)
!        
!    end if    

end subroutine hugginsmodel


! ---------------------------------------------------------------------------------------------

subroutine CJS_estim(np, algorithm, cov_meth, parameters, loglik, covariance, exit_code, cov_npd, &
    df, max_fn, beta_tol_vec )
!
!    Purporse: to maximize the capture-recapture likelihood
!
!    Inputs:
!    np = number of parameters = nx+ny
!    parameters = (np) vector of initial values
!    algorithm = integer designating the algorithm to use
!    cov_meth = method for computing covariance matrix. 1=2nd derivatives, 2=Hessian.
!
!    Output:
!    parameters = (np) vector of maximized parameters
!    loglik = maximized log likelihood value
!    covaraince = (np) x (np) matrix of covariances
!    exit_code = integer exit code from the minimization routine, see comments there for codes
!    cov_npd = 0 if covariance is positive definite, 1 if not.
!    mess = text message for exit_code  REMOVED BECAUSE R .FORTRAN HANDLES CHARACTERS DIFFERENT THAN S-PLUS .FORTRAN
!
!    A bunch of pointers to global parameters are used to calculate the likelihood, but are
!    not passed into here.

    use constants
    use globevars
    implicit none

    integer, intent(inout) :: np, algorithm, exit_code, cov_npd, cov_meth, df
    integer, intent(inout) :: max_fn  ! upon output, max_fn is the actual number of function evaluations

    double precision, intent(inout) :: loglik
    double precision, intent(inout), dimension(np) :: parameters
    double precision, intent(inout), dimension(np) :: beta_tol_vec ! required precision in each element of parameters
    double precision, intent(inout), dimension(np, np) :: covariance
    

    ! Local variables
    double precision, dimension(np) :: g ! gradient vector
    double precision, dimension(np*(np+1)/2) :: h ! hessian, see VA09AD for description of format
    double precision, dimension(3*np) :: W
    double precision, external :: CJS_loglik
    integer, external :: matrank
    external CJS_obj  ! objective function for minimization
    integer :: ij, i, j

    ! Local copies of some parameters. Added to get around the following warning: 
    ! Warning: Non-variable expression in variable definition context (actual argument to INTENT = OUT/INOUT)    
    !   You can't associate a number, like '1', with an intent(inout) variable in function calls
    double precision :: dfn = -2.0
    integer :: mode = 1


    ! Do the minimization using specified algorithm
    select case (algorithm)

        case (:1, 3:)

            ! The algorithm used by Mark
            call VA09AD(CJS_obj,np,parameters,loglik,g,h,W, dfn, beta_tol_vec, mode, max_fn, trace, exit_code)

        case (2)

            ! The algorithm used by MRAWIN4 = DFPMIN
            exit_code = -1

    end select

    if( .not.(exit_code == 1)) then
        ! Maximization did not work for some reason, do not run any of the hessian or covariance routines
!        if ( trace  /= 0) then
!            write(logfile,*)
!            write(logfile,*) " *** Likelihood maximization failed: Exit_code = ", exit_code
!            if (exit_code == 0) then
!                write(logfile,*) "     Hessian not positive definite."
!            else if (exit_code == 2) then
!                write(logfile,*) "     Underflow/rounding error."
!            else if (exit_code == 3) then
!                write(logfile,*) "     Did not converge - too many function evaluations."
!            end if
!        end if
        parameters = 0
        covariance = -1
        cov_npd = 1
        df = 0
    else
        ! Maximization did converge, compute covariances

        ! When return from minimization algorithm, logliklihood is positive (so can minimize).
        ! Change it back to negative.
        loglik = -loglik


        ! ---- Hessian and Covariance Section
        if( cov_meth == 1 ) then

!            if( trace /= 0) write(logfile,*) "Calling comp_hessian..."
            ! compute hessian matrix by numeric 2nd derivatives.  Covariance is actually the
            ! Hessian.  I'm just using covarinace as storage here
            call comp_hessian(CJS_loglik, np, parameters, loglik, covariance)
            cov_npd = 0

        else
            ! Use the returned Hessian
            ! Routine MC11FD inverts h and converts from factorized form to regular (triangular) form, all at once
            ij = np
            call MC11FD(h,np,ij)
            cov_npd = 0

            ! Put into a square
            IJ = 0
            DO j = 1,np
                IJ = ij + 1
                covariance(j,j) = h(ij)
                DO i = j+1,np
                    IJ = ij + 1
                    covariance(i,j) = h(IJ)
                    covariance(j,i) = h(ij)
                end do
            end do

            ! I realize this is inefficient, but invert the covariance back to hessian
            ! in order to compute number of parameters.
            cov_npd = 0
            call syminv(covariance, np, cov_npd)
        end if

        ! ------ Compute number of parameters = rank of hessian matrix.
        !   Gary White does this, computes rank, on the hessian prior to inverting.
        !   I do not fully understand this, but the singular values are much smaller and better
        !   behaved when computed on hessian than when computed on inverted hessian. So,
        !   we compute singular value decomposition on hessian, then count number of
        !   non-zero SV's.  Later, we invert the hessian and check for non-positive definitness.
        !
        !       cov_npd > 0 if inversion for cov.meth = 2 fails
        !       df = 1 if we want to estimate number of parameters (the norm)
        !       df = 0 means don't bother, R sets df to user specified value or np

        if( cov_npd == 0 .and. df > 0 ) then
            !df = np
            df = matrank(covariance, np, np)
!            if( trace /= 0 ) write(logfile,*) "Number of estimated parameters = ", df
        else
!            if( trace /= 0 ) write(logfile,*) "Number of estimated parameters =  USER OVERRIDE"
        end if


        ! ---- Now invert the negative matrix of 2nd derivatives
        call syminv(covariance, np, cov_npd)

!        if( trace /= 0 ) then 
!            write(logfile,*) " ----- Coefficient covariance matrix -----"
!            do i = 1, np
!                write(logfile,"(1000(g20.10,','))") (covariance(i,j), j=1,np)
!            end do
!        end if

        !   Could check for positive definiteness here


        if( cov_npd == 1 ) then
!            if( trace /= 0 ) write(logfile,*) "COVARIANCE MATRIX IS SINGULAR.  Error code= ", cov_npd
            !df = 0
        else if (cov_npd == 2) then
!            if( trace /= 0 ) write(logfile,*) "COVARIANCE IS NOT POSITIVE DEFINITE. Error code= ", cov_npd
        end if
    end if


end subroutine

! ---------------------------------------------------------------------------------------------

subroutine Huggins_estim(np, algorithm, cov_meth, parameters, loglik, covariance, exit_code, cov_npd, df, &
    maxfn, beta_tol_vec )
!
!    Purporse: to maximize the Huggins closed popln capture-recapture likelihood
!
!    Inputs:
!    np = number of parameters = nx
!    parameters = (np) vector of initial values
!    algorithm = integer designating the algorithm to use
!    cov_meth = method for computing covariance matrix. 1=2nd derivatives, 2=Hessian.
!
!    Output:
!    parameters = (np) vector of maximized parameters
!    loglik = maximized log likelihood value
!    covaraince = (np) x (np) matrix of covariances
!    exit_code = integer exit code from the minimization routine, see comments there for codes
!    cov_npd = 0 if covariance is positive definite, 1 if not.
!
!    A bunch of pointers to global parameters are used to calculate the likelihood, but are
!    not passed into here.
!
!    Note: this is largely redundance code because 'CJS_estim' (for CJS) is almost identical.
!    You could write a wrapper function to these to, and pass in the objective function
!    and likelihood function. But, for now repeat the code...

    use constants
    use globevars
    implicit none

    integer, intent(inout) :: np, algorithm, exit_code, cov_npd, cov_meth, df, maxfn

    double precision, intent(inout) :: loglik
    double precision, intent(inout), dimension(np) :: parameters
    double precision, intent(inout), dimension(np, np) :: covariance
    double precision, dimension(np) :: beta_tol_vec ! required precision in each element of parameters

    ! Local variables
    double precision, dimension(np) :: g ! gradient vector
    double precision, dimension(np*(np+1)/2) :: h ! hessian, see VA09AD for description of format
    double precision, dimension(3*np) :: W

    double precision, external :: Huggins_loglik  ! Needed for covariance estimation
    integer, external :: matrank

    external Huggins_obj  ! objective function for minimization

    integer :: ij, i, j

    double precision :: dfn = -2.0
    integer :: mode = 1


    ! Set required precision
    !beta_tol_vec = beta_tol

    ! Do the minimization using specified algorithm
    select case (algorithm)

        case (:1, 3:)

            ! The algorithm used by Mark
            call VA09AD(Huggins_obj,np,parameters,loglik,g,h,W, dfn, beta_tol_vec, mode, maxfn, trace, exit_code)

        case (2)

            ! The algorithm used by MRAWIN4 = DFPMIN
            exit_code = -1

    end select

    if( .not.(exit_code == 1)) then
        ! Maximization did not work for some reason, do not run any of the hessian or covariance routines
!        if( trace /= 0 ) then
!            write(logfile,*)
!            write(logfile,*) " *** Likelihood maximization failed: Exit_code = ", exit_code
!            if (exit_code == 0) then
!                write(logfile,*) "     Hessian not positive definite."
!            else if (exit_code == 2) then
!                write(logfile,*) "     Underflow/rounding error."
!            else if (exit_code == 3) then
!                write(logfile,*) "     Did not converge - too many function evaluations."
!            end if
!        end if 
        parameters = 0
        covariance = -1
        cov_npd = 1
        df = 0
    else
        ! Maximization did converge, compute covariances

        ! When return from minimization algorithm, logliklihood is positive (so can minimize).
        ! Change it back to negative.
        loglik = -loglik


        ! ---- Hessian and Covariance Section
        if( cov_meth == 1 ) then

            ! compute hessian matrix by numeric 2nd derivatives.  What is stored in covariance is actually the
            ! Hessian.  I'm just using covarinace as storage here
            call comp_hessian(Huggins_loglik, np, parameters, loglik, covariance)
            cov_npd = 0

        else
            ! Use the returned Hessian
            ! Routine MC11FD inverts h and converts from factorized form to regular (triangular) form, all at once
            ij = np
            call MC11FD(h,np,ij)
            cov_npd = 0

            ! Put into a square
            IJ = 0
            DO j = 1,np
                IJ = ij + 1
                covariance(j,j) = h(ij)
                DO i = j+1,np
                    IJ = ij + 1
                    covariance(i,j) = h(IJ)
                    covariance(j,i) = h(ij)
                end do
            end do

            ! I realize this is inefficient, but invert the covariance back to hessian
            ! in order to compute number of parameters.
            cov_npd = 0
            call syminv(covariance, np, cov_npd)
        end if

        ! ------ Compute number of parameters = rank of hessian matrix.
        !   Gary White does this, computes rank, on the hessian prior to inverting.
        !   I do not fully understand this, but the singular values are much smaller and better
        !   behaved when computed on hessian than when computed on inverted hessian. So,
        !   we compute singular value decomposition on hessian, then count number of
        !   non-zero SV's.  Later, we invert the hessian and check for non-positive definitness.
        !
        !       cov_npd > 0 if inversion for cov.meth = 2 fails
        !       df = 1 if we want to estimate number of parameters (the norm)
        !       df = 0 means don't bother, R sets df to user specified value or np

        if( cov_npd == 0 .and. df > 0 ) then
            df = matrank(covariance, np, np)
!            if( trace /= 0 ) write(logfile,*) "Number of estimated parameters = ", df
        else
!            if( trace /= 0 ) write(logfile,*) "Number of estimated parameters =  USER OVERRIDE"
        end if


        ! ---- Now invert the negative matrix of 2nd derivatives
        call syminv(covariance, np, cov_npd)

!        if( trace /= 0 ) then
!            write(logfile,*)
!            write(logfile,*) " ----- Coefficient covariance matrix -----"
!            do i = 1, np
!                write(logfile,"(1000(g20.10,','))") (covariance(i,j), j=1,np)
!            end do
!        end if
        
        !   Could check for positive definiteness here


        if( cov_npd == 1 ) then
!            if( trace /= 0 ) write(logfile,*) "COVARIANCE MATRIX IS SINGULAR.  Error code= ", cov_npd
            df = 0
        else if (cov_npd == 2) then
!            if( trace /= 0 ) write(logfile,*) "COVARIANCE IS NOT POSITIVE DEFINITE. Error code= ", cov_npd
        end if
    end if


end subroutine

! ----------------------------------------------------------------------------------

SUBROUTINE CJS_obj(p,beta,lnlik,grad)
!
!    Purpose: to calculate the Capture-recapture loglikelihood given data and values of
!    parameters.  This is the routine that will be minimized.
!
!    Inputs:
!    p = number of parameters/coefficients in beta
!    beta = 1 x p vector of current coefficients
!
!    Output:
!    lnlik = the negative of the CJS log-likelihood (since we are minimizing)
!    grad = 1 x p vector containing the gradient.  I.e., partial derivatives of the
!        likelihood for parameters at thier current location.
!
!
    use constants
    implicit none

    integer, intent(inout) :: p
    double precision, intent(inout) :: lnlik
    double precision, intent(inout), dimension(p) :: grad, beta

    double precision, external :: CJS_loglik

    !integer :: i   ! debugging var


!    Calculate the log-likelihood
    lnlik = -1.0D0 * CJS_loglik(p, beta)

!    Calculate the gradient
    call CJS_gradient(p, beta, lnlik, grad)

end subroutine



! ------------------------------------------------------------------------------------
double precision function CJS_loglik(p, beta)
!double precision function CJS_loglik(p, beta)  THIS DOES NOT WORK IN GFORTRAN, but worked in G95 and Lahey
!
!         Purpose: to compute the log-likelihood, given parameters
!
!    Input:
!    p = number of parameters
!    beta = coefficients in the logistic equation
!
!    Global variables, accessed through pointers:
!    ptr_nan = pointer to number of animals
!    ptr_nx = pointer to number of capture covariates
!    ptr_ny = pointer to number of survival covariates
!    ptr_ns = pointer to number of capture occasions
!    ptr_first = pointer to vector containing first capture occasion, nan X 1
!    ptr_last = pointer to vector containing last capture occasion, nan X 1
!    ptr_hist = pointer to capture histories matrix, which is nan X ns
!    ptr_dead = pointer to vector containing death on capture indicators, nan X 1
!    ptr_capX = pointer to capture covariate 3-D array, nan X ns X nx
!    ptr_survX = pointer to survival covariate 3-D array, nan X ns X ny
!
!    All of these global variable pointers are stored in module globevars and set in subroutine
!    mrawin.

    use constants
    use globevars
    implicit none

    !double precision :: CJS_loglik
    integer, intent(inout) :: p
    double precision, intent(inout), dimension(p) :: beta


    double precision, dimension(ptr_nx) :: cap_beta
    double precision, dimension(ptr_ny) :: surv_beta
    double precision :: pij, phiij, sum1, sum2, prod, xlnlik
    double precision, dimension(ptr_ns) :: vpij, vphiij
    !double precision, dimension(ptr_ns,p) :: W
    integer :: i, j, init1, init2, ii, jj

    cap_beta = beta(1:ptr_nx)
    surv_beta = beta( (ptr_nx+1):p )

    xlnlik=0.0D0
    do i=1,ptr_nan

        ! ---- Debugging: ptr_hist(i,j) should never be 2, only 0 and 1
        !do j=1,ptr_ns
        !    if( ptr_hist(i,j) == 0 .or. ptr_hist(i,j) == 1 ) cycle
        !    write(logfile,*) "**** Invalid capture history with a 2 ***: Animal", i, 
        !       " Occasion", j, " Hist value=", ptr_hist(i,j)
        !end do

        ! ---- First part of the log-likelihood function, between first and last 1
        sum1=0.0D0
        sum2=0.0D0
        vpij = 0.0D0
        vphiij = 0.0D0


        ! First(i) is never 1, because subroutine location sets first to 1 occasion
        ! after the first encounter = first estimable capture probability.
        ! Note first(i) = 0 for histories with first capture at the last occasion
        if (ptr_first(i) == 0) then
            init1=ptr_ns+1
            init2=ptr_ns+1
        else if (ptr_first(i) > 0) then
            init1=ptr_first(i)
            init2=ptr_first(i)-1
        end if


        ! Compute all probabilities of capturing bird i, from first occasion to end.
        if (init1 <= ptr_ns) then
            do j=init1,ptr_ns
                call procap(pij, i, j, cap_beta, ptr_nx)
                vpij(j)=pij
            end do
        end if


        ! compute probability of bird i surviving to time j
        if (init2 < ptr_ns) then
            do j=init2,ptr_ns-1
                call prosur(phiij, i, j, surv_beta, ptr_ny)
                vphiij(j)=phiij
            end do
        end if


        ! Compute log-likelihood contribution for animal i
        if ((ptr_first(i) > 0) .and. (ptr_first(i) <= ptr_last(i))) then
            do j=ptr_first(i),ptr_last(i)
                sum1=sum1 + ptr_hist(i,j)*log(vpij(j)) + &
                    (1-ptr_hist(i,j))*log(1.0D0-vpij(j)) + &
                    log(vphiij(j-1))
            end do
        end if


        ! ---- Find second part of likelihood, after last 1 = Chi parameters
        ! Chi = probability of animal i not being seen again
        ! If animal died on capture before release, prob of not seeing again is 1
        if (ptr_dead(i,ptr_last(i)) == 1) then
            sum2 = 0.0D0
        else if ((ptr_last(i) > 0) .and. (ptr_last(i) < ptr_ns)) then
            sum2=1.0D0-vphiij(ptr_last(i))
            do ii=ptr_last(i)+1, ptr_ns
                prod=1.0D0
                do jj=ptr_last(i),ii-1
                    prod=prod*vphiij(jj)*(1.0D0-vpij(jj+1))
                end do
                if (ii < ptr_ns) then
                    prod=prod*(1.0D0 - vphiij(ii))
                endif
                sum2=sum2+prod
            end do
            sum2=log(sum2)
        endif


        ! inserted debugging code. Set PRNT in constants module
!        if( trace >= 2 ) then
!            write(logfile,'(1x,a,5(1x,a,i5))') "---------", "i=", i, "init1=", init1, "init2=", init2, "first(i)=", &
!                ptr_first(i), "last(i)=", ptr_last(i)
!            write(logfile,*) "cap_coefs=", (cap_beta(j), j=1,ptr_nx)
!            write(logfile,*) "surv_coefs=", (surv_beta(j), j=1,ptr_ny)
!            write(logfile,62) "vpij=", (vpij(j), j=2,ptr_ns)
!            write(logfile,62) "vphiij=", (vphiij(j), j=1,ptr_ns-1)
!            write(logfile,63) "sum1=", sum1, " sum2=", sum2
!            !61  format(1x,a,25(1x,f4.2))
!            62  format(1x,a,25(1x,f6.3))
!            63  format(1x,a, f12.7,a, f12.7)
!        end if

        xlnlik=xlnlik+sum1+sum2
    end do

!    if (trace >= 2) then
!        write(logfile,*)
!        write(logfile,*) "Log likelihood = ", xlnlik
!    end if

    CJS_loglik = xlnlik

end function

! ------------------------------------------------------------------------------------

SUBROUTINE CJS_gradient(p, beta, f, grad)
!
!    Purpose: compute the derivative of the likelihood at beta.
!
!    Input:
!    p = number of coefficients
!    beta = coefficients
!    f = value of loglikelihood at beta
!
!    Output:
!    grad = pX1 vector of partial derivatives of log likelihood w.r.t. all parameters.
!        Derivatives computed by either central differences or one-sided differences.
!        Method is controled by the constant central_diffs.  Central differences are
!        more accurate, but you have to compute the likelihood twice for each parameter,
!        instead of just once.
!
    use constants
    use globevars
    implicit none

    integer, intent(inout) :: p
    double precision, intent(inout), dimension(p) :: beta
    double precision, intent(inout), dimension(p) :: grad
    double precision, intent(inout) :: f

    double precision, external :: CJS_loglik
    double precision :: f1, f2, deltai, tmp_b
    double precision, dimension(p) :: beta2
    integer :: i

    beta2 = beta

    ! central_diffs and delta set in module constants

    do i=1,p
        deltai=(deltax/2.0D0)*(1.0D0 + abs(beta2(i)))*1.D5

        tmp_b = beta2(i)
        beta2(i) = beta2(i) + deltai
        f1 = -1.0D0 * CJS_loglik(p, beta2)

        if (central_diffs) then
            beta2(i) = tmp_b - deltai
            f2 = -1.0D0 * CJS_loglik(p, beta2)
            grad(i)=(f1-f2)/(2.0D0 * deltai)
        else
            grad(i)=(f1-f)/deltai
        end if

        beta2(i) = tmp_b

    end do

!    if(trace >= 3) then
!        write(logfile,*) "Gradient vector:"
!        write(logfile,*) "    i         coef     Gradient"
!        write(logfile,*) "----- ------------ ------------"
!        write(logfile,10) (i, beta2(i), grad(i), i=1,p)
!        10 format(1x,i5,1x,2F12.7)
!    end if


end subroutine



! ----------------------------------------------------------------------------------

SUBROUTINE Huggins_obj(p,beta,lnlik,grad)
!
!    Purpose: to calculate the Huggins loglikelihood given data and values of
!    parameters.  This is the routine that will be minimized.
!
!    Inputs:
!    p = number of parameters/coefficients in beta
!    beta = 1 x p vector of current coefficients
!
!    Output:
!    lnlik = the negative of the CJS log-likelihood (since we are minimizing)
!    grad = 1 x p vector containing the gradient.  I.e., partial derivatives of the
!        likelihood for parameters at thier current location.
!
!
    use globevars
    implicit none

    integer, intent(inout) :: p
    double precision, intent(inout) :: lnlik
    double precision, intent(inout), dimension(p) :: grad, beta

    double precision, external :: Huggins_loglik

    !debugging
    integer :: i

!    Calculate the log-likelihood
!    if( trace >= 2 ) then 
!        write(logfile,*) "*****CALCULATING LOG LIKE******"
!        write(logfile,*) "p=", p, "beta=", (beta(i), i=1,p)
!
!        close(logfile)
!        open(logfile,file='mra.log',status='old',access='append')
!        
!    end if
    
    lnlik = -1.0D0 * Huggins_loglik(p, beta)

!    Calculate the gradient
!    if( trace >= 3 ) then
!        write(logfile,*) "*****CALLING GRADIENT******"
!        close(logfile)
!        open(logfile,file='mra.log',status='old',access='append')
!    end if
    
    call Huggins_gradient(p, beta, lnlik, grad)

end subroutine



! ------------------------------------------------------------------------------------
double precision function Huggins_loglik(p, beta)
!
!         Purpose: to compute the Huggins model log-likelihood, given parameters
!
!    Input:
!    p = number of parameters
!    beta = coefficients in the logistic equation
!
!    Global variables, accessed through pointers:
!    ptr_nan = pointer to number of animals
!    ptr_ny = pointer to number of recapture covariates
!    ptr_nx = pointer to number of initial capture covariates
!    ptr_ns = pointer to number of capture occasions
!    ptr_first = pointer to vector containing first capture occasion, nan X 1
!    ptr_hist = pointer to capture histories matrix, which is nan X ns
!    ptr_capX = pointer to capture covariate 3-D array, nan X ns X nx
!
!    All of these global variable pointers are stored in module globevars and set in subroutine
!    mrawin.

    use constants
    use globevars
    implicit none

    !double precision :: Huggins_loglik
    integer, intent(inout) :: p    ! number of parameters
    double precision, intent(inout), dimension(p) :: beta    ! parameters


    ! Note: 'p' is for initial capture probability.
    !       'c' is for recapture probability

    !double precision, dimension(ptr_nx) :: p_beta  ! Coefficients in model for initial captures
    !double precision, dimension(ptr_ny) :: c_beta  ! Coefficients in model for subsequent recaptures
    double precision :: sum_ppart, sum_cpart, xlnlik, denom
    double precision, dimension(ptr_ns) :: vpij, vcij
    integer :: i, j

    !p_beta = beta(1:ptr_nx)
    !c_beta = beta( (ptr_nx+1):p )

    xlnlik=0.0D0
    do i=1,ptr_nan

        ! ---- Debugging: ptr_hist(i,j) should never be 2, only 0 and 1
        !do j=1,ptr_ns
        !    if( ptr_hist(i,j) == 0 .or. ptr_hist(i,j) == 1 ) cycle
        !    write(logfile,*) "**** Invalid capture history with a 2 ***: Animal", i, " Occasion", j, " Hist value=", ptr_hist(i,j)
        !end do
        !if (i==1) write(logfile,*) "animal=", i, "----------------------------------------"

        ! ---- First part of the log-likelihood function, between initial occasion and first capture
        sum_ppart=0.0D0
        sum_cpart=0.0D0
        vpij = 0.0D0
        vcij = 0.0D0
        denom = 1.0D0

        ! Compute probabilities of initial capture for animal i.
        ! We need capture probability from first occasion to last occasion for denominator.
        ! For numerator, we only need probability of capture from first occasion to first capture.
        ! Recall: ptr_first(i) is occasion of first capture
        !   Covariates for capture model are in capX or ptr_capX
        !   If ptr_ny > 0, covariates for recapture model are in capY or ptr_capY
        do j=1,ptr_ns
            call procap(vpij(j), i, j, beta, ptr_nx)  ! beta is at least ptr_nx long, so this works
            denom = denom * (1-vpij(j))
            if( j <= ptr_first(i) ) then
                sum_ppart=sum_ppart + ptr_hist(i,j)*log(vpij(j)) + &
                    (1-ptr_hist(i,j))*log(1.0D0-vpij(j))
            end if
        end do
        
        denom = log(1-denom)

        ! compute probability of recaptures after initial capture.
        ! If ptr_ny == 0, this is same as procap(...)
        if (ptr_first(i) < ptr_ns) then
            do j=(ptr_first(i)+1),ptr_ns
                call prorecap(vcij(j), i, j, beta, ptr_nx, ptr_ny, ptr_remove )
                sum_cpart=sum_cpart + ptr_hist(i,j)*log(vcij(j)) + &
                    (1-ptr_hist(i,j))*log(1.0D0-vcij(j))
            end do
        end if

        ! Compute log-likelihood contribution for animal i
        !write(logfile, *) "Done with animal ", i
        !write(logfile, *) "sum_part = ", sum_ppart, "sum_cpart = ", sum_cpart, "denom =", denom
        !write(logfile, *) "---------------------------------------------------------------------------------------"
        
        xlnlik = xlnlik + sum_ppart + sum_cpart - denom

    end do

!    if (trace >= 2) then
!        write(logfile,*)
!        write(logfile,*) "Log likelihood = ", xlnlik
!    end if

    Huggins_loglik = xlnlik

end function

! ------------------------------------------------------------------------------------

SUBROUTINE Huggins_gradient(p, beta, f, grad)
!
!    Purpose: compute the derivative of the likelihood at beta.
!
!    Input:
!    p = number of coefficients
!    beta = coefficients
!    f = value of loglikelihood at beta
!
!    Output:
!    grad = pX1 vector of partial derivatives of log likelihood w.r.t. all parameters.
!        Derivatives computed by either central differences or one-sided differences.
!        Method is controled by the constant central_diffs.  Central differences are
!        more accurate, but you have to compute the likelihood twice for each parameter,
!        instead of just once.
!
    use constants
    use globevars
    implicit none

    integer, intent(inout) :: p
    double precision, intent(inout), dimension(p) :: beta
    double precision, intent(inout), dimension(p) :: grad
    double precision, intent(inout) :: f

    double precision, external :: Huggins_loglik
    double precision :: f1, f2, deltai, tmp_b, tmp_g
    double precision, dimension(p) :: beta2
    integer :: i

    beta2 = beta

    ! central_diffs and delta set in module constants

    do i=1,p
        deltai=(deltax/2.0D0)*(1.0D0 + abs(beta2(i)))*1.D5

        tmp_b = beta2(i)
        beta2(i) = beta2(i) + deltai
        f1 = -1.0D0 * Huggins_loglik(p, beta2)

        if (central_diffs) then
            beta2(i) = tmp_b - deltai
            f2 = -1.0D0 * Huggins_loglik(p, beta2)
                    grad(i)=(f1-f2)/(2.0D0 * deltai)
            tmp_g = (f1-f)/deltai  ! for debug printing only
        else
                    grad(i)=(f1-f)/deltai
            tmp_g = 0.0D0
        end if

        beta2(i) = tmp_b

           end do

!    if(trace >= 3) then
!        write(logfile,*) "Gradient vector:"
!        write(logfile,*) "    i         coef     Gradient    1-sided G"
!        write(logfile,*) "----- ------------ ------------ ------------"
!        write(logfile,10) (i, beta2(i), grad(i), i=1,p)
!        10 format(1x,i5,1x,3F12.7)
!    end if


end subroutine



! --------------------------------------------------------------------------

subroutine procap(pij, i, j, coef, nx)
!
!     Subroutine to evaluate the probability of capture for bird i
!
!    Input:
!    i = animial number
!    j = occasion number
!    coef = 1 x nx vector of coefficients
!    nx = number of coefficients
!
!    Output:
!    pij = probability of capture for animal i at occasion j
!
!    Uses the ptr_capX array in module globevars

    use constants
    use globevars
    implicit none

    integer, intent(in) :: nx, i, j
    double precision, intent(in), dimension(nx) :: coef
    double precision, intent(out) :: pij

    double precision, external :: logit_link, sine_link, hazard_link
    double precision :: sum
    integer :: k

    sum = 0.0D0


    do k = 1, nx
        sum = sum + coef(k)*ptr_capX(i,j,k)
    end do

    !   Apply the link function.  Correct link to use was specified by user and stored in 'link' in globevars
    if( link == 1 ) then
        pij = logit_link( sum )
        
    else if (link == 2) then
        pij = sine_link( sum )
        
    else if (link == 3) then
        pij = hazard_link( sum )
        
    else  ! unknown link function.  This will not bomb gracefully.
        pij = -1.
    end if
    

!    if( i <= 2 ) then
!        write(logfile,*) "  Cap", i, j, ":", (coef(k), "*", ptr_capX(i,j,k), " + ", k=1,nx), "p=", pij
!    end if

end subroutine

! -------------------------------------------------------------------

subroutine prorecap(cij, i, j, coef, nx, ny, remove)
!
!     Subroutine to evaluate the probability of RE-capture for animal i
!
!    Input:
!    i = animial number
!    j = occasion number
!    coef = 1 x (nx+ny) vector of coefficients
!    nx = number of coefficients in capture equation
!    ny = number of coefficients in recapture equation, could be 0.
!    remove = vector of 0's and 1's indicating which capture coefficients to remove from recapture equation
!
!    Output:
!    pij = probability of capture for animal i at occasion j
!
!    Note: this uses ptr_capY array in module globevars for covariates, rather than capX.
!    Otherwise, this and procap routine are the same.

    use constants
    use globevars
    implicit none

    integer, intent(in) :: nx, ny, i, j
    double precision, intent(in), dimension(nx+ny) :: coef
    double precision, intent(out) :: cij
    integer, intent(in), dimension(nx) :: remove

    double precision, external :: logit_link, sine_link, hazard_link
    double precision :: sum, z
    integer :: k

    sum = 0.0D0



    do k = 1, nx+ny
        if (k <= nx) then
            if (remove(k) == 0) then
                sum = sum + coef(k)*ptr_capX(i,j,k)
            end if
        else
            sum = sum + coef(k)*ptr_capY(i,j,k-nx)
        end if
    end do

    !   Apply the link function.  Correct link to use was specified by user and stored in 'link' in globevars
    if( link == 1 ) then
        cij = logit_link( sum )
        
    else if (link == 2) then
        cij = sine_link( sum )
        
    else if (link == 3) then
        cij = hazard_link( sum )
        
    else  ! unknown link function.  This will not bomb gracefully.
        cij = -1.
    end if


!    if (i==1) then
!        write(logfile,*) "Recap", i, j, ":", (ptr_capY(i,j,k), k=1,ny), "p=", cij
!        write(logfile,*) "     nx=", nx, "ny=", ny, "coef=", (coef(k),k=1,nx+ny)
!    end if

end subroutine

! -------------------------------------------------------------------

subroutine prosur(sij, i, j, coef, ny)
!
!     Subroutine to evaluate the probability of survival for bird i
!
!    Input:
!    i = animial number
!    j = occasion number
!    coef = 1 x ny vector of coefficients
!    ny = number of coefficients
!
!    Output:
!    sij = probability of survival for animal i between occasion j and j+1
!
!    Uses the ptr_survX array in module globevars

    use constants
    use globevars
    implicit none

    integer, intent(in) :: ny, i, j
    double precision, intent(in), dimension(ny) :: coef
    double precision, intent(out) :: sij

    double precision, external :: logit_link, sine_link, hazard_link
    double precision :: sum, z
    integer :: k

    sum = 0.0D0

    do k = 1, ny
        sum = sum + coef(k)*ptr_survX(i,j,k)
    end do

    !   Apply the link function.  Correct link to use was specified by user and stored in 'link' in globevars
    if( link == 1 ) then
        sij = logit_link( sum )
        
    else if (link == 2) then
        sij = sine_link( sum )
        
    else if (link == 3) then
        sij = hazard_link( sum )
        
    else  ! unknown link function.  This will not bomb gracefully.
        sij = -1.
    end if

    ! Account for the time interval here.  This is only place in whole program where
    ! we use the interval information.
    sij = sij ** ptr_intervals(j)

end subroutine

!-------------------------------------------

double precision function logit_link( eta )
!
!   Compute the inverse of the logistic link. 
!   eta is the linear predictor
!
    use constants
    implicit none
    
    !double precision :: logit_link
    double precision :: eta
    double precision :: z

    ! Careful, these two statements set boundaries where gradients fail (not differentiable)
    ! With double precision, we should have at least 13 digits of accuracy.  See module
    ! constants
    z=min(eta,max_e_able)
    z=max(z,-max_e_able)

    z = exp( z )
    logit_link = z / (1.0D0 + z)

end function

!-------------------------------------------

double precision function sine_link( eta )
!
!   Compute the inverse of the cosine link. 
!   eta is the linear predictor
!
    use nrtype ! for PI_D constant = value of pie.
    use constants
    implicit none
    
    !double precision :: sine_link
    double precision :: eta
    
    if( eta < -pi_mult ) then
        sine_link = 0.0D0
    else if( eta > pi_mult ) then
        sine_link = 1.0D0
    else 
        sine_link = 0.5D0 + 0.5D0 * sin( (eta * PI_D) / (2.0D0 * pi_mult) ) 
    end if

end function

!-------------------------------------------

double precision function hazard_link( eta )
!
!   Compute the inverse of the hazard link. 
!   eta is the linear predictor
!
    use constants 
    implicit none
    
    !double precision :: hazard_link
    double precision :: eta
    double precision :: z

    z=min(eta,max_e_able)
    z=max(z,-max_e_able)

    hazard_link = 1.0D0 - exp( -exp( z ))

end function

!-------------------------------------------

subroutine location(nan,ns,ic,first,last)
!
!    Purpose: compute first and last capture for every animal
!
!    Input:
!    nan = number of animals
!    ns = number of sample occasions
!    ic = capture history matrix, nan X ns
!
!    Output:
!    first = nan X1 vector containing occasion of first capture
!    last = nan X1 vector containing occasion of last capture

    implicit none

    integer, intent(in) :: nan, ns
    integer, intent(in), dimension(nan,ns) :: ic
    integer, intent(out), dimension(nan) :: first, last

    logical :: findic
    integer :: i, j

    first = 0
    last = 0
    do i=1,nan
        findic=.true.
        do j=1,ns
            if (ic(i,j) >= 1) then
                if (findic .and. (j < ns)) then
                    first(i)=int(j+1,1)
                    findic=.false.
                end if
                last(i)=j
            end if
        end do
    end do

end subroutine

!-------------------------------------------

subroutine first_capture(nan,ns,ic,first)
!
!    Purpose: compute occasion of first capture for every animal.
!    Used by Huggins model
!
!    Input:
!    nan = number of animals
!    ns = number of sample occasions
!    ic = capture history matrix, nan X ns
!
!    Output:
!    first = nan X1 vector containing occasion of first capture
!
!    We assume that at least 1 appears in each row of histories. This
!    condition is checked in R, prior to getting here.

    implicit none

    integer, intent(in) :: nan, ns
    integer, intent(in), dimension(nan,ns) :: ic
    integer, intent(out), dimension(nan) :: first

    integer :: i, j

    first = 0
    do i=1,nan
        do j=1,ns
            if (ic(i,j) == 1) then   ! Careful here, Must be sure only 0's and 1's
                first(i)=j
                exit
            end if
        end do
    end do

end subroutine


! ---------------------------------------------------------------------------------------------

SUBROUTINE VA09AD(FUNCT,N,X,F,G,H,W,DFN,EPS,MODE,MAXFN,IPRINT,IEXIT)
! Parameters in call from CJS: CJS_obj,    np,parameters,loglik,g,h,W, -2.0D0, beta_tol_vec, 1, max_fn, trace, exit_code
! Parameters in call from Hug: Huggins_obj,np,parameters,loglik,g,h,W, -2.0D0, beta_tol_vec, 1, maxfn, trace, exit_code
!
!    Purpose: to maximize the capture-recapture log likelihood using the VA09AD algorithm
!    from HSL.  This is the same routine used by MARK.  See hsl.rl.ac.uk for the HSL library.
!    In general, this will minimize the FUNCT(x).
!
!    Input:
!    funct = name of the function that will calculate value of the log-likelihood
!    n = number of parameters
!    x = current values of the parameters
!    f = value of the log-likelihood (at return)
!    g = gradient vector of size n
!    h = array of size n(n+1)/2 containing the Hessian.  This is represented in the product LDL',
!        where L is lower triangular matrix with unit diagonals and D is a diagonal matrix.
!        The lower triangel of L is stored by columns in H, excepting that the unit
!        diagonal elements are replaced by the corresponding elements of D.  The setting of H
!        on entry is controlled by the parameter MODE.
!    w = array size 3n used for work space
!    dfn = "a number which is set so as to give VA09 an estimate of the likely reduction to be obtained in F(x)."
!        This is used only on the first iteration, so an order of magnitude estimate will suffice.
!        If dfn > 0, the setting of dfn itself will be taken as the likely reduction to be obtained in F(x).
!        If dfn = 0, it will be assumed that an estimate of the minimum value of F(x) has been set in argument
!        f, and the likely reduction in F(x) will be computed according to the initial function value.
!        If dfn < 0, a multiple |dfn| of the modulus of the initial function value will be taken as
!        an estimate of the likely reduction.
!    eps = array of size n containing accuracy required in each element of X.
!    mode = controls setting of the initial estimate of the hessian.
!        If mode = 1, an estimate corresponding to a unit matrix is set in H by VA09.
!        If mode = 2, VA09 assumes that the hessian matrix itself has been set in H by columns of its
!        lower triangle, and the coversion to LDL' form is carried out by VA09.  H must be positive definite.
!        If mode = 3, VA09 assumes that the hessian has been set in H in product form, in which case the
!        the contents of h are passed on unchanged.
!    maxfn = maximum number of calls to funct allowed.
!    iprint = integer controlling printing.  Printing occurs every |iprint| interations and also on exit.  Output
!        is of the form:
!        Function value, x(1), x(2),...x(n), G(1), ...G(n).
!        Nothing is printed when iprint = 0.  Intermediate printing is surpressed with iprint > maxfn+1.
!        values of X and G are surpressed if iprint < 0
!    iexit = integer reason for exiting.
!        if iexit = 0 (mode = 2 only) the estimate of the hessian matrix is not positive definite.
!        if iexit = 1, a normal exit has been made in which |dx(i)| < eps(i) for all i, where dx(i) is
!        change in parameter i on an iteration.
!        if iexit = 2, G'dx >= 0. This is not possible without rounding error.  The probably cause is that
!        eps is set too small for the computer word length.
!        if iexit = 3, funct has been called maxfn times.
!
!
!    Modified by Trent McDonald for Fortran 95 conventions.
    use constants, only: logfile
    implicit none

    double precision, intent(inout) :: DFN,F
    integer, intent(inout) :: IEXIT,IPRINT,MAXFN,MODE, N
    double precision, intent(inout), dimension(n) :: EPS,G,X
    double precision, intent(inout), dimension(n*(n+1)/2) :: H
    double precision, intent(inout), dimension(3*n) :: W
    external FUNCT

    double precision :: ALPHA,DF,DGS,EPSMCH,FY,GS,GS0,GYS,SIG,TOT,Z,ZZ
    INTEGER :: I,ICON,IFN,IG,IGG,IJ,IR,IS,ITN,J,NN
    
    integer :: One = 1
    integer :: Zero = 0
    double precision :: ZeroD = 0.0
    

    double precision, external :: FD05AD

!    EXTERNAL MC11AD,MC11BD,MC11ED
!    INTRINSIC DABS,DSQRT,MOD

!     Following, down to the end of the routine, is unaltered code from HSL, except that I
!    changed all the "0.D0" and the like to "0.0_dbl" and the like.

        
        EPSMCH = FD05AD(1)*10.0D0
!        IF (IPRINT.NE.0) WRITE (logfile,FMT=1000)
 1000   FORMAT (' ENTRY TO VA09AD',/)
        NN = N* (N+1)/2
        IG = N
        IGG = N + N
        IS = IGG
        IEXIT = 0
        IR = N
        IF (MODE.EQ.3) GO TO 15
        IF (MODE.EQ.2) GO TO 10
        IJ = NN + 1
        DO 5 I = 1,N
          DO 6 J = 1,I
            IJ = IJ - 1
            H(IJ) = 0.0D0
    6     CONTINUE     
          H(IJ) = 1.0D0
    5   CONTINUE   
	    
        GO TO 15
   10   CONTINUE
        CALL MC11BD(H,N,IR)
        IF (IR.LT.N) RETURN
   15   CONTINUE
        IF (DFN.EQ.0.0D0) Z = F
        ITN = 0
        CALL FUNCT(N,X,F,G)
        IFN = 1
        DF = DFN
        IF (DFN.EQ.0.0D0) DF = F - Z
        IF (DFN.LT.0.0D0) DF = ABS(DF*F)
        IF (DF.LE.0.0D0) DF = 1.0D0
   20   CONTINUE
        IF (IPRINT.EQ.0) GO TO 21
        IF (MOD(ITN,IPRINT).NE.0) GO TO 21
!        WRITE (logfile,FMT=1001) ITN,IFN,IEXIT,F
! 1001   FORMAT (1x,"Iteration=",I5,1x,"Function evals=",I5,1x,"Exit code=",I5,1x,"Loglik=",D24.16)
!        WRITE (logfile,FMT=1002) "F=", F
! 1002   FORMAT (1x,A,1x,500(1x,5D24.16,/))
!        IF (IPRINT.LT.0) GO TO 21
!        WRITE (logfile,FMT=1002) "X=", (X(I),I=1,N)
!        WRITE (logfile,FMT=1002) "G=", (G(I),I=1,N)
   21   CONTINUE
        ITN = ITN + 1
        DO 22 I = 1,N
          W(IG+I) = G(I)
   22   CONTINUE
        CALL MC11ED(H,N,G,W,IR)
        GS = 0.0D0
        DO 29 I = 1,N
          W(IS+I) = -G(I)
          GS = GS - G(I)*W(IG+I)
   29   CONTINUE
        
        IEXIT = 2
        IF (GS.GE.0.0D0) GO TO 92
        GS0 = GS
        ALPHA = -2.0D0*DF/GS
        IF (ALPHA.GT.1.0D0) ALPHA = 1.0D0
        DF = F
        TOT = 0.0D0

   30   CONTINUE
        IEXIT = 3
        IF (IFN.EQ.MAXFN) GO TO 92
        ICON = 0
        IEXIT = 1
        DO 31 I = 1,N
          Z = ALPHA*W(IS+I)
          IF (ABS(Z).GE.EPS(I)) ICON = 1
          X(I) = X(I) + Z
   31   CONTINUE
        
        CALL FUNCT(N,X,FY,G)
        IFN = IFN + 1
        GYS = 0.0D0
        DO 32 I = 1,N
          GYS = GYS + G(I)*W(IS+I)
   32   CONTINUE
        
        IF (FY.GE.F) GO TO 40
        IF (ABS(GYS/GS0).LE..9D0) GO TO 50
        IF (GYS.GT.0.0D0) GO TO 40
        TOT = TOT + ALPHA
        Z = 10.0D0
        IF (GS.LT.GYS) Z = GYS/ (GS-GYS)
        IF (Z.GT.10.0D0) Z = 10.0D0
        ALPHA = ALPHA*Z
        F = FY
        GS = GYS
        GO TO 30
   40   CONTINUE
        DO 41 I = 1,N
          X(I) = X(I) - ALPHA*W(IS+I)
   41   CONTINUE
        IF (ICON.EQ.0) GO TO 92
        Z = 3.0D0* (F-FY)/ALPHA + GYS + GS
        ZZ = SQRT(Z**2-GS*GYS)
        Z = 1.0D0 - (GYS+ZZ-Z)/ (2.0D0*ZZ+GYS-GS)
        ALPHA = ALPHA*Z
        GO TO 30

   50   CONTINUE
        ALPHA = TOT + ALPHA
        F = FY
        IF (ICON.EQ.0) GO TO 90
        DF = DF - F
        DGS = GYS - GS0
        DO 51 I = 1,N
          W(IGG+I) = G(I)
          G(I) = -W(IG+I)
   51   CONTINUE
        IF (DGS+ALPHA*GS0.GT.0.0D0) GO TO 60
        SIG = 1.0D0/GS0
        IR = -IR
        
        CALL MC11AD(H,N,G,SIG,W,IR,One,ZeroD)
        DO 52 I = 1,N
          G(I) = W(IGG+I) - W(IG+I)
   52   CONTINUE
        SIG = 1.0D0/ (ALPHA*DGS)
        IR = -IR
        CALL MC11AD(H,N,G,SIG,W,IR,Zero,ZeroD)
        GO TO 70
   60   CONTINUE
        ZZ = ALPHA/ (DGS-ALPHA*GS0)
        SIG = -ZZ
        CALL MC11AD(H,N,G,SIG,W,IR,One,EPSMCH)
        Z = DGS*ZZ - 1.0D0
        DO 61 I = 1,N
          G(I) = W(IGG+I) + Z*W(IG+I)
   61   CONTINUE
        SIG = 1.0D0/ (ZZ*DGS**2)
        CALL MC11AD(H,N,G,SIG,W,IR,Zero,ZeroD)
   70   CONTINUE
        DO 71 I = 1,N
          G(I) = W(IGG+I)
   71   CONTINUE
        GO TO 20
   92   CONTINUE
        DO 91 I = 1,N
          G(I) = W(IG+I)
   91   CONTINUE
   90   CONTINUE
!       Trent added the following line so that number of function evaluations is returned. 
!       Careful: this means MAXFN must be a variable in the calling routine, using constant, like 1000, will not work.   
        MAXFN = IFN    
!        IF (IPRINT.EQ.0) RETURN
!        WRITE (logfile,FMT=1001) ITN,IFN,IEXIT,F
!        !WRITE (logfile,FMT=1002) "Final F=", F
!        WRITE (logfile,FMT=1002) "Final X=", (X(I),I=1,N)
!        WRITE (logfile,FMT=1002) "Final G=", (G(I),I=1,N)
        RETURN

end subroutine

! ---------------------------------------------------------------------------------

double precision FUNCTION FD05AD(INUM)
!
!    Purpose: return a real constant for IEEE double precision arithmetic.
!    Code from H.S.L. subroutine ZE02AM, modified by TLM where needed.
    use constants
    implicit none

    !double precision :: FD05AD
    integer, intent(in) :: INUM

    double precision, save, dimension(5) :: DC

    !write(logfile,*) "In fd05ad: dbl=", dbl, "kind(fd05ad)=", kind(fd05ad)

!      SAVE :: DC
!      DC(1) THE SMALLEST POSITIVE NUMBER: 1.0 + DC(1) > 1.0.
!      DC(2) THE SMALLEST POSITIVE NUMBER: 1.0 - DC(2) < 1.0.
!      DC(3) THE SMALLEST NONZERO +VE REAL NUMBER.
!    DC(4) THE SMALLEST FULL PRECISION +VE REAL NUMBER.
!    DC(5) THE LARGEST FINITE +VE REAL NUMBER.

!    I commented these DATA lines out, thinking that huge(), tiny(), etc. would do the job.
!      DATA DC(1)/2.2204460492504D-16/
!      DATA DC(2)/1.1102230246253D-16/
!C     DATA DC(3)/4.9406564584126D-324/  ! this one was commented out in the origninal code from HSL
!      DATA DC(4)/2.2250738585073D-308/
!      DATA DC(5)/1.7976931348622D+308/

      DC(1) = epsilon( DC(1) )
      DC(2) = epsilon( DC(1) )   ! Not sure dc(2) and dc(3) are right. but never use here, so okay.
      DC(3) = tiny( DC(1) )
      DC(4) = tiny( DC(1) )
      DC(5) = huge( DC(1) )

      IF ( INUM .LE. 0 ) THEN
         FD05AD = DC( 1 )
      ELSE IF ( INUM .GE. 6 ) THEN
         FD05AD = DC( 5 )

!      I commented these out because I have a valid dc(3) value
!      ELSE IF ( INUM .EQ. 3 ) THEN
!         FD05AD = DC(4)/2.0D0**52

      ELSE
         FD05AD = DC( INUM )
      ENDIF

end function

! ---------------------------------------------------------------------------------

SUBROUTINE MC11AD(A,N,Z,SIG,W,IR,MK,EPS)
!
!    Purpose: No idea.  Something to do with hessian and gradient and checking for convergence
!
      implicit none

      !integer, parameter :: dbl = selected_real_kind(p=13,r=200)

      integer, intent(inout) :: IR, n, MK
      double precision, intent(inout), dimension(n) :: Z
      double precision, intent(inout), dimension(n*(n+1)/2) :: A
      double precision, intent(inout), dimension(3*n) :: W
      double precision, intent(inout) :: sig, eps


      double precision :: AL,B,GM,R,TI,TIM,V,Y
      INTEGER :: I,IJ,IP,J,MM,NP

      IF (N.GT.1) GO TO 1
      A(1) = A(1) + SIG*Z(1)**2
      IR = 1
      IF (A(1).GT.0.0D0) RETURN
      A(1) = 0.D0
      IR = 0
      RETURN
    1 CONTINUE
      NP = N + 1
      IF (SIG.GT.0.0D0) GO TO 40
      IF (SIG.EQ.0.0D0 .OR. IR.EQ.0) RETURN
      TI = 1.0D0/SIG
      IJ = 1
      IF (MK.EQ.0) GO TO 10
      DO 7 I = 1,N
        IF (A(IJ).NE.0.0D0) TI = TI + W(I)**2/A(IJ)
        IJ = IJ + NP - I
    7 CONTINUE
      GO TO 20
   10 CONTINUE
      DO 11 I = 1,N
        W(I) = Z(I)
   11 CONTINUE
      DO 15 I = 1,N
        IP = I + 1
        V = W(I)
        IF (A(IJ).GT.0.0D0) GO TO 12
        W(I) = 0.D0
        IJ = IJ + NP - I
        GO TO 15
   12   CONTINUE
        TI = TI + V**2/A(IJ)
        IF (I.EQ.N) GO TO 14
        DO 13 J = IP,N
          IJ = IJ + 1
          W(J) = W(J) - V*A(IJ)
   13   CONTINUE
   14   IJ = IJ + 1
   15 CONTINUE
   20 CONTINUE
      IF (IR.LE.0) GO TO 21
      IF (TI.GT.0.0D0) GO TO 22
      IF ((MK-1).GT.0.0D0) THEN
        GO TO 23
      ELSE 
        GO TO 40
      END IF
   21 TI = 0.D0
      IR = -IR - 1
      GO TO 23
   22 TI = EPS/SIG
      IF (EPS.EQ.0.0D0) IR = IR - 1
   23 CONTINUE
      MM = 1
      TIM = TI
      DO 30 I = 1,N
        J = NP - I
        IJ = IJ - I
        IF (A(IJ).NE.0.0D0) TIM = TI - W(J)**2/A(IJ)
        W(J) = TI
        TI = TIM
   30 CONTINUE
      GO TO 41
   40 CONTINUE
      MM = 0
      TIM = 1.0D0/SIG
   41 CONTINUE
      IJ = 1
      DO 66 I = 1,N
        IP = I + 1
        V = Z(I)
        IF (A(IJ).GT.0.0D0) GO TO 53
        IF (IR.GT.0 .OR. SIG.LT.0.0D0 .OR. V.EQ.0.0D0) GO TO 52
        IR = 1 - IR
        A(IJ) = V**2/TIM
        IF (I.EQ.N) RETURN
        DO 51 J = IP,N
          IJ = IJ + 1
          A(IJ) = Z(J)/V
   51   CONTINUE
        RETURN
   52   CONTINUE
        TI = TIM
        IJ = IJ + NP - I
        GO TO 66
   53   CONTINUE
        AL = V/A(IJ)
        IF (MM.GT.0.0D0) THEN
          GO TO 55
        ELSE
          GO TO 54
        END IF
   54   TI = TIM + V*AL
        GO TO 56
   55   TI = W(I)
   56   CONTINUE
        R = TI/TIM
        A(IJ) = A(IJ)*R
        IF (R.EQ.0.0D0) GO TO 70
        IF (I.EQ.N) GO TO 70
        B = AL/TI
        IF (R.GT.4.0D0) GO TO 62
        DO 61 J = IP,N
          IJ = IJ + 1
          Z(J) = Z(J) - V*A(IJ)
          A(IJ) = A(IJ) + B*Z(J)
   61   CONTINUE
        GO TO 64
   62   GM = TIM/TI
        DO 63 J = IP,N
          IJ = IJ + 1
          Y = A(IJ)
          A(IJ) = B*Z(J) + Y*GM
          Z(J) = Z(J) - V*Y
   63   CONTINUE
   64   CONTINUE
        TIM = TI
        IJ = IJ + 1
   66 CONTINUE
   70 CONTINUE
      IF (IR.LT.0) IR = -IR
      RETURN
end subroutine

! -------------------------------------------------------------------------------------------


SUBROUTINE MC11BD(A,N,IR)
!
!     Purpose: Factorize the matrix A into L'DL format.  Upon input, the lower triangle of a
!    symetric matrix H is in stored
!    by columns in A, that is A(1) = H(1,1); A(2) = H(2,1); A(3) = H(3,1), ..., A(N) = H(N,1)
!    A(N+1)=H(2,2), A(N+2) = H(3,2), ..., A(2N-1) = H(N,2), A(2N)=H(3,3), ..., A(N(N+1)/2) = H(N,N).
!
!    Upon output, A contains the lower triagular matrix L and the diagonal matrix D. That is,
!    D(i,i) = A( (i-1)N - (i-2)(i-1)/2 + 1 ), and L is stored everywhere else
!
!    MC11BD and MC11CD are inverses of one another
!
!    Parameters:
!    A = matrix to hold the hessian, converted to L'DL format
!    n = size of hessian if it were square, A holds upper half of hessian, which is size n(n+1)/2
!    ir = a returned code, as near as I can figure.
!
      use constants
      implicit none

      INTEGER, intent(inout) :: IR, n
      double precision, intent(inout), dimension(n*(n+1)/2) :: A

      double precision ::  AA,V
      INTEGER :: I,II,IJ,IK,IP,JK,NI,NP

      IR = N
      IF (N.GT.1) GO TO 100
      IF (A(1).GT.0.0D0) RETURN
      A(1) = 0.D0
      IR = 0
      RETURN
  100 CONTINUE
      NP = N + 1
      II = 1
      DO 104 I = 2,N
        AA = A(II)
        NI = II + NP - I
        IF (AA.GT.0.0D0) GO TO 101
        A(II) = 0.D0
        IR = IR - 1
        II = NI + 1
        GO TO 104
  101   CONTINUE
        IP = II + 1
        II = NI + 1
        JK = II
        DO 103 IJ = IP,NI
          V = A(IJ)/AA
          DO 102 IK = IJ,NI
            A(JK) = A(JK) - A(IK)*V
            JK = JK + 1
  102     CONTINUE
          A(IJ) = V
  103   CONTINUE
  104 CONTINUE
      IF (A(II).GT.0.0D0) RETURN
      A(II) = 0.D0
      IR = IR - 1
      RETURN
end subroutine

! -------------------------------------------------------------------------------------------

SUBROUTINE MC11CD(A,N)
!
!    Purpose: to take matrix A, which is in LDL' format, and muliply it out to get A.
!    That is, on input, we have LDL'.  This will compute A = LDL'.
!
!
      use constants
      implicit none

      integer, intent(inout) :: N
      double precision, intent(inout), dimension(N*(N+1)/2) :: A
      double precision :: AA, V
      INTEGER :: II,IJ,IK,IP,JK,NI,NIP,NP

      IF (N.EQ.1) RETURN
      NP = N + 1
      II = N*NP/2
      DO 202 NIP = 2,N
        JK = II
        NI = II - 1
        II = II - NIP
        AA = A(II)
        IP = II + 1
        IF (AA.GT.0.0D0) GO TO 203
        DO 204 IJ = IP,NI
          A(IJ) = 0.D0
  204   CONTINUE
        GO TO 202
  203   CONTINUE
        DO 201 IJ = IP,NI
          V = A(IJ)*AA
          DO 200 IK = IJ,NI
            A(JK) = A(JK) + A(IK)*V
            JK = JK + 1
  200     CONTINUE
          A(IJ) = V
  201   CONTINUE
  202 CONTINUE
      RETURN
end subroutine

! -------------------------------------------------------------------------------------------

SUBROUTINE MC11ED(A,N,Z,W,IR)
!
!     Purpose: I don't really know, but has something to do with the gradient and hessian
!
!    Parameters:
!    A = matrix to hold the hessian, in weird format
!    n = size of hessian if it were square, A holds upper half of hessian, which is size n(n+1)/2
!    z = gradient vector
!    w = work vector
!    ir = a returned code, as near as I can figure.
!
      implicit none

      !integer, parameter :: dbl = selected_real_kind(p=13,r=200)

      INTEGER, intent(inout) :: IR,n
      double precision, intent(inout), dimension(n) :: Z
      double precision, intent(inout), dimension(n*(n+1)/2) :: A
      double precision, intent(inout), dimension(3*n) :: W

      double precision :: V
      INTEGER :: I,I1,II,IJ,IP,J,NIP,NP

      IF (IR.LT.N) RETURN
      W(1) = Z(1)
      IF (N.GT.1) GO TO 400
      Z(1) = Z(1)/A(1)
      RETURN
  400 CONTINUE
      DO 402 I = 2,N
        IJ = I
        I1 = I - 1
        V = Z(I)
        DO 401 J = 1,I1
          V = V - A(IJ)*Z(J)
          IJ = IJ + N - J
  401   CONTINUE
        W(I) = V
        Z(I) = V
  402 CONTINUE
      Z(N) = Z(N)/A(IJ)
      NP = N + 1
      DO 411 NIP = 2,N
        I = NP - NIP
        II = IJ - NIP
        V = Z(I)/A(II)
        IP = I + 1
        IJ = II
        DO 410 J = IP,N
          II = II + 1
          V = V - A(II)*Z(J)
  410   CONTINUE
      Z(I) = V
  411 CONTINUE
      RETURN
end subroutine

! ---------------------------------------------------------------------------------------------

SUBROUTINE MC11FD(A,N,IR)
!
!    Compute the inverse of a factorized matrix
!
      use constants
      implicit none

      INTEGER, intent(inout) :: IR,N
      double precision, intent(inout), dimension(n*(n+1)/2) :: A
      double precision :: AA,V
      INTEGER :: I,I1,II,IJ,IK,IP,J,JK,K,N1,NI,NIP,NP
      IF (IR.LT.N) RETURN
      A(1) = 1.0D0/A(1)
      IF (N.EQ.1) RETURN
      NP = N + 1
      N1 = N - 1
      II = 2
      DO 511 I = 2,N
        A(II) = -A(II)
        IJ = II + 1
        IF (I.EQ.N) GO TO 502
        DO 501 J = I,N1
          IK = II
          JK = IJ
          V = A(IJ)
          DO 500 K = I,J
            JK = JK + NP - K
            V = V + A(IK)*A(JK)
            IK = IK + 1
  500     CONTINUE
          A(IJ) = -V
          IJ = IJ + 1
  501   CONTINUE
  502   CONTINUE
        A(IJ) = 1.0D0/A(IJ)
        II = IJ + 1
        AA = A(IJ)
        IJ = I
        IP = I + 1
        NI = N - I
        DO 509 J = 2,I
          V = A(IJ)*AA
          IK = IJ
          K = IJ - IP + J
          I1 = IJ - 1
          NIP = NI + IJ
          DO 510 JK = K,I1
            A(JK) = A(JK) + V*A(IK)
            IK = IK + NIP - JK
  510     CONTINUE
          A(IJ) = V
          IJ = IJ + NP - J
  509   CONTINUE
  511 CONTINUE
      RETURN
end subroutine


! ---------------------------------------------------------------------------------------------

subroutine CJS_probs_and_vars(nan,ns,np,parameters,covariance,p_hat,s_hat, se_p_hat,se_s_hat)
!
!    Purpose: Compute variance of estimated prob of capture and survival for all animals
!
!    Input
!    nan,ns,nx,ny,parameters,capX,surX,covariance,p_hat,s_hat
!
!    Output:
!    se_p_hat,se_s_hat
    use constants
    use globevars
    implicit none


    integer, intent(inout) :: np, nan, ns

    double precision, intent(inout), dimension(np) :: parameters
    double precision, intent(inout), dimension(np, np) :: covariance
    double precision, intent(inout), dimension(nan,ns) :: p_hat, s_hat, se_p_hat, se_s_hat

    integer :: i,j,k,l
    double precision :: sum, x1, x2, x
    double precision, dimension(ptr_nx) :: cap_beta
    double precision, dimension(ptr_ny) :: surv_beta


    cap_beta = parameters(1:ptr_nx)
    surv_beta = parameters( (ptr_nx+1):np )


    ! Capture probability first
    do i = 1,nan
        do j = 1,ns
            !call procap(p_hat(i,j), i, j, cap_beta, ptr_nx)

            x = 0.0
            SUM=0.0
            DO K=1,ptr_nx
                x1 = ptr_capX(i,j,k)
                x  = x + cap_beta(k)*x1  ! this is eta = linear predictor
                DO L=1,ptr_nx
                    x2 = ptr_capX(i,j,L)
                    SUM=SUM + (X1 * X2 * covariance(K,L))
                end do
            end do

            call ilink_n_se( x, SUM, p_hat(i,j), se_p_hat(i,j) )
                            
        end do
    end do

    ! Survival probability next
    do i = 1,nan
        do j = 1,ns

            x = 0.0
            SUM=0.0
            DO K=1,ptr_ny
                x1 = ptr_survX(i,j,k)
                x  = x + surv_beta(k)*x1  ! this is eta = linear predictor
                DO L=1, ptr_ny
                    x2 = ptr_survX(i,j,L)
                    SUM=SUM + (X1 * X2 * covariance(k+ptr_nx,L+ptr_nx))
                end do
            end do

            call ilink_n_se( x, SUM, s_hat(i,j), se_s_hat(i,j) )

        end do
    end do

end subroutine

! ----------------------------------------------------------------------------------------------

subroutine est_n_hat(nan,ns,np,covariance,p_hat,se_p_hat,nhat_v_meth, n_hat,se_n_hat)
!
!    Routine to compute CJS Horvitz-Thompson estimates of N hat and varainces of N hats.
!

    use constants
    use globevars
    implicit none

    integer, intent(inout) :: np, nan, ns, nhat_v_meth

    double precision, intent(inout), dimension(np, np) :: covariance
    double precision, intent(inout), dimension(nan,ns) :: p_hat, se_p_hat
    double precision, intent(inout), dimension(ns) :: n_hat, se_n_hat

    integer :: i,j, i1, i2
    double precision :: sum1, sum2, vp, fact
    double precision, external :: phat_cov


    ! --- Compute N hats
    n_hat = 0.0
    do j = 1, ns
        do i = 1, nan
            if (p_hat(i,j) <= 0.0 ) then
                n_hat(j) = -1.0
                exit
            else if ( ptr_hist(i,j) >= 1 )then
                n_hat(j) = n_hat(j) + (1.0D0 / p_hat(i,j))
            end if
        end do
    end do

    ! ---- Compute standard errors
    se_n_hat = 0.0

    select case( nhat_v_meth )


        case (3)
            ! McDonald and Amstrup method
            do j = 1, ns
                sum1 = 0.0
                do i = 1,nan
                    if( ptr_hist(i,j) <= 0 ) cycle
                    vp = se_p_hat(i,j)*se_p_hat(i,j)
                    sum1 = sum1 + ((1.0D0 - p_hat(i,j)) / p_hat(i,j)**2) + (vp / p_hat(i,j)**4)
                end do

                se_n_hat(j) = sqrt( sum1 )

            end do


        case (:2, 4:) ! method == 1 or 2, or anything besides 3

            do j = 1, ns
                sum1 = 0.0
                sum2 = 0.0
                do i = 1,nan
                    if (ptr_hist(i,j) <= 0) cycle
                    if( nhat_v_meth == 2 ) then
                        vp = se_p_hat(i,j)*se_p_hat(i,j)
                        fact = 1.0/p_hat(i,j) + 3*vp/p_hat(i,j)**3 + vp**2/p_hat(i,j)**5
                    else
                        fact = 1.0/p_hat(i,j)
                    end if
                    sum1 = sum1 + fact*(1-p_hat(i,j))/p_hat(i,j)
                end do

                do i1 = 1,nan
                    if (ptr_hist(i1,j) <= 0) cycle
                    do i2 = 1,nan
                        if (ptr_hist(i2,j) <= 0) cycle
                        vp = phat_cov(nan, ns, np, p_hat, covariance, j,i1,i2)
                        sum2 = sum2 + vp / (p_hat(i1,j)*p_hat(i1,j)*p_hat(i2,j)*p_hat(i2,j))
                    end do
                end do

                se_n_hat(j) = sqrt( sum1 + sum2 )
            end do

    end select

end subroutine


! ----------------------------------------------------------------------------------------------
subroutine huggins_pc_hat(nan,ns,nx,ny,np,parameters,covariance,p_hat,se_p_hat,c_hat,se_c_hat)
!
!    Calculate probability of initial and subsequent captures, and SE for both
!
!    Input
!    nan,ns,nx,ny,parameters,capX,surX,covariance,p_hat,s_hat
!
!    Output:
!    se_p_hat,se_s_hat,p_hat,c_hat
!
    use constants
    use globevars
    implicit none


    integer, intent(inout) :: np, nan, ns, nx, ny

    double precision, intent(inout), dimension(np) :: parameters
    double precision, intent(inout), dimension(np, np) :: covariance
    double precision, intent(inout), dimension(nan,ns) :: p_hat, c_hat, se_p_hat, se_c_hat

    integer :: i,j,k,l
    double precision :: sum, x1, x2, x


    ! Initial capture probability first
    do i = 1,nan
        do j = 1,ns
            !call procap(p_hat(i,j), i, j, parameters, nx)

            x = 0.0
            SUM=0.0
            DO K=1,nx
                x1 = ptr_capX(i,j,k)
                x = x + parameters(k)*x1
                DO L=1,nx
                    x2 = ptr_capX(i,j,L)
                    SUM = SUM + (X1 * X2 * covariance(K,L))
                end do
            end do

            call ilink_n_se( x, SUM, p_hat(i,j), se_p_hat(i,j) )

        end do
    end do

    ! Recapture probability next
    do i = 1,nan
        do j = 1,ns

            !call prorecap(c_hat(i,j), i, j, parameters, nx, ny, ptr_remove)

            x = 0.0
            SUM=0.0
            DO K=1,nx+ny

                if (k <= nx) then
                    if (ptr_remove(k) == 1) cycle
                    x1 = ptr_capX(i,j,k)
                else
                    x1 = ptr_capY(i,j,k-nx)
                end if

                x = x + parameters(k)*x1
 
                DO L=1, nx+ny

                    if (L <= nx) then
                        if (ptr_remove(L) == 1) cycle
                        x2 = ptr_capX(i,j,L)
                    else
                        x2 = ptr_capY(i,j,L-nx)
                    end if

                    SUM=SUM + (X1 * X2 * covariance(k,L))
                end do
            end do

            call ilink_n_se( x, SUM, c_hat(i,j), se_c_hat(i,j) )

        end do
    end do



end subroutine

! ----------------------------------------------------------------------------------------------

subroutine ilink_n_se( xbeta, x_vbeta_x, prob, se_prob )
!
!   Invert the link function, and compute variance of inverse link
!
!   xbeta = linear predictor
!   x_vbeta_x = x' V(beta) x quadratic form. assumed 1x1 scaler
!   prob = output probability = answer resulting from inverting link
!   se_prob = output standard error of probability
!
    use nrtype
    use globevars
    implicit none
    
    double precision, intent(in) :: xbeta, x_vbeta_x
    double precision, intent(out) :: prob, se_prob
    double precision, external :: logit_link, sine_link, hazard_link
    
    double precision :: sq_xvx

    sq_xvx = sqrt( max( 0.0D0, x_vbeta_x ) )
    
    if( link == 1 ) then
        prob = logit_link( xbeta )
        se_prob = prob * (1.0D0 - prob) * sq_xvx
    else if( link == 2 ) then
        prob = sine_link( xbeta )
        if( xbeta < -pi_mult .or. xbeta > pi_mult ) then
            se_prob = 0.0D0
        else
            se_prob = sq_xvx * PI_D * cos( xbeta * PI_D / (2.0D0 * pi_mult) ) / (4*pi_mult)
        end if
    else if( link == 3 ) then
        prob = hazard_link( xbeta )
        se_prob = prob * exp(xbeta) * sq_xvx
    else 
        prob = -1.
        se_prob = -1.
    end if
    
end subroutine

! ----------------------------------------------------------------------------------------------

subroutine huggins_n_hat(nan,ns,np,nx,beta,covariance,p_hat,nhat_v_meth, n_hat,se_n_hat,lci,uci)
!
!    Calculate Huggins estimate of N-hat and variance.
!
!
!    Routine to compute estimates of N hat and varainces of N hats.
!

    use constants
    use globevars
    implicit none

    integer, intent(inout) :: np, nan, ns, nx, nhat_v_meth

    double precision, intent(inout), dimension(np) :: beta
    double precision, intent(inout), dimension(np, np) :: covariance
    double precision, intent(inout), dimension(nan,ns) :: p_hat
    double precision, intent(inout) :: n_hat, se_n_hat, lci, uci

    integer :: i,j, p
    double precision :: p_hat_j, p_cap_1, p_cap_2, var_nhat, part1, c, se, deltai
    double precision, dimension(nan) :: p_cap
    double precision, external :: phat_cov
    double precision, dimension(nx) :: p_beta, g
    double precision, dimension(nx,nx) :: v

    ! --- Compute N hat
    n_hat = 0.0
    p_cap = 1.0
    do i = 1, nan
        do j = 1, ns
            p_cap(i) = p_cap(i) * (1.0D0 - p_hat(i,j))
        end do
        p_cap(i) = 1.0D0 - p_cap(i)
        n_hat = n_hat + (1.0D0 / p_cap(i))
    end do

    ! ---- Compute standard error of n_hat
    se_n_hat = -1.0

    p_beta = beta(1:nx)

    select case( nhat_v_meth )

        case (2)

            ! Future bootstrap variance estimate of variance, but no method 2 yet.

        case (:1, 3:) ! method == 1 or anything besides 2

            ! I had to look at the MARK code to do this.  See 'estmat.f', line 4702 where Gary does this
            var_nhat = 0.D0
            g = 0.D0

            do i=1,nan

              var_nhat = var_nhat + (1.D0 - p_cap(i)) / (p_cap(i) * p_cap(i))

              ! Compute gradient of p* wrt each beta
              do p=1,nx


                deltai = (deltax / 2.D0) * (1.D0 + abs(beta(p)))*1.D5

                ! Positive offset
                p_beta(p) = beta(p) + deltai
                p_cap_1 = 1.D0
                do j = 1,ns
                    call procap(p_hat_j, i, j, p_beta, nx)
                    p_cap_1 = p_cap_1 * (1.0D0 - p_hat_j)
                end do
                p_cap_1 = 1.D0 - p_cap_1

                ! Negative offset
                p_beta(p) = beta(p) - deltai
                p_cap_2 = 1.D0
                do j = 1,ns
                    call procap(p_hat_j, i, j, p_beta, nx)
                    p_cap_2 = p_cap_2 * (1.0D0 - p_hat_j)
                end do
                p_cap_2 = 1.D0 - p_cap_2

                g(p) = g(p) + ((p_cap_1 - p_cap_2) / (2.D0 * deltai)) / (p_cap(i) * p_cap(i))

              end do
            end do
            v = covariance(1:nx, 1:nx)
            var_nhat = var_nhat + dot_product(matmul(g,v),g)

            se_n_hat = sqrt( max(0.D0, var_nhat) )

            !  LOG-Based CI's
            part1=max(min_log_able, n_hat-nan )
            se = sqrt(max(0.D0, log(1.D0 + (se_n_hat/part1)**2)))
            c=exp(max(min_e_able, min(max_e_able, 1.96*se)))

            lci=part1/c
            uci=part1*c
            lci=lci + nan
            uci=uci + nan

    end select


end subroutine

! ----------------------------------------------------------------------------------------------
!
!   OKAY, I ADMIT I CANNOT REMEMBER WHY I WROTE THESE GOF ROUTINES HERE.  THEY ARE 
!   ALL IMPLEMENTED IN PURE R CODE IN F.CJS.GOF().  I MUST HAVE STARTED HERE THEN 
!   SWITCHED TO ALL R.  iN ANY EVENT, I AM COMMENTING THEM ALL OUT, JUST IN CASE I NEED 
!   THEM LATER.
!subroutine mragof( nan, ns, hist, p_hat, s_hat, resid_type, &
!        t4_table, t4_chi, t4_df, &
!        t5_table, t5_chi, t5_df, &
!        HL_table, HL_chi, HL_df, &
!        marray_table, marray_chi, marray_df, &
!        roc, residuals, released, cell_expected, &
!        input_trace )
!!
!!   Compute various measures of goodness-of-fit for a capture recapture model.
!!   This routine is designed to be exposed in the DLL, and called after estimation
!!   in mrawin.
!!
!!   Input:
!!   nan = number of animals = number of rows in capture histories
!!   ns = number of capture occasions/samples
!!   hist = capture history matrix (size nan X ns)
!!   p_hat = matrix of probability of capture estimates (size nanXns, first column not used)
!!   s_hat = matrix of survival probabilities (size nanXns, last column never used)
!!   resid_type = Type of residuals to compute.  0 (default) produces Pearson residuals.
!!           1 produces deviance residuals.
!!   input_trace = whether to write results to log file. See defn of trace in globevars module.
!!
!!   Output:
!!   t4_chi = Test 4 chi square value.  Test 4 computes an expected number of captures each occasion.
!!   t4_df = Test 4 degrees of freedom.
!!   t5_chi = Test 5 Chi square value.  Test 5 computes expected number of captures for each animal, then sums over animals.
!!   t5_df = Test 5 degrees of freedom.
!!   HL_chi = Hosmer-Lemeshow chi square value.  Computed over all animals and occasions
!!   HL_df = Hosmer-Lemeshow degrees of freedom.
!!   ROC = overall ROC statistic.
!!   residuals = matrix of residuals (nanXns)
!!   marray_table = observed, expected, contribution and use cell for marry_test
!!   marray_chi = chi square statistic for marray test
!!   marray_df = df for the marray test
!!   released = number released each occasion.  Part of the m-array
!!   cell_expected = expected value of each capture-indicator in the capture history, following
!!       initial encounter.  Return this because its possble to pool over different things
!!       and construct your own test.
!!
!!   Some other gof test, such as the Osius-Rojek, are computed elsewhere (in R or S)
!!   because it was easier there.
!!
!    use constants
!    use globevars
!    implicit none
!
!!   Can't have this if compiling with G95.  Must use if compiling with Layhey lf95.
!!    dll_export mragof
!
!!    Input variables
!    integer, intent(inout), target :: nan, ns
!    integer, intent(inout), dimension(nan,ns), target :: hist
!    double precision, intent(inout), dimension(nan,ns) :: p_hat, s_hat
!    integer, intent(inout) :: resid_type, input_trace
!
!!    Output variables
!    double precision, intent(inout) :: t4_chi, t4_df, t5_chi, t5_df, HL_chi, HL_df, roc, marray_chi, marray_df
!    double precision, intent(inout), dimension(nan,ns) :: residuals, cell_expected
!    double precision, intent(inout), dimension(chi_tab_rows,nan) :: t5_table
!    double precision, intent(inout), dimension(chi_tab_rows,ns)  :: t4_table
!    double precision, intent(inout), dimension(chi_tab_rows,HL_nbins) :: HL_table
!    double precision, intent(inout), dimension(chi_tab_rows,(ns-1)*ns)  :: marray_table
!    integer, intent(inout), dimension(ns)  :: released
!
!!   Local varaibles
!    integer, dimension(nan), target :: first, last
!    integer :: ioerr
!
!!   Store input trace value in globevars module
!    trace = input_trace 
!    
!!    ---- Open a log file to store intermediate tables that may be of use
!    if( trace /= 0 ) then 
!        OPEN(logfile,FILE="mra_gof.log",status="replace",iostat=ioerr)
!        if ( ioerr /= 0 ) then
!            ! Cannot open log file
!            trace = 0
!        end if
!    end if
!
!!   ---- Assign pointers
!    ptr_hist => hist
!
!!   ---- Compute locations of first and last capture.  Needed for tests.
!    call location( nan, ns, hist, first, last)
!    ptr_first => first
!    ptr_last => last
!
!!   ---- Run tests
!
!    call marray_gof(nan, ns, p_hat, s_hat, marray_table, marray_chi, marray_df, released )
!
!    call t4_gof(nan, ns, p_hat, s_hat, resid_type, t4_table, t4_chi, t4_df, residuals, cell_expected )
!
!    call t5_gof(nan, ns, cell_expected, t5_table, t5_chi, t5_df )
!
!    call HL_gof(nan, ns, cell_expected, HL_table, HL_chi, HL_df )
!
!    call roc_gof(nan, ns, cell_expected, roc)
!
!    close(logfile)
!
!end subroutine
!
!! ----------------------------------------------------------------------------------------------
!
!subroutine marray_gof(nan, ns, p_hat, s_hat, fit_table, fit_chisq, fit_chidf, relese )
!!
!!    Purpose: Compute a global goodness of fit statistic using the m-array for
!!   counts.
!!
!!    Input
!!    nan,ns,p_hat,s_hat
!!
!!    Output:
!!    fit_chisq (scalar), fit_chidf (scalar)
!!   fit_table = 4X((ns-1)*ns) array containing observed, expected, contribution, and whether
!!       to use the cell for computing Chi square
!!   relese = number of releases each occasion
!
!    use constants
!    use globevars
!    implicit none
!
!    integer, intent(inout) :: nan, ns
!
!    double precision, intent(inout), dimension(nan,ns) :: p_hat, s_hat
!    double precision, intent(inout) :: fit_chisq, fit_chidf
!    double precision, intent(inout), dimension(chi_tab_rows,(ns-1)*ns) :: fit_table
!    integer, intent(inout), dimension(ns) :: relese
!
!    integer :: i,j,k, l
!    double precision :: contrib, o, e
!
!    fit_table = 0.D0
!    fit_chidf = 0.D0
!    fit_chisq = 0.D0
!    relese = 0
!
!    ! ---- Compute the M array
!    DO i=1,nan
!        if (ptr_first(i) > 0 .and. ptr_first(i) <= ns) then  ! don't really need the second condition
!            do j = ptr_first(i)-1, ns-1
!                if (ptr_hist(i,j) >= 1) then
!                    relese(j) = relese(j) + 1
!                    do k = j+1, ns
!                        if (ptr_hist(i,k) >= 1) then
!                               fit_table(orow,(j-1)*ns+k-1) = fit_table(orow,(j-1)*ns+k-1) + 1.D0
!                               exit
!                           end if
!                    end do
!                end if
!             end do
!        end if
!    end do
!
!    ! ---- Compute expected values for counts in M array
!    DO i=1,nan
!        if (ptr_first(i) > 0 .and. ptr_first(i) <= ns) then  ! don't really need the second condition
!            do j = ptr_first(i)-1, ns-1
!                do k = j+1, ns
!                       ! compute contribution to expected value for count m_(jk)
!                       contrib = s_hat(i,j)
!                       do l = j+1, k-1
!                           contrib = contrib * (1.D0 - p_hat(i,l)) * s_hat(i,l)
!                       end do
!                       contrib = contrib * p_hat(i,k)
!                       fit_table(erow,(j-1)*ns+k-1) = fit_table(erow,(j-1)*ns+k-1) + contrib
!                end do
!             end do
!        end if
!    end do
!
!    ! ---- Compute contribution to chi-square = (obs-exp)^2/exp
!    do j = 1,ns-1
!        do k = j+1,ns
!            o = fit_table(orow,(j-1)*ns+k-1)
!            e = fit_table(erow,(j-1)*ns+k-1)
!            fit_table(oerow,(j-1)*ns+k-1) = (o - e)*(o - e) / e
!            if( e >= chi_ruleofthumb ) then
!                fit_table(userow, (j-1)*ns+k-1) = 1.D0
!                fit_chisq = fit_chisq + fit_table(oerow,(j-1)*ns+k-1)
!                fit_chidf = fit_chidf + 1.D0
!            end if
!        end do
!    end do
!
!    if( fit_chidf > 0.D0 ) then
!        fit_chidf = fit_chidf - 1.D0
!    end if
!
!    ! Write table to the log
!    if(trace /= 0) then
!            write(logfile, 100) (i, i=2,ns)
!        100 format(1x, /" Overall Goodness of Fit Based on M-array"/ &
!                        " ========================================"/ &
!                        " OBSERVED", (1x,10(i8,1x)))
!        300 format(1x, i8,(1x,10(f8.1,1x)))
!            do j = 1,ns-1
!                write(logfile, 300) j,(fit_table(orow,(j-1)*ns + k-1), k=2,ns)
!            end do
!        
!            write(logfile, 400) (i, i=2,ns)
!        400 format(1x, /" EXPECTED", (1x,10(i8,1x)))
!            do j = 1,ns-1
!                write(logfile, 300) j,(fit_table(erow,(j-1)*ns + k-1), k=2,ns)
!            end do
!        
!            write(logfile, 600) (i, i=2,ns)
!        600 format(1x, /" CONTRIB ", (1x,10(i8,1x)))
!            do j = 1,ns-1
!                write(logfile, 300) j,(fit_table(oerow,(j-1)*ns + k-1), k=2,ns)
!            end do
!        
!            write(logfile, 800) (i, i=2,ns)
!        800 format(1x, /" USE CELL", (1x,10(i8,1x)))
!            do j = 1,ns-1
!                write(logfile, 300) j,(fit_table(erow,(j-1)*ns + k-1), k=2,ns)
!            end do
!        
!        
!            WRITE(logfile,700) fit_chisq, fit_chidf
!        700 FORMAT( " ======================================="// &
!                     " Overall GOF Test on M-array Chisq =",F14.5, " on ", F7.0, " df")
!    end if
!
!end subroutine
!
!
!! ----------------------------------------------------------------------------------------------
!
!subroutine HL_gof(nan, ns, cell_expected, HL_table, HL_chi, HL_df )
!!
!!    Purpose: Compute Hosmer-Lemeshow goodness of fit statistic.
!!
!!    Input
!!    nan,ns,p_hat,s_hat
!!
!!    Output:
!!   HL_chi = Chi square value of HL statistic, See Hos-Lem, p.147
!!   HL_df = degrees of freedome for the HL stat
!!
!    use constants
!    use globevars
!    implicit none
!
!    integer, intent(in) :: nan, ns
!
!    !double precision, intent(inout), dimension(nan,ns) :: p_hat, s_hat
!    double precision, intent(in), dimension(nan,ns) :: cell_expected
!    double precision, intent(inout) :: HL_chi, HL_df
!    double precision, intent(inout), dimension(chi_tab_rows,HL_nbins) :: HL_table
!
!    integer :: i, j, n, n_per_bin
!    !double precision :: ee, s
!    double precision, dimension(nan*ns) :: p
!    integer, dimension(nan*ns) :: y
!
!!   ---- Reformat expected values
!    n = 0
!    p = 0.D0
!    y = -1
!    do i = 1, nan
!
!        ! This code recalculates cell expected values, but we have already done that
!        ! in routine t4_gof.
!!        s = 1.D0
!!        if (ptr_first(i) > 0 .and. ptr_first(i) <= ns) then  ! don't really need the second condition
!!            do j = ptr_first(i), ns
!!                s = s*s_hat(i,j-1)  ! note last survival is never used
!!                ee = s*p_hat(i,j)   ! note first capture probability never used
!!                n = n + 1
!!                p(n) = ee
!!                y(n) = ptr_hist(i,j)
!!            end do
!!        end if
!
!        ! Find the non-negative expected values, and put them into p and y vectors so can sort properly.
!        do j = 1, ns
!            if( cell_expected(i,j) >= 0 ) then
!                n = n + 1
!                p(n) = cell_expected(i,j)
!                y(n) = ptr_hist(i,j)
!            end if
!        end do
!    end do
!
!!   ---- Sort the p's and compute HL statistic
!    call bubble_sort(n,p,y)
!
!    n_per_bin = nint( real(n) / real(HL_nbins) )
!    HL_table = 0.D0
!    do i = 1,HL_nbins-1
!        do j = 1, n_per_bin
!            HL_table(orow, i) = HL_table(orow, i) + y( (i-1)*n_per_bin + j )
!            HL_table(erow, i) = HL_table(erow, i) + p( (i-1)*n_per_bin + j )
!        end do
!    end do
!    do j = (HL_nbins-1)*n_per_bin + 1, n  ! last cell
!        HL_table(orow, HL_nbins) = HL_table(orow, HL_nbins) + y( j )
!        HL_table(erow, HL_nbins) = HL_table(erow, HL_nbins) + p( j )
!    end do
!
!    HL_chi = 0.D0
!    HL_df  = 0.D0
!    do i = 1,HL_nbins
!        if( HL_table(erow,i) >= chi_ruleofthumb ) then
!            HL_table(oerow,i) = (HL_table(orow,i)-HL_table(erow,i))*(HL_table(orow,i)-HL_table(erow,i)) / HL_table(erow,i)
!            HL_table(userow,i) = 1.D0
!            HL_chi = HL_chi + HL_table(oerow,i)
!            HL_df  = HL_df  + 1.D0
!        end if
!    end do
!    if( HL_df > 0.D0 ) then
!        HL_df = HL_df - 1.D0
!    end if
!
!
!    ! Write table to the log
!    if (trace /= 0) then 
!            write(logfile, 100) (i, i=1,HL_nbins)
!        100 format(1x, /" Hosmer-Lemeshow Goodness of Fit Contingency Table"/ &
!                        " ================================================="/ &
!                        "Prob.Bin", (1x,10(i8,1x)))
!            write(logfile, 300) (HL_table(orow,i), i=1,HL_nbins)
!        300 format(1x, "Observed",(1x,10(f8.1,1x)))
!            write(logfile, 400) (HL_table(erow,i), i=1,HL_nbins)
!        400 format(1x, "Expected",(1x,10(f8.1,1x)))
!            write(logfile, 500) (HL_table(oerow,i), i=1,HL_nbins)
!        500 format(1x, "Contrib.",(1x,10(f8.1,1x)))
!            write(logfile, 600) (HL_table(userow,i), i=1,HL_nbins)
!        600 format(1x, "Use Cell",(1x,10(f8.1,1x)))
!        
!            WRITE(logfile,700) HL_chi, HL_df
!        700 FORMAT(     " =============================================="// &
!                         " Hosmer-Lemeshow GOF Chisq =",F14.5, " on ", F7.0, " df")
!    end if
!
!end subroutine
!
!! ----------------------------------------------------------------------------------------------
!
!subroutine bubble_sort(n,x,y)
!!
!!   sort the vector x, with parallel vector y
!!
!      use constants
!      implicit none
!
!      INTEGER, intent(in) ::  n
!      double precision, intent(inout), dimension(n) :: x
!      integer, intent(inout), dimension(n) :: y
!
!      integer :: tempi
!      double precision :: tempr, smallest
!      integer :: smallpos, cur, i
!
!!     This is an implementation of bubble sort, which is slow, but
!!     easy to code.
!!     EVENTUALLY, IT WOULD BE NICE TO IMPLEMENT QUICKSORT HERE.
!
!      do cur=1,(n-1)
!
!         smallest = x(cur)
!         smallpos = cur
!
!         do i = (cur+1),n
!            if( x(i) < smallest ) then
!               smallest = x(i)
!               smallpos = i
!            endif
!         end do
!
!!        Swap cur and smallest value
!         if( smallpos .ne. cur ) then
!            tempr = x(cur)
!            x(cur) = x(smallpos)
!            x(smallpos) = tempr
!
!            tempi = y(cur)
!            y(cur) = y(smallpos)
!            y(smallpos) = tempi
!         endif
!
!      end do
!
!end subroutine
!
!! ----------------------------------------------------------------------------------------------
!
!subroutine roc_gof(nan, ns, cell_expected, roc )
!!
!!    Purpose: Compute ROC (area under curve) goodness of fit statistic.
!!
!!    Input
!!    nan,ns,p_hat,s_hat
!!
!!    Output:
!!   roc = value of roc statistic, See Hos-Lem, p.160
!!
!    use constants
!    use globevars
!    implicit none
!
!    integer, intent(in) :: nan, ns
!
!!    double precision, intent(inout), dimension(nan,ns) :: p_hat, s_hat
!    double precision, intent(in), dimension(nan,ns) :: cell_expected
!    double precision, intent(inout) :: roc
!
!    integer :: i, j, n_ones, n_zeros
!    double precision :: ee
!    double precision, dimension(nan*ns) :: p_ones, p_zeros
!
!    double precision, parameter :: near_zero = 1e-6
!
!!   ---- Compute expected values
!    n_ones = 0
!    n_zeros = 0
!    p_ones = 0.D0
!    p_zeros = 0.D0
!    do i = 1, nan
!
!!        s = 1.D0
!!        if (ptr_first(i) > 0 .and. ptr_first(i) <= ns) then  ! don't really need the second condition
!            do j = ptr_first(i), ns
!                if (cell_expected(i,j) > 0) then
!
!!                s = s*s_hat(i,j-1)  ! note last survival is never used
!!                ee = s*p_hat(i,j)   ! note first capture probability never used
!
!                    ee = cell_expected(i,j)
!
!                    if ( ptr_hist(i,j) >= 1 ) then
!                        n_ones = n_ones + 1
!                        p_ones( n_ones ) = ee
!                    else
!                        n_zeros = n_zeros + 1
!                        p_zeros( n_zeros ) = ee
!                    end if
!                end if
!            end do
!!        end if
!    end do
!
!    ! --- Compute ROC statistic, this is the same as the Mann-Whitney U
!    roc = 0.D0
!    do i = 1, n_ones
!        do j = 1, n_zeros
!
!            if( abs(p_ones(i) - p_zeros(j)) <= near_zero ) then
!                roc = roc + 0.5D0
!            else if( p_ones(i) > p_zeros(j) ) then
!                roc = roc + 1.D0
!            end if
!
!        end do
!    end do
!    roc = roc / (n_ones * n_zeros)
!
!end subroutine
!
!! ----------------------------------------------------------------------------------------------
!
!subroutine t5_gof(nan, ns, cell_expected, fit_table, fit_chisq, fit_chidf )
!!
!!    Purpose: Compute goodness of fit statistic across individuals for the model,
!!
!!    Input
!!    nan,ns,
!!   cell_expected = nanXns matrix containing cell expected values.  Computed in t4_gof.
!!       cell_expected < 0 for cells that should not be used.
!!
!!    Output:
!!    fit_chisq (scalar), fit_chidf (scalar)
!!   fit_table = 3Xnan array containing observed, expected, and contribution to Chi square
!
!    use constants
!    use globevars
!    implicit none
!
!    integer, intent(in) :: nan, ns
!
!    !double precision, intent(inout), dimension(nan,ns) :: p_hat, s_hat
!    double precision, intent(in), dimension(nan,ns) :: cell_expected
!    double precision, intent(inout) :: fit_chisq, fit_chidf
!    double precision, intent(inout), dimension(chi_tab_rows,ns) :: fit_table
!
!    integer :: i,j
!    !double precision :: s, ee, se_ee, lee
!    double precision :: o, e
!
!    fit_table = 0.D0
!    fit_chidf = 0.D0
!    fit_chisq = 0.D0
!    do i = 1, nan
!
!        ! This code re-computes the expected value of each cell.  But, already did this
!        ! in t4_gof, so don't do it again.
!!        s = 1.D0
!!        if (ptr_first(i) > 0 .and. ptr_first(i) <= ns) then  ! don't really need the second condition
!!            do j = ptr_first(i), ns
!               ! Recall that ptr_first(i) is set to first estimable capture probability,
!                   ! which is never 1 and is the first capture occasion AFTER first capture.
!                   ! Note also that ptr_first(i) == 0 for first capture at last occasion.  Why I did it this way I don't know.
!!                s = s*s_hat(i,j-1)  ! note last survival is never used
!!                ee = s*p_hat(i,j)   ! note first capture probability never used
!!                fit_table(erow,i) = fit_table(erow,i) + ee  ! note first p_hat is never used
!                ! Could do another test here, and assign expected values to e(j-ptr_first(i)+2)
!                ! this would compile obs and expected based on # occasions since first capture, not just capture occasion
!!                if ( ptr_hist(i,j) >= 1 ) then
!!                    fit_table(orow,i) = fit_table(orow,i) + 1.D0
!!                end if
!!            end do
!
!        ! Compute sums over occasions
!        o = 0.D0
!        e = 0.D0
!        do j = 1, ns
!            if( cell_expected(i,j) >= 0 ) then
!                e = e + cell_expected(i,j)
!                if( ptr_hist(i,j) >= 1 ) then
!                    o = o + 1.D0
!                end if
!            end if
!        end do
!        fit_table(orow,i) = o
!        fit_table(erow,i) = e
!        fit_table(oerow,i) = (o - e)*(o - e) / e
!
!        ! Compute the chi-square statistic
!        if (fit_table(erow,i) >= chi_ruleofthumb) then
!            fit_table(userow,i) = 1.D0
!            fit_chisq = fit_chisq + fit_table(oerow,i)
!            fit_chidf = fit_chidf + 1.D0
!        end if
!
!    end do
!
!    if( fit_chidf > 0.D0 ) then
!        fit_chidf = fit_chidf - 1.D0
!    end if
!
!    ! Write table to the log
!    if (trace /= 0) then
!            write(logfile, 100) (i, i=1,nan)
!        100 format(1x, /" Test 5 Goodness of Fit Contingency Table"/ &
!                        " ========================================"/ &
!                        " Animal", (1x,10(i8,1x)))
!            write(logfile, 300) (fit_table(orow,i), i=1,nan)
!        300 format(1x, "Observed",(1x,10(f8.1,1x)))
!            write(logfile, 400) (fit_table(erow,i), i=1,nan)
!        400 format(1x, "Expected",(1x,10(f8.1,1x)))
!            write(logfile, 500) (fit_table(oerow,i), i=1,nan)
!        500 format(1x, "Contrib.",(1x,10(f8.1,1x)))
!            write(logfile, 600) (fit_table(userow,i), i=1,nan)
!        600 format(1x, "Use Cell",(1x,10(f8.1,1x)))
!        
!            WRITE(logfile,700) fit_chisq, fit_chidf
!        700 FORMAT( " ======================================="// &
!                     " Test 5 GOF Chisq =",F14.5, " on ", F7.0, " df")
!    end if
!    
!end subroutine
!
!
!
!! ----------------------------------------------------------------------------------------------
!
!subroutine t4_gof(nan,ns,p_hat,s_hat, resid_type, fit_table, fit_chisq, fit_chidf, &
!    residuals, cell_expected )
!!
!!    Purpose: Compute goodness of fit statistic for this model, given estimated p's and s's
!!   for all animals.
!!
!!    Input
!!    nan,ns,p_hat,s_hat
!!
!!    Output:
!!    fit_chisq (scalar), fit_chidf (scalar), residuals (nan X ns matrix)
!!   fit_table = 3Xns array containing observed, expected, and contribution to Chi square
!
!    use constants
!    use globevars
!    implicit none
!
!    integer, intent(inout) :: nan, ns, resid_type
!    double precision, intent(in), dimension(nan,ns) :: p_hat, s_hat
!
!    double precision, intent(inout), dimension(nan,ns) :: residuals, cell_expected
!    double precision, intent(inout) :: fit_chisq, fit_chidf
!    double precision, intent(inout), dimension(chi_tab_rows,ns) :: fit_table
!
!    integer :: i,j,k, end_occasion
!    double precision :: s, ee, se_ee, lee
!
!    residuals = missing
!    fit_table = 0.D0
!    cell_expected=-1.D0
!
!    do i = 1, nan
!        s = 1.D0
!
!        !write(logfile,*) i, ptr_first(i), (ptr_hist(i,k), k=1,ns)
!
!        if (ptr_first(i) > 0 .and. ptr_first(i) <= ns) then  ! don't reall need the second condition
!            ! Find last occasion we should compute expected value for.  This is ns if animal
!            ! did not die.  It is occasion with the '2' if animal was death on capture.
!            end_occasion = ns
!            do j = ptr_first(i), ns
!                if( ptr_hist(i,j) >= 2 )then
!                    end_occasion = j
!                    exit
!                end if
!            end do
!
!            ! Compute expected value for each capture indicator between first and end_occasion
!            do j = ptr_first(i), end_occasion
!                   ! Recall that ptr_first(i) is set to first estimable capture probability,
!                   ! which is never 1 and is the first capture occasion AFTER first capture.
!                   ! Note also tht ptr_first(i) == 0 for first capture at last occasion.  Why I did it this way I don't know.
!
!                s = s*s_hat(i,j-1)  ! note last survival is never used
!                ee = s*p_hat(i,j)   ! note first capture probability never used
!
!                ! Save the expected value for the cell
!                cell_expected(i,j) = ee
!
!                ! Error check for the log and sqrt in the residual calc
!                if ((0.D0 <= ee) .and. (ee <= 1.D0)) then
!                    se_ee = sqrt(ee*(1.D0 - ee))
!                else
!                    se_ee = tiny(ee)   ! This will cause the residual to be huge
!                end if
!
!                if (ee <= 0.D0) then
!                    lee = tiny(ee)
!                else if (ee >= 1.D0) then
!                    lee = 1.D0 - tiny(ee)
!                else
!                    lee = ee
!                end if
!
!                fit_table(erow,j) = fit_table(erow,j) + ee  ! note first p_hat is never used
!
!                ! Could do another test here, and assign expected values to e(j-ptr_first(i)+2)
!                ! this would compile obs and expected based on # occasions since first capture, not just capture occasion
!
!                if ( ptr_hist(i,j) >= 1 ) then
!                    fit_table(orow,j) = fit_table(orow,j) + 1.D0
!                    if (resid_type == 1) then
!                        residuals(i,j) = sqrt( 2.D0 * abs( log(lee) ))  ! Deviance residual
!                    else
!                        residuals(i,j) = (1.D0 - ee) / se_ee  ! Pearson residual
!                    end if
!                else
!                    if (resid_type == 1) then
!                        residuals(i,j) = -sqrt( 2.D0 * abs( log(1.D0 - lee) ))  ! Deviance residual
!                    else
!                        residuals(i,j) = (0.D0 - ee) / se_ee  ! Pearson residual
!                    end if
!                end if
!
!            end do
!        end if
!    end do
!
!
!    !   Combine cells if needed to prep for the chi square test
!    do j=ns,2,-1
!        if (fit_table(erow,j) < chi_ruleofthumb) then  ! This is the expected value rule-of-thumb cut off
!
!            if ( j > 2 ) then
!                fit_table(erow,j-1) = fit_table(erow,j) + fit_table(erow,j-1)
!                fit_table(orow,j-1) = fit_table(orow,j) + fit_table(orow,j-1)
!                do k = j+1,ns
!                    fit_table(erow,k-1) = fit_table(erow,k)
!                    fit_table(orow,k-1) = fit_table(orow,k)
!                end do
!                fit_table(erow,ns) = -1.D0
!                fit_table(orow,ns) = -1.D0
!            else  ! j == 2
!                if( ns > 2 ) then
!                    fit_table(erow,j) = fit_table(erow,j) + fit_table(erow,j+1)
!                    fit_table(orow,j) = fit_table(orow,j) + fit_table(orow,j+1)
!                    do k = j+1,ns-1
!                        fit_table(erow,k) = fit_table(erow,k+1)
!                        fit_table(orow,k) = fit_table(orow,k+1)
!                    end do
!                    fit_table(erow,ns) = -1.D0
!                    fit_table(orow,ns) = -1.D0
!                else ! ns <= 2
!                    ! do nothing.  df will be < 0
!                end if
!            end if
!        end if
!    end do
!
!
!    ! Compute the chi-square statistic
!    fit_chidf = 0.D0
!    fit_chisq = 0.D0
!    do j=2,ns
!        if (fit_table(erow,j) > 0.D0) then
!            fit_table(userow,j) = 1.D0
!            fit_table(oerow,j) = (fit_table(orow,j) - fit_table(erow,j))*(fit_table(orow,j) - fit_table(erow,j))/fit_table(erow,j)
!            fit_chisq = fit_chisq + fit_table(oerow,j)
!            fit_chidf = fit_chidf + 1.D0
!        else
!            fit_table(userow,j) = 0.D0    ! don't really need these two assignmenst, redundant give initiation above
!            fit_table(oerow,j) = 0.D0
!        end if
!    end do
!
!    if( fit_chidf > 0.D0 ) then
!        fit_chidf = fit_chidf - 1.D0
!    end if
!
!    ! Write table to the log
!    if (trace /= 0) then
!            write(logfile, 100) (j, j=2,ns)
!        100 format(1x, /" Test 4 Goodness of Fit Contingency Table"/ &
!                        " ========================================"/ &
!                        " Occasion",(1x,10(i8,1x)))
!            write(logfile, 300) (fit_table(orow,j), j=2,ns)
!        300 format(1x, "Observed",(1x,10(f8.1,1x)))
!            write(logfile, 400) (fit_table(erow,j), j=2,ns)
!        400 format(1x, "Expected",(1x,10(f8.1,1x)))
!            write(logfile, 500) (fit_table(oerow,j), j=2,ns)
!        500 format(1x, "Contrib.",(1x,10(f8.1,1x)))
!            write(logfile, 600) (fit_table(userow,j), j=2,ns)
!        600 format(1x, "Use Cell",(1x,10(f8.1,1x)))
!        
!            WRITE(logfile,700) fit_chisq, fit_chidf
!        700 FORMAT( " ======================================="// &
!                     " Test 4 GOF Chisq =",F14.5, " on ", F7.0, " df")
!    end if 
!
!end subroutine


! ----------------------------------------------------------------------------------------------

double precision function phat_cov(nan, ns, np, p, cov, j, i1, i2)
!
!    compute the covariance between p(i1,j) and p(i2,j), where p is a probability computed
!    using the logisitic link.
!
!    Inputs:
!    nan = number of animals
!    ns = number of sampleing occasions
!    np = number of parameters
!    p = probability of capture for all animals and occasions
!    cov = covariance matrix for coefficient vector beta
!    j = occasion we are working on
!    i1 = animal 1 that we are working on
!    i2 = animal 2 that we are working on
!

    use constants
    use globevars
    implicit none

    !double precision :: phat_cov
    integer, intent(inout) :: nan, np, ns, j, i1, i2

    double precision, intent(inout), dimension(np, np) :: cov
    double precision, intent(inout), dimension(nan,ns) :: p

    double precision :: cov_eta12
    integer :: a, b

    !    Compute covariance of the eta's
    cov_eta12 = 0.0
    do a = 1,ptr_nx
        do b = 1,ptr_nx
            cov_eta12 = cov_eta12 + ptr_capX(i1,j,a)*ptr_capX(i2,j,b)*cov(a,b)
        end do
    end do

    !    Compute cov of backtransformed probabilities.  This was derived using the delta method, p. 9 of Seber
    phat_cov = cov_eta12 * p(i1,j) * (1.0D0 - p(i1,j)) * p(i2,j) * (1.0D0 - p(i2,j))

end function

! ----------------------------------------------------------------------------------------------

subroutine comp_hessian(FUNCT, np, beta, f, hess)
!
!     purpose: to compute the hessian matrix by differentiation.
!       I don't know whether this is proper terminology, but here I call the
!       matrix of second partial derivatives evaluated at the maximum "the hessian".
!       To get the covariance matrix, you need to multiply this "hessian" matrix
!       by -1, and then invert it.
!
!    input:
!       FUNCT = function to compute covariance of. (e.g., CJS_loglik)
!    np = number of parameters
!    beta = final coefficients in the model
!    f = value of log likelihood at beta
!
!    Output:
!    hess = np X np output hessian matrix
!    ierr = error code returned from inverse procedure

    use constants
    use globevars, only: trace
    implicit none

    integer, intent(inout) :: np
    double precision, intent(inout), dimension(np) :: beta
    double precision, intent(inout) :: f
    double precision, intent(inout), dimension(np,np) :: hess

    double precision, dimension(np) :: h, b
    double precision :: f1, f2, f3, f4
    integer :: i, j
    double precision, external :: FUNCT

    ! I emailed Gary White and asked him where he got the algorithm to compute numeric
    ! 2nd derivatives. He pointed me to the SAS documentation for NLP procedure.
    ! That is where I got the formulas.
    !   Diagonal element is:
    !       hess[i,i] <- (-FUN(x+2*h*ei) + 16*FUN(x+h*ei) - 30*FUN(x) + 16*FUN(x-h*ei) - FUN(x-2*h*ei)) / (12*h[i]*h[i])
    !   Off-Diagonal element is:
    !       hess[i,j] <- (FUN(x+h*ei+h*ej) - FUN(x+h*ei-h*ej) - FUN(x-h*ei+h*ej) + FUN(x-h*ei-h*ej)) / (4*h[i]*h[j])

    ! ---- Compute h vector, different step size in each dimension
    !   Note that eps is stored in module constants, and set in routine set_constants

!    if (trace /= 0) then
!        write(logfile,*)
!        write(logfile,*) "----- Computing Matrix of 2nd derivatives ----"
!    end if

    ! Note: h can be zero when beta is zero. Beta can be exactly zero when
    !       there is a singularity (overparameterization).

    ! This loop computes additive bits using my method from the SAS manuals. Good agreement
    ! with MARK for the right eps.  Wrong eps, much different results.
    !write(logfile,*) "EPS=", eps
    !do i=1,np
    !    h(i) = (eps**(0.25)) * beta(i)
    !    if( abs(h(i)) < eps ) then
    !        if (beta(i) < 0) then
    !            h(i) = -eps
    !        else
    !            h(i) =  eps
    !        end if
    !    end if
    !end do

    ! This loop computes additive bits same way as MARK
    ! if (trace /= 0) write(logfile,*) "DELTAX=", deltax
    do i = 1, np
        h(i) = (deltax / 2.D0) * (1.D0 + abs(beta(i)))*1.D5
    end do

    !if (trace /= 0) write(logfile,*) " ----- h (=delta) vector -----"
    !if (trace /= 0) write(logfile,"(1000(g20.10,','))") (h(i), i=1,np)
    !write(logfile, *) "(input) f= ", f, " CJS_loglik(beta)=", FUNCT(np, beta)

    ! ---- Compute hessian using "symetric" derivatives
    do i = 1, np

        b = beta

        ! ---- Compute diagonal element
        b(i) = beta(i) + 2.0D0 * h(i)
        f1 = FUNCT(np, b)

        b(i) = beta(i) + h(i)
        f2 = FUNCT(np, b)

        b(i) = beta(i) - h(i)
        f3 = FUNCT(np, b)

        b(i) = beta(i) - 2.0D0 * h(i)
        f4 = FUNCT(np, b)

        hess(i,i) = (-f1 + 16.D0 * f2 - 30.D0 * f + 16.D0 * f3 - f4 ) / (12.D0 * h(i) * h(i) )

        ! ----- Compute off-diagonal elements

        if ( (i+1) <= np ) then
            b = beta
            do j = (i+1), np
                b(i) = beta(i) + h(i)
                b(j) = beta(j) + h(j)
                f1 = FUNCT(np, b)

                ! debugging
                !write(logfile,*) " ------- "
                !write(logfile,*) b
                !write(logfile,*) f1

                b(i) = beta(i) + h(i)
                b(j) = beta(j) - h(j)
                f2 = FUNCT(np, b)

                ! debugging
                !write(logfile,*) " ------- "
                !write(logfile,*) b
                !write(logfile,*) f2

                b(i) = beta(i) - h(i)
                b(j) = beta(j) + h(j)
                f3 = FUNCT(np, b)

                ! debugging
                !write(logfile,*) " ------- "
                !write(logfile,*) b
                !write(logfile,*) f3

                b(i) = beta(i) - h(i)
                b(j) = beta(j) - h(j)
                f4 = FUNCT(np, b)

                ! debugging
                !write(logfile,*) " ------- "
                !write(logfile,*) b
                !write(logfile,*) f4

                hess(i,j) = (f1 - f2 - f3 + f4) / (4.D0 * h(i) * h(j))

                ! ----- Assume symetric
                hess(j,i) = hess(i,j)
            end do
        end if

    end do

    ! ---- Multiply hessian by -1.
    hess = -1.0D0 * hess

    ! debugging
    !if (trace >= 1) then
    !    write(logfile,*) " ----- Matrix of 2nd derivatives, prior to inversion -----"
    !    do i = 1, np
    !        write(logfile,"(1000(g20.10,','))") (hess(i,j), j=1,np)
    !    end do
    !end if

end subroutine


! -------------------------------------------------------------------------------------

subroutine syminv(b,n,ierr)
!
!    Purpose: invert a symmetric matrix
!
    use constants
    use globevars, only: trace
    implicit none

    integer, intent(inout) :: n
    double precision, intent(inout), dimension(n,n) :: b
    integer, intent(inout) :: ierr

    integer :: q, i, j, ka, k, m
    double precision, dimension(n*(n+1)/2) :: wa
    double precision, dimension(n) :: wb
    double precision :: s,t

    ! (dimension of wa must be at least n*(n+1)/2)

    ierr = 0
    if (n > 1) goto 10

    !**** quick exit if n=1
    if (b(1,1) == 0.0) goto 100
    b(1,1)=1.0D0 / b(1,1)
    return

    !**** normal calculations with n>1
    10 continue
    do i=1,n
        do j=1,i
            wa(i*(i-1)/2+j)=b(i,j)
        end do
    end do

    !**** start invertion
    do ka=1,n
        k=n+1-ka
        s=wa(1)
        m=1
        do i=2,n
            q=m
            m=m+i
            t=wa(q+1)
            if (s == 0.0) goto 100
            wb(i)=-t/s
            if (i > k) wb(i)=-wb(i)
            do j=q+2,m
                wa(j-i)=wa(j)+t*wb(j-q)
            end do
        end do
        q=q-1
        wa(m)=1.0D0 / s
        do i=2,n
            wa(q+i)=wb(i)
        end do
    end do

    !**** store answer in matrix b
    do i=1,n
        do j=1,i
            b(i,j)=wa(i*(i-1)/2+j)
            b(j,i)=b(i,j)
        end do
    end do
    return

    !**** singular matrix
    100 continue
    !if (trace /= 0) write(logfile,"(/a)") " Singular matrix of 2nd derivatives:-"
    b = 0.0
    ierr=1
    return

end subroutine

! -----------------------------------------------------------------------------------

integer function matrank( x, m, n )
!
!   Return the rank of matrix x, which is m by n.
!
!   This computes SV decomposition, then counts number of SV's > 0.
!
!   Key parameter in this is where to cut off SV's.
    use constants, only : SVD_ZERO, logfile
    use globevars, only: trace

    implicit none

    integer, intent(in) :: n,m
    double precision, intent(in), dimension(m,n) :: x

    ! Locals
    double precision, dimension(n) :: s_vals
    double precision, dimension(m,n) :: v, a
    double precision :: max_s
    integer :: i

    !   Do the singular value decomposition.  All we need is the singular values,
    !   but I don't know how to shut off computation of a and v.  Contents of a
    !   are destroyed, so must work on a copy to prevent changes to var-covar in calling routine.
    a = x
    call svdcmp_dp(a, s_vals, v, m, n)

    !   Check for convergence
    if ( s_vals(1) <= -9 ) then
        matrank = 0
        return
    end if

    !if (trace /= 0) then 
    !    write(logfile,*)
    !    write(logfile,*) " ----- Singular values of Hessian matrix -----"
    !    write(logfile,"(1000(g20.10,','))") (s_vals(i), i=1,n)
    !end if

    !   Compute maximum of singular values. Scale singular values by it
    max_s = maxval( s_vals )
    s_vals = s_vals / max_s

    !if (trace /= 0) then
    !    write(logfile,*) " ----- Conditioned Singular values of Hessian matrix -----"
    !    write(logfile,"(1000(g20.10,','))") (s_vals(i), i=1,n)
    !    write(logfile,*)
    !end if

    !   Count number of singular values > zero.
    matrank = 0
    do i = 1,n
        if( s_vals(i) > SVD_ZERO ) then
            matrank = matrank + 1
        end if
    end do

end function
! ----------------------------------------------------------------------------------

SUBROUTINE svdcmp_dp(a,w,v,m,n)
!
!   Do single value decompsition of a matrix.  Here, we are decomposing the
!   covariance matrix.
!
!   I got this routine straight out of Numerical Recipies, Chapter 2.
!   Here are thier comments.
!       Given an M x N matrix a, this routine computes its singular value decomposition,
!       A = UWV'.  The matrix U replaces a on output.  The diagonal matrix of singular
!       values W is stored in w as a vector of size N.  The NxN matrix V is output as v.
!
    !use constants, only : DP, I4B
    USE nrutil, ONLY : outerprod
    USE nr, ONLY : pythag
    IMPLICIT NONE
    integer :: i,its,j,k,l,m,n,nm
    double precision, DIMENSION(m,n), INTENT(INOUT) :: a
    double precision, DIMENSION(n), INTENT(OUT) :: w
    double precision, DIMENSION(m,n), INTENT(OUT) :: v
    double precision :: anorm,c,f,g,h,s,scale,x,y,z
    double precision, DIMENSION(m) :: tempm
    double precision, DIMENSION(n) :: rv1,tempn
    !m=size(a,1)
    !n=assert_eq(size(a,2),size(v,1),size(v,2),size(w),'svdcmp_dp')
    g=0.0
    scale=0.0
    do i=1,n
        l=i+1
        rv1(i)=scale*g
        g=0.0
        scale=0.0
        if (i <= m) then
            scale=sum(abs(a(i:m,i)))
            if (scale /= 0.0) then
                a(i:m,i)=a(i:m,i)/scale
                s=dot_product(a(i:m,i),a(i:m,i))
                f=a(i,i)
                g=-sign(sqrt(s),f)
                h=f*g-s
                a(i,i)=f-g
                tempn(l:n)=matmul(a(i:m,i),a(i:m,l:n))/h
                a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
                a(i:m,i)=scale*a(i:m,i)
            end if
        end if
        w(i)=scale*g
        g=0.0
        scale=0.0
        if ((i <= m) .and. (i /= n)) then
            scale=sum(abs(a(i,l:n)))
            if (scale /= 0.0) then
                a(i,l:n)=a(i,l:n)/scale
                s=dot_product(a(i,l:n),a(i,l:n))
                f=a(i,l)
                g=-sign(sqrt(s),f)
                h=f*g-s
                a(i,l)=f-g
                rv1(l:n)=a(i,l:n)/h
                tempm(l:m)=matmul(a(l:m,l:n),a(i,l:n))
                a(l:m,l:n)=a(l:m,l:n)+outerprod(tempm(l:m),rv1(l:n))
                a(i,l:n)=scale*a(i,l:n)
            end if
        end if
    end do
    anorm=maxval(abs(w)+abs(rv1))
    do i=n,1,-1
        if (i < n) then
            if (g /= 0.0) then
                v(l:n,i)=(a(i,l:n)/a(i,l))/g
                tempn(l:n)=matmul(a(i,l:n),v(l:n,l:n))
                v(l:n,l:n)=v(l:n,l:n)+outerprod(v(l:n,i),tempn(l:n))
            end if
            v(i,l:n)=0.0
            v(l:n,i)=0.0
        end if
        v(i,i)=1.0
        g=rv1(i)
        l=i
    end do
    do i=min(m,n),1,-1
        l=i+1
        g=w(i)
        a(i,l:n)=0.0
        if (g /= 0.0) then
            g=1.0/g
            tempn(l:n)=(matmul(a(l:m,i),a(l:m,l:n))/a(i,i))*g
            a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
            a(i:m,i)=a(i:m,i)*g
        else
            a(i:m,i)=0.0
        end if
        a(i,i)=a(i,i)+1.0
    end do
    do k=n,1,-1
        do its=1,30
            do l=k,1,-1
                nm=l-1
                if ((abs(rv1(l))+anorm) == anorm) exit
                if ((abs(w(nm))+anorm) == anorm) then
                    c=0.0
                    s=1.0
                    do i=l,k
                        f=s*rv1(i)
                        rv1(i)=c*rv1(i)
                        if ((abs(f)+anorm) == anorm) exit
                        g=w(i)
                        h=pythag(f,g)
                        w(i)=h
                        h=1.0/h
                        c= (g*h)
                        s=-(f*h)
                        tempm(1:m)=a(1:m,nm)
                        a(1:m,nm)=a(1:m,nm)*c+a(1:m,i)*s
                        a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
                    end do
                    exit
                end if
            end do
            z=w(k)
            if (l == k) then
                if (z < 0.0) then
                    w(k)=-z
                    v(1:n,k)=-v(1:n,k)
                end if
                exit
            end if
            if (its == 30) then
                !call nrerror('svdcmp_dp: no convergence in svdcmp')
                w(1) = -999
                return
            end if
            x=w(l)
            nm=k-1
            y=w(nm)
            g=rv1(nm)
            h=rv1(k)
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
            g=pythag(f,1.0D0)
            f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
            c=1.0
            s=1.0
            do j=l,nm
                i=j+1
                g=rv1(i)
                y=w(i)
                h=s*g
                g=c*g
                z=pythag(f,h)
                rv1(j)=z
                c=f/z
                s=h/z
                f= (x*c)+(g*s)
                g=-(x*s)+(g*c)
                h=y*s
                y=y*c
                tempn(1:n)=v(1:n,j)
                v(1:n,j)=v(1:n,j)*c+v(1:n,i)*s
                v(1:n,i)=-tempn(1:n)*s+v(1:n,i)*c
                z=pythag(f,h)
                w(j)=z
                if (z /= 0.0) then
                    z=1.0/z
                    c=f*z
                    s=h*z
                end if
                f= (c*g)+(s*y)
                x=-(s*g)+(c*y)
                tempm(1:m)=a(1:m,j)
                a(1:m,j)=a(1:m,j)*c+a(1:m,i)*s
                a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
            end do
            rv1(l)=0.0
            rv1(k)=f
            w(k)=x
        end do
    end do
    END SUBROUTINE svdcmp_dp


! ----------------------------------------------------------------------------------------------
FUNCTION pythag_dp(a,b)
!
!   This function is needed by svdcmp_dp.  Taken from Numerical Recipies
!   This computes calculate (a^2+b^2)^\{1/2\} without overflow.
    USE nrtype
    IMPLICIT NONE
    double precision, INTENT(IN) :: a,b
    double precision :: pythag_dp
    double precision :: absa,absb
    absa=abs(a)
    absb=abs(b)
    if (absa > absb) then
        pythag_dp=absa*sqrt(1.0+(absb/absa)**2)
    else
        if (absb == 0.0) then
            pythag_dp=0.0
        else
            pythag_dp=absb*sqrt(1.0+(absa/absb)**2)
        end if
    end if
    END FUNCTION pythag_dp

! ----------------------------------------------------------------------------------------------


SUBROUTINE TESTS(NAN,NS,IC,NG,IG,VIF,CHIGT,IDFGT)
!
!     Subroutine to do TEST 2 & TEST 3 from Burnham et al. (1987)
!
!    Note that Bryan Manly and I do these tests slightly differently than RELEASE.
!    We have different rules about Df and when to pool various tables.
!
use constants, only : chat_rot, logfile
use globevars, only : trace
implicit none

    !integer, parameter :: real_dbl=selected_real_kind(p=13,r=200)

    integer, intent(in) :: nan, ns, ng
    INTEGER, intent(in), dimension(nan,ns) :: ic
    INTEGER, intent(in), dimension(nan) :: ig
    double precision, intent(out) :: vif, chigt
    integer, intent(out) :: idfgt

    integer, dimension(ns) :: iuse, idf, relese
    integer, dimension(ns,ns) :: m
    double precision, dimension(ns) :: chisq
    integer :: i, j, k, k1, k2, L, idftot
    integer :: nang   ! number of animals in group
    double precision :: CHITOT, compr, comps, tot
    double precision :: e11, e12, e21, e22
    integer :: N11, N12, N21, N22, M11, M12, M21, M22, R1, R2, C1, C2, ib, ia
    integer :: iuser, iuses, idfr, idfs


!     -------
!    Start TEST 2, carried out group by group
!    Only chi-squared components for which all row and column
!    are 5 or more are used for calculating totals and the
!    estimated variance inflation factor (VIF).


!    if (trace /= 0) then 
!            WRITE(logfile,9005) chat_rot
!        9005     FORMAT(/"                   TESTS OF ASSUMPTIONS"/ &
!                   " Only Test2 and Test 3 chi-square component tables"/ &
!                   " for which all cells are >= ",i2," are used for "/ &
!                   " calculating totals and the estimate variance"/ &
!                   " inflation factor (C-hat)."/)
!            write(logfile,*) "Number of Groups=", ng
!    end if
    CHIGT=0.0
    IDFGT=0.0


    DO L = 1, NG

        ! Initializing
              NANG=0
        RELESE = 0
        M = 0

        ! Computing the M array
              DO I=1,NAN
                  IF (IG(I) == L) THEN
                       NANG=NANG+1
                       DO j=1,NS-1
                           IF (IC(I,J) == 1) THEN
                            RELESE(J)=RELESE(J)+1
                                DO K = J+1, NS
                                IF (IC(I,K) == 1) THEN
                                     M(J,K) = M(J,K)+1
                                     exit
                                END IF
                           end do
                           END IF
                     end do
                  END IF
         end do

        IF (NANG > 0) THEN

!            if (trace /= 0) then
!                WRITE(logfile,9010)
!    9010              FORMAT(/," FIRST RECAPTURES FROM RELEASES (i.e., m-array)")
!                write(logfile,9015) "   i  R(i)", (i, i=2,ns)
!    9015            format(a,1000i4:/(10x,1000i4))
!                write(logfile,"(1000a1)") ("=", i=1,(4*(ns-1)+10))
!                       DO I=1,NS-1
!                           WRITE(logfile,"(I4,I6,1000a4:/(10X,1000a4))",advance="no") I,RELESE(I),(" ", J=2,i)
!                    WRITE(logfile,"(1000I4:/(10X,1000I4))") (M(I,J), J=i+1,NS)
!                end do
!            end if

            CALL TEST2(NS,M,CHISQ,IDF,CHITOT,IDFTOT,IUSE)

!            if (trace /= 0) then
!                WRITE(logfile,9020) (I,CHISQ(I),IDF(I),IUSE(I), I=2,NS-2)
!    9020              FORMAT(/"                                Used for"/ &
!                                    "  Release    Chi-Sq        df      C-hat (1 = yes)"/ &
!                                    " ======================================="/ &
!                                    (I10,F10.2,2I10))
    
    
!                WRITE(logfile,9030) CHITOT,IDFTOT
!    9030              FORMAT( " ======================================="/ &
!                                    "     Total",F10.2,I10)
!            end if

            CHIGT=CHIGT+CHITOT
            IDFGT=IDFGT+IDFTOT
              ENDIF
    end do   ! End of L loop

!    -----
!    TEST 3, carried out group by group

          DO L = 1, NG

!        if (trace /= 0) WRITE(logfile,"(/A,I3)") " TEST 3 FOR GROUP", L

        ! Calculate R and S components for releases 2 to NS-1
        DO J=2,NS-1
            N11=0
            N12=0
            N21=0
            N22=0

            M11=0
            M12=0
            M21=0
            M22=0

            DO I = 1, NAN
                      IF ((IG(I) == L) .AND. (IC(I,J) == 1)) THEN
                    IB=0
                    DO K1=1,J-1
                        IB=IB+IC(I,K1)
                    end do

                    IA=0
                    DO K2=J+1,NS
                        IA=IA+IC(I,K2)
                    end do

                    ! Counts for R component
                    IF ((IB > 0) .AND. (IA > 0)) THEN
                        N11=N11+1
                    ELSEIF ((IB > 0) .AND. (IA == 0)) THEN
                        N12=N12+1
                    ELSEIF ((IB == 0) .AND. (IA > 0)) THEN
                        N21=N21+1
                    ELSE
                        N22=N22+1
                    END IF

                    ! Counts for S component
                    IF (J < NS-1) THEN
                        IF ((IB > 0) .AND. (IC(I,J+1) == 1)) THEN
                            M11=M11+1
                        ELSEIF ((IB > 0) .AND. (IA > 0)) THEN
                            M12=M12+1
                        ELSEIF ((IB == 0) .AND. (IC(I,J+1) == 1)) THEN
                            M21=M21+1
                        ELSEIF ((IB == 0) .AND. (IA > 0)) THEN
                            M22=M22+1
                        END IF
                    end if
                end if
            end do

            ! Chi-squared for R component
            IUSER=1
            R1=N11+N12
            R2=N21+N22
            C1=N11+N21
            C2=N12+N22

            IF ((R1 < chat_rot) .OR. (R2 < chat_rot) .OR. (C1 < chat_rot) .OR. (C2 < chat_rot)) then
                IUSER=0
            end if

            TOT=R1+R2
            IF ((R1 == 0) .OR. (R2 == 0) .OR. (C1 == 0) .OR. (C2 == 0)) THEN
                COMPR=0.0
                IDFR=0
            ELSE
                E11=R1*C1/TOT
                E12=R1*C2/TOT
                E21=R2*C1/TOT
                E22=R2*C2/TOT
                COMPR=(N11-E11)**2/E11+(N12-E12)**2/E12+ &
                    (N21-E21)**2/E21+(N22-E22)**2/E22
                IDFR=1
            END IF

            ! Chi-squared for S component
            IUSES=1
            R1=M11+M12
            R2=M21+M22
            C1=M11+M21
            C2=M12+M22

            IF ((R1 < chat_rot) .OR. (R2 < chat_rot) .OR. (C1 < chat_rot) .OR. (C2 < chat_rot)) then
                IUSES=0
            end if

            TOT=R1+R2
            IF ((R1 == 0).OR.(R2 == 0).OR.(C1 == 0).OR.(C2 == 0)) THEN
                COMPS=0.0
                IDFS=0
            ELSE
                E11=R1*C1/TOT
                E12=R1*C2/TOT
                E21=R2*C1/TOT
                E22=R2*C2/TOT
                COMPS=(M11-E11)**2/E11+(M12-E12)**2/E12+ &
                    (M21-E21)**2/E21+(M22-E22)**2/E22
                IDFS=1
            END IF

            ! Print tests for Jth release
!            if (trace /= 0) WRITE(logfile,9040) J,N11,N12,M11,M12,N21,N22,M21,M22,COMPR,COMPS, &
!                IDFR,IDFS,IUSER,IUSES
!9040             FORMAT(/" Chi-Squared Tests for Animals Captured in Sample",I3/ &
!                                "                          R Test       |        S Test"/ &
!                                "                                       |            Not seen"/ &
!                                "                        Seen  Not seen | Seen next      next"/ &
!                                "                       again     again |    sample    sample"/ &
!                                "============================================================"/ &
!                                "     Seen before",I12,I10,I12,I10/ &
!                                " Not seen before",I12,I10,I12,I10/ &
!                                "============================================================"/ &
!                                "     Chi-squared",12X,F10.2,12X,F10.2/ &
!                                "              df",2I22/ &
!                                " Used in Chat(1)",2I22 )

            CHIGT=CHIGT+COMPR*IUSER+COMPS*IUSES
            IDFGT=IDFGT+IDFR*IUSER+IDFS*IUSES

        end do ! end of J loop
    end do  ! end Test 3 L loop


    ! Print total chi-squared and variance inflation factor
    IF (IDFGT > 0) THEN
        VIF=CHIGT/IDFGT
        VIF=MAX(VIF,1.0D0)
    ELSE
        VIF=1.0D0
    ENDIF

!    if (trace /= 0) WRITE(logfile,9050) CHIGT,IDFGT,VIF
!9050     FORMAT(/" Total chi-squared from used components =",F8.2," with",I3," df"/&
!                " Variance inflation factor (C-hat) =",F7.2)

end subroutine


! ------------------------------------------------------------------------------------

SUBROUTINE TEST2(NS,M,CHISQ,IDF,CHITOT,IDFTOT,IUSE)
!
!     Subroutine to find components 2 to NS-2 of TEST 2 *
!
use constants, only : chat_rot, logfile
use globevars, only : trace

implicit none

    integer, intent(in) :: ns
    integer, intent(in), dimension(ns,ns) :: m

    integer, intent(out), dimension(ns) :: iuse, idf
    double precision, intent(out), dimension(ns) :: chisq
    integer, intent(out) :: idftot
    double precision, intent(out) :: chitot

    integer, dimension(2,ns) :: N
    double precision, dimension(2) :: R
    double precision, dimension(ns) :: C
    double precision :: tot, exp
    integer :: i, j, l

    ! Exit if test not possible
    IF (NS < 4) THEN
!        if (trace /= 0) WRITE(logfile,"(/A)") " TEST 2 not possible!"
        DO I=2,NS-2
            IUSE(I)=0
            CHISQ(I)=0.0
            IDF(I)=0.0
        end do
        RETURN
    END IF

    ! Compute components one by one
    CHITOT=0.0
    IDFTOT=0
    DO L=2,NS-2
        IUSE(L)=1

        ! Find frequencies for contingency table
        DO J=L+1,NS
            N(1,J)=0
            DO I=1,L-1
                N(1,J)=N(1,J)+M(I,J)
            end do
            N(2,J)=M(L,J)
        end do

        ! Find row and column totals
        DO J=L+1, NS
            C(J)=0.0
        end do

        DO I=1,2
            R(I)=0.0
            DO J=L+1,NS
                R(I)=R(I)+N(I,J)
                C(J)=C(J)+N(I,J)
            end do
        end do

        TOT=R(1)+R(2)
        IF ((R(1) < chat_rot) .OR. (R(2) < chat_rot)) then
            IUSE(L)=0
        end if

        DO J=L+1,NS
            IF (C(J) < chat_rot) then
                IUSE(L)=0
            end if
        end do

        ! Find chi-squared
        IF ((R(1) <= 0) .OR. (R(2) <= 0)) THEN
            CHISQ(L)=0.0
            IDF(L)=0
        ELSE
            CHISQ(L)=0.0
            IDF(L)=NS-L-1
            DO J=L+1,NS
                IF (C(J).LE.0) THEN
                    IDF(L)=IDF(L)-1
                ELSE
                    DO I=1,2
                        EXP=R(I)*C(J)/TOT
                        CHISQ(L)=CHISQ(L)+(N(I,J)-EXP)**2/EXP
                    end do
                END IF
            end do

            IF (IDF(L) <= 0) THEN
                IDF(L)=0
                CHISQ(L)=0.0
                IUSE(L)=0
            END IF
        ENDIF

        CHITOT=CHITOT+CHISQ(L)*IUSE(L)
        IDFTOT=IDFTOT+IDF(L)*IUSE(L)
    end do

end subroutine
