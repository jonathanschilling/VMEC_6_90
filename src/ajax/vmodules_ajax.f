      MODULE SPEC_KIND_MOD
!-------------------------------------------------------------------------------
!For REAL variables:
!  Use p=6, r=35   for single precision on 32 bit machine
!  Use p=12,r=100  for double precision on 32 bit machine
!                  for single precision on 64 bit machine
!Parameters for SELECTED_REAL_KIND:
!  p                   -number of digits
!  r                   -range of exponent from -r to r
!-------------------------------------------------------------------------------
!      INTEGER, PARAMETER ::                                                    &
!     & rspec = SELECTED_REAL_KIND(p=6,r=35),                                   &
!     & ispec = SELECTED_INT_KIND(r=8)
      use kind_spec
        integer, parameter :: rspec = rprec
        integer, parameter :: ispec = iprec
      END MODULE SPEC_KIND_MOD


      MODULE ajax_mod
!-------------------------------------------------------------------------------
!AJAX-(Take your pick)
!    -Algebraic Jacobians for Advanced eXperiments
!    -The son of Telamon of Salamis and a warrior of great stature and prowess
!     who fought against Troy
!    -The son of Ileus of Locris and a warrior of small stature and arrogant
!     character who fought against Troy
!    -Brand name of a cleaning agent
!
!AJAX_MOD is an F90 module of routines that serve as an interface between a
!  transport code and MHD equilibrium solutions.
!
!References:
!
!  W.A.Houlberg, P.I.Strand 11/2001
!
!Contains PUBLIC routines:
!
!  Input:
!
!    AJAX_LOAD_RZBDY   -loads approx to 2D MHD equilibria from boundary values
!    AJAX_LOAD_RZLAM   -loads 2D/3D MHD equilibria in inverse coordinate form
!    AJAX_LOAD_MAGFLUX -loads magnetic flux data
!
!  Coordinate conversions:
!
!    AJAX_CAR2CYL      -converts from Cartesian to cylindrical coordinates
!    AJAX_CYL2CAR      -converts from cylindrical to Cartesian coordinates
!    AJAX_CYL2FLX      -converts from cylindrical to flux coordinates
!    AJAX_FLX2CYL      -converts from flux to cylindrical coordinates
!
!  Local magnetics:
!
!    AJAX_B            -gets components of B in various forms
!    AJAX_LAMBDA       -gets stream function and its derivatives
!
!  Flux surface quantities:
!
!    AJAX_FLUXAV_B     -gets flux surface quantities that depend on B
!    AJAX_FLUXAV_G     -gets flux surface quantities that depend on geometry
!    AJAX_I            -gets enclosed toroidal and external poloidal currents
!                       and maximum and minimum values
!    AJAX_MAGFLUX      -gets poloidal and toroidal magnetic fluxes
!    AJAX_SHAPE        -gets shift, elongation, triangularity, Rmax, Rmin, etc
!
!  Miscellaneous:
!
!    AJAX_GLOBALS      -gets global characteristics of the data
!    AJAX_MINMAX_RZ    -gets maximum or minimum of R or Z
!
!Contains PRIVATE routines:
!
!    AJAX_INIT         -initializes or resets private data and its allocations
!                      -called from AJAX_LOAD_RZLAM
!    AJAX_INIT_FLUXAV_B-initializes magnetic field arrays
!                      -called from AJAX_FLUXAV_B
!    AJAX_INIT_FLUXAV_G-initializes geometry arrays
!                      -called from AJAX_FLUXAV_B
!                      -called from AJAX_FLUXAV_G
!    AJAX_INIT_LAMBDA  -calculates the magnetic stream function from R,Z data
!                      -version is only valid for axisymmetric plasmas
!                      -called from AJAX_LOAD_RZLAM
!    AJAX_LOAD_LAMBDA  -loads the magnetic stream function
!                      -called from AJAX_LOAD_RZLAM
!
!Comments:
!
!  This module is designed for use in a transport code that calculates MHD
!    equilibria (flux surface geometry) at infrequent intervals and updates the
!    the magnetic flux (e.g., safety factor or rotational transform) more
!    frequently (e.g., every time step) in a Grad-Hogan approach:
!    1) After a new equilibrium (geometry and magnetic flux) call:
!       AJAX_LOAD_RZLAM or AJAX_LOAD_RZBDY
!       AJAX_LOAD_MAGFLUX
!    2) After each additional timestep (magnetic flux) call:
!       AJAX_LOAD_MAGFLUX
!
!  The flux surface averaging routines are split to provide efficient
!    evaluation:
!    1) After a new equilibrium (geometry dependent) call:
!       AJAX_FLUXAV_G
!       AJAX_SHAPE
!    2) After each timestep (geometry and magnetic flux dependent) call:
!       AJAX_FLUXAV_B
!
!  There is extensive use of optional arguments to:
!    1) Provide flexibility of input variables
!    2) Give the user control over the fineness of the internal grids
!
!  The modernization of the code structure into an F90 module takes advantage of
!    some of the more attractive features of F90:
!    -use of KIND for precision declarations
!    -optional arguments for I/O
!    -generic names for all intrinsic functions
!    -compilation using either free or fixed form
!    -no common blocks or other deprecated Fortran features
!    -dynamic and automatic alocation of variables
!    -array syntax for vector operations
!-------------------------------------------------------------------------------
      USE SPEC_KIND_MOD
      USE LINEAR1_MOD
      USE SPLINE1_MOD
      IMPLICIT NONE

!-------------------------------------------------------------------------------
! Private procedures
!-------------------------------------------------------------------------------
      PRIVATE ::                                                               &
     & AJAX_INIT,                                                              &
     & AJAX_INIT_FLUXAV_B,                                                     &
     & AJAX_INIT_FLUXAV_G,                                                     &
     & AJAX_INIT_LAMBDA,                                                       &
     & AJAX_LOAD_LAMBDA

!-------------------------------------------------------------------------------
! Private data
!-------------------------------------------------------------------------------
!Radial, poloidal and toroidal grids
      INTEGER, PRIVATE, SAVE ::                                                & 
     & nrho_3d,ntheta_3d,nzeta_3d

      REAL(KIND=rspec), PRIVATE, SAVE, ALLOCATABLE ::                          &
     & rho_3d(:),                                                              &
     & wtheta_3d(:),wzeta_3d(:)

      REAL(KIND=rspec), PRIVATE, SAVE ::                                       &
     & dtheta_3d,dzeta_3d

!Poloidal and toroidal mode expansions
      INTEGER, PRIVATE, SAVE ::                                                &
     & krz_3d,klam_3d,                                                         &
     & km0n0_3d,km1n0_3d,                                                      &
     & nper_3d

      INTEGER, PRIVATE, SAVE, ALLOCATABLE ::                                   &
     & m_3d(:),n_3d(:),mabs_3d(:)

!Toroidal flux and magnetic field directions
      REAL(KIND=rspec), PRIVATE, SAVE ::                                       &
     & phitot_3d,signbp_3d,signbt_3d

!Splined quantities in radial coordinate:
      REAL(KIND=rspec), PRIVATE, SAVE, ALLOCATABLE ::                          &
     & iotabar_3d(:,:),                                                        &
     & r_3d(:,:,:),z_3d(:,:,:),lam_3d(:,:,:)

!Arrays for flux surface averaging
      REAL(KIND=rspec), PRIVATE, SAVE, ALLOCATABLE ::                          &
     & az_3d(:),gt_3d(:),vp_3d(:),                                             &
     & rcyl_3d(:,:,:),zcyl_3d(:,:,:),gsqrt_3d(:,:,:),                          &
     & eltheta_3d(:,:,:),elzeta_3d(:,:,:),                                     &
     & b_3d(:,:,:),btheta_3d(:,:,:),                                           &
     & gcyl_3d(:,:,:,:)

!Major and minor radius quantities
      REAL(KIND=rspec), PRIVATE, SAVE ::                                       &
     & r000_3d,                                                                &
     & rhomin_3d,rhomax_3d,rhores_3d

!Logical switches
      LOGICAL, PRIVATE, SAVE ::                                                &
     & l_fluxavb_3d,l_fluxavg_3d,l_mfilter_3d

!Physical and conversion constants
      REAL(KIND=rspec), PRIVATE, SAVE ::                                       &
     & z_large,z_mu0,z_pi,z_precision,z_small

!-------------------------------------------------------------------------------
! Public procedures
!-------------------------------------------------------------------------------
      CONTAINS

      SUBROUTINE AJAX_LOAD_RZBDY(r0,a0,s0,e0,e1,d1,iflag,message,              &
     &                           NRHO_AJAX,NTHETA_AJAX,NKLAM_AJAX,             &
     &                           RHOMAX_AJAX)
!-------------------------------------------------------------------------------
!AJAX_LOAD_RZBDY converts geometric boundary and axial constraints to values of
!  the R,Z coordinates then calls AJAX_LOAD_RZLAM to fill in the profiles
!  W.A.Houlberg, P.I.Strand 8/2001
!Input:
!  r0                  -major radius of geometric center [m]
!  a0                  -minor radius in midplane [m]
!  s0                  -axis shift normalized to a0 [-]
!  e0                  -axis elongation normalized to a0 [-]
!  e1                  -edge elongation normalized to a0 [-]
!  d1                  -edge triangularity normalized to a0 [-]
!Output:
!  iflag               -error and warning flag [-]
!                      =-1 warning
!                      =0 no warnings or errors
!                      =1 error
!  message             -warning or error message (character)
!Optional input:
!  NRHO_AJAX           -no. of radial nodes in internal data [-]
!                      =21 default
!  NTHETA_AJAX         -no. of poloidal nodes in internal data [odd]
!                      =21 default
!  NKLAM_AJAX          -no. of lambda modes in nternal data [-]
!                      =5 default
!  RHOMAX_AJAX         -value of internal radial grid at R,Z boundary [rho]
!                      =1.0 default
!-------------------------------------------------------------------------------

!Declaration of input variables
      REAL(KIND=rspec), INTENT(IN) ::                                          &
     & r0,a0,s0,e0,e1,d1

!Declaration of output variables
      INTEGER, INTENT(OUT) ::                                                  &
     & iflag

      CHARACTER(len=*), INTENT(OUT) ::                                         &
     & message

!Declaration of optional input variables
      INTEGER, INTENT(IN), OPTIONAL ::                                         &
     & NRHO_AJAX,NTHETA_AJAX,NKLAM_AJAX

      REAL(KIND=rspec), INTENT(IN), OPTIONAL ::                                &
     & RHOMAX_AJAX

!Declaration of local variables
      INTEGER ::                                                               &
     & i,k,                                                                    &
     & nk_rz,nk_lam,                                                           &
     & nr_rz,nthetal

      INTEGER, ALLOCATABLE ::                                                  &
     & m(:),n(:)

      REAL(KIND=rspec) ::                                                      &
     & c,dm1,em1,                                                              &
     & r2,rhomaxl

      REAL(KIND=rspec), ALLOCATABLE ::                                         &
     & rho_rz(:),                                                              &
     & r(:,:),z(:,:)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Set radial grid normalization
      IF(PRESENT(RHOMAX_AJAX)) THEN

       !Use input value
        rhomaxl=RHOMAX_AJAX

      ELSE

        !Use default value
        rhomaxl=1.0_rspec

      ENDIF

!Set number of radial points
      IF(PRESENT(NRHO_AJAX)) THEN

        !Use input value
        nr_rz=NRHO_AJAX

      ELSE

        !Use default value
        nr_rz=31

      ENDIF

!Set number of poloidal modes
      IF(PRESENT(NKLAM_AJAX)) THEN

        !Use input value
        nk_lam=NKLAM_AJAX

      ELSE

        !Use default value
        nk_lam=8

      ENDIF

!Set number of poloidal points
      IF(PRESENT(NTHETA_AJAX)) THEN

        !Use input value
        nthetal=NTHETA_AJAX

      ELSE

        !Use default value
        nthetal=4*nk_lam+1

      ENDIF

!Allocate poloidal mode information
      ALLOCATE(m(nk_lam),                                                      &
     &         n(nk_lam))

!Set number of (R,Z) modes and mode numbers
      nk_rz=3
      n(:)=0
      m(:)=(/ (k-1,k=1,nk_lam) /)

!-------------------------------------------------------------------------------
!Convert from geometric to moments quantities     
!-------------------------------------------------------------------------------
!Moments triangularity ~geometric triangularity/4
      dm1=d1/4.0_rspec
      c=0.0_rspec

!Iterate for triangularity value, convergence is fast and robust
      DO i=1,10 !Over iteration

        c=4.0_rspec*dm1/(SQRT(1.0_rspec+32.0_rspec*dm1**2)+1.0_rspec)
        dm1=d1/(4.0_rspec-6.0_rspec*c**2)

      ENDDO !Over iteration

!Elongation
      em1=e1/(SQRT(1.0_rspec-c**2)*(1.0_rspec+2.0_rspec*dm1*c))

!-------------------------------------------------------------------------------
!Fill in radial values of R,Z expansion coefficients
!-------------------------------------------------------------------------------
!Allocate radial grid and R,Z arrays
      ALLOCATE(rho_rz(nr_rz),                                                  &
     &         r(nr_rz,nk_rz),                                                 &
     &         z(nr_rz,nk_rz))

!Initialization
      rho_rz(:)=0.0_rspec
      r(:,:)=0.0_rspec
      z(:,:)=0.0_rspec

!Set radial grid
      rho_rz(:)=(/ (REAL(i-1,rspec)/REAL(nr_rz-1,rspec),i=1,nr_rz) /)

!Radial variation
      DO i=1,nr_rz !Over radial nodes

        r2=rho_rz(i)**2
        r(i,1)=r0+a0*(s0*(1.0_rspec-r2)-dm1*r2)
        r(i,2)=a0*rho_rz(i)
        z(i,2)=-r(i,2)*(e0*(1.0_rspec-r2)+em1*r2)
        r(i,3)=a0*r2*dm1
        z(i,3)=r(i,3)*em1
 
      ENDDO !Over radial nodes

!Load private data
      CALL AJAX_LOAD_RZLAM(nr_rz,nk_rz,rho_rz,m,n,r,z,iflag,                   &
     &                     message,                                            &
     &                     NRHO_AJAX=nr_rz,                                    &
     &                     NTHETA_AJAX=nthetal,                                &
     &                     RHOMAX_AJAX=rhomaxl,                                &
     &                     NK_LAM=nk_lam)

      !Check messages
      IF(iflag /= 0) message='AJAX_LOAD_RZBDY/'//message

 9999 CONTINUE

!-------------------------------------------------------------------------------
!Cleanup
!-------------------------------------------------------------------------------

      IF(ALLOCATED(m)) THEN

        !Deallocate mode information
        DEALLOCATE(m,                                                          &
     &             n)

      ENDIF

      IF(ALLOCATED(rho_rz)) THEN

        !Deallocate grid and R,Z arrays
        DEALLOCATE(rho_rz,                                                     &
     &             r,                                                          &
     &             z)

      ENDIF

      END SUBROUTINE AJAX_LOAD_RZBDY

      SUBROUTINE AJAX_LOAD_RZLAM(nr_rz,nk_rz,rho_rz,m,n,r,z,iflag,             &
     &                           message,                                      &
     &                           K_GRID,L_MFILTER_AJAX,NRHO_AJAX,              &
     &                           NTHETA_AJAX,NZETA_AJAX,RHOMAX_AJAX,           &
     &                           NR_LAM,NK_LAM,RHO_LAM,LAM)
!-------------------------------------------------------------------------------
!AJAX_LOAD_RZLAM loads 2D/3D MHD equilibria in inverse coordinate form
!References:
!  W.A.Houlberg, P.I.Strand 11/2001
!Input:
!  nr_rz               -no. of radial nodes in the input R,Z [-]
!  nk_rz               -no. of poloidal & toroidal modes in the input R,Z [-]
!  rho_rz(nr_rz)       -radial nodes in the input R,Z [arb]
!  m(MAX(nk_rz,NK_LAM))-poloidal mode numbers [-]
!  n(MAX(nk_rz,NK_LAM))-toroidal mode numbers [-]
!  r(nr_rz,nk_rz)      -expansion coeffs for R [m]
!  z(nr_rz,nk_rz)      -expansion coeffs for Z [m]
!Output:
!  iflag               -error and warning flag [-]
!                      =-1 warning
!                      =0 none
!                      =1 error
!  message             -warning or error message [character]
!Optional input:
!  K_GRID              -option designating type of input radial grid [-]
!                      =0 proportional to sqrt(toroidal flux)
!                      =1 proportional to toroidal flux
!                      =else not allowed
!  L_MFILTER_AJAX      -option for mode filtering [logical]
!  NRHO_AJAX           -no. of radial nodes in internal data [-]
!                      =21 default
!  NTHETA_AJAX         -no. of poloidal nodes in internal data [odd]
!                      =21 default
!  NZETA_AJAX          -no. of toroidal nodes per period in internal data [odd]
!                      =11 default for non-axisymmetric plasma
!                      =0 default for axisymmetric plasma
!  RHOMAX_AJAX         -value of internal radial grid at R,Z boundary [rho]
!                      =1.0 default
!  NR_LAM              -no. of radial nodes in the input lambda [-]
!  NK_LAM              -no. of poloidal & toroidal modes in the input lambda [-]
!  RHO_LAM(NR_LAM)     -radial nodes in the input lambda [arb]
!  LAM(NR_LAM,NK_LAM)  -expansion coeffs for lambda [-]
!Comments:
!  Knowledge of the form of the input radial grid (see K_GRID) is necessary for
!    accurate mapping to the internal grid
!  The internal radial grid, rho_3d, is proportional to the square root of the
!    toroidal flux
!  The maximum value of the R,Z radial grid defines rho/rhomax=1
!  The internal radial grid will be scaled to rhomax if specified, otherwise it
!    will default to the domain [0,1]
!  The input poloidal and toroidal mode numbers are the same for the R,Z and
!    lambda expansions and dimensioned to the greater of the two; e.g., if there
!    are more poloidal and toroidal modes for the input lambda they should be
!    appended to the end of the sequence of values for the R,Z expansions
!  Mode filtering may be tried to improve calculations near the axis and the
!    outer boundary
!-------------------------------------------------------------------------------

!Declaration of input variables
      INTEGER, INTENT(IN) ::                                                   &
     & nr_rz,nk_rz,                                                            &
     & m(:),n(:)
     
      REAL(KIND=rspec), INTENT(IN) ::                                          & 
     & rho_rz(:),                                                              &
     & r(:,:),z(:,:)

!Declaration of output variables
      INTEGER, INTENT(OUT) ::                                                  &
     & iflag

      CHARACTER(len=*), INTENT(OUT) ::                                         &
     & message
  
!Declaration of optional input variables
      LOGICAL, INTENT(IN), OPTIONAL ::                                         &
     & L_MFILTER_AJAX

      INTEGER, INTENT(IN), OPTIONAL ::                                         &
     & K_GRID,                                                                 &
     & NRHO_AJAX,NTHETA_AJAX,NZETA_AJAX,                                       &
     & NR_LAM,NK_LAM

      REAL(KIND=rspec), INTENT(IN), OPTIONAL ::                                &
     & RHOMAX_AJAX,                                                            &
     & RHO_LAM(:),                                                             &
     & LAM(:,:)

!Local variables
      INTEGER ::                                                               &
     & i,j,k,                                                                  &
     & kmax,                                                                   &
     & nset1,                                                                  &
     & k_grid_l,                                                               &
     & k_vopt(3)=(/1,0,0/),                                                    &
     & k_bc1=3,k_bcn=0
  
      REAL(KIND=rspec), ALLOCATABLE ::                                         &
     & rho(:),                                                                 &
     & rmn(:,:),zmn(:,:),                                                      &
     & fspl(:,:),                                                              &
     & values(:,:)

!-------------------------------------------------------------------------------
!Check data for number of angular modes and periodicity
!-------------------------------------------------------------------------------
!Set number of angular modes
      krz_3d=nk_rz

      IF(PRESENT(NK_LAM)) THEN

        klam_3d=NK_LAM

      ELSE

        klam_3d=MAX(krz_3d,8)

      ENDIF

      kmax=MAX(krz_3d,klam_3d)

!Set the number of toroidal field periods, < 100 periods
      i=100

      DO k=1,kmax !Over modes

        !Find lowest non-zero toroidal mode number
        IF(n(k) /= 0) i=MIN(i,ABS(n(k)))

      ENDDO !Over modes

      IF(i == 100) THEN

        !Found no peiodicity in data, set to axisymmetric
        nper_3d=0

      ELSE

        !Non-axisymmetric 
        nper_3d=i

        !Check periodicity
        DO k=1,kmax !Over modes

          IF(n(k) /= 0) THEN

            j=MOD(n(k),nper_3d)

            IF(j /= 0) THEN

              !Some mode does not fit the periodicity
              iflag=1
              message='AJAX_LOAD_RZLAM/ERROR(1):periodicity not found'
              GOTO 9999

            ENDIF

          ENDIF

        ENDDO !Over modes

      ENDIF

!-------------------------------------------------------------------------------
!Check optional input for internal grid initialization
!-------------------------------------------------------------------------------
!Set the number of internal radial nodes
      IF(PRESENT(NRHO_AJAX)) THEN

        !Use input value
        nrho_3d=NRHO_AJAX

      ELSE

        !Use default value
        nrho_3d=31

      ENDIF

!Set the number of internal poloidal nodes
      IF(PRESENT(NTHETA_AJAX)) THEN

        !Use input value
        ntheta_3d=NTHETA_AJAX

        !Make sure the number of modes is odd for 4th order Simpson integration
        IF(MOD(ntheta_3d,2) == 0) THEN

          ntheta_3d=ntheta_3d+1
          iflag=-1
          message='AJAX_LOAD_RZLAM/WARNING(2):poloidal nodes reset'

        ENDIF

      ELSE

        !Use default value
        IF(nper_3d == 0) THEN
          !Axisymmetric
          !4k+1 typically yields <0.1% error in B_zeta = RB_t ~ const
          ntheta_3d=4*klam_3d+1

        ELSE

          !Non-axisymmetric
          !Error estimate yet to be quantified, typically |m| <= 8
          ntheta_3d=33

        ENDIF

      ENDIF

!Set the number of internal toroidal nodes
      IF(PRESENT(NZETA_AJAX)) THEN

        !Use input value
        nzeta_3d=NZETA_AJAX

        IF(nper_3d == 0 .AND. nzeta_3d /= 1) THEN

          nzeta_3d=1
          iflag=-1
          message='AJAX_LOAD_RZLAM/WARNING(3):toroidal nodes reset to 1'

        ENDIF

        !Make sure the number of modes is odd for 4th order Simpson integration
        IF(MOD(nzeta_3d,2) == 0) THEN

          nzeta_3d=nzeta_3d+1
          iflag=-1
          message='AJAX_LOAD_RZLAM/WARNING(4):toroidal nodes reset'

        ENDIF

      ELSEIF(nper_3d == 0) THEN

        !Use axisymmetric value
        nzeta_3d=1

      ELSE

        !Use non-axisymmetric default value
        !Error estimate yet to be quantified, typically |n/nper| <= 5
        nzeta_3d=21

      ENDIF

!Set the option for mode filtering
      IF(PRESENT(L_MFILTER_AJAX)) THEN

        !Use input value
        l_mfilter_3d=L_MFILTER_AJAX

      ELSE

        !Use default value
        l_mfilter_3d=.TRUE.

      ENDIF

!-------------------------------------------------------------------------------
!Set general information and initialize arrays
!-------------------------------------------------------------------------------
!Scale the outer radial boundary of R,Z to rhomax
      IF(PRESENT(RHOMAX_AJAX)) THEN

        !Use input value
        rhomax_3d=RHOMAX_AJAX

      ELSE

        !Use default value
        rhomax_3d=1.0_rspec

      ENDIF

!Set the mode numbers of the (m,n)=(0,0) and (m,n)=(1,0) modes
      km0n0_3d=0
      km1n0_3d=0

      DO k=1,kmax

        IF(ABS(m(k)) == 0 .AND. n(k) == 0) km0n0_3d=k
        IF(ABS(m(k)) == 1 .AND. n(k) == 0) km1n0_3d=k
        IF(km0n0_3d /= 0 .AND. km1n0_3d /= 0) EXIT

      ENDDO !Over modes

!Call initialization routine
      CALL AJAX_INIT

!-------------------------------------------------------------------------------
!Poloidal and toroidal mode data
!-------------------------------------------------------------------------------
!Initialization
      m_3d(:)=0
      n_3d(:)=0

!Copy data
      m_3d(1:kmax)=m(1:kmax)
      n_3d(1:kmax)=n(1:kmax)

!Set |m| with mode filtering for rho**|m| normalization
      DO k=1,kmax !Over modes

        mabs_3d(k)=ABS(m_3d(k))

        !Check for filtering
        IF(l_mfilter_3d) THEN

          IF(mabs_3d(k) > 3) THEN

            !Filtering for |m| > 3
            IF(MOD(mabs_3d(k),2) == 0) THEN

              !Even mode, remove rho**2
              mabs_3d(k)=2

            ELSE

              !Odd mode, remove rho**3
              mabs_3d(k)=3

            ENDIF

          ENDIF

        ENDIF

      ENDDO !Over modes

!-------------------------------------------------------------------------------
!Load R and Z expansion data
!-------------------------------------------------------------------------------
!Initialization
      r_3d(:,:,:)=0.0_rspec
      z_3d(:,:,:)=0.0_rspec

!Make copy of radial grid and expansion coefficients
!If the user does not provide an axial value for R,Z data, need to fill in
      IF(rho_rz(1) > rhores_3d) THEN

        !Allocate radial grid and R,Z arrays, add radial node at axis
        nset1=nr_rz+1
        ALLOCATE(rho(1:nset1),                                                 &
     &           rmn(1:nset1,1:nk_rz),                                         &
     &           zmn(1:nset1,1:nk_rz))
        rho(:)=0.0_rspec
        rmn(:,:)=0.0_rspec
        zmn(:,:)=0.0_rspec
        rho(2:nset1)=rho_rz(1:nr_rz)/rho_rz(nr_rz)
        rmn(2:nset1,1:nk_rz)=r(1:nr_rz,1:nk_rz)
        zmn(2:nset1,1:nk_rz)=z(1:nr_rz,1:nk_rz)

      ELSE

        !Allocate radial grid and R,Z arrays, use input grid
        nset1=nr_rz
        ALLOCATE(rho(1:nset1),                                                 &
     &           rmn(1:nset1,1:nk_rz),                                         &
     &           zmn(1:nset1,1:nk_rz))
        rho(:)=0.0_rspec
        rmn(:,:)=0.0_rspec
        zmn(:,:)=0.0_rspec
        rho(1:nset1)=rho_rz(1:nset1)/rho_rz(nset1)
        rmn(1:nset1,1:nk_rz)=r(1:nset1,1:nk_rz)
        zmn(1:nset1,1:nk_rz)=z(1:nset1,1:nk_rz)

      ENDIF

!Convert rho to sqrt(toroidal flux) if necesssary and scale
      k_grid_l=0
      IF(PRESENT(K_GRID)) k_grid_l=K_GRID

      IF(k_grid_l == 1) THEN

        !~toroidal flux
        rho(:)=SQRT(rho(:))*rhomax_3d

      ELSEIF(k_grid_l == 0) THEN

        !~sqrt(toroidal flux)
        rho(:)=rho(:)*rhomax_3d

      ELSE

        !Unallowed choice of input grid
        iflag=1
        message='AJAX_LOAD_RZLAM/ERROR(5):unallowed choice of K_GRID'
        GOTO 9999

      ENDIF

!Set the inner radial boundary of R,Z for checking axial extrapolation
      rhomin_3d=rho(2)

!Allocate spline arrays
      ALLOCATE(fspl(4,nset1),                                                  &
     &         values(3,nrho_3d))
      fspl(:,:)=0.0_rspec
      values(:,:)=0.0_rspec

!Normalize the expansion coefficients to rho**m
      DO k=1,krz_3d !Over modes

        DO i=2,nset1 !Over radial nodes

          rmn(i,k)=rmn(i,k)/rho(i)**mabs_3d(k)
          zmn(i,k)=zmn(i,k)/rho(i)**mabs_3d(k)

        ENDDO !Over radial nodes

!Parabolic extrapolation to axis (user values ignored)
        rmn(1,k)=(rmn(2,k)*rho(3)**2-rmn(3,k)*rho(2)**2)                       &
     &           /(rho(3)**2-rho(2)**2)
        zmn(1,k)=(zmn(2,k)*rho(3)**2-zmn(3,k)*rho(2)**2)                       &
     &           /(rho(3)**2-rho(2)**2)

!Map R_mn to internal radial grid
        fspl(1,1:nset1)=rmn(1:nset1,k)
        iflag=0
        message=''
        CALL SPLINE1_INTERP(k_vopt,nset1,rho,fspl,nrho_3d,rho_3d,values,       &
     &                      iflag,message,                                     &
     &                      K_BC1=k_bc1,                                       &
     &                      K_BCN=k_bcn)

        !Check messages
        IF(iflag > 0) THEN

          message='AJAX_LOAD_RZLAM(6)/'//message
          GOTO 9999

        ENDIF

        !Respline the R coeffs for internal storage
        r_3d(1,1:nrho_3d,k)=values(1,1:nrho_3d)
        CALL SPLINE1_FIT(nrho_3d,rho_3d,r_3d(:,:,k),                           &
     &                   K_BC1=k_bc1,                                          &
     &                   K_BCN=k_bcn)

!Map Z_mn to internal radial grid
        fspl(1,1:nset1)=zmn(1:nset1,k)            ! Copy input 
        iflag=0
        message=''
        CALL SPLINE1_INTERP(k_vopt,nset1,rho,fspl,nrho_3d,rho_3d,              &
     &                      values,iflag,message,                              &
     &                      K_BC1=k_bc1,                                       &
     &                      K_BCN=k_bcn)

        !Check messages
        IF(iflag > 0) THEN

          message='AJAX_LOAD_RZLAM(7)/'//message
          GOTO 9999

        ENDIF

        !Respline the Z coeffs for internal storage
        z_3d(1,1:nrho_3d,k)=values(1,1:nrho_3d)        
        CALL SPLINE1_FIT(nrho_3d,rho_3d,z_3d(:,:,k),                           &
     &                   K_BC1=k_bc1,                                          &
     &                   K_BCN=k_bcn)

      ENDDO !Over modes

!Set the length scale factors
      r000_3d=r_3d(1,1,km0n0_3d)

!-------------------------------------------------------------------------------
!Magnetic stream function
!-------------------------------------------------------------------------------
      IF(PRESENT(NR_LAM) .AND.                                                 &
     &   PRESENT(NK_LAM) .AND.                                                 &
     &   PRESENT(RHO_LAM) .AND.                                                &
     &   PRESENT(LAM)) THEN

        !Load lambda from input
        iflag=0
        message=''
        CALL AJAX_LOAD_LAMBDA(k_grid_l,NR_LAM,NK_LAM,rho_rz(nr_rz),            &
     &                        RHO_LAM,LAM,iflag,message)

        !Check messages
        IF(iflag /= 0) message='AJAX_LOAD_RZLAM(8)/'//message

      ELSEIF(nper_3d == 0) THEN

        !Calculate lambdas for axisymmetric plasma
        !Set up 3d grid for flux surface integrals
        iflag=0
        message=''
        CALL AJAX_INIT_FLUXAV_G(iflag,message)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_LOAD_RZLAM(9)/'//message
          IF(iflag > 0) GOTO 9999

        ENDIF


        !Calculate lambda
        CALL AJAX_INIT_LAMBDA

      ELSE

        !Need to specify stream function for non-axisymmetric plasmas
        iflag=1
        message='AJAX_LOAD_RZLAM/ERROR(10):need lambdas'

      ENDIF

 9999 CONTINUE

!-------------------------------------------------------------------------------
!Cleanup
!-------------------------------------------------------------------------------

      IF(ALLOCATED(rho)) THEN

        !Deallocate radial grid and R,Z arrays
        DEALLOCATE(rho,                                                        &
     &             rmn,                                                        &
     &             zmn)

      ENDIF

      IF(ALLOCATED(fspl)) THEN

        !Deallocate spline arrays
        DEALLOCATE(fspl,                                                       &
     &             values)

      ENDIF

      END SUBROUTINE AJAX_LOAD_RZLAM

      SUBROUTINE AJAX_LOAD_MAGFLUX(phitot,k_pflx,nr_pflx,rho_pflx,pflx,        &
     &                             iflag,message)
!-------------------------------------------------------------------------------
!AJAX_LOAD_MAGFLUX loads magnetic flux data
!References:
!  W.A.Houlberg, P.I.Strand 8/2001
!Input:
!  phitot              -total toroidal flux [Wb]
!  k_pflx              -option for represention of poloidal flux [-]
!                      =1 iotabar (rotational transform)
!                      =2 d(Psi)/d(rho)
!                      =else q (safety factor)
!  nr_pflx             -no. of radial nodes in input flux [-]
!  rho_pflx(nr_pflx)   -radial nodes in input flux [rho]
!  pflx(nr_pflx)       -poloidal magnetic flux function
!                      =q [-]
!                      =iotabar [-]
!                      =d(Psi)/d(rho) [Wb/rho]
!Output:
!  iflag               -error and warning flag [-]
!                      =-1 warning
!                      =0 none
!                      =1 error
!  message             -warning or error message [character]
!-------------------------------------------------------------------------------

!Declaration of input variables
      INTEGER, INTENT(IN) ::                                                   &
     & k_pflx,                                                                 &
     & nr_pflx

      REAL(KIND=rspec), INTENT(IN) ::                                          &
     & rho_pflx(:),pflx(:),                                                    &
     & phitot

!Declaration of output variables
      INTEGER, INTENT(OUT) ::                                                  &
     & iflag

      CHARACTER(len=*), INTENT(OUT) ::                                         &
     & message

!Declaration of local variables
      INTEGER, PARAMETER ::                                                    &
     & k_vopt(3)=(/1,0,0/),                                                    &
     & k_bc1=3,k_bcn=0 

      INTEGER ::                                                               &
     & nset1
 
      REAL(KIND=rspec) ::                                                      &
     & fspl(4,nr_pflx),                                                        &
     & values(3,nrho_3d)

!-------------------------------------------------------------------------------
!Iotabar allocation
!-------------------------------------------------------------------------------
      IF(ALLOCATED(iotabar_3d)) THEN

        !Radial dimension
        nset1=SIZE(iotabar_3d,2)

        !If storage requirements have changed, reallocate
        IF(nset1 /= nrho_3d) THEN

          !Reallocate iotabar
          DEALLOCATE(iotabar_3d)
          ALLOCATE(iotabar_3d(4,nrho_3d))
          iotabar_3d(:,:)=0.0_rspec

        ENDIF

      ELSE

        !Allocate iotabar
        ALLOCATE(iotabar_3d(4,nrho_3d))
        iotabar_3d(:,:)=0.0_rspec

      ENDIF

!-------------------------------------------------------------------------------
!Set total toroidal flux in private data
!-------------------------------------------------------------------------------
      phitot_3d=phitot

!-------------------------------------------------------------------------------
!Set iotabar in private data
!-------------------------------------------------------------------------------
!Initialize
      fspl(:,:)=0.0_rspec
      values(:,:)=0.0_rspec

      IF(k_pflx == 1) THEN

        !iotabar input
        fspl(1,1:nr_pflx)=pflx(1:nr_pflx)

      ELSEIF(k_pflx == 2) THEN

        !d(Psi)/d(rho) input
        fspl(1,1:nr_pflx)=pflx(1:nr_pflx)/rho_pflx(1:nr_pflx)
        fspl(1,1:nr_pflx)=0.5_rspec*rhomax_3d**2/phitot_3d                     &
     &                    *fspl(1,1:nr_pflx)

      ELSE

        !q input
        fspl(1,1:nr_pflx)=1.0_rspec/pflx(1:nr_pflx)

      ENDIF

!Interpolate iotabar onto internal grid
      iflag=0
      message=''
      CALL SPLINE1_INTERP(k_vopt,nr_pflx,rho_pflx,fspl,nrho_3d,rho_3d,         &
     &                    values,iflag,message,                                &
     &                    K_BC1=k_bc1,                                         &
     &                    K_BCN=k_bcn)

      !Check messages
      IF(iflag > 0) THEN

        message='AJAX_LOAD_MAGFLUX(2)/'//message
        GOTO 9999

      ENDIF

!Get spline coefficients on internal grid
      iotabar_3d(:,:)=0.0_rspec
      iotabar_3d(1,1:nrho_3d)=values(1,1:nrho_3d)

      !If value at origin not supplied, overwrite with approximation
      IF(rho_pflx(1) > rhores_3d) THEN

        iotabar_3d(1,1)=iotabar_3d(2,1)-rho_3d(2)                              &
     &                  *(iotabar_3d(3,1)-iotabar_3d(2,1))                     &
     &                  /(rho_3d(3)-rho_3d(2))

      ENDIF

      CALL SPLINE1_FIT(nrho_3d,rho_3d,iotabar_3d,                              &
     &                 K_BC1=k_bc1,                                            &
     &                 K_BCN=k_bcn)

!-------------------------------------------------------------------------------
!Set signs of poloidal and toroidal magnetic fields
!-------------------------------------------------------------------------------
      signbt_3d=SIGN(1.0_rspec,phitot_3d)
      signbp_3d=SIGN(1.0_rspec,iotabar_3d(1,nrho_3d))*signbt_3d

!-------------------------------------------------------------------------------
!Set logical switches for flux surface averaging
!-------------------------------------------------------------------------------
      l_fluxavb_3d=.FALSE.

 9999 CONTINUE

      END SUBROUTINE AJAX_LOAD_MAGFLUX

      SUBROUTINE AJAX_CAR2CYL(r_car,r_cyl)
!-------------------------------------------------------------------------------
!AJAX_CAR2CYL converts from Cartesian to cylindrical coordinates
!References:
!  W.A.Houlberg, P.I.Strand 8/2001
!Input:
!  r_car(3)            -Cartesian coordinates (x,y,z) [m,m,m]
!Output:
!  r_cyl(3)            -cylindrical coordinates (R,phi,Z) [m,rad,m]
!-------------------------------------------------------------------------------

!Declaration of input variables
      REAL(KIND=rspec), INTENT(IN) ::                                          &
     & r_car(3)

!Declaration of output variables
      REAL(KIND=rspec), INTENT(OUT) ::                                         &
     & r_cyl(3)

!-------------------------------------------------------------------------------
!Calculate cylindrical coordinates
!-------------------------------------------------------------------------------
!Initialize
      r_cyl(:)=0.0_rspec

!Cartesian to cylindrical conversion
      r_cyl(1)=SQRT(r_car(1)**2+r_car(2)**2) !R
      r_cyl(2)=ATAN2(r_car(2),r_car(1))      !phi
      r_cyl(3)=r_car(3)                      !Z

!Ensure 0 <= phi <= 2*pi
      IF(r_cyl(2) < 0.0_rspec) r_cyl(2)=r_cyl(2)+2.0_rspec*z_pi

      END SUBROUTINE AJAX_CAR2CYL

      SUBROUTINE AJAX_CYL2CAR(r_cyl,r_car)
!-------------------------------------------------------------------------------
!AJAX_CYL2CAR converts from cylindrical to Cartesian coordinates
!References:
!  W.A.Houlberg, P.I.Strand 8/2001
!Input:
!  r_cyl(3)            -cylindrical coordinates (R,phi,Z) [m,rad,m]
!Output:
!  r_car(3)            -Cartesian coordinates (x,y,z) [m,m,m]
!-------------------------------------------------------------------------------

!Declaration of input variables
      REAL(KIND=rspec), INTENT(IN) ::                                          &
     & r_cyl(3)

!Declaration of output variables
      REAL(KIND=rspec), INTENT(OUT) ::                                         &
     & r_car(3)

!-------------------------------------------------------------------------------
!Calculate the Cartesian coordinates
!-------------------------------------------------------------------------------
!Initialization
      r_car(:)=0.0_rspec

!Cylindrical to Cartesian conversion
      r_car(1)=r_cyl(1)*COS(r_cyl(2))   !x
      r_car(2)=r_cyl(1)*SIN(r_cyl(2))   !y
      r_car(3)=r_cyl(3)                 !z

      END SUBROUTINE AJAX_CYL2CAR

      SUBROUTINE AJAX_CYL2FLX(r_cyl,r_flx,iflag,message,                       &
     &                        G_CYL,GSQRT,TAU)
!-------------------------------------------------------------------------------
!AJAX_CYL2FLX converts from cylindrical to flux coordinates
!References:
!  S.E.Attenberger, W.A.Houlberg, S.P.Hirshman J Comp Phys 72 (1987) 435
!  W.A.Houlberg, P.I.Strand 8/2001
!Input:
!  r_cyl(3)            -cylindrical coordinates (R,phi,Z) [m,rad,m]
!Input/output:
!  r_flx(3)            -flux coordinates (rho,theta,zeta) [rho,rad,rad]
!Output:
!  iflag               -error and warning flag [-]
!                      =-1 warning
!                      =0 none
!                      =1 error
!  message             -warning or error message [character]
!Optional output:
!  G_CYL(6)            -R,Z derivatives
!                      =(R_rho,R_theta,R_zeta,Z_rho,Z_theta,Z_zeta)
!                       [m/rho,m,      m,     m/rho,m,      m     ]
!  GSQRT               -3D Jacobian [m**3/rho]
!  TAU                 -2D Jacobian in phi=zeta=constant plane [m**2/rho]
!Comments:
!  The basic method is a Newton iteration in 2 dimensions
!  Input values of r_flx are used as an initial guess
!  Convergence is typically one order of magnitude reduction in the
!    error in R and Z per iteration - should rarely exceed 5 iterations
!  The toroidal angle in flux coordinates is the same as in cylindrical
!    coordinates giving a right-handed system with positive Jacobian
!  Point 0 is the previous best guess and point 1 is a trial point
!  If point 1 is further from (R,Z) than point 0, a new trial point is
!    generated by halving the step and using interpolated derivatives
!-------------------------------------------------------------------------------

!Declaration of input variables
      REAL(KIND=rspec), INTENT(IN) ::                                          &
     & r_cyl(3)

!Declaration of input/output variables
      REAL(KIND=rspec), INTENT(INOUT) ::                                       &
     & r_flx(3)
      
!Declaration of output variables
      INTEGER, INTENT(OUT) ::                                                  &
     & iflag

      CHARACTER(len=*), INTENT(OUT) ::                                         &
     & message

!Declaration of optional output variables
      REAL(KIND=rspec), INTENT(OUT), OPTIONAL ::                               & 
     & G_CYL(6),                                                               &
     & GSQRT,TAU
     
!Declaration of local variables
      INTEGER ::                                                               &
     & it,itmax,                                                               &
     & nh
 
      REAL(KIND=rspec) ::                                                      &
     & dr,dr0,dr1,dz,dz0,dz1,                                                  &
     & drho,dtheta,                                                            &
     & err0,err1,                                                              &            
     & gsqrt1,tau0,tau1,taut,                                                  &                          
     & tol

      REAL(KIND=rspec) ::                                                      &
     & r_flx0(3),g_cyl0(6),                                                    &
     & r_flx1(3),r_cyl1(3),g_cyl1(6),                                          &
     & g_cylt(6)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
      r_flx0(:)=0.0_rspec
      g_cyl0(:)=0.0_rspec
      r_flx1(:)=0.0_rspec
      r_cyl1(:)=0.0_rspec
      g_cyl1(:)=0.0_rspec
      g_cylt(:)=0.0_rspec

!Tolerance for convergence
      tol=0.1_rspec*rhores_3d/rhomax_3d

!Maximum iterations
      itmax=20

!Grid halving index
      nh=1

!Set flux coordinate toroidal angle equal to cylindrical angle
      r_flx(3)=r_cyl(2)                      ! phi

!Set starting points and parameters
      tau0=0.0_rspec
      err0=0.1_rspec*z_large
      r_flx0(1:3)=r_flx(1:3)
      r_flx1(1:3)=r_flx(1:3)

!Move rho away from axis if necessary
      IF(r_flx1(1) < rhores_3d) r_flx1(1)=rhores_3d

!-------------------------------------------------------------------------------
!Iterate to find flux coordinates
!-------------------------------------------------------------------------------
      DO it=1,itmax !Over iteration

!Get cylindrical coordinates at point 1
        iflag=0
        message=''
        CALL AJAX_FLX2CYL(r_flx1,r_cyl1,iflag,message,                         &
     &                    G_CYL=g_cyl1,                                        &
     &                    GSQRT=gsqrt1,                                        &
     &                    TAU=tau1)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_CYL2FLX(1)/'//message
          IF(iflag > 0) GOTO 9999

        ENDIF

        dr1=(r_cyl(1)-r_cyl1(1))
        dz1=(r_cyl(3)-r_cyl1(3))
        err1=(dr1**2+dz1**2)/r000_3d**2

!Check convergence
        IF((ABS(dr1) <= tol*r000_3d .AND. ABS(dz1) <= tol*r000_3d)             &
     &     .OR. (ABS(err1-err0) < 5.0e-6_rspec)) THEN

          !Converged, but first check if solution is in the R,Z domain
          IF(r_flx1(1) > rhomax_3d) THEN

            iflag=-1
            message='AJAX_CYL2FLX(2)/WARNING:outside R,Z domain'

          ENDIF

          !Solution is at point 1, copy to output arrays and exit
          r_flx(1:3)=r_flx1(1:3)
          IF(r_flx(2) < 0.0_rspec) r_flx(2)=r_flx(2)+2.0_rspec*z_pi
          r_flx(2)=MOD(r_flx(2),2.0_rspec*z_pi)
          IF(PRESENT(GSQRT)) GSQRT=gsqrt1
          IF(PRESENT(G_CYL)) G_CYL(:)=g_cyl1(:)
          IF(PRESENT(TAU)) TAU=tau1
          GOTO 9999

        ENDIF

!Need improved estimate for rho and theta and get new Jacobian
        IF(r_flx1(1) < rhores_3d) r_flx1(1)=rhores_3d

        IF(ABS(tau1) < 10.0_rspec*z_small) THEN

          !No solution exists for the Newton method
          !Assume the solution is at the origin where tau=0
          iflag=-1
          message='AJAX_CYL2FLX(3)/WARNING:solution at origin'
          r_flx(1)=0.0_rspec
          r_flx(2)=0.0_rspec

        ENDIF

!Check consistency of sign of Jacobian
        IF(SIGN(1.0_rspec,tau1) /= SIGN(1.0_rspec,tau0) .AND.                  &
     &          tau0 /= 0.0_rspec) THEN

          !Bad data
          iflag=1
          message='AJAX_CYL2FLX(4)/ERROR:bad Jacobian'
          GOTO 9999

        ENDIF

!Check whether convergence is improving
        IF(err1 < err0) THEN

          !Converging, double step if possible and set point 0 = point 1
          IF(nh > 1) nh=nh/2
          r_flx0(:)=r_flx1(:)
          g_cyl0(:)=g_cyl1(:)
          err0=err1
          tau0=tau1
          dr0=dr1
          dz0=dz1

        ELSE

          !Not converging - halve step
          nh=nh*2

        ENDIF

!Calculate new rho and theta steps
        !Use empirical combination of last two steps to avoid oscillation
        !0 and 1 may coincide here
        g_cylt(:)=0.75_rspec*g_cyl0(:)+0.25_rspec*g_cyl1(:)
        taut=0.75_rspec*tau0+0.25_rspec*tau1
        dr=dr0/REAL(nh,rspec)
        dz=dz0/REAL(nh,rspec)
        !Project to desired point using Jacobian information
        !drho=(R_theta*dZ-Z_theta*dR)/tau
        drho=(g_cylt(2)*dz-g_cylt(5)*dr)/taut
        !dtheta=(Z_rho*dR-R_rho*dZ)/tau
        dtheta=(g_cylt(4)*dr-g_cylt(1)*dz)/taut
        IF(ABS(dtheta) > 0.25_rspec*z_pi)                                      &
     &    dtheta=SIGN(0.25_rspec*z_pi,dtheta)

!Set flux coordinates for new point 1
        r_flx1(1)=r_flx0(1)+drho
        r_flx1(2)=r_flx0(2)+dtheta

        IF(r_flx1(1) < 0.0_rspec) THEN

          !Stepped past minor axis, flip flux coordinates
          r_flx1(1)=-r_flx1(1)
          r_flx1(2)=r_flx1(2)+z_pi-2.0_rspec*dtheta

        ENDIF

      ENDDO !Over iteration

!Exceeded maximum iterations
      iflag=1
      message='AJAX_CYL2FLX(5)/ERROR:max iterations exceeded'

 9999 CONTINUE

      END SUBROUTINE AJAX_CYL2FLX

      SUBROUTINE AJAX_FLX2CYL(r_flx,r_cyl,iflag,message,                       &
     &                        G_CYL,GSQRT,TAU)
!-------------------------------------------------------------------------------
!AJAX_FLX2CYL converts from flux to cylindrical coordinates
!References:
!  S.E.Attenberger, W.A.Houlberg, S.P.Hirshman J Comp Phys 72 (1987) 435
!  W.A.Houlberg, P.I.Strand 8/2001
!Input:
!  r_flx(3)            -flux coordinates (rho,theta,zeta) [rho,rad,rad]
!Output:
!  r_cyl(3)            -cylindrical coordinates (R,phi,Z) [m,rad,m]
!  iflag               -error and warning flag [-]
!                      =-1 warning
!                      =0 none
!                      =1 error
!  message             -warning or error message [character]
!Optional output:
!  G_CYL(6)            -R,Z derivatives
!                      =(R_rho,R_theta,R_zeta,Z_rho,Z_theta,Z_zeta)
!                       [m/rho,m,      m,     m/rho,m,      m     ]
!  GSQRT               -3D Jacobian [m**3/rho]
!  TAU                 -2D Jacobian in phi=zeta=constant plane [m**2/rho]
!Comments:
!  The representation is extended beyond the rho_3d grid by linear
!    extrapolation of the m=1, n=0 term in rho, with all other terms
!    held fixed at the edge value
!  This permits unique representation of all space for flux surfaces
!    that are everywhere convex (does not include bean-shapes)
!  rho**m is factored out of the Fourier coefficients prior to spline
!    fitting for increased accuracy near the origin
!  Mode filtering is used to improve calculations near the axis and the outer
!    boundary
!-------------------------------------------------------------------------------

!Declaration of input variables
      REAL(KIND=rspec), INTENT(IN) ::                                          &
     & r_flx(3)

!Declaration of output variables
      INTEGER, INTENT(OUT) ::                                                  &
     & iflag

      CHARACTER(len=*), INTENT(OUT) ::                                         &
     & message

      REAL(KIND=rspec), INTENT(OUT) ::                                         &
     & r_cyl(3)

!Declaration of optional output variables
      REAL(KIND=rspec), INTENT(OUT), OPTIONAL ::                               &
     & G_CYL(6),                                                               &
     & GSQRT,TAU
    
!Declaration of local variables
      INTEGER ::                                                               &
     & k

      INTEGER, SAVE ::                                                         &
     & i=1

      INTEGER, PARAMETER ::                                                    &
     & k_vopt(3)=(/1,1,0/)

      REAL(KIND=rspec) ::                                                      &
     & ct,st,                                                                  &
     & rho,rhom,drhom,                                                         &
     & rmnx,drmnx,zmnx,dzmnx,                                                  &
     & g_cylt(6),                                                              &
     & value(3)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Error messages
      iflag=0
      message=''

!Cylindrical coordinates and metrics
      r_cyl(:)=0.0_rspec
      g_cylt(:)=0.0_rspec

!Spline values
      value(:)=0.0_rspec

!phi = zeta
      r_cyl(2)=r_flx(3)

!Limit to R,Z domain
      rho=MIN(r_flx(1),rhomax_3d)

!Axial resolution
      rho=MAX(rho,rhores_3d)

!-------------------------------------------------------------------------------
!Evaluate R,Z and derivatives
!-------------------------------------------------------------------------------
!Loop over modes
      DO k=1,krz_3d !Over modes

!Set sine and cosine values
        ct=COS(m_3d(k)*r_flx(2)-n_3d(k)*r_flx(3))
        st=SIN(m_3d(k)*r_flx(2)-n_3d(k)*r_flx(3))

!Calculate rho**m and its derivative
        IF(mabs_3d(k) == 0) THEN

          !rho**0
          rhom=1.0_rspec
          drhom=0.0_rspec

        ELSEIF(r_flx(1) < rhomax_3d+rhores_3d) THEN

          !Inside last closed surface, rho**|m|
          rhom=rho**mabs_3d(k)
          drhom=mabs_3d(k)*rho**(mabs_3d(k)-1)

        ELSE

          !Outside last closed surface, rhomax**|m|
          rhom=rhomax_3d**mabs_3d(k)
          drhom=0.0_rspec

        ENDIF

!R_mn and dR_mn/drho from spline fits
        CALL SPLINE1_EVAL(k_vopt,nrho_3d,rho,rho_3d,                           &
     &                    r_3d(1:4,1:nrho_3d,k),i,value)
        rmnx=value(1)
        drmnx=value(2)

!Z_mn and dZ_mn/drho from spline fits
        CALL SPLINE1_EVAL(k_vopt,nrho_3d,rho,rho_3d,                           &
     &                    z_3d(1:4,1:nrho_3d,k),i,value)
        zmnx=value(1)
        dzmnx=value(2)

!R = sum_mn [ R_mn * cos(m*theta-n*zeta) * rho**m ]
        r_cyl(1)=r_cyl(1)+rmnx*ct*rhom

!Z = sum_mn [ Z_mn * sin(m*theta-n*zeta) * rho**m ]
        r_cyl(3)=r_cyl(3)+zmnx*st*rhom

!dR/dtheta = sum_mn [ -m * R_mn * sin(m*theta-n*zeta) * rho**m ]
        g_cylt(2)=g_cylt(2)-m_3d(k)*rmnx*st*rhom

!dR/dzeta = sum_mn [ n * R_mn * sin(m*theta-n*zeta) * rho**m ]
        g_cylt(3)=g_cylt(3)+n_3d(k)*rmnx*st*rhom

!dZ/dtheta = sum_mn [ m * Z_mn * cos(m*theta-n*zeta) * rho**m ]
        g_cylt(5)=g_cylt(5)+m_3d(k)*zmnx*ct*rhom

!dZ/dtheta = sum_mn [ -n * Z_mn * cos(m*theta-n*zeta) * rho**m ]
        g_cylt(6)=g_cylt(6)-n_3d(k)*zmnx*ct*rhom

!Radial derivatives inside R,Z domain
        IF(r_flx(1) <= rhomax_3d+rhores_3d) THEN

          !dR/drho = sum_mn [ rho**m * dR_mn/drho + R_mn *d(rho**m)/drho ]
          !                   * cos(m*theta-n*zeta)
          g_cylt(1)=g_cylt(1)+(drmnx*rhom+rmnx*drhom)*ct

          !dZ/drho = sum_mn [ rho**m * dZ_mn/drho + Z_mn *d(rho**m)/drho ]
          !                   * sin(m*theta-n*zeta)
          g_cylt(4)=g_cylt(4)+(dzmnx*rhom+zmnx*drhom)*st

        ENDIF

      ENDDO !Over modes

!Radial derivatives outside R,Z domain
      IF(r_flx(1) > rhomax_3d+rhores_3d) THEN

        !This point is off the grid toward the wall
        ct=COS(m_3d(km1n0_3d)*r_flx(2))
        st=SIN(m_3d(km1n0_3d)*r_flx(2))
        r_cyl(1)=r_cyl(1)                                                      &
     &           +(r_flx(1)-rhomax_3d)*r_3d(1,nrho_3d,km1n0_3d)*ct
        r_cyl(3)=r_cyl(3)                                                      &
     &           +(r_flx(1)-rhomax_3d)*z_3d(1,nrho_3d,km1n0_3d)*st
        g_cylt(1)=r_3d(1,nrho_3d,km1n0_3d)*ct
        g_cylt(2)=g_cylt(2)                                                    &
     &            -(r_flx(1)-rhomax_3d)*r_3d(1,nrho_3d,km1n0_3d)*st            &
     &            *m_3d(km1n0_3d)
        g_cylt(4)=z_3d(1,nrho_3d,km1n0_3d)*st
        g_cylt(5)=g_cylt(5)                                                    &
     &            +(r_flx(1)-rhomax_3d)*z_3d(1,nrho_3d,km1n0_3d)*ct            &
     &            *m_3d(km1n0_3d)

      ENDIF

!-------------------------------------------------------------------------------
!Optional output
!-------------------------------------------------------------------------------
!R,Z derivatives
      IF(PRESENT(G_CYL)) G_CYL(:)=g_cylt(:)

!2D Jacobian = dR/dtheta * dZ/drho - dR/drho * dZ/dtheta
      IF(PRESENT(TAU)) TAU=g_cylt(2)*g_cylt(4)-g_cylt(1)*g_cylt(5)

!3D Jacobian = R * tau
      IF(PRESENT(GSQRT))                                                       &
     &  GSQRT=r_cyl(1)*(g_cylt(2)*g_cylt(4)-g_cylt(1)*g_cylt(5))

      END SUBROUTINE AJAX_FLX2CYL

      SUBROUTINE AJAX_B(r_flx,r_cyl,g_cyl,iflag,message,                       &
     &                  B_CON,B_CO,B_CYL,B_CAR,B_POL,B_TOR,B_MOD)
!-------------------------------------------------------------------------------
!AJAX_B gets components of B in various forms
!References:
!  S.E.Attenberger, W.A.Houlberg, S.P.Hirshman, J Comp Phys 72 (1987) 435
!  W.A.Houlberg, P.I.Strand 8/2001
!Input:
!  r_flx(3)            -flux coordinates (rho,theta,zeta) [rho,rad,rad]
!  r_cyl(3)            -cylindrical coordinates (R,phi,Z) [m,rad,m]
!  g_cyl(6)            -R,Z derivatives
!                      =(R_rho,R_theta,R_zeta,Z_rho,Z_theta,Z_zeta)
!                       [m/rho,m,      m,     m/rho,m,      m     ]
!Output:
!  iflag               -error and warning flag [-]
!                      =-1 warning
!                      =0 none
!                      =1 error
!  message             -warning or error message [character]
!Optional output:
!  B_CON(3)            -contravariant (B^rho,B^theta,B^zeta) [T*rho/m,T/m,T/m]
!  B_CO(3)             -covariant (B_rho,B_theta,B_zeta) [T*m/rho,T*m,T*m]
!  B_CYL(3)            -cylindrical (B_R,B_phi,B_Z) [T,T,T]
!  B_CAR(3)            -Cartesian (B_x,B_y,B_z) [T,T,T]
!  B_POL               -poloidal field [T]
!  B_TOR               -toroidal field [T]
!  B_MOD               -|B| [T]
!-------------------------------------------------------------------------------

!Declaration of input variables
      REAL(KIND=rspec), INTENT(IN) ::                                          &
     & r_flx(3),r_cyl(3),                                                      &
     & g_cyl(6)

!Declaration of output variables
      INTEGER, INTENT(OUT) ::                                                  &
     & iflag

      CHARACTER(len=*), INTENT(OUT) ::                                         &
     & message

!Declaration of optional output variables
      REAL(KIND=rspec), INTENT(OUT), OPTIONAL ::                               &
     & B_CON(3),B_CO(3),B_CYL(3),B_CAR(3),                                     &
     & B_POL,B_TOR,B_MOD

!Declaration of local variables
      REAL(KIND=rspec) ::                                                      &
     & gsqrt,                                                                  &
     & iotabar(1),phiprm(1),                                                   &
     & lam_theta,lam_zeta,                                                     &
     & b_cont(3),b_cylt(3)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!3D Jacobian = R * (dR/dtheta * dZ/drho - dR/drho * dZ/dtheta)
      gsqrt=r_cyl(1)*(g_cyl(2)*g_cyl(4)-g_cyl(1)*g_cyl(5))

!-------------------------------------------------------------------------------
!Get magnetic stream function derivatives
!-------------------------------------------------------------------------------
      CALL AJAX_LAMBDA(r_flx,lam_theta,lam_zeta)

!-------------------------------------------------------------------------------
!Get rotational transform and derivative of toroidal flux
!-------------------------------------------------------------------------------
      CALL AJAX_MAGFLUX(1,r_flx,iflag,message,                                 &
     &                  IOTABAR_R=iotabar,                                     &
     &                  PHIPRM_R=phiprm)

      !Check messages
      IF(iflag /= 0) THEN

        message='AJAX_B/'//message
        IF (iflag > 0 ) GOTO 9999

      ENDIF

!-------------------------------------------------------------------------------
!Calculate the contravariant field components
!-------------------------------------------------------------------------------
!B^rho
      b_cont(1)=0.0_rspec

!B^theta
      b_cont(2)=phiprm(1)*(iotabar(1)-lam_zeta)/(2.0_rspec*z_pi*gsqrt)

!B^zeta
      b_cont(3)=phiprm(1)*(1.0_rspec+lam_theta)/(2.0_rspec*z_pi*gsqrt)

      IF(PRESENT(B_CON)) B_CON(1:3)=b_cont(1:3)

!-------------------------------------------------------------------------------
!Calculate the cylindrical field components
!-------------------------------------------------------------------------------
      IF(PRESENT(B_CYL) .OR.                                                   &
     &   PRESENT(B_CO) .OR.                                                    &
     &   PRESENT(B_CAR) .OR.                                                   &
     &   PRESENT(B_POL) .OR.                                                   &
     &   PRESENT(B_TOR) .OR.                                                   &
     &   PRESENT(B_MOD)) THEN

!B_R
        b_cylt(1)=g_cyl(2)*b_cont(2)+g_cyl(3)*b_cont(3)

!B_phi
        b_cylt(2)=r_cyl(1)*b_cont(3)

!B_Z
        b_cylt(3)=g_cyl(5)*b_cont(2)+g_cyl(6)*b_cont(3)

        IF(PRESENT(B_CYL)) B_CYL(1:3)=b_cylt(1:3)

!-------------------------------------------------------------------------------
!Calculate the covariant field components
!-------------------------------------------------------------------------------
        IF(PRESENT(B_CO)) THEN

!B_rho
          B_CO(1)=b_cylt(1)*g_cyl(1)+b_cylt(3)*g_cyl(4)

!B_theta
          B_CO(2)=b_cylt(1)*g_cyl(2)+b_cylt(3)*g_cyl(5)

!B_zeta
          B_CO(3)=b_cylt(1)*g_cyl(3)+b_cylt(2)*r_cyl(1)                        &
     &            +b_cylt(3)*g_cyl(6)

        ENDIF

!-------------------------------------------------------------------------------
!Calculate the Cartesian field components
!-------------------------------------------------------------------------------
        IF(PRESENT(B_CAR)) THEN

!B_x
          B_CAR(1)=b_cylt(1)*COS(r_cyl(2))-b_cylt(2)*SIN(r_cyl(2))

!B_y
          B_CAR(2)=b_cylt(1)*SIN(r_cyl(2))+b_cylt(2)*COS(r_cyl(2))

!B_z
          B_CAR(3)=b_cylt(3)

        ENDIF

!-------------------------------------------------------------------------------
!Calculate the poloidal field
!-------------------------------------------------------------------------------
        IF(PRESENT(B_POL)) THEN

!B_pol
          B_POL=signbp_3d*SQRT(b_cylt(1)**2+b_cylt(3)**2)

        ENDIF

!-------------------------------------------------------------------------------
!Calculate the toroidal field
!-------------------------------------------------------------------------------
        IF(PRESENT(B_TOR)) THEN

!B_tor
          B_TOR=b_cylt(2)

        ENDIF

!-------------------------------------------------------------------------------
!Calculate |B|
!-------------------------------------------------------------------------------
        IF(PRESENT(B_MOD)) THEN

!B_mod
          B_MOD=SQRT(b_cylt(1)**2+b_cylt(2)**2+b_cylt(3)**2)

        ENDIF

      ENDIF

 9999 CONTINUE

      END SUBROUTINE AJAX_B

      SUBROUTINE AJAX_LAMBDA(r_flx,lam_theta,lam_zeta)
!-------------------------------------------------------------------------------
!AJAX_LAMBDA gets the stream function and its derivatives
!References:
!  S.E.Attenberger, W.A.Houlberg, S.P.Hirshman, J Comp Phys 72 (1987) 435
!  W.A.Houlberg, P.I.Strand 8/2001
!Input:
!  r_flx(3)-flux coordinates (rho,theta,zeta) [rho,rad,rad]
!Output:
!  lam_theta           -d(lambda)/d(theta) [-]
!  lam_zeta            -d(lambda)/d(zeta) [-]
!Comments:
!  Mode filtering is used to improve calculations near the axis and the outer
!    boundary
!-------------------------------------------------------------------------------

!Declaration of input variables
      REAL(KIND=rspec), INTENT(IN) ::                                          &
     & r_flx(3)

!Declaration of output variables
      REAL(KIND=rspec), INTENT(OUT) ::                                         &
     & lam_theta,lam_zeta
     
!Declaration of local variables
      INTEGER, PARAMETER ::                                                    &
     & k_vopt(3)=(/1,0,0/)

      INTEGER, SAVE ::                                                         &
     & i=1

      INTEGER  ::                                                              &
     & k

      REAL(KIND=rspec) ::                                                      &
     & value(3),                                                               &
     & cosk,                                                                   &
     & lam,                                                                    &
     & rho

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Magnetic stream function derivatives
      lam_theta=0.0_rspec
      lam_zeta=0.0_rspec

!Spline values
      value(:)=0.0_rspec

!Limit to R,Z domain
      rho=MIN(r_flx(1),rhomax_3d)

!Axial resolution
      rho=MAX(rho,rhores_3d)

!-------------------------------------------------------------------------------
!Calculate derivatives of the magnetic stream function
!-------------------------------------------------------------------------------
!Loop over the modes     
      DO k=1,klam_3d !Over modes

        cosk=COS(m_3d(k)*r_flx(2)-n_3d(k)*r_flx(3))
        CALL SPLINE1_EVAL(k_vopt,nrho_3d,rho,rho_3d,                           &
     &                    lam_3d(1:4,1:nrho_3d,k),i,value)
        !Reintroduce normalization
        lam=value(1)*rho**mabs_3d(k)
        lam_theta=lam_theta+lam*m_3d(k)*cosk
        lam_zeta=lam_zeta-lam*n_3d(k)*cosk

      ENDDO !Over modes

      END SUBROUTINE AJAX_LAMBDA

      SUBROUTINE AJAX_FLUXAV_B(nrho_r,rho_r,iflag,message,                     &
     &                         B2_R,BM2_R,FM_R,FTRAP_R,GR2BM2_R,GRTH_R,        &
     &                         SUS11_R,SUS12_R,SUS21_R,SUS22_R)
!-------------------------------------------------------------------------------
!AJAX_FLUXAV_B gets flux surface quantities that depend on B
!References:
!  W.A.Houlberg, P.I.Strand 11/2001
!Input:
!  nrho_r              -no. of radial nodes [-]
!  rho_r(nrho_r)       -radial nodes [rho]
!Output:
!  iflag               -error and warning flag [-]
!                      =-1 warning
!                      =0 no warnings or errors
!                      =1 error
!  message             -warning or error message [character]
!Optional output:
!  B2_R(nrho_r)        -<B**2> [T**2]
!  BM2_R(nrho_r)       -<1/B**2> [/T**2]
!  FM_R(m,nrho_r)      -poloidal moments of <[n.grad(B)]**2>/<B**2> [-]
!  FTRAP_R(nrho_r)     -trapped particle fraction [-]
!  GR2BM2_R(nrho_r)    -<grad(rho)**2/B**2> [rho**2/m**2/T**2]
!  GRTH_R(nrho_r)      -n.grad(Theta) [/m]
!  SUS11_R(nrho_r)     -suceptance matrix element 11 (pol-pol) [-]
!  SUS12_R(nrho_r)     -suceptance matrix element 11 (pol-tor) [-]
!  SUS11_R(nrho_r)     -suceptance matrix element 11 (tor-pol) [-]
!  SUS22_R(nrho_r)     -suceptance matrix element 11 (tor-tor) [-]
!-------------------------------------------------------------------------------

!Declaration of input variables
      INTEGER, INTENT(IN) ::                                                   &
     & nrho_r

      REAL(KIND=rspec), INTENT(IN) ::                                          &
     & rho_r(:)

!Declaration of output variables
      INTEGER, INTENT(OUT) ::                                                  &
     & iflag

      CHARACTER(len=*), INTENT(OUT) ::                                         &
     & message

!Declaration of optional output variables
      REAL(KIND=rspec), INTENT(OUT), OPTIONAL ::                               &
     & B2_R(:),BM2_R(:),FM_R(:,:),FTRAP_R(:),GR2BM2_R(:),GRTH_R(:),            &
     & SUS11_R(:),SUS12_R(:),SUS21_R(:),SUS22_R(:)

!Declaration of local variables
      INTEGER ::                                                               &
     & i,j,k,m,                                                                &
     & nr

      REAL(KIND=rspec) ::                                                      &
     & dbdt,                                                                   &
     & h

      REAL(KIND=rspec) ::                                                      &
     & b(nrho_3d),b2(nrho_3d),bmax(nrho_3d),                                   &
     & captheta(nrho_3d,ntheta_3d),                                            &
     & gam(nrho_3d),v1(nrho_3d),v2(nrho_3d)

      REAL(KIND=rspec) ::                                                      &
     & v_r(nrho_r)

!-------------------------------------------------------------------------------
!Check whether flux surface averaging arrays have been set
!-------------------------------------------------------------------------------
      IF(.NOT. l_fluxavg_3d) THEN

        iflag=0
        message=''
        CALL AJAX_INIT_FLUXAV_G(iflag,message)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_FLUXAV_B(1)/'//message
          GOTO 9999

        ENDIF


      ENDIF

      IF(.NOT. l_fluxavb_3d) CALL AJAX_INIT_FLUXAV_B

!-------------------------------------------------------------------------------
!Check whether values outside R,Z domain are requested (nr is last point inside)
!-------------------------------------------------------------------------------
      DO i=nrho_r,1,-1

        nr=i
        IF(rho_r(nr) < rhomax_3d+rhores_3d) EXIT

      ENDDO

!-------------------------------------------------------------------------------
!<B**2> on internal grid for b2_r, ftrap_r, and fm_r
!-------------------------------------------------------------------------------
      IF(PRESENT(B2_R) .OR.                                                    &
     &   PRESENT(FM_R) .OR.                                                    &
     &   PRESENT(FTRAP_R)) THEN

!Initialization
        b2(:)=0.0_rspec
        gt_3d(:)=0.0_rspec
        az_3d(:)=0.0_rspec

!Calculate on internal grid 
        DO i=2,nrho_3d !Over radial nodes

          DO k=1,nzeta_3d !Over toroidal nodes

            gt_3d(1:ntheta_3d)=gsqrt_3d(i,1:ntheta_3d,k)                       &
     &                         *b_3d(i,1:ntheta_3d,k)**2
            az_3d(k)=SUM(wtheta_3d*gt_3d) !Over poloidal nodes

          ENDDO !Over toroidal nodes

          b2(i)=SUM(wzeta_3d*az_3d)/vp_3d(i) !Over toroidal nodes

        ENDDO !Over radial nodes

!Extrapolate to axis
        b2(1)=b2(2)-rho_3d(2)*(b2(3)-b2(2))/(rho_3d(3)-rho_3d(2))

      ENDIF

!-------------------------------------------------------------------------------
!Theta and gam=n.grad(Theta) for Fm and n.grad(Theta) in axisymmetric plasmas
!-------------------------------------------------------------------------------
      IF((PRESENT(FM_R) .OR.                                                   &
     &    PRESENT(GRTH_R)) .AND.                                               &
     &   nper_3d == 0) THEN

!Initialization
        captheta(:,:)=0.0_rspec
        gam(:)=0.0_rspec

!Calculate on internal grid 
        DO i=2,nrho_3d !Over radial nodes

          DO j=2,ntheta_3d !Over poloidal nodes

            captheta(i,j)=captheta(i,j-1)+0.5_rspec*dtheta_3d                  &
     &                    *(b_3d(i,j-1,1)/btheta_3d(i,j-1,1)                   &
     &                     +b_3d(i,j,1)/btheta_3d(i,j,1))

          ENDDO !Over poloidal nodes

          gam(i)=2.0_rspec*z_pi/captheta(i,ntheta_3d)
          captheta(i,1:ntheta_3d)=gam(i)*captheta(i,1:ntheta_3d)

        ENDDO !Over radial nodes

!Extrapolate to axis
        gam(1)=gam(2)-rho_3d(2)*(gam(3)-gam(2))/(rho_3d(3)-rho_3d(2))

      ENDIF

!-------------------------------------------------------------------------------
!<B**2>
!-------------------------------------------------------------------------------
      IF(PRESENT(B2_R)) THEN

!Initialization
        B2_R(1:nrho_r)=0.0_rspec

!Interpolate to user grid inside R,Z domain
        iflag=0
        message=''
        CALL LINEAR1_INTERP(nrho_3d,rho_3d,b2,nr,rho_r,b2_r,iflag,             &
     &                      message)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_FLUXAV_B(2)/'//message
          IF(iflag > 0) GOTO 9999

        ENDIF

        IF(nr < nrho_r) THEN

!Extrapolate to user grid outside R,Z domain
          B2_R(nr+1:nrho_r)=B2_R(nr)

        ENDIF

      ENDIF

!-------------------------------------------------------------------------------
!<1/B**2>
!-------------------------------------------------------------------------------
      IF(PRESENT(BM2_R)) THEN

!Initialization
        BM2_R(1:nrho_r)=0.0_rspec
        gt_3d(:)=0.0_rspec
        az_3d(:)=0.0_rspec
        v1(:)=0.0_rspec

!Calculate on internal grid 
        DO i=2,nrho_3d !Over radial nodes

          DO k=1,nzeta_3d !Over toroidal nodes

            gt_3d(1:ntheta_3d)=gsqrt_3d(i,1:ntheta_3d,k)                       &
     &                         /b_3d(i,1:ntheta_3d,k)**2
            az_3d(k)=SUM(wtheta_3d*gt_3d) !Over poloidal nodes

          ENDDO !Over toroidal nodes

          v1(i)=SUM(wzeta_3d*az_3d)/vp_3d(i) !Over toroidal nodes

        ENDDO !Over radial nodes

!Extrapolate to axis
        v1(1)=v1(2)-rho_3d(2)*(v1(3)-v1(2))/(rho_3d(3)-rho_3d(2))

!Interpolate to user grid inside R,Z domain
        iflag=0
        message=''
        CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,BM2_R,iflag,            &
     &                      message)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_FLUXAV_B(3)/'//message
          IF(iflag > 0) GOTO 9999

        ENDIF

        IF(nr < nrho_r) THEN

!Extrapolate to user grid outside R,Z domain
          BM2_R(nr+1:nrho_r)=BM2_R(nr)

        ENDIF

      ENDIF

!-------------------------------------------------------------------------------
!<grad(rho)**2/B**2>
!-------------------------------------------------------------------------------
      IF(PRESENT(GR2BM2_R)) THEN

!Initialization
        GR2BM2_R(1:nrho_r)=0.0_rspec
        gt_3d(:)=0.0_rspec
        az_3d(:)=0.0_rspec
        v1(:)=0.0_rspec

        DO i=2,nrho_3d !Over radial nodes

          DO k=1,nzeta_3d !Over toroidal nodes

            gt_3d(1:ntheta_3d)=(( gcyl_3d(3,i,1:ntheta_3d,k)                   &
     &                           *gcyl_3d(5,i,1:ntheta_3d,k)                   &
     &                           -gcyl_3d(2,i,1:ntheta_3d,k)                   &
     &                           *gcyl_3d(6,i,1:ntheta_3d,k))**2               &
     &                           +rcyl_3d(i,1:ntheta_3d,k)**2                  &
     &                          *(gcyl_3d(2,i,1:ntheta_3d,k)**2                &
     &                           +gcyl_3d(5,i,1:ntheta_3d,k)**2))              &
     &                           /gsqrt_3d(i,1:ntheta_3d,k)                    &
     &                           /b_3d(i,1:ntheta_3d,k)**2
            az_3d(k)=SUM(wtheta_3d*gt_3d) !Over poloidal nodes

          ENDDO !Over toroidal nodes

          v1(i)=SUM(wzeta_3d*az_3d)/vp_3d(i) !Over toroidal nodes

        ENDDO !Over radial nodes

!Extrapolate to axis
        v1(1)=v1(2)-rho_3d(2)*(v1(3)-v1(2))/(rho_3d(3)-rho_3d(2))

!Interpolate to user grid inside R,Z domain
        iflag=0
        message=''
        CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,GR2BM2_r,iflag,         &
     &                      message)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_FLUXAV_B(4)/'//message
          IF(iflag > 0) GOTO 9999

        ENDIF

        IF(nr < nrho_r) THEN

!Extrapolate to user grid outside R,Z domain
          GR2BM2_R(nr+1:nrho_r)=GR2BM2_R(nr)

        ENDIF

      ENDIF

!-------------------------------------------------------------------------------
!f_trap
!-------------------------------------------------------------------------------
      IF(PRESENT(FTRAP_R)) THEN

!Initialization
        FTRAP_R(1:nrho_r)=0.0_rspec
        gt_3d(:)=0.0_rspec
        az_3d(:)=0.0_rspec
        b(:)=0.0_rspec
        bmax(:)=0.0_rspec
        v1(:)=0.0_rspec

!<|B|> and B_max
        DO i=2,nrho_3d !Over radial nodes

          DO k=1,nzeta_3d !Over toroidal nodes

            DO j=1,ntheta_3d !Over poloidal nodes

              gt_3d(j)=b_3d(i,j,k)*gsqrt_3d(i,j,k)
              IF(b_3d(i,j,k) > bmax(i)) bmax(i)=b_3d(i,j,k)

            ENDDO !Over poloidal nodes

            az_3d(k)=SUM(wtheta_3d*gt_3d) !Over poloidal nodes

          ENDDO !Over toroidal nodes

          b(i)=SUM(wzeta_3d*az_3d)/vp_3d(i) !Over toroidal nodes

        ENDDO !Over radial nodes

!Extrapolate to axis
        b(1)=b(2)-rho_3d(2)*(b(3)-b(2))/(rho_3d(3)-rho_3d(2))
        bmax(1)=bmax(2)-rho_3d(2)*(bmax(3)-bmax(2))                            &
     &                           /(rho_3d(3)-rho_3d(2))

!Upper and lower trapped fractions bound solution
        DO i=2,nrho_3d !Over radial nodes

  !Flux surface average for upper estimate of trapped fraction
          DO k=1,nzeta_3d !Over toroidal nodes

            DO j=1,ntheta_3d !Over poloidal nodes

              h=b_3d(i,j,k)/bmax(i)
              gt_3d(j)=(1.0_rspec-SQRT(1.0_rspec-h)                            &
     &                 *(1.0_rspec+0.5_rspec*h))/h**2*gsqrt_3d(i,j,k)

            ENDDO !Over poloidal nodes

            az_3d(k)=SUM(wtheta_3d*gt_3d) !Over poloidal nodes

          ENDDO !Over toroidal nodes

          v1(i)=SUM(wzeta_3d*az_3d)/vp_3d(i) !Over toroidal nodes

!Trapped fraction normalized to sqrt(rho) = 0.25*f_trap_low + 0.75*f_trap_upper
          h=b(i)/bmax(i)
          v1(i)=(0.25_rspec*(1.0_rspec-b2(i)/bmax(i)**2*v1(i))                 &
     &          +0.75_rspec*(1.0_rspec-b2(i)/b(i)**2                           &
     &                     *(1.0_rspec-SQRT(1.0_rspec-h)                       &
     &          *(1.0_rspec+0.5_rspec*h))))/SQRT(rho_3d(i))

        ENDDO !Over radial nodes


!Extrapolate to axis using constant normalized trapped fraction 
        !Find index of first point away from axis in R,Z domain
        DO i=1,nrho_r !Over radial nodes

          j=i
          IF(rho_3d(j) > rhomin_3d) EXIT

        ENDDO !Over radial nodes

        IF(j > 1) v1(1:j-1)=v1(j)

!Interpolate to user grid inside R,Z domain
        iflag=0
        message=''
        CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,FTRAP_R,iflag,          &
     &                      message)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_FLUXAV_B(5)/'//message
          IF(iflag > 0) GOTO 9999

        ENDIF

        IF(nr < nrho_r) THEN

!Extrapolate to user grid outside R,Z domain
          FTRAP_R(nr+1:nrho_r)=FTRAP_R(nr)

        ENDIF

!Reintroduce dominant radial dependence
        FTRAP_R(1:nrho_r)=FTRAP_R(1:nrho_r)*SQRT(rho_r(1:nrho_r))

      ENDIF

!-------------------------------------------------------------------------------
!F_m
!-------------------------------------------------------------------------------
      IF(PRESENT(FM_R)) THEN

!Initialization
        FM_R(:,:)=0.0_rspec
        gt_3d(:)=0.0_rspec
        v1(:)=0.0_rspec
        v2(:)=0.0_rspec
        v_r(:)=0.0_rspec

        IF(nper_3d /= 0) THEN

!Non-axisymmetric plasma
          !q(rho) for fm_1 interpolation
          DO i=2,nrho_3d !Over radial nodes

            v1(i)=1.0_rspec/ABS(iotabar_3d(1,i))

          ENDDO !Over radial nodes

!Extrapolate to axis
          v1(1)=v1(2)-rho_3d(2)*(v1(3)-v1(2))/(rho_3d(3)-rho_3d(2))

!Interpolate to user grid inside R,Z domain
          iflag=0
          message=''
          CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,v_r,iflag,            &
     &                        message)

          !Check messages
          IF(iflag /= 0) THEN

            message='AJAX_FLUXAV_B(6)/'//message
            IF(iflag > 0) GOTO 9999

          ENDIF

          IF(nr < nrho_r) THEN

!Extrapolate to user grid outside R,Z domain
            v_r(nr+1:nrho_r)=v_r(nr)

          ENDIF

!Reintroduce dominant radial dependence
          !If first node is not at the axis, include in calculation
          j=2
          IF(rho_r(1) > rhores_3d) j=1

          DO i=j,nrho_r !Over radial nodes

            !fm_r(2)=eps**1.5
            FM_R(2,i)=(rho_r(i)/r000_3d)**1.5

            !fm_r(1)=R_0*q/eps**1.5
            FM_R(1,i)=r000_3d*v_r(i)/fm_r(2,i)

          ENDDO !Over radial nodes

        ELSE

!Axisymmetric plasma
          DO m=1,SIZE(FM_R,1) !Over poloidal moments

            DO i=2,nrho_3d !Over radial nodes

              !sin(Theta) terms
              DO j=1,ntheta_3d !Over poloidal nodes

                IF(j == 1 .OR. j == ntheta_3d) THEN

                  !d(B)/d(theta) at end points noting periodicity
                  dbdt=(b_3d(i,2,1)-b_3d(i,ntheta_3d-1,1))                     &
     &                 /2.0_rspec/dtheta_3d

                ELSE

                  dbdt=(b_3d(i,j+1,1)-b_3d(i,j-1,1))/2.0_rspec/dtheta_3d

                ENDIF

                !sin(m*Theta)*B^theta/B*(dB/dtheta)
                gt_3d(j)=gsqrt_3d(i,j,1)*btheta_3d(i,j,1)/b_3d(i,j,1)          &
     &                   *SIN(m*captheta(i,j))*dbdt

              ENDDO !Over poloidal nodes

              v1(i)=2.0_rspec*z_pi*SUM(wtheta_3d*gt_3d)                        &
     &              /vp_3d(i) !Over poloidal nodes

              DO j=1,ntheta_3d !Over poloidal nodes

                !sin(m*Theta)*B^theta*(dB/dtheta)
                gt_3d(j)=gt_3d(j)*b_3d(i,j,1)*gam(i)

              ENDDO !Over poloidal nodes

              v2(i)=v1(i)*2.0_rspec*z_pi*SUM(wtheta_3d*gt_3d)                  &
     &              /vp_3d(i) !Over poloidal nodes

              !cos(Theta) terms
              DO j=1,ntheta_3d !Over poloidal nodes

                IF(j == 1 .OR. j == ntheta_3d) THEN

                  !d(B)/d(theta) at end points noting periodicity
                  dbdt=(b_3d(i,2,1)-b_3d(i,ntheta_3d-1,1))                     &
     &                 /(2.0_rspec*dtheta_3d)

                ELSE

                  dbdt=(b_3d(i,j+1,1)-b_3d(i,j-1,1))/2.0_rspec/dtheta_3d

                ENDIF

                !cos(m*Theta)*B^theta/B*(dB/dtheta)
                gt_3d(j)=gsqrt_3d(i,j,1)*btheta_3d(i,j,1)/b_3d(i,j,1)          &
     &                   *COS(m*captheta(i,j))*dbdt

              ENDDO !Over poloidal nodes

              v1(i)=2.0_rspec*z_pi*SUM(wtheta_3d*gt_3d)                        &
     &              /vp_3d(i) !Over poloidal nodes

              !cos(m*Theta)*B^theta*(dB/dtheta)
              gt_3d(:)=gt_3d(:)*b_3d(i,:,1)*gam(i)
              v2(i)=v2(i)+v1(i)*2.0_rspec*z_pi*SUM(wtheta_3d*gt_3d)            &
     &              /vp_3d(i) !Over poloidal nodes

              !<B^theta>
              gt_3d(:)=gsqrt_3d(i,:,1)*btheta_3d(i,:,1)
              v1(i)=2.0_rspec*z_pi*SUM(wtheta_3d*gt_3d)                        &
     &              /vp_3d(i) !Over poloidal nodes

              !Multiply by 2/<B^theta>/<B**2> and normalize to rho**1.5
              v2(i)=v2(i)*2.0_rspec/v1(i)/b2(i)/rho_3d(i)**1.5

            ENDDO !Over radial nodes

!Extrapolate to axis
            v2(1)=v2(2)-rho_3d(2)*(v2(3)-v2(2))/(rho_3d(3)-rho_3d(2))

!Interpolate to user grid inside R,Z domain
            iflag=0
            message=''
            CALL LINEAR1_INTERP(nrho_3d,rho_3d,v2,nr,rho_r,v_r,iflag,          &
     &                          message)

            !Check messages
            IF(iflag /= 0) THEN

              message='AJAX_FLUXAV_B(7)/'//message
              IF(iflag > 0) GOTO 9999

            ENDIF

            IF(nr < nrho_r) THEN

!Extrapolate to user grid outside R,Z domain
              v_r(nr+1:nrho_r)=v_r(nr)

            ENDIF

!Reintroduce dominant radial dependence
            FM_R(m,1:nrho_r)=v_r(1:nrho_r)*rho_r(1:nrho_r)**1.5

          ENDDO !Over poloidal moments

        ENDIF

      ENDIF

!-------------------------------------------------------------------------------
!n.grad(Theta)
!-------------------------------------------------------------------------------
      IF(PRESENT(GRTH_R)) THEN

!Initialization
        GRTH_R(1:nrho_r)=0.0_rspec

!Interpolate to user grid inside R,Z domain
        IF(nper_3d == 0) THEN

          iflag=0
          message=''
          CALL LINEAR1_INTERP(nrho_3d,rho_3d,gam,nr,rho_r,GRTH_R,iflag,        &
     &                        message)

          !Check messages
          IF(iflag /= 0) THEN

            message='AJAX_FLUXAV_B(8)/'//message
            IF(iflag > 0) GOTO 9999

          ENDIF

        ENDIF

        IF(nr < nrho_r) THEN

!Extrapolate to user grid outside R,Z domain
          GRTH_R(nr+1:nrho_r)=GRTH_R(nr)

        ENDIF

      ENDIF

!-------------------------------------------------------------------------------
!Susceptance matrix element S11
!-------------------------------------------------------------------------------
      IF(PRESENT(SUS11_R)) THEN

!Initialization
        SUS11_R(1:nrho_r)=0.0_rspec
        gt_3d(:)=0.0_rspec
        az_3d(:)=0.0_rspec
        v1(:)=0.0_rspec

!Calculate on internal grid
        DO i=2,nrho_3d !Over radial nodes

          DO k=1,nzeta_3d !Over toroidal nodes

            gt_3d(1:ntheta_3d)=( gcyl_3d(2,i,1:ntheta_3d,k)**2                 &
     &                          +gcyl_3d(5,i,1:ntheta_3d,k)**2 )               &
     &                         /gsqrt_3d(i,1:ntheta_3d,k)
            az_3d(k)=SUM(wtheta_3d*gt_3d) !Over poloidal nodes

          ENDDO !Over toroidal nodes

          !Remove dominant radial dependence
          v1(i)=SUM(wzeta_3d*az_3d)/rho_3d(i) !Over toroidal nodes

        ENDDO !Over radial nodes

!Extrapolate to axis
        v1(1)=v1(2)-rho_3d(2)*(v1(3)-v1(2))/(rho_3d(3)-rho_3d(2))

!Interpolate to user grid inside R,Z domain
        iflag=0
        message=''
        CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,SUS11_R,iflag,          &
     &                      message)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_FLUXAV_B(9)/'//message
          IF(iflag > 0) GOTO 9999

        ENDIF

        IF(nr < nrho_r) THEN

!Extrapolate to user grid outside R,Z domain
          SUS11_R(nrho_r+1:nr)=SUS11_R(nr)

        ENDIF

!Reintroduce dominant radial dependence
        SUS11_R(1:nrho_r)=SUS11_R(1:nrho_r)*rho_r(1:nrho_r)/4.0_rspec          &
     &                    /z_pi**2

      ENDIF

!-------------------------------------------------------------------------------
!Susceptance matrix element S12
!-------------------------------------------------------------------------------
      IF(PRESENT(SUS12_R)) THEN

!Initialization
        SUS12_R(1:nrho_r)=0.0_rspec
        gt_3d(:)=0.0_rspec
        az_3d(:)=0.0_rspec
        v1(:)=0.0_rspec

!Calculate on internal grid
        DO i=2,nrho_3d !Over radial nodes

          DO k=1,nzeta_3d !Over toroidal nodes

            gt_3d(1:ntheta_3d)=( ( gcyl_3d(2,i,1:ntheta_3d,k)                  &
     &                           *gcyl_3d(3,i,1:ntheta_3d,k)                   &
     &                           +gcyl_3d(5,i,1:ntheta_3d,k)                   &
     &                           *gcyl_3d(6,i,1:ntheta_3d,k))                  &
     &               *(1.0_rspec+eltheta_3d(i,1:ntheta_3d,k))                  &
     &                          -( gcyl_3d(2,i,1:ntheta_3d,k)**2               &
     &                            +gcyl_3d(5,i,1:ntheta_3d,k)**2 )             &
     &                          *elzeta_3d(i,1:ntheta_3d,k) )                  &
     &                        /gsqrt_3d(i,1:ntheta_3d,k)
            az_3d(k)=SUM(wtheta_3d*gt_3d) !Over poloidal nodes

          ENDDO !Over toroidal nodes

          !Remove dominant radial dependence
          v1(i)=SUM(wzeta_3d*az_3d)/rho_3d(i) !Over toroidal nodes

        ENDDO !Over radial nodes

!Extrapolate to axis
        v1(1)=v1(2)-rho_3d(2)*(v1(3)-v1(2))/(rho_3d(3)-rho_3d(2))

!Interpolate to user grid inside R,Z domain
        iflag=0
        message=''
        CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,SUS12_R,iflag,          &
     &                      message)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_FLUXAV_B(10)/'//message
          IF(iflag > 0) GOTO 9999

        ENDIF

        IF(nr < nrho_r) THEN

!Extrapolate to user grid outside R,Z domain
          sus12_r(nr+1:nrho_r)=SUS12_R(nr)

        ENDIF

!Reintroduce dominant radial dependence
        SUS12_R(1:nrho_r)=SUS12_R(1:nrho_r)*rho_r(1:nrho_r)/4.0_rspec          &
     &                    /z_pi**2

      ENDIF

!-------------------------------------------------------------------------------
!Susceptance matrix element S21
!-------------------------------------------------------------------------------
      IF(PRESENT(SUS21_R)) THEN

!Initialization
        SUS21_R(1:nrho_r)=0.0_rspec
        gt_3d(:)=0.0_rspec
        az_3d(:)=0.0_rspec
        v1(:)=0.0_rspec

!Calculate on internal grid
        DO i=2,nrho_3d !Over radial nodes

          DO k=1,nzeta_3d !Over toroidal nodes

            gt_3d(1:ntheta_3d)=( gcyl_3d(2,i,1:ntheta_3d,k)                    &
     &                         *gcyl_3d(3,i,1:ntheta_3d,k)                     &
     &                         +gcyl_3d(5,i,1:ntheta_3d,k)                     &
     &                         *gcyl_3d(6,i,1:ntheta_3d,k) )                   &
     &                        /gsqrt_3d(i,1:ntheta_3d,k)
            az_3d(k)=SUM(wtheta_3d*gt_3d) !Over poloidal nodes

          ENDDO !Over toroidal nodes

          !Remove dominant radial dependence
          v1(i)=SUM(wzeta_3d*az_3d)/rho_3d(i) !Over toroidal nodes

        ENDDO !Over radial nodes

!Extrapolate to axis
        v1(1)=v1(2)-rho_3d(2)*(v1(3)-v1(2))/(rho_3d(3)-rho_3d(2))

!Interpolate to user grid inside R,Z domain
        iflag=0
        message=''
        CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,SUS21_R,iflag,          &
     &                      message)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_FLUXAV_B(11)/'//message
          IF(iflag > 0) GOTO 9999

        ENDIF

        IF(nr < nrho_r) THEN

!Extrapolate to user grid outside R,Z domain
          SUS21_R(nr+1:nrho_r)=SUS21_R(nr)

        ENDIF

!Reintroduce dominant radial dependence
        SUS21_R(1:nrho_r)=SUS21_R(1:nrho_r)*rho_r(1:nrho_r)/4.0_rspec          &
     &                    /z_pi**2

      ENDIF

!-------------------------------------------------------------------------------
!Susceptance matrix element S22
!-------------------------------------------------------------------------------
      IF(PRESENT(SUS22_R)) THEN

!Initialization
        SUS22_R(1:nrho_r)=0.0_rspec
        gt_3d(:)=0.0_rspec
        az_3d(:)=0.0_rspec
        v1(:)=0.0_rspec

!Calculate on internal grid
        DO i=2,nrho_3d !Over radial nodes

          DO k=1,nzeta_3d !Over toroidal nodes

            gt_3d(1:ntheta_3d)=( ( rcyl_3d(i,1:ntheta_3d,k)**2                 &
     &                           +gcyl_3d(3,i,1:ntheta_3d,k)**2                &
     &                           +gcyl_3d(6,i,1:ntheta_3d,k)**2 )              &
     &                          *(1.0_rspec+eltheta_3d(i,1:ntheta_3d,k))       &
     &                          -( gcyl_3d(2,i,1:ntheta_3d,k)                  &
     &                            *gcyl_3d(3,i,1:ntheta_3d,k)                  &
     &                            +gcyl_3d(5,i,1:ntheta_3d,k)                  &
     &                            *gcyl_3d(6,i,1:ntheta_3d,k) )                &
     &                           *elzeta_3d(i,1:ntheta_3d,k) )                 &
     &                        /gsqrt_3d(i,1:ntheta_3d,k)
            az_3d(k)=SUM(wtheta_3d*gt_3d) !Over poloidal nodes

          ENDDO !Over toroidal nodes

!Remove dominant radial dependence
          v1(i)=SUM(wzeta_3d*az_3d)*rho_3d(i) !Over toroidal nodes

        ENDDO !Over radial nodes

!Extrapolate to axis
        v1(1)=v1(2)-rho_3d(2)*(v1(3)-v1(2))/(rho_3d(3)-rho_3d(2))

!Interpolate to user grid inside R,Z domain
        iflag=0
        message=''
        CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,SUS22_R,iflag,          &
     &                      message)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_FLUXAV_B(12)/'//message
          IF(iflag > 0) GOTO 9999

        ENDIF

        IF(nr < nrho_r) THEN

!Extrapolate to user grid outside R,Z domain
          SUS22_R(nr+1:nrho_r)=SUS22_R(nr)

        ENDIF

!Reintroduce dominant radial dependence
        j=2
        !If first node is not at the axis, include in calculation
        IF(rho_r(1) > rhores_3d) j=1

        SUS22_R(j:nrho_r)=SUS22_R(j:nrho_r)/rho_r(j:nrho_r)/4.0_rspec          &
     &                    /z_pi**2

      ENDIF

      IF(j == 2) SUS22_R(1)=SUS22_R(2)

 9999 CONTINUE

      END SUBROUTINE AJAX_FLUXAV_B

      SUBROUTINE AJAX_FLUXAV_G(nrho_r,rho_r,iflag,message,                     &
     &                         AREA_R,DVOL_R,GRHO1_R,GRHO2_R,GRHO2RM2_R,       &
     &                         RM2_R,VOL_R,VP_R)
!-------------------------------------------------------------------------------
!AJAX_FLUXAV_G gets flux surface quantities that depend on geometry
!References:
!  W.A.Houlberg, P.I.Strand 8/2001
!Input:
!  nrho_r              -no. of radial nodes [-]
!  rho_r(nrho_r)       -radial nodes [rho]
!Output:
!  iflag               -error and warning flag [-]
!                      =-1 warning
!                      =0 no warnings or errors
!                      =1 error
!  message             -warning or error message [character]
!Optional output:
!  AREA_R(nrho_r)      -surface area [m**2]
!  DVOL_R(nrho_r)      -cell volume between rho_r(i) and rho_r(i+1) [m**3]
!  GRHO1_R(nrho_r)     -<|grad(rho)|> [rho/m]  
!  GRHO2_R(nrho_r)     -<grad(rho)**2> [rho**2/m**4]
!  GRHO2RM2_R(nrho_r)  -<grad(rho)**2/R**2>) [rho**2/m**2]
!  RM2_R(nrho_r)       -<1/R**2> [/m**2]
!  VOL_R(nrho_r)       -enclosed volume [-]
!  VP_R(nrho_r)        -dV/drho [m**3/rho]
!-------------------------------------------------------------------------------

!Declaration of input variables
      INTEGER, INTENT(IN) ::                                                   &
     & nrho_r

      REAL(KIND=rspec), INTENT(IN) ::                                          &
     & rho_r(:)

!Declaration of output variables
      INTEGER, INTENT(OUT) ::                                                  &
     & iflag

      CHARACTER(len=*), INTENT(OUT) ::                                         &
     & message

!Declaration of optional output variables
      REAL(KIND=rspec), INTENT(OUT), OPTIONAL ::                               &
     & AREA_R(:),DVOL_R(:),GRHO1_R(:),GRHO2_R(:),GRHO2RM2_R(:),RM2_R(:),       &
     & VOL_R(:),VP_R(:)

!Declaration of local variables
      INTEGER ::                                                               &
     & i,k,                                                                    &
     & nr

      REAL(KIND=rspec) ::                                                      &
     & vp(1:nrho_r),gr(1:nrho_r),                                              &
     & v1(1:nrho_3d)

!-------------------------------------------------------------------------------
!Check whether flux surface averaging arrays have been set
!-------------------------------------------------------------------------------
      IF(.NOT. l_fluxavg_3d) THEN

        iflag=0
        message=''
        CALL AJAX_INIT_FLUXAV_G(iflag,message)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_FLUXAV_G(1)/'//message
          GOTO 9999

        ENDIF

      ENDIF

!-------------------------------------------------------------------------------
!Check whether values outside R,Z domain are requested (nr is last point inside)
!-------------------------------------------------------------------------------
      DO i=nrho_r,1,-1

        nr=i
        IF(rho_r(nr) < rhomax_3d+rhores_3d) EXIT

      ENDDO

!-------------------------------------------------------------------------------
!d(V)/d(rho) on user grid for area_r, dvol_r, vol_r and vp_r
!-------------------------------------------------------------------------------
      IF(PRESENT(AREA_R) .OR.                                                  &
     &   PRESENT(DVOL_R) .OR.                                                  &
     &   PRESENT(VOL_R) .OR.                                                   &
     &   PRESENT(VP_R)) THEN

!Initialization
        vp(:)=0.0_rspec
        v1(:)=0.0_rspec

!Remove dominant radial dependence
        v1(2:nrho_3d)=vp_3d(2:nrho_3d)/rho_3d(2:nrho_3d)

!Extrapolate to axis
        v1(1)=v1(2)-rho_3d(2)*(v1(3)-v1(2))/(rho_3d(3)-rho_3d(2))

!Interpolate to user grid inside R,Z domain
        iflag=0
        message=''
        CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,vp,iflag,               &
     &                      message)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_FLUXAV_G(2)/'//message
          IF(iflag > 0) GOTO 9999

        ENDIF

        IF(nr < nrho_r) THEN

!Extrapolate to user grid outside R,Z domain
          vp(nr+1:nrho_r)=vp(nr)

        ENDIF

!Reintroduce dominant radial dependence
        vp(1:nrho_r)=vp(1:nrho_r)*rho_r(1:nrho_r)

      ENDIF

!-------------------------------------------------------------------------------
!<|grad(rho)|> on user grid for area_r and grho1_r
!-------------------------------------------------------------------------------
      IF(PRESENT(AREA_R) .OR.                                                  &
     &   PRESENT(GRHO1_R)) THEN

!Initialization
        gr(:)=0.0_rspec
        gt_3d(:)=0.0_rspec
        az_3d(:)=0.0_rspec
        v1(:)=0.0_rspec

!Calculate on internal grid
        DO i=2,nrho_3d !Over radial nodes

          DO k=1,nzeta_3d !Over toroidal nodes

            gt_3d(1:ntheta_3d)=SQRT(( ( gcyl_3d(3,i,1:ntheta_3d,k)             &
     &                                *gcyl_3d(5,i,1:ntheta_3d,k)              &
     &                                -gcyl_3d(2,i,1:ntheta_3d,k)              &
     &                                *gcyl_3d(6,i,1:ntheta_3d,k))**2          &
     &                               +rcyl_3d(i,1:ntheta_3d,k)**2              &
     &                               *(gcyl_3d(2,i,1:ntheta_3d,k)**2           &
     &                                +gcyl_3d(5,i,1:ntheta_3d,k)**2) ))
            az_3d(k)=SUM(wtheta_3d*gt_3d) !Over poloidal nodes

          ENDDO !Over toroidal nodes

          v1(i)=SUM(wzeta_3d*az_3d)/vp_3d(i) !Over toroidal nodes

        ENDDO !Over radial nodes

        v1(1)=v1(2)

!Interpolate to user grid inside R,Z domain
        iflag=0
        message=''
        CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,gr,iflag,               &
     &                      message)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_FLUXAV_G(3)/'//message
          IF(iflag > 0) GOTO 9999

        ENDIF

        IF(nr < nrho_r) THEN

!Extrapolate to user grid outside R,Z domain
          gr(nr+1:nrho_r)=gr(nr)

        ENDIF

      ENDIF

!-------------------------------------------------------------------------------
!Area
!-------------------------------------------------------------------------------
      IF(PRESENT(AREA_R)) THEN

!Initialization
        AREA_R(:)=0.0_rspec

!Use temporary arrays
        AREA_R(1:nrho_r)=vp(1:nrho_r)*gr(1:nrho_r)

      ENDIF

!-------------------------------------------------------------------------------
!delta(Volume)
!-------------------------------------------------------------------------------
      IF(PRESENT(DVOL_R)) THEN

!Initialization
        DVOL_R(:)=0.0_rspec

!Calculate from integral rho*drho*(V'/rho), where (V'/rho) is a weak function
        !Allow for first node to be away from the axis
        IF(rho_r(1) > rhores_3d) THEN

          !First node is off axis
          k=1

        ELSE

          !First node is on axis
          k=2
          DVOL_R(1)=0.5_rspec*vp(2)*rho_r(2)

        ENDIF

        DO i=k,nrho_r-1 !Over radial nodes

          DVOL_R(i)=0.5_rspec*(vp(i)/rho_r(i)+vp(i+1)/rho_r(i+1))              &
     &              *0.5_rspec*(rho_r(i+1)**2-rho_r(i)**2)

        ENDDO !Over radial nodes

!Set ghost node value equal to last node
        DVOL_R(nrho_r)=DVOL_R(nrho_r-1)

      ENDIF

!-------------------------------------------------------------------------------
!<|grad(rho)|>
!-------------------------------------------------------------------------------
      IF(PRESENT(GRHO1_R)) THEN

!Initialization
        GRHO1_R(:)=0.0_rspec

!Copy from temporary array
        GRHO1_R(1:nrho_r)=gr(1:nrho_r)

      ENDIF

!-------------------------------------------------------------------------------
!<grad(rho)**2>
!-------------------------------------------------------------------------------
      IF(PRESENT(GRHO2_R)) THEN

!Initialization
        GRHO2_R(:)=0.0_rspec
        gt_3d(:)=0.0_rspec
        az_3d(:)=0.0_rspec
        v1(:)=0.0_rspec

!Calculate on internal grid
        DO i=2,nrho_3d !Over radial nodes

          DO k=1,nzeta_3d !Over toroidal nodes

            gt_3d(1:ntheta_3d)=( (gcyl_3d(3,i,1:ntheta_3d,k)                   &
     &                          *gcyl_3d(5,i,1:ntheta_3d,k)                    &
     &                          -gcyl_3d(2,i,1:ntheta_3d,k)                    &
     &                          *gcyl_3d(6,i,1:ntheta_3d,k))**2                &
     &                          +rcyl_3d(i,1:ntheta_3d,k)**2                   &
     &                          *(gcyl_3d(2,i,1:ntheta_3d,k)**2                &
     &                          +gcyl_3d(5,i,1:ntheta_3d,k)**2) )              &
     &                        /gsqrt_3d(i,1:ntheta_3d,k)
            az_3d(k)=SUM(wtheta_3d*gt_3d) !Over poloidal nodes

          ENDDO !Over toroidal nodes

          v1(i)=SUM(wzeta_3d*az_3d)/vp_3d(i) !Over toroidal nodes

        ENDDO !Over radial nodes

!Extrapolate to axis
        v1(1)=v1(2)

!Interpolate to user grid inside R,Z domain
        iflag=0
        message=''
        CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,GRHO2_R,iflag,          &
     &                      message)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_FLUXAV_G(4)/'//message
          IF(iflag > 0) GOTO 9999

        ENDIF

        IF(nr < nrho_r) THEN

!Extrapolate to user grid outside R,Z domain
          GRHO2_R(nr+1:nrho_r)=GRHO2_R(nr)

        ENDIF

      ENDIF

!-------------------------------------------------------------------------------
!<grad(rho)**2/R**2>
!-------------------------------------------------------------------------------
      IF(PRESENT(GRHO2RM2_R)) THEN

!Initialization
        GRHO2RM2_r(:)=0.0_rspec
        gt_3d(:)=0.0_rspec
        az_3d(:)=0.0_rspec
        v1(:)=0.0_rspec

!Calculate on internal grid
        DO i=2,nrho_3d !Over radial nodes

          DO k=1,nzeta_3d !Over toroidal nodes

            gt_3d(1:ntheta_3d)=( (gcyl_3d(3,i,1:ntheta_3d,k)                   &
     &                          *gcyl_3d(5,i,1:ntheta_3d,k)                    &
     &                          -gcyl_3d(2,i,1:ntheta_3d,k)                    &
     &                          *gcyl_3d(6,i,1:ntheta_3d,k))**2                &
     &                          /rcyl_3d(i,1:ntheta_3d,k)**2                   &
     &                          +(gcyl_3d(2,i,1:ntheta_3d,k)**2                &
     &                          +gcyl_3d(5,i,1:ntheta_3d,k)**2) )              &
     &                        /gsqrt_3d(i,1:ntheta_3d,k)
            az_3d(k)=SUM(wtheta_3d*gt_3d) !Over poloidal nodes

          ENDDO !Over toroidal nodes

          v1(i)=SUM(wzeta_3d*az_3d)/vp_3d(i) !Over toroidal nodes

        ENDDO !Over radial nodes

!Extrapolate to axis
        v1(1)=v1(2)

!Interpolate to user grid inside R,Z domain
        iflag=0
        message=''
        CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,GRHO2RM2_R,iflag,       &
     &                      message)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_FLUXAV_G(5)/'//message
          IF(iflag > 0) GOTO 9999

        ENDIF

        IF(nr < nrho_r) THEN

!Extrapolate to user grid outside R,Z domain
          GRHO2RM2_R(nr+1:nrho_r)=GRHO2RM2_R(nr)

        ENDIF

      ENDIF

!-------------------------------------------------------------------------------
!<1/R**2>
!-------------------------------------------------------------------------------
      IF(PRESENT(RM2_R)) THEN

!Initialization
        RM2_R(:)=0.0_rspec
        gt_3d(:)=0.0_rspec
        az_3d(:)=0.0_rspec
        v1(:)=0.0_rspec

!Calculate on internal grid
        DO i=2,nrho_3d !Over radial nodes

          DO k=1,nzeta_3d !Over toroidal nodes

            gt_3d(1:ntheta_3d)=gsqrt_3d(i,1:ntheta_3d,k)                       &
     &                         /rcyl_3d(i,1:ntheta_3d,k)**2
            az_3d(k)=SUM(wtheta_3d*gt_3d) !Over poloidal nodes
          ENDDO !Over toroidal nodes

          v1(i)=SUM(wzeta_3d*az_3d)/vp_3d(i) !Over toroidal nodes

        ENDDO !Over radial nodes

!Extrapolate to axis
        v1(1)=v1(2)

!Interpolate to user grid inside R,Z domain
        iflag=0
        message=''
        CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,RM2_R,iflag,            &
     &                      message)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_FLUXAV_G(6)/'//message
          IF(iflag > 0) GOTO 9999

        ENDIF

        IF(nr < nrho_r) THEN

!Extrapolate to user grid outside R,Z domain
          RM2_R(nr+1:nrho_r)=RM2_R(nr)

        ENDIF

      ENDIF

!-------------------------------------------------------------------------------
!Volume
!-------------------------------------------------------------------------------
      IF(PRESENT(VOL_R)) THEN

!Initialization
        VOL_R(:)=0.0_rspec

!Sum up dvol elements
        !Allow for first node to be away from the axis

        IF(rho_r(1) > rhores_3d) THEN

          !First node is off axis
          k=1
          vol_r(1)=0.5_rspec*vp(1)*rho_r(1)

        ELSE

          !First node is on axis
          k=2
          VOL_R(1)=0.0_rspec
          VOL_R(2)=0.5_rspec*vp(2)*rho_r(2)

        ENDIF

        DO i=k+1,nrho_r !Over radial nodes

          VOL_R(i)=VOL_R(i-1)                                                  &
     &             +0.5_rspec*(vp(i-1)/rho_r(i-1)+vp(i)/rho_r(i))              &
     &             *0.5_rspec*(rho_r(i)**2-rho_r(i-1)**2)

        ENDDO !Over radial nodes

      ENDIF

!-------------------------------------------------------------------------------
!d(Volume)/d(rho)
!-------------------------------------------------------------------------------
      IF(PRESENT(VP_R)) THEN

!Initialization
        VP_R(:)=0.0_rspec

!Copy from temporary array
        VP_R(1:nrho_r)=vp(1:nrho_r)

      ENDIF

 9999 CONTINUE

      END SUBROUTINE AJAX_FLUXAV_G

      SUBROUTINE AJAX_I(nrho_r,rho_r,iflag,message,                            &
     &                  CUR_I_R,CUR_F_R,CUR_IMN_R,CUR_IMX_R,CUR_FMN_R,         &
     &                  CUR_FMX_R)
!-------------------------------------------------------------------------------
!AJAX_I gets the enclosed toroidal and external poloidal currents for a set of
!  flux surfaces from the covariant components of B, as well as maximum and
!  minimum values for use as a diagnostic on the accuracy of the stream function
!References:
!  W.A.Houlberg, P.I.Strand 8/2001
!Input:
!  nrho_r              -no. of radial nodes [-]
!  rho_r(nrho_r)       -radial nodes [rho]
!Output:
!  iflag               -error and warning flag [-]
!                      =-1 warning
!                      =0 none
!                      =1 error
!  message             -warning or error message [character]
!Optional output:
!  CUR_I_R(nrho_r)     -enclosed toroidal current [A]
!  CUR_F_R(nrho_r)     -external poloidal current [A]
!  CUR_IMN_R(nrho_r)   -minimum toroidal current on AJAX toroidal grid [A]
!  CUR_IMX_R(nrho_r)   -maximum toroidal current on AJAX toroidal grid [A]
!  CUR_FMN_R(nrho_r)   -minimum poloidal current on AJAX poloidal grid [A]
!  CUR_FMX_R(nrho_r)   -maximum poloidal current on AJAX poloidal grid [A]
!-------------------------------------------------------------------------------

!Declaration of input variables
      INTEGER, INTENT(IN) ::                                                   &
     & nrho_r

      REAL(KIND=rspec), INTENT(IN) ::                                          &
     & rho_r(:)

!Declaration of output variables
      INTEGER, INTENT(OUT) ::                                                  &
     & iflag

      CHARACTER(len=*), INTENT(OUT) ::                                         &
     & message

!Declaration of optional output variables
      REAL(KIND=rspec), INTENT(OUT), OPTIONAL ::                               &
     & CUR_I_R(:),CUR_F_R(:),CUR_IMN_R(:),CUR_IMX_R(:),CUR_FMN_R(:),           &
     & CUR_FMX_R(:)

!Declaration of local variables
      INTEGER ::                                                               &
     & i,j,k,imin

      REAL(KIND=rspec) ::                                                      &
     & cmax,cmin

      REAL(KIND=rspec) ::                                                      &
     & r_flx(3),r_cyl(3),g_cyl(6)

      REAL(KIND=rspec), ALLOCATABLE ::                                         &
     & b_co(:,:,:),gt(:),gz(:)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Allocate covariant B field array
      ALLOCATE(b_co(1:3,1:ntheta_3d,1:nzeta_3d))

      IF(PRESENT(CUR_I_R) .OR.                                                 &
     &   PRESENT(CUR_IMN_R) .OR.                                               &
     &   PRESENT(CUR_IMX_R)) ALLOCATE(gz(1:nzeta_3d))

      IF(PRESENT(CUR_F_R) .OR.                                                 &
     &   PRESENT(CUR_FMN_R) .OR.                                               &
     &   PRESENT(CUR_FMX_R)) ALLOCATE(gt(1:ntheta_3d))

      IF(PRESENT(CUR_I_R)) CUR_I_R(1:nrho_r)=0.0_rspec
      IF(PRESENT(CUR_IMN_R)) CUR_IMN_R(1:nrho_r)=0.0_rspec
      IF(PRESENT(CUR_IMX_R)) CUR_IMX_R(1:nrho_r)=0.0_rspec
      IF(PRESENT(CUR_F_R)) CUR_F_R(1:nrho_r)=0.0_rspec
      IF(PRESENT(CUR_FMN_R)) CUR_FMN_R(1:nrho_r)=0.0_rspec
      IF(PRESENT(CUR_FMX_R)) CUR_FMX_R(1:nrho_r)=0.0_rspec

      imin=1
      IF(rho_r(1) <= rhomin_3d) imin=2

!-------------------------------------------------------------------------------
!Perform integrals of covariant B components on each surface
!-------------------------------------------------------------------------------
      DO i=imin,nrho_r !Over radial nodes

        r_flx(1)=rho_r(i)

        DO k=1,nzeta_3d !Over toroidal nodes

          r_flx(3)=REAL(k-1,rspec)*dzeta_3d

          DO j=1,ntheta_3d !Over poloidal nodes

            r_flx(2)=REAL(j-1,rspec)*dtheta_3d

!Get cylindrical coordinates and metrics
            iflag=0
            message=''
            CALL AJAX_FLX2CYL(r_flx,r_cyl,iflag,message,                       &
     &                        G_CYL=g_cyl)

            !Check messages
            IF(iflag /= 0) THEN

              message='AJAX_I(1)/'//message
              IF(iflag > 0) GOTO 9999

            ENDIF

!Get covariant B components
            iflag=0
            message=''
            CALL AJAX_B(r_flx,r_cyl,g_cyl,iflag,message,                       &
     &                  B_CO=b_co(:,j,k))

            !Check messages
            IF(iflag /= 0) THEN

              message='AJAX_I(2)/'//message
              IF(iflag > 0) GOTO 9999

            ENDIF

          ENDDO !Over poloidal nodes

        ENDDO !Over toroidal nodes

!Toroidal current
        IF(PRESENT(CUR_I_R) .OR.                                               &
     &     PRESENT(CUR_IMN_R) .OR.                                             &
     &     PRESENT(CUR_IMX_R)) THEN

          cmax=-z_large
          cmin=z_large

          DO k=1,nzeta_3d !Over toroidal nodes

            gz(k)=SUM(wtheta_3d*b_co(2,1:ntheta_3d,k))/z_mu0
            IF(gz(k) < cmin) cmin=gz(k)
            IF(gz(k) > cmax) cmax=gz(k)

            IF(PRESENT(CUR_I_R)) CUR_I_R(i)=SUM(wzeta_3d*gz)                   &
     &                                      /2.0_rspec/z_pi
            IF(PRESENT(CUR_IMN_R)) CUR_IMN_R(i)=cmin
            IF(PRESENT(CUR_IMX_R)) CUR_IMX_R(i)=cmax

          ENDDO !Over toroidal nodes

        ENDIF

!Poloidal current
        IF(PRESENT(CUR_F_R) .OR.                                               &
     &     PRESENT(CUR_FMN_R) .OR.                                             &
     &     PRESENT(CUR_FMX_R)) THEN

          cmax=-z_large
          cmin=z_large

          DO j=1,ntheta_3d !Over poloidal nodes

            gt(j)=SUM(wzeta_3d*b_co(3,j,1:nzeta_3d))/z_mu0
            IF(gt(j) < cmin) cmin=gt(j)
            IF(gt(j) > cmax) cmax=gt(j)

            IF(PRESENT(CUR_F_R)) CUR_F_R(i)=SUM(wtheta_3d*gt)                  &
     &                                      /2.0_rspec/z_pi
            IF(PRESENT(CUR_FMN_R)) CUR_FMN_R(i)=cmin
            IF(PRESENT(CUR_FMX_R)) CUR_FMX_R(i)=cmax

          ENDDO !Over toroidal nodes

        ENDIF

      ENDDO !Over radial nodes

!Check whether axial values are requested
      IF(imin == 2) THEN

        !Use value of second node
        IF(PRESENT(CUR_I_R)) CUR_I_R(1)=CUR_I_R(2)
        IF(PRESENT(CUR_IMN_R)) CUR_IMN_R(1)=CUR_IMN_R(2)
        IF(PRESENT(CUR_IMX_R)) CUR_IMX_R(1)=CUR_IMX_R(2)
        IF(PRESENT(CUR_F_R)) CUR_F_R(1)=CUR_F_R(2)
        IF(PRESENT(CUR_FMN_R)) CUR_FMN_R(1)=CUR_FMN_R(2)
        IF(PRESENT(CUR_FMX_R)) CUR_FMX_R(1)=CUR_FMX_R(2)

      ENDIF

 9999 CONTINUE

!-------------------------------------------------------------------------------
!Cleanup
!-------------------------------------------------------------------------------
!Deallocate arrays
      DEALLOCATE(b_co)
      IF(ALLOCATED(gt)) DEALLOCATE(gt)
      IF(ALLOCATED(gz)) DEALLOCATE(gz)

      END SUBROUTINE AJAX_I

      SUBROUTINE AJAX_MAGFLUX(nrho_r,rho_r,iflag,message,                      &
     &                        IOTABAR_R,Q_R,PHIPRM_R,PSIPRM_R,PHI_R,           &
     &                        PSI_R)
!-------------------------------------------------------------------------------
!AJAX_MAGFLUX gets poloidal and toroidal magnetic fluxes
!References:
!  W.A.Houlberg, P.I.Strand 8/2001
!Input:
!  nrho_r              -no. of radial nodes [-]
!  rho_r(nrho_r)       -radial nodes [rho]
!Output:
!  iflag               -error and warning flag [-]
!                      =-1 warning
!                      =0 no warnings or errors
!                      =1 error
!  message             -warning or error message [character]
!Optional output:
!  IOTABAR_R(nrho_r)   -rotational transform = d(Psi)/d(Phi) [-]
!  Q_R(i)              -safety factor = d(Phi)/d(Psi) [-]
!  PHIPRM_R(nrho_r)    -d(Phi)/d(rho) [Wb/rho]
!  PSIPRM_R(nrho_r)    -d(Psi)/d(rho) [Wb/rho]
!  PHI_R(nrho_r)       -toroidal flux [Wb]
!  PSI_R(nrho_r)       -poloidal flux [Wb]
!-------------------------------------------------------------------------------

!Declaration of input variables
      INTEGER, INTENT(IN) ::                                                   &
     & nrho_r

      REAL(KIND=rspec), INTENT(IN) ::                                          &
     & rho_r(:)

!Declaration of output variables
      INTEGER, INTENT(OUT) ::                                                  &
     & iflag

      CHARACTER(len=*), INTENT(OUT) ::                                         &
     & message

!Declaration of output variables
      REAL(KIND=rspec), INTENT(OUT), OPTIONAL ::                               &
     & IOTABAR_R(:),Q_R(:),PHIPRM_R(:),PSIPRM_R(:),PHI_R(:),PSI_R(:)

!Declaration of local variables             
      INTEGER ::                                                               &
     & i,j,                                                                    &
     & nr

      INTEGER, PARAMETER ::                                                    &
     & k_vopt(3)=(/1,0,0/)

      REAL(KIND=rspec) ::                                                      &
     & value(3),                                                               &
     & iotabar_t(1:nrho_r),phiprm_t(1:nrho_r),values(1:nrho_r)

!-------------------------------------------------------------------------------
!Check whether values outside R,Z domain are requested (nr is last point inside)
!-------------------------------------------------------------------------------
      DO i=nrho_r,1,-1

        nr=i
        IF(rho_r(nr) < rhomax_3d+rhores_3d) EXIT

      ENDDO

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
      phiprm_t(:)=0.0_rspec
      iotabar_t(:)=0.0_rspec

!-------------------------------------------------------------------------------
!Get temporary base quantities for phiprm and iotabar
!-------------------------------------------------------------------------------
      phiprm_t(1:nrho_r)=2.0_rspec*rho_r(1:nrho_r)*phitot_3d                   &
     &                   /rhomax_3d**2
      IF(rho_r(1) < rhores_3d) phiprm_t(1)=2.0_rspec*rhores_3d*phitot_3d       &
     &                                     /rhomax_3d**2

!Inside MHD solution domain
      j=1

      DO i=1,nr !Over radial nodes

        CALL SPLINE1_EVAL(k_vopt,nrho_3d,rho_r(i),rho_3d,iotabar_3d,j,         &
     &                    value)
        iotabar_t(i)=value(1)

      ENDDO !Over radial nodes

      IF(nr < nrho_r) THEN

          !Extrapolate outside with constant slope in iotabar
          iotabar_t(nr+1:nrho_r)=iotabar_3d(1,nrho_3d)                         &
     &                           +(rho_r(nr+1:nrho_r)-rhomax_3d)               &
     &                           *(iotabar_3d(1,nrho_3d)                       &
     &                            -iotabar_3d(1,nrho_3d-1))                    &
     &                           /(rhomax_3d-rho_3d(nrho_3d-1))

      ENDIF

!-------------------------------------------------------------------------------
!Return the appropriate quantity
!-------------------------------------------------------------------------------
!Rotational transform
      IF(PRESENT(IOTABAR_R)) THEN

        IOTABAR_R(:)=0.0_rspec
        IOTABAR_R(1:nrho_r)=iotabar_t(1:nrho_r)

      ENDIF

!Safety factor
      IF(PRESENT(Q_R)) THEN

        Q_R(:)=0.0_rspec
        Q_R(1:nrho_r)=1.0_rspec/iotabar_t(1:nrho_r)

      ENDIF

!d(Phi)/d(rho)
      IF(PRESENT(PHIPRM_R)) THEN

        PHIPRM_R(:)=0.0_rspec
        PHIPRM_R(1:nrho_r)=phiprm_t(1:nrho_r)

      ENDIF

!d(Psi)/d(rho)
      IF(PRESENT(PSIPRM_R)) THEN

        PSIPRM_R(:)=0.0_rspec
        PSIPRM_R(1:nrho_r)=iotabar_t(1:nrho_r)*phiprm_t(1:nrho_r)

      ENDIF

!Total toroidal flux
      IF(PRESENT(PHI_R)) THEN

        PHI_R(:)=0.0_rspec
        PHI_R(1:nrho_r)=phitot_3d*(rho_r(1:nrho_r)/rhomax_3d)**2

      ENDIF

!Total poloidal flux
      IF(PRESENT(PSI_R)) THEN

        PSI_R(:)=0.0_rspec

        iflag=0
        message=''
        CALL SPLINE1_INTEG(1,nrho_3d,rho_3d,iotabar_3d,nr,rho_r,values,        &
     &                     iflag,message)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_MAGFLUX/'//message
          IF(iflag > 0) GOTO 9999

        ENDIF

        PSI_R(1:nr)=values(1:nr)*2.0_rspec*phitot_3d/rhomax_3d**2

        IF(nr < nrho_r) THEN

          DO i=nr+1,nrho_r !Over radial nodes

            PSI_R(i)=PSI_R(i-1)+phitot_3d/rhomax_3d**2                         &
     &               *(rho_r(i-1)*iotabar_t(i-1)+rho_r(i)*iotabar_t(i))        &
     &               *(rho_r(i)-rho_r(i-1))

          ENDDO !Over radial nodes

        ENDIF

      ENDIF

 9999 CONTINUE

      END SUBROUTINE AJAX_MAGFLUX

      SUBROUTINE AJAX_SHAPE(nrho_r,rho_r,iflag,message,                        &
     &                      ZETA,SHIFT_R,ELONG_R,TRIANG_R,RMAX_R,RMIN_R,       &
     &                      ZMAX_R,ZMIN_R,RBOT_R,RTOP_R)
!-------------------------------------------------------------------------------
!AJAX_SHAPE gets shift, elongation, and triangularity
!References:
!  W.A.Houlberg, P.I.Strand 8/2001
!Input:
!  nrho_r              -no. of radial nodes [-]
!  rho_r(nrho_r)       -radial nodes [rho]
!Output:
!  iflag               -error and warning flag [-]
!                      =-1 warning
!                      =0 no warnings or errors
!                      =1 error
!  message             -warning or error message [character]
!Optional input:
!  ZETA                -toroidal plane for the shape [-]
!                      =0.0 by default
!Optional output:
!  SHIFT_R(nrho_r)     -shift [-]
!  ELONG_R(nrho_r)     -elongation [-]
!  TRIANG_R(nrho_r)    -triangularity [-]
!  RMAX_R(nrho_r)      -maximum R of plasma in zeta plane [m]
!  RMIN_R(nrho_r)      -minimum R of plasma in zeta plane [m]
!  ZMAX_R(nrho_r)      -maximum Z of plasma in zeta plane [m]
!  ZMIN_R(nrho_r)      -minimum Z of plasma in zeta plane [m]
!  RBOT_R(nrho_r)      -R at ZMIN_R in zeta plane [m]
!  RTOP_R(nrho_r)      -R at ZMAX_R in zeta plane [m]
!-------------------------------------------------------------------------------

!Declaration of input variables
      INTEGER, INTENT(IN) ::                                                   &
     & nrho_r

      REAL(KIND=rspec), INTENT(IN) ::                                          &
     & rho_r(:)

!Declaration of output variables
      INTEGER, INTENT(OUT) ::                                                  &
     & iflag

      CHARACTER(len=*), INTENT(OUT) ::                                         &
     & message

!Declaration of optional input variables
      REAL(KIND=rspec), INTENT(IN), OPTIONAL ::                                &
     & ZETA

!Declaration of optional output variables
      REAL(KIND=rspec), INTENT(OUT), OPTIONAL ::                               &
     & ELONG_R(:),SHIFT_R(:),TRIANG_R(:),                                      &
     & RMAX_R(:),RMIN_R(:),ZMAX_R(:),ZMIN_R(:),RBOT_R(:),RTOP_R(:)

!Declaration of local variables
      LOGICAL ::                                                               &
     & l_rminmax

      INTEGER ::                                                               &
     & i,                                                                      &
     & nr

      REAL(KIND=rspec) ::                                                      &
     & a0,rg0

      REAL(KIND=rspec) ::                                                      &
     & r_cyl(3),r_flx(3)

      REAL(KIND=rspec) ::                                                      &
     & v1(1:nrho_3d),rbot(1:nrho_3d),rtop(1:nrho_3d),                          &
     & rmax(1:nrho_3d),rmin(1:nrho_3d),zmax(1:nrho_3d),zmin(1:nrho_3d)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
      r_flx(:)=0.0_rspec
      r_cyl(:)=0.0_rspec
      rbot(:)=0.0_rspec
      rtop(:)=0.0_rspec
      rmax(:)=0.0_rspec
      rmin(:)=0.0_rspec
      zmax(:)=0.0_rspec
      zmin(:)=0.0_rspec
      v1(:)=0.0_rspec

!Set the toroidal angle for the calculations
      IF(PRESENT(ZETA)) r_flx(3)=ZETA

!-------------------------------------------------------------------------------
!Check whether values outside R,Z domain are requested (nr is last point inside)
!-------------------------------------------------------------------------------
      DO i=nrho_r,1,-1

        nr=i
        IF(rho_r(nr) < rhomax_3d+rhores_3d) EXIT

      ENDDO

!-------------------------------------------------------------------------------
!Half diameter and center of plasma boundary
!-------------------------------------------------------------------------------
!Find coordinates where dR/dtheta=0 on inside for R_min
      l_rminmax=.TRUE.
      r_flx(1)=rhomax_3d
      r_flx(2)=z_pi
      iflag=0
      message=''
      CALL AJAX_MINMAX_RZ(l_rminmax,r_flx,r_cyl,iflag,message)

      !Check messages
      IF(iflag /= 0) THEN

        message='AJAX_SHAPE(1)/'//message
        IF(iflag > 0) GOTO 9999

      ENDIF

      rmin(nrho_3d)=r_cyl(1)

!Find coordinates where dR/dtheta=0 on outside for R_max
      l_rminmax=.TRUE.
      r_flx(1)=rhomax_3d
      r_flx(2)=0.0_rspec
      iflag=0
      message=''
      CALL AJAX_MINMAX_RZ(l_rminmax,r_flx,r_cyl,iflag,message)

      !Check messages
      IF(iflag /= 0) THEN

        message='AJAX_SHAPE(2)/'//message
        IF(iflag > 0) GOTO 9999

      ENDIF

      rmax(nrho_3d)=r_cyl(1)

!Major radius of center of plasma boundary
      rg0=0.5_rspec*(rmin(nrho_3d)+rmax(nrho_3d))

!Minor radius = half diameter
      a0=0.5_rspec*(rmax(nrho_3d)-rmin(nrho_3d))

!-------------------------------------------------------------------------------
!Find inside, outside (dR/theta=0) and top, bottom (dZ/dtheta=0) of each surface
!-------------------------------------------------------------------------------
      DO i=2,nrho_3d !Over radial nodes

!dR/dtheta=0 on inside
        l_rminmax=.TRUE.
        r_flx(1)=rho_3d(i)
        r_flx(2)=z_pi
        iflag=0
        message=''
        CALL AJAX_MINMAX_RZ(l_rminmax,r_flx,r_cyl,iflag,message)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_SHAPE(3)/'//message
          IF(iflag > 0) GOTO 9999

        ENDIF

        rmin(i)=r_cyl(1)

!dR/dtheta=0 on outside
        l_rminmax=.TRUE.
        r_flx(1)=rho_3d(i)
        r_flx(2)=0.0_rspec
        iflag=0
        message=''
        CALL AJAX_MINMAX_RZ(l_rminmax,r_flx,r_cyl,iflag,message)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_SHAPE(4)/'//message
          IF(iflag > 0) GOTO 9999

        ENDIF

        rmax(i)=r_cyl(1)

!dZ/dtheta=0 on bottom
        l_rminmax=.FALSE.
        r_flx(1)=rho_3d(i)
        r_flx(2)=0.5_rspec*z_pi
        iflag=0
        message=''
        CALL AJAX_MINMAX_RZ(l_rminmax,r_flx,r_cyl,iflag,message)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_SHAPE(5)/'//message
          IF(iflag > 0) GOTO 9999

        ENDIF

        rbot(i)=r_cyl(1)
        zmin(i)=r_cyl(3)

!dZ/dtheta=0 on top
        l_rminmax=.FALSE.
        r_flx(1)=rho_3d(i)
        r_flx(2)=-0.5_rspec*z_pi
        iflag=0
        message=''
        CALL AJAX_MINMAX_RZ(l_rminmax,r_flx,r_cyl,iflag,message)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_SHAPE(6)/'//message
          IF(iflag > 0) GOTO 9999

        ENDIF

        rtop(i)=r_cyl(1)
        zmax(i)=r_cyl(3)

      ENDDO !Over radial nodes

!-------------------------------------------------------------------------------
!Shift
!-------------------------------------------------------------------------------
      IF(PRESENT(SHIFT_R)) THEN

!Initialization
        SHIFT_R(1:nrho_r)=0.0_rspec

!Shift relative to center of outer surface   
        v1(2:nrho_3d)=(0.5_rspec*(rmax(2:nrho_3d)+rmin(2:nrho_3d))-rg0)        &
     &                /a0

!Set axial value
        v1(1)=v1(2)-rho_3d(2)*(v1(3)-v1(2))/(rho_3d(3)-rho_3d(2))

!Interpolate to external grid points inside R,Z domain
        iflag=0
        message=''
        CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,SHIFT_R,iflag,          &
     &                      message)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_SHAPE(7)/'//message
          GOTO 9999

        ENDIF

        IF(nr < nrho_r) THEN

          !Extrapolate to external grid points outside R,Z domain
          SHIFT_R(nr+1:nrho_r)=SHIFT_R(nr)

        ENDIF

      ENDIF

!-------------------------------------------------------------------------------
!Elongation
!-------------------------------------------------------------------------------
      IF(PRESENT(ELONG_R)) THEN

!Initialization
        ELONG_R(1:nrho_r)=0.0_rspec

!Elongation is total height over total width
        v1(2:nrho_3d)=(zmax(2:nrho_3d)-zmin(2:nrho_3d))                        &
     &                /(rmax(2:nrho_3d)-rmin(2:nrho_3d))

!Set axial value
        v1(1)=v1(2)-rho_3d(2)*(v1(3)-v1(2))/(rho_3d(3)-rho_3d(2))

!Interpolate to external grid points inside R,Z domain
        iflag=0
        message=''
        CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,ELONG_R,iflag,          &
     &                      message)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_SHAPE(8)/'//message
          GOTO 9999

        ENDIF

        IF(nr < nrho_r) THEN

          !Extrapolate to external grid points outside R,Z domain
          ELONG_R(nr+1:nrho_r)=ELONG_R(nr)

        ENDIF

      ENDIF

!-------------------------------------------------------------------------------
!Triangularity
!-------------------------------------------------------------------------------
      IF(PRESENT(TRIANG_R)) THEN

!Initialization
        TRIANG_R(1:nrho_r)=0.0_rspec

!Triangularity - average of upper and lower
        v1(2:nrho_3d)=((rmax(2:nrho_3d)+rmin(2:nrho_3d))                       &
     &                -(rbot(2:nrho_3d)+rtop(2:nrho_3d)))                      &
     &                /(rmax(2:nrho_3d)-rmin(2:nrho_3d))

!Set axial value
        v1(1)=0.0_rspec

!Interpolate to external grid points inside R,Z domain
        iflag=0
        message=''
        CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,TRIANG_R,iflag,         &
     &                      message)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_SHAPE(9)/'//message
          GOTO 9999

        ENDIF

        IF(nr < nrho_r) THEN

          !Extrapolate to external grid points outside R,Z domain
          TRIANG_R(nr+1:nrho_r)=TRIANG_R(nr)

        ENDIF

      ENDIF

!-------------------------------------------------------------------------------
!Maximum R
!-------------------------------------------------------------------------------
      IF(PRESENT(RMAX_R)) THEN

!Initialization
        RMAX_R(1:nrho_r)=0.0_rspec

!To cover odd-shaped plasmas, look for maxima near pi/4 and 7pi/4
        DO i=2,nrho_3d !Over radial nodes

          l_rminmax=.TRUE.
          r_flx(1)=rho_3d(i)
          r_flx(2)=0.25_rspec*z_pi
          iflag=0
          message=''
          CALL AJAX_MINMAX_RZ(l_rminmax,r_flx,r_cyl,iflag,message)

          !Check messages
          IF(iflag /= 0) THEN

            message='AJAX_SHAPE(10)/'//message
            GOTO 9999

          ENDIF

          IF(r_cyl(1) > rmax(i)) rmax(i)=r_cyl(1)
          r_flx(2)=1.75_rspec*z_pi
          iflag=0
          message=''
          CALL AJAX_MINMAX_RZ(l_rminmax,r_flx,r_cyl,iflag,message)

          !Check messages
          IF(iflag /= 0) THEN

            message='AJAX_SHAPE(11)/'//message
            GOTO 9999

          ENDIF

          IF(r_cyl(1) > rmax(i)) rmax(i)=r_cyl(1)

        ENDDO !Over radial nodes

!Set axial value
        rmax(1)=rmax(2)-rho_3d(2)*(rmax(3)-rmax(2))                            &
     &          /(rho_3d(3)-rho_3d(2))

!Interpolate to external grid points inside R,Z domain
        iflag=0
        message=''
        CALL LINEAR1_INTERP(nrho_3d,rho_3d,rmax,nr,rho_r,RMAX_R,iflag,         &
     &                      message)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_SHAPE(12)/'//message
          GOTO 9999

        ENDIF

      ENDIF

!-------------------------------------------------------------------------------
!Minimum R
!-------------------------------------------------------------------------------
      IF(PRESENT(RMIN_R)) THEN

!Initialization
        RMIN_R(1:nrho_r)=0.0_rspec

!To cover odd-shaped plasmas, look for minima near 3pi/4 and 5pi/4
        DO i=2,nrho_3d !Over radial nodes

          l_rminmax=.TRUE.
          r_flx(1)=rho_3d(i)
          r_flx(2)=0.75_rspec*z_pi
          iflag=0
          message=''
          CALL AJAX_MINMAX_RZ(l_rminmax,r_flx,r_cyl,iflag,message)

          !Check messages
          IF(iflag /= 0) THEN

            message='AJAX_SHAPE(13)/'//message
            GOTO 9999

          ENDIF

          IF(r_cyl(1) < rmin(i)) rmin(i)=r_cyl(1)
          r_flx(2)=1.25_rspec*z_pi
          iflag=0
          message=''
          CALL AJAX_MINMAX_RZ(l_rminmax,r_flx,r_cyl,iflag,message)

          !Check messages
          IF(iflag /= 0) THEN

            message='AJAX_SHAPE(14)/'//message
            GOTO 9999

          ENDIF

          IF(r_cyl(1) < rmin(i)) rmin(i)=r_cyl(1)

        ENDDO !Over radial nodes

!Set axial value
        rmin(1)=rmin(2)-rho_3d(2)*(rmin(3)-rmin(2))                            &
     &          /(rho_3d(3)-rho_3d(2))

!Interpolate to external grid points inside R,Z domain
        iflag=0
        message=''
        CALL LINEAR1_INTERP(nrho_3d,rho_3d,rmin,nr,rho_r,RMIN_R,iflag,         &
     &                      message)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_SHAPE(15)/'//message
          GOTO 9999

        ENDIF

      ENDIF

!-------------------------------------------------------------------------------
!Maximum Z
!-------------------------------------------------------------------------------
      IF(PRESENT(ZMAX_R)) THEN

!Initialization
        ZMAX_R(1:nrho_r)=0.0_rspec

!Set axial value
        zmax(1)=zmax(2)-rho_3d(2)*(zmax(3)-zmax(2))                            &
     &          /(rho_3d(3)-rho_3d(2))

!Interpolate to external grid points inside R,Z domain
        iflag=0
        message=''
        CALL LINEAR1_INTERP(nrho_3d,rho_3d,zmax,nr,rho_r,ZMAX_R,iflag,         &
     &                      message)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_SHAPE(16)/'//message
          GOTO 9999

        ENDIF

      ENDIF

!-------------------------------------------------------------------------------
!Minimum Z
!-------------------------------------------------------------------------------
      IF(PRESENT(ZMIN_R)) THEN

!Initialization
        ZMIN_R(1:nrho_r)=0.0_rspec

!Set axial value
        zmin(1)=zmin(2)-rho_3d(2)*(zmin(3)-zmin(2))                            &
     &          /(rho_3d(3)-rho_3d(2))

!Interpolate to external grid points inside R,Z domain
        iflag=0
        message=''
        CALL LINEAR1_INTERP(nrho_3d,rho_3d,zmin,nr,rho_r,ZMIN_R,iflag,         &
     &                      message)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_SHAPE(17)/'//message
          GOTO 9999

        ENDIF

      ENDIF

!-------------------------------------------------------------------------------
!R at maximum Z
!-------------------------------------------------------------------------------
      IF(PRESENT(RTOP_R)) THEN

!Initialization
        RTOP_R(1:nrho_r)=0.0_rspec

!Set axial value
        rtop(1)=rtop(2)-rho_3d(2)*(rtop(3)-rtop(2))                            &
     &          /(rho_3d(3)-rho_3d(2))

!Interpolate to external grid points inside R,Z domain
        iflag=0
        message=''
        CALL LINEAR1_INTERP(nrho_3d,rho_3d,rtop,nr,rho_r,RTOP_R,iflag,         &
     &                      message)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_SHAPE(18)/'//message
          GOTO 9999

        ENDIF

      ENDIF

!-------------------------------------------------------------------------------
!R at minimum Z
!-------------------------------------------------------------------------------
      IF(PRESENT(RBOT_R)) THEN

!Initialization
        RBOT_R(1:nrho_r)=0.0_rspec

!Set axial value
        rbot(1)=rbot(2)-rho_3d(2)*(rbot(3)-rbot(2))                            &
     &          /(rho_3d(3)-rho_3d(2))

!Interpolate to external grid points inside R,Z domain
        iflag=0
        message=''
        CALL LINEAR1_INTERP(nrho_3d,rho_3d,rbot,nr,rho_r,RBOT_R,iflag,         &
     &                      message)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_SHAPE(19)/'//message
          GOTO 9999

        ENDIF

      ENDIF

 9999 CONTINUE

      END SUBROUTINE AJAX_SHAPE

      SUBROUTINE AJAX_GLOBALS(iflag,message,                                   &
     &                        R000,PHITOT,RHOMAX,RHOMIN,SIGNBP,SIGNBT,         &
     &                        NPER,NRHO,NTHETA,NZETA)
!-------------------------------------------------------------------------------
!AJAX_GLOBALS gets global characteristics of the data
!References:
!  W.A.Houlberg, P.I.Strand 8/2001
!Output:
!  iflag               -error and warning flag [-]
!                      =-1 warning
!                      =0 none
!                      =1 error
!  message             -warning or error message [character]
!Optional output:
!  R000                -coefficient of (m,n)=(0,0) mode for R at axis [m]
!  PHITOT              -total toroidal flux [Wb]
!  RHOMAX              -outer radial boundary of R,Z domain [rho]
!  RHOMIN              -inner radial boundary of R,Z domain [rho]
!  SIGNBP              -sign of the poloidal magnetic field [-]
!                       positive is down on outside
!  SIGNBT              -sign of the toroidal magnetic field [-]
!                       positive is counterclockwise from top view
!  NPER                -number of field periods, =0 for axisymmetry [-]
!  NRHO                -number of radial nodes in internal data [-]
!  NTHETA              -number of poloidal nodes in internal data [-]
!  NZETA               -number of toroidal nodes in internal data [-]
!-------------------------------------------------------------------------------
      
!Declaration of output variables
      INTEGER, INTENT(OUT) ::                                                  &
     & iflag

      CHARACTER(len=*), INTENT(OUT) ::                                         &
     & message

!Declaration of optional output variables
      INTEGER, INTENT(OUT), OPTIONAL ::                                        &
     & NPER,NRHO,NTHETA,NZETA

      REAL(KIND=rspec), INTENT(OUT), OPTIONAL ::                               &
     & R000,                                                                   &
     & PHITOT,                                                                 &
     & RHOMIN,RHOMAX,SIGNBP,SIGNBT

      iflag=0
      message=''

      IF(PRESENT(R000)) R000=r000_3d
      IF(PRESENT(PHITOT)) PHITOT=phitot_3d
      IF(PRESENT(RHOMAX)) RHOMAX=rhomax_3d
      IF(PRESENT(RHOMIN)) RHOMIN=rhomin_3d
      IF(PRESENT(SIGNBP)) SIGNBP=signbp_3d
      IF(PRESENT(SIGNBT)) SIGNBT=signbt_3d
      IF(PRESENT(NPER)) NPER=nper_3d
      IF(PRESENT(NRHO)) NRHO=nrho_3d
      IF(PRESENT(NTHETA)) NTHETA=ntheta_3d
      IF(PRESENT(NZETA)) NZETA=nzeta_3d

      END SUBROUTINE AJAX_GLOBALS

      SUBROUTINE AJAX_MINMAX_RZ(l_rminmax,r_flx,r_cyl,iflag,message)
!-------------------------------------------------------------------------------
!AJAX_MINMAX_RZ gets maximum or minimum of R or Z
!References:
!  W.A.Houlberg, P.I.Strand 8/2001
!Input:
!  l_rminmax           -switch to select search for R or Z [logical]
!                      =.TRUE. find maximum or minimum R
!                      =.FALSE. find maximum or minimum Z
!Input/output:
!  r_flx(3)            -flux coordinates (rho,theta,zeta) [rho,rad,rad]
!Output:
!  r_cyl(3)            -cylindrical coordinates (R,phi,Z) [m,rad,m]
!  iflag               -error and warning flag [-]
!                      =-1 warning
!                      =0 none
!                      =1 error
!  message             -warning or error message [character]
!-------------------------------------------------------------------------------
!Declaration of input variables
      LOGICAL, INTENT(IN) ::                                                   &
     & l_rminmax
!Declaration of input/output variables
      REAL(KIND=rspec), INTENT(INOUT) ::                                       &
     & r_flx(3)

!Declaration of output variables
      INTEGER, INTENT(OUT) ::                                                  &
     & iflag

      CHARACTER(len=*), INTENT(OUT) ::                                         &
     & message

      REAL(KIND=rspec), INTENT(OUT) ::                                         &
     & r_cyl(3)

!Declaration of local variables
      INTEGER ::                                                               &
     & i,ig

      REAL(KIND=rspec) ::                                                      &
     & g_cyl0(6),g_cyl1(6),                                                    &
     & r_flx0(3),r_flx1(3)

!Set index for derivative
      IF(l_rminmax) THEN

        !d(R)/d(theta)
        ig=2

      ELSE

        !d(Z)/d(theta)
        ig=5

      ENDIF

!Initialize
      r_flx0(:)=0.0_rspec
      r_flx1(:)=0.0_rspec
      g_cyl0(:)=0.0_rspec
      g_cyl1(:)=0.0_rspec

      r_flx0(1:3)=r_flx(1:3)
      r_flx0(2)=r_flx0(2)-0.05_rspec*z_pi
      iflag=0
      message=''
      CALL AJAX_FLX2CYL(r_flx0,r_cyl,iflag,message,                            &
     &                  G_CYL=g_cyl0)

      !Check messages
      IF(iflag /= 0) THEN

        message='AJAX_MINMAX_RZ(1)/'//message
        IF(iflag > 0) GOTO 9999

      ENDIF

      r_flx1(1:3)=r_flx(1:3)
      r_flx1(2)=r_flx1(2)+0.05_rspec*z_pi
      iflag=0
      message=''
      CALL AJAX_FLX2CYL(r_flx1,r_cyl,iflag,message,                            &
     &                  G_CYL=g_cyl1)

      !Check messages
      IF(iflag /= 0) THEN

        message='AJAX_MINMAX_RZ(2)/'//message
        IF(iflag > 0) GOTO 9999

      ENDIF

      !Normally 1 or 2 iterations are required
      DO i=1,10 !Over iteration

        r_flx(2)=r_flx0(2)-g_cyl0(ig)*(r_flx1(2)-r_flx0(2))                    &
     &                               /(g_cyl1(ig)-g_cyl0(ig))
        r_flx0(1:3)=r_flx1(1:3)
        g_cyl0(:)=g_cyl1(:)
        r_flx1(1:3)=r_flx(1:3)
        iflag=0
        message=''
        CALL AJAX_FLX2CYL(r_flx,r_cyl,iflag,message,                           &
     &                    G_CYL=g_cyl1)

        !Check messages
        IF(iflag /= 0) THEN

          message='AJAX_MINMAX_RZ(3)/'//message
          IF(iflag > 0) GOTO 9999

        ENDIF

        IF(ABS(g_cyl1(ig)) < 1.0e-5_rspec*r000_3d/z_pi) EXIT

      ENDDO !Over iteration

 9999 CONTINUE

      END SUBROUTINE AJAX_MINMAX_RZ

      SUBROUTINE AJAX_INIT
!-------------------------------------------------------------------------------
!AJAX_INIT initializes or resets private data for the radial, poloidal and
!  toroidal grids, and allocates storage for the R, Z, and stream function
!  expansions
!References:
!  W.A.Houlberg, P.I.Strand 11/2001
!-------------------------------------------------------------------------------

!Declaration of local variables
      INTEGER ::                                                               &
     & i,                                                                      &
     & nset1,nset2,                                                            &
     & kmax

!-------------------------------------------------------------------------------
!Constants
!-------------------------------------------------------------------------------
!Machine and physical constants
      z_large=HUGE(1.0_rspec)
      z_pi=ACOS(-1.0_rspec)
      z_precision=EPSILON(1.0_rspec)
      z_small=TINY(1.0_rspec)
      z_mu0=4.0e-7_rspec*z_pi

!Resolution of rho
      rhores_3d=1.0e-3_rspec*rhomax_3d

!-------------------------------------------------------------------------------
!Radial grid and R,Z allocation
!-------------------------------------------------------------------------------
      IF(ALLOCATED(rho_3d)) THEN

!Check allocation
        !Radial dimension
        nset1=SIZE(r_3d,2)

        !Mode dimension
        nset2=SIZE(r_3d,3)

!If storage requirements have changed, reallocate
        IF(nset1 /= nrho_3d .OR.                                               &
     &     nset2 /= krz_3d) THEN

          !Reallocate radial grid and R,Z
          DEALLOCATE(rho_3d,                                                   &
     &               r_3d,                                                     &
     &               z_3d)
          ALLOCATE(rho_3d(nrho_3d),                                            &
     &              r_3d(4,nrho_3d,krz_3d),                                    &
     &             z_3d(4,nrho_3d,krz_3d))

          !Initialize
          rho_3d(:)=0.0_rspec
          r_3d(:,:,:)=0.0_rspec
          z_3d(:,:,:)=0.0_rspec

        ENDIF

      ELSE

        !Allocate radial grid and R,Z
        ALLOCATE(rho_3d(nrho_3d),                                              &
     &           r_3d(4,nrho_3d,krz_3d),                                       &
     &           z_3d(4,nrho_3d,krz_3d))

        !Initialization
        rho_3d(:)=0.0_rspec
        r_3d(:,:,:)=0.0_rspec
        z_3d(:,:,:)=0.0_rspec

      ENDIF

!Set radial grid
      rho_3d(:)=REAL((/ (i-1,i=1,nrho_3d) /),rspec)                            &
     &          /REAL(nrho_3d-1,rspec)*rhomax_3d

!-------------------------------------------------------------------------------
!Lambda allocation
!-------------------------------------------------------------------------------
!Check allocation
      IF(ALLOCATED(lam_3d)) THEN

        !Radial dimension
        nset1=SIZE(lam_3d,2)

        !Mode dimension
        nset2=SIZE(lam_3d,3)

!If storage requirements have changed, reallocate
        IF(nset1 /= nrho_3d .OR.                                               &
     &     nset2 /= klam_3d) THEN

          !Realloate lambda
          DEALLOCATE(lam_3d)
          ALLOCATE(lam_3d(4,nrho_3d,klam_3d))

          !Initialization
          lam_3d(:,:,:)=0.0_rspec

        ENDIF

      ELSE

          !Allocate lambda
          ALLOCATE(lam_3d(4,nrho_3d,klam_3d))

          !Initialization
          lam_3d(:,:,:)=0.0_rspec

      ENDIF

!-------------------------------------------------------------------------------
!Toroidal and poloidal modes
!-------------------------------------------------------------------------------
!Check allocation of m,n
      kmax=MAX(krz_3d,klam_3d)

      IF(ALLOCATED(m_3d)) THEN

        !Radial dimension
        nset1=SIZE(m_3d)

!If storage requirements have changed, reallocate
        IF(nset1 /= kmax) THEN

          !Reallocate toroidal and poloidal modes
          DEALLOCATE(m_3d,                                                     &
     &               n_3d,                                                     &
     &               mabs_3d)
          ALLOCATE(m_3d(kmax),                                                 &
     &             n_3d(kmax),                                                 &
     &             mabs_3d(kmax))

          !Initialization
          m_3d(:)=0
          n_3d(:)=0
          mabs_3d(:)=0

        ENDIF

      ELSE

        !Allocate toroidal and poloidal modes
        ALLOCATE(m_3d(kmax),                                                   &
     &           n_3d(kmax),                                                   &
     &           mabs_3d(kmax))

          !Initialization
          m_3d(:)=0
          n_3d(:)=0
          mabs_3d(:)=0

      ENDIF

!Set the poloidal grid spacing
      dtheta_3d=2.0_rspec*z_pi/(ntheta_3d-1)

!Check allocation of poloidal weighting
      IF(ALLOCATED(wtheta_3d)) THEN

        !Poloidal dimension
        nset1=SIZE(wtheta_3d)

        IF(nset1 /= ntheta_3d) THEN

          !Reallocate poloidal weighting
          DEALLOCATE(wtheta_3d)
          ALLOCATE(wtheta_3d(ntheta_3d))

          !Set poloidal weights
          wtheta_3d(1)=1.0_rspec
          wtheta_3d(ntheta_3d)=1.0_rspec

          DO i=2,ntheta_3d-1 !Over poloidal nodes

            IF(MOD(i,2) == 0) THEN

              !Even node
              wtheta_3d(i)=4.0_rspec

            ELSE

              !Odd mode
              wtheta_3d(i)=2.0_rspec

            ENDIF

          ENDDO !Over poloidal nodes

          wtheta_3d(:)=wtheta_3d(:)*2.0_rspec*z_pi/3.0_rspec                   &
     &                 /REAL(ntheta_3d-1,rspec)

        ENDIF

      ELSE

        !Allocate poloidal weighting
        ALLOCATE(wtheta_3d(ntheta_3d))

        !Set poloidal weights
        wtheta_3d(1)=1.0_rspec
        wtheta_3d(ntheta_3d)=1.0_rspec

        DO i=2,ntheta_3d-1 !Over poloidal nodes

          IF(MOD(i,2) == 0) THEN

            !Even node
            wtheta_3d(i)=4.0_rspec

          ELSE

            !Odd mode
            wtheta_3d(i)=2.0_rspec

          ENDIF

        ENDDO !Over poloidal nodes

        wtheta_3d(:)=wtheta_3d(:)*2.0_rspec*z_pi/3.0_rspec                     &
     &               /REAL(ntheta_3d-1,rspec)

      ENDIF

!Set the toroidal grid spacing
      IF(nzeta_3d == 1) THEN

        dzeta_3d=2.0_rspec*z_pi

      ELSE

        dzeta_3d=2.0_rspec*z_pi/REAL(nzeta_3d-1,rspec)                        &
     &           /REAL(nper_3d,rspec)

      ENDIF

!Check allocation of toroidal weighting
      IF(ALLOCATED(wzeta_3d)) THEN

        !Toroidal dimension
        nset1=SIZE(wzeta_3d)

        IF(nset1 /= nzeta_3d) THEN

          !Reallocate toroidal weighting
          DEALLOCATE(wzeta_3d)
          ALLOCATE(wzeta_3d(nzeta_3d))

          !Set toroidal weights
          wzeta_3d(1)=dzeta_3d

          IF(nzeta_3d > 1) THEN

            !Non-axisymmetric plasma
            wzeta_3d(1)=1.0_rspec
            wzeta_3d(nzeta_3d)=1.0_rspec

            DO i=2,nzeta_3d-1 !Over toroidal nodes

              IF(MOD(i,2) == 0) THEN

                !Even node
                wzeta_3d(i)=4.0_rspec

              ELSE

                !Odd mode
                wzeta_3d(i)=2.0_rspec

              ENDIF

            ENDDO !Over toroidal nodes

            wzeta_3d(:)=wzeta_3d(:)*REAL(nper_3d,rspec)*dzeta_3d               &
     &                  /3.0_rspec

          ENDIF

        ENDIF

      ELSE

        !Allocate toroidal weighting
        ALLOCATE(wzeta_3d(nzeta_3d))

        !Set toroidal weights
        wzeta_3d(1)=dzeta_3d

        IF(nzeta_3d > 1) THEN

          !Non-axisymmetric plasma
          wzeta_3d(1)=1.0_rspec
          wzeta_3d(nzeta_3d)=1.0_rspec

          DO i=2,nzeta_3d-1 !Over toroidal nodes

            IF(MOD(i,2) == 0) THEN

              !Even node
              wzeta_3d(i)=4.0_rspec

            ELSE

              !Odd mode
              wzeta_3d(i)=2.0_rspec

            ENDIF

          ENDDO !Over toroidal nodes

          wzeta_3d(:)=wzeta_3d(:)*REAL(nper_3d,rspec)*dzeta_3d/3.0_rspec

        ENDIF

      ENDIF

!-------------------------------------------------------------------------------
!Set logical switches for flux surface averaging
!-------------------------------------------------------------------------------
      l_fluxavb_3d=.FALSE.
      l_fluxavg_3d=.FALSE.

      END SUBROUTINE AJAX_INIT

      SUBROUTINE AJAX_INIT_FLUXAV_B
!-------------------------------------------------------------------------------
!AJAX_INIT_FLUXAV_B initializes magnetic field arrays
!References:
!  W.A.Houlberg, P.I.Strand 8/2001
!-------------------------------------------------------------------------------

!Declaration of local variables
      INTEGER ::                                                               &
     & i,j,k,                                                                  &
     & nset1,nset2,nset3

      REAL(KIND=rspec) ::                                                      &
     & br,bz,bzeta,                                                            &
     & phiprm

      REAL(KIND=rspec) ::                                                      &
     & r_flx(3)

!-------------------------------------------------------------------------------
!Load lambda derivatives for flux surface averaging
!-------------------------------------------------------------------------------
!Check allocation
        IF(ALLOCATED(eltheta_3d)) THEN

          !Radial dimension
          nset1=SIZE(eltheta_3d,1)

          !Poloidal dimension
          nset2=SIZE(eltheta_3d,2)

          !Toroidal dimension
          nset3=SIZE(eltheta_3d,3)

!If storage requirements have changed, reallocate
          IF(nset1 /= nrho_3d .OR.                                             &
     &       nset2 /= ntheta_3d .OR.                                           &
     &       nset3 /= nzeta_3d) THEN

            !Reallocate lambda derivatives
            DEALLOCATE(eltheta_3d,                                             &
     &                 elzeta_3d)
            ALLOCATE(eltheta_3d(1:nrho_3d,1:ntheta_3d,1:nzeta_3d),             &
     &               elzeta_3d(1:nrho_3d,1:ntheta_3d,1:nzeta_3d))

          ENDIF

        ELSE

          !Allocate lambda derivatives
          ALLOCATE(eltheta_3d(1:nrho_3d,1:ntheta_3d,1:nzeta_3d),               &
     &             elzeta_3d(1:nrho_3d,1:ntheta_3d,1:nzeta_3d))

        ENDIF

!Initialization
        eltheta_3d(:,:,:)=0.0_rspec
        elzeta_3d(:,:,:)=0.0_rspec
        r_flx(:)=0.0_rspec

!Fill arrays
        DO i=2,nrho_3d !Over radial nodes

          r_flx(1)=rho_3d(i)

          DO j=1,ntheta_3d !Over poloidal nodes

            r_flx(2)=REAL(j-1,rspec)*dtheta_3d

            DO k=1,nzeta_3d !Over toroidal nodes

              r_flx(3)=REAL(k-1,rspec)*dzeta_3d

              CALL AJAX_LAMBDA(r_flx,eltheta_3d(i,j,k),elzeta_3d(i,j,k))

            ENDDO !Over toroidal nodes

          ENDDO !Over poloidal nodes

        ENDDO !Over radial nodes

!-------------------------------------------------------------------------------
!Load B data for flux surface averaging
!-------------------------------------------------------------------------------
!Check allocation      
      IF(ALLOCATED(b_3d)) THEN

        !Radial dimension
        nset1=SIZE(b_3d,1)

        !Poloidal dimension
        nset2=SIZE(b_3d,2)

        !Toroidal dimension
        nset3=SIZE(b_3d,3)

!If storage requirements have changed, reallocate
        IF(nset1 /= nrho_3d .OR.                                               &
     &     nset2 /= ntheta_3d .OR.                                             &
     &     nset3 /= nzeta_3d) THEN

          !Reallocate B data
          DEALLOCATE(b_3d,                                                     &
     &               btheta_3d)
          ALLOCATE(b_3d(1:nrho_3d,1:ntheta_3d,1:nzeta_3d),                     &
     &             btheta_3d(1:nrho_3d,1:ntheta_3d,1:nzeta_3d))
        ENDIF

      ELSE

        !Allocate B data
        ALLOCATE(b_3d(1:nrho_3d,1:ntheta_3d,1:nzeta_3d),                       &
     &           btheta_3d(1:nrho_3d,1:ntheta_3d,1:nzeta_3d))

      ENDIF

!Initialization
      btheta_3d(:,:,:)=0.0_rspec
      b_3d(:,:,:)=0.0_rspec

!Evaluate 3D arrays
      DO i=2,nrho_3d !Over radial nodes

        phiprm=2.0_rspec*rho_3d(i)*phitot_3d/rhomax_3d**2

        DO j=1,ntheta_3d !Over poloidal nodes

          DO k=1,nzeta_3d !Over toroidal nodes

            btheta_3d(i,j,k)=(iotabar_3d(1,i)-elzeta_3d(i,j,k))
            bzeta=(1.0_rspec+eltheta_3d(i,j,k))
            br=gcyl_3d(2,i,j,k)*btheta_3d(i,j,k)+gcyl_3d(3,i,j,k)*bzeta
            bz=gcyl_3d(5,i,j,k)*btheta_3d(i,j,k)+gcyl_3d(6,i,j,k)*bzeta
            b_3d(i,j,k)=SQRT(br**2+bz**2+(rcyl_3d(i,j,k)*bzeta)**2)            &
     &                  /gsqrt_3d(i,j,k)
            btheta_3d(i,j,k)=btheta_3d(i,j,k)/gsqrt_3d(i,j,k)

          ENDDO !Over toroidal nodes

        ENDDO !Over poloidal nodes
        btheta_3d(i,:,:)=phiprm*btheta_3d(i,:,:)
        b_3d(i,:,:)=ABS(phiprm)*b_3d(i,:,:)

      ENDDO !Over radial nodes

!Poloidal values at origin are linear extrapolation of averages at 2 and 3
      DO k=1,nzeta_3d

        btheta_3d(1,1:ntheta_3d,k)=(SUM(btheta_3d(2,1:ntheta_3d-1,k))          &
     &                              *rho_3d(3)                                 &
     &                             -SUM(btheta_3d(3,1:ntheta_3d-1,k))          &
     &                              *rho_3d(2))                                &
     &                             /REAL(ntheta_3d-1,rspec)                    &
     &                             /(rho_3d(3)-rho_3d(2))
        b_3d(1,1:ntheta_3d,k)=(SUM(b_3d(2,1:ntheta_3d-1,k))                    &
     &                         *rho_3d(3)                                      &
     &                        -SUM(b_3d(3,1:ntheta_3d-1,k))                    &
     &                         *rho_3d(2))                                     &
     &                        /REAL(ntheta_3d-1,rspec)                         &
     &                        /(rho_3d(3)-rho_3d(2))

      ENDDO

!Divide by 2*pi
      btheta_3d(:,:,:)=0.5_rspec*btheta_3d(:,:,:)/z_pi
      b_3d(:,:,:)=0.5_rspec*b_3d(:,:,:)/z_pi

!-------------------------------------------------------------------------------
!Set initialization flag
!-------------------------------------------------------------------------------
      l_fluxavb_3d=.TRUE.

      END SUBROUTINE AJAX_INIT_FLUXAV_B

      SUBROUTINE AJAX_INIT_FLUXAV_G(iflag,message)
!-------------------------------------------------------------------------------
!AJAX_INIT_FLUXAV_G initializes geometry arrays
!References:
!  W.A.Houlberg, P.I.Strand 8/2001
!Output:
!  iflag               -error and warning flag [-]
!                      =-1 warning
!                      =0 none
!                      =1 error
!  message             -warning or error message [character]
!-------------------------------------------------------------------------------

!Declaration of output variables
      INTEGER, INTENT(OUT) ::                                                  &
     & iflag

      CHARACTER(len=*), INTENT(OUT) ::                                         &
     & message

!Declaration of local variables
      INTEGER ::                                                               &
     & i,j,k,                                                                  &
     & nset1,nset2,nset3

      REAL(KIND=rspec) ::                                                      &
     & r_cyl(3),r_flx(3)

!-------------------------------------------------------------------------------
!Load R,Z and other geometric data for flux surface averaging
!-------------------------------------------------------------------------------
!Check allocation
      IF(ALLOCATED(rcyl_3d)) THEN

        !Radial dimension
        nset1=SIZE(rcyl_3d,1)

        !Poloidal dimension
        nset2=SIZE(rcyl_3d,2)

        !Toroidal dimension
        nset3=SIZE(rcyl_3d,3)

!If storage requirements have changed, reallocate
        IF(nset1 /= nrho_3d .OR.                                               &
     &     nset2 /= ntheta_3d .OR.                                             &
     &     nset3 /= nzeta_3d) THEN

          !Reallocate R,Z and other geometric data
          DEALLOCATE(rcyl_3d,                                                  &
     &               zcyl_3d,                                                  &
     &               gsqrt_3d,                                                 &
     &               gcyl_3d,                                                  &
     &               vp_3d,                                                    &
     &               gt_3d,                                                    &
     &               az_3d)
          ALLOCATE(rcyl_3d(1:nrho_3d,1:ntheta_3d,1:nzeta_3d),                  &
     &             zcyl_3d(1:nrho_3d,1:ntheta_3d,1:nzeta_3d),                  &
     &             gsqrt_3d(1:nrho_3d,1:ntheta_3d,1:nzeta_3d),                 &
     &             gcyl_3d(1:6,1:nrho_3d,1:ntheta_3d,1:nzeta_3d),              &
     &             vp_3d(1:nrho_3d),                                           &
     &             gt_3d(1:ntheta_3d),                                         &
     &             az_3d(1:nzeta_3d))

        ENDIF

      ELSE

        !Allocate R,Z and other geometric data
        ALLOCATE(rcyl_3d(1:nrho_3d,1:ntheta_3d,1:nzeta_3d),                    &
     &           zcyl_3d(1:nrho_3d,1:ntheta_3d,1:nzeta_3d),                    &
     &           gsqrt_3d(1:nrho_3d,1:ntheta_3d,1:nzeta_3d),                   &
     &           gcyl_3d(1:6,1:nrho_3d,1:ntheta_3d,1:nzeta_3d),                &
     &           vp_3d(1:nrho_3d),                                             &
     &           gt_3d(1:ntheta_3d),                                           &
     &           az_3d(1:nzeta_3d))

      ENDIF

!Initialization
      rcyl_3d(:,:,:)=0.0_rspec
      zcyl_3d(:,:,:)=0.0_rspec
      gsqrt_3d(:,:,:)=0.0_rspec
      gcyl_3d(:,:,:,:)=0.0_rspec
      vp_3d(:)=0.0_rspec
      gt_3d(:)=0.0_rspec
      az_3d(:)=0.0_rspec
      r_flx(:)=0.0_rspec
      r_cyl(:)=0.0_rspec

!Evaluate 3D cylindrical coordinates, metrics, and Jacobian
      DO i=2,nrho_3d !Over radial nodes

        DO j=1,ntheta_3d !Over poloidal nodes

          DO k=1,nzeta_3d !Over toroidal nodes

            r_flx(1)=rho_3d(i)
            r_flx(2)=REAL(j-1,rspec)*dtheta_3d
            r_flx(3)=REAL(k-1,rspec)*dzeta_3d
            iflag=0
            message=''
            CALL AJAX_FLX2CYL(r_flx,r_cyl,iflag,message,                       &
     &                        G_CYL=gcyl_3d(1:6,i,j,k),                        &
     &                        GSQRT=gsqrt_3d(i,j,k))

            !Check messages
            IF(iflag /= 0) THEN

              message='AJAX_INIT_FLUXAV_G/'//message
              IF(iflag > 0) GOTO 9999

            ENDIF

            rcyl_3d(i,j,k)=r_cyl(1)
            zcyl_3d(i,j,k)=r_cyl(3)

          ENDDO !Over toroidal nodes

        ENDDO !Over poloidal nodes

      ENDDO !Over radial nodes

!Calculate V' on internal grid 
      DO i=2,nrho_3d !Over radial nodes

        DO k=1,nzeta_3d !Over toroidal nodes

          gt_3d(1:ntheta_3d)=gsqrt_3d(i,1:ntheta_3d,k)
          az_3d(k)=SUM(wtheta_3d*gt_3d) !Over poloidal nodes

        ENDDO !Over toroidal nodes

        vp_3d(i)=SUM(wzeta_3d*az_3d) !Over toroidal nodes

      ENDDO !Over radial nodes

!-------------------------------------------------------------------------------
!Set logical switches for flux surface averaging
!-------------------------------------------------------------------------------
      l_fluxavg_3d=.TRUE.

 9999 CONTINUE

      END SUBROUTINE AJAX_INIT_FLUXAV_G

      SUBROUTINE AJAX_INIT_LAMBDA
!-------------------------------------------------------------------------------
!AJAX_INIT_LAMBDA calculates the magnetic stream function expansion coefficients
!  for an axisymmetric plasma
!References:
!  W.A.Houlberg, P.I.Strand 8/2001
!-------------------------------------------------------------------------------

!Declaration of local variables
      INTEGER ::                                                               &
     & i,j,k,                                                                  &
     & k_bc1=3,k_bcn=0

      REAL(KIND=rspec) ::                                                      &
     & theta,                                                                  &
     & g1(1:ntheta_3d),g2(1:ntheta_3d)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
      lam_3d(:,:,:)=0.0_rspec
      g1(:)=0.0_rspec
      g2(:)=0.0_rspec

!-------------------------------------------------------------------------------
!Calculate lambda values
!-------------------------------------------------------------------------------
!Loop over harmonics
      DO k=1,klam_3d !Over modes

        IF(k /= km0n0_3d) THEN

!Loop over radial grid
          DO i=2,nrho_3d !Over radial nodes

!Loop over theta grid
            DO j=1,ntheta_3d !Over poloidal nodes

              g1(j)=gsqrt_3d(i,j,1)/rcyl_3d(i,j,1)**2
              theta=(j-1)*dtheta_3d
              g2(j)=g1(j)*COS(m_3d(k)*theta)

            ENDDO !Over poloidal nodes

!Theta averages to find lambda coefficients
            lam_3d(1,i,k)=2.0_rspec*SUM(wtheta_3d*g2)/SUM(wtheta_3d*g1)        &
     &                    /REAL(m_3d(k),rspec) !Over poloidal nodes

          ENDDO !Over radial nodes

          lam_3d(1,2:nrho_3d,k)=lam_3d(1,2:nrho_3d,k)                          &
     &                          /rho_3d(2:nrho_3d)**mabs_3d(k)

!Make a parabolic fit to axis
          lam_3d(1,1,k)=(lam_3d(1,2,k)*rho_3d(3)**2                            &
     &                  -lam_3d(1,3,k)*rho_3d(2)**2)                           &
     &                  /(rho_3d(3)**2-rho_3d(2)**2)
 
!-------------------------------------------------------------------------------
!Spline the Lmn coefficients for internal storage (not-a-knot edge BC)
!-------------------------------------------------------------------------------
          CALL SPLINE1_FIT(nrho_3d,rho_3d,lam_3d(:,:,k),                       &
     &                     K_BC1=k_bc1,                                        &
     &                     K_BCN=k_bcn)

        ENDIF

      ENDDO !Over modes

      END SUBROUTINE AJAX_INIT_LAMBDA

      SUBROUTINE AJAX_LOAD_LAMBDA(k_grid,nr_lam,nk_lam,rhorz,rho_lam,          &
     &                            lam,iflag,message)
!-------------------------------------------------------------------------------
!AJAX_LOAD_LAMBDA loads the magnetic stream function
!References:
!  W.A.Houlberg, P.I.Strand 8/2001
!Input:
!  k_grid              -option designating type of input radial grid [-]
!                      =0 proportional to sqrt(toroidal flux)
!                      =1 proportional to toroidal flux
!                      =else not allowed
!  nr_lam              -no. of radial nodes in the input lambda [-]
!  nk_lam              -no. of poloidal & toroidal modes in the input lambda [-]
!  rhorz               -max rho of the R,Z data [arb]
!  rho_lam(nr_lam)     -radial nodes in the input lambda [arb]
!  lam(nr_lam,nk_lam)  -expansion coeffs for lambda [-]
!Output:
!  iflag               -error and warning flag [-]
!                      =-1 warning
!                      =0 none
!                      =1 error
!  message             -warning or error message [character]
!Comments:
!  The radial grid is assumed to be the same as for the R,Z data, so the scale
!    factor, rhorz, and grid type, k_grid are needed input
!-------------------------------------------------------------------------------

!Declaration of input variables
      INTEGER, INTENT(IN) ::                                                   &
     & k_grid,                                                                 &
     & nr_lam,nk_lam

      REAL(KIND=rspec), INTENT(IN) ::                                          &
     & rhorz,                                                                  &
     & rho_lam(:),                                                             &
     & lam(:,:)

!Declaration of output variables
      INTEGER, INTENT(OUT) ::                                                  &
     & iflag

      CHARACTER(len=*), INTENT(OUT) ::                                         &
     & message

!Declaration of local variables
      LOGICAL ::                                                               &
     & l_edge

      INTEGER ::                                                               &
     & j,k,                                                                    &
     & nset1,nset2,                                                            &
     & k_vopt(3)=(/1,0,0/),                                                    &
     & k_bc1=3,k_bcn=0

      REAL(KIND=rspec), ALLOCATABLE ::                                         &
     & rho(:),                                                                 &
     & lmn(:,:),                                                               &
     & fspl(:,:),                                                              &
     & values(:,:)

!-------------------------------------------------------------------------------
!Load lambda expansion data
!-------------------------------------------------------------------------------
!Initialization
      lam_3d(:,:,:)=0.0_rspec

!Copy temporary radial grid and expansion coefficients
      nset2=nk_lam

!Check whether inner boundary is at axis
      IF(rho_lam(1) > 10.0_rspec*z_precision*rhorz) THEN

        !Radial node needs to be added at axis
        !Check whether boundary node is present
        IF(rho_lam(nr_lam) < (1.0_rspec+10.0_rspec*z_precision)*rhorz)         &
     &    THEN

          !Add radial nodes at center and edge
          nset1=nr_lam+2
          ALLOCATE(rho(1:nset1),                                               &
     &             lmn(1:nset1,1:nset2))
          rho(1)=0.0_rspec
          rho(2:nset1-1)=rho_lam(1:nr_lam)/rhorz
          rho(nset1)=1.0_rspec
          lmn(2:nset1-1,1:nset2)=lam(1:nr_lam,1:nset2)
          l_edge=.TRUE.

        ELSE

          !Add radial node at center
          nset1=nr_lam+1
          ALLOCATE(rho(1:nset1),                                               &
     &             lmn(1:nset1,1:nset2))
          rho(1)=0.0_rspec
          rho(2:nset1)=rho_lam(1:nr_lam)/rhorz
          lmn(2:nset1,1:nset2)=lam(1:nr_lam,1:nset2)
          l_edge=.FALSE.

        ENDIF

      ELSE

        !Radial node is present at axis
        !Check whether boundary node is present
        IF(rho_lam(nr_lam) < (1.0_rspec+10.0_rspec*z_precision)*rhorz)         &
     &    THEN

          !Add radial node at edge
          nset1=nr_lam+1
          ALLOCATE(rho(1:nset1),                                               &
     &             lmn(1:nset1,1:nset2))
          rho(1:nset1-1)=rho_lam(1:nset1-1)/rhorz
          rho(nset1)=1.0_rspec
          lmn(1:nset1-1,1:nset2)=lam(1:nset1-1,1:nset2)
          l_edge=.TRUE.

        ELSE

          !Use input grid
          nset1=nr_lam
          ALLOCATE(rho(1:nset1),                                               &
     &             lmn(1:nset1,1:nset2))
          rho(1:nset1)=rho_lam(1:nset1)/rhorz
          lmn(1:nset1,1:nset2)=lam(1:nset1,1:nset2)
          l_edge=.FALSE.

        ENDIF

      ENDIF

!Convert rho to sqrt(toroidal flux) if necesssary and scale
      IF(k_grid == 1) THEN

        !~toroidal flux
        rho(:)=SQRT(rho(:))*rhomax_3d

      ELSE

        !~sqrt(toroidal flux)
        rho(:)=rho(:)*rhomax_3d

      ENDIF

!Allocate temporary work arrays
      ALLOCATE(fspl(4,nset1),                                                  &
     &         values(3,nrho_3d))
      fspl(:,:)=0.0_rspec
      values(:,:)=0.0_rspec

!Normalize the expansion coeffs to rho**m
      DO k=1,klam_3d !Over modes

        IF(l_edge) THEN

          !Extrapolate to edge
          DO j=2,nset1-1 !Over radial nodes

            lmn(j,k)=lmn(j,k)/rho(j)**mabs_3d(k)

          ENDDO !Over radial nodes

          lmn(nset1,k)=(lmn(nset1-1,k)*(rho(nset1)-rho(nset1-2))               &
     &                 -lmn(nset1-2,k)*(rho(nset1)-rho(nset1-1)))              &
     &                 /(rho(nset1-1)-rho(nset1-2))
        ELSE

          DO j=2,nset1 !Over radial nodes

            lmn(j,k)=lmn(j,k)/rho(j)**mabs_3d(k)

          ENDDO !Over radial nodes

        ENDIF

!Make a parabolic fit to axis
        lmn(1,k)=(lmn(2,k)*rho(3)**2-lmn(3,k)*rho(2)**2)                       &
     &           /(rho(3)**2-rho(2)**2)
 
!Map lambda_mn to internal radial grid
        fspl(1,1:nset1)=lmn(1:nset1,k)
        iflag=0
        message=''
        CALL SPLINE1_INTERP(k_vopt,nset1,rho,fspl,nrho_3d,rho_3d,values,       &
     &                      iflag,message,                                     &
     &                      K_BC1=k_bc1,                                       &
     &                      K_BCN=k_bcn)

        !Check messages
        IF(iflag > 0) THEN

          message='AJAX_LOAD_LAM/'//message
          GOTO 9999

        ENDIF

          !Respline the Lmn coefficients for internal storage (not-a-knot now)
          lam_3d(1,1:nrho_3d,k)=values(1,1:nrho_3d)
          CALL SPLINE1_FIT(nrho_3d,rho_3d,lam_3d(:,:,k),                       &
     &                     K_BC1=k_bc1,                                        &
     &                     K_BCN=k_bcn)

        ENDDO !Over modes

 9999 CONTINUE

!-------------------------------------------------------------------------------
!Cleanup
!-------------------------------------------------------------------------------

      IF(ALLOCATED(rho)) THEN

        !Deallocate radial grid and lambda arrays
        DEALLOCATE(rho,                                                        &
     &             lmn)

      ENDIF

      IF(ALLOCATED(fspl)) THEN

        !Deallocate spline arrays
        DEALLOCATE(fspl,                                                       &
     &             values)

      ENDIF

      END SUBROUTINE AJAX_LOAD_LAMBDA

      END MODULE ajax_mod

      MODULE LINEAR1_MOD
!-------------------------------------------------------------------------------
!LINEAR1-LINEAR interpolation in 1d
!
!LINEAR1_MOD is an F90 module of linear interpolating routines in 1d
!
!References:
!
!  W.A.Houlberg, P.I. Strand 6/2001
!
!Contains PUBLIC routines:
!
!  LINEAR1_INTERP      -interpolate from one grid to another
!  LINEAR1_INTEG       -integrate the linear fit
!
!Comments:
!
!  Linear interpolation routines are C0 (only f is continuous)
!
!  The modernization of the code structure into an F90 module takes advantage of
!    some of the more attractive features of F90:
!    -use of KIND for precision declarations
!    -optional arguments for I/O
!    -generic names for all intrinsic functions
!    -compilation using either free or fixed form
!    -no common blocks or other deprecated Fortran features
!    -dynamic and automatic alocation of variables
!    -array syntax for vector operations
!-------------------------------------------------------------------------------
      USE SPEC_KIND_MOD
      IMPLICIT NONE

!-------------------------------------------------------------------------------
! Public procedures
!-------------------------------------------------------------------------------
      CONTAINS

      SUBROUTINE LINEAR1_INTERP(n0,x0,y0,n1,x1,y1,iflag,message)
!-------------------------------------------------------------------------------
!W_LIN_INTERP performs a linear interpolation from the x0 mesh to the x1 mesh
!References:
!  W.A.Houlberg, P.I. Strand 6/2001
!Input:
!  n0                  -number of source abscissas [-]
!  x0(n0)              -source abscissas (in increasing order) [arb]
!  y0(n0)              -source values [arb]
!  n1                  -number of target abscissas [-]
!  x1(n1)              -target abscissas [arb]
!Output:
!  y1(n1)              -target values
!  iflag               -error and warning flag
!                      =-1 warning
!                      =0 no warnings or errors
!                      =1 error
!  message             -warning or error message (character)
!-------------------------------------------------------------------------------

!Declaration of input variables
      INTEGER, INTENT(IN) ::                                                   &
     &  n0,n1

      REAL(KIND=rspec), INTENT(IN) ::                                          &
     &  x0(:),x1(:),                                                           &
     &  y0(:)

!Declaration of output variables
      INTEGER, INTENT(OUT) ::                                                  &
     & iflag
     
      CHARACTER(len=*), INTENT(OUT) ::                                         &
     & message

      REAL(KIND=rspec), INTENT(OUT) ::                                         &
     & y1(:)

!Declaration of local variables
      INTEGER ::                                                               &
     & i,il

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!If no target values are requested, return
      IF(n1 <= 0) THEN

        iflag=-1
        message='LINEAR1_INTERP/WARNING(1):no points in output array'
        GOTO 9999

      ENDIF

!Make sure there are at least two nodes in the input grid
      IF(n0 < 2) THEN

        iflag=1
        message='LINEAR1_INTERP/ERROR(2):<2 points in source array'
        GOTO 9999

      ENDIF

!Set starting index of source grid
      il=1

!-------------------------------------------------------------------------------
!Interpolate from x0 to x1
!-------------------------------------------------------------------------------
      DO i=1,n1 !Over index of target grid

   10   IF(x1(i) < x0(1)) THEN

          !Target is below data range, use innermost data value
          y1(i)=y0(1)
          iflag=-1
          message='LINEAR1_INTERP(3)/WARNING:x<x(1), use end point'

        ELSEIF(x1(i) == x0(1)) THEN

          !Target and source nodes coincide
          y1(i)=y0(1)

        ELSEIF(x1(i) > x0(il+1)) THEN

          !Beyond next source node
          !Step x0 grid forward and loop
          IF(il < n0-1) THEN

            il=il+1
            GOTO 10

          ELSE

            !Target is above data range, set to last value
            y1(i)=y0(n0)
            iflag=-1
            message='LINEAR1_INTERP(4)/WARNING:x>x(n0), use end point'

          ENDIF

        ELSE

          !Between the proper set of source nodes, interpolate
          y1(i)=y0(il)                                                         &
     &          +(y0(il+1)-y0(il))*(x1(i)-x0(il))/(x0(il+1)-x0(il))

        ENDIF

      ENDDO !Over index of target grid

 9999 CONTINUE
      
      END SUBROUTINE LINEAR1_INTERP

      SUBROUTINE LINEAR1_INTEG(k_order,n0,x0,f,n1,x1,value,iflag,              &
     &                         message)
!-------------------------------------------------------------------------------
!LINEAR1_INTEG evaluates the integral f(x)*x**k_order, where f(x) is a linear
!  function and k_order is 0 or 1
!References:
!  W.A.Houlberg, P.I. Strand 6/2001
!Input:
!  k_order             -exponent of weighting factor (s(x)*x**k_order) [-]
!                      =0 or 1 only
!  n0                  -number of source nodes [-]
!  x0(n0)              -source abcissas [arb]
!  f(n0)               -source ordinates [arb]
!  n1                  -number of output nodes (may be 1) [-]
!  x1(n1)              -output abcissas at which the integral is wanted [arb]
!Output:
!  value(n1)           -integral of f(x)*x**k_order from x0(1) to x1(i) [arb]
!  iflag               -error and warning flag [-]
!                      =-1 warning
!                      =0 no warnings or errors
!                      =1 error
!  message             -warning or error message (character)
!-------------------------------------------------------------------------------

!Declaration of input variables
      INTEGER, INTENT(IN) ::                                                   &
     & k_order,                                                                &
     & n0,n1

      REAL(KIND=rspec), INTENT(IN) ::                                          &
     & x0(:),x1(:),                                                            &
     & f(:)

!Declaration of output variables
      INTEGER, INTENT(OUT) ::                                                  &
     & iflag
     
      CHARACTER(len=*), INTENT(OUT) ::                                         &
     & message

      REAL(KIND=rspec), INTENT(OUT) ::                                         &
     & value(:)

!Declaration of local variables
      INTEGER ::                                                               &
     & i,j,                                                                    &
     & ido,jdo

      REAL(KIND=rspec) ::                                                      &
     & add,                                                                    &
     & dx,f2,                                                                  &
     & sum,                                                                    &
     & xnew,xold

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!If no target values are requested, return
      IF(n1 <= 0) GOTO 9999

      value(1:n1)=0.0_rspec

!Integral value
      sum=0.0_rspec

!Source node, value and derivative
      j=1
      xold=x0(j)
      f2=(f(2)-f(1))/(x0(2)-x0(1))

!-------------------------------------------------------------------------------
!Integrate over target nodes
!-------------------------------------------------------------------------------
      DO i=1,n1 !Over target nodes

!Find dx to next source or target nodes
        ido=0

        DO WHILE(ido == 0) !Over source nodes

          IF(x1(i) < x0(j+1)) THEN

            !Hit target nodes
            xnew=x1(i)
            ido=1
            jdo=0

          ELSEIF(x1(i) == x0(j+1)) THEN

            !Hit both source and target nodes
            xnew=x1(i)
            ido=1
            jdo=1

          ELSEIF(x1(i) > x0(j+1)) THEN

            !Hit source nodes
            xnew=x0(j+1)
            ido=0
            jdo=1

          ENDIF

!Integrate over dx
          dx=xnew-xold

          IF(k_order == 1) THEN

            !Integral x s(x)dx
            add=dx*(xold*f(j)+dx*((xold*f2+f(j))/2.0_rspec))

          ELSE

            !Integral s(x)dx
            add=dx*(f(j)+dx*f2/2.0_rspec)

          ENDIF

!Add increment and update endpoint
          sum=sum+add
          xold=xnew

          IF(jdo == 1) THEN

            !Increment source node and derivative
            j=j+1
            IF(j < n0) f2=(f(j+1)-f(j))/(x0(j+1)-x0(j))

          ENDIF

          IF(j > n0) THEN

            !Target node is out of range
            iflag=1
            message='SPLINE1_INTEG/ERROR:target node out of range'
            GOTO 9999

          ENDIF

!Set integral value
          value(i)=sum

        ENDDO !Over source nodes

      ENDDO !Over target nodes
   
 9999 CONTINUE

      END SUBROUTINE LINEAR1_INTEG

      END MODULE LINEAR1_MOD


      MODULE SPLINE1_MOD
!-------------------------------------------------------------------------------
!SPLINE1-SPLINE interpolation in 1d
!
!SPLINE1_MOD is an F90 module of cubic interpolating spline routines in 1d
!
!References:
!
!  Forsythe, Malcolm, Moler, Computer Methods for Mathematical
!    Computations, Prentice-Hall, 1977, p.76
!  Engeln-Muellges, Uhlig, Numerical Algorithms with Fortran, Springer,
!    1996, p.251
!  W.A.Houlberg, P.I.Strand, D.McCune 8/2001
!
!Contains PUBLIC routines:
!
!  SPLINE1_FIT         -get the coefficients
!  SPLINE1_EVAL        -evaluate the spline
!  SPLINE1_INTERP      -interpolate from one grid to another
!  SPLINE1_INTEG       -integrate the spline fit
!
!Contains PRIVATE routine:
!
!  SPLINE1_SEARCH-find the indices that bracket an abscissa value
!
!Comments:
!
!  Spline interpolation routines are C2 (f, f', and f'' are continuous)
!
!  The modernization of the code structure into an F90 module takes advantage of
!    some of the more attractive features of F90:
!    -use of KIND for precision declarations
!    -optional arguments for I/O
!    -generic names for all intrinsic functions
!    -compilation using either free or fixed form
!    -no common blocks or other deprecated Fortran features
!    -dynamic and automatic alocation of variables
!    -array syntax for vector operations
!-------------------------------------------------------------------------------
      USE SPEC_KIND_MOD
      IMPLICIT NONE

!-------------------------------------------------------------------------------
! Private procedures
!-------------------------------------------------------------------------------
      PRIVATE ::                                                               &
     & SPLINE1_SEARCH

!-------------------------------------------------------------------------------
! Public procedures
!-------------------------------------------------------------------------------
      CONTAINS

      SUBROUTINE SPLINE1_FIT(n,x,f,                                            &
     &                       K_BC1,K_BCN)
!-------------------------------------------------------------------------------
!SPLINE1_FIT gets the coefficients for a 1d cubic interpolating spline
!References:
!  Forsythe, Malcolm, Moler, Computer Methods for Mathematical
!    Computations, Prentice-Hall, 1977, p.76
!  Engeln-Muellges, Uhlig, Numerical Algorithms with Fortran, Springer,
!    1996, p.251
!  W.A.Houlberg, P.I.Strand, D.McCune 8/2001
!Input:
!  n                   -number of data points or knots [>=2]
!  x(n)                -abscissas of the knots in increasing order [arb]
!  f(1,n)              -ordinates of the knots [arb]
!Input/Output:
!  f(2,1)              -input value of s'(x1) for K_BC1=1 [arb]
!  f(2,n)              -input value of s'(xn) for K_BCN=1 [arb]
!  f(3,1)              -input value of s''(x1) for K_BC1=2 [arb]
!  f(3,n)              -input value of s''(xn) for K_BCN=2 [arb]
!Output:
!  f(2,i)              =s'(x(i)) [arb]
!  f(3,i)              =s''(x(i)) [arb]
!  f(4,i)              =s'''(x(i)) [arb]
!Optional input:
!  K_BC1               -option for BC at x(1) [-]
!                      =-1 periodic, ignore K_BCN
!                      =0 not-a-knot (default)
!                      =1 s'(x1) = input value of f(2,1)
!                      =2 s''(x1) = input value of f(3,1)
!                      =3 s'(x1) = 0.0
!                      =4 s''(x1) = 0.0
!                      =5 match first derivative to first 2 points
!                      =6 match second derivative to first 3 points
!                      =7 match third derivative to first 4 points
!                      =else use not-a-knot
!  K_BCN               -option for boundary condition at x(n) [-]
!                      =0 not-a-knot (default)
!                      =1 s'(x1) = input value of f(2,1)
!                      =2 s''(x1) = input value of f(3,1)
!                      =3 s'(x1) = 0.0
!                      =4 s''(x1) = 0.0
!                      =5 match first derivative to first 2 points
!                      =6 match second derivative to first 3 points
!                      =7 match third derivative to first 4 points
!                      =else use knot-a-knot
!Comments:
!  For x(i).le.x.le.x(i+1)
!    s(x) = f(1,i) + f(2,i)*(x-x(i)) + f(3,i)*(x-x(i))**2/2!
!          +f(4,i)*(x-x(i))**3/3!
!  The cubic spline is twice differentiable (C2)
!  The BCs default to not-a-knot conditions, K_BC1=0 and K_BCN=0
!-------------------------------------------------------------------------------

!Declaration of input variables
      INTEGER, INTENT(IN) ::                                                   &
     & n

      REAL(KIND=rspec), INTENT(IN) ::                                          &
     & x(:)

!Declaration of input/output variables
      REAL(KIND=rspec), INTENT(INOUT) ::                                       &
     & f(:,:)

!Declaration of optional input variables
      INTEGER, INTENT(IN), OPTIONAL ::                                         &
     & K_BC1,K_BCN

!Declaration of local variables
      INTEGER ::                                                               &
     & i,ib,                                                                   &
     & imax,imin,                                                              &
     & kbc1,kbcn

      REAL(KIND=rspec) ::                                                      &
     & a1,an,                                                                  &
     & b1,bn,                                                                  &
     & q,t,                                                                    &
     & hn,                                                                     &
     & wk(1:n)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Set left BC option, and leftmost node
      !Default to not-a-knot
      kbc1=0

      !Otherwise use requested value
      IF(PRESENT(K_BC1)) THEN

        IF(K_BC1 >= -1 .AND.                                                   &
     &     K_BC1 <= 7) kbc1=K_BC1

      ENDIF

      !Include first node for all but not-a-knot
      imin=1

      !Not-a-knot condition removes first node
      IF(kbc1 == 0) imin=2

!Set left BC values
      !Default for not-a-knot
      a1=0.0_rspec
      b1=0.0_rspec

      IF(kbc1 == 1) THEN

        !First derivative specified
        a1=f(2,1)

      ELSEIF(kbc1 == 2) THEN

        !Second derivative specified
        b1=f(3,1)

      ELSEIF(kbc1 == 5) THEN

        !Match first derivative to first two points
        a1=(f(1,2)-f(1,1))/(x(2)-x(1))

      ELSEIF(kbc1 == 6) THEN

        !Match second derivative to first three points
        b1=2.0_rspec*((f(1,3)-f(1,2))/(x(3)-x(2))                              &
     &     -(f(1,2)-f(1,1))/(x(2)-x(1)))/(x(3)-x(1))

      ENDIF

!Set right BC option, and rightmost node
      !Default to not-a-knot
      kbcn=0

      !Otherwise use requested value
      IF(PRESENT(K_BCN)) THEN

        IF(K_BCN >= -1 .AND.                                                   &
     &     K_BCN <= 7) kbcn=K_BCN

      ENDIF

      !Include last node for all but not-a-knot
      imax=n

      !Not-a-knot condition removes last node
      IF(kbcn == 0) imax=n-1

!Set right BC values
      !Default for not-a-knot
      an=0.0_rspec
      bn=0.0_rspec

      IF(kbcn == 1) THEN

        !First derivative specified
        an=f(2,n)

      ELSEIF(kbcn == 2) THEN

        !Second derivative specified
        bn=f(3,n)

      ELSEIF(kbcn == 5) THEN

        !Match first derivative to last two points
        an=(f(1,n)-f(1,n-1))/(x(n)-x(n-1))

      ELSEIF(kbcn == 6) THEN

        !Match second derivative to first three points
        bn=2.0_rspec*((f(1,n)-f(1,n-1))/(x(n)-x(n-1))                          &
     &     -(f(1,n-1)-f(1,n-2))/(x(n-1)-x(n-2)))/(x(n)-x(n-2))

      ENDIF

!Clear derivatives, f(2:4,1:n), and work array
      f(2:4,1:n)=0.0_rspec
      wk(1:n)=0.0_rspec

!-------------------------------------------------------------------------------
!Evaluate coefficients
!-------------------------------------------------------------------------------
!Two and three nodes are special cases
      IF(n == 2) THEN

        !Coefficients for n=2
        f(2,1)=(f(1,2)-f(1,1))/(x(2)-x(1))
        f(3,1)=0.0_rspec
        f(4,1)=0.0_rspec
        f(2,2)=f(2,1)
        f(3,2)=0.0_rspec
        f(4,2)=0.0_rspec

!Set up tridiagonal system for A*y=B where y(i) are the second
!  derivatives at the knots
!  f(2,i) are the diagonal elements of A
!  f(4,i) are the off-diagonal elements of A
!  f(3,i) are the B elements/3, and will become c/3 upon solution

      ELSEIF(n > 2) THEN

        f(4,1)=x(2)-x(1)
        f(3,2)=(f(1,2)-f(1,1))/f(4,1)

        DO i=2,n-1 !Over nodes

          f(4,i)=x(i+1)-x(i)
          f(2,i)=2.0_rspec*(f(4,i-1)+f(4,i))
          f(3,i+1)=(f(1,i+1)-f(1,i))/f(4,i)
          f(3,i)=f(3,i+1)-f(3,i)

        ENDDO !Over nodes

!Apply left BC
        IF(kbc1 == -1) THEN

          !Periodic
          f(2,1)=2.0_rspec*(f(4,1)+f(4,n-1))
          f(3,1)=(f(1,2)-f(1,1))/f(4,1)-(f(1,n)-f(1,n-1))/f(4,n-1)
          wk(1)=f(4,n-1)
          wk(2:n-3)=0.0_rspec
          wk(n-2)=f(4,n-2)
          wk(n-1)=f(4,n-1)

        ELSEIF((kbc1 == 1) .OR. (kbc1 == 3) .OR. (kbc1 == 5)) THEN

          !First derivative conditions
          f(2,1)=2.0_rspec*f(4,1)
          f(3,1)=(f(1,2)-f(1,1))/f(4,1)-a1

        ELSEIF((kbc1 == 2) .OR. (kbc1 == 4) .OR. (kbc1 == 6)) THEN

          !Second derivative conditions
          f(2,1)=2.0_rspec*f(4,1)
          f(3,1)=f(4,1)*b1/3.0_rspec
          f(4,1)=0.0_rspec

        ELSEIF(kbc1 == 7) THEN

          !Third derivative condition
          f(2,1)=-f(4,1)
          f(3,1)=f(3,3)/(x(4)-x(2))-f(3,2)/(x(3)-x(1))
          f(3,1)=f(3,1)*f(4,1)**2/(x(4)-x(1))

        ELSE

          !Not-a-knot condition
          f(2,2)=f(4,1)+2.0_rspec*f(4,2)
          f(3,2)=f(3,2)*f(4,2)/(f(4,1)+f(4,2))

        ENDIF

!Apply right BC
        IF((kbcn == 1) .OR. (kbcn == 3) .OR. (kbcn == 5)) THEN

          !First derivative conditions
          f(2,n)=2.0_rspec*f(4,n-1)
          f(3,n)=-(f(1,n)-f(1,n-1))/f(4,n-1)+an

        ELSEIF((kbcn == 2) .OR. (kbcn == 4) .OR. (kbcn == 6)) THEN

          !Second derivative conditions
          f(2,n)=2.0_rspec*f(4,n-1)
          f(3,n)=f(4,n-1)*bn/3.0_rspec
          f(4,n-1)=0.0_rspec

        ELSEIF(kbcn == 7) THEN

          !Third derivative condition
          f(2,n)=-f(4,n-1)
          f(3,n)=f(3,n-1)/(x(n)-x(n-2))-f(3,n-2)/(x(n-1)-x(n-3))
          f(3,n)=-f(3,n)*f(4,n-1)**2/(x(n)-x(n-3))

        ELSEIF(kbc1 /= -1) THEN

          !Not-a-knot condition
          f(2,n-1)=2.0_rspec*f(4,n-2)+f(4,n-1)
          f(3,n-1)=f(3,n-1)*f(4,n-2)/(f(4,n-1)+f(4,n-2))

        ENDIF

!Limit solution for only three points in domain
        IF(n == 3) THEN

          f(3,1)=0.0_rspec
          f(3,n)=0.0_rspec

        ENDIF

!Solve system of equations for second derivatives at the knots
        IF(kbc1 == -1) THEN

          !Periodic BC - requires special treatment at ends
          !Forward elimination
          DO i=2,n-2 !Over nodes in forward elimination

            t=f(4,i-1)/f(2,i-1)
            f(2,i)=f(2,i)-t*f(4,i-1)
            f(3,i)=f(3,i)-t*f(3,i-1)
            wk(i)=wk(i)-t*wk(i-1)
            q=wk(n-1)/f(2,i-1)
            wk(n-1)=-q*f(4,i-1)
            f(2,n-1)=f(2,n-1)-q*wk(i-1)
            f(3,n-1)=f(3,n-1)-q*f(3,i-1)

          ENDDO !Over nodes in forward elimination

          !Correct the n-1 element
          wk(n-1)=wk(n-1)+f(4,n-2)

          !Complete the forward elimination
          !wk(n-1) and wk(n-2) are the off-diag elements of the lower corner
          t=wk(n-1)/f(2,n-2)
          f(2,n-1)=f(2,n-1)-t*wk(n-2)
          f(3,n-1)=f(3,n-1)-t*f(3,n-2)

          !Back substitution
          f(3,n-1)=f(3,n-1)/f(2,n-1)
          f(3,n-2)=(f(3,n-2)-wk(n-2)*f(3,n-1))/f(2,n-2)

          DO ib=3,n-1 !Over nodes in back substitution

            i=n-ib
            f(3,i)=(f(3,i)-f(4,i)*f(3,i+1)-wk(i)*f(3,n-1))/f(2,i)

          ENDDO !Over nodes in back substitution

          f(3,n)=f(3,1)

        ELSE

          !Non-periodic BC
          !Forward elimination
          !For not-a-knot BC the off-diagonal end elements are not equal
          DO i=imin+1,imax !Over nodes in forward elimination

            IF((i == n-1) .AND. (imax == n-1)) THEN

              t=(f(4,i-1)-f(4,i))/f(2,i-1)

            ELSE

              t=f(4,i-1)/f(2,i-1)

            ENDIF

            IF((i == imin+1) .AND. (imin == 2)) THEN

              f(2,i)=f(2,i)-t*(f(4,i-1)-f(4,i-2))

            ELSE

              f(2,i)=f(2,i)-t*f(4,i-1)

            ENDIF

            f(3,i)=f(3,i)-t*f(3,i-1)

          ENDDO !Over nodes in forward elimination

          !Back substitution
          f(3,imax)=f(3,imax)/f(2,imax)

          DO ib=1,imax-imin !Over nodes in back substitution

            i=imax-ib

            IF((i == 2) .AND. (imin == 2)) THEN

              f(3,i)=(f(3,i)-(f(4,i)-f(4,i-1))*f(3,i+1))/f(2,i)

            ELSE

              f(3,i)=(f(3,i)-f(4,i)*f(3,i+1))/f(2,i)

            ENDIF

          ENDDO !Over nodes in back substitution

          !Reset d array to step size
          f(4,1)=x(2)-x(1)
          f(4,n-1)=x(n)-x(n-1)

          !Set f(3,1) for not-a-knot
          IF((kbc1 <= 0) .OR. (kbc1 > 5)) THEN

            f(3,1)=(f(3,2)*(f(4,1)+f(4,2))-f(3,3)*f(4,1))/f(4,2)

          ENDIF

          !Set f(3,n) for not-a-knot
          IF((kbcn <= 0) .OR. (kbcn > 5)) THEN

            f(3,n)=f(3,n-1)+(f(3,n-1)-f(3,n-2))*f(4,n-1)/f(4,n-2)

          ENDIF

        ENDIF

!f(3,i) is now the sigma(i) of the text and f(4,i) is the step size
!Compute polynomial coefficients
        DO i=1,n-1 !Over nodes

          f(2,i)=(f(1,i+1)-f(1,i))/f(4,i)-f(4,i)*(f(3,i+1)                     &
     &           +2.0_rspec*f(3,i))
          f(4,i)=(f(3,i+1)-f(3,i))/f(4,i)
          f(3,i)=6.0_rspec*f(3,i)
          f(4,i)=6.0_rspec*f(4,i)

        ENDDO !Over nodes

        IF(kbc1 == -1) THEN

          !Periodic BC
          f(2,n)=f(2,1)
          f(3,n)=f(3,1)
          f(4,n)=f(4,1)

        ELSE

          !All other BCs
          hn=x(n)-x(n-1)
          f(2,n)=f(2,n-1)+hn*(f(3,n-1)+0.5_rspec*hn*f(4,n-1))
          f(3,n)=f(3,n-1)+hn*f(4,n-1)
          f(4,n)=f(4,n-1)

          IF((kbcn == 1) .OR. (kbcn == 3) .OR. (kbcn == 5)) THEN

            !First derivative BC
            f(2,n)=an

          ELSEIF((kbcn == 2) .OR. (kbcn == 4) .OR. (kbcn == 6)) THEN

            !Second derivative BC
            f(3,n)=bn

          ENDIF

        ENDIF

      ENDIF

      END SUBROUTINE SPLINE1_FIT

      SUBROUTINE SPLINE1_EVAL(k_vopt,n,u,x,f,i,value)
!-------------------------------------------------------------------------------
!SPLINE1_EVAL evaluates the cubic spline function and its derivatives
!  W.A.Houlberg, P.I. Strand, D.McCune 8/2001
!Input:
!  k_vopt(1)           -calculate the function [0=off, else=on]
!  k_vopt(2)           -calculate the first derivative [0=off, else=on]
!  k_vopt(3)           -calculate the second derivative [0=off, else=on]
!  n                   -number of data points [-]
!  u                   -abscissa at which the spline is to be evaluated [arb]
!  x                   -array containing the data abcissas [arb]
!  f(4,n)              -array containing the data ordinates [arb]
!Input/Output:
!  i                   -guess for target lower bound input if 1<i<n [-]
!  i                   -target lower bound on output [-]
!Output:
!  value( )            -ordering is a subset of the sequence under k_vopt [arb]
!  value(1)            =function or lowest order derivative requested [arb]
!  value(2)            =next order derivative requested [arb]
!  value(3)            =second derivative if all k_vopt are non-zero [arb]
!Comments:
!  s=f(1,i)+f(2,i)*dx+f(3,i)*dx**2/2!+f(4,i)*dx**3/3!
!  s'=f(2,i)+f(3,i)*dx+f(4,i)*dx**2/2!
!  s''=f(3,i)+f(4,i)*dx
!    where dx=u-x(i) and x(i).lt.u.lt.x(i+1)
!  If u <= x(1) then i=1 is used
!  If u >= x(n) then i=n is used
!-------------------------------------------------------------------------------

!Declaration of input variables
      INTEGER, INTENT(IN) ::                                                   &
     & k_vopt(:),                                                              &
     & n

      REAL(KIND=rspec), INTENT(IN) ::                                          &
     & u,                                                                      &
     & f(:,:),                                                                 &
     & x(:)

!Declaration of input/output variables
      INTEGER, INTENT(INOUT) ::                                                &
     & i

!Declaration of output variables
      REAL(KIND=rspec), INTENT(OUT) ::                                         &
     & value(:)

!Declaration of local variables
      INTEGER ::                                                               &
     & j

      REAL(KIND=rspec) ::                                                      &
     & dx

!-------------------------------------------------------------------------------
!Find target cell
!-------------------------------------------------------------------------------
      CALL SPLINE1_SEARCH(x,n,u,i)

!-------------------------------------------------------------------------------
!Set values outside of range to endpoints
!-------------------------------------------------------------------------------
      IF(i <= 0) THEN

        i=1
        dx=0.0_rspec

      ELSEIF(i >= n) THEN

        i=n
        dx=0.0_rspec

      ELSE

       dx=u-x(i)

      ENDIF

!-------------------------------------------------------------------------------
!Evaluate spline
!-------------------------------------------------------------------------------
      j=0

!Value
      IF(k_vopt(1) /= 0) THEN

        j=j+1
        value(j)=f(1,i)+dx*(f(2,i)+0.5_rspec*dx*(f(3,i)                        &
     &           +dx*f(4,i)/3.0_rspec))

      ENDIF

!First derivative
      IF(k_vopt(2) /= 0) THEN

        j=j+1
        value(j)=f(2,i)+dx*(f(3,i)+0.5_rspec*dx*f(4,i))

      ENDIF

!Second derivative
      IF(k_vopt(3) /= 0) THEN

        j=j+1
        value(j)=f(3,i)+dx*f(4,i)

      ENDIF

      END SUBROUTINE SPLINE1_EVAL

      SUBROUTINE SPLINE1_INTERP(k_vopt,n0,x0,f,n1,x1,value,iflag,              &
     &                          message,                                       &
     &                          K_BC1,K_BCN)
!-------------------------------------------------------------------------------
!SPLINE1_INTERP performs a spline interpolation from the x0 mesh to the x1 mesh
!References:
!  W.A.Houlberg, P.I.Strand, D.McCune 8/2001
!Input:
!  k_vopt(1)           -calculate the function [0=off, else=on]
!  k_vopt(2)           -calculate the first derivative [0=off, else=on]
!  k_vopt(3)           -calculate the second derivative [0=off, else=on]
!  n0                  -number of source abscissas [-]
!  x0(n0)              -source abscissas (in increasing order) [arb]
!  f(1,n0)             -source values [arb]
!  n1                  -number of target abscissas [-]
!  x1(n1)              -target abscissas (in increasing order) [arb]
!Input/Output:
!  f(2,1)              -input value of s'(x1) for K_BC1=1 [arb]
!  f(2,n0)             -input value of s'(xn) for K_BCN=1 [arb]
!  f(3,1)              -input value of s''(x1) for K_BC1=2 [arb]
!  f(3,n0)             -input value of s''(xn) for K_BCN=2 [arb]
!  f(4,n)              -arrays of n0 spline coefficients
!                       f(2,i)=s'(x0(i))/1!
!                       f(3,i)=s''(x0(i))/2!
!                       f(4,i)=s'''(x0(i))/3!
!Output:
!  value( ,j)          -ordering is a subset of the sequence under k_vopt [arb]
!  value(1,j)          =function or lowest order derivative requested [arb]
!  value(2,j)          =next order derivative requested [arb]
!  value(3,j)          =second derivative if all k_vopt are non-zero [arb]
!  iflag               -error and warning flag [-]
!                      =-1 warning
!                      =0 no warnings or errors
!                      =1 error
!  message             -warning or error message (character)
!Optional input:
!  K_BC1               -option for BC at x(1) [-]
!                      =-1 periodic, ignore K_BCN
!                      =0 not-a-knot (default)
!                      =1 s'(x1) = input value of f(2,1)
!                      =2 s''(x1) = input value of f(3,1)
!                      =3 s'(x1) = 0.0
!                      =4 s''(x1) = 0.0
!                      =5 match first derivative to first 2 points
!                      =6 match second derivative to first 3 points
!                      =7 match third derivative to first 4 points
!                      =else use not-a-knot
!  K_BCN               -option for boundary condition at x(n) [-]
!                      =0 not-a-knot (default)
!                      =1 s'(x1) = input value of f(2,1)
!                      =2 s''(x1) = input value of f(3,1)
!                      =3 s'(x1) = 0.0
!                      =4 s''(x1) = 0.0
!                      =5 match first derivative to first 2 points
!                      =6 match second derivative to first 3 points
!                      =7 match third derivative to first 4 points
!                      =else use knot-a-knot
!-------------------------------------------------------------------------------

!Declaration of input variables
      INTEGER, INTENT(IN) ::                                                   &
     & k_vopt(:),                                                              &
     & n0,n1

      REAL(KIND=rspec), INTENT(IN) ::                                          &
     & x0(:),x1(:)

!Declaration of input/output variables
      REAL(KIND=rspec), INTENT(INOUT) ::                                       &
     & f(:,:)

!Declaration of output variables
      INTEGER, INTENT(OUT) ::                                                  &
     & iflag
     
      CHARACTER(len=*), INTENT(OUT) ::                                         &
     & message

      REAL(KIND=rspec), INTENT(OUT) ::                                         &
     & value(:,:)

!Declaration of optional input variables
      INTEGER, INTENT(IN), OPTIONAL ::                                         &
     & K_BC1,K_BCN

!Declaration of local variables
      INTEGER ::                                                               &
     & i,j

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!If no target values are requested, return
      IF(n1 <= 0) THEN

        iflag=-1
        message='SPLINE1_INTERP/WARNING(1):no points in output array'
        GOTO 9999

      ENDIF

!Make sure there are at least two nodes in the input grid
      IF(n0 < 2) THEN

        iflag=1
        message='SPLINE1_INTERP/ERROR(2):<2 points in source array'
        GOTO 9999

      ENDIF

!-------------------------------------------------------------------------------
!Get spline coefficients
!-------------------------------------------------------------------------------
      IF(PRESENT(K_BC1) .AND.                                                  &
     &   PRESENT(K_BCN)) THEN

        CALL SPLINE1_FIT(n0,x0,f,                                              &
     &                   K_BC1=K_BC1,                                          &
     &                   K_BCN=K_BCN)

      ELSEIF(PRESENT(K_BC1)) THEN

        CALL SPLINE1_FIT(n0,x0,f,                                              &
     &                   K_BC1=K_BC1)

      ELSEIF(PRESENT(K_BCN)) THEN

        CALL SPLINE1_FIT(n0,x0,f,                                              &
     &                   K_BCN=K_BCN)

      ELSE

        CALL SPLINE1_FIT(n0,x0,f)

      ENDIF

!-------------------------------------------------------------------------------
!Interpolate onto x1 mesh
!-------------------------------------------------------------------------------
      i=1

      DO j=1,n1 !Over nodes

        CALL SPLINE1_EVAL(k_vopt,n0,x1(j),x0,f,i,value(1:3,j))

      ENDDO !Over nodes

 9999 CONTINUE

      END SUBROUTINE SPLINE1_INTERP

      SUBROUTINE SPLINE1_INTEG(k_order,n0,x0,f,n1,x1,value,iflag,              &
     &                         message)
!-------------------------------------------------------------------------------
!SPLINE1_INTEG evaluates the integral f(x)*x**k_order, where f(x) is a cubic
!  spline function and k_order is 0 or 1
!References:
!  Forsythe, Malcolm, Moler, Computer Methods for Mathematical
!    Computations, Prentice-Hall (1977) 76
!  W.A.Houlberg, P.I.Strand, D.McCune 8/2001
!Input:
!  k_order             -exponent in integral (s(x)*x**k_order) [-]
!                      =0 or 1 only
!  n0                  -number of source nodes [-]
!  x0(n0)              -source abcissas [arb]
!  f(1,n0)             -source ordinates [arb]
!  f(2-4,n0)           -spline coefficients computed by SPLINE1_FIT [arb]
!  n1                  -number of output nodes (may be 1) [-]
!  x1(n1)              -output abcissas at which the integral is wanted [arb]
!Output:
!  value(n1)           -integral of f(x)*x**k_order from x0(1) to x1(i) [arb]
!  iflag               -error and warning flag [-]
!                      =-1 warning
!                      =0 no warnings or errors
!                      =1 error
!  message             -warning or error message (character)
!-------------------------------------------------------------------------------

!Declaration of input variables
      INTEGER, INTENT(IN) ::                                                   &
     & k_order,                                                                &
     & n0,n1

      REAL(KIND=rspec), INTENT(IN) ::                                          &
     & x0(:),x1(:),                                                            &
     & f(:,:)

!Declaration of output variables
      INTEGER, INTENT(OUT) ::                                                  &
     & iflag
     
      CHARACTER(len=*), INTENT(OUT) ::                                         &
     & message

      REAL(KIND=rspec), INTENT(OUT) ::                                         &
     & value(:)

!Declaration of local variables
      INTEGER ::                                                               &
     & i,j,                                                                    &
     & ido,jdo

      REAL(KIND=rspec) ::                                                      &
     & add,                                                                    &
     & dx,                                                                     &
     & sum,                                                                    &
     & xnew,xold

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!If no target values are requested, return
      IF(n1 <= 0) GOTO 9999

      value(1:n1)=0.0_rspec

!Integral value
      sum=0.0_rspec

!Source index and abscissa
      j=1
      xold=x0(j)

!-------------------------------------------------------------------------------
!Integrate over target abscissas
!-------------------------------------------------------------------------------
      DO i=1,n1 !Over target abscissas

!Find dx to next source or target abscissa
        ido=0

        DO WHILE(ido == 0) !Over source abscissas

          IF(x1(i) < x0(j+1)) THEN

            !Hit target abscissa
            xnew=x1(i)
            ido=1
            jdo=0

          ELSEIF(x1(i) == x0(j+1)) THEN

            !Hit both source and target abscissas
            xnew=x1(i)
            ido=1
            jdo=1

          ELSEIF(x1(i) > x0(j+1)) THEN

            !Hit source abscissa
            xnew=x0(j+1)
            ido=0
            jdo=1

          ENDIF

!Integrate over dx
          dx=xnew-xold

          IF(k_order == 1) THEN

            !Integral x s(x)dx
            add=dx*(xold*f(1,j)+dx*((xold*f(2,j)+f(1,j))/2.0_rspec             &
     &          +dx*((xold*f(3,j)/2.0_rspec+f(2,j))/3.0_rspec                  &
     &          +dx*((xold*f(4,j)/3.0_rspec+f(3,j))/8.0_rspec                  &
     &          +dx*f(4,j)/30.0_rspec))))

          ELSE

            !Integral s(x)dx
            add=dx*(f(1,j)+dx*(f(2,j)+dx*(f(3,j)                               &
     &          +dx*f(4,j)/4.0_rspec)/3.0_rspec)/2.0_rspec)

          ENDIF

!Add increment and update endpoint
          sum=sum+add
          xold=xnew

!Check whether to increment source index
          IF(jdo == 1) j=j+1

!Check whether target index is in range
          IF(j > n0) THEN

            iflag=1
            message='LINEAR1_INTEG/ERROR:target node out of range'
            GOTO 9999

          ENDIF

!Set integral value
          value(i)=sum

        ENDDO !Over source abscissas

      ENDDO !Over target abscissas
   
 9999 CONTINUE
      
      END SUBROUTINE SPLINE1_INTEG

      SUBROUTINE SPLINE1_SEARCH(x,n,xl,jlo)
!-------------------------------------------------------------------------------
!SPLINE1_SEARCH is a correlated table search routine to find the indices of the
!  array x that bound xl
!References:
!  W.A.Houlberg, P.I.Strand, D.McCune 8/2001
!Input:
!  x(n)                -monotonically increasing array of abscissas [arb]
!  n                   -number of abscissas [-]
!  xl                  -target value [arb]
!Input/Output:
!  jlo                 -input starting lower index [-]
!                      <1     binary search
!                      =1,n-1 use value
!                      >n-1   binary search
!  jlo                 -output starting lower index [-]
!                      =0     xl < x(1) 
!                      =1     x(1) <= xl <= x(2)
!                      =2,n-1 x(jlo) < xl <= x(jlo+1)
!                      =n     x(jlo) > x(n)
!Comments:
!  This is similar to the Numerical Recipes routine HUNT
!-------------------------------------------------------------------------------

!Declaration of input variables
      INTEGER, INTENT(IN) ::                                                   &
     & n

      REAL(KIND=rspec), INTENT(IN) ::                                          &
     & xl,                                                                     &
     & x(:)

!Declaration of input/output variables
      INTEGER, INTENT(INOUT) ::                                                &
     & jlo

!Declaration of local variables
      INTEGER ::                                                               &
     & inc,                                                                    &
     & jhi,                                                                    &
     & jmid

!-------------------------------------------------------------------------------
!Check lower end of array, first two points
!-------------------------------------------------------------------------------
      IF(xl < x(1)) THEN

        !Target is below node 1
        jlo=0

      ELSEIF(xl <= x(2)) THEN

        !Target is between nodes 1 and 2 inclusive
        jlo=1
                     
!-------------------------------------------------------------------------------
!Check middle range
!-------------------------------------------------------------------------------
      ELSEIF(xl <= x(n)) THEN

        !Target is between nodes 2 and n
        IF(jlo < 1 .OR. jlo > (n-1)) THEN

          !jlo from previous call is unusable
          jlo=2
          jhi=n

        !Bracket target value
        ELSE

          !Start with jlo from previous call
          inc=1

          IF(xl > x(jlo)) THEN

            !Search up
            jhi=jlo+1

            DO WHILE(xl > x(jhi))

              inc=2*inc
              jlo=jhi
              jhi=MIN(jlo+inc,n)

            ENDDO

          ELSE

            !Search down
            jhi=jlo
            jlo=jlo-1

            DO WHILE(xl <= x(jlo))

              inc=inc+inc
              jhi=jlo
              jlo=MAX(jlo-inc,1)

            ENDDO

          ENDIF

        ENDIF

!Bisection
        DO WHILE(jhi-jlo > 1)

          jmid=(jhi+jlo)/2

          IF(xl > x(jmid)) THEN

            jlo=jmid

          ELSE

            jhi=jmid

          ENDIF

        ENDDO

!-------------------------------------------------------------------------------
!Target is above node n
!-------------------------------------------------------------------------------
      ELSE

        jlo=n

      ENDIF

      END SUBROUTINE SPLINE1_SEARCH
      
      END MODULE SPLINE1_MOD

