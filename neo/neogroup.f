
PROGRAM neo
! **********************************************************************
! Program: NEO
! Authors: Winfried Kernbichler and Sergei Kasilov
! E-mail:  kernbichler@itp.tu-graz.ac.at
! Version: 3.02 
! Date:    February 2, 2001
! 
! ===================================================================
! Changes from 3.0
!  All routines have explicit interfaces in modules
!  lab_swi tested for ORNL (lab_swi 2 does not read sqrtg00
!   (has no influence at the moment - sqrtg00 not used)
!  Compiles without errors on 
!    IBM rs_aix42 Gaching
!    DEC ALPHA
!    Linux Lahey  f95 compiler
!    Linux Absoft f90 compiler
!
!  Runs without run time errors on all 4 machines
!
!  Results differ on the IBM because integraton along field lines
!    stops earlier - fixed in 3.02 (see below)
! ===================================================================
! Changes in 3.02
!
!    initialization of nstep_max_c = nstep_max added in flint_bo.f
!
!  Reading of ORNL-Data (Don Spong)
! ===================================================================
!
! This code is freely distributed to colleagues working in the field of
! stellarators. Users of NEO should acknowledge the authors and any 
! modifications made when publishing or presenting results. Problems should 
! be reported to the authors first. Any bugfixes or improvements carried out 
! by users should be made available to authors, so that these improvements 
! can be incorporated into the distribution.
! 
! The main purpose of the code is to calculate the effective helical 
! ripple $\eps_{eff}^{3/2}$ as defined in Phys. Plasmas 6(12) 4622-4632.
!
! The part of the code which computes parallel current densities (EPS 2000
! Budapest, APS 2000 Quebec City, VIII Ukrainean Conference on Plasma Physics
! is highly "experimental" and in development. Please use it only with care!
! 
! This version works in Boozer coordinates and was extensively tested 
! with VMEC output for PPPL NCSX equilibria. For this purpose the VMEC
! output was converted with a mapping program to Boozer coordinates as
! used at PPPL (Bmns-Files). Other input formats will be included in 
! future versions.
!
! There are three different input files distributed:
!  neo.in.acc   - high accuracy (slow)
!  neo.in.ref   - good accuracy (medium speed)
!  neo.in.opt   - input variables suitable for PPPL optimizer
!                 fast (380 sec on PPPL Alpha SATURN for 47 flux surfaces)
!
! There are some important input quantities responsible for accuracy versus
! speed of computation:
!
! Recommendations to date are based on the quasi-axisymmetric PPPL cases
!
! theta_n, phi_n: Grid size in poloidal and toroidal direction
!                 Determines accuracy of B representation for the double-
!                  periodic spline functions
!                 Reference: 200x200
!                 It does not make sense to go below 100x100
!
! npart:          Number of test particles for $J_{\perp}$ integration
!                 Determines the accuracy of $J_{\perp}$ integration (sum)
!                 Important to catch the "fine-structure" of B along the mfl
!                 Reference: 100
!                 To go below 50 seems to be dangerous
!
! acc_req:        Required accuracy for each individual integration along
!                  a magnetic field line (good guess!)
!                 Reference: 0.01
!                 Good choice, one should no go much lower
!
! nstep_per:      Number of integration steps per field period
!                 Determines the accuracy of the Runge-Kutta (fixed step size)
!                  solver for the integration along a mfl.
!                 Reference: 50
!                 Good choice, dangerous to go to lower values
!
! nstep_min:      Minimum number of field periods after which a decision is
!                  taken about continuing until the required accuracy is
!                  reached (less steps than nstep_max) or switching to the
!                  computational mode for near-rational surfaces
!                 Reference: 2000
!                 Good choice: 500, lower values do not make much sense
!
! nstep_max:      Maximum number of field periods to be followed
!                 See nstep_min
!                 Reference: 10000
!                 Good choice 2000, lower values do not make much sense
!
! nbins:          Number of bins in poloidal direction to be filled on a
!                  toroidal cut when hit by a mfl. This bins should be 
!                  filled evenly to ensure good coverage of the surface
!                  by following the mfl.
!                 Determines the number of starting points in theta for the
!                  field line integration if the flux surface is close to
!                  a rational one.
!                 Reference: 200
!                 Good choice: 100 (keeps time spend on rational surfaces
!                  lower with still acceptable accuracy)
!
! If one is interested in additional quantities, like
!
!   epspar(i)-  partial contributions to effective ripple
!               from different classes of trapped particles
!               index i=1,...,multra  means single-trapped,
!               double-trapped,...,up to multra-1,
!               epspar(multra) contains the contributions
!               from all upper classes 
!   ctrone   -  fraction of particles trapped once (see Eq.(36))
!   ctrtot   -  fraction of all trapped particles  (see Eq.(36))
!   bareph   -  $\bar \epsilon_h$ - 'ripple' amplitude defined
!               through fraction of single-trapped particles
!   barept   -  $\bar \epsilon_t$ - 'toroidal' amplitude
!               defined through full fraction of trapped prts
!                                 
! one has to set
!   multra to values greater than one, and
!   eout_swi to 2 (detailed output)
! To get a good resolution in this quantities mainly the number of particles
! has to be set to higher values (npart = 1000).
!
! There are some PPPL related things in neo_sub. If there is need for 
! things in other laboratories one could use lab_swi to add something there.
! If one needs other input programs, one should use inp_swi to make the choice.
!
! We agreed with PPPL to use the B(0,0)- and the R(0,0)-components of the 
! innermost flux surface as reference values for the $\eps_{eff}^{3/2}$
! calculation (ref_swi=1). Other choices (like using $B_{max}$ on each flux 
! surface) should be done in neo.f90 using ref_swi.
!
! **********************************************************************
! Modules
! **********************************************************************
  USE neo_precision
  USE neo_units
  USE neo_parameters
  USE neo_control
  USE neo_input
  USE neo_work
  USE neo_exchange
  USE neo_output
  USE sizey_bo
  USE safe_open_mod                          ! SPH
! **********************************************************************
! Local Definitions
! **********************************************************************
  IMPLICIT NONE

  INTEGER       ::                     npsi, istat
  INTEGER       ::                     fluxs_arr_i
  REAL(kind=dp) ::                     reff
  REAL(kind=dp) ::                     psi,dpsi
  REAL(kind=dp) ::                     b_ref,r_ref
! **********************************************************************
! Read input from control file
! **********************************************************************
  CALL neo_read_control
  IF (write_progress .NE. 0) WRITE (w_us,*) 'after  neo_read_control'
! **********************************************************************
! Initialize general data 
! **********************************************************************
  IF (write_progress .NE. 0) WRITE (w_us,*) 'before neo_init'
  CALL neo_init(npsi) 
  IF (write_progress .NE. 0) WRITE (w_us,*) 'after  neo_init'
!
  call safe_open(w_u3, istat, out_file, 'replace', 'formatted')
  if (istat .ne. 0) stop 'Error opening NEO output file'
  IF (calc_cur .EQ. 1) THEN
     OPEN(unit=w_u9,file=cur_file)
  END IF
! *******************************************************************
! Allocate and fill fluxs_arr if not done in neo_control
! ******************************************************************* 
  IF (no_fluxs .LE. 0) THEN
     no_fluxs = npsi
     ALLOCATE ( fluxs_arr(no_fluxs) )
     DO fluxs_arr_i = 1, no_fluxs
        fluxs_arr(fluxs_arr_i) = fluxs_arr_i
     END DO
  ENDIF
! *******************************************************************
! Loop For Magnetic Surfaces
! ******************************************************************* 
  reff=0
  DO fluxs_arr_i = 1, no_fluxs
!    psi_ind = fluxs_arr(fluxs_arr_i)
     psi_ind = fluxs_arr_i                                ! LPK

     IF (psi_ind .GE. 1 .AND. psi_ind .LE. npsi) THEN
! **********************************************************************
! Initialize data for a specific flux surface index psi_ind
! **********************************************************************
        IF (write_progress .NE. 0) WRITE (w_us,*)                         & 
             'before neo_init_s, psi_ind: ',psi_ind
        CALL neo_init_s(psi,dpsi)
        IF(psi_ind.EQ.1) dpsi=psi
        IF (write_progress .NE. 0) WRITE (w_us,*) 'after  neo_init_s'
! *******************************************************************
! Do eps_eff Calculation Using the Spline Arrays
! ******************************************************************* 
        IF (write_progress .NE. 0) WRITE (w_us,*) 'before flint_bo'
        CALL flint_bo()
        IF (write_progress .NE. 0) WRITE (w_us,*) 'after  flint_bo',ierr
        reff=reff+drdpsi*dpsi
! *******************************************************************
! Do parallel current Calculations
! ******************************************************************* 
        IF (calc_cur .EQ. 1) THEN
           IF (write_progress .NE. 0) WRITE (w_us,*) 'before flint_cur'
           CALL flint_cur()
           IF (write_progress .NE. 0) WRITE (w_us,*) 'after  flint_cur',ierr
        END IF
! *******************************************************************
! Rescale $e_{eff}^{3/2}$
! *******************************************************************
        IF (ref_swi .EQ. 1) THEN
           b_ref = bmref_g
           r_ref = rt0_g
        ELSEIF (ref_swi .EQ. 2) THEN
           b_ref = bmref
           r_ref = rt0
        ELSE
           WRITE (w_us,*) 'FATAL: This ref_swi ',ref_swi,' is not implemented!'
           STOP
        END IF
        epstot = epstot * (b_ref/bmref)**2 * (r_ref/rt0)**2
        epspar = epspar * (b_ref/bmref)**2 * (r_ref/rt0)**2
! *******************************************************************
! Write Output
! *******************************************************************
        IF (eout_swi .EQ. 1) THEN
           WRITE(w_u3,'(1(1x,i8),5(1x,e17.10))')                    &
                fluxs_arr(fluxs_arr_i),                             &   
                epstot,reff,iota(psi_ind),b_ref,r_ref
        ELSEIF (eout_swi .EQ. 2) THEN
           WRITE(w_u3,'(1(1x,i8),12(1x,e17.10))')                   &
                fluxs_arr(fluxs_arr_i),                             &   
                epstot,reff,iota(psi_ind),b_ref,r_ref,              &
                epspar(1),epspar(2),ctrone,ctrtot,bareph,barept,yps
        ELSEIF (eout_swi .EQ. 10) THEN                              !LPK
           WRITE(w_u3,*) b_ref, r_ref, epstot                       !LPK
        ELSE
           WRITE(w_us,*) 'FATAL: This eout_swi ',eout_swi,' is not implemented!'
           STOP
        END IF
        IF (calc_cur .EQ. 1) THEN
           WRITE(w_u9,'(1(1x,i8),5(1x,e17.10))')                       &
                psi_ind,                                               &   
                lambda_b,                                              &
                lambda_ps1,lambda_ps2,                                 &
                lambda_b1,lambda_b2 
        END IF
     ELSE
        WRITE (w_us,*) 'Flux surface ',psi_ind,' does not exist!'
     END IF
  END DO
  CLOSE(w_u3)
  IF (calc_cur .EQ. 1) THEN
     CLOSE(w_u9)
  END IF
! *******************************************************************
! DeAllocate Storage Arrays
! *******************************************************************
  IF (write_progress .NE. 0) WRITE (w_us,*) 'before neo_dealloc'
  CALL neo_dealloc
  IF (write_progress .NE. 0) WRITE (w_us,*) 'after neo_dealloc'
! *******************************************************************
  STOP
END PROGRAM neo


SUBROUTINE flint_bo()
! eps_eff Integration
! and additional output specified in the paper
!
! Input variables:  bmref    -  $B_0$ - reference value of mod-B
!                   rt0      -  $R$   - reference value of big radius 
!                   phi0     -  starting value of $\phi$
!                   theta0   -  starting value of $\theta$
!
! Output variables: epspar(i)-  partial contributions to effective ripple
!                               from different classes of trapped particles
!                               index i=1,...,multra  means single-trapped,
!                               double-trapped,...,up to multra-1,
!                               epspar(multra) contains the contributions
!                               from all upper classes 
!                   epstot   -  total effective ripple - $\epsilon_{eff}^{3/2}$,
!                               the main result of computation
!                   ctrone   -  fraction of particles trapped once (see Eq.(36))
!                   ctrtot   -  fraction of all trapped particles  (see Eq.(36))
!                   bareph   -  $\bar \epsilon_h$ - 'ripple' amplitude defined
!                               through fraction of single-trapped particles
!                   barept   -  $\bar \epsilon_t$ - 'toroidal' amplitude
!                               defined through full fraction of trapped prts
!                                 
! **********************************************************************
! Modules
! **********************************************************************
  USE neo_precision
  USE neo_units
  USE neo_exchange
  USE neo_output
  USE partpa_bo
  USE neo_rk4dbo
  USE safe_open_mod                                   ! SPH
! **********************************************************************
! Local definitions
! **********************************************************************
  IMPLICIT NONE
!
  REAL(kind=dp) ::   etamin,etamax,hphi
  REAL(kind=dp) ::   phi
  REAL(kind=dp) ::   heta
  REAL(kind=dp) ::   coeps
  REAL(kind=dp) ::   adimax,aditot,adimax_s,aditot_s
  REAL(kind=dp) ::   y2,y3
  REAL(kind=dp) ::   y2_s,y3_s,y4,y4_s,y3npart,y3npart_s
  REAL(kind=dp) ::   bmod,gval,geodcu,qval
  REAL(kind=dp) ::   acc_d, acc_m
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE, TARGET ::  y
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE ::  bigint, bigint_s
  INTEGER                                  ::  i,m_cl,j1,n, istat
  INTEGER,       DIMENSION(:), ALLOCATABLE ::  iswst

  REAL(kind=dp)                            ::  epstot_check
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE ::  epspar_check

! INTEGER,       PARAMETER              ::  ntheta_bins = 100
! INTEGER                               ::  theta_bind
! INTEGER,       DIMENSION(ntheta_bins) ::  theta_bins
! INTEGER,       DIMENSION(ntheta_bins) ::  theta_bins_c
! REAL(kind=dp)                         ::  delta_theta
  REAL(kind=dp)                         ::  theta
!
  REAL(kind=dp)                         ::  theta0, phi0
  REAL(kind=dp)                         ::  theta_rs, theta_iota
  REAL(kind=dp)                         ::  theta_d,theta_d_min
  REAL(kind=dp)                         ::  theta_gap
  REAL(kind=dp)                         ::  iota_bar_fp
  
  INTEGER                               ::  istepc
  INTEGER                               ::  m, m_iota, n_iota
  INTEGER                               ::  nfl
  INTEGER                               ::  n_gap
  INTEGER                               ::  nstep_max_c
  INTEGER                               ::  exist_first_ratfl
  INTEGER                               ::  max_class
  REAL(kind=dp)                         ::  j_iota_i
  REAL(kind=dp)                         ::  add_on

  REAL(kind=dp),               POINTER  ::  p_theta, p_bm2
  REAL(kind=dp),               POINTER  ::  p_bm2gv, p_bm3ge
  REAL(kind=dp), DIMENSION(:), POINTER  ::  p_i, p_h
! **********************************************************************
! Allocation
! **********************************************************************
  ndim = npq + 2 * npart
  ALLOCATE( isw(npart) )
  ALLOCATE( ipa(npart) )
  ALLOCATE( icount(npart) )
  ALLOCATE( iswst(npart) )
  ALLOCATE( eta(npart) )

  ALLOCATE( y(ndim) )

  IF (.NOT. ALLOCATED(epspar) ) ALLOCATE( epspar(multra) )
  ALLOCATE( bigint(multra) )
  ALLOCATE( bigint_s(multra) )
  ALLOCATE( epspar_check(multra) )
! **********************************************************************
! Pointer
! **********************************************************************
  p_theta => y(1)
  p_bm2   => y(2)
  p_bm2gv => y(3)
  p_bm3ge => y(4)
  p_i     => y(npq+1:npq+npart)
  p_h     => y(npq+npart+1:npq+2*npart)
! **********************************************************************
! Initial settings
! **********************************************************************
  ierr   = 0
  hphi   = twopi/(nstep_per*nper)
  bmod0  = bmref
!
  exist_first_ratfl = 0
  hit_rat = 0
!
  theta0 = theta_bmax
  phi0   = phi_bmax
! **********************************************************************
! Definition of eta-values for particles

! **********************************************************************
  etamin = b_min / bmref
  etamax = b_max / bmref
  heta=(etamax-etamin)/(npart-1)
!  heta=(etamax-etamin)/(npart+1)
!  etamin = etamin + heta
!  etamax = etamax - heta
  etamin = etamin + heta / 2
  DO i=1,npart
     eta(i)=etamin+heta*(i-1)
  ENDDO
  coeps=pi*rt0**2*heta/(8*sqrt(2._dp))
  j_iota_i = curr_pol(psi_ind) + iota(psi_ind) * curr_tor(psi_ind)
! **********************************************************************
  isw    = 0
  ipa    = 0
  iswst  = 0
  ipmax  = 0
  icount = 0
  bigint = 0
  adimax = 0
  aditot = 0
! **********************************************************************
  nfp_rat = 0
  nfl_rat = 0
! **********************************************************************
  phi     = phi0
  y       = 0
  p_theta = theta0
! **********************************************************************
! Bins for theta at periodicity boundary
! **********************************************************************
! delta_theta = twopi / no_bins
! theta_bins  = 0
! theta_bind  = MOD(FLOOR(theta0/delta_theta),no_bins)+1
! theta_bins(theta_bind) = 1
!
  theta_d_min = twopi
! **********************************************************************
  IF (write_integrate == 1) THEN
     OPEN(unit=w_u5,file='conver.dat',status='replace',form='formatted')
  ENDIF
  IF (write_diagnostic == 1) THEN
     OPEN(unit=w_u7,file='diagnostic.dat',status='replace',form='formatted')
  ENDIF
! **********************************************************************
! Evaluation of Splines
! **********************************************************************
  CALL neo_eval(theta0,phi0,bmod,gval,geodcu,pard0,qval)
! **********************************************************************
! Intergration steps (Summation)
! **********************************************************************
  nstep_max_c = nstep_max
  istepc = 0
  max_class = 0
  DO n=1,nstep_max
     DO j1=1,nstep_per
        CALL RK4D_BO(phi,y,hphi)
        DO i=1,npart
           IF(isw(i).EQ.2) THEN
              m_cl = MIN(multra,ipa(i))
              m_cl = MAX(1,m_cl)
              IF(ipa(i).EQ.1) adimax = p_i(i)
              add_on = p_h(i)**2 / p_i(i) * iswst(i)
              bigint(m_cl) = bigint(m_cl) + add_on

              IF ( write_diagnostic .EQ. 1 .AND. iswst(i) .EQ. 1) THEN
                 istepc = istepc + 1
                 WRITE (w_u7,'(1x,i8,1x,i8,1x,i8,1x,e20.10)')          &
                      i,icount(i),ipa(i),add_on
                 IF ( ipa(i) .GT. max_class ) max_class = ipa(i)
              END IF

              iswst(i) = 1
              p_h(i) = 0
              p_i(i) = 0
              isw(i) = 0
              icount(i) = 0
              ipa(i) = 0
           END IF
        END DO
!
        IF(ipmax.EQ.1) THEN
           ipmax=0
           aditot=aditot+adimax
        END IF
     END DO
! **********************************************************************
! Checks after each field-period
! **********************************************************************
     IF (write_integrate == 1) THEN
        epstot_check = 0
        DO m_cl=1,multra
           epspar_check(m_cl) = coeps*bigint(m_cl)*p_bm2/p_bm2gv**2
           epstot_check = epstot_check + epspar_check(m_cl) 
        ENDDO
        WRITE(w_u5,format220)                                          &
             DBLE(n), epstot_check, p_bm3ge,                           &
             p_i(npart)/p_bm2, aditot/p_bm2
     END IF
! **********************************************************************
! Bins for theta at periodicity boundary
! **********************************************************************
     theta = p_theta
!    theta_bind = MOD(FLOOR(theta/delta_theta),ntheta_bins)+1
!    theta_bins(theta_bind) = theta_bins(theta_bind) + 1
!    acc_d = (MAXVAL(theta_bins)-MINVAL(theta_bins)) *                 &
!         SQRT(DBLE(ntheta_bins)) / n / 2
!    theta_bins_c = 0
!    WHERE (theta_bins .NE. 0)
!       theta_bins_c = 1
!    END WHERE
!    acc_m = (DBLE(SUM(theta_bins))/DBLE(SUM(theta_bins_c)) -          &
!             MINVAL(theta_bins)) *                                    &
!             SQRT(DBLE(ntheta_bins)) / n
! **********************************************************************
! Check for (near) rational surfaces
! **********************************************************************
     IF (n .LE. nstep_min) THEN
        theta_rs = theta - theta0
        IF (n .EQ. 1) THEN
           theta_iota = theta_rs
           iota_bar_fp = theta_iota / twopi
        END IF
        m = FLOOR(theta_rs / twopi)
        theta_rs = theta_rs - m * twopi
        IF (theta_rs .LE. pi) THEN
           theta_d = theta_rs
        ELSE
           theta_d = theta_rs - twopi
        END IF
        IF (ABS(theta_d) .LT. ABS(theta_d_min)) THEN
           theta_d_min = theta_d
           n_iota = n
           IF (theta_d .GE. ZERO) THEN
              m_iota = m
           ELSE
              m_iota = m+1
           END IF
        END IF
     END IF
     IF (n .EQ. nstep_min) THEN
        theta_gap = twopi / n_iota
        n_gap = n_iota * INT(ABS(theta_gap / theta_d_min))
        IF (n_gap .GT. nstep_min) THEN
           nstep_max_c = n_gap
        ELSE
           nstep_max_c = n_gap * CEILING( real(nstep_min,kind=dp) / n_gap )
        END IF

!!$        PRINT *, 'theta_iota:      ',theta_iota
!!$        PRINT *, 'theta_rs:        ',theta_rs
!!$        PRINT *, 'theta_d:         ',theta_d
!!$        PRINT *, 'theta_d_min:     ',theta_d_min
!!$        PRINT *, 'theta_gap:       ',theta_gap
!!$        PRINT *, 'n_iota:          ',n_iota
!!$        PRINT *, 'n_gap:           ',n_gap
!!$        PRINT *, 'nstep_min:       ',nstep_min
!!$        PRINT *, 'n_gap:           ',n_gap 
!!$        PRINT *, 'DBLE(nstep_min): ',DBLE(nstep_min) 
!!$        PRINT *, 'DBLE(n_gap):     ',DBLE(n_gap) 
!!$        PRINT *, 'DBLE/DBLE:       ',DBLE(nstep_min)/DBLE(n_gap)
!!$        PRINT *, 'CEILING:         ',CEILING(DBLE(nstep_min) / DBLE(n_gap)) 
!!$        PRINT *, 'nstep_max_c:     ',nstep_max_c 

! **********************************************************************
! Decision about Rational Surfaces
! **********************************************************************
        IF (nstep_max_c .GT. nstep_max) THEN
           hit_rat = 1
           nfp_rat = CEILING(ONE / acc_req / iota_bar_fp)
           IF (MODULO(nfp_rat,n_iota) .NE. 0) THEN
              nfp_rat = nfp_rat + n_iota - MODULO(nfp_rat,n_iota)
           END IF
           IF (nfp_rat .GE. nstep_min) THEN
              exist_first_ratfl = 1
              nstep_max_c = nfp_rat
           END IF
           nfl_rat = CEILING(real(no_bins,kind=dp) / n_iota)
           delta_theta_rat =  theta_gap / (nfl_rat + 1)
           IF ( calc_nstep_max .EQ. 1 ) hit_rat = 0
           IF ( hit_rat .EQ. 1 .AND. exist_first_ratfl .EQ. 0 ) EXIT
        END IF
     END IF
! **********************************************************************
! Accuracy
! **********************************************************************
!     IF (acc_d .LT. acc_req .AND. n .GT. nstep_min) THEN
!        EXIT
!     END IF
     IF ( calc_nstep_max .EQ. 0 .AND. n .EQ. nstep_max_c ) EXIT
  END DO
  nintfp = n        
!  
  y2      = p_bm2
  y3      = p_bm2gv
  y4      = p_bm3ge
  y3npart = p_i(npart)
! **********************************************************************
  IF (write_integrate == 1) THEN
     CLOSE(unit=w_u5)
  END IF
! **********************************************************************
! Calculations for rational surfaces
! **********************************************************************
  IF (hit_rat .EQ. 1) THEN

     IF (exist_first_ratfl .EQ. 0) THEN
        IF (write_integrate == 1) THEN
           OPEN(unit=w_u5,file='conver.dat',status='replace',form='formatted')
        ENDIF
        bigint  = 0
        adimax  = 0
        aditot  = 0
        y2      = 0
        y3      = 0
        y4      = 0
        y3npart = 0
     ELSE
        IF (write_integrate == 1) THEN
           OPEN(unit=w_u5,file='conver.dat',status='old',                   &
                position='append',form='formatted')
        END IF
     END IF
     DO nfl = exist_first_ratfl,nfl_rat
        bigint_s  = 0
        adimax_s  = 0
        aditot_s  = 0
        isw       = 0
        ipa       = 0
        iswst     = 0        
        ipmax     = 0
        icount    = 0
        phi       = phi0
        y         = 0
        theta     = theta0 + nfl * delta_theta_rat
        p_theta   = theta
        DO n=1,nfp_rat
           CALL neo_eval(theta,phi,bmod,gval,geodcu,pard0,qval)
           DO j1=1,nstep_per
              CALL RK4D_BO(phi,y,hphi)
              DO i=1,npart
                 IF(isw(i).EQ.2) THEN
                    isw(i) = 0
                    icount(i) = 0
                    m_cl = MIN(multra,ipa(i))
                    m_cl = MAX(1,m_cl)
                    IF(ipa(i).EQ.1) adimax_s = p_i(i)
                    ipa(i) = 0
                    bigint_s(m_cl) = bigint_s(m_cl)+p_h(i)**2/p_i(i)*iswst(i)
                    iswst(i) = 1
                    p_h(i) = 0
                    p_i(i) = 0
                 END IF
              END DO
              IF(ipmax.EQ.1) THEN
                 ipmax    = 0
                 aditot_s = aditot_s + adimax_s
              END IF
           END DO

           IF (write_integrate == 1) THEN
              epstot_check = 0
              DO m_cl=1,multra
                 epspar_check(m_cl) = coeps*bigint_s(m_cl)*p_bm2/p_bm2gv**2
                 epstot_check = epstot_check + epspar_check(m_cl) 
              ENDDO
              WRITE(w_u5,format220)                                         &
                   DBLE(nfl*nfp_rat + n),                                   &
                   epstot_check,                                            &
                   p_bm3ge,p_i(npart)/p_bm2,aditot_s/p_bm2
           END IF
       
        END DO

        y2_s      = p_bm2
        y3_s      = p_bm2gv
        y4_s      = p_bm3ge
        y3npart_s = p_i(npart)

        bigint  = bigint + bigint_s
        aditot  = aditot + aditot_s
        y2      = y2 + y2_s
        y3      = y3 + y3_s
        y4      = y4 + y4_s
        y3npart = y3npart + y3npart_s
     END DO
     n = nfp_rat * (nfl_rat + 1)
  END IF
  IF (write_integrate == 1) THEN
     CLOSE(unit=w_u5)
  END IF
  IF (write_diagnostic == 1) THEN
     CLOSE(unit=w_u7)
     OPEN(unit=w_u7,file='diagnostic_add.dat',status='replace',form='formatted')
     WRITE(w_u7,'(4(1x,i8),6(1x,e20.10))')      &
          psi_ind,istepc,npart,max_class,b_min,b_max,bmref,coeps,y2,y3
     CLOSE(unit=w_u7)
  END IF
! **********************************************************************
! Final results
! **********************************************************************
  epstot=0
  DO m_cl=1,multra
     epspar(m_cl) = coeps*bigint(m_cl)*y2/y3**2
     epstot = epstot + epspar(m_cl) 
  ENDDO
!
  ctrone = aditot / y2
  ctrtot = y3npart / y2 
!
  bareph=(pi*ctrone)**2/8
  barept=(pi*ctrtot)**2/8
!
  drdpsi=y2/y3
!
  yps = y4 * j_iota_i
! **********************************************************************
! Log File
! **********************************************************************
  IF (w_u6_open .EQ. 0) THEN
     call safe_open(w_u6, istat, 'neolog.' // trim(extension), 'replace', 'formatted')
     w_u6_open = 1
  ELSE
     OPEN(unit=w_u6,file='neolog.' // trim(extension), status='old',   &
          iostat=istat,position='append',form='formatted')
     if (istat .ne. 0) stop 'NEOLOG.DAT CANNOT BE OPENED IN NEO FLINT_BO'
  END IF
  WRITE(w_u6,'(7(1x,i6),1x,1(d16.8))')                                 &
       psi_ind,n_iota,m_iota,n_gap,nfp_rat,nfl_rat+1,n,epstot 
  CLOSE(w_u6)
! **********************************************************************
! Deallocation
! **********************************************************************
  DEALLOCATE( isw )
  DEALLOCATE( ipa )
  DEALLOCATE( icount )
  DEALLOCATE( iswst )
  DEALLOCATE( eta )

  DEALLOCATE( y )

  DEALLOCATE( bigint )
  DEALLOCATE( bigint_s )
  DEALLOCATE( epspar_check )

  RETURN
END SUBROUTINE flint_bo


SUBROUTINE flint_cur()
! **********************************************************************
! Modules
! **********************************************************************
  USE neo_precision
  USE neo_units
  USE neo_exchange
  USE neo_output
  USE partpa_cur
  USE neo_rk4dcur
! **********************************************************************
! Local definitions
! **********************************************************************
  IMPLICIT NONE
! **********************************************************************
  INTEGER                                          ::  i,n,j1,nfl
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE, TARGET ::  y
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE         ::  y_s
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE         ::  t

  REAL(kind=dp)                                    ::  theta0, phi0
  REAL(kind=dp)                                    ::  theta, phi
  REAL(kind=dp)                                    ::  hphi
  REAL(kind=dp)                                    ::  tmin, tmax, ht

  REAL(kind=dp),               POINTER             ::  p_theta, p_bm2
  REAL(kind=dp),               POINTER             ::  p_bm2gv
  REAL(kind=dp),               POINTER             ::  p_lamps
  REAL(kind=dp),               POINTER             ::  p_lamb1n, p_lamb1d
  REAL(kind=dp),               POINTER             ::  p_lamb2n, p_lamb2d
  REAL(kind=dp), DIMENSION(:), POINTER             ::  p_l, p_k1, p_k
! **********************************************************************
! Allocation
! **********************************************************************
  ndim_cur = npq_cur + 3 * npart_cur
  ALLOCATE( t(npart_cur) )
  ALLOCATE( y_part(npart_cur) )
  ALLOCATE( yfac(npart_cur) )
  ALLOCATE( sqyfac(npart_cur) )
  ALLOCATE( y(ndim_cur) )
  ALLOCATE( y_s(ndim_cur) )
  ALLOCATE( k_fac1(npart_cur) )
  ALLOCATE( k_fac2(npart_cur) )
! **********************************************************************
! Pointer
! **********************************************************************
  p_theta  => y(1)
  p_bm2    => y(2)
  p_bm2gv  => y(3)
  p_lamps  => y(4)
  p_lamb1n => y(5)
  p_lamb1d => y(6)
  p_lamb2n => y(7)
  p_lamb2d => y(8)
  p_l      => y(npq_cur+1:npq_cur+npart_cur)
  p_k1     => y(npq_cur+npart_cur+1:npq_cur+2*npart_cur)
  p_k      => y(npq_cur+2*npart_cur+1:npq_cur+3*npart_cur)
! **********************************************************************
! Initial settings
! **********************************************************************
  ierr   = 0
  hphi   = twopi/(nstep_per*nper)
  bmod0  = bmref
!
  theta0 = theta_bmax
  phi0   = phi_bmax
! **********************************************************************
! Definition of t-values for particles
! **********************************************************************
  tmin = 1.0e-3_dp
  tmax = 1
  ht = (tmax-tmin)/(npart_cur-1)
  DO i=1,npart_cur
     t(i)=tmin+ht*(i-1)
  ENDDO
  y_part = ONE - t**alpha_cur
!!$  PRINT *, 'hphi, ht ',hphi, ht
!!$  PRINT *, 'bmod0 ', bmod0
!!$  PRINT *, 'ndim_cur ', ndim_cur
!!$  PRINT *, 'theta0 ', theta0
!!$  PRINT *, 'phi0 ', phi0  
!!$  PAUSE
!!$  PRINT *, 't ',t
!!$  PAUSE
!!$  PRINT *, 'y_part ',y_part
!!$  PAUSE
! **********************************************************************
  IF (write_cur_inte == 1) THEN
     OPEN(unit=w_u8,file='current.dat',status='replace',form='formatted')
  ENDIF
! **********************************************************************
  phi     = phi0
  y       = 0
  p_theta = theta0
! **********************************************************************
! Integration steps (Summation)
! **********************************************************************
! ATTENTION JUST FOR TEST
  hit_rat = 0
  nintfp  = 10000
! ATTENTION JUST FOR TEST
  IF (hit_rat .EQ. 0) THEN
     DO n=1,nintfp  
        DO j1=1,nstep_per
           CALL RK4D_CUR(phi,y,hphi)
        END DO
        IF (write_cur_inte == 1) THEN
           CALL calccur(n,t,ht,y)
        ENDIF
     END DO
     IF (write_cur_inte == 0) THEN
        CALL calccur(n,t,ht,y)
     END IF
  ELSE
    y_s = 0
    DO nfl = 0, nfl_rat
       phi     = phi0
       y       = 0
       theta   = theta0 + nfl * delta_theta_rat
       p_theta = theta
       DO n=1,nfp_rat  
          DO j1=1,nstep_per
             CALL RK4D_CUR(phi,y,hphi)
          END DO
          IF (write_cur_inte == 1) THEN
             CALL calccur(nfl*nfp_rat+n,t,ht,y)
          ENDIF
       END DO
       y_s = y_s + y
    END DO
    y = y_s
    CALL calccur(nfl*nfp_rat+n+1,t,ht,y)
  END IF
! **********************************************************************
  lambda_b   = lambda_b   / avnabpsi
  lambda_b1  = lambda_b1  / avnabpsi
  lambda_b2  = lambda_b2  / avnabpsi
  lambda_ps1 = lambda_ps1 / avnabpsi
  lambda_ps2 = lambda_ps2 / avnabpsi
! **********************************************************************
  IF (write_cur_inte == 1) THEN
     CLOSE(unit=w_u8)
  ENDIF
! **********************************************************************
! Deallocation
! **********************************************************************
  DEALLOCATE( t )
  DEALLOCATE( y_part )
  DEALLOCATE( yfac )
  DEALLOCATE( sqyfac )
  DEALLOCATE( y )
  DEALLOCATE( y_s )
  DEALLOCATE( k_fac1 )
  DEALLOCATE( k_fac2 )
!
  RETURN
END SUBROUTINE flint_cur


SUBROUTINE calccur(n,t,ht,y)
! **********************************************************************
! Modules
! **********************************************************************
  USE neo_precision
  USE neo_exchange
  USE neo_output
  USE neo_units
  USE partpa_cur
! **********************************************************************
! Local definitions
! **********************************************************************
  IMPLICIT NONE
! **********************************************************************
  REAL(kind=dp), DIMENSION(ndim_cur),  INTENT(in), TARGET ::  y
  REAL(kind=dp), DIMENSION(npart_cur), INTENT(in)         ::  t
  REAL(kind=dp),                       INTENT(in)         ::  ht
  INTEGER                                                 ::  i,n
  REAL(kind=dp),               POINTER                    ::  p_theta, p_bm2
  REAL(kind=dp),               POINTER                    ::  p_bm2gv
  REAL(kind=dp),               POINTER                    ::  p_lamps
  REAL(kind=dp),               POINTER                    ::  p_lamb1n, p_lamb1d
  REAL(kind=dp),               POINTER                    ::  p_lamb2n, p_lamb2d
  REAL(kind=dp), DIMENSION(:), POINTER                    ::  p_l, p_k1, p_k
! **********************************************************************
! Pointer
! **********************************************************************
  p_theta  => y(1)
  p_bm2    => y(2)
  p_bm2gv  => y(3)
  p_lamps  => y(4)
  p_lamb1n => y(5)
  p_lamb1d => y(6)
  p_lamb2n => y(7)
  p_lamb2d => y(8)
  p_l      => y(npq_cur+1:npq_cur+npart_cur)
  p_k1     => y(npq_cur+npart_cur+1:npq_cur+2*npart_cur)
  p_k      => y(npq_cur+2*npart_cur+1:npq_cur+3*npart_cur)
! **********************************************************************
! Initial settings
! **********************************************************************
  lambda_b = 0
! **********************************************************************
  avnabpsi = p_bm2gv / p_bm2
  DO i = 1, npart_cur
     lambda_b = lambda_b + t(i)**(alpha_cur-1) * y_part(i)*y_part(i) *   &
                p_k(i) / p_l(i)
  END DO
  lambda_b   = - 3.0_dp / 8.0_dp * lambda_b * alpha_cur * ht
!
  lambda_ps1 = 2.0_dp * p_lamb1n / p_lamb1d
  lambda_ps2 = 2.0_dp * p_lamb2n / p_lamb2d
!
  lambda_b1  = lambda_b + lambda_ps1
  lambda_b2  = lambda_b + lambda_ps2
! **********************************************************************
  IF (write_cur_inte .EQ. 1) THEN
     WRITE(w_u8,'(1(1x,i8),23(1x,e17.10))')                            &
                 n,                                                    & ! 1   
                 avnabpsi,lambda_b,                                    & ! 2-3
                 lambda_ps1,lambda_ps2,                                & ! 4-5
                 lambda_b1,lambda_b2,                                  & ! 6-7
                 p_lamps,p_lamb1n,p_lamb1d,p_lamb2n,p_lamb2d,          & ! 8-12
                 p_k(1)/p_l(1), p_l(1),p_k(1),p_k1(1),                 & ! 13-16
                 p_k(50)/p_l(50), p_l(50),p_k(50),p_k1(50),            & ! 17-20
                 p_k(100)/p_l(100), p_l(100),p_k(100),p_k1(100)          ! 21-24
  END IF
END SUBROUTINE calccur


SUBROUTINE neo_init(npsi)
! Initialization Routine
! **********************************************************************
! Modules
! **********************************************************************
  USE neo_precision
  USE neo_input
  USE neo_work
  USE neo_exchange
  USE neo_units
  USE neo_parameters
  USE neo_control
  USE neo_spline
! **********************************************************************
! Local Definitions
! **********************************************************************
  IMPLICIT NONE
  INTEGER, INTENT(out)       :: npsi
  INTEGER                    :: imn
! **********************************************************************
! Read input from data file and allocate necessary arrays
! **********************************************************************
  IF (write_progress .NE. 0) WRITE (w_us,*) 'before neo_read'
  CALL neo_read
  IF (write_progress .NE. 0) WRITE (w_us,*) 'after  neo_read'
! **********************************************************************
  npsi = ns
! **********************************************************************
! Allocate and prepare necessary arrays
! **********************************************************************
  IF (write_progress .NE. 0) WRITE (w_us,*) 'before neo_prep'
  CALL neo_prep
  IF (write_progress .NE. 0) WRITE (w_us,*) 'after  neo_prep'
! **********************************************************************
! Calculation of rt0 and bmref
! **********************************************************************
  rt0=0
  bmref=0
  DO imn=1,mnmax
     IF(ixm(imn).EQ.0 .AND. ixn(imn).EQ.0) THEN
        rt0 = rmnc(1,imn)
        bmref = bmnc(1,imn)

        rt0_g = rt0
        bmref_g = bmref
     ENDIF
  ENDDO
  IF(rt0.EQ.ZERO .OR. bmref.EQ.ZERO) THEN
    WRITE (w_us,*) ' NEO_INIT: Fatal problem setting rt0 or bmref'
    STOP
  ENDIF
!
  nper = nfp
! **********************************************************************
  w_u6_open = 0
! **********************************************************************
  RETURN
END SUBROUTINE neo_init


SUBROUTINE neo_init_s(psi,dpsi)
! Initialization for Specific Magnetic Surface
! **********************************************************************
! Modules
! **********************************************************************
  USE neo_precision
  USE neo_input
  USE neo_work
  USE neo_exchange
  USE neo_units
  USE neo_parameters
  USE neo_control
  USE neo_spline
! **********************************************************************
! Local Definitions
! **********************************************************************
  IMPLICIT NONE

  REAL(kind=dp),                INTENT(out)     :: psi, dpsi 
  INTEGER,       DIMENSION(2)                   :: b_minpos, b_maxpos 

  REAL(kind=dp), PARAMETER                      :: eps_newt = 1.0e-10_dp
  INTEGER                                       :: iter, error
  INTEGER,       PARAMETER                      :: iterma_newt = 100
  REAL(kind=dp)                                 :: gval_bmin
  REAL(kind=dp)                                 :: kval_bmin,pval_bmin
  REAL(kind=dp)                                 :: gval_bmax
  REAL(kind=dp)                                 :: kval_bmax,pval_bmax
  REAL(kind=dp)                                 :: qval_bmin,qval_bmax


  REAL(kind=dp)  :: tht, pht, f, g, dfdx, dfdy, dgdx, dgdy
  REAL(kind=dp)  :: thi, phi
  REAL(kind=dp)  :: iot,bval,gval,kval,pval
  INTEGER        :: i
! **********************************************************************
! Calculate Fourier sums and derived quantities
! **********************************************************************
  IF (write_progress .NE. 0) WRITE (w_us,*) 'before neo_fourier'
  CALL neo_fourier
  IF (write_progress .NE. 0) WRITE (w_us,*) 'after  neo_fourier'
! **********************************************************************
! Calculation of dpsi
! **********************************************************************
  psi = es(psi_ind)
  IF (psi_ind .EQ. 1) THEN
     dpsi = 0
  ELSE
     dpsi = psi - es(psi_ind-1)
  ENDIF
! **********************************************************************
! Initilaze spline arrays
! **********************************************************************
! IF (write_progress .NE. 0) WRITE (w_us,*) 'before neo_spline2d'
  CALL neo_spline2d
! IF (write_progress .NE. 0) WRITE (w_us,*) 'after  neo_spline2d'
! **********************************************************************
! Calculate absolute minimum and maximum of b and its location (theta, phi)
! **********************************************************************
  b_minpos   = MINLOC(b)
  b_min      = b(b_minpos(1),b_minpos(2))
  theta_bmin = theta_arr(b_minpos(1))
  phi_bmin   = phi_arr(b_minpos(2))

  b_maxpos   = MAXLOC(b)
  b_max      = b(b_maxpos(1),b_maxpos(2))
  theta_bmax = theta_arr(b_maxpos(1))
  phi_bmax   = phi_arr(b_maxpos(2))

! IF (write_progress .NE. 0) WRITE (w_us,*) 'before neo_zeros2d (min)'
  CALL neo_zeros2d(theta_bmin, phi_bmin, eps_newt, iterma_newt, iter, error)
  CALL neo_eval(theta_bmin,phi_bmin,b_min,gval_bmin,kval_bmin,&
                pval_bmin,qval_bmin)
! IF (write_progress .NE. 0) WRITE (w_us,*) 'after  neo_zeros2d'

! IF (write_progress .NE. 0) WRITE (w_us,*) 'before neo_zeros2d (max)'
  CALL neo_zeros2d(theta_bmax, phi_bmax, eps_newt, iterma_newt, iter, error)
  CALL neo_eval(theta_bmax,phi_bmax,b_max,gval_bmax,kval_bmax,&
                pval_bmax,qval_bmax)
! IF (write_progress .NE. 0) WRITE (w_us,*) 'after  neo_zeros2d'
! **********************************************************************
! ATTENTION
! Set bmref to the absolute maximum of b on flux surface
! This is absolutely necessary for internal routines
! Rescaling is done at the end of the main program neo.f90
! **********************************************************************
  bmref = b_max
  RETURN
END SUBROUTINE neo_init_s


SUBROUTINE neo_read_control
! Read Control File
!***********************************************************************
! Modules
!***********************************************************************
  USE neo_units
  USE neo_control
  USE neo_input
  USE neo_exchange
  USE sizey_bo
  USE sizey_cur
  USE safe_open_mod                                      !LPK/SPH
!***********************************************************************
! Local definitions
!***********************************************************************
  IMPLICIT NONE
  CHARACTER(1) :: dummy
  INTEGER      :: istat
  EXTERNAL GETCARG
!***********************************************************************
! Open input-unit and read data. Check for input file extension on command line
! This no longer supports SEPARATE neo_params.in and neo_in.ext files.
! If neo_params.ext or neo_params.in exists, use it first (check in that order).
! The PREFERRED new file name is neo_in.ext, which REPLACES the neo_params file
! and will work in a multi-tasking environment.
!***********************************************************************
  CALL GETCARG(1, arg1, numargs)                        ! LPK/SPH
  if (numargs .gt. 0) then
     extension = trim(arg1)                             ! LPK/SPH
! First look for new-style neo_in.ext (or neo_param.ext) control file 
! This style is needed for NEO to function correctly in a multi-tasking environment
     arg1 = trim("neo_param." // extension)             ! Prefer neo_in.ext, but this is older style
  else
     arg1 = "neo.in"
  end if

  call safe_open(r_u1, istat, arg1, 'old', 'formatted') ! SPH

! First, check for old-style neo_param.extension and then new-style neo_in.ext control file
  if (istat .ne. 0) then
     arg1 = "neo_param.in" 
     call safe_open(r_u1, istat, arg1, 'old', 'formatted')
     if (istat .ne. 0) then
        arg1 = trim("neo_in." // extension)
        call safe_open(r_u1, istat, arg1, 'old', 'formatted')
     end if
  end if
  if (istat .ne. 0) stop 'NEO control file (neo_param/neo.in) cannot be opened' 

  READ (r_u1,*) dummy
  READ (r_u1,*) dummy
  READ (r_u1,*) dummy
  READ (r_u1,*) in_file
  READ (r_u1,*) out_file
  READ (r_u1,*,iostat=istat) no_fluxs
  if (istat .ne. 0) stop 'Error reading NEO control file'
  IF (no_fluxs .LE. 0) THEN
     READ (r_u1,*) dummy
  ELSE
     ALLOCATE ( fluxs_arr(no_fluxs), stat=istat )
     if (istat .ne. 0) stop 'NEO ALLOCATION ERROR'
     fluxs_arr = 0                                     ! SPH, if user supplied arrary not long enough
     READ (r_u1,*,iostat=istat) fluxs_arr
  END IF

  READ (r_u1,*) theta_n
  READ (r_u1,*) phi_n
  READ (r_u1,*) max_m_mode
  READ (r_u1,*) max_n_mode
  READ (r_u1,*) npart 
  READ (r_u1,*) multra
  READ (r_u1,*) acc_req
  READ (r_u1,*) no_bins
  READ (r_u1,*) nstep_per
  READ (r_u1,*) nstep_min
  READ (r_u1,*) nstep_max
  READ (r_u1,*) calc_nstep_max
  READ (r_u1,*) eout_swi
  READ (r_u1,*) lab_swi
  READ (r_u1,*) inp_swi
  READ (r_u1,*) ref_swi
  READ (r_u1,*) write_progress
  READ (r_u1,*) write_output_files
  READ (r_u1,*) spline_test
  READ (r_u1,*) write_integrate
  READ (r_u1,*) write_diagnostic
  READ (r_u1,*) dummy
  READ (r_u1,*) dummy
  READ (r_u1,*) dummy
  READ (r_u1,*) calc_cur
  READ (r_u1,*) cur_file
  READ (r_u1,*) npart_cur
  READ (r_u1,*) alpha_cur
  READ (r_u1,*) write_cur_inte
 
  CLOSE (unit=r_u1)
! **********************************************************************
  RETURN

END SUBROUTINE neo_read_control


SUBROUTINE neo_read
! Read in Boozer Coordinate and Magnetic Data
!***********************************************************************
! Modules
!***********************************************************************
  USE neo_input
  USE neo_units
  USE neo_control
  USE neo_work
  USE neo_exchange
  use safe_open_mod                                   ! SPH
!***********************************************************************
! Local definitions
!***********************************************************************
  IMPLICIT NONE

  INTEGER :: i,j,j_m,j_n,istat
  INTEGER :: m,n,num_m,num_n,m_found,n_found
  INTEGER :: i_alloc
  INTEGER :: id1,id2,id3,id4,id5,id6,id7
  INTEGER :: k, ns_all                                      ! LPK (used to pack arrays)
  CHARACTER(5) :: dummy
  CHARACTER(45) :: cdum
!***********************************************************************
! Open input-unit and read first quantities
!***********************************************************************
  IF( inp_swi .eq. 0 ) THEN
     CALL READ_BOOZ_IN

  ELSE IF (inp_swi .EQ. 1) THEN            !PPPL-style Bmns file
     call safe_open(r_u1, istat, trim(in_file)//"."//extension, 'old',    &
         'formatted')
     IF (istat .ne. 0) THEN
        write (w_us, *) 'IN_FILE.' // trim(extension) //' NOT FOUND IN NEO_READ'
        stop
     END IF
     READ (r_u1,*) dummy
     READ (r_u1,*) m0b,n0b,ns,nfp,flux
     m_max = m0b+1
     n_max = 2*n0b+1
     mnmax = m_max*n_max

     ns_all = ns                                         ! LPK
     if (no_fluxs>0 .and. no_fluxs<=ns ) ns = no_fluxs    ! LPK/SPH

     if (max_m_mode .le. 0) max_m_mode = m0b             ! SPH
     if (max_n_mode .le. 0) max_n_mode = n0b * nfp       ! SPH
     
! **********************************************************************
! Allocate storage arrays
! **********************************************************************
     ALLOCATE(ixm(mnmax), ixn(mnmax), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'

     ALLOCATE(pixm(mnmax), pixn(mnmax), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays pointers failed!'

     ALLOCATE(i_m(m_max), i_n(n_max), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'

     ALLOCATE(es(ns), iota(ns), curr_pol(ns), curr_tor(ns),               & 
          pprime(ns), sqrtg00(ns), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for real arrays failed!'

     ALLOCATE(rmnc(ns,mnmax), zmnc(ns,mnmax), lmnc(ns,mnmax),             &
          bmnc(ns,mnmax),                                                 &
          stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for fourier arrays (1) failed!'
!***********************************************************************
! Read input arrays
! Only store REQUIRED arrays (ns), but need to scan all of them          ! LPK
!***********************************************************************
     DO k = 1, ns_all                                                    ! LPK
        IF (allocated (fluxs_arr)) THEN                                  ! SPH
           i = 0
           DO j = 1, no_fluxs
              IF (fluxs_arr(j) .eq. k) THEN
                  i = j
                  exit
              END IF
           END DO
        ELSE
           i = k
        END IF                                                           ! SPH END IF
           
        IF (i .gt. 0) THEN
           READ(r_u1,*) dummy
           READ(r_u1,*) es(i),iota(i),curr_pol(i),curr_tor(i),               &
                pprime(i),sqrtg00(i)
           READ(r_u1,*) dummy
           DO j=1,mnmax
              READ(r_u1,*) ixm(j),ixn(j),                                    &
                   rmnc(i,j),zmnc(i,j),lmnc(i,j),                            &
                   bmnc(i,j)
           END DO
        
        ELSE
           DO j=1,mnmax+3                             ! LPK
              READ(r_u1,*) dummy                      ! LPK
           END DO                                     ! LPK
        
        END IF   

     END DO                                                             ! LPK NS_ALL LOOP

  ELSE IF (inp_swi .EQ. 2) THEN        !ORNL-style Boozer file
     call safe_open(r_u1, istat, trim(in_file), 'old', 'formatted')
     IF (istat .ne. 0) THEN
        write (w_us, *) trim(in_file) // ' NOT FOUND IN NEO_READ'
        stop
     END IF

     READ(r_u1,'(a)') cdum
     READ(r_u1,'(a)') cdum
     READ(r_u1,'(a45,10i5)') cdum,ns,id1,id2,id3,id4,m0b,n0b,id5,nfp,mnmax
     READ(r_u1,'(a44,10i5)') cdum,id1,id2,id3,id4,id5,id6,id7

     ns_all = ns                                         ! LPK
     if (no_fluxs>0 .and. no_fluxs<=ns ) ns = no_fluxs    ! LPK/SPH
     if (max_m_mode .le. 0) max_m_mode = m0b
     if (max_n_mode .le. 0) max_n_mode = n0b * nfp

     DO i=1,4
        READ(r_u1,'(a)')cdum
     END DO
     m_max = m0b+1
     n_max = 2*n0b+1
!     mnmax = m_max*n_max - n0b
! **********************************************************************
! Allocate storage arrays
! **********************************************************************
     ALLOCATE(ixm(mnmax), ixn(mnmax), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'

     ALLOCATE(pixm(mnmax), pixn(mnmax), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays pointers failed!'

     ALLOCATE(i_m(m_max), i_n(n_max), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'

     ALLOCATE(es(ns), iota(ns), curr_pol(ns), curr_tor(ns),               & 
          pprime(ns), sqrtg00(ns), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for real arrays failed!'

     ALLOCATE(rmnc(ns,mnmax), zmnc(ns,mnmax), lmnc(ns,mnmax),             &
          bmnc(ns,mnmax),                                                 &
          stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for fourier arrays (1) failed!'
!***********************************************************************
! Read input arrays
!***********************************************************************
     DO  k = 1,ns_all
        IF (allocated (fluxs_arr)) THEN                                  ! SPH
           i = 0
           DO j = 1, no_fluxs
              IF (fluxs_arr(j) .eq. k) THEN
                  i = j
                  exit
              END IF
           END DO
        ELSE
           i = k
        END IF                                                           ! SPH END IF
           
        IF (i .gt. 0) THEN
           READ(r_u1,'(a)') cdum
           READ(r_u1,'(5e12.4)') es(i),iota(i),curr_pol(i),curr_tor(i),flux
           pprime(i) = 0; sqrtg00(i) = 0
           READ(r_u1,'(a)') cdum
           READ(r_u1,"(2i5,1p,4e16.8)") (ixm(j),ixn(j),rmnc(i,j),zmnc(i,j),  &
              lmnc(i,j),bmnc(i,j),j=1,mnmax)
           READ(r_u1,'(a)') cdum
        ELSE
           DO j = 1, mnmax+4
              READ(r_u1,'(a)') dummy                      ! SPH-EAT LINES ON UNWANTED SURFACES
           END DO
        END IF
     END DO
!   do i=1,ns                              !DAS test
!    curr_pol(i) = curr_pol(i)*TWOPI/nfp
!    curr_tor(i) = curr_tor(i)*TWOPI
!   end do
  ELSE
     WRITE (w_us,*) 'FATAL: There is yet no other input type defined'
     STOP
  END IF
!
! Filling of i_m and i_n 
! and pointers pixm from ixm to i_m, and pixn from ixn to i_n
  DO j = 1,mnmax
     m = ixm(j)
     n = ixn(j)
     IF (j .EQ. 1) THEN
        num_m = 1
        i_m(num_m) = m
        pixm(j) = num_m
        num_n = 1
        i_n(num_n) = n
        pixn(j) = num_n
     ELSE
        m_found = 0
        DO j_m = 1, num_m
           IF (m .EQ. i_m(j_m)) THEN
              pixm(j) = j_m
              m_found = 1
           END IF
        END DO
        IF (m_found .EQ. 0) THEN
           num_m = num_m + 1
           i_m(num_m) = m
           pixm(j) = num_m
        END IF
        n_found = 0
        DO j_n = 1, num_n
           IF (n .EQ. i_n(j_n)) THEN
              pixn(j) = j_n
              n_found = 1
           END IF
        END DO
        IF (n_found .EQ. 0) THEN
           num_n = num_n + 1
           i_n(num_n) = n
           pixn(j) = num_n
        END IF
     END IF
  END DO
  IF (lab_swi .EQ. 1) THEN
! 
! ATTENTION: Switch n TO -n
!            Toroidal mode numbers have to multiplied by number of field periods
!            Change iota to iota*nfp (PRINCETON)
!            it is named iota but actually it is iotabar
     ixn = - ixn * nfp
     i_n = - i_n * nfp
     max_n_mode = max_n_mode * nfp
     iota = iota*nfp
  ELSE IF  (lab_swi .EQ. 2) THEN         !ORNL Boozer file
     DO i=1,ns
        DO j=1,mnmax
           lmnc(i,j) = -nfp*lmnc(i,j)/TWOPI   
        END DO
     END DO
  ELSE IF (lab_swi .NE. 0) THEN         
     WRITE(w_us,*) 'LAB_SWI = ', lab_swi,' NOT HANDLED YET!'
     STOP
  END IF
!
! ATTENTION THIS IS JUST FOR TESTING
!
! For scaling of B change the following three with the same factor
! eps_eff should then stay unchanged if the reference for B and R is the same!
!
! bmnc = bmnc * 2.0_dp
! curr_pol = curr_pol * 2.0_dp
! curr_tor = curr_tor * 2.0_dp
!
! For scaling of R change the following four with the same factor
! eps_eff should then stay unchanged if the reference for B and R is the same!
!
! rmnc = rmnc * 2.0_dp
! zmnc = zmnc * 2.0_dp
! curr_pol = curr_pol * 2.0_dp
! curr_tor = curr_tor * 2.0_dp
! 
  CLOSE (unit=r_u1)
! **********************************************************************
! Write optional output for Plotting
! **********************************************************************
  IF (write_output_files .NE. 0) THEN
     IF (write_progress .NE. 0) WRITE (w_us,*) 'write dimension.dat'
     OPEN(unit=w_u1,file='dimension.dat',status='replace',form='formatted')
     WRITE (w_u1,*) ns
     WRITE (w_u1,*) mnmax
     WRITE (w_u1,*) nfp
     WRITE (w_u1,*) theta_n
     WRITE (w_u1,*) phi_n
     WRITE (w_u1,*) s_ind_in
     CLOSE(unit=w_u1)
  ENDIF
  IF (write_output_files .NE. 0) THEN
     IF (write_progress .NE. 0) WRITE (w_us,*) 'write es_arr.dat'
     OPEN(unit=w_u1,file='es_arr.dat',status='replace',form='formatted')
     DO j=1,ns
        WRITE(w_u1,format220) es(j),iota(j),curr_pol(j),curr_tor(j),       &
             pprime(j),sqrtg00(j)
     ENDDO
     CLOSE(unit=w_u1)

     IF (write_progress .NE. 0) WRITE (w_us,*) 'write mn_arr.dat'
     OPEN(unit=w_u1,file='mn_arr.dat',status='replace',form='formatted')
     DO j = 1,mnmax
        WRITE(w_u1,*)ixm(j),ixn(j)
     ENDDO
     CLOSE(unit=w_u1)

     IF (write_progress .NE. 0) WRITE (w_us,*) 'write rmnc_arr.dat'
     OPEN(unit=w_u1,file='rmnc_arr.dat',status='replace',form='formatted')
     DO i=1,ns
        DO j=1,mnmax
           WRITE(w_u1,*) rmnc(i,j)
        END DO
     END DO
     CLOSE(unit=w_u1)

     IF (write_progress .NE. 0) WRITE (w_us,*) 'write zmnc_arr.dat'
     OPEN(unit=w_u1,file='zmnc_arr.dat',status='replace',form='formatted')
     DO i=1,ns
        DO j=1,mnmax
           WRITE(w_u1,*) zmnc(i,j)
        END DO
     END DO
     CLOSE(unit=w_u1)

     IF (write_progress .NE. 0) WRITE (w_us,*) 'write lmnc_arr.dat'
     OPEN(unit=w_u1,file='lmnc_arr.dat',status='replace',form='formatted')
     DO i=1,ns
        DO j=1,mnmax
           WRITE(w_u1,*) lmnc(i,j)
        END DO
     END DO
     CLOSE(unit=w_u1)

     IF (write_progress .NE. 0) WRITE (w_us,*) 'write bmnc_arr.dat'
     OPEN(unit=w_u1,file='bmnc_arr.dat',status='replace',form='formatted')
     DO i=1,ns
        DO j=1,mnmax
           WRITE(w_u1,*) bmnc(i,j)
        END DO
     END DO
     CLOSE(unit=w_u1)
  ENDIF
! **********************************************************************
  RETURN
END SUBROUTINE neo_read


SUBROUTINE neo_prep
! Preparation of Arrays
!***********************************************************************
! Modules
!***********************************************************************
  USE neo_input
  USE neo_work
  USE neo_exchange
  USE neo_parameters
  USE neo_control
  USE neo_units
  USE neo_spline
!***********************************************************************
! Local definitions
!***********************************************************************
  IMPLICIT NONE

  INTEGER :: i_alloc
  INTEGER :: imn, ip, it, k, j
  INTEGER :: ixm_i, ixn_i
  INTEGER :: im, in
  INTEGER :: m, n
! **********************************************************************
! Allocate Storage Arrays
! **********************************************************************
  ALLOCATE(cosmth(theta_n,m_max),                                      &
           sinmth(theta_n,m_max),                                      &
           cosnph(phi_n,  n_max),                                      &
           sinnph(phi_n,  n_max),                                      &
           stat = i_alloc)
  IF(i_alloc /= 0) STOP 'Allocation for cos/sin-arrays failed!'
  ALLOCATE(theta_arr(theta_n),                                         &
           phi_arr(phi_n),                                             &
           stat = i_alloc)
  IF(i_alloc /= 0) STOP 'Allocation for theta/phi-arrays failed!'
! **********************************************************************
! Allocation for arrays for output quantities 
! **********************************************************************
  ALLOCATE(b(theta_n,phi_n),stat = i_alloc)
  IF(i_alloc /= 0) STOP 'Allocation for b-array failed!'
  ALLOCATE(sqrg11(theta_n,phi_n),stat = i_alloc)
  IF(i_alloc /= 0) STOP 'Allocation for sqrg11-array failed!'
  ALLOCATE(kg(theta_n,phi_n),stat = i_alloc)
  IF(i_alloc /= 0) STOP 'Allocation for kg-array failed!'
  ALLOCATE(pard(theta_n,phi_n),stat = i_alloc)
  IF(i_alloc /= 0) STOP 'Allocation for pard-array failed!'
  IF (calc_cur .EQ. 1) THEN
     ALLOCATE(bqtphi(theta_n,phi_n),stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for bqtphi-array failed!'
  END IF
! **********************************************************************
! Allocation for spline arrays
! **********************************************************************
  ALLOCATE(b_spl(4,4,theta_n,phi_n),                                   &
           k_spl(4,4,theta_n,phi_n),                                   &
           g_spl(4,4,theta_n,phi_n),                                   &
           p_spl(4,4,theta_n,phi_n),                                   &
           stat = i_alloc)
  IF(i_alloc /= 0) STOP 'Allocation for spline-arrays failed!'
  IF (calc_cur .EQ. 1) THEN
     ALLOCATE(q_spl(4,4,theta_n,phi_n),                                 &
              stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for q_spl--array failed!'
  END IF
! **********************************************************************
! Some initial work
! **********************************************************************
  theta_start = 0
  theta_end   = twopi
  theta_int   = (theta_end-theta_start)/(theta_n-1)
  phi_start   = 0
  phi_end     = twopi / nfp
  phi_int     = (phi_end-phi_start)/(phi_n-1)
! **********************************************************************
! Preparation of arrays
! **********************************************************************
  DO it=1,theta_n
    theta_arr(it) = theta_start + theta_int*(it-1)
  ENDDO

  DO ip=1,phi_n
    phi_arr(ip) = phi_start + phi_int*(ip-1)
  ENDDO

  DO im = 1,m_max
     m = i_m(im)
     IF (ABS(m) .LE. max_m_mode) THEN
        DO it=1,theta_n
           sinmth(it,im) = SIN( m * theta_arr(it) )
           cosmth(it,im) = COS( m * theta_arr(it) )
        ENDDO
     END IF
  ENDDO
  DO in = 1,n_max
     n = i_n(in)
     IF (ABS(n) .LE. max_n_mode) THEN
        DO ip=1,phi_n
           sinnph(ip,in) = SIN( n * phi_arr(ip) )
           cosnph(ip,in) = COS( n * phi_arr(ip) )
        ENDDO
     END IF
  ENDDO
! **********************************************************************
! Write optional output
! **********************************************************************
  IF (write_output_files .NE. 0) THEN
     IF (write_progress .NE. 0) WRITE (w_us,*) 'write theta_arr.dat'
     OPEN(unit=w_u1,file='theta_arr.dat',status='replace',form='formatted')
     DO j=1,theta_n
        WRITE(w_u1,*) theta_arr(j)
     ENDDO
     CLOSE(unit=w_u1)

     IF (write_progress .NE. 0) WRITE (w_us,*) 'write phi_arr.dat'
     OPEN(unit=w_u1,file='phi_arr.dat',status='replace',form='formatted')
     DO k=1,phi_n
        WRITE(w_u1,*) phi_arr(k)
     ENDDO
     CLOSE(unit=w_u1)
  ENDIF
! **********************************************************************
  RETURN
END SUBROUTINE neo_prep


SUBROUTINE neo_fourier
! Summation of Fourier Sums and Computation of Derived Quantities
!***********************************************************************
! Modules
!***********************************************************************
  USE neo_precision
  USE neo_input
  USE neo_work
  USE neo_exchange
  USE neo_parameters
  USE neo_control
  USE neo_units
!***********************************************************************
! Local definitions
!***********************************************************************
  IMPLICIT NONE

  INTEGER       :: i_alloc
  INTEGER       :: im, in, m, n, i, j
  INTEGER       :: it, ip, imn
  REAL(kind=dp) :: ri, zi, li, bi
  REAL(kind=dp) :: cosv, sinv
!***********************************************************************
! Allocation of arrays
!***********************************************************************
!  write(*,*) theta_n,phi_n
  ALLOCATE(r(theta_n,phi_n),                                           &
           z(theta_n,phi_n),                                           &
           l(theta_n,phi_n),                                           &
           stat = i_alloc)
  IF(i_alloc /= 0) STOP 'Allocation for sum-arrays failed!'
  ALLOCATE(r_tb(theta_n,phi_n),                                        &
           z_tb(theta_n,phi_n),                                        &
           p_tb(theta_n,phi_n),                                        &
           b_tb(theta_n,phi_n),                                        &
           stat = i_alloc)
  IF(i_alloc /= 0) STOP 'Allocation for sum-arrays failed!'
  ALLOCATE(r_pb(theta_n,phi_n),                                        &
           z_pb(theta_n,phi_n),                                        &
           p_pb(theta_n,phi_n),                                        &
           b_pb(theta_n,phi_n),                                        &
           stat = i_alloc)
  IF(i_alloc /= 0) STOP 'Allocation for sum-arrays failed!'
  ALLOCATE(gtbtb(theta_n,phi_n),                                       &
           gpbpb(theta_n,phi_n),                                       &
           gtbpb(theta_n,phi_n),                                       &
           isqrg(theta_n,phi_n),                                       &
           stat = i_alloc)
  IF(i_alloc /= 0) STOP 'Allocation for sum-arrays failed!'
!***********************************************************************
! Summation of Fourier components
!***********************************************************************
  r = 0
  z = 0
  l = 0
  b = 0

  r_tb = 0
  z_tb = 0
  p_tb = 0
  b_tb = 0
  r_pb = 0
  z_pb = 0
  p_pb = 0
  b_pb = 0
  
  DO imn=1,mnmax
     ri = rmnc(psi_ind,imn)
     zi = zmnc(psi_ind,imn)
     li = lmnc(psi_ind,imn)
     bi = bmnc(psi_ind,imn)
     m = ixm(imn)
     n = ixn(imn)
     im = pixm(imn)
     in = pixn(imn)
     IF (ABS(m) .LE. max_m_mode .AND. ABS(n) .LE. max_n_mode) THEN
        DO ip=1,phi_n
           DO it=1,theta_n
              cosv = cosmth(it,im) * cosnph(ip,in) + sinmth(it,im) * sinnph(ip,in)
              sinv = sinmth(it,im) * cosnph(ip,in) - cosmth(it,im) * sinnph(ip,in)

              r(it,ip) = r(it,ip) + ri*cosv
              z(it,ip) = z(it,ip) + zi*sinv
              l(it,ip) = l(it,ip) + li*sinv
              b(it,ip) = b(it,ip) + bi*cosv

              r_tb(it,ip) = r_tb(it,ip) - m*ri*sinv
              r_pb(it,ip) = r_pb(it,ip) + n*ri*sinv
              z_tb(it,ip) = z_tb(it,ip) + m*zi*cosv
              z_pb(it,ip) = z_pb(it,ip) - n*zi*cosv
              p_tb(it,ip) = p_tb(it,ip) - m*li*cosv
              p_pb(it,ip) = p_pb(it,ip) + n*li*cosv
              b_tb(it,ip) = b_tb(it,ip) - m*bi*sinv
              b_pb(it,ip) = b_pb(it,ip) + n*bi*sinv
           END DO
        END DO
     END IF
  END DO

! Compute Boozer-theta_b (tb), phi_b (pb) derivatives of the cylindrical phi angle (p)

  p_tb = p_tb * twopi / nfp
  p_pb = ONE + p_pb * twopi / nfp
! **********************************************************************
! Ensure periodicity boundaries to be the same
! **********************************************************************
  r(theta_n,:) = r(1,:)
  r(:,phi_n)   = r(:,1)
  z(theta_n,:) = z(1,:)
  z(:,phi_n)   = z(:,1)
  l(theta_n,:) = l(1,:)
  l(:,phi_n)   = l(:,1)
  b(theta_n,:) = b(1,:)
  b(:,phi_n)   = b(:,1)
  r_tb(theta_n,:) = r_tb(1,:)
  r_tb(:,phi_n)   = r_tb(:,1)
  r_pb(theta_n,:) = r_pb(1,:)
  r_pb(:,phi_n)   = r_pb(:,1)
  z_tb(theta_n,:) = z_tb(1,:)
  z_tb(:,phi_n)   = z_tb(:,1)
  z_pb(theta_n,:) = z_pb(1,:)
  z_pb(:,phi_n)   = z_pb(:,1)
  p_tb(theta_n,:) = p_tb(1,:)
  p_tb(:,phi_n)   = p_tb(:,1)
  p_pb(theta_n,:) = p_pb(1,:)
  p_pb(:,phi_n)   = p_pb(:,1)
  b_tb(theta_n,:) = b_tb(1,:)
  b_tb(:,phi_n)   = b_tb(:,1)
  b_pb(theta_n,:) = b_pb(1,:)
  b_pb(:,phi_n)   = b_pb(:,1)
! **********************************************************************
! Derived quantities
! NOTE: The radial coordinate used in these formulae is psi = (TOROIDAL FLUX)/TWOPI
! so that PHIP == d psi / ds = 1
! **********************************************************************
! metric tensor
  gtbtb = r_tb*r_tb + z_tb*z_tb + r*r*p_tb*p_tb  
  gpbpb = r_pb*r_pb + z_pb*z_pb + r*r*p_pb*p_pb  
  gtbpb = r_tb*r_pb + z_tb*z_pb + r*r*p_tb*p_pb  
! $1/sqrt(g)$
  fac = curr_pol(psi_ind) + iota(psi_ind)*curr_tor(psi_ind)
  isqrg  = b*b / fac
! $sqrt(g^{11})$  == |grad-psi|
  sqrg11 = SQRT( gtbtb*gpbpb - gtbpb*gtbpb ) * isqrg
! geodesic curvature term $k_G |\nabla \psi|$ 
! proportional to v_drift * grad(psi) (radial drift velocity)         !!SPH
  kg = (curr_tor(psi_ind)*b_pb - curr_pol(psi_ind)*b_tb) / fac
! parallel derivative of mod-B
  pard = b_pb + iota(psi_ind)*b_tb
! quasi-toroidal phi component of b (only for parallel current)
  IF (calc_cur .EQ. 1) THEN
     bqtphi = isqrg * (p_pb + iota(psi_ind)*p_tb)
  END IF
! **********************************************************************
! Optional Output
! **********************************************************************
  IF (write_output_files .NE. 0) THEN
     IF (write_progress .NE. 0) WRITE (w_us,*) 'write b_s_arr.dat'
     OPEN(unit=w_u1,file='b_s_arr.dat',status='replace',form='formatted')
     DO i=1,theta_n
        DO j=1,phi_n
           WRITE(w_u1,*)  b(i,j)
        END DO
     ENDDO
     CLOSE(unit=w_u1)

     IF (write_progress .NE. 0) WRITE (w_us,*) 'write r_s_arr.dat'
     OPEN(unit=w_u1,file='r_s_arr.dat',status='replace',form='formatted')
     DO i=1,theta_n
        DO j=1,phi_n
           WRITE(w_u1,*)  r(i,j)
        END DO
     ENDDO
     CLOSE(unit=w_u1)

     IF (write_progress .NE. 0) WRITE (w_us,*) 'write z_s_arr.dat'
     OPEN(unit=w_u1,file='z_s_arr.dat',status='replace',form='formatted')
     DO i=1,theta_n
        DO j=1,phi_n
           WRITE(w_u1,*)  z(i,j)
        END DO
     ENDDO
     CLOSE(unit=w_u1)

     IF (write_progress .NE. 0) WRITE (w_us,*) 'write l_s_arr.dat'
     OPEN(unit=w_u1,file='l_s_arr.dat',status='replace',form='formatted')
     DO i=1,theta_n
        DO j=1,phi_n
           WRITE(w_u1,*)  l(i,j)
        END DO
     ENDDO
     CLOSE(unit=w_u1)

     IF (write_progress .NE. 0) WRITE (w_us,*) 'write isqrg_arr.dat'
     OPEN(unit=w_u1,file='isqrg_arr.dat',status='replace',form='formatted')
     DO i=1,theta_n
        DO j=1,phi_n
           WRITE(w_u1,*)  isqrg(i,j)
        END DO
     ENDDO
     CLOSE(unit=w_u1)

     IF (write_progress .NE. 0) WRITE (w_us,*) 'write sqrg11_arr.dat'
     OPEN(unit=w_u1,file='sqrg11_arr.dat',status='replace',form='formatted')
     DO i=1,theta_n
        DO j=1,phi_n
           WRITE(w_u1,*)  sqrg11(i,j)
        END DO
     ENDDO
     CLOSE(unit=w_u1)

     IF (write_progress .NE. 0) WRITE (w_us,*) 'write kg_arr.dat'
     OPEN(unit=w_u1,file='kg_arr.dat',status='replace',form='formatted')
     DO i=1,theta_n
        DO j=1,phi_n
           WRITE(w_u1,*)  kg(i,j)
        END DO
     ENDDO
     CLOSE(unit=w_u1)

     IF (write_progress .NE. 0) WRITE (w_us,*) 'write pard_arr.dat'
     OPEN(unit=w_u1,file='pard_arr.dat',status='replace',form='formatted')
     DO i=1,theta_n
        DO j=1,phi_n
           WRITE(w_u1,*)  pard(i,j)
        END DO
     ENDDO
     CLOSE(unit=w_u1)
  ENDIF
! **********************************************************************
!  Deallocation of unnecessary arrays
! **********************************************************************
  DEALLOCATE (r,z,l)
  DEALLOCATE (r_tb,z_tb,p_tb,b_tb)
  DEALLOCATE (r_pb,z_pb,p_pb,b_pb)
  DEALLOCATE (gtbtb,gpbpb,gtbpb,isqrg)
! **********************************************************************
  RETURN
END SUBROUTINE neo_fourier


SUBROUTINE neo_spline2d
! Creation of Spline Arrays
! **********************************************************************
! Modules
! **********************************************************************
  USE neo_precision
  USE neo_input
  USE neo_work
  USE neo_exchange
  USE neo_units
  USE neo_parameters
  USE neo_control
  USE neo_spline
! **********************************************************************
! Local Definitions
! **********************************************************************
  IMPLICIT NONE
! **********************************************************************
! Allocation of spline arrays is done in neo_prep
! **********************************************************************
! **********************************************************************
! Double periodic splines (parameter mt=1 and mp=1) 
! **********************************************************************
! Spline for mod b
! IF (write_progress .NE. 0) WRITE (w_us,*) 'before spl2d'
  CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,b,b_spl)
! IF (write_progress .NE. 0) WRITE (w_us,*) 'after spl2d'
! Spline for sqrg11
! IF (write_progress .NE. 0) WRITE (w_us,*) 'before spl2d'
  CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,sqrg11,g_spl)
! IF (write_progress .NE. 0) WRITE (w_us,*) 'after spl2d'
! Spline for geodesic curviture
! IF (write_progress .NE. 0) WRITE (w_us,*) 'before spl2d'
  CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,kg,k_spl)
! IF (write_progress .NE. 0) WRITE (w_us,*) 'after spl2d'
! Spline for parallel derivative
! IF (write_progress .NE. 0) WRITE (w_us,*) 'before spl2d'
  CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,pard,p_spl)
! IF (write_progress .NE. 0) WRITE (w_us,*) 'after spl2d'
! Spline for quasi-toroidal phi component of b
  IF (calc_cur .EQ. 1) THEN
!    IF (write_progress .NE. 0) WRITE (w_us,*) 'before spl2d'
     CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,bqtphi,q_spl)
!    IF (write_progress .NE. 0) WRITE (w_us,*) 'after spl2d'
  END IF
! **********************************************************************
! Spline test
! **********************************************************************
  IF (spline_test .GT. 0) THEN
     IF (write_progress .NE. 0) WRITE (w_us,*) 'before neo_spline_test'
     CALL neo_spline_test
     IF (write_progress .NE. 0) WRITE (w_us,*) 'after neo_spline_test'
  ENDIF
! **********************************************************************
  RETURN
END SUBROUTINE neo_spline2d


SUBROUTINE neo_eval(theta,phi,bval,gval,kval,pval,qval)
! Evaluation of Splines
! **********************************************************************
! Modules
! **********************************************************************
  USE neo_precision
  USE neo_input
  USE neo_work
  USE neo_exchange
  USE neo_units
  USE neo_parameters
  USE neo_control
  USE neo_spline
! **********************************************************************
! Local Definitions
! **********************************************************************
  IMPLICIT NONE

  REAL(kind=dp), INTENT(in)  ::   theta, phi 
  REAL(kind=dp), INTENT(out) ::   bval, gval, kval, pval, qval
! **********************************************************************
! Evaluation of pointer
! **********************************************************************
  CALL poi2d(theta_int,phi_int,mt,mp,                                  &
             theta_start,theta_end,phi_start,phi_end,                  &
             theta,phi,theta_ind,phi_ind,theta_d,phi_d,ierr)
! **********************************************************************
! Evaluation of 2d-splines
! **********************************************************************
  CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,            &
             b_spl,bval)
  CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,            &
             g_spl,gval)
  CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,            &
             k_spl,kval)
  CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,            &
             p_spl,pval)
  IF (calc_cur .EQ. 1) THEN
     CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,         &
                q_spl,qval)
  END IF
  RETURN
END SUBROUTINE neo_eval


SUBROUTINE neo_bderiv(theta, phi, f, g, dfdx, dfdy, dgdx, dgdy)
!
! Calculates first and second derivatives of b using the 2d-splines
!
! Input:  theta, phi
! Output: f      db/dt              (t = theta)
!         g      db/dp              (p = phi)
!         dfdx   d^2b/dt^2
!         dfdy   d^2b/(dt dp)
!         dgdx   d^2b/(dt dp)
!         dgdy   d^2b/dp^2
! 
! Input/output consistent for neo_zeros2d
! **********************************************************************
! Modules
! **********************************************************************
  USE neo_precision
  USE neo_input
  USE neo_work
  USE neo_units
  USE neo_parameters
  USE neo_control
  USE neo_spline
! **********************************************************************
! Local Definitions
! **********************************************************************
  IMPLICIT NONE

  REAL(kind=dp), INTENT(in)    ::   theta, phi 
  REAL(kind=dp), INTENT(out)   ::   f, g, dfdx, dfdy, dgdx, dgdy

  REAL(kind=dp), DIMENSION(2)  ::   fderiv 
  REAL(kind=dp), DIMENSION(3)  ::   sderiv 
! **********************************************************************
! Evaluation of pointer
! **********************************************************************
  CALL poi2d(theta_int,phi_int,mt,mp,                                  &
             theta_start,theta_end,phi_start,phi_end,                  &
             theta,phi,theta_ind,phi_ind,theta_d,phi_d,ierr)
! **********************************************************************
! Evaluation of 2d-splines (first and second derivatives)
! **********************************************************************
!  print *, 'before eva2d_fd'
  CALL eva2d_fd(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,         &
             b_spl,fderiv)
!  print *, 'before eva2d_sd'
  CALL eva2d_sd(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,         &
             b_spl,sderiv)  
!  print *, 'after eva2d_sd'
! **********************************************************************
! Outputvalues (for neo_zeros2d)
! **********************************************************************
  f    = fderiv(1)
  g    = fderiv(2)
  dfdx = sderiv(1)
  dfdy = sderiv(2)
  dgdx = sderiv(2)
  dgdy = sderiv(3)
END SUBROUTINE neo_bderiv


SUBROUTINE neo_dealloc
! **********************************************************************
! Modules
! **********************************************************************
  USE neo_precision
  USE neo_input
  USE neo_work
  USE neo_exchange
  USE neo_units
  USE neo_parameters
  USE neo_control
  USE neo_spline
! **********************************************************************
! Local Definitions
! **********************************************************************
  IMPLICIT NONE
! *******************************************************************
! DeAllocate Storage Arrays
! *******************************************************************
  DEALLOCATE (es,iota,curr_pol,curr_tor,pprime,sqrtg00)
  DEALLOCATE (theta_arr,phi_arr)
  DEALLOCATE (rmnc,zmnc,lmnc,bmnc)
  DEALLOCATE (ixm, ixn)
  DEALLOCATE (pixm, pixn)
  DEALLOCATE (i_m, i_n)
  DEALLOCATE (b_spl,k_spl,g_spl,p_spl)
  DEALLOCATE (b,sqrg11,kg,pard)
  IF (calc_cur .EQ. 1) THEN
     DEALLOCATE (bqtphi)
     DEALLOCATE (q_spl)
  END IF
  DEALLOCATE (fluxs_arr)
! *******************************************************************
  RETURN
END SUBROUTINE neo_dealloc


SUBROUTINE neo_spline_test
! Test of Spline Routine
! **********************************************************************
! Modules
! **********************************************************************
  USE neo_precision
  USE neo_input
  USE neo_work
  USE neo_units
  USE neo_parameters
  USE neo_control
! **********************************************************************
! Local Definitions
! **********************************************************************
  IMPLICIT NONE
! **********************************************************************
  INTEGER            :: i, j, n
  REAL(kind=dp)      :: theta, phi, bval, gval, kval, pval, qval
  INTEGER, PARAMETER :: div = 5, ip = 43, it = 88
  REAL(kind=dp)      :: td, pd 
  n = MIN(theta_n,phi_n)

  IF (spline_test .EQ. 1) THEN ! along given phi
     OPEN(unit=w_u1,file='sptest1.dat',status='replace',form='formatted')
     OPEN(unit=w_u2,file='sptest2.dat',status='replace',form='formatted')
     td = theta_end - theta_start
     phi = phi_arr(ip)
     theta = theta_start-td 
     DO WHILE (theta < theta_end+td)
        CALL neo_eval(theta,phi,bval,gval,kval,pval,qval)
        WRITE(w_u1,*) theta,phi,bval,gval,kval
        theta = theta + theta_int/div
     END DO
     DO j = -1,1
        DO i = 1, theta_n-1
           theta = theta_arr(i) + j*td
           bval = b(i,ip)
           gval = sqrg11(i,ip)
           kval = kg(i,ip)
           WRITE(w_u2,*) theta,phi,bval,gval,kval
        END DO
     END DO
     CLOSE(unit=w_u2)
     CLOSE(unit=w_u1)
  ELSEIF (spline_test .EQ. 2) THEN ! along given theta
     OPEN(unit=w_u1,file='sptest1.dat',status='replace',form='formatted')
     OPEN(unit=w_u2,file='sptest2.dat',status='replace',form='formatted')
     pd = phi_end - phi_start
     theta = theta_arr(it)
     phi = phi_start-pd
     DO WHILE (phi < phi_end+pd)
        CALL neo_eval(theta,phi,bval,gval,kval,pval,qval)
        WRITE(w_u1,*) theta,phi,bval,gval,kval
        phi = phi + phi_int/div
     END DO
     DO j = -1,1
        DO i = 1, phi_n-1
           phi = phi_arr(i) + j*pd
           bval = b(it,i)
           gval = sqrg11(it,i)
           kval = kg(it,i)
           WRITE(w_u2,*) theta,phi,bval,gval,kval
        END DO
     END DO
     CLOSE(unit=w_u2)
     CLOSE(unit=w_u1)
  ELSEIF (spline_test .EQ. 3) THEN ! diagonal
     OPEN(unit=w_u1,file='sptest1.dat',status='replace',form='formatted')
     DO i = 1, n*div
        theta = theta_start + i*theta_int/div
        phi   = phi_start + i*phi_int/div
        CALL neo_eval(theta,phi,bval,gval,kval,pval,qval)
        WRITE(w_u1,*) theta,phi,bval,gval,kval
     END DO
     CLOSE(unit=w_u1)
     OPEN(unit=w_u1,file='sptest2.dat',status='replace',form='formatted')
     DO i = 1, n
        theta = theta_arr(i)
        phi   = phi_arr(i)
        bval = b(i,i)
        gval = sqrg11(i,i)
        kval = kg(i,i)
        WRITE(w_u1,*) theta,phi,bval,gval,kval
     END DO
     CLOSE(unit=w_u1)
  ENDIF
! **********************************************************************
  RETURN
END SUBROUTINE neo_spline_test


SUBROUTINE neo_zeros2d(x, y, eps, iter_ma, iter, error)

  USE neo_precision

  IMPLICIT NONE

  INTEGER,        INTENT(out)   :: error, iter
  INTEGER,        INTENT(in)    :: iter_ma
  REAL (kind=dp), INTENT(in)    :: eps
  REAL (kind=dp), INTENT(inout) :: x, y

  REAL (kind=dp) :: x_n, y_n
  REAL (kind=dp) :: f, dfdx, dfdy, g, dgdx, dgdy
  REAL (kind=dp) :: f_n,g_n
  REAL (kind=dp) :: det
  REAL (kind=dp) :: x_err, y_err

  error = 0

! compute f(x,y), g(x,y) and all first derivatives
  CALL neo_bderiv(x, y, f, g, dfdx, dfdy, dgdx, dgdy)
  DO iter = 1, iter_ma

     det = dfdx * dgdy - dfdy * dgdx
!!$     PRINT *, 'f,    g    ', f, g
!!$     PRINT *, 'dfdx, dgdx ', dfdx, dgdx
!!$     PRINT *, 'dfdy, dgdy ', dfdy, dgdy
!!$     PRINT *, 'det        ', det
     x_n = x + ( dfdy *  g   -  f   * dgdy ) / det
     y_n = y + (  f   * dgdx - dfdx *  g   ) / det
!!$     PRINT *, 'x,    y    ', x, y
!!$     PRINT *, 'x_n,  y_n  ', x_n, y_n

!    compute f(x,y), g(x,y) and all first derivatives
!    at the new positions
!     PRINT *, x_n, y_n
     CALL neo_bderiv(x_n, y_n, f_n, g_n, dfdx, dfdy, dgdx, dgdy)

!    remaining relatve errors
!     IF (x_n .NE. ZERO) THEN
     IF (ABS(x_n) .GT. eps) THEN
       x_err = ABS ( (x_n - x) / x_n )
     ELSE
       x_err = ABS ( x_n - x )
     END IF
!     IF (y_n .NE. ZERO) THEN
     IF (ABS(y_n) .GT. eps) THEN
       y_err = ABS ( (y_n - y) / y_n )
     ELSE
       y_err = ABS ( y_n - y )
     END IF

!    new values
     f = f_n
     g = g_n
     x = x_n
     y = y_n

!    exit if error is small enough
     IF ( MAX ( x_err, y_err ) < eps ) RETURN

  END DO

  error = 1

  RETURN
END SUBROUTINE neo_zeros2d


subroutine read_booz_in

      use read_boozer_mod, dp1 => dp
      use neo_input
      use neo_units
      use neo_control
      use neo_exchange
      implicit none
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real (rprec) ::  hs
      integer :: k, i, j, i_surf, mn0
      integer :: i_alloc
      character*(120) :: booz_input_file
!
!     Read data from the boozmn file and allocate storage:
!
      booz_input_file = extension
     
      call read_boozer_file (booz_input_file,k)
      if (k .ne. 0) stop 'Error reading boozmn file'

      m0b = mboz_b-1
      n0b = nboz_b
      mnmax = mnboz_b

      ns = ns_b
      nfp = nfp_b
      flux = abs( phi_b(ns_b) )
      max_n_mode = max_n_mode * nfp

!     truncate modes
!     assume m=0, mboz
!            n=-nboz, nboz
!     discard m> mb_retained
!             iabs(n)> nb_retained
!     perhaps better to retain largest mnmax modes
!             for a reference surfaces

      if (mnboz_b .ne. (mboz_b-1)*(2*nboz_b+1)+nboz_b+1 )   &
           stop "boozer modes inconsistent"

      k = 0
      if (max_m_mode .le. 0) max_m_mode = m0b
      if (max_n_mode .le. 0) max_n_mode = n0b * nfp
      do j = 1, mnboz_b
        if( ixm_b(j) .gt. max_m_mode ) cycle
        if( iabs(ixn_b(j) ) .gt. max_n_mode ) cycle
        k = k + 1
      enddo

      mnmax = k
      m_max = min0(m0b, max_m_mode)
      m_max = m_max + 1
      n_max = min0(n0b, max_n_mode/nfp)
      n_max = 2*n_max+1
      if ((no_fluxs.gt.0) .and. (no_fluxs.le.ns_b)) ns = no_fluxs

! **********************************************************************
! Allocate storage arrays
! **********************************************************************

      ALLOCATE(ixm(mnmax), ixn(mnmax), stat = i_alloc)
      IF(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'

      ALLOCATE(pixm(mnmax), pixn(mnmax), stat = i_alloc)
      IF(i_alloc /= 0)                                              &
           STOP 'Allocation for integer arrays with pointers failed!'

      ALLOCATE(i_m(m_max), i_n(n_max), stat = i_alloc)
      IF(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'

      ALLOCATE(es(ns), iota(ns), curr_pol(ns), curr_tor(ns),        &      
           pprime(ns), sqrtg00(ns), stat = i_alloc)
      IF(i_alloc /= 0) STOP 'Allocation for real arrays failed!'

      ALLOCATE(rmnc(ns,mnmax), zmnc(ns,mnmax), lmnc(ns,mnmax),      & 
           bmnc(ns,mnmax), stat = i_alloc)
      IF(i_alloc /= 0) STOP 'Allocation for fourier arrays (1) failed!'
! **********************************************************************

      k = 0;   mn0 = 0
      do j = 1, mnboz_b
         if( ixm_b(j) .gt. max_m_mode ) cycle
         if( iabs(ixn_b(j) ) .gt. max_n_mode ) cycle
         if (ixm_b(j).eq.0 .and. ixn_b(j).eq.0) mn0 = j
         k = k + 1
         ixm(k) = ixm_b(j)
         ixn(k) = ixn_b(j)
      enddo

      if (mn0 .eq. 0) stop 'M=0, N=0 MUST BE INCLUDE IN BOOZER SPECTRUM!'
      
      hs = ONE/(ns_b -1)
      do i = 1, ns                  !!NEED TO INCLUDE IN LOOP CASE FOR NO_FLUXS < 0 (SPH)

         if (no_fluxs > 0) then
            i_surf = fluxs_arr(i)   !!NOT ALLOCATED YET IF NO_FLUXS<0
         else
            i_surf = i
         end if
         k = 0
         do j = 1, mnboz_b

           if( ixm_b(j) .gt. max_m_mode ) cycle
           if( abs(ixn_b(j)) .gt. max_n_mode ) cycle
           k = k + 1

           rmnc(i,k) = rmnc_b(j,i_surf)
           zmnc(i,k) = zmns_b(j,i_surf)
!
!       note: zeta - zeta_b from booz_xform,
!             Note: in BOOZMN file, zeta_b - zeta is stored
!             Also, lmnc is normalized to one field period
!

           lmnc(i,k) = -pmns_b(j,i_surf)*nfp_b/twopi
           bmnc(i,k) = bmn_b(j,i_surf)
   
           if (ixm_b(j).eq.0 .and. ixn_b(j).eq.0)            &
               sqrtg00(i) = gmn_b(j,i_surf)

         enddo

         if (rmnc(i,mn0) .eq. ZERO) then
            write (w_us, *)' The surface i = ', i_surf, ' is absent in the BOOZMN file'
            stop
         end if
         
         es(i) = (i_surf-1.5_dp)*hs
         iota(i) = iota_b(i_surf)
         if (i_surf .lt. ns_b) then
            pprime(i) = (pres_b(i_surf+1) - pres_b(i_surf))/hs
         else
            pprime(i) = 0
         end if
!
!       these are the F and I functions, respectively (currents divided by 2*pi, sign)
!
         curr_pol(i) = bvco_b(i_surf)
         curr_tor(i) = buco_b(i_surf)

      enddo

      
      call read_boozer_deallocate
  
END SUBROUTINE READ_BOOZ_IN


SUBROUTINE spl2d(nx,ny,hx,hy,mx,my,f,spl)

! Makes a 2-dimensional cubic spline of function f(x,y)
!
! Input:  nx, ny              number of values in x and y  
!         hx, hy              step size in x and y (aequidistant)
!         mx, my              spline mode (0: standard, 1: periodic)
!         f(nx,ny)            f(x,y)-values
! Output: spl                 Array with spline parameters

  USE neo_precision

  IMPLICIT NONE

  INTEGER,                             INTENT(in)  :: nx, ny, mx, my
  REAL(kind=dp),                       INTENT(in)  :: hx, hy
  REAL(kind=dp), DIMENSION(nx,ny)    , INTENT(in)  :: f
  REAL(kind=dp), DIMENSION(4,4,nx,ny), INTENT(out) :: spl
  
  REAL(kind=dp), DIMENSION(:),         ALLOCATABLE :: bi, ci, di, s
  INTEGER                                          :: i, j, k, l

  ALLOCATE ( bi(nx), ci(nx), di(nx), s(nx) )
  DO j = 1,ny
     DO i = 1,nx
        s(i) = f(i,j)
     END DO
     IF (mx .EQ. 0) THEN
        CALL splreg(nx,hx,s,bi,ci,di)
     ELSE
        CALL splper(nx,hx,s,bi,ci,di)
     ENDIF
     DO i = 1,nx
        spl(1,1,i,j) = s(i)
        spl(2,1,i,j) = bi(i)
        spl(3,1,i,j) = ci(i)
        spl(4,1,i,j) = di(i)
     END DO
  END DO
  DEALLOCATE ( bi, ci, di, s )

  ALLOCATE ( bi(ny), ci(ny), di(ny), s(ny) )
  DO k = 1,4
     DO i = 1,nx
        DO j = 1,ny
           s(j) = spl(k,1,i,j)
        END DO
        IF (my .EQ. 0) THEN
           CALL splreg(ny,hy,s,bi,ci,di)
        ELSE
           CALL splper(ny,hy,s,bi,ci,di)
        ENDIF
        DO j=1,ny
           spl(k,2,i,j)=bi(j)
           spl(k,3,i,j)=ci(j)
           spl(k,4,i,j)=di(j)
        END DO
     END DO
  END DO
  DEALLOCATE ( bi, ci, di, s )

  RETURN
END SUBROUTINE spl2d
!=====================================================
SUBROUTINE eva2d(nx,ny,ix,iy,dx,dy,spl,spval)

! Evaluates a 2-dimensional cubic spline of function f(x,y)
!
! Input:  nx, ny              number of values in x and y  
!         ix, iy              pointer into the spline array spl
!         dx, dy              distance from x(ix) and y(iy)
!         spl                 array with spline data
! Output: spval               evaluated function value

  USE neo_precision
  
  IMPLICIT NONE
  
  INTEGER,                             INTENT(in)  :: nx, ny, ix, iy
  REAL(kind=dp),                       INTENT(in)  :: dx, dy
  REAL(kind=dp), DIMENSION(4,4,nx,ny), INTENT(in)  :: spl
  REAL(kind=dp),                       INTENT(out) :: spval

  REAL(kind=dp), DIMENSION(4)                      :: a
  INTEGER                                          :: l

  DO l=1,4
     a(l) = spl(1,l,ix,iy) + dx*(spl(2,l,ix,iy) +              & 
          dx*(spl(3,l,ix,iy) + dx* spl(4,l,ix,iy)))
  END DO
  spval = a(1)+dy*(a(2)+dy*(a(3)+dy*a(4)))

  RETURN
END SUBROUTINE eva2d
!=====================================================
SUBROUTINE eva2d_fd(nx,ny,ix,iy,dx,dy,spl,spval)

! Evaluates the first derivatives of 2-dimensional cubic spline of function f(x,y)
!
! Input:  nx, ny              number of values in x and y  
!         ix, iy              pointer into the spline array spl
!         dx, dy              distance from x(ix) and y(iy)
!         spl                 array with spline data
! Output: spval(2)            evaluated function value
!                             spval(1) = df/dx
!                             spval(2) = df/dy

  USE neo_precision
  
  IMPLICIT NONE
  
  INTEGER,                             INTENT(in)  :: nx, ny, ix, iy
  REAL(kind=dp),                       INTENT(in)  :: dx, dy
  REAL(kind=dp), DIMENSION(4,4,nx,ny), INTENT(in)  :: spl
  REAL(kind=dp), DIMENSION(2),         INTENT(out) :: spval

  INTEGER                                          :: i,j
  REAL(kind=dp)                                    :: muli, mulj

  spval = 0

! df/dx
  DO i=2,4
     IF (i == 2) THEN
        muli = 1
     ELSE
        muli = dx**(i-2)
     END IF
     muli = muli * (i-1)
     DO j=1,4
       IF (j == 1) THEN
           mulj = 1
        ELSE
           mulj = dy**(j-1)
        END IF
        spval(1) = spval(1) + spl(i,j,ix,iy) * muli * mulj
     END DO
  END DO

! df/dy
  DO i=1,4
     IF (i == 1) THEN
        muli = 1
     ELSE
        muli = dx**(i-1)
     END IF
     DO j=2,4
        IF (j == 2) THEN
           mulj = 1
        ELSE
           mulj = dy**(j-2)
        END IF
        mulj = mulj * (j-1)
        spval(2) = spval(2) + spl(i,j,ix,iy) * muli * mulj 
     END DO
  END DO

  RETURN
END SUBROUTINE eva2d_fd
!=====================================================
SUBROUTINE eva2d_sd(nx,ny,ix,iy,dx,dy,spl,spval)

! Evaluates the second derivatives of 2-dimensional cubic spline of function f(x,y)
!
! Input:  nx, ny              number of values in x and y  
!         ix, iy              pointer into the spline array spl
!         dx, dy              distance from x(ix) and y(iy)
!         spl                 array with spline data
! Output: spval(3)            evaluated function values
!                             spval(1) = d^2f/dx^2
!                             spval(2) = d^2f/(dxdy)
!                             spval(3) = d^2f/dy^2

  USE neo_precision
  
  IMPLICIT NONE
  
  INTEGER,                             INTENT(in)  :: nx, ny, ix, iy
  REAL(kind=dp),                       INTENT(in)  :: dx, dy
  REAL(kind=dp), DIMENSION(4,4,nx,ny), INTENT(in)  :: spl
  REAL(kind=dp), DIMENSION(3),         INTENT(out) :: spval

  INTEGER                                          :: i,j
  REAL(kind=dp)                                    :: muli, mulj

  spval = 0

! d^2f/dx^2
  DO i=3,4
     IF (i == 3) THEN
        muli = 1
     ELSE
        muli = dx**(i-3)
     END IF
     muli = muli * (i-1) * (i-2)
     DO j=1,4
        IF (j == 1) THEN
           mulj = 1
        ELSE
           mulj = dy**(j-1)
        END IF
        spval(1) = spval(1) + spl(i,j,ix,iy) * muli * mulj
     END DO
  END DO

! d^2f/(dxdy)
  DO i=2,4
     IF (i == 2) THEN
        muli = 1
     ELSE
        muli = dx**(i-2)
     END IF
     muli = muli * (i-1)
     DO j=2,4
        IF (j == 2) THEN
           mulj = 1
        ELSE
           mulj = dy**(j-2)
        END IF
        mulj = mulj * (j-1)
        spval(2) = spval(2) + spl(i,j,ix,iy) * muli * mulj
     END DO
  END DO

! d^2f/dy^2
  DO i=1,4
     IF (i == 1) THEN
        muli = 1
     ELSE
        muli = dx**(i-1)
     END IF
     DO j=3,4
        IF (j == 3) THEN
           mulj = 1
        ELSE
           mulj = dy**(j-3)
        END IF
        mulj = mulj * (j-1) * (j-2)
        spval(3) = spval(3) + spl(i,j,ix,iy) * muli * mulj
     END DO
  END DO

  RETURN
END SUBROUTINE eva2d_sd
!=====================================================
SUBROUTINE splreg(n,h,y,bi,ci,di)

! Makes a cubic spline of function y(x)
!
! Input:  n                   number of values in y  
!         h                   step size in x (aequidistant)
!         y(n)                y-values
! Output: bi(n),ci(n),di(n)   Spline parameters

  USE neo_precision

  IMPLICIT NONE

  INTEGER,                     INTENT(in)  :: n
  REAL(kind=dp),               INTENT(in)  :: h
  REAL(kind=dp), DIMENSION(n), INTENT(in)  :: y
  REAL(kind=dp), DIMENSION(n), INTENT(out) :: bi, ci, di

  REAL(kind=dp)                            :: ak1, ak2, am1, am2, c, e, c1
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: al, bt
  INTEGER                                  :: k, n2, i, i5

  ALLOCATE ( al(n), bt(n) )

  ak1 = 0
  ak2 = 0
  am1 = 0
  am2 = 0
  k = n-1
  al(1) = ak1
  bt(1) = am1
  n2 = n-2
  c = -4*h
  DO i = 1,n2
     e = -3*((y(i+2)-y(i+1))-(y(i+1)-y(i)))/h
     c1 = c-al(i)*h
     al(i+1) = h/c1
     bt(i+1) = (h*bt(i)+e)/c1
  END DO
  ci(n) = (am2+ak2*bt(k))/(1-al(k)*ak2)
  DO i = 1,k
     i5 = n-i
     ci(i5) = al(i5)*ci(i5+1)+bt(i5)
  END DO
  n2 = n-1
  DO i = 1,n2
     bi(i) = (y(i+1)-y(i))/h-h*(ci(i+1)+2*ci(i))/3
     di(i) = (ci(i+1)-ci(i))/h/3
  END DO
  DEALLOCATE ( al, bt )

  RETURN
END SUBROUTINE splreg
!====================================================
SUBROUTINE splper(n,h,y,bi,ci,di)

! Makes a cubic spline of periodic function y(x)
!
! Input:  n                   number of values in y
!         h                   step size in x (aequidistant)
!         y(n)                y-values
! Output: bi(n),ci(n),di(n)   Spline parameters

  USE neo_precision

  IMPLICIT NONE

  INTEGER,                     INTENT(in)  :: n
  REAL(kind=dp),               INTENT(in)  :: h
  REAL(kind=dp), DIMENSION(n), INTENT(in)  :: y
  REAL(kind=dp), DIMENSION(n), INTENT(out) :: bi, ci, di

  REAL(kind=dp)                            :: psi, ss
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: bmx, yl
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: amx1, amx2, amx3
  INTEGER                                  :: nmx, n1, n2, i, i1

  ALLOCATE ( bmx(n), yl(n), amx1(n), amx2(n), amx3(n) )
  
  bmx(1) = 1.e30_dp

  nmx=n-1
  n1=nmx-1
  n2=nmx-2
  psi=3.0_dp/h/h

  CALL spfper(n,amx1,amx2,amx3)

  bmx(nmx) = (y(nmx+1)-2*y(nmx)+y(nmx-1))*psi
  bmx(1)   =(y(2)-y(1)-y(nmx+1)+y(nmx))*psi
  DO i = 3,nmx
     bmx(i-1) = (y(i)-2*y(i-1)+y(i-2))*psi
  END DO
  yl(1) = bmx(1)/amx1(1)
  DO i = 2,n1
     i1 = i-1
     yl(i) = (bmx(i)-yl(i1)*amx2(i1))/amx1(i)
  END DO
  ss = 0
  DO i = 1,n1
     ss = ss+yl(i)*amx3(i)
  END DO
  yl(nmx) = (bmx(nmx)-ss)/amx1(nmx)
  bmx(nmx) = yl(nmx)/amx1(nmx)
  bmx(n1) = (yl(n1)-amx2(n1)*bmx(nmx))/amx1(n1)
  DO i = n2,1,-1
     bmx(i) = (yl(i)-amx3(i)*bmx(nmx)-amx2(i)*bmx(i+1))/amx1(i)
  END DO
  DO i = 1,nmx
     ci(i) = bmx(i)
  END DO

  DO i = 1,n1
     bi(i) = (y(i+1)-y(i))/h-h*(ci(i+1)+2*ci(i))/3
     di(i) = (ci(i+1)-ci(i))/h/3
  END DO
  bi(nmx) = (y(n)-y(n-1))/h-h*(ci(1)+2*ci(nmx))/3
  di(nmx) = (ci(1)-ci(nmx))/h/3
!
! Fix of problems at upper periodicity boundary
!
  bi(n) = bi(1)
  ci(n) = ci(1)
  di(n) = di(1)

  DEALLOCATE ( bmx, yl, amx1, amx2, amx3 )

  RETURN
END SUBROUTINE splper
!=====================================================
SUBROUTINE spfper(np1,amx1,amx2,amx3)

! Helper routine for splfi

  USE neo_precision

  IMPLICIT NONE

  INTEGER,                       INTENT(in)  :: np1
  REAL(kind=dp), DIMENSION(np1), INTENT(out) :: amx1, amx2, amx3
  REAL(kind=dp)                              :: beta, ss
  INTEGER                                    :: n, n1, i, i1

  n = np1-1

  n1 = n-1
  amx1(1) = 2
  amx2(1) = 0.5_dp
  amx3(1) = 0.5_dp
  amx1(2) = sqrt(15._dp)/2
  amx2(2) = ONE/amx1(2)
  amx3(2) = -.25_dp/amx1(2)
  beta = 3.75_dp
  DO i = 3,n1
     i1 = i-1
     beta = 4-ONE/beta
     amx1(i) = sqrt(beta)
     amx2(i) = ONE/amx1(i)
     amx3(i) = -amx3(i1)/amx1(i)/amx1(i1)
  END DO
  amx3(n1) = amx3(n1)+ONE/amx1(n1)
  amx2(n1) = amx3(n1)
  ss = 0
  DO i = 1,n1
     ss = ss+amx3(i)*amx3(i)
  END DO
  amx1(n) = sqrt(4-ss)

  RETURN
END SUBROUTINE spfper
!=====================================================
SUBROUTINE poi2d(hx,hy,mx,my,                         &
     xmin,xmax,ymin,ymax,                             &
     x,y,ix,iy,dx,dy,ierr)
! Creates Pointers for eva2d
!
! Input:  hx, hy              increment in x and y  
!         mx, my              standard (0) or periodic (1) spline
!         xmin, xmax          Minimum and maximum x
!         ymin, ymax          Minimum and maximum y
!         x, y                x and y values for spline avaluation
! Output: spval               evaluated function value
!         ix, iy              pointer into the spline array spl
!         dx, dy              distance from x(ix) and y(iy)
!         ierr                error (> 0)

  USE neo_precision

  IMPLICIT NONE

  REAL(kind=dp),                       INTENT(in)  :: hx, hy
  INTEGER,                             INTENT(in)  :: mx, my
  REAL(kind=dp),                       INTENT(in)  :: xmin, xmax, ymin, ymax
  REAL(kind=dp),                       INTENT(in)  :: x, y

  INTEGER,                             INTENT(out) :: ix, iy
  REAL(kind=dp),                       INTENT(out) :: dx, dy
  INTEGER,                             INTENT(out) :: ierr

  REAL(kind=dp)                                    :: dxx, x1, dyy, y1
  REAL(kind=dp)                                    :: dxmax, dymax

  ierr = 0

  dxx = x-xmin
  IF (mx .EQ. 0) THEN
     IF (dxx .LT. ZERO) THEN
        ierr = 1
        RETURN
     END IF
     IF (x .GT. xmax) THEN
        ierr = 2
        RETURN
     END IF
  ELSE
     dxmax = xmax - xmin
     IF(dxx .LT. ZERO) THEN
        dxx = dxx+(1+INT(abs(dxx/dxmax)))*dxmax
     ELSE IF(dxx .GT. dxmax) THEN
        dxx = dxx-(INT(abs(dxx/dxmax)))*dxmax
     END IF
  END IF
  x1 = dxx/hx
  ix = INT(x1)
  dx = hx*(x1-ix)
  ix = ix+1

  dyy = y-ymin
  IF (my .EQ. 0) THEN
     IF (dyy .LT. ZERO) THEN
        ierr = 3
        RETURN
     END IF
     IF (y .GT. ymax) THEN
        ierr = 4
        RETURN
     END IF
  ELSE
     dymax = ymax - ymin
     IF(dyy .LT. ZERO) THEN
        dyy = dyy+(1+INT(abs(dyy/dymax)))*dymax
     ELSE IF(dyy .GT. dymax) THEN
        dyy = dyy-(INT(abs(dyy/dymax)))*dymax
     END IF
  END IF
  y1 = dyy/hy
  iy = INT(y1)
  dy = hy*(y1-iy)
  iy = iy+1

  RETURN
END SUBROUTINE poi2d

  SUBROUTINE rhs_bo1(phi,y,dery)
! Right Hand Side of Differential Equation
!
!   y(1)                     - $\theta$
!   y(2)                     - $\int \rd \phi / B^2$
!   y(3)                     - $\int \rd \phi |\nabla \psi| / B^2$
!   y(4)                     - $\int \rd \phi K_G / B^3$
!   y(4+1)-y(4+npart)        - $\int \rd \phi \frac{\rd I_{fj}}{rd \phi}$ 
!   y(4+npart)-y(4+2*npart)  - $\int \rd \phi \frac{\rd H_{fj}}{rd \phi}$ 
!  
    USE neo_precision
    USE partpa_bo
    USE neo_exchange
!
    IMPLICIT NONE
!
    REAL(kind=dp),               INTENT(in)          ::  phi
    REAL(kind=dp)                                    ::  theta,bmod,gval,qval
    REAL(kind=dp)                                    ::  bmodm2,bmodm3
    REAL(kind=dp)                                    ::  geodcu,pardeb,bra,subsq,sq,sqeta
    REAL(kind=dp), DIMENSION(:),         INTENT(in)  ::  y
    REAL(kind=dp), DIMENSION(SIZE(y)), TARGET, INTENT(out) ::  dery  
    INTEGER                                          ::  ipass,i
    REAL(kind=dp),               POINTER             ::  p_iota, p_bm2
    REAL(kind=dp),               POINTER             ::  p_bm2gv, p_bm3ge
    REAL(kind=dp), DIMENSION(:), POINTER             ::  p_i, p_h
!
    p_iota  => dery(1)
    p_bm2   => dery(2)
    p_bm2gv => dery(3)
    p_bm3ge => dery(4)
    p_i     => dery(npq+1:npq+npart)
    p_h     => dery(npq+npart+1:npq+2*npart)

    theta=y(1)
! Evaluation of Splines
    CALL neo_eval(theta,phi,bmod,gval,geodcu,pardeb,qval)
!
    bmodm2 = ONE / bmod**2
    bmodm3 = bmodm2 / bmod
    bra = bmod / bmod0
!
    IF(pardeb*pard0.LE.ZERO .AND. pardeb.GT.ZERO) THEN
       ipass=1
    ELSE
       ipass=0
    ENDIF
    IF(ipmax.EQ.0.AND.pardeb*pard0.LE.ZERO.AND.pardeb.LT.ZERO) ipmax=1
    pard0=pardeb
!
    p_iota  = iota(psi_ind)
    p_bm2   = bmodm2
    p_bm2gv = bmodm2 * gval
    p_bm3ge = geodcu * bmodm3
!
    DO i=1,npart
       subsq = ONE - bra/eta(i)
       sqeta = sqrt(eta(i))
       IF(subsq.GT.ZERO) THEN
          isw(i) = 1
          icount(i) = icount(i) + 1
          ipa(i) = ipa(i)+ipass
          sq = sqrt(subsq)*bmodm2
          p_i(i) = sq 
          p_h(i) = sq*(4.0_dp/bra-ONE/eta(i))*geodcu/sqeta 
       ELSE
          sq = 0
          IF(isw(i).EQ.1) THEN
             isw(i) = 2
          ELSEIF(isw(i).EQ.2) THEN
             CONTINUE
          ELSE
             isw(i) = 0 
          ENDIF
          p_i(i) = 0
          p_h(i) = 0
       ENDIF
    ENDDO
!
    RETURN
  END SUBROUTINE rhs_bo1

SUBROUTINE rhs_cur1(phi,y,dery)
! Right Hand Side of Differential Equation
!
!   y(1)                     - $\theta$
!   y(2)                     - $\int \rd \phi / B^2$
!   y(3)                     - $\int \rd \phi |\nabla \psi| / B^2$
!   y(4)                     - $\int \rd \phi K_G / B^3$
!   y(4+1)-y(4+npart)        - $\int \rd \phi \frac{\rd I_{fj}}{rd \phi}$ 
!   y(4+npart)-y(4+2*npart)  - $\int \rd \phi \frac{\rd H_{fj}}{rd \phi}$ 
!  
  USE neo_precision
  USE partpa_cur
  USE neo_exchange
!
  IMPLICIT NONE
!
  REAL(kind=dp),               INTENT(in)          ::  phi
  REAL(kind=dp)                                    ::  theta,bmod,gval,qval
  REAL(kind=dp)                                    ::  bmodm2,bmodm3
  REAL(kind=dp)                                    ::  geodcu,pardeb,bra
  REAL(kind=dp)                                    ::  curfac
  REAL(kind=dp), DIMENSION(:)      , TARGET, INTENT(in)  ::  y
  REAL(kind=dp), DIMENSION(SIZE(y)), TARGET, INTENT(out) ::  dery  
  REAL(kind=dp),               POINTER             ::  pd_iota, pd_bm2
  REAL(kind=dp),               POINTER             ::  pd_bm2gv
  REAL(kind=dp),               POINTER             ::  pd_lamps
  REAL(kind=dp),               POINTER             ::  pd_lamb1n, pd_lamb1d
  REAL(kind=dp),               POINTER             ::  pd_lamb2n, pd_lamb2d
  REAL(kind=dp), DIMENSION(:), POINTER             ::  pd_l, pd_k1, pd_k
  REAL(kind=dp), DIMENSION(:), POINTER             ::  p_l, p_k1, p_k
  INTEGER      , PARAMETER                         ::  k_meth = 2
!
  pd_iota   => dery(1)
  pd_bm2    => dery(2)
  pd_bm2gv  => dery(3)
  pd_lamps  => dery(4)
  pd_lamb1n => dery(5)
  pd_lamb1d => dery(6)
  pd_lamb2n => dery(7)
  pd_lamb2d => dery(8)
  pd_l      => dery(npq_cur+1:npq_cur+npart_cur)
  pd_k1     => dery(npq_cur+npart_cur+1:npq_cur+2*npart_cur)
  pd_k      => dery(npq_cur+2*npart_cur+1:npq_cur+3*npart_cur)

  p_l       => y(npq_cur+1:npq_cur+npart_cur)
  p_k1      => y(npq_cur+npart_cur+1:npq_cur+2*npart_cur)
  p_k       => y(npq_cur+2*npart_cur+1:npq_cur+3*npart_cur)

  theta=y(1)
! Evaluation of Splines
  CALL neo_eval(theta,phi,bmod,gval,geodcu,pardeb,qval)
!
  bra = bmod / bmod0
  bmodm2 = ONE / bra**2
  bmodm3 = bmodm2 / bra
  curfac = geodcu * fac / bmod0

!!$  PRINT *, 'phi ',phi
!!$  PRINT *, 'theta ',theta
!!$  PRINT *, 'bra ',bra
!!$  PRINT *, 'bmodm2 ',bmodm2
!!$  PRINT *, 'bmodm3 ',bmodm3
!!$  PRINT *, 'curfac ',curfac
!!$  PRINT *, 'geodcu ',geodcu
!!$  PRINT *, 'fac ',fac
!!$  PRINT *, 'bmod0 ',bmod0
!!$  PRINT *, 'bmod ',bmod
!!$  PRINT *, 'gvql ',gval
!!$  PRINT *, 'qvql ',qval

  yfac   = ONE - y_part*bra
  sqyfac = SQRT(yfac)
!
  pd_iota   = iota(psi_ind)
  pd_bm2    = bmodm2
  pd_bm2gv  = bmodm2 * gval
  pd_lamps  = curfac * bmodm3
  pd_lamb1n = qval * bmodm2 * pd_lamps
  pd_lamb1d = qval * bmodm2
  pd_lamb2n = pd_lamps
  pd_lamb2d = 1

  yfac      = ONE - y_part*bra
  sqyfac    = SQRT(yfac)
  pd_l      = bmodm2 * sqyfac
  IF (k_meth .EQ. 1) THEN
     pd_k1     = curfac / yfac / sqyfac
     pd_k      = pd_l * p_k1
  ELSE
     k_fac1 = -(2.0_dp/bra + y_part / 2.0_dp / yfac  ) * pardeb / bmod0
     k_fac2 = bmodm2 / yfac * curfac
     pd_k1  = k_fac1 * p_k1 + k_fac2
     pd_k   = p_k1
  END IF
!
  RETURN
END SUBROUTINE rhs_cur1


SUBROUTINE rk4d_bo1(x,y,h)
! Differential Equation Solver
  USE neo_precision
  USE sizey_bo
  USE neo_rhsbo
!
  IMPLICIT NONE
!
  REAL(kind=dp),               INTENT(inout) ::  x
  REAL(kind=dp),               INTENT(in)    ::  h
  REAL(kind=dp)                              ::  hh,h6,xh
  REAL(kind=dp), DIMENSION(:), INTENT(inout) ::  y
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE   ::  dydx,yt,dyt,dym
  INTEGER                                    ::  i
!
  ALLOCATE(dydx(ndim))
  ALLOCATE(yt(ndim))
  ALLOCATE(dyt(ndim))
  ALLOCATE(dym(ndim))
!
  hh = h / 2
  h6 = h/6
  xh = x+hh
  CALL rhs_bo(x,y,dydx)
  DO i=1,ndim
     yt(i)=y(i)+hh*dydx(i)
  END DO
  CALL rhs_bo(xh,yt,dyt)
  DO i=1,ndim
     yt(i)=y(i)+hh*dyt(i)
  END DO
  CALL rhs_bo(xh,yt,dym)
  DO i=1,ndim
     yt(i)=y(i)+h*dym(i)
     dym(i)=dyt(i)+dym(i)
  END DO
  x=x+h
  CALL rhs_bo(x,yt,dyt)
  DO i=1,ndim
     y(i)=y(i)+h6*(dydx(i)+dyt(i)+2*dym(i))
  END DO

  DEALLOCATE (dydx, yt, dyt, dym)

  RETURN
END SUBROUTINE rk4d_bo1


SUBROUTINE rk4d_cur1(x,y,h)
! Differential Equation Solver
  USE neo_precision
  USE sizey_cur
  USE neo_rhscur
!
  IMPLICIT NONE
!
  REAL(kind=dp),               INTENT(inout) ::  x
  REAL(kind=dp),               INTENT(in)    ::  h
  REAL(kind=dp)                              ::  hh,h6,xh
  REAL(kind=dp), DIMENSION(:), INTENT(inout) ::  y
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE   ::  dydx,yt,dyt,dym
  INTEGER                                    ::  i
!
  ALLOCATE(dydx(ndim_cur))
  ALLOCATE(yt(ndim_cur))
  ALLOCATE(dyt(ndim_cur))
  ALLOCATE(dym(ndim_cur))
!
  hh = h/2
  h6 = h/6
  xh = x+hh
  CALL rhs_cur(x,y,dydx)
  DO i=1,ndim_cur
     yt(i)=y(i)+hh*dydx(i)
  END DO
  CALL rhs_cur(xh,yt,dyt)
  DO i=1,ndim_cur
     yt(i)=y(i)+hh*dyt(i)
  END DO
  CALL rhs_cur(xh,yt,dym)
  DO i=1,ndim_cur
     yt(i)=y(i)+h*dym(i)
     dym(i)=dyt(i)+dym(i)
  END DO
  x=x+h
  CALL rhs_cur(x,yt,dyt)
  DO i=1,ndim_cur
     y(i)=y(i)+h6*(dydx(i)+dyt(i)+2*dym(i))
  END DO

  DEALLOCATE (dydx, yt, dyt, dym)

  RETURN
END SUBROUTINE rk4d_cur1
