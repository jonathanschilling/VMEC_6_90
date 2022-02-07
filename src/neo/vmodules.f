MODULE nrtype
! Definition of types taken from Numerical Recipes
  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(12,100) 
  REAL(DP), PARAMETER :: PI=3.141592653589793238462643383279502884197_dp
  REAL(DP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_dp
  REAL(DP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_dp
  REAL(DP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_dp
  REAL(DP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_dp
  REAL(DP), PARAMETER :: ONE = 1,  ZERO = 0
END MODULE nrtype

MODULE neo_precision
  USE nrtype
END MODULE neo_precision

MODULE neo_parameters
  USE neo_precision
  USE nrtype
END MODULE neo_parameters

MODULE neo_input
! Input from data files (Boozer)
  USE neo_precision
  INTEGER, DIMENSION(:),     ALLOCATABLE :: ixm, ixn
  INTEGER, DIMENSION(:),     ALLOCATABLE :: pixm, pixn
  INTEGER, DIMENSION(:),     ALLOCATABLE :: i_m, i_n

  REAL(kind=dp),    DIMENSION(:),     ALLOCATABLE :: es
  REAL(kind=dp),    DIMENSION(:),     ALLOCATABLE :: pprime
  REAL(kind=dp),    DIMENSION(:),     ALLOCATABLE :: sqrtg00

  REAL(kind=dp),    DIMENSION(:,:),   ALLOCATABLE :: rmnc,  zmnc,  lmnc
  REAL(kind=dp),    DIMENSION(:,:),   ALLOCATABLE :: bmnc

  REAL(kind=dp) :: flux

  INTEGER  :: m0b, n0b
  INTEGER  :: ns, mnmax, nfp
  INTEGER  :: m_max, n_max
END MODULE neo_input

MODULE neo_work
! Working parameters
  USE neo_precision

  REAL(kind=dp) ::   theta_start
  REAL(kind=dp) ::   theta_end
  REAL(kind=dp) ::   theta_int

  REAL(kind=dp) ::   phi_start
  REAL(kind=dp) ::   phi_end
  REAL(kind=dp) ::   phi_int

  REAL(kind=dp),    DIMENSION(:,:),     ALLOCATABLE :: cosmth,sinmth
  REAL(kind=dp),    DIMENSION(:,:),     ALLOCATABLE :: cosnph,sinnph
  REAL(kind=dp),    DIMENSION(:),       ALLOCATABLE :: theta_arr,phi_arr
  REAL(kind=dp),    DIMENSION(:,:),     ALLOCATABLE :: r,z,l,b
  REAL(kind=dp),    DIMENSION(:,:),     ALLOCATABLE :: r_tb,z_tb,p_tb,b_tb
  REAL(kind=dp),    DIMENSION(:,:),     ALLOCATABLE :: r_pb,z_pb,p_pb,b_pb
  REAL(kind=dp),    DIMENSION(:,:),     ALLOCATABLE :: gtbtb,gpbpb,gtbpb
  REAL(kind=dp),    DIMENSION(:,:),     ALLOCATABLE :: isqrg,sqrg11,kg,pard
  REAL(kind=dp),    DIMENSION(:,:),     ALLOCATABLE :: bqtphi
END MODULE neo_work

MODULE neo_exchange
! Working parameters
  USE neo_precision

  REAL(kind=dp)              :: b_min, b_max
  INTEGER                    :: nper
  REAL(kind=dp)              :: rt0, rt0_g 
  REAL(kind=dp)              :: bmref, bmref_g
  INTEGER                    :: nstep_per, nstep_min, nstep_max
  INTEGER                    :: write_integrate
  INTEGER                    :: write_diagnostic
  INTEGER                    :: write_cur_inte
  REAL(kind=dp)              :: acc_req
  INTEGER                    :: no_bins
  INTEGER                    :: psi_ind
  INTEGER                    :: calc_nstep_max
  REAL(kind=dp)              :: theta_bmin, phi_bmin
  REAL(kind=dp)              :: theta_bmax, phi_bmax
  REAL(kind=dp),    DIMENSION(:),     ALLOCATABLE :: iota
  REAL(kind=dp),    DIMENSION(:),     ALLOCATABLE :: curr_pol
  REAL(kind=dp),    DIMENSION(:),     ALLOCATABLE :: curr_tor
  REAL(kind=dp)              :: fac
  INTEGER                    :: calc_cur
  INTEGER                    :: hit_rat, nfp_rat, nfl_rat
  REAL(kind=dp)              :: delta_theta_rat 
END MODULE neo_exchange

MODULE neo_output
  USE neo_precision
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE ::  epspar
  REAL(kind=dp)                            ::  epstot,ctrone,ctrtot
  REAL(kind=dp)                            ::  bareph,barept,drdpsi
  REAL(kind=dp)                            ::  yps
  INTEGER                                  ::  nintfp
  INTEGER                                  ::  ierr  
  REAL(kind=dp)                            ::  lambda_b
  REAL(kind=dp)                            ::  lambda_b1, lambda_b2
  REAL(kind=dp)                            ::  lambda_ps1, lambda_ps2
  REAL(kind=dp)                            ::  avnabpsi
END MODULE neo_output

MODULE neo_spline
! Spline arrays
  USE neo_precision
  INTEGER, PARAMETER       ::   mt = 1
  INTEGER, PARAMETER       ::   mp = 1
  INTEGER                  ::   theta_ind, phi_ind
  INTEGER                  ::   ierr

  REAL(kind=dp) ::  theta_d, phi_d 

! Spline array for modb
  REAL(kind=dp),    DIMENSION(:,:,:,:), ALLOCATABLE :: b_spl 
! Spline array for geodesic curviture
  REAL(kind=dp),    DIMENSION(:,:,:,:), ALLOCATABLE :: k_spl
! Spline array for sqrg11
  REAL(kind=dp),    DIMENSION(:,:,:,:), ALLOCATABLE :: g_spl
! Spline array for parallel derivative
  REAL(kind=dp),    DIMENSION(:,:,:,:), ALLOCATABLE :: p_spl
! Spline array for quasi-toroidal phi component of b
  REAL(kind=dp),    DIMENSION(:,:,:,:), ALLOCATABLE :: q_spl
END MODULE neo_spline

MODULE neo_control
! Control parameters from inpt file
  USE neo_precision
  CHARACTER(120)                     :: in_file                             !SPH
  CHARACTER(120)                     :: out_file                            !SPH
  CHARACTER(120)                      :: cur_file
  INTEGER                            :: theta_n
  INTEGER                            :: phi_n
  INTEGER                            :: s_ind_in
  INTEGER                            :: write_progress
  INTEGER                            :: write_output_files
  INTEGER                            :: spline_test
  INTEGER                            :: max_m_mode, max_n_mode
  INTEGER                            :: lab_swi, inp_swi, ref_swi, eout_swi
  INTEGER                            :: no_fluxs
  INTEGER, DIMENSION(:), ALLOCATABLE :: fluxs_arr
END MODULE neo_control

MODULE neo_units
! Units and Formats
  USE neo_precision
  INTEGER, PARAMETER ::   r_u2  = 4
  INTEGER, PARAMETER ::   r_us  = 5
  INTEGER, PARAMETER ::   r_u23 = 23
  INTEGER, PARAMETER ::   r_ua  = 21
  INTEGER, PARAMETER ::   w_u4  = 10
  INTEGER, PARAMETER ::   r_u10 = 3
  INTEGER, PARAMETER ::   w_us  = 6
  INTEGER, PARAMETER ::   w_u10 = 7
  INTEGER, PARAMETER ::   w_u20 = 8
  INTEGER, PARAMETER ::   w_u30 = 9
  INTEGER, PARAMETER ::   w_u50 = 11
  INTEGER, PARAMETER ::   w_u60 = 12
  INTEGER, PARAMETER ::   w_u70 = 13
  INTEGER, PARAMETER ::   w_u80 = 14
  INTEGER, PARAMETER ::   w_u90 = 15

  INTEGER            ::   w_u6_open, w_u1 = w_u10, w_u2 = w_u20, w_u3 = w_u30,      &
                          w_u5 = w_u50,w_u6 = w_u60, w_u7 = w_u70, w_u8 = w_u80,    &
                          w_u9 = w_u90, r_u1 = r_u10

  CHARACTER(20),PARAMETER :: format220="(500d18.5)"

  CHARACTER*(120) :: arg1, extension    ! LPK
  INTEGER  :: numargs                   ! LPK

END MODULE neo_units

MODULE sizey_bo
  USE neo_precision
! Definition for rk4d_bo also used in main routine neo
  INTEGER            ::  npart
  INTEGER            ::  multra
  INTEGER            ::  ndim
  INTEGER, PARAMETER ::  npq = 4
END MODULE sizey_bo

MODULE partpa_bo
  USE neo_precision
! Exchange between flint_bo and rhs_bo
  USE sizey_bo
  INTEGER                                  :: ipmax
  INTEGER,       DIMENSION(:), ALLOCATABLE :: isw,ipa,icount
  REAL(kind=dp)                            :: pard0,bmod0
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: eta
END MODULE partpa_bo

MODULE sizey_cur
  USE neo_precision
! Definition for rk4d_bo also used in main routine neo
  INTEGER            ::  npart_cur
  INTEGER            ::  ndim_cur
  INTEGER, PARAMETER ::  npq_cur = 8
  INTEGER            ::  alpha_cur
END MODULE sizey_cur

MODULE partpa_cur
  USE neo_precision
! Exchange between flint_cur and rhs_cur
  USE sizey_cur
  REAL(kind=dp)                            :: bmod0
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: y_part, yfac, sqyfac
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: k_fac1, k_fac2
END MODULE partpa_cur

MODULE neo_rk4dbo
  INTERFACE rk4d_bo 
     SUBROUTINE rk4d_bo1(x,y,h)
       USE neo_precision
       REAL(dp),               INTENT(inout) ::  x
       REAL(dp),               INTENT(in)    ::  h
       REAL(dp), DIMENSION(:), INTENT(inout) ::  y
     END SUBROUTINE rk4d_bo1
  END INTERFACE
END MODULE neo_rk4dbo

MODULE neo_rk4dcur
  INTERFACE rk4d_cur 
     SUBROUTINE rk4d_cur1(x,y,h)
       USE neo_precision
       REAL(dp),               INTENT(inout) ::  x
       REAL(dp),               INTENT(in)    ::  h
       REAL(dp), DIMENSION(:), INTENT(inout) ::  y
     END SUBROUTINE rk4d_cur1
  END INTERFACE
END MODULE neo_rk4dcur

MODULE neo_rhsbo
  INTERFACE rhs_bo 
     SUBROUTINE rhs_bo1(phi,y,dery)
       USE neo_precision
       REAL(dp),                             INTENT(in)  ::  phi
       REAL(dp), DIMENSION(:),               INTENT(in)  ::  y
       REAL(dp), DIMENSION(SIZE(y)), TARGET, INTENT(out) ::  dery
     END SUBROUTINE rhs_bo1
  END INTERFACE
END MODULE neo_rhsbo

MODULE neo_rhscur
  INTERFACE rhs_cur 
     SUBROUTINE rhs_cur1(phi,y,dery)
       USE neo_precision
       REAL(dp),                             INTENT(in)  ::  phi
       REAL(dp), DIMENSION(:), TARGET,       INTENT(in)  ::  y
       REAL(dp), DIMENSION(SIZE(y)), TARGET, INTENT(out) ::  dery
     END SUBROUTINE rhs_cur1
  END INTERFACE
END MODULE neo_rhscur
