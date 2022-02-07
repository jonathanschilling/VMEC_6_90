      module real_ptr_type
      use kind_spec
      
      type real_ptr
            real(rprec), dimension(:), pointer :: x
      end type real_ptr

      end module real_ptr_type

      module optim
      use optim_params
      implicit none
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: nxc = 1000
      integer, parameter :: iout = 54
      real(rprec) :: bigno = 1.e10_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ierr_vmec, ns_booz_max, ns_surf_max, iunit_dkes
     1 , ns_ball_max, ns_neo_max, unit_outdata                          !!VMECCOBRA, NEO 
      integer, dimension(nxc) :: nbrho_opt, mbrho_opt
      integer, dimension(2*ntor1d) :: nrz0_opt
      integer :: irho_bdy, irm0_bdy, izm0_bdy,
     1   nrad, num_ai
      integer :: nfp_opt, mpol_opt, mpol1_opt, ntor_opt, 
     1      ntor1_opt, mnmax_opt, iunit_opt, iunit_opt_local,
     2      nextcur_opt, nextcur_vmec, iter_min, min_count
      real(rprec), dimension(-ntord:ntord,0:mpol1d) :: rhobc
      real(rprec), dimension(-ntord:ntord,-mpol1d:mpol1d) :: delta_mn
      real(rprec) :: wp_opt, wb_opt, rmax_opt, rmin_opt, zmax_opt, 
     1      aminor_opt
      real(rprec) :: aspect_opt, coil_complex_opt, rbtor_opt
      real(rprec) :: chisq_min
      real(rprec), dimension(nsd) :: vp_opt, iota_opt, jcurv_opt,
     1     phip_opt, buco_opt, Dmerc_opt, jdotb_opt
      real(rprec), dimension(nsd) ::  pres_opt                           !!COBRA
      real(rprec) :: version_opt, am0_9, am10                            !!COBRA
      real(rprec) :: tm0_9, tm10                                         !!LEGENDRE
      real(rprec) :: raxis_old(0:ntord), zaxis_old(0:ntord)

      integer, allocatable :: ns_booz(:), ns_surf(:), ns_neo(:),        
     1  ns_ball(:)                                                       !!VMECCOBRA, NEO

      logical, allocatable, dimension(:) :: lneed_modB
      logical :: lone_step, lscale_only, lrestart, lniter1, lajax

      character*120 :: input_file, home_dir, min_ext
      character*130 :: output_file, min_input_file, min_wout_file
c-----------------------------------------------
      end module optim

      module boozer_params
      use kind_spec
      integer :: mboz, nboz, nu_boz, nv_boz, nunv, mnboz, nu2_b
      real(rprec), dimension(:), allocatable :: 
     1  rmnc_bdy, zmns_bdy, xm_bdy, xn_bdy
      real(rprec), dimension(:,:), allocatable :: 
     1  rmnc_opt, zmns_opt,lmns_opt
      end module boozer_params
      
      module chisq_mod
      use kind_spec
      integer, parameter :: ntargets=65
      integer, parameter :: unit_dkes=29, unit_cobra=28, 
     1                      unit_outdata=29, unit_diagno=29
      integer, parameter :: ivar_aspect=1, ivar_maxj=2, ivar_beta=3,
     1  ivar_pedge=4,    ivar_pgrad=5,   ivar_iota=6,   ivar_iota_p=7,  
     2  ivar_bmn=8,      ivar_bmin=9,    ivar_bmax=10,  ivar_jstar=11,
     3  ivar_jinvar=12,  ivar_ripple=13, ivar_dkes=14,  ivar_neo=15, 
     4  ivar_pseudo=16,  ivar_orbit=17,  ivar_well=18,  ivar_mercier=19,
     5  ivar_balloon=20, ivar_kink=21,   ivar_jac=22,   ivar_dsubr=23,
     6  ivar_bootsj=24,  ivar_fluxp=25,  ivar_rbtor=26, ivar_curve=27,
     7  ivar_center=28,  ivar_ellipticity=29,           ivar_kappa=30,
     8  ivar_zmax=31,    ivar_vv=32,     ivar_vv_rms=33, 
     9  ivar_coil_complexity=34,         ivar_coil_jmax=35, 
     A  ivar_berr_ave=36,  ivar_coil_len=37,   ivar_coil_sep=38,
     B  ivar_coil_curv=39, ivar_coil_acc=40,   ivar_coil_icurv=41,
     C  ivar_coil_ymin=42, ivar_coil_vfreg=43, ivar_extcur=44, 
     D  ivar_oh=45, ivar_vac_island=46, ivar_coil_polcur=47,
     E  ivar_coil_pldist=48, ivar_coil_plcross=49, ivar_coil_reg=50,
     F  ivar_coil_curdens=51, ivar_iota_bounds=52, ivar_coil_rmax=53,
     G  ivar_coil_vfrmax=54, ivar_vv_max=55, ivar_coil_aBerr=56,
     H  ivar_coil_mBerr=57, ivar_Jconf=58, ivar_p_prof=59, 
     I  ivar_emis=60, ivar_pdamp=61, ivar_emisdamp=62, ivar_diagno=63,
     J  ivar_curtor=64, ivar_jedge=65

      integer, dimension(:), allocatable :: index_array 
      real(rprec), dimension(:), allocatable :: wegt,
     1   chisq_match, chisq_target
      character*2 :: comm1
      character*(18), dimension(:), allocatable :: chisq_descript

      character*(*), dimension(ntargets), parameter :: descript =
     1   ( /
     1   'Aspect Ratio      ','Max. Current      ','<Beta>            ',
     2   'Edge Pressure     ','Pressure Gradient ',     
     3   'Iota (1/q)        ','d-Iota/ds         ',
     4   'Bmn Spectrum      ','Bmin Spectrum     ','Bmax Spectrum     ',
     5   'J*  Spectrum      ','J-Invariant       ','<Magnetic Ripple> ',
     6   'DKES L11 Transport','Ripple Diffusion  ','Pseudo Symmetry   ',     
     7   'Orbit Confinement ','Magnetic Well     ','Mercier           ',
     8   'Ballooning Stab.  ','External Kink     ','Resonant Jacobian ',
     9   'Resistive D_R     ',
     A   'Bootstrap Current ','Poloidal Flux     ','R-Btor            ',     
     B   'Surface Curvature ','Radial Centering  ','Ellipticity       ',
     C   '<Elongation>      ','Zmax Boundary     ', 
     D   '3D Boundary       ','3D Shape-RMS      ',
     E   'Nescoil Complexity','Nescoil Sheet Jmax','Nescoil Berr-ave  ',
     F   'Coil Length       ','Coil Separation   ','Coil Curvature    ',
     G   'Coil-Access Reqts ','Coil Integ. Curv. ','Coil ymin         ',
     H   'Coil VF Regulariz.','Coil Currents     ','OH Constraint     ',
     I   'Vac Island Width  ','Coil I-Poloidal   ','Coil-Plasma Dist. ',
     J   'Coil-Pl X CoilCoil','Coil MC Regulariz.','Coil Current Dens.',
     K   'Iota Bounds       ','Coil MC max radius','Coil VF max radius',
     L   '3D Shape-Max dev  ','Coil Bnorm avg err','Coil Bnorm max err',
     M   'J-contour confine.','Pressure Profile  ','SXR Emissivity    ',
     N   'Press Prof Damping','Emiss Prof Damping','DIAGNO mag. diag. ',
     O   'Total tor. current','Edge current dens.'
     Z   / )
 
      end module chisq_mod

      module cmnf1
      use kind_spec
      integer :: nmn1
      real(rprec), allocatable, dimension(:) ::  xm1, xn1, rmn1, zmn1
      real(rprec) :: q1
      end module cmnf1

      module cmnf2
      use kind_spec
      integer :: nmn2
      real(rprec), allocatable, dimension(:) ::  xm2, xn2, rmn2, zmn2
      real(rprec) :: q2
      end module cmnf2
      
      module geom
      use kind_spec
      real(rprec) ::  xp, yp, zp, pi2, alp
      end module geom

      module newtv
      use kind_spec
      integer, parameter :: np=40
      integer :: nn
      real(rprec) :: fvec(np)
      end module newtv

      module legendre_params                                                 !! LEGENDRE
      use kind_spec                                                          !! LEGENDRE
      integer :: n_leg                                                       !! LEGENDRE
      real(rprec), allocatable, dimension(:,:) :: a_leg, b_leg,              !! LEGENDRE
     1   a_leg_inv, b_leg_inv                                                !! LEGENDRE
      real(rprec), allocatable, dimension(:) :: tc, ti, tm                   !! LEGENDRE
      end module legendre_params                                             !! LEGENDRE
      
      module mpi_params                                                      ! MPI
      integer :: myid=0, numprocs=0, ierr_mpi                                  ! MPI
      integer, parameter :: master=0                                         ! MPI
      end module mpi_params                                                  ! MPI
