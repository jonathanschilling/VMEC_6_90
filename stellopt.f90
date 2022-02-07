#!/bin/sh
#10/02:  Added new targets for experimental pressure profile (CHISQ_PPROF) and
#        SXR chord integrals  (CHISQ_SXR). Also, increased CHISQ_JCONF speed.
#8/02:   Added new target for confining J contours (CHISQ_JCONF)
#7/02:   Added new Piecewise limiter constraints (PHI_LIM, R_LIM, Z_LIM), and
#        targets for plasma energy (TARGET_EPLASMA & sigma)
#6&7/02: more improvements to Levenberg alg.
#early/02: improvements to Levenberg by adding Stepopt routines
#1/02:   Added new targets for max, min values of iota (TARGET_IOTA_MAX, TARGET_IOTA_MIN) and
#        corresponding SIGMAs. Added target for finite ballooning eigenvalue (TARGET_BALLOON array).
#        Added COIL_GEOM constraints on modular coil current regularization and current density
#        to subroutine CHISQ_COILGEOM.
#        Added LKEEP_MINS logical, to allow user to restore history of input.min files generated during
#        optimization run.
#        Modified subroutine READ_COILGEOM_ELEM to include COIL_GEOM weights properly.
#8/01:   Replaced R00_OPT with R00_SCALE in OPTIMUM namelist. Now, R00_SCALE is the fraction by
#        which the initial equilibrium geometry is to be scaled. Thus, R00_OPT = R00_SCALE*RBC(0,0),
#        where RBC(0,0) is the initial equilibrium m=0, n=0 component.
#        Added "data_harvest" subroutine to gather output summaries from various files and combine
#        the results in a "data_summary" file
#
#        General "help" can now be flashed to the screen by typing 
#                 xstellopt -h   OR   xstellopt
#
#        Merger of ORNL and PPPL optimizer versions:
#        The Levenburg-Marquardt algorithm has been modified to follow the lowest chi-square value found, 
#        even when it is generated as part of the Jacobian evaluation. This speeds convergence in many cases.
#
#        Two new optimization algorithms are now available: Genetic (Long-Poe) 
#        and Differential-Evolution (MCZ).  The choice of algorithm is controlled by a new switch
#        NOPT_ALG = 0    Levenburg-Marquandt  (LM)       (default)
#                 = 1    Genetic  (GA)
#                 = 2    Differential Evolution  (DE)
#        Additional control parameters for GA and DE are specified in a new name list  (GA_DE) 
#        appended at the end of the input file control file.  See the comments in 
#        ga_opt.f90 and  de_opt.f90  for the meaning of the namelist controls
#
#        Both the GA and DE algorithms have a resume capability and write a "resume" file at 
#        the end of each generation. In order to restart an existing run, invoke xstellopt with the "-r" 
#        switch from within the stellopt_* subdirectory, e.g.  
#                 xstellopt -r  mycase
#        The resume files can also be used to preload an initial (sub-)population, particularly for DE.
#
#        Targetting of fast-ion orbit confinement, via Monte-Carlo simulation has been added (LPK). To implement,
#        set  LORBIT_OPT = .T, SIGMA_ORBIT, and see CHISQ_ORBIT and refer to the comments at the begining of LOAD_TARGET
#
#        Targetting of D_R (resistive stability), for improvement of flux-surface quality, has been added (LPK)
#        set  LDSUBR_OPT = .T, SIGMA_DSUBR, and see CHISQ_DSUBR and the comments at the beginning of LOAD_TARGET
#
#        An option to allow use of P. Garabedian's delta representation as the internal optimization representation 
#        (for fixed boundary optimization) was added (LPK), controlled by a new switch:
#            NOPT_BOUNDARY = 0 Hirshman-Breslau rep.  (default)
#                          = 1  Garabedian representation
#
#        Targetting of resonant components of the Boozer Jacobian, in order to
#        encourage good flux-surface quality.  Set N_JAC and M_JAC arrays for
#        desired resonances, and SIGMA_JAC for penalty.  See CHISQ_JAC.
#
#        Direct targetting of the effective ripple via the NEO program (LPK).
#        see SIGMA_NEO and CHISQ_NEO.
#
#        Regularization of coil currents.  See CHISQ_EXTCUR and CHISQ_OH.
#
#7/01:   Added generate_mgrid to generate the mgrid file from the coilsindata namelist.
#        Added chisq_coilgeom to compute targets based on coils       
#5/01:   Added LCOIL_GEOM for coupling to coilgeom code (varies coil geometry).
#2/01:   Modified chisq target for Bmn to be proportional to |B|**2 (modes we do NOT want)/
#        |B|**2 (modes that are or the correct helicity). 
#1/01:   Merged MPI version into Forked multiprocessed code. Still require separate
#        lev_opt.f90 libraries.
#8/00:   VMEC version of COBRA replaced older Boozer-coordinate version. Legendre representation
#        for the iota, current and pressure profile included. New format for in_cobra implemented
#        to facilitate dynamic memory allocation in xcobravmec code.
#8/00:   Fixed logic associated with rescaling the equilibrium (R00_opt, B00_scale) to be
#        consistent with free-boundary restarts.
#6/00:   MCZ & RH  added penalty on separation to a vacuum vessel surface
#4/00:   Added LEXTCUR array to optimizer namelist, to control which external currents will be 
#        varied (=T) and which currents will be fixed (=F) during a free-boundary optimization.
#3/00:   Merged ORNL/PPPL optimizer load_target functions. Added logicals laspect_max (=T, then
#        chisq_aspect = 0 if A < Target) and lbeta_min (=T, then chisq_beta = 0 if beta > Target). 
#        Both are defaulted to be F.
#2/00:   Added free-boundary optimization option. This allows user to vary coil
#        current (extcur array) rather than boundary coefficients, and is active
#        only if LFREEB = T.
#4/99:   Added R00_opt, B00_scale to control scaling of R so that R00 = R00_opt, and
#        scaling of magnetic fields (through phiedge) by the factor B00_scale. Physics parameters, 
#        such as beta and iota, are held fixed, so physics chi-sq should be invariant during scaling.
#9/98:   Separated optimization code from VMEC source
#        Added SYSTEM subroutine (for CRAY, others have it already) as a prerequisite
#        for spawning processes from within VMEC (such as executing nescoil, bnorm CODES)
#        Added specific system calls to bnorm and nescoil from LOAD_TARGET routine
#09/98:  Added optimizer namelist parameter COIL_SEPARATION (the winding-surface-to-plasma 
#        boundary separation, in [M] or whatever units R,Z are in) to OPTIMUM namelist. 
#        Eventually, this may be an array.
#        Modified READ_WOUT_OPT to write out COIL_SEPARATION to be read by BNORM code
#        (besonderes dankbar am P. Merkel for allowing us to merge (parts) of this 
#        outstanding code with VMEC!) Benchmarked successfully (9/12/98)
#08/98:  Added LCOIL_COMPLEX logical to optimizer namelist. If = T, it allows computation of a 
#        coil-complexity contribution to chi-sq in the optimizer. Also added SIGMA_COIL_COMPLEX 
#        to weight the complexity measure. REPLACED WITH LCOIL_OPT (11/00) AND LNESCOIL_OPT (1/02).
#03/97:  Included a flux surface curvature limiting term in Chi-Sq.
#        This prevents cusp-like regions from forming in the 
#        flux surfaces.  Subroutine flux_surf_curv() was added.  This
#        calculates the local curvature around outer flux surfaces at four
#        toroidal planes (phi = 0, 90, 180, 270) using finite differencing
#        methods.  Statistical measures of the variation of this curvature
#        from mean values are then used in subroutine LOAD_TARGET to prevent
#        the large variations of curvature which would be associated with
#        cusp-like regions.
#        Added a subroutine WRITE_INDATA to make periodic dumps of the
#        input, boundary shape and optimization parameters so that
#        restarts can be made from these dump files.  WRITE_INDATA
#        writes out a formatted restart file instead of the earlier
#        namelist write.  This makes the restart file easier to read
#        and allows restarting with a different number of modes
#        (i.e., using the mpold, ntord parameters in vmec0.inc).         
#02/97:  Fixed up optimization routine LOAD_PARAMS to check IER return
#        error flag from VMEC. Fixed VMEC routine EQSOLVE to check for
#        jacobian resets and to restart with increased/decreased DELT
#        if iterations are stuck. Added NITER_OPT to INPUT file 
#        (namelist $OPTIMUM) to control maximum number of optimization
#        iterations. Also added LRESET_OPT to control reset style of
#        vmec during optimization run.
#---------------------------------------------------------
cat > temp.c << "EOC"
cat > vmodules.f << "EOF"
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
EOF
cat > voptimize.f << "EOF"
!********************************************************************
!                              D * I * S * C * L * A * I * M * E * R
!
!       You are using a BETA version of the program STELLOPT, which is currently
!       under development by S. P. Hirshman at the Fusion Energy Division,
!       Oak Ridge National Laboratory.  Please report any problems or comments
!       to him.  As a BETA version, this program is subject to change
!       and improvement without notice.
!
!********************************************************************
      program stellarator_optimizer
!
!     Driver routine (front end)
!
!     Invoke optimizer as follows:
!
!     xstellopt input_file_name  (input_file_name == name of vmec input_file)
!
      use optim, only: lone_step, lscale_only, lrestart
      use mpi_params                                                            ! MPI

      implicit none
#ifdef MPI_OPT
      include 'mpif.h'                                                          ! MPI
#endif
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, parameter :: arg_len = 200
      integer :: numargs, iseq, ierr
      character*(arg_len) :: arg1_input
      character*(arg_len), allocatable, dimension(:) :: args
C-----------------------------------------------
#ifdef MPI_OPT
!
!     Call MPI Initialization routines:                                          ! MPI
!

      call MPI_INIT( ierr_mpi )                                                  ! MPI

      if (ierr_mpi .ne. 0) stop 'MPI_INIT error in stellopt'
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr_mpi )                       ! MPI

      if (ierr_mpi .ne. 0) stop 'MPI_COMM_RANK error in stellopt'
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr_mpi )                   ! MPI

      if (ierr_mpi .ne. 0) stop 'MPI_COMM_SIZE error in stellopt'
#endif
!     read input file name from command line
!      
      lscale_only = .false.
      lone_step = .false.
      lrestart = .false.
      
#if defined(LINUX) && defined(MPI_OPT)
      if( myid .eq. master) then
#endif
         call getcarg(1, arg1_input, numargs)
         allocate( args(numargs))
         args(1) = arg1_input

         do iseq = 2, numargs
            call getcarg(iseq, args(iseq), numargs)
         enddo

         if (numargs .lt. 1 .or. args(1) .eq. '-h' .or. 
     1      args(1) .eq. '/h') then
            print *,
     1      ' ENTER INPUT FILE NAME (OR LIST OF NAMES) OR INPUT-FILE',
     2      ' SUFFIX(ES) ON COMMAND LINE'
            print *,' For example: '
            print *,' xstellopt input.tftr (or tftr or ../input.tftr)'
            print *,' Optional control parameters may be entered:'
            print *,' xstellopt {(-o)ne_step} {(-s)cale_only}',
     1      '{(-r)estart} [list of files]'
            print *,' where:'
            print *,' (-o)ne_step:   runs input file(s) one step only'
            print *,'          (effect is same as setting NITER_OPT=1'
            print *,'               in the input file)'           
            print *,
     1      ' (-s)cale_only: scales entries in input files according'
            print *,'        to r00_scale, b00_scale values, writes'
            print *,'        new input file, but does NOT execute'
            print *,'        any optimization steps'
            print *,' (-r)esume: resumes genetic algorithm or'
            print *,' diff. evolution from previous population'
            stop 'Invalid command line'
         endif

#if defined(LINUX) && defined(MPI_OPT)
      endif

      call MPI_BCAST(numargs,1,MPI_INTEGER, master, MPI_COMM_WORLD, 
     1               ierr)

      if( myid .ne. master) then
         allocate( args(numargs))
      endif

      do iseq = 1, numargs
         call MPI_BCAST(args(iseq), arg_len, MPI_CHARACTER,master,
     1                  MPI_COMM_WORLD, ierr)
      enddo
#endif


      do iseq = 1, numargs
         if (args(iseq)(1:2).eq."-s" .or. args(iseq)(1:2).eq."-S"
     1   .or. args(iseq)(1:2).eq."/s" .or. args(iseq)(1:2).eq."/S")
     2   lscale_only = .true.          
         if (args(iseq)(1:2).eq."-o" .or. args(iseq)(1:2).eq."-O"
     1   .or. args(iseq)(1:2).eq."/o" .or. args(iseq)(1:2).eq."/O")
     2   lone_step = .true.          
         if (args(iseq)(1:2).eq."-r" .or. args(iseq)(1:2).eq."-R"
     1   .or. args(iseq)(1:2).eq."/r" .or. args(iseq)(1:2).eq."/R")
     2   lrestart = .true.
      end do

!
!     read each input file name from command line & optimize it
!      
      do iseq = 1, numargs
         if (args(iseq)(1:1).ne.'-' .and. args(iseq)(1:1).ne.'/')
     1       call optimize(args(iseq))
      end do

      deallocate(args)

#ifdef MPI_OPT
      call MPI_FINALIZE(ierr_mpi)      !Close out MPI                             !MPI
      if (ierr_mpi .ne. 0) stop 'MPI_FINALIZE error in stellopt'
#endif
      end program stellarator_optimizer
      

      subroutine optimize(in_file)
      use optim
      use boozer_params, ONLY: xm_bdy, xn_bdy, rmnc_bdy, zmns_bdy,
     1                         rmnc_opt, zmns_opt, lmns_opt
      use vparams, only: zero
      use system_mod
      use optim_params, ONLY: niter_opt
      use legendre_params
      use mpi_params   
      implicit none
#ifdef MPI_OPT
      include 'mpif.h'                                                   !MPI
#endif
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      character*(*) :: in_file
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: k, nvar, nopt, lwa, istat, info
      real(rprec), dimension(nxc) :: xc_opt
      real(rprec) :: tstart, tstop
      character :: scratch_dir*100, temp*200
      logical :: lprint
C----------------------------------------------- 
      iter_min = 0

      call get_extension(in_file)
      scratch_dir = "stellopt" // "_" // trim(seq_ext)
!
!     DO ALL CALCULATIONS IN SCRATCH DIRECTORY
!
      if(myid .eq. master) then                                          !START MPI
         temp = '/bin/rm -Rf ' // scratch_dir
         call system(temp)
         temp = '/bin/mkdir -m 755 ' // scratch_dir
         call system(temp)
         k = chdir(scratch_dir)
         if (k .ne. 0) then
            print *,
     1      'ierr = ',k,': Unable to chdir to working directory!'
            stop
         end if                                                             

!
!     COPY INPUT FILE FROM ORIGINAL DIRECTORY
!
         temp = "/bin/cp ../" // trim(in_file) // " ./" // trim(in_file)
         call system(temp)
      end if                                                             !END MPI

#ifdef MPI_OPT

      call MPI_BARRIER(MPI_COMM_WORLD, ierr_mpi)                         !MPI

      if (ierr_mpi .ne. 0) stop 'MPI_BARRIER error in OPTIMIZE'
      if (myid .ne. master) k = chdir(scratch_dir)                       !MPI
#endif
!
!     Initialize variables, read input file. On exit, in_file is the input file name
!     Set up input, output file extensions
!
      xc_opt = zero
      call osetup(in_file)
!
!     SETUP OPTIMIZATION ROUTINE PARAMETERS
!
      if (.not.lscale_only) then
         call initialize_opt (xc_opt, nvar, nopt)
         lwa = (nopt+5)*nvar + nopt
!
!        RUN OPTIMIZE ROUTINE
!
         call second0 (tstart)
      
         call run_optimizer (xc_opt, nopt, nvar, lwa, info)
      else
         call write_indata (min_input_file, info)
      end if

      call second0 (tstop)

      deallocate (rmnc_bdy, zmns_bdy, xm_bdy, xn_bdy)
      deallocate (rmnc_opt, zmns_opt, lmns_opt)
      deallocate (ns_surf, ns_booz, lneed_modB)
      deallocate (ns_ball, ns_neo)                                       !VMECCOBRA, NEO
      if(l_legendre .and. lcur_prof_opt)                                 !LEGENDRE
     1  deallocate (a_leg, b_leg, a_leg_inv, b_leg_inv, tc)              !LEGENDRE
      if(l_legendre .and. liota_prof_opt)                                !LEGENDRE
     1  deallocate (a_leg, b_leg, a_leg_inv, b_leg_inv, ti)              !LEGENDRE
      if(l_legendre .and. lpres_prof_opt)                                !LEGENDRE
     1  deallocate (a_leg, b_leg, a_leg_inv, b_leg_inv, tm)              !LEGENDRE
      
      if (myid .eq. master) then                                         !START MPI      
         lprint = .not.lniter1
         if (lprint) print 90, chisq_min, iter_min, trim(min_input_file)

         print *, 'All input/output files for this run are stored in',
     1          ' the directory: '
         print *, trim(scratch_dir)
!     
!        APPEND TIME, CHI-SQ TO END OF OUTPUT FILE
!
         open (unit=iunit_opt, file=output_file, status='old',
     1         position='append', action='readwrite', iostat=istat)
         if (istat .eq. 0) then
            write(iunit_opt, '(/,a,i4,a,i6)')' Number of parameters = ',
     1      nvar,' Number of constraints = ',nopt 
            if (lprint) write (iunit_opt, 90)
     1          chisq_min, iter_min, trim(min_input_file)
            write (iunit_opt, 100) tstop-tstart
         end if

         close (iunit_opt, iostat=istat)

         print  100, tstop - tstart
!
!        CHANGE FILE EXTENSIONS IN INPUT FILES, OUTPUT FILES TO ELIMINATE _OPT
!        (OTHERWISE, INPUT FILES WILL NOT RUN CORRECTLY WITH THE WOUT FILE, WHICH
!         HAS EITHER THE .MIN EXTENSION OR NONE AT ALL)
!
         temp = ""
         if (lprint) temp = ".min"
         call clean_up_extensions (min_ext, temp)

!
!        CONCATENATE INFO FROM SEVERAL FILES INTO A DATA_SUMMARY FILE
!
         call data_harvest (min_ext, output_file, nvar, nopt)         
         
      endif                                                              !END MPI

 90   format(/,' Minimum Chi-sq = ',1pe10.3,' at iteration ',i6,/,
     1  ' Minimum Chi-sq state stored in file ',a)
 100  format(/,' SECONDS IN OPTIMIZATION CODE:  ',1pe12.3,' SEC')

      end subroutine optimize


      subroutine get_extension(in_file)
      use optim_params, ONLY: seq_ext
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      character*(*) :: in_file
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      character*6, parameter :: input_ext = 'input.'
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: index_end, index_dat
C----------------------------------------------- 
      index_dat = index(in_file,INPUT_EXT)
      if (index_dat .eq. 0) then
         in_file = input_ext // trim(in_file)
         index_dat = 1
      end if   
      index_dat = index_dat + len(input_ext)
      index_end = len_trim(in_file)
      seq_ext   = in_file(index_dat:index_end)
      
      end subroutine get_extension
      
      
      subroutine osetup(in_file)
      use optim
      use vmec_input
c     use vmec_input, ONLY: ns_array, mpol, ntor, gamma, spres_ped, 
c    1    mgrid_file, lmac, lfreeb, rbc, extcur, raxis, zaxis, 
c    2    imatch_phiedge, bloat, ldiagno, ftol_array, loldout
      use boozer_params, ONLY: xm_bdy, xn_bdy, rmnc_bdy, zmns_bdy,
     1    rmnc_opt, zmns_opt, lmns_opt
      use vparams, ONLY: one, zero, ntord, mpol1d
      use system_mod
      use bootsj_input, only: damp_bs
      use mpi_params                                                     !MPI
      implicit none
#ifdef MPI_OPT
      include 'mpif.h'                                                   !MPI
#endif
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      character*(*) :: in_file
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: mpol_default = 6
      integer, parameter :: ntor_default = 0

      real(rprec), parameter :: p360 = 360
      character*6, parameter :: input_ext = 'input.'
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: js, k, iunit, ierr
      integer :: istat
      logical :: lexist, ljacobian
      integer, dimension(ini_max):: zetain, thetain                      !COBRA
      character*200 :: stripped_input_file, temp
      character*20 :: ch_myid
C----------------------------------------------- 

      call getenv('HOME', home_dir)
#if defined(FALCON)
      home_dir = trim(home_dir) // '/falcon'
#else
      home_dir = trim(home_dir) // '/bin'
#endif      
      
!
!     DELETE OUTPUT FILE IF IT ALREADY EXISTS
!
      
      output_file = 'output.' // trim(seq_ext)
      inquire (file=trim(output_file), exist=lexist)
      if (lexist) then
         open(unit=iout, file=output_file)
         close (iout, status='delete', iostat=istat)
      end if   
         
!
!     Default values for INDATA parameters
!
      mgrid_file = ' '
      lfreeb = .false.
      delt = 0
      tcon0 = 0
      nfp = 0
      ncurr = 0
      mpol = mpol_default
      ntor = ntor_default
      ntheta = 0
      nzeta = 0
      niter = 0
      nstep = 0
      nvacskip = 0
      imatch_phiedge = 1
      ns_array = 0
      ftol_array = 0
      gamma = 0
      spres_ped = 1
      curtor = 0
      phiedge = 0
      extcur = 0
      am = 0
      ai = 0
      ac = 0
      ac_form = 0
      aphi = 0
      raxis = 0
      zaxis = 0
      rbc = 0
      zbs = 0
      nextcur_opt = 0
      nextcur_vmec = 0
      lmac = .false.
     
      bloat = 1
      ldiagno = .false.
      loldout = .false.

!
!     SET DEFAULTS FOR OPTIMUM PARAMETERS
!
      lbeta_min = .false.
      laspect_max = .false.
      lcoil_geom = .false.
      sigma_berr_avg = bigno
      sigma_berr_max = bigno
      coil_separation = zero
      NumJstar = 4
      NumJinvariant = 4
      NPitch_JConf = 4
      NS_JConf_Src = 0
      NS_JConf_Tgt = 0
      target_iota_max = 1
      target_iota_max_min = 1
      target_iota_min = 1
      sigma_iota_max = bigno
      sigma_iota_max_min = bigno
      sigma_iota_min = bigno
      sigma_iota = bigno
      sigma_jedge = bigno
      target_jedge = 0
      sigma_vp   = bigno
      sigma_mercier = bigno
      target_kink = zero
      sigma_kink = bigno
      sigma_curv = bigno
      sigma_aspect = bigno
      Target_AspectRatio = 3
      sigma_coil_complex = bigno
      sigma_MaxCurrent = bigno
      target_MaxCurrent = 0
      sigma_beta = bigno
      target_beta = 0
      sigma_centering = bigno
      sigma_rmin = bigno
      sigma_rmax = bigno
      sigma_zmax = bigno
      target_rmin = 0
      target_rmax = 0
      target_zmax = 0
      sigma_ellipticity = bigno
      target_ellipticity = 0
      sigma_bmin = bigno
      sigma_bmax = bigno
      sigma_bmn  = bigno
      sigma_jstar = bigno
      sigma_jinvariant = bigno
      sigma_jconf = bigno
      sigma_ripple = bigno
      sigma_bootsj = bigno
      sigma_jac = bigno
      n_jac = 0
      m_jac =0
      sigma_vac_island = bigno
      n_vac_island = 0
      m_vac_island = 0
#ifdef MCURIE
#define MAX_PROCESS 32
#elif defined(CRAY)
#define MAX_PROCESS 8
#elif defined(OSF1)
#define MAX_PROCESS 1
#else
#define MAX_PROCESS 1
#endif
      num_levmar_params = max(3, MAX_PROCESS)           !!Number of lm parameter estimates/jacobian evaluation
      num_processors = MAX_PROCESS                        !!Number of processors to distribute jacobian evaluations over
#ifdef MPI_OPT
      num_levmar_params = numprocs         !!Number of lm parameter estimates/jacobian evaluation
      num_processors = numprocs-1                        !!Number of processors to distribute j
#endif
!
!  Experiment matching  
!
      sigma_eplasma = bigno           ! MCZ plasma kinetic stored energy
      target_eplasma = 0
      sigma_curtor = bigno            ! MCZ plasma current
      target_curtor = 0

      sigma_p_prof = bigno            ! MCZ experimental pressure profile
      sigma_p_damp = bigno
      p_prof = 0
      np_prof = 0
      factor_p_prof = 1
      lp_prof_incl_edge = .true.
      r_p_prof = 0
      z_p_prof = 0
      phi_p_prof = 0

      sigma_emis = bigno              ! MCZ SXR emissivity chords
      sigma_emis_damp = bigno
      n_emis = 0
      emis_file = " "
      aemis = 0

      ldiagno_opt = .false.           ! diagno simulations of mag. diags.
      diagno_control = 'diagno.control'
      sigma_diagno_seg = bigno
      target_diagno_seg = 0
      sigma_diagno_flx = bigno
      target_diagno_flx = 0

! 
!     LPK ADDITIONS
! 
      sigma_fluxp = bigno
      target_fluxp = 0
      sigma_pseudo = bigno
      sigma_pseudo2 = bigno
      lpseudo_sin = .false.
      nproc = 1                       
      nbmn = 0                        
      lkink_opt = .false.             
      lbal_opt = .false.              
      lbootstrap = .false.
      jboot = 0
      fboot = 0
      fboot(0) = 1
      zeff_boot = 1
      lseedcur = .false.              
      lpress_opt = .false.            
      lcoil_complex = .false.
      lnescoil_opt = .false.
      lcoil_opt = .false.      ! obsolete !
      sigma_coil_jmax = bigno
      sigma_berr_ave = bigno
      target_coil_complex = 1 
      target_coil_jmax = 0
      ldkes_opt = .false.
      ldkes_mask = .false.
      sigma_dkes = bigno
      ndkes_mask = 0
      dkes_efield = 0
      dkes_nu = 0.01_dp
      aseedcur = 0
      Target_RBtor = 1
      Sigma_RBtor = bigno
!!
!!    NEO    
!!
      sigma_neo = bigno
      lneo_opt = .false.
      nneo_mask = zero
      lneo_mask = .false.
      ldiag_opt = .false.
      lkeep_mins = .false.
      ldsubr_opt = .false.
      sigma_dsubr = bigno
      lorbit_opt = .false.
      sigma_orbit = bigno
!!
!!    END LPK ADDITIONS
!!
      Target_Well = (/ 0.0_dp, -0.39_dp, 0.19_dp, (0.0_dp,istat=1,8) /)
      niter_opt = 1
      nopt_alg = 0
      nopt_boundary = 0
      mboz_opt = 0
      nboz_opt = 0
      lreset_opt = .true. 
      lcur_prof_opt = .false.
      lcurprof_opt = .false.       ! obsolete !
      lcur_opt_edge0 = .false.
      liota_prof_opt = .false.
      lbmn = .false.
      lj_star = .false.
      lj_invariant = .false.
      lprof_opt = .false.
      lextcur = .false.
      lbootsj_opt = .false.
      lsurf_mask = .false.
      lfix_ntor = .false.
      nsurf_mask = 0
      epsfcn = -1
      helicity = cmplx(0.0_dp, 0.0_dp)
      sym_type = 'NONE'

!
!     Default values for COILSIN parameters
!
      call initialize_coilsin


!
!     Default values for BOOTIN parameters
!
      damp_bs = -1
!---------------------------------------------------------------------------------
!   CODE added by R.SANCHEZ (01/19/99): ballooning related variables and sigmas.

      target_balloon = 0
      lballoon_flip = .false.
      sigma_balloon = bigno
      sigma_bal = bigno                                                  !Old style - LPK
      bal_zeta0 = 0
      bal_theta0 = 0
      sigma_pgrad = bigno
      sigma_pedge = bigno
      lballoon_opt = .false.
      lballoon_mask = .false.
      lpres_prof_opt = .false.
      lpres_opt_edge0 = .false.
      lpres_opt_edgegr0 = .false.
      pres_opt_nmax = 11

      nballoon_mask = 0                                                  !VMECCOBRA
      l_legendre = .false.                                               !LEGENDRE

!------------------------------------------------------------------------------
!  initialize 3D boundary limiter parameters,  RH & MZ  June 2000

      lvv_tgt_min = .false.
      target_vv = 0
      target_vv_rms = 0
      sigma_vv = bigno
      sigma_vv_rms = bigno
      sigma_vv_max = bigno
      rbc_vv = 0
      zbs_vv = 0
      mpol_vv = 0
      ntor_vv = 0
      nu_vv = 5
      nv_vv = 2

      target_bd = 0
      target_bd_rms = 0
      sigma_bd = bigno
      sigma_bd_rms = bigno
      sigma_bd_max = bigno
      rbc_bd = 0
      zbs_bd = 0
      mpol_bd = 0
      ntor_bd = 0
      nu_bd = 5
      nv_bd = 2

      r00_opt = -1
      r00_scale = 1
      b00_scale = 1
      rmax_opt = 0; rmin_opt = 0; zmax_opt = 0
      rgrid_max = 0; rgrid_min = 0; zgrid_max = 0; zgrid_min = 0

      shapeweight = .FALSE.
      theta0_bw(1) = 0 ; theta0_bw(2) = 0 ; theta0_bw(3) = 0
      phi0_bw = 0 ; wtheta_bw = 0.392699 ; wphi_bw = 0.392699
      planes_bw(1) = 0 ; planes_bw(2) = 1.5707 ; planes_bw(3) = 3.14159

      phi_lim = -400  ! piecewise linear limiters
      r_lim = 0
      z_lim = 0
!------------------------------------------------------------------------------

      sigma_kappa = bigno
      target_kappa = 0
      Target_Iota = 0
      target_iota_p = 0
      Target_Iota_min = 0
      sigma_iota_pmax = bigno
      sigma_iota_pmin = bigno
      sigma_extcur = bigno
      target_extcur = 0
      sigma_oh = bigno
      oh_coefs = 0
!------------------------------------------------------------------------------

      call GA_preset           !  initialize the GA parameters
      call DE_preset           !  initialize the DE parameters

!------------------------------------------------------------------------------
!
!     read input (optimum) namelist
!     first, strip any comments from input file (F95-not necessary)
!

!     Produce clean file input_file//'.stripped'      
      if (myid .eq. master) call strip_comments(in_file)                 !MPI
#ifdef MPI_OPT

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)                             !MPI
#endif
      stripped_input_file = trim(in_file) // '.stripped' 
      write (ch_myid, *) myid
      temp = "/bin/cp " // trim(stripped_input_file) // " " //
     1   trim(stripped_input_file)//adjustl(ch_myid)
      call system(trim(temp))

      temp = stripped_input_file
      stripped_input_file = trim(stripped_input_file)//adjustl(ch_myid)
      call read_input_opt (stripped_input_file, istat)  !!read in_file and deletes it
      if (istat .ne. 0) then
         print *,'Error reading input namelist in STELLOPT',
     1           ' osetup routine: istat = ', istat
         stop
      end if   
#ifdef MPI_OPT

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif      
      if (myid .eq. master) then                                         !MPI
         temp = "/bin/rm " // trim(temp)
         call system (temp)
         if (lfreeb .and. .not.lcoil_geom) then
            k = index(mgrid_file,'/')          !!If NOT full path, must link to scratch directory
            if (k .eq. 0) then
               temp = "/bin/ln -s ../" // trim(mgrid_file) // " ./" 
     1           // trim(mgrid_file)
               call system(temp)
            end if
         end if
      end if                                                             !MPI
#ifdef MPI_OPT

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif 
      min_ext = trim(seq_ext)
      lniter1 = ((niter_opt.eq.1 .and. nopt_alg.eq.0) .or. 
     1            lone_step .or. lscale_only)
      if (.not.lniter1) min_ext = trim(min_ext) // '.min'
      min_input_file = input_ext // trim(min_ext)
      min_wout_file = 'wout.' // trim(min_ext)

      if (r00_opt.gt.0 .and. rbc(0,0).ne.zero) 
     1    r00_scale = r00_opt/rbc(0,0)
      if (abs(r00_scale*b00_scale) .ne. one)
     1   call init_scale(r00_scale, b00_scale)

      raxis_old = raxis(:,1)                                            !Save for writing into restart (.min) file
      zaxis_old = zaxis(:,1)
      jcurv_opt = 0;   iota_opt = 0
!
!     COMPATABILITY WITH LPK ADDITIONS
!
      if (lcoil_complex .or. lcoil_opt) lnescoil_opt = .true.
      if (lbootstrap) lbootsj_opt = .true.
      if (lbal_opt) lballoon_opt = lbal_opt
      if (nproc .ge. num_processors) num_processors = nproc
      target_coil_complex = max(one, target_coil_complex)
      if (damp_bs .eq. zero) damp_bs = -1                                !Old-style file: did not mean 0...
      if (mpol_vv.gt.mpol1d .or. ntor_vv.gt.ntord) stop
     1    'mpol_vv > mpol1d or ntor_vv > ntord exceed bounds'
      if (mpol_bd.gt.mpol1d .or. ntor_bd.gt.ntord) stop
     1    'mpol_bd > mpol1d or ntor_bd > ntord exceed bounds'

      NumJstar = min(NumJstard, NumJstar)
      NumJinvariant = min(NumJstard, NumJinvariant)
      mpol1_opt = mpol - 1                                               !need in initialize_opt
      ntor_opt  = ntor
      mnmax_opt = ntor + 1 + (mpol1_opt)*(1 + 2*ntor)  

      if (lcoil_geom) lfreeb = .true.

      if (lfreeb) then
         inquire(file=trim(mgrid_file), exist=lexist)
         if (.not.lexist .and. myid.eq.master .and. 
     1       .not.lcoil_geom) then
            print *,' Could not locate mgrid file!'
            lfreeb = .false.
         end if
         nextcur_opt = count(lextcur)
         if (lfreeb .and. lnescoil_opt .and. myid.eq.master) print *,
     1   ' NESCOIL Constraints are being used in free-boundary mode!'         
      end if      

!
!     Backward compatability
!
      select case(trim(sym_type))
        case ('QA')
           helicity = cmplx(one,zero)
        case ('QH')
           helicity = cmplx(one,-one)
        case ('QP')
           helicity = cmplx(zero,one)
      end select

!
!     Check consistency of logical variables with sigmas
!
      if (all(sigma_bmn .ge. bigno) .and. lbmn) then
         if (myid .eq. master) 
     1   print *,' LBMN is being set to FALSE; specify SIGMA_BMN array'
         lbmn = .false.
      endif
      if (all(sigma_bootsj .ge. bigno) .and. lbootsj_opt) then
         if (myid .eq. master) 
     1   print *,' LBOOTSJ_OPT is being set to FALSE; ',
     2           'specify SIGMA_BOOTSJ array'
         lbootsj_opt = .false.
      endif
      if (all(sigma_balloon .ge. bigno) .and. lballoon_opt) then
         if (myid .eq. master) 
     1   print *,' LBALLOON_OPT is being set to FALSE; ',
     2           'specify SIGMA_BALLOON array'
         lballoon_opt = .false.      
      endif         
      if (all(sigma_dkes .ge. bigno) .and. ldkes_opt) then 
         if (myid .eq. master)
     1   print *,' LDKES_OPT is being set to FALSE; ',
     2           'specify SIGMA_DKES array'
         ldkes_opt = .false.
      endif
         
      
      if (real(helicity) .lt. zero) helicity = -helicity
      if (lbmn .and. helicity.eq.cmplx(zero,zero)) then
      if (myid .eq. master) then                                         !MPI      
        print *,' Must specify helicity (symmetry type) for lbmn = T'
        print *,' Allowable values for helicity are:'
        print *,'    (1,0)   =>   quasi-axisymmetry'
        print *,'    (0,1)   =>   quasi-poloidal symmetry'
        print *,'    use with quasi-omnigeneity to reduce 1/R drift'
        print *,'    (l,k)   =>   quasi-helical symmetry'
        print *,'    |B| = F(lu + kv)'
      endif                                                              !MPI
        stop 'helicity MUST be supplied in OPTIMUM namelist'
      end if

      nrad = maxval(ns_array)

      allocate (rmnc_bdy(mnmax_opt), zmns_bdy(mnmax_opt), 
     1          xm_bdy(mnmax_opt), xn_bdy(mnmax_opt), 
     2          rmnc_opt(mnmax_opt,nrad), zmns_opt(mnmax_opt,nrad), 
     3          lmns_opt(mnmax_opt,nrad),
     4          ns_surf(nrad), ns_booz(nrad), ns_ball(nrad), 
     5          ns_neo(nrad), lneed_modB(nrad), stat = istat)
      if (istat.ne.0 ) stop 'Allocation error in OSETUP'
!
!     Make old-style, new-style masks conform
!
      call update_style (lsurf_mask, nsurf_mask, nrad)
      call update_style (ldkes_mask, ndkes_mask, nrad)
      call update_style (lballoon_mask, nballoon_mask, nrad)             !VMEC COBRA (RS)
      call update_style (lneo_mask, nneo_mask, nrad)                     !NEO

!     if tailoring the pressure profile: all surfaces must be included to compute gradient
      if(lpres_prof_opt .and. lballoon_opt) lballoon_mask = .true.     

      lneed_modB = .false.
      lballoon_mask(1) = .false.                !!MUST ignore axis point in COBRA   (RS)
      lballoon_mask(nrad) = .false.             !!MUST ignore last point in COBRA   (RS)

      ns_surf_max = 0                           !!Boozer surfaces based on lsurf_mask
      ns_booz_max = 0                           !!TOTAL NUMBER BOOZER SURFACES
      ns_ball_max = 0                           !!COBRAVMEC   (RS)
      ns_neo_max  = 0 
      ljacobian = any(abs(sigma_jac(:)) < bigno)
      
      do js = 2, nrad
         if (lsurf_mask(js)) then
           ns_surf_max = ns_surf_max + 1
           ns_surf(ns_surf_max) = js
         end if  
!!    NEO
         if (lneo_mask(js)) then
           ns_neo_max = ns_neo_max + 1
           ns_neo(ns_neo_max) = js
         endif

         if (lsurf_mask(js) .or. ldkes_mask(js) .or.
     1       ljacobian .or. lneo_mask(js) .or. 
     2       (sigma_jconf<bigno .and. js >= ns_jconf_src .and. 
     3                         js <= ns_jconf_tgt) ) then
           ns_booz_max = ns_booz_max + 1
           ns_booz(ns_booz_max) = js
         end if 
 
         if(lballoon_mask(js)) then                                      !COBRAVMEC   (RS)
           ns_ball_max = ns_ball_max + 1                                 !COBRAVMEC   (RS)
           ns_ball(ns_ball_max) = js                                     !COBRAVMEC   (RS)
         endif                                                           !COBRAVMEC   (RS)         
      end do                                  

!
!     IF USER FORGOT TO SPECIFY LBALLOON_MAKS ARRAY, DEFAULT TO NS_SURF...
!
      if (lballoon_opt .and. ns_ball_max.eq.0) then
         ns_ball_max = ns_surf_max
         ns_ball(1:ns_ball_max) = ns_surf(1:ns_surf_max)
         if (myid .eq. master) then                                      !MPI      
            print *,'**********************************************'
            print *,'WARNING:'
            print *,' LBALLOON_OPT = TRUE but no nonzero values for ',
     1              'NBALLOON_MASK were found'
            print *,' Using NSURF_MASK for ballooning surfaces'
            print *,' If this is wrong, set LBALLOON_OPT = FALSE,'
            print *,' or set the NBALLOON_MASK array values in ',
     1              'the OPTIMUM namelist'
            print *,'**********************************************'
         endif                                                           !MPI
      end if
!
!   if feeding back on pressure profile, must vary some pressure coefs
!
      if( lpres_prof_opt) then
         if( pres_opt_nmax < 0) then
            print *,'**********************************************'
            print *,'WARNING:'
            print *,' LPRES_PROF_OPT = TRUE but pres_opt_nmax < 0 '
            print *,' Incompatible!!'
            print *,' PRES_OPT_NMAX set to 1'
            print *,'**********************************************'

            pres_opt_nmax = 1
         else if( pres_opt_nmax > 10) then
            pres_opt_nmax = 10
            
         endif
      endif

!
!   set np_prof on basis of sigma_p_prof
!

      if( np_prof == 0 .and. any(sigma_p_prof<bigno)) then
         np_prof = 2

         do while(any(sigma_p_prof(np_prof:)<bigno) .and. 
     1            np_prof < size(sigma_p_prof))
            np_prof = np_prof + 1
         enddo
         np_prof = np_prof - 1
      endif

!
!  set ndiagno_seg and ndiagno_flx from respective sigmas
!

      if( ldiagno_opt .and. any(sigma_diagno_seg < bigno)) then
         ndiagno_seg = count(sigma_diagno_seg < bigno)

         do while( any(sigma_diagno_seg(ndiagno_seg:)<bigno) .and.
     1             ndiagno_seg < size(sigma_diagno_seg))
            ndiagno_seg = ndiagno_seg + 1
         enddo
         ndiagno_seg = ndiagno_seg - 1
      else
         ndiagno_seg = 0
      endif


      if( ldiagno_opt .and. any(sigma_diagno_flx < bigno)) then
         ndiagno_flx = count(sigma_diagno_flx < bigno)

         do while( any(sigma_diagno_flx(ndiagno_flx:)<bigno) .and.
     1             ndiagno_flx < size(sigma_diagno_flx))
            ndiagno_flx = ndiagno_flx + 1
         enddo
         ndiagno_flx = ndiagno_flx - 1

         if((.not. lpres_opt_edge0) .and. myid==master) then
            print *,'************************************************'
            print *,'WARNING:'
            print *,'Matching flux loops, but LPRES_OPT_EDGE0 = FALSE'
            print *,'=> will not see all of the plasma pressure'
            print *,'************************************************'
         endif

      else
         ndiagno_flx = 0
      endif

      if( ldiagno_opt .and. ndiagno_seg + ndiagno_flx > 0) then
         ldiagno = .true.
      endif


!---------------------------------------------------------------------------------
!   CODE added by R.SANCHEZ (02/01/99):  count number of initial theta values in degrees
!              (NINI_THETA), initial zeta values in degrees (NINI_ZETA) and total number
!              of initial positions (NINI_TOT) in namelist OPTIMUM where
!              ballooning growth rates are to be evaluated.

      if (lballoon_opt .or. lpres_prof_opt) then
         if(MINVAL(bal_theta0).lt.zero  .or. MAXVAL(bal_theta0).ge.p360 
     1   .or. MINVAL(bal_zeta0).lt.zero .or. MAXVAL(bal_zeta0) .ge.p360)
     2   then
      if (myid .eq. master) then                                            ! MPI      
           print *, 'ALL initial angles (in degrees) MUST be',
     1            ' 0 <= angle < 360 in STELLOPT input file!'
      endif                                                                 ! MPI
           stop
         end if  
         nini_theta=1
         nini_zeta=1
         thetain(1)=bal_theta0(1)
         zetain(1)=bal_zeta0(1)
         do k=2, ini_max
            if (bal_theta0(k) .ge. zero) then
               do js=1, k-1
                  if(bal_theta0(k).eq.bal_theta0(js)) exit
                  if(js .eq. k-1) then
                     nini_theta=nini_theta+1
                     thetain(nini_theta)=bal_theta0(k)
                  endif
               enddo
            endif
            if (bal_zeta0(k) .ge. zero) then
               do js=1, k-1
                  if (bal_zeta0(k).eq.bal_zeta0(js)) exit
                  if (js .eq. k-1) then
                     nini_zeta=nini_zeta+1
                     zetain(nini_zeta)=bal_zeta0(k)
                  endif
              enddo
            endif
         enddo          
         bal_theta0(1:nini_theta)=thetain(1:nini_theta)
         bal_zeta0(1:nini_zeta)=zetain(1:nini_zeta)
         nini_tot=nini_zeta*nini_theta
      endif
!--------------------------------------------------------------------------------- 
!
!     Maintain Backwards compatability
!
      if (sigma_rmin .ge. bigno) sigma_rmin = sigma_centering
      if (sigma_rmax .ge. bigno) sigma_rmax = sigma_centering
      
      if (epsfcn .eq. -1 ) then
         if (lreset_opt) then
            epsfcn = 1.e-4_dp
         else
            epsfcn = 1.e-3_dp
         endif
      endif
 
      end subroutine osetup   


      subroutine update_style (lmask, surf_mask, nrad)
      use kind_spec, only: rprec
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: nrad
      real(rprec) :: surf_mask(nrad)
      logical :: lmask(nrad)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: k, js
      logical :: l1
C-----------------------------------------------
      l1 = all(.not.lmask)               !!Check if user is using old-style input
      if (l1) then                       !!New style, find lsurf_mask array
         do k = 1, nrad
            do js = 2, nrad
            if (js .eq. (1+ nint(surf_mask(k)*(nrad-1))))
     1            lmask(js) = .true.
            end do
         end do
      else                               !!Old style, convert lsurf to nsurf
         k = 1
         do js = 2, nrad
            if (lmask(js)) then
              surf_mask(k) = real(js - 1,rprec)/(nrad-1)
              k = k+1
            end if
         end do
      end if
      lmask(1) = .false.           !!MUST ignore axis point in boozer code

      end subroutine update_style


      subroutine read_input_opt (file, istat)
      use optim_params, only: lcoil_geom, sigma_berr_avg, sigma_berr_max
      use optim, only: bigno
      use mpi_params    
      use read_namelist_mod
      use safe_open_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(out) :: istat
      character*(*), intent(in) :: file
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: iunit=7
C-----------------------------------------------
      call safe_open (iunit, istat, file, 'old', 'formatted')
      if (istat .ne. 0) goto 100

      call read_namelist (iunit, istat, 'indata')
      if (istat .ne. 0) then
         if (myid .eq. master) then                                      !MPI      
         print *,' indata namelist read error in READ_INPUT_OPT,',
     1       ' istat = ', istat         
         endif                                                           !MPI
         goto 100
      end if    

      call read_namelist (iunit, istat, 'optimum')
      if (istat .ne. 0) then
         if (myid .eq. master) then                                      !MPI
            print *,' optimum namelist read error in READ_INPUT_OPT,',
     1      ' istat = ', istat            
         endif                                                           !MPI
         go to 100
      end if    
     
      call read_namelist (iunit, istat, 'bootin')

      call read_namelist (iunit, istat, 'ga_de')

      istat = 0
      if (lcoil_geom .or. sigma_berr_avg < bigno .or.
     1    sigma_berr_max < bigno ) call read_namelist (iunit, 
     2        istat, 'coilsin')
      if (istat .ne. 0) then
         if (myid .eq. master) then                                      !MPI
            print *,' coilsin namelist read error in READ_INPUT_OPT,',
     1      ' istat = ', istat            
         endif                                                           !MPI
         stop
      end if    

 100  continue

      close (iunit, status='delete')

      end subroutine read_input_opt


      subroutine init_scale(r_scale, b_scale)
      use vmec_input, ONLY: rprec, curtor, am, raxis, zaxis, phiedge,
     1     extcur, rbc, zbs, dp, lfreeb
      use optim_params, ONLY: Target_MaxCurrent, coil_separation
      use mpi_params                                                     !MPI
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec) :: r_scale, b_scale
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, one = 1
C-----------------------------------------------
      if (rbc(0,0) .le. zero) return
      if (r_scale .le. zero) r_scale = 1
      
      if (r_scale.ne.one .and. myid.eq.master) then
         print *,'  R  being scaled by ', r_scale
         if (lfreeb) print *,
     1   ' ***CAUTION*** Free-boundary mode: dimensions in mgrid file',
     2   ' may no longer be invalid!'
      end if
         
      if (b_scale.ne.one .and. myid.eq.master)                           !MPI      
     1   print *,' |B| being scaled by ', b_scale

      phiedge = r_scale*r_scale*b_scale * phiedge                        !Scale |B| by b_scale

      rbc = r_scale*rbc
      zbs = r_scale*zbs

      raxis = r_scale*raxis
      zaxis = r_scale*zaxis
!     coil_separation = r_scale*coil_separation

      Target_MaxCurrent = r_scale*b_scale*Target_MaxCurrent
      curtor = r_scale*b_scale*curtor
      extcur = r_scale*b_scale*extcur
      am = b_scale*b_scale * am
      
!
!     Reset scale factors to unity after scaling performed
!
      r_scale = 1
      b_scale = 1

      end subroutine init_scale
      

      subroutine initialize_opt(xc_opt, nvar, nopt)
      use optim
      use legendre_params                                                !LEGENDRE
      use vmec_input, ONLY : rbc, zbs, ai, am, ac, ncurr, lfreeb,
     1     curtor, extcur, phiedge, nfp, ac_form
      use vparams, ONLY: one, zero
      use coilsnamin, only: lmodular, lsaddle, lbcoil, lvf, lsurfv,
     1    lsadsfv, ltfc, ltfcv, lmodcur, lsadcur, bcoil_file, lbnorm
      use mpi_params                                                     !MPI
      implicit none
#ifdef MPI_OPT
      include "mpif.h"
#endif
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(out) :: nvar, nopt
      real(rprec), dimension(*) :: xc_opt
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nb, mb, ik, nvariables, nvar_pre
      real(rprec), dimension(-ntord:ntord,0:mpol1d) :: 
     1   rbc_in, zbs_in
      real(rprec), dimension(1) :: fvec = (/zero/)
      real(rprec), dimension(:), allocatable :: work
      real(rprec) :: sum1, delta 
      external convert_boundary, convert_boundary_PG,
     1         unique_boundary, unique_boundary_PG
C-----------------------------------------------
      chisq_min = huge(chisq_min)
      min_count = 0
!
!     SET UP CORRESPONDENCE BETWEEN XC_OPT ARRAY (BOUNDARY COEFFICIENTS
!     USED BY OPTIMIZATION ROUTINES) AND XC ARRAY. FIRST CONVERT TO A
!     UNIQUE ANGLE REPRESENTATION AT THE BOUNDARY (BASED ON HIRSHMAN/BRESLAU)*
! 
      irm0_bdy = 0;      izm0_bdy = 0;      irho_bdy = 0;      nvar = 0
      
      if (rbc(0,1).eq.zero .or. zbs(0,1).eq.zero) then
         if (myid .eq. master)                                           !MPI      
     1       print *, 'Neither RBC nor ZBS for m=0, n=1 can be ZERO',
     2                ' in STELLOPT input file!'      
          stop
      end if    
      
      if (.not.lfreeb) then

         rbc_in = rbc
         zbs_in = zbs                                                    !copy and store m=0 components

!
!        check if conversion was made or if original bdy already in proper form
!        IF NOPT_BOUNDARY=0, USE HIRSHMAN/BRESLAU REPRESENTATION
!                        =1, USE PG (PAUL GARABEDIAN) DELTA_MN REPRESENTATION
!
         if (nopt_boundary .eq. 0 ) then
            call convert_boundary(rbc, zbs, rhobc, mpol1d, ntord)
            if (.not.lniter1) then
               call unique_boundary(rbc_in, zbs_in, rhobc, ntord, 
     1                        mpol1d, mpol1d, ntord)
               delta = sum((rbc - rbc_in)**2)/rbc(0,1)**2 
     1               + sum((zbs - zbs_in)**2)/zbs(0,1)**2
     
               if (delta .gt. 1.e-8_dp) write(*,10) 100*(one-delta)
            end if
 10   format(' Input boundary representation was converted!',/,
     1       ' Reliability of conversion = ',f7.2,'%')
         
         do nb = -ntord, ntord
            if (rbc(nb,0).ne.zero .and. (.not.lfix_ntor(nb)) .and.
     1         (nb.ne.0 .or. abs(sigma_rmax).le.1.e6_dp .or.
     2                       abs(sigma_rmin).le.1.e6_dp)) then
               irm0_bdy = irm0_bdy + 1
               nrz0_opt(irm0_bdy) = nb
               xc_opt(irm0_bdy) = rbc(nb,0)
            endif
         end do
         
         izm0_bdy = irm0_bdy
         do nb = -ntord, ntord
            if (zbs(nb,0).ne.zero .and. (.not.lfix_ntor(nb))
     1         .and. (nb.ne.0)) then
               izm0_bdy = izm0_bdy + 1
               nrz0_opt(izm0_bdy) = nb
               xc_opt(izm0_bdy) = zbs(nb,0)
            endif
         end do

         if (izm0_bdy .gt. ntord+ntor1d) then
            if (myid .eq. master)                                        !MPI
     1          print *,' Only one sign of n allowed for m=0'
            stop ' i > 2*(1+ntord) in STELLOPT initialization'
         end if
           
         izm0_bdy = izm0_bdy - irm0_bdy

         do nb = -ntord, ntord
            do mb = 0, mpol1d
               if (rhobc(nb,mb) .ne. zero .and. .not.lfix_ntor(nb)
     1             .and. (mb .ne. 0 .or. nb .ge. 0)) then
                  irho_bdy = irho_bdy + 1
                  ik = irho_bdy + irm0_bdy + izm0_bdy
                  nbrho_opt(irho_bdy) = nb
                  mbrho_opt(irho_bdy) = mb
                  xc_opt(ik) = rhobc(nb,mb)
               endif
            end do
         end do

         else if (nopt_boundary.eq.1 .and. niter_opt.gt.1) then
           call convert_boundary_PG(rbc, zbs, delta_mn, mpol1d, ntord)

           do nb = -ntord, ntord
             do mb = -mpol1d, mpol1d
               if (delta_mn(nb,mb) .ne. zero .and. .not.lfix_ntor(nb)
     1            .and. .not.(nb .eq.0 .and. mb .eq. 0)) then
                  irho_bdy = irho_bdy + 1
                  ik = irho_bdy + irm0_bdy + izm0_bdy
                  nbrho_opt(irho_bdy) = nb
                  mbrho_opt(irho_bdy) = mb
                  xc_opt(ik) = delta_mn(nb,mb)
               endif
             end do
           end do

         end if

         nvar = nvar + irm0_bdy + izm0_bdy + irho_bdy 
         if (myid .eq. master) print *, irm0_bdy+izm0_bdy+irho_bdy,
     1           ' boundary moments'
!
!     Initialize coils data (if berr targeting and lbcoil = true)
!
         if ((sigma_berr_avg < bigno .or. sigma_berr_max < bigno) 
     1           ) then
            allocate( work(10000) )
            nvar_pre = nvar

            call init_coilgeom (nfp)                                    ! sets nfp in boundary module
!
!     Initialize and count modular coil variables
!
            if (lmodular) then
               call init_modular_coils (nvariables, work)
            end if
!
!     Initialize and count saddle coil variables (if lsaddle = true)
!
            if (lsaddle) then
               call init_saddle_coils (nvariables, work)
            end if
!
!     Initialize and count background coil variables (if lbcoil = true)
!
            if (lbcoil) then
               if (myid .eq. master) then
                  call system("/bin/ln -s ../" // trim(bcoil_file) // 
     1                      " ./" // trim(bcoil_file))
               end if
#ifdef MPI_OPT
               call MPI_BARRIER (MPI_COMM_WORLD, ierr_mpi)
#endif
               call init_bg_currents (nvariables, work)
            end if
!
!     Initialize and count vf coil variables (if lvf = true)
!
            if (lvf) then
               call init_vf_currents (nvariables, work)
            end if

            lbnorm = .true.

            deallocate(work)

            if (myid .eq. master) print *, nvar-nvar_pre,
     1           ' coil shape coefs(I)'
         end if


!     Allow coils to determine free-boundary shape
      else if (lfreeb) then
         if( .not. lcoil_geom) then
            nb = 0
            do ik = 1, size(lextcur)
               if (lextcur(ik)) then
                  nb = nb + 1
                  xc_opt(nb) = extcur(ik)
               end if
            end do
            if (nb .ne. nextcur_opt) 
     1         stop 'Error counting external coils!'
            nvar = nvar + nextcur_opt
            if (myid .eq. master .and. nb .ne. 0) 
     1           print *, nextcur_opt,' coil currents'

         else                                                            ! lcoil_geom = T

            nvar_pre = nvar
            call init_coilgeom (nfp)                                    ! sets nfp in boundary module
!
!     Initialize and count modular coil variables
!
            if (lmodular) then
               call init_modular_coils (nvariables, xc_opt(nvar+1))
               nvar = nvar + nvariables
               if (lmodcur) then 
                 call init_modular_currents (nvariables, xc_opt(nvar+1))
                 nvar = nvar + nvariables
               end if
            end if
!
!     Initialize and count saddle coil variables (if lsaddle = true)
!
            if (lsaddle) then
               call init_saddle_coils (nvariables, xc_opt(nvar+1))
               nvar = nvar + nvariables
               if (lsadcur) then 
                  call init_saddle_currents (nvariables, xc_opt(nvar+1))
                  nvar = nvar + nvariables
               end if
            end if
!
!     Initialize and count background coil variables (if lbcoil = true)
!
            if (lbcoil) then
               if (myid .eq. master) then
                  call system("/bin/ln -s ../" // trim(bcoil_file) // 
     1                      " ./" // trim(bcoil_file))
               end if
#ifdef MPI_OPT

               call MPI_BARRIER (MPI_COMM_WORLD, ierr_mpi)
#endif
               call init_bg_currents (nvariables, xc_opt(nvar+1))
               nvar = nvar + nvariables
            end if
!
!     Initialize and count vf coil variables (if lvf = true)
!
            if (lvf) then
               call init_vf_currents (nvariables, xc_opt(nvar+1))
               nvar = nvar + nvariables
            end if

!
!     Load tf coil currents (if ltfc = true)
!
            if (ltfc .and. ltfcv) then
               call init_tf_coils (nvariables, xc_opt(nvar+1))
               nvar = nvar + nvariables
            end if
!
!     Initialize modular winding surface variables (if lsurfv = true)
!
            if (lsurfv) then
               call init_modular_wsurf (nvariables, xc_opt(nvar+1))
               nvar = nvar + nvariables
            end if
!
!     Initialize saddle winding surface variables (if lsadsfv = true)
!
            if (lsadsfv) then
               call init_saddle_wsurf (nvariables, xc_opt(nvar+1))
               nvar = nvar + nvariables
            end if

            lbnorm = sigma_berr_avg < bigno .or. sigma_berr_max < bigno

            if (myid .eq. master) print *, nvar-nvar_pre,
     1           ' coil shape coefs(II)'
         end if               !LCOIL_GEOM

      endif                   !LFREEB
      

!
!     COMPUTE NUMBER OF AI (NCURR=0) or AC (NCURR=1) COEFFICIENTS 
!     TO VARY FOR ACHIEVING OPTIMIZED PROFILES
!     IF NCURR = 0, AI VARIES. IF NCURR = 1, AC VARIES. 
!     THE TARGET_IOTA AND/OR TARGET_CURRENT FUNCTIONS WILL CONTRIBUTE 
!     TO CHISQ, DEPENDING ON THE RESPECTIVE SIGMAS.
!
      num_ai = 0
      if (lprof_opt) then                                                !Backwards compatability
         if (ncurr .eq. 1) lcur_prof_opt = .true.
         if (ncurr .eq. 0) liota_prof_opt = .true.
      end if   

      if (lbootsj_opt .and. (.not.liota_prof_opt) .and.
     1    (.not.lcur_prof_opt)) then
         if (myid .eq. master)                                           !MPI      
     1      print *,' LCUR_PROF_OPT was set = TRUE in initialize_opt'
         lcur_prof_opt = .true.
      end if     

      if (liota_prof_opt) then
         if (ncurr .ne. 0) then
            if (myid .eq. master) print *,' NCURR IS BEING SET TO 0',    !MPI      
     1      ' BECAUSE LIOTA_PROF_OPT = TRUE'
            ncurr = 0
         end if   
         lcur_prof_opt = .false.
!
!     VARY AI COEFFICIENTS
!      
         do ik = 0,10
            if (ai(ik) .ne. zero) num_ai = ik + 1
         end do  

         if (l_legendre) then                                             !LEGENDRE
           n_leg = num_ai-1                                               
           allocate(a_leg(0:n_leg,0:n_leg), b_leg(0:n_leg,0:n_leg),       
     1      a_leg_inv(0:n_leg,0:n_leg), b_leg_inv(0:n_leg,0:n_leg),       
     2      ti(0:n_leg))                                                  
          
           call build_matrices_legendre(n_leg, a_leg, b_leg,              
     1       a_leg_inv, b_leg_inv)                                        
           call power_to_legendre(n_leg, a_leg, b_leg, ai(0), ti)         
         end if                                                            
                                        
         do ik = 1,num_ai
            nvar = nvar + 1
            if(l_legendre) then                                          !LEGENDRE
              xc_opt(nvar) = ti(ik-1)                                     
            else                                                          
              xc_opt(nvar) = ai(ik-1)
            endif                                                         
         end do
         if (myid .eq. master) print *, num_ai,
     1           ' iota profile coefs'

      else if (lcur_prof_opt) then     
         if (ncurr .ne. 1) then
            if (myid .eq. master) print *,' NCURR IS BEING SET TO 1',    !MPI      
     1      ' BECAUSE LCUR_PROF_OPT = TRUE'
            ncurr = 1
         end if   
         liota_prof_opt = .false.
         if (l_legendre) ac_form = 0                                 !LEGENDRE
!
!     VARY AC COEFFICIENTS
!      
         do ik = 0,10                     !!Find index of last non-zero ac
           if (ac(ik) .ne. zero) num_ai = ik + 1
         end do  

         if (num_ai .eq. 0) then
            num_ai = 2
            ac(0) = 1
            ac(1) =-1
         else if (num_ai .eq. 11 .and. lcur_opt_edge0) then
c           num_ai = 10
            num_ai = 9
         end if
!
!     COMPUTE EDGE (INTEGRATED) CURRENT = SUM(AC(I)/(I+1)) FOR NORMALIZATION
!     WILL MULTIPLY IN PROFIL1D ROUTINE - BY CURTOR - TO GET PHYSICAL CURRENT           
!     SET EDGE CURRENT DENSITY, SUM(AC(I)), TO ZERO TO AVOID SPIKES THERE
!
         sum1 = zero
         select case (ac_form)
         case(1)
            sum1 = ac(0) + 2._rprec * ac(1)/3
            do ik = 2,num_ai-1
               sum1 = sum1 + ac(ik)/ik
            end do
         case default
            do ik = 0,num_ai-1
               sum1 = sum1 + ac(ik)/(ik+1)
            end do
         end select

         if (sum1 .ne. zero)
     1       ac(0:num_ai-1) = (curtor/sum1) * ac(0:num_ai-1)             

         if (l_legendre) then                                             !LEGENDRE
           n_leg = num_ai-1                                               
           allocate(a_leg(0:n_leg,0:n_leg), b_leg(0:n_leg,0:n_leg),       
     1      a_leg_inv(0:n_leg,0:n_leg), b_leg_inv(0:n_leg,0:n_leg),       
     2      tc(0:n_leg))                                                  
          
           call build_matrices_legendre(n_leg, a_leg, b_leg,              
     1       a_leg_inv, b_leg_inv)                                        
           call power_to_legendre(n_leg, a_leg, b_leg, ac(0), tc)         
         end if                                                            

         if( num_ai >= 0) then
            nvar = nvar + 1

            if(l_legendre) then
               xc_opt(nvar) = tc(0)
            else
!
!    make the first element be curtor, so that the feedback on the profile
!    shape is decoupled from the magnitude
!
               xc_opt(nvar) = curtor   
            endif
         endif

         do ik = 1,num_ai-1
!!$         do ik = 0,num_ai-1

            nvar = nvar + 1
            if(l_legendre) then                                           
              xc_opt(nvar) = tc(ik)                                     
            else                                                          
              xc_opt(nvar) = ac(ik) 
            endif                                                         
         end do
         if (myid .eq. master) print *, num_ai,
     1           ' current profile coefs'
         
      else if (ncurr .eq. 0) then        !do NOT match to iota target (since ai is fixed!)
         sigma_iota = bigno
         sigma_iota_max = bigno
         sigma_iota_max_min = bigno
         sigma_iota_min = bigno
      end if                           !End lprof_opt test

!
!        ALLOW MAGNITUDE OF THE CURRENT (CURTOR) TO VARY (IF NCURR=0)
!        TO MATCH THE TARGET_MAXCURRENT VALUE
!
       if (ncurr.eq.1 .and. (.not.lcur_prof_opt) .and.
     1      (abs(sigma_maxcurrent).lt.bigno .or. 
     2       sigma_curtor < bigno) ) then
          nvar = nvar + 1
          xc_opt(nvar) = curtor

          if (myid .eq. master) print *, 1,
     1           ' total current'
       end if   


!
!     scale factor for beta match (multiply mass profile coefficients)
!
      if (.not.lpres_prof_opt .and. 
     1    (abs(sigma_beta).lt.bigno .or. abs(sigma_eplasma) < bigno) 
     2    .and. (am(0).gt.zero)) then
         nvar = nvar + 1
         xc_opt(nvar) = am(0)

          if (myid .eq. master) print *, 1,
     1           ' pressure scale factor (I)'
      end if   
      
!---------------------------------------------------------------------------------
!   CODE added by R.SANCHEZ (01/19/99) to allow the mass profile to vary when
!      tailoring pressure profile for ballooning stability. Notice that when
!      tailoring (LPRES_PROF_OPT=.T.), ONLY the mass profile is varied).
!      NOTE: AM(10) is forced to remain equal -SUM(AM(i=0,9)) to guarantee p_edge=0.
!      It is ASSUMED that this condition is satisfied for the initial equilibrium!!
!   modified by MCZ, Mar 2003, normalize to ease optimizer's job
!       only constrain am(10) if not matching experimental data
       
      if (lpres_prof_opt) then
        nvar_pre = nvar
        if(.not. l_legendre) then                                        !LEGENDRE
          if( am(0) == 0) stop "am(0) = 0, no central pressure!"

!          am10 = am(10) / am(0)
!          am0_9 = sum(am(0:9)) / am(0)

          do ik=1, min(8, pres_opt_nmax)           
            nvar = nvar + 1
            xc_opt(nvar) = am(ik) / am(0)
          enddo

          if( .not.(lpres_opt_edge0 .and. lpres_opt_edgegr0) .and.
     1         pres_opt_nmax >= 9) then
             nvar = nvar + 1
             xc_opt(nvar) = am(9) / am(0)
          endif

          if( .not.(lpres_opt_edge0 .or. lpres_opt_edgegr0) .and.
     1         pres_opt_nmax >= 10) then
             nvar = nvar + 1
             xc_opt(nvar) = am(10) / am(0)
          end if

!   normalizing factor is s-integral of pressure
          nvar = nvar + 1
          xc_opt(nvar) = sum(am(0:10) / 
     1                       (/ 1,2,3,4,5,6,7,8,9,10,11 /) )

        else                                                             !LEGENDRE
          n_leg = 10                                                     
          allocate(a_leg(0:n_leg,0:n_leg), b_leg(0:n_leg,0:n_leg),       
     1      a_leg_inv(0:n_leg,0:n_leg), b_leg_inv(0:n_leg,0:n_leg),      
     2      tm(0:n_leg))                                                 
          
          call build_matrices_legendre(n_leg, a_leg, b_leg,              
     1       a_leg_inv, b_leg_inv)                                       
          call power_to_legendre(n_leg, a_leg, b_leg, am(0), tm)         

          tm10 = tm(10)                                                  
          tm0_9 = zero                                                   
          do ik=0, 9                                                     
            nvar = nvar + 1                                              
            xc_opt(nvar) = tm(ik)                                        
            tm0_9 = tm0_9 + tm(ik)                                       
          enddo                                                          
        endif

        if (myid .eq. master) print *, nvar-nvar_pre,
     1           ' pressure profile coefs.'

      end if   

c---------------------------------------------------------------------
! M. Zarnstorff  Mar 2003
!
! if matching discrete pressure profile and total energy or beta, need to
! adjust scale factor
c---------------------------------------------------------------------

      if ((abs(sigma_beta).lt.bigno .or. abs(sigma_eplasma) < bigno) 
     2    .and. np_prof > 0 .and. any(sigma_p_prof < bigno)) then
         nvar = nvar + 1
         xc_opt(nvar) = factor_p_prof*am(0)

         if (myid .eq. master) print *, 1,
     1           ' factor_p_prof'
      end if   

!------------------------------------------------------------------------------
! M.Zarnstorff  July 2000
! if this is a free-boundary optimization and matching to a 3D Vacuum Vessel constraint,
! then vary the enclosed toroidal flux to allow overall size of plasma to change

      if( (lfreeb .and. (sigma_vv<bigno .or. sigma_vv_rms<bigno 
     1              .or. sigma_vv_max<bigno .or. sigma_bd<bigno
     2              .or. sigma_bd_rms<bigno .or. sigma_bd_max<bigno )) 
     3     .or. (sigma_rbtor < bigno .and. .not.lcoil_geom) ) then
         nvar = nvar + 1
         xc_opt(nvar) = phiedge

         if (myid .eq. master) print *, 1,' phiedge'
      endif
!-------------------------------------------------------------------------------
      if (nvar .gt. nxc) stop 'nvar>nxc --> xc_opt out of bounds'
  
!
!     COMPUTE NOPT (NO. OF OPTIMIZATION FUNCTIONALS)
!
      nopt = -1
      call load_target (xc_opt, fvec, opt_ext, nopt, 0,  ik, .true.)
 
      if (myid .eq. master) print 100, nvar, nopt                        !MPI            
  100 format(/' No. Independent Variables = ',i5/,
     1   ' No. Dependent Constraints = ',i5/)
 
      end subroutine initialize_opt


      subroutine run_optimizer(xc_opt, nopt, nvar, lwa, info)
      use kind_spec
      use optim, ONLY: home_dir, lone_step, lrestart
      use optim_params, ONLY: epsfcn, niter_opt, seq_ext,
     1   num_processors, num_levmar_params, nopt_alg
      use vparams, ONLY: zero
      use mpi_params                                                            ! MPI
      implicit none
C-----------------------------------------------
C     D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nopt, nvar, lwa
      real(rprec), dimension(nvar) :: xc_opt
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, parameter :: info_size = 33
      integer :: info, niter
      integer :: mode                                                        
      real(rprec), dimension(:), allocatable :: fvec, diag
      real(rprec) :: tol
      character*120, dimension(1:info_size) :: info_array
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      external lsfun1
C-----------------------------------------------
      data info_array/
     1     'error in read_wout_opt opening wout file',

     2     'error in load_physics_codes in SYSTEM call',
     
     3     'error opening indata file in write_indata',
     
     4     'error opening output file in lsfun',
     
     5     'error writing new input file in load_params',
     
     6     'vmec2000 executable file not found',
     
     7     'i/o error in clean_up routine',

     8     'error reading wout file in call to read_wout_opt',
     
     9     'allocation error in load_target',
     
     A     'Boozer array dimension mismatch in load_target',
     
     B     'nvar and nopt do not match in load_target',
     
     C     'i/o error opening cobra file in chisq_ballooning',
     
     D     'i/o error opening bootstrap file in chisq_bootsj',
     
     E     'i/o error in open_comm_files',
     
     F     'error in chk_rzmnb',
     
     G     'boozer transform module - xbooz_xform - not found',
     
     H     'could not locate executable in load_physics_codes',

     I     'error reading output file in chisq_jinvar subroutine',
     
     J     'error in external kink computation',
     
     K     'error running xdkes code in chisq_dkes subroutine',
     
     L     'error in vacuum vessel matching subroutine',
     
     M     'system call to XCOILGEOM failed in generate_mgrid',
     
     N     'coils data file was not produced by XCOILGEOM',
     
     O     'error opening EXTCUR file in generate_mgrid',
     
     P     'error opening COIL_TARGETS file in chisq_coilgeom',

     Q     'error reading boozmn file in call to read_boozer_file',

     R     'error opening NEO code input file neo_in',

     S     'trouble running EQ3D in chisq_orbit',

     T     'trouble running MKJMC in chisq_orbit',

     U     'trouble running ORBIT in chisq_orbit',

     V     'error opening ORBSUM in chisq_orbit',

     W     'error opening FT79JMC in chisq_dsubr',
     
     X     'error in chisq_vac_island'
     Z    /

      allocate (fvec(nopt), diag(nvar), stat = info)      
      if (info .ne. 0)stop 'Allocation error in STELLOPT run-optimizer!'

      tol = 1.e-6_dp
      mode = 1
      niter = niter_opt
      diag = 0
      if (lone_step) niter = 1
      if(NOPT_ALG .eq. 0) then
#ifdef MPI_OPT
        call lmdif1 (lsfun1, nopt, nvar, xc_opt, fvec, tol, epsfcn,     
     1     niter, diag, mode, info, lwa)
#else
        call lmdif1 (lsfun1, nopt, nvar, xc_opt, fvec, tol, epsfcn,
     1     niter, diag, mode, info, lwa, num_processors,
     1     num_levmar_params)
#endif
      else if(NOPT_ALG .eq. 1) then
         call ga_driver (lsfun1, nopt, nvar, xc_opt, fvec, tol, epsfcn,
     1   niter, num_processors, seq_ext, info, lwa, lrestart )

      else if(NOPT_ALG .eq. 2) then
         call de_driver (lsfun1, nopt, nvar, xc_opt, fvec, tol, epsfcn,
     1   niter, num_processors, seq_ext, info, lwa, lrestart )

      else
         write(6,*) "NOPT_ALG undefined, unable to proceed"
         stop
      endif
      deallocate (fvec, diag)
      
      if (myid.eq.master .and. info.lt.0) write (*, '(/,1x,a,a)')
     1     'Stellopt status: ',trim(info_array(-info))
 
      end subroutine run_optimizer
      

      subroutine lsfun1(nopt, nvar, xc_opt, fvec, iflag, niter)
      use optim
      use safe_open_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nopt, nvar, iflag, niter
      real(rprec), dimension(nvar), intent(in) :: xc_opt
      real(rprec), dimension(nopt), intent(inout) :: fvec
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: iclean = -100
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat, iunit, i, ncnt
      character*200 :: temp_input
      logical :: ex, lscreen
C----------------------------------------------- 
!     ON ENTRY:
!        NITER:        NUMBER OF LEVENBERG ITERATIONS (=0 INITIALLY)
!        FVEC:         CONTAINS NOTHING THAT SHOULD BE USED
!        XC_OPT:       CONTAINS PERTURBED X VALUES
!        IFLAG = -1,   THIS IS A SINGLE-TASKED CALL, PRINT TO SCREEN
!        IFLAG 1-NOPT  MULTI-TASKED CALL, DO NOT PRINT TO SCREEN
!        IFLAG = -100, SPECIAL CLEAN UP CODE CALLED     
!
!     ON EXIT:
!        FVEC:         CONTAINS F(X)
!        XC_OPT:       UNCHANGED
!        IFLAG = 0     FINISHED WITHOUT A PROBLEM DETECTED
!        IFLAG < 0     ERROR SENT TO CALLER
!
!     SUGGESTIONS FOR IMPROVING CONVERGENCE
!     1. USE A REASONABLY SMALL VALUE FOR EPSFCN (3.E-4)
!     2. TRY 'LOOSENING' THE TOLERANCE IN FTOL ARRAY, i.e., RAISE
!        FTOL TO 1-5*E-9.
!     3. TRY 'TIGHTENING' THE TOLERANCE, i.e., FTOL <= 5.E-10
!     4. TRY PERTURBING R(m=1,n=0) and/or Z(m=1,n=0) A LITTLE
!
      if (iflag == iclean) then
         call clean_up (nvar, .false., iflag)
         return
      endif

      fvec(:nopt) = 0

!
!     COMPUTE UNIQUE EXTENSION FOR PARALLELIZED OPERATION
!
      istat = iflag
      if (iflag .eq. -1) istat = 0
      write (temp_input,'(i5)') istat
      opt_ext = trim(seq_ext)//'_opt'//trim(adjustl(temp_input))
      
      input_file = 'input.' // trim(opt_ext)
          
      if (iflag .ge. 0) then 
         ncnt = niter + iflag
         lscreen = .false.
      else
         ncnt = niter
         lscreen = .true.
      end if      
         
!
!     RUN VMEC TO COMPUTE TARGET CRITERIA.
!
 
      call load_params (xc_opt, nvar, iflag, lscreen)      

!
!     LOAD USER-DEFINED TARGET FUNCTIONALS
!
      if (ierr_vmec.eq.0 .and. iflag.eq.0) then
         iunit_opt_local = iout + 1
         call safe_open(iunit_opt_local, istat, 'output.' // 
     1                  trim(opt_ext), 'replace', 'formatted')
         if (istat .ne. 0) then
            iflag = -4
            return
         end if                  
         call load_target (xc_opt, fvec, opt_ext, nopt, ncnt, 
     1        iflag, lscreen)
         close (iunit_opt_local) 
      end if
!
!     IFLAG MAY HAVE FAILED IN LOAD_TARGET CALL, MUST CHECK IFLAG AGAIN
!
      if (ierr_vmec.ne.0 .or. iflag.ne.0) then

!        IF VMEC CONVERGENGE CONDITION VIOLATED, DISCOURAGE THIS
!        SEARCH DIRECTION BY IMPOSING A "BARRIER"
!        IF LSCREEN (FIRST CALL), THEN STOP: VMEC COULD NOT GET STARTED
!
         if( nopt_alg == 0) then
            fvec(:nopt) = 10*sqrt(chisq_min/nopt)
         else 
            fvec(:nopt) = 10*sqrt(bigno/nopt)
         endif
         if (ierr_vmec.ne.0 .and. .not.lscreen) iflag = 0
      end if
 

      end subroutine lsfun1


      subroutine clean_up (nvar, lastgo, iflag)
      use optim
      use optim_params, ONLY: num_levmar_params      
      use vmec_input, ONLY: raxis, zaxis, mgrid_file
      use safe_open_mod
      use system_mod
      use mpi_params   
      use gade_mod, only: npopsiz
      implicit none
#ifdef MPI_OPT
      include 'mpif.h'                                                   !MPI
#endif
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: nvar, iflag
      logical :: lastgo
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: iunit, icount, istat, istat2, itemp, icmax, 
     1           iteration, icount2
      real(rprec) :: aspect, chisq_tot
      logical :: lexist, lnew_min
      character*200 :: temp_input, save_min_ext, testfile
      character*25 :: ch1, ch2, ch3
C-----------------------------------------------
!
!     NOTE: BE SURE TO CATENATE -iflag output here, too 
!    (i starts at 0 for that case in do loop below)
!
      iunit_opt = iout + 3
      iflag = 0
      save_min_ext = " "
      lnew_min = .false.
      
      if (myid .ne. master) go to 2000
      
      inquire (file=trim(output_file), exist=lexist)
      if (lexist) then
         open(unit=iunit_opt, file=trim(output_file), status='old', 
     1        position='append', action='readwrite', iostat=istat)
      else
         open(unit=iunit_opt, file=trim(output_file), status='new', 
     1        iostat=istat)
      end if      

      if (istat .ne. 0) then
         print *,'Problem opening outputfile'
         call vmec_flush(6)

         iflag = -7
         go to 1999   
      end if   
         
      chisq_tot = chisq_min
      if( nopt_alg == 0) then
         icmax = max(nvar, num_levmar_params,numprocs)
      else if (nopt_alg > 0) then
         icmax = npopsiz
      endif

      findmin: do icount = 0, icmax
         if (lniter1 .and. icount.ne.0) exit
         write (temp_input,'(i5)') icount
         opt_ext = 
     1      trim(seq_ext)//'_opt'//trim(adjustl(temp_input))
         temp_input = 'output.' // trim(opt_ext)
         if (lniter1) save_min_ext = trim(opt_ext)

         iunit = iout+4
         call safe_open(iunit, istat, temp_input, 'old', 
     1             'formatted')
    
         do while (istat .eq. 0)
            read  (iunit,'(a)', iostat=istat) temp_input
            write (iunit_opt, '(a)') trim(temp_input)
            if (index(temp_input,'Aspect Ratio =').ne.0 .and.
     1          index(temp_input,'Chi-Sq =').ne.0) exit
         end do

         if (istat .eq. 0) then
            read  (temp_input,'(a16,f10.3,a10,f10.3,a23,i6)',
     1             iostat=istat)
     2             ch1, aspect, ch2, chisq_tot, ch3, iteration

!           Save in output file, delete temp file
            do 
               read  (iunit, '(a)', iostat=istat) temp_input
               if (istat .ne. 0) exit
               write (iunit_opt,'(a)') trim(temp_input)
            end do   

!
!        STORE PARAMETERS CORRESPONDING TO MINIMUM CHI-SQ STATE 
!        IN INPUT.OPT_EXT.MIN FILE. ALSO STORE WOUT.OPT_EXT.MIN FOR 
!        LRESET_OPT=F FAST RESTART OPTION. SAVE ITERATION IN ITER_MIN
!

            if (chisq_tot .lt. chisq_min) then
               chisq_min = chisq_tot
               iter_min = iteration   
               save_min_ext = trim(opt_ext)
               lnew_min = .true.
            end if

         end if

         close (iunit, iostat=istat)
 
      end do findmin      

!
!     STORE MIN STATE FILES, BUT ONLY IF ONE WAS FOUND
!
      if (lnew_min) then
         min_count = min_count + 1
         ch1 = ' '

!
!     Save the input.*_opt*  file ->  input.*.min
! 
         temp_input = "/bin/cp  input."//trim(save_min_ext) // " " //
     1                trim(min_input_file)
         call system(temp_input,istat)

!
!     If keeping all the .min files, store this one in sequence
!
         if(lkeep_mins) then
            ch1 = ' '
            write(ch1,'(".",i3.3)') min_count
            temp_input = "/bin/cp  input."//trim(save_min_ext) // " " //
     1                   trim(min_input_file) // trim(ch1)
            call system(temp_input,istat)
         endif

!
!     Construct the .min0 file containing the new boundary and modified axis
!
         temp_input = 'input.' // trim(save_min_ext)
         call read_input_opt (temp_input, istat2)    !!opens, deletes min_infile and loads modules
      
         if (istat2 .eq. 0) then
!           read in data (axis) corresponding to this wout file prior to writing out
            call read_wout_opt(.true., save_min_ext, istat, istat2)

            call write_indata (trim(min_input_file) // "_fb", iflag)
            if (iflag .ne. 0) goto 1999 
            
!           Blend old and new axis positions
!           This should be optional (SPH-01/10/01). Does NOT work well to update(blend) axis
!           UNLESS low enough ftol values are used...otherwise keep this option OFF
c           raxis_old = (3*raxis_old + raxis(:,1))/4
c           zaxis_old = (3*zaxis_old + zaxis(:,1))/4
            
         else if (save_min_ext(1:1) .ne. " ") then
            iflag = -7
            print *,'I/O Error reading file ', trim(temp_input),
     1              ' in STELLOPT routine CLEAN_UP, istat = ', istat2
            call vmec_flush(6)
            go to 1999
         end if

         if( save_min_ext(1:1) .ne. " ") then
!           Save (with min_input_file extension) all files for this minimum state 
!

            temp_input = "/bin/ls -1 *." // trim(save_min_ext) 
     1                // " > min_file"
            call system(temp_input, istat)


            itemp = iout+5
            call safe_open(itemp, istat, "min_file", "old", "formatted")
            do while (.true.) 
               read (itemp, '(a)', end=50) testfile
               icount = index(testfile, trim(save_min_ext)) - 1

               temp_input = "/bin/mv -f " // trim(testfile) 
     1                   // " " // testfile(1:icount) // trim(min_ext)
               call system(temp_input, istat)
            end do
 50         continue         
            close (itemp, status='delete')   
         end if
      end if              ! lnew_min
      

!     Delete input, wout files that are NOT minimum state (i.e, without .min extension)
!     Keep output file record and mgrid file
!
      if (.not.lniter1) then

!
!     make sure not to exceed the filename argument count for rm
!
         do icount = 10, icmax, 10
            write (temp_input,'(i5)') icount/10

            temp_input = 'rm -f *_opt'//trim(adjustl(temp_input))//'? &'
            call system(temp_input, istat)
         enddo
         call system('rm -f *_opt?', istat)
!         call system('rm -f fort.1???', istat)


         if( lastgo) then
!
!     pick up any stray files
!
            temp_input = "/bin/ls -1 > all_files"
            call system(temp_input)

            itemp = iout+5
            call safe_open(itemp, istat, "all_files", "old", 
     1                    "formatted")
            do while (.true.) 
               read (itemp, '(a)', end=60) testfile
!      Avoid wiping out necessary files for continuing run
            if ((nopt_alg.eq.1 .and. index(testfile,"ga_restart").gt.0)
     1     .or. (nopt_alg.eq.2 .and. index(testfile,"de_restart").gt.0)
     2     .or. (index(testfile, trim(output_file)) .gt. 0) 
     3     .or. (index(testfile, trim(mgrid_file)) .gt. 0) 
     4     .or. (index(testfile,"param") .gt. 0) 
     5     .or. (trim(testfile) .eq. "all_files")
     6     .or. (index(testfile,".dat") .gt. 0)
     5     .or. (index(testfile,"ft5") .gt. 0)
     6     .or. (index(testfile,"diagno.control") .gt. 0)) cycle
 
               icount = index(trim(testfile), trim(min_ext), 
     1                  back=.true.)
               if (icount .gt. 0) cycle
               temp_input = "rm -f " // trim(testfile)
               call system(temp_input, istat)
            end do

 60         close (itemp, status='delete', iostat=istat)

         end if     ! lastgo
      end if

 1999 close (iunit_opt, iostat=istat)     

 2000 continue


#ifdef MPI_OPT 
      call MPI_BCAST(chisq_min, 1, MPI_REAL8, master, MPI_COMM_WORLD, 
     1     istat)
      if (istat .ne. 0) stop 'MPI_BCAST error in STELLOPT CLEAN_UP'
      call MPI_BCAST(iflag, 1, MPI_INTEGER, master, MPI_COMM_WORLD,
     1     istat)
      if (istat .ne. 0) stop 'MPI_BCAST error in STELLOPT CLEAN_UP'
#endif

c    read the minimum wout file to get in the neighborhood
      call read_wout_opt(.true., min_ext, istat, istat2)

      end subroutine clean_up

      

      subroutine clean_up_extensions (extension, minext)
      use system_mod
      use safe_open_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      character*(*) :: extension, minext
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: iunit=20, istat, indx
      character*300 :: line, file_name, temp
C-----------------------------------------------
!
!     AT END OF RUN, REPLACES ALL THE WRONG FILE EXTENSIONS WITH THE .MIN EXTENSION
!     FIRST FIND (POSSIBLY) INCORRECT FILE EXTENSION IN THREED1 FILE
!
      file_name = 'threed1.' // trim(extension)
      call safe_open(iunit, istat, file_name, 'old', 'formatted')
      if (istat .ne. 0) return

      indx = 0

      do while (istat.eq.0)
         read (iunit, '(a)', iostat=istat, end=100) line
         indx = index(line, 'SHOT ID')
         if (indx .ne. 0) exit
      end do

 100  continue
      close (unit=iunit)
      if (indx .eq. 0) return

      temp = line(indx+10:)
      indx = index(temp, '_opt', back=.true.)
      if (indx .eq. 0) return

      indx = indx - 1
      line = temp(indx:)
      istat = index(temp,' ')
      if (istat .ne. 0) line = temp(indx:istat-1)

!     Make list of (text) file names to check and replace _opt with .min extension
!     Skip boozmn file, which is a binary file

      temp = "/bin/ls -1 *." // trim(extension) // " > ls_file"
      call system(temp, istat)
      call safe_open(iunit, istat, "ls_file", "old", "formatted")
      do while (istat.eq.0 .or. istat.ge.127)
         read (iunit, '(a)', end=150) file_name
!        skip binary files
         if (index(file_name, "boozmn.") .ne. 0) cycle
         temp = '/bin/sed -e "s/' // trim(line) // '/'
     1          // line(1:1) // trim(minext) // '/g" '
     2          // trim(file_name) // ' > tempxyz; /bin/mv -f tempxyz '
     3          // trim(file_name)
         call system(temp, istat)
      end do
 150  continue

      close (iunit, status='delete')

      end subroutine clean_up_extensions
      

      subroutine data_harvest (extension, output_file, nvar, nopt)
      use kind_spec
      use safe_open_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in)       :: nvar, nopt
      character*(*), intent(in) :: extension, output_file
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      character*(*), parameter :: file_list(5) =  (/
     1   "output  ", "threed1 ", "nescin  ", "nescout ", "dkes_opt" /)
      character*(*), parameter :: header(5) = (/
     1   "OPTIMIZATION PARAMETERS      ",
     2   "GEOMETRIC/MAGNETIC PARAMETERS",
     3   "COIL (NESCOIL) PARAMETERS    ",
     4   "COIL (NESCOIL) PARAMETERS    ",
     5   "TRANSPORT (DKES) PARAMETERS  "/)
      character*(*), parameter :: fmt(5) = (/
     1   "(a)   ", "(a)   ","(1x,a)", "(1x,a)", "(a)   " /)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: icount, iread=8, istat, iscratch, iwrite, jndex, 
     1      lcount, istat2
      real (rprec) :: chisq_min, chisq, dum, iota(3), es(3), dum1, dum2
      logical :: lwrite, lrewind, lheader
      character*200 :: line, filename
C-----------------------------------------------
!
!     CONTATENATES INFORMATION FROM SEVERAL DIFFERENT FILES INTO A
!     SINGLE "DATA_SUMMARY" FILE
!
      chisq_min = huge (chisq_min)
      lwrite = .false.;   lrewind = .false.
      es = -1

      iscratch = iread+1
      jndex = index(extension,".min",BACK=.TRUE.)
      if (jndex .eq. 0) jndex = len_trim(extension)+1
      call safe_open (iscratch, istat, "NULL", 'scratch', 'formatted')
      if (istat .ne. 0) return
      iwrite = iscratch+1
      call safe_open (iwrite, istat, "data_summary" // "." //
     1      trim(extension(1:jndex-1)), 'replace', 'formatted')
      if (istat .ne. 0) return


      FILE_LOOP: do icount = 1, size(file_list)
         if (icount .eq. 1) then
            filename = trim(output_file)
         else
            filename = trim(file_list(icount)) //"."// trim(extension)
         end if
         call safe_open(iread, istat, filename, 'old', 'formatted')
         if (istat .ne. 0) cycle

         lcount = 0
         if (icount .ne. 4) lheader = .true.

         PARSE_LOOP: do while (istat .eq. 0)
            read (iread, '(a)', end=200) line

!              parse line, depending on which one it is

            select case (icount)
            case (1)
               jndex = index (line, "Chi-Sq =")
               if (jndex .ne. 0) then
                  lrewind = .true.
                  read (line(jndex+8:), *) chisq
                  lwrite = (chisq .lt. chisq_min)
                  if (lwrite) chisq_min = chisq
               end if

               case (2)
               jndex = index (line, "S      <RADIAL    TOROIDAL")
               if (jndex .ne. 0) then
                  do lcount = 1,3
                     read (iread, '(a)') line
                  end do
                  lcount = 1
                  istat2 = 0
                  do while (istat2 .eq. 0)
                    read (iread, '(a)') line
                    read (line, *, iostat=istat2) dum1, dum, dum, dum2
                    if ((lcount.eq.1 .and. dum1.ge.0._dp) .or.
     1                 (lcount.eq.2 .and. dum1.ge.0.5_dp) .or.
     2                 (lcount.eq.3 .and. dum1.ge.0.9999_dp)) then
                       es(lcount) = dum1
                       iota(lcount) = dum2
                       lcount = lcount+1
                    else if (lcount .gt. 3) then
                       exit
                    end if
                  end do
                  lcount = 0
               end if
               jndex = index (line, "Aspect Ratio")
               if (jndex .ne. 0) then
                  lwrite = .true.; lrewind = .true.
               end if
               if (lwrite) lcount = lcount+1
               if (lcount .gt. 18) lwrite = .false.

               case (3)
               jndex = index (line, "Coil-Plasma separation =")
               if (jndex .ne. 0) then
                  read (line(jndex+24:), *) chisq
                  write(line,'(a,1pe10.3)') 
     1                  "Coil-Plasma Separation = ", chisq
                  lwrite = .true.;   lrewind = .true.
               end if
               if (lwrite) lcount = lcount+1
               if (lcount .gt. 1) lwrite = .false.
     
               case (4)
               jndex = index (line, "Complexity")
               if (jndex .ne. 0) lwrite = .true.
               if (lwrite) lcount = lcount+1
               if (lcount .gt. 8) lwrite = .false.

               case (5)
               lwrite = .true.
               lcount = lcount+1
               if (lcount .eq. 1) lrewind = .true.

            end select

            if (lwrite) then
               if (lrewind) then
                  rewind (iscratch, iostat=istat)
                  lrewind = .false.
               end if
               if ((icount .le. 3) .or.
     1             (icount .eq. 4 .and. (lcount==1 .or. lcount==4 .or.
     2              lcount==5 .or. lcount==8)) .or.
     3             (icount .eq. 5 .and. (lcount>2)))
     4             write (iscratch, fmt(icount)) trim(line)
            end if

         end do PARSE_LOOP

 200     continue

         close (iread)
         lwrite = .false.
         rewind (iscratch, iostat=istat)
         if (istat .eq. 0) then
            if (lheader) write (iwrite, '(1x,a)') header(icount)
!
!     WRITE NO. VARIABLES, PARAMETERS TO HEAD OF FILE
!
            if (icount .eq. 1)
     1      write(iwrite, '(/,a,i4,a,i6,/)')' Number of parameters = ',
     2      nvar,' Number of constraints = ',nopt 

            if (icount .eq. 3) lheader = .false.
            do while (istat .eq. 0)
               read (iscratch, '(a)', iostat=istat, end=300) line
               if (icount .eq. 1 .and. len_trim(line) .eq. 0) exit
               write (iwrite, '(a)') trim(line)
            end do
 300     continue
         if (icount .eq. 2) then
            write (iwrite,'(a,2(f3.2,a),f4.2,a,1p3e10.3)')
     1       ' iota(s=',es(1), ',', es(2), ',', es(3), ')  =      ',
     2       iota(1), iota(2), iota(3)
         end if
         if (icount .ne. 3) write (iwrite, *)
         rewind (iscratch)
         end if

      end do FILE_LOOP


      close (iscratch)
      close (iwrite)

      end subroutine data_harvest

      
      subroutine load_params(xc_opt, nvar, iflag, lscreen)
      use optim
      use bootsj_input
      use optim_params, ONLY: optimum
      use safe_open_mod
      use legendre_params
      use vmec_input
      use vparams, ONLY: zero
      use coilsnamin, only: lmodular, lsaddle, lbcoil, lvf, lsurfv,
     1    lsadsfv, ltfc, ltfcv, lmodcur, lsadcur
      use mpi_params                                                     !MPI
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(inout) :: iflag
      integer, intent(in) :: nvar
      real(rprec), dimension(nvar), intent(in) :: xc_opt
      logical :: lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: p96 = 0.96_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nb, ik, nvar_in, ireset, j1, istat, iunit
      integer :: nvariables, nzero
      real(rprec) :: sum0, sum1, newmass, hs
      real(rprec) :: am0_9_opt                                           !COBRA
      real(rprec) :: tm0_9_opt                                           !LEGENDRE
      character(len=len_trim(home_dir)+20) :: version
      character(len=len_trim(min_wout_file)+20) :: temp
      character :: screen*8
      logical :: lreset0, ex
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      real(rprec) , external :: eval_prof
      external unique_boundary, unique_boundary_PG
C-----------------------------------------------

      nvar_in = 0

!
!     SET UP RBC, ZBS BOUNDARY ARRAYS FROM ARRAY XC_OPT USED INTERNALLY
!     BY OPTIMIZATION CODE (THIS IS NOT THE VMEC XC ARRAY!)
!     NOTE: IF NITER_OPT = 1 AND THIS COMES FROM A FREE-BDY RUN, THEN
!     THE BOUNDARY MAY NOT BE IN THE EXPECTED FORM, SO KEEP RBC, ZBS 
!     UNCHANGED
!

      if (.not.lfreeb) then
         if (nopt_boundary .eq. 0) then
            do ik = 1, irm0_bdy
               j1 = ik + nvar_in
               rbc(nrz0_opt(ik),0) = xc_opt(j1)
            end do
            do ik = 1+irm0_bdy, irm0_bdy+izm0_bdy
               j1 = ik + nvar_in
               zbs(nrz0_opt(ik),0) = xc_opt(j1)
            end do      
            nvar_in = nvar_in + irm0_bdy + izm0_bdy
            do ik = 1, irho_bdy
               j1 = ik + nvar_in 
               rhobc(nbrho_opt(ik),mbrho_opt(ik)) = xc_opt(j1)
            end do
         else if (nopt_boundary .eq. 1) then
            do ik = 1, irho_bdy
              j1 = ik + nvar_in
              delta_mn(nbrho_opt(ik),mbrho_opt(ik)) = xc_opt(j1)
            end do
         endif

         nvar_in = nvar_in + irho_bdy

         if (.not.lniter1 .and. nopt_boundary.eq.0) 
     1       call unique_boundary(rbc, zbs, rhobc, ntord, mpol1d, 
     2                            mpol1_opt, ntor_opt)
         if (.not.lniter1 .and. nopt_boundary.eq.1) 
     1       call unique_boundary_PG(rbc, zbs, delta_mn, ntord, mpol1d, 
     2                            mpol1_opt, ntor_opt)

      else if (lfreeb) then                                              !Free-boundary optimization
         if( .not. lcoil_geom) then
            nb = 0
            do ik = 1, size(lextcur)
               if (lextcur(ik)) then
                  nb = nb + 1
                  extcur(ik) = xc_opt(nb)
               end if
            end do
            if (nb .ne. nextcur_opt) 
     1         stop 'Error counting external coils!'
            nvar_in = nvar_in + nextcur_opt

         else                                                            !lcoil_geom=T
!
!     Load and count modular coil variables
!
            if (lmodular) then
               call load_modular_coils (nvariables, xc_opt(nvar_in+1))
               nvar_in = nvar_in + nvariables
               if (lmodcur) then
                  call load_modular_currents
     1                 (nvariables, xc_opt(nvar_in+1))
                  nvar_in = nvar_in + nvariables
               end if
            end if
!
!     Load and count saddle coil variables (if lsaddle = true)
!
            if (lsaddle) then
               call load_saddle_coils (nvariables, xc_opt(nvar_in+1))
               nvar_in = nvar_in + nvariables
               if (lsadcur) then 
                  call load_saddle_currents 
     1                 (nvariables, xc_opt(nvar_in+1))
                  nvar_in = nvar_in + nvariables
               end if
            end if
!
!     Load and count background coil variables (if lbcoil = true)
!
            if (lbcoil) then
               call load_bg_currents (nvariables, xc_opt(nvar_in+1))
               nvar_in = nvar_in + nvariables
            end if
!
!     Load and count vf coil variables (if lvf = true)
!
            if (lvf) then
               call load_vf_currents (nvariables, xc_opt(nvar_in+1))
               nvar_in = nvar_in + nvariables
            end if

!
!     Load tf coil currents (if ltfc = true)
!
            if (ltfc .and. ltfcv) then
               call load_tf_coils (nvariables, xc_opt(nvar_in+1))
               nvar_in = nvar_in + nvariables
            end if
!
!     Load modular winding surface variables (if lsurfv = true)
!
            if (lsurfv) then
               call load_modular_wsurf (nvariables, xc_opt(nvar_in+1))
               nvar_in = nvar_in + nvariables
            end if
!
!     Load saddle winding surface variables (if lsadsfv = true)
!
            if (lsadsfv) then
               call load_saddle_wsurf (nvariables, xc_opt(nvar_in+1))
               nvar_in = nvar_in + nvariables
            end if

         end if               ! LCOIL_GEOM

      end if      ! LFREEB


!
!     LOAD AI or AC PROFILE COEFFICIENTS (IF LIOTA_PROF_OPT = TRUE or LCUR_PROF_OPT = TRUE)
!     PICK LAST COEFFICIENT SO THAT FOR NCURR=1, <JZETA(s=1)> = 0, TO AVOID
!     SURFACE CURRENT FORMATION
!

      if (num_ai .gt. 11) stop 'num_ai>11 in STELLOPT load_params'
      nzero = 10  ! num_ai    !  ac() entry to use to force j(a)=0.

      if (lcur_prof_opt) then
         nzero = 10  ! num_ai    !  ac() entry to use to force j(a)=0.

         if(l_legendre) then                                             !LEGENDRE
           tc(0:num_ai-1) = xc_opt(1+nvar_in:num_ai+nvar_in)              
           call legendre_to_power(num_ai-1, a_leg_inv, b_leg_inv,         
     1         tc, ac(0))                                                 
         else
           ac(0:10) = 0
           ac(0:num_ai-1) = xc_opt(1+nvar_in:num_ai+nvar_in)

!   first element is curtor
           sum0 = 0; sum1 = 0
           select case(ac_form)
           case(1)
              sum0 = ac(1)
              sum1 = 2._rprec * ac(1)/3
              do j1 = 2, num_ai-1
                 sum0 = sum0 + ac(j1)
                 sum1 = sum1 + ac(j1)/j1
              end do

           case default
              do j1 = 1, num_ai-1
                 sum0 = sum0 + ac(j1)
                 sum1 = sum1 + ac(j1)/(j1+1)
              end do
           end select

           if(.not. lcur_opt_edge0) then
              ac(0) = ac(0) - sum1          ! make sure the integral is right
             
           else
!
!     num_ai'th & 10th elements will be used to force j(edge)=0
!     without affecting total current
              ac(0) = ac(0) - sum1          ! make sure the integral is right

              select case(ac_form)
              case(1)
                 ac(num_ai) = num_ai*(ac(0)+sum0) / (10-num_ai)

              case default
                 ac(num_ai) = (num_ai+1)*(ac(0)+sum0) / (10-num_ai)

              end select
              
              ac(nzero) = -( ac(0) + sum0 + ac(num_ai))  ! force edge to zero

!!$              select case(ac_form)
!!$              case(1)
!!$                 ac(0) = ((nzero)/(nzero-1)) * (ac(0) -sum1)
!!$     1                 + sum0 / (nzero-1)
!!$              case default
!!$                 ac(0) = ((nzero+1)/nzero) * (ac(0) -sum1)
!!$     1                 + sum0 / nzero
!!$              end select
           endif
         endif                                                            

!
!    to keep edge current density == 0.
!
c         if( lcur_opt_edge0) ac(nzero) = -sum(ac(0:num_ai-1))

         sum1 = ac(0)
         select case(ac_form)
         case(1)
            sum1 = sum1 + 2._rprec * ac(1)/3
            do j1 = 2, 10
               sum1 = sum1 + ac(j1)/j1
            end do
         case default
            do j1 = 1, 10
               sum1 = sum1 + ac(j1)/(j1+1)
            end do
         end select

         curtor = sum1
         nvar_in = nvar_in + num_ai

      else if (liota_prof_opt) then
         if(l_legendre) then                                             !LEGENDRE
           ti(0:num_ai-1) = xc_opt(1+nvar_in:num_ai+nvar_in)              
           call legendre_to_power(num_ai-1, a_leg_inv, b_leg_inv,         
     1         ti, ai(0))                                                 
         else                                                             
           ai(0:num_ai-1) = xc_opt(1+nvar_in:num_ai+nvar_in)
         endif                                                            
              
         nvar_in = nvar_in + num_ai
      end if  
      
      if (ncurr.eq.1 .and. (.not.lcur_prof_opt) .and.
     1     (abs(sigma_maxcurrent).lt.bigno .or.
     2      sigma_curtor < bigno) ) then
         nvar_in = nvar_in + 1
         curtor = xc_opt(nvar_in)
      end if   
 
      if (.not.lpres_prof_opt .and. 
     1    (abs(sigma_beta).lt.bigno .or. abs(sigma_eplasma) < bigno)
     2    .and. (am(0).gt.zero)) then
         nvar_in = nvar_in + 1
         newmass = max(abs(xc_opt(nvar_in))/am(0), epsilon(am(0))) 
         am = newmass * am
      end if   

!-----------------------------------------------------------------------------------------------
!       CODE added by R.SANCHEZ (01/19/99).
!       COBRA: Vary mass coefficients to taylor pressure profile to ballooning modes.
 
      if (lpres_prof_opt) then
        if(.not. l_legendre) then                                        !LEGENDRE
          am0_9_opt = zero

          am(0) = 1
          am(1:10) = 0

          do j1 = 1, min(8, pres_opt_nmax)                  
            nvar_in = nvar_in + 1
            am(j1) = xc_opt(nvar_in)
          enddo

          if( lpres_opt_edge0 .and. lpres_opt_edgegr0) then

             am(9) = - sum( am(1:8)*
     1                    (/ 9,8,7,6,5,4,3,2 /) )      !use am(9 ) to guarantee p_edge gradient=0.  
             am0_9_opt = sum(am(0:9))  
             am(10) = - am0_9_opt                       !use am(10) to guarantee p_edge=0.  

          else if( lpres_opt_edge0 ) then
             if(pres_opt_nmax >= 9) then
                nvar_in = nvar_in + 1
                am(9) = xc_opt(nvar_in)
             endif

             am0_9_opt = sum(am(0:9))  
             am(10) = - am0_9_opt         !use am(10) to guarantee p_edge=0.  

          else if( lpres_opt_edgegr0 ) then 
             if(pres_opt_nmax >= 9) then
                nvar_in = nvar_in + 1
                am(9) = xc_opt(nvar_in)
             endif

             am(10) = - sum( am(1:9)*
     1                       (/ 1,2,3,4,5,6,7,8,9 /) ) / 10      !use am(10) to guarantee p_edge gradient=0.  

          else if(pres_opt_nmax >= 9) then
             nvar_in = nvar_in + 1
             am(9) = xc_opt(nvar_in)

             if(pres_opt_nmax >= 10) then
                nvar_in = nvar_in + 1
                am(10) = xc_opt(nvar_in)
             endif
          end if

!    normalize to s-integral of pressure profile
          nvar_in = nvar_in+1
          am(0:10) = am(0:10) * xc_opt(nvar_in) / 
     1              sum(am(0:10) / 
     2                  (/ 1,2,3,4,5,6,7,8,9,10,11 /) ) 


        else                                                             !LEGENDRE

          tm0_9_opt =zero                                                
          do j1 = 0, 9                                                   
            nvar_in = nvar_in + 1                                        
            tm(j1) = xc_opt(nvar_in)                                     
            tm0_9_opt = tm0_9_opt + tm(j1)                               
          enddo                                                          
          tm(10) = tm10 - (tm0_9_opt - tm0_9)                            
          call legendre_to_power(10, a_leg_inv, b_leg_inv, tm,           
     1        am(0))                                                     
        endif 
      endif

c---------------------------------------------------------------------
! M. Zarnstorff  Mar 2003
!
! if matching discrete pressure profile and total energy or beta, need to
! adjust scale factor
c---------------------------------------------------------------------

      if ((abs(sigma_beta).lt.bigno .or. abs(sigma_eplasma) < bigno) 
     2    .and. np_prof > 0 .and. any(sigma_p_prof < bigno)) then
         nvar_in = nvar_in + 1
         factor_p_prof = xc_opt(nvar_in) / am(0)
      end if   

!------------------------------------------------------------------------------
c M.Zarnstorff  July 2000
c if a free-boundary optimization and matching to a 3D Vacuum Vessel constraint
c then vary the enclosed toroidal flux to allow overall size of plasma to change

      if( (lfreeb .and. (sigma_vv<bigno .or. sigma_vv_rms<bigno
     1              .or. sigma_vv_max<bigno .or. sigma_bd<bigno
     2              .or. sigma_bd_rms<bigno .or. sigma_bd_max<bigno )) 
     3     .or. (sigma_rbtor < bigno .and. .not.lcoil_geom) ) then
         nvar_in = nvar_in + 1
         phiedge = xc_opt(nvar_in)
      endif
!-----------------------------------------------------------------------------------------------

      if (nvar_in .ne. nvar) then
         if (myid .eq. master)                                           !MPI      
     1      print *,' nvar_in = ', nvar_in,' nvar = ', nvar
         stop 'nvar_in!=nvar in STELLOPT load_params'
      end if   
      
      lreset0 = lreset_opt
      if (iflag .eq. -1) lreset0 = .true.
      ireset = 0

!************************************************************************
!
!     if varying pressure profile, test to make sure pressure is always
!     positive...

      if (lpres_prof_opt) then
         hs = 1._rprec / (nrad-1)
         do j1 = 2, nrad
            if( eval_prof(am, (j1-1.5_rprec)*hs ) < 0) then
               ierr_vmec = 1
               iflag = 0
               return       ! can't run vmec, so get out of here
            endif
         enddo
      endif

!************************************************************************
!     WRITE OUT DATA TO INPUT FILE (WILL BE READ BY XVMEC)
!
      lrecon = .false.
      if (lreset0) then
         raxis(:,1) = raxis_old
         zaxis(:,1) = zaxis_old
      end if        

      if (lcoil_geom) mgrid_file = 'mgrid.' // trim(opt_ext)

      call write_indata (input_file, iflag)
      if (iflag .ne. 0) return
      
!***********************************************************************
!     GENERATE MGRID FILE AND PARSE EXTERNAL CURRENTS ARRAY.
!     INPUT_FILE MUST BE CLOSED (ABOVE) TO BE READ IN GENERATE MGRID CALL
!     MUST RE-WRITE INPUT FILE TO REFLECT UPDATED CURRENTS.
!***********************************************************************
      if (lcoil_geom) then
         call generate_mgrid (opt_ext, (lscreen .and. (myid.eq.master)),
     1        iflag)
         if (iflag .ne. 0) return
        
!***********************************************************************
!        re-WRITE DATA TO INPUT FILE (WILL BE READ BY XVMEC), now that
!        we have extcur values. Also recomputes nextcur_vmec value.
!***********************************************************************
         call write_indata (input_file, iflag)
         if (iflag .ne. 0) return
      end if   

!********************************************************************
!     Execute vmec code and load vmec_input modules values
!***********************************************************************
      version = trim(home_dir) // '/xvmec2000'
      inquire(file=trim(version), exist=ex)
      if (ex) then
         screen = 'noscreen'
         if (lscreen) screen = 'screen'
         write (temp,'(1x,a,1x,l1)') screen, lreset0
         if (.not.lreset0) write
     1       (temp,'(a,1x,a)') trim(temp), trim(min_wout_file)
         iunit = 99
         call load_physics_codes (version, 'input', temp, 'wout',
     1        opt_ext, iunit, iflag)
      else
         iflag = -6
         return
      endif      

      if (lcoil_geom) call system("rm " // trim(mgrid_file))
 
      end subroutine load_params


      subroutine load_target(xc_opt, fvec, extension, nopt, ncnt, 
     1                       iflag, lscreen)
      use read_boozer_mod                                                !Loads bmn_b, etc.
      use optim
      use boozer_params
      use chisq_mod
      use safe_open_mod
      use vmec_input, ONLY: lfreeb, extcur, phiedge
      use vparams, ONLY: zero, one, twopi, dmu0
      use mpi_params                                       !MPI
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: nopt, iflag
      integer, intent(in) :: ncnt
      real(rprec), dimension(*) :: fvec, xc_opt
      character*(*), intent(in)  :: extension
      logical, intent(in) :: lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: p5 = 0.5_dp
      logical :: lprint_bmin = .false.
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: icontrol, jcount, no_elements, nsval, 
     1   n, k, i, nmodB, iunit, nweight(ntargets),
     2   mboz_surf, mskip
      real(rprec) :: chisq_tot, epsm
      real(rprec), dimension(ntargets) :: chisq, mean_weight
      real(rprec), dimension(:,:), allocatable :: bmod_b
      real(rprec), dimension(:), allocatable :: bmin_b, bmax_b
      real(rprec) :: hs, dmax_j, bmin_0, bmin_pi2, bmin_pi,  
     2  bmax_0, bmax_pi2, bmax_pi, dtheta, sqrt_nsurf
      character(len=len_trim(home_dir)+20) :: version
      logical :: ex, lscreen_master
C----------------------------------------------- 
!
!     DESCRIPTION OF LOCAL INPUT AND CONTROL VARIABLES
!
!
!  General Controls:
!--------------------
!       nopt_alg             Algorithm for optimization
!                            0 = Levenberg-Marquandt (default)
!                            1 = Genetic (GA)
!                            2 = Differential Evolution (DE)
!
!       lcoil_geom           =T, do coil geometry optimization 
!                            (coil harmonics are independent variables)
!                            Note that in this case, the initial coil shapes
!                            and the weights for the coil evaluations targets
!                            (e.g. coil curvatures and separations) are 
!                            specified in the COILSIN namelist.  
!                            See COILSOPT4.f90 for a description of this
!                            namelist.
!
!       epsfcn               'Annealing' parameter, 1.E-2 or less. If
!                            lreset_opt=F, epsfcn <= 1.E-3; if lreset_opt=T,
!                            epsfcn <= 1.E-4 is recommended
!
!       niter_opt            Maximum number of optimization iterations
!       num_processors       Maximum number of processors to try and use simultaneously
!       (nproc)                     
!
!       lreset_opt           =F, VMEC runs without resetting (after convergence) 
!                            to a coarse radial grid
!                            =T, VMEC resets each new iteration
!
!       mboz_opt             User-input number of poloidal boozer modes
!       nboz_opt             User-input number of toroidal boozer modes
!       nopt_boundary        internal representation used for optimization for 
!                            fixed-boundary runs (LFREEB=.f)
!                            0 = Hirshman-Breslau
!                            1 = Garabedian delta_mn
!       nsurf_mask           Array of reals (0 <= x <= 1) giving
!                            fraction of flux at which matching is to be done
!       lprof_opt            =F, the ai (iota, for ncurr=0, current, for ncurr=1) 
!                            profile expansion coefficients are NOT varied explicitly
!                            =T, the ai coefficient array is varied to optimize chi-sq
!                            (when ncurr=1, the TOTAL current is given by Target_MaxCurrent)
!       lpress_opt           =T  optimize pressure profile shape
!                            =F  keep pressure profile shape fixed
! 
!
!  Equilibrium & Geometry related
!------------------------------
!       Target_AspectRatio   Aspect Ratio value to match
!       Target_MaxCurrent    Maximum integrated toroidal current to match
!       Target_Beta          Volume averaged beta to match
!       Target_Eplasma       Total plasma energy to match
!       Target_Iota          Coefficients of Power series (in flux s)
!                            for iota profile to be matched
!       Target_Iota_P        Coefficients of Power series (in s) for iota-prime
!                               being d-iota/ds
!       Target_iota_min      minimum iota value (as lower bound)
!       Target_iota_max      maximum iota value (as upper bound)
!       Target_iota_max_min  maximum iota value (as lower bound)
!       Target_rmin          minimum major radius over the boundary
!       Target_rmax          maximum major radius over the boundary
!       Target_zmax          maximum height over the boundary
!       Target_ellipticity   desired elongation of phi=0 cross section
!       Target_kappa         desired <kappa>, calculated from n=0 terms
!       Target_fluxp         desired minimum poloidal flux (W)
!       sigma_aspect         Sigma for Asp.Ratio
!       sigma_MaxCurrent     Sigma for Maximum Integrated Current (used if ncurr=0)
!       sigma_beta           Sigma for volume-averaged beta
!       sigma_Eplasma        Sigma for total plasma energy
!       sigma_iota           Array (size nrad) of sigmas for iota, rel to avg iota
!       sigma_iota_pmax      Array (size nrad) of sigmas for target_iotap as
!                                  upper bound on d-iota/ds as
!       sigma_iota_pmin      Array (size nrad) of sigmas for d-iota/ds as
!                                  lower bound on d-iota/ds as
!       sigma_iota_min       sigma for iota below target_iota_min
!       sigma_iota_max       sigma for iota exceeding target_iota_max
!       sigma_iota_max_min   sigma for iota being below target_iota_max_min
!       sigma_curv           Sigma for forcing curve 'smoothness', *not relative*
!       sigma_rmin           sigma for forcing boundary > Target_rmin
!       sigma_rmax           Sigma for forcing boundary < Target_rmax
!       sigma_zmax           Sigma for forcing boundary height < target_zmax
!       sigma_kappa          sigma for <kappa>, calc. from n=0 terms
!       sigma_fluxp          Sigma for matching poloidal flux
!       laspect_max          =T, then sigma_aspect = bigno if Aspect Ratio <= Target_AspectRatio. 
!                            Used to guarantee A no larger than the target.
!       lbeta_min            =T, then sigma_beta = bigno if beta >= Target_Beta. Used to
!                            guarantee a minimum beta
!
!  Stability Related
!-------------------
!
!       Target_Well          Coefficients Power series (in flux s)
!                            for the magnetic well to be matched,
!                            [Vp(s) - Vp(0)]/Vp(0) [Replaced with Mercier criterion]
!       sigma_vp             Array (size nrad) of sigmas for well
!
!       sigma_Mercier        Sigma for approximate Mercier stability, *not relative*
!
!       lkink_opt            =T, do global kink stability calc.
!       sigma_kink           Array of Sigmas for full kink stability calc., 
!                            The first sigma value is used in a call to xtprp/ft5tpr
!                            the second value is used for xtprp2/ft5tpr2,
!                            and so-on for successive values.  This allows 
!                            targetting ofmultiple mode families.
!                            The sigma values are  *not relative*
!
!       target_kink          Array of target values for the kink eigenvalue
!                            this can be used to bias the chi-square to achieve
!                            stability
!
!       lbal_opt             =T, do ballooning calc. and add to chi-sq
!
!       lballoon_mask        1D array (dim=nrad) designating surfaces where 
!                            ballooning growth rates are to be computed (=T)
!       nballoon_mask        same function as lballoon_mask, but contains fractional
!                            s-values where ballooning growth rates are to be computed
!       lballoon_opt         If =T, ballooning evaluation is included in the optimization 
!                            If =F, no ballooning evaluation used in the optimization
!       lpres_prof_opt       If =T, ONLY the pressure profile is modified in the ballooning 
!                            optimization loop. When set to T, ALL surfaces must be included
!                            to be able to compute the pressure gradient.
!                            Must be =F to run optimizer in normal mode
!       lpres_opt_edge0      if T, and lpres_prof_opt=T, then adjust am(10) (and possibly am(9))
!                            to force the edge pressure to zero
!       lpres_opt_edgegr0    if T, and lpres_prof_opt=T, then adjust am(10) (and possibly am(9))
!                            to force the edge pressure gradient to zero
!       target_balloon       target ballooning eigenvalue (=0 for marginal stability)
!       lballoon_flip        If = F (default), nonzero penalty if eigenvalue is
!                                   more unstable than target_balloon
!                               = T  nonzero penalty if eigenvalue is more 
!                                    stable than target_balloon
!       sigma_balloon        1D array (dim=nrad) containing sigmas for ballooning stabilization
!       sigma_pgrad          1D array (dim=nrad) containing sigmas for local pressure gradient
!       sigma_pedge          1D array (dim=1) containing sigma for matching pressure value at edge
!
!
!  Coil Current - Related (for LFREEB=.T)
!----------------------------------------
!       lextcur              if lextcur(i) = T and lfreeb=True,
!                            extcur(i) will be varied during optimization
!                            Otherwise, extcur(i) is fixed during optimization
!       sigma_extcur         sigmas for weighting the coil-currents (free bdry)
!                            used for regularizing the coil currents
!       target_extcur        target values for the coil-currents
!
!       oh_coefs             array of coefficients for imposing a linear sum
!                            constraint on the coil currents, e.g. for ensuring
!                            that a PF set does not generate net poloidal flux
!                            constraint is sum_i (oh_coefs(i)*extcur(i))
!       sigma_oh             sigma for the OH flux constraint (weighted linear sum
!                            coil currents
!
!  Coil-Related
!-----------------------------------
!       sigma_berr_avg       sigma for the rms avg matching of Bnorm by Vmec
!       sigma_berr_max       sigma for the max. mismatch of Bnorm by Vmec
!
!
!  NESCOIL current-sheet coil estimate
!------------------------------------
!       lnescoil_opt         =T, evaluate NESCOIL coil optimization targets
!
!       Target_coil_complex  desired maximum coil complexity measure
!       Target_coil_jmax     desired maximum coil current density
!       sigma_coil_complex   sigma for the coil complexity measure (NESCOIL)
!       sigma_coil_jmax      sigma for the coil maximum current density (NESCOIL)
!       sigma_berr_ave       sigma for the average berr match (NESCOIL)
!
!   Transport optimization related
!---------------------------------
!
!       lfix_ntor            =F, rmn's, zmn's ARE NOT fixed in variation (free to vary)
!                            for the index 'n' of this array 
!                            =T, rmn's, zmn's ARE fixed (not varied) in optimization
!                            Useful if certain n-s are to be kept fixed (such as an
!                            axi-symmetric boundary: lfix_ntor(0) = T)
!       lbmn                 =T, add Boozer spectra to chi-sq (for QA, QH, QO depending on nbmn)
!
!       nbmn                 if lbmn=T then choose type of constrain on Boozer spectra
!                            =-1  suppress Bmn for m=0 or n != 0  (QO)
!                            = 0  suppress Bmn for n>0  (QA)
!                            = 1  suppress Bmn for m|=n (QH)
!       sigma_bmn            1D Array (nrad) of sigmas for matching BMNs
!                            relative to largest Bmn with programmed spatial weight
!       lprint_bmn (local var)  controls if boozer spectra is printed to a file (=T)
!                               only true if niter_opt = 1 and lbmn = T or 
!                               pseudo-sym is on.
!
!       sigma_jstar          2D Array (nrad by 4) of sigmas for Jstar
!                            at 4 values of ep/mu on each surface
!
!       ldkes_opt            =T optimize diffusion coefficients from DKES
!       sigma_dkes           =weights for DKES diffusion coefficients 
!
!       lneo_opt             logical, T= run NEO to calculate effective ripple
!                            default false
!       sigma_neo            array of sigmas for effective ripple for each surface
!
!       lorbit_opt           =T optimize fast ion confinement vis MC simulation
!       sigma_orbit          sigma for the fast ion loss
!
!       sigma_pseudo         Sigma for pseudo-symmetry 'water volume' of wells
!                            relative to toroidal well volume
!       sigma_pseudo2        Sigma for second 'water' measure: width of wells
!                            relative to field-line length
!       lpseudo_sin          weight pseudo water measures by abs(sin(theta))
!
!       sigma_bmin           1D Array (nrad) of sigmas for BMIN
!       sigma_bmax           1D Array (nrad) of sigmas for BMAX
!       sigma_ripple         1D Array (nrad) of sigmas for ripple
!
!
!   Bootstrap current related
!---------------------------
!       lbootsj_opt          =T, add match to bootstrap current to chi-sq
!
!       jboot                Bootstrap model to use:
!                               =0  Tolliver 3D, asymptotic collisionless
!                               =1  Simple 2D, asymptotic collisionless
!       sigma_bootsj         Sigma for matching bootstrap current self-consistently
!       sigma_boot           Sigma for bootstrap current density, rel to max of
!                            bootstrap or equil. current, w/prog. spatial weight
!
!       fboot                array of s**j coefficients for multiplicative factor
!                            on bootstrap current to mimic nu-star effects
!       zeff_boot            Zeff value for use in bootstrap current calculation
!       at                   polynomial coefficients for specifying plasma temperature
!                            Te=Ti for bootstrap, (keV)
!       lseedcur             =T, add seed current to bootstrap current before matching
!       aseedcur             polynomial coefficients for seed current, similar to AC coefs. 
!       lseed_opt            =F, the aseedcur coefficients are fixed
!                            =T, the non-zero aseedcur cofficients are optimized
!                            ***not yet implemented***
!
!
!    Flux surface quality related
!--------------------------------
!
!       n_jac                   array of n-number for resonantBoozer jacobean components
!       m_jac                   array of m-numbers for resonant Boozer jacobean components
!       sigma_jac               sigmas for supressing resonant Boozer jacobean components
!
!       n_vac_island            array of n-number for vacuum island targetting
!       m_vac_island            array of m-numbers for vacuum island targetting
!       sigma_vac_island        sigmas for supressing vacuum islands
!
!       ldsubr_opt              logical, T= run JMC to calculate Dr
!       sigma_dsubr             array of sigmas for resistive interchange (Dr) stability
!
!
!    Vacuum vessel and limiter matching related
!--------------------------------
!    ability to target two boundary shapes is provided (_vv and _bd)
!    normally _VV is used to target an enclosing PFC shape, and _BD is used to
!    target a desired plasma shape
!
!       target_vv               target for closest approach distance
!       lvv_tgt_min             if T, then there is no penalty if closest
!                               approach distance is larger than target_vv
!
!       target_vv_rms           target for RMS average distance
!       sigma_vv                sigma for closest distance to specified -VV
!                                  3D boundary shape (if <~1.e10, else ignore)
!                                  and limiters
!       sigma_vv_rms            sigma for RMS average distance from plasma to
!                                _VV 3D boundary (if <~ 1.e10, else ignore)
!       sigma_vv_max            sigma for Maximum deviation distance from 
!                                  plasma to _VV 3D boundary
!       mpol_vv                 maximum m for _VV 3D boundary
!       ntor_vv                 maximum n for _VV 3D boundary
!       rbc_vv                  R-cos(theta) array for _VV 3D boundary
!       zbs_vv                  Z-sin(theta) array for _VV 3D boundary 
!
!
!       target_bd               target for closest approach distance
!       target_bd_rms           target for RMS average distance
!       sigma_bd                sigma for closest distance to specified _BD
!                                  3D boundary shape (if <~1.e10, else ignore)
!                                  and limiters
!       sigma_bd_rms            sigma for RMS average distance from plasma to
!                                 _BD 3D boundary (if <~ 1.e10, else ignore)
!       sigma_bd_max            sigma for Maximum deviation distance from 
!                                  plasma to _BD 3D boundary
!       mpol_bd                 maximum m for _BD 3D boundary
!       ntor_bd                 maximum n for _BD 3D boundary
!       rbc_bd                  R-cos(theta) array for _BD 3D boundary
!       zbs_bd                  Z-sin(theta) array for _BD 3D boundary
!
!  Piece-wise linear limiter specification, used with _VV targets
!
!       phi_lim                 toroidal angles of limiters (1-D array)
!                               limiters are specified at each toroidal angle
!                               as an array of piecewise-linear segments, 
!                               going counter-clockwise around the plasma.
!                               For each phi_lim(iphi), the i-th segment
!                               goes from (r_lim(i,iphi),z_lim(i,iphi)) to
!                               (r_lim(i+1,iphi),zlim(i+1,iphi))
!       r_lim                   2D-array of major radial locations of limiter
!                               segments  (i, iphi)
!       z_lim                   2D-array of vertical locations of limiter
!                               segments  (i, iphi)
!
! -----------------------------------------------------------------------
!       target_RBtor            Target value for R-Btor at plasma edge
!       sigma_RBtor             Sigma for matching R-Btor
*
!       shapeweight             If F boundary weighting is turned off and
!                               the VV routine is used. Default = .T.
!       planes_bw(3) (radians)  planes in which shape is evaluated
!                               Do not normalize to Nfp
!                               0, Pi/2 and Pi are recommended
!       amplw_bw(3) (radians)   relative weighting of planes(3)
!       theta0_bw(3) (radians)  poloidal angle at which weight function maiimizes
!                               Gaussian Weigts: 15 terms
!                               theta0 centers are at theta0, theta0+Pi/2, 
!                               theta0+Pi, theta0+3Pi/2 &  theta0+2 Pi
!       phi0_bw (radians)       toroidal angles at which weights maximize
!                               centers at phi0/Nfp, (phi0+Pi/2/)Nfp, & (phi0+ Pi)/Nfp
!       wtheta_bw (radians)     poloidal half-width of shaping weights; e.g.,Pi/8
!       wphi_bw (radians)       toroidal half-width of shaping weights
!
!
!
!    Reconstruction related
!--------------------------------
!       ldiagno_opt             simulate & match magnetic diagnostics with diagno
!       diagno_control          name of diagno control file to use 
!       sigma_diagno_seg        array of sigmas for segmented coils
!       sigma_diagno_flx        array of sigmas for flux loops
!       target_diagno_seg       array of target values for segmented coils
!       target_diagno_flx       array of target values for flux loops
!
!       np_prof                 number of points to match in pressure profile
!       p_prof                  array of pressure values
!       r_p_prof                array of major radius values for pressure points (m)
!       z_p_prof                array of height values for pressure points (m)
!       phi_p_prof              array of toroidal angles for pressure points (degrees)
!       sigma_p_prof            array of sigmas for pressure points
!
!       target_curtor           target total current
!       sigma_curtor            sigma for total current
!
!       target_jedge            target for edge current density
!       sigma_jedge             sigma for edge current density
!
!-----------------------------------------------------------------------------------------------
!     Load information from wout file
!-----------------------------------------------------------------------------------------------
      if (nopt .le. 0) then
         call read_wout_opt(.false., extension, jcount, iflag)    !!loads boozer parameters 1st time thru
         if (iflag .ne. 0) then
            iflag = -1
         else if (jcount .ne. 0) then
            iflag = -8
         end if      
         if (iflag .ne. 0) return

      else
         dtheta = one
         if (nu2_b .gt. 1) dtheta = one/real(nu2_b - 1,rprec)
         allocate (index_array(nopt),
     1             wegt(nopt), chisq_match(nopt), chisq_target(nopt),
     2             chisq_descript(nopt), stat=n)

         if (n .ne. 0) then
            iflag = -9
            return
         end if   

         index_array = 0
         chisq_target = 0
         chisq_match = 0
         wegt = bigno

      endif

      allocate (bmin_b(nu2_b), bmax_b(nu2_b), bmod_b(nv_boz,nu_boz), 
     1         stat=n)
      
      if (nrad .gt. 1) hs = one/real((nrad - 1),rprec)
      nmodB = nu2_b
 
!-----------------------------------------------------------------------------------------------
!     COMPUTE BOOZER TRANSFORMATION HERE ON RADIAL SURFACES
!     DETERMINED BY LSURF_MASK LOGICAL ARRAY. IF NO RADIAL DERIVATIVES
!     OF Bmn REQUIRED, ONLY NEED ONE SURFACE AT A TIME (OTHERWISE, COMPUTE
!     AND STORE ADJACENT SURFACES).
!-----------------------------------------------------------------------------------------------

      mskip     = nu2_b/mboz
      mboz_surf = 0
      do n = 1, nu2_b, mskip
         mboz_surf = mboz_surf + 1
      end do

      sqrt_nsurf = sqrt(real(mboz_surf,rprec))
      write (comm1,'(1x,l1)') lscreen            !!Needed for command line arg
      version = trim(home_dir) // '/xbooz_xform'
      lscreen_master = lscreen .and. (myid.eq.master)

      if (nopt .gt. 0) then

!-----------------------------------------------------------------------------------------------
!     OPEN INPUT FILES NEEDED FOR COMMUNICATING WITH EXECUTABLES
!-----------------------------------------------------------------------------------------------
         call open_comm_files(mboz, nboz, ns_booz, ns_booz_max, 
     1      ns_surf, ns_surf_max, lbootsj_opt, extension, iflag)

!-----------------------------------------------------------------------------------------------
!        Write and read in bmn from boozmn.extension file (refer to MISCEL.F90 file)
!-----------------------------------------------------------------------------------------------
         call load_physics_codes (version, 'in_booz', comm1, 'boozmn', 
     1        extension, iunit, iflag)
         if (iflag .ne. 0) return
         
         if (mnboz .ne. mnboz_b .or.  mboz.ne.mboz_b .or. 
     1       nboz.ne.nboz_b) then
             iflag = -10
             return
         end if    
       
!-----------------------------------------------------------------------------------------------
!     open DKES_OPT file to accumulate OPT_DKES data for all surfaces
!-----------------------------------------------------------------------------------------------
         if (ldkes_opt) then
            iunit_dkes = 19
            call safe_open(iunit_dkes, k, 'dkes_opt.'//trim(extension),
     1                  'replace', 'formatted')
            if (k .ne. 0) then
               stop 'Error opening dkes_opt file in load_target'
            else
               write(iunit_dkes,'(2a,/)')
     1           'SUMMARY OF DKES TRANSPORT COEFFICIENTS FOR CASE: ',
     2           trim(extension)
               write(iunit_dkes,"(78('='),/,3x,3(9x,a,9x),/,78('='))")
     1           'L11(+/-)', 'L13(+/-)', 'L33(+/-)'
            end if
         end if
!-----------------------------------------------------------------------------------------------
!     initialize AJAX for inverse coordinate mapping, if needed
!-----------------------------------------------------------------------------------------------
         if( lajax) then
            call init_ajax
         endif

      else
         inquire(file=trim(version), exist=ex)
         if (.not.ex .and. myid.eq.master) then
            print *, 'xbooz_xform file not found in ' // trim(home_dir)
            iflag = -16
            return
         end if   
         allocate (bmn_b(1,nrad), stat = i)              !!Need for dummy arg calls to count nopt 1st time
         lajax = .false.                                 ! initialize flag on AJAX inverse mapping
      endif      

!-----------------------------------------------------------------------------------------------
!     BEGIN MATCHING TARGETS
!-----------------------------------------------------------------------------------------------
      no_elements = 0
      iflag = 0
!-----------------------------------------------------------------------------------------------
!
!     MATCH ASPECT RATIO
!
      call chisq_aspect (aspect_opt, target_aspectratio,
     1     sigma_aspect, no_elements, nopt)

!-----------------------------------------------------------------------------------------------
!
!     MATCH MAXIMUM TOROIDAL CURRENT 
!    (DO NOT USE ABSOLUTE VALUE: OPTIMIZER UNSTABLE THAT WAY)
! 
      call chisq_maxcurrent (target_maxcurrent, sigma_maxcurrent, 
     1                       no_elements, nrad, nopt)

!-----------------------------------------------------------------------------------------------
!
!     Match desired poloidal flux               MCZ  Sep. 98
!
      if (sigma_fluxp .lt. bigno) call chisq_polflux (target_fluxp,
     1    sigma_fluxp, hs, nrad, no_elements, nopt)      
      
!------------------------------------------------------------------------------
!
!     MATCH VOLUME-AVERAGED BETA
!
      call chisq_beta (Target_Beta, sigma_beta, 
     1                 Target_eplasma, sigma_eplasma,
     2                 no_elements, nopt) 

!------------------------------------------------------------------------------
!
!     MATCH Total current
!
      call chisq_curtor (Target_curtor, sigma_curtor, 
     1                   no_elements, nopt) 

!------------------------------------------------------------------------------
!
!     MATCH edge current density
!
      call chisq_jedge (Target_jedge, sigma_jedge, 
     1                   no_elements, nopt) 

!-----------------------------------------------------------------------------------------------
!
!     MATCH R-Btor AT S=1 (USED IN FREE-BDY MODE TO SET PHIEDGE)             MCZ  July 00
!
      call chisq_rbtor (rbtor_opt, target_rbtor, 
     1     sigma_rbtor, no_elements, nopt)

!-----------------------------------------------------------------------------------------------
!
!     MATCH IOTA PROFILE (IF NCURR=1) 
!
      call chisq_iota(iota_opt, sigma_iota, hs, no_elements, nrad, nopt)

! -----------------------------------------------------------------------
!     iota-prime                             M.Zarnstorff
!
     
      call chisq_iota_p(iota_opt, sigma_iota_pmax, sigma_iota_pmin, 
     1                 hs, no_elements, nrad, nopt)

!-----------------------------------------------------------------------------------------------
!
!     MERCIER CRITERION (DETERMINES STABLE MAGNETIC WELL = vp(s)-vp(0), normed to vp(0))
!     Dmerc < 0 is UNSTABLE
!
      call chisq_mercier (Dmerc_opt, sigma_mercier, no_elements, 
     1     nrad, nopt)

!-----------------------------------------------------------------------------------------------
!
!     EXPLICITLY DETERMINE MAGNETIC WELL = vp(s)-vp(0), normed to vp(0)
!
      if (any(abs(sigma_vp(2:nrad)) .lt. bigno)) call chisq_magwell 
     1   (hs, sigma_vp, no_elements, nrad, nopt)

!-----------------------------------------------------------------------------------------------
!     CODE added by R.SANCHEZ (01/19/99). Modified (02/01/99) to include
!     multiple initial toroidal and poloidal positions.
!-----------------------------------------------------------------------------------------------
 
      if (lballoon_opt) then
         call chisq_ballooning (sigma_balloon, target_balloon, 
     1       lballoon_flip, no_elements,  ns_ball, ns_ball_max,     !!VMEC COBRA (RS)  
     2       nopt, iflag, extension, comm1)
         if (iflag .ne. 0) return
      end if

      if (lpres_prof_opt) then
         call chisq_bpres (pres_opt, sigma_pedge, ivar_pedge, 
     1        no_elements, nrad, nopt, 0)
         call chisq_bpres (pres_opt, sigma_pgrad, ivar_pgrad, 
     1        no_elements, nrad, nopt, 1)

!-----------------------------
!    constrain pressure profile at discrete points, e.g. to match experimental data
!    MCZ 2002
!
         if(np_prof > 0 .and. any(sigma_p_prof<bigno))
     1      call chisq_p_prof(pres_opt, ivar_p_prof, no_elements, nrad,  
     2                        nopt, extension)
      end if

!-----------------------------------------------------------------------------------------------
!
!     CODE added by R. SANCHEZ, L. BERRY, and S. HIRSHMAN (01/19/99) to MATCH 
!     BOOTSTRAP CURRENT <J dot B> PROFILE COMPUTED FROM VMEC
 
      if(lbootsj_opt) then
         call chisq_bootsj (sigma_bootsj, no_elements, 
     1      ns_surf, ns_surf_max, nrad, nopt, iflag, extension, 
     2      lscreen_master)
         if (iflag .ne. 0) return
      end if


!----------------------------------------------------------------------------------------
!     DIAGNO targetting magnetic diagnostics
!-----------------------------------------------------------------------------------------

      if( ldiagno_opt .and. ndiagno_seg+ndiagno_flx > 0) then
          call chisq_diagno (no_elements, nopt, iflag, 
     1                            extension, lscreen_master)

          if (iflag .ne. 0) return
      endif


!-----------------------------------------------------------------------------------------------
!     EXTERNAL KINK MODE CRITERION (LPK, MCZ, PPPL - Feb 1999, SPH Feb 2000)
!-----------------------------------------------------------------------------------------------
      if (lkink_opt) call chisq_kink (no_elements, 
     1      nopt, iflag, extension, lscreen_master)
      if (iflag .ne. 0) return
 

!-----------------------------------------------------------------------------------------------
!
!     Coil targets from XCOILGEOM, as controled in the COILSIN namelist
!
      if (lcoil_geom .or. sigma_berr_avg < bigno .or.
     1    sigma_berr_max < bigno  ) 
     2       call chisq_coilgeom (no_elements, nopt, iflag, 
     3                            extension, lscreen_master)
      if (iflag .ne. 0) return

!-----------------------------------------------------------------------------------------------
!     COMPUTE TARGET CRITERIA THAT DEPEND ON BOOZER-TRANSFORMED QUANTITIES
!     SUCH AS: J*, Jinvariant, Bmin, Bmax, Ripple, Bmn
!-----------------------------------------------------------------------------------------------
      if (lprint_bmin) then
         iunit = 18
         call safe_open (iunit, k, 'bmin_bmax.' // trim(extension), 
     1                  'replace', 'formatted')
         write (iunit, *) ns_surf_max, nu2_b
      end if

      icontrol = 0
      SURF: do jcount = 1, ns_booz_max
         nsval = ns_booz(jcount)
!-----------------------------------------------------------------------------------------------
!        COMPUTE BOOZER-SPACE MOD-B vs theta-booz, zeta-booz IF NECESSARY
!-----------------------------------------------------------------------------------------------
         if (nopt .gt. 0) then
            if (jcount .eq. ns_booz_max) icontrol = -1    !!release memory here
            if (lneed_modB(nsval) .or. (icontrol .eq. -1)) then
               call modbooz(bmn_b(1,nsval), bmod_b, ixm_b, 
     1                      ixn_b, nfp_opt, icontrol)
               call bextrema (bmod_b, bmin_b, bmax_b, nv_boz, nu2_b)
            end if

            if (lprint_bmin) then
               write (iunit, *) nsval
               do k = 1,nu2_b
                    write (iunit, *) k, bmin_b(k), bmax_b(k)
               end do
               if (jcount .eq. ns_booz_max) close (iunit)
            end if
         end if                !!End nopt > 0

!-----------------------------------------------------------------------------------------------
!     added by R. FOWLER (03/03/00)
!-----------------------------------------------------------------------------------------------
         if (ldkes_opt .and. ldkes_mask(nsval)) then                                              
            call chisq_dkes(sigma_dkes(nsval), no_elements, 
     1          nopt, iflag, nsval, extension, comm1)                            
            if (iflag .ne. 0) return                                      
         end if                                                           

         if (.not.lsurf_mask(nsval)) cycle

!-----------------------------------------------------------------------------------------------
!     MATCH BMN-BOOZER CRITERION. DUE TO THE VARIETY OF POSSIBLE
!     WEIGHTS (SIGMAS) AND TARGET COMBINATIONS (BMN_B=0, BMN_B(m,n) = BMN_B(m1,n1), etc),
!     THE USER MAY DECIDE TO HAND-CODE THE WEIGHTS AND TARGETS IN THE FOLLOWING LOOP
!-----------------------------------------------------------------------------------------------
         if (lbmn) call chisq_bmn(bmn_b(1,nsval), sigma_bmn(nsval),
     1      hs, ixm_b, ixn_b, no_elements, mnboz, nsval, nopt)

!-----------------------------------------------------------------------------------------------
!   Pseudo-symmetry criteria; try and discourage local wells. Modified by MCZ, Feb. 99, SPH Feb. 00
!-----------------------------------------------------------------------------------------------
         if (sigma_pseudo.lt.bigno .or. sigma_pseudo2.lt.bigno) 
     1      call chisq_pseudo (sigma_pseudo, sigma_pseudo2, hs, 
     2         bmn_b(1,nsval), ixn_b, ixm_b, no_elements, nopt, nsval, 
     3         lpseudo_sin, lscreen_master)

!-----------------------------------------------------------------------------------------------
!        MINIMIZE THE VARIATION OF THE SECOND ADIABATIC INVARIANT, J-INVARIANT,
!        AROUND A FLUX SURFACE FOR SELECTED RADII AND VALUES
!        OF epsilon/mu (=Bturning). THIS CALLS xj_invariant WHICH COMPUTES |B|
!        ITSELF FROM THE BOOZMN FILE, SO IT DOES NOT REQUIRE A CALL HERE TO MODBOOZ
!        Added by D. A. Spong, 1/00
!-----------------------------------------------------------------------------------------------
         if (lj_invariant) call chisq_jinvar 
     1     (sigma_jinvariant(nsval,1:NumJinvariant), 
     2      nsval, no_elements, nopt, nu2_b/mskip, sqrt_nsurf, iflag, 
     3      extension, comm1, lscreen_master)
         if (iflag .ne. 0) return  

!-----------------------------------------------------------------------------------------------
!
!        I M P O R T A N T    D E V E L O P E R    N O T E
!
!        <<ALL>> ROUTINES HEREAFTER (IN SURF LOOP) REQUIRE BMOD_B, BMIN_B, AND/OR BMAX_B.
!        IF ANY ARE USED (FINITE SIGMA), THEN THE CALL TO MODBOOZ ABOVE MUST HAVE BEEN
!        MADE. THIS IS DETERMINED BY THE LOGICAL LNEED_MODB, WHICH IS COMPUTE THE FIRST TIME
!        THROUGH (FOR NOPT <= 0).
!
!-----------------------------------------------------------------------------------------------
         if (nopt .le. 0) then
            n = no_elements
         else if (.not.lneed_modB(nsval)) then
            cycle
         end if
         
!-----------------------------------------------------------------------------------------------
!        Minimize (Bmin - <Bmin>)
!-----------------------------------------------------------------------------------------------
         call chisq_bmin (sigma_bmin(nsval), ivar_bmin, bmin_b, 
     1          no_elements, nopt, mskip, sqrt_nsurf, dtheta)

!-----------------------------------------------------------------------------------------------
!        Minimize (Bmax - <Bmax>)
!-----------------------------------------------------------------------------------------------
         call chisq_bmin (sigma_bmax(nsval), ivar_bmax, bmax_b, 
     1          no_elements, nopt, mskip, sqrt_nsurf, dtheta)
 
!-----------------------------------------------------------------------------------------------
!        Minimize magnetic ripple on each surface
!        7/98: added sin(u/2) weight to emphasize outboard ripple reduction
!-----------------------------------------------------------------------------------------------
         call chisq_bripple (sigma_ripple(nsval), bmax_b,
     1          bmin_b, no_elements, nopt, mskip, sqrt_nsurf, dtheta)

!-----------------------------------------------------------------------------------------------
!        MINIMIZE THE THE VARIATION OF THE TRAPPED PARTICLE J*
!        AROUND A FLUX SURFACE FOR SELECTED RADII AND VALUES
!        OF epsilon/mu (=Bturning).
!-----------------------------------------------------------------------------------------------
         if (lj_star) call chisq_jstar (sigma_jstar(nsval,1:NumJstar), 
     1      bmax_b, bmin_b, bmod_b, no_elements, nopt, 
     2      mskip, sqrt_nsurf)
!-----------------------------------------------------------------------------------------------
         
         if (nopt.le.0 .and. no_elements.ne.n) lneed_modB(nsval)=.true.

      end do SURF
 
!
!     Clean up: close any open files, etc.
!
      if (ldkes_opt .and. nopt.gt.0) close(unit=iunit_dkes)

!-----------------------------------------------------------------------------------------------
!!
!!    NEO  LPK 01-19-01
!!
         if (lneo_opt) call chisq_neo(sigma_neo, ivar_neo, no_elements,
     1          ns_neo, ns_neo_max, nopt, iflag, extension, 
     2          lscreen_master)
         if (iflag .ne. 0) return  

         if (ldsubr_opt) call chisq_dsubr(sigma_dsubr, ivar_dsubr, 
     1          no_elements, nopt, iflag, extension, lscreen_master)
         if (iflag .ne. 0) return  

         if (lorbit_opt) call chisq_orbit(sigma_orbit, ivar_orbit,
     1          no_elements, nopt, iflag, extension, lscreen_master)
         if (iflag .ne. 0) return  
 
!------------------------------------------------------------------------------
!
!    JConf    M.Zarnstorff  Aug. 2002
!
!------------------------------------------------------------------------------

         call chisq_jconf(sigma_jconf, nrad, no_elements, nopt, 
     1      nu2_b/mskip, iflag, extension, comm1, lscreen_master)
         if (iflag .ne. 0) return  
!-----------------------------------------------------------------------------------------------
!
!     Minimize peak curvature of outer flux surfaces
!     at four different toroidal planes
!     (phi = 0, 90, 180, and 270 degrees)
!     added by d.a. spong
      if( sigma_curv < bigno) then
         call chisq_curvature(sigma_curv, no_elements, nopt)
      endif
!-----------------------------------------------------------------------------------------------
!
!     Maintain boundary between Target_rmax, Target_rmin, and below Target_zmax
!
      call chisq_rzbdy(Target_Rmax, rmax_opt, sigma_rmax, 
     1         ivar_center, no_elements, nopt, .true.)
      call chisq_rzbdy(Target_Rmin, rmin_opt, sigma_rmin, 
     1         ivar_center, no_elements, nopt, .false.)
      call chisq_rzbdy(Target_Zmax, zmax_opt, sigma_zmax, 
     1         ivar_zmax, no_elements, nopt, .true.)

     
! -----------------------------------------------------------------------
!
!     3D boundary limiting                      R.Hatcher & M.Zarnstorff
!
     
      if (sigma_vv < bigno .or. sigma_vv_rms < bigno .or.
     1    sigma_vv_max < bigno ) then
         call chisq_3dbd(sigma_vv, ivar_vv, sigma_vv_rms,  
     1        sigma_vv_max, target_vv, target_vv_rms, lvv_tgt_min,
     2        ivar_vv_rms, ivar_vv_max, no_elements, nopt, iflag, 
     3        lscreen_master)
      end if   

      if (sigma_bd < bigno .or. sigma_bd_rms < bigno .or.
     1    sigma_bd_max < bigno ) then
         call chisq_3dbd_bd(sigma_bd, ivar_vv, sigma_bd_rms,  
     1        sigma_bd_max, target_bd, target_bd_rms,
     2        ivar_vv_rms, ivar_vv_max, no_elements, nopt, iflag, 
     3        lscreen_master)
      end if   
! -----------------------------------------------------------------------
!     Regularize Free-Bdry Coil Currents            M.Zarnstorff
!
      call chisq_extcur(sigma_extcur, target_extcur, 
     1                  ivar_extcur, no_elements, nopt)
     
! -----------------------------------------------------------------------
!     Linear constrain on Coil Currents             M.Zarnstorff
!
      call chisq_OH(sigma_oh, ivar_oh, no_elements, nopt)
     
! -----------------------------------------------------------------------
!     <kappa> : surface-averaged elongation         M.Zarnstorff
!
      call chisq_kappa(target_kappa, sigma_kappa, ivar_kappa, 
     1                 no_elements, nopt)
     
!------------------------------------------------------------------------------
!
!     Suppress targetted resonant Boozer Jacobian components on their
!     respective resonant surfaces.
!     The resonances are indicated by each pair of entries (n_jac, m_jac)
!     with the sigmas being the respective sigma_jac entry
!
!     MCZ  Jan. 01

      call chisq_jac (iota_opt, hs, ivar_jac, 
     1   no_elements, ns_booz, ns_booz_max, mnboz, nrad, nopt,
     2   lscreen)

! -----------------------------------------------------------------------
!
!     Limit ellipticity of the phi=0 cross section     A. Ware, D. Spong, S. Hirshman 
!
      call chisq_ellipticity (Target_ellipticity, sigma_ellipticity,
     1        no_elements, nopt)

!-----------------------------------------------------------------------------------------------
!
!     Minimize coil-complexity measure (>1, the minimum value)
!     and maximum current density
!
      if (lnescoil_opt) call chisq_nescoil 
     1    (no_elements, nopt, iflag, extension, lscreen_master)

!-----------------------------------------------------------------------------------------------
!
!     Island due to 3D Shaping                       L-P. Ku
!
      call chisq_vac_island(ivar_vac_island, no_elements, nopt, iflag,
     1     extension, lscreen_master)
!-----------------------------------------------------------------------------------------------
!
!     Free memory associated with boozer transform arrays
!
      call read_boozer_deallocate
      
!-----------------------------------------------------------------------------------------------
      if (nopt .le. 0) then
         nopt = no_elements
         deallocate (bmin_b, bmax_b, bmod_b)
         return 
      else if (nopt .ne. no_elements) then
         if (myid .eq. master)                                           !MPI      
     1       print *,' NOPT = ',nopt,'!= NO_ELEMENTS = ',no_elements
         iflag = -11
         return
      else if (.not.lfreeb) then
!        Check that boundary data are correctly read in
         call chk_rzmnb(xc_opt, rmnc_bdy, zmns_bdy, xm_bdy, xn_bdy, 
     1        iflag)
         if (iflag .ne. 0) return
      endif

!
!
!     COMPUTE TOTAL CHISQ
!
      epsm = epsilon(wegt(1))**2
      wegt(:nopt) = max (abs(wegt(:nopt)), epsm)
      where (abs(wegt(:nopt)) .gt. epsm)
        fvec(:nopt) =
     1      (chisq_target(:nopt) - chisq_match(:nopt))/wegt(:nopt)
      elsewhere
        fvec(:nopt) = 0
      end where
        
      
      chisq_tot = sum(fvec(:nopt)**2)
 
!-----------------------------------------------------------------------------------------------
!     COMPUTE INDIVIDUAL CHISQs
!-----------------------------------------------------------------------------------------------      
   
      do i = 1, ntargets
         nweight(i) = count(index_array.eq.i .and. abs(wegt).gt.epsm)
         if (nweight(i) .eq. 0) then
            chisq(i) = 0
            mean_weight(i) = 0
         else
            chisq(i) = sum (fvec(:nopt)**2, 
     1           mask = (index_array.eq.i .and. abs(wegt).gt.epsm))
            mean_weight(i) = sum (abs(one/wegt(:nopt)**2),  
     1           mask = (index_array.eq.i .and. abs(wegt).gt.epsm))
         end if
      end do
 
 
      bmin_0 = 0
      bmin_pi2 = 0
      bmin_pi = 0
      bmax_0 = 0
      bmax_pi2 = 0
      bmax_pi = 0

      if( ns_booz_max > 0) then
         nsval = ns_booz(ns_booz_max)

         if( any(lneed_modB(1:ns_booz_max))) then
            bmin_0 = bmin_b(1)
            bmin_pi2 = bmin_b(nu2_b/2)
            bmin_pi = bmin_b(nu2_b)
            bmax_0 = bmax_b(1)
            bmax_pi2 = bmax_b(nu2_b/2)
            bmax_pi = bmax_b(nu2_b)
         endif
      else
         nsval = 0
      endif

      dmax_j = twopi * maxval(abs(buco_opt(2:nrad))) / dmu0
 
!
!     Begin writing to output record file
!
      if (lscreen_master) then 
         write (*, 544) aspect_opt, chisq_tot, ncnt, nsval, 
     1                  nrad, bmin_0, bmin_pi2, bmin_pi, bmax_0, 
     2                  bmax_pi2, bmax_pi,  dmax_j, wp_opt/wb_opt               
         if (sigma_vv < bigno .or. sigma_vv_rms < bigno .or.
     1      sigma_vv_max < bigno ) 
     1      write (*, 545) vv_dist, vv_dist_rms, vv_dist_max
         if (sigma_bd < bigno .or. sigma_bd_rms < bigno .or.
     1      sigma_bd_max < bigno ) 
     1      write (*, 547) bd_dist, bd_dist_rms, bd_dist_max
      endif                                                                 
      write (iunit_opt_local, 544) aspect_opt, chisq_tot, ncnt, nsval, 
     1    nrad, bmin_0, bmin_pi2, bmin_pi, bmax_0, bmax_pi2, bmax_pi, 
     2    dmax_j, wp_opt/wb_opt                                    
      if (sigma_vv < bigno .or. sigma_vv_rms < bigno)
     1   write (iunit_opt_local, 545) vv_dist,vv_dist_rms,vv_dist_max
      if (sigma_bd < bigno .or. sigma_bd_rms < bigno)
     1   write (iunit_opt_local, 547) bd_dist,bd_dist_rms,bd_dist_max
 
!
!     output either coil current or boundary coefficients for intermediate restart
!
      if (LDIAG_OPT) then
        write(iunit_opt_local,*)
        if( lfreeb ) then
            write(iunit_opt_local,'(a)') " EXTCUR ="
            write(iunit_opt_local,*) extcur
            write(iunit_opt_local,'(a,es14.7)') " PHIEDGE = ",phiedge
        else
          call write_rbzb(iunit_opt_local, i)
        endif
      end if

      do i = 1, ntargets
         if (nweight(i) .gt. 0)
     1      write (iunit_opt_local, 646) 
     2      descript(i), chisq(i), sqrt(chisq(i)/mean_weight(i)),
     3      sqrt(mean_weight(i)/nweight(i)), nweight(i)
      end do
 646  format(' Chi-Sq(',a,')=',1pe10.3, ' RMS Err=', 1pe10.3, 
     1     ' <Wgt>=', 1pe10.3, ' nw=',i4)
 
!     do i = 1, nopt
!        write (iunit_opt_local, 546) chisq_descript(i), i, fvec(i),
!    1      chisq_target(i), chisq_match(i), wegt(i)
!     end do
 
 544  format(/' Aspect Ratio = ',f10.3,' Chi-Sq = ',1pe10.3,
     1   ' Function Evaluations = ',i5,/,' At js = ',i3,
     2   ' out of ',i3,' radial points:',/,' Bmin_0 = ',1pe10.3,
     3   ' Bmin_pi/2 = ',1pe10.3,' Bmin_pi = ',1pe10.3,/,' Bmax_0 = ',
     4   1pe10.3,' Bmax_pi/2 = ',1pe10.3,' Bmax_pi = ',1pe10.3,/,
     5   ' Maximum Toroidal Current = ',1pe10.3,' Average Beta = ',
     6   1pe10.3)

 545  format(' Min Plasma-VV Sep =', 1pe11.3, 
     1          '  RMS dist =', 1pe11.3,
     2          '  Max dist =', 1pe11.3)
 
     
 546  format(a,3x,' fvec(',i3,') = ',d10.2,' Target = ',d10.2,
     1   ' Match = ',d10.2,' Sigma(with units) = ',d10.2)
 
 547  format(' Min Plasma-BD Sep =', 1pe11.3, 
     1          '  RMS dist =', 1pe11.3,
     2          '  Max dist =', 1pe11.3)

      deallocate (bmin_b, bmax_b, index_array, wegt,      
     1            chisq_match, chisq_target, chisq_descript, 
     2            bmod_b, stat=n)


      end subroutine load_target

      
!==========================================================
!     HELPER AND PHYSICS ROUTINES
!==========================================================
      subroutine load_physics_codes(executable, input, command_line,
     1    output, extension, iunit, iflag)
      use kind_spec
      use read_boozer_mod
      use safe_open_mod
      use system_mod
      use mpi_params                                                     !MPI
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: iunit, iflag
      character*(*), intent(in) :: executable, input, output, extension,
     1    command_line      
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: num_exec = 10
      character*12, dimension(num_exec), parameter :: exec_string = 
     1 ( / 'xvmec2000   ', 'xbootsj     ', 'xbooz_xform ', 
     2     'xj_invariant', 'xcobravmec  ', 'xbnorm      ',
     3     'xnescoil    ', 'xdkes       ', 'xtprp       ',
     4     'xneo        ' / )
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat, index_exec, iexec, i, ierror
      integer :: iprint
      character(LEN=LEN_TRIM(INPUT) + LEN_TRIM(EXTENSION) 
     1   + LEN_TRIM(COMMAND_LINE) + 1) :: input_ext 
      character(LEN=LEN_TRIM(OUTPUT) + LEN_TRIM(EXTENSION) + 1) :: 
     1    output_ext
C-----------------------------------------------
!     THIS ROUTINE EXECUTES A PHYSICS CODE (CALLED executable) USING A
!     SYSTEM CALL. IT HANDLES POSSIBLE ERRORS DUE TO FILE NON-EXISTENCE, 
!     OUT OF PROCESSORS, ETC. IT READS/OPENS (FOR IUNIT) APPROPRIATE FILES
!     PRODUCED BY THE SUCESSFUL EXECUTION OF THE PHYSICS CODE.
!     
!     INPUT VARIABLES
!
!     EXECUTABLE:   physics code to execute using system call
!     INPUT:        input file name, without extension
!     COMMAND_LINE: screen, noscreen, lreset
!     OUTPUT:       data output file (without extension) resulting
!                   from execution of the executable
      
      input_ext = trim(input)
      if (input(1:1) .ne. ' ')
     1    input_ext = trim(input) // '.' // trim(extension) 
      output_ext = trim(output) // '.' // trim(extension)
      iflag = 0
!
!     SPECIAL CASES
!
      iexec = -1
      do i = 1, num_exec
         index_exec = index(executable, trim(exec_string(i)))         
         if (index_exec .gt. 0) then
            iexec = i
            exit
         end if   
      end do
      
!!$      if (iexec .eq. -1) then
!!$         print *,'Executable not found: ', trim(executable)
!!$         iflag = -17
!!$         return
!!$      end if
      
      input_ext = trim(input_ext) // trim(command_line)

      iprint = -1
 10   call system(trim(executable) // ' ' // trim(input_ext), ierror)

!
!     Parallel processing control - wait for process to execute in system call
!     
      if (ierror.lt.127 .and. ierror.ne.0) then
         iprint = iprint+1
         if (iprint .eq. 0) then         
            print *,' Processor ', myid,': System error = ', ierror, 
     1        ' EXECUTING  ',
     2        trim(executable),' ', trim(input_ext)
         end if   
         if (iprint .le. 10) then
            go to 10
         else
            iflag = -2
            return
         end if   
      end if   

!
!     FOR XVMEC2000, MUST READ IN DATA FROM WOUT FILE GENERATED BY SUCCESSFUL RUN
!     FOR XBOOZ_XFORM, READ IN DATA FROM BOOZMN FILE
!     FOR ALL OTHER EXECUTABLES, OPEN OUTPUT FILE PRODUCED BY EXECUTABLE
!
      if (iexec == 1) then
         call read_wout_opt(.true., extension, ierror, istat)
         if (ierror .ne. 0) then
            iflag = -8
         else if (istat .ne. 0) then
            iflag = -1
         end if      
         if (iflag .ne. 0) return         
      else if (iexec == 3) then
         call read_boozer_file (extension, ierror, istat)
         if (ierror .ne. 0) then
            iflag = -26
            return
         end if
      else
         call safe_open(iunit, istat, output_ext, 'old', 
     1      'formatted')
      end if   
      
!
!     ADDITIONAL ERROR CONTROL (necessary for parallel processing control)
!     
      if (istat .ne. 0) then
         print * ,'Processor ',myid,': Error opening ',trim(output_ext), 
     1         ' in STELLOPT load_physics_codes: istat = ', istat
         iflag = -2
      end if
      
      end subroutine load_physics_codes
      
      subroutine init_ajax 
c-----------------------------------------------
c   M o d u l e s
c-----------------------------------------------
      use kind_spec
      use AJAX_MOD
      use vmec_input, only: phiedge
      use optim_params, only: rgrid_min, rgrid_max, zgrid_min, zgrid_max
      use optim
      use boozer_params
      use safe_open_mod
      use system_mod
      implicit none      
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
 
c-----------------------------------------------
c   L o c a l   V a r i a b l e s
c-----------------------------------------------
      integer :: i, iflag
      integer, dimension(mnmax_opt) :: m, n
      real(rprec), dimension(nrad, mnmax_opt) :: r_aj, z_aj
      real(rprec), dimension(nrad-1, mnmax_opt) :: lam_aj
      real(rprec) :: hs
      real(rprec), dimension(nrad) :: s_aj, sl_aj
      character*120 :: message

c-----------------------------------------------
      r_aj = transpose(rmnc_opt)
      z_aj = transpose(zmns_opt)
      lam_aj = transpose(lmns_opt(:,2:nrad))

      m = nint(xm_bdy)
      n = nint(xn_bdy) * nfp_opt

      hs = 1/real((nrad - 1),rprec)
      s_aj(1) = 0
      do i=2, nrad
         s_aj(i) = (i-1)*hs
         sl_aj(i-1) = (i-1.5_rprec)*hs
      enddo

      call ajax_load_rzlam(nrad, mnmax_opt, s_aj, m, n, r_aj, z_aj,  
     1         iflag, message, 
     1         K_GRID=1, L_MFILTER_AJAX=.false., NRHO_AJAX=nrad, 
     2         NTHETA_AJAX=4*mpol_opt+1,
     2         NZETA_AJAX=4*ntor_opt+1, RHOMAX_AJAX=1._rprec, 
     2         NR_LAM=nrad-1, NK_LAM=mnmax_opt, 
     3         RHO_LAM=sl_aj, LAM=lam_aj)

      if(iflag < 0) then
         print *,' Init_Ajax Load_RZLam warning:'
         print *,message

      else if(iflag > 0) then
         print *,' Init_Ajax Load_RZLam Error:'
         print *,message
         stop
      endif

      call ajax_load_magflux(phiedge, 1, nrad, s_aj, iota_opt, iflag, 
     1                       message)

      if(iflag < 0) then
         print *,' Init_Ajax Load_MagFlux warning:'
         print *,message

      else if(iflag > 0) then
         print *,' Init_Ajax Load_MagFlux Error:'
         print *,message
         stop
      endif
      end subroutine init_ajax


      subroutine generate_mgrid (extension, lscreen, iflag)
c-----------------------------------------------
c   M o d u l e s
c-----------------------------------------------
      use kind_spec
      use vmec_input, only: nzeta, extcur
      use optim_params, only: rgrid_min, rgrid_max, zgrid_min, zgrid_max
      use optim, only: rmax_opt, rmin_opt, zmax_opt, aminor_opt, 
     1    home_dir
      use safe_open_mod
      use system_mod
      use coilsnamin, only: lmodular, lsaddle, lbcoil, lvf, lsurfv,
     1    lsadsfv
      implicit none      
c-----------------------------------------------
c   D u m m y   A r g u m e n t s
c-----------------------------------------------
      character*(*), intent(in) :: extension
      integer, intent(out) :: iflag
      logical :: lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0
c-----------------------------------------------
c   L o c a l   V a r i a b l e s
c-----------------------------------------------
      integer :: ierr, iunit=10, i1, j1
      real(rprec) :: rmin, rmax, zmin, zmax, delta
      character*200 :: parse_string, output_file, line
      logical :: lexist
c-----------------------------------------------
!
!     Generates an mgrid file (for free-boundary calculation) corresponding to present state of coils
!     Parse the extcur file to read in extcur array
!
      
!     Step 1: generate coils.extension file by calling xcoilgeom (the reduced coilopt executable)
!     It reads the COILSIN namelist in the INPUT.EXTENSION file. Redirects screen output from
!     xcoilgeom to a file, COIL_TARGETS.EXTENSION, for future parsing to get engineering target 
!     chi-sq values
!
      output_file = 'coil_targets.' // trim(extension)
      parse_string = trim(home_dir) // '/xcoilgeom ' // trim(extension) 
     1       // ' > ' // trim(output_file)
      if (lscreen) print *, 'Creating MGRID file'
      call system (parse_string, ierr)
      if (ierr.lt.127 .and. ierr.ne.0) then
         if (lscreen) print *, 'XCOILGEOM failed in generate_mgrid: ',
     1      'ierr = ',ierr
         iflag = -22
         return
      end if
      
!     b. Make sure coils.extension file was successfully created

      inquire (file='coils.' // trim(extension), exist=lexist)
      if (.not.lexist) then
         iflag = -23
         return
      end if
      
!     Step 2: generate actual mgrid.extension file from coils.extension
!     a. Estimate grid box dimensions from wout file dimension (previous run of vmec)
!        or rgrid_min,max, zgrid_min,max in namelist
!        use plasma max, min +- minor radius, and adjust to make delr ~ delz
      if (rmax_opt .gt. rmin_opt) then
         rmin = max(rmin_opt - aminor_opt, epsilon(aminor_opt))
         zmin =-(zmax_opt + aminor_opt)
         rmax = rmax_opt + aminor_opt
         zmax =-zmin
         delta = ((rmax - rmin) - (zmax - zmin))/2
         if (delta .gt. zero) then
            zmax = zmax + delta
            zmin = zmin - delta
         else if (delta .lt. zero) then
            rmax = rmax - delta
            rmin = rmin + delta
            if (rmin .lt. epsilon(aminor_opt)) then
               rmax = rmax - rmin
               rmin = epsilon(aminor_opt)
            end if   
         end if
         rgrid_min = rmin;  rgrid_max = rmax
         zgrid_min = zmin;  zgrid_max = zmax
      end if
      if (rgrid_min.ge.rgrid_max .or. zgrid_min.ge.zgrid_max) then
            print *, 'rgrid_min MUST BE < rgrid_max'
            print *, 'zgrid_min MUST BE < zgrid_max'
            stop            
      end if

!     b. Parse command line for xgrid file

      write (parse_string, '(3a,1x,a1,1x,4(1pe20.12,1x),i4)') 
     1   trim(home_dir),'/xgrid ', trim(extension), 'T', 
     2   rgrid_min, rgrid_max, zgrid_min, zgrid_max, nzeta
     
      if (.not.lscreen) parse_string = trim(parse_string) // 
     1   ' > /dev/null'
      

!     c. Make system call to run xgrid, which generates mgrid file
      call system (parse_string, ierr)

!     d. Clean up by removing coils file
      call system ('rm -f coils.' // trim(extension))

!     e. Parse extcur file
      parse_string = "extcur." // trim(extension)
      call safe_open(iunit, ierr, parse_string, 'old', 'formatted')
      if (ierr .ne. 0) then
         print *,
     1   ' Error opening extcur file in STELLOPT generate_mgrid: ',
     2   ' ierr = ', ierr         
         iflag = -24
         return
      end if
      
      extcur = 0
      
      do while (ierr .eq. 0)
         read (iunit,'(a)', iostat=ierr) line
         if (ierr .ne. 0) exit
         i1 = index(line, '(')
         if (i1 .eq. 0) exit
         j1 = index(line,')') - 1
         read (line(i1+1:j1), *) j1
         if (j1.le.0 .or. j1.gt.size(extcur)) cycle
         i1 = index(line,'=')
         if (i1 .eq. 0) exit
         read (line(i1+1:), *) extcur(j1)
      end do
      
      close (iunit)
      
      iflag = 0

      end subroutine generate_mgrid


      subroutine chisq_aspect(match, target, sigma, num, nopt)
      use kind_spec
      use chisq_mod
      use optim, only: bigno, laspect_max
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: nopt
      integer :: num
      real(rprec), intent(in) :: match, target, sigma
C-----------------------------------------------
      if (abs(sigma) .ge. bigno) return 

      num = num + 1
      if (nopt .gt. 0) then
         index_array(num) = ivar_aspect
         wegt(num) = sigma
         chisq_target(num) = target
!
!        if (laspect_max), only contribute to chisq if match > target 
!       (match < target is permitted without penalty)
!
         if (laspect_max) then
            chisq_match(num) = max(match, target)
         else
            chisq_match(num) = match
         end if
          
         chisq_descript(num) = descript(ivar_aspect)
      end if

      end subroutine chisq_aspect
      
      
      subroutine chisq_rbtor(match, target, sigma, num, nopt)
      use kind_spec
      use chisq_mod
      use optim, ONLY: bigno
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: nopt
      integer :: num
      real(rprec) :: match, target, sigma
C-----------------------------------------------
      if (abs(sigma) .ge. bigno) return 

      num = num + 1
      if (nopt .gt. 0) then
         index_array(num) = ivar_rbtor
         wegt(num) = sigma
         chisq_target(num) = target
         chisq_match(num) = match
         chisq_descript(num) = descript(ivar_rbtor)
      end if

      end subroutine chisq_rbtor
      
      
      subroutine chisq_ballooning (sigma, target, lflip, num, ns_surf,  
     1     nrad, nopt, iflag, extension, command)
      use kind_spec
      use chisq_mod
      use optim_params, ONLY: bal_theta0, bal_zeta0,
     1    lpres_prof_opt, lballoon_opt, nini_theta, nini_zeta, nini_tot
      use optim, ONLY: home_dir, version_opt
      use safe_open_mod
      use mpi_params                                                     !MPI
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: nrad, nopt, ns_surf(*)      
      integer :: num, iflag
      real(rprec) :: sigma(*), target(*)
      character*(*) :: extension, command
      logical :: lflip
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ij, jk, jrad, j, k, iunit, jcount, ns_ball, num1
      real(rprec) :: grate_ball(nrad), init_zeta, init_theta
      character(len=len_trim(home_dir)+20) :: version
      logical :: ex
C-----------------------------------------------
!
!     BALLOONING STABILITY CRITERION (GRATE>0, unstable)
!
      version = trim(home_dir) // '/xcobravmec'
      if (nopt. gt. 0) then
         iunit = unit_cobra
         call safe_open(iunit, k, 'in_cobra.'//trim(extension),
     1          'replace', 'formatted')
         if (k .ne. 0) then
            iflag = -12
            return
         endif
         write(iunit, *,iostat=k) nini_zeta, nini_theta, trim(extension)
         write(iunit, *,iostat=k) (bal_zeta0(jk), jk = 1, nini_zeta)
         write(iunit, *,iostat=k) (bal_theta0(ij), ij = 1, nini_theta)
         write(iunit, *,iostat=k) nrad
         write(iunit, *,iostat=k) ns_surf(1:nrad)
         close (iunit)
         if (k .ne. 0) then
            print *,' chisq_balloon error for file: ', trim(extension)
            iflag = -12
            return
         end if   

!        Read in COBRA results: 
!        ns_ball    =  index of VMEC radial surface at which growth rate is computed
!                      (should be same as ns_balloon_max=nrad, so no need to save it ...)
!        grate_ball(j) = growth rate at initial position in ij, lk loop at surface "j"
!                      (note that index-j in list <=> ns_ball in VMEC indexing)

         call load_physics_codes (version, 'in_cobra', command, 
     1             'cobra_grate', extension, iunit, iflag)
         if (iflag .ne. 0) return

         do ij = 1, nini_theta
            do jk = 1, nini_zeta
               read (iunit, *, iostat=k) init_zeta, init_theta, jcount
               if (jcount.ne.nrad .and. myid.eq.master)                 !MPI 
     1            print *,' JCOUNT = ', jcount,' != NS_SURF_MAX (= ',
     2                    nrad,' READING COBRA_GRATE FILE'         
               if (k .eq. 0) read (iunit, *, iostat=k) (ns_ball, 
     1              grate_ball(jrad), jrad = 1, jcount)
               if (k .ne. 0) then
                  if (myid .eq. master)                                   !MPI      
     1               write (6, *) 'Error reading cobra_grate in ',
     2               ' load_targe for file: ', trim(extension)
                  iflag = -12
                  return
               end if     
 
               num1 = num+1
               num = num + nrad

               index_array(num1:num) = ivar_balloon
               wegt(num1:num) = 0.1_dp * sigma(:nrad)
               chisq_target(num1:num) = target(:nrad)
               if( .not. lflip) then
                  chisq_match(num1:num) = 
     1                           max(grate_ball(:nrad), target(:nrad))
               else
                  chisq_match(num1:num) = 
     1                           min(grate_ball(:nrad), target(:nrad))
               endif
               chisq_descript(num1:num) = descript(ivar_balloon)

            end do
         end do
         
         close(unit=iunit)                                             !Keep grate file 
 
      else
         inquire(file=trim(version), exist=ex)
         if (.not.ex) then
            if (myid .eq. master)                                       !MPI
     1          print *,'xcobra file not found in ' // trim(home_dir)
            lpres_prof_opt = .false.
            lballoon_opt = .false.
         else if (nrad .ge. 1) then
            num = num + nrad*nini_tot
         endif
         iflag = 0
      endif

      end subroutine chisq_ballooning


      subroutine chisq_beta (Target, sigma, targetW, sigmaW, 
     1                       num, nopt)  
      use kind_spec
      use chisq_mod
      use optim, ONLY: wp_opt, wb_opt, bigno, lbeta_min, 
     1                 vp_opt, pres_opt, nrad
      use vparams, ONLY: twopi
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: nopt
      integer :: num
      real(rprec) :: target, sigma, targetW, sigmaW
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0
      real(rprec), parameter :: min_beta = 1.0e-4_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec) :: match, Ekin
C-----------------------------------------------
      if (abs(sigma) < bigno) then
         num = num + 1
         if (nopt .gt. 0) then
            index_array(num) = ivar_beta
            match = zero
            if (wb_opt .gt. zero) match = wp_opt/wb_opt
            chisq_target(num) = Target
!
!        matched value of beta can exceed Target value with no penalty
!         
            if (lbeta_min) then
               chisq_match(num) = min(match, Target)
            else
               chisq_match(num) = match
            end if
            wegt(num) = abs(sigma*max(Target, min_beta))
            chisq_descript(num) = descript(ivar_beta)
         end if
      end if

      if (abs(sigmaw) < bigno) then
         num = num + 1
         if (nopt .gt. 0) then
            index_array(num) = ivar_beta
            ekin = 1.5_rprec * twopi**2 *
     1           sum(vp_opt(2:nrad)*pres_opt(2:nrad))/(nrad-1)

            chisq_match(num) = ekin
!            print *, 'Chisq_beta: ',ekin, targetw
            chisq_target(num) = Targetw
            wegt(num) = abs(sigmaw)
            chisq_descript(num) = descript(ivar_beta)
         end if
      end if

      end subroutine chisq_beta
      
      
      subroutine chisq_curtor (Target, sigma, num, nopt)  
      use kind_spec
      use chisq_mod
      use optim, ONLY: bigno
      use vmec_input, ONLY: curtor
      use vparams, ONLY: twopi
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: nopt
      integer :: num
      real(rprec) :: target, sigma
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------

C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------

C-----------------------------------------------
      if (abs(sigma) < bigno) then
         num = num + 1
         if (nopt .gt. 0) then
            index_array(num) = ivar_curtor
            chisq_target(num) = Target
            chisq_match(num) = curtor
            wegt(num) = abs(sigma)
            chisq_descript(num) = descript(ivar_curtor)
         end if
      end if

      end subroutine chisq_curtor
      
      
      subroutine chisq_jedge (Target, sigma, num, nopt)  
      use kind_spec
      use chisq_mod
      use optim, ONLY: bigno
      use vmec_input, ONLY: ac
      use vparams, ONLY: twopi
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: nopt
      integer :: num
      real(rprec) :: target, sigma
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------

C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec) :: sum1
C-----------------------------------------------
      if (abs(sigma) < bigno) then
         num = num + 1
         if (nopt .gt. 0) then
            sum1 = sum(ac(0:10))

            index_array(num) = ivar_jedge
            chisq_target(num) = Target
            chisq_match(num) = sum1
            wegt(num) = abs(sigma)
            chisq_descript(num) = descript(ivar_jedge)
         end if
      end if

      end subroutine chisq_jedge
      
      
      subroutine chisq_dkes(sigma, num, nopt, iflag, nsval, 
     1   extension, command)
      use kind_spec
      use chisq_mod
      use optim_params, ONLY: ldkes_opt, ndkes_mask, dkes_nu, 
     1                        dkes_efield, nsurf_mask
      use optim, ONLY: home_dir, iunit_dkes
      use mpi_params                                                     !MPI
      implicit none
c---------------------------------------------------------
c    D u m m y  A r g u m e n t s
c---------------------------------------------------------
      integer, intent(in) :: nopt, nsval
      integer, intent(inout) :: num, iflag
      real(rprec), intent(in) :: sigma
      character*(*), intent(in) :: extension, command
c
c---------------------------------------------------------
c    L o c a l s  P a r a m e t e r s
c---------------------------------------------------------
      integer :: itest
      real (rprec) :: zero = 0
      character*(*), parameter :: input_file = 'boozer', 
     1   output_file = 'opt_dkes'
c---------------------------------------------------------
c    L o c a l  V a r i a b l e s
c---------------------------------------------------------
      character*30 :: dkes_args
      character(len=len_trim(home_dir)+20) :: version
      character(len=len_trim(command)+ 40) :: comdline
      integer :: isurf_dkes, iunit, ierror, istat
      real (rprec) ::  L11p, L33p, L31p
      real (rprec) ::  L11m, L33m, L31m
      real (rprec) ::  scal11, scal33, scal13
      real (rprec) ::  small
      logical, parameter :: ldebug = .false.
c---------------------------------------------------------
 
       if (ldebug .and. myid.eq.master) then

        print *, "********** Beginning chisq_dkes **************"
        print *, "nopt = ", nopt
        print *, "num = ", num
        print *, "ldkes_opt = ", ldkes_opt
        print *, "surf_dkes = ", ndkes_mask(nsval)
        print *, "dkes_nu = ", dkes_nu(nsval)
        print *, "dkes_efield = ", dkes_efield(nsval)
        print *, "sigma_dkes = ", sigma
        print *, "extension = ", trim(extension)
        print *, "command = ", trim(command)
        print *, "unit_dkes = ", unit_dkes

      end if

      iflag = 0

      if (nopt .gt. 0) then
 
         iunit = unit_dkes
 
!        VERSION: EXECUTABLE FOR XDKES
         version = trim(home_dir) // '/xdkes'
 
!        USE INTERNAL FILE TO WRITE DKES NAMELIST INPUT INTO A CHARACTER VARIABLE
!        NAMELIST INPUT = SURFACE NUMBER, COLLISIONALITY, EFIELD
         write(dkes_args,'(1x,i3,1x,f10.5,1x,f10.5)') 
     1        nsval, dkes_nu(nsval), dkes_efield(nsval)
 
!        ADD DKES NAMELIST INPUT TO COMMAND LINE
!        ADD LOGICAL SCREEN OUTPUT VARIABLE TO COMMAND LINE
!        CALL TO DKES WILL AUTOMATICALLY ACTIVATE THE INTERNAL DKES_INPUT_PREPARE CALL
         comdline = trim(dkes_args) // ' ' // trim(command)

         call load_physics_codes (version, input_file, comdline,
     1         output_file, extension, iunit, iflag)
         if (iflag .ne. 0) return

!        READ OPT_DKES FILE AND STORE DATA IN DKES_OPT FILE FOR ALL SURFACES
         read(iunit,'(3(2x,e24.13))', iostat=ierror) L11p, L33p, L31p
         read(iunit,'(3(2x,e24.13))', iostat=ierror) L11m, L33m, L31m
         read(iunit,'(3(2x,e24.13))', iostat=ierror) 
     1        scal11, scal33, scal13
     
         write(iunit_dkes,'(/,a,i5)') 'RADIAL SURFACE JS =', nsval
         write(iunit_dkes,'(3(2x,e24.13))') L11p, L33p, L31p
         write(iunit_dkes,'(3(2x,e24.13))') L11m, L33m, L31m
         write(iunit_dkes,'(3(2x,e24.13))') scal11, scal33, scal13
     
         close (iunit, status='delete')
         call system("rm -f dkesout." // trim(extension))

         if (ldebug .and. myid.eq.master) then
            print *, "OPT_DKES FILE"
            print *,  L11p, L33p, L31p
            print *,  L11m, L33m, L31m
            print *,  scal11, scal33, scal13
         endif 
         if (ierror .ne. 0) then
            iflag = -20
            return
         end if
 
!        LOAD CHISQ ARRAYS: NOTE THERE IS NO NORMALIZATION FOR L11
!
         num = num + 1
         index_array(num) = ivar_dkes
         wegt(num) = sigma
         chisq_target(num) = zero
         chisq_match(num) = 0.5_dp*(L11p + L11m)
         chisq_descript(num) = descript(ivar_dkes)
 
      else
         num = num + 1
      end if
 
      if (ldebug .and. myid.eq.master) then
         print *, "at end of chisq_dkes"
         print *, "nopt = ", nopt
         print *, "num = ", num
         print *, "********** Ending chisq_dkes **************"
      end if
 
      end subroutine chisq_dkes


      subroutine chisq_pseudo (sigma, sigma2, hs, bmn_b, ixn_b, ixm_b, 
     1           num, nopt, nsval, lwgt, lscreen)
      use kind_spec
      use boozer_params, only: mnboz
      use optim, only: nfp_opt, iota_opt, bigno
      use chisq_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: num, nopt, nsval, ixn_b(*), ixm_b(*)
      real(rprec) :: sigma, sigma2, hs, bmn_b(mnboz)
      logical :: lwgt, lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, parameter :: nu = 496, nu1 = nu+1, n0 = nu/2
      integer :: mn, jcount, m, n, i
      real(rprec) :: pi, du, otaper, u, v, b, BBB, wrip, wrip2, rhoj,
     1               wrip_t, wrip2_t, BB_max, wrip_w, wrip2_w, fctr
      real(rprec), dimension(nu1) :: BB
C-----------------------------------------------
      if (nopt .gt. 0)then

         wrip = 0
         wrip2 = 0
         wrip_t = 0
         wrip2_t = 0
         wrip_w = 0
         wrip2_w = 0

         rhoj = hs*(nsval - 1.5_dp)               

!
!  "WATER" - CALCULATION OF THE MEASURE OF RIPPLES. A.SUBBOTIN 11/98
!
         pi = 4*atan(1._dp)
         du = (2*pi)/nu
         otaper = iota_opt(nsval)/nfp_opt

         do i = n0, nu1
            u = -pi + (i-1)*du
            v = u/otaper
            BB(i) = 0

            do mn = 1,mnboz
               n = ixn_b(mn)/nfp_opt
               m = ixm_b(mn)
               b = bmn_b(mn)
               BB(i) = BB(i) + b*cos(m*u-n*v)
            enddo
         enddo

         BBB = -1
         BB_max = maxval(bb(n0:nu1))

         do i = n0, nu1
            u = -pi + (i-1)*du
            fctr = abs(sin(u))

            if(BB(i).gt.BB(i-1) .and. BB(i).ge.BB(i+1) .and. 
     1         BB(i).ge.BBB) BBB = BB(i)

            wrip_t = wrip_t + (BB_max-BB(i))*du
            wrip2_t = wrip2_t + du
   
            if(BB(i).lt.BBB)  then
               wrip=wrip+(BBB-BB(i))*du      ! ripple depth
               wrip2 = wrip2 + du            ! ripple width,  MCZ Feb.99

               wrip_w = wrip+(BBB-BB(i))*du*fctr      ! weigted ripple depth     MCZ Apr 06
               wrip2_w = wrip2 + du*fctr            ! weigted ripple width,
            endif
         enddo
         
         if( lscreen) then
            print *,' '
            print *,' Surf ',nsval, 'Water fraction:', wrip/wrip_t,
     1                            ' width fraction:', wrip2/wrip2_t
            print *,'    weighted Water fraction:', wrip/wrip_t,
     1                            ' width fraction:', wrip2/wrip2_t
         endif

         if (sigma .lt. bigno) then
            num = num + 1
            index_array(num) = ivar_pseudo
            wegt(num) = sigma * wrip_t
            chisq_match(num) = wrip
            if( lwgt ) chisq_match(num) = wrip_w
            chisq_descript(num) = descript(ivar_pseudo)
            chisq_target(num) = 0
         end if

         if (sigma2 .lt. bigno) then
            num = num + 1
            index_array(num) = ivar_pseudo
            wegt(num) = sigma2 * wrip2_t
            chisq_match(num) = wrip2
            if( lwgt ) chisq_match(num) = wrip2_w
            chisq_descript(num) = descript(ivar_pseudo)
            chisq_target(num) = 0
         end if   

      else
         if (sigma .lt. bigno) num = num+1
         if (sigma2.lt. bigno) num = num+1
      endif

      end subroutine chisq_pseudo      
      
      
      subroutine chisq_neo (sigma, ivar, num, ns_neo, ns_max, 
     1                      nopt, iflag, extension, lscreen)
      use kind_spec
      use chisq_mod
      use optim, ONLY: home_dir, nrad
      use optim_params, ONLY: lneo_opt
      use system_mod, ONLY: system
      use safe_open_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: ivar, nopt, iflag, num, ns_max
      real(rprec), dimension(nrad) :: sigma
      integer, dimension(ns_max) :: ns_neo
      character*(*) :: extension
      logical :: lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, parameter :: unit_neo = 27
      integer :: iunit, ifail, istat, jcount, nsval, isurf, write_prog
      real(rprec), parameter :: zero = 0
      real(rprec) :: ripple_eff, dummy_neo
      character*200 :: version
      logical :: ex
C-----------------------------------------------
      version = trim(home_dir) // '/xneo'
      
      if (nopt > 0) then

         iunit = unit_neo
         write_prog = 0
            
         if (lscreen) then
           write_prog = 1
           write(6,*) 
     1           'Running NEO for effective ripple calculations'
         endif

!        write out input control filed neo_param.ext for surface control
         call write_neoparam (ns_max, ns_neo, extension, 
     1        write_prog, istat)
         if (istat .ne. 0) then
            iflag = -27
            return
         endif

         call load_physics_codes(version, " ", trim(extension),
     1         'neo_out', extension, iunit, iflag)
         if (iflag .ne. 0) return

         if (lscreen) write(6,*) '  ns    ripple diffusion (eps_h**1.5)'

         do jcount = 1, ns_max

              nsval = ns_neo(jcount)
              read(unit_neo,*,iostat=istat) isurf, ripple_eff, 
     1            dummy_neo, dummy_neo, dummy_neo, dummy_neo
              if (istat .ne. 0) then
                 iflag = -27
                 return
              end if
              if( lscreen) write (6, 1000) nsval, ripple_eff
              num = num + 1
              index_array(num)  = ivar
              chisq_target(num) = zero
              chisq_match(num)  = ripple_eff
              chisq_descript(num) = descript(ivar)
              wegt(num) = sigma(nsval)

          enddo
          close(unit_neo)

      else
         inquire(file=version, exist=ex)
         if (.not.ex) then
            print *, 'xneo file not found in ' // trim(home_dir)
            lneo_opt = .false.
         else
            num = num + ns_max
         end if 
      endif

 1000 format(2x,i3,3x,1pe12.3)

      end subroutine chisq_neo


      subroutine chisq_orbit (sigma, ivar, num,   
     1                      nopt, iflag, extension, lscreen)
      use kind_spec
      use chisq_mod
      use optim, ONLY: home_dir
      use optim_params, ONLY: lorbit_opt
      use system_mod, ONLY: system
      use safe_open_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: ivar, nopt, iflag, num
      real(rprec)  :: sigma
      character*(*) :: extension
      logical :: lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, parameter :: unit_orbit = 79
      integer :: iunit, ierr, istat, i 
      real(rprec) :: transit_completed, transit_asked,
     1               average_loss_time
      integer :: nzterm, np_start, np_escape
      real(rprec), parameter :: zero = 0
      character*200 :: version
      logical :: ex
C-----------------------------------------------
      
      if (nopt > 0) then

         iunit = unit_orbit
            
         if (lscreen) then
           write(6,*) 
     1           'Running ORBIT for Fast Ion Loss Calculations'
         endif

         ierr = 0
         version = trim(home_dir) // '/xmkjmc'
         version = trim(version) // ' ' // extension
         call system(version,ierr)                              
         if( ierr .lt. 127 .and. ierr .ne. 0 ) then 
            iflag = -28
            return
         endif

         version = trim(home_dir) // '/xeq3d'
         version = trim(version) // ' -b ' // extension
         call system(version,ierr)                              
         if( ierr .lt. 127 .and. ierr .ne. 0 ) then
            iflag = -29
            return
         endif

         version = trim(home_dir) // '/xorbit3d'
         version = trim(version) // ' -b ' // extension
         call system(version,ierr)                              
         if( ierr .lt. 127 .and. ierr .ne. 0 ) then
            iflag = -30
            return
         endif

         call safe_open(iunit, istat,'orbsum.'//trim(extension), 
     1                  'unknown', 'formatted')
         if( istat .ne. 0 ) then
            iflag = -31
            return
         else
            read(unit_orbit,*) nzterm, np_start, np_escape,
     >                   transit_completed, transit_asked,
     >                   average_loss_time

            close(unit_orbit)
         endif
!
             num = num + 1
             index_array(num) = ivar
             chisq_descript(num) = descript(ivar)
             chisq_match(num) =  average_loss_time

c            chisq_match(num) = float(np_escape)
c     >                         /float(np_start)
c     >                         /transit_completed
c     >                         *transit_asked

             chisq_target(num) =  zero 
             wegt(num) = sigma

      else

         version = trim(home_dir) // '/xmkjmc'
         inquire(file=version, exist=ex)
         if (.not.ex) then
            print *, 'xmkjmc not found in ' // trim(home_dir)
            lorbit_opt = .false.
         endif

         ierr=0
         version = "ln -s ../mkjmc_param.in ./mkjmc_param.in"
         call system(trim(version), ierr)
         if (ierr .lt. 127 .and. ierr .ne. 0) then
            print *,' MKJMC control file mkjmc_param.in not linked'
            lorbit_opt = .false.
         endif

         version = trim(home_dir) // '/xeq3d'
         inquire(file=version, exist=ex)
         if (.not.ex) then
            print *, 'xeq3d not found in ' // trim(home_dir)
            lorbit_opt = .false.
         endif

         version = trim(home_dir) // '/xorbit3d'
         inquire(file=version, exist=ex)
         if (.not.ex) then
            print *, 'xorbit3d not found in ' // trim(home_dir)
            lorbit_opt = .false.
         endif

         version = "ln -s ../orbit_param.in ./orbit_param.in"
         call system(trim(version), ierr)
         if (ierr .lt. 127 .and. ierr .ne. 0) then
         print *,' ORBIT control file orbit_param.in not linked'
         lorbit_opt = .false.
         endif

         if( lorbit_opt ) num = num + 1
      endif

      end subroutine chisq_orbit


      subroutine chisq_bmin (sigma, ivar, 
     1          bmin_b, num, nopt, mskip, sqrt_nsurf, dtheta)
      use kind_spec
      use chisq_mod
      use boozer_params, ONLY: nu2_b
      use optim, ONLY: bigno
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: ivar, nopt, mskip
      integer :: num
      real(rprec) :: bmin_b(*), sigma, sqrt_nsurf, dtheta
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: p5 = 0.5_dp
      real(rprec), parameter :: max_variation = 2.e-2_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: n
      real(rprec) :: avg_bmin
C-----------------------------------------------
      if (abs(sigma) .ge. bigno) return

      if (nopt .gt. 0) then

         avg_bmin = sum(bmin_b(:nu2_b))
     1            - p5*(bmin_b(1) + bmin_b(nu2_b))
         avg_bmin = avg_bmin * dtheta

         do n = 1, nu2_b, mskip
           num = num + 1
           index_array(num) = ivar
           wegt(num) = max_variation*sigma*avg_bmin*sqrt_nsurf
           chisq_target(num) = avg_bmin
           chisq_descript(num) = descript(ivar)
           chisq_match(num) = bmin_b(n)
         end do              

      else                   
         do n = 1, nu2_b, mskip
            num = num + 1
         end do   
      end if   

      end subroutine chisq_bmin


      subroutine chisq_jac (iota, hs, 
     1   ivar, num, nsurf, nsurf_max, mnboz, nrad, nopt,
     2   lscreen)
      use kind_spec
      use chisq_mod
      use read_boozer_mod, ONLY: bmn_b, gmn_b, ixm=>ixm_b, ixn=>ixn_b
      use optim, ONLY: nfp_opt, bigno
      use optim_params, ONLY: sigma=>sigma_jac, n_jac, m_jac
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: ivar, nrad, nopt, nsurf_max, nsurf(*),
     1    mnboz
      integer :: num
      logical :: lscreen
      real(rprec) :: iota(*), hs
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, one = 1,
     1   sjmax = one, c1p5 = 1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: jcount, m, n, nsval, mn
      integer :: j, k, nradii, knorm, ns_old

      real(rprec) :: sj, gtotl, gvalue, f, g, ginterp
      real(rprec) :: iota_tgt, iota_old, sj_old, sj_interp 
C-----------------------------------------------

      if (all(abs(sigma(:)) >= bigno)) return

      if( nopt > 0 ) then

        do k=1, mnboz
          if( ixm(k) == 0 .and. ixn(k) == 0) then
            knorm = k
            exit
          endif
        enddo

        if( lscreen) then
           write(6,*) ' Resonant Jacobian Scan'
           write(6,*) ' n  m   s     gmn       gmn_norm'
        endif

        do j=1, size(sigma)
          if( sigma(j) >= bigno .or. n_jac(j) == 0 
     1                          .or. m_jac(j) <= 0 ) cycle
          num = num + 1
          wegt(num) = sigma(j)
          index_array(num) = ivar
          chisq_descript(num) = descript(ivar)
          chisq_match(num) = zero
          chisq_target(num) = zero

c          if( lscreen) print *,'jac', n_jac(j), m_jac(j), sigma(j)


          mn = 0
          do k=1, mnboz
            if( ixm(k) == m_jac(j) .and. 
     1          ixn(k) == n_jac(j)*nfp_opt) then
              mn = k
              exit
            endif
          enddo

          if( mn == 0) then
            wegt(num) = bigno
            cycle
          endif

          iota_tgt = real(n_jac(j)*nfp_opt, rprec)/m_jac(j)
          nradii = 0
          gtotl = zero

          iota_old = iota(nsurf(1))
          ns_old = nsurf(1)
          sj_old = 0

          do jcount = 2, nsurf_max
            nsval = nsurf(jcount)
            sj = hs*(real(nsval,rprec) - c1p5)            !!This is correct (SPH)
            if (sj .gt. sjmax) cycle

c            write(6,'(a,i3,f5.2,4(1pe13.5))') ' ns,s,gmn,bmn=',
c     1          nsval,sj,gmn_b(mn,nsval), gmn_b(knorm,nsval),
c     1          bmn_b(mn,nsval), bmn_b(knorm,nsval)


            if( (iota_tgt > iota_old .and. iota_tgt <= iota(nsval)) .or.
     1        (iota_tgt < iota_old .and. iota_tgt >= iota(nsval))) then
              nradii = nradii + 1
              f = (iota_tgt-iota_old) / (iota(nsval)-iota_old)
              g = one - f

              ginterp = gmn_b(mn,nsval)*f + gmn_b(mn,ns_old)*g

              gvalue = ginterp /
     1                 (gmn_b(knorm,nsval)*f + gmn_b(knorm,ns_old)*g)

              sj_interp = sj*f + sj_old*g
              
              gtotl = gtotl + gvalue**2

              if( lscreen) then
                 write(6,'(2i3, f6.2, (1pe11.2))') 
     1                 n_jac(j), m_jac(j), sj_old, gmn_b(mn,ns_old)
                 write(6,'(2i3, f6.2, 2(1pe11.2))') 
     1                 n_jac(j), m_jac(j), sj_interp, ginterp, gvalue
                 write(6,'(2i3, f6.2, (1pe11.2))') 
     1                 n_jac(j), m_jac(j), sj, gmn_b(mn,nsval)
              endif
            endif

            iota_old = iota(nsval)
            ns_old = nsval
            sj_old = sj
          enddo

          if( gtotl == zero) then
             chisq_match(num) = zero
          else
             chisq_match(num) = sqrt(gtotl)             
          endif
        enddo
      else
        num = num + count( abs(sigma(:)) < bigno  .and. n_jac(:) /= 0
     1                      .and. m_jac(:)> 0)

      endif

      end subroutine chisq_jac


      subroutine chisq_vac_island (ivar, num, nopt, iflag, extension,
     1   lscreen)

      use kind_spec
      use chisq_mod
      use optim, ONLY: home_dir, chisq_min, bigno
      use optim_params, ONLY: sigma=>sigma_vac_island, 
     1                        n_vac_island, m_vac_island
      use system_mod, ONLY: system
      use safe_open_mod
      use vmec_input, ONLY: mgrid_file

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: ivar, nopt  
      integer :: num, iflag
      logical :: lscreen
      character*(*) :: extension
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(kind=rprec), parameter :: zero = 0, one = 1,
     1   sjmax = one, c1p5 = 1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, parameter :: unit_vac = 27
      integer :: jcount, kcount, ierr, istat
      integer :: j, k, iunit, m, n
      real(kind=rprec) :: ds, bnavg, ds_sum
      character(len=len_trim(home_dir)+50) :: version
      character*200 :: tmp
      character*6, parameter :: mgrid_ext='mgrid.'
      logical :: ex
      save jcount
C-----------------------------------------------

      if (.not. any(sigma(:) < bigno)) return

      version=trim(home_dir) // '/xvacisld'
      iunit=unit_vac
      
      if (nopt > 0) then
!
!         if coils.extension does not exist, copy from parent
!         using extension from mgrid

          tmp='coils.'//trim(extension)
          inquire(file=tmp, exist=ex)
          if(.not.ex) then
             j = index(mgrid_file,mgrid_ext)
             if (j .eq. 0 ) then
                iflag = -33
                return
             else
                j = len(mgrid_ext) + 1
                k = len_trim(mgrid_file)
                tmp = "ln -s ../coils." // trim(mgrid_file(j:k)) // 
     1               " ./coils." // trim(extension)
                print *, tmp
             endif
             call system(tmp, ierr)
             if (ierr .lt. 127 .and. ierr .ne. 0) then
                iflag = -33
                return
             endif
          endif

          version=trim(version)// ' ' // extension
          call system(version, ierr)

          if( ierr .lt. 127 .and. ierr .ne. 0 ) then
             print *, "Trouble running XVACISLD, ierr = ", ierr
             iflag = -33     ! temporary flag setting
             return
          endif

          call safe_open(iunit, istat,'visldout.'//trim(extension),
     1                  'unknown', 'formatted')
          if( istat .ne. 0 ) then
            print *, "In Open VISLDOUT istat = ", istat
            iflag = -33
            return
          else
 
            if( lscreen) then
                write(6,*) ' Vac-Island Estimate'
                write(6,*) ' n  m   ds     bnavg'
            endif


            do j=1, jcount

              ds_sum = zero
              read(iunit,*) kcount
              read(iunit,*)
              read(iunit,*)
              do k=1, kcount
                read(iunit,*) m, n, ds, bnavg
                if(m .eq. m_vac_island(j) .and. 
     1             n .eq. n_vac_island(j) ) then

                   ds_sum=ds_sum+ds

                   if( lscreen) write(6,'(2i3, 2(1pe11.2))') 
     1                 n, m, ds, bnavg 

                endif
              enddo
              num = num + 1
              wegt(num) = sigma(j)
              index_array(num) = ivar
              chisq_descript(num) = descript(ivar)
              chisq_match(num) = ds_sum
              chisq_target(num) = zero
              rewind iunit
            enddo

            close(iunit)
          endif

      else

          inquire(file=version, exist=ex)

          if(.not.ex) then
            print *, 'xvacisld not found in ' // trim(home_dir)
            sigma(:) = bigno
          endif
          

          jcount= count( sigma(:) < bigno  .and. 
     1                   n_vac_island(:) /= 0
     1             .and. m_vac_island(:)> 0)
          num = num + jcount
          call safe_open(iunit, istat,'vacisld.param',
     1                  'unknown', 'formatted')
          write(iunit,*) jcount
          do j=1, jcount
            write(iunit,*) m_vac_island(j), n_vac_island(j)
          enddo
          close(iunit)

      end if
 
      end subroutine chisq_vac_island


      subroutine chisq_bmn (bmn_b, sigma, hs, ixm, ixn, num, 
     1       mnboz, nsval, nopt)
      use kind_spec
      use chisq_mod
      use optim, ONLY: nfp_opt, bigno
      use optim_params, ONLY: helicity
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: mnboz, nsval, nopt, 
     1   ixm(mnboz), ixn(mnboz) 
      integer, intent(inout) :: num
      real(rprec), intent(in) :: bmn_b(mnboz), hs, sigma
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, one = 1, c1p5 = 1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: jcount, m, n, mn, num0
      integer :: l_helicity, k_helicity
      real(rprec) :: bnorm, bmn, bmax, rad_wegt, sj
      logical  :: l_symmetry
C-----------------------------------------------
      if (abs(sigma) .ge. bigno) return
!
!     Computes the ratio of energy (B**2) in modes with the undesirable helicities (match)
!     to the energy in modes with the desired helicity (bnorm) and minimizes that...
!

      l_helicity = nint(real(helicity))
      k_helicity = nint(aimag(helicity))
      
      if (nopt .gt. 0) then
         num0 = num+1
         bnorm = zero
         bmax  = maxval(abs(bmn_b(1:mnboz)))
         sj = hs*(real(nsval,rprec) - c1p5)            !!This is correct (SPH)

         do mn = 1,mnboz
            num = num + 1
            n = ixn(mn)/nfp_opt
            m = ixm(mn)
            wegt(num) = bigno
            index_array(num) = ivar_bmn
            chisq_descript(num) = descript(ivar_bmn)
            chisq_target(num) = zero
            chisq_match(num)  = zero
!           if (n.eq.0 .and. m.eq.0) cycle    !NEED FOR CONTINUOUS NORM
            bmn = bmn_b(mn)

!
!           Target for minimization Bmn-s with helicities other than the one desired 
!           General Helical Symmetry: mu - nv ~ Y(lu + kv) for integers Y != 0 (n,k in fp units)
!
            l_symmetry = .false.
            if (k_helicity .eq. 0) then                                  !!quasi-axisymmetry
               if (n .eq. 0) l_symmetry = .true.
            else if (l_helicity .eq. 0) then                             !!quasi-poloidal symmetry
               if (m .eq. 0) l_symmetry = .true.
            else if (mod(m,l_helicity) .eq. 0) then                      !!quasi-helical symmetry (lu + kv)
               if ((m*k_helicity+n*l_helicity).eq.0) l_symmetry = .true.
            endif

            if (l_symmetry) then
               bnorm = bnorm + bmn*bmn
               cycle
            end if

            wegt(num) = one
            chisq_match(num) = bmn

         end do

         bnorm = sqrt(bnorm)
         if (bnorm .eq. zero) bnorm = bmax

         if( sigma < zero) then
             rad_wegt = 1
         else if (m .lt. 3) then           ! apply radial weighting
             rad_wegt = sj
         else if (m .eq. 3) then
             rad_wegt = sj**c1p5
         else
             rad_wegt = sj**2
         end if

         chisq_match(num0:num) = chisq_match(num0:num)/bnorm            !!Norm this way to get correct rms error
         wegt(num0:num) = abs(sigma)*rad_wegt*wegt(num0:num)
      
      else
         num = num + mnboz
      end if   


      end subroutine chisq_bmn
      
      
      subroutine chisq_bootsj (sigma, num, nsurf, nsurf_max, 
     1    nrad, nopt, iflag, extension, lscreen)
      use kind_spec
      use chisq_mod
      use optim_params, ONLY: lbootsj_opt, jboot, lseedcur
      use optim, ONLY: home_dir, jdotb_opt
      use safe_open_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: nrad, nopt, nsurf_max, nsurf(*)      
      integer :: num, iflag
      real(rprec) :: sigma(*)
      character*(*) :: extension
      logical lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: unit_boot=26
      real(rprec), parameter :: zero=0, p5=0.5_dp, one=1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: jrad, iunit, k, nsval, nsin, iwrite
      real(rprec) :: bsnorm, jdotb_norm, jdotb_targ, sjmax, frac_nustar
      character(len=len_trim(home_dir)+20+len_trim(extension)) :: 
     1         version
      logical :: ex
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      real(rprec) , external :: pseed_t, pfboot
C-----------------------------------------------

!
!     BOOTSTRAP CURRENT, <J dot B>, READ IN FROM BOOTSJ FILE
!
      version = trim(home_dir) // '/xbootsj'
      if (nopt. gt. 0) then
!
!     RUN BOOTSTRAP CODE TO CREATE OUTPUT FILE
!
         iunit = unit_boot
   
         call load_physics_codes (version, 'in_bootsj', comm1,
     1         'jBbs', extension, iunit, iflag)
         if (iflag .ne. 0) return

!     READ IN jdotb from xbootsj - produced data file, jBbs
         read(iunit,'(a)', iostat=k)                                     !Discard header

         if (nrad .gt. 0) then
            jdotb_norm = sum(abs(jdotb_opt(:nrad)))/nrad
         else
            jdotb_norm = one
         end if
         if (jdotb_norm .eq. zero) jdotb_norm = one

         iwrite = iunit+1
         version = 'jboot.'//trim(extension)
         call safe_open (iwrite, k, version, 'replace', 'formatted')
         if (k .ne. 0) then
            print *,' Error opening file: ', trim(version)
            iflag = -13
            return
         end if
         write (iwrite, "(2a,/,14x,a,42x,11('-'),/,61x,a)")
     1     '      S    <J*B> (xbootsj)   <J*B> (vmec)',
     2     '    FRAC(NU-STAR)   <J*B>(boot)', '*FRAC','<J*B>(tok)'

         do jrad = 1, nsurf_max
            read (iunit, *, iostat=k) nsin, jdotb_targ, bsnorm
            if (k .ne. 0) then
               iflag = -13
               return
            end if
            nsval = nsurf(jrad)
            if (nsin .ne. nsval) stop 'nsin != nsval in chisq_bootsj'
            num = num + 1
            index_array(num) = ivar_bootsj
            wegt(num) = sigma(nsval)*jdotb_norm

            sjmax = real(nsval - 1,rprec)/(nrad - 1)
            frac_nustar = pfboot(sjmax)
            jdotb_targ = jdotb_targ * frac_nustar                ! mimic nu-star effect: pfboot = 1/(1+nu-star)

            if (jboot.eq.1 .and. abs(bsnorm).gt.zero) 
     1         jdotb_targ = jdotb_targ / bsnorm                  ! use axisymmetric value

            if (lseedcur) 
     1         jdotb_targ = jdotb_targ + pseed_t(sjmax)          ! Seed current (MCZ: Sept 98)

            chisq_target(num) = jdotb_targ                       !!Half radial grid
            chisq_match(num) = p5*(jdotb_opt(nsval)
     1                       +     jdotb_opt(nsval-1))           !!VMEC on full grid
            chisq_descript(num) = descript(ivar_bootsj)
            write(iwrite,'(f8.2,1pe15.2,2x,3e15.2,6x,3e15.2,
     1                     6x,3e15.2)') sjmax, chisq_target(num), 
     2                     chisq_match(num), frac_nustar, bsnorm
         end do

         close (iwrite)
         if (lscreen) call system ("cat -s " // trim(version))
         close(iunit, status='delete')                           !!Remove jBbs.ext file
         version = "answers_plot" // "." // trim(extension)
         inquire (file=version, exist=ex)
         if (ex) call system ("rm -f " // trim(version))

      else
        inquire(file=trim(version), exist=ex)
        if (.not.ex .and. lscreen) then
           print *,'xbootsj file not found in ' // trim(home_dir)
           lbootsj_opt = .false.
        else if (nsurf_max .ge. 1) then
           num = num + nsurf_max
        end if
      endif

      end subroutine chisq_bootsj
      

      subroutine chisq_kink (num, nopt, iflag, 
     1    extension, lscreen)
      use kind_spec
      use chisq_mod
      use optim, ONLY: home_dir, chisq_min, bigno
      use optim_params, ONLY: lkink_opt, sigma=>sigma_kink,
     1                        target=>target_kink
      use system_mod, ONLY: system
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: nopt, iflag, num
      character*(*) :: extension
      logical :: lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, parameter :: unit_kink = 27
      integer :: iunit, ifail, istat, n, i
      real(rprec), parameter :: zero = 0
      real(rprec) :: growth_rate
      character(len=len_trim(home_dir)+100) :: version, ft5_ver
      logical :: ex
      character*1 :: tag
C-----------------------------------------------
      version = trim(home_dir) // '/xtprp'
      iflag = 0

      if (nopt > 0) then
        do n=1, size(sigma)
          if( sigma(n) >= bigno) cycle
          if( n == 1) then
             tag = ' '
          else
             write(tag,'(i1)') n
          endif

          iunit = unit_kink
          iflag = 0
            
          if (lscreen) write(6,*)'Running TERPSICHORE('//tag//
     1                 ') for kink stability'

          do i=1, 3   ! we will try 3 times to get a terpsichore output
             call load_physics_codes(trim(version)//tag, " ", "-b " // 
     1         trim(extension),
     1         'ft79tpr'//tag, extension, iunit, iflag)

             if( iflag .ne. -2 ) exit      ! success or bad failure
             write (6,*) 'No Terpsichore output file!'
          enddo

          if (iflag .ne. 0 .and. iflag .ne. -2) then
             write(6,*) 'Problems with load_physics(terp.), iflag=',
     1                   iflag
             return

          else if (iflag .ne. -2) then

             ifail = 0
             read(iunit, *, iostat=istat) ifail, growth_rate

             if (istat.ne.0 .or. ifail.eq.1) then
c              if (lscreen)  
                 write(6,*) 'Error running Terpsichore'//tag//
     1               '/kink, status=',istat, ' Ifail = ', ifail
               if( ifail .eq. 1 .or. istat.eq.-4001) then

c        Terpsichore could not construct vacuum region. 
c                 if (lscreen) 
                  write(6,*) 'Terpsichore'//tag//
     1                  ' error construction vacuum, ',
     1                  'Deprecate this direction!'

                  growth_rate = -100*sqrt(chisq_min*sigma(n))
                  iflag = 0
               else
                  iflag = -19
                  return
               endif
             endif

          else
             write(6,*) 'Terpsichore'//tag//
     1                  ' error, never wrote output files. ',
     1                  'Deprecate this direction!'

             growth_rate = -100*sqrt(chisq_min*sigma(n))
             iflag = 0
             
          endif

          num = num + 1
          index_array(num) = ivar_kink
          chisq_match(num) = 0

          chisq_target(num) = target(n)
!         chisq_target(num) = sigma/2

          chisq_match(num) = min (chisq_target(num), growth_rate)

          close(iunit, status='delete')

          chisq_descript(num) = descript(ivar_kink)
          wegt(num) = sigma(n)
        enddo
      else
        do n=1, size(sigma)
          if( sigma(n) >= bigno) cycle
          if( n == 1) then
             tag = ' '
          else
             write(tag,'(i1)') n
          endif

          inquire(file=trim(version)//tag, exist=ex)
          if (.not.ex) then
            if(lscreen) then
              print *, 'xtprp'//tag//' file not found in '// 
     1                trim(home_dir)
              print *, '*** Kink evaluation disabled !'
            endif
            lkink_opt = .false.
          else
            inquire(file='../ft5tpr'//tag, exist=ex)
            if (.not.ex) then
              if(lscreen) then
                print *,' TERPSICHORE control file ft5tpr'//tag//
     1                 ' is missing'
                print *, '*** Kink evaluation disabled !'
              endif
              lkink_opt = .false.
            else   
               ifail = 0
               ft5_ver = "/bin/ln -s ../ft5tpr"//tag//" ./ft5tpr"
     1                    //tag
               call system(trim(ft5_ver), ifail)
               if (ifail.lt.127 .and. ifail.ne.0) then
                 if(lscreen) then
                   print *,' TERPSICHORE control file ft5tpr'
     1                   //tag//' unlinked, status=', ifail
                   print *, '*** Kink evaluation disabled !'
                 endif
                 lkink_opt = .false.
               else   
                  num = num + 1
               end if   
            end if 
          end if 
        enddo
      endif

      end subroutine chisq_kink
      
      
      subroutine chisq_dsubr (sigma, ivar, num,   
     1                      nopt, iflag, extension, lscreen)
      use kind_spec
      use chisq_mod
      use optim, ONLY: home_dir, nrad
      use optim_params, ONLY: ldsubr_opt
      use system_mod, ONLY: system
      use safe_open_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: ivar, nopt, iflag, num
      real(rprec), dimension(*) :: sigma
      character*(*) :: extension
      logical :: lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, parameter :: unit_jmc = 79
      integer :: iunit, ierr, istat, i, ns_jmc_max
      integer :: nsval
      real(rprec), parameter :: zero = 0
      real(rprec), dimension(nrad)  :: dsubr 
      integer, dimension(nrad) :: ns_jmc
      character*200 :: version
      logical :: ex
C-----------------------------------------------
      version = trim(home_dir) // '/xjmc'
      
      if (nopt > 0) then

         iunit = unit_jmc
            
         if (lscreen) then
           write(6,*) 
     1           'Running JMC for D_sub_r calculations'
         endif

         ierr = 0
         version = trim(version) // ' ' // extension
         call system(version,ierr)                              
         if( ierr .lt. 127 .and. ierr .ne. 0 ) then
             print *, "Trouble running JMC, ierr = ", ierr
             return
         endif

         call safe_open(iunit, istat, 'ft79jmc.'//trim(extension), 
     1                  'old', 'formatted')
         if( istat .ne. 0 ) then
            iflag = -32
            return
         else
            read(unit_jmc,*) ns_jmc_max
            read(unit_jmc,*) (ns_jmc(i),dsubr(i),i=1,ns_jmc_max)
            if( ns_jmc_max .ne. nrad - 3 ) stop "Error in JMC"
            close(unit_jmc)
         endif
!
!         calculated on all surfaces from 2 to nrad-2
!
          do i=1, ns_jmc_max
             num = num + 1
             index_array(num) = ivar
             chisq_descript(num) = descript(ivar)
             nsval = ns_jmc(i)
             chisq_match(num) =  -dsubr(i)
             chisq_target(num) =  min(chisq_match(num), zero)
             wegt(num) = sigma(nsval)
          enddo

      else
         inquire(file=version, exist=ex)
         if (.not.ex) then
            print *, 'xjmc file not found in ' // trim(home_dir)
            ldsubr_opt = .false.
         else
           num = num + nrad - 3
         end if 
      endif

      end subroutine chisq_dsubr


      subroutine chisq_nescoil (number, nopt, iflag, extension, lscreen)
      use kind_spec
      use chisq_mod
      use optim, ONLY: home_dir, bigno
      use optim_params, ONLY: lnescoil_opt, coil_separation, 
     1     target_coil_jmax, sigma_coil_complex, sigma_coil_jmax,
     2     sigma_berr_ave, target_coil_complex
      use system_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: number, nopt, iflag 
      character*(*) :: extension
      logical :: lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: nescoil0 = 28
      real(rprec), parameter :: zero = 0, one = 1
      integer, parameter :: ivals = 3
      character*10, parameter :: search_string(ivals) =
     1   (/ 'Complexity', 'J Surf max', 'Berr ave, ' /)
      integer, parameter :: search_index(ivals) =  (/ 14, 24, 22 /)
      real(rprec) :: search_value(ivals)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: iunit, istat, index1, itype, m, n
      real(rprec) :: coil_complex_opt, jsheet_max, berr_ave, berr_max
      character(len=len_trim(home_dir)+len_trim(extension)+30) 
     1    :: version
      character*200 :: line
      logical :: ex, ex1
C-----------------------------------------------
!
!     COMPUTE COIL COMPLEXITY FROM NESCOIL POTCOEFFS FILE
!
      version = trim(home_dir) // '/xbnorm'

      if (nopt .gt. 0) then
         iunit = nescoil0

            
!     Write coil_separation for use by BNORM (NESCOIL)
         version = trim(version) // ' wout.' // trim(extension)
         write (version,'(a,1pe12.3)') trim(version), coil_separation
         if (lscreen) then
           print '(/,a)',' Running BNORM code...'
           print 100,' Coil separation = ', coil_separation
         end if   
         call system(version)                              !!Produce bnorm, nescin files

         version = trim(home_dir) // '/xnescoil'  
         if (lscreen) 
     1      print '(/,a)',' Running NESCOIL code ...'
         call load_physics_codes (version, 'nescin', ' ',
     1        'nescout', extension, iunit, iflag)
         if (iflag .ne. 0) return
!
!     Read in current complexity,amaximum current density, and berr from nescout file (iunit)
!
         do itype = 1, ivals
            index1 = -1
            do while (index1 .le. 0)
               read(iunit, '(a)', iostat=istat) line
               if (istat .ne. 0) stop ' ISTAT != 0 in chisq_nescoil'
               index1 = index(line, search_string(itype))
               if (index1 .gt. 0) then
                  if (itype .eq. 3) then
                     read (line(search_index(itype):), *)
     1               search_value(itype), berr_max
                  else 
                     read (line(search_index(itype):), *) 
     1               search_value(itype)
                  end if
               end if
            end do
         end do

         coil_complex_opt = search_value(1)
         jsheet_max       = search_value(2)
         berr_ave         = search_value(3)
         
         if (lscreen) then
            print 100,' Coil Complexity = ', coil_complex_opt
            print 100,' Maximum Sheet Current Density = ', jsheet_max
            print 100,' Berr-ave = ', berr_ave
         end if

         close (iunit)

!        COIL COMPLEXITY (NO PENALTY IF ACTUAL COMPLEXITY IS BELOW TARGET)
         if (sigma_coil_complex .lt. bigno) then
            number = number + 1
            index_array(number) = ivar_coil_complexity
            wegt(number) = sigma_coil_complex
            chisq_target(number) = target_coil_complex
            chisq_match(number)  = 
     1         max(target_coil_complex, coil_complex_opt)
            chisq_descript(number) = descript(ivar_coil_complexity)
         end if

!        MAXIMUM CURRENT DENSITY
         if (sigma_coil_jmax .lt. bigno) then
            number = number + 1
            index_array(number) = ivar_coil_jmax
            wegt(number) = sigma_coil_jmax
            chisq_target(number) = min(abs(jsheet_max),target_coil_jmax)
            chisq_match(number)  = abs(jsheet_max)
            chisq_descript(number) = descript(ivar_coil_jmax)
         end if

!        AVERAGE BERR
         if (sigma_berr_ave .lt. bigno) then
            number = number + 1
            index_array(number) = ivar_berr_ave
            wegt(number) = sigma_berr_ave
            chisq_target(number) = zero
            chisq_match(number)  = abs(berr_ave)
            chisq_descript(number) = descript(ivar_berr_ave)
         end if

      else
         inquire(file=trim(version), exist=ex)
         inquire(file=trim(home_dir) // '/xnescoil', exist=ex1)
         if (.not.ex) then
            if (lscreen)
     1          print *, 'xbnorm file not found in ' // trim(home_dir)
            lnescoil_opt = .false.
         else if (.not.ex1) then
            if (lscreen)
     1          print *, 'xnescoil file not found in ' // trim(home_dir)
            lnescoil_opt = .false.
         else   
            if (sigma_coil_complex .lt. bigno) number = number+1
            if (sigma_coil_jmax .lt. bigno)    number = number+1
            if (sigma_berr_ave .lt. bigno)     number = number+1
         end if   
      endif

 100  format(a,1pe12.3)

      end subroutine chisq_nescoil
      
      
      subroutine chisq_curvature (sigma, num, nopt)
      use kind_spec
      use chisq_mod
      use boozer_params
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: nopt
      integer :: num
      real(rprec) :: sigma
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: test_pts = 4
      real(rprec), parameter :: p5 = 0.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: lzth
      real(rprec) :: test_fcn
      real(rprec), dimension(4) :: curv_kur      
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      external flux_surf_curv
C-----------------------------------------------

      if (nopt .gt. 0) then
         call flux_surf_curv (curv_kur, rmnc_bdy, zmns_bdy, 
     1        xm_bdy, xn_bdy)
         do lzth = 1, test_pts
            num = num + 1
            index_array(num) = ivar_curve
c
c       As a target to prevent flux surface cusping, use a tanh functional of
c       the curvature kurtosis around 0,90,180, and 270 degree phi planes.
c       This is designed to allow moderate elliptical and triangular elongations,
c       but not strong cusps. Note: the numbers in this functional may need
c       readjustment for particular situations.
c
            test_fcn = 0.44_dp + p5*
     1        tanh((curv_kur(lzth)-20.0_dp)/15.0_dp)
            wegt(num) = sigma
            chisq_target(num) = 1.e-3_dp
            chisq_match(num) = test_fcn
            chisq_descript(num) = descript(ivar_curve)
         end do
      else
         num = num + test_pts
      end if   


      end subroutine chisq_curvature
      
            
      subroutine chisq_iota(match, sigma, hs, num, nrad, nopt)
      use kind_spec
      use chisq_mod
      use optim, ONLY: bigno
      use optim_params, ONLY: 
     1    target_iota_max, target_iota_min, target_iota_max_min,
     2    sigma_iota_max, sigma_iota_min, sigma_iota_max_min
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: nrad, nopt
      integer :: num
      real(rprec) :: match(*), sigma(*), hs
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, c1p5 = 1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: jrad
      real(rprec) :: TargetIota, AvgIota, sj
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      real(rprec) , external :: piota_t
C-----------------------------------------------
      if (any(abs(sigma(2:nrad)) .lt. bigno)) then
        if (nopt .gt. 0) then
          avgiota = zero
          do jrad = 2, nrad
            sj = hs*(jrad - c1p5)
            AvgIota = AvgIota + abs(piota_t(sj))
          end do
          AvgIota = AvgIota/(nrad - 1)
          do jrad = 2, nrad
            num = num + 1
            index_array(num) = ivar_iota
            sj = hs*(jrad - c1p5)
            TargetIota = piota_t(sj)
            wegt(num) = sigma(jrad)*AvgIota
            chisq_target(num) = TargetIota
            chisq_match(num) = match(jrad)
            chisq_descript(num) = descript(ivar_iota)
          end do

        else
          if (nrad .gt. 1) num = num + nrad - 1 
        endif
      endif

      if( sigma_iota_max < bigno) then
         num = num+1
         if( nopt > 0) then
            index_array(num) = ivar_iota_bounds
            chisq_descript(num) = descript(ivar_iota_bounds)
            wegt(num) = sigma_iota_max
            chisq_target(num) = target_iota_max
            chisq_match(num) = max(target_iota_max, 
     1                             maxval(match(2:nrad)))
         endif
      endif

      if( sigma_iota_max_min < bigno) then
         num = num+1
         if( nopt > 0) then
            index_array(num) = ivar_iota_bounds
            chisq_descript(num) = descript(ivar_iota_bounds)
            wegt(num) = sigma_iota_max_min
            chisq_target(num) = target_iota_max_min
            chisq_match(num) = min(target_iota_max_min, 
     1                             maxval(match(2:nrad)))
         endif
      endif

      if( sigma_iota_min < bigno) then
         num = num+1
         if( nopt > 0) then
            index_array(num) = ivar_iota_bounds
            chisq_descript(num) = descript(ivar_iota_bounds)
            wegt(num) = sigma_iota_min
            chisq_target(num) = target_iota_min
            chisq_match(num) = min(target_iota_min, 
     1                             minval(match(2:nrad)))
         endif
      endif
      end subroutine chisq_iota
            

      subroutine chisq_iota_p(iota, sigma_max, sigma_min, hs, 
     1           num, nrad, nopt)
      use kind_spec
      use chisq_mod
      use optim, ONLY: bigno
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: nrad, nopt
      integer :: num
      real(rprec), intent(in) :: iota(nrad), sigma_max(nrad), 
     1    sigma_min(nrad), hs
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(kind=rprec), parameter :: zero = 0, c1p5 = 1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: jrad, ncount
      real(kind=rprec) :: Target, sj, match
     
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      real(rprec) , external :: piota_prime
C-----------------------------------------------
      ncount = count(abs(sigma_max(2:nrad-1)) < bigno)
     1       + count(abs(sigma_min(2:nrad-1)) < bigno)

C      print *,'iota_p ncount=',ncount

      if (ncount == 0) return

      if (nopt .gt. 0) then
         do jrad = 2, nrad-1
           if( sigma_max(jrad)<bigno .or. sigma_min(jrad)<bigno) then
             sj = hs*(jrad - c1p5)

             Target = piota_prime(sj)
             match = (iota(jrad+1)-iota(jrad-1))/(2*hs)

c             print *,'iota_p, request, sig_mx, sig_mn=',match,target,
c     1             sigma_max(jrad), sigma_min(jrad)
           endif

           if( sigma_max(jrad)<bigno ) then
             num = num + 1
             index_array(num) = ivar_iota_p
             wegt(num) = sigma_max(jrad)
             chisq_target(num) = Target
             chisq_match(num) = max(match, Target)
             chisq_descript(num) = descript(ivar_iota_p)
           endif

           if( sigma_min(jrad)<bigno ) then
             num = num + 1
             index_array(num) = ivar_iota_p
             wegt(num) = sigma_min(jrad)
             chisq_target(num) = Target
             chisq_match(num) = min(match, Target)
             chisq_descript(num) = descript(ivar_iota_p)
           endif

         end do
      else
         if (nrad .gt. 1) then
            num = num + ncount
         endif
      endif

      end subroutine chisq_iota_p
            

      subroutine chisq_jstar (sigma, bmax_b, bmin_b,
     1          bmod_b, num, nopt, mskip, sqrt_nsurf)
      use kind_spec
      use chisq_mod
      use boozer_params, ONLY: nu2_b, nv_boz
      use vparams, ONLY: twopi
      use optim, ONLY: bigno
      use optim_params, ONLY: NumJstar
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: nopt, mskip
      integer :: num
      real(rprec), intent(in) :: bmax_b(*), bmin_b(*), bmod_b(*), 
     1    sigma(1,*), sqrt_nsurf
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: p5 = 0.5_dp, zero = 0, one = 1
      real(rprec), parameter :: max_variation = 2.e-2_dp
      real(rprec), parameter :: epl = 0.05_dp, epu = 0.05_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: iep, nJstar, n
      real(rprec), dimension(:), allocatable :: trapJS      
      real(rprec) :: bmin_global, bmax_global, epsmu, sj, avg_Jstar
      logical, dimension(:), allocatable  :: ljstar
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      external j_star
C-----------------------------------------------
      if (all(abs(sigma(1,1:NumJstar)) .ge. bigno)) return

!        CODE GENERATED BY D. SPONG, 6/97
!        COMPUTE J_star AT NUMJSTAR VALUES OF ep/mu RANGING FROM SLIGHTLY ABOVE
!        THE TRAPPED-PASSING BOUNDARY TO SLIGHTLY BELOW THE
!        DEEPLY TRAPPED-FORBIDDEN BOUNDARY.  THE PARAMETERS epl AND epu
!        DETERMINE DISTANCE TO THESE BOUNDARIES.

      if (nopt .gt. 0) then
         allocate (ljstar(nu2_b), trapJS(nu2_b))

         bmin_global = minval(bmin_b(:nu2_b))
         bmax_global = maxval(bmax_b(:nu2_b))

         do iep = 1,NumJstar

            if (abs(sigma(1,iep)) .ge. bigno) cycle

            sj = real(iep - 1,rprec)/(NumJstar - 1)
            epsmu = bmin_global*(one + epl) + sj*(
     1         bmax_global*(one - epu) - bmin_global*(one + epl))
            call j_star (bmod_b, bmin_b, bmax_b, epsmu, trapJs, 
     1          nv_boz, nu2_b)

            ljstar = trapJs(:nu2_b) .gt. zero
            nJstar = count(ljstar)
            avg_Jstar = sum(trapJs, mask=ljstar)
            if (nJstar .gt. 0) avg_Jstar = avg_Jstar/nJstar
            where (.not.ljstar) trapJs = avg_Jstar
!
!          Target all non-zero Jstars to (d Jstar/du) = 0
!
           do n = 1, nu2_b, mskip
              num = num+1
              index_array(num) = ivar_jstar
              wegt(num) = max_variation * avg_Jstar
     1                   * sqrt_nsurf * sigma(1,iep)
              chisq_target(num) = avg_Jstar
              chisq_descript(num) = descript(ivar_jstar)
              chisq_match(num) = trapJs(n)
            end do
         end do

         deallocate (ljstar, trapJS)

      else                   
         do iep = 1,NumJstar
            if (abs(sigma(1,iep)) .ge. bigno) cycle
            do n = 1, nu2_b, mskip
               num = num + 1
            end do   
         end do   
      end if   

      end subroutine chisq_jstar


      subroutine chisq_jinvar (sigma, nrad, num, nopt, nu,
     1    sqrt_nsurf, iflag, extension, command, lscreen)
      use kind_spec
      use chisq_mod
      use optim, ONLY: bigno, home_dir
      use optim_params, ONLY: NumJinvariant, lj_invariant
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: nopt, nu, nrad
      integer :: num, iflag
      character *(*) :: extension, command
      real(rprec) :: sigma(1,NumJinvariant), sqrt_nsurf
      logical, intent(in) :: lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer :: unit_jinvar = 30
      character*20 :: temp
      real(rprec), parameter :: p5 = 0.5_dp, zero = 0, one = 1
      real(rprec), parameter :: max_variation = 2.e-2_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: iunit, ieps, n, istat, NJinvar
      character(len=len_trim(home_dir)+20) :: version
      logical :: ex
      logical, allocatable, dimension(:) :: ljinvar
      real(rprec), dimension(:), allocatable :: J_Invariant
      real(rprec) :: avg_Jinvar
      real(rprec) :: sqrt_nusurf
C-----------------------------------------------
      if (all(abs(sigma(1,1:NumJinvariant)) .ge. bigno)) return
!
!        COMPUTE J-INVARIANT AT NUMJINVARIANT VALUES OF ep/mu RANGING FROM SLIGHTLY ABOVE
!        THE TRAPPED-PASSING BOUNDARY TO SLIGHTLY BELOW THE
!        DEEPLY TRAPPED-FORBIDDEN BOUNDARY.  THE PARAMETERS epl AND epu
!        DETERMINE DISTANCE TO THESE BOUNDARIES.
!
      version = trim(home_dir) // '/xj_invariant'

      if (nopt .gt. 0) then
!
!     RUN J-INVARIANT CODE TO CREATE OUTPUT FILE
!
         iunit = unit_jinvar
   
         write (temp,'(1x,i3,1x,i3,1x,i3,a2)') 
     1       nrad, NumJinvariant, nu, command  
     
         call load_physics_codes (version, 'boozmn', temp,
     1         'j_invar_out', extension, iunit, iflag)
         if (iflag .ne. 0) return

         if (nu .le. 0) return
         sqrt_nusurf = sqrt(real(nu, rprec))
         allocate (J_Invariant(nu), ljinvar(nu))
!
!        Target all non-zero Jinvars to (d Jinvar/du) = 0
!
         do ieps = 1, NumJinvariant
            read (iunit, *, iostat=istat) (J_Invariant(n), n = 1, nu)
            if (istat .ne. 0) then
               iflag = -18
               return
            end if   

            ljinvar = J_Invariant .gt. zero
            avg_Jinvar = sum(J_Invariant, mask=ljinvar)
            nJinvar = count(ljinvar)
            if (nJInvar .gt. 0) avg_Jinvar = avg_Jinvar/nJinvar
            where (.not.ljinvar) J_Invariant = avg_Jinvar

            if (abs(sigma(1,ieps)) .ge. bigno) cycle
            do n = 1, nu                                 !!Read in nu from file...
              num = num+1
              index_array(num) = ivar_jinvar
              wegt(num) = max_variation * sqrt_nsurf * avg_Jinvar
     1                   * sigma(1,ieps)
              chisq_target(num) = avg_Jinvar
              chisq_descript(num) = descript(ivar_jinvar)
              chisq_match(num) = J_Invariant(n)
            end do
         end do

         close (iunit)                                   !!Opened in call to load_physics....
         deallocate (J_Invariant, ljinvar)
      else                   
         inquire(file=trim(version), exist=ex)
         if (.not.ex) then
            if (lscreen) print *,    
     1          'xj_invariant file not found in ' // trim(home_dir)
            lj_invariant = .false.
         else 
            do ieps = 1, NumJinvariant
               if (abs(sigma(1,ieps)) .lt. bigno) num = num + nu
            end do   
         end if   
      end if   

      end subroutine chisq_jinvar


      subroutine chisq_jconf(sigma, nrad, num, nopt, nu,
     1    iflag, extension, command, lscreen)
c
c   Added  Aug. 2002  M. Zarnstorff
c   Targets confinement of J-contours, starting from a specified sourcesurface,
c   Penalize the fraction of the J-range for the trapped particles that make it
c   to outer surfaces, weighted by sigma_conf
c
      use kind_spec
      use chisq_mod
      use optim, ONLY: bigno, home_dir
      use optim_params, ONLY: NS_JConf_Src, NS_JConf_tgt, NPitch_Jconf
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: nopt, nu, nrad
      integer :: num, iflag
      character *(*) :: extension, command
      real(rprec) :: sigma
      logical, intent(in) :: lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer :: unit_jinvar = 30
      character*120 :: temp
      real(rprec), parameter :: p5 = 0.5_dp, zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: iunit, ieps, n, istat, NJinvar, ns, is, nu2
      character(len=len_trim(home_dir)+20) :: version
      logical :: ex
      real(rprec) :: epmu, Jmin, Jmax, fract, dum1, dum2, avg

C-----------------------------------------------
      if (abs(sigma) .ge. bigno .or. NPitch_Jconf .le. 0 .or.
     1    nu .le. 0 .or. NS_JConf_Src .lt. 2 .or.
     1    NS_JConf_Src >= min(NS_Jconf_Tgt, nrad) ) return
!
!        COMPUTE J-INVARIANT AT NUMJINVARIANT VALUES OF ep/mu RANGING FROM SLIGHTLY ABOVE
!        THE TRAPPED-PASSING BOUNDARY TO SLIGHTLY BELOW THE
!        DEEPLY TRAPPED-FORBIDDEN BOUNDARY.  THE PARAMETERS epl AND epu
!        DETERMINE DISTANCE TO THESE BOUNDARIES.
!
      version = trim(home_dir) // '/xj_invariant'
!      nu2 = max(nu, NPitch_Jconf)
      nu2 = nu

      if (nopt .gt. 0) then
!
!     RUN J-INVARIANT CODE on Target surface to get range of Jinvariant
!
         iunit = unit_jinvar
   
         write (temp,'(1x,i3,1x,i3,1x,i3,a2,i4)') 
     1       NS_JConf_Src, NPitch_Jconf, nu2, command,
     2       min(nrad, NS_JConf_Tgt)
     
         call load_physics_codes (version, 'boozmn', temp,
     1         'j_invar_sum', extension, iunit, iflag)
         if (iflag .ne. 0) return

         do ieps = 1, NPitch_Jconf
            read (iunit, *, iostat=istat) fract
            if (istat .ne. 0) then
               iflag = -18
               return
            end if   

            num = num+1
            index_array(num) = ivar_jconf
            wegt(num) = sigma * sqrt(real(NPitch_Jconf,rprec))
            chisq_target(num) = 0
            chisq_descript(num) = descript(ivar_jconf)
            chisq_match(num) = sqrt(fract)
         end do

         close (iunit)                     !!Opened in call to load_physics....

      else                   
         inquire(file=trim(version), exist=ex)
         if (.not.ex) then
            if (lscreen) print *,    
     1          'xj_invariant file not found in ' // trim(home_dir)
            stop
         else 
            num = num + NPitch_Jconf
         end if   
      end if   

      end subroutine chisq_jconf


      subroutine chisq_maxcurrent(target, sigma, num, nrad, nopt)
      use kind_spec
      use chisq_mod
      use optim, ONLY: buco_opt, bigno
      use vparams, ONLY: dmu0, twopi
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: nrad, nopt
      integer :: num
      real(rprec), intent(in) :: target, sigma
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------  
      integer :: jrad, jcount, i
      real(rprec) :: dnorm, dmax_j
C-----------------------------------------------
      if (abs(sigma) .ge. bigno) return

      jcount = 0

      if (nopt .gt. 0) then
        do jrad = 2,nrad
          num = num + 1
          index_array(num) = ivar_maxj
          dmax_j = twopi * buco_opt(jrad) / dmu0
          wegt(num) = abs(sigma*target)
          chisq_match(num) = dmax_j
          chisq_descript(num) = descript(ivar_maxj)
          if (abs(dmax_j) .gt. abs(target) + wegt(num)) then
            chisq_target(num) = sign(target, dmax_j)
            jcount = jcount + 1
          else
            chisq_target(num) = chisq_match(num)
          end if  
        end do 
      else
          if (nrad .gt. 1) num = num + nrad - 1
      endif

 
      if (jcount .gt. 1) then
        i = num - (nrad-2)
        dnorm = sqrt(real(jcount,rprec))
        wegt(i:num) = wegt(i:num) * dnorm
      end if  

      end subroutine chisq_maxcurrent
      
      
      subroutine chisq_polflux (target, sigma, hs, nrad, num, nopt)
      use kind_spec
      use chisq_mod
      use optim, only: iota_opt, phip_opt
      use vparams, ONLY: twopi
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: nrad, nopt
      integer, intent(inout) :: num
      real(rprec), intent(in) :: target, sigma, hs
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec) :: psi_a
C-----------------------------------------------
      num = num + 1
      if (nopt .gt. 0)then
         psi_a = abs(-twopi * hs * 
     1               sum(iota_opt(2:nrad)*phip_opt(2:nrad)) )

         index_array(num) = ivar_fluxp
         chisq_descript(num) = descript(ivar_fluxp)
         chisq_match(num) = psi_a
         chisq_target(num) = max(psi_a, target)
         wegt(num) = sigma * abs(target)
      end if

      end subroutine chisq_polflux
      
      
      subroutine chisq_mercier(match, sigma,  
     1           num, nrad, nopt)
      use kind_spec
      use chisq_mod
      use optim, ONLY: bigno
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: nrad, nopt
      integer :: num
      real(rprec) :: match(*), sigma(*)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0.0_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: jrad
      integer, dimension(nrad):: mercier_flag, pre_flag, post_flag
C-----------------------------------------------
      if (all(abs(sigma(3:nrad-1)) .ge. bigno)) return

!
!     LOGIC TO AVOID ISOLATED RESONANCES
!
      mercier_flag = 0; pre_flag = 0; post_flag = 0                            !! MERCIERFLAG (RS)
      where(match(3:nrad-1)< zero) mercier_flag(3:nrad-1) = 1                  !! MERCIERFLAG (RS)
      pre_flag(2:nrad-2) = mercier_flag(3:nrad-1)                              !! MERCIERFLAG (RS)
      post_flag(4:nrad) = mercier_flag(3:nrad-1)                               !! MERCIERFLAG (RS)
      mercier_flag = mercier_flag + pre_flag + post_flag                       !! MERCIERFLAG (RS)

      if (nopt. gt. 0) then
         do jrad = 3, nrad-1
            num = num + 1
            index_array(num) = ivar_mercier
            wegt(num) = sigma(jrad)
            chisq_target(num) = zero
            if(mercier_flag(jrad) > 1) then                                    !! MERCIERFLAG (RS)
              chisq_match(num) = min(match(jrad), zero)
            else
              chisq_match(num) = zero                                          !! MERCIERFLAG (RS)
            endif
            chisq_descript(num) = descript(ivar_mercier)
         end do
      else
         if (nrad .gt. 3) num = num + nrad - 3
      endif

      end subroutine chisq_mercier
      

      subroutine chisq_magwell(hs, sigma, num, nrad, nopt)
      use kind_spec
      use chisq_mod
      use optim, ONLY: bigno, vp_opt, phip_opt
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: nrad, nopt
      integer :: num
      real(rprec) :: hs, sigma(*)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0.0_dp, one = 1.0_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: jrad
      real(rprec) :: TargetWell, AvgWell, vpf(nrad), sj, vpp
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      real(rprec) , external :: pwell_t
C-----------------------------------------------
      if (nopt .gt. 0) then

         vpf(1) = 1.5_dp*vp_opt(2)/phip_opt(2)
     1          - 0.5_dp*vp_opt(3)/phip_opt(3)
         vpf(nrad) = 1.5_dp*vp_opt(nrad)/phip_opt(nrad) -
     1               0.5_dp*vp_opt(nrad-1)/phip_opt(nrad-1)
         vpf(2:nrad-1) = 0.5_dp*(vp_opt(3:nrad)/phip_opt(3:nrad) +
     1                           vp_opt(2:nrad-1)/phip_opt(2:nrad-1))

         AvgWell = 0._dp
         do jrad = 2, nrad
           sj = hs * (jrad - one)
           AvgWell = AvgWell + abs(pwell_t(sj)) * hs
         end do
         do jrad = 2, nrad
            num = num + 1
            index_array(num) = ivar_well
            vpp = (vpf(jrad)-vpf(1))/vpf(1)
            if( AvgWell .lt. 1.0e-20_dp )
     1          vpp = (vpf(jrad)-vpf(jrad-1))/vpf(1)
            sj = hs*(jrad - one)
            TargetWell = pwell_t(sj)
            if( AvgWell .lt. 1.0e-20_dp .and. vpp .lt. 0.0_dp)
     1         TargetWell = vpp
            wegt(num) = sigma(jrad)
            if(AvgWell .ge. 1.0e-20_dp)
     1            wegt(num) = sigma(jrad)*avgwell
            chisq_target(num) = targetwell
            chisq_match(num) = vpp
            chisq_descript(num) = descript(ivar_well)
         end do
      else
         if (nrad .gt. 1) num = num + nrad - 1
      endif

      end subroutine chisq_magwell
      

      subroutine chisq_rzbdy(target, x_opt, sigma, ivar, num, nopt, 
     1    lmin)
      use kind_spec
      use chisq_mod
      use optim, ONLY: bigno
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: ivar, nopt
      integer :: num
      real(rprec) :: target, x_opt, sigma
      logical :: lmin
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0.0_dp, one = 1.0_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: jrad
C-----------------------------------------------
      if (sigma .ge. bigno) return
      
      num = num + 1

      if (nopt .gt. 0) then
         index_array(num) = ivar
         wegt(num) = sigma*Target
         if (wegt(num) .eq. zero) wegt(num) = one
         chisq_target(num) = Target
         if (lmin) then
            chisq_match(num) = max(Target, x_opt)
         else
            chisq_match(num) = min(Target, x_opt)
         end if   
         chisq_descript(num) = descript(ivar)
      end if   

      end subroutine chisq_rzbdy


      subroutine chisq_bpres (pres_opt, sigma, ivar, num, nrad, 
     1   nopt, nflag)
      use kind_spec
      use chisq_mod
      use optim, ONLY: home_dir
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: ivar, nrad, nopt, nflag
      integer :: num
      real(rprec) :: pres_opt(*), sigma(*)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: jrad
      real(rprec) :: balpgrad_opt
C-----------------------------------------------
!
!     MATCH EDGE PRESSURE=0
!

      if (nflag .eq. 0) then
         num = num + 1
         if (nopt .gt. 0) then
            index_array(num) = ivar
            wegt(num) = sigma(1)
            chisq_target(num) = zero
            chisq_match(num) = pres_opt(nrad)
            chisq_descript(num) = descript(ivar)
      end if

      else
!
!     MATCH MONOTONICALLY DECREASING PRESSURE PROFILE  (dP/ds<0)
!
         if (nopt. gt. 0) then
            do jrad = 3, nrad
               num = num + 1
               index_array(num) = ivar
               wegt(num) = sigma(jrad)
               balpgrad_opt = (pres_opt(jrad) - pres_opt(jrad-1))
     1             /real(nrad-1, rprec)
               chisq_target(num) = zero
               chisq_match(num) = max(balpgrad_opt, zero)
               chisq_descript(num) = descript(ivar)
            end do
         else
            if (nrad .gt. 2) num = num + nrad - 2
         endif
      end if   

      end subroutine chisq_bpres


      subroutine chisq_p_prof (pres_opt, ivar, num, nrad, nopt, 
     1                         extension)
      use kind_spec
      use chisq_mod
      use optim_params
      use safe_open_mod
      use AJAX_MOD
      use vmec_input, ONLY: am
      use optim, only: lajax, bigno
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: ivar, nrad, nopt
      integer :: num
      real(rprec) :: pres_opt(*)
      character*(*) :: extension
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, j, k, num0, iflag, navg, iunit
      real(rprec) :: s, avg_data, avg_vmec, sj, fract, factor, match,
     1               chisq_fit, int_am, int_fit_am, p_fit
      real(rprec), dimension(3) :: r_cyl, r_flx
      real(rprec), dimension(np_prof) :: s_prof, sigma_fit
      real(rprec), dimension(1:11) :: fit_am
      character*120 :: message
      real(rprec) :: pi
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      real(rprec) , external :: eval_prof
C-----------------------------------------------

      pi = 4*atan(1._rprec)

      if (nopt > 0) then

         iunit = unit_outdata
         call safe_open(iunit, k, 'p_prof.'//trim(extension),
     1          'replace', 'formatted')
         if (k .ne. 0) then
            iflag = -12
            return
         endif

         factor = 1
         if( (abs(sigma_beta).lt.bigno .or. abs(sigma_eplasma) < bigno)
     1                   ) then               ! do we need to normalize data?
            factor = factor_p_prof 
         endif

         write(iunit, '(a,a)', iostat=k) 
     1   '  r     z    phi   s    p-data    p-vmec   p-fit      sigma',
     2   'wgted dev.' 

         sigma_fit(1:np_prof) = sigma_p_prof(1:np_prof)

         do i=1, np_prof
            r_cyl(1) = r_p_prof(i)
            r_cyl(2) = phi_p_prof(i)*pi/180
            r_cyl(3) = z_p_prof(i)
            call ajax_cyl2flx(r_cyl, r_flx, iflag, message)

            s = r_flx(1)**2                     ! back to s
            s_prof(i) = s
            if( s > 1 ) sigma_fit(i) = bigno
         enddo

         i = 11
         call polyfit(s_prof, p_prof, sigma_fit, 
     1                np_prof, i, fit_am, chisq_fit)

         do i=1, np_prof
            num = num + 1

            index_array(num) = ivar
            chisq_descript(num) = descript(ivar)

            wegt(num) = sigma_fit(i) * factor
            chisq_target(num) = p_prof(i) * factor

            s = s_prof(i)
            sj = s * real(nrad-1, rprec) + 1.5   ! figure out where we are on half-grid
            j = floor(sj)
            fract = sj - j

            if( s > 1) then
               if( lp_prof_incl_edge) then
                  chisq_match(num) = 0 
               else
                  chisq_match(num) = chisq_target(num) 
               endif
            else if( j == nrad) then
               chisq_match(num) = pres_opt(nrad)
            else if( sj < 2) then
               chisq_match(num) = pres_opt(2)
            else
               chisq_match(num) = pres_opt(j) + 
     1                            fract*(pres_opt(j+1) - pres_opt(j))
            endif

            p_fit = eval_prof(fit_am, s)

            write(iunit, '(2f6.3,2f5.2,5es10.2)', iostat=k) 
     1            r_p_prof(i), z_p_prof(i), phi_p_prof(i), s, p_prof(i), 
     2            chisq_match(num), p_fit, sigma_p_prof(i),
     3            (chisq_match(num)-chisq_target(num))/wegt(num)

         end do

         write(iunit,*)
         write(iunit,*) ' Normalization factor = ',factor

         write(iunit,*)
         write(iunit,'(a)') 'Fitted am coefficients'
         do i=1, 11
            write(iunit,*) i-1, fit_am(i)
         enddo

         int_am = sum(am(0:10) / (/ 1,2,3,4,5,6,7,8,9,10,11 /) )
         int_fit_am = sum(fit_am(1:11) / 
     1                    (/ 1,2,3,4,5,6,7,8,9,10,11 /) )
         write(iunit,*)
         write(iunit,*) ' Normalization factor for fit= ',
     1                  int_am/int_fit_am

         close(iunit)
         if (k .ne. 0) then
            print *,'chisq_p_prof error writing file:',trim(extension)
            iflag = -12
            return
         end if  

      else
         num = num + np_prof
         lajax = .true.             ! we need to have AJAX initialized for the calculations
      end if   
      
      end subroutine chisq_p_prof


      subroutine polyfit(x, y, sigma, npts, nterms, a, chisq)
      use kind_spec
      implicit none
C-----------------------------------------------
c   fit (x(i), y(i)) points to polynomial in x
c   based on Eddington, p.140
c
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: npts, nterms
      real(rprec), intent(in) :: x(npts), y(npts), sigma(npts)
      real(rprec), intent(out) :: a(nterms)
      real(rprec) ::  chisq
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, j, k, l, n, nmax
      real(rprec) :: weight, xterm, yterm, xi, yi, delta
      real(rprec) :: sumx(2*nterms-1), sumy(nterms), 
     1               array(nterms, nterms) 
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      real(rprec) , external :: determ
C-----------------------------------------------
      nmax = 2*nterms - 1

      sumx = 0
      sumy = 0
      chisq = 0

c
c   accumulate weighted sums
c
      do i=1, npts
         xi = x(i)
         yi = y(i)
         weight = 1/sigma(i)**2
         xterm = weight
         yterm = weight * yi

         do n=1, nmax
            sumx(n) = sumx(n) + xterm
            xterm = xterm * xi
         enddo

         do n=1, nterms
            sumy(n) = sumy(n) + yterm
            yterm = yterm * xi
         enddo

         chisq = chisq + weight * yi**2

      enddo
c
c   construct matricies and calculate coefficients
c
      do j=1, nterms
         do k=1, nterms
            n = j + k - 1
            array(j,k) = sumx(n)
         enddo
      enddo

      delta = determ(array, nterms)

      if( delta == 0 ) then
c   singular!
         a = 0
         chisq = 0
         return

      else
         do l=1, nterms
            do j=1, nterms
               do k=1, nterms
                  n = j + k - 1
                  array(j,k) = sumx(n)
               enddo

               array(j,l) = sumy(j)
            enddo

            a(l) = determ(array, nterms) / delta
         enddo
      endif

c
c   calculate chisq
c
      do j=1, nterms
         chisq = chisq - 2 * a(j) * sumy(j)

         do k=1, nterms
            n = j + k - 1
            chisq = chisq + a(j) * a(k) * sumx(n)
         enddo
      enddo

      chisq = chisq / (npts - nterms)

      end subroutine polyfit


      function determ(a, n) result (d)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(n,n), intent(in) :: a
      integer, intent(in) :: n
      real(rprec) :: d
      
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, j, m, m2
      real(rprec) :: minor, minor2
C-----------------------------------------------

c      n = size(a,1)
      if( n > 3) then
         d = 0
         do i=1, n
            minor = a(1,i)
            minor2 = a(1,i)

            do j = 1, n-1
               m = mod(i+j-1, n) + 1
               m2 = mod(n+i-j-1, n) + 1     ! make sure it is positive

               minor = minor * a(j+1, m)
               minor2 = minor2 * a(j+1, m2)
            enddo

            d = d + minor - minor2
         enddo

      else if( n == 3) then
         d =   a(1,1) * (a(2,2)*a(3,3) - a(2,3)*a(3,2))
     1       + a(1,2) * (a(2,3)*a(3,1) - a(2,1)*a(3,3))  
     2       + a(1,3) * (a(2,1)*a(3,2) - a(2,2)*a(3,1))  

      else if( n == 2) then
         d = a(1,1) * a(2,2) - a(1,2) * a(2,1)

      else if( n == 1) then
         d = a(1,1)

      endif
      end function determ


      subroutine chisq_diagno(num, nopt, iflag, 
     1                            extension, lscreen)
      use kind_spec
      use chisq_mod
      use optim, ONLY: home_dir, chisq_min, bigno
      use optim_params, ONLY: sigma_diagno_seg,sigma_diagno_flx,
     1    target_diagno_seg, target_diagno_flx, ndiagno_seg, 
     2    ndiagno_flx, ldiagno_opt, diagno_control
      use system_mod, ONLY: system
      use safe_open_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: nopt, iflag, num
      character*(*) :: extension
      logical :: lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      character(len=len_trim(home_dir)+100) :: version
      character*48, dimension(:), allocatable :: seg_name, flx_name
      character*256 :: line
      integer :: i, k, istat, iunit
      real(rprec), dimension(:), allocatable :: segs, loops
      logical :: ex
C-----------------------------------------------
      version = trim(home_dir) // '/xdiagno'
      iflag = 0
      iunit = unit_diagno

      if (nopt > 0) then
          if (lscreen) write(6,*)'Running DIAGNO'

          call load_physics_codes(trim(version), " ", trim(extension),
     1         'diagno_seg', extension, iunit, iflag)

          if( iflag .ne. 0 ) return      ! failure

          allocate( segs(ndiagno_seg), loops(ndiagno_flx),
     1              seg_name(ndiagno_seg), flx_name(ndiagno_flx),
     2              stat=istat)
          if( istat .ne. 0) stop 'Allocation error in CHISQ_DIAGNO'

          read(iunit,*, iostat=istat) (segs(i), i=1, ndiagno_seg)
          if( istat .ne. 0)
     1              stop 'Read error in CHISQ_DIAGNO: segs data'

          do i=1, ndiagno_seg
             read(iunit,'(a)', iostat=istat) seg_name(i)
             if( istat .ne. 0) 
     1              stop 'Read error in CHISQ_DIAGNO: segs names'
          enddo

!          close(iunit, status='delete')
          close(iunit)

          call safe_open(iunit, istat, 
     1                   'diagno_flux.'//trim(extension), 
     1                   'old', 'formatted')
          if( istat .ne. 0 ) then
             deallocate( segs, loops)
             stop 'No flux-loop output file from Diagno!'
          endif

          read(iunit,*, iostat=istat) (loops(i), i=1, ndiagno_flx)
          if( istat .ne. 0) then
             deallocate( segs, loops)
             stop 'Read error in CHISQ_DIAGNO: loops'
          endif

          do i=1, ndiagno_flx
             read(iunit,'(a)', iostat=istat) flx_name(i)
             if( istat .ne. 0) 
     1              stop 'Read error in CHISQ_DIAGNO: flx names'
          enddo
!          close(iunit, status='delete')
          close(iunit)


          call safe_open(iunit, istat, 'mag_diags.'//trim(extension),
     1         'replace', 'formatted')
          if (istat .ne. 0) then
             iflag = -12
             deallocate( segs, loops)
             return
          endif

          write(iunit, '(a)')  
     1    "  calc.    measured  sigma    wgt'd dev.   name"

          do i=1, ndiagno_seg
             num = num + 1
             index_array(num) = ivar_diagno
             chisq_match(num) = segs(i)

             chisq_target(num) = target_diagno_seg(i)

             chisq_descript(num) = descript(ivar_diagno)
             wegt(num) = sigma_diagno_seg(i)

             if( sigma_diagno_seg(i) < bigno )
     1          write(iunit, '(4es10.2,a)') 
     1             segs(i),  chisq_target(num), wegt(num),
     2             (segs(i)-chisq_target(num))/wegt(num),
     3             trim(seg_name(i))
          enddo

          write(iunit, '(a)')  ' '

          do i=1, ndiagno_flx
             num = num + 1
             index_array(num) = ivar_diagno
             chisq_match(num) = loops(i)

             chisq_target(num) = target_diagno_flx(i)

             chisq_descript(num) = descript(ivar_diagno)
             wegt(num) = sigma_diagno_flx(i)
             
             if( sigma_diagno_flx(i) < bigno )
     1          write(iunit, '(4es10.2,a)') 
     1             loops(i),  chisq_target(num), wegt(num),
     2             (loops(i)-chisq_target(num))/wegt(num),
     3             trim(flx_name(i))
          enddo

          if( lscreen) then
             rewind(iunit)

             print *,' '
             print *,'Magnetic diagnostics matching'

             do while (.true.)
                read(iunit,'(a)', iostat=istat) line
                if( istat .ne. 0) exit
                print *,trim(line)
             enddo
          endif

          close(iunit)
          deallocate( segs, loops)

      else
          inquire(file=trim(version), exist=ex)
          if (.not.ex) then
            if(lscreen) then
              print *, 'xdiagno file not found in '// 
     1                trim(home_dir)
              print *, '*** Diagno evaluation disabled !'
            endif
            ldiagno_opt = .false.
            iflag = 1
          endif


          k = index(diagno_control, '/')
          if( k.eq.0) then               ! not full path
             inquire(file='../'//trim(diagno_control), exist=ex)
          else                           ! full path specified
             inquire(file=trim(diagno_control), exist=ex)
          endif

          if (.not.ex) then
            if(lscreen) then
              print *, "diagno control '",trim(diagno_control),
     1                 "' file not found" 
              print *, '*** Diagno evaluation disabled !'
            endif
            ldiagno_opt = .false.
            iflag = 1

          else if (k.eq. 0) then         ! not full path
            call system("/bin/ln -s ../" // trim(diagno_control) //
     1                  "  ./diagno.control" )
          else                           ! full path specified
            call system("/bin/ln -s " // trim(diagno_control) //
     1                  "  ./diagno.control" )
          endif

          num = num + ndiagno_seg + ndiagno_flx
      endif

      end subroutine chisq_diagno


      subroutine chisq_bripple (sigma, bmax_b,
     1          bmin_b, num, nopt, mskip, sqrt_nsurf, dtheta)
      use kind_spec
      use chisq_mod
      use boozer_params, ONLY: nu2_b
      use vparams, ONLY: twopi
      use optim, ONLY: bigno
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: nopt, mskip
      integer :: num
      real(rprec) :: bmax_b(*), bmin_b(*), sigma, sqrt_nsurf,dtheta
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: p5 = 0.5_dp, zero = 0.0_dp
      real(rprec), parameter :: max_variation = 2.e-2_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: n
      real(rprec) :: avg_bmin, cosweight
C-----------------------------------------------
      if (abs(sigma) .ge. bigno) return

      if (nopt .gt. 0) then

         avg_bmin = sum(bmin_b(:nu2_b))
     1            - p5*(bmin_b(1) + bmin_b(nu2_b))
         avg_bmin = avg_bmin * dtheta

         do n = 1, nu2_b, mskip
           num = num + 1
           index_array(num) = ivar_ripple
           wegt(num) = sigma*avg_bmin*sqrt_nsurf
           chisq_target(num) = zero
           chisq_descript(num) = descript(ivar_ripple)
           cosweight = sin (p5*twopi*real(n-1,rprec)/
     1                 real(2*(nu2_b-1),rprec))
           chisq_match(num) =
     1        cosweight*(bmax_b(n) - bmin_b(n))
           
         end do              

      else                   
         do n = 1, nu2_b, mskip
            num = num + 1
         end do   
      end if   

      end subroutine chisq_bripple


      subroutine chisq_3dbd(sigma, ivar, sigma_rms, sigma_max, 
     1           target, target_rms, lvv_tgt, ivar_rms, ivar_max,
     2           num, nopt, iflag, lscreen)
      use kind_spec
      use chisq_mod
      use optim_params, ONLY: ntor_vv, mpol_vv, rbc_vv, zbs_vv, 
     1    vv_dist, vv_dist_rms, vv_dist_max, shapeweight,
     2    phi_lim, r_lim, z_lim, nu_vv, nv_vv
      use optim, ONLY: bigno, nfp_opt, mnmax_opt, iunit_opt_local
      use boozer_params, ONLY: rmnc_bdy, zmns_bdy, xm_bdy, xn_bdy
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: ivar, ivar_rms, ivar_max, nopt
      integer :: num, iflag
      real(rprec) :: sigma, sigma_rms, sigma_max, target, target_rms
      logical, intent(in) :: lscreen, lvv_tgt
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i,m,n, ierr, mnmax_vv
      integer, parameter :: ngrid = 64
      real(rprec) :: dist, drms, lim_dist, adj_sigma
      real(rprec), dimension(mpol_vv*(2*ntor_vv+1)+ntor_vv+1) :: 
     1                    rmnc_vv, zmns_vv, xm_vv, xn_vv
      real(rprec), dimension(mnmax_opt) :: xn_per
      real(rprec), dimension(:,:), allocatable  :: grid_dist
      logical :: deviance
C-----------------------------------------------
      iflag = 0
      deviance = .false.
      vv_dist = 1.e30

      if (abs(sigma) .ge. bigno .and. abs(sigma_rms) .ge.  bigno
     1    .and. abs(sigma_max) .ge. bigno) return 

      if(shapeweight .and. 
     1   (.not. abs(sigma_rms) >=  bigno) ) deviance=.true.

      if (nopt .gt. 0) then
         mnmax_vv = mpol_vv*(2*ntor_vv+1)+ntor_vv+1

         allocate( grid_dist(ngrid, ngrid) )

         do n = 1, mnmax_opt
            xn_per(n) = nint(xn_bdy(n)/nfp_opt)       ! convert back to per period
         end do

!
!     first, calculated distance to a smooth boundary
!
         if( mnmax_vv > 1) then
c          write(6,*) '3dbd: mpol,ntor', mpol_vv, ntor_vv
           i = 0
           do m=0, mpol_vv
             do n=-ntor_vv, ntor_vv
               if (m.eq.0 .and. n.lt.0) cycle
               i = i+1
               xm_vv(i) = m
               xn_vv(i) = n
               rmnc_vv(i) = rbc_vv(n,m)
               zmns_vv(i) = zbs_vv(n,m)
             end do
           end do

           if( i .ne. mnmax_vv) then
             if (lscreen) print *, 
     1         'mnmax_vv mis-match in chisq_3dbd!:', mnmax_vv, i
             deallocate(grid_dist)
             stop
           endif

           ierr=0
           if(.not.deviance) then
             call surfsep(mnmax_opt, rmnc_bdy, zmns_bdy, xm_bdy,  
     1               xn_per, mnmax_vv, nu_vv, nv_vv, rmnc_vv, zmns_vv, 
     2               xm_vv, xn_vv, nfp_opt, ngrid, grid_dist, 
     3               vv_dist, vv_dist_rms, vv_dist_max, ierr)
           else
             call weighted_surfsep(rmnc_bdy, zmns_bdy, xm_bdy, xn_bdy, 
     1          nfp_opt, mnmax_opt, vv_dist_rms, vv_dist, vv_dist_max)
      
           endif
           if(ierr .ne. 0) then
             print *, 'SurfSep returned IERR = ',ierr, '!!!'
             iflag = -21

             deallocate(grid_dist)
             return
           endif

           if( lscreen) then
             write(6,'(a,3es12.3)')  'VV Distances: sep, rms, max:',  
     1                   vv_dist, vv_dist_rms, vv_dist_max
           endif
         endif

!
!    Calculate distance to a piecewise-linear limiter set
!
         if( r_lim(1,1) .ne. 0) then
            call limsep(mnmax_opt, rmnc_bdy, zmns_bdy, xm_bdy, xn_per, 
     1          size(phi_lim), phi_lim, size(r_lim,1), r_lim, z_lim, 
     2          nfp_opt, lim_dist, ierr)

            if(ierr .ne. 0) then
               print *, 'LimSep returned IERR = ',ierr, '!!!'
               iflag = -21

               deallocate(grid_dist)
               return
            endif

            if( lscreen) then
               write(6,'(a,3es12.3)')  
     1             'Min. Limiter Separation:', lim_dist
            endif

            vv_dist = min(vv_dist, lim_dist)
         endif

         if( sigma .lt. bigno ) then
            num = num + 1
            index_array(num) = ivar
            wegt(num) = sigma
            if( lvv_tgt) then
               chisq_target(num) = max(target, vv_dist)
            else
               chisq_target(num) = target
            endif
            chisq_descript(num) = descript(ivar)
            chisq_match(num) = vv_dist
         endif

         if( sigma_rms .lt. bigno .and. .not. deviance) then
            adj_sigma = sigma_rms * 
     1                  sqrt(float(count(grid_dist(:,:)<bigno)))

            do n=1, ngrid
               do m=1, ngrid
                  num = num + 1
                  index_array(num) = ivar_rms
                  chisq_target(num) = target_rms
                  chisq_descript(num) = descript(ivar_rms)
                  if( grid_dist(n,m) < bigno) then
                     wegt(num) = adj_sigma
                     chisq_match(num) = grid_dist(n,m)
                  else
                     wegt(num) = bigno
                     chisq_match(num) = target_rms
                  endif
               enddo
            enddo

         else if( sigma_rms .lt. bigno .and. deviance) then
            num = num + 1
            index_array(num) = ivar_rms
            chisq_target(num) = target_rms
            chisq_descript(num) = descript(ivar_rms)
            wegt(num) = sigma_rms
            chisq_match(num) = vv_dist_rms
         endif

         if( sigma_max .lt. bigno ) then
            num = num + 1
            index_array(num) = ivar_max
            wegt(num) = sigma_max
            chisq_target(num) = 0
            chisq_descript(num) = descript(ivar_max)
            chisq_match(num) = vv_dist_max
         endif

         deallocate(grid_dist)

      else
         if( sigma .lt. bigno) num = num + 1
         if( sigma_rms .lt. bigno) then
            if( .not. deviance) then
               num = num + ngrid*ngrid
            else
               num = num + 1
            endif
         endif
         if( sigma_max .lt. bigno) num = num + 1
      end if


      end subroutine chisq_3dbd


      subroutine chisq_3dbd_bd(sigma, ivar, sigma_rms, sigma_max, 
     1           target, target_rms, ivar_rms, ivar_max,
     2           num, nopt, iflag, lscreen)
      use kind_spec
      use chisq_mod
      use optim_params, ONLY: ntor_bd, mpol_bd, rbc_bd, zbs_bd, 
     1    bd_dist, bd_dist_rms, bd_dist_max, shapeweight,
     2    phi_lim, r_lim, z_lim, nu_bd, nv_bd
      use optim, ONLY: bigno, nfp_opt, mnmax_opt, iunit_opt_local
      use boozer_params, ONLY: rmnc_bdy, zmns_bdy, xm_bdy, xn_bdy
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: ivar, ivar_rms, ivar_max, nopt
      integer :: num, iflag
      real(rprec) :: sigma, sigma_rms, sigma_max, target, target_rms
      logical, intent(in) :: lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i,m,n, ierr, mnmax_bd
      integer, parameter :: ngrid = 64
      real(rprec) :: dist, drms, lim_dist, adj_sigma
      real(rprec), dimension(mpol_bd*(2*ntor_bd+1)+ntor_bd+1) :: 
     1                    rmnc_bd, zmns_bd, xm_bd, xn_bd
      real(rprec), dimension(mnmax_opt) :: xn_per
      real(rprec), dimension(:,:), allocatable :: grid_dist
      logical :: deviance
C-----------------------------------------------
      iflag = 0
      deviance = .false.
      bd_dist = 1.e30

      if (abs(sigma) .ge. bigno .and. abs(sigma_rms) .ge.  bigno
     1    .and. abs(sigma_max) .ge. bigno) return 

      if (nopt .gt. 0) then
         mnmax_bd = mpol_bd*(2*ntor_bd+1)+ntor_bd+1

         allocate( grid_dist(ngrid, ngrid) )

         do n = 1, mnmax_opt
            xn_per(n) = nint(xn_bdy(n)/nfp_opt)       ! convert back to per period
         end do

!
!     first, calculated distance to a smooth boundary
!
         if( mnmax_bd > 1) then
c          write(6,*) '3dbd_bd: mpol,ntor', mpol_bd, ntor_bd
           i = 0
           do m=0, mpol_bd
             do n=-ntor_bd, ntor_bd
               if (m.eq.0 .and. n.lt.0) cycle
               i = i+1
               xm_bd(i) = m
               xn_bd(i) = n
               rmnc_bd(i) = rbc_bd(n,m)
               zmns_bd(i) = zbs_bd(n,m)
             end do
           end do

           if( i .ne. mnmax_bd) then
             if (lscreen) print *, 
     1         'mnmax_bd mis-match in chisq_3dbd_bd!:', mnmax_bd, i
             deallocate(grid_dist)
             stop
           endif

           ierr=0
           call surfsep(mnmax_opt, rmnc_bdy, zmns_bdy, xm_bdy,  
     1               xn_per, mnmax_bd, nu_bd, nv_bd, rmnc_bd, zmns_bd, 
     2               xm_bd, xn_bd, nfp_opt, ngrid, grid_dist, 
     3               bd_dist, bd_dist_rms, bd_dist_max, ierr)
      
           if(ierr .ne. 0) then
             print *, 'SurfSep returned IERR = ',ierr, '!!! in PL'
             iflag = -21

             deallocate(grid_dist)
             return
           endif

           if( lscreen) then
             write(6,'(a,3es12.3)')  'BD Distances: sep, rms, max:',  
     1                   bd_dist, bd_dist_rms, bd_dist_max
           endif
         endif

         if( sigma .lt. bigno ) then
            num = num + 1
            index_array(num) = ivar
            wegt(num) = sigma
            chisq_target(num) = target
            chisq_descript(num) = descript(ivar)
            chisq_match(num) = bd_dist
         endif

         if( sigma_rms .lt. bigno ) then
            adj_sigma = sigma_rms * 
     1                  sqrt(float(count(grid_dist(:,:)<bigno)))

            do n=1, ngrid
               do m=1, ngrid
                  num = num + 1
                  index_array(num) = ivar_rms
                  chisq_target(num) = target_rms
                  chisq_descript(num) = descript(ivar_rms)
                  if( grid_dist(n,m) < bigno) then
                     wegt(num) = adj_sigma
                     chisq_match(num) = grid_dist(n,m)
                  else
                     wegt(num) = bigno
                     chisq_match(num) = bd_dist_rms
                  endif
               enddo
            enddo
         endif

         if( sigma_max .lt. bigno ) then
            num = num + 1
            index_array(num) = ivar_max
            wegt(num) = sigma_max
            chisq_target(num) = 0
            chisq_descript(num) = descript(ivar_max)
            chisq_match(num) = bd_dist_max
         endif

         deallocate(grid_dist)
      else
         if( sigma .lt. bigno) num = num + 1
         if( sigma_rms .lt. bigno) num = num + ngrid*ngrid
         if( sigma_max .lt. bigno) num = num + 1
      end if

      end subroutine chisq_3dbd_bd


      subroutine chisq_extcur(sigma, target, ivar, num, nopt)
      use kind_spec
      use chisq_mod
      use optim, ONLY: bigno, lextcur, nextcur_vmec
      use vmec_input, ONLY: lfreeb, extcur
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: ivar, nopt
      integer :: num
      real(rprec) :: sigma(*), target(*)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, c1p5 = 1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: j, n
C-----------------------------------------------
      n = size(lextcur)

      if (all(abs(sigma(:n)) >= bigno .or. .not.lextcur(:n)) 
     1       .or. .not. lfreeb) return

      if (nopt .gt. 0) then
         do j = 1, n
          if( abs(sigma(j)) < bigno .and. lextcur(j)) then
            num = num + 1
            index_array(num) = ivar

!               Must count it this way or else num value will be (possibly) thrown off

            if (j <= nextcur_vmec) then
               wegt(num) = sigma(j)
               chisq_match(num) = extcur(j)
               chisq_target(num) = target(j)
            else
               chisq_match(num) = zero
               chisq_target(num) = zero
            end if   
            chisq_descript(num) = descript(ivar)
          endif
         end do
           
      else
         !This MAY be an overestimate, but we may NOT know nextcur_vmec at this first call
         num = num + count(lextcur(:n).and.(abs(sigma(:n))<bigno))
      endif

      end subroutine chisq_extcur
            

      subroutine chisq_oh(sigma, ivar, num, nopt)
      use kind_spec
      use chisq_mod
      use optim, ONLY: bigno, lextcur, oh_coefs, nextcur_vmec
      use vmec_input, ONLY: lfreeb, extcur
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: ivar, nopt
      integer :: num
      real(rprec) :: sigma
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, c1p5 = 1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: j, n
      real(rprec) :: sum0
C-----------------------------------------------

      if (abs(sigma) >= bigno .or. .not. lfreeb) return

      if (nopt .gt. 0) then
         
         sum0 = sum(oh_coefs(:nextcur_vmec)*extcur(:nextcur_vmec))

         num = num + 1
         index_array(num) = ivar

         wegt(num) = sigma
         chisq_target(num) = zero
         chisq_match(num) = sum0
         chisq_descript(num) = descript(ivar)

      else
         num = num + 1
      endif

      end subroutine chisq_oh
            

      subroutine chisq_kappa(target, sigma, ivar, num, nopt)
      use kind_spec
      use chisq_mod
      use optim, ONLY: bigno, lextcur, mpol1_opt, mnmax_opt
      use boozer_params, ONLY: rmnc_bdy, zmns_bdy, xm_bdy, xn_bdy
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: ivar, nopt
      integer :: num
      real(kind=rprec) :: sigma, target
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(kind=rprec), parameter :: zero = 0, c1p5 = 1.5_dp,
     1    xtpi=3.141592654
      integer, parameter :: nthsrch = 50
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, j, jmax, k
      real(rprec) rmin, rmax, zmax, znew, thsrch, th0, dth
      real(rprec), dimension(0:mpol1_opt) :: rbc, zbs
C-----------------------------------------------
      if ( abs(sigma) >= bigno) return

      if (nopt .gt. 0) then
         num = num + 1
         index_array(num) = ivar

         wegt(num) = sigma
         chisq_target(num) = target
         chisq_descript(num) = descript(ivar)

!    find the n=0 terms
         rbc = 0
         zbs = 0

         do j=1, mnmax_opt
            if( nint(xn_bdy(j)) == 0 ) then
               i = nint(xm_bdy(j))
               rbc(i) = rmnc_bdy(j)
               zbs(i) = zmns_bdy(j)
            endif
         enddo

         rmax = sum(rbc(0:mpol1_opt))
         rmin = sum(rbc(0:mpol1_opt:2)) - sum(rbc(1:mpol1_opt:2))
         zmax = 0

         th0 = 0
         dth = xtpi/real(nthsrch+1,kind=rprec)
         jmax = 0

         do k=1, 2
           do j=1,nthsrch
             thsrch=th0 + j*dth

             znew = sum(zbs(0:mpol1_opt)*
     1             cos(thsrch*(/ (i, i=0, mpol1_opt) /) ))

             if( znew > zmax ) then
                jmax = j - 1
                zmax = znew
             endif
           enddo

           th0 = th0 + jmax*dth
           dth = (2*dth)/(nthsrch+1)
         enddo

         chisq_match(num) = 2*zmax/(rmax - rmin)

      else
         num = num + 1
      endif

      end subroutine chisq_kappa


      subroutine chisq_ellipticity(target, sigma, num, nopt)
      use kind_spec
      use chisq_mod
      use optim, ONLY: bigno, mnmax => mnmax_opt
      use vparams, ONLY: zero, one, twopi
      use boozer_params, ONLY: rmnc_bdy, zmns_bdy, xm_bdy, xn_bdy
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: nopt
      integer :: num
      real(rprec), intent(in) :: target, sigma
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: nlen = 100
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, m, mn
      real(rprec), dimension(nlen) :: theta, z1
      real(rprec) :: ellipticity, rmax0, rmin0, zmax0, zmin0
C-----------------------------------------------

      if (sigma .ge. bigno) return

      num = num + 1
      
      if (nopt .gt. 0) then

         do i = 1,nlen
            theta(i) = (i-1)*(twopi/nlen)
         enddo
!
!        Limits waist thickness in symmetry plane
! 
         rmax0 = sum(rmnc_bdy(:mnmax))
         rmin0 = sum(rmnc_bdy(:mnmax) * (-one)**nint(xm_bdy(:)))
         z1 = zero
         do mn = 1, mnmax
            m = nint(xm_bdy(mn))
            z1 = z1 + zmns_bdy(mn)*sin(m*theta)
         end do

         zmax0 = maxval(z1(:))
         zmin0 = minval(z1(:))
         ellipticity = abs(zmax0 - zmin0)/abs(rmax0 - rmin0)

         index_array(num) = ivar_ellipticity
         wegt(num) = sigma
         if (wegt(num) .eq. zero) wegt(num) = one
         chisq_descript(num) = descript(ivar_ellipticity)
         chisq_target(num) = target
         chisq_match(num) = ellipticity
      end if

      end subroutine chisq_ellipticity


      subroutine chisq_coilgeom (num, nopt, iflag, extension, lscreen)
      use kind_spec
      use real_ptr_type
      use chisq_mod
      use coilsnamin
      use safe_open_mod
      use optim, ONLY: home_dir, sigma_berr_avg, sigma_berr_max, bigno
      use system_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: nopt
      integer :: num, iflag
      character*(*) :: extension
      logical :: lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: n_list=9
      integer, parameter :: coilgeom0 = 28
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: j, k, iunit, icount, istat, ierr
      type (REAL_PTR), dimension(n_list) :: weights, weightm
      real(rprec), dimension(1), target :: scalar1
      real(rprec), dimension(2) :: no_weights
      logical :: access, vf_reg, polcur, mplcoil, splcoil, vf_rmax,
     1           berrflag
      logical, dimension(n_list) :: mflag, sflag
      integer, dimension(n_list) :: ivars
      character*30 :: header
      character*200 :: output_file, cmd

!****************************************************************************
!
!   **** NOTE CAREFULLY ****
!
!  This routine must be maintained in concert with COILGEOM/lsfun1
!
!****************************************************************************
!
!   Setup the flags indicating what is to be done
!   initialize

      mflag = .false.
      sflag = .false.
      berrflag = .false.
      access = .false.
      vf_reg = .false.
      vf_rmax = .false.
      polcur = .false.
      mplcoil = .false.
      splcoil = .false.

      no_weights = 1

!   setup ivars for the initial 5 types of penalties

      ivars(1) = ivar_coil_len
      ivars(2) = ivar_coil_ymin
      ivars(3) = ivar_coil_sep
      ivars(4) = ivar_coil_plcross
      ivars(5) = ivar_coil_curv
      ivars(6) = ivar_coil_icurv
      ivars(7) = ivar_coil_reg
      ivars(8) = ivar_coil_curdens
      ivars(9) = ivar_coil_rmax
      

!   here we set up flags indicating whether there are any values we want
!   to include for the initial 5 types of penalties, for both modular and
!   saddle representations.   Note that we must carefully parallel the logic
!   used in COILGEOM in writing these out.

!   'modular' representation

      if (lmodular .and. nmid >= 1) then
         mflag(1) = (.not. lsaddle) .and. any(lmod_wgt(1:nmid).ne.zero)
         mflag(2) = (.not. lsaddle) .and. any(ymin_wgt(1:nmid).ne.zero)
         mflag(3) = any( dcc_wgt(1:nmid) .ne. zero)
         mflag(4) = .false.                                      ! someday
         mflag(5) = any( rc_wgt(1:nmid) .ne. zero)
         mflag(6) = any( cu_wgt(1:nmid) .ne. zero)
         mflag(7) = .false.
         mflag(8) = .false.
         mflag(9) = .false.

         weightm(1)%x => lmod_wgt
         weightm(2)%x => ymin_wgt
         weightm(3)%x => dcc_wgt
         weightm(5)%x => rc_wgt
         weightm(6)%x => cu_wgt

         mplcoil = dcp_wgt .ne. zero
      endif

!   'saddle' representation

      if (lsaddle .and. nsmid >= 1) then
         sflag(1) = any( lsad_wgt(1:nsmid) .ne. zero)
         sflag(2) = any( ymin_wgt(1:nsmid) .ne. zero)
         sflag(3) = any( dsc_wgt(1:nsmid) .ne. zero)
         sflag(4) = any( dscxp_wgt(1:nsmid) .ne. zero)
         sflag(5) = any( rs_wgt(1:nsmid) .ne. zero)
         sflag(6) = any( cs_wgt(1:nsmid) .ne. zero)
         sflag(7) = any( csc_wgt(1:nsmid) .ne. zero)
         sflag(8) = any( scd_wgt(1:nsmid) .ne. zero)
         sflag(9) = any( rmax_wgt(1:nsmid) .ne. zero)

         weights(1)%x => lsad_wgt
         weights(2)%x => ymin_wgt
         weights(3)%x => dsc_wgt
         weights(4)%x => dscxp_wgt
         weights(5)%x => rs_wgt
         weights(6)%x => cs_wgt
         weights(7)%x => csc_wgt
         weights(8)%x => scd_wgt
         weights(9)%x => rmax_wgt

         splcoil = dcp_wgt .ne. zero
      endif

!   set flags for the penalties that are independent of representation

!   Berr
      berrflag = lbnorm .and. (sigma_berr_avg < bigno .or.
     1                       sigma_berr_max < bigno)

!   access

      if (laccess .and. n_access > 0)
     1   access = any (dac_wgt(1:n_access) .ne. zero)

!   PF specific: regularization and rmax

      if (lvf .and. num_vf > 0 ) then
         vf_reg = any(cvf_wgt(1:num_vf) .ne. zero)
         vf_rmax = any(rvf_wgt(1:num_vf) .ne. zero) .and. lvfvar
      endif 

!   total poloidal current

      if( lpolcur) polcur = dpc_wgt .ne. zero

! ------------------------------------------------

      if (nopt >= 0) then

!   if necessary, calculate Bnorm

         if( lbnorm) then
            cmd = trim(home_dir) // '/xbnorm wout.' // trim(extension)
            write (cmd,'(a,1pe12.3)') trim(cmd), 0.20
            if (lscreen) then
               print '(/,a)',' Running BNORM code...'
            end if
            call system(cmd)                     !!Produce bnorm file
         endif


!   execute xcoilgeom again (called previously from generate_mgrid) to generate plasmas-dependent penalties

         output_file = 'coil_targets.' // trim(extension)
         cmd = trim(home_dir) // '/xcoilgeom ' // trim(extension)
     1                    // ' -P > ' // trim(output_file)
         call system(cmd, ierr)
         if (ierr.lt.127 .and. ierr.ne.0) then
            if (lscreen) print *,
     1                   'XCOILGEOM failed in chisq_coilgeom call: ',
     2                           'ierr = ',ierr
            iflag = -22
            return
         end if

!
!  Open and read coil_targets file generated when coils file written
!
         iunit = coilgeom0

         call safe_open(iunit, k, output_file, 'old', 'formatted')
         if (k .ne. 0) then
            iflag = -25
            return
         endif

         header = ' '
         do while (index(header,'coilgeom penalties') == 0)
            read(iunit,'(a)',iostat=istat) header
            if (istat .ne. 0) stop ' ISTAT != 0 in chisq_coilgeom(1)'
         enddo

!   the first penalties are for Berr avg. & max

         call read_coilgeom_elem(iunit, num, 2, no_weights,
     1                              berrflag, ivar_coil_aBerr)
         if( berrflag) then
            index_array(num) = ivar_coil_mBerr
            chisq_descript(num) = descript(ivar_coil_mBerr)
            wegt(num) = sigma_berr_max
            wegt(num-1) = sigma_berr_avg
         endif

!   the next  n_list  penalties occur in pairs, first for the modular representation,
!   then for the saddle representation.  This pairing is broken due to some
!   penalties not being implemented for one representation or the other in
!   COILGEOM/COILOPT, causing the extra tests

         do j=1, n_list
            if (lmodular .and. (j .ne. 2 .or. (.not. lsaddle))
     1                   .and. (j.ne.4) .and. (j.ne.7) 
     2                   .and. (j.ne.8) .and. (j.ne.9) ) then
               call read_coilgeom_elem(iunit, num, nmid, weightm(j)%x,
     1                                 mflag(j), ivars(j))
            end if

            if (lsaddle) then
               call read_coilgeom_elem(iunit, num, nsmid, weights(j)%x,
     1                                 sflag(j), ivars(j))
            end if
         end do

!   Coil-Plasma distances
         scalar1(1) = dcp_wgt;  weights(1)%x => scalar1
         if (lmodular) call read_coilgeom_elem(iunit, num, 1,
     1                    weights(1)%x, mplcoil, ivar_coil_pldist)

         if (lsaddle) call read_coilgeom_elem(iunit, num, 1,
     1                          weights(1)%x, splcoil, ivar_coil_pldist)


!   Now we have penalties that do not depend on representations used

         if (laccess) then
            weights(1)%x => dac_wgt
            call read_coilgeom_elem(iunit, num, n_access, weights(1)%x,
     1                              access, ivar_coil_acc)
         endif

         if (lvf) then
            weights(1)%x => cvf_wgt
            call read_coilgeom_elem(iunit, num, num_vf, weights(1)%x,
     1                              vf_reg, ivar_coil_vfreg)
         endif

         if (lvfvar) then
            weights(1)%x => rvf_wgt
            call read_coilgeom_elem(iunit, num, num_vf, weights(1)%x,
     1                              vf_rmax, ivar_coil_vfrmax)
         endif

         if (lpolcur) then
            scalar1(1) = dpc_wgt;    weights(1)%x => scalar1
            call read_coilgeom_elem(iunit, num, 1, weights(1)%x, polcur,
     1                              ivar_coil_polcur)
         endif

         close (iunit)

      else
         num = num + nmid * count(mflag(:))
     1             + nsmid * count(sflag(:))

         if (mplcoil) num = num + 1
         if (splcoil) num = num + 1
         if (berrflag) num = num + 2
         if (access) num = num + n_access
         if (vf_reg) num = num + num_vf
         if (vf_rmax) num = num + num_vf
         if (polcur) num = num + 1

         iflag = 0
      end if

      end subroutine chisq_coilgeom


      subroutine read_coilgeom_elem(iunit, num, nwant, weight, flag, 
     1           ivars)
      use kind_spec
      use chisq_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: iunit, nwant, ivars      
      integer :: num
      real(rprec), dimension(nwant) :: weight
      logical, intent(in) :: flag
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------

C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: icount, istat, nums
      real(rprec), parameter :: one = 1, epsm = epsilon(one)
      real(rprec), dimension (nwant) :: values
C-----------------------------------------------
!
!     NOTE: THE VALUES ALREADY HAVE THE WEIGHTS ON THEM 
!     TO BE CONSISTENT WITH REST OF STELLOPT, WE DIVIDE THEM OUT HERE
!     AND MULTIPLY THEM BACK LATER
!
      read (iunit,*, iostat=istat) icount
      if (istat .ne. 0) stop 'ISTAT != 0 in read_coilgeom_elem(1)'
      if (icount .ne. nwant) 
     1          stop 'Count mismatch in read_coilgeom_elem'

      read (iunit,*, iostat=istat) values(1:nwant)
      if (istat .ne. 0) stop 'ISTAT != 0 in read_coilgeom_elem(2)'

      if (flag) then
         nums = num+nwant
         wegt(num+1:nums) = one/(weight(:nwant) + epsm)
         chisq_match(num+1:nums) = values(:nwant) * wegt(num+1:nums)
         chisq_target(num+1:nums) = 0
         index_array(num+1:nums) = ivars
         chisq_descript(num+1:nums) = descript(ivars)
         num = num + nwant
      endif

      end subroutine read_coilgeom_elem


      subroutine open_comm_files (mboz, nboz, ns_booz, ns_booz_max,
     1  ns_surf, ns_surf_max, lbootsj_opt, extension, iflag)
      use kind_spec
      use safe_open_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: mboz, nboz, ns_booz(*), ns_booz_max, 
     1   ns_surf(*), ns_surf_max
      integer :: iflag
      logical :: lbootsj_opt      
      character*(*), intent(in)  :: extension
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: unit_booz=24, unit_boot=26
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: iunit, k

      iunit = unit_booz
      call safe_open(iunit, k, 'in_booz.'//trim(extension), 'replace', 
     1    'formatted')
      if (k .ne. 0) then
         iflag = -14
         return
      end if   
      write(iunit, *) mboz, nboz           !!User may have input these through &OPTIMUM namelist
      write(iunit, *) trim(extension)
      write(iunit, *) ns_booz(1:ns_booz_max)
      close(iunit)

      if (lbootsj_opt) then
         iunit = unit_boot
         call safe_open(iunit, k, 'in_bootsj.'//trim(extension),
     1             'replace', 'formatted')
         if (k .ne. 0) then
            iflag = -14
            return
         end if   
         write(iunit, *) trim(extension)
         write(iunit, *) ns_surf(1:ns_surf_max)
         close(iunit)
      end if  

      end subroutine open_comm_files
      
      
      subroutine read_wout_opt(lread, extension, ierr, iopen)
      use read_wout_mod, version_wout => version_
      use vmec_input, only: nfp_in => nfp, mpol_in => mpol, 
     1   ntor_in => ntor, raxis_in => raxis, zaxis_in => zaxis,
     2   lfreeb, rbc, zbs
      use vparams, only: mpold, ntord, zero, one
      use optim, ierr_vmec_opt => ierr_vmec
      use boozer_params
      use mpi_params                                                     !MPI
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: ierr, iopen
      character*(*) :: extension
      logical lread
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: p03=3.e-2_dp, p5=0.5_dp, two=2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: mn, js, n1, mb, nb
      real(rprec) :: dum
C-----------------------------------------------
 
!     COMPUTE ACTUAL NO. THETA, PHI POINTS FOR INTEGRATIONS
!     NEEDED FOR DYNAMIC MEMORY ALLOCATION
 
      if (.not.lread) then
         mboz = max(6*mpol_in,2,mboz_opt)        !USER MAY SPECIFY mboz_opt
         nboz = max(2*ntor_in - 1,1,nboz_opt)    !USER MAY SPECIFY nboz_opt
         nu_boz   = 6*mboz+2
         nv_boz   = 4*nboz+2
         nu_boz   = nu_boz + mod(nu_boz,2)       !nu_boz, nv_boz MUST be even
         nv_boz   = nv_boz + mod(nv_boz,2)
         mnboz = nboz + 1 + (mboz - 1)*(1 + 2*nboz)
         nu2_b = nu_boz/2 + 1
         nunv = nu2_b*nv_boz
!        nunv = nu_boz*nv_boz
         ierr = 0
         iopen = 0
         return
      end if   
 
!
!     READ_WOUT_FILE: IN READ_WOUT_MODULE (MISCEL.F90 FILE)
!
      call read_wout_file('wout.' // trim(extension), ierr, iopen)
      
      if (iopen .ne. 0) return
      
      version_opt = version_wout                                   !!COBRA
      wb_opt   = wb
      wp_opt   = wp
      rmax_opt = rmax_surf
      rmin_opt = rmin_surf
      zmax_opt = zmax_surf
      aminor_opt = aminor
      nfp_opt  = nfp
      mpol_opt = mpol
      ntor_opt = ntor
      mnmax_opt= mnmax
      ierr_vmec_opt = ierr_vmec

      if (ierr_vmec.ne.0 .or. ierr.ne.0) then
         print *,' For processor ', myid+1,
     1   ' ierr_vmec = ', ierr_vmec,' in read_wout_opt'
         print *,' optimization will continue,',
     1   ' but this search direction will be deprecated!'
         if (ierr_vmec_opt .eq. 0) ierr_vmec_opt = ierr
         goto 1000
      end if   
!
!     Confirm that WOUT file contains same data read in from input file
!
      if (nfp_opt.ne.nfp_in .or. ns.ne.nrad .or. mpol_opt.ne.mpol_in
     1    .or. ntor_opt.ne.ntor_in) then
         if (myid .eq. master) then                                      !START MPI
            print *,' NFP_IN  = ',nfp_in, ' NFP_WOUT  = ',nfp_opt
            print *,' NS_IN   = ',nrad,   ' NS_WOUT   = ',ns
            print *,' MPOL_IN = ',mpol_in,' MPOL_WOUT = ',mpol_opt
            print *,' NTOR_IN = ',ntor_in,' NTOR_WOUT = ',ntor_opt
         endif                                                           !END MPI
         ierr = -1
         return
      end if   
      if (mpol_opt.gt.mpold .or. ntor_opt.gt.ntord)
     1  stop 'mpold or ntord too small in read_wout'
        mpol1_opt = mpol_opt - 1
        ntor1_opt = ntor_opt + 1
 
      do mn = 1, mnmax
         xm_bdy(mn)   = xm(mn)
         xn_bdy(mn)   = xn(mn)
         if (nint(xm(mn)) .eq. 0) then
            n1 = abs(nint(xn(mn)))/nfp
            raxis_in(n1,1) = rmnc(mn,1) 
            zaxis_in(n1,1) = zmns(mn,1)
         end if   
         rmnc_bdy(mn) = rmnc(mn,nrad)
         zmns_bdy(mn) = zmns(mn,nrad)
      end do

      rmnc_opt(:mnmax,:nrad) = rmnc(:mnmax,:nrad)
      zmns_opt(:mnmax,:nrad) = zmns(:mnmax,:nrad)
      lmns_opt(:mnmax,:nrad) = lmns(:mnmax,:nrad)
      
!     Account for free boundary possibly changing boundary coefficients        
!     need to write out correctly in write_rbzb routine
      if (lfreeb) then
         rbc = 0; zbs = 0
         do mn = 1, mnmax
            nb = nint(xn_bdy(mn)/nfp)
            mb = nint(xm_bdy(mn))
            rbc(nb,mb) = rmnc_bdy(mn)
            zbs(nb,mb) = zmns_bdy(mn)
         end do
      end if
      
      pres_opt(2:nrad) = pres(2:nrad)                              !!COBRA
      iota_opt(2:nrad) = iotas(2:nrad)
      phip_opt(2:nrad) = phip(2:nrad)
      buco_opt(2:nrad) = buco(2:nrad)
      vp_opt(2:nrad)   = vp(2:nrad)
      jcurv_opt(1:nrad) = jcurv(1:nrad)                            !!local <current density>, NOT integrated in s
      jdotb_opt(1:nrad) = jdotb(1:nrad)
      aspect_opt = aspect
      rbtor_opt  = rbtor

!
!     MERCIER CRITERION (NORM IT FOR USE IN OPTIMIZER)
!
      Dmerc_opt(2:nrad-1) = Dmerc(2:nrad-1)
      Dmerc_opt(1) = 0
      Dmerc_opt(nrad) = 0
      
      do js = 2,nrad-1
         dum = Dshear(js) + abs(Dwell(js)) + abs(Dcurr(js))
         dum = 0.1_dp*max(dum, abs(Dgeod(js)))
         if (dum .gt. zero) Dmerc_opt(js) = Dmerc_opt(js)/dum
      end do
      
!
!     IF POOR FORCE BALANCE, DO NOT BELIEVE MERCIER...
!
      equif(2:nrad-1) = abs(equif(2:nrad-1))/p03
      where (equif(2:nrad-1) .gt. one) 
     1    Dmerc_opt(2:nrad-1) = Dmerc_opt(2:nrad-1) / equif(2:nrad-1)

!
!     PUT BUCO_OPT ON FULL MESH
!
      do js = 2,nrad-1
        buco_opt(js) = p5*(buco_opt(js) + buco_opt(js+1))
      end do
      
      buco_opt(1) = zero
      buco_opt(nrad) = two*buco_opt(nrad) - buco_opt(nrad-1)
 
!
!     Deallocate memory
!
 1000 call read_wout_deallocate

      end subroutine read_wout_opt
 

      subroutine write_indata(input_file, istat)
      use optim, input_file_opt=>input_file
      use vmec_input
      use vparams, ONLY: zero
      use bootsj_input
      use safe_open_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(out) :: istat
      character*(*) :: input_file
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, dimension(1) :: ins
      integer :: iftol, m, n, i, j, k, ns0, iunit,  
     1    nextcur_max
      character*100, dimension(9), save :: comments
      character*100 :: temp_write
      character*(*), parameter :: form1 = "(a,i1,a)",
     1    form2 = "(a,i2,a)"
      character*8 :: form
      data comments/
     1  'OPTIMIZER RUN CONTROL PARAMETERS (INCLUDING MULTI-PROCESSING)',
     2  'LOGICAL VARIABLES CONTROLLING ACCESS TO PHYSICS MODULES',
     3  'SIZE AND FIELD STRENGTH SCALE PARAMETERS',
     4  'BOOZER COORDINATE TRANSFORMATION PARAMETERS',
     5  'SCALAR OPTIMIZATION PARAMETERS',
     6  'TRANSPORT OPTIMIZATION PARAMETERS',
     7  'BALLOONING MODE OPTIMIZATION PARAMETERS',
     8  'BOOTSTRAP CURRENT OPTIMIZATION PARAMETERS',
     9  'DKES AND NEO OPTIMIZATION PARAMETERS'
     A   /      
C-----------------------------------------------
      iunit = iout+2
      call safe_open(iunit, istat, input_file, 'replace', 'formatted')
      if (istat .ne. 0) then
         istat = -3
         return
      end if
      
      ins = maxloc(ns_array)
      ns0 = ns_array(ins(1))

      iftol = 1
      do while(ftol_array(iftol).ne.zero .and. iftol.lt.100)
         iftol = iftol + 1
      end do
 
!
!      COMPUTE nextcur_vmec, which could have changed if lcoil_geom = .true.
!
      do nextcur_vmec = size(extcur), 1, -1
         if (extcur(nextcur_vmec) .ne. zero) exit
      end do
      nextcur_vmec = max(0, nextcur_vmec)

      do nextcur_max = size(lextcur), 1, -1
         if (lextcur(nextcur_max)) exit
      end do
      nextcur_max = max(0, nextcur_max)
!
!     Start writing out formatted VMEC data with &INDATA heading
!     so that it can be read in under the indata namelist
!
      write (iunit, '(a7)') '&INDATA'
      write (iunit, '(2x,3a)') "MGRID_FILE = '", trim(mgrid_file), "'"
      write (iunit, 1007) 'LFREEB = ', lfreeb
      if( loldout) write (iunit, 1007) 'LOLDOUT = ', loldout 
      if( ldiagno) write (iunit, 1007) 'LDIAGNO = ', ldiagno 
      write (iunit, '(2x,a,1pe10.2)') 'DELT = ', delt
      write (iunit, '(2x,a,1pe10.2)') 'TCON0 = ', tcon0
      write (iunit, '(2x,a,i4)') 'NFP = ', nfp
      write (iunit, '(2x,a,i4)') 'NCURR = ', ncurr
      write (iunit, '(2(2x,a,i4))') 'MPOL = ', mpol, 'NTOR = ', ntor
      if (ntheta.gt.0) write (iunit, '(2x,a,i4)') 'NTHETA = ',ntheta 
      if (nzeta.gt.0)  write (iunit, '(2x,a,i4)') 'NZETA  = ',nzeta 
      write (iunit, '(2x,a)') 'NS_ARRAY = '
      write (iunit, 980, err=2000) (ns_array(i),i=1,ins(1))
      write (iunit, '(2x,a,i4)') 'NITER = ', niter
      write (iunit, '(2x,a,i4)') 'NSTEP = ', nstep
      write (iunit, '(2x,a,i4)') 'NVACSKIP = ', nvacskip
      write (iunit, 994) gamma
      write (iunit, '(2x,a)') 'FTOL_ARRAY ='
      write (iunit, 995, err=2000) (ftol_array(i),i=1,iftol - 1)
      write (iunit, 996) phiedge
      write (iunit, '(2x,a,1pe14.6)') 'BLOAT = ', bloat
      if (nextcur_vmec .gt. 0)
     1    call write_array(iunit, 'EXTCUR', extcur, nextcur_vmec)

      call write_rbzb (iunit, istat)
      if (istat .ne. 0) goto 2000
      
  980 format(2x,16i5)
  981 format(2x,a,10(i5))
  990 format(2x,a,6(1pe22.14))
  991 format(x,a,5(1pe122.14))
  994 format(2x,'GAMMA = ',1pe14.6)
  995 format(2x,1p4e14.6)
  996 format(2x,'PHIEDGE = ',1pe22.14)
 1000 format(5(1pe22.14))
 1005 format('  LFIX_NTOR(',i1,') = T')
 1006 format('  LFIX_NTOR(',i2,') = T')
 1007 format(2x,a,L1)

!
!     Finish by writing out formatted optimization data with &OPTIMUM
!     so that it can be read in under the optimum namelist
!     Option to add comments (starting with '!') now available
!
      write (iunit, '(a)') '/'
      write (iunit, '(a8)') '&OPTIMUM'
      write (iunit, 100) comments(1)
      write (iunit, 990) 'EPSFCN = ', epsfcn
      write (iunit, '(2x,a,i5)') 'NITER_OPT = ', niter_opt
      write (iunit, '(2x,a,i5)') 'NUM_PROCESSORS = ', num_processors
      write (iunit, '(2x,a,i5)') 'NUM_LEVMAR_PARAMS = ', 
     1   num_levmar_params
      write (iunit, '(2x,a)') 'NSURF_MASK = '
      write (iunit, '(10f7.3)', err=2000) (nsurf_mask(i),i=1,ns0)
      
      write (iunit, '(2x,a,i5)') 'NOPT_ALG = ', nopt_alg
      write (iunit, '(2x,a,i5)') 'NOPT_BOUNDARY = ', nopt_boundary
      write (iunit, 1007) 'LRESET_OPT = ', lreset_opt
      write (iunit, 1007) 'LDIAG_OPT = ', ldiag_opt
      write (iunit, 1007) 'LKEEP_MINS = ', lkeep_mins
      write (iunit, 100) comments(2)
      write (iunit, 1007) 'LBOUNDARY = ', lboundary
      write (iunit, 1007) 'LBMN = ', lbmn
      write (iunit, 1007) 'LJ_STAR = ', lj_star
      write (iunit, 1007) 'LASPECT_MAX = ', laspect_max
      write (iunit, 1007) 'LBETA_MIN = ', lbeta_min
      write (iunit, 1007) 'LJ_INVARIANT = ', lj_invariant
      write (iunit, 1007) 'LIOTA_PROF_OPT = ', liota_prof_opt
      write (iunit, 1007) 'LCUR_PROF_OPT = ', lcur_prof_opt
      write (iunit, 1007) 'LCUR_OPT_EDGE0 = ', lcur_opt_edge0
      write (iunit, 1007) 'LBOOTSJ_OPT = ', lbootsj_opt
      write (iunit, 1007) 'LKINK_OPT = ', lkink_opt
      write (iunit, 1007) 'LBALLOON_OPT = ', lballoon_opt                          
      write (iunit, 1007) 'L_LEGENDRE = ', l_legendre                 !! LEGENDRE
      write (iunit, 1007) 'LDKES_OPT = ', ldkes_opt                   !! RHF
      write (iunit, 1007) 'LNEO_OPT = ', lneo_opt                      
      write (iunit, 1007) 'LDSUBR_OPT = ', ldsubr_opt                      
      write (iunit, 1007) 'LORBIT_OPT = ', lorbit_opt                      
      write (iunit, 1007) 'LPRES_PROF_OPT = ', lpres_prof_opt
      write (iunit, 1007) 'LPRES_OPT_EDGE0 = ', lpres_opt_edge0
      write (iunit, 1007) 'LPRES_OPT_EDGEGR0 = ', lpres_opt_edgegr0
      write (iunit, 1007) 'LDIAGNO_OPT = ', ldiagno_opt
      write (iunit, 1007) 'LNESCOIL_OPT = ', lnescoil_opt
      write (iunit, 1007) 'LCOIL_GEOM = ', lcoil_geom
      write (iunit, 990) 'SIGMA_BERR_AVG = ', sigma_berr_avg
      write (iunit, 990) 'SIGMA_BERR_MAX = ', sigma_berr_max
      write (iunit, '(2x,a,i3,a)', err=2000) 
     1   'LFIX_NTOR = ',ntord+ntor1d,'*F'
      do i = -ntord,ntord
        if (lfix_ntor(i)) then
          if (i.ge.-9 .and. i.le.9) then
            write(iunit,1005) i
          else
            write(iunit,1006) i
          end if
        end if
      end do  

      if (nextcur_opt .gt. 0) then
         write (iunit, '(2x,a10,30L2)', err=2000) 
     1      'LEXTCUR = ', lextcur(1:nextcur_max)

         call write_array(iunit, 'SIGMA_EXTCUR', 
     1                     sigma_extcur, nextcur_vmec)
         call write_array(iunit, 'TARGET_EXTCUR', 
     1                     target_extcur, nextcur_vmec)
         call write_array(iunit, 'OH_COEFS', 
     1                     oh_coefs, nextcur_vmec)
      endif

      write (iunit, 990) 'SIGMA_OH = ', sigma_oh

      write (iunit, 100) comments(3)
      write (iunit, '(2(2x,a,1pe14.6))')
     1    'R00_SCALE = ', r00_scale, 'B00_SCALE = ', b00_scale

      write (iunit, '(2(2x,a,1pe14.6))')
     1    'RGRID_MIN = ', rgrid_min, 'RGRID_MAX = ', rgrid_max
      write (iunit, '(2(2x,a,1pe14.6))')
     1    'ZGRID_MIN = ', zgrid_min, 'ZGRID_MAX = ', zgrid_max
      write (iunit, 100) comments(4)
      write (iunit, '(2(2x,a,i4))') 'MBOZ_OPT = ', mboz_opt, 
     1                              'NBOZ_OPT = ',nboz_opt

      write (iunit, 100) comments(5)
      write (iunit, 990) 'COIL_SEPARATION = ', coil_separation
      write (iunit, 990) 'TARGET_ASPECTRATIO = ', target_aspectratio
      write (iunit, 990) 'TARGET_BETA = ', target_beta
      write (iunit, 990) 'TARGET_EPLASMA = ', target_eplasma
      write (iunit, 990) 'TARGET_CURTOR = ', target_curtor
      write (iunit, 990) 'TARGET_JEDGE = ', target_jedge
      call write_array(iunit, 'TARGET_KINK', target_kink,
     1                 size(target_kink))
      write (iunit, 990) 'TARGET_MAXCURRENT = ', target_maxcurrent
      write (iunit, 990) 'TARGET_RMAX = ', target_rmax
      write (iunit, 990) 'TARGET_RMIN = ', target_rmin
      write (iunit, 990) 'TARGET_ZMAX = ', target_zmax
      write (iunit, 990) 'TARGET_ELLIPTICITY = ', target_ellipticity
      write (iunit, 990) 'TARGET_KAPPA = ', target_kappa
      write (iunit, 990) 'TARGET_FLUXP = ', target_fluxp
      write (iunit, 990) 'TARGET_RBTOR = ', target_RBtor
      write (iunit, 990) 'TARGET_COIL_COMPLEX = ', target_coil_complex
      write (iunit, 990) 'TARGET_COIL_JMAX = ', target_coil_jmax
      call write_array(iunit, 'TARGET_IOTA', target_iota, 11)
      call write_array(iunit, 'TARGET_IOTA_P', target_iota_p, 11)
      write (iunit, 990) 'TARGET_IOTA_MIN = ', target_iota_min
      write (iunit, 990) 'TARGET_IOTA_MAX = ', target_iota_max
      write (iunit, 990) 'TARGET_IOTA_MAX_MIN = ', target_iota_max_min

      write (iunit, 990) 'SIGMA_ASPECT = ', sigma_aspect
      write (iunit, 990) 'SIGMA_BETA = ', sigma_beta
      write (iunit, 990) 'SIGMA_EPLASMA = ', sigma_eplasma
      write (iunit, 990) 'SIGMA_CURTOR = ', sigma_curtor
      write (iunit, 990) 'SIGMA_JEDGE = ', sigma_jedge
      write (iunit, 990) 'SIGMA_CURV = ', sigma_curv
      write (iunit, 990) 'SIGMA_COIL_COMPLEX = ', sigma_coil_complex
      write (iunit, 990) 'SIGMA_COIL_JMAX = ', sigma_coil_jmax
      write (iunit, 990) 'SIGMA_BERR_AVE = ', sigma_berr_ave
      call write_array(iunit, 'SIGMA_KINK', sigma_kink,size(sigma_kink))
      write (iunit, 990) 'SIGMA_MAXCURRENT = ', sigma_maxcurrent       
      write (iunit, 990) 'SIGMA_RMAX = ', sigma_rmax
      write (iunit, 990) 'SIGMA_RMIN = ', sigma_rmin
      write (iunit, 990) 'SIGMA_ZMAX = ', sigma_zmax
      write (iunit, 990) 'SIGMA_ELLIPTICITY = ', sigma_ellipticity
      write (iunit, 990) 'SIGMA_KAPPA = ', sigma_kappa
      write (iunit, 990) 'SIGMA_FLUXP = ', sigma_fluxp
      write (iunit, 990) 'SIGMA_RBTOR = ', sigma_RBtor
      write (iunit, 990) 'SIGMA_PSEUDO = ', sigma_pseudo
      write (iunit, 990) 'SIGMA_PSEUDO2 = ', sigma_pseudo2
      write (iunit, 1007) 'LPSEUDO_SIN = ', lpseudo_SIN
      call write_array(iunit, 'SIGMA_IOTA', sigma_iota, ns0)
      call write_array(iunit,'SIGMA_IOTA_PMAX',sigma_iota_pmax, ns0)
      call write_array(iunit,'SIGMA_IOTA_PMIN',sigma_iota_pmin, ns0)
      write (iunit, 990) 'SIGMA_IOTA_MIN = ', sigma_iota_min
      write (iunit, 990) 'SIGMA_IOTA_MAX = ', sigma_iota_max
      write (iunit, 990) 'SIGMA_IOTA_MAX_MIN = ', sigma_iota_max_min
      call write_array(iunit, 'SIGMA_MERCIER', sigma_mercier, ns0)
      call write_array(iunit, 'SIGMA_JAC', sigma_jac, size(sigma_jac))
      write (iunit,981) 'N_JAC = ', n_jac
      write (iunit,981) 'M_JAC = ', m_jac
      call write_array(iunit, 'SIGMA_VAC_ISLAND', sigma_vac_island, 
     1                 size(sigma_vac_island))
      write (iunit,981) 'N_VAC_ISLAND = ', n_vac_island
      write (iunit,981) 'M_VAC_ISLAND = ', m_vac_island

      write (iunit, 100) comments(6)
      write (iunit, '(2x,a12,i3,a2,i3,a1)') 'HELICITY = (',
     1  nint(real(helicity)),', ',nint(aimag(helicity)),')'
      call write_array(iunit,'SIGMA_BMIN', sigma_bmin, ns0)
      call write_array(iunit,'SIGMA_BMAX', sigma_bmax, ns0)
      call write_array(iunit,'SIGMA_BMN',sigma_bmn, ns0)
      call write_array(iunit,'SIGMA_RIPPLE',sigma_ripple, ns0)

      write (iunit, '(2(2x,a,i4))') 'NUMJSTAR = ',NumJstar, 
     1   ' NUMJINVARIANT = ', NumJinvariant
      do k = 1, NumJstar
         if (k .lt. 10) form = form1
         if (k .ge. 10) form = form2
         write (temp_write, form, err=2000) 'SIGMA_JSTAR(1:,',k,')'
         call write_array(iunit,trim(temp_write),sigma_jstar(1,k),ns0)
      end do
      do k = 1, NumJinvariant
         if (k .lt. 10) form = form1
         if (k .ge. 10) form = form2
         write(temp_write, form, err=2000)'SIGMA_JINVARIANT(1:,',k,')'
         call write_array(iunit, trim(temp_write),
     1        sigma_jinvariant(1,k), ns0)
      end do

      write (iunit, '(2(2x,a,i4))') 'NS_JCONF_SRC = ',NS_JConf_Src,
     1   'NS_JCONF_TGT = ',  NS_JConf_Tgt
      write(iunit,981) 'NPITCH_JCONF = ', NPitch_JConf
      write(iunit,990) 'SIGMA_JCONF = ',sigma_jconf

!----------------------------------------------------------------------------------------
!     CODE added by R.SANCHEZ (01/19/99) to write out ballooning-related info.
!     Modified (02/01/99) to include multiple initial position output.
!
      write (iunit,100) comments(7)
      if(lballoon_opt .or. lpres_prof_opt) then
        call write_array(iunit,'NBALLOON_MASK',nballoon_mask, ns0)
        call write_array(iunit,'TARGET_BALLOON',target_balloon, ns0)
        write (iunit, 1007) 'LBALLOON_FLIP = ', lballoon_flip
        call write_array(iunit,'SIGMA_BALLOON',sigma_balloon, ns0)
        call write_array(iunit,'SIGMA_PGRAD',sigma_pgrad, ns0)
        write (iunit, 990) 'SIGMA_PEDGE = ', sigma_pedge
        call write_array(iunit,'BAL_ZETA0',bal_zeta0, nini_zeta)
        call write_array(iunit,'BAL_THETA0', bal_theta0, nini_theta)
      endif
!-----------------------------------------------------------------------------------------      

      write (iunit, 100) comments(8)
      if (lbootsj_opt) then
         call write_array(iunit,'FBOOT',fboot, 11)
         if( jboot .ne. 0) write (iunit, '(2x,a,i1)') 'JBOOT = ',jboot
         if( lseedcur ) 
     1      call write_array(iunit,'ASEEDCUR',aseedcur, 11)
         call write_array(iunit,'SIGMA_BOOTSJ',sigma_bootsj, ns0)
      end if

      write (iunit, 100) comments(9)                                     !RHF
      call write_array(iunit, 'NDKES_MASK', ndkes_mask, ns0)        !RHF
      call write_array(iunit, 'DKES_NU', dkes_nu, ns0)              !RHF
      call write_array(iunit, 'DKES_EFIELD', dkes_efield, ns0)      !RHF
      call write_array(iunit,'SIGMA_DKES', sigma_dkes, ns0)         !RHF

c ---------------------------------------------------------------------------
!  NEO
      if (lneo_opt) then
        call write_array(iunit, 'NNEO_MASK', nneo_mask, ns0)
        call write_array(iunit, 'SIGMA_NEO', sigma_neo, ns0)
      end if

      if (ldsubr_opt) 
     1   call write_array(iunit, 'SIGMA_DSUBR', sigma_dsubr, ns0)

      if (lorbit_opt) 
     1   write (iunit, 990) 'SIGMA_ORBIT = ', sigma_orbit

c ---------------------------------------------------------------------------
c     write (iunit, 990) 'SIGMA_P_DAMP = ', sigma_p_damp
      write(iunit,981) 'NP_PROF = ', NP_PROF
      if( np_prof > 0 .and. any(sigma_p_prof > 0)) then
         write (iunit, 990) 'FACTOR_P_PROF = ', factor_p_prof
         write (iunit, 1007) 'LP_PROF_INCL_EDGE = ', lp_prof_incl_edge
         call write_array(iunit,'P_PROF', p_prof, np_prof)
         call write_array(iunit,'SIGMA_P_PROF', sigma_p_prof, np_prof)
         call write_array(iunit,'R_P_PROF', r_p_prof, np_prof)
         call write_array(iunit,'Z_P_PROF', z_p_prof, np_prof)
         call write_array(iunit,'PHI_P_PROF', phi_p_prof, np_prof)
c        write (iunit, 1007) 'LP_PROF_ABS = ', lp_prof_abs
      endif
c ---------------------------------------------------------------------------
c  MCZ  Diagno simulations of magnetic diagnostics
c
      if( ldiagno_opt) then
         write (iunit,'(2x,3a)') "DIAGNO_CONTROL = '",
     1                            trim(diagno_control),"'"

         if( ndiagno_seg > 0) then
            call write_array(iunit,'SIGMA_DIAGNO_SEG', sigma_diagno_seg, 
     1                       ndiagno_seg)
            call write_array(iunit,'TARGET_DIAGNO_SEG', 
     1                       target_diagno_seg, ndiagno_seg)
         endif
         
         if( ndiagno_flx > 0) then
            call write_array(iunit,'SIGMA_DIAGNO_FLX', sigma_diagno_flx, 
     1                       ndiagno_flx)
            call write_array(iunit,'TARGET_DIAGNO_FLX', 
     1                       target_diagno_flx, ndiagno_flx)
         endif
         
         
      endif

c ---------------------------------------------------------------------------
c     write(iunit,981) 'N_EMIS = ', NP_EMIS
c     if( n_emis > 0 .and. any(sigma_emis > 0)) then
c        call write_array(iunit,'SIGMA_EMIS', sigma_emis, n_emis)
c        write (iunit, 990) 'SIGMA_EMIS_DAMP = ', sigma_emis_damp
c        call write_array(iunit,'AEMIS', aemis, 11)
c        write (iunit,'(2x,3a)') "EMIS_FILE = '",trim(emis_file),"'"
c        call write_array(iunit,'EMIS_CHORD', emis_chord, n_emis)
c     endif
c ---------------------------------------------------------------------------
      if (sigma_vv < bigno .or. sigma_vv_rms < bigno .or.
     1    sigma_vv_max < bigno ) then
         write (iunit, '(a,1pe13.6)') '  TARGET_VV = ',target_vv
         write (iunit, 1007) 'LVV_TGT_MIN = ', lvv_tgt_min
         write (iunit,'(a,1pe13.6)' ) '  TARGET_VV_RMS = ',target_vv_rms
         write (iunit, '(a,1pe13.6)') '  SIGMA_VV = ',sigma_vv
         write (iunit,'(a,1pe13.6)' ) '  SIGMA_VV_RMS = ',sigma_vv_rms
         write (iunit,'(a,1pe13.6)' ) '  SIGMA_VV_MAX = ',sigma_vv_max

         write (iunit,'(a,i3)' ) '  MPOL_VV = ', mpol_vv
         write (iunit,'(a,i3)' ) '  NTOR_VV = ', ntor_vv
         if( mpol_vv > 0 .or. ntor_vv > 0) then
            write (iunit,'(a,i3)' ) '  NU_VV = ', nu_vv
            write (iunit,'(a,i3)' ) '  NV_VV = ', nv_vv
            do m = 0, mpol_vv
               do n = -ntor_vv, ntor_vv
                 if (rbc_vv(n,m)==zero .and. zbs_vv(n,m)==zero) cycle
                 if (n < -9 .and. m >= 10) then
                   write (iunit,'(2(A,i3,A,i2,A,1pe13.6,3x))',err=2000) 
     1                   '  RBC_VV(', n, ',', m, ') = ', rbc_vv(n,m),
     2                   '  ZBS_VV(', n, ',', m, ') = ', zbs_vv(n,m)
                 else if (n < -9 .and. m < 10) then
                   write (iunit,'(2(A,i3,A,i1,A,1pe13.6,3x))',err=2000) 
     1                   '  RBC_VV(', n, ',', m, ') = ', rbc_vv(n,m),
     2                   '  ZBS_VV(', n, ',', m, ') = ', zbs_vv(n,m)
                 else if ((n < 0 .or. n > 9).and. m >= 10) then
                   write (iunit,'(2(A,i2,A,i2,A,1pe13.6,3x))',err=2000) 
     1                   '  RBC_VV(', n, ',', m, ') = ', rbc_vv(n,m),
     2                   '  ZBS_VV(', n, ',', m, ') = ', zbs_vv(n,m)
                 else if ((n < 0 .or. n > 9).and. m < 10) then
                   write (iunit,'(2(A,i2,A,i1,A,1pe13.6,3x))',err=2000) 
     1                   '  RBC_VV(', n, ',', m, ') = ', rbc_vv(n,m),
     2                   '  ZBS_VV(', n, ',', m, ') = ', zbs_vv(n,m)
                 else if ( m < 10 ) then
                   write (iunit,'(2(A,i1,A,i1,A,1pe13.6,3x))',err=2000) 
     1                   '  RBC_VV(', n, ',', m, ') = ', rbc_vv(n,m),
     2                   '  ZBS_VV(', n, ',', m, ') = ', zbs_vv(n,m)
                 else
                   write (iunit,'(2(A,i1,A,i2,A,1pe13.6,3x))',err=2000) 
     1                   '  RBC_VV(', n, ',', m, ') = ', rbc_vv(n,m),
     2                   '  ZBS_VV(', n, ',', m, ') = ', zbs_vv(n,m)
               
                  endif
               end do
            end do

            write (iunit,'(a,L1)' ) '  SHAPEWEIGHT = ', shapeweight
            if( shapeweight ) then
                write (iunit, '(a,1p3e13.6)') '  PLANES_BW = ',
     1                    (planes_bw(n),n=1,3)
                write (iunit, '(a,1p3e13.6)') '  AMPLW_BW = ',
     1                    (amplw_bw(n),n=1,3)
                write (iunit, '(a,1p3e13.6)') '  THETA0_BW = ',
     1                    (theta0_bw(n),n=1,3)
                write (iunit, '(a,1pe13.6)') '  PHI0_BW = ',phi0_bw
                write (iunit, '(a,1pe13.6)') '  WTHETA_BW = ', wtheta_bw
                write (iunit, '(a,1pe13.6)') '  WPHI_BW = ', wphi_bw
            endif
         endif

         do n = 1, size(phi_lim)
            if( phi_lim(n) < -360 .or. phi_lim(n) > 360 .or.
     1          r_lim(1,n) == 0 .or. r_lim(2,n) == 0) exit

            m = 2
            do i = 3, size(r_lim, 1)
                if( r_lim(i,n) == 0) exit
                m = i
            enddo

            write(temp_write,*) n
            temp_write = adjustl(temp_write)
            write(iunit,'(a,a,a,1pe13.6)') ' PHI_LIM(', 
     1               trim(temp_write),') = ',phi_lim(n)
            temp_write = '(1:,'//trim(temp_write)//')'
            call write_array(iunit, 'R_LIM'//trim(temp_write),
     1                       r_lim(1,n), m)
            call write_array(iunit, 'Z_LIM'//trim(temp_write),
     1                       z_lim(1,n), m)
         enddo
      endif

c ---------------------------------------------------------------------------
      if (sigma_bd < bigno .or. sigma_bd_rms < bigno .or.
     1    sigma_bd_max < bigno ) then
         write (iunit, '(a,1pe13.6)') '  TARGET_BD = ',target_bd
         write (iunit,'(a,1pe13.6)' ) '  TARGET_BD_RMS = ',target_bd_rms
         write (iunit, '(a,1pe13.6)') '  SIGMA_BD = ',sigma_bd
         write (iunit,'(a,1pe13.6)' ) '  SIGMA_BD_RMS = ',sigma_bd_rms
         write (iunit,'(a,1pe13.6)' ) '  SIGMA_BD_MAX = ',sigma_bd_max

         write (iunit,'(a,i3)' ) '  MPOL_BD = ', mpol_bd
         write (iunit,'(a,i3)' ) '  NTOR_BD = ', ntor_bd
         if( mpol_bd > 0 .or. ntor_bd > 0) then
            write (iunit,'(a,i3)' ) '  NU_BD = ', nu_bd
            write (iunit,'(a,i3)' ) '  NV_BD = ', nv_bd
            do m = 0, mpol_bd
               do n = -ntor_bd, ntor_bd
                 if (rbc_bd(n,m)==zero .and. zbs_bd(n,m)==zero) cycle
                 if (n < -9 .and. m >= 10) then
                   write (iunit,'(2(A,i3,A,i2,A,1pe13.6,3x))',err=2000) 
     1                   '  RBC_BD(', n, ',', m, ') = ', rbc_bd(n,m),
     2                   '  ZBS_BD(', n, ',', m, ') = ', zbs_bd(n,m)
                 else if (n < -9 .and. m < 10) then
                   write (iunit,'(2(A,i3,A,i1,A,1pe13.6,3x))',err=2000) 
     1                   '  RBC_BD(', n, ',', m, ') = ', rbc_bd(n,m),
     2                   '  ZBS_BD(', n, ',', m, ') = ', zbs_bd(n,m)
                 else if ((n < 0 .or. n > 9).and. m >= 10) then
                   write (iunit,'(2(A,i2,A,i2,A,1pe13.6,3x))',err=2000) 
     1                   '  RBC_BD(', n, ',', m, ') = ', rbc_bd(n,m),
     2                   '  ZBS_BD(', n, ',', m, ') = ', zbs_bd(n,m)
                 else if ((n < 0 .or. n > 9).and. m < 10) then
                   write (iunit,'(2(A,i2,A,i1,A,1pe13.6,3x))',err=2000) 
     1                   '  RBC_BD(', n, ',', m, ') = ', rbc_bd(n,m),
     2                   '  ZBS_BD(', n, ',', m, ') = ', zbs_bd(n,m)
                 else if ( m < 10 ) then
                   write (iunit,'(2(A,i1,A,i1,A,1pe13.6,3x))',err=2000) 
     1                   '  RBC_BD(', n, ',', m, ') = ', rbc_bd(n,m),
     2                   '  ZBS_BD(', n, ',', m, ') = ', zbs_bd(n,m)
                 else
                   write (iunit,'(2(A,i1,A,i2,A,1pe13.6,3x))',err=2000) 
     1                   '  RBC_BD(', n, ',', m, ') = ', rbc_bd(n,m),
     2                   '  ZBS_BD(', n, ',', m, ') = ', zbs_bd(n,m)
               
                  endif
               end do
            end do
         endif
      endif

      write (iunit, '(a)') '/'
      
  100 format('!',40('-'),/,'!',5x,a,/,'!',40('-'))
 1010 format(8(1pe16.6))

      write(iunit, nml=bootin, err=2000)
      call write_gade_nml(iunit)

      if (lcoil_geom .or. sigma_berr_avg < bigno .or. 
     1    sigma_berr_max < bigno) call write_coilsin (iunit, istat)
      if (istat .ne. 0) goto 2000

      close (iunit)
      
      return
      
!
!     Handle errors here      
!
 2000 istat = -5
      print *,'write_indata io error '
      call vmec_flush(6)
     
      end subroutine write_indata


      subroutine write_rbzb(iunit, istat)
      use optim
      use vmec_input
      use vparams, only: zero
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: iunit, istat
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      character*(*), dimension(7), parameter :: outfmt =
     1   (/ '(2(A6,i1,A1,i1,A4,1pe22.14,3x))', 
     2      '(2(A6,i1,A1,i2,A4,1pe21.14,3x))',
     3      '(2(A6,i2,A1,i1,A4,1pe21.14,3x))',
     4      '(2(A6,i2,A1,i2,A4,1pe21.14,3x))',
     5      '(2(A6,i3,A1,i1,A4,1pe21.14,3x))',
     6      '(2(A6,i3,A1,i2,A4,1pe21.14,3x))',
     7      '(2(A6,i ,A1,i ,A4,1pe21.14,3x))' /)
      integer, parameter :: n_max_leg = 6                               !! LEGENDRE (uses 6+1coeff.)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, n, m
      integer(kind=iprec):: js                                          !! LEGENDRE
      real(rprec):: sum1                                                !! LEGENDRE
      real(rprec), allocatable, dimension(:) :: ac_leg, ai_leg          !! LEGENDRE
      real(rprec), dimension(nrad) :: x_rad                             !! LEGENDRE
      character*(len(outfmt(1))) :: outcfmt
C-----------------------------------------------
!
!     CONVERT MESH QUANTITIES (iota_opt, jcurv_opt) TO S-EXPANSION COEFFICIENTS
!
      if (lcur_prof_opt .and. nrad.gt.3) then                                              
        x_rad(2:nrad) = (/( (real(js, rprec)-0.5_dp)                     
     1    /real(nrad-1, rprec), js=1,nrad-1)/)                           
        allocate(ai_leg(0:n_max_leg), stat=istat)
        if (istat .ne. 0) goto 2000                                    
        call point_wise_to_power(nrad-1, x_rad(2), iota_opt(2),          
     1    n_max_leg, ai_leg, 1, 1000)                                    
                                                                         
        ai = 0                                                            
        do js = 0, n_max_leg                                             
          ai(js) = ai_leg(js)                                            
        enddo                                                            
        deallocate(ai_leg)                                               
                                                                         
      else if (liota_prof_opt .and. nrad.gt.3) then                     !(JCURV_OPT on full-mesh, NOT integrated in s)                 
                                                                         
        x_rad(1:nrad) = (/( real(js-1, rprec)/(nrad-1), js=1,nrad)/)
        allocate(ac_leg(0:n_max_leg), stat=istat)
        if (istat .ne. 0) goto 2000                                    
        call point_wise_to_power(nrad, x_rad, jcurv_opt,            
     1    n_max_leg, ac_leg, 1, 1000)                                    
                                                                         
        ac = 0                                                           
        do js = 0, n_max_leg                                             
          ac(js) = ac_leg(js)                                            
        enddo
        ac_form = 0
        deallocate(ac_leg)                                               
      endif                                                              
!
!     BE SURE TO WRITE OUT ORIGINAL RAXIS, ZAXIS (THEY MAY BE OVERRIDDEN IN READ_WOUT)
!     THIS WILL PRESERVE CHI-SQ IN .MIN FILE
!
      write (iunit, 100) '  CURTOR = ',curtor
      write (iunit, 100) '  SPRES_PED = ',spres_ped
      write (iunit, 100, err=2000) '  AM = ', (am(n-1), n=1,size(am))
      write (iunit, 100, err=2000) '  AI = ', (ai(n-1), n=1,size(ai))
      write (iunit, 101) '  AC_FORM = ',ac_form
      write (iunit, 100, err=2000) '  AC = ', (ac(n-1), n=1,size(ac))
!     if (any(aphi .ne. zero))
!    1  write (iunit, 100, err=2000) '  APHI = ', (aphi(n), n=0,10)
! write out both the new axis and the original one, so the user can easily
! try either in case of convergence problems.  The original axis is written
! second, so it will have precedence. NOTE: RBC, ZBS WILL BE REASSIGNED TO
! THE BOUNDARY VALUES IN READ_WOUT_OPT IF LFREEB = T

      write (iunit, 100, err=2000) '  RAXIS = ',(raxis(n,1), n=0,ntor)
      write (iunit, 100, err=2000) '  ZAXIS = ',(zaxis(n,1), n=0,ntor)
      write (iunit, 100, err=2000) '  RAXIS = ',(raxis_old(n), n=0,ntor)
      write (iunit, 100, err=2000) '  ZAXIS = ',(zaxis_old(n), n=0,ntor)
      do m = 0, mpol - 1
         do n = -ntor, ntor
            if ((rbc(n,m).ne.zero) .or. (zbs(n,m).ne.zero)) then
!     
!     handle formatting for up to 2 digit n by 3 digit m.
!     while this is probably overkill, we at least have to handle up
!     thru 2 digit n by 2 digit m to handle known cases.
!
               outcfmt = outfmt(1)
               if( m > 9) outcfmt = outfmt(2)
               if(( n>-10 .and. n<0 ) .or. (n > 9 .and. n < 100)) then
                  outcfmt = outfmt(3)
                  if( m > 9) outcfmt = outfmt(4)
               else if( n>-100 .and. n< -9 ) then
                  outcfmt = outfmt(5)
                  if( m > 9) outcfmt = outfmt(6)
               endif

               if( n>= 100 .or. n <= -100 .or. m >=100) then
                  outcfmt = outfmt(7)
               endif
                  
               write (iunit, outcfmt, err=2000) 
     1            '  RBC(', n, ',', m, ') = ', rbc(n,m),
     2            '  ZBS(', n, ',', m, ') = ', zbs(n,m)
            endif
         end do
      end do

 100  format(a,(1p4e22.14))
 101  format(a,i)

      istat = 0
      return

 2000 istat = -5

      end subroutine write_rbzb
            

      subroutine write_neoparam (no_fluxs, fluxs_arr, extension, 
     1    write_progress, istat)
      use safe_open_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: no_fluxs, fluxs_arr(*), istat, write_progress
      character*(*) :: extension
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, parameter :: max_m_mode = 0, max_n_mode = 0     !LET NEO COMPUTE INTERNALLY
      integer, parameter :: write_output_files = 0, write_integrate = 0,
     1     write_diagnostic = 0, calc_cur = 0, write_cur_inte = 0      
      integer :: w_neo = 10
      character*200 :: arg1
C-----------------------------------------------
      
      arg1 = "neo_in." // trim(extension)

      call safe_open(w_neo, istat, arg1, 'replace', 'formatted')
      if (istat .ne. 0) return

      WRITE (w_neo,*) "#"
      WRITE (w_neo,*) "#"
      WRITE (w_neo,*) "#"
      WRITE (w_neo,'(a)') "boozmn." // trim(extension)
      WRITE (w_neo,'(a)') "neo_out." // trim(extension)

      WRITE (w_neo,*) no_fluxs

      IF (no_fluxs .LE. 0) THEN    
         WRITE (w_neo,*) "#"
      ELSE
         WRITE (w_neo,*) fluxs_arr(1:no_fluxs)
      END IF

      WRITE (w_neo,*) 100
      WRITE (w_neo,*) 100
      WRITE (w_neo,*) max_m_mode
      WRITE (w_neo,*) max_n_mode
      WRITE (w_neo,*) 50
      WRITE (w_neo,*) 1
      WRITE (w_neo,*) 0.01
      WRITE (w_neo,*) 100
      WRITE (w_neo,*) 50
      WRITE (w_neo,*) 500
      WRITE (w_neo,*) 3000
      WRITE (w_neo,*) 0
      WRITE (w_neo,*) 1
      WRITE (w_neo,*) 0
      WRITE (w_neo,*) 0
      WRITE (w_neo,*) 2
      WRITE (w_neo,*) write_progress
      WRITE (w_neo,*) write_output_files
      WRITE (w_neo,*) 0
      WRITE (w_neo,*) write_integrate
      WRITE (w_neo,*) write_diagnostic
      WRITE (w_neo,*) "#"
      WRITE (w_neo,*) "#"
      WRITE (w_neo,*) "#"
      WRITE (w_neo,*) calc_cur
      WRITE (w_neo,*) "Bmns_cur.dat"
      WRITE (w_neo,*) 200
      WRITE (w_neo,*) 2
      WRITE (w_neo,*) write_cur_inte

      CLOSE (unit=w_neo)

      end subroutine write_neoparam


      subroutine bextrema(modb, bmin, bmax, nzeta, ntheta)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: nzeta, ntheta
      real(rprec), intent(in) :: modb(nzeta,ntheta)
      real(rprec), intent(out) :: bmin(ntheta), bmax(ntheta)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ku
C-----------------------------------------------
*
*     Computes max, min of |B| along v (zeta) between two angle lines (theta = 0,pi)
*
      do ku = 1,ntheta
         bmin(ku)  = minval(modb(:,ku))
         bmax(ku)  = maxval(modb(:,ku))
      enddo

      end subroutine bextrema


      subroutine j_star(modb, bmin, bmax, ep_mu, jstar, nzeta, ntheta)
      use kind_spec
      use mpi_params                                                     !MPI
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nzeta, ntheta
      real(rprec), intent(in) :: ep_mu
      real(rprec), dimension(nzeta,ntheta), intent(in) :: modb
      real(rprec), dimension(ntheta), intent(in) :: bmin, bmax
      real(rprec), dimension(ntheta), intent(out) :: jstar
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, one = 1
      integer :: ku_lw, ku_up, ku, istat
      real(rprec) :: dzeta
      real(rprec), dimension(ntheta) :: test_bmin, test_bmax
      real(rprec) , dimension(ntheta) :: test1, test2
      real(rprec) , allocatable,  dimension(:,:) :: vpar1
      logical, allocatable, dimension(:,:) :: l3v
C-----------------------------------------------
*
*     Computes the trapped branch of Jstar on a single flux surface
*     and for a single value of ep/mu.  Jstar is only non-zero for
*     values of theta where trapped particle orbits can exist.  At
*     theta values which are in the passing particle regime (ep/mu > Bmax)
*     or the forbidden regime (ep/mu < Bmin) Jstar is set to 0.  A Jstar
*     topology is assumed here such that as theta runs from 0 to Pi,
*     one first encounters the ep/mu = Bmax point and then
*     the ep/mu = Bmin point.
 
      ku_lw = 1
      ku_up = ntheta
      jstar = zero
      dzeta = one
      if (nzeta .gt. 1) dzeta = one/real((nzeta-1),rprec)
      test_bmax = ep_mu - bmax(:ntheta)
      test_bmin = ep_mu - bmin(:ntheta)

      test1(2:ntheta) = test_bmax(2:ntheta)*test_bmax(:ntheta-1)
      test2(2:ntheta) = test_bmin(2:ntheta)*test_bmin(:ntheta-1)
      do ku = 2,ntheta
          if(test1(ku) .le. zero) ku_lw = ku
          if(test2(ku) .le. zero) ku_up = ku
      end do
 
      if (ku_lw .ge. ku_up) return
      
      allocate (vpar1(nzeta,ku_up-ku_lw+1), l3v(nzeta,ku_up-ku_lw+1),
     1   stat=istat)
      if (istat .ne. 0) then
         if (myid .eq. master) print *,' Allocation error in J_STAR!'    !MPI               
         return
      end if   

      vpar1 = one - modb(:,ku_lw:ku_up)/ep_mu
      l3v = vpar1 .gt. zero
      where (l3v) vpar1 = sqrt(vpar1)/modb(:,ku_lw:ku_up)
      jstar(ku_lw:ku_up) = dzeta * sum(vpar1, mask=l3v, dim=1)

      deallocate (vpar1, l3v)

      end subroutine j_star


      subroutine modbooz(bmnb, bmod, ixmb, ixnb, nfp, icontrol)
      use boozer_params
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: nfp, icontrol, ixmb(*), ixnb(*)
      real(rprec), dimension(mnboz), intent(in) :: bmnb
!    1 , rmnb, zmnb, pmnb
      real(rprec), dimension(nunv), intent(out) :: bmod
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: mn, m, n, i, lk, lz, lt, isgn
      real(rprec), parameter :: zero = 0.0_dp, one = 1.0_dp
      real(rprec), dimension(:), allocatable :: cost, sint
      real(rprec), dimension(:), allocatable :: thgrd, ztgrd
      real(rprec), dimension(:,:), allocatable, save :: 
     1     cosmm, sinmm, cosnn, sinnn 
      real(rprec) :: twopi, dth, dzt
C-----------------------------------------------
!
!     PERFORMS INVERSE TRANSFORM OF |BMN_B| INTO BOOZER ANGLE SPACE
!     CAN BE EASILY EXPANDED TO INVERSE TRANSFORM OTHER QUANTITIES (R, Z, p)
!
      twopi = 8*atan(one)

      if (.not. allocated(cosmm)) then
         allocate (cosmm(nunv,0:mboz), sinmm(nunv,0:mboz), 
     1             cosnn(nunv,0:nboz), sinnn(nunv,0:nboz), 
     2             thgrd(nunv), ztgrd(nunv), stat = i)
      if (i .ne. 0) stop 'allocation error in modbooz'
      
!
!     COMPUTE POLOIDAL (thgrd) AND TOROIDAL (ztgrd) ANGLES
!
!        dth = twopi/real(nu_boz,rprec)         !USE THIS FOR FULL 2-pi
         dth = twopi/real(2*(nu2_b-1),rprec)    !Half-around in theta
         dzt = twopi/real(nv_boz,rprec)
         lk = 0

         do lt = 1, nu2_b
            do lz = 1, nv_boz
               lk = lk + 1
               thgrd(lk) = (lt-1)*dth
               ztgrd(lk) = (lz-1)*dzt
            end do
         end do

         if (lk .ne. nunv) stop 'lk != nunv in modbooz'

         cosmm(:,0) = one
         sinmm(:,0) = zero
         cosmm(:,1) = cos(thgrd)
         sinmm(:,1) = sin(thgrd)

         cosnn(:,0) = one
         sinnn(:,0) = zero
         cosnn(:,1) = cos(ztgrd)
         sinnn(:,1) = sin(ztgrd)
          
         deallocate (thgrd, ztgrd) 

         do m = 2,mboz
            cosmm(:,m) = cosmm(:,m-1)*cosmm(:,1)
     1                 - sinmm(:,m-1)*sinmm(:,1)
            sinmm(:,m) = sinmm(:,m-1)*cosmm(:,1)
     1                 + cosmm(:,m-1)*sinmm(:,1)
         end do

         do n = 2,nboz
            cosnn(:,n) = cosnn(:,n-1)*cosnn(:,1)
     1                 - sinnn(:,n-1)*sinnn(:,1)
            sinnn(:,n) = sinnn(:,n-1)*cosnn(:,1)
     1                 + cosnn(:,n-1)*sinnn(:,1)
         end do
      end if
      
      allocate (cost(nunv), sint(nunv), stat = i) 

      bmod = zero

      do mn = 1,mnboz
        m = ixmb(mn)
        n = abs(ixnb(mn))/nfp
        if (m .gt. mboz) stop 'm > mboz in modbooz'
        if (n .gt. nboz) stop 'n > nboz in modbooz'
        isgn = sign(1,ixnb(mn))
        cost = cosmm(:,m)*cosnn(:,n)
     1       + sinmm(:,m)*sinnn(:,n)*isgn
c       sint = sinmm(:,m)*cosnn(:,n)
c    1       - cosmm(:,m)*sinnn(:,n)*isgn
        bmod = bmod + bmnb(mn)*cost
c       r12  = r12  + rmnb(mn)*cost
c       z12  = z12  + zmnb(mn)*sint
c       p12  = p12  + pmnb(mn)*sint
      end do

      deallocate (cost, sint)
      if (icontrol .eq. -1) deallocate
     1   (cosmm, sinmm, cosnn, sinnn)
       
      end subroutine modbooz


      subroutine chk_rzmnb(xcb_opt, rmnc_a, zmns_a, xm, xn, iflag)
      use optim
      use vmec_input, ONLY : rbc, zbs
      use vparams, ONLY: zero
      use mpi_params                                                     !MPI
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(out) :: iflag
      real(rprec), dimension(*) :: xcb_opt
      real(rprec), dimension(mnmax_opt), intent(in) :: 
     1   rmnc_a, zmns_a, xm, xn
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: eps = 1.0e-10_dp
      character*(*), parameter :: error_message =
     1   'Boundary array inconsistency in CHK_RZMNB'
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nb, mb, ikt, mn, ikrun, j1
      real(rprec) :: rzero, tol, compit
      logical :: lerror
C-----------------------------------------------
      iflag = 0;  lerror = .false.
      if (irm0_bdy.eq.0 .or. lniter1) return
      rzero = eps
      
!     rhobc = zero
      zbs(0,0) = zero

!
!     NOTE: THIS ASSUMES THAT R,Z ARE THE FIRST ELEMENTS IN XCB_OPT ARRAY
!
      do ikt = 1, irm0_bdy
         rbc(nrz0_opt(ikt),0) = xcb_opt(ikt)
      end do
      do ikt = 1, izm0_bdy
         j1 = ikt+irm0_bdy
         zbs(nrz0_opt(j1),0) = xcb_opt(j1)
      end do      
      do ikt = 1, irho_bdy
         j1 = ikt + irm0_bdy + izm0_bdy
         rhobc(nbrho_opt(ikt),mbrho_opt(ikt)) = xcb_opt(j1)
      end do
      
      call unique_boundary(rbc, zbs, rhobc, ntord, mpol1d, 
     1     mpol1_opt, ntor_opt)

      tol = sqrt(rzero)            
      compit = rzero + abs(rbc(0,0)) * tol

      do mn = 1,mnmax_opt
         nb = nint(xn(mn))/nfp_opt
         mb = nint(xm(mn))
         if(abs(rbc(nb,mb) - rmnc_a(mn)).gt.compit) then
            if (myid .eq. master)                                        !MPI
     1         write(*,57)nb,mb,rbc(nb,mb),rmnc_a(mn)
            lerror = .true.
         end if   
         if(abs(zbs(nb,mb) - zmns_a(mn)).gt.compit) then
           if (myid .eq. master)                                         !MPI      
     1        write(*,58)nb,mb,zbs(nb,mb),zmns_a(mn)
           lerror = .true.
         end if   
      enddo

 57   format(' nb = ',i5,' mb = ',i5,' rmnc diff: ',2e15.5)
 58   format(' nb = ',i5,' mb = ',i5,' zmns diff: ',2e15.5)

      if (lerror) then
         print *, trim(error_message),' in process ', myid
         write (iunit_opt_local,*) trim(error_message),
     1      ' in process ', myid
         iflag = -15
      end if   

      end subroutine chk_rzmnb


      function piota_t (x)
      use optim_params, only: rprec, target_iota
      real(rprec) :: x, piota_t
      piota_t = target_iota(0) + x*(target_iota(1)+x*(target_iota(2)+x*(
     1   target_iota(3)+x*(target_iota(4)+x*(target_iota(5)+x*(
     2   target_iota(6)+x*(target_iota(7)+x*(target_iota(8)+x*(
     3   target_iota(9)+x*target_iota(10))))))))))

      end function piota_t
      
      function piota_prime (x)
      use optim_params, only: rprec, target_iota_p
      real(rprec) :: x, piota_prime
      piota_prime = target_iota_p(0) + x*(target_iota_p(1)+
     1   x*(target_iota_p(2)+x*(
     1   target_iota_p(3)+x*(target_iota_p(4)+x*(target_iota_p(5)+x*(
     2   target_iota_p(6)+x*(target_iota_p(7)+x*(target_iota_p(8)+x*(
     3   target_iota_p(9)+x*target_iota_p(10))))))))))

      end function piota_prime
      
      function pwell_t (x)
      use optim_params, only: rprec, target_well
      real(rprec) :: x, pwell_t
      pwell_t = target_well(0) + x*(target_well(1)+x*(target_well(2)+x*(
     1   target_well(3)+x*(target_well(4)+x*(target_well(5)+x*(
     2   target_well(6)+x*(target_well(7)+x*(target_well(8)+x*(
     3   target_well(9)+x*target_well(10))))))))))

      end function pwell_t


      function pseed_t (x)
      use optim_params, only: rprec, aseedcur
      real(rprec) :: x, pseed_t
C-----------------------------------------------
      pseed_t = aseedcur(0) + x*(aseedcur(1) + x*(aseedcur(2) +
     1    x*(aseedcur(3) + x*(aseedcur(4) + x*(aseedcur(5) +
     2    x*(aseedcur(6) + x*(aseedcur(7) + x*(aseedcur(8) +
     3    x*(aseedcur(9) + x*aseedcur(10))))))))))
      end function pseed_t
      
      function pfboot (x)
      use optim_params, only: rprec, fboot
      real(rprec) :: x, pfboot
C-----------------------------------------------
      pfboot = fboot(0) + x*(fboot(1) + x*(fboot(2) +
     1      x*(fboot(3) + x*(fboot(4) + x*(fboot(5) +
     2      x*(fboot(6) + x*(fboot(7) + x*(fboot(8) +
     3      x*(fboot(9) + x* fboot(10))))))))))
      end function pfboot

      function eval_prof (a, x)
      use optim_params, only: rprec
      real(rprec) :: x, eval_prof
      real(rprec), dimension(0:10) :: a
C-----------------------------------------------
      eval_prof = a(0) + x*(a(1) + x*(a(2) +
     1      x*(a(3) + x*(a(4) + x*(a(5) +
     2      x*(a(6) + x*(a(7) + x*(a(8) +
     3      x*(a(9) + x* a(10))))))))))
      end function eval_prof

      subroutine flux_surf_curv(curv_kur, rmnc_a, zmns_a, xm, xn)
      use optim
      use safe_open_mod
      use vparams, only: zero, one
      use mpi_params                                                     !MPI
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(4) :: curv_kur
      real(rprec), dimension(mnmax_opt) :: rmnc_a, zmns_a, xm, xn
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: four = 4.0_dp
      logical, parameter :: lprint = .false.         !!Set .true. to get debugging info
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: lzth, l, ith, i, iunit1, iunit2, iunit3, istat
      real(rprec), dimension(100) :: r_major, z_elev
      real(rprec), dimension(100,4) :: curv
      real(rprec), dimension(4) :: 
     1   avg_curv, curv_dev2, curv_dev4, curv_max
      real(rprec) :: xpi, zeta_1, theta_1, arg, s_all, sa2, 
     1   sb2, s_x2, s_y2, delt_1, ratio
      character*15 :: filename
C-----------------------------------------------
 
!
!     Open up several diagnostic data files.
!
      if (lprint) then
         iunit1 = 34
         iunit2 = iunit1+1
         iunit3 = iunit2+1
      
         call safe_open(iunit1, istat, 'curv_dev', 'replace', 
     1      'formatted')
         call safe_open(iunit2, istat, 'outer_surf', 'replace', 
     1      'formatted')
      endif

      xpi = four*atan(one)
      curv_kur = zero
!
!     check curvature at phi=0,90,180,270 degree planes (lzth = 1,2,3,4):
!
      do 100 lzth = 1, 4
         avg_curv(lzth) = zero
         if (lprint) then
            if (lzth .eq. 1) filename = 'test_out_0'
            if (lzth .eq. 2) filename = 'test_out_90'
            if (lzth .eq. 3) filename = 'test_out_180'
            if (lzth .eq. 4) filename = 'test_out_270'
            call safe_open(iunit3, istat, filename, 'replace', 
     1         'formatted')
            if (istat.ne.0 .and. myid.eq.master)                         !MPI      
     1      print *, ' Error opening ',trim(filename),
     2               ' in FLUX_SURF_CURV'
         end if
         zeta_1 = ((lzth - 1)*xpi)/real((2*nfp_opt),rprec)
         if (lzth.eq.1 .and. lprint) then
            do l = 1, mnmax_opt
               write (iunit2, *) xm(l), xn(l), rmnc_a(l), zmns_a(l)
            end do
            close(iunit2)
         endif
c
c        Construct real space R, Z coordinates of outer flux surface:
c
         do ith = 1,100
            theta_1 = (2.0_dp*xpi*ith)/100.0_dp
            R_major(ith) = zero
            z_elev(ith) = zero
            do l=1,mnmax_opt
               arg = xm(l)*theta_1 - xn(l)*zeta_1
               R_major(ith) = R_major(ith) + rmnc_a(l)*cos(arg)
               z_elev(ith) = z_elev(ith) + zmns_a(l)*sin(arg)
            end do
         end do
c
c        Calculate local curvature around flux surface by finite
c        differencing the local unit tangent vectors:
c
         do i=2,99
            s_all = zero
            sa2 = sqrt((R_major(i)-R_major(i-1))**2 +
     1            (z_elev(i)-z_elev(i-1))**2)
            sb2 = sqrt((R_major(i+1)-R_major(i))**2 +
     1            (z_elev(i+1)-z_elev(i))**2)
            s_x2 =  ((R_major(i) - R_major(i-1))/sa2
     1          -   (R_major(i+1) - R_major(i))/sb2)**2
            s_y2 =  ((z_elev(i) - z_elev(i-1))/sa2
     1          - (z_elev(i+1) - z_elev(i))/sb2)**2
            s_all = sqrt(s_x2 + s_y2)
            delt_1 = .5_dp*sqrt((R_major(i-1) - R_major(i+1))**2
     1             + (z_elev(i-1) - z_elev(i+1))**2)

            curv(i,lzth) = zero
            if (delt_1 .ne. zero) curv(i,lzth) = s_all/delt_1
            avg_curv(lzth) = avg_curv(lzth) + curv(i,lzth)/98.0_dp
            if (lprint) write (iunit3,998) R_major(i),
     1         z_elev(i),curv(i,lzth)
         end do
         if (lprint) close(iunit3) 
 100  continue

      do lzth = 1,4
         curv_dev2(lzth) = zero
         curv_dev4(lzth) = zero
         curv_max(lzth) = -1.e30_dp
c        Calculate average curvature, variance of curvature, and kurtosis
c        (4th moment/2nd moment) around flux surface.
         do i=2,99
            curv_max(lzth) = max(curv_max(lzth),curv(i,lzth))
            curv_dev2(lzth) = curv_dev2(lzth) + ((curv(i,lzth)
     1                     - avg_curv(lzth))**2)/98.0_dp
            curv_dev4(lzth) = curv_dev4(lzth) + ((curv(i,lzth)
     1                     - avg_curv(lzth))**4)/98.0_dp
         end do
         curv_kur(lzth) = curv_dev4(lzth)/(curv_dev2(lzth)**2)
         if(lprint) then
            ratio = curv_max(lzth)/avg_curv(lzth)
            write(iunit1,997) lzth, avg_curv(lzth),
     1          curv_dev2(lzth), ratio, curv_kur(lzth)
         end if
      end do

  998 format(3(1x,e15.6))
  997 format(1x,i2,4(1x,e15.6))
      if (lprint) close(iunit1)

      end subroutine flux_surf_curv
 
! *************************************************************************
!     SURFSEP SUBS SUITE - nested surface separation subroutines
!     modified from SURFSEP4 by A. Brooks to compute the sign of
!     the distance to determine if the point on the free surface
!     lies inside the fixed surface Note that this could only be
!     compiled with optimization at -O 1 when it was a standalone
!     routine. This problem does not exist in this version so we
!     compile with -O (= -O 2) In this set of routines 2 refers to
!     the plasma surface while 1 refers to the vacuum vessel or
!     "fixed" surface.
!     R.E. Hatcher - May 2000
! *************************************************************************

      subroutine surfsep(mnmax_pl, r_pl, z_pl, xm_pl, xn_pl, mnmax_vv,
     1    nu_vv, nv_vv, r_vv, z_vv, xm_vv, xn_vv, nper, igrid, 
     2    grid_dist, distance, drms, dmax, ierr)
      use kind_spec
      use cmnf1
      use cmnf2
      use geom

      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer mnmax_pl, mnmax_vv, nper, igrid, nu_vv, nv_vv, ierr
      real(rprec) distance, drms, dmax
      real(rprec), dimension(mnmax_pl) :: r_pl, z_pl, xm_pl, xn_pl
      real(rprec), dimension(mnmax_vv) :: r_vv, z_vv, xm_vv, xn_vv
      real(rprec), dimension(igrid, igrid) :: grid_dist
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iperr,ngood2,i,isurf,ip,it,np,nt,ngood
      real(rprec) :: xun, yun, zun, davg, xvec, yvec, zvec, xf, yf,
     1   zf, r1, z1, d, dotp, du, dv, phi, u1, v1, u2, v2, r2, z2,
     2   dmean, std
      real(rprec), dimension(:), allocatable :: dmat
 
      integer :: jp, jt, kp, kt
      real(rprec) :: u1_sav, v1_sav, dist_test
      integer, dimension(1) :: maxlocation

!-----------------------------------------------
!
      ierr = 0                                   ! presume success to start
!
!     define some angular constants
      q1 = nper
      q2 = nper
      pi2 = 4*asin(1._dp)
      alp = pi2/nper
!
!     set np and nt to igrid
      np = igrid
      nt = igrid

      grid_dist = 1.e30_rprec
!
!     allocate an array for distances
      allocate(dmat(np*nt), stat=ierr)
!
!     compute the total number of modes for each surface
      nmn2 = mnmax_pl
      nmn1 = mnmax_vv
!
!     equate the common block variables with the arguments
!     Surface 1 - Vacuum Vessel
!
      allocate (xm1(nmn1), xn1(nmn1), rmn1(nmn1), zmn1(nmn1),
     1          xm2(nmn2), xn2(nmn2), rmn2(nmn2), zmn2(nmn2),stat=ierr)
      if (ierr .ne. 0) return
!
!
      xm1(:nmn1) = xm_vv(:nmn1)
      xn1(:nmn1) = xn_vv(:nmn1)
      rmn1(:nmn1) = r_vv(:nmn1)
      zmn1(:nmn1) = z_vv(:nmn1)
!
!     equate the common block variables with the arguments
!     Surface 2 - Plasma
!
      xm2(:nmn2) = xm_pl(:nmn2)
      xn2(:nmn2) = xn_pl(:nmn2)
      rmn2(:nmn2) = r_pl(:nmn2)
      zmn2(:nmn2) = z_pl(:nmn2)
!
!     write (6,*) 'done copying data', nmn2,nmn1
!
!     foreach each point on plasfree,calculate dist to plas fix
      distance = 1.0e30_dp
      ngood = 0

      du = 1.0_dp/np
      dv = 1.0_dp/nt
      do ip = 1, np
         do it = 1, nt
            u2 = (ip - 1)*du
            v2 = (it - 1)*dv
!           ! get point (r2,z2) on plasma surface corrseponding to u2,v2
            call dmnf2 (u2, v2, r2, z2)
            phi = alp*v2
            xp = r2*cos(phi)
            yp = r2*sin(phi)
            zp = z2
!           ! assume that u2,v2 is a good starting estimate for u1,v1
            u1 = u2
            v1 = v2

!           un-vectorized loop to find a better starting u1, v1       LPK-011002

            u1_sav = u1
            v1_sav = v1
            dist_test=1.0e30_dp

!           The following neighborhood size should be made as part of input
!          
            do jt=-nv_vv, nv_vv
               do jp=-nu_vv, nu_vv
                  kt = it + jt
                  kp = ip + jp
                  if ( kt .le. 0 ) kt = kt + nt
                  if ( kp .le. 0 ) kp = kp + np
                  u1 = (kp -1)*du
                  v1 = (kt -1)*dv
                  call dmnf1( u1, v1, r1, z1 )
                  phi = alp*v1
                  xf = r1*cos(phi)
                  yf = r1*sin(phi)
                  zf = z1
                  d  = (xf-xp)**2+(yf-yp)**2+(zf-zp)**2
                  if( d .lt. dist_test ) then
                     u1_sav = u1
                     v1_sav = v1
                     dist_test = d
                  end if
               end do
            end do
            u1 = u1_sav
            v1 = v1_sav
            d  = dist_test
!
!           end initial search of minimum distance

!           ! find the point of a minimum distance from the plasma
!           ! to the vacuum vessel
            iperr = 0
            call p2surf (u1, v1, d, iperr)
            if (iperr == 1) cycle                ! could not find a minimum distance
!           ! check for whether free surface lies inside
!           ! the fixed surface by checking the sign of the
!           ! dot product between the normal and the vector
!           ! from (u2,v2) to (u1,v1)
!           ! get r1 and z1
            call dmnf1 (u1, v1, r1, z1)
            phi = alp*v1
            xf = r1*cos(phi)
            yf = r1*sin(phi)
            zf = z1
!           ! compute vector from surface 2 to surface 1
            xvec = xf - xp
            yvec = yf - yp
            zvec = zf - zp
!           ! compute the normal at (u1,v1) on the fixed surface
            isurf = 1
            call surnorm (isurf, u1, v1, xun, yun, zun)
!           ! compute the dot product of the normal and the vector
            dotp = xvec*xun + yvec*yun + zvec*zun
            if (dotp < 0._dp) then
!              write(*,*) 'ip =',ip,'it = ',it
               d = -d
            else if (dotp == 0._dp) then
               distance = 0
               ierr = 1
               return
            endif
            distance = min(distance,d)
            ngood = ngood + 1
            dmat(ngood) = d

            grid_dist(ip,it) = d
         end do
      end do
!
c     compute the mean and the standard deviation
      dmean = sum(dmat(1:ngood)) / ngood
      std = sqrt(sum((dmat(1:ngood) - dmean)**2)/ngood)
      maxlocation = maxloc(abs(dmat(1:ngood)))
      dmax = dmat(maxlocation(1))

c     davg is the mean of the absolute values
      davg = sum(abs(dmat(1:ngood))) / ngood
      ngood2 = 0
      drms = 0
      do i = 1, ngood
         if(abs(dmat(i) - dmean) .le. 2*std) then
            ngood2 = ngood2 + 1
            drms = drms + dmat(i)**2
         endif
      enddo
      drms = sqrt(drms/ngood2)
!
      deallocate(dmat, xm1, xn1, rmn1, zmn1, xm2, xn2, rmn2, zmn2)

      end subroutine surfsep


      subroutine dmnf1(u, v, sumr, sumz) 
      use kind_spec
      use cmnf1
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(rprec) u, v, sumr, sumz 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(rprec) :: pi2 
!-----------------------------------------------
!
!     constants
      pi2 = 4*asin(1._dp) 
!
      sumr = sum(rmn1(:nmn1)*cos(pi2*(xm1(:nmn1)*u+xn1(:nmn1)*v))) 
      sumz = sum(zmn1(:nmn1)*sin(pi2*(xm1(:nmn1)*u+xn1(:nmn1)*v))) 
!
      end subroutine dmnf1
 
 
      subroutine dmnf2(u, v, sumr, sumz) 
      use kind_spec
      use cmnf2
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(rprec) u, v, sumr, sumz 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(rprec) :: pi2 
!-----------------------------------------------
!
!
!     constants
      pi2 = 4*asin(1._dp) 
!
      sumr = sum(rmn2(:nmn2)*cos(pi2*(xm2(:nmn2)*u+xn2(:nmn2)*v))) 
      sumz = sum(zmn2(:nmn2)*sin(pi2*(xm2(:nmn2)*u+xn2(:nmn2)*v))) 
 
      end subroutine dmnf2
 
 
      subroutine fdjac(n, x, fvec, np, df) 
!  (c) copr. 1986-92 numerical recipes software
 
      use kind_spec
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer n, np 
      real(rprec), dimension(n) :: x, fvec 
      real(rprec), dimension(np,np) :: df 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: nmax = 40 
      real(rprec), parameter :: eps = 1.e-4_dp 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j 
      real(rprec) :: h, temp 
      real(rprec), dimension(nmax) :: f 
!-----------------------------------------------
!
      do j = 1, n 
         temp = x(j) 
         h = eps*abs(temp) 
         if (h == 0._dp) h = eps 
         x(j) = temp + h 
         h = x(j) - temp 
         call funcv (n, x, f) 
         x(j) = temp 
         df(:n,j) = (f(:n)-fvec(:n))/h 
      end do 

      end subroutine fdjac
 
 
      function fmin_nr (x) 
!  (c) copr. 1986-92 numerical recipes software
 
      use kind_spec
      use newtv
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(rprec), dimension(*) :: x 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i 
      real(rprec) :: fmin_nr, r1
!-----------------------------------------------
!
!     uses funcv
!
      call funcv (nn, x, fvec) 
      r1 = dot_product(fvec(:nn),fvec(:nn)) 
      fmin_nr = r1/2 

      end function fmin_nr
 
 
      subroutine funcv(n, xx, fvec) 
 
      use kind_spec
      use cmnf1
      use geom
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer n 
      real(rprec), dimension(n) :: xx, fvec 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i 
      real(rprec) :: zv, cm, rv, u, v, 
     1   ru,  zu, r, z, cn, dy, dz, dx, yy, 
     2   coh, sih, xy, y, xv, yu, x, si, co, yv, xu 
!-----------------------------------------------
!
!
!
      u = xx(1) 
      v = xx(2) 
!
!     initializations
      r = 0; z = 0; ru = 0; zu = 0; rv = 0; zv = 0 
 
      do i = 1, nmn1 
         cm = xm1(i)*pi2 
         cn = xn1(i)*pi2 
         co = cos(cm*u + cn*v) 
         si = sin(cm*u + cn*v) 
         r = r + rmn1(i)*co 
         z = z + zmn1(i)*si 
         ru = ru - cm*rmn1(i)*si 
         rv = rv - cn*rmn1(i)*si 
         zu = zu + cm*zmn1(i)*co 
         zv = zv + cn*zmn1(i)*co 
      end do 
!
      coh = cos(alp*v) 
      sih = sin(alp*v) 
      x = coh*r 
      y = sih*r 
      xu = coh*ru 
      yu = sih*ru 
      xv = coh*rv - alp*y 
      yv = sih*rv + alp*x 
!
      dx = x - xp 
      dy = y - yp 
      dz = z - zp 
!
      fvec(1) = 2*(dx*xu + dy*yu + dz*zu) 
      fvec(2) = 2*(dx*xv + dy*yv + dz*zv) 
!
      end subroutine funcv
 
 
      subroutine lnsrch(n, xold, fold, g, p, x, f, stpmax, check, func) 
!  (c) copr. 1986-92 numerical recipes software
 
      use kind_spec
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer n 
      real(rprec) fold, f, stpmax, func 
      logical check 
      real(rprec), dimension(n) :: xold, g, p, x 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(rprec), parameter :: alf = 1.e-4_dp
      real(rprec), parameter :: tolx = 1.e-7_dp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i 
      real(rprec) :: a, alam, alam2, alamin, b, disc, f2, fold2, 
     1    rhs1, rhs2, slope, sumn, temp, test, tmplam 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      EXTERNAL func 
!-----------------------------------------------
!
!     uses func
!
      check = .FALSE. 
!
      sumn = dot_product(p,p) 
      sumn = sqrt(sumn) 
      if (sumn > stpmax)  p = p*stpmax/sumn 
      slope = dot_product(g,p) 
      test = maxval(abs(p)/max(abs(xold),1._dp))
      alamin = tolx/test 
      alam = 1 
    1 continue 
      x = xold + alam*p 
      f = func(x) 
      if (alam < alamin) then 
         x = xold 
         check = .TRUE. 
         return  
      else if (f <= fold + alf*alam*slope) then 
         return  
      else 
         if (alam == 1._dp) then 
            tmplam = -slope/(2*(f - fold - slope)) 
         else 
            rhs1 = f - fold - alam*slope 
            rhs2 = f2 - fold2 - alam2*slope 
            a = (rhs1/alam**2 - rhs2/alam2**2)/(alam - alam2) 
            b = ((-alam2*rhs1/alam**2) + alam*rhs2/alam2**2)
     1           /(alam - alam2) 
            if (a == 0._dp) then 
               tmplam = -slope/(2*b) 
            else 
               disc = b*b - 3*a*slope 
               tmplam = ((-b) + sqrt(disc))/(3*a) 
            endif 
            tmplam = min(.5_dp*alam,tmplam) 
         endif 
      endif 
      alam2 = alam 
      f2 = f 
      fold2 = fold 
      alam = max(tmplam,.1_dp*alam) 
      go to 1 

      end subroutine lnsrch
 
 
      subroutine lubksb(a, n, np, indx, b) 
!  (c) copr. 1986-92 numerical recipes software
 
      use kind_spec
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer n, np 
      integer, dimension(n) :: indx 
      real(rprec), dimension(np,np) :: a 
      real(rprec), dimension(n) :: b 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, ii, j, ll 
      real(rprec) :: sumn 
!-----------------------------------------------
!
!
      ii = 0 
      do i = 1, n 
         ll = indx(i) 
         sumn = b(ll) 
         b(ll) = b(i) 
         if (ii .ne. 0) then 
            sumn = sumn - sum(a(i,ii:i-1)*b(ii:i-1)) 
         else if (sumn .ne. 0._dp) then 
            ii = i 
         endif 
         b(i) = sumn 
      end do 
      do i = n, 1, -1 
         sumn = b(i) 
         sumn = sumn - sum(a(i,i+1:n)*b(i+1:n)) 
         b(i) = sumn/a(i,i) 
      end do 

      end subroutine lubksb


      subroutine ludcmp(a, n, np, indx, d) 
!  (c) copr. 1986-92 numerical recipes software
 
      use kind_spec
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer n, np 
      real(rprec) d 
      integer, dimension(n) :: indx 
      real(rprec), dimension(np,np) :: a 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: nmax = 500 
      real(rprec), parameter :: tiny = 1.0e-20_dp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, imax, j, k 
      real(rprec) :: aamax, dum, sumn 
      real(rprec), dimension(nmax) :: vv 
!-----------------------------------------------
!
      d = 1 
      vv = 0
      do i = 1, n 
         aamax = maxval(abs(a(i,:n)))
         if (aamax == 0._dp) then 
            write (*, '(2a)') 'PAUSE ', 'singular matrix in ludcmp' 
         else 
            vv(i) = 1._dp/aamax 
         endif 
      end do 
      do j = 1, n 
         do i = 1, j - 1 
            sumn = a(i,j) 
            sumn = sumn - sum(a(i,:i-1)*a(:i-1,j)) 
            a(i,j) = sumn 
         end do 
         aamax = 0 
         do i = j, n 
            sumn = a(i,j) 
            sumn = sumn - sum(a(i,:j-1)*a(:j-1,j)) 
            a(i,j) = sumn 
            dum = vv(i)*abs(sumn) 
            if (dum >= aamax) then 
               imax = i 
               aamax = dum 
            endif 
         end do 
         if (j /= imax) then 
            do k = 1, n 
               dum = a(imax,k) 
               a(imax,k) = a(j,k) 
               a(j,k) = dum 
            end do 
            d = -d 
            vv(imax) = vv(j) 
         endif 
         indx(j) = imax 
         if (a(j,j) == 0._dp) a(j,j) = tiny 
         if (j /= n) then 
            dum = 1._dp/a(j,j) 
            a(j+1:n,j) = a(j+1:n,j)*dum 
         endif 
      end do 

      end subroutine ludcmp

 
      subroutine newt(x, n, check) 
!  (c) copr. 1986-92 numerical recipes software
 
!     uses fdjac,fmin_nr,lnsrch,lubksb,ludcmp
 
      use kind_spec
      use newtv
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer :: n 
      logical :: check 
      real(rprec), dimension(n) :: x 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: maxits = 250 
      real(rprec), parameter :: tolf = 1.e-4_dp 
      real(rprec), parameter :: tolmin = 1.e-6_dp
      real(rprec), parameter :: tolx = 1.e-7_dp 
      real(rprec), parameter :: stpmx = 100 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, its, j 
      integer, dimension(np) :: indx 
      real(rprec) :: d, den, f, fold, stpmax, temp, test, sumn 
      real(rprec), dimension(np,np) :: fjac 
      real(rprec), dimension(np) :: g, p, xold 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      real(rprec) , EXTERNAL :: fmin_nr 
!-----------------------------------------------
!
      if (n > np) stop 'N > NP in NEWT' 
      nn = n 
      f = fmin_nr(x) 
      test = 0 
      do i = 1, n 
         if (abs(fvec(i)) > test) test = abs(fvec(i)) 
      end do 
      if (test < tolf/100) return  
!
      sumn = dot_product(x(:n),x(:n)) 
!
      stpmax = stpmx*max(sqrt(sumn),real(n,rprec)) 
!
      do its = 1, maxits 
         call fdjac (n, x, fvec, np, fjac) 
         do i = 1, n 
            sumn = sum(fjac(:n,i)*fvec(:n)) 
            g(i) = sumn 
         end do 
         fold = f 
         xold(:n) = x(:n) 
         p(:n) = -fvec(:n) 
         call ludcmp (fjac, n, np, indx, d) 
         call lubksb (fjac, n, np, indx, p) 
         call lnsrch (n, xold, fold, g, p, x, f, stpmax, check, fmin_nr) 
         test = maxval(abs(fvec(:n)))
         if (test < tolf) then 
            check = .FALSE. 
            return  
         endif 
         if (check) then 
            den = max(f,.5_dp*n) 
            test = maxval(abs(g(:n))*max(abs(x(:n)),1._dp)/den)
            if (test < tolmin) then 
               check = .TRUE. 
            else 
               check = .FALSE. 
            endif 
            return  
         endif 
         test = maxval( abs(x(:n)-xold(:n))/max(abs(x(:n)),1._dp) )
         if (test < tolx) return  
      end do 
!      write(*,*) 'maxits exceeded in newt'
      check = .TRUE. 

      end subroutine newt
 
 
      subroutine p2surf(u, v, d, ierr)
      use kind_spec
      use geom
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer ierr
      real(dp) u, v, d
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: n = 2
      integer, parameter :: ldfjac = 2
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: tries
      real(dp), dimension(n) :: xx
      real(dp) :: step2, r, step1, x, y, z, phi, factor
      logical :: check
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      EXTERNAL newt
!-----------------------------------------------
!
      ierr = 0
      call random_seed                           ! initialize the random number generator
      factor = 100
      xx(1) = u
      xx(2) = v
!
      check = .FALSE.                            ! presume success to start
      do tries = 1, 20
!
         call newt (xx, n, check)
!
         if (.not.check) then
            exit
         else
!            write(*,*) 'check = ', check
            call random_number (step1)
            call random_number (step2)
            if (step1 <= 0.5_dp) step1 = -step1
            if (step2 <= 0.5_dp) step2 = -step2
            xx(1) = (1 + step1)*xx(1)
            xx(2) = (1 + step2)*xx(2)
         endif
      end do
!
!
!     if successful get point r,z on plasfree at u,v
      if (.not.check) then
         call dmnf1 (xx(1), xx(2), r, z)

         u = xx(1)
         v = xx(2)

         phi = alp*v
         x = r*cos(phi)
         y = r*sin(phi)
         d = sqrt((x - xp)**2 + (y - yp)**2 + (z - zp)**2)
         ierr = 0
      else
         ierr = 1
      endif

      end subroutine p2surf

 
      subroutine surnorm(isurf, u, v, xnorm, ynorm, znorm) 
      use kind_spec
      use cmnf1
      use cmnf2
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer isurf 
      real(rprec) u, v, xnorm, ynorm, znorm 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i 
      real(rprec) :: si, coh, sih, co, cm, cn, xv, yv, crmag, yu, x, y,
     1     xu, pi2, zu, rv, zv, ru, alp, r, z 
!-----------------------------------------------
!
!     compute the normal to the surface at (u,v)
!
!
      pi2 = 4*asin(1._dp) 
      if (isurf == 1) then 
         alp = pi2/q1 
      else if (isurf == 2) then 
         alp = pi2/q1 
      else 
         stop 'ISURF must be surface 1 or 2' 
      endif 
!
      r = 0;  z = 0;  ru = 0;  zu = 0;  rv = 0;  zv = 0 
!
      if (isurf == 1) then 
         do i = 1, nmn1 
            cm = xm1(i)*pi2 
            cn = xn1(i)*pi2 
            co = cos(cm*u + cn*v) 
            si = sin(cm*u + cn*v) 
            r = r + rmn1(i)*co 
            z = z + zmn1(i)*si 
            ru = ru - cm*rmn1(i)*si 
            rv = rv - cn*rmn1(i)*si 
            zu = zu + cm*zmn1(i)*co 
            zv = zv + cn*zmn1(i)*co 
         end do 
!
         coh = cos(alp*v) 
         sih = sin(alp*v) 
!
         x = coh*r 
         y = sih*r 
         xu = coh*ru 
         yu = sih*ru 
         xv = coh*rv - alp*y 
         yv = sih*rv + alp*x 
!
      else if (isurf == 2) then 
!
         do i = 1, nmn2 
            cm = xm2(i)*pi2 
            cn = xn2(i)*pi2 
            co = cos(cm*u + cn*v) 
            si = sin(cm*u + cn*v) 
            r = r + rmn2(i)*co 
            z = z + zmn2(i)*si 
            ru = ru - cm*rmn2(i)*si 
            rv = rv - cn*rmn2(i)*si 
            zu = zu + cm*zmn2(i)*co 
            zv = zv + cn*zmn2(i)*co 
         end do 
!
         coh = cos(alp*v) 
         sih = sin(alp*v) 
!
         x = coh*r 
         y = sih*r 
         xu = coh*ru 
         yu = sih*ru 
         xv = coh*rv - alp*y 
         yv = sih*rv + alp*x 
!
      else 
!
         stop 'ISURF must be surface 1 or 2' 
!
      endif 
!
!     compute the components of (Xu x Xv)
      xnorm = yu*zv - yv*zu 
      ynorm = -(xu*zv - xv*zu) 
      znorm = xu*yv - xv*yu 
!     ! magnitude of (Xu x Xv)
      crmag = sqrt(xnorm*xnorm + ynorm*ynorm + znorm*znorm) 
!     ! normalize the components and multiply by -1
!     ! to get right sign of the normal
      xnorm = -xnorm/crmag 
      ynorm = -ynorm/crmag 
      znorm = -znorm/crmag 
 
      end subroutine surnorm


      subroutine weighted_surfsep(rmnc,zmns,xm,xn,nfp,mnmax,rms,unrms,
     1         dist_max)
C-----------------------------------------------
C EAL 9/00
C The boundary is a function of theta, the target a function of angle
C d=sqrt( (rb-rt)**2 + (zb-zt)**2)
C Dd= d (d)/d(angle)
C search angle ~ theta for a zero crossing of Dd, find nearest theta values
C symmetric about theta0
C find the root Dd=0 via quadratic interpolation
C if quadratic formula fails use linear interpolation
C on exit angle contains nu values which make the nu points (Rt,Zt) nearest
C to the nu points (Rb,Zb) evaluated at theta. (For each point)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use kind_spec
      use optim, only: rbc_vv, zbs_vv, mpol_vv, ntor_vv, 
     1  theta0=>theta0_bw, phi0=>phi0_bw, wtheta=>wtheta_bw, 
     2  wphi=>wphi_bw, amplw=>amplw_bw, planes=>planes_bw
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(*), intent(in) :: 
     1   rmnc, zmns, xm, xn
      real(rprec),  intent(out) :: rms, unrms, dist_max
      integer, intent(in) :: mnmax, nfp
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      real(rprec) , external :: weightfcn
C-----------------------------------------------


C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: mu = 12     !       modes - poloidal
      integer, parameter :: nu = 1024   !       points - poloidal
      integer, parameter :: nv = 16     !       slices - toroidal
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: niter, nstep, nfp_d, nphi_d, nphi2, itype,
     1   ntheta_d, ii, i, k, j, m, n, mexp, mpol_d, mpol_d1,
     2  k1, k2, ijk
      real(rprec) :: twopi, zeta, sumwgt, wgtampl(3),
     1   err, aiota, dumphi, arg, tol, rnow, znow, pit, v, u, xxm, xxn,
     2   cosa, sina, a, b, c, fa, fb, fc, r, s, t, p, q, sol, 
     3  r1, r2, z1, z2, rp2, zp2, sp, val, rpnow, zpnow
      real(rprec), dimension(:), allocatable :: dummy1, dummy2
      integer, dimension(:), allocatable :: index, mask
      integer, dimension(1) :: iloc, jloc
      real(rprec), dimension(:), allocatable :: rin, zin,
     1  trin, tzin, tprin,tpzin, theta, angle, wgt
      real(rprec) :: rmin, rmax, zmin, zmax, siz
      integer :: kz, nplanes1
      character*100 lt 
      logical first
      data first/.true./
      data nplanes1/4/ ! 1 + # planes in weight function
C-----------------------------------------------
      sp(r1,r2,z1,z2,rp2,zp2)=1/sqrt((r1-r2)*(r1-r2)+(z1-z2)*(z1-z2))*
     .(rp2*(r2-r1)+zp2*(z2-z1))
C-----------------------------------------------
C-----------------------------------------------
      twopi = 8*atan(1._dp)
      dist_max = 0
      wgtampl=amplw
      ntheta_d=nu;nphi_d=nv
      mpol_d = mu
      mpol_d1 = mpol_d - 1
      nphi2 = 1 + nphi_d/2
      if(.not.allocated(rin))
     .allocate(rin(ntheta_d),zin(ntheta_d))
      if(.not.allocated(trin))
     .allocate(trin(ntheta_d),tzin(ntheta_d))
      if(.not.allocated(tprin))
     .allocate(tprin(ntheta_d),tpzin(ntheta_d))
      if(.not.allocated(theta))allocate(wgt(ntheta_d),
     .theta(ntheta_d),angle(ntheta_d),index(ntheta_d),mask(ntheta_d))
      if(.not.allocated(dummy1))allocate(dummy1(ntheta_d))
      if(.not.allocated(dummy2))allocate(dummy2(ntheta_d))
 
10001 continue
CEAL      v=twopi*(nfp-1)*(kz-1)/nplanes1/2
      rms=0
      unrms=0
      sumwgt=0
      do ijk=1,nplanes1-1
      v=planes(ijk)/nfp
      trin=0;tzin=0;rin=0;zin=0;tprin=0;tpzin=0
         call totbdy(ntheta_d,v,
     .          rin(1:ntheta_d),
     .          zin(1:ntheta_d),
     .          rmnc,zmns,xm,xn,mnmax)
       do j=1,ntheta_d
        u=twopi*(j-1)/(ntheta_d-1)
        theta(j)=u
        wgt(j)= weightfcn(theta(j), v, theta0(ijk), phi0,
     .   wtheta, wphi, amplw, Nfp)
        rnow=0;znow=0;rpnow=0;zpnow=0
        do n=-ntor_vv,ntor_vv
         do m=0,mpol_vv
          xxm=m;xxn=n
          arg=xxm*u-xxn*planes(ijk)
          cosa=cos(arg);sina=sin(arg)
          rnow=rnow+rbc_vv(n,m)*cosa
          znow=znow+zbs_vv(n,m)*sina
          rpnow=rpnow-rbc_vv(n,m)*xxm*sina
          zpnow=zpnow+zbs_vv(n,m)*xxm*cosa
         enddo
        enddo
        trin(j)=rnow
        tzin(j)=znow
        tprin(j)=rpnow
        tpzin(j)=zpnow
       enddo                                 
 ! grid search
      index=(/(k,k=1,ntheta_d)/)
      i=0
      angle=0
      do while (i < size(theta))
c        k=nint(10.*ntheta_d/360.)
c        k=nint(20.*ntheta_d/360.)      !       seems better for li383
        k=nint(30.*ntheta_d/360.)       !       still better
        i=i+1
        mask=cshift(index,-k)
        k1=mask(i)
        mask=cshift(index,k)
        k2=mask(i)
        fa=sp(rin(i),trin(k2),zin(i),tzin(k2),tprin(k2),tpzin(k2))
        fb=sp(rin(i),trin(k1),zin(i),tzin(k1),tprin(k1),tpzin(k1))
        if(fb.eq.0.)angle(i)=theta(k1)
        if(fb.eq.0.)cycle
        val=fa/fb
        do while(val > 0 .and. k>0)
         k=k-1
         mask=cshift(index,-k)
         k1=mask(i)
         mask=cshift(index,k)
         k2=mask(i)
         fa=sp(rin(i),trin(k2),zin(i),tzin(k2),tprin(k2),tpzin(k2))
         fb=sp(rin(i),trin(k1),zin(i),tzin(k1),tprin(k1),tpzin(k1))
         val=fa/fb
        enddo   !       while(val < 0 .and. k>0)
          if(val.gt.0)then
           angle(j)=theta(j)
           cycle
          endif
        do while(val < 0 .and. k>0)
         k=k-1
         mask=cshift(index,-k)
         k1=mask(i)
         mask=cshift(index,k)
         k2=mask(i)
         fa=sp(rin(i),trin(k2),zin(i),tzin(k2),tprin(k2),tpzin(k2))
         fb=sp(rin(i),trin(k1),zin(i),tzin(k1),tprin(k1),tpzin(k1))
         if(fa.eq.0.)angle(i)=theta(k2)
         if(fb.eq.0.)angle(i)=theta(k1)
         if(fb.eq.0.)cycle            
         if(fa.eq.0.)cycle            
         val=fa/fb
        enddo   !       while(val < 0 .and. k>0)
        if(fb.eq.0.)cycle
        k=k+1
        mask=cshift(index,-k)
        k1=mask(i)
        mask=cshift(index,k)
        k2=mask(i)
!       find zero crossing
        a=theta(k1)
        b=theta(i)
        c=theta(k2)
        if(b < a) a=a-twopi
        if(c < b) c=c+twopi
        if(c < a) c=c+twopi

        fa=sp(rin(i),trin(k1),zin(i),tzin(k1),tprin(k1),tpzin(k1))
        fb=sp(rin(i),trin(i),zin(i),tzin(i),tprin(i),tpzin(i))
        fc=sp(rin(i),trin(k2),zin(i),tzin(k2),tprin(k2),tpzin(k2))
        if(fa.eq.0.)angle(i)=a
        if(fb.eq.0.)angle(i)=b
        if(fc.eq.0.)angle(i)=c
        if((fa.eq.0.).or.(fb.eq.0.).or.(fc.eq.0.))cycle

10002   if(fc/fb < 0)then 
            sol=-(c*fb-b*fc)/(fc-fb)
        elseif(fa/fb < 0)then
            sol=(a*fb-b*fa)/(fb-fa)
        else
            sol=theta(i)
        endif   !       if(fc/fb < 0)
10003   angle(i)=sol
      enddo     !       while(val < 0 .and. k>0)
      where(angle < 0.)angle=angle+twopi
      where(angle > twopi)angle=angle-twopi

! evaluate trin,tzin on new mesh
       do j=1,ntheta_d
        u=angle(j)
        rnow=0;znow=0
        do n=-ntor_vv,ntor_vv
         do m=0,mpol_vv
          xxm=m;xxn=n
          arg=xxm*u-xxn*planes(ijk)
          cosa=cos(arg);sina=sin(arg)
          rnow=rnow+rbc_vv(n,m)*cosa
          znow=znow+zbs_vv(n,m)*sina
         enddo
        enddo
        trin(j)=rnow
        tzin(j)=znow
       enddo
        dummy1(1:ntheta_d)=wgt(1:ntheta_d)*sqrt(
     &  (rin(1:ntheta_d)-trin(1:ntheta_d))**2 +
     &  (zin(1:ntheta_d)-tzin(1:ntheta_d))**2
     &)
        dummy2(1:ntheta_d)=sqrt(
     &  (rin(1:ntheta_d)-trin(1:ntheta_d))**2 +
     &  (zin(1:ntheta_d)-tzin(1:ntheta_d))**2
     &)
!      rms=rms+sum(dummy1(1:ntheta_d))/sum(wgt)
!       need normalization to get lengths, but do not remove planar
!        weighting
      rms=rms+sum(dummy1(1:ntheta_d))/sum(wgt)*wgtampl(ijk)/sum(wgtampl)

      unrms=unrms+sum(dummy2(1:ntheta_d))/ntheta_d
      sumwgt=sumwgt+sum(wgt)
      enddo     !       ijk
      ijk=ijk-1
      unrms=unrms/ijk
      rmin = minval(rin(1:ntheta_d))
      rmax = maxval(rin(1:ntheta_d))
      rmin=min(rmin,minval(trin(1:ntheta_d)))
      rmax=max(rmax,maxval(trin(1:ntheta_d)))
      zmin = minval(zin(1:ntheta_d))
      zmax = maxval(zin(1:ntheta_d))
      zmin=min(zmin,minval(tzin(1:ntheta_d)))
      zmax=max(zmax,maxval(tzin(1:ntheta_d)))
! deallocate
       if(allocated(dummy1))deallocate(dummy1)
       if(allocated(dummy2))deallocate(dummy2)
       if(allocated(tzin))deallocate(tzin)
       if(allocated(trin))deallocate(trin)
       if(allocated(tzin))deallocate(tpzin)
       if(allocated(trin))deallocate(tprin)
       if(allocated(rin))deallocate(rin)
       if(allocated(zin))deallocate(zin)
       if(allocated(zin))deallocate(theta)
       if(allocated(zin))deallocate(angle)
       if(allocated(zin))deallocate(index)
      return
 2100 format(' RBC(',i3,',',i2,') = ',1pe12.4,3x,
     1       ' ZBS(',i3,',',i2,') = ',1pe12.4)               
      end subroutine weighted_surfsep
      subroutine totbdy(ntheta, phiin, r, z, 
     1   rmnc, zmns, xm, xn, mnmax )
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ntheta, nplanes1, mnmax
      real(rprec), dimension(*), intent(out) :: r, z
      real(rprec), dimension(*), intent(in) :: 
     1   rmnc, zmns, xm, xn 
      real(rprec) :: phiin

C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nrth, l, j, mn, kt, m, iarg
      real(rprec) :: pit, piz, xm0, xn0, arg, twopi
C-----------------------------------------------
      twopi = atan(1.)*8. 
      pit = 1.d0/(ntheta - 1)
      nrth = 1*ntheta
      r(:nrth) = 0.d0
      z(:nrth) = 0.d0
         do mn = 1, mnmax
            m = mn 
            if(rmnc(m).eq.0.)cycle
            xm0 = xm(mn)*pit
            xn0 = xn(mn)
            do kt = 1, ntheta
               arg = xm0*(kt - 1) - xn0*phiin/twopi
               iarg = arg
               arg = twopi*(arg - iarg)
               r(kt) = r(kt) + rmnc(m)*cos(arg) 
               z(kt) = z(kt) + zmns(m)*sin(arg) 
            end do
         end do
      end subroutine totbdy


      function weightfcn(thetaw, phiw, theta0, phi0, wtheta, wphi,
     .   amplw, Nfp)
C EAL 9/00      Form sum of gaussian weights, in theta, phi plane
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: nfp
      real(rprec) :: thetaw, phiw, theta0, phi0, 
     .wtheta, wphi, amplw(3)
C-----------------------------------------------
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec) :: e, pi, weightfcn
      real(rprec) :: c1, c2, aweight, a(3)
      integer :: i, j
C-----------------------------------------------
      e=0
      Pi=4*atan(1._dp)
      do j = 1,3
       c2=phi0 + float(j-1)/2.*Pi/nfp
       do i = 1,5
        c1=theta0 + float(i-1)/4.*(2.*pi)
        e=e+Amplw(j)*aweight(thetaw,phiw,c1,c2,wtheta,wphi)
       enddo
      enddo
      weightfcn = e
      return
      end function weightfcn


      function aweight(x, y, a, b, v, w)
CEAL    9/00    !       avoid underflow in evaluating a 2D gaussian
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: nfp
      real(rprec) ::thetaw, phiw, theta0, phi0, wtheta, wphi, amplw(3)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec) :: e, pi, aweight, x, y
      real(rprec) :: n12, n3, n4, a, b, w, v
      data n3/-1.e10_dp/
C-----------------------------------------------
      Pi=4*atan(1._dp)
      n12=(-((-a + x)**2/v**2) - (-b + y)**2/w**2)
      if(n12 > n3)then
         n4=dexp(n12)
      else
         aweight=0
         return
      endif
      aweight = n4

      end function aweight


! *************************************************************************
!  LIMSEP
!
!  This subroutine computes minimum separation distance (signed) between
!  the plasma and a set of limiting surfaces, described as an array of
!  piecewise linear curves at specified toroidal angles.
!
!  This subroutine is modified from SURFSEP, which does a similar calculation
!  to a bounding surface specified by Fourier Coefficients.
!
!  M. Zarnstorff   July 2002
! *************************************************************************

      subroutine limsep(mnmax_pl, r_pl, z_pl, xm_pl, xn_pl, nphi_lim,
     1    phi_lim, nrz_max, r_lim, z_lim, nper, distance, ierr)
      use kind_spec
!      use cmnf1
      use cmnf2
!      use geom

      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer mnmax_pl, nphi_lim, nrz_max, nper, ierr
      real(rprec) distance
      real(rprec), dimension(mnmax_pl) :: r_pl, z_pl, xm_pl, xn_pl
      real(rprec), dimension(nphi_lim) :: phi_lim
      real(rprec), dimension(nrz_max, nphi_lim) :: r_lim, z_lim
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer, parameter :: igrid = 64

      integer :: i, ip, it, il, iu, np
      real(rprec) :: d, du, u2, v2, u2_new, d_min
      real(rprec), dimension(igrid) :: r2, z2, dist

      integer :: nlim, it_min, ip_min, il_min
      integer, dimension(1) :: iresult, iresult2

!-----------------------------------------------
!
      ierr = 0                                   ! presume success to start

!
!     set np to igrid
      np = igrid

!
!     compute the total number of modes for each surface
      nmn2 = mnmax_pl
!
!     equate the common block variables with the arguments
!
      allocate (xm2(nmn2), xn2(nmn2), rmn2(nmn2), zmn2(nmn2),
     1          stat=ierr)
      if (ierr .ne. 0) return
!
!
!     equate the common block variables with the arguments
!     Surface 2 - Plasma
!
      xm2(:nmn2) = xm_pl(:nmn2)
      xn2(:nmn2) = xn_pl(:nmn2)
      rmn2(:nmn2) = r_pl(:nmn2)
      zmn2(:nmn2) = z_pl(:nmn2)
!
!     write (6,*) 'done copying data', nmn2,nmn1
!
      distance = 1.0e30_dp

      du = 1.0_dp/np

      d_min = 1.0e30_dp

      it_min = -1
      ip_min = -1
      il_min = -1

!
!  grided search to find nearest segment and plasma location
!
      do it = 1, nphi_lim
         if( phi_lim(it) > 360 .or. phi_lim(it) < -360 ) exit
         if( r_lim(1,it) == 0 .or. r_lim(2,it) == 0) exit

         do ip = 3, nrz_max
            if( r_lim(ip, it) == 0 ) then
                nlim = ip - 1
                exit
            endif
         enddo

         v2 = modulo((nper*phi_lim(it))/360, 1._rprec)

!
!   get (r,z) points on plasma surface
! 
         do ip = 1, np
            u2 = (ip - 1)*du
            call dmnf2 (u2, v2, r2(ip), z2(ip))
         enddo

         do il = 2, nlim
            do ip = 1, np
               call lim_dist( r2(ip), z2(ip), r_lim(il-1, it), 
     1                        z_lim(il-1, it), dist(ip))
            enddo
            
            iresult = minloc(dist)

            if( minval(dist) < d_min) then
               d_min = minval(dist)
               it_min = it
               iresult = minloc(dist)
               ip_min = iresult(1)
               il_min = il
            endif
         enddo
      enddo

!
!  if no valid limiter segments specified, we are done
!
      if( it == -1) then
         distance = 0
         ierr = 1
         return
      endif

!
!     refine the point of minimum distance from the plasma
!     to the limiter segment by binary search
!
!      print *,'min found:',phi_lim(it_min),d_min,il_min
      v2 = modulo((nper*phi_lim(it_min))/360, 1._rprec)

      u2 = (ip_min - 1)*du
      u2_new = u2

      do ip = 1, 10
         du = du/2

         do iu = -1, 1, 2
            call dmnf2( u2+du*iu, v2, r2(1), z2(1))

            call lim_dist( r2(1), z2(1), r_lim(il_min-1, it_min), 
     1                        z_lim(il_min-1, it_min), dist(1))
         
            if( dist(1) < d_min ) then
               u2_new = u2 + du*iu
               d_min = dist(1)
            endif
         enddo

         u2 = u2_new
      enddo

      distance = d_min

      deallocate(xm2, xn2, rmn2, zmn2)

      end subroutine limsep


      subroutine lim_dist( rp, zp, r_lim, z_lim, d)
      use kind_spec
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(rprec) :: rp, zp, d
      real(rprec), dimension(2) :: r_lim, z_lim
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(rprec) :: dlr, dlz, dl, dr, dz, d1, d2

!-----------------------------------------------

      dlr = r_lim(2) - r_lim(1)
      dlz = z_lim(2) - z_lim(1)
      dl = (dlr**2 + dlz**2)

      dr = (r_lim(1)*dlz**2 + rp*dlr**2 - dlr*dlz*(z_lim(1)-zp)) 
     1     / dl

      dz = z_lim(1) + (dr - r_lim(1)) * dlz / dlr

      d = (rp-dr)**2 + (zp-dz)**2
      d1 = (dr-r_lim(1))**2 + (dz-z_lim(1))**2
      d2 = (dr-r_lim(2))**2 + (dz-z_lim(2))**2

      if( d1 > dl .or. d2 > dl) then
         d = 1.e30
      else
         d = sign(sqrt(d), dlr*(zp-dz) - dlz*(rp-dr))
      endif

      end subroutine lim_dist

EOF
EOC




