

      module date_and_computer
      character*(3), dimension(12), parameter :: months =
     1  ( / 'Jan','Feb','Mar','Apr','May','Jun',
     2      'Jul','Aug','Sep','Oct','Nov','Dec' / )
      character*(*), parameter :: computer =
#if defined(RISC)
     1  ' IBM RISC-6000 WORKSTATION. '
#elif defined(IRIX64)
     1  ' ORIGIN 2000. '
#elif defined(IRIX)
     1  ' PPPL CCF. '
#elif defined(HPUX)
     1  ' HP-UX WORKSTATION. '
#elif defined(OSF1)
     1  ' DECSTATION AXP. '
#elif defined(SUN)
     1  ' SUN WORKSTATION. '
#elif defined(VAX)
     1  ' DEC VAX. '
#elif defined(WINNT)
     1  ' WINDOWS NT SYSTEM. '
#elif defined(LINUX)
     1  ' LINUX SYSTEM. '
#elif defined(AXPVMS)
     1  ' AXP/VMS. '
#elif defined(CRAY)
     1  ' CRAY SUPER-COMPUTER. '
#elif defined(SX5)
     1  ' SX-5 SUPER-COMPUTER. '
#else
     1  ' '
#endif
      end module date_and_computer

      module vparams
      use kind_spec
      implicit none
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!
!     MAXIMUM PARAMETERS FOR VMEC CODE (FOR READING INPUT)
!     USER SHOULD NOT ALTER THESE
!
      integer, parameter :: nsd = 3001      !maximum number of radial nodes
      integer, parameter :: mpold = 61     !maximum number of poloidal harmonics (in r,z,lam fourier series)
      integer, parameter :: ntord = 61     !maximum number of toroidal harmonics
      integer, parameter :: ndatafmax = 101
      integer, parameter :: nstore_seq = 100

!
!     CONSTANTS
!
      integer, parameter :: nthreed0=9, nmac0=nthreed0+1,
     1                      indata0=nthreed0+2, nwout0 = nthreed0+3,
     2                      jxbout0 = nthreed0+4
      integer, parameter :: nfort8 = 28, nfort18 = 18,  nbfld0 = 50
      integer, parameter :: nfort21 = 21
      integer, parameter :: nlog0 = 51, nmercier0 = 52
      integer :: nthreed, nmac, nlog

!
!     DERIVED (FROM FUNDAMENTAL) PARAMETERS FOR VMEC CODE
!
      integer, parameter :: mpol1d = mpold - 1
      integer, parameter :: ntor1d = 1 + ntord

!
!     MISCELLANEOUS PARAMETERS
!
      real(rprec), parameter :: c1pm2  =  1.e-2_dp
      real(rprec), parameter :: cp15   =  0.15_dp
      real(rprec), parameter :: cp25   = 0.25_dp
      real(rprec), parameter :: cp5    = 0.50_dp
      real(rprec), parameter :: c1pm8  = 1.0e-8_dp
      real(rprec), parameter :: one    = 1.0_dp
      real(rprec), parameter :: cbig   = 0.9e30_dp
      real(rprec), parameter :: c2p0   = 2.0_dp
      real(rprec), parameter :: c3p0   = 3.0_dp
      real(rprec), parameter :: zero   = 0.0_dp
      real(rprec), parameter :: cp05   = 0.05_dp
      real(rprec), parameter :: c1pm13 = 1.0e-13_dp
      real(rprec), parameter :: twopi  = 6.28318530717958623_dp
      real(rprec), parameter :: dmu0   = 2.0e-7_dp*twopi
      real(rprec), parameter :: osqrt2 = 0.707106781186547462_dp
      real(rprec), parameter :: epstan = 1.e-10_dp

      end module vparams

      module vsvd0
      use kind_spec
      implicit none
      real(rprec), parameter :: fturnon_axis = 3.e-9_dp
      real(rprec), parameter :: fopt_axis = 3.e-2_dp*fturnon_axis
      integer, parameter :: isamecoil = -2
      integer, parameter :: needit = -1
      integer, parameter :: idontneed = 0
      integer, parameter :: isymcoil = 1
      integer, parameter :: ithom0 = 1
      integer, parameter :: istark0 = 2
      integer, parameter :: islope0 = 3
      integer, parameter :: icurr0 = 4
      integer, parameter :: idiam0 = 5
      integer, parameter :: iflxs0 = 6
      integer, parameter :: ibrzfld = 7
      integer, parameter :: natur = 0
      integer, parameter :: ideriv = 1
      integer, parameter :: intder = 1
      integer, parameter :: intfun = 2
      integer, parameter :: nmse = 100        !number of mse measurements
      integer, parameter :: ntse = 100        !number of thompson scattering measurements
      integer, parameter :: nfloops = 100     !number of external poloidal flux loops
      integer, parameter :: nbsetsp = 5       !number of external b-field loop sets allowed
      integer, parameter :: nbcoilsp = 100    !number of external b-field coils per set
      integer, parameter :: nbctotp = nbsetsp*nbcoilsp
      integer, parameter :: jngrn = 1001      !number of "greens function" points
      integer, parameter :: jchix = 7         !number of data types contributing to rms error match
      integer, parameter :: mstp = 100        !number of time steps to store chisq error
      integer, parameter :: jchix1 = jchix + 1
      integer, parameter :: nlimset = 2       !number of different limiters
      integer, parameter :: nparts = 4        !number of items needed to specify pf coils
      integer, parameter :: npfcoil = 40      !number of filaments in pf coil pack (for plotting)
      integer, parameter :: nigroup = 100     !number of external current groups
      integer, parameter :: ipedsvd = 8
      end module vsvd0

      module vmec_input
      use vparams, only: rprec, dp, mpol1d, ntord, ndatafmax
      use vsvd0
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nfp, ncurr, nsin, niter, nstep, nvacskip, mpol, ntor,
     1           ntheta, nzeta
      integer, dimension(100) :: ns_array
      integer :: imse, isnodes, itse, ipnodes, iopt_raxis,
     1   imatch_phiedge, nflxs, ac_form
      integer, dimension(nbsetsp) :: nbfld
      integer, dimension(nfloops) :: indxflx
      integer, dimension(nbcoilsp,nbsetsp) :: indxbfld
      real(rprec), dimension(-ntord:ntord,0:mpol1d) ::
     1   rbs, zbc, rbc, zbs
      real(rprec) :: time_slice, curtor, delt, ftol, tcon0,
     1   gamma, phiedge, phidiam, sigma_current, sigma_delphid, tensi,
     2   tensp, tensi2, fpolyi, presfac, mseangle_offset, pres_offset,
     3   mseangle_offsetm, spres_ped, bloat
      real(rprec), dimension(0:10) :: am, ai, ac, aphi
      real(rprec), dimension(0:ntord,2) :: raxis, zaxis
      real(rprec), dimension(100) :: ftol_array
      real(rprec), dimension(nigroup) :: extcur
      real(rprec), dimension(nmse) :: mseprof
      real(rprec), dimension(ntse) :: rthom, datathom, sigma_thom
      real(rprec), dimension(nmse) :: rstark, datastark,
     1    sigma_stark
      real(rprec), dimension(nfloops) :: dsiobt, sigma_flux
      real(rprec), dimension(nbcoilsp,nbsetsp) :: bbc, sigma_b
      real(rprec), dimension(ndatafmax) :: psa, pfa, isa, ifa
      logical :: lpofr, lmac, lfreeb, lrecon, loldout, ledge_dump,
     1           lasym, ldiagno
     2  , lspectrum_dump, loptim           !!Obsolete
      character*(100) :: mgrid_file
      character*(120) :: arg1
      character*(100) :: input_extension

      namelist /indata/ mgrid_file, time_slice, nfp, ncurr, nsin,
     1   lasym, niter, nstep, nvacskip, delt, ftol, gamma, am, ai, ac,
     2   aphi, rbc, zbs, rbs, zbc, raxis, zaxis, spres_ped,
     3   mpol, ntor, ntheta, nzeta, psa, pfa, isa, ifa,
     4   ns_array, ftol_array, tcon0, curtor, sigma_current, extcur,
     5   phiedge, imatch_phiedge, iopt_raxis, tensi, tensp,
     6   mseangle_offset, mseangle_offsetm, imse, isnodes, rstark,
     7   datastark, sigma_stark, itse, ipnodes, presfac, pres_offset,
     8   rthom, datathom, sigma_thom, phidiam, sigma_delphid, tensi2,
     9   fpolyi, nflxs, indxflx, dsiobt, sigma_flux, nbfld, indxbfld,
     A   bbc, sigma_b, lpofr, lfreeb, lrecon, lmac, loldout, ldiagno,
     B   ledge_dump, lspectrum_dump, loptim, bloat, ac_form

      namelist /mseprofile/ mseprof

      contains

      subroutine read_indata_namelist (iunit, istat)
      integer :: iunit, istat

      read (iunit, nml=indata, iostat=istat)

      end subroutine read_indata_namelist

      subroutine read_mse_namelist (iunit, istat)
      integer :: iunit, istat

      read (iunit, nml=mseprofile, iostat=istat)

      end subroutine read_mse_namelist

      end module vmec_input

      module vmec_seq
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: nseqmax = 100
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nseq
      integer, dimension(nseqmax) :: nseq_select
      character*120, dimension(nseqmax) :: extension

      namelist /vseq/ nseq, nseq_select, extension

      end module vmec_seq

      module bootsj_input
      use kind_spec
      integer :: nrho, mbuse, nbuse, isymm0
      real(rprec) :: tempres, ate(0:11), ati(0:11)
      real(rprec) :: zeff1, dens0, teti, damp, damp_bs

      namelist /bootin/ nrho, mbuse, nbuse, zeff1, dens0, teti, tempres,
     1   damp, damp_bs, isymm0, ate, ati

      contains

      subroutine read_boot_namelist (iunit, istat)
      integer :: iunit, istat

      read (iunit, nml=bootin, iostat=istat)

      end subroutine read_boot_namelist

      end module bootsj_input

      module optim_params
      use vparams, only: rprec, dp, nsd, ntord, ntor1d, mpol1d
      use vsvd0, only: nigroup
      implicit none
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: NumJstard = 10
      integer, parameter :: ini_max=100                                  !! COBRA
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: niter_opt, NumJstar, NumJinvariant,
     1   mboz_opt, nboz_opt,
     2   nbmn, nproc, jboot,                                              !LPK
     3   NS_JConf_Src, NS_JConf_Tgt, NPitch_Jconf
      integer :: num_processors, num_levmar_params
      integer, dimension(20) :: n_jac, m_jac
      integer, dimension(20) :: n_vac_island, m_vac_island
      real(rprec) :: r00_scale, b00_scale, epsfcn, rgrid_min, rgrid_max,
     1               zgrid_min, zgrid_max, r00_opt
      real(rprec) :: sigma_jstar(nsd, NumJstard),
     1   sigma_jinvariant(nsd, NumJstard), sigma_jconf
      real(rprec), dimension(nsd) ::
     1   sigma_iota, sigma_mercier, sigma_vp
      real(rprec), dimension(nsd) :: sigma_bmin,
     1   sigma_bmax, sigma_ripple, sigma_bmn, nsurf_mask
      real(rprec) :: sigma_aspect, sigma_ellipticity, sigma_maxcurrent,
     1   sigma_coil_complex, sigma_curv, sigma_beta, target_aspectratio,
     2   target_maxcurrent, target_beta, target_rmax, target_rmin,
     3   target_ellipticity, sigma_rmax, sigma_rmin, sigma_iota_max,
     4   sigma_iota_min, target_iota_max, target_iota_min,
     5   sigma_iota_max_min, target_iota_max_min,
     6   sigma_berr_avg, sigma_berr_max, target_eplasma, sigma_eplasma,
     7   target_curtor, sigma_curtor,
     9   sigma_centering                                               !!Obsolete, replaced by sigma_rmax, sigma_rmin
      real(rprec) :: coil_separation
      real(rprec), dimension(0:10) :: target_iota, target_well,
     1   at, aseedcur, fboot                                             !!LPK
      real(rprec) :: sigma_bal, sigma_boot, zeff_boot,                   !!LPK
     1    sigma_fluxp, target_fluxp, sigma_zmax, target_zmax,
     2    sigma_pseudo, sigma_pseudo2, sigma_rbtor, target_rbtor
      real(rprec), dimension(9) :: sigma_kink, target_kink
      real(rprec), dimension(20) :: sigma_jac, sigma_vac_island
      real(rprec), dimension(60) :: phi_lim
      real(rprec), dimension(40,60) :: r_lim, z_lim

      complex :: helicity
      logical, dimension(-ntord:ntord) :: lfix_ntor
      logical, dimension(nsd) :: lsurf_mask, ldkes_mask
      logical :: lextcur(nigroup)
      logical :: lreset_opt, lbmn, lcoil_complex, lboundary,
     1   lcur_prof_opt, liota_prof_opt, lbootsj_opt, lj_star,
     2   lj_invariant, lcoil_opt, lcoil_geom, laspect_max, lbeta_min,
     3   lbal_opt, lbootstrap, lseedcur, lpress_opt, lkink_opt,          !!LPK
     4   lnescoil_opt, lvv_tgt_min, lballoon_flip, lpseudo_sin
     5  ,lcurprof_opt, lprof_opt                                         !!Obsolete, replaced by lcur_prof_opt, liota_prof_opt
      character*(100) :: seq_ext, opt_ext
      character*(4) :: sym_type

      real(rprec) :: target_coil_complex, target_coil_jmax               !!LPK
      real(rprec) :: sigma_coil_jmax, sigma_berr_ave                     !!LPK, SPH

      real(rprec), dimension(nsd) :: sigma_neo, nneo_mask                !! NEO
      logical, dimension(nsd) :: lneo_mask                               !! NEO
      logical :: lneo_opt                                                !! NEO
      real(rprec), dimension(nsd) :: sigma_dsubr
      logical :: ldsubr_opt
      real(rprec) :: sigma_orbit
      logical :: lorbit_opt
      integer :: nopt_alg, nopt_boundary
      logical :: ldiag_opt, lkeep_mins                                   !! Diagnostic output

      real(rprec) :: sigma_kappa, target_kappa, sigma_oh
      real(rprec), dimension(nigroup) :: sigma_extcur, target_extcur,
     1                                   oh_coefs
      real(rprec), dimension(0:10) :: target_iota_p
      real(rprec), dimension(nsd) :: sigma_iota_pmax, sigma_iota_pmin
      real(rprec) :: sigma_jedge, target_jedge

      integer :: nini_theta, nini_zeta, nini_tot                         !! COBRA
      real(rprec), dimension(nsd) :: sigma_balloon, sigma_pgrad,
     1                               target_balloon                      !! COBRA
      real(rprec), dimension(ini_max):: bal_theta0, bal_zeta0            !! COBRA

      real(rprec) :: sigma_pedge(1)                                      !! COBRA
      real(rprec) :: sigma_bootsj(nsd)
      logical :: lballoon_mask(nsd)                                      !! COBRA
      logical :: lballoon_opt, lpres_prof_opt, lpres_opt_edge0,          !! COBRA
     1           lpres_opt_edgegr0, lcur_opt_edge0
      integer :: pres_opt_nmax

      real(rprec) :: nballoon_mask(nsd)                                  !! VMECCOBRA (RS)
      logical :: l_legendre                                              !! LEGENDRE (RS)

      real(rprec), dimension(nsd) :: ndkes_mask, dkes_nu, dkes_efield,
     1     sigma_dkes                                                    !! RHF
      logical :: ldkes_opt                                               !! RHF

      real(rprec) :: sigma_vv, sigma_vv_rms, vv_dist, vv_dist_rms,
     1           target_vv, target_vv_rms, sigma_vv_max, vv_dist_max     !! RH & MZ
      integer :: mpol_vv, ntor_vv, nu_vv, nv_vv
      real(rprec), dimension(-ntord:ntord,0:mpol1d) :: rbc_vv, zbs_vv

      real(rprec) :: sigma_bd, sigma_bd_rms, bd_dist, bd_dist_rms,
     1           target_bd, target_bd_rms, sigma_bd_max, bd_dist_max     !! MZ
      integer :: mpol_bd, ntor_bd, nu_bd, nv_bd
      real(rprec), dimension(-ntord:ntord,0:mpol1d) :: rbc_bd, zbs_bd

      logical :: shapeweight
C     !       shape deviation weighting   (EAL)
      real(rprec) :: theta0_bw(3), phi0_bw,
     1   wtheta_bw, wphi_bw, amplw_bw(3),  planes_bw(3)

      logical :: ldiagno_opt                                             !! DIAGNO mcz
      character(120) :: diagno_control
      real(rprec), dimension(100) :: sigma_diagno_seg,                   !! diagno mcz
     1         target_diagno_seg, sigma_diagno_flx, target_diagno_flx    !! diagno mcz
      integer :: ndiagno_seg, ndiagno_flx

      integer :: np_prof                        ! MCZ  pressure profile matching
      real(rprec), dimension(50) :: sigma_p_prof, r_p_prof, z_p_prof,
     1           phi_p_prof, p_prof
      real(rprec) :: sigma_p_damp, factor_p_prof
      logical :: lp_prof_incl_edge

      integer :: n_emis                                                  !! SXR chords
      real(rprec), dimension(0:10) :: aemis                              !! mcz
      character(120) :: emis_file
      real(rprec) :: sigma_emis_damp
      real(rprec), dimension(1000) :: sigma_emis, emis_chord


      namelist /optimum/ lreset_opt, lcur_prof_opt, liota_prof_opt,
     1   lbootsj_opt, lbmn, lj_star, lj_invariant, lfix_ntor, lboundary,
     2   lextcur, nopt_alg, nopt_boundary,
     3   niter_opt, num_processors, num_levmar_params,
     4   helicity, epsfcn, r00_scale, b00_scale, r00_opt,
     5   target_aspectratio, target_iota, mboz_opt, nboz_opt, NumJstar,
     6   NumJInvariant, nsurf_mask, target_kink, target_eplasma,
     7   target_maxcurrent, target_beta, target_rmax, target_rmin,
     8   target_ellipticity, sigma_jstar, sigma_jinvariant, sigma_rmax,
     9   sigma_rmin, sigma_iota, sigma_aspect, sigma_ellipticity,
     A   sigma_mercier, sigma_kink, sigma_bmin, sigma_bmax, sigma_bmn,
     B   sigma_ripple, sigma_beta, sigma_maxcurrent, sigma_eplasma,
     C   sigma_curv, sigma_coil_complex, lcoil_complex, coil_separation,
     D   laspect_max, lbeta_min, sigma_vp, sym_type, sigma_bootsj,
     E   sigma_balloon, sigma_pgrad, sigma_pedge, lballoon_opt,          !! COBRA
     F   target_balloon, lpres_prof_opt, bal_theta0, bal_zeta0,          !! COBRA
     G   lpres_opt_edge0, lpres_opt_edgegr0,nballoon_mask,lballoon_mask, !! VMECCOBRA (RS)
     H   l_legendre, rgrid_min, rgrid_max, zgrid_min, zgrid_max,         !! LEGENDRE (RS)
     I   nproc, sigma_zmax, target_zmax, sigma_pseudo, sigma_pseudo2,    !! LPK
     J   lpseudo_sin, lbal_opt, sigma_bal, nbmn, zeff_boot,              !! LPK
     K   lbootstrap, jboot, lseedcur, fboot, sigma_boot, at, aseedcur,   !! LPK
     L   sigma_fluxp, target_fluxp, lpress_opt, lkink_opt, lkeep_mins,   !! LPK
     M   lnescoil_opt, lcoil_geom, target_coil_complex,target_coil_jmax,
     N   sigma_coil_jmax, sigma_berr_ave, ldkes_mask,                    !! LPK, SPH
     J   sigma_neo, lneo_mask, nneo_mask, lneo_opt, ldiag_opt,           !! NEO
     J   sigma_dsubr, ldsubr_opt, sigma_orbit, lorbit_opt,
     O   ldkes_opt, ndkes_mask, sigma_dkes, dkes_nu, dkes_efield,        !! RHF
     P   sigma_vv, sigma_vv_rms, mpol_vv, ntor_vv, rbc_vv, zbs_vv,       !! MZ
     Q   target_rbtor, sigma_rbtor, sigma_kappa, target_kappa,           !! MZ
     R   sigma_extcur, target_iota_p, sigma_iota_pmax,sigma_iota_pmin,   !! MZ
     S   shapeweight, wtheta_bw, wphi_bw, theta0_bw, phi0_bw, amplw_bw,  !!EAL
     T   planes_bw, sigma_jac, m_jac, n_jac,sigma_oh, oh_coefs,
     U   m_vac_island, n_vac_island, sigma_vac_island, sigma_vv_max,
     V   target_vv, target_vv_rms, nu_vv, nv_vv, sigma_iota_max,
     W   sigma_iota_min, target_iota_max, target_iota_min, lvv_tgt_min,
     X   sigma_berr_avg, sigma_berr_max, lballoon_flip, sigma_p_damp,
     Y   sigma_iota_max_min, target_iota_max_min, phi_lim, r_lim, z_lim,
     1   sigma_jconf, ns_jconf_src, ns_jconf_tgt, npitch_jconf,
     2   np_prof, p_prof, r_p_prof, z_p_prof, phi_p_prof, sigma_p_prof,
     3   factor_p_prof,n_emis, aemis, emis_chord, emis_file, sigma_emis,
     4   sigma_emis_damp, target_extcur, ldiagno_opt, sigma_diagno_seg,
     5   target_diagno_seg, sigma_diagno_flx, target_diagno_flx,
     6   diagno_control, lcur_opt_edge0, sigma_curtor, target_curtor,
     7   pres_opt_nmax, target_jedge, sigma_jedge, lp_prof_incl_edge,
     8   sigma_bd, sigma_bd_rms, mpol_bd, ntor_bd, rbc_bd, zbs_bd,       !! MZ
     9   target_bd, target_bd_rms, nu_bd, nv_bd, sigma_bd_max
     Z   ,lsurf_mask, lcurprof_opt, lprof_opt, target_well, lcoil_opt    !!Obsolete, retain for consistency
     Z   ,sigma_centering
      contains

      subroutine read_optimum_namelist (iunit, istat)
      integer :: iunit, istat

      read (iunit, nml=optimum, iostat=istat)

      end subroutine read_optimum_namelist

      end module optim_params


      module system_mod
      interface system
         subroutine vmec_system(cmd, error)
         character*(*), intent(in) :: cmd
         integer, optional :: error
         end subroutine vmec_system
      end interface   !system
      interface chdir
         integer function vmec_chdir(path)
         character*(*), intent(in) :: path
         end function vmec_chdir
      end interface   !chdir
      interface getenv
         subroutine vmec_getenv(ename, evalue)
         character*(*) :: ename, evalue
         end subroutine vmec_getenv
      end interface   !getenv
      interface putenv
         subroutine vmec_putenv(ename, evalue, ierror)
         character*(*) :: ename, evalue
         integer :: ierror
         end subroutine vmec_putenv
      end interface   !putenv

#if !defined(CRAY) && !defined(IRIX64)
      interface
         subroutine pxffork(ipid, ierror)
         integer :: ipid, ierror
         end subroutine pxffork
      end interface
      interface
         subroutine pxfgetpid(ipid, ierror)
         integer :: ipid, ierror
         end subroutine pxfgetpid
      end interface
      interface
         subroutine pxfwait(istat, iretpid, ierror)
         integer :: istat, iretpid, ierror
         end subroutine pxfwait
      end interface
#endif
      end module system_mod

      module fdjac_mod
      use kind_spec
      integer :: ncnt, max_processors, num_lm_params, ix_min,
     1           jac_count
c      real(rprec), dimension(:), allocatable :: xp, wap
      logical, dimension(:), allocatable :: flip
      integer, dimension(:), allocatable :: jac_order
      real(rprec) :: eps
      end module fdjac_mod

      module lmpar_mod
      use kind_spec
      integer :: nscan, ldfjac, lev_state, m, n, myid, numprocs
      integer, dimension(:), allocatable :: ipvt
      real(rprec) :: pnorm, fnorm1, delta, par, spread_ratio
      real(rprec), dimension(:), allocatable :: x, wa1, wa2, wa3, wa4
      real(rprec), dimension(:), allocatable :: diag, qtf, fvec
      real(rprec), dimension(:,:), allocatable :: fjac
      logical :: first
      end module lmpar_mod

      module ga_mod
      use kind_spec
      use gade_mod
      implicit none

      integer, dimension(nparmax) :: ig2
      integer, dimension(nchrmax,indmax) :: iparent, ichild
      integer, dimension(nchrmax) :: ibest
      real(rprec), dimension(nparmax,indmax) :: parent, child
      real(rprec), dimension(indmax) :: fitness
      real(rprec), dimension(nparmax) :: g0, g1,
     +                     pardel, par_max, par_min
      real(rprec), dimension(max_gen) :: geni, genavg, genmax
      integer :: nparam, nchrome, num_obj
      integer :: jbest,irestrt
      integer :: maxgen,nfit_eval,
     +                 kountmx
      integer :: iunit_ga_restart, iunit_ga_out
      real(rprec), dimension(:), pointer :: f_obj

      end module ga_mod

      module de_mod
      use kind_spec
      use gade_mod
      implicit none
      integer :: n_pop, n_free, nopt, nfev
      real(rprec), dimension(:,:), allocatable :: ui_XC

      end module DE_mod

      module read_boozer_mod
!
!     USER NOTES:
!
!     USE READ_BOOZ_MOD to include variables dynamically allocated
!     in the module
!     CALL DEALLOCATE_READ_BOOZER to free this memory when it is no longer needed
!
      use kind_spec
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: mnboz_b, mboz_b, nboz_b, nfp_b, ns_b
      integer, dimension(:), allocatable :: idx_b, ixm_b, ixn_b
      real(rprec), dimension(:), allocatable :: iota_b, pres_b,
     1    phip_b, phi_b, beta_b, buco_b, bvco_b
      real(rprec) :: aspect_b, rmax_b, rmin_b, betaxis_b
      real(rprec), dimension(:,:), allocatable ::
     1   bmn_b, rmnc_b, zmns_b, pmns_b, gmn_b

      interface
         subroutine read_boozer_file(extension, ierr, iopen)
         integer :: ierr
         integer, optional :: iopen
         character*(*) :: extension
         end subroutine read_boozer_file
      end interface

      contains

      subroutine read_boozer_deallocate
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: istat
!-----------------------------------------------

      if (allocated(iota_b)) deallocate (iota_b, pres_b, beta_b,
     1    phip_b, phi_b, bvco_b, buco_b, idx_b, stat = istat)

      if (allocated(bmn_b)) deallocate (bmn_b, rmnc_b,
     1   zmns_b, pmns_b, gmn_b, ixm_b, ixn_b, stat = istat)

      end subroutine read_boozer_deallocate

      end module read_boozer_mod

      module read_wout_mod
!
!     USE READ_WOUT_MOD to include variables dynamically allocated
!     in the module
!     CALL DEALLOCATE_READ_WOUT to free this memory when it is no longer needed
!
      use kind_spec
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nfp, ns, mpol, ntor, mnmax, itfsq, niter, iasym,
     1    ireconstruct, ierr_vmec, imse, itse, nbsets, nobd, nextcur,
     2    nstore_seq, isnodes, ipnodes, nsets, nparts, nlim, nrgrid,
     3    nzgrid, imatch_phiedge, nbfldn, isigng
      integer, allocatable, dimension(:) :: limitr, nbfld, nsetsn
      real(rprec) :: wb, wp, gamma, pfac, rmax_surf, rmin_surf,
     1    zmax_surf, aspect, betatot, betapol, betator, betaxis, b0,
     2    tswgt, msewgt, flmwgt, bcwgt, phidiam, version_,
     3    delphid, rx1, rx2, zy1, zy2, condif, IonLarmor, VolAvgB,
     4    Aminor, Rmajor, Volume, RBtor, RBtor0, Itor
      real(rprec), dimension(:,:), allocatable ::
     1    rmnc, zmns, lmns, rmns, zmnc, lmnc, bmn, gmn, bsubumn,
     2    bsubvmn,bsubsmn, bsupumn, bsupvmn, currvmn, bcoil, plbfld,
     3    bbc, rlim, zlim, raxis, zaxis
      real(rprec), dimension(:), allocatable ::
     1   iotas, mass, pres, beta_vol, xm, xn,
     2   phip, buco, bvco, phi, vp, overr, jcuru, jcurv, specw,
     3   jdotb, bdotgradv, fsqt, wdot,
     3   Dmerc, Dshear, Dwell, Dcurr, Dgeod, equif, extcur,
     4   sknots, ystark, y2stark, pknots, ythom, y2thom,
     5   anglemse, rmid, qmid, shear, presmid, alfa, curmid, rstark,
     6   qmeas, datastark, rthom, datathom, dsiext, plflux, dsiobt
      real(rprec), dimension(:,:,:), allocatable :: pfcspec
      character :: mgrid_file*100, tokid*60, input_extension*100
      character*(8), dimension(:), allocatable :: curlabel

!     OVERLOAD SUBROUTINE READ_WOUT_FILE TO ACCEPT BOTH UNIT NO. (OPENED EXTERNALLY)
!     OR FILENAME (HANDLE OPEN/CLOSE HERE)

      interface read_wout_file
         subroutine readw_and_open(filename, ierr, iopen)
            integer :: ierr
            integer, optional :: iopen
            character*(*) :: filename
         end subroutine readw_and_open
         subroutine readw_only(iunit, ierr, iopen)
            integer :: iunit, ierr
            integer, optional :: iopen
         end subroutine readw_only
      end interface

      contains

      subroutine read_wout_deallocate
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: istat
!-----------------------------------------------

      if (allocated(xm)) deallocate (nbfld, xm, xn, rmnc, zmns, lmns,
     1  rmns, zmnc, lmnc, bmn, gmn, bsubumn,
     2  bsubvmn, bsubsmn, bsupumn, bsupvmn, currvmn, iotas, mass,
     3  pres, beta_vol, phip, buco, bvco, phi, vp, overr, jcuru,
     4  jcurv, specw, Dmerc, Dshear, Dwell, Dcurr, Dgeod, equif, jdotb,
     5  bdotgradv,extcur, curlabel, raxis, zaxis, fsqt, wdot,
     6  stat = istat)

      if (ireconstruct.gt.0 .and. allocated(sknots)) deallocate (
     1    ystark, y2stark, pknots, anglemse, rmid, qmid, shear,
     2    presmid, alfa, curmid, rstark, datastark, rthom, datathom,
     3    ythom, y2thom, plflux, dsiobt, bcoil, plbfld, bbc, sknots,
     4    pfcspec, limitr, rlim, zlim, nsetsn, stat = istat)


      end subroutine read_wout_deallocate

      end module read_wout_mod

      module read_namelist_mod

      interface read_namelist
         subroutine read_namelist(iunit, io_stat, lc_name)
         use vmec_input
         use vmec_seq
         use bootsj_input
         use optim_params
         integer iunit, io_stat
         character*(*) :: lc_name
         end subroutine read_namelist
      end interface


      interface autofill
         subroutine autofill(X, DESX, NX)
         use kind_spec
         integer :: nx
         real(rprec) :: X(NX)
         character*(*) :: DESX
         end subroutine autofill
      end interface

      end module read_namelist_mod

      module safe_open_mod
         interface safe_open
            subroutine safe_open_it(iunit, istat, filename, filestat,
     1      fileform, record, access)
            integer, intent(inout) :: iunit
            integer, intent(out) :: istat
            integer, intent(in), optional :: record
            character*(*), intent(in) :: filename, filestat, fileform
            character*(*), intent(in), optional :: access
            end subroutine safe_open_it
         end interface
      end module safe_open_mod
