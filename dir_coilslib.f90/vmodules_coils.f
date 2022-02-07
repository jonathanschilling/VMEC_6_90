!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     COILOPT4 RELATED MODULES 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module cparms
      use kind_spec
      real(rprec), parameter :: zero = 0, one = 1,
     1 dpi = 3.141592653589793238462643_dp, twopi = 2*dpi,
     2 dmu0 = 4.0e-7_dp*dpi, dtwopi = 2*dpi
      end module cparms

      module control_mod
      integer :: nopt_alg
      logical :: lrestart
      end module control_mod

      module Vcoilpts
      use cparms
      integer, parameter :: ncdim = 100
      integer, parameter :: nwdim = 400
      integer, parameter :: nwdim1 = nwdim + 1
      integer, parameter :: nfdim = 5
      end module Vcoilpts

      module modular_coils
      use Vcoilpts
      integer, parameter :: nfourier=20, nsurf_p=100
      integer :: nf_rho, nf_phi, nstep, niter_opt, ncoils,
     1    nmod_coils, nmid, nodd, nfper, nmod_coils_per_period,
     2    m_in, n_in, nmod_unique_coils, nmod_coeffs, nmod_currents
      integer :: numsurf, nall
      integer, dimension(100) :: niter_array
      real(rprec), target :: dcp_wgt
      real(rprec) :: epsfcn, dcp_exp, dcp_tgt
      real(rprec), dimension(100) :: epsfcn_array
      real(rprec), dimension(ncdim), target ::
     1                  dcc_wgt, dcc_exp, dcc_tgt, cc_min
      real(rprec), dimension(ncdim), target ::
     1                  rc_wgt, rc_exp, rc_tgt, rc_min
      real(rprec), dimension(ncdim), target ::
     1                  lmod_wgt, lmod_tgt, cmod_scl
      real(rprec), dimension(ncdim), target ::
     1                  ymin_wgt, ymin_tgt, ymin_cls
      real(rprec), dimension(ncdim), target ::
     1                  cu_wgt, cu_tgt, cu_sum
      real(rprec), dimension(ncdim), target :: r_ext
      real(rprec), dimension(ncdim), target :: curmod
      real(rprec), dimension(ncdim) :: phimin, phimax
      real(rprec), dimension(ncdim) :: mod_length
      real(rprec), dimension(nwdim1,3,ncdim) :: x_mod, y_mod, z_mod
      real(rprec), dimension(ncdim) :: curcon
      real(rprec), dimension(ncdim,0:nfourier) :: phic, phis
      real(rprec), dimension(ncdim,0:nfourier) :: rhoc, rhos
      real(rprec), dimension(ncdim,nwdim1) :: rho , phi
      real(rprec), dimension(ncdim,nwdim) :: rcoil, zcoil
      integer, dimension(nsurf_p) :: m_num,n_num
      real(rprec), dimension(nsurf_p) :: rmn_sf, zmn_sf
      real(rprec), dimension(ncdim) :: phi_full
      real(rprec) :: p_d_min, p_d_max
      real (rprec), dimension (ncdim) :: b_max
      logical :: lmodular
      logical :: lmodcur
      logical :: lsurfv
      logical :: lncsx
      logical :: lsymm
      end module modular_coils

      module vf_coils
      use Vcoilpts
      integer :: num_vf, nvf, nrvf_c, nvf_coeffs
      real(rprec), dimension(ncdim) :: rc_vf, zc_vf, cc_vf
      real(rprec), dimension(ncdim,ncdim) :: rcfc_vf, rcfs_vf
      real(rprec), dimension(ncdim) :: rvf, zvf, cvf, cvf_tgt
      real(rprec), dimension(ncdim), target :: cvf_wgt, rvf_wgt
      real(rprec), dimension(ncdim) :: rvf_tgt, rvf_max
      real(rprec), dimension(nwdim1,3,ncdim) :: x_vf, y_vf, z_vf
      real(rprec) :: cvf_ssq
      logical :: lvf, lvfc, lvfvar, lvfr, lvfz
      logical, dimension(ncdim) :: lcc_vf
      end module vf_coils

      module boundary
      use cparms
      integer :: mpol, ntor, nfp, mnmax, ns, nu, nv, nuv, nedge
      integer, allocatable, dimension(:) :: ixm, ixn
      real(rprec) :: TotalArea, sum_d_area, rbphi_avg
      real(rprec) :: iota_b, phip_b
      real(rprec), dimension(:), allocatable :: phib, thetab
      real(rprec), dimension(:), allocatable :: rb, zb
      real(rprec), dimension(:), allocatable :: rb_th, rb_ph
      real(rprec), dimension(:), allocatable :: zb_th, zb_ph
      real(rprec), dimension(:), allocatable :: n_r, n_phi, n_z
      real(rprec), dimension(:), allocatable :: d_area
      real(rprec), dimension(:), allocatable :: x_p, y_p, z_p
      real(rprec), dimension(:), allocatable :: phi_d, theta_d
      real(rprec), dimension(:), allocatable :: xm, xn
      real(rprec), dimension(:), allocatable :: rmnc_b, zmns_b
      real(rprec), dimension(:), allocatable :: gmn_b, lmns_b
      real(rprec), dimension(:), allocatable :: bnormal_match
      real(rprec), dimension(:), allocatable :: b_error, b_mod
      character boundary_mn_file*30
      end module boundary

      module tor_field
      use cparms
      use Vcoilpts
      real(rprec) :: i_pol, i_tfc, pol_cur
      real(rprec), target :: dpc_wgt
      integer :: mtfcoil, mtfwire
      real(rprec), dimension(ncdim) :: tfc_cur
      real(rprec), dimension(ncdim, nwdim1) ::
     1  tfc_x, tfc_y, tfc_z
      logical :: ltfc, ltfcv
      logical :: lqos
      logical :: lpolcur
      end module tor_field

      module bnorm_mod
      use cparms
      use Vcoilpts
      integer, parameter :: bnorm_dim=1000
      integer :: mnbn_max
      real(rprec), dimension(bnorm_dim) :: xbn_m, xbn_n, bn_coef
      integer :: mmax_bmn, nmax_bmn, mnmax_bmn
      real(rprec), dimension(:), allocatable :: xm_bmn, xn_bmn
      real(rprec), dimension(:), allocatable :: bmn, bmn_error
      real(rprec), dimension(:), allocatable :: luv, bsl_error
      character*200 :: bnorm_file
      logical :: lbnorm
      end module bnorm_mod

      module bcoils_mod
      use Vcoilpts
      integer :: mbcoils, mc_max
      integer, dimension(ncdim) :: mbwires, mc_bg
      real(rprec) :: mxb_wgt
      real(rprec), dimension(ncdim) :: bcoil_cur, cc_bg
      real(rprec), dimension(ncdim, nwdim1) ::
     1  bcoil_x, bcoil_y, bcoil_z
      logical :: lbcoil, lbcoil_cur
      logical :: lp_bg(ncdim)
      integer :: n_access, n_access_pts
      real(rprec), dimension(ncdim) :: x0_access, y0_access, z0_access
      real(rprec), dimension(ncdim) :: x1_access, y1_access, z1_access
      real(rprec), dimension(ncdim) :: dac_exp, dac_tgt, acc_min
      real(rprec), dimension(ncdim), target :: dac_wgt
      real(rprec), dimension(:, :), allocatable ::
     1  x_access, y_access, z_access
      logical :: laccess
      character*200 :: bcoil_file
      end module bcoils_mod


      module saddle_coils
      use Vcoilpts
      integer, parameter :: mfourier = 100
      integer :: nsad_v, nsad_u, nsad_coeffs, nfils
      integer :: nsad_coils_per_period, nsad_coils, nsmid, nsodd
      integer :: nsad_unique_coils
      integer :: nsad_currents, num_cursad
      integer, dimension(ncdim) :: nsad_group
      logical :: ls_cur(ncdim)
      real(rprec), dimension(ncdim) :: csad_scl
      real(rprec), dimension(ncdim,0:mfourier) :: sad_v_c,
     1                  sad_v_s, sad_u_c, sad_u_s
      real(rprec), dimension(ncdim) :: sad_u0, sad_v0,
     1                  sad_phi0, sad_theta0
      real(rprec), dimension(ncdim) :: cursad, c_sad, sad_length
      real(rprec), dimension(ncdim), target ::
     1                  dsc_wgt, dsc_exp, dsc_tgt, sc_min
      real(rprec), dimension(ncdim), target ::
     1                  rs_wgt, rs_exp, rs_tgt, rs_min
      real(rprec), dimension(ncdim), target ::
     1                  cs_wgt, cs_tgt, cs_sum
      real(rprec), dimension(ncdim), target ::
     1                  dscxp_wgt, dscxp_exp, dscxp_tgt, scxp_min
      real(rprec), dimension(ncdim), target ::
     1                  lsad_wgt, lsad_tgt
      real(rprec), dimension(ncdim), target ::
     1                  rmax_sad, rmax_wgt, rmax_tgt
      real(rprec), dimension(ncdim), target :: ymin_sad
      real(rprec), dimension(:,:), allocatable :: u_sad, v_sad
      real(rprec), dimension(:,:,:), allocatable :: x_sad, y_sad,
     1                  z_sad
      real(rprec) :: deln, delt, p_s_min, bkp_min, bkp_wgt, bkp_tgt
      real(rprec), dimension(ncdim), target :: csc_wgt, scd_wgt
      real(rprec), dimension(ncdim) :: csc_tgt, scd_tgt
      logical :: lsaddle, lsadsfv, lsadcur, lsmod, lsadshape, 
     1           lspline, lsplbkp
      logical, dimension(0:mfourier) :: lsad_uv_m
      integer, dimension(ncdim,0:mfourier) :: nvar_vc, nvar_uc
      end module saddle_coils

      module saddle_surface
      use cparms
      integer, parameter :: nsurf = 500
      integer :: numsurf_sad, nopt_wsurf
      integer, dimension(nsurf) :: m_sad, n_sad
      integer :: mpol_opt, ntor_opt, irho_bdy, irm0_bdy, izm0_bdy
      real(rprec), dimension(nsurf) :: rmn_sad, zmn_sad
      integer, dimension(:), allocatable :: nbrho_opt, mbrho_opt
      integer, allocatable :: nrz0_opt(:)
      real(rprec), dimension(:,:), allocatable :: rbc,zbs,rhobc,delta_mn
      end module saddle_surface


      module coilsnamin
      use modular_coils
      use saddle_coils
      use saddle_surface
      use vf_coils
      use tor_field
      use bcoils_mod
      use bnorm_mod
      use control_mod
 
      namelist /coilsin/ nmod_coils_per_period, nf_phi, nf_rho, epsfcn,
     1   lvf, lmodcur, lsurfv, lbnorm, rhoc, rhos, phic, phis, curmod, 
     2   dcc_wgt, dcc_exp, dcc_tgt, dcp_wgt, dcp_exp, dcp_tgt, rc_wgt,
     3   rc_exp, rc_tgt, lmod_wgt, lmod_tgt, niter_opt, nstep, i_pol,
     4   numsurf, m_num, n_num, rmn_sf, zmn_sf, num_vf, rc_vf, zc_vf,
     5   cc_vf, lbcoil, lncsx, lsymm, lsaddle, nsad_coils_per_period,
     6   ltfc, ltfcv, i_tfc, lsadsfv, nsad_u, nsad_v, nfils, sad_v_c,
     7   sad_v_s, sad_u_c, sad_u_s, sad_v0, sad_u0, lsadcur, lpolcur,
     8   cursad, numsurf_sad, m_sad, n_sad, rmn_sad, zmn_sad, dsc_wgt,
     9   dsc_exp, dsc_tgt, r_ext, ymin_wgt, ymin_tgt, lmodular, lqos,
     1   rs_wgt, rs_exp, rs_tgt, lsmod, laccess, n_access, x0_access,
     2   y0_access, z0_access, x1_access, y1_access, z1_access,
     3   dac_wgt, dac_exp, dac_tgt, cvf_wgt, cvf_tgt, cs_wgt, cs_tgt,
     4   cu_wgt, cu_tgt, dpc_wgt, mc_bg, lp_bg, bcoil_cur, dscxp_wgt,
     5   dscxp_exp, dscxp_tgt, deln, delt, lspline, lsplbkp, nvar_vc,
     6   nvar_uc, bkp_wgt, bkp_tgt, mxb_wgt, lvfvar, nrvf_c, rcfc_vf,
     7   rcfs_vf, lvfr, lvfz, nopt_alg, lrestart, nopt_wsurf, rvf_wgt,
     8   rvf_tgt, nsad_group, ls_cur, csad_scl, lsad_wgt, lsad_tgt,
     9   rmax_wgt, rmax_tgt, bcoil_file, lbcoil_cur, lsadshape,
     1   lvfc, lcc_vf, csc_wgt, csc_tgt, scd_wgt, scd_tgt,
     2   niter_array, epsfcn_array, lsad_uv_m

      contains

      subroutine read_coils_namelist (iunit, istat)
      integer :: iunit, istat

      read (iunit, nml=coilsin, iostat=istat)

      end subroutine read_coils_namelist
      
      end module coilsnamin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     END COIL OPT MODULES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
