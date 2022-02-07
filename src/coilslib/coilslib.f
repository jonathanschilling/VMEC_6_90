      subroutine init_coilgeom (nfper) 
      use boundary, ONLY: nfp
      use modular_coils, ONLY: nmod_coils
      use saddle_coils, ONLY: nsad_coils
      implicit none
      
      integer :: nfper
      
      nfp = nfper
      nmod_coils = -1
      nsad_coils = -1
      
      end subroutine init_coilgeom


      subroutine initialize_coilsin
      use kind_spec
      use coilsnamin

      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
!-----------------------------------------------
      bcoil_file = "bcoilxyz.dat"
      lvf = .false.
      lvfc = .true.
      lcc_vf = .true.
      lvfvar = .false.
      lvfr = .false.
      lvfz = .false.
      ltfc = .false.
      ltfcv = .false.
      lbcoil = .false.
      lbcoil_cur = .true.
      lbnorm = .false.
      lmodcur = .false.
      lsadcur = .false.
      lsadshape = .true.
      lmodular = .false.
      lsaddle = .true.
      lsmod = .false.
      lspline = .false.
      lsplbkp = .false.
      lncsx = .false.
      lqos = .false.
      lsymm = .true.
      laccess = .false.
      lpolcur = .false.
      lp_bg = .true.
      lrestart = .false.
      ls_cur = .true.
 
      n_access = 0
      nf_phi = 0
      nf_rho = 0
      nmod_coils_per_period = 0
      nsad_coils_per_period = 0
      nsad_u = 0
      nsad_v = 0
      numsurf = 0
      numsurf_sad = 0
      nstep  = 10
      niter_opt = 1
      nopt_alg = 0
      nopt_wsurf = -1
      nfils = 1
      nvar_vc = 1
      nvar_uc = 1
  
      i_pol = 0
      i_tfc = 0
      rhoc = 0
      rhos = 0
      phic = 0
      phis = 0
      curmod = 0
      dcc_wgt = 0
      dcc_exp = 0
      dcc_tgt = 0
      rc_wgt = 0
      rc_exp = 0
      rc_tgt = 0
      lmod_wgt = 0
      lmod_tgt = 0
      lsad_wgt = 0
      lsad_tgt = 0
      sad_v_c = 0
      sad_v_s = 0
      sad_u_c = 0
      sad_u_s = 0
      sad_v0 = 0
      sad_u0 = 0
      lsad_uv_m = .true.
      cursad = 0
      r_ext = 0
      ymin_wgt = 0
      ymin_tgt = 0
      rmax_wgt = 0
      rmax_tgt = 0
      rs_wgt = 0
      rs_exp = 0
      rs_tgt = 0
      dac_wgt = 0
      dac_exp = 0
      dac_tgt = 0
      cvf_wgt = 0
      cvf_tgt = 0
      rvf_wgt = 0
      rvf_tgt = 0
      rc_vf = 0
      zc_vf = 0
      rcfc_vf = 0
      rcfs_vf = 0
      cs_wgt = 0
      cs_tgt = 0
      csc_wgt = 0
      csc_tgt = 0
      scd_wgt = 0
      scd_tgt = 0
      dpc_wgt = 0
      mc_bg = 0
      bcoil_cur = 0
      dscxp_wgt = 0
      dscxp_exp = 0
      dscxp_tgt = 0
      delt = 0.032_dp
      deln = 0.055_dp           ! NCSX value
      bkp_wgt = 10000
      bkp_tgt = 0.02_dp
      mxb_wgt = 0
      csad_scl = 1
      nmid = 0
      nsmid = 0
      n_access = 0
      num_vf = 0
      nrvf_c = 0

      do i = 1, ncdim
         nsad_group(i) = i
      end do
 
      end subroutine initialize_coilsin


      subroutine init_modular_coils (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use boundary, only: nfp
      use modular_coils
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nc, i, j, k, n
      integer :: nvariables, modes
      real(rprec) :: xvariables(*)
!-----------------------------------------------

      nc = nmod_coils_per_period

      if (nmod_coils .lt. 0) then
         nfper = nfp
         nmod_coils = nc*nfper
         ncoils = nmod_coils
         nmid = (nc+1)/2
         nodd = mod(nc,2)
         if ((nodd.eq.0) .and. (.not.lsymm)) nmid = nmid + 1
      end if

      nvariables = 0
      nmod_coeffs = 0
 
      if (nc .le. 0) return
 
!     Initialize the variables to values of unique coil parameters
!     and count the number of variables
 
      n = 0
 
!     First consider the case with coils on both symmetry planes at
!     phi = 0 and phi = pi/nfp (lsymm = F). This implies that the number
!     of coils per field period must be even (nodd = 0).
 
      if ((nodd .eq. 0) .and. (.not.lsymm)) then
 
!     Symmetry coil at phi = 0.
 
        i = 1
        do modes = 1, nf_phi
          n = n + 1
          xvariables(n) = phis(i,modes)
        end do
 
        do modes = 1, nf_rho
          n = n + 1
          xvariables(n) = rhos(i,modes)
        end do
 
!     Coils 2 through nmid-1.
        do i = 2, nmid-1
          modes = 0
          n = n + 1
          xvariables(n) = phic(i,modes)
          do modes = 1,nf_phi
            n = n + 1
            xvariables(n) = phic(i,modes)
            n = n + 1
            xvariables(n) = phis(i,modes)
          end do
 
          do modes = 1,nf_rho
            n = n + 1
            xvariables(n) = rhos(i,modes)
          end do
        end do
 
!     Symmetry coil at phi = pi/nfp.
        i = nmid
        do modes = 1, nf_phi
          n = n + 1
          xvariables(n) = phis(i,modes)
        end do
 
        do modes = 1, nf_rho
          n = n + 1
          xvariables(n) = rhos(i,modes)
        end do
 
!     Next consider the cases with a coil on phi = 0 (lsymm = T), or
!     on phi = pi/nfp (lsymm = F), but not both. Then there may be an
!     even number of coils per period (nodd = 0) or an odd number of
!     coils per period (nodd = 1).
 
      else

        do i = 1, nmid-nodd
          modes = 0
          n = n + 1
          xvariables(n) = phic(i,modes)
          do modes = 1,nf_phi
            n = n + 1
            xvariables(n) = phic(i,modes)
            n = n + 1
            xvariables(n) = phis(i,modes)
          end do
 
          do modes = 1,nf_rho
            n = n + 1
            xvariables(n) = rhos(i,modes)
          end do
        end do
 
        if (nodd .eq. 1) then
          i = nmid
          do modes = 1, nf_phi
            n = n + 1
            xvariables(n) = phis(i,modes)
          end do
 
          do modes = 1, nf_rho
            n = n + 1
            xvariables(n) = rhos(i,modes)
          end do
        end if
 
      end if  ! end if ((nodd .eq. 0) .and. (lsymm .eqv. .false.))
 
      nmod_coeffs = n
      nvariables = n

      end subroutine init_modular_coils


      subroutine init_modular_currents (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use modular_coils
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nc, i, n
      integer :: nvariables
      real(rprec) :: xvariables(*)
!-----------------------------------------------
 
      nvariables = 0
      nmod_currents = 0
      nc = nmod_coils_per_period
 
      if (.not.lmodcur .or. nc.le.0) return
 
 
!     Initialize the variables to values of coil currents
!     and count the number of variables
 
      n = 0
      do i = 1, nmid
         n = n + 1
         xvariables(n) = curmod(i)
      end do
 
      nmod_currents = n
      nvariables = n

      end subroutine init_modular_currents


      subroutine load_modular_coils (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use boundary, only: nfp
      use modular_coils
      use tor_field
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, k, n
      integer :: nvariables, modes
      real(rprec) :: xvariables(*)
      real(rprec) :: csum, ctfc
!-----------------------------------------------
 
!     load the unique coil parameters with values from variables in
!     optimization
 
      nvariables = 0
      n = 0
 
!     Case with coils on both symmetry planes phi = 0 and phi = pi/nfp
!     (lsymm = F). No. of coils per field period must be even (nodd = 0).
 
      if ((nodd.eq.0) .and. (.not.lsymm)) then
 
!        Symmetry coil at phi = 0.
         i = 1
         do modes = 1, nf_phi
            n = n + 1
            phis(i,modes) = xvariables(n)
         end do
 
         rhos(i,0) = 0
         do modes = 1, nf_rho
            n = n + 1
            rhos(i,modes) = xvariables(n)
         end do
 
!        Coils 2 through nmid - 1
         do i = 2, nmid-1
            modes = 0
            n = n + 1
            phic(i,modes) = xvariables(n)
            do modes = 1,nf_phi
               n = n + 1
               phic(i,modes) = xvariables(n)
               n = n + 1
               phis(i,modes) = xvariables(n)
            end do
 
            rhos(i,0) = 0
            do modes = 1,nf_rho
               n = n + 1
               rhos(i,modes) = xvariables(n)
            end do
         end do
 
!        Symmetry coil at phi = pi/nfp
         i = nmid
         do modes = 1, nf_phi
            n = n + 1
            phis(i,modes) = xvariables(n)
         end do
 
         rhos(i,0) = 0
         do modes = 1, nf_rho
            n = n + 1
            rhos(i,modes) = xvariables(n)
         end do

!     Cases with a coil on phi = 0 (lsymm = T), or on phi = pi/nfp 
!     (lsymm = F), but not both. There may be an even number of coils
!     per period (nodd = 0), or an odd number of coils per period
!     (nodd = 1).

      else
 
         do i = 1, nmid-nodd
            modes = 0
            n = n + 1
            phic(i,modes) = xvariables(n)
            do modes = 1,nf_phi
               n = n + 1
               phic(i,modes) = xvariables(n)
               n = n + 1
               phis(i,modes) = xvariables(n)
            end do
 
            rhos(i,0) = 0
            do modes = 1,nf_rho
               n = n + 1
               rhos(i,modes) = xvariables(n)
            end do
         end do
         if (nodd .eq. 1) then
            i = nmid
            do modes = 1, nf_phi
               n = n + 1
               phis(i,modes) = xvariables(n)
            end do
 
            rhos(i,0) = 0
            do modes = 1, nf_rho
               n = n + 1
               rhos(i,modes) = xvariables(n)
            end do
         end if
 
      endif   ! end if ((nodd .eq. 0) .and. (lsymm .eqv. .false.))
 
      nvariables = n

      end subroutine load_modular_coils


      subroutine load_modular_currents (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use modular_coils
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, n
      integer :: nvariables
      real(rprec) :: xvariables(*)

!-----------------------------------------------
 
!     Load the coil currents with values from optimization variables
 
      n = 0

      if (lmodcur) then
         do i = 1, nmid
            n = n + 1
            curmod(i) = xvariables(n)
         end do
      end if
 
      nvariables = n

      end subroutine load_modular_currents


      subroutine init_saddle_coils (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use boundary, only: nfp
      use saddle_coils
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, k, n, nsv, nc, modes
      integer :: nvariables
      real(rprec) :: xvariables(*)
!-----------------------------------------------
 
      nc = nsad_coils_per_period
      nvariables = 0
      nsad_coeffs = 0

      if (nsad_coils .lt. 0) then
         nsad_coils = nc * nfp
         nsmid = nc/2
         nsodd = mod(nc,2)
         nsad_unique_coils = nsmid
      end if
 
      if (nsad_coils .le. 0) return
 
!     initialize the variables to values of unique coil parameters
!     and count the number of variables
 
      nsv = 0
 
      if (lspline) then

         ! initialize bspline series representation
         do i = 1, nsmid
            ! coefficients for v
            do n = 1, nsad_v - 3
               if (nvar_vc(i,n).ne.0) then
                  nsv = nsv + 1
                  xvariables(nsv) = sad_v_c(i,n)
               end if
            end do
            if (lsplbkp) then
               ! breakpoints for v
               do n = 5, nsad_v
                  nsv = nsv + 1
                  xvariables(nsv) = sad_v_s(i,n)
               end do
            end if
            if (nsad_u .gt. 0) then
               ! coefficients for u
               do n = 1, nsad_u - 3
                  if (nvar_uc(i,n).ne.0) then
                     nsv = nsv + 1
                     xvariables(nsv) = sad_u_c(i,n)
                  end if
               end do
               if (lsplbkp) then
                  ! breakpoints for u
                  do n = 5, nsad_u
                     nsv = nsv + 1
                     xvariables(nsv) = sad_u_s(i,n)
                  end do
               end if
            end if
         end do

      else

         ! initialize fourier series representation
         do i = 1, nsmid
            do modes = 0,nsad_v
               if( .not. lsad_uv_m(modes)) cycle
               nsv = nsv + 1
               xvariables(nsv) = sad_v_c(i,modes)
               if( modes .ne. 0) then
                  nsv = nsv + 1
                  xvariables(nsv) = sad_v_s(i,modes)
               endif
            end do

            if (nsad_u .le. 0) cycle
            do modes = 0,nsad_u
               if( .not. lsad_uv_m(modes)) cycle
               nsv = nsv + 1
               xvariables(nsv) = sad_u_c(i,modes)
               if( modes .ne. 0) then
                  nsv = nsv + 1
                  xvariables(nsv) = sad_u_s(i,modes)
               endif
            end do
         end do

      end if                 ! if (lspline)
 
      nsad_coeffs = nsv

      if (lsadshape) then
         nvariables = nsv
      else
         nvariables = 0
      end if

      end subroutine init_saddle_coils


      subroutine init_saddle_currents (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use saddle_coils
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, nsv, nc
      integer :: nvariables
      real(rprec) :: xvariables(*)
!-----------------------------------------------
 
      nc = nsad_coils_per_period
      nsad_currents = 0
      nvariables = 0
 
      if (nc .le. 0) return
 
!     initialize the variables to values of unique coil parameters
!     and count the number of variables
 
      nsv = 0
 
      if (lsadcur) then
         do i=1, num_cursad
            if (ls_cur(i)) then
               nsv = nsv + 1
               xvariables(nsv) = cursad(i)
            end if
         end do
      end if
 
      nsad_currents = nsv
      nvariables = nsv

      end subroutine init_saddle_currents


      subroutine load_saddle_coils (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use saddle_coils
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, n, nsv, modes
      integer :: nvariables
      real(rprec) :: xvariables(*)
!-----------------------------------------------
 
!     load the unique coil parameters with values from variables in
!     optimization
 
      nvariables = 0
      nsv = 0

      if (.not.lsadshape) return
 
      if (lspline) then

         ! load bspline series representation
         do i = 1, nsmid
            ! coefficients for v
            do n = 1, nsad_v - 3
               if (nvar_vc(i,n).ne.0) then
                  nsv = nsv + 1
                  sad_v_c(i,n) = xvariables(nsv)
               end if
            end do
            if (lsplbkp) then
               ! breakpoints for v
               do n = 5, nsad_v
                  nsv = nsv + 1
                  sad_v_s(i,n) = xvariables(nsv)
                  if (sad_v_s(i,n) .le. sad_v_s(i,n-1)) then
                     stop 'v-spline knots decreasing'
                  end if
               end do
            end if
            if (nsad_u .gt. 0) then
               ! coefficients for u
               do n = 1, nsad_u - 3
                  if (nvar_uc(i,n).ne.0) then
                     nsv = nsv + 1
                     sad_u_c(i,n) = xvariables(nsv)
                  end if
               end do
               if (lsplbkp) then
                  ! breakpoints for u
                  do n = 5, nsad_u
                     nsv = nsv + 1
                     sad_u_s(i,n) = xvariables(nsv)
                     if (sad_u_s(i,n) .le. sad_u_s(i,n-1)) then
                        stop 'u-spline knots decreasing'
                     end if
                  end do
               end if
            end if
         end do

      else

         ! load fourier series representation
         do i = 1, nsmid
            do modes = 0,nsad_v
               if( .not. lsad_uv_m(modes)) cycle
               nsv = nsv + 1
               sad_v_c(i,modes) = xvariables(nsv)
               if( modes .ne. 0) then
                  nsv = nsv + 1
                  sad_v_s(i,modes) = xvariables(nsv)
               else
                  sad_v_s(i,modes) = 0
               endif
            end do
 
            if (nsad_u .le. 0) cycle
            do modes = 0,nsad_u
               if( .not. lsad_uv_m(modes)) cycle
               nsv = nsv + 1
               sad_u_c(i,modes) = xvariables(nsv)
               if( modes .ne. 0) then
                  nsv = nsv + 1
                  sad_u_s(i,modes) = xvariables(nsv)
               else
                  sad_u_s(i,modes) = 0
               end if
            end do
         end do

      end if                 ! if (lspline)
 
      nvariables = nsv
 
      end subroutine load_saddle_coils


      subroutine load_saddle_currents (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use saddle_coils
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, nsv
      integer :: nvariables
      real(rprec) :: xvariables(*)
!-----------------------------------------------
 
!     Load the coil currents with values from variables in
!     optimization
 
      nsv = 0
 
      if (lsadcur) then
!     Vary saddle currents
         do i = 1, num_cursad
            if (ls_cur(i)) then
               nsv = nsv + 1
               cursad(i) = xvariables(nsv)
            end if
         end do
      end if
 
      nvariables = nsv
 
      end subroutine load_saddle_currents


      subroutine init_tf_coils (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use tor_field
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nvariables, nx_max, ny_max, nx, ny, ntf
      real(rprec) :: xvariables(*)
      real(rprec) :: tfc_xmin, tfc_xmax
      real(rprec) :: tfc_ymin, tfc_ymax
      real(rprec) :: tfc_zmin, tfc_zmax
      real(rprec) :: tfc_dphi, tfc_phi, tfc_rad
      real(rprec) :: tfc_dx, tfc_dy, tfc_scl
!-----------------------------------------------
      mtfcoil = 1
      mtfwire = 2
 
!     Model QPS tf coil with return legs ...
      nx_max = 6
      ny_max = 2
      tfc_dphi = twopi/(nx_max*ny_max)
      tfc_phi = 0.5_dp*tfc_dphi
!     tfc_scl = 0.9_dp ! since the following are engineering values
      tfc_scl = 1.0_dp ! for vmec input based on actual plasma size
      tfc_rad = 2.15_dp/tfc_scl
      tfc_xmin = -0.375_dp/tfc_scl
      tfc_xmax =  0.375_dp/tfc_scl
      tfc_dx = 0.150_dp/tfc_scl
      tfc_ymax = 0.013_dp/tfc_scl
      tfc_dy = 0.026_dp/tfc_scl
      tfc_zmin = -2.48_dp/tfc_scl
      tfc_zmax =  2.48_dp/tfc_scl
      if (lqos) then
         mtfwire = 5
!        distribute straight tf filaments in a rectangle in x-y plane
         mtfcoil = nx_max*ny_max
         ntf = 0
         do ny = 1, ny_max
            do nx = 1, nx_max
               ntf = ntf + 1
               if (ny .eq. 1) then
                  tfc_x(ntf,1) = (tfc_xmax - (nx-1)*tfc_dx)
               else
                  tfc_x(ntf,1) = (tfc_xmin + (nx-1)*tfc_dx)
               end if
               tfc_x(ntf,2) = tfc_x(ntf,1)
               tfc_x(ntf,3) = tfc_rad*cos(tfc_phi)
               tfc_x(ntf,4) = tfc_x(ntf,3)
               tfc_x(ntf,5) = tfc_x(ntf,1)
               tfc_y(ntf,1) = (tfc_ymax - (ny-1)*tfc_dy)
               tfc_y(ntf,2) = tfc_y(ntf,1)
               tfc_y(ntf,3) = tfc_rad*sin(tfc_phi)
               tfc_y(ntf,4) = tfc_y(ntf,3)
               tfc_y(ntf,5) = tfc_y(ntf,1)
               tfc_z(ntf,1) = tfc_zmax
               tfc_z(ntf,2) = tfc_zmin
               tfc_z(ntf,3) = tfc_zmin
               tfc_z(ntf,4) = tfc_zmax
               tfc_z(ntf,5) = tfc_zmax
               tfc_cur(ntf) = i_tfc/mtfcoil
               tfc_phi = tfc_phi + tfc_dphi
            end do
         end do
      else
!        single filament at x=0, y=0
         tfc_zmin = -1000
         tfc_zmax =  1000
         tfc_x(1,1) = 0
         tfc_x(1,2) = 0
         tfc_y(1,1) = 0
         tfc_y(1,2) = 0
         tfc_z(1,1) = tfc_zmax
         tfc_z(1,2) = tfc_zmin
         tfc_cur(1) = i_tfc
      end if
 
      xvariables (1) = i_tfc
      nvariables = 1
 
      end subroutine init_tf_coils


      subroutine load_tf_coils (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use tor_field
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, nvariables
      real(rprec) :: xvariables(*)
!-----------------------------------------------
      nvariables = 1
      i_tfc = xvariables(1)
      tfc_cur(1) = xvariables(1)
 
      if (lqos) then
         do i = 1, mtfcoil
            tfc_cur(i) = xvariables(1)/mtfcoil
         end do
      end if
 
      end subroutine load_tf_coils


      subroutine init_vf_coils (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vf_coils
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, k, nsv
      integer :: nvariables
      real(rprec) :: xvariables(*)
!-----------------------------------------------

      nsv = 0

      if (lvfvar) then
         do i=1, num_vf
            if (lvfz) then
               nsv = nsv + 1
               xvariables(nsv) = zc_vf(i)
            end if
            if (lvfr) then
               nsv = nsv + 1
               xvariables(nsv) = rc_vf(i)
            end if
            if (nrvf_c .gt. 0) then
               do k=1, nrvf_c
                  nsv = nsv + 1
                  xvariables(nsv) = rcfc_vf(i,k)
                  nsv = nsv + 1
                  xvariables(nsv) = rcfs_vf(i,k)
               end do
            end if
         end do
      end if

      nvariables = nsv
      nvf_coeffs = nsv

      end subroutine init_vf_coils


      subroutine load_vf_coils (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vf_coils
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, k, nsv
      integer :: nvariables
      real(rprec) :: xvariables(*)
!-----------------------------------------------

      nsv = 0

      if (lvfvar) then
         do i=1, num_vf
            if (lvfz) then
               nsv = nsv + 1
               zc_vf(i) = xvariables(nsv)
            end if
            if (lvfr) then
               nsv = nsv + 1
               rc_vf(i) = xvariables(nsv)
            end if
            if (nrvf_c .gt. 0) then
               do k=1, nrvf_c
                  nsv = nsv + 1
                  rcfc_vf(i,k) = xvariables(nsv)
                  nsv = nsv + 1
                  rcfs_vf(i,k) = xvariables(nsv)
               end do
            end if
         end do
      end if

      nvariables = nsv

      end subroutine load_vf_coils


      subroutine init_vf_currents (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vf_coils
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      integer :: nvariables
      real(rprec) :: xvariables(*)
!-----------------------------------------------
      nvariables = 0

      if (lvfc) then
        do i=1, num_vf
          if (lcc_vf(i)) then
             nvariables = nvariables + 1
             xvariables(nvariables) = cc_vf(i)
          endif
        end do
      end if

      end subroutine init_vf_currents


      subroutine load_vf_currents (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vf_coils
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      integer :: nvariables
      real(rprec) :: xvariables(*)
!-----------------------------------------------
!     Load the coil currents with values from optimization variables
 
      nvariables = 0

      if (lvfc) then
        do i=1, num_vf
          if (lcc_vf(i)) then
             nvariables = nvariables + 1
             cc_vf(i) = xvariables(nvariables)
          endif
        end do
      end if
 
      end subroutine load_vf_currents


      subroutine init_bg_currents (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use bcoils_mod
      use safe_open_mod
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, m, iunit, istat
      integer, intent(out) :: nvariables
      real(rprec) :: xvariables(*)
      character*200 :: temp
!-----------------------------------------------
      mbcoils = 0
      iunit = 77
      call safe_open(iunit, istat, trim(bcoil_file), 'old',
     1                       'formatted')
      if (istat .eq. 0) then
         read (iunit, *) mbcoils
         close (iunit)
      else
         print *, 'Background file - ', trim(bcoil_file),
     1            ' - could not be opened'
      end if
 
!     Set variable background currents based on index mc_bg
 
      mc_max = 0
      do i = 1, mbcoils
         m = mc_bg(i)
         if (m .gt. mc_max) then
            cc_bg(m) = bcoil_cur(i)
            mc_max = m
         end if
      end do

!
!     lbcoil_cur == .false. is equivalent to mc_max = 0
!
      if (mc_max == 0) lbcoil_cur = .false. 

      if (lbcoil_cur) then
         do m = 1, mc_max
            xvariables(m) = cc_bg(m)
         end do
      end if
 
      nvariables = mc_max

      end subroutine init_bg_currents


      subroutine load_bg_currents (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use bcoils_mod
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, m
      integer :: nvariables
      real(rprec) :: xvariables(*)
!-----------------------------------------------
 
!     Load variable background currents cc_bg into bcoil_cur
!     based on index mc_bg
 
      if (lbcoil_cur) then

         do m = 1, mc_max
            cc_bg(m) = xvariables(m)
         end do
 
         do i = 1, mbcoils
            m = mc_bg(i)
            if (m .gt. 0) then
               bcoil_cur(i) = cc_bg(m)
            end if
         end do

         nvariables = mc_max
      
      else

         nvariables = 0

      end if

      end subroutine load_bg_currents


      subroutine init_modular_wsurf (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use modular_coils
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nvariables
      real(rprec) :: xvariables(*)
!-----------------------------------------------
      xvariables (1:numsurf) = rmn_sf (1:numsurf)
      xvariables (numsurf+1:2*numsurf) = zmn_sf (1:numsurf)
      nvariables = 2*numsurf
 
      end subroutine init_modular_wsurf


      subroutine load_modular_wsurf (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use modular_coils
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nvariables
      real(rprec) :: xvariables(*)
!-----------------------------------------------
      rmn_sf (1:numsurf) = xvariables (1:numsurf)
      zmn_sf (1:numsurf) = xvariables (numsurf+1:2*numsurf)
      nvariables = 2*numsurf
 
      end subroutine load_modular_wsurf


      subroutine init_saddle_wsurf (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use coilsnamin, ONLY : nopt_wsurf
      use saddle_surface
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nb, mb, ik
      real(rprec), allocatable, dimension(:,:) :: rbc_in, zbs_in
      real(rprec) :: delta
      integer :: nvariables
      real(rprec) :: xvariables(*)
      external unique_boundary, convert_boundary, 
     1   unique_boundary_PG, convert_boundary_PG
!-----------------------------------------------

      if (nopt_wsurf .eq. -1) then
         xvariables(1:numsurf_sad) = rmn_sad(:numsurf_sad)
         xvariables(numsurf_sad+1:2*numsurf_sad) = zmn_sad(:numsurf_sad)
         nvariables = 2*numsurf_sad
         return
      end if

      irm0_bdy = 0;      izm0_bdy = 0;      irho_bdy = 0;

      ntor_opt = max(maxval(n_sad), abs(minval(n_sad)))
      mpol_opt = max(maxval(m_sad), abs(minval(m_sad)))
      if (nopt_wsurf .eq. 1) mpol_opt = mpol_opt + 1
      ik = (2*ntor_opt+1)*(mpol_opt+1)
      
      allocate (nbrho_opt(ik), mbrho_opt(ik), 
     1  rbc(-ntor_opt:ntor_opt,0:mpol_opt), 
     2  zbs(-ntor_opt:ntor_opt,0:mpol_opt),
     3  rbc_in(-ntor_opt:ntor_opt,0:mpol_opt), 
     4  zbs_in(-ntor_opt:ntor_opt,0:mpol_opt),
     5  rhobc(-ntor_opt:ntor_opt,0:mpol_opt), 
     6  nrz0_opt(2*ntor_opt),
     7  delta_mn(-ntor_opt:ntor_opt,-mpol_opt:mpol_opt))
      
      rbc = 0
      zbs = 0
      nvariables = 0
      do ik = 1, numsurf_sad
         mb = m_sad(ik)
 
         if (mb .eq. 0) then
            nb = n_sad(ik)
            if (nb .ge. 0) then
               rbc(nb, mb) = rmn_sad(ik)
               zbs(nb, mb) = -zmn_sad(ik)
            else
               nb = -nb
               rbc(nb, mb) = rmn_sad(ik)
               zbs(nb, mb) = zmn_sad(ik)
            end if                              ! mb=0

         else
            nb = -n_sad(ik)
            rbc(nb, mb) = rmn_sad(ik)
            zbs(nb, mb) = zmn_sad(ik)
         end if

      end do


      if (nopt_wsurf .eq. 0) then

         rbc_in = rbc
         zbs_in = zbs                                                    !copy and store m=0 components

!
!        check if conversion was made or if original bdy already in proper form
!        IF NOPT_BOUNDARY=0, USE HIRSHMAN/BRESLAU REPRESENTATION
!                        =1, USE PG (PAUL GARABEDIAN) DELTA_MN REPRESENTATION
!
            call convert_boundary(rbc, zbs, rhobc, mpol_opt, ntor_opt)

               call unique_boundary(rbc_in, zbs_in, rhobc, ntor_opt,
     1                        mpol_opt, mpol_opt, ntor_opt)
               delta = sum((rbc - rbc_in)**2)/rbc(0,1)**2
     1               + sum((zbs - zbs_in)**2)/zbs(0,1)**2

               if (delta .gt. 1.e-8_dp) write(*,10) 100*(one-delta)

 10   format(' Input boundary representation was converted!',/,
     1       ' Reliability of conversion = ',f7.2,'%')

         do nb = -ntor_opt, ntor_opt
            if (rbc(nb,0).ne.zero .and. (nb.ne.0)) then   !mb=0, nb=0 not varied
               nvariables = nvariables + 1
               irm0_bdy = irm0_bdy + 1
               nrz0_opt(irm0_bdy) = nb
               xvariables(nvariables) = rbc(nb,0)
            endif
         end do

         do nb = -ntor_opt, ntor_opt
            if (zbs(nb,0).ne.zero .and. (nb.ne.0)) then
               nvariables = nvariables + 1
               izm0_bdy = izm0_bdy + 1
               nrz0_opt(nvariables) = nb
               xvariables(nvariables) = zbs(nb,0)
            endif
         end do

         do nb = -ntor_opt, ntor_opt
            do mb = 0, mpol_opt
               if (rhobc(nb,mb) .ne. zero 
     1             .and. (mb .ne. 0 .or. nb .ge. 0)) then
                  nvariables = nvariables + 1
                  irho_bdy = irho_bdy + 1
                  nbrho_opt(irho_bdy) = nb
                  mbrho_opt(irho_bdy) = mb
                  xvariables(nvariables) = rhobc(nb,mb)
               endif
            end do
         end do

         end if

         if (nopt_wsurf.eq.1) then
           call convert_boundary_PG(rbc,zbs,delta_mn,mpol_opt,ntor_opt)

           nvariables = 0
           do nb = -ntor_opt, ntor_opt
             do mb = -mpol_opt, mpol_opt
               if (delta_mn(nb,mb) .ne. zero 
     1            .and. .not.(nb .eq.0 .and. mb .eq. 0)) then
                  irho_bdy = irho_bdy + 1
                  nvariables = nvariables + 1
                  nbrho_opt(irho_bdy) = nb
                  mbrho_opt(irho_bdy) = mb
                  xvariables(nvariables) = delta_mn(nb,mb)
               endif
             end do
           end do
          
         end if

      deallocate (rbc_in, zbs_in)      

      end subroutine init_saddle_wsurf


      subroutine load_saddle_wsurf (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use saddle_surface
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nvariables, mb, nb, ik, j1
      real(rprec) :: xvariables(*)
      external unique_boundary, unique_boundary_PG
!-----------------------------------------------

      if (nopt_wsurf .eq. -1) then
         rmn_sad(:numsurf_sad) = xvariables(1:numsurf_sad)
         zmn_sad(:numsurf_sad) = xvariables(numsurf_sad+1:2*numsurf_sad)
         nvariables = 2*numsurf_sad
         return
      end if

      if (nopt_wsurf .eq. 0) then
         do ik = 1, irm0_bdy
            j1 = ik
            rbc(nrz0_opt(ik),0) = xvariables(j1)
         end do
         do ik = 1+irm0_bdy, irm0_bdy+izm0_bdy
            j1 = ik
            zbs(nrz0_opt(ik),0) = xvariables(j1)
         end do
         do ik = 1, irho_bdy
            j1 = ik + irm0_bdy + izm0_bdy
            rhobc(nbrho_opt(ik),mbrho_opt(ik)) = xvariables(j1)
         end do

         call unique_boundary(rbc, zbs, rhobc, ntor_opt, mpol_opt,
     2                            mpol_opt, ntor_opt)
         nvariables = 0
         do mb = 0, mpol_opt
            do nb = -ntor_opt, ntor_opt
               nvariables = nvariables + 1
               m_sad(nvariables) = mb
               if (mb .eq. 0) then
                  n_sad(nvariables) = nb
                  rmn_sad(nvariables) = rbc(nb,mb)
                  zmn_sad(nvariables) = -zbs(nb,mb)
               else
                  n_sad(nvariables) = -nb
                  rmn_sad(nvariables) = rbc(nb,mb)
                  zmn_sad(nvariables) = zbs(nb,mb)
               end if
            end do
         end do

         nvariables = irm0_bdy + izm0_bdy + irho_bdy
      end if

      if (nopt_wsurf .eq. 1) then
         do ik = 1, irho_bdy
            j1 = ik
            delta_mn(nbrho_opt(ik),mbrho_opt(ik)) = xvariables(j1)
         end do

         call unique_boundary_PG(rbc, zbs, delta_mn, ntor_opt, mpol_opt,
     1                           mpol_opt, ntor_opt)

         nvariables = 0
         do mb = 0, mpol_opt
            do nb = -ntor_opt, ntor_opt
               if( abs(rbc(nb,mb)) .gt. zero .or.
     1            abs(zbs(nb,mb)) .gt. zero ) then
                  nvariables = nvariables + 1
                  m_sad(nvariables) = mb
                  if (mb .eq. 0) then
                     n_sad(nvariables) = nb
                     rmn_sad(nvariables) = rbc(nb,mb)
                     zmn_sad(nvariables) = -zbs(nb,mb)
                  else
                     n_sad(nvariables) = -nb
                     rmn_sad(nvariables) = rbc(nb,mb)
                     zmn_sad(nvariables) = zbs(nb,mb)
                  end if
               end if
            end do
         end do

         if (nvariables .ne. numsurf_sad) then
            write(6,*) "nvariables != numsurf_sad in load_saddle_surf"
            stop
         end if

         nvariables = irho_bdy

      end if

      end subroutine load_saddle_wsurf


      subroutine write_coilsin (iunit, istat)
      use modular_coils
      use saddle_coils
      use saddle_surface
      use tor_field
      use coilsnamin
      use bcoils_mod
      use cparms
      use control_mod
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer :: iunit, istat
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i,j,ierr, iter
      real(rprec) :: tempc(nmid), temps(nmid)
      real(rprec) :: temp2c(nsmid), temp2s(nsmid)
      character*(*), parameter :: fmt(4) = ( /
     1  '(A8,I1,A1) ','(A8,I2,A1) ','(A11,I1,A1)','(A11,I2,A1)' / )
      character :: cfmt*15, temp_write*100
!-----------------------------------------------
      istat = 0
      do iter = 2, 100
         if(niter_array(iter) .le. 0 .or. 
     1      epsfcn_array(iter) .eq. 0) exit
      enddo
      iter = min(iter, 100)

!     Write out two-dimensional arrays with indices explicitly
!     so that if NFPDIM changes, the namelist can still be read
!     correctly
 
      write(iunit,'(1x,a8)')'&COILSIN'
      write(iunit,'(1x,3a)') "BCOIL_FILE = '", trim(bcoil_file),"'"
      write(iunit,100) 'LRESTART = ', lrestart
      write(iunit,100) 'LSURFV = ', lsurfv
      write(iunit,100) 'LSADSFV = ', lsadsfv
      write(iunit,100) 'LVF = ', lvf
      write(iunit,100) 'LVFC = ', lvfc
      write(iunit,100) 'LVFVAR = ', lvfvar
      write(iunit,100) 'LVFR = ', lvfr
      write(iunit,100) 'LVFZ = ', lvfz
      write(iunit,100) 'LTFC = ', ltfc
      write(iunit,100) 'LTFCV = ', ltfcv
      write(iunit,100) 'LSADDLE = ', lsaddle
      write(iunit,100) 'LSMOD = ', lsmod
      write(iunit,100) 'LSPLINE = ', lspline
      write(iunit,100) 'LSPLBKP = ', lsplbkp
      write(iunit,100) 'LMODULAR = ', lmodular
      write(iunit,100) 'LMODCUR = ', lmodcur
      write(iunit,100) 'LSADCUR = ', lsadcur
      write(iunit,100) 'LSADSHAPE = ', lsadshape
      write(iunit,100) 'LPOLCUR = ', lpolcur
      write(iunit,100) 'LBNORM = ', lbnorm
      write(iunit,100) 'LBCOIL = ', lbcoil
      write(iunit,100) 'LBCOIL_CUR = ', lbcoil_cur
      write(iunit,100) 'LNCSX = ', lncsx
      write(iunit,100) 'LQOS = ', lqos
      write(iunit,100) 'LSYMM = ', lsymm
      write(iunit,100) 'LACCESS = ', laccess
      write(iunit,200) 'NOPT_ALG = ', nopt_alg 
      write(iunit,200) 'NOPT_WSURF = ', nopt_wsurf 
      write(iunit,200) 'NITER_OPT = ',niter_opt, 'NSTEP = ', nstep
      write(iunit,201) 'NITER_ARRAY = ',(niter_array(i),i=1,iter)
      write(iunit,200) 'NF_PHI = ',nf_phi
      if (nf_rho .gt. 0) then
         write(iunit,200) 'NF_RHO = ',nf_rho
      end if
      write(iunit,300) 'EPSFCN = ', epsfcn
      call write_array(iunit,'EPSFCN_ARRAY', epsfcn_array, iter)
      write(iunit,300) 'I_POL = ', i_pol
      write(iunit,300) 'I_TFC = ', i_tfc
      write(iunit,200) 'NMOD_COILS_PER_PERIOD = ',
     1                  nmod_coils_per_period
      write(iunit,200) 'NUM_VF = ', num_vf
!     Use following weights, etc., for both modular and saddle minimum
!     coil-plasma distance penalties
      write(iunit,300) 'DCP_WGT = ', dcp_wgt, 'DCP_EXP = ', dcp_exp,
     1                 'DCP_TGT = ', dcp_tgt
      write(iunit,450) 'DPC_WGT = ', dpc_wgt
      write(iunit,450) 'MXB_WGT = ', mxb_wgt
      if (nmod_coils_per_period .gt. 0) then
         call write_array(iunit,'DCC_WGT', dcc_wgt, nmid)
         call write_array(iunit,'DCC_EXP', dcc_exp, nmid)
         call write_array(iunit,'DCC_TGT', dcc_tgt, nmid)
         write(iunit,'(a)') ' LMOD_WGT = '
         write(iunit,400) (lmod_wgt(i),i=1,nmid)
         write(iunit,'(a)') ' LMOD_TGT = '
         write(iunit,400) (lmod_tgt(i),i=1,nmid)
         call write_array(iunit,'RC_WGT', rc_wgt, nmid)
         call write_array(iunit,'RC_EXP', rc_exp, nmid)
         call write_array(iunit,'RC_TGT', rc_tgt, nmid)
         call write_array(iunit,'CU_WGT', cu_wgt, nmid)
         call write_array(iunit,'CU_TGT', cu_tgt, nmid)
      end if
      if (lmodular .and. (.not.lsaddle)) then
         call write_array(iunit,'YMIN_WGT', ymin_wgt, nmid)
         call write_array(iunit,'YMIN_TGT', ymin_tgt, nmid)
      end if
      if (lsaddle) then
         call write_array(iunit,'YMIN_WGT', ymin_wgt, nsmid)
         call write_array(iunit,'YMIN_TGT', ymin_tgt, nsmid)
      end if
      if (lncsx) then
         call write_array(iunit,'R_EXT', r_ext, nmid)
      end if
      if (laccess) then
         write(iunit,200) 'N_ACCESS = ', n_access
         call write_array(iunit,'X0_ACCESS', x0_access, n_access)
         call write_array(iunit,'Y0_ACCESS', y0_access, n_access)
         call write_array(iunit,'Z0_ACCESS', z0_access, n_access)
         call write_array(iunit,'X1_ACCESS', x1_access, n_access)
         call write_array(iunit,'Y1_ACCESS', y1_access, n_access)
         call write_array(iunit,'Z1_ACCESS', z1_access, n_access)
         call write_array(iunit,'DAC_WGT', dac_wgt, n_access)
         call write_array(iunit,'DAC_EXP', dac_exp, n_access)
         call write_array(iunit,'DAC_TGT', dac_tgt, n_access)
      end if
      if (lbcoil) then
         write(iunit,'(a)') ' MC_BG = '
         write(iunit,140) (mc_bg(i), i=1,mbcoils)
         write(iunit,'(a)') ' LP_BG = '
         write(iunit,150) (lp_bg(i), i=1,mbcoils)
         call write_array(iunit,'BCOIL_CUR', bcoil_cur, mbcoils)
      end if
      if (num_vf .gt. 0) then
         write(iunit,'(a)') ' LCC_VF = '
         write(iunit,150) (lcc_vf(i), i=1,num_vf)
         call write_array(iunit,'CC_VF', cc_vf, num_vf)
         call write_array(iunit,'RC_VF', rc_vf, num_vf)
         call write_array(iunit,'ZC_VF', zc_vf, num_vf)
         call write_array(iunit,'CVF_WGT', cvf_wgt, num_vf)
         call write_array(iunit,'CVF_TGT', cvf_tgt, num_vf)
         call write_array(iunit,'RVF_WGT', rvf_wgt, num_vf)
         call write_array(iunit,'RVF_TGT', rvf_tgt, num_vf)
         write(iunit,200) 'NRVF_C = ', nrvf_c
         if (nrvf_c .gt. 0) then
            do j = 1, nrvf_c
              cfmt = fmt(3)
              if( j.gt.9 ) cfmt = fmt(4)
              do i = 1, num_vf
                temp2c(i) = rcfc_vf(i,j)
                if (abs(temp2c(i)) .lt. 1.e-10) temp2c(i) = 0
                temp2s(i) = rcfs_vf(i,j)
                if (abs(temp2s(i)) .lt. 1.e-10) temp2s(i) = 0
              end do
              write(temp_write, cfmt) 'RCFC_VF(1:,',J,')'
              call write_array(iunit,trim(temp_write), temp2c, num_vf)
              write(temp_write, cfmt) 'RCFS_VF(1:,',J,')'
              call write_array(iunit,trim(temp_write), temp2s, num_vf)
            end do
         end if
      end if
      if (nmod_coils_per_period .gt. 0) then
         call write_array(iunit,'CURMOD', curmod, nmid)
      end if
      if (nf_phi .gt. 0) then
         do j = 0,nf_phi
           cfmt = fmt(1)
           if( j.gt.9 ) cfmt = fmt(2)
           do i = 1,nmid
             tempc(i) = phic(i,j)                        !modular(i)%phic(j)
             temps(i) = phis(i,j)                        !modular(i)%phis(j)
             if (abs(tempc(i)) .lt. 1.e-10) tempc(i) = 0
             if (abs(temps(i)) .lt. 1.e-10) temps(i) = 0
           end do
           write(temp_write, cfmt) 'PHIC(1:,',J,')'
           call write_array(iunit,trim(temp_write), tempc, nmid)
           write(temp_write, cfmt) 'PHIS(1:,',J,')'
           call write_array(iunit,trim(temp_write), temps, nmid)
         end do
      end if
      if (nf_rho .gt. 0) then
         do j = 0,nf_rho
           cfmt = fmt(1)
           if( j.gt.9 ) cfmt = fmt(2)
           do i = 1,nmid
             temps(i) = rhos(i,j)                        !modular(i)%rhos(j)
             if (abs(temps(i)) .lt. 1.e-10) temps(i) = 0
           end do
           write(temp_write, cfmt) 'RHOS(1:,',J,')'
           call write_array(iunit,trim(temp_write), rhos, nmid)
         end do
      end if
 
! 2/20/98 WHM Write out surface coefficients
  
      if (numsurf .gt. 0) then
         write(iunit,200)'numsurf = ',numsurf
         write(iunit,'(a)')' m_num = '
         write(iunit,210)(m_num(i),i=1,numsurf)
         write(iunit,'(a)')' n_num = '
         write(iunit,210)(n_num(i),i=1,numsurf)
         call write_array(iunit,'rmn_sf', rmn_sf, numsurf)
         call write_array(iunit,'zmn_sf', zmn_sf, numsurf)
      end if
 
!     Saddle coil input parameters
 
      write(iunit,200) 'NSAD_COILS_PER_PERIOD = ',nsad_coils_per_period
      write(iunit,200) 'NSAD_U = ',nsad_u
      write(iunit,200) 'NSAD_V = ',nsad_v
      if (nsad_coils_per_period .gt. 0) then
         write(iunit,200) 'NFILS = ',nfils
         write(iunit,300) 'DELN = ', deln, 'DELT = ', delt
         write(iunit,'(a)') ' NSAD_GROUP = '
         write(iunit,140) (nsad_group(i), i=1,nsmid)
         write(iunit,'(a)') ' LS_CUR = '
         write(iunit,150) (ls_cur(i), i=1,nsmid)
         call write_array(iunit,'CURSAD', cursad, nsmid)
         call write_array(iunit,'CSC_WGT', csc_wgt, nsmid)
         call write_array(iunit,'CSC_TGT', csc_tgt, nsmid)
         call write_array(iunit,'CSAD_SCL', csad_scl, nsmid)
         call write_array(iunit,'DSC_WGT', dsc_wgt, nsmid)
         call write_array(iunit,'DSC_EXP', dsc_exp, nsmid)
         call write_array(iunit,'DSC_TGT', dsc_tgt, nsmid)
         call write_array(iunit,'RS_WGT', rs_wgt, nsmid)
         call write_array(iunit,'RS_EXP', rs_exp, nsmid)
         call write_array(iunit,'RS_TGT', rs_tgt, nsmid)
         call write_array(iunit,'CS_WGT', cs_wgt, nsmid)
         call write_array(iunit,'CS_TGT', cs_tgt, nsmid)
         call write_array(iunit,'LSAD_WGT', lsad_wgt, nsmid)
         call write_array(iunit,'LSAD_TGT', lsad_tgt, nsmid)
         call write_array(iunit,'RMAX_WGT', rmax_wgt, nsmid)
         call write_array(iunit,'RMAX_TGT', rmax_tgt, nsmid)
         call write_array(iunit,'SAD_V0', sad_v0, nsmid)
         call write_array(iunit,'SAD_U0', sad_u0, nsmid)
         call write_array(iunit,'DSCXP_WGT', dscxp_wgt, nsmid)
         call write_array(iunit,'DSCXP_EXP', dscxp_exp, nsmid)
         call write_array(iunit,'DSCXP_TGT', dscxp_tgt, nsmid)
         call write_array(iunit,'SCD_WGT', scd_wgt, nsmid)
         call write_array(iunit,'SCD_TGT', scd_tgt, nsmid)
         if (lspline .and. lsplbkp) then
            write(iunit,350) 'BKP_WGT = ',bkp_wgt,'BKP_TGT = ',bkp_tgt
         end if
      end if
      if( max(nsad_v, nsad_u) .gt. 0) then
         j = max(nsad_v, nsad_u)
         call write_larray(iunit,'LSAD_UV_M', lsad_uv_m(0:j), j+1)
      endif

      if (nsad_v .gt. 0) then
         do j = 0,nsad_v
           cfmt = fmt(3)
           if( j.gt.9 ) cfmt = fmt(4)
           do i = 1,nsmid
             temp2c(i) = sad_v_c(i,j)                    !saddle(i)%v_c(j)
             temp2s(i) = sad_v_s(i,j)                    !saddle(i)%v_s(j)
             if (abs(temp2c(i)) .lt. 1.e-10) temp2c(i) = 0
             if (abs(temp2s(i)) .lt. 1.e-10) temp2s(i) = 0
           end do
           write(temp_write, cfmt) 'SAD_V_C(1:,',J,')'
           call write_array(iunit,trim(temp_write), temp2c, nsmid)
           if (lspline) then
              write(temp_write, cfmt) 'NVAR_VC(1:,',J,')'
              write(iunit,160) trim(temp_write)
              write(iunit,140) (nvar_vc(i,j), i=1,nsmid)
           end if
           write(temp_write, cfmt) 'SAD_V_S(1:,',J,')'
           call write_array(iunit,trim(temp_write), temp2s, nsmid)
         end do
      end if
      if (nsad_u .gt. 0) then
         do j = 0,nsad_u
           cfmt = fmt(3)
           if( j.gt.9 ) cfmt = fmt(4)
           do i = 1,nsmid
             temp2c(i) = sad_u_c(i,j)                   !saddle(i)%u_c(j)
             temp2s(i) = sad_u_s(i,j)                   !saddle(i)%u_s(j)
             if (abs(temp2c(i)) .lt. 1.e-10) temp2c(i) = 0
             if (abs(temp2s(i)) .lt. 1.e-10) temp2s(i) = 0
           end do
           write(temp_write, cfmt) 'SAD_U_C(1:,',J,')'
           call write_array(iunit,trim(temp_write), temp2c, nsmid)
           if (lspline) then
              write(temp_write, cfmt) 'NVAR_UC(1:,',J,')'
              write(iunit,160) trim(temp_write)
              write(iunit,140) (nvar_uc(i,j), i=1,nsmid)
           end if
           write(temp_write, cfmt) 'SAD_U_S(1:,',J,')'
           call write_array(iunit,trim(temp_write), temp2s, nsmid)
         end do
      end if
 
!     Saddle surface coefficients
  
      if (numsurf_sad .gt. 0) then
         write(iunit,200)'numsurf_sad = ', numsurf_sad
         write(iunit,'(a)')' m_sad = '
         write(iunit,210)(m_sad(i),i=1,numsurf_sad)
         write(iunit,'(a)')' n_sad = '
         write(iunit,210)(n_sad(i),i=1,numsurf_sad)
         call write_array(iunit,'rmn_sad', rmn_sad, numsurf_sad)
         call write_array(iunit,'zmn_sad', zmn_sad, numsurf_sad)
      end if
 
      write(iunit,'(a)')'/'
 
 100  format(4(1x,a,l2,','))
 140  format(10(i2,','))
 150  format(10(l2,','))
 160  format(2x,a,' = ')
 200  format(4(1x,a,i6,','))
 201  format(1x,a,(' ', 10i7))
 210  format(10(1x,i4,','))
 300  format(3(1x,a,1pe12.4,','))
 350  format(2(1x,a,1pe12.4,','))
 400  format(3(1pe15.8,','))
 450  format(1x,a,1pe12.4,',')
 1000 format(3(1pe25.18,','))
 
  
      end subroutine write_coilsin
