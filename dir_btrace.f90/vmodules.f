      module btrace_mod
        use kind_spec
!**********************************************************************
!
!  namelist for controlling btrace
!
!
!  nfp  number of field periods.
!       0 assume no toroidal periodicity
!
!  lsym T  assume stellarator symmetry
!       F  assume no symmetry
!
!  field_source
!       0 Biot-Savart from coils-file
!       1 3D cubic splined from grid-file
!       2 dynamically generated splined-grid using Biot-Savart from coils
!
!  bspline
!       F field is integrated directly
!       T field is integrated on a grid, fit to bspline as a map.  
!         bspline coeffs are interpolated to subsequent tracing
!
!  rmin, rmax, zmin, zmax
!       bounds of computational domain
!
!  ntransits
!       number of period transits (or toroidal transits for nfp=0) to trace
!
!  phi_pl
!       toroidal angle (degrees) of poloidal planes to make plots at
!
!  step
!       horizontal step size in initial positions
!**********************************************************************
        integer, parameter :: max_planes = 100
        integer, parameter :: max_currents = 100

        integer :: nfp, ntransits, nplanes, field_source, 
     1             n_starts, max_starts, n_evals
        real(rprec) :: rmax, rmin, zmax, zmin, step, tol
        real(rprec), dimension(max_planes) :: phi_pl
        real(rprec), dimension(max_planes) :: phi_poin
        real(rprec), dimension(max_currents) :: coil_current
        character*(100) :: coils_file, mgrid_file 
        logical :: bspline, lsym

        real(rprec), dimension(:,:,:), allocatable :: r_poin,z_poin
        
        namelist /btrnml/ nfp, ntransits, field_source, 
     1    coils_file, mgrid_file, lsym, coil_current, tol,
     1    rmax, rmin, zmax, zmin, step, phi_pl, bspline, max_starts

      end module btrace_mod

!**********************************************************************
      module mgrid_mod
      use kind_spec 
!
!     grid dimensions:
!              ir  = no. radial (r) points in box
!              jz  = no. z points in box
!
!              suggest choosing hr == (rmax-rmin)/(ir-1) equal to
!                               hz == (zmax-zmin)/(jz-1)
!
      integer, parameter :: maxgroups = 100
      integer, parameter :: ir = 101, jz = 101
!                          ,kp = 16
!     real(rprec), parameter :: rmin = 0.09, rmax = 2.51,
!    1                          zmin =-1.21, zmax = 1.21
      integer :: nspul(maxgroups), ngroup(maxgroups)
      integer :: unit_parsed = 200
      integer :: nextcur
      real(rprec) :: extcur(maxgroups)
      character*8 :: dsilabel, bloopnames
      end module mgrid_mod

      module bcoils_mod
      use kind_spec
      integer :: nall, isteu
      real(rprec), dimension(:), allocatable :: xw, yw, zw, curre
      real(rprec) :: curcoils
      real(rprec), dimension(:), allocatable :: vx, vy, vz, dx, dy, dz
      end module bcoils_mod
