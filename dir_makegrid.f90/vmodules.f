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
      integer :: ir = 101, jz = 101

      integer :: nspul(maxgroups), ngroup(maxgroups)
      integer :: unit_parsed = 200
      real(rprec) :: extcur(maxgroups)
      character*8 :: dsilabel, bloopnames
      end module mgrid_mod

      module bpol_mod
      use kind_spec
      integer :: nall, isteu
      real(rprec), dimension(:), allocatable :: xw, yw, zw, curre
      real(rprec) :: curcoils
      end module bpol_mod

      module becoil1_mod
      use kind_spec
      real(rprec), dimension(:), allocatable :: vx, vy, vz, dx, dy, dz
      end module becoil1_mod
