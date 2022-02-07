      module booz_params
      use kind_spec
      implicit none
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: unit_booz = 20
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, dimension(0:1) :: ntorsum
      integer :: nfp, ns, mpol, mpol1, ntor, mnmax
      integer :: nu_boz, nv_boz, nunv, mboz, nboz, mnboz
      real(rprec), dimension(:), allocatable ::
     1  hiota, phip, gpsi, ipsi
      real(rprec), dimension(:), allocatable :: pmn, bmnb,
     1   rmnb, zmnb, pmnb, gmnb
      real(rprec), dimension(:,:), allocatable :: bmod_b
      real(rprec) :: ohs
      logical :: lscreen
      end module booz_params
          
      module booz_extern
      use kind_spec
C     VARIABLES NEEDED OUTSIDE THE BOOZER-COORDINATE CODE ARE USED HERE
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec), dimension(:,:), allocatable ::
     1   bsubtmn, bsubzmn
c-----------------------------------------------
      end module booz_extern
 
      module booz_persistent
      use kind_spec
C     ONLY ARRAYS SMALL ENOUGH TO SAVE BETWEEN INTERNAL CALLS SHOULD BE STORED HERE
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nu2_b
      integer, dimension(:), allocatable :: ixm, ixn
      real(rprec), dimension(:), allocatable :: 
     1   sfull, scl, thgrd, ztgrd, xm, xn, xnb, xmb
      real(rprec), dimension(:,:), allocatable ::
     1   cosm_b, cosn_b, sinm_b, sinn_b, rmnc, zmns, lmns 
C-----------------------------------------------
      end module booz_persistent
