      module NumParams
      use kind_spec
      integer, parameter :: inescoil = 6
      integer :: inesc = inescoil
      real(rprec), parameter :: zero = 0,
     1    one = 1, two = 2, three = 3,
     2    pi = 3.14159265358979312_dp, pi2 = two*pi, mu0 = 4.0e-7_dp*pi

      contains
        logical function is_zero(x)
        use kind_spec
        real(rprec) :: x
        is_zero = abs(x) < tiny(x)
        end function is_zero

      end module NumParams

      module LoopCtrl
      integer :: iloop
      end module LoopCtrl

      module SvdCtrl
      use kind_spec
      integer :: mstrt, mstep, mkeep, mdspw, noaccuracy
      real(rprec) ::  curwt, trgwt
      end module SvdCtrl

      module OutCtrl
      use kind_spec
      integer :: w_psurf, w_csurf, w_bnuv, w_jsurf, w_xerr, w_svd
      end module OutCtrl

      module Vmeshes
      use kind_spec
      integer :: nu, nv, nuv, nuvh, nu1, nv1, nuv1, nuvh1,
     1    mf, nf, npol, ntor, md, nd, nmax, mnd
      end module Vmeshes

      module Vbiopre
      use Vmeshes
      real(rprec), allocatable, dimension(:) ::
     1   vx, vy, vz, dx, dy, dz
      end module Vbiopre

      module Vbnorm
      use Vmeshes
      real(rprec), dimension(:), allocatable :: bn_ext
      character*200 :: extension
      end module Vbnorm

      module Vcen
      use kind_spec
      real(rprec) ::  xp, yp, zp, bmod, bmod1, bx, by, bz,
     1  dbxx, dbyx, dbzx, dbxy, dbyy, dbzy, dbxz, dbyz, dbzz,
     2  xgb, ygb, zgb, esm, esm1, rl, xmu, xmu2, vits2
      end module Vcen

      module Vculine1
      use Vmeshes
      integer :: nw
      real(rprec), allocatable, dimension(:) :: xw, yw, zw, curre
      end module Vculine1

      module Vdiagno2
      use Vmeshes
      real(rprec), dimension(:,:), allocatable :: pot
      end module Vdiagno2

      module Vdiagno3
      use Vmeshes
      real(rprec), dimension(:), allocatable :: pote
      end module Vdiagno3

      module Vfourier2
      use kind_spec
      integer, dimension(13) :: ifaxv, ifax
      real(rprec), dimension(:), allocatable :: trigsv, trigs
      end module Vfourier2

      module Voptim1
      use Vmeshes
      real(rprec) :: averro, ermax, ermin
      end module Voptim1

      module Vprecal1
      use kind_spec
      real(rprec) :: alp, alu, alv, alvp, fnuv, fnv
      integer :: np
      end module Vprecal1

      module Vprecal10
      use Vmeshes
      real(rprec), dimension(:,:), allocatable :: conv1, sinv1
      end module Vprecal10

      module Vprecal2
      use kind_spec
      real(rprec) :: alu1, alv1, alvp1, fnuv1, fnv1
      end module Vprecal2

      module Vprecal3
      use Vmeshes
      real(rprec), dimension(:), allocatable :: factor
      real(rprec), dimension(:), allocatable :: eps
      end module Vprecal3

      module Vprecal4
      use Vmeshes
      real(rprec), dimension(:), allocatable :: cok, sik
      end module Vprecal4

      module Vprecal5
      use Vmeshes
      real(rprec), dimension(:), allocatable :: coh, sih
      end module Vprecal5

      module Vprecal6
      use Vmeshes
      real(rprec), dimension(:,:), allocatable :: comu, simu
      end module Vprecal6

      module Vprecal7
      use Vmeshes
      real(rprec), dimension(:,:), allocatable :: conv, sinv
      end module Vprecal7

      module Vprecal8
      use Vmeshes
      real(rprec), dimension(:), allocatable :: coh1, sih1
      end module Vprecal8

      module Vprecal9
      use Vmeshes
      real(rprec), dimension(:,:), allocatable :: comu1, simu1
      end module Vprecal9

      module Vsolver1
      use Vmeshes
      real(rprec), dimension(:), allocatable :: a, work,
     1   result, dumv, bpx, bpy, bpz, svdw, daccu, bmag
      real(rprec), dimension(:,:), allocatable :: ab,
     1      bfx, bfy, bfz, bfn, svdv
      end module Vsolver1

      module Vsurfac10
      use Vmeshes
      real(rprec), dimension(:), allocatable :: x1u, y1u, z1u
      end module Vsurfac10

      module Vsurfac11
      use Vmeshes
      real(rprec), dimension(:), allocatable :: x1v, y1v, z1v
      end module Vsurfac11

      module Vsurfac12
      use Vmeshes
      real(rprec), dimension(:), allocatable :: r1uu, z1uu, r1u
      end module Vsurfac12

      module Vsurfac13
      use Vmeshes
      real(rprec), dimension(:), allocatable :: snx1, sny1, snz1,
     1   dsur1
      real(rprec) :: fla
      end module Vsurfac13

      module Vsurfac14
      use Vmeshes
      real(rprec), dimension(:), allocatable :: bex, bey, bez, ben
      end module Vsurfac14

      module Vsurface1
      use Vmeshes
      real(rprec), dimension(:), allocatable :: x, y, z, r
      end module Vsurface1

      module Vsurface2
      use Vmeshes
      real(rprec), dimension(:), allocatable :: xu, yu, xv, ru, rv
      end module Vsurface2

      module Vsurface3
      use Vmeshes
      real(rprec), dimension(:), allocatable :: yv, zu, zv
      end module Vsurface3

      module Vsurface4
      use Vmeshes
      real(rprec), dimension(:), allocatable :: xuu, yuu, zuu, ruu
      end module Vsurface4

      module Vsurface5
      use Vmeshes
      real(rprec), dimension(:), allocatable :: xuv, yuv, zuv, ruv
      end module Vsurface5

      module Vsurface6
      use Vmeshes
      real(rprec), dimension(:), allocatable :: xvv, yvv, zvv, rvv
      end module Vsurface6

      module Vsurface7
      use Vmeshes
      real(rprec), dimension(:), allocatable :: snx, sny, snz
      end module Vsurface7

      module Vsurface8
      use Vmeshes
      real(rprec), dimension(:), allocatable :: xcur, ycur, zcur
      end module Vsurface8

      module Vsurface9
      use Vmeshes
      real(rprec), dimension(:), allocatable :: x1, y1, z1, r1
      end module Vsurface9

      module Vsurfcurdiag
      use Vmeshes
      real(rprec) :: curvra, curvri, cudema, cudemi, avcude
      end module Vsurfcurdiag

      module Vvacuum1
      use Vmeshes
      integer   ms, ns
      real(rprec), dimension(:,:), allocatable :: cr, cz
      end module Vvacuum1

      module Vvacuum2
      use Vmeshes
      integer   ms1, ns1
      real(rprec), dimension(:,:), allocatable :: cr1, cz1, cl1
      real(rprec) :: iota_edge, phip_edge, curpol
      end module Vvacuum2

      module Vvacuum3
      use kind_spec
      integer   ibex
      real(rprec) ::  cup, cut
      end module Vvacuum3

      module Vvacuum4
      use Vmeshes
      integer   ms2, ns2
      real(rprec), dimension(:,:), allocatable :: cr2, cz2
      end module Vvacuum4

      module Vvacuum5
      use Vmeshes
      integer   ms3, ns3
      real(rprec), dimension(:,:), allocatable :: cr3, cz3
      end module Vvacuum5

      module Vvacuum6
      use Vmeshes
      integer   msf, nsf
      real(rprec), dimension(:,:), allocatable :: cf, sf
      end module Vvacuum6

