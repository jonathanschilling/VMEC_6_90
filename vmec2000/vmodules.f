      module vmec_params
      use kind_spec
      use vparams, only: mpold
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: meven = 0, modd = 1 
      integer, parameter :: ndamp = 10 
      integer, parameter :: ns4 = 25 

      integer, private :: ink
      integer, parameter, dimension(0:mpold) ::
     1  jmin1 = (/ 1,1,(2,ink=2,mpold) /),        !starting js(m) values where R,Z are non-zero
     2  jmin2 = (/ 1,2,(2,ink=2,mpold) /),        !starting js(m) values for which R,Z are evolved
     3  jlam  = (/ 2,2,(2,ink=2,mpold) /)         !starting js(m) values for which Lambda is evolved

      character*(*), parameter :: version_ = '6.90'
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ntmax
      real(rprec), allocatable :: mscale(:), nscale(:)
      real(rprec) :: signgs

      end module vmec_params
      
      module vacmod0 
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: mf, nf, nu, nv, mf1, nf1, mnpd, mnpd2, 
     1           nvp, nuv, nu2, nu3, nuv2
!-----------------------------------------------
      end module vacmod0 


      module vmec_dim
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: mpol1, ntor1, mnmax, ntheta1, ntheta2, ntheta3,
     1           nznt, nrzt, mns, mnsize, ns, ns1
c-----------------------------------------------
      end module vmec_dim
 
 
      module vmec_io
      use kind_spec
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec) :: VolAvgB, IonLarmor, Aminor_p, Rmajor_p,
     1  betatot, betapol, betator, betaxis, b0, volume_plasma,
     2  cross_area
      real(rprec) :: rmax_surf, rmin_surf, zmax_surf
      end module vmec_io


      module vmec_persistent
      use kind_spec
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, dimension(:), allocatable :: ixm, jmin3
      real(rprec), dimension(:,:), allocatable :: cosmu, sinmu,
     1   cosmum, sinmum, cosmui, cosmumi, sinmui, sinmumi,
     2   cosnv, sinnv, cosnvn, sinnvn
      real(rprec), dimension(:), allocatable ::
     1      xm, xn
c-----------------------------------------------
      end module vmec_persistent

      module xstuff
      use kind_spec
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec), dimension(:), allocatable :: 
     1    xc, gc, xcdot, xstore, scalxc
C-----------------------------------------------
      end module xstuff
 
      module vmec_main
      use vmec_dim
      use vmec_input
      use vmec_persistent
      use vmec_params, only: ndamp
      use vparams
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec), dimension(:,:), allocatable ::
     1    ard, arm, brd, brm, azd, azm, bzd, bzm, bmin, bmax
      real(rprec), dimension(:), allocatable ::
     1    crd, iotaf, mass, phi, presf, beta_vol, jcuru, jcurv, jdotb,
     2    buco, bvco, bdotgradv, equif, specw, tcon, fpsi, psi, yellip,
     3    yinden, ytrian, yshift, ygeo, overr, sm, sp, iotas, phips,
     4    pres, vp, jpar2, jperp2, bdotb, blam, clam, dlam
      real(rprec), dimension(:,:,:,:), allocatable :: faclam
      real(rprec), dimension(0:mpol1d,3) :: xmpq
      real(rprec), dimension(0:mpol1d) :: faccon
      real(rprec) :: dcon, currv, aspect, hs, ohs, voli, r01,
     1   signiota, rc0mse, r00, r0scale, z00, dkappa, fsqsum0,
     2   pressum0, fnorm,
     3   fsqr, fsqz, fsql, fnorm1, fsqr1, fsqz1, fsql1, fsq, fedge,
     4   wb, wp
      real(rprec), dimension(nstore_seq) :: fsqt, wdot
      real(rprec) :: ftolv, otav
      real(rprec), dimension(ndamp) :: otau
      real(rprec), dimension(0:10) :: timer
      real(rprec), dimension(:,:,:), allocatable, target ::
     1    rmn_bdy, zmn_bdy
      real(rprec), dimension(:,:), allocatable :: bsqsav
      real(rprec), dimension(:), allocatable :: bsubu0, dbsq, rbsq
      real(rprec) :: rbtor, rbtor0, ctor, delbsq
      real(rprec), dimension(ndatafmax) ::
     1  spfa, spfa2, hp, sifa, sifa2, hi
      logical :: lthreed
      integer, dimension(:), allocatable :: ireflect
      integer :: multi_ns_grid, iequi, irst,
     1    iter1, iter2, ijacob, itfsq, iresidue, neqs, neqs1,
     2    neqs2, irzloff, ivac, ndatap, ndatai
c-----------------------------------------------
      end module vmec_main

      module realspace
      use kind_spec
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec), dimension(:,:), allocatable ::
     1   r1, ru, rv, z1, zu, zv, rcon, zcon
      real(rprec), dimension(:), allocatable :: guu, guv,
     1   gvv, ru0, zu0, rcon0, zcon0, phip, shalf, sqrts, wint
C-----------------------------------------------
      end module realspace
 
      module vforces
      use kind_spec
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec), dimension(:), allocatable, target ::
     1    armn, azmn, brmn, bzmn, crmn, czmn, blmn, clmn
      real(rprec), pointer, dimension(:) ::
     1    armn_e, armn_o, azmn_e, azmn_o,
     2    brmn_e, brmn_o, bzmn_e, bzmn_o,
     3    crmn_e, crmn_o, czmn_e, czmn_o, blmn_e,
     4    blmn_o, clmn_e, clmn_o
c-----------------------------------------------
      end module vforces

      module vac_persistent
      use kind_spec
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, dimension(:), allocatable :: imirr
      real(rprec), dimension(:), allocatable :: sinper, cosper,
     1   sinuv, cosuv, sin2v, tanu, tanv, xmpot, xnpot, csign
      real(rprec), dimension(:,:), allocatable :: sinu, cosu,
     1 sinv, cosv, sinui, cosui, sinu1, cosu1, sinv1, cosv1
      real(rprec), dimension(:,:,:), allocatable :: cmns
      end module vac_persistent


      module vacwires
      use kind_spec
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec), dimension(:), allocatable ::
     1   xw, yw, zw, vx, vy, vz, dx, dy, dz
      end module vacwires

      module vacmod
      use vacmod0
      use vac_persistent
      use vmec_input, only: lasym
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nr0b, np0b, nfper0, nz0b, nfper
      real(rprec), dimension(:), allocatable, target :: potvac
      real(rprec), dimension(:), allocatable :: bvecsav, amatsav,
     1   bexni, brv, bphiv, bzv, bpolvac, bsqvac, r1b, rub, rvb, z1b,
     2   zub, zvb, bexu, bexv, bexn, auu, auv, avv, snr, snv, snz, drv,
     3   guu_b, guv_b, gvv_b, rb2, rcosuv, rsinuv,
     4   bredge, bpedge, bzedge
      real(rprec) :: bsubvvac, rminb, zminb, rmaxb, zmaxb,
     1   delrb, delzb, pi2,
     2   pi3, pi4, alp, alu, alv, alvp, onp, onp2
      end module vacmod


      module mgrid_mod
      use kind_spec
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nbvac
      real(rprec), dimension(:,:,:), allocatable :: btemp
      real(rprec), dimension(:,:), allocatable :: bvac
      character*256 :: mgrid_path
      end module mgrid_mod

      module vmercier
      use vparams, only: nsd, rprec, dp
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec), dimension(nsd) ::
     1   Dshear, Dwell, Dcurr, Dmerc, Dgeod
      end module vmercier

      module vsvd
      use vmec_input
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, dimension(:), allocatable ::
     1   indexr, imid, needflx, nbcoils
      integer, dimension(:,:), allocatable :: iconnect, needbfld
      integer, dimension(:), allocatable :: nk_ia, nk_pa, nk_ib,
     1   nk_pb, indexs2, indexu2, indexs1, indexu1, isortr, isorts
      integer, dimension(nigroup) :: nsetsn
      integer :: nmeasurements, nbcoil_max, nsets, nlim_max, imovephi, 
     1   imse2, icurrout, islope, itse2, iphidiam, ipresin, ipresout, 
     2   nchistp
      integer, dimension(mstp) :: nchi2
      integer :: nchisaddle, nchitot, nrgrid, nzgrid, nchidia, 
     1   nchipres, nchimse, nlim
      integer, dimension(nlimset) :: limitr
      integer :: nobd, nobser, nextcur, nbfldn, nbsets, nbcoilsn
      real(rprec), dimension(:), allocatable, target ::
     1      datamse, qmid, shear, presmid, alfa, curmid,
     2      curint, psimid, ageo, volpsi, phimid
      real(rprec), dimension(:), allocatable ::
     1   current, rm2, vrm2, ovrm2, ochip, presph, 
     2   presint, w_ia, w1_ia, u_ia, u1_ia, w_pa, w1_pa, u_pa, u1_pa, 
     3   w_ib, w1_ib, u_ib, u1_ib, w_pb, w1_pb, u_pb, u1_pb, rmid, 
     4   isplinef, isplineh, psplinef, psplineh, sthom, delse2, delso2, 
     5   pcalc, delse1, delso1, starkcal, qmeas, qcalc, fpsical, 
     6   stark_weight, rsort, rsort0, xobser, xobsqr, zobser, dsiext, 
     7   psiext, plflux, b_chi
      real(rprec), dimension(:,:), allocatable :: 
     1   pm, im, unpsiext, plbfld, rbcoil, zbcoil, 
     2   abcoil, bcoil, rbcoilsqr, rlim, zlim, reslim, seplim
      real(rprec), dimension(:,:,:,:), allocatable :: pmb, imb
      real(rprec), dimension(:,:,:), allocatable :: 
     1   dbcoil, pfcspec
      real(rprec), dimension(jngrn) :: 
     1   yf, dyf, qsq, yek, yeq, dyek, dyeq
      real(rprec) :: odqk2, rstarkmin, rstarkmax, torflux, ppeak, 
     1   pfac, phifac, phifsave, rthompeak, pthommax, rthommax, 
     2   rthommin, delphid, dlim_min, rlim_min, zlim_min, router, 
     3   rinner, apres, aminor, grmse, gphifac, rstepx0,
     4   rsfac, raxmse, rwidth, errsvd
      real(rprec), dimension(jchix1) :: chisqerr
      real(rprec), dimension(jchix1,mstp) :: chi2
      real(rprec) :: total_chi_square_n, total_chisq_n0, 
     1   total_chi_square, total_saddle_chi, total_b_chi, 
     2   total_pres_chi, total_mse_chi, total_chi_cur, total_chi_dia,
     3   scstark, scthom, rx1, rx2, zy1, 
     4   zy2, condif, flmwgt, bcwgt, tswgt, msewgt
      logical :: lpprof
      character*(20) :: tokid
      character*8, dimension(:), allocatable :: dsilabel, bloopnames
      character*30, dimension(:), allocatable ::  curlabel
C-----------------------------------------------
      end module vsvd

      module vspline
      use vmec_input
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, dimension(1) :: jspmin
      integer :: iknots
      real(rprec), dimension(:), allocatable :: hthom, ythom,
     1   y2thom, pknots, hstark, y2stark, ystark, sknots
C-----------------------------------------------
      end module vspline

      module csplinx
      use kind_spec
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nptsx
      real(rprec), dimension(:), allocatable :: rmidx,hmidx,
     1   wmidx,qmidx,tenmidx,ymidx,y2midx
c-----------------------------------------------
      end module csplinx
