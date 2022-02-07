      module Vpname1
      use kind_spec
      integer, parameter :: nthpts = 99
      integer, parameter :: nfloops = 40
      integer, parameter :: nbloops = 40
      integer, parameter :: nbsetsp = 2
      integer, parameter :: nbcoilsp = 50
      integer, parameter :: nbrbz_id = 1
      integer, parameter :: nmirnov_id = 2
      integer, parameter :: nthreed2 = 39
      integer, parameter :: max_no_angles_p = 10
      integer :: nrt, max_no_angles
      real(rprec), parameter :: dmu0 = 1.256637e-06_dp
      end module Vpname1

      module Vpname2
      use kind_spec
      use Vpname1
      integer, dimension(:), allocatable :: ixm, ixn, ns_ball
      integer :: nbrbz, nmirnov,  nmirnovset
      real(rprec), dimension(:), allocatable :: dbzcods, ffp, ub, 
     1   iotaf, darea, iotazb, psi, itors, ipols
      real(rprec), dimension(:), allocatable :: brcoil, plbrfld, 
     1   brbc, bzcoil, plbzfld, bzbc
      real(rprec), dimension(:), allocatable :: sqrt_phimod, 
     1   phimod, dbsubudstz, bsubutz, bsubvtz, jboots
      real(rprec), dimension(:,:), allocatable :: balloon_grate
      real(rprec), dimension(:), allocatable :: init_theta_v,init_zeta_v
      real(rprec) :: itor, ipol
      character   limpos_file*120, gmeta_file*40
      character   cdate*30, runlabel*50
      character*100 :: threed2_file, spreadsheet, answers_file, 
     1     bal_grate_file, bal_input_file
      end module Vpname2

      module Vpname3
      use kind_spec
      integer, parameter :: nlimd = 50
      integer, parameter :: nlimset = 2
      integer, parameter :: nltot = nlimd*nlimset
      integer, parameter :: ndim1 = 100
      integer, parameter :: n2max = 40
      integer, parameter :: nsetsmax = 13
      integer :: initnml = 0, ltouch
      integer, dimension(nsetsmax) :: ipset
      real(rprec) :: rx1sv, rx2sv, zy1sv, zy2sv, condifsv
      character :: tokid1*20
      integer :: nrgridsv, nzgridsv
      real(rprec), dimension(ndim1) :: rgrid, zgrid
      real(rprec), dimension(:,:), allocatable :: pfltot
      real(rprec), dimension(:,:,:), allocatable :: pflcoilgr
      real(rprec), dimension(nlimd,nlimset) :: xlim, ylim
      real(rprec), dimension(nltot) :: xlimt, ylimt, seplimt
      integer :: nrgrid_17 = 60, nzgrid_17 = 60, 
     1           nrgrid_18 = 60, nzgrid_18 = 60
      real(rprec) :: rx1_17 = 1.5, rx2_17 = 4.0, zy1_17 = -1.40, 
     1        zy2_17 = 1.40,
     1        rx1_18 = 0.6, rx2_18 = 5.0, zy1_18 = -2.88, zy2_18 = 2.88
      end module Vpname3

      module Vindat2
      use kind_spec
      integer   mn0
      real(rprec) :: hs, ohs, twopi
      end module Vindat2

      module Vmagaxis
      use kind_spec
      real(rprec) :: rmagaxis, zmagaxis
      end module Vmagaxis

      module Vplotdata
      use kind_spec
      integer :: ndata, j1
      integer, parameter :: iselect = 2 
      real, parameter :: lx1 = 0.10
      real, parameter :: lx2 = 0.50
      real, parameter :: lx3 = 0.50
      real, parameter :: lx4 = 0.90
      real, parameter :: ly1 = 0.10
      real, parameter :: ly2 = 0.50
      real, parameter :: ly3 = 0.50
      real, parameter :: ly4 = 0.90
      real, dimension(4,5) :: gwnd
      real, dimension(:), allocatable :: rdata, pldata

      data (gwnd(j1,1),j1=1,4) / lx1, lx2, ly3, ly4 /
      data (gwnd(j1,2),j1=1,4) / lx3, lx4, ly3, ly4 /
      data (gwnd(j1,3),j1=1,4) / lx1, lx2, ly1, ly2 /
      data (gwnd(j1,4),j1=1,4) / lx3, lx4, ly1, ly2 /
      data (gwnd(j1,5),j1=1,4) / lx1, lx4, ly1, ly4 / 

      end module Vplotdata

      module Vpltcn2
      use kind_spec
      real(rprec) :: rmin, rmax, zmin, zmax
      end module Vpltcn2

      module Vpltcn6
      integer   lxc, lyc, noxc, noyc
      logical   newframe
      end module Vpltcn6

      module Vrzarray
      use kind_spec
      real(rprec), dimension(:), allocatable :: r, z
      real(rprec), dimension(:), allocatable :: ru, zu
      end module Vrzarray

      module Vtraneq
      use kind_spec
      integer, parameter :: maxmom = 20
      integer, parameter :: maxpsi = 502
      integer, parameter :: maxtht = 502
      integer, parameter :: maxord = 10
      integer   itype, ieqedg, numtht, intord, ieqax, mombnd, nequil, 
     1   nstep, njav, mom, ncycle, irun, isyms, ipest, npsit, kmax
      real(rprec), dimension(maxmom) :: rmb, ymb
      real(rprec), dimension(maxpsi) :: psi, press, dpdpsi, q, 
     1 dqdpsi, g, ggp, ainvrs
      real(rprec), dimension(maxtht) :: xbound, zbound
      real(rprec) :: zlowp, zlowq, zlowps, zcycp, zcycq, zcycps, 
     1   pfilt1, pfilt2, pfilt3, r0b, shotnm, runum, time, btor, rtor, 
     2   eqcamp, times, xaxes, zmags, gzeros, apls, betas, betaps, 
     3   ali2s, qsaws, psimins, psilims
      character   thdfil*80, flnmi1*80, flnmo1*80
      real(rprec), dimension(maxpsi) :: jdotbc
      real(rprec) :: batot, bapol, bator, baxis
      end module Vtraneq


      module Vcpmpcm1
      integer, parameter :: ntmax = 101
      real, dimension(ntmax) :: comxi, comyj
      end module Vcpmpcm1

      module Vcpmpcm2
      real, dimension(:,:), allocatable :: comxij, comyij
      end module Vcpmpcm2

      module Vcpmpinf
      integer   mxdim, nydim
      end module Vcpmpinf

      module Vthrint
      integer   ithrmj, ithrmn, ithrtx
      end module Vthrint
