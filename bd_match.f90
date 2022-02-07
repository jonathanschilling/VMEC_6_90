#!/bin/sh
#---------------------------------------------------------------------
cat > temp.c << "EOC"
cat > vmodules.f << "EOF"
      module boundaries
      use kind_spec
      integer ::  mnmax_vv, dimnow, mnmax_bd, mnmax_pl, mpol_pl, ntor_pl
      integer ::  nrad
      real(rprec), dimension(:),  allocatable  :: 
     1  rmnc_vv, zmns_vv, xm_vv, xn_vv,
     1  rmnc_bd, zmns_bd, xm_bd, xn_bd, 
     1  rmns_vv, zmnc_vv, rmns_bd, zmnc_bd,
     1  rmnc_pl, zmns_pl, xm_pl, xn_pl, rmns_pl, zmnc_pl
      logical vv , bd
      end module boundaries
      module Vpname1
      use kind_spec
      integer, parameter :: nthpts = 99
      integer, parameter :: nfloops = 40
      integer, parameter :: nbloops = 40
      integer, parameter :: nbsetsp = 2
      integer, parameter :: nbcoilsp = 50
      integer, parameter :: nbrbz_id = 1
      integer, parameter :: nmirnov_id = 2
      integer :: nrt
      real(kind=rprec), parameter :: dmu0 = 1.256637d-06
      end module Vpname1

      module Vpname2
      use kind_spec
      use Vpname1
      integer, dimension(:), allocatable :: ixm, ixn
      integer :: nbrbz, nmirnov,  nmirnovset
      real(kind=rprec), dimension(:), allocatable :: dbzcods, ffp, ub, 
     1   iotaf, darea, iotazb, psi, itors, ipols
      real(kind=rprec), dimension(:), allocatable :: brcoil, plbrfld, 
     1   brbc, bzcoil, plbzfld, bzbc
      real(kind=rprec), dimension(:), allocatable :: sqrt_phimod, 
     1   phimod, dbsubudstz, bsubutz, bsubvtz
      real(kind=rprec) :: itor, ipol
      character   limpos_file*40, gmeta_file*40
      character   cdate*30, runlabel*60
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
      real(kind=rprec) :: rx1sv, rx2sv, zy1sv, zy2sv, condifsv
      character :: tokid1*20
      integer :: nrgridsv, nzgridsv
      real(kind=rprec), dimension(ndim1) :: rgrid, zgrid
      real(kind=rprec), dimension(:,:), allocatable :: pfltot
      real(kind=rprec), dimension(:,:,:), allocatable :: pflcoilgr
      real(kind=rprec), dimension(nlimd,nlimset) :: xlim, ylim
      real(kind=rprec), dimension(nltot) :: xlimt, ylimt, seplimt
      integer :: nrgrid_17 = 60, nzgrid_17 = 60, 
     1           nrgrid_18 = 60, nzgrid_18 = 60
      real(kind=rprec) :: rx1_17 = 1.5, rx2_17 = 4.0, zy1_17 = -1.40, 
     1        zy2_17 = 1.40,
     1        rx1_18 = 0.6, rx2_18 = 5.0, zy1_18 = -2.88, zy2_18 = 2.88
      end module Vpname3

      module Vindat2
      use kind_spec
      integer   mn0
      real(kind=rprec) :: hs, ohs, twopi
      end module Vindat2

      module Vmagaxis
      use kind_spec
      real(kind=rprec) :: rmagaxis, zmagaxis
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
      real, dimension(100) :: rdata, pldata

      data (gwnd(j1,1),j1=1,4) / lx1, lx2, ly3, ly4 /
      data (gwnd(j1,2),j1=1,4) / lx3, lx4, ly3, ly4 /
      data (gwnd(j1,3),j1=1,4) / lx1, lx2, ly1, ly2 /
      data (gwnd(j1,4),j1=1,4) / lx3, lx4, ly1, ly2 /
      data (gwnd(j1,5),j1=1,4) / lx1, lx4, ly1, ly4 / 

      end module Vplotdata

      module Vpltcn2
      use kind_spec
      real(kind=rprec) :: rmin, rmax, zmin, zmax
      end module Vpltcn2
      module pgplotinc
C-----------------------------------------------------------------------
C This module is simply a copy of pgplot.inc from version 5.2
C This module may need updating with future releases of PGPLOT
C Contents used for basic coordinate conversions, e.g., xyitoxys
C Ed Lazarus
C-----------------------------------------------------------------------
C PGPLOT: common block definition.
C-----------------------------------------------------------------------
C Maximum number of concurrent devices (should match GRIMAX).
C-----------------------------------------------------------------------
      INTEGER PGMAXD
      PARAMETER (PGMAXD=8)
C-----------------------------------------------------------------------
C Indentifier of currently selected device.
C-----------------------------------------------------------------------
      INTEGER PGID
C-----------------------------------------------------------------------
C Device status (indexed by device identifier).
C-----------------------------------------------------------------------
C PGDEVS  =0 if device is not open; 1 if device is open.
C PGADVS  Set to 0 by PGBEGIN, set to 1 by PGPAGE; used to suppress
C         the prompt for the first page.
C PROMPT  If .TRUE., ask user before clearing page; set by PGASK
C         and (indirectly) by PGBEGIN, used in PGENV.
C PGBLEV  Buffering level: incremented by PGBBUF, decremented by
C         PGEBUF.
C PGPFIX  TRUE if PGPAP has been called, FALSE otherwise.
C
      INTEGER PGDEVS(PGMAXD), PGADVS(PGMAXD), PGBLEV(PGMAXD)
      LOGICAL PGPRMP(PGMAXD), PGPFIX(PGMAXD)
C-----------------------------------------------------------------------
C Panel parameters (indexed by device identification).
C-----------------------------------------------------------------------
C NX      Number of panels in x direction
C NY      Number of panels in y direction
C NXC     Ordinal number of current X panel
C NYC     Ordinal number of current Y panel
C XSZ     X dimension of panel (device units)
C YSZ     Y dimension of panel (device units)
C PGROWS  TRUE if panels are used in row order, FALSE for column
C         order.
C
      INTEGER PGNX  (PGMAXD), PGNY  (PGMAXD)
      INTEGER PGNXC (PGMAXD), PGNYC (PGMAXD)
      REAL    PGXSZ (PGMAXD), PGYSZ (PGMAXD)
      LOGICAL PGROWS(PGMAXD)
C-----------------------------------------------------------------------
C Attributes (indexed by device identification).
C-----------------------------------------------------------------------
C PGCLP   clipping enabled/disabed
C PGFAS   fill-area style
C PGCHSZ  character height
C PGAHS   arrow-head fill style
C PGAHA   arrow-head angle
C PGAHV   arrow-head vent
C PGTBCI  text background color index
C PGMNCI  lower range of color indices available to PGGRAY/PGIMAG
C PGMXCI  upper range of color indices available to PGGRAY/PGIMAG
C PGITF   type of transfer function used by PGGRAY/PGIMAG
C PGHSA   hatching line angle
C PGHSS   hatching line separation
C PGHSP   hatching line phase
C
      INTEGER PGCLP (PGMAXD)
      INTEGER PGFAS (PGMAXD)
      REAL    PGCHSZ(PGMAXD)
      INTEGER PGAHS (PGMAXD)
      REAL    PGAHA (PGMAXD)
      REAL    PGAHV (PGMAXD)
      INTEGER PGTBCI(PGMAXD)
      INTEGER PGMNCI(PGMAXD)
      INTEGER PGMXCI(PGMAXD)
      INTEGER PGITF (PGMAXD)
      REAL    PGHSA (PGMAXD)
      REAL    PGHSS (PGMAXD)
      REAL    PGHSP (PGMAXD)
C-----------------------------------------------------------------------
C Viewport parameters (indexed by device identification); all are device
C coordinates:
C-----------------------------------------------------------------------
C PGXOFF  X coordinate of blc of viewport.
C PGYOFF  Y coordinate of blc of viewport.
C PGXVP   X coordinate of blc of viewport, relative to blc of subpage.
C PGYVP   Y coordinate of blc of viewport, relative to blc of subpage.
C PGXLEN  Width of viewport. 
C PGYLEN  Height of viewport.
C
      REAL   PGXOFF(PGMAXD), PGYOFF(PGMAXD)
      REAL   PGXVP (PGMAXD), PGYVP (PGMAXD)
      REAL   PGXLEN(PGMAXD), PGYLEN(PGMAXD)
C-----------------------------------------------------------------------
C Scaling parameters (indexed by device identification):
C-----------------------------------------------------------------------
C PGXORG  device coordinate value corresponding to world X=0
C PGYORG  device coordinate value corresponding to world Y=0
C PGXSCL  scale in x (device units per world coordinate unit)
C PGYSCL  scale in y (device units per world coordinate unit)
C PGXPIN  device x scale in device units/inch
C PGYPIN  device y scale in device units/inch
C PGXSP   Character X spacing (device units)
C PGYSP   Character Y spacing (device units)
C
      REAL   PGXORG(PGMAXD), PGYORG(PGMAXD)
      REAL   PGXSCL(PGMAXD), PGYSCL(PGMAXD)
      REAL   PGXPIN(PGMAXD), PGYPIN(PGMAXD)
      REAL   PGXSP (PGMAXD), PGYSP (PGMAXD)
C-----------------------------------------------------------------------
C Window parameters (indexed by device identification); all are world
C coordinate values:
C-----------------------------------------------------------------------
C PGXBLC  world X at bottom left corner of window
C PGXTRC  world X at top right corner of window
C PGYBLC  world Y at bottom left corner of window
C PGYTRC  world Y at top right corner of window
C
      REAL    PGXBLC(PGMAXD), PGXTRC(PGMAXD)
      REAL    PGYBLC(PGMAXD), PGYTRC(PGMAXD)
C-----------------------------------------------------------------------
C The following parameters are used in the contouring routines to pass
C information to the action routine. They do not need to be indexed.
C-----------------------------------------------------------------------
C TRANS   Transformation matrix for contour plots; copied
C         from argument list by PGCONT and used by PGCP.
C
      INTEGER PGCINT, PGCMIN
      REAL    TRANS(6)
      CHARACTER*32 PGCLAB
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      COMMON/PGPLT1/ PGID,PGDEVS,PGADVS,PGNX,  PGNY,  PGNXC, PGNYC ,
     1        PGXPIN,PGYPIN,PGXSP, PGYSP, PGXSZ, PGYSZ,
     2        PGXOFF,PGYOFF,PGXVP, PGYVP, PGXLEN,PGYLEN,PGXORG,PGYORG,
     3        PGXSCL,PGYSCL,PGXBLC,PGXTRC,PGYBLC,PGYTRC,TRANS,
     4        PGPRMP,PGCLP, PGFAS, PGCHSZ,PGBLEV,PGROWS,
     5        PGAHS, PGAHA, PGAHV, PGTBCI,PGMNCI,PGMXCI,PGCINT,PGCMIN,
     6        PGPFIX,PGITF, PGHSA, PGHSS, PGHSP
      COMMON/PGPLT2/ PGCLAB
      SAVE    /PGPLT1/
      SAVE    /PGPLT2/
C-----------------------------------------------------------------------
      end module pgplotinc
      module grpckg1inc

C-----------------------------------------------------------------------
C      From PGPLOT V5.2 -- Be sure this is up-to-date
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C            Include file for GRPCKG
C Modifications:
C   29-Jan-1985 - add HP2648 (KS/TJP).
C   16-Sep-1985 - remove tabs (TJP).
C   30-Dec-1985 - add PS, VPS (TJP).
C   27-May-1987 - remove ARGS, NULL, PS, VPS, QMS, VQMS, HIDMP,
C                 HP7221, GRINL (TJP).
C    6-Jun-1987 - remove PRTX, TRILOG, VERS, VV (TJP).
C   11-Jun-1987 - remove remaining built-in devices (TJP).
C    5-Jul-1987 - replace GRINIT, GRPLTD by GRSTAT.
C   16-Aug-1987 - remove obsolete variables.
C    9-Sep-1989 - add SAVE statement.
C   26-Nov-1990 - remove GRCTYP.
C    5-Jan-1993 - add GRADJU.
C    1-Sep-1994 - add GRGCAP.
C   21-Dec-1995 - increase GRIMAX to 8.
C   30-Apr-1997 - remove GRC{XY}SP
C-----------------------------------------------------------------------
C
C Parameters:
C   GRIMAX : maximum number of concurrent devices
C   GRFNMX : maximum length of file names
C   GRCXSZ : default width of chars (pixels)
C   GRCYSZ : default height of chars (pixels)
C
      INTEGER   GRIMAX, GRFNMX
      REAL      GRCXSZ, GRCYSZ
      PARAMETER (GRIMAX = 8)
      PARAMETER (GRFNMX = 90)
      PARAMETER (GRCXSZ =  7.0, GRCYSZ =  9.0)
C
C Common blocks:
C   GRCIDE : identifier of current plot
C   GRGTYP : device type of current plot
C The following are qualified by a plot id:
C   GRSTAT : 0 => workstation closed
C            1 => workstation open
C            2 => picture open
C   GRPLTD :
C   GRDASH : software dashing in effect?
C   GRUNIT : unit associated with id
C   GRFNLN : length of filename
C   GRTYPE : device type
C   GRXMXA : x size of plotting surface
C   GRYMXA : y size of plotting surface
C   GRXMIN : blc of plotting window
C   GRYMIN : ditto
C   GRXMAX : trc of plotting window
C   GRYMAX : ditto
C   GRSTYL : line style (integer code)
C   GRWIDT : line width (integer code)
C   GRCCOL : current color index (integer code)
C   GRMNCI : minimum color index on this device
C   GRMXCI : maximum color index on this device
C   GRCMRK : marker number
C   GRXPRE : previous (current) pen position (x)
C   GRYPRE : ditto (y)
C   GRXORG : transformation variables (GRTRAN)
C   GRYORG : ditto
C   GRXSCL : ditto
C   GRYSCL : ditto
C   GRCSCL : character scaling factor
C   GRCFAC :
C   GRCFNT : character font
C   GRFILE : file name (character)
C   GRGCAP : device capabilities (character)
C   GRPXPI : pixels per inch in x
C   GRPYPI : pixels per inch in y
C   GRADJU : TRUE if GRSETS (PGPAP) has been called
C
      INTEGER   GRCIDE, GRGTYP
      LOGICAL   GRPLTD(GRIMAX), GRDASH(GRIMAX), GRADJU(GRIMAX)
      INTEGER   GRSTAT(GRIMAX)
      INTEGER   GRUNIT(GRIMAX), GRFNLN(GRIMAX), GRTYPE(GRIMAX),
     1          GRXMXA(GRIMAX), GRYMXA(GRIMAX), 
     2          GRSTYL(GRIMAX), GRWIDT(GRIMAX), GRCCOL(GRIMAX),
     3          GRCMRK(GRIMAX), GRIPAT(GRIMAX), GRCFNT(GRIMAX),
     4          GRMNCI(GRIMAX), GRMXCI(GRIMAX)
      REAL      GRXMIN(GRIMAX), GRYMIN(GRIMAX),
     1          GRXMAX(GRIMAX), GRYMAX(GRIMAX)
      REAL      GRXPRE(GRIMAX), GRYPRE(GRIMAX), GRXORG(GRIMAX),
     1          GRYORG(GRIMAX), GRXSCL(GRIMAX), GRYSCL(GRIMAX),
     2          GRCSCL(GRIMAX), GRCFAC(GRIMAX), GRPOFF(GRIMAX),
     3          GRPATN(GRIMAX,8),GRPXPI(GRIMAX),GRPYPI(GRIMAX)
      COMMON /GRCM00/ GRCIDE, GRGTYP, GRSTAT, GRPLTD, GRUNIT,
     1                GRFNLN, GRTYPE, GRXMXA, GRYMXA, GRXMIN, GRYMIN,
     2                GRXMAX, GRYMAX, GRWIDT, GRCCOL, GRSTYL,
     3                GRXPRE, GRYPRE, GRXORG, GRYORG, GRXSCL, GRYSCL,
     4                GRCSCL, GRCFAC, GRDASH, GRPATN, GRPOFF,
     5                GRIPAT, GRCFNT, GRCMRK, GRPXPI, GRPYPI, GRADJU,
     6                GRMNCI, GRMXCI
C
      CHARACTER*(GRFNMX) GRFILE(GRIMAX)
      CHARACTER*11       GRGCAP(GRIMAX)
      COMMON /GRCM01/ GRFILE, GRGCAP
      SAVE /GRCM00/, /GRCM01/
C-----------------------------------------------------------------------
      end module grpckg1inc


      module Vpltcn6
      integer   lxc, lyc, noxc, noyc
      contains
       subroutine newframe
       use pgplotinc, only: pgnxc,pgnyc,pgnx,pgny
       pgnxc=pgnx
       pgnyc=pgny
       return
       end subroutine newframe
      end module Vpltcn6

      module Vrzarray
      use kind_spec
      real(kind=rprec), dimension(:), allocatable :: r, z
      real(kind=rprec), dimension(:), allocatable :: ru, zu
      end module Vrzarray

      module Vtraneq
      use kind_spec
      integer, parameter :: maxmom = 20
      integer, parameter :: maxpsi = 502
      integer, parameter :: maxtht = 502
      integer, parameter :: maxord = 10
      integer   itype, ieqedg, numtht, intord, ieqax, mombnd, nequil, 
     1   nstep, njav, mom, ncycle, irun, isyms, ipest, npsit, kmax
      real(kind=rprec), dimension(maxmom) :: rmb, ymb
      real(kind=rprec), dimension(maxpsi) :: psi, press, dpdpsi, q, 
     1 dqdpsi, g, ggp, ainvrs
      real(kind=rprec), dimension(maxtht) :: xbound, zbound
      real(kind=rprec) :: zlowp, zlowq, zlowps, zcycp, zcycq, zcycps, 
     1   pfilt1, pfilt2, pfilt3, r0b, shotnm, runum, time, btor, rtor, 
     2   eqcamp, times, xaxes, zmags, gzeros, apls, betas, betaps, 
     3   ali2s, qsaws, psimins, psilims
      character   thdfil*80, flnmi1*80, flnmo1*80
      real(kind=rprec), dimension(maxpsi) :: jdotbc
      real(kind=rprec) :: batot, bapol, bator, baxis
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


      module optim
      use optim_params
      implicit none
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: nxc = 300
      integer, parameter :: iout = 54
      real(kind=rprec) :: bigno = 1.e10_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ierr_vmec, ns_booz_max, ns_surf_max
      integer, dimension(nxc) :: nbrho_opt, mbrho_opt
      integer, dimension(2*ntor1d) :: nrz0_opt
      integer :: irho_bdy, irm0_bdy, izm0_bdy,
     1   nrad, num_ai
      integer :: nfp_opt, mpol_opt, mpol1_opt, ntor_opt, 
     1      ntor1_opt, mnmax_opt, iunit_opt, iunit_opt_local,
     2      nextcur_opt, nextcur_vmec
      real(kind=rprec), dimension(-ntord:ntord,0:mpol1d) :: rhobc
      real(kind=rprec) :: wp_opt, wb_opt, rmax_opt, rmin_opt, zmax_opt
      real(kind=rprec) :: aspect_opt, coil_complex_opt
      real(kind=rprec) :: chisq_min, rbtor_opt
      real(kind=rprec), dimension(nsd) :: 
     1  vp_opt, iota_opt, phip_opt, buco_opt, Dmerc_opt, jdotb_opt
      real(kind=rprec), dimension(nsd) ::  pres_opt                            !!COBRA
      real(kind=rprec) :: version_opt, am0_9, am10                             !!COBRA
      real(kind=rprec) :: raxis_old(0:ntord), zaxis_old(0:ntord)
 
      integer, allocatable :: ns_booz(:), ns_surf(:)
      character*120 :: input_file, home_dir, min_ext
      character*130 :: output_file, min_input_file, min_wout_file
c-----------------------------------------------
      end module optim

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

EOF
cat > bd_match.f << "EOF"

      program bd_match

C-----------------------------------------------
c  Ed Lazarus Jan. 2000

c  A version of prout which uses PGPLOT
c  http://astro.caltech.edu/~tjp/pgplot/

c  Aside from a new module in avmodules.f, subroutine and file
c  names are the same as prout with "pgp" prepended.

c  The plots having to do with reconstruction of an experiment are
c  likely to fail at this time (1/26/00)
C-----------------------------------------------

C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      use read_wout_mod, itor_w=>itor
      use boundaries
      use Vpname1
      use Vpname2
      use Vpname3
      use system_mod, only: system
      use Vindat2
      use optim, only: rbc_vv, zbs_vv, mpol_vv, ntor_vv,
     &    rbc_bd, zbs_bd, mpol_bd, ntor_bd,
     1    vv_dist,  vv_dist_rms,  vv_dist_max,
     2    nu_vv, nv_vv,target_rbtor,target_vv,
     1    bd_dist, bd_dist_rms, bd_dist_max,
     2    nu_bd, nv_bd, target_beta, target_cur

      use vmec_input, only: rbc,rbs,zbc,zbs, lfreeb
      use read_namelist_mod
      use safe_open_mod
      use  boundaries
      implicit none
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: lastplot = 18
      integer, parameter :: maxlist = 100
      logical, parameter :: lwstbf = .TRUE.
      character*(*), parameter :: version =
     1   ' VMEC Plotter Version 6.10pg  y2k '

      character*(*), parameter :: tempfile = 'QqZzXLftemp'

C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: id1, id2, pgopen
      integer :: numargs, numchars, lenlist, i,
     1   js, lj, mn, l, j, n, ii, k, m, iunit,
     2   ntor0, lp, ierr
      real(kind=rprec), dimension(:), allocatable :: szc, szb, dummy
      real(kind=rprec) :: r0,tvvon,tvvoff,tbdon,tbdoff
      character(len=100) :: input_id, explist, device
      character(len=11) ::device_type
      character(len=50), dimension(0:lastplot) :: pagedir
      character*120 :: stripped_input_file, in_file

      character*30 timeloc

C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      integer , EXTERNAL :: iargc
C-----------------------------------------------
*
*                 THIS PROGRAM - PROUT - ACCEPTS THE OUTPUT
*                 FROM THE EQUILIBRIUM PROGRAM VMEC (WOUT FILE)
*                 AND PRINTS/PLOTS THE APPROPRIATE DATA
*
*                GET FILE ID EXTENSION FROM COMMAND LINE
*
*                 Modified (Aug, 1993) ... by M. Vachharajan, R Wieland
*                 Removed DISSPLA Graphics, replaced by NCAR Graphics (V3.0),
*                 using a modified version (ncarplts1.)of the PPPL "ncarplts"
*                 graphics package by Derek Fox.
*                 Graphics calls in this program are a mix of AUTOGRAPH
*                 calls (for 2d plots) and NCARPLTS calls for contours
*
*                To generate all plots %xprout fileid
*                To generate selected plots %xprout fileid "1,5 7 12 14,$"
*                To generate no plots %xprout fileid 0
*
*                 Modified Jan 1994 by R. Wieland to generate fileid.EQI
*                 file for stability code analysis. (Cf lwstbf vbl in WRFCN)
*************
*               Bsubs,Bsubu,Bsubv,Gmn  come in on the half-mesh
*               Vp,Mass,Pres(*),phip   come in on the half-mesh
*               Iotas(*), Lmns         come in on the half-mesh
*               Rmnc,Zmns              come in on the full-mesh
*               Jcuru,Jcurv            come in on the full-mesh
*               Gmn ==> Gsqrt                  on the half mesh
*               sknots,ystark,y2stark  spline knots come in on the sqrt(s) grid
*               pknots,ythom,y2thom    spline knots come in on the sqrt(s) grid
*************
*               (*)  Use Iota Spline Knots instead
*               (*)  Use Pressure Spline Knots instead

      timeloc = 'date >> '//tempfile
      call system (timeloc)
      open(unit=99, file=tempfile, status='old', err=100)
      read (99, 5) cdate
    5 format(a)
  100 continue
      close(unit=99, status='delete')


      numchars = 0

*       Get the Command Line Arguments (Case & Plot Selection)

      numargs = iargc()
      if (input_id .eq. '-h' ) then
         write (*, 104)
  104    format(/' ** To generate all plots %xprout fileid'/
     2' ** choose plotting device %xprout fileid  /device'/
     3' ** typically choose /xw (default) or /cps for/device')
         stop
      endif
      if (numargs >= 1) call getarg (1, input_id)
      if (numargs .eq. 2) call getarg (2, device)
      if (numargs .lt. 2) device='/xw'
      if(numargs .eq. 0) then
         write (6, 104)
         write(6,*)
     .'enter fileid <cr> , device <cr>'
         read(5,fmt='(a)')input_id
         read(5,fmt='(a)')device
         if(device .eq. ' ') device='/xw'
         write(6,*)input_id
         write(6,*)device
       endif


*
*       Does the Input File Exist ?
*
      numchars = len_trim(input_id)
      if (numchars .eq. 0) then
         write(6,*) ' MUST ENTER FILE SUFFIX ON COMMAND LINE'
         stop
      endif
      if (index(input_id,'wout.') .eq. 1) then
         input_id = input_id(6:)
         numchars = numchars - 6
      endif
      runlabel = input_id(1:numchars+1)//' '//cdate(1:len_trim(cdate)+1)
      gmeta_file = 'gmeta.'//input_id(1:numchars)

      call read_wout_file('wout.'//trim(input_id),ierr)
      if (ierr .ne. 0 ) stop "ERROR: call read_wout_file"

      if (ierr .eq. 1) stop 'could not read wout file: check extension'
      if (ierr.ne.0 .and. ierr.le.10) stop 'error in plotter read_wout'
      if (niter .le. 0) stop 'VMEC CODE DID NOT RUN PROPERLY!'


      if (imse .ne. (-1)) then
         i = index(mgrid_file,'.')
         if (i .eq. 0) then
            write (*, *)
     1      'MGRID_FILE (in WOUT & INPUT) has incorrect format!'
             stop
         endif
         tokid = mgrid_file(i+1:)
      else
         tokid = ' '
         mgrid_file = ' '
      endif

*****************************
* Read input file for RBC, etc.
* Needed for LFREEB=T
*****************************
C  need to strip comments from input file producved by stellopt
*****************************
!------------------------------------------------------------------------------

!
!     read input (optimum) namelist
!     first, strip any comments from input file (F95-not necessary)
!
      in_file='input.'//trim(input_id)
      call strip_comments(trim(in_file))
!!Produces clean file input_file//'.stripped'
      stripped_input_file = trim(in_file) // '.stripped'
      id1 = 99
      call safe_open (id1, id2, trim(stripped_input_file)  ,
     .  'old', 'formatted')
      if (id2 .ne. 0) then
         print *,  trim(stripped_input_file)  ,
     .   ': input file open error: iostat = ',
     .      id2
      else
         call read_namelist (id1, id2, 'indata')
      endif
      close(id1) !  close and reopen to avoid error on hecate
      if (id2 .ne. 0) then
         write (*, *) ' indata namelist READ error: iostat = ', id2
         ierr = 3
      endif
      ierr=0; id1 = 99; id2=0
      call safe_open (id1, id2, trim(stripped_input_file)  ,
     .  'old', 'formatted')
      if (id2 .ne. 0) then
         print *,  trim(stripped_input_file)  ,
     .   ': input file open error: iostat = ',
     .      id2
      else
         call read_namelist (id1, id2, 'optimum')
      endif
      if (id2 .ne. 0) then
         write (*, *) ' optimum namelist READ error: iostat = ', id2
         ierr = ierr+6
      endif
*****************************
      mpol_pl=mpol-1; ntor_pl=ntor
      n=size(rmnc(1,:))
      ntor_pl=nint(maxval(abs(xn)))/nfp
      mnmax_pl=size(rmnc(:,1))
      mpol_pl=nint(maxval(xm))
      dimnow=(mpol_pl+1)*(2*ntor_pl+1)
      if(.not.allocated(rmnc_pl))allocate(xm_pl(dimnow),
     &xn_pl(dimnow),rmnc_pl(dimnow),zmns_pl(dimnow),rmns_pl(dimnow),
     &zmnc_pl(dimnow))
      rmns_pl=0;zmns_pl=0
      do mn = 1, mnmax_pl
         xm_pl(mn)   = xm(mn)
         xn_pl(mn)   = xn(mn)/nfp
         rmnc_pl(mn) = rmnc(mn,n)
         zmns_pl(mn) = zmns(mn,n)
      end do
      vv_dist = 0; vv_dist_rms = 0; vv_dist_max = 0
      bd_dist = 0; bd_dist_rms = 0; bd_dist_max = 0
      mnmax_vv=0
      tvvon=0;tvvoff=0;tbdon=0;tbdoff=0
      vv=lfreeb .and. rbc_vv(0,0).gt.0
       print *,'vv=',vv,' R_00 =',real(rbc_vv(0,0)),' wait for surfsep'
      if(id2.eq.0) close(id1,status='delete')
      if(id2.ne.0) close(id1)
      if(vv) then
        mnmax_vv = mpol_vv*(2*ntor_vv+1)+ntor_vv+1
        mnmax_vv = (mpol_vv+1)*(2*ntor_vv+1)-ntor_vv
        dimnow=(mpol_vv+1)*(2*ntor_vv+1)-ntor_vv
        if(.not.allocated(xm_vv)) allocate(rmnc_vv(dimnow),
     &   zmns_vv(dimnow), xm_vv(dimnow), xn_vv(dimnow),
     &   rmns_vv(dimnow), zmnc_vv(dimnow))
            rmns_vv=0 ; zmnc_vv=0; rmnc_vv=0;  zmns_vv=0
           i = 0
           do m=0, mpol_vv
             do n=-ntor_vv, ntor_vv
               if (m.eq.0 .and. n.lt.0) cycle
               i = i+1
               xm_vv(i) = m
               xn_vv(i) = n  !  NOT *nfp
               rmnc_vv(i) = rbc_vv(n,m)
               zmns_vv(i) = zbs_vv(n,m)
             end do
           end do
           call second0(tvvoff)
             call surfsep(mnmax_pl, rmnc_pl, zmns_pl, xm_pl,
     1               xn_pl, mnmax_vv, nu_vv, nv_vv, rmnc_vv, zmns_vv,
     2               xm_vv, xn_vv, nfp, 64, vv_dist, vv_dist_rms,
     2               vv_dist_max, ierr)
           call second0(tvvon)
           print *,nint(tvvon-tvvoff),' seconds in surfsep'
      endif  !  if(vv)

*****************************

      mnmax_bd=0
      bd=lfreeb .and. rbc_bd(0,0).gt.0
       print *,'bd=',bd,' R_00 =',real(rbc_bd(0,0)),' wait for surfsep'
      if(id2.eq.0) close(id1,status='delete')
      if(id2.ne.0) close(id1)
      if(bd) then
        mnmax_bd = mpol_bd*(2*ntor_bd+1)+ntor_bd+1
        mnmax_bd = (mpol_bd+1)*(2*ntor_bd+1)-ntor_bd
        dimnow=  (mpol_bd+1)*(2*ntor_bd+1)-ntor_bd
        if(.not.allocated(xm_bd)) allocate(rmnc_bd(dimnow),
     &   zmns_bd(dimnow), xm_bd(dimnow), xn_bd(dimnow),
     &   rmns_bd(dimnow), zmnc_bd(dimnow))
           rmns_bd=0 ; zmnc_bd=0; rmnc_bd=0 ; zmns_bd=0
           i = 0
           do m=0, mpol_bd
             do n=-ntor_bd, ntor_bd
               if (m.eq.0 .and. n.lt.0) cycle
               i = i+1
               xm_bd(i) = m
               xn_bd(i) = n  !  NOT *nfp
               rmnc_bd(i) = rbc_bd(n,m)
               zmns_bd(i) = zbs_bd(n,m)
             end do
           end do

            call second0(tbdoff)
             call surfsep(mnmax_pl, rmnc_pl, zmns_pl, xm_pl,
     1               xn_pl, mnmax_bd, nu_bd, nv_bd, rmnc_bd, zmns_bd,
     2               xm_bd, xn_bd, nfp, 64, bd_dist, bd_dist_rms,
     2               bd_dist_max, ierr)
           call second0(tbdon)
           print *,nint(tbdon-tbdoff),' seconds in surfsep'
      endif  !  if(bd)

*****************************

      allocate (sqrt_phimod(mnmax*ns), phimod(mnmax*ns), dbzcods(ns),
     1   ffp(ns), ub(ns), iotaf(ns), darea(ns), iotazb(ns), psi(ns),
     2   itors(ns), ipols(ns), szc(ns), szb(ns), dummy(ns),
     3   ixm(mnmax), ixn(mnmax), stat=mn)
      if (mn .ne. 0) stop 'Allocation error in plotout'

      ixm(:mnmax) = nint(xm(:mnmax))
      ixn(:mnmax) = nint(xn(:mnmax))

      do mn = 1, mnmax
         if (ixm(mn).eq.0 .and. ixn(mn).eq.0) mn0 = mn
      enddo

      ntor0 = 1 + ntor
      nrt = ns*nthpts

*****************************
*        COMPUTE FOURIER COEFFICIENTS OF d Bu/ds and d Bs/du
*        ON ZONE BNDRY GRID
*        Bsubu,v on 1/2 grid; dBsubu,v on full grid
*****************************
      hs = 1./(ns - 1)
      phip(1) = phip(2)


      do mn = 1, ntor0
         bmn(mn,1) = 1.5*bmn(mn,2) - 0.5*bmn(mn,3)
      enddo

      do mn = 1, ntor
         gmn(mn,1) = 1.5*gmn(mn,2) - 0.5*gmn(mn,3)
      enddo

      currvmn(1+ntor0:mnmax,1) = 0.
      bmn(1+ntor0:mnmax,1) = 0.
      gmn(1+ntor0:mnmax,1) = 0.
      darea(2:ns) = vp(2:ns)*overr(2:ns)
      overr(2:ns-1) = .5*(overr(2:ns-1)+overr(3:ns))
      iotaf(2:ns-1) = .5*(iotas(2:ns-1)+iotas(3:ns))
      jcuru(1) = 2.*jcuru(2) - jcuru(3)
      overr(1) = 2.*overr(2) - overr(3)
      jcuru(ns) = 2.*jcuru(ns-1) - jcuru(ns-2)
      overr(ns) = 2.*overr(ns-1) - overr(ns-2)
      iotas(1) = 1.5*iotas(2) - 0.5*iotas(3)
      iotaf(1) = iotas(1)
      iotaf(ns) = 1.5*iotas(ns) - 0.5*iotas(ns-1)

*     Forming <J-dot-GradPhi> / < R**-1 > ...
cj      where (overr .ne. 0.) jcurv = jcurv/overr  !  obsolete
      hs = 1.0/(ns - 1)
      ohs = 1.0/hs
      r0 = rmnc(mn0,ns)
      call wrfcn (input_id)
      id2 = pgopen (trim(device))
      if (id2.le.0)
     .  write(*,*)'Fail to Open ',trim(device),' with code ',id2
      if (id2.le.0) stop 'id2'
      call pgscf(2)
      call bd_plotter (r0, explist, lenlist, input_id, trim(device))
      call read_wout_deallocate
      deallocate (sqrt_phimod, phimod, dbzcods, ffp, ub, iotaf,
     1   darea, iotazb, psi, itors, ipols, szc, szb, dummy,
     2   ixm, ixn)
      if (allocated(brcoil)) deallocate (brcoil, plbrfld, brbc,
     1    bzcoil, plbzfld, bzbc)
8033  call pgclos

      end program bd_match
      function polyval(a, na, xv)
        use kind_spec
        implicit none
        real(kind=rprec) :: polyval
        real(kind=rprec), intent(in) :: a(*), xv
        integer :: na ,j
         polyval=xv*0.+a(1)
         do j=2,na
              polyval = polyval+a(j)*xv**(j-1)
         enddo
      end function polyval





      
 
      integer function gen_find_first_in_set (string, set)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      character string*(*), set*(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, j, lstring, lset
C-----------------------------------------------
 
 
 
c     Mimic Str$Find_First_In_Set
 
      gen_find_first_in_set = 0
 
      lstring = len_trim(string)
      lset = len_trim(set)
 
      if (lstring==0 .or. lset==0) return 
 
      do i = 1, lstring
         do j = 1, lset
            if (string(i:i) .eq. set(j:j)) goto 100
         enddo
      enddo
 
      return 
 
  100 continue
 
c     Found a character in SET
 
      gen_find_first_in_set = i
 
      end function !gen_find_first_in_set

      subroutine  get_lun(iu)
      integer*4  iun(20),ist(20)
      save    iun,ist
      data iun/  119,118,117,116,115,114,113,112,111,110,
     &             109,108,107,106,105,104,103,102,101,100/
      data ist/  20*0/
      do j = 1,20
        if (ist(j) .eq. 0) then
          iu = iun(j)
          ist(j) = 1
          return
        endif
      enddo
      iu = -1
      return
      entry free_lun(iu)
      ist(120-iu) = 0
      return
      end subroutine get_lun

        subroutine graf1 (x,y,n,lx,ly,lt,runlbl)
         implicit none
         integer n
         real , dimension(n), intent(in) :: x, y
         real ymin, ymax, xmin, xmax, siz
         character*(*) lx,ly,lt,runlbl
         call pgsave
         call pgbbuf
         call pgsci(1)
         xmin=minval(x);xmax=maxval(x)
         ymin=minval(y);ymax=maxval(y)
         if(ymax.eq.ymin)then
          if(ymax.gt.0.)ymin=.99*ymax
          if(ymax.lt.0.)ymax=.99*ymin
          if(ymax.eq.0.)ymax=1.e-3;ymin=-1.e-3
         endif
         call pgenv(minval(x),maxval(x),ymin,ymax,0,0)
         call pglab(trim(lx),trim(ly),trim(lt))
         call pgqch(siz)
         call pgsch(siz-.5)
         call pgmtxt('B',2.85,0.5,0.5,trim(runlbl))
         call pgsch(siz)
         call pgsci(2)
         call pgslw(3)
         call pgline(n,x,y)
         call pgebuf
         call pgunsa
        return
        end subroutine graf1
        subroutine graf1x (x,y,n,lx,ly,lt,runlbl)
         implicit none
         integer n
         real , dimension(n), intent(in) :: x, y
         real ymin, ymax, xmin, xmax, siz
         character*(*) lx,ly,lt,runlbl
         call pgsave
         call pgbbuf
         call pgsci(1)
         xmin=minval(x);xmax=maxval(x)
         ymin=minval(y);ymax=maxval(y)
         if(ymax.eq.ymin)then
          if(ymax.gt.0.)ymin=.99*ymax
          if(ymax.lt.0.)ymax=.99*ymin
          if(ymax.eq.0.)ymax=1.e-3;ymin=-1.e-3
         endif
         call pgenv(minval(x),maxval(x),ymin,ymax,0,1)
         call pglab(trim(lx),trim(ly),trim(lt))
         call pgqch(siz)
         call pgsch(siz-.5)
         call pgmtxt('B',2.85,0.5,0.5,trim(runlbl))
         call pgsch(siz)
         call pgsci(2)
         call pgslw(3)
         call pgline(n,x,y)
         call pgebuf
         call pgunsa
        return
        end subroutine graf1x

        subroutine graf1pt (x,y,n,lx,ly,lt,runlbl)
         implicit none
         integer n
         real , dimension(n), intent(in)  :: x, y
         real ymin, ymax, siz
         character*(*) lx,ly,lt,runlbl
         call pgsave
         call pgsci(1)
         ymin=minval(y);ymax=maxval(y)
         if(ymax.eq.ymin)then
          if(ymax.gt.0.)ymin=.99*ymax
          if(ymax.lt.0.)ymax=.99*ymin
          if(ymax.eq.0.)ymax=1.e-3;ymin=-1.e-3
         endif
         call pgpage
         call pgvstd
         call pgswin(minval(x),maxval(x),ymin,ymax)
         call pgbox('BCNT',0.,0,'BCNTP1',0.,0)
c         call pgenv(minval(x),maxval(x),ymin,ymax,0,0)
         call pglab(trim(lx),trim(ly),trim(lt))
         call pgqch(siz)
         call pgsch(siz-.5)
         call pgmtxt('B',2.85,0.5,0.5,trim(runlbl))
         call pgsch(siz)
         call pgsci(2)
         call pgpt(n,x,y,2)
         call pgunsa
        return
        end subroutine graf1pt
        subroutine midplpt (x1,y1,nx,ny,lx,ly,lt,runlbl)
         USE Vmagaxis
         implicit none
         integer nx,ny
         real , dimension(nx,ny), intent(in)  :: x1, y1
         real ymin, ymax, siz
         real, dimension(:), allocatable :: x,y
         character*(*) lx,ly,lt,runlbl
         call pgsave
         call pgsci(1)
         allocate(x(2*ny),y(2*ny))
         x=(/x1(1,1:ny),x1((nx+1)/2,1:ny)/)
         y=(/y1(1,1:ny),y1((nx+1)/2,1:ny)/)
         ymin=minval(y);ymax=maxval(y)
         if(ymax.eq.ymin)then
          if(ymax.gt.0.)ymin=.99*ymax
          if(ymax.lt.0.)ymax=.99*ymin
          if(ymax.eq.0.)ymax=1.e-3;ymin=-1.e-3
         endif
         call pgenv(minval(x),maxval(x),ymin,ymax,0,0)
         call pglab(trim(lx),trim(ly),trim(lt))
         call pgqch(siz)
         call pgsch(siz-.5)
         call pgmtxt('B',2.85,0.5,0.5,trim(runlbl))
         call pgsch(siz)
         call pgsci(2)
         call pgpt(2*ny,x,y,1)
         call pgsci(4)
         call pgpt1 ( real(rmagaxis), ymin+(ymax-ymin)/3.,-4)
         call pgunsa
         deallocate(x,y)
        return
        end subroutine midplpt

        subroutine graf2 (x,y,yfit,n,lx,ly,lt,runlbl)
c  two lines on a single scale
         implicit none
         integer i, n
         real ymin, ymax, siz
         real *4, dimension(n) :: x,y,yfit
         character*(*) lx,ly,lt,runlbl
         call pgsave
         call pgbbuf
         call pgsci(1)
         ymin=min(minval(y),minval(yfit))
         ymax=max(maxval(y),maxval(yfit))
         if(ymax.eq.ymin)then
          if(ymax.gt.0.)ymin=.99*ymax
          if(ymax.lt.0.)ymax=.99*ymin
          if(ymax.eq.0.)ymax=1.e-3;ymin=-1.e-3
         endif
         call pgenv(minval(x),maxval(x),ymin,ymax,0,0)
         call pglab(trim(lx),trim(ly),trim(lt))
         call pgqch(siz)
         call pgsch(siz-.5)
         call pgmtxt('B',2.85,0.5,0.5,trim(runlbl))
         call pgsch(siz)
         call pgsci(3)
         call pgpt(n,x,y,-2)
         call pgslw(3)
         call pgsls(1)
         call pgsci(2)
         call pgline(n,x,yfit)
         call pgsci(3)
         call pgpt(n,x,y,2)
         call pgebuf
         call pgunsa
        return
        end subroutine graf2

        subroutine graf2pt (x1,x2,y1,y2,n,lx,ly1,ly2,lt,runlbl)
         integer n
         real , dimension(n), intent(in)  :: x1, x2, y1, y2
         real ymin, ymax, siz
         character*(*) lx,ly1,ly2,lt,runlbl
         call pgsave
         call pgslw(1)
         call pgsci(1)
         xmin=minval(x1);xmax=maxval(x1)
         xmin=min(xmin,minval(x2))
         xmax=max(xmax,maxval(x2))
         ymin=minval(y1);ymax=maxval(y1)
         ymin=min(ymin,minval(y2))
         ymax=max(ymax,maxval(y2))
         if(ymax.eq.ymin)then
          if(ymax.gt.0.)ymin=.99*ymax
          if(ymax.lt.0.)ymax=.99*ymin
          if(ymax.eq.0.)ymax=1.e-3;ymin=-1.e-3
         endif
         call pgpage
         call pgvstd
         call pgswin(xmin,xmax,ymin,ymax)
         call pgbox('BCNT',0.,0,'BCNTP1',0.,0)
c         call pgenv(xmin,xmax,ymin,ymax,0,0)
         call pglab(trim(lx),'(m)',trim(lt))
         call pgqch(siz)
         call pgsch(siz-.5)
         call pgmtxt('B',2.85,0.5,0.5,trim(runlbl))
         call pgsch(siz)
         call pgsci(4)
         call pgslw(1)
         call pgpt(n,x1,y1,-1)
         call pgmtxt('L',2.65,0.,0.,trim(ly1))
         call pgsci(2)
         call pgslw(3)
         call pgpt(n,x2,y2,-2)
         call pgmtxt('L',2.65,1.0,1.0,trim(ly2))
         call pgunsa
        return
        end subroutine graf2pt

        subroutine graf3pt (x,y1,y2,y3,n,lx,ly1,ly2,ly3,lt,runlbl)
         integer n
         real , dimension(n), intent(in)  :: x, y1, y2, y3
         real ymin, ymax, siz
         character*(*) lx,ly1,ly2,ly3,lt,runlbl
         call pgsave
         call pgsci(1)
         ymin=minval(y1);ymax=maxval(y1)
         ymin=min(ymin,minval(y2))
         ymax=max(ymax,maxval(y2))
         ymin=min(ymin,minval(y3))
         ymax=max(ymax,maxval(y3))
         if(ymax.eq.ymin)then
          if(ymax.gt.0.)ymin=.99*ymax
          if(ymax.lt.0.)ymax=.99*ymin
          if(ymax.eq.0.)ymax=1.e-3;ymin=-1.e-3
         endif
         call pgenv(minval(x),maxval(x),ymin,ymax,0,0)
         call pglab(trim(lx),' ',trim(lt))
         call pgqch(siz)
         call pgsch(siz-.5)
         call pgmtxt('B',2.85,0.5,0.5,trim(runlbl))
         call pgsch(siz)
         call pgsci(1)
         call pgpt(n,x,y1,2)
         call pgmtxt('L',2.65,0.,0.,trim(ly1))
         call pgsci(2)
         call pgpt(n,x,y2,3)
         call pgmtxt('L',2.65,0.5,0.5,trim(ly2))
         call pgsci(4)
         call pgpt(n,x,y3,4)
         call pgmtxt('L',2.65,1.,1.,trim(ly3))
         call pgunsa
        return
        end subroutine graf3pt

        subroutine graf2x (x,y,yfit,n,lx,lyl,lyr,lt,runlbl)
c two lines in a box with y axis left and right
         implicit none
         integer i, n
         real ymin, ymax, dx, siz
         real *4, dimension(n) :: x,y,yfit
         character*(*) lx,lyl,lyr,lt,runlbl
         character*8 xopt,yopt
         dx=x(3)-x(1)
         if(n.lt.30)dx=x(2)-x(1)
         call pgsave
         call pgbbuf
         call pgsci(1)
         call pgpage
         call pgvstd
         call pgswin(minval(x),maxval(x),minval(y),maxval(y))
         xopt='BCNST';yopt='B'
         call pgbox(trim(xopt),0.,0,trim(yopt),0.,0)
c         call pgmtxt('B',1.5,0.5,0.5,trim(lx))
         call pglab(trim(lx),' ',trim(lt))
c         call pgmtxt('T',1.5,0.5,0.5,trim(lt))
         call pgqch(siz)
         call pgsch(siz-.5)
         call pgmtxt('B',2.85,0.5,0.5,trim(runlbl))
         call pgsch(siz)
         call pgsci(2)
         yopt='BNST'
         call pgbox(' ',0.,0,trim(yopt),0.,0)
         call pgmtxt('L',2.25,0.5,0.5,trim(lyl))
         call pgsls(1)
         call pgslw(3)
         call pgline(n,x,y)
c         call pgarro(x(n/3)-dx,y(n/3),x(n/3-n/5)-dx,y(n/3))
c         call pgunsa
         xopt=' ';yopt='CMST'
         call pgsci(4)
         call pgswin(minval(x),maxval(x),minval(yfit),maxval(yfit))
         call pgbox(' ',0.,0,trim(yopt),0.,0)
         call pgmtxt('R',2.25,0.5,0.5,trim(lyr))
c         call pgsave
         call pgsls(4)
         call pgslw(8)
         call pgline(n,x,yfit)
c         call pgarro(x(2*n/3)+dx,yfit(2*n/3),x(2*n/3+n/5)+dx,yfit(2*n/3))
         call pgebuf
         call pgunsa
        return
        end subroutine graf2x



        subroutine loggraf1 (x,y,n,lx,ly,lt,runlbl)
         implicit none
         integer i, n
         real ymin, ymax, siz
         real *4, dimension(n) :: x,y
         real *4, dimension(:), allocatable :: logx,logy
         character*(*) lx,ly,lt,runlbl
         character*100 loglabel
         allocate(logy(n))
         logy=alog10(y)
         call pgsave
         call pgbbuf
         call pgsci(1)
         ymin=minval(logy);ymax=maxval(logy)
         if(ymax.eq.ymin)then
          if(ymax.gt.0.)ymin=.99*ymax
          if(ymax.lt.0.)ymax=.99*ymin
          if(ymax.eq.0.)ymax=1.e-3;ymin=-1.e-3
         endif
         call pgenv(minval(x),maxval(x),ymin,ymax,0,20)
         loglabel='Log\d10\u\(2223)'//trim(ly)//'\(2224)'
         call pglab(
     &     trim(lx),trim(loglabel),trim(lt))
         call pgqch(siz)
         call pgsch(siz-.5)
         call pgmtxt('B',2.85,0.5,0.5,trim(runlbl))
         call pgsch(siz)
         call pgslw(3)
         call pgline(n,x,logy)
         call pgebuf
         call pgunsa
         deallocate(logy)
        return
        end subroutine loggraf1


        subroutine loggraf3pt (x,y1,y2,y3,n,lx,ly1,ly2,ly3,lt,runlbl)
         integer n
         real , dimension(n), intent(in)  :: x, y1, y2, y3
         real *4, dimension(:), allocatable :: logx,logy1,logy2,logy3
         real ymin, ymax, siz
         character*(*) lx,ly1,ly2,ly3,lt,runlbl
         character*100 loglabel
         allocate(logy1(n),logy2(n),logy3(n)   )
         logy1=alog10(y1)
         logy2=alog10(y2)
         logy3=alog10(y3)
         call pgsave
         call pgsci(1)
         ymin=minval(logy1);ymax=maxval(logy1)
         ymin=min(ymin,minval(logy2))
         ymax=max(ymax,maxval(logy2))
         ymin=min(ymin,minval(logy3))
         ymax=max(ymax,maxval(logy3))
         if(ymax.eq.ymin)then
          if(ymax.gt.0.)ymin=.99*ymax
          if(ymax.lt.0.)ymax=.99*ymin
          if(ymax.eq.0.)ymax=1.e-3;ymin=-1.e-3
         endif
         call pgenv(minval(x),maxval(x),ymin,ymax,0,20)
         loglabel='Log\d10\u\(2223)'//trim(ly2)//'\(2224)'
         call pglab(trim(lx),' ',trim(lt))
         call pgqch(siz)
         call pgsch(siz-.5)
         call pgmtxt('B',2.85,0.5,0.5,trim(runlbl))
         call pgsch(siz)
         call pgsci(1)
         call pgpt(n,x,logy1,2)
         call pgmtxt('L',2.65,0.,0.,trim(ly1))
         call pgsci(2)
         call pgpt(n,x,logy2,3)
         call pgmtxt('L',2.65,0.5,0.5,trim(loglabel))
         call pgsci(4)
         call pgpt(n,x,logy3,4)
         call pgmtxt('L',2.65,1.,1.,trim(ly3))
         call pgunsa
        return
        end subroutine loggraf3pt
        subroutine loggrafpt (x,y,n,lx,ly,lt,runlbl)
         implicit none
         integer i, n
         real ymin, ymax, dx, siz
         real *4, dimension(n) :: x,y
         real *4, dimension(:), allocatable :: logx,logy
         character*(*) lx,ly,lt,runlbl
         character*100 loglabel
         allocate(logy(n))
         logy=alog10(y)
         call pgsave
         call pgbbuf
         call pgsci(1)
         ymin=minval(logy);ymax=maxval(logy)
         if(ymax.eq.ymin)then
          if(ymax.gt.0.)ymin=.99*ymax
          if(ymax.lt.0.)ymax=.99*ymin
          if(ymax.eq.0.)ymax=1.e-3;ymin=-1.e-3
         endif
         call pgenv(minval(x),maxval(x),ymin,ymax,0,20)
         loglabel='Log\d10\u\(2223)'//trim(ly)//'\(2224)'
         call pglab(
     &     trim(lx),trim(loglabel),trim(lt))
         call pgqch(siz)
         call pgsch(siz-.5)
         call pgmtxt('B',2.85,0.5,0.5,trim(runlbl))
         call pgsch(siz)
         call pgsci(2)
         call pgpt(n,x,logy,2)
         call pgsci(15)
         call pgslw(1)
         call pgline(n,x,logy)
         call pgebuf
         call pgunsa
         deallocate(logy)
        return
        end subroutine loggrafpt

        subroutine logint(contr,zmin,zmax,n,irr)
        implicit none
c logarithmic interpolation
        real zmin,zmax,zlmin,zlmax,amin, amax
        logical misign,allpos,zerol,zerou
        complex c,d,bmin,bmax
        real cs, conv, pi, el, da
        integer i,n,n1, irr
        real*4 contr(*)
        pi=4.*atan(1.)
        if(zmax-zmin .lt. 10.)goto 1001  !  linear
        zlmax=zmax
        zlmin=zmin
        zerol=zmin.eq.0.
        zerou=zmax.eq.0.
        if(zerol)zlmin=-zmax/1.e4
        if(zerou)zlmax=abs(zmin/1.e4)
        conv=log(10.)
        allpos=(zlmax.gt.0.).and.(zlmin.gt.0.)
        misign=zlmin.lt.0.
        bmin=cmplx(zlmin,.0)
        bmin=(log(bmin))
        if(misign)cs=real(bmin)
        if(misign)bmin=cmplx(cs,pi)
        amin=real(bmin)
        if(zerol)amin=0
        if(misign)amin=-amin
        misign=zlmax.lt.0.
        bmax=cmplx(zlmax,.0)
        bmax=(log(bmax))
        if(misign)cs=real(bmax)
        if(misign)bmax=cmplx(cs,pi)
        amax=real(bmax)
        if(zerou)amax=0
        if(amin.gt.amax)then
          cs=amin
          amin=amax
          amax=cs
        endif
        da=real(amax-amin)/(n-1)
        if(allpos) then
          do i=1,n
          c=cmplx(amin+(i-1)*da)
          contr(i)=real(exp(c)) 
          enddo
        else
          el=abs(zlmax)+abs(zlmin)
          n1=-(zlmin/el)*n
          n1=max(1,n1)
          do i=1,n1
          c=cmplx(amin+(i-1)*da,pi)
          contr(i)=real(exp(-c)) 
          enddo
          i=n1+1
          contr(i)=0
          do i=n1+2,n
          c=cmplx(amin+(i-1)*da)
          contr(i)=real(exp(c)) 
          enddo
        endif
        irr=0
        return

 1001   continue     
        da=real(zmax-zmin)/(n-1)
        do i=1,n
           contr(i)=zmin +da*(i-1)
        enddo
        irr= 1001
        return
        end subroutine logint

      subroutine bd_plotter(r0, explist, lenlist, input_id, device)
C-----------------------------------------------
CSave as working, has axes3d, offset corrected
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use read_wout_mod, itor_w=>itor
      use boundaries
      use Vpname2
      USE Vpltcn2
      USE Vpltcn6
      USE Vrzarray
      USE Vplotdata
      USE Vmagaxis
      USE Vindat2
      use system_mod, only: system
      use vmec_input, only: rbc,rbs,zbc,zbs, lfreeb
      use safe_open_mod
      use optim, only: rbc_vv, zbs_vv, mpol_vv, ntor_vv,
     &    rbc_bd, zbs_bd, mpol_bd, ntor_bd, bigno,
     1    vv_dist,  vv_dist_rms,  vv_dist_max, sigma_iota,
     2    nu_vv, nv_vv,target_rbtor,target_vv, 
     1    bd_dist, bd_dist_rms, bd_dist_max,target_iota,
     2    nu_bd, nv_bd, target_beta, target_cur

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer lenlist
      real(kind=rprec) :: r0
      character explist*(*), input_id*(*), device*(*)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(kind=rprec) :: zero = 0.0_dp , one = 1.0_dp
      integer, parameter :: ntdu = 31, ntdv = 21
      integer, parameter :: nplot1 = 25, ntheta1 = 2, 
     1  nox = 1, noy = 3, lx = -1, ly = 0, iend = 3,
     2  ilog = 1, icart = 0, oneppg = 2, fourppg = 1
      integer, dimension(4) :: ivar = (/0, 1, 2, 3/),
     1  ncon = (/10, 25, 15, 20/)      
      real(kind=rprec), parameter :: grsize1 = 3.5, grsize2 = 2.5, 
     1  grfull = 7.5
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: id1, id2, np, id, pgopen, nzeta, ntheta, ntheta2
      integer :: i, ii, nplots, les, kz, k, kt, j, l, js, kz1, nplt1,
     1   ndeg, ns1, mn, mp, nt, ichar, istat, nznt, nznt1, ipag
      integer :: nsurf , ntht , nphi ,jxbcode
      real(kind=rprec) :: small 
      real(kind=rprec), dimension(ns) :: br, bz, phin
      real(kind=rprec), dimension(2*ns) :: oqmid
      real(kind=rprec), dimension(100) :: time
      real(kind=rprec), dimension(nplot1) :: raxis_v, zaxis_v
      real(kind=rprec), dimension(nfloops + 1) :: delflm
      real(kind=rprec), dimension(nbloops) :: delbc
      real(kind=rprec), dimension(nfloops + 1) :: indflm
      real(kind=rprec), dimension(nbloops) :: indbc
      real(kind=rprec), dimension(:,:), allocatable :: dummy1
      real(kind=rprec), dimension(:), allocatable :: rbndy,zbndy
      real(kind=rprec), dimension(:), allocatable :: modb, sqrt_phim, 
     1  phim, gsqrt, torcur, dummy2, r12, z12, ru12, zu12
      real(kind=rprec) :: dth, denom, bdotgradv0, phiangle, offset,
     1   bdotgradvn, phinorm, t1, a0, an, cosphi, sinphi
      character :: page*10, pagedesc*50, nchar*100, mchar*100, 
     1   pchar*100, ititlet*100, fname*80, fdum*80, lt*100,
     1   xlabel*100, ylabel*100, zlabel*100, mtitle*40
      character ::device_type*11
      real xspread,xxspread,yspread,yyspread,deltx,siz,
     1 xmin,xmax,ymin,ymax,yzero,ymmax,ymmin,xmmin,xmmax
      integer maxc,ic,nmarks,nspread
      real(kind=rprec), dimension(:), allocatable :: hiota, presb, 
     1   beta_volb, phipb,  phib, bvcob,  bucob , xmb, xnb
      integer , dimension(:), allocatable :: jlistb
      real(kind=rprec), dimension(:,:), allocatable :: 
     1   bmnb, rmnb, zmnb, pmnb
      logical first
      data first/.true./
      real(kind=rprec) :: aspectb, rmax_surfb, rmin_surfb, betaxisb
      integer :: nsb, nfpb, nboz, mboz, mnboz, versionb, 
     1 iunit, ierr, iread, jrad, mb, nb

      integer :: nx, ny
      integer, dimension(:,:), allocatable :: lka
      real(kind=rprec), dimension(:,:), allocatable :: 
     1  bsupu, bsupv, bsubu, bsubv, bsubs, residual12, xfres12
      real(kind=rprec), dimension(:,:), allocatable :: 
     1  jsupu, jsupv, jcrossb, pprime, bdotj, jsups
      real(kind=rprec), dimension(:), allocatable :: pang, tang,
     1  avforce, hshere, amaxfor, aminfor
      real, dimension(:), allocatable :: xplot,yplot
      logical, dimension(:), allocatable :: iplot
      real x1,x2,y1,y2
      integer itheta,izeta
      logical addon, alldone
      real(kind=rprec) :: polyval

C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
C-----------------------------------------------
      ipag=0


      mn = mnmax*ns
      allocate (dummy1(mnmax,ns), dummy2(nrt), r(nrt), z(nrt), 
     1   ru(nrt), zu(nrt), modb(nrt), sqrt_phim(nrt), phim(nrt), 
     2   gsqrt(nrt),torcur(nrt), dbsubudstz(nrt), bsubutz(nrt), 
     3   bsubvtz(nrt),r12(nrt), z12(nrt), ru12(nrt), zu12(nrt), 
     4   stat = istat)
      if (istat .ne. 0) stop 'allocation error in plotter'
      allocate(rbndy(nthpts),zbndy(nthpts), stat = istat)
      if (istat .ne. 0) stop 'allocation error in plotter'


      dummy1 = 0.
      dummy2 = 0.
 
      noxc = 0
      noyc = noy - 1
      lxc = lx
      lyc = ly
 
*
*                 SETUP PGPLOT GRAPHICS
*
      call pgpap(7.5,1.)
      call pgask(.true.)
      call pgsch(1.85)
      call pgslw(3)
      call pgsubp(1,1)

 
*
*                 COMPUTE R,Z SCALES
*

      nplots = 1
      if (ntor .ne. 0) nplots = 4
      les = 1 + nthpts*(ns - 1)
      do kz = 1, nplots
         call totz(nthpts,ns,nplots,kz,r,z,rmnc,zmns,rmns,zmnc)
         if (kz .eq. 1) then
            rmax = r(les)
            rmin = r(les)
            zmax = z(les)
            zmin = z(les)
         endif
         rmax = max(rmax,maxval(r(les:nrt)))
         rmin = min(rmin,minval(r(les:nrt)))
         zmax = max(zmax,maxval(z(les:nrt)))
         zmin = min(zmin,minval(z(les:nrt)))
      enddo

      call pgsubp(2,2)
      if(allocated(rbndy))deallocate(rbndy,zbndy)
      if(allocated(ru12))deallocate(ru12,zu12)
      if(allocated(r12))deallocate(r12,z12)
      allocate(rbndy(nthpts),zbndy(nthpts),ru12(nthpts),zu12(nthpts),
     & r12(nthpts),z12(nthpts), stat = istat)
      if (istat .ne. 0) stop 'allocation error for boundary figure'
      do kz = 1, nplots
         call totb2(nthpts, ns, nplots, kz,
     .    rbndy, zbndy, rmnc_pl ,zmns_pl ,rmns_pl ,zmnc_pl ,
     .    nfp, mnmax_pl, xm_pl, xn_pl,  mpol_pl+1, ntor_pl )
         if(vv)call totb2(nthpts, ns, nplots, kz,
     .    r12, z12, rmnc_vv ,zmns_vv ,rmns_vv ,zmnc_vv ,
     .    nfp, mnmax_vv, xm_vv, xn_vv,  mpol_vv+1, ntor_vv )
         if(bd)call totb2(nthpts, ns, nplots, kz,
     .    ru12, zu12, rmnc_bd, zmns_bd, rmns_bd, zmnc_bd,
     .    nfp, mnmax_bd, xm_bd, xn_bd,  mpol_bd+1, ntor_bd )
        rmin=minval(rbndy)
        rmax=maxval(rbndy)
        zmin=minval(zbndy)
        zmax=maxval(zbndy)
        if(vv)rmin=min(rmin,minval(r12))
        if(vv)rmax=max(rmax,maxval(r12))
        if(vv)zmin=min(zmin,minval(z12))
        if(vv)zmax=max(zmax,maxval(z12))
        if(bd)rmin=min(rmin,minval(ru12))
        if(bd)rmax=max(rmax,maxval(ru12))
        if(bd)zmin=min(zmin,minval(zu12))
        if(bd)zmax=max(zmax,maxval(zu12))
      call pgsci(1)
      call  pgenv (real(rmin), real(rmax), real(zmin), real(zmax), 1, 0)      
           xlabel='R (m)'
           ylabel='Z (m)'
           call pgsci(3)
           zlabel= 'Boundary(g), '
           call pgsci(4)
           zlabel=trim(zlabel)//' target(b)'
           call pgsci(2)
           zlabel=trim(zlabel)//', VV(r)'
           call pgsci(1)
      call pglab(trim(xlabel),trim(ylabel),trim(zlabel))
      write(zlabel,'(f4.0)')real(720.d0/(nfp*nfp*nplots)*(kz-1))
      lt='N\dfp\u\(0647)= '//adjustl(zlabel)
      call pgmtxt('T',1.0,0.5,0.5,trim(lt))
      lt='vv_dist_min, vv_dist_rms, vv_dist_max (in mm):'
      write(zlabel,'(3(x,i4))')nint(1000*vv_dist), 
     &nint(1000*vv_dist_rms),nint(1000*abs(vv_dist_max))
      lt=trim(lt)//adjustl(zlabel)
      if(kz==3)print*, lt
      call pgsci(2)
      if(kz==3)call pgmtxt('T',3.15,0.,0.,trim(lt))
      lt='bd_dist_min, bd_dist_rms, bd_dist_max (in mm):'
      write(zlabel,'(3(x,i4))')nint(1000*bd_dist), 
     &nint(1000*bd_dist_rms),nint(1000*abs(bd_dist_max))
      lt=trim(lt)//adjustl(zlabel)
      if(kz==1)print*, lt
      call pgsci(4)
      if(kz==1)call pgmtxt('T',3.15,0.,0.,trim(lt))
      call pgqch(siz)
      call pgsch(siz-.5)
      call pgmtxt('B',3.25,0.5,0.5,trim(runlabel))
      call pgsch(siz)
      
c         Plot the Boundary -- Last Closed Flux Surface
      call  pgscf (2  )                                               
      call pgsci(3)
      call pgslw(1)
      call pgpt(nthpts,real(rbndy),real(zbndy),5)
      call pgsci(2)
      call pgslw(3)
      if(vv)call pgline(nthpts,real(r12),real(z12))
      call pgsci(4)
      call pgslw(4)
      if(bd)call pgline(nthpts,real(ru12),real(zu12))

      enddo  !kz  
      call newframe
*
*                 COMPUTE EQUILIBRIUM PROFILES
*
       phin(:ns) = abs(phi(:ns))
 

            ndata = 0
         jcurv(1) = 2.*jcurv(2) - jcurv(3)
         zlabel='FLUX-AV. JTOR (A/m\u2\d)'
         ylabel= '<J\.\(2266)\(0647)>\(2770)<R\u-1\d>'
         xlabel='\gF'
         call graf1
     &     (real(phin(1:)),real(jcurv(1:)),
     &     ns,xlabel,ylabel,zlabel,runlabel)
         xlabel="('R\(0210)B\dt\u/R=',f5.2,1hT,' A=',f5.2)"
         write(ylabel,fmt=xlabel)RBtor/ Rmajor,  Rmajor/ Aminor
         call pgmtxt('T',0.9,0.5,0.5,trim(ylabel))
      lt='TARGET_BETA (%), TARGET_RBTOR (m-T) : '
      write(zlabel,'(2(3x,f4.2))')100.*target_beta, target_rbtor
      lt=trim(lt)//(zlabel)
      call pgsci(2)
      print*,trim(lt)
      call pgmtxt('T',3.175,0.,0.,trim(lt))
      call pgsci(1)


         zlabel='PRESSURE'
         ylabel= 'P(\gF)'
         xlabel='\gF'
         call graf1
     &     (real(phin(1:)),real(pres(1:)),
     &     ns,xlabel,ylabel,zlabel,runlabel)
         xlabel=
     &  "('R\(0210)B\d1\u=',f5.2,x,'R\(0210)B\d0\u=',f5.2,4h m-T)"
         write(ylabel,fmt=xlabel)RBtor, RBtor0
         call pgmtxt('T',0.9,0.5,0.5,trim(ylabel))
      addon=.false.
      do i=1,49
       if(sigma_iota(i) .lt. bigno) addon = .true.
      enddo
      zlabel='TRANSFORM'
      ylabel= '\gi'
      xlabel='\gF'
      if (.not. addon ) call graf1
     &     (real(phin(1:)),real(iotaf(1:)),
     &     ns,xlabel,ylabel,zlabel,runlabel)
         xlabel="('<\gb>=',f5.2,1h%,x,'<B>=',f5.2,2h T)"
      if (addon ) then
       dummy2=0.
       do i=1,ns
        dummy2(i) = polyval(target_iota,11,phin(i))
       enddo
       call graf2
     &     (real(phin(1:)),real(dummy2(1:)),real(iotaf(1:)),
     &     ns,xlabel,ylabel,zlabel,runlabel)
         xlabel="('<\gb>=',f5.2,1h%,x,'<B>=',f5.2,2h T)"
      endif
         write(ylabel,fmt=xlabel)wp/wb*100.,volavgb
         call pgmtxt('T',0.9,0.5,0.5,trim(ylabel))
      lt='TARGET_CUR, TARGET_VV (mm) : '
      write(zlabel,'(1(5x,1pe10.3,5x,i4))')
     &  target_cur,nint(1000*target_vv)
      lt=trim(lt)//(zlabel)
      call pgsci(4)
      print*,trim(lt)
      call pgmtxt('T',3.175,0.,0.,trim(lt))
      call pgsci(1)
 
         ndata = 0
         zlabel='TOROIDAL CURRENT'
         ylabel= 'I\dTOR\u'
         xlabel='\gF'
         call graf1
     &     (real(phin(1:)),real(itors(1:)),
     &     ns,xlabel,ylabel,zlabel,runlabel)
         xlabel="('Plasma Current=',1pe10.3,' A')"
         write(ylabel,fmt=xlabel)itors(size(itors))
         call pgmtxt('T',0.9,0.5,0.5,trim(ylabel))


      if(allocated( dummy1))deallocate( dummy1)
      if(allocated( dummy2))deallocate( dummy2)
      if(allocated(  r))deallocate(  r)
      if(allocated( z))deallocate( z)
      if(allocated( ru))deallocate( ru)
      if(allocated( modb))deallocate( modb)
      if(allocated( sqrt_phim))deallocate( sqrt_phim)
      if(allocated( phim))deallocate( phim)
      if(allocated( gsqrt))deallocate( gsqrt)
      if(allocated( torcur))deallocate( torcur)
      if(allocated( dbsubudstz))deallocate( dbsubudstz)
      if(allocated( bsubutz))deallocate( bsubutz)
      if(allocated( bsubvtz))deallocate( bsubvtz)
      if(allocated(  r12))deallocate(  r12)
      if(allocated( z12))deallocate( z12)
      if(allocated( ru12))deallocate( ru12)
      if(allocated( zu12))deallocate( zu12)

 8033 return

  150 format(a)
      end subroutine bd_plotter

      subroutine splaan(n, x, y, b, c, d)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer n
      real(kind=rprec), dimension(n) :: x, y, b, c, d
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nm1, ib, i
      real(kind=rprec) :: t
C-----------------------------------------------
c
c  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
c  for a cubic interpolating spline for which s-prime(x1)=0.
c
c    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
c
c    for  x(i) .le. x .le. x(i+1)
c
c  input..
c
c    n = the number of data points or knots (n.ge.2)
c    x = the abscissas of the knots in strictly increasing order
c    y = the ordinates of the knots
c
c  output..
c
c    b, c, d  = arrays of spline coefficients as defined above.
c
c  using  p  to denote differentiation,
c
c    y(i) = s(x(i))
c    b(i) = sp(x(i))
c    c(i) = spp(x(i))/2
c    d(i) = sppp(x(i))/6  (derivative from the right)
c
ccccccccccccccc
c  the accompanying function subprogram  seval  can be used
c  to evaluate the spline.
c
c
c
      nm1 = n - 1
      if (n < 2) return 
      if (n >= 3) then
c
c  set up tridiagonal system
c
c  b = diagonal, d = offdiagonal, c = right hand side.
c
         d(1) = x(2) - x(1)
         c(2) = (y(2)-y(1))/d(1)
         d(2:nm1) = x(3:nm1+1) - x(2:nm1)
         b(2:nm1) = 2.*(d(:nm1-1)+d(2:nm1))
         c(3:nm1+1) = (y(3:nm1+1)-y(2:nm1))/d(2:nm1)
         c(2:nm1) = c(3:nm1+1) - c(2:nm1)
c
c  end conditions.  third derivatives at  x(1)  and  x(n)
c  obtained from divided differences
c
         b(1) = 2.*d(1)
         b(n) = -d(n-1)
         c(1) = 0.
         c(n) = 0.
         if (n .ne. 3) then
            c(1) = (y(2)-y(1))/d(1)
            c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
            c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
c
c  forward elimination
c
         endif
         do i = 2, n
            t = d(i-1)/b(i-1)
            b(i) = b(i) - t*d(i-1)
            c(i) = c(i) - t*c(i-1)
         enddo
c
c  back substitution
c
         c(n) = c(n)/b(n)
         do ib = 1, nm1
            i = n - ib
            c(i) = (c(i)-d(i)*c(i+1))/b(i)
         enddo
c
c  c(i) is now the sigma(i) of the text
c
c  compute polynomial coefficients
c
         b(n) = (y(n)-y(nm1))/d(nm1) + d(nm1)*(c(nm1)+2.*c(n))
         b(:nm1) = (y(2:nm1+1)-y(:nm1))/d(:nm1) - d(:nm1)*(c(2:nm1+1)+2.
     1      *c(:nm1))
         d(:nm1) = (c(2:nm1+1)-c(:nm1))/d(:nm1)
         c(:nm1) = 3.*c(:nm1)
         c(n) = 3.*c(n)
         d(n) = d(n-1)
         return 
c
      endif
      b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.
      d(1) = 0.
      b(2) = b(1)
      c(2) = 0.
      d(2) = 0.

      end subroutine splaan
 
      subroutine spline(n, x, y, b, c, d)
c       the codes (spline & seval) are taken from:
c       forsythe,malcolm and moler,
c       "computer methods for mathematical computations",
c       prentice-hall, 1977.
c
c       the codes (spleen,splaan & speval) are adaptations
c       by r.m. wieland for special cases ... see comments
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer n
      real(kind=rprec), dimension(n) :: x, y, b, c, d
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nm1, ib, i
      real(kind=rprec) :: t
C-----------------------------------------------
c
c  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
c  for a cubic interpolating spline
c
c    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
c
c    for  x(i) .le. x .le. x(i+1)
c
c  input..
c
c    n = the number of data points or knots (n.ge.2)
c    x = the abscissas of the knots in strictly increasing order
c    y = the ordinates of the knots
c
c  output..
c
c    b, c, d  = arrays of spline coefficients as defined above.
c
c  using  p  to denote differentiation,
c
c    y(i) = s(x(i))
c    b(i) = sp(x(i))
c    c(i) = spp(x(i))/2
c    d(i) = sppp(x(i))/6  (derivative from the right)
c
ccccccccccccccc
c  the accompanying function subprogram  seval  can be used
c  to evaluate the spline.
c
c
c
      nm1 = n - 1
      if (n < 2) return 
      if (n >= 3) then
c
c  set up tridiagonal system
c
c  b = diagonal, d = offdiagonal, c = right hand side.
c
         d(1) = x(2) - x(1)
         c(2) = (y(2)-y(1))/d(1)
         d(2:nm1) = x(3:nm1+1) - x(2:nm1)
         b(2:nm1) = 2.*(d(:nm1-1)+d(2:nm1))
         c(3:nm1+1) = (y(3:nm1+1)-y(2:nm1))/d(2:nm1)
         c(2:nm1) = c(3:nm1+1) - c(2:nm1)
c
c  end conditions.  third derivatives at  x(1)  and  x(n)
c  obtained from divided differences
c
         b(1) = -d(1)
         b(n) = -d(n-1)
         c(1) = 0.
         c(n) = 0.
         if (n .ne. 3) then
            c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
            c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
            c(1) = c(1)*d(1)**2/(x(4)-x(1))
            c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
c
c  forward elimination
c
         endif
         do i = 2, n
            t = d(i-1)/b(i-1)
            b(i) = b(i) - t*d(i-1)
            c(i) = c(i) - t*c(i-1)
         enddo
c
c  back substitution
c
         c(n) = c(n)/b(n)
         do ib = 1, nm1
            i = n - ib
            c(i) = (c(i)-d(i)*c(i+1))/b(i)
         enddo
c
c  c(i) is now the sigma(i) of the text
c
c  compute polynomial coefficients
c
         b(n) = (y(n)-y(nm1))/d(nm1) + d(nm1)*(c(nm1)+2.*c(n))
         b(:nm1) = (y(2:nm1+1)-y(:nm1))/d(:nm1) - d(:nm1)*(c(2:nm1+1) + 
     1      2.d0*c(:nm1))
         d(:nm1) = (c(2:nm1+1)-c(:nm1))/d(:nm1)
         c(:nm1) = 3.d0*c(:nm1)
         c(n) = 3.d0*c(n)
         d(n) = d(n-1)
         return 
c
      endif
      b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.
      d(1) = 0.
      b(2) = b(1)
      c(2) = 0.
      d(2) = 0.

      end subroutine spline

      subroutine splint(xa, ya, y2a, n, x, y, yp, ndim)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer n, ndim
      real(kind=rprec), dimension(*) :: xa, ya, y2a, x, y, yp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, klo, khi, k
      real(kind=rprec) :: c1o6, deriv, a, b, h, h2, a2, b2, y26lo, y26hi
C-----------------------------------------------
*
*       SPLINE INTERPOLATION ROUTINE (Numerical Recipes, pg. 89)
*       XA: ordered array of length N of ordinates at which function YA=F(XA)
*           is tabulated
*       YA: array of length N , = F(XA)
*       Y2A: array of second derivatives at XA points
*       computed from call to SPLINE
*       X : value at which Y = F(X) is to be computed from splines
*       YP = dY/dX at X
*       NDIM: dimension of X, Y, YP arrays
 
 
      c1o6 = 1.0/6.0
      deriv = yp(1)
      klo = 1
      khi = n
      do i = 1, ndim
 
         do while(khi - klo > 1)
            k = (khi + klo)/2
            if (xa(k) > x(i)) then
               khi = k
            else
               klo = k
            endif
         enddo
 
         h = xa(khi) - xa(klo)
         a = xa(khi) - x(i)
         b = x(i) - xa(klo)
         h2 = h*h
         a2 = a*a
         b2 = b*b
         y26lo = c1o6*y2a(klo)
         y26hi = c1o6*y2a(khi)
         y(i) = (a*(ya(klo)+(a2-h2)*y26lo)+b*(ya(khi)+(b2-h2)*y26hi))/h
         if (deriv .ne. 0.0) yp(i) = (ya(khi)-ya(klo)+y26hi*(3.0*b2-h2)-
     1      y26lo*(3.0*a2-h2))/h
         if (i<ndim .and. x(i+1)>x(i)) then
            khi = n
         else
            klo = 1
         endif
      enddo
      end subroutine splint

      subroutine str_strip(string, inblen)
C
C...       This routine compresses a character string, deleting all blanks
c...       the compressed string is returned in the pkg it came in [STRING]],
c...       together with a label [INBLEN] telling how long it is. (Dick
c
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer inblen
      character string*(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ilen, ind, j
C-----------------------------------------------
c
 
      ilen = len(string)
      inblen = 0
      if (ilen .eq. 0) return 
c
      ind = 0
      do j = 1, ilen
         if (string(j:j) .ne. ' ') then
            ind = ind + 1
            string(ind:ind) = string(j:j)
         endif
      enddo
      if (ind < ilen) string(ind+1:) = ' '
 
      inblen = ind
 
      end subroutine str_strip


      subroutine totz(ntheta, ns1, nplots, kz, r, z, 
     1   rmnc, zmns, rmns, zmnc)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use read_wout_mod, 
     .only: nfp, mnmax, xm, xn, rprec, mpol, ntor
      use Vpname1, only: nrt
      use Vindat2, only: twopi
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ntheta, ns1, nplots, kz
      real(kind=rprec), dimension(*), intent(out) :: r, z
      real(kind=rprec), dimension(*), intent(in) :: 
     1   rmnc, zmns, rmns, zmnc

C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nrth, l, j, jes, mes, mn, kt, m, iarg
      real(kind=rprec) :: pit, piz, xm0, xn0, arg
C-----------------------------------------------
      pit = 1.d0/(ntheta - 1)
      piz=2.d0/(nfp*nfp*nplots)
      nrth = ns1*ntheta
      if (nrth .gt. nrt) stop 'nrth > nrt in totz'
      r(:nrth) = 0.d0
      z(:nrth) = 0.d0
      do j = 1, ns1
         jes = ntheta*(j - 1)
         mes = mnmax*(j - 1)
         do mn = 1, mnmax
            xm0 = xm(mn)*pit
            xn0 = xn(mn)*piz
            do kt = 1, ntheta
               l = kt + jes
               m = mn + mes
               arg = xm0*(kt - 1) - xn0*(kz - 1)
               iarg = arg
               arg = twopi*(arg - iarg)
               r(l) = r(l) + rmnc(m)*cos(arg) + rmns(m)*sin(arg)
               z(l) = z(l) + zmns(m)*sin(arg) + zmnc(m)*cos(arg)
            enddo
         enddo
      enddo
      end subroutine totz

      subroutine totzu(ntheta, ns1, nplots, kz, r, z, 
     1   rmnc, zmns, rmns, zmnc)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use read_wout_mod, only: rprec, nfp, mnmax, xm, xn
      use Vindat2, only: twopi
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ntheta, ns1, nplots, kz
      real(kind=rprec), dimension(*), intent(out) :: r, z
      real(kind=rprec), dimension(*), intent(in) :: 
     1  rmnc, zmns, rmns, zmnc
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nrth, l, j, jes, mes, mn, kt, m, iarg
      real(kind=rprec) :: pit, piz, xm0, xn0, arg
C-----------------------------------------------
      pit = 1.d0/(ntheta - 1)
      piz=2.d0/(nfp*nfp*nplots)
      nrth = ns1*ntheta
      r(:nrth) = 0.
      z(:nrth) = 0.
      do j = 1, ns1
         jes = ntheta*(j - 1)
         mes = mnmax*(j - 1)
         do mn = 1, mnmax
            xm0 = xm(mn)*pit
            xn0 = xn(mn)*piz
            do kt = 1, ntheta
               l = kt + jes
               m = mn + mes
               arg = xm0*(kt - 1) - xn0*(kz - 1)
               iarg = arg
               arg = twopi*(arg - iarg)
               r(l) =r(l) + xm(mn)*(rmns(m)*cos(arg) - rmnc(m)*sin(arg))
               z(l) =z(l) - xm(mn)*(zmnc(m)*sin(arg) - zmns(m)*cos(arg))
            enddo
         enddo
      enddo
      end subroutine totzu

      subroutine totb2(ntheta, ns1, nplots, kz, r, z, 
     1   rmnc, zmns, rmns, zmnc,
     2   nfp, mnmax, xm, xn,  mpol, ntor)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use read_wout_mod,only: rprec
      use Vindat2, only: twopi
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ntheta, ns1, nplots, kz
      real(kind=rprec), dimension(*), intent(out) :: r, z
      integer, intent(in) :: nfp, mnmax, mpol, ntor
      real(kind=rprec), dimension(*), intent(in) :: 
     1   rmnc, zmns, rmns, zmnc, xm, xn

C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nrth, l, j, jes, mes, mn, kt, m, iarg
      real(kind=rprec) :: pit, piz, xm0, xn0, arg
C-----------------------------------------------
      pit = 1.d0/(ntheta - 1)
      piz=2.d0/(nfp*nfp*nplots)
      piz=2.d0/(nfp*nplots)  !  xn does not have NFP here
      nrth = 1*ntheta
      r(1:nrth) = 0.d0
      z(1:nrth) = 0.d0
      do kt=1,ntheta
       do l=1,mnmax
          arg=twopi*(xm(l)*pit*(kt-1)-xn(l)*piz*(kz-1))
          r(kt)=r(kt)+rmnc(l)*cos(arg)
          z(kt)=z(kt)+zmns(l)*sin(arg)
       enddo
      enddo
      end subroutine totb2
      subroutine totb(ntheta, ns1, nplots, kz, r, z, 
     1   rmnc, zmns, rmns, zmnc)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use read_wout_mod, 
     .only: nfp, mnmax, xm, xn, rprec, mpol, ntor
      use Vpname1, only: nrt
      use Vindat2, only: twopi
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ntheta, ns1, nplots, kz
      real(kind=rprec), dimension(*), intent(out) :: r, z
      real(kind=rprec), dimension(*), intent(in) :: 
     1   rmnc, zmns, rmns, zmnc

C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nrth, l, j, jes, mes, mn, kt, m, iarg
      real(kind=rprec) :: pit, piz, xm0, xn0, arg
C-----------------------------------------------
      pit = 1.d0/(ntheta - 1)
      piz=2.d0/(nfp*nfp*nplots)
      piz=2.d0/(nfp*nfp*nplots)
      nrth = 1*ntheta
      if (nrth .gt. nrt) stop 'nrth > nrt in totz'
      r(:nrth) = 0.d0
      z(:nrth) = 0.d0
      j = ns1
         jes = ntheta*(j - 1)
         mes = mnmax*(j - 1)
         do mn = 1, mnmax
            xm0 = xm(mn)*pit
            xn0 = xn(mn)*piz
            do kt = 1, ntheta
               l = kt + jes
               m = mn + mes
               arg = xm0*(kt - 1) - xn0*(kz - 1)
               iarg = arg
               arg = twopi*(arg - iarg)
               r(kt) = r(kt) + rmnc(m)*cos(arg) + rmns(m)*sin(arg)
               z(kt) = z(kt) + zmns(m)*sin(arg) + zmnc(m)*cos(arg)
            enddo
         enddo
      end subroutine totb


      subroutine wrfcn(input_id)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use read_wout_mod, itor_w=>itor
      use Vpname1
      use Vpname2
      use Vindat2
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      character input_id*(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: iloop, nmax1, mn, n1, j, jp, i
      real(kind=rprec) :: es, tb1, tm1, tv1, tp1, ti1, ub1
      character, dimension(3) :: ichar1*5, ichar2*5
C-----------------------------------------------
 
      data ichar1/'RmnC(', 'ZmnS(', 'LmnS('/
      data ichar2/'RmnS(', 'ZmnC(', 'LmnC('/
      twopi = 8.*atan(1.0)
*
*                 PRINT SYMMETRIC TERMS -- RmnC and ZmnS
*                 INTERPOLATE LAMBDA ONTO FULL MESH FIRST
*

      do mn = 1, mnmax
         if (ixm(mn) .ne. 0) then
            lmns(mn,1) = 0.
         else
            lmns(mn,1) = 1.5*lmns(mn,2) - 0.5*lmns(mn,3)
         endif
      enddo
               
      do j = 2, ns-1
         do mn = 1,mnmax
            lmns(mn,j) = 0.5*(lmns(mn,j) + lmns(mn,j+1))
         enddo         
      enddo   
      
      do mn = 1,mnmax
         lmns(mn,ns) = 2.0*lmns(mn,ns-1) - lmns(mn,ns-2)
      enddo      

      write (39, 200)
      do iloop = 1, 3
         nmax1 = 5
         do mn = 1, mnmax, 6
            if (mn > mnmax - 6) nmax1 = mnmax - mn
            write (39, 210) (ichar1(iloop),ixm(mn+n1),ixn(mn+n1),n1=0,
     1         nmax1)
            do j = 1, ns
               es = (j - 1)*hs
               select case (iloop) 
               case default
                  write (39, 220) es, (rmnc(mn+n1,j),n1=0,nmax1)
                  cycle 
               case (2) 
                  write (39, 220) es, (zmns(mn+n1,j),n1=0,nmax1)
                  cycle 
               case (3) 
                  write (39, 220) es, (lmns(mn+n1,j),n1=0,nmax1)
               end select
            enddo
         enddo
      enddo
 
*
*                 PRINT ASYMMETRIC TERMS -- RmnS and ZmnC
*
      if (iasym .eq. 1) then
         write (39, 200)
         do iloop = 1, 3
            nmax1 = 5
            do mn = 1, mnmax, 6
               if (mn > mnmax - 6) nmax1 = mnmax - mn
               write (39, 210) (ichar2(iloop),ixm(mn+n1),ixn(mn+n1),
     1            n1=0,nmax1)
               do j = 1, ns
                  es = (j - 1)*hs
                  select case (iloop) 
                  case default
                     write (39, 220) es, (rmns(mn+n1,j),n1=0,nmax1)
                     cycle 
                  case (2) 
                     write (39, 220) es, (zmnc(mn+n1,j),n1=0,nmax1)
                     cycle 
                  case (3) 
                     write (39, 220) es, (lmnc(mn+n1,j),n1=0,nmax1)
                  end select
               enddo
            enddo
         enddo
      endif
*
*     DETERMINE RADIAL BETA PROFILE (SURFACE AVERAGED)
*
      phi(1) = 0.
      sqrt_phimod = 0.
      phimod = 0.
      do j = 1, ns
         sqrt_phimod(mn0+mnmax*(j-1)) = sqrt(abs(phi(j)))
         phimod(mn0+mnmax*(j-1)) = phi(j)
      enddo
*     d(BVCO)/ds
      dbzcods(2:ns-1) = (bvco(3:ns)-bvco(2:ns-1))/hs
      dbzcods(ns) = 2*dbzcods(ns-1) - dbzcods(ns-2)
      dbzcods(1) = 2*dbzcods(2) - dbzcods(3)
      beta_vol(1) = 1.5*beta_vol(2) - 0.5*beta_vol(3)
      tb1 = 1.5*beta_vol(ns) - 0.5*beta_vol(ns-1)
      bvco(1) = 1.5*bvco(2) - 0.5*bvco(3)
      buco(1) = 0.
      ub(2:ns) = vp(2:ns)/phip(2:ns)
      mass(2:ns) = mass(2:ns)/abs(phip(2:ns))**gamma
      mass(1) = 1.5*mass(2) - 0.5*mass(3)
      pres(1) = 1.5*pres(2) - 0.5*pres(3)
      vp(1)   = 1.5*vp(2) - 0.5*vp(3)
      tm1 = 1.5*mass(ns) - .5*mass(ns-1)
      tv1 = 1.5*vp(ns) - .5*vp(ns-1)
      tp1 = 1.5*pres(ns) - .5*pres(ns-1)
      ti1 = 1.5*iotas(ns) - .5*iotas(ns-1)
      ub1 = 1.5*ub(ns) - .5*ub(ns-1)
      ub(1) = ub(2)*twopi
CDIR$   IVDEP
      do j = 2, ns - 1
         jp = j + 1
         buco(j) = .5*(buco(jp)+buco(j))
         ub(j) = .5*(ub(jp)+ub(j))*twopi
         vp(j) = .5*(vp(jp)+vp(j))
         mass(j) = .5*(mass(jp)+mass(j))
         pres(j) = 0.5*(pres(jp) + pres(j))
         beta_vol(j) = .5*(beta_vol(jp)+beta_vol(j))
         bvco(j) = .5*(bvco(jp)+bvco(j))
      enddo
      ub(ns) = twopi*ub1
      vp(ns) = tv1
      mass(ns) = tm1
      pres(ns) = tp1
      beta_vol(ns) = tb1
      buco(ns) = 2.*buco(ns) - buco(ns-1)
      bvco(ns) = 2.*bvco(ns) - bvco(ns-1)
      write (39, 230)
*     Normalized FF-Prime : ffp(j)
      do j = 1, ns
         es = (j - 1)*hs
         itor = -twopi*buco(j)
         itors(j) = itor/dmu0
         ipol = twopi*(bvco(j)-bvco(ns))
         ipols(j) = ipol/dmu0
         ffp(j) = bvco(j)*dbzcods(j)/phip(j)/iotaf(j)/
     1      bvco(ns)**2
         write (39, 220) es, bvco(j), buco(j), itor, ipol, ffp(j)
      enddo
      write (39, 240)
      do i = 1, ns
         es = (i - 1)*hs
         write (39, 220) es, twopi**2*vp(i), ub(i), mass(i), pres(i), 
     1      iotaf(i), beta_vol(i)
      enddo
 
      return 
 
*  Do the above before the call to PLOTTER; do the below after
 
      entry wrfcn2 (input_id)
 
      write (39, 260)
      do i = 1, ns
         es = (i - 1)*hs
         write (39, 220) es, jcurv(i), jdotb(i)
      enddo
 
      write (39, 250) beta_vol(1)
 
  200 format(//,40x,'FOURIER COEFFICIENTS X(m,n)',/)
  210 format(//,9x,' S ',4x,6(4x,a5,i1,',',i3,')'),/)
  220 format(1p8e15.3)
  230 format(//,25x,'COVARIANT COMPONENTS OF B',
     1   ' AND INTEGRATED CURRENTS',2/,9x,' S ',10x,'<BZETA>',8x,
     2   '<BTHETA>',8x,'ITOR',11x,'IPOL',8x,'FF''',/)
  240 format(//,9x,' S ',11x,'VP',12x,'dV/dPHI',10x,'MASS',13x,'P',10x,
     1   'IOTA',12x,'BETA',/)
  250 format(//,'  BETA ON AXIS (SUM OVER MODES) = ',1pe10.3)
  260 format(//,9x,' S ',11x,'<JTOR>',2x,'<JdotB>/<bdotgradv>',/,9x,3x
     1   ,11x,'[A/M2]',10x,'[A/M]',/)
 
      end subroutine wrfcn


      subroutine xywtoxyan(xw,yw,xso,yso)
       use pgplotinc, only: pgxscl,pgxorg,pgxsz,pgyscl,pgyorg,pgysz
c       include '/usr/local/pgplot/pgplot.inc' ! might not be there
       integer idx idc
       real xw,yw,xso,yso,xa,ya
         xa(xw,idc)=(xw*pgxscl(idc)+pgxorg(idc))
         ya(yw,idc)=(yw*pgyscl(idc)+pgyorg(idc))
         call pgqid(idx)
         xso=xa(xw,idx)/pgxsz(idx)
         yso=ya(yw,idx)/pgysz(idx)
         return
      end subroutine xywtoxyan

      real function zeroin(ax,bx,f,tol)
      real ax,bx,f,tol
c
c      a zero of the function  f(x)  is computed in the interval ax,bx .
c
c  input..
c
c  ax     left endpoint of initial interval
c  bx     right endpoint of initial interval
c  f      function subprogram which evaluates f(x) for any x in
c         the interval  ax,bx
c  tol    desired length of the interval of uncertainty of the
c         final result ( .ge. 0.0)
c
c
c  output..
c
c  zeroin abcissa approximating a zero of  f  in the interval ax,bx
c
c
c      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
c  without  a  check.  zeroin  returns a zero  x  in the given interval
c  ax,bx  to within a tolerance  4*macheps*abs(x) + tol, where macheps
c  is the relative machine precision.
c      this function subprogram is a slightly  modified  translation  of
c  the algol 60 procedure  zero  given in  richard brent, algorithms for
c  minimization without derivatives, prentice - hall, inc. (1973).
c
c
      real  a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s
c
c  compute eps, the relative machine precision
c
      eps = 1.0
   10 eps = eps/2.0
      tol1 = 1.0 + eps
      if (tol1 .gt. 1.0) go to 10
c
c initialization
c
      a = ax
      b = bx
      fa = f(a)
      fb = f(b)
c
c begin step
c
   20 c = a
      fc = fa
      d = b - a
      e = d
   30 if (abs(fc) .ge. abs(fb)) go to 40
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
c
c convergence test
c
   40 tol1 = 2.0*eps*abs(b) + 0.5*tol
      xm = .5*(c - b)
      if (abs(xm) .le. tol1) go to 90
      if (fb .eq. 0.0) go to 90
c
c is bisection necessary
c
      if (abs(e) .lt. tol1) go to 70
      if (abs(fa) .le. abs(fb)) go to 70
c
c is quadratic interpolation possible
c
      if (a .ne. c) go to 50
c
c linear interpolation
c
      s = fb/fa
      p = 2.0*xm*s
      q = 1.0 - s
      go to 60
c
c inverse quadratic interpolation
c
   50 q = fa/fc
      r = fb/fc
      s = fb/fa
      p = s*(2.0*xm*q*(q - r) - (b - a)*(r - 1.0))
      q = (q - 1.0)*(r - 1.0)*(s - 1.0)
c
c adjust signs
c
   60 if (p .gt. 0.0) q = -q
      p = abs(p)
c
c is interpolation acceptable
c
      if ((2.0*p) .ge. (3.0*xm*q - abs(tol1*q))) go to 70
      if (p .ge. abs(0.5*e*q)) go to 70
      e = d
      d = p/q
      go to 80
c
c bisection
c
   70 d = xm
      e = d
c
c complete step
c
   80 a = b
      fa = fb
      if (abs(d) .gt. tol1) b = b + d
      if (abs(d) .le. tol1) b = b + sign(tol1, xm)
      fb = f(b)
      if ((fb*(fc/abs(fc))) .gt. 0.0) go to 20
      go to 30
c
c done
c
   90 zeroin = b
      return
      end function zeroin

      subroutine surfsep(mnmax_pl, r_pl, z_pl, xm_pl, xn_pl, mnmax_vv,
     1    nu_vv, nv_vv, r_vv, z_vv, xm_vv, xn_vv, nper, igrid, distance, 
     2    drms, dmax, ierr)
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
      dist_max=0
      twopi = 8*atan(1._dp)
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
!      print *,'nper =',nper
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
!         print *, phi_lim(it), v2
!
!   get (r,z) points on plasma surface
! 
         do ip = 1, np
            u2 = (ip - 1)*du
            call dmnf2 (u2, v2, r2(ip), z2(ip))
!            print *,ip, r2(ip), z2(ip)
         enddo


!            print *,'lim: il,mindist,minseg,maxdist,maxseg'

         do il = 2, nlim
            do ip = 1, np
               call lim_dist( r2(ip), z2(ip), r_lim(il-1, it), 
     1                        z_lim(il-1, it), dist(ip))
            enddo
            
            iresult = minloc(dist)
!            print *, phi_lim(it),v2,il,minval(dist),iresult(1)

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

!      print *,r_lim(1),r_lim(2),z_lim(1), z_lim(2)
!      print *,rp, zp, dr, dz, d
      end subroutine lim_dist


EOF


EOC
