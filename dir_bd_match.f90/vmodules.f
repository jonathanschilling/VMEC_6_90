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

