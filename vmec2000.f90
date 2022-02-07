#!/bin/sh
# ************************************************************************************
#        
#        WELCOME TO VMEC2000!
#
# ************************************************************************************
#08/04:  Add LDIAGNO to write out edge current potential for use by Diagno
#03/02:  Add LASYM to indata namelist to control running VMEC in symmetric (default) or
#        asymmetric mode (stellarator symmetric...)
#03/01:  Introduced some damping on the free-boundary motion into SCALFOR in order to
#        accelerate (and improve) the convergence of the free-bdy code.
#03/01:  Modified GUESS_AXIS subroutine to make better guess for axis. A grid search
#        over the plasma cross section is made to find the axis value which maximizes
#        the minimum value of the jacobian for each value of toroidal angle. This algorithm
#        was implemented to improve convergence of VMEC2000 inside STELLOPT optimization code.
#
#03/00:  Begin modification of MHD forces to improve (post-processed) force balance at
#        low aspect ratio. The first step is to divide out the factor of R in FR, FZ 
#        forces. At low aspect ratio this produces mode-coupling which impacts the force
#        balance in the grad-s direction.
#01/00:  VMEC2000 is born!
#        The latest modifications were made to improve the convergence with increasing
#        mode numbers and finer radial meshes.
#        Converted code to internal full mesh differencing for lambda, with a hybrid full/half
#        differencing (blend in radius) for the lambda force differencing. The EXTRAP
#        subroutine was eliminated, as well as the radial damping (in LAMCAL) on the
#        lambda force. WROUT was modified so it continues to write out lambda on the
#        half grid, consistent with previous versions of VMEC.
#
#        Added (but temporarily disabled for backwards compatibility) a radial mesh zoning
#        option to put the toroidal flux on a separate grid
#        from the "s" grid. User can input aphi expansion coefficients for phi = phi(s)
#        in the namelist: phi = aphi(1)*s + aphi(2)*s**2 ... The current, mass and iota
#        expansions are now input as functions of the normalized toroidal flux (same as "s"
#        in the old code). 
#
#        Added an optional command line argument (look for "command line" below for further
#        clarification) to specify the restart file name. The default, if lreset = F on the
#        command line, is to look for the wout file with the input file name extension.
#
#08/99:  Fixed scaling of tcon with radial mesh (1000*hs**2),  to keep it invariant
#04/99:  Version 6.20: WROUT writes out signgs (VMEC now accepts pos jacobian)
#        Additional diagnostic info is also output
#01/99:  Version 6.10: WROUT writes out pres, mass in pascals (1/mu0), 
#        and currents in A(/m**2), NOT internal VMEC units.
#09/98:  Version 6.00. Start splitting out optimization routines.
#        Decoupled BOOZER COORDINATE transformation from vmec/optimizer. Replaced with
#        system call to xbooz
#        Removed lspectrum_dump option. User can run xbooz code separately now.
#        Removed loptim option. User now runs xstelopt code for optimization.
#07/98:  Added LOLDOUT logical to print out fort.8, fort.18 old-style output
#        Added LEDGE_DUMP logical to print out fort.99 (edge-bfield) dump
#        Added (1-cos(u))/2 weight to ripple in optimizer
#05/98:  Version 5.20 introduced as new marker. Added version ID in WOUT file to monitor changes.
#        Fixed various subroutines to avoid 'segmentation fault' on DEC alphas,
#        due to low stack space on those machines (eliminated large automatic
#        vectors, replaced with dynamically allocated arrays). Users of DEC machines should
#        check for adequate stack space using uname -a (or -s). If stack < 30000K, ask system
#        administrator to increase it (although now VMEC will run with stack >= 2048K)
#        Added Mercier criterion to optimization
#        Added approximate external kink mode criterion to load_target
#        Added Mercier condition calculation to jxbforce
#        Fixed damping in lamcal: power = 0.5*m
#        Added logical variable LMAC. If lmac = F (default), the
#        mac-file (unit=nmac0) will be deleted (contains reconstruction data).
#04/98:  Added |B|mn-boozer spectral targets for to optimization chi-sq. (Version 5.10)
#        Added LSPECTRUM_DUMP logical to indata file. If = T, will call load_target to dump
#        boozer spectra, even for free/fixed boundary NOT being optimized.
#        Writes to file bmn_spectrum.file-extension
#        Fixed some 'bugs' associated with converting to a unique boundary representation
#        used by the optimizer. Now the new input file written by the optimizer
#        is consistent with the boundary representation
#03/98:  Improved boundary representation (unique_boundary, convert_boundary) for optimization
#        code, as per Hirshman/Breslau representation.
#01/98:  Fixed default MGRID-FILE so user can specify ANY name, not necessarily
#        prefixed by mgrid.####. (mgrid_file = #### now acceptable in input file.)
#11/97:  Version 5.00: Fixed a number of F90-related bugs associated with
#        the time-stepper and optimization loops. In particular, in EVOLVE, 
#        when irst = 2, we now ALWAYS return before evolving XC, since GC 
#        will NOT be the true force array coming out of FUNCT3D unless irst = 1.
#        Also, the namelists have now been modularized (vmec_input, vmec_seq, optim_params modules).
#        Added spres_ped: (0 < spres_ped <= 1) is a normed flux value. For
#        s > spres_ped, the pressure profile is assumed flat. Default = 1 (no pedestal)
#08/97:  A separate sequence file, SEQ.EXT, is no longer supported. Now, to
#        run a sequence of input files, add a namelist section &VSEQ to an
#        input file, and the filenames (or extensions to input.ext..) will
#        be sequentially executed.
#        Added ntheta, nzeta to namelist (indata).
#08/97:  Added LPROF_OPT to optimizer namelist. When LPROF_OPT=T, then the AI (or AC)
#        coefficients are allowed to vary. Thus, for ncurr=1, the current profile
#        is varied, while for ncurr=0, the iota profile is varied (even though
#        iota is matched in a chi-sq sense). When LPROF_OPT=F, AI (AC) coefficient array
#        is fixed. LPROF_OPT replaces the older LCURPROF_OPT variable (which
#        can still be read in, but is obsolete).
#07/97:  Beginning to phase out gcc/cc (c-precompilation). To this end, have introduced
#        three new logical flags: LFREEB, LRECON, LOPTIM, which all default to
#        F (fixed boundary, no reconstruction, no optimization) unless:
#        (a) mgrid_file is specified in indata file, then LFREEB=T (unless LOPTIM=T is 
#        specified);
#        (b) itse or imse are nonzero in the indata file, then LRECON=T and
#        LFREEB=T (unless LOPTIM=T or no mgrid_file exists)
#        (c) if LOPTIM=T in indata, then LFREEB=LRECON=F regardless of their input
#        values 
#07/97:  Established new makefile structure for VMEC. Now, the user can
#        run vmec.lsqh script on a UNIX machine and only newer vmec files
#        will be updated. To 'make' the vmec executable, type make debug (release)
#        for a debug (release) version.
#04/97:  Began conversion of VMEC to F90 standard (Version 4.00). Making changes to
#        argument passing conventions so as to be compatible with
#        parallel (CRAY) processors. The code will no longer run under
#        F77 with these changes. Converted all includes and common blocks
#        to modules.
#07/96:  Added iresidue=3 condition to turn on FSQR(0,0) for fsq<FOPT_AXIS
#        and turn off RADFOR-pfac time-variation
#07/96:  Moved pressum0 stuff into radfor routine
#07/96:  Fixed GETDIAM routine to use equilibrium pressure balance
#        to compute diamagnetic flux correctly
#05/96:  Adjusted spatial damping parameter in FACLAM to improve convegence
#05/96:  Added constraint weighting-parameter, TCON0, to INDATA input file
#05/96:  WROUT modified: LMN Output on HALF-RADIAL mesh (same as internal VMEC mesh)
#05/96:  Added multigrid capability, NS_ARRAY and FTOL_ARRAY
#05/96:  Removed testing for MSE points changing sign in FIXRECON
#05/96:  Improved algorithm for moving axis to minimize CHISQ in AXISOPT
#04/96:  Upgraded fixed-boundary discrete p & iota profiles in profil1d
#04/96:  Remove Input_Update routine
#04/96:  Removed match to slope of MSE data in GETMSE routine
#04/96:  Added imatch_ip as imatch_edge=3 (for fixed boundary mode)
#03/96:  Made np (number field periods) a true parameter and added the
#        variable nfper to vacuum include file.
#02/96:  Allow Pressure vs. s (s=normalized toroidal flux) input
#        instead of Pressure vs. R (experimental input), when flag
#        LPOFR = .FALSE. in namelist file. This is useful when TRANSP
#        output, for example, is being used and the boundary shape is 
#        uncertain.
#01/96:  ADDED OPTIONS FOR WINDOWS-NT COMPATABILITY AND DEBUGGING
#        SPLIT FILES INTO SEPARATE .F MODULES
#11/95:  ASYMMETRIC THOUGHTS:
#        We need to match Bpol = sqrt(BZ**2 + BR**2) in the
#        MSE measurement, rather than just BZ!!! (Check with
#        S.Batha,F.Levinton) 
#10/95:  Added Variable tension option for weighting tension
#        on iota spline knots independently, i.e., namelist
#        variables TENSI2, FPOLYI
#10/95:  Added NAMELIST variables ISNODES, IPNODES allowing
#        user to pick number of iota,pressure knots
#10/95:  Accelerate the feedback in the imatch_phidege=2 loop
#09/95:  Merged INITSPLINE call into SETSPLINE routine
#06/95:  Converted input pressure coefficients (AM) to 
#        be read in MKS units, NWT/M**2
#05/95:  Changed over to arbitrarily oriented B loops in GETBFLD
#        Added namelist VSEQ for sequential running
#        Added GETLIM, CAUCHY Subroutines to match to limiter position
#02/95:  Added Pres_Offset in Namelist for Thomson Data
#12/94:  Re-wrote GETTHOM, GETMSE subroutines to compute
#        CHI-SQ directly from data (not SMOOTHED and SPLINED data).
#        Added indexing arrays INDEX_I, INDEX_P in NEWPROFIL so
#        calling order of GET... routines is now irrelevant.
#        Looped around ALL GET... routines every IPEDSVD times.
#        Computed extrapolated current in GETMSE correctly.
#        Eliminated user options for fixed (in R) spline nodes
#        in favor of fixed, uniformly distributed number of nodes
#        in SQRT-S (or S) space.
#        Rewrote CHISQ subroutine so all CHI-SQ's are computed
#        from matrix elements, rather than from defining relations.
#        Added VERSION NO. to THREED1 file output.
#10/94:  Improved pfac, phifac algorithms (a little)
#09/94:  Pass J-dot-B material through WOUT
#     :  User should customize mgrid_defarea in subr readin for his site
#04/94:  LAMBDA DIFFERENCING ON HALF-MESH IMPLEMENTED
#03/94:  ELIMINATED RMNSS FOR PUSHING AXIS, REPLACED WITH RMNCC(JS=1)
#        THIS WAS NECESSITATED BY DISCREPANCIES FOR SMALL ITSE...
#01/94:  IMPLEMENTED CHANGES FOR FIXED, FREE BOUNDARY
#        FOR UP-DOWN NON-SYMMETRIC PLASMAS (symmetry_mode qualifier)
#11/93:  REPLACED WEIGHTS WITH STANDARD DEVIATIONS
#10/93:  IMPLEMENTED PHIEDGE MATCHING BASED ON WIDTH OF PRESSURE PROFILE DATA
#VERSION WHICH USES SQRT(S) MESH FOR SPLINES AND USES BOTH INNER AND
#OUTER EDGES FOR CURRENT MATCH (OPTIONAL RADIAL REDISTRIBUTION OF DATA)
#---------------------------------------------------------------------
cat > temp.c << "EOC"
cat > vmodules.f << "EOF"
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
EOF

cat > meq.f << "EOF"
      program vmec
      use vmec_input
      use vmec_seq
      use read_namelist_mod
      use safe_open_mod
      use vparams, only: nlog, nlog0
      implicit none
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: nseq0 = 12
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: numargs, ierr_vmec, index_end,
     1   iopen, isnml, iread, iseq, index_seq, 
     2   index_dat, ireset, iunit
      character*120 :: input_file, seq_ext, reset_file_name, arg
      character*120 :: log_file
      character*120, dimension(10) :: command_arg
      logical :: lfirst=.true., lreseta, lscreen
C-----------------------------------------------
!***
!                              D   I   S   C   L   A   I   M   E   R
!
!       You are using a BETA version of the program VMEC, which is currently
!       under development by S. P. Hirshman at the Fusion Energy Division,
!       Oak Ridge National Laboratory.  Please report any problems or comments
!       to him.  As a BETA version, this program is subject to change
!       and improvement without notice.
!
!       1. CODE SYNOPSIS
!
!       THIS PROGRAM - VMEC (Variational Moments Equilibrium Code)  -
!       SOLVES THREE-DIMENSIONAL MHD EQUILIBRIUM EQUATIONS USING
!       FOURIER SPECTRAL (MOMENTS) METHODS. A CYLINDRICAL COORDINATE
!       REPRESENTATION IS USED (R-Z COORDINATES). THE POLOIDAL
!       ANGLE VARIABLE IS RENORMALIZED THROUGH THE STREAM FUNCTION
!       LAMBDA, WHICH IS SELF-CONSISTENTLY DETERMINED AND DIFFERENCED
!       VARIATIONALLY ON THE HALF-RADIAL MESH. THE POLOIDAL ANGLE IS
!       DETERMINED BY MINIMIZING <M> = m**2 S(m) , WHERE S(m) =
!       Rm**2 + Zm**2 . AN EVEN-ODD DECOMPOSITION IN THE POLOIDAL MODE
!       NO. OF R,Z, AND LAMDA IS USED TO IMPROVE RADIAL RESOLUTION.
!       A FREE-BOUNDARY OPTION IS AVAILABLE (FOR lfreeb=T), WITH A
!       USER-SUPPLIED DATA-FILE "MGRID" NEEDED TO COMPUTE THE PLASMA
!       VACUUM FIELD COMPONENTS BR, BPHI, BZ (see subroutine BECOIL)
!
!       THE MAGNETIC FIELD IS REPRESENTED INTERNALLY AS FOLLOWS:
!
!       B(s,u,v) = grad(phiT) X ( grad(u) + grad(lambda) ) +
!
!                  iota(s) * grad(v) X grad(phiT)
!
!       where phiT is the toroidal flux (called phi in code) and
!       u,v are the poloidal, toroidal angles, respectively.
!
!       2. ADDITIONAL CODES REQUIRED
!       For the fixed boundary calculation, the user must provide the Fourier
!       coefficients for the plasma boundary (the last surface outside of which
!       the pressure gradient vanishes). For all but the simplest geometry, the
!       SCRUNCH code (available from R. Wieland), based on the DESCUR curve-fitting
!       code, can be used to produce the optimized VMEC Fourier representation for
!       an arbritrary closed boundary (it need not be a 'star-like' domain, nor
!       need it possess vertical, or 'stellarator', symmetry).
!   
!       For the free boundary calculation, the MAKEGRID code (available upon
!       request) is needed to create a binary Green''s function table for the
!       vacuum magnetic field(s) and, if data analysis is to be done, flux and
!       field loops as well. The user provides a subroutine (BFIELD) which can be
!       called at an arbitrary spatial location and which should return the three
!       cylindrical components of the vacuum field at that point. (Similary,
!       locations of diagnostic flux loops, Rogowski coils, etc. are required if
!       equilibrium reconstruction is to be done.)
!   
!       Plotting is handled by a stand-alone package, PROUT.NCARG (written by
!       R. M. Wieland). It uses NCAR-graphics calls and reads the primary VMEC output
!       file, WOUT.EXT, where 'EXT' is the command-line extension of the INPUT file.
!
!   
!       3. UNIX SCRIPT SETUP PARAMETERS
!       The VMEC source code (vmec.lsqh) is actually a UNIX script file which uses
!       the C-precompiler to produce both the machine-specific Fortran source and a
!       make-file specific to any one of the following platforms:
!
!       IBM-RISC6000, CRAY, ALPHA (DEC-STATION), HP-UX WORKSTATION,
!       WINDOWS-NT, DEC-VMS
!
!       Additional platforms are easy to add to the existing script as required.
!   
!
!       4. FORTRAN PARAMETER STATEMENTS set by user
!       In the Fortran-90 version of VMEC these parameter statements have
!       been replaced by dynamic memory allocation. So the user should set the
!       run-time parameters ns (through ns_array), mpol, ntor in the namelist INDATA.
!
!
!       Added features since last edition
!       1. Implemented preconditioning algorithm for R,Z
!       2. The physical (unpreconditioned) residuals are used
!          to determine the level of convergence
!       3. The original (MOMCON) scaling of lambda is used, i.e.,
!          Bsupu = phip*(iota - lamda[sub]v)/sqrt(g). This is needed to
!          maintain consistency with the time-stepper for arbitrary PHIP.
!
!       WRITTEN BY S. P. HIRSHMAN (8/28/85 - REVISED 3/1/86) BASED ON
!       1. S. P. Hirshman and J. C. Whitson, Phys. Fluids 26, 3553 (1983).
!       2. S. P. Hirshman and H. K. Meier, Phys. Fluids 28, 1387 (1985).
!       3. S. P. Hirshman and D. K. Lee, Comp. Phys. Comm. 39, 161 (1986).
!***

!    
!     Read in command-line arguments to get input file or sequence file,
!     screen display information, and restart information
!
      call getcarg(1, command_arg(1), numargs)
      do iseq = 2, numargs
         call getcarg(iseq, command_arg(iseq), numargs)
      end do   

      lreseta = .true.            !!Default value: runvmec MUST be called this way the first time
      lscreen = .true.
      
      if (numargs .lt. 1) then
         stop 'Invalid command line'
      else if (command_arg(1).eq.'-h' .or. command_arg(1).eq.'/h') then
         print *,
     1   ' ENTER INPUT FILE NAME OR INPUT-FILE SUFFIX ON COMMAND LINE'
         print *
         print *,' For example: '
         print *,'    xvmec input.tftr OR xvmec tftr ',
     1           'OR xvmec ../input.tftr'
         print *
         print *,' Sequence files, containing a LIST of input files',
     1           ' are also allowed: '
         print *,'    xvmec input.tftr_runs'
         print *
         print *,' Here, input.tftr_runs contains a &VSEQ namelist',
     1           ' entry'
         print *
         print *,' Additional (optional) command arguments are',
     1           ' allowed:'
         print *
         print *,'    xvmec <filename> noscreen F'
         print *
         print *,' noscreen: supresses all output to screen ',
     1           ' (default, or "screen", displays output)'    
         print *,' F (or T): if "T", forces reset on',
     1           ' a coarse mesh (used for sequencing control)'   
         
         stop
      else if (numargs .gt. 1) then
         arg = command_arg(2)
         if (trim(arg).eq.'noscreen' .or. trim(arg).eq.'NOSCREEN')
     1      lscreen = .false.
      end if
      if (numargs .gt. 2) then
          arg = command_arg(3)
          if (arg(1:1).eq.'f' .or. arg(1:1).eq.'F') lreseta = .false.   
      end if    
      if (numargs .gt. 3) then
          reset_file_name = command_arg(4)
      end if    
 

!
!     Determine type of file opened (sequential or input-data)      
!     ARG1 (char var)
!          By default, ARG1 obtained from the command
!          line is parsed as follows to determine the input data file(s):
!               a. Attempt to open file ARG1 (full path + file name).
!                  Look for the VSEQ namelist to obtain nseq, nseq_select, and
!                  extension array. If they exist and nseq>0, VMEC will run
!                  sequentially using input determined from the array EXTENSION[i] 
!                  or input.EXTENSION[i]
!               b. If the command argument is not a sequence namelist, then the data file 
!                  ARG1 or input.ARG1 is read directly, with NSEQ=1.
!
      arg = command_arg(1)
      index_dat = index(arg,'.')
      index_end = len_trim(arg)
      if (index_dat .gt. 0) then
         seq_ext  = arg(index_dat:index_end)
         input_file = trim(arg)
      else
         seq_ext = trim(arg)
         input_file = 'input.'//trim(seq_ext)
      end if   
 
      if (numargs .le. 3) reset_file_name = 'wout.' // seq_ext
      
      nseq = 1
      nseq_select(1) = 1
      extension(1) = input_file
!
!     READ IN NAMELIST VSEQ TO GET ARRAY
!     OF INPUT FILE EXTENSIONS AND INDEXING ARRAY, NSEQ_SELECT
!
      nlog = nlog0
      iunit = nseq0
      do iseq = 1, 2
         if (iseq .eq. 1) then
           arg = input_file
         else
           arg = seq_ext
         end if
         call safe_open(iunit, iopen, trim(arg), 'old', 'formatted')
         if (iopen .eq. 0) then
           call read_namelist (iunit, isnml, 'vseq')
           if (isnml.eq.0 .and. nseq .gt. nseqmax) stop 'NSEQ>NSEQMAX'
 
!
!       OPEN FILE FOR STORING SEQUENTIAL RUN HISTORY
!
           if (isnml .eq. 0) then
              log_file = 'log.'//seq_ext
 
              call safe_open(nlog, iread, log_file, 'replace', 
     1           'formatted')
              if (iread .ne. 0) then
                 print *, log_file, 
     1           ' LOG FILE IS INACCESSIBLE: IOSTAT= ',iread
                 stop 3
              else
                 exit        !!Break out of loop
              end if
           endif
        endif  

        close (iunit)

      end do  

!
!     CALL EQUILIBRIUM SOLVER 
!
!     nseq_select:      if sequence file (VSEQ namelist given with nseq >0)
!                       array giving indices into EXTENSION array prescribing
!                       the order in which the input files are run by VMEC
!     nseq:             number of sequential VMEC runs to make
!
!
!     CALL VMEC WITH POSSIBLE SEQUENCE EXTENSION (SEQ_EXT)
!     AND ARRAY OF INPUT FILE EXTENSIONS (EXTENSION)
!
      do iseq = 1, nseq 
         index_seq = nseq_select(iseq)
         ireset = 0
         ierr_vmec = 0
         if (iseq .gt. 1) reset_file_name = 
     1       'wout.' // trim(extension(index_seq))
 100     continue
         call runvmec (extension(index_seq), iseq-1, lreseta, ierr_vmec, 
     1                 lfirst, lscreen, reset_file_name)
         lfirst = .false. 
         select case (ierr_vmec) 
!        case (1:2)    !BAD JACOBIAN AFTER 75 ITERATIONS...
!          ireset = ireset + 1
!          lreseta = .true.
!          if (ireset .le. 2) go to 100
         case (4)                                !Try a few more iterations
           ireset = ireset + 1
           lreseta = .false.
           if (ireset .le. 1) then
              if (lscreen) write (*, '(/,1x,a)') 
     1           'RUNNING A FEW MORE ITERATIONS THAN REQUESTED'
              go to 100
           else if (lscreen) then
              print *, 'DECREASE DELT OR INCREASE NITER'
           endif
         case (6)    !BAD JACOBIAN AFTER AXIS RESET: TRY DECREASING TO NS=3
           ireset = ireset + 1
           lreseta = .true.
           if (ireset .le. 1) go to 100             
         case default
           lreseta = .false.
         end select
      end do

!
!     FREE ANY LONG-TERM (PERSISTENT THROUGH ISEQ > 1, OR XC, SCALXC FOR
!     ITERATIVE OPTIMIZATION) POINTERS
!
      call free_persistent_mem
 
      close (nlog)
 
      end program vmec
      

      subroutine bcovar (lu, lv)
      use vmec_main
      use vmec_params, only: ns4, signgs
      use realspace, weight_f => sqrts, sqrts => sqrts
      use vforces, r12 => armn_o, ru12 => azmn_e, gsqrt => azmn_o,
     1   rs => bzmn_e, zs => brmn_e, zu12 => armn_e,
     2   bsubu_e => clmn_e, bsubv_e => blmn_e, bsubu_o => clmn_o,
     3   bsubv_o => blmn_o, bsq => bzmn_o, phipog => brmn_o
      use vsvd, only: phifac, phifsave, imovephi
      use xstuff, only: xc
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(nrzt,0:1) :: lu, lv
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
!     GENERALLY, IF TEMPORAL CONVERGENCE IS POOR, TRY TO INCREASE PDAMP (< 1)
      real(rprec), parameter :: c1p5 = 1.5_dp, p5 = 0.5_dp,
     1    pdamp = 0.05_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: l, is, js, ndim
      real(rprec) :: r2, r3
      real(rprec), dimension(:), pointer :: r12sq, bsupu, bsupv,
     1    bsubuh, bsubvh
      real(rprec) :: luu12, luv12, lvv12
      real(rprec) :: arnorm, aznorm, volume, tcon_mul
      real(rprec), allocatable, dimension(:) :: luu, luv, lvv
      real(rprec), external :: dot_g
C-----------------------------------------------
      ndim = 1+nrzt
      r12sq => bsubu_o

      allocate (luu(ndim), luv(ndim), lvv(ndim), stat = is)

      if (is .ne. 0) stop 'allocation error in bcovar'

!
!     FOR OPTIMIZATION ON CRAY, MUST USE COMPILER DIRECTIVES TO
!     GET VECTORIZATION OF LOOPS INVOLVING (MORE THAN ONE) POINTER!
!
      guu = 0                           ! mcz
      guv = 0;  gvv = 0

!
!     COMPUTE METRIC ELEMENTS GIJ ON HALF MESH
!     FIRST, GIJ = EVEN PART (ON FULL MESH), LIJ = ODD PART (ON FULL MESH)
!     THEN, GIJ(HALF) = < GIJ(even)> + SHALF < GIJ(odd) >
!

CDIR$ IVDEP
      do l = 1, nrzt
         r2 = sqrts(l)*sqrts(l)
         guu(l)   =  ru(l,0)*ru(l,0) + r2*ru(l,1)*ru(l,1)
     1            +  zu(l,0)*zu(l,0) + r2*zu(l,1)*zu(l,1)
         luu(l)   = (ru(l,0)*ru(l,1) + zu(l,0)*zu(l,1))*2
         r12sq(l) = r1(l,0)*r1(l,0) + r2*r1(l,1)*r1(l,1)
         phipog(l)= 2*r1(l,0)*r1(l,1)
      end do

      if (lthreed) then
CDIR$ IVDEP
         do l = 1, nrzt
            r2 = sqrts(l)*sqrts(l)
            guv(l)   =  ru(l,0)*rv(l,0) + r2*ru(l,1)*rv(l,1)
     1               +  zu(l,0)*zv(l,0) + r2*zu(l,1)*zv(l,1)
            luv(l)   =  ru(l,0)*rv(l,1) + ru(l,1)*rv(l,0)
     1               +  zu(l,0)*zv(l,1) + zu(l,1)*zv(l,0)
            gvv(l)   =  rv(l,0)*rv(l,0) + r2*rv(l,1)*rv(l,1)
     1               +  zv(l,0)*zv(l,0) + r2*zv(l,1)*zv(l,1)
            lvv(l)   = (rv(l,0)*rv(l,1) + zv(l,0)*zv(l,1))*2
         end do
      end if
CDIR$ IVDEP
      do l = nrzt, 2, -1
         guu(l) = p5*(guu(l) + guu(l-1) + shalf(l)*(luu(l) + luu(l-1)))
         r12sq(l) = p5*(r12sq(l) + r12sq(l-1) + shalf(l)*
     1                (phipog(l) + phipog(l-1)))
      end do
      if (lthreed) then
CDIR$ IVDEP
         do l = nrzt, 2, -1
            guv(l) = p5*(guv(l) + guv(l-1) +
     1         shalf(l)*(luv(l) + luv(l-1)))
            gvv(l) = p5*(gvv(l) + gvv(l-1) +
     1         shalf(l)*(lvv(l) + lvv(l-1)))
         end do
      end if

      do js = 2, ns
         vp(js) = signgs*DOT_G(nznt,gsqrt(js),ns,wint(js),ns)
      end do
      if (iter2 .eq. 1) voli = twopi*twopi*hs*sum(vp(2:ns))

      gvv(2:nrzt) = gvv(2:nrzt) + r12sq(2:nrzt)
      phipog(2:nrzt) = phip(2:nrzt)/gsqrt(2:nrzt)
      phipog(1)    = 0
      phipog(ndim) = 0
!
!     RECONSTRUCT IOTA, PRESSURE PROFILE FROM DATA
!
      if (lrecon) call newprofil (phipog)
      if (.not.lrecon .and. imovephi.gt.0) call newphi (phipog)

!
!     STORE CURRENT(0), MAGNETIC PITCH, PRESSURE DATA
!     STORE HALF-MESH VALUES OF LU, LV FOR ACCURATE CURRENT CALCULATION
!
      if (iequi .eq. 1) call storesvd (r1(1,0), r1(1,1), lu(1,0),
     1   lu(1,1), lv(1,0), phipog, zu0)

!
!     COMPUTE COVARIANT COMPONENTS OF B ON RADIAL FULL-MESH
!     NOTE: LU = 1+LAMU, LV = -LAMV COMING INTO THIS ROUTINE
!     WILL ADD IOTAF, PHIP/GSQRT FACTOR LATER...
!          LV == BSUPU = (PHIP/GSQRT)*(iota - LAMV),
!          LU == BSUPV = (PHIP/GSQRT)*(1    + LAMU)
!

      bsupu => bsubu_o                  !!Save these later in lv(l,0), lu(l,0)
      bsupv => bsubv_o

!
!     FIRST, PUT LAMBDA DERIVATIVES ON RADIAL HALF-MESH
!
CDIR$ IVDEP
      do l = 2, nrzt
         bsupv(l) = p5*phipog(l)*(lu(l,0) + lu(l-1,0) + shalf(l)*
     1                           (lu(l,1) + lu(l-1,1)))
         bsupu(l) = p5*phipog(l)*(lv(l,0) + lv(l-1,0) + shalf(l)*
     1                           (lv(l,1) + lv(l-1,1)))
      end do

      bsupv(1) = 0
      bsupu(1) = 0

!
!     COMPUTE UPDATED IOTA PROFILE
!
      call getiota(phipog, bsupu, bsupv)


!     ADD PRESENT VALUE OF IOTAF HERE.
      do js = 1, ns
         lv(js:nrzt:ns,0) = lv(js:nrzt:ns,0) + iotaf(js)
      end do

!
!     NEXT COMPUTE LAMBDA FORCES ON FULL MESH BY AVERAGING HALF-MESH METRICS
!
      luu = phipog(:ndim)*guu
      luv = phipog(:ndim)*guv
      lvv = phipog(:ndim)*gvv

CDIR$ IVDEP
      do l = 1,nrzt
         luu12 = p5*(luu(l) + luu(l+1))
         luv12 = p5*(luv(l) + luv(l+1))
         lvv12 = p5*(lvv(l) + lvv(l+1))
         bsubu_e(l) = luu12*lv(l,0) + luv12*lu(l,0)
         bsubv_e(l) = luv12*lv(l,0) + lvv12*lu(l,0)
      end do

      luu = luu*shalf
      luv = luv*shalf
      lvv = lvv*shalf

CDIR$ IVDEP
      do l = 1,nrzt
         luu12 = p5*(luu(l) + luu(l+1))
         luv12 = p5*(luv(l) + luv(l+1))
         lvv12 = p5*(lvv(l) + lvv(l+1))
         bsubu_e(l) = bsubu_e(l)+ luu12*lv(l,1) + luv12*lu(l,1)
         bsubv_e(l) = bsubv_e(l)+ luv12*lv(l,1) + lvv12*lu(l,1)
      end do

!
!     FINALLY, COMPUTE LAMBDA FORCES ON RADIAL HALF-MESH
!
      lu(:nrzt,0) = bsupv(:nrzt)
      lv(:nrzt,0) = bsupu(:nrzt)

      bsubuh => bsubu_o
      bsubvh => bsubv_o

      bsubuh(:nrzt) = guu(:nrzt)*lv(:nrzt,0) + guv(:nrzt)*lu(:nrzt,0)
      bsubvh(:nrzt) = guv(:nrzt)*lv(:nrzt,0) + gvv(:nrzt)*lu(:nrzt,0)

      bsubuh(ndim) = 0
      bsubvh(ndim) = 0

      do js = 1,ns
         fpsi(js) = DOT_G(nznt,bsubvh(js),ns,wint(js),ns)
      enddo

!
!     COMPUTE TOROIDAL AND POLOIDAL CURRENTS
!
      rbtor = c1p5*fpsi(ns) - cp5*fpsi(ns-1)
      rbtor0= c1p5*fpsi(2)  - cp5*fpsi(3)
      ctor = signgs*twopi*sum((c1p5*bsubuh(ns:nrzt:ns) -
     1   cp5*bsubuh(ns-1:nrzt:ns))*wint(ns:nrzt:ns))

!
!     COMPUTE KINETIC PRESSURE ON HALF-GRID
!     IF R*dp/ds NONVARIATIONAL FORM TO BE USED IN FORCES, UNCOMMENT ALL THE "RPRES" 
!     COMMENTS IN THIS SUBROUTINE, AS WELL AS FORCES, FUNCT3D SUBROUTINES
!
      pres(2:ns) = mass(2:ns)/vp(2:ns)**gamma
      wp = hs*sum(vp(2:ns)*pres(2:ns))
      do js = 2,ns
         bsq(js:nrzt:ns) = pres(js)
      end do   
!RPRES     bsq = 0
!RPRES     do js = 2,ns
!RPRES        lu(js:nrzt:ns,1) = pres(js)
!RPRES     end do

!
!     COMPUTE MAGNETIC PRESSURE + KINETIC PRESSURE
!
      bsq(:nrzt) = bsq(:nrzt) + p5*(lv(:nrzt,0)*bsubuh(:nrzt) +
     1      lu(:nrzt,0)*bsubvh(:nrzt))
      wb = -wp + hs*abs(sum(wint(:nrzt)*gsqrt(:nrzt)*bsq(:nrzt)))
!RPRES      wb = hs*abs(sum(wint(:nrzt)*gsqrt(:nrzt)*bsq(:nrzt)))


!
!     AVERAGE LAMBDA FORCES ONTO FULL MESH
!
CDIR$ IVDEP
      do l = 1, nrzt
         r2 = (1 - weight_f(l)*weight_f(l))**2 * pdamp
         r3 = p5*(1 - r2)
         bsubu_e(l) = bsubu_e(l) * r2
     1              + r3*(bsubuh(l) + bsubuh(l+1))
         bsubv_e(l) = bsubv_e(l) * r2
     1              + r3*(bsubvh(l) + bsubvh(l+1))
      end do

!
!     COMPUTE COVARIANT BSUBU,V (EVEN, ODD) ON HALF RADIAL MESH
!     FOR FORCE BALANCE AND RETURN (IEQUI=1)
!
      if (iequi .eq. 1) then

!DBG     luu = bsubuh
!DBG     lvv = bsubvh

!RPRES         bsq(:nrzt) = bsq(:nrzt) + lu(:nrzt,1)

         do js = ns-1,2,-1
         do l = js, nrzt, ns
            bsubuh(l) = 2*bsubu_e(l) - bsubuh(l+1)
            bsubvh(l) = 2*bsubv_e(l) - bsubvh(l+1)
         end do
         end do

         bsubu_e(:nrzt) = bsubuh(:nrzt)
         bsubv_e(:nrzt) = bsubvh(:nrzt)

         bsubu_o(:nrzt) = shalf(:nrzt)*bsubu_e(:nrzt)
         bsubv_o(:nrzt) = shalf(:nrzt)*bsubv_e(:nrzt)


!DBG     do js = 2, ns
!DBG        print *, 'JS = ', js
!DBG        print *,
!DBG 1      '    BSUBU(old)     BSUBU(new)    BSUBV(old)     BSUBV(new)'
!DBG        do l = js, nrzt, ns
!DBG           write(*,1223) luu(l), bsubu_e(l), lvv(l), bsubv_e(l)
!DBG        end do
!DBG     end do

 1223    format(1p4e14.5)

         go to 1000

      end if

      bsubu_o(:nrzt) = sqrts(:nrzt)*bsubu_e(:nrzt)
      bsubv_o(:nrzt) = sqrts(:nrzt)*bsubv_e(:nrzt)

!
!     COMPUTE R,Z AND LAMBDA PRE-CONDITIONING MATRIX
!     ELEMENTS AND FORCE NORMS EVERY NS4 STEPS.
!
      if(mod(iter2-iter1,ns4).eq.0)then
         phifsave = phifac
         phipog(:nrzt) = phipog(:nrzt)*wint(:nrzt)
         call lamcal(phipog, guu, guv, gvv)
         call precondn(lu,bsq,gsqrt,r12,zs,zu12,zu,zu(1,1),
     1                 z1(1,1),arm,ard,brm,brd,crd)
         call precondn(lu,bsq,gsqrt,r12,rs,ru12,ru,ru(1,1),
     1                r1(1,1),azm,azd,bzm,bzd,crd)

         guu(:ndim) = guu(:ndim)*r12(:ndim)**2
         volume = hs*sum(vp(2:ns))
         r2 = max(wb,wp)/volume
         fnorm = one/(sum(guu(1:nrzt)*wint(1:nrzt))*(r2*r2))
         fnorm1 = one/sum(xc(1+ns:2*irzloff)**2)

!
!        COMPUTE CONSTRAINT FORCE SCALING FACTOR (TCON)
!
!        OVERRIDE USER INPUT VALUE HERE
!
         r2 = ns
         tcon0 = min(abs(tcon0), one)                      !!ignore large tcon0 value from old-style file
         tcon_mul = tcon0*(1 + r2*(one/60 + r2/(200*120)))

         do js = 2, ns-1
           arnorm = sum(wint(js:nrzt:ns)*ru0(js:nrzt:ns)**2)
           aznorm = sum(wint(js:nrzt:ns)*zu0(js:nrzt:ns)**2)
           if (arnorm .eq. zero .or. aznorm .eq. zero)
     1     stop 'arnorm or aznorm=0'

           tcon(js) = min(abs(ard(js,1)/arnorm),abs(azd(js,1)/
     1          aznorm))                 * tcon_mul * (32*hs)**2
         end do
         tcon(ns) = cp5*tcon(ns-1)

      endif

!
!     STORE LU * LV COMBINATIONS USED IN FORCES
!
      do l=2,nrzt
         guu(l) = lv(l,0)*lv(l,0)*gsqrt(l)
         guv(l) = lv(l,0)*lu(l,0)*gsqrt(l)
         gvv(l) = lu(l,0)*lu(l,0)*gsqrt(l)
         lv(l,0)  = bsq(l)*gsqrt(l)/r12(l)
         lu(l,0)  = bsq(l)*r12(l)
      enddo

 1000 continue

      deallocate (luu, luv, lvv, stat = l)

      end subroutine bcovar


      subroutine getiota(phipog, bsupu, bsupv)
      use vmec_main
      use realspace
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(nrzt), intent(in) :: phipog, bsupv
      real(rprec), dimension(nrzt), intent(inout) :: bsupu
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: p5=0.5_dp, c1p5=1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: js, l
      real(rprec) :: top, bot
C-----------------------------------------------

      if (ncurr .eq. 0) goto 100
      do js = 2, ns
         top = jcurv(js)
         bot = 0
         do l = js, nrzt, ns
            top = top - wint(l)*(guu(l)*bsupu(l) + guv(l)*bsupv(l))
            bot = bot + wint(l)*phipog(l)*guu(l)
         end do
         iotas(js) = top/bot
      end do

!     Do not compute iota too near origin
      iotaf(1)  = c1p5*iotas(2) - p5*iotas(3)           !!zero gradient near axis
      iotaf(ns) = c1p5*iotas(ns) - p5*iotas(ns-1)
      do js = 2, ns-1
         iotaf(js) = p5*(iotas(js) + iotas(js+1))
      end do

 100  continue

      do js = 2, ns
         bsupu(js:nrzt:ns) = bsupu(js:nrzt:ns) + phipog(js:nrzt:ns)
     1                     * iotas(js)
      end do

      end subroutine getiota


      subroutine lamcal(phipog, guu, guv, gvv)
      use vmec_main
      use vmec_params, only: ntmax, jlam
      use realspace, only: sqrts
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,nznt), intent(in) ::
     1   phipog, guu, guv, gvv
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: damping_fac = -2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: m,n,js
      real(rprec) :: tnn, tnm, tmm, power
C-----------------------------------------------

      blam(:ns) = sum(guu *phipog, dim=2)
      clam(:ns) = sum(gvv *phipog, dim=2)
      dlam(:ns) = sum(guv *phipog, dim=2)
      blam(1) =  blam(2)
      clam(1) =  clam(2)
      dlam(1) =  dlam(2)
      blam(ns+1) =  0
      clam(ns+1) =  0
      dlam(ns+1) =  0
      do js = 2, ns
        blam(js) = cp5*(blam(js) + blam(js+1))
        clam(js) = cp5*(clam(js) + clam(js+1))
        dlam(js) = cp5*(dlam(js) + dlam(js+1))
      end do

!
!       REDUCE FACLAM AT SOME FSQ THRESHOLD TO IMPROVE 3D CONVERGENCE
!
      do m = 0, mpol1
         tmm = m*m
         power = min(tmm/256, 8._dp)
         do n = 0, ntor
            if (m.eq.0 .and. n.eq.0) cycle
            tnn = (n*nfp)**2
            tnm = 2*m*n*nfp
            do js = jlam(m), ns
               faclam(js,n,m,1) = 2*damping_fac/
     1         ((blam(js)+blam(js+1))*tnn
     2         + sign((dlam(js)+dlam(js+1)),blam(js))*tnm
     3         + (clam(js) + clam(js+1))*tmm)
     4         * sqrts(js)**power                                   !Damps m > 16 modes
            end do
         end do
      end do

      do n = 2, ntmax
         faclam(:ns,0:ntor,0:mpol1,n) = faclam(:ns,0:ntor,0:mpol1,1)
      end do

      end subroutine lamcal


      subroutine precondn(lu1, bsq, gsqrt, r12, xs, xu12, xue, xuo,
     1   xodd, axm, axd, bxm, bxd, cx)
      use vmec_main
      use realspace
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(nrzt), intent(in) ::
     1  lu1, bsq, gsqrt, r12, xs, xu12, xue, xuo, xodd
      real(rprec), dimension(ns+1,2), intent(out) ::
     1  axm, axd, bxm, bxd
      real(rprec), dimension(ns+1), intent(out) :: cx
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: js, l, lk
      real(rprec), dimension(:,:), allocatable :: ax, bx
      real(rprec), dimension(:), allocatable :: ptau
      real(rprec) t1, t2, t3
C-----------------------------------------------
!
!     COMPUTE PRECONDITIONING MATRIX ELEMENTS FOR R,Z
!     FORCE (ALL ARE MULTIPLIED BY 0.5)
!

      allocate (ax(ns+1,4), bx(ns+1,4), ptau(nznt))

      ax = 0
      bx = 0
      cx = 0

      do 20 js = 2,ns
!
!     COMPUTE DOMINANT (1/DELTA-S)**2 PRECONDITIONING
!     MATRIX ELEMENTS
!
        lk = 0
        do l = js,nrzt,ns
          lk = lk + 1
          ptau(lk) = r12(l)*r12(l)*bsq(l)*wint(l)/gsqrt(l)
          t1 = xu12(l)*ohs
          t2 = cp25*(xue(l)/shalf(js) + xuo(l))/shalf(js)
          t3 = cp25*(xue(l-1)/shalf(js) + xuo(l-1))/shalf(js)
          ax(js,1) = ax(js,1) + ptau(lk)*t1*t1
          ax(js,2) = ax(js,2) + ptau(lk)*(-t1+t3)*(t1+t2)
          ax(js,3) = ax(js,3) + ptau(lk)*(t1+t2)*(t1+t2)
          ax(js,4) = ax(js,4) + ptau(lk)*(-t1+t3)*(-t1+t3)
        end do
!
!       COMPUTE ORDER M**2 PRECONDITIONING MATRIX ELEMENTS
!
        lk = 0
        do l = js,nrzt,ns
          lk = lk+1
          t1 = cp5*(xs(l) + cp5*xodd(l)/shalf(js))
          t2 = cp5*(xs(l) + cp5*xodd(l-1)/shalf(js))
          bx(js,1) = bx(js,1) + ptau(lk)*t1*t2
          bx(js,2) = bx(js,2) + ptau(lk)*t1*t1
          bx(js,3) = bx(js,3) + ptau(lk)*t2*t2
          cx(js) = cx(js) + cp25*lu1(l)**2*gsqrt(l)*wint(l)
        end do
 20   continue

      do js = 1,ns
        axm(js,1) =-ax(js,1)
        axd(js,1) = ax(js,1) + ax(js+1,1)
        axm(js,2) = ax(js,2) * sm(js) * sp(js-1)
        axd(js,2) = ax(js,3)*sm(js)**2 + ax(js+1,4)*sp(js)*sp(js)
        bxm(js,1) = bx(js,1)
        bxm(js,2) = bx(js,1) * sm(js) * sp(js-1)
        bxd(js,1) = bx(js,2) + bx(js+1,3)
        bxd(js,2) = bx(js,2)*sm(js)**2 + bx(js+1,3)*sp(js)*sp(js)
        cx(js)    = cx(js) + cx(js+1)
      end do

      axm(ns+1,:) = 0
      bxm(ns+1,:) = 0
      axd(ns+1,:) = 0
      bxd(ns+1,:) = 0
      cx(ns+1) = 0

      deallocate (ax, bx, ptau)

      end subroutine precondn

      subroutine forces
      use vmec_main
      use realspace, z1 => z1, gcon => z1
      use vforces, crmn_e => crmn_e, lv_e => crmn_e, 
     1   czmn_e => czmn_e, lu_e => czmn_e,
     2   czmn_o => czmn_o, lu_o => czmn_o
!RPRES , presg => czmn_o
      use vsvd, only: torflux_edge => torflux
      implicit none
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: p5 = 0.5_dp, p25 = p5*p5
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: l, lk, ndim
      real(rprec), dimension(:), allocatable :: 
     1    bsqr, gvvs, guvs, guus
      real(rprec) :: rcon1, zcon1, dphids
C-----------------------------------------------

      ndim = 1+nrzt 
      allocate (bsqr(ndim), gvvs(ndim), guvs(ndim), guus(ndim), 
     1           stat=l)
      if (l .ne. 0) stop 'Allocation error in VMEC FORCES'
      
!
!     ON ENTRY, ARMN=ZU,BRMN=ZS,AZMN=RU,BZMN=RS,LU=R*BSQ,LV = BSQ*SQRT(G)/R12
!     HERE, XS (X=Z,R) DO NOT INCLUDE DERIVATIVE OF EXPLICIT SQRT(S)
!     BSQ = |B|**2/2 + p        
!     GIJ = (BsupI * BsupJ) * SQRT(G)  (I,J = U,V)
!     IT IS ESSENTIAL THAT LU,LV AT j=1 ARE ZERO INITIALLY
!
!     SOME OF THE BIGGER LOOPS WERE SPLIT TO FACILITATE CACHE
!     HITS, PIPELINING ON RISCS
!
!     FOR OPTIMIZATION ON CRAY, MUST USE COMPILER DIRECTIVES TO
!     GET VECTORIZATION OF LOOPS INVOLVING POINTERS!
!
      dphids = p25/torflux_edge
      
      lu_e(1:ndim:ns) = 0; lv_e(1:ndim:ns) = 0
!RPRES      lu_o(1:ndim:ns) = 0
      guu(1:ndim:ns)  = 0; guv(1:ndim:ns)  = 0; gvv(1:ndim:ns) = 0
      guus = guu*shalf;    guvs = guv*shalf;    gvvs = gvv*shalf

CDIR$ IVDEP
      do l = 1, ndim
         armn_e(l)  = ohs*armn_e(l) * lu_e(l)
         azmn_e(l)  =-ohs*azmn_e(l) * lu_e(l)
         brmn_e(l)  = brmn_e(l) * lu_e(l)
         bzmn_e(l)  =-bzmn_e(l) * lu_e(l)
         bsqr(l)    = phip(l)*lu_e(l)/shalf(l)
      end do
      
CDIR$ IVDEP
      do l = 1, ndim
         armn_o(l)  = armn_e(l) *shalf(l)
         azmn_o(l)  = azmn_e(l) *shalf(l)
         brmn_o(l)  = brmn_e(l) *shalf(l)
         bzmn_o(l)  = bzmn_e(l) *shalf(l)
      end do   
! 
!     CONSTRUCT CYLINDRICAL FORCE KERNELS
!     NOTE: presg(ns+1) == 0, AND WILL BE "FILLED IN" AT EDGE
!     FOR FREE-BOUNDARY BY RBSQ
!
CDIR$ IVDEP
      do l = 1, nrzt
         guu(l) = p5*(guu(l) + guu(l+1))
         gvv(l) = p5*(gvv(l) + gvv(l+1))
         bsqr(l) = dphids*(bsqr(l) + bsqr(l+1))
         guus(l) = p5*(guus(l) + guus(l+1))
         gvvs(l) = p5*(gvvs(l) + gvvs(l+1))
!RPRES         presg(l)= ohs*(presg(l+1) - presg(l))*
!RPRES     1             (r1(l,0) + sqrts(l)*r1(l,1))                         !R*dp/ds
      end do

CDIR$ IVDEP
      do l = 1, nrzt
         armn_e(l) = armn_e(l+1) - armn_e(l) + p5*(lv_e(l) + lv_e(l+1))
     1             - gvv(l)*r1(l,0)
!RPRES     2             + presg(l)*zu0(l)
         azmn_e(l) = azmn_e(l+1) - azmn_e(l)
!RPRES                   - presg(l)*ru0(l)
         brmn_e(l) = p5*(brmn_e(l) + brmn_e(l+1))
         bzmn_e(l) = p5*(bzmn_e(l) + bzmn_e(l+1))
      end do

CDIR$ IVDEP
      do l = 1, nrzt
        armn_e(l) = armn_e(l) - gvvs(l)*r1(l,1)
        brmn_e(l) = brmn_e(l) + bsqr(l)*z1(l,1) 
     1            - guus(l)*ru(l,1) - guu(l)*ru(l,0) 
        bzmn_e(l) = bzmn_e(l) - bsqr(l)*r1(l,1) 
     1            - guus(l)*zu(l,1) - guu(l)*zu(l,0)
      end do  
!
!     ORIGIN OF VARIOUS TERMS
!
!     LU :  VARIATION OF DOMINANT .5*(RU-odd*Zodd - ZU-odd*Rodd) TERM
!           IN JACOBIAN
      
CDIR$ IVDEP
      do l = 1, nrzt
         armn_o(l) = armn_o(l+1) - armn_o(l) - zu(l,0)*bsqr(l) 
     1             + p5*(lv_e(l)*shalf(l) + lv_e(l+1)*shalf(l+1)) 
         azmn_o(l) = azmn_o(l+1) - azmn_o(l) + ru(l,0)*bsqr(l) 
         brmn_o(l) = p5*(brmn_o(l) + brmn_o(l+1))
         bzmn_o(l) = p5*(bzmn_o(l) + bzmn_o(l+1))         
      end do

!RPRES      presg(:nrzt)  = presg(:nrzt)*sqrts(:nrzt)
!RPRES      armn_o(:nrzt) = armn_o(:nrzt) + presg(:nrzt)*zu0(:nrzt)
!RPRES      azmn_o(:nrzt) = azmn_o(:nrzt) - presg(:nrzt)*ru0(:nrzt)

CDIR$ IVDEP
      do l = 1, nrzt
         lu_o(l)   = dphids*(lu_e(l)*phip(l) + lu_e(l+1)*phip(l+1))
         guu(l)    = guu(l)*sqrts(l)*sqrts(l)
         bsqr(l)   = gvv(l)*sqrts(l)*sqrts(l)
         armn_o(l) = armn_o(l) - zu(l,1)*lu_o(l)
     1             - bsqr(l)*r1(l,1) - gvvs(l)*r1(l,0)
         azmn_o(l) = azmn_o(l) + ru(l,1)*lu_o(l)
         brmn_o(l) = brmn_o(l) + z1(l,1)*lu_o(l)
     1             - guu(l)*ru(l,1) - guus(l)*ru(l,0)
         bzmn_o(l) = bzmn_o(l) - r1(l,1)*lu_o(l)
     1             - guu(l)*zu(l,1) - guus(l)*zu(l,0)
      end do      

      if (lthreed) then
CDIR$ IVDEP
         do l = 1, nrzt
            guv(l)  = p5*(guv(l) + guv(l+1))
            guvs(l) = p5*(guvs(l) + guvs(l+1))
            brmn_e(l) = brmn_e(l) - guv(l)*rv(l,0) - guvs(l)*rv(l,1)
            bzmn_e(l) = bzmn_e(l) - guv(l)*zv(l,0) - guvs(l)*zv(l,1)
            crmn_e(l) = guv(l) *ru(l,0) + gvv(l) *rv(l,0) 
     1                + gvvs(l)*rv(l,1) + guvs(l)*ru(l,1)
            czmn_e(l) = guv(l) *zu(l,0) + gvv(l) *zv(l,0)
     1                + gvvs(l)*zv(l,1) + guvs(l)*zu(l,1)
         end do

CDIR$ IVDEP
         do l = 1, nrzt
            guv(l) = guv(l) *sqrts(l)*sqrts(l)
            brmn_o(l) = brmn_o(l) - guvs(l)*rv(l,0) - guv(l)*rv(l,1)
            bzmn_o(l) = bzmn_o(l) - guvs(l)*zv(l,0) - guv(l)*zv(l,1)
            crmn_o(l) = guvs(l)*ru(l,0) + gvvs(l)*rv(l,0) 
     1                + bsqr(l)*rv(l,1) + guv(l) *ru(l,1) 
            czmn_o(l) = guvs(l)*zu(l,0) + gvvs(l)*zv(l,0) 
     1                + bsqr(l)*zv(l,1) + guv(l) *zu(l,1) 
         end do
      endif
!
!     ASSIGN EDGE FORCES (JS = NS) FOR FREE BOUNDARY CALCULATION
!
      if (ivac .ge. 1) then
        lk = 0
CDIR$ IVDEP
        do l = ns,nrzt,ns
          lk = lk+1
          armn_e(l) = armn_e(l) + zu0(l)*rbsq(lk)
          azmn_e(l) = azmn_e(l) - ru0(l)*rbsq(lk)
          armn_o(l) = armn_o(l) + zu0(l)*rbsq(lk)
          azmn_o(l) = azmn_o(l) - ru0(l)*rbsq(lk)
        end do  
      endif

      deallocate (bsqr, gvvs, guvs, guus, stat=l)
!
!     COMPUTE CONSTRAINT FORCE KERNELS
!
CDIR$ IVDEP
      do l = 1,nrzt
         rcon1   = (rcon(l,0) - rcon0(l)) * gcon(l,0)
         zcon1   = (zcon(l,0) - zcon0(l)) * gcon(l,0)
         brmn_e(l) = brmn_e(l) + rcon1
         bzmn_e(l) = bzmn_e(l) + zcon1
         brmn_o(l) = brmn_o(l)+ rcon1*sqrts(l)
         bzmn_o(l) = bzmn_o(l)+ zcon1*sqrts(l)
         rcon(l,0) =  ru0(l) * gcon(l,0)
         zcon(l,0) =  zu0(l) * gcon(l,0)
         rcon(l,1) = rcon(l,0) * sqrts(l)
         zcon(l,1) = zcon(l,0) * sqrts(l)
      end do   
 
      end subroutine forces

      subroutine tomnspa(frzl_array, armn, brmn, crmn, azmn, bzmn, 
     1   czmn, blmn, clmn, arcon, azcon)
      use vmec_main
      use vmec_params, only: jlam, jmin2, ntmax
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,0:ntor,0:mpol1,3*ntmax),
     1   intent(inout) :: frzl_array
      real(rprec), dimension(ns*nzeta,ntheta3,0:1), intent(in) :: 
     1   armn, brmn, crmn, azmn, bzmn, czmn, blmn, clmn, arcon, azcon
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: frcs = 3, frsc = 4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: fzcc, fzss,  flcc, flss
      integer :: jmax, m, mparity, i, jk, n, k, js
      real(rprec), dimension(ns,nzeta,12) :: work1
      real(rprec), dimension(ns*nzeta) :: temp1, temp3
C-----------------------------------------------
 
      fzcc = frcs + ntmax
      fzss = frsc + ntmax
      flcc = frcs + 2*ntmax
      flss = frsc + 2*ntmax
 
      jmax = ns
      if (ivac .lt. 1) jmax = ns1
!
!     BEGIN INVERSE FOURIER TRANSFORM
!     DO THETA (U) INTEGRATION FIRST
!
      do m = 0, mpol1
         mparity = mod(m,2)
         work1 = zero
         do i = 1, ntheta2
            do jk = 1, ns*nzeta
               temp1(jk) = armn(jk,i,mparity) + xmpq(m,1)*arcon(jk,i,
     1            mparity)
               temp3(jk) = azmn(jk,i,mparity) + xmpq(m,1)*azcon(jk,i,
     1            mparity)
               work1(jk,1,3) = work1(jk,1,3) + temp1(jk)*sinmui(i,m) + 
     1            brmn(jk,i,mparity)*cosmumi(i,m)
               work1(jk,1,5) = work1(jk,1,5) + temp3(jk)*cosmui(i,m) + 
     1            bzmn(jk,i,mparity)*sinmumi(i,m)
               work1(jk,1,9) = work1(jk,1,9) + blmn(jk,i,mparity)*
     1            sinmumi(i,m)
            end do
            if (lthreed) then
               do jk = 1, ns*nzeta
                  work1(jk,1,1) = work1(jk,1,1) + temp1(jk)*cosmui(i,m)
     1                + brmn(jk,i,mparity)*sinmumi(i,m)
                  work1(jk,1,7) = work1(jk,1,7) + temp3(jk)*sinmui(i,m)
     1                + bzmn(jk,i,mparity)*cosmumi(i,m)
                  work1(jk,1,11) = work1(jk,1,11) + blmn(jk,i,mparity)*
     1               cosmumi(i,m)
                  work1(jk,1,2) = work1(jk,1,2) - crmn(jk,i,mparity)*
     1               cosmui(i,m)
                  work1(jk,1,4) = work1(jk,1,4) - crmn(jk,i,mparity)*
     1               sinmui(i,m)
                  work1(jk,1,6) = work1(jk,1,6) - czmn(jk,i,mparity)*
     1               cosmui(i,m)
                  work1(jk,1,8) = work1(jk,1,8) - czmn(jk,i,mparity)*
     1               sinmui(i,m)
                  work1(jk,1,10) = work1(jk,1,10) - clmn(jk,i,mparity)*
     1               cosmui(i,m)
                  work1(jk,1,12) = work1(jk,1,12) - clmn(jk,i,mparity)*
     1               sinmui(i,m)
               end do
            endif
         end do
!
!        NEXT, DO ZETA (V) INTEGRATION
!
         do n = 0, ntor
            do k = 1, nzeta
 
               do js = jmin2(m), jmax
                  frzl_array(js,n,m,frsc) = frzl_array(js,n,m,frsc) + 
     1               work1(js,k,3)*cosnv(k,n)
                  frzl_array(js,n,m,fzcc) = frzl_array(js,n,m,fzcc) + 
     1               work1(js,k,5)*cosnv(k,n)
               end do
 
               do js = jlam(m), ns
                  frzl_array(js,n,m,flcc) = frzl_array(js,n,m,flcc) + 
     1               work1(js,k,9)*cosnv(k,n)
               end do
 
               if (lthreed) then
                  do js = jmin2(m), jmax
                     frzl_array(js,n,m,frsc) = frzl_array(js,n,m,frsc)
     1                   + work1(js,k,4)*sinnvn(k,n)
                     frzl_array(js,n,m,fzcc) = frzl_array(js,n,m,fzcc)
     1                   + work1(js,k,6)*sinnvn(k,n)
                     frzl_array(js,n,m,frcs) = frzl_array(js,n,m,frcs)
     1                   + work1(js,k,1)*sinnv(k,n) + work1(js,k,2)*
     2                  cosnvn(k,n)
                     frzl_array(js,n,m,fzss) = frzl_array(js,n,m,fzss)
     1                   + work1(js,k,7)*sinnv(k,n) + work1(js,k,8)*
     2                  cosnvn(k,n)
                  end do
                  do js = jlam(m), ns
                     frzl_array(js,n,m,flcc) = frzl_array(js,n,m,flcc)
     1                   + work1(js,k,10)*sinnvn(k,n)
                     frzl_array(js,n,m,flss) = frzl_array(js,n,m,flss)
     1                   + work1(js,k,11)*sinnv(k,n) + work1(js,k,12)*
     2                  cosnvn(k,n)
                  end do
               endif
 
            end do
         end do
      end do

      end subroutine tomnspa

      subroutine tomnsps(frzl_array, armn, brmn, crmn, azmn, bzmn, 
     1   czmn, blmn, clmn, arcon, azcon)
      use vmec_main
      use vmec_params, only: jlam, jmin2, ntmax
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,0:ntor,0:mpol1,3*ntmax),
     1   intent(out) :: frzl_array
      real(rprec), dimension(ns*nzeta,ntheta3,0:1), intent(in) :: 
     1   armn, brmn, crmn, azmn, bzmn, czmn, blmn, clmn, arcon, azcon
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: frcc = 1, frss = 2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: fzcs, fzsc, flcs, flsc
      integer :: jmax, m, mparity, i, n, k, js, l
      real(rprec), dimension(ns*nzeta,12) :: work1
      real(rprec), dimension(ns*nzeta) :: temp1, temp3
C-----------------------------------------------

      fzcs = frcc+ntmax
      fzsc = frss+ntmax
      flcs = frcc+2*ntmax
      flsc = frss+2*ntmax

      frzl_array = zero
      
      jmax = ns
      if (ivac .lt. 1) jmax = ns1

!
!     BEGIN INVERSE FOURIER TRANSFORM
!     DO THETA (U) INTEGRATION FIRST ON HALF INTERVAL (0 < U < PI)
!
!       FRmn = ARmn - d(BRmn)/du + d(CRmn)/dv
!       FZmn = AZmn - d(BZmn)/du + d(CZmn)/dv
!       FLmn =      - d(BLmn)/du + d(CLmn)/dv
!
!       NOTE: sinmumi = -m sin(mu),  sinnvn = -n sin(nv)
!
      do m = 0, mpol1
         mparity = mod(m,2)
         work1 = zero
         do i = 1, ntheta2
            temp1(:) = armn(:,i,mparity) + xmpq(m,1)*
     1            arcon(:,i,mparity)
            temp3(:) = azmn(:,i,mparity) + xmpq(m,1)*
     1            azcon(:,i,mparity)
            work1(:,1) = work1(:,1) + temp1(:)*cosmui(i,m) + 
     1            brmn(:,i,mparity)*sinmumi(i,m)
            work1(:,7) = work1(:,7) + temp3(:)*sinmui(i,m) + 
     1            bzmn(:,i,mparity)*cosmumi(i,m)
            work1(:,11) = work1(:,11) + blmn(:,i,mparity)*
     1            cosmumi(i,m)
            if (lthreed) then
               work1(:,2) = work1(:,2) - crmn(:,i,mparity)*
     1               cosmui(i,m)
               work1(:,4) = work1(:,4) - crmn(:,i,mparity)*
     1               sinmui(i,m)
               work1(:,3) = work1(:,3) + temp1(:)*sinmui(i,m)
     1                + brmn(:,i,mparity)*cosmumi(i,m)
               work1(:,5) = work1(:,5) + temp3(:)*cosmui(i,m)
     1                + bzmn(:,i,mparity)*sinmumi(i,m)
               work1(:,6) = work1(:,6) - czmn(:,i,mparity)*
     1               cosmui(i,m)
               work1(:,8) = work1(:,8) - czmn(:,i,mparity)*
     1               sinmui(i,m)
               work1(:,9) = work1(:,9) + blmn(:,i,mparity)*
     1               sinmumi(i,m)
               work1(:,10) = work1(:,10) - clmn(:,i,mparity)*
     1               cosmui(i,m)
               work1(:,12) = work1(:,12) - clmn(:,i,mparity)*
     1               sinmui(i,m)
            endif
         end do
!
!                NEXT, DO ZETA (V) INTEGRATION
!
         do n = 0, ntor
            do k = 1, nzeta
               l = ns*(k-1)
               do js = jmin2(m), jmax
                  frzl_array(js,n,m,frcc) = frzl_array(js,n,m,frcc) + 
     1               work1(js+l,1)*cosnv(k,n)
                  frzl_array(js,n,m,fzsc) = frzl_array(js,n,m,fzsc) + 
     1               work1(js+l,7)*cosnv(k,n)
               end do
               do js = jlam(m), ns
                  frzl_array(js,n,m,flsc) = frzl_array(js,n,m,flsc) + 
     1               work1(js+l,11)*cosnv(k,n)
               end do
 
               if (lthreed) then
                  do js = jmin2(m), jmax
                     frzl_array(js,n,m,frcc) = frzl_array(js,n,m,frcc)
     1                   + work1(js+l,2)*sinnvn(k,n)
                     frzl_array(js,n,m,fzsc) = frzl_array(js,n,m,fzsc)
     1                   + work1(js+l,8)*sinnvn(k,n)
                     frzl_array(js,n,m,frss) = frzl_array(js,n,m,frss)
     1                   + work1(js+l,3)*sinnv(k,n) + work1(js+l,4)*
     2                  cosnvn(k,n)
                     frzl_array(js,n,m,fzcs) = frzl_array(js,n,m,fzcs)
     1                   + work1(js+l,5)*sinnv(k,n) + work1(js+l,6)*
     2                  cosnvn(k,n)
                  end do
                  do js = jlam(m), ns
                     frzl_array(js,n,m,flsc) = frzl_array(js,n,m,flsc)
     1                   + work1(js+l,12)*sinnvn(k,n)
                     frzl_array(js,n,m,flcs) = frzl_array(js,n,m,flcs)
     1                   + work1(js+l,9)*sinnv(k,n) + work1(js+l,10)*
     2                  cosnvn(k,n)
                  end do
               endif
 
            end do
         end do
      end do

      end subroutine tomnsps

      subroutine jacobian
      use vmec_main, only: ohs, nrzt, irst
      use vmec_params, only: meven, modd
      use realspace
      use vmec_dim, only: ns
      use vforces, r12 => armn_o, ru12 => azmn_e, z12 => blmn_e,
     1    zu12 => armn_e, rs => bzmn_e, zs => brmn_e, gsqrt => azmn_o
      use vsvd, only: torflux_edge => torflux
      implicit none
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, p5 = 0.5_dp, p25 = p5*p5
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: l
      real(rprec) :: taumax, taumin, dphids
C-----------------------------------------------
 
!
!     (RS, ZS)=(R, Z) SUB S, (RU12, ZU12)=(R, Z) SUB THETA(=U)
!     AND GSQRT=SQRT(G) ARE DIFFERENCED ON HALF MESH
!     NOTE: LOOPS WERE SPLIT TO ALLOW EFFICIENT MEMORY BUS USAGE
!
!     FOR OPTIMIZATION ON CRAY, MUST USE COMPILER DIRECTIVES TO
!     GET VECTORIZATION OF LOOPS INVOLVING MORE THAN ONE POINTER!
!
!
!     HERE, SQRT(G) = R( Ru * Zs - Rs * Zu ). THE DERIVATIVES OF SHALF = SQRT(PHI(s))
!     WERE COMPUTED EXPLICITLY AS: d Shalf/ds = .5/shalf * d(PHI)/ds, WHERE
!     d(PHI)/ds = phip(s)/torflux_edge
!
!
      dphids = p25/torflux_edge
      
CDIR$ IVDEP
      do l = 2,nrzt
        ru12(l) = p5*(ru(l,meven) + ru(l-1,meven) + 
     1       shalf(l)*(ru(l,modd)  + ru(l-1,modd)))
        zs(l)   = ohs*(z1(l,meven) - z1(l-1,meven) +
     1       shalf(l)*(z1(l,modd)  - z1(l-1,modd)))
        z12(l)  = p5*(z1(l,meven) + z1(l-1,meven) +
     1       shalf(l)*(z1(l,modd)  + z1(l-1,modd)))
        gsqrt(l) = ru12(l)*zs(l) + dphids*phip(l)*
     1  (ru(l,modd) *z1(l,modd) + ru(l-1,modd) *z1(l-1,modd) +
     2  (ru(l,meven)*z1(l,modd) + ru(l-1,meven)*z1(l-1,modd))/shalf(l))
      enddo


CDIR$ IVDEP
      do l = 2,nrzt
        zu12(l) = p5*(zu(l,meven) + zu(l-1,meven) +
     1       shalf(l)*(zu(l,modd)  + zu(l-1,modd)))
        rs(l)   = ohs*(r1(l,meven) - r1(l-1,meven) +
     1       shalf(l)*(r1(l,modd)  - r1(l-1,modd)))
        r12(l)  = p5*(r1(l,meven) + r1(l-1,meven) +
     1       shalf(l)*(r1(l,modd)  + r1(l-1,modd)))
        gsqrt(l) = r12(l)*(gsqrt(l)- rs(l)*zu12(l) - dphids*phip(l)*
     1    (zu(l,modd) *r1(l,modd)+zu(l-1,modd) *r1(l-1,modd)
     2  + (zu(l,meven)*r1(l,modd)+zu(l-1,meven)*r1(l-1,modd))/shalf(l)))
      end do

!
!     TEST FOR SIGN CHANGE IN JACOBIAN
!
      gsqrt(1:nrzt:ns) = gsqrt(2:nrzt:ns)
      taumax = maxval(gsqrt(2:nrzt))
      taumin = minval(gsqrt(2:nrzt))
      if (taumax*taumin .lt. zero) irst = 2

      end subroutine jacobian
      

      subroutine allocate_funct3d
      use vmec_main
      use realspace
      use vforces
      use vacmod
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat1, ndim, ndim2
C-----------------------------------------------

      ndim  = 1+nrzt
      ndim2 = 2*ndim

      call free_mem_funct3d

      allocate( armn(ndim2), azmn(ndim2), brmn(ndim2), bzmn(ndim2), 
     1   crmn(ndim2), czmn(ndim2), blmn(ndim2), clmn(ndim2), 
     2   r1(nrzt,0:1), ru(nrzt,0:1), rv(nrzt,0:1),
     3   z1(nrzt,0:1), zu(nrzt,0:1), zv(nrzt,0:1),
     4   rcon(nrzt,0:1), zcon(nrzt,0:1), ru0(ndim), zu0(ndim),
     6   rcon0(ndim), zcon0(ndim), guu(ndim), guv(ndim), gvv(ndim),
     7   stat=istat1 )
      if (istat1.ne.0) stop 'allocation error #1 in funct3d'

      if (lfreeb) then
      allocate (brv(nznt), bphiv(nznt), bzv(nznt), 
     1   bpolvac(nznt), bsqvac(nznt), stat=istat1) 
      if (istat1.ne.0) stop 'allocation error #2 in funct3d'
      end if
!
!     Pointer alias assignments 
!     NOTE: In FORCES, X_e(nrzt+1) overlaps X_o(1), which should never be used...
!      
      armn_e => armn(:ndim)
      armn_o => armn(ndim:)
      armn(:ndim2) = zero
      brmn_e => brmn(:ndim)
      brmn_o => brmn(ndim:)
      brmn(:ndim2) = zero
      azmn_e => azmn(:ndim)
      azmn_o => azmn(ndim:)
      azmn(:ndim2) = zero
      bzmn_e => bzmn(:ndim)
      bzmn_o => bzmn(ndim:)
      bzmn(:ndim2) = zero
      crmn_e => crmn(:ndim)
      crmn_o => crmn(ndim:)
      crmn(:ndim2) = zero
      czmn_e => czmn(:ndim)
      czmn_o => czmn(ndim:)
      czmn(:ndim2) = zero
      blmn_e => blmn(:ndim)
      blmn_o => blmn(ndim:)
      blmn(:ndim2) = zero
      clmn_e => clmn(:ndim)
      clmn_o => clmn(ndim:)
      clmn(:ndim2) = zero
      rcon0(:ndim) = zero
      zcon0(:ndim) = zero
                  
      end subroutine allocate_funct3d
      

      subroutine allocate_ns (lreset, linterp, neqs2_old)
      use vmec_main
      use vmec_params, only: ntmax
      use realspace
      use vsvd
      use vspline
      use vforces
      use xstuff
      use csplinx
      implicit none
C-----------------------------------------------
C   D u m m y   V a r i a b l e s
C-----------------------------------------------
      integer neqs2_old
      logical :: lreset, linterp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ndim, nsp1, istat1
      real(rprec), dimension(:), allocatable :: xc_old, scalxc_old
      real(rprec) delr_mse
C-----------------------------------------------
!
!     FIRST STORE COARSE-MESH XC FOR INTERPOLATION
!

      ndim  = 1 + nrzt
      nsp1  = 1 + ns
      delr_mse = zero

!
!     Make consistency checks (just for extra protection)
!
      if (linterp .and. (neqs2_old .eq. neqs2)) then
         print *,' Setting linterp = F in allocate_ns. '
         linterp = .false.
      end if

!
!     Save old xc, scalxc for possible interpolation or if iterations restarted on same mesh...
!            
      if (neqs2_old .gt. 0 .and. allocated(scalxc) .and. linterp) then
         allocate(xc_old(neqs2_old), scalxc_old(neqs2_old), stat=istat1)
         if (istat1.ne.0) stop 'allocation error #1 in allocate_ns'
         xc_old(:neqs2_old) = xc(:neqs2_old)
         scalxc_old(:neqs2_old) = scalxc(:neqs2_old)
         if (lrecon) delr_mse = xc(neqs2_old)
      end if

!
!     ALLOCATES MEMORY FOR NS-DEPENDENT ARRAYS
!     FIRST BE SURE TO FREE MEMORY PREVIOUSLY ALLOCATED
!
      call free_mem_ns (lreset)

      allocate (phip(ndim), shalf(ndim), sqrts(ndim), wint(ndim), 
     1  stat=istat1)
      if (istat1.ne.0) stop 'allocation error #2 in allocate_ns'
      allocate( ireflect(ns*nzeta), indexr(2*ns) ,imid(2*ns),  
     1  stat=istat1)
      if (istat1.ne.0) stop 'allocation error #3 in allocate_ns'
      allocate( current(ns),rm2(ns),vrm2(ns),
     1  ovrm2(ns), ochip(ns), presph(ns), presint(ns),
     2  w_ia(ns), w1_ia(ns), u_ia(ns), u1_ia(ns),
     3  w_pa(ns), w1_pa(ns), u_pa(ns), u1_pa(ns),
     4  w_ib(ns), w1_ib(ns), u_ib(ns), u1_ib(ns),
     5  w_pb(ns), w1_pb(ns), u_pb(ns), u1_pb(ns),
     6  rmid(2*ns),datamse(2*ns),qmid(2*ns),
     7  shear(2*ns),presmid(2*ns),alfa(2*ns),curmid(2*ns),
     8  curint(2*ns),psimid(2*ns),ageo(2*ns),volpsi(2*ns),
     9  isplinef(ns),isplineh(ns),psplinef(ns),psplineh(ns),
     A  phimid(2*ns),pm(ns,0:nobser+nobd),
     B  im(ns,0:nobser+nobd),stat=istat1)
      if (istat1.ne.0) stop 'allocation error #4 in allocate_ns'

      if (nbsets.gt.0)
     1  allocate( pmb(ns,0:nbcoil_max,nbsets,2),
     2            imb(ns,0:nbcoil_max,nbsets,2),stat=istat1)
      if (istat1.ne.0) stop 'allocation error #5 in allocate_ns'

      allocate( ard(nsp1,2),arm(nsp1,2),brd(nsp1,2),brm(nsp1,2),
     1          azd(nsp1,2),azm(nsp1,2),bzd(nsp1,2), bzm(nsp1,2),
     2          sm(ns), sp(0:ns), bmin(ntheta2,ns), bmax(ntheta2,ns),
     3          stat=istat1)
      if (istat1.ne.0) stop 'allocation error #6 in allocate_ns'

      allocate( iotaf(nsp1), crd(nsp1), mass(ns), phi(ns), presf(ns),
     1          jcuru(ns), jcurv(ns), jdotb(ns), buco(ns), bvco(ns),
     2          bdotgradv(ns), equif(ns), specw(ns), tcon(ns), fpsi(ns),
     3          psi(ns),yellip(ns),yinden(ns), ytrian(ns),yshift(ns),
     4          ygeo(ns),overr(ns), faclam(ns,0:ntor,0:mpol1,ntmax),
     5          iotas(nsp1), phips(nsp1), pres(nsp1), vp(nsp1), 
     6          beta_vol(ns), jperp2(ns), jpar2(ns), bdotb(ns),
     7          blam(nsp1), clam(nsp1), dlam(nsp1), 
     7          stat=istat1)
      if (istat1.ne.0) stop 'allocation error #7 in allocate_ns'

      allocate( rmidx(2*ns), hmidx(2*ns), wmidx(2*ns), qmidx(2*ns),
     1          tenmidx(2*ns), ymidx(2*ns), y2midx(2*ns), stat=istat1)
      if (istat1.ne.0) stop 'allocation error #8 in allocate_ns'
     
      allocate (gc(neqs2), xcdot(neqs2), xstore(neqs2), stat=istat1)
      if (istat1.ne.0) stop 'allocation error #9 in allocate_ns'
      
      if (lreset) then
         allocate (xc(neqs2), scalxc(neqs2), stat=istat1)
      if (istat1.ne.0) stop 'allocation error #10 in allocate_ns'
         xc(:neqs2) = zero
      end if         
   
      if (allocated(xc_old)) then
         xstore(1:neqs2_old) = xc_old(1:neqs2_old)
         scalxc(1:neqs2_old) = scalxc_old(1:neqs2_old)
         deallocate (xc_old, scalxc_old)         
      end if
      

      xc(neqs2) = delr_mse

!
!     Allocate nrzt-dependent arrays (persistent) for funct3d
!
      call allocate_funct3d

      end subroutine allocate_ns
      

      subroutine allocate_nunv
      use vmec_main
      use vmec_params, only: ntmax
      use vacmod
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat1
C-----------------------------------------------

      call free_mem_nunv

      allocate (bsubu0(nznt), rbsq(nznt), dbsq(nznt), stat=istat1) 
      if (istat1.ne.0) stop 'allocation error #1 in allocate_nunv'

      allocate (rmn_bdy(0:ntor,0:mpol1,ntmax),
     1          zmn_bdy(0:ntor,0:mpol1,ntmax), stat=istat1)
      if (istat1.ne.0) stop 'allocation error #2 in allocate_nunv'

!     PERSISTENT ARRAYS (DURATION OF PROGRAM)
      if (lfreeb) 
     1   allocate (amatsav(mnpd2*mnpd2),bvecsav(mnpd2), 
     2          bsqsav(nznt,3), potvac(2*mnpd), stat=istat1)
         if (istat1.ne.0) stop 'allocation error #3 in allocate_nunv'

      end subroutine allocate_nunv
        

      function aspectratio ()
      use vmec_main
      use realspace
      use vmec_io
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: lk, l
      real(rprec) :: rb, zub, pi, t1, aspectratio
C-----------------------------------------------
!
!     routine for computing aspect-ratio
!     A = <R>/<a>
!
!     where pi <a>**2 = Area (toroidally averaged)
!           2*pi * <R> * Area = Volume
!     Use integration by parts to compute as surface integral (Stoke''s theorem)
!
 
      pi = 4*atan(one)
 
!
!     Compute Volume and Area
!
      volume_plasma = 0
      cross_area = 0
      do lk = 1, nznt
         l = ns*lk
         rb  = r1(l,0) + r1(l,1)
         zub = zu(l,0) + zu(l,1)
         t1  = rb*zub*wint(l)
         volume_plasma = volume_plasma + rb*t1
         cross_area = cross_area + t1
      end do
 
      volume_plasma = 2*pi*pi*abs(volume_plasma)
      cross_area = 2*pi*abs(cross_area)
 
      Rmajor_p = volume_plasma/(2*pi*cross_area)
      Aminor_p = sqrt(cross_area/pi)

      aspectratio = Rmajor_p/Aminor_p

      end function aspectratio 

                 
      subroutine open_output_files (extension, iseq, lmac, lscreen,
     1           lfirst)
      use safe_open_mod
      use vparams, only: nmac, nthreed, nmac0, nthreed0
      implicit none      
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: iseq
      character*(*) :: extension
      logical :: lmac, lscreen, lfirst
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: iread, inthreed=0, imac0=0
      character*120 :: mac_file, threed1_file
C-----------------------------------------------
 
!
!     OPEN FILES FOR READING, WRITING
!
      threed1_file = 'threed1.'//extension
      mac_file = 'mac.'//extension
 
      if (lscreen .and. lfirst) write (*, '(33('' -''))')
      nthreed = nthreed0
      call safe_open(nthreed, iread, threed1_file, 'new', 'formatted')
      if (iread .ne. 0) then
         if (iseq .eq. 0 .and. lscreen .and. lfirst) print *, 
     1   ' VMEC OUTPUT FILES ALREADY EXIST: OVERWRITING THEM ...'
         call safe_open(nthreed, inthreed, threed1_file, 'replace',
     1     'formatted')
      endif
   
 
      nmac = max(nmac0, nthreed)
      if (lmac) then
         call safe_open(nmac, imac0, mac_file, 'replace', 'formatted')
      end if
      if (inthreed.ne.0 .or. imac0.ne.0) then
         print *,' nthreed = ', nthreed, ' istat_threed = ', inthreed,
     1           ' nmac0   = ', nmac,' istat_mac0 = ', imac0 
         print *, 'Error opening output file in VMEC OPEN_OUTPUT_FILES'
         stop 10
      endif

      end subroutine open_output_files


      subroutine close_all_files
      use vparams, only: nmac, nthreed
      implicit none
C-----------------------------------------------

      if (nthreed .gt. 0) close (nthreed)
      if (nmac .gt. 0) close (nmac)

      end subroutine close_all_files
      

      subroutine storesvd(re, ro, lue, luo, lve, phipog, zu00)
      use vmec_main
      use realspace
      use vsvd
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,nzeta,*) :: re, ro, lue
      real(rprec), dimension(ns,*) :: luo
      real(rprec), dimension(*) :: lve, phipog
      real(rprec), dimension(ns,nzeta,*) :: zu00
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: lk, l, noff, jmin, lt, js, j, ks, j1, j2
      real(rprec), dimension(ns,ntheta2) :: luef, luof
      real(rprec) :: signi, phipog0, fpsie, w1,
     1    fpsi2, fpsi1, diota
C-----------------------------------------------
!
!     COMPUTE FULL MESH VALUES OF LAMBDA AT THETA=0,PI
!
      call lamfull (luef, luof, lue, luo)
!
!     COMPUTE PHIP*GUU*(IOTA-LV)/GSQRT AS S->0 FOR CURRENT(0)
!
        do lk = 1,nznt
          l = 2+ns*(lk-1)
          phipog0 = 1.5_dp*phipog(l) - 0.5_dp*phipog(l+1)
          bsubu0(lk) = phipog0*( (iotas(2) + lve(l))*guu(l)
     1    + shalf(2)*luo(2,lk)*guv(l) )
        enddo
 
      if (.not.lrecon) return
      
      signi = sign(one,iotaf(2))
      noff = ntheta2
      jmin = 2
      do lt = 1, 2
         j2   = ns*nzeta*(noff-1)
         if (lt .eq. 1) then
            l = ns+1-jmin
            imid(l:1:(-1)) = (/(j1,j1=jmin+j2,ns+j2)/)
            rmid(l:1:(-1)) = re(jmin:ns,1,noff) + sqrts(jmin:ns)
     1         *ro(jmin:ns,1,noff)
            datamse(l:1:(-1)) = atan(iotaf(jmin:ns)*zu00(jmin:ns
     1         ,1,noff)/(rmid(l:1:(-1))*(luef(jmin:ns,noff)+
     2         sqrts(jmin:ns)*luof(jmin:ns,noff))))/dcon
            qmid(l:1:(-1)) = signi/iotaf(jmin:ns)
            presmid(l:1:(-1)) = presf(jmin:ns)/dmu0
         else
            l = ns+jmin-1
            imid(l:2*ns-1)=(/(j1,j1=jmin+j2,ns+j2)/)
            rmid(l:2*ns-1) = re(jmin:ns,1,noff) + sqrts(jmin:ns)
     1         *ro(jmin:ns,1,noff)
            datamse(l:2*ns-1) = atan(iotaf(jmin:ns)*zu00(jmin:ns
     1         ,1,noff)/(rmid(l:2*ns-1)*(luef(jmin:ns,noff)+
     2         sqrts(jmin:ns)*luof(jmin:ns,noff))))/dcon
            qmid(l:2*ns-1) = signi/iotaf(jmin:ns)
            presmid(l:2*ns-1) = presf(jmin:ns)/dmu0
         endif
         noff = 1
         jmin = 1
      end do
 
      call findphi (re, ro, rthom, delse2, delso2, rmid, indexs2, 
     1   indexu2, indexr, itse)
      pcalc(:itse) = (presf(indexs2(:itse))) + abs(delse2(:itse))*(presf
     1   (indexs2(:itse)+1)-(presf(indexs2(:itse))))
 
!
!       SORT ON RSORT(KS) ARRAY
!       INCLUDE EDGE PITCH MATCH TO TOTAL CURRENT
!
      fpsie = 1.5_dp*fpsi(ns) - 0.5_dp*fpsi(ns1)
 
      do j = 1, imse2
         ks = isorts(j)
         js = indexs1(ks)
         w1 = delse1(ks)
         if (js .eq. ns) then
            fpsi2 = fpsie
            fpsi1 = fpsie
         else if (js .eq. ns1) then
            fpsi2 = fpsie
         else
            fpsi2 = 0.5_dp*(fpsi(js+1)+fpsi(js+2))
         endif
         if (js.lt.ns .and. js.gt.1) fpsi1 = .5_dp*(fpsi(js)+fpsi(js+1))
         if (js .eq. 1) fpsi1 = fpsi(2)
         diota = (one - w1)*iotaf(js) + w1*iotaf(js+1)
         fpsical(ks) = (one - w1)*fpsi1 + w1*fpsi2
         if (stark_weight(ks) .ne. zero) qmeas(ks) = datastark(ks)/
     1      stark_weight(ks)
         qcalc(ks) = diota
         starkcal(ks) = diota*stark_weight(ks)
         rsort0(ks) = sqrt(hs*(js - 1 + w1))
      end do

      end subroutine storesvd

      
      subroutine bextrema(modb, bmin, bmax, nzeta, ntheta)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: nzeta, ntheta
      real(rprec) :: modb(nzeta,ntheta), bmin(ntheta), bmax(ntheta)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer ku
C-----------------------------------------------
!
!     Computes max, min of |B| along v (zeta) between two angle lines (theta = 0,pi)
!
      do ku = 1,ntheta
         bmin(ku)  = minval(modb(:,ku))
         bmax(ku)  = maxval(modb(:,ku))
      enddo

      end subroutine bextrema


      subroutine convert(rmnc,zmns,lmns,rmns,zmnc,lmnc,rzl_array,js)
      use vmec_main
      use vmec_params, only: mscale, nscale, ntmax
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer js
      real(rprec), dimension(mnmax), intent(out) ::
     1    rmnc, zmns, lmns, rmns, zmnc, lmnc
      real(rprec), dimension(ns,0:ntor,0:mpol1,3*ntmax),
     1    intent(in) :: rzl_array
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: rmncc=1
      integer, parameter :: rmnss=2
      integer, parameter :: rmncs=3
      integer, parameter :: rmnsc=4
      real(rprec), parameter :: p5 = 0.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: zmncs, zmnsc, zmncc, zmnss, lmncs, lmnsc, lmncc, lmnss
      integer :: mn, m, n, n1
      real(rprec) :: t1, sign0
C-----------------------------------------------
!
!     CONVERTS INTERNAL MODE REPRESENTATION TO STANDARD
!     FORM FOR OUTPUT (COEFFICIENTS OF cos(mu-nv), sin(mu-nv))
!
      zmncs = rmncc + ntmax
      lmncs = rmncc + 2*ntmax
      zmnsc = rmnss + ntmax
      lmnsc = rmnss + 2*ntmax
      zmncc = rmncs + ntmax
      lmncc = rmncs + 2*ntmax
      zmnss = rmnsc + ntmax
      lmnss = rmnsc + 2*ntmax

      mn = 0
!
!     DO M = 0 MODES SEPARATELY (ONLY KEEP N >= 0 HERE: COS(-NV), SIN(-NV))
!
      m = 0
      do n = 0, ntor
         t1 = mscale(m)*nscale(n)
         mn = mn + 1
         rmnc(mn) = t1*rzl_array(js,n,m,rmncc)
         zmns(mn) =-t1*rzl_array(js,n,m,zmncs)
         lmns(mn) =-t1*rzl_array(js,n,m,lmncs)
         if (lasym) then
         rmns(mn) = t1*rzl_array(js,n,m,rmncs)
         zmnc(mn) = t1*rzl_array(js,n,m,zmncc)
         lmnc(mn) = t1*rzl_array(js,n,m,lmncc)
         end if
      end do

      if (js .eq. 1) then
         mn = 0
         do n = 0, ntor
            t1 = mscale(m)*nscale(n)
            mn = mn + 1
            lmns(mn) =-t1*(2._dp*rzl_array(2,n,m,lmncs)
     1               -           rzl_array(3,n,m,lmncs))            
         end do   
      end if
       
      do m = 1, mpol1
         do n = -ntor, ntor
            n1 = abs(n)
            t1 = mscale(m)*nscale(n1)
            mn = mn + 1
            if (n .eq. 0) then
               rmnc(mn) = t1*rzl_array(js,n,m,rmncc)
               zmns(mn) = t1*rzl_array(js,n,m,zmnsc)
               lmns(mn) = t1*rzl_array(js,n,m,lmnsc)
               if (lasym) then
               rmns(mn) = t1*rzl_array(js,n,m,rmnsc)
               zmnc(mn) = t1*rzl_array(js,n,m,zmncc)
               lmnc(mn) = t1*rzl_array(js,n,m,lmncc)
               end if
            else if (js .gt. 1) then
               sign0 = n/n1
               rmnc(mn) = p5*t1*(rzl_array(js,n1,m,rmncc)+sign0*
     1            rzl_array(js,n1,m,rmnss))
               zmns(mn) = p5*t1*(rzl_array(js,n1,m,zmnsc)-sign0*
     1            rzl_array(js,n1,m,zmncs))
               lmns(mn) = p5*t1*(rzl_array(js,n1,m,lmnsc)-sign0*
     1            rzl_array(js,n1,m,lmncs))
               if (lasym) then
               rmns(mn) = p5*t1*(rzl_array(js,n1,m,rmncs)-sign0*
     1            rzl_array(js,n1,m,rmnsc))
               zmnc(mn) = p5*t1*(rzl_array(js,n1,m,zmncc)+sign0*
     1            rzl_array(js,n1,m,zmnss))
               lmnc(mn) = p5*t1*(rzl_array(js,n1,m,lmncc)+sign0*
     1            rzl_array(js,n1,m,lmnss))
               end if
            else if (js .eq. 1) then
               rmnc(mn) = zero
               zmns(mn) = zero
               lmns(mn) = zero
               if (lasym) then
               rmns(mn) = zero
               zmnc(mn) = zero
               lmnc(mn) = zero
               end if
            end if
         end do
      end do

      if( .not. lasym) then
         rmns = 0
         zmnc = 0
         lmnc = 0
      endif
 
      end subroutine convert
      

      subroutine eqsolve(nsval, interp_flag, ier_flag, 
     1    lreset, lfirst, lscreen, reset_file_name)
      use vmec_main
      use vmec_params, only: ntmax, ns4
      use realspace
      use vsvd
      use xstuff
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: nsval, ier_flag
      character*(*) :: reset_file_name
      logical :: interp_flag, lreset, lfirst, lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: bad_jacobian_flag = 6
      real(rprec), parameter :: p98 = 0.98_dp, p96 = 0.96_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nsold=0, neqs2_old=0
      real(rprec) :: w1, r00s, w0, res0, wdota, r0dot
      real(rprec) :: delt0
      logical :: liter_flag, lreset_internal
C-----------------------------------------------
!
!                INDEX OF LOCAL VARIABLES
!
!        hs      radial mesh size increment
!        iequi   counter used to call -EQFOR- at end of run
!        irzloff offset in xc array between R,Z,L components
!        ijacob  counter for number of times jacobian changes sign
!        irst    counter monitoring sign of jacobian; resets R, Z, and
!                Lambda when jacobian changes sign and decreases time step
!        signgs  sign of Jacobian : must be =1 (right-handed) or =-1 (left-handed)
!        iterj   stores position in main iteration loop (j=1,2)
!        itfsq   counter for storing FSQ into FSQT for plotting
!        ivac    counts number of free-boundary iterations
!        ndamp   number of iterations over which damping is averaged
!        meven   parity selection label for even poloidal modes of R and Z
!        modd    parity selection label for odd poloidal modes of R and
!        gc      stacked array of R, Z, Lambda Spectral force coefficients (see readin for stack order)
!        xc      stacked array of scaled R, Z, Lambda Fourier coefficients
 
!
!       INITIALIZE MESH-DEPENDENT SCALARS
!
      ns = nsval
      ns1 = ns - 1
      hs = one/ns1
      ohs = one/hs
      mns = ns*mnsize
      irzloff = ntmax*mns
      nrzt = nznt*ns
      neqs = 3*irzloff
      neqs1 = neqs + 1
      neqs2 = neqs1 + 1
      ijacob = 0
      delt0   = delt

      fsqt = 0                              ! mcz
      wdot = 0

      write (nthreed, 10) ns, mnmax, ftolv
      if (lscreen) print 10, ns, mnmax, ftolv
   10 format(/'  NS = ',i4,' NO. FOURIER MODES = ',i4,' FTOLV = ',
     1   1pe10.3)
      if (imovephi .gt. 0 .and. lscreen) print *, 
     1   'Turning on Ip Matching by varying Phi-Edge.'
 
!
!     ALLOCATE NS-DEPENDENT ARRAYS
!
      call allocate_ns(lreset, interp_flag, neqs2_old)

      lreset_internal = (ns .gt. nsold)

      if (lreset_internal .and. neqs2_old.gt.0) gc(1:neqs2_old) =
     1  scalxc(1:neqs2_old)*xstore(1:neqs2_old)

      lreset_internal = lreset_internal .or. lreset

 1000 continue
 
      itfsq = 0
      fsq     = one
      rsfac   = one
      w1      = zero
      r00s    = zero
      gphifac = zero
      grmse   = zero
 
!
!     COMPUTE INITIAL R, Z AND MAGNETIC FLUX PROFILES
!
   20 continue
      iter2 = 1
      irst = 1
      call profil1d (xc, xcdot, lreset)
      if (lfirst .and. .not.lreset)               !!Command line with lreset=F
     1    call load_xc_from_wout(xc(1), xc(1+irzloff), xc(1+2*irzloff),
     2         ntor, mpol1, ns, reset_file_name, lreset_internal)     
      call profil3d (xc(1), xc(1+irzloff), lreset_internal)
 
!
!     INTERPOLATE FROM COARSE TO NEXT FINER RADIAL GRID
!
      if (interp_flag) call interp (xc, gc, scalxc, ns, nsold)
      nsold = ns
      neqs2_old = neqs2
 
!
!     Store XC,XCDOT for possible restart
!
      call restart(delt0)
      iter1 = iter2
      liter_flag = .true.
      ier_flag = 0

!
!     ENTER FORCE ITERATION LOOP
!
      iter_loop: do while (liter_flag)
!
!     ADVANCE FOURIER AMPLITUDES OF R, Z, AND LAMBDA
!
         call evolve (delt0, ier_flag, liter_flag, lscreen)
         if (ijacob.eq.0 .and. ier_flag.eq.1) then
            if (lscreen .and. irst.eq.2) print *,
     1      ' INITIAL JACOBIAN CHANGED SIGN: IMPROVING THE',
     2      ' GUESS FOR THE MAGNETIC AXIS'
            if (lscreen .and. irst.eq.4) print *,
     1      ' INITIAL FORCES LARGE: IMPROVING JACOBIAN ESTIMATE'
            call guess_axis (r1, z1, ru0, zu0, lscreen)
            lreset_internal = .true.
            ijacob = 1
            go to 20
         else if (ier_flag .ne. 0) then
            if (ier_flag .eq. 1) ier_flag = bad_jacobian_flag
            return
         endif

         w0 = wb + wp/(gamma - one)
 
!
!     ADDITIONAL STOPPING CRITERION (set liter_flag to FALSE)
!
         if (ijacob .eq. 25) then
            irst = 2
            call restart(delt0)
            delt0 = p98*delt
            if (lscreen) print 120, delt0
            go to 1000
         else if (ijacob .eq. 50) then
            irst = 2
            call restart(delt0)
            delt0 = p96*delt
            if (lscreen) print 120, delt0
            go to 1000
         else if (ijacob .ge. 75) then
            ier_flag = 2
            liter_flag = .false.
         else if (iter2.ge.niter .and. liter_flag) then
            ier_flag = 4
            liter_flag = .false.
         endif
 
!
!       TIME STEP CONTROL
!
         if (iter2 .eq. iter1) res0 = fsq
         res0 = min(res0,fsq)
!       Store current state (irst=1)
         if (fsq .le. res0 .and. iter2-iter1 .gt. 10) then
            call restart(delt0)
!       Residuals are growing in time, reduce time step
         else if (fsq .gt. 100.0_dp*res0 .and. iter2 .gt. iter1) then
            irst = 2
         else if (iter2 - iter1 .gt. ns4/2 .and. iter2 .gt. 2*ns4 
     1        .and. fsqr+fsqz .gt. c1pm2) then
            irst = 3
         endif
 
         if (irst .ne. 1) then
!       Retrieve previous good state
            call restart(delt0)
            iter1 = iter2
         else
!       Increment time step and Printout every nstep iterations
            if (mod(iter2,nstep).eq.0 .or. iter2.eq.1 .or.
     1         .not.liter_flag) call printout(iter2, delt0, w0, lscreen)
            iter2 = iter2 + 1
         endif
 
!       Store force residual, wdot for plotting
         wdota = abs(w0 - w1)/w0
         r0dot = abs(r00 - r00s)/r00
         r00s = r00
         w1 = w0
         if (ivac.eq.1 .and. lreset) then
            if (lscreen) print 110, iter2
            write (nthreed, 110) iter2
            ivac = ivac + 1
         endif
!
!       STORE FSQ FOR PLOTTING. EVENTUALLY, STORE FOR EACH RADIAL MESH
!
         if (mod(iter2,niter/nstore_seq + 1).eq.0 .and. ns.eq.
     1      ns_array(multi_ns_grid)) then
            if (itfsq .lt. nstore_seq) then
              itfsq = itfsq + 1
              fsqt(itfsq) = fsqr + fsqz
              wdot(itfsq) = max(wdota,c1pm13)
            end if  
         end if
 
      end do iter_loop
 
      write (nthreed, 60) wdota, r0dot
      write (nthreed, 70) 1.e-6_dp*ctor/dmu0, rbtor, rbtor0
      if (lrecon) write (nthreed, 65) r00*fsqsum0/wb
 
   60 format(/,' d(ln W)/dt = ',1pe10.3,' d(ln R0)/dt = ',1pe10.3,/)
   65 format(' Average radial force balance: Int[FR(m=0)]',
     1   '/Int(B**2/R) = ',1pe12.5,' (should tend to zero)'/)
   70 format(' TOROIDAL CURRENT = ',1pe10.2,' [MA] ',' R-BTOR(s=1) = ',
     1   1pe10.2,' R-BTOR(s=0) = ',1pe10.2)
  110 format(/,2x,'VACUUM PRESSURE TURNED ON AT ',i4,' ITERATIONS'/)
  120 format(2x,'HAVING A CONVERGENCE PROBLEM: RESETTING DELT TO ',f8.3,
     1  /,2x,'If this does NOT resolve problem, inform vendor ',
     2       'that lamcal scaling needs to be checked')

      end subroutine eqsolve


      subroutine guess_axis(r1, z1, ru0, zu0, lscreen)
      use vmec_main
      use vmec_params, only: nscale, signgs
      use realspace, only: sqrts
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,nzeta,ntheta3,0:1),
     1     intent(in) :: r1, z1
      real(rprec), dimension(ns,nzeta,ntheta3), intent(in) :: ru0, zu0
      logical :: lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: limpts = 61
      real(rprec), parameter :: p5 = 0.5_dp, two = 2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: iv, iu, iu_r, ivminus, nlim, ns12, klim, n
      real(rprec), dimension(nzeta) :: rcom, zcom
      real(rprec), dimension(ntheta1) :: r1b, z1b, rub, zub
      real(rprec), dimension(ntheta1) :: r12, z12
      real(rprec), dimension(ntheta1) :: rs, zs, tau, ru12, zu12, tau0
      real(rprec) :: rlim(limpts), zlim
      real(rprec) :: rmax, rmin, zmax, zmin, dzeta
      real(rprec) :: ds, mintau, mintemp
C-----------------------------------------------
!
!     COMPUTES GUESS FOR MAGNETIC AXIS IF USER GUESS
!     LEADS TO INITIAL SIGN CHANGE OF JACOBIAN. DOES A GRID
!     SEARCH (irgrid, izgrid) IN EACH PHI-PLANE FOR POINTS WHICH
!     YIELD A VALUE FOR THE JACOBIAN WITH THE CORRECT SIGN (SIGNGS)
!     CHOOSES THE AXIS POSITION SO THE MIN VALUE OF THE JACOBIAN IS MAXIMIZED
!
      ns12 = (ns+1)/2

      planes: do iv = 1, nzeta
         if (.not.lasym .and. iv .gt. nzeta/2+1) then
            rcom(iv) = rcom(nzeta+2-iv)
            zcom(iv) =-zcom(nzeta+2-iv)
            cycle
         end if
         r1b(:ntheta3) = r1(ns,iv,:,0) + r1(ns,iv,:,1)
         z1b(:ntheta3) = z1(ns,iv,:,0) + z1(ns,iv,:,1)
         r12(:ntheta3) = r1(ns12,iv,:,0) + r1(ns12,iv,:,1)*sqrts(ns12)
         z12(:ntheta3) = z1(ns12,iv,:,0) + z1(ns12,iv,:,1)*sqrts(ns12)
         rub(:ntheta3) = ru0(ns,iv,:)
         zub(:ntheta3) = zu0(ns,iv,:)
         ru12(:ntheta3) =  p5*(ru0(ns,iv,:) + ru0(ns12,iv,:))
         zu12(:ntheta3) =  p5*(zu0(ns,iv,:) + zu0(ns12,iv,:))

         if (.not.lasym) then
!
!     USE Z(v,-u) = -Z(twopi-v,u), R(v,-u) = R(twopi-v,u)
!     TO DO EXTEND R,Z, etc. OVER ALL THETA (NOT JUST 0,PI)
!
         ivminus = mod(nzeta + 1 - iv,nzeta) + 1           !!(twopi-v)
         do iu = 1+ntheta2, ntheta1
            iu_r = ntheta1 + 2 - iu
            r1b(iu) = r1(ns,ivminus,iu_r,0) + r1(ns,ivminus,iu_r,1)
            z1b(iu) =-(z1(ns,ivminus,iu_r,0) + z1(ns,ivminus,iu_r,1))
            r12(iu) = r1(ns12,ivminus,iu_r,0) +
     1                r1(ns12,ivminus,iu_r,1)*sqrts(ns12)
            z12(iu) =-(z1(ns12,ivminus,iu_r,0) +
     1                z1(ns12,ivminus,iu_r,1)*sqrts(ns12))
            rub(iu) =-ru0(ns,ivminus,iu_r)
            zub(iu) = zu0(ns,ivminus,iu_r)
            ru12(iu)=-p5*(ru0(ns,ivminus,iu_r) + ru0(ns12,ivminus,iu_r))
            zu12(iu)= p5*(zu0(ns,ivminus,iu_r) + zu0(ns12,ivminus,iu_r))
         end do

         end if
!
!        Scan over r-z grid for interior point
!
         rmin = minval(r1b);  rmax = maxval(r1b)
         zmin = minval(z1b);  zmax = maxval(z1b)
         rcom(iv) = (rmax + rmin)/2; zcom(iv) = (zmax + zmin)/2

!
!        Estimate jacobian based on boundary and 1/2 surface
!
         ds = (ns - ns12)*hs
         do iu = 1, ntheta1
            rs(iu) = (r1b(iu) - r12(iu))/ds + r1(1,iv,1,0)
            zs(iu) = (z1b(iu) - z12(iu))/ds + z1(1,iv,1,0)
            tau0(iu) = ru12(iu)*zs(iu) - zu12(iu)*rs(iu)
         end do

         mintau = 0

         do nlim = 1, limpts
            zlim = zmin + ((zmax - zmin)*(nlim-1))/(limpts-1)
            if (.not.lasym .and. (iv.eq.1 .or. iv.eq.nzeta/2+1)) then
               zlim = 0
               if (nlim .gt. 1) exit
            end if
            rlim(:) = ( / ((rmin + ((klim-1)*(rmax - rmin))/(limpts-1)),
     1              klim = 1, limpts) / )

!
!           Find value of magnetic axis that maximizes the minimum jacobian value
!
            do klim = 1, limpts
               tau = tau0 - ru12(:)*zlim + zu12(:)*rlim(klim)
               mintemp = minval(tau*signgs)
               if (mintemp .gt. mintau) then
                  mintau = mintemp
                  rcom(iv) = rlim(klim)
                  zcom(iv) = zlim
!           If up-down symmetric and lasym=T, need this to pick z = 0
               else if (mintemp .eq. mintau) then
                  if (abs(zcom(iv)).gt.abs(zlim)) zcom(iv) = zlim
               end if
            end do
         end do

      end do planes

!
!     FOURIER TRANSFORM RCOM, ZCOM
!
      dzeta = two/nzeta
      do n = 0, ntor
         raxis(n,1) = dzeta*sum(cosnv(:,n)*rcom(:))/nscale(n)
         zaxis(n,1) =-dzeta*sum(sinnv(:,n)*zcom(:))/nscale(n)
         raxis(n,2) =-dzeta*sum(sinnv(:,n)*rcom(:))/nscale(n)
         zaxis(n,2) = dzeta*sum(cosnv(:,n)*zcom(:))/nscale(n)
         if (n.eq.0 .or. n.eq.nzeta/2) then
            raxis(n,1) = p5*raxis(n,1)
            zaxis(n,2) = p5*zaxis(n,2)
         end if
!        if (lscreen) print 100, n, raxis(n), zaxis(n)
      end do

!  100 format(' n = ',i4,' raxis = ',1pe10.3,' zaxis = ',1pe10.3)

      end subroutine guess_axis
      
      subroutine heading(extension, time_slice, iseq_count, lmac, 
     1     lscreen, lfirst)
      use vmec_main, only: rprec
      use vparams, only: nthreed, nmac
      use vmec_params, only: version_
      use date_and_computer
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: iseq_count
      real(rprec) :: time_slice
      character*(*) :: extension
      logical :: lmac, lscreen, lfirst
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      character*(50), parameter ::
     1   banner = ' THIS IS VMEC2000, A 3D EQUILIBRIUM CODE, VERSION '
      character*(*), parameter :: VersionID1 =
     1   ' Lambda: Full Radial Mesh. L-Force: hybrid full/half.',
     2   VersionID2 = ' Forces Are Conservative'
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: imon, nout
      character*(10) :: date0, time0, zone0
      character*(50) :: dateloc, Version
C----------------------------------------------- 
!
!     open output files
!
      call open_output_files (extension, iseq_count, lmac, lscreen, 
     1     lfirst)

c     FORTRAN-90 ROUTINE
      call date_and_time(date0,time0,zone0)
      read(date0(5:6),'(i2)')imon
      write(dateloc,100)months(imon),date0(7:8),date0(1:4),
     1  time0(1:2),time0(3:4),time0(5:6)
 100  format('DATE = ',a3,' ',a2,',',a4,' ',' TIME = ',2(a2,':'),a2)

      if (lscreen .and. lfirst) write (*,'(a,i4,a,1pe12.4/2a)')
     1  '  SEQ = ', iseq_count+1,
     2  ' TIME SLICE',time_slice,'  PROCESSING INPUT.', trim(extension)

      Version = trim(adjustl(version_)) 
      write(nthreed,'(a,1x,a,/,2a,//,2a,2x,a)') trim(banner), 
     1     trim(Version), trim(VersionID1), trim(VersionID2),
     2     ' COMPUTER: ', trim(computer), trim(dateloc)
      if (lscreen .and. lfirst) 
     1   write (*,'(1x,a,1x,a,/,1x,2a,//,1x,2a,2x,a)') trim(banner),
     2   trim(Version), trim(VersionID1), trim(VersionID2), 
     3   ' COMPUTER: ', trim(computer), trim(dateloc)

      do nout = nthreed, nthreed+1
        imon = nout
        if (imon .eq. nthreed+1) imon = nmac
        if (imon.eq.nmac .and. .not.lmac) cycle
        write (imon,3) trim(extension),iseq_count,time_slice
      enddo

 3    format(' SHOT ID.: ',a,2x,'SEQ. NO.:',i4,/,
     1       ' TIME SLICE = ',f5.0,' ms')

      end subroutine heading      

      subroutine interp(xnew, xold, scalxc, nsnew, nsold)
      use vmec_main, only: dp, rprec, mnsize
      use vmec_params, only: ntmax
      use vmec_persistent, only: ixm
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nsnew, nsold
      real(rprec), dimension(nsnew,mnsize,3*ntmax),
     1  intent(out) :: xnew
      real(rprec), dimension(nsnew,mnsize,3*ntmax),
     1  intent(in) :: scalxc
      real(rprec), dimension(nsold,mnsize,3*ntmax) :: xold
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: p5 = 0.5_dp, two = 2.0_dp
     1   ,zero = 0.0_dp, one = 1.0_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ntype, js, js1, js2
      real(rprec) :: hsold, sj, s1, xint
C-----------------------------------------------
 
      if (nsold .le. 0) return 
      hsold = one/(nsold - 1)
!
!       INTERPOLATE R,Z AND LAMBDA ON FULL GRID
!       (EXTRAPOLATE M=1 MODES,OVER SQRT(S), TO ORIGIN)
!       ON ENTRY, XOLD = X(COARSE MESH) * SCALXC(COARSE MESH)
!       ON EXIT,  XNEW = X(NEW MESH)   [ NOT SCALED BY 1/SQRTS ]
!
 
      do ntype = 1, 3*ntmax
 
         where (mod(ixm(:mnsize),2) .eq. 1) xold(1,:,ntype) =
     1       two*xold(2,:,ntype) - xold(3,:,ntype)

         do js = 1, nsnew
            sj = real(js - 1,rprec)/(nsnew - 1)
            js1 = 1 + ((js - 1)*(nsold - 1))/(nsnew - 1)
            js2 = min(js1 + 1,nsold)
            s1 = (js1 - 1)*hsold
            xint = (sj - s1)/hsold
            xint = min(one,xint)
            xint = max(zero,xint)
            xnew(js,:,ntype) = ((one - xint)*xold(js1,:,ntype)+xint*
     1         xold(js2,:,ntype))/scalxc(js,:,1)
         end do
 
 
!        Zero M=1 modes at origin
         where (mod(ixm(:mnsize),2) .eq. 1) xnew(1,:,ntype) = 0
 
      end do
 
      end subroutine interp
             

      subroutine load_xc_from_wout(rmn, zmn, lmn, ntor_in, mpol1_in,
     1    ns_in, reset_file, lreset)
      use read_wout_mod, only: rmnc, zmns, lmns, rmns, zmnc, lmnc, 
     1    xm, xn, ntor, ns, 
     2    nfp, mnmax, read_wout_file, read_wout_deallocate
      use vmec_params, only: rprec, mscale, nscale, ntmax
      use vmec_dim, only: mpol1
      use vparams, only: one, zero
      use vmec_input, only: lasym 
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: ns_in, mpol1_in, ntor_in
      real(rprec), dimension(ns_in,0:ntor_in,0:mpol1_in,ntmax), 
     1   intent(out) :: rmn, zmn, lmn
      character*(*) :: reset_file
      logical lreset
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: rbcc=1, rbss=2, rbcs=3, rbsc=4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ierr, mn, m, n, n1, zbcs, zbsc, zbcc, zbss
      real(rprec) :: t1, t2
C-----------------------------------------------


!
!     THIS ALLOWS SEQUENTIAL RUNNING OF VMEC FROM THE COMMAND LINE
!     i.e., WHEN VMEC INTERNAL ARRAYS ARE NOT KEPT IN MEMORY (compared to sequence file input)
!     THIS IS THE CASE WHEN VMEC IS CALLED FROM, SAY, THE OPTIMIZATION CODE
!
      call read_wout_file(reset_file, ierr)
      
      if (ierr .ne. 0) then
         print *,' Error opening/reading wout file in VMEC load_xc!'
         lreset = .true.
         return
      end if
      
      lreset = .false.   

      if (ns_in .ne. ns) then
         lreset = .true.
         print *, 'ns_in != ns in load_xc'
         return
      end if

      if (ntor_in  .ne. ntor ) stop 'ntor_in != ntor in load_xc'
      if (mpol1_in .ne. mpol1) stop 'mpol1_in != mpol1 in load_xc'
      if (nfp .eq. 0) stop 'nfp = 0 in load_xc'
      
      rmn = zero
      zmn = zero
      lmn = zero

      zbcs = rbcc
      zbsc = rbss
      zbcc = rbcs
      zbss = rbsc
      
      do mn = 1, mnmax
         m = nint(xm(mn))
         n = nint(xn(mn))/nfp
         n1 = abs(n)
         t1 = one/(mscale(m)*nscale(n1))
         t2 = t1
         if (n .lt. 0) t2 = -t2
         if (n .eq. 0) t2 = zero
         rmn(:ns, n1, m, rbcc) = rmn(:ns, n1, m, rbcc) + t1*rmnc(mn,:ns)
         zmn(:ns, n1, m, zbsc) = zmn(:ns, n1, m, zbsc) + t1*zmns(mn,:ns)
         lmn(:ns, n1, m, zbsc) = lmn(:ns, n1, m, zbsc) + t1*lmns(mn,:ns)
         rmn(:ns, n1, m, rbss) = rmn(:ns, n1, m, rbss) + t2*rmnc(mn,:ns)
         zmn(:ns, n1, m, zbcs) = zmn(:ns, n1, m, zbcs) - t2*zmns(mn,:ns)
         lmn(:ns, n1, m, zbcs) = lmn(:ns, n1, m, zbcs) - t2*lmns(mn,:ns)
         if (lasym) then
         rmn(:ns, n1, m, rbcs) = rmn(:ns, n1, m, rbcs) - t2*rmns(mn,:ns)
         zmn(:ns, n1, m, zbcc) = zmn(:ns, n1, m, zbcc) + t1*zmnc(mn,:ns)
         lmn(:ns, n1, m, zbcc) = lmn(:ns, n1, m, zbcc) + t1*lmnc(mn,:ns)
         rmn(:ns, n1, m, rbsc) = rmn(:ns, n1, m, rbsc) + t1*rmns(mn,:ns)
         zmn(:ns, n1, m, zbss) = zmn(:ns, n1, m, zbss) + t2*zmnc(mn,:ns)
         lmn(:ns, n1, m, zbss) = lmn(:ns, n1, m, zbss) + t2*lmnc(mn,:ns)
         end if
         if (m .eq. 0) then
            rmn(:ns, n1, m, rbss) = zero
            zmn(:ns, n1, m, zbsc) = zero
            lmn(:ns, n1, m, zbsc) = zero
            if (lasym) then
            rmn(:ns, n1, m, rbsc) = zero
            zmn(:ns, n1, m, zbss) = zero
            lmn(:ns, n1, m, zbss) = zero
            end if
         end if   
      end do   


      call read_wout_deallocate

      end subroutine load_xc_from_wout


      subroutine profil1d(xc, xcdot, lreset)
      use vmec_main
      use vmec_params, only: signgs
      use vsvd, torflux_edge => torflux
      use vspline
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(neqs2), intent(out) :: xc, xcdot
      logical, intent(in) :: lreset
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: c1p5 = 1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i
      real(rprec) :: Itor, si, pedge, tflux, tfluxd
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      real(rprec), external :: pcurr, pmass, piota, torflux, 
     1    torflux_deriv
C-----------------------------------------------
!
!                INDEX OF LOCAL VARIABLES
!
!        ai       array of coefficients in phi-series for iota (ncurr=0)
!        ac       array of coefficients in phi-series for the quantity d(jcurv)/ds = toroidal
!                 current density * Vprime, so jcurv = Itor(s) (ncurr=1)
!        am       array of coefficients in phi-series for mass (NWT/m**2)
!        iotas    rotational transform , on half radial mesh
!        jcurv    (-)toroidal current inside flux surface (vanishes like s)
!        mass     mass profile on half-grid
!        phiedge  value of real toroidal flux at plasma edge (s=1)
!        phips    same as phip , one-dimensional array
!        presf    pressure profile on full-grid, mass/phip**gamma
!        spres_ped value of s beyond which pressure profile is flat (pedestal)
 
!
!       COMPUTE PHIP, IOTA PROFILES ON FULL-GRID
!       COMPUTE MASS PROFILE ON HALF-GRID
!       BY READING INPUT COEFFICIENTS. PRESSURE CONVERTED TO
!       INTERNAL UNITS BY MULTIPLICATION BY mu0 = 4*pi*10**-7
!
      torflux_edge = signgs * phifac * phiedge / twopi
      r00 = rmn_bdy(0,0,1)

      phips(1) = 0
      jcurv(1) = 0

      do i = 2,ns
         si = hs*(i - c1p5)
         tflux = torflux(si)
         phips(i) = torflux_edge * torflux_deriv(si)
         iotas(i) = piota(tflux)
         jcurv(i) = pcurr(tflux)
      end do

      do i = 1,ns
         si = hs*(i-1)
         tflux = torflux(si)
         iotaf(i) = piota(tflux)
      enddo
!
!     SCALE CURRENT TO MATCH INPUT EDGE VALUE, CURTOR
!
      pedge = pcurr(one)
      Itor = 0
      if (abs(pedge) .gt. abs(epsilon(pedge)*curtor)) 
     1   Itor = signgs*currv/(twopi*pedge)
      jcurv(2:ns) = -signgs*Itor*jcurv(2:ns)

!
!     POSSIBLE PRESSURE PEDESTAL FOR S >= SPRES_PED
!
      spres_ped = abs(spres_ped)
      if (.not.lrecon) then
        do i = 2,ns
          si = hs*(i - c1p5)
          tflux = torflux(si)
          tfluxd = torflux_edge * torflux_deriv(si)
          if (si .gt. spres_ped) then 
             pedge = pmass(spres_ped)
          else
             pedge = pmass(tflux)
          end if   
          mass(i) = pedge*(abs(tfluxd)*r00)**gamma
        end do
 
      else
        iotas(:ns) = 0
        iotaf(:ns) = 0
        mass (:ns) = 0
        presf(:ns) = 0
      end if

      pres(:ns+1) = 0
      xcdot(:neqs2) = 0

      if (lreset) then
c        xc(:neqs1) = 0
        xc(:) = 0
        if (lrecon) iresidue = 0
      end if  
      if (lrecon) then
        if (iresidue .gt. 1) iresidue = 1        
!
!       COMPUTE INDEX ARRAY FOR FINDPHI ROUTINE
!
        do i = 1,ns
          indexr(i)    = ns + 1 - i                    !FINDPHI
          indexr(i+ns) = i + 1                         !FINDPHI
        enddo
        indexr(2*ns) = ns
      end if 

      end subroutine profil1d
      

      function pmass (xx)
      use kind_spec
      use vmec_input, only: am, bloat
      use vparams, only: dmu0
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer     :: im
      real(rprec) :: xx, pmass, x
C-----------------------------------------------
!     NOTE: On entry, am is in pascals. pmass internal units are mu0*pascals (B**2 units)
      x = min (abs(xx * bloat), 1._dp)
      pmass = 0
      do im = size(am), 1, -1
         pmass = x*pmass + am(im-1)
      end do
      pmass = dmu0*pmass

      end function pmass

      
      function piota (x)
      use kind_spec
      use vmec_input, only: ai
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer     :: iq
      real(rprec) :: x, piota
C-----------------------------------------------
      piota = 0
      do iq = size(ai), 1, -1
         piota = x*piota + ai(iq-1)
      end do

      end function piota

      
      function pcurr (xx)
      use kind_spec
      use vmec_input, only: ac, bloat, ac_form
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer     :: ic
      real(rprec) :: xx, pcurr, x
C-----------------------------------------------
      x = min (abs(xx * bloat), 1._dp)
      pcurr = 0
      select case (ac_form)
      case(1)
         pcurr = ac(0)*x + (2.*ac(1)*x**1.5)/3
         do ic = 2, size(ac)
            pcurr = pcurr + (ac(ic)*x**ic)/ic
         end do

      case default
         do ic = size(ac), 1, -1
            pcurr = x*pcurr + ac(ic-1)/ic
         end do
         pcurr = x*pcurr
      end select
   
      end function pcurr
      

      function torflux (x)
      use kind_spec
!     use vmec_input, only: af => aphi
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec) :: x, torflux
C-----------------------------------------------
!     TEMPORARILY DISABLED....
!     torflux = x*(af(1) + x*(af(2) + x*(af(3) + x*(af(4) + 
!    1     x*(af(5) + x*(af(6) + x*(af(7) + x*(af(8) + x*(af(9) + 
!    2     x*af(10))))))))))
      torflux = x

      end function torflux


      function torflux_deriv (x)
      use kind_spec
!     use vmec_input, only: af => aphi
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec) :: x, torflux_deriv
C-----------------------------------------------
!     TEMPORARILY DISABLED....
!     torflux_deriv = af(1) + x*(2*af(2) + x*(3*af(3) + x*(4*af(4) + 
!    1           x*(5*af(5) + x*(6*af(6) + x*(7*af(7) + x*(8*af(8) + 
!    2           x*(9*af(9) + x*10*af(10)))))))))
      torflux_deriv = 1

      end function torflux_deriv
      

      subroutine profil3d(rmn, zmn, lreset)
      use vmec_main
      use vmec_params, only: mscale, nscale, ntmax
      use vsvd, torflux_svd => torflux
      use vspline
      use realspace
      use xstuff
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,0:ntor,0:mpol1,ntmax),
     1    intent(inout) ::  rmn, zmn
      logical, intent(in) :: lreset
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: js, l, lk, lt, lz, ntype, m, n, mn, ntype2
      real(rprec), dimension(0:ntor,ntmax) :: rold, zold
      real(rprec) :: sm0, t1, facj, si
      integer :: jcount, jk, k
      real(rprec), external :: torflux
C-----------------------------------------------

!
!                INDEX OF LOCAL VARIABLES
!
!        phip     radial derivative of phi/(2*pi) on half-grid
!        shalf    sqrt(s) ,two-dimensional array on half-grid
!        sqrts    sqrt(s), two-dimensional array on full-grid
!        wint     two-dimensional array for normalizing angle integrations
!        ireflect two-dimensional array for computing 2pi-v angle

      do js = 1, ns
         phip(js:nrzt:ns) = phips(js)
      end do

      phip(nrzt+1) = zero

      faclam(:ns,0:ntor,0:mpol1,:ntmax) = zero


      do js = 1, ns
         si = hs*abs(js - 1.5_dp)
         si = torflux(si)
         shalf(js:nrzt:ns) = sqrt(si)
         si = hs*(js - 1)
         si = torflux(si)
         sqrts(js:nrzt:ns) = sqrt(si)
      end do
      sqrts(ns:nrzt:ns) = one     !!Avoid round-off

      lk = 0
      do lt = 1, ntheta3
         do lz = 1, nzeta
            lk = lk + 1
            if (lasym) then
               do js = 2, ns
                  wint(js+ns*(lk-1)) = one/(nzeta*ntheta1)
               end do
            else
               do js = 2, ns
                  wint(js+ns*(lk-1)) = cosmui(lt,0)/mscale(0)
               end do
            end if
         end do
      end do

      shalf(nrzt+1) = 1
      wint(1:nrzt:ns) = 0

!
!       COMPUTE ARRAY FOR REFLECTING v = -v (only needed for lasym)
!
      jcount = 0
      do k = 1, nzeta
         jk = nzeta + 2 - k
         if (k .eq. 1) jk = 1
         do js = 1,ns
           jcount = jcount+1
           ireflect(jcount) = js+ns*(jk-1)           !Index for -zeta[k]
         enddo
      end do

      do js = 2,ns
         sm(js) = shalf(js)/sqrts(js)
         sp(js) = shalf(js+1)/sqrts(js)
      enddo
      sm(1) = 0
      sp(0) = 0
      sp(1) = sm(2)


!
!     COMPUTE INITIAL R AND Z FOURIER COEFFICIENTS,
!     FROM SCALED BOUNDARY VALUES, AND SCALXC ARRAY
!     (1/SQRTS FACTOR FOR ODD M VALUES)
!
!
!     ALLOW FOR lreset == F RESTART PERTURBATION OF FIXED-BDY (LFREEB=F)
!

         do js = 1, ns
            si = hs*(js - 1)
            sm0 = one - torflux(si)
            do ntype = 1, ntmax
               do m = 0, mpol1
                  do n = 0, ntor
                     t1 = one/(mscale(m)*nscale(n))
                     mn = n + ntor1*m
                     l = js + ns*mn + (ntype - 1)*mns
                     if (mod(m,2) .eq. 0) then
                        scalxc(l) = one
                        scalxc(l+2*irzloff) = one   !Lambda on full mesh
                     else
                        scalxc(l) = one/max(sqrts(js),sqrts(2))
                        scalxc(l+2*irzloff)=one/max(sqrts(js),sqrts(2))
                     endif
                     if (.not.lreset .and. lfreeb) cycle
                     if (m .eq. 0) then
                        rmn(js,n,m,ntype) = rmn(js,n,m,ntype) +
     1                     (one - sm0)*(rmn_bdy(n,m,ntype)*t1
     2                    - rmn(ns,n,m,ntype))
                        zmn(js,n,m,ntype) = zmn(js,n,m,ntype) +
     1                    (one - sm0)*(zmn_bdy(n,m,ntype)*t1
     2                    - zmn(ns,n,m,ntype))
                        if (mod(ntype,2).eq.1 .and. lreset) then
                           if (js .eq. 1) then
                              rold(n,ntype) = rmn(1,n,0,ntype)
                              zold(n,ntype) = zmn(1,n,0,ntype)
                           endif
                           ntype2 = (ntype+1)/2
                           rmn(js,n,m,ntype) = rmn(js,n,m,ntype) +
     1                       sm0*(raxis(n,ntype2)*t1 - rold(n,ntype))
                           zmn(js,n,m,ntype) = zmn(js,n,m,ntype) +
     1                       sm0*(-zaxis(n,ntype2)*t1 - zold(n,ntype))
                        endif

                     else
                        facj = sqrts(js)**m        !!TURN OFF NEXT 3 LINES IF THIS ONE ACTIVATED
!                       if (mod(m,2) .eq. 0) then
!                           facj = sqrts(js)*sqrts(js)
!                       else if (mod(m,2) .eq. 1) then
!                          facj = sqrts(js)**min(m,3)
!                       end if
                        rmn(js,n,m,ntype) = rmn(js,n,m,ntype) + (rmn_bdy
     1                     (n,m,ntype)*t1 - rmn(ns,n,m,ntype))*facj
                        zmn(js,n,m,ntype) = zmn(js,n,m,ntype) + (zmn_bdy
     1                     (n,m,ntype)*t1 - zmn(ns,n,m,ntype))*facj
                     endif

                  end do
               end do
            end do
         end do

         scalxc(1+irzloff:2*irzloff) = scalxc(:irzloff)              !Z-components

!
!     STORE PHIFAC IN XC(NEQS1) ARRAY ELEMENT
!     STORE DEL-RMSE IN XC(NEQS2) ARRAY ELEMENT
!
      xc(neqs1) = phifac
      scalxc(neqs1) = 1
      scalxc(neqs2) = 1

      if (lrecon) then
!
!       COMPUTE SPLINE SEGMENTS FOR INTEGRALS OF FUNCTIONS IN SPLININT
!       THESE ARE FIXED DURING THE ITERATION SEQUENCE FOR FIXED SK,PK-NOTS,
!       BUT MAY CHANGE IF SKNOTS, PKNOTS CHANGE FOR DIFFERENT DATA FILES
!
        call setup_int (sknots, shalf(2), hstark, w_ia, w1_ia, u_ia,
     1       u1_ia, nk_ia, isnodes, ns1)
c-08-96 call setup_int(pknots,shalf(2),hthom,w_pa,w1_pa,u_pa,u1_pa,
c-08-96 >  nk_pa,ipnodes,ns1)
!
!       COMPUTE SPLINE SEGMENTS FOR INTEGRALS OF DERIVATIVES IN SPLININT
!
        call setup_int (sknots, sqrts, hstark, w_ib, w1_ib, u_ib,
     1       u1_ib, nk_ib, isnodes, ns)
        call setup_int (pknots, sqrts, hthom, w_pb, w1_pb, u_pb,
     1       u1_pb, nk_pb, ipnodes, ns)
      endif

      end subroutine profil3d

      
      subroutine evolve(time_step, ier_flag, liter_flag, lscreen)
      use vmec_main
      use vsvd
      use xstuff
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec) :: time_step
      integer ier_flag
      logical liter_flag, lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i
      real(rprec) :: fsq1, dtau, b1, fac
C-----------------------------------------------
 
!     COMPUTE MHD FORCES
 
      call funct3d (lscreen, ier_flag)
 
!     COMPUTE ABSOLUTE STOPPING CRITERION
 
      if (irst .ne. 1) then
         if (iter2 .eq. 1) ier_flag = 1
         return 
      else if (fsqr .le. ftolv .and. fsqz .le. ftolv .and. 
     1    fsql .le. ftolv) then
         liter_flag = .false.
         return 
      endif
 
!     COMPUTE DAMPING PARAMETER (DTAU) AND EVOLVE
!     R, Z, AND LAMBDA ARRAYS IN FOURIER SPACE
 
      if (iter2 .eq. iter1) otau(:ndamp) = cp15/time_step

      fsq1 = fsqr1 + fsqz1 + fsql1
      if (iter2 .gt. iter1) dtau = min(abs(log(fsq/fsq1)),cp15)
      fsq = fsq1
      if (iter2 .le. 1) return 
 
      otau(1:ndamp-1) = otau(2:ndamp)
 
      if (iter2 .gt. iter1) otau(ndamp) = dtau/time_step
      otav = sum(otau(:ndamp))/ndamp
      dtau = time_step*otav
      b1  = one - cp5*dtau
      fac = one/(one + cp5*dtau)
 
      do i = 1,neqs2
        xcdot(i) = fac*(xcdot(i)*b1 + time_step*gc(i))
        xc(i) = xc(i) + xcdot(i)*time_step
      end do  
 
      end subroutine evolve

      
      subroutine fileout(iseq, ier_flag, lscreen)
      use vmec_main
      use vac_persistent
      use realspace
      use vmec_params, only: mscale, nscale, signgs
      use vforces, czmn => czmn, lu => czmn, crmn => crmn, lv => crmn
      use vsvd
      use xstuff, ONLY: xc
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer iseq, ier_flag
      logical :: lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      logical, parameter :: lreset_xc = .false.
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: js, istat1=0
      real(rprec) :: twouton, twoutoff, teqon, teqoff      
      character*(92), dimension(0:10), save :: werror
      character*(35), dimension(0:6), save :: form
      character*(*), parameter :: 
     1    Warning = " global memory deallocation error "
      logical log_open
C-----------------------------------------------
      data werror/'EXECUTION TERMINATED NORMALLY', 
     1   'INITIAL JACOBIAN CHANGED SIGN (IMPROVE BOUNDARY GUESS)', 
     2   'MORE THAN 75 JACOBIAN ITERATIONS (DECREASE DELT)', 
     3   'VMEC INDATA ERROR: NCURR.ne.1 but BLOAT.ne.1.',
     4   'FORCE RESIDUALS EXCEED FTOL: INCREASING NUMBER ITERATIONS', 
     5   'ERROR READING INPUT FILE OR NAMELIST',
     6   'NEW AXIS GUESS STILL FAILED TO GIVE GOOD JACOBIAN',
     7   'PHIEDGE HAS WRONG SIGN IN VACUUM SUBROUTINE',
     8   'NS ARRAY MUST NOT BE ALL ZEROES',
     9   'ERROR READING MGRID FILE',
     A   'VACUUM-VMEC MISMATCH IN TOROIDAL CURRENT: INITIAL BOUNDARY MAY
     A BE ENCLOSING EXTERNAL CURRENT'/
      data form/'  TOTAL COMPUTATIONAL TIME :       ',
     1   '  TIME IN VACUUM LOOP :            ',
     2   '  TIME TO READ IN DATA:            ',
     3   '  TIME TO WRITE DATA TO WOUT:      ',
     4   '  TIME IN EQFORCE                  ',
     5   '  TIME (REMAINDER) IN FUNCT3D:     ',
     6   '  TIME IN PROFILE RECONSTRUCTION:  '/


      INTERFACE

      subroutine eqfor(bsubu, bsubv, tau, rz_array)
      use vmec_main
      use vmec_params, only: mscale, nscale, ntmax, signgs
      use realspace, br=>rcon, bz=>zcon
      use vforces, r12=>armn_o, bsupu => crmn_e, bsupv => czmn_e,
     1   gsqrt => azmn_o, bsq => bzmn_o, izeta => azmn_e,
     2   lu => czmn, lv => crmn, bphi => czmn_o
      use vacmod
      use vsvd
      use vspline
      use csplinx
      use vmec_io
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,nznt,0:1), intent(inout) :: 
     1  bsubu, bsubv
      real(rprec), dimension(nrzt), intent(out) :: tau
      real(rprec), dimension(ns,0:ntor,0:mpol1,3*ntmax),
     1  target, intent(in) :: rz_array

      end subroutine eqfor

      END INTERFACE

C-----------------------------------------------
!
!     COMPUTE REMAINING COVARIANT COMPONENT OF B (BSUBS),
!     CYLINDRICAL COMPONENTS OF B (BR, BPHI, BZ), AND
!     AVERAGE EQUILIBRIUM PROPERTIES AT END OF RUN
!
      iequi = 1

      if (ier_flag.eq.0 .or. ier_flag.eq.4) then
!
!     The sign of the jacobian MUST multiply phi to get the PHYSICALLY 
!     correct toroidal flux
!      
         phi(1) = zero
         do js = 2, ns
            phi(js) = phi(js-1) + phip(js)
         end do
         phi = (signgs*twopi*hs)*phi
      
         call funct3d (lscreen, ier_flag)
      end if

      call second0 (teqon)
      if (ier_flag .eq. 0) call eqfor (clmn, blmn, rcon(1,1), xc)
      call second0 (teqoff)
      timer(4) = timer(4) + teqoff - teqon

!
!     MUST Call WROUT To Write Error, if nothing else
!
      call second0 (twouton)

      call wrout (bzmn_o, azmn_o, clmn, blmn, crmn_o,
     1      czmn_e, crmn_e, azmn_e, ier_flag)
 
      call second0 (twoutoff)

      timer(3) = timer(3) + twoutoff - twouton
      timer(0) = timer(0) + timer(3) + timer(4)

      if (ier_flag .ne. 0) goto 1000

      timer(5) = timer(5) - timer(1) - timer(6)

      if (lscreen) print 10, ijacob
      write (nthreed, 10) ijacob
 10   format(/,'  NUMBER OF JACOBIAN RESETS = ',i4,/)

      if (lfreeb .and. lrecon) then
         if (lscreen) print 20, form(0),timer(0),form(2),timer(2),
     1            form(3),timer(3), form(4),timer(4),
     1            form(1),timer(1), form(6),timer(6),form(5),timer(5)
         write (nthreed, 20) form(0),timer(0),form(2),timer(2),
     1            form(3),timer(3),form(4),timer(4), form(1),timer(1),
     2            form(4),timer(4), form(6),timer(6),form(5),timer(5)
      else if (lfreeb) then
         if (lscreen) print 20, form(0),timer(0),form(2),timer(2),
     1               form(3),timer(3),form(4),timer(4), form(1),
     2               timer(1),form(5),timer(5)
         write (nthreed, 20) form(0),timer(0),form(2),timer(2), 
     1        form(3),timer(3),form(4),timer(4), form(1),timer(1),
     2        form(5),timer(5)
      else if (lrecon) then
         if (lscreen) print 20, form(0),timer(0),form(2),timer(2),
     1               form(3),timer(3),form(4),timer(4), 
     2               form(6),timer(6),form(5),timer(5)
         write (nthreed, 20) form(0),timer(0),form(2),timer(2),
     1        form(3),timer(3),form(4),timer(4), form(6),timer(6),
     2        form(5),timer(5)
      else
         if (lscreen) print 20,form(0),timer(0),form(2),timer(2),
     1              form(3),timer(3),form(4),timer(4), form(5),timer(5)
         write (nthreed, 20) form(0),timer(0),form(2),timer(2),
     1              form(3),timer(3),form(4),timer(4), form(5),timer(5)
      end if 
   20 format(a35,f12.2,' SECONDS')

         inquire(unit=nlog,opened=log_open)
         if (lrecon .and. log_open) then
!
!     WRITE SEQUENCE HISTORY FILE
!
         if (iseq .eq. 0) write (nlog, 100)
         write (nlog, 110) iseq + 1, iter2, total_chi_square_n,
     1   1.e-6_dp*ctor/dmu0, 1.e-3_dp*ppeak/dmu0, torflux, r00,
     2   timer(0), input_extension
      
      endif

  100 format(' SEQ ITERS  CHISQ/N',
     1   '  TORCUR  PRESMAX  PHIEDGE     R00 CPU-TIME  EXTENSION')
  110 format(i4,i6,f8.2,3f9.2,f8.2,f9.2,2x,a20)

 1000 continue

      if (lscreen) print 120, trim(werror(ier_flag)), input_extension
      write (nthreed, 120) trim(werror(ier_flag)), input_extension
  120 format(2x,a,/,2x,'FILE : ',a/)

!
!     DEALLOCATE GLOBAL MEMORY
!
      if (allocated(cosmu))
     1  deallocate(cosmu, sinmu, cosmum, sinmum, cosmui, cosmumi,
     2  sinmui, sinmumi, cosnv, sinnv, cosnvn, sinnvn, stat=istat1)
      if (istat1 .ne. 0) print *, Warning // "#1"

      if (allocated(xm)) deallocate (xm, xn, ixm, jmin3,
     1   mscale, nscale, stat=istat1)
      if (istat1 .ne. 0) print *, Warning // "#2"
     
      if (allocated(tanu))
     1  deallocate(tanu, tanv, sin2v, sinper, cosper, sinuv, cosuv, 
     2  sinu, cosu, sinv, cosv, sinui, cosui, cmns, csign, sinu1,
     3  cosu1, sinv1, cosv1, imirr, xmpot, xnpot, stat=istat1)
      if (istat1 .ne. 0) print *, Warning // "#3"

      call free_mem_funct3d
      call free_mem_ns (lreset_xc)
      call free_mem_nunv

!
!     CLOSE OPENED FILES
!
      if (ier_flag .ne. 4) call close_all_files
 
      end subroutine fileout


      subroutine fixaray
      use vmec_main
      use vmec_params, only: jmin2, mscale, nscale, signgs
      implicit none
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: 
     1    two=2, three=3, five=5, pexp = 4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, m, j, n, mn, mn1, nmin0, istat1, istat2, mnyq, nnyq
      real(rprec):: argi, arg, argj, dnorm
C-----------------------------------------------
!
!                INDEX OF LOCAL VARIABLES
!
!       mscale   array for norming theta-trig functions (internal use only)
!                so that the discrete SUM[cos(mu)*cos(m'u)] = .5 delta(m,m')
!       nscale   array for norming zeta -trig functions (internal use only)
 
 
!
!    COMPUTE TRIGONOMETRIC FUNCTION ARRAYS
!    NOTE: ARRAYS ALLOCATED HERE ARE GLOBAL AND ARE DEALLOCATED IN FILEOUT
!
      mnyq = max(0, ntheta1/2)           !!Nyquist (max) frequency: needed in jxbforce aliasing
      nnyq = max(0, nzeta/2)
      allocate(cosmu(ntheta2,0:mnyq), sinmu(ntheta2,0:mnyq),
     1           cosmum(ntheta2,0:mnyq), sinmum(ntheta2,0:mnyq), 
     2           cosmui(ntheta2,0:mnyq), cosmumi(ntheta2,0:mnyq),
     3           sinmui(ntheta2,0:mnyq), sinmumi(ntheta2,0:mnyq),
     4           cosnv(nzeta,0:nnyq),     sinnv(nzeta,0:nnyq), 
     5           cosnvn(nzeta,0:nnyq),    sinnvn(nzeta,0:nnyq),
     6           stat=istat1 )
      allocate( xm(mnmax),xn(mnmax), ixm(mnsize), jmin3(0:mnsize),
     1          mscale(0:mnyq), nscale(0:nnyq), stat=istat2)

      if (istat1.ne.0) stop 'allocation error in fixaray: istat1'
      if (istat2.ne.0) stop 'allocation error in fixaray: istat2'
      
      dnorm = one/(nzeta*(ntheta2 - 1))

      mscale(0) = osqrt2
      nscale(0) = osqrt2
      mscale(1:mnyq) = 1
      nscale(1:nnyq) = 1
      r0scale = mscale(0)*nscale(0)
 
      do i = 1, ntheta2
         argi = twopi*(i - 1)/ntheta1
         do m = 0, mnyq
            arg = argi*m
            cosmu(i,m) = cos(arg)*mscale(m)
            sinmu(i,m) = sin(arg)*mscale(m)
            cosmui(i,m) = dnorm*cosmu(i,m)
            sinmui(i,m) = dnorm*sinmu(i,m)
            if (i.eq.1 .or. i.eq.ntheta2) cosmui(i,m)=0.5_dp*cosmui(i,m)
            cosmum(i,m) = cosmu(i,m)*(m)
            sinmum(i,m) = -sinmu(i,m)*(m)
            cosmumi(i,m) = cosmui(i,m)*(m)
            sinmumi(i,m) = -sinmui(i,m)*(m)
         end do
      end do
 
      do j = 1, nzeta
         argj = twopi*(j - 1)/nzeta
         do n = 0, nnyq
            arg = argj*(n)
            cosnv(j,n) = cos(arg)*nscale(n)
            sinnv(j,n) = sin(arg)*nscale(n)
            cosnvn(j,n) = cosnv(j,n)*(n*nfp)
            sinnvn(j,n) = -sinnv(j,n)*(n*nfp)
         end do
      end do
 
!
!     R,Z,L / s**(m/2) ARE LINEAR NEAR ORIGIN
!
      mn = 0
      mn1 = 0
      jmin3 = 0                       ! mcz
      do m = 0, mpol1
         xmpq(m,1) = (m*(m - 1))
         xmpq(m,2) = (m**pexp)
         xmpq(m,3) = (m**(pexp+1))
         do n = 0, ntor
            jmin3(mn) = jmin2(m)
            mn = mn + 1
            ixm(mn) = m
         end do
         nmin0 = -ntor
         if (m .eq. 0) nmin0 = 0
         do n = nmin0, ntor
            mn1 = mn1 + 1
            xm(mn1) = (m)
            xn(mn1) = (n*nfp)
         end do
      end do

      if (mn1.ne.mnmax) stop 'mn1 != mnmax'
       
      faccon(0) = zero
      faccon(mpol1) = zero
c05-96    faccon(m) = -0.5_dp*signgs/((1+m)**4)
!!!         formfac = real(m**2,rprec)/(m+1)**2
c                                                !!!* formfac
      faccon(1:mpol1-1) = -0.25_dp*signgs/xmpq(2:mpol1,1)**2
 
      if (lrecon) call getgreen
      if (lfreeb) call precal                               !Fixed arrays for VACUUM

      end subroutine fixaray


      subroutine free_mem_funct3d
      use vmec_main
      use realspace
      use vforces
      use vacmod
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat1 = 0
C-----------------------------------------------
 
      if (allocated(armn))
     1   deallocate (armn, azmn, brmn, bzmn, crmn, czmn, blmn, clmn,
     1   r1, ru, rv, z1, zu, zv, rcon, zcon, ru0, zu0,
     2   rcon0, zcon0, guu, guv, gvv, stat=istat1)
      if (istat1 .ne. 0) stop 'deallocation error in funct3d'

      if (allocated(brv))
     1   deallocate (brv, bphiv, bzv, bpolvac, bsqvac, stat=istat1)
      if (istat1 .ne. 0) stop 'deallocation error in funct3d'

      end subroutine free_mem_funct3d

      
      subroutine free_mem_ns(lreset)
      use vmec_main
      use realspace
      use vforces
      use vsvd
      use xstuff
      use csplinx
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      logical, intent(in) :: lreset
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat1 = 0, istat2 = 0, istat3 = 0, istat4 = 0,
     >           istat5 = 0, istat6 = 0, istat7 = 0, istat8 = 0,
     >           istat9 = 0, istat10 = 0
C-----------------------------------------------
      if (allocated(phip)) 
     1  deallocate (phip, shalf, sqrts, wint, stat=istat3)
     
      if (allocated(ireflect))
     1  deallocate (ireflect,indexr,imid, stat=istat4)

      if (allocated(current))
     1  deallocate (current, rm2, vrm2, ovrm2, ochip, presph, presint,
     2      w_ia, w1_ia, u_ia, u1_ia, w_pa, w1_pa, u_pa, u1_pa,
     3      w_ib, w1_ib, u_ib, u1_ib, w_pb, w1_pb, u_pb, u1_pb,
     4      rmid, datamse, qmid, shear, presmid, alfa, curmid,
     5      curint, psimid, ageo, volpsi, isplinef, isplineh,
     6      psplinef, psplineh, phimid, pm, im, stat=istat5)


      if (allocated(pmb)) deallocate (pmb,imb,stat=istat6)

      if (allocated(ard))
     1  deallocate (ard,arm,brd,brm,crd,azd,azm,bzd,bzm, sm,sp,
     2        bmin, bmax,stat=istat7)

      if (allocated(iotaf))
     1  deallocate (iotaf,mass,phi,presf,jcuru,jcurv,jdotb,buco,bvco,
     2     bdotgradv,equif,specw,tcon,fpsi,psi,yellip,yinden,
     3     ytrian,yshift,ygeo,overr,faclam,iotas,phips,pres,vp,
     4     beta_vol, jperp2, jpar2, bdotb, clam, blam, dlam, 
     5     stat=istat8)

      if (allocated(rmidx))
     1  deallocate (rmidx,hmidx,wmidx,qmidx,tenmidx,ymidx,y2midx,
     2     stat=istat9)

      if (allocated(gc))
     1  deallocate (gc, xstore, xcdot, stat=istat10)
      if (allocated(xc) .and. lreset) deallocate (xc, scalxc)
     
      if (istat1.ne.0 .or. istat2.ne.0 .or. istat3.ne.0 .or.
     1      istat4.ne.0 .or. istat5.ne.0 .or. istat6.ne.0 .or.
     2      istat7.ne.0 .or. istat8.ne.0 .or. istat9.ne.0 .or.
     3      istat10.ne.0) then
          print *,' deallocation problem in free_mem_ns'
          print *,' istat1 = ',istat1,' istat2 = ',istat2
          print *,' istat3 = ',istat3,' istat4 = ',istat4
          print *,' istat5 = ',istat5,' istat6 = ',istat6
          print *,' istat7 = ',istat7,' istat8 = ',istat8
          print *,' istat9 = ',istat9,' istat10= ',istat10
       endif  

      end subroutine free_mem_ns

        
      subroutine free_mem_nunv
      use vmec_main
      use vacmod
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat1 = 0, istat2 = 0, istat3 = 0
C-----------------------------------------------

      if (allocated(bsubu0))
     1    deallocate (bsubu0, rbsq, dbsq, stat=istat1)
      if (allocated(rmn_bdy))
     1    deallocate (rmn_bdy, zmn_bdy, stat=istat2)

      if (allocated(amatsav))
     1    deallocate (amatsav, bvecsav, potvac, bsqsav, stat=istat3)
     

      if (istat1.ne.0 .or. istat2.ne.0 .or. istat3.ne.0) then
          print *,' deallocation problem in free_mem_nunv'
          print *,' istat1 = ',istat1,' istat2 = ',istat2
          print *,' istat3 = ',istat3
      endif  

      end subroutine free_mem_nunv

        
      subroutine free_mem_recon
      use vsvd
      use vspline
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat1 = 0, istat2 = 0, istat3 = 0
C-----------------------------------------------
      if (allocated(nk_ia))
     1  deallocate (nk_ia,nk_ib, nk_pa, nk_pb,
     2      indexs2,indexu2, indexs1, indexu1,
     3      isortr, isorts, stat=istat1)
      if (allocated(hthom))
     1  deallocate( hthom, ythom, y2thom, pknots,
     2  hstark,y2stark,ystark, sknots, stat=istat2 )
      if (allocated(sthom))
     1  deallocate( sthom, delse2, delso2, pcalc,
     2      delse1, delso1, starkcal, qmeas, qcalc,
     3      fpsical, stark_weight, rsort,rsort0, stat=istat3)
      if ((istat1 .ne. 0) .or. (istat2 .ne. 0)
     1                    .or. (istat3 .ne. 0) )then
          print *,' in free_mem_recon, istat1 = ',istat1
          print *,' istat2 = ',istat2,' istat3 = ',istat3
        end if

      end subroutine free_mem_recon


      subroutine free_persistent_mem()
      use mgrid_mod
      use vmec_main
      use vsvd
      use vspline
      use xstuff
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat0 = 0
      integer :: istat1 = 0, istat2 = 0, istat3 = 0, istat4 = 0
c-----------------------------------------------
      if (allocated(xc)) deallocate (xc, scalxc, stat=istat0)
      if (allocated(btemp)) deallocate(btemp, stat=istat1)
      if (allocated(bvac)) deallocate( bvac,stat=istat2 )
      if (allocated(xobser))
     1  deallocate( xobser, xobsqr, zobser, unpsiext, dsiext,
     2      psiext,plflux, iconnect, needflx, needbfld, plbfld, 
     3      nbcoils, rbcoil, zbcoil, abcoil, bcoil, rbcoilsqr, dbcoil,
     4      pfcspec,dsilabel, bloopnames, curlabel, b_chi,stat=istat3 )
      if (allocated(rlim))
     1  deallocate( rlim,zlim, reslim,seplim,stat=istat4 )
      if (lrecon) call free_mem_recon

      if (istat1.ne.0 .or. istat2.ne.0 .or.
     1    istat3.ne.0 .or. istat4.ne.0) then
          print *,'problem in free_persistent_mem'
          print *,' istat0 = '
          print *,' istat1 = ',istat1,' istat2 = ',istat2
          print *,' istat3 = ',istat3,' istat4 = ',istat4
      endif
     
      end subroutine free_persistent_mem
        

      subroutine funct3d (lscreen, ier_flag)
      use vmec_main
      use vacmod
      use realspace, z1 => z1, gcon => z1
      use vforces, czmn => czmn, lu => czmn, crmn => crmn, lv => crmn
      use vsvd
      use vspline
      use xstuff
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(inout) :: ier_flag
      logical, intent(in) :: lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: l0pi, l, lk, ivacskip
      real(rprec), dimension(mnmax) :: 
     1   rmnc, zmns, lmns, rmns, zmnc, lmnc
      real(rprec), dimension(nznt) :: rax, zax
      real(rprec), allocatable, dimension(:,:) :: extra1
     1   , extra2, extra3, extra4
      real(rprec) :: tfunon, tvacon, presf_ns, tvacoff,
     1   tfunoff, delr_mse
C-----------------------------------------------
 
      call second0 (tfunon)

      allocate (extra1(nrzt,0:1), stat=l)
      if (lasym) allocate (
     1   extra2(nrzt,0:1), extra3(nrzt,0:1), extra4(nrzt,0:1),
     2   stat=l)
      if (l .ne. 0) stop 'Allocation error in funct3d'
!
!     CONVERT ODD M TO 1/SQRT(S) INTERNAL REPRESENTATION
!
      call extrap(xc)
      gc(:neqs2) = xc(:neqs2)*scalxc(:neqs2)
 
!
!     RIGID BODY SHIFT OF RMNCC(JS.GT.1,0,0) BY DELR_MSE= R00-RAXMSE
!
      if (lrecon) then
         delr_mse = xc(neqs2)
         gc(1:ns) = gc(1:ns) + delr_mse
      endif

!
!     INVERSE FOURIER TRANSFORM TO S,THETA,ZETA SPACE
!     R, Z, AND LAMBDA ARRAYS IN FOURIER SPACE
!     FIRST, DO SYMMETRIC [ F(u,v) = F(-u,-v) ] PIECES
!     ON THE RANGE u = 0,pi  and v = 0,2*pi
!
      call totzsps (gc, r1, ru, rv, z1, zu, zv, lu, lv, rcon, zcon)

      if (lasym) then
!
!     NEXT, DO ANTI-SYMMETRIC PIECES
!
      call totzspa (gc, armn, brmn, extra3, azmn, bzmn, extra4, blmn, 
     1     clmn, extra1, extra2)

!     NOW SUM SYMMETRIC, ANTISYMMETRIC PIECES APPROPRIATELY
!     TO GET "R's, Z's, L's" ON FULL RANGE OF u (0 to 2*pi)
 
      call symrzl (r1, ru, rv, z1, zu, zv, lu, lv, rcon, zcon, armn, 
     1   brmn, extra3, azmn, bzmn, extra4, blmn, clmn, extra1, extra2)

      endif
 
      l0pi = ns*(1 + nzeta*(ntheta2 - 1))        !u = pi, v = 0, js = ns
      router = r1(ns,0) + r1(ns,1)
      rinner = r1(l0pi,0) + r1(l0pi,1)
      r00 = r1(1,0)
      z00 = z1(1,0)

!
!     COMPUTE CONSTRAINT RCON, ZCON
!
      do l = 1,nrzt
         rcon(l,0) = rcon(l,0) + rcon(l,1)*sqrts(l)
         zcon(l,0) = zcon(l,0) + zcon(l,1)*sqrts(l)
         ru0(l) = ru(l,0) + ru(l,1)*sqrts(l)
         zu0(l) = zu(l,0) + zu(l,1)*sqrts(l)
      end do 

! 
!     COMPUTE S AND THETA DERIVATIVE OF R AND Z AND JACOBIAN ON HALF-GRID
!
      call jacobian
      if (irst.eq.2 .and. iequi.eq.0) goto 100 

!
!     COMPUTE RCON0, ZCON0 FOR FIXED BOUNDARY BY SCALING EDGE VALUES
!     SCALE BY POWER OF SQRTS SO THAT RESTART FOR FIXED BOUNDARY DOES NOT
!     HAVE A DISCONTINUITY DUE TO NEW RCON0....
!
      if (iter2.eq.1 .and. .not.lfreeb) then
         do l = 1, ns
            rcon0(l:nrzt:ns) = rcon(ns:nrzt:ns,0)*sqrts(l:nrzt:ns)**2
            zcon0(l:nrzt:ns) = zcon(ns:nrzt:ns,0)*sqrts(l:nrzt:ns)**2
         end do
!        rcon0(:nrzt) = rcon(:nrzt,0)
!        zcon0(:nrzt) = zcon(:nrzt,0)
      endif

! 
!     COMPUTE COVARIANT COMPONENTS OF B, MAGNETIC AND KINETIC
!     PRESSURE, AND METRIC ELEMENTS ON HALF-GRID
!
      call bcovar (lu, lv)

!     COMPUTE VACUUM MAGNETIC PRESSURE AT PLASMA EDGE
!     NOTE: FOR FREE BOUNDARY RUNS, THE VALUE OF RBTOR=R*BTOR
!     AT THE PLASMA EDGE SHOULD BE ADJUSTED TO APPROXIMATELY
!     EQUAL THE VACUUM VALUE. THIS CAN BE DONE BY CHANGING
!     EITHER PHIEDGE OR THE INITIAL CROSS SECTION ACCORDING
!     TO THE SCALING LAW  R*BTOR .EQ. PHIEDGE/(R1 * Z1).
 
      if (lfreeb .and. iter2.gt.1) then
         if (fsqr + fsqz .le. 1.e2_dp) ivac = ivac + 1
         if (ivac .ge. 0) then
            call second0 (tvacon)
            ivacskip = mod(iter2 - iter1,nvacskip)
            if (ivac .le. 2) ivacskip = 0
            do lk = 1, nznt
               rax(lk) = r1(1 + ns*(lk - 1),0)
               zax(lk) = z1(1 + ns*(lk - 1),0)
            end do
            call convert (rmnc, zmns, lmns, rmns, zmnc, lmnc, gc, ns)
            call vacuum (rmnc, rmns, zmns, zmnc, xm, xn, rax, zax, 
     1         ctor, rbtor, wint, ns, ivacskip, ivac, mnmax, 
     2         ier_flag, lscreen)
            if (ier_flag .ne. 0) return
            lk = 0
!
!           IN CASE PRESSURE IS NOT ZERO AT EXTRAPOLATED EDGE...
!           UNCOMMENT ALL "RPRES" COMMENTS HERE AND IN BCOVAR, FORCES ROUTINES
!           IF NON-VARIATIONAL FORCES ARE DESIRED 
!
            presf_ns = 1.5_dp*pres(ns) - 0.5_dp*pres(ns1)
!RPRES      if (iequi .ne. 1) presf_ns = 0

            do l = ns, nrzt, ns
               lk = lk + 1
               bsqsav(lk,3) = 1.5_dp*bzmn_o(l) - 0.5_dp*bzmn_o(l-1)
               rbsq(lk) = (bsqvac(lk) + presf_ns)*ohs*(r1(l,0)+r1(l,1))
               dbsq(lk) = abs(bsqvac(lk) + presf_ns - bsqsav(lk,3))
            end do
            if (ivac .eq. 1) then
               bsqsav(:nznt,1) = bzmn_o(ns:nrzt:ns)
               bsqsav(:nznt,2) = bsqvac(:nznt)
            endif
            call second0 (tvacoff)
            timer(1) = timer(1) + (tvacoff - tvacon)
         endif
      endif

!
!     COMPUTE CONSTRAINT FORCE (GCON => Z1)
!
      extra1(:nrzt,1) = (rcon(:nrzt,0) - rcon0(:nrzt))*ru0(:nrzt)
     1    + (zcon(:nrzt,0) - zcon0(:nrzt))*zu0(:nrzt)
      call alias (gcon, extra1(:,0), extra1(:,1), gc, gc(1+mns), 
     1   gc(1+2*mns), gc(1+3*mns))


      if (iequi .eq. 1) then
         if (lrecon) xc(:ns) = xc(:ns) + delr_mse
         goto 100
      end if

!
!     COMPUTE MHD FORCES ON INTEGER-MESH
!
      call forces

!
!     SYMMETRIZE FORCES (in u-v space) 
!
      if (lasym) call symforce (armn, brmn, crmn, azmn, bzmn,
     1     czmn, blmn, clmn, rcon, zcon, r1, ru, rv, z1, zu, zv, 
     2     extra3, extra4, extra1, extra2)
 
!
!     FOURIER-TRANSFORM MHD FORCES TO (M,N)-SPACE
!
      call tomnsps (gc, armn, brmn, crmn, azmn, bzmn, czmn, blmn, clmn, 
     1     rcon, zcon)

      if (lasym) call tomnspa (gc, r1, ru, rv, z1, zu, zv,
     1     extra3, extra4, extra1, extra2)

!
!     COMPUTE FORCE RESIDUALS
!
      gc = gc * scalxc
      call residue (gc, gc(1+irzloff), gc(1+2*irzloff), faclam)
 
      gc(neqs1) = gphifac
      if (iopt_raxis .gt. 0) gc(neqs2) = grmse

 100  continue

      deallocate (extra1)
      if (lasym) deallocate (extra2, extra3, extra4)
      call second0 (tfunoff)
      timer(5) = timer(5) + (tfunoff - tfunon)
 
      end subroutine funct3d

              
      subroutine extrap(rzl_array)
      use vmec_main
      use vmec_params, only: jmin2, jlam, ntmax
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,0:ntor,0:mpol1,3*ntmax),
     1   intent(inout) :: rzl_array
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: rmnmin, rmnmax
      integer :: zmnmax, lmnmin, lmnmax, m
      real(rprec) :: power
C-----------------------------------------------
      rmnmin = 1
      if (lasym) then
         rmnmax = 4
      else
         rmnmax = 2
      end if

      zmnmax = rmnmax + ntmax
      lmnmin = rmnmin + 2*ntmax
      lmnmax = rmnmax + 2*ntmax

!
!     EXTRAPOLATION AT JS=2 (PUT NONSYMMETRIC HERE, TOO)
!
      do m = 0, mpol1
         if (jmin2(m) .lt. 3) cycle
         power = osqrt2**m
         rzl_array(2,:,m,rmnmin:zmnmax) = 
     1   rzl_array(3,:,m,rmnmin:zmnmax)*power
      end do
      do m = 0, mpol1
         if (jlam(m) .lt. 3) cycle
         power = osqrt2**m
         rzl_array(2,:,m,lmnmin:lmnmax) = 
     1   rzl_array(3,:,m,lmnmin:lmnmax)*power
      end do

      end subroutine extrap




      subroutine alias(gcons, gcona, ztemp, gcs, gsc, gcc, gss)
      use vmec_main
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), parameter :: p5 = 0.5_dp
      real(rprec), dimension(ns*nzeta,ntheta3) ::
     1   gcons, gcona, ztemp
      real(rprec), dimension(ns,0:ntor,0:mpol1) ::
     1   gcs, gsc, gcc, gss
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: m, i, ir, jk, jka, n, k, js, l
      real(rprec), dimension(:,:), allocatable :: work
C-----------------------------------------------

      allocate (work(ns*nzeta,4))

      gcons = 0
      gcona = 0

      gcs(:,:ntor,:mpol1) = 0;  gsc(:,:ntor,:mpol1) = 0
      gcc(:,:ntor,:mpol1) = 0;  gss(:,:ntor,:mpol1) = 0

      do m = 1, mpol1 - 1
         work = 0
         do i = 1, ntheta2
            ir = ntheta1 + 2 - i
            if (i .eq. 1) ir = 1
            do jk = 1, ns*nzeta
               work(jk,1) = work(jk,1) + ztemp(jk,i)*cosmui(i,m)
               work(jk,2) = work(jk,2) + ztemp(jk,i)*sinmui(i,m)
            end do
            if (lasym) then
            do jk = 1, ns*nzeta
               jka = ireflect(jk)
               work(jk,3) = work(jk,3) + ztemp(jka,ir)*cosmui(i,m)
               work(jk,4) = work(jk,4) + ztemp(jka,ir)*sinmui(i,m)
            end do
            end if
         end do

         do n = 0, ntor
            do k = 1, nzeta
               l = ns*(k-1)
               if (.not.lasym) then
               do js = 2,ns
                  gcs(js,n,m) = gcs(js,n,m) + tcon(js)*work(js+l,1)*
     1               sinnv(k,n)
                  gsc(js,n,m) = gsc(js,n,m) + tcon(js)*work(js+l,2)*
     1               cosnv(k,n)
               end do
               else
               do js = 2,ns
                  gcs(js,n,m) = gcs(js,n,m) + p5*tcon(js)*sinnv(k,n)*
     1               (work(js+l,1)-work(js+l,3))
                  gsc(js,n,m) = gsc(js,n,m) + p5*tcon(js)*cosnv(k,n)*
     1               (work(js+l,2)-work(js+l,4))
                  gss(js,n,m) = gss(js,n,m) + p5*tcon(js)*sinnv(k,n)*
     1               (work(js+l,2)+work(js+l,4))
                  gcc(js,n,m) = gcc(js,n,m) + p5*tcon(js)*cosnv(k,n)*
     1               (work(js+l,1)+work(js+l,3))
               end do
               end if
            end do
         end do
!
!        INVERSE FOURIER TRANSFORM DE-ALIASED GCON
!
         work = 0

         do n = 0, ntor
            do k = 1, nzeta
               l = ns*(k-1)
               do js = 2, ns
                  work(js+l,3) = work(js+l,3) + gcs(js,n,m)*sinnv(k,n)
                  work(js+l,4) = work(js+l,4) + gsc(js,n,m)*cosnv(k,n)
               end do
               if (lasym) then
               do js = 2, ns
                  work(js+l,1) = work(js+l,1) + gcc(js,n,m)*cosnv(k,n)
                  work(js+l,2) = work(js+l,2) + gss(js,n,m)*sinnv(k,n)
               end do
               end if
            end do
         end do
         do i = 1, ntheta2
            do jk = 1, ns*nzeta
               gcons(jk,i) = gcons(jk,i) + (work(jk,3)*cosmu(i,m)
     1                     + work(jk,4)*sinmu(i,m))*faccon(m)
            end do
            if (lasym) then
            do jk = 1, ns*nzeta
               gcona(jk,i) = gcona(jk,i) + (work(jk,1)*cosmu(i,m)
     1                     + work(jk,2)*sinmu(i,m))*faccon(m)
            end do
            end if
         end do
      end do

      if (lasym) then
!
!     EXTEND GCON INTO THETA = PI,2*PI DOMAIN
!
      do i = 1 + ntheta2, ntheta1
         ir = ntheta1 + 2 - i
         do jk = 1, ns*nzeta
            jka = ireflect(jk)
            gcons(jk,i) = -gcons(jka,ir) + gcona(jka,ir)
         end do
      end do
!
!     ADD SYMMETRIC, ANTI-SYMMETRIC PIECES IN THETA = 0,PI DOMAIN
!
      gcons(:,:ntheta2) = gcons(:,:ntheta2) + gcona(:,:ntheta2)

      end if

      deallocate (work)

      end subroutine alias

      
      subroutine bss(r12, rs, zs, ru12, zu12, bsubs, bsupu, bsupv,
     1               br, bphi, bz)
      use vmec_main
      use realspace
      use vsvd, only: torflux_edge => torflux
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(nrzt), intent(in) :: r12, rs, zs, 
     1     ru12, zu12, bsupu, bsupv 
      real(rprec), dimension(nrzt), intent(out) :: 
     1     br, bphi, bz, bsubs
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: p5 = 0.5_dp, p25 = p5*p5
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: l
      real(rprec) :: rv12, zv12, gsu, gsv, dphids, t1,
     1    rs12, zs12      
C-----------------------------------------------
      dphids = p25/torflux_edge
 
      do l = 2, nrzt
         t1 = dphids*phip(l)
         rv12 = p5*(rv(l,0)+rv(l-1,0) + shalf(l)*(rv(l,1) + rv(l-1,1)))
         zv12 = p5*(zv(l,0)+zv(l-1,0) + shalf(l)*(zv(l,1) + zv(l-1,1)))
         rs12 = rs(l) + t1*(r1(l,1) + r1(l-1,1))/shalf(l)
         zs12 = zs(l) + t1*(z1(l,1) + z1(l-1,1))/shalf(l)
         gsu  = rs12*ru12(l) + zs12*zu12(l)
         gsv  = rs12*rv12    + zs12*zv12
         br(l)    = bsupu(l)*ru12(l) + bsupv(l)*rv12
         bphi(l)  = bsupv(l)*r12(l)
         bz(l)    = bsupu(l)*zu12(l) + bsupv(l)*zv12
         bsubs(l) = bsupu(l)*gsu + bsupv(l)*gsv
      end do
      end subroutine bss
      

      subroutine eqfor(bsubu, bsubv, tau, rz_array)
      use vmec_main
      use vmec_params, only: mscale, nscale, ntmax, signgs
      use realspace, br=>rcon, bz=>zcon
      use vforces, r12=>armn_o, bsupu => crmn_e, bsupv => czmn_e,
     1   gsqrt => azmn_o, bsq => bzmn_o, izeta => azmn_e,
     2   lu => czmn, lv => crmn, bphi => czmn_o
      use vacmod
      use vsvd
      use vspline
      use csplinx
      use vmec_io
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,nznt,0:1), intent(inout) :: 
     1  bsubu, bsubv
      real(rprec), dimension(nrzt), intent(out) :: tau
      real(rprec), dimension(ns,0:ntor,0:mpol1,3*ntmax),
     1  target, intent(in) :: rz_array
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: p5=0.5_dp, c1p5=1.5_dp, two=2.0_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, icount, itheta, js, js1, l, loff,
     1   lpi, lt, n, n1, nchicur, nchiiota0, ncol, noff,
     2   nout, nsort, iv, iu, lk
      real(rprec), dimension(:,:), pointer ::
     1   rmags, zmags, rmaga, zmaga
      real(rprec), dimension(:,:,:), pointer :: rmncc,zmnsc
      real(rprec), dimension(ns) :: phipf, phi1, chi1, t3,
     1   bvcof, t11u, t21u, jPS2
      real(rprec) :: modb(nznt)
      real(rprec), dimension(:), allocatable ::
     1   bsup1u, dbts1u, dint1u, t12u, guu_1u, guus1u, r3v,
     2   redg1u, rbps1u
      real(rprec) :: aminr1, aminr2, aminr2in, anorm,
     1   area, aspectratio, betai, betpol, betstr,
     2   bminz2, bminz2in, bsq1, btor, iotamax, musubi,
     3   btorsq, bzcalc, bzin, chisq, chiwgt, circum, cur0,
     4   delphid_exact, delta1, delta2, delta3, denwgt, lambda,
     5   dlogidsi, dmusubi_meas, er, es, fac, facnorm, factor, fgeo,
     6   fmax, fmin, flao, fpsi0, pavg, pitchc, pitchm,
     7   pprime, qedge, qmin1, qmin2, qmin3, qzero,
     8   raxis0, rcalc, rcen, rcenin, rgeo, rjs,
     9   rjs1, rlao, rqmin1, rqmin2, rshaf, rshaf1, rshaf2, s11, s12,
     A   s13, s2, s3, sarea, sigr0, sigr1, sigz1, smaleli,
     B   splintx, splints, sqmin, sumbpo, sumbto, sumbtr, sump,
     C   sump2, sump20, t1, jpar_perp=zero, jparPs_perp=zero,
     D   tol, vnorm, volf, vprime, wght0, xmax,
     E   xmida, xmidb, xmin, rzmax, rzmin, zxmax, zxmin, zaxis0,
     F   zmax, zmin, yr1u, yz1u, waist(2), height(2)
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
       external splintx,splints
C-----------------------------------------------

!
!     POINTER ASSOCIATIONS
!
      rmags => rz_array(:,:,0,1)
      zmags => rz_array(:,:,0,1+ntmax)
      rmaga => rz_array(:,:,0,3)
      zmaga => rz_array(:,:,0,3+ntmax)
      zmnsc => rz_array(:,:,:,2+ntmax)
      rmncc => rz_array(:,:,:,1)

!
!     NOTE: JXBFORCE ROUTINE MUST BE CALLED FIRST TO COMPUTE IZETA, JDOTB
!           ON OUTPUT, J, IZETA, JDOTB ARE IN MKS UNITS (HAVE 1/MU0)
!
      call bss (r12, bzmn, brmn, azmn, armn, lv(1+nrzt), bsupu, bsupv,
     1          br, bphi, bz)

!
!     STORE EDGE VALUES OF B-FIELD
!
      if (lfreeb .or. ledge_dump) then
         allocate (bredge(2*nznt), bpedge(2*nznt), bzedge(2*nznt))
         do iv = 1,nzeta
            do iu = 1,ntheta3
               lk = iv + nzeta*(iu-1)
               n1 = ns*lk
               bredge(lk) = 1.5_dp*br(n1,0)   - p5*br(n1-1,0)
               bpedge(lk) = 1.5_dp*bphi(n1)   - p5*bphi(n1-1)
               bzedge(lk) = 1.5_dp*bz(n1,0)   - p5*bz(n1-1,0)
            end do
         end do
      end if

      call jxbforce (bsupu, bsupv, bsubu, bsubv, lv(1+nrzt),
     1               br, bz, rcon0, zcon0, gsqrt, izeta, bsq)

!
!     HALF-MESH VOLUME-AVERAGED BETA
!
      tau(1) = zero
      tau(2:nrzt) = signgs*wint(2:nrzt)*gsqrt(2:nrzt)
      do i = 2, ns
         s2 = sum(bsq(i:nrzt:ns)*tau(i:nrzt:ns))/vp(i) - pres(i)
         beta_vol(i) = pres(i)/s2
      end do

      betaxis = c1p5*beta_vol(2) - p5*beta_vol(3)


      write (nthreed, 5)
    5 format(/,' NOTE: <RADIAL FORCE> = d(Ipol)/dPHI',
     1         ' - IOTA*d(Itor)/dPHI - dp/dPHI * dV/dPHI',/,
     1         ' (NORMED TO SUM OF INDIVIDUAL TERMS)',//,
     1         '     S      <RADIAL    TOROIDAL     IOTA     ',
     1         ' <JPOL>     <JTOR>     d(VOL)/',
     2         '   d(PRES)/    <M>     PRESF     <BSUBV> ',
     3         '     <J.B>     <B.B>',/,
     4         '             FORCE>      FLUX                ',
     5         ' (X mu0)    (X mu0)    d(PHI) ',
     6         '    d(PHI)                             ',/,137('-'),/)

      phipf(1) = twopi*signgs*(c1p5*phip(2) - p5*phip(3))
      presf(1) = c1p5*pres(2) - p5*pres(3)
      do i = 2,ns1
         presf(i) = p5*(pres(i) + pres(i+1))
         phipf(i) = p5*twopi*signgs*(phip(i) + phip(i+1))
      end do
      presf(ns) = c1p5*pres(ns)- p5*pres(ns-1)
      phipf(ns) = twopi*signgs*(c1p5*phip(ns) - p5*phip(ns1))

      phi1(1) = zero
      chi1(1) = zero
      do i = 2, ns
         buco(i) = sum(bsubu(i,:,0)*wint(i:nrzt:ns))
         bvco(i) = sum(bsubv(i,:,0)*wint(i:nrzt:ns))
         phi1(i) = phi1(i-1) + hs*phip(i)
         chi1(i) = chi1(i-1) + hs*(phip(i)*iotas(i))
      end do

      bvcof(1) = c1p5*bvco(2) - p5*bvco(3)

!
!     NOTE:  jcuru, jcurv on FULL radial mesh now.
!            They are local (surface-averaged) current densities (NOT integrated in s)
!
      do i = 2,ns1
         jcurv(i) = (signgs*ohs)*(buco(i+1) - buco(i))
         jcuru(i) =-(signgs*ohs)*(bvco(i+1) - bvco(i))
         t11u(i)  = jcurv(i)*iotaf(i)
         t21u(i)  = ohs*(pres(i+1) - pres(i))
         t3(i)    = p5*(vp(i+1)/phip(i+1) + vp(i)/phip(i))
         equif(i) = (jcuru(i) - t11u(i) - (t21u(i)*t3(i)))/(abs(t11u(i))
     1            +  abs(jcuru(i))+abs((t21u(i)*t3(i))))
         bvcof(i) = p5*(bvco(i) + bvco(i+1))
      end do

      bvcof(ns) = c1p5*bvco(ns) - p5*bvco(ns1)

      equif(1) = two*equif(2) - equif(3)
      jcuru(1) = two*jcuru(2) - jcuru(3)
      jcurv(1) = two*jcurv(2) - jcurv(3)
      t21u(1)  = two*t21u(2)  - t21u(3)
      t21u(ns) = two*t21u(ns1) - t21u(ns1-1)
      t3(1)  = two*t3(2) - t3(3)
      t3(ns) = two*t3(ns1) - t3(ns1-1)
      equif(ns) = two*equif(ns1) - equif(ns1-1)
      jcuru(ns) = two*jcuru(ns1) - jcuru(ns1-1)
      jcurv(ns) = two*jcurv(ns1) - jcurv(ns1-1)
      fac = twopi*signgs
      do js = 1, ns
         es = (js - 1)*hs
         cur0 = signgs*fac*t3(js)*phipf(js)
         write (nthreed, 30) es, equif(js), fac*phi1(js),
     1     iotaf(js), fac*jcuru(js)/cur0, fac*jcurv(js)/cur0,
     2     fac*t3(js), t21u(js)/phipf(js), specw(js), presf(js)/dmu0,
     3     bvcof(js), jdotb(js), bdotb(js)
      end do
 30   format(1p2e10.2,1p6e11.3,0pf7.3,1p4e11.3)

!
!     Calculate poloidal cross section area & toroidal flux
!
      anorm = twopi*hs
      vnorm = twopi*anorm
      area = sum(tau(2:nrzt)/r12(2:nrzt))
      do i = 2, ns
         overr(i) = sum(tau(i:(nznt-1)*ns+i:ns)
     1      /r12(i:(nznt-1)*ns+i:ns)) / vp(i)
      end do
      area = area*anorm
      volf = vnorm*sum(vp(2:ns))
      torflux = anorm * sum(bsupv(2:nrzt)*tau(2:nrzt))
      write (nthreed, 3) torflux, phifac
    3 format(/,' TOROIDAL FLUX    = ',f10.2,7x,'PHIFAC = ',f10.2)

!
!
!     OUTPUT BETAS, INDUCTANCES, SAFETY FACTORS, ETC.
!     (EXTRACTED FROM FQ-CODE, 9-10-92)
!
!     b poloidals (cylindrical estimates)
!
      rcen = p5*(router + rinner)               !geometric center
      n = 0
      n1 = n + 1
      rcenin = dot_product(rmncc(ns,n1,:mpol1+1:2),
     1                     mscale(:mpol1:2)*nscale(n))

      l = (mpol1+1)/2
      allocate (t12u(l))
      t12u(:l) = mscale(1:mpol1:2)*nscale(n)
      aminr2in = dot_product(rmncc(ns,n1,2:mpol1+1:2),t12u(:l))
      bminz2in = dot_product(zmnsc(ns,n1,2:mpol1+1:2),t12u(:l))
      bminz2 = dot_product(zmnsc(ns,n1,2:mpol1+1:2),t12u(:l))
      deallocate (t12u)
      aminr1 = sqrt(two*volf/(twopi*twopi*r00))  !vol av minor radius
      aminr2 = p5*(router - rinner)             !geometric plasma radius
!
!       cylindrical estimates for beta poloidal
      sump = vnorm*sum(vp(2:ns)*pres(2:ns))
      pavg = sump/volf
      ppeak = presf(1)
      factor = two*pavg
!
!       delphid_exact = Integral[ (Bvac - B) * dSphi ]
!
      delphid_exact = zero                         !Eq. 20 in Shafranov
      musubi = zero

      rshaf1 = zero
      rshaf2 = zero
      allocate (bsup1u(nznt), dbts1u(nznt), dint1u(nznt))
      do js = 2, ns
         bsup1u(:nznt) = bsubvvac/(r12(js:nrzt:ns)*
     1      r12(js:nrzt:ns))
         delphid_exact = delphid_exact + sum((bsup1u(:nznt) -
     1      bsupv(js:nrzt:ns))*tau(js:nrzt:ns))
         dbts1u(:nznt) = bsup1u(:nznt)*bsubvvac -
     1      bsupv(js:nrzt:ns)*bsubv(js,:nznt,0)
         musubi = musubi + sum(dbts1u(:nznt)*tau(js:nrzt:ns))
         dint1u(:nznt) = (bsupu(js:nrzt:ns)*bsubu(js,:nznt,0)
     1      + dbts1u(:nznt) + two*pres(js))*tau(js:nrzt:ns)
         rshaf1 = rshaf1 + sum(dint1u(:nznt))
         rshaf2 = rshaf2 + sum(dint1u(:nznt)/r12(js:nrzt:ns))
      end do
      deallocate (bsup1u, dbts1u, dint1u)

      delphid_exact = anorm*delphid_exact
      rshaf = rshaf1/rshaf2

      if (lrecon) then
         if (apres .ne. aminr1) then
            write (*, 50) (apres/aminr1)**2
            write (nthreed, 50) (apres/aminr1)**2
         endif
   50 format(/'  Multiply final PHIEDGE by ',f7.3,
     1   ' to make apres = aminr1')
!
!       Output some of measured data
!

         raxis0 = sum(raxis(0:ntor,1))
         zaxis0 = sum(zaxis(0:ntor,1))
         fpsi0 = c1p5*bvco(2) - p5*bvco(3)
         b0 = fpsi0/r00
         cur0 = signgs*iotaf(1)*fpsi0/r00**2*(dkappa + one/dkappa)
         if (lpprof) write (nthreed, 60) presfac*pfac, (-100.*(r00 -
     1      rthompeak))
         write (nthreed, 65) b0, cur0/dmu0
   60    format(//,' Input pressure scaled by ',f6.3/,
     1      ' Pressure profile shifted ',f6.2,' cm.',
     2      ' relative to magnetic axis')
   65    format(//' B-PHI(R=Raxis,Z=0) = ',f8.3,' [Wb/M**2]',/,
     1   ' J-PHI(R=Raxis,Z=0) = iota(0)*Bt0/(mu0*R0)*(1+k**2)/k = ',1p
     2   e10.3,' [A/M**2]',2/,' Comparison of Input vs Output',8x,
     3   'Input',9x,'VMEC Output',/,1x,71('-'))
         write (nthreed, 70) 1.e-6_dp*currv/dmu0, 1.e-6_dp*ctor/dmu0,
     1      phiedge, torflux, 1.e3_dp*phidiam, 
     1      1.e3_dp*delphid, 1.e-3_dp*pthommax,
     2      1.e-3_dp*ppeak/dmu0, rcenin, rcen, 
     2      raxis0, r00, zaxis0, z00, rc0mse,
     3      apres, aminr1, aminr2in, aminr2, bminz2in, bminz2, rinner,
     4      router, 1.e3_dp*delphid_exact
   70    format(' Toroidal Current     =',2(10x,f10.3),'  [MA]',/,
     1   ' Edge Toroidal Flux   =',2(10x,f10.3),'  [Wb]',/,
     2   ' Diamagnetic Flux(*)  =',2(10x,f10.3),'  [mWb]',/,
     3   ' Peak Pressure        =',2(10x,f10.3),'  [KPa]',/,
     4   ' Geometric Center     =',2(10x,f10.3),'  [M]',/,
     5   ' Magnetic R Axis      =',2(10x,f10.3),'  [M]',/,
     6   ' Magnetic Z Axis      =',2(10x,f10.3),'  [M]',/,
     7   ' MSE R Axis           =',10x,f10.3,20x,'  [M]',/,
     8   ' Minor Radius (apres) =',2(10x,f10.3),'  [M]',/,
     9   ' Minor Radius (a)     =',2(10x,f10.3),'  [M]',/,
     .   ' Minor Radius (b)     =',2(10x,f10.3),'  [M]',/,
     1   ' Inboard  Midplane R  =',30x,f10.3,'  [M]',/,
     2   ' Outboard Midplane R  =',30x,f10.3,'  [M]',/,50('-')/,
     3   ' * Exact diamagnetic flux (not based on equilibrium) = ',f8.3,
     4   '  [mWb]')

!       Calculate components of total chi-squared
         total_chi_cur = (currv - ctor)**2/sigma_current**2
         nchicur = 1
         if (iphidiam .eq. 1) total_chi_dia = (phidiam - delphid)**2/
     1      sigma_delphid**2
         nchidia = 1

      else                               ! mcz
         b0 = 0
      end if   !!if(lrecon)

      rmax_surf = maxval(r1(ns:nrzt:ns,0)+r1(ns:nrzt:ns,1))
      rmin_surf = minval(r1(ns:nrzt:ns,0)+r1(ns:nrzt:ns,1))
      zmax_surf = maxval(z1(ns:nrzt:ns,0)+z1(ns:nrzt:ns,1))

!
!     Calculate poloidal circumference and normal surface area
!
      allocate (guu_1u(nznt), guus1u(nznt))
      guu_1u(:nznt) = ru0(ns:nrzt:ns)*ru0(ns:nrzt:ns) +
     1   zu0(ns:nrzt:ns)*zu0(ns:nrzt:ns)

c     *** old version ***
c     guus1u(:nznt) = twopi*wint(ns:nrzt:ns)*sqrt(guu_1u(:nznt))
c     circum = sum(guus1u(:nznt))
c     sarea = twopi*dot_product((r1(ns:nrzt:ns,0) + r1(ns:nrzt:ns,1)),
c    1   guus1u(:nznt))

c    *** new version ***
      guus1u(:nznt) = wint(ns:nrzt:ns)*sqrt(guu_1u(:nznt))
      circum = twopi*sum(guus1u(:nznt))
      guus1u(:nznt) = wint(ns:nrzt:ns)*sqrt(
     1                (r1(ns:nrzt:ns,0)+r1(ns:nrzt:ns,1))**2
     2                    *guu_1u(:nznt) +
     3                 ((rv(ns:nrzt:ns,0)+rv(ns:nrzt:ns,1))*
     3                   zu0(ns:nrzt:ns) -
     4                  (zv(ns:nrzt:ns,0)+zv(ns:nrzt:ns,1))*
     5                   ru0(ns:nrzt:ns))**2 )
      sarea = twopi**2*sum(guus1u(:nznt))

      deallocate (guu_1u, guus1u)

      aspect = aspectratio()
      do js = 2, ns
         modb(:nznt) = sqrt(two*(bsq(js:nrzt:ns)-pres(js)))
         call bextrema (modb, bmin(1,js), bmax(1,js), nzeta, ntheta2)
      end do

!
!     output geometrical, |B| quantities
!
      call elongation (r1, z1, waist, height)

      write (nthreed, 75) bmin(1,ns), bmax(1,ns), bmin(ntheta2,ns), bmax
     1   (ntheta2,ns)
   75 format(/
     1   ' Magnetic field modulation (averaged over toroidal angle)',/,
     2   1x,71('-')/,' Bmin(u=0)             = ',f14.6/
     3   ' Bmax(u=0)             = ',f14.6/' Bmin(u=pi)            = ',
     4   f14.6/' Bmax(u=pi)            = ',f14.6/)

      sumbto = two*(vnorm*sum(bsq(:nrzt)*tau(:nrzt)) - sump)
      VolAvgB = sqrt(abs(sumbto/volf))
      IonLarmor = 0.0032_dp/VolAvgB
      jPS2(2:ns1) = (jpar2(2:ns1) - jdotb(2:ns1)**2/bdotb(2:ns1))
      s2 = sum(abs(jperp2(2:ns1))*(vp(2:ns1) + vp(3:ns)))
      jpar_perp = sum(jpar2(2:ns1)*(vp(2:ns1) + vp(3:ns)))
      jparPS_perp = sum(jPS2(2:ns1)*(vp(2:ns1) + vp(3:ns)))
      s2 = sum(jperp2(2:ns1)*(vp(2:ns1) + vp(3:ns)))
      if (s2 .ne. zero) then
         jpar_perp = jpar_perp/s2
         jparPS_perp = jparPS_perp/s2
      end if
      if (ntor .gt. 1) then
      write (nthreed, 80) aspect, volf, area, sarea, circum, Rmajor_p,
     1   Aminor_p, rmin_surf, rmax_surf, zmax_surf, waist(1), height(1),
     2   waist(2), height(2)
      else
      write (nthreed, 80) aspect, volf, area, sarea, circum, Rmajor_p,
     1   Aminor_p, rmin_surf, rmax_surf, zmax_surf, waist(1), height(1)
      end if
 80   format(/,' Geometric and Magnetic Quantities',/,1x,71('-')/,
     1   ' Aspect Ratio          = ',f14.6,/' Plasma Volume         = ',
     2   f14.6,' [M**3]',/' Cross Sectional Area  = ',f14.6,' [M**2]',/
     3   ' Normal Surface Area   = ',f14.6,' [M**2]',/
     4   ' Poloidal Circumference= ',f14.6,' [M]',/
     5   ' Major Radius          = ',f14.6,' [M]',
     6   ' (from Volume and Cross Section)',/
     7   ' Minor Radius          = ',f14.6,' [M]',
     8   ' (from Cross Section)',/
     9   ' Minimum (inboard)  R  = ',f14.6,' [M]',/
     A   ' Maximum (outboard) R  = ',f14.6,' [M]',/
     A   ' Maximum height     Z  = ',f14.6,' [M]',/
     B   ' Waist (v = 0)   in R  = ',f14.6,' [M]',/
     B   ' Full Height(v = 0)    = ',f14.6,' [M]',:,/
     B   ' Waist (v = pi)  in R  = ',f14.6,' [M]',:,/
     B   ' Full Height(v = pi)   = ',f14.6,' [M]')
      write (nthreed, 85) VolAvgB, IonLarmor, jpar_perp, jparPS_perp
 85   format(
     1   ' Volume Average B      = ',f14.6,' [T]',/
     2   ' Ion Larmor Radius     = ',f14.6,' [M] X Ti(keV)**0.5',/
     3   ' <J||**2>/<J-perp**2>  = ',f14.6,' (Vol. Averaged)',/
     4   ' <JPS**2>/<J-perp**2>  = ',f14.6,' (Vol. Averaged)',/)

      write (nthreed, 90)
   90 format(//' More Geometric and Physics Quantities',1x,71('-')/,5x,
     1   'j',3x,'psi-psiaxis',9x,'a [M]',3x,'ellipticity',3x,
     2   'indentation',7x,'d-shape',9x,'shift',6x,'<J||**2>/',4x,
     3   '<JPS**2>/',/,95x,
     4   '<J-perp**2>',3x,'<J-perp**2>'/,' -----',8(2x,12('-')))

      fac = twopi*hs*signgs
      psi(1) = zero
      allocate (r3v(ns-1))
      r3v(:ns-1) = fac*phip(2:ns)*iotas(2:ns)
      do i = 1, ns - 1
         psi(1+i) = psi(i) + r3v(i)
      end do
      deallocate (r3v)

      ygeo(1) = zero
      do js = 2, ns
         zmin =  HUGE(zmin)
         zmax = -HUGE(zmax)
         xmin =  HUGE(xmin)
         xmax = -HUGE(xmax)
         rzmax = zero

c                            !Theta = 0 to pi in upper half of X-Z plane
         noff = 1            !!nphi-plane, noff = 1,....,nzeta
         do icount = 1,2
            n1 = noff        !!nphi-plane, n1 = noff,...,nzeta
            if (icount .eq. 2)
     1      n1 = mod(nzeta + 1 - noff,nzeta) + 1           !!(twopi-v)
            loff = js + ns*(n1-1)
            t1 = one
            if (icount .eq. 2) t1 = -one
            do itheta = 1,ntheta2
               yr1u = r1(loff,0) + sqrts(js)*r1(loff,1)
               yz1u = z1(loff,0) + sqrts(js)*z1(loff,1)
               yz1u = t1*yz1u
               if (yz1u .ge. zmax) then
                  zmax = abs(yz1u)
                  rzmax = yr1u
               else if (yz1u .le. zmin) then
                  zmin = yz1u
                  rzmin = yr1u
               end if
               if (yr1u .ge. xmax) then
                  xmax = yr1u
                  zxmax = yz1u
               else if (yr1u .le. xmin) then
                  xmin = yr1u
                  zxmin = yz1u
               end if
               loff = loff + ns*nzeta
            end do
         end do


         lpi = ns*nzeta*(ntheta2 - 1)
         xmida = r1(js+lpi,0) + sqrts(js)*r1(js+lpi,1)
         xmidb = r1(js,0)     + sqrts(js)*r1(js,1)
!
         rgeo = p5*(xmidb + xmida)              !Geometric major radius
         ygeo(js) = p5*(xmidb - xmida)          !Geometric minor radius
c
         yinden(js) = (xmida - xmin)/(xmax - xmin) !Geometric indentation
         yellip(js) = (zmax - zmin)/(xmax - xmin)  !Geometric ellipticity
c
         ytrian(js) = (rgeo - rzmax)/(xmax - xmin) !Geometric triangularity
         yshift(js) = (r1(1,0)-rgeo)/(xmax - xmin) !Geometric shift
c
         if (jperp2(js) .eq. zero) jperp2(js) = epsilon(jperp2(js))
         jpar_perp = jpar2(js)/jperp2(js)
         if (js .lt. ns) then
            jparPS_perp = jPS2(js)/jperp2(js)
         else
            jparPS_perp = zero
         end if
         write (nthreed, 120) js, psi(js), ygeo(js), yellip(js),
     1      yinden(js), ytrian(js), yshift(js), jpar_perp, jparPS_perp

      end do
  120 format(1x,i5,6f14.5,1p3e14.2)

      write (nthreed, 130)
  130 format(//,' Magnetic Fields and Pressure',/,1x,71('-'))
      sumbpo = zero
      sumbtr = zero
      do i = 2, nrzt
         js = mod(i - 1,ns) + 1
         ncol = (i - 1)/ns + 1
         btorsq = (r12(i)*bsupv(i))**2
         bsq1 = bsupu(i)*bsubu(js,ncol,0) + bsupv(i)*bsubv(js,ncol,0)
         sumbpo = sumbpo + vnorm*tau(i)*(bsq1 - btorsq)
         sumbtr = sumbtr + vnorm*tau(i)*btorsq
      end do
      fac = p5/dmu0
      write (nthreed, 140) sump/dmu0, pavg/dmu0, fac*sumbpo, fac*sumbpo/
     1   volf, fac*sumbtr, fac*sumbtr/volf, fac*sumbto, fac*sumbto/volf,
     2   c1p5*sump/dmu0, c1p5*pavg/dmu0
  140 format(' Volume Integrals (Joules) and Volume ',
     1   'Averages (Pascals)',/,24x,'Integral',6x,'Average',/,
     2   ' pressure         = ',1p2e14.6,/,' bpol**2 /(2 mu0) = ',
     3   1p2e14.6,/,' btor**2/(2 mu0)  = ',1p2e14.6,/,
     4   ' b**2/(2 mu0)     = ',1p2e14.6,/,' EKIN (3/2p)      = ',
     5   1p2e14.6,/)

      write (nthreed, 800)
  800 format(/,' MAGNETIC AXIS COEFFICIENTS'/,
     1   '    n     rmag       zmag        rmag        zmag',/,
     2   '          (cc)       (cs)        (cs)        (cc)',/)
      n1 = 1
      zmags(1,n1) = zero                         !Used for pushing PHIEDGE
      do n = 0, ntor
         n1 = n + 1
         t1 = mscale(0)*nscale(n)
         if (lasym) then
            write (nthreed, 820) n, t1*rmags(1,n1), (-t1*zmags(1,n1)),
     1                             -t1*rmaga(1,n1),   t1*zmaga(1,n1)
         else
            write (nthreed, 820) n, t1*rmags(1,n1), (-t1*zmags(1,n1))
         end if
      end do
  820 format(i5,1p4e12.4)

      betpol = two*sump/sumbpo
      sump20 = two*sump
      sump2 = zero
      sump2 = sum(pres(2:ns)*pres(2:ns)*vp(2:ns)*vnorm)

      betstr = two*sqrt(sump2/volf)/(sumbto/volf)
      betatot = sump20/sumbto
      betapol = betpol
      betator = sump20/sumbtr

      write (nthreed, 150) betatot, betapol, betator
  150 format(' From volume averages over plasma, betas are',/,
     1   ' beta total    = ',f14.6,/,' beta poloidal = ',f14.6,/,
     2   ' beta toroidal = ',f14.6,/)

      write (nthreed, 160) 1.e-6_dp*ctor/dmu0, rbtor, betaxis, betstr
  160 format(' Toroidal Current     = ',f14.6,'  [MA]',/
     1   ' R * Btor-vac         = ',f14.6,' [Wb/M]',/,
     2   ' Peak Beta            = ',f14.6,/,' Beta-star            = ',
     3   f14.6,/)

      if (lrecon) then
!
!
!     Shafranov surface integrals s1,s2
!     Plasma Physics vol 13, pp 757-762 (1971)
!     Also, s3 = .5*S3, defined in Lao, Nucl. Fusion 25, p.1421 (1985)
!     Note: if ctor = 0, use Int(Bsupu*Bsubu dV) for ctor*ctor/R
!
      if (lfreeb) then
        factor = zero                             !Compute current-like norm
        do l = ns, nrzt, ns
           js = mod(l - 1,ns) + 1
           ncol = (l - 1)/ns + 1
           factor = factor + twopi*wint(l)*abs(bsubu(js,ncol,0))
        end do
        factor = one/factor**2
        facnorm = factor*twopi*twopi

        allocate (redg1u(nznt), rbps1u(nznt))
        redg1u(:nznt) = r1(ns:nznt*ns:ns,0) + r1(ns:nznt*ns:ns,1)
        rbps1u(:nznt) = two*facnorm*redg1u(:nznt)*(bpolvac(:nznt) +
     1    presf(ns))*wint(ns:nznt*ns:ns)
        sigr0 = dot_product(rbps1u(:nznt),zu0(ns:nznt*ns:ns))
        sigr1 = dot_product(rbps1u(:nznt)*zu0(ns:nznt*ns:ns),
     1                    redg1u(:nznt))
        sigz1 = -dot_product(rbps1u(:nznt)*ru0(ns:nznt*ns:ns),
     1           z1(ns:nznt*ns:ns,0) + z1(ns:nznt*ns:ns,1))
        deallocate (redg1u, rbps1u)

        er = sigr1 + sigz1
        rlao = volf/(twopi*area)               !LAO, NUCL.FUS.25(1985)1421
        flao = rshaf/rlao
        fgeo = rshaf/rcen
        factor = two*factor/rshaf

        smaleli = factor*sumbpo
        betai = two*factor*sump
        musubi = vnorm*factor*musubi
        dmusubi_meas = two*twopi*factor*delphid*rbtor
        lambda = p5*smaleli + betai
        s11 = (er - rshaf*sigr0)/rshaf         !Shafranov def. based on RT
        s12 = (er - rcen*sigr0)/rcen               !R = Rgeometric
        s13 = (er - rlao*sigr0)/rlao               !R = RLao
        s2 = sigr0
        s3 = sigz1/rshaf
        delta1 = zero
        delta2 = one - fgeo
        delta3 = one - flao
        write (nthreed, 170) rshaf, rcen, rlao, dmusubi_meas, delta1,
     1   delta2, delta3, s11, s12, s13, s2, s2, s2, s3, s3*fgeo, s3*flao
     2   , smaleli, smaleli*fgeo, smaleli*flao, musubi, musubi*fgeo,
     3   musubi*flao, betai, betai*fgeo, betai*flao, musubi + s11,
     4   musubi*fgeo + s12 + s2*(one - fgeo), musubi*flao + s13 + s2*(
     5   one - flao), lambda, fgeo*dlam, flao*dlam, 0.5*s11 + s2, 0.5*(
     6   s12 + s2*(one + fgeo)), 0.5*(s13 + s2*(one + flao)), (3*betai
     7    + smaleli - musubi)*0.5/(s11 + s2) - one, fgeo*(3*betai +
     8   smaleli - musubi)*0.5/(s12 + s2) - one, flao*(3*betai + smaleli
     9    - musubi)*0.5/(s13 + s2) - one, (betai + smaleli + musubi)*0.5
     .   /s2 - one, fgeo*(betai + smaleli + musubi)*0.5/s2 - one,
     .   flao*(betai + smaleli + musubi)*0.5/s2 - one
  170 format(' Integrals of Shafranov',/,1x,22('-'),/,
     1   ' RT (Flux-weighted)   = ',f14.6,' [M]',/,
     2   ' RG (Geometric)       = ',f14.6,' [M]',/,
     3   ' RL (Vol/2*pi*Area)   = ',f14.6,' [M]',/,
     4   ' Mui (diamagnetism)   = ',f14.6,2/,32x,'R = RT',12x,'R = RG',
     5   12x,'R = RL',/,20x,3(10x,8('-')),/,' delta = 1 - RT/R     = ',3
     6   (f14.6,4x),/,' s1                   = ',3(f14.6,4x),/,
     7   ' s2                   = ',3(f14.6,4x),/,
     8   ' s3                   = ',3(f14.6,4x),/,
     9   ' Li                   = ',3(f14.6,4x),/,
     .   ' Mui                  = ',3(f14.6,4x),/,
     1   ' Betai (Calculated)   = ',3(f14.6,4x),/,
     2   ' Betai (Mui + s1)     = ',3(f14.6,4x),/,
     3   ' Lambda (Calculated)  = ',3(f14.6,4x),/,
     4   ' Lambda (s1/2 + s2)   = ',3(f14.6,4x),/,
     5   ' 1st Shafr''v relation = ',3(f14.6,4x),/,
     6   ' (3*Betai + Li - Mui)/[2*(s1+s2)] - 1',/,
     7   ' Radial force balance = ',3(f14.6,4x),/,
     8   ' (Betai + Li + Mui)/(2*s2) - 1',/)

      end if
!
!     Safety Factors (q)
!
      qzero = one/abs(iotaf(1))
      qedge = one/abs(iotaf(ns))
      write (nthreed, 180) qzero, qedge, qzero/qedge
  180 format(' Safety Factors',/,1x,14('-'),/,' q (on axis) = ',f12.4,/,
     1   ' q (at edge) = ',f12.4,/,' q(0)/qedge  = ',f12.4,/)

!
!     PRINT OUT IOTA, PRESSURE SPLINE COEFFICIENTS
!     (PRESSURE IN MKS UNITS, NWT/M**2)
!
      write (nthreed, 190)
  190 format(/,' SUMMARY OF IOTA AND PRESSURE SPLINES'/,
     1   '  K   Spline Node       IOTA(K)      IOTA"(K)'/,
     2   '        sqrt(s)',/,3('-'),3(4x,10('-')))
      do i = 1, isnodes
         write (nthreed, 200) i, sknots(i), ystark(i), y2stark(i)
      end do
  200 format(i3,1p3e14.3)
      write (nthreed, 210)
  210 format(/,'  K   Spline Node       PRES(K)      PRES"(K)'/,
     1   '        sqrt(s)',/,3('-'),3(4x,10('-')))
      factor = pthommax
      do i = 1, ipnodes
         write (nthreed, 200) i, pknots(i), factor*ythom(i), factor*
     1      y2thom(i)
      end do

!
!     PRINT-OUT MAGNETIC PITCH DATA
!
      nchimse = imse
      total_mse_chi = zero

      if (imse .ne. 0) then

         write (nthreed, 220)
  220    format(//,4x,'N     R(data)      Spline',3x,
     1      'atan[BZ/BT] atan[BZ/BT] Chi-Sq-Err',4x,
     2      'Sigma       BZ(in)     BZ(calc)       BT        q(in)',/,
     3      12x,'[m]',5x,'sqrt(s)-node',2x,'(in-deg)',2x,'(calc-deg)',
     4      20x,3(9x,'[T]')/,2x,3('-'),1x,2(3x,9('-')),8(2x,10('-'))/)

         msewgt = zero
         denwgt = zero
         do n = 1, imse + 1
            nsort = isortr(n)
            pitchc = atan(starkcal(n))/dcon
            pitchm = atan(datastark(nsort))/dcon
            wght0 = atan(sigma_stark(nsort))/dcon
            js = indexs1(n)
            lt = indexu1(n)
            noff = ns*(lt - 1)
            rjs = r1(js+noff,0) + sqrts(js)*r1(js+noff,1)
            js1 = js + 1
            rjs1 = r1(js1+noff,0) + sqrts(js1)*r1(js1+noff,1)
            rcalc = (one - delso1(n))*rjs + delso1(n)*rjs1
            if (delse1(n) .lt. zero) rcalc = zero
            if (rcalc .ne. 0.) btor = fpsical(n)/rcalc
            bzin = tan(dcon*pitchm)*btor
            bzcalc = tan(dcon*pitchc)*btor
            chisq = zero
            if (abs(rcalc - router) .lt. epstan) then
               write (nthreed, 230) n, rcalc, rsort0(n), pitchm,
     1            pitchc, wght0, one/(qcalc(n) + epstan)
            else
               if (abs(rsort(n) - rcalc) .le. epstan) then
                  chisq = ((pitchc - pitchm)/wght0)**2
                  msewgt = msewgt + chisq
                  denwgt = denwgt + (pitchm/wght0)**2
               endif
               write (nthreed, 240) n, rcalc, rsort0(n), pitchm,
     1            pitchc, chisq, wght0, bzin, bzcalc, btor,
     2            one/(qmeas(n) + epstan)
            endif
         end do

         chiwgt = msewgt/(imse)
         msewgt = sqrt(msewgt/denwgt)
c                                             !total chi-squared for mse
         total_mse_chi = (imse)*chiwgt
         write (nthreed, 250) chiwgt, msewgt

  230    format(i5,' ',1p2e12.3,1p2e12.3,12x,1pe12.3,6x,
     1      '- Outer (phantom) Edge -',6x,1pe12.3)
  240    format(i5,' ',1p2e12.3,1p8e12.3)
  250    format(/' MEAN CHI-SQ ERROR IN STARK DATA MATCH : ',1pe10.3,/,
     1      ' RMS ERROR IN STARK DATA MATCH : ',1pe10.3,2/)

      endif

!
!     PRINT-OUT PRESSURE DATA
!
      nchipres = itse
      total_pres_chi = zero

      if (lpprof) then
         tswgt = zero
         denwgt = zero
         write (nthreed, 300) presfac*pfac
  300    format(4x,'N     R(in)      R(calc)',4x,f6.2,
     1      ' X Pres(in)      Pres(calc)','  Chi-Sq Err       Sigma'/,2x
     2      ,3('-'),2(2x,10('-')),2(6x,12('-')),2(3x,9('-')),/)

         do n = 1, itse
            js = indexs2(n)
            lt = indexu2(n)
            noff = ns*(lt - 1)
            rjs = r1(js+noff,0) + sqrts(js)*r1(js+noff,1)
            js1 = js + 1
            rjs1 = r1(js1+noff,0) + sqrts(js1)*r1(js1+noff,1)
            if (delso2(n) .eq. (-one)) delso2(n) = one
            rcalc = (one - delso2(n))*rjs + delso2(n)*rjs1
            wght0 = sigma_thom(n)
            chisq = ((datathom(n)*pfac-pcalc(n)/dmu0)/wght0)**2
            tswgt = tswgt + chisq
            denwgt = denwgt + (datathom(n)/wght0)**2
            write (nthreed, 310) n, rthom(n), rcalc, presfac*pfac*
     1         datathom(n), pcalc(n)/dmu0, chisq, wght0
         end do

         chiwgt = tswgt/(itse)
         tswgt = sqrt(tswgt/denwgt)
         total_pres_chi = (itse)*chiwgt !total
         write (nthreed, 320) chiwgt, tswgt
  310    format(i5,1p2e12.3,1p2e18.3,1p2e12.3)
  320    format(/' MEAN CHI-SQ ERROR IN PRESSURE DATA MATCH: ',1pe10.3,/
     1      ,' RMS ERROR IN PRESSURE DATA MATCH: ',1pe10.3/)
      endif

!
!     SUMMARIZE MAGNETICS DATA AND MATCH
!
      call magnetics_data

!
!     COMPUTE REAL TOROIDAL CURRENT ALONG MIDPLANE (MULTIPLY IZETA BY R/SQRT(G))
!     IN PHI = 0 PLANE
!
      call getcurmid (curmid, izeta, gsqrt, r12)

      do nout = nthreed, nmac, (nmac - nthreed)
         if (nout .eq. nmac) then
            if (.not.lmac) cycle
            write (nout, *)
     1      'FOLLOWING DATA EQUALLY SPACED IN TOROIDAL FLUX'
         end if
         write (nout, 700)
         if (nout .eq. nthreed) write (nout, 710)
         iotas(1) = two*iotas(2) - iotas(3)
         iotas(ns+1) = two*iotas(ns) - iotas(ns1)
         vp(1) = two*vp(2) - vp(3)
         vp(ns+1) = two*vp(ns) - vp(ns1)
         pres(1) = two*pres(2) - pres(3)
         pres(ns+1) = two*pres(ns) - pres(ns1)
         do icount = 1, 2*ns - 1
            js = mod(imid(icount) - 1,ns) + 1
            ageo(icount) = ygeo(js)
            phimid(icount) = torflux*(js - 1)/(ns1)
            psimid(icount) = psi(js)
            volpsi(icount) = vnorm*sum(vp(2:js))
            if (js .eq. 1) volpsi(icount) = zero
            dlogidsi = (iotas(js+1)-iotas(js))*ohs/iotaf(js)
            vprime = p5*vnorm*ohs*(vp(js)+vp(js+1))
            pprime = (pres(js+1)-pres(js))*ohs
            rgeo = rmid(icount) + ygeo(js)
            if (icount .gt. ns) rgeo = rgeo - two*ygeo(js)
            alfa(icount) = -two*pprime*vprime*sqrt(two*volpsi(icount)/
     1         rgeo/twopi)/(iotaf(js)*phipf(js))**2
            shear(icount) = -two*volpsi(icount)*dlogidsi/vprime
            write (nout, 720) rmid(icount), ygeo(js), psi(js),
     1         volpsi(icount), qmid(icount), shear(icount),
     2         presmid(icount),
     3         alfa(icount), curmid(icount), datamse(icount)
         end do
      end do
  700 format(/3x,'RMID[M]',6x,'AMID[M]',6x,'PSI[Wb]',4x,'VOL[M**3]',9x,
     1   'q(r)',2x,'SHEAR(Spsi)',6x,'P(PASC)',6x,'ALF(P'')',2x,
     2   'JTOR(A/M**2)',1x,'ATAN(Bz/Bt)[deg]')
  710 format(10('-'),9(3x,10('-')))
  720 format(1pe10.3,9(3x,1pe10.3))

!     Determine q-min on the s grid by min-splining on the sqrt-s grid
      nptsx = 2*ns - 1

      wmidx(:nptsx) = one
      tenmidx(:nptsx) = 0.1_dp
      rmidx(:nptsx) = rmid(:nptsx)
      qmidx(:nptsx) = qmid(:nptsx)

      tol = .005_dp
      sqmin = fmax(zero,one,splints,tol)
      iotamax = splints(sqmin)
      qmin3 = -99999.0_dp
      if (iotamax .ne. zero) qmin3 = one/iotamax
      sqmin = sqmin**2

!     Determine q-min on a fine r-midplane grid
      tol = .005_dp
!     outboard side
      rqmin1 = fmin(rmid(ns),rmid(nptsx),splintx,tol)
      qmin1 = splintx(rqmin1)
      rqmin2 = fmin(rmid(1),rmid(ns),splintx,tol)!outboard side only
      qmin2 = splintx(rqmin2)
      write (nthreed, 730) qmin1, rqmin1, qmin2, rqmin2, qmin3, sqmin
  730 format(//' MINIMUM Q :     OUTBOARD SIDE: QMIN= ',f6.3,'  AT R= ',
     1   f6.3,/,'                  INBOARD SIDE: QMIN= ',f6.3,'  AT R= '
     2   ,f6.3,/,'                    IN S SPACE: QMIN= ',f6.3,
     3   '  AT S= ',f6.3,/)

      nchiiota0 = 0
      total_chi_square = total_b_chi + total_saddle_chi + total_pres_chi
     1    + total_mse_chi + total_chi_cur + total_chi_dia + chisqerr(
     2   islope0)                        !Total CHISQ of equilibrium fit
      nchitot = nbfldn + nchiiota0 + nchisaddle + nchipres + nchimse +
     1   nchicur + nchidia
      total_chi_square_n = total_chi_square/max(1,nchitot)
      write (nthreed, 900)
      write (nthreed, 901) (bloopnames(n),nbfld(n),b_chi(n),b_chi(n)/
     1   max(1,nbfld(n)),n=1,nbsets)
      write (nthreed, 902) nchisaddle, total_saddle_chi,
     1   total_saddle_chi/max(1,nchisaddle), nchipres, total_pres_chi,
     2   total_pres_chi/max(1,nchipres), nchimse, total_mse_chi,
     3   total_mse_chi/max(1,nchimse), nchicur, total_chi_cur,
     4   total_chi_cur/max(1,nchicur), nchidia, total_chi_dia,
     5   total_chi_dia, nchitot, total_chi_square, total_chi_square_n

  900 format(/,' CHI-SQUARED SUMMARY:',t25,'Type',t50,'Nchi',t65,'ChiSq'
     1   ,t84,'ChiSq/N'/,t25,'----',t50,'----',t65,'-----',t84,'-------'
     2   )
  901 format(t25,'B-Loops-',a,t50,i3,t60,1pe12.4,t80,1pe12.4)
  902 format(t25,'Saddle',t50,i3,t60,1pe12.4,t80,1pe12.4,/,t25,
     1   'Pressure',t50,i3,t60,1pe12.4,t80,1pe12.4,/,t25,'MSE',t50,i3,
     2   t60,1pe12.4,t80,1pe12.4,/,t25,'Ip',t50,i3,t60,1pe12.4,t80,1p
     3   e12.4,/,t25,'Diamagnetic Flux',t50,i3,t60,1pe12.4,t80,1pe12.4,/
     4   ,t25,'TOTAL',t50,i3,t60,1pe12.4,t80,1pe12.4)

      end if          !!IF(LRECON)


      end subroutine eqfor


      subroutine elongation (r1, z1, waist, height)
      use vmec_main
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), intent(out) :: waist(2), height(2)
      real(rprec), intent(in), dimension(ns,nzeta,ntheta3,0:1) :: r1, z1
      integer :: nv, n1
C-----------------------------------------------
!
!     Compute Waist thickness, Height in phi = 0, pi symmetry planes
!
      n1 = 0
      do nv = 1, nzeta/2+1
         if (nv.ne.1 .and. nv.ne.nzeta/2+1) cycle
         n1 = n1+1
         waist(n1) = (r1(ns,nv,1,0)       + r1(ns,nv,1,1)) -
     1               (r1(ns,nv,ntheta2,0) + r1(ns,nv,ntheta2,1))
         height(n1) = 2*maxval(z1(ns,nv,:,0) + z1(ns,nv,:,1))
      end do

      end subroutine elongation

      
      subroutine getcurmid (curmid, izeta, gsqrt, r12)
      use vmec_input, only: rprec, dp, nzeta
      use vmec_dim, only: ns, ns1, ntheta2
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec) :: curmid(2*ns)
      real(rprec) :: izeta(ns,nzeta,*), gsqrt(ns,nzeta,*), 
     1    r12(ns,nzeta,*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec) :: midcur(ns)
C-----------------------------------------------
!     THETA = pi, PHI = 0      
      midcur(2:ns) = r12(2:ns,1,ntheta2)/gsqrt(2:ns,1,ntheta2)
      
      curmid(1) = izeta(ns,1,ntheta2)*midcur(ns)
      curmid(2:ns1) = 0.5_dp*izeta(ns1:2:-1,1,ntheta2)*
     1                   (midcur(ns1:2:-1) + midcur(ns:3:-1))
  
!     THETA = 0, PHI = 0      
      midcur(2:ns) = r12(2:ns,1,1)/gsqrt(2:ns,1,1)

      curmid(ns+1:2*ns-1) = 0.5_dp*izeta(2:ns1,1,1)*
     1                   (midcur(2:ns1) + midcur(3:ns))

      curmid(ns) = 0.5_dp*(curmid(ns-1) + curmid(ns+1))
      curmid(2*ns) = 2*curmid(2*ns-1) - curmid(2*ns-2)

      end subroutine getcurmid
            

      subroutine getfsq(gcr, gcz, gnormr, gnormz, gnorm, mprecon)
      use vmec_main
      use vmec_params, only: ntmax
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer mprecon
      real(rprec) gnormr, gnormz, gnorm
      real(rprec), dimension(ns,mnsize,ntmax) :: gcr, gcz
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: jsmax
C-----------------------------------------------
      jsmax = ns1 + mprecon
      gnormr = gnorm * sum(gcr(:jsmax,:,:)**2)
      gnormz = gnorm * sum(gcz(:jsmax,:,:)**2)

      end subroutine getfsq
      

      subroutine jxbforce(bsupu, bsupv, bsubu, bsubv, bsubs, bsubsu,
     1   bsubsv, itheta, brho, gsqrt, izeta, bsq)
      use safe_open_mod
      use vmec_main
      use vmec_params, only: mscale, nscale, signgs
      use realspace
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,nznt), intent(in) ::
     1  bsupu, bsupv, bsq, gsqrt
      real(rprec), dimension(ns,nznt,0:1), intent(inout) ::
     1  bsubu, bsubv
      real(rprec), dimension(ns,nznt), intent(inout) :: bsubs
      real(rprec), dimension(ns,nznt), intent(out) ::
     1  itheta, brho, izeta
      real(rprec), dimension(ns,nznt,0:1), intent(out) ::
     1  bsubsu, bsubsv
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      logical, parameter :: lbsubs = .true.       !!True  to use NEW bsubs calculation (from mag. diff. eq.)
                                                  !!False to use OLD bsubs calculation (from metrics)
      logical, parameter :: lprint = .false.      !!Prints out bsubs spectrum to fort.33
      real(rprec), parameter :: two=2, p5=0.5_dp, c1p5=1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer lk, lz, lt, k, m, js, j, n, injxbout, mparity
      integer :: njxbout = jxbout0, mnyq, nnyq, nmin
      integer, parameter :: ns_skip = 1, nu_skip = 1, nv_skip = 1
      real(rprec), dimension(:,:), allocatable ::
     1    bdotj, bsubuv, bsubvu, brhomn
      real(rprec), dimension(:,:,:), allocatable :: bsubsmn
      real(rprec), dimension(:), allocatable     :: jperpu, jperpv, 
     2    sqgb2, jp2, jxb, jxb2, bsupu1, bsupv1, bsubu1, bsubv1
        real(rprec), dimension(:,:), allocatable   :: bsubua, bsubva
      real(rprec) ::
     1    bsubsmn1, bsubsmn2, bsubvmn1, bsubvmn2, bsubumn1, bsubumn2,
     1    bsubsmn3, bsubsmn4, bsubvmn3, bsubvmn4, bsubumn3, bsubumn4,
     2    dnorm1, tsini1, tsini2, tcosi1, tcosi2, tcosm1, tcosm2,
     3    tcosn1, tcosn2, tsinm1, tsinm2, tcos1, tcos2, tsin1, tsin2,
     4    tsinn1, tsinn2, avforce, pprime, amaxfor, aminfor, tjnorm,
     5    ovp, pnorm, brho00(ns)
      real(rprec), dimension(:,:), allocatable ::
     1    bsubs_s, bsubs_a, bsubu_s, bsubu_a, bsubv_s, bsubv_a
      character :: jxbout_file*100
      real(rprec), external :: dot_g
C-----------------------------------------------
      jxbout_file = 'jxbout.'//input_extension
      call safe_open(njxbout, injxbout, jxbout_file, 'replace',
     1    'formatted')
      if (injxbout .ne. 0) then
         print *,' Error opening JXBOUT file in jxbforce'
         return
      end if

      write (njxbout,6) (ns1-1)/ns_skip, ntheta2/nu_skip, nzeta/nv_skip,
     1    mpol, ntor
 6    format(/,' Radial surfaces = ',i3, ' Poloidal grid points = ',i3,
     1         ' Toroidal grid points = ',i3,/,
     2         ' Poloidal modes = ',i3,' Toroidal modes = ', i3)

!
!     PROGRAM FOR COMPUTING LOCAL JXB = grad-p FORCE BALANCE
!
!     Compute u (=theta), v (=zeta) derivatives of B sub s
!
      mnyq = max(0, ntheta1/2)
      nnyq = max(0, nzeta/2)
      write (njxbout, 5)
 5    format(' LEGEND:',/,
     1  " U = VMEC poloidal angle, V = VMEC (geometric) toroidal angle"/
     2  " SQRT(g') = SQRT(g-VMEC) / V': Jacobian based on VOLUME",/,
     3  " V' = dV/ds: differential volume element",/,
     4  " Es = SQRT(g') [grad(U) X grad(V)] : covariant radial",
     4  " unit vector based on volume",/,
     5  " BSUP{U,V} = B DOT GRAD{U,V}:  contravariant components of B"/,
     6  " JSUP{U,V} = SQRT(g') J DOT GRAD{U,V}",/,
     7  " J X B = Es DOT [J X B]: covariant component of J X B force",/,
     8  " J * B = J DOT B * SQRT(g')",/,
     9  " p' = dp/dV: pressure gradient based on volume derivative",//)

      lz = nzeta*ntheta2
      allocate (bdotj(ns,nznt), bsubuv(ns,nznt),
     1          bsubvu(ns,nznt), jperpu(nznt), jperpv(nznt),
     2          sqgb2(nznt), brhomn(0:mnyq,-nnyq:nnyq),jp2(nznt),
     3          bsubua(nznt,0:1), bsubva(nznt,0:1), jxb(nznt), 
     4          jxb2(nznt), bsupu1(nznt), bsupv1(nznt), bsubu1(nznt), 
     5          bsubv1(nznt), bsubsmn(ns,0:mnyq,-nnyq:nnyq),
     6          bsubs_s(lz,0:1), bsubs_a(lz,0:1),
     7          bsubu_s(lz,0:1), bsubu_a(lz,0:1),
     8          bsubv_s(lz,0:1), bsubv_a(lz,0:1), stat = j)
      if (j .ne. 0) stop 'allocation error in jxbforce'


!
!     NOTE: bsubuv, bsubvu are used to compute the radial current (should be zero)
!
      bsubsu = 0; bsubsv = 0; bsubuv = 0; bsubvu = 0; bdotj  = 0

!
!     IT IS ASSUMED THAT MSCALE, NSCALE = [1/sqrt(2), 1, ....] HERE
!     IF NOT, MUST PUT THESE M,N DEPENDENT NORMALIZATION FACTORS INTO DNORM1
!
      radial: do js = 1, ns
         if (js.gt.1 .and. js.lt.ns) then     !!Put on full mesh
            bsubs(js,:) = p5*(bsubs(js,:) + bsubs(js+1,:))
         end if
         bsubu(js,:,1) = bsubu(js,:,1)/shalf(js)
         bsubv(js,:,1) = bsubv(js,:,1)/shalf(js)
         bsubua = 0;   bsubva = 0

         if (lasym)  then
            call fsym_fft (bsubs(js,:), bsubu(js,:,:), bsubv(js,:,:),
     1             bsubs_s, bsubu_s, bsubv_s, bsubs_a, bsubu_a, bsubv_a)
         else
            bsubs_s(:,0) = bsubs(js,:); bsubu_s = bsubu(js,:,:)
            bsubv_s = bsubv(js,:,:)
         end if

         do m = 0, mnyq
            mparity = mod(m, 2)
            do n = 0, nnyq
               dnorm1 = 4*mscale(m)*nscale(n)        !! One factor mscale*nscale cancels in trig...
               if (m.eq.mnyq) dnorm1 = p5*dnorm1
               if (n.eq.nnyq .and. n.ne.0) dnorm1 = p5*dnorm1
               bsubsmn1 = 0;  bsubsmn2 = 0
               if (lasym) bsubsmn3 = 0;  bsubsmn4 = 0
               if (m.gt.mpol1 .or. n.gt.ntor) goto 222
               bsubumn1 = 0;  bsubumn2 = 0;  bsubvmn1 = 0;  bsubvmn2 = 0
               if (lasym)
     1         bsubumn3 = 0;  bsubumn4 = 0;  bsubvmn3 = 0;  bsubvmn4 = 0

               do k = 1, nzeta
                  lk = k
                  do j = 1, ntheta2
                     tsini1 = sinmui(j,m)*cosnv(k,n)*dnorm1
                     tsini2 = cosmui(j,m)*sinnv(k,n)*dnorm1
                     tcosi1 = cosmui(j,m)*cosnv(k,n)*dnorm1
                     tcosi2 = sinmui(j,m)*sinnv(k,n)*dnorm1
                     bsubsmn1 = bsubsmn1 + tsini1*bsubs_s(lk,0)
                     bsubsmn2 = bsubsmn2 + tsini2*bsubs_s(lk,0)
                     bsubvmn1 = bsubvmn1 + tcosi1*bsubv_s(lk, mparity)
                     bsubvmn2 = bsubvmn2 + tcosi2*bsubv_s(lk, mparity)
                     bsubumn1 = bsubumn1 + tcosi1*bsubu_s(lk, mparity)
                     bsubumn2 = bsubumn2 + tcosi2*bsubu_s(lk, mparity)

                     if (lasym) then
                     bsubsmn3 = bsubsmn3 + tcosi1*bsubs_a(lk,0)
                     bsubsmn4 = bsubsmn4 + tcosi2*bsubs_a(lk,0)
                     bsubvmn3 = bsubvmn3 + tsini1*bsubv_a(lk, mparity)
                     bsubvmn4 = bsubvmn4 + tsini2*bsubv_a(lk, mparity)
                     bsubumn3 = bsubumn3 + tsini1*bsubu_a(lk, mparity)
                     bsubumn4 = bsubumn4 + tsini2*bsubu_a(lk, mparity)
                     end if

                     lk = lk + nzeta
                  end do
               end do

               dnorm1 = one/(mscale(m)*nscale(n))

!
!              Compute on half u grid (must add symmetric, antisymmetric parts for lasym)
! 
               do k = 1, nzeta
                  lk = k
                  do j = 1, ntheta2
                     tcos1 = cosmu(j,m)*cosnv(k,n)*dnorm1
                     tcos2 = sinmu(j,m)*sinnv(k,n)*dnorm1
!
!                    MUST DEALIAS BSUBU,V IF LBSUBS = TRUE
!
                     bsubua(lk,0) = bsubua(lk,0) + tcos1*bsubumn1 +
     1                  tcos2*bsubumn2
                     bsubva(lk,0) = bsubva(lk,0) + tcos1*bsubvmn1 +
     1                  tcos2*bsubvmn2

                     tcosm1 = cosmum(j,m)*cosnv(k,n)*dnorm1
                     tcosm2 = sinmum(j,m)*sinnv(k,n)*dnorm1
                     bsubsu(js,lk,0) = bsubsu(js,lk,0) +
     1                  tcosm1*bsubsmn1 + tcosm2*bsubsmn2
                     tcosn1 = sinmu(j,m)*sinnvn(k,n)*dnorm1
                     tcosn2 = cosmu(j,m)*cosnvn(k,n)*dnorm1
                     bsubsv(js,lk,0) = bsubsv(js,lk,0) +
     1                  tcosn1*bsubsmn1 + tcosn2*bsubsmn2
                     bsubvu(js,lk) = bsubvu(js,lk) + dnorm1*(
     1                               sinmum(j,m)*cosnv(k,n)*bsubvmn1 +
     2                               cosmum(j,m)*sinnv(k,n)*bsubvmn2)
                     bsubuv(js,lk) = bsubuv(js,lk) + dnorm1*(
     1                               cosmu(j,m)*sinnvn(k,n)*bsubumn1 +
     2                               sinmu(j,m)*cosnvn(k,n)*bsubumn2)

                     if (lasym) then
                     tsin1 = sinmui(j,m)*cosnv(k,n)*dnorm1
                     tsin2 = cosmui(j,m)*sinnv(k,n)*dnorm1
                     bsubua(lk,1) = bsubua(lk,1) + tsin1*bsubumn3 +
     1                  tsin2*bsubumn4
                     bsubva(lk,1) = bsubva(lk,1) + tsin1*bsubvmn3 +
     1                  tsin2*bsubvmn4

                     tsinm1 = sinmum(j,m)*cosnv(k,n)*dnorm1
                     tsinm2 = cosmum(j,m)*sinnv(k,n)*dnorm1
                     bsubsu(js,lk,1) = bsubsu(js,lk,1) +
     1                   tsinm1*bsubsmn3 + tsinm2*bsubsmn4
                     tsinn1 = cosmu(j,m)*sinnvn(k,n)*dnorm1
                     tsinn2 = sinmu(j,m)*cosnvn(k,n)*dnorm1
                     bsubsv(js,lk,1) = bsubsv(js,lk,1) +
     1                   tsinn1*bsubsmn3 + tsinn2*bsubsmn4
                     end if
                     lk = lk + nzeta
                  end do
               end do

 222           continue

!
!        FIX MULTIPLIER FOR M != 0, N != 0 (p5)
!
               if (m .eq. 0) then
                  bsubsmn(js,m,n) =-bsubsmn2
                  if (n.ne.0) bsubsmn(js,m,-n)= 0
               else if (n .eq. 0) then
                  bsubsmn(js,m,n) = bsubsmn1
               else
                  bsubsmn(js,m,n)  = p5*(bsubsmn1 - bsubsmn2)
                  bsubsmn(js,m,-n) = p5*(bsubsmn1 + bsubsmn2)
               end if

            end do
         end do

         if (lasym) call fsym_invfft (bsubua, bsubva, 1)
         bsubu(js,:,0) = bsubua(:,0)
         bsubv(js,:,0) = bsubva(:,0)

      end do radial

      if (lasym) call fsym_invfft (bsubsu, bsubsv, ns)

      if (.not.lbsubs) go to 1500          !!SKIPS Bsubs Correction - uses Bsubs from metric elements

!
!     Compute corrected B-sub-s (impacts currents)
!
      correct_bsubs: do js = 2, ns-1
         jxb(:) = p5*(gsqrt(js,:) + gsqrt(js+1,:))
         bsupu1(:) = p5*(bsupu(js,:)*gsqrt(js,:)
     1    + bsupu(js+1,:)*gsqrt(js+1,:))
         bsupv1(:) = p5*(bsupv(js,:)*gsqrt(js,:)
     1    + bsupv(js+1,:)*gsqrt(js+1,:))
         brho(js,:) = ohs*
     1   ( bsupu1(:)*(bsubu(js+1,:,0) - bsubu(js,:,0))
     2   + bsupv1(:)*(bsubv(js+1,:,0) - bsubv(js,:,0)))
     3   + (pres(js+1) - pres(js))*ohs*jxb(:)
         brho00(js) = DOT_G(nznt,brho(js,1),ns,wint(js),ns)
         brho(js,:) = brho(js,:) - signgs*jxb(:)*brho00(js)/
     1      (p5*(vp(js) + vp(js+1)))

         jxb(:) = brho(js,:)
         call getbrho (brhomn, jxb, bsupu1, bsupv1, mnyq, nnyq)
         if (lprint) then
            write (33, *) ' JS = ', js
            write (33, *) '  M    N       BSUBS(old)        BSUBS(new)'
            do m = 0, mnyq
               nmin = -nnyq
               if (m .eq. 0) nmin = 0
               do n = nmin, nnyq
                  write(33,1223) m, n, bsubsmn(js,m,n), brhomn(m,n)
               end do
            end do
         end if
 1223    format (i4,1x,i4,2(6x,1pe12.3))

!
!        RECOMPUTE bsubsu,v now using corrected bsubs
!
         itheta(js,:) = bsubsu(js,:,0)        !!Store old values here
         izeta (js,:) = bsubsv(js,:,0)
         bsubsu(js,:,:) = 0
         bsubsv(js,:,:) = 0

         do m = 0, mnyq
            do n = 0, nnyq
               dnorm1 = one/(mscale(m)*nscale(n))
               if (n .eq. 0) then
                  bsubsmn1 = brhomn(m,n)
                  bsubsmn2 = 0
               else
                  bsubsmn1 = brhomn(m,n) + brhomn(m,-n)
                  bsubsmn2 =-brhomn(m,n) + brhomn(m,-n)
               end if
               do k = 1, nzeta
                  lk = k
                  do j = 1, ntheta2
                     tcosm1 = cosmum(j,m)*cosnv(k,n)*dnorm1
                     tcosm2 = sinmum(j,m)*sinnv(k,n)*dnorm1
                     bsubsu(js,lk,0) = bsubsu(js,lk,0) +
     1                  tcosm1*bsubsmn1 + tcosm2*bsubsmn2
                     tcosn1 = sinmu(j,m)*sinnvn(k,n)*dnorm1
                     tcosn2 = cosmu(j,m)*cosnvn(k,n)*dnorm1
                     bsubsv(js,lk,0) = bsubsv(js,lk,0) +
     1                  tcosn1*bsubsmn1 + tcosn2*bsubsmn2
                     lk = lk + nzeta
                  end do
               end do
            end do
         end do

      end do correct_bsubs

      if (lasym) call fsym_invfft (bsubsu, bsubsv, ns)

!
!     CHECK FORCE BALANCE: sqrt(g)*(bsupu*bsubsu + bsupv*bsubsv) = brho
!
      check_fb: do js = 2, ns-1
         bsupu1(:) = p5*(bsupu(js,:)*gsqrt(js,:)
     1    + bsupu(js+1,:)*gsqrt(js+1,:))
         bsupv1(:) = p5*(bsupv(js,:)*gsqrt(js,:)
     1    + bsupv(js+1,:)*gsqrt(js+1,:))
!        jp2(:) = bsupu1(:)*bsubsu(js,:,0) + bsupv1(:)*bsubsv(js,:,0)
!        jxb(:) = bsupu1(:)*itheta(js,:) + bsupv1(:)*izeta(js,:)

!        print *,'JS = ',js
!        do lk = 1, nznt
!           write(*,1224) lk, brho(js,lk), jxb(lk), jp2(lk)
!        end do

!        pause

      end do check_fb

 1224    format ('lk = ',i4,' brho(rhs) = ', 1pe12.4,
     1   ' B dot grad Bs(old) = ', 1pe12.4, ' B dot grad Bs(new) = ',
     2   1pe12.4)

 1500 continue

      deallocate (bsubs_s, bsubs_a, bsubu_s, bsubu_a, bsubv_s, bsubv_a)
!
!     Now compute currents on the FULL radial mesh
!     Itheta = sqrt(g) * Jsupu, Izeta = sqrt(g) * Jsupv,
!     where Jsupx = J dot grad(x)
!     jxb = J X B   bdotj = sqrt(g)*J dot B
!     jperp-x = (B X gradp) dot grad(x) / |B|**2, x=(u,v)
!     Here, we compute |j-perp|**2 = (j-perp-u)**2 * guu + ...
!     This was compared to the alternative expression (agreed very well):
!     |j-perp|**2 = |grad-s|**2 * (dp/ds)**2 / |B|**2
!
!     Note: Multiply currents, pressure by 1/mu0 to get in mks units!
!     Note: IZETA <=> GVV OVERWRITES GVV!!
!
      do js = 2, ns1
         ovp = two/(vp(js+1) + vp(js))
         tjnorm = ovp*signgs
         pprime = ovp*ohs*(pres(js+1)-pres(js))/dmu0
         sqgb2(:nznt) = (gsqrt(js+1,:nznt)*(bsq(js+1,:nznt)- pres(js+1))
     1                +  gsqrt(js,:nznt)  *(bsq(js,:nznt) - pres(js)))
         jperpu(:nznt) = p5*(bsubv(js+1,:nznt,0) +
     1                       bsubv(js,:nznt,0))*pprime/sqgb2
         jperpv(:nznt) =-p5*(bsubu(js+1,:nznt,0) +
     1                       bsubu(js,:nznt,0))*pprime/sqgb2
         jp2(:nznt)=p5*(jperpu**2*(guu(js+1:nrzt:ns) + guu(js:nrzt:ns))
     1        + two*jperpu*jperpv*(guv(js+1:nrzt:ns) + guv(js:nrzt:ns))
     2        +         jperpv**2*(gvv(js+1:nrzt:ns) + gvv(js:nrzt:ns)))
         itheta(js,:nznt) = bsubsv(js,:nznt,0)
     1                   - ohs*(bsubv(js+1,:nznt,0) - bsubv(js,:nznt,0))
         izeta(js,:nznt) = -bsubsu(js,:nznt,0)
     1                   + ohs*(bsubu(js+1,:nznt,0) - bsubu(js,:nznt,0))
         itheta(js,:nznt) = itheta(js,:nznt)/dmu0
         izeta(js,:nznt)  = izeta(js,:nznt)/dmu0
         jxb(:) = p5*(gsqrt(js,:) + gsqrt(js+1,:))
         bsupu1(:nznt) = p5*(bsupu(js+1,:nznt)*gsqrt(js+1,:)
     1                 +     bsupu(js,:nznt)  *gsqrt(js,:))  / jxb(:)
         bsupv1(:nznt) = p5*(bsupv(js+1,:nznt)*gsqrt(js+1,:)
     1                 +     bsupv(js,:nznt)  *gsqrt(js,:))  / jxb(:)
         bsubu1(:nznt) = p5*(bsubu(js+1,:nznt,0) + bsubu(js,:nznt,0))
         bsubv1(:nznt) = p5*(bsubv(js+1,:nznt,0) + bsubv(js,:nznt,0))
         jxb(:nznt) = ovp*(itheta(js,:nznt) * bsupv1(:nznt)
     1              -      izeta (js,:nznt) * bsupu1(:nznt))
         bdotj(js,:nznt) = itheta(js,:nznt) * bsubu1(:nznt) +
     1                     izeta (js,:nznt) * bsubv1(:nznt)
         pnorm = one/(abs(pprime) + epsilon(pprime))
         amaxfor = max(maxval(jxb(:nznt)-pprime)*pnorm, zero)
         aminfor = min(minval(jxb(:nznt)-pprime)*pnorm, zero)
         avforce = sum(wint(2:nrzt:ns)*(jxb(:nznt) - pprime))


!        Compute <j dot B>, <B sup v> = signgs*phip
!        jpar2 = <j||**2>, jperp2 = <j-perp**2>

         jdotb(js) = tjnorm*sum(bdotj(js,:nznt)*wint(2:nrzt:ns))
         bdotb(js) = tjnorm*sum(sqgb2(:nznt)*wint(2:nrzt:ns))
         bdotgradv(js) = p5*(phip(js) + phip(js+1))*tjnorm
         jpar2(js) = tjnorm*sum(bdotj(js,:nznt)**2*wint(2:nrzt:ns)
     1                         /sqgb2(:nznt))
         jperp2(js) = tjnorm*sum(jp2(:nznt)*wint(2:nrzt:ns)*
     1                        p5*(gsqrt(js+1,:nznt) + gsqrt(js,:nznt)))

         if (mod(js,ns_skip) .eq. 0) then
            amaxfor = min(amaxfor,9.999_dp)
            aminfor = max(aminfor,-9.999_dp)
            write (njxbout, 200) phi(js)/phi(ns), avforce, jdotb(js),
     1         bdotgradv(js), 100.0_dp*amaxfor, 100.0_dp*aminfor
            write (njxbout, 90)
            do lz = 1, nzeta, nv_skip
               write (njxbout, 100) 360._dp*(lz-1)/nzeta, lz
               do lt = 1, ntheta2, nu_skip
                  lk = lz + nzeta*(lt - 1)
                  write (njxbout, 110) lt, ovp*itheta(js,lk),
     1              ovp*izeta(js,lk), ovp*(bsubuv(js,lk) -
     2              bsubvu(js,lk))/dmu0, bsupu1(lk), bsupv1(lk),
     3              jxb(lk), pprime,
     4        (jxb(lk) - pprime), ovp*bdotj(js,lk), bsubu(js,lk,0),
     5              bsubv(js,lk,0), bsubs(js,lk)
               end do
            end do
         endif
      end do

      close (njxbout)

      izeta(1,:nznt) = two*izeta(2,:nznt) - izeta(3,:nznt)           !!For print out in wrout
      izeta(ns,:nznt)= two*izeta(ns-1,:nznt) - izeta(ns-2,:nznt)     !!For print out in wrout
      jdotb(1) = two*jdotb(2) - jdotb(3)
      jdotb(ns) = two*jdotb(ns-1) - jdotb(ns-2)
      bdotb(1) = two*bdotb(3) - bdotb(2)
      bdotb(ns) = two*bdotb(ns-1) - bdotb(ns-2)
      bdotgradv(1) = two*bdotgradv(2) - bdotgradv(3)
      bdotgradv(ns) = two*bdotgradv(ns-1) - bdotgradv(ns-2)
      jpar2(1)   = 0; jpar2(ns)  = 0; jperp2(1)  = 0; jperp2(ns) = 0

      deallocate (jperpu, jperpv, sqgb2, jp2, brhomn, bsubsmn, bsubua,
     1    bsubva, jxb, jxb2, bsupu1, bsupv1, bsubu1, bsubv1, stat = j)
!
!     COMPUTE MERCIER CRITERION
!
      bdotj = dmu0*bdotj
      call Mercier(gsqrt,bsq,bdotj,iotas,wint,r1,ru,rv,zu,zv,bsubu,
     1             vp,phips,pres,ns,nznt)


   90 format(/"   LU      JSUPU      JSUPV      JSUPS      BSUPU",
     1   "      BSUPV      J X B       p'    J X B - p'     J * B",
     2   "      BSUBU      BSUBV      BSUBS   "/)
  100 format( " TOROIDAL ANGLE (PER PERIOD) = ", f8.3," DEGREES",
     1        " (PLANE #", i3,")")
  110 format(i5,1p12e11.3)
  200 format(/" TOROIDAL FLUX = ",1pe12.3,3x,"<J X B - p'> = ",
     1   1pe12.3,3x,"<J DOT B> = ",1pe12.3,3x,
     2   "<B DOT GRAD(V)> = ",1pe12.3,/,
     3   " MAXIMUM FORCE DEVIATIONS (RELATIVE TO p'): ",sp,0pf7.2,"%",
     4     3x,f7.2,"%")

      deallocate (bdotj, bsubuv, bsubvu, stat = j)

      end subroutine jxbforce


      subroutine fsym_fft (bs, bu, bv, bs_s, bu_s, bv_s, 
     1    bs_a, bu_a, bv_a)
      use vmec_main
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(nzeta,ntheta3), intent(in) :: bs
      real(rprec), dimension(nzeta,ntheta3,0:1), intent(in) :: bu, bv
      real(rprec), dimension(nzeta,ntheta2,0:1), intent(out) :: 
     1   bs_s, bu_s, bv_s, bs_a, bu_a, bv_a
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ir, i, jk, jka, mpar, kz, kr
C-----------------------------------------------
 
!
!       SYMMETRIZE CONTRAVARIANT COMPONENTS OF B
!       SO COS,SIN INTEGRALS CAN BE PERFORMED ON HALF-THETA INTERVAL
!
!       bs_s(v,u) = .5*( bs_s(v,u) - bs_s(-v,-u) )     ! * sin(mu - nv)
!       bs_a(v,u) = .5*( bs_s(v,u) + bs_s(-v,-u) )     ! * cos(mu - nv)
!
!     
      do mpar = 0, 1
         do i = 1, ntheta2
            ir = ntheta1 + 2 - i                 !-theta
            if (i == 1) ir = 1
            do kz = 1, nzeta
               kr = ireflect(ns*kz)/ns           !-zeta
               bs_a(kz,i,mpar) = cp5*(bs(kz,i)+bs(kr,ir))
               bs_s(kz,i,mpar) = cp5*(bs(kz,i)-bs(kr,ir))
               bu_a(kz,i,mpar) = cp5*(bu(kz,i,mpar)-bu(kr,ir,mpar))
               bu_s(kz,i,mpar) = cp5*(bu(kz,i,mpar)+bu(kr,ir,mpar))
               bv_a(kz,i,mpar) = cp5*(bv(kz,i,mpar)-bv(kr,ir,mpar))
               bv_s(kz,i,mpar) = cp5*(bv(kz,i,mpar)+bv(kr,ir,mpar))
            end do
         end do
      end do
     
      end subroutine fsym_fft


      subroutine fsym_invfft (bsubsu, bsubsv, ns)
      use vmec_main, only: rprec, nzeta, ntheta1, ntheta2, 
     1    ntheta3, ireflect
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: ns
      real(rprec), dimension(ns*nzeta,ntheta3,0:1),
     1   intent(inout) :: bsubsu, bsubsv
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ir, i, jkz, jkr
C-----------------------------------------------

      do i = 1 + ntheta2, ntheta1
         ir = ntheta1 + 2 - i                 !-theta
         do jkz= 1, ns*nzeta
            jkr = ireflect(jkz)               !-zeta
            bsubsu(jkz,i,0) = bsubsu(jkr,ir,0) - bsubsu(jkr,ir,1)
            bsubsv(jkz,i,0) = bsubsv(jkr,ir,0) - bsubsv(jkr,ir,1)
         end do
      end do

      bsubsu(:,:ntheta2,0)=bsubsu(:,:ntheta2,0) + bsubsu(:,:ntheta2,1)
      bsubsv(:,:ntheta2,0)=bsubsv(:,:ntheta2,0) + bsubsv(:,:ntheta2,1)

      end subroutine fsym_invfft


      subroutine getbrho (brhomn, frho, bsupu, bsupv, mmax, nmax)
      use kind_spec
      use vmec_input, only: nfp, nzeta
      use vmec_dim, only: ntheta1, ntheta2, ntheta3
      use vmec_persistent, only: cosmu, sinmu, cosnv, sinnv
      use vmec_params, only: mscale, nscale
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: mmax, nmax
      real(rprec), intent(out) :: brhomn(0:mmax, -nmax:nmax)
      real(rprec), dimension(nzeta, ntheta3), intent(in) ::
     1    bsupu, bsupv, frho
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: p5 = 0.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, j, m, n, nmax1, itotal, ijtot, mntot
      real(rprec) :: amn, dnorm, ccmn, ssmn, dm, dn, termsc, termcs
      real(rprec), allocatable :: amatrix(:,:), save_matrix(:,:),
     1   brhs(:)
      logical :: lpior0
C-----------------------------------------------
      external solver
C-----------------------------------------------
!
!     Solves the radial force balance B dot Brho = Fs for Brho in Fourier space
!     Here, Fs = frho(mn) is the Fourier transform sqrt(g)*F (part of radial force
!     balance sans the B dot Brho term)
!

      nmax1 = max(0,nmax-1)
      itotal = ntheta2*nzeta - 2*nmax1
      allocate (amatrix(itotal, itotal),
     1      brhs(itotal), save_matrix(itotal, itotal), stat=m)
      if (m .ne. 0) stop 'Allocation error in getbrho'
      

      amatrix = 0

!
!     BRHO = BSC(M,N)*SIN(MU)COS(NV) + BCS(M,N)*COS(MU)SIN(NV) 
!          + BCC(M,N)*COS(MU)COS(NV) + BSS(M,N)*SIN(MU)SIN(NV)   (ONLY IF ASYMMETRIC MODE ALLOWED)
!

      ijtot = 0
      brhs = 0
      do i = 1, ntheta2
         do j = 1, nzeta
!           IGNORE u=0,pi POINTS FOR v > pi: REFLECTIONAL SYMMETRY         
            lpior0 = ((i.eq.1 .or. i.eq.ntheta2) .and. (j.gt.nzeta/2+1))
            if (lpior0) cycle
            ijtot = ijtot + 1
            brhs(ijtot) = frho(j,i)
            mntot = 0
            do m = 0, mmax
               do n = 0, nmax
                  dnorm = 1._dp/(mscale(m)*nscale(n))
                  ccmn = cosmu(i,m)*cosnv(j,n)*dnorm
                  ssmn = sinmu(i,m)*sinnv(j,n)*dnorm
                  dm = m * bsupu(j,i)
                  dn = n * bsupv(j,i) * nfp
                  termsc = dm*ccmn - dn*ssmn
                  termcs =-dm*ssmn + dn*ccmn
                  if (n.eq.0 .or. n.eq.nmax) then      
                     mntot = mntot + 1
                     if (m .gt. 0) then
                        amatrix(ijtot,mntot) = termsc     !!only bsc != 0 for n=0, nmax1
                     else if (n .eq. 0) then
                        amatrix(ijtot,mntot) = bsupv(j,i) !!pedestal for m=0,n=0 mode, which should = 0
                     else 
                        amatrix(ijtot,mntot) = termcs     !!bcs(m=0,n=nmax)
                     end if
                  else if (m.eq.0 .or. m.eq.mmax) then
                     mntot = mntot + 1
                     amatrix(ijtot,mntot) = termcs        !!only bcs != 0 for m=0,mmax
                  else 
                     amatrix(ijtot,mntot+1) = termsc
                     amatrix(ijtot,mntot+2) = termcs
                     mntot = mntot + 2
                  end if   
               end do
            end do   
         end do
      end do   

      save_matrix = amatrix
      
      if (ijtot .ne. itotal .or. mntot .ne. itotal) then
         print *,' itotal = ', itotal,' ijtot = ', ijtot,
     1   ' mntot = ', mntot
         stop
      end if
      if (mmax+1 .ne. ntheta2) stop 'Error 1 in getbrho'
      if (nmax   .ne. nzeta/2) stop 'Error 2 in getbrho'

      call solver (amatrix, brhs, itotal)
      
!
!     CHECK SOLUTION FROM SOLVER
!
      
      ijtot = 0
      do i = 1, ntheta2
         do j = 1, nzeta
            lpior0 = ((i.eq.1 .or. i.eq.ntheta2) .and. (j.gt.nzeta/2+1))
            if (lpior0) cycle
            ijtot = ijtot + 1
            amn = sum(save_matrix(ijtot,:)*brhs(:))
            if (abs(frho(j,i) - amn) .gt. 1.e-8_dp*abs(amn))
     1      print *,' i = ',i,' j = ',j,' Original force = ',
     2      frho(j,i),' Final force = ', amn
         end do
      end do   
      
!
!     CONVERT TO BS*SIN(MU - NV) REPRESENTATION
!
      mntot = 0
      brhomn = 0
      do m = 0, mmax
         do n = 0, nmax
            if (n.eq.0 .or. n.eq.nmax) then     
               mntot = mntot + 1
               if (m .gt. 0) then
                  if (n .eq. 0) then
                     brhomn(m,n) = brhs(mntot)         !!bcs(m,0), m > 0
                  else   
                     brhomn(m,n) = p5*brhs(mntot)      !!bsc(m,nmax), m > 0
                     brhomn(m,-n)= p5*brhs(mntot)
                  end if   
               else if (n .ne. 0) then
                  brhomn(m,n) = -brhs(mntot)            !!bcs(0,nmax): needed for derivative
               end if   
            else if (m.eq.0 .or. m.eq.mmax) then
               mntot = mntot + 1
               if (m .eq. 0) then
                  brhomn(m,n)  = -brhs(mntot)          !!bcs(0,n) for n !=0, nmax
               else 
                  brhomn(m,n)  = -p5*brhs(mntot)       !!bcs(mmax,n)
                  brhomn(m,-n) =  p5*brhs(mntot)
               end if
            else
               brhomn(m,n)  = p5*(brhs(mntot+1) - brhs(mntot+2))
               brhomn(m,-n) = p5*(brhs(mntot+1) + brhs(mntot+2))
               mntot = mntot + 2
            end if
         end do
      end do      

      if (mntot .ne. ijtot) stop 'mntot != ijtot at end of getbrho'


      deallocate (amatrix, save_matrix, brhs)

      end subroutine getbrho
      

      subroutine mercier(gsqrt, bsq, bdotj, iotas, wint,
     1           r1, rt, rz, zt, zz, bsubu, vp, phips, pres, ns, nznt)
      use safe_open_mod
      use vmercier
      use vmec_input, ONLY: input_extension
      use vparams, only: one, zero, twopi, nmercier0
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: ns, nznt
      real(rprec), dimension(ns,nznt), intent(in) ::
     1     gsqrt, bsq
      real(rprec), dimension(ns,nznt), intent(inout) :: bdotj
      real(rprec), dimension(ns*nznt), intent(in) :: wint, bsubu
      real(rprec), dimension(ns,nznt,0:1), intent(in) ::
     1     r1, rt, rz, zt, zz
      real(rprec), dimension(ns), intent(in) ::
     1     iotas, vp, phips, pres
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: p5 = 0.5_dp, two = 2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ns1, i, imercier0, nmerc = nmercier0
      real(rprec) :: sign_jac, hs, sqs, denom
      real(rprec), dimension(:,:), allocatable :: 
     1     gpp, gsqrt_full, b2
      real(rprec), dimension(nznt) :: gtt, rtf, ztf,
     1   rzf, zzf, r1f, jdotb, ob2, b2i
      real(rprec), dimension(ns) :: vp_real, phip_real, 
     1   shear, vpp, presp, torcur, ip, sj, tpp, tjb, tbb, tjj
      character*120 mercier_file
      real(rprec), external :: dot_g
C-----------------------------------------------

      mercier_file = 'mercier.'//trim(input_extension)
      call safe_open (nmerc, imercier0, mercier_file, 'replace',
     1   'formatted')
      if (imercier0 .ne. 0) return


      allocate (gpp(ns,nznt), gsqrt_full(ns,nznt), b2(ns,nznt), stat=i)
      if (i .ne. 0) stop 'allocation error in Mercier'

!
!     SCALE VP, PHIPS TO REAL UNITS (VOLUME, TOROIDAL FLUX DERIVATIVES)
!     AND PUT GSQRT IN ABS UNITS (SIGNGS MAY BE NEGATIVE)      
!     NOTE: VP has (coming into this routine) the sign of the jacobian multiplied out
!          i.e., vp = signgs*<gsqrt>
!     THE SHEAR TERM MUST BE MULTIPLIED BY THE SIGN OF THE JACOBIAN
!     (OR A BETTER SOLUTION IS TO RETAIN THE JACOBIAN SIGN IN ALL TERMS, INCLUDING
!      VP, THAT DEPEND EXPLICITLY ON THE JACOBIAN. WE CHOOSE THIS LATTER METHOD...)
!
!     COMING INTO THIS ROUTINE, THE JACOBIAN(gsqrt) = 1./(grad-s . grad-theta X grad-zeta)
!     WE CONVERT THIS FROM grad-s to grad-phi DEPENDENCE BY DIVIDING gsqrt by PHIP_REAL
!
!     NOTE: WE ARE USING 0 < s < 1 AS THE FLUX VARIABLE, BEING CAREFUL
!     TO KEEP d(phi)/ds == PHIP_REAL FACTORS WHERE REQUIRED
!     THE V'' TERM IS d2V/d(PHI)**2, PHI IS REAL TOROIDAL FLUX
!
!     SHEAR = d(iota)/d(phi)   :  FULL MESH
!     VPP   = d(vp)/d(phi)     :  FULL MESH
!     PRESP = d(pres)/d(phi)   :  FULL MESH  (PRES IS REAL PRES*mu0)
!     IP    = d(Itor)/d(phi)   :  FULL MESH
!
!     ON ENTRY, BDOTJ = Jacobian * J*B  ON THE FULL RADIAL GRID
!               BSQ = 0.5*|B**2| + p IS ON THE HALF RADIAL GRID
!

      ns1 = ns - 1 
      if (ns1 .le. 0) return
      hs = one/ns1
      sign_jac = zero
      if (gsqrt(ns,1) .ne. zero)
     1    sign_jac = abs(gsqrt(ns,1))/gsqrt(ns,1)

      if (sign_jac .eq. zero) return
      phip_real = twopi * phips * sign_jac
!
!     NOTE: phip_real should be > 0 to get the correct physical sign of real-space gradients
!     For example, grad-p, grad-Ip, etc. However, with phip_real defined this way,
!     Mercier will be correct
!
      vp_real(2:ns) = sign_jac*(twopi*twopi)*vp(2:ns)/phip_real(2:ns)  !!dV/d(PHI) on half mesh
      
!
!     COMPUTE INTEGRATED TOROIDAL CURRENT
!
      do i = 2,ns
         torcur(i)=sign_jac*twopi*dot_g(nznt,bsubu(i),ns,wint(i),ns)
      end do

!
!     COMPUTE SURFACE AVERAGE VARIABLES ON FULL RADIAL MESH
!
      do i = 2,ns1
        phip_real(i) = p5*(phip_real(i+1) + phip_real(i))
        denom     = one/(hs*phip_real(i))
        shear(i)  = (iotas(i+1) - iotas(i))*denom       !!d(iota)/d(PHI)
        vpp(i)    = (vp_real(i+1) - vp_real(i))*denom   !!d(VP)/d(PHI)
        presp(i)  = (pres(i+1) - pres(i))*denom         !!d(p)/d(PHI)
        ip(i)     = (torcur(i+1) - torcur(i))*denom     !!d(Itor)/d(PHI)
      end do  
      
!
!     COMPUTE GPP == |grad-phi|**2 = PHIP**2*|grad-s|**2           (on full mesh)
!             GSQRT_FULL = JACOBIAN/PHIP == jacobian based on flux (on full mesh)
!

      do i = 2, ns1
        gsqrt_full(i,:) = p5*(gsqrt(i,:) + gsqrt(i+1,:))
        bdotj(i,:) = bdotj(i,:)/gsqrt_full(i,:)
        gsqrt_full(i,:) = gsqrt_full(i,:)/phip_real(i)
        sj(i) = hs*(i-1)
        sqs = sqrt(sj(i))
        rtf(:) = rt(i,:,0) + sqs*rt(i,:,1)
        ztf(:) = zt(i,:,0) + sqs*zt(i,:,1)
        gtt(:) = rtf(:)*rtf(:) + ztf(:)*ztf(:)
        rzf(:) = rz(i,:,0) + sqs*rz(i,:,1)
        zzf(:) = zz(i,:,0) + sqs*zz(i,:,1)
        r1f(:) = r1(i,:,0) + sqs*r1(i,:,1)
        gpp(i,:) = gsqrt_full(i,:)**2/(gtt(:)*r1f(:)**2 + 
     1             (rtf(:)*zzf(:) - rzf(:)*ztf(:))**2)     !!1/gpp
      end do
      
!
!     COMPUTE SURFACE AVERAGES OVER dS/|grad-PHI|**3 => |Jac| du dv / |grad-PHI|**2      
!     WHERE Jac = gsqrt/phip_real
!
      do i = 2,ns
        b2(i,:) = two*(bsq(i,:) - pres(i))
      end do
        
      do i = 2,ns1
        b2i(:) = p5*(b2(i+1,:) + b2(i,:))
        ob2(:) = gsqrt_full(i,:)/b2i(:)
        tpp(i) = DOT_G(nznt,ob2,1,wint(i),ns)                !<1/B**2>
        ob2(:) = b2i(:) * gsqrt_full(i,:) * gpp(i,:)
        tbb(i) = DOT_G(nznt,ob2,1,wint(i),ns)                !<b*b/|grad-phi|**3>
        jdotb(:) = bdotj(i,:) * gpp(i,:) * gsqrt_full(i,:)
        tjb(i) = DOT_G(nznt,jdotb,1,wint(i),ns)              !<j*b/|grad-phi|**3>
        jdotb(:) = jdotb(:) * bdotj(i,:) / b2i(:)
        tjj(i) = DOT_G(nznt,jdotb,1,wint(i),ns)              !<(j*b)2/b**2*|grad-phi|**3>
      end do
            
      deallocate (gpp, gsqrt_full, b2, stat=i)

!
!     REFERENCE: BAUER, BETANCOURT, GARABEDIAN, MHD Equilibrium and Stability of Stellarators
!     We break up the Omega-subs into a positive shear term (Dshear) and a net current term, Dcurr
!     Omega_subw == Dwell and Omega-subd == Dgeod (geodesic curvature, Pfirsch-Schluter term)
!
!     Include (eventually) Suydam for reference (cylindrical limit)
!
       
      write(nmerc,90)
 90   format(6x,'S',10x,'PHI',9x,'IOTA',8x,'SHEAR',7x,' VP ',8x,'WELL',
     1       8x,'ITOR',7x,'ITOR''',7x,'PRES',7x,'PRES''',/,120('-'))
      
      do i = 2,ns1
         sqs = p5*(vp_real(i) + vp_real(i+1))*sign_jac
         if (sqs .eq. zero) cycle
         write(nmerc,100) sj(i), hs*sum(phip_real(2:i)),
     1   p5*(iotas(i+1)+iotas(i)), shear(i)/sqs,
     2   sqs, -vpp(i)*sign_jac, 
     3   p5*(torcur(i) + torcur(i+1)), ip(i)/sqs,
     4   p5*(pres(i) + pres(i+1)), presp(i)/sqs
      end do

 100  format(1p10e12.4)      
      
      write(nmerc,190)
 190  format(/,6x,'S',8x,'DMerc',8x,'DShear',7x,'DCurr',7x,'DWell',
     1     7x,'Dgeod',/,100('-'))      
     
      do i = 2,ns1
         tpp(i) = (twopi*twopi)*tpp(i)
         tjb(i) = (twopi*twopi)*tjb(i)
         tbb(i) = (twopi*twopi)*tbb(i)
         tjj(i) = (twopi*twopi)*tjj(i)
         Dshear(i) = shear(i) * shear(i)/4
         Dcurr(i)  =-shear(i) * (tjb(i) - ip(i) *tbb(i))
         Dwell(i)  = presp(i) * (vpp(i) - presp(i) *tpp(i))*tbb(i) 
         Dgeod(i)  = tjb(i) *tjb(i)  - tbb(i) *tjj(i) 
         DMerc(i)  = Dshear(i) + Dcurr(i) + Dwell(i) + Dgeod(i)
         write(nmerc,100) sj(i), Dmerc(i), Dshear(i),
     1         Dcurr(i), Dwell(i), Dgeod(i)
      end do

      close (nmerc)
     
      end subroutine mercier


      subroutine residue(gcr, gcz, gcl, facmul)
      use vmec_main
      use vmec_params, only: meven, modd, ntmax, signgs
      use vsvd
      use xstuff
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,0:ntor,0:mpol1,ntmax) :: 
     1  gcr, gcz, gcl, facmul
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: n0 = 0, m0 = 0, m1 = 1
      real(rprec), parameter :: two = 2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: n, nsfix
      real(rprec) :: r1, r2
      real(rprec) :: tnorm, t1, t2, fac, bz1
C-----------------------------------------------
 
!
!       IMPOSE M=1 MODE CONSTRAINT TO MAKE THETA ANGLE
!       INVARIANT TO PHI-SHIFTS (AND THETA SHIFTS FOR ASYMMETRIC CASE)
!       ( ZCS = RSS, ZCC = RSC ARE THE CORRECT POLAR RELATIONS; HERE WE USE
!         THE SIMPLER RELATION ZCS = 0 INSTEAD)
!
 
      if (lthreed .and. fsqz.lt.c1pm8 .and. iter1.ne.1) then
         gcz(2:ns,1:ntor,m1,1) = zero
      endif

      if (lasym) then
      gcr(2:,:,m1,4) = gcr(2:,:,m1,4) + gcz(2:,:,m1,3)
      gcz(2:,:,m1,3) = 0
      do n = 0, ntor
         t1 = azd(ns,2) + bzd(ns,2) + (n*nfp)**2*crd(ns)
         t2 = ard(ns,2) + brd(ns,2) + (n*nfp)**2*crd(ns)
         t1 = t2/(t1 + t2)
         gcr(2:,n,m1,4) = t1*gcr(2:,n,m1,4) 
      end do

      n = max (3, mpol1)
      gcr(2:, :, n:,4) = 0              !!NEED THIS TO IMPROVE CONVERGENCE....
      gcz(2:, :, 4:,3) = 0

      end if
!
!       PUT FORCES INTO PHIFSAVE UNITS THAT PRECONDITIONERS, FNORM ARE IN
!
      if (phifac .eq. zero) then
         stop 'phifac = 0 in residue'
      else
         tnorm = phifsave/phifac           !put all forces into phifac=phifsave units
      end if   
 
 
      if (lrecon) then
!
!       MOVE R(n=0,m=0) TO SATISFY LIMITER OR AXIS POSITION
!       USE XC(NEQS2) TO STORE THIS PERTURBATION
!       TO SATISFY FORCE BALANCE AT JS=1, ADJUST PFAC IN RADFOR
!       ALSO, SCALE TOROIDAL FLUX AT EDGE TO MATCH MINOR RADIUS
 
        r1 = sum(gcr(:ns,n0,m0,1))
        if (r0scale .eq. zero) stop 'r0scale = 0'
        fsqsum0 = signgs*hs*r1/r0scale
        nsfix = 1                   !fix origin for reconstruction mode
        gcr = gcr * tnorm**2
        gcz = gcz * tnorm**2
        gcl = gcl * tnorm
        if (iopt_raxis.gt.0 .and. iresidue.eq.2 
     1     .and. fsq.lt.fopt_axis) iresidue = 3
        if (iresidue .lt. 3) gcr(nsfix,n0,m0,1) = zero
      else
!
!     ADJUST PHIEDGE
!
         if (imovephi .gt. 0) call movephi1 (gphifac)
      endif
      gc(neqs1) = gphifac
!
!       CONSTRUCT INVARIANT RESIDUALS
!
 
      call getfsq (gcr, gcz, fsqr, fsqz, fnorm, meven)
 
      r2 = sum(gcr(ns,:ntor,:mpol1,:)*gcr(ns,:ntor,:mpol1,:)
     1   +     gcz(ns,:ntor,:mpol1,:)*gcz(ns,:ntor,:mpol1,:))
      fedge = fnorm*r2
!
!       PERFORM PRECONDITIONING AND COMPUTE RESIDUES
!
      call scalfor (gcr, arm, brm, ard, brd, crd)
      call scalfor (gcz, azm, bzm, azd, bzd, crd)
!
!     APPLY ZCC = RSC CONSTRAINT FOR M=1 MODES
!
      if (lasym) gcz(2:,:,m1,3) = gcr(2:,:,m1,4)

      call getfsq (gcr, gcz, fsqr1, fsqz1, fnorm1, modd)
      fac = one/(one + (fsqr + fsqz))
      gcr = fac*gcr
      gcz = fac*gcz

!
!     CONSTRUCT INVARIANT AND PRECONDITIONED LAMBDA FORCES
!     NOTE: IF PHIP**2 USED IN BCOVAR FOR BU, BV COMPONENTS (RATHER THAN PHIP)
!     MUST MULTIPLY RBTOR IN BZ1 BY ADDITIONAL HS*SUM(ABS(PHIPS(2:NS))) FACTOR
!
      bz1 = two*hs/(rbtor*tnorm)**2
      fsql = bz1*sum(gcl*gcl)
      gcl = facmul*gcl
      fsql1 = hs*sum(gcl*gcl)
 
      if (fsqr.gt.1.e+4_dp .or. fsqz.gt.1.e+4_dp) then
         if (iter2.eq.iter1) gcl = 2.e-1_dp*gcl
         if (iter2.eq.1) irst = 4
      end if   
          
      end subroutine residue
      

      subroutine scalfor(gcx, axm, bxm, axd, bxd, cx)
      use vmec_main
      use vmec_params, only: jmin2, ntmax
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,0:ntor,0:mpol1,ntmax), 
     1  intent(inout) :: gcx
      real(rprec), dimension(ns + 1,2), intent(in) :: 
     1  axm, bxm, axd, bxd
      real(rprec), dimension(ns), intent(in) :: cx
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec), parameter :: ftol_edge = 1.e-9_dp
      integer :: m , mp, n, js, mmin, nmin, jmax, jmin4(0:mnsize)
      real(rprec), dimension(:,:,:), allocatable :: ax, bx, dx
      logical, save :: ledge
C-----------------------------------------------
      allocate (ax(ns,0:ntor,0:mpol1), bx(ns,0:ntor,0:mpol1),
     1    dx(ns,0:ntor,0:mpol1))

      jmax = ns
      if (ivac .lt. 1) jmax = ns1
!
!     ACCELERATE (IMPROVE) CONVERGENCE OF FREE BOUNDARY. THIS WAS ADDED
!     TO DEAL WITH CASES WHICH MIGHT OTHERWISE DIVERGE. BY DECREASING THE
!     FSQ TOLERANCE LEVEL WHERE THIS KICKS IN (FTOL_EDGE), THE USER CAN
!     TURN-OFF THIS FEATURE
!
      if (fsqr + fsqz .lt. ftol_edge) ledge = .true.
      if (iter2.lt.400 .or. ivac.lt.1) ledge = .false.

      do m = 0, mpol1
         mp = mod(m,2) + 1
         do n = 0, ntor
            do js = jmin2(m), jmax
               ax(js,n,m) = axm(js+1,mp) + bxm(js+1,mp)*m**2
               bx(js,n,m) = axm(js,mp) + bxm(js,mp)*m**2
               dx(js,n,m) = axd(js,mp) + bxd(js,mp)*m**2
     1                    + cx(js)*(n*nfp)**2
            end do

            if (m .eq. 1) dx(2,n,m) = dx(2,n,m) + bx(2,n,m)

            if (jmax .lt. ns) then
               dx(ns,n,m) = zero
            else if (m .le. 1) then
               dx(ns,n,m) = 1.05_dp*dx(ns,n,m)   !May need to raise 1.05 for 3D
            else
               dx(ns,n,m) = 1.10_dp*dx(ns,n,m)
            endif
         end do
      end do
      
!     FOR DATA MATCHING MODE (0 <= IRESIDUE < 3),
!     MAGNETIC AXIS IS FIXED SO JMIN3(0) => 2 FOR M=0,N=0
 
      jmin4 = jmin3       
      if (iresidue.ge.0 .and. iresidue.lt.3) jmin4(0) = 2
      
!     DIAGONALIZE (DX DOMINANT) AND REDUCE FORCE (DX ENHANCED) AT EDGE TO IMPROVE CONVERGENCE
 
      if (ledge) then
         bx(ns,:,:) = 0*bx(ns,:,:)
         dx(ns,:,:) = 3*dx(ns,:,:)
      end if
      
      call tridslv (ax, dx, bx, gcx, jmin4, 
     1      jmax, mnsize - 1, ns, ntmax)

      deallocate (ax, bx, dx)

      end subroutine scalfor

      subroutine symforce(ars, brs, crs, azs, bzs, czs, bls, cls, rcs, 
     1   zcs, ara, bra, cra, aza, bza, cza, bla, cla, rca, zca)
      use vmec_main
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns*nzeta,ntheta3,0:1), 
     1   intent(inout) :: ars, brs, crs, azs, bzs, czs, 
     2   bls, cls, rcs, zcs
      real(rprec), dimension(ns*nzeta,ntheta3,0:1), intent(out) :: 
     1   ara, bra, cra, aza, bza, cza, bla, cla, rca, zca
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: mpar, ir, i, jk, jka
C-----------------------------------------------
 
!
!       SYMMETRIZE FORCES ON RESTRICTED THETA INTERVAL (0 <= u <= pi)
!       SO COS,SIN INTEGRALS CAN BE PERFORMED. FOR EXAMPLE,
!
!       ARS(v,u) = .5*( ARS(v,u) + ARS(-v,-u) )     ! * cos(mu - nv)
!       ARA(v,u) = .5*( ARS(v,u) - ARS(-v,-u) )     ! * sin(mu - nv)
!
!
      do mpar = 0, 1
         do i = 1, ntheta2
            ir = ntheta1 + 2 - i                 !-theta
            if (i == 1) ir = 1
            do jk = 1, ns*nzeta
               jka = ireflect(jk)                !-zeta
               ara(jk,i,mpar) = cp5*(ars(jk,i,mpar)-ars(jka,ir,mpar))
               ars(jk,i,mpar) = cp5*(ars(jk,i,mpar)+ars(jka,ir,mpar))
               bra(jk,i,mpar) = cp5*(brs(jk,i,mpar)+brs(jka,ir,mpar))
               brs(jk,i,mpar) = cp5*(brs(jk,i,mpar)-brs(jka,ir,mpar))
               aza(jk,i,mpar) = cp5*(azs(jk,i,mpar)+azs(jka,ir,mpar))
               azs(jk,i,mpar) = cp5*(azs(jk,i,mpar)-azs(jka,ir,mpar))
               bza(jk,i,mpar) = cp5*(bzs(jk,i,mpar)-bzs(jka,ir,mpar))
               bzs(jk,i,mpar) = cp5*(bzs(jk,i,mpar)+bzs(jka,ir,mpar))
               bla(jk,i,mpar) = cp5*(bls(jk,i,mpar)-bls(jka,ir,mpar))
               bls(jk,i,mpar) = cp5*(bls(jk,i,mpar)+bls(jka,ir,mpar))
               rca(jk,i,mpar) = cp5*(rcs(jk,i,mpar)-rcs(jka,ir,mpar))
               rcs(jk,i,mpar) = cp5*(rcs(jk,i,mpar)+rcs(jka,ir,mpar))
               zca(jk,i,mpar) = cp5*(zcs(jk,i,mpar)+zcs(jka,ir,mpar))
               zcs(jk,i,mpar) = cp5*(zcs(jk,i,mpar)-zcs(jka,ir,mpar))
            end do
            if (lthreed) then
               do jk = 1, ns*nzeta
                  jka = ireflect(jk)
                  cra(jk,i,mpar) = cp5*(crs(jk,i,mpar)+crs(jka,ir,mpar))
                  crs(jk,i,mpar) = cp5*(crs(jk,i,mpar)-crs(jka,ir,mpar))
                  cza(jk,i,mpar) = cp5*(czs(jk,i,mpar)-czs(jka,ir,mpar))
                  czs(jk,i,mpar) = cp5*(czs(jk,i,mpar)+czs(jka,ir,mpar))
                  cla(jk,i,mpar) = cp5*(cls(jk,i,mpar)-cls(jka,ir,mpar))
                  cls(jk,i,mpar) = cp5*(cls(jk,i,mpar)+cls(jka,ir,mpar))
               end do
            endif
         end do
      end do
 
      end subroutine symforce

        
      subroutine symrzl(r1s, rus, rvs, z1s, zus, zvs, lus, lvs, rcons, 
     1   zcons, r1a, rua, rva, z1a, zua, zva, lua, lva, rcona, zcona)
      use vmec_main
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns*nzeta,ntheta3,0:1), intent(inout) :: 
     1   r1s, rus, rvs, z1s, zus, zvs, lus, lvs, rcons, zcons
      real(rprec), dimension(ns*nzeta,ntheta3,0:1), intent(in) :: 
     1   r1a, rua, rva, z1a, zua, zva, lua, lva, rcona, zcona
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: mpar, ir, i, jk, jka, n2
C-----------------------------------------------
 
!
!     FIRST SUM SYMMETRIC, ANTISYMMETRIC PIECES ON EXTENDED INTERVAL, THETA = [PI,2*PI]
!
      do mpar = 0, 1
         do i = 1 + ntheta2, ntheta1
            ir = ntheta1 + 2 - i                 !-theta
            do jk = 1, ns*nzeta
               jka = ireflect(jk)                !-zeta
               r1s(jk,i,mpar) = r1s(jka,ir,mpar) - r1a(jka,ir,mpar)
               rus(jk,i,mpar) =-rus(jka,ir,mpar) + rua(jka,ir,mpar)
               z1s(jk,i,mpar) =-z1s(jka,ir,mpar) + z1a(jka,ir,mpar)
               zus(jk,i,mpar) = zus(jka,ir,mpar) - zua(jka,ir,mpar)
               lus(jk,i,mpar) = lus(jka,ir,mpar) - lua(jka,ir,mpar)
               rcons(jk,i,mpar)= rcons(jka,ir,mpar)-rcona(jka,ir,mpar)
               zcons(jk,i,mpar)=-zcons(jka,ir,mpar)+zcona(jka,ir,mpar)
            end do
            if (lthreed) then
               do jk = 1, ns*nzeta
                  jka = ireflect(jk)
                  rvs(jk,i,mpar)=(-rvs(jka,ir,mpar))+rva(jka,ir,mpar)
                  zvs(jk,i,mpar) = zvs(jka,ir,mpar) - zva(jka,ir,mpar)
                  lvs(jk,i,mpar) = lvs(jka,ir,mpar) - lva(jka,ir,mpar)
               end do
            endif
         end do
!
!        NOW SUM SYMMETRIC, ANTISYMMETRIC PIECES FOR THETA = [0,PI]
!
         n2 = ntheta2
         r1s(:,:n2,mpar) = r1s(:,:n2,mpar) + r1a(:,:n2,mpar)
         rus(:,:n2,mpar) = rus(:,:n2,mpar) + rua(:,:n2,mpar)
         z1s(:,:n2,mpar) = z1s(:,:n2,mpar) + z1a(:,:n2,mpar)
         zus(:,:n2,mpar) = zus(:,:n2,mpar) + zua(:,:n2,mpar)
         lus(:,:n2,mpar) = lus(:,:n2,mpar) + lua(:,:n2,mpar)
         rcons(:,:n2,mpar) = rcons(:,:n2,mpar) + rcona(:,:n2,mpar)
         zcons(:,:n2,mpar) = zcons(:,:n2,mpar) + zcona(:,:n2,mpar)
         if (lthreed) then
            rvs(:,:n2,mpar) = rvs(:,:n2,mpar) + rva(:,:n2,mpar)
            zvs(:,:n2,mpar) = zvs(:,:n2,mpar) + zva(:,:n2,mpar)
            lvs(:,:n2,mpar) = lvs(:,:n2,mpar) + lva(:,:n2,mpar)
         endif
      end do
 
      end subroutine symrzl


      subroutine totzspa(rzl_array, r11, ru1, rv1, z11, zu1, zv1, lu1, 
     1   lv1, rcn1, zcn1)
      use vmec_main
      use vmec_params, only: jmin1, jlam, ntmax
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,0:ntor,0:mpol1,3*ntmax), 
     1   intent(inout) :: rzl_array
      real(rprec), dimension(ns*nzeta,ntheta3,0:1), intent(out) ::
     1   r11, ru1, rv1, z11, zu1, zv1, lu1, lv1, rcn1, zcn1
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: rmncs = 3, rmnsc = 4
      integer, parameter :: m0 = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: zmncc, zmnss, lmncc, lmnss
      integer :: m, n, mparity, k, js, i, jk, j1, l
      real(rprec), dimension(:,:), allocatable :: work1
      real(rprec) :: cosmux, sinmux
C-----------------------------------------------

      allocate (work1(ns*nzeta,12))

      zmncc = rmncs + ntmax
      zmnss = rmnsc + ntmax
      lmncc = rmncs + 2*ntmax
      lmnss = rmnsc + 2*ntmax

!
!                INITIALIZATION BLOCK
!
      r11 = 0;  ru1 = 0;  rv1 = 0;  z11 = 0;  zu1 = 0
      zv1 = 0;  lu1 = 0;  lv1 = 0;  rcn1 = 0; zcn1 = 0
  
      if (jlam(m0) .gt. 1) then
      rzl_array(1,1:ntor,m0,lmncc) = 2*rzl_array(2,1:ntor,m0,lmncc)
     1                             -   rzl_array(3,1:ntor,m0,lmncc)
      end if

      do m = 0, mpol1
         mparity = mod(m,2)
         work1 = 0
         j1 = jmin1(m)
         do n = 0, ntor
            do k = 1, nzeta
               l = ns*(k-1)
               do js = j1, ns
                  work1(js+l,1) = work1(js+l,1) + 
     1               rzl_array(js,n,m,rmnsc)*cosnv(k,n)
                  work1(js+l,6) = work1(js+l,6) + 
     1               rzl_array(js,n,m,zmncc)*cosnv(k,n)
                  work1(js+l,10) = work1(js+l,10) + 
     1               rzl_array(js,n,m,lmncc)*cosnv(k,n)
               end do
               if (lthreed) then
                  do js = j1, ns
                     work1(js+l,2) = work1(js+l,2) + 
     1                  rzl_array(js,n,m,rmncs)*sinnv(k,n)
                     work1(js+l,3) = work1(js+l,3) + 
     1                  rzl_array(js,n,m,rmnsc)*sinnvn(k,n)
                     work1(js+l,4) = work1(js+l,4) + 
     1                  rzl_array(js,n,m,rmncs)*cosnvn(k,n)
                     work1(js+l,5) = work1(js+l,5) + 
     1                  rzl_array(js,n,m,zmnss)*sinnv(k,n)
                     work1(js+l,7) = work1(js+l,7) + 
     1                  rzl_array(js,n,m,zmnss)*cosnvn(k,n)
                     work1(js+l,8) = work1(js+l,8) + 
     1                  rzl_array(js,n,m,zmncc)*sinnvn(k,n)
                     work1(js+l,9) = work1(js+l,9) + 
     1                  rzl_array(js,n,m,lmnss)*sinnv(k,n)
                     work1(js+l,11) = work1(js+l,11) + 
     1                  rzl_array(js,n,m,lmnss)*cosnvn(k,n)
                     work1(js+l,12) = work1(js+l,12) + 
     1                  rzl_array(js,n,m,lmncc)*sinnvn(k,n)
                  end do
               endif
            end do
         end do
 
!
!        INVERSE TRANSFORM IN M-THETA
!
         do i = 1, ntheta2
            cosmux = xmpq(m,1)*cosmu(i,m)
            sinmux = xmpq(m,1)*sinmu(i,m)
            r11(:,i,mparity) = r11(:,i,mparity) + work1(:,1)*
     1            sinmu(i,m)
            ru1(:,i,mparity) = ru1(:,i,mparity) + work1(:,1)*
     1            cosmum(i,m)
            z11(:,i,mparity) = z11(:,i,mparity) + work1(:,6)*
     1            cosmu(i,m)
            zu1(:,i,mparity) = zu1(:,i,mparity) + work1(:,6)*
     1            sinmum(i,m)
            lu1(:,i,mparity) = lu1(:,i,mparity) + work1(:,10)*
     1            sinmum(i,m)
            rcn1(:,i,mparity) = rcn1(:,i,mparity) + work1(:,1)*
     1            sinmux
            zcn1(:,i,mparity) = zcn1(:,i,mparity) + work1(:,6)*
     1            cosmux
 
            if (lthreed) then
               r11(:,i,mparity) = r11(:,i,mparity) + work1(:,2)*
     1               cosmu(i,m)
               ru1(:,i,mparity) = ru1(:,i,mparity) + work1(:,2)*
     1               sinmum(i,m)
               z11(:,i,mparity) = z11(:,i,mparity) + work1(:,5)*
     1               sinmu(i,m)
               zu1(:,i,mparity) = zu1(:,i,mparity) + work1(:,5)*
     1               cosmum(i,m)
               lu1(:,i,mparity) = lu1(:,i,mparity) + work1(:,9)*
     1               cosmum(i,m)
               rcn1(:,i,mparity) = rcn1(:,i,mparity) + work1(:,2)*
     1               cosmux
               zcn1(:,i,mparity) = zcn1(:,i,mparity) + work1(:,5)*
     1               sinmux
               rv1(:,i,mparity) = rv1(:,i,mparity) + work1(:,3)*
     1               sinmu(i,m) + work1(:,4)*cosmu(i,m)
               zv1(:,i,mparity) = zv1(:,i,mparity) + work1(:,7)*
     1               sinmu(i,m) + work1(:,8)*cosmu(i,m)
               lv1(:,i,mparity) = lv1(:,i,mparity) - (work1(:,11)*
     1               sinmu(i,m)+work1(:,12)*cosmu(i,m))
            endif
 
         end do
      end do
 
      deallocate (work1)

      end subroutine totzspa


      subroutine totzsps(rzl_array, r11, ru1, rv1, z11, zu1, zv1, lu1, 
     1   lv1, rcn1, zcn1)
      use vmec_main
      use vmec_params, only: jmin1, jlam, ntmax
      use vmec_persistent
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,0:ntor,0:mpol1,3*ntmax),
     1   intent(inout) :: rzl_array
      real(rprec), dimension(ns*nzeta,ntheta3,0:1),
     1   intent(out) :: r11, ru1, 
     1   rv1, z11, zu1,  zv1, lu1, lv1, rcn1, zcn1
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: rmncc = 1, rmnss = 2
      integer, parameter :: m0 = 0, m1 = 1, n0 = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: zmncs, zmnsc, lmncs, lmnsc
      integer :: n, m, mparity, k, js, i, j1, l
      real(rprec), dimension(:,:), allocatable :: work1
      real(rprec) :: cosmux, sinmux
C-----------------------------------------------

      zmncs = rmncc + ntmax
      zmnsc = rmnss + ntmax
      lmncs = rmncc + 2*ntmax
      lmnsc = rmnss + 2*ntmax

      allocate (work1(ns*nzeta,12))

      r11 = 0;  ru1 = 0;  rv1 = 0;  z11 = 0; zu1 = 0
      zv1 = 0;  rcn1 = 0; zcn1 = 0
      lu1(:,:,0) = 1;  lu1(:,:,1) = 0;  lv1 = 0

!
!     EXTRAPOLATION AT JS=1 FOR M=ODD MODES AND M=0 MODES OF LAMBDA
!
      rzl_array(1,:,m1,:) = 2*rzl_array(2,:,m1,:)
     1                    -   rzl_array(3,:,m1,:)

      if (jlam(m0) .gt. 1) then
      rzl_array(1,1:ntor,m0,lmncs) = 2*rzl_array(2,1:ntor,m0,lmncs)
     1                               - rzl_array(3,1:ntor,m0,lmncs)
      end if

!
!     COMPUTE R, Z, AND LAMBDA IN REAL SPACE
!     BEGIN INVERSE TRANSFORM IN N-ZETA
!     NOTE: LU = d(Lam)/du, LV = -d(Lam)/dv    
!
      do m = 0, mpol1
         mparity = mod(m,2)
         work1 = 0
         j1 = jmin1(m)
         do n = 0, ntor
            do k = 1, nzeta
               l = ns*(k-1)
               do js = j1,ns
               work1(js+l,1) = work1(js+l,1) + 
     1            rzl_array(js,n,m,rmncc)*cosnv(k,n)
               work1(js+l,6) = work1(js+l,6) + 
     1            rzl_array(js,n,m,zmnsc)*cosnv(k,n)
               work1(js+l,10) = work1(js+l,10) +
     1            rzl_array(js,n,m,lmnsc)*cosnv(k,n)
               end do
               if (lthreed) then
                  do js = j1,ns
                     work1(js+l,2) = work1(js+l,2) + 
     1                  rzl_array(js,n,m,rmnss)*sinnv(k,n)
                     work1(js+l,4) = work1(js+l,4) + 
     1                  rzl_array(js,n,m,rmnss)*cosnvn(k,n)
                     work1(js+l,9) = work1(js+l,9) + 
     1                  rzl_array(js,n,m,lmncs)*sinnv(k,n)
                     work1(js+l,11) = work1(js+l,11) + 
     1                  rzl_array(js,n,m,lmncs)*cosnvn(k,n)
                  end do
                  do js = j1,ns
                     work1(js+l,5) = work1(js+l,5) + 
     1                  rzl_array(js,n,m,zmncs)*sinnv(k,n)
                     work1(js+l,7) = work1(js+l,7) + 
     1                  rzl_array(js,n,m,zmncs)*cosnvn(k,n)
                     work1(js+l,3) = work1(js+l,3) + 
     1                  rzl_array(js,n,m,rmncc)*sinnvn(k,n)
                     work1(js+l,8) = work1(js+l,8) + 
     1                  rzl_array(js,n,m,zmnsc)*sinnvn(k,n)
                     work1(js+l,12) = work1(js+l,12) + 
     1                  rzl_array(js,n,m,lmnsc)*sinnvn(k,n)
                  end do
               endif
            end do
         end do
!
!        INVERSE TRANSFORM IN M-THETA
!
         do i = 1, ntheta2
            cosmux = xmpq(m,1)*cosmu(i,m)
            sinmux = xmpq(m,1)*sinmu(i,m)
            r11(:,i,mparity) = r11(:,i,mparity) + work1(:,1)*
     1            cosmu(i,m)
            ru1(:,i,mparity) = ru1(:,i,mparity) + work1(:,1)*
     1            sinmum(i,m)
            rcn1(:,i,mparity) = rcn1(:,i,mparity) + work1(:,1)*
     1            cosmux
            z11(:,i,mparity) = z11(:,i,mparity) + work1(:,6)*
     1            sinmu(i,m)
            zu1(:,i,mparity) = zu1(:,i,mparity) + work1(:,6)*
     1            cosmum(i,m)
            zcn1(:,i,mparity) = zcn1(:,i,mparity) + work1(:,6)*
     1            sinmux
            lu1(:,i,mparity) = lu1(:,i,mparity) + work1(:,10)*
     1            cosmum(i,m)
 
            if (lthreed) then
               r11(:,i,mparity) = r11(:,i,mparity) + work1(:,2)*
     1               sinmu(i,m)
               ru1(:,i,mparity) = ru1(:,i,mparity) + work1(:,2)*
     1               cosmum(i,m)
               rcn1(:,i,mparity) = rcn1(:,i,mparity) + work1(:,2)*
     1               sinmux
               z11(:,i,mparity) = z11(:,i,mparity) + work1(:,5)*
     1               cosmu(i,m)
               zu1(:,i,mparity) = zu1(:,i,mparity) + work1(:,5)*
     1               sinmum(i,m)
               zcn1(:,i,mparity) = zcn1(:,i,mparity) + work1(:,5)*
     1               cosmux
               lu1(:,i,mparity) = lu1(:,i,mparity) + work1(:,9)*
     1               sinmum(i,m)
               rv1(:,i,mparity) = rv1(:,i,mparity) + work1(:,3)*
     1               cosmu(i,m) + work1(:,4)*sinmu(i,m)
               zv1(:,i,mparity) = zv1(:,i,mparity) + work1(:,7)*
     1               cosmu(i,m) + work1(:,8)*sinmu(i,m)
               lv1(:,i,mparity) = lv1(:,i,mparity) - (work1(:,11)*
     1               cosmu(i,m) + work1(:,12)*sinmu(i,m))
            endif
 
         end do
      end do
 
      deallocate (work1)

      if (rzl_array(1,n0,m1,rmncc) .eq. zero)
     1   stop 'r01(0) = 0 in totzsps'      
!     dkappa = rzl_array(1,n0,m1,zmnsc)/rzl_array(1,n0,m1,rmncc)
!     r01 = rzl_array(ns,n0,m1,rmncc)

      end subroutine totzsps


      subroutine wrout(bsq, gsqrt, bsubu, bsubv, bsubs, bsupv, bsupu,
     1   izeta, ier_flag)
      use vmec_main
      use vmec_params, only: mscale, nscale, signgs, version_
      use vmercier
      use vmec_persistent
      use vsvd
      use vspline
      use xstuff
      use vmec_io
      use realspace
      use safe_open_mod
      use vacmod, ONLY: potvac,xmpot,xnpot,mnpd   !added for diagno, J.Geiger
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: ier_flag
      real(rprec), dimension(ns,nznt), intent(in) ::
     1   bsq, gsqrt, bsubu, bsubv, bsubs, bsupv, bsupu, izeta
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: two = 2, p5 = 0.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: iasym, ireconstruct, i, j, js, lj, mn, lk, isgn,
     1   m, nmin0, n, k, iwout0, l, ntotal, n1, nwout, nfort
      real(rprec), dimension(nznt) :: bmod
      real(rprec), dimension(mnmax) :: rmnc, zmns, lmns,
     1   rmns, zmnc, lmnc, bmodmn, bmodmn1, lmns_half, lmnc_half
      real(rprec) :: dmult, gmn, bmn, currvmn, bsubumn, bsubvmn,
     1   bsubsmn, bsupumn, bsupvmn, tcosi, tsini, presfactor, zeta,
     2   lmnsh, lmnch
!     4  ,potsin, potcos                            !added for diagno, J.Geiger
      character*120 :: wout_file, fort_file
C-----------------------------------------------

      rmns = 0
      zmnc = 0
      lmnc = 0
!
!     THIS SUBROUTINE CREATES THE TEXT FILE WOUT WHICH
!     CONTAINS THE CYLINDRICAL COORDINATE SPECTRAL
!     COEFFICIENTS RMN,ZMN (FULL MESH), LMN (HALF MESH - CONVERTED FROM
!              INTERNAL FULL MESH REPRESENTATION)
!
!     IZETA (FULL), BSQ, BSUPU,V, BSUBU,V, GSQRT (HALF)
!

      wout_file = 'wout.'//input_extension
      nwout = nwout0
      call safe_open(nwout, iwout0, wout_file, 'replace', 'formatted')
      if (iwout0 .ne. 0) stop 'Error writing WOUT file in VMEC WROUT'

      if (loldout) then
         nfort = nfort8
         fort_file = 'fort.8.'//input_extension
         call safe_open(nfort,iwout0,fort_file,'replace','formatted')
         if (iwout0 .ne. 0) then
             print *,'Problem opening fort.8 file'
             print *,'iwout0=',iwout0, 'for ',input_extension
             stop 'Error writing fort.8. file in VMEC WROUT '
         endif

         nfort = nfort18
         fort_file = 'fort.18.'//input_extension
         call safe_open(nfort,iwout0,fort_file,'replace','formatted')
         if (iwout0 .ne. 0) 
     1       stop 'Error writing fort.18. file in VMEC WROUT'
      endif

      if (lasym) then
         iasym = 1                                  ! asymmetric mode
      else
         iasym = 0
      end if
      if (lrecon) then
         ireconstruct = 1
      else
         itse = 0
         imse2 = 0
         ireconstruct = 0
      end if

!
!     Insert version information into wout file. This will be parsed in
!     read_wout_file to return the real value version_ to check the version number.
!
      write (nwout, '(a15,a)') 'VMEC VERSION = ', version_
      write (nwout, *) wb, wp, gamma, pfac,
     1  rmax_surf, rmin_surf, zmax_surf

      write (nwout, *) nfp, ns, mpol, ntor, mnmax, itfsq, niter,
     1  iasym, ireconstruct, ier_flag

      write (nwout, *) imse2 - 1, itse, nbsets, nobd, nextcur,
     1   nstore_seq
      if (nbsets .gt. 0) write (nwout, *) (nbfld(i),i=1,nbsets)
      write (nwout, '(a)') mgrid_file

!-----------------------------------------------
!     Modification to obtain old fort.8 file
!-----------------------------------------------
      if (loldout) then
        write (nfort8, 710) gamma, nfp, ns, mpol, ntor, mnmax, itfsq,
     1         niter/100+1
!       the variables helsym (4th) and kpres (last) were set to 0.
        write (nfort18, 40)voli,gamma,1.0/nfp,0.,mnmax,ns,mpol,ntor,
     1              ntor+1,1,itfsq,niter/100+1,0
      end if
 40   format(1x,1p,e22.12,1x,3e12.5,8i4,i1)
 710  format(f10.3,7i6)

      if (ier_flag.ne.0 .and. ier_flag.ne.4) goto 1000

      pres(1) = pres(2)
      ntotal = 0

      do js = 1, ns
         lj = (js - 1)*mnmax

         call convert (rmnc, zmns, lmns, rmns, zmnc, lmnc, xc, js)

         mn = 0
         bmod = sqrt(two*abs(bsq(js,:nznt)-pres(js)))
         do m = 0, mpol1
            nmin0 = -ntor
            if (m .eq. 0) nmin0 = 0
            do n = nmin0, ntor
               mn = mn + 1
               ntotal = max(mn,ntotal)
               dmult = two/(mscale(m)*nscale(abs(n)))
               if (m .eq. 0 .and. n .eq. 0) dmult = p5*dmult
               n1 = abs(n)
               isgn = sign(1, n)
               gmn = 0
               bmn = 0
               currvmn = 0
               bsubumn = 0
               bsubvmn = 0
               bsubsmn = 0
               bsupumn = 0
               bsupvmn = 0
               do j = 1, ntheta2
                  do k = 1, nzeta
                     lk = k + nzeta*(j - 1)
                     tcosi = dmult*(cosmui(j,m)*cosnv(k,n1) +
     1                         isgn*sinmui(j,m)*sinnv(k,n1))
                     tsini = dmult*(sinmui(j,m)*cosnv(k,n1) -
     1                         isgn*cosmui(j,m)*sinnv(k,n1))
                     bmn = bmn + tcosi*bmod(lk)
                     gmn = gmn + tcosi*gsqrt(js,lk)
                     currvmn = currvmn + tcosi*izeta(js,lk)
                     bsubumn = bsubumn + tcosi*bsubu(js,lk)
                     bsubvmn = bsubvmn + tcosi*bsubv(js,lk)
                     bsubsmn = bsubsmn + tsini*bsubs(js,lk)
                     bsupumn = bsupumn + tcosi*bsupu(js,lk)
                     bsupvmn = bsupvmn + tcosi*bsupv(js,lk)
                  end do
               end do
               if (js .eq. ns/2) bmodmn(mn) = bmn
               if (js .eq. ns) bmodmn1(mn) = bmn
               if(js .eq. 1) then
                  if (m .ne. nint(xm(mn))) stop 'm != nint(xm)'
                  if (n*nfp .ne. nint(xn(mn))) stop ' n = nint(xn)'
                  write (nwout, *) m, n*nfp
                  raxis(0:ntor,1) = rmnc(1:ntor+1)
                  zaxis(0:ntor,1) = zmns(1:ntor+1)
!                 ZERO HALF-MESH QUANTITIES
                  bsubumn = 0
                  bsubvmn = 0
                  bsubsmn = 0
                  bsupumn = 0
                  bsupvmn = 0
                  bmn = 0
                  gmn = 0
               end if
!
!       PUT LMNS, LMNC ONTO HALF-MESH FOR BACKWARDS CONSISTENCY
!
               if (js .eq. 1) then
                  lmnsh = 0
                  lmnch = 0
               else if (mod(m,2) .eq. 0) then
                  lmnsh = cp5*(lmns(mn) + lmns_half(mn))
                  lmnch = cp5*(lmnc(mn) + lmnc_half(mn))
               else
                  if (js.eq.2 .and. m.eq.1) then
                     lmns_half(mn) = lmns(mn)
                     lmnc_half(mn) = lmnc(mn)
                  end if
                  lmnsh = cp5*(sm(js)*lmns(mn) + sp(js-1)*lmns_half(mn))
                  lmnch = cp5*(sm(js)*lmnc(mn) + sp(js-1)*lmnc_half(mn))
               end if

               write (nwout, *) rmnc(mn), zmns(mn), lmnsh,
     1            bmn, gmn, bsubumn,
     2            bsubvmn, bsubsmn, bsupumn, bsupvmn, currvmn
               if (lasym) then
                  write (nwout, *) rmns(mn), zmnc(mn), lmnch
               endif

!-----------------------------------------------
!   Modification to obtain old fort.8 file
!-----------------------------------------------
               if (loldout) then
                 if (js .eq. 1)
     1              write (nfort8, 721) nint(xm(mn)),nint(xn(mn))
                 write (nfort18,50)  xm(mn),xn(mn),rmnc(mn),zmns(mn),gmn
                 write (nfort8, 731) rmnc(mn), zmns(mn), lmns(mn), bmn,
     1            gmn, bsubumn, bsubvmn, bsubsmn, bsupumn, bsupvmn
               end if
            end do
         end do
         lmns_half = lmns        !!Store previous full point values for averaging
         lmnc_half = lmnc
      end do
 50   format(1x,1p,5e14.6)
 721  format(2i10)
 731  format(5e20.13)

!
!     HALF-MESH QUANTITIES (except phi, jcuru, jcurv which are FULL MESH - computed in eqfor)
!     NOTE: jcuru are local current densities, NOT integrated in s
!     NOTE: In version <= 6.00, mass, press are written out in INTERNAL units
!     and should be multiplied by 1/mu0 to transform to pascals. In version > 6.00,
!     the pressure, mass are in correct (physical) units
!
      write (nwout, *) (iotas(js), mass(js)/dmu0, pres(js)/dmu0,
     1   beta_vol(js), phip(js), buco(js), bvco(js), phi(js), vp(js),
     2   overr(js), jcuru(js)/dmu0, jcurv(js)/dmu0, specw(js),js=2,ns)
!-----------------------------------------------
!   Modification to obtain old fort.8 file
!-----------------------------------------------
      if (loldout) then
         write (nfort8, 731) (iotas(js),mass(js),pres(js),phips(js),
     1          buco(js),bvco(js),phi(js),vp(js),jcuru(js)/dmu0,
     2          jcurv(js)/dmu0,specw(js),js=2,ns)
        write (nfort18,50)(iotas(js),mass(js),pres(js),-phips(js),
     1          vp(js),js=1,ns)
      end if
!-----------------------------------------------

      write (nwout, *) aspect, betatot, betapol, betator, betaxis, b0

!-----------------------------------------------
!     New output added to version 6.20
!-----------------------------------------------
      write (nwout, *) nint(signgs)
      write (nwout, '(a)') input_extension
      write (nwout, *) IonLarmor, VolAvgB, rbtor0, rbtor, ctor/dmu0,
     1  Aminor_p, Rmajor_p, volume_plasma
!-----------------------------------------------
!     MERCIER CRITERION
!-----------------------------------------------
      write (nwout, *) (Dmerc(js), Dshear(js), Dwell(js), Dcurr(js),
     1       Dgeod(js), equif(js), js=2,ns-1)

      if (nextcur .gt. 0) then
         write (nwout, *) (extcur(i),i=1,nextcur)
         write (nwout, *) (curlabel(i),i=1,nextcur)
      endif

!-----------------------------------------------
!     NOTE: jdotb is in units of A (1/mu0 incorporated in jxbforce...)
!     prior to version 6.00, this was output in internal VMEC units...
!-----------------------------------------------
      write (nwout, *) (fsqt(i),wdot(i),i=1,nstore_seq)
      write (nwout, *) (jdotb(js),bdotgradv(js),js=1,ns)

!-----------------------------------------------
!   Modification to obtain old fort.8 file
!-----------------------------------------------
      if (loldout) write (nfort8, 731) (fsqt(i),wdot(i),i=1,100)

!-----------------------------------------------
!   Write diagno file (J.Geiger)
!-----------------------------------------------
      IF(ldiagno)THEN         !added for diagno, J.Geiger start
         IF(lfreeb .and. .not.lasym)THEN
            nfort = nfort21
            fort_file = 'diagno_in.'//input_extension
            call safe_open(nfort,iwout0,fort_file,'replace',
     1                     'formatted')
            if (iwout0 .ne. 0) 
     1          stop 'Error writing diagno_in. file in VMEC WROUT'

            write(nfort21,'(a)') "vmec2000"
            write(nfort21,*) "nfp  mpol  ntor"
            write(nfort21,*) nfp, mpol, ntor

            CALL convert (rmnc, zmns, lmns, rmns, zmnc, lmnc, xc, ns)
            write(nfort21,*) "rmnc"
            write(nfort21,'(1p,e24.16)') (rmnc(mn),mn=1,mnmax)
            write(nfort21,*) "zmns"
            write(nfort21,'(1p,e24.16)') (zmns(mn),mn=1,mnmax)

            write(nfort21,*) "potsin"
            DO i = 1, mnpd
               write(nfort21,'(1p,e24.16)') potvac(i)
            END DO
            write(nfort21,*) "phiedge"
            write(nfort21,*) phiedge
            write(nfort21,*) "nextcur"
            write(nfort21,*) nextcur
            write(nfort21,*) "external currents"
            write(nfort21,*) extcur(1:nextcur)

            write(nfort21,*) "plasma current"
            write(nfort21,*) ctor
            CALL convert (rmnc, zmns, lmns, rmns, zmnc, lmnc, xc,1)
            write(nfort21,*) "plasma current filament fc R"
            write(nfort21,*) rmnc(1:ntor+1)
            write(nfort21,*) "plasma current filament fc z"
            write(nfort21,*) zmns(1:ntor+1)

            close(unit=nfort21)
         ELSE
            write(6,*)"Diagno-file request not completed!"
            write(6,*)"VMEC2000 not running in free-boundary mode!"
            write(6,*)"-or- LASYM = .true. !"
            write(6,*)"Check mgrid-file and namelist!"
         ENDIF
      ENDIF                   !added for diagno, J.Geiger end

!-----------------------------------------------
!     DATA AND MSE FITS
!-----------------------------------------------
      if (.not.lrecon) goto 900

      if (imse2 - 1.gt.2 .or. itse.gt.0) then
         write (nwout, *) tswgt, msewgt
         call smoothdata(nwout)

!       These knot values are on sqrt(s) grid
         presfactor = dmu0*pthommax             !!*pfac moved to getthom
         write (nwout, *) isnodes, (sknots(i),ystark(i),y2stark(i),i=
     1      1,isnodes)
         write (nwout, *) ipnodes, (pknots(i),presfactor*ythom(i),
     1      presfactor*y2thom(i),i=1,ipnodes)
         write (nwout, *)(datamse(i),rmid(i),qmid(i),shear(i),
     1      presmid(i),alfa(i),curmid(i),i=1,2*ns-1)
         write (nwout, *)(rstark(i),datastark(i),qmeas(i),i=1,imse)
         write (nwout, *)(rthom(i),datathom(i),i=1,itse)
      endif
      if (nobd .gt. 0) then
         write (nwout, *) (dsiext(i),plflux(i),dsiobt(i),i=1,nobd)
         write (nwout, *) flmwgt
      endif
      if (nbfldn .gt. 0) then
         do n = 1, nbsets
            write (nwout, *) (bcoil(i,n),plbfld(i,n),bbc(i,n),
     1         i=1,nbfld(n))
         end do
         write (nwout, *) bcwgt
      endif

      write (nwout, *) phidiam, delphid
!
!     Write Limiter & Prout plotting specs
!
      write (nwout, *) nsets, nparts, nlim
      write (nwout, *) (nsetsn(i),i=1,nsets)
      write (nwout, *) (((pfcspec(i,j,k),i=1,nparts),j=1,nsetsn(k)),
     1   k=1,nsets)
      write (nwout, *) (limitr(i), i=1,nlim)
      write (nwout, *) ((rlim(i,j),zlim(i,j),i=1,limitr(j)),j=1,nlim)
      write (nwout, *) nrgrid, nzgrid
      write (nwout, *) tokid
      write (nwout, *) rx1, rx2, zy1, zy2, condif
      write (nwout, *) imatch_phiedge
      if (imatch_phiedge .eq. 2) call wroutlim(nwout)

 900  continue
!-----------------------------------------------
!     FREE BOUNDARY DATA
!-----------------------------------------------
      if (ier_flag .eq. 0) call freeb_data(rmnc, zmns, bmodmn, bmodmn1)


 1000 continue

      close (nwout)
      if( loldout) then
         close(nfort8)
         close(nfort18)
      endif

      end subroutine wrout
            
      
      subroutine freeb_data (rmnc, zmns, bmodmn, bmodmn1)
      use vmec_main
      use vacmod
      use realspace, only: r1, z1
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(mnmax) :: rmnc, zmns, bmodmn, bmodmn1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: iprint, nzskip, l, k, lk, mnlim, mnres, mn,
     1           mn0, n, nedge, nedge0 = 99, iu, iv, nl, lkr
      real(rprec), allocatable, dimension(:) :: rb, phib, zb
C-----------------------------------------------
!
!     WRITE OUT EDGE VALUES OF FIELDS TO FORT.NEDGE0 (INCLUDE REFLECTED POINT)
!
!     NOTE: BR, BPHI, BZ WERE COMPUTED IN BSS, CALLED FROM EQFOR
!
      if (.not.lfreeb .and. .not.ledge_dump) return

      allocate (rb(2*nznt), phib(2*nznt), zb(2*nznt), stat=l)
      if (l .ne. 0) stop 'allocation error in freeb_data'

      nedge = 0
      lkr = nznt
      do iv = 1,nzeta
         zeta = (twopi*(iv-1))/(nzeta*nfp)
         do iu = 1,ntheta3
            lk = iv + nzeta*(iu-1)
            nl = ns*lk
            nedge = nedge+1
            rb(lk)   = r1(nl,0) + r1(nl,1)
            phib(lk) = zeta
            zb(lk)   = z1(nl,0) + z1(nl,1)
!
!      INCLUDE -u,-v POINTS HERE BY STELLARATOR SYMMETRY
!
            if (.not.lasym .and. (iu.ne.1 .and. iu.ne.ntheta2)) then
              lkr = lkr + 1
              nedge = nedge+1
              rb(lkr)   = rb(lk)
              phib(lkr) =-phib(lk)
              zb(lkr)   =-zb(lk)  
              bredge(lkr) = -bredge(lk)
              bpedge(lkr) =  bpedge(lk)         
              bzedge(lkr) =  bzedge(lk)
            endif  
         end do
      end do  

      if (ledge_dump) then
        write(NEDGE0,*) 'INPUT FILE = ',arg1
        write(NEDGE0,*) 'NEDGE = ',nedge
        write(NEDGE0,*) 'RB = ',  (rb(i), i=1,nedge)
        write(NEDGE0,*) 'PHIB = ',(phib(i), i=1,nedge)
        write(NEDGE0,*) 'ZB = ',  (zb(i), i=1,nedge)
        write(NEDGE0,*) 'BREDGE = ', (bredge(i), i=1,nedge)
        write(NEDGE0,*) 'BPEDGE = ', (bpedge(i), i=1,nedge)
        write(NEDGE0,*) 'BZEDGE = ', (bzedge(i), i=1,nedge)
      end if

!
!     WRITE OUT (TO THREED1 FILE) VACUUM INFORMATION
!

      if (.not.lfreeb) then
         deallocate (rb, phib, zb, stat=l)
         return
      end if
      
      do iprint = 1, 2
         if (iprint .eq. 1) write (nthreed, 750)
         if (iprint .eq. 2) write (nthreed, 760)
         nzskip = 1 + nzeta/6
         do l = 1, nzeta, nzskip
            zeta = (360.0_dp*(l - 1))/nzeta
            if (iprint .eq. 1) then
               do k = 1, ntheta2
                  lk = l + nzeta*(k - 1)
                  write (nthreed, 770) zeta, rb(lk),
     1            zb(lk), (bsqsav(lk,n),n=1,3), bsqvac(lk)
               end do
            else
               do k = 1, ntheta2
                  lk = l + nzeta*(k - 1)
                  write (nthreed, 780) zeta, rb(lk), zb(lk), 
     1            bredge(lk), bpedge(lk), bzedge(lk),
     2            brv(lk), bphiv(lk), bzv(lk)
               end do
            endif
         end do
      end do

      deallocate (rb, phib, zb, bredge, bpedge, bzedge, stat=l)

      write (nthreed, 800)
      mnlim = (mnmax/2)*2
      do mn = 1, mnlim, 2
         write (nthreed, 810) nint(xn(mn)/nfp), nint(xm(mn)), rmnc(mn),
     1      zmns(mn), bmodmn(mn), bmodmn1(mn), nint(xn(mn+1)/nfp),
     2      nint(xm(mn+1)), rmnc(mn+1), zmns(mn+1), bmodmn(mn+1),
     3      bmodmn1(mn+1)
      end do

      write (nthreed, 820)

      mnlim = mnpd/3
      do mn = 1, mnlim
         mn0 = 3*mn - 2
         write (nthreed, 830) nint(xnpot(mn0)/nfp), nint(xmpot(mn0)),
     1      potvac(mn0), nint(xnpot(mn0+1)/nfp), nint(xmpot(mn0+1)),
     2      potvac(mn0+1), nint(xnpot(mn0+2)/nfp), 
     3      nint(xmpot(mn0+2)), potvac(mn0+2)
      end do
      mnres = mnpd - mnlim*3
      mn = mnlim*3 + 1
      if (mnres .eq. 2) then
         write (nthreed, 830) nint(xnpot(mn)/nfp), nint(xmpot(mn)),
     1      potvac(mn), nint(xnpot(mn+1)/nfp), nint(xmpot(mn+1)),
     2      potvac(mn+1)
      else if (mnres .eq. 1) then
         write (nthreed, 830) nint(xnpot(mn)/nfp), nint(xmpot(mn)),
     1      potvac(mn)
      endif

  750 format(/,3x,'NF*PHI',7x,' Rb ',8x,' Zb ',6x,'BSQMHDI',5x,'BSQVACI'
     1   ,5x,'BSQMHDF',5x,'BSQVACF',/)
  760 format(/,3x,'NF*PHI',7x,' Rb ',8x,' Zb ',6x,'BR',8x,'BPHI',6x,'BZ'
     1   ,8x,'BRv',7x,'BPHIv',5x,'BZv',/)
  770 format(1pe10.2,1p6e12.4)
  780 format(1pe10.2,1p2e12.4,1p6e10.2)
  790 format(i5,/,(1p3e12.4))
  800 format(//,3x,'nb',2x,'mb',9x,'rbc',9x,'zbs',3x,'|B|(s=.5)',3x,
     1   '|B|(s=1.)',6x,'nb',2x,'mb',9x,'rbc',9x,'zbs',3x,'|B|(s=.5)',3x
     2   ,'|B|(s=1.)'/)
  810 format(i5,i4,1p4e12.4,3x,i5,i4,1p4e12.4)
  820 format(/,3x,'nf',2x,'mf',5x,'potvacs',6x,'nf',2x,'mf',5x,'potvacs'
     1   ,6x,'nf',2x,'mf',5x,'potvacs'/)
  830 format(i5,i4,1pe12.4,3x,i5,i4,1pe12.4,3x,i5,i4,1pe12.4)

      end subroutine freeb_data
      

      subroutine lamfull(luef, luof, lue, luo)
      use vmec_main
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,*), intent(out) :: luef, luof
      real(rprec), dimension(ns,nzeta,*), intent(in) :: lue, luo
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: p5 = 0.5_dp, c1p5 = 1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: lt, lt1
C-----------------------------------------------

C
C     COMPUTES LAMBDA ON FULL RADIAL MESH AT THETA=0, PI
C 
      do lt = 1, 2
         lt1 = 1
         if (lt .eq. 2) lt1 = ntheta2
         luof(1,lt1) = p5*(c1p5*luo(2,1,1)-p5*luo(3,1,1)+c1p5*luo(2,1,
     1      ntheta2)-p5*luo(3,1,ntheta2))
         luef(1,lt1) = p5*(c1p5*lue(2,1,1)-p5*lue(3,1,1)+c1p5*lue(2,1,
     1      ntheta2)-p5*lue(3,1,ntheta2))
         luof(ns,lt1) = c1p5*luo(ns,1,lt1) - p5*luo(ns1,1,lt1)
         luef(ns,lt1) = c1p5*lue(ns,1,lt1) - p5*lue(ns1,1,lt1)
         luof(2:ns1,lt1) = p5*(luo(3:ns1+1,1,lt1)+luo(2:ns1,1,lt1))
         luef(2:ns1,lt1) = p5*(lue(3:ns1+1,1,lt1)+lue(2:ns1,1,lt1))
      end do
 
      end subroutine lamfull

      
      subroutine magnetics_data
      use vmec_main
      use vacmod
      use vsvd
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: iobset, j, i, n, ibfld
      real(rprec) :: denwgt, dsitotal,     
     1   chierr, wght0, abs_flux, chiwgt, raddeg, denbwgt, btotal, 
     2   relerr_b, wghtb, chisqx, chibwgt
C-----------------------------------------------
!
!       ALGEBRAIC SUM OF FLUXES FOR GROUPS OF OBSERVER POSITIONS
!
      if (nobd .ne. 0) then
         flmwgt = 0
         denwgt = 0
         write (nthreed, 400)
  400    format(/' External Poloidal Flux Loop Measurements (Wb)',/,
     1      '   I     Connection          Coils       Plasma    ',
     2      '    Total     Measured   Chi-Sq Err        Sigma',
     3      '   Designation',/,1x,3('-'),1x,16('-'),7(1x,12('-')),'-')
 
         do iobset = 1, nobd
            dsitotal = dsiext(iobset) + plflux(iobset)
            chierr = dsiobt(iobset) - dsitotal
            wght0 = cbig + one
            if (any(iobset.eq.indxflx(:nflxs)))
     1          wght0 = sigma_flux(iobset)
            flmwgt = flmwgt + (chierr/wght0)**2
            denwgt = denwgt + (dsitotal/wght0)**2
            chierr = (chierr/wght0)**2
            if (wght0 .lt. cbig) then
               write (nthreed, 410) iobset, (iconnect(j,iobset),j=1,4), 
     1            dsiext(iobset), plflux(iobset), dsitotal, dsiobt(
     2            iobset), chierr, wght0, dsilabel(iobset)
            else
               write (nthreed, 420) iobset, (iconnect(j,iobset),j=1,4), 
     1            dsiext(iobset), plflux(iobset), dsitotal, dsiobt(
     2            iobset), chierr, dsilabel(iobset)
            endif
         end do
 
         write (nthreed, 600)
         do i = 1, nobser
            abs_flux = plflux(nobd+i)
            write (nthreed, 610) i, xobser(i), zobser(i), psiext(i), 
     1         abs_flux, psiext(i) + abs_flux
         end do
  600    format(//,' Poloidal Flux at Individual Coil Locations',/,
     1      '   I         R [M]        Z [M]        Coils       Plasma',
     2      '        Total',/,2x,3('-'),5(1x,12('-')))
  610    format(i4,1x,5f13.4)
      endif
 
      nchisaddle = nflxs
      total_saddle_chi = 0.
      if (nflxs .GT. 0) then
         chiwgt = flmwgt/(nflxs)
         flmwgt = sqrt(flmwgt/denwgt)
         total_saddle_chi = (nflxs)*chiwgt
         write (nthreed, 430) chiwgt, flmwgt
      endif
  410 format(i4,1x,4i4,6f13.4,7x,a8)
  420 format(i4,1x,4i4,5f13.4,5x,'No Match',7x,a8)
  430 format(/' MEAN CHI-SQ ERROR IN SADDLE LOOP MATCH: ',1pe10.3,/,
     1   ' RMS ERROR IN SADDLE LOOP MATCH: ',1pe10.3)
 
!
!       BR, BZ FIELD COMPONENTS FOR MATCHING AT OBSERVER POSITIONS
!
      total_b_chi = 0.
      if (nbfldn .ne. 0) then
  500    format(//' BFIELD GROUP: ',a,' // ',
     1      'External B Loop Magnetic Field Measurements (T)',/,
     2      '  I  Rcoil[M]  Zcoil[M]    Degree    B[coil]',
     3      '  B[plasma]   B[Total]  ','  B[Meas]    Sigma  Chi-Sq Err',
     4      /,1x,2('-'),3(1x,9('-')),4(1x,10('-')),1(1x,8('-')),1(1x,11(
     5      '-')))
 
         raddeg = 360.0_dp/twopi
         do n = 1, nbsets                        ! over all allowed sets
            write (nthreed, 500) bloopnames(n)
            bcwgt = 0.
            denbwgt = 0.
            b_chi(n) = 0.
            ibfld = 0
            do j = 1, nbcoils(n)
               ibfld = ibfld + 1
               btotal = bcoil(j,n) + plbfld(j,n)
               relerr_b = btotal - bbc(j,n)
               wghtb = cbig + 1.0_dp
c                                                ! over namelist sets
               if (any(j.eq.indxbfld(:nbfld(n),n))) wghtb = sigma_b(j,n)
               bcwgt = bcwgt + (relerr_b/wghtb)**2
               denbwgt = denbwgt + (btotal/wghtb)**2
               chisqx = (relerr_b/wghtb)**2
               if (wghtb .LT. cbig) then
                  write (nthreed, 520) ibfld, rbcoil(j,n), zbcoil(j,n), 
     1               abcoil(j,n)*raddeg, bcoil(j,n), plbfld(j,n), btotal
     2               , bbc(j,n), wghtb, chisqx
               else
                  write (nthreed, 530) ibfld, rbcoil(j,n), zbcoil(j,n), 
     1               abcoil(j,n)*raddeg, bcoil(j,n), plbfld(j,n), btotal
     2               , bbc(j,n)
               endif
            end do
            if (nbfld(n) .gt. 0) then
               chibwgt = bcwgt/(nbfld(n))
               bcwgt = sqrt(bcwgt/denbwgt)
               b_chi(n) = (nbfld(n))*chibwgt
               write (nthreed, 540) bloopnames(n), chibwgt, 
     1            bloopnames(n), bcwgt
               total_b_chi = total_b_chi + b_chi(n)
            endif
         end do
 
      endif
  520 format(i3,3f10.3,4f11.3,f9.4,f12.4)
  530 format(i3,3f10.3,4f11.3,3x,'No Match')
  540 format(/' MEAN CHI-SQ ERROR IN B-COIL GROUP ',a,': ',1pe10.3,/,
     1   ' RMS ERROR IN B-COIL GROUP ',a,': ',1pe10.3)
 
      end subroutine magnetics_data


      subroutine movephi1 (gphifacx)   !wiecode
      use vmec_main
      use vsvd
      implicit none 
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec) :: gphifacx
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec) :: phifac0
C-----------------------------------------------
!
!     UPDATE PHIEDGE SCALE FACTOR BY THE FORMULA
!     FDOT/F = -2*(1 - Ipexp/Ipcode), WHERE F = PHISCALE
!
      if (ctor .eq. zero) stop 'ctor = 0 in movephi1'
      phifac0 = phifac*(currv/ctor)
      gphifacx = rsfac*c1pm2*(phifac0 - phifac)
 
      end subroutine movephi1


      subroutine newphi(phipog)
      use kind_spec
      use vmec_main
      use realspace
      use vsvd
      use xstuff
c-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(nrzt) :: phipog
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec) :: phifold, tnorm
C-----------------------------------------------

!
!     UPDATE PHI SCALE FACTOR (PHIFAC)
!
      phifold = phifac
      phifac = xc(neqs1)

      if (phifold .eq. zero) stop 'phifac = 0 in newphi'
      tnorm = phifac/phifold
      phip = tnorm * phip
      phipog = tnorm * phipog
 
      end subroutine newphi


      subroutine printout(i0, delt0, w0, lscreen)
      use vmec_main
      use realspace
      use vsvd
!     use vspline
      use xstuff
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer i0
      real(rprec) :: delt0, w0
      logical lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec) :: betav, w, avm, den
C-----------------------------------------------
      betav = wp/wb
      w = w0*twopi*twopi
      den = zero
      specw(1) = one
      call spectrum (xc(:irzloff), xc(1+irzloff:2*irzloff))
      den = sum(phip(2:ns))
      avm = dot_product(phip(2:ns),specw(2:ns)+specw(:ns-1))
      avm = 0.5_dp*avm/den
      if (ivac .ge. 1 .and. iter2.gt.1) delbsq =
     1     sum(dbsq(:nznt)*wint(2:nrzt:ns))/
     1     sum(bsqsav(:nznt,3)*wint(2:nrzt:ns))
      if (i0.eq.1 .and. lfreeb) then
         if (lscreen) print 20
         if (imatch_phiedge .eq. 1) then
            write (nthreed, 15)
         else
            write (nthreed, 16)
         endif
      else if (i0.eq.1 .and. .not.lfreeb) then
         if (lscreen) print 30
         write (nthreed, 25)
      endif
   15 format(/,' ITER    FSQR      FSQZ       FSQL      fsqr      ',
     1   'fsqz      DELT    RAX(v=0)       WMHD      BETA',
     2   '      <M>   DEL-BSQ   FEDGE',/)
   16 format(/,' ITER    FSQR      FSQZ      FSQL      fsqr      ',
     1   'fsqz      DELT      RAX(v=0)     WMHD       BETA',
     2   '     PHIEDGE  DEL-BSQ  FEDGE',/)
   20 format(/,' ITER    FSQR      FSQZ      FSQL    ',
     1   'RAX(v=0)      WMHD      DEL-BSQ',/)
   25 format(/,' ITER    FSQR      FSQZ      FSQL      fsqr      ',
     1   'fsqz      DELT     RAX(v=0)      WMHD',
     2   '       BETA     <M>        ',/)
   30 format(/,' ITER    FSQR      FSQZ      FSQL      fsqr      ',
     1   'fsqz     RAX(v=0)     WMHD',/)
      if (.not.lfreeb) then
         if (lscreen) print 45, i0, fsqr, fsqz, fsql, fsqr1, fsqz1,r00,w
         write (nthreed, 40) i0, fsqr, fsqz, fsql, fsqr1, fsqz1, delt0, 
     1      r00, w, betav, avm
         return
      endif
      if (lscreen) print 50, i0, fsqr, fsqz, fsql, r00, w, delbsq
      if (imatch_phiedge .eq. 1) then
         write (nthreed, 40) i0, fsqr, fsqz, fsql, fsqr1, fsqz1, delt0, 
     1      r00, w, betav, avm, delbsq, fedge
      else
         write (nthreed, 40) i0, fsqr, fsqz, fsql, fsqr1, fsqz1, delt0, 
     1      r00, w, betav, abs(phiedge*phifac), delbsq, fedge
      endif
   40 format(i5,1p6e10.2,1pe11.3,1pe12.4,1pe11.3,1pe9.2,0pf7.3,1p2e9.2)
   45 format(i5,1p5e10.2,1pe11.3,1pe12.4)
   50 format(i5,1p3e10.2,1pe10.3,1pe11.3,1pe12.4,1pe11.3)
!     if (lrecon) then
!        if (ystark(1)*ystark(isnodes).lt.zero .or. 
!    1       ystark(2)*ystark(isnodes-1).lt.zero) then
!           if (lscreen) print 200
!           write (nthreed, 200)
!        endif
! 200    format(' CHECK THE SIGN OF MSE DATA: INCONSISTENT WITH ',
!    1      'CURRENT, TOROIDAL FIELD!')
!       if (pressum0 .ne. zero .and. apres .ne. zero .and. lscreen)
!    1     print 60, errsvd, pfac, phifac, 
!    2               abs(fsqsum0/pressum0), aminor/apres
!  60    format(2x,'chisq-err   = ',f8.3,3x,'pres-scale = ',f8.3,3x,
!    1      'phifac    = ',f8.3,/,2x,'<F00> force = ',f8.3,3x,
!    2      '<a>/apres  = ',f8.3)
!        if (imatch_phiedge.eq.2 .and. ivac>1) then
!           if (lscreen) print 65, dlim_min, rlim_min, zlim_min
!  65       format(2x,'dlim (min)  = ',f8.3,3x,'at rlim    = ',f8.3,
!    1         '   zlim = ',f8.3/)
!        else
!           if (lscreen) print *, ' '
!        endif
!     endif
      
      end subroutine printout


      subroutine spectrum(rmn, zmn)
      use vmec_main
      use vmec_params, only: mscale, nscale, ntmax
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,0:ntor,0:mpol1,ntmax) :: 
     1   rmn, zmn
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: js, ntype, n, m
      real(rprec), dimension(ns) :: t1, dnumer, denom
      real(rprec) :: scale
C-----------------------------------------------
 
      dnumer(2:ns) = zero
      denom(2:ns) = zero
      do ntype = 1,ntmax
        do n = 0,ntor
          do m = 1,mpol1
             scale = (mscale(m)*nscale(n))**2
             do js = 2,ns
                t1(js) =(rmn(js,n,m,ntype)**2 + zmn(js,n,m,ntype)**2)
     2               *scale
             end do
             dnumer(2:ns) = dnumer(2:ns) + t1(2:ns)*xmpq(m,3)
             denom (2:ns) = denom (2:ns) + t1(2:ns)*xmpq(m,2)
          end do
        end do
      enddo

      specw(2:ns) = dnumer(2:ns)/denom(2:ns)

      end subroutine spectrum


      subroutine readin(input_file, iseq_count, lfirst, ier_flag,
     1      lscreen)
      use vmec_main
      use vmec_params, only: ntmax, signgs
      use vacmod
      use vsvd
      use vspline
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer iseq_count, ier_flag
      logical lfirst, lscreen
      character*(*) :: input_file
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: ns_default = 31
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: iexit, ipoint, n, iunit, ier_flag_init,
     1   i, n1, m, nsmin, igrid, m1, isgn
      real(rprec), dimension(:,:), pointer ::
     1  rbcc, rbss, rbcs, rbsc, zbcs, zbsc, zbcc, zbss
      real(rprec) rtest, ztest, treadon, treadoff, tzc, trc, delta
      character*(80) :: line
C-----------------------------------------------
!
!       LOCAL VARIABLES
!
!       rbcc,rbss,rbcs,rbsc
!                boundary Fourier coefficient arrays for R (of cosu*cosv, etc)
!       zbcc,zbss,zbcs,zbsc
!                boundary Fourier coefficient arrays for Z

!
!       STACKING ORDER (DEPENDS ON LASYM):
!       XCC*COS(MU)COS(NV), XCS*COS(MU)SIN(NV), ETC
!               LASYM = T                             LASYM = F
!          1      ,  mns => rmncc                1      ,  mns => rmncc
!          1+  mns,2*mns => rmnss                1+  mns,2*mns => rmnss
!          1+2*mns,3*mns => rmncs
!          1+3*mns,4*mns => rmnsc
!          1+4*mns,5*mns => zmncs                1+2*mns,3*mns => zmncs
!          1+5*mns,6*mns => zmnsc                1+3*mns,4*mns => zmnsc
!          1+6*mns,7*mns => zmncc
!          1+7*mns,8*mns => zmnss
!          1+8*mns,9*mns => lmncs                1+4*mns,5*mns => lmncs
!          1+9*mns,10*mns=> lmnsc                1+5*mns,neqs  => lmnsc
!          1+10*mns,11*mns=> lmncc
!          1+11*mns,neqs  => lmnss
!          neqs1          => phifac
!          neqs2          => MSE axis correction
!


!
!                STANDARD INPUT DATA AND RECOMMENDED VALUES
!
!   Plasma parameters (MKS units)
!          ai:   iota (ncurr=0) or toroidal current density (=ac, ncurr=1)
!                expansion coefficients (series in s)
!          am:   mass or pressure (gamma=0) expansion coefficients (series in s)
!                in MKS units (NWT/M**2)
!      curtor:   value of toroidal current (A). Used if ncurr = 1 to specify
!                current profile, or if in data reconstruction mode.
!      extcur:   array of currents in each external current group. Used to
!                multiply Green''s function for fields and loops read in from
!                MGRID file. Should use real current units (A).
!       gamma:   value of compressibility index (gamma=0 => pressure prescribed)
!         nfp:   number of toroidal field periods ( =1 for Tokamak)
!         rbc:   boundary coefficients of cos(m*theta-n*zeta) for R
!         zbs:   boundary coefficients of sin(m*theta-n*zeta) for Z
!         rbs:   boundary coefficients of sin(m*theta-n*zeta) for R
!         zbc:   boundary coefficients of cos(m*theta-n*zeta) for Z
!
!
!   Numerical parameters
!  mgrid_file:   full path for vacuum green''s function data
!       ncurr:   flux conserving (=0) or prescribed toroidal current (=1)
!    ns_array:   array of radial mesh sizes to be used in multigrid sequence
!
!   Convergence control parameters
!  ftol_array:   array of value of residual(s) at which each multigrid
!                iteration ends
!       niter:   number of iterations (used to terminate run)
!       nstep:   number of timesteps between printouts on screen
!    nvacskip:   iterations skipped between full update of vacuum solution
!       tcon0:   weight factor for constraint force (=1 by default)
!
!   Equilibrium reconstruction parameters
!      phifac:   factor scaling toroidal flux to match apres or limiter
!   datastark:   pitch angle data from stark measurement
!    datathom:   pressure data from Thompson, CHEERS (Pa)
!     imatch_         = 1 (default),match value of PHIEDGE in input file
!     phiedge:   = 0, use pressure profile width to determine PHIEDGE
!                = 2, use LIMPOS data (in mgrid file) to find PHIEDGE
!                = 3, use Ip to find PHIEDGE (fixed-boundary only)
!        imse:   number of Motional Stark effect data points
!                >0, use mse data to find iota; <=0, fixed iota profile ai
!        itse:   number of pressure profile data points
!                = 0, no thompson scattering data to read
!     isnodes:   number of iota spline points (computed internally unless specified explicitly)
!     ipnodes:   number of pressure spline points (computed internally unless specified explicitly)
!       lpofr:   logical variable. =.true. if pressure data are
!                prescribed in real space. =.false. if data in flux space.
!      pknots:   array of pressure knot values in sqrt(s) space
!      sknots:   array of iota knot values in sqrt(s) space
!       tensp:   spline tension for pressure profile
!
!       tensi:   spline tension for iota
!      tensi2:   vbl spline tension for iota
!      fpolyi:   vbl spline tension form factor (note: if tensi!=tensi2
!               then tension(i-th point) = tensi+(tensi2-tensi)*(i/n-1))**fpolyi
!               - - - - - - - - - - - - - - - - - -
!    mseangle_   uniform experimental offset of MSE data
!     offset:    (calibration offset) ... PLUS ...
!    mseangle_   multiplier on mseprof offset array
!     offsetM:   (calibration offset)
!     mseprof:   offset array from namelist MSEPROFIL
!                so that the total offset on the i-th MSE data point is
!                taken to be
!                = mseangle_offset+mseangle_offsetM*mseprof(i)
!               - - - - - - - - - - - - - - - - - -
! pres_offset:   uniform arbitrary  radial offset of pressure data
!     presfac:   number by which Thomson scattering data is scaled
!                to get actual pressure
!     phidiam:   diamagnetic toroidal flux (Wb)
!      dsiobt:   measured flux loop signals corresponding to the
!                combination of signals in iconnect array
!     indxflx:   array giving index of flux measurement in iconnect array
!    indxbfld:   array giving index of bfield measurement used in matching
!        nobd:   number of connected flux loop measurements
!      nobser:   number of individual flux loop positions
!      nbsets:   number of B-coil sets defined in mgrid file
!  nbcoils(n):   number of bfield coils in each set defined in mgrid file
!    nbcoilsn:   total number of bfield coils defined in mgrid file
!    bbc(m,n):   measured magnetic field at rbcoil(m,n),zbcoil(m,n) at
!                the orientation br*cos(abcoil) + bz*sin(abcoil)
! rbcoil(m,n):   R position of the m-th coil in the n-th set from mgrid file
! zbcoil(m,n):   Z position of the m-th coil in the n-th set from mgrid file
! abcoil(m,n):   orientation (surface normal wrt R axis; in radians)
!                of the m-th coil in the n-th set from mgrid file.
!       nflxs:   number of flux loop measurements used in matching
!    nbfld(n):   number of selected external bfield measurements in set n from nml file
!      nbfldn:   total number of external bfield measurements used in matching
!               - - - - - - - - - - - - - - - - - -
!             NOTE: FOR STANDARD DEVIATIONS (sigma''s) < 0, INTERPRET
!             AS PERCENT OF RESPECTIVE MEASUREMENT
!  sigma_thom:   standard deviation (Pa) for pressure profile data
! sigma_stark:   standard deviation (degrees) in MSE data
!  sigma_flux:   standard deviaton (Wb) for external poloidal flux data
!     sigma_b:   standard deviation (T) for external magnetic field data
!sigma_current:  standard deviation (A) in toroidal current
!sigma_delphid:  standard deviation (Wb) for diamagnetic match
!
!                - - - - - - - - - - - - - - - - - -
!               INPUTS IN MAKEGRID CODE WHICH ARE EXPORTED
!               IN MGRID FILE
!                - - - - - - - - - - - - - - - - - -
!               FOR SPECIFYING BFIELD GREEN''S FUNCTIONS
!        rmin:   minimum R position for green''s function data
!        rmax:   maximum R position for green''s function data
!        zmin:   minimum Z position for green''s function data
!        zmax:   maximum Z position for green''s function data
!          ir:   parameter specifying number of points in R mesh
!                for green''s function box
!          jz:   parameter specifying number of points in Z mesh
!                for green''s function box
!          kp:   parameter specifying number of toroidal planes
!                within one field period used for green''s
!                function box. MUST agree with NV parameter in
!                VMEC for now
!     nextcur:         no. of external current groups (eg., TF, PF, helical)
!    curlabel:   array of labels describing each current group
!                     included in green''s function BFIELD response
!
!               - - - - - - - - - - - - - - - - - -
!               FOR DIAGNOSTICS AND DATA ANALYSIS
!               (HERE,COILS ARE FOR MEASURING FIELDS, FLUXES)
!    iconnect:   two-dimensional array describing electrical
!                      connection of up to four flux loops. Specifies
!                     the sign and flux loop number of (up to) four
!                     connected individual loops (indexing based on
!                     xobser,zobser arrays).
!      nobser:   no. of distinct flux loops
!        nobd:   no. of distinct connections of flux loops
!      nbsets:   no. of b-field coil sets
!     nbcoils:   array specifying no. of b-field coils in each set
!  bloopnames:   array of labels describing b-field sets
!    dsilabel:   array of labels describing connected flux loops
!      xobser:   array of flux loop R-positions
!      zobser:   array of flux loop Z-positions
!      rbcoil:   two-dimensional array of b-field coil R-positions
!                (no. of coil groups, no. coils in a group)
!      zbcoil:   two-dimensional array of b-field coil Z-positions
!      abcoil:   two-dimensional array of angles (in degrees) of
!                       each b-field coil
!
!
!
!       THE (ABSOLUTE) CHI-SQ ERROR IS DEFINED AS FOLLOWS:
!
!          2
!       CHI      =     SUM [ EQ(K,IOTA,PRESSURE)  -  DATA(K) ] ** 2
!                     (K) -----------------------------------
!                                   SIGMA(K)**2
!
!       HERE, SIGMA IS THE STANDARD DEVIATION OF THE MEASURED DATA, AND
!       EQ(IOTA,PRESSURE) IS THE EQUILIBRIUM EXPRESSION FOR THE DATA TO BE
!       MATCHED:
!
!       EQ(I)   =    SUM [ W(I,J)*X(J) ]
!                   (J)
!
!       WHERE W(I,J) ARE THE (LINEAR) MATRIX ELEMENTS AND X(J) REPRESENT
!       THE KNOT VALUES OF IOTA (AND/OR PRESSURE). THE RESULTING LEAST-SQUARES
!       MATRIX ELEMENTS AND DATA ARRAY CAN BE EXPRESSED AS FOLLOWS:
!
!       ALSQ(I,J) = SUM [ W(K,I) * W(K,J) / SIGMA(K) ** 2]
!                   (K)
!
!       BLSQ(I)   = SUM [ W(K,I) * DATA(K)/ SIGMA(K) ** 2]
!                   (K)
!
!       THEREFORE, INTERNALLY IT IS CONVENIENT TO WORK WITH THE 'SCALED'
!       W'(K,I) = W(K,I)/SIGMA(K) AND DATA'(K) = DATA(K)/SIGMA(K)
!
!       ****!   I - M - P - O - R - T - A - N - T     N - O - T - E   *****
!
!       THE INPUT DATA FILE WILL ACCEPT BOTH POSITIVE AND NEGATIVE
!       SIGMAS, WHICH IT INTERPRETS DIFFERENTLY. FOR SIGMA > 0, IT
!       TAKES SIGMA TO BE THE STANDARD DEVIATION FOR THAT MEASUREMENT
!       AS DESCRIBED ABOVE. FOR SIGMA < 0, SIGMA IS INTERPRETED AS
!       THE FRACTION OF THE MEASURED DATA NEEDED TO COMPUTE THE ABSOLUTE
!       SIGMA, I.E., (-SIGMA * DATA) = ACTUAL SIGMA USED IN CODE.
!
      ier_flag_init = ier_flag
      ier_flag = 0
      call second0(treadon)
      if (ier_flag_init .eq. 4) goto 1000

!
!     READ IN DATA FROM INDATA FILE
!
      call read_indata(input_file, iunit, ier_flag)
      if (ier_flag .ne. 0) return


      if (tensi2 .eq. zero ) tensi2 = tensi

!
!     open output files here, print out heading to threed1 file
!
      call heading(input_extension, time_slice,
     1      iseq_count, lmac, lscreen, lfirst)

!
!     READ IN COMMENTS DEMARKED BY "!"
!
      rewind (iunit)
      iexit = 0
      do while( iexit.eq.0 )
         read (iunit, '(a)') line
         iexit = index(line,'INDATA') + index(line,'indata')
         ipoint = index(line,'!')
         if (ipoint .eq. 1) write (nthreed,*) line
      enddo
      close (iunit)

!
!     READ IN AND STORE (FOR SEQUENTIAL RUNNING) MAGNETIC FIELD DATA
!     FROM MGRID_FILE FIRST TIME (lfirst = T) ONLY
!     SET LOGICAL FLAGS FOR ALL SUBSEQUENT RUNS
!
      if (lfirst .and. lfreeb) call read_mgrid (lscreen, ier_flag)
      if (ier_flag .ne. 0) return

!
!     PARSE NS_ARRAY
!
      nsin = max (3, nsin)
      multi_ns_grid = 1
      if (ns_array(1) .eq. 0) then                   !!Old input style
          ns_array(1) = min(nsin,nsd)
          multi_ns_grid = 2
          ns_array(multi_ns_grid) = ns_default        !!Run on 31-point mesh
      else
          nsmin = 1
          do while (ns_array(multi_ns_grid) .gt. nsmin .and.
     1             multi_ns_grid .lt. 100)
             nsmin = max(nsmin, ns_array(multi_ns_grid))
            if (nsmin.le.nsd) then
               multi_ns_grid = multi_ns_grid + 1
            else                                      !!Optimizer, Boozer code overflows otherwise
               ns_array(multi_ns_grid) = nsd
               nsmin = nsd
               print *,' NS_ARRAY ELEMENTS CANNOT EXCEED ',nsd
               print *,' CHANGING NS_ARRAY(',multi_ns_grid,') to ', nsd
            end if
          end do
          multi_ns_grid = multi_ns_grid - 1
      endif
      if (ftol_array(1) .eq. zero) then
         ftol_array(1) = 1.e-8_dp
         if (multi_ns_grid .eq. 1) ftol_array(1) = ftol
         do igrid = 2, multi_ns_grid
            ftol_array(igrid) = 1.e-8_dp * (1.e8_dp * ftol)**
     1        ( real(igrid-1,rprec)/(multi_ns_grid-1) )
         end do
      endif

!
!     WRITE OUT DATA TO THREED1 FILE
!
      write (nthreed,100)
     1  ns_array(multi_ns_grid),ntheta1,nzeta,mpol,ntor,nfp,
     2  gamma,spres_ped,phiedge,curtor
 100  format(/,' COMPUTATION PARAMETERS: (u = theta, v = zeta)'/,
     1  1x,45('-'),/,
     2  '     ns     nu     nv     mu     mv',/,
     3  5i7,//,' CONFIGURATION PARAMETERS:',/,1x,25('-'),/,
     4  '    nfp      gamma      spres_ped    phiedge(wb)     curtor(A)'
     5  ,/,i7,1pe11.3,1p2e15.3,1pe14.3,/)

      if (nvacskip.le.0) nvacskip = nfp
      write (nthreed,110) ncurr,niter,ns_array(1),nstep,nvacskip,
     1  ftol_array(multi_ns_grid),tcon0,lasym
 110    format(' RUN CONTROL PARAMETERS:',/,1x,23('-'),/,
     1  '  ncurr  niter   nsin  nstep  nvacskip      ftol     tcon0',
     2  '   lasym',/, 4i7,i10,1p2e10.2,L8/)
      if (nextcur.ne.0) then
        write(nthreed,120)
        do i = 1,nextcur,5
          write (nthreed, 125) trim(curlabel(i)), trim(curlabel(i+1)),
     1    trim(curlabel(i+2)), trim(curlabel(i+3)), trim(curlabel(i+4)),
     2    extcur(i), extcur(i+1), extcur(i+2), extcur(i+3), extcur(i+4)
        enddo
      endif
 120  format(' EXTERNAL CURRENTS',/,1x,17('-'))
 125  format(5a12,/,1p5e12.4,/)

      if (bloat .ne. one) then
          write (nthreed,'(" Profile Bloat Factor: ",f7.5)') bloat
          phiedge = phiedge*bloat
      endif

      write(nthreed,130)
 130  format(' MASS PROFILE COEFFICIENTS am - newton/m**2',
     1  ' (EXPANSION IN TOROIDAL FLUX):',/,1x,35('-'))
      write(nthreed,135)(am(i-1),i=1, size(am))
      if (ncurr.eq.0) then
          write(nthreed,140)
          write(nthreed,135)(ai(i-1),i=1, size(ai))
      else
          select case(ac_form)
          case (1)
             write(nthreed,1145)
          case default
             write(nthreed,145)
          end select
             
          write(nthreed,135)(ac(i-1),i=1, size(ac))
      endif
!     write(nthreed,150)
!     write(nthreed,135)(aphi(i-1),i=1, size(aphi))    !!Temporarily disabled

 135  format(1p6e12.3)
 140  format(/' IOTA PROFILE COEFFICIENTS ai',
     1   ' (EXPANSION IN TOROIDAL FLUX):',/,1x,35('-'))
 145  format(/' TOROIDAL CURRENT DENSITY (*V'') COEFFICIENTS',
     1        ' ac (EXPANSION IN TOROIDAL FLUX):',/,1x,38('-'))
 1145 format(/' TOROIDAL CURRENT DENSITY (*V'') COEFFICIENTS',
     1        ' ac (EXPANSION IN TOROIDAL FLUX w/ SQRT Term):',
     2          /,1x,38('-'))
 150  format(/' NORMALIZED TOROIDAL FLUX COEFFICIENTS aphi',
     1   ' (EXPANSION IN S):',/,1x,35('-'))
      write(nthreed,180)
 180  format(/,' R-Z FOURIER BOUNDARY COEFFICIENTS',/,
     1  ' R = RBC*cos(m*u - n*v) + RBS*sin(m*u - n*v),',
     2  ' Z = ZBC*cos(m*u - n*v) + ZBS*sin(m*u-n*v)'/1x,86('-'),
     3  /,'   nb  mb     rbc         rbs         zbc         zbs   ',
     4   '     raxis(cc)   raxis(cs)   zaxis(cc)   zaxis(cs)')

 1000  continue

      if (.not.lasym) then
!
!       CONVERT TO REPRESENTATION WITH RBS(m=1) = ZBC(m=1)
!

      delta = atan( (rbs(0,1) - zbc(0,1))/(rbc(0,1) + zbs(0,1)) )
      if (delta .ne. zero) then
        do m = 0,mpol1
          do n = -ntor,ntor
            trc = rbc(n,m)*cos(m*delta) + rbs(n,m)*sin(m*delta)
            rbs(n,m) = rbs(n,m)*cos(m*delta) - rbc(n,m)*sin(m*delta)
            rbc(n,m) = trc
            tzc = zbc(n,m)*cos(m*delta) + zbs(n,m)*sin(m*delta)
            zbs(n,m) = zbs(n,m)*cos(m*delta) - zbc(n,m)*sin(m*delta)
            zbc(n,m) = tzc
          enddo
        enddo
      endif

      endif

!
!     ALLOCATE MEMORY FOR NU, NV, MPOL, NTOR SIZED ARRAYS
!
      call allocate_nunv

!
!     CONVERT TO INTERNAL REPRESENTATION OF MODES
!
!     R = RBCC*COS(M*U)*COS(N*V) + RBSS*SIN(M*U)*SIN(N*V)
!         + RBCS*COS(M*U)*SIN(N*V) + RBSC*SIN(M*U)*COS(N*V)
!     Z = ZBCS*COS(M*U)*SIN(N*V) + ZBSC*SIN(M*U)*COS(N*V)
!         + ZBCC*COS(M*U)*COS(N*V) + ZBSS*SIN(M*U)*SIN(N*V)
!

!     POINTER ASSIGNMENTS (NOTE: INDICES START AT 1, NOT 0, FOR POINTERS, EVEN THOUGH
!                          THEY START AT ZERO FOR RMN_BDY)

      rbcc => rmn_bdy(:,:,1)
      rbss => rmn_bdy(:,:,2)
      zbcs => zmn_bdy(:,:,1)
      zbsc => zmn_bdy(:,:,2)

      if (lasym) then
      rbcs => rmn_bdy(:,:,ntmax-1)
      rbsc => rmn_bdy(:,:,ntmax)
      zbcc => zmn_bdy(:,:,ntmax-1)
      zbss => zmn_bdy(:,:,ntmax)
      endif

      rmn_bdy(0:ntor,0:mpol1,:ntmax) = zero
      zmn_bdy(0:ntor,0:mpol1,:ntmax) = zero

      do 190 m=0,mpol1
         m1 = m+1
         do 190 n=-ntor,ntor
            n1 = abs(n) + 1
            if (n .eq. 0) then
               isgn = 0
            else if (n .gt. 0) then
               isgn = 1
            else
               isgn = -1
            end if
            rbcc(n1,m1) = rbcc(n1,m1) + rbc(n,m)
            rbss(n1,m1) = rbss(n1,m1) + isgn*rbc(n,m)
            zbcs(n1,m1) = zbcs(n1,m1) - isgn*zbs(n,m)
            zbsc(n1,m1) = zbsc(n1,m1) + zbs(n,m)

            if (lasym) then
            rbcs(n1,m1) = rbcs(n1,m1) - isgn*rbs(n,m)
            rbsc(n1,m1) = rbsc(n1,m1) + rbs(n,m)
            zbcc(n1,m1) = zbcc(n1,m1) + zbc(n,m)
            zbss(n1,m1) = zbss(n1,m1) + isgn*zbc(n,m)
            if (m .eq. 0) zbss(n1,m1) = zero
            if (m .eq. 0) rbsc(n1,m1) = zero
            end if

            if (ier_flag_init .ne. 0) cycle
            trc = abs(rbc(n,m)) + abs(rbs(n,m))
     1          + abs(zbc(n,m)) + abs(zbs(n,m))
            if (m .eq. 0) then
               rbss(n1,m1) = zero
               zbsc(n1,m1) = zero
               if (n .lt. 0) cycle
               if (trc.eq.zero .and. abs(raxis(n,1)).eq.zero .and.
     1             abs(zaxis(n,1)).eq.zero) cycle
               write (nthreed,195) n, m, rbc(n,m), rbs(n,m),
     1                   zbc(n,m), zbs(n,m), raxis(n,1), raxis(n,2),
     2                   zaxis(n,2), zaxis(n,1)
            else
               if (trc .eq. zero) cycle
               write (nthreed,195) n, m, rbc(n,m), rbs(n,m),
     1                   zbc(n,m), zbs(n,m)
            end if
 190  continue
 195  format(i5,i4,1p8e12.4)


!
!     CHECK SIGN OF JACOBIAN (SHOULD BE SAME AS SIGNGS)
!
      m = 1
      m1 = m+1
      rtest = sum(rbcc(1:ntor1,m1))
      ztest = sum(zbsc(1:ntor1,m1))
      signgs = one
      if (rtest*ztest .gt. zero) signgs = -one


      iresidue = -1
      if (lrecon) then
!
!       DETERMINE CURRENT-FLUX CONSISTENCY CHECK
!
        signiota = one
        if (signgs*curtor*phiedge .lt. zero)signiota = -one
        if (sigma_current .eq. zero) then
          write (*,*) 'Sigma_current cannot be zero!'
          ier_flag = 1
          return
        end if

!
!       SET UP RECONSTRUCTION FIXED PROFILES
!
        dcon = atan(one)/45
        call readrecon                   !Setup for reconstruction mode
        call fixrecon(ier_flag)          !Fixed arrays for reconstruction
        if (ier_flag .ne. 0) return
      end if

      currv = dmu0*curtor              !Convert to Internal units

      call second0(treadoff)
      timer(2) = timer(2) + (treadoff-treadon)

      end subroutine readin


      subroutine read_indata(in_file, iunit, ier_flag)
      use vmec_main
      use vmec_input, only: bloat, ncurr
      use vmec_params, only: ntmax
      use vacmod
      use read_namelist_mod
      use safe_open_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ier_flag, iunit
      character*(*) :: in_file
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ireadseq, lt, iosnml
      real(rprec) :: denom
C-----------------------------------------------
!     INDATA0 MUST BE CLOSED EXTERNAL TO THIS RETURN

      iunit = indata0
      call safe_open (iunit, ireadseq, in_file, 'old', 'formatted')
      if (ireadseq .ne. 0) then
         print *, ' In VMEC, error opening input file: ',
     1   trim(in_file), '. Iostat = ', ireadseq
         ier_flag = 5
         return
      endif


      call read_namelist (iunit, iosnml, 'indata')
      if (iosnml .ne. 0) then
         print *, ' In VMEC, indata namelist error: iostat = ', iosnml
         ier_flag = 5
         return
      endif

!
!     This is temporarily disabled in torflux, torflux_deriv routines...
!
      aphi(0) = zero
      denom = sum (aphi)
      if (denom .ne. zero) then
         aphi = aphi/denom
      else
         aphi = ( / zero, one, (zero, lt=2,10) / )
         denom = sum (aphi)
         aphi = aphi/denom
      end if   
 
      call read_namelist (iunit, iosnml, 'mseprofile')

      if (lrecon .and. itse.le.0 .and. imse.le.0) lrecon = .false.
      if (lfreeb .and. mgrid_file.eq.'NONE') lfreeb = .false.

      if (bloat .eq. zero) bloat = one
      if ((bloat.ne.one) .and. (ncurr.ne.1)) then
         ier_flag = 3
         return 
      endif
!
!     COMPUTE NTHETA, NZETA VALUES
!
      mpol = abs(mpol)
      ntor = abs(ntor)
      if (mpol .gt. mpold) stop 'mpol>mpold: lower mpol'
      if (ntor .gt. ntord) stop 'ntor>ntord: lower ntor'
      mpol1 = mpol - 1
      ntor1 = ntor + 1
      if (ntheta .le. 0) ntheta = 2*mpol + 6    !number of theta grid points (>=2*mpol+6)
      ntheta1 = 2*(ntheta/2)
      ntheta2 = 1 + ntheta1/2                   !u = pi
      if (ntor .eq. 0) lthreed = .false.
      if (ntor .gt. 0) lthreed = .true.

      if (nzeta .le. 0) nzeta = 2*ntor + 4      !number of zeta grid points (=1 if ntor=0)
      if (ntor .eq. 0) nzeta = 1
      mnmax = ntor1 + mpol1*(1 + 2*ntor)        !size of rmnc,  rmns,  ...
      mnsize = mpol*ntor1                       !size of rmncc, rmnss, ...

      mf = mpol + 1
      nf = ntor
      nu = ntheta1
      nv = nzeta
      mf1 = 1 + mf
      nf1 = 2*nf + 1
      mnpd = mf1*nf1
      nfper = nfp
!     NEED FOR INTEGRATION IN BELICU FOR AXISYMMETRIC PLASMA
      if (nf.eq.0 .and. nv.eq.1) nfper = 64

      if (lasym) then
         ntmax = 4
         ntheta3 = ntheta1
         mnpd2 = 2*mnpd
      else
         ntmax = 2
         ntheta3 = ntheta2
         mnpd2 = mnpd
      end if
      
      nvp = nv*nfper
      nuv = nu*nv
      nu2 = nu/2 + 1
      nu3 = ntheta3
      nznt = nzeta*ntheta3
      nuv2 = nznt
!     if (nuv2 < mnpd) then
!        print *, ' nuv2 < mnpd: not enough integration points'
!        stop 11
!     endif

      if (ncurr.eq.1 .and. all(ac.eq.cbig)) ac = ai            !!Old format: may not be reading in ac
      where (ac .eq. cbig) ac = zero
      
      end subroutine read_indata

 
      subroutine read_mgrid(lscreen, ier_flag)
      use vmec_main
      use vacmod
      use vsvd
      use vspline
      use mgrid_mod
      use safe_open_mod
      use system_mod
      implicit none
C-----------------------------------------------
C   D u m m y   V a r i a b l e s
C-----------------------------------------------
      integer :: ier_flag
      logical :: lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
#ifdef VMS
      character*(*), parameter :: mgrid_defarea = 'vmec$:[makegrid]'
#else
      character*(*), parameter :: mgrid_defarea = '$HOME/vmec/MAKEGRID'
#endif
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ipath, ierr, ierr1, ierr2, istat0, istat, ii, 
     1   i, j, n, n1, m, nsets_max, k, iunit = nbfld0
      character*200 :: home_dir
      logical :: lgrid_exist
C-----------------------------------------------
      mgrid_path = trim(mgrid_file)
      inquire (file=mgrid_path,exist=lgrid_exist)
      if (.not.lgrid_exist) then
          if (lscreen) print *,' MGRID FILE NOT FOUND IN SPECIFIED ',
     1       'PATH: TRYING TO FIND MGRID-FILE IN DEFAULT AREA'
          ii = index(mgrid_file,'/',back=.true.)
          istat = index(mgrid_defarea, '$HOME')
          if (istat .ne. 0) then
             call getenv('HOME', home_dir)
             if (istat .gt. 1) then
                home_dir = mgrid_defarea(1:istat-1) // trim(home_dir)
     1                 // trim(mgrid_defarea(istat+5:))
             else
                home_dir = trim(home_dir) // 
     1                     trim(mgrid_defarea(istat+5:))
             end if
          else
             home_dir = mgrid_defarea
          end if      
          mgrid_path = trim(home_dir) // mgrid_file(ii+1:)
          inquire (file=mgrid_path,exist=lgrid_exist)
      end if

      if (lgrid_exist) then
          call safe_open(iunit, ipath, mgrid_path, 'old', 'unformatted')
          if (ipath .ne. 0)then
             if (lscreen) then
                print *,' ',trim(mgrid_file(ii+1:)),' NOT FOUND IN ', 
     1                      trim(home_dir) 
                print *, ' MUST SUPPLY VACUUM BFIELD ON GRID TO ',
     1                   'RUN FREE-BOUNDARY'
                print *,' WILL PROCEED TO RUN VMEC IN',
     1                  ' FIXED BOUNDARY MODE'
             end if
             lfreeb = .false.
          endif
      else
         lfreeb = .false.
      end if

      if (.not.lfreeb) then
        lrecon = .false.
        return
      end if

! 
!     PARTITION B-LOOPS INTO SETS
!
      if (lscreen)
     1    print '(2x,2a)','Opening vacuum field file: ',mgrid_file

      read(iunit,iostat=ierr) nr0b, nz0b, np0b, nfper0, nextcur
      if (ierr .ne. 0) ier_flag = 9
      read(iunit,iostat=ierr) rminb, zminb, rmaxb, zmaxb
      if (ierr .ne. 0) ier_flag = 9

      if (np0b .ne. nv) then
        print *,' NV=',nv,' NOT EQUAL TO NP0B=',np0b,' IN MGRID FILE'
        ier_flag = 9
      else if ((nf.ne.0) .and. (nfper0.ne.nfp)) then
        print *,' NFP(read in) = ',nfp,' DOES NOT AGREE WITH ',
     1    'NFPER (in vacuum field file) = ',nfper0
        ier_flag = 9
      else if (nextcur .le. 0) then
        print *,' NEXTCUR = ',nextcur,' IN READING MGRID FILE'
        ier_flag = 9
      end if

      if (ier_flag .eq. 10) return
      
      allocate (curlabel(5*(nextcur/5+1)), stat=istat0)    !min of 5 for printing
      curlabel = " "
      read(iunit,iostat=istat) (curlabel(n),n=1,nextcur)
      if (istat.ne.0 .or. istat0.ne.0) then
         print *,' reading mgrid-files (curlabel) failed'
         ier_flag = 9
         return
      end if

!
!       NOTE: IF THE BTEMP ARRAY GETS TOO LARGE, USER CAN
!           ADD UP BVAC DIRECTLY FOR EACH EXTERNAL CURRENT GROUP
!           IN LOOP 50 WITHOUT STORING BTEMP(...,II)

      nbvac = nr0b*nz0b*nv
      allocate (btemp(nbvac,3,nextcur), stat=istat)
      if (istat .ne. 0) then
        print *,' allocation for b-vector storage failed'
        ier_flag = 9
        return
      end if

      do ii = 1,nextcur
         read(iunit) (btemp(i,1,ii),btemp(i,2,ii),btemp(i,3,ii),
     1              i=1,nbvac)
      enddo
      
!
!       READ IN EXTERNAL POLOIDAL FLUX, FIELD MEASURMENT
!       LOOP COORDINATES AND LABELS
!
      read(iunit,iostat=ierr1) nobser, nobd, nbsets
      if (ierr1.ne.0) then
         nobser = 0
         nobd   = 0
         nbsets = 0
         if (lscreen) print *,' No observation data in mgrid data'
         go to 900
      end if

      nbfldn = sum(nbfld(:nbsets))
      allocate (nbcoils(nbsets), stat=istat0)
      read(iunit) (nbcoils(n),n=1,nbsets)

      nbcoil_max = maxval(nbcoils(:nbsets))

      allocate (xobser(nobser), zobser(nobser), dsilabel(nobd),
     1       iconnect(4,nobser+nobd), unpsiext(nobser,nextcur),
     2       xobsqr(nobser), needflx(nobser), plflux(nobser+nobd),
     3       dsiext(nobd), psiext(nobser), bloopnames(nbsets),
     4       needbfld(nbcoil_max,nbsets), plbfld(nbcoil_max,nbsets),
     5       rbcoil(nbcoil_max,nbsets), zbcoil(nbcoil_max,nbsets),
     6       abcoil(nbcoil_max,nbsets), bcoil(nbcoil_max,nbsets),
     7       rbcoilsqr(nbcoil_max,nbsets), b_chi(nbsets),
     8       dbcoil(nbcoil_max,nbsets,nextcur), stat = istat)
      if (istat .ne. 0) then
          if (lscreen) 
     1       print *,' allocation error for xobser: istat = ',istat
          ier_flag = 9
          return
      end if
     
      if (nobser .gt. nfloops) then
         print *, 'NOBSER>NFLOOPS'
         ier_flag = 9
      end if   
      if (nobd .gt. nfloops) then
         print *, 'NOBD>NFLOOPS'
         ier_flag = 9
      end if
      if (nflxs .gt. nfloops) then
         print *, 'NFLXS>NFLOOPS'
         ier_flag = 9
      end if   
      if (nbfldn .gt. nbctotp) then
         print *, 'NBFLDN>NBCTOTP'
         ier_flag = 9
      end if   
      if (nbcoil_max .gt. nbcoilsp) then
         print *, 'NBCOIL_MAX>NBCOILSP'
         ier_flag = 9
      end if   

      if (ier_flag .eq. 10) return

      if (nobser+nobd .gt. 0) iconnect(:4,:nobser+nobd) = 0

      read(iunit) (xobser(n), zobser(n),n=1,nobser)
      read(iunit) (dsilabel(n),n=1,nobd)
      read(iunit) ((iconnect(j,n),j=1,4),n=1,nobd)

      if (nbcoil_max.gt.0 .and. nbsets.gt.0) then
         rbcoil(:nbcoil_max,:nbsets) = 0
         zbcoil(:nbcoil_max,:nbsets) = 0
         abcoil(:nbcoil_max,:nbsets) = 0

         do n=1,nbsets
           if (nbcoils(n).gt.0) then
           read(iunit) n1,bloopnames(n1)
           read(iunit)(rbcoil(m,n),zbcoil(m,n),abcoil(m,n),
     1             m=1,nbcoils(n))
           endif
         enddo

         dbcoil(:nbcoil_max,:nbsets,:nextcur) = 0
      end if
      do ii = 1,nextcur
        !un-connected coil fluxes
         read(iunit) (unpsiext(n,ii),n=1,nobser)
         do n = 1,nbsets
            read(iunit) (dbcoil(m,n,ii),m=1,nbcoils(n))
         enddo
      enddo

!
!     READ LIMITER & PROUT PLOTTING SPECS
!
      read(iunit,iostat=ierr2) nlim,(limitr(i),i=1,nlim)
      if (ierr2 .ne. 0)then
        nlim = 0
        if (lscreen) print *,' No limiter data in mgrid file'
        go to 900
      end if

      nlim_max = maxval(limitr)

      if (nlim .gt. nlimset) then
         print *, 'nlim>nlimset'
         ier_flag = 9
         return
      end if   
 
      allocate( rlim(nlim_max,nlim),   zlim(nlim_max,nlim),
     1          reslim(nlim_max,nlim) ,seplim(nlim_max,nlim),
     2          stat=istat)
      if (istat .ne. 0) then
         print *, 'rlim istat!=0'
         ier_flag = 9
         return
      end if   

      read(iunit, iostat=ierr2)
     1   ((rlim(i,j),zlim(i,j),i=1,limitr(j)),j=1,nlim)
      read(iunit, iostat=ierr2) nsets,(nsetsn(i), i=1,nsets)
      
      if (nsets .gt. nigroup) then
         print *, 'nsets>nigroup'
         ier_flag = 9
         return
      else if (ierr2 .ne. 0) then
         ier_flag = 9
         return
      end if

      nsets_max = maxval(nsetsn)
 
      if (nsets_max .gt. npfcoil) then
         print *, 'nsetsn>npfcoil'
         ier_flag = 9
         return
      end if   

      allocate (pfcspec(nparts,nsets_max,nsets), stat=istat)

!     NOTE TO RMW: SHOULD READ IN NPARTS HERE (PUT INTO MGRID FILE)

      read(iunit, iostat=ierr2) (((pfcspec(i,j,k),i=1,nparts),
     1        j=1,nsetsn(k)), k=1,nsets)
      read(iunit, iostat=ierr2) rx1,rx2,zy1,zy2,condif,
     1  nrgrid,nzgrid,tokid

      if (ierr2.ne.0 .or. istat.ne.0) then
         ier_flag = 9
         return
      end if
      
      if (nobser .gt. 0) xobsqr(:nobser) = sqrt(xobser(:nobser))
! 
!       PARTITION MGRID B-LOOPS INTO SETS
!
      nbcoilsn = sum(nbcoils(:nbsets))

      do n = 1,nbsets
        rbcoilsqr(:nbcoils(n),n) = sqrt(rbcoil(:nbcoils(n),n))
      enddo

 900  continue

      close (iunit)

      delrb = (rmaxb-rminb)/(nr0b-1)
      delzb = (zmaxb-zminb)/(nz0b-1)

!
!     SUM UP CONTRIBUTIONS FROM INDIVIDUAL COIL GROUPS
!
      if (lfreeb) then
         write (nthreed,20) nr0b, nz0b, np0b, rminb, rmaxb,
     1        zminb, zmaxb, trim(mgrid_file)
 20      format(//,' VACUUM FIELD PARAMETERS:',/,1x,24('-'),/,
     1  '  nr-grid  nz-grid  np-grid      rmin      rmax      zmin',
     2  '      zmax     input-file',/,3i9,4f10.3,5x,a)

         istat = 0
         if (.not.allocated(bvac)) allocate (bvac(nbvac,3), stat=istat)
         if (istat.ne.0) then
           print *,' bvac allocation failed'
           ier_flag = 9
           return
         end if

         bvac(:nbvac,:3) = 0
         if (nobser .gt. 0) psiext(:nobser) = 0
         if (nbcoil_max.gt.0 .and. nbsets.gt.0)
     1       bcoil(:nbcoil_max, :nbsets) = 0
 
         do ii = 1,nextcur
            do i = 1,3
               bvac(:nbvac, i) = bvac(:nbvac, i) + 
     1                           extcur(ii)*btemp(:nbvac,i,ii)
            enddo
            if (nobser .gt. 0)
     1      psiext(:nobser) = psiext(:nobser) +
     2                        extcur(ii)*unpsiext(:nobser,ii)
            do n=1,nbsets
               n1 = nbcoils(n)
               bcoil(:n1,n) = bcoil(:n1,n) + 
     1                        extcur(ii)*dbcoil(:n1,n,ii)
            enddo
         enddo
      endif                   !!IF LFREEB
      
      end subroutine read_mgrid


      subroutine restart(time_step)
      use vmec_main
      use xstuff
      implicit none
C-----------------------------------------------
C   D u m m y   V a r i a b l e s
C-----------------------------------------------
      real(rprec) :: time_step
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: c1p03 = 1.03_dp, cp90 = 0.90_dp
c-----------------------------------------------
 
      select case (irst) 
      case default
         xstore(:neqs2) = xc(:neqs2)
         return 
      case (2:3) 
         xcdot(:neqs2) = zero
         xc(:neqs2) = xstore(:neqs2)
         time_step = time_step*((irst-2)/c1p03 + cp90*(3-irst))
         if (irst .eq. 2) ijacob = ijacob + 1
         irst = 1
         return 
      end select

      end subroutine restart
      

      subroutine runvmec(input_file, iseq_count, lreseta, ier_flag, 
     1   lfirst_call, lscreen, reset_file_name)
      use vmec_main
      use vsvd
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer iseq_count, ier_flag
      logical lreseta, lfirst_call, lscreen
      character*(*) :: input_file, reset_file_name
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: bad_jacobian_flag = 6
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: igrid0, igrid, nsval, index_end, index_dat, 
     1   ns_min, ier_init, jacob_off
      real(rprec) :: timeon, timeoff
      logical :: interp_flag, lreset
C-----------------------------------------------
 
!
!
!     INDEX OF LOCAL VARIABLES
!
!     ier_flag   specifies error condition if nonzero
!     interp_flag 
!                = FALSE,compute xc array from scaled boundary data
!                = TRUE, compute xc array by interpolation
!     lfirst_call
!                = T, initial call to runvmec
!                = F, runvmec has been previously called
!
      call second0 (timeon)

!
!     PARSE input_file into path/input.ext
!
      index_dat = index(input_file,'input.')
      index_end = len_trim(input_file)
      if (index_dat .gt. 0) then
         input_extension  = input_file(index_dat+6:index_end)
      else
         input_extension = input_file(1:index_end)
         input_file = 'input.'//trim(input_extension)
      end if   

!
!     INITIALIZE PARAMETERS
!
      lreset = .false.
      ier_init = ier_flag
      if (lfirst_call .or. lreseta) lreset = .true.
      if (ier_init.ne.4 .and. ier_init.ne.bad_jacobian_flag) then
         call vsetup (lreset, iseq_count)
      else
         iequi = 0
         if (lfreeb) ivac = 1    !!Must restart vacuum calculations if free boundary
      end if   

!
!     READ INPUT FILES INDATA, MGRID_FILE
!     USE ISEQ_COUNT TO AVOID REREADING MGRID FILE
!
      call readin 
     1    (input_file, iseq_count, lfirst_call, ier_flag, lscreen)
      if (ier_flag .ne. 0) goto 1000      
 
!
!     COMPUTE INVARIANT ARRAYS 
!
      call fixaray
    
      if (ier_init .ne. 4) then
         write (nthreed, 230)
  230 format(' FSQR, FSQZ = Normalized Physical Force Residuals',/,
     1   ' fsqr, fsqz = Preconditioned Force Residuals',/,1x,23('-'),/,
     2   ' BEGIN FORCE ITERATIONS',/,1x,23('-'),/)
      end if

!
!     COMPUTE INITIAL SOLUTION ON COARSE GRID
!     IF PREVIOUS SEQUENCE DID NOT CONVERGE, DO COARSE RESTART
!

      if (lreseta) then
        igrid0 = 1
      else
        igrid0 = multi_ns_grid
      endif


      imovephi = 0
      ns_min = 0
      jacob_off = 0
      ier_flag = ier_init
      if (ier_flag .eq. bad_jacobian_flag) jacob_off = 1
      if (all(ns_array .eq. 0)) then
         ier_flag = 8
         goto 1000
      end if
         
      do igrid = igrid0, multi_ns_grid + jacob_off
         if (jacob_off.eq.1 .and. igrid.eq.igrid0) then
!           TRY TO GET NON-ZERO JACOBIAN ON A 3 PT RADIAL MESH         
            nsval = 3
            ftolv = 1.e-7_dp
         else
            nsval = ns_array(igrid-jacob_off)
            if (nsval .le. ns_min) cycle
            ns_min = nsval
            ftolv = ftol_array(igrid-jacob_off)
         end if
         if (igrid .eq. igrid0) then
            interp_flag = .false.
         else
            interp_flag = .true.
         endif
         call eqsolve (nsval, interp_flag, ier_flag, lreseta, 
     1                 lfirst_call, lscreen, reset_file_name)
         if (imatch_phiedge .eq. 3) imovephi = 1
         if (ier_flag .ne. 0 .and. ier_flag .ne. 4) exit
      end do
 
  100 continue
 
      call second0 (timeoff)
      timer(0) = timer(0) + timeoff - timeon
 
!
!     WRITE OUTPUT TO THREED1, WOUT FILES; FREE MEMORY ALLOCATED GLOBALLY
!
 1000 call fileout (iseq_count, ier_flag, lscreen) 
      
      end subroutine runvmec
      
      subroutine vsetup (lreset, iseq_count)
      use vmec_main
      use vacmod
      use realspace
      use vsvd
      use vspline
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: iseq_count
      logical :: lreset
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: mpol_default = 6
      integer, parameter :: ntor_default = 0
C----------------------------------------------- 
!
!     Default values
!
      lasym = .false.
      mpol = mpol_default
      ntor = ntor_default
      ntheta = 0
      nzeta  = 0
      lfreeb = .true.              !reset later after input file read
      lrecon = .true.
      loldout = .false.
      ldiagno = .false.
      lmac = .false.
      ledge_dump = .false.

      z00 = zero
      mgrid_file = 'NONE'

!
!     ZERO ARRAYS WHICH MAY BE READ IN
!
      rbc = zero
      rbs = zero
      zbc = zero
      zbs = zero
 
      am = zero
      ai = zero
      ac = cbig
      ac_form = 0
      ac_form = 0
      aphi = zero
      extcur = zero
      bloat = one

      ns_array = 0
      nsin = 0
      ftol_array = zero
      raxis = zero
      zaxis = zero
      gamma = zero
      spres_ped = one
       
      iequi = 0
      ivac  = -1
      delbsq = one
      delt = 1.1
      tcon0 = one
      curtor = 1.e30_dp
      time_slice = zero
      if (lreset) then
         pfac   = one
         phifac = one
         timer = zero
      endif
      fsqr = one
      fsqz = one
      ftolv = fsqr
!
!     Reconstruction stuff
!
      lpofr = .true.
      lpprof = .true.
      icurrout = 0
      total_chi_square_n = one
      total_chisq_n0 = total_chi_square_n
      nchistp = 0
      imse = -1
      itse = 0
      isnodes = 0
      ipnodes = 0
      iopt_raxis = 1
      imatch_phiedge = 1
      nflxs = 0
      nbfld = 0
      mseangle_offset = zero
      mseangle_offsetm = zero
      pres_offset = zero
      sigma_current = 1.e30_dp
      sigma_delphid = 1.e30_dp
      tensi = one
      tensp = one
      tensi2 = zero
      fpolyi = one
      presfac = one
      phidiam = 1.e30_dp
 
      mseprof = one
      indxflx = 0
      indxbfld = 0
      sigma_stark = 1.1*cbig
      sigma_thom = 1.1*cbig
      sigma_flux = 1.1*cbig
      sigma_b = 1.1*cbig
!
!     FREE-BOUNDARY STUFF, READ IN FIRST TIME ONLY
!
      if (iseq_count .eq. 0) then
        nbcoil_max = 0
        nlim_max = 0
      end if

      end subroutine vsetup

      
      subroutine wroutlim(nwout)
      use vmec_main
      use vsvd
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: nwout
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer l1, l2, ltouch, limpts, n, i
      integer, dimension(nlim_max*nlim) :: l1min,l2min
      real(rprec) :: sepmin
      character, dimension(nlim_max,nlim) :: closest
C-----------------------------------------------
 
 
      if (nlim .le. 0) return 
 
      do l1 = 1, nlim
         do l2 = 1, limitr(l1)
            closest(l2,l1) = ' '
            seplim(l2,l1) = sqrt(abs(seplim(l2,l1)))
         end do
      end do

      sepmin = seplim(1,1)
      do l1 = 1, nlim
         do l2 = 1, limitr(l1)
            if (seplim(l2,l1) .lt. sepmin)
     1      sepmin = seplim(l2,l1)
         end do
      end do

!     Do it this way to catch multiple mins
      ltouch = 0
      do l1 = 1, nlim
         do l2 = 1, limitr(l1)
            if (abs(seplim(l2,l1)-sepmin) .lt. 1.e-5_dp) then
               closest(l2,l1) = '*'
               ltouch = ltouch + 1
               l2min(ltouch) = l2
               l1min(ltouch) = l1
            endif
         end do
      end do
      write (nwout, 702) ltouch
      write (nwout, 704) (l1min(i),l2min(i),rlim(l2min(i),l1min(i)),
     1   zlim(l2min(i),l1min(i)),seplim(l2min(i),l1min(i)),i=1,ltouch)
  702 format(8i10)
  704 format(2i5,3e20.13)
      if (lrecon) then
      write (nthreed, 705)
  705 format(/,' PLASMA BOUNDARY-LIMITER PROXIMITY'/
     1   ' POINTS OUTSIDE PLASMA HAVE RESIDUE < 0.5')
      write (nthreed, 707)
  707 format(/,13x,' nlim    n       R         Z        Residue',
     1   '          min |d|        nearest'/13x,' ----',4x,'-',7x,'-',9x
     2   ,'-',8x,'-------',10x,'-------',8x,'-------')
 
      do n = 1, nlim
         limpts = limitr(n)
         do i = 1, limpts
            write (nthreed, 708) n, i, rlim(i,n), zlim(i,n), reslim(i,n)
     1         , seplim(i,n), closest(i,n)
         end do
      end do
  708 format(13x,i5,i5,2f10.3,1pe15.4,5x,1pe12.4,8x,a)
      end if

      end subroutine wroutlim

EOF

cat > vvacuum.f << "EOF"
!     THIS ROUTINE COMPUTES .5 * B**2 ON THE VACUUM / PLASMA SURFACE
!     BASED ON THE PROGRAM BY P. MERKEL [J. Comp. Phys. 66, 83 (1986)]
!     AND MODIFIED BY W. I. VAN RIJ AND S. P. HIRSHMAN (1987)
 
!     THE USER MUST SUPPLY THE FILE << MGRID >> WHICH INCLUDES THE MAGNETIC
!     FIELD DATA TO BE READ BY THE SUBROUTINE BECOIL
!     THE "VACUUM.INC" FILE IS DEFINED IN VMEC.UNIX
!
      subroutine vacuum(rmnc, rmns, zmns, zmnc, xm, xn, rcurr, zcurr, 
     1   plascur, rbtor, wint, ns, ivac_skip, ivac, mnmax, ier_flag, 
     2   lscreen)
      use vacmod
      use vparams, only: nthreed, zero, one, dmu0
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: ns, ivac_skip, ivac, mnmax, ier_flag
      real(rprec) :: plascur, rbtor
      real(rprec), dimension(mnmax), intent(in) :: 
     1   rmnc, rmns, zmns, zmnc, xm, xn
      real(rprec), dimension(nv), intent(in) :: rcurr, zcurr
      real(rprec), dimension(*), intent(in) :: wint
      logical :: lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: mn, n, n1, m, i, istore_max
      real(rprec), dimension(:), pointer :: potcos, potsin
      real(rprec), allocatable :: bsubu(:), bsubv(:), potu(:), potv(:)
      real(rprec), allocatable :: amatrix(:)
      real(rprec):: dn2, dm2, cosmn, sinmn, huv, hvv,
     1    det, bsupu, bsupv, bsubuvac, fac
C-----------------------------------------------
      if (.not.allocated(potvac)) stop 'POTVAC not allocated in VACCUM'

      allocate (amatrix(mnpd2*mnpd2), bsubu(nuv2), bsubv(nuv2),
     1    potu(nuv2), potv(nuv2), stat = i)
      if (i .ne. 0) stop 'Allocation error in vacuum'      

      potsin => potvac(1:mnpd)
      potcos => potvac(1+mnpd:)

      allocate (bexu(nuv2), bexv(nuv2), bexn(nuv2),
     1     bexni(nuv2), r1b(nuv), rub(nuv2), rvb(nuv2),
     2     z1b(nuv), zub(nuv2), zvb(nuv2), auu(nuv2), auv(nuv2), 
     3     avv(nuv2), snr(nuv2), snv(nuv2), snz(nuv2), drv(nuv2),
     4     guu_b(nuv2), guv_b(nuv2), gvv_b(nuv2), rb2(nuv),
     5     rcosuv(nuv), rsinuv(nuv), stat=i)
      if (i .ne. 0) stop 'Allocation error in vacuum'

!
!       INDEX OF LOCAL VARIABLES
!
!       rmnc,rmns,zmns,zmnc:     Surface Fourier coefficients (m,n) of R,Z
!       xm,xn:     m, n values corresponding to rc,zs array
!       rcurr,zcurr:
!                  Position of magnetic axis (toroidal current filament)
!       bsqvac:    B**2/2 at the vacuum interface
!       plascur:   net toroidal current
!       mnmax:     number of R, Z modes in Fourier series of R,Z
!       ivac_skip: regulates whether full (=0) or incremental (>0)
!                 update of matrix elements is necessary
!
!
!       compute and store mean magnetic fields (due to
!       toroidal plasma current and external tf-coils)
!       note: these are fixed for a constant current iteration
!       bfield = rbtor*grad(zeta) + plascur*grad("theta") - grad(potential)
!            where "theta" is computed using Biot-Savart law for filaments
!
!       Here, the potential term is needed to satisfy B ! dS = 0 and has the form:
!
!       potential = sum potsin*sin(mu - nv) + potcos*cos(mu - nv)
!
 
      call surface (rmnc, rmns, zmns, zmnc, xm, xn, mnmax)
      call bextern (rcurr, zcurr, plascur, rbtor, wint, ns)
!
!       Determine scalar magnetic potential
!
      istore_max = min(64,nuv2)
      call scalpot (potvac, amatrix, wint, ns, istore_max, ivac_skip)
      call solver (amatrix, potvac, mnpd2)
!
!       compute tangential covariant (sub u,v) and contravariant
!       (super u,v) magnetic field components on the plasma surface
!
      potu(:nuv2) = zero;  potv(:nuv2) = zero
 
      mn = 0
      do n = -nf, nf
         dn2 = -(n*nfper)
         n1 = abs(n)
         do m = 0, mf
            mn = mn + 1
            dm2 = (m)
            do i = 1, nuv2
               cosmn = cosu1(i,m)*cosv1(i,n1) + csign(n)*sinu1(i,m)*
     1            sinv1(i,n1)
               potu(i) = potu(i) + dm2*potsin(mn)*cosmn 
               potv(i) = potv(i) + dn2*potsin(mn)*cosmn 
            end do
            if (.not.lasym) cycle
            do i = 1, nuv2
               sinmn = sinu1(i,m)*cosv1(i,n1) - csign(n)*cosu1(i,m)*
     1            sinv1(i,n1)
               potu(i) = potu(i) - dm2*potcos(mn)*sinmn
               potv(i) = potv(i) - dn2*potcos(mn)*sinmn
            end do
         end do
      end do
      do i = 1, nuv2
         bsubu(i) = potu(i) + bexu(i)
         bsubv(i) = potv(i) + bexv(i)
         huv = 0.5_dp*guv_b(i)*(nfper)
         hvv = gvv_b(i)*(nfper*nfper)
         det = one/(guu_b(i)*hvv-huv*huv)
         bsupu = (hvv*bsubu(i)-huv*bsubv(i))*det
         bsupv = ((-huv*bsubu(i))+guu_b(i)*bsubv(i))*det
         bpolvac(i) = 0.5_dp*bsubu(i)*bsupu
         bsqvac(i) = bpolvac(i) + 0.5_dp*bsubv(i)*bsupv
         brv(i) = rub(i)*bsupu + rvb(i)*bsupv
         bphiv(i) = r1b(i)*bsupv
         bzv(i) = zub(i)*bsupu + zvb(i)*bsupv
      end do

!
!       PRINT OUT VACUUM PARAMETERS
!
      if (ivac .eq. 0) then
         ivac = ivac + 1
         if (lscreen) write (*, 200) nfper, mf, nf, nu, nv
         write (nthreed, 200) nfper, mf, nf, nu, nv
  200    format(/,2x,'In VACUUM, np =',i3,2x,'mf =',i3,2x,'nf =',i3,
     1      ' nu =',i3,2x,'nv = ',i4)
         bsubuvac = sum(bsubu(:nuv2)*wint(ns:ns*nuv2:ns))
         bsubvvac = sum(bsubv(:nuv2)*wint(ns:ns*nuv2:ns))
         fac = 1.e-6_dp/dmu0
         if (lscreen )write (*,1000) (-pi2*bsubuvac*fac),
     1       plascur*fac, bsubvvac, rbtor
         write (nthreed, 1000) (-pi2*bsubuvac*fac), plascur*fac, 
     1       bsubvvac, rbtor
 1000    format(2x,'2*pi * a * -BPOL(vac) = ',1pe10.2,
     1      ' TOROIDAL CURRENT = ',1pe10.2,/,2x,'R * BTOR(vac) = ',
     2      1pe10.2,' R-BTOR = ',1pe10.2)
         if (rbtor*bsubvvac .lt. zero) ier_flag = 7
         if (abs((plascur+pi2*bsubuvac)/rbtor) .gt. 1.e-3_dp)
     1      ier_flag = 10
      endif

      if (allocated(bexu))
     1    deallocate (bexu, bexv, bexn, bexni, r1b, rub, rvb, z1b, zub, 
     2    zvb, auu, auv, avv, snr, snv, snz, drv, guu_b, guv_b, gvv_b, 
     3    rb2, rcosuv, rsinuv, stat=i)
      if (i .ne. 0) stop 'Deallocation error in vacuum'
       
      deallocate (amatrix, bsubu, bsubv, potu, potv, stat = i)

      end subroutine vacuum

      
      subroutine analysum(grpmn, bvec, sl, tl, m, n, l, ivacskip)
      use vacmod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer m, n, l, ivacskip
      real(rprec), intent(inout) :: grpmn(nuv2,0:mf,-nf:nf,*)
      real(rprec), intent(inout) :: bvec(0:mf,-nf:nf,*)
      real(rprec), dimension(nuv2), intent(in) :: sl, tl
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat
      real(rprec), allocatable :: sinp(:), cosp(:)
C-----------------------------------------------
      allocate (sinp(nuv2), stat=istat)
      if (istat .ne. 0) stop 'Allocation error in analysum'
      
      if (n .ne. abs(n)) stop 'error calling analysum!'
 
      sinp = (sinu1(:,m)*cosv1(:,n) - sinv1(:,n)*
     1            cosu1(:,m))*cmns(l,m,n)
      bvec(m,n,1) = bvec(m,n,1) + sum(tl*bexni*sinp)

      if (ivacskip .eq. 0) grpmn(:,m,n,1) = grpmn(:,m,n,1) + sl*sinp
      deallocate (sinp)

      if (.not.lasym) return

      allocate (sinp(nuv2), stat=istat)
      if (istat .ne. 0) stop 'Allocation error in analysum'
      cosp = (cosu1(:,m)*cosv1(:,n) + cosv1(:,n)*
     1        cosu1(:,m))*cmns(l,m,n)
      bvec(m,n,2) = bvec(m,n,2) + sum(tl*bexni*cosp)
      if (ivacskip .eq. 0) grpmn(:,m,n,2) = grpmn(:,m,n,2) + sl*cosp
      deallocate (cosp)
 
      end subroutine analysum

      
      subroutine analysum2(grpmn, bvec, slp, tlp, slm, tlm, 
     1    m, n, l, ivacskip)
      use vacmod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer m, n, l, ivacskip
      real(rprec), intent(inout) :: grpmn(nuv2,0:mf,-nf:nf,*)
      real(rprec), intent(inout) :: bvec(0:mf,-nf:nf,*)
      real(rprec), dimension(nuv2), intent(in) :: 
     1   slp, tlp, slm, tlm
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat
      real(rprec), allocatable :: sinp(:), sinm(:), temp(:)
      real(rprec), allocatable :: cosp(:), cosm(:)
C-----------------------------------------------
      allocate (sinp(nuv2), sinm(nuv2), temp(nuv2), stat=istat)
      if (istat .ne. 0) stop 'Allocation error in analysum2'
      if (n .ne. abs(n)) stop 'error calling analysum2!'
 
      sinp =  sinu1(:,m)*cosv1(:,n)*cmns(l,m,n)
      temp = -sinv1(:,n)*cosu1(:,m)*cmns(l,m,n)
      sinm = sinp - temp                !sin(mu + |n|v) * cmns (l,m,|n|)
      sinp = sinp + temp                !sin(mu - |n|v) * cmns
      bvec(m,n,1)  = bvec(m,n,1)  + sum(tlp*bexni*sinp)
      bvec(m,-n,1) = bvec(m,-n,1) + sum(tlm*bexni*sinm)

      if (ivacskip .eq. 0) then 
         grpmn(:,m,n,1)  = grpmn(:,m,n,1)  + slp*sinp
         grpmn(:,m,-n,1) = grpmn(:,m,-n,1) + slm*sinm
      end if   

      deallocate (sinp, sinm, temp)

      if (.not.lasym) return
      
      allocate (cosp(nuv2), cosm(nuv2), stat=istat)  
      if (istat .ne. 0) stop 'Allocation error in analysum2'

      cosp = cosu1(:,m)*cosv1(:,n)*cmns(l,m,n)     
      temp = cosv1(:,n)*cosu1(:,m)*cmns(l,m,n)
      cosm = cosp - temp                !cos(mu + |n|v) * cmns (l,m,|n|) 
      cosp = cosp + temp                !cos(mu - |n|v) * cmns (l,m,|n|) 
      bvec(m,n,2)  = bvec(m,n,2)  + sum(tlp*bexni*cosp)
      bvec(m,-n,2) = bvec(m,-n,2) + sum(tlm*bexni*cosm)

      if (ivacskip .eq. 0) then 
         grpmn(:,m,n,2)  = grpmn(:,m,n,2)  + slp*cosp
         grpmn(:,m,-n,2) = grpmn(:,m,-n,2) + slm*cosm
      end if   

      deallocate (cosp, cosm)
 
      end subroutine analysum2


      subroutine analyt(grpmn, bvec, ivacskip)
      use vacmod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ivacskip
      real(rprec), intent(inout) :: grpmn(nuv2,0:mf,-nf:nf,*)
      real(rprec), intent(inout) :: bvec(0:mf,-nf:nf,*)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero=0, one=1, two=2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: l, n, m
      real(rprec), dimension(:), allocatable :: 
     1   r0p, r1p, r0m, r1m, sqrtc, sqrta, tlp2, tlp1, tlp, tlm2,
     2    tlm1, tlm, adp, adm, cma, ra1p, ra1m, slm, slp, tlpm, slpm,
     3    delt1u, azp1u, azm1u, cma11u, sqad1u, sqad2u
      real(rprec) :: fl, sign1
C-----------------------------------------------
      allocate (r0p(nuv2), r1p(nuv2), r0m(nuv2), r1m(nuv2), 
     1          sqrtc(nuv2), sqrta(nuv2), tlp2(nuv2), tlp1(nuv2), 
     2          tlp(nuv2), tlm2(nuv2), tlm1(nuv2), tlm(nuv2), adp(nuv2),
     3          adm(nuv2), cma(nuv2), ra1p(nuv2), ra1m(nuv2), slm(nuv2), 
     4          slp(nuv2), tlpm(nuv2), slpm(nuv2), delt1u(nuv2), 
     5          azp1u(nuv2), azm1u(nuv2), cma11u(nuv2), sqad1u(nuv2), 
     6          sqad2u(nuv2), stat = l)
      if (l .ne. 0) stop 'Allocation error in analyt'

!
!     ALL EQUATIONS REFER TO THE APPENDIX OF THE PAPER BY P. MERKEL
!     IN J. COMPUT. PHYSICS 66, p83 (1986)
!
!     THE REQUIRED INTEGRALS ARE:
!
!     BVECS(m,n) = Int< sin(mu - nv) Int<BNORM`*Gan(u`,v`,u-u`,v-v`)>
!     BVECC(m,n) = Int< cos(mu - nv) Int<BNORM`*Gan(u`,v`,u-u`,v-v`)>
!
!     Where Int<...> means integration of u (theta) and v (zeta) and
!     summation over field periods. In terms of Merkels Imn integral,
!     a = guu (g theta-theta), etc., we have
!
!     BVECS(m,n) = (2*pi/nfp) * Int<BNORM *sin(mu` - nv`)[Im,-n](a,b,c)>
!     BVECC(m,n) = (2*pi/nfp) * Int<BNORM *cos(mu` - nv`)[Im,-n](a,b,c)>
!
!     Similarly, the analytic part of the matrix A(m,n;m`,n`) can be written:
!
!     A(m,n;m`,n`) = (2*pi/nfp) * Int<sin(mu` - nv`)*sin(m`u` - n`v`)
!                              [Km,-n](a`,b`,c`;A`,B`,C`)>
!
!     On exit, GRPMN(ip,m,n) = ALP * SIN(ip,m,n) * K[m,-n](ip)
!
!
!     COMPUTE ALL QUANTITIES INDEPENDENT OF THE MODE INDICES L,M,N
!
!     ADP(M): a +(-)2b + c
!     CMA:    c - a
!     DELTA:  4*(ac - b**2)
!     AZP(M): A +(-)2*B + C
!     CMA1:   C - A
!     R1P(M): Coefficient of l*Tl+(-) in eq (A17)
!     R0P(M): Coefficient of l*T(l-1)+(-) in eq (A17)
!     RA1P(M):Coefficient of Tl+(-) in eq (A17)
!
      adp  = guu_b  + guv_b  + gvv_b 
      adm  = guu_b  - guv_b  + gvv_b 
      cma  = gvv_b  - guu_b 
      sqrtc  = two*sqrt(gvv_b)
      sqrta  = two*sqrt(guu_b)
 
      if (ivacskip .eq. 0) then
         delt1u  = adp *adm  - cma *cma 
         azp1u  = auu  + auv  + avv 
         azm1u  = auu  - auv  + avv 
         cma11u  = avv  - auu 
         r1p  = (azp1u *(delt1u  - (cma *cma ))/adp 
     1        - (azm1u *adp ) + ((two*cma11u )*cma ))/delt1u 
         r1m  = (azm1u *(delt1u  - (cma *cma ))/adm 
     1        - (azp1u *adm ) + ((two*cma11u )*cma ))/delt1u 
         r0p  = ((-(azp1u *adm )*cma /adp ) - azm1u *cma 
     1        + (two*cma11u )*adm )/delt1u 
         r0m  = ((-(azm1u *adp )*cma /adm ) - azp1u *cma 
     1        + (two*cma11u )*adp )/delt1u 
         ra1p = azp1u /adp 
         ra1m = azm1u /adm 
      endif
!
!     INITIALIZE T0+ and T0-
!
!     TLP(M): TL+(-)
!     TLP(M)1:T(L-1)+(-)
!     TLP(M)2:T(L-2)+(-)
!
      sqad1u  = sqrt(adp )
      sqad2u  = sqrt(adm )
      tlp1 = zero
      tlm1 = zero
      tlp  = one/sqad1u *log((sqad1u *sqrtc  + adp 
     1     + cma )/(sqad1u *sqrta  - adp  + cma ))
      tlm  = one/sqad2u *log((sqad2u *sqrtc  + adm  
     1     + cma )/(sqad2u *sqrta  - adm  + cma ))
      tlpm = tlp  + tlm 
!
!     BEGIN L-SUM IN EQ (A14) TO COMPUTE Cmn COEFFICIENTS
!
      do l = 0, mf + nf
         fl = l
         sign1 = one - two*mod(l,2)          !!(-1)**l
!
!       COMPUTE SL+ and SL- , Eq (A17)
!       SLP(M): SL+(-)
!
         if (ivacskip .eq. 0) then
            slp = (r1p*fl + ra1p)*tlp + r0p*fl*tlp1 - (r1p + r0p)/sqrtc
     1          + sign1*(r0p - r1p)/sqrta
            slm = (r1m*fl + ra1m)*tlm + r0m*fl*tlm1 - (r1m + r0m)/sqrtc
     1          + sign1*(r0m - r1m)/sqrta
            slpm = slp + slm
         endif
!
!       BEGIN MODE NUMBER (m,n) LOOP
!
         do n = 0, nf
            do m = 0, mf
               if (cmns(l,m,n) .eq. zero) cycle 

               if (n.eq.0 .or. m.eq.0) then
!
!       1. n = 0 and  m >= 0  OR n > 0 and m = 0
!
                 call analysum (grpmn, bvec, slpm, tlpm, m, n, l, 
     1               ivacskip)
 
               else
!
!       2. n>=1  and  m>=1
!
                 call analysum2 (grpmn, bvec, slm, tlm, slp, tlp, 
     1               m, n, l, ivacskip)
               endif
            end do
         end do
!
!       UPDATE "TL's" (FOR L -> L+1) USING EQ (A15)
!
         tlp2 = tlp1
         tlm2 = tlm1
         tlp1 = tlp
         tlm1 = tlm
         tlp = ((sqrtc - (sign1*sqrta)) - (((two*fl) + one)*cma)*tlp1 - 
     1      fl*adm*tlp2)/(adp*(fl + one))
         tlm = ((sqrtc - (sign1*sqrta)) - (((two*fl) + one)*cma)*tlm1 - 
     1      fl*adp*tlm2)/(adm*(fl + one))
         tlpm = tlp + tlm
      end do
 
      deallocate (r0p, r1p, r0m, r1m, sqrtc, sqrta, tlp2, tlp1, 
     1          tlp, tlm2, tlm1, tlm, adp, adm, cma, ra1p, ra1m, slm, 
     2          slp, tlpm, slpm, delt1u, azp1u, azm1u, cma11u, sqad1u,
     3          sqad2u, stat = l)

      end subroutine analyt

      
      subroutine becoil(rad, zee, br, bp, bz, brvac, bpvac, bzvac)
      use vparams, only: one, nthreed
      use vacmod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(*), intent(in) :: rad, zee
      real(rprec), dimension(*), intent(out) :: br, bp, bz
      real(rprec), dimension(nr0b,nz0b,nv), intent(in) :: 
     1   brvac, bpvac, bzvac
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      character*(50), parameter :: warning =
     1   'Plasma Boundary exceeded Vacuum Grid Size'
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, save :: icount = 0
      integer :: igrid, i, kv, ir, jz, ir1, jz1
      real(rprec) :: rad0, zee0, ri, zj, 
     1   dri, dzj, f22, f21, f12, f11
C-----------------------------------------------
!
!       THIS ROUTINE FINDS THE CYLINDRICAL COMPONENTS OF THE EXTERNAL
!       MAGNETIC FIELD AT A FIXED PHI PLANE BY 2-D INTERPOLATION
!
      igrid = 0
      icount = icount + 1

      do i = 1, nuv2
!
!       CHECK THAT BOUNDARY POINTS ARE INSIDE VACUUM GRID.  IF NOT,
!       SET THEM EQUAL TO LARGEST (OR SMALLEST) VACUUM GRID POINTS
!
         if (rad(i) .gt. rmaxb) then
            rad0 = rmaxb
            igrid = 1
         else if (rad(i) .lt. rminb) then
            igrid = 1
            rad0 = rminb
         else
            rad0 = rad(i)
         endif
         if (zee(i) .gt. zmaxb) then
            igrid = 1
            zee0 = zmaxb
         else if (zee(i) .lt. zminb) then
            igrid = 1
            zee0 = zminb
         else
            zee0 = zee(i)
         endif
!
!       DETERMINE PHI-PLANE, KV (MUST LIE IN FIRST FIELD PERIOD)
!
         kv = 1 + mod(i - 1,nv)
!
!
!       DETERMINE R, Z CORNER INDICES (IR,IZ)
!
         ir = int((rad0 - rminb)/delrb) + 1
         jz = int((zee0 - zminb)/delzb) + 1
         ir1 = min0(nr0b,ir + 1)
         jz1 = min0(nz0b,jz + 1)
!
!       COMPUTE R , Z AND DR , DZ AT MESH POINT (IR , IZ)
!
         ri = rminb + (ir - 1)*delrb
         zj = zminb + (jz - 1)*delzb
         dri = (rad0 - ri)/delrb
         dzj = (zee0 - zj)/delzb
         f22 = dri*dzj
         f21 = dri - f22
         f12 = dzj - f22
         f11 = one + f22 - (dri + dzj)
!
!       COMPUTE INTERPOLATED B FIELD
!
         br(i) = f11*brvac(ir,jz,kv) + f22*brvac(ir1,jz1,kv) + 
     1      f21*brvac(ir1,jz,kv) + f12*brvac(ir,jz1,kv)
         bz(i) = f11*bzvac(ir,jz,kv) + f22*bzvac(ir1,jz1,kv) + 
     1      f21*bzvac(ir1,jz,kv) + f12*bzvac(ir,jz1,kv)
         bp(i) = f11*bpvac(ir,jz,kv) + f22*bpvac(ir1,jz1,kv) + 
     1      f21*bpvac(ir1,jz,kv) + f12*bpvac(ir,jz1,kv)
 
      end do
 
      if (igrid.eq.1 .and. mod(icount,25).eq.0) then
         print *, warning
         write (nthreed, *) warning
      endif
 
      end subroutine becoil
      

      subroutine belicu(bx, by, bz, cos1, sin1, rp, zp, torcur)
      use vacmod
      use vacwires
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), intent(in) :: torcur
      real(rprec), dimension(*), intent(in)  :: cos1, sin1, rp, zp
      real(rprec), dimension(*), intent(out) :: bx, by, bz
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
      integer :: i, j
      real(rprec), allocatable, dimension(:) :: x0, y0, z0, rw, fa
      real(rprec) :: sumx, sumy, sumz, ax, ay, az, curpi, 
     1   xp, yp, r12
C-----------------------------------------------
 
      curpi = torcur/pi4
      if (nvp .ne. nv*nfper) stop 'nvp != nv*nfper in belicu'
      allocate (x0(nvp+1), y0(nvp+1), z0(nvp+1), rw(nvp+1), fa(nvp+1))

      do j = 1,nuv2
        xp   = rp(j) * cos1(j)
        yp   = rp(j) * sin1(j)

        x0   = xp - xw
        y0   = yp - yw
        z0   = zp(j) - zw
        rw   = sqrt(x0*x0 + y0*y0 + z0*z0 )

        do 20 i = 1,nvp
          r12    = rw(i+1)*rw(i)
          fa(i)  = (rw(i+1)+rw(i))/
     1            (r12*(r12 +x0(i+1)*x0(i)+y0(i+1)*y0(i)+z0(i+1)*z0(i)))
 20     continue

        ax = 0;    ay = 0;      az = 0
        sumx = 0;  sumy = 0;    sumz = 0

        do 30 i = 1,nvp
          ax = ax + fa(i)*dx(i)
          ay = ay + fa(i)*dy(i)
          az = az + fa(i)*dz(i)
          sumx = sumx + fa(i)*vx(i)
          sumy = sumy + fa(i)*vy(i)
          sumz = sumz + fa(i)*vz(i)
 30     continue
        bx(j)  = curpi * (sumx - yp*az    + zp(j)*ay)
        by(j)  = curpi * (sumy - zp(j)*ax + xp*az)
        bz(j)  = curpi * (sumz - xp*ay    + yp*ax)
      enddo

      deallocate (x0, y0, z0, rw, fa)

      end subroutine belicu

      
      subroutine bextern(rcurr, zcurr, plascur, rbtor, wint, ns)
      use vacmod
      use vacwires
      use mgrid_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: ns
      real(rprec), intent(in) :: plascur, rbtor
      real(rprec), dimension(nv), intent(in) :: rcurr, zcurr
      real(rprec), dimension(*), intent(in) :: wint
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i
      real(rprec), allocatable :: brad(:), bphi(:), bz(:)
C-----------------------------------------------
      if (.not.allocated(bvac)) stop 'BVAC Unallocated in bextern'
      allocate (brad(nuv2), bphi(nuv2), bz(nuv2), stat=i)
      if (i .ne. 0) stop 'allocation error in bextern'

!
!     THIS ROUTINE COMPUTES THE EXTERNAL, COIL PRODUCED B DOT DS
!     NOTE THAT BEXN = - BEX * DS IS THE EFFECTIVE SOURCE TERM
!
      call becoil(r1b,z1b,brad,bphi,bz,bvac(1,1),bvac(1,3),bvac(1,2))

!
!     COMPUTE CONTRIBUTION FROM NET TOROIDAL PLASMA CURRENT
!
      if (abs(plascur) > 1.e-8_dp*abs(rbtor)) then
         allocate (xw(nvp+1),yw(nvp+1),zw(nvp+1),vx(nvp+1),
     1     vy(nvp+1),vz(nvp+1),dx(nvp+1),dy(nvp+1),dz(nvp+1),stat=i)
         if (i .ne. 0) stop 'allocation error in bextern subroutine'
         call tolicu (rcurr, zcurr)
         call belicu (bexu, bexv, bexn, cosuv, sinuv, r1b, z1b, plascur)
         deallocate (xw, yw, zw, vx, vy, vz, dx, dy, dz)
         do i = 1, nuv2
            brad(i) = brad(i) + bexu(i)*cosuv(i) + bexv(i)*sinuv(i)
            bphi(i) = bphi(i) - bexu(i)*sinuv(i) + bexv(i)*cosuv(i)
            bz(i) = bz(i) + bexn(i)
         end do
      endif
 
      do i = 1, nuv2
        bexu(i) = rub(i)*brad(i) + zub(i)*bz(i)
        bexv(i) = rvb(i)*brad(i) + zvb(i)*bz(i) + r1b(i)*bphi(i)
        bexn(i) = -(brad(i)*snr(i)+bphi(i)*snv(i)+bz(i)*snz(i))
      end do

      deallocate (brad, bphi, bz)
      bexni(:nuv2) = wint(ns:nuv2*ns:ns)*bexn*pi2*pi2
       
      end subroutine bextern
      

      subroutine fouri(grpmn, gsource, amatrix, amatsq, bvec, wint, ns)
      use vacmod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ns
      real(rprec), dimension(nv,nu3,mnpd,*), intent(in) :: grpmn
      real(rprec), dimension(nuv), intent(in) :: gsource
      real(rprec), dimension(mnpd,mnpd,*) :: amatrix
      real(rprec), dimension(mnpd2,mnpd2), intent(out) :: amatsq
      real(rprec), dimension(0:mf,-nf:nf,*) :: bvec
      real(rprec), dimension(*) :: wint
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: k, i, j, n, kvi, kui, mn, m
      real(rprec), allocatable, dimension(:,:,:) :: bcos, bsin, source
      real(rprec), allocatable :: actemp(:,:,:,:), astemp(:,:,:,:)
      real(rprec) :: cosn, sinn, cosm, sinm
C-----------------------------------------------
!
!       AMATRIX(,1) = A(Sin)(Sin);  AMATRIX(,2) = A(Sin)(Cos);
!       AMATRIX(,3) = A(Cos)(Sin);  AMATRIX(,4) = A(Cos)(Cos)
!
!       SYMMETRIZE SOURCE TERMS
!
      allocate (bcos(nu2,-nf:nf,2), bsin(nu2,-nf:nf,2),
     1   actemp(nu3,-nf:nf,mnpd,2), astemp(nu3,-nf:nf,mnpd,2),
     2   source(nv,nu2,2), stat = i)
      if (i .ne. 0) stop 'allocation error in fouri'

      k = 0
      do i = 1, nu2
         do j = 1, nv
            k = k + 1
            source(j,i,1) = 0.5_dp*onp*(gsource(k)-gsource(imirr(k)))
            if (lasym)
     1      source(j,i,2) = 0.5_dp*onp*(gsource(k)+gsource(imirr(k)))
         end do
      end do
!
!       INITIALIZE SUMMED VECTORS TO ZERO
!
      bcos(:,(-nf):nf,:) = 0
      bsin(:,(-nf):nf,:) = 0
      actemp(:,(-nf):nf,:,:) = 0
      astemp(:,(-nf):nf,:,:) = 0
!
!       PERFORM KV (TOROIDAL ANGLE) TRANSFORM
!
      do n = 0, nf
         do kvi = 1, nv
            cosn = cosv(n,kvi)
            sinn = sinv(n,kvi)
            bcos(:,n,1) = bcos(:,n,1) + cosn*source(kvi,:,1)
            bsin(:,n,1) = bsin(:,n,1) + sinn*source(kvi,:,1)
            actemp(:,n,:,1) = actemp(:,n,:,1) + cosn*grpmn(kvi,:,:,1)
            astemp(:,n,:,1) = astemp(:,n,:,1) + sinn*grpmn(kvi,:,:,1)

            if (lasym) then
               bcos(:,n,2) = bcos(:,n,2) + cosn*source(kvi,:,2)
               bsin(:,n,2) = bsin(:,n,2) + sinn*source(kvi,:,2)
               actemp(:,n,:,2) = actemp(:,n,:,2) + cosn*grpmn(kvi,:,:,2)
               astemp(:,n,:,2) = astemp(:,n,:,2) + sinn*grpmn(kvi,:,:,2)
            end if
   
            if (n .ne. 0) then
               bcos(:,(-n),1) = bcos(:,n,1)
               bsin(:,(-n),1) = -bsin(:,n,1)
               actemp(:,(-n),:,1) = actemp(:,n,:,1)
               astemp(:,(-n),:,1) = -astemp(:,n,:,1)

               if (lasym) then
                  bcos(:,(-n),2) = bcos(:,n,2)
                  bsin(:,(-n),2) = -bsin(:,n,2)
               actemp(:,(-n),:,2) = actemp(:,n,:,2)
               astemp(:,(-n),:,2) = -astemp(:,n,:,2)
               end if

            endif
         end do
      end do
!
!       PERFORM KU (POLOIDAL ANGLE) TRANSFORM
!
      do m = 0, mf
         do kui = 1, nu2
            cosm = -cosui(m,kui)
            sinm = sinui(m,kui)
            bvec(m,-nf:nf,1) = bvec(m,-nf:nf,1) + 
     1       bcos(kui,-nf:nf,1)*sinm + bsin(kui,-nf:nf,1)*cosm
            if (lasym) then
            bvec(m,-nf:nf,2) = bvec(m,-nf:nf,2) - 
     1         bcos(kui,-nf:nf,2)*cosm + bsin(kui,-nf:nf,2)*sinm
            end if
         end do
!
!        NOTE: TRANSPOSE KUI, MN INDICES HERE ...
!
         do kui = 1, nu3
            cosm = -cosu(m,kui)*wint(kui*ns*nv)*pi2*pi2
            sinm = sinu(m,kui)*wint(kui*ns*nv)*pi2*pi2
            amatrix(:,m+1:nf*2*mf1+m+1:mf1,1) = amatrix(:,m+1:nf*2*mf1+m
     1         +1:mf1,1) + sinm*transpose(actemp(kui,-nf:nf,:,1)) + 
     2         cosm*transpose(astemp(kui,-nf:nf,:,1))
            if (lasym) then
            amatrix(:,m+1:nf*2*mf1+m+1:mf1,2) = amatrix(:,m+1:nf*2*mf1+m
     1         +1:mf1,2) + sinm*transpose(actemp(kui,-nf:nf,:,2)) + 
     2         cosm*transpose(astemp(kui,-nf:nf,:,2))
            amatrix(:,m+1:nf*2*mf1+m+1:mf1,3) = amatrix(:,m+1:nf*2*mf1+m
     1         +1:mf1,3) - cosm*transpose(actemp(kui,-nf:nf,:,1)) + 
     2         sinm*transpose(astemp(kui,-nf:nf,:,1))
            amatrix(:,m+1:nf*2*mf1+m+1:mf1,4) = amatrix(:,m+1:nf*2*mf1+m
     1         +1:mf1,4) - cosm*transpose(actemp(kui,-nf:nf,:,2)) + 
     2         sinm*transpose(astemp(kui,-nf:nf,:,2))
            end if
         end do
      end do

      deallocate (bcos, bsin, actemp, astemp, source, stat = i)

!
!       ZERO BVEC(0,n) AND AMATRIX(0,n,m`,n`) FOR n < 0
!
      bvec(0,:0,1) = 0
      if (lasym) bvec(0,:0,2) = 0
!
!     M = 0 MODES
!
      amatrix(:nf*mf1+1:mf1,:,1) = 0
!
!     ADD DIAGONAL ELEMENT TO AMATRIX
!
      do mn = 1, mnpd
         amatrix(mn,mn,1) = amatrix(mn,mn,1) + pi3
      end do

      if (lasym) then
         amatrix(:nf*mf1+1:mf1,:,2) = 0
         amatrix(:nf*mf1+1:mf1,:,3) = 0
         amatrix(:nf*mf1+1:mf1,:,4) = 0
         do mn = 1, mnpd
            amatrix(mn,mn,4) = amatrix(mn,mn,4) + pi3
         end do
      end if
 
!
!       PUT ELEMENTS INTO SQUARE MATRIX
!
      amatsq(:mnpd,:mnpd) = amatrix(:,:,1)                      !Sin-Sin
 
      if (lasym) then
         amatsq(1+mnpd:mnpd*2,:mnpd) = amatrix(:,:,2)           !Sin-Cos
         amatsq(:mnpd,1+mnpd:mnpd*2) = amatrix(:,:,3)           !Cos-Sin
         amatsq(1+mnpd:mnpd*2,1+mnpd:mnpd*2) = amatrix(:,:,4)   !Cos-Cos
      end if
 
      end subroutine fouri

      
      subroutine fourp (grpmn, grp, istore, istart, iend)
      use vacmod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: istart, iend, istore
      real(rprec), intent(in) :: grp(nuv,istore)
      real(rprec) :: grpmn(nuv2,0:mf,-nf:nf,*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: n, kv, ku, ip, iu, m
      real(rprec), allocatable, dimension(:,:,:,:) :: g1, g2
      real(rprec), allocatable :: kernels(:), kernelc(:), 
     1     gcos(:,:),gsin(:,:)
      real(rprec) :: cosm, sinm, cosn, sinn
C-----------------------------------------------
!
!     PERFORM KV (TOROIDAL ANGLE) TRANSFORM
!     NOTE: THE m,n INDICES HERE CORRESPOND TO THE FIRST INDEX OF AMATRIX
!     NOTE: THE .5 FACTOR (IN COSN,SINN) ACCOUNTS FOR THE SUM IN KERNELM
!     ON ENTRY THE FIRST TIME, GRPMN IS SIN,COS * Kmn(analytic)
!
      allocate (g1(istore,nu2,0:nf,2), g2(istore,nu2,0:nf,2),
     1   kernels(istore), kernelc(istore), gcos(istore,2), 
     2   gsin(istore,2), stat = m)
      if (m .ne. 0) stop 'Allocation error in fourp'
      
      g1 = 0
      g2 = 0
      
      do 10 n = 0,nf
        do 10 kv = 1,nv
          cosn = 0.5_dp*onp*cosv(n,kv)
          sinn = 0.5_dp*onp*sinv(n,kv)
          iu = kv
          do ku = 1,nu2
            do ip = 1,istore
              kernels(ip) = grp(iu,ip) - grp(imirr(iu),ip)        !sin symmetry
              g1(ip,ku,n,1) = g1(ip,ku,n,1) + cosn*kernels(ip)
              g2(ip,ku,n,1) = g2(ip,ku,n,1) + sinn*kernels(ip)
              if (lasym) then
              kernelc(ip) = grp(iu,ip) + grp(imirr(iu),ip)        !cos symmetry
              g1(ip,ku,n,2) = g1(ip,ku,n,2) + cosn*kernelc(ip)
              g2(ip,ku,n,2) = g2(ip,ku,n,2) + sinn*kernelc(ip)
              end if
            end do
            iu = iu + nv
          end do   
 10   continue

!
!     PERFORM KU (POLOIDAL ANGLE) TRANFORM
!
      do 30 m = 0,mf
        do 30 ku = 1,nu2
          cosm = -cosui(m,ku)
          sinm =  sinui(m,ku)
          do 30 n= 0,nf
            do ip = 1,istore
              gcos(ip,1) = g1(ip,ku,n,1)*sinm
              gsin(ip,1) = g2(ip,ku,n,1)*cosm
              grpmn(ip+istart,m,n,1) = grpmn(ip+istart,m,n,1)
     1        + gcos(ip,1) + gsin(ip,1)
            end do

            if (n .ne. 0) then
              do ip = 1,istore
                grpmn(ip+istart,m,-n,1) = grpmn(ip+istart,m,-n,1)
     1        + gcos(ip,1) - gsin(ip,1)
              end do
            endif  

            if (.not.lasym) cycle

            do ip = 1, istore
              gcos(ip,2) =-g1(ip,ku,n,2)*cosm
              gsin(ip,2) = g2(ip,ku,n,2)*sinm
              grpmn(ip+istart,m,n,2) = grpmn(ip+istart,m,n,2)
     1           + gcos(ip,2) + gsin(ip,2)
            end do

            if (n .ne. 0) then
              do ip = 1,istore
                 grpmn(ip+istart,m,-n,2) = grpmn(ip+istart,m,-n,2)
     1           + gcos(ip,2) - gsin(ip,2)
              end do
            end if  
 30   continue

      istart = iend

      deallocate (g1, g2, kernels, kernelc, gcos, gsin, stat = m)

      end subroutine fourp

      
      subroutine greenf(delgr, delgrp, ip)
      use vacmod
      use vparams, only: one, zero
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ip
      real(rprec), dimension(nuv), intent(out) :: delgr, delgrp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, dimension(2) :: ilow, ihigh
      integer :: ivoff, iskip, iuoff, i, kp, nloop
      real(rprec), dimension(:), allocatable ::
     1    ftemp, gsave, htemp, ga1, ga2, dsave
      real(rprec):: z2p, rcosip, rsinip, cosp, sinp,
     1    sxsave, sysave
C-----------------------------------------------
!
!       ON EXIT, DELGR IS THE DIFFERENCE OF "GREEN'S FUNCTION"
!       AND ANALYTIC APPROXIMATION, SUMMED OVER FIELD PERIODS
!       DELGRP IS DIFFERENCE OF DERIVATIVE OF "GREEN'S FUNCTION"
!       AND ANALYTIC APPROXIMATION
!
!       COMPUTE OFFSETS FOR U,V ANGLE DIFFERENCES AND CONSTANTS
!
      allocate (ftemp(nuv), gsave(nuv), htemp(nuv), ga1(nuv), ga2(nuv),
     1          dsave(nuv), stat=i)
      if (i .ne. 0) stop 'allocation error in greenf'
            
      ilow(1) = 1
      ilow(2) = ip + 1
      ihigh(1) = ip - 1
      ihigh(2) = nuv
      ivoff = nuv + 1 - ip
      iskip = (ip - 1)/nv
      iuoff = nuv - nv*iskip
      z2p = -2*z1b(ip)
      rcosip = -2*rcosuv(ip)
      rsinip = -2*rsinuv(ip)
      delgr(ip)  = zero
      delgrp(ip) = zero
!
!     INITIALIZE ANALYTIC APPROXIMATIONS AND COMPUTE FIELD-PERIOD
!     INVARIANT VECTORS
!
      do i = 1, nuv
         ga1(i) = tanu(i + iuoff)*(guu_b(ip)*tanu(i+iuoff)+guv_b(ip)*
     1      tanv(i+ivoff)) + gvv_b(ip)*tanv(i+ivoff)*tanv(i+ivoff)
         ga2(i) = tanu(i + iuoff)*(auu(ip)*tanu(i+iuoff)+auv(ip)*
     1      tanv(i+ivoff)) + avv(ip)*tanv(i+ivoff)*tanv(i+ivoff)
         gsave(i) = rb2(ip) + rb2(i) + z1b(i)*z2p
         dsave(i) = drv(ip) + z1b(i)*snz(ip)
      end do
!
!     SUM OVER FIELD-PERIODS
!
      do kp = 1, nfper
         cosp = rcosip*cosper(kp) + rsinip*sinper(kp)
         sinp = rsinip*cosper(kp) - rcosip*sinper(kp)
         sxsave = -0.5_dp*(snr(ip)*cosp-snv(ip)*sinp)/r1b(ip)
         sysave = -0.5_dp*(snr(ip)*sinp+snv(ip)*cosp)/r1b(ip)
         if (kp .le. 1) then
            do nloop = 1, 2
               do i = ilow(nloop), ihigh(nloop)
                 ga2(i) = ga2(i)/ga1(i)
                 ga1(i) = one/sqrt(ga1(i))
                 ftemp(i) = one/(gsave(i) + cosp*rcosuv(i) 
     1                    + sinp*rsinuv(i))
                 htemp(i) = sqrt(ftemp(i))
                 delgrp(i) = -ga2(i)*ga1(i) + ftemp(i)*htemp(i)*
     1              (rcosuv(i)*sxsave + rsinuv(i)*sysave + dsave(i))
                 delgr(i) = htemp(i) - ga1(i)
               end do
            end do
         else
            do i = 1,nuv        
              ftemp(i) =one/(gsave(i) + cosp*rcosuv(i) + sinp*rsinuv(i))
              htemp(i) = sqrt(ftemp(i))
              delgrp(i) = delgrp(i) + ftemp(i)*htemp(i)*
     1           (rcosuv(i)*sxsave + rsinuv(i)*sysave + dsave(i))
              delgr(i) = delgr(i) + htemp(i)
           end do  
         endif
      end do

      deallocate (ftemp, gsave, htemp, ga1, ga2, dsave, stat=i)

      end subroutine greenf

      
      subroutine scalpot(bvec, amatrix, wint, ns, istore_max, ivacskip)
      use vacmod
      use vparams, only: zero
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: ns, ivacskip, istore_max
      real(rprec), intent(out) :: bvec(mnpd2), amatrix(mnpd2*mnpd2)
      real(rprec), intent(in) :: wint(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ip, istore, istart
      real(rprec), allocatable :: grpmn(:), green(:), gstore(:)
      real(rprec), allocatable :: greenp(:,:)
C-----------------------------------------------
      if (.not.allocated(amatsav)) 
     1   stop 'AMATSAV: Allocation error in scalpot' 

      allocate (grpmn(nuv2*mnpd2), stat=ip)
      if (ip .ne. 0) stop 'GRPMN: Allocation error in scalpot'
!
!     INITIALIZE VECTORS
!
      if (ivacskip .eq. 0) then
         amatrix = zero
         grpmn   = zero
      endif

      bvec = zero
!
!       COMPUTE TRANFORM OF ANALYTIC SOURCE AND KERNEL
!       ON EXIT, BVEC CONTAINS THE TRANSFORM OF THE ANALYTIC SOURCE
!       AND GRPMN CONTAINS SIN,COS * TRANSFORM OF NORMAL DERIVATIVE
!       OF THE "GREEN'S FUNCTION"
!
      call analyt (grpmn, bvec, ivacskip)

      if (ivacskip .ne. 0) then
         bvec = bvec + bvecsav
      else
         allocate (green(nuv), gstore(nuv), greenp(nuv,istore_max))
         bvecsav = bvec
         gstore  = zero
!
!       COMPUTE SURFACE INTEGRALS OF SOURCE, "GREEN'S FUNCTION" NEEDED
!       FOR SPECTRAL DECOMPOSITION OF POTENTIAL INTEGRAL EQUATION
!       NOTE: SOURCE IS THE RHS, KERNEL IS THE LHS OF EQN.
!
         istart = 0
         do ip = 1, nuv2
            istore = 1 + mod(ip-1,istore_max)
!
!       COMPUTE EXACT AND APPROXIMATE "GREEN'S" FUNCTION AND GRADIENT
!
            call greenf (green, greenp(1,istore), ip)
            gstore = gstore + bexni(ip)*green
!
!       COMPUTE FOURIER INTEGRAL OF GRADIENT KERNEL ON UNPRIMED MESH
!
            if (istore.eq.istore_max .or. ip.eq.nuv2) 
     1        call fourp (grpmn, greenp, istore, istart, ip)

         end do
 
!
!       COMPUTE FOURIER INTEGRAL OF GRADIENT (SOURCE) PRIMED (UNPRIMED)
!
         call fouri (grpmn, gstore, amatrix, amatsav, bvec, wint, ns)
         deallocate (green, greenp, gstore)
      endif

      deallocate (grpmn)
      amatrix = amatsav

      if (ivacskip .ne. 0) return 
!
!     SAVE NON-SINGULAR CONTRIBUTION TO BVEC IN BVECSAV
!
      bvecsav(:mnpd2) = bvec - bvecsav(:mnpd2)

      end subroutine scalpot

      
      subroutine surface(rc, rs, zs, zc, xm, xn, mnmax)
      use vacmod
      use vparams, only: zero
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer mnmax
      real(rprec), dimension(mnmax) :: rc, rs, zs, zc, xm, xn
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
      integer :: i, mn, m, n, n1
      real(rprec), allocatable, dimension(:) :: 
     1   ruu, ruv, rvv, zuu, zuv, zvv, cosmn1, sinmn1
C-----------------------------------------------
!
!       THIS ROUTINE COMPUTES THE SURFACE VALUES OF R,Z AND DERIVATIVES
!
!
!       Compute R & Z (and their derivatives) on surface
!
      allocate (ruu(nuv2), ruv(nuv2), rvv(nuv2), zuu(nuv2), zuv(nuv2),
     1          zvv(nuv2), cosmn1(nuv2), sinmn1(nuv2), stat = i)
      if (i .ne. 0) stop 'Allocation error in SURFACE'
       
      r1b = 0;   rub = 0;   rvb = 0;  ruu = 0; ruv = 0; rvv = 0
      z1b = 0;   zub = 0;   zvb = 0;  zuu = 0; zuv = 0; zvv = 0
      do mn = 1, mnmax
         m = nint(xm(mn))
         n = nint(xn(mn)/(nfper))
         n1 = abs(n)
         cosmn1(:) = cosu1(:,m)*cosv1(:,n1) + csign(n)*sinu1(:,m)*
     1               sinv1(:,n1)
         sinmn1(:) = sinu1(:,m)*cosv1(:,n1) - csign(n)*cosu1(:,m)*
     1               sinv1(:,n1)
         do i = 1, nuv2
            r1b(i) = r1b(i) + rc(mn) * cosmn1(i)
            rub(i) = rub(i) - xm(mn) * rc(mn) * sinmn1(i)
            rvb(i) = rvb(i) + xn(mn) * rc(mn) * sinmn1(i)
            z1b(i) = z1b(i) + zs(mn) * sinmn1(i)
            zub(i) = zub(i) + xm(mn) * zs(mn) * cosmn1(i)
            zvb(i) = zvb(i) - xn(mn) * zs(mn) * cosmn1(i)
            ruu(i) = ruu(i) - xm(mn)*xm(mn)*rc(mn) * cosmn1(i)
            ruv(i) = ruv(i) + xm(mn)*xn(mn)*rc(mn) * cosmn1(i)
            rvv(i) = rvv(i) - xn(mn)*xn(mn)*rc(mn) * cosmn1(i)
            zuu(i) = zuu(i) - xm(mn)*xm(mn)*zs(mn) * sinmn1(i)
            zuv(i) = zuv(i) + xm(mn)*xn(mn)*zs(mn) * sinmn1(i)
            zvv(i) = zvv(i) - xn(mn)*xn(mn)*zs(mn) * sinmn1(i)
         end do
 
         if (.not.lasym) cycle

         do i = 1, nuv2
            r1b(i) = r1b(i) + rs(mn) * sinmn1(i)
            rub(i) = rub(i) + xm(mn) * rs(mn) * cosmn1(i)
            rvb(i) = rvb(i) - xn(mn) * rs(mn) * cosmn1(i)
            z1b(i) = z1b(i) + zc(mn) * cosmn1(i)
            zub(i) = zub(i) - xm(mn) * zc(mn) * sinmn1(i)
            zvb(i) = zvb(i) + xn(mn) * zc(mn) * sinmn1(i)
            ruu(i) = ruu(i) - xm(mn)*xm(mn)*rs(mn) * sinmn1(i)
            ruv(i) = ruv(i) + xm(mn)*xn(mn)*rs(mn) * sinmn1(i)
            rvv(i) = rvv(i) - xn(mn)*xn(mn)*rs(mn) * sinmn1(i)
            zuu(i) = zuu(i) - xm(mn)*xm(mn)*zc(mn) * cosmn1(i)
            zuv(i) = zuv(i) + xm(mn)*xn(mn)*zc(mn) * cosmn1(i)
            zvv(i) = zvv(i) - xn(mn)*xn(mn)*zc(mn) * cosmn1(i)
         end do
      end do

!
!       COMPUTE METRIC COEFFICIENTS AND AREA ELEMENTS
!       NOTE: guv = .5*np GUV; gvv = np*np* GVV, where GUV, GVV are the
!            real metric elements
!
      do i = 1,nuv2
        guu_b(i) = rub(i)*rub(i) + zub(i)*zub(i)
        guv_b(i) = (rub(i)*rvb(i)+ zub(i)*zvb(i))*onp*2.0_dp
        gvv_b(i) = (rvb(i)*rvb(i)+ zvb(i)*zvb(i)+(r1b(i)*r1b(i)))*onp2
        rb2(i) = (r1b(i)*r1b(i)) + z1b(i)*z1b(i)
        snr(i) = -r1b(i)*zub(i)
        snv(i) = rvb(i)*zub(i) - rub(i)*zvb(i)
        snz(i) = r1b(i)*rub(i)
        drv(i) = (-r1b(i)*snr(i)) - z1b(i)*snz(i)
        auu(i) = (0.5_dp*r1b(i))*(zuu(i)*rub(i) - ruu(i)*zub(i))
        auv(i) = (snv(i)*rub(i)+(zuv(i)*rub(i) 
     1         - ruv(i)*zub(i))*r1b(i))*onp
        avv(i) = (snv(i)*rvb(i)+(0.5_dp*r1b(i))*(zvv(i)*rub(i)
     1         - zub(i)*rvv(i) - snr(i)))*onp2
      end do

      if (.not.lasym) then
         do i = 1 + nv, nuv2 - nv
            rb2(imirr(i)) = rb2(i)
            r1b(imirr(i)) = r1b(i)
            z1b(imirr(i)) =-z1b(i)
         end do
      end if

      do i = 1,nuv
        rcosuv(i) = r1b(i)*cosuv(i)
        rsinuv(i) = r1b(i)*sinuv(i)
      end do  

      deallocate (ruu, ruv, rvv, zuu, zuv, zvv, cosmn1, sinmn1, stat=i)

      end subroutine surface

      
      subroutine tolicu(rcurr, zcurr)
      use vacmod
      use vacwires
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(nv) :: rcurr, zcurr
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, kper, kv
C-----------------------------------------------
 
      if (nv*nfper .ne. nvp) stop ' nvp!=nv*nfper in tolicu'
      i = 0
      do kper = 1, nfper
         do kv = 1, nv
            i = i + 1
            xw(i) = rcurr(kv)*(cosper(kper)*cosuv(kv) - sinper(kper)*
     1         sinuv(kv))
            yw(i) = rcurr(kv)*(sinper(kper)*cosuv(kv) + cosper(kper)*
     1         sinuv(kv))
            zw(i) = zcurr(kv)
         end do
      end do
      xw(nvp+1) = xw(1)
      yw(nvp+1) = yw(1)
      zw(nvp+1) = zw(1)
      do i = 1,nvp
        dx(i) = xw(i+1) - xw(i)
        dy(i) = yw(i+1) - yw(i)
        dz(i) = zw(i+1) - zw(i)
        vx(i) = yw(i)*dz(i) - zw(i)*dy(i)
        vy(i) = zw(i)*dx(i) - xw(i)*dz(i)
        vz(i) = xw(i)*dy(i) - yw(i)*dx(i)
      end do  

      end subroutine tolicu

      
      subroutine precal
      use vparams, only: zero, one, epstan
      use vacmod
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: kp, ku, kuminus, kv, kvminus, i, m, n, mn, 
     1   jmn, kmn, l, istat1
      real(rprec), dimension(0:mf + nf,0:mf,0:nf) :: cmn
      real(rprec) :: argu, argv, dn1, smn, f1, f2, f3
C-----------------------------------------------
!
!     THIS ROUTINE COMPUTES INITIAL CONSTANTS AND ARRAYS
!
      pi2 = 8*atan(one)
      pi3 = 0.5_dp*pi2**3
      pi4 = 2*pi2
      alp = pi2/(nfper)
      alu = pi2/(nu)
      alv = pi2/(nv)
      onp = one/(nfper)
      onp2 = onp*onp
      alvp = onp*alv

!
!     ALLOCATE PERSISTENT ARRAYS. DEALLOCATED IN FILEOUT ROUTINE
!
      allocate( tanu(2*nuv), tanv(2*nuv), sin2v(2*nuv),
     1     sinper(nfper), cosper(nfper), sinuv(nuv), cosuv(nuv),
     2     sinu(0:mf,nu), cosu(0:mf,nu), sinv(-nf:nf,nv),
     3     cosv(-nf:nf,nv), sinui(0:mf,nu), cosui(0:mf,nu),
     4     cmns(0:(mf+nf),0:mf,0:nf), csign(-nf:nf),
     5     sinu1(nuv2,0:mf), cosu1(nuv2,0:mf),
     6     sinv1(nuv2,0:nf), cosv1(nuv2,0:nf), imirr(nuv), 
     7     xmpot(mnpd), xnpot(mnpd), stat=istat1)
      if (istat1.ne.0) stop 'allocation error in precal'


!
!       IMIRR(I) GIVES THE INDEX OF THE POINT TWOPI-THETA(I),TWOPI-ZETA(I)
!
      do kp = 1, nfper
         cosper(kp) = cos(alp*(kp - 1))
         sinper(kp) = sin(alp*(kp - 1))
      end do
      do ku = 1, nu
         kuminus = mod(nu + 1 - ku,nu) + 1
         do kv = 1, nv
            kvminus = mod(nv + 1 - kv,nv) + 1
            i = kv + nv*(ku - 1)
            imirr(i) = kvminus + nv*(kuminus - 1)
            cosuv(i) = cos(alvp*(kv - 1))
            sinuv(i) = sin(alvp*(kv - 1))
         end do
      end do
!
!       NOTE: ACTUAL ANGLE DIFFERENCE IS (KUP-1) - (KU-1)
!
      i = 0
      do ku = 1, 2*nu
         argu = 0.5_dp*alu*(ku - 1)
         do kv = 1, nv
            i = i + 1
            argv = 0.5_dp*alv*(kv - 1)
            sin2v(i) = sin(argv)*sin(argv)
            if (abs(argu - 0.25_dp*pi2)<epstan .or.
     1      abs(argu - 0.75_dp*pi2) < epstan) then
               tanu(i) = 0.9e30_dp
            else
               tanu(i) = 2.0_dp*tan(argu)
            endif
            if (abs(argv - 0.25_dp*pi2) < epstan) then
               tanv(i) = 0.9e30_dp
            else
               tanv(i) = 2.0*tan(argv)
            endif
         end do
      end do
      do m = 0, mf
         l40: do ku = 1, nu
            cosu(m,ku) = cos(alu*(m*(ku - 1)))
            sinu(m,ku) = sin(alu*(m*(ku - 1)))
            cosui(m,ku) = cosu(m,ku)*alu*alv*2.0
            sinui(m,ku) = sinu(m,ku)*alu*alv*2.0
            if (ku.eq.1 .or. ku.eq.nu2) cosui(m,ku) = 0.5_dp*cosui(m,ku)
            do kv = 1, nv
               i = kv + nv*(ku - 1)
               if (i > nuv2) cycle  l40
               cosu1(i,m) = cosu(m,ku)
               sinu1(i,m) = sinu(m,ku)
            end do
         end do l40
      end do
      do n = -nf, nf
         dn1 = alvp*(n*nfper)
         csign(n) = sign(one,dn1)
         l50: do ku = 1, nu
            do kv = 1, nv
               i = kv + nv*(ku - 1)
               cosv(n,kv) = cos(dn1*(kv - 1))
               sinv(n,kv) = sin(dn1*(kv - 1))
               if (i.gt.nuv2 .or. n.lt.0) cycle  l50
               cosv1(i,n) = cosv(n,kv)
               sinv1(i,n) = sinv(n,kv)
            end do
         end do l50
      end do
      mn = 0
      do n = -nf, nf
         do m = 0, mf
            mn = mn + 1
            xmpot(mn) = (m)
            xnpot(mn) = (n*nfper)
         end do
      end do
!
!       COMPUTE "CMN'S" AND THEIR SUMS , EQ (A14 AND A13) IN J.COMP.PHYS PAPER
!
      do m = 0, mf
         do n = 0, nf
            jmn = m + n
            kmn = abs(m - n)
            smn = 0.5_dp*(jmn + kmn)
            f1 = 1.0_dp
            f2 = 1.0_dp
            f3 = 1.0_dp
            do i = 1, kmn
               f1 = f1*(smn + (1 - i))
               f2 = f2*(i)
            end do
            cmn(0:mf+nf,m,n) = zero
            do l = kmn, jmn, 2
               cmn(l,m,n) = f1/(f2*f3)*((-1)**((l - m + n)/2))
               f1 = f1*0.25_dp*((jmn + l + 2)*(jmn - l))
               f2 = f2*0.5_dp*(l + 2 + kmn)
               f3 = f3*0.5_dp*(l + 2 - kmn)
            end do
         end do
      end do
!
!       The ALP factor comes from integral over field periods
!
        do m = 1,mf
           do n = 1,nf
              cmns(0:mf+nf,m,n) = 0.5_dp*alp*(cmn(0:mf+nf,m,n) +
     1        cmn(0:mf+nf,m-1,n) + cmn(0:mf+nf,m,n-1) + 
     2        cmn(0:mf+nf,m-1,n-1))
           end do
        end do   
      cmns(0:mf+nf,1:mf,0) = (0.5_dp*alp)*(cmn(0:mf+nf,1:mf,0)
     1                       + cmn(0:mf+nf,:mf-1,0))
      cmns(0:mf+nf,0,1:nf) = (0.5_dp*alp)*(cmn(0:mf+nf,0,1:nf)
     1                    + cmn(0:mf+nf,0,:nf-1))
      cmns(0:mf+nf,0,0)    = (0.5_dp*alp)*(cmn(0:mf+nf,0,0)
     1                    + cmn(0:mf+nf,0,0)) 

      end subroutine precal
EOF

cat > reconstruct.f << "EOF"
      subroutine axisopt(fsq, r00, iresidue, ivac)
      use vsvd
      use vparams, only: zero, one, nthreed
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: iresidue, ivac
      real(rprec) :: fsq, r00
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: smax = 0.998_dp
      real(rprec), parameter :: smin = 0.985_dp
      character*(60), parameter :: optbegin =
     1   'Begin variation of Raxis to minimize total RMS error'
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec) :: delstep, dedrmax, factor, delerr,
     1   dedr, rstepx1
      real(rprec), save :: delstep_old, errmax, errold, rstepx,
     1   raxold, scale
C-----------------------------------------------
 
 
      if (iresidue.lt.1 .or. errsvd*1.e6_dp.lt.one .or.
     1    fsq.gt.fturnon_axis .or. ivac .le.2) return 
!
!     MOVE R-AXIS BASED ON dR/dt = (-dEsvd/dR)
!     LIMIT MAXIMUM RSTEPX TO RSTEPX0
!     TRY TO FIND ZERO-CROSSING IN dEsvd/dR (ESTIMATED NUMERICALLY)
!
 
      if (iresidue .eq. 1) then                    !First time through
         iresidue = 2
         raxold = r00
         errold = errsvd
         errmax = zero
         rstepx = rstepx0
         scale = smax
         if (iopt_raxis .gt. 0) then
            write (*, 115) optbegin
            write (nthreed, 115) optbegin
         endif
      else
         delerr = errsvd - errold                !delta E-svd
         delstep = r00 - raxold                  !delta R-axis
         if (delerr.ne.zero .and. abs(delstep).gt.1.e-3_dp*rstepx0) then
            dedr = delerr/delstep
            errmax = max(errmax,errsvd)
            dedrmax = 2.0*errmax/rwidth
            rstepx1 = min(one,abs(dedr)/dedrmax)*rsfac*rstepx0
            factor = sign(one,(-dedr))        !Move in -dE/dR direction
            rstepx = rstepx1*factor
            scale = smax
            if (delstep*delstep_old .le. zero) scale = smin
            delstep_old = delstep
            raxold = r00
            errold = errsvd
         endif
      endif
      rsfac = scale*rsfac
c-5/1/96 raxmse = raxmse + rstepx
      raxmse = raxold + rstepx
  115 format(2x,a)
 
      end subroutine axisopt


      subroutine chisq(amat_i, amat_p, data, idata, isize, itotal)
      use vmec_main
      use vsvd
      use vspline
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: itotal
      integer, dimension(*) :: idata, isize
      real(rprec), dimension(isnodes,*) :: amat_i
      real(rprec), dimension(ipnodes,*) :: amat_p
      real(rprec), dimension(*) :: data
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, n, i1, ispec, is(5), ip(5), j2
      real(rprec) :: delsq, delp, dels
      character*(130) label
C-----------------------------------------------
!
!       COMPUTES CHI**2 FROM DIFFERENT SUBROUTINES AT VARIOUS TIME-STEPS
!       WRITTEN BY D.K. LEE (3/93)
!
 
         chisqerr(:jchix) = zero
 
         do i = 1, itotal
            delsq = (sum(ystark(:isnodes)*amat_i(:isnodes,i)) +
     1               sum(ythom(:ipnodes)*amat_p(:ipnodes,i))-data(i))**2
            if (i.ge.idata(ithom0) .and. i<idata(ithom0)+isize(ithom0)) 
     1         then
               chisqerr(ithom0) = chisqerr(ithom0) + delsq
            else if (i.ge.idata(istark0) .and. i<idata(istark0)
     1            +isize(istark0)) then
               i1 = i - idata(istark0) + 1
               if (i1 .eq. islope) then
                  chisqerr(islope0) = delsq
               else if (i1 .eq. icurrout) then
                  chisqerr(icurr0) = delsq
               else
                  chisqerr(istark0) = chisqerr(istark0) + delsq
               endif
            else if (i.ge.idata(idiam0) .and. i<idata(idiam0)
     1            +isize(idiam0)) then
               chisqerr(idiam0) = chisqerr(idiam0) + delsq
            else if (i.ge.idata(iflxs0) .and. i<idata(iflxs0)
     1            +isize(iflxs0)) then
               chisqerr(iflxs0) = chisqerr(iflxs0) + delsq
            else if (i.ge.idata(ibrzfld) .and. i<idata(ibrzfld)
     1            +isize(ibrzfld)) then
               chisqerr(ibrzfld) = chisqerr(ibrzfld) + delsq
            endif
         end do
 
!
         errsvd = sum(chisqerr(:jchix))
         if (.not.lpprof) errsvd = errsvd - chisqerr(ithom0)
 
      if (iequi.ne.1 .or. .not.lrecon) return
         if (.not.lpprof) then
            write (nthreed, 15)
         else
            write (nthreed, 10)
         endif
         if (lpprof) then
            do n = 1, nchistp
               write (nthreed, 20) nchi2(n), chi2(ithom0,n), 
     1            chi2(istark0,n), chi2(icurr0,n), chi2(idiam0,n), 
     2            chi2(iflxs0,n), chi2(ibrzfld,n), chi2(jchix1,n)
            end do
         else
            do n = 1, nchistp
               write (nthreed, 20) nchi2(n), chi2(istark0,n), 
     1            chi2(icurr0,n), chi2(idiam0,n), chi2(iflxs0,n), 
     2            chi2(ibrzfld,n), chi2(jchix1,n)
            end do
         endif
 
!
!       PRINT OUT MATRIX ELEMENTS (5 EACH FOR PRESSURE, IOTA)
!
         write (nthreed, 200)
         delp = (ipnodes - 1)/4.
         dels = (isnodes - 1)/4.
         ip(1) = 1
         is(1) = 1
         ip(2:4) = ip(1) + int((((/(j2,j2=2,4)/)) - 1)*delp)
         is(2:4) = is(1) + int((((/(j2,j2=2,4)/)) - 1)*dels)
         ip(5) = ipnodes
         is(5) = isnodes
         write (label, 210) ip(1), ip(2), ip(3), ip(4), ip(5), is(1), 
     1      is(2), is(3), is(4), is(5)
         write (nthreed, 220) label
         ispec = 0
         do i = 1, itotal
            if (i.ge.idata(ithom0) .and. 
     1         i.lt.idata(ithom0)+isize(ithom0)) then
               i1 = i - idata(ithom0) + 1
               call printmatrix (amat_p(1,i), amat_i(1,i), data(i), i, 
     1            i1, ip, is, '   PRES  (')
            else if (i.ge.idata(istark0) .and. i<idata(istark0)
     1            +isize(istark0)) then
               i1 = i - idata(istark0) + 1
               if (i1 .eq. islope) then
                  ispec = ispec - 1
                  call printmatrix (amat_p(1,i), amat_i(1,i), data(i), 
     1               i, 1, ip, is, '  IOTA0  (')
               else if (i1 .eq. icurrout) then
                  ispec = ispec - 1
                  call printmatrix (amat_p(1,i), amat_i(1,i), data(i), 
     1               i, 1, ip, is, ' CURRENT (')
               else
                  call printmatrix (amat_p(1,i), amat_i(1,i), data(i), 
     1               i, i1 + ispec, ip, is, '   MSE   (')
               endif
            else if (i.ge.idata(idiam0) .and. i<idata(idiam0)
     1            +isize(idiam0)) then
               i1 = i - idata(idiam0) + 1
               call printmatrix (amat_p(1,i), amat_i(1,i), data(i), i, 
     1            i1, ip, is, ' DIAMAG  (')
            else if (i.ge.idata(iflxs0) .and. i<idata(iflxs0)
     1            +isize(iflxs0)) then
               i1 = i - idata(iflxs0) + 1
               call printmatrix (amat_p(1,i), amat_i(1,i), data(i), i, 
     1            i1, ip, is, ' FLUXES  (')
            else if (i.ge.idata(ibrzfld) .and. i<idata(ibrzfld)
     1            +isize(ibrzfld)) then
               i1 = i - idata(ibrzfld) + 1
               call printmatrix (amat_p(1,i), amat_i(1,i), data(i), i, 
     1            i1, ip, is, '  BR-BZ  (')
            endif
         end do
   20 format(i6,1p8e12.4)
   10 format(/,30x,'ABSOLUTE CHI-SQ ERROR BY DATA TYPE'/,30x,
     1   '(NOT NORMED BY NUMBER DATA POINTS)'/,20x,
     2   'NOTE: STARK CHISQ MAY BE EVALUATED AT REDISTRIBUTED KNOTS'/,
     3   '  ITER   Thomscat      Stark',5x,'Current',5x,'Diamag.',5x,
     4   'Saddle',6x,' B-Loops',6x,'TOTAL',/,1x,5('-'),7(2x,10('-')))
   15 format(/,30x,'ABSOLUTE CHI-SQ ERROR BY DATA TYPE'/,30x,
     1   '(NOT NORMED BY NUMBER DATA POINTS)'/,20x,
     2   'NOTE: STARK CHISQ MAY BE EVALUATED AT REDISTRIBUTED KNOTS'/,
     3   '  ITER      Stark',5x,'Current',5x,'Diamag.',5x,'Saddle',6x,
     4   ' B-Loops',6x,'TOTAL',/,1x,5('-'),7(2x,10('-')))
  200 format(//,38x,'SPLINE MATRIX ELEMENTS BY DATA TYPE'/,30x,
     1   ' AI(i,j)*iota(j) + AP(i,j)*[mu0*pres(j)] = DATA(i)'/,30x,
     2   ' NOTE: DATA(I) IS THE RAW DATA NORMED TO SIGMA(I)'//)
  210 format('   I   TYPE         DATA(I)','  AP(I,',i2,')  AP(I,',i2,
     1   ')  AP(I,',i2,')  AP(I,',i2,')','  AP(I,',i2,')  AI(I,',i2,
     2   ')  AI(I,',i2,')  AI(I,',i2,')','  AI(I,',i2,')  AI(I,',i2,')')
  220 format(a,/,3x,'-',3x,4('-'),9x,7('-'),10(2x,8('-')))
 
      end subroutine chisq

      
      subroutine printmatrix(amatp, amati, data, i, i1, ip, is, type)
      use vparams, only: rprec, nthreed
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer i, i1
      real(rprec) data
      character*(10) type
      integer, dimension(*) :: ip, is
      real(rprec), dimension(*) :: amatp, amati
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: k
C-----------------------------------------------
 
      write (nthreed, 10) i, type, i1, data, (amatp(ip(k)),k=1,5), (
     1   amati(is(k)),k=1,5)
 
   10 format(1x,i3,a10,i2,')',1p11e10.2)
 
      end subroutine printmatrix

       
      subroutine store_chisq
      use vmec_main
      use vsvd
      use vspline
      implicit none
C-----------------------------------------------
!
!       COMPUTES CHI**2 FROM DIFFERENT SUBROUTINES AT VARIOUS TIME-STEPS
!       WRITTEN BY D.K. LEE (3/93)
!
 
      if (mod(iter2,nstep).ne.10 .and. iequi.eq.0) return
         chisqerr(jchix1) = sum(chisqerr(:jchix))
         if (.not.lpprof) chisqerr(jchix1) = chisqerr(jchix1) 
     1      - chisqerr(ithom0)
         nchistp = nchistp + 1
         if (nchistp .gt. mstp) return
         chi2(:,nchistp) = chisqerr
         nchi2(nchistp) = iter2 - 10
         if (iequi .eq. 1) nchi2(nchistp) = iter2
         if (iter2 .eq. 10) nchi2(nchistp) = 1

      end subroutine store_chisq

      
      subroutine findphi(reven, rodd, rmeas, dse, dso, rmid, ismeas, 
     1   iumeas, indexr, npts)
      use vmec_main
      use realspace
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer npts
      integer, dimension(npts) :: ismeas, iumeas
      integer, dimension(2*ns) :: indexr
      real(rprec), dimension(ns,nzeta,*) :: reven, rodd
      real(rprec), dimension(npts) :: rmeas, dse, dso
      real(rprec), dimension(2*ns) :: rmid
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, k, itemp, jtemp, i1, j1, ns2
C-----------------------------------------------
 
 
!
!       THIS ROUTINE FINDS THE LOWER FLUX INDEX [=INDEXS] CORRESPONDING
!       TO THE MEASURED VALUE [=RMEAS] OF (R,Z=0) ALONG THE MIDPLANE
!       (THETA=0 OR PI).
!       THE QUANTITIES RMEAS(K=1,NPTS) ARE INTERPOLATED AS FOLLOWS:
!
!       RMEAS(K) = R[ISMEAS(K),IUMEAS(K)]*[1-DSO(K)] +
!                  R[ISMEAS(K)+1,IUMEAS(K)]*DSO(K)
!
!       BECAUSE OF THE SQRT(S) BEHAVIOUR OF R IN THE FIRST RADIAL ZONE,
!       THE ACTUAL S-INTERPOLAND IN THE FIRST ZONE IS DSO(K)**2 = DSE(K).
!       IN ALL OTHER ZONES, DSE(K) = DSO(K).
!
 
      if (npts .le. 0) return 
      ns2 = 2*ns
!
!     COMPUTE THE GRID VALUES (S-COORDINATE) OF R ALONG THE MIDPLANE,
!     STARTING AT THETA=PI (I=NTHETHA2) AND ENDING AT THETA=0 (I=1)
!
 
      rmid(:ns) = reven(indexr(:ns),1,ntheta2) + sqrts(indexr(:ns))
     1   *rodd(indexr(:ns),1,ntheta2)
      rmid(ns+1:ns2) = reven(indexr(ns+1:ns2),1,1) + 
     1   sqrts(indexr(ns+1:ns2))*rodd(indexr(ns+1:ns2),1,1)
 
!
!     FIND THE RADIAL ZONE INDEX [=ITEMP], WHICH BRACKETS THE MEASURED R-VALUE
!
!     RMID(ITEMP-1) .le. RMEAS .le. RMID(ITEMP)
!
 
      do k = 1, npts
         itemp = 0
         do i = 1, ns2 - 1
            if (rmeas(k) .lt. rmid(i)) then
               itemp = i
               go to 100
            endif
         end do
         itemp = ns2
!
!         FIND FLUX-COORDINATE S-INDEX [=ISMEAS], POLOIDAL ANGLE
!         INDEX [=IUMEAS], AND INTERPOLAND [=DSO]
!
 
  100    continue
         if (itemp.gt.1 .and. itemp.lt.ns2) then
            i1 = itemp - 1
            jtemp = indexr(itemp)
            j1 = indexr(i1)
            dso(k) = (rmeas(k)-rmid(i1))/(rmid(itemp)-rmid(i1))
            if (j1 .lt. jtemp) then                 !THETA = 0
               ismeas(k) = j1
               iumeas(k) = 1
            else
               ismeas(k) = jtemp
               dso(k) = 1.0 - dso(k)
               iumeas(k) = ntheta2
            endif
         else
            dso(k) = 1.0
            ismeas(k) = indexr(itemp) - 1
!           IGNORE MEASURED POINTS OUTSIDE GRID
            if (itemp.eq.1 .or. rmeas(k).gt.rmid(ns2-1)) dso(k) = -1.0
            if (itemp .eq. ns2) iumeas(k) = 1
         endif
!
!        ACCOUNT FOR SQRT(S) SINGULARITY IN 1st ZONE
!        DSE IS THE S-FRACTIONAL INTERPOLAND
!
         dse(k) = dso(k)
         if (ismeas(k) .eq. 1) dse(k) = dso(k)*dso(k)
      end do
 
      end subroutine findphi

      
      subroutine fixrecon(ier)
      use vmec_main
      use vsvd
      use vspline
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ier
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: angle_variance  = 0.3_dp
      real(rprec), parameter :: radial_variance = 0.3_dp
      real(rprec), parameter :: p_threshold = 1.e-3_dp
      real(rprec), parameter :: c1p5 = 1.5_dp
      integer, parameter :: inode_max = 15
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, dimension(itse) :: isortp
      integer, dimension(1) :: isamax
      integer :: istat1 = 0, istat2 = 0, istat3 = 0, ii, i, ineg
     1   , ipos, ipmax, ioff, ileft, n, n1, icount, index, ind, m, m1
      real(rprec) :: delstark, datapos, dataneg
      real(rprec), dimension(imse+1) :: datalsq_s, stark_temp
      real(rprec) :: datamax, datamin, t1, rneg,
     1   rpos, datanorm, presmax, presin, presout, presmin, tsign
      real(rprec), dimension(itse) ::
     1   datalsq_p, ythom0, y2thom0, qtemp
      logical :: l1v(imse), l4v(itse)
C-----------------------------------------------
!
!
!                INDEX OF LOCAL VARIABLES
!
!         needflx    =NEEDIT, loop required for flux match
!                    >=ISYMCOIL, loop required but flux computed by
!                               invoking symmetry in Z
!                    =IDONTNEED, loop not required for flux match
!        needbfld    =NEEDIT, loop required for B-field match
!                    =ISAMECOIL, loop at same position as previous loop
!                    >=ISYMCOIL, loop required but B-field computed
!                               by invoking symmetry in Z
!                    =IDONTNEED, loop not required in B-field match
!          dsiext    connected flux loop signals due to external coils
!          plflux    array of measured (inferred) plasma contrib. to flux loops
!          plbfld    array of measured (inferred) plasma contrib. to B-loops
!

      call free_mem_recon

      ier = 0
!
!       ONLY QUANTITIES INDEPENDENT OF RADIAL MESH CAN GO HERE
!
!       STARK-DATA CONSISTENCY CHECK
!
      delstark = tan(dcon*angle_variance)

!
!       SORT STARK DATA IN ASCENDING ORDER IN R-SPACE
!       AND RE-INDEX DATASTARK, SIGMA_STARK ARRAYS ON ASCENDING RSTARK ARRAY
!
!       SCALE MOTIONAL STARK DATA TO RADIANS
!       RSTARK = R position along Z=0 plane of measurement
!       DATASTARK = ARCTAN ( Bpol/Btor ) at RSTARK, in degrees
!       AND IS CONVERTED TO Bpol/Btor
!

      ii = 0
      l1v(:imse) = (rstark(:imse) .gt. zero)
      do i = 1, imse
         if (l1v(i)) ii = ii + 1
      end do
      if (ii .ne. imse) then
         print *, 'There is a zero in the RSTARK array ?!'
         ier = 1
         return
      endif

      allocate (indexs1(imse+3), indexu1(imse+3), isortr(imse+3),
     1   isorts(imse+3), delse1(imse+3), delso1(imse+3),
     2   starkcal(imse+3), qmeas(imse+3), qcalc(imse+3),
     3   fpsical(imse+3), stark_weight(imse+3),
     4   rsort(imse+3), rsort0(imse+3), stat=istat1)
      if (istat1 .ne. 0) then
         print *, ' ISTAT1 = ', istat1, ' ALLOCATING INDEXS1'
         stop
      endif

      datalsq_s(:imse) = datastark(:imse)
      stark_temp(:imse) = sigma_stark(:imse)
      call sort_data (rstark, isorts, imse)
      do i = 1, imse
         datastark(i) = datalsq_s(isorts(i))
         sigma_stark(i) = stark_temp(isorts(i))
         if (sigma_stark(i) .ge. cbig) then
            print *, 'SIGMA_STARK missing'
            ier = 1
            return
         endif
         if (sigma_stark(i) .lt. zero) sigma_stark(i) =
     1     abs(sigma_stark(i)*datastark(i))
          !CONVERT TO Bpol/Btor ... applying profile offsets along the way!
         datastark(i) = tan(dcon*(datastark(i)+mseangle_offset+
     1      mseangle_offsetm*mseprof(i)))
         sigma_stark(i) = tan(dcon*sigma_stark(i))
         if (sigma_stark(i) .eq. zero) sigma_stark(i) = 0.1*delstark
      end do
      rstarkmin = rstark(1)
      rstarkmax = rstark(imse)

!
!     NEED FOR SCALING EDGE IOTA
!     HERE, SINCE RSTARK IS ORDERED, DATAMAX -> RIGHT OF AXIS
!     AND DATAMIN -> LEFT OF AXIS
!
!       GET RC0MSE (ESTIMATE FOR MAGNETIC AXIS)
!
      rwidth = radial_variance*(rbc(0,1)+abs(rbs(0,1)))

      if (imse .gt. 0) then
         datamax = datastark(imse)
         datamin = datastark(1)
         rc0mse = 0.
         ineg = 1
c                                     !NO ZERO CROSSING: FIND MIN ANYHOW
          if (datamax * datamin.gt.zero) then      !NO ZERO CROSSING: FIND MIN ANYHOW
            datamin = abs(datamin)
            do i = 2,imse
              t1 = abs(datastark(i))
              if (t1.lt.datamin) then
                datamin = t1
                ineg = i
              endif
            end do
             if (ineg.eq.imse) ineg = ineg-1
            ipos = ineg + 1
            goto 310
          else if (( datamax*signiota.lt.zero ) .or.
     >             (datamin*signiota.gt.zero)) then
            datastark = -datastark
          endif

!
!       ALLOW FOR POSSIBLE MULTIPLE ZERO CROSSINGS (WIGGLES) IN DATA
!

         do i = 1, imse
            if (datastark(i)*signiota .le. zero) then
               ineg = i                          !LEFT OF MAGNETIC AXIS
            else
               exit                              !RIGHT OF MAGNETIC AXIS
            endif
         end do
         do i = imse, ineg + 1, -1
            if (datastark(i)*signiota .le. zero) then
               exit                              !LEFT OF MAGNETIC AXIS
            else
               ipos = i                          !RIGHT OF MAGNETIC AXIS
            endif
         end do

  310    continue

         rneg = rstark(ineg)
         rpos = rstark(ipos)
         dataneg = datastark(ineg)
         datapos = datastark(ipos)
      endif                                      !End of if(imse>0)
      if (datapos .ne. dataneg) then
         rc0mse = (datapos*rneg - dataneg*rpos)/(datapos - dataneg)
         rwidth=delstark*abs((rneg-rpos)/(datapos-dataneg))+rwidth
      endif
      if (ipos .gt. ineg + 1) rc0mse = 0.5_dp*(rpos + rneg)

!
!       ESTIMATE MAGNETIC AXIS FROM RAXIS
!
      raxmse = rc0mse
      if (rc0mse.eq.0.0_dp .or. iopt_raxis.ne.1) raxmse = raxis(0,1)
      rstepx0 = 0.005_dp*rwidth

!
!       COMPUTE SPLINES IN R-SPACE FOR MATCHING IOTA(0)
!

      datanorm = zero
      delse1(:imse) = one
      qcalc(:imse) = one/sigma_stark(:imse)**2
      datalsq_s(:imse) = datastark(:imse)*qcalc(:imse)
      scstark = sum(abs(datastark(:imse)/sigma_stark(:imse)))
      datanorm = sum(abs(datastark(:imse)))
      scstark = datanorm/scstark
c04-96        call setspline(rstark,qcalc,datalsq_s,stark_temp,ystark0,
c04-96     >  y2stark0,delse1,0.1*tensi/scstark**2,imse,NATUR)

!
!       DETERMINE NUMBER OF IOTA KNOTS IN SQRT(S)-SPACE
!
      if (isnodes .le. 0) then
         isnodes = min(inode_max,imse + 1)
         isnodes = max(5,isnodes)               !At least 5 spline knots
      endif
      if (isnodes .lt. 5) stop 'MUST PICK ISNODES > 4'
      write (nthreed, *)
     1   'Number of iota-spline knots (in s-space):     ', isnodes

      allocate (nk_ia(isnodes), nk_ib(isnodes), hstark(isnodes),
     1  y2stark(isnodes), ystark(isnodes), sknots(isnodes), stat=istat1)
      if (istat1 .ne. 0) then
         print *, ' ISTAT1 = ', istat1, ' ALLOCATING NK_IA'
         stop
      endif

!
!       COMPUTES NODES IN SQRT(S) SPACE FOR SPLINES
!       THIS ASSUMES A FIXED NO - ISNODES - OF KNOTS
!       THIS MAY NOT PRESERVE ISNODES.
!       ALSO, IT IS NOT NECESSARY TO TAKE EQUALLY SPACED KNOTS IN
!       SQRT(S) SPACE. INDEED, THE FOLLOWING CHOICES ARE POSSIBLE:
!
!       SKNOTS(I) = HNODES*(I-1)   .eq.>   EQUAL-SPACED IN SQRT(S)
!
!       SKNOTS(I) = SQRT(HNODES*(I-1)) .eq.>  EQUAL-SPACED IN S
!
!       DO NOT - UNDER ANY CIRCUMSTANCES - CHANGE THE ARGUMENTS TO
!       THE SPLINT, GETSPLINE, SETUP_INT ROUTINES FROM SQRTS,SHALF
!       TO SQRTS**2, SHALF**2 TO DO S-INTERPOLATION. RATHER, CHANGE
!       SKNOTS (AND PKNOTS) ACCORDING TO THE ABOVE FORMULA. THIS IS
!       ABSOLUTELY CRUCIAL, SINCE ONLY IN SQRT(S) SPACE DO THE
!       FIRST DERIVATIVE BOUNDARY CONDITIONS, d IOTA/d SQRT(S) = 0
!       (SIMILAR FOR P) APPLY AT THE AXIS, S=0.
!
!
      do i = 1, isnodes
         sknots(i) = real(i - 1,rprec)/(isnodes - 1)
      end do

      hstark(:isnodes-1) = sknots(2:isnodes) - sknots(:isnodes-1)


!
!       SET UP DATA ARRAY SCALE FACTORS
!       ACTUAL PRESSURE = PRESPEAK * PFAC * P(INTERNAL)
!       IF (LPOFR) THEN DATA ARE INPUT vs R (REAL SPACE)
!       IF (.NOT.LPOFR),DATA ARE INPUT vs S (FLUX SPACE)
!
      if (itse .eq. 0) call getpresprofile         !!Simulate 'data'
      if (itse .gt. 0) then

         allocate (sthom(itse), delse2(itse), delso2(itse), pcalc(itse),
     1      indexs2(itse), indexu2(itse), stat=istat1)
         if (istat1 .ne. 0) then
            print *, ' ISTAT1 = ', istat1, ' ALLOCATING STHOM'
            stop
         endif


         datathom(:itse) = datathom(:itse)*presfac
         presmax = maxval(datathom(:itse))

!
!       SORT DATA IN ASCENDING ORDER IN R-SPACE (LPOFR) OR S-SPACE(.NOT.LPOFR)
!       AND RE-INDEX DATATHOM, SIGMA_THOM ARRAYS ON ASCENDING RTHOM ARRAY

         datalsq_p(:itse) = datathom(:itse)
         qtemp(:itse) = sigma_thom(:itse)
         call sort_data (rthom, isortp, itse)
         if (lpofr) then
            do i = 1, itse
               datathom(i) = datalsq_p(isortp(i))
               rthom(i) = rthom(i) + pres_offset
               if (rthom(i) .le. zero) then
                  print *, 'Check units of PRES_OFFSET: rthom < 0!'
                  ier = 1
                  return
               endif
               if (datathom(i) .eq. presmax) ipmax = i
               sigma_thom(i) = qtemp(isortp(i))
               if (sigma_thom(i) .ge. cbig) then
                  print *, 'SIGMA_THOM missing'
                  ier = 1
                  return
               endif
               if (sigma_thom(i) .lt. zero) then
                  sigma_thom(i) = abs(sigma_thom(i)*datathom(i))
               else
                  if (sigma_thom(i) .gt. zero) then
                     sigma_thom(i) = presfac*sigma_thom(i)
                  else
                     sigma_thom(i) = p_threshold*presmax
                  endif
               endif
            end do
         else
            do i = 1, itse
               datathom(i) = datalsq_p(isortp(i))
               sthom(i) = rthom(i)
               if (datathom(i) .eq. presmax) ipmax = i
               sigma_thom(i) = qtemp(isortp(i))
               if (sigma_thom(i) .ge. cbig) then
                  print *, 'SIGMA_THOM missing'
                  ier = 1
                  return
               endif
               if (sigma_thom(i) .lt. zero) then
                  sigma_thom(i) = abs(sigma_thom(i)*datathom(i))
               else
                  if (sigma_thom(i) .gt. zero) then
                     sigma_thom(i) = presfac*sigma_thom(i)
                  else
                     sigma_thom(i) = p_threshold*presmax
                  endif
               endif
            end do
         endif

!
!       THROW AWAY NOISY (SMALL) PRESSURE DATA BELOW P_THRESHOLD
!       STARTING FROM PEAK WORKING TO LARGER, SMALLER R
!
         ineg = ipmax
         ipos = ipmax
         do while(ineg.gt.1 .and.
     1            datathom(ineg-1).ge.p_threshold*presmax)
            ineg = ineg - 1
         end do
         do while(ipos.lt.itse .and.
     1            datathom(ipos+1).ge.p_threshold*presmax)
            ipos = ipos + 1
         end do
         itse = ipos - ineg + 1
         ioff = ineg - 1
         do i = 1, itse
            datathom(i) = datathom(ioff+i)
         end do
         do i = 1, itse
            rthom(i) = rthom(ioff+i)
         end do
         do i = 1, itse
            sigma_thom(i) = sigma_thom(ioff+i)
         end do
!
!       COMPUTE PRESSURE AND 1/SIGMA SPLINES IN R-SPACE (OR S-SPACE)
!       a. PRESSURE SPLINE
!
         datanorm = zero
         delse2(:itse) = one
         pcalc(:itse) = one/sigma_thom(:itse)**2
         datalsq_p(:itse) = datathom(:itse)*pcalc(:itse)
         scthom = sum(abs(datathom(:itse)/sigma_thom(:itse)))
         datanorm = sum(abs(datathom(:itse)))
         scthom = datanorm/scthom
         call setspline (rthom, pcalc, datalsq_p, qtemp, ythom0,
     1      y2thom0, delse2, 0.1*tensp/scthom**2, itse, natur)


!
!       FIND PRESSURE PEAK USING SMOOTHED DATA
!
         isamax = maxloc(ythom0(:itse))
         i = isamax(1)
         pthommax = ythom0(i)
         rthompeak = rthom(i)

         ileft = 0                    !Count data points to left of peak
         l4v(:itse) = rthom(:itse) < rthompeak
         do i = 1, itse
            if (l4v(i)) ileft = ileft + 1
         end do

         if (ipnodes .le. 0) then
            ipnodes = max(ileft + 1,itse - ileft)
            ipnodes = min(inode_max,ipnodes)
            ipnodes = max(5,ipnodes)            !At least 5 spline knots
            if (.not.lpprof) ipnodes = 7
         endif
         if (ipnodes < 5) stop 'MUST PICK IPNODES > 4'
         write (nthreed, *)
     1      'Number of pressure-spline knots (in s-space): ', ipnodes

         allocate( nk_pa(ipnodes), nk_pb(ipnodes), ythom(ipnodes),
     1      y2thom(ipnodes), hthom(ipnodes), pknots(ipnodes) )
         if (istat1 .ne. 0) then
            print *, ' ISTAT1 = ', istat1, ' ALLOCATION NK_PA'
            stop
         endif

!
!       COMPUTE NODES IN SQRT(S) SPACE FOR SPLINES
!       (SEE COMMENTS ABOVE PRECEDING SKNOTS(I) CALCULATION)
!
         do i = 1, ipnodes
            pknots(i) = real(i - 1,rprec)/(ipnodes - 1)
         end do
         hthom(:ipnodes-1) = pknots(2:ipnodes) - pknots(:ipnodes-1)
!
!       COMPUTE MINOR RADII FOR DETERMINING PHIEDGE
!
         if (lpofr) then
            rthommax = rthom(itse)
            rthommin = rthom(1)
            presin = datathom(1)
            presout = datathom(itse)
            presmin = min(presin,presout)
            ipresin = 0
            ipresout = 0
            if (presin .eq. presmin) then
               ipresin = 1
               if (presout.le.c1p5*presmin .or. presout<=0.1_dp*presmax)
     1            ipresout = 1
            else
               ipresout = 1
               if(presin.le.c1p5*presmin .or. presin<=0.1_dp*presmax)
     1            ipresin=1
            endif
         else
            ipresin = 0                        !Only use theta=0 in pofs
            ipresout = 1
         endif
      endif                                   !End of if(itse.gt.0) test

!
!       COMPUTE INDICES OF FLUX LOOP MEASUREMENTS NEEDED FOR
!       LOOP SIGNAL MATCHING
!       ALSO COMPUTE CONNECTED EXTERNAL POLOIDAL FLUX ARRAY

      if (.not.lfreeb) then
         nflxs = 0
         nobser = 0
         nobd  = 0
         nbsets = 0
      end if

      nmeasurements = imse + itse + 2 + nflxs   !!1 for diamag, 1 for edge MSE
      do n = 1, nobser
         needflx(n) = idontneed
         iconnect(1,nobd+n) = n            !For outputting absolute flux

         if (lasym) cycle

!        Save Index of up-down symmetric spatial observation points
         do n1 = 1,n-1
           if ((xobser(n1).eq. xobser(n)) .and.
     >         (zobser(n1).eq.-zobser(n))) needflx(n) = n1
         enddo
      end do

      do n = 1, nobd
         dsiext(n) = zero
         if (sigma_flux(n) .lt. zero) sigma_flux(n) =
     1     abs(sigma_flux(n)*dsiobt(n))
         if (sigma_flux(n) .eq. zero) sigma_flux(n) = 0.0001
         do icount = 1, 4
            index = iconnect(icount,n)
            tsign = sign(1,index)
            if(index.ne.0)dsiext(n) = dsiext(n)+psiext(abs(index))*tsign
         end do
      end do

      do n = 1,nflxs
        index = indxflx(n)
        if (index.gt.0) then
          plflux(index) = dsiobt(index) - dsiext(index)        !n-th connected PLASMA flux
          do icount = 1,4
            ind = abs(iconnect(icount,index))
            if ((ind.gt.0).and.(needflx(ind).eq.IDONTNEED))
     1      needflx(ind) = NEEDIT
          enddo
        endif
      enddo

!
!       COMPUTE INDICES OF EXTERNAL BFIELD MEASUREMENTS NEEDED
!       FOR SIGNAL MATCHING
!       FOR MULTIPLE ANGLES AT THE SAME R,Z,PHI LOCATION, IT IS
!       ASSUMED THE LOOP DATA ARE CONSECUTIVELY ORDERED
!
      do n = 1, nbsets
         nmeasurements = nmeasurements + nbfld(n)
         do m = 1,nbcoils(n)
            needbfld(m,n) = IDONTNEED
            if( (m.gt.1).and.(rbcoil(m,n).eq.rbcoil(m-1,n)).and.
     1        (zbcoil(m,n).eq.zbcoil(m-1,n)) )needbfld(m,n)=ISAMECOIL
            if( sigma_b(m,n).lt.zero )sigma_b(m,n) =
     1      abs(sigma_b(m,n) * bbc(m,n))
            if( sigma_b(m,n).eq.zero )sigma_b(m,n) = 0.0001

            if (lasym) cycle
!       CHECK FOR ANTISYMMETRIC SITUATED COIL FOR M1 < M, N <= NSETS
            do n1 = 1, nbsets
               do m1 = 1, m-1
               if( (rbcoil(m1,n1).eq.rbcoil(m,n)) .and.
     1            (zbcoil(m1,n1).eq.-zbcoil(m,n)).and.
     2            (abcoil(m1,n1).eq.abcoil(m,n)) )
     3         needbfld(m,n) = n1 + nbsets*(m1-1)
               enddo
            enddo

         enddo
      enddo

      do n1 = 1, nbsets
        do m1 = 1,nbfld(n1)
          index = indxbfld(m1,n1)
          if( index.gt.0 )then
!m-th PLASMA B-field
            plbfld(index,n1) = bbc(index,n1) - bcoil(index,n1)
            if( needbfld(index,n1).eq.IDONTNEED )
     1      needbfld(index,n1) = NEEDIT
          endif
        enddo
      enddo

      end subroutine fixrecon

      
      subroutine getpresprofile
      use vmec_main
      use vsvd
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i
      real(rprec) :: sigmin
C-----------------------------------------------
 
!
!       WRITE OVER THOMPSON DATA
!
      lpofr = .false.                   !!these data are at s-half nodes
      lpprof = .false.
      itse = 10
      pthommax = datathom(1)                  !!compute in final version
      sigmin = 0.03*pthommax
 
      do i = 1, itse
         rthom(i) = real(i - 1,rprec)/(itse - 1)
         datathom(i) = (pthommax - sigmin)*(1. - rthom(i))**2 + sigmin
         sigma_thom(i) = 0.2*pthommax
      end do
 
      end subroutine getpresprofile

      
      subroutine getgreen
      use vsvd
      use vparams, only: twopi
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: jnm, i
      real(rprec), dimension(10) :: ak, bk
      real(rprec), dimension(9) :: ae, be
      real(rprec), dimension(jngrn) :: ye, yk, sqrt1u
      real(rprec):: dqk2, qk2, eta, alg, sum1, suma, sum2,
     1     sumb, sum3, sumc, sum4, sumd
C-----------------------------------------------
 
      data ak/3.0072519903686507E-04_dp, 3.9684709020989806E-03_dp,
     1   1.0795990490591656E-02_dp, 1.0589953620989356E-02_dp,
     2   7.5193867218083799E-03_dp, 8.9266462945564728E-03_dp,
     3   1.4942029142282098E-02_dp, 3.0885173001899746E-02_dp,
     4   9.6573590301742396E-02_dp, 1.3862943611198872e+0_dp/
      data bk/6.6631752464607272E-05_dp, 1.7216147097986537E-03_dp,
     1   9.2811603829686118E-03_dp, 2.0690240005100891E-02_dp,
     2   2.9503729348688723E-02_dp, 3.7335546682286003E-02_dp,
     3   4.8827155048118076E-02_dp, 7.0312495459546653E-02_dp,
     4   1.2499999999764055e-1_dp, 5.0000000000000000e-1_dp/
      data ae/3.2519201550638976E-04_dp, 4.3025377747931137E-03_dp,
     1   1.1785841008733922E-02_dp, 1.1841925995501268E-02_dp,
     2   9.0355277375409049E-03_dp, 1.1716766944657730E-02_dp,
     3   2.1836131405486903E-02_dp, 5.6805223329308374E-02_dp,
     4   4.4314718058336844E-1_dp/
      data be/7.2031696345715643E-05_dp, 1.8645379184063365E-03_dp,
     1   1.0087958494375104E-02_dp, 2.2660309891604169E-02_dp,
     2   3.2811069172721030E-02_dp, 4.2672510126591678E-02_dp,
     3   5.8592707184265347E-02_dp, 9.3749995116366946E-02_dp,
     4   2.4999999999746159E-1_dp/
 
!
!       Compute "Green's Functions" for Poloidal Flux, 2*pi*R*A-sub-phi,
!       BR, and BZ at point (XT,ZT) due to unit current (mu0*I = 1) at (XS,ZS) ...
!       modified to interpolate on k**2 - 3-34-92 - sph
!
      jnm = jngrn - 1
      odqk2 = (jnm)
      dqk2 = 1.0_dp/odqk2
      do i = 2, jnm
         qk2 = dqk2*(i - 1)
         qsq(i) = qk2
         eta = 1 - qk2
         alg = log(eta)
         sum1 = ((((ak(1)*eta+ak(2))*eta+ak(3))*eta+ak(4))*eta+ak(5))*
     1      eta + ak(6)
         suma = (((sum1*eta + ak(7))*eta+ak(8))*eta+ak(9))*eta + ak(10)
         sum2 = ((((bk(1)*eta+bk(2))*eta+bk(3))*eta+bk(4))*eta+bk(5))*
     1      eta + bk(6)
         sumb = (((sum2*eta + bk(7))*eta+bk(8))*eta+bk(9))*eta + bk(10)
         yk(i) = suma - alg*sumb
         sum3 = (((ae(1)*eta+ae(2))*eta+ae(3))*eta+ae(4))*eta
         sumc = (((((sum3 + ae(5))*eta+ae(6))*eta+ae(7))*eta+ae(8))*eta+
     1      ae(9))*eta
         sum4 = (((be(1)*eta+be(2))*eta+be(3))*eta+be(4))*eta
         sumd = (((((sum4 + be(5))*eta+be(6))*eta+be(7))*eta+be(8))*eta+
     1      be(9))*eta
         ye(i) = sumc - alg*sumd + 1
         yf(i) = ((1 + eta)*yk(i)-2*ye(i))/qk2
      end do
      ye(1) = 0.25_dp*twopi
      ye(jngrn) = 1
      yk(1) = ye(1)
      yk(jngrn) = 2*yk(jnm) - yk(jngrn-2)
      yf(1) = 0.
      yf(jngrn) = 2*yf(jnm) - yf(jngrn-2)
      qsq(1) = 0
      qsq(jngrn) = 1
 
      sqrt1u = sqrt(qsq(:jngrn))/twopi
c                                      !Factor of 1/2 from sqrt(4*xs*xt)
      yek(:jngrn) = 0.5_dp*sqrt1u*(ye(:jngrn)-yk(:jngrn))
      yeq(:jngrn) = 0.25_dp*qsq(:jngrn)*sqrt1u*ye(:jngrn)
c                                 !Factor of 2 absorbed by sqrt(4 xt xs)
      yf(:jngrn) = twopi*sqrt1u*yf(:jngrn)
      dyek(:jnm) = (yek(2:jnm+1)-yek(:jnm))*odqk2
      dyeq(:jnm) = (yeq(2:jnm+1)-yeq(:jnm))*odqk2
      dyf(:jnm) = (yf(2:jnm+1)-yf(:jnm))*odqk2
 
      end subroutine getgreen

      
      subroutine getlim
      use vmec_main
      use realspace
      use vsvd
      implicit none
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: resup = 3.0_dp
      real(rprec), parameter :: resdn = 0.5_dp
      real(rprec), parameter :: eps = 0.005_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ntheta_2pi, nphi_plane, i,  
     1   nthtot, iexpand, ishrink, ionlim, n,
     2   limpts, nonlim, nexpand, nshrink, ilim0, nlim0
      real(rprec), dimension(2*ntheta1) ::
     1   rbdy, zbdy, rubdy, zubdy
      real(rprec) :: fshrink, distmax, fexpand
C-----------------------------------------------
 
c
c     DETERMINES WHEN PLASMA TOUCHES LIMITER
c     USE DOUBLE THE NO. OF THETA POINTS FOR INCREASED RESOLUTION
c
      ntheta_2pi = ntheta1
      nphi_plane = 1                 !Pick a phi plane (phi = 0 for now)
 
      rbdy(:ntheta3*2-1:2) = r1(ns*nphi_plane:ns*ntheta3*nphi_plane:ns*
     1   nphi_plane,0) + r1(ns*nphi_plane:ns*ntheta3*nphi_plane:ns*
     2   nphi_plane,1)
      zbdy(:ntheta3*2-1:2) = z1(ns*nphi_plane:ns*ntheta3*nphi_plane:ns*
     1   nphi_plane,0) + z1(ns*nphi_plane:ns*ntheta3*nphi_plane:ns*
     2   nphi_plane,1)
      rubdy(:ntheta3*2-1:2) = ru(ns*nphi_plane:ns*ntheta3*nphi_plane:ns*
     1   nphi_plane,0) + ru(ns*nphi_plane:ns*ntheta3*nphi_plane:ns*
     2   nphi_plane,1)
      zubdy(:ntheta3*2-1:2) = zu(ns*nphi_plane:ns*ntheta3*nphi_plane:ns*
     1   nphi_plane,0) + zu(ns*nphi_plane:ns*ntheta3*nphi_plane:ns*
     2   nphi_plane,1)
 
 
      if (.not.lasym) then
!     FOR NOW, THIS ONLY WORKS FOR NZETA=1 (PHI=0 PLANE)
!     TO EXTEND TO OTHER PHI PLANES, MUST USE IREFLECT(NZETA)
      do i = 1, ntheta_2pi - ntheta2
         rbdy(2*(ntheta2+i)-1) = rbdy(2*(ntheta1-ntheta2-i)+3)
      end do
      do i = 1, ntheta_2pi - ntheta2
         zbdy(2*(ntheta2+i)-1) = -zbdy(2*(ntheta1-ntheta2-i)+3)
      end do
      do i = 1, ntheta_2pi - ntheta2
         rubdy(2*(ntheta2+i)-1) = -rubdy(2*(ntheta1-ntheta2-i)+3)
      end do
      do i = 1, ntheta_2pi - ntheta2
         zubdy(2*(ntheta2+i)-1) = zubdy(2*(ntheta1-ntheta2-i)+3)
      end do
      end if
 
!
!     FIND EVEN INDEXED POINTS BY INTERPOLATION
!
      nthtot = 2*ntheta_2pi
      rbdy(nthtot) = .5_dp*(rbdy(1)+rbdy(nthtot-1))
      zbdy(nthtot) = .5_dp*(zbdy(1)+zbdy(nthtot-1))
      rubdy(nthtot) = .5_dp*(rubdy(1)+rubdy(nthtot-1))
      zubdy(nthtot) = .5_dp*(zubdy(1)+zubdy(nthtot-1))
      rbdy(2:(ntheta_2pi-1)*2:2) = .5_dp*(rbdy(3:ntheta_2pi*2-1:2) +
     1      rbdy(:ntheta_2pi*2-3:2))
      zbdy(2:(ntheta_2pi-1)*2:2) = .5_dp*(zbdy(3:ntheta_2pi*2-1:2) +
     1      zbdy(:ntheta_2pi*2-3:2))
      rubdy(2:(ntheta_2pi-1)*2:2) = .5_dp*(rubdy(3:ntheta_2pi*2-1:2)+
     1   rubdy(:ntheta_2pi*2-3:2))
      zubdy(2:(ntheta_2pi-1)*2:2) = .5_dp*(zubdy(3:ntheta_2pi*2-1:2)+
     1   zubdy(:ntheta_2pi*2-3:2))
 
      fshrink = 0.0_dp
      distmax = sum(rbdy(:nthtot)**2) + sum(zbdy(:nthtot)**2)
      fexpand = distmax
      iexpand = 0
      ishrink = 0
      ionlim = 0
 
      do n = 1, nlim
         limpts = limitr(n)
         call cauchy (rbdy, zbdy, rubdy, zubdy, rlim(:,n), zlim(:,n),
     1      reslim(:,n), seplim(:,n), distmax, nthtot, limpts)

           do i = 1,limpts
c       LIMITER POINT ON PLASMA
            if( (abs(reslim(i,n)-resdn).lt.eps) )then
!    .gt.      .or. (abs(reslim(i,n)).gt.resup) )then
              ionlim = i
              nonlim = n
c       LIMITER POINT OUTSIDE PLASMA
            else if( reslim(i,n).lt.RESDN )then
              if( seplim(i,n).le.fexpand )then
                fexpand = seplim(i,n)
                iexpand = i
                nexpand = n
              endif
c       LIMITER POINT INSIDE PLASMA
            else if( reslim(i,n).ge.RESDN )then
              if( seplim(i,n).gt.fshrink )then
                fshrink =  seplim(i,n) 
                ishrink = i
                nshrink = n
              endif
            endif
          enddo
      end do
 
c
c       LOGIC: IF THERE IS A LIMITER POINT INSIDE PLASMA, THEN MUST
c       SHRINK CONTOUR. OTHERWISE, IF THERE IS AT LEAST ONE LIMITER
c       POINT ON CONTOUR, AND ALL THE REST OUTSIDE, DO NOT CHANGE ANYTHING.
c       FINALLY, IF ALL LIMITER POINTS ARE OUTSIDE PLASMA, EXPAND PLASMA
c       TO OSCULATE WITH LIMITER
c
      if (ishrink .gt. 0) then
         gphifac = -sqrt(fshrink)
         ilim0 = ishrink
         nlim0 = nshrink
      else if (ionlim .gt. 0) then
         gphifac = 0.
         ilim0 = ionlim
         nlim0 = nonlim
      else if (iexpand .gt. 0) then
         gphifac = sqrt(fexpand)
         ilim0 = iexpand
         nlim0 = nexpand
      endif
 
      dlim_min = gphifac
      rlim_min = rlim(ilim0,nlim0)
      zlim_min = zlim(ilim0,nlim0)
 
c     overall damping in time, rsfac/sqrt(rsfac) = sqrt(rsfac)
      gphifac = gphifac/r01
      if (abs(gphifac) .gt. 0.04*one)
     1   gphifac = 0.04*gphifac/abs(gphifac)
      gphifac = 0.20*gphifac/sqrt(rsfac)
 
      end subroutine getlim


      subroutine cauchy(rbdy, zbdy, rubdy, zubdy, rlim, zlim, residue, 
     1   sep, distmax, ntheta, nlim)
      use vparams, only: twopi, dp, rprec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: ntheta, nlim
      real(rprec), intent(in) :: distmax
      real(rprec), dimension(ntheta), intent(in) :: 
     1     rbdy, zbdy, rubdy, zubdy
      real(rprec), dimension(nlim) :: rlim, zlim, residue, sep
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero=0, p5=0.5_dp, two=2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: n, i, imin, imax
      real(rprec), dimension(ntheta) :: dsq, dsepdu
      real(rprec), dimension(ntheta) :: x1u, y1u
      real(rprec) :: delu, dmin, delta_d, alpha, gam0
C-----------------------------------------------
c       Check that the points (rlim(i),zlim(i)) are inside boundary surface
c       using Cauchys theorem in "complex"-plane (for a fixed
c       toroidal plane, nphi=const. It is assumed that rbdy, zbdy are
c       extended around the full interval, 0-2pi, in theta.
c       with rbdy(1) = rbdy(ntheta+1) (i.e., ntheta intervals)
c
c       Because of numerical inaccuracies, instead of testing on
c       res = 0 (outside), res = 1 (inside), we use the test:
c       res >= .5, inside;  res < .5, outside
c**********************************************************************
c
c       LOCAL VARIABLE ARRAYS
c
c       dsq:    Distance squared between limiter point and plasma boundary
c       sep:    Minimum dsq for each limiter point
c    dsepdu:    .5* d(dsq)/d(theta)
c   residue:    Contour integral of 1/(X-rlim)in complex X=(R,Z) plane
c
      delu = twopi/ntheta
      dmin = 1.E-20_DP*distmax
 
      do n = 1, nlim
         residue(n) = zero
         x1u = rbdy(:ntheta) - rlim(n)
         y1u = zbdy(:ntheta) - zlim(n)
         dsq(:ntheta) = x1u*x1u + y1u*y1u
         dsepdu(:ntheta) = x1u*rubdy(:ntheta) + y1u*zubdy(:ntheta)
         residue(n) = residue(n) + sum((x1u*zubdy(:ntheta) - 
     1      y1u*rubdy(:ntheta))/(dsq(:ntheta)+dmin))
 
         residue(n) = residue(n)/ntheta
!
!        Find actual minimum distance from nth limiter point to boundary
!
         sep(n) = distmax
         do i = 1,ntheta
            if( dsq(i).le.sep(n) )then
               imin = i
               sep(n) = dsq(i)
            endif
         enddo
 
!        gamu = two*abs(dsepdu(imin))*delu
 
         if (dsepdu(imin) .le. zero) then
            imax = 1 + mod(imin,ntheta)
         else
            imax = imin
            imin = imax - 1
            if (imin .eq. 0) imin = ntheta
         endif
 
         delta_d = two*(dsepdu(imax)-dsepdu(imin))
         alpha = delta_d/delu
!        gamu = gamu/delta_d
!        sep(n) = sep(n) - p5*alpha*gamu**2
         if (alpha .ne. zero) 
     1      gam0 = 0.5_dp - (dsq(imax)-dsq(imin))/(alpha*delu**2)
         sep(n) = dsq(imin) - p5*alpha*(gam0*delu)**2
         if (sep(n) .lt. zero) sep(n) = zero
 
      end do
 
      end subroutine cauchy

           
      subroutine newprofil(phipog)
      use vmec_main
      use vacmod
      use realspace
      use vforces, lu => czmn, lv => crmn
      use vsvd
      use vspline
      use xstuff
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(*) :: phipog
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: icount, inodes, js
      integer, save :: idata(jchix), isize(jchix)       !INCREMENTAL UPDATES
      real(rprec), dimension(ipnodes + isnodes) :: datalsq, wten
      real(rprec), dimension(1+nrzt) :: r12sqr
      real(rprec) :: amat_lsq(isnodes+ipnodes,isnodes+ipnodes),
     1   djac_p(ipnodes,nmeasurements), djac_i(isnodes,nmeasurements)
      real(rprec), dimension(nmeasurements) :: datainput
      real(rprec) :: treconon, delt1, pfac0, aminout,
     1   aminin, ymin, pfactor, treconoff

C----------------------------------------------- 
!     IDATA: GIVES STARTING INDEX IN DATAINPUT ARRAY
!            FOR STORING EACH OF THE DATA TYPES
!     ISIZE: GIVES NUMBER OF DATA STARTING AT IDATA
!     INDEXING OF DJAC ARRAY (ASSUMES K STARTS AT 0)
 
      inodes = isnodes + ipnodes
 
      call second0 (treconon)
 
!     Unfreeze magnetic axis 
      if (iresidue.eq.0 .and. fsq*1.e6_dp.lt.one) iresidue = 1
      delt1 = one/real(ipedsvd,rprec)
      if (iresidue .eq. 0) delt1 = one
 
!
!       COMPUTE AVERAGE RADIAL FORCE BALANCE CONSTRAINT
!
      call radfor (pfac0)
      pfac = pfac + delt1*(pfac0 - pfac)
!
!       UPDATE PHI SCALE FACTOR (PHIFAC)
!       SCALE TOROIDAL FLUX TO MATCH PRESSURE WIDTH OR LIMITER
!
      aminout = max(rthommax,rstarkmax) - r00
      aminin = r00 - min(rthommin,rstarkmin)
      apres = (aminin*ipresin + aminout*ipresout)/(ipresin +ipresout)
      aminor = ((r00 - rinner)*ipresin + (router - r00)*ipresout)/
     1       (ipresin + ipresout)
      if (imatch_phiedge.ne.1 .and. ivac.gt.1 .or. imatch_phiedge.eq.3 
     1   .and. (.not.lfreeb)) then
         call newphi (phipog)
         call gettflux
      endif
 
      icount = 0
 
      if (.not.(mod(iter2 - iter1,ipedsvd).ne.0 .and. iequi.eq.0
     1     .and. iresidue.gt.0)) then
 
!
!       SETUP COMMON BLOCKS FOR FLUX-MATCHING ROUTINES
!

         if (iphidiam + nflxs + nbfldn.gt.0 .or. iequi.gt.0) then
            r12sqr(2:nrzt) = sqrt(armn_o(2:nrzt))
            call flux_init (phipog)
         endif
 
!
!       COMPUTE MATRIX ELEMENTS FOR THOMPSON SCATTERING DATA
!
         idata(ithom0) = icount + 1
         call getthom(djac_i(1,idata(ITHOM0)), djac_p(1,idata(ITHOM0)),
     1     datainput(idata(ITHOM0)), r1(1:,0), r1(1:,1), isize(ITHOM0))
         icount = icount + isize(ithom0)
 
!
!       COMPUTE MOTIONAL STARK EFFECT. THIS CALL ALSO INITIALIZES
!       THE ALSQ, DATALSQ ARRAYS AND SETS UP THE SPLINE NODES.
!
         idata(istark0) = icount + 1
         call getmse(djac_i(1,idata(ISTARK0)), djac_p(1,idata(ISTARK0)),
     1     datainput(idata(ISTARK0)), r1(1:,0), r1(1:,1),lu,
     2     lu(1+nrzt), zu(1:,0), zu(1:,1), phipog, isize(ISTARK0))
         icount = icount + isize(istark0)
 
!
!       COMPUTE MATRIX ELEMENTS FOR DIAMAGNETIC FLUX LOOP
!
         idata(idiam0) = icount + 1
         call getdiam (djac_i(1,idata(idiam0)), djac_p(1,idata(idiam0))
     1      , datainput(idata(idiam0)), isize(idiam0))
         icount = icount + isize(idiam0)
 
!
!       COMPUTE MATRIX ELEMENTS FOR EXTERNAL POLOIDAL FLUXES
!
         idata(iflxs0) = icount + 1
         call getflux (djac_i(1,idata(iflxs0)), djac_p(1,idata(iflxs0))
     1      , datainput(idata(iflxs0)), r12sqr, clmn_e(1), clmn_o(1), 
     2      blmn_o(1), armn_o(1), blmn_e(1), azmn_o(1), isize(iflxs0))
         icount = icount + isize(iflxs0)
 
!
!       COMPUTE MATRIX ELEMENTS FOR EXTERNAL MAGNETIC FIELD MATCHING
!
         idata(ibrzfld) = icount + 1
         call getbfld (djac_i(1,idata(ibrzfld)), djac_p(1,idata(ibrzfld)
     1      ), datainput(idata(ibrzfld)), r12sqr, azmn_o(1), blmn_o(1), 
     2      clmn_e(1), clmn_o(1), armn_o(1), blmn_e(1), isize(ibrzfld))
         icount = icount + isize(ibrzfld)
 
!
!       SQUARE DATA MATRIX ELEMENTS AND STORE IN ALSQ
!
 
         if (icount .gt. nmeasurements) stop 'icount>nmeasurements'
         if (iequi .eq. 0) then
 
            call sgemvmm (djac_i, djac_p, amat_lsq, datainput, datalsq, 
     1         wten, icount, isnodes, ipnodes, inodes)
 
!
!       COMPUTE IOTA, PRESSURE SPLINE COEFFICIENTS
!
            call set_dual (datalsq, hstark, ystark, y2stark, hthom, 
     1         ythom, y2thom, wten, amat_lsq, isnodes, ipnodes, inodes)
 
            if (.not.lpprof) then
               ymin = minval(ythom(1:ipnodes))
               ythom(:ipnodes) = ythom(:ipnodes) - ymin
            endif
 
!
!       COMPUTE IOTA, PRESSURE AT R(js) FROM SPLINED INPUT
!       DATA ALONG THE MIDPLANE
!
            call splint (sknots, ystark, y2stark, isnodes, sqrts, 
     1         isplinef, zero, ns)
            call splint (sknots, ystark, y2stark, isnodes, shalf(2), 
     1         isplineh(2), zero, ns1)
            call splint (pknots, ythom, y2thom, ipnodes, sqrts, psplinef
     1         , zero, ns)
            call splint (pknots, ythom, y2thom, ipnodes, shalf(2), 
     1         psplineh(2), zero, ns1)
 
            pfactor = dmu0*pthommax             !!*pfac moved to getthom
            do js = 1,ns
              isplinef(js) = isplinef(js) - iotaf(js)
              isplineh(js) = isplineh(js) - iotas(js)
              psplinef(js) = pfactor*psplinef(js) - presf(js)
              psplineh(js) = pfactor*psplineh(js) - mass(js)
            end do  
 
         endif                     ! iequi>0
!
!       COMPUTE CHISQ
!
         call chisq (djac_i, djac_p, datainput, idata, isize, icount)

      endif                        ! mod(iter2-iter1,ipedsvd) == 0

      if (iequi .eq. 0) then
!
!       UPDATE PRESSURE SPLINE AND ESTABLISH INTERNAL
!       (CODE) UNITS OF PRESSURE. P(real) = dmu0 * pthommax * P(splined)
!       WHERE P(code-units) = mu0 * P(real)
!       SMOOTH TIME VARIATION OF PROFILES
!
         do js = 1,ns
           iotaf(js) = iotaf(js) + delt1*isplinef(js)
           iotas(js) = iotas(js) + delt1*isplineh(js)
           presf(js) = presf(js) + delt1*psplinef(js)
           mass(js)  = mass(js)  + delt1*psplineh(js)
         end do
      endif
 
!
!     STORE CHISQ
!
      call store_chisq
 
!
!       OPTIMIZE MAGNETIC AXIS POSITION BY MINIMIZING RMS ERROR
!       RETURNS RAXMSE AS UPDATED GUESS FOR NEW AXIS POSITION
!       TO BE USED IN SUBROUTINE RESIDUE
!
      call axisopt (fsq, r00, iresidue, ivac)
 
!       Compute force to fix axis at RAXMSE
      grmse = -0.05*(r00 - raxmse)
 
      call second0 (treconoff)
      timer(6) = timer(6) + (treconoff - treconon)
 
      end subroutine newprofil
      

      subroutine flux_init(phipog)
      use vmec_main
      use vmec_params, only: signgs
      use vforces, only : r12=>armn_o, gsqrt=>azmn_o, orsq=>blmn_o
      use vsvd
      use realspace
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(*) :: phipog
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: p5 = 0.5_dp, c1p5 = 1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: l, js
      real(rprec), external :: dot_g
C-----------------------------------------------
 
!
!     COMPUTE OBSERVATION POINT - INVARIANT FUNCTIONS OF RADIUS
!     CURRENT = PHIP * INTEGRAL(0,2pi)[ guu / gsqrt]   (on full mesh)
!     RM2     = < R**(-2) >
!     VRM2    = V` * RM2    (V` = 2pi * VP)
!     ORSQ    = SQRT(G) * R**(-2)   (on half mesh)
!
!     MUST HAVE GONE THROUGH NEWPROFILE DETERMINATION OF IOTA AT
!     LEAST ONCE, OTHERWISE IOTAS IS UNDEFINED!
!
      if (iresidue .le. 0) return 
      current(1) = zero
      presint(1) = one
 
      do l = 2,nrzt-1
        orsq(l) = p5*( phipog(l) + phipog(l+1) ) *
     1  (ru0(l)*ru0(l) + zu0(l)*zu0(l))
      enddo
      do l = ns,nrzt,ns
        orsq(l) = ( c1p5*phipog(l) - p5*phipog(l-1) ) *
     1  (ru0(l)*ru0(l) + zu0(l)*zu0(l))
      enddo

      do js = 2, ns
         current(js) = twopi*DOT_G(nznt,orsq(js),ns,wint(js),ns)
         presint(js) = one
      end do
 
      do l = 2, nrzt
         orsq(l) = gsqrt(l)/r12(l)**2
      end do
      do js = 2, ns
         vrm2(js) = twopi*DOT_G(nznt,orsq(js),ns,wint(js),ns)
         rm2(js) = vrm2(js)/(twopi*signgs*vp(js))
         ovrm2(js) = one/vrm2(js)
         ochip(js) = one/(phip(js)*iotas(js))
         presph(js) = presf(js) - presf(js - 1)
      end do
 
      end subroutine flux_init

      
      subroutine getbfld(amat_i, amat_p, data_array, r12sqr, 
     1  gsqrt, orsq, gobsr1, gobsz1, r12, z12, kcbfld)
      use vmec_main
      use vmec_params, only: signgs
      use vsvd
      use realspace, only: wint
      use vspline, only: hthom, hstark
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer kcbfld
      real(rprec), dimension(isnodes,*) :: amat_i
      real(rprec), dimension(ipnodes,*) :: amat_p
      real(rprec), dimension(*) :: data_array
      real(rprec), dimension(nrzt) :: r12sqr,
     1  gsqrt, orsq, gobsr1, gobsz1, r12, z12
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: js, n1, m1, iobs, isym, iobsmax, iloop, 
     1   msym, nsym, indexbfld, l, lk
      real(rprec), dimension(:), allocatable :: gobsz2, gobsr2
      real(rprec):: tpisg, sumir, sumiz, wscaleb, sumpr,
     1    sumpz, deltab, plasbfld, coscoil, sincoil, t2, t1
C----------------------------------------------
       
!     IRESIDUE > 0, OTHERWISE FLUX_INIT NOT CALLED YET
      kcbfld = 0
      if (iresidue.le.0 .or. (nbfldn.eq.0 .and. iequi.eq.0))return
      
      allocate( gobsz2(nrzt), gobsr2(nrzt), stat=l)
      if (l .ne. 0) stop 'allocation problem in getbfld'
!
!       COMPUTE "GREEN'S FUNCTION" KERNEL ONLY FOR NEEDED OBSERVATION POINTS
!       (OR FOR FINAL OUTPUT IF IEQUI=1)
!       GOBSR1 = BR, GOBSZ1 = BZ 
!       IF A (BR,BZ) OR (BRHO,BTHETA) PAIR, JUST CALL GRNBFLD ONCE
!
      if (any(rm2(2:ns) .eq. zero)) stop 'rm2 = 0'
      tpisg = twopi * signgs                !Positive volume integral
      do n1 = 1, nbsets
         do m1 = 1, nbcoils(n1)
            isym = needbfld(m1,n1)
            if (isym.eq.needit .or. iequi.eq.1) then
               call grnbfld(r12sqr,r12,z12,gobsr1,gobsz1,nrzt,m1,n1)
               do l = 2,nrzt
                 gobsr2(l) = gobsr1(l)*orsq(l)
                 gobsr1(l) = gobsr1(l)*gsqrt(l)
                 gobsz2(l) = gobsz1(l)*orsq(l)
                 gobsz1(l) = gobsz1(l)*gsqrt(l)
               end do
!
!       DO INTEGRAL OVER ANGLES (ALL INTEGRALS ARE FROM THETA=0,TWOPI)
!
               do js = 2, ns
                  sumir = zero
                  sumiz = zero
                  sumpr = zero
                  sumpz = zero
                  do lk = js,nrzt,ns
                    sumir = sumir + gobsr2(lk)*wint(lk)
                    sumpr = sumpr + gobsr1(lk)*wint(lk)
                    sumiz = sumiz + gobsz2(lk)*wint(lk)
                    sumpz = sumpz + gobsz1(lk)*wint(lk)
                  enddo
                  imb(js,m1,n1,1) = tpisg*sumir
                  pmb(js,m1,n1,1) = (-tpisg*sumpr) + imb(js,m1,n1,1)
     1               /rm2(js)
                  imb(js,m1,n1,2) = tpisg*sumiz
                  pmb(js,m1,n1,2) = (-tpisg*sumpz) + imb(js,m1,n1,2)
     1               /rm2(js)
               end do

            else if (isym .eq. ISAMECOIL) then       !Same coil position as previous coil
              do js = 2,ns
                imb(js,m1,n1,1) = imb(js,m1-1,n1,1)
                pmb(js,m1,n1,1) = pmb(js,m1-1,n1,1)
                imb(js,m1,n1,2) = imb(js,m1-1,n1,2)
                pmb(js,m1,n1,2) = pmb(js,m1-1,n1,2)
              enddo
            endif
          enddo   !m1
        enddo   !n1

!
!       CHECK FOR SYMMETRIC COIL (MAY BE IN DIFFERENT COIL SET,
!       SO HAD TO MOVE OUT OF M1,N1 LOOP ABOVE)
!
        do n1 = 1,nbsets
          do m1 = 1,nbcoils(n1)
            isym = needbfld(m1,n1)
            if (isym .ge. ISYMCOIL) then
              msym = 1 + (isym-1)/nbsets
              nsym = 1 + mod(isym-1,nbsets)
              do js = 2,ns                     !BR(-Z) = -BR(Z), BZ(-Z) = BZ(Z)
                imb(js,m1,n1,1) =-imb(js,msym,nsym,1)
                pmb(js,m1,n1,1) =-pmb(js,msym,nsym,1)
                imb(js,m1,n1,2) = imb(js,msym,nsym,2)
                pmb(js,m1,n1,2) = pmb(js,msym,nsym,2)
              enddo
            endif
          enddo
        enddo

!
!       COMPUTE SPLINE MATRIX ELEMENTS BY INTEGRATING OVER RADIUS
!
        do 2000 iloop = 0,iequi                !iequi = 0 normally, = 1 at end
          do n1 = 1, nbsets
            iobsmax = nbfld(n1)
            if (iloop .eq. 1) iobsmax = nbcoils(n1)
            if (iobsmax .gt. 0) then
              do 1000 iobs = 1, iobsmax
                indexbfld = indxbfld(iobs,n1)
                if (iloop .eq. 1) indexbfld = iobs
                if (indexbfld .le. 0) goto 1000
                coscoil = cos( abcoil(indexbfld,n1) )
                sincoil = sin( abcoil(indexbfld,n1) )
                do js = 2,ns
                  pmb(js,0,n1,1) = ochip(js) *
     >            (pmb(js,indexbfld,n1,1)*coscoil +
     >             pmb(js,indexbfld,n1,2)*sincoil)
                  imb(js,0,n1,1) = ovrm2(js) *
     >            (imb(js,indexbfld,n1,1)*coscoil +
     >             imb(js,indexbfld,n1,2)*sincoil) 
                end do

                     if (iloop .eq. 0) then
                        deltab = plbfld(indexbfld,n1)
                        kcbfld = kcbfld + 1
 
                        call splinint (imb(1,0,n1,1), current, 
     1                     amat_i(1,kcbfld), hstark, u_ib, u1_ib, 
     2                     w_ib, w1_ib, nk_ib, isnodes, intder, ns)
 
                        call splinint (pmb(1,0,n1,1), presint, 
     1                     amat_p(1,kcbfld), hthom, u_pb, u1_pb, w_pb, 
     2                     w1_pb, nk_pb, ipnodes, intder, ns)
 
                        wscaleb = one/sigma_b(indexbfld,n1)
                        data_array(kcbfld) = wscaleb*deltab
                        t2 = dmu0*pthommax      !!*pfac moved to getthom
 
                        amat_i(:,kcbfld) = wscaleb*amat_i(:,kcbfld)
                        wscaleb = wscaleb*t2
                        amat_p(:,kcbfld) = wscaleb*amat_p(:,kcbfld)
 
                     else      !Store plasma fluxes in EXTFLX for output
                        plasbfld = zero
                        do js = 2, ns
                           t1 = current(js)*iotaf(js) - current(js-1)*
     1                        iotaf(js - 1)
                           plasbfld = plasbfld + pmb(js,0,n1,1)*
     1                        presph(js) + imb(js,0,n1,1)*t1
                        end do
                        plbfld(iobs,n1) = plasbfld
                     endif
 1000         continue
            endif
          enddo       !n1
 2000   continue

      deallocate( gobsz2, gobsr2, stat=l)

      end subroutine getbfld

      
      subroutine getdiam(amat_i, amat_p, data_array, kcdiam)
      use vmec_main
      use vmec_params, only: signgs
      use realspace
      use vforces, only : r12=>armn_o, ru12=>azmn_e
      use vsvd
      use vspline
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer kcdiam
      real(rprec), dimension(isnodes,*) :: amat_i
      real(rprec), dimension(ipnodes,*) :: amat_p
      real(rprec), dimension(*) :: data_array
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: ilimit = 2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: js, lk, l
      real(rprec), dimension(ns) :: gp, gi, gip
      real(rprec), dimension(isnodes) :: amat2_i
      real(rprec) :: wdiam, z12, tv, ti, t2, sum1
C-----------------------------------------------
 
 
      kcdiam = 0
      if (iphidiam.eq.0 .or. iresidue.lt.ilimit) return 
!
!       COMPUTE FIT TO DIAMAGNETIC SIGNAL, USING EQUILIBRIUM RELATION
!       (modified 7/96 by SPH)
!
!       PHI-DIAMAG = 2*pi*INT[ Gp dp/ds + Gi d(<Bu>)/ds ]
!
!       where
!
!       Gp = Ru * Z * <sqrt(g)> /(R * phip)
!       Gi = Ru * Z * iota / R
!
 
      kcdiam = kcdiam + 1
c-7/96  dNewPhiedge = signgs*twopi*hs*SSUM_1(ns1,phip(2),1)
c-7/96  VacPhiedge  = signgs*bsubvvac*hs*SSUM_1(ns1,vrm2(2),1)
c-7/96  delphid0    = VacPhiedge - dNewPhiedge
 
      wdiam = one/sigma_delphid
      gp(1) = zero
      gi(1) = zero
      do js = 2, ns
         gp(js) = zero
         do lk = 1, nznt
            l = js + ns*(lk - 1)
            z12 = .5_dp*(z1(l,0)+z1(l-1,0)+shalf(l)*(z1(l,1)+z1(l-1,1)))
            gp(js) = gp(js) + ru12(l)*z12/r12(l)*wint(l)
         end do
      end do
!
!       NOTE: gip terms comes from linearizing the iota*d/ds[current*iota]
!             terms
!
      do js = 2, ns
         tv = twopi*vp(js)/phip(js)
         ti = -gp(js)*signgs*wdiam
         gi(js) = ti*iotas(js)
         gp(js) = -gp(js)*tv*wdiam
         gip(js) = ti*(current(js)*iotaf(js)-current(js-1)*iotaf(js-1))
      end do
 
      call splinint (gi, current, amat_i(1,kcdiam), hstark, u_ib, u1_ib
     1   , w_ib, w1_ib, nk_ib, isnodes, intder, ns)
 
      call splinint (gip(2), current(2), amat2_i, hstark, u_ia, u1_ia, 
     1   w_ia, w1_ia, nk_ia, isnodes, intfun, ns1)
 
      call splinint (gp, presint, amat_p(1,kcdiam), hthom, u_pb, u1_pb, 
     1   w_pb, w1_pb, nk_pb, ipnodes, intder, ns)
 
      amat_i(:isnodes,kcdiam) = amat_i(:isnodes,kcdiam) + amat2_i(:
     1   isnodes)
 
      t2 = dmu0*pthommax                        !!*pfac moved to getthom
      amat_p(:ipnodes,kcdiam) = t2*amat_p(:ipnodes,kcdiam)
      sum1 = sum(iotas(2:ns)*gip(2:ns))
      data_array(kcdiam) = wdiam*phidiam + sum1
      if (iequi .eq. 0) then
!
!       Eliminate p variation until well-converged
!
!@        do i = 1,ipnodes
!@          data_array(kcdiam) = data_array(kcdiam) -
!@     >    amat_p(i,kcdiam)*ythom(i)
!@          amat_p(i,kcdiam) = 0.
!@        end do
 
!
!       FINAL OUTPUT (ALSO USE FOR DEBUGGING)
!
      else
!
!       Integrate by parts
!
         delphid = gp(ns)*presf(ns) + gi(ns)*iotaf(ns)*current(ns) - 
     1      gp(2)*presf(1) - gi(2)*current(1)*iotaf(1)
         do js = 2, ns1
            delphid = delphid - presf(js)*(gp(js+1)-gp(js)) - iotaf(js)*
     1         current(js)*(gi(js+1)-gi(js))
         end do
         delphid = delphid/wdiam
      endif
!@        do js = 2,ns
!@        end do
!@
!@        sumi = sum(amat_i(:isnodes,kcdiam)*ystark(:isnodes))
!@     >       - sum(amat2_i(isnodes)*ystark(:isnodes))
!@        sump = DOT_G(ipnodes,amat_p(1,kcdiam),1,ythom,1)
!@
!@        write(*,1212)delphid,(sumi+sump)/wdiam,delphid0
!@ 1212   format(' DelPhid = ',1pe10.3,' PhiD = ',1pe10.3,
!@     >   ' DelPhid0 = ',1pe10.3)
 
      end subroutine getdiam

      
      subroutine getflux(amat_i, amat_p, data_array, r12sqr, 
     1  gobser1, gobser2, orsq, r12, z12, gsqrt, kcflux)
      use vmec_main
      use vmec_params, only: signgs
      use realspace
      use vsvd
      use vspline
      use vparams, only: zero
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer kcflux
      real(rprec), dimension(isnodes,*) :: amat_i
      real(rprec), dimension(ipnodes,*) :: amat_p
      real(rprec), dimension(*) :: data_array
      real(rprec), dimension(nrzt) :: r12sqr,
     1  gobser1, gobser2, orsq, r12, z12, gsqrt     
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: js, l, iloop, iobsmax, iobs, index, indexflx, 
     1   isym, n1, lk
      real(rprec) :: t1, t2, tpisg, sign0, sumi, sump, delta
C-----------------------------------------------
 
!       IRESIDUE > 0, OTHERWISE FLUX_INIT NOT CALLED YET!
        kcflux = 0
        if (iresidue.le.0 .or. (nflxs.eq.0 .and. iequi.eq.0)) return

!
!       COMPUTES MATRIX ELEMENTS NEEDED TO RELATE OBSERVED FLUX
!       MEASUREMENTS TO CURRENT AND PRESSURE EXPANSION COEFFICIENTS
!       R12,Z12 ARE THE PLASMA R,Z COORDINATES AT THE HALF
!       RADIAL NODE POINTS
!


!
!       COMPUTE SYMMETRIZED PSI(R,Z)+PSI(R,-Z) FLUX "GREEN'S FUNCTION"
!
!
!       COMPUTE "GREEN'S FUNCTION" KERNEL ONLY FOR NEEDED OBSERVATION POINTS
!       (OR FOR FINAL OUTPUT IF IEQUI=1)
!
      tpisg = twopi*signgs                     !Positive volume integral
      do n1 = 1, nobser
         isym = needflx(n1)
         if (isym.eq.needit .or. iequi.eq.1) then
            call grnflx (r12sqr, r12, z12, gobser1, nrzt, n1)
            do l = 2,nrzt
              gobser2(l) = gobser1(l)*orsq(l)
              gobser1(l) = gobser1(l)*gsqrt(l)
            end do  
!
!       DO INTEGRAL OVER ANGLES (ALL INTEGRALS ARE FROM THETA=0,TWOPI)
!       IM = <G/R**2>, PM = <G(1/R**2/<R**-2> - 1)>
!
            do js = 2, ns
              sumi = zero
              sump = zero
              do lk = js ,nrzt, ns
                sumi = sumi + gobser2(lk)*wint(lk)
                sump = sump + gobser1(lk)*wint(lk)
              enddo
              im(js,n1) = tpisg*sumi
              pm(js,n1) = (-tpisg*sump) + im(js,n1)/rm2(js)
            end do

          else if( isym.ge.ISYMCOIL )then    !only for up-down symmetric plasmas
            do js = 2,ns
              im(js,n1) = im(js,isym)
              pm(js,n1) = pm(js,isym)
            enddo
          endif
        enddo    !n1 loop
!
!       COMPUTE SPLINE MATRIX ELEMENTS BY INTEGRATING OVER RADIUS
!
        do 2000 iloop = 0,iequi                !iequi = 0 normally, = 1 at end
          iobsmax = nflxs
          if( iloop.eq.1 )iobsmax = nobd + nobser
          do 1000 iobs = 1,iobsmax
            indexflx = indxflx(iobs)
            if( iloop.eq.1 )indexflx = iobs
            if( indexflx.le.0 )go to 1000
            do js = 2,ns
              pm(js,0) = zero
              im(js,0) = zero
            enddo

            do l = 1,4                !This could be halved by using symmetry
              index = iconnect(l,indexflx)
              if( index.ne.0 )then
                sign0 = 1.0
                if( index.lt.0 )then
                  sign0 = -sign0
                  index = -index
                endif
                do js = 2,ns
                  pm(js,0) = pm(js,0) + sign0*pm(js,index)
                  im(js,0) = im(js,0) + sign0*im(js,index)
                enddo
              endif
            enddo
          
            do js = 2,ns
              pm(js,0) = pm(js,0)*ochip(js)
              im(js,0) = im(js,0)*ovrm2(js)
            enddo

               if (iloop .eq. 0) then
                  kcflux = kcflux + 1
                  delta = plflux(indexflx)
 
                  call splinint (im, current, amat_i(1,kcflux), hstark, 
     1               u_ib, u1_ib, w_ib, w1_ib, nk_ib, isnodes, intder, 
     2               ns)
 
                  call splinint (pm, presint, amat_p(1,kcflux), hthom, 
     1               u_pb, u1_pb, w_pb, w1_pb, nk_pb, ipnodes, intder, 
     2               ns)
 
                  t1 = one/sigma_flux(indexflx)
                  data_array(kcflux) = t1*delta
                  amat_i(:,kcflux) = t1*amat_i(:,kcflux)
                  t2 = t1*dmu0*pthommax         !!*pfac moved to getthom
                  amat_p(:,kcflux) = t2*amat_p(:,kcflux)
 
 
               else            !Store plasma fluxes in PLFLUX for output
                  plflux(indexflx) = zero
                  do js = 2, ns
                     plflux(indexflx) = plflux(indexflx) + pm(js,0)*
     1                  presph(js) + im(js,0)*(current(js)*iotaf(js)-
     2                  current(js-1)*iotaf(js-1))
                  end do
               endif
 1000     continue
 2000   continue

      end subroutine getflux

      
      subroutine getmse(amat_i, amat_p, data_array, re, ro, lue, luo, 
     1   zue, zuo, phipog, kcstark)
      use vmec_main
      use vmec_params, only: signgs
      use realspace
      use vsvd
      use vspline
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer kcstark
      real(rprec), dimension(isnodes,*) :: amat_i
      real(rprec), dimension(ipnodes,*) :: amat_p
      real(rprec), dimension(*) :: data_array
      real(rprec), dimension(ns,nzeta,*) :: 
     1   re, ro, lue, luo, zue, zuo
      real(rprec), dimension(*) :: phipog
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: c1p5 = 1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: lt, i, js, ks, j, irnodes
      real(rprec), dimension(ns,ntheta2) :: luef, luof
      real(rprec) :: dlu, dzu, guu1, edgeiota, guu2, edgefactor, 
     1   wedge, t1
      real(rprec), save :: facedge 
C-----------------------------------------------
 
!
!       THIS SUBROUTINE COMPUTES THE LEAST-SQUARES AMATRIX AND DATA
!       ARRAYS FOR MATCHING TO THE MSE DATA AT EQUALLY-SPACED KNOTS
!       IN SQRT(PHI-FLUX) SPACE (ISNODES, INCLUDING S=0 AND S=1)
!
!       THE RANGE OF MSE DATA IS EXTENDED TO INCLUDE S=1 BY USING THE
!       CURRENT MATCHING CONDITION TO CONSTRAIN IOTA(S=1)
!
!       COMING INTO THIS ROUTINE, THE RSTARK, DATASTARK HAVE BEEN
!       PREVIOUSLY ORDERED SO RSTARK(1) < RSTARK(2) < ....
!
 
!
!       COMPUTE FULL MESH VALUES OF LAMBDA
!
      call lamfull (luef, luof, lue, luo)
!
!       COMPUTE OUTER EDGE PITCH (IOTA) TO MATCH TOTAL CURRENT
!       IOTA(EDGE) = MU0 * IPLASMA/ 2 * PI *< Guu/SQRTg > PHIP(s=1)
!       NEED THIS TO SPLINE IOTA OVER FULL S-RANGE ( 0 .le. S .le.1 )
!
      guu1 = dot_product(c1p5*wint(ns:nrzt:ns)*guu(ns:nrzt:ns),
     1   phipog(ns:nrzt:ns))
      guu2 = dot_product(cp5*wint(ns-1:nrzt-1:ns)*guu(ns-1:nrzt-1:ns),
     1   phipog(ns-1:nrzt-1:ns))
 
      if (iresidue.eq.0 .or. iotaf(ns).eq.zero) then
         facedge = one
      else if (mod(iter2 - iter1,ipedsvd) .eq. 0) then
         facedge = (guu1*iotas(ns) - guu2*iotas(ns1))/(iotaf(ns)*(guu1
     1       - guu2))
      endif
      edgefactor = facedge*(guu1 - guu2)*signgs*twopi
      edgeiota = currv/edgefactor
 
      irnodes = max(0,imse) + 1
      lt = 1                                     !Outer R edge
      dlu = luef(ns,lt) + luof(ns,lt)
      wedge = (zue(ns,1,lt)+zuo(ns,1,lt))/(dlu*router)
      rstark(irnodes) = router
      datastark(irnodes) = wedge*edgeiota        !Edge pitch
      sigma_stark(irnodes) = abs(sigma_current*wedge/edgefactor)
 
!
!       THROW AWAY POINTS OUTSIDE GRID
!       NOTE: IF ONLY OUTER POINT KEPT, THE CALL TO SORT IS UNNECESSARY
!
      rsort0(:irnodes) = rstark(:irnodes)
      call sort_data (rsort0,isortr,irnodes)
      kcstark = 0
      do i = 1,irnodes
        j = isortr(i)
        if( ((rsort0(i).gt.rinner) .and.
     >       (rsort0(i).le.router)) .or. (iequi.ne.0) )then
           kcstark = kcstark+1
           rsort(kcstark) = rsort0(i)                        !kcstark <= i
           starkcal(kcstark) = datastark(j)                 !sorted data array
           qcalc(kcstark) = one/sigma_stark(j)                !qcalc = sorted weight array
        endif
      enddo
 
!
!       COMPUTE IOTA(0) FROM SLOPE AT RSTARK=R00
!
c04-96        kcstark = kcstark+1
c04-96        rsort(kcstark) = r00                !Magnetic axis (s=0)
c04-96        qcalc(kcstark) = 1.0/scstark
c04-96
c04-96        if( imse.gt.0 )then
c04-96        slope0 = 1.0
c04-96        call splint(rstark,ystark0,y2stark0,
c04-96     >  imse,r00,dum,slope0,1)
c04-96        starkcal(kcstark) = r00*slope0*luef(1,1)/dkappa
c04-96        else
c04-96c       EXTEND BOUNDARY POINTS TO INCLUDE FULL RANGE IN THETA
c04-96        starkcal(kcstark) = ai(0)
c04-96        endif
 
!
!       FIND S,THETA INDICES FOR RSORT ARRAY
!
      call findphi (re, ro, rsort, delse1, delso1, rmid, indexs1, 
     1   indexu1, indexr, kcstark)
 
!
!       COMPUTE MATRIX ELEMENTS FOR IOTA SPLINE NODES CORRESPONDING
!       TO ORDERED RSORT ARRAY ( = RSORT S )
!
      if (kcstark .gt. nmse) stop 'kcstark>nmse'
      call getspline (amat_i, sknots, hstark, delse1, hs, indexs1, 
     1   isorts, kcstark, isnodes)
 
!
!       MATCH TO STARK MSE DATA ACCORDING TO THE FORMULA:
!
!       Bz/Btor = IOTA*Zu/[ R*(1+LAMu) ]
!
!       NOTE: QCALC, DATA = STARKCAL CORRESPOND TO RSORT_S(I)
!       WITH INDEX KS = ISORTS(I) (INDEXED ON RSORT BEFORE IT WAS SORTED)
!       SAME IS TRUE FOR DELSE,O1, INDEXS1, INDEXU1
!
 
      islope = 0
      do i = 1, kcstark
c                     !Index BEFORE sorting on sknots (indexed on rsort)
         ks = isorts(i)
         js = indexs1(ks)
         lt = indexu1(ks)
!
!       COMPUTE WEIGHT(J) = Zu / (R * [1 + LAMu]), WHICH IS THE FACTOR
!       RELATING MSE PITCH = WEIGHT(J) * IOTA(J) TO ROTATIONAL TRANSFORM.
!
!       ON ENTRY INTO THIS LOOP,
!       QCALC = 1/SIGMA_STARK
!
         dlu = (one - delse1(ks))*luef(js,lt) + delse1(ks)*luef(js+1,lt
     1      ) + (one - delso1(ks))*sqrts(js)*luof(js,lt) + delso1(ks)*
     2      sqrts(js+1)*luof(js+1,lt)
         dzu = (one - delse1(ks))*zue(js,1,lt) + delse1(ks)*zue(js+1,1,
     1      lt) + (one - delso1(ks))*sqrts(js)*zuo(js,1,lt) + delso1(ks
     2      )*sqrts(js+1)*zuo(js+1,1,lt)
         stark_weight(ks) = dzu/(rsort(ks)*dlu)
 
         if (rsort(ks) .eq. router) icurrout = i
 
c04-96        if( rsort(ks).eq.r00 )then                        !IOTA(0)
c04-96          islope = i
c04-96          stark_weight(ks) = abs(wedge)                  !Need in
c04-96          starkcal(ks) = weight(ks)*starkcal(ks)        !Need for
c04-96          data_array(i) = starkcal(ks) * qcalc(ks)
c04-96          amat_i(1,i) = amat_i(1,i) * stark_weight(ks) * qcalc(ks)
c04-96         else
         data_array(i) = starkcal(ks)*qcalc(ks)
         t1 = stark_weight(ks)*qcalc(ks)
         amat_i(:isnodes,i) = t1*amat_i(:isnodes,i)
c04-96        endif
         if (i.eq.icurrout .and. qcalc(ks).eq.zero) stop 'CURR ERR'
      end do
 
      imse2 = kcstark
      amat_p(:,:kcstark) = zero
 
      end subroutine getmse

      
      subroutine gettflux
      use vsvd
      implicit none
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec) :: p01=1.e-2_dp, zero = 0.0_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec) :: tfac, phifac0
C-----------------------------------------------
 
!
!       UPDATE PHIEDGE SCALE FACTOR BY THE FORMULA
!       FDOT/F = -2*(1 - apres/aminor), WHERE F = PHISCALE
!
      tfac = p01*rsfac
 
      if (imatch_phiedge .eq. 0) then
         if (aminor .eq. zero) then
            phifac0 = phifac
         else
           phifac0 = phifac*(apres/aminor)
         end if  
         gphifac = tfac*(phifac0 - phifac)
 
      else if (imatch_phiedge .eq. 2) then
         call getlim
         gphifac = tfac*phifac*gphifac
      endif
 
      end subroutine gettflux

      
      subroutine getthom(amat_i, amat_p, data_array, re, ro, kcthom)
      use vmec_main
      use vsvd
      use vspline
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer kcthom
      real(rprec), dimension(isnodes,*) :: amat_i
      real(rprec), dimension(ipnodes,*) :: amat_p
      real(rprec), dimension(*) :: data_array
      real(rprec), dimension(ns,nzeta,*) :: re, ro
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, ks
      integer, dimension(itse + 2) :: isortp
      real(rprec), dimension(itse) :: datalsq_p
      real(rprec), dimension(itse + 2) :: rgrid
      real(rprec):: sig_fac, t1
      real(rprec) :: rmat1u(itse)
      logical :: l1v(itse)
C-----------------------------------------------
 
!
!       THIS SUBROUTINE COMPUTES THE LEAST-SQUARES AMATRIX AND DATA
!       ARRAYS FOR MATCHING TO THE PRESSURE DATA AT EQUALLY-SPACED KNOTS
!       IN SQRT(PHI-FLUX) SPACE (IPNODES, INCLUDING S=0 AND S=1)
!
!       COMING INTO THIS ROUTINE, THE RTHOM, DATATHOM HAVE BEEN
!       PREVIOUSLY ORDERED SO RTHOM(1) < RTHOM(2) < ....
!
 
!       IF (LPOFR), user has input P(R,Z)
!       If (.NOT.LPOFR),  then user has input P(s), not P(R)
!
 
      if (.not.lpofr) then                       !p(R) or p(s) ?
 
!       CONSTRUCT RTHOM BASED ON P(s)
         rthom(:itse) = sthom(:itse)
 
         call pofs (re, ro, ns, rthom, itse)
         rthommax = rthom(itse)
         rthommin = rthom(1)
 
      endif
 
!
!       IF NO PRESSURE DATA, MAKE SURE CHISQ-THOM <= CHI_TARGET
!
      sig_fac = one
 
!
!       CONSTRUCT EVERYTHING BASED ON P(R)
!       FOR POINTS OUTSIDE GRID, SET R = either rmin,rmax
!
      kcthom = 0
      if (itse .gt. 0) then
         l1v(:itse) = .false.
         datalsq_p(:itse) = datathom(:itse)*pfac   !sorted data array
         sigma_thom(:itse) = sigma_thom(:itse)/sig_fac
         pcalc(:itse) = 1.0/sigma_thom(:itse)      !pcalc = sorted sigma array
         rmat1u(:itse) = rthom(:itse)
         where (rmat1u(:itse) .lt. rinner) 
            rmat1u(:itse) = rinner
         elsewhere
            l1v(:itse) = rmat1u(:itse) .gt. router
         end where
         where (l1v(:itse)) rmat1u(:itse) = router
         rgrid(kcthom+1:itse+kcthom) = rmat1u(:itse)
         kcthom = itse + kcthom
      endif
 
 
!
!       FIND S,THETA INDICES FOR GRIDDED R-THOM ARRAY (RGRID)
!
      call findphi (re, ro, rgrid, delse2, delso2, rmid, indexs2, 
     1   indexu2, indexr, kcthom)
 
!
!       COMPUTE MATRIX ELEMENTS FOR PRESSURE SPLINE NODES CORRESPONDING
!       TO ORDERED RGRID ARRAY
!
      call getspline (amat_p, pknots, hthom, delse2, hs, indexs2,
     1   isortp, kcthom, ipnodes)
 
!
!       MATCH PRESSURE SPLINE KNOTS TO THOMSON SCATTERING DATA
!
!       ON ENTRY INTO THIS LOOP, PCALC = 1/SIGMA_THOM
!
 
      do i = 1, kcthom
         ks = isortp(i)         !Index BEFORE sorting on pknots (indexed
 
         data_array(i) = datalsq_p(ks)*pcalc(ks)
         t1 = pthommax*pcalc(ks)
         amat_p(:ipnodes,i) = t1*amat_p(:ipnodes,i)
      end do
      if (.not.lpofr) rthompeak = rgrid(isortp(1))
 
      itse2 = kcthom
      amat_i(:,:itse2) = zero
 
      end subroutine getthom

      
      subroutine grnbfld(xsqr, xs, zs, br, bz, idim, nobs1, nobs2)
      use vsvd
      use vparams, only: one
      implicit none
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: four=4.0_dp, p5=0.5_dp
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer idim, nobs1, nobs2
      real(rprec), dimension(idim) :: xsqr, xs, zs, br, bz
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer j, i1, i2
      real(rprec) :: xt, zt, xtsqr, xt4, oxt, zrp, zrm, xvv, 
     1 fxu, sqrxu, qqp, qqm, delqp, delqm, yeqp, yeqm, brp, brm
C-----------------------------------------------
!
!       COMPUTE BR = (ZT-ZS)/RT/SQRT(4*RT*RS) * F1(k)
!               BZ = 1/SQRT(4*RT*RS)*[RS/RT * F2(k) - F1(k)]
!       WHERE   F1 = k/2pi[ (E(k) - K(k)) + q1(k)*E(k) ]
!               F2 = k/(2pi)  [ q1(k)*E(k) ]
!               q1 = .5*k**2/(1. - k**2)   [Most singular piece near k=1]
!               k**2 = 4*RT*RS/[(RT+RS)**2 + (ZT-ZS)**2]
!
      xt = rbcoil(nobs1,nobs2)
      zt = zbcoil(nobs1,nobs2)
      xtsqr = p5/rbcoilsqr(nobs1,nobs2)            !1/2 from symmetrizing
      xt4 = four*xt
      oxt = one/xt

      do j = 2,idim
        zrp = zt - zs(j)
        zrm = zt + zs(j)
        xvv =(xt + xs(j))**2
        fxu = xs(j)*xt4
        sqrxu = xtsqr/xsqr(j)
        qqp = fxu/(xvv + zrp*zrp)
        qqm = fxu/(xvv + zrm*zrm)
!
!       WHICH INDEX LIES BELOW ?
!
        i1 = int(qqp*odqk2) + 1
        i2 = int(qqm*odqk2) + 1
!
!       LINEAR INTERPOLATION
!
        delqp = qqp - qsq(i1)
        delqm = qqm - qsq(i2)
        yeqp = (yeq(i1) + delqp*dyeq(i1))/(one - qqp)
        yeqm = (yeq(i2) + delqm*dyeq(i2))/(one - qqm)
        brp = yek(i1) + delqp*dyek(i1)
        brm = yek(i2) + delqm*dyek(i2)
        br(j) = sqrxu*oxt*(zrp*(brp+yeqp) + zrm*(brm+yeqm))
        bz(j) = sqrxu*((xs(j)*oxt-1.0)*(yeqp + yeqm) - (brp + brm))
      enddo

      end subroutine grnbfld

      
      subroutine grnflx(xsqr, xs, zs, ansp, idim, nobs)
      use vsvd
      implicit none
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: p5 =0.5_dp, four=4.0_dp
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer idim, nobs
      real(rprec), dimension(idim) :: xsqr, xs, zs, ansp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i1, i2, j
      real(rprec) :: xt, zt, xtsqr, xt4, 
     4   zrp, zrm, xvv, fxu, sqrxu, qqp, qqm
C-----------------------------------------------
!
!       EVALUATES "GREEN'S" FUNCTION FOR POLOIDAL FLUX AT INTERIOR
!       POINTS XS,ZS AND OBSERVATION POINT XT,ZT (ANSP)
!       (RECALL THETA INTEGRATION ONLY FROM ZERO TO PI, SO NEED
!        TO REFLECT ZS TO -ZS, AT LEAST IN UP-DOWN SYMMETRIC CASE)
!
      xt = xobser(nobs)
      zt = zobser(nobs)
      xtsqr = p5*xobsqr(nobs)        !1/2 factor from averaging up,down
      xt4 = four*xt
      do j = 2,idim
        zrp = zt - zs(j)
        zrm = zt + zs(j)
        xvv =(xt + xs(j))**2
        fxu = xs(j)*xt4
        sqrxu = xsqr(j)*xtsqr
        qqp = fxu/(xvv + zrp*zrp)                !k**2 for zplasma > 0
        qqm = fxu/(xvv + zrm*zrm)                !k**2 for zplasma < 0
!
!       WHICH INDEX LIES BELOW ?
!
        i1 = int(qqp*odqk2) + 1
        i2 = int(qqm*odqk2) + 1
!
!       LINEAR INTERPOLATION
!
        ansp(j)  =  sqrxu *( ( yf(i1)+(qqp-qsq(i1))*dyf(i1) )
     >   + ( yf(i2)+(qqm-qsq(i2))*dyf(i2) ) )
      enddo
 
      end subroutine grnflx

      
      subroutine pofs(re, ro, ns, rthom, itse)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ns, itse
      real(rprec), dimension(ns) :: re, ro
      real(rprec), dimension(*) :: rthom
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: j, js, i
      real(rprec), dimension(ns) :: s
      real(rprec) :: sqrjs, rt2, ds
C-----------------------------------------------
!
!       Interpolate Rmid = Re + Ro to get R(s)
!
!       THIS ROUTINE INTERPOLATES THE RTHOM "S-SPACE" ARRAY
!       ONTO THE INSTANTANEOUS RMID ARRAY
!       ON INPUT, RTHOM IS THE S-VALUE ARRAY FOR THOMPSON DATA
!       ON OUTPUT,RTHOM IS THE CORRESPONDING (THETA=0) RMID ARRAY
!
      do j = 1, ns
         s(j) = real(j - 1,rprec)/(ns - 1)
      end do
 
      js = 1
      do i = 1, itse
         rt2 = rthom(i)
  100    continue
         if (rt2.ge.s(js) .and. rt2.le.s(js+1)) then
            ds = (rt2 - s(js))/(s(js+1)-s(js))
            sqrjs = sqrt(rt2)
            rthom(i) = re(js) + (re(js+1)-re(js))*ds + sqrjs*
     1         (ro(js)+(ro(js+1)-ro(js))*ds)
         else
            js = js + 1
            if (js < ns) go to 100
         endif
      end do
 
      end subroutine pofs

      
      subroutine radfor(pfac0)
      use vmec_main
      use vmec_params, only: signgs
      use vforces, only : r12=>armn_o, gsqrt=>azmn_o, gor=>clmn_o
      use realspace
      use vsvd
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec) pfac0
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec) :: p05 = 0.05, p5 = 0.5_dp, c1p5 = 1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: js
      real(rprec), dimension(ns) :: vpres
      real(rprec) :: delpres, pedge, t1
      real(rprec), external :: dot_g
C-----------------------------------------------
 
 
!
!       COMPUTE VPRES, NEEDED FOR F00 PRESSURE BALANCE
!
      gor(2:nrzt) = gsqrt(2:nrzt) / r12(2:nrzt)
      do js = 2, ns
         vpres(js) =signgs*DOT_G(nznt,gor(js),ns,wint(js),ns)
      end do
 
      pedge = c1p5*pres(ns) - p5*pres(ns1)
      pressum0 = dot_product(wint(ns:nrzt:ns)*zu0(ns:nrzt:ns),
     1   r1(ns:nrzt:ns,0)+r1(ns:nrzt:ns,1))
      pressum0 = signgs*pedge*pressum0
      pressum0 = pressum0 + hs*dot_product(vpres(2:ns),pres(2:ns))
 
      if (pressum0 .eq. zero) pressum0 = one
 
      pfac0 = pfac
      if (iresidue .ge. 3) return                    !!AXIS MOVED BY FSQR IN RESIDUE
!
!       COMPUTE AVERAGE FORCE BALANCE CONSTRAINT FOR FIXING R(0)
!       (INTEGRAL OF RADIAL FORCE BALANCE,M=0,N=0, OVER S)
!
 
      if (1.e6_dp*fsq .le. one) then
         delpres = 0.
         delpres = -fsqsum0/pressum0
         t1 = abs(delpres)
         if (t1 .gt. p05) delpres = p05*delpres/t1   !!Wait til close
         pfac0 = pfac*(one + delpres)
      endif
 
      end subroutine radfor

      
      subroutine readrecon
      use vmec_main
      use vsvd
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      character*(50), dimension(2), save :: raxis_message
      character*(50), dimension(4), save :: phiedge_message
C-----------------------------------------------
      data raxis_message/'Magnetic axis position fixed', 
     1   'Magnetic axis position optimized'/
      data phiedge_message/
     1   'Phiedge determined to match pressure minor radius', 
     2   'Phiedge matched to input value', 
     3   'Phiedge determined by limiter position', 
     4   'Phiedge determined by Ip'/
 
      iphidiam = 0
      if (imse > nmse) stop 'IMSE>NMSE'
      if (itse > ntse) stop 'ITSE>NTSE'
      if ((imse>0 .or. nflxs>0 .or. nbfldn>0) .and. itse.ge.0) then
         iresidue = 0
      else
         lrecon = .false.
      end if   
 
      if (.not.lrecon) return                   !No reconstruction matching
      ncurr = 0                                  !Just to be safe
      if (sigma_current .ge. cbig) stop 'SIGMA_CURRENT missing'
      if (sigma_delphid .ge. cbig) print *, ' SIGMA_DELPHID missing'
      if (sigma_current < zero) sigma_current = abs(sigma_current*
     1   curtor)
      if (sigma_delphid < zero) sigma_delphid = abs(sigma_delphid*
     1   phidiam)
      write (nthreed, 150)
  150 format(/' DATA MATCHING PARAMETERS: ',/,1x,35('-'))
      write (nthreed, 155) imse, itse, nflxs, nobser, nobd, nbfldn, 
     1   nbcoilsn, sigma_current, 1.e3_dp*sigma_delphid, tensp, tensi,
     2   tensi2, fpolyi, mseangle_offset, presfac, pres_offset, lpofr
      write (nthreed, 152) mseangle_offsetm
  152 format('mse-angleM offset',/,f13.3)
  155 format('   imse       itse      nflxs     nobser       nobd',
     1'     nbfldn    nbcoils  sigma_current(A)   sigma_delphid(mWb)',/,
     2   i7,6i11,3x,1pe15.3,4x,0pf17.3,/,
     3   '    tension(p)   tension(i)  tension2(i)  fpolyi  ',
     4   'mse-angle offset  pres scale factor pressure offset  lpofr',/,
     5   3f13.3,f9.3,f18.3,f19.3,f16.3,6x,l1)
      write (nthreed, 200)
  200 format(/,' LEGEND',/,1x,6('-'))
 
      if (curtor < cbig) then
         write (nthreed, 210) 1.e-6_dp*curtor
      else
         write (nthreed, *) 'Need toroidal plasma current'
         stop 15
      endif
  210 format(' Matching to toroidal current = ',f10.3,' [MA]')
      sigma_current = dmu0*sigma_current
      if (nflxs > 0) then
         write (nthreed, *) 'Fitting ', nflxs, 
     1      ' external flux loop measurements'
      else
         write (nthreed, *) 
     1      'Not fitting external flux loop measurements.'
      endif
      if (phidiam<cbig .and. sigma_delphid<cbig) then
         iphidiam = 1
         write (nthreed, 220) 1.e3_dp*phidiam
      else
         write (nthreed, *) 'No fit to diamagnetic flux'
      endif
  220 format(' Fitting diamagnetic flux     = ',f10.3,' [mWb]')
      write (nthreed, *) raxis_message(iopt_raxis+1)
      write (nthreed, *) phiedge_message(imatch_phiedge+1)

      end subroutine readrecon

      
      subroutine sgemvmm(amat_i, amat_p, amatsq, b, bsq, wten,
     1   mdata, niota, npres, nots)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer mdata, niota, npres, nots
      real(rprec), dimension(niota,mdata) :: amat_i
      real(rprec), dimension(npres,mdata) :: amat_p
      real(rprec), dimension(nots,nots) :: amatsq
      real(rprec), dimension(mdata) :: b
      real(rprec), dimension(nots) :: bsq, wten
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0.0_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, j, ioff, joff
C-----------------------------------------------

      amatsq = zero
      bsq    = zero
 
!
!       INITIALIZE IOTA, PRESSURE DIAGONAL ELEMENTS (ALREADY 'SQUARED')
!
!
!       COMPUTE LOWER DIAGONAL HALF OF SQUARE OF MATRIX
!       A(trans)*A and A(trans)*B
!       BY ADDING CONTRIBUTIONS FROM EXTERNAL MAGNETICS SIGNALS
!
 
!
!       FIRST UPPER LEFT NIOTA X NIOTA BLOCK
!
      do i = 1, niota
         bsq(i) = bsq(i) + sum(b*amat_i(i,:))
         do j = 1, i
            amatsq(i,j) = amatsq(i,j) + sum(amat_i(i,:)*amat_i(j,:))
         end do
      end do
 
!
!       LOWER NPRES X NIOTA BLOCK, NPRES X NPRES BLOCK
!
      do i = 1, npres
         ioff = i + niota
         bsq(ioff) = bsq(ioff) + sum(b*amat_p(i,:))
         do j = 1, niota
            amatsq(ioff,j) = amatsq(ioff,j) + 
     1                       sum(amat_p(i,:)*amat_i(j,:))
         end do
         do j = 1, i
            joff = j + niota
            amatsq(ioff,joff) = amatsq(ioff,joff) +
     1                          sum(amat_p(i,:)*amat_p(j,:))
         end do
      end do
 
      do i = 1, nots
         wten(i) = amatsq(i,i)
         amatsq(1:i-1,i) = amatsq(i,1:i-1)
      end do
 
 
      end subroutine sgemvmm
      

      subroutine smoothdata(nwout)
      use vmec_main
      use vsvd
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: nwout
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: ndata_elems = 11
c-----------------------------------------------
c   L o c a l   V a r i a b l e s
c-----------------------------------------------
      integer :: npts, i, ndata1
      real(rprec), dimension(2*ns-1) :: hmid,ymid,y2mid,
     1  wmid,tenmid,requal
      real(rprec), dimension(:), pointer :: data

c-----------------------------------------------
!
!       spline output data onto equally-spaced r-mesh for plotting
!
      npts = 2*ns - 1
      if (npts .le. 1) return
      
      do i = 1, npts
         if (i .le. ns) curint(i) = twopi/dmu0*buco(ns + 1 - i)
         if (i .gt. ns) curint(i) = twopi/dmu0*buco(i - ns + 1)
      end do
 
      do i = 1, npts
         wmid(i) = 1.0
         tenmid(i) = 0.1
         requal(i) = rmid(1) + ((i - 1)*(rmid(npts)-rmid(1)))/
     1      (npts - 1)
      end do
 
      do ndata1 = 1, ndata_elems
         select case (ndata1) 
         case (1) 
            data => datamse
         case (2) 
            data => qmid
         case (3) 
            data => shear
         case (4) 
            data => presmid
         case (5) 
            data => alfa
         case (6) 
            data => curmid
         case (7) 
            data => curint
         case (8) 
            data => psimid
         case (9) 
            data => ageo
         case (10) 
            data => volpsi
         case (11) 
            data => phimid
         end select
         call setspline (rmid, wmid, data, hmid, ymid, y2mid, tenmid, 
     1      tenmid(1), npts, natur)
         call splint (rmid, ymid, y2mid, npts, requal, data(1), 
     1      zero, npts)
      end do
!
!     write out splined data
!
      write (nwout, 703) (datamse(i),requal(i),qmid(i),shear(i),presmid
     1   (i),alfa(i),curmid(i),i=1,npts)
      write (nwout, 703) (rsort(i),atan(datastark(isortr(i)))/dcon,abs(
     1   qmeas(i)),i=1,imse2 - 1)
      write (nwout, 703) (rthom(i),datathom(i),i=1,itse)
  703 format(5e20.13)
 
      if (lmac) then
        write (nmac, 705)
  705 FORMAT(//,' FOLLOWING DATA EQUALLY SPACED IN R-MIDPLANE'//,
     1   '        RMID       J-PHI       SHEAR        QMID',
     2   '   MSE-PITCH     PRESMID         PSI        AMID',
     3   '      VOLUME         PHI',/,
     4   '         [M]    [A/M**2]                        ',
     5   '       [Deg]        [Pa]        [Wb]         [M]',
     6   '      [M**3]        [Wb]',/)
        write (nmac, 707) (requal(i),curmid(i),shear(i),qmid(i),
     1    datamse(i),presmid(i),psimid(i),ageo(i),volpsi(i),
     2    phimid(i),i=1,npts)
        write (nmac, 709) phimid(2*ns-1), psimid(2*ns-1)
      end if
  707 format(1p10e12.3)
  709 format(/,' phi-edge =',t16,1pe12.3,t40,'psi-edge =',t56,1pe12.3)
       
      end subroutine smoothdata
EOF
cat > splines.f << "EOF"
      subroutine add_tension(amat, wten, hx, tens, tensv, fpoly, n, nb, 
     1   ioff, nmat)
      use vspline
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer n, nb, ioff, nmat
      real(rprec) tens, tensv, fpoly
      real(rprec), dimension(nmat,*) :: amat
      real(rprec), dimension(nmat) :: wten
      real(rprec), dimension(n) :: hx
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: k, i, koff
      real(rprec), dimension(n) :: wten0, tenshx, work
      real(rprec) :: delta_x, tension
C-----------------------------------------------
 
!
!       THIS SUBROUTINE ADDS THE TENSION TERM WTEN*(DEL(K) - DEL(K-1))
!       TO THE REST OF THE CHI-SQ AMAT, WHERE DEL(K) = (G(K+1)-G(K))/H(K)
!       AND G(K) IS THE SECOND DERIVATIVE AT THE K-TH KNOT
!
!       IOFF ALLOWS FOR ADDING TENSION SEPARATELY TO IOTA (IOFF=0)
!       AND PRESSURE (IOFF=NIOTA) SPLINES
!
!       tens:   spline tension (optionally, at the 1st pt only;see note)
!       tensv:  vbl spline tension for n-1th point (see note)
!       fpoly:  vbl spline tension form factor (note: if tens1<>tens2
!               then tension(i-th point) = tens+(tensv-tens)*(i/n-1))**fpoly)
 
!
!       BOUNDS CHECKING
!
      if (n + ioff > nmat) stop '(n+ioff>nmat)'
      if (fpoly < 0.) stop '(fpoly<0)'
      if (n < 1) stop '(n < 1)'
 
!
!       COMPUTE TENSION COEFFICIENTS
!
      
      delta_x = sum(hx(:n-1))
      tension = 0.5*(delta_x/(n))**3
      if (fpoly.eq.0. .or. tens.eq.tensv) then
         tenshx(:n-1) = tens
      else
         do i = 1, n - 1
            tenshx(i) = tens + (tensv - tens)*(real(i - 1,rprec)/
     1         (n - 1))**fpoly
         end do
      endif
 
      do i = 1,n-1
        tenshx(i) = tension * tenshx(i) / hx(i)
        work(i) = hx(i)*(wten(i+ioff) + wten(i+ioff+1))
      enddo
      do i = 2,n-1
        wten0(i) = 0.5 * ( work(i) + work(i-1) )/(hx(i) + hx(i-1))
      enddo
      wten0(1) = wten0(2)
      wten0(n) = wten(n+ioff)
!
!       COMPUTE, FOR K = 1,N, B(K,L)*JACOBIAN(L,I) = W(K,I),
!       WHERE JACOBIAN = D[G]/D[F] and B is TRIDIAGONAL
!       SEE EQN(27) IN PHYS.PLASMAS 1, p 2277.
!
      do k = 1, n
         koff = k + ioff
         work(:n-1) = 0
!       SET UP COEFFICIENTS IN [G(K+1)-G(K)]/h(k) - [G(K)-G(K-1)]/h(k-1)
         if (k .eq. 1) then
            work(2) = tenshx(1)*wten0(2)
            work(1) = -work(2)
         else if (k .eq. n) then
            work(n-1) = tenshx(n-1)*wten0(n-1)
         else
            work(k-1) = tenshx(k-1)*wten0(k-1)
            work(k) = -(tenshx(k)+tenshx(k-1))*wten0(k)
            work(k+1) = tenshx(k)*wten0(k+1)
         endif
         if (nb .eq. natur) work(1) = 0
         work(n) = 0
!
!       COMPUTE work(j) = work(i)*Jacobian(i,j) and add to amat(k,j)
!
         call jacprod (work, hx, n, nb)
         amat(koff,1+ioff:n+ioff) = amat(koff,1+ioff:n+ioff) + work(:n)
      end do
 
      end subroutine add_tension
      
      subroutine getspline(amat, splnots, hk, delse, hs, indexs, isort, 
     1   ndata0, nots)
      use vsvd0
      use vparams, only: zero
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ndata0, nots
      real(rprec) hs
      integer, dimension(ndata0) :: indexs, isort
      real(rprec), dimension(nots,ndata0) :: amat
      real(rprec), dimension(nots) :: splnots, hk
      real(rprec), dimension(ndata0) :: delse
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, dimension(nots) :: nk
      integer :: i, js, ia, k, k1, ib, j, nb
      real(rprec), dimension(ndata0) :: w, w1, u, u1, snodes
      real(rprec), dimension(nots,ndata0) :: bmat
C-----------------------------------------------
 
!
!       ON EXIT, AMAT = AMAT + BMAT, WHERE BMAT IS THE 2nd
!       DERIVATIVE COEFFICIENT MATRIX ARRAY, MULTIPLIED BY JACOBIAN
!       AND AMAT (ON RHS) WAS FUNCTION COEFFICIENT MATRIX ARRAY
!
 
!
!       SORT KNOT POSITIONS IN ASCENDING ORDER IN S-SPACE
!       USE SQRT(S) KNOT POSITIONS FOR IMPROVED RESOLUTION
!       NOTE: SNODES(I) IS THE VALUE OF SQRT(S) CORRESPONDING TO
!       THE MESH VALUES CORRESPONDING TO DELSE, INDEXS (COMPUTED OUTSIDE
!       THIS PROGRAM)
!
 
      do i = 1, ndata0
         js = indexs(i)
         snodes(i) = sqrt(hs*((js - 1) + delse(i)))
!        if (snodes(i) .le. zero) snodes(i) = epstan
      end do

!     Avoid roundoff error in SPLININT
      if( snodes(ndata0) .gt. splnots(nots) )
     1    snodes(ndata0) = splnots(nots)
 
      call sort_data (snodes, isort, ndata0)
 
!
!       COMPUTE MATRIX COEFFICIENTS RELATING SPLINE AT SPLNOTS
!       TO REAL-SPACE FUNCTION AT SORTED MESH POINTS RMESH
!
      amat(:nots,:ndata0) = zero
      bmat(:nots,:ndata0) = zero
 
!
!       SETUP SPLINE PARAMETERS AT EACH TIME STEP, SINCE SNODES
!       MAY BE CHANGING DYNAMICALLY IN S-SPACE
!
      call setup_int(splnots,snodes,hk,w,w1,u,u1,nk,nots,ndata0)
 
      ia = 1
      do k = 1, nots - 1
         if (nk(k) .gt. 0) then
            k1 = k + 1
            ib = ia + nk(k) - 1
            amat(k,ia:ib) = amat(k,ia:ib) + w(ia:ib)
            bmat(k,ia:ib) = bmat(k,ia:ib) + u(ia:ib)
            amat(k1,ia:ib) = amat(k1,ia:ib) + w1(ia:ib)
            bmat(k1,ia:ib) = bmat(k1,ia:ib) + u1(ia:ib)
            ia = ib + 1
         endif
      end do
 
      if (ib .ne. ndata0) stop 'ib!=ndat'
      nb = ideriv
 
      do j = 1, ndata0
         bmat(nots,j) = 0.
         call jacprod (bmat(1,j), hk, nots, nb)
         amat(:nots,j) = amat(:nots,j) + bmat(:nots,j)
      end do
 
      end subroutine getspline
      
      subroutine gety2(y, y2, h, nots, nb)
      use vspline
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nots, nb
      real(rprec), dimension(*) :: y, y2, h
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer jmax
      real(rprec), dimension(nots) :: aspline, bspline, dspline
C-----------------------------------------------
 
      jspmin(1) = 2
      if (nb .eq. ideriv) jspmin(1) = 1
      jmax = nots - 1
      aspline(1) = h(1)
      dspline(1) = 2.0*h(1)
      y2(1) = 0.
      y2(nots) = 0.
      if (nb .eq. ideriv) y2(1) = 6.0*(y(2)-y(1))/h(1)
      aspline(2:jmax) = h(2:jmax)
      bspline(2:jmax) = h(:jmax-1)
      dspline(2:jmax) = 2.0*(h(2:jmax)+h(:jmax-1))
      y2(2:jmax) = 6.0*((y(3:jmax+1)-y(2:jmax))/h(2:jmax)
     1   -(y(2:jmax)-y(:jmax-1))/h(:jmax-1))
 
      call tridslv(aspline,dspline,bspline,y2,jspmin,jmax,0,nots,1)
 
      end subroutine gety2
      
      subroutine initspline(amat, splnot, h, weight, nots)
      use vspline
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nots
      real(rprec), dimension(nots,nots) :: amat
      real(rprec), dimension(*) :: splnot, h, weight
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer i
      real(rprec) :: eps
C-----------------------------------------------
 
      if (nots .lt. 3) stop 'nots<3'
      eps = 1.0/(splnot(nots)-splnot(1))
 
      amat = 0.
      do i = 1, nots
         amat(i,i) = weight(i)
      end do
 
      do i = 1, nots - 1
         h(i) = splnot(i+1) - splnot(i)
         if (eps*h(i) .le. 1.e-8_dp) stop 'h(i)<1.e-8'
      end do
 
      end subroutine initspline
      
      subroutine jacprod(c, h, nots, nb)
      use vspline
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nots, nb, jmax
      real(rprec), dimension(*) :: c, h
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec), dimension(nots) :: 
     1   aspline, bspline, dspline, dum1
C-----------------------------------------------
 
!
!       THIS ROUTINE COMPUTES THE INNER PRODUCT COUT(I) = CIN(J)*JACOBIAN(J,I)
!       WHERE JACOBIAN(J,I) = D[G(J)]/D[F(I)]
!       HERE, G(J) ARE SECOND-DERIVATIVE KNOTS, F(I) FUNCTION KNOTS
!
!       COMPUTE COEFFICIENT ARRAY ELEMENTS A*X(I+1) + D*X(I) + B*X(I-1)
!       (TO BE SAFE, RECOMPUTE EACH TIME, SINCE IOTA, P SPLINES MAY
!        DIFFER FROM CALL TO CALL)
!
      aspline(1) = h(1)
      dspline(1) = 2.0*h(1)
      aspline(2:nots-1) = h(2:nots-1)
      bspline(2:nots-1) = h(:nots-2)
      dspline(2:nots-1) = 2.0*(h(2:nots-1)+h(:nots-2))
 
      jspmin(1) = 2
      if (nb .eq. ideriv) jspmin(1) = 1
      jmax = nots - 1
      call tridslv(aspline,dspline,bspline,c,jspmin,jmax,0,nots,1)
      dum1(1) = 6.0*(c(2)-c(1))/h(1)
      dum1(2:nots) = 6.0*(c(:nots-1)-c(2:nots))/h(:nots-1)
      c(2:nots-1) = dum1(2:nots-1) - dum1(3:nots)
      c(1) = dum1(1)
      c(nots) = dum1(nots)

      end subroutine jacprod
      
      subroutine set_dual(data, hi, yi, yi2, hp, yp, yp2, wten, alsq, 
     1   niota, npres, nots)
      use vspline
      use vparams, only: zero
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer niota, npres, nots
      real(rprec), dimension(nots) :: data
      real(rprec), dimension(niota) :: hi, yi, yi2
      real(rprec), dimension(npres) :: hp, yp, yp2
      real(rprec), dimension(nots) :: wten
      real(rprec), dimension(nots,nots) :: alsq
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer nb, ioff
C-----------------------------------------------
!
!       ADD TENSION TO DIAGONAL BLOCKS
!
      nb = ideriv
      ioff = 0
      call add_tension (alsq, wten, hi, tensi, tensi2, fpolyi, niota, nb
     1   , ioff, nots)
      call add_tension (alsq, wten, hp, tensp, zero, zero, npres, nb, 
     1   niota, nots)
 
!
!       FREEZE EDGE PRESSURE IF NO PRESSURE SPECIED
!       COMPUTE SOLUTION FOR SPLINES
!
      call solver (alsq, data, nots)
      yi(:niota) = data(:niota)
      yp(:npres) = data(1+niota:npres+niota)
 
!
!       COMPUTE SECOND DERIVATIVES
!
      call gety2 (yi, yi2, hi, niota, nb)
      call gety2 (yp, yp2, hp, npres, nb)
 
      end subroutine set_dual
      
      subroutine setspline(x,weight,y,h,yfit,y2,wten,tens,nots,nb)
      use vspline
      use vparams, only: zero
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nots, nb
      real(rprec) tens
      real(rprec), dimension(*) :: x, weight
      real(rprec), dimension(nots) :: y
      real(rprec), dimension(*) :: h
      real(rprec), dimension(nots) :: yfit
      real(rprec), dimension(*) :: y2, wten
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ioff, istat
      real(rprec) :: alsq(nots,nots)
C-----------------------------------------------
!
!       x:        independent coordinate array
!       y:        dependent y(x) array
!       yfit:     fitted values (under tension) to y array
!       h:        x(i+1) - x(i) array
!       y2:       y'' array used for splines
!       wten:        weight array for tension (changed on exit)
!       alsq:        matrix elements for least squares fit (from s-integrations)
!       nots:        number of independent coordinates (knots)
!       nb:        = NATUR, use natural boundary condition at left knot
!                  = IDERIV, use derivative (dy/dx =0) boundary condition at left knot
 
!
!       IT IS ASSUMED THAT X,Y,WTEN ARE ALL SORTED (ON X(I) < X(I+1))
!
 
!
!       INITIALIZE ALSQ TO ZERO, COMPUTE H ELEMENTS
!
      call initspline (alsq, x, h, weight, nots)
 
!
!       SET UP SPLINE MATRIX ASPLINE AND NON-DIMENSIONLIZE TENSION
!
      ioff = 0
      call add_tension (alsq, wten, h, tens, zero, zero, nots, nb, 
     1   ioff, nots)
 
!
!       SOLVE FOR COEFFICIENTS
!
      yfit(:nots) = y(:nots)
      call solver (alsq, yfit, nots)
!
!       OBTAIN Y'' COEFFICIENTS AND STORE IN Y2
!
      call gety2 (yfit, y2, h, nots, nb)
 
      end subroutine setspline
      
      subroutine setup_int(xknots,smesh,hx,w,w1,u,u1,nk,nots,nmesh)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nots, nmesh
      integer, dimension(nots) :: nk
      real(rprec), dimension(nots) :: xknots
      real(rprec), dimension(nmesh) :: smesh
      real(rprec), dimension(nots) :: hx
      real(rprec), dimension(nmesh) :: w, w1, u, u1
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: epstan = 1.d-10
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ksp1, k, i, k1
      real(rprec) :: smesh1, hk6
C-----------------------------------------------
!
!       FOR THE SPLINED FUNCTIONS CORRESPONDING TO XKNOTS (PRESSURE,
!       IOTA), THIS ROUTINE COMPUTES THE A AND B MATRIX ELEMENTS
!       (STORED IN W,U,W1,U1) WHICH ARE NEEDED TO EVALUATE THE FUNCTIONS
!       IN REAL-SPACE IN TERMS OF THEIR (VARIABLE) SPLINE KNOT VALUES.
!       THIS 'UNDOES' WHAT SPLINT ROUTINE DOES. LET Y(I) DENOTE THE
!       FUNCTION AT THE POINT SMESH(I) SUCH THAT
!
!                   XKNOTS(K) < SMESH(I) <= XKNOTS(K)
!
!       THEN,  Y(I) = W(I)*YK  + U(I)*GK  + W1(I)*YK1   + U1(I)*GK1
!
!       WHERE YK, GK ARE THE SPLINE AND 2ND DERIVATIVES AT KNOT K
!             YK1,GK1 ARE THE SAME AT KNOT K+1
!
      ksp1 = nots - 1
      smesh1 = smesh(1)
      if (smesh1 .le. xknots(1)) smesh(1) = xknots(1) + epstan
 
      nk = 0
 
      k = 1
      do i = 1, nmesh
  140    continue
         k1 = k + 1
!
!       XKNOTS = SQRT(HS*(JS-1)) DEFINED IN STARK,PRESSURE ROUTINE
!       (THIS CORRESPONDS TO APPROXIMATELY EQUAL SPACING ALONG MIDPLANE)
!
         if (smesh(i).gt.xknots(k) .and. smesh(i).le.xknots(k1)) then
            nk(k) = nk(k) + 1
            hk6 = hx(k)*hx(k)/6.0
            w1(i) = (smesh(i)-xknots(k))/hx(k)
            if (w1(i)<(-epstan) .or. w1(i)>1.0+epstan) stop 'w1(i)'
            w(i) = 1.0 - w1(i)
            u(i) = hk6*w(i)*(w(i)*w(i)-1.0)
            u1(i) = hk6*w1(i)*(w1(i)*w1(i)-1.0)
         else
            k = k + 1
            if (k .gt. ksp1) stop 'K>KSP1'
            go to 140
         endif
      end do
 
      smesh(1) = smesh1
 
      end subroutine setup_int
      
      subroutine sort_data (x, index_array, n)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
c-----------------------------------------------
      integer n
      integer, dimension(n) :: index_array
      real(rprec), dimension(n) :: x
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, j
      integer, dimension(1) :: isamax
      real(rprec), dimension(n) :: dumx
C-----------------------------------------------
!
!       RETURNS INDEX(I) ARRAY, SO THAT X(INDEX(I)) IS SORTED, I=1,N
!       RETURNS Xin(INDEX(I)) = Xout(I)
!
 
      do i = n, 1, -1
         isamax = maxloc(abs(x))
         j = isamax(1)
         dumx(i) = x(j)
         x(j) = 0.
         index_array(i) = j
      end do
 
      x = dumx
 
      end subroutine sort_data
      
      subroutine splinint(grn, cm, jacob, h, u, u1, w, w1, nk, nots, 
     1   ifunc, nmesh)
      use vparams
      use vsvd
      use vspline
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nots, ifunc, nmesh
      integer, dimension(*) :: nk
      real(rprec), dimension(nmesh) :: grn, cm
      real(rprec), dimension(nots) :: jacob, h
      real(rprec), dimension(nmesh) :: u, u1, w, w1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: j, k, nb, ia, ib, k1, ksp1, nmesh1
      real(rprec), dimension(nmesh) :: func
      real(rprec), dimension(nots) :: af, bs
C-----------------------------------------------
 
!
!       COMPUTES af,bs FACTORS IN
!       (ifunc=INTFUN)     Int[ GRN(X) * F(X) ds ]      or
!       (ifunc=INTDER)     Int[ GRN(X) * {d(CmF)/ds} ds ]
!                 =   af(k)*f(k) + bs(k)*f''(k)
!       WHERE f(k),f''(k) ARE SPLINE COEFFICIENTS OF F(X)
!
!       NOTE: FOR ifunc = INTDER, the OHS factor in CmF cancels
!             THE HS factor in the Integral.
!             FOR ifunc = INTFUN, GRN is assumed to be pre-multiplied
!             OUTSIDE this routine by HS factor
!       ALSO, COMPUTES af(k) + (SUM on i)bs(i)*J(i,k) = jacob(k),
!       WHERE J(i,k) = d[g(i)]/d[f(k)], g = f''
!
!       nk(k): Number of smesh-pts in k-th spline interval
!       xknots(k) < smesh <= xknots(k+1),   k = 1,nots-1
!
!       NOTE: The ifunc=INTDER case is done by integrating by parts,
!             so that the half-point integration (GRN at half mesh pts)
!             becomes a full-point integration in Cm*F.
!
      nb = ideriv           !Pressure, iota derivatives vanish at origin
      ksp1 = nots - 1
      nmesh1 = nmesh - 1
 
      if (ifunc .eq. intder) then
!
!       Integrate by parts (in s), func(1) and func(nmesh) are 'surface terms'
!
         func(1) = -cm(1)*grn(2)
         func(2:nmesh1) = cm(2:nmesh1)*(grn(2:nmesh1)-grn(3:nmesh1+1))
         func(nmesh) = cm(nmesh)*grn(nmesh)
      else
         func = grn
      endif
 
      af(:nots) = zero
      bs(:nots) = zero
 
      ia = 1
      do k = 1, ksp1
         if (nk(k) .ne. 0) then
            k1 = k + 1
            ib = ia + nk(k) - 1
            do j = ia,ib
              af(k)  = af(k)  + func(j)*w(j)
              bs(k)  = bs(k)  + func(j)*u(j)
              af(k1) = af(k1) + func(j)*w1(j)
              bs(k1) = bs(k1) + func(j)*u1(j)
            enddo
            ia = ib + 1
         endif
      end do
 
      if (ib .ne. nmesh) stop 'ib!=nmesh'
      if (nb .eq. natur) bs(1) = 0.           !Natural boundary conditions
      bs(nots) = 0.
      call jacprod (bs, h, nots, nb)         !Returns bs(i)=bs(j)*J(j,i)
      jacob(:nots) = af(:nots) + bs(:nots)
 
      end subroutine splinint
      
      subroutine splint(xa, ya, y2a, n, x, y, yp, ndim)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer n, ndim
      real(rprec), dimension(*) :: xa, ya, y2a, x, y, yp
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, c1o6 = 1._dp/6._dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: klo, khi, i, k
      real(rprec) :: h, a, a2, b, b2, h2, y26lo, y26hi, deriv
C-----------------------------------------------
!
!       SPLINE INTERPOLATION ROUTINE (Numerical Recipes, pg. 89)
!       XA: ordered array of length N of ordinates at which function YA=F(XA)
!           is tabulated
!       YA: array of length N , = F(XA)
!       Y2A: array of second derivatives at XA points
!       computed from call to SPLINE
!       X : value at which Y = F(X) is to be computed from splines
!       YP = dY/dX at X
!       NDIM: dimension of X, Y, YP arrays
 
 
      deriv = yp(1)
      klo = 1
      khi = n
      do i = 1, ndim
         do while(khi - klo .gt. 1)
            k = (khi + klo)/2
            if (xa(k) .gt. x(i)) then
               khi = k
            else
               klo = k
            endif
         end do
 
         h = xa(khi) - xa(klo)
         a = xa(khi) - x(i)
         b = x(i) - xa(klo)
         h2 = h*h
         a2 = a*a
         b2 = b*b
         y26lo = c1o6*y2a(klo)
         y26hi = c1o6*y2a(khi)
         y(i) = (a*(ya(klo)+(a2-h2)*y26lo)+b*(ya(khi)+(b2-h2)*y26hi))/h
         if (deriv .ne. zero) yp(i) = (ya(khi)-ya(klo)+y26hi*(3.0*b2-h2)
     1      -y26lo*(3.0*a2-h2))/h
         if (i.lt.ndim .and. x(i+1).gt.x(i)) then
            khi = n
         else
            klo = 1
         endif
      end do
 
      end subroutine splint
      
      function splints (x)
      use vspline
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec) x
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: c1o6 = 1._dp/6._dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: klo, khi, k
      real(rprec) :: h, a, a2, b, b2, h2, y26lo, y26hi, yx, splints
C-----------------------------------------------
 
      klo = 1
      khi = isnodes
 
    1 continue
      if (khi - klo .gt. 1) then
         k = (khi + klo)/2
         if (sknots(k) .gt. x) then
            khi = k
         else
            klo = k
         endif
         go to 1
      endif
 
      h = sknots(khi) - sknots(klo)
      a = sknots(khi) - x
      b = x - sknots(klo)
      h2 = h*h
      a2 = a*a
      b2 = b*b
      y26lo = c1o6*y2stark(klo)
      y26hi = c1o6*y2stark(khi)
      yx = (a*(ystark(klo)+(a2-h2)*y26lo)+b*(ystark(khi)+(b2-h2)*y26hi))
     1   /h
      splints = yx

      end function splints
       
      function splintx(x)
      use vparams
      use vsvd
      use csplinx
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec) x
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: c1o6 = 1._dp/6._dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: klo, khi, k
      real(rprec) :: h, a, b, h2, a2, b2, y26lo, y26hi, qmidx0,
     1   splintx      
C-----------------------------------------------
 
 
       call setspline (rmidx, wmidx, qmidx, hmidx, ymidx, y2midx, 
     1    tenmidx, tenmidx(1), nptsx, natur)
 
      klo = 1
      khi = nptsx
 
    1 continue
      if (khi - klo .gt. 1) then
         k = (khi + klo)/2
         if (rmidx(k) .gt. x) then
            khi = k
         else
            klo = k
         endif
         go to 1
      endif
 
      h = rmidx(khi) - rmidx(klo)
      if( h.eq.zero )then
        splintx = zero
        return
      end if  
      a = rmidx(khi) - x
      b = x - rmidx(klo)
      h2 = h*h
      a2 = a*a
      b2 = b*b
      y26lo = c1o6*y2midx(klo)
      y26hi = c1o6*y2midx(khi)
      qmidx0 = (a*(ymidx(klo)+(a2-h2)*y26lo)+b*(ymidx(khi)+
     1   (b2-h2)*y26hi))/h
      splintx = qmidx0

      end function splintx
EOF
cat > input.demo << "EOF"
!-  5-MAR-95 15:41:59 wieland 83738s03t290 83738 2.90003   sq_all:s_v18s.ang
--------------------------------------------------------------------------------
c  In general, for "on/off" switches, think "1=on", "0=off"
c
c iopt_raxis:     = 1 (default) finds optimum Raxis value automatically
c                 = 0 freezes Raxis at the namelist value
c
c imatch_phiedge: = 1 (default) matches PHIEDGE to namelist value
c                 = 0 attempts to find phiedge automatically so as
c                     to match the minor radius based on the outboard
c                     pressure
c
c *_sigma         > 0 means use this as an absolute sigma
c                 < 0 means apply the abs of this as the percent error
c
--------------------------------------------------------------------------------
 &INDATA
 MGRID_FILE = 'mgrid.tftr                                                  ',
 TIME_SLICE =  0.000E+00, DELT =  1.100E+00,
 FTOL_ARRAY = 1.E-6  5.000E-11, NITER =   600, NSTEP =   100, NVACSKIP =    12,
 MPOL = 6, NS_ARRAY = 16,31, NFP =     1, NCURR =     0, GAMMA =  0.000E+00,
 AI =      2.354498895432447, -9.907517672615517, 42.86570171033211,    
           -122.8693254810872, 188.6838680668763, -143.1773314731515,    
           42.19838185779675,    
 AM =      0.2762993401771816E6, -1.398637555311701E6,
           6.435929178901299E6, -19.13351221515564E6,
           30.46522063279099E6, -23.99696824024823E6,
           7.363942466650127E6,
 PSA=      0.00000E+00,1.00000E-01,2.00000E-01,3.00000E-01,4.00000E-01,
           5.00000E-01, 6.00000E-01, 7.00000E-01, 8.00000E-01, 9.00000E-01,
           1.00000E+00,
 PFA=      3.47209E-01, 2.31819E-01, 1.79043E-01, 1.42202E-01, 1.08563E-01,
           7.98048E-02, 5.91538E-02, 4.51807E-02, 3.22540E-02, 1.76674E-02,
           1.54225E-02,
 ISA=      0.00000E+00, 1.00000E-01, 2.00000E-01, 3.00000E-01, 4.00000E-01,
           5.00000E-01, 6.00000E-01, 7.00000E-01, 8.00000E-01, 9.00000E-01,
           1.00000E+00,
 IFA=      2.35450E+00, 1.68701E+00, 1.36345E+00, 1.13387E+00, 9.23383E-01,
           7.36300E-01, 5.90636E-01, 4.83033E-01, 3.84059E-01, 2.63870E-01,
           1.48272E-01,
 RBC(00,00) =  2.599E+00, ZBS(00,00) =  0.000E+00, RBS(00,00) =  0.000E+00, ZBC(00,00) =  0.000E+00,
 RBC(00,01) =  9.370E-01, ZBS(00,01) =  9.590E-01, RBS(00,01) =  0.000E+00, ZBC(00,01) =  0.000E+00,
 RAXIS(00) =  2.772E+00, ZAXIS(00) =  0.000E+00,
 CURTOR =  1.582E+06,  SIGMA_CURRENT = -7.000E-03, 
 EXTCUR = -6.800E+01, -1.070E+01,  1.435E+01, -1.405E+01,  2.498E-02, 
  0.000E+00,  0.000E+00,  0.000E+00,  0.000E+00,  0.000E+00, 
 PHIEDGE =-1.3534E+01, IMATCH_PHIEDGE =   2, IOPT_RAXIS =   1,
 TENSI =  1.000E+00, TENSP =  1.000E+00,
 MSEANGLE_OFFSET =  0.000E+00, IMSE =  12,
 RSTARK(01) =  2.563E+00, DATASTARK(01) = -1.631E+00, SIGMA_STARK(01) =  2.099E-01,
 RSTARK(02) =  2.619E+00, DATASTARK(02) = -4.520E-01, SIGMA_STARK(02) =  1.209E-01,
 RSTARK(03) =  2.726E+00, DATASTARK(03) =  4.920E-01, SIGMA_STARK(03) =  1.066E-01,
 RSTARK(04) =  2.776E+00, DATASTARK(04) =  6.870E-01, SIGMA_STARK(04) =  1.056E-01,
 RSTARK(05) =  2.825E+00, DATASTARK(05) =  1.188E+00, SIGMA_STARK(05) =  1.056E-01,
 RSTARK(06) =  2.916E+00, DATASTARK(06) =  2.047E+00, SIGMA_STARK(06) =  1.121E-01,
 RSTARK(07) =  2.960E+00, DATASTARK(07) =  2.903E+00, SIGMA_STARK(07) =  1.156E-01,
 RSTARK(08) =  3.003E+00, DATASTARK(08) =  3.657E+00, SIGMA_STARK(08) =  1.149E-01,
 RSTARK(09) =  3.086E+00, DATASTARK(09) =  4.940E+00, SIGMA_STARK(09) =  1.184E-01,
 RSTARK(10) =  3.240E+00, DATASTARK(10) =  6.397E+00, SIGMA_STARK(10) =  2.189E-01,
 RSTARK(11) =  3.383E+00, DATASTARK(11) =  6.974E+00, SIGMA_STARK(11) =  2.250E-01,
 RSTARK(12) =  3.453E+00, DATASTARK(12) =  7.655E+00, SIGMA_STARK(12) =  2.462E-01,
 PRESFAC =  1.000E+00,  PRES_OFFSET =  0.000E+00,  ITSE =  23,
 RTHOM(01) =  2.701E+00, DATATHOM(01) =  1.238E+05, SIGMA_THOM(01) = -6.000E-02,
 RTHOM(02) =  2.749E+00, DATATHOM(02) =  1.260E+05, SIGMA_THOM(02) = -9.000E-02,
 RTHOM(03) =  2.772E+00, DATATHOM(03) =  1.271E+05, SIGMA_THOM(03) = -1.200E-01,
 RTHOM(04) =  2.796E+00, DATATHOM(04) =  1.260E+05, SIGMA_THOM(04) = -9.000E-02,
 RTHOM(05) =  2.842E+00, DATATHOM(05) =  1.238E+05, SIGMA_THOM(05) = -6.000E-02,
 RTHOM(06) =  2.888E+00, DATATHOM(06) =  1.153E+05, SIGMA_THOM(06) = -6.000E-02,
 RTHOM(07) =  2.933E+00, DATATHOM(07) =  1.041E+05, SIGMA_THOM(07) = -6.000E-02,
 RTHOM(08) =  2.977E+00, DATATHOM(08) =  8.969E+04, SIGMA_THOM(08) = -6.000E-02,
 RTHOM(09) =  3.021E+00, DATATHOM(09) =  7.668E+04, SIGMA_THOM(09) = -6.000E-02,
 RTHOM(10) =  3.064E+00, DATATHOM(10) =  6.410E+04, SIGMA_THOM(10) = -6.000E-02,
 RTHOM(11) =  3.106E+00, DATATHOM(11) =  5.245E+04, SIGMA_THOM(11) = -6.000E-02,
 RTHOM(12) =  3.146E+00, DATATHOM(12) =  4.180E+04, SIGMA_THOM(12) = -6.000E-02,
 RTHOM(13) =  3.186E+00, DATATHOM(13) =  3.320E+04, SIGMA_THOM(13) = -6.000E-02,
 RTHOM(14) =  3.224E+00, DATATHOM(14) =  2.617E+04, SIGMA_THOM(14) = -6.000E-02,
 RTHOM(15) =  3.262E+00, DATATHOM(15) =  1.970E+04, SIGMA_THOM(15) = -6.000E-02,
 RTHOM(16) =  3.298E+00, DATATHOM(16) =  1.511E+04, SIGMA_THOM(16) = -6.000E-02,
 RTHOM(17) =  3.332E+00, DATATHOM(17) =  1.083E+04, SIGMA_THOM(17) = -8.000E-02,
 RTHOM(18) =  3.366E+00, DATATHOM(18) =  7.903E+03, SIGMA_THOM(18) = -8.000E-02,
 RTHOM(19) =  3.397E+00, DATATHOM(19) =  5.482E+03, SIGMA_THOM(19) = -1.000E-01,
 RTHOM(20) =  3.428E+00, DATATHOM(20) =  3.802E+03, SIGMA_THOM(20) = -1.100E-01,
 RTHOM(21) =  3.457E+00, DATATHOM(21) =  2.289E+03, SIGMA_THOM(21) = -1.200E-01,
 RTHOM(22) =  3.485E+00, DATATHOM(22) =  1.047E+03, SIGMA_THOM(22) = -1.400E-01,
 RTHOM(23) =  3.511E+00, DATATHOM(23) =  2.544E+02, SIGMA_THOM(23) = -4.000E-01,
 PHIDIAM =  1.307E-02,  SIGMA_DELPHID =  1.000E-03, 
 NFLXS =   6, INDXFLX = 01,02,03,04,05,06,
 DSIOBT(01) =  3.794E-01, SIGMA_FLUX(01) =  3.000E-02,
 DSIOBT(02) =  1.305E+00, SIGMA_FLUX(02) =  3.000E-02,
 DSIOBT(03) = -8.790E-03, SIGMA_FLUX(03) =  3.000E-02,
 DSIOBT(04) = -1.322E+00, SIGMA_FLUX(04) =  3.000E-02,
 DSIOBT(05) = -3.722E-01, SIGMA_FLUX(05) =  3.000E-02,
 DSIOBT(06) = -2.940E-03, SIGMA_FLUX(06) =  3.000E-02,
 NBFLD =  20,   0, 
 INDXBFLD(1,01) = 01, 02, 03, 04, 05, 06, 07, 08, 09, 10, 11, 12, 13, 14, 15, 
 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 
 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 
 BBC(01,01) =  0.000E+00, SIGMA_B(01,01) =  1.200E-03,
 BBC(02,01) = -2.790E-01, SIGMA_B(02,01) =  4.300E-03,
 BBC(03,01) =  5.420E-02, SIGMA_B(03,01) =  2.100E-03,
 BBC(04,01) = -2.770E-01, SIGMA_B(04,01) =  5.500E-03,
 BBC(05,01) =  1.030E-01, SIGMA_B(05,01) =  3.600E-03,
 BBC(06,01) = -2.580E-01, SIGMA_B(06,01) =  6.100E-03,
 BBC(07,01) =  1.490E-01, SIGMA_B(07,01) =  5.700E-03,
 BBC(08,01) = -2.250E-01, SIGMA_B(08,01) =  6.400E-03,
 BBC(09,01) =  1.950E-01, SIGMA_B(09,01) =  9.100E-03,
 BBC(10,01) = -1.820E-01, SIGMA_B(10,01) =  6.400E-03,
 BBC(11,01) =  2.290E-01, SIGMA_B(11,01) =  1.220E-02,
 BBC(12,01) = -1.360E-01, SIGMA_B(12,01) =  4.800E-03,
 BBC(13,01) =  2.500E-01, SIGMA_B(13,01) =  1.490E-02,
 BBC(14,01) = -8.260E-02, SIGMA_B(14,01) =  2.500E-03,
 BBC(15,01) =  2.500E-01, SIGMA_B(15,01) =  1.270E-02,
 BBC(16,01) =  5.140E-02, SIGMA_B(16,01) =  5.000E-03,
 BBC(17,01) =  2.200E-01, SIGMA_B(17,01) =  9.300E-03,
 BBC(18,01) =  1.110E-01, SIGMA_B(18,01) =  6.500E-03,
 BBC(19,01) =  1.740E-01, SIGMA_B(19,01) =  6.400E-03,
 BBC(20,01) =  1.800E-01, SIGMA_B(20,01) =  7.200E-03,
 BBC(21,01) =  1.230E-01, SIGMA_B(21,01) =  2.700E-03,
 BBC(22,01) =  2.380E-01, SIGMA_B(22,01) =  4.800E-03,
 BBC(23,01) =  6.350E-02, SIGMA_B(23,01) =  1.600E-03,
 BBC(24,01) =  2.790E-01, SIGMA_B(24,01) =  4.200E-03,
 BBC(25,01) =  5.600E-03, SIGMA_B(25,01) =  1.500E-03,
 BBC(26,01) =  2.950E-01, SIGMA_B(26,01) =  4.100E-03,
 BBC(27,01) = -5.820E-02, SIGMA_B(27,01) =  1.800E-03,
 BBC(28,01) =  2.830E-01, SIGMA_B(28,01) =  4.200E-03,
 BBC(29,01) = -1.180E-01, SIGMA_B(29,01) =  2.700E-03,
 BBC(30,01) =  2.440E-01, SIGMA_B(30,01) =  4.800E-03,
 BBC(31,01) = -1.710E-01, SIGMA_B(31,01) =  6.400E-03,
 BBC(32,01) =  1.820E-01, SIGMA_B(32,01) =  7.200E-03,
 BBC(33,01) = -2.210E-01, SIGMA_B(33,01) =  9.300E-03,
 BBC(34,01) =  1.150E-01, SIGMA_B(34,01) =  6.500E-03,
 BBC(35,01) = -2.480E-01, SIGMA_B(35,01) =  1.270E-02,
 BBC(36,01) =  4.980E-02, SIGMA_B(36,01) =  5.000E-03,
 BBC(37,01) = -2.550E-01, SIGMA_B(37,01) =  1.820E-02,
 BBC(38,01) = -1.060E-02, SIGMA_B(38,01) =  2.700E-03,
 BBC(39,01) = -2.660E-01, SIGMA_B(39,01) =  1.630E-02,
 BBC(40,01) = -7.570E-02, SIGMA_B(40,01) =  2.500E-03,
 BBC(41,01) = -2.380E-01, SIGMA_B(41,01) =  1.280E-02,
 BBC(42,01) = -1.370E-01, SIGMA_B(42,01) =  5.100E-03,
 BBC(43,01) = -1.930E-01, SIGMA_B(43,01) =  9.200E-03,
 BBC(44,01) = -1.810E-01, SIGMA_B(44,01) =  6.500E-03,
 BBC(45,01) = -1.480E-01, SIGMA_B(45,01) =  5.700E-03,
 BBC(46,01) = -2.220E-01, SIGMA_B(46,01) =  6.300E-03,
 BBC(47,01) = -1.020E-01, SIGMA_B(47,01) =  3.500E-03,
 BBC(48,01) = -2.560E-01, SIGMA_B(48,01) =  6.000E-03,
 BBC(49,01) = -5.550E-02, SIGMA_B(49,01) =  2.000E-03,
 BBC(50,01) = -2.740E-01, SIGMA_B(50,01) =  5.400E-03,
 /
 &END
EOF
cat > input.demo3d << "EOF"
 &INDATA
 MGRID_FILE = 'NONE',
 DELT = 1.0,
 NFP     =          12,
 NCURR   =           0,
 NITER   =        3000,
 NSTEP   =        200,
 NVACSKIP        =          12,
 FTOL_ARRAY    = 1.e-8,  1.E-11,
 NS_ARRAY    =  26, 51
 MPOL = 6  NTOR = 6
 GAMMA   =   1.667,
 PHIEDGE = 0.2300000    ,
 CURTOR  =  0.0000000E+00,
 EXTCUR = -3.9385E2, 29.4795, -61.3872, -3.133E-2,
 AM      =  3.6500000E+04, -7.300000E+04,  3.650000E+04, 8*0.0000000E+00,
 AI      =  0.3000000,  0.4200000,  0.2500000, 8*0.0000000E+00,
 RAXIS   =   1.7200000,
 ZAXIS   = 0.0000000E+00,
 RBC(0,0) =   1.7183E+00    ZBS(0,0) =   0.0000E+00
 RBC(0,1) =   2.1404E-01    ZBS(0,1) =   2.5040E-01
 RBC(1,1) =  -5.5576E-02    ZBS(1,1) =   6.2457E-02
 RBC(-1,2)=   9.8162E-04    ZBS(-1,2)=   1.3940E-03
 RBC(0,2) =   4.0247E-03    ZBS(0,2) =  -2.6518E-03
 RBC(1,2) =  -6.1395E-03    ZBS(1,2) =   1.0371E-02
 RBC(0,3) =  -2.2444E-04    ZBS(0,3) =  -2.2444E-04
 RBC(1,3) =  -7.2860E-04    ZBS(1,3) =  -7.2860E-04
 RBC(0,4) =   2.0864E-04    ZBS(0,4) =   2.0864E-04
 RBC(1,4) =  -5.1596E-04    ZBS(1,4) =  -5.1596E-04
 /
 &END
EOF
cat > input.sequence << "EOF"
!SAMPLE SEQUENCE FILE FOR RUNNING VMEC SEQUENTIALLY
 &VSEQ
 NSEQ = 2,
 NSEQ_SELECT = 1,2,3,
 EXTENSION = 't01', 't02', 't03'
 /
 &END
EOF
EOC
