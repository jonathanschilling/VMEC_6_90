!!!!
!                              D   I   S   C   L   A   I   M   E   R
!
!       You are using a BETA version of the program DKES, which is currently
!       under development by S. P. Hirshman at the Fusion Energy Division,
!       Oak Ridge National Laboratory.  Please report any problems or comments
!       to him.  As a BETA version, this program is subject to change
!       and improvement without notice.
!
!!!!
      program dkes2
c-----------------------------------------------------------------------
c     SPH:         implemented in UNIX, F77                 February 21, 1994
c     SPH:         upgraded to F90, dynamic memory allocation,
c                  matmul routines                          May 1999
c-----------------------------------------------------------------------
c
c
c   Example usage:
c
c   OR   xdkes INPUT_FILE_NAME
c
c        xdkes BOOZ_FILE_NAME NSURF CMUL EFIELD LSCREEN
c
c
c  WHERE
c
c        INPUT_FILE_NAME:   name of the input file prepared by dkes calling
c                           subroutine dkes_input_prepare (the old way of running DKES)
c
c        BOOZ_FILE_NAME:    name of Boozermn file; see dkes_input_prepare subroutine for
c                           description of this and nsurf, cmul, efield, lscreen parameters
c
c
c  DKES.VAR (Drift Kinetic Equation Solver, Variational) solves a set
c  of 3-D drift kinetic equations to obtain upper and lower bounds for
c  the diffusion coefficients of a prescribed toroidal plasma
c  equilibrium. The 3 dimensions are theta (poloidal angle),
c  zeta (toroidal angle), and alpha (pitch angle). Straight-line flux
c  coordinates are used to describe the equilibrium, which satisfies
c  the stellarator symmetry conditions R(theta,zeta) = R(-theta,-zeta),
c  Z(theta,zeta) = - Z(-theta,-zeta).
c
c  Reference: W. I. van Rij and S. P. Hirshman, Variational Bounds for
c  Transport Coefficients in Three-Dimensional Plasmas,
c  Phys. Fluids B 1,563(1989).
c
c  Boozer Coordinate Version
c
c-----------------------------------------------------------------------
c
c  mpnt      : The total number of distribution modes is computed from
c              the mmnn matrix in the MNSET subroutine
c  iswpm = 1 : "Plus" sources, distributions used  (Maximizing)
c  iswpm = 2 : "Minus" sources,distributions used  (Minimizing)
c
c-----------------------------------------------------------------------
c
c  subroutines required:
c
c       name:     purpose:
c
c  DKES2 essential routines:
c       blk5d     solves the block-pentadiagonal system of
c                 equations
c       blox      forms the l-row block matrices and rows
c       cescale   cmul and efield scaling and
c                 the dominant-diagonal scaling arrays
c       ftconv    computes magnetic field spectral arrays and Fourier transforms
c       lcalc     degree of Legendre polynomial arrays
c       printout  calculate diffusion coefficients and output
c       residue   calculates the solution residuals
c       reverse   returns solution arrays to their original
c                 order with respect to l
c       scalel    dominant-diagonal l matrix scaling
c       wrout     output of magnetic field and of final
c                 distribution functions
c
c  DKES2 auxiliary routines:
c       rddisk    auxilary routine for fast internal disk i/o
c       wrdisk    auxilary routine for fast internal disk i/o
c
c  LINPACK (BLAS) routines:
c       GETRF_G, GETRS_G                    GENERIC (_G) LAPACK ROUTINES WITH CORRECT PRECISION
c
c-----------------------------------------------------------------------
c
c     main program
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      use Vimatrix
      use Vnamecl2
      use Vcmain
      use dkes_input
      use dkes_realspace
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, ir, irun, istat, neqs
      real(rprec), dimension(:), pointer :: f0p1, f0p2, f0m1, f0m2
      real(rprec) :: tcpu0, tcpu1, tcpui, tcput, tcpu, tcpua
C-----------------------------------------------------------------------
      call second0 (tcpu0)
!
!     OPEN FILES FOR IO AND READ INPUT FILE DATA
!
      call read_dkes_input

!
!     set-up : calculate spectral arrays, matrix elements, and l arrays
!
      call set_mndim

      neqs = mpnt*(lalpha + 1)                                           !add l = -1 "constraint"
      allocate (fzerop(2*neqs), fzerom(2*neqs), stat=istat)
      if (istat .ne. 0) stop 'allocation error(1) in dkes2!'


      f0p1 => fzerop(1:neqs);   f0p2 => fzerop(neqs+1:)                  !+ functions
      f0m1 => fzerom(1:neqs);   f0m2 => fzerom(neqs+1:)                  !- functions

      call ftconv

      call lcalc

      call second0 (tcpu1)
      tcpui = tcpu1 - tcpu0

c  header for file DKESOUT; output magnetic field

      call header

      call second0 (tcpu0)
      tcput = zero

c  cmul, efield loops
      nrun = max(nrun, 1)
      nrun = min(nrun, krun)

      do ir = 1, nrun
         irun = ir
         efield1 = efield(irun)
         cmul1 = cmul(irun)

c  cmul and efield scaling; scaling arrays

         call cescale (srces0)

         write (ioout, 950) dashes, cmul1, efield1, weov, wtov,
     1          wcyclo, vthermi
  950    format(/9x,'CMUL',7x,'EFIELD',4x,'OMEGA-E/v',4x,'OMEGA-T/v',
     1         5x,'OMEGA-ci',3x,'VI-THERMAL'/
     2         7x,'(nu/v)',7x,'(Es/v)',2x,'(ExB drift)',4x,'(transit)',
     3         3x, '(H+, B=B00)',2x,'(H+, Ti=1keV)'/,a,/
     4         1p6e13.4/)

c  block-pentadiagonal solutions and residuals

         if (ipmb < 2) then
            iswpm = 1
            call blk5d (blk1, blk2, blk3, blk4, blk5, blk6, 
     1               blk7, f0p1, f0p2, srces0)

            if (ier .ne. 0) then
               write (ioout, 1000) ier
 1000          format(/' blk5d error in block = ',i5)
               stop
            endif
            call residue (blk1, blk2, blk3, blk4, blk5, blk6, 
     1           f0p1, f0p2, srces0, rsd1p, rsd3p, g11p, g33p, 
     2           g31p, g13p, crs1p, crs3p)
         endif

         if (ipmb .ne. 1) then
            iswpm = 2
            call blk5d (blk1, blk2, blk3, blk4, blk5, blk6, 
     1         blk7, f0m1, f0m2, srces0)

            if (ier .ne. 0) then
               write (ioout, 1050) ier
 1050          format(/' blk5d error : ierm = ',i5)
               stop
            endif
            call residue (blk1, blk2, blk3, blk4, blk5, blk6, 
     1           f0m1, f0m2, srces0, rsd1m, rsd3m, g11m, g33m, 
     2           g31m, g13m, crs1m, crs3m)
         endif

c  calculate and output diffusion coefficients

         call printout (f0p1, f0m1, f0p2, f0m2, srces0)

c  timing and check remaining run time

         call second0 (tcpu1)
         tcpu = tcpu1 - tcpu0
         tcpu0 = tcpu1
         tcput = tcput + tcpu
         tcpua = tcput/irun
         write (ioout, 1100) tcpu
!        if (irun<nrun .and. tcpu1<1.05*tcpua+3.0) exit
 1100    format(/' time used:    tcpu =',1p,e10.2,'  sec'/)
      end do


c  output magnetic field and final distributions

      if (lfout .ne. 0) call wrout (f0p1, f0m1, f0p2, f0m2,
     1   srces0)



c  clean-up memory
      call free_mndim
      deallocate (cols, al1, al2, al3, al4, bl1, bl2, bl3, bl4, cl1,
     1   cl2, cl3, cl4, cols0, omgl, al01, al02, al03, al04, bl01, 
     2   bl02, bl03, bl04, cl01, cl02, cl03, cl04, fzerop, fzerom)

c  time and date

      tcput = tcput + tcpui
      write (ioout, 1200) dashes, tcpui, tcpua, tcput
      if (irun < nrun) then
         write (ioout, '(a,i3,a,i3)')
     1      'DKES2 run not completed (time limit):   irun = ',
     2      irun, '   nrun =', nrun
      else
         write (ioout, '(a,i3)') ' DKES2 run completed: nrun = ', nrun
      endif

      write (ioout, '(1x,a)') dashes
 1200 format(1x,a/1p,' tcpui =',e9.2,' s',5x,'tcpua =',e9.2,' s',5x,
     1    'tcput =',e9.2,' s')

      end program dkes2


      subroutine blk5d(a, bm1, bp1, cp2, cm2back, cm2, pl2back,
     1   fz1, fz3, srces)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      use kind_spec
      use Vimatrix
      use Vnamecl2
      use dkes_input, only: idisk, lalpha
      use dkes_realspace, only: mn0, mpnt, mpntsq, diagle, diagl
      use safe_open_mod
      implicit none
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: bytes_per_rprec = 8
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(mpnt,mpnt), intent(out) ::
     1     a, bp1, bm1, cp2, cm2back, cm2, pl2back
      real(rprec), dimension(mpnt,0:lalpha), intent(out) :: fz1, fz3
      real(rprec), dimension(mpnt,lsource,2,2), intent(in) :: srces
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ibuph, incnow, irecl, incbu, iunit
      integer :: kbot, k, ll, mn, mnp, mblk2
      integer, allocatable, dimension(:) :: ipiv
      real(rprec) :: fac
      real(rprec), dimension(:,:), pointer :: ql, pl, ql2back
      real(rprec), dimension(:,:,:), pointer :: transf
      logical :: ldisk
C-----------------------------------------------
C   E x t e r n a l  S u b r o u t i n e s
C-----------------------------------------------
      external GETRF_G, GETRS_G
C-----------------------------------------------
c
c  original version:                   W.I. van Rij
c                         in the variational version of DKES code
c
c  modified for CRAY-XMP (Garching):   H. Maassberg       May 89
c                  ---->  FORTRAN subroutine MC32AD
c                  ---->  simulation of BLAS routines SGEMM/SGEMV
c                         by routine SAXPY
c                  ---->  disk i/o by routines WRDISK / RDDISK
c
c  modified (May, 1999, ORNL):         S. P. Hirshman
c                  Removed MC32AD, replaced with F90 matmul routines
c
c  significantly modified (June, 2001, ORNL) S. P. Hirshman
c                  fixed bug associated with particle conservation
c                  rewrote block pentadiagonal solver in readable fashion
c                  and to incorporate l=0 constraints as a "ghost" l=-1 component
c                  rewrote blox routine
c
c-----------------------------------------------------------------------
c
c  This subroutine solves the block-pentadiagonal system of equations.
c  It is called once with "plus" indexing for fz1,fz3 and once with "minus"
c  indexing, corresponding to maximizing, minimizing bounding distributions
c  for use in the variational equations
c
c-----------------------------------------------------------------------
c
c  iunit       : unit number for block-pentadiagonal solution disk file.
c  fz1(mn,l)   : distribution function in response to density gradients (1)
c  fz3(mn,l)   : distribution function in response to parallel electric field (3)
c                Note: the l=0 (real l=-1) component is the associated F(-,+)(l=0)
c                component needed for particle conservation constraint, previously
c                call fo1 and fo3
c
c  srces(mn,l,itype,plus-min)
c              : sources, mn=Fourier index;       l=Legendre index (+1);
c                itype=1,3 (uses 2) (n'',<E.B>);   plus(+) or minus(-) type (max,min)
c
c  Distributions are indexed in m-n Fourier-space, Legendre-space. The penta-diagonal
c  equation is:
c
c  cm2 * f(l-2) + bm1 * f(l-1) + a * f(l) + bp1 * f(l+1) + cp2 * f(l+2) = source(l)
c
c  where:   bp1(l) = bm1(l+1) (transpose);     cp1(l) = cm1(l+2) (transpose)
c
c
c
c     GENERAL SOLUTION SCHEME APPLIED TO EACH BLOCK (L)
c
c     1. CALL GETRF_G:   Perform LU factorization of diagonal block (A)
c     2. CALL GETRS_G:   With multiple (mpnt) right-hand sides, to do block inversion
c                        operation, A X = B  (store result in X; here B is a matrix)
c     3. CALL GETRS_G:   With single right hand side (source) to solve A x = b (b a vector)
c

      mblk2 = 2*mpntsq
      allocate (transf(mpnt,mpnt,2), ql2back(mpnt,mpnt), ipiv(mpnt), 
     1          stat=ier)
      if (ier .ne. 0) stop 'Allocation error in blk5'

      pl => transf(:,:,1)
      ql => transf(:,:,2)

      kbot = idisk*lam6 + 1


c  create disk file for doing direct access i/o.

      incnow = mblk2
      irecl  = bytes_per_rprec*incnow
      incbu = 1 + (mblk2 - 1)/incnow
      ibuph = 0

      iunit = 10
      call safe_open(iunit, ier, 'NULL', 'scratch', 'unformatted', 
     1     irecl, 'DIRECT')
      if (ier .ne. 0) stop 'Error opening scratch file in blk5d DKES'
      
c  load and process first block-row. -----------------------------------

      k = 1
      ll = lalpha

      call blox (ll, a, bm1, cm2, fz1(1,ll), fz3(1,ll), srces)
      ql = bm1
      pl = cm2
      cm2back = cm2

!
!     Compute (and save) ql = A-1 ql,  pl = A-1 pl, and source terms A-1 fz1,3
!
      call GETRF_G (mpnt, mpnt, a, mpnt, ipiv, ier)
      if (ier .ne. 0) go to 200

      call GETRS_G ('n', mpnt, mpnt, a, mpnt, ipiv, ql, mpnt, ier)
      call GETRS_G ('n', mpnt, mpnt, a, mpnt, ipiv, pl, mpnt, ier)

      call GETRS_G ('n', mpnt, 1, a, mpnt, ipiv, fz1(1,ll), mpnt, ier)
      call GETRS_G ('n', mpnt, 1, a, mpnt, ipiv, fz3(1,ll), mpnt, ier)

!
!     save pl as pl2back (will use at l+2 iteration)
!
      pl2back = pl
      ql2back = ql

      ldisk = (idisk.eq.0) .or. ((idisk.eq.1) .and. (lalpha.le.7))

      if (ldisk) call wrdisk(iunit, transf, mblk2, 
     1          incnow, ibuph, incbu, ier)
      if (ier .ne. 0) go to 302


c  load and process second block-row. ----------------------------------

      k = k+1
      ll = ll-1
      bp1 = transpose(bm1)

      call blox (ll, a, bm1, cm2, fz1(1,ll), fz3(1,ll), srces)

      a   = a - matmul(bp1, ql)
      ql  = bm1 - matmul(bp1, pl)
      pl  = cm2

      fz1(:,ll) = fz1(:,ll) - matmul(bp1, fz1(:,ll+1))
      fz3(:,ll) = fz3(:,ll) - matmul(bp1, fz3(:,ll+1))

      call GETRF_G (mpnt, mpnt, a, mpnt, ipiv, ier)
      if (ier .ne. 0) go to 200

      call GETRS_G ('n', mpnt, mpnt, a, mpnt, ipiv, ql, mpnt, ier)
      call GETRS_G ('n', mpnt, mpnt, a, mpnt, ipiv, pl, mpnt, ier)

      call GETRS_G ('n', mpnt, 1, a, mpnt, ipiv, fz1(1,ll), mpnt, ier)
      call GETRS_G ('n', mpnt, 1, a, mpnt, ipiv, fz3(1,ll), mpnt, ier)

      if (ldisk) call wrdisk(iunit, transf, mblk2, 
     1          incnow, ibuph, incbu, ier)
      if (ier .ne. 0) go to 302

c  main loop. load and process block-rows 3 to lalpha+1. The last row (k=lalph+1) corresponds
c  to the constraint function F- which has been added to the f vector to satisfy particle
c  conservation

      fac = iswpm - 1

      BLOCKS: do k = 3, lalpha + 1
         ll = ll - 1

         if (k .le. lalpha) then
            bp1 = transpose(bm1)
            cp2 = transpose(cm2back)
            cm2back = cm2                                  !!stored 2 blocks back

            call blox (ll, a, bm1, cm2, fz1(1,ll), fz3(1,ll), srces)
         else                                               
!
!         particle conservation constraint; note V(l=0)F changes parity (+,-)
!         which is why the -transpose is used here
!
            a = 0
            a(mn0, mn0) = 1
            bp1 =-transpose(diagle(:,:,iswpm))
            cp2 =-transpose(diagl(:,:,iswpm))                        !f(l=1) part of V(l=0)
            fz1(:,0) = fac*srces(:,1,1,1)
            fz3(:,0) = 0
         end if

!
!      Update diagonal "a" matrix and source terms and store pl,ql 2 l-steps back
!
         bp1 = bp1 - matmul(cp2, ql2back)
         a   = a - matmul(bp1, ql) - matmul(cp2, pl2back)
         ql2back = ql
         pl2back = pl

         fz1(:,ll) = fz1(:,ll) - matmul(bp1,fz1(:,ll+1))
     1                         - matmul(cp2,fz1(:,ll+2))
         fz3(:,ll) = fz3(:,ll) - matmul(bp1,fz3(:,ll+1))
     1                         - matmul(cp2,fz3(:,ll+2))

         if (k .gt. lam1) then
            cm2 = 0
            bm1 = 0
         end if
!
!        Compute (-,+)V[F(l=0)] contributions from l = 1, l = 0 Legendre moments (of V)
!         
         if (k .eq. lam1) cm2 = diagl(:,:,iswpm)
         if (k .eq. lalpha) bm1 = diagle(:,:,iswpm)
         
!
!        Compute a-1; pl = A-1 * pl,  ql = A-1 * ql; sources = A-1 * sources
!
         call GETRF_G (mpnt, mpnt, a, mpnt, ipiv, ier)
         if (ier .ne. 0) go to 200

         if (k .le. lalpha) then
            ql = bm1 - matmul(bp1, pl)
            pl = cm2
            call GETRS_G('n', mpnt, mpnt, a, mpnt, ipiv, ql, mpnt, ier)
            call GETRS_G('n', mpnt, mpnt, a, mpnt, ipiv, pl, mpnt, ier)

            call GETRS_G('n',mpnt, 1, a, mpnt, ipiv, fz1(1,ll),mpnt,ier)
            call GETRS_G('n',mpnt, 1, a, mpnt, ipiv, fz3(1,ll),mpnt,ier)

            ldisk = (idisk.eq.0) .or. ((idisk.eq.1) .and. (k.gt.lam6))
            if (ldisk) call wrdisk(iunit, transf, mblk2, 
     1          incnow, ibuph, incbu, ier)
            if (ier .ne. 0) go to 302
         else
            call GETRS_G('n',mpnt, 1, a, mpnt, ipiv, fz1(1,ll),mpnt,ier)
            call GETRS_G('n',mpnt, 1, a, mpnt, ipiv, fz3(1,ll),mpnt,ier)
         end if
         
      end do BLOCKS

 
c  backward solution sweep for block-rows ll = 1 (l=0) to ll = lap1-kbot

      do k = lalpha, kbot, -1
         ll = lap1 - k
c  read blocks transf => (pl,ql) from disk.

         call rddisk (iunit, transf, mblk2, incnow, ibuph, ier)
         if (ier .ne. 0) then
            write (ioout, '(a)')
     1         ' BLK5D:   error in I/O routine RDDISK'
            go to 303
         endif
         ibuph = ibuph - incbu
         

         fz1(:,ll) = fz1(:,ll) - matmul(ql,fz1(:,ll-1)) 
         fz3(:,ll) = fz3(:,ll) - matmul(ql,fz3(:,ll-1))
         
         if (ll .ge. 2) then
            fz1(:,ll) = fz1(:,ll) - matmul(pl,fz1(:,ll-2)) 
            fz3(:,ll) = fz3(:,ll) - matmul(pl,fz3(:,ll-2)) 
         end if
      end do

      go to 400

c  error returns. ------------------------------------------------------

  200 continue
      ier = ll-1
      go to 400
  301 continue
      write (ioout, '(a,i8)') ' BLK5D:   error in opening file:  ',
     1   'RECL = ', irecl
  302 continue
      write (ioout, '(a)') ' BLK5D:   error in I/O routine WRDISK'
  303 continue
      ier = -2
  305 continue
      write (ioout, '(2/a,i4,2/)') ' BLK5D:   error detected:   ier =',
     1   ier
      stop

c  destroy disk file and return. ---------------------------------------

  400 continue

      close (iunit)

      deallocate (transf, ql2back, ipiv)

      end subroutine blk5d


      subroutine blox(ll, al, bl, cl, s1, s3, srces)

c  This subroutine forms the l-row block matrices and sources.
c
c  On exit, s1 and s3 contain the sources for Legendre index ll
c           a, b, c are the MPNT X MPNT Fourier blocks for this index
c
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      use Vnamecl2
      use dkes_realspace
      implicit none
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: l0 = 1, l1 = 2, l2 = 3,
     1   pgrad = 1, epar = 2
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: ll
      real(rprec), intent(out), dimension(mpnt,mpnt) :: al, bl, cl
      real(rprec), intent(out), dimension(mpnt) :: s1, s3
      real(rprec), intent(in), dimension(mpnt,lsource,2,2) :: srces
C-----------------------------------------------
      if (ll .ge. l2) then
         cl = cl1(ll)*bmat1(:,:) + cl2(ll)*bmat2(:,:,iswpm) 
     1      + cl3(ll)*bmat4(:,:,iswpm)
     2      + cl4(ll)*transpose(bmat4(:,:,iswpm))
      end if

      if (ll .ge. l1) then
         bl = bl1(ll)*bmat5(:,:,iswpm)
     1      + bl2(ll)*transpose(bmat5(:,:,iswpm)) + bl3(ll)
     2       * bmat6(:,:,iswpm) + bl4(ll)*transpose(bmat6(:,:,iswpm))
      end if

      al = al1(ll)*bmat1(:,:) + al2(ll)*bmat2(:,:,iswpm)  
     1      + al3(ll)*bmat3(:,:,iswpm) + cols(ll)*matjac(:,:,iswpm)
     2      + al4(ll)*(bmat4(:,:,iswpm) + transpose(bmat4(:,:,iswpm)))

      if (iswpm.eq.1 .or. ll.eq.l0) al(mn0, mn0) = 1

c  sources (ll <= lsource only)

      if (ll .le. lsource) then
         s1 = srces(:,ll,pgrad,iswpm)
         s3 = srces(:,ll,epar, iswpm)
      else
         s1 = 0
         s3 = 0
      end if

      end subroutine blox


      subroutine cescale(srces)

c  This subroutine performs cmul and efield scaling, and calculates the
c  dominant - diagonal scaling arrays.

C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      use Vnamecl2
      use dkes_input, only: lalpha, psip
      use dkes_realspace
      implicit none
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: l0 = 1, l1 = 2, l2 = 3, l3 = 4,
     1   epar = 2, pgrad = 1, plus = 1, minus = 2
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(mpnt,lsource,2,2), intent(out) :: srces
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: mn, n1
      real(rprec) :: cmul2, eoc, e2oc, sqrt10, bmod
C-----------------------------------------------
      sqrt10 = sqrt(10._dp)
      
      if (psip .ne. zero) then
         weov = efield1/psip                               !!EXB drift frequency (normed to 1/v)
      else
         weov = efield1/epsilon(psip)
      end if

      cmul2 = cmul1*cmul1
      wcyclo = 9.58e7_dp * b00                          !!b00 in [Tesla]
      vthermi = 9.79e3_dp * sqrt(2._dp) * sqrt(1.e3_dp) !!vi for Ti=1keV, [m/s]

c  Spitzer function contribution to conductivity
c  Note that fspitzer = qB/nu/sqrt(bsqav), S3 = qB*jacobian/sqrt(bsqav), so that
c  g33s ~ int(q**2)/nu * vp.
      g33s = one/(cmul1*rt3o2**2)                       !!1./rt3o2**2 = 1/(2/3) from pitch integral of q**2

      eoc  = -efield1/cmul1
      e2oc = -efield1*eoc

      cols = cmul1*cols0(1:lalpha)
      al1 = al01/cmul1
      al2 = al02/cmul1
      al3 = e2oc*al03
      al4 = al04/cmul1
      bl1 = eoc*bl01
      bl2 = eoc*bl02
      bl3 = eoc*bl03
      bl4 = eoc*bl04
      cl1 = cl01/cmul1
      cl2 = cl02/cmul1
      cl3 = cl03/cmul1
      cl4 = cl04/cmul1

      bmod = sqrt(bsqav)
      s1cs1 = s1cs10/cols(l2)

c     e-field dependent particle conservation matrix elements
      diagle = 0
      do mn = 1,mpnt
         diagle(mn,mn,1) = efield1*exbgrad(mn)
         diagle(mn,mn,2) =-efield1*exbgrad(mn)
      end do

c     itype = pgrad  (density gradient sources)
      srces(:,l1,pgrad,minus) = omgl(l2)*auxs1(:,1)/cols(l2)
      srces(:,l2,pgrad,minus) = auxs1(:,2)*eoc/cols0(l2)
      srces(:,l3,pgrad,minus) = omgl(l3)*auxs1(:,3)/cols(l2)


c     itype = epar (parallel E-field sources)
      srces(:,l1,epar,plus) = auxs3p(:,1)*eoc/bmod                       !!l=1 component of s3
      srces(:,l2,epar,plus) = 3*omgl(l2)*auxs3p(:,2)/cmul1/bmod          !!l=2 component of s3
      srces(:,l1,epar,minus) = auxs3m(:)/bmod

      end subroutine cescale


      subroutine header
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      use Vnamecl2
      use Vimatrix, ONLY: ioout
      use dkes_input
      use dkes_realspace, only: mvalue, nvalue, mpnt, borbi1
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: m, n, mm, mn, mnmax
      integer :: index(mpnt)
C-----------------------------------------------
      write (ioout, 10) dashes, chip, psip, btheta, bzeta, vp,
     1   nzperiod, mpolb, ntorb, dashes, blabl(ibbi), blabl(ibbi),
     2   blabl(ibbi), blabl(ibbi)
 10   format(50x,'VARIATIONAL DKES-II CODE - 06/2001'/,
     1       54x,'UNITS ARE ASSUMED TO BE MKS'//1x,a/8x,
     1   'CHIP',8x,'PSIP',6x,'BTHETA',7x,'BZETA',10x,"V'",4x,
     2   'NZPERIOD',7x,'MPOLB',7x,'NTORB',/
     3   1p5e12.4,3i12/1x,a/4x,4('N',4x,'M',14x,a3,9x))

      mnmax = 0
      do mn = 1, mpnt
         if (borbi1(mn) .eq. zero) cycle
         mnmax = mnmax + 1
         index(mnmax) = mn
      end do
      write (ioout, 20) (nvalue(index(mn)), mvalue(index(mn)),
     1   borbi1(index(mn)), mn = 1,mnmax)
 20   format(1p,4(2i5,5x,e12.4,5x))
      write (ioout, '(1x,a)') dashes

c  output numerical parameters and Fourier spectrum

      write (ioout, 30) mpol, ntor, lalpha, mpnt, meshtz, idisk,
     1   ipmb, dashes
 30   format(8x,'MPOL',8x,'NTOR',6x,'LALPHA',2x,'BLOCK SIZE',
     1   6x,'MESHTZ',7x,'IDISK',9x,'IPMB'/7i12/1x,a/4x,'N')
      do n = 1, ntor
         mm = mmnn(2,n)
         write (ioout, 40) mmnn(1,n), (mmnn(2+m,n),m=1,mm)
      end do
 40   format(i5,'    M =',20i6/12x,20i6/12x,20i6)
      write (ioout, '(1x,a)') dashes
      end subroutine header


      subroutine lcalc
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      use Vnamecl2
      use dkes_realspace, only: diagl
      use dkes_input, only: lalpha
      implicit none
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: l1 = 2
      real(rprec), parameter :: half = 0.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: l, ll
      real(rprec) :: el, ap, am
      real(rprec), dimension(:), allocatable :: cinv
C-----------------------------------------------
c  This subroutine calculates l (degree of pitch-angle Legendre
c  polynomial) arrays that are stored once and reused each cmul,efield run.
      lap1 = lalpha + 1
      lam1 = lalpha - 1
      lam3 = lalpha - 3
      lam6 = lalpha - 6
 
      allocate (omgl(0:lap1), cols(lalpha), al1(lalpha), al2(lalpha),
     1   al3(lalpha), al4(lalpha), bl1(lalpha), bl2(lalpha),
     2   bl3(lalpha), bl4(lalpha), cl1(lalpha), cl2(lalpha),
     3   cl3(lalpha), cl4(lalpha), cols0(lap1), al01(lalpha),
     4   al02(lalpha), al03(lalpha), al04(lalpha), bl01(lalpha),
     5   bl02(lalpha), bl03(lalpha), bl04(lalpha), cl01(lalpha),
     6   cl02(lalpha), cl03(lalpha), cl04(lalpha), cinv(0:lap1),
     9   stat=l)
      if (l .ne. 0) stop 'allocation error in LCALC'

      omgl(0) = 0;   omgl(1) = 0;   cols0(1) = 0
      cinv(0) = 0;   cinv(1) = 0
      do l = 2, lap1
         ll = l - 1
         omgl(l)  = half*ll/sqrt(4*ll*ll - one)
         cols0(l) = half*ll*l                                            !.5*l*(l+1)
         cinv(l)  = one/cols0(l)                                         !1/nu(l), l>0
      end do

      do l = 1, lalpha
         ll = l - 1
         ap = cinv(l+1)*omgl(l+1)**2
         am = cinv(l-1)*omgl(l)**2
!
!        coefficients of BMAT1 - BMAT6 comprising A(l) [multiplies f(l)]
!
         al01(l) = 4*(ap + am)
         al02(l) = ap*ll**2 + am*l**2
         al03(l) = cinv(l)
         al04(l) = 2*(am*l - ap*ll)
         bl01(l) = 2*cinv(l-1)*omgl(l)
         bl02(l) = 2*cinv(l)*omgl(l)
         bl03(l) = cinv(l-1)*l*omgl(l)
         bl04(l) =-cinv(l)*(ll - 1)*omgl(l)
         cl01(l) = cinv(l-1)*omgl(l-1)*omgl(l)                           !w(l) * w(l-1) /-nu(l-1)
!
!        coefficients of BMAT1 - BMAT6 comprising C-(l) [multiplies f(l-2)]
!
         cl02(l) =-(ll - 2)*l*cl01(l)
         cl03(l) = 2*l*cl01(l)
         cl04(l) =-2*(ll - 2)*cl01(l)

         cl01(l) = 4*cl01(l)

      end do

      deallocate (cinv)

c   Matrix elements for particle conservation. These are the coefficients
c   of Vf(l=0) (l=0,1 contributions), and are related by -transpose to V(l=0) f
c   Recall that -efield (Bsubv d/du - Bsubu d/dv) is the electric drift term

      diagl(:,:,1)  =-2*omgl(l1) * diagl(:,:,1)
      diagl(:,:,2)  = 2*omgl(l1) * diagl(:,:,2)

      end subroutine lcalc


      subroutine printout(fz1p, fz1m, fz3p, fz3m, srces)

c  This subroutine calculates the diffusion coefficients, parallel
c  viscous stress, and banana-plateau flux, and writes them to DKESOUT.
c
c      fz1p:      Distribution function in response to source type = 1,+
c      fz1m:      Distribution function in response to source type = 1,-
c      fz3p:      Distribution function in response to source type = 3,+
c      fz3m:      Distribution function in response to source type = 3,-
c
c      Uses Eqs (24) and (33) in W. I. van Rij, et. al. paper
c
c      See subroutine RESIDUE for definitions of residuals
c
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      use Vimatrix
      use Vnamecl2
      use dkes_input, only: ipmb, lalpha, lscreen
      use dkes_realspace, only: bstrs, mpnt, mpnt2, mpnt3, mpnt4,
     1    diagl
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(mpnt,0:lalpha) :: fz1p, fz1m, fz3p, fz3m
      real(rprec), dimension(mpnt,4,2,2) :: srces
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: l0 = 1, l1 = 2, l2 = 3, l3 = 4,
     1   pgrad = 1, epar = 2, plus = 1, minus = 2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: mn
      real(rprec) :: L11p, L33p, rs11p, rs33p, L13p, L31p,
     1   rs13p, rs31p,
     1   stress1, bpflux1, L11m, L33m, rs11m, rs33m,
     2   L13m, L31m, rs13m, rs31m, g31max, g31min,
     3   del, av, scal11, scal33, scal13, srcp, srcm, t1, t2, t3, t4
      real(rprec), external :: DOT_G
C-----------------------------------------------------------------------

      srcp = DOT_G(mpnt3,fz1p(1,l0),1,srces(1,l0,1,plus),1)

      t1 = -DOT_G(mpnt4,fz1m(1,l0),1,srces(1,l0,pgrad,minus),1)/vp
      t2 =  DOT_G(mpnt,fz1p(1,l0),1,srces(1,l0,1,plus),1)/vp
      t3 =  -s1cs1
      t4 = sum(fz1m(:,0)*srces(:,l0,1,1))/vp

      srcm = t1 + t3 + t4

      if (lscreen) then
      print *,'  {F+,sigma+}     = ', srcp/vp
      print *,' -{F-,sigma-}+... = ', srcm

      print *,' -{F-,s-}    = ', t1
      print *,'  {F+,s+}    = ', t2,' {F0+,s+} = ',   t4
      print *,'  {s+,C-1s+} = ', t3
      end if

      if (ipmb .eq. 2) then
         L11p = 0;         L33p = 0;         L31p = 0;        L13p = 0
         rs11p = 0;        rs33p = 0;        rs13p = 0;       rs31p = 0
         rsd1p = 0;        rsd3p = 0;        crs1p = 0;       crs3p = 0
         stress1 = 0;      bpflux1 = 0
      else
c   See expressions following Eq.(24) in W. I. van Rij, et. al.
         L11p = sum(fz1p(:,l0)*srces(:,l0,pgrad,plus) +
     1              fz1p(:,l2)*srces(:,l2,pgrad,plus))/vp                !!{f1+, sigma1+}
         rs11p = abs(g11p/L11p)
         L11p = L11p - g11p                                              !!g11p = {f1+,(sigma1+ -Wf1+)} => 0
         L33p = sum(fz3p(:,l1)*srces(:,l1,epar,plus) +
     1              fz3p(:,l2)*srces(:,l2,epar,plus))/vp                 !!{f3+, sigma3+}
         rs33p = abs(g33p/L33p)
         L33p = L33p - g33p
         L13p = sum(fz1p(:,l1)*srces(:,l1,epar,plus) +
     1              fz1p(:,l2)*srces(:,l2,epar,plus))/vp                 !!{f1+, sigma3+}
         rs13p = abs(g13p/L13p)
         L13p  = L13p - g13p
         L31p = sum(fz3p(:,l0)*srces(:,l0,pgrad,plus) +
     1              fz3p(:,l2)*srces(:,l2,pgrad,plus))/vp                !!{f3+, sigma1+}
         rs31p = abs(g31p/L31p)
         L31p  = L31p - g31p
         stress1 = sum((fz1p(:,l2)+fz1m(:,l2))*bstrs)                    !!l=2 component of distribution
         bpflux1 = bpfac*stress1
      endif

      if (ipmb .eq. 1) then
         L11m = 0;     L33m = 0;      L13m = 0;     L31m = 0
         rs11m = 0;    rs33m = 0;     rs13m = 0;    rs31m = 0
         rsd1m  = 0;   rsd3m  = 0;    crs1m  = 0;   crs3m  = 0
      else
         L11m = sum(fz1m(:,0)*srces(:,l0,pgrad,plus))/vp                 !!{f01+, sigma1+}
         do mn = l1, l3
            L11m = L11m - sum(fz1m(:,mn)*srces(:,mn,pgrad,minus))/vp     !!-{f1-, sigma1-}
         end do
         rs11m = abs(g11m/L11m)
         L11m  = L11m + g11m - s1cs1                                     !!-{sigma1+, C-1(sigma1+)}

         L33m  =-sum(fz3m(:,l1)*srces(:,l1,epar,minus))/vp               !!-{f1-, sigma1-}
         rs33m = abs(g33m/L33m)
         L33m  = L33m + g33m + g33s                                      !!+{Fs, sigma_s}

         L13m  = sum(fz3m(:,0)*srces(:,l0,pgrad,plus))/vp
         do mn = l1, l3
            L13m = L13m - sum(fz3m(:,mn)*srces(:,mn,pgrad,minus))/vp
         end do
         rs13m = abs(g13m/L13m)
         L13m  = L13m + g13m

         L31m  =-sum(fz1m(:,l1)*srces(:,l1,epar,minus))/vp
         rs31m = abs(g31m/L31m)
         L31m  = L31m + g31m
      endif

      if (ipmb .ne. 0) then
         g31max = 0
         g31min = 0
      else                                                               !!Eq. (26)
         del = .5_dp*sqrt(abs((L11m - L11p)*(L33m - L33p)))
         av = .25_dp*(L13p + L31p + L13m + L31m)
         g31max = av + del
         g31min = av - del
      endif

c  scale factors relating DIJ in Eq. (36) of W. I. van Rij, et. al.
c  with output from mono-energetic code (these results for LIJ). Note the factor
c  of 0.5 arises due to the conversion from v (velocity) to K (kinetic energy)

      scal11 = 0.5_dp*(b00*(vthermi/wcyclo))**2 * vthermi
      scal33 = 0.5_dp*vthermi
      scal13 = 0.5_dp*(b00*(vthermi/wcyclo)) * vthermi

c  output results summary

      write (ioout, 10) '+', L11p, L33p, L13p, L31p, g31min,
     1                 '-', L11m, L33m, L13m, L31m, g31max,
     2                 scal11, scal33, scal13,
     2                 stress1, bpflux1, g33s,
     3                 'F(1,-)', rsd1m, crs1m, rs11m, rs13m,
     4                 'F(1,+)', rsd1p, crs1p, rs11p, rs13p,
     5                 'F(3,-)', rsd3m, crs3m, rs33m, rs31m,
     6                 'F(3,+)', rsd3p, crs3p, rs33p, rs31p

 10   format(/1x,'NEOCLASSICAL TRANSPORT MATRIX ELEMENTS: ',
     1       /1x,'DIJ(Eq.36,K=1) = LIJ * [(MKS FACT) * (Ti SCALE)]',//,
     1      6x,'PARITY',10x,'L11',12x,'L33',12x,'L13',12x,'L31',7x,
     2      'L13(min/max)',/,1x,131('-'),/,2(9x,a,3x,5(3x,1pe12.4),/),
     3      4x,'MKS FACT.',3(3x,1pe12.4),/,
     4      4x,'Ti SCALE',8x,'TI**1.5',8x,'TI**0.5',8x,'TI**1.0',
     3      //,4x,'PARALLEL',7x,'BAN-PLAT',10x,'L33'/,
     4      6x,'STRESS',9x,'FLUX',9x,'(SPITZER)',/,
     5      1pe12.4,4x,1pe12.4,3x,1pe12.4///,1x,
     5      'EQUATION RESIDUALS (L = KINETIC EQUATION OPERATOR)'//,
     6      3x,' F(I,+/-)',7x,'KIN. EQ.',4x,'PART. CONSERV.',3x,
     7      'RES(FI,LI)',5x,'RES(FJ,LI)',/,
     8      18x,'{L[FI]**2}',20x,'{FI,L[FI]}',5x,'{FJ,L[FI]}',/,
     9      1x,131('-'),/,4(5x,a,2x,4(3x,1pe12.4)/)/)

!
!     WRITE SUMMARY OPT_FILE FOR USE BY OPTIMIZER
!
      write(ioout_opt,'(3(2x,e24.13))') L11p, L33p, L31p
      write(ioout_opt,'(3(2x,e24.13))') L11m, L33m, L31m
      write(ioout_opt,'(3(2x,e24.13))') scal11, scal33, scal13

      close(unit=ioout_opt)

      end subroutine printout


      subroutine rddisk(iunit, a, now, incnow, irec, ierr)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer iunit, now, incnow, irec, ierr
      real(rprec), dimension(now) :: a
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ig, irech, il
C-----------------------------------------------
c
c     Author:    H. Maassberg       Sept. 1988
c
c     auxilary routine for routines BLKTRD / BLK5D (DKES code)
c     for disk i/o on CRAY at Garching
c
c-----------------------------------------------------------------------
c
c     read NOW words of vector A from disk (fortran IUNIT) with
c     direct access in record IREC
c
c-----------------------------------------------------------------------
      ig = 0
      irech = irec - 1
   
      do while (ig < now)
         il = ig + 1
         irech = irech + 1
         ig = min(now,ig + incnow)
         read (iunit, rec=irech, err=10) a(il:ig)
      end do
      ierr = 0
      return
   10 continue
      ierr = 1
      print *, ' error detected in disk i/o  (routine RDDISK)'

      end subroutine rddisk


      subroutine residue(a, bm1, bp1, cm2, cp2, csave, fz1, fz3,
     1   srces, rsd1, rsd3, g11, g33, g31, g13, crs1, crs3)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      use Vnamecl2
      use dkes_input, only: idisk, lalpha
      use dkes_realspace, only: diagl, diagle, srces0, mpnt
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), intent(out) ::
     1     rsd1, rsd3, g11, g33, g31, g13, crs1, crs3
      real(rprec), dimension(mpnt,0:lalpha), intent(in) :: fz1, fz3
      real(rprec), dimension(mpnt,mpnt) :: a, bm1, bp1, cm2, cp2, csave
      real(rprec), intent(in), dimension(mpnt,4,2,2) :: srces
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: k, kkount, mnp, j, ll
      real(rprec) :: fnrm1, fnrm3, fac, fnx1, fnx3, gnx1
      real(rprec), dimension(:), allocatable :: src1, src3, xrc1, xrc3
C-----------------------------------------------------------------------
!
!  This subroutine calculates the solution residuals.
!  See Eq.(24) and following in W. I. van Rij, et. al.
!
!  cm2 * f(l-2) + bm1 * f(l-1) + a * f(l) + bp1 * f(l+1) + cp2 * f(l+2) = source(l)
!
!  where:   bp1(l) = bm1(l+1) (transpose);     cp2(l) = cm1(l+2) (transpose)
!
!
!     read in f for tri-diagonal code 
!
!     read (33, iostat=k) fz1(1:mpnt,1:), fz3(1:mpnt,1:)
!     read (33, iostat=k) fz1(1:mpnt,0), fz3(1:mpnt,0)
!     if (k .ne. 0) stop 'Unable to read in data from tri-diag unit33'

      allocate (src1(mpnt), src3(mpnt), xrc1(mpnt), xrc3(mpnt), stat=k)
      if (k .ne. 0) stop 'allocation error in DKES residue!'

      fnrm1 = 0;  fnrm3 = 0;  rsd1 = 0;  rsd3 = 0
      g11 = 0;    g33 = 0;    g31 = 0;   g13 = 0

      k = idisk*lam6 + 1
      ll = lap1 - k
      kkount = lam3 - k

c  l >= 2

      call blox (ll, a, bm1, cm2, src1, src3, srces0)
      if (idisk .eq. 0) then                                !!ll = lalpha block row
         fnrm1 = fnrm1 + sum(src1*src1)
         fnrm3 = fnrm3 + sum(src3*src3)
         do mnp = 1, mpnt
            src1(:) = src1(:) - a(:,mnp)*fz1(mnp,ll) - bm1(:,mnp)
     1         *fz1(mnp,ll-1) - cm2(:,mnp)*fz1(mnp,ll-2)
            src3(:) = src3(:) - a(:,mnp)*fz3(mnp,ll) - bm1(:,mnp)
     1         *fz3(mnp,ll-1) - cm2(:,mnp)*fz3(mnp,ll-2)
         end do
         rsd1 = rsd1 + sum(src1*src1)
         rsd3 = rsd3 + sum(src3*src3)
         g11 = g11 - sum(src1*fz1(:,ll))
         g33 = g33 - sum(src3*fz3(:,ll))
         g31 = g31 - sum(src1*fz3(:,ll))
         g13 = g13 - sum(src3*fz1(:,ll))
      endif

      bp1 = bm1
      cp2 = cm2

      k = k + 1
      ll = ll - 1
      call blox (ll, a, bm1, cm2, src1, src3, srces0)
      if (idisk .eq. 0) then
         fnrm1 = fnrm1 + sum(src1*src1)
         fnrm3 = fnrm3 + sum(src3*src3)
         do mnp = 1, mpnt
            src1(:) = src1(:) - bp1(mnp,:)*fz1(mnp,ll+1) -
     1         a(:,mnp)*fz1(mnp,ll) - bm1(:,mnp)*fz1(mnp,ll-1) -
     2         cm2(:,mnp)*fz1(mnp,ll-2)
            src3(:) = src3(:) - bp1(mnp,:)*fz3(mnp,ll+1) -
     1         a(:,mnp)*fz3(mnp,ll) - bm1(:,mnp)*fz3(mnp,ll-1) -
     2         cm2(:,mnp)*fz3(mnp,ll-2)
         end do
         rsd1 = rsd1 + sum(src1*src1)
         rsd3 = rsd3 + sum(src3*src3)
         g11 = g11 - sum(src1*fz1(:,ll))
         g33 = g33 - sum(src3*fz3(:,ll))
         g31 = g31 - sum(src1*fz3(:,ll))
         g13 = g13 - sum(src3*fz1(:,ll))
      endif

      bp1 = bm1
      csave = cm2
      
      do j = 1, kkount
         k = k + 1
         ll = ll -1 
         call blox (ll, a, bm1, cm2, src1, src3, srces0)
         fnrm1 = fnrm1 + sum(src1*src1)
         fnrm3 = fnrm3 + sum(src3*src3)
         do mnp = 1, mpnt
            src1(:) = src1(:) - cp2(mnp,:)*fz1(mnp,ll+2) -
     1         bp1(mnp,:)*fz1(mnp,ll+1) - a(:,mnp)*fz1(mnp,ll) -
     2         bm1(:,mnp)*fz1(mnp,ll-1) - cm2(:,mnp)*fz1(mnp,ll-2)
            src3(:) = src3(:) - cp2(mnp,:)*fz3(mnp,ll+2) -
     1         bp1(mnp,:)*fz3(mnp,ll+1) - a(:,mnp)*fz3(mnp,ll) -
     2         bm1(:,mnp)*fz3(mnp,ll-1) - cm2(:,mnp)*fz3(mnp,ll-2)
         end do
         rsd1 = rsd1 + sum(src1*src1)
         rsd3 = rsd3 + sum(src3*src3)
         g11 = g11 - sum(src1*fz1(:,ll))
         g33 = g33 - sum(src3*fz3(:,ll))
         g31 = g31 - sum(src1*fz3(:,ll))
         g13 = g13 - sum(src3*fz1(:,ll))
         bp1 = bm1
         cp2 = csave
         csave = cm2
      end do

c  l = 1

      k = k + 1
      ll = ll - 1

      call blox (ll, a, bm1, cm2, src1, src3, srces0)
      fnrm1 = fnrm1 + sum(src1*src1)
      fnrm3 = fnrm3 + sum(src3*src3)

      do mnp = 1, mpnt
         src1(:) = src1(:) - diagl(:,mnp,iswpm)*fz1(mnp,0) -
     1      cp2(mnp,:)*fz1(mnp,ll+2) -
     2      bp1(mnp,:)*fz1(mnp,ll+1) - a(:,mnp)*fz1(mnp,ll) -
     3      bm1(:,mnp)*fz1(mnp,ll-1)
         src3(:) = src3(:) - diagl(:,mnp,iswpm)*fz3(mnp,0) -
     1      cp2(mnp,:)*fz3(mnp,ll+2) -
     2      bp1(mnp,:)*fz3(mnp,ll+1) - a(:,mnp)*fz3(mnp,ll) -
     3      bm1(:,mnp)*fz3(mnp,ll-1)
      end do
      rsd1 = rsd1 + sum(src1*src1)
      rsd3 = rsd3 + sum(src3*src3)
      g11 = g11 - sum(src1*fz1(:,ll))
      g33 = g33 - sum(src3*fz3(:,ll))
      g31 = g31 - sum(src1*fz3(:,ll))
      g13 = g13 - sum(src3*fz1(:,ll))
      bp1 = bm1
      cp2 = csave

c  l = 0

      k = k + 1
      ll = ll - 1
      call blox (ll, a, bm1, cm2, src1, src3, srces0)
      fnrm1 = fnrm1 + sum(src1*src1)
      fnrm3 = fnrm3 + sum(src3*src3)
      fac = iswpm - 1
      fnx1 = 0
      fnx3 = 0

      do mnp = 1, mpnt
         src1(:) = src1(:) - diagle(:,mnp,iswpm)*fz1(mnp,0)
     1           - cp2(mnp,:)*fz1(mnp,ll+2) - bp1(mnp,:)
     2           * fz1(mnp,ll+1) - a(:,mnp)*fz1(mnp,ll)
         src3(:) = src3(:) - diagle(:,mnp,iswpm)*fz3(mnp,0)
     1           - cp2(mnp,:)*fz3(mnp,ll+2) - bp1(mnp,:)
     2           *fz3(mnp,ll+1) - a(:,mnp)*fz3(mnp,ll)
      end do


c  Particle conservation. Here, The diagl, diagl coefficients (for V(l=0)) are
c  -transpose of the ones used in the Vf(l=0) contributions above.

      xrc1(:) = matmul(transpose(diagle(:,:,iswpm)),fz1(:,ll)) 
     1        + matmul(transpose(diagl(:,:,iswpm)), fz1(:,ll+1)) 
     2        + fac*srces(:,1,1,1)
      xrc3(:) = matmul(transpose(diagle(:,:,iswpm)),fz3(:,ll)) 
     1        + matmul(transpose(diagl(:,:,iswpm)), fz3(:,ll+1))
      fnx1 = sum((abs(matmul(transpose(diagle(:,:,1)),fz1(:,ll))) 
     1        + abs(matmul(transpose(diagl(:,:,1)),fz1(:,ll+1))))**2)
      fnx3 = sum((abs(matmul(transpose(diagle(:,:,1)),fz3(:,ll))) 
     1        + abs(matmul(transpose(diagl(:,:,1)),fz3(:,ll+1))))**2)

      crs1 = sum(xrc1*xrc1)
      crs3 = sum(xrc3*xrc3)
      gnx1 = fac*sum(srces(:,1,1,1)*srces(:,1,1,1))

      if (gnx1 .ne. zero) fnx1 = gnx1
      if (fnx1 .eq. zero) fnx1 = one
      if (fnx3 .eq. zero) fnx3 = one
      if (fnrm1 .eq. zero) fnrm1 = epsilon(fnrm1)
      if (fnrm3 .eq. zero) fnrm3 = epsilon(fnrm3)
      crs1 = sqrt(crs1/fnx1)
      crs3 = sqrt(crs3/fnx3)
      if (efield1.eq.zero .and. iswpm.eq.2) crs3 = zero
      rsd1 = sqrt((rsd1 + sum(src1*src1))/fnrm1)
      rsd3 = sqrt((rsd3 + sum(src3*src3))/fnrm3)
      g11 = (g11 - sum(src1*fz1(1,ll)) - sum(xrc1*fz1(:,0)))/vp
      g33 =  g33 - sum(src3*fz3(1,ll)) - sum(xrc3*fz3(:,0))/vp
      g31 = (g31 - sum(src1*fz3(1,ll)) - sum(xrc1*fz3(:,0)))/vp
      g13 = (g13 - sum(src3*fz1(1,ll)) - sum(xrc3*fz1(:,0)))/vp

      deallocate (src1, src3, xrc1, xrc3, stat=k)

      end subroutine residue


      subroutine wrdisk(iunit, a, now, incnow, irec, incb, ierr)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: iunit, now, incnow, irec, incb, ierr
      real(rprec), dimension(now) :: a
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ig, irech, il
C-----------------------------------------------
c
c-----------------------------------------------------------------------
c
c     Author:    H. Maassberg       Sept. 1988
c     Modified:  S. Hirshman        June, 2001
c
c     auxilary routine for routines BLKTRD / BLK5D (DKES code)
c     for disk i/o on CRAY at Garching
c
c-----------------------------------------------------------------------
c
c     write NOW words of vector A on disk (fortran IUNIT) with
c     direct access in record IREC
c
c-----------------------------------------------------------------------
      ig = 0
      irec = irec + incb
      irech = irec - 1
 
      do while (ig < now)
         il = ig + 1
         irech = irech + 1
         ig = min(now,ig + incnow)
         write (iunit, rec=irech, err=10) a(il:ig)
      end do
      ierr = 0
      return
   10 continue
      ierr = 1
      print *, ' error detected in disk i/o  (routine WRDISK)'

      end subroutine wrdisk


      subroutine wrout(fz1s, fz1c, fz3s, fz3c, srces)

c  This subroutine creates the file DKESF containing the magnetic field
c  and the final solution distributions.

C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      use Vnamecl2
      use safe_open_mod
      use dkes_input
      use dkes_realspace, only: mvalue, nvalue, mpnt
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(mpnt,0:lalpha) :: fz1s, fz1c, fz3s, fz3c
      real(rprec), dimension(mpnt,4,2,2) :: srces
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: idkes = 7
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: m, n, mn, ll, l, iodkes = idkes
C-----------------------------------------------

      call safe_open(iodkes, mn, 'DKESF', 'unknown', 'formatted')
      if (mn .ne. 0) stop 'ERROR OPENING DKESF FILE'

      write (iodkes, 100) chip, psip, nzperiod, cmul1, efield1, wtov,
     1   weov, mpolb, ntorb, lfout, mpnt, blabl(ibbi)
  100 format(/8x,'CHIP',8x,'PSIP',4x,'NZPERIOD',8x,'CMUL',6x,'EFIELD',4x
     1   ,'OMEGAT/V',4x,'OMEGAE/V',7x,'MPOLB',7x,'NTORB',7x,'LFOUT',2x,
     2   'BLOCK SIZE'//1p2e12.4,i12,4e12.4,4i12,2/4x,'M',4x,'N',14x,a3,
     3   18x,'BTHETA',7x,'BZETA'/)
      write (iodkes, 200) 0, nvalsb(1), borbi(1,1), btheta, bzeta
  200 format(2i5,5x,1pe12.4,12x,2e12.4)
      write (iodkes, 300) (m - 1,nvalsb(1),borbi(1,m),m=2,mpolb)
  300 format(2i5,5x,1pe12.4)
      do n = 2, ntorb
         write (iodkes, 300) (m - 1,nvalsb(n),borbi(n,m),m=1,mpolb)
      end do

      if (ipmb .eq. 1) then
         fz1c(:,:lfout) = zero
         fz3c(:,:lfout) = zero
      endif
      if (ipmb .eq. 2) then
         fz1s(:,:lfout) = zero
         fz3s(:,:lfout) = zero
      endif

      n = nzperiod
      if (n .eq. 0) n = 1

      write (iodkes, 400)
  400 format(/4x,'M',4x,'N',13x,'F11C',8x,'F11S',8x,'F13C',8x,'F13S'/)
      write (iodkes, 500) (mvalue(mn),nvalue(mn),fz1s(mn,0),fz1c(mn,0),
     1   fz3s(mn,0),fz3c(mn,0),mn=1,mpnt)
  500 format(2i5,5x,1p4e12.4)

      write (iodkes, 600)
  600 format(/4x,'M',4x,'N',4x,'L',8x,'F01C',8x,'F01S',8x,'F03C',8x,
     1   'F03S',20x,'S01C',8x,'S01S',8x,'S03C',8x,'S03S'/)

      l = -1
      do ll = 1, lfout
         l = l + 1
         if (ll > 4) write (iodkes, 700) (mvalue(mn),nvalue(mn),l,
     1      fz1c(mn,ll),fz1s(mn,ll),fz3c(mn,ll),fz3s(mn,ll),mn=1,mpnt)
  700    format(3i5,1p4e12.4)
         if (ll <= 4)
     1      write (iodkes, 800) (mvalue(mn),nvalue(mn),l,
     2      fz1c(mn,ll),fz1s(mn,ll),fz3c(mn,ll),fz3s(mn,ll),
     3      srces(mn,ll,1,2),srces(mn,ll,1,1),srces(mn,ll,2,2),
     4      srces(mn,ll,2,1),mn=1,mpnt)
  800    format(3i5,1p4e12.4,12x,4e12.4)
      end do

      close(iodkes)


      end subroutine wrout
