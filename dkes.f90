#!/bin/sh
#---------------------------------------------------------------------
cat > temp.c << "EOC"
cat > vmodules.f << "EOF"
      module Vimatrix
      integer :: ioout, ioout_opt
      end module Vimatrix

      module dkes_input
!
!     CONTAINS INPUT PARAMETERS AND ARRAYS READ IN FROM THE DATAIN NAMELIST
!
!     mpol      : Number of Fourier harmonics used to represent the theta
!                 dependence of the distribution functions; 1 < mpol <= mpold
!     mpolb     : Number of Fourier harmonics used to represent the theta
!                 dependence of the equilibrium magnetic field; 1 < mpolb <= mpolbd
!     ntor      : Number of Fourier harmonics used to represent the zeta
!                 dependence of the distributions; 0 < ntor <= ntord
!     ntorb     : Number of Fourier harmonics used to represent the zeta
!                 dependence of the equilibrium field; 0 < ntorb <= ntorbd
!     mmnn      : Matrix containing the distribution toroidal mode numbers
!                 (expressed in units of nzperiod) in column 1, and the
!                 number of poloidal modes associated with each toroidal
!                 mode in column 2. The poloidal mode numbers associated
!                 with each toroidal mode are stored in the same row,
!                 starting in column 3. WARNING : this spectrum must
!                 encompass the spectrum of the input magnetic field.

!     nrun      : Number of cmul = nu/v  and  efield = Erad/v values
!                 for which solutions are to be calculated; must be < krun
!     cmul      : Corresponding array of cmul values; should be non-zero.
!     efield    : Corresponding array of efield values.
!     nzperiod  : Number of toroidal field periods.
!     chip      : chip = d(chi)/d(rho), radial derivative of the poloidal flux.
!     psip      : psip = d(psi)/d(rho), radial derivative of the toroidal flux
!     btheta    : covariant poloidal component of B. In Boozer coordinates,
!                 a flux function (I, the TOROIDAL current density)
!     bzeta     : covariant toroidal component of B. In Boozer coordinates,
!                 a flux function (J, the POLOIDAL current density)
!     borbi     : Array of nonzero Fourier coefficients of |B| (or 1/|B|) defined by
!                 BORBI_REALSPACE (theta,zeta) = SUM borbi(n,m)*cos[m*theta - nzperiod*n*zeta],
!                 where m = 0 to mpolb, n = -ntorb to ntorb
!                 Note: borbi represents "|B|" or 1/|B| ("binverse"), depending on ibbi)
!                 OBSOLETE (see nvalsb array below):
!                 BORBI_REALSPACE(theta,zeta) = borbi(n,m)*cos[(m-1)*theta - nzperiod*nvalsb(n)*zeta].
!     ibbi  = 1 : borbi represents B (default, comes from boozmn file, for example)
!     ibbi  = 2 : borbi represents 1/B
!     nvalsb    : (NOW OBSOLETE) Array of ntorb toroidal mode numbers
!                 expressed in units of nzperiod. For an axisymmetric
!                 device, ntorb = 1 and nvalsb(1) = 0.

!     ipmb  = 0 : DKES computes both sine and cosine solutions (default)
!     ipmb  = 1 : DKES computes only sine            solutions.
!     ipmb  = 2 : DKES computes only          cosine solutions.

!     idisk = 0 : DKES computes all l values.
!     idisk = 1 : DKES computes only l = 0,1,2,3,4,5 (default)

!     lfout = 0 : Final distributions not written in the file DKESF. (default)
!     lfout > 0 : Final distributions     written in the file DKESF,
!                 for l+1 = 1,2,.....,lfout <= lalpha.

!     meshtz    : Determines theta and zeta meshes for Fourier transforms.
!                 See ntheta and nzeta in the mnset subroutine. Used to
!                 ensure accuracy of spectral matrix elements.


!     lalpha    : Number of orthonormalized Legendre polynomials used to
!                 represent the dependence of the distributions on the
!                 pitch-angle variable cos(alpha); lalpha >= 6



      use kind_spec
      integer, parameter :: mpold = 100
      integer, parameter :: ntord = 100
      integer, parameter :: mpolbd = 25
      integer, parameter :: ntorbd = 25
      integer, parameter :: bigint = 2**16
      integer, parameter :: krun = 100

      integer :: mpol, ntor, lalpha, ipmb, meshtz, idisk,
     1   lfout, ifscl, nrun, nzperiod, mpolb, ntorb, ibbi

      integer, dimension(mpold + 2,ntord) :: mmnn

      real(rprec), dimension(krun) :: cmul, efield
      real(rprec) :: chip, psip, btheta, bzeta
      real(rprec), dimension(-ntorbd:ntorbd,0:mpolbd) :: borbi
      integer, dimension(ntorbd) :: nvalsb                              !OBSOLETE
      character*(131) :: dashes
      logical :: lscreen
      
      namelist /datain/mpol, ntor, lalpha, ipmb, mmnn, meshtz, idisk,
     1   lfout, ifscl, nrun, cmul, efield, nzperiod, chip, psip, mpolb,
     2   ntorb, nvalsb, borbi, btheta, bzeta, ibbi

      contains

      subroutine read_dkes_input
      use safe_open_mod
      use Vimatrix, only: ioout, ioout_opt
      implicit none
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: idata = 7
      integer, parameter :: iout = 20
      integer, parameter :: iout_opt = 14
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: numargs, icount, istat, index_blank, iodata
      character*120 :: output_file, opt_file, input_file
      character*120 :: arg1(5)
      logical :: lexist
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      external getcarg, dkes_input_prepare
C-----------------------------------------------
      lscreen = .true.

c  read input data from "datain" namelist
      call getcarg(1, arg1(1), numargs)

      do icount = 2, numargs
         call getcarg(icount, arg1(icount), istat)
      end do
      if (numargs .ge. 4) then
         if (numargs .eq. 4) arg1(5) = 'T'
         call dkes_input_prepare (arg1, input_file)
         if (arg1(5)(1:1).eq.'f' .or. arg1(5)(1:1).eq.'F') 
     1      lscreen = .false.
      else if (numargs .eq. 1) then
         input_file = trim(arg1(1))
      else
         print *, ' MUST ENTER INPUT FILENAME ON COMMAND LINE'
         stop
      end if

      index_blank = index(input_file,'input.')

      if (index_blank .eq. 0) then
         inquire (file=trim(input_file), exist=lexist)        
         if (.not.lexist) then
            input_file = 'input.' // trim(input_file)
            index_blank = 1
         else
            index_blank = index(input_file,'.') - 5
         end if
      end if

      index_blank = max(index_blank + 5,1)
      output_file= 'dkesout' // input_file(index_blank:)
      opt_file= 'opt_dkes' // input_file(index_blank:)             !DAS 2/21/2000

!
!     OPEN INPUT AND OUTPUT FILES FOR READING (AND WRITING OUTPUT)
!
      iodata = idata
      call safe_open(iodata, istat, input_file, 'old', 'formatted')
      if (istat .ne. 0) stop 'Error reading input file in DKES'
      ioout = iout
      call safe_open(ioout, istat, output_file, 'replace', 'formatted')
      if (istat .ne. 0) stop 'Error writing output file'
      ioout_opt = iout_opt
      call safe_open(ioout_opt, istat, opt_file, 'replace','formatted')
      if (istat .ne. 0) stop 'Error writing opt_output file'

!     Read namelist (datain) input
      nvalsb(1) = -bigint-1
      idisk = 1
      lfout = 0
      ibbi = 1
      read (iodata, datain)
      close (iodata)

!     perform error checking on input data

      if (nzperiod .le. 0) nzperiod = 1
      if (ibbi<1 .or. ibbi >2) then
         write (ioout,'(a)') ' ibbi must =1 or =2'
         stop ' ibbi <1 or ibbi >2 in DKES input'
      end if

      if (mpol < 1) then
         write (ioout, 10) mpol
  10     format(' mpol = ',i5,'  is less than 1')
         stop ' mpol < 1 in DKES input'
      endif

      if (mpol > mpold) then
         write (ioout, 15) mpol, mpold
  15     format(' mpol = ',i5,'  is greater than mpold = ',i5)
         stop ' mpol > mpold in DKES input'
      endif

      if (mpolb < 2) then
         write (ioout, 20) mpolb
  20     format(' mpolb = ',i5,'  is less than 2')
         stop ' mpolb < 2 in DKES input'
      endif

      if (mpolb > mpolbd) then
         write (ioout, 25) mpolb, mpolbd
  25     format(' mpolb = ',i5,'  is greater than mpolbd = ',i5)
         stop ' mpolb > mpolbd in DKES input'
      endif

      if (ntor < 1) then
         write (ioout, 30) ntor
  30     format(' ntor = ',i5,'  is less than 1')
         stop ' ntor < 1 in DKES input'
      endif

      if (ntor > ntord) then
         write (ioout, 35) ntor, ntord
  35     format(' ntor = ',i5,'  is greater than ntord = ',i5)
         stop ' ntor > ntord in DKES input'
      endif

      if (ntorb < 1) then
         write (ioout, 40) ntorb
  40     format(' ntorb = ',i5,'  is less than 1')
         stop ' ntorb < 1 in DKES input'
      endif

      if (ntorb > ntorbd) then
         write (ioout, 45) ntorb, ntorbd
  45     format(' ntorb = ',i5,'  is greater than ntorbd = ',i5)
         stop ' ntorb > ntorbd in DKES input'
      endif

      if (lalpha < 6) then
         write (ioout, 50) lalpha
  50     format(' lalpha = ',i5,'  is less than 6')
         stop ' lalpha < 6 in DKES input'
      endif

      if (ipmb<0 .or. ipmb>2) ipmb = 0

      meshtz = max(0,meshtz)

      do istat = 1, len(dashes)
         dashes(istat:istat) = "-"
      end do

      end subroutine read_dkes_input

      end module dkes_input


      module Vcmain
      use kind_spec
      integer, parameter :: lsource = 4
      real(rprec), dimension(:), allocatable, target :: fzerop, fzerom
      end module Vcmain

      module Vnamecl2
      use Vcmain, only: lsource, rprec, dp
      integer :: lap1, lam1, lam3, lam6, iswpm, ier
      real(rprec), parameter :: zero = 0, one = 1, two = 2
      character*3, dimension(2), parameter :: blabl = (/'  B','1/B'/)

      real(rprec), save :: p2, root2, rthalf, rt3o2, rt5o2
      real(rprec) :: b00, bsqav
      real(rprec), dimension(:), allocatable ::
     1   cols, al1, al2, al3, al4, bl1,
     2   bl2, bl3, bl4, cl1, cl2, cl3, cl4
      real(rprec), dimension(:), allocatable ::
     1   cols0, omgl, al01, al02, al03, al04, bl01,
     2   bl02, bl03, bl04, cl01, cl02, cl03, cl04
      real(rprec) :: cmul1, efield1, wtov, weov, wcyclo, vthermi, vp,
     1   bpfac, s1cs10, s1cs1, 
     2   rsd1p, rsd1m, rsd3p, rsd3m, g11p, g11m, g33p, g33m, g31p,
     3   g33s, g31m, g13p, g13m, crs1p, crs3p, crs1m, crs3m

      end module Vnamecl2

      module dkes_realspace
      use kind_spec
      integer, dimension(:), allocatable :: mvalue, nvalue
      real(rprec), dimension(:), allocatable :: exbgrad, bgrad, bstrs,
     1   srces0, borbi1, auxs3m
      real(rprec), dimension(:,:), allocatable :: auxs1, auxs3p
      real(rprec), dimension(:,:,:), allocatable :: diagl, diagle
      real(rprec), dimension(:,:,:), allocatable :: bmat2, bmat3,
     1   bmat4, bmat5, bmat6, matjac
      real(rprec), dimension(:,:), allocatable :: bmat1
      real(rprec), allocatable, dimension(:) :: blk1, blk2, blk3, blk4,
     1   blk5, blk6, blk7
      integer :: mn0, mpnt, mpntp1, mpnt2, mpnt3, mpnt4,
     1   mpnt5, mpntsq, mmax, nmax, ntheta, nzeta

      contains

      subroutine set_mndim
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      use Vimatrix
      use Vnamecl2
      use dkes_input
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: mn, mtop, ntop, nn, mbig, mm, istat, m1, n1
C-----------------------------------------------
c  This subroutine sets up Fourier mode number arrays and block dimensions
c  and allocates memory for spatial block arrays

!
!     Constants for Legendre polynomial norms
!
      p2 = one/sqrt(20._dp)
      root2 = sqrt(two);          rthalf = one/root2
      rt3o2 = sqrt(3._dp/2._dp);  rt5o2 = sqrt(5._dp/2._dp)
      
c  stack Fourier modes [(m=0, n<0) components are discarded]
c  only m < mpol modes retained. only ntor toroidal modes retained.

      mn = 0

c  first compute mpnt from mmnn array in namelist

      ntor = max(1, abs(ntor))
      ntor = min(ntor, ntord)
      do nn = 1, ntor
         mbig = 2 + abs(mmnn(2,nn))
         do mm = 3, mbig
            m1 = mmnn(mm,nn)
            n1 = mmnn(1,nn)
            if ((m1 .lt. mpol) .and. (m1.ne.0 .or. n1.ge.0) .and.
     1          all(m1 .ne. mmnn(3:mm-1,nn))) mn = mn + 1
         end do
      end do

      mpnt = mn

      allocate (mvalue(mpnt), nvalue(mpnt), stat=istat)
      if (istat .ne. 0) stop 'Allocation error(1) in DKES2'

      mvalue = -100000;   nvalue = -100000
      mn = 0
      mtop = 0
      ntop = 0
      do nn = 1, ntor
         mbig = 2 + abs(mmnn(2,nn))
         do mm = 3, mbig
            m1 = mmnn(mm,nn)
            n1 = mmnn(1,nn)
            if ((m1.lt.mpol) .and. (m1.ne.0 .or. n1.ge.0)) then
               if (any(m1.eq.mvalue .and. n1.eq.nvalue)) cycle           !!do not read same m,n more than once
               mn = mn + 1
               mvalue(mn) = m1
               nvalue(mn) = n1                                           !!Per/period n
               mtop = max(mtop, abs(m1))
               ntop = max(ntop, abs(n1))
               if (m1.eq.0 .and. n1.eq.0) mn0 = mn
            endif
         end do
      end do

      if (mn .ne. mpnt)
     1   stop 'Error counting mmnn modes in DKES set_mndim'

      mpnt2 = 2*mpnt
      mpnt3 = 3*mpnt
      mpnt4 = 4*mpnt
      mpnt5 = 5*mpnt
      mpntp1 = mpnt + 1
      mpntsq = mpnt*mpnt
      mmax = 2*mtop
      nmax = 2*ntop

      allocate (exbgrad(mpnt), bgrad(mpnt), bstrs(mpnt), auxs3p(mpnt,2),
     1         auxs1(mpnt,3), auxs3m(mpnt), diagl(mpnt,mpnt,2),
     2         diagle(mpnt,mpnt, 2), srces0(16*mpnt),
     3         blk1(mpntsq), blk2(mpntsq), blk3(mpntsq),
     4         blk4(mpntsq), blk5(mpntsq), blk6(mpntsq),
     5         blk7(mpntsq), borbi1(mpnt),
     6         bmat2(mpnt,mpnt,2), bmat4(mpnt,mpnt,2), bmat1(mpnt,mpnt),
     7         bmat3(mpnt,mpnt,2), bmat5(mpnt,mpnt,2), 
     8         bmat6(mpnt,mpnt,2), matjac(mpnt,mpnt,2), stat = istat)
      if (istat .ne. 0) stop 'allocation error(1) in set_mndim!'


      end subroutine set_mndim

      subroutine free_mndim

      deallocate (exbgrad, bgrad, bstrs, borbi1,
     1   auxs1, auxs3p, auxs3m, diagl, diagle, 
     2   srces0, blk1, blk2, blk3, blk4, blk5, blk6, blk7,
     3   bmat2, bmat4, bmat1, matjac, bmat3, bmat5, bmat6)
      deallocate (mvalue, nvalue)

      end subroutine free_mndim

      subroutine ftconv
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      use Vnamecl2
      use dkes_input
      implicit none
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: plus = 1, minus = 2
      real(rprec), parameter :: specac = 1.e-12_dp
      real(rprec), parameter :: pt1 = 0.1_dp, half = 0.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, dimension(:), allocatable :: nval_bfld, mval_bfld
      integer :: nznt, mn, n, m, nn, mm, k, j, l, mnp, mp, np, 
     1          ms, ns, md, nd, is, imask(1)
      real(rprec), dimension(:), allocatable :: trigs, trigc, 
     1    rmn, bmn, qmn, q1mn, q2mn, dnorm
      real(rprec), dimension(:,:,:), allocatable :: blank
      real(rprec), dimension(:,:), allocatable :: overb1, bsqi1, 
     1   trigtc, trigts, trigzc, trigzs
      real(rprec) :: dnrm0, dnrm, twopi, dangt, dangz, ang,
     1   bi0, bb1, bi1, bisq1, fac, bi20, btc0, bzc0, facjac,
     2   vpnz, jac, bgradb, blebi, qs, bi,
     3   jaci, gradbgi, bgbi2gi, bigi, bgbibigi, blbi, qbi, qlbi,
     4   b1, trc, trs, facc, c2, c3, stfac, bi2gi
C-----------------------------------------------
!     This subroutine calculates magnetic field spectral arrays (Fourier
!     transforms that are needed for convolutions).
!
!     overb1        :  2D array of Fourier coefficients of 1/|B|
!                      1/|B|  = sum [ overb1(m,n) * cos(mu - nv) ]
!
!     bsqi1         :  2D array of Fourier coefficients of 1/|B|*2
!                      1/|B|**2  = sum [ bsqi1(m,n) * cos(mu - nv) ]
!
!     bgrad         :  m*chip - n*phip (Fourier representation of B dot grad operator)
!
!     exbgrad       : (m*bzeta + n*btheta)/|B00|**2 (multiplier of -Erad X B drift term; multiplied
!                                                   later by -efield to get correct ovarall sign)
!
!     srces0(1:mpnt):  (m*bzeta + n*btheta)*(1/B**2)(m,n)/3 (Vdrift source for l=0)
!     srces0(mpnt2+1: 3*mpnt)
!                   :  (1/sqrt(20)*sreces0(1:mpnt)  (Vdrift source for l=2 Legendre)
!
!     bgbi2gi       :  (jacobian * B dot grad 1/|B|)**2 / jacobian
!                      (Note mneumonics: bgbi -> jacobian * B dot grad 1/|B|, gi -> 1/jacobian...)
!
!     bgbibigi      :  (jacobian * B dot grad (1/|B|)] / (jacobian |B|)  (bi -> 1/|B|...)
!
!     bi2gi         :  1/(jacobian*|B|**2)
!

!     Prepare m,n arrays for stacking magnetic field coefficients

      allocate (nval_bfld(-ntorb:ntorb), mval_bfld(0:mpolb))

      if (nvalsb(1) .gt. -bigint) then
         if (lscreen) then
         print *,' This is an old-style input file.'
         print *,' Convert to new-style by eliminating nvalsb array'
         print *,' and then index borbi array with actual (n,m) indices'
         print *,' Note: n is in units of field period'
         end if
         nval_bfld = -bigint
         do l = 1, ntorb
            nval_bfld(l) = nvalsb(l)
         end do
         mval_bfld(0) = -bigint
         do l = 1, mpolb
            mval_bfld(l) = l - 1
         end do
      else
         do l = -ntorb, ntorb
            nval_bfld(l) = l
         end do
         do l = 0, mpolb
            mval_bfld(l) = l
         end do
      end if


      mmax = max(mmax, maxval(mval_bfld))
      nmax = max(nmax, maxval(nval_bfld))
      imask = minloc(nval_bfld, nval_bfld > -bigint) - 1 - ntorb
      nn  = min(-nmax, nval_bfld(imask(1)))
      nmax = max(nmax,abs(nn))

      ntheta = (3 + meshtz)*(mmax/2 + 1) - 2 - meshtz
      nzeta  = (3 + meshtz)*(nmax/2 + 1) - 2 - meshtz

      allocate (trigtc(0:mmax,ntheta), trigts(0:mmax,ntheta),
     1          trigzc(-nmax:nmax,nzeta), trigzs(-nmax:nmax,nzeta),
     2          blank(0:mmax,-nmax:nmax,7), overb1(0:mmax,-nmax:nmax),
     3          bsqi1(0:mmax,-nmax:nmax), rmn(mpnt), bmn(mpnt), 
     4          trigs(mpnt), trigc(mpnt), 
     5          qmn(mpnt), q1mn(mpnt), q2mn(mpnt), dnorm(mpnt), stat=j)
      if (j .ne. 0) stop 'Allocation error in ftconv!'


c  initialize to zero; Fourier transform normalization arrays

      nznt = ntheta*nzeta
      dnrm0 = sqrt(one/nznt)
      dnrm = root2*dnrm0

      srces0 = 0;      borbi1 = 0;      overb1 = 0;      bsqi1 = 0
      rmn    = 0;      bmn    = 0;      qmn    = 0
      q1mn   = 0;      q2mn   = 0
      s1cs10 = 0;      dnorm  = dnrm;   blank  = 0
      dnorm(mn0) = dnrm0


!     stack magnetic field Fourier coefficients into 1D array borbi1
      
      do nn = -ntorb, ntorb
         mloop: do mm = 0, mpolb
            do mn = 1, mpnt
               if (mval_bfld(mm).eq.mvalue(mn) .and.
     1             nval_bfld(nn).eq.nvalue(mn)) then
                  borbi1(mn) = borbi(nn,mm)
                  cycle mloop
               endif
            end do
         end do mloop
      end do


c  theta and zeta trig functions over one field period

      twopi = 8*atan(one)
      dangt = twopi/ntheta
      dangz = twopi/nzeta

      do k = 1, ntheta
         trigtc(0,k) = 1
         trigts(0,k) = 0
         do m = 1, mmax
            ang = dangt*mod(m*(k - 1), ntheta)
            trigtc(m,k) = cos(ang)
            trigts(m,k) = sin(ang)
         end do
      end do
      do j = 1, nzeta
         trigzc(0,j) = 1
         trigzs(0,j) = 0
         do n = 1, nmax
            ang = dangz*mod(n*(j - 1), nzeta)
            trigzc(n,j) = cos(ang)
            trigzs(n,j) = sin(ang)
            trigzc(-n,j) =  trigzc(n,j)
            trigzs(-n,j) = -trigzs(n,j)
         end do
      end do

!     form bi1 = 1/B and bisq1 = 1/B**2 cosine series in REAL space to compute
!     respective Fourier coefficients for 1/B (overb1) and 1/B**2 (bsqi1).
!     need for ALL m = 0, mmax, n= -nmax,nmax, not JUST mpnt space

      do k = 1, ntheta
         do j = 1, nzeta
            bb1 = 0
            do n = -ntorb, ntorb
               np = nval_bfld(n)
               if (np .le. -bigint) cycle
               do m = 0, mpolb
                  mp = mval_bfld(m)
                  if (mp .le. -bigint) cycle              
                  trc = trigtc(mp,k)*trigzc(np,j)
     1                + trigts(mp,k)*trigzs(np,j)
                  bb1 = bb1 + trc*borbi(n,m)
                end do
            end do
            if (bb1 .le. zero)
     1          stop ' 1/|B| <= 0: Check BORBI array input'
            bi1 = bb1
            if (ibbi .eq. 1) bi1 = one/bi1
            bisq1 = bi1*bi1
            do n = -nmax, nmax
               do m = 0, mmax
                  trc=trigtc(m,k)*trigzc(n,j)+trigts(m,k)*trigzs(n,j)
                  overb1(m,n) = overb1(m,n) + trc*bi1
                  bsqi1(m,n)  = bsqi1(m,n)  + trc*bisq1
               end do
            end do
         end do
      end do

      overb1(0,-nmax:-1) = 0
      bsqi1(0,-nmax:-1) = 0
      
      deallocate (nval_bfld, mval_bfld)

      overb1 = dnrm**2 * overb1
      bsqi1  = dnrm**2 * bsqi1
      overb1(0,0) = overb1(0,0)/2
      bsqi1(0,0)  = bsqi1(0,0)/2
      if (overb1(0,0) .eq. zero) then
         stop 'Error: |B|(0,0) = 0'
      else
         b00 = one/overb1(0,0)
      end if

      if (ibbi .eq. 2) then
         fac = specac*bsqi1(0,0)
         where (abs(bsqi1) < fac) bsqi1 = 0
         fac = specac*overb1(0,0)
         where (abs(overb1) < fac) overb1 = 0
      end if

!     compute source S1+(for l=0,2 Legendre pitch harmonics) sine series, Jacobian cosine series,
!     and parallel-streaming and electric (EXB-drift) matrix elements
!     NOTE: S1 source should be multiplied by B00*Rhoi*VTi (VTi=ion thermal speed),
!        and f1 (distribution) by B00*Rhoi to be in proper units. The implied
!        scaling of the transport matrix elements (to get into real units) is:
!
!        L11(real) = L11(DKES2) * (B00*Rhoi)**2 * VTi
!        L13(real) = L13(DKES2) * (B00*Rhoi) * VTi
!        L33(real) = L33(DKES2) * VTi
!

      bi0 = overb1(0,0)                                                  !m=0, n=0 1/|B| coefficient
      bi20 = bsqi1(0,0)                                                  !m=0, n=0 1/|B|**2 coefficient
      facjac = btheta*chip + bzeta*psip                                  !jacobian = facjac / |B|**2
      if (facjac .eq. zero) stop 'facjac = 0 in ftconv'
      bsqav = one/bi20                                                   !<B**2> in Boozer coordinates
      vp = facjac*bi20                                                   !<<jacobian>>, <<...>> u,v average
      vpnz = vp*nznt
      wtov = bi0*abs(chip/vp)                                            !transit freq, divided by 1/v
      if (chip .ne. zero) then
         bpfac = bzeta/(chip*bsqav)
      else 
         bpfac = bzeta/(epsilon(chip)*bsqav)
      end if
      
      do mn = 1, mpnt
         m = mvalue(mn);   n = nvalue(mn)
         srces0(mn) = (m*bzeta + nzperiod*n*btheta)
     1               * bsqi1(m,n)/3                                      !l=0 source (sine coefficients): S1
         srces0(mn+mpnt2) = p2*srces0(mn)                                !l=2 source ~ S1
         auxs3m(mn) = overb1(m,n)
         bgrad(mn) = m*chip - nzperiod*n*psip                            !jacobian*(B dot grad)
         exbgrad(mn) = bi20*(m*bzeta + nzperiod*n*btheta)                !Eradial X B term
      end do

      bi2gi = one/facjac

!     Compute required Fourier transforms

      THETA: do k = 1, ntheta
         ZETA: do j = 1, nzeta
            bi1 = 0;  jac = 0;  bgradb = 0;  blebi = 0; qs = 0
            do n = -nmax, nmax
               do m = 0, mmax
                  trc=trigtc(m,k)*trigzc(n,j)+trigts(m,k)*trigzs(n,j)
                  trs=trigts(m,k)*trigzc(n,j)-trigtc(m,k)*trigzs(n,j)
                  bi1 = bi1 + overb1(m,n)*trc                            !1/B(nu,nv)
                  jac = jac + bsqi1(m,n)*trc                             !jacobian(nu,nv)
                  bgradb = bgradb - 
     1                   (m*chip - nzperiod*n*psip)*overb1(m,n)*trs      !jacobian B dot grad (1/|B|)
                  trs = (m*bzeta + nzperiod*n*btheta)*trs
                  blebi = blebi - overb1(m,n)*trs                        !grad of 1/B in E X B direction
                  qs  = qs - bsqi1(m,n)*trs                              !-source (l=0)  
               end do
            end do

            jac = jac*facjac;  qs = qs/3;  blebi = blebi*bi20
            b1 = one/bi1                      !|B|(nu,nv)
            jaci = one/jac                    !1/jacobian(nu,nv)
            bigi = bi1*jaci                   !1/(jac*B)(nu,nv)

            gradbgi = bgradb*jaci             !B dot grad (1/|B|))
            bgbi2gi = bgradb*gradbgi          ![jacobian * B dot grad (1/|B|)]**2 / jacobian
            bgbibigi = bgradb*bigi            ![jacobian * B dot grad (1/|B|)] / (jacobian |B|)
            blbi = bgradb*b1                  !|B| (jac * B dot grad (1/|B|))

            qbi = qs*bigi                     !-source /(B*jac)
            qs = qs*jaci                      !-source /jac
            qlbi = qs*bgradb                  !-source * B dot grad (1/|B|)
            blebi = blebi*b1*b1


c  source terms {sigma+(i), C[-1] sigma+(j): s1cs1

            s1cs10 = s1cs10 - jac*qs**2       !-jac*(source/jac)**2

c  series for parallel stress, S3+(l=1,2), S1-(l=0,1,2,3), and
c  S3-(l=0,1,2,3)
 
            do mn = 1, mpnt
               n = nvalue(mn)
               m = mvalue(mn)
               trigc(mn) = trigtc(m,k)*trigzc(n,j) 
     1                    + trigts(m,k)*trigzs(n,j)
               trigs(mn) = trigts(m,k)*trigzc(n,j) 
     1                    - trigtc(m,k)*trigzs(n,j)
            end do

            rmn(:) = rmn(:) + trigs(:)*blbi
            bmn(:) = bmn(:) + trigs(:)*blebi
            qmn(:) = qmn(:) + trigs(:)*qs
            q1mn(:) = q1mn(:) + trigs(:)*qbi
            q2mn(:) = q2mn(:) + trigc(:)*qlbi
                        
c  cosine and sine integrals for matrix elements

            do n = -nmax, nmax
               do m = 0, mmax
                  trc=trigtc(m,k)*trigzc(n,j)+trigts(m,k)*trigzs(n,j)
                  trs=trigts(m,k)*trigzc(n,j)-trigtc(m,k)*trigzs(n,j)
                  blank(m,n,1) = blank(m,n,1) + trc*jac
                  blank(m,n,2) = blank(m,n,2) + trc*jaci
                  blank(m,n,3) = blank(m,n,3) + trc*bigi
                  blank(m,n,4) = blank(m,n,4) + trc*bgbi2gi
                  blank(m,n,5) = blank(m,n,5) + trs*gradbgi
                  blank(m,n,6) = blank(m,n,6) + trs*bgbibigi
                  blank(m,n,7) = blank(m,n,7) + trc*bi1
               end do
            end do

         end do ZETA
      end do THETA

      deallocate (trigs, trigc, trigtc, trigts, trigzc, trigzs)

      rmn = rmn*dnorm;   bmn = bmn*dnorm;   qmn = qmn*dnorm
      q1mn = q1mn*dnorm; q2mn = q2mn*dnorm

      dnorm = rthalf * dnorm

      l = 0
      do mnp = 1, mpnt
         mp = mvalue(mnp)
         np = nvalue(mnp)
         do mn = 1, mpnt
            facc = dnorm(mn)*dnorm(mnp)
            ms = mvalue(mn) + mp
            ns = nvalue(mn) + np
            md = mvalue(mn) - mp
            nd = nvalue(mn) - np
            is = 1
            if (md < 0) is = -1
            md = is*md
            nd = is*nd

c  matrix elements <ei M ej> where ei, ej are e+ = sin, e- = cos
c  symmetric sine-sine and cosine-cosine matrix elements

            if (mn .le. mnp) then
               matjac(mn,mnp,1) = facc*(blank(md,nd,1)-blank(ms,ns,1))   !<e+ M e+>, M = jacobian
               matjac(mn,mnp,2) = facc*(blank(md,nd,1)+blank(ms,ns,1))   !<e- M e->
               bmat2(mn,mnp,1)  = facc*(blank(md,nd,4)-blank(ms,ns,4))   !<e+ M e+>, M = (bgrad(1/B))**2/jacobian
               bmat2(mn,mnp,2)  = facc*(blank(md,nd,4)+blank(ms,ns,4))   !<e- M e->

               bmat3(mn,mnp,1) = facc*(blank(md,nd,2)+blank(ms,ns,2))    !<e- M e->, M = 1/jacobian
               bmat3(mn,mnp,2) = facc*(blank(md,nd,2)-blank(ms,ns,2))    !<e+ M e+>
               bmat5(mn,mnp,1) = facc*(blank(md,nd,3)+blank(ms,ns,3))    !<e- M e->, M = 1/(jacobian*B)
               bmat5(mn,mnp,2) = facc*(blank(md,nd,3)-blank(ms,ns,3))    !<e+ M e+>

               diagl(mn,mnp,1) = facc*(blank(md,nd,7)-blank(ms,ns,7))
               diagl(mn,mnp,2) = facc*(blank(md,nd,7)+blank(ms,ns,7))
            endif

c  sine-cosine matrix elements

            bmat4(mn,mnp,1) = facc*(is*blank(md,nd,6) + blank(ms,ns,6))
            bmat6(mn,mnp,1) = facc*(is*blank(md,nd,5) + blank(ms,ns,5))

         end do
      end do

      bmat4(:,:,2) = -transpose(bmat4(:,:,1))
      bmat6(:,:,2) = -transpose(bmat6(:,:,1))

      deallocate (blank)

c  complete integrals, sources, and fill out symmetric matrix elements
      c2 = one/(rt3o2**2*vpnz)
      s1cs10 = (root2*p2)**2*s1cs10/vpnz

      c2 = dnrm0/rt3o2                                                   !rt3o2 comes from norm of p(l=1) source
      auxs3m(:) = auxs3m(:)*facjac*rthalf/rt3o2                           !assumes jac = facjac/B**2
      auxs3m(mn0) = auxs3m(mn0)/rthalf

      bmat1 = 0
      do mn = 1, mpnt
         bmat1(mn,mn) = bgrad(mn)*bgrad(mn)*bi2gi
      end do

      do l = plus, minus
         do mn = 1, mpnt
            diagl(mn:mpnt,mn,l) = diagl(mn,mn:mpnt,l)
            bmat2(mn:mpnt,mn,l) = bmat2(mn,mn:mpnt,l)
            bmat3(mn:mpnt,mn,l) = bmat3(mn,mn:mpnt,l)
            bmat5(mn:mpnt,mn,l) = bmat5(mn,mn:mpnt,l)
            matjac(mn:mpnt,mn,l) = matjac(mn,mn:mpnt,l)
         end do

         do mnp = 1, mpnt
            bmat3(:,mnp,l) = exbgrad(:)*exbgrad(mnp)*bmat3(:,mnp,l)
            bmat4(:,mnp,l) = bgrad(mnp)*bmat4(:,mnp,l)
            bmat5(:,mnp,l) = bgrad(:)*exbgrad(mnp)*bmat5(:,mnp,l)
            bmat6(:,mnp,l) = exbgrad(mnp)*bmat6(:,mnp,l)
            diagl(:,mnp,l) = bgrad(mnp) * diagl(:,mnp,l)
         end do
      end do

      
      stfac = -dnrm0/(rt5o2*vp)
      c3 = p2*root2*dnrm0                   !p2 arises from l=2 component of source(+)
      
      bstrs(:)    = stfac*rmn(:)
!
!     note: these sources are all negative, since they get divided
!     in cescale by -1/nu(l) so the (-) signs cancel..
!
      auxs3p(:,1) = c2*bmn(:)
      auxs3p(:,2) = c2*rmn(:)
      auxs1(:,1)  = c3*(2*q1mn(:)*bgrad(:) + q2mn(:))
      auxs1(:,2)  = c3*qmn(:)*exbgrad(:)
      auxs1(:,3)  = 2*c3*(q1mn(:)*bgrad(:)-2*q2mn(:))


      deallocate (overb1, rmn, bmn, qmn, q1mn, q2mn, bsqi1, dnorm)

      end subroutine ftconv

      end module dkes_realspace
EOF
cat > dkes.f << "EOF"
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
EOF
cat > dkes_input_prepare.f << "EOF"
      subroutine dkes_input_prepare (arg, dkes_input_file)
c
c   This code prepares an input file for DKES based on data
c   in the boozmn.* file which is periodically written out
c   by the VMEC optimizer.  After compiling (see instructions
c   given below) and running a file called input.dkes is written.
c   DKES can then be run by typing "xdkes input.dkes".
c   This input file is based on a fixed set of parameters
c   specified below that is intended to provide
c   a rapidly evaluated (high collisionality) optimization target.
c   These can be modified to consider other parameter ranges and
c   a more complete mode spectrum as the need arises.  Parameters
c   which may be of interest to vary are:
c
c   The variables in the call to dkes_input_prepare are:
c
c     booz_file_name (=arg(1))
c             = name of boozmn file containing |B| data at various surfaces
c
c     nsurf (=arg(2))
c             = flux surface index (in boozer file) where
c                   DKES is to evaluate transport coefficients
c
c     cmul (arg(3))
c             = DKES collisionality parameter: nu/v  (in meter**-1)
c
c     efield (arg(4))
c             = DKES electric field parameter: E_s/v
c
c     lscreen (arg(5))
c             = logical variable which controls screen output
c            (= .true. screen output on, = .false. screen output off)
c
c
c
c     max_bmns = Number of Bmns which are retained for the
c                  DKES spectrum (this also influences the
c                  Fourier spectrum used for the distribution
c                  function in DKES)
c
c     legendre_modes = number of Legendre polynomials used in
c                  DKES to represent the pitch angle variation
c                  of the distribution function
c
c     coupling_order = Parameter for mode coupling order.
c                  This is the number of iterations used to
c                  generate the optimal m,n Fourier spectrum for
c                  the DKES distribution function.  At high
c                  collisionalities (cmul > 0.1), 2 iterations seems
c                  to be adequate (here adequate means that the upper
c                  and lower bounds coming out of DKES are close to
c                  to each other for the transport coefficient of
c                  interest).  As one goes to lower collisionalities,
c                  progressively more iterations must be used to get
c                  convergence of the upper/lower bounds.  This will
c                  imply increasingly longer run times for DKES.
c
c
c  Authors:  W.I. van Rij  --->  wrote original version of SPECFIL code
c            H. Maassberg  --->  Adapted SPECFIL to work within an interactive
c                                code to prepare input for DKES
c            D. Spong      --->  converted to f90, made to read new
c                                    form of Bmn file, different spectrum
c                                    generation rules tried, new output files,
c                                    interactive input turned off.  Adapted
c                                    for running with VMEC optimizer. Added
c                                    comments to the spectrum optimization.
c
c
      use kind_spec
      use date_and_computer
      use read_boozer_mod
      use read_wout_mod, rmnc_w=>rmnc, zmns_w=>zmns, lmns_w=>lmns,
     &    xm_w=>xm, xn_w=>xn, phip_w=>phip, mpol_w=>mpol,
     &    ntor_w=>ntor, nfp_w=>nfp, ns_w=>ns, mnmax_w=>mnmax
      use safe_open_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      character*(*) :: arg(5), dkes_input_file
c-----------------------------------------------
c   L o c a l   P a r a m e t e r s
c-----------------------------------------------
c     real(rprec), parameter :: radial_position = 0.5_dp
c     real(rprec), parameter :: cmul = 0.05_dp
c     real(rprec), parameter :: efield = 0._dp
      integer, parameter :: max_bmns = 20
      integer, parameter :: legendre_modes = 100
      integer, parameter :: coupling_order = 3
c
      integer, parameter :: mhi = 50
      integer, parameter :: ibbi = 1
      integer, parameter :: nhi = 50
      integer, parameter :: nfcd = 100
      integer, parameter :: mtd = 15
      integer, parameter :: nzd = 15
      integer, parameter :: mtd1 = mtd + 1
      integer, parameter :: nzd1 = 2*nzd + 1
      real(rprec), parameter :: aval1 = 1.e-4_dp
      real(rprec), parameter :: aval2 = 1.e-4_dp
      real(rprec), parameter :: one = 1, two = 2,
     1  half = 0.5_dp, zero = 0
      character*(50), parameter ::
     1   banner = 'THIS IS THE DKES input prep CODE Version 1.0'
c
c    aval1 = lower limit for the |B| Fourier coefficients
c    aval2 = lower limit for dominant zeta dependent spectrum
c         (note: either one or the other of these is active, but not
c          both.- DAS)
c    The parameters nhi and mhi are described in comments to follow.
c-----------------------------------------------
c   L o c a l   V a r i a b l e s
c-----------------------------------------------
      real(rprec), dimension(:), allocatable :: bmn_local,
     >   abs_blocal, blocal_ordered
      integer, dimension(:), allocatable :: m_ordered, n_ordered
      integer, dimension(:), allocatable :: nsurf
      integer :: iunit, js, iloc, n_norm, i, j, k, ierr, istat
      real(rprec) btheta, bzeta, check, iota_local, psip,
     >  chip, TWOPI
      real(rprec) phi_edge, cmul, efield, bmn_test
      integer, dimension(mhi,nhi,3) :: nplus, nmins
      integer, dimension(2*nhi,mhi) :: mns_out = 0
      integer, dimension(2*nhi) :: n_out = 0
      integer, dimension(nfcd) :: mtheta, ntheta, mzeta, nzeta, ivec
      integer, dimension(nhi) :: nvalsb
      integer, dimension(nfcd) :: intor, impol
      integer m, n, mp, nt, mp0, nt0, mpu, ntu, ntorb, mpolb,
     >  nszt, nsth, kk, jj, ii, ntor, mpol, ij, ik, im, in, iout,
     >  ier, nord, nmin, nmax, mmin, mmax, numargs
      real(rprec), dimension(nfcd) :: afc
      real(rprec), dimension(mtd1,nzd1) :: alf
      real(rprec) a, alfm, alfu, eps, au, aa, alfz
      real(rprec) hs, vnorm, vol_inner, vol_outer,
     >  rminor_i, rminor_o
      integer index_dat, index_end
      real(rprec) time_begin, time_end
      character*(120) :: booz_input_file, wout_input_file
      character*(10) :: date0, time0, zone0
      character*(40) :: dateloc
!     character*120 dkes_input_file
      integer :: imon, index_surf
      logical :: lscreen = .true.
c-----------------------------------------------
      TWOPI = 8*atan(one)

      call second0(time_begin)
c
      read(arg(2),'(i20)') js
      read(arg(3),'(f10.5)') cmul
      read(arg(4),'(f10.5)') efield
      if (arg(5)(1:1).eq.'f' .or. arg(5)(1:1).eq.'F') lscreen = .false.
c
      if (lscreen) then
         print 48
         call date_and_time(date0,time0,zone0)
         read (date0(5:6),'(i2)') imon
         write (dateloc,100) months(imon),date0(7:8),date0(1:4),
     1         time0(1:2),time0(3:4),time0(5:6)
         write (*,'(1x,a,/,1x,a)') banner, dateloc
         print *
         write(*,'(" js = ",i4,", ",i4," Bmns")') js,max_bmns
         write(*,'(" cmul = ",f10.6,", efield = ",f10.6)')cmul,efield
         write(*,'(" coupling order = ",i2,
     1    ", ",i4," Legendre modes")') coupling_order,legendre_modes
         print *
      end if
 100  format('DATE = ',a3,' ',a2,',',a4,' ',' TIME = ',2(a2,':'),a2)
 120  format(/,' TIME IN DKES input prep CODE:',1pe12.2,' SEC')
  48  format('====================================================')
c
      index_dat = index(arg(1),'.')
      index_end = len_trim(arg(1))
      booz_input_file  = arg(1)(index_dat+1:index_end)
      call read_boozer_file (booz_input_file, k)
      if (k .ne. 0) stop 'Error reading boozmn file in DKES_INPUT_PREP'
      allocate(bmn_local(mnboz_b), abs_blocal(mnboz_b),
     1    blocal_ordered(mnboz_b), stat=ierr)
      allocate(m_ordered(mnboz_b), n_ordered(mnboz_b), stat=ierr)
      allocate(nsurf(ns_b), stat=ierr)
c
c     Read data from the wout file and allocate storage:
c
      wout_input_file = 'wout.' // booz_input_file
      call readw_and_open(wout_input_file,k)
c
c     Find out which surfaces are not included in the boozmn file
c         (i.e., by checking if Bmn's are zero)
c
      index_surf = 0
      do j = 1, ns_b
         bmn_test = zero
         do i = 1,mnboz_b
            abs_blocal(i) = abs(bmn_b(i,j))
         end do
         bmn_test = sum(abs_blocal)
         if(bmn_test .gt. 1.e-8_dp) then
            index_surf = index_surf + 1
            nsurf(index_surf) = j
         endif
      end do


      do i = 1,mnboz_b
         bmn_local(i) = bmn_b(i,js)
         abs_blocal(i) = abs(bmn_local(i))
      end do
      bmn_test = sum(abs_blocal)
      if(bmn_test .lt. 1.e-8_dp .and. lscreen) then
         write(*,'(" This surface is not in the boozmn file")')
         write(*,'("  - need to pick a different surface")')
         write(*,'(" The following surfaces are available:")')
         write(*,'(6(3x,i2))') (nsurf(j), j=1,index_surf-1)
c
         call second0(time_end)
         if (lscreen) then
            write (*,120) time_end - time_begin
            write (*,48)
         end if
c
         stop 23
      endif
c
      iota_local = iota_b(js)
c      btheta = buco_b(js)
c      bzeta = bvco_b(js)
      btheta = half*(buco(js) + buco(js-1))
      bzeta = half*(bvco(js) + bvco(js-1))
      phi_edge = abs(phi_b(ns_b))
c
      hs = one/real((ns_w - 1), rprec)
      vnorm = TWOPI*TWOPI*hs
      vol_inner = sum(vp(2:js-1))
      vol_outer = vol_inner + vp(js)
      vol_inner = vol_inner*vnorm
      vol_outer = vol_outer*vnorm
      rminor_i = sqrt(two*vol_inner/(TWOPI*TWOPI*Rmajor))
      rminor_o = sqrt(two*vol_outer/(TWOPI*TWOPI*Rmajor))
      psip = -phip_w(js)*hs/(rminor_o - rminor_i)
      chip = psip*iota_local
c       write(*,'(" ns_w = ",i3)') ns_w

c
c     Sort the Bmn's (along with associated m's and n's in order
c     of increasing abs(Bmn):
c
      do i=1,mnboz_b
         check = maxval(abs_blocal)
c    Find location of this value of check in the abs_blocal array
         do j=1,mnboz_b
            if (abs(check - abs_blocal(j)) .lt. 1.e-6_dp*abs_blocal(j))
     1      iloc = j
         end do
         blocal_ordered(i) = bmn_local(iloc)
         m_ordered(i) = ixm_b(iloc)
         n_ordered(i) = ixn_b(iloc)
         abs_blocal(iloc) = zero
      end do
!
!     Find min/max m and n out of the first max_bmns of Bmn
!     modes. Compute ntorb and mpolb for DKES.
!
      nmin = 100
      mmin = 100
      nmax = -100
      mmax = -100
      do i=1,max_bmns
         n_norm = n_ordered(i)/nfp_b
         if(m_ordered(i) .lt. mmin) mmin = m_ordered(i)
         if(m_ordered(i) .gt. mmax) mmax = m_ordered(i)
         if(n_norm .lt. nmin) nmin = n_norm
         if(n_norm .gt. nmax) nmax = n_norm
      end do
      ntorb = nmax - nmin + 1
      mpolb = mmax - mmin + 1
!
!     Write out the initial part of the DKES input file.
!
      dkes_input_file = 'input_dkes.' // booz_input_file
      iunit = 15
      call safe_open(iunit, istat, dkes_input_file, 'replace', 
     1    'formatted')
      write (iunit,'(1x,"&datain")')
      write (iunit,'(1x,"nzperiod= ",i2,",")') nfp_b
      write (iunit,'(1x,"lalpha= ",i3,", nrun = 1,")') legendre_modes
      write (iunit,'(1x,"cmul = ",e12.4,",")') cmul
      write (iunit,'(1x,"efield = ",f7.4,",")') efield
      write (iunit,'(1x,"mpolb = ",i2,",",2x,"ntorb = ",
     1     i2,",",2x,"ibbi = 1,")') mpolb, ntorb
      write (iunit,'(1x,"chip = ",f7.4,",","  psip = ",f7.4,",")') 
     1     chip, psip
      write (iunit,'(1x,"btheta = ",f7.4,","," bzeta = ",f7.4,",")')
     1     btheta, bzeta
c
c     Write out only the top max_bmns largest Bmn's:
c
      do i=1,max_bmns
         n_norm = n_ordered(i)/nfp_b
         if(n_norm .ge. 0 .and. m_ordered(i) .le. 9) then
            write (iunit,'(" borbi(",i1,",",i1,")= ",1x,e11.5,",")')
     1      n_norm,m_ordered(i),blocal_ordered(i)
         else if (n_norm .lt. 0 .and. m_ordered(i) .le. 9) then
            write (iunit,'(" borbi(",i2,",",i1,")= ",1x,e11.5,",")')
     1      n_norm,m_ordered(i),blocal_ordered(i)
         else if(n_norm .ge. 0 .and. m_ordered(i) .gt. 9) then
            write (iunit,'(" borbi(",i1,",",i2,")= ",1x,e11.5,",")')
     1      n_norm,m_ordered(i),blocal_ordered(i)
         else if(n_norm .lt. 0 .and. m_ordered(i) .gt. 9) then
            write (iunit,'(" borbi(",i2,",",i2,")= ",1x,e11.5,",")')
     1      n_norm,m_ordered(i),blocal_ordered(i)
         end if
      end do
c
c     Translate m, n, and Bmn arrays into arrays used in the
c       GIDKES2 code:
c
      nord = coupling_order
      eps = max(zero,min(aval1,0.01_dp))
      alfm = zero
      do i=1,max_bmns
         mp = m_ordered(i)
         nt = n_ordered(i)/nfp_b
         a = blocal_ordered(i)
         impol(i) = mp
         intor(i) = nt
         afc(i) = a
         if (mp<=mtd .and. iabs(nt)<=nzd) alf(mp+1,nzd+1+nt) = a
         if (mp.eq.0 .and. nt.eq.0) then
            mp0 = mp
            nt0 = nt
         else if (abs(a) > alfm) then
            alfm = abs(a)
            alfu = a
            mpu = mp
            ntu = nt
         endif
      end do
      alfz = aval2
c     alfz = min(half*alfm,max(0.2_dp*alfm,10*eps))
c
c     The following loop assigns modes to S(zeta) (non-zero n modes through
c     the mzeta and nzeta arrays) and S(theta) (non-zero m modes through
c     the mtheta and ntheta arrays).  There is some flexibility as how this
c     is done through the alfz parameter defined above.  If alfz is
c     relatively large (i.e., if above alfz is used) not too many modes
c     make it into the mzeta and nzeta arrays, more go into the mtheta,
c     ntheta array.  If alfz is small (comment out the above alfz line)
c     then more modes go into the mzeta and nzeta arrays and not so many
c     into the mtheta and ntheta arrays.
c
c
c-----------------------------------------------------------------------
c
c     Description of the SPECFIL code:      (W.I. van Rij)
c
c     Purpose:
c     Estimate the Fourier spectrum for the DKES2 distribution function
c
c     Let S denote the fourier spectrum of 1/B [except for (0,0) term]
c     Let S(zeta) denote the dominant zeta-dependent spectrum of 1/B
c     Define S(theta) = S - S(zeta)
c
c     Let nszt (<=nhi) be the number of Fourier modes in S(zeta)
c     Let nsth (<=mhi) be the number of Fourier modes in S(theta)
c              (nhi, mhi are defined in parameter statement)
c
c     For i-th mode of S(zeta) load m value in mzeta(i) and n value in
c     nzeta(i) [i=1,2,....,nszt]
c     For i-th mode of S(theta) load m value in mtheta(i) and n value in
c     ntheta(i) [i=1,2,....,nsth]
c
c     n values are expressed in units of nzperiod
c
c     n values must satisfy   -nhi <= n <= nhi-1
c     m values must satisfy      0 <= m <= mhi-1
c     These bounds may be violated depending on the value of the
c     coupling order NORD (error message is printed on terminal).
c     In this case, you are prompted for a reduced value of NORD.
c
c-----------------------------------------------------------------------
c     subroutines required:
c     FILIT, SORTI
c-----------------------------------------------------------------------
c
      nszt = 0
      nsth = 0
      ntorb = 0
      do m = 0, mtd
         au = zero
         l12: do n = -nzd, nzd
            a = alf(m+1,nzd+1+n)
            aa = abs(a)
            if (aa > eps) then
               if (n.ne.0 .and. aa>alfz) then
                  nszt = nszt + 1
                  mzeta(nszt) = m
                  nzeta(nszt) = n
               else if (m.ne.0 .or. n.ne.0) then
                  nsth = nsth + 1
                  mtheta(nsth) = m
                  ntheta(nsth) = n
               endif
c
c       Keep track of the number of modes (mpolb, ntorb); also set
c       up the array of n's (nvalsb) as used in the old dkes input
c       file format (DAS).
c
               mpolb = m + 1
               if (ntorb .eq. 0) then
                  ntorb = 1
                  nvalsb(1) = n
               else
                  do i = 1, ntorb
                     if (nvalsb(i) .eq. n) cycle  l12
                  end do
                  ntorb = ntorb + 1
                  nvalsb(ntorb) = n
               endif
            endif
         end do l12
      end do
c
      call sorti (nvalsb, ntorb)
c
c     estimate the Fourier spectrum for the distribution function
c     -----   original SPECFIL version   -----
c
c     The nplus and nmins arrays are used to indicate where m,n
c     modes are used for the distribution function. nplus holds the
c     n .ge. 0 modes while nmins holds the n .lt. 0 modes. These
c     arrays contain the current 3 interates of S through the third
c     index.
c
      nord = max(2,nord)

      nplus(:mhi,:nhi,:) = 0
      nmins(:mhi,:nhi,:) = 0

      nplus(1,1,1) = 1
c
c     Generate S(0) = S(theta) + S(zeta)
c
      do i = 1, nsth
         m = mtheta(i)
         n = ntheta(i)
         call filit (m, n, 1, nplus, nmins, mhi, nhi, ier)
      end do
      if (ier .ne.0) stop

      do i = 1, nszt
         m = mzeta(i)
         n = nzeta(i)
         call filit (m, n, 1, nplus, nmins, mhi, nhi, ier)
      end do
      if (ier .ne.0) stop

      nplus(:mhi,:nhi,2) = nplus(:mhi,:nhi,1)
      nmins(:mhi,:nhi,2) = nmins(:mhi,:nhi,1)
c
c     Generate S(1) = S(0)*S(zeta)
c

      do i = 1, nszt
         do in = 1, nhi
            do im = 1, mhi
               if (nplus(im,in,1) .eq. 1) then
                  m = mzeta(i) + im - 1
                  n = nzeta(i) + in - 1
                  call filit (m, n, 2, nplus, nmins, mhi, nhi, ier)
                  if (ier .ne.0) stop
                  m = mzeta(i) - im + 1
                  n = nzeta(i) - in + 1
                  call filit (m, n, 2, nplus, nmins, mhi, nhi, ier)
                  if (ier .ne.0) stop
               endif
               if (nmins(im,in,1) .eq. 1) then
                  m = mzeta(i) + im - 1
                  n = nzeta(i) - in
                  call filit (m, n, 2, nplus, nmins, mhi, nhi, ier)
                  if (ier .ne.0) stop
                  m = mzeta(i) - im + 1
                  n = nzeta(i) + in
                  call filit (m, n, 2, nplus, nmins, mhi, nhi, ier)
                  if (ier .ne.0) stop
               endif
            end do
         end do
      end do

      kk = 2
      if (nord .ne.2) then
c
c
c     Generate S(n) = S(zeta)*S(n-1) + S(theta)*S(n-2):
c
c     Note: For the S(n) = S(zeta)*S(n-1) + S(theta)*S(n-2) recurrence
c      relation (this orders eps_tor << eps_hel) set ii = 0 and jj = 1 initially.
c      For the S(n) = S(zeta)*S(n-1) + S(theta)*S(n-1) recurrence
c      relation set ii = 1 and jj = 0 initially.
c
         ii = 1
         jj = 0
         do k = 3, nord
            ii = ii + 1
            if (ii .eq. 4) ii = 1
            jj = jj + 1
            if (jj .eq. 4) jj = 1
            kk = kk + 1
            if (kk .eq. 4) kk = 1
            nplus(:mhi,:nhi,kk) = nplus(:mhi,:nhi,jj)
            nmins(:mhi,:nhi,kk) = nmins(:mhi,:nhi,jj)
c
c        Generate the S(zeta)S(n-1) convolution:
c
            do i = 1, nszt
               do in = 1, nhi
                  do im = 1, mhi
                     if (nplus(im,in,jj) .eq. 1) then
                        m = mzeta(i) + im - 1
                        n = nzeta(i) + in - 1
                        call filit(m,n,kk,nplus,nmins,mhi,nhi,ier)
                        if (ier .ne.0) go to 72
                        m = mzeta(i) - im + 1
                        n = nzeta(i) - in + 1
                        call filit(m,n,kk,nplus,nmins,mhi,nhi,ier)
                        if (ier .ne.0) go to 72
                     endif
                     if (nmins(im,in,jj) .eq. 1) then
                        m = mzeta(i) + im - 1
                        n = nzeta(i) - in
                        call filit(m,n,kk,nplus,nmins,mhi,nhi,ier)
                        if (ier .ne.0) go to 72
                        m = mzeta(i) - im + 1
                        n = nzeta(i) + in
                        call filit(m,n,kk,nplus,nmins,mhi,nhi,ier)
                        if (ier .ne.0) go to 72
                     endif
                  end do
               end do
            end do
c
c        Generate the S(theta)S(n-2) convolution:
c

            do i = 1, nsth
               do in = 1, nhi
                  do im = 1, mhi
                     if (nplus(im,in,ii) .eq. 1) then
                        m = mtheta(i) + im - 1
                        n = ntheta(i) + in - 1
                        call filit(m,n,kk,nplus,nmins,mhi,nhi,ier)
                        if (ier .ne.0) go to 72
                        m = mtheta(i) - im + 1
                        n = ntheta(i) - in + 1
                        call filit(m,n,kk,nplus,nmins,mhi,nhi,ier)
                        if (ier .ne.0) go to 72
                     endif
                     if (nmins(im,in,ii) .eq. 1) then
                        m = mtheta(i) + im - 1
                        n = ntheta(i) + in
                        call filit(m,n,kk,nplus,nmins,mhi,nhi,ier)
                        if (ier .ne.0) go to 72
                        m = mtheta(i) - im + 1
                        n = ntheta(i) - in
                        call filit(m,n,kk,nplus,nmins,mhi,nhi,ier)
                        if (ier .ne.0) go to 72
                     endif
                  end do
               end do
            end do
         end do
      endif

      ntor = 0
      mpol = 0
      ii = 0
      do n = 1, nhi
         jj = 0
         do m = 1, mhi
            if (nplus(m,n,kk) .eq. 1) then
               jj = jj + 1
               ivec(jj) = m - 1
            endif
         end do
         if (jj .ne.0) then
            ii = ii + jj
            ntor = ntor + 1
            mpol = max(mpol,jj)
c         Capture the n's and m's into arrays n_out and mns_out
c           for later printout:
            n_out(n+nhi) = n-1
            do iout=1,jj
              mns_out(n+nhi,iout) = ivec(iout)
            end do
c
            ij = min(13,jj)
            if (ntor <= 9 .and. ntor >= -9) then
               write (iunit, 101) ntor, n - 1, jj, (ivec(i),i=1,ij)
            else
               write (iunit, 201) ntor, n - 1, jj, (ivec(i),i=1,ij)
            endif
            if (jj > 13) then
               ik = -4
   61          continue
               ik = ik + 18
               ij = min(jj,ik + 17)
               write (iunit, 202) (ivec(i),i=ik,ij)
               if (ij < jj) go to 61
            endif
c            endif
         endif
      end do

      do n = 1, nhi
         jj = 0
         do m = 1, mhi
            if (nmins(m,n,kk) .eq. 1) then
               jj = jj + 1
               ivec(jj) = m - 1
            endif
         end do
         if (jj .ne.0) then
            ii = ii + jj
            ntor = ntor + 1
            mpol = max(mpol,jj)
c         Capture the n's and m's into arrays n_out and mns_out
c           for later printout:
            n_out(nhi-n+1) = -n
            do iout=1,jj
              mns_out(nhi-n+1,iout) = ivec(iout)
            end do
c
            ij = min(13,jj)
            if (ntor <= 9 .and. ntor >= -9) then
               write (iunit, 101) ntor, (-n), jj, (ivec(i),i=1,ij)
            else
               write (iunit, 201) ntor, (-n), jj, (ivec(i),i=1,ij)
            endif
            if (jj > 13) then
               ik = -4
   64          continue
               ik = ik + 18
               ij = min(jj,ik + 17)
               write (iunit, 202) (ivec(i),i=ik,ij)
               if (ij < jj) go to 64
            endif
         endif
      end do
c
      write (iunit, 203) mpol, ntor
  101 format(' mmnn(1,',i1,')=',15(i3,','))
  201 format(' mmnn(1,',i2,')=',15(i3,','))
  202 format(1x,18(i3,','))
  203 format(' mpol=',i2,', ntor=',i3,' /')

      go to 73
  72  write(*,'("error in call to filit")')
  73  continue
      close(unit= iunit)
c
       call second0(time_end)
       if (lscreen) then
          write (*,120) time_end - time_begin
          write (*,48)
       end if
c

      deallocate(bmn_local, abs_blocal, blocal_ordered,
     >   stat=istat)
      deallocate(m_ordered, n_ordered, stat=istat)
      deallocate(nsurf, stat=istat)

      call read_boozer_deallocate
      call read_wout_deallocate

      end subroutine dkes_input_prepare


      subroutine filit(m, n, k, nplus, nmins, mhi, nhi, ier)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer m, n, k, mhi, nhi, ier
      integer, dimension(mhi,nhi,3) :: nplus, nmins
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: mp1, np1
C-----------------------------------------------

c-----------------------------------------------------------------------
c     auxilary routine for program GIDKES2
c     Author:     W.I. van Rij
c     modified:   HNM   STOP replaced by error parameter IER and RETURN
c-----------------------------------------------------------------------


      ier = 0
      if (m < 0) then
         m = -m
         n = -n
      endif

      mp1 = m + 1
      if (mp1 > mhi) then
         write (6, '(a,i3)') ' error in FILIT:   m+1 exceeds mhi = ',
     1      mhi
         ier = 1
         return
      endif

      if (n<0 .and. m.eq.0) n = -n

      if (n >= 0) then
         np1 = n + 1
         if (np1 > nhi) then
            write (6, '(a,i3)') ' error in FILIT:   n+1 exceeds nhi = '
     1         , nhi
            ier = 2
            return
         endif
         nplus(mp1,np1,k) = 1
         return

      else

         if ((-n) > nhi) then
            write (6, '(a,i3)') ' error in FILIT:   -n  exceeds nhi = '
     1         , nhi
            ier = 2
            return
         endif
         nmins(mp1,(-n),k) = 1
         return

      endif

      end subroutine filit


      SUBROUTINE SORTI(IX, NM)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER NM
      INTEGER, DIMENSION(*) :: IX
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ML, MG, IIL, IIG, ML1, MG1, NL, NG, N, II, IX0
C-----------------------------------------------

Ccom*    Sorting of an Integer Array (size)
C----------------------------------------------------------------------
C     Sorting of an integer array. On exit, the elements of IX increase
C     with index i.
C
C     IX      integer array
C     NM      no. of elements
C----------------------------------------------------------------------
Cend
C
      IF (NM <= 1) RETURN
      ML = 1
      MG = NM
      IF (NM .ne.2) THEN
C
    1    CONTINUE
         IIL = IX(ML)
         IIG = IIL
         ML1 = ML + 1
         MG1 = MG - 1
         NL = ML
         NG = ML
C
         DO N = ML1, MG
            II = IX(N)
            IF (II < IIL) THEN
               NL = N
               IIL = II
            ENDIF
            IF (II > IIG) THEN
               NG = N
               IIG = II
            ENDIF
         END DO
C
         IF (NL .ne.ML) THEN
            IX0 = IX(NL)
            IX(NL) = IX(ML)
            IX(ML) = IX0
         ENDIF
         IF (NG .ne.MG) THEN
            IF (NG .eq. ML) NG = NL
            IX0 = IX(NG)
            IX(NG) = IX(MG)
            IX(MG) = IX0
         ENDIF
         ML = ML1
         MG = MG1
         IF (MG - ML >= 2) GO TO 1
C
         IF (MG <= ML) RETURN
      ENDIF
      IF (IX(MG) >= IX(ML)) RETURN
      IX0 = IX(MG)
      IX(MG) = IX(ML)
      IX(ML) = IX0
      RETURN
      END SUBROUTINE SORTI
EOF
EOC
