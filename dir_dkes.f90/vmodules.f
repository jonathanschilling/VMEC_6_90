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
