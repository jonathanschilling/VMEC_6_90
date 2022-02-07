      subroutine autofill(X, DESX, NX)
      use kind_spec
      implicit none
!-----------------------------------------------
!   D u m m y   V a r i a b l e s
!-----------------------------------------------
      integer :: nx
      real(rprec) :: X(NX)
      character*(*) :: DESX     
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(rprec), parameter :: big_no = 1.e10_dp, em6 = 1.e-6_dp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(rprec) :: value
!-----------------------------------------------
      if (all(X(1:nx) .ge. big_no)) return
      if (any(X(1:nx) .ge. big_no)) then
         print *,' Auto-filling (in js) sigma_', trim(desx),' array!'
         value = maxval (x(1:nx), mask=x .lt. big_no)
         value = max(em6, value)
         where (x(1:nx) .ge. big_no) x = value
      endif

      end subroutine autofill


      subroutine safe_open_it(iunit, istat, filename, filestat, 
     1           fileform, record, access)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(inout) :: iunit
      integer, intent(out) :: istat
      integer, intent(in), optional :: record
      character*(*), intent(in) :: filename, filestat, fileform
      character*(*), intent(in), optional :: access
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      character*(*), parameter :: cdelim = 'apostrophe'
      character*(10) :: acc_type
      logical :: lopen, lexist, linvalid
C-----------------------------------------------
!
!     Make sure unit is not already opened
!
      linvalid = .true.
      do while (linvalid)
         inquire(iunit, exist=lexist, opened=lopen)
         linvalid = (.not.lexist) .or. lopen
         if (.not.linvalid) exit
         iunit = iunit+1
      end do

      lexist = (filestat(1:1).eq.'s') .or. (filestat(1:1).eq.'S')        !Scratch file

      if (present(access)) then
         acc_type = trim(access)
      else
         acc_type = 'SEQUENTIAL'
      end if
    
      if (fileform(1:1).eq.'u' .or. fileform(1:1).eq.'U') then
         if (present(record)) then
            if (lexist) then
            open(unit=iunit, form="unformatted", status="scratch", 
     1           recl=record, access=acc_type, iostat=istat)
            else
            open(unit=iunit, file=trim(filename), form="unformatted", 
     1           status=trim(filestat), recl=record,
     2           access=acc_type, iostat=istat)
            end if
         else
            if (lexist) then
            open(unit=iunit, form="unformatted", status="scratch", 
     1           access=acc_type, iostat=istat)
            else
            open(unit=iunit, file=trim(filename), form="unformatted", 
     1           status=trim(filestat), access=acc_type, iostat=istat)
            end if
         end if
      else
         if (present(record)) then
            if (lexist) then
            open(unit=iunit, form="formatted", status="scratch", 
     1           delim=trim(cdelim), recl=record, access=acc_type, 
     2           iostat=istat)
            else
            open(unit=iunit, file=trim(filename), form="formatted", 
     1           status=trim(filestat), delim=trim(cdelim), 
     2           recl=record, access=acc_type, iostat=istat)
            end if
         else
            if (lexist) then
            open(unit=iunit, form="formatted", status="scratch",
     1           delim=trim(cdelim), access=acc_type, iostat=istat)
            else
            open(unit=iunit, file=trim(filename), form="formatted", 
     1          status=trim(filestat), delim=trim(cdelim), 
     2          access=acc_type, iostat=istat)
            end if
         end if
      end if

      end subroutine safe_open_it


      subroutine second0(stime)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: count, count_rate
      real(rprec) :: stime
C-----------------------------------------------
      call system_clock(count, count_rate)
      if (count_rate .ne. 0) then
         stime = real(count, rprec)/count_rate
      else
         stime = 0
      end if
      
      end subroutine second0


      subroutine getcarg(index, arg, numargs)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: index, numargs
      character*(*), intent(out) :: arg
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: numchars, ier
C-----------------------------------------------
#if defined(WINNT)
      integer nargs
      numargs = nargs() - 1
      call getarg(index, arg, numchars)
#elif defined(LINUX)
      integer iargc
      numargs = iargc()
      call getarg(index, arg)
#elif defined(VMS)
      call lib$get_foreign(arg,,numchars)
      numargs = min(1,numchars)
#elif defined(LONESTAR) || defined(MCURIE)
      integer iargc
      numargs = iargc()
      call pxfgetarg(index, arg, numchars, ier)
#else
      integer iargc, getarg
      numargs = iargc()
      numchars = getarg(index, arg)
#endif

      end subroutine getcarg


      subroutine vmec_getenv(ename, evalue)
      implicit none
      character*(*) :: ename, evalue
#if defined(LONESTAR) || defined(MCURIE)
      INTEGER :: lenname=0, lenval, ierror
      call pxfgetenv(ename, lenname, evalue, lenval, ierror)
#else     
      call getenv(ename, evalue)
#endif      
      end subroutine vmec_getenv


      subroutine vmec_putenv(ename, evalue, ierror)
      implicit none
      integer :: ierror, len1
      character*(*), intent(in) :: ename, evalue
      character(len=200) :: temp
      integer, external :: putenv

      len1=len_trim(ename) + len_trim(evalue) + 2
      if (len1 .gt. 200) then
          print *,' error in vmec_putenv: ename+evalue too long'
          return
      end if
                 
      temp = trim(ename) // "=" // trim(evalue) // char(0)

      ierror = putenv(temp)

      end subroutine vmec_putenv


      subroutine vmec_system(cmd, ierror)
      integer, optional :: ierror
      integer :: ireturn
      character*(*), intent(in) :: cmd
      
#if defined(CRAY)
      integer, external :: ishell
      ireturn = ishell(cmd)
#elif defined(RISC)
      call system(cmd, ireturn)
#elif defined(LINUX) || defined(OSF1)
      integer, external :: system
      ireturn = system(trim(cmd))
#else
      integer, external :: system
      ireturn = system(trim(cmd) // char(0))
#endif
      if (present(ierror)) ierror = ireturn

      end subroutine vmec_system

      
      integer function vmec_chdir(new_path)
      implicit none
      character*(*), intent(in) :: new_path
      
#if defined(CRAY)
      integer :: ilen
      ilen = 0
      call pxfchdir(new_path, ilen, vmec_chdir)
#else
      integer, external :: chdir
      vmec_chdir = chdir(trim(new_path) // char(0))
#endif
      end function vmec_chdir


      subroutine vmec_flush(n)
      use kind_spec
      implicit none

C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: n


#ifdef RISC
      call flush_(n)
#else
      call flush(n)
#endif      
      end subroutine vmec_flush


      subroutine strip_comments(input_file)
      use safe_open_mod
      implicit none
!
!     strips comment lines (starting with '!') from input_file
!     renames clean file input_file // '.stripped'
!
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      character*(*) :: input_file
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, parameter :: unit_strip = 20
      integer :: istat, iustrip, iunew
      character*(200) :: line
      logical :: lex
C-----------------------------------------------
      inquire (file=input_file, exist = lex)
      if (.not. lex) then
         print *,trim(input_file),' does not exist!'
         stop
      end if
      iustrip = unit_strip
      call safe_open(iustrip, istat, input_file, 'old','formatted')
      if (istat .ne. 0) then       
        line = 'error opening ' // trim(input_file) // 
     1         ' in STRIP_COMMENTS'
        print *, line
        print *,'istat = ', istat
        stop
      end if  
      iunew = iustrip + 1
      call safe_open(iunew, istat, trim(input_file) // '.stripped',
     1     'replace', 'formatted')
      if (istat .ne. 0) then       
        line = 'error opening ' // trim(input_file) // 
     1      '.stripped  in STRIP_COMMENTS'
        print *, line
        print *,'istat = ', istat
        stop
      end if  
      do 
         read(iustrip, '(a)', end=100) line
         line = adjustl(line)
         if (line(1:1) == '!') cycle
         write(iunew, '(a)') trim(line)
      end do

 100  continue
       
      close(iustrip)
      close(iunew)
      
      end subroutine strip_comments


      subroutine multxy(x,y,n)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer n
      real(rprec), dimension(n) :: x, y
c-----------------------------------------------

      x = x * y

      end subroutine multxy


      subroutine solver(amat, b, m)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: m
      real(rprec), intent(inout) :: amat(m,*), b(m,*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, parameter :: nrhs = 1
      integer :: info 
      integer, allocatable :: ipiv(:)
C-----------------------------------------------
      EXTERNAL GESV_G
C-----------------------------------------------
      info = 0
      allocate (ipiv(m))

!     Compute the solution to a real system of linear equations
!       A * X = B,
!     where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
!
!     FACTOR AMATRIX INTO LU FORM
!     AND SOLVE BY GAUSSIAN ELIMINATION
!
      CALL GESV_G (m, nrhs, amat, m, ipiv, b, m, info)
      if (info .ne. 0) print *, ' Condition No. = 0 in Solver'
      
      deallocate (ipiv)
 
      end subroutine solver


      function fmin (ax, bx, f, tol)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec) ax, bx, f, tol
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec),parameter :: zero = 0,one = 1, two = 2,
     1   three = 3, five = 5, p5 = 0.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec):: a,b,c,d,e,eps,xm,p,q,r,tol1,tol2,
     1   u,v,w,fu,fv,fw,fx,x,fmin
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      external f
C-----------------------------------------------
c
c  c is the squared inverse of the golden ratio
c
      c = p5*(three - sqrt(five))
c
c  eps is approximately the square root of the relative machine
c  precision.
c
      eps = one
      eps = eps/two
      tol1 = one + eps
 1003 continue
      if (tol1 .le. one) go to 1002
      eps = eps/two
      tol1 = one + eps
      if (tol1 .le. one) go to 1002
      eps = eps/two
      tol1 = one + eps
      if (tol1 .le. one) go to 1002
      eps = eps/two
      tol1 = one + eps
      if (tol1 .le. one) go to 1002
      eps = eps/two
      tol1 = one + eps
      go to 1003
 1002 continue
      eps = sqrt(eps)
c
c  initialization
c
      a = ax
      b = bx
      v = a + c*(b - a)
      w = v
      x = v
      e = zero
      fx = f(x)
      fv = fx
      fw = fx
c
c  main loop starts here
c
   20 continue
      xm = p5*(a + b)
      tol1 = eps*abs(x) + tol/three
      tol2 = two*tol1
c
c  check stopping criterion
c
      if (abs(x - xm) .le. tol2 - p5*(b - a)) go to 90
c
c is golden-section necessary
c
      if (abs(e) .le. tol1) go to 40
c
c  fit parabola
c
      r = (x - w)*(fx - fv)
      q = (x - v)*(fx - fw)
      p = (x - v)*q - (x - w)*r
      q = two*(q - r)
      if (q .gt. zero) p = -p
      q = abs(q)
      r = e
      e = d
c
c  is parabola acceptable
c
      if (abs(p) .ge. abs(p5*q*r)) go to 40
      if (p .le. q*(a - x)) go to 40
      if (p .ge. q*(b - x)) go to 40
c
c  a parabolic interpolation step
c
      if( q.ne.zero )d = p/q
      u = x + d
c
c  f must not be evaluated too close to ax or bx
c
      if (u - a < tol2) d = sign(tol1,xm - x)
      if (b - u < tol2) d = sign(tol1,xm - x)
      go to 50
c
c  a golden-section step
c
   40 continue
      if (x .ge. xm) then
         e = a - x
      else
         e = b - x
      endif
      d = c*e
c
c  f must not be evaluated too close to x
c
   50 continue
      if (abs(d) .ge. tol1) then
         u = x + d
      else
         u = x + sign(tol1,d)
      endif
      fu = f(u)
c
c  update  a, b, v, w, and x
c
      if (fu .gt. fx) go to 60
      if (u .ge. x) then
         a = x
      else
         b = x
      endif
      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu
      go to 20
   60 continue
      if (u .lt. x) then
         a = u
      else
         b = u
      endif
      if (fu .le. fw) go to 70
      if (w .eq. x) go to 70
      if (fu .le. fv) go to 80
      if (v .eq. x) go to 80
      if (v .eq. w) go to 80
      go to 20
   70 continue
      v = w
      fv = fw
      w = u
      fw = fu
      go to 20
   80 continue
      v = u
      fv = fu
      go to 20
c
c  end of main loop
c
   90 continue
      fmin = x

      end function fmin

       
      function fmax (ax, bx, f, tol)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec) ax, bx, f, tol
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero=0, one=1, two=2,
     1   three = 3, five = 5, p5 = 0.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec):: a,b,c,d,e,eps,xm,p,q,r,tol1,tol2,
     1   u,v,w,fu,fv,fw,fx,x, fmax
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      external f
C-----------------------------------------------
c
c  c is the squared inverse of the golden ratio
c
      c = p5*(three - sqrt(five))
c
c  eps is approximately the square root of the relative machine
c  precision.
c
      eps = one
      eps = eps/two
      tol1 = one + eps
 1003 continue
      if (tol1 .le. one) go to 1002
      eps = eps/two
      tol1 = one + eps
      if (tol1 .le. one) go to 1002
      eps = eps/two
      tol1 = one + eps
      if (tol1 .le. one) go to 1002
      eps = eps/two
      tol1 = one + eps
      if (tol1 .le. one) go to 1002
      eps = eps/two
      tol1 = one + eps
      go to 1003
 1002 continue
      eps = sqrt(eps)
c
c  initialization
c
      a = ax
      b = bx
      v = a + c*(b - a)
      w = v
      x = v
      e = zero
      fx = f(x)
      fv = fx
      fw = fx
c
c  main loop starts here
c
   20 continue
      xm = p5*(a + b)
      tol1 = eps*abs(x) + tol/three
      tol2 = two*tol1
c
c  check stopping criterion
c
      if (abs(x - xm) .le. tol2 - p5*(b - a)) go to 90
c
c is golden-section necessary
c
      if (abs(e) .le. tol1) go to 40
c
c  fit parabola
c
      r = (x - w)*(fx - fv)
      q = (x - v)*(fx - fw)
      p = (x - v)*q - (x - w)*r
      q = two*(q - r)
cwie      if (q.gt.0.0) p = -p
      if (q < zero) p = -p
      q = abs(q)
      r = e
      e = d
c
c  is parabola acceptable
c
cwie   30 if (abs(p).ge.abs(0.5*q*r)) goto 40
cwie      if (p.le.q*(a - x)) goto 40
cwie      if (p.ge.q*(b - x)) goto 40
      if (abs(p) .le. abs(p5*q*r)) go to 40
      if (p .ge. q*(a - x)) go to 40
      if (p .le. q*(b - x)) go to 40
c
c  a parabolic interpolation step
c
      d = p/q
      u = x + d
c
c  f must not be evaluated too close to ax or bx
c
      if (u - a < tol2) d = sign(tol1,xm - x)
      if (b - u < tol2) d = sign(tol1,xm - x)
      go to 50
c
c  a golden-section step
c
   40 continue
      if (x .ge. xm) then
         e = a - x
      else
         e = b - x
      endif
      d = c*e
c
c  f must not be evaluated too close to x
c
   50 continue
      if (abs(d) .ge. tol1) then
         u = x + d
      else
         u = x + sign(tol1,d)
      endif
      fu = f(u)
c
c  update  a, b, v, w, and x
c
cwie      if (fu.gt.fx) goto 60
      if (fu < fx) go to 60
      if (u .ge. x) then
         a = x
      else
         b = x
      endif
      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu
      go to 20
   60 continue
      if (u < x) then
         a = u
      else
         b = u
      endif
cwie      if (fu.le.fw) goto 70
      if (fu .ge. fw) go to 70
      if (w .eq. x) go to 70
cwie      if (fu.le.fv) goto 80
      if (fu .ge. fv) go to 80
      if (v .eq. x) go to 80
      if (v .eq. w) go to 80
      go to 20
   70 continue
      v = w
      fv = fw
      w = u
      fw = fu
      go to 20
   80 continue
      v = u
      fv = fu
      go to 20
c
c  end of main loop
c
   90 continue
      fmax = x

      end function fmax


      subroutine tridslv(a, d, b, c, jmin, jmax, mnd1, ns, nrhs)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: jmax, mnd1, ns, nrhs
      integer, dimension(0:mnd1), intent(in) :: jmin
      real(rprec), dimension(ns,0:mnd1) :: a, d, b
      real(rprec), dimension(ns,0:mnd1, nrhs) :: c
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: mn, in, i0, in1, jrhs
      real(rprec), allocatable, dimension(:,:) :: alf
      real(rprec), dimension(0:mnd1) :: psi0
C-----------------------------------------------
!
!     SOLVES B(I)*X(I-1)+D(I)*X(I)+A(I)*X(I+1)=C(I), I=IN,JMAX
!     AND RETURNS ANSWER IN C(I)
!     ADDED VECTORIZATION ON FOURIER MODE ARGUMENT (1/2000)
!     AND NEW ARGUMENT (NRHS) TO DO MULTIPLE RIGHT SIDES SIMULTANEOUSLY
!
      if (jmax .gt. ns) stop 'jmax>ns in tridslv'

      allocate (alf(ns,0:mnd1), stat = in)
      if (in .ne. 0) stop 'Allocation error in tridslv'

      in = minval(jmin)
!
!      FILL IN MN BELOW MAX(JMIN) WITH DUMMY VALUES
!      TO ALLOW VECTORIZATION ON MN INDEX
!
      do mn = 0, mnd1
         in1 = jmin(mn)-1
         if (in1 .ge. in) then
            d(in:in1, mn) = 1
            c(in:in1, mn, 1:nrhs) = 0
            b(in:in1, mn) = 0
            a(in:in1, mn) = 0
         end if
      end do

      in1 = in + 1

      psi0(:)= d(in,:)
!     if (any(psi0 .eq. zero)) stop 'error in tridslv'
      psi0 = one/psi0
      do jrhs = 1, nrhs
         c(in,:,jrhs) = c(in,:,jrhs)*psi0(:)
      end do

      do i0 = in1,jmax
         alf(i0-1,:) = a(i0-1,:)*psi0(:)
         psi0(:)  = one/(d(i0,:) - b(i0,:)*alf(i0-1,:))
         do jrhs = 1, nrhs
            c(i0,:,jrhs) = (c(i0,:,jrhs) - b(i0,:)*c(i0-1,:,jrhs))
     1                   * psi0(:)
         end do
      end do
        
      do i0 = jmax - 1, in, -1
         do jrhs = 1,nrhs
            c(i0,:,jrhs) = c(i0,:,jrhs) - alf(i0,:)*c(i0+1,:,jrhs)
         end do
      end do

      deallocate (alf)

      end subroutine tridslv


      subroutine read_boozer_file(extension, ierr, iopen)
      use read_boozer_mod, except_this => read_boozer_file    !!Avoid redundant interface
      use safe_open_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: ierr
      integer, optional :: iopen
      character*(*) :: extension
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: unit_booz = 14
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nsval, mn, num, iunit
      character*38 :: version
C-----------------------------------------------
!
!     Read in bmn''s from boozmn.extension file
!
      iunit = unit_booz
      
      call safe_open(iunit, ierr, 'boozmn.'//extension, 'old',
     1     'unformatted')
      if (present(iopen)) iopen = ierr
       
      if (ierr .ne. 0) return

      read(iunit, iostat=ierr, err=100) nfp_b, ns_b, aspect_b, 
     1   rmax_b, rmin_b, betaxis_b

      if (allocated(iota_b)) call read_boozer_deallocate
      allocate (iota_b(ns_b), pres_b(ns_b), beta_b(ns_b), phip_b(ns_b),
     1  phi_b(ns_b), bvco_b(ns_b), buco_b(ns_b), idx_b(ns_b), stat=ierr)
      if (ierr .ne. 0) then 
        print *,' Allocation error in read_boozer_file' 
        return
      end if     
      iota_b(1) = 0; pres_b(1) = 0; beta_b(1) = 0
      phip_b(1) = 0; phi_b(1) = 0; bvco_b(1) = 0
      buco_b(1) = 0
      
      do nsval = 2, ns_b
         read(iunit, iostat=ierr, err=100) iota_b(nsval), 
     1   pres_b(nsval), beta_b(nsval), phip_b(nsval), phi_b(nsval),
     2   bvco_b(nsval), buco_b(nsval)
      end do

      read(iunit, iostat=ierr, err=100) mboz_b, nboz_b, mnboz_b
      read(iunit, iostat=ierr, err=100) version
              
      allocate (bmn_b(mnboz_b,ns_b), rmnc_b(mnboz_b,ns_b),
     1  zmns_b(mnboz_b,ns_b), pmns_b(mnboz_b,ns_b), gmn_b(mnboz_b,ns_b),
     2  ixm_b(mnboz_b), ixn_b(mnboz_b), stat = ierr)   
      if (ierr .ne. 0) then
         print *,' Allocation error in read_boozer_file'
         return
      end if

      idx_b = 0
      ixm_b = 0
      rmnc_b = 0; zmns_b = 0; pmns_b = 0; bmn_b = 0; gmn_b = 0
      num = 1
      
      do while (ierr .eq. 0)
        read(iunit, iostat=ierr, end=200, err=100) nsval
        if (nsval .gt. ns_b) then
           print *,' nsval > ns_b in read_boozer'
           cycle
        end if   
        idx_b(nsval) = 1
        
        if (num .eq. 1) then
           do mn = 1, mnboz_b
              read(iunit, iostat=ierr, err=100) ixn_b(mn), ixm_b(mn)
           end do   
           num = num + 1
        end if

        do mn = 1, mnboz_b
           read(iunit, iostat=ierr, err=100, end=200) bmn_b(mn,nsval), 
     1       rmnc_b(mn,nsval), zmns_b(mn,nsval), pmns_b(mn,nsval),
     2       gmn_b(mn,nsval)
        end do
      end do

 100  continue
      if (ierr .gt. 0)
     1    print *,' Error reading in subroutine read_boozer_file:',
     2            ' ierr = ', ierr
 200  continue      
      if (ierr .lt. 0) ierr = 0       !End-of-file, ok
      close(iunit)
      
      end subroutine read_boozer_file
 

      subroutine readw_and_open(filename, ierr, iopen)
      use safe_open_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(out) :: ierr
      integer, optional :: iopen
      character*(*), intent(in) :: filename
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, parameter :: iunit_init = 10
      integer :: iunit
C-----------------------------------------------
*
*     THIS SUBROUTINE READS THE TEXT FILE WOUT CREATED BY THE VMEC CODE
*     AND STORES THE INFORMATION IN THE READ_WOUT MODULE
*
      iunit = iunit_init
      call safe_open (iunit, ierr, filename, 'old', 'formatted')
      if (present(iopen)) iopen = ierr
      if (ierr .ne. 0) return
      
      call readw_only_priv(iunit, ierr)
      
      close(unit=iunit)

      end subroutine readw_and_open
 

      subroutine readw_only(iunit, ierr, iopen)
      use read_wout_mod, except_this => readw_only       !!Avoid redundant interface
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: iunit, ierr
      integer, optional :: iopen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      character*256 :: vmec_version
      logical exfile
C-----------------------------------------------
*
      ierr = 0

      inquire(unit=iunit, exist=exfile, name=vmec_version)
      if (.not.exfile) then
        print *,' In READ_WOUT_FILE, Unit = ',iunit,
     1          ' File = ',trim(vmec_version),' DOES NOT EXIST'
        if (present(iopen)) iopen = -1
        ierr = -1
        return
      else
        if (present(iopen)) iopen = 0
      end if
      
      call readw_only_priv(iunit, ierr)

      end subroutine readw_only

      
      subroutine readw_only_priv(iunit, ierr)
      use read_wout_mod, except_this => readw_only       !!Avoid redundant interface
      use vparams, ONLY: dmu0
      use vsvd0, only: nigroup
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: iunit, ierr
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: eps_w = 1.e-4_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat(15), i, j, k, js, m, n, n1, mn
      character*(256) :: vmec_version
C-----------------------------------------------
*
*     THIS SUBROUTINE READS THE TEXT FILE WOUT CREATED BY THE VMEC CODE
*     AND STORES THE INFORMATION IN THE READ_WOUT MODULE
*   
*     CALL READ_WOUT_FILE - GENERIC INTERFACE - CAN BE CALLED WITH EITHER UNIT NO. OR FILENAME
*
*     RMNC, ZMNS: FULL-GRID
*     LMNS      : HALF-GRID
*
      istat = 0
      ierr = 0
      nextcur = 0

      read (iunit, '(a)', iostat=istat(2), err=1000) vmec_version

      i = index(vmec_version,'=')
      if (i .ge. 0) then
         read(vmec_version(i+1:len_trim(vmec_version)),*) version_
      else 
         version_ = -1.0
      end if      

      ierr_vmec = 0    
 
      if (version_ .le. (5.10 + eps_w)) then
         read (iunit, *, iostat=istat(2), err=1000) wb, wp, gamma, 
     1      pfac, nfp, ns, 
     1      mpol, ntor, mnmax, itfsq, niter, iasym, ireconstruct
      else 
         if (version_ .lt. 6.54) then
            read (iunit, *, iostat=istat(2), err=1000) wb, wp, gamma, 
     1         pfac, rmax_surf, rmin_surf
         else
            read (iunit, *, iostat=istat(2), err=1000) wb, wp, gamma, 
     1         pfac, rmax_surf, rmin_surf, zmax_surf
         end if
         read (iunit, *, iostat=istat(2), err=1000) nfp, ns, mpol, ntor,
     1     mnmax, itfsq, niter, iasym, ireconstruct, ierr_vmec
      end if  
     

      if (version_ .gt. (6.20+eps_w)) then      
         read (iunit, *, iostat=istat(1), err=1000)imse, itse, nbsets, 
     1      nobd, nextcur, nstore_seq
      else
         read (iunit, *, iostat=istat(1), err=1000)imse, itse, nbsets, 
     1        nobd, nextcur
         nstore_seq = 100
      end if   

      if (ierr_vmec.ne.0 .and. ierr_vmec.ne.4) go to 1000

      if (nextcur .gt. nigroup) istat(15) = -1
      

      if (allocated(nbfld) .or. allocated(xm)) call read_wout_deallocate

      allocate (nbfld(nbsets), stat = istat(3))

      allocate (xm(mnmax), xn(mnmax), rmnc(mnmax,ns), zmns(mnmax,ns),
     1  lmns(mnmax,ns), rmns(mnmax,ns), zmnc(mnmax,ns), lmnc(mnmax,ns),
     2  bmn(mnmax,ns), gmn(mnmax,ns), bsubumn(mnmax,ns),
     2  bsubvmn(mnmax,ns), bsubsmn(mnmax,ns), bsupumn(mnmax,ns),
     3  bsupvmn(mnmax,ns), currvmn(mnmax,ns), iotas(ns), mass(ns), 
     4  pres(ns), beta_vol(ns), phip(ns), buco(ns), bvco(ns), phi(ns),
     5  vp(ns), overr(ns), jcuru(ns), jcurv(ns), specw(ns), Dmerc(ns),
     6  Dshear(ns), Dwell(ns), Dcurr(ns), Dgeod(ns), equif(ns),
     7  extcur(nextcur), curlabel(nextcur), raxis(0:ntor,2), 
     8  zaxis(0:ntor,2), jdotb(ns), bdotgradv(ns), 
     9  fsqt(nstore_seq), wdot(nstore_seq), stat = istat(6))
     

      rmns = 0; zmnc = 0; lmnc = 0; fsqt = 0; wdot = 0
      
      if (nbsets .gt. 0) read (iunit, *, iostat=istat(4), err=1000) 
     1   (nbfld(i),i=1,nbsets)
      read (iunit, '(a)', iostat=istat(5), err=1000) mgrid_file

      do js = 1, ns
         do mn = 1, mnmax
            if(js .eq. 1) then
               read (iunit, *, iostat=istat(7), err=1000) m, n
               xm(mn) = real(m,rprec)
               xn(mn) = real(n,rprec)
            end if   
            if (version_ .le. (6.20+eps_w)) then
              read (iunit, 730, iostat=istat(8), err=1000) 
     1        rmnc(mn,js), zmns(mn,js), lmns(mn,js),
     2        bmn(mn,js), gmn(mn,js), bsubumn(mn,js), bsubvmn(mn,js),
     3        bsubsmn(mn,js), bsupumn(mn,js), bsupvmn(mn,js), 
     4        currvmn(mn,js)
            else
              read (iunit, *, iostat=istat(8), err=1000) 
     1        rmnc(mn,js), zmns(mn,js), lmns(mn,js),
     2        bmn(mn,js), gmn(mn,js), bsubumn(mn,js), bsubvmn(mn,js),
     3        bsubsmn(mn,js), bsupumn(mn,js), bsupvmn(mn,js), 
     4        currvmn(mn,js)
            end if
            if (iasym .gt. 0) then
               read (iunit, *, iostat=istat(8), err=1000) 
     1         rmns(mn,js), zmnc(mn,js), lmnc(mn,js)
            end if
            if (js.eq.1 .and. m.eq.0) then
               n1 = abs(n/nfp)
               if (n1 .le. ntor) then
                  raxis(n1,1) = rmnc(mn,1)
                  zaxis(n1,1) = zmns(mn,1)
               end if   
            end if
         end do
      end do

!
!     Read HALF-MESH QUANTITIES (except jcuru, jcurv which are FULL MESH)
!
!     NOTE: In version_ <= 6.00, mass, press were written out in INTERNAL (VMEC) units
!     and are therefore multiplied here by 1/mu0 to transform to pascals. Same is true
!     for all the currents (jcuru, jcurv, jdotb). Also, in version_ = 6.10 and
!     above, PHI is the true (physical) toroidal flux (has the sign of jacobian correctly
!     built into it)
!
      iotas(1) = 0; mass(1) = 0; pres(1) = 0; phip(1) = 0; 
      buco(1) = 0; bvco(1) = 0; phi(1) = 0; vp(1) = 0; overr(1) = 0
      jcuru(1) = 0; jcurv(1) = 0; specw(1) = 0

      if (version_ .le. (6.05+eps_w)) then
         read (iunit, 730, iostat=istat(9), err=1000) 
     1     (iotas(js), mass(js), pres(js), 
     2      phip(js), buco(js), bvco(js), phi(js), vp(js), overr(js),
     3      jcuru(js), jcurv(js), specw(js),js=2,ns)
         read (iunit, 730, iostat=istat(10), err=1000) 
     1      aspect, betatot, betapol, betator, betaxis, b0
      else if (version_ .le. (6.20+eps_w)) then
         read (iunit, 730, iostat=istat(9), err=1000) 
     1     (iotas(js), mass(js), pres(js), beta_vol(js),
     2      phip(js), buco(js), bvco(js), phi(js), vp(js), overr(js),
     3      jcuru(js), jcurv(js), specw(js),js=2,ns)
         read (iunit, 730, iostat=istat(10), err=1000) 
     1      aspect, betatot, betapol, betator, betaxis, b0
      else
         read (iunit, *, iostat=istat(9), err=1000) 
     1     (iotas(js), mass(js), pres(js), beta_vol(js),
     2      phip(js), buco(js), bvco(js), phi(js), vp(js), overr(js),
     3      jcuru(js), jcurv(js), specw(js),js=2,ns)
         read (iunit, *, iostat=istat(10), err=1000) 
     1      aspect, betatot, betapol, betator, betaxis, b0
      end if


      if (version_ .gt. (6.10+eps_w)) then
         read (iunit, *, iostat=istat(10), err=1000) isigng
         read (iunit, *, iostat=istat(10), err=1000) input_extension
         read (iunit, *, iostat=istat(10), err=1000) IonLarmor, 
     1     VolAvgB, RBtor0, RBtor, Itor, Aminor, Rmajor, Volume
      end if


!-----------------------------------------------
!     MERCIER CRITERION
!-----------------------------------------------
      if (version_.gt.(5.10+eps_w) .and. version_.lt.(6.20-eps_w)) then
         read (iunit, 730, iostat=istat(11), err=1000) 
     1      (Dmerc(js), Dshear(js), Dwell(js), Dcurr(js),
     2       Dgeod(js), equif(js), js=2,ns-1)
      else if (version_ .ge. (6.20-eps_w)) then        
         read (iunit, *, iostat=istat(11), err=1000) 
     1      (Dmerc(js), Dshear(js), Dwell(js), Dcurr(js),
     2       Dgeod(js), equif(js), js=2,ns-1)
      end if     

      if (nextcur .gt. 0) then
         if (version_ .le. (6.20+eps_w)) then
            read (iunit, 730, iostat=istat(12), err=1000)
     1      (extcur(i),i=1,nextcur)
         else
            read (iunit, *, iostat=istat(12), err=1000)
     1      (extcur(i),i=1,nextcur)
         end if
         read (iunit, *, iostat=istat(13), err=1000)
     1      (curlabel(i),i=1,nextcur)
      endif

      if (version_ .le. (6.20+eps_w)) then
         read (iunit, 730, iostat=istat(14)) 
     1     (fsqt(i), wdot(i), i=1,nstore_seq)
      else
         read (iunit, *, iostat=istat(14)) 
     1     (fsqt(i), wdot(i), i=1,nstore_seq)
      end if

      if ((version_.ge.6.20-eps_w) .and. (version_ .lt. (6.50-eps_w))
     1   .and. (istat(14).eq.0)) then
         read (iunit, 730, iostat=istat(14), err=1000) 
     1     (jdotb(js), bdotgradv(js),js=1,ns)
      else if (version_ .ge. (6.50-eps_w)) then
         read (iunit, *, iostat=istat(14), err=1000) 
     1     (jdotb(js), bdotgradv(js),js=1,ns)
      else
        istat(14) = 0
      endif  
!
!     CONVERT FROM INTERNAL UNITS TO PHYSICAL UNITS IF NEEDED
!
      if (version_ .le. (6.05+eps_w)) then
         mass = mass/dmu0
         pres = pres/dmu0
         jcuru = jcuru/dmu0
         jcurv = jcurv/dmu0
         jdotb = jdotb/dmu0
         phi   = -phi
      end if          

!-----------------------------------------------
!     DATA AND MSE FITS
!-----------------------------------------------
      if (ireconstruct .gt. 0) then

        n1 = maxval(nbfld(:nbsets))
        allocate (sknots(isnodes), ystark(isnodes), y2stark(isnodes), 
     1     pknots(ipnodes), ythom(ipnodes), y2thom(ipnodes), 
     2     anglemse(2*ns), rmid(2*ns), qmid(2*ns), shear(2*ns), 
     3     presmid(2*ns), alfa(2*ns), curmid(2*ns), rstark(imse),
     4     datastark(imse), rthom(itse), datathom(itse),
     5     dsiext(nobd), plflux(nobd), dsiobt(nobd), bcoil(n1,nbsets), 
     6     plbfld(n1,nbsets), bbc(n1,nbsets))
         if (imse.ge.2 .or. itse.gt.0) then
            read (iunit, *) tswgt, msewgt
            read (iunit, *) isnodes, (sknots(i),ystark(i),y2stark(i),
     1         i=1,isnodes)
            read (iunit, *) ipnodes, (pknots(i), ythom(i), 
     1         y2thom(i),i=1,ipnodes)
            read(iunit, *)(anglemse(i),rmid(i),qmid(i),shear(i),
     1      presmid(i),alfa(i),curmid(i),i=1,2*ns-1)
            read(iunit, *)(rstark(i),datastark(i),qmeas(i),i=1,imse)
            read(iunit, *)(rthom(i),datathom(i),i=1,itse)
         endif
  
         if (nobd .gt. 0) then
            read (iunit, *) (dsiext(i),plflux(i),dsiobt(i),i=1,nobd)
            read (iunit, *) flmwgt
         endif

         nbfldn = sum(nbfld(:nbsets))
         if (nbfldn .gt. 0) then
            do n = 1, nbsets
               read (iunit, *) (bcoil(i,n),plbfld(i,n),bbc(i,n),
     1            i=1,nbfld(n))
            end do
            read (iunit, *) bcwgt
         endif

         read (iunit, *) phidiam, delphid
!
!     read Limiter & Prout plotting specs
!
         read (iunit, *) nsets, nparts, nlim

         allocate (nsetsn(nsets))
         read (iunit, *) (nsetsn(i),i=1,nsets)

         n1 = maxval(nsetsn(:nsets))
         allocate (pfcspec(nparts,n1,nsets), limitr(nlim))

         read (iunit, *) (((pfcspec(i,j,k),i=1,nparts),
     1      j=1,nsetsn(k)),k=1,nsets)

         read (iunit, *) (limitr(i), i=1,nlim)

         m  = maxval(limitr(:nlim))
         allocate (rlim(m,nlim), zlim(m,nlim)) 

         read (iunit, *) ((rlim(i,j),zlim(i,j),i=1,limitr(j)),
     1      j=1,nlim)
         read (iunit, *) nrgrid, nzgrid
         read (iunit, *) tokid
         read (iunit, *) rx1, rx2, zy1, zy2, condif
         read (iunit, *) imatch_phiedge

      end if
      
 1000 continue

      if (istat(2) .ne. 0) ierr_vmec = 1

      do m = 1,15
        if (istat(m) .ne. 0) then
           print *,' Error No. ',m,' in READ_WOUT, iostat = ',istat(m)
           ierr = m
           exit
        end if  
      end do  


  720 format(8i10)
  730 format(5e20.13)
  740 format(a)
  790 format(i5,/,(1p3e12.4))

      end subroutine readw_only_priv


      subroutine write_array(iunit, name, array, n)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: iunit, n
      character*(*), intent(in) :: name
      real(rprec), dimension(n), intent(in) :: array
C-----------------------------------------------
C   L o c a l  V a r i a b l e s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0
      integer :: start_count, end_count, length, line
      real(rprec) :: temp
      character*(*), parameter :: adyes = "yes", adno = "no"
      character*(3) :: adv
C-----------------------------------------------
!
!      This subroutine writes repeated namelist entries in N*Value format to save space
!      in namelist file connected to "iunit"
!
  980 format(2x,a,' = ')
  982 format(2x,i4,'*',1pe20.14)
  983 format(2x,i4,'*',1pe21.14)
  984 format(2x,1pe20.14)
  985 format(2x,1pe21.14)
  990 format(2x,a,' = ',3(1pe22.14))
  992 format(2x,a,' = ',i4,'*',1pe20.14)
  993 format(2x,a,' = ',i4,'*',1pe21.14)

      if( n == 1) then
         write(iunit,990) name, array(1)

      else if( all(array(1)==array(2:n))) then
        if( array(1) >= zero ) then
          write(iunit,992) name, n, array(1)
        else
          write(iunit,993) name, n, array(1)
        endif

      else
!       Look for repeats within the array
         start_count = 1
         line = 1
         write(iunit, 980, advance="no") name
         do while (start_count .le. n)
            temp = array(start_count)
            end_count = start_count
            do while (end_count .lt. n)
               if (array(end_count+1) .ne. temp) exit
               end_count = end_count+1
            end do
!              Limit no. records/line AND start new record for each array
!              IF maximum packing desired, eliminate the end_count==n in following test            
            if (line==2 .or. end_count.eq.n) then
               adv = "yes"
               line = 0
            else
               adv = "no"
               line = line+1
            end if
            length = end_count - start_count + 1
            if (length > 1) then
               if (temp >= zero) then
                  write (iunit, 982, advance=trim(adv)) length, temp
               else
                  write (iunit, 983, advance=trim(adv)) length, temp
               end if
            else
               if (temp >= zero) then
                  write (iunit, 984, advance=trim(adv)) temp
               else
                  write (iunit, 985, advance=trim(adv)) temp
               end if
            end if
            start_count = end_count+1
         end do
      endif

      end subroutine write_array
      

      subroutine write_larray(iunit, name, array, n)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: iunit, n
      character*(*), intent(in) :: name
      logical, dimension(n), intent(in) :: array
C-----------------------------------------------
C   L o c a l  V a r i a b l e s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0
      integer :: start_count, end_count, length, line
      logical :: temp
      character*(*), parameter :: adyes = "yes", adno = "no"
      character*(3) :: adv
C-----------------------------------------------
!
!      This subroutine writes repeated namelist entries in N*Value format to save space
!      in namelist file connected to "iunit"
!
  980 format(2x,a,' = ')
  982 format(2x,i4,'*',l1)
  984 format(2x,l1)
  990 format(2x,a,' =',35(1x,l1))
  992 format(2x,a,' = ',i4,'*',l1)

      if( n == 1) then
         write(iunit,990) name, array(1)

      else if( all(array(1) .eqv. array(2:n))) then
         write(iunit,992) name, n, array(1)

      else
!       Look for repeats within the array
         start_count = 1
         line = 1
         write(iunit, 980, advance="no") trim(name)
         line = 5 + len_trim(name)
         do while (start_count .le. n)
            temp = array(start_count)
            end_count = start_count

            do while (end_count .lt. n)
               if (array(end_count+1) .neqv. temp) exit
               end_count = end_count+1
            end do
            length = end_count - start_count + 1

            if( length > 1) then
               line = line + 8
            else
               line = line + 3
            endif

!              Limit no. records/line AND start new record for each array
!              IF maximum packing desired, eliminate the end_count==n in following test            
            if (line>=70 .or. end_count.eq.n) then
               adv = "yes"
               line = 0
            else
               adv = "no"
            end if

            if (length > 1) then
               write (iunit, 982, advance=trim(adv)) length, temp
            else
               write (iunit, 984, advance=trim(adv)) temp
            end if
            start_count = end_count+1
         end do
      endif

      end subroutine write_larray


      subroutine read_namelist(iunit, io_stat, lc_name)
      use vmec_input, only: read_indata_namelist, 
     1   read_mse_namelist, ns_array
      use vmec_seq, only: vseq
      use bootsj_input, only: read_boot_namelist
      use optim_params, only: read_optimum_namelist, sigma_bmax, 
     1    sigma_bmin, sigma_bmn, sigma_ripple, sigma_jstar, lprof_opt,
     2    lcurprof_opt, numjstar, nsd    
      use coilsnamin, only: read_coils_namelist
      use gade_mod, only: read_gade_namelist
      use read_namelist_mod, only: autofill
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer iunit, io_stat
      character*(*) :: lc_name
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ifind, ieof, ilc, iuc, nsmax
      character*(132) :: line, nameuc, namelc
!-----------------------------------------------
 
      io_stat = -1
      rewind (iunit)
      namelc = '&'//adjustl(lc_name)
#if defined(CRAY) || defined(SUNOS)
!
!     work-around Cray, Sun bug: cannot read multiple namelists in same file!
!
      nameuc = namelc
      iuc = ichar('A')
      ilc = ichar('a')
      do ifind = 1,len_trim(namelc)
        ieof = ichar(namelc(ifind:ifind))
        if (ieof.ge.ilc .and. ieof.le.ichar('z')) then
           nameuc(ifind:ifind) = char(ieof+iuc-ilc)
        else
           nameuc(ifind:ifind) = char(ieof)
        end if      
      end do
      
      ifind = 0
      ieof = 0
 
      do while(ifind.eq.0 .and. ieof.eq.0)
         read (iunit, '(a)', iostat=ieof) line
         ifind = index(line,trim(nameuc)) + index(line,trim(namelc))
      end do
 
      if (ifind .ne. 0) then
         backspace (iunit)
      else
         return 
      endif
#endif
      ifind = min(len_trim(namelc), 132)
      if (namelc(1:ifind) .eq. '&indata') then
         call read_indata_namelist (iunit, io_stat)

      else if (namelc(1:ifind) .eq. '&optimum') then
         call read_optimum_namelist (iunit, io_stat)
         if (io_stat .gt. 0) return
!
!        obsolete assignments
!
         if (lcurprof_opt) lprof_opt = .true.
         nsmax = min(nsd, maxval(ns_array))
         do iuc = 1,NumJstar
            call autofill(sigma_Jstar(1:nsmax,iuc), 'Jstar', nsmax)
         end do   
         call autofill(sigma_bmin, 'bmin', nsmax)
         call autofill(sigma_bmax, 'bmax', nsmax)
         call autofill(sigma_ripple, 'ripple', nsmax)
         call autofill(sigma_bmn, 'bmn', nsmax)

      else if (namelc(1:ifind) .eq. '&bootin') then
         call read_boot_namelist (iunit, io_stat)
      else if (namelc(1:ifind) .eq. '&mseprofile') then
         call read_mse_namelist (iunit, io_stat)
      else if (namelc(1:ifind) .eq. '&vseq') then
         read (iunit, nml=vseq, iostat=io_stat)
      else if (namelc(1:ifind) .eq. '&coilsin') then
         call read_coils_namelist (iunit, io_stat) 
      else if (namelc(1:ifind) .eq. '&ga_de') then
         call read_gade_namelist (iunit, io_stat)
      end if
      
#if defined(IRIX64)      
      if (io_stat .eq. -1) io_stat = 0                  !Catches EOF on ORIG2000 Machine
#endif      
      end subroutine read_namelist

#if !defined(CRAY) && !defined(IRIX64)
      subroutine pxffork(ipid, ierror)
      implicit none
      integer :: ipid, ierror
      integer, external :: fork

      ierror = 0

      ipid = fork()

      if (ipid < 0) ierror = -ipid

      end subroutine pxffork
      

      subroutine pxfgetpid(ipid, ierror)
      implicit none
      integer :: ipid, ierror
      integer, external :: getpid

      ierror = 0

      ipid = getpid()

      if (ipid < 0) ierror = -ipid

      end subroutine pxfgetpid
      

      subroutine pxfwait(istat, iretpid, ierror)
      implicit none
      integer :: istat, iretpid, ierror
      integer, external :: wait

      iretpid = 0
      istat = 0

      ierror = wait(0)

      end subroutine pxfwait
#endif
      subroutine convert_boundary(rbc, zbs, rhobc, mpol, ntor)
      use kind_spec
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: mpol, ntor
      real(rprec), dimension(-ntor:ntor,0:mpol), intent(in) ::
     1   rbc, zbs
      real(rprec), dimension(-ntor:ntor,0:mpol), intent(out) ::
     1   rhobc
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: pexp = 4
      integer :: mcount, ncount, m_bdy, n_bdy
      real(rprec), parameter :: p25 = 0.25_dp, p50 = 0.50_dp,
     1  zero = 0, one = 1
      real(rprec), dimension(0:mpol) :: t1, t2
!-----------------------------------------------
!
!     GIVEN A BOUNDARY REPRESENTATION OF THE FORM
!     
!     R = RBC(n,m)*COS(mu-nv)
!     Z = ZBS(n,m)*SIN(mu-nv)
!
!     CONVERTS (APPROXIMATELY) TO FIND POLAR RADIUS ARRAY RHOBC
!     REPRESENTATION (FOR M>0 MODES) USING HIRSHMAN/BRESLAU
!     PRESCRIPTION WITH EXPONENT = PEXP
!
!     THIS WOULD MAKE A GOOD INITIAL GUESS FOR DESCUR CODE
!     CHECKED THAT IF THE BOUNDARY IS IN THE DESIRED FORM, IT
!     WILL NOT BE CHANGED BY A CALL TO THIS CODE!
!

!     First determine maximum m-number in boundary representation
!     DO NOT exceed this mmax-1 in rhobc
!
      m_bdy = 0
      n_bdy = 0
      do mcount = 1, mpol
         do ncount = -ntor, ntor
            if (rbc(ncount,mcount).ne.zero .or.
     1          zbs(ncount,mcount).ne.zero) then
                m_bdy = max(m_bdy,mcount)
                n_bdy = max(n_bdy,abs(ncount))
            end if
         end do
      end do                          

      rhobc = zero
      
      do mcount = 1, mpol
        t1(mcount) = ( real(mcount-1,rprec)
     1        /real(mcount,rprec) )**pexp
        t2(mcount) = ( real(mcount+1,rprec)
     1        /real(mcount,rprec) )**pexp
      end do
      t1(1) = one
      
!
!     NOTE: Rhobc(n,m=0) is different for n>0 and n<0, since
!           it is a linear combination of k|n| +- rho(0,|n|)
!

      do ncount = -n_bdy, n_bdy
        do mcount = 0,1
          rhobc(ncount,mcount) = p50*(rbc(ncount,mcount+1) +
     1                               zbs(ncount,mcount+1))/t1(mcount+1)
        end do   

        do mcount = 2,m_bdy-1
          rhobc(ncount,mcount) = p25*(
     1    (rbc(ncount,mcount+1) + zbs(ncount,mcount+1))/t1(mcount+1)
     2  + (rbc(ncount,mcount-1) - zbs(ncount,mcount-1))/t2(mcount-1))
        end do     
      end do  
      
      end subroutine convert_boundary
 

      subroutine convert_boundary_PG(rbc, zbs, rhobc, mpol, ntor)
      use kind_spec
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: mpol, ntor
      real(rprec), dimension(-ntor:ntor,0:mpol) ::
     1   rbc, zbs
      real(rprec), dimension(-ntor:ntor,-mpol:mpol), intent(out) ::
     1   rhobc
      real(rprec) :: rnorm, r00
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: pexp = 4
      integer :: mcount, ncount, m_bdy, n_bdy
      real(rprec), parameter :: p25 = 0.25_dp, p50 = 0.50_dp
     1  , zero = 0.0_dp, one = 1.0_dp
      real(rprec), dimension(0:mpol) :: t1, t2
!-----------------------------------------------
!
!     GIVEN A BOUNDARY REPRESENTATION OF THE FORM
!     
!     R = RBC(n,m)*COS(mu-nv)
!     Z = ZBS(n,m)*SIN(mu-nv)
!
!
!     CONVERT TO GARABEDIAN DELTA REPRESENTATION
!
      m_bdy = 0
      n_bdy = 0
      do mcount = 1, mpol
         do ncount = -ntor, ntor
            if (rbc(ncount,mcount).ne.zero .or.
     1          zbs(ncount,mcount).ne.zero) then
                m_bdy = max(m_bdy,mcount)
                n_bdy = max(n_bdy,abs(ncount))
            end if
         end do
      end do

      if(m_bdy+1 .gt. mpol) then
      write(6,*) "In Conversion to Delta-mn, mpol too small"
      stop
      endif

      rhobc = zero
      r00 = rbc(0,0)
      rnorm = rbc(0,1)+zbs(0,1)
      rnorm = 2/rnorm

      do mcount = 0, mpol
         do ncount = -ntor, ntor
            rbc(ncount,mcount)=rbc(ncount,mcount)*rnorm
            zbs(ncount,mcount)=zbs(ncount,mcount)*rnorm
         enddo
      enddo

       do mcount = 0, m_bdy
       do ncount = -n_bdy, n_bdy
         rhobc(ncount,mcount+1) =
     >                  0.5*(rbc(ncount,mcount)-zbs(ncount,mcount))
     >                + rhobc(ncount,mcount+1)
         rhobc(-ncount,-mcount+1) =
     >                  0.5*(rbc(ncount,mcount)+zbs(ncount,mcount))
     >                + rhobc(-ncount,-mcount+1)
       enddo
       enddo
!
!      restore r(0,0) and keep fixed
!
       rhobc(0,0) = r00

!      do mcount = -mmax, mmax
!      do ncount = -nmax, nmax
!      if( abs(rhobc(ncount,mcount)) .gt. 1.0e-10)
!    > write(*,*) mcount, ncount, rhobc(ncount,mcount)
!      enddo
!      enddo
      
      end subroutine convert_boundary_PG


      subroutine unique_boundary(rbc, zbs, rhobc, nmax, mmax, mpol,ntor)
      use kind_spec
      implicit none
!-----------------------------------------------
C   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: nmax, mmax, mpol, ntor
      real(rprec), dimension(-nmax:nmax,0:mmax), intent(inout) ::
     1   rhobc
      real(rprec), dimension(-nmax:nmax,0:mmax), intent(out) ::
     1   rbc, zbs
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: pexp = 4
      integer :: mcount, ncount, nrho
      real(rprec), parameter :: zero = 0, one = 1
      real(rprec), dimension(0:mpol) :: t1, t2
!-----------------------------------------------
!
!     GIVEN A RADIUS ARRAY RHOBC, COMPUTES A UNIQUE BOUNDARY
!     REPRESENTATION (FOR M>0 MODES) USING HIRSHMAN/BRESLAU
!     PRESCRIPTION
!
      do mcount = 1, mpol
        t1(mcount) = ( real(mcount-1,rprec)/
     1                 real(mcount,rprec) )**pexp
        t2(mcount) = ( real(mcount+1,rprec)/
     1                 real(mcount,rprec) )**pexp
      end do
      t1(1) = one
      
      rbc(:,1:mmax) = zero
      zbs(:,1:mmax) = zero
      
!
!     NOTE: RHOBC(n,0) includes both signs of n, since
!     it is a linear combination of k(n) and rho(|n|,0)
!     We need to impose the constraint rbc(n,1) - rbc(-n,1) = -(zbs(n,1) - zbs(-n,1))
!     corresponding to rhobs(n,0) = 0. In terms of rhobc, this becomes
!     rhobc(0,n) = rhobc(0,-n)
!
      do ncount = -ntor,-1
          rhobc(ncount,0) = rhobc(-ncount,0)      !m=1 rbc,zbs constraint
      end do

      do ncount = -ntor, ntor
        do mcount = 1,mpol-2
          rbc(ncount,mcount) = t1(mcount)*rhobc(ncount,mcount-1)
     1                       + t2(mcount)*rhobc(ncount,mcount+1)           
          zbs(ncount,mcount) = t1(mcount)*rhobc(ncount,mcount-1)
     1                       - t2(mcount)*rhobc(ncount,mcount+1)           
        end do
        do mcount = max(2,mpol-1),mpol
          rbc(ncount,mcount) = t1(mcount)*rhobc(ncount,mcount-1)
          zbs(ncount,mcount) = rbc(ncount,mcount)
        end do
      end do  

      end subroutine unique_boundary


      subroutine unique_boundary_PG
     >                          (rbc, zbs, rhobc, nmax, mmax, mpol,ntor)
      use kind_spec
      implicit none
!-----------------------------------------------
C   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: nmax, mmax, mpol, ntor
      real(rprec), dimension(-nmax:nmax,-mmax:mmax), 
     >                  intent(inout) ::  rhobc
      real(rprec), dimension(-nmax:nmax,0:mmax), intent(out) ::
     1   rbc, zbs
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer :: mcount, ncount, nrho
      real(rprec), parameter :: zero = 0, one = 1
      real(rprec), dimension(0:mpol) :: t1, t2
      integer :: mb, nb, m_bdyn, m_bdyp, n_bdy
      real(rprec) :: rnorm, r00
!-----------------------------------------------

      rbc(:,0:mpol) = zero
      zbs(:,0:mpol) = zero

      m_bdyn = 0
      m_bdyp = 0
      n_bdy = 0
      do mcount = -mmax, mmax
         do ncount = -nmax, nmax
            if (rhobc(ncount,mcount) .ne. zero ) then
                m_bdyn = min(m_bdyn,mcount)
                m_bdyp = max(m_bdyp,mcount)
                n_bdy = max(n_bdy,abs(ncount))
            end if
         end do
      end do
!
       rnorm = rhobc(0,0)/rhobc(0,1)
!
!      note: rhobc(0,0)=1.0 in delta representation
!            we use it to temporarily store the
!            actual major radius
!
       r00 = rhobc(0,0)
       rhobc(0,0) = one

       do mcount = m_bdyn, 0
       mb = iabs(mcount) + 1
       do ncount = -n_bdy, n_bdy
       nb = -ncount
       rbc( nb, mb) = rbc( nb, mb) + rhobc( ncount, mcount)
       zbs( nb, mb) = zbs( nb, mb) + rhobc( ncount, mcount)
       enddo
       enddo
       do mcount = 1, m_bdyp
       mb = mcount - 1
       do ncount = -n_bdy, n_bdy
       nb =  ncount
       rbc( nb, mb) = rbc( nb, mb) + rhobc( ncount, mcount)
       zbs( nb, mb) = zbs( nb, mb) - rhobc( ncount, mcount)
       enddo
       enddo
       do ncount = 1, ntor
       rbc( ncount, 0 ) = rbc( ncount, 0 ) + rbc( -ncount, 0 )
       rbc( -ncount, 0 ) = zero
       zbs( ncount, 0 ) = zbs( ncount, 0 ) - zbs( -ncount, 0 )
       zbs( -ncount, 0 ) = zero
       enddo
       zbs( 0, 0 ) = zero
 
       do mcount = 0, mpol
          do ncount = -ntor, ntor
             rbc(ncount,mcount)=rnorm*rbc(ncount,mcount)
             zbs(ncount,mcount)=rnorm*zbs(ncount,mcount)
          enddo
       enddo
       rbc(0,0) = r00
       rhobc(0,0) = r00


      end subroutine unique_boundary_PG
