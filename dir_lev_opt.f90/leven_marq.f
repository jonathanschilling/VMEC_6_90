!-----------------------------------------------
!
!     THE FOLLOWING ROUTINES COMPRISE THE LEVENBERG-MARQUARDT OPTIMIZATION SUITE
!
!-----------------------------------------------
#ifdef MPI_OPT
      subroutine lmdif1(fcn, m, n, x, fvec, tol, epsfcn,
     1   nfev_end, diag, mode, info, lwa)
!
!!    THIS IS THE MULTIPLE PROCESSOR VERSION OF LEV_OPT USING MPI CALLS
!!          (D. A. Spong 10/27/00)
      
!!!   epsfcn was added as a dummy argument, eliminating the circular
!!!   module reference. 
!!    diag and mode were added as arguments, so that different choices
!!    for these can be tried in the future (as suggested by L. A. Berry)
!!    Currently mode is set = 1 in the calling routine (stellopt) so that
!!    diag is internally generated.  See definitions of mode and diag below.
!!

#else
      subroutine lmdif1(fcn, m, n, x, fvec, tol, epsfcn, nfev_end,
     1           diag, mode, info, lwa, max_processors, num_lm_params)
!!    ADDED EPSFCN TO ARG LIST: BY SPH (2/97)
!!    ADDED NFEV_END == MAXFEV TO ARG LIST (6/31/99)
!!    ADDED MAX_PROCESSORS, NUM_LM_PARAMS TO ARG LIST (11/23/99)
!!    (NEED MAX_PROCESSORS, NUM_LM_PARAMS FOR MULTI-PROCESSOR APPLICATIONS)

      use fdjac_mod, ONLY: maxj_processors=>max_processors,
     1     numj_lm_params=>num_lm_params      
#endif
      use lmpar_mod, diag_mod=>diag, x_mod=>x, fvec_mod=>fvec, 
     1               m_mod=>m, n_mod=>n

      use kind_spec
      implicit none


C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: m, n, lwa, nfev_end, mode
      integer, intent(out) :: info
      real(rprec), intent(in) :: tol, epsfcn
      real(rprec), dimension(n), intent(inout) :: x, diag
      real(rprec), dimension(m), intent(out) :: fvec
#ifndef MPI_OPT
      integer, intent(in) :: max_processors, num_lm_params
#endif
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, factor = 100
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: maxfev, mp5n, nfev, nprint
      real(rprec) :: ftol, gtol, xtol
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      external fcn
C-----------------------------------------------
c
c     subroutine lmdif1
c
c     the purpose of lmdif1 is to minimize the sum of the squares of
c     m nonlinear functions in n variables by a modification of the
c     levenberg-marquardt algorithm. this is done by using the more
c     general least-squares solver lmdif. the user must provide a
c     subroutine which calculates the functions. the jacobian is
c     then calculated by a forward-difference approximation.
c
c     the subroutine statement is
c
c       subroutine lmdif1(fcn,m,n,x,fvec,tol,info,lwa)
c
c     where
c
c       fcn is the name of the user-supplied subroutine which
c         calculates the functions. fcn must be declared
c         in an external statement in the user calling
c         program, and should be written as follows.
c
c         subroutine fcn(m, n, x, fvec, iflag, ncnt)
c         integer m,n,iflag
c         real(rprec) x(n),fvec(m)
c         ----------
c         calculate the functions at x and
c         return this vector in fvec.
c         ----------
c         return
c         end
c
c         the value of iflag should not be changed by fcn unless
c         the user wants to terminate execution of lmdif1.
c         in this case set iflag to a negative integer. On a multi-processor
c         machine, iflag will be initialized to the particular processor id.
c
c
c       m is a positive integer input variable set to the number
c         of functions.
c
c       n is a positive integer input variable set to the number
c         of variables. n must not exceed m.
c
c       x is an array of length n. on input x must contain
c         an initial estimate of the solution vector. on output x
c         contains the final estimate of the solution vector.
c
c       fvec is an output array of length m which contains
c         the functions evaluated at the output x.
c
c       ncnt is a positive integer input variable set to the current
c         iteration count (added by SPH - 7/99)
c
c       tol is a nonnegative input variable. termination occurs
c         when the algorithm estimates either that the relative
c         error in the sum of squares is at most tol or that
c         the relative error between x and the solution is at
c         most tol.
c
c       info is an integer output variable. if the user has
c         terminated execution, info is set to the (negative)
c         value of iflag. see description of fcn.
c
c       lwa is a positive integer input variable not less than
c         m*n+5*n+m.
c
c     subprograms called
c
c       user-supplied ...... fcn
c
c       minpack-supplied ... lmdif
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     modified to accept improvements from jacobian calc. MZarnstorff Oct 2001
c     modified to flip sign of jacobian offsets for more efficient search
c     and start with an exponential levenberg spread to settle on scale
c             M. Zarnstorff                          Jan 2002
c
c
c     modified to simplify storage and arguments, to promote clarity, debugging
c     and reduce issues with multiple paths and names for the same variables.
c     Also, MPI version improved to have similar convergence as fork version,
c     and to not need pure-supervisor process..
c             M.Zarnstorff                           April 2006
c     **********

      info = 0
c
c     check the input parameters for errors.
c
      if (n.le.0 .or. m<n .or. tol<zero .or. lwa<m*n+5*n+m) then
         print *,' Input parameter error in LMDIF1 (levenburg/marq)'
         return
      end if
c
c     call lmdif.
c
      allocate (x_mod(n), wa1(n), wa2(n), wa3(n), wa4(m), ipvt(n),  
     1          diag_mod(n), qtf(n), fvec_mod(m), fjac(m,n), 
     2          stat = info)      
      if (info .ne. 0) stop 'Allocation error in lmdif1!' 

#ifndef MPI_OPT
!
!     Load fdjac module values
!
      maxj_processors = max(max_processors,1)
      numj_lm_params  = max(num_lm_params,1)
#endif      

      maxfev = 200*(n + 1)
      maxfev = min (maxfev, nfev_end)             !!SPH-Added 7/99
      ftol = tol
      xtol = tol
      gtol = zero

      m_mod = m
      n_mod = n
      x_mod = x
      diag_mod = diag

!!    ADDED BY SPH -- PASSED IN ARG LIST(2/97)
!!    epsfcn = zero
!!    mode = 1       (DAS, passed through arg list 9/13/00)

      nprint = 0
      mp5n = m + 5*n
      call lmdif (fcn, ftol, xtol, gtol, maxfev, 
     1     epsfcn, mode, factor, nprint, info, nfev)

#ifndef MPI_OPT
      if (info .eq. 8) info = 4
#endif

      x = x_mod
      diag = diag_mod
      fvec = fvec_mod

      deallocate (x_mod, wa1, wa2, wa3, wa4, ipvt,  
     1          diag_mod, qtf, fvec_mod, fjac )      

      end subroutine lmdif1

       
      subroutine lmdif(fcn, ftol, xtol, gtol, maxfev,
     1                 epsfcn, mode, factor, nprint, info, nfev)
      use kind_spec
      use lmpar_mod
#ifdef MPI_OPT
      use fdjac_mod, ONLY: flip, ix_min, jac_order, jac_count
#else
      use fdjac_mod, ONLY: max_processors, flip, ix_min, jac_order,
     1                     jac_count, num_lm_params
#endif
      implicit none
#ifdef MPI_OPT
      include 'mpif.h'                                       !mpi stuff
#endif
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: maxfev, mode, nprint, info, nfev
      real(rprec), intent(in) ::  ftol, xtol, gtol, epsfcn, factor
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: FLAG_SINGLETASK = -1, FLAG_CLEANUP = -100
      integer, parameter :: master = 0                       !only mpi uses this...
      real(rprec), parameter :: zero = 0, one = 1,
     1   p1=0.1_dp, p5=0.5_dp, p25=0.25_dp, p75=0.75_dp, p0001=1.e-4_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: iflag, iter, j, l, istat, lev_step_range
      real(rprec) :: actred, dirder, epsmch, fnorm,
     1   gnorm, prered, ratio, sum0, temp,
     2   temp1, temp2, xnorm, delta_old, actred_lev, par_old
      real(rprec), dimension(:,:), allocatable :: fjac_save
#ifdef MPI_OPT
      integer :: ierr
#else
      real(rprec) :: wall_time, wall_time_lev
#endif
      real(rprec) :: fnorm_min
      real(rprec), allocatable, dimension(:) :: x_min, fvec_min
      integer :: cycle_count, subcycle
      character*130, dimension(0:9) :: info_array
      logical :: step_improved, first_jacobian
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      external fcn
      real(rprec), external :: dpmpar, enorm
C-----------------------------------------------
      data info_array/
     1     'improper input parameters (number constraints MUST be greate
     1r than number variables (>0).',

     1     'algorithm estimates that the relative error in the sum of sq
     1uares is at most ftol.',

     2     'algorithm estimates that the relative error between x and th
     2e solution is at most xtol.', 

     3     'algorithm estimates that the relative error in the sum of sq
     3uares and between x and the solution is at most ftol and xtol.',
     
     4     'the cosine of the angle between fvec and any column of the j
     4acobian is at most gtol in absolute value.',

     5     'number of calls to fcn has reached or exceeded maxfev.',

     6     'ftol is too small. no further reduction in the sum of square
     6s is possible.',

     7     'xtol is too small. no further improvement in the approximate 
     7 solution x is possible.',
     
     8     'gtol is too small. fvec is orthogonal to the columns of the 
     8jacobian to machine precision.',
     
     9     'levenberg-marquardt optimizer terminated properly.' /

c
c     subroutine lmdif
c
c     the purpose of lmdif is to minimize the sum of the squares of
c     m nonlinear functions in n variables by a modification of
c     the levenberg-marquardt algorithm. the user must provide a
c     subroutine which calculates the functions. the jacobian is
c     then calculated by a forward-difference approximation.
c
c     the subroutine statement is
c
c       subroutine lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,
c                        diag,mode,factor,nprint,info,nfev,fjac,
c                        ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
c
c     where
c
c       fcn is the name of the user-supplied subroutine which
c         calculates the functions. fcn must be declared
c         in an external statement in the user calling. See LMDIF1 for
c         documentation and should be written as follows.
c
c         subroutine fcn(m, n, x, fvec, iflag, ncnt)
c         integer m,n,iflag
c         real(rprec) x(n),fvec(m)
c         ----------
c         calculate the functions at x and
c         return this vector in fvec.
c         ----------
c         return
c         end
c
c
c       ftol is a nonnegative input variable. termination
c         occurs when both the actual and predicted relative
c         reductions in the sum of squares are at most ftol.
c         therefore, ftol measures the relative error desired
c         in the sum of squares.
c
c       xtol is a nonnegative input variable. termination
c         occurs when the relative error between two consecutive
c         iterates is at most xtol. therefore, xtol measures the
c         relative error desired in the approximate solution.
c
c       gtol is a nonnegative input variable. termination
c         occurs when the cosine of the angle between fvec and
c         any column of the jacobian is at most gtol in absolute
c         value. therefore, gtol measures the orthogonality
c         desired between the function vector and the columns
c         of the jacobian.
c
c       maxfev is a positive integer input variable. termination
c         occurs when the number of calls to fcn is at least
c         maxfev by the end of an iteration.
c
c       epsfcn is an input variable used in determining a suitable
c         step length for the forward-difference approximation. this
c         approximation assumes that the relative errors in the
c         functions are of the order of epsfcn. if epsfcn is less
c         than the machine precision, it is assumed that the relative
c         errors in the functions are of the order of the machine
c         precision.
c
c       diag is an array of length n. if mode = 1 (see
c         below), diag is internally set. if mode = 2, diag
c         must contain positive entries that serve as
c         multiplicative scale factors for the variables.
c
c       mode is an integer input variable. if mode = 1, the
c         variables will be scaled internally. if mode = 2,
c         the scaling is specified by the input diag. other
c         values of mode are equivalent to mode = 1.
c
c       factor is a positive input variable used in determining the
c         initial step bound. this bound is set to the product of
c         factor and the euclidean norm of diag*x if nonzero, or else
c         to factor itself. in most cases factor should lie in the
c         interval (.1,100.). 100. is a generally recommended value.
c
c       nprint is an integer input variable that enables controlled
c         printing of iterates if it is positive. in this case,
c         fcn is called with iflag = 0 at the beginning of the first
c         iteration and every nprint iterations thereafter and
c         immediately prior to return, with x and fvec available
c         for printing. if nprint is not positive, no special calls
c         of fcn with iflag = 0 are made.
c
c       info is an integer output variable. if the user has
c         terminated execution, info is set to the (negative)
c         value of iflag. see description of fcn. otherwise,
c         info is set as follows.
c
c         info = 0  improper input parameters.
c
c         info = 1  both actual and predicted relative reductions
c                   in the sum of squares are at most ftol.
c
c         info = 2  relative error between two consecutive iterates
c                   is at most xtol.
c
c         info = 3  conditions for info = 1 and info = 2 both hold.
c
c         info = 4  the cosine of the angle between fvec and any
c                   column of the jacobian is at most gtol in
c                   absolute value.
c
c         info = 5  number of calls to fcn has reached or
c                   exceeded maxfev.
c
c         info = 6  ftol is too small. no further reduction in
c                   the sum of squares is possible.
c
c         info = 7  xtol is too small. no further improvement in
c                   the approximate solution x is possible.
c
c         info = 8  gtol is too small. fvec is orthogonal to the
c                   columns of the jacobian to machine precision.
c
c       nfev is an integer output variable set to the number of
c         calls to fcn.
c
c       fjac is an output m by n array. the upper n by n submatrix
c         of fjac contains an upper triangular matrix r with
c         diagonal elements of nonincreasing magnitude such that
c
c                t     t           t
c               p *(jac *jac)*p = r *r,
c
c         where p is a permutation matrix and jac is the final
c         calculated jacobian. column j of p is column ipvt(j)
c         (see below) of the identity matrix. the lower trapezoidal
c         part of fjac contains information generated during
c         the computation of r.
c
c       ldfjac is a positive integer input variable not less than m
c         which specifies the leading dimension of the array fjac.
c
c       ipvt is an integer output array of length n. ipvt
c         defines a permutation matrix p such that jac*p = q*r,
c         where jac is the final calculated jacobian, q is
c         orthogonal (not stored), and r is upper triangular
c         with diagonal elements of nonincreasing magnitude.
c         column j of p is column ipvt(j) of the identity matrix.
c
c       qtf is an output array of length n which contains
c         the first n elements of the vector (q transpose)*fvec.
c
c       wa1, wa2, and wa3 are work arrays of length n.
c
c       wa4 is a work array of length m.
c
c     subprograms called
c
c       user-supplied ...... fcn
c
c       minpack-supplied ... dpmpar,enorm,fdjac2,lmpar,qrfac
c
c       fortran-supplied ... abs,max,min,sqrt,mod
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
c      INTERFACE
c#ifdef MPI_OPT
c      subroutine fdjac2_mp(fcn, fnorm,
c     1           iflag, ncnt, epsfcn, fvec_min, fnorm_min, x_min)
c      use kind_spec
c      use lmpar_mod
c      use fdjac_mod, ONLY: flip, ix_min, jac_order, jac_count
c      implicit none
c      include 'mpif.h'                                       !mpi stuff
c
c      external fcn
c      integer, intent(in) :: ncnt
c      integer, intent(out) :: iflag
c      real(rprec), intent(in) :: epsfcn, fnorm
c      real(rprec), intent(out) :: fnorm_min
c      real(rprec), intent(out), dimension(n) :: x_min
c      real(rprec), intent(out), dimension(m) :: fvec_min
c
c      end subroutine fdjac2_mp
c
c#else
c      subroutine fdjac2(fcn, fnorm, iflag, ncnt, epsfcn, time, 
c     1                  fnorm_min, x_min, fvec_min)
c      use kind_spec
c      use lmpar_mod
c      use fdjac_mod, eps1=>eps, ncnt1=>ncnt
c      implicit none
c
c      external fcn
c      integer, intent(in) :: ncnt
c      integer, target :: iflag
c      real(rprec) epsfcn, time, fnorm
c      real(rprec), intent(out) :: fnorm_min, x_min(n), fvec_min(m)
c      end subroutine fdjac2
c
c
c#endif
c      END INTERFACE



#ifdef MPI_OPT
c     Get mpi parameters and allocate iflag array
      call MPI_COMM_RANK (MPI_COMM_WORLD, myid, ierr)       !mpi stuff
      call MPI_COMM_SIZE (MPI_COMM_WORLD, numprocs, ierr)   !mpi stuff

      if (numprocs > (n+1)) then
         if (myid .eq. master) then
            write (6, *)'Warning: more processors have been requested',
     1      ' than the maximum required = ',n+1
         end if
      else if (numprocs < 1) then   
         if (myid .eq. master)
     1      write (6, *)'Must request at least ONE processor'
         goto 401
      end if
#endif

      allocate (x_min(n), fvec_min(m), flip(n), jac_order(n), 
     1          stat=istat)
      if (istat .ne. 0) stop 'Allocation error in lmdif'

!
!     epsmch is the machine precision.
!     flip is control for direction flipping in fdjac2!

      epsmch = dpmpar(1)
      flip = .false.;      jac_order = 0           
      info = 0;      iflag = 0;      nfev = 0;      cycle_count = 0
      lev_state = 0
      delta = 0

#ifndef MPI_OPT
      myid = 0;      wall_time = 0;      wall_time_lev = 0
#endif
!
!     ASSIGN MODULE POINTERS (FACILITATES PASSING TO SUBROUTINES)
!
c     ldfjac_mod = ldfjac
c     ipvt_mod => ipvt
c     fjac_mod => fjac
c     diag_mod => diag
c     qtf_mod => qtf
!      
!     check the input parameters for errors.
!
      if (n.le.0 .or. m.lt.n .or. ftol.lt.zero 
     1   .or. xtol.lt. zero .or. gtol.lt.zero .or. maxfev.le.0 
     2   .or. factor.le.zero) goto 400

      if (mode .eq. 2) then
         do j = 1, n
            if (diag(j) .le. zero) go to 300
         end do
      endif

!
!     evaluate the function at the starting point
!     and calculate its norm.
!
      if (myid .eq. master) then
         iflag = FLAG_SINGLETASK          
         call fcn (m, n, x, fvec, iflag, nfev)
      endif

      iflag = FLAG_CLEANUP                !!Clean-up and set everyone to result
      call fcn (m, n, x, fvec, iflag, nfev)

#ifdef MPI_OPT
      call MPI_BCAST(fvec, m, MPI_REAL8, master, MPI_COMM_WORLD, ierr)
      if (ierr .ne. 0) goto 3000
      call MPI_BCAST(iflag,1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      if (ierr .ne. 0) goto 3000
#endif

      if (iflag .lt. 0) go to 300
      fnorm = enorm(m,fvec)
      if (nfev.ge.maxfev .or. maxfev.eq.1) info = 5
      if (info .ne. 0) go to 300
!
!     initialize levenberg-marquardt parameter (par) and iteration counter.
!
      par = 0
      iter = 1
#ifdef MPI_OPT
      if (myid .eq. master) write (6, 1000) numprocs
 1000 format (/,' Beginning Levenberg-Marquardt Iterations',/,
     1        ' Number of Processors: ',i4,' (1 controller proc)',//,
     2        70('='),/,2x,'Iteration',3x,'Processor',7x,'Chi-Sq',7x,
     3       'LM Parameter',6x,'Delta Tol'/,70('='))
#else
      write (6, 1000) max_processors 
 1000 format (/,' Beginning Levenberg-Marquardt Iterations',/,
     1        ' Number processors requested: ', i4,//,
     1        59('='),/,2x,'Iteration',8x,'Chi-Sq',7x,
     2       'LM Parameter',6x,'Delta Tol',/,59('='))
#endif

!
!     beginning of the outer loop.
!
      first_jacobian = .true.

      outerloop: do while (nfev .lt. maxfev)
           delta_old = delta
           par_old = par
!
!        calculate the jacobian matrix.
!
#ifdef MPI_OPT
         call fdjac2_mp(fcn, fnorm, iflag, nfev, epsfcn, 
     1                  fvec_min, fnorm_min, x_min)
#else
         iflag = 2
         call fdjac2(fcn, fnorm, iflag, nfev, epsfcn, wall_time, 
     2               fnorm_min, x_min, fvec_min)
#endif     
         nfev = nfev + n
         if (iflag .lt. 0) exit outerloop

!
!        compute the qr factorization of the jacobian.
!
         call qrfac(m, n, fjac, m, .true., ipvt, 
     1              n, wa1, wa2, wa3)

!
!        on the first iteration and if mode is 1, scale according
!        to the norms of the columns of the initial jacobian.
!
         if (iter .eq. 1) then
            if (mode .ne. 2) then
               diag = wa2
               where (wa2 .eq. zero) diag = one
            endif

!
!        also on the first iteration, calculate the norm of the scaled x
!        and initialize the step bound delta.
!
            wa3 = diag*x
            xnorm = enorm(n,wa3)
            delta = factor*xnorm
            if (delta .eq. zero) delta = factor
         endif

!
!        form (q transpose)*fvec and store the first n components in qtf.
!
         wa4 = fvec
         do j = 1, n
            if (fjac(j,j) .ne. zero) then
               sum0 = sum(fjac(j:m,j)*wa4(j:m))
               temp = -sum0/fjac(j,j)
               wa4(j:m) = wa4(j:m) + fjac(j:m,j)*temp
            endif
            fjac(j,j) = wa1(j)
            qtf(j) = wa4(j)
         end do

!
!        compute the norm of the scaled gradient.
!
         gnorm = zero
         if (fnorm .ne. zero) then
            do j = 1, n
               l = ipvt(j)
               if (wa2(l) .ne. zero) then
                  sum0 = sum(fjac(:j,j)*(qtf(:j)/fnorm))
                  gnorm = max(gnorm,abs(sum0/wa2(l)))
               endif
            end do
         endif

!
!        test for convergence of the gradient norm.
!
         if (gnorm .le. gtol) info = 4
         if (info .ne. 0) exit outerloop

!
!        rescale if necessary.
!
         if (mode .ne. 2) diag = max(diag,wa2)

!
!        set up for inner loop (levmarqloop) to determine x update.
!
         subcycle = 0
         ratio = 0

         allocate (fjac_save(n,n), stat=istat)
         if (istat .ne. 0) stop 'Fjac_save allocation error'

         fjac_save(:n,:n) = fjac(:n,:n)

         levmarqloop: do while (ratio .lt. p0001)
 
           subcycle = subcycle + 1
           fjac(:n,:n) = fjac_save(:n,:n)
           first = (iter == 1) .and. (subcycle == 1)
           spread_ratio = abs(epsfcn)/delta/10

!        Determine the levenberg-marquardt parameter.  
#ifdef MPI_OPT
!        for parallel processing, scan a range of values for this
!        parameter to find the 'optimal' choice
!
           call levmarq_param_mp (nfev, iflag, fcn, lev_step_range)
#else
           call levmarq_param(wall_time_lev, nfev, iflag, fcn, 
     1                        lev_step_range)
#endif
           if (iflag .lt. 0) exit

!
!        on the first iteration, adjust the initial step bound.
!
c          if (iter .eq. 1) delta = min(delta, pnorm)

!
!        compute the scaled actual reduction.
!
           actred = -one
           if (p1*fnorm1 .lt. fnorm) actred = one - (fnorm1/fnorm)**2

!
!        compute the scaled predicted reduction (prered) and
!        the scaled directional derivative (dirder).
!
           do j = 1, n
              wa3(j) = zero
              l = ipvt(j)
              temp = wa1(l)
              wa3(:j) = wa3(:j) + fjac(:j,j)*temp
           end do
           temp1 = enorm(n,wa3)/fnorm
           temp2 = (sqrt(par)*pnorm)/fnorm
           prered = temp1**2 + temp2**2/p5
           dirder = -(temp1**2 + temp2**2)

!
!        compute the ratio of the actual to the predicted reduction.
!
           actred_lev = actred
           ratio = zero
           if (prered .ne. zero) ratio = actred/prered

!
!        test if Jacobian calculation gave best improvement
!
           step_improved = fnorm_min < min(fnorm,fnorm1)
           if ( step_improved ) then

!  Jacobian gave best improvement => having problems following the
!  Gradient adequately using the single-sided differences.   Hunt around 
!  a bit guided by the improvements found by the Jacobian calculation.

#ifdef MPI_OPT
              if( myid == master) 
     1           jac_count = min(count(jac_order(:) .ne. 0), 
     2                        int(sqrt(real(numprocs))))
!
!  Need to broadcast the jac_order array from the last jacobean calc.
!
              call MPI_BCAST(jac_count, 1, MPI_INTEGER, master,
     1                       MPI_COMM_WORLD, ierr)
              if (ierr .ne. 0) goto 3000

              if(myid .ne. master) jac_order = 0
              call MPI_BCAST(jac_order, jac_count, MPI_INTEGER, master,
     1                       MPI_COMM_WORLD, ierr)
              if (ierr .ne. 0) goto 3000

              call stepopt_mp(fcn, x_min, fvec_min, fnorm_min, iflag, 
     1                        nfev, epsfcn)
#else
              jac_count = min(count(jac_order(:) .ne. 0), 
     2                        int(sqrt(real(num_lm_params))))

              iflag = 2
              call stepopt(fcn, x_min, fvec_min, fnorm_min, iflag, nfev,  
     1                     epsfcn, wall_time_lev )
#endif     


              if (myid .eq. master) then
                 write(6,'(a,i6,a)') 
     1           ' Using minimum from Jacobian/step improvement (',
     2           ix_min, ')'
                 call vmec_flush(6)
              endif

              wa2 = x_min
              wa4 = fvec_min
              fnorm1 = fnorm_min
              if( cycle_count > 1) then
                 delta = max(delta, delta_old)
                 par = min(par, par_old)
              endif

              actred = 1 - (fnorm1/fnorm)**2
              ratio = actred/prered

           else if( first_jacobian .and. actred < 0) then
!
!  if there was no direction of improvement found
!  try recomputing jacobian with signs of displacements flipped
!
              first_jacobian = .false.
              deallocate (fjac_save)
              cycle outerloop

           else

!
!        update the step bound, if Levenberg Step gave best improvement
!

#ifdef MPI_OPT
             if( numprocs > 8 ) then
#else
             if( max_processors > 8 ) then
#endif 

!
!        In the multiprocessor case, if we have enough samples, defer to 
!        the best delta parameter found.  Because of the asymmetry in the
!        'factors' (in levmarq_param), need to goose delta when it is at the
!        top end of the range
!
               if( actred < 0) then    ! if no success this iteration

                 delta = min(delta, delta_old) / 4
                 par = par*4

               else if( lev_step_range == 1 ) then
                 delta = 2*delta
                 par = par/2
               endif

             else
               if (ratio .le. p25) then
                 if (actred .ge. zero) then
                    temp = p5
                 else
                    temp = p5*dirder/(dirder + p5*actred)
                 endif
                 if (p1*fnorm1.ge.fnorm .or. temp.lt.p1) temp = p1
                 delta = max( temp*min(delta,pnorm/p1), p1*delta)
                 par = par/temp
               else if (par.eq.zero .or. ratio.ge.p75) then
                 delta = max(pnorm/p5, p1*delta)
                 par = p5*par
               endif
             endif
           endif

!
!        test for successful iteration.
!        update x, fvec, and their norms.
!
#ifdef MPI_OPT
           if( n > numprocs .and. subcycle < 2 .and.
     1         ratio .lt. p0001 ) cycle levmarqloop
#else
           if( n > max_processors .and. subcycle < 2 .and.
     1         ratio .lt. p0001) cycle levmarqloop
#endif 
           if (ratio .ge. p0001) then
              x = wa2
              wa2 = diag*x
              fvec = wa4
              xnorm = enorm(n,wa2)
              fnorm = fnorm1
              iter = iter + 1
              if (myid .eq. master) then
                 write(6,'(/,i6,3(2x,a,es10.3)/)') 
     1            nfev, 'new minimum =', fnorm**2,'lm-par =', par, 
     2            'delta-tol =', delta
                 call vmec_flush(6)
              endif
           endif

           cycle_count = cycle_count + 1

!
!        tests for convergence.
!
           if (abs(actred).le.ftol .and. prered.le.ftol 
     1        .and. p5*ratio.le.one) info = 1
!
!        next test made more stringent to avoid premature declaration of 
!        completion observed on some problems.
!
           if (delta .le. xtol*xnorm/10 .and. cycle_count>2 
     1         .and. .not. step_improved ) info = 2
!          if (delta .le. xtol*xnorm) info = 2

           if (abs(actred).le.ftol .and. prered.le.ftol 
     1        .and. p5*ratio.le.one .and. info.eq.2) info = 3
           if (info .ne. 0) exit levmarqloop

!
!        tests for termination and stringent tolerances.
!
           if (nfev .ge. maxfev) info = 5
           if (abs(actred).le.epsmch .and. prered.le.epsmch 
     1        .and. p5*ratio.le.one) info = 6
           if (delta .le. epsmch*xnorm) info = 7
           if (gnorm .le. epsmch) info = 8
           if (info .ne. 0) exit levmarqloop
c
c        end of the inner loop. repeat if iteration unsuccessful.
c
         end do levmarqloop
 
         deallocate (fjac_save)
         if (info.ne.0 .or. iflag.ne.0) exit outerloop

         first_jacobian = .true.
      end do outerloop
c
c     termination, either normal or user imposed.
c
 300  continue

      if (iflag.eq.0 .and. info.eq.0) then
         info = 9
      else if (iflag .lt. 0) then
         info = iflag
      end if


 400  deallocate (x_min, fvec_min, flip, jac_order)

 401  if (myid .eq. master) then                                         ! MPI      
         if (info.ge.-1 .and. info.le.8) write (*, '(/,1x,a,/,1x,a)')
     1  'Levenberg-Marquardt optimizer status: ',trim(info_array(info))
      endif                                                              ! MPI

      if (nfev .le. 1) then
         nfev = 1
         iflag = FLAG_CLEANUP                !!Clean-up last time through for master
         call fcn (m, n, x, fvec, iflag, nfev)
      end if
         
#ifdef MPI_OPT
      return
 3000 continue
      print *, 'MPI_BCAST error in LMDIF, ierr = ', ierr
#else      
      write(*, '(2(/,a, f10.2))') 
     1     ' Total wall clock time in jacobian multi-process call  = ',
     2     wall_time,
     3     ' Total wall clock time in lev param multi-process call = ',
     4     wall_time_lev
#endif
      end subroutine lmdif

 
      function enorm (n, x) 
      use kind_spec
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer n 
      real(rprec), dimension(n) :: x
      real(rprec) :: enorm
!-----------------------------------------------
!
!     function enorm
!
!     given an n-vector x, this function calculates the
!     euclidean norm of x.
!
!
!     the function statement is
!
!       function enorm(n,x)
!
!     where
!
!       n is a positive integer input variable.
!
!       x is an input array of length n.
!
      enorm = sqrt (sum(x(:n)*x(:n)))

      end function enorm

#ifdef MPI_OPT
      subroutine fdjac2_mp(fcn, fnorm,
     1           iflag, ncnt, epsfcn, fvec_min, fnorm_min, x_min)
      use kind_spec
      use lmpar_mod
      use fdjac_mod, ONLY: flip, ix_min, jac_order, jac_count
      implicit none
      include 'mpif.h'                                       !mpi stuff
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: ncnt
      integer, intent(out) :: iflag
      real(rprec), intent(in) :: epsfcn, fnorm
c     real(rprec), dimension(n) :: x
c     real(rprec), dimension(m), intent(in) :: fvec
c     real(rprec), dimension(ldfjac,n), intent(out) :: fjac
      real(rprec), intent(out) :: fnorm_min
      real(rprec), intent(out), dimension(n) :: x_min
      real(rprec), intent(out), dimension(m) :: fvec_min

C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, parameter :: FLAG_CLEANUP = -100
      integer, parameter :: master = 0                       !mpi stuff
      integer :: j, ix_temp, nleft
      integer, dimension(1) :: isort
      integer, dimension(numprocs) :: iflag_array
      logical, dimension(n) :: lmask
      real(rprec) :: eps, epsmch, dpmpar, temp
      real(rprec), parameter :: one = 1, zero = 0
      real(rprec), dimension(n) :: h, fnorm_array            !mpi stuff
      real(rprec), dimension(m) :: buffer
      integer :: status(MPI_STATUS_SIZE)                     !mpi stuff
      integer :: i                                           !mpi stuff
      integer :: numsent, sender, ierr, istat                !mpi stuff
      integer :: anstype, column, ierr_flag                  !mpi stuff
      integer :: ic1, ic2, ic3, irate, count_max
      real(rprec) :: time, time2
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      external fcn, dpmpar, enorm
      real(rprec) enorm
C-----------------------------------------------
c
c     subroutine fdjac2
c
c     this subroutine computes a forward-difference approximation
c     to the m by n jacobian matrix associated with a specified
c     problem of m functions in n variables.
c
c     the subroutine statement is
c
c       subroutine fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)
c
c     where
c
c       fcn is the name of the user-supplied subroutine which
c         calculates the functions. fcn must be declared
c         in an external statement in the user calling
c         program, and should be written as follows.
c
c         subroutine fcn(m,n,x,fvec,iflag)
c         integer m,n,iflag
c         real(rprec) x(n),fvec(m)
c         ----------
c         calculate the functions at x and
c         return this vector in fvec.
c         ----------
c         return
c         end
c
c         the value of iflag should not be changed by fcn unless
c         the user wants to terminate execution of fdjac2.
c         in this case set iflag to a negative integer.
c
c       m is a positive integer input variable set to the number
c         of functions.
c
c       n is a positive integer input variable set to the number
c         of variables. n must not exceed m.
c
c       x is an input array of length n.
c
c       fvec is an input array of length m which must contain the
c         functions evaluated at x.
c
c       fjac is an output m by n array which contains the
c         approximation to the jacobian matrix evaluated at x.
c
c       ldfjac is a positive integer input variable not less than m
c         which specifies the leading dimension of the array fjac.
c
c       iflag is an integer variable which can be used to terminate
c         the execution of fdjac2. see description of fcn.
c
c       epsfcn is an input variable used in determining a suitable
c         step length for the forward-difference approximation. this
c         approximation assumes that the relative errors in the
c         functions are of the order of epsfcn. if epsfcn is less
c         than the machine precision, it is assumed that the relative
c         errors in the functions are of the order of the machine
c         precision.
c
c       wa is a work array of length m.
c
c     subprograms called
c
c       user-supplied ...... fcn
c
c       minpack-supplied ... dpmpar
c
c       fortran-supplied ... abs,max,sqrt
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********

      ierr_flag = 0
c
c     epsmch is the machine precision.
c
      epsmch = dpmpar(1)
c
      eps = sqrt(max(abs(epsfcn),epsmch))
      if( epsfcn < 0 ) eps = - eps

c     call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )       !mpi stuff
c     call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )   !mpi stuff
c     call MPI_BARRIER(MPI_COMM_WORLD, ierr)                 !mpi stuff
c
c     Begin parallel MPI master/worker version of the function call.  A
c     bank queue or master/worker MPI algorithm is used here to achieve
c     parallel load balancing in the somewhat uneven work load involved
c     in calculating the Jacobian function needed in the Levenberg-Marquardt
c     optimization.  This model is based on the example given in Chapt. 5 of
c     "The LAM Companion to Using MPI" by Zdzislaw Meglicki (online see:
c     http://carpanta.dc.fi.udc.es/docs/mpi/mpi-lam/mpi.html).  These
c     modifications were made by D.A. Spong 8/23/2000.
c
c     ****Master portion of the code****
c
      if (myid .eq. master) then
         numsent = 0    !numsent is a counter used to track how many
                        !jobs have been sent to workers
         fnorm_min = 0
c        cur_norm = enorm(m,fvec)
c        cur_norm = fnorm
c
c        calculate the displacements
c
         h(1:n) = eps*abs(x(1:n))
         where (h == zero) h=eps
         where (flip) h = -h

c     Send forward difference displacements from master to each
c           worker process. Tag with these with the column number.
c
         nleft = n

         do j = 1,min(numprocs-1,n)
            temp = x(j)
            x(j) = temp + h(j)
            call MPI_SEND(x, n, MPI_REAL8, j, 
     1                  j, MPI_COMM_WORLD, ierr)
            if (ierr .ne. 0) stop 'MPI_SEND error(1) in fdjac2' 
            x(j) = temp
            numsent = numsent+1
         end do          !j = 1,min(numprocs-1,n)

c      should we do one of the cases? 
c      Only if needed and will make a difference in # cycles

         if( numsent < n .and. mod(n-1, numprocs-1) == 0) then
            numsent = numsent+1
            temp = x(numsent)
            x(numsent) = temp + h(numsent)
               
c           Call the chisq fcn 
            iflag = numsent
            call fcn(m, n, x, wa4, iflag, ncnt)
            ierr_flag = iflag

            x(numsent) = temp
            fjac(:m,numsent) = (wa4(:m) - fvec(:m))/h(numsent)
!
!           STORE FNORM OF PERTURBED STATE (X + H)
!
            temp = enorm(m,wa4)
            fnorm_array(numsent) = temp
            if (temp > fnorm) flip(numsent) = .not. flip(numsent)

            write (6, '(2x,i6,8x,i3,7x,1es12.4)') ncnt+numsent,
     1             myid, temp**2

            nleft = nleft - 1
         endif

c
c      Looping through the columns, collect answers from the workers.
c      As answers are received, new uncalculated columns are sent
c      out to these same workers.
c
         do j = 1,nleft
            call MPI_RECV(wa4, m, MPI_REAL8, 
     1           MPI_ANY_SOURCE, MPI_ANY_TAG, 
     2           MPI_COMM_WORLD, status, ierr)
            if (ierr .ne. 0) stop 'MPI_RECV error(1) in fdjac2'
            sender     = status(MPI_SOURCE)    
            anstype    = status(MPI_TAG)       ! column is tag value
            if (anstype .gt. n) stop 'ANSTYPE > N IN FDJAC2'
            
            fjac(:m,anstype) = (wa4(:m) - fvec(:m))/h(anstype)
!
!           STORE FNORM OF PERTURBED STATE (X + H)
!
            temp = enorm(m,wa4)
            fnorm_array(anstype) = temp
            if (temp > fnorm) flip(anstype) = .not. flip(anstype)

            write (6, '(2x,i6,8x,i3,7x,1es12.4)') ncnt+anstype,
     1             sender, temp**2
c
c           If more columns are left, then send another column to the worker(sender)
c           that just sent in an answer
c
            if (numsent .lt. n) then
               numsent = numsent+1
               temp = x(numsent)
               x(numsent) = temp + h(numsent)
               
               call MPI_SEND(x, n, MPI_REAL8, 
     1                       sender, numsent, MPI_COMM_WORLD, ierr)
               if (ierr .ne. 0) stop 'MPI_SEND error(2) in fdjac2'
               x(numsent) = temp

            else                ! Tell worker that there is no more work to do
               
               call MPI_SEND(MPI_BOTTOM, 0, MPI_REAL8,
     1                       sender, 0, MPI_COMM_WORLD, ierr)
               if (ierr .ne. 0) stop 'MPI_END error(3) in fdjac2'
            endif      ! if( myid .eq. master ) then
         end do     ! do j = 1,n

c         ierr_flag = 0
         call vmec_flush(6)
c
c     ****Worker portion of the code****
c        Skip this when processor id exceeds work to be done
c
      else if (myid .le. n) then        ! i.e., if( myid .ne. master )
c
c        Otherwise accept the next available column, check the tag,
c        and if the tag is non-zero call subroutine fcn.
c        If the tag is zero, there are no more columns
c        and worker skips to the end.
c
 90      call MPI_RECV(x, n, MPI_REAL8, master, 
     1                 MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
         if (ierr .ne. 0) stop 'MPI_RECV error(2) in fdjac2'
        
         column = status(MPI_TAG)
         if (column .eq. 0) then
            go to 200
         else
            iflag = column
c           Call the chisq fcn for the portion of displacement vector which
c           was just received. Note that WA4 stores the local fvec_min array
            call fcn(m, n, x, wa4, iflag, ncnt)

            if (iflag.ne.0 .and. ierr_flag.eq.0) ierr_flag = iflag
c
c           Send this function evaluation back to the master process tagged
c           with the column number so the master knows where to put it
c
            call MPI_SEND(wa4, m, MPI_REAL8, master, 
     1                    column, MPI_COMM_WORLD, ierr)
            if (ierr .ne. 0) stop 'MPI_SEND error(4) in fdjac2'
            go to 90    !Return to 90 and check if master process has sent any more jobs
         endif
 200     continue

      else        ! this is an extra processor
         ierr_flag = 0             ! make sure

      endif       ! if( myid .ne. master )


!
!     Gather iflag information to ALL processors and check for iflag < 0
!
      call MPI_ALLGATHER(ierr_flag, 1, MPI_INTEGER, iflag_array, 1,
     1     MPI_INTEGER, MPI_COMM_WORLD, ierr)
      if (ierr .ne. 0) stop 'MPI_ALLGATHER failed in FDJAC2'

      iflag = minval(iflag_array)
      if (iflag .lt. 0) return


!
!     Broadcast the fjac matrix
!
c     do j=1,n
c        if (myid .eq. master) buffer(:m) = fjac(:m,j)
c        call MPI_BCAST(buffer, m, MPI_REAL8, master,
c    1        MPI_COMM_WORLD, ierr)
c        if (ierr .ne. 0) go to 100
c        if (myid .ne. master) fjac(:m,j) = buffer(:m)
c     end do   


      call MPI_BCAST(fjac, m*n, MPI_REAL8, master,
     1        MPI_COMM_WORLD, ierr)
      if (ierr .ne. 0) go to 100

!
!     Find processor with minimum fnorm_min value and broadcast fvec_min, x_min, fnorm_min
!
      if (myid .eq. master) then
         jac_order = 0
         ix_temp = 1
         temp = 0
         lmask = .true.

         do while (ix_temp <= n)
            isort = minloc(fnorm_array, MASK=lmask)

            temp = fnorm_array(isort(1))
            jac_order(ix_temp) = isort(1)

            if(isort(1) <= 0 .or. isort(1) > n) then
               exit
            else if(fnorm_array(isort(1)) > fnorm) then
               exit
            else
               lmask(isort(1)) = .false.
               ix_temp = ix_temp + 1
            endif
         enddo

         ix_min = jac_order(1)
         if(ix_temp <= n) jac_order(ix_temp:) = 0
         jac_count = ix_temp - 1

         print *,jac_order(1:jac_count)

         if (ix_min .le. 0 .or. ix_min .gt. n) then
            print *,' IX_MIN = ',ix_min,' out of range'
            stop
         end if
         fnorm_min = fnorm_array(ix_min)
         fvec_min(:) = fjac(:, ix_min)*h(ix_min) + fvec(:)
         x_min(:) = x(:) 
         x_min(ix_min) = x(ix_min) + h(ix_min)
      end if

      call MPI_BCAST(ix_min, 1, MPI_INTEGER, master,
     1               MPI_COMM_WORLD, ierr)
      call MPI_BCAST(fnorm_min, 1, MPI_REAL8, master,
     1               MPI_COMM_WORLD, ierr)
      if (ierr .ne. 0) go to 100
      call MPI_BCAST(fvec_min, m, MPI_REAL8, master, MPI_COMM_WORLD, 
     1               ierr)
      if (ierr .ne. 0) go to 100
      call MPI_BCAST(x_min, n, MPI_REAL8, master, MPI_COMM_WORLD, 
     1               ierr)
      if (ierr .ne. 0) go to 100 
c
c     Original serial version of the Jacobian calculation:
c
c      do j = 1, n
c         temp = x(j)
c         h = eps*abs(temp)
c         if (h .eq. zero) h = eps
c         x(j) = temp + h
c         call fcn (m, n, x, wa, iflag)
c         if (iflag .lt. 0) exit 
c         x(j) = temp
c         fjac(:m,j) = (wa - fvec)/h
c      end do

!
!     Do any special cleanup now for IFLAG = FLAG_CLEANUP
!
      column = FLAG_CLEANUP
      call fcn(m, n, x, wa4, column, ncnt)

!
!     Reassign initial x value to all processors and perform error handling
!     IS THIS NECESSARY? CHECK...         
      call MPI_BCAST(x, n, MPI_REAL8, master, MPI_COMM_WORLD, ierr)
      if (ierr .ne. 0) goto 100

      
      return

 100  continue
      print *,' MPI_BCAST error in FDJAC2_MP: IERR=', ierr

      end subroutine fdjac2_mp


      subroutine levmarq_param_mp(nfev, iflag, fcn, lev_step_range)
      use lmpar_mod
      implicit none
      include 'mpif.h'                                       !mpi stuff
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: nfev, iflag, lev_step_range
c      real(rprec) :: x(n), wa1(n), wa2(n), wa3(n), wa4(m)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: FLAG_CLEANUP = -100, master = 0
      real(rprec), parameter :: zero = 0
      real(rprec), dimension(11), parameter :: factors =
     1  (/ 1.0_dp, 0.5_dp, 0.25_dp, 0.128_dp, 2.1_dp, 0.75_dp,
     2      1.25_dp, 1.5_dp, 0.9_dp, 1.1_dp, 1.75_dp /)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: iproc, iproc_min, nfact, num_lev, istat, ierr, j, len
      integer, dimension(numprocs) :: iflag_array
      real(rprec) :: scale_factor
      real(rprec), dimension(numprocs) :: 
     1                           fnorm_array, delta_array, par_array
      character*1 :: ext, low_mark
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      external fcn
      real(rprec), external :: enorm
C-----------------------------------------------      
!     Perform Numprocs different function calls (in parallel) find the minimum norm (chi-sq).
!     MPI calls are used to determine which processor has the minimum norm and then the 
!     information is sent to all processors using MPI Broadcasts (modifications made by D. A. Spong 8/27/2000).

      nfact = size(factors)
      num_lev = numprocs
      iproc = myid+1

      if( first .and. num_lev > 2) then            
c
c       do an exponential spread the first time to see where we are
c
         scale_factor = exp((iproc-1)*log(spread_ratio)/num_lev)
      else if (num_lev > 2*nfact) then
         scale_factor = (iproc*maxval(factors))/num_lev
      else if (iproc .le. nfact) then
         scale_factor = factors(iproc)
      else
         scale_factor = ((iproc-nfact)*minval(factors)*0.7)/
     1                   (num_lev-nfact)
      end if

      delta = scale_factor*delta

      call lmpar (n, fjac, m, ipvt, diag, qtf,
     1            delta, par, wa1, wa2, wa3, wa4)
!
!     store the direction p and x + p. calculate the norm of p.
!
      if (par .eq. zero) wa1 = wa1*scale_factor
      wa1 = -wa1
c        wa1 = -wa1*scale_factor
      wa2 = x + wa1
      wa3 = diag*wa1
      pnorm = enorm(n,wa3)

c     delta = scale_factor*delta
c     par = par / scale_factor
!
!     evaluate the function at x + p and calculate its norm.
!     Only do for 0 <= myid < n processors to avoid clean-up problems (in lsfun1)
!

      iflag = iproc
      call fcn (m, n, wa2, wa4, iflag, nfev)
      fnorm1 = enorm(m,wa4)

!
!     Gather iflag information to ALL processors and check for iflag < 0
!
      call MPI_ALLGATHER(iflag, 1, MPI_INTEGER, iflag_array, 1,
     1     MPI_INTEGER, MPI_COMM_WORLD, ierr)
      if (ierr .ne. 0) stop 'MPI_ALLGATHER failed in LMDIF'

      iflag = minval(iflag_array)
      if (iflag .lt. 0) return


!
!     Find processor with minimum fnorm1 value
!
      call MPI_ALLGATHER(fnorm1, 1, MPI_REAL8, fnorm_array, 1,
     1     MPI_REAL8, MPI_COMM_WORLD, ierr)
      if (ierr .ne. 0) stop 'MPI_ALLGATHER of fnorm failed in LMDIF'           
      iflag_array(1:1) = minloc(fnorm_array) 
      iproc_min = iflag_array(1) - 1
      call MPI_ALLGATHER(delta, 1, MPI_REAL8, delta_array, 1,
     1     MPI_REAL8, MPI_COMM_WORLD, ierr)
      if (ierr .ne. 0) stop 'MPI_ALLGATHER of delta failed in LMDIF'           

      call MPI_ALLGATHER(par, 1, MPI_REAL8, par_array, 1,
     1     MPI_REAL8, MPI_COMM_WORLD, ierr)

      if( myid .eq. master) then
         do j=1, num_lev
            ext = ' '
            if (j == 0) ext = '*'
            low_mark = ' '
            if (j == iproc_min+1) low_mark = '*'

            write(6, '(2x,i6,8x,i3,4x,2(3x,es12.4,a),(3x,es12.4))') 
     1         j+nfev, j, fnorm_array(j)**2, low_mark, 
     2         par_array(j), ext, delta_array(j)
         enddo
         write(6, '(a)') '  '

         call vmec_flush(6)
      endif

      fnorm1 = fnorm_array(iproc_min+1)
      delta = delta_array(iproc_min+1)
      par = par_array(iproc_min+1)

      lev_step_range = 0
      iflag_array(1:1) = maxloc(delta_array)
      if( iproc_min == iflag_array(1) - 1 ) then
         lev_step_range = 1
      else if( iproc_min == iflag_array(numprocs) - 1 ) then
         lev_step_range = -1
      endif

!
!     Broadcast all relevant scalars and arrays from the
!     processor with minimum fnorm1 to the other processors,
!     overwriting their data. Note: diag, ipvt are same already on
!     all processors. wa3 is overwritten...
!
      call MPI_BCAST(pnorm,1,MPI_REAL8,iproc_min,
     1     MPI_COMM_WORLD,ierr)
      if (ierr .ne. 0) goto 3000

c     call MPI_BCAST(par,1,MPI_REAL8,iproc_min,
c    1     MPI_COMM_WORLD,ierr)
c     if (ierr .ne. 0) goto 3000
c     call MPI_BCAST(delta,1,MPI_REAL8,iproc_min,
c    1     MPI_COMM_WORLD,ierr)
c     if (ierr .ne. 0) goto 3000

      call MPI_BCAST(wa1,n,MPI_REAL8,iproc_min,
     1     MPI_COMM_WORLD,ierr)
      if (ierr .ne. 0) goto 3000
      call MPI_BCAST(wa2,n,MPI_REAL8,iproc_min,
     1     MPI_COMM_WORLD,ierr)
      if (ierr .ne. 0) goto 3000
      call MPI_BCAST(wa4,m,MPI_REAL8,iproc_min,
     1     MPI_COMM_WORLD,ierr)
      if (ierr .ne. 0) goto 3000

c   this is not actually needed...mcz
!
!     BROADCAST JACOBIAN fjac(:n,j), j=1,n ROW BY ROW (CHANGED IN LMPAR) TO OTHER PROCESSES
!
c     do j = 1, n
c        if (myid .eq. iproc_min) wa3(:n) = fjac(:n,j)
c        call MPI_BCAST(wa3, n, MPI_REAL8, iproc_min,
c    1        MPI_COMM_WORLD, ierr)
c        if (ierr .ne. 0) goto 3000
c        if (myid .ne. iproc_min) fjac(:n,j) = wa3(:n)
c     end do

!
!     CLEANUP AFTER LEVENBERG-MARQUARDT LOOP AS NEEDED (WA4 IS NOT CHANGED)
!

2999  iflag = FLAG_CLEANUP                        
      call fcn (m, n, x, wa4, iflag, nfev)             !Contains Bcast Barrier

      nfev = nfev + num_lev

      return
      
 3000 continue
       
      print *, 'MPI_BCAST error in LEVMARQ_PARAM_MP, ierr = ', ierr
 
      end subroutine levmarq_param_mp



      subroutine stepopt_mp(fcn, x, fvec, fnorm, iflag, nfev, epsfcn )

      use kind_spec
      use fdjac_mod, ONLY: flip, ix_min, jac_order, jac_count
      use lmpar_mod, ONLY: m, n, wa2, wa4, myid, numprocs
      implicit none
      include 'mpif.h'                                       !mpi stuff
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
c     integer, intent(in) :: m, n, myid, numprocs
      integer :: nfev
      integer, intent(out) :: iflag
      real(rprec), intent(in) :: epsfcn
      real(rprec) :: fnorm
      real(rprec), dimension(n) :: x
      real(rprec), dimension(m) :: fvec
c     real(rprec) :: wa2(n), wa4(m)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, parameter :: FLAG_CLEANUP = -100
      integer, parameter :: master = 0                       !mpi stuff
      integer :: i, j, k, l, istat, iread, ic1, ic2, irate, count_max,
     1     jmin, ll, kk, mm, jj, ierr, iproc_min, num_lev  
      integer, dimension(numprocs) :: iflag_array
      real(rprec) :: fnorm1, temp, h
      real(rprec), dimension(numprocs) :: fnorm_array
      character*1 :: ext, low_mark
C-----------------------------------------------
C   L o c a l   Parameters
C-----------------------------------------------
      integer, parameter :: skip = 4
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      external fcn, dpmpar, enorm
      real(rprec) :: enorm
C-----------------------------------------------
c
c     subroutine stepopt
c
c     this subroutine uses the ordered array of jacobian-evaluated 
c     improvements to attempt ot find incremental improvements in fnorm
c     in the vicinity of the previous best case.  This routine is called
c     when the optimizer is having difficulty interpretting the local
c     gradient, as syptified by the levenberg-marquardt step not giving the
c     best improvement.  Rather than immediately going for another Jacobian 
c     evaluation, we use the information from the last one to hunt a bit.
c
c     M. Zarnstorff  March 2002
C-----------------------------------------------

      j = myid 
      num_lev = numprocs-1

      if(j .gt. 0 ) then

!
!     propagate the ordered list of directions to bump
!


!
!     store the direction p and x + p. calculate the norm of p.
!
!     l = ((j-1) / jac_count)*skip + 1

         l = ((j-1) / jac_count)
         if( l <= skip/2 ) then
            l = l+1
         else
            l = (l-skip/2) * skip + 1
         endif

         k = mod((j-1), jac_count) + 1
         if( k == 1) l = l+1
         
         wa2 = x
         
         do i=1, k
            jj = jac_order(i)
            temp = wa2(jj)
            h = epsfcn * abs(temp)
            if (h .eq. 0 ) h = epsfcn
            if( flip(jj)) h = -h
            if( i .ne. 1) then
               wa2(jj) = temp + l*h
            else
               wa2(jj) = temp + (l-1)*h
            endif
         enddo
      
c
c     evaluate the function at x + p and calculate its norm.
c
c     if(j .le. num_lev) then
         iflag = j
         call fcn (m, n, wa2, wa4, iflag, nfev)
         fnorm1 = enorm(m, wa4)
      else
         iflag = 0
         fnorm1 = huge(fnorm1)
      endif
!
!     Gather iflag information to ALL processors and check for iflag < 0
!
      call MPI_ALLGATHER(iflag, 1, MPI_INTEGER, iflag_array, 1,
     1     MPI_INTEGER, MPI_COMM_WORLD, ierr)
      if (ierr .ne. 0) stop 'MPI_ALLGATHER failed in STEPOPT_MP'

      iflag = minval(iflag_array)
      if (iflag .lt. 0) return

!
!     Find processor with minimum fnorm1 value
!
      call MPI_ALLGATHER(fnorm1, 1, MPI_REAL8, fnorm_array, 1,
     1     MPI_REAL8, MPI_COMM_WORLD, ierr)
      if (ierr .ne. 0) stop 'MPI_ALLGATHER failed in LMDIF'           
      iflag_array(1:1) = minloc(fnorm_array) 
      iproc_min = iflag_array(1) - 1
      fnorm1 = fnorm_array(iproc_min+1)

      if( myid == master) then
         jmin = iproc_min
         low_mark = ' '
         if (fnorm < fnorm1) then
            low_mark = '*'
            jmin = numprocs+1
         endif


        write (6, '(2x,6x,4x,1es12.4,a,3x,a,i4,2x,i4,2x,i4,a)')
     1      fnorm**2, low_mark, '(',1,1,jac_order(1), ')' 
        do j=1, num_lev
!           ll = ((j-1) / jac_count)*skip + 1

           ll = ((j-1) / jac_count)
           if( ll <= skip/2 ) then
              ll = ll+1
           else
              ll = (ll-skip/2) * skip + 1
           endif

           kk = mod((j-1), jac_count) + 1
           if( kk == 1) ll = ll+1
           jj = jac_order(kk)
           low_mark = ' '
           if( j == jmin) low_mark = '*'
           
           write(6, '(2x,i6,4x,1es12.4,a,3x,a,i4,2x,i4,2x,i4,a)')
     1           j+nfev, fnorm_array(j+1)**2, low_mark, 
     2           '(', ll,kk,jj, ')'
        enddo
      endif

      if( fnorm1 < fnorm) then
         fnorm = fnorm1

!
!     Broadcast all relevant scalars and arrays from the
!     processor with minimum fnorm1 to the other processors,
!     overwriting their data. Note: diag, ipvt are same already on
!     all processors. wa3 is overwritten...
!
         call MPI_BCAST(wa2,n,MPI_REAL8,iproc_min,
     1     MPI_COMM_WORLD,ierr)
         if (ierr .ne. 0) goto 3000
         call MPI_BCAST(wa4,m,MPI_REAL8,iproc_min,
     1     MPI_COMM_WORLD,ierr)
         if (ierr .ne. 0) goto 3000

         x = wa2
         fvec = wa4
      endif

!
!     CLEANUP 
!
      iflag = FLAG_CLEANUP                        
      call fcn (m, n, x, wa4, iflag, nfev)             !Contains Bcast Barrier

      nfev = nfev + num_lev

      return
      
 3000 continue
       
      print *, 'MPI_BCAST error in STEPOPT_MP, ierr = ', ierr
 
      end subroutine stepopt_mp

#else
      subroutine fdjac2(fcn, fnorm, iflag, ncnt, epsfcn, time, 
     1                  fnorm_min, x_min, fvec_min)
      use kind_spec
      use lmpar_mod
      use fdjac_mod, eps1=>eps, ncnt1=>ncnt
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: ncnt
      integer, target :: iflag
      real(rprec) epsfcn, time, fnorm
      real(rprec), intent(out) :: fnorm_min, x_min(n), fvec_min(m)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: FLAG_CLEANUP = -100
      real(rprec), parameter :: zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: j, k, istat, iread, ic1, ic2, irate, count_max,
     1           ix_temp, ix_temp2
      integer, dimension(1) :: isort
      logical, dimension(n) :: lmask
      real(rprec) :: eps, epsmch, h, dpmpar, temp
      real(rprec), dimension(n) :: fnorm_array
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      external fcn, dpmpar, multiprocess, fdjac_parallel
      real(rprec), external :: enorm
C-----------------------------------------------
c
c     subroutine fdjac2
c
c     this subroutine computes a forward-difference approximation
c     to the m by n jacobian matrix associated with a specified
c     problem of m functions in n variables.
c
c     Here
c
c       fcn is the name of the user-supplied subroutine which
c         calculates the functions. fcn must be declared
c         in an external statement in the user calling
c         program (see LMDIF1 for documentation), and should be written as follows:
c
c         subroutine fcn(m,n,x,fvec,iflag,ncnt)
c         integer m,n,iflag
c         real(rprec) x(n),fvec(m)
c         ----------
c         calculate the functions at x and
c         return this vector in fvec.
c         ----------
c         return
c         end
c
c       fjac is an output m by n array which contains the
c         approximation to the jacobian matrix evaluated at x.
c
c       ldfjac is a positive integer input variable not less than m
c         which specifies the leading dimension of the array fjac.
c
c       iflag is an integer variable which can be used to terminate
c         the execution of fdjac2. see description of fcn.
c
c       epsfcn is an input variable used in determining a suitable
c         step length for the forward-difference approximation. this
c         approximation assumes that the relative errors in the
c         functions are of the order of epsfcn. if epsfcn is less
c         than the machine precision, it is assumed that the relative
c         errors in the functions are of the order of the machine
c         precision.
c
c       wa is a work array of length m.
c
c     subprograms called
c
c       user-supplied ...... fcn
c
c       minpack-supplied ... dpmpar
c
c       fortran-supplied ... abs,max,sqrt
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
 
c
c     epsmch is the machine precision.
c
      epsmch = dpmpar(1)
c
      eps = sqrt(max(abs(epsfcn),epsmch))
      if( epsfcn < 0 ) eps = -eps
!
!     Load module values. Pointers will automatically update targets...
!     Prepare for multi-processing...
!

      ncnt1 = ncnt
      eps1 = eps

!     Find min chisq = fnorm**2 state for this jacobian evaluation
!     (Do NOT retain from previous evaluations, or could get into a non-converging loop...)
      fnorm_min = -1                                                     
      ix_min = 0
            
      call system_clock(ic1, irate)
      call multiprocess(n, max_processors, fdjac_parallel, fcn)
      call system_clock(ic2, irate, count_max)
      if (ic2 .lt. ic1) ic2 = ic2 + count_max
      
c     do j = 1, n
c         temp = x(j)
c         h = eps*abs(temp)
c         if (h .eq. zero) h = eps
c         x(j) = temp + h
c         call fcn (m, n, x, wa4, iflag, ncnt)
c         if (iflag .lt. 0) exit 
c         x(j) = temp
c         fjac(:m,j) = (wa4 - fvec)/h
c      end do


c     cur_norm = enorm(m,fvec)      ! where are we now?
c     cur_norm = fnorm

      do j = 1, n

        read (j+1000, iostat=iread) istat, iflag, h, temp
        if (iread .ne. 0) then           
           write (6, *) 'Error reading from file fort.', j+1000,
     1     ' in fdjac2: IOSTAT = ', iread           
           iflag = -14
        else if (j .ne. istat) then
           write (6, *) 'Wrong value for index j read in fdjac2'
           iflag = -14
        end if   
          
        if (iflag .ne. 0) exit
        
#ifdef CRAY
        do k = 1, m
           read (j+1000) wa4(k)
        end do   
#else
        read (j+1000) wa4
#endif        
        fjac(:m,j) = (wa4 - fvec)/h

        fnorm_array(j) = temp
        if( temp > fnorm) flip(j) = .not. flip(j)  ! flip for next time

        if( fnorm_min < 0) fnorm_min = temp * 1.1_dp
        if( temp < fnorm_min) then
           fnorm_min = temp
           ix_min = j
           fvec_min = wa4
#ifdef CRAY
           do k = 1, n
              read (j+1000) x_min(k)
           end do   
#else
           read (j+1000) x_min
#endif        
        endif

        close (j+1000, status='delete')                        !!Needed to run correctly in multi-tasking...

      end do
      
      jac_order = 0
      ix_temp = 1
      temp = 0
      lmask = .true.
        
      do while (ix_temp <= n)
         isort = minloc(fnorm_array, MASK=lmask)

         temp = fnorm_array(isort(1))
         jac_order(ix_temp) = isort(1)

         if(isort(1) <= 0 .or. isort(1) > n) then
            exit
         else if(fnorm_array(isort(1)) > fnorm) then
            exit
         else
            lmask(isort(1)) = .false.
            ix_temp = ix_temp + 1
         endif
      enddo

      ix_min = jac_order(1)
      if(ix_temp <= n) jac_order(ix_temp:) = 0
      jac_count = ix_temp - 1

      print *,jac_order(1:jac_count)

      if (ix_min .le. 0 .or. ix_min .gt. n) then
         print *,' IX_MIN = ',ix_min,' out of range'
         stop
      end if
      fnorm_min = fnorm_array(ix_min)
!
!     Do any special cleanup now for IFLAG = -100
!
      iflag = FLAG_CLEANUP
      call fcn(m, n, x, wa4, iflag, ncnt)

      time = time + real(ic2 - ic1)/real(irate)                !!Time in multi-process call
      
      end subroutine fdjac2

 
      subroutine fdjac_parallel(j, fcn)
      use lmpar_mod
      use fdjac_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: j
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer istat, iflag, i
      real(rprec) :: temp, temp2, h, enorm
      external fcn, enorm
C-----------------------------------------------      
!
!     THIS ROUTINE IS PASSED TO THE MULTI-PROCESSOR HANDLING
!     ROUTINE


      temp = x(j)
      h = eps*abs(temp)
      if (h .eq. zero) h = eps
      if( flip(j)) h = -h
      x(j) = temp + h
      iflag = j

      call fcn (m, n, x, wa4, iflag, ncnt)

      temp2 = enorm(m, wa4)
      write(6, '(2x,i6,7x,1es12.4)') ncnt+j, temp2**2
      
!
!     WRITE TO A UNIQUE FILE FOR I/O IN MULTI-PROCESSOR SYSTEM
!
      write (j+1000) j, iflag, h, temp2
#ifdef CRAY
      do i = 1,m
         write (j+1000) wa4(i)
      end do   
      do i = 1,n
         write (j+1000) x(i)
      end do   
#else
      write (j+1000) wa4
      write (j+1000) x
#endif                  
      close (j+1000)                      !!Needed to run correctly in multi-tasking...

      x(j) = temp

      end subroutine fdjac_parallel


      subroutine levmarq_param(time, nfev, iflag, fcn, lev_step_range)
      use fdjac_mod
      use lmpar_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: nfev, iflag, lev_step_range
      real(rprec) :: time
c     real(rprec), target :: x(n), wa1(n), wa2(n), wa3(n),
c    1    wa4(m)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: FLAG_CLEANUP = -100
      real(rprec), parameter :: zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: j, k, istat, iread, ic1, ic2, irate, count_max,
     1     jmin      
      real(rprec), dimension(num_lm_params) :: 
     1      fnorm_min, pnorm_min, delta_min, par_min
      character*1 :: ext, low_mark
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      external fcn, multiprocess, lmpar_parallel
C-----------------------------------------------      
!
!     Initialize module variables for use in lmpar_parallel
!
      
c     xp => x
c     wap => wa1
c     wa2p => wa2
c     wa3p => wa3
c     wa4p => wa4
c     np = n
c     mp = m
      ncnt = nfev    

      call system_clock(ic1, irate)

      call multiprocess(num_lm_params, max_processors, 
     1    lmpar_parallel, fcn)

      call system_clock(ic2, irate, count_max)
      if (ic2 .lt. ic1) ic2 = ic2 + count_max

      nfev = nfev + num_lm_params

!
!     Read in optimal wa1, wa2, wa4, par, delta, pnorm, fnorm1 value from file
!
      do j = 1, num_lm_params

        read (j+1000, iostat=iread) istat, iflag, pnorm_min(j), 
     1        fnorm_min(j), par_min(j), delta_min(j)
        if (iread .ne. 0) then
           write (6, *) 'Error reading from file fort.', j+1000,
     1       ' in levmarq_param', ' IOSTAT = ', iread           
           iflag = -15
        else if (j .ne. istat) then
           write (6, *)
     1        'Incorrect value read in for index j in levmarq_param'
           iflag = -15
        end if    

        if (iflag .ne. 0) return
        
        if (j .eq. 1) fnorm1 = fnorm_min(j)
        if (fnorm_min(j) .le. fnorm1) then
           jmin = j
           fnorm1 = fnorm_min(jmin)
           pnorm  = pnorm_min(jmin)
           par    = par_min(jmin)
           delta  = delta_min(jmin)
#ifdef CRAY
           do k = 1, n
              read (j+1000) wa1(k), wa2(k)
              do istat = 1, n
                 read (j+1000) fjac(k, istat)
              end do
           end do
           do k = 1, m   
              read (j+1000) wa4(k)
           end do   
#else
           read (j+1000) wa1, wa2, wa4, fjac(1:n, 1:n)
#endif           
        end if   

        close (j+1000, status='delete')                        !!Needed to run correctly in multi-tasking...

      end do
      
      do j = 1, num_lm_params
         ext = ' '
         low_mark = ' '
         if (j .eq. 1) ext = '*'
         if (j .eq. jmin) low_mark = '*'
         write (6, '(2x,i6,4x,2(3x,1es12.4,a),3x,1es12.4)') j+ncnt,
     1         fnorm_min(j)**2, low_mark, par_min(j), ext, delta_min(j)
      end do

      call vmec_flush(6)

      lev_step_range = 0
      if( delta == maxval(delta_min)) then
         lev_step_range = 1
      else if (delta == minval(delta_min)) then
         lev_step_range = -1
      endif

!
!     Do any special cleanup now for IFLAG = FLAG_CLEANUP. WA4 LEFT UNCHANGED
!
      iflag = FLAG_CLEANUP
      call fcn(m, n, x, wa4, iflag, ncnt)

      time = time + real(ic2 - ic1)/real(irate)                !!Time in multi-process call

      end subroutine levmarq_param


      subroutine lmpar_parallel(j, fcn)
      use fdjac_mod
      use lmpar_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: j
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), dimension(11), parameter :: factors =
     1  (/ 1.0_dp, 0.5_dp, 0.25_dp, 0.128_dp, 2.1_dp, 0.75_dp,
     2      1.25_dp, 1.5_dp, 0.9_dp, 1.1_dp, 1.75_dp /)
      real(rprec) :: sjp1 = 2, sjnp = 10      
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat, iflag, k, nfact
      real(rprec) :: deltain, parin, fnorm_in, pnorm_in, scale_factor
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      external fcn
      real(rprec), external :: enorm
C-----------------------------------------------      
!
!     THIS ROUTINE IS PASSED TO THE MULTI-PROCESSOR HANDLING
!     ROUTINE


c ***************************************************
c  stepping algorithm similar to that used in the orignal parallel optimizer
c  by M.Zarnstorff and S. Ethier,  Feb. 1999
c
c  Re-implemented,  MCZ  July 2000
c ***************************************************
      nfact = size(factors)

      if (first .and. num_lm_params > 2) then            
c
c       do an exponential spread the first time to see where we are
c
         scale_factor = exp((j-1)*log(spread_ratio)/num_lm_params)
      else if (num_lm_params > 2*nfact) then
         scale_factor = (j*maxval(factors))/num_lm_params
      else if (j .le. nfact) then
         scale_factor = factors(j)
      else
         scale_factor =((j-nfact)*minval(factors)*0.7)/
     1                  (num_lm_params-nfact)
      endif

      deltain = delta * scale_factor

!
!     Compute perturbation vector (wa1 ) and Lev/Marq parameter (par)
!     for different tolerances, delta
!      
      parin = par

      call lmpar (n, fjac, m, ipvt, diag, qtf, deltain, parin,
     1            wa1, wa2, wa3, wa4)

!
!     store the direction p and x + p. calculate the norm of p.
!
      wa1 = -wa1
      if (parin .eq. 0._dp) wa1 = wa1*scale_factor
      wa2 = x + wa1
      wa3 = diag*wa1
      pnorm_in = enorm(n, wa3)

c
c     evaluate the function at x + p and calculate its norm.
c
      iflag = j
      call fcn (m, n, wa2, wa4, iflag, ncnt)

      fnorm_in = enorm(m, wa4)

!
!     OPEN A UNIQUE FILE FOR I/O IN MULTI-PROCESSOR SYSTEM
!
      write (j+1000) j, iflag, pnorm_in, fnorm_in, parin, deltain
#ifdef CRAY
      do k = 1, n
         write (j+1000) wa1(k), wa2(k)
         do istat = 1, n
            write (j+1000) fjac(k, istat)
         end do
      end do
      do k = 1, m      
         write (j+1000) wa4(k)
      end do   
#else
      write (j+1000) wa1, wa2, wa4, fjac(1:n, 1:n)
#endif            
      close (j+1000)                      !!Needed to run correctly in multi-tasking...

      end subroutine lmpar_parallel


      subroutine stepopt(fcn, x, fvec, fnorm, iflag, nfev, epsfcn, time)

      use kind_spec
c     use fdjac_mod, ONLY: flip, ix_min, jac_order
      use fdjac_mod
      use lmpar_mod, ONLY: m, n, wa2, wa4
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
c     integer, intent(in) :: m, n
      integer :: nfev
      integer, intent(out) :: iflag
      real(rprec), intent(in) :: epsfcn
      real(rprec) :: fnorm, time
      real(rprec), dimension(n) :: x
      real(rprec), dimension(m) :: fvec
c     real(rprec), target :: wa2(n), wa4(m)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, parameter :: FLAG_CLEANUP = -100
      integer :: j, k, istat, iread, ic1, ic2, irate, count_max,
     1     jmin, ll, kk, mm    
      real(rprec), dimension(num_lm_params) :: fnorm_array
      real(rprec) fnorm_min
      integer, dimension(num_lm_params) :: ll, kk, mm
      character*1 :: ext, low_mark

C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      external fcn, dpmpar, enorm, stepopt_parallel
C-----------------------------------------------
C   L o c a l   Parameters
C-----------------------------------------------
      integer, parameter :: skip = 4
      integer, parameter :: master = 0                       !mpi stuff
C-----------------------------------------------
c
c     subroutine stepopt
c
c     this subroutine uses the ordered array of jacobian-evaluated 
c     improvements to attempt ot find incremental improvements in fnorm
c     in the vicinity of the previous best case.  This routine is called
c     when the optimizer is having difficulty interpretting the local
c     gradient, as syptified by the levenberg-marquardt step not giving the
c     best improvement.  Rather than immediately going for another Jacobian 
c     evaluation, we use the information from the last one to hunt a bit.
c
c     M. Zarnstorff  March 2002
C-----------------------------------------------

C-----------------------------------------------      
!
!     Initialize module variables for use in stepopt_parallel
!
      
c     xp => x
c     wa2p => wa2
c     wa4p => wa4
c     np = n
c     mp = m

      wa2 = x
      ncnt = nfev

      call system_clock(ic1, irate)

      call multiprocess(num_lm_params, max_processors, 
     1    stepopt_parallel, fcn)

      call system_clock(ic2, irate, count_max)
      if (ic2 .lt. ic1) ic2 = ic2 + count_max

      nfev = nfev + num_lm_params

!
!     Read in optimal wa1, wa2, wa4, fnorm, and indexing  values from file
!
      do j = 1, num_lm_params

        read (j+1000, iostat=iread) istat, iflag,fnorm_array(j),
     1                              ll(j),kk(j),mm(j) 
        if (iread .ne. 0) then
           write (6, *) 'Error reading from file fort.', j+1000,
     1       ' in stepopt', ' IOSTAT = ', iread           
           iflag = -15
        else if (j .ne. istat) then
           write (6, *)
     1        'Incorrect value read in for index j in stepopt'
           iflag = -15
        end if    

        if (iflag .ne. 0) return
        
        if (j .eq. 1) fnorm_min = 2*fnorm_array(j)
        if (fnorm_array(j) < fnorm_min) then
           jmin = j
           fnorm_min = fnorm_array(jmin)
#ifdef CRAY
           do k = 1, n
              read (j+1000) wa2(k)
           end do
           do k = 1, m   
              read (j+1000) wa4(k)
           end do   
#else
           read (j+1000) wa2, wa4
#endif           
        end if   

        close (j+1000, status='delete')                        !!Needed to run correctly in multi-tasking...

      end do
      
      low_mark = ' '
      if (fnorm < fnorm_min) then
         low_mark = '*'
         jmin = num_lm_params+1
      endif

      write (6, '(2x,6x,4x,1es12.4,a,3x,a,i4,2x,i4,2x,i4,a)')
     1      fnorm**2, low_mark, '(', 
     2      1, 1, mm(1), ')'

      do j = 1, num_lm_params
         low_mark = ' '
         if (j .eq. jmin) low_mark = '*'
         write (6, '(2x,i6,4x,1es12.4,a,3x,a,i4,2x,i4,2x,i4,a)')
     1      j+ncnt, fnorm_array(j)**2, low_mark, '(', 
     2      ll(j), kk(j), mm(j), ')'
      end do

      if( fnorm_min < fnorm) then
         x = wa2
         fvec = wa4
         fnorm = fnorm_min
      end if

!
!     Do any special cleanup now for IFLAG = FLAG_CLEANUP. WA4 LEFT UNCHANGED
!
      iflag = FLAG_CLEANUP
      call fcn(m, n, x, wa4, iflag, ncnt)

      time = time + real(ic2 - ic1)/real(irate)              !!Time in multi-process call

      end subroutine stepopt


      subroutine stepopt_parallel(j, fcn)
      use fdjac_mod
      use lmpar_mod, ONLY: m, n, wa2, wa4
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: j
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: skip = 4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat, iflag, k, l, i, jj
      real(rprec) :: temp, fnorm_in, h
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      external fcn
      real(rprec), external :: enorm
C-----------------------------------------------      
!
!     THIS ROUTINE IS PASSED TO THE MULTI-PROCESSOR HANDLING
!     ROUTINE


!
!     store the direction p and x + p. calculate the norm of p.
!
!      l = ((j-1) / jac_count)*skip + 1         ! original version

      l = ((j-1) / jac_count)
      if( l <= skip/2 ) then
         l = l+1
      else
         l = (l-skip/2) * skip + 1
      endif

      k = mod((j-1), jac_count) + 1
      if( k == 1) l = l+1

c     wa2  = x

      do i=1, k
         jj = jac_order(i)
         temp = wa2(jj)
         h = eps*abs(temp)
         if (h .eq. 0) h = eps
         if( flip(jj)) h = -h
         if( i .ne. 1) then
            wa2(jj) = temp + l*h
         else
            wa2(jj) = temp + (l-1)*h
         endif
      enddo
      
c
c     evaluate the function at x + p and calculate its norm.
c
      iflag = j
      call fcn (m, n, wa2, wa4, iflag, ncnt)

      fnorm_in = enorm(m, wa4)

!
!     OPEN A UNIQUE FILE FOR I/O IN MULTI-PROCESSOR SYSTEM
!
      write (j+1000) j, iflag, fnorm_in, l, k, jac_order(k)

#ifdef CRAY
      do k = 1, n
         write (j+1000) wa2(k)
      end do
      do k = 1, m      
         write (j+1000) wa4(k)
      end do   
#else
      write (j+1000) wa2, wa4
#endif            
      close (j+1000)                      !!Needed to run correctly in multi-tasking...

      end subroutine stepopt_parallel

#endif      

      subroutine lmpar (n, r, ldr, ipvt, diag, qtb, delta, par, x,
     1                 sdiag, wa1, wa2)
      use kind_spec 
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: n, ldr, ipvt(n)
      real(rprec) :: delta, par
      real(rprec), dimension(ldr,n) :: r
      real(rprec), dimension(n) :: wa1, wa2
      real(rprec), dimension(n), intent(in) :: diag, qtb
      real(rprec), dimension(n), intent(out) :: x, sdiag
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: p1=0.1_dp, zero = 0,
     1     p001 = 0.001_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: iter, j, jm1, jp1, k, l, nsing
      real(rprec) :: dxnorm, dwarf, fp, gnorm, parc, parl, paru,
     1   sum0, enorm, dpmpar, temp, epsmch
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      external dpmpar, enorm
C-----------------------------------------------
c
c     subroutine lmpar
c
c     given an m by n matrix a, an n by n nonsingular diagonal
c     matrix d, an m-vector b, and a positive number delta,
c     the problem is to determine a value for the parameter
c     par such that if x solves the system
c
c           a*x = b ,     sqrt(par)*d*x = 0 ,
c
c     in the least squares sense, and dxnorm is the euclidean
c     norm of d*x, then either par is zero and
c
c           (dxnorm-delta) .le. 0.1*delta ,
c
c     or par is positive and
c
c           abs(dxnorm-delta) .le. 0.1*delta .
c
c     this subroutine completes the solution of the problem
c     if it is provided with the necessary information from the
c     qr factorization, with column pivoting, of a. that is, if
c     a*p = q*r, where p is a permutation matrix, q has orthogonal
c     columns, and r is an upper triangular matrix with diagonal
c     elements of nonincreasing magnitude, then lmpar expects
c     the full upper triangle of r, the permutation matrix p,
c     and the first n components of (q transpose)*b. on output
c     lmpar also provides an upper triangular matrix s such that
c
c            t   t                   t
c           p *(a *a + par*d*d)*p = s *s .
c
c     s is employed within lmpar and may be of separate interest.
c
c     only a few iterations are generally needed for convergence
c     of the algorithm. if, however, the limit of 10 iterations
c     is reached, then the output par will contain the best
c     value obtained so far.
c
c     the subroutine statement is
c
c       subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,
c                        wa1,wa2)
c
c     where
c
c       n is a positive integer input variable set to the order of r.
c
c       r is an n by n array. on input the full upper triangle
c         must contain the full upper triangle of the matrix r.
c         on output the full upper triangle is unaltered, and the
c         strict lower triangle contains the strict upper triangle
c         (transposed) of the upper triangular matrix s.
c
c       ldr is a positive integer input variable not less than n
c         which specifies the leading dimension of the array r.
c
c       ipvt is an integer input array of length n which defines the
c         permutation matrix p such that a*p = q*r. column j of p
c         is column ipvt(j) of the identity matrix.
c
c       diag is an input array of length n which must contain the
c         diagonal elements of the matrix d.
c
c       qtb is an input array of length n which must contain the first
c         n elements of the vector (q transpose)*b.
c
c       delta is a positive input variable which specifies an upper
c         bound on the euclidean norm of d*x.
c
c       par is a nonnegative variable. on input par contains an
c         initial estimate of the levenberg-marquardt parameter.
c         on output par contains the final estimate.
c
c       x is an output array of length n which contains the least
c         squares solution of the system a*x = b, sqrt(par)*d*x = 0,
c         for the output par.
c
c       sdiag is an output array of length n which contains the
c         diagonal elements of the upper triangular matrix s.
c
c       wa1 and wa2 are work arrays of length n.
c
c     subprograms called
c
c       minpack-supplied ... dpmpar,enorm,qrsolv
c
c       fortran-supplied ... abs,max,min,sqrt
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
 
c
c     dwarf is the smallest positive magnitude.
c
      dwarf = dpmpar(2)
      epsmch = dpmpar(1)
      epsmch = 100*epsmch                  !!SPH: gives more reliable results
c
c     compute and store in x the gauss-newton direction. if the
c     jacobian is rank-deficient, obtain a least squares solution.
c
      nsing = n
      do j = 1, n
         wa1(j) = qtb(j)
!SPH     if (r(j,j).eq.zero .and. nsing.eq.n) nsing = j - 1
         if (abs(r(j,j)).le.epsmch .and. nsing.eq.n) nsing = j - 1
         if (nsing .lt. n) wa1(j) = zero
      end do
      if (nsing .ge. 1) then
         do k = 1, nsing
            j = nsing - k + 1
            wa1(j) = wa1(j)/r(j,j)
            temp = wa1(j)
            jm1 = j - 1
            if (jm1 .ge. 1) then
               wa1(:jm1) = wa1(:jm1) - r(:jm1,j)*temp
            endif
         end do
      endif
      do j = 1, n
         l = ipvt(j)
         x(l) = wa1(j)
      end do
c
c     initialize the iteration counter.
c     evaluate the function at the origin, and test
c     for acceptance of the gauss-newton direction.
c
      wa2 = diag * x
      dxnorm = enorm(n,wa2)
      fp = dxnorm - delta

      if (fp .le. p1*delta) then
         par = zero
         return
      end if

c
c     BEGIN GAUSS-NEWTON STEP
c
c     if the jacobian is not rank deficient, the newton
c     step provides a lower bound, parl, for the zero of
c     the function. otherwise set this bound to zero.
c
      parl = zero
      if (nsing .ge. n) then
         wa1 = diag(ipvt)*(wa2(ipvt)/dxnorm)
         do j = 1, n
            sum0 = zero
            jm1 = j - 1
            if (jm1 .ge. 1) sum0 = sum(r(:jm1,j)*wa1(:jm1))
            wa1(j) = (wa1(j)-sum0)/r(j,j)
         end do
         temp = enorm(n,wa1)
         parl = ((fp/delta)/temp)/temp
      endif
c
c     calculate an upper bound, paru, for the zero of the function.
c
      do j = 1, n
         sum0 = sum(r(:j,j)*qtb(:j))
         l = ipvt(j)
         wa1(j) = sum0/diag(l)
      end do
      gnorm = enorm(n,wa1)
      paru = gnorm/delta
      if (paru .eq. zero) paru = dwarf/min(delta,p1)
c
c     if the input par lies outside of the interval (parl,paru),
c     set par to the closer endpoint.
c
      par = max(par,parl)
      par = min(par,paru)
      if (par .eq. zero) par = gnorm/dxnorm

c
c     beginning of an iteration loop.
c
      do iter = 1, 10
c
c        evaluate the function at the current value of par.
c
         if (par .le. zero) par = max(dwarf,p001*paru)
         temp = sqrt(par)

         wa1 = temp*diag
         call qrsolv (n, r, ldr, ipvt, wa1, qtb, x, sdiag, wa2)
         wa2 = diag*x
         dxnorm = enorm(n,wa2)
         temp = fp
         fp = dxnorm - delta
c
c        if the function is small enough, accept the current value
c        of par. also test for the exceptional cases where parl
c        is zero or the number of iterations has reached 10.
c
         if (abs(fp).le.p1*delta .or. parl.eq.zero .and. 
     1       fp.le.temp .and. temp.lt.zero .or. iter.eq.10) exit
c
c        compute the newton correction.
c
         wa1 = diag(ipvt)*(wa2(ipvt)/dxnorm)
         do j = 1, n
            wa1(j) = wa1(j)/sdiag(j)
            temp = wa1(j)
            jp1 = j + 1
            if (n .ge. jp1) wa1(jp1:n) = wa1(jp1:n) - r(jp1:n,j)*temp
         end do
         temp = enorm(n,wa1)
         parc = ((fp/delta)/temp)/temp
c
c        depending on the sign of the function, update parl or paru.
c
         if (fp .gt. zero) parl = max(parl,par)
         if (fp .lt. zero) paru = min(paru,par)
c
c        compute an improved estimate for par.
c
         par = max(parl,par + parc)

      enddo                        !!end of an iteration.
 
      end subroutine lmpar
 

      subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: m, n, lda, lipvt
      logical, intent(in) :: pivot
      integer, dimension(lipvt), intent(out) :: ipvt
      real(rprec), dimension(lda,n), intent(inout) :: a
      real(rprec), dimension(n), intent(out) :: rdiag, acnorm, wa
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: one = 1, p05 = 0.05_dp,
     1    zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: j, jp1, k, kmax, minmn
      real(rprec) :: ajnorm, epsmch, sum0, temp, enorm, dpmpar
      real(rprec) :: temp1u(m)
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      external dpmpar,enorm
C-----------------------------------------------
c
c     subroutine qrfac
c
c     this subroutine uses householder transformations with column
c     pivoting (optional) to compute a qr factorization of the
c     m by n matrix a. that is, qrfac determines an orthogonal
c     matrix q, a permutation matrix p, and an upper trapezoidal
c     matrix r with diagonal elements of nonincreasing magnitude,
c     such that a*p = q*r. the householder transformation for
c     column k, k = 1,2,...,min(m,n), is of the form
c
c                           t
c           i - (1/u(k))*u*u
c
c     where u has zeros in the first k-1 positions. the form of
c     this transformation and the method of pivoting first
c     appeared in the corresponding linpack subroutine.
c
c     the subroutine statement is
c
c       subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
c
c     where
c
c       m is a positive integer input variable set to the number
c         of rows of a.
c
c       n is a positive integer input variable set to the number
c         of columns of a.
c
c       a is an m by n array. on input a contains the matrix for
c         which the qr factorization is to be computed. on output
c         the strict upper trapezoidal part of a contains the strict
c         upper trapezoidal part of r, and the lower trapezoidal
c         part of a contains a factored form of q (the non-trivial
c         elements of the u vectors described above).
c
c       lda is a positive integer input variable not less than m
c         which specifies the leading dimension of the array a.
c
c       pivot is a logical input variable. if pivot is set true,
c         then column pivoting is enforced. if pivot is set false,
c         then no column pivoting is done.
c
c       ipvt is an integer output array of length lipvt. ipvt
c         defines the permutation matrix p such that a*p = q*r.
c         column j of p is column ipvt(j) of the identity matrix.
c         if pivot is false, ipvt is not referenced.
c
c       lipvt is a positive integer input variable. if pivot is false,
c         then lipvt may be as small as 1. if pivot is true, then
c         lipvt must be at least n.
c
c       rdiag is an output array of length n which contains the
c         diagonal elements of r.
c
c       acnorm is an output array of length n which contains the
c         norms of the corresponding columns of the input matrix a.
c         if this information is not needed, then acnorm can coincide
c         with rdiag.
c
c       wa is a work array of length n. if pivot is false, then wa
c         can coincide with rdiag.
c
c     subprograms called
c
c       minpack-supplied ... dpmpar,enorm
c
c       fortran-supplied ... max,sqrt,min
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
c
c     epsmch is the machine precision.
c
      epsmch = dpmpar(1)
c
c     compute the initial column norms and initialize several arrays.
c
      do j = 1, n
         acnorm(j) = enorm(m,a(1,j))
         rdiag(j) = acnorm(j)
         wa(j) = rdiag(j)
         if (pivot) ipvt(j) = j
      end do
c
c     reduce a to r with householder transformations.
c
      minmn = min(m,n)

      do j = 1, minmn
         if (pivot) then
c
c        bring the column of largest norm into the pivot position.
c
            kmax = j
            do k = 1, n - j + 1
               if (rdiag(j+k-1) .gt. rdiag(kmax)) kmax = j + k - 1
            end do
            if (kmax .ne. j) then
               temp1u(:m) = a(:m,j)
               a(:m,j) = a(:m,kmax)
               a(:m,kmax) = temp1u(:m)
               rdiag(kmax) = rdiag(j)
               wa(kmax) = wa(j)
               k = ipvt(j)
               ipvt(j) = ipvt(kmax)
               ipvt(kmax) = k
            endif
         endif
c
c        compute the householder transformation to reduce the
c        j-th column of a to a multiple of the j-th unit vector.
c
         ajnorm = enorm(m - j + 1,a(j,j))
         if (ajnorm .ne. zero) then
            if (a(j,j) .lt. zero) ajnorm = -ajnorm
            a(j:m,j) = a(j:m,j)/ajnorm
            a(j,j) = a(j,j) + one
c
c        apply the transformation to the remaining columns
c        and update the norms.
c
            jp1 = j + 1
            if (n .ge. jp1) then
               do k = jp1, n
                  sum0 = sum(a(j:m,j)*a(j:m,k))
                  temp = sum0/a(j,j)
                  a(j:m,k) = a(j:m,k) - temp*a(j:m,j)
                  if (pivot .and. rdiag(k).ne.zero) then
                     temp = a(j,k)/rdiag(k)
                     rdiag(k) = rdiag(k)*sqrt(max(zero,one - temp*temp))
                     if (p05*(rdiag(k)/wa(k))**2 .le. epsmch) then
                        rdiag(k) = enorm(m - j,a(jp1,k))
                        wa(k) = rdiag(k)
                     endif
                  endif
               end do
            endif
         endif
         rdiag(j) = -ajnorm
      end do

      end subroutine qrfac
 

      subroutine qrsolv(n, r, ldr, ipvt, diag, qtb, x, sdiag, wa)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: n, ldr
      integer, dimension(n), intent(in) :: ipvt
      real(rprec), dimension(ldr,n) :: r
      real(rprec), dimension(n) :: diag, qtb, x, sdiag, wa
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: one = 1.0_dp, zero = 0.0_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, j, jp1, k, kp1, l, nsing, l1
      real(rprec) :: cos, cotan, qtbpj, sin, sum0, tan, temp
      real(rprec) , allocatable :: temp1u(:)
      logical l1s(n)
C-----------------------------------------------
c
c     subroutine qrsolv
c
c     given an m by n matrix a, an n by n diagonal matrix d,
c     and an m-vector b, the problem is to determine an x which
c     solves the system
c
c           a*x = b ,     d*x = 0 ,
c
c     in the least squares sense.
c
c     this subroutine completes the solution of the problem
c     if it is provided with the necessary information from the
c     qr factorization, with column pivoting, of a. that is, if
c     a*p = q*r, where p is a permutation matrix, q has orthogonal
c     columns, and r is an upper triangular matrix with diagonal
c     elements of nonincreasing magnitude, then qrsolv expects
c     the full upper triangle of r, the permutation matrix p,
c     and the first n components of (q transpose)*b. the system
c     a*x = b, d*x = 0, is then equivalent to
c
c                  t       t
c           r*z = q *b ,  p *d*p*z = 0 ,
c
c     where x = p*z. if this system does not have full rank,
c     then a least squares solution is obtained. on output qrsolv
c     also provides an upper triangular matrix s such that
c
c            t   t               t
c           p *(a *a + d*d)*p = s *s .
c
c     s is computed within qrsolv and may be of separate interest.
c
c     the subroutine statement is
c
c       subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)
c
c     where
c
c       n is a positive integer input variable set to the order of r.
c
c       r is an n by n array. on input the full upper triangle
c         must contain the full upper triangle of the matrix r.
c         on output the full upper triangle is unaltered, and the
c         strict lower triangle contains the strict upper triangle
c         (transposed) of the upper triangular matrix s.
c
c       ldr is a positive integer input variable not less than n
c         which specifies the leading dimension of the array r.
c
c       ipvt is an integer input array of length n which defines the
c         permutation matrix p such that a*p = q*r. column j of p
c         is column ipvt(j) of the identity matrix.
c
c       diag is an input array of length n which must contain the
c         diagonal elements of the matrix d.
c
c       qtb is an input array of length n which must contain the first
c         n elements of the vector (q transpose)*b.
c
c       x is an output array of length n which contains the least
c         squares solution of the system a*x = b, d*x = 0.
c
c       sdiag is an output array of length n which contains the
c         diagonal elements of the upper triangular matrix s.
c
c       wa is a work array of length n.
c
c     subprograms called
c
c       fortran-supplied ... abs,sqrt
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
c
c     copy r and (q transpose)*b to preserve input and initialize s.
c     in particular, save the diagonal elements of r in x.
c
      do j = 1, n
         r(j:n,j) = r(j,j:n)
         x(j) = r(j,j)
      end do
      wa(1:n) = qtb(1:n)
c
c     eliminate the diagonal matrix d using a givens rotation.
c
      do j = 1, n
c
c        prepare the row of d to be eliminated, locating the
c        diagonal element using p from the qr factorization.
c
         l = ipvt(j)
         if (diag(l) .ne. zero) then
            sdiag(j:n) = zero
            sdiag(j) = diag(l)
c
c        the transformations to eliminate the row of d
c        modify only a single element of (q transpose)*b
c        beyond the first n, which is initially zero.
c
            qtbpj = zero
            do k = j, n
c
c           determine a givens rotation which eliminates the
c           appropriate element in the current row of d.
c
               if (sdiag(k) .ne. zero) then
                  if (abs(r(k,k)) .lt. abs(sdiag(k))) then
                     cotan = r(k,k)/sdiag(k)
                     sin = one/sqrt(one + cotan*cotan)
                     cos = sin*cotan
                  else
                     tan = sdiag(k)/r(k,k)
                     cos = one/sqrt(one + tan*tan)
                     sin = cos*tan
                  endif
c
c           compute the modified diagonal element of r and
c           the modified element of ((q transpose)*b,0).
c
                  r(k,k) = cos*r(k,k) + sin*sdiag(k)
                  temp = cos*wa(k) + sin*qtbpj
                  qtbpj = (-sin*wa(k)) + cos*qtbpj
                  wa(k) = temp
c
c           accumulate the tranformation in the row of s.
c
                  kp1 = k + 1
                  if (n .ge. kp1) then
                     l1 = n-kp1+1
                     allocate (temp1u(l1))
                     temp1u(:l1) = cos*r(kp1:n,k) + sin*sdiag(kp1:n)
                     sdiag(kp1:n) = (-sin*r(kp1:n,k)) + cos*sdiag(kp1:n)
                     r(kp1:n,k) = temp1u(:l1)
                     deallocate (temp1u)
                  endif
               endif
            end do
         endif
c
c        store the diagonal element of s and restore
c        the corresponding diagonal element of r.
c
         sdiag(j) = r(j,j)
         r(j,j) = x(j)
      end do
c
c     solve the triangular system for z. if the system is
c     singular, then obtain a least squares solution.
c
      nsing = n
      do j = 1, n
         if (sdiag(j).eq.zero .and. nsing.eq.n) nsing = j - 1
         l1s = nsing .lt. n
      end do
      where (l1s) wa = zero
      if (nsing .ge. 1) then
         do k = 1, nsing
            j = nsing - k + 1
            sum0 = zero
            jp1 = j + 1
            if (nsing .ge. jp1) then
               sum0 = sum(r(jp1:nsing,j)*wa(jp1:nsing))
            endif
            wa(j) = (wa(j)-sum0)/sdiag(j)
         end do
      endif
c
c     permute the components of z back to components of x.
c
      do j = 1, n
         l = ipvt(j)
         x(l) = wa(j)
      end do

      end subroutine qrsolv

       
      function dpmpar (i)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: i
      real(rprec) :: dpmpar
C-----------------------------------------------
c
c     function dpmpar
c
c     this function provides single (double) precision machine parameters
c     when the appropriate set of data statements is activated (by
c     removing the c from column 1) and all other data statements are
c     rendered inactive. most of the parameter values were obtained
c     from the corresponding bell laboratories port library function.
c
c     the function statement is
c
c       function dpmpar(i)
c
c     where
c
c       i is an integer input variable set to 1, 2, or 3 which
c         selects the desired machine parameter. if the machine has
c         t base b digits and its smallest and largest exponents are
c         emin and emax, respectively, then these parameters are
c
c         dpmpar(1) = b**(1 - t), the machine precision,
c
c         dpmpar(2) = b**(emin - 1), the smallest magnitude,
c
c         dpmpar(3) = b**emax*(1 - b**(-t)), the largest magnitude.
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     SINGLE PRECISION
c     data rmach(1)/1.490116100E-8/
c     data rmach(2)/1.4693679000E-30/
c     data rmach(3)/1.701411800E+30/
c     DOUBLE PRECISION - IEEE
c     data dmach(1) /2.22044604926d-16/
c     data dmach(2) /2.22507385852d-308/
c     data dmach(3) /1.79769313485d+308/
c     modified for f90 (sph, august 1997)
c
      select case(i)
      case(:1)
        dpmpar = epsilon(dpmpar)      !2.22044604926e-16_dp
      case(2)
        dpmpar = tiny(dpmpar)         !2.22507385852e-308_dp
      case(3:)
        dpmpar = huge(dpmpar)         !1.79769313485e+308_dp
      end select
 
      end function dpmpar
