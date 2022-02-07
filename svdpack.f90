#!/bin/sh
#---------------------------------------------------------------------
cat > temp.c << "EOC"
cat > svdpack.f << "EOF"     
      subroutine svd_solve(m, n, mp, np, a, b, v, w, nw)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer m, n, mp, np, nw
      real(rprec), dimension(mp,np) :: a
      real(rprec), dimension(n,n) :: v
      real(rprec), dimension(n) :: w
      real(rprec), dimension(m) :: b
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec) :: small = 1.0e-10_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat, i, j
      real(rprec), allocatable :: u(:,:)
C-----------------------------------------------
 
c       Solves Matrix equation A * V = B for V using SVD method
c       Uses svd routines from numerical recipes
c       Prashant Valanju (Dec 1998) pvalanju@mail.utexas.edu
 
c       Almost same as SvdSolveB except this one finds
c       nw solutions V(:,i) by keeping i = 1 to nw vectors, where
c       nw is the number of vectors with nonzero weights > small

c       It also uses V internally for SVD decomposition to save space, so
 
c       Inputs:
c       A(M,N) - Matrix with physical dimensions (Mp,Np)
c       B(M)   - R.H.S. of A X = B
c       M      - Number of rows of A to use
c       N      - Number of columns of A to use
c       Mp,NP  - Storage dimensions of A(Mp,Np)
c
c       Output:
c       nw      - number of vectors with normalized weights above small
c          Eventually, this will be the number of optimum weights
c          after we decide on a criterion for optimization
c       V(N,nw) - nw Solutions of A V = B, each a vector of length N
c          V(:N,iw) = solution with top iw weights kept
c          Physical dimensions of V are (N,N)
c       w(n)    - Weights in decreasing order
c
 
c  Allocate local arrays. 
c  Note u(m,n) is enough because needed part of a(mp,np) is copied to
c  u(m,n), and a(mp,np) is never used directly. This saves space.
c  It is essential to use a local u since svdcmp changes it

      allocate (u(m,n),  stat=istat)
      if(istat.ne.0) stop 'Stop: No memory in svd_nesc'  
 
c.......................................
c  Initialize all to zero to wipe out effects of old call
      w(:n) = 0                                !Zero all weights
      do j = 1, n
         u(:m,j) = a(:m,j)                  !Because U will be changed by svdcmp
         v(:n,j) = 0
      end do

c  Do the SVD decomposition of a, i.e, of u into u, v and w 
       call svdcmp (u, m, n, m, n, w, v)
 
c  Sort weights and matrices with DECREASING weights so w(1) is biggest
c  Permute weight w(i) AND column vectors U(*,i), V(*,i) at the same time
       call sortsvd (m, n, m, n, w, u, v)
 
c  Find nw = number of large weights (dcreasing ordered by sortsvd)
         do nw = n, 1, -1          !Find first large weight and get out
            if ( abs(w(nw)/w(1)) .gt. small) exit
         end do

c  Find nw solutions by successively adding up eigenvectors V 
c     with correct weight given by (U(i) dot b)/w(i) in eq 14.3.17 (NR). 
c  The coeff of ith vector V(i) is the dot product of the i-th column
c     of U with the rhs vector B, divided by w(i)
c  This does the svdbksb directly, faster than the NR svdbksub routine
c     and uses less memory due to the dual role of V
c  Note: any optimization scheme to find the 'best' nw will 

c      First set the 1st vector with largest weight w(1)
       v(:n,1) = sum(u(:m,1)*b(:m)) *v(:n,1) /w(1)
c      Next add the vectors with successive weights (in decreasing order)
         do i = 2, nw
            j = i - 1
            v(:n,i) = v(:n,j) + sum(u(:m,i)*b(:m)) *v(:n,i) / w(i)
         end do
c................................................

      deallocate (u, stat=istat)

!     end subroutine svd_solve
      end


      subroutine svdsolveb(m, n, mp, np, wcut, file, a, b, x)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer m, n, mp, np
      real(rprec) wcut
      real(rprec), dimension(mp,np) :: a
      real(rprec), dimension(n) :: b, x
      character file*(*)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec) :: zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: n1, nonzero, nlast, istat, i, j, iunit
      real(rprec), allocatable :: u(:,:), v(:,:), w(:)
      real(rprec) :: wmax, wmin, wcuta, udotb, err
C-----------------------------------------------
 
c       Solves Matrix equation A X = B for X using SVD method
c       Allows supression of weights, storage and recall of answers
c       Uses svd routines from numerical recipes
c       Prashant Valanju (Sept 1998) pvalanju@mail.utexas.edu
 
c       Almost same as SvdSolveUV except this one stores basis vectors
c       rather than U, w, V. This one is a bit faster than SvdSolveUV, so
c       use this one when you want to change only weights, not B or A
 
c       Inputs:
c       A(M,N) - Matrix with physical dimensions (Mp,Np)
c       B(N)   - R.H.S. of A X = B
c       M      - Number of rows of A to use
c       N      - Number of columns of A to use
c       Mp,NP  - Storage dimensions of A(Mp,Np)
c       Wcut  - Cutoff parameter
c                If Wcut = 0, writes answers into "file", all weights are kept
c                If Wcut < 0, reads previous answers from "file",
c                              and uses -Wcut for weights to keep
c                If 0 < Wcut < 1 : Cuts off weights with w(i)/wmax < wcut
c                      Use this to cut off small weights, does calculation
c                If Wcut = integer > 1, keeps Wcut weights, does calculation
c                      Use this to cut off fixed number of weights
c       Output:
c       X(N)   - Solution of A X = B
c
      n1 = n
 
c     allocate local arrays
      allocate (u(mp,np), w(np), v(np,np), stat=istat)

      if (istat .eq. 0) then                       !have enough memory
 
         iunit = 41

         x(:n1) = zero                    !zero answer from previous call
         wcuta = abs(wcut)
c...............................
c     If Wcut < 0 use precalculated weights to get answer fast
c        use wcuta = abs(wcut) for everything
c     If any error happens while reading "file", do full svd calculation
         if (wcut .lt. 0) then         !Read file from previous calculation
            call safe_open(iunit, istat, file, 'old', 'formatted')
            if (istat .ne. 0) goto 98
            read (iunit, *, err=98)
            read (iunit, *, err=98)
            read (iunit, *, err=98) m, n, nonzero
c                           !Read exactly the same way they were written
            do i = 1, nonzero
               read (iunit, *, err=98) w(i)         !read i-th weight
               do j = 1, n                 !read i-th basis vector/ w(i)
c                             !read j-th entry (row) of i-th column of V
                  read (iunit, *, err=98) v(j,i)
               end do
            end do
            close(unit=iunit)
 
c        Decide how many weights to keep
            if (wcuta .gt. 1) nlast = wcuta         !Keep Nlast weights
            if (wcuta .lt. 1) then                  !cutoff small weights
               wmax = w(1)                          !weights are already ordered
               nlast = 0
               do i = 1, nonzero
                  if (w(i) .gt. wmax*wcuta) then
                     nlast = nlast + 1              !accept this weight
                  else
                     go to 96
                  endif
               end do
            endif
 
c        Just find linear sum of Nlast vectors (weights were already included)
   96       continue
            if (nlast.le.0 .or. nlast.gt.nonzero) nlast = nonzero
            do i = 1, nlast             !add ith column of V(*,i) to X(*)
               x(:n) = x(:n) + v(:n,i)  !add to the j-th entry of X
            end do
            go to 999                   !all done
         endif                          !End of precalculated branch Wcut < 0
   98    continue
c.......................................
c     If Wcut .ge. 0, do full svd calculation
c     Initialize all to zero to wipe out effects of old call
         do j = 1, n
            w(j) = zero                           !Zero all weights
            u(:m,j) = a(:m,j)  !Because U will be changed by svdcmp
            v(:m,j) = zero
         end do
 
         call svdcmp (u, m, n, mp, np, w, v)     !Do SVD decomposition
 
c       Sort weights and matrices with DECREASING weights
c       Permute weight w(i) AND column vectors U(*,i), V(*,i) at the same time
         call sortsvd (m, n, mp, np, w, u, v)
 
c       Find the number of nonzero weights (already dcreasing ordered)
         do nonzero = n, 1, -1
c                                   !Found first nonzero weight, get out
            if (w(nonzero) .ne. 0) exit 
         end do
         if (nonzero .gt. 0) then
            if (wcuta .gt. 1) then
 
c        Decide how many weights to keep
               nlast = wcuta                     !Keep Nlast weights
            else                                 !cutoff small weights
               wmax = w(1)                  !weights are already ordered
               nlast = 0
               do i = 1, nonzero
                  if (w(i) .gt. wmax*wcuta) then
                     nlast = nlast + 1           !accept this weight
                  else
                     go to 97
                  endif
               end do
            endif
 
c       Find solution X and basis vectors (put into columns of V)
   97       continue
            if (nlast.le.0 .or. nlast.gt.nonzero) nlast = nonzero
            do i = 1, nlast  !Calc i-th coeff (U.b/w) for noNonzero w(i)
c        First calculate i-th coeff (U(i)*b/w(i)) in the sum in eq 14.3.17 (NR)
c                                Dot product of i-th column of U with B
c                                summed over all plasma points
               udotb = sum(u(:m,i)*b(:m))
c                        !=(Ui.b/wi) in eq 14.3.17, saves many divisions
               udotb = udotb/w(i)
c         Now run down the i-th column of the vector V(j,i)
c                             j-th entry (row) of vectors in eq 14.3.17
c                             Now V is the basis vector times weight
               v(:n,i) = udotb*v(:n,i)
c                             add it to the j-th entry of X
               x(:n) = x(:n) + v(:n,i)
                                  !This did the svdbksb, faster than NR
            end do                !and also calculated the basis vectors
 
c       Write all weights and weighted basis functions into SVD file (Wcut=0)
            if (wcut .eq. 0) then
               call safe_open(iunit, istat, file, 'unknown','formatted')
               write (iunit, *) 'Max w = ', w(1), ', Min w = ', 
     1                          w(nonzero)
               write (iunit, *) 'Ratio = ', w(1)/w(nonzero)
               write (iunit, *) m, n, nonzero
               do i = 1, nonzero             !write in DECREASING order of weights
                  write (iunit, *) w(i)      !write i-th weight
                  do j = 1, n                !write i-th basis vector/ w(i)
                                             !write j-th entry (row) of i-th column of V
                     write (iunit, *) v(j,i)
                  end do
               end do
               close(unit=iunit)
            endif
c................................................
         endif
      endif
  999 continue

      deallocate (u, w, v, stat=istat)

!     end subroutine svdsolveb
      end


      subroutine svdsolveuv(m, n, mp, np, wcut, file, a, b, x)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer m, n, mp, np
      real(rprec) wcut
      character file*(*)
      real(rprec), dimension(mp,np) :: a
      real(rprec), dimension(n) :: b, x
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec) :: zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: n1, nonzero, nlast, istat, i, j, iunit
      real(rprec), allocatable :: u(:,:), w(:), v(:,:)
      real(rprec) :: wmax, wmin, wcuta, udotb, err
C-----------------------------------------------
 
c       Solves Matrix equation A X = B for X using SVD method
c       Allows supression of weights, storage and recall of answers
c       Uses svd routines from numerical recipes
c       Prashant Valanju (Sept 1998) pvalanju@mail.utexas.edu
 
c       Almost same as SvdSolveB except this one stores U, w, V
c       rather than basis vectots.
c       Use this one when you want to change both B and weights for fixed A
 
c       Inputs:
c       A(M,N) - Matrix with physical dimensions (Mp,Np)
c       B(N)   - R.H.S. of A X = B
c       M      - Number of rows of A to use
c       N      - Number of columns of A to use
c       Mp,NP  - Storage dimensions of A(Mp,Np)
c       Wcut   - Cutoff parameter
c                If Wcut = 0, writes answers into "file", all weights are kept
c                If Wcut < 0, reads previous answers from "file",
c                              and uses -Wcut for weights to keep
c                If 0 < Wcut < 1 : Cuts off weights with w(i)/wmax < wcut
c                      Use this to cut off small weights, does calculation
c                If Wcut = integer > 1, keeps Wcut weights, does calculation
c                      Use this to cut off fixed number of weights
c       Output:
c       X(N)   - Solution of A X = B
c
      n1 = n
 
c     allocate local arrays
      allocate (u(mp,np), w(np), v(np,np), stat=istat)
      if (istat .eq. 0) then                       !have enough memory
         iunit = 41
 
         x(:n1) = zero                    !zero answer from previous call
         wcuta = abs(wcut)
c...............................
c     If Wcut < 0 use precalculated weights to get answer fast
c        use wcuta = abs(wcut) for everything
c     If any error happens while reading "file", do full svd calculation
         if (wcut .lt. 0) then         !Read file from previous calculation
            call safe_open(iunit, istat, file, 'old', 'formatted')
            if (istat .ne. 0) goto 98
            read (iunit, *, err=98)
            read (iunit, *, err=98)
            read (iunit, *, err=98) m, n, nonzero
            do i = 1, n     !Read exactly the same way they were written
               read (iunit, *, err=98) w(i)         !read i-th weight
               do j = 1, n
                  read (iunit, *, err=98) v(j,i)
               end do
               do j = 1, m
                  read (iunit, *, err=98) u(j,i)
               end do
            end do
            close(unit=iunit)
 
            go to 99                !Use U, w, V matrices to calculate X
         endif                     !End of precalculated branch Wcut < 0
   98    continue
c.......................................
c     If Wcut .ge. 0, do full svd calculation
c     Initialize all to zero to wipe out effects of old call
         do j = 1, n
            w(j) = zero                           !Zero all weights
            u(:m,j) = a(:m,j)  !Because U will be changed by svdcmp
            v(:m,j) = zero
         end do
 
         call svdcmp (u, m, n, mp, np, w, v)     !Do SVD decomposition
 
c       Sort weights and matrices with DECREASING weights
c       Permute weight w(i) AND column vectors U(*,i), V(*,i) at the same time
         call sortsvd (m, n, mp, np, w, u, v)
 
c       Find the number of nonzero weights (already dcreasing ordered)
         do nonzero = n, 1, -1
c                                   !Found first nonzero weight, get out
            if (w(nonzero) .ne. 0) exit 
         end do
         if (nonzero .le. 0) go to 999
 
c............................................
c        Decide how many weights to keep
   99    continue
         if (wcuta .gt. 1) then
            nlast = wcuta
         else                                    !cutoff small weights
            wmax = w(1)                     !weights are already ordered
            nlast = 0
            do i = 1, nonzero
               if (w(i) .gt. wmax*wcuta) then
                  nlast = nlast + 1              !accept this weight
               else
                  go to 96
               endif
            end do
         endif
 
c       Find solution X and basis vectors (put into columns of V)
   96    continue
         if (nlast.le.0 .or. nlast.gt.nonzero) nlast = nonzero
         do i = 1, nlast     !Calc i-th coeff (U.b/w) for noNonzero w(i)
c        First calculate i-th coeff (U(i)*b/w(i)) in the sum in eq 14.3.17 (NR)
c                            !Dot product of i-th column of U with B
c                            !summed over all plasma points
            udotb = sum(u(:m,i)*b(:m))
c                        !=(Ui.b/wi) in eq 14.3.17, saves many divisions
            udotb = udotb/w(i)
c         Now run down the i-th column of the vector V(j,i)
c                         !j-th entry (row) of vectors in eq 14.3.17
c                         !add it to the j-th entry of X
            x(:n) = x(:n) + udotb*v(:n,i)
c                                  !This did the svdbksb, faster than NR
         end do                   !and also calculated the basis vectors
 
c............................................
c       Write all weights and U, V into SVD file (if Wcut=0)
         if (wcut .eq. 0) then
            call safe_open(iunit, istat, file, 'unknown', 'formatted')
            write (iunit, *) 'Max w = ', w(1), ', Min w = ', w(nonzero)
            write (iunit, *) 'Ratio = ', w(1)/w(nonzero)
            write (iunit, *) m, n, nonzero
            do i = 1, n
               write (iunit, *) w(i)                !write i-th weight
               do j = 1, n
                  write (iunit, *) v(j,i)
               end do
               do j = 1, m
                  write (iunit, *) u(j,i)
               end do
            end do
            close(unit=iunit)
         endif
c................................................
      endif
  999 continue

      deallocate (u, w, v, stat=istat)

!     end subroutine svdsolveuv
      end


      subroutine svdinv(m, a, wcut, file)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer m
      real(rprec) wcut
      character file*(*)
      real(rprec), dimension(m,m) :: a
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec) :: zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nonzero, nlast, istat, i, j, k, iunit
      real(rprec), allocatable :: u(:,:), w(:), v(:,:)
      real(rprec) :: wmax, wmin, wcuta, udotb, err
C-----------------------------------------------
 
c       Inverts Square Matrix A using SVD method
c       Allows supression of weights, storage and recall of answers
c       Uses svd routines from numerical recipes
c       Prashant Valanju (Sept 1998) pvalanju@mail.utexas.edu
 
c       A(M,M) - Matrix on input, Inverse on output
c       M      - Size of A
c       Wcut   - Cutoff parameter
c                If Wcut = 0, writes answers into "file", all weights are kept
c                If Wcut < 0, reads previous answers from "file",
c                              and uses -Wcut for weights to keep
c                If 0 < Wcut < 1 : Cuts off weights with w(i)/wmax < wcut
c                      Use this to cut off small weights, does calculation
c                If Wcut = integer > 1, keeps Wcut weights, does calculation
c                      Use this to cut off fixed number of weights
c
c       Note: M is same as N, Mp, and Np in this routine
 
c     allocate local arrays
      allocate (u(m,m), w(m), v(m,m), stat=istat)
      if (istat .eq. 0) then                       !have enough memory

         iunit = 41
         wcuta = abs(wcut)
c...............................
c     If Wcut < 0 use precalculated weights to get answer fast
c        use wcuta = abs(wcut) for everything
c     If any error happens while reading "file", do full svd calculation
         if (wcut .lt. 0) then         !Read file from previous calculation
            call safe_open(iunit, istat, file, 'old', 'formatted')
            if (istat .ne. 0) goto 98
            read (iunit, *, err=98)
            read (iunit, *, err=98)
            read (iunit, *, err=98) m, nonzero
c                           !Read exactly the same way they were written
            do i = 1, nonzero
               read (iunit, *, err=98) w(i)         !read i-th weight
               do j = 1, m
                  read (iunit, *, err=98) v(j,i), u(j,i)
               end do
            end do
            close(unit=iunit)
 
c        Decide how many weights to keep
            nlast = wcuta                        !remember Wcut is < 0
            nlast = min(nonzero,nlast)      !Do not use nonzero weights
 
            go to 99           !bypass svdcmp and go to matrix inversion
         endif               !End of precalculated U,w,V branch Wcut < 0
   98    continue
c.......................................
c     If Wcut .ge. 0, do full svd calculation and save V,w,U
c     Initialize all to zero to wipe out effects of old call
         do j = 1, m
            w(j) = zero                           !Zero all weights
            u(:m,j) = a(:m,j)  !Because U will be changed by svdcmp
            v(:m,j) = zero
         end do
 
         call svdcmp (u, m, m, m, m, w, v)       !Do SVD decomposition
 
c       Sort weights and matrices with DECREASING weights
c       Permute weight w(i) AND column vectors U(*,i), V(*,i) at the same time
         call sortsvd (m, m, m, m, w, u, v)
 
c       Find the number of nonzero weights (already dcreasing ordered)
         do nonzero = m, 1, -1
c                                   !Found first nonzero weight, get out
            if (w(nonzero) .ne. 0) exit 
         end do
         if (nonzero .le. 0) go to 999
 
c       Write all weights and U, V into SVD file (if Wcut=0)
         if (wcut .eq. 0) then    !no nonzero weights, you must be kidding
            call safe_open(iunit, istat, file, 'unknown', 'formatted')
            write (iunit, *) 'Max w = ', w(1), ', Min w = ', w(nonzero)
            write (iunit, *) 'Ratio = ', w(1)/w(nonzero)
            write (iunit, *) m, nonzero
c                          !write exactly the same way they were written
            do i = 1, nonzero
               write (iunit, *) w(i)                !write i-th weight
               do j = 1, m
                  write (iunit, *) v(j,i), u(j,i)
               end do
            end do
            close(unit=iunit)
         endif
c.......................................
 
c       Following gets done in both cases (using saved or recalculating)
c        Decide how many weights to keep
   99    continue
         if (wcuta .gt. 1) then
            nlast = wcuta
         else                                    !cutoff small weights
            wmax = w(1)                     !weights are already ordered
            nlast = 0
            do i = 1, nonzero
               if (w(i) .gt. wmax*wcuta) then
                  nlast = nlast + 1              !accept this weight
               else
                  go to 96
               endif
            end do
         endif
 
c       Find Inverse of A = V [diag(1/wi)] Utranspose, return in A
c       First do [diag(1/wi)] Utranspose, store answer in U again
   96    continue
         if (nlast.le.0 .or. nlast.gt.nonzero) nlast = nonzero
         do i = 1, nlast
c                                         !divide ith row of Utr by w(i)
            u(:m,i) = u(:m,i)/w(i)
         end do
c        Zero the infinite 1/weights
         do i = nlast + 1, m
            u(:m,i) = zero                !multiply ith row of Utr by zero
         end do
c       Next multiply by V matrix to get A inverse, put it in A
         do i = 1, m                    !multiplt ith row vector of V by
            do j = 1, m                         !jth column vector of wU
                                                !vector multiply V and wU
               a(i,j) = sum(v(i,:m)*u(j,:m))    !Put inverse back in A
            end do
         end do
c................................................
      endif
  999 continue

      deallocate (u, w, v, stat=istat)

!     end subroutine svdinv
      end
 

      SUBROUTINE SVDCMP(A, M, N, MP, NP, W, V)
      use kind_spec
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER M, N, MP, NP
      REAL(rprec), DIMENSION(MP,NP) :: A
      REAL(rprec), DIMENSION(NP) :: W
      REAL(rprec), DIMENSION(NP,NP) :: V
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec) :: zero = 0, one = 1, two = 2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I, L, K, J, ITS, NM
      REAL(rprec), DIMENSION(MP*NP) :: RV1
      REAL(rprec) :: G, SCALE, ANORM, S, F, H, C, Y, Z, X
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      REAL(rprec) , EXTERNAL :: PYTHAG
C-----------------------------------------------
      IF (M .GT. MP) STOP 'M > MP IN SVDCMP'
      IF (N .GT. NP) STOP 'N > NP IN SVDCMP'

      G = zero
      SCALE = zero
      ANORM = zero
      DO I = 1, N
         L = I + 1
         RV1(I) = SCALE*G
         G = zero
         S = zero
         SCALE = zero
         IF (I .le. M) THEN
            SCALE = SUM(ABS(A(I:M,I)))
            IF (SCALE .ne. zero) THEN
               A(I:M,I) = A(I:M,I)/SCALE
               S = SUM(A(I:M,I)*A(I:M,I))
               F = A(I,I)
               G = -SIGN(SQRT(S),F)
               H = F*G - S
               A(I,I) = F - G
               IF (I .ne. N) THEN
                  DO J = L, N
                     S = SUM(A(I:M,I)*A(I:M,J))
                     F = S/H
                     A(I:M,J) = A(I:M,J) + F*A(I:M,I)
                  END DO
               ENDIF
               A(I:M,I) = SCALE*A(I:M,I)
            ENDIF
         ENDIF
         W(I) = SCALE*G
         G = zero
         S = zero
         SCALE = zero
         IF (I.le.M .AND. I.ne.N) THEN
            SCALE = SUM(ABS(A(I,L:N)))
            IF (SCALE .ne. zero) THEN
               A(I,L:N) = A(I,L:N)/SCALE
               S = SUM(A(I,L:N)*A(I,L:N))
               F = A(I,L)
               G = -SIGN(SQRT(S),F)
               H = F*G - S
               A(I,L) = F - G
               RV1(L:N) = A(I,L:N)/H
               IF (I .ne. M) THEN
                  DO J = L, M
                     S = SUM(A(J,L:N)*A(I,L:N))
                     A(J,L:N) = A(J,L:N) + S*RV1(L:N)
                  END DO
               ENDIF
               A(I,L:N) = SCALE*A(I,L:N)
            ENDIF
         ENDIF
         ANORM = MAX(ANORM,ABS(W(I))+ABS(RV1(I)))
      END DO
      DO I = N, 1, -1
         IF (I .lt. N) THEN
            IF (G .ne. zero) THEN
               V(L:N,I) = (A(I,L:N)/A(I,L))/G
               DO J = L, N
                  S = SUM(A(I,L:N)*V(L:N,J))
                  V(L:N,J) = V(L:N,J) + S*V(L:N,I)
               END DO
            ENDIF
            V(I,L:N) = zero
            V(L:N,I) = zero
         ENDIF
         V(I,I) = one
         G = RV1(I)
         L = I
      END DO
      DO I = N, 1, -1
         L = I + 1
         G = W(I)
         IF (I .lt. N) THEN
            A(I,L:N) = zero
         ENDIF
         IF (G .ne. zero) THEN
            G = one/G
            IF (I .ne. N) THEN
               DO J = L, N
                  S = SUM(A(L:M,I)*A(L:M,J))
                  F = (S/A(I,I))*G
                  A(I:M,J) = A(I:M,J) + F*A(I:M,I)
               END DO
            ENDIF
            A(I:M,I) = A(I:M,I)*G
         ELSE
            A(I:M,I) = zero
         ENDIF
         A(I,I) = A(I,I) + one
      END DO
      L49: DO K = N, 1, -1
         DO ITS = 1, 50
            DO L = K, 1, -1
               NM = L - 1
               IF (ABS(RV1(L)) + ANORM .eq. ANORM) GO TO 2
               IF (ABS(W(NM)) + ANORM .eq. ANORM) EXIT 
            END DO
            C = zero
            S = one
            DO I = L, K
               F = S*RV1(I)
               IF (ABS(F) + ANORM .ne. ANORM) THEN
                  G = W(I)
                  H = PYTHAG(F,G)
c              H=SQRT(F*F+G*G)
                  W(I) = H
                  H = one/H
                  C = G*H
                  S = -F*H
                  DO J = 1, M
                     Y = A(J,NM)
                     Z = A(J,I)
                     A(J,NM) = Y*C + Z*S
                     A(J,I) = (-Y*S) + Z*C
                  END DO
               ENDIF
            END DO
    2       CONTINUE
            Z = W(K)
            IF (L .eq. K) THEN
               IF (Z .lt. zero) THEN
                  W(K) = -Z
                  V(:N,K) = -V(:N,K)
               ENDIF
               CYCLE  L49
            ENDIF
            IF (ITS .eq. 50) THEN
               WRITE (*, '(2A)') 'PAUSE ', 
     1            'No convergence in 150 iterations'
               READ *
            ENDIF
            X = W(L)
            NM = K - 1
            Y = W(NM)
            G = RV1(NM)
            H = RV1(K)
            F = ((Y - Z)*(Y + Z) + (G - H)*(G + H))/(TWO*H*Y)
            G = PYTHAG(F,one)
c          G=SQRT(F*F + one)
            F = ((X - Z)*(X + Z) + H*(Y/(F + SIGN(G,F)) - H))/X
            C = one
            S = one
            DO J = L, NM
               I = J + 1
               G = RV1(I)
               Y = W(I)
               H = S*G
               G = C*G
               Z = PYTHAG(F,H)
c            Z=SQRT(F*F+H*H)
               RV1(J) = Z
               C = F/Z
               S = H/Z
               F = X*C + G*S
               G = (-X*S) + G*C
               H = Y*S
               Y = Y*C
               DO NM = 1, N
                  X = V(NM,J)
                  Z = V(NM,I)
                  V(NM,J) = X*C + Z*S
                  V(NM,I) = (-X*S) + Z*C
               END DO
               Z = PYTHAG(F,H)
c            Z=SQRT(F*F+H*H)
               W(J) = Z
               IF (Z .ne. zero) THEN
                  Z = one/Z
                  C = F*Z
                  S = H*Z
               ENDIF
               F = C*G + S*Y
               X = (-S*G) + C*Y
               DO NM = 1, M
                  Y = A(NM,J)
                  Z = A(NM,I)
                  A(NM,J) = Y*C + Z*S
                  A(NM,I) = (-Y*S) + Z*C
               END DO
            END DO
            RV1(L) = zero
            RV1(K) = F
            W(K) = X
         END DO
      END DO L49

!     END SUBROUTINE SVDCMP
      END


      SUBROUTINE SVBKSB(U, W, V, M, N, MP, NP, B, X)
      use kind_spec
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER M, N, MP, NP
      REAL(rprec), DIMENSION(MP,NP) :: U
      REAL(rprec), DIMENSION(NP) :: W
      REAL(rprec), DIMENSION(NP,NP) :: V
      REAL(rprec), DIMENSION(MP) :: B
      REAL(rprec), DIMENSION(NP) :: X
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec) :: zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: J, I, JJ
      REAL(rprec), DIMENSION(MP*NP) :: TMP
      REAL(rprec) :: S
C-----------------------------------------------
      IF (M .GT. MP) STOP 'M > MP IN SVBKSB'
      IF (N .GT. NP) STOP 'N > NP IN SVBKSB'

      DO J = 1, N
         S = zero
         IF (W(J) .ne. zero) THEN
            S = SUM(U(:M,J)*B(:M))
            S = S/W(J)
         ENDIF
         TMP(J) = S
      END DO
      DO J = 1, N
         S = SUM(V(J,:N)*TMP(:N))
         X(J) = S
      END DO

!     END SUBROUTINE SVBKSB
      END

      SUBROUTINE SORTSVD(M, N, MP, NP, W, U, V)
c------------------------------
c     Sorts the weights and U, V matrices from SVD decomposition
c     Sorts weights in DECREASING order
c     Based on routine SORT2 from Numerical Recipes, pg 231
c     Created by Prashant Valanju, (Sept 1998)
c------------------------------
      use kind_spec
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER M, N, MP, NP
      REAL(rprec), DIMENSION(NP) :: W
      REAL(rprec), DIMENSION(MP,NP) :: U
      REAL(rprec), DIMENSION(NP,NP) :: V
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: L, IR, II, I, J
      REAL(rprec), DIMENSION(MP) :: RU
      REAL(rprec), DIMENSION(NP) :: RV
      REAL(rprec) :: RW
C-----------------------------------------------
      IF (M .GT. MP) STOP 'M > MP IN SORTSVD'
      IF (N .GT. NP) STOP 'N > NP IN SORTSVD'

      L = N/2 + 1
      IR = N
   10 CONTINUE
      IF (L .gt. 1) THEN
         L = L - 1
         RW = W(L)
         RU(:M) = U(:M,L)
         RV(:N) = V(:N,L)
      ELSE
         RW = W(IR)
         RU(:M) = U(:M,IR)
         RV(:N) = V(:N,IR)
         W(IR) = W(1)
         U(:M,IR) = U(:M,1)
         V(:N,IR) = V(:N,1)
         IR = IR - 1
         IF (IR .eq. 1) THEN
            W(1) = RW
            U(:M,1) = RU(:M)
            V(:N,1) = RV(:N)
            RETURN 
         ENDIF
      ENDIF
      I = L
      J = L + L
   20 CONTINUE
      IF (J .le. IR) THEN
         IF (J .lt. IR) THEN
C                                    !Just change these 2 for increasing
            IF (W(J) .gt. W(J+1)) J = J + 1
         ENDIF
         IF (RW .gt. W(J)) THEN         !Just change these 2 for increasing
            W(I) = W(J)
            U(:M,I) = U(:M,J)
            V(:N,I) = V(:N,J)
            I = J
            J = J + J
         ELSE
            J = IR + 1
         ENDIF
         GO TO 20
      ENDIF
      W(I) = RW
      U(:M,I) = RU(:M)
      V(:N,I) = RV(:N)
      GO TO 10

!     END SUBROUTINE SORTSVD
      END


      function pythag (a, b)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec) a, b
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec) :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec) :: absa, absb, pythag
C-----------------------------------------------
c computes sqrt(a^2+b^2) without destructive underflow or overflow
      absa = abs(a)
      absb = abs(b)
      if (absa .gt. absb) then
         pythag = absa*sqrt(one + (absb/absa)**2)
      else
         if (absb .eq. zero) then
            pythag = zero
         else
            pythag = absb*sqrt(one + (absa/absb)**2)
         endif
      endif

!     end function pythag
      end
EOF
EOC
