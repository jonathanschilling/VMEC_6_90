      SUBROUTINE SCAL_G (N, SA, SX, INCX)
      USE KIND_SPEC
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER, INTENT(IN) :: N, INCX
      REAL(RPREC), INTENT(IN) :: SA
      REAL(RPREC), INTENT(INOUT) :: SX(*)
C-----------------------------------------------
C   E X T E R N A L   F U N C T I O N S
C-----------------------------------------------
      EXTERNAL SCAL_F90
C-----------------------------------------------
!     GENERIC INTERFACE TO SINGLE OR DOUBLE PRECISION SCALE ROUTINES

      CALL SCAL_F90 (N, SA, SX, INCX)

      END SUBROUTINE SCAL_G


      SUBROUTINE AXPY_G(N, SA, SX, INCX, SY, INCY)
      USE KIND_SPEC
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER :: N, INCX, INCY
      REAL(RPREC) :: SA, SX(*), SY(*)
C-----------------------------------------------
C   E X T E R N A L   F U N C T I O N S
C-----------------------------------------------
      EXTERNAL AXPY_F90
C-----------------------------------------------
!     GENERIC INTERFACE TO SINGLE, DOUBLE PRECISION AX + Y ROUTINES

      CALL AXPY_F90 (N, SA, SX, INCX, SY, INCY)

      END SUBROUTINE AXPY_G


      FUNCTION DOT_G (N, SX, INCX, SY, INCY)
      USE KIND_SPEC
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER :: N, INCX, INCY
      REAL(RPREC) :: SX(*), SY(*)
      REAL(RPREC) :: DOT_G
C-----------------------------------------------
C   E X T E R N A L   F U N C T I O N S
C-----------------------------------------------
      REAL(RPREC), EXTERNAL :: DOT_F90
C-----------------------------------------------
!     GENERIC INTERFACE TO SINGLE, DOUBLE PRECISION DOT-PRODUCT ROUTINES

      DOT_G = DOT_F90 (N, SX, INCX, SY, INCY)

      END FUNCTION DOT_G


      SUBROUTINE COPY_G (N, SX, INCX, SY, INCY)
      USE KIND_SPEC
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER N, INCX, INCY
      REAL(RPREC), DIMENSION(*) :: SX, SY
C-----------------------------------------------
C   E X T E R N A L   F U N C T I O N S
C-----------------------------------------------
      EXTERNAL COPY_F90
C-----------------------------------------------

!     GENERIC INTERFACE TO SINGLE, DOUBLE PRECISION COPY ROUTINES

      CALL COPY_F90 (N, SX, INCX, SY, INCY)

      END SUBROUTINE COPY_G


      SUBROUTINE SWAP_G (N, SX, INCX, SY, INCY)
      USE KIND_SPEC
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER :: N, INCX, INCY
      REAL (RPREC) :: SX(*), SY(*)
C-----------------------------------------------
C   E X T E R N A L   F U N C T I O N S
C-----------------------------------------------
      EXTERNAL SWAP_F90
C-----------------------------------------------

!     GENERIC INTERFACE TO SINGLE, DOUBLE PRECISION SWAP ROUTINES

      CALL SWAP_F90 (N, SX, INCX, SY, INCY)

      END SUBROUTINE SWAP_G


      SUBROUTINE TRSM_G(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA,
     1    A, LDA, B, LDB)
      USE KIND_SPEC
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER M, N, LDA, LDB
      REAL(RPREC) ALPHA
      CHARACTER SIDE, UPLO, TRANSA, DIAG
      REAL(RPREC), DIMENSION(LDA,*) :: A
      REAL(RPREC), DIMENSION(LDB,*) :: B
C-----------------------------------------------
C   E X T E R N A L   F U N C T I O N S
C-----------------------------------------------
      EXTERNAL TRSM_F90
C-----------------------------------------------

!     GENERIC INTERFACE TO SINGLE, DOUBLE PRECISION TRSM ROUTINES

      CALL TRSM_F90 (SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA,
     1    A, LDA, B, LDB)

      END SUBROUTINE TRSM_G


      SUBROUTINE GETF2_G(M, N, A, LDA, IPIV, INFO)
      USE KIND_SPEC
      IMPLICIT NONE
!-----------------------------------------------
!   D U M M Y   A R G U M E N T S
!-----------------------------------------------
      INTEGER M, N, LDA, INFO
      INTEGER, DIMENSION(*) :: IPIV
      REAL(RPREC), DIMENSION(LDA,*) :: A
C-----------------------------------------------
C   E X T E R N A L   F U N C T I O N S
C-----------------------------------------------
      EXTERNAL GETF2_F90
C-----------------------------------------------
!     GENERIC INTERFACE TO SINGLE, DOUBLE PRECISION GETF2 ROUTINES

      CALL GETF2_F90 (M, N, A, LDA, IPIV, INFO)

      END SUBROUTINE GETF2_G


      SUBROUTINE GER_G(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
      USE KIND_SPEC
      IMPLICIT NONE
!-----------------------------------------------
!   D U M M Y   A R G U M E N T S
!-----------------------------------------------
      INTEGER :: M,N,INCX,INCY,LDA
      REAL(RPREC) :: ALPHA, X(*), Y(*), A(LDA,*)
C-----------------------------------------------
C   E X T E R N A L   F U N C T I O N S
C-----------------------------------------------
      EXTERNAL GER_F90
C-----------------------------------------------
!     GENERIC INTERFACE TO SINGLE, DOUBLE PRECISION GETF2 ROUTINES

      CALL GER_F90 (M,N,ALPHA,X,INCX,Y,INCY,A,LDA)

      END SUBROUTINE GER_G


      SUBROUTINE LASWP_G(N, A, LDA, K1, K2, IPIV, INCX)
      USE KIND_SPEC
      IMPLICIT NONE
!-----------------------------------------------
!   D U M M Y   A R G U M E N T S
!-----------------------------------------------
      INTEGER N, LDA, K1, K2, INCX
      INTEGER, DIMENSION(*) :: IPIV
      REAL(RPREC), DIMENSION(LDA,*) :: A
C-----------------------------------------------
C   E X T E R N A L   F U N C T I O N S
C-----------------------------------------------
      EXTERNAL LASWP_F90
C-----------------------------------------------

      CALL LASWP_F90(N, A, LDA, K1, K2, IPIV, INCX)

      END SUBROUTINE LASWP_G


      FUNCTION LAMCH_G( CMACH )
      USE KIND_SPEC
      IMPLICIT NONE
!-----------------------------------------------
!   D U M M Y   A R G U M E N T S
!-----------------------------------------------
      REAL(RPREC) :: LAMCH_G
      CHARACTER(LEN=1), INTENT(IN) :: CMACH
C-----------------------------------------------
C   E X T E R N A L   F U N C T I O N S
C-----------------------------------------------
      REAL(RPREC), EXTERNAL :: LAMCH_F90
C-----------------------------------------------

      LAMCH_G = LAMCH_F90(CMACH)

      END FUNCTION LAMCH_G


      SUBROUTINE GEMM_G( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA,
     1                    B, LDB, BETA, C, LDC )
      USE KIND_SPEC
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      CHARACTER*1         TRANSA, TRANSB
      INTEGER             M, N, K, LDA, LDB, LDC
      REAL(RPREC) :: ALPHA, BETA
      REAL(RPREC) :: A( LDA, * ), B( LDB, * ), C( LDC, * )
C-----------------------------------------------
C   E X T E R N A L   F U N C T I O N S
C-----------------------------------------------
      EXTERNAL GEMM_F90
C-----------------------------------------------
!     GENERIC INTERFACE TO SINGLE, DOUBLE PRECISION GEMM ROUTINES

      CALL GEMM_F90 (TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
     1                BETA, C, LDC)

      END SUBROUTINE GEMM_G


      SUBROUTINE GETRF_G(M, N, A, LDA, IPIV, INFO)
      USE KIND_SPEC
      IMPLICIT NONE
!-----------------------------------------------
!   D U M M Y   A R G U M E N T S
!-----------------------------------------------
      INTEGER M, N, LDA, INFO
      INTEGER, DIMENSION(*) :: IPIV
      REAL(RPREC), DIMENSION(LDA,*) :: A
C-----------------------------------------------
C   E X T E R N A L   F U N C T I O N S
C-----------------------------------------------
      EXTERNAL GETRF_F90
C-----------------------------------------------

!     GENERIC INTERFACE TO SINGLE, DOUBLE PRECISION GETRF ROUTINES

      CALL GETRF_F90(M, N, A, LDA, IPIV, INFO)

      END SUBROUTINE GETRF_G


      SUBROUTINE GETRS_G(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)
      USE KIND_SPEC
      IMPLICIT NONE
!-----------------------------------------------
!   D U M M Y   A R G U M E N T S
!-----------------------------------------------
      INTEGER N, NRHS, LDA, LDB, INFO
      CHARACTER TRANS
      INTEGER, DIMENSION(*) :: IPIV
      REAL(RPREC), DIMENSION(LDA,*) :: A
      REAL(RPREC), DIMENSION(LDB,*) :: B
C-----------------------------------------------
C   E X T E R N A L   F U N C T I O N S
C-----------------------------------------------
      EXTERNAL GETRS_F90
C-----------------------------------------------

!     GENERIC INTERFACE TO SINGLE, DOUBLE PRECISION GETRS ROUTINES

      CALL GETRS_F90 (TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)

      END SUBROUTINE GETRS_G


      SUBROUTINE GEFA_G (A, LDA, N, IPVT, INFO)
      USE KIND_SPEC
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER LDA, N, INFO
      INTEGER, DIMENSION(N) :: IPVT
      REAL(RPREC), DIMENSION(LDA,N) :: A
C-----------------------------------------------

      CALL GEFA_F90 (A, LDA, N, IPVT, INFO)

      END SUBROUTINE GEFA_G


      SUBROUTINE GESL_G (A, LDA, N, IPVT, B, JOB)
      USE KIND_SPEC
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER LDA, N, JOB
      INTEGER, DIMENSION(N) :: IPVT
      REAL(RPREC), DIMENSION(LDA,N) :: A
      REAL(RPREC), DIMENSION(N) :: B
C-----------------------------------------------

      CALL GESL_F90 (A, LDA, N, IPVT, B, JOB)

      END SUBROUTINE GESL_G


      SUBROUTINE GESV_G (N, NRHS, A, LDA, IPIV, B, LDB, INFO)
      USE KIND_SPEC
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER :: INFO, LDA, LDB, N, NRHS
      INTEGER :: IPIV(N)
      REAL(RPREC) :: A(LDA, N), B(LDB, NRHS)
C-----------------------------------------------
      CALL GESV_F90 (N, NRHS, A, LDA, IPIV, B, LDB, INFO)

      END SUBROUTINE GESV_G


      SUBROUTINE GBFA_G (ABD, LDA, N, ML, MU, IPVT, INFO)
      USE KIND_SPEC
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER LDA, N, ML, MU, INFO
      INTEGER, DIMENSION(N) :: IPVT
      REAL(RPREC), DIMENSION(LDA,N) :: ABD
C-----------------------------------------------
      CALL GBFA_F90(ABD, LDA, N, ML, MU, IPVT, INFO)

      END SUBROUTINE GBFA_G


      SUBROUTINE GBSL_G (ABD, LDA, N, ML, MU, IPVT, B, JOB)
      USE KIND_SPEC
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER LDA, N, ML, MU, JOB
      INTEGER, DIMENSION(N) :: IPVT
      REAL(RPREC), DIMENSION(LDA,N) :: ABD
      REAL(RPREC), DIMENSION(N) :: B
C-----------------------------------------------
      CALL GBSL_F90 (ABD, LDA, N, ML, MU, IPVT, B, JOB)

      END SUBROUTINE GBSL_G


      SUBROUTINE GEFA_F90(A, LDA, N, IPVT, INFO)
      USE KIND_SPEC
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER LDA, N, INFO
      INTEGER, DIMENSION(N) :: IPVT
      REAL(RPREC), DIMENSION(LDA,N) :: A
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      REAL(RPREC), PARAMETER :: ZERO = 0, ONE = 1
      INTEGER :: J, K, KP1, L, NM1
      INTEGER, DIMENSION(1) :: ISAMAX
      REAL(RPREC) :: ELEMENT
C-----------------------------------------------
c
c     sgefa factors a real matrix by gaussian elimination.
c
c     sgefa is usually called by sgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for sgeco) = (1 + 9/n)*(time for sgefa) .
c
c     on entry
c
c        a       real(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k).eq.0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that sgesl or sgedi will divide by zero
c                     if called.  use  rcond  in sgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .ge. 1) then
         do k = 1, nm1
            kp1 = k + 1
c
c         find l = pivot index
c
            isamax = maxloc(abs(a(k:n,k)))
            l = isamax(1) + k - 1
C           l = IDAMAX(n-k+1,a(k,k),1) + k - 1
            ipvt(k) = l
c

c         zero pivot implies this column already triangularized
c
            if (a(l,k) .ne. zero) then
c
c           interchange if necessary
c
               if (l .ne. k) then
                  element = a(l,k)
                  a(l,k) = a(k,k)
                  a(k,k) = element
               endif
c
c         compute multipliers
c
               element = -one/a(k,k)
!              call scal_g(n-k,element,a(k+1,k),1)
               a(k+1:n,k) = element*a(k+1:n,k)
c
c           row elimination with column indexing
c
               do j = kp1, n
                 element = a(l,j)
                 if (l .ne. k) then
                    a(l,j) = a(k,j)
                    a(k,j) = element
                 end if
!                call axpy_g(n-k,element,a(k+1,k),1,a(k+1,j),1)
                 a(k+1:n,j) = a(k+1:n,j) + element*a(k+1:n,k)
               end do

            else
               info = k
            endif
         end do
      endif
      ipvt(n) = n
      if (a(n,n) .eq. zero) info = n

      END SUBROUTINE GEFA_F90


      SUBROUTINE GESL_F90(A, LDA, N, IPVT, B, JOB)
      USE KIND_SPEC
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER LDA, N, JOB
      INTEGER, DIMENSION(N) :: IPVT
      REAL(RPREC), DIMENSION(LDA,N) :: A
      REAL(RPREC), DIMENSION(N) :: B
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      INTEGER :: K, L, NM1
      REAL(RPREC) :: ELEMENT
C-----------------------------------------------
c
c     sgesl solves the real system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by sgeco or sgefa.
c
c     on entry
c
c        a       real(lda, n)
c                the output from sgeco or sgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from sgeco or sgefa.
c
c        b       real(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if sgeco has set rcond.gt.0.0
c        or sgefa has set info.eq.0 .
c
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c
      nm1 = n - 1
      if (job .eq. 0) then
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .ge. 1) then
            do k = 1, nm1
               l = ipvt(k)
               element = b(l)
               if (l .ne. k) then
                  b(l) = b(k)
                  b(k) = element
               endif
!              call axpy_g(n-k,element,a(k+1,k),1,b(k+1),1)
               b(k+1:n) = b(k+1:n) + element*a(k+1:n,k)
            end do
         endif
c
c        now solve  u*x = y
c
         do k = n, 1, -1
            b(k) = b(k)/a(k,k)
            element = -b(k)
!           call axpy_g(k-1,element,a(1,k),1,b,1)
            b(1:k-1) = b(1:k-1) + element*a(1:k-1,k)
         end do
      else
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do k = 1, n
            element = sum(a(:k-1,k)*b(:k-1))
            b(k) = (b(k)-element)/a(k,k)
         end do
c
c        now solve trans(l)*x = y
c
         if (nm1 .ge. 1) then
            do k = nm1, 1, -1
               b(k) = b(k) + sum(a(k+1:n,k)*b(k+1:n))
               l = ipvt(k)
               if (l .ne. k) then
                  element = b(l)
                  b(l) = b(k)
                  b(k) = element
               endif
            end do
         endif
      endif

      END SUBROUTINE GESL_F90


      SUBROUTINE GBFA_F90 (ABD, LDA, N, ML, MU, IPVT, INFO)
      USE KIND_SPEC
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER LDA, N, ML, MU, INFO
      INTEGER, DIMENSION(N) :: IPVT
      REAL(RPREC), DIMENSION(LDA,N) :: ABD
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      REAL(RPREC), PARAMETER :: ZERO = 0, ONE = 1
      INTEGER :: I0, J, JU, JZ, J0, J1, K, KP1, L, LM, M, MM, NM1
      REAL(RPREC) :: T
      INTEGER :: ISAMAX1(1)
C-----------------------------------------------
c
c     sgbfa factors a real band matrix by elimination.
c
c     sgbfa is usually called by sgbco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c
c     on entry
c
c        abd     real(lda, n)
c                contains the matrix in band storage.  the columns
c                of the matrix are stored in the columns of  abd  and
c                the diagonals of the matrix are stored in rows
c                ml+1 through 2*ml+mu+1 of  abd .
c                see the comments below for details.
c
c        lda     integer
c                the leading dimension of the array  abd .
c                lda must be .ge. 2*ml + mu + 1 .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c                0 .le. ml .lt. n .
c
c        mu      integer
c                number of diagonals above the main diagonal.
c                0 .le. mu .lt. n .
c                more efficient if  ml .le. mu .
c     on return
c
c        abd     an upper triangular matrix in band storage and
c                the multipliers which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that sgbsl will divide by zero if
c                     called.  use  rcond  in sgbco for a reliable
c                     indication of singularity.
c
c     band storage
c
c           if  a  is a band matrix, the following program segment
c           will set up the input.
c
c                   ml = (band width below the diagonal)
c                   mu = (band width above the diagonal)
c                   m = ml + mu + 1
c                   do 20 j = 1, n
c                      i1 = max(1, j-mu)
c                      i2 = min(n, j+ml)
c                      do 10 i = i1, i2
c                         k = i - j + m
c                         abd(k,j) = a(i,j)
c                10    continue
c                20 continue
c
c           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
c           in addition, the first  ml  rows in  abd  are used for
c           elements generated during the triangularization.
c           the total number of rows needed in  abd  is  2*ml+mu+1 .
c           the  ml+mu by ml+mu  upper left triangle and the
c           ml by ml  lower right triangle are not referenced.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas saxpy ,sscal
c     fortran max,min
c
c     internal variables
c
c
c
      m = ml + mu + 1
      info = 0
c
c     zero initial fill-in columns
c
      j0 = mu + 2
      j1 = min(n,m) - 1
      do jz = j0, j1
         i0 = m + 1 - jz
         abd(i0:ml,jz) = zero
      end do
      jz = j1
      ju = 0
c
c     gaussian elimination with partial pivoting
c
      nm1 = n - 1
      if (nm1 .ge. 1) then
         do k = 1, nm1
            kp1 = k + 1
c
c        zero next fill-in column
c
            jz = jz + 1
            if (jz .le. n) then
               abd(:ml,jz) = zero
            endif
c
c        find l = pivot index
c
            lm = min(ml,n - k)
            isamax1 = maxloc(abs(abd(m:m+lm,k)))
            l = isamax1(1) + m - 1
!           l = isamax(lm + 1,abd(m,k),1) + m - 1
            ipvt(k) = l + k - m
c
c        zero pivot implies this column already triangularized
c
            if (abd(l,k) .ne. zero) then
c
c           interchange if necessary
c
               if (l .ne. m) then
                  t = abd(l,k)
                  abd(l,k) = abd(m,k)
                  abd(m,k) = t
               endif
c
c           compute multipliers
c
               t = -one/abd(m,k)
               call scal_g (lm, t, abd(m+1,k), 1)
c
c           row elimination with column indexing
c
               ju = min(max(ju,mu + ipvt(k)),n)
               mm = m
               do j = kp1, ju
                  l = l - 1
                  mm = mm - 1
                  t = abd(l,j)
                  if (l .ne. mm) then
                     abd(l,j) = abd(mm,j)
                     abd(mm,j) = t
                  endif
                  call axpy_g (lm, t, abd(m+1,k), 1, abd(mm+1,j), 1)
               end do
            else
               info = k
            endif
         end do
      endif
      ipvt(n) = n
      if (abd(m,n) .eq. zero) info = n

      END SUBROUTINE GBFA_F90


      SUBROUTINE GBSL_F90 (ABD, LDA, N, ML, MU, IPVT, B, JOB)
      USE KIND_SPEC
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER LDA, N, ML, MU, JOB
      INTEGER, DIMENSION(N) :: IPVT
      REAL(RPREC), DIMENSION(LDA,N) :: ABD
      REAL(RPREC), DIMENSION(N) :: B
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      INTEGER :: K, KB, L, LA, LB, LM, M, NM1
      REAL(RPREC) :: T
C-----------------------------------------------
C   E X T E R N A L   F U N C T I O N S
C-----------------------------------------------
      REAL(RPREC) , EXTERNAL :: DOT_G
      EXTERNAL AXPY_G
C-----------------------------------------------
c
c     sgbsl solves the real band system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by sgbco or sgbfa.
c
c     on entry
c
c        abd     real(lda, n)
c                the output from sgbco or sgbfa.
c
c        lda     integer
c                the leading dimension of the array  abd .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c
c        mu      integer
c                number of diagonals above the main diagonal.
c
c        ipvt    integer(n)
c                the pivot vector from sgbco or sgbfa.
c
c        b       real(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b , where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if sgbco has set rcond .gt. 0.0
c        or sgbfa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call sgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call sgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas saxpy ,sdot
c     fortran min
c
c     internal variables
c
c
      m = mu + ml + 1
      nm1 = n - 1
      if (job .eq. 0) then
c
c        job = 0 , solve  a * x = b
c        first solve l*y = b
c
         if (ml .ne. 0) then
            do k = 1, nm1
               lm = min(ml,n - k)
               l = ipvt(k)
               t = b(l)
               if (l .ne. k) then
                  b(l) = b(k)
                  b(k) = t
               endif
               call axpy_g (lm, t, abd(m+1,k), 1, b(k+1), 1)
            end do
         endif
c
c        now solve  u*x = y
c
         do kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/abd(m,k)
            lm = min(k,m) - 1
            la = m - lm
            lb = k - lm
            t = -b(k)
            call axpy_g (lm, t, abd(la,k), 1, b(lb), 1)
         end do
      else
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do k = 1, n
            lm = min(k,m) - 1
            la = m - lm
            lb = k - lm
            t = dot_g(lm,abd(la,k),1,b(lb),1)
            b(k) = (b(k)-t)/abd(m,k)
         end do
c
c        now solve trans(l)*x = y
c
         if (ml .ne. 0) then
            do kb = 1, nm1
               k = n - kb
               lm = min(ml,n - k)
               b(k) = b(k) + dot_g(lm,abd(m+1,k),1,b(k+1),1)
               l = ipvt(k)
               if (l .ne. k) then
                  t = b(l)
                  b(l) = b(k)
                  b(k) = t
               endif
            end do
         endif
      endif

      END SUBROUTINE GBSL_F90



      SUBROUTINE GETRF_F90(M, N, A, LDA, IPIV, INFO)
      USE KIND_SPEC
      IMPLICIT NONE
!-----------------------------------------------
!   D U M M Y   A R G U M E N T S
!-----------------------------------------------
      INTEGER, INTENT(IN) :: LDA, M, N
      INTEGER, INTENT(OUT) :: INFO
      INTEGER, INTENT( OUT ) :: IPIV( * )
      REAL(RPREC), INTENT( INOUT ) :: A( LDA, * )
!-----------------------------------------------
!   L O C A L   P A R A M E T E R S
!-----------------------------------------------
      REAL(RPREC), PARAMETER :: ONE = 1
!-----------------------------------------------
!   L O C A L   V A R I A B L E S
!-----------------------------------------------
      INTEGER :: IINFO, J, JB, NB
!-----------------------------------------------
!   E X T E R N A L   F U N C T I O N S
!-----------------------------------------------
      INTEGER , EXTERNAL :: ILAENV
      EXTERNAL GEMM_G, GETF2_G, LASWP_G, TRSM_G
!-----------------------------------------------
!
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  SGETRF computes an LU factorization of a general M-by-N matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 3 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) REAL array, dimension (LDA,N)
!          On entry, the M-by-N matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!                has been completed, but the factor U is exactly
!                singular, and division by zero will occur if it is used
!                to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF (M < 0) THEN
         INFO = -1
      ELSE IF (N < 0) THEN
         INFO = -2
      ELSE IF (LDA < MAX(1,M)) THEN
         INFO = -4
      ENDIF
      IF (INFO .NE. 0) THEN
         PRINT *,'INFO = ', INFO,' IN SGETRF'
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
!
!     Determine the block size for this environment.
!
      NB = ILAENV(1,'SGETRF',' ',M,N,-1,-1)
!
      IF (NB<=1 .OR. NB>=MIN(M,N)) THEN
!
!        Use unblocked code.
!
         CALL GETF2_G (M, N, A, LDA, IPIV, INFO)
      ELSE
!
!        Use blocked code.
!
         DO J = 1, MIN(M,N), NB
            JB = MIN(MIN(M,N) - J + 1,NB)
!
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
!
            CALL GETF2_G (M - J + 1, JB, A(J,J), LDA, IPIV(J), IINFO)
!
!           Adjust INFO and the pivot indices.
!
            IF (INFO.EQ.0 .AND. IINFO>0) INFO = IINFO + J - 1
            IPIV(J:MIN(M,J+JB-1)) = J - 1 + IPIV(J:MIN(M,J+JB-1))
!
!           Apply interchanges to columns 1:J-1.
!
            CALL LASWP_G (J - 1, A, LDA, J, J + JB - 1, IPIV, 1)
!
            IF (J + JB <= N) THEN
!
!              Apply interchanges to columns J+JB:N.
!
               CALL LASWP_G (N-J-JB+1,A(1,J+JB),LDA,J,J+JB-1,IPIV,1)
!
!              Compute block row of U. (S/DTRSM) and update trailing
!              submatrix (S/DGEMM)
!
               CALL TRSM_G ('Left', 'Lower', 'No transpose', 'Unit',
     1                        JB, N - J - JB + 1, ONE, A(J,J),
     2                        LDA, A(J,J+JB), LDA)


               IF (J + JB <= M) CALL GEMM_G ('No transpose',
     1                'No transpose', M - J - JB + 1,
     2                 N - J - JB + 1, JB, (-ONE), A(J+JB,J),
     3                 LDA, A(J,J+JB), LDA, ONE, A(J+JB,J+JB), LDA)

            ENDIF
         END DO
      ENDIF

      END SUBROUTINE GETRF_F90


      SUBROUTINE GETRS_F90(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)
      USE KIND_SPEC
      IMPLICIT NONE
!-----------------------------------------------
!   D U M M Y   A R G U M E N T S
!-----------------------------------------------
      CHARACTER(LEN=1), INTENT(IN) :: TRANS
      INTEGER, INTENT(IN) :: LDA, LDB, NRHS, N
      INTEGER, INTENT(OUT) :: INFO
      INTEGER, INTENT(IN) :: IPIV(N)
      REAL(RPREC), INTENT(IN) :: A(LDA,N)
      REAL(RPREC), INTENT(INOUT) :: B(LDB,NRHS)
!-----------------------------------------------
!   L O C A L   P A R A M E T E R S
!-----------------------------------------------
      REAL(RPREC), PARAMETER :: ONE = 1
!-----------------------------------------------
!   L O C A L   V A R I A B L E S
!-----------------------------------------------
      LOGICAL :: NOTRAN
!-----------------------------------------------
!   E X T E R N A L   F U N C T I O N S
!-----------------------------------------------
      LOGICAL , EXTERNAL :: LSAME
      EXTERNAL LASWP_G, TRSM_G
!-----------------------------------------------
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  GETRS solves a system of linear equations
!     A * X = B  or  A'' * X = B
!  with a general N-by-N matrix A using the LU factorization computed
!  by GETRF.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER*1
!          Specifies the form of the system of equations:
!          = 'N':  A * X = B  (No transpose)
!          = 'T':  A''* X = B  (Transpose)
!          = 'C':  A''* X = B  (Conjugate transpose = Transpose)
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input) REAL array, dimension (LDA,N)
!          The factors L and U from the factorization A = P*L*U
!          as computed by SGETRF.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices from SGETRF; for 1<=i<=N, row i of the
!          matrix was interchanged with row IPIV(i).
!
!  B       (input/output) REAL array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      NOTRAN = LSAME(TRANS,'N')
      IF ( .NOT.NOTRAN .AND.  .NOT.LSAME(TRANS,'T') .AND.
     1     .NOT.LSAME(TRANS,'C')) THEN
         INFO = -1
      ELSE IF (N < 0) THEN
         INFO = -2
      ELSE IF (NRHS < 0) THEN
         INFO = -3
      ELSE IF (LDA < MAX(1,N)) THEN
         INFO = -5
      ELSE IF (LDB < MAX(1,N)) THEN
         INFO = -8
      ENDIF
      IF (INFO .NE. 0) THEN
         PRINT *, 'INFO = ', INFO, 'IN GETRS'
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF (N.EQ.0 .OR. NRHS.EQ.0) RETURN
!

      IF (NOTRAN) THEN
!
!           Solve A * X = B.
!
!           Apply row interchanges to the right hand sides.
            CALL LASWP_G (NRHS, B, LDB, 1, N, IPIV, 1)
!
!           Solve L*X = B, overwriting B with X.
            CALL TRSM_G('Left', 'Lower', 'No transpose', 'Unit',
     1                  N, NRHS, ONE, A, LDA, B, LDB)
!
!           Solve U*X = B, overwriting B with X.
            CALL TRSM_G('Left', 'Upper', 'No transpose', 'Non-unit',
     1                  N, NRHS, ONE, A, LDA, B, LDB)
      ELSE
!
!           Solve A'' * X = B.
!
!           Solve U''*X = B, overwriting B with X.
            CALL TRSM_G('Left', 'Upper', 'Transpose', 'Non-unit', N,
     1                  NRHS, ONE, A, LDA, B, LDB)
!
!           Solve L''*X = B, overwriting B with X.
            CALL TRSM_G('Left', 'Lower', 'Transpose', 'Unit', N, NRHS,
     1                  ONE, A, LDA, B, LDB)
!
!           Apply row interchanges to the solution vectors.
            CALL LASWP_G (NRHS, B, LDB, 1, N, IPIV, -1)
      END IF

      END SUBROUTINE GETRS_F90


      SUBROUTINE GESV_F90 (N, NRHS, A, LDA, IPIV, B, LDB, INFO)
      USE KIND_SPEC
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER, INTENT(IN) :: LDA, LDB, NRHS, N
      INTEGER, INTENT(OUT) :: INFO
      INTEGER, INTENT(OUT) :: IPIV(*)
      REAL(RPREC), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
C-----------------------------------------------
*  Purpose
*  =======
*
*  SGESV computes the solution to a real system of linear equations
*     A * X = B,
*  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
*
*  The LU decomposition with partial pivoting and row interchanges is
*  used to factor A as
*     A = P * L * U,
*  where P is a permutation matrix, L is unit lower triangular, and U is
*  upper triangular.  The factored form of A is then used to solve the
*  system of equations A * X = B.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of linear equations, i.e., the order of the
*          matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input/output) REAL array, dimension (LDA,N)
*          On entry, the N-by-N coefficient matrix A.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (output) INTEGER array, dimension (N)
*          The pivot indices that define the permutation matrix P;
*          row i of the matrix was interchanged with row IPIV(i).
*
*  B       (input/output) REAL array, dimension (LDB,NRHS)
*          On entry, the N-by-NRHS matrix of right hand side matrix B.
*          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
*                has been completed, but the factor U is exactly
*                singular, so the solution could not be computed.
*
*  =====================================================================
*
*     .. External Subroutines ..
      EXTERNAL           GETRF_G, GETRS_G, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGESV ', -INFO )
         RETURN
      END IF
*
*     Compute the LU factorization of A.
*
      CALL GETRF_G( N, N, A, LDA, IPIV, INFO )
      IF( INFO.EQ.0 ) THEN
*
*        Solve the system A*X = B, overwriting B with X.
*
         CALL GETRS_G( 'No transpose', N, NRHS, A, LDA, IPIV, B, LDB,
     $                INFO )
      END IF
      END SUBROUTINE GESV_F90


      SUBROUTINE LASWP_F90(N, A, LDA, K1, K2, IPIV, INCX)
      USE KIND_SPEC
      IMPLICIT NONE
!-----------------------------------------------
!   D U M M Y   A R G U M E N T S
!-----------------------------------------------
      INTEGER N, LDA, K1, K2, INCX
      INTEGER, DIMENSION(*) :: IPIV
      REAL(RPREC), DIMENSION(LDA,*) :: A
!-----------------------------------------------
!   L O C A L   V A R I A B L E S
!-----------------------------------------------
      INTEGER :: I, IP, IX
!-----------------------------------------------
!   E X T E R N A L   F U N C T I O N S
!-----------------------------------------------
      EXTERNAL SWAP_G
!-----------------------------------------------
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  SLASWP performs a series of row interchanges on the matrix A.
!  One row interchange is initiated for each of rows K1 through K2 of A.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.
!
!  A       (input/output) REAL array, dimension (LDA,N)
!          On entry, the matrix of column dimension N to which the row
!          interchanges will be applied.
!          On exit, the permuted matrix.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!
!  K1      (input) INTEGER
!          The first element of IPIV for which a row interchange will
!          be done.
!
!  K2      (input) INTEGER
!          The last element of IPIV for which a row interchange will
!          be done.
!
!  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
!          The vector of pivot indices.  Only the elements in positions
!          K1 through K2 of IPIV are accessed.
!          IPIV(K) = L implies rows K and L are to be interchanged.
!
!  INCX    (input) INTEGER
!          The increment between successive values of IPIV.  If IPIV
!          is negative, the pivots are applied in reverse order.
!
! =====================================================================
!
!     .. Local Scalars ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Executable Statements ..
!
!     Interchange row I with row IPIV(I) for each of rows K1 through K2.
!
      IF (INCX .EQ. 0) RETURN
      IF (INCX > 0) THEN
         IX = K1
      ELSE
         IX = 1 + (1 - K2)*INCX
      ENDIF
      SELECT CASE (INCX)
      CASE (1)
         DO I = K1, K2
            IP = IPIV(I)
            IF (IP .NE. I) CALL SWAP_G (N, A(I,1), LDA, A(IP,1), LDA)
         END DO
      CASE (2:)
         DO I = K1, K2
            IP = IPIV(IX)
            IF (IP .NE. I) CALL SWAP_G (N, A(I,1), LDA, A(IP,1), LDA)
            IX = IX + INCX
         END DO
      CASE (:(-1))
         DO I = K2, K1, -1
            IP = IPIV(IX)
            IF (IP .NE. I) CALL SWAP_G (N, A(I,1), LDA, A(IP,1), LDA)
            IX = IX + INCX
         END DO
      END SELECT

      END SUBROUTINE LASWP_F90


      FUNCTION LAMCH_F90( CMACH )
      USE KIND_SPEC, WP => RPREC
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER*1,  INTENT(IN)  :: CMACH

      REAL(WP) :: LAMCH_F90, RMACH
*     ..
*
*  Purpose
*  =======
*
*  SLAMCH determines single precision machine parameters.
*
*  Arguments
*  =========
*
*  CMACH   (input) CHARACTER*1
*          Specifies the value to be returned by SLAMCH:
*          = 'E' or 'e',   SLAMCH := eps
*          = 'S' or 's',   SLAMCH := sfmin
*          = 'B' or 'b',   SLAMCH := base
*          = 'P' or 'p',   SLAMCH := eps*base
*          = 'N' or 'n',   SLAMCH := t
*          = 'R' or 'r',   SLAMCH := rnd
*          = 'M' or 'm',   SLAMCH := emin
*          = 'U' or 'u',   SLAMCH := rmin
*          = 'L' or 'l',   SLAMCH := emax
*          = 'O' or 'o',   SLAMCH := rmax
*
*          where
*
*          eps   = relative machine precision                           EPSILON/RADIX
*          sfmin = safe minimum, such that 1/sfmin does not overflow    TINY
*          base  = base of the machine                                  RADIX
*          prec  = eps*base                                             EPSILON
*          t     = number of (base) digits in the mantissa              DIGITS
*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*          emin  = minimum exponent before (gradual) underflow          MINEXPONENT
*          rmin  = underflow threshold - base**(emin-1)                 TINY
*          emax  = largest exponent before overflow                     MAXEXPONENT
*          rmax  = overflow threshold  - (base**emax)*(1-eps)           HUGE
*
* =====================================================================
*
*     .. Parameters ..
      REAL(WP), PARAMETER :: ONE = 1.0_WP
*     ..
*     .. External Functions ..
      LOGICAL, EXTERNAL :: LSAME
*     ..
*     .. Executable Statements ..
*
*
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPSILON(ONE)/RADIX(ONE)
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         RMACH = TINY(ONE)
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = RADIX(ONE)
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = EPSILON(ONE)
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = DIGITS(ONE)
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = ONE
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = MINEXPONENT(ONE)
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = TINY(ONE)
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = MAXEXPONENT(ONE)
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = HUGE(ONE)
      END IF
*
      LAMCH_F90 = RMACH

      END FUNCTION LAMCH_F90


      SUBROUTINE GETF2_F90(M, N, A, LDA, IPIV, INFO)
      USE KIND_SPEC
      IMPLICIT NONE
!-----------------------------------------------
!   D U M M Y   A R G U M E N T S
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: M, N, LDA
      INTEGER, INTENT(OUT) :: INFO
      INTEGER, INTENT(OUT) :: IPIV(*)
      REAL(RPREC), INTENT(INOUT) :: A(LDA,*)
!-----------------------------------------------
!   L O C A L   P A R A M E T E R S
!-----------------------------------------------
      REAL(RPREC), PARAMETER :: ONE = 1, ZERO = 0
!-----------------------------------------------
!   L O C A L   V A R I A B L E S
!-----------------------------------------------
      INTEGER :: ISMAX (1)
      INTEGER :: J, JP, I
!-----------------------------------------------
!   E X T E R N A L   F U N C T I O N S
!-----------------------------------------------
      EXTERNAL SWAP_G, GER_G
!-----------------------------------------------
!   I N T R I N S I C  F U N C T I O N S
!-----------------------------------------------
      INTRINSIC MAX, MIN
!-----------------------------------------------
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1992
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  SGETF2 computes an LU factorization of a general m-by-n matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 2 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) REAL array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF (M < 0) THEN
         INFO = -1
      ELSE IF (N < 0) THEN
         INFO = -2
      ELSE IF (LDA < MAX(1,M)) THEN
         INFO = -4
      ENDIF
      IF (INFO .NE. 0) THEN
         PRINT *, 'INFO = ', INFO,' IN SGETF2'
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
!
      DO J = 1, MIN(M,N)
!
!        Find pivot and test for singularity.
!
!        JP = J - 1 + IDAMAX(M - J + 1,A(J,J),1)
         ISMAX = MAXLOC(ABS(A(J:M,J)))
         JP = J - 1 + ISMAX(1)
         IPIV(J) = JP
         IF (A(JP,J) .NE. ZERO) THEN
!
!           Apply the interchange to columns 1:N.
!
            IF (JP .NE. J) CALL SWAP_G (N, A(J,1), LDA, A(JP,1), LDA)
!
!           Compute elements J+1:M of J-th column.
!
!            IF (J < M) CALL DSCAL (M - J, ONE/A(J,J), A(J+1,J), 1)
            A(J+1:M,J) = ONE/A(J,J)*A(J+1:M,J)
!
         ELSE IF (INFO .EQ. 0) THEN
!
            INFO = J
         ENDIF
!
!           IF (J < MIN(M,N)) CALL GER_G (M - J, N - J, (-ONE), A(J+1,J),
!     1          1, A(J,J+1), LDA, A(J+1,J+1), LDA)


!        Update trailing submatrix.
            IF (J < MIN(M,N)) THEN
                DO i=J+1, N
                   A(J+1:M,i) = A(J+1:M,i) - A(J+1:M,J)*A(J,i)
                END DO
            END IF

      END DO

      END SUBROUTINE GETF2_F90


      SUBROUTINE GER_F90(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
      USE KIND_SPEC
      IMPLICIT NONE

      INTEGER :: M,N,INCX,INCY,LDA
      REAL(RPREC) :: ALPHA, X(*), Y(*), A(LDA,*)
*
*  Purpose
*  =======
*
*  SGER   performs the rank 1 operation
*
*     A := alpha*x*y'' + A,
*
*  where alpha is a scalar, x is an m element vector, y is an n element
*  vector and A is an m by n matrix.
*
*  Parameters
*  ==========
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least
*           ( 1 + ( m - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the m
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*  Y      - REAL             array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients. On exit, A is
*           overwritten by the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least max(1,m).
*           Unchanged on exit.
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 30-August-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      INTEGER :: I,IX,J,JY,KX
      REAL(RPREC), PARAMETER :: ZERO = 0
      REAL(RPREC) :: TEMP
      LOGICAL :: OK

      OK = (M.GT.0) .AND. (N.GT.0) .AND. (LDA.GE.M)
*
*
*     Quick return if possible.
*
      IF ( .NOT. OK .OR. (ALPHA.EQ.ZERO)) RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
          DO 20,J = 1,N
             IF (Y(J).NE.ZERO) THEN
                 TEMP = ALPHA*Y(J)
                 DO 10,I = 1,M
                    A(I,J) = A(I,J) + X(I)*TEMP
   10            CONTINUE
             END IF
*
   20     CONTINUE
*
      ELSE
          IF (INCX.GT.0) THEN
              KX = 1
*
          ELSE
              KX = 1 - (M-1)*INCX
          END IF
*
          IF (INCY.GT.0) THEN
              JY = 1
*
          ELSE
              JY = 1 - (N-1)*INCY
          END IF
*
          DO 40,J = 1,N
             IF (Y(JY).NE.ZERO) THEN
                 TEMP = ALPHA*Y(JY)
                 IX = KX
                 DO 30,I = 1,M
                    A(I,J) = A(I,J) + X(IX)*TEMP
                    IX = IX + INCX
   30            CONTINUE
             END IF
*
             JY = JY + INCY
   40     CONTINUE
      END IF
*
      RETURN
*
*     End of SGER  .
*
      END SUBROUTINE GER_F90


      INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
*     ..
*
*  Purpose
*  =======
*
*  ILAENV is called from the LAPACK routines to choose problem-dependent
*  parameters for the local environment.  See ISPEC for a description of
*  the parameters.
*
*  This version provides a set of parameters which should give good,
*  but not optimal, performance on many of the currently available
*  computers.  Users are encouraged to modify this subroutine to set
*  the tuning parameters for their particular machine using the option
*  and problem size information in the arguments.
*
*  This routine will not function correctly if it is converted to all
*  lower case.  Converting it to all upper case is allowed.
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies the parameter to be returned as the value of
*          ILAENV.
*          = 1: the optimal blocksize; if this value is 1, an unblocked
*               algorithm will give the best performance.
*          = 2: the minimum block size for which the block routine
*               should be used; if the usable block size is less than
*               this value, an unblocked routine should be used.
*          = 3: the crossover point (in a block routine, for N less
*               than this value, an unblocked routine should be used)
*          = 4: the number of shifts, used in the nonsymmetric
*               eigenvalue routines
*          = 5: the minimum column dimension for blocking to be used;
*               rectangular blocks must have dimension at least k by m,
*               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
*          = 6: the crossover point for the SVD (when reducing an m by n
*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
*               this value, a QR factorization is used first to reduce
*               the matrix to a triangular form.)
*          = 7: the number of processors
*          = 8: the crossover point for the multishift QR and QZ methods
*               for nonsymmetric eigenvalue problems.
*          = 9: maximum size of the subproblems at the bottom of the
*               computation tree in the divide-and-conquer algorithm
*               (used by xGELSD and xGESDD)
*          =10: ieee NaN arithmetic can be trusted not to trap
*          =11: infinity arithmetic can be trusted not to trap
*
*  NAME    (input) CHARACTER*(*)
*          The name of the calling subroutine, in either upper case or
*          lower case.
*
*  OPTS    (input) CHARACTER*(*)
*          The character options to the subroutine NAME, concatenated
*          into a single character string.  For example, UPLO = 'U',
*          TRANS = 'T', and DIAG = 'N' for a triangular routine would
*          be specified as OPTS = 'UTN'.
*
*  N1      (input) INTEGER
*  N2      (input) INTEGER
*  N3      (input) INTEGER
*  N4      (input) INTEGER
*          Problem dimensions for the subroutine NAME; these may not all
*          be required.
*
* (ILAENV) (output) INTEGER
*          >= 0: the value of the parameter specified by ISPEC
*          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  The following conventions have been used when calling ILAENV from the
*  LAPACK routines:
*  1)  OPTS is a concatenation of all of the character options to
*      subroutine NAME, in the same order that they appear in the
*      argument list for NAME, even if they are not used in determining
*      the value of the parameter specified by ISPEC.
*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
*      that they appear in the argument list for NAME.  N1 is used
*      first, N2 second, and so on, and unused problem dimensions are
*      passed a value of -1.
*  3)  The parameter value returned by ILAENV is checked for validity in
*      the calling subroutine.  For example, ILAENV is used to retrieve
*      the optimal blocksize for STRTRI as follows:
*
*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
*      IF( NB.LE.1 ) NB = MAX( 1, N )
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
*     ..
*     .. External Functions ..
      INTEGER, EXTERNAL :: IEEECK
*     ..
*     .. Executable Statements ..
*
      GO TO ( 100, 100, 100, 400, 500, 600, 700, 800, 900, 1000,
     $        1100 ) ISPEC
*
*     Invalid value for ISPEC
*
      ILAENV = -1
      RETURN
*
  100 CONTINUE
*
*     Convert NAME to upper case if the first character is lower case.
*
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
*
*        ASCII character set
*
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
*
*        EBCDIC character set
*
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $             ( IC.GE.162 .AND. IC.LE.169 ) )
     $            SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
*
*        Prime machines:  ASCII+128
*
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF
*
      C1 = SUBNAM( 1:1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
*
      GO TO ( 110, 200, 300 ) ISPEC
*
  110 CONTINUE
*
*     ISPEC = 1:  block size
*
*     In these examples, separate code is provided for setting NB for
*     real and complex.  We assume that NB will take the same value in
*     single or double precision.
*
      NB = 1
*
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
*
  200 CONTINUE
*
*     ISPEC = 2:  minimum block size
*
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
*
  300 CONTINUE
*
*     ISPEC = 3:  crossover point
*
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
*
  400 CONTINUE
*
*     ISPEC = 4:  number of shifts (used by xHSEQR)
*
      ILAENV = 6
      RETURN
*
  500 CONTINUE
*
*     ISPEC = 5:  minimum column dimension (not used)
*
      ILAENV = 2
      RETURN
*
  600 CONTINUE
*
*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
*
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
*
  700 CONTINUE
*
*     ISPEC = 7:  number of processors (not used)
*
      ILAENV = 1
      RETURN
*
  800 CONTINUE
*
*     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
*
      ILAENV = 50
      RETURN
*
  900 CONTINUE
*
*     ISPEC = 9:  maximum size of the subproblems at the bottom of the
*                 computation tree in the divide-and-conquer algorithm
*                 (used by xGELSD and xGESDD)
*
      ILAENV = 25
      RETURN
*
 1000 CONTINUE
*
*     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
*
C     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 0, 0.0, 1.0 )
      END IF
      RETURN
*
 1100 CONTINUE
*
*     ISPEC = 11: infinity arithmetic can be trusted not to trap
*
C     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 1, 0.0, 1.0 )
      END IF
*
      END FUNCTION ILAENV


      INTEGER FUNCTION IEEECK( ISPEC, ZERO, ONE )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1998
*
*     .. Scalar Arguments ..
      INTEGER            ISPEC
      REAL               ONE, ZERO
*     ..
*
*  Purpose
*  =======
*
*  IEEECK is called from the ILAENV to verify that Infinity and
*  possibly NaN arithmetic is safe (i.e. will not trap).
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies whether to test just for inifinity arithmetic
*          or whether to test for infinity and NaN arithmetic.
*          = 0: Verify infinity arithmetic only.
*          = 1: Verify infinity and NaN arithmetic.
*
*  ZERO    (input) REAL
*          Must contain the value 0.0
*          This is passed to prevent the compiler from optimizing
*          away this code.
*
*  ONE     (input) REAL
*          Must contain the value 1.0
*          This is passed to prevent the compiler from optimizing
*          away this code.
*
*  RETURN VALUE:  INTEGER
*          = 0:  Arithmetic failed to produce the correct answers
*          = 1:  Arithmetic produced the correct answers
*
*     .. Local Scalars ..
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF,
     $                   NEGZRO, NEWZRO, POSINF
*     ..
*     .. Executable Statements ..
      IEEECK = 1
*
      POSINF = ONE / ZERO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = -ONE / ZERO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGZRO = ONE / ( NEGINF+ONE )
      IF( NEGZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = ONE / NEGZRO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEWZRO = NEGZRO + ZERO
      IF( NEWZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      POSINF = ONE / NEWZRO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = NEGINF*POSINF
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      POSINF = POSINF*POSINF
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
*
*
*
*     Return if we were only asked to check infinity arithmetic
*
      IF( ISPEC.EQ.0 )
     $   RETURN
*
      NAN1 = POSINF + NEGINF
*
      NAN2 = POSINF / NEGINF
*
      NAN3 = POSINF / POSINF
*
      NAN4 = POSINF*ZERO
*
      NAN5 = NEGINF*NEGZRO
*
      NAN6 = NAN5*0.0
*
      IF( NAN1.EQ.NAN1 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN2.EQ.NAN2 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN3.EQ.NAN3 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN4.EQ.NAN4 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN5.EQ.NAN5 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN6.EQ.NAN6 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      END FUNCTION IEEECK



      SUBROUTINE SCAL_F90 (N, SA, SX, INCX)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     GENERIC F90 VERSIONS (THESE ROUTINES TAKE ANY RPREC PARAMETER)
!     MOST OF THESE ARE IN BLAS, SO WE PROBABLY DO NOT NEED THEM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      USE KIND_SPEC
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER, INTENT(IN) :: N, INCX
      REAL(RPREC), INTENT(IN) :: SA
      REAL(RPREC), INTENT(INOUT) :: SX(*)
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      INTEGER :: I, IX
C-----------------------------------------------
c
c     scales a vector by a constant.
c     jack dongarra, linpack, 3/11/78.
c     modified to correct problem with negative increments, 9/29/88.
c
      if(n.le.0)return
      if (incx .ne. 1) then
c
c     code for unequal increments or equal increments
c     not equal to 1
c
         ix = 1
         if(incx.lt.0)ix = (-n+1)*incx + 1
         sx(ix:(n-1)*incx+ix:incx) = sa*sx(ix:(n-1)*incx+ix:incx)
      else
c
c     code for both increments equal to 1
c
         sx(:n) = sa * sx(:n)
      endif

      END SUBROUTINE SCAL_F90


      SUBROUTINE AXPY_F90 (N, SA, SX, INCX, SY, INCY)
      USE KIND_SPEC, WP => RPREC
      IMPLICIT NONE
C
C     CONSTANT TIMES A VECTOR PLUS A VECTOR.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER, INTENT(IN) :: N, INCX, INCY
      REAL(WP), INTENT(IN) :: SA, SX(*)
      REAL(WP), INTENT(INOUT) :: SY(*)
C-----------------------------------------------
C   L O C A L   P A R A M E T E R S
C-----------------------------------------------
      REAL(WP), PARAMETER :: ZERO = 0.0_WP
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      INTEGER :: IX, IY
C-----------------------------------------------
C
      if (n .le. 0) return
      if (sa .eq. zero) return

      if (incx.ne.1 .or. incy.ne.1) then
c
c     code for unequal increments or equal increments
c     not equal to 1
c
         ix = 1
         iy = 1
         if (incx .lt. 0)ix = (-n+1)*incx + 1
         if (incy .lt. 0)iy = (-n+1)*incy + 1

         sy(iy:(n-1)*incy+iy:incy) = sy(iy:(n-1)*incy+iy:incy)
     1            + sa * sx(ix:(n-1)*incx+ix:incx)
      else
c
c     code for both increments equal to 1
c
         sy(:n) = sy(:n) + sa*sx(:n)
      end if

      END SUBROUTINE AXPY_F90


      FUNCTION DOT_F90 (N, SX, INCX, SY, INCY)
      USE KIND_SPEC
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER :: N, INCX, INCY
      REAL(RPREC) :: SX(*), SY(*)
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      INTEGER :: IX, IY
      REAL(RPREC) :: DOT_F90
C-----------------------------------------------
c
c     forms the dot product of two vectors.
c     jack dongarra, linpack, 3/11/78.
c
c
      if (n .le. 0) return
      if (incx.ne.1 .or. incy.ne.1) then
c
c     code for unequal increments or equal increments
c     not equal to 1
c
         ix = 1
         iy = 1
         if (incx .lt. 0) ix = (-n+1)*incx + 1
         if (incy .lt. 0) iy = (-n+1)*incy + 1
         DOT_F90 = sum (sx(ix:(n-1)*incx+ix:incx)*
     1                   sy(iy:(n-1)*incy+iy:incy))
      else
c
c     code for both increments equal to 1
c
         DOT_F90 = sum(sx(:n)*sy(:n))
      end if

      END FUNCTION DOT_F90


      SUBROUTINE COPY_F90 (N, SX, INCX, SY, INCY)
      USE KIND_SPEC
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER N, INCX, INCY
      REAL(RPREC), DIMENSION(*) :: SX, SY
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      INTEGER :: IX, IY
C-----------------------------------------------
c
c     copies a vector, x, to a vector, y.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
c
      if (n .le. 0) return
      if (incx.ne.1 .or. incy.ne.1) then
c
c        code for unequal increments or equal increments
c          not equal to 1
c
         ix = 1
         iy = 1
         if (incx .lt. 0) ix = ((-n) + 1)*incx + 1
         if (incy .lt. 0) iy = ((-n) + 1)*incy + 1
         sy(iy:(n-1)*incy+iy:incy) = sx(ix:(n-1)*incx+ix:incx)
      else
c
c     code for both increments equal to 1
c
         sy(:n) = sx(:n)
      end if

      END SUBROUTINE COPY_F90


      SUBROUTINE SWAP_F90(N, SX, INCX, SY, INCY)
      USE KIND_SPEC
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER :: N, INCX, INCY
      REAL(RPREC), DIMENSION(*) :: SX, SY
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      INTEGER :: IX, IY
      REAL(RPREC), DIMENSION(:), ALLOCATABLE :: STEMP
C-----------------------------------------------
c
c     interchanges two vectors.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
c
      if (n .le. 0) return

      allocate (stemp(n), stat=ix)
      if (ix .ne. 0) stop 'Allocation error in swap'

      if (incx.ne.1 .or. incy.ne.1) then
c
c       code for unequal increments or equal increments not equal
c         to 1
c
         ix = 1
         iy = 1
         if (incx .lt. 0) ix = ((-n) + 1)*incx + 1
         if (incy .lt. 0) iy = ((-n) + 1)*incy + 1
         stemp(:n) = sx(ix:(n-1)*incx+ix:incx)
         sx(ix:(n-1)*incx+ix:incx) = sy(iy:(n-1)*incy+iy:incy)
         sy(iy:(n-1)*incy+iy:incy) = stemp(:n)

      else
c
c     code for both increments equal to 1
c
         stemp(:n) = sx(:n)
         sx(:n) = sy(:n)
         sy(:n) = stemp(:n)
      endif

      deallocate (stemp)

      END SUBROUTINE SWAP_F90


      SUBROUTINE TRSM_F90(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA,
     1    A, LDA, B, LDB)
      USE KIND_SPEC
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER :: M, N, LDA, LDB
      REAL(RPREC) :: ALPHA
      CHARACTER(LEN=1) :: SIDE, UPLO, TRANSA, DIAG
      REAL(RPREC), DIMENSION(LDA,*) :: A
      REAL(RPREC), DIMENSION(LDB,*) :: B
C-----------------------------------------------
C   L O C A L   P A R A M E T E R S
C-----------------------------------------------
      REAL(RPREC), PARAMETER :: ONE = 1, ZERO = 0
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      INTEGER :: I, INFO, J, K, NROWA
      REAL(RPREC) :: TEMP
      LOGICAL :: LSIDE, NOUNIT, UPPER
C-----------------------------------------------
C   E X T E R N A L   F U N C T I O N S
C-----------------------------------------------
      LOGICAL , EXTERNAL :: LSAME
C-----------------------------------------------
C   I N T R I N S I C  F U N C T I O N S
C-----------------------------------------------
      INTRINSIC MAX
C-----------------------------------------------

*     .. ARRAY ARGUMENTS ..
*     ..
*
*  Purpose
*  =======
*
*  STRSM  solves one of the matrix equations
*
*     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
*
*  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*     op( A ) = A   or   op( A ) = A''.
*
*  The matrix X is overwritten on B.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry, SIDE specifies whether op( A ) appears on the left
*           or right of X as follows:
*
*              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
*
*              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n'   op( A ) = A.
*
*              TRANSA = 'T' or 't'   op( A ) = A''.
*
*              TRANSA = 'C' or 'c'   op( A ) = A''.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit triangular
*           as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of B. M must be at
*           least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of B.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            .
*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero then  A is not referenced and  B need not be set before
*           entry.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, k ), where k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*           then LDA must be at least max( 1, n ).
*           Unchanged on exit.
*
*  B      - REAL             array of DIMENSION ( LDB, n ).
*           Before entry,  the leading  m by n part of the array  B must
*           contain  the  right-hand  side  matrix  B,  and  on exit  is
*           overwritten by the solution matrix  X.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
*     .. External Subroutines ..
*     .. Intrinsic Functions ..
*     .. Local Scalars ..
*     .. Parameters ..
      LSIDE = LSAME(SIDE,'L')
      IF (LSIDE) THEN
         NROWA = M
      ELSE
         NROWA = N
      ENDIF
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')
*
      INFO = 0
      IF ( .NOT.LSIDE .AND.  .NOT.LSAME(SIDE,'R')) THEN
         INFO = 1
      ELSE IF ( .NOT.UPPER .AND.  .NOT.LSAME(UPLO,'L')) THEN
         INFO = 2
      ELSE IF ( .NOT.LSAME(TRANSA,'N') .AND.  .NOT.LSAME(TRANSA,'T')
     1       .AND.  .NOT.LSAME(TRANSA,'C')) THEN
         INFO = 3
      ELSE IF ( .NOT.LSAME(DIAG,'U') .AND.  .NOT.LSAME(DIAG,'N')) THEN
         INFO = 4
      ELSE IF (M < 0) THEN
         INFO = 5
      ELSE IF (N < 0) THEN
         INFO = 6
      ELSE IF (LDA < MAX(1,NROWA)) THEN
         INFO = 9
      ELSE IF (LDB < MAX(1,M)) THEN
         INFO = 11
      ENDIF
      IF (INFO .NE. 0) THEN
         PRINT *,'INFO = ', INFO,' IN STRSM '
         RETURN
      ENDIF
*
*     Quick return if possible.
*
      IF (N .EQ. 0) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (ALPHA .EQ. ZERO) THEN
         B(:M,:N) = ZERO
         RETURN
      ENDIF
*
*     Start the operations.
*
      IF (LSIDE) THEN
         IF (LSAME(TRANSA,'N')) THEN
*
*           Form  B := alpha*inv( A )*B.
*
            IF (UPPER) THEN
               DO J = 1, N
                  IF (ALPHA .NE. ONE) THEN
                     B(:M,J) = ALPHA*B(:M,J)
                  ENDIF
                  DO K = M, 1, -1
                     IF (B(K,J) .NE. ZERO) THEN
                        IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                        B(:K-1,J) = B(:K-1,J) - B(K,J)*A(:K-1,K)
                     ENDIF
                  END DO
               END DO
            ELSE
               DO J = 1, N
                  IF (ALPHA .NE. ONE) THEN
                     B(:M,J) = ALPHA*B(:M,J)
                  ENDIF
                  DO K = 1, M
                     IF (B(K,J) .NE. ZERO) THEN
                        IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                        B(K+1:M,J) = B(K+1:M,J) - B(K,J)*A(K+1:M,K)
                     ENDIF
                  END DO
               END DO
            ENDIF
         ELSE
*
*           Form  B := alpha*inv( A'' )*B.
*
            IF (UPPER) THEN
               DO J = 1, N
                  DO I = 1, M
                     TEMP = ALPHA*B(I,J)
                     TEMP = TEMP - SUM(A(:I-1,I)*B(:I-1,J))
                     IF (NOUNIT) TEMP = TEMP/A(I,I)
                     B(I,J) = TEMP
                  END DO
               END DO
            ELSE
               DO J = 1, N
                  DO I = M, 1, -1
                     TEMP = ALPHA*B(I,J)
                     TEMP = TEMP - SUM(A(I+1:M,I)*B(I+1:M,J))
                     IF (NOUNIT) TEMP = TEMP/A(I,I)
                     B(I,J) = TEMP
                  END DO
               END DO
            ENDIF
         ENDIF
      ELSE
         IF (LSAME(TRANSA,'N')) THEN
*
*           Form  B := alpha*B*inv( A ).
*
            IF (UPPER) THEN
               DO J = 1, N
                  IF (ALPHA .NE. ONE) THEN
                     B(:M,J) = ALPHA*B(:M,J)
                  ENDIF
                  DO K = 1, J - 1
                     IF (A(K,J) .NE. ZERO) THEN
                        B(:M,J) = B(:M,J) - A(K,J)*B(:M,K)
                     ENDIF
                  END DO
                  IF (NOUNIT) THEN
                     TEMP = ONE/A(J,J)
                     B(:M,J) = TEMP*B(:M,J)
                  ENDIF
               END DO
            ELSE
               DO J = N, 1, -1
                  IF (ALPHA .NE. ONE) THEN
                     B(:M,J) = ALPHA*B(:M,J)
                  ENDIF
                  DO K = J + 1, N
                     IF (A(K,J) .NE. ZERO) THEN
                        B(:M,J) = B(:M,J) - A(K,J)*B(:M,K)
                     ENDIF
                  END DO
                  IF (NOUNIT) THEN
                     TEMP = ONE/A(J,J)
                     B(:M,J) = TEMP*B(:M,J)
                  ENDIF
               END DO
            ENDIF
         ELSE
*
*           Form  B := alpha*B*inv( A'' ).
*
            IF (UPPER) THEN
               DO K = N, 1, -1
                  IF (NOUNIT) THEN
                     TEMP = ONE/A(K,K)
                     B(:M,K) = TEMP*B(:M,K)
                  ENDIF
                  DO J = 1, K - 1
                     IF (A(J,K) .NE. ZERO) THEN
                        TEMP = A(J,K)
                        B(:M,J) = B(:M,J) - TEMP*B(:M,K)
                     ENDIF
                  END DO
                  IF (ALPHA .NE. ONE) THEN
                     B(:M,K) = ALPHA*B(:M,K)
                  ENDIF
               END DO
            ELSE
               DO K = 1, N
                  IF (NOUNIT) THEN
                     TEMP = ONE/A(K,K)
                     B(:M,K) = TEMP*B(:M,K)
                  ENDIF
                  DO J = K + 1, N
                     IF (A(J,K) .NE. ZERO) THEN
                        TEMP = A(J,K)
                        B(:M,J) = B(:M,J) - TEMP*B(:M,K)
                     ENDIF
                  END DO
                  IF (ALPHA .NE. ONE) THEN
                     B(:M,K) = ALPHA*B(:M,K)
                  ENDIF
               END DO
            ENDIF
         ENDIF
      ENDIF

      END SUBROUTINE TRSM_F90


      SUBROUTINE GEMM_F90(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B,
     1   LDB, BETA, C, LDC)
      USE KIND_SPEC
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: M, N, K, LDA, LDB, LDC
      REAL(RPREC) :: ALPHA, BETA
      CHARACTER(LEN=1) :: TRANSA, TRANSB
      REAL(RPREC), DIMENSION(LDA,*) :: A
      REAL(RPREC), DIMENSION(LDB,*) :: B
      REAL(RPREC), DIMENSION(LDC,*) :: C
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(RPREC), PARAMETER :: ONE = 1, ZERO = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I, INFO, J, L, NCOLA, NROWA, NROWB
      REAL(RPREC) :: TEMP
      LOGICAL :: NOTA, NOTB
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      LOGICAL , EXTERNAL :: LSAME
C-----------------------------------------------
*
*  Purpose
*  =======
*
*  SGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X'',
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A''.
*
*              TRANSA = 'C' or 'c',  op( A ) = A''.
*
*           Unchanged on exit.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B''.
*
*              TRANSB = 'C' or 'c',  op( B ) = B''.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - REAL             array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - REAL            .
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - REAL             array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
*     .. External Subroutines ..
*     .. Intrinsic Functions ..
*     .. Local Scalars ..
*     .. Parameters ..
      NOTA = LSAME(TRANSA,'N')
      NOTB = LSAME(TRANSB,'N')
      IF (NOTA) THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      ENDIF
      IF (NOTB) THEN
         NROWB = K
      ELSE
         NROWB = N
      ENDIF
*
*     Test the input parameters.
*
      INFO = 0
      IF ( .NOT.NOTA .AND.  .NOT.LSAME(TRANSA,'C') .AND.  .NOT.LSAME(
     1   TRANSA,'T')) THEN
         INFO = 1
      ELSE IF ( .NOT.NOTB .AND.  .NOT.LSAME(TRANSB,'C') .AND.  .NOT.
     1      LSAME(TRANSB,'T')) THEN
         INFO = 2
      ELSE IF (M < 0) THEN
         INFO = 3
      ELSE IF (N < 0) THEN
         INFO = 4
      ELSE IF (K < 0) THEN
         INFO = 5
      ELSE IF (LDA < MAX(1,NROWA)) THEN
         INFO = 8
      ELSE IF (LDB < MAX(1,NROWB)) THEN
         INFO = 10
      ELSE IF (LDC < MAX(1,M)) THEN
         INFO = 13
      ENDIF
      IF (INFO .NE. 0) THEN
         PRINT *,'INFO = ', INFO, ' IN GEMM '
         RETURN
      ENDIF
*
*     Quick return if possible.
*
      IF(M.EQ.0.OR.N.EQ.0.OR.(ALPHA.EQ.ZERO.OR.K.EQ.0).AND.BETA.EQ.ONE)
     1   RETURN
*
*     And if  alpha.eq.zero.
*
      IF (ALPHA .EQ. ZERO) THEN
         IF (BETA .EQ. ZERO) THEN
            C(:M,:N) = ZERO
         ELSE
            C(:M,:N) = BETA*C(:M,:N)
         ENDIF
         RETURN
      ENDIF
*
*     Start the operations.
*
      IF (NOTB) THEN
         IF (NOTA) THEN
*
*           Form  C := alpha*A*B + beta*C.
*
            DO J = 1, N
               IF (BETA .EQ. ZERO) THEN
                  C(:M,J) = ZERO
               ELSE IF (BETA .NE. ONE) THEN
                  C(:M,J) = BETA*C(:M,J)
               ENDIF
               DO L = 1, K
                  IF (B(L,J) .NE. ZERO) THEN
                     TEMP = ALPHA*B(L,J)
                     C(:M,J) = C(:M,J) + TEMP*A(:M,L)
                  ENDIF
               END DO
            END DO
         ELSE
*
*           Form  C := alpha*A''*B + beta*C
*
            DO J = 1, N
               DO I = 1, M
                  TEMP = SUM(A(:K,I)*B(:K,J))
                  IF (BETA .EQ. ZERO) THEN
                     C(I,J) = ALPHA*TEMP
                  ELSE
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  ENDIF
               END DO
            END DO
         ENDIF
      ELSE
         IF (NOTA) THEN
*
*           Form  C := alpha*A*B'' + beta*C
*
            DO J = 1, N
               IF (BETA .EQ. ZERO) THEN
                  C(:M,J) = ZERO
               ELSE IF (BETA .NE. ONE) THEN
                  C(:M,J) = BETA*C(:M,J)
               ENDIF
               DO L = 1, K
                  IF (B(J,L) .NE. ZERO) THEN
                     TEMP = ALPHA*B(J,L)
                     C(:M,J) = C(:M,J) + TEMP*A(:M,L)
                  ENDIF
               END DO
            END DO
         ELSE
*
*           Form  C := alpha*A''*B'' + beta*C
*
            DO J = 1, N
               DO I = 1, M
                  TEMP = SUM(A(:K,I)*B(J,:K))
                  IF (BETA .EQ. ZERO) THEN
                     C(I,J) = ALPHA*TEMP
                  ELSE
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  ENDIF
               END DO
            END DO
         ENDIF
      ENDIF
*
      END SUBROUTINE GEMM_F90


      LOGICAL FUNCTION LSAME( CA, CB )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          CA, CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME returns .TRUE. if CA is the same letter as CB regardless of
*  case.
*
*  Arguments
*  =========
*
*  CA      (input) CHARACTER*1
*  CB      (input) CHARACTER*1
*          CA and CB specify the single characters to be compared.
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
*     ..
*     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
*     ..
*     .. Executable Statements ..
*
*     Test if the characters are equal
*
      LSAME = CA.EQ.CB
      IF( LSAME ) RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      ZCODE = ICHAR( 'Z' )
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
*
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
*
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     1       INTA.GE.145 .AND. INTA.LE.153 .OR.
     2       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     1       INTB.GE.145 .AND. INTB.LE.153 .OR.
     2       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
*
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
*
      END FUNCTION LSAME

      SUBROUTINE XERBLA( SRNAME, INFO )
*
*  -- LAPACK auxiliary routine (preliminary version) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
*     ..
*
*  Purpose
*  =======
*
*  XERBLA  is an error handler for the LAPACK routines.
*  It is called by an LAPACK routine if an input parameter has an
*  invalid value.  A message is printed and execution stops.
*
*  Installers may consider modifying the STOP statement in order to
*  call system-specific exception-handling facilities.
*
*  Arguments
*  =========
*
*  SRNAME  (input) CHARACTER*6
*          The name of the routine which called XERBLA.
*
*  INFO    (input) INTEGER
*          The position of the invalid parameter in the parameter list
*          of the calling routine.
*
*
      WRITE( *, FMT = 9999 )SRNAME, INFO
*
      STOP
*
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ',
     $      'an illegal value' )
*
*     End of XERBLA
*
      END SUBROUTINE XERBLA


