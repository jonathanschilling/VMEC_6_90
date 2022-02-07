!
! ===================================
! NIST Guide to Available Math Software.
! Fullsource for module TVAL from package NAPACK.
! Retrieved from NETLIB on Fri Aug 7 17:41:59 1998.
! ===================================
!
      SUBROUTINE TVAL(E, K, L, D, U, N, W) 

!      ________________________________________________________
!     |                                                        |
!     |   FIND THE K-TH SMALLEST EIGENVALUE OF A TRIDIAGONAL   |
!     |  MATRIX WHOSE CROSS-DIAGONAL PRODUCTS ARE NONNEGATIVE  |
!     |USING BOTH THE BISECTION METHOD AND A NEWTON-LIKE METHOD|
!     |                                                        |
!     |    INPUT:                                              |
!     |                                                        |
!     |         K     --INDEX OF DESIRED EIGENVALUE (REPLACE K |
!     |                 BY -K TO GET K-TH LARGEST EIGENVALUE)  |
!     |                                                        |
!     |         L     --SUBDIAGONAL (CAN BE IDENTIFIED WITH U) |
!     |                                                        |
!     |         D     --DIAGONAL                               |
!     |                                                        |
!     |         U     --SUPERDIAGONAL                          |
!     |                                                        |
!     |         N     --MATRIX DIMENSION                       |
!     |                                                        |
!     |         W     --WORK ARRAY (LENGTH AT LEAST N)         |
!     |                                                        |
!     |    OUTPUT:                                             |
!     |                                                        |
!     |         E     --EIGENVALUE                             |
!     |                                                        |
!     |    BUILTIN FUNCTIONS: ABS,MAX                        |
!     |    PACKAGE SUBROUTINES: CP,EQL,INP,STM                 |
!     |________________________________________________________|
!
      use kind_spec
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: K, N 
      REAL(rprec):: E 
      REAL(rprec), DIMENSION(N) :: L, D, U, W 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(rprec), parameter :: zero = 0, one = 1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: I, M 
      REAL(rprec) :: S, T, Y, Z 
!-----------------------------------------------
      IF (N > 1) THEN 
!     -------------------------------------------------------
!     |*** BOUND EIGENVALUES USING GERSCHGORIN'S THEOREM ***|
!     -------------------------------------------------------
         M = N - 1 
         Y = D(1) 
         Z = D(1) 
         S = 0 
         DO I = 1, M 
            W(I) = L(I)*U(I) 
            IF (W(I) < ZERO) GO TO 50 
            S = S + ABS(U(I)) 
            T = D(I) + S 
            Z = MAX(T,Z) 
            T = D(I) - S 
            Y = MIN(T,Y) 
            S = ABS(L(I)) 
         END DO 
         T = D(N) + S 
         Z = MAX(T,Z) 
         T = D(N) - S 
         Y = MIN(T,Y) 
         M = K 
         IF (K .LE. 0) M = N + K + 1 
         T = ONE 
         T = T/2 
         S = ONE + T 
 1003    CONTINUE 
         IF (S .LE. ONE) GO TO 1002 
         T = T/2 
         S = ONE + T 
         IF (S .LE. ONE) GO TO 1002 
         T = T/2 
         S = ONE + T 
         IF (S .LE. ONE) GO TO 1002 
         T = T/2 
         S = ONE + T 
         IF (S .LE. ONE) GO TO 1002 
         T = T/2 
         S = ONE + T 
         GO TO 1003 
 1002    CONTINUE 
         T = 8*N*T*MAX(ABS(Y),ABS(Z)) 
         IF (T .NE. ZERO) THEN 
            CALL STM (E, Y, Z, M, T, D, W, N) 
            RETURN  
         ENDIF 
      ENDIF 
      E = D(1) 
      RETURN  

   50 CONTINUE 
      WRITE (6, *) 'ERROR: SUBROUTINE TVAL CAN ONLY BE APPLIED' 
      WRITE (6, *) 'WHEN THE CROSS-DIAGONAL PRODUCTS ARE NONNEGATIVE' 
      STOP  

      END SUBROUTINE TVAL 


      SUBROUTINE STM(E, Y, Z, K, T, D, P, N) 
!      ________________________________________________________
!     |                                                        |
!     |   FIND THE K-TH SMALLEST EIGENVALUE OF A TRIDIAGONAL   |
!     |  MATRIX WHOSE CROSS-DIAGONAL PRODUCTS ARE NONNEGATIVE  |
!     |USING BOTH THE BISECTION METHOD AND A NEWTON-LIKE METHOD|
!     |                                                        |
!     |    INPUT:                                              |
!     |                                                        |
!     |         Y,Z   --ENDS OF AN INTERVAL THAT CONTAINS THE  |
!     |                 DESIRED EIGENVALUE                     |
!     |                                                        |
!     |         K     --INDEX OF DESIRED EIGENVALUE            |
!     |                                                        |
!     |         T     --TOLERANCE (ITERATIONS CONTINUE UNTIL   |
!     |                 THE ERROR IN THE EIGENVALUE .LE. T)    |
!     |                                                        |
!     |         D     --DIAGONAL OF THE COEFFICIENT MATRIX A   |
!     |                                                        |
!     |         P     --CROSS-DIAGONAL PRODUCTS A SUB I+1,I    |
!     |                 TIMES A SUB I,I+1                      |
!     |                                                        |
!     |         N     --MATRIX DIMENSION                       |
!     |                                                        |
!     |    OUTPUT:                                             |
!     |                                                        |
!     |         E     --EIGENVALUE                             |
!     |                                                        |
!     |    BUILTIN FUNCTIONS: ABS,MAX,SIGN,SQRT              |
!     |    PACKAGE SUBROUTINES: CP,EQL,INP                     |
!     |________________________________________________________|
!
      use kind_spec
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer:: K, N 
      REAL(rprec):: E, Y, Z, T 
      REAL(rprec), DIMENSION(N) :: D, P 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(rprec), parameter :: one = 1, zero = 0
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: I, J, M 
      REAL(rprec) :: A,B,FL,FR,L,Q,R,S,U,V,
     1   W,D2,D3,D4,P2,P3,P4,D34,D42,D23 
      REAL(rprec), DIMENSION(4) :: F, H, X 
      REAL(rprec), DIMENSION(4) :: G 
      REAL(rprec) :: C, GL, GR 
!-----------------------------------------------
      L = Y
      R = Z
      IF ( L .LT. R ) GOTO 10
      L = Z
      R = Y
10    IF ( N .EQ. 1 ) GOTO 210
      S = T + T
      U = ONE
20    U = U/2
      A = ONE + U
      IF (A .GT. ONE) GOTO 20
      U = 5*U
      V = U/2
      CALL CP(L,FL,GL,I,D,P,N)
      CALL CP(R,FR,GR,J,D,P,N)
      IF ( I .GE. K ) GOTO 190
      IF ( J .LT. K ) GOTO 200
C     --------------------------------
C     |*** ISOLATE THE EIGENVALUE ***|
C     --------------------------------
30    E = R - L
      IF ( E .LE. S ) GOTO 180
      IF ( E .LE. U*(ABS(L)+ABS(R)) ) GOTO 180
      IF ( J .EQ. I+1 ) GOTO 70
40    A = .5*(L+R)
      CALL CP(A,B,C,M,D,P,N)
      IF ( K .LE. M ) GOTO 50
      L = A
      I = M
      FL = B
      GL = C
      GOTO 30
50    R = A
      J = M
      FR = B
      GR = C
      GOTO 30
60    E = R - L
      IF ( E .LE. S ) GOTO 180
      IF ( E .LE. U*(ABS(L)+ABS(R)) ) GOTO 180
70    X(1) = L
      F(1) = FL
      G(1) = GL
      X(2) = R
      F(2) = FR
      G(2) = GR
      CALL EQL(X,F,G,H,2)
      IF ( H(1) .EQ. H(2) ) GOTO 160
C     ---------------------
C     |*** SECANT STEP ***|
C     ---------------------
      A = X(1) - H(1)*(X(1)-X(2))/(H(1)-H(2))
      Q = A
      W = max(T,V*(ABS(L)+ABS(R)))
      IF ( ABS(A-L) .LT. W ) A = L + W
      IF ( ABS(A-R) .LT. W ) A = R - W
      CALL CP(A,B,C,J,D,P,N)
      IF ( I .GE. J ) GOTO 80
      R = A
      FR = B
      GR = C
      GOTO 90
80    L = A
      FL = B
      GL = C
90    X(3) = A
      F(3) = B
      G(3) = C
      W = R - L
      IF ( W .LE. S ) GOTO 220
      IF ( W .LE. U*(ABS(L)+ABS(R)) ) GOTO 220
      CALL EQL(X,F,G,H,3)
C     --------------------------------------
C     |*** QUADRATIC INTERPOLATION STEP ***|
C     --------------------------------------
      CALL INP(A,X(1),X(2),X(3),H(1),H(2),H(3),L,R)
      B = L
      IF ( ABS(A-L) .GT. ABS(A-R) ) B = R
C     ------------------------------------
C     |*** APPLY PSEUDO-NEWTON METHOD ***|
C     ------------------------------------
100   Q = A
      W = max(T,V*(ABS(L)+ABS(R)))
      IF ( ABS(A-L) .LT. W ) GOTO 110
      IF ( ABS(A-R) .GT. W ) GOTO 130
110   IF ( A+A .GT. L+R ) GOTO 120
      A = L + W
      GOTO 130
120   A = R - W
130   IF ( A .LE. L ) GOTO 160
      IF ( A .GE. R ) GOTO 160
      E = .5*E
      IF ( E .LT. ABS(B-A) ) GOTO 160
      CALL CP(A,B,C,J,D,P,N)
      IF ( I .GE. J ) GOTO 140
      R = A
      FR = B
      GR = C
      GOTO 150
140   L = A
      FL = B
      GL = C
150   W = R - L
      IF ( W .LE. S ) GOTO 220
      IF ( W .LE. U*(ABS(L)+ABS(R)) ) GOTO 220
      X(4) = A
      F(4) = B
      G(4) = C
      CALL EQL(X,F,G,H,4)
      IF ( X(1) .LT. L ) GOTO 160
      IF ( X(1) .GT. R ) GOTO 160
      B = X(1)
      D4 = X(4) - B
      D3 = X(3) - B
      D2 = X(2) - B
      D34 = X(3) - X(4)
      D42 = X(4) - X(2)
      D23 = X(2) - X(3)
      P2 = D2*(ONE+((D2/D3)*D42+(D2/D4)*D23)/D34)
      P3 = D3*(ONE+((D3/D2)*D34+(D3/D4)*D23)/D42)
      P4 = D4*(ONE+((D4/D2)*D34+(D4/D3)*D42)/D23)
      IF (P2 .NE. ZERO) P2 = (H(2)-H(1))/P2
      IF (P3 .NE. ZERO) P3 = (H(3)-H(1))/P3
      IF (P4 .NE. ZERO) P4 = (H(4)-H(1))/P4
      P2 = P2 + P3 + P4
      IF (P2 .EQ. ZERO) GOTO 160
      A = B - H(1)/P2
      GOTO 100
C     --------------------------
C     |*** BISECTION METHOD ***|
C     --------------------------
160   A = (L+R)/2
      CALL CP(A,B,C,J,D,P,N)
      IF ( I .GE. J ) GOTO 170
      R = A
      FR = B
      GR = C
      GOTO 60
170   L = A
      FL = B
      GL = C
      GOTO 60
180   E = (L+R)/2
      RETURN
190   E = L
      RETURN
200   E = R
      RETURN
210   E = D(1)
      IF ( L .GT. E ) E = L
      IF ( R .LT. E ) E = R
      RETURN
220   E = Q
      RETURN
      END SUBROUTINE STM


      SUBROUTINE INP(A, X, Y, Z, U, V, W, L, R) 
      use kind_spec
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec) :: A, X, Y, Z, U, V, W, L, R 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(rprec), parameter :: zero = 0, one = 1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec) :: B, C, P, Q, S, T 
!-----------------------------------------------
      S = Z - X 
      T = (Y - X)/S 
      A = (U - V)/T + (W - V)/(1. - T) 
      IF (A .NE. ZERO) THEN 
         B = .5_dp*(W - U)/A - .5_dp 
         C = U/A 
         T = SQRT(ABS(C)) 
         IF (ABS(B) < SIGN(T,C)) GO TO 60 
         T = MAX(T,ABS(B)) 
         IF (T .EQ. ZERO) GO TO 50 
         Q = ONE/T 
         P = SQRT((Q*B)**2 - Q*C*Q) 
         P = T*P 
         IF (ABS(P + B) .LE. ABS(P - B)) THEN 
            Q = P - B 
         ELSE 
            Q = -(B + P) 
         ENDIF 
         P = C/Q 
         Q = X + S*Q 
         P = X + S*P 
         IF (Q .GE. L) THEN 
            IF (Q .LE. R) THEN 
               A = Q 
               RETURN  
            ENDIF 
         ENDIF 
         A = P 
         RETURN  
      ENDIF 
      IF (U .NE. W) THEN 
         A = X + S*U/(U - W) 
         RETURN  
      ENDIF 
   50 CONTINUE 
      A = L 
      RETURN  
   60 CONTINUE 
      A = X - S*B 
      RETURN  

      END SUBROUTINE INP 


      SUBROUTINE CP(T, B, V, L, D, P, N) 
!------------------------------------------------
!*** EVALUATE CHARACTERISTIC POLYNOMIAL AND ***
!              SIGN ALTERNATIONS               
!------------------------------------------------
      use kind_spec
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer :: L, N 
      REAL(rprec) ::  T, B 
      REAL(rprec) :: V 
      REAL(rprec), DIMENSION(N) :: D, P 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(rprec), parameter :: zero = 0, one = 1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: I, J, M 
      REAL(rprec) :: A, C, E, F, U, Z 
!-----------------------------------------------
      L = 0 
      Z = 0 
      V = Z 
      IF (N .LE. 1) THEN 
         B = D(1) - T 
         IF (B .LE. Z) L = 1 
         RETURN  
      ENDIF 
      F = 65536
      F = F**4
      E = ONE/F 
      U = 16 
      M = 0 
      I = 1 
   20 CONTINUE 
      C = ONE 
      B = D(I) - T 
   30 CONTINUE 
      IF (B .LE. Z) GO TO 70 
      IF (I .GE. N) GO TO 130 
   40 CONTINUE 
      J = I 
      I = I + 1 
      A = (D(I)-T)*B - P(J)*C 
      C = B 
      B = A 
   50 CONTINUE 
      A = ABS(B) 
      IF (A > F) GO TO 60 
      IF (A > E) GO TO 30 
      IF (A .EQ. Z) GO TO 70 
      C = C*F 
      B = B*F 
      V = V - U 
      GO TO 50 
   60 CONTINUE 
      C = C*E 
      B = B*E 
      V = V + U 
      GO TO 50 
   70 CONTINUE 
      L = L + 1 
      IF (I .GE. N) GO TO 130 
      IF (B < Z) GO TO 90 
      IF (P(I) > ZERO) GO TO 90 
      I = I + 1 
      M = 1 
      V = Z 
      GO TO 20 
   80 CONTINUE 
      IF (B .GE. Z) GO TO 120 
      IF (I .GE. N) GO TO 130 
   90 CONTINUE 
      J = I 
      I = I + 1 
      A = (D(I)-T)*B - P(J)*C 
      C = B 
      B = A 
  100 CONTINUE 
      A = ABS(B) 
      IF (A > F) GO TO 110 
      IF (A > E) GO TO 80 
      IF (A .EQ. Z) GO TO 120 
      C = C*F 
      B = B*F 
      V = V - U 
      GO TO 100 
  110 CONTINUE 
      C = C*E 
      B = B*E 
      V = V + U 
      GO TO 100 
  120 CONTINUE 
      L = L + 1 
      IF (I .GE. N) GO TO 130 
      IF (B > Z) GO TO 40 
      IF (P(I) > ZERO) GO TO 40 
      I = I + 1 
      M = 1 
      V = Z 
      GO TO 20 
  130 CONTINUE 
      IF (M .NE. 1) THEN 
         IF (B .NE. ZERO) THEN 
            A = ONE/U 
            IF (ABS(B) .GE. A) THEN 
  140          CONTINUE 
               IF (ABS(B) < ONE) RETURN  
               B = B*A 
               V = V + ONE 
               GO TO 140 
            ENDIF 
            B = B*U 
            V = V - ONE
 1003       CONTINUE 
            IF (ABS(B) .GE. A) GO TO 1002 
            B = B*U 
            V = V - ONE
            IF (ABS(B) .GE. A) GO TO 1002 
            B = B*U 
            V = V - ONE
            IF (ABS(B) .GE. A) GO TO 1002 
            B = B*U 
            V = V - ONE
            IF (ABS(B) .GE. A) GO TO 1002 
            B = B*U 
            V = V - ONE
            GO TO 1003 
 1002       CONTINUE 
            RETURN  
         ENDIF 
      ENDIF 
      V = Z 
      B = Z 
      RETURN  

      END SUBROUTINE CP 


      SUBROUTINE EQL(X, F, G, H, N) 
      use kind_spec
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer N 
      REAL(rprec), DIMENSION(N) :: X, F, H 
      REAL(rprec), DIMENSION(N) :: G 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: I1, J2, I, J, K, L, M, NM1 
      REAL(rprec) :: R, S, Z, U 
!-----------------------------------------------
      Z = 0 
      I = 1 
      R = X(N) 
      S = F(N) 
      U = G(N) 
      IF (S .EQ. Z) GO TO 30 
      I1 = I 
      J2 = MAX(N - 1,I1) 
      DO I = I1, J2 
         IF (F(I) .EQ. Z) CYCLE  
         IF (U < G(I)) GO TO 30 
         IF (U > G(I)) CYCLE  
         IF (ABS(F(I)) .GE. ABS(S)) GO TO 30 
      END DO 
      GO TO 50 
   30 CONTINUE 
      M = N + I 
      NM1 = N - 1 
      DO J = I, NM1 
         K = M - J 
         L = K - 1 
         X(K) = X(L) 
         F(K) = F(L) 
         G(K) = G(L) 
      END DO 
      X(I) = R 
      F(I) = S 
      G(I) = U 
   50 CONTINUE 
      U = G(N) 
      WHERE ((G(:N)-U).LE.(-99._dp) .OR. F(:N).EQ.Z) H(:N) = Z 
      WHERE (F(:N).NE.Z .AND. (G(:N)-U)>(-99._dp))
     1    H(:N) = F(:N)*16._dp**(G(:N)-U) 

      RETURN  
      END SUBROUTINE EQL 


      SUBROUTINE TVECT(E, X, L, D, U, N, W) 
!      ________________________________________________________
!     |                                                        |
!     |     COMPUTE EIGENVECTOR CORRESPONDING TO GIVEN REAL    |
!     |        EIGENVALUE FOR A REAL TRIDIAGONAL MATRIX        |
!     |                                                        |
!     |    INPUT:                                              |
!     |                                                        |
!     |         E     --EIGENVALUE                             |
!     |                                                        |
!     |         L     --SUBDIAGONAL                            |
!     |                                                        |
!     |         D     --DIAGONAL                               |
!     |                                                        |
!     |         U     --SUPERDIAGONAL                          |
!     |                                                        |
!     |         N     --MATRIX DIMENSION                       |
!     |                                                        |
!     |         W     --WORK ARRAY (LENGTH AT LEAST 4N)        |
!     |                                                        |
!     |    OUTPUT:                                             |
!     |                                                        |
!     |         E     --IMPROVED ESTIMATE FOR EIGENVALUE       |
!     |                                                        |
!     |         X     --EIGENVECTOR                            |
!     |                                                        |
!     |    BUILTIN FUNCTIONS: ABS,SQRT                         |
!     |________________________________________________________|
!
      use kind_spec
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer:: N 
      REAL(rprec):: E 
      REAL(rprec), DIMENSION(N) :: X, L, D, U
      REAL(rprec), DIMENSION(4*N) :: W
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(rprec), parameter :: zero = 0, one = 1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: G1, J1, J2, J4, J3, F, G, H, I, J, K, M 
      REAL(rprec) :: O, P, Q, R, S, T, V, Y, Z 
!-----------------------------------------------
      IF (N .LE. 1) THEN 
         E = D(1) 
         X(1) = ONE 
         RETURN  
      ENDIF 
      M = N - 1 
      J = 2 
!     ---------------------------------------------------------
!     |*** STORE MATRIX IN W ARRAY AND SUBTRACT EIGENVALUE ***|
!     ---------------------------------------------------------
      X(1:M) = 0
      DO I = 1,M
         W(J) = U(I)
         W(J+1) = D(I) - E
         W(J+2) = L(I)
         J = J + 4
      END DO   
      W(J) = 0 
      W(J+1) = D(N) - E 
      O = 65536
      O = O**(-4) 
      T = O/2
      S = T 
      T = T/2
      P = S 
      S = S + T 
      IF (S .GE. O) GO TO 40 
      DO WHILE(S + T > S) 
         T = T/2
         P = S 
         S = S + T 
         IF (S .GE. O) GO TO 40 
      END DO 
   40 CONTINUE 
      R = W(3) 
      V = ABS(R) + ABS(W(4)) 
      F = 4 
      M = 4*N - 3 
      G = -2 
!     ---------------------------
!     |*** FACTOR THE MATRIX ***|
!     ---------------------------
      G = G + 4 
      G1 = G 
      DO G = G1, M, 4 
         H = G - 1 
         I = G + 2 
         J = G + 5 
!     --------------------
!     |*** FIND PIVOT ***|
!     --------------------
         Q = W(I) 
         Y = ABS(Q) 
         Z = ABS(R) 
         IF (Z < Y) THEN 
!     -------------------
!     |*** SWAP ROWS ***|
!     -------------------
            IF (V < Y) THEN 
               V = Y 
               F = I 
            ENDIF 
            T = W(G) 
            W(G) = W(J) 
            W(J) = T 
            T = R/Q 
            K = G + 1 
            W(K) = Q 
            K = J - 1 
            S = W(K) 
            IF (S .NE. ZERO) THEN 
               IF (S .EQ. O) S = P 
               W(K) = -S*T 
               W(H) = S 
               GO TO 100 
            ENDIF 
            W(K) = S 
            W(H) = O 
            GO TO 100 
         ENDIF 
         W(H) = 0 
         IF (V < Z) THEN 
            V = Z 
            F = I 
         ENDIF 
         IF (R .EQ. ZERO) GO TO 120 
         T = Q/R 
!     -------------------
!     |*** ELIMINATE ***|
!     -------------------
  100    CONTINUE 
         R = W(J) - T*W(G) 
         W(J) = R 
         W(I) = T 
      END DO 
      IF (ABS(R) < V) THEN 
         V = R 
         F = J + 1 
!     ---------------------------------------------------
!     |*** COMPUTE INITIAL EIGENVECTOR APPROXIMATION ***|
!     ---------------------------------------------------
      ENDIF 
  120 CONTINUE 
      J = F/4 
      X(J) = ONE 
      IF (J .NE. 1) THEN 
         K = F - 5 
         J = J - 1 
         X(J) = (X(J)-W(K-1)*X(J+1))/W(K) 
         J1 = J 
         IF (J1 - 1 > 0) THEN 
            DO J = J1, 2, -1 
               K = K - 4 
               T = W(K-2) 
               IF (T .EQ. O) T = 0 
               X(J-1) = (X(J-1)-W(K-1)*X(J)-T*X(J+1))/W(K) 
            END DO 
         ENDIF 
         GO TO 140 
      ENDIF 
  140 CONTINUE 
      IF (V .NE. ZERO) THEN 
         S = MAXVAL(ABS(X(:N)))
         S = ONE/S 
         X(:N) = S*X(:N)
         R = SUM(X(:N)*X(:N))
         K = 0 
         J = 1 
         Y = X(1) 
!     -----------------------------------------------------
!     |*** APPLY ONE ITERATION OF INVERSE POWER METHOD ***|
!     -----------------------------------------------------
         J2 = J 
         J4 = MAX(N - 1,J2) 
         DO J = J2, J4 
            K = K + 4 
            I = J 
            S = W(K) 
            W(K) = Y 
            Y = X(J+1) 
            IF (W(K-3) .NE. ZERO) THEN 
               T = X(J+1) 
               X(J+1) = X(I) 
               X(I) = T 
            ENDIF 
            X(J+1) = X(J+1) - S*X(I) 
         END DO 
!     ---------------------------
!     |*** BACK SUBSTITUTION ***|
!     ---------------------------
         S = X(J)/W(K+3) 
         X(J) = S 
         T = ABS(S) 
         V = S*Y 
         J = J - 1 
         K = K - 1 
         S = (X(J)-W(K-1)*S)/W(K) 
         X(J) = S 
         T = MAX(ABS(S),T) 
         V = V + S*W(K+1) 
         J3 = J 
         IF (J3 - 1 > 0) THEN 
            DO J = J3, 2, -1 
               K = K - 4 
               Z = W(K-2) 
               IF (Z .EQ. O) Z = 0 
               S = (X(J-1)-W(K-1)*S-Z*X(J+1))/W(K) 
               X(J-1) = S 
               T = MAX(ABS(S),T) 
               V = V + S*W(K+1) 
            END DO 
         ENDIF 
         IF (V .NE. ZERO) V = R/V 
         T = ONE/T 
         S = SUM((X(:N)*V)**2) 
         Z = SUM((T*X(:N))**2) 
         T = T/SQRT(Z) 
         X(:N) = T*X(:N) 
!     --------------------------------------------------------------
!     |*** USE RAYLEIGH QUOTIENT TO IMPROVE EIGENVALUE ESTIMATE ***|
!     --------------------------------------------------------------
         IF (R + R .GE. S) E = E + V 
         RETURN  
      ENDIF 
      T = MAXVAL(ABS(X(:N)))
      T = ONE/T 
      Z = SUM((T*X(:N))**2) 
      T = T/SQRT(Z) 
      X(:N) = T*X(:N) 

      END SUBROUTINE TVECT 
