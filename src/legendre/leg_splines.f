      SUBROUTINE PCHEZ(N, X, F, D, SPLINE, WK, LWK, IERR)
      USE kind_spec
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: N, LWK, IERR
      LOGICAL SPLINE
      REAL(rprec), DIMENSION(N) :: X, F, D
      REAL(rprec), DIMENSION(LWK) :: WK
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(2) :: IC
      INTEGER :: INCFD
      REAL(rprec), DIMENSION(2) :: VC
!-----------------------------------------------
!***BEGIN PROLOGUE  PCHEZ
!***DATE WRITTEN   870821   (YYMMDD)
!***REVISION DATE  870908   (YYMMDD)
!***CATEGORY NO.  E1B
!***KEYWORDS  CUBIC HERMITE MONOTONE INTERPOLATION, SPLINE
!             INTERPOLATION, EASY TO USE PIECEWISE CUBIC INTERPOLATION
!***AUTHOR  KAHANER, D.K., (NBS)
!             SCIENTIFIC COMPUTING DIVISION
!             NATIONAL BUREAU OF STANDARDS
!             GAITHERSBURG, MARYLAND 20899
!             (301) 975-3808
!***PURPOSE  Easy to use spline or cubic Hermite interpolation.
!***DESCRIPTION
!
!          PCHEZ:  Piecewise Cubic Interpolation, Easy to Use.
!
!     From the book "Numerical Methods and Software"
!          by  D. Kahaner, C. Moler, S. Nash
!               Prentice Hall 1988
!
!     Sets derivatives for spline (two continuous derivatives) or
!     Hermite cubic (one continuous derivative) interpolation.
!     Spline interpolation is smoother, but may not "look" right if the
!     data contains both "steep" and "flat" sections.  Hermite cubics
!     can produce a "visually pleasing" and monotone interpolant to
!     monotone data. This is an easy to use driver for the routines
!     by F. N. Fritsch in reference (4) below. Various boundary
!     conditions are set to default values by PCHEZ. Many other choices
!     are available in the subroutines PCHIC, PCHIM and PCHSP.
!
!     Use PCHEV to evaluate the resulting function and its derivative.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:   CALL  PCHEZ (N, X, F, D, SPLINE, WK, LWK, IERR)
!
!     INTEGER  N, IERR,  LWK
!     REAL  X(N), F(N), D(N), WK(*)
!     LOGICAL SPLINE
!
!   Parameters:
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!           If N=2, simply does linear interpolation.
!
!     X -- (input) real array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real array of dependent variable values to be inter-
!           polated.  F(I) is value corresponding to X(I).
!
!     D -- (output) real array of derivative values at the data points.
!
!     SPLINE -- (input) logical variable to specify if the interpolant
!           is to be a spline with two continuous derivaties
!           (set SPLINE=.TRUE.) or a Hermite cubic interpolant with one
!           continuous derivative (set SPLINE=.FALSE.).
!        Note: If SPLINE=.TRUE. the interpolating spline satisfies the
!           default "not-a-knot" boundary condition, with a continuous
!           third derivative at X(2) and X(N-1). See reference (3).
!              If SPLINE=.FALSE. the interpolating Hermite cubic will be
!           monotone if the input data is monotone. Boundary conditions
!           computed from the derivative of a local quadratic unless this
!           alters monotonicity.
!
!     WK -- (scratch) real work array, which must be declared by the calling
!           program to be at least 2*N if SPLINE is .TRUE. and not used
!           otherwise.
!
!     LWK -- (input) length of work array WK. (Error return if
!           LWK.LT.2*N and SPLINE is .TRUE., not checked otherwise.)
!  
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning error:
!              IERR.GT.0  (can only occur when SPLINE=.FALSE.) means that
!                 IERR switches in the direction of monotonicity were detected.
!                 When SPLINE=.FALSE.,  PCHEZ guarantees that if the input
!                 data is monotone, the interpolant will be too. This warning
!                 is to alert you to the fact that the input data was not
!                 monotone.
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -7  if LWK is less than 2*N and SPLINE is .TRUE.
!             (The D-array has not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!  
!----------------------------------------------------------------------
!***REFERENCES  1. F.N.FRITSCH AND R.E.CARLSON, 'MONOTONE PIECEWISE CUBIC INTERPOLATION,'
!                  SIAM J.NUMER.ANAL. 17, 2 (APRIL 1980), 238-246.
!               2. F.N.FRITSCH AND J.BUTLAND, 'A METHOD FOR CONSTRUCTING LOCAL MONOTONE PIECEWISE CUBIC INTERPOLANTS,'
!                  LLNL PREPRINT UCRL-87559 (APRIL 1982).
!               3. CARL DE BOOR, A PRACTICAL GUIDE TO SPLINES, SPRINGER-
!                 VERLAG (NEW YORK, 1978).  (ESP. CHAPTER IV, PP.49-62.)
!               4. F.N.FRITSCH, 'PIECEWISE CUBIC HERMITE INTERPOLATION PACKAGE, FINAL SPECIFICATIONS',
!                  LAWRENCE LIVERMORE NATIONAL LABORATORY, COMPUTER DOCUMENTATION UCID-30194, AUGUST 1982.
!***ROUTINES CALLED  PCHIM,PCHSP
!***END PROLOGUE  PCHEZ
!
!  DECLARE LOCAL VARIABLES.
!
      DATA IC(1)/0/
      DATA IC(2)/0/
      DATA INCFD/1/
!
!
!***FIRST EXECUTABLE STATEMENT  PCHEZ
!
      IF (SPLINE) THEN
         CALL PCHSP (IC, VC, N, X, F, D, INCFD, WK, LWK, IERR)
      ELSE
         CALL PCHIM (N, X, F, D, INCFD, IERR)
      ENDIF
!
!  ERROR CONDITIONS ALREADY CHECKED IN PCHSP OR PCHIM
 
      RETURN 
!------------- LAST LINE OF PCHEZ FOLLOWS ------------------------------
      END SUBROUTINE PCHEZ


      SUBROUTINE PCHIM(N, X, F, D, INCFD, IERR)
      USE kind_spec
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: N, INCFD, IERR
      REAL(rprec), DIMENSION(N) :: X
      REAL(rprec), DIMENSION(INCFD,N) :: F, D
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, NLESS1
      REAL(rprec) :: R1, DEL1, DEL2, DMAX, DMIN 
      REAL(rprec) :: DRAT1, DRAT2, DSAVE, H1, H2 
      REAL(rprec) :: HSUM, HSUMT3, THREE, W1, W2, ZERO
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      REAL(rprec) , EXTERNAL :: PCHST
!-----------------------------------------------
!***BEGIN PROLOGUE  PCHIM
!***DATE WRITTEN   811103   (YYMMDD)
!***REVISION DATE  870707   (YYMMDD)
!***CATEGORY NO.  E1B
!***KEYWORDS  LIBRARY=SLATEC(PCHIP),
!             TYPE=SINGLE PRECISION(PCHIM-S DPCHIM-D),
!             CUBIC HERMITE INTERPOLATION,MONOTONE INTERPOLATION,
!             PIECEWISE CUBIC INTERPOLATION
!***AUTHOR  FRITSCH, F. N., (LLNL)
!             MATHEMATICS AND STATISTICS DIVISION
!             LAWRENCE LIVERMORE NATIONAL LABORATORY
!             P.O. BOX 808  (L-316)
!             LIVERMORE, CA  94550
!             FTS 532-4275, (415) 422-4275
!***PURPOSE  Set derivatives needed to determine a monotone piecewise
!            cubic Hermite interpolant to given data.  Boundary values
!            are provided which are compatible with monotonicity.  The
!            interpolant will have an extremum at each point where mono-
!            tonicity switches direction.  (See PCHIC if user control is
!            desired over boundary or switch conditions.)
!***DESCRIPTION
!
!          PCHIM:  Piecewise Cubic Hermite Interpolation to
!                  Monotone data.
!
!     Sets derivatives needed to determine a monotone piecewise cubi!
!     Hermite interpolant to the data given in X and F.
!
!     Default boundary conditions are provided which are compatible
!     with monotonicity.  (See PCHIC if user control of boundary con-
!     ditions is desired.)
!
!     If the data are only piecewise monotonic, the interpolant will
!     have an extremum at each point where monotonicity switches direc-
!     tion.  (See PCHIC if user control is desired in such cases.)
!
!     To facilitate two-dimensional applications, includes an increment
!     between successive values of the F- and D-arrays.
!
!     The resulting piecewise cubic Hermite function may be evaluated
!     by PCHFE or PCHFD.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  N, IERR
!        REAL  X(N), F(INCFD,N), D(INCFD,N)
!
!        CALL  PCHIM (N, X, F, D, INCFD, IERR)
!
!   Parameters:
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!           If N=2, simply does linear interpolation.
!
!     X -- (input) real array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real array of dependent variable values to be inter-
!           polated.  F(1+(I-1)*INCFD) is value corresponding to X(I).
!           PCHIM is designed for monotonic data, but it will work for
!           any F-array.  It will force extrema at points where mono-
!           tonicity switches direction.  If some other treatment of
!           switch points is desired, PCHIC should be used instead.
!                                     -----
!     D -- (output) real array of derivative values at the data points.
!           If the data are monotonic, these values will determine a
!           a monotone cubic Hermite function.
!           The value corresponding to X(I) is stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!           No other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in F and D.
!           This argument is provided primarily for 2-D applications.
!           (Error return if  INCFD.LT.1 .)
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning error:
!              IERR.GT.0  means that IERR switches in the direction
!                 of monotonicity were detected.
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -2  if INCFD.LT.1 .
!              IERR = -3  if the X-array is not strictly increasing.
!             (The D-array has not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
!***REFERENCES  1. F.N.FRITSCH AND R.E.CARLSON, 'MONOTONE PIECEWISE CUBIC INTERPOLATION,'
!                 SIAM J.NUMER.ANAL. 17, 2 (APRIL 1980), 238-246.
!              2. F.N.FRITSCH AND J.BUTLAND, 'A METHOD FOR CONSTRUCTING LOCAL MONOTONE PIECEWISE CUBIC INTERPOLANTS,'
!                 LLNL PREPRINT UCRL-87559 (APRIL 1982).
!***ROUTINES CALLED  PCHST
!***END PROLOGUE  PCHIM
!
! ----------------------------------------------------------------------
!
!  Change record:
!     82-02-01   1. Introduced  PCHST  to reduce possible over/under-
!                   flow problems.
!                2. Rearranged derivative formula for same reason.
!     82-06-02   1. Modified end conditions to be continuous functions
!                   of data when monotonicity switches in next interval.
!                2. Modified formulas so end conditions are less prone
!                   of over/underflow problems.
!     82-08-03   Minor cosmetic changes for release 1.
!
! ----------------------------------------------------------------------
!
!  Programming notes:
!
!     1. The function  PCHST(ARG1,ARG2)  is assumed to return zero if
!        either argument is zero, +1 if they are of the same sign, and
!        -1 if they are of opposite sign.
!     2. To produce a double precision version, simply:
!        a. Change PCHIM to DPCHIM wherever it occurs,
!        b. Change PCHST to DPCHST wherever it occurs,
!        c. Change all references to the Fortran intrinsics to their
!           double precision equivalents,
!        d. Change the real declarations to double precision, and
!        e. Change the constants ZERO and THREE to double precision.
!
!  DECLARE ARGUMENTS.
!
!
!  DECLARE LOCAL VARIABLES.
!
      DATA ZERO/0/
      DATA THREE/3/
!
!  VALIDITY-CHECK ARGUMENTS.
!
!***FIRST EXECUTABLE STATEMENT  PCHIM
      IF (N >= 2) THEN
         IF (INCFD < 1) GO TO 5002
         DO I = 2, N
            IF (X(I) <= X(I-1)) GO TO 5003
         END DO
!
!  FUNCTION DEFINITION IS OK, GO ON.
!
         IERR = 0
         NLESS1 = N - 1
         H1 = X(2) - X(1)
         DEL1 = (F(1,2)-F(1,1))/H1
         DSAVE = DEL1
!
!  SPECIAL CASE N=2 -- USE LINEAR INTERPOLATION.
!
         IF (NLESS1 <= 1) THEN
            D(1,1) = DEL1
            D(1,N) = DEL1
         ELSE
!
!  NORMAL CASE  (N .GE. 3).
!
            H2 = X(3) - X(2)
            DEL2 = (F(1,3)-F(1,2))/H2
!
!  SET D(1) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
!  SHAPE-PRESERVING.
!
            HSUM = H1 + H2
            W1 = (H1 + HSUM)/HSUM
            W2 = -H1/HSUM
            D(1,1) = W1*DEL1 + W2*DEL2
            IF (PCHST(D(1,1),DEL1) <= ZERO) THEN
               D(1,1) = ZERO
            ELSE IF (PCHST(DEL1,DEL2) < ZERO) THEN
!        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
               DMAX = THREE*DEL1
               IF (ABS(D(1,1)) > ABS(DMAX)) D(1,1) = DMAX
            ENDIF
!
!  LOOP THROUGH INTERIOR POINTS.
!
            DO I = 2, NLESS1
               IF (I /= 2) THEN
!
                  H1 = H2
                  H2 = X(I+1) - X(I)
                  HSUM = H1 + H2
                  DEL1 = DEL2
                  DEL2 = (F(1,I+1)-F(1,I))/H2
               ENDIF
!
!        SET D(I)=0 UNLESS DATA ARE STRICTLY MONOTONIC.
!
               D(1,I) = ZERO
               R1 = PCHST(DEL1,DEL2)
               IF (R1 >= ZERO) THEN
                  IF (R1 > ZERO) GO TO 45
                  IF (DEL2 == ZERO) CYCLE 
                  IF (PCHST(DSAVE,DEL2) < ZERO) IERR = IERR + 1
                  DSAVE = DEL2
                  CYCLE 
!
               ENDIF
               IERR = IERR + 1
               DSAVE = DEL2
               CYCLE 
!
!        USE BRODLIE MODIFICATION OF BUTLAND FORMULA.
!
   45          CONTINUE
               HSUMT3 = HSUM + HSUM + HSUM
               W1 = (HSUM + H1)/HSUMT3
               W2 = (HSUM + H2)/HSUMT3
               DMAX = MAX(ABS(DEL1),ABS(DEL2))
               DMIN = MIN(ABS(DEL1),ABS(DEL2))
               DRAT1 = DEL1/DMAX
               DRAT2 = DEL2/DMAX
               D(1,I) = DMIN/(W1*DRAT1 + W2*DRAT2)
!
            END DO
!
!  SET D(N) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
!     SHAPE-PRESERVING.
!
            W1 = -H2/HSUM
            W2 = (H2 + HSUM)/HSUM
            D(1,N) = W1*DEL1 + W2*DEL2
            IF (PCHST(D(1,N),DEL2) <= ZERO) THEN
               D(1,N) = ZERO
            ELSE IF (PCHST(DEL1,DEL2) < ZERO) THEN
!        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
               DMAX = THREE*DEL2
               IF (ABS(D(1,N)) > ABS(DMAX)) D(1,N) = DMAX
            ENDIF
!
!  NORMAL RETURN.
!
         ENDIF
         RETURN 
!
!  ERROR RETURNS.
!
      ENDIF
!     N.LT.2 RETURN.
      IERR = -1
      STOP 'PCHIM -- NUMBER OF DATA POINTS LESS THAN TWO'
      RETURN 
!
 5002 CONTINUE
!     INCFD.LT.1 RETURN.
      IERR = -2
      STOP 'PCHIM -- INCREMENT LESS THAN ONE'
      RETURN 
!
 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      STOP 'PCHIM -- X-ARRAY NOT STRICTLY INCREASING'
      RETURN 
!------------- LAST LINE OF PCHIM FOLLOWS ------------------------------
      END SUBROUTINE PCHIM


      SUBROUTINE PCHSP(IC, VC, N, X, F, D, INCFD, WK, NWK, IERR)
      use kind_spec
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: N, INCFD, NWK, IERR
      INTEGER, DIMENSION(2) :: IC
      REAL(rprec), DIMENSION(2) :: VC
      REAL(rprec), DIMENSION(N) :: X
      REAL(rprec), DIMENSION(INCFD,N) :: F, D
      REAL(rprec), DIMENSION(2,N) :: WK
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IBEG, IEND, INDEX, J, NM1
      REAL(rprec) :: G, HALF, ONE
      REAL(rprec), DIMENSION(3) :: STEMP
      REAL(rprec) :: THREE, TWO
      REAL(rprec), DIMENSION(4) :: XTEMP
      REAL(rprec) :: ZERO
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      REAL(rprec) , EXTERNAL :: PCHDF
!-----------------------------------------------
!***BEGIN PROLOGUE  PCHSP
!***DATE WRITTEN   820503   (YYMMDD)
!***REVISION DATE  870707   (YYMMDD)
!***CATEGORY NO.  E1B
!***KEYWORDS  LIBRARY=SLATEC(PCHIP),
!             TYPE=SINGLE PRECISION(PCHSP-S DPCHSP-D),
!             CUBIC HERMITE INTERPOLATION,PIECEWISE CUBIC INTERPOLATION,
!             SPLINE INTERPOLATION
!***AUTHOR  FRITSCH, F. N., (LLNL)
!             MATHEMATICS AND STATISTICS DIVISION
!             LAWRENCE LIVERMORE NATIONAL LABORATORY
!             P.O. BOX 808  (L-316)
!             LIVERMORE, CA  94550
!             FTS 532-4275, (415) 422-4275
!***PURPOSE  Set derivatives needed to determine the Hermite represen-
!            tation of the cubic spline interpolant to given data, with
!            specified boundary conditions.
!***DESCRIPTION
!
!          PCHSP:   Piecewise Cubic Hermite Spline
!
!     Computes the Hermite representation of the cubic spline inter-
!     polant to the data given in X and F satisfying the boundary
!     conditions specified by IC and VC.
!
!     To facilitate two-dimensional applications, includes an increment
!     between successive values of the F- and D-arrays.
!
!     The resulting piecewise cubic Hermite function may be evaluated
!     by PCHFE or PCHFD.
!
!     NOTE:  This is a modified version of C. de Boor''S cubic spline
!            routine CUBSPL.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  IC(2), N, NWK, IERR
!        REAL  VC(2), X(N), F(INCFD,N), D(INCFD,N), WK(NWK)
!
!        CALL  PCHSP (IC, VC, N, X, F, D, INCFD, WK, NWK, IERR)
!
!   Parameters:
!
!     IC -- (input) integer array of length 2 specifying desired
!           boundary conditions:
!           IC(1) = IBEG, desired condition at beginning of data.
!           IC(2) = IEND, desired condition at end of data.
!
!           IBEG = 0  to set D(1) so that the third derivative is con-
!              tinuous at X(2).  This is the "not a knot" condition
!              provided by de Boor''S cubic spline routine CUBSPL.
!              < This is the default boundary condition. >
!           IBEG = 1  if first derivative at X(1) is given in VC(1).
!           IBEG = 2  if second derivative at X(1) is given in VC(1).
!           IBEG = 3  to use the 3-point difference formula for D(1).
!                     (Reverts to the default b.c. if N.LT.3 .)
!           IBEG = 4  to use the 4-point difference formula for D(1).
!                     (Reverts to the default b.c. if N.LT.4 .)
!          NOTES:
!           1. An error return is taken if IBEG is out of range.
!           2. For the "natural" boundary condition, use IBEG=2 and
!              VC(1)=0.
!
!           IEND may take on the same values as IBEG, but applied to
!           derivative at X(N).  In case IEND = 1 or 2, the value is
!           given in VC(2).
!
!          NOTES:
!           1. An error return is taken if IEND is out of range.
!           2. For the "natural" boundary condition, use IEND=2 and
!              VC(2)=0.
!
!     VC -- (input) real array of length 2 specifying desired boundary
!           values, as indicated above.
!           VC(1) need be set only if IC(1) = 1 or 2 .
!           VC(2) need be set only if IC(2) = 1 or 2 .
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!
!     X -- (input) real array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real array of dependent variable values to be inter-
!           polated.  F(1+(I-1)*INCFD) is value corresponding to X(I).
!
!     D -- (output) real array of derivative values at the data points.
!           These values will determine the cubic spline interpolant
!           with the requested boundary conditions.
!           The value corresponding to X(I) is stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!           No other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in F and D.
!           This argument is provided primarily for 2-D applications.
!           (Error return if  INCFD.LT.1 .)
!
!     WK -- (scratch) real array of working storage.
!
!     NWK -- (input) length of work array.
!           (Error return if NWK.LT.2*N .)
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -2  if INCFD.LT.1 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -4  if IBEG.LT.0 or IBEG.GT.4 .
!              IERR = -5  if IEND.LT.0 of IEND.GT.4 .
!              IERR = -6  if both of the above are true.
!              IERR = -7  if NWK is too small.
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!             (The D-array has not been changed in any of these cases.)
!              IERR = -8  in case of trouble solving the linear system
!                         for the interior derivative values.
!             (The D-array may have been changed in this case.)
!             (             Do **NOT** use it!                )
!
!***REFERENCES  CARL DE BOOR, A PRACTICAL GUIDE TO SPLINES, SPRINGER-
!                 VERLAG (NEW YORK, 1978), PP. 53-59.
!***ROUTINES CALLED  PCHDF
!***END PROLOGUE  PCHSP
!
! ----------------------------------------------------------------------
!
!  Change record:
!     82-08-04   Converted to SLATEC library version.
!     87-07-07   Minor cosmetic changes to prologue.
!
! ----------------------------------------------------------------------
!
!  Programming notes:
!
!     To produce a double precision version, simply:
!        a. Change PCHSP to DPCHSP wherever it occurs,
!        b. Change the real declarations to double precision, and
!        c. Change the constants ZERO, HALF, ... to double precision.
!
!  DECLARE ARGUMENTS.
!
!
!  DECLARE LOCAL VARIABLES.
!
!
      DATA ZERO/0/
      DATA HALF/0.5_DP/
      DATA ONE/1/
      DATA TWO/2/
      DATA THREE/3/
!
!  VALIDITY-CHECK ARGUMENTS.
!
!***FIRST EXECUTABLE STATEMENT  PCHSP
      IF (N >= 2) THEN
         IF (INCFD < 1) GO TO 5002
         DO J = 2, N
            IF (X(J) <= X(J-1)) GO TO 5003
         END DO
!
         IBEG = IC(1)
         IEND = IC(2)
         IERR = 0
         IF (IBEG<0 .OR. IBEG>4) IERR = IERR - 1
         IF (IEND<0 .OR. IEND>4) IERR = IERR - 2
         IF (IERR < 0) GO TO 5004
!
!  FUNCTION DEFINITION IS OK -- GO ON.
!
         IF (NWK < 2*N) GO TO 5007
!
!  COMPUTE FIRST DIFFERENCES OF X SEQUENCE AND STORE IN WK(1,.). ALSO,
!  COMPUTE FIRST DIVIDED DIFFERENCE OF DATA AND STORE IN WK(2,.).
         WK(1,2:N) = X(2:N) - X(:N-1)
         WK(2,2:N) = (F(1,2:N)-F(1,:N-1))/WK(1,2:N)
!
!  SET TO DEFAULT BOUNDARY CONDITIONS IF N IS TOO SMALL.
!
         IF (IBEG > N) IBEG = 0
         IF (IEND > N) IEND = 0
!
!  SET UP FOR BOUNDARY CONDITIONS.
!
         IF (IBEG==1 .OR. IBEG==2) THEN
            D(1,1) = VC(1)
         ELSE IF (IBEG > 2) THEN
!        PICK UP FIRST IBEG POINTS, IN REVERSE ORDER.
            DO J = 1, IBEG
               INDEX = IBEG - J + 1
!           INDEX RUNS FROM IBEG DOWN TO 1.
               XTEMP(J) = X(INDEX)
               IF (J < IBEG) STEMP(J) = WK(2,INDEX)
            END DO
!                 --------------------------------
            D(1,1) = PCHDF(IBEG,XTEMP,STEMP,IERR)
!                 --------------------------------
            IF (IERR /= 0) GO TO 5009
            IBEG = 1
         ENDIF
!
         IF (IEND==1 .OR. IEND==2) THEN
            D(1,N) = VC(2)
         ELSE IF (IEND > 2) THEN
!        PICK UP LAST IEND POINTS.
            DO J = 1, IEND
               INDEX = N - IEND + J
!           INDEX RUNS FROM N+1-IEND UP TO N.
               XTEMP(J) = X(INDEX)
               IF (J < IEND) STEMP(J) = WK(2,INDEX+1)
            END DO
!                 --------------------------------
            D(1,N) = PCHDF(IEND,XTEMP,STEMP,IERR)
!                 --------------------------------
            IF (IERR /= 0) GO TO 5009
            IEND = 1
         ENDIF
!
! --------------------( BEGIN CODING FROM CUBSPL )--------------------
!
!  **** A TRIDIAGONAL LINEAR SYSTEM FOR THE UNKNOWN SLOPES S(J) OF
!  F  AT X(J), J=1,...,N, IS GENERATED AND THEN SOLVED BY GAUSS ELIM-
!  INATION, WITH S(J) ENDING UP IN D(1,J), ALL J.
!     WK(1,.) AND WK(2,.) ARE USED FOR TEMPORARY STORAGE.
!
!  CONSTRUCT FIRST EQUATION FROM FIRST BOUNDARY CONDITION, OF THE FORM
!             WK(2,1)*S(1) + WK(1,1)*S(2) = D(1,1)
!
         IF (IBEG == 0) THEN
            IF (N == 2) THEN
!           NO CONDITION AT LEFT END AND N = 2.
               WK(2,1) = ONE
               WK(1,1) = ONE
               D(1,1) = TWO*WK(2,2)
            ELSE
!           NOT-A-KNOT CONDITION AT LEFT END AND N .GT. 2.
               WK(2,1) = WK(1,3)
               WK(1,1) = WK(1,2) + WK(1,3)
               D(1,1) = ((WK(1,2)+TWO*WK(1,1))*WK(2,2)
     1            *WK(1,3)+WK(1,2)**2*WK(2,3))/WK(1,1)
            ENDIF
         ELSE IF (IBEG == 1) THEN
!        SLOPE PRESCRIBED AT LEFT END.
            WK(2,1) = ONE
            WK(1,1) = ZERO
         ELSE
!        SECOND DERIVATIVE PRESCRIBED AT LEFT END.
            WK(2,1) = TWO
            WK(1,1) = ONE
            D(1,1) = THREE*WK(2,2) - HALF*WK(1,2)*D(1,1)
         ENDIF
!
!  IF THERE ARE INTERIOR KNOTS, GENERATE THE CORRESPONDING EQUATIONS AND
!  CARRY OUT THE FORWARD PASS OF GAUSS ELIMINATION, AFTER WHICH THE J-TH
!  EQUATION READS    WK(2,J)*S(J) + WK(1,J)*S(J+1) = D(1,J).
!
         NM1 = N - 1
         IF (NM1 > 1) THEN
            DO J = 2, NM1
               IF (WK(2,J-1) == ZERO) GO TO 5008
               G = -WK(1,J+1)/WK(2,J-1)
               D(1,J) = G*D(1,J-1) + THREE*(WK(1,J)
     1            *WK(2,J+1)+WK(1,J+1)*WK(2,J))
               WK(2,J) = G*WK(1,J-1) + TWO*(WK(1,J)+WK(1,J+1))
            END DO
         ENDIF
!
!  CONSTRUCT LAST EQUATION FROM SECOND BOUNDARY CONDITION, OF THE FORM
!           (-G*WK(2,N-1))*S(N-1) + WK(2,N)*S(N) = D(1,N)
!
!     IF SLOPE IS PRESCRIBED AT RIGHT END, ONE CAN GO DIRECTLY TO BACK-
!     SUBSTITUTION, SINCE ARRAYS HAPPEN TO BE SET UP JUST RIGHT FOR IT
!     AT THIS POINT.
         IF (IEND /= 1) THEN
!
            IF (IEND == 0) THEN
               IF (N==2 .AND. IBEG==0) THEN
!           NOT-A-KNOT AT RIGHT ENDPOINT AND AT LEFT ENDPOINT AND N = 2.
                  D(1,2) = WK(2,2)
                  GO TO 30
               ELSE IF (N==2 .OR. N==3 .AND. IBEG==0) THEN
!           EITHER (N=3 AND NOT-A-KNOT ALSO AT LEFT) OR (N=2 AND *NOT*
!           NOT-A-KNOT AT LEFT END POINT).
                  D(1,N) = TWO*WK(2,N)
                  WK(2,N) = ONE
                  IF (WK(2,N-1) == ZERO) GO TO 5008
                  G = -ONE/WK(2,N-1)
               ELSE
!           NOT-A-KNOT AND N .GE. 3, AND EITHER N.GT.3 OR  ALSO NOT-A-
!           KNOT AT LEFT END POINT.
                  G = WK(1,N-1) + WK(1,N)
!           DO NOT NEED TO CHECK FOLLOWING DENOMINATORS (X-DIFFERENCES).
                  D(1,N) = ((WK(1,N)+TWO*G)*WK(2,N)*WK(1,N-1)
     1            +WK(1,N)**2*(F(1,N-1)-F(1,N-2))/WK(1,N-1))/G
                  IF (WK(2,N-1) == ZERO) GO TO 5008
                  G = -G/WK(2,N-1)
                  WK(2,N) = WK(1,N-1)
               ENDIF
            ELSE
!        SECOND DERIVATIVE PRESCRIBED AT RIGHT ENDPOINT.
               D(1,N) = THREE*WK(2,N) + HALF*WK(1,N)*D(1,N)
               WK(2,N) = TWO
               IF (WK(2,N-1) == ZERO) GO TO 5008
               G = -ONE/WK(2,N-1)
            ENDIF
!
!  COMPLETE FORWARD PASS OF GAUSS ELIMINATION.
!
            WK(2,N) = G*WK(1,N-1) + WK(2,N)
            IF (WK(2,N) == ZERO) GO TO 5008
            D(1,N) = (G*D(1,N-1)+D(1,N))/WK(2,N)
!
!  CARRY OUT BACK SUBSTITUTION
!
         ENDIF
   30    CONTINUE
         DO J = NM1, 1, -1
            IF (WK(2,J) == ZERO) GO TO 5008
            D(1,J) = (D(1,J)-WK(1,J)*D(1,J+1))/WK(2,J)
         END DO
! --------------------(  END  CODING FROM CUBSPL )--------------------
!
!  NORMAL RETURN.
!
         RETURN 
!
!  ERROR RETURNS.
!
      ENDIF
!     N.LT.2 RETURN.
      IERR = -1
      STOP 'PCHSP -- NUMBER OF DATA POINTS LESS THAN TWO'
      RETURN 
!
 5002 CONTINUE
!     INCFD.LT.1 RETURN.
      IERR = -2
      STOP 'PCHSP -- INCREMENT LESS THAN ONE'
      RETURN 
!
 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      STOP 'PCHSP -- X-ARRAY NOT STRICTLY INCREASING'
      RETURN 
!
 5004 CONTINUE
!     IC OUT OF RANGE RETURN.
      IERR = IERR - 3
      STOP 'PCHSP -- IC OUT OF RANGE'
      RETURN 
!
 5007 CONTINUE
!     NWK TOO SMALL RETURN.
      IERR = -7
      STOP 'PCHSP -- WORK ARRAY TOO SMALL'
      RETURN 
!
 5008 CONTINUE
!     SINGULAR SYSTEM.
!   *** THEORETICALLY, THIS CAN ONLY OCCUR IF SUCCESSIVE X-VALUES   ***
!   *** ARE EQUAL, WHICH SHOULD ALREADY HAVE BEEN CAUGHT (IERR=-3). ***
      IERR = -8
      STOP 'PCHSP -- SINGULAR LINEAR SYSTEM'
      RETURN 
!
 5009 CONTINUE
!     ERROR RETURN FROM PCHDF.
!   *** THIS CASE SHOULD NEVER OCCUR ***
      IERR = -9
      STOP 'PCHSP -- ERROR RETURN FROM PCHDF'
      RETURN 
!------------- LAST LINE OF PCHSP FOLLOWS ------------------------------
      END SUBROUTINE PCHSP


      FUNCTION PCHST (ARG1, ARG2)
      USE kind_spec
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec) :: ARG1, ARG2
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec) :: ONE, ZERO, PCHST
!-----------------------------------------------
!***BEGIN PROLOGUE  PCHST
!***REFER TO  PCHCE,PCHCI,PCHCS,PCHIM
!***ROUTINES CALLED  (NONE)
!***DESCRIPTION
!
!         PCHST:  PCHIP Sign-Testing Routine.
!
!
!     Returns:
!        -1. if ARG1 and ARG2 are of opposite sign.
!         0. if either argument is zero.
!        +1. if ARG1 and ARG2 are of the same sign.
!
!     The object is to do this without multiplying ARG1*ARG2, to avoid
!     possible over/underflow problems.
!
!  Fortran intrinsics used:  SIGN.
!
! ----------------------------------------------------------------------
!
!  Programmed by:  Fred N. Fritsch,  FTS 532-4275, (415) 422-4275,
!                  Mathematics and Statistics Division,
!                  Lawrence Livermore National Laboratory.
!
!  Change record:
!     82-08-05   Converted to SLATEC library version.
!
! ----------------------------------------------------------------------
!
!  Programming notes:
!
!     To produce a double precision version, simply:
!        a. Change PCHST to DPCHST wherever it occurs,
!        b. Change all references to the Fortran intrinsics to their
!           double presision equivalents,
!        c. Change the real declarations to double precision, and
!        d. Change the constants  ZERO  and  ONE  to double precision.
!***END PROLOGUE  PCHST
!
!  DECLARE LOCAL VARIABLES.
!
      DATA ZERO/0/
      DATA ONE/1/
!
!  PERFORM THE TEST.
!
!***FIRST EXECUTABLE STATEMENT  PCHST
      PCHST = SIGN(ONE,ARG1)*SIGN(ONE,ARG2)
      IF (ARG1==ZERO .OR. ARG2==ZERO) PCHST = ZERO
!
      RETURN 
!------------- LAST LINE OF PCHST FOLLOWS ------------------------------
      END FUNCTION PCHST


      FUNCTION PCHDF (K, X, S, IERR)
      USE kind_spec
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: K, IERR
      REAL(rprec), DIMENSION(K) :: X, S
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, J
      REAL(rprec) :: VALUE, ZERO, PCHDF
!-----------------------------------------------
!***BEGIN PROLOGUE  PCHDF
!***REFER TO  PCHCE,PCHSP
!***DESCRIPTION
!
!          PCHDF:   PCHIP Finite Difference Formula
!
!     Uses a divided difference formulation to compute a K-point approx-
!     imation to the derivative at X(K) based on the data in X and S.
!
!     Called by  PCHCE  and  PCHSP  to compute 3- and 4-point boundary
!     derivative approximations.
!
! ----------------------------------------------------------------------
!
!     On input:
!        K      is the order of the desired derivative approximation.
!               K must be at least 3 (error return if not).
!        X      contains the K values of the independent variable.
!               X need not be ordered, but the values **MUST** be
!               distinct.  (Not checked here.)
!        S      contains the associated slope values:
!                  S(I) = (F(I+1)-F(I))/(X(I+1)-X(I)), I=1(1)K-1.
!               (Note that S need only be of length K-1.)
!
!     On return:
!        S      will be destroyed.
!        IERR   will be set to -1 if K.LT.2 .
!        PCHDF  will be set to the desired derivative approximation if
!               IERR=0 or to zero if IERR=-1.
!
! ----------------------------------------------------------------------
!
!  Reference:  Carl de Boor, A Practical Guide to Splines, Springer-
!              Verlag (New York, 1978), pp. 10-16.
!
! ----------------------------------------------------------------------
!
!  Programmed by:  Fred N. Fritsch,  FTS 532-4275, (415) 422-4275,
!                  Mathematics and Statistics Division,
!                  Lawrence Livermore National Laboratory.
!
!  Change record:
!     82-08-05   Converted to SLATEC library version.
!
! ----------------------------------------------------------------------
!
!  Programming notes:
!
!     To produce a double precision version, simply:
!        a. Change PCHDF to DPCHDF wherever it occurs,
!        b. Change the real declarations to double precision, and
!        c. Change the constant ZERO to double precision.
!***END PROLOGUE  PCHDF
!
!  DECLARE LOCAL VARIABLES.
!
      DATA ZERO/0/
!
!  CHECK FOR LEGAL VALUE OF K.
!
!***FIRST EXECUTABLE STATEMENT  PCHDF
      IF (K >= 3) THEN
!
!  COMPUTE COEFFICIENTS OF INTERPOLATING POLYNOMIAL.
!
         DO J = 2, K - 1
            S(:K-J) = (S(2:K-J+1)-S(:K-J))/(X(1+J:K)-X(:K-J))
         END DO
!
!  EVALUATE DERIVATIVE AT X(K).
!
         VALUE = S(1)
         DO I = 2, K - 1
            VALUE = S(I) + VALUE*(X(K)-X(I))
         END DO
!
!  NORMAL RETURN.
!
         IERR = 0
         PCHDF = VALUE
         RETURN 
!
!  ERROR RETURN.
!
      ENDIF
!     K.LT.3 RETURN.
      IERR = -1
      STOP 'PCHDF -- K LESS THAN THREE'
      PCHDF = ZERO
      RETURN 
!------------- LAST LINE OF PCHDF FOLLOWS ------------------------------
      END FUNCTION PCHDF


      SUBROUTINE PCHEV(N, X, F, D, NVAL, XVAL, FVAL, DVAL, IERR)
      USE kind_spec
      IMPLICIT NONE
!- ----------------------------------------------
!   D u m m y   A r g u m e n t s
!- ----------------------------------------------
      INTEGER :: N, NVAL, IERR
      REAL(rprec), DIMENSION(N) :: X, F, D
      REAL(rprec), DIMENSION(NVAL) :: XVAL, FVAL, DVAL
!- ----------------------------------------------
!   L o c a l   V a r i a b l e s
!- ----------------------------------------------
      INTEGER :: INCFD
      LOGICAL :: SKIP
!- ----------------------------------------------
!*  **BEGIN PROLOGUE  PCHEV
!*  **DATE WRITTEN   870828   (YYMMDD)
!*  **REVISION DATE  870828   (YYMMDD)
!*  **CATEGORY NO.  E3,H1
!*  **KEYWORDS  CUBIC HERMITE OR SPLINE DIFFERENTIATION,CUBIC HERMITE
!             EVALUATION,EASY TO USE SPLINE OR CUBIC HERMITE EVALUATOR
!*  **AUTHOR  KAHANER, D.K., (NBS)
!             SCIENTIFIC COMPUTING DIVISION
!             NATIONAL BUREAU OF STANDARDS
!             ROOM A161, TECHNOLOGY BUILDING
!             GAITHERSBURG, MARYLAND 20899
!             (301) 975-3808
!*  **PURPOSE  Evaluates the function and first derivative of a piecewise
!            cubic Hermite or spline function at an array of points XVAL,
!            easy to use.
!*  **DESCRIPTION
!
!          PCHEV:  Piecewise Cubic Hermite or Spline Derivative Evaluator,
!                  Easy to Use.
!
!     From the book "Numerical Methods and Software"
!          by  D. Kahaner, C.  Moler, S. Nash
!                 Prentice Hall 1988
!
!     Evaluates the function and first derivative of the cubic Hermite
!     or spline function defined by  N, X, F, D, at the array of points
!
!     This is an easy to use driver for the routines by F.N. Fritsch
!     described in reference (2) below. Those also have other capabilities.
!
! ----------------------------------------------------------------------
!
!  Calling sequence: CALL  PCHEV (N, X, F, D, NVAL, XVAL, FVAL, DVAL, IERR)
!
!     INTEGER N, NVAL, IERR
!     REAL X(N), F(N), D(N), XVAL(NVAL), FVAL(NVAL), DVAL(NVAL)
!
!   Parameters:
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!
!     X -- (input) real array of independent variable values.  The
!           elements of X must be strictly increasing:
!             X(I-1) .LT. X(I),  I = 2(1)N. (Error return if not.)
!
!     F -- (input) real array of function values.  F(I) is
!           the value corresponding to X(I).
!
!     D -- (input) real array of derivative values.  D(I) is
!           the value corresponding to X(I).
!
!  NVAL -- (input) number of points at which the functions are to be
!           evaluated. ( Error return if NVAL.LT.1 )
!
!  XVAL -- (input) real array of points at which the functions are to
!           be evaluated.
!
!          NOTES:
!           1. The evaluation will be most efficient if the elements
!              of XVAL are increasing relative to X;
!              that is,   XVAL(J) .GE. X(I)
!              implies    XVAL(K) .GE. X(I),  all K.GE.J .
!           2. If any of the XVAL are outside the interval [X(1),X(N)],
!              values are extrapolated from the nearest extreme cubic,
!              and a warning error is returned.
!
!  FVAL -- (output) real array of values of the cubic Hermite function
!           defined by  N, X, F, D  at the points  XVAL.
!
!  DVAL -- (output) real array of values of the first derivative of
!           the same function at the points  XVAL.
!
!  IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning error:
!              IERR.GT.0  means that extrapolation was performed at
!                 IERR points.
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -4  if NVAL.LT.1 .
!           (Output arrays have not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT*! been validated.
!              IERR = -5  if an error has occurred in the lower-level
!                         routine CHFDV.  NB: this should never happen.
!                         Notify the author **IMMEDIATELY*! if it does.
!
! ----------------------------------------------------------------------
!*  **REFERENCES  1. F.N.FRITSCH AND R.E.CARLSON, 'MONOTONE PIECEWISE CUBIC INTERPOLATION,'
!                   SIAM J.NUMER.ANAL. 17, 2 (APRIL 1980), 238-246.
!                2. F.N.FRITSCH, 'PIECEWISE CUBIC HERMITE INTERPOLATION PACKAGE, FINAL SPECIFICATIONS',
!                   LAWRENCE LIVERMORE NATIONAL LABORATORY, COMPUTER DOCUMENTATION UCID-30194, AUGUST 1982.
!*  **ROUTINES CALLED  PCHFD
!*  **END PROLOGUE  PCHEV
!
!  DECLARE LOCAL VARIABLES.
!
      DATA SKIP/.TRUE./
      DATA INCFD/1/
 
!
!
!***FIRST EXECUTABLE STATEMENT  PCHEV
!
      CALL PCHFD (N, X, F, D, INCFD, SKIP, NVAL, XVAL, FVAL,
     1 DVAL, IERR)
!
!
      RETURN 
!
!- ------------ LAST LINE OF PCHEV FOLLOWS ------------------------------
      END SUBROUTINE PCHEV


      SUBROUTINE PCHFD(N, X, F, D, INCFD, SKIP, NE, XE, FE, 
     1 DE, IERR)
      USE kind_spec
      IMPLICIT NONE
!- ----------------------------------------------
!   D u m m y   A r g u m e n t s
!- ----------------------------------------------
      INTEGER :: N, INCFD, NE, IERR
      LOGICAL :: SKIP
      REAL(rprec), DIMENSION(N) :: X
      REAL(rprec), DIMENSION(INCFD,N) :: F, D
      REAL(rprec), DIMENSION(NE) :: XE, FE, DE
!- ----------------------------------------------
!   L o c a l   V a r i a b l e s
!- ----------------------------------------------
      INTEGER :: I, IERC, IR, J, JFIRST
      INTEGER, DIMENSION(2) :: NEXT
      INTEGER :: NJ
!- ----------------------------------------------
!*  **BEGIN PROLOGUE  PCHFD
!*  **DATE WRITTEN   811020   (YYMMDD)
!*  **REVISION DATE  870707   (YYMMDD)
!*  **CATEGORY NO.  E3,H1
!*  **KEYWORDS  LIBRARY=SLATEC(PCHIP),
!             TYPE=SINGLE PRECISION(PCHFD-S DPCHFD-D),
!             CUBIC HERMITE DIFFERENTIATION,CUBIC HERMITE EVALUATION,
!             HERMITE INTERPOLATION,PIECEWISE CUBIC EVALUATION
!*  **AUTHOR  FRITSCH, F. N., (LLNL)
!             MATHEMATICS AND STATISTICS DIVISION
!             LAWRENCE LIVERMORE NATIONAL LABORATORY
!             P.O. BOX 808  (L-316)
!             LIVERMORE, CA  94550
!             FTS 532-4275, (415) 422-4275
!*  **PURPOSE  Evaluate a piecewise cubic hermite function and its first
!            derivative at an array of points.  May be used by itself
!            for Hermite interpolation, or as an evaluator for PCHIM
!            or PCHI!.   If only function values are required, use
!            PCHFE instead.
!*  **DESCRIPTION
!
!          PCHFD:  Piecewise Cubic Hermite Function and Derivative
!                  evaluator
!
!     Evaluates the cubic Hermite function defined by  N, X, F, D,  to-
!     gether with its first derivative, at the points  XE(J), J=1(1)NE.
!
!     If only function values are required, use PCHFE, instead.
!
!     To provide compatibility with PCHIM and PCHIC, includes an
!     increment between successive values of the F- and D-arrays.
!
C ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  N, NE, IERR
!        REAL X(N), F(INCFD,N), D(INCFD,N), XE(NE), FE(NE), DE(NE)
!        LOGICAL  SKIP
!
!        CALL  PCHFD (N, X, F, D, INCFD, SKIP, NE, XE, FE, DE, IERR)
!
!   Parameters:
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!
!     X -- (input) real array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real array of function values.  F(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     D -- (input) real array of derivative values.  D(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     INCFD -- (input) increment between successive values in F and D.
!           (Error return if  INCFD.LT.1 .)
!
!     SKIP -- (input/output) logical variable which should be set to
!           .TRUE. if the user wishes to skip checks for validity of
!           preceding parameters, or to .FALSE. otherwise.
!           This will save time in case these checks have already
!           been performed (say, in PCHIM or PCHIC).
!           SKIP will be set to .TRUE. on normal return.
!
!     NE -- (input) number of evaluation points.  (Error return if
!           NE.LT.1 .)
!
!     XE -- (input) real array of points at which the functions are to
!           be evaluated.
!
!
!          NOTES:
!           1. The evaluation will be most efficient if the elements
!              of XE are increasing relative to X;
!              that is,   XE(J) .GE. X(I)
!              implies    XE(K) .GE. X(I),  all K.GE.J .
!           2. If any of the XE are outside the interval [X(1),X(N)],
!              values are extrapolated from the nearest extreme cubic,
!              and a warning error is returned.
!
!     FE -- (output) real array of values of the cubic Hermite function
!           defined by  N, X, F, D  at the points  XE.
!
!     DE -- (output) real array of values of the first derivative of
!           the same function at the points  XE.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning error:
!              IERR.GT.0  means that extrapolation was performed at
!                 IERR points.
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -2  if INCFD.LT.1 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -4  if NE.LT.1 .
!           (Output arrays have not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT*! been validated.
!              IERR = -5  if an error has occurred in the lower-level
!                         routine CHFDV.  NB: this should never happen.
!                         Notify the author **IMMEDIATELY*! if it does.
!
!*  **REFERENCES  (NONE)
!*  **ROUTINES CALLED  CHFDV
!*  **END PROLOGUE  PCHFD
!
! ----------------------------------------------------------------------
!
!  Change record:
!     82-08-03   Minor cosmetic changes for release 1.
!     87-07-07   Minor cosmetic changes to prologue.
!
! ----------------------------------------------------------------------
!
!  Programming notes:
!
!     1. To produce a double precision version, simply:
!        a. Change PCHFD to DPCHFD, and CHFDV to DCHFDV, wherever they
!           occur,
!        b. Change the real declaration to double precision,
!
!     2. Most of the coding between the call to CHFDV and the end of
!        the IR-loop could be eliminated if it were permissible to
!        assume that XE is ordered relative to X.
!
!     3. CHFDV does not assume that X1 is less than X2.  thus, it would
!        be possible to write a version of PCHFD that assumes a strict-
!        ly decreasing X-array by simply running the IR-loop backwards
!        (and reversing the order of appropriate tests).
!
!     4. The present code has a minor bug, which I have decided is not
!        worth the effort that would be required to fix it.
!        If XE contains points in [X(N-1),X(N)], followed by points .LT.
!        X(N-1), followed by points .GT.X(N), the extrapolation points
!        will be counted (at least) twice in the total returned in IERR.
!
!  DECLARE ARGUMENTS.
!
!
!  DECLARE LOCAL VARIABLES.
!
!
!  VALIDITY-CHECK ARGUMENTS.
!
!*  **FIRST EXECUTABLE STATEMENT  PCHFD
      IF (.NOT.SKIP) THEN
!
         IF (N < 2) GO TO 5001
         IF (INCFD < 1) GO TO 5002
         DO I = 2, N
            IF (X(I) <= X(I-1)) GO TO 5003
         END DO
!
!  FUNCTION DEFINITION IS OK, GO ON.
!
      ENDIF
      IF (NE < 1) GO TO 5004
      IERR = 0
      SKIP = .TRUE.
!
!  LOOP OVER INTERVALS.        (   INTERVAL INDEX IS  IL = IR-1  . )
!                              ( INTERVAL IS X(IL).LE.X.LT.X(IR) . )
      JFIRST = 1
      IR = 2
   10 CONTINUE
!
!     SKIP OUT OF LOOP IF HAVE PROCESSED ALL EVALUATION POINTS.
!
      IF (JFIRST > NE) GO TO 5000
!
!     LOCATE ALL POINTS IN INTERVAL.
!
      DO J = JFIRST, NE
         IF (XE(J) >= X(IR)) GO TO 30
      END DO
      J = NE + 1
      GO TO 40
!
!     HAVE LOCATED FIRST POINT BEYOND INTERVAL.
!
   30 CONTINUE
      IF (IR == N) J = NE + 1
!
   40 CONTINUE
      NJ = J - JFIRST
!
!     SKIP EVALUATION IF NO POINTS IN INTERVAL.
!
      IF (NJ /= 0) THEN
!
!     EVALUATE CUBIC AT XE(I),  I = JFIRST (1) J-1 .
!
!       ----------------------------------------------------------------
         CALL CHFDV (X(IR-1), X(IR), F(1,IR-1), F(1,IR), 
     1      D(1,IR-1), D(1,IR), NJ, XE(JFIRST), FE(JFIRST),
     2      DE(JFIRST), NEXT, IERC)
!       ----------------------------------------------------------------
         IF (IERC < 0) GO TO 5005
!
         IF (NEXT(2) /= 0) THEN
!        IF (NEXT(2) .GT. 0)  THEN
!           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(2) TO THE
!           RIGHT OF X(IR).
!
            IF (IR >= N) THEN
!           IF (IR .EQ. N)  THEN
!              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
               IERR = IERR + NEXT(2)
            ELSE
!           ELSE
!              WE SHOULD NEVER HAVE GOTTEN HERE.
               GO TO 5005
!           ENDIF
!        ENDIF
            ENDIF
         ENDIF
!
         IF (NEXT(1) /= 0) THEN
!        IF (NEXT(1) .GT. 0)  THEN
!           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(1) TO THE
!           LEFT OF X(IR-1).
!
            IF (IR <= 2) THEN
!           IF (IR .EQ. 2)  THEN
!              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
               IERR = IERR + NEXT(1)
            ELSE
!           ELSE
!              XE IS NOT ORDERED RELATIVE TO X, SO MUST ADJUST
!              EVALUATION INTERVAL.
!
!              FIRST, LOCATE FIRST POINT TO LEFT OF X(IR-1).
               DO I = JFIRST, J - 1
                  IF (XE(I) < X(IR-1)) GO TO 45
               END DO
!              NOTE-- CANNOT DROP THROUGH HERE UNLESS THERE IS AN ERROR
!                     IN CHFDV.
               GO TO 5005
!
   45          CONTINUE
!              RESET J.  (THIS WILL BE THE NEW JFIRST.)
               J = I
!
!              NOW FIND OUT HOW FAR TO BACK UP IN THE X-ARRAY.
               DO I = 1, IR - 1
                  IF (XE(J) < X(I)) EXIT 
               END DO
!              AT THIS POINT, EITHER  XE(J) .LT. X(1)
!                 OR      X(I-1) .LE. XE(J) .LT. X(I) .
!              RESET IR, RECOGNIZING THAT IT WILL BE INCREMENTED BEFORE
!              CYCLING.
               IR = MAX(1,I - 1)
!           ENDIF
!        ENDIF
            ENDIF
         ENDIF
!
         JFIRST = J
!
!     END OF IR-LOOP.
!
      ENDIF
      IR = IR + 1
      IF (IR <= N) GO TO 10
!
!  NORMAL RETURN.
!
 5000 CONTINUE
      RETURN 
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     N.LT.2 RETURN.
      IERR = -1
      STOP 'PCHFD -- NUMBER OF DATA POINTS LESS THAN TWO'
      RETURN 
!
 5002 CONTINUE
!     INCFD.LT.1 RETURN.
      IERR = -2
      STOP 'PCHFD -- INCREMENT LESS THAN ONE'
      RETURN 
!
 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      STOP 'PCHFD -- X-ARRAY NOT STRICTLY INCREASING'
      RETURN 
!
 5004 CONTINUE
!     NE.LT.1 RETURN.
      IERR = -4
      STOP 'PCHFD -- NUMBER OF EVALUATION POINTS LESS THAN ONE'
      RETURN 
!
 5005 CONTINUE
!     ERROR RETURN FROM CHFDV.
!   **! THIS CASE SHOULD NEVER OCCUR ***
      IERR = -5
      STOP 'PCHFD -- ERROR RETURN FROM CHFDV -- FATAL'
      RETURN 
!- ------------ LAST LINE OF PCHFD FOLLOWS ------------------------------
      END SUBROUTINE PCHFD


      SUBROUTINE CHFDV(X1,X2,F1,F2,D1,D2,NE,XE,FE,DE,NEXT,IERR)
      USE kind_spec
      IMPLICIT NONE
!- ----------------------------------------------
!   D u m m y   A r g u m e n t s
!- ----------------------------------------------
      INTEGER :: NE, IERR
      REAL(rprec):: X1, X2, F1, F2, D1, D2
      INTEGER, DIMENSION(2) :: NEXT
      REAL(rprec), DIMENSION(NE) :: XE, FE, DE
!- ----------------------------------------------
!   L o c a l   V a r i a b l e s
!- ----------------------------------------------
      INTEGER :: I
      REAL(rprec) :: C2,C2T2,C3,C3T3,DEL1,DEL2,DELTA,H,
     1 X,XMI,XMA,ZERO
!- ----------------------------------------------
!*  **BEGIN PROLOGUE  CHFDV
!*  **DATE WRITTEN   811019   (YYMMDD)
!*  **REVISION DATE  870707   (YYMMDD)
!*  **CATEGORY NO.  E3,H1
!*  **KEYWORDS  LIBRARY=SLATEC(PCHIP),
!             TYPE=SINGLE PRECISION(CHFDV-S DCHFDV-D),
!             CUBIC HERMITE DIFFERENTIATION,CUBIC HERMITE EVALUATION,
!             CUBIC POLYNOMIAL EVALUATION
!*  **AUTHOR  FRITSCH, F. N., (LLNL)
!             MATHEMATICS AND STATISTICS DIVISION
!             LAWRENCE LIVERMORE NATIONAL LABORATORY
!             P.O. BOX 808  (L-316)
!             LIVERMORE, CA  94550
!             FTS 532-4275, (415) 422-4275
!*  **PURPOSE  Evaluate a cubic polynomial given in Hermite form and its
!            first derivative at an array of points.  While designed for
!            use by PCHFD, it may be useful directly as an evaluator for
!            a piecewise cubic Hermite function in applications, such as
!            graphing, where the interval is known in advance.
!            If only function values are required, use CHFEV instead.
!*  **DESCRIPTION
!
!        CHFDV:  Cubic Hermite Function and Derivative Evaluator
!
!     Evaluates the cubic polynomial determined by function values
!     F1,F2 and derivatives D1,D2 on interval (X1,X2), together with
!     its first derivative, at the points  XE(J), J=1(1)NE.
!
!     If only function values are required, use CHFEV, instead.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        INTEGER  NE, NEXT(2), IERR
!        REAL X1, X2, F1, F2, D1, D2, XE(NE), FE(NE), DE(NE)
!
!        CALL  CHFDV (X1,X2, F1,F2, D1,D2, NE, XE, FE, DE, NEXT, IERR)
!
!   Parameters:
!
!     X1,X2 -- (input) endpoints of interval of definition of cubic.
!           (Error return if  X1.EQ.X2 .)
!
!     F1,F2 -- (input) values of function at X1 and X2, respectively.
!
!     D1,D2 -- (input) values of derivative at X1 and X2, respectively.
!
!     NE -- (input) number of evaluation points.  (Error return if
!           NE.LT.1 .)
!
!     XE -- (input) real array of points at which the functions are to
!           be evaluated.  If any of the XE are outside the interval
!           [X1,X2], a warning error is returned in NEXT.
!
!     FE -- (output) real array of values of the cubic function defined
!           by  X1,X2, F1,F2, D1,D2  at the points  XE.
!
!     DE -- (output) real array of values of the first derivative of
!           the same function at the points  XE.
!
!     NEXT -- (output) integer array indicating number of extrapolation
!           points:
!            NEXT(1) = number of evaluation points to left of interval.
!            NEXT(2) = number of evaluation points to right of interval.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -1  if NE.LT.1 .
!              IERR = -2  if X1.EQ.X2 .
!                (Output arrays have not been changed in either case.)
!
!*  **REFERENCES  (NONE)
!*  **END PROLOGUE  CHFDV
!
! ----------------------------------------------------------------------
!
!  Change record:
!     82-08-03   Minor cosmetic changes for release 1.
!
! ----------------------------------------------------------------------
!
!  Programming notes:
!
!     To produce a double precision version, simply:
!        a. Change CHFDV to DCHFDV wherever it occurs,
!        b. Change the real declaration to double precision,
!        c. Change the constant ZERO to double precision, and
!        d. Change the names of the Fortran functions:  MAX, MIN.
!
!  DECLARE ARGUMENTS.
!
!
!  DECLARE LOCAL VARIABLES.
!
      DATA ZERO/0/
!
!  VALIDITY-CHECK ARGUMENTS.
!
!*  **FIRST EXECUTABLE STATEMENT  CHFDV
      IF (NE >= 1) THEN
         H = X2 - X1
         IF (H == ZERO) GO TO 5002
!
!  INITIALIZE.
!
         IERR = 0
         NEXT(1) = 0
         NEXT(2) = 0
         XMI = MIN(ZERO,H)
         XMA = MAX(ZERO,H)
!
!  COMPUTE CUBIC COEFFICIENTS (EXPANDED ABOUT X1).
!
         DELTA = (F2 - F1)/H
         DEL1 = (D1 - DELTA)/H
         DEL2 = (D2 - DELTA)/H
!                                           (DELTA IS NO LONGER NEEDED.)
         C2 = -(DEL1 + DEL1 + DEL2)
         C2T2 = C2 + C2
         C3 = (DEL1 + DEL2)/H
!                               (H, DEL1 AND DEL2 ARE NO LONGER NEEDED.)
         C3T3 = C3 + C3 + C3
!
!  EVALUATION LOOP.
!
         DO I = 1, NE
            X = XE(I) - X1
            FE(I) = F1 + X*(D1 + X*(C2 + X*C3))
            DE(I) = D1 + X*(C2T2 + X*C3T3)
!          COUNT EXTRAPOLATION POINTS.
            IF (X < XMI) NEXT(1) = NEXT(1) + 1
            IF (X > XMA) NEXT(2) = NEXT(2) + 1
!        (NOTE REDUNDANCY--IF EITHER CONDITION IS TRUE, OTHER IS FALSE.)
         END DO
!
!  NORMAL RETURN.
!
         RETURN 
!
!  ERROR RETURNS.
!
      ENDIF
!     NE.LT.1 RETURN.
      IERR = -1
      STOP 'CHFDV -- NUMBER OF EVALUATION POINTS LESS THAN ONE'
      RETURN 
!
 5002 CONTINUE
!     X1.EQ.X2 RETURN.
      IERR = -2
      STOP 'CHFDV -- INTERVAL ENDPOINTS EQUAL'
      RETURN 
!- ------------ LAST LINE OF CHFDV FOLLOWS ------------------------------
      END SUBROUTINE CHFDV
